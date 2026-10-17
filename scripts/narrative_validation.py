"""Conservative narrative checks; semantic faithfulness still requires human review."""
import re

DISCLAIMER = 'This is a research prototype and must not be interpreted as a clinical diagnosis.'

def factual_anchor(packet):
    p = packet['prediction']
    if p['decision'] == 'abstain':
        return (f"KeloidBench abstained because model confidence ({p['confidence']:.3f}) "
                f"did not reach the operational threshold ({p['confidence_threshold']:.2f}).")
    return (f"KeloidBench returned a {p['decision'].replace('_', ' ')} research-triage result "
            f"with P(keloid)={p['prob_keloid']:.3f} and model confidence={p['confidence']:.3f}.")

TOWARD_NON = re.compile(r'\b(?:non[- ]keloid|unaffected|against|away|oppos\w*|decreas\w*|lower\w*|reduc\w*|negative\w*)\b', re.I)
TOWARD_KELOID = re.compile(r'\b(?:towards?|favou?r\w*|support\w*|increas\w*|rais\w*|positive\w*)\b', re.I)
CLAUSE = re.compile(r'[,;:]|\b(?:while|whereas|but|however|although|in contrast)\b|(?<=[.!?])\s+', re.I)

def contribution_roles_ok(text, packet):
    """Clause-level check that each named program keeps its supplied direction.

    Added in narrative_local_v2 after v1 showed the toward/against phrase check
    was too weak. Contributions are oriented toward P(keloid). A clause naming
    programs of both directions together with a direction word is ambiguous and fails.
    """
    roles = {}
    for key, role in [('top_supporting_contributions', 'keloid'), ('top_opposing_contributions', 'non')]:
        for item in packet.get(key, []):
            name = item.get('display_name', '').lower()
            for alias in {name, name.split(' ')[0], item.get('feature', '').lower().replace('_', ' ')}:
                if alias: roles.setdefault(alias, set()).add(role)
    for clause in CLAUSE.split(text.lower().replace('_', ' ')):
        named = set().union(*[r for alias, r in roles.items() if re.search(r'\b' + re.escape(alias) + r'\b', clause)] or [set()])
        if not named:
            continue
        direction = 'non' if TOWARD_NON.search(clause) else 'keloid' if TOWARD_KELOID.search(clause) else None
        if direction and (len(named) > 1 or direction not in named):
            return False
    return True

def validate_narrative(narrative, packet):
    anchor = factual_anchor(packet)
    text = narrative.strip()
    anchor_ok = text.startswith(anchor)
    # Strip only the exact allowed factual anchor and disclaimer, not arbitrary
    # sentences mentioning research use (which could also contain harmful claims).
    body = text[len(anchor):] if anchor_ok else text
    body = body.replace(DISCLAIMER, '')
    qc_note = packet.get('quality_control',{}).get('transfer_note','')
    stripped = body.replace(qc_note,'') if qc_note else body
    prohibited = re.search(r'\b(diagnos\w*|prescrib\w*|therapy|therapeutic|treat(?:ment)? with|causes?|drives?|leads to|recurrence risk|prognos\w*|doi|PMID|et al\.)\b', stripped, re.I)
    # A contradictory class assertion is disallowed even when both class names occur.
    contradiction = re.search(r'\b(?:is|was|has|indicates|confirms|diagnosed as|classified as|predicts|prediction is)\s+(?:a\s+)?(?:non[- ]keloid|unaffected|keloid)\b', stripped, re.I)
    assertive_abstention = packet['prediction']['decision']=='abstain' and bool(re.search(r'\b(?:definitely|certainly|confirmed|confidently)\b',stripped,re.I))
    # No additional numbers: the prompt requires the labeled anchor, so swapped
    # numbers cannot pass merely because both values occurred somewhere in JSON.
    extra_numbers = bool(re.search(r'(?<![A-Za-z])[-+]?\d+(?:\.\d+)?',stripped))
    direction_ok = True
    contributions = packet.get('top_supporting_contributions',[])+packet.get('top_opposing_contributions',[])
    for item in contributions:
        name = item.get('display_name','')
        value = float(item.get('contribution',0))
        # Clause-level (v2): a correct sentence may name programs in both directions.
        for sentence in CLAUSE.split(stripped):
            if name and name.lower() in sentence.lower():
                toward = bool(re.search(r'(toward|support\w*|increas\w*)\s+(?:the\s+)?keloid',sentence,re.I))
                against = bool(re.search(r'(against|away from|decreas\w*)\s+(?:the\s+)?keloid',sentence,re.I))
                if (toward and value < 0) or (against and value > 0): direction_ok=False
    checks = {'factual_anchor':anchor_ok, 'numeric_fields_grounded':anchor_ok and not extra_numbers,
              'decision_consistent':anchor_ok and not contradiction and not assertive_abstention,
              'abstention_disclosed':packet['prediction']['decision']!='abstain' or (anchor_ok and 'abstained' in anchor),
              'no_prohibited_claims':not bool(prohibited), 'research_use_stated':text.endswith(DISCLAIMER),
              'qc_disclosed':bool(qc_note) and qc_note in text, 'contribution_direction':direction_ok,
              'nonempty_and_within_word_limit':bool(text) and len(text.split())<=130,
              'contribution_role':contribution_roles_ok(stripped, packet)}
    return {'checks':checks,'automated_pass':all(checks.values()),
            'failed_checks':[k for k,v in checks.items() if not v],
            'limitation':'Automated screening is conservative and incomplete; it is not a semantic faithfulness guarantee.'}
