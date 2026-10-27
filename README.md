# DermTrace GitHub Pages

This branch is served at <https://willibmwang.github.io/DermTrace/>.

- `index.html`: the DermTrace landing page. Its source is `website/` on `main`.
- `favicon.svg`: the DermTrace mark.
- `404.html`: not-found page; also forwards legacy root-level `/DermTrace/assets/...`
  links into `keloidbench/`, where the old KeloidBench site now lives.
- `keloidbench/`: the archived KeloidBench site (NeurIPS GenAI4Health 2026),
  KeloidBench being the keloid-only predecessor of DermTrace. It is a snapshot
  of `website/` at commit `be60f64`, with its repository links pinned to that
  commit. Treat it as immutable — do not regenerate it from `main`.
  See `keloidbench/README.md`.

Earlier versions of the KeloidBench site, which were served from the branch
root, remain in this branch's history.
