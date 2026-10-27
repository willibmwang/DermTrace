# DermTrace site source

Source of the project site served at <https://willibmwang.github.io/DermTrace/>.

| File | Purpose |
|---|---|
| `index.html` | The landing page. Self-contained: inline CSS, one inline SVG, no build step. |
| `404.html` | Not-found page. Also forwards legacy `/DermTrace/assets/...` deep links to the KeloidBench archive. |
| `favicon.svg` | DermTrace mark. |

Preview locally:

```bash
python3 -m http.server 8765 --directory website
```

## Publishing

The site is served from the **`gh-pages` branch**, not from `main`, and nothing in CI
deploys it. Copy these files to the root of `gh-pages` and commit there.

The page deliberately keeps its assets at the branch root (`favicon.svg`) rather than
under `assets/`, so it does not shadow the legacy KeloidBench deep links that
`404.html` forwards into `keloidbench/`.

## Do not regenerate the archive

`gh-pages:keloidbench/` is a frozen snapshot of `keloidbench/website/` taken at commit
`be60f64`, with an added archive banner and every repository link rewritten to
`willibmwang/DermTrace` pinned at that commit. The source on `main` still points at the
pre-rename `willibmwang/KeloidBench` repository, so re-syncing it from `main` would
silently break the archive. Treat that folder as immutable.
