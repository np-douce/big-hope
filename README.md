# NP-douce Ising Lab

Standalone offline browser app for weighted Ising ground-state experiments using Taylor/cumulant partition-function approximations through kappa5.

Open `ising.html` directly from the local folder, or open `index.html` to launch it from a repository root. The app has no backend, no API keys, and no CDN dependencies.

With the bundled Codex Node runtime, this folder can be served with:

```powershell
node server.cjs
```

## Files

- `index.html` - GitHub Pages/local launcher that opens the app.
- `ising.html` - UI shell modeled after the `np-douce/website` static app.
- `styles.css` - local dark teal/gold NP-douce styling.
- `ising-app.js` - single-file browser runtime so `ising.html` works from a local folder.
- `ising-math.js` - Hamiltonian evaluation, ghost-spin graph, cycle enumeration, cumulants, and Taylor ln Z.
- `ising-search.js` - spin conditioning, importance scoring, Taylor-guided search, and optional alternate-branch backtracking.
- `ising-precheck.js` - exact preprocessing reductions and diagnostics.
- `ising-examples.js` - built-in pasteable examples, including the default 50-node spin glass.
- `ising-tests.js` - optional developer sanity checks outside the app.
- `sw.js` and `manifest.webmanifest` - offline/PWA support.
- `examples/` - plain text inputs for the required built-in test problems.

## Local developer sanity checks

```powershell
node ising-tests.js
```

Use the `examples/` files for pasteable test data. External tools can be used for independent validation.

## GitHub Pages

GitHub Pages supports one user or organization site named `<owner>.github.io`, plus one project site for each repository. This app can be published as a project page without replacing another website. Enable Pages from the repository root on the main branch, and `index.html` will launch `ising.html`.
