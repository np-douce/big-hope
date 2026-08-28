# NP-DOUCE Ising Lab

Offline Ising spin-glass ground-state solver using exact cumulants through fifth order, a Taylor expansion of the partition function, Boltzmann-based spin importance, exact preprocessing reductions, and optional backtracking.

## Live app

**GitHub Pages:** https://np-douce.github.io/big-hope/

The app runs entirely in the browser with no backend, API keys, or CDN dependencies.

## What this project explores

The solver works with weighted Ising models of the form

\[
H(\mathbf{s})=-\sum_{i<j}J_{ij}s_i s_j-\sum_i h_i s_i,
\qquad s_i\in\{-1,+1\}.
\]

It approximates the partition function using a Taylor/cumulant expansion through order five:

\[
\ln Z^{(5)}(\beta)
=
N\ln 2
-\beta\kappa_1
+\frac{\beta^2}{2}\kappa_2
-\frac{\beta^3}{6}\kappa_3
+\frac{\beta^4}{24}\kappa_4
-\frac{\beta^5}{120}\kappa_5.
\]

The cumulants through \(\kappa_5\) are calculated exactly from the interaction graph at \(\beta=0\), using the ghost-spin construction for local fields. For fixed Taylor order five, this calculation remains polynomial in the number of spins.

For each unresolved spin, the search compares the two conditioned approximations

\[
I_i^{(r)}
=
\ln Z^{(r)}(s_i=+1)
-
\ln Z^{(r)}(s_i=-1),
\]

and uses the score to guide the next spin assignment.

## Exact preprocessing

Before each Taylor decision, the app applies fast exact reductions when they are mathematically guaranteed, including:

- zero-coupling removal;
- duplicate-coupling merging;
- self-coupling conversion to a constant energy offset;
- disconnected-component detection;
- isolated-spin fixing;
- zero-field global spin symmetry reduction;
- dominant local-field fixing;
- propagation of fixed spins into effective local fields;
- bipartite and unfrustrated-gauge diagnostics.

These checks simplify the model but do not replace the Taylor/cumulant search for genuinely unresolved spins.

## Ground-state search

The solver supports:

- Taylor orders 2 through 5;
- weighted positive and negative couplings;
- nonzero local fields;
- arbitrary sparse or nonplanar interaction graphs;
- spin-glass instances;
- candidate selection by maximum \(|I_i|\);
- exact final-energy evaluation;
- optional alternate-branch backtracking;
- brute-force validation for small instances;
- interactive visualization of the spin graph and search state.

The fifth-order cumulants are exact, but the truncated \(\ln Z^{(5)}\) is an approximation to the full partition function. The project therefore treats the Taylor-guided ground-state method as an experimental search algorithm rather than a proof of polynomial-time exact ground-state optimization.

## Files

- `index.html` - GitHub Pages/local launcher.
- `ising.html` - main offline Ising solver UI.
- `styles.css` - local NP-DOUCE styling.
- `ising-app.js` - browser runtime.
- `ising-math.js` - Hamiltonian evaluation, ghost-spin graph, cycle enumeration, cumulants, and Taylor \(\ln Z\).
- `ising-search.js` - conditioning, importance scoring, Taylor-guided search, and backtracking.
- `ising-precheck.js` - exact preprocessing reductions and diagnostics.
- `ising-examples.js` - built-in examples, including the 50-node spin glass.
- `ising-tests.js` - developer sanity checks.
- `sw.js` and `manifest.webmanifest` - offline/PWA support.
- `examples/` - plain-text input examples.

## Run locally

Open `ising.html` directly from the local folder, or open `index.html` to launch it from the repository root.

With the bundled Node runtime, the folder can also be served with:

```powershell
node server.cjs
```

Developer checks:

```powershell
node ising-tests.js
```

## Suggested GitHub topics

`ising-model` · `spin-glass` · `ground-state` · `qubo` · `cumulant-expansion` · `high-temperature-expansion` · `statistical-mechanics` · `partition-function` · `boltzmann-distribution` · `combinatorial-optimization` · `optimization` · `ising-solver` · `taylor-expansion`

## Project keywords

Ising model, spin glass, ground state, statistical mechanics, partition function, cumulants, Taylor expansion, Boltzmann distribution, QUBO, combinatorial optimization, nonplanar graphs, offline solver, NP-DOUCE.
