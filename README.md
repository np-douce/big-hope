# NP-DOUCE Ising Lab

Offline Ising spin-glass ground-state solver using cumulants through fifth order, a Taylor expansion of the partition function, Boltzmann-based spin importance, exact preprocessing reductions, and optional backtracking.

## Live app

**GitHub Pages:** https://np-douce.github.io/big-hope/

The app runs entirely in the browser with no backend, API keys, or CDN dependencies. Open `ising.html` directly from the local folder, or open `index.html` to launch it from the repository root.

## What this project explores

The solver works with weighted Ising models of the form

```text
H(s) = - sum_{i<j} J_ij s_i s_j - sum_i h_i s_i,
where s_i is in {-1, +1}.
```

It approximates the partition function using a Taylor/cumulant expansion through order five:

```text
ln Z^(5)(beta)
= N ln 2
- beta kappa1
+ beta^2 kappa2 / 2
- beta^3 kappa3 / 6
+ beta^4 kappa4 / 24
- beta^5 kappa5 / 120.
```

For each unresolved spin, max mode compares both conditioned approximations:

```text
I_i = ln Z^(5)(s_i = +1) - ln Z^(5)(s_i = -1)
```

and scores every remaining spin before choosing the next assignment. There is no candidate cap or separate confirmation path in the app.

## Speedups

The math surface is unchanged, but the implementation avoids unnecessary work:

- direct fixed-order sums replace generic cycle canonicalization;
- sparse anchored cycle walks avoid scanning absent dense tuples on sparse graphs;
- cycle products are accumulated immediately instead of stored in arrays;
- each model keeps cached edge, adjacency, degree, and radius metadata;
- conditioned scoring updates local fields in an augmented graph instead of rebuilding full remapped models twice per candidate;
- dense prechecks skip structural diagnostics that cannot help dense connected graphs.

## Exact Preprocessing

Before each Taylor decision, the app applies exact reductions when they are mathematically guaranteed, including:

- zero-coupling removal;
- duplicate-coupling merging;
- self-coupling conversion to a constant energy offset;
- isolated-spin fixing;
- zero-field global spin symmetry reduction;
- dominant local-field fixing;
- propagation of fixed spins into effective local fields;
- bipartite and unfrustrated-gauge diagnostics when useful.

These checks simplify the model but do not replace the Taylor/cumulant search for genuinely unresolved spins.

## Files

- `index.html` - GitHub Pages/local launcher.
- `ising.html` - main offline Ising solver UI.
- `styles.css` - local NP-DOUCE styling.
- `ising-app.js` - browser runtime.
- `ising-math.js` - Hamiltonian evaluation, ghost-spin graph, cumulants, and Taylor ln Z.
- `ising-search.js` - conditioning, importance scoring, Taylor-guided search, and backtracking.
- `ising-precheck.js` - preprocessing reductions and diagnostics.
- `ising-examples.js` - built-in examples, including the 50-node spin glass.
- `sw.js` and `manifest.webmanifest` - offline/PWA support.
- `examples/` - plain-text benchmark-style inputs.

## Run Locally

Open `ising.html` directly from the local folder, or serve the folder with:

```powershell
node server.cjs
```

## Suggested GitHub Topics

`ising-model` · `spin-glass` · `ground-state` · `qubo` · `cumulant-expansion` · `high-temperature-expansion` · `statistical-mechanics` · `partition-function` · `boltzmann-distribution` · `combinatorial-optimization` · `optimization` · `ising-solver` · `taylor-expansion`

## Project Keywords

Ising model, spin glass, ground state, statistical mechanics, partition function, cumulants, Taylor expansion, Boltzmann distribution, QUBO, combinatorial optimization, nonplanar graphs, offline solver, NP-DOUCE.
