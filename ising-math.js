export const EPS = 1e-12;

export function edgeKey(a, b) {
  return a < b ? `${a}:${b}` : `${b}:${a}`;
}

export function cloneModel(model) {
  return {
    n: model.n,
    labels: [...model.labels],
    fields: [...model.fields],
    couplings: new Map(model.couplings),
    offset: model.offset || 0,
    normalizationStats: { ...(model.normalizationStats || {}) }
  };
}

export function createModel(n, couplings = [], fields = []) {
  const model = {
    n,
    labels: Array.from({ length: n }, (_, i) => i),
    fields: Array(n).fill(0),
    couplings: new Map(),
    offset: 0,
    normalizationStats: {
      zeroCouplingsRemoved: 0,
      duplicateEdgesMerged: 0,
      selfCouplingsConverted: 0
    }
  };
  fields.forEach((value, i) => {
    model.fields[i] = Number(value) || 0;
  });
  for (const c of couplings) {
    const i = c.i ?? c[0];
    const j = c.j ?? c[1];
    const value = c.J ?? c.value ?? c[2];
    if (Math.abs(value) < EPS) {
      model.normalizationStats.zeroCouplingsRemoved++;
      continue;
    }
    if (i === j) {
      model.offset -= Number(value);
      model.normalizationStats.selfCouplingsConverted++;
      continue;
    }
    const key = edgeKey(i, j);
    if (model.couplings.has(key)) model.normalizationStats.duplicateEdgesMerged++;
    model.couplings.set(key, (model.couplings.get(key) || 0) + Number(value));
  }
  for (const [key, value] of [...model.couplings.entries()]) {
    if (Math.abs(value) < EPS) {
      model.couplings.delete(key);
      model.normalizationStats.zeroCouplingsRemoved++;
    }
  }
  return model;
}

export function parseIsingText(text) {
  const lines = text.split(/\r?\n/).map((line) => line.replace(/#.*/, "").trim()).filter(Boolean);
  let n = 0;
  let section = "";
  const fields = [];
  const couplings = [];

  for (const line of lines) {
    const upper = line.toUpperCase();
    if (upper === "FIELDS" || upper === "COUPLINGS") {
      section = upper;
      continue;
    }
    const parts = line.split(/\s+/);
    if (parts[0].toUpperCase() === "N") {
      n = Number(parts[1]);
      continue;
    }
    if (section === "FIELDS") {
      const i = Number(parts[0]) - 1;
      fields[i] = Number(parts[1]);
    } else if (section === "COUPLINGS") {
      couplings.push({ i: Number(parts[0]) - 1, j: Number(parts[1]) - 1, J: Number(parts[2]) });
    }
  }

  if (!Number.isInteger(n) || n < 1) throw new Error("Input must include a positive `N` line.");
  return createModel(n, couplings, Array.from({ length: n }, (_, i) => fields[i] || 0));
}

export function serializeModel(model) {
  const fields = model.fields
    .map((h, i) => `${i + 1} ${formatNumber(h)}`)
    .join("\n");
  const couplings = [...model.couplings.entries()]
    .map(([key, value]) => {
      const [i, j] = key.split(":").map(Number);
      return `${i + 1} ${j + 1} ${formatNumber(value)}`;
    })
    .join("\n");
  return `N ${model.n}\n\nFIELDS\n${fields}\n\nCOUPLINGS\n${couplings}`;
}

export function buildAugmentedGraph(model) {
  const size = model.n + 1;
  const weights = Array.from({ length: size }, () => Array(size).fill(0));
  const adjacency = Array.from({ length: size }, () => []);

  const addEdge = (a, b, value) => {
    if (Math.abs(value) < EPS) return;
    weights[a][b] = value;
    weights[b][a] = value;
    adjacency[a].push(b);
    adjacency[b].push(a);
  };

  model.fields.forEach((h, i) => addEdge(0, i + 1, h));
  for (const [key, value] of model.couplings.entries()) {
    const [i, j] = key.split(":").map(Number);
    addEdge(i + 1, j + 1, value);
  }
  adjacency.forEach((neighbors) => neighbors.sort((a, b) => a - b));
  return { size, weights, adjacency };
}

function canonicalCycle(path) {
  const n = path.length;
  let best = null;
  const variants = [path, [...path].reverse()];
  for (const seq of variants) {
    for (let shift = 0; shift < n; shift++) {
      const rotated = [...seq.slice(shift), ...seq.slice(0, shift)];
      const key = rotated.join(",");
      if (best === null || key < best) best = key;
    }
  }
  return best;
}

export function enumerateCycles(graph, length) {
  const cycles = new Map();
  const { size, adjacency, weights } = graph;
  const walk = (start, path, used) => {
    const last = path[path.length - 1];
    if (path.length === length) {
      if (weights[last][start] !== 0) cycles.set(canonicalCycle(path), [...path]);
      return;
    }
    for (const next of adjacency[last]) {
      if (next === start || used.has(next)) continue;
      used.add(next);
      path.push(next);
      walk(start, path, used);
      path.pop();
      used.delete(next);
    }
  };

  for (let start = 0; start < size; start++) walk(start, [start], new Set([start]));
  return [...cycles.values()];
}

function cycleProduct(cycle, weights) {
  let product = 1;
  for (let i = 0; i < cycle.length; i++) {
    product *= weights[cycle[i]][cycle[(i + 1) % cycle.length]];
  }
  return product;
}

export function calculateCumulants(model, maxOrder = 5) {
  const graph = buildAugmentedGraph(model);
  const edges = [];
  for (let i = 0; i < graph.size; i++) {
    for (let j = i + 1; j < graph.size; j++) {
      if (graph.weights[i][j] !== 0) edges.push(graph.weights[i][j]);
    }
  }

  const triangles = maxOrder >= 3 || maxOrder >= 5 ? enumerateCycles(graph, 3) : [];
  const fourCycles = maxOrder >= 4 ? enumerateCycles(graph, 4) : [];
  const fiveCycles = maxOrder >= 5 ? enumerateCycles(graph, 5) : [];
  const kappa = [0, 0, 0, 0, 0, 0];
  kappa[2] = edges.reduce((sum, k) => sum + k ** 2, 0);
  if (maxOrder >= 3) kappa[3] = -6 * triangles.reduce((sum, c) => sum + cycleProduct(c, graph.weights), 0);
  if (maxOrder >= 4) {
    kappa[4] = -2 * edges.reduce((sum, k) => sum + k ** 4, 0)
      + 24 * fourCycles.reduce((sum, c) => sum + cycleProduct(c, graph.weights), 0);
  }
  if (maxOrder >= 5) {
    kappa[5] = 40 * triangles.reduce((sum, c) => {
      const product = cycleProduct(c, graph.weights);
      const squareSum = c.reduce((inner, v, idx) => inner + graph.weights[v][c[(idx + 1) % c.length]] ** 2, 0);
      return sum + product * squareSum;
    }, 0) - 120 * fiveCycles.reduce((sum, c) => sum + cycleProduct(c, graph.weights), 0);
  }

  return {
    kappa,
    counts: {
      edges: edges.length,
      triangles: triangles.length,
      cycles4: fourCycles.length,
      cycles5: fiveCycles.length
    }
  };
}

export function approximateLnZ(model, beta, order) {
  const { kappa, counts } = calculateCumulants(model, order);
  let value = -beta * (model.offset || 0) + model.n * Math.log(2);
  for (let r = 2; r <= order; r++) {
    value += ((-beta) ** r / factorial(r)) * kappa[r];
  }
  return { value, kappa, counts };
}

export function exactEnergy(originalModel, spinsByOriginalIndex) {
  let energy = originalModel.offset || 0;
  for (let i = 0; i < originalModel.n; i++) energy -= originalModel.fields[i] * spinsByOriginalIndex[i];
  for (const [key, value] of originalModel.couplings.entries()) {
    const [i, j] = key.split(":").map(Number);
    energy -= value * spinsByOriginalIndex[i] * spinsByOriginalIndex[j];
  }
  return energy;
}

export function formatNumber(value) {
  if (!Number.isFinite(value)) return String(value);
  if (Math.abs(value) < 1e-10) return "0";
  return Number(value.toPrecision(8)).toString();
}

function factorial(n) {
  return n <= 1 ? 1 : n * factorial(n - 1);
}
