const EPS = 1e-12;

function edgeKey(a, b) {
  return a < b ? `${a}:${b}` : `${b}:${a}`;
}

function cloneModel(model) {
  return {
    n: model.n,
    labels: [...model.labels],
    fields: [...model.fields],
    couplings: new Map(model.couplings),
    offset: model.offset || 0,
    normalizationStats: { ...(model.normalizationStats || {}) }
  };
}

function createModel(n, couplings = [], fields = []) {
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
    const value = Number(c.J ?? c.value ?? c[2]);
    if (Math.abs(value) < EPS) {
      model.normalizationStats.zeroCouplingsRemoved++;
      continue;
    }
    if (i === j) {
      model.offset -= value;
      model.normalizationStats.selfCouplingsConverted++;
      continue;
    }
    const key = edgeKey(i, j);
    if (model.couplings.has(key)) model.normalizationStats.duplicateEdgesMerged++;
    model.couplings.set(key, (model.couplings.get(key) || 0) + value);
  }
  for (const [key, value] of [...model.couplings.entries()]) {
    if (Math.abs(value) < EPS) {
      model.couplings.delete(key);
      model.normalizationStats.zeroCouplingsRemoved++;
    }
  }
  return model;
}

function parseIsingText(text) {
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
      fields[Number(parts[0]) - 1] = Number(parts[1]);
    } else if (section === "COUPLINGS") {
      couplings.push({ i: Number(parts[0]) - 1, j: Number(parts[1]) - 1, J: Number(parts[2]) });
    }
  }
  if (!Number.isInteger(n) || n < 1) throw new Error("Input must include a positive N line.");
  return createModel(n, couplings, Array.from({ length: n }, (_, i) => fields[i] || 0));
}

function buildAugmentedGraph(model) {
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
  for (const seq of [path, [...path].reverse()]) {
    for (let shift = 0; shift < n; shift++) {
      const key = [...seq.slice(shift), ...seq.slice(0, shift)].join(",");
      if (best === null || key < best) best = key;
    }
  }
  return best;
}

function enumerateCycles(graph, length) {
  const cycles = new Map();
  const walk = (start, path, used) => {
    const last = path[path.length - 1];
    if (path.length === length) {
      if (graph.weights[last][start] !== 0) cycles.set(canonicalCycle(path), [...path]);
      return;
    }
    for (const next of graph.adjacency[last]) {
      if (next === start || used.has(next)) continue;
      used.add(next);
      path.push(next);
      walk(start, path, used);
      path.pop();
      used.delete(next);
    }
  };
  for (let start = 0; start < graph.size; start++) walk(start, [start], new Set([start]));
  return [...cycles.values()];
}

function cycleProduct(cycle, weights) {
  let product = 1;
  for (let i = 0; i < cycle.length; i++) product *= weights[cycle[i]][cycle[(i + 1) % cycle.length]];
  return product;
}

function calculateCumulants(model, maxOrder = 5) {
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
  return { kappa, counts: { edges: edges.length, triangles: triangles.length, cycles4: fourCycles.length, cycles5: fiveCycles.length } };
}

function approximateLnZ(model, beta, order) {
  const { kappa, counts } = calculateCumulants(model, order);
  let value = -beta * (model.offset || 0) + model.n * Math.log(2);
  for (let r = 2; r <= order; r++) value += ((-beta) ** r / factorial(r)) * kappa[r];
  return { value, kappa, counts };
}

function exactEnergy(originalModel, spins) {
  let energy = originalModel.offset || 0;
  for (let i = 0; i < originalModel.n; i++) energy -= originalModel.fields[i] * spins[i];
  for (const [key, value] of originalModel.couplings.entries()) {
    const [i, j] = key.split(":").map(Number);
    energy -= value * spins[i] * spins[j];
  }
  return energy;
}

function fixSpin(model, localSpinIndex, value) {
  const nextLabels = [];
  const nextFields = [];
  const nextCouplings = new Map();
  const localMap = new Map();
  const originalField = model.fields[localSpinIndex];
  let offset = model.offset - originalField * value;
  const fieldUpdates = [];
  for (let i = 0; i < model.n; i++) {
    if (i === localSpinIndex) continue;
    localMap.set(i, nextLabels.length);
    nextLabels.push(model.labels[i]);
    nextFields.push(model.fields[i]);
  }
  for (const [key, coupling] of model.couplings.entries()) {
    const [a, b] = key.split(":").map(Number);
    if (a === localSpinIndex || b === localSpinIndex) {
      const other = a === localSpinIndex ? b : a;
      if (localMap.has(other)) {
        const mapped = localMap.get(other);
        nextFields[mapped] += coupling * value;
        fieldUpdates.push({ spin: model.labels[other] + 1, delta: coupling * value });
      }
    } else {
      nextCouplings.set(edgeKey(localMap.get(a), localMap.get(b)), coupling);
    }
  }
  return {
    originalSpin: model.labels[localSpinIndex],
    value,
    constantChange: -originalField * value,
    remaining: nextLabels.length,
    fieldUpdates,
    model: {
      n: nextLabels.length,
      labels: nextLabels,
      fields: nextFields,
      couplings: nextCouplings,
      offset,
      normalizationStats: { ...(model.normalizationStats || {}) }
    }
  };
}

function connectedComponents(model) {
  const adj = adjacency(model);
  const seen = new Set();
  const components = [];
  for (let start = 0; start < model.n; start++) {
    if (seen.has(start)) continue;
    const queue = [start];
    const vertices = [];
    seen.add(start);
    for (let p = 0; p < queue.length; p++) {
      const v = queue[p];
      vertices.push(v);
      for (const next of adj[v]) {
        if (!seen.has(next)) {
          seen.add(next);
          queue.push(next);
        }
      }
    }
    const vertexSet = new Set(vertices);
    let edges = 0;
    for (const [key] of model.couplings.entries()) {
      const [a, b] = key.split(":").map(Number);
      if (vertexSet.has(a) && vertexSet.has(b)) edges++;
    }
    const gauge = bipartiteGauge(model, vertices, adj);
    components.push({ vertices, edges, bipartite: gauge.bipartite, unfrustratedGauge: gauge.unfrustratedGauge });
  }
  return components;
}

function preprocessState(model, spins, options = {}) {
  let current = cloneModel(model);
  const assignments = [...spins];
  const log = [];
  const summary = {
    originalSpins: current.n,
    originalCouplings: current.couplings.size,
    zeroCouplingsRemoved: options.countNormalization === false ? 0 : current.normalizationStats?.zeroCouplingsRemoved || 0,
    duplicateEdgesMerged: options.countNormalization === false ? 0 : current.normalizationStats?.duplicateEdgesMerged || 0,
    selfCouplingsConverted: options.countNormalization === false ? 0 : current.normalizationStats?.selfCouplingsConverted || 0,
    disconnectedComponents: 0,
    isolatedSpinsFixed: 0,
    dominantFieldSpinsFixed: 0,
    symmetrySpinsFixed: 0,
    degenerateSpinsFixed: 0,
    treeComponentsSolved: 0,
    bipartiteComponents: 0,
    unfrustratedGaugeComponents: 0,
    remainingSpins: current.n
  };
  let changed = true;
  while (changed) {
    changed = false;
    normalizeCurrent(current, summary);
    const components = connectedComponents(current);
    summary.disconnectedComponents = Math.max(summary.disconnectedComponents, components.length);
    summary.bipartiteComponents = Math.max(summary.bipartiteComponents, components.filter((c) => c.bipartite).length);
    summary.unfrustratedGaugeComponents = Math.max(summary.unfrustratedGaugeComponents, components.filter((c) => c.unfrustratedGauge).length);

    const isolated = current.fields.findIndex((_, i) => degree(current, i) === 0);
    if (isolated >= 0) {
      const h = current.fields[isolated];
      const value = h < -EPS ? -1 : 1;
      const reason = Math.abs(h) < EPS ? "DEGENERATE" : "FORCED";
      const result = fixSpin(current, isolated, value);
      assignments[result.originalSpin] = value;
      current = result.model;
      summary.isolatedSpinsFixed++;
      if (reason === "DEGENERATE") summary.degenerateSpinsFixed++;
      log.push(makePrecheckLog(result, reason, Math.abs(h) < EPS ? "isolated zero-field spin" : "isolated field spin", `h = ${formatNumber(h)}`));
      changed = true;
      continue;
    }

    const symmetricComponent = components.find((component) => component.vertices.length > 1 && component.vertices.every((i) => Math.abs(current.fields[i]) < EPS));
    if (symmetricComponent) {
      const local = [...symmetricComponent.vertices].sort((a, b) => current.labels[a] - current.labels[b])[0];
      const result = fixSpin(current, local, 1);
      assignments[result.originalSpin] = 1;
      current = result.model;
      summary.symmetrySpinsFixed++;
      log.push(makePrecheckLog(result, "SYMMETRY", "zero-field global spin symmetry", "lowest-numbered spin fixed to +1"));
      changed = true;
      continue;
    }

    const dominant = findDominantField(current);
    if (dominant) {
      const value = dominant.h > 0 ? 1 : -1;
      const result = fixSpin(current, dominant.local, value);
      assignments[result.originalSpin] = value;
      current = result.model;
      summary.dominantFieldSpinsFixed++;
      log.push(makePrecheckLog(result, "FORCED", "dominant local field", `|h| = ${formatNumber(Math.abs(dominant.h))}, sum |J| = ${formatNumber(dominant.radius)}`));
      changed = true;
    }
  }
  summary.remainingSpins = current.n;
  return { model: current, spins: assignments, log, summary };
}

function normalizeCurrent(model, summary) {
  for (const [key, value] of [...model.couplings.entries()]) {
    if (Math.abs(value) < EPS) {
      model.couplings.delete(key);
      summary.zeroCouplingsRemoved++;
    }
  }
}

function degree(model, local) {
  let count = 0;
  for (const key of model.couplings.keys()) {
    const [a, b] = key.split(":").map(Number);
    if (a === local || b === local) count++;
  }
  return count;
}

function findDominantField(model) {
  for (let i = 0; i < model.n; i++) {
    const h = model.fields[i];
    let radius = 0;
    for (const [key, value] of model.couplings.entries()) {
      const [a, b] = key.split(":").map(Number);
      if (a === i || b === i) radius += Math.abs(value);
    }
    if (Math.abs(h) > radius + EPS) return { local: i, h, radius };
  }
  return null;
}

function adjacency(model) {
  const adj = Array.from({ length: model.n }, () => []);
  for (const [key, value] of model.couplings.entries()) {
    if (Math.abs(value) < EPS) continue;
    const [a, b] = key.split(":").map(Number);
    adj[a].push(b);
    adj[b].push(a);
  }
  return adj;
}

function bipartiteGauge(model, vertices, adj) {
  const color = new Map();
  const sign = new Map();
  let bipartite = true;
  let unfrustratedGauge = true;
  const vertexSet = new Set(vertices);
  for (const start of vertices) {
    if (!color.has(start)) {
      color.set(start, 1);
      sign.set(start, 1);
    }
    const queue = [start];
    for (let p = 0; p < queue.length; p++) {
      const v = queue[p];
      for (const n of adj[v]) {
        if (!vertexSet.has(n)) continue;
        if (!color.has(n)) {
          color.set(n, -color.get(v));
          queue.push(n);
        } else if (color.get(n) === color.get(v)) {
          bipartite = false;
        }
        const coupling = model.couplings.get(edgeKey(v, n));
        const wanted = Math.sign(coupling || 1) * sign.get(v);
        if (!sign.has(n)) sign.set(n, wanted);
        else if (sign.get(n) !== wanted) unfrustratedGauge = false;
      }
    }
  }
  return { bipartite, unfrustratedGauge };
}

function makePrecheckLog(result, reason, rule, detail) {
  return {
    type: "precheck",
    spin: result.originalSpin + 1,
    value: result.value,
    reason,
    rule,
    detail,
    constantChange: result.constantChange,
    remaining: result.remaining,
    fieldUpdates: result.fieldUpdates
  };
}

function conditionModel(model, localSpinIndex, value) {
  return fixSpin(model, localSpinIndex, value).model;
}

function scoreSpin(model, localSpinIndex, beta, order) {
  const plus = approximateLnZ(conditionModel(model, localSpinIndex, 1), beta, order);
  const minus = approximateLnZ(conditionModel(model, localSpinIndex, -1), beta, order);
  return {
    localSpinIndex,
    originalSpin: model.labels[localSpinIndex],
    Lplus: plus.value,
    Lminus: minus.value,
    importance: plus.value - minus.value,
    preferred: plus.value - minus.value >= 0 ? 1 : -1,
    alternate: plus.value - minus.value >= 0 ? -1 : 1
  };
}

function solveTaylor(originalModel, options = {}) {
  const beta = Number(options.beta ?? 0.5);
  const order = 5;
  const mode = options.mode || "max";
  const maxBacktracks = Math.max(0, Number(options.maxBacktracks ?? 0));
  const started = performance.now();
  let bestEnergy = Infinity;
  let bestSpins = null;
  let bestPrecheckSummary = {};
  let backtracks = 0;
  let decisions = 0;
  let firstComplete = null;
  const decisionLog = [];
  const initial = {
    model: cloneModel(originalModel),
    spins: Array(originalModel.n).fill(null),
    depth: 0,
    trace: [],
    deviations: 0
  };
  const stack = [initial];
  while (stack.length) {
    const node = stack.pop();
    const prechecked = preprocessState(node.model, node.spins, { countNormalization: !node.precheckSummary });
    node.model = prechecked.model;
    node.spins = prechecked.spins;
    node.trace = [...node.trace, ...prechecked.log];
    node.precheckSummary = mergeSummaries(node.precheckSummary, prechecked.summary);
    if (node.model.n === 0) {
      const energy = exactEnergy(originalModel, node.spins);
      if (!firstComplete) firstComplete = { energy, spins: [...node.spins], trace: node.trace };
      if (energy < bestEnergy - 1e-10) {
        bestEnergy = energy;
        bestSpins = [...node.spins];
        bestPrecheckSummary = node.precheckSummary || {};
        decisionLog.splice(0, decisionLog.length, ...node.trace);
      }
      continue;
    }
    const scores = [];
    if (mode === "sequential") {
      scores.push(scoreSpin(node.model, 0, beta, order));
    } else {
      for (let i = 0; i < node.model.n; i++) scores.push(scoreSpin(node.model, i, beta, order));
      scores.sort((a, b) => Math.abs(b.importance) - Math.abs(a.importance));
    }
    const chosen = scores[0];
    decisions++;
    for (const branch of [chosen.alternate, chosen.preferred]) {
      const isAlternate = branch === chosen.alternate;
      const deviations = node.deviations + (isAlternate ? 1 : 0);
      if (isAlternate && deviations > maxBacktracks) continue;
      if (isAlternate) backtracks++;
      const spins = [...node.spins];
      spins[chosen.originalSpin] = branch;
      const fixed = fixSpin(node.model, chosen.localSpinIndex, branch);
      const step = {
        type: "taylor",
        step: node.depth + 1,
        spin: chosen.originalSpin + 1,
        Lplus: chosen.Lplus,
        Lminus: chosen.Lminus,
        importance: chosen.importance,
        chosen: branch,
        preferred: chosen.preferred,
        isAlternate,
        candidateScores: mode === "max" ? scores : [chosen],
        backtracksSoFar: backtracks,
        fieldUpdates: fixed.fieldUpdates,
        constantChange: fixed.constantChange
      };
      stack.push({ model: fixed.model, spins, depth: node.depth + 1, trace: [...node.trace, step], deviations, precheckSummary: node.precheckSummary });
    }
  }
  const fallback = firstComplete || { energy: Infinity, spins: [], trace: [] };
  return {
    energy: Number.isFinite(bestEnergy) ? bestEnergy : fallback.energy,
    spins: bestSpins || fallback.spins,
    bestEnergy: Number.isFinite(bestEnergy) ? bestEnergy : fallback.energy,
    bestSpins: bestSpins || fallback.spins,
    decisions,
    backtracks,
    runtimeMs: performance.now() - started,
    decisionLog: decisionLog.length ? decisionLog : fallback.trace,
    precheckSummary: bestPrecheckSummary
  };
}

function mergeSummaries(a = {}, b = {}) {
  const next = { ...a };
  for (const [key, value] of Object.entries(b)) {
    if (typeof value !== "number") {
      next[key] = value;
    } else if (key === "remainingSpins") {
      next[key] = value;
    } else if (key === "originalSpins" || key === "originalCouplings") {
      next[key] ??= value;
    } else {
      next[key] = (next[key] || 0) + value;
    }
  }
  return next;
}

function factorial(n) {
  return n <= 1 ? 1 : n * factorial(n - 1);
}

function buildSpinGlass50Text() {
  let seed = 501337;
  const rand = () => {
    seed = (1664525 * seed + 1013904223) >>> 0;
    return seed / 2 ** 32;
  };
  const lines = ["N 50", "", "FIELDS"];
  for (let i = 1; i <= 50; i++) lines.push(`${i} ${Math.round((rand() * 1.2 - 0.6) * 100) / 100}`);
  lines.push("", "COUPLINGS");
  const used = new Set();
  for (let i = 1; i <= 49; i++) {
    const j = i + 1;
    const J = Math.round((rand() * 2 - 1) * 100) / 100 || 0.25;
    used.add(`${i}:${j}`);
    lines.push(`${i} ${j} ${J}`);
  }
  while (used.size < 120) {
    const i = 1 + Math.floor(rand() * 50);
    const j = 1 + Math.floor(rand() * 50);
    if (i === j) continue;
    const a = Math.min(i, j);
    const b = Math.max(i, j);
    const key = `${a}:${b}`;
    if (used.has(key)) continue;
    used.add(key);
    const J = Math.round((rand() * 2.4 - 1.2) * 100) / 100 || -0.35;
    lines.push(`${a} ${b} ${J}`);
  }
  return lines.join("\n");
}

const examples = {
  spinGlass50: {
    name: "50-node spin glass",
    benchmark: "No built-in certificate; use for visual/search behavior.",
    text: buildSpinGlass50Text()
  },
  ferroRing8: {
    name: "8-spin ferromagnetic ring",
    benchmark: "Known ground energy -8; all spins equal.",
    text: "N 8\n\nFIELDS\n1 0\n2 0\n3 0\n4 0\n5 0\n6 0\n7 0\n8 0\n\nCOUPLINGS\n1 2 1\n2 3 1\n3 4 1\n4 5 1\n5 6 1\n6 7 1\n7 8 1\n8 1 1"
  },
  antiferroRing8: {
    name: "8-spin antiferro ring",
    benchmark: "Known ground energy -8; alternating spins.",
    text: "N 8\n\nFIELDS\n1 0\n2 0\n3 0\n4 0\n5 0\n6 0\n7 0\n8 0\n\nCOUPLINGS\n1 2 -1\n2 3 -1\n3 4 -1\n4 5 -1\n5 6 -1\n6 7 -1\n7 8 -1\n8 1 -1"
  },
  frustratedRing9: {
    name: "9-spin frustrated AF ring",
    benchmark: "Known ground energy -7; one unsatisfied bond is unavoidable.",
    text: "N 9\n\nFIELDS\n1 0\n2 0\n3 0\n4 0\n5 0\n6 0\n7 0\n8 0\n9 0\n\nCOUPLINGS\n1 2 -1\n2 3 -1\n3 4 -1\n4 5 -1\n5 6 -1\n6 7 -1\n7 8 -1\n8 9 -1\n9 1 -1"
  },
  plantedField12: {
    name: "12-spin planted fields",
    benchmark: "Known planted assignment: + + - + - - + - + + - -.",
    text: "N 12\n\nFIELDS\n1 2.4\n2 2.1\n3 -2.6\n4 2.2\n5 -2.5\n6 -2.3\n7 2.7\n8 -2.4\n9 2.2\n10 2.5\n11 -2.1\n12 -2.8\n\nCOUPLINGS\n1 2 0.3\n2 3 -0.2\n3 4 0.25\n4 5 -0.3\n5 6 0.2\n6 7 -0.25\n7 8 0.35\n8 9 -0.2\n9 10 0.3\n10 11 -0.25\n11 12 0.2\n12 1 -0.3\n1 7 0.15\n3 9 -0.15\n5 11 0.1"
  }
};

const $ = (id) => document.getElementById(id);
const visualState = {
  model: null,
  result: null,
  nodes: [],
  edges: [],
  assignments: new Map(),
  reasons: new Map(),
  importance: new Map(),
  selected: null,
  step: 0,
  timer: null
};

function options() {
  return {
    beta: Number($("beta").value),
    order: 5,
    mode: $("mode").value,
    maxBacktracks: Number($("backtracks").value),
    visualizationOn: $("visualizationOn").checked,
    maxRenderedSpins: Number($("maxRenderedSpins").value)
  };
}

function write(text) {
  $("output").textContent = text;
}

function setSummary(metrics) {
  $("summaryGrid").innerHTML = metrics.map(([label, value]) => (
    `<div class="metric"><span>${escapeHtml(label)}</span><strong>${escapeHtml(String(value))}</strong></div>`
  )).join("");
}

function currentModel() {
  return parseIsingText($("isingInput").value);
}

function modelStats(model) {
  return {
    fields: model.fields.filter((h) => Math.abs(h) > EPS).length,
    couplings: model.couplings.size
  };
}

function runSolver() {
  try {
    write("Running fifth-order Taylor/cumulant solver. For 50 nodes this can take a while in the browser.");
    setTimeout(() => {
      try {
        const model = currentModel();
        const opts = options();
        const cumulants = calculateCumulants(model, 5);
        const lnz = approximateLnZ(model, opts.beta, 5);
        const result = solveTaylor(model, opts);
        if (opts.visualizationOn) prepareVisualization(model, result);
        const stats = modelStats(model);
        setSummary([
          ["spins", model.n],
          ["couplings", stats.couplings],
          ["fields", stats.fields],
          ["energy", formatNumber(result.energy)],
          ["backtracks", result.backtracks],
          ["runtime", `${formatNumber(result.runtimeMs)} ms`]
        ]);
        const lines = [
          "Input summary",
          `N: ${model.n}`,
          `Couplings: ${stats.couplings}`,
          `Nonzero fields: ${stats.fields}`,
          `Beta: ${opts.beta}`,
          "Taylor order: 5",
          "",
          "Cumulants",
          ...[1, 2, 3, 4, 5].map((r) => `kappa${r}: ${formatNumber(cumulants.kappa[r])}`),
          `Cycles: triangles=${cumulants.counts.triangles}, C4=${cumulants.counts.cycles4}, C5=${cumulants.counts.cycles5}`,
          `ln Z approximation: ${formatNumber(lnz.value)}`,
          "",
          "Search result",
          `Final spins: ${formatSpins(result.spins)}`,
          `Exact final energy: ${formatNumber(result.energy)}`,
          `Best energy found: ${formatNumber(result.bestEnergy)}`,
          `Decisions: ${result.decisions}`,
          `Backtracks: ${result.backtracks}`,
          `Runtime: ${formatNumber(result.runtimeMs)} ms`,
          "",
          "Precheck summary",
          ...formatPrecheckSummary(result.precheckSummary),
          "",
          "Decision trace",
          ...formatDecisionLog(result.decisionLog)
        ];
        write(lines.join("\n"));
      } catch (error) {
        write(`Error: ${error.message}`);
      }
    }, 20);
  } catch (error) {
    write(`Error: ${error.message}`);
  }
}

function runBetaBenchmark() {
  try {
    const model = currentModel();
    const opts = options();
    const betas = $("betaList").value.split(/\s+/).map(Number).filter(Number.isFinite);
    const rows = ["Beta | Energy found | Backtracks | Runtime ms", "---- | ------------ | ---------- | ----------"];
    for (const beta of betas) {
      const result = solveTaylor(model, { ...opts, beta });
      rows.push(`${beta} | ${formatNumber(result.energy)} | ${result.backtracks} | ${formatNumber(result.runtimeMs)}`);
    }
    setSummary([["benchmark", "beta sweep"], ["order", 5], ["mode", opts.mode], ["cases", betas.length]]);
    write(rows.join("\n"));
  } catch (error) {
    write(`Error: ${error.message}`);
  }
}

function formatDecisionLog(log) {
  if (!log.length) return ["No decisions recorded."];
  return log.flatMap((step) => {
    if (step.type === "precheck") {
      const updates = step.fieldUpdates.length ? step.fieldUpdates.map((u) => `  h${u.spin} += ${formatNumber(u.delta)}`).join("\n") : "  none";
      return [
        "Precheck reduction",
        `Spin: ${step.spin}`,
        `Assigned value: ${step.value > 0 ? "+1" : "-1"}`,
        `Reason: ${step.reason}`,
        `Rule: ${step.rule}`,
        `Data: ${step.detail}`,
        `Constant-energy change: ${formatNumber(step.constantChange)}`,
        `Field updates:\n${updates}`,
        `Remaining spins: ${step.remaining}`,
        ""
      ];
    }
    const candidates = step.candidateScores.map((score) => `  s${score.originalSpin + 1}: I=${formatNumber(score.importance)} preferred=${score.preferred > 0 ? "+1" : "-1"}`).join("\n");
    return [
      `Step ${step.step}`,
      `Spin tested: ${step.spin}`,
      `lnZ(s${step.spin} = +1): ${formatNumber(step.Lplus)}`,
      `lnZ(s${step.spin} = -1): ${formatNumber(step.Lminus)}`,
      `Importance: ${formatNumber(step.importance)}`,
      `Chosen: s${step.spin} = ${step.chosen > 0 ? "+1" : "-1"}${step.isAlternate ? " (alternate branch)" : ""}`,
      `Backtracks so far: ${step.backtracksSoFar}`,
      candidates,
      ""
    ];
  });
}

function formatPrecheckSummary(summary = {}) {
  return [
    ["Original spins", summary.originalSpins],
    ["Original couplings", summary.originalCouplings],
    ["Zero couplings removed", summary.zeroCouplingsRemoved],
    ["Duplicate edges merged", summary.duplicateEdgesMerged],
    ["Self-couplings converted", summary.selfCouplingsConverted],
    ["Disconnected components observed", summary.disconnectedComponents],
    ["Isolated spins fixed", summary.isolatedSpinsFixed],
    ["Dominant-field spins fixed", summary.dominantFieldSpinsFixed],
    ["Symmetry spins fixed", summary.symmetrySpinsFixed],
    ["Degenerate spins fixed", summary.degenerateSpinsFixed],
    ["Bipartite components detected", summary.bipartiteComponents],
    ["Unfrustrated gauge components detected", summary.unfrustratedGaugeComponents],
    ["Remaining unresolved spins", summary.remainingSpins]
  ].map(([label, value]) => `${label}: ${value ?? 0}`);
}

function prepareVisualization(model, result) {
  visualState.model = model;
  visualState.result = result;
  visualState.assignments = new Map();
  visualState.reasons = new Map();
  visualState.importance = new Map();
  visualState.selected = null;
  visualState.step = 0;
  const cap = Number($("maxRenderedSpins").value || 150);
  const count = Math.min(model.n, cap);
  const nodes = Array.from({ length: count }, (_, i) => {
    const angle = (i / count) * Math.PI * 2;
    const radius = 285 + 76 * Math.sin(i * 1.7);
    return { index: i, x: 450 + Math.cos(angle) * radius, y: 380 + Math.sin(angle) * radius, vx: 0, vy: 0 };
  });
  const rendered = new Set(nodes.map((n) => n.index));
  visualState.edges = [...model.couplings.entries()].map(([key, J]) => {
    const [a, b] = key.split(":").map(Number);
    return { a, b, J };
  }).filter((e) => rendered.has(e.a) && rendered.has(e.b));
  visualState.nodes = nodes;
  settleLayout();
  drawVisualization();
}

function settleLayout() {
  const nodes = visualState.nodes;
  const edges = visualState.edges;
  for (let tick = 0; tick < 140; tick++) {
    for (const node of nodes) {
      node.vx += (450 - node.x) * 0.0007;
      node.vy += (380 - node.y) * 0.0007;
    }
    for (let i = 0; i < nodes.length; i++) {
      for (let j = i + 1; j < nodes.length; j++) {
        const a = nodes[i];
        const b = nodes[j];
        const dx = a.x - b.x;
        const dy = a.y - b.y;
        const force = 56 / Math.max(64, dx * dx + dy * dy);
        a.vx += dx * force;
        a.vy += dy * force;
        b.vx -= dx * force;
        b.vy -= dy * force;
      }
    }
    for (const edge of edges) {
      const a = nodes[edge.a];
      const b = nodes[edge.b];
      if (!a || !b) continue;
      const dx = b.x - a.x;
      const dy = b.y - a.y;
      const dist = Math.max(1, Math.hypot(dx, dy));
      const force = (dist - 122) * 0.0028;
      a.vx += dx / dist * force;
      a.vy += dy / dist * force;
      b.vx -= dx / dist * force;
      b.vy -= dy / dist * force;
    }
    for (const node of nodes) {
      node.vx *= 0.86;
      node.vy *= 0.86;
      node.x = Math.min(870, Math.max(30, node.x + node.vx));
      node.y = Math.min(730, Math.max(30, node.y + node.vy));
    }
  }
}

function drawVisualization() {
  const canvas = $("spinCanvas");
  if (!canvas || !visualState.model) return;
  const ctx = canvas.getContext("2d");
  const rect = canvas.getBoundingClientRect();
  const scale = Math.max(1, window.devicePixelRatio || 1);
  canvas.width = Math.max(640, Math.floor(rect.width * scale));
  canvas.height = Math.max(520, Math.floor(760 * scale));
  ctx.setTransform(canvas.width / 900, 0, 0, canvas.height / 760, 0, 0);
  ctx.clearRect(0, 0, 900, 760);
  const maxJ = Math.max(0.01, ...visualState.edges.map((e) => Math.abs(e.J)));
  for (const edge of visualState.edges) {
    const a = visualState.nodes[edge.a];
    const b = visualState.nodes[edge.b];
    if (!a || !b) continue;
    ctx.globalAlpha = visualState.assignments.has(edge.a) || visualState.assignments.has(edge.b) ? 0.16 : 0.48;
    ctx.strokeStyle = edge.J > 0 ? "#36d3c2" : "#f4b44f";
    ctx.lineWidth = 0.55 + 3.2 * Math.abs(edge.J) / maxJ;
    ctx.setLineDash(edge.J < 0 ? [7, 5] : []);
    ctx.beginPath();
    ctx.moveTo(a.x, a.y);
    ctx.lineTo(b.x, b.y);
    ctx.stroke();
  }
  ctx.setLineDash([]);
  ctx.globalAlpha = 1;
  for (const node of visualState.nodes) {
    const value = visualState.assignments.get(node.index);
    const reason = visualState.reasons.get(node.index) || "UNRESOLVED";
    const h = visualState.model.fields[node.index] || 0;
    const fieldRadius = 13 + Math.min(8, Math.abs(h) * 9);
    ctx.globalAlpha = 0.24;
    ctx.fillStyle = h >= 0 ? "#36d3c2" : "#f4b44f";
    ctx.beginPath();
    ctx.arc(node.x, node.y, fieldRadius, 0, Math.PI * 2);
    ctx.fill();
    ctx.globalAlpha = 1;
    ctx.fillStyle = nodeColor(value, reason);
    ctx.strokeStyle = visualState.selected === node.index ? "#edf6f5" : "#283443";
    ctx.lineWidth = visualState.selected === node.index ? 3 : 1.5;
    ctx.beginPath();
    ctx.arc(node.x, node.y, 12, 0, Math.PI * 2);
    ctx.fill();
    ctx.stroke();
    if (value !== undefined) {
      ctx.fillStyle = "#041c19";
      ctx.font = "bold 14px system-ui";
      ctx.textAlign = "center";
      ctx.textBaseline = "middle";
      ctx.fillText(value > 0 ? "+" : "-", node.x, node.y);
      ctx.textBaseline = "alphabetic";
    }
    if (visualState.selected === node.index || node.index % 5 === 0) {
      ctx.fillStyle = "#edf6f5";
      ctx.font = visualState.selected === node.index ? "12px system-ui" : "10px system-ui";
      ctx.textAlign = "center";
      ctx.fillText(String(node.index + 1), node.x, node.y - 19);
    }
  }
  updateVisualInfo();
}

function nodeColor(value, reason) {
  if (reason === "BACKTRACKED") return "#ff7d7d";
  if (reason === "FORCED") return "#77e0a3";
  if (reason === "SYMMETRY") return "#b893ff";
  if (reason === "DEGENERATE") return "#98a6b3";
  if (reason === "TAYLOR_DECISION") return value > 0 ? "#36d3c2" : "#f4b44f";
  return "#17212c";
}

function runVisualAnimation() {
  if (!visualState.result) runSolver();
  pauseVisualAnimation();
  visualState.timer = setInterval(() => stepVisualization(1), 650);
}

function pauseVisualAnimation() {
  if (visualState.timer) clearInterval(visualState.timer);
  visualState.timer = null;
}

function resetVisualization() {
  pauseVisualAnimation();
  prepareVisualization(currentModel(), visualState.result);
}

function stepVisualization(count) {
  if (!visualState.result) return;
  const log = visualState.result.decisionLog || [];
  for (let c = 0; c < count && visualState.step < log.length; c++) {
    const entry = log[visualState.step++];
    const spin = entry.spin - 1;
    if (spin >= visualState.nodes.length) continue;
    visualState.selected = spin;
    visualState.assignments.set(spin, entry.type === "precheck" ? entry.value : entry.chosen);
    visualState.reasons.set(spin, entry.type === "precheck" ? entry.reason : (entry.isAlternate ? "BACKTRACKED" : "TAYLOR_DECISION"));
    if (entry.importance !== undefined) visualState.importance.set(spin, entry.importance);
  }
  if (visualState.step >= log.length) pauseVisualAnimation();
  drawVisualization();
}

function showFinalAssignment() {
  if (!visualState.result) {
    runSolver();
    return;
  }
  const log = visualState.result.decisionLog || [];
  while (visualState.step < log.length) {
    const entry = log[visualState.step++];
    const spin = entry.spin - 1;
    if (spin >= visualState.nodes.length) continue;
    visualState.selected = spin;
    visualState.assignments.set(spin, entry.type === "precheck" ? entry.value : entry.chosen);
    visualState.reasons.set(spin, entry.type === "precheck" ? entry.reason : (entry.isAlternate ? "BACKTRACKED" : "TAYLOR_DECISION"));
    if (entry.importance !== undefined) visualState.importance.set(spin, entry.importance);
  }
  for (let i = 0; i < visualState.result.spins.length && i < visualState.nodes.length; i++) {
    if (!visualState.assignments.has(i)) {
      visualState.assignments.set(i, visualState.result.spins[i]);
      visualState.reasons.set(i, "TAYLOR_DECISION");
    }
  }
  drawVisualization();
}

function selectCanvasSpin(event) {
  const canvas = $("spinCanvas");
  const rect = canvas.getBoundingClientRect();
  const x = (event.clientX - rect.left) * 900 / rect.width;
  const y = (event.clientY - rect.top) * 620 / rect.height;
  let best = null;
  let bestDist = Infinity;
  for (const node of visualState.nodes) {
    const dist = Math.hypot(node.x - x, node.y - y);
    if (dist < bestDist) {
      best = node;
      bestDist = dist;
    }
  }
  if (best && bestDist < 24) {
    visualState.selected = best.index;
    drawVisualization();
  }
}

function updateVisualInfo() {
  const model = visualState.model;
  const box = $("visualInfo");
  if (!model || !box) return;
  if (model.n > visualState.nodes.length) {
    box.textContent = `Model contains ${model.n} spins. Visualization capped at ${visualState.nodes.length} spins. Solver is still using the full model.`;
    return;
  }
  if (visualState.selected === null) {
    box.textContent = `Showing ${visualState.nodes.length} spins and ${visualState.edges.length} couplings. Assigned ${visualState.assignments.size}/${visualState.nodes.length}. Select a spin to inspect it.`;
    return;
  }
  const i = visualState.selected;
  const degree = visualState.edges.filter((e) => e.a === i || e.b === i).length;
  const value = visualState.assignments.has(i) ? (visualState.assignments.get(i) > 0 ? "+1" : "-1") : "unresolved";
  const reason = visualState.reasons.get(i) || "UNRESOLVED";
  const importance = visualState.importance.has(i) ? formatNumber(visualState.importance.get(i)) : "n/a";
  box.textContent = `Spin ${i + 1} | value ${value} | h = ${formatNumber(model.fields[i] || 0)} | degree ${degree} | reason ${reason} | I = ${importance}`;
}

function formatNumber(value) {
  if (!Number.isFinite(value)) return String(value);
  if (Math.abs(value) < 1e-10) return "0";
  return Number(value.toPrecision(8)).toString();
}

function formatSpins(spins) {
  return `[${spins.map((s) => s > 0 ? "+1" : "-1").join(" ")}]`;
}

function escapeHtml(value) {
  return value.replace(/[&<>"']/g, (ch) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[ch]));
}

function init() {
  $("startupHelp").style.display = "none";
  $("startupHelp").textContent = "JavaScript loaded. Buttons are active.";
  $("startupHelp").style.display = "block";
  for (const [key, example] of Object.entries(examples)) {
    const button = document.createElement("button");
    button.type = "button";
    button.className = "miniButton";
    button.innerHTML = `<strong>${escapeHtml(example.name)}</strong><span>${escapeHtml(example.benchmark)}</span>`;
    button.addEventListener("click", () => {
      $("isingInput").value = example.text;
      prepareVisualization(currentModel(), null);
      setSummary([["example", example.name], ["visual", "ready"]]);
      write(`${example.name} loaded.\nBenchmark note: ${example.benchmark}`);
    });
    $("exampleButtons").append(button);
  }
  $("isingInput").value = examples.spinGlass50.text;
  document.querySelectorAll(".tab").forEach((tab) => {
    tab.addEventListener("click", () => {
      document.querySelectorAll(".tab, .panel").forEach((node) => node.classList.remove("active"));
      tab.classList.add("active");
      $(tab.dataset.panel).classList.add("active");
      if (tab.dataset.panel === "visualizer") drawVisualization();
    });
  });
  $("previewGraph").addEventListener("click", () => {
    prepareVisualization(currentModel(), null);
    document.querySelectorAll(".tab, .panel").forEach((node) => node.classList.remove("active"));
    document.querySelector('[data-panel="visualizer"]').classList.add("active");
    $("visualizer").classList.add("active");
    setSummary([["visual", "ready"], ["spins", currentModel().n]]);
    write("Graph preview ready. This does not run the solver.");
  });
  $("runSolver").addEventListener("click", runSolver);
  $("runBetaBenchmark").addEventListener("click", runBetaBenchmark);
  $("clearOutput").addEventListener("click", () => {
    setSummary([]);
    write("Ready.");
  });
  $("inputFile").addEventListener("change", async (event) => {
    const file = event.target.files[0];
    if (file) {
      $("isingInput").value = await file.text();
      prepareVisualization(currentModel(), null);
    }
  });
  $("visualRun").addEventListener("click", runVisualAnimation);
  $("visualPause").addEventListener("click", pauseVisualAnimation);
  $("visualNext").addEventListener("click", () => stepVisualization(1));
  $("visualFinal").addEventListener("click", showFinalAssignment);
  $("visualReset").addEventListener("click", resetVisualization);
  $("spinCanvas").addEventListener("click", selectCanvasSpin);
  window.addEventListener("resize", drawVisualization);
  prepareVisualization(currentModel(), null);
  setSummary([["spins", 50], ["example", "50-node spin glass"], ["visual", "ready"]]);
  write("Ready. Preview shows unresolved spins. Press Run Taylor solver to compute assignments, then use Next or Final in the Visualizer.");
}

init();
