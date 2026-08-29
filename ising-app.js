const EPS = 1e-12;

function edgeKey(a, b) {
  return a < b ? `${a}:${b}` : `${b}:${a}`;
}

function hydrateModelMetadata(model) {
  const edges = [];
  const adjacency = Array.from({ length: model.n }, () => []);
  const degrees = Array(model.n).fill(0);
  const radii = Array(model.n).fill(0);
  for (const [key, value] of model.couplings.entries()) {
    if (Math.abs(value) < EPS) continue;
    const [i, j] = key.split(":").map(Number);
    const edge = { i, j, J: value };
    edges.push(edge);
    adjacency[i].push({ to: j, J: value });
    adjacency[j].push({ to: i, J: value });
    degrees[i]++;
    degrees[j]++;
    radii[i] += Math.abs(value);
    radii[j] += Math.abs(value);
  }
  model.edges = edges;
  model.adjacency = adjacency;
  model.degrees = degrees;
  model.radii = radii;
  return model;
}

function cloneModel(model) {
  return hydrateModelMetadata({
    n: model.n,
    labels: [...model.labels],
    fields: [...model.fields],
    couplings: new Map(model.couplings),
    offset: model.offset || 0,
    normalizationStats: { ...(model.normalizationStats || {}) }
  });
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
  return hydrateModelMetadata(model);
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
  for (const edge of model.edges || edgesFromCouplings(model)) {
    addEdge(edge.i + 1, edge.j + 1, edge.J);
  }
  adjacency.forEach((neighbors) => neighbors.sort((a, b) => a - b));
  return { size, weights, adjacency };
}

function edgesFromCouplings(model) {
  return [...model.couplings.entries()].map(([key, value]) => {
    const [i, j] = key.split(":").map(Number);
    return { i, j, J: value };
  });
}

function calculateCumulants(model, maxOrder = 5) {
  const graph = buildAugmentedGraph(model);
  const sums = fixedCycleSums(graph, maxOrder);
  const kappa = [0, 0, 0, 0, 0, 0];
  kappa[2] = sums.edgeSquareSum;
  if (maxOrder >= 3) kappa[3] = -6 * sums.triangleProductSum;
  if (maxOrder >= 4) kappa[4] = -2 * sums.edgeFourthSum + 24 * sums.fourCycleProductSum;
  if (maxOrder >= 5) kappa[5] = 40 * sums.triangleSquareProductSum - 120 * sums.fiveCycleProductSum;
  return { kappa, counts: { edges: sums.edges, triangles: sums.triangles, cycles4: sums.cycles4, cycles5: sums.cycles5 } };
}

function fixedCycleSums(graph, maxOrder) {
  const weights = graph.weights;
  const size = graph.size;
  const sums = {
    edges: 0,
    triangles: 0,
    cycles4: 0,
    cycles5: 0,
    edgeSquareSum: 0,
    edgeFourthSum: 0,
    triangleProductSum: 0,
    triangleSquareProductSum: 0,
    fourCycleProductSum: 0,
    fiveCycleProductSum: 0
  };
  for (let i = 0; i < graph.size; i++) {
    for (let j = i + 1; j < graph.size; j++) {
      const value = weights[i][j];
      if (value === 0) continue;
      sums.edges++;
      sums.edgeSquareSum += value ** 2;
      if (maxOrder >= 4) sums.edgeFourthSum += value ** 4;
    }
  }
  const possibleEdges = size * (size - 1) / 2;
  if (size > 12 && possibleEdges > 0 && sums.edges / possibleEdges < 0.35) {
    accumulateSparseCycleSums(graph, maxOrder, sums);
    return sums;
  }
  if (maxOrder >= 3 || maxOrder >= 5) {
    for (let i = 0; i < size - 2; i++) {
      for (let j = i + 1; j < size - 1; j++) {
        const ij = weights[i][j];
        if (ij === 0) continue;
        for (let k = j + 1; k < size; k++) {
          const jk = weights[j][k];
          const ki = weights[k][i];
          if (jk === 0 || ki === 0) continue;
          const product = ij * jk * ki;
          sums.triangles++;
          sums.triangleProductSum += product;
          if (maxOrder >= 5) sums.triangleSquareProductSum += product * (ij ** 2 + jk ** 2 + ki ** 2);
        }
      }
    }
  }
  if (maxOrder >= 4) {
    for (let a = 0; a < size - 3; a++) {
      for (let b = a + 1; b < size - 2; b++) {
        for (let c = b + 1; c < size - 1; c++) {
          for (let d = c + 1; d < size; d++) {
            sums.fourCycleProductSum += addFourCycle(weights, a, b, c, d, sums);
            sums.fourCycleProductSum += addFourCycle(weights, a, b, d, c, sums);
            sums.fourCycleProductSum += addFourCycle(weights, a, c, b, d, sums);
          }
        }
      }
    }
  }
  if (maxOrder >= 5) {
    const orders = [
      [1, 2, 3, 4], [1, 2, 4, 3], [1, 3, 2, 4], [1, 3, 4, 2],
      [1, 4, 2, 3], [1, 4, 3, 2], [2, 1, 3, 4], [2, 1, 4, 3],
      [2, 3, 1, 4], [2, 4, 1, 3], [3, 1, 2, 4], [3, 2, 1, 4]
    ];
    for (let a = 0; a < size - 4; a++) {
      for (let b = a + 1; b < size - 3; b++) {
        for (let c = b + 1; c < size - 2; c++) {
          for (let d = c + 1; d < size - 1; d++) {
            for (let e = d + 1; e < size; e++) {
              const vertices = [a, b, c, d, e];
              for (const order of orders) {
                sums.fiveCycleProductSum += addCycle(weights, [
                  a,
                  vertices[order[0]],
                  vertices[order[1]],
                  vertices[order[2]],
                  vertices[order[3]]
                ], sums, "cycles5");
              }
            }
          }
        }
      }
    }
  }
  return sums;
}

function accumulateSparseCycleSums(graph, maxOrder, sums) {
  const weights = graph.weights;
  const adjacency = graph.adjacency;
  const used = Array(graph.size).fill(false);
  const path = [];
  const lengths = [];
  if (maxOrder >= 3 || maxOrder >= 5) lengths.push(3);
  if (maxOrder >= 4) lengths.push(4);
  if (maxOrder >= 5) lengths.push(5);
  for (const length of lengths) {
    for (let start = 0; start < graph.size; start++) {
      path.length = 0;
      path.push(start);
      used[start] = true;
      walkSparseCycle(start, start, length, 1, 0);
      used[start] = false;
    }
  }
  function walkSparseCycle(start, last, length, product, squareSum) {
    if (path.length === length) {
      const close = weights[last][start];
      if (close === 0 || path[1] > last) return;
      const cycleProduct = product * close;
      if (length === 3) {
        sums.triangles++;
        sums.triangleProductSum += cycleProduct;
        if (maxOrder >= 5) sums.triangleSquareProductSum += cycleProduct * (squareSum + close ** 2);
      } else if (length === 4) {
        sums.cycles4++;
        sums.fourCycleProductSum += cycleProduct;
      } else {
        sums.cycles5++;
        sums.fiveCycleProductSum += cycleProduct;
      }
      return;
    }
    for (const next of adjacency[last]) {
      if (next <= start || used[next]) continue;
      const value = weights[last][next];
      used[next] = true;
      path.push(next);
      walkSparseCycle(start, next, length, product * value, squareSum + value ** 2);
      path.pop();
      used[next] = false;
    }
  }
}

function addFourCycle(weights, a, b, c, d, sums) {
  return addCycle(weights, [a, b, c, d], sums, "cycles4");
}

function addCycle(weights, cycle, sums, countKey) {
  let product = 1;
  for (let i = 0; i < cycle.length; i++) {
    const value = weights[cycle[i]][cycle[(i + 1) % cycle.length]];
    if (value === 0) return 0;
    product *= value;
  }
  sums[countKey]++;
  return product;
}

function approximateLnZ(model, beta, order) {
  const { kappa, counts } = calculateCumulants(model, order);
  let value = -beta * (model.offset || 0) + model.n * Math.log(2);
  for (let r = 2; r <= order; r++) value += ((-beta) ** r / factorial(r)) * kappa[r];
  return { value, kappa, counts };
}

function approximateConditionedLnZ(model, localSpinIndex, spinValue, beta, order) {
  const graph = buildConditionedAugmentedGraph(model, localSpinIndex, spinValue);
  const sums = fixedCycleSums(graph, order);
  const kappa = [0, 0, 0, 0, 0, 0];
  kappa[2] = sums.edgeSquareSum;
  if (order >= 3) kappa[3] = -6 * sums.triangleProductSum;
  if (order >= 4) kappa[4] = -2 * sums.edgeFourthSum + 24 * sums.fourCycleProductSum;
  if (order >= 5) kappa[5] = 40 * sums.triangleSquareProductSum - 120 * sums.fiveCycleProductSum;
  let value = -beta * (model.offset - model.fields[localSpinIndex] * spinValue) + (model.n - 1) * Math.log(2);
  for (let r = 2; r <= order; r++) value += ((-beta) ** r / factorial(r)) * kappa[r];
  return { value, kappa, counts: { edges: sums.edges, triangles: sums.triangles, cycles4: sums.cycles4, cycles5: sums.cycles5 } };
}

function buildConditionedAugmentedGraph(model, localSpinIndex, spinValue) {
  const size = model.n;
  const weights = Array.from({ length: size }, () => Array(size).fill(0));
  const adjacency = Array.from({ length: size }, () => []);
  const oldToAugmented = Array(model.n).fill(-1);
  const fields = [];
  for (let i = 0; i < model.n; i++) {
    if (i === localSpinIndex) continue;
    oldToAugmented[i] = fields.length + 1;
    fields.push(model.fields[i]);
  }
  const addEdge = (a, b, value) => {
    if (Math.abs(value) < EPS) return;
    weights[a][b] = value;
    weights[b][a] = value;
    adjacency[a].push(b);
    adjacency[b].push(a);
  };
  for (const edge of model.edges || edgesFromCouplings(model)) {
    if (edge.i === localSpinIndex || edge.j === localSpinIndex) {
      const other = edge.i === localSpinIndex ? edge.j : edge.i;
      const mapped = oldToAugmented[other];
      if (mapped > 0) fields[mapped - 1] += edge.J * spinValue;
    } else {
      addEdge(oldToAugmented[edge.i], oldToAugmented[edge.j], edge.J);
    }
  }
  fields.forEach((h, i) => addEdge(0, i + 1, h));
  adjacency.forEach((neighbors) => neighbors.sort((a, b) => a - b));
  return { size, weights, adjacency };
}

function exactEnergy(originalModel, spins) {
  let energy = originalModel.offset || 0;
  for (let i = 0; i < originalModel.n; i++) energy -= originalModel.fields[i] * spins[i];
  for (const edge of originalModel.edges || edgesFromCouplings(originalModel)) {
    energy -= edge.J * spins[edge.i] * spins[edge.j];
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
    model: hydrateModelMetadata({
      n: nextLabels.length,
      labels: nextLabels,
      fields: nextFields,
      couplings: nextCouplings,
      offset,
      normalizationStats: { ...(model.normalizationStats || {}) }
    })
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
    const graph = graphStats(current);
    const denseFastPath = useDensePrecheckFastPath(current, graph);
    const components = denseFastPath ? [] : connectedComponents(current, graph.adj);
    if (denseFastPath) {
      summary.disconnectedComponents = Math.max(summary.disconnectedComponents, current.n > 0 ? 1 : 0);
    } else {
      summary.disconnectedComponents = Math.max(summary.disconnectedComponents, components.length);
      summary.bipartiteComponents = Math.max(summary.bipartiteComponents, components.filter((c) => c.bipartite).length);
      summary.unfrustratedGaugeComponents = Math.max(summary.unfrustratedGaugeComponents, components.filter((c) => c.unfrustratedGauge).length);
    }

    const isolated = graph.degrees.findIndex((d) => d === 0);
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

    const symmetricComponent = denseFastPath
      ? (current.n > 1 && current.fields.every((h) => Math.abs(h) < EPS) ? { vertices: Array.from({ length: current.n }, (_, i) => i) } : null)
      : components.find((component) => component.vertices.length > 1 && component.vertices.every((i) => Math.abs(current.fields[i]) < EPS));
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

    const dominant = findDominantField(current, graph.radii);
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
  let changed = false;
  for (const [key, value] of [...model.couplings.entries()]) {
    if (Math.abs(value) < EPS) {
      model.couplings.delete(key);
      summary.zeroCouplingsRemoved++;
      changed = true;
    }
  }
  if (changed) hydrateModelMetadata(model);
}

function findDominantField(model, radii = null) {
  for (let i = 0; i < model.n; i++) {
    const h = model.fields[i];
    const radius = radii ? radii[i] : localRadius(model, i);
    if (Math.abs(h) > radius + EPS) return { local: i, h, radius };
  }
  return null;
}

function localRadius(model, local) {
  let radius = 0;
  for (const [key, value] of model.couplings.entries()) {
    const [a, b] = key.split(":").map(Number);
    if (a === local || b === local) radius += Math.abs(value);
  }
  return radius;
}

function graphStats(model) {
  const adj = (model.adjacency || []).map((neighbors) => neighbors.map((edge) => edge.to));
  const degrees = model.degrees ? [...model.degrees] : adj.map((neighbors) => neighbors.length);
  const radii = model.radii ? [...model.radii] : Array(model.n).fill(0);
  const edges = model.edges ? model.edges.length : model.couplings.size;
  const possible = model.n * (model.n - 1) / 2;
  const density = possible > 0 ? edges / possible : 0;
  return { adj, degrees, radii, edges, density };
}

function useDensePrecheckFastPath(model, graph) {
  return model.n > 12 && graph.density >= 0.45 && graph.degrees.every((degree) => degree > 0);
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
  const plus = approximateConditionedLnZ(model, localSpinIndex, 1, beta, order);
  const minus = approximateConditionedLnZ(model, localSpinIndex, -1, beta, order);
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
  const maxBacktracks = Math.max(0, Number(options.maxBacktracks ?? 0));
  const started = performance.now();
  let bestEnergy = Infinity;
  let bestSpins = null;
  let bestPrecheckSummary = {};
  let backtracks = 0;
  let decisions = 0;
  let firstComplete = null;
  const decisionLog = [];
  const branchLimit = maxBacktracks + 1;
  const queue = [{
    model: cloneModel(originalModel),
    spins: Array(originalModel.n).fill(null),
    depth: 0,
    trace: [],
    remaining: originalModel.n,
    penalty: 0,
    branchBacktracks: 0,
    label: "greedy"
  }];
  let explored = 0;
  let queued = 1;
  while (queue.length && explored < branchLimit) {
    queue.sort(compareBranches);
    followPath(queue.shift());
    explored++;
  }
  const fallback = firstComplete || { energy: Infinity, spins: [], trace: [] };
  return {
    energy: Number.isFinite(bestEnergy) ? bestEnergy : fallback.energy,
    spins: bestSpins || fallback.spins,
    bestEnergy: Number.isFinite(bestEnergy) ? bestEnergy : fallback.energy,
    bestSpins: bestSpins || fallback.spins,
    decisions,
    backtracks: Math.max(0, explored - 1),
    branchesExplored: explored,
    branchesQueued: queued,
    runtimeMs: performance.now() - started,
    decisionLog: decisionLog.length ? decisionLog : fallback.trace,
    precheckSummary: bestPrecheckSummary
  };

  function followPath(startNode) {
    let node = startNode;
    while (true) {
      const prechecked = preprocessState(node.model, node.spins, { countNormalization: !node.precheckSummary });
      node = {
        ...node,
        model: prechecked.model,
        spins: prechecked.spins,
        trace: [...node.trace, ...prechecked.log],
        precheckSummary: mergeSummaries(node.precheckSummary, prechecked.summary)
      };
      if (node.model.n === 0) {
        recordComplete(node);
        return;
      }
      const scores = [];
      for (let i = 0; i < node.model.n; i++) scores.push(scoreSpin(node.model, i, beta, order));
      scores.sort((a, b) => Math.abs(b.importance) - Math.abs(a.importance));
      const chosen = scores[0];
      decisions++;
      enqueueAlternatives(node, scores);
      node = applyBranch(node, chosen, chosen.preferred, false, scores);
    }
  }

  function enqueueAlternatives(node, scores) {
    if (scores.length <= 1) return;
    const best = scores[0];
    const alternatives = [];
    for (let optionIndex = 1; optionIndex < scores.length && alternatives.length < 2; optionIndex++) {
      alternatives.push({
        chosen: scores[optionIndex],
        branch: scores[optionIndex].preferred,
        optionIndex,
        penalty: scoreRegret(best, scores[optionIndex])
      });
    }
    alternatives.push({
      chosen: best,
      branch: best.alternate,
      optionIndex: scores.length,
      penalty: Math.abs(best.importance)
    });
    for (const alternative of alternatives.sort((a, b) => a.penalty - b.penalty || a.optionIndex - b.optionIndex)) {
      queue.push({
        ...applyBranch(snapshotNode(node), alternative.chosen, alternative.branch, true, scores),
        penalty: node.penalty + alternative.penalty + alternative.optionIndex * 1e-9,
        branchBacktracks: (node.branchBacktracks || 0) + 1,
        label: `regret ${formatPenalty(node.penalty + alternative.penalty)}`
      });
      queued++;
      pruneQueue();
    }
  }

  function pruneQueue() {
    if (queue.length + explored > branchLimit) {
      queue.sort(compareBranches);
      queue.splice(branchLimit - explored);
    }
  }

  function compareBranches(a, b) {
    return a.remaining - b.remaining || a.penalty - b.penalty;
  }

  function applyBranch(node, chosen, branch, isAlternate, scores) {
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
      candidateScores: scores,
      backtracksSoFar: (node.branchBacktracks || 0) + (isAlternate ? 1 : 0),
      fieldUpdates: fixed.fieldUpdates,
      constantChange: fixed.constantChange
    };
    return { model: fixed.model, spins, depth: node.depth + 1, remaining: fixed.model.n, trace: [...node.trace, step], penalty: node.penalty || 0, branchBacktracks: node.branchBacktracks || 0, precheckSummary: node.precheckSummary };
  }

  function snapshotNode(node) {
    return {
      model: node.model,
      spins: [...node.spins],
      depth: node.depth,
      remaining: node.model.n,
      trace: [...node.trace],
      penalty: node.penalty || 0,
      branchBacktracks: node.branchBacktracks || 0,
      precheckSummary: node.precheckSummary
    };
  }

  function recordComplete(node) {
    const energy = exactEnergy(originalModel, node.spins);
    if (!firstComplete) firstComplete = { energy, spins: [...node.spins], trace: node.trace };
    if (energy < bestEnergy - 1e-10) {
      bestEnergy = energy;
      bestSpins = [...node.spins];
      bestPrecheckSummary = node.precheckSummary || {};
      decisionLog.splice(0, decisionLog.length, ...node.trace);
    }
  }
}

function scoreRegret(best, candidate) {
  return Math.max(0, Math.abs(best.importance) - Math.abs(candidate.importance));
}

function formatPenalty(value) {
  if (!Number.isFinite(value)) return String(value);
  if (Math.abs(value) < 1e-10) return "0";
  return Number(value.toPrecision(8)).toString();
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

function buildSpinGlass100Text() {
  let seed = 501337;
  const rand = () => {
    seed = (1664525 * seed + 1013904223) >>> 0;
    return seed / 2 ** 32;
  };
  const n = 100;
  const targetEdges = 240;
  const lines = [`N ${n}`, "", "FIELDS"];
  for (let i = 1; i <= n; i++) lines.push(`${i} ${Math.round((rand() * 1.2 - 0.6) * 100) / 100}`);
  lines.push("", "COUPLINGS");
  const used = new Set();
  for (let i = 1; i < n; i++) {
    const j = i + 1;
    const J = Math.round((rand() * 2 - 1) * 100) / 100 || 0.25;
    used.add(`${i}:${j}`);
    lines.push(`${i} ${j} ${J}`);
  }
  while (used.size < targetEdges) {
    const i = 1 + Math.floor(rand() * n);
    const j = 1 + Math.floor(rand() * n);
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
  spinGlass100: {
    name: "100-node spin glass",
    benchmark: "No built-in certificate; use for visual/search behavior.",
    text: buildSpinGlass100Text()
  },
  frustratedRing9: {
    name: "9-spin frustrated AF ring",
    benchmark: "Known ground energy -7; one unsatisfied bond is unavoidable.",
    text: "N 9\n\nFIELDS\n1 0\n2 0\n3 0\n4 0\n5 0\n6 0\n7 0\n8 0\n9 0\n\nCOUPLINGS\n1 2 -1\n2 3 -1\n3 4 -1\n4 5 -1\n5 6 -1\n6 7 -1\n7 8 -1\n8 9 -1\n9 1 -1"
  }
};

const $ = (id) => document.getElementById(id);
const VISUAL_SPIN_LIMIT = 100;
let previewTimer = null;
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
  sourceText: ""
};

function options() {
  return {
    beta: Number($("beta").value),
    order: 5,
    maxBacktracks: Number($("backtracks").value)
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
    write("Running fifth-order Taylor/cumulant solver. For 100 nodes this can take a while in the browser.");
    setTimeout(() => {
      try {
        const model = currentModel();
        const opts = options();
        const cumulants = calculateCumulants(model, 5);
        const lnz = approximateLnZ(model, opts.beta, 5);
        const result = solveTaylor(model, opts);
        prepareVisualization(model, result);
        showFinalAssignment();
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
          "Answer",
          `Final spins: ${formatSpins(result.spins)}`,
          `Exact final energy: ${formatNumber(result.energy)}`,
          "",
          "Run summary",
          `N: ${model.n}`,
          `Couplings: ${stats.couplings}`,
          `Nonzero fields: ${stats.fields}`,
          `Beta: ${opts.beta}`,
          "Taylor order: 5",
          `Decisions: ${result.decisions}`,
          `Backtrack branches explored: ${result.backtracks}`,
          `Total branches explored: ${result.branchesExplored || 1}`,
          `Runtime: ${formatNumber(result.runtimeMs)} ms`,
          `ln Z approximation: ${formatNumber(lnz.value)}`,
          `Cycles: triangles=${cumulants.counts.triangles}, C4=${cumulants.counts.cycles4}, C5=${cumulants.counts.cycles5}`,
          "",
          "Precheck summary",
          ...formatPrecheckSummary(result.precheckSummary)
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
  visualState.sourceText = $("isingInput")?.value || "";
  visualState.assignments = new Map();
  visualState.reasons = new Map();
  visualState.importance = new Map();
  visualState.selected = null;
  visualState.step = 0;
  if (model.n > VISUAL_SPIN_LIMIT) {
    visualState.nodes = [];
    visualState.edges = [];
    drawVisualization();
    return;
  }
  const nodes = layoutNodes(model.n, 450, 380, 325);
  const rendered = new Set(nodes.map((n) => n.index));
  visualState.edges = (model.edges || [...model.couplings.entries()].map(([key, J]) => {
    const [a, b] = key.split(":").map(Number);
    return { i: a, j: b, J };
  })).map((edge) => ({ a: edge.i, b: edge.j, J: edge.J }))
    .filter((e) => rendered.has(e.a) && rendered.has(e.b));
  visualState.nodes = nodes;
  drawVisualization();
}

function previewCurrentGraph() {
  try {
    const text = $("isingInput").value;
    if (text !== visualState.sourceText) prepareVisualization(parseIsingText(text), null);
    else drawVisualization();
  } catch (error) {
    const box = $("visualInfo");
    if (box) box.textContent = `Input error: ${error.message}`;
  }
}

function scheduleGraphPreview() {
  clearTimeout(previewTimer);
  previewTimer = setTimeout(previewCurrentGraph, 180);
}

function layoutNodes(count, centerX, centerY, maxRadius) {
  if (count <= 0) return [];
  if (count === 1) return [{ index: 0, x: centerX, y: centerY }];
  const rings = Math.max(1, Math.ceil(Math.sqrt(count / 8)));
  const nodes = [];
  let placed = 0;
  for (let ring = 1; ring <= rings && placed < count; ring++) {
    const remaining = count - placed;
    const ringCapacity = ring === rings ? remaining : Math.max(6, Math.round((count * ring) / ((rings * (rings + 1)) / 2)));
    const radius = Math.max(42, (maxRadius * ring) / rings);
    const phase = (ring % 2) * Math.PI / Math.max(1, ringCapacity);
    for (let k = 0; k < ringCapacity && placed < count; k++, placed++) {
      const angle = phase + (k / ringCapacity) * Math.PI * 2;
      nodes.push({
        index: placed,
        x: centerX + Math.cos(angle) * radius,
        y: centerY + Math.sin(angle) * radius
      });
    }
  }
  return nodes;
}

function drawVisualization() {
  const canvas = $("spinCanvas");
  if (!canvas || !visualState.model) return;
  const ctx = canvas.getContext("2d");
  const rect = canvas.getBoundingClientRect();
  const scale = Math.max(1, window.devicePixelRatio || 1);
  canvas.width = Math.max(640, Math.floor(rect.width * scale));
  const logicalHeight = 760;
  canvas.height = Math.max(620, Math.floor(rect.height * scale));
  ctx.setTransform(canvas.width / 900, 0, 0, canvas.height / logicalHeight, 0, 0);
  ctx.clearRect(0, 0, 900, logicalHeight);
  if (visualState.model.n > VISUAL_SPIN_LIMIT) {
    ctx.fillStyle = "#edf6f5";
    ctx.font = "18px system-ui";
    ctx.textAlign = "center";
    ctx.fillText("Visualization unavailable for more than 100 spins.", 450, 365);
    updateVisualInfo();
    return;
  }
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

function resetVisualization() {
  prepareVisualization(currentModel(), visualState.result);
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
  if (model.n > VISUAL_SPIN_LIMIT) {
    box.textContent = `Visualization unavailable for ${model.n} spins. The solver can still run on the full input.`;
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
  $("isingInput").value = examples.spinGlass100.text;
  $("isingInput").addEventListener("input", scheduleGraphPreview);
  $("runSolver").addEventListener("click", runSolver);
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
  $("spinCanvas").addEventListener("click", selectCanvasSpin);
  window.addEventListener("resize", drawVisualization);
  prepareVisualization(currentModel(), null);
  setSummary([["spins", 100], ["example", "100-node spin glass"], ["visual", "ready"]]);
  write("Ready. Preview shows unresolved spins. Press Run Taylor solver to compute assignments and draw the final Visualizer state.");
}

init();
