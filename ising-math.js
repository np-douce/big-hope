export const EPS = 1e-12;

export function edgeKey(a, b) {
  return a < b ? `${a}:${b}` : `${b}:${a}`;
}

export function hydrateModelMetadata(model) {
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

export function cloneModel(model) {
  return hydrateModelMetadata({
    n: model.n,
    labels: [...model.labels],
    fields: [...model.fields],
    couplings: new Map(model.couplings),
    offset: model.offset || 0,
    normalizationStats: { ...(model.normalizationStats || {}) }
  });
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
  return hydrateModelMetadata(model);
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

export function calculateCumulants(model, maxOrder = 5) {
  const graph = buildAugmentedGraph(model);
  const sums = fixedCycleSums(graph, maxOrder);
  const kappa = [0, 0, 0, 0, 0, 0];
  kappa[2] = sums.edgeSquareSum;
  if (maxOrder >= 3) kappa[3] = -6 * sums.triangleProductSum;
  if (maxOrder >= 4) kappa[4] = -2 * sums.edgeFourthSum + 24 * sums.fourCycleProductSum;
  if (maxOrder >= 5) kappa[5] = 40 * sums.triangleSquareProductSum - 120 * sums.fiveCycleProductSum;

  return {
    kappa,
    counts: {
      edges: sums.edges,
      triangles: sums.triangles,
      cycles4: sums.cycles4,
      cycles5: sums.cycles5
    }
  };
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

export function approximateLnZ(model, beta, order) {
  const { kappa, counts } = calculateCumulants(model, order);
  let value = -beta * (model.offset || 0) + model.n * Math.log(2);
  for (let r = 2; r <= order; r++) {
    value += ((-beta) ** r / factorial(r)) * kappa[r];
  }
  return { value, kappa, counts };
}

export function approximateConditionedLnZ(model, localSpinIndex, spinValue, beta, order) {
  const graph = buildConditionedAugmentedGraph(model, localSpinIndex, spinValue);
  const sums = fixedCycleSums(graph, order);
  const kappa = [0, 0, 0, 0, 0, 0];
  kappa[2] = sums.edgeSquareSum;
  if (order >= 3) kappa[3] = -6 * sums.triangleProductSum;
  if (order >= 4) kappa[4] = -2 * sums.edgeFourthSum + 24 * sums.fourCycleProductSum;
  if (order >= 5) kappa[5] = 40 * sums.triangleSquareProductSum - 120 * sums.fiveCycleProductSum;

  let value = -beta * (model.offset - model.fields[localSpinIndex] * spinValue) + (model.n - 1) * Math.log(2);
  for (let r = 2; r <= order; r++) {
    value += ((-beta) ** r / factorial(r)) * kappa[r];
  }
  return {
    value,
    kappa,
    counts: { edges: sums.edges, triangles: sums.triangles, cycles4: sums.cycles4, cycles5: sums.cycles5 }
  };
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

export function exactEnergy(originalModel, spinsByOriginalIndex) {
  let energy = originalModel.offset || 0;
  for (let i = 0; i < originalModel.n; i++) energy -= originalModel.fields[i] * spinsByOriginalIndex[i];
  for (const edge of originalModel.edges || edgesFromCouplings(originalModel)) {
    energy -= edge.J * spinsByOriginalIndex[edge.i] * spinsByOriginalIndex[edge.j];
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
