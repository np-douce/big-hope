import { EPS, cloneModel, edgeKey } from "./ising-math.js";

export function preprocessState(model, spins, options = {}) {
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
      const rule = Math.abs(h) < EPS ? "isolated zero-field spin" : "isolated field spin";
      const result = fixSpin(current, isolated, value);
      assignments[result.originalSpin] = value;
      current = result.model;
      summary.isolatedSpinsFixed++;
      if (reason === "DEGENERATE") summary.degenerateSpinsFixed++;
      log.push(makeLog(result, reason, rule, `h = ${format(h)}`));
      changed = true;
      continue;
    }

    const treeComponent = options.solveTreesExactly ? components.find((c) => c.vertices.length > 1 && c.edges === c.vertices.length - 1) : null;
    if (treeComponent) {
      const exact = solveTreeComponent(current, treeComponent.vertices);
      for (const originalSpin of exact.order) {
        const local = current.labels.indexOf(originalSpin);
        if (local < 0) continue;
        const value = exact.assignment.get(originalSpin);
        const result = fixSpin(current, local, value);
        assignments[result.originalSpin] = value;
        current = result.model;
        log.push(makeLog(result, "FORCED", "exact tree component DP", `component energy = ${format(exact.energy)}`));
      }
      summary.treeComponentsSolved++;
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
      log.push(makeLog(result, "SYMMETRY", "zero-field global spin symmetry", "lowest-numbered spin fixed to +1"));
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
      log.push(makeLog(result, "FORCED", "dominant local field", `|h| = ${format(Math.abs(dominant.h))}, sum |J| = ${format(dominant.radius)}`));
      changed = true;
    }
  }

  summary.remainingSpins = current.n;
  return { model: current, spins: assignments, log, summary, components: connectedComponents(current) };
}

export function connectedComponents(model) {
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

export function fixSpin(model, localSpinIndex, value) {
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
    field: originalField,
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
  for (const start of vertices) {
    if (!color.has(start)) {
      color.set(start, 1);
      sign.set(start, 1);
    }
    const queue = [start];
    for (let p = 0; p < queue.length; p++) {
      const v = queue[p];
      for (const n of adj[v]) {
        if (!vertices.includes(n)) continue;
        if (!color.has(n)) {
          color.set(n, -color.get(v));
          queue.push(n);
        } else if (color.get(n) === color.get(v)) {
          bipartite = false;
        }
        const coupling = model.couplings.get(edgeKey(v, n));
        const wanted = Math.sign(coupling || 1) * sign.get(v);
        if (!sign.has(n)) {
          sign.set(n, wanted);
        } else if (sign.get(n) !== wanted) {
          unfrustratedGauge = false;
        }
      }
    }
  }
  return { bipartite, unfrustratedGauge };
}

function solveTreeComponent(model, vertices) {
  const vertexSet = new Set(vertices);
  const adj = new Map(vertices.map((v) => [v, []]));
  for (const [key, value] of model.couplings.entries()) {
    const [a, b] = key.split(":").map(Number);
    if (vertexSet.has(a) && vertexSet.has(b)) {
      adj.get(a).push({ node: b, J: value });
      adj.get(b).push({ node: a, J: value });
    }
  }
  const root = vertices[0];
  const parent = new Map([[root, -1]]);
  const order = [root];
  for (let p = 0; p < order.length; p++) {
    const v = order[p];
    for (const edge of adj.get(v)) {
      if (edge.node === parent.get(v)) continue;
      parent.set(edge.node, v);
      order.push(edge.node);
    }
  }

  const dp = new Map();
  const choice = new Map();
  for (const v of [...order].reverse()) {
    const row = new Map();
    const pick = new Map();
    for (const s of [-1, 1]) {
      let energy = -model.fields[v] * s;
      const childChoices = new Map();
      for (const edge of adj.get(v)) {
        if (parent.get(edge.node) !== v) continue;
        const child = dp.get(edge.node);
        const minus = child.get(-1) - edge.J * s * -1;
        const plus = child.get(1) - edge.J * s * 1;
        const childSpin = plus < minus ? 1 : -1;
        energy += Math.min(minus, plus);
        childChoices.set(edge.node, childSpin);
      }
      row.set(s, energy);
      pick.set(s, childChoices);
    }
    dp.set(v, row);
    choice.set(v, pick);
  }

  const rootSpin = dp.get(root).get(1) <= dp.get(root).get(-1) ? 1 : -1;
  const assignment = new Map();
  const stack = [{ node: root, spin: rootSpin }];
  while (stack.length) {
    const { node, spin } = stack.pop();
    assignment.set(model.labels[node], spin);
    for (const [child, childSpin] of choice.get(node).get(spin).entries()) {
      stack.push({ node: child, spin: childSpin });
    }
  }
  return {
    energy: Math.min(dp.get(root).get(-1), dp.get(root).get(1)),
    assignment,
    order: vertices.map((v) => model.labels[v])
  };
}

function makeLog(result, reason, rule, detail) {
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

function format(value) {
  if (Math.abs(value) < 1e-10) return "0";
  return Number(value.toPrecision(8)).toString();
}
