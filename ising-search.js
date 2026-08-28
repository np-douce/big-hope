import { approximateLnZ, cloneModel, exactEnergy } from "./ising-math.js";
import { fixSpin, preprocessState } from "./ising-precheck.js";

export function conditionModel(model, localSpinIndex, value) {
  return fixSpin(model, localSpinIndex, value).model;
}

export function scoreSpin(model, localSpinIndex, beta, order) {
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

export function solveTaylor(originalModel, options = {}) {
  const beta = Number(options.beta ?? 0.5);
  const order = Number(options.order ?? 5);
  const mode = options.mode || "max";
  const maxBacktracks = Math.max(0, Number(options.maxBacktracks ?? 0));
  const started = performance.now();
  let bestEnergy = Infinity;
  let bestSpins = null;
  let backtracks = 0;
  let decisions = 0;
  let firstComplete = null;
  let bestPrecheckSummary = {};
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
    const prechecked = preprocessState(node.model, node.spins, {
      solveTreesExactly: Boolean(options.solveTreesExactly),
      countNormalization: !node.precheckSummary
    });
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
      const nextModel = fixed.model;
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
      stack.push({ model: nextModel, spins, depth: node.depth + 1, trace: [...node.trace, step], deviations, precheckSummary: node.precheckSummary });
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
    if (typeof value === "number") {
      if (key === "remainingSpins") next[key] = value;
      else if (key === "originalSpins" || key === "originalCouplings") next[key] ??= value;
      else next[key] = (next[key] || 0) + value;
    } else {
      next[key] = value;
    }
  }
  return next;
}
