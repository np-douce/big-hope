import { approximateConditionedLnZ, cloneModel, exactEnergy } from "./ising-math.js";
import { fixSpin, preprocessState } from "./ising-precheck.js";

export function conditionModel(model, localSpinIndex, value) {
  return fixSpin(model, localSpinIndex, value).model;
}

export function scoreSpin(model, localSpinIndex, beta, order) {
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

export function solveTaylor(originalModel, options = {}) {
  const beta = Number(options.beta ?? 0.5);
  const order = Number(options.order ?? 5);
  const maxBacktracks = Math.max(0, Number(options.maxBacktracks ?? 0));
  const started = performance.now();
  let bestEnergy = Infinity;
  let bestSpins = null;
  let backtracks = 0;
  let decisions = 0;
  let firstComplete = null;
  let bestPrecheckSummary = {};
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
      const prechecked = preprocessState(node.model, node.spins, {
        solveTreesExactly: Boolean(options.solveTreesExactly),
        countNormalization: !node.precheckSummary
      });
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
