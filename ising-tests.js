import {
  calculateCumulants,
  createModel,
  edgeKey,
  exactEnergy
} from "./ising-math.js";
import { fixSpin, preprocessState } from "./ising-precheck.js";

export function runCumulantVerification({ tolerance = 1e-8 } = {}) {
  const checks = [];
  const add = (name, actual, expected) => {
    checks.push({ name, passed: Math.abs(actual - expected) < tolerance, actual, expected });
  };

  const singleEdge = calculateCumulants(createModel(2, [{ i: 0, j: 1, J: 2 }], [0, 0])).kappa;
  add("single edge kappa2", singleEdge[2], 4);
  add("single edge kappa3", singleEdge[3], 0);
  add("single edge kappa4", singleEdge[4], -32);
  add("single edge kappa5", singleEdge[5], 0);

  const triangle = calculateCumulants(createModel(3, [
    { i: 0, j: 1, J: 1 },
    { i: 1, j: 2, J: 1 },
    { i: 0, j: 2, J: 1 }
  ], [0, 0, 0])).kappa;
  add("ferromagnetic triangle kappa2", triangle[2], 3);
  add("ferromagnetic triangle kappa3", triangle[3], -6);
  add("ferromagnetic triangle kappa4", triangle[4], -6);
  add("ferromagnetic triangle kappa5", triangle[5], 120);

  const ghostTriangle = calculateCumulants(createModel(2, [{ i: 0, j: 1, J: 2 }], [0.5, -0.25])).kappa;
  add("field ghost triangle kappa2", ghostTriangle[2], 4.3125);
  add("field ghost triangle kappa3", ghostTriangle[3], 1.5);

  return { passed: checks.every((check) => check.passed), checks };
}

export function runBuiltInChecks() {
  const checks = [
    {
      name: "Two-spin ferromagnet",
      model: createModel(2, [{ i: 0, j: 1, J: 1 }], [0, 0]),
      spins: [1, 1],
      expected: -1
    },
    {
      name: "Two-spin antiferromagnet",
      model: createModel(2, [{ i: 0, j: 1, J: -1 }], [0, 0]),
      spins: [1, -1],
      expected: -1
    },
    {
      name: "Energy evaluator sanity",
      model: createModel(2, [{ i: 0, j: 1, J: 1 }], [0.25, -0.5]),
      spins: [1, -1],
      expected: 0.25
    }
  ];
  return checks.map((check) => {
    const actual = exactEnergy(check.model, check.spins);
    return { name: check.name, passed: Math.abs(actual - check.expected) < 1e-10, expected: check.expected, actual };
  });
}

export function runPrecheckVerification({ tolerance = 1e-9 } = {}) {
  const checks = [];
  const add = (name, passed, detail = {}) => checks.push({ name, passed, ...detail });

  const zero = createModel(2, [{ i: 0, j: 1, J: 0 }], [0.4, -0.1]);
  add("Zero coupling is removed", zero.couplings.size === 0, { removed: zero.normalizationStats.zeroCouplingsRemoved });

  const duplicate = createModel(2, [{ i: 0, j: 1, J: 0.7 }, { i: 1, j: 0, J: -0.2 }], [0, 0]);
  add("Duplicate coupling merges by summation", Math.abs(duplicate.couplings.get(edgeKey(0, 1)) - 0.5) < tolerance, {
    merged: duplicate.normalizationStats.duplicateEdgesMerged
  });

  const self = createModel(1, [{ i: 0, j: 0, J: 1.7 }], [0]);
  add("Self coupling converts to constant", self.couplings.size === 0 && Math.abs(self.offset + 1.7) < tolerance, {
    offset: self.offset
  });

  const isoPlus = preprocessState(createModel(1, [], [2]), [null]);
  add("Isolated positive field fixes +1", isoPlus.spins[0] === 1 && isoPlus.model.n === 0);

  const isoMinus = preprocessState(createModel(1, [], [-2]), [null]);
  add("Isolated negative field fixes -1", isoMinus.spins[0] === -1 && isoMinus.model.n === 0);

  const symmetryModel = createModel(3, [{ i: 0, j: 1, J: 1 }, { i: 1, j: 2, J: -0.7 }], [0, 0, 0]);
  const spin = [1, -1, 1];
  const flipped = spin.map((s) => -s);
  const symmetryReduced = preprocessState(symmetryModel, [null, null, null]);
  add("Zero-field symmetry preserves flipped energy", Math.abs(exactEnergy(symmetryModel, spin) - exactEnergy(symmetryModel, flipped)) < tolerance
    && symmetryReduced.log.some((l) => l.reason === "SYMMETRY"));

  const dominant = createModel(3, [{ i: 0, j: 1, J: -1 }, { i: 0, j: 2, J: 0.5 }], [2, 0, 0]);
  const dominantPrep = preprocessState(dominant, [null, null, null]);
  add("Dominant strict field fixes predicted sign", dominantPrep.spins[0] === 1);

  const equality = createModel(2, [{ i: 0, j: 1, J: -1 }], [1, 0]);
  const equalityPrep = preprocessState(equality, [null, null]);
  add("Dominant equality case is not forced", !equalityPrep.log.some((l) => l.rule === "dominant local field" && l.spin === 1));

  const propagation = createModel(3, [{ i: 0, j: 1, J: 1.2 }, { i: 0, j: 2, J: -0.4 }], [0.3, -0.2, 0.7]);
  const fixed = fixSpin(propagation, 0, -1);
  const original = exactEnergy(propagation, [-1, 1, -1]);
  const reduced = exactEnergy(fixed.model, [1, -1]);
  add("Propagation preserves a conditioned energy", Math.abs(original - reduced) < tolerance);

  return { passed: checks.every((check) => check.passed), checks };
}

if (typeof process !== "undefined" && process.argv[1]?.endsWith("ising-tests.js")) {
  const builtIns = runBuiltInChecks();
  const cumulants = runCumulantVerification();
  const prechecks = runPrecheckVerification();
  console.log(JSON.stringify({ builtIns, cumulants, prechecks }, null, 2));
  if (!builtIns.every((test) => test.passed) || !cumulants.passed || !prechecks.passed) process.exit(1);
}
