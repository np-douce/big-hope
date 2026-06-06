const DEFAULT_ZERO_EDGES_TO_TRY = 50;
const MAX_VISIBLE_LABEL_TERMS = 65535;
const MAX_VISIBLE_GROUPS = 200;

const edgeInput = document.getElementById("edgeInput");
const zeroEdgeCapInput = document.getElementById("zeroEdgeCapInput");
const calculateBtn = document.getElementById("calculateBtn");
const sampleBtn = document.getElementById("sampleBtn");
const clearBtn = document.getElementById("clearBtn");
const summary = document.getElementById("summary");
const subtractionResult = document.getElementById("subtractionResult");

class Rational {
  constructor(numerator, denominator = 1n) {
    if (denominator === 0n) throw new Error("Rational denominator cannot be zero.");
    if (denominator < 0n) {
      numerator = -numerator;
      denominator = -denominator;
    }
    const divisor = gcd(absBigInt(numerator), denominator);
    this.n = numerator / divisor;
    this.d = denominator / divisor;
  }

  static zero() {
    return new Rational(0n);
  }

  static fromInt(value) {
    return new Rational(BigInt(value));
  }

  add(other) {
    return new Rational((this.n * other.d) + (other.n * this.d), this.d * other.d);
  }

  sub(other) {
    return new Rational((this.n * other.d) - (other.n * this.d), this.d * other.d);
  }

  mulInt(value) {
    return new Rational(this.n * BigInt(value), this.d);
  }

  divInt(value) {
    return new Rational(this.n, this.d * BigInt(value));
  }

  equalsInt(value) {
    return this.n === BigInt(value) * this.d;
  }

  toDecimal(places = 10) {
    if (this.d === 1n) return this.n.toString();
    const sign = this.n < 0n ? "-" : "";
    const numerator = absBigInt(this.n);
    const integerPart = numerator / this.d;
    let remainder = numerator % this.d;
    let fraction = "";
    for (let i = 0; i < places && remainder !== 0n; i++) {
      remainder *= 10n;
      fraction += (remainder / this.d).toString();
      remainder %= this.d;
    }
    fraction = fraction.replace(/0+$/, "");
    return fraction ? `${sign}${integerPart}.${fraction}` : `${sign}${integerPart}`;
  }

  toDisplay() {
    if (this.d === 1n) return this.n.toString();
    return `${this.n}/${this.d} ~= ${this.toDecimal(10)}`;
  }
}

function absBigInt(value) {
  return value < 0n ? -value : value;
}

function gcd(a, b) {
  while (b !== 0n) {
    const next = a % b;
    a = b;
    b = next;
  }
  return a || 1n;
}

function edgeKey(u, v) {
  return u < v ? `${u}-${v}` : `${v}-${u}`;
}

function parseEdges(text) {
  const edges = [];
  let maxVertex = 0;
  const seen = new Set();

  text.split(/\r?\n/).forEach((line, index) => {
    const trimmed = line.trim();
    if (!trimmed) return;

    const numbers = trimmed.match(/\d+/g);
    if (!numbers || numbers.length !== 2) {
      throw new Error(`Line ${index + 1} must contain exactly two vertex numbers.`);
    }

    let u = Number(numbers[0]);
    let v = Number(numbers[1]);
    if (!Number.isInteger(u) || !Number.isInteger(v) || u < 1 || v < 1) {
      throw new Error(`Line ${index + 1} has an invalid vertex number.`);
    }
    if (u === v) {
      throw new Error(`Line ${index + 1} is a self-edge. Use two different vertices.`);
    }
    if (u > v) [u, v] = [v, u];

    const key = edgeKey(u, v);
    if (!seen.has(key)) {
      seen.add(key);
      edges.push({ u, v });
      maxVertex = Math.max(maxVertex, u, v);
    }
  });

  return { edges, n: maxVertex, seen };
}

function buildWeight(n, realEdges) {
  const weight = Array.from({ length: n + 1 }, () => Array(n + 1).fill(0));
  for (const edge of realEdges) {
    weight[edge.u][edge.v] = -1;
    weight[edge.v][edge.u] = -1;
  }
  return weight;
}

function allTourAverageExact(n, weight) {
  let edgeWeightSum = 0n;
  for (let i = 1; i < n; i++) {
    for (let j = i + 1; j <= n; j++) {
      edgeWeightSum += BigInt(weight[i][j]);
    }
  }
  return new Rational(2n * edgeWeightSum, BigInt(n - 1));
}

function factorialBigInt(value) {
  if (value < 0) return 0n;
  let result = 1n;
  for (let i = 2; i <= value; i++) result *= BigInt(i);
  return result;
}

function omegaAllTours(n) {
  return factorialBigInt(n - 1) / 2n;
}

function zeroEdgesFromGraph(n, realEdgeKeys) {
  const zeroEdges = [];
  for (let i = 1; i < n; i++) {
    for (let j = i + 1; j <= n; j++) {
      if (!realEdgeKeys.has(edgeKey(i, j))) zeroEdges.push({ u: i, v: j });
    }
  }
  return zeroEdges;
}

function readZeroEdgeTryCap() {
  const rawValue = zeroEdgeCapInput ? zeroEdgeCapInput.value.trim() : "";
  const value = rawValue === "" ? DEFAULT_ZERO_EDGES_TO_TRY : Number(rawValue);
  if (!Number.isInteger(value) || value < 0) {
    throw new Error("Zero-edge try cap must be a nonnegative whole number.");
  }
  return value;
}

function orderZeroEdgesForPruning(n, zeroEdges) {
  const zeroDegree = Array(n + 1).fill(0);
  for (const edge of zeroEdges) {
    zeroDegree[edge.u]++;
    zeroDegree[edge.v]++;
  }

  const remaining = canonicalEdges(zeroEdges);
  const ordered = [];
  const orderedIncident = Array(n + 1).fill(0);

  while (remaining.length > 0) {
    let bestIndex = 0;
    let bestScore = -Infinity;
    let bestTie = "";

    for (let index = 0; index < remaining.length; index++) {
      const edge = remaining[index];
      const sharedWithOrder = orderedIncident[edge.u] + orderedIncident[edge.v];
      const zeroDegreeScore = zeroDegree[edge.u] + zeroDegree[edge.v];
      const pressureScore = Math.max(orderedIncident[edge.u], orderedIncident[edge.v]);
      const score = (1000 * sharedWithOrder) + (100 * pressureScore) + zeroDegreeScore;
      const tie = edgeKey(edge.u, edge.v);

      if (score > bestScore || (score === bestScore && tie < bestTie)) {
        bestIndex = index;
        bestScore = score;
        bestTie = tie;
      }
    }

    const [nextEdge] = remaining.splice(bestIndex, 1);
    ordered.push(nextEdge);
    orderedIncident[nextEdge.u]++;
    orderedIncident[nextEdge.v]++;
  }

  return ordered;
}

function certifiesHamiltonianAverage(n, result) {
  return result.omegaValid > 0n && result.averageValid && result.averageValid.equalsInt(-n);
}

function attachSearchMetadata(result, metadata) {
  return Object.assign(result, metadata);
}

function runSmartSubtractionSearch(n, weight, zeroEdges, zeroEdgeTryCap = DEFAULT_ZERO_EDGES_TO_TRY) {
  const orderedZeroEdges = orderZeroEdgesForPruning(n, zeroEdges);
  const forcedSubsetDepth = n;
  const maxPrefix = Math.min(orderedZeroEdges.length, zeroEdgeTryCap);
  const realEdges = [];
  for (let i = 1; i < n; i++) {
    for (let j = i + 1; j <= n; j++) {
      if (weight[i][j] === -1) realEdges.push({ u: i, v: j });
    }
  }
  let termsVisited = 0n;
  let nonzeroTerms = 0n;
  let zeroOmegaTerms = 0n;
  let prunedTerms = 0n;
  const groups = new Map();
  const statsBySignature = new Map();
  let structuralCacheHits = 0n;

  function countCompletions(remainingChoices) {
    return 2n ** BigInt(remainingChoices);
  }

  function statsFromSignature(signature) {
    if (signature.kind === "zero") return zeroStats();

    if (signature.kind === "empty") {
      const omega = omegaAllTours(n);
      const average = allTourAverageExact(n, weight);
      return { omega, average, total: average.mulInt(omega), impossible: false };
    }

    if (signature.kind === "cycle") {
      const omega = 1n;
      const average = Rational.zero();
      return { omega, average, total: average.mulInt(omega), impossible: false };
    }

    if (signature.denominator === 0) {
      const omega = 2n ** BigInt(signature.ceCount - 1);
      const average = new Rational(signature.closingWeight);
      return { omega, average, total: average.mulInt(omega), impossible: false };
    }

    const d = BigInt(signature.denominator);
    const average = Rational.zero()
      .add(new Rational(signature.cceSum, 2n * d))
      .add(new Rational(signature.neSum, d))
      .add(new Rational(2n * signature.ieSum, d));
    const omega = (2n ** BigInt(signature.ceCount - 1)) * factorialBigInt(signature.denominator);

    return {
      omega,
      average,
      total: average.mulInt(omega),
      impossible: false
    };
  }

  function getStatsForSignature(signature) {
    const cached = statsBySignature.get(signature.key);
    if (cached) {
      structuralCacheHits++;
      return cached;
    }

    const stats = statsFromSignature(signature);
    statsBySignature.set(signature.key, stats);
    return stats;
  }

  function addSignature(signature, subsetSize) {
    termsVisited++;
    const sign = subsetSize % 2 === 0 ? 1n : -1n;
    const key = signature.key;
    let group = groups.get(key);

    if (group) {
      nonzeroTerms++;
      structuralCacheHits++;
      group.coefficient += sign;
      group.subsets++;
      return;
    }

    const stats = getStatsForSignature(signature);
    if (stats.omega === 0n) {
      zeroOmegaTerms++;
      return;
    }

    nonzeroTerms++;
    group = {
      coefficient: sign,
      omega: stats.omega,
      average: stats.average,
      total: stats.total,
      subsets: 1n,
      signature: signature.label
    };
    groups.set(key, group);
  }

  function skipZeroTerms(count) {
    termsVisited += count;
    zeroOmegaTerms += count;
    prunedTerms += count;
  }

  function buildResult(usedZeroEdges, stopReason) {
    let omegaValid = 0n;
    let totalValid = Rational.zero();
    let cancelledGroups = 0;
    const activeGroups = [];

    for (const group of groups.values()) {
      if (group.coefficient === 0n) {
        cancelledGroups++;
      } else {
        omegaValid += group.coefficient * group.omega;
        totalValid = totalValid.add(group.total.mulInt(group.coefficient));
        activeGroups.push(group);
      }
    }

    activeGroups.sort((a, b) => {
      const absA = absBigInt(a.coefficient);
      const absB = absBigInt(b.coefficient);
      if (absA !== absB) return absA > absB ? -1 : 1;
      if (a.subsets === b.subsets) return 0;
      return a.subsets > b.subsets ? -1 : 1;
    });

    const averageValid = omegaValid !== 0n
      ? totalValid.divInt(omegaValid)
      : null;

    return attachSearchMetadata({
      omegaValid,
      totalValid,
      averageValid,
      termsVisited,
      nonzeroTerms,
      zeroOmegaTerms,
      prunedTerms,
      groupedTerms: groups.size,
      structuralGroups: statsBySignature.size,
      structuralCacheHits,
      activeGroupedTerms: activeGroups.length,
      cancelledGroups,
      activeGroups
    }, {
      orderedZeroEdges,
      usedZeroEdges,
      totalZeroEdges: orderedZeroEdges.length,
      zeroEdgeTryCap,
      forcedSubsetDepth,
      stoppedEarly: usedZeroEdges < orderedZeroEdges.length && (stopReason === "average" || stopReason === "empty"),
      stopReason,
      prefixLimited: stopReason === "cap"
    });
  }

  const degree = Array(n + 1).fill(0);
  const touched = Array(n + 1).fill(false);
  const adjacency = Array.from({ length: n + 1 }, () => []);
  const visitMarks = Array(n + 1).fill(0);
  const chosen = [];
  let touchedCount = 0;
  let visitStamp = 0;

  function areConnected(u, v) {
    if (!touched[u] || !touched[v]) return false;

    visitStamp++;
    const stack = [u];
    visitMarks[u] = visitStamp;

    while (stack.length > 0) {
      const vertex = stack.pop();
      if (vertex === v) return true;

      for (const neighbor of adjacency[vertex]) {
        if (visitMarks[neighbor] !== visitStamp) {
          visitMarks[neighbor] = visitStamp;
          stack.push(neighbor);
        }
      }
    }

    return false;
  }

  function tryAddEdge(edge) {
    if (chosen.length + 1 > n) {
      return { possible: false, completeHamiltonian: false };
    }

    if (degree[edge.u] >= 2 || degree[edge.v] >= 2) {
      return { possible: false, completeHamiltonian: false };
    }

    const createsCycle = areConnected(edge.u, edge.v);
    if (!createsCycle) {
      return { possible: true, completeHamiltonian: false };
    }

    const nextTouchedCount =
      touchedCount +
      (touched[edge.u] ? 0 : 1) +
      (touched[edge.v] ? 0 : 1);
    const completeHamiltonian = chosen.length + 1 === n && nextTouchedCount === n;

    return {
      possible: completeHamiltonian,
      completeHamiltonian
    };
  }

  function pushEdge(edge) {
    const token = {
      edge,
      uWasTouched: touched[edge.u],
      vWasTouched: touched[edge.v]
    };

    chosen.push(edge);

    if (!touched[edge.u]) {
      touched[edge.u] = true;
      touchedCount++;
    }
    if (!touched[edge.v]) {
      touched[edge.v] = true;
      touchedCount++;
    }

    degree[edge.u]++;
    degree[edge.v]++;
    adjacency[edge.u].push(edge.v);
    adjacency[edge.v].push(edge.u);

    return token;
  }

  function removeNeighbor(vertex, neighbor) {
    const neighbors = adjacency[vertex];
    const index = neighbors.indexOf(neighbor);
    if (index !== -1) neighbors.splice(index, 1);
  }

  function popEdge(token) {
    const { edge } = token;

    chosen.pop();

    removeNeighbor(edge.u, edge.v);
    removeNeighbor(edge.v, edge.u);
    degree[edge.u]--;
    degree[edge.v]--;

    if (!token.uWasTouched) {
      touched[edge.u] = false;
      touchedCount--;
    }
    if (!token.vWasTouched) {
      touched[edge.v] = false;
      touchedCount--;
    }
  }

  function currentForcedSetSignature(completeHamiltonian) {
    if (chosen.length === 0) {
      return {
        key: "empty",
        label: "empty tour ensemble",
        kind: "empty"
      };
    }

    if (completeHamiltonian) {
      return {
        key: "cycle|zero-forced",
        label: "Hamiltonian forced zero-edge cycle",
        kind: "cycle"
      };
    }

    const visited = Array(n + 1).fill(false);
    const endpointComponent = Array(n + 1).fill(-1);
    const internalVertex = Array(n + 1).fill(false);
    let ceCount = 0;

    for (let start = 1; start <= n; start++) {
      if (!touched[start] || visited[start]) continue;

      const vertices = [];
      const stack = [start];
      visited[start] = true;

      while (stack.length > 0) {
        const vertex = stack.pop();
        vertices.push(vertex);

        for (const neighbor of adjacency[vertex]) {
          if (!visited[neighbor]) {
            visited[neighbor] = true;
            stack.push(neighbor);
          }
        }
      }

      for (const vertex of vertices) {
        if (degree[vertex] === 1) {
          endpointComponent[vertex] = ceCount;
        } else {
          internalVertex[vertex] = true;
        }
      }

      ceCount++;
    }

    const denominator = n - 1 + ceCount - touchedCount;
    if (denominator < 0) {
      return {
        key: "zero|negative-denominator",
        label: "impossible denominator",
        kind: "zero"
      };
    }

    if (denominator === 0) {
      const endpoints = [];
      for (let vertex = 1; vertex <= n; vertex++) {
        if (endpointComponent[vertex] !== -1) endpoints.push(vertex);
      }

      const closingWeight = endpoints.length === 2
        ? BigInt(weight[endpoints[0]][endpoints[1]])
        : 0n;
      const key = [
        "path-close",
        `ce=${ceCount}`,
        `touch=${touchedCount}`,
        `close=${closingWeight}`
      ].join("|");

      return {
        key,
        label: `ce=${ceCount}, touched=${touchedCount}, D=0, close=${closingWeight}`,
        kind: "path-close",
        ceCount,
        touchedCount,
        denominator,
        closingWeight
      };
    }

    let cceSum = 0n;
    let neSum = 0n;
    let ieSum = 0n;

    for (const realEdge of realEdges) {
      const i = realEdge.u;
      const j = realEdge.v;
      if (internalVertex[i] || internalVertex[j]) continue;

      const iComponent = endpointComponent[i];
      const jComponent = endpointComponent[j];
      const iIsEndpoint = iComponent !== -1;
      const jIsEndpoint = jComponent !== -1;

      if (iIsEndpoint && jIsEndpoint) {
        if (iComponent !== jComponent) cceSum -= 1n;
      } else if (iIsEndpoint || jIsEndpoint) {
        neSum -= 1n;
      } else if (!touched[i] && !touched[j]) {
        ieSum -= 1n;
      }
    }

    const key = [
      "paths",
      `ce=${ceCount}`,
      `touch=${touchedCount}`,
      `D=${denominator}`,
      `cce=${cceSum}`,
      `ne=${neSum}`,
      `ie=${ieSum}`
    ].join("|");

    return {
      key,
      label: `ce=${ceCount}, touched=${touchedCount}, D=${denominator}, cce=${cceSum}, ne=${neSum}, ie=${ieSum}`,
      kind: "paths",
      ceCount,
      touchedCount,
      denominator,
      cceSum,
      neSum,
      ieSum
    };
  }

  function addCurrentStats(completeHamiltonian = false) {
    addSignature(currentForcedSetSignature(completeHamiltonian), chosen.length);
  }

  addCurrentStats();

  if (orderedZeroEdges.length === 0) {
    return buildResult(0, "full");
  }

  for (let newEdgeIndex = 0; newEdgeIndex < maxPrefix; newEdgeIndex++) {
    function visitPreviousEdges(previousIndex, completeHamiltonian) {
      if (completeHamiltonian) {
        addCurrentStats(true);
        const skippedSupersets = countCompletions(newEdgeIndex - previousIndex) - 1n;
        skipZeroTerms(skippedSupersets);
        return;
      }

      if (previousIndex === newEdgeIndex) {
        addCurrentStats();
        return;
      }

      visitPreviousEdges(previousIndex + 1, false);

      const edge = orderedZeroEdges[previousIndex];
      const attempt = tryAddEdge(edge);
      if (!attempt.possible) {
        skipZeroTerms(countCompletions(newEdgeIndex - previousIndex - 1));
        return;
      }

      const token = pushEdge(edge);
      visitPreviousEdges(previousIndex + 1, attempt.completeHamiltonian);
      popEdge(token);
    }

    const firstAttempt = tryAddEdge(orderedZeroEdges[newEdgeIndex]);
    if (!firstAttempt.possible) {
      skipZeroTerms(countCompletions(newEdgeIndex));
    } else {
      const token = pushEdge(orderedZeroEdges[newEdgeIndex]);
      visitPreviousEdges(0, firstAttempt.completeHamiltonian);
      popEdge(token);
    }

    const usedZeroEdges = newEdgeIndex + 1;
    const candidate = buildResult(usedZeroEdges, "prefix");
    if (certifiesHamiltonianAverage(n, candidate)) {
      return buildResult(usedZeroEdges, "average");
    }
    if (candidate.omegaValid === 0n) {
      return buildResult(usedZeroEdges, "empty");
    }
  }

  return buildResult(maxPrefix, maxPrefix === orderedZeroEdges.length ? "full" : "cap");
}

function find(parent, x) {
  if (parent[x] !== x) parent[x] = find(parent, parent[x]);
  return parent[x];
}

function union(parent, a, b) {
  const rootA = find(parent, a);
  const rootB = find(parent, b);
  if (rootA !== rootB) parent[rootB] = rootA;
}

function canonicalEdges(edges) {
  const unique = new Map();
  for (const edge of edges) {
    const u = Math.min(edge.u, edge.v);
    const v = Math.max(edge.u, edge.v);
    unique.set(edgeKey(u, v), { u, v });
  }
  return Array.from(unique.values()).sort((a, b) => a.u === b.u ? a.v - b.v : a.u - b.u);
}

function zeroStats() {
  return {
    omega: 0n,
    average: Rational.zero(),
    total: Rational.zero(),
    impossible: true
  };
}

function analyzeForcedSetPrefix(n, forcedEdges) {
  const edges = canonicalEdges(forcedEdges);
  if (edges.length === 0) return { impossible: false, completeHamiltonian: false };

  const parent = Array.from({ length: n + 1 }, (_, index) => index);
  const degree = Array(n + 1).fill(0);
  const touched = new Set();
  let cycleDetected = false;

  for (const edge of edges) {
    const { u, v } = edge;
    touched.add(u);
    touched.add(v);
    degree[u]++;
    degree[v]++;
    if (degree[u] > 2 || degree[v] > 2) {
      return { impossible: true, completeHamiltonian: false };
    }

    const rootU = find(parent, u);
    const rootV = find(parent, v);
    if (rootU === rootV) {
      cycleDetected = true;
    } else {
      union(parent, u, v);
    }
  }

  if (!cycleDetected) return { impossible: false, completeHamiltonian: false };

  const roots = new Set(Array.from(touched, vertex => find(parent, vertex)));
  const allDegreeTwo = Array.from(touched).every(vertex => degree[vertex] === 2);
  const completeHamiltonian =
    touched.size === n &&
    edges.length === n &&
    roots.size === 1 &&
    allDegreeTwo;

  return {
    impossible: !completeHamiltonian,
    completeHamiltonian
  };
}

function forcedSetStatsExact(n, weight, forcedEdges) {
  const edges = canonicalEdges(forcedEdges);
  if (edges.length === 0) {
    const omega = omegaAllTours(n);
    const average = allTourAverageExact(n, weight);
    return { omega, average, total: average.mulInt(omega), impossible: false };
  }

  const parent = Array.from({ length: n + 1 }, (_, index) => index);
  const degree = Array(n + 1).fill(0);
  const touched = new Set();
  const chosenKeys = new Set();
  let ceSum = 0n;
  let cycleDetected = false;

  for (const edge of edges) {
    const { u, v } = edge;
    chosenKeys.add(edgeKey(u, v));
    touched.add(u);
    touched.add(v);
    degree[u]++;
    degree[v]++;
    if (degree[u] > 2 || degree[v] > 2) return zeroStats();

    const rootU = find(parent, u);
    const rootV = find(parent, v);
    if (rootU === rootV) {
      cycleDetected = true;
    } else {
      union(parent, u, v);
    }
    ceSum += BigInt(weight[u][v]);
  }

  if (cycleDetected) {
    const roots = new Set(Array.from(touched, vertex => find(parent, vertex)));
    const allDegreeTwo = Array.from(touched).every(vertex => degree[vertex] === 2);
    const isHamiltonianCycle = touched.size === n && edges.length === n && roots.size === 1 && allDegreeTwo;
    if (!isHamiltonianCycle) return zeroStats();

    const omega = 1n;
    const average = new Rational(ceSum);
    return { omega, average, total: average.mulInt(omega), impossible: false };
  }

  const components = new Map();
  for (const vertex of touched) {
    const root = find(parent, vertex);
    if (!components.has(root)) components.set(root, []);
    components.get(root).push(vertex);
  }

  const endpointToComponent = new Map();
  const internalVertices = new Set();
  for (const [componentId, vertices] of components.entries()) {
    for (const vertex of vertices) {
      if (degree[vertex] === 1) {
        endpointToComponent.set(vertex, componentId);
      } else {
        internalVertices.add(vertex);
      }
    }
  }

  const ceCount = components.size;
  const vceCount = touched.size;
  const denominator = n - 1 + ceCount - vceCount;
  if (denominator < 0) return zeroStats();

  if (denominator === 0) {
    const endpoints = Array.from(endpointToComponent.keys());
    const closingWeight = endpoints.length === 2 ? BigInt(weight[endpoints[0]][endpoints[1]]) : 0n;
    const omega = 2n ** BigInt(ceCount - 1);
    const average = new Rational(ceSum + closingWeight);
    return { omega, average, total: average.mulInt(omega), impossible: false };
  }

  let cceSum = 0n;
  let neSum = 0n;
  let ieSum = 0n;

  for (let i = 1; i < n; i++) {
    for (let j = i + 1; j <= n; j++) {
      if (chosenKeys.has(edgeKey(i, j))) continue;
      if (internalVertices.has(i) || internalVertices.has(j)) continue;

      const iComponent = endpointToComponent.get(i);
      const jComponent = endpointToComponent.get(j);
      const iIsEndpoint = iComponent !== undefined;
      const jIsEndpoint = jComponent !== undefined;
      const edgeWeight = BigInt(weight[i][j]);

      if (iIsEndpoint && jIsEndpoint) {
        if (iComponent !== jComponent) cceSum += edgeWeight;
      } else if (iIsEndpoint || jIsEndpoint) {
        neSum += edgeWeight;
      } else if (!touched.has(i) && !touched.has(j)) {
        ieSum += edgeWeight;
      }
    }
  }

  const d = BigInt(denominator);
  const average = new Rational(ceSum)
    .add(new Rational(cceSum, 2n * d))
    .add(new Rational(neSum, d))
    .add(new Rational(2n * ieSum, d));
  const omega = (2n ** BigInt(ceCount - 1)) * factorialBigInt(denominator);

  return {
    omega,
    average,
    total: average.mulInt(omega),
    impossible: false
  };
}

function runExactSubtraction(n, weight, zeroEdges) {
  let termsVisited = 0;
  let nonzeroTerms = 0;
  let zeroOmegaTerms = 0;
  let prunedTerms = 0;
  const chosen = [];
  const groups = new Map();

  function addStatsToGroups(stats, subsetSize) {
    termsVisited++;
    if (stats.omega === 0n) {
      zeroOmegaTerms++;
      return;
    }

    nonzeroTerms++;
    const sign = subsetSize % 2 === 0 ? 1n : -1n;
    const key = `${stats.omega}|${stats.average.n}/${stats.average.d}`;
    let group = groups.get(key);
    if (!group) {
      group = {
        coefficient: 0n,
        omega: stats.omega,
        average: stats.average,
        total: stats.total,
        subsets: 0
      };
      groups.set(key, group);
    }
    group.coefficient += sign;
    group.subsets++;
  }

  function countCompletions(index) {
    return Number(2n ** BigInt(zeroEdges.length - index));
  }

  function visit(index) {
    if (chosen.length > n) {
      const skippedTerms = countCompletions(index);
      termsVisited += skippedTerms;
      zeroOmegaTerms += skippedTerms;
      prunedTerms += skippedTerms;
      return;
    }

    const prefix = analyzeForcedSetPrefix(n, chosen);
    if (prefix.impossible) {
      const skippedTerms = countCompletions(index);
      termsVisited += skippedTerms;
      zeroOmegaTerms += skippedTerms;
      prunedTerms += skippedTerms;
      return;
    }

    if (prefix.completeHamiltonian) {
      const stats = forcedSetStatsExact(n, weight, chosen);
      addStatsToGroups(stats, chosen.length);
      const skippedSupersets = countCompletions(index) - 1;
      termsVisited += skippedSupersets;
      zeroOmegaTerms += skippedSupersets;
      prunedTerms += skippedSupersets;
      return;
    }

    if (index === zeroEdges.length) {
      const stats = forcedSetStatsExact(n, weight, chosen);
      addStatsToGroups(stats, chosen.length);
      return;
    }

    visit(index + 1);
    chosen.push(zeroEdges[index]);
    visit(index + 1);
    chosen.pop();
  }

  visit(0);

  let omegaValid = 0n;
  let totalValid = Rational.zero();
  let cancelledGroups = 0;
  const activeGroups = [];

  for (const group of groups.values()) {
    if (group.coefficient === 0n) {
      cancelledGroups++;
      continue;
    }
    omegaValid += group.coefficient * group.omega;
    totalValid = totalValid.add(group.total.mulInt(group.coefficient));
    activeGroups.push(group);
  }

  activeGroups.sort((a, b) => {
    const absA = absBigInt(a.coefficient);
    const absB = absBigInt(b.coefficient);
    if (absA === absB) return b.subsets - a.subsets;
    return absA > absB ? -1 : 1;
  });

  const averageValid = omegaValid !== 0n
    ? totalValid.divInt(omegaValid)
    : null;

  return {
    omegaValid,
    totalValid,
    averageValid,
    termsVisited,
    nonzeroTerms,
    zeroOmegaTerms,
    prunedTerms,
    groupedTerms: groups.size,
    activeGroupedTerms: activeGroups.length,
    cancelledGroups,
    activeGroups
  };
}

function formatBigInt(value) {
  return value.toString();
}

function formatRational(value) {
  return value ? value.toDisplay() : "n/a";
}

function renderTermLabels(zeroEdges) {
  if (zeroEdges.length === 0) {
    return `
      <div class="term-labels">
        <div class="term-title">Inclusion-exclusion labels</div>
        <div class="term-note">No zero edges, so there are no correction labels.</div>
      </div>
    `;
  }

  const edgeDictionary = zeroEdges
    .map((edge, index) => `<span>e${index + 1}=(${edge.u},${edge.v})</span>`)
    .join("");
  const termCount = (2n ** BigInt(zeroEdges.length)) - 1n;

  if (termCount > BigInt(MAX_VISIBLE_LABEL_TERMS)) {
    return `
      <div class="term-labels">
        <div class="term-title">Inclusion-exclusion labels</div>
        <div class="edge-dictionary">${edgeDictionary}</div>
        <div class="term-note">
          There are ${termCount.toString()} non-empty unique subset labels. The exact arithmetic ran, but the label display is capped at ${MAX_VISIBLE_LABEL_TERMS} terms so the browser stays responsive.
        </div>
      </div>
    `;
  }

  const lines = [];
  const combo = [];

  function buildCombinations(size, start) {
    if (combo.length === size) {
      const sign = size % 2 === 0 ? "+" : "-";
      lines[lines.length - 1].push(`${sign}${combo.map(index => `e${index + 1}`).join("")}`);
      return;
    }

    for (let index = start; index < zeroEdges.length; index++) {
      combo.push(index);
      buildCombinations(size, index + 1);
      combo.pop();
    }
  }

  for (let size = 1; size <= zeroEdges.length; size++) {
    lines.push([]);
    buildCombinations(size, 0);
  }

  const expansion = lines
    .map((terms, index) => `size ${index + 1}: ${terms.join(" ")}`)
    .join("\n");

  return `
    <div class="term-labels">
      <div class="term-title">Inclusion-exclusion labels</div>
      <div class="edge-dictionary">${edgeDictionary}</div>
      <pre class="term-list">${expansion}</pre>
    </div>
  `;
}

function renderGroupedSummary(result) {
  const hiddenGroups = Math.max(0, result.activeGroups.length - MAX_VISIBLE_GROUPS);
  const rows = result.activeGroups.slice(0, MAX_VISIBLE_GROUPS).map((group, index) => `
    <tr>
      <td>${index + 1}</td>
      <td>${group.coefficient.toString()}</td>
      <td>${group.subsets}</td>
      <td>${formatBigInt(group.omega)}</td>
      <td>${formatRational(group.average)}</td>
      <td>${group.signature || "same Omega/Average"}</td>
    </tr>
  `).join("");

  return `
    <div class="group-summary">
      <div class="term-title">Exact grouped terms</div>
      <div class="term-note">
        Groups combine subsets with the same exact Equation 3 signature. Forced CE weights are zero, so only real -1 edges affect the cce/ne/ie sums.
      </div>
      <div class="table-wrap compact">
        <table>
          <thead>
            <tr>
              <th>#</th>
              <th>Signed coefficient</th>
              <th>Subsets grouped</th>
              <th>Omega</th>
              <th>Average</th>
              <th>Equation 3 signature</th>
            </tr>
          </thead>
          <tbody>${rows}</tbody>
        </table>
      </div>
      ${hiddenGroups ? `<div class="term-note">${hiddenGroups} more active groups are hidden.</div>` : ""}
    </div>
  `;
}

function renderEraseOrderSummary(result) {
  const orderedZeroEdges = result.orderedZeroEdges || [];
  const usedZeroEdges = result.usedZeroEdges ?? orderedZeroEdges.length;
  if (orderedZeroEdges.length === 0) {
    return `<div class="term-note">No zero edges need to be erased.</div>`;
  }

  const shownEdges = orderedZeroEdges
    .slice(0, usedZeroEdges)
    .map((edge, index) => `<span>e${index + 1}=(${edge.u},${edge.v})</span>`)
    .join("");

  const unusedCount = Math.max(0, orderedZeroEdges.length - usedZeroEdges);
  return `
    <div class="term-note">
      Zero edges are ordered to share vertices early, so forced unions hit degree/cycle impossibilities sooner.
    </div>
    <div class="edge-dictionary">${shownEdges}</div>
    ${unusedCount ? `<div class="term-note">${unusedCount} later zero edges were not needed for this exact result.</div>` : ""}
  `;
}

function renderSummary(n, realEdgeCount, zeroEdgeCount, allAverage, subsetCountText) {
  summary.className = "summary";
  summary.innerHTML = `
    <div class="stat"><span>Vertices inferred</span><strong>${n}</strong></div>
    <div class="stat"><span>Real edges</span><strong>${realEdgeCount}</strong></div>
    <div class="stat"><span>Zero edges</span><strong>${zeroEdgeCount}</strong></div>
    <div class="stat"><span>Subsets</span><strong>${subsetCountText}</strong></div>
    <div class="stat"><span>Average of all tours</span><strong>${formatRational(allAverage)}</strong></div>
  `;
}

function renderExactResult(n, zeroEdges, result) {
  const usedZeroEdges = result.usedZeroEdges ?? zeroEdges.length;
  const totalZeroEdges = result.totalZeroEdges ?? zeroEdges.length;
  const zeroEdgeTryCap = result.zeroEdgeTryCap ?? DEFAULT_ZERO_EDGES_TO_TRY;
  const forcedSubsetDepth = result.forcedSubsetDepth ?? n;
  const hasTours = result.omegaValid > 0n;
  const isHamiltonian = certifiesHamiltonianAverage(n, result);
  const hcTourCount = isHamiltonian ? formatBigInt(result.omegaValid) : "n/a";
  const message = !hasTours
    ? "NO HAMILTONIAN CYCLE"
    : (
      isHamiltonian
        ? `HAMILTONIAN CYCLE EXISTS (${hcTourCount} tours)`
        : (
          result.prefixLimited
            ? "PARTIAL EXACT SEARCH: ZERO-EDGE CAP REACHED"
            : "SUBTRACTION COMPLETE, CHECK AVERAGE"
        )
    );
  const searchStatus = !hasTours
    ? "Stopped: Omega after is 0"
    : (
      isHamiltonian
        ? (usedZeroEdges < totalZeroEdges ? "Stopped early: Average after is -n" : "Certified: Average after is -n")
        : (result.prefixLimited ? `Stopped at ${zeroEdgeTryCap} tried zero edges` : "Full erase set processed")
    );

  subtractionResult.className = "subtraction-card";
  subtractionResult.innerHTML = `
    <div class="subtraction-title">Exact zero-edge subtraction</div>
    <div class="subtraction-grid">
      <div><span>Zero edges used</span><strong>${usedZeroEdges}</strong></div>
      <div><span>Total zero edges</span><strong>${totalZeroEdges}</strong></div>
      <div><span>Zero-edge try cap</span><strong>${zeroEdgeTryCap}</strong></div>
      <div><span>Forced subset depth</span><strong>n = ${forcedSubsetDepth}</strong></div>
      <div><span>Search status</span><strong>${searchStatus}</strong></div>
      <div><span>Terms visited</span><strong>${result.termsVisited}</strong></div>
      <div><span>Nonzero terms</span><strong>${result.nonzeroTerms}</strong></div>
      <div><span>Impossible terms</span><strong>${result.zeroOmegaTerms}</strong></div>
      <div><span>Pruned zero terms</span><strong>${result.prunedTerms}</strong></div>
      <div><span>Grouped terms</span><strong>${result.groupedTerms}</strong></div>
      <div><span>Equation 3 structures</span><strong>${result.structuralGroups}</strong></div>
      <div><span>Structure reuses</span><strong>${result.structuralCacheHits}</strong></div>
      <div><span>Active groups</span><strong>${result.activeGroupedTerms}</strong></div>
      <div><span>Cancelled groups</span><strong>${result.cancelledGroups}</strong></div>
      <div><span>Omega after</span><strong>${formatBigInt(result.omegaValid)}</strong></div>
      <div><span>HC tours inferred</span><strong>${hcTourCount}</strong></div>
      <div><span>Total after</span><strong>${formatRational(result.totalValid)}</strong></div>
      <div><span>Average after</span><strong>${formatRational(result.averageValid)}</strong></div>
      <div><span>Result</span><strong>${message}</strong></div>
    </div>
    <div class="subtraction-formula">
      Omega after = sum (-1)^|S| Omega(S)<br>
      Total after = sum (-1)^|S| Omega(S) * Average(S)<br>
      Average after = Total after / Omega after
    </div>
    ${renderEraseOrderSummary(result)}
    ${renderGroupedSummary(result)}
  `;
}

function showError(message) {
  summary.className = "summary error";
  summary.textContent = message;
  subtractionResult.className = "subtraction-card empty";
  subtractionResult.textContent = "The exact inclusion-exclusion result will appear here.";
}

function calculate() {
  try {
    const { edges, n, seen } = parseEdges(edgeInput.value);
    if (n < 3) {
      showError("Enter enough edges to infer at least 3 vertices.");
      return;
    }

    const weight = buildWeight(n, edges);
    const zeroEdges = zeroEdgesFromGraph(n, seen);
    const zeroEdgeTryCap = readZeroEdgeTryCap();
    const allAverage = allTourAverageExact(n, weight);
    const subsetCount = 2n ** BigInt(zeroEdges.length);

    renderSummary(n, edges.length, zeroEdges.length, allAverage, subsetCount.toString());

    const result = runSmartSubtractionSearch(n, weight, zeroEdges, zeroEdgeTryCap);
    renderExactResult(n, zeroEdges, result);
  } catch (error) {
    showError(error.message);
  }
}

calculateBtn.addEventListener("click", calculate);
sampleBtn.addEventListener("click", () => {
  edgeInput.value = "1 2\n2 3\n3 4\n4 1\n1 3";
  calculate();
});
clearBtn.addEventListener("click", () => {
  edgeInput.value = "";
  summary.className = "summary empty";
  summary.textContent = "Enter edges and run exact subtraction.";
  subtractionResult.className = "subtraction-card empty";
  subtractionResult.textContent = "The exact inclusion-exclusion result will appear here.";
});
