import {
  approximateLnZ,
  calculateCumulants,
  formatNumber,
  parseIsingText
} from "./ising-math.js";
import { solveTaylor } from "./ising-search.js";
import { examples } from "./ising-examples.js";

const $ = (id) => document.getElementById(id);
const output = $("output");
const summaryGrid = $("summaryGrid");

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
  output.textContent = text;
}

function setSummary(metrics) {
  summaryGrid.innerHTML = metrics.map(([label, value]) => (
    `<div class="metric"><span>${escapeHtml(label)}</span><strong>${escapeHtml(String(value))}</strong></div>`
  )).join("");
}

function currentModel() {
  return parseIsingText($("isingInput").value);
}

function modelStats(model) {
  return {
    fields: model.fields.filter((h) => Math.abs(h) > 1e-12).length,
    couplings: model.couplings.size
  };
}

function runSolver() {
  try {
    const model = currentModel();
    const opts = options();
    const cumulants = calculateCumulants(model, opts.order);
    const lnz = approximateLnZ(model, opts.beta, opts.order);
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
      `Taylor order: ${opts.order}`,
      "",
      "Cumulants",
      ...Array.from({ length: opts.order }, (_, i) => i + 1).map((r) => `kappa${r}: ${formatNumber(cumulants.kappa[r])}`),
      `Cycles: triangles=${cumulants.counts.triangles}, C4=${cumulants.counts.cycles4}, C5=${cumulants.counts.cycles5}`,
      `ln Z approximation: ${formatNumber(lnz.value)}`,
      "",
      "Search result",
      `Final spins: ${formatSpins(result.spins)}`,
      `Exact final energy: ${formatNumber(result.energy)}`,
      `Best energy found: ${formatNumber(result.bestEnergy)}`,
      `Decisions: ${result.decisions}`,
      `Backtracks: ${result.backtracks}`,
      `Runtime: ${formatNumber(result.runtimeMs)} ms`
    ];

    lines.push("", "Precheck summary", ...formatPrecheckSummary(result.precheckSummary));

    lines.push("", "Decision trace", ...formatDecisionLog(result.decisionLog));
    write(lines.join("\n"));
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
    setSummary([["benchmark", "beta sweep"], ["order", opts.order], ["mode", opts.mode], ["cases", betas.length]]);
    write(rows.join("\n"));
  } catch (error) {
    write(`Error: ${error.message}`);
  }
}

function formatDecisionLog(log) {
  if (!log.length) return ["No decisions recorded."];
  return log.flatMap((step) => {
    if (step.type === "precheck") {
      const updates = step.fieldUpdates.length
        ? step.fieldUpdates.map((u) => `  h${u.spin} += ${formatNumber(u.delta)}`).join("\n")
        : "  none";
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
    const candidates = step.candidateScores
      .map((score) => `  s${score.originalSpin + 1}: I=${formatNumber(score.importance)} preferred=${score.preferred > 0 ? "+1" : "-1"}`)
      .join("\n");
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
  const entries = [
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
    ["Tree components solved", summary.treeComponentsSolved],
    ["Bipartite components detected", summary.bipartiteComponents],
    ["Unfrustrated gauge components detected", summary.unfrustratedGaugeComponents],
    ["Remaining unresolved spins", summary.remainingSpins]
  ];
  return entries.map(([label, value]) => `${label}: ${value ?? 0}`);
}

function formatSpins(spins) {
  return `[${spins.map((s) => s > 0 ? "+1" : "-1").join(" ")}]`;
}

function escapeHtml(value) {
  return value.replace(/[&<>"']/g, (ch) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[ch]));
}

function init() {
  $("startupHelp").style.display = "none";
  for (const [key, example] of Object.entries(examples)) {
    const option = document.createElement("option");
    option.value = key;
    option.textContent = example.name;
    $("exampleSelect").append(option);

    const button = document.createElement("button");
    button.type = "button";
    button.className = "miniButton";
    button.textContent = example.name;
    button.addEventListener("click", () => {
      $("exampleSelect").value = key;
      $("isingInput").value = example.text;
      runSolver();
    });
    $("exampleButtons").append(button);
  }
  $("isingInput").value = examples.fields5.text;
  $("isingInput").value = examples.spinGlass50.text;

  document.querySelectorAll(".tab").forEach((tab) => {
    tab.addEventListener("click", () => {
      document.querySelectorAll(".tab, .panel").forEach((node) => node.classList.remove("active"));
      tab.classList.add("active");
      $(tab.dataset.panel).classList.add("active");
      if (tab.dataset.panel === "visualizer") drawVisualization();
    });
  });

  $("loadExample").addEventListener("click", () => {
    $("isingInput").value = examples[$("exampleSelect").value].text;
    prepareVisualization(currentModel(), null);
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
  $("visualReset").addEventListener("click", resetVisualization);
  $("spinCanvas").addEventListener("click", selectCanvasSpin);
  window.addEventListener("resize", drawVisualization);

  prepareVisualization(currentModel(), null);
  setSummary([["spins", 50], ["example", "50-node spin glass"], ["visual", "ready"]]);
  write("Ready. The 50-node spin glass is loaded. Open the Visualizer tab to inspect it, or press Run Taylor solver when you want to compute decisions.");
}

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
    const radius = 210 + 48 * Math.sin(i * 1.7);
    return {
      index: i,
      x: 450 + Math.cos(angle) * radius,
      y: 310 + Math.sin(angle) * radius,
      vx: 0,
      vy: 0
    };
  });
  const rendered = new Set(nodes.map((n) => n.index));
  const edges = [...model.couplings.entries()].map(([key, J]) => {
    const [a, b] = key.split(":").map(Number);
    return { a, b, J };
  }).filter((e) => rendered.has(e.a) && rendered.has(e.b));
  visualState.nodes = nodes;
  visualState.edges = edges;
  settleLayout();
  drawVisualization();
}

function settleLayout() {
  const nodes = visualState.nodes;
  const edges = visualState.edges;
  for (let tick = 0; tick < 180; tick++) {
    for (const node of nodes) {
      node.vx += (450 - node.x) * 0.0008;
      node.vy += (310 - node.y) * 0.0008;
    }
    for (let i = 0; i < nodes.length; i++) {
      for (let j = i + 1; j < nodes.length; j++) {
        const a = nodes[i];
        const b = nodes[j];
        const dx = a.x - b.x;
        const dy = a.y - b.y;
        const d2 = Math.max(64, dx * dx + dy * dy);
        const force = 56 / d2;
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
      const force = (dist - 92) * 0.003;
      a.vx += dx / dist * force;
      a.vy += dy / dist * force;
      b.vx -= dx / dist * force;
      b.vy -= dy / dist * force;
    }
    for (const node of nodes) {
      node.vx *= 0.86;
      node.vy *= 0.86;
      node.x = Math.min(870, Math.max(30, node.x + node.vx));
      node.y = Math.min(590, Math.max(30, node.y + node.vy));
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
  canvas.height = Math.max(420, Math.floor(620 * scale));
  ctx.setTransform(canvas.width / 900, 0, 0, canvas.height / 620, 0, 0);
  ctx.clearRect(0, 0, 900, 620);

  const maxJ = Math.max(0.01, ...visualState.edges.map((e) => Math.abs(e.J)));
  for (const edge of visualState.edges) {
    const a = visualState.nodes[edge.a];
    const b = visualState.nodes[edge.b];
    if (!a || !b) continue;
    const resolved = visualState.assignments.has(edge.a) || visualState.assignments.has(edge.b);
    ctx.globalAlpha = resolved ? 0.22 : 0.72;
    ctx.strokeStyle = edge.J > 0 ? "#36d3c2" : "#f4b44f";
    ctx.lineWidth = 0.75 + 4 * Math.abs(edge.J) / maxJ;
    if (edge.J < 0) ctx.setLineDash([7, 5]);
    else ctx.setLineDash([]);
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
    ctx.fillStyle = nodeColor(value, reason);
    ctx.strokeStyle = visualState.selected === node.index ? "#edf6f5" : "#283443";
    ctx.lineWidth = visualState.selected === node.index ? 3 : 1.5;
    ctx.beginPath();
    ctx.arc(node.x, node.y, 12, 0, Math.PI * 2);
    ctx.fill();
    ctx.stroke();
    ctx.fillStyle = "#edf6f5";
    ctx.font = "10px system-ui";
    ctx.textAlign = "center";
    ctx.fillText(String(node.index + 1), node.x, node.y - 17);
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
    box.textContent = `Showing ${visualState.nodes.length} spins and ${visualState.edges.length} couplings. Select a spin to inspect it.`;
    return;
  }
  const i = visualState.selected;
  const degree = visualState.edges.filter((e) => e.a === i || e.b === i).length;
  const value = visualState.assignments.has(i) ? (visualState.assignments.get(i) > 0 ? "+1" : "-1") : "unresolved";
  const reason = visualState.reasons.get(i) || "UNRESOLVED";
  const importance = visualState.importance.has(i) ? formatNumber(visualState.importance.get(i)) : "n/a";
  box.textContent = `Spin ${i + 1} | value ${value} | h = ${formatNumber(model.fields[i] || 0)} | degree ${degree} | reason ${reason} | I = ${importance}`;
}

init();
