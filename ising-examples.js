function buildSpinGlass100Text() {
  let seed = 501337;
  const rand = () => {
    seed = (1664525 * seed + 1013904223) >>> 0;
    return seed / 2 ** 32;
  };
  const n = 100;
  const targetEdges = 240;
  const lines = [`N ${n}`, "", "FIELDS"];
  for (let i = 1; i <= n; i++) {
    const h = Math.round((rand() * 1.2 - 0.6) * 100) / 100;
    lines.push(`${i} ${h}`);
  }
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

export const examples = {
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
