function buildSpinGlass50Text() {
  let seed = 501337;
  const rand = () => {
    seed = (1664525 * seed + 1013904223) >>> 0;
    return seed / 2 ** 32;
  };
  const lines = ["N 50", "", "FIELDS"];
  for (let i = 1; i <= 50; i++) {
    const h = Math.round((rand() * 1.2 - 0.6) * 100) / 100;
    lines.push(`${i} ${h}`);
  }
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

export const examples = {
  spinGlass50: {
    name: "50-node spin glass",
    text: buildSpinGlass50Text()
  },
  ferro2: {
    name: "Two-spin ferromagnet",
    text: "N 2\n\nFIELDS\n1 0\n2 0\n\nCOUPLINGS\n1 2 1"
  },
  antiferro2: {
    name: "Two-spin antiferromagnet",
    text: "N 2\n\nFIELDS\n1 0\n2 0\n\nCOUPLINGS\n1 2 -1"
  },
  triangleFerro: {
    name: "Triangle ferromagnet",
    text: "N 3\n\nFIELDS\n1 0\n2 0\n3 0\n\nCOUPLINGS\n1 2 1\n2 3 1\n1 3 1"
  },
  frustratedTriangle: {
    name: "Frustrated triangle",
    text: "N 3\n\nFIELDS\n1 0\n2 0\n3 0\n\nCOUPLINGS\n1 2 1\n2 3 1\n1 3 -1"
  },
  fields5: {
    name: "Five-spin field example",
    text: "N 5\n\nFIELDS\n1 0.5\n2 -0.2\n3 0\n4 1.0\n5 -0.4\n\nCOUPLINGS\n1 2 -1\n1 3 0.5\n2 3 1\n2 4 -0.8\n3 5 1.2\n4 5 -1"
  }
};
