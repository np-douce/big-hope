const CACHE = "np-douce-ising-lab-v17";
const ASSETS = [
  "index.html",
  "ising.html",
  "styles.css",
  "ising-app.js",
  "ising.js",
  "ising-examples.js",
  "ising-math.js",
  "ising-precheck.js",
  "ising-search.js",
  "ising-tests.js",
  "manifest.webmanifest",
  "examples/spin-glass-50.txt",
  "examples/ferromagnetic-ring-8.txt",
  "examples/antiferromagnetic-ring-8.txt",
  "examples/frustrated-antiferro-ring-9.txt",
  "examples/planted-field-12.txt"
];

self.addEventListener("install", (event) => {
  event.waitUntil(caches.open(CACHE).then((cache) => cache.addAll(ASSETS)));
});

self.addEventListener("activate", (event) => {
  event.waitUntil(caches.keys().then((keys) => Promise.all(keys.filter((key) => key !== CACHE).map((key) => caches.delete(key)))));
});

self.addEventListener("fetch", (event) => {
  if (event.request.method !== "GET") return;
  event.respondWith(caches.match(event.request).then((cached) => cached || fetch(event.request)));
});
