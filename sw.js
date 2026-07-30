/* LAMPrime — service worker (offline-first PWA).
   Caches the app shell so LAMPrime works with no network after the first visit. */
const CACHE = 'lamprime-v2.0.0';
const CORE = [
  './', './index.html', './LAMPrime.html', './LAMPrime_en.html',
  './app.js', './styles.css', './install_LAMPrime.html', './manifest.webmanifest',
  './icons/icon-192.png', './icons/icon-512.png',
  './icons/apple-touch-icon-180.png', './icons/favicon-32.png'
];

self.addEventListener('install', (e) => {
  self.skipWaiting();
  e.waitUntil(caches.open(CACHE).then((c) => c.addAll(CORE)).catch(() => {}));
});

self.addEventListener('activate', (e) => {
  e.waitUntil(
    caches.keys()
      .then((keys) => Promise.all(keys.filter((k) => k !== CACHE).map((k) => caches.delete(k))))
      .then(() => self.clients.claim())
  );
});

self.addEventListener('fetch', (e) => {
  const req = e.request;
  if (req.method !== 'GET') return;
  e.respondWith(
    caches.match(req).then((hit) => hit || fetch(req).then((res) => {
      const copy = res.clone();
      caches.open(CACHE).then((c) => { try { c.put(req, copy); } catch (_) {} });
      return res;
    }).catch(() => caches.match('./LAMPrime_en.html')))
  );
});
