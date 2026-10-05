// Radar replay: reflectivity frames with flight tracks; seeding aircraft's air is carried
// downwind with its drift so it can be seen arriving on its snow streak.
(function () {
  const data = JSON.parse(document.getElementById("replay-data").textContent);
  const canvas = document.getElementById("replay");
  const ctx = canvas.getContext("2d");
  const slider = document.getElementById("replay-time");
  const play = document.getElementById("replay-play");
  const clock = document.getElementById("replay-clock");
  const listEl = document.getElementById("replay-seeders");
  const half = data.half_km;
  // frames are 8-bit codes (0 = no echo, 1..255 = 0..35 dBZ); colour them per theme
  const frames = data.frames.map((f) => {
    const img = new Image();
    img.src = "data:image/png;base64," + f.png;
    return { t: f.t, img, cache: {} };
  });

  function hexToRgb(h) {
    const v = h.replace("#", "");
    return [parseInt(v.slice(0, 2), 16), parseInt(v.slice(2, 4), 16), parseInt(v.slice(4, 6), 16)];
  }

  function lut() {
    const stops = [1, 2, 3, 4, 5].map((i) => hexToRgb(css("--r" + i)));
    const out = new Uint8ClampedArray(256 * 4);
    for (let v = 1; v < 256; v++) {
      const s = (v - 1) / 254 * (stops.length - 1);
      const i = Math.min(Math.floor(s), stops.length - 2), w = s - i;
      for (let c = 0; c < 3; c++) out[v * 4 + c] = stops[i][c] * (1 - w) + stops[i + 1][c] * w;
      out[v * 4 + 3] = 255;
    }
    return out;
  }

  function coloured(f) {
    const key = css("--r1") + css("--r5");
    if (f.cache[key]) return f.cache[key];
    if (!f.img.complete || !f.img.naturalWidth) return null;
    const c = document.createElement("canvas");
    c.width = f.img.naturalWidth;
    c.height = f.img.naturalHeight;
    const g = c.getContext("2d");
    g.drawImage(f.img, 0, 0);
    const im = g.getImageData(0, 0, c.width, c.height);
    const L = lut();
    for (let p = 0; p < im.data.length; p += 4) {
      const v = im.data[p];
      im.data[p] = L[v * 4]; im.data[p + 1] = L[v * 4 + 1]; im.data[p + 2] = L[v * 4 + 2];
      im.data[p + 3] = v ? 255 : 0;
    }
    g.putImageData(im, 0, 0);
    f.cache = { [key]: c };
    return c;
  }
  const seeders = data.tracks.filter((t) => t.seeder);
  const others = data.tracks.filter((t) => !t.seeder);
  const info = data.seeder_info || {};
  let k = 0, timer = null, hover = null;

  const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

  function size() {
    const r = canvas.getBoundingClientRect();
    const dpr = window.devicePixelRatio || 1;
    canvas.width = Math.round(r.width * dpr);
    canvas.height = Math.round(r.width * dpr);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    draw();
  }

  function px(xkm) { return (xkm + half) / (2 * half) * canvas.clientWidth; }
  function py(ykm) { return (half - ykm) / (2 * half) * canvas.clientWidth; }

  function hhmm(t, offsetH) {
    const d = new Date((t + offsetH * 3600) * 1000);
    return String(d.getUTCHours()).padStart(2, "0") + ":" + String(d.getUTCMinutes()).padStart(2, "0");
  }

  function drawTrackWindow(tr, T, minAge, maxAge, colour, width, alpha) {
    ctx.strokeStyle = colour;
    ctx.lineWidth = width;
    ctx.globalAlpha = alpha;
    ctx.beginPath();
    let open = false;
    for (let i = 0; i < tr.t.length; i++) {
      const age = T - tr.t[i];
      if (age < minAge || age > maxAge) { open = false; continue; }
      const X = px(tr.x[i]), Y = py(tr.y[i]);
      if (!open) { ctx.moveTo(X, Y); open = true; } else { ctx.lineTo(X, Y); }
    }
    ctx.stroke();
    ctx.globalAlpha = 1;
  }

  function drawDrifted(tr, T, colour, width) {
    const s = tr.seed;
    if (!s) return false;
    ctx.strokeStyle = colour;
    ctx.lineWidth = width;
    ctx.lineCap = "round";
    ctx.lineJoin = "round";
    ctx.beginPath();
    let open = false, any = false;
    const alt0 = s.seed_alt_ft;
    for (let i = 0; i < s.t.length; i++) {
      const age = T - s.t[i];
      // only the part of the flight inside the seeding layer, within 1000 ft of the altitude
      // at which it matched its streak
      const inLayer = alt0 == null || Math.abs(s.alt[i] - alt0) <= 1000;
      if (!inLayer || age < 5 * 60 || age > 75 * 60) { open = false; continue; }
      const X = px(s.x[i] + s.u[i] * age / 1000), Y = py(s.y[i] + s.v[i] * age / 1000);
      if (!open) { ctx.moveTo(X, Y); open = true; } else { ctx.lineTo(X, Y); }
      any = true;
    }
    ctx.stroke();
    return any;
  }

  function draw() {
    const W = canvas.clientWidth;
    ctx.fillStyle = css("--radar-bg");
    ctx.fillRect(0, 0, W, W);
    const f = frames[k];
    ctx.imageSmoothingEnabled = false;
    const img = coloured(f);
    if (img) ctx.drawImage(img, 0, 0, W, W);
    else f.img.onload = draw;
    const T = f.t;
    // range rings every 20 km
    ctx.strokeStyle = css("--rule");
    ctx.lineWidth = 1;
    for (let r = 20; r <= half; r += 20) {
      ctx.beginPath();
      ctx.arc(px(0), py(0), r / (2 * half) * W, 0, 2 * Math.PI);
      ctx.stroke();
    }
    for (const tr of others) drawTrackWindow(tr, T, 0, 30 * 60, css("--traffic"), 0.8, 0.55);
    const active = [];
    for (const tr of seeders) {
      drawTrackWindow(tr, T, 0, 75 * 60, css("--seed-path"), 1.2, 0.7);
      const on = drawDrifted(tr, T, css("--seed"), hover === tr.id ? 4 : 2.6);
      if (on) active.push(tr.id);
    }
    // markers: radar and airport
    ctx.fillStyle = css("--ink");
    ctx.beginPath();
    ctx.moveTo(px(0), py(0) - 6); ctx.lineTo(px(0) + 5, py(0) + 4); ctx.lineTo(px(0) - 5, py(0) + 4);
    ctx.closePath(); ctx.fill();
    if (data.airport) {
      const ax = px(data.airport.x), ay = py(data.airport.y);
      ctx.fillRect(ax - 4, ay - 4, 8, 8);
      ctx.font = "600 12px " + css("--font-ui");
      ctx.fillText(data.airport.name, ax + 8, ay - 6);
    }
    ctx.font = "12px " + css("--font-ui");
    ctx.fillStyle = css("--ink-2");
    ctx.fillText("radar", px(0) + 8, py(0) + 14);
    // scale bar 10 km
    const x0 = 14, y0 = W - 16, L = 10 / (2 * half) * W;
    ctx.strokeStyle = css("--ink");
    ctx.lineWidth = 2;
    ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x0 + L, y0); ctx.stroke();
    ctx.fillStyle = css("--ink");
    ctx.fillText("10 km", x0, y0 - 6);
    clock.textContent = hhmm(T, 0) + " UTC · " + hhmm(T, data.local_offset_h) + " " + data.local_tz;
    for (const li of listEl.querySelectorAll("li")) {
      li.classList.toggle("is-active", active.includes(Number(li.dataset.id)));
    }
  }

  function setFrame(i) { k = Math.max(0, Math.min(frames.length - 1, i)); slider.value = k; draw(); }

  function toggle() {
    if (timer) { clearInterval(timer); timer = null; play.textContent = "Play"; play.setAttribute("aria-pressed", "false"); return; }
    if (k >= frames.length - 1) setFrame(0);
    play.textContent = "Pause";
    play.setAttribute("aria-pressed", "true");
    timer = setInterval(() => {
      if (k >= frames.length - 1) { toggle(); return; }
      setFrame(k + 1);
    }, 450);
  }

  // seeder list
  const rows = seeders
    .map((tr) => ({ tr, i: info[tr.id] || {} }))
    .sort((a, b) => (a.i.t_pass || 0) - (b.i.t_pass || 0));
  for (const { tr, i } of rows) {
    const li = document.createElement("li");
    li.dataset.id = tr.id;
    li.tabIndex = 0;
    const temp = i.seed_temp_c != null ? Math.round(i.seed_temp_c) + " °C" : "";
    const alt = i.seed_alt_ft != null ? Math.round(i.seed_alt_ft / 100) * 100 + " ft" : "";
    const when = i.t_pass ? hhmm(i.t_pass, 0) : "";
    li.innerHTML = `<span class="cs">${tr.callsign}</span><span class="ty">${tr.type}</span>` +
      `<span class="meta">${when} UTC · ${alt} · ${temp}</span>`;
    const jump = () => {
      if (!i.t_match) return;
      let best = 0;
      frames.forEach((f, j) => { if (Math.abs(f.t - i.t_match) < Math.abs(frames[best].t - i.t_match)) best = j; });
      setFrame(best);
    };
    li.addEventListener("mouseenter", () => { hover = tr.id; draw(); });
    li.addEventListener("mouseleave", () => { hover = null; draw(); });
    li.addEventListener("focus", () => { hover = tr.id; draw(); });
    li.addEventListener("blur", () => { hover = null; draw(); });
    li.addEventListener("click", jump);
    li.addEventListener("keydown", (e) => { if (e.key === "Enter") jump(); });
    listEl.appendChild(li);
  }

  slider.max = frames.length - 1;
  slider.addEventListener("input", () => setFrame(Number(slider.value)));
  play.addEventListener("click", toggle);
  window.addEventListener("resize", size);
  const mq = window.matchMedia("(prefers-color-scheme: dark)");
  if (mq.addEventListener) mq.addEventListener("change", draw);
  new MutationObserver(draw).observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
  // start on the frame with the most seeding aircraft in view
  setFrame(Math.round(frames.length * 0.6));
  size();
})();
