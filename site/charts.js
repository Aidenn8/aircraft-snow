// Charts and tables for the page, drawn as inline SVG from the page data. Colours come from
// the page's CSS tokens so both themes work; text stays in ink tokens, marks carry identity.
(function () {
  const data = JSON.parse(document.getElementById("page-data").textContent);
  const css = (n) => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
  const NS = "http://www.w3.org/2000/svg";
  const fmt = (v, d = 0) => (v == null ? "–" : Number(v).toFixed(d));
  const minus = (s) => String(s).replace(/^-/, "−");

  function el(name, attrs, parent) {
    const e = document.createElementNS(NS, name);
    for (const k in attrs) e.setAttribute(k, attrs[k]);
    if (parent) parent.appendChild(e);
    return e;
  }

  // one shared tooltip
  const tip = document.createElement("div");
  tip.className = "tip";
  tip.hidden = true;
  document.body.appendChild(tip);
  function showTip(evt, html) {
    tip.innerHTML = html;
    tip.hidden = false;
    const r = tip.getBoundingClientRect();
    let x = evt.clientX + 14, y = evt.clientY + 14;
    if (x + r.width > window.innerWidth - 8) x = evt.clientX - r.width - 14;
    if (y + r.height > window.innerHeight - 8) y = evt.clientY - r.height - 14;
    tip.style.left = x + "px";
    tip.style.top = y + "px";
  }
  const hideTip = () => { tip.hidden = true; };

  function frame(svg, W, H, m, xs, ys, xTicks, yTicks, xLabel, yLabel, xfmt, yfmt) {
    const g = el("g", {}, svg);
    for (const t of yTicks) {
      el("line", { x1: m.l, x2: W - m.r, y1: ys(t), y2: ys(t), stroke: css("--rule"), "stroke-width": 1 }, g);
      const tx = el("text", { x: m.l - 8, y: ys(t) + 4, "text-anchor": "end", class: "tick" }, g);
      tx.textContent = yfmt(t);
    }
    for (const t of xTicks) {
      const tx = el("text", { x: xs(t), y: H - m.b + 18, "text-anchor": "middle", class: "tick" }, g);
      tx.textContent = xfmt(t);
    }
    el("line", { x1: m.l, x2: W - m.r, y1: H - m.b, y2: H - m.b, stroke: css("--muted"), "stroke-width": 1 }, g);
    const xl = el("text", { x: (m.l + W - m.r) / 2, y: H - 6, "text-anchor": "middle", class: "axis-label" }, g);
    xl.textContent = xLabel;
    const yl = el("text", { x: 14, y: (m.t + H - m.b) / 2, "text-anchor": "middle", class: "axis-label",
                            transform: `rotate(-90 14 ${(m.t + H - m.b) / 2})` }, g);
    yl.textContent = yLabel;
  }

  // ---- life cycle: echo on the drifted track minus the sides, with a 90% band ----
  function lifecycle() {
    const box = document.getElementById("chart-lifecycle");
    const L = data.lifecycle;
    if (!box || !L) return;
    box.innerHTML = "";
    const W = 640, H = 300, m = { l: 56, r: 16, t: 14, b: 44 };
    const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img",
      "aria-label": "Reflectivity on the drifted track minus 5 to 8 km either side, against minutes since the aircraft passed" }, box);
    const x0 = -30, x1 = 120, y0 = -4, y1 = 16;
    const xs = (v) => m.l + (v - x0) / (x1 - x0) * (W - m.l - m.r);
    const ys = (v) => m.t + (y1 - v) / (y1 - y0) * (H - m.t - m.b);
    // the matching window, where the composite is shaped by selection
    el("rect", { x: xs(10), y: m.t, width: xs(50) - xs(10), height: H - m.t - m.b, fill: css("--rule"), opacity: 0.5 }, svg);
    const wl = el("text", { x: (xs(10) + xs(50)) / 2, y: m.t + 14, "text-anchor": "middle", class: "note-text" }, svg);
    wl.textContent = "matching window";
    frame(svg, W, H, m, xs, ys, [-30, 0, 30, 60, 90, 120], [-4, 0, 4, 8, 12, 16],
          "minutes since the aircraft passed", "dB above the sides", (t) => minus(t), (t) => minus(t));
    el("line", { x1: xs(0), x2: xs(0), y1: m.t, y2: H - m.b, stroke: css("--muted"), "stroke-width": 1 }, svg);
    const pts = L.age_min.map((a, i) => [a, L.mean_db[i], L.lo_db[i], L.hi_db[i]]).filter((p) => p[1] != null);
    const band = pts.map((p) => `${xs(p[0])},${ys(p[3])}`).concat(pts.slice().reverse().map((p) => `${xs(p[0])},${ys(p[2])}`));
    el("polygon", { points: band.join(" "), fill: css("--seed"), opacity: 0.16 }, svg);
    el("polyline", { points: pts.map((p) => `${xs(p[0])},${ys(p[1])}`).join(" "), fill: "none", stroke: css("--seed"),
                     "stroke-width": 2, "stroke-linejoin": "round" }, svg);
    // hover: nearest age bin
    const cross = el("line", { y1: m.t, y2: H - m.b, stroke: css("--ink-2"), "stroke-width": 1, visibility: "hidden" }, svg);
    const dot = el("circle", { r: 4, fill: css("--seed"), stroke: css("--panel"), "stroke-width": 2, visibility: "hidden" }, svg);
    const hit = el("rect", { x: m.l, y: m.t, width: W - m.l - m.r, height: H - m.t - m.b, fill: "transparent" }, svg);
    hit.addEventListener("mousemove", (e) => {
      const r = svg.getBoundingClientRect();
      const vx = (e.clientX - r.left) / r.width * W;
      const a = x0 + (vx - m.l) / (W - m.l - m.r) * (x1 - x0);
      let best = pts[0];
      for (const p of pts) if (Math.abs(p[0] - a) < Math.abs(best[0] - a)) best = p;
      cross.setAttribute("x1", xs(best[0])); cross.setAttribute("x2", xs(best[0]));
      cross.setAttribute("visibility", "visible");
      dot.setAttribute("cx", xs(best[0])); dot.setAttribute("cy", ys(best[1])); dot.setAttribute("visibility", "visible");
      showTip(e, `<b>${minus(fmt(best[0] - 2.5))} to ${minus(fmt(best[0] + 2.5))} min</b><br>` +
                 `${minus(fmt(best[1], 1))} dB above the sides<br><span class="muted">90%: ${minus(fmt(best[2], 1))} to ${minus(fmt(best[3], 1))}</span>`);
    });
    hit.addEventListener("mouseleave", () => { hideTip(); cross.setAttribute("visibility", "hidden"); dot.setAttribute("visibility", "hidden"); });
  }

  // ---- efficiency: share of first passages that left a streak, by temperature ----
  function efficiency() {
    const box = document.getElementById("chart-efficiency");
    const E = (data.efficiency || []).filter((r) => r.engine_class === "jet" && r.n >= 8);
    if (!box || !E.length) return;
    box.innerHTML = "";
    const W = 640, H = 300, m = { l: 56, r: 16, t: 14, b: 44 };
    const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img",
      "aria-label": "Share of jets crossing the cloud layer that left a streak, by air temperature, with the chance level" }, box);
    const x0 = -8, x1 = -23, y0 = 0, y1 = 0.5;            // colder to the right
    const xs = (v) => m.l + (v - x0) / (x1 - x0) * (W - m.l - m.r);
    const ys = (v) => m.t + (y1 - v) / (y1 - y0) * (H - m.t - m.b);
    frame(svg, W, H, m, xs, ys, [-8, -11, -14, -17, -20, -23], [0, 0.1, 0.2, 0.3, 0.4, 0.5],
          "air temperature where the aircraft crossed the layer (°C)", "share that left a streak",
          (t) => minus(t), (t) => `${Math.round(t * 100)}%`);
    for (const r of E) {
      const cx = xs((r.t_lo + r.t_hi) / 2);
      el("line", { x1: cx, x2: cx, y1: ys(r.ci_lo), y2: ys(r.ci_hi), stroke: css("--seed"), "stroke-width": 2, "stroke-linecap": "round" }, svg);
      el("circle", { cx, cy: ys(r.chance), r: 5, fill: "none", stroke: css("--muted"), "stroke-width": 1.5 }, svg);
      el("circle", { cx, cy: ys(r.share), r: 6, fill: css("--seed"), stroke: css("--panel"), "stroke-width": 2 }, svg);
      const lab = el("text", { x: cx + 10, y: ys(r.share) + 4, class: "tick" }, svg);
      lab.textContent = `${r.made_streak} of ${r.n}`;
      const hit = el("rect", { x: cx - 22, y: m.t, width: 44, height: H - m.t - m.b, fill: "transparent" }, svg);
      hit.addEventListener("mousemove", (e) => showTip(e,
        `<b>${minus(r.t_hi)} to ${minus(r.t_lo)} °C</b><br>${r.made_streak} of ${r.n} jets left a streak ` +
        `(${Math.round(r.share * 100)}%)<br><span class="muted">90% interval ${Math.round(r.ci_lo * 100)}–${Math.round(r.ci_hi * 100)}%; ` +
        `by chance ${Math.round(r.chance * 100)}%</span>`));
      hit.addEventListener("mouseleave", hideTip);
    }
  }

  // ---- the matched aircraft ----
  function aircraftTable() {
    const body = document.getElementById("aircraft-rows");
    if (!body) return;
    body.innerHTML = "";
    for (const a of data.case.aircraft) {
      const tr = document.createElement("tr");
      const cells = [
        `<span class="mono">${a.callsign}</span>`, `<span class="mono">${a.icao_type}</span>`, a.phase,
        `<span class="num">${a.pass_utc}</span>`, `<span class="num">${Math.round(a.seed_alt_ft / 100) * 100}</span>`,
        `<span class="num">${minus(fmt(a.seed_temp_c, 0))}</span>`, `<span class="num">${fmt(a.age_min, 0)}</span>`,
        `<span class="num">${fmt(a.v_fall, 1)}</span>`, `<span class="num">${fmt(a.length_km, 0)}</span>`,
        `<span class="num">${fmt(a.snow_mm_h_max, 1)}</span>`];
      tr.innerHTML = cells.map((c) => `<td>${c}</td>`).join("");
      body.appendChild(tr);
    }
  }

  function fillNumbers() {
    const c = data.case;
    const set = (id, v) => { const e = document.getElementById(id); if (e) e.textContent = v; };
    set("n-flights", c.flights_p01);
    set("n-expected", fmt(c.expected_false_p01, 1));
    set("n-objects", c.n_streak_objects);
    set("e-real", fmt(c.E, 0));
    set("e-null", fmt(c.E_null_median, 0));
    set("e-q95", fmt(c.E_null_q95, 0));
    set("n-rot", c.n_rotations);
  }

  // ---- census: event days by radar and winter ----
  function censusRates() {
    const box = document.getElementById("chart-census");
    const C = data.census;
    if (!box || !C || !C.rates) return;
    box.innerHTML = "";
    const rates = C.rates;
    const winters = [...new Set(rates.map((r) => r.winter))].sort();
    const sites = [...new Set(rates.map((r) => r.site))];
    const total = (s) => rates.filter((r) => r.site === s).reduce((a, r) => a + r.event_days, 0);
    sites.sort((a, b) => total(b) - total(a));
    const colours = [css("--seed"), css("--accent-ink"), css("--muted")];
    const W = 640, H = 300, m = { l: 56, r: 16, t: 14, b: 44 };
    const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Event days per radar and winter" }, box);
    const ymax = Math.max(1, ...rates.map((r) => r.event_days));
    const top = Math.ceil(ymax * 1.15);
    const band = (W - m.l - m.r) / sites.length;
    const bw = Math.min(18, (band * 0.8) / winters.length);
    const ys = (v) => m.t + (top - v) / top * (H - m.t - m.b);
    const step = top <= 6 ? 1 : Math.ceil(top / 5);
    const ticks = []; for (let t = 0; t <= top; t += step) ticks.push(t);
    frame(svg, W, H, m, (v) => v, ys, [], ticks, "", "event days", (t) => t, (t) => t);
    sites.forEach((s, i) => {
      const cx = m.l + band * (i + 0.5);
      const lab = el("text", { x: cx, y: H - m.b + 18, "text-anchor": "middle", class: "tick" }, svg);
      lab.textContent = s;
      winters.forEach((w, j) => {
        const r = rates.find((q) => q.site === s && q.winter === w);
        if (!r) return;
        const x = cx + (j - (winters.length - 1) / 2) * (bw + 2) - bw / 2;
        const y = ys(r.event_days);
        const h = Math.max(ys(0) - y, r.event_days > 0 ? 2 : 0);
        el("rect", { x, y: ys(0) - h, width: bw, height: h, rx: 2, fill: colours[j % colours.length] }, svg);
        const hit = el("rect", { x: x - 2, y: m.t, width: bw + 4, height: H - m.t - m.b, fill: "transparent" }, svg);
        hit.addEventListener("mousemove", (e) => showTip(e, `<b>${s}, winter ${w}</b><br>${r.event_days} event days ` +
          `of ${r.tested} tested (${r.radar_days} radar-days)`));
        hit.addEventListener("mouseleave", hideTip);
      });
    });
    const leg = document.getElementById("legend-census");
    if (leg) leg.innerHTML = winters.map((w, j) =>
      `<span><i class="sw" style="border-radius:2px;background:${colours[j % colours.length]}"></i>winter ${w}</span>`).join("");
  }

  // ---- census: seeding temperature by engine class ----
  function censusTemps() {
    const box = document.getElementById("chart-temps");
    const S = (data.census && data.census.seeders) || [];
    if (!box || !S.length) return;
    box.innerHTML = "";
    const classes = ["jet", "turboprop", "piston"].filter((c) => S.some((r) => r.engine_class === c));
    const W = 640, H = 260, m = { l: 90, r: 16, t: 14, b: 44 };
    const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Seeding temperature of each flight by engine class" }, box);
    const x0 = 0, x1 = -30;
    const xs = (v) => m.l + (v - x0) / (x1 - x0) * (W - m.l - m.r);
    const row = (H - m.t - m.b) / classes.length;
    frame(svg, W, H, m, xs, (v) => v, [0, -5, -10, -15, -20, -25, -30], [], "air temperature where the aircraft seeded (°C)", "",
          (t) => minus(t), (t) => t);
    classes.forEach((c, i) => {
      const cy = m.t + row * (i + 0.5);
      const v = S.filter((r) => r.engine_class === c && r.seed_temp_c != null);
      const lab = el("text", { x: m.l - 10, y: cy + 4, "text-anchor": "end", class: "axis-label" }, svg);
      lab.textContent = `${c}s (${v.length})`;
      v.forEach((r, k) => {
        const jitter = ((k * 37) % 17) / 17 - 0.5;
        const dot = el("circle", { cx: xs(r.seed_temp_c), cy: cy + jitter * row * 0.5, r: 4,
                                   fill: c === "jet" ? css("--seed") : css("--accent-ink"), opacity: 0.8 }, svg);
        dot.addEventListener("mousemove", (e) => showTip(e, `<b>${r.callsign} (${r.icao_type})</b><br>${r.site} ${r.day}<br>` +
          `${minus(fmt(r.seed_temp_c, 1))} °C at ${Math.round(r.seed_alt_ft / 100) * 100} ft`));
        dot.addEventListener("mouseleave", hideTip);
      });
      if (v.length) {
        const med = v.map((r) => r.seed_temp_c).sort((a, b) => a - b)[Math.floor(v.length / 2)];
        el("line", { x1: xs(med), x2: xs(med), y1: cy - row * 0.35, y2: cy + row * 0.35, stroke: css("--ink"), "stroke-width": 2 }, svg);
      }
    });
  }

  function censusNumbers() {
    const C = data.census;
    const sec = document.getElementById("census");
    if (!C || !sec) return;
    sec.hidden = false;
    const set = (id, v) => { const e = document.getElementById(id); if (e) e.textContent = v; };
    const r = C.rates || [];
    set("c-radar-days", r.reduce((a, q) => a + q.radar_days, 0).toLocaleString());
    set("c-tested", r.reduce((a, q) => a + q.tested, 0).toLocaleString());
    set("c-events", r.reduce((a, q) => a + q.event_days, 0));
    set("c-seeders", (C.seeders || []).length);
    if (C.null) {
      set("c-null-days", C.null.tested.toLocaleString());
      set("c-null-selected", C.null.selected);
    }
  }

  function drawAll() { lifecycle(); efficiency(); censusRates(); censusTemps(); }
  censusNumbers();
  aircraftTable();
  fillNumbers();
  drawAll();
  // redraw when the theme changes so the token colours follow
  window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", drawAll);
  new MutationObserver(drawAll).observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
})();
