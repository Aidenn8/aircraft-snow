// Charts and tables for the page, drawn as inline SVG from the page data. Colours come from
// the page's CSS tokens so both themes work; text stays in ink tokens, marks carry identity.
(function () {
  const data = JSON.parse(document.getElementById("page-data").textContent);
  const css = (n) => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
  const NS = "http://www.w3.org/2000/svg";
  const fmt = (v, d = 0) => (v == null ? "–" : Number(v).toFixed(d));
  const minus = (s) => String(s).replace(/^-/, "−");

  // charts are drawn in a 640-unit-wide coordinate system, or 420 when the box is narrow
  // (phones), so that their text is not scaled down to unreadable sizes
  const chartW = (box) => (box && box.clientWidth && box.clientWidth < 520 ? 420 : 640);

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
    const W = chartW(box), H = 300, m = { l: 56, r: 16, t: 14, b: 44 };
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
    const W = chartW(box), H = 300, m = { l: 56, r: 16, t: 14, b: 44 };
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

  // ---- census: days beyond chance per radar, against the no-seeding expectation ----
  const SITE_ORDER = ["KFTG", "KMPX", "KLOT", "KIND", "KDTX", "KCLE", "KBUF", "KPBZ"];
  function censusSites() {
    const box = document.getElementById("chart-census-sites");
    const C = data.census;
    if (!box || !C || !C.sites) return;
    box.innerHTML = "";
    const rows = SITE_ORDER.map((s) => C.sites.find((r) => r.site === s)).filter(Boolean);
    const W = chartW(box), H = 300, m = { l: 150, r: 24, t: 14, b: 44 };
    const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img",
      "aria-label": "Days with p at or below 0.01 at each radar, against the number expected with no seeding" }, box);
    const top = Math.max(2, ...rows.map((r) => Math.max(r.p01, r.expected_p01))) * 1.12;
    const xs = (v) => m.l + v / top * (W - m.l - m.r);
    const step = top <= 8 ? 1 : top <= 20 ? 2 : 5;
    const ticks = []; for (let t = 0; t <= top; t += step) ticks.push(t);
    for (const t of ticks) el("line", { x1: xs(t), x2: xs(t), y1: m.t, y2: H - m.b, stroke: css("--rule"), "stroke-width": 1 }, svg);
    frame(svg, W, H, m, xs, (v) => v, ticks, [], "radar-days with p ≤ 0.01", "", (t) => t, (t) => t);
    const band = (H - m.t - m.b) / rows.length;
    rows.forEach((r, i) => {
      const cy = m.t + band * (i + 0.5);
      const lab = el("text", { x: m.l - 10, y: cy + 4, "text-anchor": "end", class: "axis-label" }, svg);
      lab.textContent = `${r.name} (${r.site})`;
      el("rect", { x: xs(0), y: cy - band * 0.28, width: xs(r.expected_p01) - xs(0), height: band * 0.56, rx: 2,
                   fill: css("--rule-strong") }, svg);
      el("circle", { cx: xs(r.p01), cy, r: 6, fill: css("--accent-ink"), stroke: css("--panel"), "stroke-width": 2 }, svg);
      const hit = el("rect", { x: m.l, y: cy - band / 2, width: W - m.l - m.r, height: band, fill: "transparent" }, svg);
      hit.addEventListener("mousemove", (e) => showTip(e, `<b>${r.name} (${r.site})</b><br>${r.p01} days with p ≤ 0.01 ` +
        `of ${r.tested} tested<br><span class="muted">${fmt(r.expected_p01, 1)} expected with no seeding; ` +
        `${r.p05} at p ≤ 0.05 (${fmt(r.expected_p05, 1)} expected)</span>`));
      hit.addEventListener("mouseleave", hideTip);
    });
    const leg = document.getElementById("legend-census-sites");
    if (leg) leg.innerHTML = `<span><i class="sw" style="border-radius:50%;background:${css("--accent-ink")}"></i>observed</span>` +
      `<span><i class="sw" style="border-radius:2px;background:${css("--rule-strong")}"></i>expected with no seeding</span>`;
  }

  // ---- census: days beyond chance by hours of supercooled cloud ----
  function wilson(k, n, z = 1.645) {
    const p = k / n, d = 1 + z * z / n;
    const c = (p + z * z / (2 * n)) / d, h = z * Math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d;
    return [c - h, c + h];
  }
  function censusSlw() {
    const box = document.getElementById("chart-census-slw");
    const S = data.census && data.census.slw;
    if (!box || !S || !S.length) return;
    box.innerHTML = "";
    const W = chartW(box), H = 300, m = { l: 56, r: 78, t: 14, b: 58 };
    const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img",
      "aria-label": "Share of days with p at or below 0.05, by hours of supercooled cloud over the radar" }, box);
    const rows = S.map((r) => ({ ...r, share: r.p05 / r.tested, ci: wilson(r.p05, r.tested) }));
    const y1 = Math.max(0.1, ...rows.map((r) => r.ci[1])) * 1.1;
    const ys = (v) => m.t + (y1 - v) / y1 * (H - m.t - m.b);
    const band = (W - m.l - m.r) / rows.length;
    const ticks = []; for (let t = 0; t <= y1 + 1e-9; t += 0.05) ticks.push(t);
    frame(svg, W, H, m, (v) => v, ys, [], ticks, "hours of supercooled cloud water over the radar (HRRR)",
          "days with p ≤ 0.05", (t) => t, (t) => `${Math.round(t * 100)}%`);
    el("line", { x1: m.l, x2: W - m.r, y1: ys(0.05), y2: ys(0.05), stroke: css("--ink-2"), "stroke-width": 1.2, "stroke-dasharray": "5 4" }, svg);
    const nl = el("text", { x: W - m.r + 6, y: ys(0.05) - 3, "text-anchor": "start", class: "note-text" }, svg);
    nl.textContent = "no seeding";
    const nl2 = el("text", { x: W - m.r + 6, y: ys(0.05) + 10, "text-anchor": "start", class: "note-text" }, svg);
    nl2.textContent = "gives 5%";
    rows.forEach((r, i) => {
      const cx = m.l + band * (i + 0.5), bw = Math.min(70, band * 0.55);
      el("rect", { x: cx - bw / 2, y: ys(r.share), width: bw, height: ys(0) - ys(r.share), rx: 3, fill: css("--accent-ink") }, svg);
      el("line", { x1: cx, x2: cx, y1: ys(r.ci[0]), y2: ys(r.ci[1]), stroke: css("--ink"), "stroke-width": 1.5 }, svg);
      const l1 = el("text", { x: cx, y: H - m.b + 18, "text-anchor": "middle", class: "tick" }, svg);
      l1.textContent = r.slw_class;
      const l2 = el("text", { x: cx, y: H - m.b + 32, "text-anchor": "middle", class: "note-text" }, svg);
      l2.textContent = `${r.tested} days`;
      const hit = el("rect", { x: cx - band / 2, y: m.t, width: band, height: H - m.t - m.b, fill: "transparent" }, svg);
      hit.addEventListener("mousemove", (e) => showTip(e, `<b>${r.slw_class} of supercooled cloud</b><br>${r.p05} of ` +
        `${r.tested} days at p ≤ 0.05 (${fmt(r.share * 100, 1)}%)<br><span class="muted">90% interval ` +
        `${fmt(r.ci[0] * 100, 1)}–${fmt(r.ci[1] * 100, 1)}%; ${r.p01} at p ≤ 0.01</span>`));
      hit.addEventListener("mouseleave", hideTip);
    });
    // the axis labels sit lower than usual to leave room for the day counts
    svg.querySelectorAll("text.axis-label")[0].setAttribute("y", H - 4);
  }

  // ---- census: share of jets leaving a streak by temperature, arrivals and departures ----
  function censusPhase() {
    const box = document.getElementById("chart-census-phase");
    const E = ((data.census && data.census.efficiency) || []).filter((r) => r.engine_class === "jet" && r.tbin !== "all" &&
      (r.phase === "arrival" || r.phase === "departure") && r.n >= 20);
    if (!box || !E.length) return;
    box.innerHTML = "";
    const mid = (s) => { const [a, b] = s.replace(/[(\]]/g, "").split(",").map(Number); return (a + b) / 2; };
    const W = chartW(box), H = 300, m = { l: 56, r: 16, t: 14, b: 44 };
    const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img",
      "aria-label": "Share of jets crossing the layer that left a streak, by temperature, for arrivals and departures" }, box);
    // the axis runs to just past the coldest bin's midpoint (the coldest bin is -30 to -20 C)
    const x0 = -4, x1 = Math.min(-20, ...E.map((r) => mid(r.tbin))) - 1.5;
    const y1 = Math.max(0.2, ...E.map((r) => r.hi90)) * 1.08;
    const xs = (v) => m.l + (v - x0) / (x1 - x0) * (W - m.l - m.r);
    const ys = (v) => m.t + (y1 - v) / y1 * (H - m.t - m.b);
    const ticks = []; for (let t = 0; t <= y1 + 1e-9; t += 0.05) ticks.push(t);
    const xticks = []; for (let t = -4; t >= x1; t -= 4) xticks.push(t);
    frame(svg, W, H, m, xs, ys, xticks, ticks, "air temperature where the aircraft crossed the layer (°C)",
          "share that left a streak", (t) => minus(t), (t) => `${Math.round(t * 100)}%`);
    const series = [["arrival", css("--accent-ink"), -0.18], ["departure", css("--seed"), 0.18]];
    for (const [ph, col, dx] of series) {
      const g = E.filter((r) => r.phase === ph).map((r) => ({ ...r, mid: mid(r.tbin) + dx })).sort((a, b) => b.mid - a.mid);
      el("polyline", { points: g.map((r) => `${xs(r.mid)},${ys(r.made_streak)}`).join(" "), fill: "none", stroke: col,
                       "stroke-width": 2, "stroke-linejoin": "round" }, svg);
      for (const r of g) {
        const cx = xs(r.mid);
        el("line", { x1: cx, x2: cx, y1: ys(r.lo90), y2: ys(r.hi90), stroke: col, "stroke-width": 1.5 }, svg);
        el("circle", { cx, cy: ys(r.chance), r: 4, fill: css("--panel"), stroke: col, "stroke-width": 1.5 }, svg);
        el("circle", { cx, cy: ys(r.made_streak), r: 5.5, fill: col, stroke: css("--panel"), "stroke-width": 2 }, svg);
        const hit = el("rect", { x: cx - 14, y: m.t, width: 28, height: H - m.t - m.b, fill: "transparent" }, svg);
        const [a, b] = r.tbin.replace(/[(\]]/g, "").split(",");
        hit.addEventListener("mousemove", (e) => showTip(e, `<b>jet ${ph}s, ${minus(b.trim())} to ${minus(a.trim())} °C</b><br>` +
          `${r.k} of ${r.n} left a streak (${fmt(r.made_streak * 100, 1)}%)<br><span class="muted">90% interval ` +
          `${fmt(r.lo90 * 100, 1)}–${fmt(r.hi90 * 100, 1)}%; by chance ${fmt(r.chance * 100, 1)}%</span>`));
        hit.addEventListener("mouseleave", hideTip);
      }
    }
    const leg = document.getElementById("legend-census-phase");
    if (leg) leg.innerHTML = series.map(([ph, col]) => `<span><i class="sw" style="border-radius:50%;background:${col}"></i>jet ${ph}s</span>`).join("") +
      `<span><i class="sw" style="border-radius:50%;background:transparent;border:1.5px solid ${css("--ink-2")}"></i>chance (rotated tracks)</span>`;
  }

  // ---- census: seeding temperature by engine class, against the quoted thresholds ----
  function censusTemps() {
    const box = document.getElementById("chart-temps");
    const S = (data.census && data.census.seeders) || [];
    if (!box || !S.length) return;
    box.innerHTML = "";
    const classes = ["jet", "turboprop", "piston"].filter((c) => S.some((r) => r.engine_class === c));
    const W = chartW(box), H = 300, m = { l: 110, r: 16, t: 30, b: 44 };
    const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Seeding temperature of each catalogue flight by engine class" }, box);
    const x0 = -4, x1 = -28;
    const xs = (v) => m.l + (v - x0) / (x1 - x0) * (W - m.l - m.r);
    // thresholds usually quoted (NCAR 2011): jets about -20 to -25 C, turboprops about -10 to -15 C
    for (const [a, b, lab] of [[-10, -15, "quoted for turboprops"], [-20, -25, "quoted for jets"]]) {
      el("rect", { x: xs(a), y: m.t, width: xs(b) - xs(a), height: H - m.t - m.b, fill: css("--rule"), opacity: 0.7 }, svg);
      const t = el("text", { x: (xs(a) + xs(b)) / 2, y: m.t - 10, "text-anchor": "middle", class: "note-text" }, svg);
      t.textContent = lab;
    }
    frame(svg, W, H, m, xs, (v) => v, [-4, -8, -12, -16, -20, -24, -28], [], "air temperature where the aircraft seeded (°C)", "",
          (t) => minus(t), (t) => t);
    const row = (H - m.t - m.b) / classes.length;
    classes.forEach((c, i) => {
      const cy = m.t + row * (i + 0.5);
      const v = S.filter((r) => r.engine_class === c && r.seed_temp_c != null);
      const lab = el("text", { x: m.l - 10, y: cy + 4, "text-anchor": "end", class: "axis-label" }, svg);
      lab.textContent = `${c}s (${v.length})`;
      v.forEach((r, k) => {
        const jitter = ((k * 37) % 17) / 17 - 0.5;
        const dot = el("circle", { cx: xs(r.seed_temp_c), cy: cy + jitter * row * 0.5, r: 4,
                                   fill: c === "jet" ? css("--accent-ink") : css("--seed"), opacity: 0.8 }, svg);
        dot.addEventListener("mousemove", (e) => showTip(e, `<b>${r.callsign} (${r.type_name})</b><br>${r.site_name}, ${r.date}<br>` +
          `${minus(fmt(r.seed_temp_c, 1))} °C at ${(Math.round(r.seed_alt_ft / 100) * 100).toLocaleString()} ft, ${r.phase}`));
        dot.addEventListener("mouseleave", hideTip);
      });
      if (v.length) {
        const sorted = v.map((r) => r.seed_temp_c).sort((a, b) => a - b);
        const med = sorted.length % 2 ? sorted[(sorted.length - 1) / 2] : (sorted[sorted.length / 2 - 1] + sorted[sorted.length / 2]) / 2;
        el("line", { x1: xs(med), x2: xs(med), y1: cy - row * 0.35, y2: cy + row * 0.35, stroke: css("--ink"), "stroke-width": 2.5 }, svg);
      }
    });
  }

  // ---- census: the catalogue ----
  function catalogueTable() {
    const body = document.getElementById("catalogue-rows");
    const S = (data.census && data.census.seeders) || [];
    if (!body || !S.length) return;
    body.innerHTML = S.map((r) => `<tr><td>${r.date}</td><td>${r.site_name}</td><td>${r.callsign}</td>` +
      `<td title="${r.type_name}">${r.icao_type}</td><td>${r.phase}</td><td class="num">${(Math.round(r.seed_alt_ft / 100) * 100).toLocaleString()}</td>` +
      `<td class="num">${minus(fmt(r.seed_temp_c, 0))}</td><td class="num">${fmt(r.age_min, 0)}</td>` +
      `<td class="num">${fmt(r.length_km, 0)}</td><td class="num">${fmt(r.snow_mm_h_max, 1)}</td></tr>`).join("");
  }

  function censusNumbers() {
    const C = data.census;
    const sec = document.getElementById("census");
    if (!C || !sec) return;
    sec.hidden = false;
    const set = (id, v) => { const e = document.getElementById(id); if (e) e.textContent = v; };
    const T = C.total || {};
    set("c-winters", C.winters_text || "");
    set("c-radar-days", (T.radar_days || 0).toLocaleString());
    set("c-tested", (T.tested || 0).toLocaleString());
    set("c-p01", T.p01);
    set("c-exp01", fmt(T.expected_p01, 1));
    const K = C.catalogue || {};
    set("c-cat-flights", K.flights);
    set("c-cat-days", K.days);
    set("c-cat-radars", K.radars);
    set("c-cat-false", fmt(K.expected_false_p01, 0));
    const note = (id, html) => { const e = document.getElementById(id); if (e && html) e.innerHTML = html; };
    note("c-null-note", C.null_note);
    note("c-band-note", C.band_note);
    note("c-phase-note", C.phase_note);
  }

  function drawAll() { lifecycle(); efficiency(); censusSites(); censusSlw(); censusPhase(); censusTemps(); }
  censusNumbers();
  catalogueTable();
  aircraftTable();
  fillNumbers();
  drawAll();
  // redraw on resize (the coordinate width depends on the box) and when the theme changes
  let resizeTimer = null;
  window.addEventListener("resize", () => { clearTimeout(resizeTimer); resizeTimer = setTimeout(drawAll, 150); });
  window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", drawAll);
  new MutationObserver(drawAll).observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
})();
