/* live.js -- the app's live design view: the same sums as the Python programs, fast enough to run on
   every keystroke. Build (the Python programs) stays the authority for the parts; tests/live_check.py
   compares the two.
     horn:    hornLive(computeFamily, opts)     section, mould regions, numbers, 3-D outline
     box:     boxLive(driver, opts)             alignments, curves, box size, port, warnings
     bracket: bracketLive(horn, box, opts)      flat pattern, numbers, inserts                    */
(function (root) {
  "use strict";
  const C = 343, RHO = 1.204, P0 = 20e-6, PI = Math.PI;
  const num = (v, dflt) => { if (v === undefined || v === null || v === "" || v === false) return dflt;
    const x = parseFloat(String(v).replace(",", ".")); return isNaN(x) ? dflt : x; };

  // ---------- complex numbers ----------
  const cx = (re, im) => ({ re, im: im || 0 });
  const add = (a, b) => cx(a.re + b.re, a.im + b.im);
  const sub = (a, b) => cx(a.re - b.re, a.im - b.im);
  const mul = (a, b) => cx(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
  const div = (a, b) => { const d = b.re * b.re + b.im * b.im; return cx((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d); };
  const sc = (a, k) => cx(a.re * k, a.im * k);
  const inv = a => div(cx(1), a);
  const abs = a => Math.hypot(a.re, a.im);

  // ======================================================================= box
  function driverParams(d) {
    const p = { Re: d.re_ohm, Le: (d.le_mh || 0) / 1000, Bl: d.bl_tm, Mms: d.mms_g / 1000, Cms: d.cms_mm_per_n / 1000,
      Sd: d.sd_cm2 / 1e4, fs: d.fs_hz, Qms: d.qms, Qts: d.qts, Vas: d.vas_l / 1000, Xmax: d.xmax_mm / 1000 };
    p.Rms = 2 * PI * p.fs * p.Mms / p.Qms;
    return p;
  }
  function logspace(a, b, n) { const out = []; for (let i = 0; i < n; i++) out.push(Math.pow(10, Math.log10(a) + (Math.log10(b) - Math.log10(a)) * i / (n - 1))); return out; }

  // the lumped circuit of boxmould.py response(): Dp is the port's (equivalent) diameter in mm
  function response(p, f, Vb, kind, fb, Dp, QL, volts, hpf, order) {
    volts = volts === undefined ? 2.83 : volts;
    const Cab = Vb / (RHO * C * C);
    const wref = 2 * PI * (kind === "vented" ? fb : p.fs * Math.sqrt(1 + p.Vas / Vb));
    const Ral = QL / (wref * Cab);
    const Sp = kind === "vented" ? PI * Math.pow(Dp / 2000, 2) : 0, Map = kind === "vented" ? 1 / (wref * wref * Cab) : 0;
    const out = { spl: [], x: [], Zin: [], vport: kind === "vented" ? [] : null, pc: [] };
    for (let i = 0; i < f.length; i++) {
      const w = 2 * PI * f[i], s = cx(0, w);
      const Ze = add(cx(p.Re), sc(s, p.Le));
      const Zm = add(add(sc(s, p.Mms), cx(p.Rms)), inv(sc(s, p.Cms)));
      let Y = add(sc(s, Cab), cx(1 / Ral)), Zp = null;
      if (kind === "vented") { Zp = add(sc(s, Map), cx(1e-3 * Map * wref)); Y = add(Y, inv(Zp)); }
      const Zn = inv(Y);
      let e = cx(volts);
      if (hpf) {
        const x = sc(s, 1 / (2 * PI * hpf)), x2 = mul(x, x);
        let H;
        if (order === 2) H = div(x2, add(add(x2, sc(x, Math.SQRT2)), cx(1)));
        else if (order === 4) H = div(mul(x2, x2), mul(add(add(x2, sc(x, 0.7654)), cx(1)), add(add(x2, sc(x, 1.8478)), cx(1))));
        else H = div(x, add(x, cx(1)));
        e = mul(e, H);
      }
      const ZmT = add(Zm, sc(Zn, p.Sd * p.Sd));
      const u = div(sc(e, p.Bl), add(mul(Ze, ZmT), cx(p.Bl * p.Bl)));
      out.x.push(abs(div(u, s)) * Math.SQRT2);
      out.Zin.push(abs(add(Ze, div(cx(p.Bl * p.Bl), ZmT))));
      const Ud = sc(u, p.Sd), pb = sc(mul(Ud, Zn), -1);
      let U = add(Ud, sc(pb, 1 / Ral));
      if (Zp) { const Up = div(pb, Zp); U = add(U, Up); out.vport.push(abs(Up) / Sp * Math.SQRT2); }
      out.spl.push(20 * Math.log10(RHO * w * abs(U) / (2 * PI) / P0));
      out.pc.push(mul(cx(0, RHO * w / (2 * PI)), U));           // the pressure itself (Pa, half space, 1 m), to add drivers
    }
    return out;
  }
  function interp(x, xs, ys) {
    if (x <= xs[0]) return ys[0];
    for (let i = 1; i < xs.length; i++) if (xs[i] >= x) { const t = (x - xs[i - 1]) / (xs[i] - xs[i - 1]); return ys[i - 1] + t * (ys[i] - ys[i - 1]); }
    return ys[ys.length - 1];
  }
  function f3Ripple(f, spl, fref) {
    fref = fref || 300;
    const ref = interp(fref, f, spl);
    let hump = -Infinity, last = -1;
    for (let i = 0; i < f.length; i++) {
      const d = spl[i] - ref;
      if (f[i] < fref) { if (d > hump) hump = d; if (d < -3) last = i; }
    }
    let f3 = f[0];
    if (last >= 0) { const d0 = spl[last] - ref, d1 = spl[last + 1] - ref; f3 = f[last] + (-3 - d0) / (d1 - d0) * (f[last + 1] - f[last]); }
    return { f3, hump, ref };
  }
  function bestFb(p, f, Vb, Deq, QL) {
    let best = null;
    for (let fb = 0.6 * p.fs; fb < 1.8 * p.fs; fb += 0.25) {
      const r = response(p, f, Vb, "vented", fb, Deq, QL);
      const k = f3Ripple(f, r.spl);
      if (k.hump <= 0.25 && (!best || k.f3 < best.f3)) best = { fb, f3: k.f3 };
    }
    return best;
  }
  function portLength(Vb, fb, DeqM) {
    const Leff = PI * Math.pow(DeqM / 2, 2) * C * C / (4 * PI * PI * fb * fb * Vb);
    return { Leff, Lphys: Leff - 1.463 * DeqM / 2 };
  }
  function chamferD(d, o) {
    return o.chamfer === false ? d.cutout_d_mm : Math.min(d.cutout_d_mm + 2 * (o.wall - d.frame_t_mm - o.straight), o.chamfer_max);
  }
  function innerWidth(d, o) { return o.width || Math.max(d.outer_d_mm + 2 * o.margin, chamferD(d, o) + 10); }
  function parseInner(s) {
    if (!s) return null;
    const m = String(s).match(/^\s*(\d+(?:[.,]\d+)?)\s*[x×X*,;/ ]\s*(\d+(?:[.,]\d+)?)\s*[x×X*,;/ ]\s*(\d+(?:[.,]\d+)?)/);
    return m ? [1, 2, 3].map(i => parseFloat(m[i].replace(",", "."))) : null;
  }
  // boxmould.py size_box(): keep the two together
  function sizeBox(d, o, VbL, kind, Lphys) {
    const t = o.wall, inner = parseInner(o.inner);
    let Wi = inner ? inner[0] : innerWidth(d, o), Hi, Di;
    const port_in = kind === "vented" ? Lphys - t : 0;
    let disp = 0;
    if (kind === "vented" && o.port_shape === "slot") disp = Wi * (o.slot_h_eff + o.shelf) * port_in / 1e6;
    else if (kind === "vented") disp = PI * Math.pow((o.port_d + 2 * o.tube_wall) / 2, 2) * port_in / 1e6;
    const gross = VbL + (d.displacement_l || 0.8) + disp;
    if (inner) [Wi, Hi, Di] = inner;
    else {   // as boxmould.py size_box()
      const hole = Math.ceil((chamferD(d, o) + 4) / 5) * 5;
      if (o.height) Hi = o.height;
      else if (o.depth) Hi = Math.max(Math.ceil(gross * 1e6 / (Wi * o.depth) / 5) * 5, hole);
      else if (o.hw_ratio) Hi = Math.max(Math.ceil(Wi * o.hw_ratio / 5) * 5, hole);
      else Hi = Math.max(Math.ceil(Math.sqrt(gross * 1e6 / (Wi * o.ratio)) / 5) * 5, hole);
      Di = o.depth || Math.ceil(gross * 1e6 / (Wi * Hi) / 5) * 5;
    }
    return { Wi, Hi, Di, disp, port_in };
  }

  function boxOpts(raw, d) {
    const o = {
      type: raw.type || "vented", vb: num(raw.vb, null), fb: num(raw.fb, null), f3: num(raw.f3, null), qtc: num(raw.qtc, 0.707),
      ql: num(raw.ql, null), port_shape: raw.port_shape || "slot", slot_h: num(raw.slot_h, null), shelf: num(raw.shelf, null),
      port_d: num(raw.port_d, 75), flare: num(raw.flare, 15), tube_wall: num(raw.tube_wall, 3), port_on: raw.port_on || "auto",
      hpf: num(raw.hpf, 32), hpf_order: num(raw.hpf_order, 2), amp_v: num(raw.amp_v, 15), wall: num(raw.wall, 35),
      inner: raw.inner && !/^auto$/i.test(raw.inner) ? raw.inner : null, width: num(raw.width, null), ratio: num(raw.ratio, 1.17),
      hw_ratio: num(raw.hw_ratio, null), height: num(raw.height, null), depth: num(raw.depth, null), driver_y: num(raw.driver_y, null), compare: raw.compare === true || raw.compare === "true",
      margin: 12, straight: 10, chamfer: true, chamfer_max: 230 };
    if (o.shelf === null) o.shelf = o.wall;          // the slot's shelf: as thick as the wall
    o.power = d.power_rms_w || 50; o.znom = d.z_nominal_ohm || 8;
    if (o.height && o.depth && !o.inner) o.inner = `${innerWidth(d, o)}x${o.height}x${o.depth}`;   // fully set: the volume follows
    return o;
  }

  function boxLive(d, raw, nf) {
    const o = boxOpts(raw, d), p = driverParams(d), notes = [];
    const f = logspace(10, 2000, nf || 1200), fs = logspace(10, 2000, 400);
    // the port: its area; a slot as wide as the box inside
    if (o.port_shape === "slot") {
      const inner = parseInner(o.inner), Wi_ = inner ? inner[0] : innerWidth(d, o);
      o.slot_h_eff = o.slot_h || Math.round(PI * Math.pow(o.port_d / 2, 2) / Wi_ * 10) / 10;
      o.port_area = Wi_ * o.slot_h_eff; o._Wi = Wi_;
    } else o.port_area = PI * Math.pow(o.port_d / 2, 2);
    const Deq = 2 * Math.sqrt(o.port_area / PI);
    const qlOf = k => o.ql || (k === "vented" ? 7 : 50);
    function design(kind) {
      if (kind === "sealed") return { Vb: (o.vb && kind === o.type) ? o.vb / 1000 : p.Vas / (Math.pow(o.qtc / p.Qts, 2) - 1), fb: null };
      if (o.vb && o.fb) return { Vb: o.vb / 1000, fb: o.fb };
      if (o.vb) { const b = bestFb(p, f, o.vb / 1000, Deq, o.ql || 7); return { Vb: o.vb / 1000, fb: b ? b.fb : 0.9 * p.fs }; }
      if (o.f3) {   // smallest box (0.5 L steps) reaching f3: bisection, then the exact step
        // as boxmould.py: the first volume from 5 L up (0.5 L steps) whose best tuning reaches f3; coarse, then fine
        const ok = (v, grid) => { const b = bestFb(p, grid, v / 1000, Deq, o.ql || 7); return b && b.f3 <= o.f3; };
        let hit = null;
        for (let v = 5; v <= 120; v += 2) if (ok(v, fs)) { hit = v; break; }
        if (hit === null) throw new Error("that f3 needs more than 120 L");
        for (let v = Math.max(5, hit - 2); v <= hit + 1; v += 0.5) if (ok(v, f)) { const b = bestFb(p, f, v / 1000, Deq, o.ql || 7); return { Vb: v / 1000, fb: b.fb }; }
        const b = bestFb(p, f, hit / 1000, Deq, o.ql || 7); return { Vb: hit / 1000, fb: b.fb };
      }
      throw new Error("vented: give a net volume (and optionally fb) or a target f3");
    }
    // an explicit inside size sets the volume
    const inner = parseInner(o.inner);
    if (inner) {
      const gross = inner[0] * inner[1] * inner[2] / 1e6, asked = o.vb;
      let vb = gross - (d.displacement_l || 0.8);
      for (let i = 0; i < 6; i++) {
        o.vb = vb;
        if (o.type !== "vented") break;
        const ds = design("vented");
        const Lp = portLength(ds.Vb, ds.fb, Deq / 1000).Lphys * 1000 - o.wall;
        const disp = (o.port_shape === "slot" ? inner[0] * (o.slot_h_eff + o.shelf) * Lp : PI * Math.pow((o.port_d + 2 * o.tube_wall) / 2, 2) * Lp) / 1e6;
        vb = gross - (d.displacement_l || 0.8) - Math.max(disp, 0);
      }
      o.vb = Math.round(vb * 100) / 100;
      if (asked === null || Math.abs(asked - o.vb) > 0.05)
        notes.push(`The inside size ${o.inner} gives ${o.vb.toFixed(2)} L net, used instead of the net volume.`);
    }
    const kinds = [o.type].concat(o.compare ? ["vented", "sealed"].filter(k => k !== o.type) : []);
    const designs = [];
    for (const kind of kinds) {
      if (kind === "vented" && kind !== o.type && !(o.vb || o.f3)) continue;
      const { Vb, fb } = design(kind);
      const QL = qlOf(kind);
      const r = response(p, f, Vb, kind, fb, Deq, QL), rh = response(p, f, Vb, kind, fb, Deq, QL, 2.83, o.hpf || null, o.hpf_order);
      const k = f3Ripple(f, r.spl), kh = f3Ripple(f, rh.spl);
      // max SPL: the least of amplifier, Xmax, driver power and port 17 m/s
      const r1 = response(p, f, Vb, kind, fb, Deq, QL, 1.0, o.hpf || null, o.hpf_order);
      const mx = [], lim = [];
      for (let i = 0; i < f.length; i++) {
        const c = { amplifier: o.amp_v, Xmax: p.Xmax / Math.max(r1.x[i], 1e-12), "driver power": Math.sqrt(o.power * o.znom) };
        if (kind === "vented") c["port 17 m/s"] = 17 / Math.max(r1.vport[i], 1e-12);
        let best = "amplifier"; for (const key in c) if (c[key] < c[best]) best = key;
        mx.push(r1.spl[i] + 20 * Math.log10(c[best])); lim.push(best);
      }
      const e = { type: kind, Vb, fb, net_volume_L: Math.round(Vb * 1e5) / 100, passband_spl: k.ref, f3: k.f3, f3_hpf: kh.f3, hump: k.hump,
        r, rh, mx, lim };
      if (kind === "sealed") { e.fc = p.fs * Math.sqrt(1 + p.Vas / Vb); e.Qtc = p.Qts * Math.sqrt(1 + p.Vas / Vb); }
      else {
        const L = portLength(Vb, fb, Deq / 1000);
        e.port = o.port_shape === "slot" ? `slot ${Math.round(o._Wi)} × ${o.slot_h_eff} mm` : `round ${o.port_d} mm`;
        e.port_eff_mm = L.Leff * 1000; e.port_phys_mm = L.Lphys * 1000; e.pipe_res = C / (2 * L.Leff);
        e.port_speed = Math.max.apply(null, rh.vport) * o.amp_v / 2.83;
        const sz = sizeBox(d, o, Vb * 1000, "vented", e.port_phys_mm);
        e.port_fits = e.port_phys_mm - o.wall + (o.port_shape === "slot" ? o.slot_h_eff : o.port_d / 2) <= sz.Di;
        if (kind === o.type && e.port_phys_mm < o.wall + 10) {   // shorter than the panel it goes through: it can't be built
          const slot = o.port_shape === "slot"; let x = 0.5;
          for (; x < 2000; x += 0.5) { const A = slot ? o._Wi * x : PI * x * x / 4;
            if (portLength(Vb, fb, 2 * Math.sqrt(A / PI) / 1000).Lphys * 1000 >= o.wall + 10) break; }
          e.port_short = true;
          notes.push(`The port only needs to be ${Math.round(e.port_phys_mm)} mm long, shorter than the ${o.wall} mm wall it goes through, so it can't be built: ` +
            `use ${slot ? `a slot at least ${x} mm tall (Slot height)` : `a port at least ${x} mm across`}, or raise fb.`);
        }
        if (!e.port_fits) notes.push(`The ${kind} port is ${Math.round(e.port_phys_mm)} mm long and won't fit in ${Math.round(sz.Di)} mm of inner ` +
          `depth: raise fb, use a ${o.port_shape === "slot" ? "lower slot" : "smaller port diameter"}, or a bigger box.`);
      }
      designs.push(e);
    }
    // the box for the design being built
    const D0 = designs[0];
    const Lp = D0.type === "vented" ? D0.port_phys_mm : 0;
    const sz = sizeBox(d, o, D0.Vb * 1000, D0.type, Lp), t = o.wall;
    const We = sz.Wi + 2 * t, He = sz.Hi + 2 * t, De = sz.Di + 2 * t;
    const box = { inner: [sz.Wi, sz.Hi, sz.Di], outer: [We, He, De], gross: sz.Wi * sz.Hi * sz.Di / 1e6,
      net: sz.Wi * sz.Hi * sz.Di / 1e6 - (d.displacement_l || 0.8) - sz.disp, hw: sz.Hi / sz.Wi, dh: sz.Di / sz.Hi };
    const rch = chamferD(d, o) / 2;
    let yd = o.driver_y !== null ? o.driver_y : Math.min(He / 2 + (sz.Hi / 2 - rch - 2), He - d.outer_d_mm / 2 - 8);
    box.driver_y = yd;
    if (yd + rch > t + sz.Hi + 0.01 || yd - rch < t - 0.01)
      notes.push(`The woofer centre must be between ${Math.round(t + rch)} and ${Math.round(t + sz.Hi - rch)} mm from the bottom for its hole to fit inside.`);
    if (D0.type === "vented") {
      if (o.port_shape === "slot") {
        box.port_on = o.port_on === "front" ? "front" : "back";
        if (box.port_on === "front" && t + o.slot_h_eff + o.shelf + 5 > yd - rch)
          notes.push("The slot and its shelf run into the woofer's hole: a lower slot, the slot on the back, or a taller box.");
      } else box.port_on = o.port_on === "auto" ? "auto" : (o.port_on === "rear" ? "back" : "front");
      box.port_inside = Lp - t;
    }
    // what the drawings need: the box, the woofer and the port, in mm (x across, y up from the bottom, z back from the front)
    const g = { We, He, De, t, Wi: sz.Wi, Hi: sz.Hi, Di: sz.Di, yd, rch, cut: d.cutout_d_mm, od: d.outer_d_mm, frame_t: d.frame_t_mm,
      depth: d.depth_mm, mag: d.magnet_d_mm, screw_n: d.screw_n, screw_bc: d.screw_circle_d_mm, straight: o.straight, cable: 24, port: null };
    if (D0.type === "vented" && o.port_shape === "slot") {
      const R = Math.min(o.flare, t - 5, o.slot_h_eff + o.shelf);
      g.port = { kind: "slot", side: o.port_on === "front" ? "front" : "back", h: o.slot_h_eff, ts: o.shelf, L: Lp, R };
    } else if (D0.type === "vented") {
      const Dout = o.port_d + 2 * o.tube_wall, rfo = Dout / 2 + o.flare, ypf = yd - (rch + Dout / 2 + 5);
      const frontOk = ypf - Dout / 2 >= t + 3 && Math.abs(yd - ypf) >= rch + Dout / 2 + 5 && ypf - rfo >= 8;
      let side = o.port_on === "auto" ? (frontOk ? "front" : "back") : (o.port_on === "front" ? "front" : "back");
      if (side === "front" && !frontOk) notes.push("The round port doesn't fit under the woofer on the front: put it on the back.");
      g.port = { kind: "round", side, yp: side === "front" ? ypf : t + Dout / 2 + 15, Dp: o.port_d, Dout, flare: o.flare, L: Lp };
      box.port_on = side;
    }
    box.geom = g;
    // concrete, roughly: the shell, the back panel, the slot's shelf
    const shell = (We * He * (De - t) - sz.Wi * sz.Hi * (sz.Di)) + We * He * t;
    const slot = D0.type === "vented" && o.port_shape === "slot";
    const shelf = slot ? sz.Wi * o.shelf * (Lp - t) - sz.Wi * o.slot_h_eff * t : 0;   // shelf in, slot opening out
    const rAvg = (d.cutout_d_mm / 2 + rch) / 2;                                        // the woofer's chamfered hole
    const hole = PI * rAvg * rAvg * t + PI * Math.pow(o.cable_d || 24, 2) / 4 * t;
    box.concrete_L = (shell + shelf - hole) / 1e6;
    box.material = raw.material === "wood" ? "wood" : "concrete";
    box.mass_kg = box.concrete_L * (box.material === "wood" ? num(raw.wood_density, 0.68) : 2.2);
    return { f, designs, box, notes, opts: o };
  }

  // ======================================================================= horn
  function hornLive(computeFamily, raw) {
    const o = { family: raw.family || "jmlc", fc: num(raw.fc, 600), T0: num(raw.T0, 0.707), trunc: num(raw.trunc, 180),
      throatD: num(raw.throatD, 50), thickness: num(raw.thickness, 10), thickness_throat: num(raw.thickness_throat, null),
      taper_power: num(raw.taper_power, 1), mode: raw.mode || "mould", skin: num(raw.skin, 1.6), flange_od: num(raw.flange_od, 170),
      flange_t: num(raw.flange_t, 12), bolt_circle: num(raw.bolt_circle, 102), bolt_d: num(raw.bolt_d, 6.5), bolt_n: num(raw.bolt_n, 4),
      bolt_angle: num(raw.bolt_angle, 45), fillet: num(raw.fillet, 6), wall: num(raw.wall, 15), base: num(raw.base, 15) };
    const P = { family: o.family, fc: o.fc, T0: o.T0, trunc: o.trunc, rt: o.throatD / 2, entryAuto: true, petf: false };
    String(raw.hs_extra || "").split(/\s+/).filter(Boolean).forEach(tok => {
      const m = tok.match(/^([A-Za-z]\w*)=(.+)$/); if (m) P[m[1]] = isNaN(+m[2]) ? m[2] : +m[2]; });
    const res = computeFamily(P), notes = [];
    if (res.wallV) throw new Error("this family/setting isn't round: round horns only");
    if (res.terminated === "area-overshoot") notes.push("Horn Studio reports an area overshoot: this geometry isn't validated.");
    const z0 = res.wall[0].z;
    const w = res.wall.map(q => [q.z - z0, q.r]);                         // (z, r)
    const n = w.length;
    // thickness along the curve
    const sArr = [0]; for (let i = 1; i < n; i++) sArr.push(sArr[i - 1] + Math.hypot(w[i][0] - w[i - 1][0], w[i][1] - w[i - 1][1]));
    const tw = sArr.map(s => o.thickness_throat === null ? o.thickness :
      o.thickness_throat + (o.thickness - o.thickness_throat) * Math.pow(s / sArr[n - 1], o.taper_power));
    // outer curve: offset on the left normal (numpy gradient: central inside, one-sided at the ends)
    const off = [];
    for (let i = 0; i < n; i++) {
      const a = w[Math.max(0, i - 1)], b = w[Math.min(n - 1, i + 1)];
      const div2 = (i === 0 || i === n - 1) ? 1 : 2;
      let dz = (b[0] - a[0]) / div2, dr = (b[1] - a[1]) / div2; const L = Math.hypot(dz, dr) || 1; dz /= L; dr /= L;
      off.push([w[i][0] - tw[i] * dr, w[i][1] + tw[i] * dz]);
    }
    if (selfIntersects(off)) notes.push("The wall is thicker than the curve's radius somewhere: its outer surface folds. Use a thinner wall there.");
    const rb = w[0][1], Rf = o.flange_od / 2, tf = o.flange_t;
    const chain = w.map(q => [q[1], q[0]]);                               // (r, z)
    const outer = off.map(q => [q[1], q[0]]);
    const shell = chain.concat(outer.slice().reverse());
    const flange = [[rb, -tf], [Rf, -tf], [Rf, 0], [rb, 0]];
    const fil = o.fillet > 0 ? fillet(outer, o.fillet, rb, Rf) : null;
    // numbers (Pappus)
    let A = 0, M = 0;
    for (const poly of [shell, flange].concat(fil ? [fil.poly] : [])) { const g = areaCentroid(poly); A += g.A; M += g.A * g.cx; }
    const vol = 2 * PI * M / 1000;
    let rmax = Rf, zmax = 0;
    for (const q of shell) { rmax = Math.max(rmax, q[0]); zmax = Math.max(zmax, q[1]); }
    const re_ = chain[n - 1][0], ze = chain[n - 1][1];
    // mould regions (painter's order: J, I, plug, resin, K reproduces the partition)
    const zb = -tf - o.base, Er = rmax + o.wall, ztop = zmax + o.wall;
    const regions = {
      J_jacket: [[Rf, -tf], [Er, -tf], [Er, ztop], [Rf, ztop]],
      I_flange_ring: [[0, -tf], [Rf, -tf], [Rf, ztop], [0, ztop]],
      plug: [[0, zb], [rb, zb], [rb, 0]].concat(chain.slice(1)).concat([[Er, ze], [Er, ztop], [0, ztop]]),
      K_rear_plate: [[0, zb], [Er, zb], [Er, -tf], [0, -tf]] };
    // the horn's outline for a 3-D lathe: plate, fillet, outer curve, lip end, inside curve back to the throat
    let k0 = 0; while (k0 < outer.length - 1 && outer[k0][1] < 0) k0++;
    let lathe = [[rb, -tf], [Rf, -tf], [Rf, 0]];
    if (fil) lathe = lathe.concat(fil.arc).concat(outer.slice(fil.iB));
    else lathe = lathe.concat(outer.slice(k0));
    lathe = lathe.concat(chain.slice().reverse()).concat([[rb, -tf]]);
    return { opts: o, res, chain, outer, shell, flange, fillet: fil, regions, lathe,
      bounds: { zb, Er, ztop, rmax, zmax },
      stats: { mouth_d: 2 * Math.max.apply(null, chain.map(q => q[0])), depth: Math.max.apply(null, chain.map(q => q[1])),
        throat_d: 2 * rb, entry: res.th0Deg, wall: [tw[0], tw[Math.floor(n / 2)], tw[n - 1]], volume_cm3: vol, resin_kg: vol * 1.15 / 1000,
        lip_end: [re_, ze] }, notes };
  }
  function areaCentroid(poly) {   // shoelace in (x = r, y = z)
    let A = 0, cx_ = 0, cy = 0;
    for (let i = 0; i < poly.length; i++) {
      const [x0, y0] = poly[i], [x1, y1] = poly[(i + 1) % poly.length], c = x0 * y1 - x1 * y0;
      A += c; cx_ += (x0 + x1) * c; cy += (y0 + y1) * c;
    }
    A /= 2; return { A: Math.abs(A), cx: A ? cx_ / (6 * A) : 0, cy: A ? cy / (6 * A) : 0 };
  }
  function selfIntersects(pts) {
    const segX = (a, b, c, d) => { const o1 = orient(a, b, c), o2 = orient(a, b, d), o3 = orient(c, d, a), o4 = orient(c, d, b); return o1 * o2 < 0 && o3 * o4 < 0; };
    const orient = (a, b, c) => Math.sign((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]));
    const step = Math.max(1, Math.floor(pts.length / 400));
    const q = pts.filter((_, i) => i % step === 0);
    for (let i = 0; i < q.length - 1; i++) for (let j = i + 2; j < q.length - 1; j++) if (segX(q[i], q[i + 1], q[j], q[j + 1])) return true;
    return false;
  }
  // a fillet of radius f between the plate's front face (z = 0) and the horn's outer curve
  function fillet(outer, f, rb, Rf) {
    const distTo = (P) => { let best = Infinity, bi = 0, bp = null;
      for (let i = 0; i < outer.length - 1; i++) {
        const A = outer[i], B = outer[i + 1], dx = B[0] - A[0], dy = B[1] - A[1], L2 = dx * dx + dy * dy || 1;
        const t = Math.max(0, Math.min(1, ((P[0] - A[0]) * dx + (P[1] - A[1]) * dy) / L2));
        const X = [A[0] + t * dx, A[1] + t * dy], dd = Math.hypot(P[0] - X[0], P[1] - X[1]);
        if (dd < best) { best = dd; bi = i + 1; bp = X; }
      } return { d: best, i: bi, p: bp }; };
    // outer radius at z = f
    let rAt = null;
    for (let i = 0; i < outer.length - 1; i++) if (outer[i][1] <= f && outer[i + 1][1] >= f) {
      const t = (f - outer[i][1]) / (outer[i + 1][1] - outer[i][1]); rAt = outer[i][0] + t * (outer[i + 1][0] - outer[i][0]); break; }
    if (rAt === null) return null;
    let lo = rAt, hi = rAt + 4 * f;
    if (distTo([hi, f]).d < f) return null;
    for (let k = 0; k < 40; k++) { const m = (lo + hi) / 2; if (distTo([m, f]).d < f) lo = m; else hi = m; }
    const rc = (lo + hi) / 2;
    if (rc > Rf || rc < rb + 2) return null;
    const T = distTo([rc, f]);
    const thB = Math.atan2(T.p[1] - f, T.p[0] - rc);   // from the centre to the tangent point on the curve
    const arc = [];
    // from the plate face point (rc, 0) (angle -pi/2) round to the curve's tangent point
    let a0 = -PI / 2, a1 = thB; if (a1 > a0) a1 -= 2 * PI;
    for (let k = 0; k <= 16; k++) { const a = a0 + (a1 - a0) * k / 16; arc.push([rc + f * Math.cos(a), f + f * Math.sin(a)]); }
    // region: plate face from the curve's foot to (rc, 0), the arc, then the curve back down to z = 0
    let k0 = 0; while (k0 < outer.length - 1 && outer[k0][1] < 0) k0++;
    const back = outer.slice(k0, T.i).reverse();
    const poly = arc.concat(back).concat([[outer[k0][0], 0]]);
    return { poly, arc, iB: T.i, rc };
  }

  // ======================================================================= bracket
  function bracketLive(H, B, raw) {
    const h = H.horn, pl = h.plate;
    const [We, He, De] = B.box.outer_mm, wall = B.box.wall_mm;
    const o = { type: raw.type || "L", t: num(raw.t, 5), bend_r: num(raw.bend_r, null), k: num(raw.k, 0.33), width: num(raw.width, null),
      clearance: num(raw.clearance, 10), setback: num(raw.setback, 0), foot_back: num(raw.foot_back, 120), foot_front: num(raw.foot_front, 60),
      hole: num(raw.hole, 6.6), density: num(raw.density, 2.66), driver_d: num(raw.driver_d, 165), driver_depth: num(raw.driver_depth, 70) };
    const t = o.t, R = o.bend_r === null ? t : o.bend_r, W = o.width || pl.od, rt = h.throat_d_mm / 2, notes = [];
    const Ya = He + h.outer_r_mm + o.clearance, Xa = num(raw.xa, We / 2), Zp = o.setback + h.front_z_mm + pl.t;
    if (Xa - W / 2 < 12 || Xa + W / 2 > We - 12) notes.push(`${raw.name || "The bracket"}'s foot runs off the box top (it needs ${Math.round(W)} mm, ` +
      `centred ${Math.round(Xa)} mm from the left of a ${Math.round(We)} mm box): move it with Horizontal shift, widen the box or narrow the bracket.`);
    const Rd = Math.max(W / 2, pl.od / 2, o.driver_d / 2), Ytop = Ya + Rd;
    if (Rd < pl.bolt_circle / 2 + pl.bolt_d) throw new Error("the round plate is too small for the driver's bolt circle");
    const minW = Math.max(2 * o.hole + 16, o.type === "C" ? 3 * (o.hole + 12) : 0);
    if (W < minW) throw new Error(`bracket too narrow for its foot bolts: at least ${Math.round(minW)} mm`);
    if (Zp + t + o.foot_back > De - wall - 10) notes.push("The foot runs onto the back panel's rim: shorten it or move the horn forward.");
    if (o.type === "C" && Zp - o.foot_front < 10) notes.push("The front tabs run off the front of the box.");
    const L1 = Ytop - (He + t + R), BA = PI / 2 * (R + o.k * t), vb0 = L1, vb1 = L1 + BA;
    const feet = o.type === "L" ? [[-W / 2, W / 2, o.foot_back - R, 1]]
      : [[-W / 2, -W / 6, o.foot_front - R, -1], [-W / 6, W / 6, o.foot_back - R, 1], [W / 6, W / 2, o.foot_front - R, -1]];
    const holes = [], inserts = [];
    for (const [u0, u1, Lf, dir] of feet) {
      const wf = u1 - u0, cols = wf > 100 ? [u0 + wf / 4, u1 - wf / 4] : [(u0 + u1) / 2];
      const rows = Lf > 70 ? [R + Math.max(18, 2 * o.hole), Lf + R - Math.max(18, 2 * o.hole)] : [R + Lf / 2 + R / 2];
      for (const c of cols) for (const dd of rows) {
        holes.push([c, vb1 + (dd - R)]);
        inserts.push([Math.round((Xa + c) * 10) / 10, Math.round((dir > 0 ? Zp + t + dd : Zp - dd) * 10) / 10]);
      }
    }
    const vc = Rd;
    let phis;
    if (W / 2 >= Rd - 1e-6) { phis = []; for (let k = 0; k <= 180; k += 3) phis.push((180 - k) * PI / 180); }
    else {
      const sq = Math.sqrt(Rd * Rd - W * W / 4), pL = Math.atan2(-sq, -W / 2) + 2 * PI, pR = Math.atan2(-sq, W / 2);
      const nn = Math.max(24, Math.floor((pL - pR) / (3 * PI / 180))); phis = [];
      for (let k = 0; k <= nn; k++) phis.push(pL - (pL - pR) * k / nn);
      if (vc + sq > L1 - t - 5) throw new Error("the leg is too short below the round plate: more mouth clearance");
    }
    const outline = phis.map(ph => [Rd * Math.cos(ph), -vc + Rd * Math.sin(ph)]);
    const slot = Math.max(t, 2), tabs = feet.slice().sort((a, b) => b[1] - a[1]);
    tabs.forEach(([u0, u1, Lf], i) => {
      const right = i === 0 ? u1 : u1 - slot / 2, left = i === tabs.length - 1 ? u0 : u0 + slot / 2;
      if (i > 0) outline.push([right, -(vb0 - t)]);
      outline.push([right, -(vb1 + Lf)]); outline.push([left, -(vb1 + Lf)]);
      if (i < tabs.length - 1) outline.push([left, -(vb0 - t)]);
    });
    const bolts = [];
    for (let k = 0; k < pl.bolt_n; k++) { const a = (pl.bolt_angle + 360 * k / pl.bolt_n) * PI / 180;
      bolts.push([pl.bolt_circle / 2 * Math.cos(a), -vc + pl.bolt_circle / 2 * Math.sin(a)]); }
    const area = areaCentroid(outline).A - PI * rt * rt - bolts.length * PI * Math.pow(pl.bolt_d / 2, 2) - holes.length * PI * Math.pow(o.hole / 2, 2);
    const minV = Math.min.apply(null, outline.map(q => q[1]));
    return { opts: o, outline, throat: [[0, -vc], rt], bolts: bolts.map(b => [b, pl.bolt_d / 2]), holes: holes.map(q => [[q[0], -q[1]], o.hole / 2]),
      bends: feet.map(([u0, u1, , dir]) => ({ u0, u1, v: -(vb0 + vb1) / 2, dir })), inserts,
      geo: { Ya, Xa, Zp, t, R, Rd, W, Ytop, He, We, De, foot_back: o.foot_back, foot_front: o.foot_front, type: o.type,
        driver_d: o.driver_d, driver_depth: o.driver_depth },
      stats: { type: o.type, t, plate_d: 2 * Rd, width: W, flat: [Math.max(W, 2 * Rd), -minV], mass_kg: area * t / 1e6 * o.density,
        axis_above_box: Ya - He, plate_from_front: Zp, driver_back: Zp + t + o.driver_depth, box_depth: De, bend_allowance: BA }, notes };
  }

  // two horns side by side on the box top: their axes across the box (x from the left, seen from the front),
  // the pair centred on the box, `gap` between the mouths' outsides. r1, r2: the horns' outer radii.
  function pairX(We, r1, r2, gap, side1) {
    const left1 = side1 !== "right", rl = left1 ? r1 : r2, rr = left1 ? r2 : r1, d = rl + gap + rr;
    const xl = (We - d - rr + rl) / 2, xr = xl + d;
    return left1 ? [xl, xr] : [xr, xl];
  }


  // ---- two sealed chambers in one box: a subwoofer below, the midbass above (as boxmould.py dual_size()) ----
  // inside width from the bigger driver; the depth so that each chamber is at least as tall as its driver needs;
  // each chamber's height from its volume. A partition of the wall's thickness between them.
  function dualSize(dm, ds, o, VmL, VsL) {
    const m = 12, Wi = o.width || Math.max(dm.outer_d_mm, ds.outer_d_mm) + 2 * m;
    const gm = VmL + (dm.displacement_l || 0.8), gs = VsL + (ds.displacement_l || 0.8);
    const hm = dm.outer_d_mm + 2 * m, hs = ds.outer_d_mm + 2 * m, need = Math.max(dm.depth_mm, ds.depth_mm) + 40;
    let Di = o.depth || Math.floor(Math.min(gm * 1e6 / (Wi * hm), gs * 1e6 / (Wi * hs)) / 5) * 5;
    Di = Math.max(Di, need);
    const Hm = Math.ceil(gm * 1e6 / (Wi * Di) / 5) * 5, Hs = Math.ceil(gs * 1e6 / (Wi * Di) / 5) * 5;
    return { Wi, Di, Hm, Hs, gm, gs, hm, hs, need };
  }
  function lr4(f, fc, high) {     // Linkwitz-Riley 4th order, complex
    return f.map(q => { const x = cx(0, q / fc), b = inv(add(add(mul(x, x), sc(x, Math.SQRT2)), cx(1)));
      const h = high ? mul(mul(x, x), b) : b; return mul(h, h); });
  }
  function dualLive(dm, ds, raw, nf) {
    const t = num(raw.wall, 18), xo = num(raw.xo_lo, 90), notes = [];
    const mid = boxLive(dm, { type: "sealed", vb: raw.vb, qtc: raw.qtc, amp_v: raw.amp_v, hpf: 0, wall: t, ql: raw.ql }, nf);
    const sub = boxLive(ds, { type: "sealed", vb: raw.sub_vb, qtc: raw.sub_qtc, amp_v: raw.sub_amp_v || 25, hpf: 0, wall: t, ql: raw.ql }, nf);
    const em = mid.designs[0], es = sub.designs[0], f = mid.f;
    // the DSP: LR4 between them at xo; the sub's gain sets its level to the midbass's (at 3 x below the crossover)
    const LP = lr4(f, xo, false), HP = lr4(f, xo, true), fl = xo / 3;
    // the sub's bass boost (a Linkwitz transform in the DSP): its sealed roll-off (fc, Qtc) moved down to f0, Q 0.707
    const f0 = num(raw.sub_lt, null);
    const LT = f.map(q => { if (!f0 || f0 >= es.fc) return cx(1);
      const s_ = cx(0, q), w1 = 2 * PI * es.fc, w2 = 2 * PI * f0, sq = z => mul(z, z);
      return div(add(add(sq(sc(s_, 2 * PI)), sc(s_, 2 * PI * w1 / es.Qtc)), cx(w1 * w1)), add(add(sq(sc(s_, 2 * PI)), sc(s_, 2 * PI * w2 / 0.707)), cx(w2 * w2))); });
    const boost_db = Math.max(...LT.map(z => 20 * Math.log10(abs(z))));
    const at = (ys, x) => { for (let i = 1; i < f.length; i++) if (f[i] >= x) return ys[i - 1] + (ys[i] - ys[i - 1]) * (x - f[i - 1]) / (f[i] - f[i - 1]); return ys[ys.length - 1]; };
    const subEq = es.r.pc.map((p, i) => mul(p, LT[i])).map(p => 20 * Math.log10(Math.max(abs(p), 1e-12) / P0));
    const gain_db = em.passband_spl - at(subEq, Math.max(fl, 1.6 * (f0 && f0 < es.fc ? f0 : es.fc)));
    const g = Math.pow(10, gain_db / 20);
    const cm = em.r.pc.map((p, i) => mul(p, HP[i])), cs = es.r.pc.map((p, i) => sc(mul(mul(p, LT[i]), LP[i]), g));
    const sum = cm.map((p, i) => add(p, cs[i]));
    const dB = p => 20 * Math.log10(Math.max(abs(p), 1e-12) / P0);
    const spl = sum.map(dB), splM = cm.map(dB), splS = cs.map(dB);
    // the most the pair can play: each driver hits its own limit where it carries its share of the sum
    // the most the pair can play: each driver's own output is its share of the sum; it stops at its own limit
    // (Xmax, amplifier, power: em.mx / es.mx, which no EQ changes)
    const mx = f.map((q, i) => { const a = abs(sum[i]) || 1e-12;
      return Math.min(em.mx[i] - 20 * Math.log10(Math.max(abs(cm[i]) / a, 1e-9)), es.mx[i] - 20 * Math.log10(Math.max(abs(cs[i]) / a, 1e-9))); });
    const k = f3Ripple(f, spl);
    // the box
    const sz = dualSize(dm, ds, { width: num(raw.width, null), depth: num(raw.depth, null) }, em.net_volume_L, es.net_volume_L);
    const We = sz.Wi + 2 * t, He = sz.Hs + t + sz.Hm + 2 * t, De = sz.Di + 2 * t;
    if (sz.Hm < sz.hm - 0.5) notes.push(`The midbass's chamber is ${Math.round(sz.Hm)} mm tall inside, less than the ${Math.round(sz.hm)} mm its driver needs: more volume, or a shallower box.`);
    if (sz.Hs < sz.hs - 0.5) notes.push(`The sub's chamber is ${Math.round(sz.Hs)} mm tall inside, less than the ${Math.round(sz.hs)} mm its driver needs: more volume, or a shallower box.`);
    const ys = t + sz.Hs / 2, yp0 = t + sz.Hs, ym = yp0 + t + sz.Hm / 2;
    const drvG = (d, yd) => ({ yd, rch: d.cutout_d_mm / 2, cut: d.cutout_d_mm, od: d.outer_d_mm, frame_t: d.frame_t_mm, depth: d.depth_mm, mag: d.magnet_d_mm,
      screw_n: d.screw_n, screw_bc: d.screw_circle_d_mm });
    const geom = Object.assign({ We, He, De, t, Wi: sz.Wi, Hi: He - 2 * t, Di: sz.Di, straight: 10, cable: 24, port: null,
      dual: { partition: [yp0, yp0 + t], mid: drvG(dm, ym), sub: drvG(ds, ys), Hm: sz.Hm, Hs: sz.Hs } }, drvG(dm, ym));
    const inner = [sz.Wi, He - 2 * t, sz.Di];
    const wood_L = (We * He * De - sz.Wi * sz.Di * (sz.Hm + sz.Hs) - PI * t / 4 * (dm.cutout_d_mm ** 2 + ds.cutout_d_mm ** 2)) / 1e6;
    const box = { inner, outer: [We, He, De], gross: (sz.Hm + sz.Hs) * sz.Wi * sz.Di / 1e6,
      net: em.net_volume_L + es.net_volume_L, hw: inner[1] / inner[0], dh: inner[2] / inner[1], driver_y: ym, geom,
      concrete_L: wood_L, material: raw.material === "wood" ? "wood" : "concrete" };
    box.mass_kg = wood_L * (box.material === "wood" ? num(raw.wood_density, 0.68) : 2.2);
    notes.push(...mid.notes.filter(n => !/woofer centre/.test(n)), ...sub.notes.filter(n => !/woofer centre/.test(n)));
    return { f, mid: em, sub: es, opts: mid.opts, designs: [Object.assign({ role: "midbass" }, em), Object.assign({ role: "sub" }, es)],
      dual: { xo, gain_db, boost_db, f0, spl, splM, splS, mx, f3: k.f3, ref: k.ref, hump: k.hump, Hm: sz.Hm, Hs: sz.Hs }, box, notes };
  }

  const api = { dualLive, dualSize, pairX, response, f3Ripple, bestFb, portLength, sizeBox, boxLive, hornLive, bracketLive, areaCentroid, logspace };
  if (typeof module !== "undefined") module.exports = api; else root.Live = api;
})(this);
