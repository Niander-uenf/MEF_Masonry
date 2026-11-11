/* 3D post-processor in plain WebGL (no external libraries).
   Reads the visual.json written by mef_masonry: mesh (6-node wedges, hexahedra or triangles),
   material and states with nodal or element fields and displacements. */
(function () {
  "use strict";

  const COLOUR_MAPS = {
    viridis: ["#440154", "#482878", "#3e4989", "#31688e", "#26828e", "#1f9e89", "#35b779", "#6ece58", "#b5de2b", "#fde725"],
    inferno: ["#000004", "#1b0c41", "#4a0c6b", "#781c6d", "#a52c60", "#cf4446", "#ed6925", "#fb9b06", "#f7d13d", "#fcffa4"],
    diverging: ["#3b4cc0", "#5977e3", "#7b9ff9", "#9ebeff", "#c0d4f5", "#dddcdc", "#f2cbb7", "#f7ac8e", "#ee8468", "#d65244", "#b40426"],
    damage: ["#fff8e6", "#fde3c8", "#fbc89c", "#f6a06d", "#ec7647", "#d9492f", "#b5222c", "#7e0b27"],
    grey: ["#2b2b2b", "#f0f0f0"],
  };
  const MATERIAL_COLOUR = { 0: "#cde2fb", 1: "#d6b58f", 2: "#eb6834", 3: "#eda100", 4: "#8a8984", 5: "#7f8fa6", 6: "#5b5f66" };
  const NO_DATA = [0.72, 0.76, 0.82];   // bluish grey: no data (for example, damage outside the ceramic)

  function hex(c) { return [parseInt(c.slice(1, 3), 16) / 255, parseInt(c.slice(3, 5), 16) / 255, parseInt(c.slice(5, 7), 16) / 255]; }
  const LUT = {};
  for (const k in COLOUR_MAPS) LUT[k] = COLOUR_MAPS[k].map(hex);
  function mapColour(name, t) {
    const s = LUT[name] || LUT.viridis;
    if (!(t >= 0)) t = 0; if (t > 1) t = 1;
    const x = t * (s.length - 1), i = Math.min(Math.floor(x), s.length - 2), f = x - i;
    return [s[i][0] + (s[i + 1][0] - s[i][0]) * f, s[i][1] + (s[i + 1][1] - s[i][1]) * f, s[i][2] + (s[i + 1][2] - s[i][2]) * f];
  }

  // --- 4x4 matrices (column-major, as WebGL expects)
  function perspective(fov, asp, near, far) {
    const f = 1 / Math.tan(fov / 2), nf = 1 / (near - far);
    return [f / asp, 0, 0, 0, 0, f, 0, 0, 0, 0, (far + near) * nf, -1, 0, 0, 2 * far * near * nf, 0];
  }
  function sub(a, b) { return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]; }
  function cross(a, b) { return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]; }
  function norm(a) { const l = Math.hypot(a[0], a[1], a[2]) || 1; return [a[0] / l, a[1] / l, a[2] / l]; }
  function lookAt(eye, target, up) {
    const z = norm(sub(eye, target)), x = norm(cross(up, z)), y = cross(z, x);
    return [x[0], y[0], z[0], 0, x[1], y[1], z[1], 0, x[2], y[2], z[2], 0,
      -(x[0] * eye[0] + x[1] * eye[1] + x[2] * eye[2]), -(y[0] * eye[0] + y[1] * eye[1] + y[2] * eye[2]),
      -(z[0] * eye[0] + z[1] * eye[1] + z[2] * eye[2]), 1];
  }
  function mult(a, b) {
    const r = new Array(16);
    for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) {
      let s = 0; for (let k = 0; k < 4; k++) s += a[k * 4 + j] * b[i * 4 + k]; r[i * 4 + j] = s;
    }
    return r;
  }

  const VS = `attribute vec3 p; attribute vec3 n; attribute vec3 c; uniform mat4 mvp; varying vec3 vn; varying vec3 vc;
    void main(){ vn = n; vc = c; gl_Position = mvp * vec4(p, 1.0); }`;
  const FS = `precision mediump float; varying vec3 vn; varying vec3 vc; uniform vec3 light;
    void main(){ float l = 0.42 + 0.58 * abs(dot(normalize(vn), light)); gl_FragColor = vec4(vc * l, 1.0); }`;
  const VSL = `attribute vec3 p; uniform mat4 mvp; void main(){ gl_Position = mvp * vec4(p, 1.0); }`;
  const FSL = `precision mediump float; uniform vec4 colour; void main(){ gl_FragColor = colour; }`;

  function compile(gl, vs, fs) {
    const pr = gl.createProgram();
    for (const [kind, src] of [[gl.VERTEX_SHADER, vs], [gl.FRAGMENT_SHADER, fs]]) {
      const s = gl.createShader(kind); gl.shaderSource(s, src); gl.compileShader(s);
      if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s));
      gl.attachShader(pr, s);
    }
    gl.linkProgram(pr);
    if (!gl.getProgramParameter(pr, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(pr));
    return pr;
  }

  const FACES_W6 = [[0, 2, 1], [3, 4, 5], [0, 1, 4, 3], [1, 2, 5, 4], [2, 0, 3, 5]];

  class Viewer {
    constructor(canvas) {
      this.canvas = canvas;
      const gl = canvas.getContext("webgl", { antialias: true, preserveDrawingBuffer: true }) ||
        canvas.getContext("experimental-webgl", { preserveDrawingBuffer: true });
      if (!gl) throw new Error("WebGL is not available in this browser.");
      this.gl = gl;
      this.pf = compile(gl, VS, FS);
      this.pl = compile(gl, VSL, FSL);
      this.buf = { p: gl.createBuffer(), n: gl.createBuffer(), c: gl.createBuffer(), l: gl.createBuffer() };
      this.op = { state: 0, field: null, cmap: "auto", auto: true, min: 0, max: 1, deformed: false, scale: 50,
        mirror: true, edges: true, cuts: { x: null, y: null, z: null }, hidden: new Set([0]) };
      this.cam = { az: -0.65, el: 0.42, dist: 1, target: [0, 0, 0], pan: [0, 0] };
      this.data = null;
      this.nTri = 0; this.nLin = 0;
      this._events();
      new ResizeObserver(() => this.draw()).observe(canvas);
    }

    load(d) {
      this.data = d;
      const n = d.nodes.length / 3;
      this.X = Float32Array.from(d.nodes);
      this.E = Int32Array.from(d.elements);
      this.npe = d.nodes_per_element;
      this.elFaces = d.faces || FACES_W6;      // faces of the element (local nodes), W6 or H8
      this.mat = Int8Array.from(d.material);
      const ne = this.mat.length, npe = this.npe;
      this.C = new Float32Array(ne * 3);
      for (let e = 0; e < ne; e++) for (let k = 0; k < npe; k++) {
        const i = this.E[e * npe + k];
        for (let a = 0; a < 3; a++) this.C[e * 3 + a] += this.X[i * 3 + a] / npe;
      }
      const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
      for (let i = 0; i < n; i++) for (let a = 0; a < 3; a++) { const v = this.X[i * 3 + a]; if (v < lo[a]) lo[a] = v; if (v > hi[a]) hi[a] = v; }
      this.lo = lo; this.hi = hi;
      this.nNodes = n;
      // fields: indices and conversion to Float32 (null -> NaN)
      this.fields = d.fields;
      this.states = d.states.map(e => {
        const v = {};
        for (const k in e.values) v[k] = Float32Array.from(e.values[k], x => x === null ? NaN : x);
        return { label: e.label, v, u: e.u ? Float32Array.from(e.u, x => x === null ? 0 : x) : null };
      });
      this.op.state = Math.max(0, this.states.length - 1);
      this.op.field = this.fields.length ? this.fields[0].id : null;
      this.op.hidden = new Set([0]);
      if (d.kind === "3d" && this.states.some(e => e.u)) this.op.hidden.add(4);   // steel outside the mechanical model
      this.view(d.kind === "2d" ? "top" : "iso");
      this.updateGeometry();
    }

    fieldValues(id, ist) {
      if (!this.states.length) return null;
      const e = this.states[ist];
      if (e && e.v[id]) return e.v[id];
      for (const o of this.states) if (o.v[id]) return o.v[id];
      return null;
    }

    fieldInfo(id) { return this.fields.find(c => c.id === id) || { id, name: id, unit: "", local: "node", cmap: "viridis" }; }

    _mirrored() { return !!(this.data && this.data.symmetry && this.op.mirror && this.data.kind === "3d"); }

    _limits() {
      const lo = this.lo.slice(), hi = this.hi.slice();
      if (this._mirrored()) hi[1] = 2 * this.data.model_width;
      return [lo, hi];
    }

    _visible(copy) {
      const ne = this.mat.length, vis = new Uint8Array(ne), c = this.op.cuts, [lo, hi] = this._limits();
      const lim = [c.x, c.y, c.z].map((f, a) => f === null ? null : lo[a] + (hi[a] - lo[a]) * f);
      const W2 = 2 * this.data.model_width;
      for (let e = 0; e < ne; e++) {
        if (this.op.hidden.has(this.mat[e])) continue;
        let ok = true;
        for (let a = 0; a < 3; a++) {
          if (lim[a] === null) continue;
          const x = (a === 1 && copy) ? W2 - this.C[e * 3 + 1] : this.C[e * 3 + a];
          if (x > lim[a] + 1e-9) { ok = false; break; }
        }
        if (ok) vis[e] = 1;
      }
      return vis;
    }

    _faces(vis) {
      const npe = this.npe, E = this.E, n = this.nNodes, faces = [];
      if (npe === 3) {
        for (let e = 0; e < vis.length; e++) if (vis[e]) faces.push([e, E[e * 3], E[e * 3 + 1], E[e * 3 + 2]]);
        return faces;
      }
      const seen = new Map();
      for (let e = 0; e < vis.length; e++) {
        if (!vis[e]) continue;
        for (const f of this.elFaces) {
          const ids = f.map(k => E[e * npe + k]);
          const s = ids.slice().sort((a, b) => a - b);
          const key = s[0] + n * (s[1] + n * s[2]);
          const prev = seen.get(key);
          if (prev === undefined) { seen.set(key, faces.length); faces.push([e, ...ids]); }
          else if (prev >= 0) { faces[prev] = null; seen.set(key, -1); }
        }
      }
      return faces.filter(f => f);
    }

    updateGeometry() {
      if (!this.data) return;
      const copies = this._mirrored() ? 2 : 1;
      const groups = [];
      for (let cp = 0; cp < copies; cp++) groups.push(this._faces(this._visible(cp)));
      let nt = 0, nl = 0;
      for (const g of groups) for (const f of g) { nt += f.length === 5 ? 2 : 1; nl += f.length - 1; }
      this.vNode = new Int32Array(nt * 3); this.vEl = new Int32Array(nt * 3); this.vSym = new Uint8Array(nt * 3);
      this.lNode = new Int32Array(nl * 2); this.lSym = new Uint8Array(nl * 2);
      let iv = 0, il = 0;
      groups.forEach((g, cp) => { for (const f of g) {
        const e = f[0], ids = f.slice(1);
        const tris = ids.length === 4 ? [[0, 1, 2], [0, 2, 3]] : [[0, 1, 2]];
        for (const t of tris) for (const k of t) { this.vNode[iv] = ids[k]; this.vEl[iv] = e; this.vSym[iv] = cp; iv++; }
        for (let k = 0; k < ids.length; k++) {
          this.lNode[il] = ids[k]; this.lSym[il++] = cp;
          this.lNode[il] = ids[(k + 1) % ids.length]; this.lSym[il++] = cp;
        }
      } });
      this.nTri = nt; this.nLin = nl;
      this.updatePositions();
      this.updateColours();
    }

    _pos(i, sym, out, o) {
      const X = this.X, d = this.op.deformed ? this.states[this.op.state] : null;
      let x = X[i * 3], y = X[i * 3 + 1], z = X[i * 3 + 2];
      if (d && d.u) { const s = this.op.scale; x += s * d.u[i * 3]; y += s * d.u[i * 3 + 1]; z += s * d.u[i * 3 + 2]; }
      if (sym) y = 2 * this.data.model_width - y;
      out[o] = x; out[o + 1] = y; out[o + 2] = z;
    }

    updatePositions() {
      if (!this.data) return;
      const gl = this.gl, nv = this.nTri * 3;
      const P = new Float32Array(nv * 3), N = new Float32Array(nv * 3);
      for (let v = 0; v < nv; v++) this._pos(this.vNode[v], this.vSym[v], P, v * 3);
      for (let t = 0; t < this.nTri; t++) {
        const o = t * 9;
        const a = [P[o + 3] - P[o], P[o + 4] - P[o + 1], P[o + 5] - P[o + 2]];
        const b = [P[o + 6] - P[o], P[o + 7] - P[o + 1], P[o + 8] - P[o + 2]];
        const nn = norm(cross(a, b));
        for (let k = 0; k < 3; k++) { N[o + k * 3] = nn[0]; N[o + k * 3 + 1] = nn[1]; N[o + k * 3 + 2] = nn[2]; }
      }
      this.P = P;
      gl.bindBuffer(gl.ARRAY_BUFFER, this.buf.p); gl.bufferData(gl.ARRAY_BUFFER, P, gl.STATIC_DRAW);
      gl.bindBuffer(gl.ARRAY_BUFFER, this.buf.n); gl.bufferData(gl.ARRAY_BUFFER, N, gl.STATIC_DRAW);
      const L = new Float32Array(this.nLin * 2 * 3);
      for (let k = 0; k < this.nLin * 2; k++) this._pos(this.lNode[k], this.lSym[k], L, k * 3);
      gl.bindBuffer(gl.ARRAY_BUFFER, this.buf.l); gl.bufferData(gl.ARRAY_BUFFER, L, gl.STATIC_DRAW);
      this.draw();
    }

    range() {
      const id = this.op.field;
      if (!id || id === "material") return [0, 1];
      const v = this.fieldValues(id, this.op.state);
      if (!v) return [0, 1];
      const info = this.fieldInfo(id), perEl = info.local === "element";
      // range of what is visible only
      let lo = Infinity, hi = -Infinity;
      const nv = this.nTri * 3;
      for (let k = 0; k < nv; k++) {
        const x = v[perEl ? this.vEl[k] : this.vNode[k]];
        if (x === x) { if (x < lo) lo = x; if (x > hi) hi = x; }
      }
      if (!isFinite(lo)) return [0, 1];
      if (this.effectiveMap() === "diverging" && lo < 0 && hi > 0) { const m = Math.max(-lo, hi); lo = -m; hi = m; }
      if (info.cmap === "damage") { lo = 0; hi = Math.max(hi, 1e-6) > 0.5 ? 1 : hi; }
      if (hi - lo < 1e-12) hi = lo + 1e-12;
      return [lo, hi];
    }

    effectiveMap() {
      if (this.op.cmap !== "auto") return this.op.cmap;
      const m = this.fieldInfo(this.op.field).cmap;
      return COLOUR_MAPS[m] ? m : "viridis";
    }

    updateColours() {
      if (!this.data) return;
      const gl = this.gl, nv = this.nTri * 3, Cc = new Float32Array(nv * 3), id = this.op.field;
      if (this.op.auto) { const [a, b] = this.range(); this.op.min = a; this.op.max = b; }
      const lo = this.op.min, hi = this.op.max, cmap = this.effectiveMap();
      const v = id && id !== "material" ? this.fieldValues(id, this.op.state) : null;
      const perEl = this.fieldInfo(id).local === "element";
      for (let k = 0; k < nv; k++) {
        let c;
        if (!v) c = hex(MATERIAL_COLOUR[this.mat[this.vEl[k]]] || "#bbbbbb");
        else {
          const x = v[perEl ? this.vEl[k] : this.vNode[k]];
          c = x === x ? mapColour(cmap, (x - lo) / (hi - lo)) : NO_DATA;
        }
        Cc[k * 3] = c[0]; Cc[k * 3 + 1] = c[1]; Cc[k * 3 + 2] = c[2];
      }
      gl.bindBuffer(gl.ARRAY_BUFFER, this.buf.c); gl.bufferData(gl.ARRAY_BUFFER, Cc, gl.STATIC_DRAW);
      this.draw();
      if (this.onScaleChange) this.onScaleChange();
    }

    cssGradient() {
      const s = COLOUR_MAPS[this.effectiveMap()] || COLOUR_MAPS.viridis;
      return "linear-gradient(90deg," + s.join(",") + ")";
    }

    view(name) {
      const c = this.cam;
      if (name === "iso") { c.az = -0.65; c.el = 0.42; }
      else if (name === "front") { c.az = -Math.PI / 2; c.el = 0.0; }
      else if (name === "side") { c.az = 0; c.el = 0.0; }
      else if (name === "top") { c.az = -Math.PI / 2; c.el = Math.PI / 2 - 1e-3; }
      c.pan = [0, 0];
      if (this.lo && this.data) {
        const [lo, hi] = this._limits();
        c.target = [(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, (lo[2] + hi[2]) / 2];
        c.dist = 2.6 * Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]);
      }
      this.draw();
    }

    _mvp() {
      const c = this.cam, cv = this.canvas;
      const asp = cv.width / Math.max(cv.height, 1);
      const dir = [Math.cos(c.el) * Math.cos(c.az), Math.cos(c.el) * Math.sin(c.az), Math.sin(c.el)];
      const up = Math.abs(c.el) > 1.5 ? [0, 1, 0] : [0, 0, 1];
      // lateral shift (pan) in the plane of the screen
      const dirn = norm(dir), dx = norm(cross(up, dirn)), dy = cross(dirn, dx);
      const target = [c.target[0] + dx[0] * c.pan[0] + dy[0] * c.pan[1], c.target[1] + dx[1] * c.pan[0] + dy[1] * c.pan[1],
        c.target[2] + dx[2] * c.pan[0] + dy[2] * c.pan[1]];
      const eye = [target[0] + dir[0] * c.dist, target[1] + dir[1] * c.dist, target[2] + dir[2] * c.dist];
      const V = lookAt(eye, target, up);
      const Pm = perspective(30 * Math.PI / 180, asp, c.dist / 100, c.dist * 10);
      return { mvp: mult(Pm, V), light: dirn };
    }

    draw() {
      const gl = this.gl, cv = this.canvas, dpr = window.devicePixelRatio || 1;
      const w = Math.max(1, Math.round(cv.clientWidth * dpr)), h = Math.max(1, Math.round(cv.clientHeight * dpr));
      if (cv.width !== w || cv.height !== h) { cv.width = w; cv.height = h; }
      gl.viewport(0, 0, w, h);
      gl.clearColor(0, 0, 0, 0); gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
      if (!this.data || !this.nTri) return;
      gl.enable(gl.DEPTH_TEST);
      const { mvp, light } = this._mvp();
      this.lastMVP = mvp;
      gl.useProgram(this.pf);
      gl.enable(gl.POLYGON_OFFSET_FILL); gl.polygonOffset(1, 1);
      for (const [name, b] of [["p", this.buf.p], ["n", this.buf.n], ["c", this.buf.c]]) {
        const loc = gl.getAttribLocation(this.pf, name);
        gl.bindBuffer(gl.ARRAY_BUFFER, b); gl.enableVertexAttribArray(loc); gl.vertexAttribPointer(loc, 3, gl.FLOAT, false, 0, 0);
      }
      gl.uniformMatrix4fv(gl.getUniformLocation(this.pf, "mvp"), false, new Float32Array(mvp));
      gl.uniform3fv(gl.getUniformLocation(this.pf, "light"), new Float32Array(light));
      gl.drawArrays(gl.TRIANGLES, 0, this.nTri * 3);
      gl.disable(gl.POLYGON_OFFSET_FILL);
      for (const name of ["n", "c"]) gl.disableVertexAttribArray(gl.getAttribLocation(this.pf, name));
      if (this.op.edges && this.nLin) {
        gl.useProgram(this.pl);
        const loc = gl.getAttribLocation(this.pl, "p");
        gl.bindBuffer(gl.ARRAY_BUFFER, this.buf.l); gl.enableVertexAttribArray(loc); gl.vertexAttribPointer(loc, 3, gl.FLOAT, false, 0, 0);
        gl.uniformMatrix4fv(gl.getUniformLocation(this.pl, "mvp"), false, new Float32Array(mvp));
        gl.uniform4fv(gl.getUniformLocation(this.pl, "colour"), new Float32Array([0.1, 0.1, 0.1, 0.35]));
        gl.enable(gl.BLEND); gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
        gl.drawArrays(gl.LINES, 0, this.nLin * 2);
        gl.disable(gl.BLEND);
      }
    }

    probe(px, py) {
      if (!this.P || !this.lastMVP) return null;
      const m = this.lastMVP, cv = this.canvas, dpr = window.devicePixelRatio || 1;
      const X = px * dpr, Y = py * dpr, P = this.P, nv = this.nTri * 3;
      let best = -1, dmin = Infinity, front = -1, zmin = Infinity;
      for (let k = 0; k < nv; k++) {
        const x = P[k * 3], y = P[k * 3 + 1], z = P[k * 3 + 2];
        const cw = m[3] * x + m[7] * y + m[11] * z + m[15];
        const sx = ((m[0] * x + m[4] * y + m[8] * z + m[12]) / cw * 0.5 + 0.5) * cv.width;
        const sy = (1 - ((m[1] * x + m[5] * y + m[9] * z + m[13]) / cw * 0.5 + 0.5)) * cv.height;
        const sz = (m[2] * x + m[6] * y + m[10] * z + m[14]) / cw;
        const d = Math.hypot(sx - X, sy - Y);
        if (d < dmin) { dmin = d; best = k; }
        if (d < 6 * dpr && sz < zmin) { zmin = sz; front = k; }
      }
      if (front >= 0) best = front;
      else if (dmin > 14 * dpr) best = -1;
      if (best < 0) return null;
      const i = this.vNode[best], e = this.vEl[best], id = this.op.field;
      const info = this.fieldInfo(id);
      const v = id && id !== "material" ? this.fieldValues(id, this.op.state) : null;
      const val = v ? v[info.local === "element" ? e : i] : null;
      const nm = this.data.material_names[this.mat[e]] || this.mat[e];
      return { node: i, element: e, material: nm, x: this.X[i * 3], y: this.X[i * 3 + 1], z: this.X[i * 3 + 2],
        value: val, field: info.name, unit: info.unit };
    }

    png() { this.draw(); return this.canvas.toDataURL("image/png"); }

    _events() {
      const cv = this.canvas; let drag = null;
      cv.addEventListener("contextmenu", e => e.preventDefault());
      cv.addEventListener("mousedown", e => { drag = { x: e.clientX, y: e.clientY, pan: e.shiftKey || e.button === 2 }; });
      window.addEventListener("mouseup", () => { drag = null; });
      window.addEventListener("mousemove", e => {
        if (!drag) return;
        const dx = e.clientX - drag.x, dy = e.clientY - drag.y; drag.x = e.clientX; drag.y = e.clientY;
        if (drag.pan) { const s = this.cam.dist / 900; this.cam.pan[0] -= dx * s; this.cam.pan[1] += dy * s; }
        else { this.cam.az -= dx * 0.008; this.cam.el = Math.max(-1.55, Math.min(1.55, this.cam.el + dy * 0.008)); }
        this.draw();
      });
      cv.addEventListener("wheel", e => { e.preventDefault(); this.cam.dist *= Math.exp(e.deltaY * 0.0012); this.draw(); }, { passive: false });
      let touch = null;
      cv.addEventListener("touchstart", e => { const t = e.touches[0]; touch = { x: t.clientX, y: t.clientY }; }, { passive: true });
      cv.addEventListener("touchmove", e => {
        if (!touch) return; const t = e.touches[0];
        this.cam.az -= (t.clientX - touch.x) * 0.008; this.cam.el = Math.max(-1.55, Math.min(1.55, this.cam.el + (t.clientY - touch.y) * 0.008));
        touch = { x: t.clientX, y: t.clientY }; this.draw();
      }, { passive: true });
    }
  }

  window.Viewer = Viewer;
  window.COLOUR_MAPS = COLOUR_MAPS;
  window.MATERIAL_COLOUR = MATERIAL_COLOUR;
})();
