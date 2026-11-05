/* Pós-processador 3D em WebGL puro (sem bibliotecas externas).
   Lê o visual.json gravado pelo mef_masonry: malha (cunhas de 6 nós ou triângulos),
   material e estados com campos nodais ou por elemento e deslocamentos. */
(function () {
  "use strict";

  const MAPAS = {
    viridis: ["#440154", "#482878", "#3e4989", "#31688e", "#26828e", "#1f9e89", "#35b779", "#6ece58", "#b5de2b", "#fde725"],
    inferno: ["#000004", "#1b0c41", "#4a0c6b", "#781c6d", "#a52c60", "#cf4446", "#ed6925", "#fb9b06", "#f7d13d", "#fcffa4"],
    divergente: ["#3b4cc0", "#5977e3", "#7b9ff9", "#9ebeff", "#c0d4f5", "#dddcdc", "#f2cbb7", "#f7ac8e", "#ee8468", "#d65244", "#b40426"],
    dano: ["#fff8e6", "#fde3c8", "#fbc89c", "#f6a06d", "#ec7647", "#d9492f", "#b5222c", "#7e0b27"],
    cinza: ["#2b2b2b", "#f0f0f0"],
  };
  const COR_MATERIAL = { 0: "#cde2fb", 1: "#d6b58f", 2: "#eb6834", 3: "#eda100", 4: "#8a8984", 5: "#7f8fa6", 6: "#5b5f66" };
  const SEM_DADO = [0.72, 0.76, 0.82];   // cinza azulado: sem dado (por exemplo, dano fora da cerâmica)

  function hex(c) { return [parseInt(c.slice(1, 3), 16) / 255, parseInt(c.slice(3, 5), 16) / 255, parseInt(c.slice(5, 7), 16) / 255]; }
  const LUT = {};
  for (const k in MAPAS) LUT[k] = MAPAS[k].map(hex);
  function corMapa(nome, t) {
    const s = LUT[nome] || LUT.viridis;
    if (!(t >= 0)) t = 0; if (t > 1) t = 1;
    const x = t * (s.length - 1), i = Math.min(Math.floor(x), s.length - 2), f = x - i;
    return [s[i][0] + (s[i + 1][0] - s[i][0]) * f, s[i][1] + (s[i + 1][1] - s[i][1]) * f, s[i][2] + (s[i + 1][2] - s[i][2]) * f];
  }

  // --- matrizes 4x4 (coluna maior, como o WebGL espera)
  function perspectiva(fov, asp, perto, longe) {
    const f = 1 / Math.tan(fov / 2), nf = 1 / (perto - longe);
    return [f / asp, 0, 0, 0, 0, f, 0, 0, 0, 0, (longe + perto) * nf, -1, 0, 0, 2 * longe * perto * nf, 0];
  }
  function sub(a, b) { return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]; }
  function cruz(a, b) { return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]; }
  function norm(a) { const l = Math.hypot(a[0], a[1], a[2]) || 1; return [a[0] / l, a[1] / l, a[2] / l]; }
  function olharPara(olho, alvo, cima) {
    const z = norm(sub(olho, alvo)), x = norm(cruz(cima, z)), y = cruz(z, x);
    return [x[0], y[0], z[0], 0, x[1], y[1], z[1], 0, x[2], y[2], z[2], 0,
      -(x[0] * olho[0] + x[1] * olho[1] + x[2] * olho[2]), -(y[0] * olho[0] + y[1] * olho[1] + y[2] * olho[2]),
      -(z[0] * olho[0] + z[1] * olho[1] + z[2] * olho[2]), 1];
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
  const FS = `precision mediump float; varying vec3 vn; varying vec3 vc; uniform vec3 luz;
    void main(){ float l = 0.42 + 0.58 * abs(dot(normalize(vn), luz)); gl_FragColor = vec4(vc * l, 1.0); }`;
  const VSL = `attribute vec3 p; uniform mat4 mvp; void main(){ gl_Position = mvp * vec4(p, 1.0); }`;
  const FSL = `precision mediump float; uniform vec4 cor; void main(){ gl_FragColor = cor; }`;

  function compilar(gl, vs, fs) {
    const pr = gl.createProgram();
    for (const [tipo, src] of [[gl.VERTEX_SHADER, vs], [gl.FRAGMENT_SHADER, fs]]) {
      const s = gl.createShader(tipo); gl.shaderSource(s, src); gl.compileShader(s);
      if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s));
      gl.attachShader(pr, s);
    }
    gl.linkProgram(pr);
    if (!gl.getProgramParameter(pr, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(pr));
    return pr;
  }

  const FACES_W6 = [[0, 2, 1], [3, 4, 5], [0, 1, 4, 3], [1, 2, 5, 4], [2, 0, 3, 5]];

  class Visualizador {
    constructor(canvas) {
      this.canvas = canvas;
      const gl = canvas.getContext("webgl", { antialias: true, preserveDrawingBuffer: true }) ||
        canvas.getContext("experimental-webgl", { preserveDrawingBuffer: true });
      if (!gl) throw new Error("WebGL indisponível neste navegador.");
      this.gl = gl;
      this.pf = compilar(gl, VS, FS);
      this.pl = compilar(gl, VSL, FSL);
      this.buf = { p: gl.createBuffer(), n: gl.createBuffer(), c: gl.createBuffer(), l: gl.createBuffer() };
      this.op = { estado: 0, campo: null, mapa: "auto", auto: true, min: 0, max: 1, deformada: false, escala: 50,
        espelhar: true, arestas: true, cortes: { x: null, y: null, z: null }, ocultos: new Set([0]) };
      this.cam = { az: -0.65, el: 0.42, dist: 1, alvo: [0, 0, 0], pan: [0, 0] };
      this.dados = null;
      this.nTri = 0; this.nLin = 0;
      this._eventos();
      new ResizeObserver(() => this.desenhar()).observe(canvas);
    }

    carregar(d) {
      this.dados = d;
      const n = d.nos.length / 3;
      this.X = Float32Array.from(d.nos);
      this.E = Int32Array.from(d.elementos);
      this.npe = d.nos_por_elemento;
      this.facesEl = d.faces || FACES_W6;      // faces do elemento (nós locais), W6 ou H8
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
      this.nNos = n;
      // campos: índices e conversão para Float32 (null -> NaN)
      this.campos = d.campos;
      this.estados = d.estados.map(e => {
        const v = {};
        for (const k in e.valores) v[k] = Float32Array.from(e.valores[k], x => x === null ? NaN : x);
        return { rotulo: e.rotulo, v, u: e.u ? Float32Array.from(e.u, x => x === null ? 0 : x) : null };
      });
      this.op.estado = Math.max(0, this.estados.length - 1);
      this.op.campo = this.campos.length ? this.campos[0].id : null;
      this.op.ocultos = new Set([0]);
      if (d.tipo === "3d" && this.estados.some(e => e.u)) this.op.ocultos.add(4);   // aço fora do modelo mecânico
      this.vista(d.tipo === "2d" ? "topo" : "iso");
      this.atualizarGeometria();
    }

    valoresCampo(id, iest) {
      if (!this.estados.length) return null;
      const e = this.estados[iest];
      if (e && e.v[id]) return e.v[id];
      for (const o of this.estados) if (o.v[id]) return o.v[id];
      return null;
    }

    infoCampo(id) { return this.campos.find(c => c.id === id) || { id, nome: id, unidade: "", local: "no", mapa: "viridis" }; }

    _espelhado() { return !!(this.dados && this.dados.simetria && this.op.espelhar && this.dados.tipo === "3d"); }

    _limites() {
      const lo = this.lo.slice(), hi = this.hi.slice();
      if (this._espelhado()) hi[1] = 2 * this.dados.largura_modelo;
      return [lo, hi];
    }

    _visiveis(copia) {
      const ne = this.mat.length, vis = new Uint8Array(ne), c = this.op.cortes, [lo, hi] = this._limites();
      const lim = [c.x, c.y, c.z].map((f, a) => f === null ? null : lo[a] + (hi[a] - lo[a]) * f);
      const W2 = 2 * this.dados.largura_modelo;
      for (let e = 0; e < ne; e++) {
        if (this.op.ocultos.has(this.mat[e])) continue;
        let ok = true;
        for (let a = 0; a < 3; a++) {
          if (lim[a] === null) continue;
          const x = (a === 1 && copia) ? W2 - this.C[e * 3 + 1] : this.C[e * 3 + a];
          if (x > lim[a] + 1e-9) { ok = false; break; }
        }
        if (ok) vis[e] = 1;
      }
      return vis;
    }

    _faces(vis) {
      const npe = this.npe, E = this.E, n = this.nNos, faces = [];
      if (npe === 3) {
        for (let e = 0; e < vis.length; e++) if (vis[e]) faces.push([e, E[e * 3], E[e * 3 + 1], E[e * 3 + 2]]);
        return faces;
      }
      const mapa = new Map();
      for (let e = 0; e < vis.length; e++) {
        if (!vis[e]) continue;
        for (const f of this.facesEl) {
          const ids = f.map(k => E[e * npe + k]);
          const s = ids.slice().sort((a, b) => a - b);
          const chave = s[0] + n * (s[1] + n * s[2]);
          const ant = mapa.get(chave);
          if (ant === undefined) { mapa.set(chave, faces.length); faces.push([e, ...ids]); }
          else if (ant >= 0) { faces[ant] = null; mapa.set(chave, -1); }
        }
      }
      return faces.filter(f => f);
    }

    atualizarGeometria() {
      if (!this.dados) return;
      const copias = this._espelhado() ? 2 : 1;
      const grupos = [];
      for (let cp = 0; cp < copias; cp++) grupos.push(this._faces(this._visiveis(cp)));
      let nt = 0, nl = 0;
      for (const g of grupos) for (const f of g) { nt += f.length === 5 ? 2 : 1; nl += f.length - 1; }
      this.vNo = new Int32Array(nt * 3); this.vEl = new Int32Array(nt * 3); this.vSim = new Uint8Array(nt * 3);
      this.lNo = new Int32Array(nl * 2); this.lSim = new Uint8Array(nl * 2);
      let iv = 0, il = 0;
      grupos.forEach((g, cp) => { for (const f of g) {
        const e = f[0], ids = f.slice(1);
        const tris = ids.length === 4 ? [[0, 1, 2], [0, 2, 3]] : [[0, 1, 2]];
        for (const t of tris) for (const k of t) { this.vNo[iv] = ids[k]; this.vEl[iv] = e; this.vSim[iv] = cp; iv++; }
        for (let k = 0; k < ids.length; k++) {
          this.lNo[il] = ids[k]; this.lSim[il++] = cp;
          this.lNo[il] = ids[(k + 1) % ids.length]; this.lSim[il++] = cp;
        }
      } });
      this.nTri = nt; this.nLin = nl;
      this.atualizarPosicoes();
      this.atualizarCores();
    }

    _pos(i, sim, out, o) {
      const X = this.X, d = this.op.deformada ? this.estados[this.op.estado] : null;
      let x = X[i * 3], y = X[i * 3 + 1], z = X[i * 3 + 2];
      if (d && d.u) { const s = this.op.escala; x += s * d.u[i * 3]; y += s * d.u[i * 3 + 1]; z += s * d.u[i * 3 + 2]; }
      if (sim) y = 2 * this.dados.largura_modelo - y;
      out[o] = x; out[o + 1] = y; out[o + 2] = z;
    }

    atualizarPosicoes() {
      if (!this.dados) return;
      const gl = this.gl, nv = this.nTri * 3;
      const P = new Float32Array(nv * 3), N = new Float32Array(nv * 3);
      for (let v = 0; v < nv; v++) this._pos(this.vNo[v], this.vSim[v], P, v * 3);
      for (let t = 0; t < this.nTri; t++) {
        const o = t * 9;
        const a = [P[o + 3] - P[o], P[o + 4] - P[o + 1], P[o + 5] - P[o + 2]];
        const b = [P[o + 6] - P[o], P[o + 7] - P[o + 1], P[o + 8] - P[o + 2]];
        const nn = norm(cruz(a, b));
        for (let k = 0; k < 3; k++) { N[o + k * 3] = nn[0]; N[o + k * 3 + 1] = nn[1]; N[o + k * 3 + 2] = nn[2]; }
      }
      this.P = P;
      gl.bindBuffer(gl.ARRAY_BUFFER, this.buf.p); gl.bufferData(gl.ARRAY_BUFFER, P, gl.STATIC_DRAW);
      gl.bindBuffer(gl.ARRAY_BUFFER, this.buf.n); gl.bufferData(gl.ARRAY_BUFFER, N, gl.STATIC_DRAW);
      const L = new Float32Array(this.nLin * 2 * 3);
      for (let k = 0; k < this.nLin * 2; k++) this._pos(this.lNo[k], this.lSim[k], L, k * 3);
      gl.bindBuffer(gl.ARRAY_BUFFER, this.buf.l); gl.bufferData(gl.ARRAY_BUFFER, L, gl.STATIC_DRAW);
      this.desenhar();
    }

    faixa() {
      const id = this.op.campo;
      if (!id || id === "material") return [0, 1];
      const v = this.valoresCampo(id, this.op.estado);
      if (!v) return [0, 1];
      const info = this.infoCampo(id), porEl = info.local === "elemento";
      // faixa só do que está visível
      let lo = Infinity, hi = -Infinity;
      const nv = this.nTri * 3;
      for (let k = 0; k < nv; k++) {
        const x = v[porEl ? this.vEl[k] : this.vNo[k]];
        if (x === x) { if (x < lo) lo = x; if (x > hi) hi = x; }
      }
      if (!isFinite(lo)) return [0, 1];
      if (this.mapaEfetivo() === "divergente" && lo < 0 && hi > 0) { const m = Math.max(-lo, hi); lo = -m; hi = m; }
      if (info.mapa === "dano") { lo = 0; hi = Math.max(hi, 1e-6) > 0.5 ? 1 : hi; }
      if (hi - lo < 1e-12) hi = lo + 1e-12;
      return [lo, hi];
    }

    mapaEfetivo() {
      if (this.op.mapa !== "auto") return this.op.mapa;
      const m = this.infoCampo(this.op.campo).mapa;
      return MAPAS[m] ? m : "viridis";
    }

    atualizarCores() {
      if (!this.dados) return;
      const gl = this.gl, nv = this.nTri * 3, Cc = new Float32Array(nv * 3), id = this.op.campo;
      if (this.op.auto) { const [a, b] = this.faixa(); this.op.min = a; this.op.max = b; }
      const lo = this.op.min, hi = this.op.max, mapa = this.mapaEfetivo();
      const v = id && id !== "material" ? this.valoresCampo(id, this.op.estado) : null;
      const porEl = this.infoCampo(id).local === "elemento";
      for (let k = 0; k < nv; k++) {
        let c;
        if (!v) c = hex(COR_MATERIAL[this.mat[this.vEl[k]]] || "#bbbbbb");
        else {
          const x = v[porEl ? this.vEl[k] : this.vNo[k]];
          c = x === x ? corMapa(mapa, (x - lo) / (hi - lo)) : SEM_DADO;
        }
        Cc[k * 3] = c[0]; Cc[k * 3 + 1] = c[1]; Cc[k * 3 + 2] = c[2];
      }
      gl.bindBuffer(gl.ARRAY_BUFFER, this.buf.c); gl.bufferData(gl.ARRAY_BUFFER, Cc, gl.STATIC_DRAW);
      this.desenhar();
      if (this.aoMudarEscala) this.aoMudarEscala();
    }

    gradienteCSS() {
      const s = MAPAS[this.mapaEfetivo()] || MAPAS.viridis;
      return "linear-gradient(90deg," + s.join(",") + ")";
    }

    vista(nome) {
      const c = this.cam;
      if (nome === "iso") { c.az = -0.65; c.el = 0.42; }
      else if (nome === "frente") { c.az = -Math.PI / 2; c.el = 0.0; }
      else if (nome === "lado") { c.az = 0; c.el = 0.0; }
      else if (nome === "topo") { c.az = -Math.PI / 2; c.el = Math.PI / 2 - 1e-3; }
      c.pan = [0, 0];
      if (this.lo && this.dados) {
        const [lo, hi] = this._limites();
        c.alvo = [(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, (lo[2] + hi[2]) / 2];
        c.dist = 2.6 * Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]);
      }
      this.desenhar();
    }

    _mvp() {
      const c = this.cam, cv = this.canvas;
      const asp = cv.width / Math.max(cv.height, 1);
      const dir = [Math.cos(c.el) * Math.cos(c.az), Math.cos(c.el) * Math.sin(c.az), Math.sin(c.el)];
      const cima = Math.abs(c.el) > 1.5 ? [0, 1, 0] : [0, 0, 1];
      // deslocamento lateral (pan) no plano da tela
      const dirn = norm(dir), dx = norm(cruz(cima, dirn)), dy = cruz(dirn, dx);
      const alvo = [c.alvo[0] + dx[0] * c.pan[0] + dy[0] * c.pan[1], c.alvo[1] + dx[1] * c.pan[0] + dy[1] * c.pan[1],
        c.alvo[2] + dx[2] * c.pan[0] + dy[2] * c.pan[1]];
      const olho = [alvo[0] + dir[0] * c.dist, alvo[1] + dir[1] * c.dist, alvo[2] + dir[2] * c.dist];
      const V = olharPara(olho, alvo, cima);
      const Pm = perspectiva(30 * Math.PI / 180, asp, c.dist / 100, c.dist * 10);
      return { mvp: mult(Pm, V), luz: dirn };
    }

    desenhar() {
      const gl = this.gl, cv = this.canvas, dpr = window.devicePixelRatio || 1;
      const w = Math.max(1, Math.round(cv.clientWidth * dpr)), h = Math.max(1, Math.round(cv.clientHeight * dpr));
      if (cv.width !== w || cv.height !== h) { cv.width = w; cv.height = h; }
      gl.viewport(0, 0, w, h);
      gl.clearColor(0, 0, 0, 0); gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
      if (!this.dados || !this.nTri) return;
      gl.enable(gl.DEPTH_TEST);
      const { mvp, luz } = this._mvp();
      this.ultimoMVP = mvp;
      gl.useProgram(this.pf);
      gl.enable(gl.POLYGON_OFFSET_FILL); gl.polygonOffset(1, 1);
      for (const [nome, b] of [["p", this.buf.p], ["n", this.buf.n], ["c", this.buf.c]]) {
        const loc = gl.getAttribLocation(this.pf, nome);
        gl.bindBuffer(gl.ARRAY_BUFFER, b); gl.enableVertexAttribArray(loc); gl.vertexAttribPointer(loc, 3, gl.FLOAT, false, 0, 0);
      }
      gl.uniformMatrix4fv(gl.getUniformLocation(this.pf, "mvp"), false, new Float32Array(mvp));
      gl.uniform3fv(gl.getUniformLocation(this.pf, "luz"), new Float32Array(luz));
      gl.drawArrays(gl.TRIANGLES, 0, this.nTri * 3);
      gl.disable(gl.POLYGON_OFFSET_FILL);
      for (const nome of ["n", "c"]) gl.disableVertexAttribArray(gl.getAttribLocation(this.pf, nome));
      if (this.op.arestas && this.nLin) {
        gl.useProgram(this.pl);
        const loc = gl.getAttribLocation(this.pl, "p");
        gl.bindBuffer(gl.ARRAY_BUFFER, this.buf.l); gl.enableVertexAttribArray(loc); gl.vertexAttribPointer(loc, 3, gl.FLOAT, false, 0, 0);
        gl.uniformMatrix4fv(gl.getUniformLocation(this.pl, "mvp"), false, new Float32Array(mvp));
        gl.uniform4fv(gl.getUniformLocation(this.pl, "cor"), new Float32Array([0.1, 0.1, 0.1, 0.35]));
        gl.enable(gl.BLEND); gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
        gl.drawArrays(gl.LINES, 0, this.nLin * 2);
        gl.disable(gl.BLEND);
      }
    }

    sondar(px, py) {
      if (!this.P || !this.ultimoMVP) return null;
      const m = this.ultimoMVP, cv = this.canvas, dpr = window.devicePixelRatio || 1;
      const X = px * dpr, Y = py * dpr, P = this.P, nv = this.nTri * 3;
      let melhor = -1, dmin = Infinity, frente = -1, zmin = Infinity;
      for (let k = 0; k < nv; k++) {
        const x = P[k * 3], y = P[k * 3 + 1], z = P[k * 3 + 2];
        const cw = m[3] * x + m[7] * y + m[11] * z + m[15];
        const sx = ((m[0] * x + m[4] * y + m[8] * z + m[12]) / cw * 0.5 + 0.5) * cv.width;
        const sy = (1 - ((m[1] * x + m[5] * y + m[9] * z + m[13]) / cw * 0.5 + 0.5)) * cv.height;
        const sz = (m[2] * x + m[6] * y + m[10] * z + m[14]) / cw;
        const d = Math.hypot(sx - X, sy - Y);
        if (d < dmin) { dmin = d; melhor = k; }
        if (d < 6 * dpr && sz < zmin) { zmin = sz; frente = k; }
      }
      if (frente >= 0) melhor = frente;
      else if (dmin > 14 * dpr) melhor = -1;
      if (melhor < 0) return null;
      const i = this.vNo[melhor], e = this.vEl[melhor], id = this.op.campo;
      const info = this.infoCampo(id);
      const v = id && id !== "material" ? this.valoresCampo(id, this.op.estado) : null;
      const val = v ? v[info.local === "elemento" ? e : i] : null;
      const nm = this.dados.nomes_materiais[this.mat[e]] || this.mat[e];
      return { no: i, elemento: e, material: nm, x: this.X[i * 3], y: this.X[i * 3 + 1], z: this.X[i * 3 + 2],
        valor: val, campo: info.nome, unidade: info.unidade };
    }

    png() { this.desenhar(); return this.canvas.toDataURL("image/png"); }

    _eventos() {
      const cv = this.canvas; let arr = null;
      cv.addEventListener("contextmenu", e => e.preventDefault());
      cv.addEventListener("mousedown", e => { arr = { x: e.clientX, y: e.clientY, pan: e.shiftKey || e.button === 2 }; });
      window.addEventListener("mouseup", () => { arr = null; });
      window.addEventListener("mousemove", e => {
        if (!arr) return;
        const dx = e.clientX - arr.x, dy = e.clientY - arr.y; arr.x = e.clientX; arr.y = e.clientY;
        if (arr.pan) { const s = this.cam.dist / 900; this.cam.pan[0] -= dx * s; this.cam.pan[1] += dy * s; }
        else { this.cam.az -= dx * 0.008; this.cam.el = Math.max(-1.55, Math.min(1.55, this.cam.el + dy * 0.008)); }
        this.desenhar();
      });
      cv.addEventListener("wheel", e => { e.preventDefault(); this.cam.dist *= Math.exp(e.deltaY * 0.0012); this.desenhar(); }, { passive: false });
      let toque = null;
      cv.addEventListener("touchstart", e => { const t = e.touches[0]; toque = { x: t.clientX, y: t.clientY }; }, { passive: true });
      cv.addEventListener("touchmove", e => {
        if (!toque) return; const t = e.touches[0];
        this.cam.az -= (t.clientX - toque.x) * 0.008; this.cam.el = Math.max(-1.55, Math.min(1.55, this.cam.el + (t.clientY - toque.y) * 0.008));
        toque = { x: t.clientX, y: t.clientY }; this.desenhar();
      }, { passive: true });
    }
  }

  window.Visualizador = Visualizador;
  window.MAPAS_CORES = MAPAS;
  window.COR_MATERIAL = COR_MATERIAL;
})();
