/* Interface of the MEF_Masonry program: forms generated from the server schema,
   run with progress, results, history and post-processor. */
(function () {
  "use strict";
  const $ = s => document.querySelector(s);
  const $$ = s => Array.from(document.querySelectorAll(s));
  const FORM_PAGES = ["geometry", "ceramic", "model", "test", "mesh", "sensitivity"];
  const WALLS = ["wall", "sensitivity"];
  // the Sensitivity page only appears with that analysis
  const withSens = () => project && project.model.analysis === "sensitivity";
  const order = () => [...FORM_PAGES.filter(p => p !== "sensitivity" || withSens()), "run", "results", "post", "history"];
  function updateMenu() { const b = $('#menu button[data-page="sensitivity"]'); if (b) b.hidden = !withSens(); }
  let SCHEMA = null, DEFAULT = null, PRESETS = {}, LIB = {}, project = null, page = "geometry";
  const libNotes = {};      // group -> note of the last library item applied
  let task = null, since = 0, timer = null, currentResult = null, currentFolder = null, viewer = null;

  const fmt = v => {
    if (typeof v !== "number") return v;
    if (v !== 0 && (Math.abs(v) < 1e-3 || Math.abs(v) >= 1e6)) return v.toExponential(3);
    return v.toLocaleString("en-GB", { maximumFractionDigits: 6, useGrouping: false });
  };
  const fmtRes = v => typeof v === "number" ? v.toLocaleString("en-GB", { maximumSignificantDigits: 4 }) : v;
  // accepts a decimal point or a decimal comma
  const num = s => { const x = parseFloat(String(s).trim().replace(/\s/g, "").replace(",", ".")); return isFinite(x) ? x : NaN; };
  const esc = s => String(s).replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));

  async function api(route, body) {
    const r = await fetch(route, body === undefined ? {} : { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) });
    const j = await r.json().catch(() => ({}));
    if (!r.ok && !j.errors) throw new Error(j.error || ("HTTP " + r.status));
    return j;
  }

  // ---------------------------------------------------------------- forms
  function visible(field) {
    if (!field.show_if) return true;
    return Object.entries(field.show_if).every(([k, vals]) => {
      const [g, c] = k.split("."); return vals.includes(project[g][c]);
    });
  }

  function buildForm() {
    const g = SCHEMA.find(x => x.id === page);
    $("#form-title").textContent = g.title;
    $("#form-description").textContent = (g.description || "") + (libNotes[g.id] ? " " + libNotes[g.id] : "");
    const f = $("#form"); f.innerHTML = "";
    for (const c of g.fields) {
      if (!visible(c)) continue;
      const div = document.createElement("div");
      div.className = "field " + c.kind + (c.line_break ? " line-break" : "");
      const val = project[g.id][c.id];
      const idc = `c-${g.id}-${c.id}`;
      if (c.kind === "number") {
        div.innerHTML = `<label for="${idc}">${esc(c.label)}</label><div class="row"><input id="${idc}" type="text" inputmode="decimal" value="${esc(fmt(val))}"><span class="unit">${esc(c.unit || "")}</span></div>` +
          (c.help ? `<div class="help">${esc(c.help)}</div>` : "");
        const inp = div.querySelector("input");
        inp.addEventListener("change", () => {
          const x = num(inp.value);
          const bad = isNaN(x) || (c.min !== undefined && x < c.min) || (c.max !== undefined && x > c.max);
          div.classList.toggle("invalid", bad);
          if (!bad) { project[g.id][c.id] = (c.step === 1) ? Math.round(x) : x; onChange(g.id, c.id); }
        });
      } else if (c.kind === "option") {
        const ops = c.options.map(([v, r]) => {
          const off = / under development/.test(r);
          return `<option value="${esc(v)}" ${v === val ? "selected" : ""} ${off ? "disabled class='option-disabled'" : ""}>${esc(r)}</option>`;
        }).join("");
        div.innerHTML = `<label for="${idc}">${esc(c.label)}</label><div class="row"><select id="${idc}">${ops}</select></div>` +
          (c.help ? `<div class="help">${esc(c.help)}</div>` : "");
        div.querySelector("select").addEventListener("change", e => { project[g.id][c.id] = e.target.value; onChange(g.id, c.id); });
      } else if (c.kind === "boolean") {
        div.innerHTML = `<input id="${idc}" type="checkbox" ${val ? "checked" : ""}><div><label for="${idc}">${esc(c.label)}</label>` +
          (c.help ? `<div class="help">${esc(c.help)}</div>` : "") + `</div>`;
        div.querySelector("input").addEventListener("change", e => { project[g.id][c.id] = e.target.checked; onChange(g.id, c.id); });
      } else if (c.kind === "table") {
        div.innerHTML = `<label>${esc(c.label)}</label>` + (c.help ? `<div class="help">${esc(c.help)}</div>` : "") +
          `<table class="mini-table"><thead><tr>${c.columns.map(x => `<th>${esc(x)}</th>`).join("")}<th></th></tr></thead><tbody></tbody></table>` +
          `<div><button type="button" class="secondary">+ row</button> <span class="help">empty = constant value</span></div>`;
        const body = div.querySelector("tbody");
        const draw = () => {
          body.innerHTML = "";
          (project[g.id][c.id] || []).forEach((line, i) => {
            const tr = document.createElement("tr");
            tr.innerHTML = line.map((v, j) => `<td><input type="text" inputmode="decimal" value="${esc(fmt(v))}" data-j="${j}"></td>`).join("") +
              `<td><button type="button" class="secondary" title="remove">×</button></td>`;
            tr.querySelectorAll("input").forEach(inp => inp.addEventListener("change", () => {
              const x = num(inp.value); if (!isNaN(x)) project[g.id][c.id][i][+inp.dataset.j] = x;
            }));
            tr.querySelector("button").addEventListener("click", () => { project[g.id][c.id].splice(i, 1); draw(); });
            body.appendChild(tr);
          });
        };
        div.querySelector("div > button").addEventListener("click", () => {
          const t = project[g.id][c.id] = project[g.id][c.id] || [];
          const last = t[t.length - 1];
          t.push(last ? [last[0] + 100, last[1]] : [20, 1]); draw();
        });
        draw();
      }
      f.appendChild(div);
    }
    $("#mesh-preview").hidden = !(page === "geometry" || page === "mesh");
    $("#bt-previous").disabled = page === FORM_PAGES[0];
  }

  function applyLibrary(group, field) {
    const tab = LIB[`${group}.${field}`];
    const item = tab && tab[project[group][field]];
    if (!item) return;
    for (const [g, d] of Object.entries(item.changes || {})) Object.assign(project[g], JSON.parse(JSON.stringify(d)));
    const note = [item.note, item.source ? "Source: " + item.source + "." : ""].filter(Boolean).join(" ");
    for (const g of new Set([group, ...Object.keys(item.changes || {})])) libNotes[g] = note ? "Library: " + note : "";
  }

  const UNIT_FIELDS = ["unit_length", "width", "course_height", "hole_diameter", "holes_per_unit", "joint_type",
    "joint_thickness", "bedding", "filled_head_joints", "hole_shape", "face_shell", "web",
    "end_shell", "wall_cells", "min_shell", "holes_width", "holes_height", "bed_shell"];
  const MATERIAL_FIELDS = ["k", "rho", "c", "E_unit_gross", "nu", "alpha"];

  function onChange(group, field) {
    // a value changed by hand is no longer the library item
    if (group === "geometry" && UNIT_FIELDS.includes(field)) { project.geometry.unit_library = "custom"; libNotes.geometry = ""; }
    if (group === "ceramic" && MATERIAL_FIELDS.includes(field)) { project.ceramic.material_library = "custom"; libNotes.ceramic = ""; }
    if (group === "ceramic" && field.startsWith("mor_")) project.ceramic.mortar_library = "custom";
    if (group === "ceramic" && field.startsWith("grt_") && field !== "grt_en1992") project.ceramic.grout_library = "custom";
    applyLibrary(group, field);
    const mo = project.model;
    if (group === "model" && field === "analysis") {
      if (mo.analysis === "section_2d") mo.element = "T3";
      else if (mo.element === "T3") mo.element = "W6";
      if (WALLS.includes(mo.analysis)) {
        if (mo.analysis === "sensitivity" || !["ambient", "iso834"].includes(project.test.kind)) project.test.kind = "iso834";
        if (project.mesh.h > 8) project.mesh.h = 5;
        if (project.test.kind === "iso834" && project.test.h_convection < 20) project.test.h_convection = 25;
      }
    }
    if (group === "test" && field === "kind") {
      const ts = project.test;
      if (ts.kind === "iso834" && ts.h_convection < 20) ts.h_convection = 25;
      if (ts.kind !== "iso834" && ts.h_convection === 25) ts.h_convection = 10;
    }
    updateMenu();
    buildForm();
    saveLocal();
    previewPending = true;
    if (["geometry", "mesh"].includes(page)) schedulePreview();
  }

  let previewPending = true, tPreview = null;
  function schedulePreview() {
    clearTimeout(tPreview);
    tPreview = setTimeout(() => { previewPending = false; meshPreview(); }, 600);
  }

  // ---------------------------------------------------------------- navigation
  function goTo(p) {
    page = p;
    $$("#menu button").forEach(b => b.classList.toggle("active", b.dataset.page === p));
    const isForm = FORM_PAGES.includes(p);
    $("#pg-form").hidden = !isForm;
    $("#pg-run").hidden = p !== "run";
    $("#pg-results").hidden = p !== "results";
    $("#pg-post").hidden = p !== "post";
    $("#pg-history").hidden = p !== "history";
    if (isForm) buildForm();
    if ((p === "geometry" || p === "mesh") && previewPending) schedulePreview();
    if (p === "run") projectSummary();
    if (p === "history") loadHistory();
    if (p === "post" && viewer) setTimeout(() => viewer.draw(), 30);
    window.scrollTo(0, 0);
  }

  // ---------------------------------------------------------------- project
  function saveLocal() { try { localStorage.setItem("mef_masonry_project", JSON.stringify(project)); } catch (e) { /* ignore */ } }
  function setProject(p) {
    project = JSON.parse(JSON.stringify(DEFAULT));
    for (const g in p) { if (typeof p[g] === "object" && p[g] && project[g]) Object.assign(project[g], p[g]); else project[g] = p[g]; }
    $("#project-name").value = project.name || "Project";
    updateMenu();
    saveLocal();
    previewPending = true;
    if (["geometry", "mesh"].includes(page)) schedulePreview();
    if (FORM_PAGES.includes(page)) buildForm();
    if (page === "run") projectSummary();
  }
  // project read from a file, from the history or from the browser: the server converts
  // the projects of version 1.x (Portuguese names) to the current format
  async function setProjectFile(p) {
    const r = await api("/api/complete", p);
    setProject(r.project || p);
  }

  async function meshPreview() {
    $("#preview-info").textContent = "Generating the mesh…";
    try {
      const r = await api("/api/mesh", project);
      if (r.errors && r.errors.length) { $("#preview-info").textContent = r.errors.join(" "); $("#preview-img").removeAttribute("src"); return; }
      $("#preview-img").src = "data:image/png;base64," + r.png;
      $("#preview-info").textContent = r.summary;
    } catch (e) { $("#preview-info").textContent = "Error: " + e.message; }
  }

  function optionLabel(group, field) {
    const c = SCHEMA.find(g => g.id === group).fields.find(x => x.id === field);
    const op = c.options.find(o => o[0] === project[group][field]);
    return op ? op[1] : project[group][field];
  }

  async function projectSummary() {
    const g = project.geometry, ts = project.test, me = project.mesh;
    const items = [
      ["Geometry", `${fmt(g.units_per_course)} unit(s) × ${g.n_courses} courses; unit ${fmt(g.unit_length)} × ${fmt(g.width)} × ${fmt(g.course_height)} mm (centre to centre); ${g.hole_diameter > 0 ? `${g.holes_per_unit} hole(s) of ${fmt(g.hole_diameter)} mm` : "solid"}; ${optionLabel("geometry", "bond")}; ${g.joint_type === "mortar" ? `mortar of ${fmt(g.joint_thickness)} mm (${g.bedding === "lateral" ? "on the face shells only" : "over the whole face"}${g.filled_head_joints ? "" : ", empty head joints"})` : `dry joint of ${fmt(g.joint_thickness)} mm`}`],
      ["Analysis", optionLabel("model", "analysis")],
      WALLS.includes(project.model.analysis)
        ? ["Wall", `${fmt(g.units_per_course * g.unit_length)} × ${fmt(g.n_courses * g.course_height)} × ${fmt(g.width)} mm (H/t = ${fmt(Math.round(g.n_courses * g.course_height / g.width * 10) / 10)}); ${ts.base_support} at the base and ${ts.top_support} at the top; eccentricities ${fmt(ts.top_eccentricity)} and ${fmt(ts.base_eccentricity)} mm; imperfection ${ts.imperfection > 0 ? "H/" + fmt(ts.imperfection) : "none"}`]
        : ["Element", optionLabel("model", "element")],
      ["Test", optionLabel("test", "kind") + (["hot", "residual", "transient"].includes(ts.kind) ? ` — ${fmt(ts.T_target)} °C, ${fmt(ts.rate)} °C/min, dwell ${fmt(ts.dwell)} min` : "") + (ts.kind === "transient" ? `, load ${fmt(ts.load)} kN` : "") + (ts.kind === "iso834" ? ` — ${fmt(ts.iso_duration)} min` : "") + (WALLS.includes(project.model.analysis) && ts.kind === "iso834" ? `, load ${fmt(ts.load)} kN` : "")],
      ["Mesh", withSens() && project.sensitivity.fast_mesh ? `fast mesh of the screening: section with h = ${fmt(Math.max(me.h, 10))} mm, ${Math.min(me.n_elements, 16)} beam elements, Δt = ${fmt(me.dt)} s`
        : WALLS.includes(project.model.analysis) ? `section with h = ${fmt(me.h)} mm, ${me.n_elements} beam elements, Δt = ${fmt(me.dt)} s`
        : project.model.analysis === "section_2d" ? `h = ${fmt(me.h)} mm, Δt = ${fmt(me.dt)} s` : `h = ${fmt(me.h)} mm, dz = ${fmt(me.dz)} mm, Δt = ${fmt(me.dt)} s`],
    ];
    if (withSens()) {
      const se = project.sensitivity;
      const k = Object.keys(se).filter(c => c.startsWith("var_") && se[c]).length;
      const n = se.trajectories * (k + 1) + 1;
      items.push(["Sensitivity", `${k} parameters, ${se.trajectories} trajectories: ${n} wall runs` +
        ` (about ${se.fast_mesh ? "10 to 30 s" : "1 to 3 min"} each per process; ${se.processes > 0 ? se.processes + " process(es)" : "processes: automatic"})`]);
    }
    if (project.model.analysis === "nonlinear") items.push(["Materials", `unit: ft = ${fmt(project.model.ft)} MPa, fc = ${fmt(project.model.fc)} MPa` + (g.joint_type === "mortar" ? `; mortar: fc = ${fmt(project.ceramic.mor_fc)} MPa, E = ${fmt(project.ceramic.mor_E)} MPa` : `; gap ${fmt(project.model.gap_std)} mm`)]);
    $("#project-summary").innerHTML = items.map(([a, b]) => `<dt>${esc(a)}</dt><dd>${esc(b)}</dd>`).join("");
    const r = await api("/api/validate", project).catch(e => ({ errors: [e.message] }));
    $("#validation-errors").innerHTML = (r.errors || []).map(e => `<li>${esc(e)}</li>`).join("");
    $("#bt-run").disabled = !!(r.errors && r.errors.length) || (task && task.state === "running");
  }

  // ---------------------------------------------------------------- run
  async function run() {
    project.name = $("#project-name").value || "Project";
    const r = await api("/api/run", project).catch(e => ({ errors: [e.message] }));
    if (r.errors) { $("#validation-errors").innerHTML = r.errors.map(e => `<li>${esc(e)}</li>`).join(""); return; }
    task = { id: r.id, state: "running", folder: r.folder }; since = 0;
    $("#log").textContent = "";
    $("#bt-run").disabled = true; $("#bt-cancel").hidden = false;
    follow();
  }

  function follow() {
    clearInterval(timer);
    timer = setInterval(async () => {
      let r;
      try { r = await api(`/api/task/${task.id}?since=${since}`); } catch (e) { return; }
      since = r.n;
      const log = $("#log");
      const atEnd = log.scrollTop + log.clientHeight >= log.scrollHeight - 20;
      if (r.log.length) { log.textContent += r.log.join("\n") + "\n"; if (atEnd) log.scrollTop = log.scrollHeight; }
      const min = Math.floor(r.time_s / 60), sec = Math.round(r.time_s % 60);
      $("#timer").textContent = `${min} min ${String(sec).padStart(2, "0")} s`;
      const s = $("#status");
      s.className = "status " + (r.state === "done" ? "done" : r.state === "error" ? "error" : "running");
      s.textContent = { "running": "Running…", "queued": "Queued…", "done": "Done.", "error": "Error: " + (r.error || ""), "cancelled": "Cancelled." }[r.state] || r.state;
      $("#menu-status").textContent = r.state === "running" ? `Analysis running (${min} min)` : "";
      task.state = r.state;
      if (r.state !== "running" && r.state !== "queued") {
        clearInterval(timer);
        $("#bt-cancel").hidden = true; $("#bt-run").disabled = false;
        if (r.state === "done") { await openResult(r.folder); goTo("results"); }
      }
    }, 1000);
  }

  // ---------------------------------------------------------------- results
  async function openResult(folder) {
    const r = await api(`/api/result/${encodeURIComponent(folder)}`);
    if (r.error) { alert(r.error); return; }
    currentResult = r; currentFolder = folder;
    showResult();
    if (r.visual) loadVisual(folder);
  }

  function showResult() {
    const r = currentResult, p = currentFolder, base = `/files/${encodeURIComponent(p)}/`;
    $("#res-source").innerHTML = `${esc(r.project && r.project.name || "")} · folder <code>${esc(p)}</code>`;
    let h = `<div class="result-summary">` + r.summary.map(([lbl, v, u]) =>
      `<div class="indicator"><div class="r">${esc(lbl)}</div><div class="v">${esc(fmtRes(v))}<small>${esc(u || "")}</small></div></div>`).join("") + `</div>`;
    if (r.notices && r.notices.length) h += `<ul class="notices">${r.notices.map(a => `<li>${esc(a)}</li>`).join("")}</ul>`;
    if (r.visual) h += `<p><button id="bt-go-post">◈ Open in the post-processor</button></p>`;
    if (r.figures.length) h += `<div class="figures">${r.figures.map(f => `<figure><img src="${base}${encodeURIComponent(f.file)}" alt="${esc(f.title)}"><figcaption>${esc(f.title)}</figcaption></figure>`).join("")}</div>`;
    h += `<div class="card" style="margin-top:16px"><h2>Tables</h2>` + r.tables.map(t =>
      `<details class="result-table"><summary>${esc(t.title)} (${t.lines.length} rows)</summary><div class="scroll"><table class="data-table"><thead><tr>${t.columns.map(c => `<th>${esc(c)}</th>`).join("")}</tr></thead><tbody>${t.lines.slice(0, 400).map(l => `<tr>${l.map(v => `<td>${esc(fmtRes(v))}</td>`).join("")}</tr>`).join("")}</tbody></table></div></details>`).join("") + `</div>`;
    h += `<div class="card"><h2>Files</h2><div class="files">${r.files.map(a => `<a href="${base}${encodeURIComponent(a.name)}?download=1" title="${esc(a.description)}">${esc(a.name)}</a>`).join("")}<a href="${base}project.json?download=1">project.json</a><a href="${base}results.json?download=1">results.json</a></div></div>`;
    $("#res-content").innerHTML = h;
    const b = $("#bt-go-post"); if (b) b.addEventListener("click", () => goTo("post"));
  }

  async function loadHistory() {
    const r = await api("/api/history");
    $("#results-folder").textContent = r.folder;
    const names = { thermal: "thermal 3D", linear: "linear 3D", nonlinear: "non-linear 3D", section_2d: "2D section", wall: "wall", sensitivity: "sensitivity (wall)" };
    const t = $("#history-table");
    if (!r.runs.length) { t.innerHTML = "<tr><td>No analysis stored yet.</td></tr>"; return; }
    t.innerHTML = `<thead><tr><th>Date and project</th><th>Analysis</th><th>Test</th><th>Status</th><th></th></tr></thead><tbody>` +
      r.runs.map(e => `<tr><td>${esc(e.folder.slice(0, 8))} ${esc(e.folder.slice(9, 11))}:${esc(e.folder.slice(11, 13))} · ${esc(e.name)}</td><td>${esc(names[e.analysis] || e.analysis)}</td><td>${esc(e.test)}</td><td>${e.done ? "done" : "no results"}</td>
        <td>${e.done ? `<button class="secondary" data-action="res" data-folder="${esc(e.folder)}">Results</button> ` : ""}${e.visual ? `<button class="secondary" data-action="post" data-folder="${esc(e.folder)}">Post</button> ` : ""}<button class="secondary" data-action="proj" data-folder="${esc(e.folder)}">Use inputs</button></td></tr>`).join("") + "</tbody>";
    t.querySelectorAll("button").forEach(b => b.addEventListener("click", async () => {
      const folder = b.dataset.folder;
      if (b.dataset.action === "proj") {
        const f = await fetch(`/files/${encodeURIComponent(folder)}/project.json`);
        const p = f.ok ? await f.json() : await (await fetch(`/files/${encodeURIComponent(folder)}/projeto.json`)).json();   // runs of version 1.x
        await setProjectFile(p); goTo("geometry");
      } else { await openResult(folder); goTo(b.dataset.action === "post" ? "post" : "results"); }
    }));
  }

  // ---------------------------------------------------------------- post-processor
  function startViewer() {
    if (viewer) return viewer;
    try { viewer = new Viewer($("#post-canvas")); }
    catch (e) { $("#post-empty").textContent = "WebGL could not be started: " + e.message; return null; }
    viewer.onScaleChange = updateLegend;
    $("#post-canvas").addEventListener("dblclick", ev => {
      const rc = ev.target.getBoundingClientRect();
      const s = viewer.probe(ev.clientX - rc.left, ev.clientY - rc.top);
      const box = $("#post-probe");
      if (!s) { box.hidden = true; return; }
      box.hidden = false;
      box.innerHTML = `<b>${esc(s.field)}:</b> ${s.value === null || s.value !== s.value ? "—" : esc(fmtRes(s.value))} ${esc(s.unit)}<br>` +
        `node ${s.node}, element ${s.element} (${esc(s.material)})<br>x = ${fmtRes(s.x)}, y = ${fmtRes(s.y)}, z = ${fmtRes(s.z)} mm`;
    });
    return viewer;
  }

  async function loadVisual(folder) {
    if (!startViewer()) return;
    $("#post-source").textContent = "Loading…";
    const d = await (await fetch(`/files/${encodeURIComponent(folder)}/visual.json`)).json();
    if (!d.nodes) { $("#post-source").textContent = "This result was written by version 1.x: run it again to see it in the post-processor."; return; }
    viewer.load(d);
    $("#post-empty").hidden = true;
    $("#post-source").textContent = `${d.description} · ${folder}`;
    const se = $("#post-state");
    se.innerHTML = d.states.map((e, i) => `<option value="${i}">${esc(e.label)}</option>`).join("");
    se.value = viewer.op.state;
    const sl = $("#post-state-slider"); sl.max = Math.max(0, d.states.length - 1); sl.value = viewer.op.state;
    $("#post-field").innerHTML = d.fields.map(c => `<option value="${c.id}">${esc(c.name)}${c.unit ? " (" + esc(c.unit) + ")" : ""}</option>`).join("") +
      `<option value="material">Material</option>`;
    $("#post-field").value = viewer.op.field || "material";
    const used = [...new Set(d.material)].sort();
    $("#post-materials").innerHTML = used.map(m => `<label class="check"><input type="checkbox" data-mat="${m}" ${viewer.op.hidden.has(m) ? "" : "checked"}> <span style="display:inline-block;width:10px;height:10px;border-radius:2px;background:${MATERIAL_COLOUR[m] || "#bbb"}"></span> ${esc(d.material_names[m] || m)}</label>`).join("");
    $$("#post-materials input").forEach(c => c.addEventListener("change", () => {
      const m = +c.dataset.mat; if (c.checked) viewer.op.hidden.delete(m); else viewer.op.hidden.add(m); viewer.updateGeometry();
    }));
    $("#post-deformed").disabled = !d.states.some(e => e.u);
    $("#post-mirror").disabled = !(d.symmetry && d.kind === "3d");
    updateLegend();
  }

  function updateLegend() {
    if (!viewer || !viewer.data) return;
    const id = viewer.op.field, info = viewer.fieldInfo(id);
    const leg = $("#post-legend"); leg.hidden = false;
    if (!id || id === "material") {
      $("#legend-title").textContent = "Material";
      $("#legend-bar").style.background = "linear-gradient(90deg," + Object.values(MATERIAL_COLOUR).join(",") + ")";
      $("#legend-min").textContent = ""; $("#legend-mid").textContent = ""; $("#legend-max").textContent = "";
      return;
    }
    $("#legend-title").textContent = info.name + (info.unit ? ` (${info.unit})` : "") + " — " + viewer.states[viewer.op.state].label;
    $("#legend-bar").style.background = viewer.cssGradient();
    $("#legend-min").textContent = fmtRes(viewer.op.min);
    $("#legend-mid").textContent = fmtRes((viewer.op.min + viewer.op.max) / 2);
    $("#legend-max").textContent = fmtRes(viewer.op.max);
    if (viewer.op.auto) { $("#post-min").value = +viewer.op.min.toPrecision(4); $("#post-max").value = +viewer.op.max.toPrecision(4); }
  }

  function bindPostControls() {
    const ok = () => viewer && viewer.data;
    $("#post-state").addEventListener("change", e => { if (!ok()) return; viewer.op.state = +e.target.value; $("#post-state-slider").value = e.target.value; if (viewer.op.deformed) viewer.updatePositions(); viewer.updateColours(); });
    $("#post-state-slider").addEventListener("input", e => { if (!ok()) return; viewer.op.state = +e.target.value; $("#post-state").value = e.target.value; if (viewer.op.deformed) viewer.updatePositions(); viewer.updateColours(); });
    $("#post-field").addEventListener("change", e => { if (!ok()) return; viewer.op.field = e.target.value; viewer.updateColours(); });
    $("#post-cmap").addEventListener("change", e => { if (!ok()) return; viewer.op.cmap = e.target.value; viewer.updateColours(); });
    $("#post-auto").addEventListener("change", e => { if (!ok()) return; viewer.op.auto = e.target.checked; viewer.updateColours(); });
    for (const k of ["min", "max"]) $("#post-" + k).addEventListener("change", e => {
      if (!ok()) return; viewer.op.auto = false; $("#post-auto").checked = false; viewer.op[k] = num(e.target.value); viewer.updateColours();
    });
    $("#post-deformed").addEventListener("change", e => { if (!ok()) return; viewer.op.deformed = e.target.checked; viewer.updatePositions(); });
    $("#post-scale").addEventListener("change", e => { if (!ok()) return; viewer.op.scale = num(e.target.value) || 0; if (viewer.op.deformed) viewer.updatePositions(); });
    $("#post-mirror").addEventListener("change", e => { if (!ok()) return; viewer.op.mirror = e.target.checked; viewer.updateGeometry(); viewer.view("iso"); });
    $("#post-edges").addEventListener("change", e => { if (!ok()) return; viewer.op.edges = e.target.checked; viewer.draw(); });
    for (const a of ["x", "y", "z"]) {
      const f = () => { if (!ok()) return; viewer.op.cuts[a] = $("#cut-" + a).checked ? +$("#cut-" + a + "-v").value / 100 : null; viewer.updateGeometry(); };
      $("#cut-" + a).addEventListener("change", f); $("#cut-" + a + "-v").addEventListener("input", f);
    }
    $$("[data-view]").forEach(b => b.addEventListener("click", () => viewer && viewer.view(b.dataset.view)));
    $("#post-png").addEventListener("click", () => {
      if (!ok()) return;
      const a = document.createElement("a"); a.href = viewer.png(); a.download = `${currentFolder || "mef"}_${viewer.op.field}.png`; a.click();
    });
  }

  // ---------------------------------------------------------------- start
  async function start() {
    const r = await api("/api/schema");
    SCHEMA = r.schema; DEFAULT = r.default; PRESETS = r.presets; LIB = r.library || {};
    $("#version").textContent = "version " + r.version;
    const sel = $("#presets");
    for (const k of Object.keys(PRESETS)) { const o = document.createElement("option"); o.value = k; o.textContent = k; sel.appendChild(o); }
    sel.addEventListener("change", () => { if (sel.value) { setProject({ ...PRESETS[sel.value], name: sel.value }); sel.value = ""; } });
    let saved = null, old = null;
    try {
      saved = JSON.parse(localStorage.getItem("mef_masonry_project"));
      if (!saved) old = JSON.parse(localStorage.getItem("mef_masonry_projeto") || localStorage.getItem("mef_bcepq_projeto"));   // version 1.x
    } catch (e) { /* ignore */ }
    if (old) await setProjectFile(old).catch(() => setProject(DEFAULT));
    else setProject(saved || DEFAULT);
    $$("#menu button").forEach(b => b.addEventListener("click", () => goTo(b.dataset.page)));
    $("#bt-next").addEventListener("click", () => { const o = order(); goTo(o[Math.min(o.length - 1, o.indexOf(page) + 1)]); });
    $("#bt-previous").addEventListener("click", () => { const o = order(); goTo(o[Math.max(0, o.indexOf(page) - 1)]); });
    $("#bt-preview").addEventListener("click", meshPreview);
    $("#bt-run").addEventListener("click", run);
    $("#bt-cancel").addEventListener("click", () => task && api(`/api/cancel/${task.id}`, {}));
    $("#project-name").addEventListener("change", e => { project.name = e.target.value; saveLocal(); });
    $("#bt-new").addEventListener("click", () => setProject(DEFAULT));
    $("#bt-save").addEventListener("click", () => {
      project.name = $("#project-name").value;
      const a = document.createElement("a");
      a.href = URL.createObjectURL(new Blob([JSON.stringify(project, null, 2)], { type: "application/json" }));
      a.download = (project.name || "project").replace(/[^\w\-]+/g, "_") + ".json"; a.click();
    });
    $("#bt-open").addEventListener("click", () => $("#file-open").click());
    $("#file-open").addEventListener("change", async e => {
      const f = e.target.files[0]; if (!f) return;
      try { await setProjectFile(JSON.parse(await f.text())); } catch (err) { alert("Invalid file: " + err.message); }
      e.target.value = "";
    });
    bindPostControls();
    goTo("geometry");
    // resumes an analysis in progress (page reloaded)
    const cur = await api("/api/current").catch(() => ({}));
    if (cur.id && (cur.state === "running" || cur.state === "queued")) {
      task = { id: cur.id, state: cur.state, folder: cur.folder }; since = 0; $("#bt-cancel").hidden = false; follow();
    }
  }
  start().catch(e => { document.body.insertAdjacentHTML("afterbegin", `<p style="color:#c00;padding:12px">Error at start-up: ${esc(e.message)}</p>`); });
})();
