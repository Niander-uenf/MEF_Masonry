/* Interface do programa MEF_Masonry: formulários gerados pelo esquema do servidor,
   execução com acompanhamento, resultados, histórico e pós-processador. */
(function () {
  "use strict";
  const $ = s => document.querySelector(s);
  const $$ = s => Array.from(document.querySelectorAll(s));
  const PAGINAS_FORM = ["geometria", "ceramica", "modelo", "ensaio", "malha", "sensibilidade"];
  const PAREDES = ["parede", "sensibilidade"];
  // a página Sensibilidade só aparece com essa análise
  const comSens = () => projeto && projeto.modelo.analise === "sensibilidade";
  const ordem = () => [...PAGINAS_FORM.filter(p => p !== "sensibilidade" || comSens()), "executar", "resultados", "pos", "historico"];
  function atualizarMenu() { const b = $('#menu button[data-pagina="sensibilidade"]'); if (b) b.hidden = !comSens(); }
  let ESQ = null, PADRAO = null, PRED = {}, BIB = {}, projeto = null, pagina = "geometria";
  const notasBib = {};      // grupo -> nota do último item da biblioteca aplicado
  let tarefa = null, desde = 0, temporizador = null, resultadoAtual = null, pastaAtual = null, vis = null;

  const fmt = v => {
    if (typeof v !== "number") return v;
    if (v !== 0 && (Math.abs(v) < 1e-3 || Math.abs(v) >= 1e6)) return v.toExponential(3).replace(".", ",");
    return v.toLocaleString("pt-BR", { maximumFractionDigits: 6, useGrouping: false });
  };
  const fmtRes = v => typeof v === "number" ? v.toLocaleString("pt-BR", { maximumSignificantDigits: 4 }) : v;
  const num = s => { const x = parseFloat(String(s).trim().replace(/\s/g, "").replace(",", ".")); return isFinite(x) ? x : NaN; };
  const esc = s => String(s).replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));

  async function api(rota, corpo) {
    const r = await fetch(rota, corpo === undefined ? {} : { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(corpo) });
    const j = await r.json().catch(() => ({}));
    if (!r.ok && !j.erros) throw new Error(j.erro || ("HTTP " + r.status));
    return j;
  }

  // ---------------------------------------------------------------- formulários
  function visivel(campo) {
    if (!campo.mostrar_se) return true;
    return Object.entries(campo.mostrar_se).every(([k, vals]) => {
      const [g, c] = k.split("."); return vals.includes(projeto[g][c]);
    });
  }

  function montarFormulario() {
    const g = ESQ.find(x => x.id === pagina);
    $("#form-titulo").textContent = g.titulo;
    $("#form-descricao").textContent = (g.descricao || "") + (notasBib[g.id] ? " " + notasBib[g.id] : "");
    const f = $("#formulario"); f.innerHTML = "";
    for (const c of g.campos) {
      if (!visivel(c)) continue;
      const div = document.createElement("div");
      div.className = "campo " + c.tipo + (c.quebra ? " quebra" : "");
      const val = projeto[g.id][c.id];
      const idc = `c-${g.id}-${c.id}`;
      if (c.tipo === "numero") {
        div.innerHTML = `<label for="${idc}">${esc(c.rotulo)}</label><div class="linha"><input id="${idc}" type="text" inputmode="decimal" value="${esc(fmt(val))}"><span class="unidade">${esc(c.unidade || "")}</span></div>` +
          (c.ajuda ? `<div class="ajuda">${esc(c.ajuda)}</div>` : "");
        const inp = div.querySelector("input");
        inp.addEventListener("change", () => {
          const x = num(inp.value);
          const ruim = isNaN(x) || (c.min !== undefined && x < c.min) || (c.max !== undefined && x > c.max);
          div.classList.toggle("invalido", ruim);
          if (!ruim) { projeto[g.id][c.id] = (c.passo === 1) ? Math.round(x) : x; aoMudar(g.id, c.id); }
        });
      } else if (c.tipo === "opcao") {
        const ops = c.opcoes.map(([v, r]) => {
          const desab = / em desenvolvimento/.test(r);
          return `<option value="${esc(v)}" ${v === val ? "selected" : ""} ${desab ? "disabled class='op-desabilitada'" : ""}>${esc(r)}</option>`;
        }).join("");
        div.innerHTML = `<label for="${idc}">${esc(c.rotulo)}</label><div class="linha"><select id="${idc}">${ops}</select></div>` +
          (c.ajuda ? `<div class="ajuda">${esc(c.ajuda)}</div>` : "");
        div.querySelector("select").addEventListener("change", e => { projeto[g.id][c.id] = e.target.value; aoMudar(g.id, c.id); });
      } else if (c.tipo === "logico") {
        div.innerHTML = `<input id="${idc}" type="checkbox" ${val ? "checked" : ""}><div><label for="${idc}">${esc(c.rotulo)}</label>` +
          (c.ajuda ? `<div class="ajuda">${esc(c.ajuda)}</div>` : "") + `</div>`;
        div.querySelector("input").addEventListener("change", e => { projeto[g.id][c.id] = e.target.checked; aoMudar(g.id, c.id); });
      } else if (c.tipo === "tabela") {
        div.innerHTML = `<label>${esc(c.rotulo)}</label>` + (c.ajuda ? `<div class="ajuda">${esc(c.ajuda)}</div>` : "") +
          `<table class="mini-tabela"><thead><tr>${c.colunas.map(x => `<th>${esc(x)}</th>`).join("")}<th></th></tr></thead><tbody></tbody></table>` +
          `<div><button type="button" class="secundario">+ linha</button> <span class="ajuda">vazia = valor constante</span></div>`;
        const corpo = div.querySelector("tbody");
        const desenhar = () => {
          corpo.innerHTML = "";
          (projeto[g.id][c.id] || []).forEach((linha, i) => {
            const tr = document.createElement("tr");
            tr.innerHTML = linha.map((v, j) => `<td><input type="text" inputmode="decimal" value="${esc(fmt(v))}" data-j="${j}"></td>`).join("") +
              `<td><button type="button" class="secundario" title="remover">×</button></td>`;
            tr.querySelectorAll("input").forEach(inp => inp.addEventListener("change", () => {
              const x = num(inp.value); if (!isNaN(x)) projeto[g.id][c.id][i][+inp.dataset.j] = x;
            }));
            tr.querySelector("button").addEventListener("click", () => { projeto[g.id][c.id].splice(i, 1); desenhar(); });
            corpo.appendChild(tr);
          });
        };
        div.querySelector("div > button").addEventListener("click", () => {
          const t = projeto[g.id][c.id] = projeto[g.id][c.id] || [];
          const ult = t[t.length - 1];
          t.push(ult ? [ult[0] + 100, ult[1]] : [20, 1]); desenhar();
        });
        desenhar();
      }
      f.appendChild(div);
    }
    $("#previa-malha").hidden = !(pagina === "geometria" || pagina === "malha");
    $("#bt-anterior").disabled = pagina === PAGINAS_FORM[0];
  }

  function aplicarBiblioteca(grupo, campo) {
    const tab = BIB[`${grupo}.${campo}`];
    const item = tab && tab[projeto[grupo][campo]];
    if (!item) return;
    for (const [g, d] of Object.entries(item.mudancas || {})) Object.assign(projeto[g], JSON.parse(JSON.stringify(d)));
    const nota = [item.nota, item.fonte ? "Fonte: " + item.fonte + "." : ""].filter(Boolean).join(" ");
    for (const g of new Set([grupo, ...Object.keys(item.mudancas || {})])) notasBib[g] = nota ? "Biblioteca: " + nota : "";
  }

  const CAMPOS_BLOCO = ["comprimento_bloco", "largura", "altura_fiada", "diametro_furo", "furos_por_bloco", "tipo_junta",
    "espessura_junta", "assentamento", "juntas_verticais_preenchidas", "forma_furo", "parede_externa", "septo",
    "parede_topo", "celulas_parede", "parede_minima", "furos_largura", "furos_altura", "parede_leito"];
  const CAMPOS_MATERIAL = ["k", "rho", "c", "E_bloco_bruto", "nu", "alfa"];

  function aoMudar(grupo, campo) {
    // um valor trocado à mão deixa de ser o item da biblioteca
    if (grupo === "geometria" && CAMPOS_BLOCO.includes(campo)) { projeto.geometria.bloco_biblioteca = "personalizado"; notasBib.geometria = ""; }
    if (grupo === "ceramica" && CAMPOS_MATERIAL.includes(campo)) { projeto.ceramica.material_biblioteca = "personalizado"; notasBib.ceramica = ""; }
    if (grupo === "ceramica" && campo.startsWith("arg_")) projeto.ceramica.argamassa_biblioteca = "personalizada";
    if (grupo === "ceramica" && campo.startsWith("gra_") && campo !== "gra_en1992") projeto.ceramica.graute_biblioteca = "personalizado";
    aplicarBiblioteca(grupo, campo);
    const mo = projeto.modelo;
    if (grupo === "modelo" && campo === "analise") {
      if (mo.analise === "secao_2d") mo.elemento = "T3";
      else if (mo.elemento === "T3") mo.elemento = "W6";
      if (PAREDES.includes(mo.analise)) {
        if (mo.analise === "sensibilidade" || !["ambiente", "iso834"].includes(projeto.ensaio.tipo)) projeto.ensaio.tipo = "iso834";
        if (projeto.malha.h > 8) projeto.malha.h = 5;
        if (projeto.ensaio.tipo === "iso834" && projeto.ensaio.h_conveccao < 20) projeto.ensaio.h_conveccao = 25;
      }
    }
    if (grupo === "ensaio" && campo === "tipo") {
      const en = projeto.ensaio;
      if (en.tipo === "iso834" && en.h_conveccao < 20) en.h_conveccao = 25;
      if (en.tipo !== "iso834" && en.h_conveccao === 25) en.h_conveccao = 10;
    }
    atualizarMenu();
    montarFormulario();
    salvarLocal();
    previaPendente = true;
    if (["geometria", "malha"].includes(pagina)) agendarPrevia();
  }

  let previaPendente = true, tPrevia = null;
  function agendarPrevia() {
    clearTimeout(tPrevia);
    tPrevia = setTimeout(() => { previaPendente = false; previaMalha(); }, 600);
  }

  // ---------------------------------------------------------------- navegação
  function irPara(p) {
    pagina = p;
    $$("#menu button").forEach(b => b.classList.toggle("ativo", b.dataset.pagina === p));
    const ehForm = PAGINAS_FORM.includes(p);
    $("#pg-formulario").hidden = !ehForm;
    $("#pg-executar").hidden = p !== "executar";
    $("#pg-resultados").hidden = p !== "resultados";
    $("#pg-pos").hidden = p !== "pos";
    $("#pg-historico").hidden = p !== "historico";
    if (ehForm) montarFormulario();
    if ((p === "geometria" || p === "malha") && previaPendente) agendarPrevia();
    if (p === "executar") resumoProjeto();
    if (p === "historico") carregarHistorico();
    if (p === "pos" && vis) setTimeout(() => vis.desenhar(), 30);
    window.scrollTo(0, 0);
  }

  // ---------------------------------------------------------------- projeto
  function salvarLocal() { try { localStorage.setItem("mef_masonry_projeto", JSON.stringify(projeto)); } catch (e) { /* ignora */ } }
  function definirProjeto(p) {
    projeto = JSON.parse(JSON.stringify(PADRAO));
    for (const g in p) { if (typeof p[g] === "object" && p[g] && projeto[g]) Object.assign(projeto[g], p[g]); else projeto[g] = p[g]; }
    $("#nome-projeto").value = projeto.nome || "Projeto";
    atualizarMenu();
    salvarLocal();
    previaPendente = true;
    if (["geometria", "malha"].includes(pagina)) agendarPrevia();
    if (PAGINAS_FORM.includes(pagina)) montarFormulario();
    if (pagina === "executar") resumoProjeto();
  }

  async function previaMalha() {
    $("#previa-info").textContent = "Gerando a malha…";
    try {
      const r = await api("/api/malha", projeto);
      if (r.erros && r.erros.length) { $("#previa-info").textContent = r.erros.join(" "); $("#previa-img").removeAttribute("src"); return; }
      $("#previa-img").src = "data:image/png;base64," + r.png;
      $("#previa-info").textContent = r.resumo;
    } catch (e) { $("#previa-info").textContent = "Erro: " + e.message; }
  }

  function rotuloOpcao(grupo, campo) {
    const c = ESQ.find(g => g.id === grupo).campos.find(x => x.id === campo);
    const op = c.opcoes.find(o => o[0] === projeto[grupo][campo]);
    return op ? op[1] : projeto[grupo][campo];
  }

  async function resumoProjeto() {
    const g = projeto.geometria, en = projeto.ensaio, ma = projeto.malha;
    const itens = [
      ["Geometria", `${fmt(g.blocos_por_fiada)} bloco(s) × ${g.n_fiadas} fiadas; bloco ${fmt(g.comprimento_bloco)} × ${fmt(g.largura)} × ${fmt(g.altura_fiada)} mm (eixo a eixo); ${g.diametro_furo > 0 ? `${g.furos_por_bloco} furo(s) de ${fmt(g.diametro_furo)} mm` : "maciço"}; ${g.amarracao}; ${g.tipo_junta === "argamassa" ? `argamassa de ${fmt(g.espessura_junta)} mm (${g.assentamento === "lateral" ? "só nas paredes externas" : "em toda a face"}${g.juntas_verticais_preenchidas ? "" : ", juntas verticais vazias"})` : `junta seca de ${fmt(g.espessura_junta)} mm`}`],
      ["Análise", rotuloOpcao("modelo", "analise")],
      PAREDES.includes(projeto.modelo.analise)
        ? ["Parede", `${fmt(g.blocos_por_fiada * g.comprimento_bloco)} × ${fmt(g.n_fiadas * g.altura_fiada)} × ${fmt(g.largura)} mm (H/t = ${fmt(Math.round(g.n_fiadas * g.altura_fiada / g.largura * 10) / 10)}); apoios ${en.apoio_base} na base e ${en.apoio_topo} no topo; excentricidades ${fmt(en.excentricidade_topo)} e ${fmt(en.excentricidade_base)} mm; imperfeição ${en.imperfeicao > 0 ? "H/" + fmt(en.imperfeicao) : "nenhuma"}`]
        : ["Elemento", rotuloOpcao("modelo", "elemento")],
      ["Ensaio", rotuloOpcao("ensaio", "tipo") + (["quente", "residual", "transiente"].includes(en.tipo) ? ` — ${fmt(en.T_alvo)} °C, ${fmt(en.taxa)} °C/min, patamar ${fmt(en.patamar)} min` : "") + (en.tipo === "transiente" ? `, carga ${fmt(en.carga)} kN` : "") + (en.tipo === "iso834" ? ` — ${fmt(en.duracao_iso)} min` : "") + (PAREDES.includes(projeto.modelo.analise) && en.tipo === "iso834" ? `, carga ${fmt(en.carga)} kN` : "")],
      ["Malha", comSens() && projeto.sensibilidade.malha_rapida ? `malha rápida da triagem: seção com h = ${fmt(Math.max(ma.h, 10))} mm, ${Math.min(ma.n_elementos, 16)} elementos de viga, Δt = ${fmt(ma.dt)} s`
        : PAREDES.includes(projeto.modelo.analise) ? `seção com h = ${fmt(ma.h)} mm, ${ma.n_elementos} elementos de viga, Δt = ${fmt(ma.dt)} s`
        : projeto.modelo.analise === "secao_2d" ? `h = ${fmt(ma.h)} mm, Δt = ${fmt(ma.dt)} s` : `h = ${fmt(ma.h)} mm, dz = ${fmt(ma.dz)} mm, Δt = ${fmt(ma.dt)} s`],
    ];
    if (comSens()) {
      const se = projeto.sensibilidade;
      const k = Object.keys(se).filter(c => c.startsWith("var_") && se[c]).length;
      const n = se.trajetorias * (k + 1) + 1;
      itens.push(["Sensibilidade", `${k} parâmetros, ${se.trajetorias} trajetórias: ${n} rodadas da parede` +
        ` (cerca de ${se.malha_rapida ? "10 a 30 s" : "1 a 3 min"} cada uma por processo; ${se.processos > 0 ? se.processos + " processo(s)" : "processos: automático"})`]);
    }
    if (projeto.modelo.analise === "nao_linear") itens.push(["Materiais", `bloco: ft = ${fmt(projeto.modelo.ft)} MPa, fc = ${fmt(projeto.modelo.fc)} MPa` + (g.tipo_junta === "argamassa" ? `; argamassa: fc = ${fmt(projeto.ceramica.arg_fc)} MPa, E = ${fmt(projeto.ceramica.arg_E)} MPa` : `; folga ${fmt(projeto.modelo.desvio_folga)} mm`)]);
    $("#resumo-projeto").innerHTML = itens.map(([a, b]) => `<dt>${esc(a)}</dt><dd>${esc(b)}</dd>`).join("");
    const r = await api("/api/validar", projeto).catch(e => ({ erros: [e.message] }));
    $("#erros-validacao").innerHTML = (r.erros || []).map(e => `<li>${esc(e)}</li>`).join("");
    $("#bt-executar").disabled = !!(r.erros && r.erros.length) || (tarefa && tarefa.estado === "rodando");
  }

  // ---------------------------------------------------------------- execução
  async function executar() {
    projeto.nome = $("#nome-projeto").value || "Projeto";
    const r = await api("/api/executar", projeto).catch(e => ({ erros: [e.message] }));
    if (r.erros) { $("#erros-validacao").innerHTML = r.erros.map(e => `<li>${esc(e)}</li>`).join(""); return; }
    tarefa = { id: r.id, estado: "rodando", pasta: r.pasta }; desde = 0;
    $("#registro").textContent = "";
    $("#bt-executar").disabled = true; $("#bt-cancelar").hidden = false;
    acompanhar();
  }

  function acompanhar() {
    clearInterval(temporizador);
    temporizador = setInterval(async () => {
      let r;
      try { r = await api(`/api/tarefa/${tarefa.id}?desde=${desde}`); } catch (e) { return; }
      desde = r.n;
      const reg = $("#registro");
      const noFim = reg.scrollTop + reg.clientHeight >= reg.scrollHeight - 20;
      if (r.log.length) { reg.textContent += r.log.join("\n") + "\n"; if (noFim) reg.scrollTop = reg.scrollHeight; }
      const min = Math.floor(r.tempo_s / 60), seg = Math.round(r.tempo_s % 60);
      $("#cronometro").textContent = `${min} min ${String(seg).padStart(2, "0")} s`;
      const s = $("#situacao");
      s.className = "situacao " + (r.estado === "concluída" ? "concluida" : r.estado === "erro" ? "erro" : "rodando");
      s.textContent = { "rodando": "Rodando…", "na fila": "Na fila…", "concluída": "Concluída.", "erro": "Erro: " + (r.erro || ""), "cancelada": "Cancelada." }[r.estado] || r.estado;
      $("#estado-menu").textContent = r.estado === "rodando" ? `Análise rodando (${min} min)` : "";
      tarefa.estado = r.estado;
      if (r.estado !== "rodando" && r.estado !== "na fila") {
        clearInterval(temporizador);
        $("#bt-cancelar").hidden = true; $("#bt-executar").disabled = false;
        if (r.estado === "concluída") { await abrirResultado(r.pasta); irPara("resultados"); }
      }
    }, 1000);
  }

  // ---------------------------------------------------------------- resultados
  async function abrirResultado(pasta) {
    const r = await api(`/api/resultado/${encodeURIComponent(pasta)}`);
    if (r.erro) { alert(r.erro); return; }
    resultadoAtual = r; pastaAtual = pasta;
    mostrarResultado();
    if (r.visual) carregarVisual(pasta);
  }

  function mostrarResultado() {
    const r = resultadoAtual, p = pastaAtual, base = `/arquivos/${encodeURIComponent(p)}/`;
    $("#res-origem").innerHTML = `${esc(r.projeto && r.projeto.nome || "")} · pasta <code>${esc(p)}</code>`;
    let h = `<div class="resumo-res">` + r.resumo.map(([rot, v, u]) =>
      `<div class="indicador"><div class="r">${esc(rot)}</div><div class="v">${esc(fmtRes(v))}<small>${esc(u || "")}</small></div></div>`).join("") + `</div>`;
    if (r.avisos && r.avisos.length) h += `<ul class="avisos">${r.avisos.map(a => `<li>${esc(a)}</li>`).join("")}</ul>`;
    if (r.visual) h += `<p><button id="bt-ir-pos">◈ Abrir no pós-processador</button></p>`;
    if (r.figuras.length) h += `<div class="figuras">${r.figuras.map(f => `<figure><img src="${base}${encodeURIComponent(f.arquivo)}" alt="${esc(f.titulo)}"><figcaption>${esc(f.titulo)}</figcaption></figure>`).join("")}</div>`;
    h += `<div class="cartao" style="margin-top:16px"><h2>Tabelas</h2>` + r.tabelas.map(t =>
      `<details class="tab-res"><summary>${esc(t.titulo)} (${t.linhas.length} linhas)</summary><div class="rolagem"><table class="tabela"><thead><tr>${t.colunas.map(c => `<th>${esc(c)}</th>`).join("")}</tr></thead><tbody>${t.linhas.slice(0, 400).map(l => `<tr>${l.map(v => `<td>${esc(fmtRes(v))}</td>`).join("")}</tr>`).join("")}</tbody></table></div></details>`).join("") + `</div>`;
    h += `<div class="cartao"><h2>Arquivos</h2><div class="arquivos">${r.arquivos.map(a => `<a href="${base}${encodeURIComponent(a.nome)}?baixar=1" title="${esc(a.descricao)}">${esc(a.nome)}</a>`).join("")}<a href="${base}projeto.json?baixar=1">projeto.json</a><a href="${base}resultados.json?baixar=1">resultados.json</a></div></div>`;
    $("#res-conteudo").innerHTML = h;
    const b = $("#bt-ir-pos"); if (b) b.addEventListener("click", () => irPara("pos"));
  }

  async function carregarHistorico() {
    const r = await api("/api/historico");
    $("#pasta-res").textContent = r.pasta;
    const nomes = { termica: "térmica 3D", linear: "linear 3D", nao_linear: "não linear 3D", secao_2d: "seção 2D", parede: "parede", sensibilidade: "sensibilidade (parede)" };
    const t = $("#tab-historico");
    if (!r.execucoes.length) { t.innerHTML = "<tr><td>Nenhuma análise gravada ainda.</td></tr>"; return; }
    t.innerHTML = `<thead><tr><th>Data e projeto</th><th>Análise</th><th>Ensaio</th><th>Situação</th><th></th></tr></thead><tbody>` +
      r.execucoes.map(e => `<tr><td>${esc(e.pasta.slice(0, 8))} ${esc(e.pasta.slice(9, 11))}:${esc(e.pasta.slice(11, 13))} · ${esc(e.nome)}</td><td>${esc(nomes[e.analise] || e.analise)}</td><td>${esc(e.ensaio)}</td><td>${e.concluida ? "concluída" : "sem resultados"}</td>
        <td>${e.concluida ? `<button class="secundario" data-acao="res" data-pasta="${esc(e.pasta)}">Resultados</button> ` : ""}${e.visual ? `<button class="secundario" data-acao="pos" data-pasta="${esc(e.pasta)}">Pós</button> ` : ""}<button class="secundario" data-acao="proj" data-pasta="${esc(e.pasta)}">Usar entradas</button></td></tr>`).join("") + "</tbody>";
    t.querySelectorAll("button").forEach(b => b.addEventListener("click", async () => {
      const pasta = b.dataset.pasta;
      if (b.dataset.acao === "proj") {
        const p = await (await fetch(`/arquivos/${encodeURIComponent(pasta)}/projeto.json`)).json();
        definirProjeto(p); irPara("geometria");
      } else { await abrirResultado(pasta); irPara(b.dataset.acao === "pos" ? "pos" : "resultados"); }
    }));
  }

  // ---------------------------------------------------------------- pós-processador
  function iniciarVisualizador() {
    if (vis) return vis;
    try { vis = new Visualizador($("#pos-canvas")); }
    catch (e) { $("#pos-vazio").textContent = "Não foi possível iniciar o WebGL: " + e.message; return null; }
    vis.aoMudarEscala = atualizarLegenda;
    $("#pos-canvas").addEventListener("dblclick", ev => {
      const rc = ev.target.getBoundingClientRect();
      const s = vis.sondar(ev.clientX - rc.left, ev.clientY - rc.top);
      const caixa = $("#pos-sonda");
      if (!s) { caixa.hidden = true; return; }
      caixa.hidden = false;
      caixa.innerHTML = `<b>${esc(s.campo)}:</b> ${s.valor === null || s.valor !== s.valor ? "—" : esc(fmtRes(s.valor))} ${esc(s.unidade)}<br>` +
        `nó ${s.no}, elemento ${s.elemento} (${esc(s.material)})<br>x = ${fmtRes(s.x)}, y = ${fmtRes(s.y)}, z = ${fmtRes(s.z)} mm`;
    });
    return vis;
  }

  async function carregarVisual(pasta) {
    if (!iniciarVisualizador()) return;
    $("#pos-origem").textContent = "Carregando…";
    const d = await (await fetch(`/arquivos/${encodeURIComponent(pasta)}/visual.json`)).json();
    vis.carregar(d);
    $("#pos-vazio").hidden = true;
    $("#pos-origem").textContent = `${d.descricao} · ${pasta}`;
    const se = $("#pos-estado");
    se.innerHTML = d.estados.map((e, i) => `<option value="${i}">${esc(e.rotulo)}</option>`).join("");
    se.value = vis.op.estado;
    const sl = $("#pos-estado-slider"); sl.max = Math.max(0, d.estados.length - 1); sl.value = vis.op.estado;
    $("#pos-campo").innerHTML = d.campos.map(c => `<option value="${c.id}">${esc(c.nome)}${c.unidade ? " (" + esc(c.unidade) + ")" : ""}</option>`).join("") +
      `<option value="material">Material</option>`;
    $("#pos-campo").value = vis.op.campo || "material";
    const usados = [...new Set(d.material)].sort();
    $("#pos-materiais").innerHTML = usados.map(m => `<label class="marcar"><input type="checkbox" data-mat="${m}" ${vis.op.ocultos.has(m) ? "" : "checked"}> <span style="display:inline-block;width:10px;height:10px;border-radius:2px;background:${COR_MATERIAL[m] || "#bbb"}"></span> ${esc(d.nomes_materiais[m] || m)}</label>`).join("");
    $$("#pos-materiais input").forEach(c => c.addEventListener("change", () => {
      const m = +c.dataset.mat; if (c.checked) vis.op.ocultos.delete(m); else vis.op.ocultos.add(m); vis.atualizarGeometria();
    }));
    $("#pos-deformada").disabled = !d.estados.some(e => e.u);
    $("#pos-espelhar").disabled = !(d.simetria && d.tipo === "3d");
    atualizarLegenda();
  }

  function atualizarLegenda() {
    if (!vis || !vis.dados) return;
    const id = vis.op.campo, info = vis.infoCampo(id);
    const leg = $("#pos-legenda"); leg.hidden = false;
    if (!id || id === "material") {
      $("#legenda-titulo").textContent = "Material";
      $("#legenda-barra").style.background = "linear-gradient(90deg," + Object.values(COR_MATERIAL).join(",") + ")";
      $("#legenda-min").textContent = ""; $("#legenda-meio").textContent = ""; $("#legenda-max").textContent = "";
      return;
    }
    $("#legenda-titulo").textContent = info.nome + (info.unidade ? ` (${info.unidade})` : "") + " — " + vis.estados[vis.op.estado].rotulo;
    $("#legenda-barra").style.background = vis.gradienteCSS();
    $("#legenda-min").textContent = fmtRes(vis.op.min);
    $("#legenda-meio").textContent = fmtRes((vis.op.min + vis.op.max) / 2);
    $("#legenda-max").textContent = fmtRes(vis.op.max);
    if (vis.op.auto) { $("#pos-min").value = +vis.op.min.toPrecision(4); $("#pos-max").value = +vis.op.max.toPrecision(4); }
  }

  function ligarControlesPos() {
    const ok = () => vis && vis.dados;
    $("#pos-estado").addEventListener("change", e => { if (!ok()) return; vis.op.estado = +e.target.value; $("#pos-estado-slider").value = e.target.value; vis.op.deformada ? vis.atualizarPosicoes() : null; vis.atualizarCores(); });
    $("#pos-estado-slider").addEventListener("input", e => { if (!ok()) return; vis.op.estado = +e.target.value; $("#pos-estado").value = e.target.value; if (vis.op.deformada) vis.atualizarPosicoes(); vis.atualizarCores(); });
    $("#pos-campo").addEventListener("change", e => { if (!ok()) return; vis.op.campo = e.target.value; vis.atualizarCores(); });
    $("#pos-mapa").addEventListener("change", e => { if (!ok()) return; vis.op.mapa = e.target.value; vis.atualizarCores(); });
    $("#pos-auto").addEventListener("change", e => { if (!ok()) return; vis.op.auto = e.target.checked; vis.atualizarCores(); });
    for (const k of ["min", "max"]) $("#pos-" + k).addEventListener("change", e => {
      if (!ok()) return; vis.op.auto = false; $("#pos-auto").checked = false; vis.op[k] = num(e.target.value); vis.atualizarCores();
    });
    $("#pos-deformada").addEventListener("change", e => { if (!ok()) return; vis.op.deformada = e.target.checked; vis.atualizarPosicoes(); });
    $("#pos-escala").addEventListener("change", e => { if (!ok()) return; vis.op.escala = num(e.target.value) || 0; if (vis.op.deformada) vis.atualizarPosicoes(); });
    $("#pos-espelhar").addEventListener("change", e => { if (!ok()) return; vis.op.espelhar = e.target.checked; vis.atualizarGeometria(); vis.vista("iso"); });
    $("#pos-arestas").addEventListener("change", e => { if (!ok()) return; vis.op.arestas = e.target.checked; vis.desenhar(); });
    for (const a of ["x", "y", "z"]) {
      const f = () => { if (!ok()) return; vis.op.cortes[a] = $("#corte-" + a).checked ? +$("#corte-" + a + "-v").value / 100 : null; vis.atualizarGeometria(); };
      $("#corte-" + a).addEventListener("change", f); $("#corte-" + a + "-v").addEventListener("input", f);
    }
    $$("[data-vista]").forEach(b => b.addEventListener("click", () => vis && vis.vista(b.dataset.vista)));
    $("#pos-png").addEventListener("click", () => {
      if (!ok()) return;
      const a = document.createElement("a"); a.href = vis.png(); a.download = `${pastaAtual || "mef"}_${vis.op.campo}.png`; a.click();
    });
  }

  // ---------------------------------------------------------------- início
  async function iniciar() {
    const r = await api("/api/esquema");
    ESQ = r.esquema; PADRAO = r.padrao; PRED = r.predefinidos; BIB = r.biblioteca || {};
    $("#versao").textContent = "versão " + r.versao;
    const sel = $("#predefinidos");
    for (const k of Object.keys(PRED)) { const o = document.createElement("option"); o.value = k; o.textContent = k; sel.appendChild(o); }
    sel.addEventListener("change", () => { if (sel.value) { definirProjeto({ ...PRED[sel.value], nome: sel.value }); sel.value = ""; } });
    let salvo = null; try { salvo = JSON.parse(localStorage.getItem("mef_masonry_projeto") || localStorage.getItem("mef_bcepq_projeto")); } catch (e) { /* ignora */ }
    definirProjeto(salvo || PADRAO);
    $$("#menu button").forEach(b => b.addEventListener("click", () => irPara(b.dataset.pagina)));
    $("#bt-proximo").addEventListener("click", () => { const o = ordem(); irPara(o[Math.min(o.length - 1, o.indexOf(pagina) + 1)]); });
    $("#bt-anterior").addEventListener("click", () => { const o = ordem(); irPara(o[Math.max(0, o.indexOf(pagina) - 1)]); });
    $("#bt-previa").addEventListener("click", previaMalha);
    $("#bt-executar").addEventListener("click", executar);
    $("#bt-cancelar").addEventListener("click", () => tarefa && api(`/api/cancelar/${tarefa.id}`, {}));
    $("#nome-projeto").addEventListener("change", e => { projeto.nome = e.target.value; salvarLocal(); });
    $("#bt-novo").addEventListener("click", () => definirProjeto(PADRAO));
    $("#bt-salvar").addEventListener("click", () => {
      projeto.nome = $("#nome-projeto").value;
      const a = document.createElement("a");
      a.href = URL.createObjectURL(new Blob([JSON.stringify(projeto, null, 2)], { type: "application/json" }));
      a.download = (projeto.nome || "projeto").replace(/[^\w\-]+/g, "_") + ".json"; a.click();
    });
    $("#bt-abrir").addEventListener("click", () => $("#arquivo-abrir").click());
    $("#arquivo-abrir").addEventListener("change", async e => {
      const f = e.target.files[0]; if (!f) return;
      try { definirProjeto(JSON.parse(await f.text())); } catch (err) { alert("Arquivo inválido: " + err.message); }
      e.target.value = "";
    });
    ligarControlesPos();
    irPara("geometria");
    // retoma uma análise em andamento (página recarregada)
    const at = await api("/api/atual").catch(() => ({}));
    if (at.id && (at.estado === "rodando" || at.estado === "na fila")) {
      tarefa = { id: at.id, estado: at.estado, pasta: at.pasta }; desde = 0; $("#bt-cancelar").hidden = false; acompanhar();
    }
  }
  iniciar().catch(e => { document.body.insertAdjacentHTML("afterbegin", `<p style="color:#c00;padding:12px">Erro ao iniciar: ${esc(e.message)}</p>`); });
})();
