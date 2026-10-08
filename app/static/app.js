// DNA Analysis web UI: a plain ES module with no build step or framework.
// All user-controlled text (FASTA headers, sequences) is inserted with
// textContent via the h() helper, never as HTML.

const SAMPLES = ["patient_01", "patient_02", "patient_03", "patient_04", "patient_05"];

const MODES = {
  align: {
    button: "Run alignment",
    hint: "Global alignment to REF_001 with SNV, insertion and deletion calls.",
  },
  translate: {
    button: "Translate",
    hint: "Standard genetic code, stops before the first stop codon. Codons containing N become X.",
  },
  analyze: {
    button: "Run full analysis",
    hint: "Alignment + frame-1 translation + effect of each variant (synonymous, missense, nonsense, frameshift…).",
  },
};

const EFFECT_STATUS = {
  synonymous: "ok",
  missense: "warn",
  inframe_insertion: "warn",
  inframe_deletion: "warn",
  nonsense: "bad",
  stop_lost: "bad",
  frameshift: "bad",
  undetermined: "neutral",
};

// `pending` holds an in-flight sample download so Run never races it.
const state = { source: "upload", file: null, sample: null, mode: "align", frame: 1, pending: null };

const $ = (selector, root = document) => root.querySelector(selector);
const results = $("#results");

// ---------------------------------------------------------------------------- DOM helpers

function h(tag, props = {}, ...children) {
  const el = document.createElement(tag);
  for (const [key, value] of Object.entries(props)) {
    if (value == null || value === false) continue;
    if (key === "class") el.className = value;
    else if (key === "dataset") Object.assign(el.dataset, value);
    else if (key.startsWith("on")) el.addEventListener(key.slice(2), value);
    else el.setAttribute(key, value === true ? "" : value);
  }
  for (const child of children.flat(Infinity)) {
    if (child == null || child === false) continue;
    el.append(child instanceof Node ? child : String(child));
  }
  return el;
}

const plural = (n, word) => `${n} ${word}${n === 1 ? "" : "s"}`;
const formatBytes = (n) => (n < 1024 ? `${n} B` : `${(n / 1024).toFixed(1)} KB`);
const formatPercent = (x) => `${(x * 100).toFixed(2)}%`;
const ungapped = (s) => s.replaceAll("-", "").length;

/** Width of one monospace character inside `container`, cached per font size. */
const charWidthCache = new Map();
function charWidth(container) {
  const font = getComputedStyle(container).font;
  if (!charWidthCache.has(font)) {
    const probe = h("span", { style: "position:absolute;visibility:hidden;white-space:pre" }, "0".repeat(50));
    container.append(probe);
    charWidthCache.set(font, probe.getBoundingClientRect().width / 50 || 8);
    probe.remove();
  }
  return charWidthCache.get(font);
}

/** Re-render `draw()` whenever the element's width changes (sequence viewers wrap to fit). */
function fitToWidth(element, draw) {
  let lastWidth = -1;
  const observer = new ResizeObserver(([entry]) => {
    const width = Math.round(entry.contentRect.width);
    if (width !== lastWidth) {
      lastWidth = width;
      draw(width);
    }
  });
  observer.observe(element);
  return element;
}

// ---------------------------------------------------------------------------- controls

function bindSegment(attribute, onSelect) {
  const buttons = document.querySelectorAll(`[data-${attribute}]`);
  for (const button of buttons) {
    button.addEventListener("click", () => {
      for (const other of buttons) other.setAttribute("aria-selected", String(other === button));
      onSelect(button.dataset[attribute]);
    });
  }
}

function setSource(source) {
  state.source = source;
  $("#upload-pane").hidden = source !== "upload";
  $("#paste-pane").hidden = source !== "paste";
  for (const b of document.querySelectorAll("[data-source]")) {
    b.setAttribute("aria-selected", String(b.dataset.source === source));
  }
}

function setMode(mode) {
  state.mode = mode;
  $("#run-btn").textContent = MODES[mode].button;
  $("#mode-hint").textContent = MODES[mode].hint;
  $("#frame-control").hidden = mode !== "translate";
}

function setFile(file, sample = null) {
  state.file = file;
  state.sample = sample;
  $("#dropzone").hidden = Boolean(file);
  $("#file-pill").hidden = !file;
  if (file) {
    $("#file-name").textContent = file.name;
    $("#file-size").textContent = formatBytes(file.size);
  } else {
    $("#file-input").value = "";
  }
  for (const chip of document.querySelectorAll("[data-sample]")) {
    chip.setAttribute("aria-pressed", String(chip.dataset.sample === sample));
  }
}

function loadSample(name) {
  state.pending = (async () => {
    try {
      const response = await fetch(`/samples/${name}.fasta`);
      if (!response.ok) throw new Error(`HTTP ${response.status}`);
      const blob = await response.blob();
      setSource("upload");
      setFile(new File([blob], `${name}.fasta`, { type: "text/plain" }), name);
    } catch (error) {
      renderError(0, { code: "sample_unavailable", message: `Could not load ${name}: ${error.message}` });
    } finally {
      state.pending = null;
    }
  })();
  return state.pending;
}

function currentFile() {
  if (state.source === "paste") {
    // Sent even when empty: the API's own validation produces the error message.
    return new File([$("#paste-input").value], "pasted.fasta", { type: "text/plain" });
  }
  return state.file; // null -> the API answers 400 missing_file
}

function setBusy(busy) {
  const button = $("#run-btn");
  button.disabled = busy;
  button.setAttribute("aria-busy", String(busy));
  results.setAttribute("aria-busy", String(busy));
}

async function run() {
  if (state.pending) await state.pending;
  const form = new FormData();
  const file = currentFile();
  if (file) form.append("file", file, file.name);
  if (state.mode === "translate") form.append("frame", String(state.frame));

  setBusy(true);
  try {
    const response = await fetch(`/api/${state.mode}`, { method: "POST", body: form });
    const body = await response.json().catch(() => null);
    if (!response.ok) {
      renderError(response.status, body?.error);
    } else {
      RENDERERS[state.mode](body);
    }
  } catch {
    renderError(0, { code: "network_error", message: "Could not reach the API. Is the server running?" });
  } finally {
    setBusy(false);
  }
}

// ---------------------------------------------------------------------------- shared views

function resultHead(body, kind) {
  const copyButton = h("button", { class: "btn btn--ghost btn--sm", type: "button" }, "Copy JSON");
  copyButton.addEventListener("click", async () => {
    try {
      await navigator.clipboard.writeText(JSON.stringify(body, null, 2));
      copyButton.textContent = "Copied";
    } catch {
      copyButton.textContent = "Copy failed";
    }
    setTimeout(() => (copyButton.textContent = "Copy JSON"), 1500);
  });
  const download = () => {
    const blob = new Blob([JSON.stringify(body, null, 2)], { type: "application/json" });
    const link = h("a", { href: URL.createObjectURL(blob), download: `${body.sequence_id}_${kind.toLowerCase()}.json` });
    link.click();
    URL.revokeObjectURL(link.href);
  };
  return h(
    "div",
    { class: "result-head" },
    h(
      "div",
      { class: "result-head__title" },
      h("span", { class: "eyebrow" }, `${kind}${body.reference_id ? ` · vs ${body.reference_id}` : ""}`),
      h("h2", { class: "h2" }, body.sequence_id),
    ),
    h(
      "div",
      { class: "result-head__actions" },
      copyButton,
      h("button", { class: "btn btn--ghost btn--sm", type: "button", onclick: download }, "Download"),
    ),
  );
}

function statGrid(stats) {
  return h(
    "div",
    { class: "stat-grid" },
    stats.map(([label, value, sub, tone]) =>
      h(
        "div",
        { class: "stat" },
        h("div", { class: "stat__label" }, label),
        h("div", { class: `stat__value${tone ? ` stat__value--${tone}` : ""}` }, value),
        sub ? h("div", { class: "stat__sub" }, sub) : null,
      ),
    ),
  );
}

function sectionCard(title, aside, ...body) {
  return h(
    "section",
    { class: "card section-card" },
    h("div", { class: "section-card__head" }, h("h3", { class: "title" }, title), aside),
    h("div", { class: "section-card__body" }, body),
  );
}

function rawJson(body) {
  return h(
    "details",
    { class: "card section-card raw" },
    h("summary", {}, "Raw JSON response"),
    h("pre", {}, JSON.stringify(body, null, 2)),
  );
}

function renderError(status, error) {
  const err = error ?? { code: "unexpected_response", message: "The server returned an unexpected response." };
  results.replaceChildren(
    h(
      "div",
      { class: "card error-card", role: "alert" },
      h("span", { class: "eyebrow" }, status ? `Request rejected · HTTP ${status}` : "Request failed"),
      h("h2", { class: "h3" }, "This input could not be processed"),
      h("p", { class: "body" }, err.message),
      h(
        "div",
        {},
        h("span", { class: "status status--bad error-card__code" }, err.code),
        err.field ? h("span", { class: "caption" }, `  field: ${err.field}`) : null,
      ),
    ),
  );
}

// ---------------------------------------------------------------------------- alignment viewer

/** Map every variant to the alignment columns it covers (for click-to-highlight). */
function variantColumns(referenceAligned, variant) {
  const columnOf = [-1]; // columnOf[p] = column holding reference base p (1-based)
  for (let i = 0; i < referenceAligned.length; i++) {
    if (referenceAligned[i] !== "-") columnOf.push(i);
  }
  const range = (start, count) => Array.from({ length: count }, (_, k) => start + k);
  if (variant.type === "insertion") return range(columnOf[variant.position] + 1, variant.alt.length);
  if (variant.type === "deletion") return range(columnOf[variant.position], variant.ref.length);
  return [columnOf[variant.position]];
}

function columnClass(r, q) {
  if (r === "-" || q === "-") return "c-g";
  if (r === "N" || q === "N") return "c-n";
  return r === q ? "" : "c-x";
}

function alignmentView({ reference_aligned: ref, query_aligned: qry }) {
  const container = h("div", { class: "aln" });
  const draw = (width) => {
    const perRow = Math.max(20, Math.min(120, Math.floor((width - 130) / charWidth(container) / 10) * 10));
    const blocks = [];
    let refPos = 0;
    let qryPos = 0;
    for (let start = 0; start < ref.length; start += perRow) {
      const end = Math.min(ref.length, start + perRow);
      const refRow = h("span", { class: "aln__seq" });
      const qryRow = h("span", { class: "aln__seq" });
      let mid = "";
      const refStart = refPos + 1;
      const qryStart = qryPos + 1;
      for (let i = start; i < end; i++) {
        const r = ref[i];
        const q = qry[i];
        if (r !== "-") refPos++;
        if (q !== "-") qryPos++;
        const cls = columnClass(r, q);
        refRow.append(h("span", { class: cls || null, dataset: { col: i } }, r));
        qryRow.append(h("span", { class: cls || null, dataset: { col: i } }, q));
        mid += cls === "" ? "|" : cls === "c-x" ? "." : " ";
      }
      blocks.push(
        h(
          "div",
          { class: "aln__block" },
          h("span", { class: "aln__label" }, "REF"),
          h("span", { class: "aln__num" }, Math.min(refStart, refPos)),
          refRow,
          h("span", { class: "aln__num" }, refPos),
          h("span", {}),
          h("span", {}),
          h("span", { class: "aln__seq aln__mid" }, mid),
          h("span", {}),
          h("span", { class: "aln__label" }, "QRY"),
          h("span", { class: "aln__num" }, Math.min(qryStart, qryPos)),
          qryRow,
          h("span", { class: "aln__num" }, qryPos),
        ),
      );
      if (end === ref.length) break;
    }
    container.replaceChildren(...blocks);
  };
  return fitToWidth(container, draw);
}

function alignmentCard(alignment) {
  const legend = h(
    "div",
    { class: "legend caption" },
    h("span", {}, h("i", { class: "swatch c-x" }), "mismatch"),
    h("span", {}, h("i", { class: "swatch c-g" }), "gap"),
    h("span", {}, h("i", { class: "swatch", style: "background:var(--ink-4)" }), "N (unknown)"),
  );
  return sectionCard("Alignment", legend, alignmentView(alignment));
}

function flashColumns(columns) {
  const cells = columns.flatMap((c) => [...results.querySelectorAll(`.aln [data-col="${c}"]`)]);
  if (!cells.length) return;
  cells[0].scrollIntoView({ behavior: "smooth", block: "center", inline: "center" });
  for (const cell of cells) {
    cell.classList.remove("flash");
    void cell.offsetWidth; // restart the animation
    cell.classList.add("flash");
  }
}

// ---------------------------------------------------------------------------- variants table

function variantsCard(variants, referenceAligned, { annotated = false } = {}) {
  if (!variants.length) {
    return sectionCard(
      "Variants",
      null,
      h("p", { class: "body" }, "No variants: every A/C/G/T position agrees with the reference."),
    );
  }
  const counts = variants.reduce((acc, v) => ({ ...acc, [v.type]: (acc[v.type] ?? 0) + 1 }), {});
  const summary = h("span", { class: "caption" }, Object.entries(counts).map(([t, n]) => plural(n, t)).join(" · "));
  const head = ["#", "Position", "Type", "Ref", "Alt"];
  if (annotated) head.push("Effect", "Protein", "Notes");

  const rows = variants.map((v, index) => {
    const cells = [
      h("td", { class: "tabular caption" }, index + 1),
      h("td", { class: "tabular" }, v.position),
      h("td", {}, v.type),
      h("td", { class: "mono", title: v.ref }, v.ref),
      h("td", { class: "mono", title: v.alt }, v.alt),
    ];
    if (annotated) {
      cells.push(
        h("td", {}, h("span", { class: `status status--${EFFECT_STATUS[v.effect] ?? "neutral"}` }, v.effect.replace("_", " "))),
        h("td", { class: "mono" }, v.protein_change ?? "–"),
        h("td", { class: "notes" }, v.notes.length ? v.notes.join(" ") : ""),
      );
    }
    const locate = () => flashColumns(variantColumns(referenceAligned, v));
    return h(
      "tr",
      {
        tabindex: "0",
        title: "Show in alignment",
        onclick: locate,
        onkeydown: (event) => {
          if (event.key === "Enter" || event.key === " ") {
            event.preventDefault();
            locate();
          }
        },
      },
      cells,
    );
  });

  return sectionCard(
    "Variants",
    summary,
    h(
      "div",
      { class: "table-wrap" },
      h("table", { class: "vtable" }, h("thead", {}, h("tr", {}, head.map((t) => h("th", {}, t)))), h("tbody", {}, rows)),
    ),
  );
}

// ---------------------------------------------------------------------------- protein viewer

function proteinView(protein, hasStop) {
  const container = h("div", { class: "protein" });
  const draw = (width) => {
    // 10-residue blocks separated by a space, numbered like UniProt.
    const blocksPerLine = Math.max(1, Math.min(6, Math.floor((width - 60) / (11 * charWidth(container)))));
    const perLine = blocksPerLine * 10;
    const lines = [];
    for (let start = 0; start < Math.max(protein.length, 1); start += perLine) {
      const seq = h("span", {});
      for (let b = start; b < Math.min(protein.length, start + perLine); b += 10) {
        if (b > start) seq.append(" ");
        for (const aa of protein.slice(b, b + 10)) seq.append(aa === "X" ? h("span", { class: "aa-x" }, aa) : aa);
      }
      if (start + perLine >= protein.length && hasStop) {
        seq.append(h("span", { class: "protein__stop", title: "stop codon (not part of the protein)" }, " stop"));
      }
      lines.push(h("div", { class: "protein__line" }, h("span", { class: "protein__num" }, start + 1), seq));
    }
    container.replaceChildren(...lines);
  };
  return fitToWidth(container, draw);
}

function emptyProteinMessage(body) {
  return body.has_stop_codon
    ? "The first codon in this frame is a stop codon, so the protein is empty."
    : "This frame leaves no complete codon to translate.";
}

// ---------------------------------------------------------------------------- page renderers

function renderAlign(body) {
  const { alignment, variants } = body;
  results.replaceChildren(
    resultHead(body, "Alignment"),
    statGrid([
      ["Score", alignment.score.toFixed(1), "+1 / −1 / −2 per gap · N = 0"],
      ["Identity", formatPercent(alignment.identity), "identical A/C/G/T ÷ all columns"],
      ["Variants", variants.length, variants.length ? "click a row to locate it" : "matches the reference", variants.length ? null : "ok"],
      [
        "Columns",
        alignment.reference_aligned.length,
        `${ungapped(alignment.reference_aligned)} ref · ${ungapped(alignment.query_aligned)} query nt`,
      ],
    ]),
    variantsCard(variants, alignment.reference_aligned),
    alignmentCard(alignment),
    rawJson(body),
  );
}

function renderTranslate(body) {
  const unknown = [...body.protein].filter((aa) => aa === "X").length;
  results.replaceChildren(
    resultHead(body, "Translation"),
    statGrid([
      ["Frame", body.frame, `starts at nucleotide ${body.frame}`],
      ["Length", body.length, "amino acids (stop excluded)"],
      ["Stop codon", body.has_stop_codon ? "Yes" : "No", body.has_stop_codon ? "translation ended at a stop" : "ran to the last complete codon"],
      ["Unknown (X)", unknown, "codons containing N", unknown ? null : "ok"],
    ]),
    sectionCard(
      "Protein",
      h("span", { class: "caption" }, plural(body.length, "residue")),
      body.length ? proteinView(body.protein, body.has_stop_codon) : h("p", { class: "body" }, emptyProteinMessage(body)),
    ),
    rawJson(body),
  );
}

function renderAnalyze(body) {
  const { alignment, variants, translation, reference_translation: reference } = body;
  const damaging = variants.filter((v) => EFFECT_STATUS[v.effect] === "bad").length;
  results.replaceChildren(
    resultHead(body, "Full analysis"),
    statGrid([
      ["Identity", formatPercent(alignment.identity), `score ${alignment.score.toFixed(1)}`],
      ["Variants", variants.length, plural(damaging, "high-impact change"), damaging ? "bad" : null],
      ["Protein", `${translation.length}`, `aa vs ${reference.length} aa in REF`, translation.length === reference.length ? null : "bad"],
      ["Stop codon", translation.has_stop_codon ? "Yes" : "No", "frame 1 of the sample"],
    ]),
    variantsCard(variants, alignment.reference_aligned, { annotated: true }),
    sectionCard(
      "Protein (frame 1)",
      null,
      h(
        "div",
        { class: "protein-pair" },
        h("div", {}, h("h4", { class: "eyebrow" }, `Reference · ${reference.length} aa`), proteinView(reference.protein, reference.has_stop_codon)),
        h("div", {}, h("h4", { class: "eyebrow" }, `${body.sequence_id} · ${translation.length} aa`),
          translation.length ? proteinView(translation.protein, translation.has_stop_codon) : h("p", { class: "body" }, emptyProteinMessage(translation))),
      ),
    ),
    alignmentCard(alignment),
    rawJson(body),
  );
}

const RENDERERS = { align: renderAlign, translate: renderTranslate, analyze: renderAnalyze };

// ---------------------------------------------------------------------------- wiring

function init() {
  bindSegment("source", setSource);
  bindSegment("mode", setMode);
  bindSegment("frame", (frame) => (state.frame = Number(frame)));
  setMode(state.mode);

  const input = $("#file-input");
  const dropzone = $("#dropzone");
  dropzone.addEventListener("click", () => input.click());
  input.addEventListener("change", () => setFile(input.files[0] ?? null));
  for (const type of ["dragenter", "dragover"]) {
    dropzone.addEventListener(type, (event) => {
      event.preventDefault();
      dropzone.classList.add("is-dragover");
    });
  }
  for (const type of ["dragleave", "drop"]) {
    dropzone.addEventListener(type, () => dropzone.classList.remove("is-dragover"));
  }
  dropzone.addEventListener("drop", (event) => {
    event.preventDefault();
    const file = event.dataTransfer?.files?.[0];
    if (file) setFile(file);
  });
  $("#file-clear").addEventListener("click", () => setFile(null));

  $("#samples").append(
    ...SAMPLES.map((name) =>
      h("button", { type: "button", class: "chip", "aria-pressed": "false", dataset: { sample: name }, onclick: () => loadSample(name) }, name),
    ),
  );

  $("#run-btn").addEventListener("click", run);
  $("#paste-input").addEventListener("keydown", (event) => {
    if (event.key === "Enter" && (event.ctrlKey || event.metaKey)) run();
  });

  $("#theme-toggle").addEventListener("click", () => {
    const next = document.documentElement.dataset.theme === "dark" ? "light" : "dark";
    document.documentElement.dataset.theme = next;
    try {
      localStorage.setItem("dna-theme", next);
    } catch {
      /* storage unavailable: theme still applies for this visit */
    }
  });
}

init();
