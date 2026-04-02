<script>
  import { onMount } from 'svelte';
  import { store } from '../stores/data.svelte.js';

  let { filters = {} } = $props();

  let heatmapData = $state(null);
  let container;
  let canvas;
  let rowDendroSvg;
  let colDendroSvg;
  let tooltip = $state({ show: false, x: 0, y: 0, text: '' });

  // Layout constants
  const ROW_DENDRO_W = 80;
  const COL_DENDRO_H = 60;
  const LABEL_W = 0;  // auto-hide labels for large matrices
  const LABEL_H = 0;

  onMount(async () => {
    try {
      let res = await fetch('/data/heatmap.json.gz');
      if (!res.ok) res = await fetch('/data/heatmap.json');
      if (res.ok) heatmapData = await res.json();
    } catch (e) {
      console.error('Failed to load heatmap data:', e);
    }
  });

  // Viridis-ish colormap (256 entries)
  function viridis(t) {
    // Simplified viridis: dark purple -> teal -> yellow
    const r = Math.round(Math.max(0, Math.min(255, -510 * t * t + 713 * t + 68)));
    const g = Math.round(Math.max(0, Math.min(255, -270 * t * t + 442 * t + 24)));
    const b = Math.round(Math.max(0, Math.min(255, 320 * t * t - 570 * t + 280)));
    return `rgb(${r},${g},${b})`;
  }

  // Viridis sampled at 9 stops, linearly interpolated to 256
  const VIRIDIS_STOPS = [
    [68, 1, 84], [72, 35, 116], [64, 67, 135], [52, 94, 141],
    [33, 145, 140], [53, 183, 121], [109, 205, 89], [180, 222, 44], [253, 231, 37],
  ];

  function buildViridisLUT() {
    const lut = [];
    const nStops = VIRIDIS_STOPS.length;
    for (let i = 0; i < 256; i++) {
      const t = i / 255 * (nStops - 1);
      const lo = Math.floor(t);
      const hi = Math.min(lo + 1, nStops - 1);
      const f = t - lo;
      lut.push([
        Math.round(VIRIDIS_STOPS[lo][0] * (1 - f) + VIRIDIS_STOPS[hi][0] * f),
        Math.round(VIRIDIS_STOPS[lo][1] * (1 - f) + VIRIDIS_STOPS[hi][1] * f),
        Math.round(VIRIDIS_STOPS[lo][2] * (1 - f) + VIRIDIS_STOPS[hi][2] * f),
      ]);
    }
    return lut;
  }

  const viridisLUT = buildViridisLUT();

  $effect(() => {
    if (!heatmapData || !canvas || !rowDendroSvg || !colDendroSvg) return;

    const { z, sampleIds, asvIds, asvLabels, rowDendro, colDendro } = heatmapData;
    const nRows = z.length;
    const nCols = z[0]?.length || 0;
    if (nRows === 0 || nCols === 0) return;

    // Find max value for color scaling
    let zMax = 0;
    for (const row of z) for (const v of row) if (v > zMax) zMax = v;
    if (zMax === 0) zMax = 1;

    // Get container size
    const rect = container.getBoundingClientRect();
    const heatW = rect.width - ROW_DENDRO_W;
    const heatH = rect.height - COL_DENDRO_H;

    const cellW = heatW / nCols;
    const cellH = heatH / nRows;

    // ── Draw heatmap on canvas ──
    canvas.width = heatW * window.devicePixelRatio;
    canvas.height = heatH * window.devicePixelRatio;
    canvas.style.width = heatW + 'px';
    canvas.style.height = heatH + 'px';

    const ctx = canvas.getContext('2d');
    ctx.scale(window.devicePixelRatio, window.devicePixelRatio);

    const imgData = ctx.createImageData(Math.ceil(heatW * window.devicePixelRatio), Math.ceil(heatH * window.devicePixelRatio));
    const pixels = imgData.data;
    const dpr = window.devicePixelRatio;

    // Direct pixel painting for speed
    for (let row = 0; row < nRows; row++) {
      const y0 = Math.floor(row * cellH * dpr);
      const y1 = Math.floor((row + 1) * cellH * dpr);
      for (let col = 0; col < nCols; col++) {
        const t = Math.min(z[row][col] / zMax, 1);
        const ci = Math.round(t * 255);
        const [r, g, b] = viridisLUT[ci];
        const x0 = Math.floor(col * cellW * dpr);
        const x1 = Math.floor((col + 1) * cellW * dpr);
        for (let py = y0; py < y1; py++) {
          for (let px = x0; px < x1; px++) {
            const idx = (py * imgData.width + px) * 4;
            pixels[idx] = r;
            pixels[idx + 1] = g;
            pixels[idx + 2] = b;
            pixels[idx + 3] = 255;
          }
        }
      }
    }
    ctx.putImageData(imgData, 0, 0);

    // ── Draw column dendrogram (SVG, above heatmap) ──
    colDendroSvg.innerHTML = '';
    colDendroSvg.setAttribute('width', heatW);
    colDendroSvg.setAttribute('height', COL_DENDRO_H);

    if (colDendro.icoord.length > 0) {
      // scipy dendrogram: icoord values are 5, 15, 25... (step 10, centered on leaves)
      // Scale to pixel positions
      const maxDist = Math.max(...colDendro.dcoord.flat());
      for (let i = 0; i < colDendro.icoord.length; i++) {
        const ix = colDendro.icoord[i];
        const iy = colDendro.dcoord[i];
        const points = [];
        for (let j = 0; j < 4; j++) {
          const px = ((ix[j] - 5) / (nCols * 10 - 10)) * heatW;
          const py = COL_DENDRO_H - (iy[j] / maxDist) * (COL_DENDRO_H - 4);
          points.push(`${px},${py}`);
        }
        const line = document.createElementNS('http://www.w3.org/2000/svg', 'polyline');
        line.setAttribute('points', points.join(' '));
        line.setAttribute('fill', 'none');
        line.setAttribute('stroke', '#64748b');
        line.setAttribute('stroke-width', '0.8');
        colDendroSvg.appendChild(line);
      }
    }

    // ── Draw row dendrogram (SVG, left of heatmap) ──
    rowDendroSvg.innerHTML = '';
    rowDendroSvg.setAttribute('width', ROW_DENDRO_W);
    rowDendroSvg.setAttribute('height', heatH);

    if (rowDendro.icoord.length > 0) {
      const maxDist = Math.max(...rowDendro.dcoord.flat());
      for (let i = 0; i < rowDendro.icoord.length; i++) {
        const ix = rowDendro.icoord[i];  // leaf positions (vertical)
        const iy = rowDendro.dcoord[i];  // heights (horizontal)
        const points = [];
        for (let j = 0; j < 4; j++) {
          const py = ((ix[j] - 5) / (nRows * 10 - 10)) * heatH;
          const px = ROW_DENDRO_W - (iy[j] / maxDist) * (ROW_DENDRO_W - 4);
          points.push(`${px},${py}`);
        }
        const line = document.createElementNS('http://www.w3.org/2000/svg', 'polyline');
        line.setAttribute('points', points.join(' '));
        line.setAttribute('fill', 'none');
        line.setAttribute('stroke', '#64748b');
        line.setAttribute('stroke-width', '0.8');
        rowDendroSvg.appendChild(line);
      }
    }

    // ── Hover handler ──
    canvas.onmousemove = (e) => {
      const cr = canvas.getBoundingClientRect();
      const mx = e.clientX - cr.left;
      const my = e.clientY - cr.top;
      const col = Math.floor(mx / cellW);
      const row = Math.floor(my / cellH);
      if (row >= 0 && row < nRows && col >= 0 && col < nCols) {
        const val = z[row][col];
        tooltip = {
          show: true,
          x: e.clientX,
          y: e.clientY,
          text: `${sampleIds[row]} | ${asvLabels[col]} | ${val.toFixed(4)}`,
        };
      } else {
        tooltip = { show: false, x: 0, y: 0, text: '' };
      }
    };
    canvas.onmouseleave = () => {
      tooltip = { show: false, x: 0, y: 0, text: '' };
    };
  });
</script>

<div class="flex h-full flex-col">
  {#if !heatmapData}
    <div class="flex-1 flex items-center justify-center">
      <div class="text-center">
        <div class="mb-4 h-8 w-8 animate-spin rounded-full border-2 border-blue-500 border-t-transparent mx-auto"></div>
        <p class="text-sm text-slate-400">Loading heatmap data...</p>
      </div>
    </div>
  {:else}
    <div class="flex-1 relative" bind:this={container}>
      <!-- Grid: [dendro-spacer][col-dendro] / [row-dendro][heatmap] -->
      <div class="absolute inset-0 grid" style="grid-template-columns: {ROW_DENDRO_W}px 1fr; grid-template-rows: {COL_DENDRO_H}px 1fr;">
        <!-- Top-left spacer -->
        <div></div>
        <!-- Column dendrogram -->
        <div class="overflow-hidden">
          <svg bind:this={colDendroSvg}></svg>
        </div>
        <!-- Row dendrogram -->
        <div class="overflow-hidden">
          <svg bind:this={rowDendroSvg}></svg>
        </div>
        <!-- Heatmap canvas -->
        <div class="overflow-hidden">
          <canvas bind:this={canvas} class="block"></canvas>
        </div>
      </div>

      <!-- Title -->
      <div class="absolute top-1 left-1/2 -translate-x-1/2 text-xs text-slate-500 pointer-events-none">
        {heatmapData.nSamples} samples × {heatmapData.nAsvs} ASVs (Bray-Curtis, avg. linkage)
      </div>

      <!-- Tooltip -->
      {#if tooltip.show}
        <div class="fixed z-50 rounded bg-slate-800/95 px-3 py-1.5 text-xs text-slate-200 shadow-lg pointer-events-none"
          style="left: {tooltip.x + 12}px; top: {tooltip.y - 20}px;">
          {tooltip.text}
        </div>
      {/if}
    </div>
  {/if}
</div>
