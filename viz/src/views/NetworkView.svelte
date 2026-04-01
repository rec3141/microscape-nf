<script>
  import { onMount } from 'svelte';
  import createScatterplot from 'regl-scatterplot';
  import {
    store, countsByAsv,
    GROUP_COLORS, GROUP_HEX,
    buildTaxColorMap, getAsvColor, hexToRegl,
  } from '../stores/data.svelte.js';

  let { filters = {} } = $props();

  let canvasContainer = $state(null);
  let canvas = $state(null);
  let scatterplot = $state(null);
  let tooltip = $state({ show: false, x: 0, y: 0, text: '' });
  let hasZoomed = false;

  let taxonRe = $derived(() => {
    try { return filters.taxonFilter ? new RegExp(filters.taxonFilter, 'i') : null; }
    catch { return null; }
  });

  let filteredAsvs = $derived.by(() => {
    const re = taxonRe();
    const gf = filters.groupFlags || {};
    return store.asvs.filter(a => {
      if ((a.prevalence ?? 0) < (filters.minPrevalence || 0)) return false;
      const group = a.group ?? 'unknown';
      if (gf[group] === false) return false;
      if (re && !(re.test(a.taxonomy ?? '') || re.test(a.id ?? ''))) return false;
      return true;
    });
  });

  let idxMap = $derived.by(() => {
    const m = new Map();
    filteredAsvs.forEach((a, fi) => {
      const oi = store.asvs.indexOf(a);
      m.set(oi, fi);
    });
    return m;
  });

  let filteredEdges = $derived.by(() => {
    if (!filters.showEdges) return [];
    return (store.network?.edges ?? store.network ?? []).filter(e => {
      if (Math.abs(e.weight ?? 0) < (filters.corrThreshold || 0.3)) return false;
      return idxMap.has(e.source) && idxMap.has(e.target);
    });
  });

  let selectedAsvObj = $derived(
    store.selectedAsv != null ? store.asvs[store.selectedAsv] : null
  );

  onMount(() => {
    return () => { if (scatterplot) scatterplot.destroy(); };
  });

  $effect(() => {
    if (canvas && !scatterplot) {
      const rect = canvasContainer.getBoundingClientRect();
      const sp = createScatterplot({
        canvas,
        width: rect.width,
        height: rect.height,
        pointSize: 200,
        opacity: 0.85,
        lassoOnLongPress: true,
        backgroundColor: [0.02, 0.06, 0.1, 1],
      });

      sp.subscribe('pointover', (idx) => {
        const a = filteredAsvs[idx];
        if (a) {
          tooltip = { show: true, x: 0, y: 0,
            text: `${a.id ?? 'ASV'} | ${a.taxonomy ?? ''} | ${(a.total_reads ?? 0).toLocaleString()} reads`,
          };
        }
      });
      sp.subscribe('pointout', () => { tooltip = { show: false, x: 0, y: 0, text: '' }; });
      sp.subscribe('select', ({ points: indices }) => {
        if (indices.length > 0) {
          const oi = store.asvs.indexOf(filteredAsvs[indices[0]]);
          store.selectedAsv = oi >= 0 ? oi : null;
        }
      });
      scatterplot = sp;
    }
  });

  $effect(() => {
    if (!scatterplot || filteredAsvs.length === 0) return;
    const xArr = filteredAsvs.map(a => a.x ?? 0);
    const yArr = filteredAsvs.map(a => a.y ?? 0);
    const sizes = filteredAsvs.map(a => Math.max(25, Math.log2((a.total_reads ?? 1) + 1) * 12));
    const cmap = filters.colorByLevel !== 'group' ? buildTaxColorMap(filters.colorByLevel).colorMap : null;
    const perPointHex = filteredAsvs.map(a => {
      if (cmap) return getAsvColor(a.id, filters.colorByLevel, cmap);
      return GROUP_HEX[a.group ?? 'prokaryote'] ?? GROUP_HEX.unknown;
    });

    // Build palette + index mapping for regl-scatterplot
    const uniqueColors = [...new Set(perPointHex)];
    const colorIdx = {};
    uniqueColors.forEach((c, i) => { colorIdx[c] = i; });
    const zArr = perPointHex.map(c => colorIdx[c]);

    const maxSize = Math.max(...sizes, 1);
    const wArr = sizes.map(s => s / maxSize);

    scatterplot.set({
      pointColor: uniqueColors, colorBy: 'valueZ',
      pointSize: [2, 40], sizeBy: 'valueW',
    });
    scatterplot.draw({ x: xArr, y: yArr, z: zArr, w: wArr }).then(() => {
      if (!hasZoomed) {
        scatterplot.zoomToPoints(Array.from({ length: xArr.length }, (_, i) => i), {
          padding: 0.2, transition: true, transitionDuration: 500,
        });
        hasZoomed = true;
      }
    });
  });

  function handleResize() {
    if (scatterplot && canvasContainer) {
      const rect = canvasContainer.getBoundingClientRect();
      scatterplot.set({ width: rect.width, height: rect.height });
    }
  }
</script>

<svelte:window onresize={handleResize} />

<div class="flex h-full flex-col">
  <div class="relative flex-1" bind:this={canvasContainer}>
    <canvas bind:this={canvas} class="absolute inset-0 h-full w-full"></canvas>

    {#if tooltip.show}
      <div class="pointer-events-none absolute left-4 top-4 rounded bg-slate-800/90 px-3 py-1.5 text-xs text-slate-200 shadow-lg">
        {tooltip.text}
      </div>
    {/if}

    <div class="absolute bottom-4 left-4 text-xs text-slate-500">
      {filteredAsvs.length} / {store.asvs.length} ASVs | {filteredEdges.length} edges
    </div>
  </div>

  {#if selectedAsvObj}
    <div class="border-t border-slate-800 bg-slate-900/80 p-4">
      <div class="flex items-center justify-between">
        <h3 class="text-sm font-semibold text-slate-200">
          {selectedAsvObj.id ?? 'ASV'} &mdash; {selectedAsvObj.group ?? ''}
        </h3>
        <button class="text-xs text-slate-500 hover:text-slate-300" onclick={() => store.selectedAsv = null}>Close</button>
      </div>
      <p class="mt-1 text-xs text-slate-400">{selectedAsvObj.taxonomy ?? 'No taxonomy'}</p>
      <p class="text-xs text-slate-500">
        {(selectedAsvObj.total_reads ?? 0).toLocaleString()} total reads |
        Prevalence: {selectedAsvObj.prevalence ?? 0}
      </p>
    </div>
  {/if}
</div>
