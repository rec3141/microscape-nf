<script>
  import { onMount } from 'svelte';
  import createScatterplot from 'regl-scatterplot';
  import {
    store, countsBySample,
    GROUP_COLORS, GROUP_HEX,
    buildTaxColorMap, getAsvColor, hexToRegl,
  } from '../stores/data.svelte.js';

  let { filters = {} } = $props();

  // ── Canvas + scatterplot ──────────────────────────────────────────────────
  let canvasContainer = $state(null);
  let canvas = $state(null);
  let scatterplot = $state(null);
  let tooltip = $state({ show: false, x: 0, y: 0, text: '' });
  let hasZoomed = false;

  // ── Derived data ──────────────────────────────────────────────────────────

  let filteredSamples = $derived.by(() => {
    let s = store.samples.filter(s => (s.total_reads ?? 0) >= (filters.minReads || 0));
    if (filters.sampleFilter) {
      try {
        const re = new RegExp(filters.sampleFilter, 'i');
        const matched = s.filter(sample => re.test(sample.id ?? ''));
        if (matched.length > 0) s = matched;
      } catch {}
    }
    return s;
  });

  let cMap = $derived(countsBySample());

  let taxonRe = $derived(() => {
    try { return filters.taxonFilter ? new RegExp(filters.taxonFilter, 'i') : null; }
    catch { return null; }
  });

  let taxColorMap = $derived.by(() => {
    if (filters.colorByLevel === 'group') return null;
    return buildTaxColorMap(filters.colorByLevel);
  });

  let overlayPoints = $derived.by(() => {
    if (!filters.showOverlay || filteredSamples.length === 0 || store.asvs.length === 0) return [];

    const re = taxonRe();
    const pts = [];
    const gf = filters.groupFlags || {};
    const cmap = taxColorMap?.colorMap;

    for (const sample of filteredSamples) {
      const sIdx = store.samples.indexOf(sample);
      const entries = cMap.get(sIdx) ?? [];
      const totalCount = entries.reduce((s, e) => s + e.count, 0) || 1;

      for (const { asv_idx, count } of entries) {
        const asv = store.asvs[asv_idx];
        if (!asv) continue;

        const group = asv.group ?? 'prokaryote';
        if (gf[group] === false) continue;
        if (re && !(re.test(asv.taxonomy ?? '') || re.test(asv.id ?? ''))) continue;

        const proportion = count / totalCount;
        let color;
        if (cmap) {
          color = hexToRegl(getAsvColor(asv.id, filters.colorByLevel, cmap));
        } else {
          color = GROUP_COLORS[group] ?? GROUP_COLORS.prokaryote;
        }

        pts.push({
          x: sample.x + (Math.random() - 0.5) * 0.3,
          y: sample.y + (Math.random() - 0.5) * 0.3,
          size: Math.max(1, proportion * 20),
          color,
          sampleIdx: sIdx,
          asvIdx: asv_idx,
        });
      }
    }
    return pts;
  });

  let topTaxa = $derived.by(() => {
    if (store.selectedSample == null) return [];
    const entries = cMap.get(store.selectedSample) ?? [];
    const total = entries.reduce((s, e) => s + e.count, 0) || 1;
    return entries
      .map(e => ({
        asv: store.asvs[e.asv_idx],
        count: e.count,
        pct: ((e.count / total) * 100).toFixed(1),
      }))
      .filter(e => e.asv)
      .sort((a, b) => b.count - a.count)
      .slice(0, 20);
  });

  let selectedSampleObj = $derived(
    store.selectedSample != null ? store.samples[store.selectedSample] : null
  );

  // ── Scatterplot lifecycle ─────────────────────────────────────────────────
  onMount(() => {
    return () => {
      if (scatterplot) scatterplot.destroy();
    };
  });

  $effect(() => {
    if (canvas && !scatterplot) {
      const rect = canvasContainer.getBoundingClientRect();
      const sp = createScatterplot({
        canvas,
        width: rect.width,
        height: rect.height,
        pointSize: 200,
        opacity: 0.8,
        lassoOnLongPress: true,
        backgroundColor: [0.02, 0.06, 0.1, 1],
      });

      sp.subscribe('pointover', (idx) => {
        const s = filteredSamples[idx];
        if (s) {
          tooltip = {
            show: true,
            x: 0, y: 0,
            text: `${s.id ?? 'Sample ' + idx} | ${(s.total_reads ?? 0).toLocaleString()} reads`,
          };
        }
      });

      sp.subscribe('pointout', () => {
        tooltip = { show: false, x: 0, y: 0, text: '' };
      });

      sp.subscribe('select', ({ points: indices }) => {
        if (indices.length > 0) {
          const sIdx = store.samples.indexOf(filteredSamples[indices[0]]);
          store.selectedSample = sIdx >= 0 ? sIdx : null;
        }
      });

      scatterplot = sp;
    }
  });

  $effect(() => {
    if (!scatterplot || filteredSamples.length === 0) return;

    const xArr = filteredSamples.map(s => s.x);
    const yArr = filteredSamples.map(s => s.y);
    const sizes = filteredSamples.map(s => Math.max(25, Math.log2((s.total_reads ?? 1) + 1) * 12));
    const colors = filteredSamples.map(() => [0.5, 0.55, 0.65, 0.5]);

    scatterplot.draw({
      x: xArr,
      y: yArr,
      size: sizes,
      color: colors,
    }).then(() => {
      if (!hasZoomed) {
        scatterplot.zoomToPoints(Array.from({ length: xArr.length }, (_, i) => i), {
          padding: 0.2,
          transition: true,
          transitionDuration: 500,
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
      {filteredSamples.length} / {store.samples.length} samples
    </div>
  </div>

  {#if selectedSampleObj}
    <div class="border-t border-slate-800 bg-slate-900/80 p-4">
      <div class="mb-2 flex items-center justify-between">
        <h3 class="text-sm font-semibold text-slate-200">
          {selectedSampleObj.id ?? 'Sample'} &mdash; {(selectedSampleObj.total_reads ?? 0).toLocaleString()} reads
        </h3>
        <button
          class="text-xs text-slate-500 hover:text-slate-300"
          onclick={() => store.selectedSample = null}
        >Close</button>
      </div>

      {#if topTaxa.length > 0}
        <div class="max-h-48 overflow-y-auto">
          <table class="w-full text-xs">
            <thead class="sticky top-0 bg-slate-900 text-left text-slate-400">
              <tr>
                <th class="py-1 pr-4">ASV</th>
                <th class="py-1 pr-4">Taxonomy</th>
                <th class="py-1 pr-4">Group</th>
                <th class="py-1 pr-4 text-right">Reads</th>
                <th class="py-1 text-right">%</th>
              </tr>
            </thead>
            <tbody class="text-slate-300">
              {#each topTaxa as row}
                <tr class="border-t border-slate-800/50 hover:bg-slate-800/30">
                  <td class="py-1 pr-4 font-mono">{row.asv.id ?? ''}</td>
                  <td class="py-1 pr-4 max-w-xs truncate">{row.asv.taxonomy ?? ''}</td>
                  <td class="py-1 pr-4">
                    <span class="inline-block h-2 w-2 rounded-full mr-1" style="background:{GROUP_HEX[row.asv.group] ?? GROUP_HEX.prokaryote}"></span>
                    {row.asv.group ?? ''}
                  </td>
                  <td class="py-1 pr-4 text-right font-mono">{row.count.toLocaleString()}</td>
                  <td class="py-1 text-right font-mono">{row.pct}</td>
                </tr>
              {/each}
            </tbody>
          </table>
        </div>
      {:else}
        <p class="text-xs text-slate-500">No taxa for this sample.</p>
      {/if}
    </div>
  {/if}
</div>
