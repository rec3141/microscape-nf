<script>
  import { onMount } from 'svelte';
  import createScatterplot from 'regl-scatterplot';
  import {
    store, countsBySample,
    GROUP_COLORS, GROUP_HEX,
    buildTaxColorMap, getAsvColor,
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

  // Overlay: per-sample ASV entries placed at sample position, sorted large-to-small
  let overlayEntries = $derived.by(() => {
    if (!filters.showOverlay || filteredSamples.length === 0 || store.asvs.length === 0) return [];

    const re = taxonRe();
    const gf = filters.groupFlags || {};
    const entries = [];

    for (const sample of filteredSamples) {
      const sIdx = store.samples.indexOf(sample);
      const counts = cMap.get(sIdx) ?? [];
      const totalCount = counts.reduce((s, e) => s + e.count, 0) || 1;

      for (const { asv_idx, count } of counts) {
        const asv = store.asvs[asv_idx];
        if (!asv) continue;

        const group = asv.group ?? 'prokaryote';
        if (gf[group] === false) continue;
        if (re && !(re.test(asv.taxonomy ?? '') || re.test(asv.id ?? ''))) continue;

        entries.push({
          x: sample.x,
          y: sample.y,
          proportion: count / totalCount,
          sampleIdx: sIdx,
          asvIdx: asv_idx,
        });
      }
    }
    // Sort largest first so they draw behind smaller ones
    entries.sort((a, b) => b.proportion - a.proportion);
    return entries;
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
        let s, asvInfo;
        if (idx < numBasePts) {
          s = filteredSamples[idx];
        } else {
          const entry = overlayEntries[idx - numBasePts];
          if (entry) {
            s = filteredSamples.find((_, i) => store.samples.indexOf(filteredSamples[i]) === entry.sampleIdx)
              || store.samples[entry.sampleIdx];
            asvInfo = store.asvs[entry.asvIdx];
          }
        }
        if (s) {
          const text = asvInfo
            ? `${s.id} | ${asvInfo.id} ${asvInfo.taxonomy ?? ''}`
            : `${s.id ?? 'Sample'} | ${(s.total_reads ?? 0).toLocaleString()} reads`;
          tooltip = { show: true, x: 0, y: 0, text };
        }
      });

      sp.subscribe('pointout', () => {
        tooltip = { show: false, x: 0, y: 0, text: '' };
      });

      sp.subscribe('select', ({ points: indices }) => {
        if (indices.length > 0) {
          const idx = indices[0];
          let sIdx;
          if (idx < numBasePts) {
            sIdx = store.samples.indexOf(filteredSamples[idx]);
          } else {
            const entry = overlayEntries[idx - numBasePts];
            sIdx = entry ? entry.sampleIdx : -1;
          }
          store.selectedSample = sIdx >= 0 ? sIdx : null;
        }
      });

      scatterplot = sp;
    }
  });

  // Track how many base points for click routing
  let numBasePts = 0;

  $effect(() => {
    if (!scatterplot || filteredSamples.length === 0) return;

    // Read colorByLevel to ensure reactivity
    const colorLevel = filters.colorByLevel;
    const cmap = colorLevel !== 'group' ? buildTaxColorMap(colorLevel).colorMap : null;

    const xArr = [];
    const yArr = [];
    const sizes = [];
    const hexArr = [];

    // Layer 1: gray base points (background)
    for (const s of filteredSamples) {
      xArr.push(s.x);
      yArr.push(s.y);
      sizes.push(Math.max(30, Math.log2((s.total_reads ?? 1) + 1) * 14));
      hexArr.push('#445566');
    }

    // Layer 2: colored overlay (taxa at sample positions, large first)
    for (const entry of overlayEntries) {
      xArr.push(entry.x);
      yArr.push(entry.y);
      sizes.push(Math.max(3, Math.pow(entry.proportion, 0.25) * 50 * (filters.pointScale ?? 1)));

      const asv = store.asvs[entry.asvIdx];
      if (cmap && asv) {
        hexArr.push(getAsvColor(asv.id, colorLevel, cmap));
      } else if (asv) {
        hexArr.push(GROUP_HEX[asv.group ?? 'prokaryote'] ?? GROUP_HEX.unknown);
      } else {
        hexArr.push('#999999');
      }
    }

    numBasePts = filteredSamples.length;

    // Build palette + z indices
    const uniqueColors = [...new Set(hexArr)];
    const colorIdx = {};
    uniqueColors.forEach((c, i) => { colorIdx[c] = i; });
    const zArr = hexArr.map(c => colorIdx[c]);

    // Build size palette + w indices (categorical: each unique size is a bin)
    const uniqueSizes = [...new Set(sizes)].sort((a, b) => a - b);
    const sizeIdx = {};
    uniqueSizes.forEach((s, i) => { sizeIdx[s] = i; });
    const wArr = sizes.map(s => sizeIdx[s]);

    scatterplot.set({
      pointColor: uniqueColors, colorBy: 'valueZ',
      pointSize: uniqueSizes, sizeBy: 'valueW',
      opacity: 1.0,
    });
    scatterplot.draw({
      x: xArr,
      y: yArr,
      z: zArr,
      w: wArr,
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
