<script>
  import { onMount } from 'svelte';
  import Plotly from 'plotly.js-dist-min';
  import {
    store, countsBySample,
    GROUP_HEX,
    buildTaxColorMap, getAsvColor, getEffectiveColorLevel, getClusterColor,
  } from '../stores/data.svelte.js';

  let { filters = {} } = $props();

  let plotDiv;
  let hasPlot = false;

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

  let cMap = $derived.by(() => countsBySample());

  let taxonRe = $derived(() => {
    try { return filters.taxonFilter ? new RegExp(filters.taxonFilter, 'i') : null; }
    catch { return null; }
  });

  let topTaxa = $derived.by(() => {
    const re = taxonRe();
    const gf = filters.groupFlags || {};

    // Merge counts: single click > lasso > all filtered
    const merged = new Map();
    let sampleIds;
    if (store.selectedSample != null) {
      sampleIds = [store.samples[store.selectedSample]?.id];
    } else if (filters.lassoSampleIds.size > 0) {
      sampleIds = [...filters.lassoSampleIds];
    } else {
      sampleIds = filteredSamples.map(s => s.id);
    }

    for (const sampleId of sampleIds) {
      for (const e of (cMap.get(sampleId) ?? [])) {
        const prev = merged.get(e.asv_idx) || 0;
        merged.set(e.asv_idx, prev + e.count);
      }
    }

    const total = [...merged.values()].reduce((s, c) => s + c, 0) || 1;

    return [...merged.entries()]
      .map(([asv_idx, count]) => ({
        asv: store.asvs[asv_idx],
        count,
        pct: ((count / total) * 100).toFixed(1),
      }))
      .filter(e => {
        if (!e.asv) return false;
        const group = e.asv.group ?? 'prokaryote';
        if (gf[group] === false) return false;
        if (re && !(re.test(e.asv.taxonomy ?? '') || re.test(e.asv.id ?? ''))) return false;
        return true;
      })
      .sort((a, b) => b.count - a.count)
      .slice(0, 50);
  });

  let selectedSampleObj = $derived(
    store.selectedSample != null ? store.samples[store.selectedSample] : null
  );

  // ── Build plotly traces ───────────────────────────────────────────────────

  $effect(() => {
    if (!plotDiv) return;
    if (filteredSamples.length === 0) {
      if (hasPlot) Plotly.react(plotDiv, [{ x: [], y: [], type: 'scattergl' }], {
        plot_bgcolor: 'rgba(2, 6, 15, 1)', paper_bgcolor: 'rgba(2, 6, 15, 1)',
        xaxis: { visible: false }, yaxis: { visible: false },
        annotations: [{ text: 'No samples match filters', showarrow: false,
          font: { color: '#64748b', size: 14 }, xref: 'paper', yref: 'paper', x: 0.5, y: 0.5 }],
      });
      return;
    }

    const colorLevel = filters.colorMode === 'group' ? 'group' : getEffectiveColorLevel(filters.colorByLevel, filters.taxonFilter);
    const cmap = colorLevel !== 'group' ? buildTaxColorMap(colorLevel, filters.taxonFilter).colorMap : null;
    const re = taxonRe();
    const gf = filters.groupFlags || {};
    const scale = filters.pointScale ?? 1;

    // Collect all points, then sort largest first so they draw behind small ones
    const allPoints = [];

    for (const sample of filteredSamples) {
      const entries = cMap.get(sample.id) ?? [];
      const totalCount = entries.reduce((s, e) => s + e.count, 0) || 1;

      for (const { asv_idx, count } of entries) {
        const asv = store.asvs[asv_idx];
        if (!asv) continue;

        const group = asv.group ?? 'prokaryote';
        if (gf[group] === false) continue;
        if ((asv.n_samples ?? 0) < (filters.minPrevalence || 0)) continue;
        if (re && !(re.test(asv.taxonomy ?? '') || re.test(asv.id ?? ''))) continue;

        const proportion = count / totalCount;
        let color;
        if (filters.colorMode === 'cluster') {
          color = getClusterColor(sample.id, 'sampleCluster', filters.sampleClusterK);
        } else if (cmap) {
          color = getAsvColor(asv.id, colorLevel, cmap);
        } else {
          color = GROUP_HEX[group] ?? GROUP_HEX.unknown;
        }

        allPoints.push({
          x: sample.x,
          y: sample.y,
          size: Math.max(2, Math.pow(proportion, 0.25) * 15 * scale),
          color,
          proportion,
          text: `${sample.id}<br>${(sample.total_reads ?? 0).toLocaleString()} reads | ${sample.n_asvs ?? 0} ASVs`,
        });
      }
    }

    // Sort largest first — plotly draws array order, so big points go first (behind)
    allPoints.sort((a, b) => b.proportion - a.proportion);

    const overlayTraces = [{
      x: allPoints.map(p => p.x),
      y: allPoints.map(p => p.y),
      mode: 'markers',
      type: 'scattergl',
      marker: {
        size: allPoints.map(p => p.size),
        color: allPoints.map(p => p.color),
        opacity: 0.7,
        line: { width: 0 },
      },
      text: allPoints.map(p => p.text),
      hoverinfo: 'text',
      showlegend: false,
    }];

    const savedZoom = filters.sampleZoom;
    const layout = {
      dragmode: 'pan',
      xaxis: { title: '', zeroline: false, showgrid: false, showticklabels: false,
               ...(savedZoom ? { range: savedZoom.xRange } : {}) },
      yaxis: { title: '', zeroline: false, showgrid: false, showticklabels: false, scaleanchor: 'x',
               ...(savedZoom ? { range: savedZoom.yRange } : {}) },
      plot_bgcolor: 'rgba(2, 6, 15, 1)',
      paper_bgcolor: 'rgba(2, 6, 15, 1)',
      font: { color: '#94a3b8' },
      legend: {
        bgcolor: 'rgba(15, 23, 42, 0.8)',
        font: { size: 10 },
      },
      margin: { l: 20, r: 20, t: 10, b: 20 },
      title: { text: `${filteredSamples.length} samples`, font: { size: 12, color: '#64748b' }, x: 0.01, y: 0.99 },
    };

    const config = { scrollZoom: true, displayModeBar: false, doubleClick: 'reset+autosize' };

    if (!hasPlot) {
      Plotly.newPlot(plotDiv, overlayTraces, layout, config);
      hasPlot = true;

      plotDiv.on('plotly_click', (data) => {
        if (data.points?.[0]) {
          const pt = data.points[0];
          let bestIdx = -1, bestDist = Infinity;
          filteredSamples.forEach((s, i) => {
            const d = (s.x - pt.x) ** 2 + (s.y - pt.y) ** 2;
            if (d < bestDist) { bestDist = d; bestIdx = i; }
          });
          const sId = bestIdx >= 0 ? filteredSamples[bestIdx]?.id : null;
          const sIdx = sId ? store.samples.findIndex(s => s.id === sId) : -1;
          store.selectedSample = sIdx >= 0 ? sIdx : null;
          filters.lassoSampleIds = new Set();
        }
      });

      plotDiv.on('plotly_selected', (data) => {
        if (data?.points?.length > 0) {
          const ids = new Set();
          for (const pt of data.points) {
            let bestIdx = -1, bestDist = Infinity;
            filteredSamples.forEach((s, i) => {
              const d = (s.x - pt.x) ** 2 + (s.y - pt.y) ** 2;
              if (d < bestDist) { bestDist = d; bestIdx = i; }
            });
            if (bestIdx >= 0) ids.add(filteredSamples[bestIdx].id);
          }
          filters.lassoSampleIds = ids;
          store.selectedSample = null;
        }
      });

      plotDiv.on('plotly_deselect', () => {
        filters.lassoSampleIds = new Set();
      });

      plotDiv.on('plotly_relayout', (update) => {
        if (update['xaxis.range[0]'] != null) {
          filters.sampleZoom = {
            xRange: [update['xaxis.range[0]'], update['xaxis.range[1]']],
            yRange: [update['yaxis.range[0]'], update['yaxis.range[1]']],
          };
        }
      });

    } else {
      // Restore persisted zoom or preserve current
      const zoom = filters.sampleZoom;
      if (zoom) {
        layout.xaxis.range = zoom.xRange;
        layout.yaxis.range = zoom.yRange;
      } else {
        const curLayout = plotDiv.layout;
        if (curLayout?.xaxis?.range) layout.xaxis.range = curLayout.xaxis.range;
        if (curLayout?.yaxis?.range) layout.yaxis.range = curLayout.yaxis.range;
      }
      Plotly.react(plotDiv, overlayTraces, layout, config);
    }
  });

  function handleKey(e) {
    if (!plotDiv || !hasPlot) return;
    if (e.key === 'Shift') {
      Plotly.relayout(plotDiv, { dragmode: e.type === 'keydown' ? 'lasso' : 'pan' });
    } else if (e.key === 'Escape' && e.type === 'keydown') {
      store.selectedSample = null;
      filters.lassoSampleIds = new Set();
      Plotly.restyle(plotDiv, { selectedpoints: [null] });
    }
  }

  onMount(() => {
    document.addEventListener('keydown', handleKey);
    document.addEventListener('keyup', handleKey);
    return () => {
      document.removeEventListener('keydown', handleKey);
      document.removeEventListener('keyup', handleKey);
      if (plotDiv && hasPlot) Plotly.purge(plotDiv);
    };
  });
</script>

<div class="flex h-full flex-col">
  <div class="flex-1 relative">
    <div bind:this={plotDiv} class="absolute inset-0"></div>
  </div>

  {#if topTaxa.length > 0}
    <div class="border-t border-slate-800 bg-slate-900/80 p-4">
      <div class="mb-2 flex items-center justify-between">
        <h3 class="text-sm font-semibold text-slate-200">
          {#if selectedSampleObj}
            {selectedSampleObj.id} &mdash; {(selectedSampleObj.total_reads ?? 0).toLocaleString()} reads
          {:else if filters.lassoSampleIds.size > 0}
            {filters.lassoSampleIds.size} selected samples &mdash; top {topTaxa.length} ASVs
          {:else}
            All samples ({filteredSamples.length}) &mdash; top {topTaxa.length} ASVs
          {/if}
        </h3>
        {#if selectedSampleObj}
          <button
            class="text-xs text-slate-500 hover:text-slate-300"
            onclick={() => store.selectedSample = null}
          >Clear selection</button>
        {/if}
      </div>

        <div class="max-h-48 overflow-y-auto">
          <table class="w-full text-xs">
            <thead class="sticky top-0 bg-slate-900 text-left text-slate-400">
              <tr>
                <th class="py-1 pr-4">ASV</th>
                <th class="py-1 pr-4">Taxonomy</th>
                <th class="py-1 pr-4 text-right">Reads</th>
                <th class="py-1 text-right">%</th>
              </tr>
            </thead>
            <tbody class="text-slate-300">
              {#each topTaxa as row}
                {@const colorLevel = filters.colorMode === 'group' ? 'group' : getEffectiveColorLevel(filters.colorByLevel, filters.taxonFilter)}
                {@const cmap = colorLevel !== 'group' ? buildTaxColorMap(colorLevel, filters.taxonFilter).colorMap : null}
                {@const rowColor = cmap ? getAsvColor(row.asv.id, colorLevel, cmap) : (GROUP_HEX[row.asv.group] ?? GROUP_HEX.prokaryote)}
                <tr class="border-t border-slate-800/50 hover:bg-slate-800/30">
                  <td class="py-1 pr-4 font-mono">
                    <span class="inline-block h-2.5 w-2.5 rounded-full mr-1.5" style="background:{rowColor}"></span>
                    {row.asv.id ?? ''}
                  </td>
                  <td class="py-1 pr-4 max-w-xs truncate">{row.asv.taxonomy ?? ''}</td>
                  <td class="py-1 pr-4 text-right font-mono">{row.count.toLocaleString()}</td>
                  <td class="py-1 text-right font-mono">{row.pct}</td>
                </tr>
              {/each}
            </tbody>
          </table>
        </div>
    </div>
  {/if}
</div>
