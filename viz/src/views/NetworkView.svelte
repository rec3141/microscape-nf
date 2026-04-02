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

  let taxonRe = $derived(() => {
    try { return filters.taxonFilter ? new RegExp(filters.taxonFilter, 'i') : null; }
    catch { return null; }
  });

  let filteredAsvs = $derived.by(() => {
    const re = taxonRe();
    const gf = filters.groupFlags || {};
    return store.asvs.filter(a => {
      if ((a.n_samples ?? 0) < (filters.minPrevalence || 0)) return false;
      const group = a.group ?? 'unknown';
      if (gf[group] === false) return false;
      if (re && !(re.test(a.taxonomy ?? '') || re.test(a.id ?? ''))) return false;
      return true;
    });
  });

  let selectedAsvObj = $derived(
    store.selectedAsv != null ? store.asvs[store.selectedAsv] : null
  );

  $effect(() => {
    if (!plotDiv || filteredAsvs.length === 0) return;

    const colorLevel = filters.colorMode === 'group' ? 'group' : getEffectiveColorLevel(filters.colorByLevel, filters.taxonFilter);
    const cmap = colorLevel !== 'group' ? buildTaxColorMap(colorLevel, filters.taxonFilter).colorMap : null;

    const colors = filteredAsvs.map(a => {
      if (filters.colorMode === 'asvCluster') return getClusterColor(a.id, 'asvCluster', filters.asvClusterK);
      if (filters.colorMode === 'sampleCluster') return GROUP_HEX[a.group ?? 'prokaryote'] ?? GROUP_HEX.unknown;
      if (cmap) return getAsvColor(a.id, colorLevel, cmap);
      return GROUP_HEX[a.group ?? 'prokaryote'] ?? GROUP_HEX.unknown;
    });

    const trace = {
      x: filteredAsvs.map(a => a.x ?? 0),
      y: filteredAsvs.map(a => a.y ?? 0),
      mode: 'markers',
      type: 'scattergl',
      marker: {
        size: (() => {
          // If samples selected, size by abundance in those samples
          const selIds = store.selectedSample != null
            ? new Set([store.samples[store.selectedSample]?.id])
            : filters.lassoSampleIds?.size > 0 ? filters.lassoSampleIds : null;
          if (selIds) {
            const cMap = countsBySample();
            const asvCounts = new Map();
            for (const sid of selIds) {
              for (const e of (cMap.get(sid) ?? [])) {
                asvCounts.set(e.asv_idx, (asvCounts.get(e.asv_idx) || 0) + e.count);
              }
            }
            return filteredAsvs.map(a => {
              const idx = store.asvs.indexOf(a);
              const count = asvCounts.get(idx) || 0;
              const s = filters.networkPointScale ?? 1;
              return count > 0 ? Math.max(4, Math.log2(count + 1) * 2 * s) : 2;
            });
          }
          const s = filters.networkPointScale ?? 1;
          return filteredAsvs.map(a => Math.max(3, Math.log2((a.total_reads ?? 1) + 1) * 1.5 * s));
        })(),
        color: colors,
        opacity: 0.7,
      },
      text: filteredAsvs.map(a =>
        `${a.id}<br>${a.taxonomy ?? ''}<br>${(a.total_reads ?? 0).toLocaleString()} reads<br>${a.n_samples ?? 0} samples`
      ),
      hoverinfo: 'text',
      showlegend: false,
    };

    const savedZoom = filters.networkZoom;
    const layout = {
      dragmode: 'pan',
      xaxis: { title: '', zeroline: false, showgrid: false, showticklabels: false,
               ...(savedZoom ? { range: savedZoom.xRange } : {}) },
      yaxis: { title: '', zeroline: false, showgrid: false, showticklabels: false, scaleanchor: 'x',
               ...(savedZoom ? { range: savedZoom.yRange } : {}) },
      plot_bgcolor: 'rgba(2, 6, 15, 1)',
      paper_bgcolor: 'rgba(2, 6, 15, 1)',
      font: { color: '#94a3b8' },
      margin: { l: 20, r: 20, t: 10, b: 20 },
      title: { text: `${filteredAsvs.length} ASVs`, font: { size: 12, color: '#64748b' }, x: 0.01, y: 0.99 },
    };

    const config = { scrollZoom: true, displayModeBar: false, doubleClick: 'reset+autosize' };

    if (!hasPlot) {
      Plotly.newPlot(plotDiv, [trace], layout, config);
      hasPlot = true;

      plotDiv.on('plotly_click', (data) => {
        if (data.points?.[0]) {
          const idx = data.points[0].pointNumber;
          const oi = store.asvs.indexOf(filteredAsvs[idx]);
          store.selectedAsv = oi >= 0 ? oi : null;
        }
      });

      plotDiv.on('plotly_relayout', (update) => {
        if (update['xaxis.range[0]'] != null) {
          filters.networkZoom = {
            xRange: [update['xaxis.range[0]'], update['xaxis.range[1]']],
            yRange: [update['yaxis.range[0]'], update['yaxis.range[1]']],
          };
        }
      });

    } else {
      const zoom = filters.networkZoom;
      if (zoom) {
        layout.xaxis.range = zoom.xRange;
        layout.yaxis.range = zoom.yRange;
      } else {
        const curLayout = plotDiv.layout;
        if (curLayout?.xaxis?.range) layout.xaxis.range = curLayout.xaxis.range;
        if (curLayout?.yaxis?.range) layout.yaxis.range = curLayout.yaxis.range;
      }
      Plotly.react(plotDiv, [trace], layout, config);
    }
  });

  function handleKey(e) {
    if (!plotDiv || !hasPlot) return;
    if (e.key === 'Shift') {
      Plotly.relayout(plotDiv, { dragmode: e.type === 'keydown' ? 'lasso' : 'pan' });
    } else if (e.key === 'Escape' && e.type === 'keydown') {
      store.selectedAsv = null;
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
        Prevalence: {selectedAsvObj.n_samples ?? 0}
      </p>
    </div>
  {/if}
</div>
