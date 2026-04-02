<script>
  import { onMount } from 'svelte';
  import Plotly from 'plotly.js-dist-min';
  import { store } from '../stores/data.svelte.js';

  let { filters = {} } = $props();

  let plotDiv;
  let hasPlot = false;
  let heatmapData = $state(null);

  // Load heatmap data
  onMount(async () => {
    try {
      const res = await fetch('/data/heatmap.json.gz');
      if (res.ok) {
        const buf = await res.arrayBuffer();
        const bytes = new Uint8Array(buf);
        let text;
        if (bytes[0] === 0x1f && bytes[1] === 0x8b) {
          const ds = new DecompressionStream('gzip');
          const reader = new Blob([buf]).stream().pipeThrough(ds).getReader();
          const chunks = [];
          while (true) {
            const { done, value } = await reader.read();
            if (done) break;
            chunks.push(value);
          }
          const combined = new Uint8Array(chunks.reduce((a, c) => a + c.length, 0));
          let offset = 0;
          for (const c of chunks) { combined.set(c, offset); offset += c.length; }
          text = new TextDecoder().decode(combined);
        } else {
          text = new TextDecoder().decode(buf);
        }
        heatmapData = JSON.parse(text);
      }
    } catch (e) {
      console.error('Failed to load heatmap data:', e);
    }

    return () => { if (plotDiv && hasPlot) Plotly.purge(plotDiv); };
  });

  $effect(() => {
    if (!plotDiv || !heatmapData) return;

    const { z, sampleIds, asvIds, asvLabels, rowDendro, colDendro } = heatmapData;

    if (!z || z.length < 2) return;

    // Scipy dendrogram coords use 5, 15, 25... for leaf positions
    // Scale to match plotly heatmap cell indices (0, 1, 2...)
    const nSamples = sampleIds.length;
    const nAsvs = asvIds.length;

    function scaleCoords(coords, n) {
      // Scipy uses 5, 15, 25... (step 10) for n leaves
      // We need 0, 1, 2... (step 1)
      return coords.map(arr => arr.map(v => (v - 5) / 10));
    }

    // Column dendrogram (above heatmap)
    const colDendroTraces = [];
    if (colDendro.icoord.length > 0) {
      const scaledI = scaleCoords(colDendro.icoord, nAsvs);
      for (let i = 0; i < colDendro.icoord.length; i++) {
        colDendroTraces.push({
          x: scaledI[i],
          y: colDendro.dcoord[i],
          type: 'scatter',
          mode: 'lines',
          line: { color: '#64748b', width: 0.8 },
          xaxis: 'x',
          yaxis: 'y3',
          hoverinfo: 'skip',
          showlegend: false,
        });
      }
    }

    // Row dendrogram (left of heatmap)
    const rowDendroTraces = [];
    if (rowDendro.icoord.length > 0) {
      const scaledI = scaleCoords(rowDendro.icoord, nSamples);
      for (let i = 0; i < rowDendro.icoord.length; i++) {
        rowDendroTraces.push({
          x: rowDendro.dcoord[i].map(v => -v),  // flip to left side
          y: scaledI[i],
          type: 'scatter',
          mode: 'lines',
          line: { color: '#64748b', width: 0.8 },
          xaxis: 'x2',
          yaxis: 'y',
          hoverinfo: 'skip',
          showlegend: false,
        });
      }
    }

    // Main heatmap
    const heatmapTrace = {
      z: z,
      x: asvLabels,
      y: sampleIds,
      type: 'heatmap',
      colorscale: 'Viridis',
      hovertemplate: 'Sample: %{y}<br>ASV: %{x}<br>Value: %{z:.4f}<extra></extra>',
      xaxis: 'x',
      yaxis: 'y',
      colorbar: {
        title: '∜(rel. abund.)',
        len: 0.4,
        thickness: 12,
        x: 1.02,
        titleside: 'right',
        tickfont: { size: 9, color: '#94a3b8' },
        titlefont: { size: 10, color: '#94a3b8' },
      },
    };

    const traces = [heatmapTrace, ...colDendroTraces, ...rowDendroTraces];

    const layout = {
      plot_bgcolor: 'rgba(2,6,15,1)',
      paper_bgcolor: 'rgba(2,6,15,1)',
      font: { color: '#94a3b8', size: 9 },
      margin: { l: 10, r: 60, t: 10, b: 10 },
      // Main heatmap axes
      xaxis: {
        domain: [0.12, 0.98],
        showticklabels: nAsvs <= 80,
        tickfont: { size: 7 },
        tickangle: -90,
      },
      yaxis: {
        domain: [0, 0.85],
        showticklabels: nSamples <= 80,
        tickfont: { size: 7 },
        autorange: true,
      },
      // Row dendrogram axis (left)
      xaxis2: {
        domain: [0, 0.11],
        showticklabels: false,
        showgrid: false,
        zeroline: false,
        autorange: 'reversed',
      },
      // Column dendrogram axis (top)
      yaxis3: {
        domain: [0.87, 1.0],
        showticklabels: false,
        showgrid: false,
        zeroline: false,
        anchor: 'x',
      },
      title: {
        text: `${nSamples} samples × ${nAsvs} ASVs (Bray-Curtis, avg. linkage)`,
        font: { size: 11, color: '#64748b' },
        x: 0.55, y: 0.995,
      },
    };

    const config = { scrollZoom: true, displayModeBar: false };

    if (!hasPlot) {
      Plotly.newPlot(plotDiv, traces, layout, config);
      hasPlot = true;
    } else {
      Plotly.react(plotDiv, traces, layout, config);
    }
  });
</script>

<div class="flex h-full flex-col">
  {#if !heatmapData}
    <div class="flex-1 flex items-center justify-center">
      <div class="text-center">
        <div class="mb-4 h-8 w-8 animate-spin rounded-full border-2 border-blue-500 border-t-transparent mx-auto"></div>
        <p class="text-sm text-slate-400">Computing heatmap clustering...</p>
      </div>
    </div>
  {:else}
    <div class="flex-1 relative">
      <div bind:this={plotDiv} class="absolute inset-0"></div>
    </div>
  {/if}
</div>
