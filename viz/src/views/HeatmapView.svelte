<script>
  import { onMount } from 'svelte';
  import Plotly from 'plotly.js-dist-min';
  import {
    store, countsBySample,
    GROUP_HEX,
    buildTaxColorMap, getAsvColor, getEffectiveColorLevel,
  } from '../stores/data.svelte.js';

  let { filters = {} } = $props();

  let plotDiv;
  let hasPlot = false;

  // ── Clustering helpers ────────────────────────────────────────────────────

  function brayCurtis(a, b) {
    let sumMin = 0, sumAB = 0;
    for (let i = 0; i < a.length; i++) {
      sumMin += Math.min(a[i], b[i]);
      sumAB += a[i] + b[i];
    }
    return sumAB === 0 ? 0 : 1 - (2 * sumMin / sumAB);
  }

  function distMatrix(rows) {
    const n = rows.length;
    const D = Array.from({ length: n }, () => new Float64Array(n));
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        const d = brayCurtis(rows[i], rows[j]);
        D[i][j] = d;
        D[j][i] = d;
      }
    }
    return D;
  }

  // Average linkage clustering — returns {left, right, dist, size, id}
  function cluster(D) {
    const n = D.length;
    if (n <= 1) return { id: 0, size: 1 };

    const active = new Set(Array.from({ length: n }, (_, i) => i));
    const nodes = Array.from({ length: n }, (_, i) => ({ id: i, size: 1 }));
    const dist = Array.from({ length: 2 * n }, () => new Float64Array(2 * n));
    for (let i = 0; i < n; i++)
      for (let j = 0; j < n; j++) dist[i][j] = D[i][j];

    let nextId = n;
    while (active.size > 1) {
      let minD = Infinity, mi = -1, mj = -1;
      const arr = [...active];
      for (let ii = 0; ii < arr.length; ii++) {
        for (let jj = ii + 1; jj < arr.length; jj++) {
          if (dist[arr[ii]][arr[jj]] < minD) {
            minD = dist[arr[ii]][arr[jj]];
            mi = arr[ii]; mj = arr[jj];
          }
        }
      }

      const newNode = { left: nodes[mi], right: nodes[mj], dist: minD, size: nodes[mi].size + nodes[mj].size, id: nextId };
      nodes[nextId] = newNode;

      // Update distances (average linkage)
      for (const k of active) {
        if (k === mi || k === mj) continue;
        dist[nextId][k] = dist[k][nextId] =
          (dist[mi][k] * nodes[mi].size + dist[mj][k] * nodes[mj].size) / (nodes[mi].size + nodes[mj].size);
      }

      active.delete(mi);
      active.delete(mj);
      active.add(nextId);
      nextId++;
    }

    return nodes[nextId - 1];
  }

  // Get leaf order from dendrogram
  function leafOrder(node) {
    if (!node.left && !node.right) return [node.id];
    return [...leafOrder(node.left), ...leafOrder(node.right)];
  }

  // Build dendrogram trace (x or y coords)
  function dendrogramTrace(node, leafPositions, isRow) {
    const x = [], y = [];

    function walk(n) {
      if (!n.left && !n.right) {
        return { pos: leafPositions[n.id], height: 0 };
      }
      const L = walk(n.left);
      const R = walk(n.right);
      const h = n.dist;

      if (isRow) {
        // Row dendrogram: height on x-axis (negative, left of heatmap), position on y-axis
        x.push(-L.height, -h, -h, -R.height, null);
        y.push(L.pos, L.pos, R.pos, R.pos, null);
      } else {
        // Col dendrogram: position on x-axis, height on y-axis (above heatmap)
        x.push(L.pos, L.pos, R.pos, R.pos, null);
        y.push(L.height, h, h, R.height, null);
      }

      return { pos: (L.pos + R.pos) / 2, height: h };
    }

    walk(node);
    return { x, y };
  }

  // ── Build heatmap ─────────────────────────────────────────────────────────

  $effect(() => {
    if (!plotDiv || store.asvs.length === 0 || store.samples.length === 0) return;

    const cMap = countsBySample();
    const colorLevel = filters.colorMode === 'group' ? 'group' : getEffectiveColorLevel(filters.colorByLevel, filters.taxonFilter);
    const cmap = colorLevel !== 'group' ? buildTaxColorMap(colorLevel, filters.taxonFilter).colorMap : null;

    // Filter samples and ASVs
    const sampleIds = store.samples
      .filter(s => (s.total_reads ?? 0) >= (filters.minReads || 0))
      .map(s => s.id);
    const asvIds = store.asvs
      .filter(a => (a.n_samples ?? 0) >= (filters.minPrevalence || 0))
      .map(a => a.id);

    if (sampleIds.length < 2 || asvIds.length < 2) {
      if (hasPlot) Plotly.react(plotDiv, [], {
        plot_bgcolor: 'rgba(2,6,15,1)', paper_bgcolor: 'rgba(2,6,15,1)',
        annotations: [{ text: 'Need ≥2 samples and ≥2 ASVs', showarrow: false,
          font: { color: '#64748b', size: 14 }, xref: 'paper', yref: 'paper', x: 0.5, y: 0.5 }],
      });
      return;
    }

    // Limit for performance
    const maxSamples = Math.min(sampleIds.length, 200);
    const maxAsvs = Math.min(asvIds.length, 200);
    const useSamples = sampleIds.slice(0, maxSamples);
    const useAsvs = asvIds.slice(0, maxAsvs);

    // Build abundance matrix: 4th-root of relative abundance
    const asvIdxMap = {};
    store.asvs.forEach((a, i) => { asvIdxMap[a.id] = i; });

    const matrix = [];
    for (const sid of useSamples) {
      const row = new Float64Array(useAsvs.length);
      const entries = cMap.get(sid) ?? [];
      const total = entries.reduce((s, e) => s + e.count, 0) || 1;
      for (const e of entries) {
        const asv = store.asvs[e.asv_idx];
        if (!asv) continue;
        const asvPos = useAsvs.indexOf(asv.id);
        if (asvPos >= 0) {
          row[asvPos] = Math.pow(e.count / total, 0.25);
        }
      }
      matrix.push(row);
    }

    // Cluster samples (rows)
    const rowDist = distMatrix(matrix);
    const rowTree = cluster(rowDist);
    const rowOrder = leafOrder(rowTree);

    // Cluster ASVs (columns) — transpose first
    const colMatrix = [];
    for (let j = 0; j < useAsvs.length; j++) {
      const col = new Float64Array(useSamples.length);
      for (let i = 0; i < useSamples.length; i++) col[i] = matrix[i][j];
      colMatrix.push(col);
    }
    const colDist = distMatrix(colMatrix);
    const colTree = cluster(colDist);
    const colOrder = leafOrder(colTree);

    // Reorder matrix
    const orderedZ = rowOrder.map(ri => colOrder.map(ci => matrix[ri][ci]));
    const orderedSampleIds = rowOrder.map(i => useSamples[i]);
    const orderedAsvIds = colOrder.map(i => useAsvs[i]);

    // Side colors
    const rowColors = orderedSampleIds.map(() => '#475569');  // could color by metadata later

    const colColors = orderedAsvIds.map(id => {
      if (cmap) return getAsvColor(id, colorLevel, cmap);
      const asv = store.asvs.find(a => a.id === id);
      return GROUP_HEX[asv?.group ?? 'prokaryote'] ?? GROUP_HEX.unknown;
    });

    // Dendrograms
    const rowLeafPos = {};
    rowOrder.forEach((id, i) => { rowLeafPos[id] = i; });
    const colLeafPos = {};
    colOrder.forEach((id, i) => { colLeafPos[id] = i; });

    const rowDendro = dendrogramTrace(rowTree, rowLeafPos, true);
    const colDendro = dendrogramTrace(colTree, colLeafPos, false);

    // ASV labels (short taxonomy)
    const db = Object.keys(store.taxonomy)[0];
    const assignments = db ? store.taxonomy[db]?.assignments : {};
    const asvLabels = orderedAsvIds.map(id => {
      const tax = assignments?.[id];
      if (tax) {
        for (let i = tax.length - 1; i >= 0; i--) {
          if (tax[i] && tax[i] !== 'unclassified') return `${tax[i]}`;
        }
      }
      return id;
    });

    const traces = [
      // Heatmap
      {
        z: orderedZ,
        x: asvLabels,
        y: orderedSampleIds,
        type: 'heatmap',
        colorscale: 'Viridis',
        hovertemplate: 'Sample: %{y}<br>ASV: %{x}<br>Value: %{z:.4f}<extra></extra>',
        xaxis: 'x',
        yaxis: 'y',
        colorbar: { title: '∜(rel. abund.)', len: 0.5, thickness: 12, titleside: 'right',
                    tickfont: { size: 9, color: '#94a3b8' }, titlefont: { size: 10, color: '#94a3b8' } },
      },
      // Column-side color bar
      {
        z: [colColors.map(() => 1)],
        x: asvLabels,
        y: ['group'],
        type: 'heatmap',
        colorscale: [[0, colColors[0] || '#475569'], [1, colColors[0] || '#475569']],
        showscale: false,
        xaxis: 'x',
        yaxis: 'y2',
        hoverinfo: 'skip',
        // Use marker colors directly
        zmin: 0, zmax: 1,
      },
      // Column dendrogram
      {
        x: colDendro.x,
        y: colDendro.y,
        type: 'scatter',
        mode: 'lines',
        line: { color: '#64748b', width: 0.5 },
        xaxis: 'x',
        yaxis: 'y3',
        hoverinfo: 'skip',
        showlegend: false,
      },
    ];

    // Replace the colorscale hack with actual per-cell colors for the side bar
    // plotly doesn't support per-cell colors on heatmap, so use a scatter instead
    traces[1] = {
      x: asvLabels,
      y: asvLabels.map(() => 'group'),
      type: 'scatter',
      mode: 'markers',
      marker: { color: colColors, size: 8, symbol: 'square' },
      xaxis: 'x',
      yaxis: 'y2',
      hoverinfo: 'skip',
      showlegend: false,
    };

    const layout = {
      plot_bgcolor: 'rgba(2,6,15,1)',
      paper_bgcolor: 'rgba(2,6,15,1)',
      font: { color: '#94a3b8', size: 9 },
      margin: { l: 150, r: 30, t: 60, b: 10 },
      xaxis: { showticklabels: false, domain: [0, 1] },
      yaxis: { domain: [0, 0.85], autorange: true, tickfont: { size: 8 } },
      yaxis2: { domain: [0.86, 0.88], showticklabels: false },
      yaxis3: { domain: [0.89, 1.0], showticklabels: false, showgrid: false, zeroline: false },
      title: { text: `${orderedSampleIds.length} samples × ${orderedAsvIds.length} ASVs`,
               font: { size: 12, color: '#64748b' }, x: 0.5, y: 0.99 },
    };

    const config = { scrollZoom: true, displayModeBar: false };

    if (!hasPlot) {
      Plotly.newPlot(plotDiv, traces, layout, config);
      hasPlot = true;
    } else {
      Plotly.react(plotDiv, traces, layout, config);
    }
  });

  onMount(() => {
    return () => { if (plotDiv && hasPlot) Plotly.purge(plotDiv); };
  });
</script>

<div class="flex h-full flex-col">
  <div class="flex-1 relative">
    <div bind:this={plotDiv} class="absolute inset-0"></div>
  </div>
</div>
