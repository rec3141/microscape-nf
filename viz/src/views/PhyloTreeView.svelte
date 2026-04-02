<script>
  import PhylocanvasTree from '../components/PhylocanvasTree.svelte';
  import {
    store,
    GROUP_COLORS, GROUP_HEX,
    buildTaxColorMap, getAsvColor, getEffectiveColorLevel, hexToRgba255,
  } from '../stores/data.svelte.js';

  let { filters = {} } = $props();

  let treeType = $derived(filters.treeLayout || 'rc');

  // ---- Taxonomy data ----
  let primaryDb = $derived(Object.keys(store.taxonomy)[0] || null);
  let taxByAsv = $derived.by(() => {
    if (!primaryDb || !store.taxonomy[primaryDb]?.assignments) return {};
    return store.taxonomy[primaryDb].assignments;
  });
  let taxLevels = $derived(
    primaryDb && store.taxonomy[primaryDb]?.levels ? store.taxonomy[primaryDb].levels : []
  );

  // ---- Taxonomy regex for filtering ----
  let taxonRe = $derived(() => {
    try { return filters.taxonFilter ? new RegExp(filters.taxonFilter, 'i') : null; }
    catch { return null; }
  });

  // ---- Filtered ASV set (by taxonomy, group, prevalence) ----
  let filteredAsvIds = $derived.by(() => {
    const re = taxonRe();
    const gf = filters.groupFlags || {};
    const minPrev = filters.treeMinPrevalence || 0;
    const ids = new Set();
    for (const asv of store.asvs) {
      if ((asv.prevalence ?? 0) < minPrev) continue;
      const group = asv.group ?? 'unknown';
      if (gf[group] === false) continue;
      if (re && !(re.test(asv.taxonomy ?? '') || re.test(asv.id ?? ''))) continue;
      ids.add(asv.id);
    }
    return ids;
  });

  // ---- Node styles using shared color-by ----
  let effectiveColorLevel = $derived(
    filters.colorMode === 'group' ? 'group' : getEffectiveColorLevel(filters.colorByLevel, filters.taxonFilter)
  );

  let taxColorMap = $derived.by(() => {
    if (effectiveColorLevel === 'group') return null;
    return buildTaxColorMap(effectiveColorLevel, filters.taxonFilter);
  });

  // Bootstrap data
  let bootByAsv = $derived.by(() => {
    if (!primaryDb || !store.taxonomy[primaryDb]?.bootstraps) return {};
    return store.taxonomy[primaryDb].bootstraps;
  });

  // Build label for an ASV based on selected label levels
  function makeLabel(asvId) {
    const levels = filters.treeLabelLevels || ['id'];
    const tax = taxByAsv[asvId];
    const boot = bootByAsv[asvId];
    const parts = [];

    for (const level of levels) {
      if (level === 'id') {
        parts.push(asvId);
      } else if (level === 'bootstrap') {
        if (boot) {
          // Show bootstrap for the deepest assigned level
          for (let i = boot.length - 1; i >= 0; i--) {
            if (tax && tax[i]) { parts.push(`(${boot[i]})`); break; }
          }
        }
      } else if (tax) {
        const idx = taxLevels.indexOf(level);
        if (idx >= 0 && tax[idx]) parts.push(tax[idx]);
      }
    }

    return parts.length > 0 ? parts.join(' | ') : asvId;
  }

  let nodeStyles = $derived.by(() => {
    const styles = {};
    const cmap = taxColorMap?.colorMap;

    for (const asv of store.asvs) {
      const id = asv.id;
      const isFiltered = filteredAsvIds.has(id);
      const label = makeLabel(id);

      if (!isFiltered) {
        styles[id] = {
          fillColour: 'rgba(71,85,105,0.2)',
          fontColour: 'rgba(71,85,105,0.2)',
          shape: 'circle',
          nodeSize: 0.3,
          label,
        };
        continue;
      }

      let hex;
      if (cmap) {
        hex = getAsvColor(id, effectiveColorLevel, cmap);
      } else {
        hex = GROUP_HEX[asv.group] || GROUP_HEX.unknown;
      }

      styles[id] = {
        fillColour: hex,
        fontColour: hex,
        shape: 'circle',
        nodeSize: 1,
        label,
      };
    }
    return styles;
  });

  // ---- ASV lookup ----
  let asvById = $derived.by(() => {
    const map = {};
    for (const a of store.asvs) map[a.id] = a;
    return map;
  });

  // ---- Click handling ----
  let clickedNode = $state(null);

  function handleNodeClick(nodeId) {
    const asv = asvById[nodeId];
    if (asv) {
      const tax = taxByAsv[nodeId];
      clickedNode = { id: nodeId, asv, tax };
    } else {
      clickedNode = { id: nodeId, asv: null, tax: null };
    }
  }

  function closeInfoPanel() { clickedNode = null; }
</script>

<div class="flex h-full flex-col">

  <!-- Tree + info panel -->
  <div class="flex flex-1 overflow-hidden">
    {#if !store.treeNewick}
      <div class="flex-1 flex items-center justify-center">
        <div class="text-center">
          <p class="text-slate-400 mb-2">No phylogenetic tree available</p>
          <p class="text-xs text-slate-500">Run with <code class="bg-slate-800 px-1.5 py-0.5 rounded text-cyan-400">--run_phylogeny</code></p>
        </div>
      </div>
    {:else}
      <div class="{clickedNode ? 'w-2/3' : 'w-full'} transition-all">
        <PhylocanvasTree
          newick={store.treeNewick}
          {treeType}
          styles={nodeStyles}
          onNodeClick={handleNodeClick}
        />
      </div>

      {#if clickedNode}
        <div class="w-1/3 border-l border-slate-800 bg-slate-900/60 overflow-y-auto p-4">
          <div class="flex items-center justify-between mb-3">
            <h3 class="text-sm font-semibold text-cyan-400">{clickedNode.id}</h3>
            <button class="text-slate-500 hover:text-slate-300 text-lg" onclick={closeInfoPanel}>&times;</button>
          </div>

          {#if clickedNode.asv}
            <dl class="space-y-2 text-sm">
              <div class="flex justify-between">
                <dt class="text-slate-400">Group</dt>
                <dd>
                  <span class="inline-block h-2 w-2 rounded-full mr-1" style="background:{GROUP_HEX[clickedNode.asv.group] || GROUP_HEX.unknown}"></span>
                  {clickedNode.asv.group || 'unknown'}
                </dd>
              </div>
              <div class="flex justify-between">
                <dt class="text-slate-400">Total Reads</dt>
                <dd class="font-mono">{(clickedNode.asv.total_reads || 0).toLocaleString()}</dd>
              </div>
              <div class="flex justify-between">
                <dt class="text-slate-400">Prevalence</dt>
                <dd class="font-mono">{clickedNode.asv.prevalence || 0} samples</dd>
              </div>

              {#if clickedNode.tax}
                <div class="border-t border-slate-700 pt-2 mt-2">
                  <dt class="text-slate-400 mb-1 font-medium">Taxonomy</dt>
                  <dd class="font-mono text-xs">
                    {#each taxLevels as level, i}
                      {#if clickedNode.tax[i]}
                        <div><span class="text-slate-500">{level}:</span> {clickedNode.tax[i]}</div>
                      {/if}
                    {/each}
                  </dd>
                </div>
              {/if}
            </dl>
          {:else}
            <p class="text-xs text-slate-500">Internal node</p>
          {/if}
        </div>
      {/if}
    {/if}
  </div>
</div>
