<script>
  import { store, GROUP_HEX, buildTaxColorMap, getEffectiveColorLevel, findTaxonLevel } from '../stores/data.svelte.js';
  import AutocompleteInput from './AutocompleteInput.svelte';

  let { activeTab = 'samples', filters = $bindable({}) } = $props();

  // When user picks from autocomplete, navigate to that taxon's level
  function navigateToTaxon(name) {
    if (!name || name === 'unclassified' || filters.colorMode !== 'taxonomy') return;
    const level = findTaxonLevel(name);
    if (level) {
      const stack = [...(filters.navStack || [])];
      stack.push({ level: filters.colorByLevel, filter: filters.taxonFilter });
      filters.navStack = stack;
      filters.colorByLevel = level;
      filters.taxonFilter = name;
    }
  }

  // Collapsible sections
  let sections = $state({
    taxonomy: true,
    samples: true,
    network: true,
    phylogeny: true,
  });

  function toggle(section) {
    sections[section] = !sections[section];
  }

  // Taxonomy levels from data
  let taxLevels = $derived(
    store.taxonomy[Object.keys(store.taxonomy)[0]]?.levels || []
  );

  // Autocomplete candidates
  let taxCandidates = $derived.by(() => {
    const db = Object.keys(store.taxonomy)[0];
    if (!db || !store.taxonomy[db]?.assignments) return [];
    const terms = new Set();
    for (const asvId in store.taxonomy[db].assignments) {
      const vals = store.taxonomy[db].assignments[asvId];
      for (const v of vals) {
        if (v && v !== 'unclassified') terms.add(v);
      }
    }
    return [...terms].sort();
  });

  let sampleCandidates = $derived(
    store.samples.map(s => s.id).filter(Boolean).sort()
  );
</script>

<aside class="flex w-64 shrink-0 flex-col overflow-y-auto border-r border-slate-800 bg-slate-900/60">
  <div class="p-3 border-b border-slate-800">
    <p class="text-xs text-slate-500">{store.samples.length} samples | {store.asvs.length} ASVs</p>
  </div>

  <!-- ══ Taxonomy Filters (shared) ══ -->
  <div class="border-b border-slate-800">
    <button class="flex w-full items-center justify-between px-3 py-2 text-xs font-semibold uppercase tracking-wider text-slate-400 hover:text-slate-200" onclick={() => toggle('taxonomy')}>
      Taxonomy
      <span class="text-[10px]">{sections.taxonomy ? '▾' : '▸'}</span>
    </button>
    {#if sections.taxonomy}
      <div class="space-y-3 px-3 pb-3">
        <AutocompleteInput
          bind:value={filters.taxonFilter}
          label="Filter (regex)"
          placeholder="e.g. Proteobacteria"
          candidates={taxCandidates}
          onPick={navigateToTaxon}
        />

        <label class="block">
          <span class="text-xs text-slate-400">Color by</span>
          <select bind:value={filters.colorMode} class="mt-1 w-full rounded border border-slate-700 bg-slate-800 px-2 py-1 text-sm text-slate-200"
            onchange={() => { filters.navStack = []; filters.taxonFilter = ''; }}>
            <option value="taxonomy">Taxonomy</option>
            <option value="group">Group (broad)</option>
          </select>
        </label>

        {#if filters.colorMode === 'group'}
          <fieldset class="space-y-1">
            <legend class="text-xs text-slate-400">Groups</legend>
            {#each Object.keys(filters.groupFlags || {}) as group}
              <label class="flex items-center gap-2 text-sm capitalize">
                <input type="checkbox" bind:checked={filters.groupFlags[group]} class="accent-blue-500" />
                <span class="inline-block h-2.5 w-2.5 rounded-full" style="background:{GROUP_HEX[group]}"></span>
                {group}
              </label>
            {/each}
          </fieldset>
        {:else}
          {@const effectiveLevel = getEffectiveColorLevel(filters.colorByLevel, filters.taxonFilter)}
          {@const taxColors = buildTaxColorMap(effectiveLevel, filters.taxonFilter)}
          <div class="space-y-0.5 max-h-48 overflow-y-auto">
            <p class="text-[10px] text-slate-500 mb-1">
              {effectiveLevel === '_asv' ? 'ASV' : effectiveLevel} ({taxColors.ranked.length})
            </p>
            {#if filters.navStack?.length > 0}
              {@const prev = filters.navStack[filters.navStack.length - 1]}
              <button
                class="flex items-center gap-1.5 w-full text-left text-xs hover:bg-slate-800 rounded px-1 py-1 text-cyan-400 border-b border-slate-700 mb-1"
                onclick={() => {
                  const stack = [...filters.navStack];
                  const popped = stack.pop();
                  filters.navStack = stack;
                  filters.colorByLevel = popped.level;
                  filters.taxonFilter = popped.filter;
                }}
              >
                &#x25B4; Up{prev.filter ? ` to ${prev.filter}` : ` to ${prev.level}`}
              </button>
            {:else if filters.taxonFilter}
              <button
                class="flex items-center gap-1.5 w-full text-left text-xs hover:bg-slate-800 rounded px-1 py-1 text-cyan-400 border-b border-slate-700 mb-1"
                onclick={() => {
                  filters.taxonFilter = '';
                }}
              >
                &#x25B4; Up to {filters.colorByLevel}
              </button>
            {/if}
            {#each taxColors.ranked as item}
              <button
                class="flex items-center gap-1.5 w-full text-left text-xs hover:bg-slate-800 rounded px-1 py-0.5"
                onclick={() => {
                  if (effectiveLevel === '_asv') return;
                  if (item.name === 'unclassified') return;
                  const stack = [...(filters.navStack || [])];
                  stack.push({ level: filters.colorByLevel, filter: filters.taxonFilter });
                  filters.navStack = stack;
                  filters.colorByLevel = effectiveLevel;
                  filters.taxonFilter = item.name;
                }}
              >
                <span class="inline-block h-2 w-2 rounded-full shrink-0" style="background:{item.color}"></span>
                <span class="truncate">{item.name}</span>
                <span class="ml-auto text-slate-500 shrink-0">{item.count}</span>
              </button>
            {/each}
          </div>
        {/if}
      </div>
    {/if}
  </div>

  <!-- ══ Sample Controls ══ -->
  {#if activeTab === 'samples'}
    <div class="border-b border-slate-800">
      <button class="flex w-full items-center justify-between px-3 py-2 text-xs font-semibold uppercase tracking-wider text-slate-400 hover:text-slate-200" onclick={() => toggle('samples')}>
        Sample Filters
        <span class="text-[10px]">{sections.samples ? '▾' : '▸'}</span>
      </button>
      {#if sections.samples}
        <div class="space-y-3 px-3 pb-3">
          <label class="block">
            <span class="text-xs text-slate-400">Min reads: {(filters.minReads || 0).toLocaleString()}</span>
            <input type="range" min="0" max="50000" step="100" bind:value={filters.minReads} class="mt-1 w-full accent-blue-500" />
          </label>

          <AutocompleteInput
            bind:value={filters.sampleFilter}
            label="Sample filter (regex)"
            placeholder="e.g. Plate1"
            candidates={sampleCandidates}
          />

          <label class="block">
            <span class="text-xs text-slate-400">Point scale: {(filters.pointScale ?? 1).toFixed(1)}x</span>
            <input type="range" min="0.1" max="5" step="0.1" bind:value={filters.pointScale} class="mt-1 w-full accent-blue-500" />
          </label>
        </div>
      {/if}
    </div>
  {/if}

  <!-- ══ Network Controls ══ -->
  {#if activeTab === 'network'}
    <div class="border-b border-slate-800">
      <button class="flex w-full items-center justify-between px-3 py-2 text-xs font-semibold uppercase tracking-wider text-slate-400 hover:text-slate-200" onclick={() => toggle('network')}>
        Network Filters
        <span class="text-[10px]">{sections.network ? '▾' : '▸'}</span>
      </button>
      {#if sections.network}
        <div class="space-y-3 px-3 pb-3">
          <label class="block">
            <span class="text-xs text-slate-400">Min prevalence: {filters.minPrevalence || 0}</span>
            <input type="range" min="0" max="100" step="1" bind:value={filters.minPrevalence} class="mt-1 w-full accent-blue-500" />
          </label>

          <label class="block">
            <span class="text-xs text-slate-400">Correlation threshold: {(filters.corrThreshold || 0.3).toFixed(2)}</span>
            <input type="range" min="0" max="1" step="0.01" bind:value={filters.corrThreshold} class="mt-1 w-full accent-blue-500" />
          </label>

          <label class="flex items-center gap-2 text-sm">
            <input type="checkbox" bind:checked={filters.showEdges} class="accent-blue-500" />
            Show edges
          </label>
        </div>
      {/if}
    </div>
  {/if}

  <!-- ══ Phylogeny Controls ══ -->
  {#if activeTab === 'phylogeny'}
    <div class="border-b border-slate-800">
      <button class="flex w-full items-center justify-between px-3 py-2 text-xs font-semibold uppercase tracking-wider text-slate-400 hover:text-slate-200" onclick={() => toggle('phylogeny')}>
        Phylogeny
        <span class="text-[10px]">{sections.phylogeny ? '▾' : '▸'}</span>
      </button>
      {#if sections.phylogeny}
        <div class="space-y-3 px-3 pb-3">
          <label class="block">
            <span class="text-xs text-slate-400">Layout</span>
            <select bind:value={filters.treeLayout} class="mt-1 w-full rounded border border-slate-700 bg-slate-800 px-2 py-1 text-sm text-slate-200">
              <option value="rc">Rectangular</option>
              <option value="rd">Radial</option>
              <option value="cr">Circular</option>
              <option value="dg">Diagonal</option>
              <option value="hr">Hierarchical</option>
            </select>
          </label>

          <label class="block">
            <span class="text-xs text-slate-400">Tip labels</span>
            <select bind:value={filters.treeLabelLevel} class="mt-1 w-full rounded border border-slate-700 bg-slate-800 px-2 py-1 text-sm text-slate-200">
              <option value="id">ASV ID</option>
              {#each taxLevels as level}
                <option value={level}>{level}</option>
              {/each}
            </select>
          </label>

          <label class="block">
            <span class="text-xs text-slate-400">Min prevalence: {filters.treeMinPrevalence || 0}</span>
            <input type="range" min="0" max="100" step="1" bind:value={filters.treeMinPrevalence} class="mt-1 w-full accent-blue-500" />
          </label>
        </div>
      {/if}
    </div>
  {/if}

  <!-- ══ Selection info ══ -->
  <div class="mt-auto border-t border-slate-800 p-3">
    {#if store.selectedSample != null}
      <p class="text-xs text-slate-400">Selected: {store.samples[store.selectedSample]?.id || ''}</p>
    {:else if store.selectedAsv != null}
      <p class="text-xs text-slate-400">Selected: {store.asvs[store.selectedAsv]?.id || ''}</p>
    {:else}
      <p class="text-xs text-slate-500">Click to select, Shift+drag to lasso</p>
      <p class="text-xs text-slate-500">Shift+double-click to deselect</p>
    {/if}
  </div>
</aside>
