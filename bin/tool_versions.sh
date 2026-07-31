#!/usr/bin/env bash
# Emit the ACTUAL version of every tool this pipeline runs, as JSON, by asking each
# binary rather than reading the environment spec.
#
# The spec says what was requested ("cutadapt", "r-base >=4.1"); only the binary knows
# what was resolved. A dataset analysed months later needs the second one — a claim that
# rests on DADA2's error model is not reproducible against "bioconductor-dada2".
#
# Never fails: a tool that is absent or does not answer is reported as null, because a
# manifest that is missing is worse than one with a hole in it.
set -uo pipefail

first_line() { head -n1 2>/dev/null | tr -d '\r' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//'; }
# Most tools print to stdout; a few (FastTree, mafft) use stderr.
probe() { { eval "$1"; } 2>&1 | first_line; }

# tool -> command, and a regex to pull the bare version out of a chatty banner
emit() {
    local name="$1" raw="$2" ver
    ver=$(printf '%s' "$raw" | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?([a-z0-9._-]*)?' | head -n1)
    if [ -z "$raw" ]; then
        printf '  "%s": null' "$name"
    else
        printf '  "%s": {"version": %s, "reported": "%s"}' "$name" \
            "$([ -n "$ver" ] && printf '"%s"' "$ver" || printf 'null')" \
            "$(printf '%s' "$raw" | sed 's/\\/\\\\/g; s/"/\\"/g' | cut -c1-160)"
    fi
}

printf '{\n'
sep=""
for spec in \
    "cutadapt|cutadapt --version" \
    "mafft|mafft --version" \
    "fasttree|FastTree -expert" \
    "vsearch|vsearch --version" \
    "python|python --version" \
    "R|R --version" \
    "nextflow|nextflow -v" \
    "dada2|Rscript -e 'cat(as.character(packageVersion(\"dada2\")))'" \
    "decipher|Rscript -e 'cat(as.character(packageVersion(\"DECIPHER\")))'" \
    "ape|Rscript -e 'cat(as.character(packageVersion(\"ape\")))'" \
    "biopython|python -c 'import Bio;print(Bio.__version__)'" \
    "scikit_learn|python -c 'import sklearn;print(sklearn.__version__)'" \
    "scipy|python -c 'import scipy;print(scipy.__version__)'" \
    "papa2|python -c 'import papa2;print(getattr(papa2,\"__version__\",\"unknown\"))'" \
    "microscape|python -c 'import microscape;print(getattr(microscape,\"__version__\",\"unknown\"))'" \
; do
    name="${spec%%|*}"; cmd="${spec#*|}"
    printf '%s' "$sep"; emit "$name" "$(probe "$cmd")"; sep=$',\n'
done
printf '\n}\n'
