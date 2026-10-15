#!/usr/bin/env bash
# WSL -> Windows R. Windows binary packages must not be loaded into Linux R.
# All arguments after SCRIPT are filesystem paths and are converted by wslpath.
set -euo pipefail
WIN_RSCRIPT="${WIN_RSCRIPT:-/mnt/c/Program Files/R/R-4.6.0/bin/Rscript.exe}"
WIN_R_LIB="${WIN_R_LIB:-C:/Users/q6321/AppData/Local/R/win-library/4.6}"
[[ -f "$WIN_RSCRIPT" ]] || { echo "[ERROR] Windows Rscript not found: $WIN_RSCRIPT" >&2; exit 1; }
command -v wslpath >/dev/null
# Pass the Windows library path literally; /p would incorrectly translate it twice.
R_ENV="${WSLENV:+$WSLENV:}R_LIBS_USER/w:LC_ALL/w"
if [[ "${1:-}" == --check ]]; then
    env WSLENV="$R_ENV" R_LIBS_USER="$WIN_R_LIB" LC_ALL=C "$WIN_RSCRIPT" --vanilla -e 'stopifnot(.Platform$OS.type == "windows", dir.exists(Sys.getenv("R_LIBS_USER"))); print(R.version.string); print(.libPaths()); p <- c("ggplot2", "dplyr", "tidyr", "purrr", "scales", "patchwork"); ok <- vapply(p, requireNamespace, logical(1), quietly=TRUE); print(ok); stopifnot(all(ok), capabilities("cairo"), capabilities("png"))'
    exit
fi
[[ $# -ge 1 && -s "$1" ]] || { echo "Usage: bash rrmg_windows_r.sh SCRIPT [PATH ...] | --check" >&2; exit 2; }
args=()
for path in "$@"; do args+=("$(wslpath -w "$(realpath -m "$path")")"); done
env WSLENV="$R_ENV" R_LIBS_USER="$WIN_R_LIB" LC_ALL=C "$WIN_RSCRIPT" --vanilla "${args[@]}"
