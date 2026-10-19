#!/usr/bin/env bash
# WSL -> Windows R. Windows binary packages must not be loaded into Linux R.
# All arguments after SCRIPT are filesystem paths and are converted by wslpath.
set -euo pipefail
command -v wslpath >/dev/null || { echo '[ERROR] This launcher requires WSL; native R can run the .R scripts directly.' >&2; exit 1; }
WIN_RSCRIPT="${WIN_RSCRIPT:-$(command -v Rscript.exe || true)}"
if [[ -z "$WIN_RSCRIPT" ]]; then
    shopt -s nullglob
    candidates=(/mnt/c/Program\ Files/R/R-*/bin/Rscript.exe)
    shopt -u nullglob
    if [[ "${#candidates[@]}" == 1 ]]; then
        WIN_RSCRIPT="${candidates[0]}"
    fi
fi
[[ -n "$WIN_RSCRIPT" && -f "$WIN_RSCRIPT" ]] || {
    echo '[ERROR] Set WIN_RSCRIPT to one Windows Rscript.exe (WSL path); no unique installation was detected.' >&2
    exit 1
}
WIN_R_LIB="${WIN_R_LIB:-}"
# Pass the Windows library path literally; /p would incorrectly translate it twice.
R_ENV="${WSLENV:+$WSLENV:}LC_ALL/w"
R_COMMAND=(env -u R_LIBS_USER)
if [[ -n "$WIN_R_LIB" ]]; then
    R_ENV="$R_ENV:R_LIBS_USER/w"
    R_COMMAND=(env "R_LIBS_USER=$WIN_R_LIB")
fi
R_COMMAND+=("WSLENV=$R_ENV" LC_ALL=C "$WIN_RSCRIPT" --vanilla)
if [[ "${1:-}" == --check ]]; then
    "${R_COMMAND[@]}" -e 'stopifnot(.Platform$OS.type == "windows"); print(R.version.string); print(.libPaths()); p <- c("ggplot2", "dplyr", "tidyr", "purrr", "scales", "patchwork"); ok <- vapply(p, requireNamespace, logical(1), quietly=TRUE); print(ok); stopifnot(all(ok), capabilities("cairo"), capabilities("png"))'
    exit
fi
[[ $# -ge 1 && -s "$1" ]] || { echo "Usage: bash rrmg_windows_r.sh SCRIPT [PATH ...] | --check" >&2; exit 2; }
args=()
for path in "$@"; do args+=("$(wslpath -w "$(realpath -m "$path")")"); done
"${R_COMMAND[@]}" "${args[@]}"
