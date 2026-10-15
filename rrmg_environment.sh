#!/usr/bin/env bash
# Source after the project config. Environment activation belongs to the caller.
# Empty overrides mean use PATH; no Conda installation or site path is assumed.
if [[ -n "${BBMAP_DIR:-}" ]]; then
  export PATH="$BBMAP_DIR:$PATH"
fi
if [[ -n "${JAVA_BIN:-}" ]]; then
  JAVA_BIN="$(command -v "$JAVA_BIN")" || { echo "[ERROR] Invalid JAVA_BIN" >&2; return 1; }
  [[ -x "$JAVA_BIN" && "$(basename "$JAVA_BIN")" == java ]] || {
    echo "[ERROR] JAVA_BIN must be an executable named java: $JAVA_BIN" >&2; return 1;
  }
  export PATH="$(dirname "$JAVA_BIN"):$PATH"
elif [[ -n "${JAVA_HOME:-}" ]]; then
  [[ -x "$JAVA_HOME/bin/java" ]] || { echo "[ERROR] Invalid JAVA_HOME: $JAVA_HOME" >&2; return 1; }
  JAVA_BIN="$JAVA_HOME/bin/java"
  export PATH="$JAVA_HOME/bin:$PATH"
else
  JAVA_BIN="$(command -v java || true)"
fi
if [[ -n "$JAVA_BIN" ]]; then
  # Conda bin/java may be a symlink into lib/jvm/bin/java.
  JAVA_HOME="$(dirname "$(dirname "$(readlink -f "$JAVA_BIN")")")"
  export JAVA_BIN JAVA_HOME
fi
RANDOMREADSMG_SH="${RANDOMREADSMG_SH:-$(command -v randomreadsmg.sh || true)}"
if [[ -n "$RANDOMREADSMG_SH" ]]; then
  RANDOMREADSMG_SH="$(readlink -f "$RANDOMREADSMG_SH")"
  export RANDOMREADSMG_SH
fi
RANDOMREADS_SH="${RANDOMREADS_SH:-$(command -v randomreads.sh || true)}"
if [[ -z "$RANDOMREADS_SH" && -n "$RANDOMREADSMG_SH" && -s "$(dirname "$RANDOMREADSMG_SH")/randomreads.sh" ]]; then
  RANDOMREADS_SH="$(dirname "$RANDOMREADSMG_SH")/randomreads.sh"
fi
if [[ -n "$RANDOMREADS_SH" ]]; then
  RANDOMREADS_SH="$(readlink -f "$RANDOMREADS_SH")"
  export RANDOMREADS_SH
fi
