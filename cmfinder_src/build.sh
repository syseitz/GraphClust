#!/bin/bash
# Build CMfinder04 from bundled source.
# Works on Linux x86_64 (SSE) and Apple Silicon/ARM64 (dummy/non-SIMD fallback).
#
# Usage:
#   bash cmfinder_src/build.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

echo "Building CMfinder04 in $SCRIPT_DIR"

# Configure if needed
if [ ! -f Makefile ] || [ configure -nt Makefile ]; then
  echo "Running configure..."
  ./configure 2>&1 | tee configure.log
fi

# Ensure -DHMMER_THREADS is set in cmfinder04/Makefile
# (needed on all platforms — cmfinder.c guards esl_threads.h behind this define)
if ! grep -q "DHMMER_THREADS" cmfinder04/Makefile 2>/dev/null; then
  echo "Patching cmfinder04/Makefile to add -DHMMER_THREADS..."
  sed -i.bak 's|^DEFAULT_INCLUDES = \(.*\)|DEFAULT_INCLUDES = \1 -DHMMER_THREADS|' cmfinder04/Makefile
  rm -f cmfinder04/Makefile.bak
fi

# Build everything (Infernal libs + ViennaRNA 1.4 + cmfinder04)
echo "Building..."
make 2>&1 | tee build.log

# Copy binary to bin/
mkdir -p "$SCRIPT_DIR/bin"
cp "$SCRIPT_DIR/cmfinder04/cmfinder04" "$SCRIPT_DIR/bin/cmfinder"
chmod +x "$SCRIPT_DIR/bin/cmfinder"

echo ""
echo "CMfinder04 built successfully: $SCRIPT_DIR/bin/cmfinder"
echo ""
"$SCRIPT_DIR/bin/cmfinder" -h 2>&1 | head -3
