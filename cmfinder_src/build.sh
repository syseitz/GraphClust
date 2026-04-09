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

# On ARM (Apple Silicon): patch the generated cmfinder04/Makefile to add -DHMMER_THREADS
# (the bundled Infernal 1.1 builds with impl_dummy on ARM, but cmfinder04 needs this define)
ARCH=$(uname -m)
if [ "$ARCH" = "arm64" ] || [ "$ARCH" = "aarch64" ]; then
  echo "ARM detected — ensuring -DHMMER_THREADS is set for cmfinder04..."
  if ! grep -q "DHMMER_THREADS" cmfinder04/Makefile 2>/dev/null; then
    sed -i.bak 's|^DEFAULT_INCLUDES = \(.*\)|DEFAULT_INCLUDES = \1 -DHMMER_THREADS|' cmfinder04/Makefile
    rm -f cmfinder04/Makefile.bak
  fi
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
