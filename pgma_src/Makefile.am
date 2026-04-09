bin_PROGRAMS = pgma
pgma_SOURCES = pgma.cc pgma_aux.cc pgma_aux.hh

# On macOS, GNU libtool strips -Wl,-rpath from the linker command, so the
# OpenMP shared library (libomp.dylib) cannot be found at runtime.
# Re-embed the RPATH with install_name_tool after installation.
install-exec-hook:
	@case `uname` in \
	  Darwin) \
	    if test -n "@OMP_LIBDIR@"; then \
	      install_name_tool -add_rpath "@OMP_LIBDIR@" \
	        "$(DESTDIR)$(bindir)/pgma" 2>/dev/null || true; \
	      echo "Embedded RPATH @OMP_LIBDIR@ into pgma"; \
	    fi \
	    ;; \
	esac
