Profiling
=========

Compile Times
-------------
1. Select a clang compiler
2. Configure and include `-ftime-trace` in your CXX FLAGS (this works with clang+cuda).
3. Clone and build https://github.com/aras-p/ClangBuildAnalyzer. Put the binary directory in your `PATH`.
4. Compile Kokkos and KokkosKernels
5. Create a directory called `ftime-trace-artifacts` in your build directory
6. Copy the json files you care about in this directory, for example:

.. code-block::

  cp ./{sparse,blas}/unit_test/CMakeFiles/*.dir/backends/*.json ftime-trace-artifacts/

7. Run `ClangBuildAnalyzer`:

.. code-block::

  ClangBuildAnalyzer --all ftime-trace-artifacts/ profile.txt
  ClangBuildAnalyzer --analyze profile.txt > analyze.txt

8. Open `analyze.txt`