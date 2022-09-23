Building Developer Documentation
================================

.. code-block::
    :caption: Installing dependencies on MacOS

        brew install doxygen
        pip install sphinx
        pip install breathe
        pip install sphinx-rtd-theme

.. code-block::
    :caption: How to build developer documentation

        cmake -DKokkosKernels_ENABLE_DOCS:BOOL=ON /path/to/kokkos-kernels
        make Doxygen
        make Sphinx
        open build/docs/docs/sphinx/index.html