SPARSE -- KokkosKernels sparse interfaces
=========================================

crsmatrix
---------
.. doxygenclass::    KokkosSparse::CrsMatrix
    :members:

ccsmatrix
---------
.. doxygenclass::    KokkosSparse::CcsMatrix
    :members:

crs2ccs
-------
.. doxygenfunction:: KokkosSparse::crs2ccs(OrdinalType nrows, OrdinalType ncols, SizeType nnz, ValViewType vals, RowMapViewType row_map, ColIdViewType col_ids)
.. doxygenfunction:: KokkosSparse::crs2ccs(KokkosSparse::CrsMatrix<ScalarType, OrdinalType, DeviceType, MemoryTraitsType, SizeType> &crsMatrix)

ccs2crs
-------
.. doxygenfunction:: KokkosSparse::ccs2crs(OrdinalType nrows, OrdinalType ncols, SizeType nnz, ValViewType vals, ColMapViewType col_map, RowIdViewType row_ids)
.. doxygenfunction:: KokkosSparse::ccs2crs(KokkosSparse::CcsMatrix<ScalarType, OrdinalType, DeviceType, MemoryTraitsType, SizeType> &ccsMatrix)

spmv
----

.. doxygenfunctions:: KokkosSparse::spmv(KokkosKernels::Experimental::Controls, const char[], const AlphaType&, const AMatrix&, const XVector&, const BetaType&, const YVector&)
.. doxygenfunctions:: KokkosSparse::spmv(KokkosKernels::Experimental::Controls controls, const char mode[], const AlphaType &alpha, const AMatrix &A, const XVector &x, const BetaType &beta, const YVector &y)
.. doxygenfunctions:: KokkosSparse::spmv(KokkosKernels::Experimental::Controls controls, const char mode[], const AlphaType &alpha, const AMatrix &A, const XVector &x, const BetaType &beta, const YVector &y, const RANK_ONE)
.. doxygenfunctions:: KokkosSparse::spmv(KokkosKernels::Experimental::Controls controls, const char mode[], const AlphaType &alpha, const AMatrix &A, const XVector &x, const BetaType &beta, const YVector &y, const RANK_TWO)
.. doxygenfunctions:: KokkosSparse::spmv(const char mode[], const AlphaType &alpha, const AMatrix &A, const XVector &x, const BetaType &beta, const YVector &y)


trsv
----
.. doxygenfunction:: KokkosSparse::trsv

spgemm
------
.. doxygenfunction:: KokkosSparse::spgemm

gauss
-----
.. doxygenfunction:: KokkosSparse::gauss
