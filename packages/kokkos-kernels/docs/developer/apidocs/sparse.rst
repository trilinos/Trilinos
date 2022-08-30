SPARSE -- KokkosKernels sparse interfaces
=========================================

crsmatrix
---------
.. doxygenclass::    KokkosSparse::CrsMatrix
    :members:

spmv
----
.. doxygenfunction:: KokkosSparse::spmv(KokkosKernels::Experimental::Controls, const char[], const AlphaType&, const AMatrix&, const XVector&, const BetaType&, const YVector&)
.. doxygenfunction:: KokkosSparse::spmv(KokkosKernels::Experimental::Controls, const char mode[], const AlphaType &alpha, const AMatrix &A, const XVector &x, const BetaType &beta, const YVector &y)
.. doxygenfunction:: KokkosSparse::spmv(KokkosKernels::Experimental::Controls controls, const char mode[], const AlphaType &alpha, const AMatrix &A, const XVector &x, const BetaType &beta, const YVector &y, const RANK_ONE)
.. doxygenfunction:: KokkosSparse::spmv(KokkosKernels::Experimental::Controls, const char mode[], const AlphaType &alpha, const AMatrix &A, const XVector &x, const BetaType &beta, const YVector &y, const RANK_TWO)
.. doxygenfunction:: KokkosSparse::spmv(const char mode[], const AlphaType &alpha, const AMatrix &A, const XVector &x, const BetaType &beta, const YVector &y)

trsv
----
.. doxygenfunction:: KokkosSparse::trsv

spgemm
------
.. doxygenfunction:: KokkosSparse::spgemm

gauss
-----
.. doxygenfunction:: KokkosSparse::gauss
