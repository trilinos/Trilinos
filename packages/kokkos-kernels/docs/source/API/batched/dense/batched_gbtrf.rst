KokkosBatched::Gbtrf
####################

Defined in header: :code:`KokkosBatched_Gbtrf.hpp`

.. code:: c++

    template <typename ArgAlgo>
    struct SerialGbtrf {
      template <typename ABViewType, typename PivViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ABViewType &Ab, const PivViewType &piv, const int kl, const int ku,
             const int m = -1);
    };

Computes a LU factorization of a general m-by-n band matrix :math:`A` using partial pivoting with row interchanges. The factorization has the format
:math:`A = P \cdot L \cdot U`. 

where 

- :math:`P` is a permutation matrix
- :math:`L` is a lower triangular matrix with unit diagonal elements
- :math:`U` is an upper triangular matrix

Parameters
==========

:Ab: On input, :math:`Ab` is a m by n general band matrix :math:`A` in the band storage. On output, the factors :math:`L` and :math:`U` from the factorization :math:`A = P \cdot L \cdot U`. The unit diagonal elements of :math:`L` are not stored. See `LAPACK reference <https://www.netlib.org/lapack/lug/node124.html>`_ for the band storage format.
:piv: On output, the pivot indices.
:kl: The number of subdiagonals within the band of :math:`A`. kl >= 0
:ku: The number of superdiagonals within the band of :math:`A`. ku >= 0
:m: The number of rows of the matrix :math:`A`. (optional, default is -1, corresponding to m == n)

Type Requirements
-----------------

- ``ArgAlgo`` must be ``KokkosBatched::Algo::Gbtrf::Unblocked`` for the unblocked algorithm
- ``ABViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 containing the general band matrix :math:`A`
- ``PivViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 containing the pivot indices

Example
=======

.. literalinclude:: ../../../../../example/batched_solve/serial_gbtrs.cpp
  :language: c++

output:

.. code::

   gbtrf/gbtrs works correctly!
