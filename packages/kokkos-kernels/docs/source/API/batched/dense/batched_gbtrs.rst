KokkosBatched::Gbtrs
####################

Defined in header: :code:`KokkosBatched_Gbtrs.hpp`

.. code:: c++

    template <typename ArgTrans, typename ArgAlgo>
    struct SerialGbtrs {
      template <typename AViewType, typename PivViewType, typename BViewType>
      KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A, const PivViewType &piv, const BViewType &b, const int kl,
                                               const int ku);
    };

Solves a system of the linear equations :math:`A \cdot X = B` or :math:`A^T \cdot X = B` or :math:`A^H \cdot X = B` with a general n-by-n band matrix :math:`A` using :math:`LU` factorization computed by ``Gbtrf``.
This operation is equivalent to the LAPACK routine ``SGBTRS`` (``CGBTRS``) or ``DGBTRS`` (``ZGBTRS``) for single or double precision for real (complex) matrix. 

Parameters
==========

:A: Input view containing the factors :math:`L` and :math:`U` from the factorization :math:`A = P \cdot L \cdot U` of the band matrix :math:`A` computed by ``Gbtrf``.  See `LAPACK reference <https://www.netlib.org/lapack/lug/node124.html>`_ for the band storage format.
:piv: The pivot indices computed by ``Gbtrf``.
:B: Input/output view containing the right-hand side on input and the solution on output.
:kl: The number of subdiagonals within the band of :math:`A`. kl >= 0
:ku: The number of superdiagonals within the band of :math:`A`. ku >= 0

Type Requirements
-----------------

- ``ArgTrans`` must be one of the following:
   - ``KokkosBatched::Trans::NoTranspose`` to solve a system :math:`A \cdot X = B`
   - ``KokkosBatched::Trans::Transpose`` to solve a system :math:`A^T \cdot X = B`
   - ``KokkosBatched::Trans::ConjTranspose`` to solve a system :math:`A^H \cdot X = B`
- ``ArgAlgo`` must be ``KokkosBatched::Algo::Gbtrs::Unblocked`` for the unblocked algorithm
- ``AViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 containing the general band matrix :math:`A`
- ``PivViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 containing the pivot indices
- ``BViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 containing the right-hand side that satisfies
  - ``std::is_same_v<typename BViewType::value_type, typename BViewType::non_const_value_type> == true``

Example
=======

.. literalinclude:: ../../../../../example/batched_solve/serial_gbtrs.cpp
  :language: c++

output:

.. code::

   gbtrf/gbtrs works correctly!
