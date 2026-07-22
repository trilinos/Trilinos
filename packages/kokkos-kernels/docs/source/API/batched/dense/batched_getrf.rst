KokkosBatched::Getrf
####################

Defined in header: :code:`KokkosBatched_Getrf.hpp`

.. code:: c++

    template <typename ArgAlgo>
    struct SerialGetrf {
      template <typename AViewType, typename PivViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const AViewType &A, const PivViewType &piv);
    };

Computes a LU factorization of a general m-by-n matrix :math:`A` using partial pivoting with row interchanges. The factorization has the format
:math:`A = P \cdot L \cdot U`

where 

- :math:`P` is a permutation matrix
- :math:`L` is a lower triangular matrix with unit diagonal elements (lower trapezoidal if m > n)
- :math:`U` is an upper triangular matrix (upper trapezoidal if m < n)

This operation is equivalent to the LAPACK routine ``SGETRF`` (``CGETRF``) or ``DGETRF`` (``ZGETRF``) for single or double precision for real (complex) matrix. 
This is the recusive version of the algorithm. It divides the matrix :math:`A` into four submatrices:

.. math::

   A =
   \begin{bmatrix}
      A_{00} & A_{01} \\
      A_{10} & A_{11}
   \end{bmatrix}


where :math:`A_{00}` is a square matrix of size n0, :math:`A_{11}` is a matrix of size n1 by n1 with n0 = min(m, n) / 2 and n1 = n - n0.
This function calls itself to factorize

.. math::

   A_{0} =
   \begin{bmatrix}
      A_{00} \\
      A_{10}
   \end{bmatrix}

do the swaps on 

.. math::

   A_{1} =
   \begin{bmatrix}
      A_{01} \\
      A_{11}
   \end{bmatrix}

solve :math:`A_{01}`, update :math:`A_{11}`, then calls itself to factorize :math:`A_{11}` and do the swaps on :math:`A_{10}`.

.. note::

   On GPUs, the maximum matrix size of :math:`A` is limited to 4096 by 4096 as we rely on static stack to achive the recursion.

Parameters
==========

:A: On input, :math:`A` is a m by n general matrix. On output, the factors :math:`L` and :math:`U` from the factorization :math:`A = P \cdot L \cdot U`. The unit diagonal elements of :math:`L` are not stored.
:piv: On output, the pivot indices.

Type Requirements
-----------------

- ``ArgAlgo`` must be ``KokkosBatched::Algo::Getrf::Unblocked`` for the unblocked algorithm
- ``AViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 containing the general matrix :math:`A`
- ``PivViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 containing the pivot indices

Example
=======

.. literalinclude:: ../../../../../example/batched_solve/serial_getrs.cpp
  :language: c++

output:

.. code::

   getrf/getrs works correctly!
