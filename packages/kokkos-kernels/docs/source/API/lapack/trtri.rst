KokkosLapack::trtri
###################

Defined in header: :code:`KokkosLapack_trtri.hpp`

.. code:: cppkokkos

  template <class AViewType>
  int trtri(const char uplo[], const char diag[], const AViewType& A);

Find the inverse of the triangular matrix A and store it in-place

.. math::

   A=A^{-1}

The function returns 0 upon success. Otherwise, it returns :math:`i` if the i-th diagonal element of :math:`A` is zero, and hence, :math:`A` is singular, and its inverse could not be completed.

Parameters
==========

:uplo: "U" or "u" indicates matrix :math:`A` is an upper triangular matrix. "L" or "l" indicates matrix :math:`A` is a lower triangular matrix.

:diag: "U" or "u" indicates the diagonal of :math:`A` is assumed to be unit. "N" or "n" indicates the diagonal of :math:`A` is assumed to be non-unit.

:A: Input matrix, as a 2-D Kokkos::View. On entry, :math:`A`, on successful exit, :math:`A^{-1}`.

Type Requirements
=================

No requirements

..
   .. note::

      We should provide an execution space instance overload of the function as well as implement some static assertion to check input is a view of rank 2 with accessible memory space...

Example
=======

TBD
