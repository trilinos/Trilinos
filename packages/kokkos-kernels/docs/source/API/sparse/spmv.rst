KokkosSparse::spmv
##################

Defined in header ``KokkosSparse_spmv.hpp``

.. code:: cppkokkos

  template <class ExecutionSpace, class Handle, class AlphaType, class AMatrix,
            class XVector, class BetaType, class YVector>
  void spmv(const ExecutionSpace& space, Handle* handle, const char mode[],
            const AlphaType& alpha, const AMatrix& A, const XVector& x,
            const BetaType& beta, const YVector& y);

  template <class ExecutionSpace, class AlphaType, class AMatrix, class XVector,
            class BetaType, class YVector,
	    typename = std::enable_if_t<Kokkos::is_execution_space_v<ExecutionSpace>>>
  void spmv(const ExecutionSpace& space, const char mode[], const AlphaType& alpha,
            const AMatrix& A, const XVector& x,
            const BetaType& beta, const YVector& y);

  template <class Handle, class AlphaType, class AMatrix, class XVector,
            class BetaType, class YVector,
            typename = std::enable_if_t<!Kokkos::is_execution_space<Handle>::value>>
  void spmv(Handle* handle, const char mode[], const AlphaType& alpha,
            const AMatrix& A, const XVector& x,
            const BetaType& beta, const YVector& y);

  template <class AlphaType, class AMatrix, class XVector, class BetaType, class YVector>
  void spmv(const char mode[], const AlphaType& alpha, const AMatrix& A,
            const XVector& x, const BetaType& beta, const YVector& y);

Kokkos sparse matrix-vector multiply. Computes y := alpha*Op(A)*x + beta*y, where Op(A) is controlled by mode (see below).

.. math::

   y = \beta*y + \alpha*A*x


1. Scale the ``y`` vector by ``beta``, and accumulate the result of the sparse matrix-vector product (``A*x``), scaled by ``alpha``, into ``y``.
2. Calls 1. using an spmv handle with ``SPMV_FAST_SETUP`` algorithm for the handle parameter.
3. Calls 1. using an instance of ``Handle::ExecutionSpaceType`` for the execution space parameter.
4. Calls 1. using an spmv handle with ``SPMV_FAST_SETUP`` algorithm and an instance of ``Handle::ExecutionSpaceType`` for the handle and execution space parameters respectively.

Parameters
==========

:space: execution space instance.

:handle: an spmv handle that stores multiple parameters for algorithm and third party libraries choices at run time.

:mode: mode to be applied to ``A``, possible values are "N" (normal), "T" (transpose), "C" (conjugate) and "H" (hermitian or conjugate-transpose).

:alpha, beta: scaling coefficents for the matrix ``A`` and left hand side ``y`` vector respectively.

:A: The matrix used to perform the matrix-vector product.

:x, y: The right and left hand side vectors used as input and output respectively.

Type Requirements
-----------------

- ``AMatrix`` must be either a :doc:`KokkosSparse::CrsMatrix <crs_matrix>` or a :doc:`KokkosSparse::BsrMatrix <bsr_matrix>` and have a memory space compatible with the ``ExecutionSpace`` type:

  - ``Kokkos::SpaceAccessibility<ExecutionSpace, typename AMatrix::memory_space>::accessible == true``

- `ExecutionSpace` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_

- ``XVector`` and ``YVector`` are two Kokkos Views of same rank, either rank 1 or rank 2 with memory spaces accessible from ``ExecutionSpace`` and ``YVector`` must store non-const data:

  - ``Kokkos::is_view_v<XVector> == true && Kokkos::is_view_v<YVector> == true``
  - ``XVector::rank() == YVector::rank() && (XVector::rank() == 1 || XVector::rank() == 2)``
  - ``Kokkos::SpaceAccessibility<ExecutionSpace, typename XVector::memory_space>::accessible``
  - ``Kokkos::SpaceAccessibility<ExecutionSpace, typename YVector::memory_space>::accessible``
  - ``!std::is_const_v<typename YVector::value_type> == true``

Example
=======

.. literalinclude:: ../../../../example/wiki/sparse/KokkosSparse_wiki_spmv.cpp
  :language: c++
  :lines: 16-

