KokkosBatched::Swap
###################

Defined in header: :code:`KokkosBatched_Swap.hpp`

.. code:: c++

    struct SerialSwap {
      template <typename XViewType, typename YViewType>
      KOKKOS_INLINE_FUNCTION static int invoke(const XViewType &x, const YViewType &y);
    };

    template <typename MemberType>
    struct TeamSwap {
      template <typename XViewType, typename YViewType>
      KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const XViewType &x, const YViewType &y);
    };

    template <typename MemberType>
    struct TeamVectorSwap {
      template <typename XViewType, typename YViewType>
      KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const XViewType &x, const YViewType &y);
    };

Swaps the elements of two vectors or matrices :math:`x` and :math:`y`.

- For a real vector or matrix :math:`x`, this operation is equivalent to the BLAS routine `SSWAP <https://www.netlib.org/blas/sswap.f>`_ or `DSWAP <https://www.netlib.org/blas/dswap.f>`_ for single or double precision.
- For a complex vector or matrix :math:`x`, this operation is equivalent to the BLAS routine `CSWAP <https://www.netlib.org/blas/cswap.f>`_ or `ZSWAP <https://www.netlib.org/blas/zswap.f>`_ for single or double precision.

Parameters
==========

:x: :math:`x` is a length :math:`n` vector or a :math:`m` by :math:`n` matrix.
:y: :math:`y` is a length :math:`n` vector or a :math:`m` by :math:`n` matrix.

Type Requirements
-----------------

- ``MemberType`` must be a Kokkos team member handle (only for ``TeamSwap`` and ``TeamVectorSwap``)

- ``XViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 or 2 containing a vector or matrix :math:`x` that satisfies ``std::is_same_v<typename XViewType::value_type, typename XViewType::non_const_value_type>``
- ``YViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 or 2 containing a vector or matrix :math:`y` that satisfies ``std::is_same_v<typename YViewType::value_type, typename YViewType::non_const_value_type>``
- The value type of ``XViewType`` and ``YViewType`` must be swappable.

.. note::

  This kernel supports both vector and matrix operations. When the input views :math:`X` and :math:`Y` are of rank 1, the kernel performs a vector operation (BLAS swap).
  When the input views :math:`X` and :math:`Y` are of rank 2, the kernel swaps matrices :math:`X` and :math:`Y`.

Example
=======

.. literalinclude:: ../../../../../example/batched_solve/serial_swap.cpp
  :language: c++
  :linenos:
  :lines: 4-

output:

.. code::

   swap works correctly!
