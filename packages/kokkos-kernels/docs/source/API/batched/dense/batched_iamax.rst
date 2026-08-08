KokkosBatched::Iamax
####################

Defined in header: :code:`KokkosBatched_Iamax.hpp`

.. code:: c++

    struct SerialIamax {
      template <typename XViewType>
      KOKKOS_INLINE_FUNCTION static typename XViewType::size_type invoke(const XViewType &x);
    };

    template <typename MemberType>
    struct TeamIamax {
      template <typename XViewType>
      KOKKOS_INLINE_FUNCTION static typename XViewType::size_type invoke(const MemberType &member, const XViewType &x);
    };

    template <typename MemberType>
    struct TeamVectorIamax {
      template <typename XViewType>
      KOKKOS_INLINE_FUNCTION static typename XViewType::size_type invoke(const MemberType &member, const XViewType &x);
    };

Finds the index of the first element of :math:`x` having maximum absolute value. As well as Blas, this returns 0 for an empty vector.

- For a real vector :math:`x`, this operation is equivalent to the BLAS routine `ISAMAX <https://www.netlib.org/blas/isamax.f>`_ or `IDAMAX <https://www.netlib.org/blas/idamax.f>`_ for single or double precision.
- For a complex vector :math:`x`, this operation is equivalent to the BLAS routine `ICAMAX <https://www.netlib.org/blas/icamax.f>`_ or `IZAMAX <https://www.netlib.org/blas/izamax.f>`_ for single or double precision.

Parameters
==========

:x: :math:`x` is a length n vector.

Type Requirements
-----------------

- ``MemberType`` must be a Kokkos team member handle (only for ``TeamIamax`` and ``TeamVectorIamax``)

- ``XViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 containing a vector :math:`x`

Example
=======

.. literalinclude:: ../../../../../example/batched_solve/serial_iamax.cpp
  :language: c++
  :linenos:
  :lines: 4-

output:

.. code::

   iamax works correctly!
