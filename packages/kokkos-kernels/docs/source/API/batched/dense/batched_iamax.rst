KokkosBatched::Iamax
####################

Defined in header: :code:`KokkosBatched_Iamax.hpp`

.. code:: c++

    struct SerialIamax {
      template <typename XViewType>
      KOKKOS_INLINE_FUNCTION static typename XViewType::size_type invoke(const XViewType &x);
    };

Finds the index of the first element of :math:`x` having maximum absolute value. As well as Blas, this returns 0 for an empty vector.

Parameters
==========

:x: :math:`x` is a length n vector.

Type Requirements
-----------------

- ``XViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 containing a vector :math:`x`

Example
=======

.. literalinclude:: ../../../../../example/batched_solve/serial_iamax.cpp
  :language: c++

output:

.. code::

   iamax works correctly!
