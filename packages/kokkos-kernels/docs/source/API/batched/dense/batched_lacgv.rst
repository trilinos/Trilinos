KokkosBatched::Lacgv
####################

Defined in header: :code:`KokkosBatched_Lacgv.hpp`

.. code:: c++

    struct SerialLacgv {
      template <typename XViewType>
      KOKKOS_INLINE_FUNCTION static int invoke(const XViewType &x);
    };

Conjugates the elements of a complex vector :math:`x`. No operation is performed if :math:`x` is real.

Parameters
==========

:x: :math:`x` is a length n vector.

Type Requirements
-----------------

- ``XViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 representing the vector :math:`x`

Example
=======

.. literalinclude:: ../../../../../example/batched_solve/serial_lacgv.cpp
  :language: c++

output:

.. code::

   lacgv works correctly!
