KokkosBatched::Laswp
####################

Defined in header: :code:`KokkosBatched_Laswp.hpp`

.. code:: c++

    template <typename ArgDirect>
    struct SerialLaswp {
      template <typename PivViewType, typename AViewType>
      KOKKOS_INLINE_FUNCTION static int invoke(const PivViewType &piv, const AViewType &A);
    };

Performs a series of row interchanges on the matrix :math:`A`. One row interchange is initiated for each of rows `K1` through `K2` of :math:`A`, where 
`piv` includes the elements from `K1` to `K2`.

Parameters
==========

:A: On entry, the matrix :math:`A` to which the row interchanges will be applied. On exit, the permuted matrix.
:piv: The vector of pivot indices.

Type Requirements
-----------------

- ``ArgDirect`` must be one of the following:
   - ``KokkosBatched::Direct::Forward`` to apply the pivot in the forward order
   - ``KokkosBatched::Direct::Backward`` to apply the pivot in the reverse order

- ``PivViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 containing the pivot indices
- ``AViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 containing the matrix :math:`A`

Example
=======

.. literalinclude:: ../../../../../example/batched_solve/serial_laswp.cpp
  :language: c++

output:

.. code::

   laswp works correctly!
