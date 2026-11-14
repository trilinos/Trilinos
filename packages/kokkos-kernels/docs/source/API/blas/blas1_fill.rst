KokkosBlas::fill
################

Defined in header: :code:`KokkosBlas1_fill.hpp`

.. code:: c++

  template <class execution_space, class XMV>
  void fill(const execution_space& space, const XMV& X, const typename XMV::non_const_value_type& val);

  template <class XMV>
  void fill(const XMV& X, const typename XMV::non_const_value_type& val);

Fills `X` with the value of `val`.

1. Overwrites the entries of `X` with `val` using ``space`` instance.
2. Overwrites the entries of `X` with `val` using the default instance of ``typename XMV::execution_space``.

The function will throw a runtime exception if ``R.extent(0) != X.extent(0) || R.extent(1) != X.extent(1)``

Parameters
==========

:space: execution space instance

:X: input vector(s)

:val: the value used to overwrite the entries of `X`

Type Requirements
-----------------

- `execution_space` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_

- `XMV` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ that satisfies

  - ``Kokkos::SpaceAccessibility<execution_space, typename XMV::memory_space>::accessible``

Example
=======

.. literalinclude:: ../../../../example/docs/blas/KokkosBlas1_docs_fill.cpp
  :language: c++
  :linenos:
  :lines: 4-

output:

.. code::

   orignal values:
      x: {0, -1, -2, 3, -4}
   final values:
      x: {7, 7, 7, 7, 7}
