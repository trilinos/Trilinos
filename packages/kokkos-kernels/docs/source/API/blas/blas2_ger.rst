KokkosBlas::ger
###############

Defined in header: :code:`KokkosBlas2_ger.hpp`

.. code:: cppkokkos

  template <class ExecutionSpace, class XViewType, class YViewType, class AViewType>
  void ger(const ExecutionSpace& space, const char trans[],
           const typename AViewType::const_value_type& alpha,
           const XViewType& x, const YViewType& y, const AViewType& A);

  template <class XViewType, class YViewType, class AViewType>
  void ger(const char trans[], const typename AViewType::const_value_type& alpha,
           const XViewType& x, const YViewType& y, const AViewType& A);

Perform a rank-1 update of matrix ``A`` by vectors ``x`` and ``y`` with scaling factor ``alpha``

.. math::

   A = A + \alpha (x * y^T).

1. Iterate of the entries of ``A`` and update them with the ``alpha`` scaled product of the entries of ``x`` and ``y`` using resources from the ``space`` instance.
2. Iterate of the entries of ``A`` and update them with the ``alpha`` scaled product of the entries of ``x`` and ``y`` using resources from the default instance of ``typename AViewType::execution_space``.

The function will throw a runtime exception if ``A.extent(0) != X.extent(0) || A.extent(1) != Y.extent(0)``.

Parameters
==========

:space: execution space instance.

:trans: "T" for transpose and "H" for conjugate transpose of ``y``, all characters after the first one are ignored.

:alpha: scaling parameter for the rank update.

:x, y: vectors used to perform the rank update.

:A: matrix being updated.

Type Requirements
-----------------

- `execution_space` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_

- `XViewType` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 that satisfies

  - ``Kokkos::SpaceAccessibility<execution_space, typename XViewType::memory_space>::accessible``

- `YViewType` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 that satisfies

  - ``Kokkos::SpaceAccessibility<execution_space, typename YViewType::memory_space>::accessible``

- `AViewType` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 that satisfies

  - ``Kokkos::SpaceAccessibility<execution_space, typename AViewType::memory_space>::accessible``
  - ``std::is_same_v<typename AViewType::value_type, typename AViewType::non_const_value_type>``


Example
=======

.. literalinclude:: ../../../../example/wiki/blas/KokkosBlas2_wiki_ger.cpp
  :language: c++
