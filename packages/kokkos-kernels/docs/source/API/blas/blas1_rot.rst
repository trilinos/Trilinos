KokkosBlas::rot
###############

Defined in header: :code:`KokkosBlas1_rot.hpp`

.. code:: c++

  template <class execution_space, class VectorView, class MagnitudeView, class ScalarView>
  void rot(
    execution_space const& space,
    VectorView const& X, VectorView const& Y,
    MagnitudeView const& c, ScalarView const& s);
  
  template <class VectorView, class MagnitudeView, class ScalarView>
  void rot(
    VectorView const& X, VectorView const& Y,
    MagnitudeView const& c, ScalarView const& s);

Applies plane rotation ``(c, s)`` to vector pair ``X, Y``. This can be used to apply a Givens rotation to ``X`` and ``Y``,
if the coefficients ``c`` and ``s`` were previously computed by :doc:`KokkosBlas::rotg(a, b, c, s) <blas1_rotg>`.

1. Replaces the values of ``X`` and ``Y`` with ``c*X + s*Y`` and ``c*Y - s*X`` respectively, using the provided ``space`` instance.
2. Replaces the values of ``X`` and ``Y`` with ``c*X + s*Y`` and ``c*Y - s*X`` respectively, using the default instance of type ``typename VectorView::execution_space``.

The function will throw a runtime exception if ``X.extent(0) != Y.extent(0)``

Parameters
==========

:space: execution space instance
:X, Y: Pair of vectors to rotate
:c, s: cosine and sine of the angle rotation.

Type Requirements
-----------------

- `execution_space` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_
- `VectorView` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 that satisfies:

  - ``Kokkos::SpaceAccessibility<execution_space, typename VectorView::memory_space>::accessible == true``
  - ``std::is_same_v<typename VectorView::non_const_value_type, typename VectorView::value_type> == true``

- `MagnitudeView` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 0 that satisfies:

  - ``Kokkos::SpaceAccessibility<execution_space, typename MagnitudeView::memory_space>::accessible == true``
  - ``!KokkosKernels::ArithTraits<typename MagnitudeView::value_type>::is_complex``

- `ScalarView` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 0 that satisfies:

  - ``Kokkos::SpaceAccessibility<execution_space, typename ScalarView::memory_space>::accessible == true``

Example
=======

This example shows how to eliminate an entry using a Givens rotation.
It uses :doc:`rotg <blas1_rotg>` to compute the rotation coefficients and :code:`rot` to apply the rotation.

.. literalinclude:: ../../../../example/wiki/blas/KokkosBlas1_wiki_rotg_rot.cpp
  :language: c++

