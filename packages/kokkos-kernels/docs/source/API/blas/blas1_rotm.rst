KokkosBlas::rotm
################

Defined in header: :code:`KokkosBlas1_rotm.hpp`

.. code:: c++

  template <class execution_space, class VectorView, class ParamView>
  void rotm(
    execution_space const& space,
    VectorView const& X, VectorView const& Y,
    ParamView const& param);
  
  template <class VectorView, class ParamView>
  void rotm(
    VectorView const& X, VectorView const& Y,
    ParamView const& param);

1. Applies the modified Givens rotations defined by the parameters in ``param`` to the pair of vectors ``X`` and ``Y``. This executes on the provided ``space`` instance.
2. Applies the modified Givens rotations defined by the parameters in ``param`` to the pair of vectors ``X`` and ``Y``. This executes on the default instance of type ``typename VectorView::execution_space``.

Please refer to the documentation of :doc:`rotmg <blas1_rotmg>` for the description of the data stored in ``param``.

The function will throw a runtime exception if ``X.extent(0) != Y.extent(0)``.

Parameters
==========

:space: execution space instance
:X, Y: pair of vectors to rotate
:params: 5-element view containing the control parameter and coefficients defining the rotation to apply

Type Requirements
-----------------

- `execution_space` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_

- `VectorView` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 that satisfies:

  - ``Kokkos::SpaceAccessibility<execution_space, typename VectorView::memory_space>::accessible == true``
  - ``std::is_same_v<typename VectorView::non_const_value_type, typename VectorView::value_type> == true``
  - ``!Kokkos::ArithTraits<typename VectorView::value_type>::is_complex``

- `ParamView` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 and compile-time extent 5 that satisfies:

  - ``Kokkos::SpaceAccessibility<execution_space, typename ParamView::memory_space>::accessible == true``
  - ``std::is_same_v<typename ParamView::non_const_value_type, typename ParamView::value_type> == true``
  - ``!Kokkos::ArithTraits<typename ParamView::value_type>::is_complex``

Example
=======

This example shows how to eliminate an entry using a modified Givens rotation.
It uses :doc:`rotmg <blas1_rotmg>` to compute the rotation parameters and ``rotm`` to apply the rotation.

.. literalinclude:: ../../../../example/wiki/blas/KokkosBlas1_wiki_rotmg_rotm.cpp
  :language: c++
