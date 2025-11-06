KokkosBlas::rotmg
#################

Defined in header: :code:`KokkosBlas1_rotmg.hpp`

.. code:: c++

  template <class execution_space, class DXView, class YView, class PView>
  void rotmg(execution_space const& space, DXView const& d1, DXView const& d2,
             DXView const& x1, YView const& y1, PView const& param);
  
  template <class DXView, class YView, class PView>
  void rotmg(DXView const& d1, DXView const& d2, DXView const& x1,
             YView const& y1, PView const& param);

1. Computes the modified Givens rotation coefficients that zeroes out the second entry of vector :math:`[x1\sqrt{d1}, y1\sqrt{d2}]`, using the provided execution space instance.
2. Computes the modified Givens rotation coefficients that zeroes out the second entry of vector :math:`[x1\sqrt{d1}, y1\sqrt{d2}]`, using the default instance of :code:`PView::execution_space`.

The output flag parameter (stored in :code:`params(0)`) along with the other four elements of :code:`param` determine the matrix :math:`H` that applies the rotation:

flag = -1.0:

.. math::

   H=\begin{bmatrix} h_{11} && h_{12} \\ h_{21} && h_{22}\end{bmatrix}

flag = 0.0:

.. math::

   H=\begin{bmatrix} 1 && h_{12} \\ h_{21} && 1\end{bmatrix}

flag = 1.0:

.. math::

   H=\begin{bmatrix} h_{11} && 1 \\ -1 && h_{22}\end{bmatrix}

flag = -2.0:

.. math::

   H=\begin{bmatrix} 1 && 0\\ 0  && 1\end{bmatrix}

where elements :math:`1 \ldots 4` of :code:`params` contain :math:`h_{11}, h_{21}, h_{12}, h_{22}`.

Note that the values -1.0, 0.0 and 1.0 implied by the flag are not stored explicitly in :code:`params`.

Parameters
==========

:space: execution space instance
:d1, d2: Rank-0 views. On input, they contain the squares of the initial scaling factors for x,y respectively. On output, they contain the squares of the scaling factors to be applied to x,y after applying the rotation with :doc:`rotm <blas1_rotm>`.
:x1, y1: Rank-0 views containing the components of the vector to rotate.
:params: A 1D vector with compile-time dimension 5. On output, :code:`params(0)` is the control flag and :code:`params(1...4)` are the coefficients of the rotation matrix :math:`H` as described above.

Type Requirements
-----------------

- `execution_space` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_

- `DXView` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 0 that satisfies:

  - ``Kokkos::SpaceAccessibility<execution_space, typename DXView::memory_space>::accessible == true``
  - ``std::is_same_v<typename DXView::non_const_value_type, typename DXView::value_type> == true``
  - ``!KokkosKernels::ArithTraits<typename DXView::value_type>::is_complex``

- `YView` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 0 that satisfies:

  - ``Kokkos::SpaceAccessibility<execution_space, typename YView::memory_space>::accessible == true``
  - ``!KokkosKernels::ArithTraits<typename YView::value_type>::is_complex``

- `PView` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 and compile-time extent 5 that satisfies:

  - ``Kokkos::SpaceAccessibility<execution_space, typename PView::memory_space>::accessible == true``
  - ``std::is_same_v<typename PView::non_const_value_type, typename PView::value_type> == true``
  - ``!KokkosKernels::ArithTraits<typename PView::value_type>::is_complex``

Example
=======

This example shows how to eliminate an entry using a modified Givens rotation.
It uses :code:`rotmg` to compute the rotation parameters and :doc:`rotm <blas1_rotm>` to apply the rotation.

.. literalinclude:: ../../../../example/wiki/blas/KokkosBlas1_wiki_rotmg_rotm.cpp
  :language: c++
