KokkosBlas::rotg
################

Defined in header: :code:`KokkosBlas1_rotg.hpp`

.. code:: c++

  template <class execution_space, class SViewType, class MViewType>
  void rotg(
    execution_space const& space,
    SViewType const& a, SViewType const& b,
    MViewType const& c, SViewType const& s);

Compute the Givens rotation coefficients ``c`` and ``s`` that align [a, b] along the direction of :math:`e_1`:

.. math::
   \begin{bmatrix} c && s\\ -s && c\end{bmatrix}\begin{bmatrix}a\\ b\end{bmatrix}=\begin{bmatrix}r\\ 0\end{bmatrix}

satisfying :math:`c^2+s^2=1` and :math:`r=\sqrt{a^2 + b^2}`

Parameters
==========

:space: execution space instance
:a, b: 0-D views. On input, the components of the 2D vector to rotate. On output, ``a`` contains the first component of the rotated vector and ``b`` contains the ``z`` parameter that provides an alternate way to define the rotation.
:c, s: cosine and sine of the rotation that rotates [a, b] onto :math:`e_1`

Type Requirements
-----------------

- `execution_space` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_

- `SViewType` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 0 that satisfies:

  - ``Kokkos::SpaceAccessibility<execution_space, typename SViewType::memory_space>::accessible == true``

- `MViewType` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 0 that satisfies:

  - ``Kokkos::SpaceAccessibility<execution_space, typename MViewType::memory_space>::accessible == true``
  - ``!KokkosKernels::ArithTraits<typename MViewType::value_type>::is_complex``

Example
=======

This example shows how to eliminate an entry using a Givens rotation.
It uses both :code:`rotg` to compute the rotation coefficients and :doc:`rot <blas1_rot>` to apply the rotation.

.. literalinclude:: ../../../../example/wiki/blas/KokkosBlas1_wiki_rotg_rot.cpp
  :language: c++
