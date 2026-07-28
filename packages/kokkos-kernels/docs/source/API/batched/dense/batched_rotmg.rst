KokkosBatched::Rotmg
####################

Defined in header: :code:`KokkosBatched_Rotmg.hpp`

.. code:: c++

   struct Rotmg {
     template <class DXViewType, class YViewType, class PViewType>
     KOKKOS_INLINE_FUNCTION static int invoke(const DXViewType &d1, const DXViewType &d2, const DXViewType &x1,
                                              const YViewType &y1, const PViewType &param);
   };

Constructs the modified Givens transformation

.. math::
  
  H = \begin{bmatrix} h_{11} & h_{12} \\ h_{21} & h_{22} \end{bmatrix}

that zeros out the second component of a 2D vector :math:`(\sqrt{d1} \times x1, \sqrt{d2} \times y1)^T`. The transformation can be applied to vectors :math:`x` and :math:`y` by calling the `Rotm` routine described in :doc:`batched_rotm`. This operation is equivalent to the BLAS routine `SROTMG <https://www.netlib.org/blas/srotmg.f>`_ or `DROTMG <https://www.netlib.org/blas/drotmg.f>`_ for single or double precision.
A :math:`flag` and matrix H are stored in a 5-element vector :math:`param` as follows:

.. math::

   \begin{align}
   param(0) &= flag \\
   param(1) &= h_{11} \\
   param(2) &= h_{21} \\
   param(3) &= h_{12} \\
   param(4) &= h_{22}
   \end{align}

The value of :math:`flag` can be -1, 0, 1, or -2, which corresponds to the following forms of H:

.. list-table::
   :widths: 25 25 25 25
   :header-rows: 1

   * - :math:`flag == -1`
     - :math:`flag == 0`
     - :math:`flag == 1`
     - :math:`flag == -2`
   * - :math:`\begin{bmatrix} h_{11} & h_{12} \\ h_{21} & h_{22} \end{bmatrix}`
     - :math:`\begin{bmatrix} 1.0 & h_{12} \\ h_{21} & 1.0 \end{bmatrix}`
     - :math:`\begin{bmatrix} h_{11} & 1.0 \\ -1.0 & h_{22} \end{bmatrix}`
     - :math:`\begin{bmatrix} 1.0 & 0.0 \\ 0.0 & 1.0 \end{bmatrix}`

Parameters
==========

:d1, d2: 0-D view. On input, the squares of the initial scaling factors for :math:`x` and :math:`y`. On output, they contains squares of the scaling factors to be applied to :math:`x` and :math:`y`, respectively.
:x1, y1: 0-D view. On input, the components of the 2D vector to rotate. On output, :math:`x1` contains the first component of the rotated vector, while :math:`y1` is unchanged on output.
:param: 1-D view. On output,  A length 5 vector, containing the parameters of the modified Givens rotation. The first element of :math:`param` is a flag that indicates the form of the rotation matrix :math:`H`, and the remaining four elements are the entries of the rotation matrix :math:`H`.

Type Requirements
-----------------

- `DXViewType` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 0 that satisfies:

  - ``std::is_same_v<typename DXViewType::non_const_value_type, typename DXViewType::value_type> == true``
  - ``!KokkosKernels::ArithTraits<typename DXViewType::value_type>::is_complex``

- `YViewType` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 0 that satisfies:

  - ``!KokkosKernels::ArithTraits<typename YViewType::value_type>::is_complex``

- `PViewType` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 that satisfies:

  - ``std::is_same_v<typename PViewType::non_const_value_type, typename PViewType::value_type> == true``
  - ``!KokkosKernels::ArithTraits<typename PViewType::value_type>::is_complex``

Example
=======

.. literalinclude:: ../../../../../example/batched_solve/serial_rotm.cpp
  :language: c++
  :linenos:
  :lines: 4-

output:

.. code::

   rotm/rotmg works correctly!
