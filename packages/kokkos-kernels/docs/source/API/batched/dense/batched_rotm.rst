KokkosBatched::Rotm
###################

Defined in header: :code:`KokkosBatched_Rotm.hpp`

.. code:: c++

  struct SerialRotm {
    template <typename XViewType, typename YViewType, typename ParamViewType>
    KOKKOS_INLINE_FUNCTION static int invoke(const XViewType &x, const YViewType &y, const ParamViewType &param);
  };

  template <typename MemberType>
  struct TeamRotm {
    template <typename XViewType, typename YViewType, typename ParamViewType>
    KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const XViewType &x, const YViewType &y,
                                             const ParamViewType &param);
  };

  template <typename MemberType>
  struct TeamVectorRotm {
    template <typename XViewType, typename YViewType, typename ParamViewType>
    KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const XViewType &x, const YViewType &y,
                                             const ParamViewType &param);
  };

Applies the modified Givens transformation

.. math::

  H = \begin{bmatrix} h_{11} & h_{12} \\ h_{21} & h_{22} \end{bmatrix}

to vectors :math:`x` and :math:`y`:

.. math::

   \begin{align}
   x(i) &= h_{11} \times x(i) + h_{12} \times y(i) \\
   y(i) &= h_{21} \times x(i) + h_{22} \times y(i)
   \end{align}

This operation is equivalent to the BLAS routine `SROTM <https://www.netlib.org/blas/srotm.f>`_ or `DROTM <https://www.netlib.org/blas/drotm.f>`_ for single or double precision. A :math:`flag` and matrix H are stored in a 5-element vector :math:`param` as follows:

.. math::

   \begin{align}
   param(0) &= flag \\
   param(1) &= h_{11} \\
   param(2) &= h_{21} \\
   param(3) &= h_{12} \\
   param(4) &= h_{22}
   \end{align}

where :math:`flag` is a scalar that indicates the form of the rotation matrix H. The value of :math:`flag` can be -1, 0, 1, or -2, which corresponds to the following forms of H:

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

.. note::

   This function does not support complex data types since the modified Givens rotation is not defined for complex numbers.

Parameters
==========

:x: On input, :math:`x` is a length :math:`n` vector. On output, :math:`x` is overwritten by the rotated vector.
:y: On input, :math:`y` is a length :math:`n` vector. On output, :math:`y` is overwritten by the rotated vector.
:param: A length 5 vector, a rank 1 view, containing the parameters of the modified Givens rotation. The first element of :math:`param` is a flag that indicates the form of the rotation matrix H, and the remaining four elements are the entries of the rotation matrix H.

Type Requirements
-----------------

- ``MemberType`` must be a Kokkos team member handle (only for ``TeamRotm`` and ``TeamVectorRotm``).
- ``XViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 containing a vector :math:`X` that satisfies:

  - ``std::is_same_v<typename XViewType::non_const_value_type, typename XViewType::value_type> == true``
  - ``!KokkosKernels::ArithTraits<typename XViewType::value_type>::is_complex``
- ``YViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 containing a vector :math:`Y` that satisfies:

  - ``std::is_same_v<typename YViewType::non_const_value_type, typename YViewType::value_type> == true``
  - ``!KokkosKernels::ArithTraits<typename YViewType::value_type>::is_complex``

- ``ParamViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 containing a vector :math:`param` that satisfies:

  - ``!KokkosKernels::ArithTraits<typename ParamViewType::value_type>::is_complex``

Example
=======

.. literalinclude:: ../../../../../example/batched_solve/serial_rotm.cpp
  :language: c++
  :linenos:
  :lines: 4-

output:

.. code::

   rotm/rotmg works correctly!
