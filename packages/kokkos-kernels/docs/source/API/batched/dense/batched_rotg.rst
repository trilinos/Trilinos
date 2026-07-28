KokkosBatched::Rotg
###################

Defined in header: :code:`KokkosBatched_Rotg.hpp`

.. code:: c++

   struct Rotg {
     template <class SViewType, class MViewType>
     KOKKOS_INLINE_FUNCTION static int invoke(const SViewType &a, const SViewType &b, const MViewType &c,
                                              const SViewType &s);
   };

Constructs a plane rotation such that:

.. math::

   \begin{aligned}
   \begin{bmatrix} c & s\\ -s & c\end{bmatrix}\begin{bmatrix}a\\ b\end{bmatrix} &= \begin{bmatrix}r\\ 0\end{bmatrix} & \text{(for real arithmetic)} \\
   \begin{bmatrix} c & s\\ -\mathrm{conj}(s) & c\end{bmatrix}\begin{bmatrix}a\\ b\end{bmatrix} &= \begin{bmatrix}r\\ 0\end{bmatrix} & \text{(for complex arithmetic)}
   \end{aligned}

1. For real arithmetic, the rotation is defined such that :math:`c^2+s^2=1` and :math:`r=\pm\sqrt{a^2 + b^2}`. This operation is equivalent to the BLAS routine ``SROTG`` or ``DROTG`` for single or double precision.

2. For complex arithmetic, the rotation is defined such that :math:`c^2+|s|^2=1` and :math:`r=a/|a| \times \sqrt{a^2 + |b|^2}`. This operation is equivalent to the BLAS routine ``CROTG`` or ``ZROTG`` for single or double precision.

The rotation can be applied to vectors :math:`x` and :math:`y` by calling the `Rot` routine described in :doc:`batched_rot`. 

Parameters
==========

:a, b: 0-D views. On input, the components of the 2D vector to rotate. On output, ``a`` contains the first component of the rotated vector. For real arithmetic, ``b`` contains the ``z`` parameter that provides an alternate way to define the rotation. For complex arithmetic, ``b`` is unchanged on output.
:c, s: 0-D views. cosine and sine of the rotation that rotates [a, b] onto :math:`e_1`

Type Requirements
-----------------

- `SViewType` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 0 that satisfies:
   - ``std::is_same_v<typename SViewType::value_type, typename SViewType::non_const_value_type>``

- `MViewType` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 0 that satisfies:
   - ``std::is_same_v<typename MViewType::value_type, typename MViewType::non_const_value_type>``
   - ``!KokkosKernels::ArithTraits<typename MViewType::value_type>::is_complex``


Example
=======

.. literalinclude:: ../../../../../example/batched_solve/serial_rot.cpp
  :language: c++
  :linenos:
  :lines: 4-

output:

.. code::

   rot/rotg works correctly!
