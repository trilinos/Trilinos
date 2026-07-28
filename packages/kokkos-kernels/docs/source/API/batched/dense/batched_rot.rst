KokkosBatched::Rot
##################

Defined in header: :code:`KokkosBatched_Rot.hpp`

.. code:: c++

  template <bool Conj = false>
  struct SerialRot {
    template <typename XViewType, typename YViewType, typename CType, typename SType>
    KOKKOS_INLINE_FUNCTION static int invoke(const XViewType &x, const YViewType &y, const CType c, const SType s);
  };

  template <typename MemberType, bool Conj = false>
  struct TeamAxpy {
    template <typename XViewType, typename YViewType, typename CType, typename SType>
    KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const XViewType &x, const YViewType &y,
                                             const CType c, const SType s);
  };

  template <typename MemberType, bool Conj = false>
  struct TeamVectorRot {
    template <typename XViewType, typename YViewType, typename CType, typename SType>
    KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const XViewType &x, const YViewType &y,
                                             const CType c, const SType s);
  }; 

Applies a plane rotation to vectors :math:`x` and :math:`y`:

.. math::

   \begin{align}
   x &= c * x + s * y \\
   y &= -s * x + c * y \: \text{(if Conj is false)} \\
   y &= -s * \mathrm{conj}(x) + c * y \: \text{(if Conj is true)}
   \end{align}

1. For real vectors :math:`X` and :math:`Y`, this operation is equivalent to the BLAS routine ``SROT`` or ``DROT`` for single or double precision.
2. For complex vectors :math:`X` and :math:`Y`, this operation is equivalent to the BLAS routine ``CSROT`` or ``ZDROT`` for single or double precision if ``Conj`` is false. If ``Conj`` is true, this operation is equivalent to the BLAS routine ``CROT`` or ``ZROT`` for single or double precision.

Parameters
==========

:x: On input, :math:`x` is a length :math:`n` vector. On output, :math:`x` is overwritten by the rotated vector.
:y: On input, :math:`y` is a length :math:`n` vector. On output, :math:`y` is overwritten by the rotated vector.
:c: A scalar of cosine of the rotation (real scalar)
:s: A scalar of sine of the rotation (real or complex scalar)

Type Requirements
-----------------

- ``Conj`` must be a boolean template parameter that indicates whether the rotation is a conjugate rotation or not.
- ``MemberType`` must be a Kokkos team member handle (only for ``TeamRot`` and ``TeamVectorRot``).
- ``XViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 containing a vector :math:`X` that satisfies ``std::is_same_v<typename XViewType::value_type, typename XViewType::non_const_value_type>``
- ``YViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 containing a vector :math:`Y` that satisfies ``std::is_same_v<typename YViewType::value_type, typename YViewType::non_const_value_type>``
- ``CType`` must be a built-in real type like ``float``, or ``double``
- ``SType`` must be a built-in arithmetic type like ``float``, ``double``, ``Kokkos::complex<float>``, or ``Kokkos::complex<double>``

Example
=======

.. literalinclude:: ../../../../../example/batched_solve/serial_rot.cpp
  :language: c++
  :linenos:
  :lines: 4-

output:

.. code::

   rot/rotg works correctly!
