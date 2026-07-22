KokkosBatched::ApplyHouseholder
###############################

Defined in header: :code:`KokkosBatched_ApplyHouseholder_Decl.hpp`

.. code:: c++

    template <typename ArgSide, typename ArgTrans = Trans::NoTranspose>
    struct SerialApplyHouseholder {
      template <typename uViewType, typename tauViewType, typename AViewType, typename wViewType>
      KOKKOS_INLINE_FUNCTION static int invoke(const uViewType &u2, const tauViewType &tau, const AViewType &A,
      const wViewType &w);
    };

    template <typename MemberType, typename ArgSide>
    struct TeamVectorApplyHouseholder {
      template <typename uViewType, typename tauViewType, typename AViewType, typename wViewType>
      KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const uViewType &u2, const tauViewType &tau,
      const AViewType &A, const wViewType &w);
    };

Apply the Householder reflection associated with vector :math:`u2` and scaling coefficient :math:`tau` to matrix :math:`A`.

.. math::
   A = (I - \frac{1}{\tau}*u2*u2^H)*A

.. note::

   that the scaling factor :math:`tau` is computed as the inverse of the corresponding LAPACK scaling factor.

1. Apply the Householder reflector :math:`u2` to :math:`A` in serial at the inner most parallel level.
2. Apply the Householder reflector :math:`u2` to :math:`A` using parallelism exposed with the `member` handle.

Parameters
==========

:member: Kokkos team member handle.
:u2: The reflection vector.
:tau: The scaling factor associated with the reflection vector.
:A: The matrix reflected about vector :math:`u2`.
:w: Memory buffer required to perform the reflection on :math:`A`.

Type Requirements
-----------------

- ``MemberType`` must be a Kokkos team member handle.
- ``ArgSide`` must be a ``KokkosBatched::Side`` that indicates if the reflection is applied on the left or right of :math:`A`.
- ``ArgTrans`` must be a ``KokkosBatched::Trans`` that indicates if the reflection or it (conjugate) transpose is applied to :math:`A`.
- ``uViewType``  must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1.
- ``tauViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 0 or 1.
- ``AViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 or 2 that satisfies ``std::is_same_v<typename AViewType::value_type, typename AViewType::non_const_value_type>``
- ``wViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1.

Example
=======

..
   .. literalinclude:: ../../../../../example/apply_householder.cpp
     :language: c++

   output:

   .. code::

      ger works correctly!
