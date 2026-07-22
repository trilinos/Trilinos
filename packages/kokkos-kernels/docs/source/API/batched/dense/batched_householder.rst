KokkosBatched::Householder
##########################

Defined in header: :code:`KokkosBatched_Householder_Decl.hpp`

.. code:: c++

    template <typename ArgSide>
    struct SerialHouseholder {
      template <typename aViewType, typename tauViewType>
      KOKKOS_INLINE_FUNCTION static int invoke(const aViewType &a, const tauViewType &tau);
    };

    template <typename ArgSide>
    struct TeamVectorHouseholder {
      template <typename MemberType, typename aViewType, typename tauViewType>
      KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const aViewType &a, const tauViewType &tau);
    };

Computes the Householder reflector associated with vector :math:`a` and stores the result in-place using the LAPACK format.
.. note::

   that the scaling factor :math:`tau` is computed as the inverse of the corresponding LAPACK scaling factor.

1. Compute the Householder reflector in serial at the inner most parallel level.
2. Compute the Householder reflector using parallelism exposed with the `member` handle.

Parameters
==========

:member: Kokkos team member handle
:a: On input a vector, on output stores the signed norm of :math:`a` in the first entry and a scaled Householder reflector in the remaining entries.
:tau: A vector containing the scaling factor for the reflector on output.

Type Requirements
-----------------

- ``MemberType`` must be a Kokkos team member handle.
- ``AViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2 that satisfies ``std::is_same_v<typename AViewType::value_type, typename AViewType::non_const_value_type>``
- ``tauViewType`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 that satisfies ``std::is_same_v<typename tauViewType::value_type, typename tauViewType::non_const_value_type>``

Example
=======

..
   .. literalinclude:: ../../../../../example/householder.cpp
     :language: c++

   output:

   .. code::

      ger works correctly!
