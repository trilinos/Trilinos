KokkosBlas::set
###############

Defined in header: :code:`KokkosBlas1_set.hpp`

.. code:: c++

  struct SerialSet {
    template <typename ScalarType, typename AViewType>
    KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A);
  };

  template <typename MemberType>
  struct TeamSet {
    template <typename ScalarType, typename AViewType>
    KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A);
  };

  template <typename MemberType>
  struct TeamVectorSet {
    template <typename ScalarType, typename AViewType>
    KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A);
  };

Replaces the entries of `A` by the value `alpha`.

1. Execute operation using a single execution thread.
2. Execute operation using a Kokkos TeamThreadRange
3. Exceute operation using a Kokkos TeamVectorRange

Parameters
==========

:member: Kokkos team handle

:A: output view set to ``alpha``

:alpha: input value to set entries of ``A`` to.

Type Requirements
-----------------

- `member` must be a Kokkos `TeamPolicy::member_type <https://kokkos.org/kokkos-core-wiki/API/core/policies/TeamPolicy.html>`_

- `A` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1 or 2

- `alpha` must be convertible to the ``value_type`` of ``A``
