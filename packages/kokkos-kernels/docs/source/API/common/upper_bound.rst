KokkosKernels::upper_bound_thread and KokkosKernels::upper_bound_team
#####################################################################

Defined in header: :code:`KokkosKernels_UpperBound.hpp`

.. code:: c++

  template <typename ViewLike,
            typename Pred = LT<typename ViewLike::non_const_value_type>>
  KOKKOS_INLINE_FUNCTION typename ViewLike::size_type
  upper_bound_thread(const ViewLike& view,
                     const typename ViewLike::non_const_value_type& value,
                     Pred pred = Pred());

  template <typename TeamMember, typename ViewLike,
            typename Pred = LT<typename ViewLike::non_const_value_type>>
  KOKKOS_INLINE_FUNCTION typename ViewLike::size_type
  upper_bound_team(const TeamMember& handle, const ViewLike& view,
                   const typename ViewLike::non_const_value_type& value,
                   Pred pred = Pred());

These functions search a rank-1 ordered input and return the first index
``i`` such that ``pred(value, view(i))`` is true. With the default predicate,
that is the first position whose value is strictly greater than ``value``.

``upper_bound_thread`` performs the search in a single thread. ``upper_bound_team``
provides the same operation cooperatively across a Kokkos team.

Parameters
==========

:handle: Kokkos team handle used by the team overload.

:view: rank-1 ordered input to search.

:value: target value used in the comparison.

:pred: binary predicate that defines the ordering relation.

Type Requirements
-----------------

- ``ViewLike`` must be a rank-1 view-like type with ``size()`` and ``operator()``.

- ``ViewLike`` must expose ``non_const_value_type`` and ``size_type``.

- ``Pred`` must be callable as ``pred(value, view(i))`` and return a Boolean-like value.

- ``TeamMember`` must be a Kokkos team-policy member type for the team overload.

Notes
=====

- The input must already be ordered with respect to the predicate.

- If no entry satisfies the upper-bound condition, the return value is
  ``view.size()``.

Example
=======

.. code:: c++

  #include <Kokkos_Core.hpp>
  #include <KokkosKernels_UpperBound.hpp>

  using view_type = Kokkos::View<int*>;

  void example(const view_type& v) {
    auto idx = KokkosKernels::upper_bound_thread(v, 7);
    (void)idx;
  }
