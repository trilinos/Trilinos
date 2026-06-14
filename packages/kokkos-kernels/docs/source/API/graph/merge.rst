KokkosGraph Merge Routines
##########################

Defined in header: :code:`KokkosGraph_Merge.hpp`

.. code:: c++

  template <typename CView, typename AView, typename BView>
  KOKKOS_INLINE_FUNCTION
  void merge_into_thread(const CView& c, const AView& a, const BView& b);

  template <typename TeamMember, typename CView, typename AView, typename BView>
  KOKKOS_INLINE_FUNCTION
  void merge_into_team(const TeamMember& handle, const CView& c,
                       const AView& a, const BView& b);

  template <typename ExecSpace, typename CView, typename AView, typename BView>
  void merge_into(const ExecSpace& space, const CView& c,
                  const AView& a, const BView& b);

  template <typename CView, typename AView, typename BView>
  void merge_into(const CView& c, const AView& a, const BView& b);

  template <typename ExecSpace, typename CView, typename AView, typename BView>
  void resize_and_merge_into(const ExecSpace& space, CView& c,
                             const AView& a, const BView& b);

  template <typename CView, typename AView, typename BView>
  void resize_and_merge_into(CView& c, const AView& a, const BView& b);

  template <typename CView, typename ExecSpace, typename AView, typename BView>
  CView merge(const ExecSpace& space, const AView& a, const BView& b);

  template <typename CView, typename AView, typename BView>
  CView merge(const AView& a, const BView& b);

These routines merge two ordered rank-1 inputs into a single ordered result.
They preserve stability with respect to the input sequences and expose thread,
team, and global interfaces.

Use ``merge_into_*`` when the output storage is already available. Use
``resize_and_merge_into`` when the output should be resized automatically. Use
``merge`` when the routine should allocate and return the merged view.

Parameters
==========

:handle: Kokkos team handle used by the team overload.

:space: execution space instance.

:a, b: rank-1 ordered input views to merge.

:c: rank-1 output view that stores the merged result.

Type Requirements
-----------------

- ``AView``, ``BView``, and ``CView`` must be rank-1 Kokkos views.

- ``c.size()`` must equal ``a.size() + b.size()`` for the ``merge_into`` overloads.

- The memory spaces of the inputs and output must be accessible from the
  selected execution space.

- The input value types must support the ordering relation used by the merge.

Notes
=====

- The inputs are assumed to be sorted in nondecreasing order.

- If one input is empty, the result is a copy of the other input.

Example
=======

.. code:: c++

  #include <Kokkos_Core.hpp>
  #include <KokkosGraph_Merge.hpp>

  using view_type = Kokkos::View<int*>;

  void example(const view_type& a, const view_type& b) {
    auto c = KokkosGraph::merge<view_type>(a, b);
    (void)c;
  }
