KokkosGraph Load-Balance Routines
#################################

Defined in header: :code:`KokkosGraph_LoadBalance.hpp`

.. code:: c++

  template <typename TaskSizesView>
  struct LoadBalanceResult {
    tasks_view_type tasks;
    ranks_view_type ranks;
  };

  template <typename ExecSpace, typename View>
  LoadBalanceResult<View>
  load_balance_exclusive(const ExecSpace& space, const View& scanTasks,
                         typename View::non_const_value_type totalWork);

  template <typename View>
  LoadBalanceResult<View>
  load_balance_exclusive(const View& scanTasks,
                         typename View::non_const_value_type totalWork);

  template <typename ExecSpace, typename View>
  LoadBalanceResult<View>
  load_balance(const ExecSpace& space, const View& taskSizes);

  template <typename View>
  LoadBalanceResult<View> load_balance(const View& taskSizes);

  template <typename TeamMember, typename View>
  LoadBalanceResult<View> load_balance_team(const TeamMember& handle,
                                            const View& taskSizes);

  template <typename TeamMember, typename OutView, typename InView>
  KOKKOS_INLINE_FUNCTION
  void inclusive_prefix_sum_team(const TeamMember& handle, const OutView& out,
                                 const InView& in);

These routines map irregular task sizes onto a dense work-item space.
The returned ``LoadBalanceResult`` stores, for each work item, the source task
index in ``tasks`` and the rank of that work item within the source task in
``ranks``.

Use ``load_balance`` when starting from raw task sizes. Use
``load_balance_exclusive`` when an exclusive prefix sum is already available.
Use ``load_balance_team`` for team-local cooperative load balancing.

Parameters
==========

:handle: Kokkos team handle used by the team-local overloads.

:space: execution space instance.

:taskSizes: rank-1 view whose entries store the amount of work in each task.

:scanTasks: rank-1 exclusive prefix sum of task sizes.

:totalWork: total number of work items represented by ``taskSizes`` or ``scanTasks``.

:out: rank-1 output view for the team-local inclusive prefix sum.

:in: rank-1 input view for the team-local inclusive prefix sum.

Type Requirements
-----------------

- ``View`` must be a rank-1 Kokkos view.

- The memory space of ``View`` must be accessible from the selected execution
  space.

- ``OutView`` and ``InView`` must be rank-1 Kokkos views with compatible extents.

- ``TeamMember`` must be a Kokkos team-policy member type for the team overloads.

Notes
=====

- For a task-size input ``{2, 1, 0, 3}``, the result encodes
  ``tasks = {0, 0, 1, 3, 3, 3}`` and ``ranks = {0, 1, 0, 0, 1, 2}``.

- ``load_balance_exclusive`` assumes the leading zero is part of the exclusive
  scan and omits it internally when building the result.

Example
=======

.. code:: c++

  #include <Kokkos_Core.hpp>
  #include <KokkosGraph_LoadBalance.hpp>

  using view_type = Kokkos::View<int*>;

  void example(const view_type& taskSizes) {
    auto result = KokkosGraph::load_balance(taskSizes);
    (void)result;
  }
