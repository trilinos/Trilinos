KokkosGraph Merge-Path Routines
###############################

Defined in header: :code:`KokkosGraph_MergePath.hpp`

.. code:: c++

  struct StepperContext {
    size_t ai;
    size_t bi;
    size_t pi;
  };

  enum class StepDirection { a, b };

  template <typename AView, typename BView, typename Stepper, typename... Ctxs>
  KOKKOS_INLINE_FUNCTION
  void merge_path_thread(const AView& a, const BView& b, size_t pathLength,
                         const Stepper& stepper, Ctxs... ctxs);

  template <typename TeamHandle, typename AView, typename BViewLike,
            typename Stepper>
  KOKKOS_INLINE_FUNCTION
  void merge_path_team(const TeamHandle& handle, const AView& a,
                       const BViewLike& b, size_t pathLength,
                       const Stepper& stepper);

These routines traverse the merge path induced by two ordered rank-1 inputs.
Rather than producing a merged output directly, they invoke a user-provided
stepper for each step along the path.

``StepperContext`` records the position in the first input, the second input,
and the overall path. ``StepDirection`` indicates whether a given step consumes
an element from ``a`` or ``b``.

Parameters
==========

:handle: Kokkos team handle used by the team overload.

:a, b: rank-1 ordered inputs that define the merge path.

:pathLength: maximum number of path steps to process.

:stepper: callable object invoked for each path step.

:ctxs: optional contexts forwarded to the stepper by the thread overload.

Type Requirements
-----------------

- ``AView`` and ``BViewLike`` must be rank-1 view-like types.

- ``Stepper`` must be invocable with ``StepDirection`` followed by one, two, or
  three ``StepperContext`` arguments, depending on the overload and optional
  contexts in use.

- ``TeamHandle`` must be a Kokkos team-policy member type for the team overload.

Notes
=====

- ``merge_path_thread`` is useful when the caller wants custom per-step work
  rather than a materialized merged output.

- ``merge_path_team`` partitions the path across a Kokkos team and forwards a
  thread-level context to each stepper invocation.

Example
=======

.. code:: c++

  #include <Kokkos_Core.hpp>
  #include <KokkosGraph_MergePath.hpp>

  template <class AView, class BView>
  void example(const AView& a, const BView& b) {
    auto stepper = [](KokkosGraph::StepDirection dir,
                      KokkosGraph::StepperContext step) {
      (void)dir;
      (void)step;
    };
    KokkosGraph::merge_path_thread(a, b, a.size() + b.size(), stepper);
  }
