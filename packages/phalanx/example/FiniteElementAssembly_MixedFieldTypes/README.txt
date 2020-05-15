This is an example for finite element assembly that mixes different
field types in Evaluators: MDField, Field, and Kokkos::View objects to
demonstrate that the various field implementations can be mixed in an
evaluation DAG.

For testing/code coverage purposes, must call the following functions
in the EvaluatorWithBaseImpl object for both Field and Kokkos::View
Types in the FieldManagerObject (other assembly example covers
MDFields):

addDependentField()
addContributedField()
addEvalautedField()

To cover this, IntegrateDiffusion and ProjectGradientToQP will use
Field objects while IntegrateSourceTerm and GatherSolution will use
Kokkos::Views.

This test will also cover manually enforcing user defined kokkos
layouts. The same objects above (for PHX::Field and Kokkos::View) will
pick a non-default layout. User specified layouts is not supported for
PHX::MDField. This will result in non-optimal performance of some
evaluators on certain devices, but we need to test this special case.
