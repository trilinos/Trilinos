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
