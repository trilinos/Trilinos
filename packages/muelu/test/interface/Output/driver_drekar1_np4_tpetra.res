verbosity = test
number of equations = 1
coarse: max size = 1000
sa: use filtered matrix = 1
filtered matrix: use lumping = 1
smoother: type = CHEBYSHEV
aggregation: drop scheme = distance laplacian
repartition: enable = 1
repartition: min rows per proc = 2000
repartition: max imbalance = 1.327
repartition: start level = 1
repartition: remap parts = 1
repartition: keep proc 0 = 1
repartition: partitioner = zoltan2
max levels = 10   [default]
debug: graph level = -1   [default]
repartition: rebalance P and R = 1   [default]
transpose: use implicit = 0   [default]
smoother: pre or post = both   [default]
aggregation: type = uncoupled   [default]
multigrid algorithm = sa   [default]
problem: symmetric = 1   [default]
aggregation: export visualization data = 0   [default]
smoother: params ->
 chebyshev: degree = 2   [unused]
 chebyshev: ratio eigenvalue = 20   [unused]
 chebyshev: min eigenvalue = 1   [unused]
 chebyshev: zero starting solution = 1   [unused]
 chebyshev: eigenvalue max iterations = 10   [unused]
repartition: params ->
 algorithm = multijagged   [unused]

Level 0
 Setup Smoother (MueLu::Ifpack2Smoother{type = CHEBYSHEV})
 chebyshev: degree = 2
 chebyshev: ratio eigenvalue = 20
 chebyshev: min eigenvalue = 1
 chebyshev: zero starting solution = 1
 chebyshev: eigenvalue max iterations = 10
 chebyshev: min diagonal value = 2.22045e-16   [default]
 chebyshev: assume matrix does not change = 0   [default]

Level 1
 Build (MueLu::RebalanceTransferFactory)
  Build (MueLu::RepartitionFactory)
   Computing Ac (MueLu::RAPFactory)
    Prolongator smoothing (MueLu::SaPFactory)
     Matrix filtering (MueLu::FilteredAFactory)
      Build (MueLu::CoalesceDropFactory)
      lightweight wrap = 1
      aggregation: drop scheme = distance laplacian
      aggregation: drop tol = 0   [default]
      aggregation: Dirichlet threshold = 0   [default]

     filtered matrix: use lumping = 1   [unused]
     filtered matrix: reuse graph = 1   [default]
     filtered matrix: reuse eigenvalue = 1   [default]

     Build (MueLu::TentativePFactory)
      Build (MueLu::UncoupledAggregationFactory)
      aggregation: mode = old   [default]
      aggregation: max agg size = 2147483647   [default]
      aggregation: min agg size = 2   [default]
      aggregation: max selected neighbors = 0   [default]
      aggregation: ordering = natural   [default]
      aggregation: enable phase 1 = 1   [default]
      aggregation: enable phase 2a = 1   [default]
      aggregation: enable phase 2b = 1   [default]
      aggregation: enable phase 3 = 1   [default]
      aggregation: preserve Dirichlet points = 0   [default]
      UseOnePtAggregationAlgorithm = 0   [default]
      UsePreserveDirichletAggregationAlgorithm = 0   [default]
      UseUncoupledAggregationAlgorithm = 1   [default]
      UseMaxLinkAggregationAlgorithm = 1   [default]
      UseIsolatedNodeAggregationAlgorithm = 1   [default]
      UseEmergencyAggregationAlgorithm = 1   [default]
      OnePt aggregate map name =    [default]

      Build (MueLu::AmalgamationFactory)
      [empty list]

      Nullspace factory (MueLu::NullspaceFactory)
      Fine level nullspace = Nullspace

      Build (MueLu::CoarseMapFactory)
      Striding info = {}   [default]
      Strided block id = -1   [default]
      Domain GID offsets = {0}   [default]

     [empty list]

    sa: damping factor = 1.33333   [default]
    sa: calculate eigenvalue estimate = 0   [default]
    sa: eigenvalue estimate num iterations = 10   [default]

    Transpose P (MueLu::TransPFactory)
    [empty list]

   Build (MueLu::CoordinatesTransferFactory)
   write start = -1   [default]
   write end = -1   [default]

   transpose: use implicit = 0
   Keep AP Pattern = 0   [default]
   Keep RAP Pattern = 0   [default]
   CheckMainDiagonal = 0   [default]
   RepairMainDiagonal = 0   [default]

  repartition: start level = 1
  repartition: min rows per proc = 2000
  repartition: max imbalance = 1.327
  repartition: keep proc 0 = 1
  repartition: remap parts = 1
  repartition: print partition distribution = 0   [default]
  repartition: remap num values = 4   [default]

 type = Interpolation
 repartition: rebalance P and R = 1
 transpose: use implicit = 0   [default]
 useSubcomm = 1   [default]
 write start = -1   [default]
 write end = -1   [default]

 Build (MueLu::RebalanceTransferFactory)
 type = Restriction
 repartition: rebalance P and R = 1
 transpose: use implicit = 0
 useSubcomm = 1   [default]
 write start = -1   [default]
 write end = -1   [default]

 Computing Ac (MueLu::RebalanceAcFactory)
 useSubcomm = 1   [default]

 Setup Smoother (MueLu::Ifpack2Smoother{type = CHEBYSHEV})
 chebyshev: degree = 2
 chebyshev: ratio eigenvalue = 20
 chebyshev: min eigenvalue = 1
 chebyshev: zero starting solution = 1
 chebyshev: eigenvalue max iterations = 10
 chebyshev: min diagonal value = 2.22045e-16   [default]
 chebyshev: assume matrix does not change = 0   [default]

Level 2
 Build (MueLu::RebalanceTransferFactory)
  Build (MueLu::RepartitionFactory)
   Computing Ac (MueLu::RAPFactory)
    Prolongator smoothing (MueLu::SaPFactory)
     Matrix filtering (MueLu::FilteredAFactory)
      Build (MueLu::CoalesceDropFactory)
      lightweight wrap = 1
      aggregation: drop scheme = distance laplacian
      aggregation: drop tol = 0   [default]
      aggregation: Dirichlet threshold = 0   [default]

     filtered matrix: use lumping = 1   [unused]
     filtered matrix: reuse graph = 1   [default]
     filtered matrix: reuse eigenvalue = 1   [default]

     Build (MueLu::TentativePFactory)
      Build (MueLu::UncoupledAggregationFactory)
      aggregation: mode = old   [default]
      aggregation: max agg size = 2147483647   [default]
      aggregation: min agg size = 2   [default]
      aggregation: max selected neighbors = 0   [default]
      aggregation: ordering = natural   [default]
      aggregation: enable phase 1 = 1   [default]
      aggregation: enable phase 2a = 1   [default]
      aggregation: enable phase 2b = 1   [default]
      aggregation: enable phase 3 = 1   [default]
      aggregation: preserve Dirichlet points = 0   [default]
      UseOnePtAggregationAlgorithm = 0   [default]
      UsePreserveDirichletAggregationAlgorithm = 0   [default]
      UseUncoupledAggregationAlgorithm = 1   [default]
      UseMaxLinkAggregationAlgorithm = 1   [default]
      UseIsolatedNodeAggregationAlgorithm = 1   [default]
      UseEmergencyAggregationAlgorithm = 1   [default]
      OnePt aggregate map name =    [default]

      Build (MueLu::AmalgamationFactory)
      [empty list]

      Nullspace factory (MueLu::NullspaceFactory)
      Fine level nullspace = Nullspace

      Build (MueLu::CoarseMapFactory)
      Striding info = {}   [default]
      Strided block id = -1   [default]
      Domain GID offsets = {0}   [default]

     [empty list]

    sa: damping factor = 1.33333   [default]
    sa: calculate eigenvalue estimate = 0   [default]
    sa: eigenvalue estimate num iterations = 10   [default]

    Transpose P (MueLu::TransPFactory)
    [empty list]

   Build (MueLu::CoordinatesTransferFactory)
   write start = -1   [default]
   write end = -1   [default]

   transpose: use implicit = 0
   Keep AP Pattern = 0   [default]
   Keep RAP Pattern = 0   [default]
   CheckMainDiagonal = 0   [default]
   RepairMainDiagonal = 0   [default]

  repartition: start level = 1
  repartition: min rows per proc = 2000
  repartition: max imbalance = 1.327
  repartition: keep proc 0 = 1
  repartition: remap parts = 1
  repartition: print partition distribution = 0   [default]
  repartition: remap num values = 4   [default]

 type = Interpolation
 repartition: rebalance P and R = 1
 transpose: use implicit = 0   [default]
 useSubcomm = 1   [default]
 write start = -1   [default]
 write end = -1   [default]

 Build (MueLu::RebalanceTransferFactory)
 type = Restriction
 repartition: rebalance P and R = 1
 transpose: use implicit = 0
 useSubcomm = 1   [default]
 write start = -1   [default]
 write end = -1   [default]

 Computing Ac (MueLu::RebalanceAcFactory)
 useSubcomm = 1   [default]

 Setup Smoother (MueLu::Ifpack2Smoother{type = CHEBYSHEV})
 chebyshev: degree = 2
 chebyshev: ratio eigenvalue = 20
 chebyshev: min eigenvalue = 1
 chebyshev: zero starting solution = 1
 chebyshev: eigenvalue max iterations = 10
 chebyshev: min diagonal value = 2.22045e-16   [default]
 chebyshev: assume matrix does not change = 0   [default]

Level 3
 Build (MueLu::RebalanceTransferFactory)
  Build (MueLu::RepartitionFactory)
   Computing Ac (MueLu::RAPFactory)
    Prolongator smoothing (MueLu::SaPFactory)
     Matrix filtering (MueLu::FilteredAFactory)
      Build (MueLu::CoalesceDropFactory)
      lightweight wrap = 1
      aggregation: drop scheme = distance laplacian
      aggregation: drop tol = 0   [default]
      aggregation: Dirichlet threshold = 0   [default]

     filtered matrix: use lumping = 1   [unused]
     filtered matrix: reuse graph = 1   [default]
     filtered matrix: reuse eigenvalue = 1   [default]

     Build (MueLu::TentativePFactory)
      Build (MueLu::UncoupledAggregationFactory)
      aggregation: mode = old   [default]
      aggregation: max agg size = 2147483647   [default]
      aggregation: min agg size = 2   [default]
      aggregation: max selected neighbors = 0   [default]
      aggregation: ordering = natural   [default]
      aggregation: enable phase 1 = 1   [default]
      aggregation: enable phase 2a = 1   [default]
      aggregation: enable phase 2b = 1   [default]
      aggregation: enable phase 3 = 1   [default]
      aggregation: preserve Dirichlet points = 0   [default]
      UseOnePtAggregationAlgorithm = 0   [default]
      UsePreserveDirichletAggregationAlgorithm = 0   [default]
      UseUncoupledAggregationAlgorithm = 1   [default]
      UseMaxLinkAggregationAlgorithm = 1   [default]
      UseIsolatedNodeAggregationAlgorithm = 1   [default]
      UseEmergencyAggregationAlgorithm = 1   [default]
      OnePt aggregate map name =    [default]

      Build (MueLu::AmalgamationFactory)
      [empty list]

      Nullspace factory (MueLu::NullspaceFactory)
      Fine level nullspace = Nullspace

      Build (MueLu::CoarseMapFactory)
      Striding info = {}   [default]
      Strided block id = -1   [default]
      Domain GID offsets = {0}   [default]

     [empty list]

    sa: damping factor = 1.33333   [default]
    sa: calculate eigenvalue estimate = 0   [default]
    sa: eigenvalue estimate num iterations = 10   [default]

    Transpose P (MueLu::TransPFactory)
    [empty list]

   Build (MueLu::CoordinatesTransferFactory)
   write start = -1   [default]
   write end = -1   [default]

   transpose: use implicit = 0
   Keep AP Pattern = 0   [default]
   Keep RAP Pattern = 0   [default]
   CheckMainDiagonal = 0   [default]
   RepairMainDiagonal = 0   [default]

  repartition: start level = 1
  repartition: min rows per proc = 2000
  repartition: max imbalance = 1.327
  repartition: keep proc 0 = 1
  repartition: remap parts = 1
  repartition: print partition distribution = 0   [default]
  repartition: remap num values = 4   [default]

 type = Interpolation
 repartition: rebalance P and R = 1
 transpose: use implicit = 0   [default]
 useSubcomm = 1   [default]
 write start = -1   [default]
 write end = -1   [default]

 Build (MueLu::RebalanceTransferFactory)
 type = Restriction
 repartition: rebalance P and R = 1
 transpose: use implicit = 0
 useSubcomm = 1   [default]
 write start = -1   [default]
 write end = -1   [default]

 Computing Ac (MueLu::RebalanceAcFactory)
 useSubcomm = 1   [default]

 Setup Smoother (MueLu::Amesos2Smoother{type = Superlu})
 presmoother ->
  [empty list]


 --------------------------------------------------------------------------------
 ---                            Multigrid Summary                             ---
 --------------------------------------------------------------------------------
 Number of levels    = 4
 Operator complexity = 1.48
 Max Coarse Size     = 1000
 Implicit Transpose  = false

 matrix rows    nnz  nnz/row procs
 A 0    9999  29995     3.00  4
 A 1    3335  10015     3.00  1
 A 2    1111   3331     3.00  1
 A 3     371   1111     2.99  1

 Smoother (level 0) both : "Ifpack2::Chebyshev": {Initialized: true, Computed: true, "Ifpack2::Details::Chebyshev":{degree: 2, lambdaMax: 1.94747, alpha: 20, lambdaMin: 0.0973734}, Global matrix dimensions: [9999, 9999], Global nnz: 29995}

 Smoother (level 1) both : "Ifpack2::Chebyshev": {Initialized: true, Computed: true, "Ifpack2::Details::Chebyshev":{degree: 2, lambdaMax: 1.94405, alpha: 20, lambdaMin: 0.0972025}, Global matrix dimensions: [3335, 3335], Global nnz: 10015}

 Smoother (level 2) both : "Ifpack2::Chebyshev": {Initialized: true, Computed: true, "Ifpack2::Details::Chebyshev":{degree: 2, lambdaMax: 1.94378, alpha: 20, lambdaMin: 0.0971891}, Global matrix dimensions: [1111, 1111], Global nnz: 3331}

 Smoother (level 3) pre  : SuperLU solver interface, direct solve
 Smoother (level 3) post : no smoother
