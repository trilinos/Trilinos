coarse: max size = 100
smoother: type = CHEBYSHEV
verbosity = test
max levels = 10   [default]
debug: graph level = -1   [default]
number of equations = 1   [default]
transpose: use implicit = 0   [default]
smoother: pre or post = both   [default]
aggregation: type = uncoupled   [default]
multigrid algorithm = sa   [default]
problem: symmetric = 1   [default]
aggregation: export visualization data = 0   [default]
repartition: enable = 0   [default]
smoother: params ->
 chebyshev: ratio eigenvalue = 2   [unused]
level 1 ->
 smoother: type = CHEBYSHEV
 smoother: params ->
  chebyshev: ratio eigenvalue = 5   [unused]
level 2 ->
 smoother: type = CHEBYSHEV
 smoother: params ->
  chebyshev: ratio eigenvalue = 4   [unused]

Level 0
 Setup Smoother (MueLu::IfpackSmoother{type = Chebyshev})
 chebyshev: ratio eigenvalue = 2

Level 1
 Prolongator smoothing (MueLu::SaPFactory)
  Build (MueLu::TentativePFactory)
   Build (MueLu::UncoupledAggregationFactory)
    Build (MueLu::CoalesceDropFactory)
    lightweight wrap = 1
    aggregation: drop tol = 0   [default]
    aggregation: Dirichlet threshold = 0   [default]
    aggregation: drop scheme = classical   [default]

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

 Computing Ac (MueLu::RAPFactory)
 transpose: use implicit = 0
 Keep AP Pattern = 0   [default]
 Keep RAP Pattern = 0   [default]
 CheckMainDiagonal = 0   [default]
 RepairMainDiagonal = 0   [default]

 Setup Smoother (MueLu::IfpackSmoother{type = Chebyshev})
 chebyshev: ratio eigenvalue = 5   [unused]

Level 2
 Prolongator smoothing (MueLu::SaPFactory)
  Build (MueLu::TentativePFactory)
   Build (MueLu::UncoupledAggregationFactory)
    Build (MueLu::CoalesceDropFactory)
    lightweight wrap = 1
    aggregation: drop tol = 0   [default]
    aggregation: Dirichlet threshold = 0   [default]
    aggregation: drop scheme = classical   [default]

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

 Computing Ac (MueLu::RAPFactory)
 transpose: use implicit = 0
 Keep AP Pattern = 0   [default]
 Keep RAP Pattern = 0   [default]
 CheckMainDiagonal = 0   [default]
 RepairMainDiagonal = 0   [default]

 Setup Smoother (MueLu::IfpackSmoother{type = Chebyshev})
 chebyshev: ratio eigenvalue = 4   [unused]

Level 3
 Prolongator smoothing (MueLu::SaPFactory)
  Build (MueLu::TentativePFactory)
   Build (MueLu::UncoupledAggregationFactory)
    Build (MueLu::CoalesceDropFactory)
    lightweight wrap = 1
    aggregation: drop tol = 0   [default]
    aggregation: Dirichlet threshold = 0   [default]
    aggregation: drop scheme = classical   [default]

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

 Computing Ac (MueLu::RAPFactory)
 transpose: use implicit = 0
 Keep AP Pattern = 0   [default]
 Keep RAP Pattern = 0   [default]
 CheckMainDiagonal = 0   [default]
 RepairMainDiagonal = 0   [default]

 Setup Smoother (MueLu::IfpackSmoother{type = Chebyshev})
 chebyshev: ratio eigenvalue = 2.99461   [unused]

Level 4
 Prolongator smoothing (MueLu::SaPFactory)
  Build (MueLu::TentativePFactory)
   Build (MueLu::UncoupledAggregationFactory)
    Build (MueLu::CoalesceDropFactory)
    lightweight wrap = 1
    aggregation: drop tol = 0   [default]
    aggregation: Dirichlet threshold = 0   [default]
    aggregation: drop scheme = classical   [default]

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

 Computing Ac (MueLu::RAPFactory)
 transpose: use implicit = 0
 Keep AP Pattern = 0   [default]
 Keep RAP Pattern = 0   [default]
 CheckMainDiagonal = 0   [default]
 RepairMainDiagonal = 0   [default]

 Setup Smoother (MueLu::IfpackSmoother{type = Chebyshev})
 chebyshev: ratio eigenvalue = 2.99194   [unused]

Level 5
 Prolongator smoothing (MueLu::SaPFactory)
  Build (MueLu::TentativePFactory)
   Build (MueLu::UncoupledAggregationFactory)
    Build (MueLu::CoalesceDropFactory)
    lightweight wrap = 1
    aggregation: drop tol = 0   [default]
    aggregation: Dirichlet threshold = 0   [default]
    aggregation: drop scheme = classical   [default]

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

 Computing Ac (MueLu::RAPFactory)
 transpose: use implicit = 0
 Keep AP Pattern = 0   [default]
 Keep RAP Pattern = 0   [default]
 CheckMainDiagonal = 0   [default]
 RepairMainDiagonal = 0   [default]

 Setup Smoother (MueLu::AmesosSmoother{type = Superlu})
 presmoother ->
  [empty list]


 --------------------------------------------------------------------------------
 ---                            Multigrid Summary                             ---
 --------------------------------------------------------------------------------
 Number of levels    = 6
 Operator complexity = 1.50
 Max Coarse Size     = 100
 Implicit Transpose  = false

 matrix rows    nnz  nnz/row procs
 A 0    9999  29995     3.00  1
 A 1    3333   9997     3.00  1
 A 2    1111   3331     3.00  1
 A 3     371   1111     2.99  1
 A 4     124    370     2.98  1
 A 5      42    124     2.95  1

 Smoother (level 0) both : MueLu::IfpackSmoother{type = Chebyshev}

 Smoother (level 1) both : MueLu::IfpackSmoother{type = Chebyshev}

 Smoother (level 2) both : MueLu::IfpackSmoother{type = Chebyshev}

 Smoother (level 3) both : MueLu::IfpackSmoother{type = Chebyshev}

 Smoother (level 4) both : MueLu::IfpackSmoother{type = Chebyshev}

 Smoother (level 5) pre  : MueLu::AmesosSmoother{type = Superlu}
 Smoother (level 5) post : no smoother
