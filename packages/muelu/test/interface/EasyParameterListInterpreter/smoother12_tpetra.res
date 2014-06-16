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
 Setup Smoother (MueLu::Ifpack2Smoother{type = CHEBYSHEV})
 chebyshev: ratio eigenvalue = 2
 chebyshev: min diagonal value = 2.22045e-16   [default]
 chebyshev: degree = 1   [default]
 chebyshev: eigenvalue max iterations = 10   [default]
 chebyshev: zero starting solution = 1   [default]
 chebyshev: assume matrix does not change = 0   [default]
 
Level 1
 Prolongator smoothing (MueLu::SaPFactory)
  Build (MueLu::TentativePFactory)
   Build (MueLu::UncoupledAggregationFactory)
    Build (MueLu::CoalesceDropFactory)
    lightweight wrap = 1
    aggregation threshold = 0
    Dirichlet detection threshold = 0
    algorithm = original
    
   mode = old
   Ordering = 0
   MaxNeighAlreadySelected = 0
   MinNodesPerAggregate = 2
   MaxNodesPerAggregate = 2147483647
   UseOnePtAggregationAlgorithm = 0
   UseSmallAggregatesAggregationAlgorithm = 0
   UsePreserveDirichletAggregationAlgorithm = 0
   UseUncoupledAggregationAlgorithm = 1
   UseMaxLinkAggregationAlgorithm = 1
   UseIsolatedNodeAggregationAlgorithm = 1
   UseEmergencyAggregationAlgorithm = 1
   aggregation: enable phase 1 = 1   [unused]
   aggregation: enable phase 2a = 1   [unused]
   aggregation: enable phase 2b = 1   [unused]
   aggregation: enable phase 3 = 1   [unused]
   OnePt aggregate map name =
   SmallAgg aggregate map name =
   
   Build (MueLu::AmalgamationFactory)
   [empty list]
   
   Nullspace factory (MueLu::NullspaceFactory)
   Fine level nullspace = Nullspace
   
   Build (MueLu::CoarseMapFactory)
   Striding info = {}   [default]
   Strided block id = -1   [default]
   Domain GID offsets = {0}   [default]
   
  [empty list]
  
 Damping factor = 1.33333
 
 Transpose P (MueLu::TransPFactory)
 [empty list]
 
 Computing Ac (MueLu::RAPFactory)
 Keep AP Pattern = 0
 Keep RAP Pattern = 0
 implicit transpose = 0
 CheckMainDiagonal = 0
 RepairMainDiagonal = 0
 
 Setup Smoother (MueLu::Ifpack2Smoother{type = CHEBYSHEV})
 chebyshev: ratio eigenvalue = 5
 chebyshev: min diagonal value = 2.22045e-16   [default]
 chebyshev: degree = 1   [default]
 chebyshev: eigenvalue max iterations = 10   [default]
 chebyshev: zero starting solution = 1   [default]
 chebyshev: assume matrix does not change = 0   [default]
 
Level 2
 Prolongator smoothing (MueLu::SaPFactory)
  Build (MueLu::TentativePFactory)
   Build (MueLu::UncoupledAggregationFactory)
    Build (MueLu::CoalesceDropFactory)
    lightweight wrap = 1
    aggregation threshold = 0
    Dirichlet detection threshold = 0
    algorithm = original
    
   mode = old
   Ordering = 0
   MaxNeighAlreadySelected = 0
   MinNodesPerAggregate = 2
   MaxNodesPerAggregate = 2147483647
   UseOnePtAggregationAlgorithm = 0
   UseSmallAggregatesAggregationAlgorithm = 0
   UsePreserveDirichletAggregationAlgorithm = 0
   UseUncoupledAggregationAlgorithm = 1
   UseMaxLinkAggregationAlgorithm = 1
   UseIsolatedNodeAggregationAlgorithm = 1
   UseEmergencyAggregationAlgorithm = 1
   aggregation: enable phase 1 = 1   [unused]
   aggregation: enable phase 2a = 1   [unused]
   aggregation: enable phase 2b = 1   [unused]
   aggregation: enable phase 3 = 1   [unused]
   OnePt aggregate map name =
   SmallAgg aggregate map name =
   
   Build (MueLu::AmalgamationFactory)
   [empty list]
   
   Nullspace factory (MueLu::NullspaceFactory)
   Fine level nullspace = Nullspace
   
   Build (MueLu::CoarseMapFactory)
   Striding info = {}   [default]
   Strided block id = -1   [default]
   Domain GID offsets = {0}   [default]
   
  [empty list]
  
 Damping factor = 1.33333
 
 Transpose P (MueLu::TransPFactory)
 [empty list]
 
 Computing Ac (MueLu::RAPFactory)
 Keep AP Pattern = 0
 Keep RAP Pattern = 0
 implicit transpose = 0
 CheckMainDiagonal = 0
 RepairMainDiagonal = 0
 
 Setup Smoother (MueLu::Ifpack2Smoother{type = CHEBYSHEV})
 chebyshev: ratio eigenvalue = 4
 chebyshev: min diagonal value = 2.22045e-16   [default]
 chebyshev: degree = 1   [default]
 chebyshev: eigenvalue max iterations = 10   [default]
 chebyshev: zero starting solution = 1   [default]
 chebyshev: assume matrix does not change = 0   [default]
 
Level 3
 Prolongator smoothing (MueLu::SaPFactory)
  Build (MueLu::TentativePFactory)
   Build (MueLu::UncoupledAggregationFactory)
    Build (MueLu::CoalesceDropFactory)
    lightweight wrap = 1
    aggregation threshold = 0
    Dirichlet detection threshold = 0
    algorithm = original
    
   mode = old
   Ordering = 0
   MaxNeighAlreadySelected = 0
   MinNodesPerAggregate = 2
   MaxNodesPerAggregate = 2147483647
   UseOnePtAggregationAlgorithm = 0
   UseSmallAggregatesAggregationAlgorithm = 0
   UsePreserveDirichletAggregationAlgorithm = 0
   UseUncoupledAggregationAlgorithm = 1
   UseMaxLinkAggregationAlgorithm = 1
   UseIsolatedNodeAggregationAlgorithm = 1
   UseEmergencyAggregationAlgorithm = 1
   aggregation: enable phase 1 = 1   [unused]
   aggregation: enable phase 2a = 1   [unused]
   aggregation: enable phase 2b = 1   [unused]
   aggregation: enable phase 3 = 1   [unused]
   OnePt aggregate map name =
   SmallAgg aggregate map name =
   
   Build (MueLu::AmalgamationFactory)
   [empty list]
   
   Nullspace factory (MueLu::NullspaceFactory)
   Fine level nullspace = Nullspace
   
   Build (MueLu::CoarseMapFactory)
   Striding info = {}   [default]
   Strided block id = -1   [default]
   Domain GID offsets = {0}   [default]
   
  [empty list]
  
 Damping factor = 1.33333
 
 Transpose P (MueLu::TransPFactory)
 [empty list]
 
 Computing Ac (MueLu::RAPFactory)
 Keep AP Pattern = 0
 Keep RAP Pattern = 0
 implicit transpose = 0
 CheckMainDiagonal = 0
 RepairMainDiagonal = 0
 
 Setup Smoother (MueLu::Ifpack2Smoother{type = CHEBYSHEV})
 chebyshev: ratio eigenvalue = 2.99461
 chebyshev: min diagonal value = 2.22045e-16   [default]
 chebyshev: degree = 1   [default]
 chebyshev: eigenvalue max iterations = 10   [default]
 chebyshev: zero starting solution = 1   [default]
 chebyshev: assume matrix does not change = 0   [default]
 
Level 4
 Prolongator smoothing (MueLu::SaPFactory)
  Build (MueLu::TentativePFactory)
   Build (MueLu::UncoupledAggregationFactory)
    Build (MueLu::CoalesceDropFactory)
    lightweight wrap = 1
    aggregation threshold = 0
    Dirichlet detection threshold = 0
    algorithm = original
    
   mode = old
   Ordering = 0
   MaxNeighAlreadySelected = 0
   MinNodesPerAggregate = 2
   MaxNodesPerAggregate = 2147483647
   UseOnePtAggregationAlgorithm = 0
   UseSmallAggregatesAggregationAlgorithm = 0
   UsePreserveDirichletAggregationAlgorithm = 0
   UseUncoupledAggregationAlgorithm = 1
   UseMaxLinkAggregationAlgorithm = 1
   UseIsolatedNodeAggregationAlgorithm = 1
   UseEmergencyAggregationAlgorithm = 1
   aggregation: enable phase 1 = 1   [unused]
   aggregation: enable phase 2a = 1   [unused]
   aggregation: enable phase 2b = 1   [unused]
   aggregation: enable phase 3 = 1   [unused]
   OnePt aggregate map name =
   SmallAgg aggregate map name =
   
   Build (MueLu::AmalgamationFactory)
   [empty list]
   
   Nullspace factory (MueLu::NullspaceFactory)
   Fine level nullspace = Nullspace
   
   Build (MueLu::CoarseMapFactory)
   Striding info = {}   [default]
   Strided block id = -1   [default]
   Domain GID offsets = {0}   [default]
   
  [empty list]
  
 Damping factor = 1.33333
 
 Transpose P (MueLu::TransPFactory)
 [empty list]
 
 Computing Ac (MueLu::RAPFactory)
 Keep AP Pattern = 0
 Keep RAP Pattern = 0
 implicit transpose = 0
 CheckMainDiagonal = 0
 RepairMainDiagonal = 0
 
 Setup Smoother (MueLu::Ifpack2Smoother{type = CHEBYSHEV})
 chebyshev: ratio eigenvalue = 2.99194
 chebyshev: min diagonal value = 2.22045e-16   [default]
 chebyshev: degree = 1   [default]
 chebyshev: eigenvalue max iterations = 10   [default]
 chebyshev: zero starting solution = 1   [default]
 chebyshev: assume matrix does not change = 0   [default]
 
Level 5
 Prolongator smoothing (MueLu::SaPFactory)
  Build (MueLu::TentativePFactory)
   Build (MueLu::UncoupledAggregationFactory)
    Build (MueLu::CoalesceDropFactory)
    lightweight wrap = 1
    aggregation threshold = 0
    Dirichlet detection threshold = 0
    algorithm = original
    
   mode = old
   Ordering = 0
   MaxNeighAlreadySelected = 0
   MinNodesPerAggregate = 2
   MaxNodesPerAggregate = 2147483647
   UseOnePtAggregationAlgorithm = 0
   UseSmallAggregatesAggregationAlgorithm = 0
   UsePreserveDirichletAggregationAlgorithm = 0
   UseUncoupledAggregationAlgorithm = 1
   UseMaxLinkAggregationAlgorithm = 1
   UseIsolatedNodeAggregationAlgorithm = 1
   UseEmergencyAggregationAlgorithm = 1
   aggregation: enable phase 1 = 1   [unused]
   aggregation: enable phase 2a = 1   [unused]
   aggregation: enable phase 2b = 1   [unused]
   aggregation: enable phase 3 = 1   [unused]
   OnePt aggregate map name =
   SmallAgg aggregate map name =
   
   Build (MueLu::AmalgamationFactory)
   [empty list]
   
   Nullspace factory (MueLu::NullspaceFactory)
   Fine level nullspace = Nullspace
   
   Build (MueLu::CoarseMapFactory)
   Striding info = {}   [default]
   Strided block id = -1   [default]
   Domain GID offsets = {0}   [default]
   
  [empty list]
  
 Damping factor = 1.33333
 
 Transpose P (MueLu::TransPFactory)
 [empty list]
 
 Computing Ac (MueLu::RAPFactory)
 Keep AP Pattern = 0
 Keep RAP Pattern = 0
 implicit transpose = 0
 CheckMainDiagonal = 0
 RepairMainDiagonal = 0
 
 Setup Smoother (MueLu::Amesos2Smoother{type = Superlu})
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
 
 Smoother (level 0) both : "Ifpack2::Chebyshev": {Initialized: true, Computed: true, "Ifpack2::Details::Chebyshev":{degree: 1, lambdaMax: 1.9506, alpha: 2, lambdaMin: 0.975299}, Global matrix dimensions: [9999, 9999], Global nnz: 29995}
 
 Smoother (level 1) both : "Ifpack2::Chebyshev": {Initialized: true, Computed: true, "Ifpack2::Details::Chebyshev":{degree: 1, lambdaMax: 1.94634, alpha: 5, lambdaMin: 0.389268}, Global matrix dimensions: [3333, 3333], Global nnz: 9997}
 
 Smoother (level 2) both : "Ifpack2::Chebyshev": {Initialized: true, Computed: true, "Ifpack2::Details::Chebyshev":{degree: 1, lambdaMax: 1.95747, alpha: 4, lambdaMin: 0.489368}, Global matrix dimensions: [1111, 1111], Global nnz: 3331}
 
 Smoother (level 3) both : "Ifpack2::Chebyshev": {Initialized: true, Computed: true, "Ifpack2::Details::Chebyshev":{degree: 1, lambdaMax: 1.93464, alpha: 2.99461, lambdaMin: 0.646043}, Global matrix dimensions: [371, 371], Global nnz: 1111}
 
 Smoother (level 4) both : "Ifpack2::Chebyshev": {Initialized: true, Computed: true, "Ifpack2::Details::Chebyshev":{degree: 1, lambdaMax: 1.95476, alpha: 2.99194, lambdaMin: 0.653344}, Global matrix dimensions: [124, 124], Global nnz: 370}
 
 Smoother (level 5) pre  : SuperLU solver interface
 Smoother (level 5) post : no smoother
 
