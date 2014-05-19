smoother: pre type = CHEBYSHEV
smoother: post type = ILUT
verbosity = test
coarse: max size = 2000   [default]
max levels = 10   [default]
debug: graph level = -1   [default]
number of equations = 1   [default]
transpose: use implicit = 0   [default]
smoother: pre or post = both   [default]
aggregation: type = uncoupled   [default]
multigrid algorithm = sa   [default]
repartition: enable = 0   [default]

Level 0
 Setup Smoother (MueLu::Ifpack2Smoother{type = CHEBYSHEV})
 Setup Smoother (MueLu::Ifpack2Smoother{type = ILUT})
 presmoother -> 
  chebyshev: ratio eigenvalue = 20
  chebyshev: min diagonal value = 2.22045e-16   [default]
  chebyshev: degree = 1   [default]
  chebyshev: eigenvalue max iterations = 10   [default]
  chebyshev: zero starting solution = 1   [default]
  chebyshev: assume matrix does not change = 0   [default]
 postsmoother -> 
  [empty list]
 
Level 1
 Prolongator smoothing (MueLu::SaPFactory)
  Build (MueLu::TentativePFactory)
   Build (MueLu::UncoupledAggregationFactory)
    Build (MueLu::CoalesceDropFactory)
    lightweight wrap = 1
    aggregation threshold = 0
    Dirichlet detection threshold = 0
    algorithm = original
    
   Ordering = 0   [default]
   MaxNeighAlreadySelected = 0   [default]
   MinNodesPerAggregate = 2   [default]
   MaxNodesPerAggregate = 2147483647   [default]
   UseOnePtAggregationAlgorithm = 0   [default]
   UseSmallAggregatesAggregationAlgorithm = 0   [default]
   UsePreserveDirichletAggregationAlgorithm = 0   [default]
   UseUncoupledAggregationAlgorithm = 1   [default]
   UseMaxLinkAggregationAlgorithm = 1   [default]
   UseIsolatedNodeAggregationAlgorithm = 1   [default]
   UseEmergencyAggregationAlgorithm = 1   [default]
   OnePt aggregate map name =    [default]
   SmallAgg aggregate map name =    [default]
   
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
 Setup Smoother (MueLu::Ifpack2Smoother{type = ILUT})
 presmoother -> 
  chebyshev: ratio eigenvalue = 20
  chebyshev: min diagonal value = 2.22045e-16   [default]
  chebyshev: degree = 1   [default]
  chebyshev: eigenvalue max iterations = 10   [default]
  chebyshev: zero starting solution = 1   [default]
  chebyshev: assume matrix does not change = 0   [default]
 postsmoother -> 
  [empty list]
 
Level 2
 Prolongator smoothing (MueLu::SaPFactory)
  Build (MueLu::TentativePFactory)
   Build (MueLu::UncoupledAggregationFactory)
    Build (MueLu::CoalesceDropFactory)
    lightweight wrap = 1
    aggregation threshold = 0
    Dirichlet detection threshold = 0
    algorithm = original
    
   Ordering = 0   [default]
   MaxNeighAlreadySelected = 0   [default]
   MinNodesPerAggregate = 2   [default]
   MaxNodesPerAggregate = 2147483647   [default]
   UseOnePtAggregationAlgorithm = 0   [default]
   UseSmallAggregatesAggregationAlgorithm = 0   [default]
   UsePreserveDirichletAggregationAlgorithm = 0   [default]
   UseUncoupledAggregationAlgorithm = 1   [default]
   UseMaxLinkAggregationAlgorithm = 1   [default]
   UseIsolatedNodeAggregationAlgorithm = 1   [default]
   UseEmergencyAggregationAlgorithm = 1   [default]
   OnePt aggregate map name =    [default]
   SmallAgg aggregate map name =    [default]
   
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
 Number of levels    = 3
 Operator complexity = 1.44
 Max Coarse Size     = 2000
 Implicit Transpose  = false
 
 matrix rows    nnz  nnz/row procs
 A 0    9999  29995     3.00  1
 A 1    3333   9997     3.00  1
 A 2    1111   3331     3.00  1
 
 Smoother (level 0) pre  : "Ifpack2::Chebyshev": {Initialized: true, Computed: true, "Ifpack2::Details::Chebyshev":{degree: 1, lambdaMax: 1.9506, alpha: 20, lambdaMin: 0.0975299}, Global matrix dimensions: [9999, 9999], Global nnz: 29995}
 Smoother (level 0) post : "Ifpack2::ILUT": {Initialized: true, Computed: true, Level-of-fill: 1, absolute threshold: 0, relative threshold: 1, relaxation value: 0, Global matrix dimensions: [9999, 9999], Global nnz: 29995}
 
 Smoother (level 1) pre  : "Ifpack2::Chebyshev": {Initialized: true, Computed: true, "Ifpack2::Details::Chebyshev":{degree: 1, lambdaMax: 1.94634, alpha: 20, lambdaMin: 0.0973171}, Global matrix dimensions: [3333, 3333], Global nnz: 9997}
 Smoother (level 1) post : "Ifpack2::ILUT": {Initialized: true, Computed: true, Level-of-fill: 1, absolute threshold: 0, relative threshold: 1, relaxation value: 0, Global matrix dimensions: [3333, 3333], Global nnz: 9997}
 
 Smoother (level 2) pre  : SuperLU solver interface
 Smoother (level 2) post : no smoother
 
