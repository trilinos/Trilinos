Level 0
 Setup Smoother (MueLu::Ifpack2Smoother{type = CHEBYSHEV})
 chebyshev: ratio eigenvalue = 20
 chebyshev: min diagonal value = 2.22045e-16   [default]
 chebyshev: degree = 1   [default]
 chebyshev: eigenvalue max iterations = 10   [default]
 chebyshev: zero starting solution = 1   [default]
 chebyshev: assume matrix does not change = 0   [default]
 
 Setup Smoother (MueLu::Ifpack2Smoother{type = RILUK})
 [empty list]
 
Level 1
 Prolongator smoothing (MueLu::SaPFactory)
  Matrix filtering (MueLu::FilteredAFactory)
   Build (MueLu::CoalesceDropFactory)
    Build (MueLu::AmalgamationFactory)
    [empty list]

   aggregation: drop tol = 0   [default]
   aggregation: Dirichlet threshold = 0   [default]
   aggregation: drop scheme = classical   [default]
   lightweight wrap = 1
   
  filtered matrix: use lumping = 1   [default]
  filtered matrix: reuse graph = 1   [default]
  filtered matrix: reuse eigenvalue = 1   [default]
  
  Build (MueLu::TentativePFactory)
   Build (MueLu::UncoupledAggregationFactory)
   aggregation: mode = old   [default]
   aggregation: max agg size = -1   [default]
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
   
   Nullspace factory (MueLu::NullspaceFactory)
   Fine level nullspace = Nullspace
   
   Build (MueLu::CoarseMapFactory)
   Striding info = {}   [default]
   Strided block id = -1   [default]
   Domain GID offsets = {0}   [default]
   
  [empty list]
  
 sa: damping factor = 1.33   [default]
 sa: calculate eigenvalue estimate = 0   [default]
 sa: eigenvalue estimate num iterations = 10   [default]
 
 Transpose P (MueLu::TransPFactory)
 [empty list]
 
 Computing Ac (MueLu::RAPFactory)
 transpose: use implicit = 0   [default]
 Keep AP Pattern = 0   [default]
 Keep RAP Pattern = 0   [default]
 CheckMainDiagonal = 0   [default]
 RepairMainDiagonal = 0   [default]
 
 Setup Smoother (MueLu::Ifpack2Smoother{type = CHEBYSHEV})
 chebyshev: ratio eigenvalue = 20
 chebyshev: min diagonal value = 2.22045e-16   [default]
 chebyshev: degree = 1   [default]
 chebyshev: eigenvalue max iterations = 10   [default]
 chebyshev: zero starting solution = 1   [default]
 chebyshev: assume matrix does not change = 0   [default]
 
 Setup Smoother (MueLu::Ifpack2Smoother{type = RILUK})
 [empty list]
 
Level 2
 Prolongator smoothing (MueLu::SaPFactory)
  Matrix filtering (MueLu::FilteredAFactory)
   Build (MueLu::CoalesceDropFactory)
    Build (MueLu::AmalgamationFactory)
    [empty list]

   aggregation: drop tol = 0   [default]
   aggregation: Dirichlet threshold = 0   [default]
   aggregation: drop scheme = classical   [default]
   lightweight wrap = 1
   
  filtered matrix: use lumping = 1   [default]
  filtered matrix: reuse graph = 1   [default]
  filtered matrix: reuse eigenvalue = 1   [default]
  
  Build (MueLu::TentativePFactory)
   Build (MueLu::UncoupledAggregationFactory)
   aggregation: mode = old   [default]
   aggregation: max agg size = -1   [default]
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
   
   Nullspace factory (MueLu::NullspaceFactory)
   Fine level nullspace = Nullspace
   
   Build (MueLu::CoarseMapFactory)
   Striding info = {}   [default]
   Strided block id = -1   [default]
   Domain GID offsets = {0}   [default]
   
  [empty list]
  
 sa: damping factor = 1.33   [default]
 sa: calculate eigenvalue estimate = 0   [default]
 sa: eigenvalue estimate num iterations = 10   [default]
 
 Transpose P (MueLu::TransPFactory)
 [empty list]
 
 Computing Ac (MueLu::RAPFactory)
 transpose: use implicit = 0   [default]
 Keep AP Pattern = 0   [default]
 Keep RAP Pattern = 0   [default]
 CheckMainDiagonal = 0   [default]
 RepairMainDiagonal = 0   [default]
 
 Setup Smoother (MueLu::Amesos2Smoother{type = Superlu})
 presmoother -> 
  [empty list]
 

--------------------------------------------------------------------------------
---                            Multigrid Summary                             ---
--------------------------------------------------------------------------------
Number of levels    = 3
Operator complexity = 1.44

matrix rows    nnz  nnz/row procs
A 0    9999  29995     3.00  1
A 1    3333   9997     3.00  1
A 2    1111   3331     3.00  1

Smoother (level 0) pre  : "Ifpack2::Chebyshev": {Initialized: true, Computed: true, "Ifpack2::Details::Chebyshev":{degree: 1, lambdaMax: 1.95138, alpha: 20, lambdaMin: 0.0975692}, Global matrix dimensions: [9999, 9999], Global nnz: 29995}
Smoother (level 0) post : "Ifpack2::RILUK": {Initialized: true, Computed: true, Level-of-fill: 0, Global matrix dimensions: [9999, 9999], Global nnz: 29995}

Smoother (level 1) pre  : "Ifpack2::Chebyshev": {Initialized: true, Computed: true, "Ifpack2::Details::Chebyshev":{degree: 1, lambdaMax: 1.94982, alpha: 20, lambdaMin: 0.0974912}, Global matrix dimensions: [3333, 3333], Global nnz: 9997}
Smoother (level 1) post : "Ifpack2::RILUK": {Initialized: true, Computed: true, Level-of-fill: 0, Global matrix dimensions: [3333, 3333], Global nnz: 9997}

Smoother (level 2) pre  : SuperLU solver interface, direct solve
Smoother (level 2) post : no smoother

