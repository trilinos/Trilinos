Level 0
 Setup Smoother (MueLu::Ifpack2Smoother{type = SCHWARZ})
 schwarz: overlap level = 1   [unused]
 schwarz: combine mode = Zero   [unused]
 schwarz: use reordering = 0   [unused]
 subdomain solver name = RILUK   [unused]
 subdomain solver parameters -> 
  fact: iluk level-of-fill = 0   [unused]
  fact: absolute threshold = 0   [unused]
  fact: relative threshold = 1   [unused]
  fact: relax value = 0   [unused]
 
Level 1
 Build (MueLu::TentativePFactory)
  Build (MueLu::UncoupledAggregationFactory)
   Build (MueLu::CoalesceDropFactory)
   aggregation: drop tol = 0   [default]
   aggregation: Dirichlet threshold = 0   [default]
   aggregation: drop scheme = classical   [default]
   lightweight wrap = 1
   
  aggregation: mode = new   [unused]
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
  
  Build (MueLu::AmalgamationFactory)
  [empty list]
  
  Nullspace factory (MueLu::NullspaceFactory)
  Fine level nullspace = Nullspace
  
  Build (MueLu::CoarseMapFactory)
  Striding info = {}   [default]
  Strided block id = -1   [default]
  Domain GID offsets = {0}   [default]
  
 [empty list]
 
 Transpose P (MueLu::TransPFactory)
 [empty list]
 
 Computing Ac (MueLu::RAPFactory)
 transpose: use implicit = 1
 Keep AP Pattern = 0   [default]
 Keep RAP Pattern = 0   [default]
 CheckMainDiagonal = 0   [default]
 RepairMainDiagonal = 0   [default]
 
 Setup Smoother (MueLu::Ifpack2Smoother{type = SCHWARZ})
 schwarz: overlap level = 1   [unused]
 schwarz: combine mode = Zero   [unused]
 schwarz: use reordering = 0   [unused]
 subdomain solver name = RILUK   [unused]
 subdomain solver parameters -> 
  fact: iluk level-of-fill = 0   [unused]
  fact: absolute threshold = 0   [unused]
  fact: relative threshold = 1   [unused]
  fact: relax value = 0   [unused]
 
Level 2
 Build (MueLu::TentativePFactory)
  Build (MueLu::UncoupledAggregationFactory)
   Build (MueLu::CoalesceDropFactory)
   aggregation: drop tol = 0   [default]
   aggregation: Dirichlet threshold = 0   [default]
   aggregation: drop scheme = classical   [default]
   lightweight wrap = 1
   
  aggregation: mode = new   [unused]
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
  
  Build (MueLu::AmalgamationFactory)
  [empty list]
  
  Nullspace factory (MueLu::NullspaceFactory)
  Fine level nullspace = Nullspace
  
  Build (MueLu::CoarseMapFactory)
  Striding info = {}   [default]
  Strided block id = -1   [default]
  Domain GID offsets = {0}   [default]
  
 [empty list]
 
 Transpose P (MueLu::TransPFactory)
 [empty list]
 
 Computing Ac (MueLu::RAPFactory)
 transpose: use implicit = 1
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

Smoother (level 0) both : "Ifpack2::AdditiveSchwarz": {Initialized: true, Computed: true, Overlap level: 0, Subdomain reordering: "none", Combine mode: "ZERO", Global matrix dimensions: [9999, 9999], Inner solver: {"Ifpack2::RILUK": {Initialized: true, Computed: true, Level-of-fill: 0, Global matrix dimensions: [9999, 9999], Global nnz: 29995}}}

Smoother (level 1) both : "Ifpack2::AdditiveSchwarz": {Initialized: true, Computed: true, Overlap level: 0, Subdomain reordering: "none", Combine mode: "ZERO", Global matrix dimensions: [3333, 3333], Inner solver: {"Ifpack2::RILUK": {Initialized: true, Computed: true, Level-of-fill: 0, Global matrix dimensions: [3333, 3333], Global nnz: 9997}}}

Smoother (level 2) pre  : SuperLU solver interface, direct solve
Smoother (level 2) post : no smoother

