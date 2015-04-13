Level 0
 Setup Smoother (MueLu::IfpackSmoother{type = Chebyshev})
 chebyshev: ratio eigenvalue = 20   [unused]
 
 Setup Smoother (MueLu::IfpackSmoother{type = ILUT})
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
 
 Setup Smoother (MueLu::IfpackSmoother{type = Chebyshev})
 chebyshev: ratio eigenvalue = 20   [unused]
 
 Setup Smoother (MueLu::IfpackSmoother{type = ILUT})
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
 
 Setup Smoother (MueLu::AmesosSmoother{type = Superlu})
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

Smoother (level 0) pre  : MueLu::IfpackSmoother{type = Chebyshev}
Smoother (level 0) post : MueLu::IfpackSmoother{type = ILUT}

Smoother (level 1) pre  : MueLu::IfpackSmoother{type = Chebyshev}
Smoother (level 1) post : MueLu::IfpackSmoother{type = ILUT}

Smoother (level 2) pre  : MueLu::AmesosSmoother{type = Superlu}
Smoother (level 2) post : no smoother

