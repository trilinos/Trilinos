coarse: max size = 100
smoother: type = CHEBYSHEV
verbosity = test
max levels = 10   [default]
debug: graph level = -1   [default]
number of equations = 1   [default]
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
    
   Ordering = 0   [default]
   MaxNeighAlreadySelected = 0   [default]
   MinNodesPerAggregate = 2   [default]
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
   [empty list]
   
   Build (MueLu::CoarseMapFactory)
   [empty list]
   
  [empty list]
  
 Damping factor = 1.33333
 
 Transpose P (MueLu::TransPFactory)
 [empty list]
 
 Computing Ac (MueLu::RAPFactory)
 Keep AP Pattern = 0   [default]
 Keep RAP Pattern = 0   [default]
 CheckMainDiagonal = 0   [default]
 RepairMainDiagonal = 0   [default]
 
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
    
   Ordering = 0   [default]
   MaxNeighAlreadySelected = 0   [default]
   MinNodesPerAggregate = 2   [default]
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
   [empty list]
   
   Build (MueLu::CoarseMapFactory)
   [empty list]
   
  [empty list]
  
 Damping factor = 1.33333
 
 Transpose P (MueLu::TransPFactory)
 [empty list]
 
 Computing Ac (MueLu::RAPFactory)
 Keep AP Pattern = 0   [default]
 Keep RAP Pattern = 0   [default]
 CheckMainDiagonal = 0   [default]
 RepairMainDiagonal = 0   [default]
 
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
    
   Ordering = 0   [default]
   MaxNeighAlreadySelected = 0   [default]
   MinNodesPerAggregate = 2   [default]
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
   [empty list]
   
   Build (MueLu::CoarseMapFactory)
   [empty list]
   
  [empty list]
  
 Damping factor = 1.33333
 
 Transpose P (MueLu::TransPFactory)
 [empty list]
 
 Computing Ac (MueLu::RAPFactory)
 Keep AP Pattern = 0   [default]
 Keep RAP Pattern = 0   [default]
 CheckMainDiagonal = 0   [default]
 RepairMainDiagonal = 0   [default]
 
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
    
   Ordering = 0   [default]
   MaxNeighAlreadySelected = 0   [default]
   MinNodesPerAggregate = 2   [default]
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
   [empty list]
   
   Build (MueLu::CoarseMapFactory)
   [empty list]
   
  [empty list]
  
 Damping factor = 1.33333
 
 Transpose P (MueLu::TransPFactory)
 [empty list]
 
 Computing Ac (MueLu::RAPFactory)
 Keep AP Pattern = 0   [default]
 Keep RAP Pattern = 0   [default]
 CheckMainDiagonal = 0   [default]
 RepairMainDiagonal = 0   [default]
 
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
    
   Ordering = 0   [default]
   MaxNeighAlreadySelected = 0   [default]
   MinNodesPerAggregate = 2   [default]
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
   [empty list]
   
   Build (MueLu::CoarseMapFactory)
   [empty list]
   
  [empty list]
  
 Damping factor = 1.33333
 
 Transpose P (MueLu::TransPFactory)
 [empty list]
 
 Computing Ac (MueLu::RAPFactory)
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
 
 Smoother (level 0) both : Ifpack2::Chebyshev{status = initialized, computed, "Ifpack2::Details::Chebyshev":{degree: 1, lambdaMax: 1.9506, alpha: 2, lambdaMin: 0.975299}, global rows = 9999, global cols = 9999, global nnz  = 29995}
 
 Smoother (level 1) both : Ifpack2::Chebyshev{status = initialized, computed, "Ifpack2::Details::Chebyshev":{degree: 1, lambdaMax: 1.94634, alpha: 5, lambdaMin: 0.389268}, global rows = 3333, global cols = 3333, global nnz  = 9997}
 
 Smoother (level 2) both : Ifpack2::Chebyshev{status = initialized, computed, "Ifpack2::Details::Chebyshev":{degree: 1, lambdaMax: 1.95747, alpha: 4, lambdaMin: 0.489368}, global rows = 1111, global cols = 1111, global nnz  = 3331}
 
 Smoother (level 3) both : Ifpack2::Chebyshev{status = initialized, computed, "Ifpack2::Details::Chebyshev":{degree: 1, lambdaMax: 1.93464, alpha: 2.99461, lambdaMin: 0.646043}, global rows = 371, global cols = 371, global nnz  = 1111}
 
 Smoother (level 4) both : Ifpack2::Chebyshev{status = initialized, computed, "Ifpack2::Details::Chebyshev":{degree: 1, lambdaMax: 1.95476, alpha: 2.99194, lambdaMin: 0.653344}, global rows = 124, global cols = 124, global nnz  = 370}
 
 Smoother (level 5) pre  : SuperLU solver interface
 Smoother (level 5) post : no smoother
 
