smoother: pre type = CHEBYSHEV
smoother: post type = ILUT
verbosity = test
coarse: max size = 2000   [default]
max levels = 10   [default]
debug: graph level = -1   [default]
number of equations = 1   [default]
smoother: pre or post = both   [default]
aggregation: type = uncoupled   [default]
multigrid algorithm = sa   [default]
repartition: enable = 0   [default]

Level 0
 Setup Smoother (MueLu::IfpackSmoother{type = Chebyshev})
 Setup Smoother (MueLu::IfpackSmoother{type = ILU})
 presmoother -> 
  chebyshev: ratio eigenvalue = 30   [default]
  chebyshev: min eigenvalue = 0   [default]
  chebyshev: max eigenvalue = -1   [default]
  chebyshev: degree = 1   [default]
  chebyshev: min diagonal value = 0   [default]
  chebyshev: zero starting solution = 1   [default]
  chebyshev: operator inv diagonal = 0   [default]
  chebyshev: eigenvalue max iterations = 10   [default]
  chebyshev: use block mode = 0   [default]
  chebyshev: solve normal equations = 0   [default]
 postsmoother -> 
  fact: relax value = 0   [default]
  fact: absolute threshold = 0   [default]
  fact: relative threshold = 1   [default]
  fact: level-of-fill = 0   [default]
 
Level 1
 Prolongator smoothing (MueLu::SaPFactory)
  Build (MueLu::TentativePFactory)
   Build (MueLu::UncoupledAggregationFactory)
    Build (MueLu::CoalesceDropFactory)
    lightweight wrap = 1
    Dirichlet detection threshold = 0
    aggregation threshold = 0
    algorithm = original
    
   Ordering = 0   [default]
   MaxNeighAlreadySelected = 0   [default]
   MinNodesPerAggregate = 2   [default]
   UseOnePtAggregationAlgorithm = 1   [default]
   UseSmallAggregatesAggregationAlgorithm = 0   [default]
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
 
 Setup Smoother (MueLu::IfpackSmoother{type = Chebyshev})
 Setup Smoother (MueLu::IfpackSmoother{type = ILU})
 presmoother -> 
  chebyshev: ratio eigenvalue = 30   [default]
  chebyshev: min eigenvalue = 0   [default]
  chebyshev: max eigenvalue = -1   [default]
  chebyshev: degree = 1   [default]
  chebyshev: min diagonal value = 0   [default]
  chebyshev: zero starting solution = 1   [default]
  chebyshev: operator inv diagonal = 0   [default]
  chebyshev: eigenvalue max iterations = 10   [default]
  chebyshev: use block mode = 0   [default]
  chebyshev: solve normal equations = 0   [default]
 postsmoother -> 
  fact: relax value = 0   [default]
  fact: absolute threshold = 0   [default]
  fact: relative threshold = 1   [default]
  fact: level-of-fill = 0   [default]
 
Level 2
 Prolongator smoothing (MueLu::SaPFactory)
  Build (MueLu::TentativePFactory)
   Build (MueLu::UncoupledAggregationFactory)
    Build (MueLu::CoalesceDropFactory)
    lightweight wrap = 1
    Dirichlet detection threshold = 0
    aggregation threshold = 0
    algorithm = original
    disable Dirichlet detection = 0   [unused]
    
   Ordering = 0   [default]
   MaxNeighAlreadySelected = 0   [default]
   MinNodesPerAggregate = 2   [default]
   UseOnePtAggregationAlgorithm = 1   [default]
   UseSmallAggregatesAggregationAlgorithm = 0   [default]
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
 
 Setup Smoother (MueLu::AmesosSmoother{type = Superlu})
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
 
 Smoother (level 0) pre  : MueLu::IfpackSmoother{type = Chebyshev}
 Smoother (level 0) post : MueLu::IfpackSmoother{type = ILU}
 
 Smoother (level 1) pre  : MueLu::IfpackSmoother{type = Chebyshev}
 Smoother (level 1) post : MueLu::IfpackSmoother{type = ILU}
 
 Smoother (level 2) pre  : MueLu::AmesosSmoother{type = Superlu}
 Smoother (level 2) post : no smoother
 
