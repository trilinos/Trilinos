repartition: enable = 1
transpose: use implicit = 1
verbosity = test
coarse: max size = 2000   [default]
max levels = 10   [default]
debug: graph level = -1   [default]
number of equations = 1   [default]
repartition: rebalance P and R = 1   [default]
smoother: pre or post = both   [default]
aggregation: type = uncoupled   [default]
multigrid algorithm = sa   [default]
repartition: partitioner = zoltan   [default]

Level 0
 Setup Smoother (MueLu::IfpackSmoother{type = point relaxation stand-alone})
 relaxation: type = symmetric Gauss-Seidel   [unused]
 relaxation: sweeps = 1   [unused]
 relaxation: damping factor = 1   [unused]
 
Level 1
 Build (MueLu::RebalanceTransferFactory)
  Build (MueLu::RepartitionFactory)
   Computing Ac (MueLu::RAPFactory)
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
    
   Build (MueLu::CoordinatesTransferFactory)
   write start = -1   [default]
   write end = -1   [default]
   
   Keep AP Pattern = 0
   Keep RAP Pattern = 0
   implicit transpose = 1
   CheckMainDiagonal = 0
   RepairMainDiagonal = 0
   
  startLevel = 2
  minRowsPerProcessor = 800
  nonzeroImbalance = 1.2
  remapPartitions = 1
  numRemapValues = 4   [unused]
  alwaysKeepProc0 = 1
  
 type = Interpolation
 implicit = 0
 implicit transpose = 0   [default]
 useSubcomm = 1   [default]
 write start = -1   [default]
 write end = -1   [default]
 
 Computing Ac (MueLu::RebalanceAcFactory)
 useSubcomm = 1   [default]
 
 Setup Smoother (MueLu::IfpackSmoother{type = point relaxation stand-alone})
 relaxation: type = symmetric Gauss-Seidel   [unused]
 relaxation: sweeps = 1   [unused]
 relaxation: damping factor = 1   [unused]
 
Level 2
 Build (MueLu::RebalanceTransferFactory)
  Build (MueLu::RepartitionFactory)
   Computing Ac (MueLu::RAPFactory)
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
    
   Build (MueLu::CoordinatesTransferFactory)
   write start = -1   [default]
   write end = -1   [default]
   
   Keep AP Pattern = 0
   Keep RAP Pattern = 0
   implicit transpose = 1
   CheckMainDiagonal = 0
   RepairMainDiagonal = 0
   
  startLevel = 2
  minRowsPerProcessor = 800
  nonzeroImbalance = 1.2
  remapPartitions = 1
  numRemapValues = 4
  alwaysKeepProc0 = 1
  
 type = Interpolation
 implicit = 0
 implicit transpose = 0   [default]
 useSubcomm = 1   [default]
 write start = -1   [default]
 write end = -1   [default]
 
 Computing Ac (MueLu::RebalanceAcFactory)
 useSubcomm = 1   [default]
 
 Setup Smoother (MueLu::AmesosSmoother{type = Superlu})
 presmoother -> 
  [empty list]
 
 
 --------------------------------------------------------------------------------
 ---                            Multigrid Summary                             ---
 --------------------------------------------------------------------------------
 Number of levels    = 3
 Operator complexity = 1.45
 Max Coarse Size     = 2000
 Implicit Transpose  = true
 
 matrix rows    nnz  nnz/row procs
 A 0    9999  29995     3.00  4
 A 1    3335  10015     3.00  4
 A 2    1112   3340     3.00  1
 
 Smoother (level 0) both : MueLu::IfpackSmoother{type = point relaxation stand-alone}
 
 Smoother (level 1) both : MueLu::IfpackSmoother{type = point relaxation stand-alone}
 
 Smoother (level 2) pre  : MueLu::AmesosSmoother{type = Superlu}
 Smoother (level 2) post : no smoother
 
