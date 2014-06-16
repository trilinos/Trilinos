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
  repartition: print partition distribution = 0   [unused]
  
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
  repartition: print partition distribution = 0   [unused]
  
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
 Operator complexity = 1.44
 Max Coarse Size     = 2000
 Implicit Transpose  = true
 
 matrix rows    nnz  nnz/row procs
 A 0    9999  29995     3.00  1
 A 1    3333   9997     3.00  1
 A 2    1111   3331     3.00  1
 
 Smoother (level 0) both : MueLu::IfpackSmoother{type = point relaxation stand-alone}
 
 Smoother (level 1) both : MueLu::IfpackSmoother{type = point relaxation stand-alone}
 
 Smoother (level 2) pre  : MueLu::AmesosSmoother{type = Superlu}
 Smoother (level 2) post : no smoother
 
