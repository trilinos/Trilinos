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
 Setup Smoother (MueLu::Ifpack2Smoother{type = RELAXATION})
 relaxation: type = Symmetric Gauss-Seidel
 relaxation: sweeps = 1
 relaxation: damping factor = 1
 relaxation: zero starting solution = 1   [default]
 relaxation: backward mode = 0   [default]
 relaxation: use l1 = 0   [default]
 relaxation: l1 eta = 1.5   [default]
 relaxation: min diagonal value = 0   [default]
 relaxation: fix tiny diagonal entries = 0   [default]
 relaxation: check diagonal entries = 0   [default]
 relaxation: local smoothing indices = Teuchos::ArrayRCP<int>{ptr=0,lowerOffset=0,upperOffset=-1,size=0,node=0,strong_count=0,weak_count=0}   [default]
 
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
 
 Setup Smoother (MueLu::Ifpack2Smoother{type = RELAXATION})
 relaxation: type = Symmetric Gauss-Seidel
 relaxation: sweeps = 1
 relaxation: damping factor = 1
 relaxation: zero starting solution = 1   [default]
 relaxation: backward mode = 0   [default]
 relaxation: use l1 = 0   [default]
 relaxation: l1 eta = 1.5   [default]
 relaxation: min diagonal value = 0   [default]
 relaxation: fix tiny diagonal entries = 0   [default]
 relaxation: check diagonal entries = 0   [default]
 relaxation: local smoothing indices = Teuchos::ArrayRCP<int>{ptr=0,lowerOffset=0,upperOffset=-1,size=0,node=0,strong_count=0,weak_count=0}   [default]
 
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
  numRemapValues = 4
  alwaysKeepProc0 = 1
  repartition: print partition distribution = 0
  
 type = Interpolation
 implicit = 0
 implicit transpose = 0   [default]
 useSubcomm = 1   [default]
 write start = -1   [default]
 write end = -1   [default]
 
 Computing Ac (MueLu::RebalanceAcFactory)
 useSubcomm = 1   [default]
 
 Setup Smoother (MueLu::Amesos2Smoother{type = Superlu})
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
 
 Smoother (level 0) both : "Ifpack2::Relaxation": {Initialized: true, Computed: true, Type: Symmetric Gauss-Seidel, sweeps: 1, damping factor: 1, Global matrix dimensions: [9999, 9999], Global nnz: 29995}
 
 Smoother (level 1) both : "Ifpack2::Relaxation": {Initialized: true, Computed: true, Type: Symmetric Gauss-Seidel, sweeps: 1, damping factor: 1, Global matrix dimensions: [3335, 3335], Global nnz: 10015}
 
 Smoother (level 2) pre  : SuperLU solver interface
 Smoother (level 2) post : no smoother
 
