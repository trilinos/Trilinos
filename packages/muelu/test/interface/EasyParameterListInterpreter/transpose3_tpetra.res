transpose: use implicit = 1
repartition: enable = 1
repartition: rebalance P and R = 0
verbosity = test
coarse: max size = 2000   [default]
max levels = 10   [default]
debug: graph level = -1   [default]
number of equations = 1   [default]
smoother: pre or post = both   [default]
aggregation: type = uncoupled   [default]
multigrid algorithm = sa   [default]
aggregation: visualize = 0   [default]
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
       
      mode = old   [unused]
      Ordering = 0   [unused]
      MaxNeighAlreadySelected = 0   [unused]
      MinNodesPerAggregate = 2   [unused]
      MaxNodesPerAggregate = 2147483647   [unused]
      UseOnePtAggregationAlgorithm = 0   [unused]
      UsePreserveDirichletAggregationAlgorithm = 0   [unused]
      UseUncoupledAggregationAlgorithm = 1   [unused]
      UseMaxLinkAggregationAlgorithm = 1   [unused]
      UseIsolatedNodeAggregationAlgorithm = 1   [unused]
      UseEmergencyAggregationAlgorithm = 1   [unused]
      aggregation: preserve Dirichlet points = 0   [unused]
      aggregation: enable phase 1 = 1   [unused]
      aggregation: enable phase 2a = 1   [unused]
      aggregation: enable phase 2b = 1   [unused]
      aggregation: enable phase 3 = 1   [unused]
      OnePt aggregate map name = 
      
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
 implicit = 1
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
       
      mode = old   [unused]
      Ordering = 0   [unused]
      MaxNeighAlreadySelected = 0   [unused]
      MinNodesPerAggregate = 2   [unused]
      MaxNodesPerAggregate = 2147483647   [unused]
      UseOnePtAggregationAlgorithm = 0   [unused]
      UsePreserveDirichletAggregationAlgorithm = 0   [unused]
      UseUncoupledAggregationAlgorithm = 1   [unused]
      UseMaxLinkAggregationAlgorithm = 1   [unused]
      UseIsolatedNodeAggregationAlgorithm = 1   [unused]
      UseEmergencyAggregationAlgorithm = 1   [unused]
      aggregation: preserve Dirichlet points = 0   [unused]
      aggregation: enable phase 1 = 1   [unused]
      aggregation: enable phase 2a = 1   [unused]
      aggregation: enable phase 2b = 1   [unused]
      aggregation: enable phase 3 = 1   [unused]
      OnePt aggregate map name = 
      
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
 implicit = 1
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
 Operator complexity = 1.44
 Max Coarse Size     = 2000
 Implicit Transpose  = true
 
 matrix rows    nnz  nnz/row procs
 A 0    9999  29995     3.00  1
 A 1    3333   9997     3.00  1
 A 2    1111   3331     3.00  1
 
 Smoother (level 0) both : "Ifpack2::Relaxation": {Initialized: true, Computed: true, Type: Symmetric Gauss-Seidel, sweeps: 1, damping factor: 1, Global matrix dimensions: [9999, 9999], Global nnz: 29995}
 
 Smoother (level 1) both : "Ifpack2::Relaxation": {Initialized: true, Computed: true, Type: Symmetric Gauss-Seidel, sweeps: 1, damping factor: 1, Global matrix dimensions: [3333, 3333], Global nnz: 9997}
 
 Smoother (level 2) pre  : SuperLU solver interface
 Smoother (level 2) post : no smoother
 
