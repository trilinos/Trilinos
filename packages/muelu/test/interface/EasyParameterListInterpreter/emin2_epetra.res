multigrid algorithm = emin
emin: iterative method = sd
emin: num iterations = 4
verbosity = test
coarse: max size = 2000   [default]
max levels = 10   [default]
debug: graph level = -1   [default]
number of equations = 1   [default]
transpose: use implicit = 0   [default]
smoother: pre or post = both   [default]
aggregation: type = uncoupled   [default]
emin: pattern = AkPtent   [default]
aggregation: visualize = 0   [default]
repartition: enable = 0   [default]

Level 0
 Setup Smoother (MueLu::IfpackSmoother{type = point relaxation stand-alone})
 relaxation: type = symmetric Gauss-Seidel   [unused]
 relaxation: sweeps = 1   [unused]
 relaxation: damping factor = 1   [unused]
 
Level 1
 Prolongator minimization (MueLu::EminPFactory)
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
   aggregation: preserve Dirichlet points = 0   [unused]
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
  
  Constraint (MueLu::ConstraintFactory)
   Ppattern (MueLu::PatternFactory)
   k = 1
   
  [empty list]
  
 emin: num iterations = 4
 emin: num reuse iterations = 1   [unused]
 P0 = Teuchos::RCP<Xpetra::Matrix<double, int, int, KokkosClassic::SerialNode, KokkosClassic::AltSparseOps<void, int, KokkosClassic::SerialNode, KokkosClassic::details::AltSparseOpsDefaultAllocator<int, KokkosClassic::SerialNode> > > >{ptr=0,node=0,strong_count=0,weak_count=0}   [unused]
 Keep P0 = 0
 Constraint0 = Teuchos::RCP<MueLu::Constraint<double, int, int, KokkosClassic::SerialNode, KokkosClassic::AltSparseOps<void, int, KokkosClassic::SerialNode, KokkosClassic::details::AltSparseOpsDefaultAllocator<int, KokkosClassic::SerialNode> > > >{ptr=0,node=0,strong_count=0,weak_count=0}   [unused]
 Keep Constraint0 = 0
 emin: iterative method = sd
 
 Transpose P (MueLu::TransPFactory)
 [empty list]
 
 Computing Ac (MueLu::RAPFactory)
 Keep AP Pattern = 0
 Keep RAP Pattern = 0
 implicit transpose = 0
 CheckMainDiagonal = 0
 RepairMainDiagonal = 0
 
 Setup Smoother (MueLu::IfpackSmoother{type = point relaxation stand-alone})
 relaxation: type = symmetric Gauss-Seidel   [unused]
 relaxation: sweeps = 1   [unused]
 relaxation: damping factor = 1   [unused]
 
Level 2
 Prolongator minimization (MueLu::EminPFactory)
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
   aggregation: preserve Dirichlet points = 0   [unused]
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
  
  Constraint (MueLu::ConstraintFactory)
   Ppattern (MueLu::PatternFactory)
   k = 1
   
  [empty list]
  
 emin: num iterations = 4
 emin: num reuse iterations = 1   [unused]
 P0 = Teuchos::RCP<Xpetra::Matrix<double, int, int, KokkosClassic::SerialNode, KokkosClassic::AltSparseOps<void, int, KokkosClassic::SerialNode, KokkosClassic::details::AltSparseOpsDefaultAllocator<int, KokkosClassic::SerialNode> > > >{ptr=0,node=0,strong_count=0,weak_count=0}   [unused]
 Keep P0 = 0
 Constraint0 = Teuchos::RCP<MueLu::Constraint<double, int, int, KokkosClassic::SerialNode, KokkosClassic::AltSparseOps<void, int, KokkosClassic::SerialNode, KokkosClassic::details::AltSparseOpsDefaultAllocator<int, KokkosClassic::SerialNode> > > >{ptr=0,node=0,strong_count=0,weak_count=0}   [unused]
 Keep Constraint0 = 0
 emin: iterative method = sd
 
 Transpose P (MueLu::TransPFactory)
 [empty list]
 
 Computing Ac (MueLu::RAPFactory)
 Keep AP Pattern = 0
 Keep RAP Pattern = 0
 implicit transpose = 0
 CheckMainDiagonal = 0
 RepairMainDiagonal = 0
 
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
 
 Smoother (level 0) both : MueLu::IfpackSmoother{type = point relaxation stand-alone}
 
 Smoother (level 1) both : MueLu::IfpackSmoother{type = point relaxation stand-alone}
 
 Smoother (level 2) pre  : MueLu::AmesosSmoother{type = Superlu}
 Smoother (level 2) post : no smoother
 
