verbosity = high
coarse: max size = 100
transpose: use implicit = 0
max levels = 4
number of equations = 1   [default]
level 1 -> 
 P = Teuchos::RCP<Xpetra::Matrix<ignored> >{ptr=,node=,strong_count=6,weak_count=0}
 Nullspace = Teuchos::RCP<Xpetra::MultiVector<double, int, int, Kokkos::Compat::KokkosDeviceWrapperNode<Kokkos::Serial> > >{ptr=,node=,strong_count=6,weak_count=0}

Clearing old data (if any)
MueLu::Amesos2Smoother: using "Superlu"
MueLu::AmesosSmoother: using "Superlu"
Using default factory (MueLu::SmootherFactory{pre = MueLu::DirectSolver{type = }, post = null}) for building 'CoarseSolver'.
Using default factory (MueLu::SmootherFactory{pre = MueLu::TrilinosSmoother{type = RELAXATION}, post = pre}) for building 'Smoother'.
Level 0
 Setup Smoother (MueLu::Ifpack2Smoother{type = RELAXATION})
  "Ifpack2::Relaxation": {Initialized: true, Computed: true, Type: Symmetric Gauss-Seidel, sweeps: 1, damping factor: 1, Global matrix dimensions: [10000, 10000], Global nnz: 49600}
MueLu::Amesos2Smoother: using "Superlu"
MueLu::AmesosSmoother: using "Superlu"
Using default factory (MueLu::SmootherFactory{pre = MueLu::DirectSolver{type = }, post = null}) for building 'CoarseSolver'.
Using default factory (MueLu::SmootherFactory{pre = MueLu::TrilinosSmoother{type = RELAXATION}, post = pre}) for building 'Smoother'.
Using default factory (MueLu::AmalgamationFactory) for building 'UnAmalgamationInfo'.
Level 1
 Transpose P (MueLu::TransPFactory)
  R size =  1700 x 10000, nnz = 24807
  R Load balancing info
  R   # active processes: 1/1
  R   # rows per proc   : avg = 1.70e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
  R   #  nnz per proc   : avg = 2.48e+04,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
 Computing Ac (MueLu::RAPFactory)
  MxM: A x P
   Matrix product nnz per row estimate = 5
  MxM: R x (AP) (explicit)
   Matrix product nnz per row estimate = 18
  Ac size =  1700 x 1700, nnz = 14928
  Ac Load balancing info
  Ac   # active processes: 1/1
  Ac   # rows per proc   : avg = 1.70e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
  Ac   #  nnz per proc   : avg = 1.49e+04,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
 Setup Smoother (MueLu::Ifpack2Smoother{type = RELAXATION})
  "Ifpack2::Relaxation": {Initialized: true, Computed: true, Type: Symmetric Gauss-Seidel, sweeps: 1, damping factor: 1, Global matrix dimensions: [1700, 1700], Global nnz: 14928}
MueLu::Amesos2Smoother: using "Superlu"
MueLu::AmesosSmoother: using "Superlu"
Using default factory (MueLu::SmootherFactory{pre = MueLu::DirectSolver{type = }, post = null}) for building 'CoarseSolver'.
Using default factory (MueLu::SmootherFactory{pre = MueLu::TrilinosSmoother{type = RELAXATION}, post = pre}) for building 'Smoother'.
Using default factory (MueLu::AmalgamationFactory) for building 'UnAmalgamationInfo'.
Level 2
 Prolongator smoothing (MueLu::SaPFactory)
  Matrix filtering (MueLu::FilteredAFactory)
   Build (MueLu::CoalesceDropFactory)
    lightweight wrap = 1
    algorithm = "classical": threshold = 0.00, blocksize = 1
    Detected 0 Dirichlet nodes
    Number of dropped entries in unamalgamated matrix graph: 0/14928 (0.00%)
   Filtered matrix is not being constructed as no filtering is being done
  Build (MueLu::TentativePFactory)
   Build (MueLu::UncoupledAggregationFactory)
    Algo "Phase - (Dirichlet)"
     BuildAggregates (Phase - (Dirichlet))
       aggregated : 0 (phase), 0/1700 [0.00%] (total)
       remaining  : 1700
       aggregates : 0 (phase), 0 (total)
    Algo "Phase 1 (main)"
     BuildAggregates (Phase 1 (main))
       aggregated : 1649 (phase), 1649/1700 [97.00%] (total)
       remaining  : 51
       aggregates : 192 (phase), 192 (total)
    Algo "Phase 2 (cleanup)"
     BuildAggregates (Phase 2 (cleanup))
       aggregated : 51 (phase), 1700/1700 [100.00%] (total)
       remaining  : 0
       aggregates : 0 (phase), 192 (total)
    Algo "Phase 3 (emergency)"
     BuildAggregates (Phase 3 (emergency))
       aggregated : 0 (phase), 1700/1700 [100.00%] (total)
       remaining  : 0
       aggregates : 0 (phase), 192 (total)
    Algo "Phase - (isolated)"
     BuildAggregates (Phase - (isolated))
       aggregated : 0 (phase), 1700/1700 [100.00%] (total)
       remaining  : 0
       aggregates : 0 (phase), 192 (total)
    "UC": MueLu::Aggregates{nGlobalAggregates = 192}
   Build (MueLu::AmalgamationFactory)
    AmalagamationFactory::Build(): found fullblocksize=1 and stridedblocksize=1 from strided maps. offset=0
   Build (MueLu::CoarseMapFactory)
    domainGIDOffset: 0 block size: 1 stridedBlockId: -1
   Column map is consistent with the row map, good.
   TentativePFactory : aggregates do not cross process boundaries
   Ptent size =  1700 x 192, nnz = 1700
   Ptent Load balancing info
   Ptent   # active processes: 1/1
   Ptent   # rows per proc   : avg = 1.70e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
   Ptent   #  nnz per proc   : avg = 1.70e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
  Eigenvalue estimate
   Calculating max eigenvalue estimate now (max iters = 10)
   Prolongator damping factor = 0.98 (1.33 / 1.35)
  Fused (I-omega*D^{-1} A)*Ptent
  P size =  1700 x 192, nnz = 4551
  P Load balancing info
  P   # active processes: 1/1
  P   # rows per proc   : avg = 1.70e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
  P   #  nnz per proc   : avg = 4.55e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
 Transpose P (MueLu::TransPFactory)
  R size =  192 x 1700, nnz = 4551
  R Load balancing info
  R   # active processes: 1/1
  R   # rows per proc   : avg = 1.92e+02,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
  R   #  nnz per proc   : avg = 4.55e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
 Computing Ac (MueLu::RAPFactory)
  MxM: A x P
   Matrix product nnz per row estimate = 10
  MxM: R x (AP) (explicit)
   Matrix product nnz per row estimate = 28
  Ac size =  192 x 192, nnz = 1674
  Ac Load balancing info
  Ac   # active processes: 1/1
  Ac   # rows per proc   : avg = 1.92e+02,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
  Ac   #  nnz per proc   : avg = 1.67e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
 Setup Smoother (MueLu::Ifpack2Smoother{type = RELAXATION})
  "Ifpack2::Relaxation": {Initialized: true, Computed: true, Type: Symmetric Gauss-Seidel, sweeps: 1, damping factor: 1, Global matrix dimensions: [192, 192], Global nnz: 1674}
MueLu::Amesos2Smoother: using "Superlu"
MueLu::AmesosSmoother: using "Superlu"
Using default factory (MueLu::SmootherFactory{pre = MueLu::DirectSolver{type = }, post = null}) for building 'CoarseSolver'.
Using default factory (MueLu::SmootherFactory{pre = MueLu::TrilinosSmoother{type = RELAXATION}, post = pre}) for building 'Smoother'.
Level 3
 Prolongator smoothing (MueLu::SaPFactory)
  Matrix filtering (MueLu::FilteredAFactory)
   Build (MueLu::CoalesceDropFactory)
    lightweight wrap = 1
    algorithm = "classical": threshold = 0.00, blocksize = 1
    Detected 0 Dirichlet nodes
    Number of dropped entries in unamalgamated matrix graph: 0/1674 (0.00%)
   Filtered matrix is not being constructed as no filtering is being done
  Build (MueLu::TentativePFactory)
   Build (MueLu::UncoupledAggregationFactory)
    Algo "Phase - (Dirichlet)"
     BuildAggregates (Phase - (Dirichlet))
       aggregated : 0 (phase), 0/192 [0.00%] (total)
       remaining  : 192
       aggregates : 0 (phase), 0 (total)
    Algo "Phase 1 (main)"
     BuildAggregates (Phase 1 (main))
       aggregated : 178 (phase), 178/192 [92.71%] (total)
       remaining  : 14
       aggregates : 24 (phase), 24 (total)
    Algo "Phase 2 (cleanup)"
     BuildAggregates (Phase 2 (cleanup))
       aggregated : 14 (phase), 192/192 [100.00%] (total)
       remaining  : 0
       aggregates : 0 (phase), 24 (total)
    Algo "Phase 3 (emergency)"
     BuildAggregates (Phase 3 (emergency))
       aggregated : 0 (phase), 192/192 [100.00%] (total)
       remaining  : 0
       aggregates : 0 (phase), 24 (total)
    Algo "Phase - (isolated)"
     BuildAggregates (Phase - (isolated))
       aggregated : 0 (phase), 192/192 [100.00%] (total)
       remaining  : 0
       aggregates : 0 (phase), 24 (total)
    "UC": MueLu::Aggregates{nGlobalAggregates = 24}
   Build (MueLu::AmalgamationFactory)
    AmalagamationFactory::Build(): found fullblocksize=1 and stridedblocksize=1 from strided maps. offset=0
   Nullspace factory (MueLu::NullspaceFactory)
   Build (MueLu::CoarseMapFactory)
    domainGIDOffset: 0 block size: 1 stridedBlockId: -1
   Column map is consistent with the row map, good.
   TentativePFactory : aggregates do not cross process boundaries
   
   ******* WARNING *******
   Level::Set: unable to store "Nullspace" generated by factory  on level 3, as it has not been requested and no keep flags were set for it
   Ptent size =  192 x 24, nnz = 192
   Ptent Load balancing info
   Ptent   # active processes: 1/1
   Ptent   # rows per proc   : avg = 1.92e+02,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
   Ptent   #  nnz per proc   : avg = 1.92e+02,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
  Eigenvalue estimate
   Calculating max eigenvalue estimate now (max iters = 10)
   Prolongator damping factor = 0.96 (1.33 / 1.38)
  Fused (I-omega*D^{-1} A)*Ptent
  P size =  192 x 24, nnz = 512
  P Load balancing info
  P   # active processes: 1/1
  P   # rows per proc   : avg = 1.92e+02,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
  P   #  nnz per proc   : avg = 5.12e+02,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
 Transpose P (MueLu::TransPFactory)
  R size =  24 x 192, nnz = 512
  R Load balancing info
  R   # active processes: 1/1
  R   # rows per proc   : avg = 2.40e+01,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
  R   #  nnz per proc   : avg = 5.12e+02,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
 Computing Ac (MueLu::RAPFactory)
  MxM: A x P
   Matrix product nnz per row estimate = 10
  MxM: R x (AP) (explicit)
   Matrix product nnz per row estimate = 29
  Ac size =  24 x 24, nnz = 192
  Ac Load balancing info
  Ac   # active processes: 1/1
  Ac   # rows per proc   : avg = 2.40e+01,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
  Ac   #  nnz per proc   : avg = 1.92e+02,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
 Setup Smoother (MueLu::Amesos2Smoother{type = <ignored>})

--------------------------------------------------------------------------------
---                            Multigrid Summary                             ---
--------------------------------------------------------------------------------
Number of levels    = 4
Operator complexity = 1.34

matrix  rows    nnz  nnz/row procs
A 0    10000  49600     4.96  1
A 1     1700  14928     8.78  1
A 2      192   1674     8.72  1
A 3       24    192     8.00  1

Smoother (level 0) both : "Ifpack2::Relaxation": {Initialized: true, Computed: true, Type: Symmetric Gauss-Seidel, sweeps: 1, damping factor: 1, Global matrix dimensions: [10000, 10000], Global nnz: 49600}

Smoother (level 1) both : "Ifpack2::Relaxation": {Initialized: true, Computed: true, Type: Symmetric Gauss-Seidel, sweeps: 1, damping factor: 1, Global matrix dimensions: [1700, 1700], Global nnz: 14928}

Smoother (level 2) both : "Ifpack2::Relaxation": {Initialized: true, Computed: true, Type: Symmetric Gauss-Seidel, sweeps: 1, damping factor: 1, Global matrix dimensions: [192, 192], Global nnz: 1674}

Smoother (level 3) pre  : <Direct> solver interface
Smoother (level 3) post : no smoother

