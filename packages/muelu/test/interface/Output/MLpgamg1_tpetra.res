========================= Aggregate option summary  =========================
min Nodes per aggregate :               1
min # of root nbrs already aggregated : 0
aggregate ordering :                    natural
=============================================================================
Clearing old data (if any)
Using default factory (MueLu::AmalgamationFactory) for building 'UnAmalgamationInfo'.
Using default factory (MueLu::CoarseMapFactory) for building 'CoarseMap'.
Level 0
 Setup Smoother (MueLu::Ifpack2Smoother{type = RELAXATION})
  "Ifpack2::Relaxation": {Initialized: true, Computed: true, Type: Symmetric Gauss-Seidel, sweeps: 2, damping factor: 1, Global matrix dimensions: [9999, 9999], Global nnz: 29995}
Using default factory (MueLu::AmalgamationFactory) for building 'UnAmalgamationInfo'.
Using default factory (MueLu::CoarseMapFactory) for building 'CoarseMap'.
Level 1
 Prolongator smoothing (PG-AMG) (MueLu::PgPFactory)
  Build (MueLu::TentativePFactory)
   Build (MueLu::UncoupledAggregationFactory)
    Build (MueLu::CoalesceDropFactory)
     lightweight wrap = 0
     CoalesceDropFactory::Build(): found blockdim=1 from strided maps. offset=0
     Build (MueLu::AmalgamationFactory)
      AmalagamationFactory::Build(): found fullblocksize=1 and stridedblocksize=1 from strided maps. offset=0
     CoalesceDropFactory::SetupAmalgamationData() # of amalgamated blocks=9999
     CoalesceDropFactory: nodeMap 9999/9999 elements
     Detected 0 Dirichlet nodes
    Algo "Phase - (Dirichlet)"
     BuildAggregates (Phase - (Dirichlet))
       aggregated : 0 (phase), 0/9999 [0.00%] (total)
       remaining  : 9999
       aggregates : 0 (phase), 0 (total)
    Algo "Phase 1 (main)"
     BuildAggregates (Phase 1 (main))
       aggregated : 9998 (phase), 9998/9999 [99.99%] (total)
       remaining  : 1
       aggregates : 3333 (phase), 3333 (total)
    Algo "Phase 2 (cleanup)"
     BuildAggregates (Phase 2 (cleanup))
       aggregated : 1 (phase), 9999/9999 [100.00%] (total)
       remaining  : 0
       aggregates : 0 (phase), 3333 (total)
    Algo "Phase 3 (emergency)"
     BuildAggregates (Phase 3 (emergency))
       aggregated : 0 (phase), 9999/9999 [100.00%] (total)
       remaining  : 0
       aggregates : 0 (phase), 3333 (total)
    Algo "Phase - (isolated)"
     BuildAggregates (Phase - (isolated))
       aggregated : 0 (phase), 9999/9999 [100.00%] (total)
       remaining  : 0
       aggregates : 0 (phase), 3333 (total)
    "UC": MueLu::Aggregates{nGlobalAggregates = 3333}
   Build (MueLu::CoarseMapFactory)
    domainGIDOffset: 0 block size: 1 stridedBlockId: -1
   Column map is consistent with the row map, good.
   TentativePFactory : aggregates do not cross process boundaries
   Ptent size =  9999 x 3333, nnz = 9999
   Ptent Load balancing info
   Ptent   # active processes: 1/1
   Ptent   # rows per proc   : avg = 1.00e+04,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
   Ptent   #  nnz per proc   : avg = 1.00e+04,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
  Utils::Multiply : Estimate for nnz per row of product matrix = 3
  PgPFactory::ComputeRowBasedOmega (MueLu::PgPFactory)
   Utils::Multiply : Estimate for nnz per row of product matrix = 3
   PgPFactory: smoothed aggregation (scheme: DinvAnorm)
   Damping parameter: min = 0.57, max = 0.80
   # negative omegas: 0 out of 3333 column-based omegas
   # NaNs: 0 out of 3333 column-based omegas
  Utils::TwoMatrixAdd : special case detected (one matrix has a fixed nnz per row), using static profiling
  P size =  9999 x 3333, nnz = 16663
  P Load balancing info
  P   # active processes: 1/1
  P   # rows per proc   : avg = 1.00e+04,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
  P   #  nnz per proc   : avg = 1.67e+04,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
 Call prolongator factory for calculating restrictor (MueLu::GenericRFactory)
  Prolongator smoothing (PG-AMG) (MueLu::PgPFactory)
   Transpose A
   Utils::Multiply : Estimate for nnz per row of product matrix = 3
   PgPFactory::ComputeRowBasedOmega (MueLu::PgPFactory)
    Utils::Multiply : Estimate for nnz per row of product matrix = 3
    PgPFactory: smoothed aggregation (scheme: DinvAnorm)
    Damping parameter: min = 0.57, max = 0.80
    # negative omegas: 0 out of 3333 column-based omegas
    # NaNs: 0 out of 3333 column-based omegas
   Utils::TwoMatrixAdd : special case detected (one matrix has a fixed nnz per row), using static profiling
   P size =  3333 x 9999, nnz = 16663
   P Load balancing info
   P   # active processes: 1/1
   P   # rows per proc   : avg = 3.33e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
   P   #  nnz per proc   : avg = 1.67e+04,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
 Computing Ac (MueLu::RAPFactory)
  MxM: A x P
   Utils::Multiply : Estimate for nnz per row of product matrix = 3
  MxM: R x (AP) (explicit)
   Utils::Multiply : Estimate for nnz per row of product matrix = 5
  Ac size =  3333 x 3333, nnz = 9997
  Ac Load balancing info
  Ac   # active processes: 1/1
  Ac   # rows per proc   : avg = 3.33e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
  Ac   #  nnz per proc   : avg = 1.00e+04,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
 Setup Smoother (MueLu::Ifpack2Smoother{type = RELAXATION})
  "Ifpack2::Relaxation": {Initialized: true, Computed: true, Type: Symmetric Gauss-Seidel, sweeps: 2, damping factor: 1, Global matrix dimensions: [3333, 3333], Global nnz: 9997}
Using default factory (MueLu::AmalgamationFactory) for building 'UnAmalgamationInfo'.
Using default factory (MueLu::CoarseMapFactory) for building 'CoarseMap'.
Level 2
 Prolongator smoothing (PG-AMG) (MueLu::PgPFactory)
  Build (MueLu::TentativePFactory)
   Build (MueLu::UncoupledAggregationFactory)
    Build (MueLu::CoalesceDropFactory)
     lightweight wrap = 0
     CoalesceDropFactory::Build(): found blockdim=1 from strided maps. offset=0
     Build (MueLu::AmalgamationFactory)
      AmalagamationFactory::Build(): found fullblocksize=1 and stridedblocksize=1 from strided maps. offset=0
     CoalesceDropFactory::SetupAmalgamationData() # of amalgamated blocks=3333
     CoalesceDropFactory: nodeMap 3333/3333 elements
     Detected 0 Dirichlet nodes
    Algo "Phase - (Dirichlet)"
     BuildAggregates (Phase - (Dirichlet))
       aggregated : 0 (phase), 0/3333 [0.00%] (total)
       remaining  : 3333
       aggregates : 0 (phase), 0 (total)
    Algo "Phase 1 (main)"
     BuildAggregates (Phase 1 (main))
       aggregated : 3332 (phase), 3332/3333 [99.97%] (total)
       remaining  : 1
       aggregates : 1111 (phase), 1111 (total)
    Algo "Phase 2 (cleanup)"
     BuildAggregates (Phase 2 (cleanup))
       aggregated : 1 (phase), 3333/3333 [100.00%] (total)
       remaining  : 0
       aggregates : 0 (phase), 1111 (total)
    Algo "Phase 3 (emergency)"
     BuildAggregates (Phase 3 (emergency))
       aggregated : 0 (phase), 3333/3333 [100.00%] (total)
       remaining  : 0
       aggregates : 0 (phase), 1111 (total)
    Algo "Phase - (isolated)"
     BuildAggregates (Phase - (isolated))
       aggregated : 0 (phase), 3333/3333 [100.00%] (total)
       remaining  : 0
       aggregates : 0 (phase), 1111 (total)
    "UC": MueLu::Aggregates{nGlobalAggregates = 1111}
   Nullspace factory (MueLu::NullspaceFactory)
   Build (MueLu::CoarseMapFactory)
    domainGIDOffset: 0 block size: 1 stridedBlockId: -1
   Column map is consistent with the row map, good.
   TentativePFactory : aggregates do not cross process boundaries
   Ptent size =  3333 x 1111, nnz = 3333
   Ptent Load balancing info
   Ptent   # active processes: 1/1
   Ptent   # rows per proc   : avg = 3.33e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
   Ptent   #  nnz per proc   : avg = 3.33e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
  Utils::Multiply : Estimate for nnz per row of product matrix = 3
  PgPFactory::ComputeRowBasedOmega (MueLu::PgPFactory)
   Utils::Multiply : Estimate for nnz per row of product matrix = 3
   PgPFactory: smoothed aggregation (scheme: DinvAnorm)
   Damping parameter: min = 0.57, max = 0.78
   # negative omegas: 0 out of 1111 column-based omegas
   # NaNs: 0 out of 1111 column-based omegas
  Utils::TwoMatrixAdd : special case detected (one matrix has a fixed nnz per row), using static profiling
  P size =  3333 x 1111, nnz = 5553
  P Load balancing info
  P   # active processes: 1/1
  P   # rows per proc   : avg = 3.33e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
  P   #  nnz per proc   : avg = 5.55e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
 Call prolongator factory for calculating restrictor (MueLu::GenericRFactory)
  Prolongator smoothing (PG-AMG) (MueLu::PgPFactory)
   Transpose A
   Utils::Multiply : Estimate for nnz per row of product matrix = 3
   PgPFactory::ComputeRowBasedOmega (MueLu::PgPFactory)
    Utils::Multiply : Estimate for nnz per row of product matrix = 3
    PgPFactory: smoothed aggregation (scheme: DinvAnorm)
    Damping parameter: min = 0.57, max = 0.78
    # negative omegas: 0 out of 1111 column-based omegas
    # NaNs: 0 out of 1111 column-based omegas
   Utils::TwoMatrixAdd : special case detected (one matrix has a fixed nnz per row), using static profiling
   P size =  1111 x 3333, nnz = 5553
   P Load balancing info
   P   # active processes: 1/1
   P   # rows per proc   : avg = 1.11e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
   P   #  nnz per proc   : avg = 5.55e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
 Computing Ac (MueLu::RAPFactory)
  MxM: A x P
   Utils::Multiply : Estimate for nnz per row of product matrix = 3
  MxM: R x (AP) (explicit)
   Utils::Multiply : Estimate for nnz per row of product matrix = 5
  Ac size =  1111 x 1111, nnz = 3331
  Ac Load balancing info
  Ac   # active processes: 1/1
  Ac   # rows per proc   : avg = 1.11e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
  Ac   #  nnz per proc   : avg = 3.33e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
 Setup Smoother (MueLu::Ifpack2Smoother{type = RELAXATION})
  "Ifpack2::Relaxation": {Initialized: true, Computed: true, Type: Symmetric Gauss-Seidel, sweeps: 2, damping factor: 1, Global matrix dimensions: [1111, 1111], Global nnz: 3331}
Using default factory (MueLu::AmalgamationFactory) for building 'UnAmalgamationInfo'.
Using default factory (MueLu::CoarseMapFactory) for building 'CoarseMap'.
Level 3
 Prolongator smoothing (PG-AMG) (MueLu::PgPFactory)
  Build (MueLu::TentativePFactory)
   Build (MueLu::UncoupledAggregationFactory)
    Build (MueLu::CoalesceDropFactory)
     lightweight wrap = 0
     CoalesceDropFactory::Build(): found blockdim=1 from strided maps. offset=0
     Build (MueLu::AmalgamationFactory)
      AmalagamationFactory::Build(): found fullblocksize=1 and stridedblocksize=1 from strided maps. offset=0
     CoalesceDropFactory::SetupAmalgamationData() # of amalgamated blocks=1111
     CoalesceDropFactory: nodeMap 1111/1111 elements
     Detected 0 Dirichlet nodes
    Algo "Phase - (Dirichlet)"
     BuildAggregates (Phase - (Dirichlet))
       aggregated : 0 (phase), 0/1111 [0.00%] (total)
       remaining  : 1111
       aggregates : 0 (phase), 0 (total)
    Algo "Phase 1 (main)"
     BuildAggregates (Phase 1 (main))
       aggregated : 1111 (phase), 1111/1111 [100.00%] (total)
       remaining  : 0
       aggregates : 371 (phase), 371 (total)
    Algo "Phase 2 (cleanup)"
     BuildAggregates (Phase 2 (cleanup))
       aggregated : 0 (phase), 1111/1111 [100.00%] (total)
       remaining  : 0
       aggregates : 0 (phase), 371 (total)
    Algo "Phase 3 (emergency)"
     BuildAggregates (Phase 3 (emergency))
       aggregated : 0 (phase), 1111/1111 [100.00%] (total)
       remaining  : 0
       aggregates : 0 (phase), 371 (total)
    Algo "Phase - (isolated)"
     BuildAggregates (Phase - (isolated))
       aggregated : 0 (phase), 1111/1111 [100.00%] (total)
       remaining  : 0
       aggregates : 0 (phase), 371 (total)
    "UC": MueLu::Aggregates{nGlobalAggregates = 371}
   Nullspace factory (MueLu::NullspaceFactory)
   Build (MueLu::CoarseMapFactory)
    domainGIDOffset: 0 block size: 1 stridedBlockId: -1
   Column map is consistent with the row map, good.
   TentativePFactory : aggregates do not cross process boundaries
   Ptent size =  1111 x 371, nnz = 1111
   Ptent Load balancing info
   Ptent   # active processes: 1/1
   Ptent   # rows per proc   : avg = 1.11e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
   Ptent   #  nnz per proc   : avg = 1.11e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
  Utils::Multiply : Estimate for nnz per row of product matrix = 3
  PgPFactory::ComputeRowBasedOmega (MueLu::PgPFactory)
   Utils::Multiply : Estimate for nnz per row of product matrix = 3
   PgPFactory: smoothed aggregation (scheme: DinvAnorm)
   Damping parameter: min = 0.00, max = 0.76
   # negative omegas: 0 out of 371 column-based omegas
   # NaNs: 0 out of 371 column-based omegas
  Utils::TwoMatrixAdd : special case detected (one matrix has a fixed nnz per row), using static profiling
  P size =  1111 x 371, nnz = 1851
  P Load balancing info
  P   # active processes: 1/1
  P   # rows per proc   : avg = 1.11e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
  P   #  nnz per proc   : avg = 1.85e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
 Call prolongator factory for calculating restrictor (MueLu::GenericRFactory)
  Prolongator smoothing (PG-AMG) (MueLu::PgPFactory)
   Transpose A
   Utils::Multiply : Estimate for nnz per row of product matrix = 3
   PgPFactory::ComputeRowBasedOmega (MueLu::PgPFactory)
    Utils::Multiply : Estimate for nnz per row of product matrix = 3
    PgPFactory: smoothed aggregation (scheme: DinvAnorm)
    Damping parameter: min = 0.00, max = 0.76
    # negative omegas: 0 out of 371 column-based omegas
    # NaNs: 0 out of 371 column-based omegas
   Utils::TwoMatrixAdd : special case detected (one matrix has a fixed nnz per row), using static profiling
   P size =  371 x 1111, nnz = 1851
   P Load balancing info
   P   # active processes: 1/1
   P   # rows per proc   : avg = 3.71e+02,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
   P   #  nnz per proc   : avg = 1.85e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
 Computing Ac (MueLu::RAPFactory)
  MxM: A x P
   Utils::Multiply : Estimate for nnz per row of product matrix = 3
  MxM: R x (AP) (explicit)
   Utils::Multiply : Estimate for nnz per row of product matrix = 5
  Ac size =  371 x 371, nnz = 1111
  Ac Load balancing info
  Ac   # active processes: 1/1
  Ac   # rows per proc   : avg = 3.71e+02,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
  Ac   #  nnz per proc   : avg = 1.11e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
 Setup Smoother (MueLu::Ifpack2Smoother{type = RELAXATION})
  "Ifpack2::Relaxation": {Initialized: true, Computed: true, Type: Symmetric Gauss-Seidel, sweeps: 2, damping factor: 1, Global matrix dimensions: [371, 371], Global nnz: 1111}
Using default factory (MueLu::AmalgamationFactory) for building 'UnAmalgamationInfo'.
Using default factory (MueLu::CoarseMapFactory) for building 'CoarseMap'.
Level 4
 Prolongator smoothing (PG-AMG) (MueLu::PgPFactory)
  Build (MueLu::TentativePFactory)
   Build (MueLu::UncoupledAggregationFactory)
    Build (MueLu::CoalesceDropFactory)
     lightweight wrap = 0
     CoalesceDropFactory::Build(): found blockdim=1 from strided maps. offset=0
     Build (MueLu::AmalgamationFactory)
      AmalagamationFactory::Build(): found fullblocksize=1 and stridedblocksize=1 from strided maps. offset=0
     CoalesceDropFactory::SetupAmalgamationData() # of amalgamated blocks=371
     CoalesceDropFactory: nodeMap 371/371 elements
     Detected 0 Dirichlet nodes
    Algo "Phase - (Dirichlet)"
     BuildAggregates (Phase - (Dirichlet))
       aggregated : 0 (phase), 0/371 [0.00%] (total)
       remaining  : 371
       aggregates : 0 (phase), 0 (total)
    Algo "Phase 1 (main)"
     BuildAggregates (Phase 1 (main))
       aggregated : 371 (phase), 371/371 [100.00%] (total)
       remaining  : 0
       aggregates : 124 (phase), 124 (total)
    Algo "Phase 2 (cleanup)"
     BuildAggregates (Phase 2 (cleanup))
       aggregated : 0 (phase), 371/371 [100.00%] (total)
       remaining  : 0
       aggregates : 0 (phase), 124 (total)
    Algo "Phase 3 (emergency)"
     BuildAggregates (Phase 3 (emergency))
       aggregated : 0 (phase), 371/371 [100.00%] (total)
       remaining  : 0
       aggregates : 0 (phase), 124 (total)
    Algo "Phase - (isolated)"
     BuildAggregates (Phase - (isolated))
       aggregated : 0 (phase), 371/371 [100.00%] (total)
       remaining  : 0
       aggregates : 0 (phase), 124 (total)
    "UC": MueLu::Aggregates{nGlobalAggregates = 124}
   Nullspace factory (MueLu::NullspaceFactory)
   Build (MueLu::CoarseMapFactory)
    domainGIDOffset: 0 block size: 1 stridedBlockId: -1
   Column map is consistent with the row map, good.
   TentativePFactory : aggregates do not cross process boundaries
   Ptent size =  371 x 124, nnz = 371
   Ptent Load balancing info
   Ptent   # active processes: 1/1
   Ptent   # rows per proc   : avg = 3.71e+02,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
   Ptent   #  nnz per proc   : avg = 3.71e+02,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
  Utils::Multiply : Estimate for nnz per row of product matrix = 3
  PgPFactory::ComputeRowBasedOmega (MueLu::PgPFactory)
   Utils::Multiply : Estimate for nnz per row of product matrix = 3
   PgPFactory: smoothed aggregation (scheme: DinvAnorm)
   Damping parameter: min = 0.57, max = 0.74
   # negative omegas: 0 out of 124 column-based omegas
   # NaNs: 0 out of 124 column-based omegas
  Utils::TwoMatrixAdd : special case detected (one matrix has a fixed nnz per row), using static profiling
  P size =  371 x 124, nnz = 617
  P Load balancing info
  P   # active processes: 1/1
  P   # rows per proc   : avg = 3.71e+02,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
  P   #  nnz per proc   : avg = 6.17e+02,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
 Call prolongator factory for calculating restrictor (MueLu::GenericRFactory)
  Prolongator smoothing (PG-AMG) (MueLu::PgPFactory)
   Transpose A
   Utils::Multiply : Estimate for nnz per row of product matrix = 3
   PgPFactory::ComputeRowBasedOmega (MueLu::PgPFactory)
    Utils::Multiply : Estimate for nnz per row of product matrix = 3
    PgPFactory: smoothed aggregation (scheme: DinvAnorm)
    Damping parameter: min = 0.57, max = 0.74
    # negative omegas: 0 out of 124 column-based omegas
    # NaNs: 0 out of 124 column-based omegas
   Utils::TwoMatrixAdd : special case detected (one matrix has a fixed nnz per row), using static profiling
   P size =  124 x 371, nnz = 617
   P Load balancing info
   P   # active processes: 1/1
   P   # rows per proc   : avg = 1.24e+02,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
   P   #  nnz per proc   : avg = 6.17e+02,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
 Computing Ac (MueLu::RAPFactory)
  MxM: A x P
   Utils::Multiply : Estimate for nnz per row of product matrix = 3
  MxM: R x (AP) (explicit)
   Utils::Multiply : Estimate for nnz per row of product matrix = 6
  Ac size =  124 x 124, nnz = 370
  Ac Load balancing info
  Ac   # active processes: 1/1
  Ac   # rows per proc   : avg = 1.24e+02,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
  Ac   #  nnz per proc   : avg = 3.70e+02,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
 Max coarse size (<= 128) achieved
 Setup Smoother (MueLu::Amesos2Smoother{type = <ignored>})

--------------------------------------------------------------------------------
---                            Multigrid Summary                             ---
--------------------------------------------------------------------------------
Number of levels    = 5
Operator complexity = 1.49

matrix rows    nnz  nnz/row procs
A 0    9999  29995     3.00  1
A 1    3333   9997     3.00  1
A 2    1111   3331     3.00  1
A 3     371   1111     2.99  1
A 4     124    370     2.98  1

Smoother (level 0) both : "Ifpack2::Relaxation": {Initialized: true, Computed: true, Type: Symmetric Gauss-Seidel, sweeps: 2, damping factor: 1, Global matrix dimensions: [9999, 9999], Global nnz: 29995}

Smoother (level 1) both : "Ifpack2::Relaxation": {Initialized: true, Computed: true, Type: Symmetric Gauss-Seidel, sweeps: 2, damping factor: 1, Global matrix dimensions: [3333, 3333], Global nnz: 9997}

Smoother (level 2) both : "Ifpack2::Relaxation": {Initialized: true, Computed: true, Type: Symmetric Gauss-Seidel, sweeps: 2, damping factor: 1, Global matrix dimensions: [1111, 1111], Global nnz: 3331}

Smoother (level 3) both : "Ifpack2::Relaxation": {Initialized: true, Computed: true, Type: Symmetric Gauss-Seidel, sweeps: 2, damping factor: 1, Global matrix dimensions: [371, 371], Global nnz: 1111}

Smoother (level 4) pre  : <Direct> solver interface
Smoother (level 4) post : no smoother

