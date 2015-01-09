========================= Aggregate option summary  =========================
min Nodes per aggregate :               1
min # of root nbrs already aggregated : 0
aggregate ordering :                    natural
=============================================================================
Clearing old data (if any)
Level 0
 Setup Smoother (MueLu::IfpackSmoother{type = point relaxation stand-alone})
  IFPACK (Local SGS, sweeps=2, damping=1)
Level 1
 Build (MueLu::RebalanceTransferFactory)
  Build (MueLu::RepartitionFactory)
   Computing Ac (MueLu::RAPFactory)
    Prolongator smoothing (MueLu::SaPFactory)
     Build (MueLu::TentativePFactory)
      Build (MueLu::UncoupledAggregationFactory)
       Build (MueLu::CoalesceDropFactory)
        lightweight wrap = 0
        CoalesceDropFactory::Build(): found blockdim=1 from strided maps. offset=0
        Build (MueLu::AmalgamationFactory)
        CoalesceDropFactory::SetupAmalgamationData() # of amalgamated blocks=9999
        CoalesceDropFactory: nodeMap 9999/9999 elements
        Detected 0 Dirichlet nodes
       BuildAggregates (Phase - (Dirichlet))
         aggregated : 0 (phase), 0/9999 [0.00%] (total)
         remaining  : 9999
         aggregates : 0 (phase), 0 (total)
       BuildAggregates (Phase 1 (main))
         aggregated : 9998 (phase), 9998/9999 [99.99%] (total)
         remaining  : 1
         aggregates : 3333 (phase), 3333 (total)
       BuildAggregates (Phase 2 (cleanup))
         aggregated : 1 (phase), 9999/9999 [100.00%] (total)
         remaining  : 0
         aggregates : 0 (phase), 3333 (total)
       BuildAggregates (Phase 3 (emergency))
         aggregated : 0 (phase), 9999/9999 [100.00%] (total)
         remaining  : 0
         aggregates : 0 (phase), 3333 (total)
       BuildAggregates (Phase - (isolated))
         aggregated : 0 (phase), 9999/9999 [100.00%] (total)
         remaining  : 0
         aggregates : 0 (phase), 3333 (total)
       "UC": MueLu::Aggregates{nGlobalAggregates = 3333}
      Build (MueLu::CoarseMapFactory)
      Ptent size =  9999 x 3333, nnz = 9999
      Ptent Load balancing info
      Ptent   # active processes: 1/1
      Ptent   # rows per proc   : avg = 1.00e+04,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
      Ptent   #  nnz per proc   : avg = 1.00e+04,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
     Calculating max eigenvalue estimate now (max iters = 10)
     Prolongator damping factor = 0.68 (1.33 / 1.95)
     P size =  9999 x 3333, nnz = 16663
     P Load balancing info
     P   # active processes: 1/1
     P   # rows per proc   : avg = 1.00e+04,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
     P   #  nnz per proc   : avg = 1.67e+04,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
    Transpose P (MueLu::TransPFactory)
     R size =  3333 x 9999, nnz = 16663
     R Load balancing info
     R   # active processes: 1/1
     R   # rows per proc   : avg = 3.33e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
     R   #  nnz per proc   : avg = 1.67e+04,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
    Ac size =  3333 x 3333, nnz = 9997
    Ac Load balancing info
    Ac   # active processes: 1/1
    Ac   # rows per proc   : avg = 3.33e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
    Ac   #  nnz per proc   : avg = 1.00e+04,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
   Repartitioning?  NO:
     current level = 1, first level where repartitioning can happen is 2
  Using original prolongator
 Build (MueLu::RebalanceTransferFactory)
  Using original restrictor
 Computing Ac (MueLu::RebalanceAcFactory)
 Setup Smoother (MueLu::IfpackSmoother{type = point relaxation stand-alone})
  IFPACK (Local SGS, sweeps=2, damping=1)
Level 2
 Build (MueLu::RebalanceTransferFactory)
  Build (MueLu::RepartitionFactory)
   Computing Ac (MueLu::RAPFactory)
    Prolongator smoothing (MueLu::SaPFactory)
     Build (MueLu::TentativePFactory)
      Build (MueLu::UncoupledAggregationFactory)
       Build (MueLu::CoalesceDropFactory)
        lightweight wrap = 0
        CoalesceDropFactory::Build(): found blockdim=1 from strided maps. offset=0
        Build (MueLu::AmalgamationFactory)
        CoalesceDropFactory::SetupAmalgamationData() # of amalgamated blocks=3333
        CoalesceDropFactory: nodeMap 3333/3333 elements
        Detected 0 Dirichlet nodes
       BuildAggregates (Phase - (Dirichlet))
         aggregated : 0 (phase), 0/3333 [0.00%] (total)
         remaining  : 3333
         aggregates : 0 (phase), 0 (total)
       BuildAggregates (Phase 1 (main))
         aggregated : 3332 (phase), 3332/3333 [99.97%] (total)
         remaining  : 1
         aggregates : 1111 (phase), 1111 (total)
       BuildAggregates (Phase 2 (cleanup))
         aggregated : 1 (phase), 3333/3333 [100.00%] (total)
         remaining  : 0
         aggregates : 0 (phase), 1111 (total)
       BuildAggregates (Phase 3 (emergency))
         aggregated : 0 (phase), 3333/3333 [100.00%] (total)
         remaining  : 0
         aggregates : 0 (phase), 1111 (total)
       BuildAggregates (Phase - (isolated))
         aggregated : 0 (phase), 3333/3333 [100.00%] (total)
         remaining  : 0
         aggregates : 0 (phase), 1111 (total)
       "UC": MueLu::Aggregates{nGlobalAggregates = 1111}
      Build (MueLu::CoarseMapFactory)
      Ptent size =  3333 x 1111, nnz = 3333
      Ptent Load balancing info
      Ptent   # active processes: 1/1
      Ptent   # rows per proc   : avg = 3.33e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
      Ptent   #  nnz per proc   : avg = 3.33e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
     Calculating max eigenvalue estimate now (max iters = 10)
     Prolongator damping factor = 0.68 (1.33 / 1.95)
     P size =  3333 x 1111, nnz = 5553
     P Load balancing info
     P   # active processes: 1/1
     P   # rows per proc   : avg = 3.33e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
     P   #  nnz per proc   : avg = 5.55e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
    Transpose P (MueLu::TransPFactory)
     R size =  1111 x 3333, nnz = 5553
     R Load balancing info
     R   # active processes: 1/1
     R   # rows per proc   : avg = 1.11e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
     R   #  nnz per proc   : avg = 5.55e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
    Ac size =  1111 x 1111, nnz = 3331
    Ac Load balancing info
    Ac   # active processes: 1/1
    Ac   # rows per proc   : avg = 1.11e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
    Ac   #  nnz per proc   : avg = 3.33e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
   Repartitioning?  NO:
     # processes with rows = 1
  Using original prolongator
 Build (MueLu::RebalanceTransferFactory)
  Using original restrictor
 Computing Ac (MueLu::RebalanceAcFactory)
 Setup Smoother (MueLu::IfpackSmoother{type = point relaxation stand-alone})
  IFPACK (Local SGS, sweeps=2, damping=1)
Level 3
 Build (MueLu::RebalanceTransferFactory)
  Build (MueLu::RepartitionFactory)
   Computing Ac (MueLu::RAPFactory)
    Prolongator smoothing (MueLu::SaPFactory)
     Build (MueLu::TentativePFactory)
      Build (MueLu::UncoupledAggregationFactory)
       Build (MueLu::CoalesceDropFactory)
        lightweight wrap = 0
        CoalesceDropFactory::Build(): found blockdim=1 from strided maps. offset=0
        Build (MueLu::AmalgamationFactory)
        CoalesceDropFactory::SetupAmalgamationData() # of amalgamated blocks=1111
        CoalesceDropFactory: nodeMap 1111/1111 elements
        Detected 0 Dirichlet nodes
       BuildAggregates (Phase - (Dirichlet))
         aggregated : 0 (phase), 0/1111 [0.00%] (total)
         remaining  : 1111
         aggregates : 0 (phase), 0 (total)
       BuildAggregates (Phase 1 (main))
         aggregated : 1111 (phase), 1111/1111 [100.00%] (total)
         remaining  : 0
         aggregates : 371 (phase), 371 (total)
       BuildAggregates (Phase 2 (cleanup))
         aggregated : 0 (phase), 1111/1111 [100.00%] (total)
         remaining  : 0
         aggregates : 0 (phase), 371 (total)
       BuildAggregates (Phase 3 (emergency))
         aggregated : 0 (phase), 1111/1111 [100.00%] (total)
         remaining  : 0
         aggregates : 0 (phase), 371 (total)
       BuildAggregates (Phase - (isolated))
         aggregated : 0 (phase), 1111/1111 [100.00%] (total)
         remaining  : 0
         aggregates : 0 (phase), 371 (total)
       "UC": MueLu::Aggregates{nGlobalAggregates = 371}
      Build (MueLu::CoarseMapFactory)
      Ptent size =  1111 x 371, nnz = 1111
      Ptent Load balancing info
      Ptent   # active processes: 1/1
      Ptent   # rows per proc   : avg = 1.11e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
      Ptent   #  nnz per proc   : avg = 1.11e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
     Calculating max eigenvalue estimate now (max iters = 10)
     Prolongator damping factor = 0.69 (1.33 / 1.95)
     P size =  1111 x 371, nnz = 1851
     P Load balancing info
     P   # active processes: 1/1
     P   # rows per proc   : avg = 1.11e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
     P   #  nnz per proc   : avg = 1.85e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
    Transpose P (MueLu::TransPFactory)
     R size =  371 x 1111, nnz = 1851
     R Load balancing info
     R   # active processes: 1/1
     R   # rows per proc   : avg = 3.71e+02,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
     R   #  nnz per proc   : avg = 1.85e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
    Ac size =  371 x 371, nnz = 1111
    Ac Load balancing info
    Ac   # active processes: 1/1
    Ac   # rows per proc   : avg = 3.71e+02,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
    Ac   #  nnz per proc   : avg = 1.11e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
   Repartitioning?  NO:
     # processes with rows = 1
  Using original prolongator
 Build (MueLu::RebalanceTransferFactory)
  Using original restrictor
 Computing Ac (MueLu::RebalanceAcFactory)
 Setup Smoother (MueLu::IfpackSmoother{type = point relaxation stand-alone})
  IFPACK (Local SGS, sweeps=2, damping=1)
Level 4
 Build (MueLu::RebalanceTransferFactory)
  Build (MueLu::RepartitionFactory)
   Computing Ac (MueLu::RAPFactory)
    Prolongator smoothing (MueLu::SaPFactory)
     Build (MueLu::TentativePFactory)
      Build (MueLu::UncoupledAggregationFactory)
       Build (MueLu::CoalesceDropFactory)
        lightweight wrap = 0
        CoalesceDropFactory::Build(): found blockdim=1 from strided maps. offset=0
        Build (MueLu::AmalgamationFactory)
        CoalesceDropFactory::SetupAmalgamationData() # of amalgamated blocks=371
        CoalesceDropFactory: nodeMap 371/371 elements
        Detected 0 Dirichlet nodes
       BuildAggregates (Phase - (Dirichlet))
         aggregated : 0 (phase), 0/371 [0.00%] (total)
         remaining  : 371
         aggregates : 0 (phase), 0 (total)
       BuildAggregates (Phase 1 (main))
         aggregated : 371 (phase), 371/371 [100.00%] (total)
         remaining  : 0
         aggregates : 124 (phase), 124 (total)
       BuildAggregates (Phase 2 (cleanup))
         aggregated : 0 (phase), 371/371 [100.00%] (total)
         remaining  : 0
         aggregates : 0 (phase), 124 (total)
       BuildAggregates (Phase 3 (emergency))
         aggregated : 0 (phase), 371/371 [100.00%] (total)
         remaining  : 0
         aggregates : 0 (phase), 124 (total)
       BuildAggregates (Phase - (isolated))
         aggregated : 0 (phase), 371/371 [100.00%] (total)
         remaining  : 0
         aggregates : 0 (phase), 124 (total)
       "UC": MueLu::Aggregates{nGlobalAggregates = 124}
      Build (MueLu::CoarseMapFactory)
      Ptent size =  371 x 124, nnz = 371
      Ptent Load balancing info
      Ptent   # active processes: 1/1
      Ptent   # rows per proc   : avg = 3.71e+02,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
      Ptent   #  nnz per proc   : avg = 3.71e+02,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
     Calculating max eigenvalue estimate now (max iters = 10)
     Prolongator damping factor = 0.69 (1.33 / 1.94)
     P size =  371 x 124, nnz = 617
     P Load balancing info
     P   # active processes: 1/1
     P   # rows per proc   : avg = 3.71e+02,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
     P   #  nnz per proc   : avg = 6.17e+02,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
    Transpose P (MueLu::TransPFactory)
     R size =  124 x 371, nnz = 617
     R Load balancing info
     R   # active processes: 1/1
     R   # rows per proc   : avg = 1.24e+02,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
     R   #  nnz per proc   : avg = 6.17e+02,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
    Ac size =  124 x 124, nnz = 370
    Ac Load balancing info
    Ac   # active processes: 1/1
    Ac   # rows per proc   : avg = 1.24e+02,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
    Ac   #  nnz per proc   : avg = 3.70e+02,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
   Repartitioning?  NO:
     # processes with rows = 1
  Using original prolongator
 Build (MueLu::RebalanceTransferFactory)
  Using original restrictor
 Computing Ac (MueLu::RebalanceAcFactory)
 Setup Smoother (MueLu::AmesosSmoother{type = Klu})

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

Smoother (level 0) both : IFPACK (Local SGS, sweeps=2, damping=1)

Smoother (level 1) both : IFPACK (Local SGS, sweeps=2, damping=1)

Smoother (level 2) both : IFPACK (Local SGS, sweeps=2, damping=1)

Smoother (level 3) both : IFPACK (Local SGS, sweeps=2, damping=1)

Smoother (level 4) pre  : MueLu::AmesosSmoother{type = Klu}
Smoother (level 4) post : no smoother

