========================= Aggregate option summary  =========================
min Nodes per aggregate :               1
min # of root nbrs already aggregated : 0
aggregate ordering :                    natural
=============================================================================
Level 0
 Setup Smoother (MueLu::IfpackSmoother{type = point relaxation stand-alone})
 relaxation: type = symmetric Gauss-Seidel   [unused]
 relaxation: sweeps = 2   [unused]
 relaxation: damping factor = 1   [unused]

Level 1
 Build (MueLu::RebalanceTransferFactory)
  Build (MueLu::RepartitionFactory)
   Computing Ac (MueLu::RAPFactory)
    Prolongator smoothing (MueLu::SaPFactory)
     Build (MueLu::TentativePFactory)
      Build (MueLu::UncoupledAggregationFactory)
       Build (MueLu::CoalesceDropFactory)
        Build (MueLu::AmalgamationFactory)
        [empty list]

       aggregation: drop tol = 0   [default]
       aggregation: Dirichlet threshold = 0   [default]
       aggregation: drop scheme = classical   [default]
       lightweight wrap = 0   [default]

      aggregation: mode = old   [default]
      aggregation: max agg size = -1   [default]
      aggregation: min agg size = 1   [unused]
      aggregation: max selected neighbors = 0   [unused]
      aggregation: ordering = natural   [unused]
      aggregation: enable phase 1 = 1   [default]
      aggregation: enable phase 2a = 1   [default]
      aggregation: enable phase 2b = 1   [default]
      aggregation: enable phase 3 = 1   [default]
      aggregation: preserve Dirichlet points = 0   [default]
      UseOnePtAggregationAlgorithm = 0   [default]
      UsePreserveDirichletAggregationAlgorithm = 0   [unused]
      UseUncoupledAggregationAlgorithm = 1   [default]
      UseMaxLinkAggregationAlgorithm = 1   [default]
      UseIsolatedNodeAggregationAlgorithm = 1   [default]
      UseEmergencyAggregationAlgorithm = 1   [default]
      OnePt aggregate map name =    [default]

      Build (MueLu::CoarseMapFactory)
      Striding info = {}   [default]
      Strided block id = -1   [default]
      Domain GID offsets = {0}   [default]

     [empty list]

    sa: damping factor = 1.33333
    sa: calculate eigenvalue estimate = 0   [default]
    sa: eigenvalue estimate num iterations = 10   [default]

    Transpose P (MueLu::TransPFactory)
    [empty list]

   transpose: use implicit = 0   [default]
   Keep AP Pattern = 0   [default]
   Keep RAP Pattern = 0   [default]
   CheckMainDiagonal = 0   [default]
   RepairMainDiagonal = 0

  repartition: start level = 2   [default]
  repartition: min rows per proc = 512
  repartition: max imbalance = 1.3
  repartition: print partition distribution = 0   [default]
  repartition: remap parts = 1   [default]
  repartition: remap num values = 4   [default]

 repartition: rebalance P and R = 0   [default]
 transpose: use implicit = 0   [default]
 repartition: use subcommunicators = 1   [default]
 type = Interpolation
 write start = -1   [default]
 write end = -1   [default]

 Build (MueLu::RebalanceTransferFactory)
 repartition: rebalance P and R = 0   [default]
 transpose: use implicit = 0   [default]
 repartition: use subcommunicators = 1   [default]
 type = Restriction
 write start = -1   [default]
 write end = -1   [default]

 Computing Ac (MueLu::RebalanceAcFactory)
 repartition: use subcommunicators = 1   [default]

 Setup Smoother (MueLu::IfpackSmoother{type = point relaxation stand-alone})
 relaxation: type = symmetric Gauss-Seidel   [unused]
 relaxation: sweeps = 2   [unused]
 relaxation: damping factor = 1   [unused]

Level 2
 Build (MueLu::RebalanceTransferFactory)
  Build (MueLu::RepartitionFactory)
   Computing Ac (MueLu::RAPFactory)
    Prolongator smoothing (MueLu::SaPFactory)
     Build (MueLu::TentativePFactory)
      Build (MueLu::UncoupledAggregationFactory)
       Build (MueLu::CoalesceDropFactory)
        Build (MueLu::AmalgamationFactory)
        [empty list]

       aggregation: drop tol = 0   [default]
       aggregation: Dirichlet threshold = 0   [default]
       aggregation: drop scheme = classical   [default]
       lightweight wrap = 0   [default]

      aggregation: mode = old   [default]
      aggregation: max agg size = -1   [default]
      aggregation: min agg size = 1   [unused]
      aggregation: max selected neighbors = 0   [unused]
      aggregation: ordering = natural   [unused]
      aggregation: enable phase 1 = 1   [default]
      aggregation: enable phase 2a = 1   [default]
      aggregation: enable phase 2b = 1   [default]
      aggregation: enable phase 3 = 1   [default]
      aggregation: preserve Dirichlet points = 0   [default]
      UseOnePtAggregationAlgorithm = 0   [default]
      UsePreserveDirichletAggregationAlgorithm = 0   [unused]
      UseUncoupledAggregationAlgorithm = 1   [default]
      UseMaxLinkAggregationAlgorithm = 1   [default]
      UseIsolatedNodeAggregationAlgorithm = 1   [default]
      UseEmergencyAggregationAlgorithm = 1   [default]
      OnePt aggregate map name =    [default]

      Build (MueLu::CoarseMapFactory)
      Striding info = {}   [default]
      Strided block id = -1   [default]
      Domain GID offsets = {0}   [default]

     [empty list]

    sa: damping factor = 1.33333
    sa: calculate eigenvalue estimate = 0   [default]
    sa: eigenvalue estimate num iterations = 10   [default]

    Transpose P (MueLu::TransPFactory)
    [empty list]

   transpose: use implicit = 0   [default]
   Keep AP Pattern = 0   [default]
   Keep RAP Pattern = 0   [default]
   CheckMainDiagonal = 0   [default]
   RepairMainDiagonal = 0

  repartition: start level = 2   [default]
  repartition: min rows per proc = 512
  repartition: max imbalance = 1.3
  repartition: print partition distribution = 0   [default]
  repartition: remap parts = 1   [default]
  repartition: remap num values = 4   [default]

 repartition: rebalance P and R = 0   [default]
 transpose: use implicit = 0   [default]
 repartition: use subcommunicators = 1   [default]
 type = Interpolation
 write start = -1   [default]
 write end = -1   [default]

 Build (MueLu::RebalanceTransferFactory)
 repartition: rebalance P and R = 0   [default]
 transpose: use implicit = 0   [default]
 repartition: use subcommunicators = 1   [default]
 type = Restriction
 write start = -1   [default]
 write end = -1   [default]

 Computing Ac (MueLu::RebalanceAcFactory)
 repartition: use subcommunicators = 1   [default]

 Setup Smoother (MueLu::IfpackSmoother{type = point relaxation stand-alone})
 relaxation: type = symmetric Gauss-Seidel   [unused]
 relaxation: sweeps = 2   [unused]
 relaxation: damping factor = 1   [unused]

Level 3
 Build (MueLu::RebalanceTransferFactory)
  Build (MueLu::RepartitionFactory)
   Computing Ac (MueLu::RAPFactory)
    Prolongator smoothing (MueLu::SaPFactory)
     Build (MueLu::TentativePFactory)
      Build (MueLu::UncoupledAggregationFactory)
       Build (MueLu::CoalesceDropFactory)
        Build (MueLu::AmalgamationFactory)
        [empty list]

       aggregation: drop tol = 0   [default]
       aggregation: Dirichlet threshold = 0   [default]
       aggregation: drop scheme = classical   [default]
       lightweight wrap = 0   [default]

      aggregation: mode = old   [default]
      aggregation: max agg size = -1   [default]
      aggregation: min agg size = 1   [unused]
      aggregation: max selected neighbors = 0   [unused]
      aggregation: ordering = natural   [unused]
      aggregation: enable phase 1 = 1   [default]
      aggregation: enable phase 2a = 1   [default]
      aggregation: enable phase 2b = 1   [default]
      aggregation: enable phase 3 = 1   [default]
      aggregation: preserve Dirichlet points = 0   [default]
      UseOnePtAggregationAlgorithm = 0   [default]
      UsePreserveDirichletAggregationAlgorithm = 0   [unused]
      UseUncoupledAggregationAlgorithm = 1   [default]
      UseMaxLinkAggregationAlgorithm = 1   [default]
      UseIsolatedNodeAggregationAlgorithm = 1   [default]
      UseEmergencyAggregationAlgorithm = 1   [default]
      OnePt aggregate map name =    [default]

      Build (MueLu::CoarseMapFactory)
      Striding info = {}   [default]
      Strided block id = -1   [default]
      Domain GID offsets = {0}   [default]

     [empty list]

    sa: damping factor = 1.33333
    sa: calculate eigenvalue estimate = 0   [default]
    sa: eigenvalue estimate num iterations = 10   [default]

    Transpose P (MueLu::TransPFactory)
    [empty list]

   transpose: use implicit = 0   [default]
   Keep AP Pattern = 0   [default]
   Keep RAP Pattern = 0   [default]
   CheckMainDiagonal = 0   [default]
   RepairMainDiagonal = 0

  repartition: start level = 2   [default]
  repartition: min rows per proc = 512
  repartition: max imbalance = 1.3
  repartition: print partition distribution = 0   [default]
  repartition: remap parts = 1   [default]
  repartition: remap num values = 4   [default]

 repartition: rebalance P and R = 0   [default]
 transpose: use implicit = 0   [default]
 repartition: use subcommunicators = 1   [default]
 type = Interpolation
 write start = -1   [default]
 write end = -1   [default]

 Build (MueLu::RebalanceTransferFactory)
 repartition: rebalance P and R = 0   [default]
 transpose: use implicit = 0   [default]
 repartition: use subcommunicators = 1   [default]
 type = Restriction
 write start = -1   [default]
 write end = -1   [default]

 Computing Ac (MueLu::RebalanceAcFactory)
 repartition: use subcommunicators = 1   [default]

 Setup Smoother (MueLu::IfpackSmoother{type = point relaxation stand-alone})
 relaxation: type = symmetric Gauss-Seidel   [unused]
 relaxation: sweeps = 2   [unused]
 relaxation: damping factor = 1   [unused]

Level 4
 Build (MueLu::RebalanceTransferFactory)
  Build (MueLu::RepartitionFactory)
   Computing Ac (MueLu::RAPFactory)
    Prolongator smoothing (MueLu::SaPFactory)
     Build (MueLu::TentativePFactory)
      Build (MueLu::UncoupledAggregationFactory)
       Build (MueLu::CoalesceDropFactory)
        Build (MueLu::AmalgamationFactory)
        [empty list]

       aggregation: drop tol = 0   [default]
       aggregation: Dirichlet threshold = 0   [default]
       aggregation: drop scheme = classical   [default]
       lightweight wrap = 0   [default]

      aggregation: mode = old   [default]
      aggregation: max agg size = -1   [default]
      aggregation: min agg size = 1   [unused]
      aggregation: max selected neighbors = 0   [unused]
      aggregation: ordering = natural   [unused]
      aggregation: enable phase 1 = 1   [default]
      aggregation: enable phase 2a = 1   [default]
      aggregation: enable phase 2b = 1   [default]
      aggregation: enable phase 3 = 1   [default]
      aggregation: preserve Dirichlet points = 0   [default]
      UseOnePtAggregationAlgorithm = 0   [default]
      UsePreserveDirichletAggregationAlgorithm = 0   [unused]
      UseUncoupledAggregationAlgorithm = 1   [default]
      UseMaxLinkAggregationAlgorithm = 1   [default]
      UseIsolatedNodeAggregationAlgorithm = 1   [default]
      UseEmergencyAggregationAlgorithm = 1   [default]
      OnePt aggregate map name =    [default]

      Build (MueLu::CoarseMapFactory)
      Striding info = {}   [default]
      Strided block id = -1   [default]
      Domain GID offsets = {0}   [default]

     [empty list]

    sa: damping factor = 1.33333
    sa: calculate eigenvalue estimate = 0   [default]
    sa: eigenvalue estimate num iterations = 10   [default]

    Transpose P (MueLu::TransPFactory)
    [empty list]

   transpose: use implicit = 0   [default]
   Keep AP Pattern = 0   [default]
   Keep RAP Pattern = 0   [default]
   CheckMainDiagonal = 0   [default]
   RepairMainDiagonal = 0

  repartition: start level = 2   [default]
  repartition: min rows per proc = 512
  repartition: max imbalance = 1.3
  repartition: print partition distribution = 0   [default]
  repartition: remap parts = 1   [default]
  repartition: remap num values = 4   [default]

 repartition: rebalance P and R = 0   [default]
 transpose: use implicit = 0   [default]
 repartition: use subcommunicators = 1   [default]
 type = Interpolation
 write start = -1   [default]
 write end = -1   [default]

 Build (MueLu::RebalanceTransferFactory)
 repartition: rebalance P and R = 0   [default]
 transpose: use implicit = 0   [default]
 repartition: use subcommunicators = 1   [default]
 type = Restriction
 write start = -1   [default]
 write end = -1   [default]

 Computing Ac (MueLu::RebalanceAcFactory)
 repartition: use subcommunicators = 1   [default]

 Setup Smoother (MueLu::IfpackSmoother{type = point relaxation stand-alone})
 relaxation: type = symmetric Gauss-Seidel   [unused]
 relaxation: sweeps = 2   [unused]
 relaxation: damping factor = 1   [unused]

Level 5
 Build (MueLu::RebalanceTransferFactory)
  Build (MueLu::RepartitionFactory)
   Computing Ac (MueLu::RAPFactory)
    Prolongator smoothing (MueLu::SaPFactory)
     Build (MueLu::TentativePFactory)
      Build (MueLu::UncoupledAggregationFactory)
       Build (MueLu::CoalesceDropFactory)
        Build (MueLu::AmalgamationFactory)
        [empty list]

       aggregation: drop tol = 0   [default]
       aggregation: Dirichlet threshold = 0   [default]
       aggregation: drop scheme = classical   [default]
       lightweight wrap = 0   [default]

      aggregation: mode = old   [default]
      aggregation: max agg size = -1   [default]
      aggregation: min agg size = 1   [unused]
      aggregation: max selected neighbors = 0   [unused]
      aggregation: ordering = natural   [unused]
      aggregation: enable phase 1 = 1   [default]
      aggregation: enable phase 2a = 1   [default]
      aggregation: enable phase 2b = 1   [default]
      aggregation: enable phase 3 = 1   [default]
      aggregation: preserve Dirichlet points = 0   [default]
      UseOnePtAggregationAlgorithm = 0   [default]
      UsePreserveDirichletAggregationAlgorithm = 0   [unused]
      UseUncoupledAggregationAlgorithm = 1   [default]
      UseMaxLinkAggregationAlgorithm = 1   [default]
      UseIsolatedNodeAggregationAlgorithm = 1   [default]
      UseEmergencyAggregationAlgorithm = 1   [default]
      OnePt aggregate map name =    [default]

      Build (MueLu::CoarseMapFactory)
      Striding info = {}   [default]
      Strided block id = -1   [default]
      Domain GID offsets = {0}   [default]

     [empty list]

    sa: damping factor = 1.33333
    sa: calculate eigenvalue estimate = 0   [default]
    sa: eigenvalue estimate num iterations = 10   [default]

    Transpose P (MueLu::TransPFactory)
    [empty list]

   transpose: use implicit = 0   [default]
   Keep AP Pattern = 0   [default]
   Keep RAP Pattern = 0   [default]
   CheckMainDiagonal = 0   [default]
   RepairMainDiagonal = 0

  repartition: start level = 2   [default]
  repartition: min rows per proc = 512
  repartition: max imbalance = 1.3
  repartition: print partition distribution = 0   [default]
  repartition: remap parts = 1   [default]
  repartition: remap num values = 4   [default]

 repartition: rebalance P and R = 0   [default]
 transpose: use implicit = 0   [default]
 repartition: use subcommunicators = 1   [default]
 type = Interpolation
 write start = -1   [default]
 write end = -1   [default]

 Build (MueLu::RebalanceTransferFactory)
 repartition: rebalance P and R = 0   [default]
 transpose: use implicit = 0   [default]
 repartition: use subcommunicators = 1   [default]
 type = Restriction
 write start = -1   [default]
 write end = -1   [default]

 Computing Ac (MueLu::RebalanceAcFactory)
 repartition: use subcommunicators = 1   [default]

 Setup Smoother (MueLu::AmesosSmoother{type = Klu})
 [empty list]


--------------------------------------------------------------------------------
---                            Multigrid Summary                             ---
--------------------------------------------------------------------------------
Number of levels    = 6
Operator complexity = 1.50

matrix rows    nnz  nnz/row procs
A 0    9999  29995     3.00  1
A 1    3333   9997     3.00  1
A 2    1111   3331     3.00  1
A 3     371   1111     2.99  1
A 4     124    370     2.98  1
A 5      42    124     2.95  1

Smoother (level 0) both : MueLu::IfpackSmoother{type = point relaxation stand-alone}

Smoother (level 1) both : MueLu::IfpackSmoother{type = point relaxation stand-alone}

Smoother (level 2) both : MueLu::IfpackSmoother{type = point relaxation stand-alone}

Smoother (level 3) both : MueLu::IfpackSmoother{type = point relaxation stand-alone}

Smoother (level 4) both : MueLu::IfpackSmoother{type = point relaxation stand-alone}

Smoother (level 5) pre  : MueLu::AmesosSmoother{type = Klu}
Smoother (level 5) post : no smoother

