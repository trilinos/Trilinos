Level 0
 Setup Smoother (MueLu::IfpackSmoother{type = Chebyshev})
 chebyshev: degree = 2
 chebyshev: ratio eigenvalue = 20
 chebyshev: min eigenvalue = 1
 chebyshev: zero starting solution = 1
 chebyshev: eigenvalue max iterations = 10
 
Level 1
 Build (MueLu::RebalanceTransferFactory)
  Build (MueLu::RepartitionFactory)
   Computing Ac (MueLu::RAPFactory)
    Prolongator smoothing (MueLu::SaPFactory)
     Matrix filtering (MueLu::FilteredAFactory)
      Build (MueLu::CoalesceDropFactory)
      aggregation: drop tol = 0.02
      aggregation: Dirichlet threshold = 0   [default]
      aggregation: drop scheme = distance laplacian
      lightweight wrap = 1
      
     filtered matrix: use lumping = 1
     filtered matrix: reuse graph = 1   [default]
     filtered matrix: reuse eigenvalue = 1   [default]
     
     Build (MueLu::TentativePFactory)
      Build (MueLu::UncoupledAggregationFactory)
      aggregation: mode = old   [default]
      aggregation: max agg size = -1   [default]
      aggregation: min agg size = 2   [default]
      aggregation: max selected neighbors = 0   [default]
      aggregation: ordering = natural   [default]
      aggregation: enable phase 1 = 1   [default]
      aggregation: enable phase 2a = 1   [default]
      aggregation: enable phase 2b = 1   [default]
      aggregation: enable phase 3 = 1   [default]
      aggregation: preserve Dirichlet points = 0   [default]
      UseOnePtAggregationAlgorithm = 0   [default]
      UsePreserveDirichletAggregationAlgorithm = 0   [default]
      UseUncoupledAggregationAlgorithm = 1   [default]
      UseMaxLinkAggregationAlgorithm = 1   [default]
      UseIsolatedNodeAggregationAlgorithm = 1   [default]
      UseEmergencyAggregationAlgorithm = 1   [default]
      OnePt aggregate map name =    [default]
      
      Build (MueLu::AmalgamationFactory)
      [empty list]
      
      Nullspace factory (MueLu::NullspaceFactory)
      Fine level nullspace = Nullspace
      
      Build (MueLu::CoarseMapFactory)
      Striding info = {}   [default]
      Strided block id = -1   [default]
      Domain GID offsets = {0}   [default]
      
     [empty list]
     
    sa: damping factor = 1.33   [default]
    sa: calculate eigenvalue estimate = 0   [default]
    sa: eigenvalue estimate num iterations = 10   [default]
    
    Transpose P (MueLu::TransPFactory)
    [empty list]
    
   Build (MueLu::CoordinatesTransferFactory)
   write start = -1   [default]
   write end = -1   [default]
   
   transpose: use implicit = 0   [default]
   Keep AP Pattern = 0   [default]
   Keep RAP Pattern = 0   [default]
   CheckMainDiagonal = 0   [default]
   RepairMainDiagonal = 0   [default]
   
  repartition: start level = 1
  repartition: min rows per proc = 2000
  repartition: max imbalance = 1.327
  repartition: print partition distribution = 0   [default]
  repartition: remap parts = 1
  repartition: remap num values = 4   [default]
  
 repartition: rebalance P and R = 0   [default]
 transpose: use implicit = 0   [default]
 useSubcomm = 1   [default]
 type = Interpolation
 write start = -1   [default]
 write end = -1   [default]
 
 Build (MueLu::RebalanceTransferFactory)
 repartition: rebalance P and R = 0   [default]
 transpose: use implicit = 0   [default]
 useSubcomm = 1   [default]
 type = Restriction
 write start = -1   [default]
 write end = -1   [default]
 
 Computing Ac (MueLu::RebalanceAcFactory)
 useSubcomm = 1   [default]
 

--------------------------------------------------------------------------------
---                            Multigrid Summary                             ---
--------------------------------------------------------------------------------
Number of levels    = 4
Operator complexity = 1.48

matrix rows    nnz  nnz/row procs
A 0    9999  29995     3.00  4
A 1    3335  10015     3.00  1
A 2    1111   3331     3.00  1
A 3     371   1111     2.99  1

Smoother (level 0) both : MueLu::IfpackSmoother{type = Chebyshev}

Smoother (level 1) both : MueLu::IfpackSmoother{type = Chebyshev}

Smoother (level 2) both : MueLu::IfpackSmoother{type = Chebyshev}

Smoother (level 3) pre  : MueLu::AmesosSmoother{type = Superlu}
Smoother (level 3) post : no smoother

