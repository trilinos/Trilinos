verbosity = high
coarse: max size = 100
transpose: use implicit = 1
max levels = 2
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
 Setup Smoother (MueLu::IfpackSmoother{type = point relaxation stand-alone})
  IFPACK (Local SGS, sweeps=1, damping=1)
MueLu::Amesos2Smoother: using "Superlu"
MueLu::AmesosSmoother: using "Superlu"
Using default factory (MueLu::SmootherFactory{pre = MueLu::DirectSolver{type = }, post = null}) for building 'CoarseSolver'.
Using default factory (MueLu::SmootherFactory{pre = MueLu::TrilinosSmoother{type = RELAXATION}, post = pre}) for building 'Smoother'.
Level 1
 Computing Ac (MueLu::RAPFactory)
  MxM: A x P
  MxM: P' x (AP) (implicit)
  Ac size =  1700 x 1700, nnz = 14928
  Ac Load balancing info
  Ac   # active processes: 1/1
  Ac   # rows per proc   : avg = 1.70e+03,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
  Ac   #  nnz per proc   : avg = 1.49e+04,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
 Setup Smoother (MueLu::AmesosSmoother{type = <ignored>})

--------------------------------------------------------------------------------
---                            Multigrid Summary                             ---
--------------------------------------------------------------------------------
Number of levels    = 2
Operator complexity = 1.30

matrix  rows    nnz  nnz/row procs
A 0    10000  49600     4.96  1
A 1     1700  14928     8.78  1

Smoother (level 0) both : IFPACK (Local SGS, sweeps=1, damping=1)

Smoother (level 1) pre  : MueLu::AmesosSmoother{type = <ignored>}
Smoother (level 1) post : no smoother

