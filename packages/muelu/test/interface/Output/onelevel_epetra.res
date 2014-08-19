max levels = 1
coarse: type = RELAXATION
verbosity = test
coarse: max size = 2000   [default]
debug: graph level = -1   [default]
number of equations = 1   [default]
transpose: use implicit = 0   [default]
smoother: pre or post = both   [default]
aggregation: type = uncoupled   [default]
multigrid algorithm = sa   [default]
problem: symmetric = 1   [default]
aggregation: export visualization data = 0   [default]
repartition: enable = 0   [default]

Level 0
 Setup Smoother (MueLu::IfpackSmoother{type = point relaxation stand-alone})
 [empty list]


 --------------------------------------------------------------------------------
 ---                            Multigrid Summary                             ---
 --------------------------------------------------------------------------------
 Number of levels    = 1
 Operator complexity = 1.00
 Max Coarse Size     = 2000
 Implicit Transpose  = false

 matrix rows    nnz  nnz/row procs
 A 0    9999  29995     3.00  1

 Smoother (level 0) pre  : MueLu::IfpackSmoother{type = point relaxation stand-alone}
 Smoother (level 0) post : no smoother
