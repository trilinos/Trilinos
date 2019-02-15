
See Trilinos PR#4320
https://github.com/trilinos/Trilinos/pull/4320

Input files generated from driver of @sehereacer 's code

mpirun -np 8 ./anasazi.exe --matrix=Laplacian --normalize --maxiterations=100000 --use1D --numtrials=1 --file=oregon1_010512.mtx --nev=2 --tol=1e-02

and Zoltan2::VectorAdapter::generateFiles() method (PR #4357)

