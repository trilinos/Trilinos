export PYTHONPATH=/home/knliege/local/trilinos/pytrilinos_test/lib/python3.8/site-packages:$PYTHONPATH
export LD_LIBRARY_PATH=/home/knliege/local/trilinos/pytrilinos_test/lib:$LD_LIBRARY_PATH
mpirun -np 4 --map-by ppr:2:socket python CG.py
