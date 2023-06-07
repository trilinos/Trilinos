export PYTHONPATH=/ascldap/users/knliege/local/trilinos_all/pytrilinos_test/lib/python3.8/site-packages
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/ascldap/users/knliege/local/trilinos_all/pytrilinos_test/lib:/home/projects/ppc64le-pwr9-nvidia/cuda/10.1.105/lib64
/ascldap/users/projects/ppc64le-pwr9-nvidia/openmpi/4.0.1/gcc/7.2.0/cuda/10.1.105/bin/mpiexec  -n 4 --map-by ppr:2:socket python CG.py &
nvidia-smi -lms &> log.txt