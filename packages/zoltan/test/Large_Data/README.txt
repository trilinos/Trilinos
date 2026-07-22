The tests in this directory create in parallel an arbitrarily large problem, run 
Zoltan on the problem, and report the results.

On Linux systems, a signal handler can print out /proc/meminfo if the test fails,
indicating whether there is a bug in Zoltan or the test, or whether the test
simply ran out of memory.  The line of interest is Committed_AS, how much memory
would be required to satisfy all of the outstanding mallocs.

stressTestRCB
=============
Create a problem, run recursive coordinate bisection, and report the results.

  usage:  mpiexec -np {num_procs} stressTestRCB {num_coords} {dim_weights} {dim_coords}

num_coords - the global number of coordinates that the test will create
dim_weights - the number of weights per coordinate that the test will create
dim_coords - the dimension (1, 2 or 3) of the coordinates.

stressTestRIB
=============
Create a problem, run recursive inertial bisection, and report the results.

  usage:  mpiexec -np {num_procs} stressTestRIB {num_coords} {dim_weights} {dim_coords}

num_coords - the global number of coordinates that the test will create
dim_weights - the number of weights per coordinate that the test will create
dim_coords - the dimension (1, 2 or 3) of the coordinates.

stressTestPHG
=============
Create a problem, run parallel hypergraph partitioning, and report the results.

  usage:  mpiexec -np {num_procs} stressTestPHG {num_vertices} 

num_vertices - the global number of vertices that the test will create and partition
