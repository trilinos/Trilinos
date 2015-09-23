# @HEADER
#
########################################################################
#
#  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
#                  Copyright 2012 Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Questions? Contact Karen Devine	kddevin@sandia.gov
#                    Erik Boman	        egboman@sandia.gov
#
########################################################################
#
# @HEADER
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
