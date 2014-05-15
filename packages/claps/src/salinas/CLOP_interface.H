//@HEADER
// ************************************************************************
// 
//         Claps: A Collection of Domain Decomposition Preconditioners
//                and Solvers
//         Copyright (2006) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Clark R. Dohrmann (crdohrm@sandia.gov) 
// 
// ************************************************************************
//@HEADER

/* interface for CLOP (CLark's Overlapping Preconditioner) 

   This interface is written assuming one substructure per processor.
   Each substructure can be decomposed into two or more subdomains 
   inside the CLOP solver

   Subdomain Variables Passed:
     Node X,Y,Z
     Element Connectivity
     K in sparse format
     Force Vector (passed in solve)
     list of boundary conditions (Node #, nodal_dof)
     List of Node communicators
     dofs per node
     local to global node map
     resultant displacement (returned in solve)
     parameters and flags

  Processor information  (for future use, but required)
     pointer into subdomains owned on processor
     subdomain number list on this processor
     (for one subdomain/proc, these are {0,1,2,3..n} and {0,1,2,..n-1} )

  CLOP will be used to solve the global equation
    K * displacement = force


  calling CLOP is performed in three calls

    CLOP_init(...);
    CLOP_solve( force, displacement );
    CLOP_cleanup();


Example

  12--------13--------14--------15--------16--------17
  |         |          |         |         |         |
  |         |          |         |         |         |
  |   0     |     1    |   2     |  3      |  4      |
  |         |          |         |         |         |
  |         |          |         |         |         |
  6---------7----------8---------9---------10-------11
  |         |          |         |         |         |
  |         |          |         |         |         |
  |   5     |     6    |   7     |  8      |  9      |
  |         |          |         |         |         |
  0---------1----------2---------3---------4---------5

  assume 4 subdomains.
    0 - contains elements 0 1 2
    1 - contains elements 3 4
    2 - contains elements 5 6 7
    3 - contains elements 8 9

  assume 2 processors
    0 - contains subdomains 0,1  (elems 0:4)
    1 - contains subdomains 2,3  (elems 5:9)

............................................
Example Using Global numbering
      proc=0, subdomain=0               proc=0, sub=1

  12--------13--------14--------15  15--------16--------17
  |         |          |         |   |         |         |
  |         |          |         |   |         |         |
  |   0     |     1    |   2     |   |  3      |  4      |
  |         |          |         |   |         |         |
  |         |          |         |   |         |         |
  6---------7----------8---------9   9---------10-------11

      proc=1, subdomain=2               proc=1, sub=3
  6---------7----------8---------9   9---------10-------11
  |         |          |         |   |         |         |
  |         |          |         |   |         |         |
  |   5     |     6    |   7     |   |  8      |  9      |
  |         |          |         |   |         |         |
  0---------1----------2---------3   3---------4---------5

............................................
Example Using Local numbering
      proc=0, subdomain=0               proc=0, sub=1

  4---------5----------6---------7   3---------4---------5
  |         |          |         |   |         |         |
  |         |          |         |   |         |         |
  |   0     |     1    |   2     |   |    0    |   1     |
  |         |          |         |   |         |         |
  |         |          |         |   |         |         |
  0---------1----------2---------3   0---------1---------2

      proc=1, subdomain=0               proc=1, sub=1
  4---------5----------6---------7   3---------4---------5
  |         |          |         |   |         |         |
  |         |          |         |   |         |         |
  |   0     |     1    |   2     |   |    0    |   1     |
  |         |          |         |   |         |         |
  0---------1----------2---------3   0---------1---------2

...................

**** in what follows, only the data for subdomain 0 are shown
****
****

 Data for subdomain 0
   (note: this contains global nodes 6 7 8 9 12 13 14 15)
   X[]={ 0.0 1.0 2.0 3.0 0.0 1.0 2.0 3.0 }
   Y[]={ 1.0 1.0 1.0 1.0 2.0 2.0 2.0 2.0 }

 Element Connectivity:
     This is determined using local node numbering. Two integer
     vectors are passed. elemptr points to the offset in the 
     elem_conn array. For this example,
    int elem_ptr[]={0 4 8 12 }
    int elem_conn[]={ 0 1 5 4  1 2 6 5   2 3 7 6 }
     the local connectivity is illustrated below
      4------5------6------7
      |      |      |      |
      |      |      |      |
      |      |      |      |
      0------1------2------3

 The stiffness matrix K, will be represented in harwell boeing
   format, the upper triangle (and diagonal) only. Three vectors
   are passed, K_ai, K_aj, K_aa. Fortran indexing is used in K_ai
   and K_aj. K_ai[nrows]=nnz+1
  int nrows, K_ai[nrows+1], K_aj[nnz];
  double K_aa[nnz];
  int neq=nrows;

 Force Vector
  the right hand side vector corresponding to the stiffness matrix, and
  dimensioned to neq. 
  The force can be thought of as composed of body forces (such as 
  those applied by pressures or gravity), and of applied forces.

 A dof map for K and F must also be provided. The standard numbering
  would have node 0 having the first degrees of freedom in the
  map (as determined by dofs_per_node below). The dimension of
  the map is neq. 

 --------------------------------------------------------------------
 boundary conditions:
     int ndofs_deleted
     int bc_node_ids[]
     unsigned short bc_node_dof[]
   both vectors above are dimensioned at "ndofs_deleted". The first
   lists the nodes where constraints are applied. The second, the
   nodal dof (a number between 0:5) where the bc is applied. In our
   example, if the left side were clamped,
     ndofs_deleted = 12
     bc_node_ids[]={ 0 0 0 0 0 0 4 4 4 4 4 4 } (local numbers used)
     bc_node_dof[]={ 0 1 2 3 4 5 0 1 2 3 4 5 }

 --------------------------------------------------------------------


 Node Communicators
    int numads      :       Number of adjacent subdomains
 
    int*           Size numads
    subdad          :       Subdomain adjacent Lists

    int*           Size numads + 1
    com_ptr          :       Pointer array to com_lis
                            com_ptr[numads] marks the end of com_lis

    int*           Size com_ptr[numads]
    com_lis          :       Shared nodes definition. A list of
			     global node numbers that are shared
			     with the neighbors


   Once again take the subdomain 0 as an example:
    numads = 3
    subdad = [1,2,3]
    com_ptr = [0,2,6,7]
    com_lis = [9,15,6,7,8,9,9]


 The dofs/node are determined using a bitwise OR of the active dofs
  in the unconstrained element. For example, if only X translation is
  available "1" is the value. "2" for Y, "4" for Z, 8 for RotX, etc.
  dofmap_on_node[]={63 63 63 63 63 63 63 63} 
      (because these are shells, with 6 dofs/node and 2^6=64)

 local to global_map[]={ 6 7 8 9 12 13 14 15 }


 --------------------------------------------------------------------
 Note: subdomain/processor information is not used by CLOP
 subdomains/processor
   int subs_ptr_proc[nproc+1]={ 0, 2, 4 }
   int  subs_per_proc[]= { 0, 1,     2, 3 }
   int local_subdomain_index
 NOTE: this subdomain information is for future use. Currently the interface
 requires one subdomain/processor. In that case
    subs_ptr_proc={0, 1, 2, .. nproc}
    subs_per_proc={0, 1, 2, .. nproc-1}
    local_subdomain_index=0
 --------------------------------------------------------------------


 Computed result
   double displacement[neq]. returned in call to solve.


 Parameters and Flags
   Because the parameters and flags are subject to change, and to provide
   default values for them, a C++ object containing the class is passed.
   The object is included below. It is expected that this class will change
   in the future, but it is also expected that all such changes can be
   included in this interface file.


 The solver returns:
   0 if successful
   positive integers for warnings
   negative integers for failure


*/
#ifndef SALINAS_CLOP_INTERFACE
#define SALINAS_CLOP_INTERFACE


// provide definition of MPI_Comm
//#ifdef SALINAS_MPI
#include "mpi.h"
//#else

//#ifndef MPI_Comm
//#define MPI_Comm int
//#endif

//#endif

class CLOPParams;
class CLOPValues;
class MpcLocal;  // defined externally

void CLOP_init(const double *X,
	       const double *Y,
	       const double *Z,
	       int nnodes_this_subdomain,
	       int *elem_ptr,
	       int *elem_conn,
	       int num_elem_this_subdomain,
	       int neq,
	       const double *K_aa,
	       const int *K_aj,
	       const int *K_ai,
	       const int *map,
	       int ndofs_deleted,
	       const int *bc_node_ids,
	       int *bc_node_dof,
	       MPI_Comm *communicator,
	       int numads,
	       int* subdad,
	       int* com_ptr,
	       int* com_lis,
	       unsigned short *dofmap_on_node,
	       int *local_to_global_node_map,
	       const int *subs_ptr_proc,
	       const int *subs_per_proc,
	       int local_subdomain_index,
	       CLOPParams* params,
	       int NumberMpc,
	       MpcLocal* MpcVector);

/**********************************************************************/

void CLOP_cleanup();
/**********************************************************************/

int CLOP_solve(double *Force, double *displacement, 
	       CLOPValues &returnvals);

void CLOP_MpcForces( double *cvals);


/**********************************************************************/

class CLOPParams {

public:
  //
  // let's just have Salinas set these variables for now
  //
  int max_iterations;      // maximum number of iterations
  double solver_tol;       // solver tolerance
  double stagnation_tol;   // stagnation tolerance
  int max_N_orthog_vecs;   // maximum number of orthogonalization vectors
  int coarse_solver;       // 1 for direct solver, 2 for 3-level method
  int overlap;             // overlap for additive Schwarz preconditioner
  int local_solver;        // subdomain solver, 1 for direct
  int prt_summary;         // 0 no summary output, 1 print summary information
  int prt_debug;           // 0 no debug output, 1 print debug information
  int krylov_method;       // 0 for pcg, 1 for gmres, 2 for pcg w/out dot
                           // product assertion test
  int scale_option;        // 0 no scaling, 1 use scaling in factorizations
  int num_rigid_mode;      // number of rigid body modes for mode acceleration
                           // approach in modal transient analysis
  int BailOut;             // controls what happens when stagnation or max
                           // iterations occurs during a CLOP solve. default is
                           // to exit if either occurs before convergence, 
                           // and to continue if bailout option is set.
  //
  // other variables
  //
  int solution_type;       // 0 = statics
                           // 1 = eigen
                           // 2 = implicit transient
                           // 3 = direct frequency response
  int atype;               // analysis type
  int ndim;                // spatial dimension
  //
  // constructor
  //
  CLOPParams() {
    max_iterations        = 400;
    solver_tol            = 1.0e-06;
    stagnation_tol        = 1.0e-06;
    max_N_orthog_vecs     = 200;
    coarse_solver         = 1;
    overlap               = 0;
    local_solver          = 1;
    prt_summary           = 1;
    prt_debug             = 0;
    krylov_method         = 0;
    scale_option          = 0;
    num_rigid_mode        = 0;
    BailOut               = 0;
    solution_type         = 0;
    atype                 = 2;
    ndim                  = 3;
    // mpc_method          = 1;
    rdCLOPParams();
  }
  void rdCLOPParams() {};
  void Echo(int flag );

}; // end of CLOPPParams definition


class CLOPValues
{
public:
  int number_iterations;      // number of iterations
  double primal_residual;     // primal residual
  int number_reorthog_used;   // # reorthogonal vectors used
  double memory_use_this_sub; // total memory usage on this subdomain  
  double time_in_CLOP;        // total time spent in CLOP
  int CLOP_status;            // 0=success, 1=stagnation, 
                              // 2=iteration count exceeded, 3=other error

  CLOPValues() 
  { 
    number_iterations = 0;
    primal_residual = 0.0;
    number_reorthog_used = 0;
    memory_use_this_sub = 0;
    time_in_CLOP = 0.0;
    CLOP_status = 0;
  };

}; // end of CLOPValues definition

/**********************************************************************/
#endif
/* revision information
  Created 03/10/2004 by Clark Dohrmann
  Revised 04/19/2006 by Clark Dohrmann
*/ 
