#include <stdio.h>
#include <iostream>
#include <assert.h>
#include <mpi.h>
#include "CLOP_interface.H"
#include "SCLOP_solver.hpp"

static SCLOP_solver *CLOP_ptr = 0;

void CLOP_init(
    const double *x,
    const double *y,
    const double *z,
    int nnode,
    int *E2,
    int *E1,
    int nelem,
    int neq,
    const double *KU,
    const int *colidx_KU,
    const int *rowbeg_KU,
    const int *mapdof2row,
    int ndofs_deleted,
    const int *bc_node_ids,
    int *bc_node_dof,
    MPI_Comm *mpicomm,
    int nadj_proc,
    int* adj_proc,
    int* H2,
    int* H1_global,
    unsigned short *dofmap_on_node,
    int *gnn,
    const int *subs_ptr_proc,
    const int *subs_per_proc,
    int local_subdomain_index,
    CLOPParams* params,
    int NumberMpc,
    MpcLocal* MpcVector
)
{
  int overlap, max_iter, atype, ndim, max_N_orthog_vecs, coarse_solver;
  int *adj_elem1 = 0, *adj_elem2 = 0, local_solver;
  double solver_tol;
  overlap           = params->overlap;
  max_iter          = params->max_iterations;
  solver_tol        = params->solver_tol;
  atype             = params->atype;
  ndim              = params->ndim;
  max_N_orthog_vecs = params->max_N_orthog_vecs;
  coarse_solver     = params->coarse_solver;
  local_solver      = params->local_solver;
  //
  // initialize CLOP_ptr
  //
  //  assert(CLOP_ptr == NULL);
  delete CLOP_ptr;
  CLOP_ptr = new SCLOP_solver(nnode, nelem, neq, nadj_proc, ndofs_deleted, 
		     NumberMpc, E1, E2, adj_proc, H1_global, H2, mpicomm, gnn,
  		     dofmap_on_node, x, y, z, rowbeg_KU, colidx_KU, KU,
  		     mapdof2row, bc_node_ids, bc_node_dof, MpcVector);
  //
  // construct H1 array from H1_global
  //
  CLOP_ptr->construct_H1();
  //
  // determine dof ownership and construct Epetra objects
  //
  CLOP_ptr->determine_ownership();
  //
  // construct full stiffness matrix K_base and construct Epetra objects
  //
  CLOP_ptr->construct_K_base();
  //
  // initialize CLOP solver
  //
  CLOP_ptr->CLOP_solver_init(overlap, solver_tol, max_iter, atype, ndim,
			     local_solver, max_N_orthog_vecs);
}

int CLOP_solve( double *f, double *u, CLOPValues &returnvals)
{
  int number_iterations, number_reorthog_used, SCLOP_status;

  CLOP_ptr->solve(f, u, number_iterations, SCLOP_status);
  returnvals.number_iterations = number_iterations;
  returnvals.CLOP_status = SCLOP_status;
  return 0;
}

void CLOP_MpcForces( double *cvals)
{
  CLOP_ptr->MpcForces(cvals);
}

void CLOP_cleanup()
{
  delete CLOP_ptr; CLOP_ptr = NULL;
}


