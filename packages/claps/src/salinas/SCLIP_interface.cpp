#define SALINAS_MPI
#include <stdio.h>
#include <iostream>
#include <assert.h>
#include <mpi.h>
#include "CLIP_interface.H"
#include "SCLIP_solver.hpp"

static SCLIP_solver *CLIP_ptr = 0;

void CLIP_init(
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
    CLIPParams* params,
    int NumberMpc,
    MpcLocal* MpcVector
)
{
  int cdof_option, max_iter, atype, ndim, max_N_orthog_vecs, coarse_solver;
  int *adj_elem1 = 0, *adj_elem2 = 0, local_solver, prt_debug, prt_summary;
  int chk_sub_singularity, krylov_method, scale_option;
  double solver_tol;
  cdof_option         = params->cdof_option;
  max_iter            = params->max_iterations;
  solver_tol          = params->solver_tol;
  atype               = params->atype;
  ndim                = params->ndim;
  max_N_orthog_vecs   = params->max_N_orthog_vecs;
  coarse_solver       = params->coarse_solver;
  local_solver        = params->local_solver;
  prt_debug           = params->prt_debug;
  prt_summary         = params->prt_summary;
  chk_sub_singularity = params->chk_sub_singularity;
  krylov_method       = params->krylov_method;
  scale_option        = params->scale_option;
  //
  // initialize CLIP_ptr
  //
  //  assert(CLIP_ptr == NULL);
  delete CLIP_ptr;
  CLIP_ptr = new SCLIP_solver(nnode, nelem, neq, nadj_proc, ndofs_deleted, 
		     NumberMpc, E1, E2, adj_proc, H1_global, H2, mpicomm, gnn,
  		     dofmap_on_node, x, y, z, rowbeg_KU, colidx_KU, KU,
  		     mapdof2row, bc_node_ids, bc_node_dof, MpcVector);
  //
  // construct H1 array from H1_global
  //
  CLIP_ptr->construct_H1();
  //
  // determine dof ownership and construct Epetra objects
  //
  CLIP_ptr->determine_ownership();
  //
  // construct full stiffness matrix K_base and construct Epetra objects
  //
  CLIP_ptr->construct_K_base();
  //
  // initialize CLIP solver
  //
  CLIP_ptr->CLIP_solver_init(cdof_option, solver_tol, max_iter, atype, ndim,
			     local_solver, max_N_orthog_vecs, prt_debug,
			     prt_summary, chk_sub_singularity, krylov_method,
			     scale_option);
}

int CLIP_solve( double *f, double *u, CLIPValues &returnvals)
{
  int number_iterations, number_reorthog_used, SCLIP_status, max_added_corner;

  CLIP_ptr->solve(f, u, number_iterations, SCLIP_status, max_added_corner);
  returnvals.number_iterations = number_iterations;
  returnvals.CLIP_status = SCLIP_status;
  returnvals.max_added_corner = max_added_corner;
  return 0;
}

void CLIP_MpcForces( double *cvals)
{
  CLIP_ptr->MpcForces(cvals);
}

void CLIP_cleanup()
{
  delete CLIP_ptr; CLIP_ptr = NULL;
}


