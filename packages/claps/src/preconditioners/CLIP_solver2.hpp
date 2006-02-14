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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Clark R. Dohrmann (crdohrm@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef CLIP_SOLVER2_HPP
#define CLIP_SOLVER2_HPP
#include <mpi.h>
#include <math.h>
#include <algorithm>
//#include "CLIP_graph.hpp"
//#include "CLIP_sub.hpp"
#include "CLOP_constraint.hpp"
#include "CRD_utils.hpp"
#include "sparse_lu2.hpp"
//#include "../include/sort_prototypes.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_MpiComm.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_MultiVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_LocalMap.h"
#include "Epetra_BLAS.h"
#include "Epetra_LAPACK.h"
#include "EpetraExtCD_MatrixMatrix.hpp"
#include "CRS_serial.hpp"
#include "preconditioner_crd.hpp"

//class CLIP_solver2 : public preconditioner_crd
class CLIP_solver2
{
 public: // functions
CLIP_solver2::CLIP_solver2(
	  CRS_serial* A_,              // problem data
	  const Epetra_Map* SubMap_,   // subdomain map
	  const Epetra_Map* OwnMap_,   // unique owner map
	  const double* clip_params_); // array of solver parameters
  ~CLIP_solver2();
  void solve(Epetra_Vector* uStand, const Epetra_Vector* fStand,
	     int & num_iter, int & solver_status, int & max_added_cor);
  void mpcforces(Epetra_Vector* uLocal, Epetra_Import* ImporterStLam2Loc);
  void tie_down_coarse(int n, int rowbeg[], int colidx[], 
		       double vals[], int ne, CLAPS_sparse_lu* & AA,
		       double* & Xvecs);
 private: // functions
  void copy_map(const Epetra_Map A, 
		Epetra_Map* & B);
  void copy_multivector(const Epetra_MultiVector* A,
			const Epetra_Map RowMap,
			Epetra_MultiVector* & B);
  void copy_intvector(const Epetra_IntVector* A,
		      const Epetra_Map RowMap,
		      Epetra_IntVector* & B);
  void construct_transpose(Epetra_CrsMatrix* & A, Epetra_CrsMatrix* & AT);
  void determine_components();
  void determine_components(int A1[], int A2[], int N, 
			    int* & comp1, int* & comp2, int & ncomp);
  void determine_dof_sets();
  void determine_corner_dofs();
  void share_corner_info();
  void modify_dof_sets();
  void factor_sub_matrices();
  void calculate_coarse();
  void pcg_init();
  void get_matrix_diag(CRS_serial* AA, double* & diag);
  int determine_shell(int node);
  bool find_sub(int dof, int sub);
  void check_two_dofs(int ii, int* adof, bool* adof_flag);
  void check_three_dofs(int ii, int* adof, bool* adof_flag);
  void get_adof(int ii, int* adof, int & nadof, bool* adof_flag);
  bool contains(int ii, int dof2);
  int find_dof1(int* adof, int nadof);
  void find_farthest1(int* adof, int nadof, int dof1, int & dof2);
  void find_farthest2(int* adof, int nadof, int & dof1, int & dof2);
  void find_dof1_dof2(int* adof, int nadof, int & dof1, int & dof2);
  int find_max_area(int* adof, int nadof, int dof1, int dof2);
  void gen_matrix(CRS_serial* A, int na, int adof[], 
		  int* & rowbeg, int* & colidx, double* &vals);
  void zero_pointers();
  void construct_subdomains();
  int initialize_subdomains();
  void calculate_coarse_stiff();
  void gather_coarse_stiff();
  void factor_coarse_stiff();
  void pcg_solve(Epetra_Vector* uStand, const Epetra_Vector* fStand,
		 int & num_iter, int & pcg_status);
  void gmres_solve(Epetra_Vector* uStand, const Epetra_Vector* fStand, 
		   int & num_iter, int & gmres_status);
  double stat_cond();
  void remove_orthog_null(Epetra_Vector *vec);
  double initial_update();
  void apply_preconditioner();
  void final_update(Epetra_Vector* uB, Epetra_Vector* rB);
  void remove_search(Epetra_Vector* v);
  void determine_AxB(Epetra_Vector* x, Epetra_Vector* b);
  void store_search(double pAp);
  void calculate_u_St();
  void calculate_Au(Epetra_Vector *a_St, Epetra_Vector *b_St);
  void calculate_AuStand(Epetra_Vector *u, Epetra_Vector *Au);
  void calculate_condition(int miter);
  void two_steps_CGS(int gmres_iter, Epetra_Vector* r);
  void hessenberg_qr(int gmres_iter);
  void gmres_givens(double a, double b, double & c, double & s);
  void construct_solution(int gmres_iter, double normb);
  void spmat_datfile(const Epetra_CrsMatrix & A, char fname[], int opt);

 private: // variables
  CRS_serial *A;
  const Epetra_Map *SubMap, *OwnMap;
  const double* clip_params;
  const Epetra_Comm & Comm;

  int cdof_option, maxiter, atype, ndim, local_solver, prt_debug, prt_summary;
  int max_orthog, chk_sub_singularity, krylov_method, scale_option, sub_solver;
  int num_rigid_mode, coarse_solver;
  double solver_tol;

  // extra to work for now
  const Epetra_IntVector *NodalDofs;
  const Epetra_CrsMatrix *ASub;
  Epetra_Map *StandardMap;
  // current stuff
  Epetra_Import *Importer, *ImporterC, *ImporterB;
  Epetra_Export *Exporter, *ExporterC, *ExporterB;
  Epetra_Map *SubMapC, *SubMapB;
  Epetra_Map *OwnMapC, *OwnMapB;
  int MyPID, NumProc, ndof_sub, ndof_global, print_flag;
  int *comp1, *comp2, ncomp, ncomp_sum, *sub1, *sub2, *dset1, *dset2;
  int ndof_set, *corner_flag, *dof2node, *local_sub;
  bool *owner_flag, *sub_flag;
  double *PhiB;

  // old stuff
  Epetra_Import *ImporterSt2B, *Importer_Kc;
  Epetra_Export *Exporter_lam, *Exporter_Kc;
  Epetra_Import *ImporterStand;
  Epetra_Export *ExporterStand;
  Epetra_Vector *Lambda, *Lambda_local, *ConError, *ApB_Sub, *ApB_St;
  Epetra_Vector *pB_Sub, *pB_St, *zB_Sub, *zB_St, *uB_St, *u_St, *vStand;
  Epetra_Vector *r_St, *r_Sub, *rB_St, *rB_Sub, *work_Sub, *work_Kc;
  Epetra_Vector *rBSub, *workB_St, *work_St, *work_St2, *vec_PhiB, *prod_PhiB;
  Epetra_Vector *vec1_ASub, *vec2_ASub, *wStand, *Rhs_null;
  Epetra_Map *RowMapMyCon, *MapB_Sub, *MapB_St;
  Epetra_CrsMatrix *CtT, *Tran, *Block_Stiff;
  Epetra_MultiVector *Phir_St;
  Epetra_IntVector *NodalDofs_red;
  solver_crd *AR, *AI, *AKc;
  int *mycdof, nmycon, cg_iter, n_orthog_used, ncon_global, max_added_corner;
  int nI, nB, nC, nR, nB_own, *dofI, *dofB, *dofC, *dofR, sub_begin;
  int num_tied_down, *tied_down;
  double *lambda, *lambda_local, *weight, *ARinvCT, *CARinvCT, *lambda_r;
  double *RHS_cg, *SOL_cg, *TEMP_cg, *SOL_Kc, *TEMP_Kc;
  double *rcurra, *rhoa, *betaa, *pApa, *Dtri, *Etri, *econa;
  double *P_ortho, *AP_ortho, *orth1, *orth2, *PhirTPhir;
  double *VV, *RR, *HH, *zz, *cc, *ss, *norms, *gmres_vec, *gmres_sum;
  ofstream fout;
};
#endif // CLIP_SOLVER2_HPP
