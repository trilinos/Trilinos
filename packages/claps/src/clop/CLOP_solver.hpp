#include <mpi.h>
#include <math.h>
#include "CLOP_graph.hpp"
#include "CLOP_sub.hpp"
#include "CLOP_constraint.hpp"
#include "CRD_utils.hpp"
//#include "../include/sort_prototypes.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_MpiComm.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
//#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_LocalMap.h"
#include "Epetra_BLAS.h"
#include "Epetra_LAPACK.h"
#include "EpetraExtCD_MatrixMatrix.hpp"

class CLOP_solver 
{
 public: // functions
  CLOP_solver(
     const Epetra_CrsMatrix* AStandard_,   // stiffness matrix
     const Epetra_IntVector* LDStandard_,  // local dof information
     const Epetra_MultiVector* CStandard_, // coordinates of dofs
     const Epetra_Map* SubMap_,            // map for subdomain dofs
     const Epetra_CrsMatrix* ConStandard_, // subdomain constraint matrix
     const Epetra_IntVector* GNStandard_,  // global nodes for subdomain dofs
     const double* clop_params_);          // array of solver parameters
  ~CLOP_solver();
  void solve(Epetra_Vector* uStand, const Epetra_Vector* fStand,
	     int & num_iter, int & solver_status);
  void mpcforces(Epetra_Vector* uLocal, Epetra_Import* ImporterStLam2Loc);
 private: // functions
  void zero_pointers();
  void process_constraints();
  int transform_constraints();
  void construct_transpose(Epetra_CrsMatrix* & A, Epetra_CrsMatrix* & AT);
  void construct_Overlap();
  void construct_subdomains();
  int initialize_subdomains();
  void calculate_coarse_stiff();
  void gather_coarse_stiff();
  void factor_coarse_stiff();
  void solve_init();
  void construct_Overlap_Subs(Epetra_CrsGraph* & Overlap_Subs);
  void flag_sub_bound(unsigned char on_sub_bound[], unsigned char nsubdof[], 
		      Epetra_CrsGraph* & Overlap_Sub);
  void correct_shape_ela(int rbm, double Edof[], Epetra_Vector & Edof_Overlap,
		      Epetra_Vector & Edof_Standard, unsigned char nsubdof[]);
  void correct_shape_dkt(int rbm, double Edof[], Epetra_Vector & Edof_Overlap,
		      Epetra_Vector & Edof_Standard, unsigned char nsubdof[]);
  void assemble_Phi();
  void pcg_solve(Epetra_Vector* uStand, const Epetra_Vector* fStand,
		 int & num_iter, int & pcg_status);
  void gmres_solve(Epetra_Vector* uStand, const Epetra_Vector* fStand, 
		   int & num_iter, int & gmres_status);
  void apply_preconditioner(Epetra_Vector* r, Epetra_Vector* z);
  void search_correction(Epetra_Vector *r, Epetra_Vector *u, double BETA2);
  void coarse_correction(const Epetra_Vector* r, Epetra_Vector* u);
  void remove_projection_search(Epetra_Vector* v);
  void store_search(double pAp);
  void mat_vec_prod(Epetra_Vector* u, Epetra_Vector* Au);
  void calculate_multipliers(Epetra_Vector* uStand, double & norm_rconstraint, 
			     double & norm_conerror);
  void calculate_condition(int miter);
  void two_steps_CGS(int gmres_iter, Epetra_Vector* r);
  void hessenberg_qr(int gmres_iter);
  void gmres_givens(double a, double b, double & c, double & s);
  void construct_solution(int gmres_iter, double normb);
  void spmat_datfile(const Epetra_CrsMatrix & A, char fname[], int opt);

 private: // variables
  const Epetra_CrsMatrix *AStandard, *ConStandard;
  const Epetra_IntVector *LDStandard, *GNStandard;
  const Epetra_MultiVector* CStandard;
  const Epetra_Map* SubMap;
  const double* clop_params;
  const Epetra_Comm & Comm;

  Epetra_CrsMatrix *AOverlap, *Tran, *ASt_red_keep;
  const Epetra_CrsMatrix *ASt_red;
  Epetra_IntVector *LDOverlap, *GNOverlap;
  Epetra_MultiVector *COverlap, *AP_matrix, *P_matrix;
  Epetra_Import *ImporterST2O, *Importer_coarse;
  Epetra_Export *ExporterO2ST, *ExporterSub2ST, *Exporter_lam;
  Epetra_CrsMatrix *Phi, *PhiT, *Kc, *Kc_gathered, *CtT;
  Epetra_Vector *PhiTr, *PhiTr_gathered, *CSol_gathered, *Lambda_local;
  Epetra_Vector *Lambda, *gSt_red, *ConError;
  Epetra_Vector *rSt_red, *ApSt_red, *zSt_red, *rOverlap, *zOverlap;
  Epetra_Vector *pSt_red, *uSt_red, *vSt_red, *vStand, *wStand;
  Epetra_Map *RowMap_coarse, *RowMapMyCon;
  Epetra_BLAS EB;
  Epetra_LAPACK EL;
  MPI_Comm mpicomm;
  int overlap, maxiter, atype, ndim, local_solver, prt_debug, prt_summary;
  int max_orthog;
  double solver_tol;
  int npart, ndof, ndof_overlap, MyPID, NumProc, gpart0, gcdof0;
  int ndof_rot, *count1, *cs_local, *csdima, ncdof_proc, max_csdim;
  int *dofpart1, *dofpart2, *imap, *cdof_proc, ncdof, max_ndof;
  int ndof_Standard, *sub_gdofs, nsub_gdofs, ncon_global, nx2, nx2_global;
  int *x2_dof, nmycon, *mycdof, gmres_flag, n_orthog_used;
  int pre_type_orthog, pre_type_coarse, orthog_option, *IPIV, ndof_global;
  int ndof_global_red, print_flag;
  double *xcent, *ycent, *zcent, *sol_coarse, *temp_coarse, *rhs_coarse;
  double *rcurra, *r_overlap, *z_overlap, *rhs_work, *sol_work, *tmp_work;
  double *rhoa, *betaa, *pApa, *Etri, *Dtri, *econa, *lambda_local;
  double *lambda, *ortho_vec, *pAp_vec, *ortho_sum, *PAP, *PAP_sum, *PAP_store;
  double *VV, *HH, *RR, *zz, *cc, *ss, *norms, *gmres_vec, *gmres_sum;
  ofstream fout;
  CLOP_sub *Asub;
  sparse_lu *Kc_fac;
};
