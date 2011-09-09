/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

#include "ml_common.h"
#include "ml_include.h"
#include "ml_agg_min_energy.h"
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include "float.h"


using namespace std;

// ============ //
// private data //
// ============ //

static int     Dinv_size = -1;
static double* Dinv      = 0;
static bool    SinglePrecision = true; // to save memory

#ifdef ICL
void *ml_void_mem_ptr;
#endif

// ====================================================================== 
// Take each column of Op, compute its 2-norm squared, and store
// in Column2Norm[]. 
//
// NOTE: Column2Norm[] must have enough space to store ghost unknowns.
//
// Note: Must a CSR matrix but it is easy to change for a generic getrow()
//
static void ML_multiply_self_all(ML_Operator* Op, double* Column2Norm)
{
  int n = Op->invec_leng;
  int n_rows = Op->getrow->Nrows;

  for (int i = 0 ; i < n ; ++i) Column2Norm[i] = 0.0;

  struct ML_CSR_MSRdata* data = 0;
  data = (struct ML_CSR_MSRdata*)ML_Get_MyGetrowData(Op);

  int nnz     = (data->rowptr)[n_rows];
  int *bindx  = data->columns;
  double *val = data->values;
  double register dtemp;

  for (int i = 0; i < nnz; i++) {
       dtemp = *val++;
       Column2Norm[*bindx++] += (dtemp*dtemp);
  }

  if (Op->getrow->pre_comm != NULL)
    ML_transposed_exchange_bdry(Column2Norm, Op->getrow->pre_comm, n, 
				Op->comm, ML_ADD);
}

// ====================================================================== 
// Take inner product between column(i) in left and column(i) in right 
// and store in InnerProd[i]. 
// 
// NOTE: InnerProd[] must have enough space to store ghost unknowns in right.
//
// Note: Must be CSR matrices but it is easy to change for a generic getrow()
//
inline static void ML_multiply_all(ML_Operator* left, ML_Operator* right,
                                double* InnerProd)
{
  int n = left->invec_leng;
  int n_rows = left->getrow->Nrows;
  int LeftGhost = 0, RightGhost = 0;
  int *LeftGlobal = NULL, *RightGlobal = NULL;
  int *LeftLocal  = NULL, *RightLocal  = NULL, *NewLeftLocal = NULL;

  
  if (n != right->invec_leng || n_rows != right->getrow->Nrows)
  {
    cerr << "Error: non-comparible operators" << endl;
    cerr << "(file " << __FILE__ << ", line " << __LINE__ << ")" << endl;
    exit(EXIT_FAILURE);
  }
  
  // compute the number of ghosts and make up some global ids

  if ( left->getrow->pre_comm != NULL) {
      ML_CommInfoOP_Compute_TotalRcvLength( left->getrow->pre_comm);
      LeftGhost =  left->getrow->pre_comm->total_rcv_length;
      ML_create_unique_id(left->invec_leng, &LeftGlobal, 
                          left->getrow->pre_comm, left->comm, -1);
   }
  if (right->getrow->pre_comm != NULL) {
      ML_CommInfoOP_Compute_TotalRcvLength(right->getrow->pre_comm);
      RightGhost = right->getrow->pre_comm->total_rcv_length;
      ML_create_unique_id(right->invec_leng, &RightGlobal, 
                          right->getrow->pre_comm, right->comm, -1);
   }
   if (  (LeftGhost > 0) && (RightGhost > 0) ) {

      LeftLocal    = (int *) ML_allocate(sizeof(int)*(n+LeftGhost));
      RightLocal   = (int *) ML_allocate(sizeof(int)*(n+RightGhost));
      NewLeftLocal = (int *) ML_allocate(sizeof(int)*(n+LeftGhost));

      // Sort global ids

      for (int i = 0; i <  n+LeftGhost; i++) LeftLocal[i] = i;
      for (int i = 0; i < n+RightGhost; i++) RightLocal[i] = i;

      ML_az_sort( LeftGlobal,  n+LeftGhost,  LeftLocal, NULL);
      ML_az_sort(RightGlobal, n+RightGhost, RightLocal, NULL);

      // Compute new local ids for the left matrix that correspond
      // to the local ids for the right matrix. If there is a left global
      // id that is not found in the right, assign it the id of n+RightGhost+1.
      // It will never be accessed and this avoids an if. Any index
      // that is not found in both matrices doesn't contribute to final result.

      for (int i = 0; i <  n+LeftGhost; i++) NewLeftLocal[i] = n+RightGhost+1;

      int j = 0;
      for (int i = 0; i <  n+LeftGhost; i++) {
         while ( (j < n+RightGhost) && ( RightGlobal[j] < LeftGlobal[i]) ) j++;
         if (RightGlobal[j] == LeftGlobal[i])
            NewLeftLocal[LeftLocal[i]] = RightLocal[j];
      }
   }

   for (int i = 0 ; i < n+RightGhost ; ++i) InnerProd[i] = 0.0;

   struct ML_CSR_MSRdata* left_data = 0;
   struct ML_CSR_MSRdata* right_data = 0;
   left_data  = (struct ML_CSR_MSRdata*)ML_Get_MyGetrowData(left);
   right_data = (struct ML_CSR_MSRdata*)ML_Get_MyGetrowData(right);

   int*    lrowptr = left_data->rowptr;
   int*    rrowptr = right_data->rowptr;
   double *temp_array = NULL;

   temp_array = (double *) ML_allocate(sizeof(double)*(n+RightGhost+1));
                                // use RightGhost because these indices are
                                // based off right matrix
  for (int i = 0; i < n+RightGhost+1; i++) temp_array[i] = 0.;

  double *rval = right_data->values;
  double *lval = left_data->values;
  int*    rbindx  = right_data->columns;
  int     rpos;
  int*    lbindx1  = left_data->columns;
  int*    lbindx2  = left_data->columns;
  for (int row = 0 ; row < n_rows ; ++row)
  {
    int     llen    = lrowptr[row + 1] - lrowptr[row];
    int     rlen    = rrowptr[row + 1] - rrowptr[row];

    if (NewLeftLocal == NULL) 
        for (int i = 0 ; i < llen ; ++i) temp_array[*lbindx1++] = *lval++;
    else
        for (int i = 0 ; i < llen ; ++i)
               temp_array[NewLeftLocal[*lbindx1++]] = *lval++;

    for (int i = 0 ; i < rlen ; ++i) {
      rpos = *rbindx++;
      InnerProd[rpos] += (*rval++) * temp_array[rpos];
    }

    if (NewLeftLocal == NULL) 
      for (int i = 0; i < llen; ++i) temp_array[*lbindx2++] = 0.;
    else
      for (int i = 0; i < llen; ++i) temp_array[NewLeftLocal[*lbindx2++]] = 0.;
  }
  if (right->getrow->pre_comm != NULL)
     ML_transposed_exchange_bdry(InnerProd, right->getrow->pre_comm, n, 
                                 right->comm, ML_ADD);

  if (temp_array  != NULL) ML_free(temp_array);
  if (NewLeftLocal!= NULL) ML_free(NewLeftLocal);
  if (  LeftLocal != NULL) ML_free(LeftLocal);
  if ( RightLocal != NULL) ML_free(RightLocal);
  if (RightGlobal != NULL) ML_free(RightGlobal);
  if ( LeftGlobal != NULL) ML_free(LeftGlobal);
}
// ====================================================================== 
// generate smooth prolongator by minimizing energy                      
//
// \author Marzio Sala and Ray Tuminaro, 9214
// 
// \date 14-Jul-05
//
//  rst: Marzio, as you can see, I have made billions of changes. 
//  I can't really remember what motivated all of them but here is
//  my memory of issues: 
//      1) I did not agree with how you converted column based omegas
//         to row based omegas. If you do not agree with what I did, we
//         need to talk.
//      2) I added a capability for handling block diagonal scaling. I
//         basically tried to keep what you had for the scalar PDE case and
//         add new stuff to handle a block diagonal.
//      3) I also had some confusion about how the restriction was done.
//         Something seemed funny to me, but now I don't remember. I basically
//         ripped out the old stuff (of course it is in CVS) and decided to
//         just take the same omegas that were used for the prolongator.
//         We should probably put back in some capability to either keep the
//         same omegas (which reduces setup cost) or recompute them. It would
//         be nice if I could remember what I did not like in the previous
//         setup.
//      4) I put in a signficant effort to reduce the run time. A lot of this
//         was motivated by the Premo milestone (that was due in early 
//         Sept.). Originally, I was hoping to show some new results ... but
//         I didn't get them done in time. The modifications for speed include:
//              a) I reworked things to avoid explicit transpose operations
//                 involving Amat. These were taking a huge amount of time on
//                 the Premo matrices (which have billions of nonzeros).
//              b) I put in some funky stuff to drop small coefficients. This
//                 was again motivated by premo. Since I don't really know how
//                 to do this robustly. I have commented this out.
//              c) I started working on an algorithm that computes only a 
//                 subset of the omegas. I believe that this worked in serial
//                 and not in parallel. It is currently turned off.
//      
// ====================================================================== 
int ML_AGG_Gen_Prolongator_MinEnergy(ML *ml,int level, int clevel, void *data)
{
  int         Ncoarse, Nfine, gNfine, gNcoarse;
  ML_Operator **prev_P_tentatives;
  ML_Aggregate * ag = (ML_Aggregate *) data;
  double dropping = ag->minimizing_energy_droptol;

#ifdef ML_TIMING
  double t0;
  t0 =  GetClock();
#endif

  ML_Operator* Amat = &(ml->Amat[level]); //already created and filled
  prev_P_tentatives = ag->P_tentative;    // for keep_P_tentative

  Amat->num_PDEs = ag->num_PDE_eqns;

  Nfine    = Amat->outvec_leng;
  gNfine   = ML_Comm_GsumInt(ml->comm, Nfine);
  ML_Aggregate_Set_CurrentLevel(ag, level);

  // =============================================== //
  // Creates the non-smoothed prolongator, P0, then  //
  // Checks the dimension of the coarse problem.     //
  // Note that I always need P0 to form R.           //
  // =============================================== //

  ML_Operator* P0 = ML_Operator_Create(ml->comm);

  if (prev_P_tentatives == 0) {
    ag->P_tentative = ML_Operator_ArrayCreate(ag->max_levels);
    prev_P_tentatives = ag->P_tentative;
    for (int jj = 0 ; jj < ag->max_levels; jj++) prev_P_tentatives[jj] = 0;
  }

  if ((prev_P_tentatives != 0) && (prev_P_tentatives[clevel] != 0)) {
    P0 = prev_P_tentatives[clevel];
    Ncoarse = P0->invec_leng;
  }
  else {
    Ncoarse  = ML_Aggregate_Coarsen(ag,Amat,&P0,ml->comm);
    prev_P_tentatives[clevel] = P0;
  }

  gNcoarse = ML_Comm_GsumInt( ml->comm, Ncoarse);
  gNcoarse = gNcoarse / P0->num_PDEs;
  gNfine   = gNfine   / Amat->num_PDEs;
  if (gNcoarse == 0 || ((1.0*gNfine) / (1.0*gNcoarse+0.1) < 1.05)) {
    ML_Operator_Destroy(&P0);
    if (prev_P_tentatives != 0) prev_P_tentatives[clevel] = NULL;
    return -1;
  }

  int row_length;
  int allocated = 128;
  int*    bindx = (int    *)  ML_allocate(allocated*sizeof(int   ));
  double* val   = (double *)  ML_allocate(allocated*sizeof(double));
  int n = Amat->getrow->Nrows;

  // ======================================================  //
  // Compute the new energy smoothed prolongator.            //
  // This computation is done a little differently depending //
  // on whether we have a scalar PDE or a PDE system.        //
  //                                                         //
  // scalar PDE: the diagonal of Amat is obtained and        //
  //   inverted. A new ML operator is created (and overwrites//
  //   the old one) containing the implicitly scaled Amat.   //
  //   In this way, we no longer have to think about scaling.//
  // system PDE: implicit scaling is no longer practical.    //
  //   instead the diagonal blocks are computed and inverted.//
  //   Whenever a D^{-1} operation is needed, it is done     //
  //   explicitly in the code.
  // ======================================================  //


  // Begin by computing the inverted matrix diagonal in the  //
  // scalar PDE case.                                        //  

  if (Dinv != 0 || Dinv_size != -1) {
    cerr << "Error: Static data Dinv is not null or Dinv_size is wrong!" << endl;
    cerr << "(file " << __FILE__ << ", line " << __LINE__ << ")" << endl;
    exit(EXIT_FAILURE);
  }

  int AmatNghost = 0;
  if (Amat->getrow->pre_comm != 0)
    AmatNghost = Amat->getrow->pre_comm->total_rcv_length;
  ML_Operator *UnscaledAmat = NULL;
  if (Amat->num_PDEs == 1 || ag->block_scaled_SA == 0) {

    // Use point scaling here instead of block scaling. Incorporate the 
    // point scaling automatically in Amat so we don't need to worry about
    // it any more.                                                       

    Dinv = (double *) ML_allocate(sizeof(double)*(Amat->outvec_leng+AmatNghost +1));
    for (int row = 0; row < Amat->outvec_leng; row++) {
      Dinv[row] = 1.;
      ML_get_matrix_row(Amat, 1, &row, &allocated, &bindx,&val,&row_length,0); 
      for (int j = 0; j < row_length; j++) 
	if ( (bindx[j] == row) && (val[j] != 0.)) Dinv[row] = 1./val[j];
    }
    UnscaledAmat = Amat;
    Amat = ML_Operator_ImplicitlyVScale(UnscaledAmat, Dinv, 0);
  }

  ML_Operator* DinvAmat = NULL;
  ML_Operator *DinvAP0_subset = NULL;
  int *Subset = NULL, Nsubset = 0, Nsubset_tot, NComputedOmegas;
  double *Numerator = NULL, *Denominator = NULL;

  int compress = ag->cheap_minimizing_energy;   
                      // compress = 1 corresponds to only computing
                      // the omegas at a subset of points          
                      // some form of this code works, but I can't 
                      // really remember so it is turned off.      

  NComputedOmegas = P0->invec_leng;

  // Compute D^{-1} A P0

  ML_Operator *DinvAP0 = ML_Operator_Create(P0->comm);
  ML_2matmult(Amat, P0, DinvAP0, ML_CSR_MATRIX);

  // Scale result. Note: if Amat corresponds to a scalar PDE, the
  // point scaling is already incorporated into Amat so there is 
  // only a need to scale explicitly if solving a PDE system.    


  if (Amat->num_PDEs != 1 && ag->block_scaled_SA == 1) {
    DinvAmat = ML_Operator_ImplicitlyBlockDinvScale(Amat);
    ML_Operator_ExplicitDinvA(Amat->num_PDEs,(MLSthing *)DinvAmat->data,
			      DinvAP0);
  }

  if (SinglePrecision)
    ML_Operator_ChangeToSinglePrecision(DinvAP0);

  int AP0Nghost = 0;
  if (DinvAP0->getrow->pre_comm != NULL) {
      ML_CommInfoOP_Compute_TotalRcvLength(DinvAP0->getrow->pre_comm);
      AP0Nghost = DinvAP0->getrow->pre_comm->total_rcv_length;
   }

  // rst: This is a routine that will drop small entries. In Premo
  //      there are billions of small values. I don't really remember
  //      the current state of this code nor do I remember how to choose
  //      good values for the dropping ... so it is commented out.
  //
  if (dropping != 0.0)
    ML_CSR_DropSmall(DinvAP0, 0.0, dropping, 0.0);

  DinvAP0_subset = DinvAP0;

  int n_0     = P0->invec_leng;
  int n_0_tot = ML_gsum_int(n_0, P0->comm);
  if (n_0_tot < 100) compress = 0;
  if (ag->minimizing_energy != 2) compress = 0;
      //
      // We don't have a subset of P0 so compress is only supported for
      // minimizing_energy = 2.
  Nsubset_tot = n_0_tot;

  if (compress == 1) {

    // Only compute omega's at a subset of coarse grid points and then
    // later approximate the omegas at the neighboring points.               

    double *dtemp;
    int    count = 0;

    dtemp    = (double *) ML_allocate((n_0+1+AP0Nghost)*sizeof(double));
    Subset   = (int    *) ML_allocate((n_0+1+AP0Nghost)*sizeof(int));

    ML_random_vec(dtemp, n_0, P0->comm);
    ML_exchange_bdry(dtemp,DinvAP0->getrow->pre_comm, n_0,
                     DinvAP0->comm, ML_OVERWRITE,NULL);
    for (int i = 0; i < n_0+AP0Nghost; i++) {
       if ( ML_dabs(dtemp[i])*((double)Amat->num_PDEs) < .5 ) {
          Subset[i] = count++;
          if (i < n_0 ) Nsubset++;
       }
       else Subset[i] = -1;
    }
    ML_free(dtemp);


    /* Here is some old code that works in serial but I'm not sure in parallel.
     * It can be used as an alternative to random() above. It essentially
     * uses root points of an MIS(P0T_AP0) for the locations of computed omegas.
     *
     * ML_Operator *P0_T = NULL, *P0_TAP0 = NULL;
     * P0_T       = ML_Operator_Create(Amat->comm);
     * P0_TAP0    = ML_Operator_Create(Amat->comm);
     *
     * ML_Operator_Transpose(P0,P0_T);
     * ML_2matmult(P0_T, DinvAP0, P0_TAP0, ML_CSR_MATRIX);
     * Nsubset = ML_Operator_MisRootPts( P0_TAP0,  Amat->num_PDEs, &Subset);
     *
     * if (P0_TAP0)  ML_Operator_Destroy(&P0_TAP0);
     * if (P0_T)     ML_Operator_Destroy(&P0_T);
     */

    Nsubset_tot = ML_gsum_int(Nsubset, P0->comm);
    if (Nsubset_tot == 0) compress = 0;
  }

  if (compress == 1) {

    NComputedOmegas = Nsubset;


    // Form the compressed matrix

    DinvAP0_subset = ML_CSRmatrix_ColumnSubset(DinvAP0, Nsubset,Subset);

    if (dropping != 0.0)
      ML_CSR_DropSmall(DinvAP0_subset, 0.0, dropping, 0.0);

  }

  ML_Operator* DinvADinvAP0 = 0;



  //
  // Compute the matrix D^{-1} A  D^{-1} A P0 which will
  // be needed when the minimization norm corresponds to
  // A + A' and A'A.

  if (ag->minimizing_energy != 1) {
    DinvADinvAP0 = ML_Operator_Create(P0->comm);
    ML_2matmult(Amat, DinvAP0_subset, DinvADinvAP0, ML_CSR_MATRIX);

    // Scale result. Note: if Amat corresponds to a scalar PDE, the 
    // point scaling is already incorporated into Amat so there is  
    // only a need to scale explicitly if solving a PDE system.     
    if (Amat->num_PDEs != 1 && ag->block_scaled_SA == 1) 
      ML_Operator_ExplicitDinvA(Amat->num_PDEs,(MLSthing *)DinvAmat->data,
				DinvADinvAP0);  

    if (dropping != 0.0)
      ML_CSR_DropSmall(DinvADinvAP0, 0.0, dropping, 0.0); 

    if (SinglePrecision)
      ML_Operator_ChangeToSinglePrecision(DinvADinvAP0);
  }
  int MaxGhost;

   MaxGhost = AP0Nghost;
  if ( (DinvADinvAP0 != NULL) && (DinvADinvAP0->getrow->pre_comm != NULL)) {
      ML_CommInfoOP_Compute_TotalRcvLength(DinvADinvAP0->getrow->pre_comm);
      if (DinvADinvAP0->getrow->pre_comm->total_rcv_length > MaxGhost)
      MaxGhost = DinvADinvAP0->getrow->pre_comm->total_rcv_length;
   }
  Numerator   = (double *) ML_allocate(sizeof(double)*(P0->invec_leng+MaxGhost));
  Denominator = (double *) ML_allocate(sizeof(double)*(P0->invec_leng+MaxGhost));
  for (int i = 0 ; i < NComputedOmegas+MaxGhost ; ++i) Numerator[i] = 0.;
  for (int i = 0 ; i < NComputedOmegas+MaxGhost ; ++i) Denominator[i] = 0.;

  vector<double> tmp(n_0+MaxGhost);
  vector<double> ColBasedOmega(n_0+ MaxGhost);

  for (int i = 0 ; i < n_0+MaxGhost ; ++i) tmp[i] = 0.;
  for (int i = 0 ; i < n_0+MaxGhost ; ++i) ColBasedOmega[i] = 0.;

  switch (ag->minimizing_energy) {
  case 1:
    // Minimize with respect to L2 norm

    //                  diag( P0' D^{-1} A P0)
    //   omega =   -----------------------------
    //             diag( P0' A' D^{-1}' D^{-1} A P0)
    //

    ML_multiply_all(P0, DinvAP0, &Numerator[0]);
    ML_multiply_self_all(DinvAP0, &Denominator[0]);
    break;

  case 2:
    // Minimize with respect to the (D^{-1} A)' D^{-1} A norm. 
    // Need to be smart here to avoid the construction of A' A
    //
    //                   diag( P0' ( A' D^{-1}' D^{-1} A) D^{-1} A P0)
    //   omega =   --------------------------------------------------------
    //             diag( P0' A' D^{-1}' ( A' D^{-1}' D^{-1} A) D^{-1} A P0)
    //
    ML_multiply_all(DinvAP0_subset, DinvADinvAP0, &Numerator[0]);
    ML_multiply_self_all(DinvADinvAP0, &Denominator[0]);

    break;

  case 3:
    //             diag( P0' ( A'D' + DA) D A P0)
    //   omega =   -----------------------------
    //             diag( P0'A'D' ( A'D' + DA) D A P0)
    //
    //             diag( DinvAP0'DinvAP0 + P0'DinvADinvAP0)
    //         =   -----------------------------
    //                2*diag( DinvADinvAP0'DinvAP0)
    //
    //
    ML_multiply_all(P0, DinvADinvAP0,  &Numerator[0]);
    ML_multiply_self_all(DinvAP0, &tmp[0]);

    for (int i = 0 ; i < n_0 ; ++i) Numerator[i] += tmp[i];
    ML_multiply_all(DinvAP0, DinvADinvAP0, &Denominator[0]);
    for (int i = 0 ; i < n_0 ; ++i) Denominator[i] *= 2.;
    break;

  default:
    // should never be here
    cerr << "Incorrect parameter (" << ag->minimizing_energy << ")" << endl;
    cerr << "(file " << __FILE__ << ", line " << __LINE__ << ")" << endl;
    exit(EXIT_FAILURE);
  }

  int zero_local   = 0;
  double min_local = DBL_MAX;
  double max_local = DBL_MIN;

  // Compute 'column-based' Omega's (with safeguard)

#ifdef EMIN_IN_PAPER
FILE *fp;
fp = fopen("colomega","w");
#endif
  for (int i = 0 ; i < NComputedOmegas ; ++i) {
    ColBasedOmega[i] = Numerator[i]/Denominator[i];
#ifdef EMIN_IN_PAPER
fprintf(fp,"%20.13e\n",ColBasedOmega[i]);
#endif
    double& val = ColBasedOmega[i];
    if (val < 0.0) {
      val = 0.0;
      ++zero_local;
    }
    if (val < min_local) min_local = val;
    if (val > max_local) max_local = val;
  }
#ifdef EMIN_IN_PAPER
fclose(fp);
#endif
  if (Denominator != NULL) ML_free(Denominator);
  if (Numerator   != NULL) ML_free(Numerator);

  double min_all  = ML_gmin_double(min_local, P0->comm);
  double max_all  = ML_gmax_double(max_local, P0->comm);
  double zero_all = ML_gsum_int(zero_local, P0->comm);

  if (ML_Get_PrintLevel() > 5 && P0->comm->ML_mypid == 0) { 
    cout << endl;
    cout << "Prolongator Smoothing: Using energy minimization (scheme = " 
         << ag->minimizing_energy << ")" << endl;
    cout << "Damping parameter: min = " << min_all <<  ", max = " << max_all 
         << " (" << zero_all << " zeros out of " << Nsubset_tot << ")" << endl;
    cout << "Dropping tolerance for DinvAP_0 = " << dropping << endl;
    cout << endl;
  }

  //  Stick the Omegas in their proper column (if they have been compressed).

  if (compress == 1) {
    double *TmpOmega = (double *) ML_allocate(sizeof(double)*(1+n_0));
    for (int i = 0; i < n_0; i++) {
       if (Subset[i] != -1) TmpOmega[i] =  ColBasedOmega[Subset[i]]; 
       else TmpOmega[i] = -666.;
    }
    for (int i = 0; i < n_0; i++) ColBasedOmega[i] = TmpOmega[i];  
    ML_free(TmpOmega);
  }

  // convert the omega's from column-based to row-based

  double* RowOmega = (double *) ML_allocate(sizeof(double)*(n + AmatNghost + 1));
  ag->old_RowOmegas = RowOmega;


  if (DinvAP0_subset->getrow->pre_comm != NULL) 
     ML_exchange_bdry(&ColBasedOmega[0],DinvAP0->getrow->pre_comm, n_0,
                      DinvAP0->comm, ML_OVERWRITE,NULL);

#ifdef EMIN_IN_PAPER
fp = fopen("rowomega","w");
#endif
  int AtLeastOneDefined, NRowOmegasSet = 0;
  for (int row = 0 ; row < n ; row++) {
     RowOmega[row] = -666.;            // RowOmega not set
     ML_get_matrix_row(DinvAP0, 1, &row, &allocated, &bindx, &val,
                      &row_length, 0);
     AtLeastOneDefined = 0;
     for  (int j = 0; j < row_length; j++) {
       double omega = ColBasedOmega[bindx[j]];
       if ( omega != -666. ) {  // ColBasedOmega not set
           AtLeastOneDefined = 1;
	   if (RowOmega[row] == -666.)  RowOmega[row] = omega;
           else if (omega < RowOmega[row]) RowOmega[row] = omega;
       }
    }
    if (AtLeastOneDefined == 1) {
       NRowOmegasSet++;
       if (RowOmega[row] < 0.) RowOmega[row] = 0.;
    }
#ifdef EMIN_IN_PAPER
fprintf(fp,"%20.13e\n",RowOmega[row]);
#endif
  }
#ifdef EMIN_IN_PAPER
fclose(fp);
#endif

  // If only a subset of column omegas are computed then we might not have
  // all of the RowOmega's defined. We need to loop a few times to get them.
  // Note: This should only be necessary if the compressed version is used.

  int Tries = 0;
  int Nrows_tot = ML_gsum_int(n, P0->comm);
  int NSet_tot  = ML_gsum_int(NRowOmegasSet, P0->comm);
  double *NewOmega = NULL;

  while ( NSet_tot < Nrows_tot) {
     if (NewOmega == NULL) 
        NewOmega = (double *) ML_allocate(sizeof(double)*(n + AmatNghost + 1));
     ML_exchange_bdry(RowOmega, Amat->getrow->pre_comm, Amat->invec_leng,
                       Amat->comm, ML_OVERWRITE, NULL);

     for (int row = 0 ; row < n ; row++) {
        if ( RowOmega[row] == -666.) {
           NewOmega[row] = -666.;
           int    NDefined = 0;
           double sum = 0.;
           ML_get_matrix_row(Amat, 1, &row, &allocated, &bindx, &val,
                             &row_length, 0);
           for  (int j = 0; j < row_length; j++) {
              double omega = RowOmega[bindx[j]];
              if (omega != -666.) {
                 sum = sum + omega; //just take average
                 NDefined++;
              }
           }
           if (NDefined != 0) {
              NewOmega[row] = sum/((double) NDefined);
              NRowOmegasSet++;
           }
        }
     } 
     for (int row = 0 ; row < n ; row++) {
        if ( (RowOmega[row] == -666.) && (NewOmega[row] != -666.))
           RowOmega[row] = NewOmega[row];
     }
     Tries++;
     if (Tries == 3) { // just sent any remaining ones to 0
        for (int row = 0 ; row < n ; row++) 
           if ( RowOmega[row] == -666.) RowOmega[row] = 0.;
        NRowOmegasSet = n;
     }
     NSet_tot  = ML_gsum_int(NRowOmegasSet, P0->comm);
  }
  if (NewOmega != NULL) ML_free(NewOmega);

  ML_Operator* OmegaDinvAP0 = 0;

  // Compute the new prolongator. Note: since the omegas are always scalar
  // quantities, we can always use the implicit scaling function here.

  OmegaDinvAP0 = ML_Operator_ImplicitlyVScale(DinvAP0, &RowOmega[0], 0);

  // The following sometimes crashes when used with MIS aggregation
  // (in parallel).
  ML_Operator_Add(P0, OmegaDinvAP0, &(ml->Pmat[clevel]), ML_CSR_MATRIX, -1.0);
  if (dropping != 0.0)
    ML_CSR_DropSmall(&(ml->Pmat[clevel]), 0.0, dropping, 0.0);

  ML_Operator_Set_1Levels(&(ml->Pmat[clevel]), &(ml->SingleLevel[clevel]),
                          &(ml->SingleLevel[level]));

#ifdef EMIN_IN_PAPER
// load in the matlab version of omega
system("/usr/local/matlab/bin/matlab -nodisplay < omegar.m >> outfile");
system("ls");
fp = fopen("newomegas","r");
for (int i = 0; i < Amat->invec_leng; i++) {
  fscanf(fp,"%lf",&(RowOmega[i]));
}
fclose(fp);
#endif

  if (bindx != NULL)    ML_free(bindx);
  if (val   != NULL)    ML_free(val);

  if (UnscaledAmat != NULL) {
     ML_Operator_Destroy(&Amat);
     Amat = UnscaledAmat; UnscaledAmat = NULL;
  }
  if (DinvAmat != NULL) ML_Operator_Destroy(&DinvAmat);
  if (DinvAP0_subset != DinvAP0) ML_Operator_Destroy(&DinvAP0_subset);
  if (DinvAP0)          ML_Operator_Destroy(&DinvAP0);
  if (DinvADinvAP0)     ML_Operator_Destroy(&DinvADinvAP0);
  if (OmegaDinvAP0)     ML_Operator_Destroy(&OmegaDinvAP0);
  if (Subset != NULL) ML_free(Subset);

#ifdef ML_TIMING
  ml->Pmat[clevel].build_time =  GetClock() - t0;
  ml->timing->total_build_time += ml->Pmat[clevel].build_time;
#endif
  return(0);
}

// ====================================================================== 
// generate smoothed restriction by minimizing energy                      
// Note: in the symmetric case we don't really need this, a transpose
// of P suffices.
//
// \author Marzio Sala and Ray Tuminaro, 9214
//
// \date 14-Jul-05
// ====================================================================== 

int ML_AGG_Gen_Restriction_MinEnergy(ML *ml,int level, int clevel, void *data)
{
  ML_Operator **prev_P_tentatives;
  ML_Aggregate * ag = (ML_Aggregate *) data;
  double       *tempOmega = NULL;

#ifdef ML_TIMING
  double t0;
  t0 =  GetClock();
#endif

  ML_Operator* Amat = &(ml->Amat[level]); //already created and filled
  prev_P_tentatives = ag->P_tentative; // for keep_P_tentative
  ML_Operator* P0 = prev_P_tentatives[clevel]; // already created and filled
  ML_Operator* P0_trans = NULL;
  ML_Operator* P0TA     = NULL;
  ML_Operator* ttt      = NULL;
  ML_Operator* Scaled_P0TA = NULL;

  double *RowOmega = ag->old_RowOmegas;
  bool NSR = false; // NonSmoothed Restriction

  if (NSR) 
  {
    // I just take the transpose of the tentative prolongator and go home.
    // The idea is nice but simply doesn't work at all -- you get the
    // same results of "pure" NSR (without energy minimization).
    if (ml->comm->ML_mypid == 0)
      printf("Using non-smoothed restriction\n\n");

    ML_Gen_Restrictor_TransP(ml, level, clevel, P0);
  }
  else
  {
    P0_trans = ML_Operator_Create(P0->comm);
    P0TA     = ML_Operator_Create(P0->comm);

    ML_Operator_Transpose_byrow(P0, P0_trans);

    if (Amat->num_PDEs == 1 || ag->block_scaled_SA == 0)
    {
      if (Amat->getrow->pre_comm != NULL) 
      {
        // communicate Dinv and RowOmega for ghost columns 
        ML_exchange_bdry(Dinv, Amat->getrow->pre_comm, Amat->invec_leng,
                         Amat->comm, ML_OVERWRITE, NULL);
        ML_exchange_bdry(RowOmega, Amat->getrow->pre_comm, Amat->invec_leng,
                         Amat->comm, ML_OVERWRITE, NULL);
        int Nghost = Amat->getrow->pre_comm->total_rcv_length;
        for (int i = 0; i < Amat->outvec_leng + Nghost; ++i)
          Dinv[i] *= RowOmega[i];
      }
      ML_Operator *Amat_Scaled;


      Amat_Scaled = ML_Operator_ImplicitlyVCScale(Amat, Dinv, 0);

      ML_2matmult(P0_trans, Amat_Scaled, P0TA, ML_CSR_MATRIX);
      ML_Operator_Destroy(&Amat_Scaled);

      ML_Operator_Add(P0_trans,P0TA,&(ml->Rmat[level]),ML_CSR_MATRIX,-1.0);
      if (SinglePrecision)
        ML_Operator_ChangeToSinglePrecision(&(ml->Rmat[level]));
    }
    else
    {
      // This code should work just fine, however, it is a bit of a waste
      // when num_PDEs = 1. 
      ML_2matmult(P0_trans, Amat, P0TA, ML_CSR_MATRIX);
      if (SinglePrecision)
        ML_Operator_ChangeToSinglePrecision(P0TA);

      ttt = ML_Operator_ImplicitlyBlockDinvScale(Amat);

      ML_AGG_DinvP(P0TA,(MLSthing *)ttt->data,Amat->num_PDEs,COL_SCALE_WITH_DT);

      int Nghost;
      Nghost = ML_CommInfoOP_Compute_TotalRcvLength(P0TA->getrow->pre_comm);

      tempOmega = (double *) ML_allocate(sizeof(double)*(P0TA->invec_leng + Nghost + 1));
      for (int i = 0; i < P0TA->invec_leng; i++) tempOmega[i] = RowOmega[i];

      ML_exchange_bdry(tempOmega,P0TA->getrow->pre_comm, P0TA->invec_leng, 
                       P0TA->comm, ML_OVERWRITE,NULL);
      Scaled_P0TA = ML_Operator_ImplicitlyVCScale(P0TA, &(tempOmega[0]), 0);

      ML_Operator_Add(P0_trans,Scaled_P0TA,&(ml->Rmat[level]),ML_CSR_MATRIX,-1.0);
    }
  }

  ML_Operator_Set_1Levels(&(ml->Rmat[level]), &(ml->SingleLevel[level]), 
                          &(ml->SingleLevel[clevel]));

  if ( ag->old_RowOmegas != NULL) ML_free(ag->old_RowOmegas);
  ag->old_RowOmegas = NULL;

  Dinv_size = -1;
  if (Dinv != NULL) {ML_free(Dinv); Dinv = NULL;}

#ifdef ML_TIMING
  ml->Rmat[level].build_time =  GetClock() - t0;
  ml->timing->total_build_time += ml->Rmat[level].build_time;
#endif

  if (ttt != NULL) ML_Operator_Destroy(&ttt);
  if (P0_trans != NULL) ML_Operator_Destroy(&P0_trans);  
  if (P0TA != NULL) ML_Operator_Destroy(&P0TA);
  if (Scaled_P0TA != NULL) ML_Operator_Destroy(&Scaled_P0TA);
  if (tempOmega != NULL) ML_free(tempOmega);

  return(0);
}

#include "ml_lapack.h"
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//		Begin Routines Written by JBS
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// ====================================================================== 
// The same as multiply_all, but use diagonal to scale one of the matrices 
// when you calculate the inner product between column(i) in left 
// and column(i) in right and store in InnerProd[i]. 
// 
// NOTE: InnerProd[] must have enough space to store ghost unknowns in right.
//
// Note: Must be CSR matrices but it is easy to change for a generic getrow()
//
void ML_multiply_all_vscale(ML_Operator* left, ML_Operator* right,
                                double* InnerProd, double* diagonal)
{
  int n = left->invec_leng;
  int n_rows = left->getrow->Nrows;
  int LeftGhost = 0, RightGhost = 0;
  int *LeftGlobal = NULL, *RightGlobal = NULL;
  int *LeftLocal  = NULL, *RightLocal  = NULL, *NewLeftLocal = NULL;

  
  if (n != right->invec_leng || n_rows != right->getrow->Nrows)
    {
      cerr << "Error: non-comparible operators" << endl;
      cerr << "(file " << __FILE__ << ", line " << __LINE__ << ")" << endl;
      exit(EXIT_FAILURE);
    }
  

  // compute the number of ghosts and make up some global ids

  if ( left->getrow->pre_comm != NULL) {
    ML_CommInfoOP_Compute_TotalRcvLength( left->getrow->pre_comm);
    LeftGhost =  left->getrow->pre_comm->total_rcv_length;
    ML_create_unique_id(left->invec_leng, &LeftGlobal, 
			left->getrow->pre_comm, left->comm, -1);
  }
  if (right->getrow->pre_comm != NULL) {
    ML_CommInfoOP_Compute_TotalRcvLength(right->getrow->pre_comm);
    RightGhost = right->getrow->pre_comm->total_rcv_length;
    ML_create_unique_id(right->invec_leng, &RightGlobal, 
			right->getrow->pre_comm, right->comm, -1);
  }
  if (  (LeftGhost > 0) && (RightGhost > 0) ) {

    LeftLocal    = (int *) ML_allocate(sizeof(int)*(n+LeftGhost));
    RightLocal   = (int *) ML_allocate(sizeof(int)*(n+RightGhost));
    NewLeftLocal = (int *) ML_allocate(sizeof(int)*(n+LeftGhost));

    // Sort global ids

    for (int i = 0; i <  n+LeftGhost; i++) LeftLocal[i] = i;
    for (int i = 0; i < n+RightGhost; i++) RightLocal[i] = i;

    ML_az_sort( LeftGlobal,  n+LeftGhost,  LeftLocal, NULL);
    ML_az_sort(RightGlobal, n+RightGhost, RightLocal, NULL);

    // Compute new local ids for the left matrix that correspond
    // to the local ids for the right matrix. If there is a left global
    // id that is not found in the right, assign it the id of n+RightGhost+1.
    // It will never be accessed and this avoids an if. Any index
    // that is not found in both matrices doesn't contribute to final result.

    for (int i = 0; i <  n+LeftGhost; i++) NewLeftLocal[i] = n+RightGhost+1;

    int j = 0;
    for (int i = 0; i <  n+LeftGhost; i++) {
      while ( (j < n+RightGhost) && ( RightGlobal[j] < LeftGlobal[i]) ) j++;
      if (RightGlobal[j] == LeftGlobal[i])
	NewLeftLocal[LeftLocal[i]] = RightLocal[j];
    }
  }

  for (int i = 0 ; i < n+RightGhost ; ++i) InnerProd[i] = 0.0;

  struct ML_CSR_MSRdata* left_data = 0;
  struct ML_CSR_MSRdata* right_data = 0;
  left_data  = (struct ML_CSR_MSRdata*)ML_Get_MyGetrowData(left);
  right_data = (struct ML_CSR_MSRdata*)ML_Get_MyGetrowData(right);

  int*    lrowptr = left_data->rowptr;
  int*    rrowptr = right_data->rowptr;
  double *temp_array = NULL;

  temp_array = (double *) ML_allocate(sizeof(double)*(n+RightGhost+1));
                                // use RightGhost because these indices are
                                // based off right matrix
  for (int i = 0; i < n+RightGhost+1; i++) temp_array[i] = 0.;

  double *rval = right_data->values;
  double *lval = left_data->values;
  int*    rbindx  = right_data->columns;
  int     rpos;
  int*    lbindx1  = left_data->columns;
  int*    lbindx2  = left_data->columns;
  for (int row = 0 ; row < n_rows ; ++row)
    {
      int     llen    = lrowptr[row + 1] - lrowptr[row];
      int     rlen    = rrowptr[row + 1] - rrowptr[row];

      if (NewLeftLocal == NULL) 
        for (int i = 0 ; i < llen ; ++i) temp_array[*lbindx1++] = *lval++;
      else
        for (int i = 0 ; i < llen ; ++i)
	  temp_array[NewLeftLocal[*lbindx1++]] = *lval++;

      for (int i = 0 ; i < rlen ; ++i) {
	rpos = *rbindx++;
	//put in the diagonal scaling of the right matrix
	InnerProd[rpos] += diagonal[row]*(*rval++) * temp_array[rpos];
      }

      if (NewLeftLocal == NULL) 
	for (int i = 0; i < llen; ++i) temp_array[*lbindx2++] = 0.;
      else
	for (int i = 0; i < llen; ++i) temp_array[NewLeftLocal[*lbindx2++]] = 0.;
    }
  if (right->getrow->pre_comm != NULL)
    ML_transposed_exchange_bdry(InnerProd, right->getrow->pre_comm, n, 
				right->comm, ML_ADD);

  if (temp_array  != NULL) ML_free(temp_array);
  if (NewLeftLocal!= NULL) ML_free(NewLeftLocal);
  if (  LeftLocal != NULL) ML_free(LeftLocal);
  if ( RightLocal != NULL) ML_free(RightLocal);
  if (RightGlobal != NULL) ML_free(RightGlobal);
  if ( LeftGlobal != NULL) ML_free(LeftGlobal);
}

// ====================================================================== 
// Begin Four routines to wrap an arbitrary CSR Mat so that
// All of its entries have an absolute value applied to them
// Written by JBS -- they are routines for ImplicitScale with minor modifications
void ML_ImplicitAbs_Destroy(void *data)
{
  struct ml_matscale *temp;

  temp = (struct ml_matscale *) data;
  if (temp != NULL) {
    if (temp->destroy_child) ML_Operator_Destroy( &(temp->Amat));
    ML_free(temp);
  }
}

/* ******************************************************************** */
/* Getrow function that is uses the wrapped mat's getrow to get the     */
/* correct sparsity pattern, and then sets everything to abs(.) of 	*/
/* stored value 							*/
/* ******************************************************************** */
int ML_ImplicitAbs_Getrow(ML_Operator *data, int N_requested_rows, 
			  int requested_rows[], int allocated_space, 
			  int columns[], double values[], 
			  int row_lengths[])
{
  struct ml_matscale *temp;
  int    i, status = 1, size = 0;
 
  if (N_requested_rows > 1) {
    printf("ML_implicitmatscale_getrow: Not implemented for > 1 row at a time\n");
    exit(1);
  }
  temp = (struct ml_matscale *) ML_Get_MyGetrowData(data);
  status = ML_Operator_Getrow(temp->Amat, N_requested_rows, requested_rows,
			      allocated_space, columns,
			      values, &size );
  if (status) {
    for (i = 0; i < size; i++) values[i] = fabs(values[i]);// old value 1.0;
    row_lengths[0] = size;
  }
  return(status);
}

//Implict absolute value MatVec
int ML_ImplicitAbs_Matvec(ML_Operator *Amat_in, int ilen, double p[], 
			  int olen, double ap[])
{
  struct ml_matscale *temp;
  int    status = 1;

  temp = (struct ml_matscale *) ML_Get_MyGetrowData(Amat_in);
  status = ML_Operator_Apply(temp->Amat, ilen, p, olen, ap);

  return(status);
}


//The ImplicitAbs wrapper
ML_Operator *ML_Operator_ImplicitAbs(ML_Operator *Amat, int OnDestroy_FreeChild)
{
  ML_Operator *matrix;
  struct ml_matscale *new_data;

  matrix = ML_Operator_Create(Amat->comm);

  new_data = (struct ml_matscale *) ML_allocate( sizeof(struct ml_matscale));
  if (new_data == NULL)
    pr_error("ML_Operator_ImplicitAbs: out of space\n");
  new_data->Amat          = Amat;
  new_data->destroy_child = 0;
  ML_Operator_Set_ApplyFuncData(matrix,Amat->invec_leng, 
				Amat->outvec_leng,new_data,
				Amat->matvec->Nrows, ML_ImplicitAbs_Matvec,
				Amat->from_an_ml_operator);

  ML_Operator_Set_Getrow(matrix,Amat->getrow->Nrows,ML_ImplicitAbs_Getrow);
  matrix->data_destroy   = ML_ImplicitAbs_Destroy;
  if (OnDestroy_FreeChild) new_data->destroy_child = 1;


  /* Note: this is need for any functions doing getrow(). The matvec  */
  /* wrapper does not invoke communication as it is already contained */
  /* in the lower level (unscaled) matvec function.                   */

  if (Amat->getrow->pre_comm != NULL) 
    ML_CommInfoOP_Clone( &(matrix->getrow->pre_comm),Amat->getrow->pre_comm); 

  return matrix;
}

//		End Four Routines to Define ABS(.) Wrapper
//=====================================================================


// =====================================================================
// Written by JBS
// Input: CSR_Data := CSR Data for a matrix
// 	  nRows := number of rows
// Output:  CSR_Data := defines identical matrix that is inputted, only the entries for each row 
//			Are sorted according to column
void ML_Sort_Cols(struct ML_CSR_MSRdata * CSR_Data, int nRows)
{
  int i, row_start, row_end;
	
  for(i = 0; i < nRows; i++)
    {
      row_start = CSR_Data->rowptr[i];
      row_end = CSR_Data->rowptr[i+1];
      ML_az_sort( &(CSR_Data->columns[row_start]),  row_end-row_start,  NULL, &(CSR_Data->values[row_start]) );
    }

  return;
}


// =====================================================================
// Written By JBS
// Input: A, in CSR format
// Output:  The maximum magnitude entry in an ML CSR mat

double ML_MaxEntry(ML_Operator * A)
{
  int i, nnzs;
  double max = 0.0;
  double value;
  struct ML_CSR_MSRdata * CSR_Data = (struct ML_CSR_MSRdata *) A->data;
	
  nnzs = CSR_Data->rowptr[A->outvec_leng];
	
  for(i = 0; i < nnzs; i++)
    {
      value = fabs(CSR_Data->values[i]);
      if(value > max)
	{	max = value; }
    }
	
  return max;
}

// =====================================================================
// Written by JBS
// Input:  A := matrix, must be in CSR format, with each row sorted according to column
//	   Pattern := A desired sparsity pattern in CSR format
// Output: A := original matrix with 0's now replacing values in A that did appear as 
//		nonzeros in Pattern
void ML_Enforce_Sparsity(ML_Operator * A, struct ML_CSR_MSRdata *Pattern)
{
  int Arow_start, Arow_end, Prow_start, Prow_end, Acurrent_col, Pcurrent_col, k, zero_flag;
  struct ML_CSR_MSRdata * Acsr_data = ( struct ML_CSR_MSRdata *) A->data;
	
  //traverse the rows of A and P, which are assumed to be the same
  for(int i = 0; i < A->outvec_leng; i++)
    {
      Arow_start = Acsr_data->rowptr[i]; 
      Arow_end =   Acsr_data->rowptr[i+1];
      Prow_start = Pattern->rowptr[i]; 
      Prow_end = Pattern->rowptr[i+1];
      k = Prow_start;
		
      //		if(i == 63) { cout << "Arow_start = " << Arow_start << ".  Arow_end = " << Arow_end << endl; cout << "Prow_start = " << Prow_start << ".  Prow_end = " << Prow_end << endl; } 

		
      for(int j = Arow_start; j < Arow_end; j++)
	{
	  //traverse Pattern until you reach or go past this entry
	  Acurrent_col = Acsr_data->columns[j];
			
	  //			if(i == 63) { cout << "j = " << j << ".  Acurrent_col = " << Acurrent_col << endl; }
			
	  zero_flag = 0;
	  if(0)	//dumb loop.  assumes colun entries aren't sorted
	    {
				//look for this element in pattern
	      for(int m = Prow_start; m < Prow_end; m++)
		{
		  Pcurrent_col = Pattern->columns[m];
		  if(Pcurrent_col == Acurrent_col)
		    {	zero_flag = 1; break;}
		}
				//If we've never found this element, zero it out
	      if(zero_flag == 0)
		{	Acsr_data->values[j] = 0.0; }

	    }
	  else	//This assumes that the column entries are always sorted.
	    {
	      for(int m = k; m < Prow_end; m++)
		{
		  Pcurrent_col = Pattern->columns[m];
					
		  //					if(i == 63) { cout << "m = " << m << ".  Pcurrent_col = " << Pcurrent_col << endl; }
				
		  //If this entry exists in Pattern, keep it
		  if(Pcurrent_col == Acurrent_col)
		    {	
		      //if we've reached the end of the row in Pattern, zero the rest of row i in A
		      if(m == (Prow_end-1) )
			{	zero_flag = 1; 
			//							if(i == 63) { cout << "Here 1" << endl;}
			break; }
		      else
			{	k = m+1; 
			//							if(i == 63) { cout << "Here 2" << endl;}
			break;	}
		    }
		  //If it doesn't, zero it out
		  else if(Pcurrent_col > Acurrent_col)
		    {	
		      Acsr_data->values[j] = 0.0; 
		      k = m; 
		      //						if(i == 63) { cout << "Here 3" << endl; }
		      break;	
		    }
		  else if( Acurrent_col > Pcurrent_col)
		    {
	
		      //						if(i == 63) { cout << "Here 4" << endl; }
		      if(m == (Prow_end-1) )
			{
			  //							if(i == 63) { cout << "Here 5" << endl; }
			  zero_flag = 1; 
			  j--;
			  break; 
			}
		    }
		}

				//Handle Special Case that this row is zero in Pattern
	      if(Prow_start == Prow_end)
		{	j = Arow_start - 1; zero_flag = 1;}
			
				//if we've reached the end of the row in Pattern, zero the rest of row i in A
	      if(zero_flag == 1)
		{
		  //					if(i == 63) {cout << "Here 6" << endl;}
		  for(j = j+1; j < Arow_end; j++)
		    {	
		      Acsr_data->values[j] = 0.0; 
		      //						if(i == 63) { cout << "j = " << j << endl;}
		    }
		  break;
		}
	    }
	}
		
    }
	
  return;
}

// ======================================================================
//Written by JBS
// Input: mat := a matrix stored in a linear vector in column major format
//	  rows
//	  cols
//	  FileName := File to delete and then write mat to
// Output: ./FileName  := a file with mat written in  "i j  value" format
//			  This is a file readable by Matlab's spconvert function
void ML_print_mat(double * mat, int rows, int cols, char FileName[])
{
  int counter = 0;
  char str[80];
	
  //open file for writing, overwriting any previous contents
  FILE * fid = fopen(FileName, "w");
	
  for(int j = 0; j < cols; j++)
    {
      for(int i = 0; i < rows; i++)
	{
	  sprintf(str, "%d   %d      %1.16e \n", i+1, j+1, mat[counter]);
	  fprintf(fid, str);  
	  counter++;
	}
    }		

  return;
}

// ====================================================================
//Written by JBS
// Input: A, in CSR format
// Output: A := with only nonzero entries.  All stored 0's in A are removed.  
//		Requires allocating space for another copy of A while entries 
//		are being removed

void ML_Squeeze_Out_Zeros(ML_Operator *A)
{
  struct ML_CSR_MSRdata * ACSRdata = (struct ML_CSR_MSRdata *) A->data;
  int counter, i, j, rowstart, rowend;
  double * values_new;
  int *columns_new, *rowptr_new;
  int numRows = A->outvec_leng;
  int NNZ_A = ACSRdata->rowptr[numRows];
	
  //We first must count the number of non-zeros
  counter = 0;
  for(i = 0; i < NNZ_A; i++)
    {	
      if(fabs(ACSRdata->values[i]) > 1e-10)
	{	counter++; }
    }

  //Allocate new arrays
  columns_new = (int *) ML_allocate( sizeof(int)*counter);
  rowptr_new = (int *) ML_allocate( sizeof(int)*(numRows + 1) );
  values_new = (double *) ML_allocate( sizeof(double)*counter );
	
  //Drop zero entries in A.  
  counter = 0;
  rowptr_new[0] = 0;
  for(i = 0; i < numRows; i++)
    {	
      rowstart = ACSRdata->rowptr[i];
      rowend = ACSRdata->rowptr[i+1];
		
      for(j = rowstart; j < rowend; j++)
	{
	  //if the value is ``non-zero"
	  if(fabs(ACSRdata->values[j]) > 1e-10)
	    {	
	      columns_new[counter] = ACSRdata->columns[j];
	      values_new[counter] = ACSRdata->values[j];
	      counter++; 
	    }
	}

      //mark where this row ends
      rowptr_new[i+1] = counter;
    }

  //free old CSR data for A so that we increase memory usage as little as possible
  if(ACSRdata->rowptr != NULL) ML_free(ACSRdata->rowptr); 
  if(ACSRdata->columns != NULL) ML_free(ACSRdata->columns); 
  if(ACSRdata->values != NULL) ML_free(ACSRdata->values); 

  //set ACSRdata to point to new CSR information
  ACSRdata->columns = columns_new;
  ACSRdata->rowptr = rowptr_new;
  ACSRdata->values = values_new;
	
  return;
}


// ====================================================================
// Written By JBS
// Input: A
//	  droptol := row-based drop tolerance
// Output: A  := Is identical to A passed in, only that entries in row, i,
//		that are less than droptol*max_entry(A(i,:)) have been set
//		to zero

void ML_Drop(ML_Operator * A, double droptol)
{
  int i, j, rowstart, rowend, numrows;
  double max, val;
  struct ML_CSR_MSRdata * ACSRdata = (struct ML_CSR_MSRdata *) A->data; 
	
  numrows = A->outvec_leng;
	
  for(i = 0; i < numrows; i++)
    {
      rowstart = ACSRdata->rowptr[i];
      rowend = ACSRdata->rowptr[i+1];

      //find maximumn POSITIVE!! off-diagonal for this row
      max = 0.0;
      for(j = rowstart; j < rowend; j++)
	{
	  //If not diagonal
	  if(ACSRdata->columns[j] != i)
	    {
				//NO ABS VALUE, i think.  JBS
	      val = ACSRdata->values[j];
	      if( val > max)
		{	max = val;}
	    }
	}
		
      //set tolerance for this row;
      max = droptol*max;

      //Apply drop-tolerance to this row
      for(j = rowstart; j < rowend; j++)
	{
	  //If not diagonal
	  if(ACSRdata->columns[j] != i)
	    {
				//do we want to take an absolute value or not?  i.e.
				//	do we want to just drop negative values and assume that 
				//	negative's imply a weak connection?
				//	Currently, we say no.
	      if(  fabs(ACSRdata->values[j]) < max)
		{	ACSRdata->values[j] = 0.0; }
	    }
	}
    }
	
  return;
}



// ====================================================================
// Written By JBS
// Input:  A
//	   num_steps  := number of times to evolve a delta function in order
//			calculate strenth information
//	   t_final    := t_final/rho(DinvA) is the stop time for the time 
//			 evolution of a delta function
//	   drop_tol   := drop tolerance to apply to the ODE based strength 
//		         matrix.  Used by the routin "Drop"
// Output: ML_Operator *   := A strength of connection matrix that is generated 
// 				by evolving delta functions in time.  Essentially, this
//				matrix is (I - deltaT*Dinv*A)^num_steps.  
//				Energy based post-processing can be turned to 
//				improve strength information
ML_Operator * ML_ODE_Strength_Matrix(ML_Operator * A, int num_steps, double t_final, double drop_tol)
{
  int i;
  ML_Operator * Atilde, *I, *S, *Sodd, *Seven, *omegaDinvA, *DinvA, *AtildeNew;
  double * diago;			
  double * diago_local =  (double *) ML_allocate(sizeof(double)*(A->outvec_leng));
	
  Atilde = ML_Operator_Create(A->comm); 
  I = ML_Operator_Create(A->comm); 
  S = ML_Operator_Create(A->comm); 
  Sodd = ML_Operator_Create(A->comm); 
  Seven = ML_Operator_Create(A->comm); 
  omegaDinvA = ML_Operator_Create(A->comm); 
  AtildeNew = ML_Operator_Create(A->comm);
  DinvA = ML_Operator_Create(A->comm); 

  //Create an Identity
  //In the words of Ray, this is a Kludge, but I borrowed this from 
  //	ml_amg_genP
  ML_Operator_Set_ApplyFuncData(I, A->invec_leng,
				A->outvec_leng, (void*) A,
				A->matvec->Nrows, NULL, 0);
  ML_Operator_Set_Getrow(I, A->getrow->Nrows, 
			 ML_AMG_Identity_Getrows);
  ML_CommInfoOP_Clone(&(I->getrow->pre_comm), A->getrow->pre_comm);


  //Create S
  //Get the inverse diagonal
  ML_Operator_Getrow_Diag(A, &diago);
  for (int row = 0; row < A->outvec_leng; row++) 
    { diago_local[row] = 1.0/diago[row];	}

  //Implicitly Calculate Dinv*A, i.e. wrap around A with a Dinv routine
  DinvA = ML_Operator_ImplicitlyVScale(A, diago_local, 0);    

  //Calculate the spectral radius of Dinv*A for the SA weight
  ML_Gimmie_Eigenvalues(DinvA, ML_NO_SCALE, ML_USE_POWER, ML_NO_SYMMETRIZE);					

  //Implicitly wrap around Dinv*A with another diagonal scaling.
  //	Scale by the time step size where t_final/rho is the real t_final and we divide
  //	by the number of desired time steps
printf("these guys are %e  %d   %e\n",t_final,num_steps,DinvA->lambda_max); 
  omegaDinvA = ML_Operator_ImplicitlyScale(DinvA, t_final/(num_steps*(DinvA->lambda_max)), 0); 
	
  // Add I to a constant times OmegaDinvA and store the result in S
  ML_Operator_Add(I, omegaDinvA, S, ML_CSR_MATRIX, -1.0 );             
	
  //Set up Sodd for the below iterations to equal S
  ML_Operator_Add(S, I, Sodd, ML_CSR_MATRIX, 0.0 ); 
	
  //free up some memory
  ML_Operator_Destroy(&omegaDinvA);
  ML_Operator_Destroy(&DinvA);
  ML_free(diago_local);

  //Calculate the raw strength matrix, note that we have already calculated the first iteration in S
  //	Could easily change to calculate S^{2*i} so that S
  //	raised to powers of two could be efficiently calculated
  for(i = 1; i < num_steps; i++)
    {
      if( (i%2) != 0 )
	{
	  ML_Operator_Destroy(&Seven);
	  Seven = ML_Operator_Create(A->comm);
	  //Seven = Sodd*S
	  ML_2matmult(Sodd, S, Seven, ML_CSR_MATRIX);	        
	}
      else
	{
	  ML_Operator_Destroy(&Sodd);
	  Sodd = ML_Operator_Create(A->comm);
	  //Sodd = Seven*S
	  ML_2matmult(Seven, S, Sodd, ML_CSR_MATRIX);
	}
    }

		
  //Free up space and Transpose.
  //If Amat is symmetric, S^{i} need not be transposed.  But in general
  //	the strength information is in columns, so this is necessary.
  if( (num_steps%2) != 0)
    {	
      ML_Operator_Destroy(&S);  ML_Operator_Destroy(&Seven);
      ML_Operator_Transpose_byrow(Sodd, Atilde);
      ML_Operator_Destroy(&Sodd);
    }
  else
    {	
      ML_Operator_Destroy(&S);  ML_Operator_Destroy(&Sodd);
      ML_Operator_Transpose_byrow(Seven, Atilde);
      ML_Operator_Destroy(&Seven);
    }
	
  //Choose which type of post-processing to use.
  if(1)	//Apply a simple drop tolerance
    {
      //Enforce sparsity of Amat, this is Kludge-Like, because A and Atilde must
      //	be in CSR format for routines to work
      S = ML_Operator_Create(A->comm);
      ML_Operator_Add(A, I, S, ML_CSR_MATRIX, 0.0);
      ML_Operator_Add(Atilde, I, Atilde, ML_CSR_MATRIX, 0.0 ); 
      ML_Sort_Cols((struct ML_CSR_MSRdata *) Atilde->data, Atilde->outvec_leng);
      ML_Sort_Cols((struct ML_CSR_MSRdata *) S->data, S->outvec_leng);
      ML_Enforce_Sparsity(Atilde, ((struct ML_CSR_MSRdata *) S->data) );
		
      //free up space
      ML_Operator_Destroy(&S);
      ML_Operator_Destroy(&I);
		

      //In case diagnostic output is desired...
      if(0)
	{
	  //print out a row of Atilde
	  int Dimen = (int) floor(sqrt((double) A->outvec_leng));
	  int x = Dimen/2;
	  int y = Dimen/2;
	  int point = (x - 1)*Dimen + y;

	  //print out row "point" from Atilde
	  struct ML_CSR_MSRdata * AtildeCSRdata = (struct ML_CSR_MSRdata *) Atilde->data;
	  cout << endl << "Dimension:  " << Dimen << endl;
	  cout << "Atilde at point (" << x << "," << y << ") is ..." << endl;
	  for(int j = (AtildeCSRdata->rowptr[point]); j < (AtildeCSRdata->rowptr[point + 1]); j++)
	    {  	cout << "(" << point << "," << AtildeCSRdata->columns[j] << ")  = " << (AtildeCSRdata->values[j]) << endl; }
	  cout << endl;
	}
		
      //Apply Drop Tolerance
      ML_Drop(Atilde, drop_tol);

      //It is important to get rid of 0 entries, as we wrap this non-zero pattern with 1's later.
      ML_Squeeze_Out_Zeros(Atilde);
		
      return Atilde;
    }
  else
    {
      //Do the energy-based post-processing of Scott's measure
      //As it stands, it uses 3 vectors of numRows length and transposes A -- Not exactly cheap.
      //Atilde should already be transposed, if necessary.
      //Transpose A now, is not necessary if A is symmetric...this is an expensive operation	
      ML_Operator *Atrans = ML_Operator_Create(A->comm); 
      ML_Operator_Transpose_byrow(A, Atrans);
	
      //Create a new Atilde with the sparsity pattern preset to A's
      //	This is Kludge-like...but works and should be too inefficient
      ML_Operator_Add(A, I, AtildeNew, ML_CSR_MATRIX, 0.0);
	
      //Declare some variables for just this algorithm
      int rowstart, rowend, rowstartA, rowendA, j, k, rowstart2, rowend2;
      int numRows = A->outvec_leng;
      double denom, val;
      double * g = (double *) ML_allocate(sizeof(double)*(numRows));
      double * Ag = (double *) ML_allocate(sizeof(double)*(numRows));
      double * AgMod = (double *) ML_allocate(sizeof(double)*(numRows));
      struct ML_CSR_MSRdata * AtildeCSRdata = (struct ML_CSR_MSRdata *) Atilde->data;
      struct ML_CSR_MSRdata * AtildeNewCSRdata = (struct ML_CSR_MSRdata *) AtildeNew->data;
		
      //We need A's sparsity pattern, but it may not be a CSR matrix...get CSR sparsity pattern in Kludge-like fashion
      S = ML_Operator_Create(A->comm);
      ML_Operator_Add(A, I, S, ML_CSR_MATRIX, 0.0);
      struct ML_CSR_MSRdata * ACSRdata = (struct ML_CSR_MSRdata *) S->data;
		
      //zero out g
      for(i = 0; i < numRows; i++)
	{	g[i] = 0.0; }
		
		
      //Begin loop to calculate AtildeNew(i,j) = [ || Atilde(i,:)_j ||_A / || Atilde(i,:) ||_A ]  - 1
      //	where Atilde(i,:)_j is that row of Atilde with entry j zeroed out
      for(i = 0; i < numRows; i++)
	{
	  rowstart = AtildeCSRdata->rowptr[i];
	  rowend = AtildeCSRdata->rowptr[i+1];
	  rowstartA = ACSRdata->rowptr[i];
	  rowendA = ACSRdata->rowptr[i+1];
			
	  //set g = Atilde(i,:)
	  for(j = rowstart; j < rowend; j++)
	    {	g[ (AtildeCSRdata->columns[j]) ] = AtildeCSRdata->values[j]; }
				
	  //set Ag = A*g
	  ML_Operator_Apply(A, A->invec_leng, &(g[0]), numRows, &(Ag[0]));
			
	  // ||g||_A
	  denom = sqrt(ML_gdot(numRows, g, Ag, A->comm));
 
	  //we calculate strength of connection only at the nonzero locations of A
	  //   each iteration of this loop calculates AtildeNew(i,A->columns[j])
	  for(j = rowstartA; j < rowendA; j++)
	    {
				//We now mimic zeroing out an element of g and multiplying by A, by simply subtracting off 
				//	( a row of Atrans * the corresponding entry of g) from Ag
				//	AgMod = Ag - g[cols[j]]*Atrans(cols[j],:)
	      val = g[ (ACSRdata->columns[j]) ];
	      memcpy( &(AgMod[0]), &(Ag[0]), sizeof(double)*numRows ); 
	      rowstart2 = ACSRdata->rowptr[ (ACSRdata->columns[j]) ];
	      rowend2 = ACSRdata->rowptr[ (ACSRdata->columns[j] + 1) ];
	      for(k = rowstart2; k < rowend2; k++)
		{	AgMod[ACSRdata->columns[k]] -= val*ACSRdata->values[k];	}
				
				//AgMod' * g
	      g[ACSRdata->columns[j]] = 0.0;
				
				//Set value AtildeNew(i, A->columns[j])
	      AtildeNewCSRdata->values[j] = ( sqrt(ML_gdot(numRows, AgMod, g, A->comm))/denom ) - 1;
				
				//set g back
	      g[ACSRdata->columns[j]] = val;
	    }
			
	  //Finally, zero out all possible non-zeros in g, as opposed to zeroing out the whole vector
	  for(j = rowstart; j < rowend; j++)
	    {	g[AtildeCSRdata->columns[j]] = 0.0; }
	}

		
      //if curious, can print out Atilde
      //char str[80];
      //sprintf(str,"temp");
      //ML_Operator_Print(Atilde,str);   

      //Clean Up
      ML_Operator_Destroy(&I);	
      ML_Operator_Destroy(&Atilde);
      ML_Operator_Destroy(&Atrans);
      ML_Operator_Destroy(&S);
		
      //Diagonal of Atilde is probably negative.
      //Set diagonal equal to maximum entry.  Don't worry...the drop 
      //tolerance function below just drops with respect to the largest
      //off diagonal entry in a row
      double Biggest = ML_MaxEntry(AtildeNew);
      for(i = 0; i < numRows; i++)
	{
	  rowstart = AtildeNewCSRdata->rowptr[i];
	  rowend = AtildeNewCSRdata->rowptr[i+1];
	
	  //loop over this row until you get to diagonal entry...
	  //	set it to max and break 
	  for(j = rowstart; j < rowend; j++)
	    {
	      if(AtildeNewCSRdata->columns[j] == i)
		{	AtildeNewCSRdata->values[j] = Biggest; }
	    }
	}
		
		
      //In case diagnostic output is desired...
      if(1)
	{
	  //print out a row of AtildeNew
	  int Dimen = (int) floor(sqrt((double) A->outvec_leng));
	  int x = Dimen/2;
	  int y = Dimen/2;
	  int point = (x - 1)*Dimen + y;

	  //print out row "point" from AtildeNew
	  cout << endl << "Dimension:  " << Dimen << endl;
	  cout << "Atilde at point (" << x << "," << y << ") is ..." << endl;
	  for(int j = (AtildeNewCSRdata->rowptr[point]); j < (AtildeNewCSRdata->rowptr[point + 1]); j++)
	    {  	cout << "(" << point << "," << AtildeNewCSRdata->columns[j] << ")  = " << (AtildeNewCSRdata->values[j]) << endl; }
	  cout << endl;
	}
		
      ML_free(g);  ML_free(Ag); ML_free(AgMod);
		
      //Apply Drop Tolerance and Return
      ML_Drop(AtildeNew, drop_tol);
      return AtildeNew;
    }
}


// =====================================================================
// Written by JBS
// This routine satisfies the constraints for the energy minimization algorithm
// for generating P.  The constraints are enforcing preservation of nullspace 
// at desired locations vectors and preservation of a sparsity pattern.
// ---- Input ----
// Update:  is the size of P0 and is added upon return to P in order to step in the 
//          direction of minimizing the energy of P inside the space of the constraints.
// Pattern: sparsity pattern that must be enforced on P and Update
// Bone:    coarse grid nullspace vectors
// BtBinv:  is a local Cholesky factor for each node (not for each dof) of B_i^T*B_i, where 
//	    B_i is B restricted to node i's allowed nonzeros in Pattern
// F:       holds a NullDim sized vector for each Node that is binary.  It switches on and off the 
//	    enforcement of constraints.  0 => off     1 => on
// numPDEs, numDOFs, numNodes are self-explanatory
//
// ---- Output----  
// Update is changed so that Update*Bone = 0 with the changes occuring only at allowed non-zeros

void ML_Satisfy_Constraints(ML_Operator *Update, ML_Operator *Pattern, double *Bone, double *BtBinv, int * F, int numPDEs, int numDOFs, int numNodes, int NullDim)
{
  int i, j, k, length, rowstart, rowend, counter, indx, nodeNullDim, node;
  int oneInt = 1;
  int numCoarseDOFs = Update->invec_leng;
  int NullDimSq = NullDim*NullDim;
  char NN = 'N';
  double *NQ, *NQlocal, *BNQlocal, *Bonelocal= NULL, *Updatelocal;
  double one = 1.0;
  double NegOne = -1.0;
  double zero = 0.0;
  struct ML_CSR_MSRdata * PatternCSRdata = (struct ML_CSR_MSRdata *) Pattern->data;
	
  NQlocal = (double *) ML_allocate( sizeof(double)*NullDim );
  BNQlocal = (double *) ML_allocate( sizeof(double)*NullDim );
	
  //Construct NQ = Update*Bone, where NQ is stored in column major
  //	order when Bone has multiple columns, i.e. multiple nullspace vects
  NQ =  (double *) ML_allocate( sizeof(double)*numDOFs*NullDim );
  for(i = 0; i < NullDim; i++)
    {
      //Update*Bone(:,i) ===> NQ(:,i)
      ML_Operator_Apply(Update, Update->invec_leng, &(Bone[i*numCoarseDOFs]), Update->outvec_leng, &(NQ[i*numDOFs]));
    }
	

  //Begin enforce constraints loop
  for(i = 0; i < numDOFs; i++)
    {	
      //Assume integer division truncates.  Calculate (current node number).
      node = ( (int) (i/numPDEs));
      nodeNullDim = node*NullDim;
		
      //Construct NQlocal
      for(j = 0; j < NullDim; j++)
	{	
	  //Use F to enforce constraints
	  if(F[nodeNullDim + j] == 1)
	    {	NQlocal[j] = NQ[j*numDOFs + i]; }
	  else
	    {	NQlocal[j] = 0.0; }
	}
	
      //if(i == 3) {for(j = 0; j < NullDim; j++) cout << NQlocal[j] << endl; }
	
      //We use the sparsity pattern for this node to decide which entries of the NullSpace vectors to use
      //	We assume a unifrom sparsity pattern for each dof on a node
      indx = i;// - i%numPDEs;
      rowstart = PatternCSRdata->rowptr[indx];
      rowend = PatternCSRdata->rowptr[indx + 1];
      length = rowend - rowstart;
		
      //Construct Bonelocal in scratch memory, only do for every node, not every dof -- 
      // commented if stmt out for now, to agree with Matlab
      // if( (i%numPDEs) == 0 )
      {
	//Clean up previous Bonelocal
	if( i > 0)
	  { if(Bonelocal != NULL) ML_free(Bonelocal); }
			
	Bonelocal = (double *) ML_allocate( sizeof(double)*length*NullDim );
	counter = 0;
	for(k = 0; k < NullDim; k++)  //loop over cols
	  {
				//Use F to enforce constraints
	    if(F[nodeNullDim + k] == 1)
	      {
		for(j = rowstart; j < rowend; j++)  //loop over rows
		  {
		    indx = PatternCSRdata->columns[j] + (k*numCoarseDOFs);
		    Bonelocal[counter] = Bone[indx];
		    counter++;
		  }
	      }
	    else // don't enforce constraint
	      {
		for(j = 0; j < length; j++)	//loop over rows
		  {	Bonelocal[counter] = 0.0; counter++; }
	      }
	  }
      }


      //DGEMM:  &(BtBinv[node*NullDimSq]) * NQlocal^T  ===> BNQlocal^T
      //BtBinv is NullDim X NullDim, NQlocal is 1 X NullDim
      // params: X,  Y,   op(X)rows, op(Y)columns, op(X)columns,   ALPHA,               X,            Xrows,        Y,         Yrows,   BETA,     Z,Zrows  
      DGEMM_F77(&NN, &NN, &NullDim,    &oneInt,     &NullDim,      &one,   &(BtBinv[node*NullDimSq]), &NullDim, &(NQlocal[0]), &NullDim, &zero,  &(BNQlocal[0]), &NullDim);
		
      //DGEMM of (-1*BNQlocal^T * Bonelocal^T) ==> Updatelocal
      //       i.e. (-1*Bonelocal*BNQlocal)^T ==> Updatelocal
      Updatelocal = (double *) ML_allocate( sizeof(double)*length );
      //params:  X'   Y', op(X)rows, op(Y)columns, op(X)columns,  ALPHA,       X,           Xrows,      Y,          Yrows,   BETA,         Z,         Zrow
      DGEMM_F77(&NN, &NN, &length,      &oneInt,     &NullDim,   &NegOne, &(Bonelocal[0]), &length, &(BNQlocal[0]), &NullDim, &zero, &(Updatelocal[0]), &length);

		
      //if(i == 3) {for(j = 0; j < length; j++) cout << Updatelocal[j] << endl; }
		
      //Here we do something tricky.  We insert the new values in Updatelocal into Pattern
      //Pattern is already an MLOperator, so it doesn't use any more space.  The values 
      //of Pattern are also not used anywhere, only its rowptr and columns entries are used.
      //We will then later do an ML_Operator_Add of Pattern and Update.
      //Add Updatelocal to Update(i,:)
      counter = 0;
      rowstart = PatternCSRdata->rowptr[i];
      rowend = PatternCSRdata->rowptr[i+1];
      for(j = rowstart; j < rowend; j++)
	{
	  PatternCSRdata->values[j] = Updatelocal[counter];
	  counter++;
	}
		
      // -----Clean Up this loop iteration-----
      if(Updatelocal != NULL) ML_free(Updatelocal);
		
    }	

  //We now update Update with an MLOperator Add
  //              1.0*A +   B    =   C
  ML_Operator_Add(Update, Pattern, Update, ML_CSR_MATRIX, 1.0   );

  ML_free(NQlocal);
  ML_free(BNQlocal);
  ML_free(NQ); 
}

// ====================================================================== 
// Written by JBS
// Generate smooth prolongator by minimizing energy according to Mandel's
// Paper, "Energy Optimization of Algebraic Multigrid Bases."  Uses the 
// above helper functions.
// Input: ml  := ml hierarchy
//	  level  := current fine grid level
//       clevel  := current coarse grid level, integer +/- 1 from level
//	 data    := holds aggregate information
// Output  P	 := inserts P into ML hierarchy with energy optimized bases

int ML_AGG_Gen_Prolongator_MandelMinEnergy(ML *ml,int level, int clevel, void *data)
{	
  int UseODEStrength = 1;

  if ( ml->comm->ML_nprocs > 1 )
    {	cerr << "ML_AGG_Gen_Prolongator_MinEnergy only works in serial" << endl;   exit(EXIT_FAILURE); }
	
  int         Ncoarse;
  ML_Operator **prev_P_tentatives;
  ML_Aggregate * ag = (ML_Aggregate *) data;
  double *Bzero;

#ifdef ML_TIMING
  double t0;
  t0 =  GetClock();
#endif

  ML_Operator* Amat = &(ml->Amat[level]); //already created and filled
  prev_P_tentatives = ag->P_tentative;    // for keep_P_tentative

  Amat->num_PDEs = ag->num_PDE_eqns;

  ML_Aggregate_Set_CurrentLevel(ag, level);

  //-----------------JBS----------------------------------------
  //Grab these values out of ag before it is overwritten during the P0 generation
  int numPDEs = ag->num_PDE_eqns;
  int NullDim = ag->nullspace_dim;

  //Save Nullspace_vect in agg from previous level, need for the 
  //satisfy constraints algorithm.  If this vector is not preset, 
  //assume it is a constant.

  Bzero =  (double *) ML_allocate(sizeof(double)*(NullDim*Amat->outvec_leng));
  if(ag->nullspace_vect == NULL) {
    if (NullDim != numPDEs) {
       cerr << "Null space not given but nullspace dimension not equal to the " <<endl;
       cerr << "number of PDEs: " << NullDim << " vs. "  << numPDEs << ". Cannot set default null space!" << endl;
       cerr << "(file " << __FILE__ << ", line " << __LINE__ << ")" << endl;
       exit(EXIT_FAILURE);
    }

    for(int i = 0; i < NullDim*Amat->outvec_leng; i++) Bzero[i] = 0.;

    for(int i = 0; i < Amat->outvec_leng; i++) {
      for (int jj=0; jj < numPDEs; jj++) {
         if ( (i%numPDEs) == jj ) 
            Bzero[i + jj*(Amat->outvec_leng)] = 1.0;
      }
    }
  }
  else {	
      memcpy( &(Bzero[0]), &(ag->nullspace_vect[0]), sizeof(double)*(ag->nullspace_dim)*(Amat->outvec_leng) ); 
  }
  //---------------End JBS------------------------------------

  // =============================================== //
  // Creates the non-smoothed prolongator, P0, then  //
  // Checks the dimension of the coarse problem.     //
  // Note that I always need P0 to form R.           //
  // =============================================== //

  ML_Operator* P0 = ML_Operator_Create(ml->comm);
  ML_Operator *Atilde;
  int time_steps = 38;
  double tfinal =  70.0;
  double theta = 0.0006;

  if (prev_P_tentatives == 0) {
    ag->P_tentative = ML_Operator_ArrayCreate(ag->max_levels);
    prev_P_tentatives = ag->P_tentative;
    for (int jj = 0 ; jj < ag->max_levels; jj++) prev_P_tentatives[jj] = 0;
  }

  if ((prev_P_tentatives != 0) && (prev_P_tentatives[clevel] != 0)) {
    P0 = prev_P_tentatives[clevel];
    Ncoarse = P0->invec_leng;
  }
  else {
    if (UseODEStrength) {
       Atilde = ML_ODE_Strength_Matrix(Amat, time_steps, tfinal, theta);
ML_Operator_Print(Atilde,"Atilde");
exit(1);
       Ncoarse  = ML_Aggregate_Coarsen(ag,Atilde,&P0,ml->comm);
    }
    else Ncoarse  = ML_Aggregate_Coarsen(ag,Amat,&P0,ml->comm);
    prev_P_tentatives[clevel] = P0;
  }

  /***************************************************************************/
  /* Begin JS Code for Mandel Energy Min Algorithm using CG                  */
  /***************************************************************************/

  //Declare Variables For The Algorithm
  // ---- Algorithm Parameters ----
  int Nits = 4;
  int UseF = 0;
  double tol = 1e-10;

  // ---- Workspace Variables ----
  char NN = 'N';
  char TT = 'T';
  char AA = 'A';
  ML_Operator *pk, *ap, *P, *rk, *minusA, *zk, *sparsity_pattern, *AbsA, *AbsP;
  int i, j, k, counter, counter2, rowstart, rowend, indx, indx2, dummy, flag;
  int numDOFs = Amat->outvec_leng;
  int numCoarseDOFs = P0->invec_leng;
  int numNodes = numDOFs/numPDEs;
  int NullDimSq = NullDim*NullDim;
  int * F = (int *) ML_allocate( sizeof(int)*numNodes*NullDim );
  int * Eye = (int *) ML_allocate( sizeof(int)*NullDim );
  double betak, newsum, oldsum = 1., alphak, resid;
  double * InnerProd = (double *) ML_allocate(sizeof(double)*Ncoarse);
  double * BtBinv = (double *) ML_allocate( sizeof(double)*numNodes*NullDim*NullDim );
  double * ABzero = (double *) ML_allocate( sizeof(double)*numDOFs*NullDim );
  double * B, *Scratch;
  double * U =  (double *) ML_allocate( sizeof(double)*NullDimSq);
  double *S =  (double *) ML_allocate( sizeof(double)*NullDim);
  double *VT =  (double *) ML_allocate( sizeof(double)*NullDimSq);
  double one = 1.0;
  double zero = 0.0;
  struct ML_CSR_MSRdata * SparsityCSRdata;

  for(i = 0; i < NullDim; i++) {	Eye[i] = 1; }

  pk = ML_Operator_Create(ml->comm); 
  ap = ML_Operator_Create(ml->comm);
  P  = ML_Operator_Create(ml->comm);
  rk = ML_Operator_Create(ml->comm);
  minusA = ML_Operator_Create(ml->comm);
  AbsA = ML_Operator_Create(ml->comm);
  AbsP = ML_Operator_Create(ml->comm);
  zk = ML_Operator_Create(ml->comm);
  sparsity_pattern = ML_Operator_Create(ml->comm);


  //Choose to use Amat or ODE Strength routine to generate strength matrix
  if (UseODEStrength) {
      //The ODE_Strength routines initializes Atilde
      //params		     Mat,  Time_Steps, t_final, drop_tol
//      Atilde = ML_ODE_Strength_Matrix(Amat, time_steps, tfinal, theta);
	
      //We want to wrap our operators inside of an abs() when calculating the sparsity pattern, 
      //	in order to mimic a symbolic mat-mat
      AbsA = ML_Operator_ImplicitAbs(Atilde, 0);
      AbsP = ML_Operator_ImplicitAbs(P0, 0);
	
      //Generate Sparsity Pattern
      ML_2matmult(AbsA, AbsP, sparsity_pattern, ML_CSR_MATRIX);
      SparsityCSRdata = (struct ML_CSR_MSRdata *) sparsity_pattern->data;
      ML_Operator_Destroy(&Atilde);
  }
  else {
      //We want to wrap our operators inside of an abs() when calculating the sparsity pattern, 
      //	in order to mimic a symbolic mat-mat
      AbsA = ML_Operator_ImplicitAbs(Amat, 0);
      AbsP = ML_Operator_ImplicitAbs(P0, 0);
	
      //Generate Sparsity Pattern
      ML_2matmult(AbsA, AbsP, sparsity_pattern, ML_CSR_MATRIX);
      SparsityCSRdata = (struct ML_CSR_MSRdata *) sparsity_pattern->data;
    }

  //Sparsity_Pattern may not have its column entries sorted per row, lets do that.
  ML_Sort_Cols(SparsityCSRdata, Amat->outvec_leng);

  //Get the inverse diagonal
  double * diagonal;			
  double * diagonal_local =  (double *) ML_allocate(sizeof(double)*(Amat->outvec_leng));
  ML_Operator_Getrow_Diag(Amat, &diagonal);
  for (int row = 0; row < Amat->outvec_leng; row++) 
    {      diagonal_local[row] = 1.0/diagonal[row];	}

  //Calculate Amat*Bzero, we need this to calculate F.  We need to know what nullspace components are 
  //preserved by the operator. 
  if(UseF) {
      //ABzero is stored in column major order
      for(i = 0; i < NullDim; i++)
	{
	  //Amat*Bzero(:,i) ===> ABzero(:,i)
	  ML_Operator_Apply(Amat, Amat->invec_leng, &(Bzero[i*numDOFs]), Amat->outvec_leng, &(ABzero[i*numDOFs]));
	}
  }

  //Calculate F.  F( node*NullDim ) is a NullDim sized binary vector that switches on and off
  //	the enforcement of constraints for each nullspace vector at current node.
  counter = 0;	//stores current dof offset
  counter2 = 0;   //stores current offset in F
  for(i = 0; i < numNodes; i++)
    {
      if(!UseF) //Don't Use F
	{	memcpy(&(F[i*NullDim]), &(Eye[0]), sizeof(int)*NullDim); }
      else	//Use F, Brave Soul.
	{
	  //memcpy(&(F[i*NullDim]), &(Eye[0]), sizeof(int)*NullDim);
	  for(j = 0; j < NullDim; j++)
	    {
	      flag = 0;
	      //for node i and nullspace vector j, we need to examine ABzero at each dof on node i
	      //if all the dofs on node i preserve this nullspace vector, we set F to 1.
	      for(k = 0; k < numPDEs; k++)
		{
		  if(fabs(ABzero[k + counter + j*numDOFs]) > 1e-10)
		    {flag = 1; break; }
		}
			
	      if(flag)
		{	F[counter2 + j] = 0; }	//Don't enforce constraints
	      else
		{	F[counter2 + j] = 1; } 	//Enforce constraints
	    }
	}
	
      counter += numPDEs;
      counter2 += NullDim;
    }
  //Free ABzero, we don't need it anymore
  ML_free(ABzero);

  //Calculate inv(B^T * B)_i using the SVD.  B is the nullspace vectors 
  //	restricted to each node i, where each node could have multiple dofs 
  //	depending on the numPDEs
  counter2 = 6*NullDim;
  Scratch =  (double *) ML_allocate( sizeof(double)*counter2);
  indx2 = 0;
  for(i = 0; i < numNodes; i++)
    {
      //We use the sparsity pattern for this node to decide which entries of the NullSpace vectors to use
      //	We assume a unifrom sparsity pattern for each dof on a node
      rowstart = SparsityCSRdata->rowptr[i*numPDEs];
      rowend = SparsityCSRdata->rowptr[i*numPDEs + 1];
      int length = rowend - rowstart;
	
      //Build B in scratch memory in column major storage
      counter = 0;
      B = (double *) ML_allocate( sizeof(double)*length*NullDim );
      for(k = 0; k < NullDim; k++)  //loop over cols
	{
	  //Use F(i) as a switch to either set this column to zero or set it to nullspace_vect, 
	  //	i.e. binary switch for enforcing constraints
	  if(F[i*NullDim + k] == 1)
	    {
	      for(j = rowstart; j < rowend; j++)  //loop over rows
		{
		  indx = SparsityCSRdata->columns[j] + (k*numCoarseDOFs);
		  B[counter] = ag->nullspace_vect[indx];
		  counter++;
		}
	    }
	  else // don't enforce constraint
	    {
	      for(j = 0; j < length; j++)	//loop over rows
		{	B[counter] = 0.0; counter++; }
	    }
	}
	
      //Calculate BtB with BLAS, set option to transpose first B
      //	   X = Y = B, and store result in Z=B
      // params: X',  Y,  op(X)rows, op(Y)columns, op(X)columns, ALPHA,  X,       Xrows,   Y,       Yrows,   BETA,   Z,                Zrows  
      DGEMM_F77(&TT, &NN, &NullDim,   &NullDim,    &length,      &one,   &(B[0]), &length, &(B[0]), &length, &zero,  &(BtBinv[indx2]), &NullDim);

      //BtB may be singular, depending on F.  So calculate a pseudo-inverse based on the SVD.
      //Allocate the Scratch space for the SVD routine, which according to the LAPACK webpage should be >= 5*(dimension of mat)
      //Calculate SVD of BtB, which is NullDim X NullDim
      //params: AA => Full-SVD     Rows      Cols           Mat          LDA     S  U    LDU     VT    LDVT             Size-Scratch      Flag
      DGESVD_F77( &AA, &AA,      &NullDim, &NullDim, &(BtBinv[indx2]), &NullDim, S, U, &NullDim, VT, &NullDim, Scratch,   &counter2,      &dummy);
	
      //Filter Sigma for numerically zero entries and invert.  LAPACK does guarantee that first entry of Sigma is the largest.
      if(fabs(S[0]) < 1e-10)
	{
	  for(k = 0; k < NullDim; k++)
	    {	S[k] = 0.0; }
	}
      else
	{
	  for(k = 1; k < NullDim; k++)
	    {
	      if(fabs(S[k]/S[0]) < 1e-10)
		{	S[k] = 0.0; }
	      else
		{	S[k] = 1/S[k]; }
	    }
	  S[0] = 1/S[0];
	}

      //SigmaInv*U^T, i.e. just scale each column of U by the corresponding entry of SigmaInv	
      counter = 0;
      for(k = 0; k < NullDim ; k++)  //loop over cols
	{
	  for(j = 0; j < NullDim; j++) //loop over rows
	    {
	      U[counter] *= S[k];
	      counter++;	
	    }
	}
	
      //V*(SigmaInv*U^T)  ==>  &(BtBinv[indx2]).   V and U are NullDim X NullDim 
      // params: X',  Y', op(X)rows, op(Y)columns, op(X)columns,   ALPHA,    X,       Xrows,     Y,      Yrows,   BETA,     Z,                Zrows  
      DGEMM_F77(&TT, &TT, &NullDim,    &NullDim,     &NullDim,     &one,   &(VT[0]), &NullDim, &(U[0]), &NullDim, &zero,  &(BtBinv[indx2]), &NullDim);
	
      indx2 += NullDimSq;
	
      if(B != NULL) ML_free(B);
    }
  ML_free(U); ML_free(VT); ML_free(S);
  if(Scratch != NULL) ML_free(Scratch);


  //Begin Preparations for CG Iterations
  minusA = ML_Operator_ImplicitlyScale(Amat, -1.0, 0);    //  wrap A so that minusA = -A
  ML_2matmult(minusA, P0, rk, ML_CSR_MATRIX);	        //  rk = -A*P0

  //Enforce constraints on rk, but first 
  //	rk may not have its column entries sorted per row, lets do that.
  //	Experimentally, the unsorted nature of some row's column indices was verified.
  ML_Sort_Cols((struct ML_CSR_MSRdata *) rk->data, Amat->outvec_leng);
  ML_Enforce_Sparsity(rk, SparsityCSRdata);
  ML_Satisfy_Constraints(rk, sparsity_pattern, ag->nullspace_vect, BtBinv, F, numPDEs, numDOFs, numNodes, NullDim);

  resid = ML_MaxEntry(rk);
  cout << printf("Iteration 0 ------ Max(Abs( R ) = %f", resid) << endl;

  //Begin CG iterations to minimize the energy of P
  i = 1;
  while( (i <= Nits) && (resid > tol) )
    {
      //Implicitly scale rk with Dinv for for the application of the preconditioner
      zk = ML_Operator_ImplicitlyVScale(rk, diagonal_local, 0);  
	
      //Calculate the innerprodct of rk and zk, which are "vectors" in the CG sense 
      //but are represented here as matrices.
      //Note that we scale rk by Dinv in the multiply_all_vscale routine, so that 
      //we are really "passing" zk and rk to multiply_all.
      ML_multiply_all_vscale(rk, rk, &(InnerProd[0]), diagonal_local);	
      newsum = 0.0;
      for(int j = 0; j < Ncoarse; j++)
	{	newsum += InnerProd[j]; }	

      //pk = zk;
      if(i == 1)
	{	
	  ML_Operator_Add(zk, zk, pk, ML_CSR_MATRIX, 0.0);
	}
      else // pk = zk + betak*pk
	{
	  betak = newsum/oldsum;
	  ML_Operator_Add(zk, pk, pk, ML_CSR_MATRIX, betak);
	}
      oldsum = newsum;	
	
      //	cout << "newsum = " << newsum << endl;
      //	cout << "pk(1,1) = " << ((struct ML_CSR_MSRdata *) pk->data)->values[0] << endl;

      //ap = Amat*pk
      if(i > 1)
	{	
	  ML_Operator_Destroy(&ap);
	  ap = ML_Operator_Create(ml->comm);
	}
      ML_2matmult(Amat, pk, ap, ML_CSR_MATRIX); 
      //	ML_Operator_Print(ap,"ap1_before");

      //	Enforce constraints on ap
      //ap may not have its column entries sorted per row.  This was 
      //	experimentally verified to be a problem
      ML_Sort_Cols( (struct ML_CSR_MSRdata *) ap->data, Amat->outvec_leng);
      ML_Enforce_Sparsity(ap, (struct ML_CSR_MSRdata *) sparsity_pattern->data);
      ML_Satisfy_Constraints(ap, sparsity_pattern, ag->nullspace_vect, BtBinv, F, numPDEs, numDOFs, numNodes, NullDim);

      //Alpha -- Calculate the innerprodct of pk and ap, which are "vectors" in the CG sense 
      //but are represented here as matrices.  pk, ap, zk and rk are all the 
      //same size, so InnerProd will just be written over.
      ML_multiply_all(pk, ap, &(InnerProd[0])); 
      alphak = 0.0;
      for(int j = 0; j < Ncoarse; j++)
	{	alphak += InnerProd[j]; }
      //	ML_Operator_Print(pk,"pk1");
      //	ML_Operator_Print(ap,"ap1");
      alphak = newsum/alphak;
	
      //P = P + alpha*pk 
      if(i == 1)
	{	ML_Operator_Add(P0, pk, P, ML_CSR_MATRIX, alphak); }
      else
	{	ML_Operator_Add(P, pk, P, ML_CSR_MATRIX, alphak); }
	
      //rk = rk - alpha*ap
      ML_Operator_Add(rk, ap, rk, ML_CSR_MATRIX, -alphak);
	
      resid = ML_MaxEntry(rk);
      cout << printf("Iteration %d ------ Max(Abs( R )) = %f", i, resid) << endl;
      i++;
    }

  ML_Operator_Add(P, ap, &(ml->Pmat[clevel]),         // Add the matrix in P to 0
		  ML_CSR_MATRIX, 0.0 );              // and store the result
  // in the prolongator hierarchy
  // array as a CSR matrix.

  ML_Operator_Set_1Levels(&(ml->Pmat[clevel]),        // This piece of code is
			  &(ml->SingleLevel[clevel]), // needed for grid transfer
			  &(ml->SingleLevel[level])); // operators. It tells ML
  // what array elements 
  // correspond to the current
  // fine level and to the 
  // current coarse level.


  //Deconstruct operators and vectors
  if(diagonal_local != NULL) ML_free(diagonal_local);
  if(InnerProd != NULL) ML_free(InnerProd);
  if(Bzero != NULL) ML_free(Bzero);
  ML_free(Eye);
  ML_free(F);
  ML_free(BtBinv);
  ML_Operator_Destroy(&pk);
  ML_Operator_Destroy(&ap);
  ML_Operator_Destroy(&P);
  ML_Operator_Destroy(&rk);
  ML_Operator_Destroy(&minusA);
  ML_Operator_Destroy(&zk);
  ML_Operator_Destroy(&AbsA);
  ML_Operator_Destroy(&AbsP);
  ML_Operator_Destroy(&sparsity_pattern);

  return(0); 
}
