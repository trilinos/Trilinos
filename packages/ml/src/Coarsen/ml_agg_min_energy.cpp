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

// ====================================================================== 
// Take each column of Op, compute its 2-norm squared, and store
// in Column2Norm[]. 
//
// NOTE: Column2Norm[] must have enough space to store ghost unknowns.
//
// Note: Must a CSR matrix but it is easy to change for a generic getrow()
//
static void multiply_self_all(ML_Operator* Op, double* Column2Norm)
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
inline static void multiply_all(ML_Operator* left, ML_Operator* right,
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

    multiply_all(P0, DinvAP0, &Numerator[0]);
    multiply_self_all(DinvAP0, &Denominator[0]);
    break;

  case 2:
    // Minimize with respect to the (D^{-1} A)' D^{-1} A norm. 
    // Need to be smart here to avoid the construction of A' A
    //
    //                   diag( P0' ( A' D^{-1}' D^{-1} A) D^{-1} A P0)
    //   omega =   --------------------------------------------------------
    //             diag( P0' A' D^{-1}' ( A' D^{-1}' D^{-1} A) D^{-1} A P0)
    //
    multiply_all(DinvAP0_subset, DinvADinvAP0, &Numerator[0]);
    multiply_self_all(DinvADinvAP0, &Denominator[0]);

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
    multiply_all(P0, DinvADinvAP0,  &Numerator[0]);
    multiply_self_all(DinvAP0, &tmp[0]);

    for (int i = 0 ; i < n_0 ; ++i) Numerator[i] += tmp[i];
    multiply_all(DinvAP0, DinvADinvAP0, &Denominator[0]);
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

  for (int i = 0 ; i < NComputedOmegas ; ++i) {
    ColBasedOmega[i] = Numerator[i]/Denominator[i];
    double& val = ColBasedOmega[i];
    if (val < 0.0) {
      val = 0.0;
      ++zero_local;
    }
    if (val < min_local) min_local = val;
    if (val > max_local) max_local = val;
  }
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
  }

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

      Amat = ML_Operator_ImplicitlyVCScale(Amat, Dinv, 0);

      ML_2matmult(P0_trans, Amat, P0TA, ML_CSR_MATRIX);

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
