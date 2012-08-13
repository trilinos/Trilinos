//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER

#include "AztecOO.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#endif
#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_MsrMatrix.h"
#include "AztecOO_Scaling.h"
#include "AztecOO_StatusTestMaxIters.h"
#include "AztecOO_StatusTestResNorm.h"
#include "AztecOO_StatusTestCombo.h"
#include "AztecOO_Operator.h"
#include "AztecOO_Version.h"

bool argument_is_present(const char* argument,
                         int argc,
                         char** argv);

double resid2norm(Epetra_CrsMatrix& A,
                  Epetra_Vector& x,
                  Epetra_Vector& b);

int test_azoo_as_precond_op(Epetra_CrsMatrix& A,
                            Epetra_Vector& x,
                            Epetra_Vector& b,
                            bool verbose);

int test_azoo_conv_anorm(Epetra_CrsMatrix& A,
                         Epetra_Vector& x,
                         Epetra_Vector& b,
                         bool verbose);

int test_azoo_conv_with_scaling(int conv_option, int scaling_option,
                                const Epetra_Comm& comm, bool verbose);

int test_azoo_with_ilut(Epetra_CrsMatrix& A,
                        Epetra_Vector& x,
                        Epetra_Vector& b,
                        bool verbose);

int test_azoo_scaling(Epetra_CrsMatrix& A,
                      Epetra_Vector& x,
                      Epetra_Vector& b,
                      bool verbose);

int call_AZ_iterate(AZ_MATRIX* Amat,
                    AZ_PRECOND* P,
                    AZ_SCALING* S,
                    double* x,
                    double* b,
                    int* options,
                    double* params,
                    double* status,
                    int* proc_config,
                    int keep_info,
                    int pre_calc,
                    bool verbose);

int create_and_transform_simple_matrix(int matrix_type,
                                    int N,
                                    double diag_term,
                                    int* proc_config,
                                    AZ_MATRIX*& Amat,
                                    int*& external,
                                    int*& update_index,
                                    int*& external_index);

int test_AZ_iterate_AZ_pre_calc_AZ_reuse(Epetra_Comm& Comm,
                                         int* options,
                                         bool verbose);

int test_AZ_iterate_then_AZ_scale_f(Epetra_Comm& Comm, bool verbose);

int test_bug2554(Epetra_Comm& Comm, bool verbose);

int test_bug2890(Epetra_Comm& Comm, bool verbose);

Epetra_CrsMatrix* create_and_fill_crs_matrix(const Epetra_Map& emap);

void destroy_matrix(AZ_MATRIX*& Amat);

int main(int argc, char *argv[])
{
#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  int numprocs = comm.NumProc();
  int localproc = comm.MyPID();

  ///////////////////////////////////////////////////
  //First figure out whether verbose output is on.
  //(and if so, only turn it on for proc 0)

  bool verbose = false;
  if (localproc == 0) {
    verbose = argument_is_present("-v", argc, argv);
  }

  ///////////////////////////////////////////////////

  int local_n = 30;
  long long global_n = ((long long)numprocs)*((long long)local_n);

  Epetra_Map emap(global_n, 0, comm);

  Epetra_CrsMatrix* A = create_and_fill_crs_matrix(emap);

  Epetra_Vector x(emap), b(emap);

  x.PutScalar(1.0);

  A->Multiply(false, x, b);
  x.PutScalar(0.0);

  double initial_norm = resid2norm(*A, x, b);

  if (verbose) {
    cout << "Initial 2-norm of b-A*x: "<<initial_norm<<endl;
  }

  int err = test_azoo_as_precond_op(*A, x, b, verbose);
  if (err != 0) {
    cout << "test_azoo_as_precond_op err, test FAILED."<<endl;
    return(err);
  }

  err = test_azoo_conv_anorm(*A, x, b, verbose);
  if (err != 0) {
    cout << "test_azoo_conv_anorm err, test FAILED."<<endl;
    return(err);
  }

  if (verbose) {
    cout << "testing AztecOO with AZ_conv = AZ_Anorm, and AZ_scaling = AZ_sym_diag" << endl;
  }

  err = test_azoo_conv_with_scaling(AZ_Anorm, AZ_sym_diag, A->Comm(), verbose);
  if (err != 0) {
    cout << "test_azoo_conv_with_scaling err, test FAILED."<<endl;
    return(err);
  }

  if (verbose) {
    cout << "testing AztecOO with AZ_conv = AZ_rhs, and AZ_scaling = AZ_sym_diag" << endl;
  }

  err = test_azoo_conv_with_scaling(AZ_rhs, AZ_sym_diag, A->Comm(), verbose);
  if (err != 0) {
    cout << "test_azoo_conv_with_scaling err, test FAILED."<<endl;
    return(err);
  }

  err = test_azoo_with_ilut(*A, x, b, verbose);
  if (err != 0) {
    cout << "test_azoo_with_ilut err, test FAILED."<<endl;
    return(err);
  }

  err = test_azoo_scaling(*A, x, b, verbose);
  if (err != 0) {
    cout << "test_azoo_scaling err="<<err<<", test FAILED."<<endl;
    return(err);
  }

  delete A;

  int* options = new int[AZ_OPTIONS_SIZE];
  options[AZ_solver] = AZ_cg;
  options[AZ_subdomain_solve] = AZ_none;
  options[AZ_precond] = AZ_Jacobi;

  if (verbose)
    std::cout << "about to call test_AZ_iterate_AZ_pre_calc_AZ_reuse"
       <<std::endl;

  err = test_AZ_iterate_AZ_pre_calc_AZ_reuse(comm, options, verbose);
  if (err != 0) {
    cout << "test_AZ_iterate_AZ_pre_calc_AZ_reuse err, test FAILED."<<endl;
    return(err);
  }

  options[AZ_solver] = AZ_cgs;
  options[AZ_subdomain_solve] = AZ_icc;
  options[AZ_precond] = AZ_dom_decomp;

  err = test_AZ_iterate_AZ_pre_calc_AZ_reuse(comm, options, verbose);
  if (err != 0) {
    cout << "test_AZ_iterate_AZ_pre_calc_AZ_reuse err, test FAILED."<<endl;
    return(err);
  }

  options[AZ_solver] = AZ_gmres;
  options[AZ_subdomain_solve] = AZ_ilut;

  err = test_AZ_iterate_AZ_pre_calc_AZ_reuse(comm, options, verbose);
  if (err != 0) {
    cout << "test_AZ_iterate_AZ_pre_calc_AZ_reuse err, test FAILED."<<endl;
    return(err);
  }

  options[AZ_solver] = AZ_tfqmr;
  options[AZ_subdomain_solve] = AZ_ilu;

  err = test_AZ_iterate_AZ_pre_calc_AZ_reuse(comm, options, verbose);
  if (err != 0) {
    cout << "test_AZ_iterate_AZ_pre_calc_AZ_reuse err, test FAILED."<<endl;
    return(err);
  }

  options[AZ_solver] = AZ_bicgstab;
  options[AZ_subdomain_solve] = AZ_rilu;

  err = test_AZ_iterate_AZ_pre_calc_AZ_reuse(comm, options, verbose);
  if (err != 0) {
    cout << "test_AZ_iterate_AZ_pre_calc_AZ_reuse err, test FAILED."<<endl;
    return(err);
  }

  err = test_AZ_iterate_then_AZ_scale_f(comm, verbose);
  if (err != 0) {
    cout << "test_AZ_iterate_then_AZ_scale_f err, test FAILED."<<endl;
    return(err);
  }

  delete [] options;

  err = test_bug2554(comm, verbose);
  if (err != 0) {
    cout << "test_bug2554 err, test FAILED."<<endl;
    return(err);
  }

  err = test_bug2890(comm, verbose);
  if (err != 0) {
    cout << "test_bug2890 err, test FAILED."<<endl;
    return(err);
  }

  cout << "********* Test passed **********" << endl;

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

  return(0);
}

bool argument_is_present(const char* argument,
                         int argc,
                         char** argv)
{
  if (argument == NULL || argc < 1) return(false);

  for(int i=0; i<argc; ++i) {
    if (strcmp(argument, argv[i]) == 0) {
      return(true);
    }
  }
  return(false);
}

double resid2norm(Epetra_CrsMatrix& A,
                  Epetra_Vector& x,
                  Epetra_Vector& b)
{
  Epetra_Vector r(x);
  A.Multiply(false, x, r);

  r.Update(1.0, b, -1.0);//r  =  b - r  =  b - A*x

  double nrm2;
  r.Norm2(&nrm2);

  return(nrm2);
}

Epetra_CrsMatrix* create_and_fill_crs_matrix(const Epetra_Map& emap)
{
  int localproc = emap.Comm().MyPID();
  int local_n = emap.NumMyElements();
  long long global_n = emap.NumGlobalElements64();
  long long myFirstGlobalRow = ((long long)localproc)*((long long)local_n);
  long long globalCols[3];
  double values[3];
  Epetra_CrsMatrix* A = new Epetra_CrsMatrix(Copy, emap, 3);

  for(int i=0; i<local_n; ++i) {
    long long globalRow = myFirstGlobalRow +i;

    int numcols = 0;
    if (globalRow > 0) {
      globalCols[numcols] = globalRow-1;
      values[numcols++] = -1.0;
    }

    globalCols[numcols] = globalRow;
    values[numcols++] = 4.0;

    if (globalRow < global_n-1) {
      globalCols[numcols] = globalRow+1;
      values[numcols++] = -1.0;
    }

    A->InsertGlobalValues(globalRow, numcols, values, globalCols);
  }

  A->FillComplete();

  return A;
}

int test_azoo_as_precond_op(Epetra_CrsMatrix& A,
                            Epetra_Vector& x,
                            Epetra_Vector& b,
                            bool verbose)
{
  AztecOO* azoo0 = new AztecOO(&A, &x, &b);

  azoo0->SetAztecOption(AZ_solver, AZ_cg);
  azoo0->SetAztecOption(AZ_output, AZ_none);
  azoo0->SetAztecOption(AZ_conv, AZ_none);

  AztecOO_Operator azooOp(azoo0, 10);

  Epetra_LinearProblem elp(&A, &x, &b);

  x.PutScalar(0.0);
  AztecOO* azoo1 = new AztecOO(elp);
  azoo1->SetPrecOperator(&azooOp);
  azoo1->SetAztecOption(AZ_solver, AZ_gmres);
  azoo1->SetAztecOption(AZ_output, AZ_none);
  azoo1->SetAztecOption(AZ_conv, AZ_none);

  if (verbose) {
    cout << "testing recursive solve (AztecOO as"
    << " preconditioner for another AztecOO)."<<endl;
  }

  int maxiters = 100;
  double tolerance = 1.e-12;
  azoo1->Iterate(maxiters, tolerance);

  double resid1 = resid2norm(A, x, b);
  if (verbose) {
    cout << "residual 2-norm after recursive solve: "
        << resid1 << endl;
  }

  if (resid1 > 1.e-6) {
    return(-1);
  }

  if (verbose) {
    cout << "now make sure the precond AztecOO instance"
      << " hasn't been corrupted."<<endl;
  }

//  AZ_manage_memory(0, -43, 0, 0, 0);

  x.PutScalar(0.0);
  azoo0->Iterate(&A, &x, &b, maxiters, tolerance);

  double resid0 = resid2norm(A, x, b);
  if (verbose) {
    cout << "residual 2-norm: " << resid0 << endl;
  }
  if (resid0 > 1.e-6) {
    return(-1);
  }

  delete azoo0;
  delete azoo1;

  return(0);
}

int test_azoo_with_ilut(Epetra_CrsMatrix& A,
                        Epetra_Vector& x,
                        Epetra_Vector& b,
                        bool verbose)
{
  x.PutScalar(0.0);

  AztecOO* azoo0 = new AztecOO(&A, &x, &b);

  azoo0->SetAztecOption(AZ_solver, AZ_gmres);
  azoo0->SetAztecOption(AZ_output, AZ_none);
  azoo0->SetAztecOption(AZ_conv, AZ_Anorm);
  azoo0->SetAztecOption(AZ_precond, AZ_dom_decomp);
  azoo0->SetAztecOption(AZ_subdomain_solve, AZ_ilut);
  azoo0->SetAztecOption(AZ_keep_info, 1);

  if (verbose) {
    cout << "testing AztecOO with GMRES and ILUT, AZ_keep_info==1" << endl;
  }

  int maxiters = 100;
  double tolerance = 1.e-12;
  int err = azoo0->Iterate(maxiters, tolerance);
  if (err != 0) {
    if (verbose) cout << "AztecOO::Iterate return err="<<err<<endl;
    return(err);
  }

  double resid = resid2norm(A, x, b);
  if (verbose) {
    cout << "residual 2-norm after GMRES/ILUT solve: "
        << resid << endl;
  }

  if (resid > 1.e-6) {
    return(-1);
  }

  if (verbose) {
    cout << "solving with GMRES/ILUT again, AZ_pre_calc==AZ_reuse"
        << endl << "(will error out if factors weren't kept from"
         << " previous solve)"<<endl;
  }

  azoo0->SetAztecOption(AZ_pre_calc, AZ_reuse);
  x.PutScalar(0.0);
  err = azoo0->Iterate(maxiters, tolerance);
  if (err != 0) {
    if (verbose) cout << "AztecOO::Iterate return err="<<err<<endl;
    return(err);
  }

  double resid2 = resid2norm(A, x, b);
  if (verbose) {
    cout << "after second GMRES/ILUT solve, residual 2-norm: "
      << resid2 <<endl;
  }
  if (resid2 > 1.e-6) {
    return(-1);
  }

  if (verbose) {
    cout << "azoo0->SolveTime(): " << azoo0->SolveTime() << endl;
  }

  delete azoo0;

  return(0);
}

int test_azoo_conv_with_scaling(int conv_option, int scaling_option,
                                const Epetra_Comm& comm, bool verbose)
{
  int localN = 20;
  int numprocs = comm.NumProc();
  long long globalN = ((long long)numprocs)*((long long)localN);
 
  Epetra_Map emap(globalN, 0, comm);
  Epetra_CrsMatrix* Acrs = create_and_fill_crs_matrix(emap);

  Epetra_Vector x_crs(emap), b_crs(emap);
  x_crs.PutScalar(1.0);

  Acrs->Multiply(false, x_crs, b_crs);
  x_crs.PutScalar(0.0);

  AztecOO azoo(Acrs, &x_crs, &b_crs);
  azoo.SetAztecOption(AZ_conv, conv_option);
  azoo.SetAztecOption(AZ_solver, AZ_cg);
  azoo.SetAztecOption(AZ_scaling, scaling_option);

  azoo.Iterate(100, 1.e-9);

  //now, do the same thing with 'old-fashioned Aztec', and compare
  //the solutions.

    int* proc_config = new int[AZ_PROC_SIZE];

#ifdef EPETRA_MPI
  AZ_set_proc_config(proc_config, MPI_COMM_WORLD);
  AZ_set_comm(proc_config, MPI_COMM_WORLD);
#else
  AZ_set_proc_config(proc_config, 0);
#endif

  int *external, *update_index, *external_index;
  int *external2, *update_index2, *external_index2;
  AZ_MATRIX* Amsr = NULL;
  AZ_MATRIX* Avbr = NULL;
  int err = create_and_transform_simple_matrix(AZ_MSR_MATRIX, localN, 4.0,
                                               proc_config, Amsr,
                                               external, update_index,
                                               external_index);

  int N_update = localN+Amsr->data_org[AZ_N_border];
  double* x_msr = new double[N_update];
  double* b_msr = new double[N_update*2];
  double* b_msr_u = b_msr+N_update;
  double* x_vbr = new double[N_update];
  double* b_vbr = new double[N_update*2];
  double* b_vbr_u = b_vbr+N_update;

  err = create_and_transform_simple_matrix(AZ_VBR_MATRIX, localN, 4.0,
                                           proc_config, Avbr,
                                           external2, update_index2,
                                           external_index2);
  for(int i=0; i<N_update; ++i) {
    x_msr[i] = 1.0;
    b_msr[i] = 0.0;
    b_msr_u[i] = 0.0;
    x_vbr[i] = 1.0;
    b_vbr[i] = 0.0;
    b_vbr_u[i] = 0.0;
  }

  Amsr->matvec(x_msr, b_msr, Amsr, proc_config);
  Avbr->matvec(x_vbr, b_vbr, Avbr, proc_config);

  for(int i=0; i<N_update; ++i) {
    x_msr[i] = 0.0;
    x_vbr[i] = 0.0;
  }

  //check that the rhs's are the same.
  double max_rhs_diff1 = 0.0;
  double max_rhs_diff2 = 0.0;
  double* bptr_crs = b_crs.Values();

  AZ_invorder_vec(b_msr, Amsr->data_org, update_index, NULL, b_msr_u);
  AZ_invorder_vec(b_vbr, Avbr->data_org, update_index2, Avbr->rpntr, b_vbr_u);
  for(int i=0; i<localN; ++i) {
    if (std::abs(bptr_crs[i] - b_msr_u[i]) > max_rhs_diff1) {
      max_rhs_diff1 = std::abs(bptr_crs[i] - b_msr_u[i]);
    }
    if (std::abs(bptr_crs[i] - b_vbr_u[i]) > max_rhs_diff2) {
      max_rhs_diff2 = std::abs(bptr_crs[i] - b_vbr_u[i]);
    }
  }

  if (max_rhs_diff1> 1.e-12) {
    cout << "AztecOO rhs not equal to Aztec msr rhs "<<max_rhs_diff1<<endl;
    return(-1);
  }

  if (max_rhs_diff2> 1.e-12) {
    cout << "AztecOO rhs not equal to Aztec vbr rhs "<<max_rhs_diff2<<endl;
    return(-1);
  }

  int* az_options = new int[AZ_OPTIONS_SIZE];
  double* params = new double[AZ_PARAMS_SIZE];
  double* status = new double[AZ_STATUS_SIZE];
  AZ_defaults(az_options, params);
  az_options[AZ_solver] = AZ_cg;
  az_options[AZ_conv] = conv_option;
  az_options[AZ_scaling] = scaling_option;

  az_options[AZ_max_iter] = 100;
  params[AZ_tol] = 1.e-9;

  AZ_iterate(x_msr, b_msr, az_options, params, status, proc_config,
             Amsr, NULL, NULL);
  AZ_iterate(x_vbr, b_vbr, az_options, params, status, proc_config,
             Avbr, NULL, NULL);

  AZ_invorder_vec(x_msr, Amsr->data_org, update_index, NULL, b_msr_u);
  AZ_invorder_vec(x_vbr, Avbr->data_org, update_index2, Avbr->rpntr, b_vbr_u);

  double max_diff1 = 0.0;
  double max_diff2 = 0.0;
  double* xptr_crs = x_crs.Values();

  for(int i=0; i<localN; ++i) {
    if (std::abs(xptr_crs[i] - b_msr_u[i]) > max_diff1) {
      max_diff1 = std::abs(xptr_crs[i] - b_msr_u[i]);
    }
    if (std::abs(xptr_crs[i] - b_vbr_u[i]) > max_diff2) {
      max_diff2 = std::abs(xptr_crs[i] - b_vbr_u[i]);
    }
  }

  if (max_diff1 > 1.e-7) {
    cout << "AztecOO failed to match Aztec msr with scaling and Anorm conv."
      << endl;
    return(-1);
  }

  if (max_diff2 > 1.e-7) {
    cout << "AztecOO failed to match Aztec vbr with scaling and Anorm conv."
      << endl;
    return(-1);
  }

  delete Acrs;
  delete [] x_msr;
  delete [] b_msr;
  delete [] x_vbr;
  delete [] b_vbr;
  destroy_matrix(Amsr);
  destroy_matrix(Avbr);
  delete [] proc_config;
  free(update_index);
  free(external);
  free(external_index);
  free(update_index2);
  free(external2);
  free(external_index2);
  delete [] az_options;
  delete [] params;
  delete [] status;

  return(0);
}

int test_azoo_conv_anorm(Epetra_CrsMatrix& A,
                         Epetra_Vector& x,
                         Epetra_Vector& b,
                         bool verbose)
{
  if (verbose) {
    cout << "testing AztecOO with AZ_conv = AZ_Anorm" << endl;
  }

  Epetra_Vector soln_Anorm(x), soln_none(x), vec1(x), rhs(x);

  //We'll put large numbers in a vector and use that to generate an rhs
  //which has norm much larger than the infinity-norm of the matrix.
  vec1.PutScalar(1000.0);
  soln_Anorm.PutScalar(0.0);
  soln_none.PutScalar(0.0);

  A.Multiply(false, vec1, rhs);

  AztecOO azoo(&A, &soln_Anorm, &rhs);
  azoo.SetAztecOption(AZ_conv, AZ_Anorm);
  azoo.SetAztecOption(AZ_solver, AZ_cg);
  //azoo.SetAztecOption(AZ_scaling, AZ_sym_diag);

  azoo.Iterate(30, 1.e-5);

  AztecOO azoo1(&A, &soln_none, &rhs);
  azoo1.SetAztecOption(AZ_conv, AZ_rhs);
  azoo1.SetAztecOption(AZ_solver, AZ_cg);

  azoo1.Iterate(30, 1.e-5);

  double rhsnorm = 0.0;
  rhs.Norm2(&rhsnorm);

  double anorm = A.NormInf();

  double rnrm_anorm = resid2norm(A, soln_Anorm, rhs);
  double rnrm_rhs = resid2norm(A, soln_none, rhs);

  //we expect the ratio rnrm_anorm/rnrm_rhs to roughly equal
  //the ratio anorm/rhsnorm, since rnrm_anorm is the residual norm
  //obtained for the solve that used Anorm scaling to determine
  //convergence, and rnrm_rhs is the residual norm obtained for
  //the solve that used rhs scaling.

  double ratio1 = rnrm_anorm/rnrm_rhs;
  double ratio2 = anorm/rhsnorm;

  cout << "ratio1: " << ratio1 << ", ratio2: " << ratio2 << endl;
  if (std::abs(ratio1 - ratio2) > 1.e-1) {
    if (verbose) {
      cout << "anorm: " << anorm << ", rhsnorm: " << rhsnorm
       << "rnrm_anorm: " << rnrm_anorm << ", rnrm_rhs: " << rnrm_rhs<<endl;
    }
    return(-1);
  }

  return(0);
}

int test_azoo_scaling(Epetra_CrsMatrix& A,
                      Epetra_Vector& x,
                      Epetra_Vector& b,
                      bool verbose)
{
  Epetra_Vector vec1(x);
  Epetra_Vector vec2(x);
  Epetra_Vector diag(x);
  Epetra_Vector vec3(x);
  Epetra_Vector vec4(x);
  Epetra_Vector rhs(x);
  Epetra_Vector soln_none(x);
  Epetra_Vector soln_jacobi(x);
  Epetra_Vector soln_rowsum(x);
  Epetra_Vector soln_symdiag(x);

  vec1.PutScalar(1.0);

  A.Multiply(false, vec1, vec2);

  A.ExtractDiagonalCopy(diag);

  double* diag_vals = NULL;
  diag.ExtractView(&diag_vals);

  int* options = new int[AZ_OPTIONS_SIZE];
  double* params = new double[AZ_PARAMS_SIZE];
  AZ_defaults(options, params);

  options[AZ_output] = verbose ? 1 : AZ_none;

  options[AZ_scaling] = AZ_Jacobi;
  AztecOO::MatrixData mdata(&A);
  AZ_MATRIX* Amat = AZ_matrix_create(vec1.Map().NumMyElements());
  AZ_set_MATFREE(Amat, (void*)(&mdata), Epetra_Aztec_matvec);

  AZ_SCALING* scaling = AZ_scaling_create();

  double* xvals = NULL, *bvals = NULL;
  x.ExtractView(&xvals);
  b.ExtractView(&bvals);

  int err = AztecOO_scale_epetra(AZ_SCALE_MAT_RHS_SOL, Amat,
                                 options, bvals, xvals, NULL, scaling);
  if (err != 0) {
    if (verbose) {
      cout << "AztecOO_scale_epetra returned err="<<err<<endl;
    }
    return(err);
  }

  A.Multiply(false, vec1, vec3);

  vec4.Multiply(1.0, diag, vec3, 0.0);

  double vec2nrm, vec4nrm;

  vec2.Norm2(&vec2nrm);
  vec4.Norm2(&vec4nrm);

  if (fabs(vec2nrm - vec4nrm) > 1.e-6) {
    return(-1);
  }

  //now call the scaling function again, just to allow for
  //testing memory-leak issues.
  err = AztecOO_scale_epetra(AZ_SCALE_MAT_RHS_SOL, Amat,
                             options, bvals, xvals, NULL, scaling);
  if (err != 0) {
    if (verbose) {
      cout << "AztecOO_scale_epetra returned err="<<err<<endl;
    }
    return(err);
  }

  AztecOO_scale_epetra(AZ_DESTROY_SCALING_DATA, Amat, options,
                       bvals, xvals, NULL, scaling);

  x.PutScalar(1.0);

  Epetra_CrsMatrix* Atmp = create_and_fill_crs_matrix(A.RowMap());
  Atmp->Multiply(false, x, rhs);

  x.PutScalar(0.0);

  AztecOO azoo(&A, &x, &b);

  azoo.SetAztecOption(AZ_scaling, AZ_Jacobi);
  if (verbose) {
    azoo.SetAztecOption(AZ_output, 1);
  }
  else {
    azoo.SetAztecOption(AZ_output, AZ_none);
  }

  azoo.Iterate(100, 1.e-6);

  delete Atmp;

  Epetra_CrsMatrix* Atmp1 = create_and_fill_crs_matrix(A.RowMap());

  x.PutScalar(1.0);

  Atmp1->Multiply(false, x, rhs);

  soln_rowsum.PutScalar(0.0);

  AztecOO azoo1(Atmp1, &soln_rowsum, &rhs);

  azoo1.SetAztecOption(AZ_scaling, AZ_row_sum);

  azoo1.Iterate(100, 1.e-8);

  delete Atmp1;

  Epetra_CrsMatrix* Atmp2 = create_and_fill_crs_matrix(A.RowMap());

  x.PutScalar(1.0);

  Atmp2->Multiply(false, x, rhs);

  soln_symdiag.PutScalar(0.0);

  AztecOO azoo2(Atmp2, &soln_symdiag, &rhs);

  azoo2.SetAztecOption(AZ_scaling, AZ_sym_diag);

  azoo2.Iterate(100, 1.e-8);

  delete Atmp2;

  Epetra_CrsMatrix* Atmp3 = create_and_fill_crs_matrix(A.RowMap());

  x.PutScalar(1.0);

  Atmp3->Multiply(false, x, rhs);

  soln_none.PutScalar(0.0);

  AztecOO azoo3(Atmp3, &soln_none, &rhs);

  azoo3.SetAztecOption(AZ_scaling, AZ_none);

  azoo3.Iterate(100, 1.e-8);

  delete Atmp3;


  Epetra_CrsMatrix* Atmp4 = create_and_fill_crs_matrix(A.RowMap());

  x.PutScalar(1.0);

  Atmp4->Multiply(false, x, rhs);

  soln_jacobi.PutScalar(0.0);

  AztecOO azoo4(Atmp4, &soln_jacobi, &rhs);

  azoo4.SetAztecOption(AZ_scaling, AZ_Jacobi);

  azoo4.Iterate(100, 1.e-8);

  delete Atmp4;

  //at this point, soln_none, soln_jacobi, soln_rowsum and soln_symdiag
  //should be the same or at least close to the same, since the
  //matrix used in the solution has well-behaved coefficients.
  
  //form vec1 = soln_none - soln_rowsum
  vec1.PutScalar(0.0);
  vec1.Update(1.0, soln_none, 0.0);
  vec1.Update(-1.0, soln_rowsum, 1.0);

  double norm_check1= 0.0;
  vec1.Norm2(&norm_check1);

  //form vec2 = soln_none - soln_symdiag
  vec2.PutScalar(0.0);
  vec2.Update(1.0, soln_none, 0.0);
  vec2.Update(-1.0, soln_symdiag, 1.0);

  double norm_check2= 0.0;
  vec2.Norm2(&norm_check2);

  //form vec3 = soln_none - soln_jacobi
  vec3.PutScalar(0.0);
  vec3.Update(1.0, soln_none, 0.0);
  vec3.Update(-1.0, soln_jacobi, 1.0);

  double norm_check3= 0.0;
  vec3.Norm2(&norm_check3);


  if (std::abs(norm_check1) > 1.e-6) {
    if (verbose) {
      cerr << "AZ_row_sum scaling produced bad soln"
      << endl;
    }
    return(-1);
  }

  if (std::abs(norm_check2) > 1.e-6) {
    if (verbose) {
      cerr << "AZ_sym_diag scaling produced bad soln"
      << endl;
    }
    return(-1);
  }

  if (std::abs(norm_check3) > 1.e-6) {
    if (verbose) {
      cerr << "AZ_Jacobi scaling produced bad soln"
      << endl;
    }
    return(-1);
  }

  options[AZ_pre_calc] = AZ_reuse;

  err = AztecOO_scale_epetra(AZ_SCALE_MAT_RHS_SOL, Amat,
                             options, bvals, xvals, NULL, scaling);
  if (err == 0) {
    if (verbose) {
      cerr << "AztecOO_scale_epetra failed to return err when"
        << " asked to reuse non-existent scaling data."<<endl;
    }
    return(-1);
  }

  options[AZ_keep_info] = 1;
  options[AZ_pre_calc] = AZ_calc;
  err = AztecOO_scale_epetra(AZ_SCALE_MAT_RHS_SOL, Amat,
                             options, bvals, xvals, NULL, scaling);
  if (err != 0) {
    if (verbose) {
      cerr << "AztecOO_scale_epetra returned err=="<<err<<endl;
    }
    return(err);
  }

  options[AZ_keep_info] = 0;
  options[AZ_pre_calc] = AZ_reuse;
  err = AztecOO_scale_epetra(AZ_SCALE_MAT_RHS_SOL, Amat,
                             options, bvals, xvals, NULL, scaling);
  if (err != 0) {
    if (verbose) {
      cerr << "AztecOO_scale_epetra returned err=="<<err
          <<" when asked to reuse scaling data"<<endl;
    }
    return(err);
  }

  options[AZ_pre_calc] = AZ_calc;
  err = AztecOO_scale_epetra(AZ_DESTROY_SCALING_DATA, Amat,
                             options, bvals, xvals, NULL, scaling);
  if (err != 0) {
    if (verbose) {
      std::cerr << "AztecOO_scale_epetra returned err=="<<err
      << " when asked to destroy scaling data."<<std::endl;
    }
    return(err);
  }

  AZ_matrix_destroy(&Amat);
  delete [] options; 
  delete [] params;
  AZ_scaling_destroy(&scaling);

  AZ_manage_memory(0, AZ_CLEAR_ALL, 0, 0, 0);

  return(0);
}

int test_AZ_iterate_AZ_pre_calc_AZ_reuse(Epetra_Comm& Comm,
                                         int* options,
                                         bool verbose)
{
  (void)Comm;
  if (verbose) {
    cout << "testing AZ_keep_info/AZ_reuse with 'old' Aztec (solver "
         <<options[AZ_solver] <<", precond "<<options[AZ_precond]<<"/"
         << options[AZ_subdomain_solve]<<")"<<endl;
  }

  int* proc_config = new int[AZ_PROC_SIZE];

#ifdef EPETRA_MPI
  AZ_set_proc_config(proc_config, MPI_COMM_WORLD);
  AZ_set_comm(proc_config, MPI_COMM_WORLD);
#else
  AZ_set_proc_config(proc_config, 0);
#endif

  //We're going to create 2 Aztec matrices, one MSR and one VBR. We're going
  //to do 2 solves with each, reusing the preconditioner for the second solve.

  int *external, *update_index, *external_index;
  int *external2, *update_index2, *external_index2;

  int i, N = 5;
  AZ_MATRIX* Amsr = NULL;
  AZ_MATRIX* Avbr = NULL;
  int err = create_and_transform_simple_matrix(AZ_MSR_MATRIX, N, 3.0,
                                               proc_config, Amsr,
                                               external, update_index,
                                               external_index);

  err += create_and_transform_simple_matrix(AZ_VBR_MATRIX, N, 3.0,
                                            proc_config, Avbr,
                                            external2, update_index2,
                                            external_index2);


  int* az_options = new int[AZ_OPTIONS_SIZE];
  double* params = new double[AZ_PARAMS_SIZE];
  double* status = new double[AZ_STATUS_SIZE];
  AZ_defaults(az_options, params);
  az_options[AZ_solver] = options[AZ_solver];
  az_options[AZ_precond] = options[AZ_precond];
  az_options[AZ_subdomain_solve] = options[AZ_subdomain_solve];
  az_options[AZ_scaling] = AZ_sym_diag;
  if (verbose) {
    az_options[AZ_output] = AZ_warnings;
  }
  else {
    az_options[AZ_output] = 0;
  }

  int N_update = N+Amsr->data_org[AZ_N_border];
  double* x = new double[N_update];
  double* b = new double[N_update];

  for(i=0; i<N_update; ++i) {
    x[i] = 0.0;
    b[i] = 1.0;
  }

  AZ_PRECOND* Pmsr = AZ_precond_create(Amsr, AZ_precondition, NULL);
  AZ_SCALING* Smsr = AZ_scaling_create();
  AZ_PRECOND* Pvbr = AZ_precond_create(Avbr, AZ_precondition, NULL);
  AZ_SCALING* Svbr = AZ_scaling_create();

 // Amsr->data_org[AZ_name] = 1;
 // Avbr->data_org[AZ_name] = 2;

  //First solve with the first matrix (Amsr).
  if (verbose)
    cout << "solve Amsr, name: "<<Amsr->data_org[AZ_name]<<endl;

  call_AZ_iterate(Amsr, Pmsr, Smsr, x, b, az_options, params, status,
                  proc_config, 1, AZ_calc, verbose);

  //First solve with the second matrix (Avbr).
  if (verbose)
    cout << "solve Avbr, name: " <<Avbr->data_org[AZ_name]<<endl;

  call_AZ_iterate(Avbr, Pvbr, Svbr, x, b, az_options, params, status,
                  proc_config, 0, AZ_calc, verbose);

  //Second solve with Amsr, reusing preconditioner
  if (verbose)
    cout << "solve Amsr (first reuse)"<<endl;

  call_AZ_iterate(Amsr, Pmsr, Smsr, x, b, az_options, params, status,
                  proc_config, 1, AZ_reuse, verbose);

  //Second solve with Avbr, not reusing preconditioner
  if (verbose)
    cout << "solve Avbr (keepinfo==0), name: " <<Avbr->data_org[AZ_name]<<endl;

  call_AZ_iterate(Avbr, Pvbr, Svbr, x, b, az_options, params, status,
                  proc_config, 0, AZ_calc, verbose);

  if (verbose)
    std::cout << "calling AZ_free_memory..."<<std::endl;

  AZ_free_memory(Amsr->data_org[AZ_name]);
  AZ_free_memory(Avbr->data_org[AZ_name]);

  //solve with Amsr again, not reusing preconditioner
  if (verbose)
    cout << "solve Amsr (keepinfo==0)"<<endl;

  call_AZ_iterate(Amsr, Pmsr, Smsr, x, b, az_options, params, status,
                  proc_config, 0, AZ_calc, verbose);

  //Second solve with Avbr, this time with keepinfo==1
  if (verbose)
    cout << "solve Avbr (keepinfo==1), name: " <<Avbr->data_org[AZ_name]<<endl;

  call_AZ_iterate(Avbr, Pvbr, Svbr, x, b, az_options, params, status,
                  proc_config, 1, AZ_calc, verbose);

  //Second solve with Amsr, not reusing preconditioner
  if (verbose)
    cout << "solve Amsr (keepinfo==0, calc)"<<endl;

  call_AZ_iterate(Amsr, Pmsr, Smsr, x, b, az_options, params, status,
                  proc_config, 0, AZ_calc, verbose);

  //Second solve with Avbr, not reusing preconditioner
  if (verbose)
    cout << "solve Avbr (keepinfo==1, reuse), name: "<<Avbr->data_org[AZ_name]<<endl;

  call_AZ_iterate(Avbr, Pvbr, Svbr, x, b, az_options, params, status,
                  proc_config, 1, AZ_reuse, verbose);

  AZ_free_memory(Amsr->data_org[AZ_name]);
  AZ_free_memory(Avbr->data_org[AZ_name]);

  AZ_scaling_destroy(&Smsr);
  AZ_precond_destroy(&Pmsr);
  AZ_scaling_destroy(&Svbr);
  AZ_precond_destroy(&Pvbr);
  destroy_matrix(Amsr);
  destroy_matrix(Avbr);

  delete [] x;
  delete [] b;

  delete [] az_options;
  delete [] params;
  delete [] status;
  delete [] proc_config;
  free(update_index);
  free(external);
  free(external_index);
  free(update_index2);
  free(external2);
  free(external_index2);

  return(0);
}

int call_AZ_iterate(AZ_MATRIX* Amat,
                    AZ_PRECOND* P,
                    AZ_SCALING* S,
                    double* x,
                    double* b,
                    int* options,
                    double* params,
                    double* status,
                    int* proc_config,
                    int keep_info,
                    int pre_calc,
                    bool verbose)
{
  options[AZ_keep_info] = keep_info;
  options[AZ_pre_calc] = pre_calc;

  std::string keepstr;
  std::string calcstr;

  if (keep_info == 1) keepstr = "true";
  else keepstr = "false";

  if (pre_calc == AZ_calc) calcstr = "AZ_calc";
  else calcstr = "AZ_reuse";

  if (verbose)
    cout << "   solve with AZ_keep_info="<<keepstr<<", AZ_pre_calc="<<calcstr<<endl;

  for(int i=0; i<Amat->N_update; ++i) x[i] = 0.0;

  AZ_iterate(x, b, options, params, status, proc_config,
             Amat, P, S);

  return(0);
}
                    
int test_AZ_iterate_then_AZ_scale_f(Epetra_Comm& Comm, bool verbose)
{
  (void)Comm;
  if (verbose) {
    cout << "testing AZ_iterate/AZ_scale_f with 'old' Aztec"<<endl;
  }

  int* proc_config = new int[AZ_PROC_SIZE];

#ifdef EPETRA_MPI
  AZ_set_proc_config(proc_config, MPI_COMM_WORLD);
  AZ_set_comm(proc_config, MPI_COMM_WORLD);
#else
  AZ_set_proc_config(proc_config, 0);
#endif

  int *external, *update_index, *external_index;

  int i, N = 5;
  AZ_MATRIX* Amat = NULL;
  int err = create_and_transform_simple_matrix(AZ_MSR_MATRIX, N, 3.0,
                                               proc_config, Amat,
                                                 external, update_index,
                                                 external_index);
  if (err != 0) {
    return(err);
  }
 
  int* options = new int[AZ_OPTIONS_SIZE];
  double* params = new double[AZ_PARAMS_SIZE];
  double* status = new double[AZ_STATUS_SIZE];
  AZ_defaults(options, params);
  options[AZ_scaling] = AZ_sym_diag;
  if (verbose) {
    options[AZ_output] = AZ_warnings;
  }
  else {
    options[AZ_output] = 0;
  }
  
  int N_update = N+Amat->data_org[AZ_N_border];
  double* x = new double[N_update];
  double* b = new double[N_update];

  for(i=0; i<N_update; ++i) {
    x[i] = 0.0;
    b[i] = 1.0;
  }

  AZ_PRECOND* Pmat = AZ_precond_create(Amat, AZ_precondition, NULL);
  AZ_SCALING* Scal = AZ_scaling_create();

  options[AZ_keep_info] = 1;

  AZ_iterate(x, b, options, params, status, proc_config,
             Amat, Pmat, Scal);

  //now set options[AZ_pre_calc] = AZ_reuse and try to call AZ_scale_f.
  options[AZ_pre_calc] = AZ_reuse;

  AZ_scale_f(AZ_SCALE_MAT_RHS_SOL, Amat, options, b, x, proc_config, Scal);

  AZ_scaling_destroy(&Scal);
  AZ_precond_destroy(&Pmat);
  destroy_matrix(Amat);

  delete [] x;
  delete [] b;

  delete [] options;
  delete [] params;
  delete [] status;
  delete [] proc_config;
  free(update_index);
  free(external);
  free(external_index);

  return(0);
}

void destroy_matrix(AZ_MATRIX*& Amat)
{
  delete [] Amat->update;
  delete [] Amat->val;
  delete [] Amat->bindx;

  if (Amat->matrix_type==AZ_VBR_MATRIX) {
    delete [] Amat->indx;
    delete [] Amat->rpntr;
    free(Amat->cpntr);
    delete [] Amat->bpntr;
  }

  AZ_matrix_destroy(&Amat);
  Amat = NULL;
}

int create_and_transform_simple_matrix(int matrix_type,
                                    int N,
                                    double diag_term,
                                    int* proc_config,
                                    AZ_MATRIX*& Amat,
                                    int*& external,
                                    int*& update_index,
                                    int*& external_index)
{
  //We're going to create a very simple tri-diagonal matrix with diag_term
  //on the diagonal, and -1.0 on the off-diagonals.

  Amat = AZ_matrix_create(N);

  int* update = new int[N];
  int i;
  int numprocs = proc_config[AZ_N_procs];
  int first_eqn = proc_config[AZ_node]*N;
  int adjustment = 0;
  for(i=0; i<N; ++i) {
    update[i] = first_eqn+i;
    if (update[i] == 0 || update[i] == numprocs*N-1) ++adjustment;
  }

  int* data_org;

  //global row 0 and global-N-1 (N*numprocs-1) will have 2 nonzeros in the first
  //and last rows, and there will be 3 nonzeros in all other rows.
  //If you are brave enough to attempt to modify any of the following code,
  //bear in mind that the number of nonzeros per row (3) is hard-coded in
  //a few places.

  int nnz = 3*N - adjustment;
  double* val = new double[nnz+2];
  int* bindx = new int[nnz+2];
  int* indx = NULL;
  int* rpntr = NULL;
  int* cpntr = NULL;
  int* bpntr = NULL;

  int offs = N+1;
  for(i=0; i<N; ++i) {
    val[i] = diag_term;
    bindx[i] = offs;
    int num_off_diagonals = 2;
    if (update[i]==0 || update[i]==numprocs*N-1) num_off_diagonals = 1;
    offs += num_off_diagonals;
  }
  bindx[N] = offs;

  for(i=0; i<N; ++i) {
    int global_row = update[i];

    if (global_row > 0) {
      int ks = bindx[i];
      val[ks] = -1.0;
      bindx[ks] = global_row-1;
    }
    if (global_row < numprocs*N-1) {
      int ke = bindx[i+1]-1;
      val[ke] = -1.0;
      bindx[ke] = global_row+1;
    }
  }

  if (matrix_type == AZ_VBR_MATRIX) {
    //AZ_transform allocates cpntr
    rpntr = new int[N+1];
    bpntr = new int[N+1];
    indx  = new int[nnz+2];

    offs = 0;
    for(i=0; i<N; ++i) {
      rpntr[i] = i;
 
      bpntr[i] = offs;

      if (update[i]==0 || update[i]==numprocs*N-1) offs += 2;
      else offs += 3;
    }
    rpntr[N] = N;
    bpntr[N] = offs;

    for(i=0; i<N; ++i) {
      int global_col = update[i] - 1;
      if (update[i]==0) ++global_col;

      for(int j=bpntr[i]; j<=bpntr[i+1]-1; ++j) {
        if (global_col == update[i]) val[j] = diag_term;
        else val[j] = -1.0;

        bindx[j] = global_col++;
      }
    }

    for(i=0; i<nnz+2; ++i) {
      indx[i] = i;
    }
  }

  AZ_transform(proc_config, &external, bindx, val, update, &update_index,
               &external_index, &data_org, N, indx, bpntr, rpntr,
               &cpntr, matrix_type);

  if (matrix_type == AZ_MSR_MATRIX) {
    AZ_set_MSR(Amat, bindx, val, data_org, N, update, AZ_LOCAL);
  }
  else {
    AZ_set_VBR(Amat, rpntr, cpntr, bpntr, indx, bindx, val, data_org,
               N, update, AZ_LOCAL);
  }

  Amat->must_free_data_org = 1;

  return(0);
}

int test_bug2554(Epetra_Comm& Comm, bool verbose)
{
//This function contains code submitted by Joe Young to
//expose bug 2554. The bug has now been fixed, so this
//function executes without problem. It will be kept as
//a regression test.

  // Construct maps that do not have consecutive indices 
  int             RowIndices[3];
  if (Comm.MyPID() == 0) {
    RowIndices[0] = 1;
    RowIndices[1] = 2;
    RowIndices[2] = 3;

  } else {
    RowIndices[0] = 4;
    RowIndices[1] = 5;
    RowIndices[2] = 6;
  }
  Epetra_Map      RowMap((long long)-1, 3, RowIndices, 0, Comm);

  // Construct a graph with two entries per line 
  Epetra_CrsGraph Graph(Copy, RowMap, 2);
  for (int i = 0; i < RowMap.NumMyElements(); i++) {
    int             ig = RowIndices[i];
    Graph.InsertGlobalIndices(ig, 1, &ig);
  }
  Graph.FillComplete();

  // Make some matrix out of this
  Epetra_FECrsMatrix *Matrix=new Epetra_FECrsMatrix(Copy, Graph);

  // Fill it up with ones 
  Matrix->PutScalar(1.0);

  // Create a rhs and lhs
  Epetra_Vector *rhs=new Epetra_Vector(RowMap);
  Epetra_Vector *lhs=new Epetra_Vector(RowMap);
  rhs->PutScalar(2.0);
  lhs->PutScalar(0.0);


  // Create a solver and problem;
  AztecOO *solver=new AztecOO();
  Epetra_LinearProblem *problem=new Epetra_LinearProblem();

  // Load the problem into the solver
  problem->SetOperator(Matrix);
  problem->SetRHS(rhs);
  problem->SetLHS(lhs);
  solver->SetProblem(*problem, true);

  // Set some options
  solver->SetAztecOption(AZ_solver,AZ_cg);
  solver->SetAztecOption(AZ_precond,AZ_ls);
  solver->SetAztecOption(AZ_poly_ord,9);

  // Solve the problem
  solver->Iterate(50,1e-12);

  // Delete the matrix, lhs, rhs
  delete Matrix;
  delete lhs;
  delete rhs;

  /* Somehow, C++ reallocates objects in the same location where
  they were previously allocated.  So, we need to trick it.  If we
  don't do this, the error will not appear. */
  int *dummy=new int[1000];
  double *dummy2=new double[1000];

  // Reallocate all of them
  Matrix=new Epetra_FECrsMatrix(Copy, Graph);
  rhs=new Epetra_Vector(RowMap);
  lhs=new Epetra_Vector(RowMap);
  Matrix->PutScalar(1.0);
  rhs->PutScalar(2.0);
  lhs->PutScalar(0.0);

  // Load the problem into the solver
  problem->SetOperator(Matrix);
  problem->SetRHS(rhs);
  problem->SetLHS(lhs);
  //For the following call to SetProblem, we must pass 'true' for
  //the optional bool argument which would otherwise default to false.
  //Passing 'true' forces it to internally reset the preconditioner.
  solver->SetProblem(*problem, true);

  // Solve the problem
  solver->Iterate(50,1e-12);

  // Clean up some memory
  solver->UnsetLHSRHS(); // Make sure this function works
  delete problem;
  delete solver;
  delete [] dummy;
  delete [] dummy2;
  delete Matrix;
  delete lhs;
  delete rhs;

  return(0);
}

int test_bug2890(Epetra_Comm& Comm, bool verbose)
{
//This function tests the AZ_random1() function in AztecOO.
//The implementation of the Park and Miller random number
//generator was incorrect and resulted in an overflow condition.
//This is *not* a complete test of AztecOO's RNG.
//
//A more robust check is to compile AztecOO with gcc -ftrapv and run
//a Krylov method that invokes AZ_random_vector().

  int seed = -127773;
  double rand_num;

  rand_num = AZ_srandom1(&seed);

  if (verbose && Comm.MyPID() == 0)
    printf("test_bug2890: rand_num = %e (should be in [0,1])\n",rand_num);

  if ( (rand_num > 1) || (rand_num < 0) )
    return 1;    // rand_num should be in [0,1]
  else
    return 0;
}
