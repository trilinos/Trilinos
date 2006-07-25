//@HEADER
// ***********************************************************************
//
//        AztecOO: An Object-Oriented Aztec Linear Solver Package
//                 Copyright (2002) Sandia Corporation
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

  int local_n = 20;
  int global_n = numprocs*local_n;

  Epetra_Map emap(global_n, 0, comm);

  Epetra_CrsMatrix A(Copy, emap, 3);
  Epetra_Vector x(emap), b(emap);

  x.PutScalar(1.0);

  int myFirstGlobalRow = localproc*local_n;
  int globalCols[3];
  double values[3];

  for(int i=0; i<local_n; ++i) {
    int globalRow = myFirstGlobalRow +i;

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

    A.InsertGlobalValues(globalRow, numcols, values, globalCols);
  }

  A.FillComplete();

  A.Multiply(false, x, b);
  x.PutScalar(0.0);

//  AztecOO azoo(&A, &x, &b);
//  azoo.PrintLinearSystem("test");

  double initial_norm = resid2norm(A, x, b);

  if (verbose) {
    cout << "Initial 2-norm of b-A*x: "<<initial_norm<<endl;
  }

//  AZ_manage_memory(0, -43, 0, 0, 0);

  int err = test_azoo_as_precond_op(A, x, b, verbose);
  if (err != 0) {
    cout << "test_azoo_as_precond_op err, test FAILED."<<endl;
    return(err);
  }

//  AZ_manage_memory(0, -43, 0, 0, 0);

  err = test_azoo_with_ilut(A, x, b, verbose);
  if (err != 0) {
    cout << "test_azoo_with_ilut err, test FAILED."<<endl;
    return(err);
  }

  err = test_azoo_scaling(A, x, b, verbose);
  if (err != 0) {
    cout << "test_azoo_scaling err="<<err<<", test FAILED."<<endl;
  }

  int* options = new int[AZ_OPTIONS_SIZE];
  options[AZ_solver] = AZ_cg;
  options[AZ_subdomain_solve] = AZ_none;
  options[AZ_precond] = AZ_Jacobi;

  if (verbose)
    std::cout << "about to call test_AZ_iterate_AZ_pre_calc_AZ_reuse"
       <<std::endl;
//  AZ_manage_memory(0, -43, 0, 0, 0);

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

  delete azoo0;

  return(0);
}

int test_azoo_scaling(Epetra_CrsMatrix& A,
                      Epetra_Vector& x,
                      Epetra_Vector& b,
                      bool verbose)
{
  Epetra_CrsMatrix Atmp(A);
  Epetra_Vector vec1(x);
  Epetra_Vector vec2(x);
  Epetra_Vector diag(x);
  Epetra_Vector vec3(x);
  Epetra_Vector vec4(x);
  Epetra_Vector rhs(x);
  Epetra_Vector soln_rowsum(x);
  Epetra_Vector soln_symdiag(x);

  vec1.PutScalar(1.0);

  Atmp.Multiply(false, vec1, vec2);

  Atmp.ExtractDiagonalCopy(diag);

  double* diag_vals = NULL;
  diag.ExtractView(&diag_vals);

  int* options = new int[AZ_OPTIONS_SIZE];
  double* params = new double[AZ_PARAMS_SIZE];
  AZ_defaults(options, params);

  options[AZ_output] = verbose ? 1 : AZ_none;

  options[AZ_scaling] = AZ_Jacobi;
  AztecOO::MatrixData mdata(&Atmp);
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

  Atmp.Multiply(false, vec1, vec3);

  vec4.Multiply(1.0, diag, vec3, 0.0);

  double vec2nrm, vec4nrm;

  vec2.Norm2(&vec2nrm);
  vec4.Norm2(&vec4nrm);

  if (vec2nrm != vec4nrm) {
    return(-1);
  }

  x.PutScalar(1.0);

  Atmp.Multiply(false, x, rhs);

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

  Epetra_CrsMatrix Atmp1(A);

  x.PutScalar(1.0);

  Atmp1.Multiply(false, x, rhs);

  soln_rowsum.PutScalar(0.0);

  AztecOO azoo1(&Atmp1, &soln_rowsum, &rhs);

  azoo1.SetAztecOption(AZ_scaling, AZ_row_sum);

  azoo1.Iterate(100, 1.e-8);

  Epetra_CrsMatrix Atmp2(A);

  x.PutScalar(1.0);

  Atmp2.Multiply(false, x, rhs);

  soln_symdiag.PutScalar(0.0);

  AztecOO azoo2(&Atmp2, &soln_symdiag, &rhs);

  azoo2.SetAztecOption(AZ_scaling, AZ_sym_diag);

  azoo2.Iterate(100, 1.e-8);

  //at this point, soln_rowsum and soln_symdiag should be the same
  //or at least very close to the same.
  
  //form vec1 = soln_symdiag - soln_rowsum
  vec1.PutScalar(0.0);
  vec1.Update(1.0, soln_symdiag, 0.0);
  vec1.Update(-1.0, soln_rowsum, 1.0);

  double norm_check = 0.0;
  vec1.Norm2(&norm_check);

  if (std::abs(norm_check) > 1.e-6) {
    if (verbose) {
      cerr << "AZ_row_sum and AZ_sym_diag produced different solns"
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
  int err = create_and_transform_simple_matrix(AZ_MSR_MATRIX, N, 2.0,
                                               proc_config, Amsr,
                                               external, update_index,
                                               external_index);

  err += create_and_transform_simple_matrix(AZ_VBR_MATRIX, N, 2.0,
                                            proc_config, Avbr,
                                            external2, update_index2,
                                            external_index2);

//  Epetra_MsrMatrix emsr(proc_config, Amsr);

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

  double* x = new double[N];
  double* b = new double[N];

  for(i=0; i<N; ++i) {
    x[i] = 0.0;
    b[i] = 1.0;
  }

  AZ_PRECOND* Pmsr = AZ_precond_create(Amsr, AZ_precondition, NULL);
  AZ_SCALING* Smsr = AZ_scaling_create();
  AZ_PRECOND* Pvbr = AZ_precond_create(Avbr, AZ_precondition, NULL);
  AZ_SCALING* Svbr = AZ_scaling_create();

 // Amsr->data_org[AZ_name] = 1;
 // Avbr->data_org[AZ_name] = 2;

//  AZ_manage_memory(0, -43, 0, 0, 0);

  //First solve with the first matrix (Amsr).
  if (verbose)
    cout << "solve Amsr, name: "<<Amsr->data_org[AZ_name]<<endl;

  call_AZ_iterate(Amsr, Pmsr, Smsr, x, b, az_options, params, status,
                  proc_config, 1, AZ_calc, verbose);

//  AZ_manage_memory(0, -43, 0, 0, 0);

  //First solve with the second matrix (Avbr).
  if (verbose)
    cout << "solve Avbr, name: " <<Avbr->data_org[AZ_name]<<endl;

//  AZ_manage_memory(0, -43, 0, 0, 0);

  call_AZ_iterate(Avbr, Pvbr, Svbr, x, b, az_options, params, status,
                  proc_config, 0, AZ_calc, verbose);

//  AZ_manage_memory(0, -43, 0, 0, 0);

  //Second solve with Amsr, reusing preconditioner
  if (verbose)
    cout << "solve Amsr (first reuse)"<<endl;

  call_AZ_iterate(Amsr, Pmsr, Smsr, x, b, az_options, params, status,
                  proc_config, 1, AZ_reuse, verbose);

//  AZ_manage_memory(0, -43, 0, 0, 0);

  //Second solve with Avbr, not reusing preconditioner
  if (verbose)
    cout << "solve Avbr (keepinfo==0), name: " <<Avbr->data_org[AZ_name]<<endl;

  call_AZ_iterate(Avbr, Pvbr, Svbr, x, b, az_options, params, status,
                  proc_config, 0, AZ_calc, verbose);

//  AZ_manage_memory(0, -43, 0, 0, 0);
  if (verbose)
    std::cout << "calling AZ_free_memory..."<<std::endl;

  AZ_free_memory(Amsr->data_org[AZ_name]);
  AZ_free_memory(Avbr->data_org[AZ_name]);

//  AZ_manage_memory(0, -43, 0, 0, 0);

  //solve with Amsr again, not reusing preconditioner
  if (verbose)
    cout << "solve Amsr (keepinfo==0)"<<endl;

  call_AZ_iterate(Amsr, Pmsr, Smsr, x, b, az_options, params, status,
                  proc_config, 0, AZ_calc, verbose);

//  AZ_manage_memory(0, -43, 0, 0, 0);

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
  int err = create_and_transform_simple_matrix(AZ_MSR_MATRIX, N, 2.0,
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
  
  double* x = new double[N];
  double* b = new double[N];

  for(i=0; i<N; ++i) {
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
  //We're going to create a very simple diagonal matrix.

  Amat = AZ_matrix_create(N);

  int* update = new int[N];
  int i;
  int first_eqn = proc_config[AZ_node]*N;
  for(i=0; i<N; ++i) {
    update[i] = first_eqn+i;
  }

  int* data_org;

  int nnz = N;
  double* val = new double[nnz+1];
  int* bindx = new int[nnz+1];
  int* indx = NULL;
  int* rpntr = NULL;
  int* cpntr = NULL;
  int* bpntr = NULL;

  for(i=0; i<nnz+1; ++i) {
    val[i] = diag_term;
    bindx[i] = N+1;
  }

  if (matrix_type == AZ_VBR_MATRIX) {
    //AZ_transform allocates cpntr
    rpntr = new int[N+1];
    bpntr = new int[N+1];
    indx  = new int[N+1];

    for(i=0; i<N+1; ++i) {
      rpntr[i] = i;
      bpntr[i] = i;
      bindx[i] = first_eqn+i;
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
