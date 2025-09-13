// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// This software is a result of the research described in the report
//
//     "A comparison of algorithms for modal analysis in the absence
//     of a sparse direct method", P. Arbenz, R. Lehoucq, and U. Hetmaniuk,
//     Sandia National Laboratories, Technical report SAND2003-1028J.
//
// It is based on the Epetra, AztecOO, and ML packages defined in the Trilinos
// framework ( http://trilinos.org/ ).

#include "AMGOperator.h"


AMGOperator::AMGOperator(const Epetra_Comm &_Comm, const Epetra_Operator *KK, int verb,
                         int nLevel, int smoother, int param,
                         int coarseSolver, int cycle,
                         int _numDofs, const Epetra_MultiVector *Z)
               : MyComm(_Comm),
                 callBLAS(),
                 callLAPACK(),
                 K(KK),
                 Prec(0),
                 Q(Z),
                 QtQ(0),
                 numDofs(_numDofs),
                 leftProjection(false),
                 rightProjection(false),
                 ml_handle(0),
                 ml_agg(0),
                 AMG_NLevels(nLevel),
                 coarseLocalSize(0),
                 coarseGlobalSize(0),
                 ZcoarseTZcoarse(0),
                 verbose(verb)
               {

  preProcess((numDofs == 1) ? 1 : 50);

  // Set a smoother for the MG method
  if (smoother == 1) {
    if (numDofs == 1) {
      for (int j = 0; j < AMG_NLevels; ++j)
        ML_Gen_Smoother_MLS(ml_handle, j, ML_BOTH, 30., param);
    }
    else {
      int nBlocks, *blockIndices;
      for (int j = 0; j < AMG_NLevels; ++j) {
        ML_Gen_Blocks_Aggregates(ml_agg, j, &nBlocks, &blockIndices);
        ML_Gen_Smoother_BlockDiagScaledCheby(ml_handle, j, ML_BOTH, 30., param,
                                             nBlocks, blockIndices);
      }
    } // if (numDofs == 1)
  }

  // Note that Symmetric Gauss Seidel does not parallelize well
  if (smoother == 2) {
    // Note: ML_DEFAULT specifies the damping parameter
    //       param specifies the number of iterations
    if (numDofs == 1)
      ML_Gen_Smoother_SymGaussSeidel(ml_handle, ML_ALL_LEVELS, ML_BOTH, param, ML_DEFAULT);
    else
      ML_Gen_Smoother_BlockGaussSeidel(ml_handle, ML_ALL_LEVELS, ML_BOTH, param, ML_DEFAULT,
                                       numDofs);
  }

  if (smoother == 3) {
    // Note: ML_DEFAULT specifies the damping parameter
    //       param specifies the number of iterations
    if (numDofs == 1)
      ML_Gen_Smoother_Jacobi(ml_handle, ML_ALL_LEVELS, ML_BOTH, param, ML_DEFAULT);
    else {
      int size = ml_handle->Amat[0].getrow->Nrows;
      int *blockIndices = new int[size];
      int j;
      for (j = 0; j < size; ++j)
        blockIndices[j] = j/numDofs;
      for (j = 0; j < AMG_NLevels; ++j) {
        size = ml_handle->Amat[j].getrow->Nrows;
        ML_Gen_Smoother_VBlockJacobi(ml_handle, j, ML_BOTH, param, ML_DEFAULT,
                                     size/numDofs, blockIndices);
      }
      delete[] blockIndices;
    } // if (numDofs == 1)
  }

  setCoarseSolver_Cycle(coarseSolver, cycle);

}


AMGOperator::AMGOperator(const Epetra_Comm &_Comm, const Epetra_Operator *KK, int verb,
                         int nLevel, int smoother, int *param,
                         int coarseSolver, int cycle,
                         int _numDofs, const Epetra_MultiVector *Z)
               : MyComm(_Comm),
                 callBLAS(),
                 callLAPACK(),
                 K(KK),
                 Prec(0),
                 Q(Z),
                 QtQ(0),
                 numDofs(_numDofs),
                 leftProjection(false),
                 rightProjection(false),
                 ml_handle(0),
                 ml_agg(0),
                 AMG_NLevels(nLevel),
                 coarseLocalSize(0),
                 coarseGlobalSize(0),
                 ZcoarseTZcoarse(0),
                 verbose(verb)
               {

  preProcess((numDofs == 1) ? 1 : 50);

  // Set a smoother for the MG method
  if (smoother == 1) {
    if (numDofs == 1) {
      for (int j = 0; j < AMG_NLevels; ++j)
        ML_Gen_Smoother_MLS(ml_handle, j, ML_BOTH, 30., (param) ? param[j] : 2);
    }
    else {
      int nBlocks, *blockIndices;
      for (int j = 0; j < AMG_NLevels; ++j) {
        ML_Gen_Blocks_Aggregates(ml_agg, j, &nBlocks, &blockIndices);
        ML_Gen_Smoother_BlockDiagScaledCheby(ml_handle, j, ML_BOTH, 30., (param) ? param[j] : 2,
                                             nBlocks, blockIndices);
      }
    } // if (numDofs == 1)
  }

  // Note that Symmetric Gauss Seidel does not parallelize well
  if (smoother == 2) {
    // Note: ML_DEFAULT specifies the damping parameter
    //       param specifies the number of iterations
    if (numDofs == 1) {
      for (int j = 0; j < AMG_NLevels; ++j) {
        ML_Gen_Smoother_SymGaussSeidel(ml_handle, j, ML_BOTH, (param) ? param[j] : 1, ML_DEFAULT);
      }
    }
    else {
      for (int j = 0; j < AMG_NLevels; ++j) {
        ML_Gen_Smoother_BlockGaussSeidel(ml_handle, j, ML_BOTH, (param) ? param[j] : 1,
                                         ML_DEFAULT, numDofs);
      }
    }
  }

  if (smoother == 3) {
    // Note: ML_DEFAULT specifies the damping parameter
    //       param specifies the number of iterations
    if (numDofs == 1) {
      for (int j = 0; j < AMG_NLevels; ++j) {
        ML_Gen_Smoother_Jacobi(ml_handle, j, ML_BOTH, (param) ? param[j] : 1, ML_DEFAULT);
      }
    }
    else {
      int size = ml_handle->Amat[0].getrow->Nrows;
      int *blockIndices = new int[size];
      int j;
      for (j = 0; j < size; ++j)
        blockIndices[j] = j/numDofs;
      for (j = 0; j < AMG_NLevels; ++j) {
        size = ml_handle->Amat[j].getrow->Nrows;
        ML_Gen_Smoother_VBlockJacobi(ml_handle, j, ML_BOTH, (param) ? param[j] : 1,
                                     ML_DEFAULT, size/numDofs, blockIndices);
      }
      delete[] blockIndices;
    } // if (numDofs == 1)
  }

  setCoarseSolver_Cycle(coarseSolver, cycle);

}


void AMGOperator::preProcess(int maxCoarseSize) {

  // Note the constness is cast away for ML
  Epetra_RowMatrix *Krow = dynamic_cast<Epetra_RowMatrix*>(const_cast<Epetra_Operator*>(K));
  if (Krow == 0) {
    if (MyComm.MyPID() == 0) {
      cerr << endl;
      cerr << " !!! For AMG preconditioner, the matrix must be 'Epetra_RowMatrix' object !!!\n";
      cerr << endl;
    }
    return;
  }

  // Generate an ML multilevel preconditioner
  ML_Set_PrintLevel(verbose);
  ML_Create(&ml_handle, AMG_NLevels);

  EpetraMatrix2MLMatrix(ml_handle, 0, Krow);

  ML_Aggregate_Create(&ml_agg);
  ML_Aggregate_Set_MaxCoarseSize(ml_agg, maxCoarseSize);
  ML_Aggregate_Set_Threshold(ml_agg, 0.0);

  // Set the aggregation scheme
  if ((Q == 0) || (numDofs == 1)) {
    ML_Aggregate_Set_CoarsenScheme_Uncoupled(ml_agg);
  }
  else {
    //--------------------------------------
    //ML_Aggregate_Set_DampingFactor(ml_agg, 1.3333);
    //ML_Aggregate_Set_CoarsenScheme_METIS(ml_agg);
    //ML_Aggregate_Set_NodesPerAggr(ml_handle, ml_agg, -1, );
    //--------------------------------------
    ML_Aggregate_Set_DampingFactor(ml_agg, 1.3333);
    ML_Aggregate_Set_CoarsenScheme_MIS(ml_agg);
    //ML_Aggregate_Set_Phase3AggregateCreationAggressiveness(ml_agg, 0.1);
    //--------------------------------------
    //ML_Aggregate_Set_CoarsenScheme_Uncoupled(ml_agg);
    //--------------------------------------
    //ML_Aggregate_Set_BlockDiagScaling(ml_agg);
    //--------------------------------------
    ML_Aggregate_Set_NullSpace(ml_agg, numDofs, Q->NumVectors(), Q->Values(), Q->MyLength());
    QtQ = new double[Q->NumVectors()*Q->NumVectors()];
    double *work = new double[Q->NumVectors()*Q->NumVectors()];
    callBLAS.GEMM('T', 'N', Q->NumVectors(), Q->NumVectors(), Q->MyLength(),
                  1.0, Q->Values(), Q->MyLength(), Q->Values(), Q->MyLength(),
                  0.0, work, Q->NumVectors());
    MyComm.SumAll(work, QtQ, Q->NumVectors()*Q->NumVectors());
    delete[] work;
    int info = 0;
    callLAPACK.POTRF('U', Q->NumVectors(), QtQ, Q->NumVectors(), &info);
    if (info) {
      cerr << endl;
      cerr << " !!! In AMG preconditioner, the null space vectors are linearly dependent !!!\n";
      cerr << endl;
      delete[] QtQ;
      QtQ = 0;
    }
  }

  AMG_NLevels = ML_Gen_MGHierarchy_UsingAggregation(ml_handle, 0, ML_INCREASING, ml_agg);

}


void AMGOperator::setCoarseSolver_Cycle(int coarseSolver, int cycle) {

  // Specifies the solver on the coarsest level
  if (coarseSolver > -1) {
    switch (coarseSolver) {
      case 0:
        // Use SuperLU on the coarsest level
        ML_Gen_CoarseSolverSuperLU(ml_handle, AMG_NLevels-1);
        break;
      case 1:
        // Use Aztec PCG on the coarsest level
        int options[AZ_OPTIONS_SIZE];
        double params[AZ_PARAMS_SIZE];
        int *proc_config = new int[AZ_PROC_SIZE];
        double *status = new double[AZ_STATUS_SIZE];
        AZ_defaults(options, params);
        options[AZ_solver] = AZ_cg;
        options[AZ_precond] = AZ_user_precond;
        options[AZ_scaling] = AZ_none;
        options[AZ_output] = 1;
//        if (verbose < 3)
//          options[AZ_output] = AZ_last;
        if (verbose < 2)
          options[AZ_output] = AZ_none;
#ifdef EPETRA_MPI
        AZ_set_proc_config(proc_config, get_global_comm());
#endif
        params[AZ_tol] = 1.0e-14;
        coarseLocalSize = ml_handle->Amat[AMG_NLevels-1].invec_leng;
        MyComm.SumAll(&coarseLocalSize, &coarseGlobalSize, 1);
        if ((verbose > 0) && (MyComm.MyPID() == 0))
          cout << "\n --> Total number of coarse grid unknowns = " << coarseGlobalSize << endl;
        // Define the data for the projection
        if (Q) {
          int zCol = ml_agg->nullspace_dim;
          double *ZZ = ml_agg->nullspace_vect;
          ZcoarseTZcoarse = new double[zCol*zCol];
          double *tmp = new double[zCol*zCol];
          int info = 0;
          callBLAS.GEMM('T', 'N', zCol, zCol, coarseLocalSize, 1.0, ZZ, coarseLocalSize,
                        ZZ, coarseLocalSize, 0.0, tmp, zCol);
          MyComm.SumAll(tmp, ZcoarseTZcoarse, zCol*zCol);
          callLAPACK.POTRF('U', zCol, ZcoarseTZcoarse, zCol, &info);
          if (info) {
            cerr << endl;
            cerr << " !!! In AMG preconditioner, for the coarse level, ";
            cerr << "the null space vectors are linearly dependent !!!\n";
            cerr << endl;
          }
#ifndef MACOSX
          singularCoarse::setNullSpace(ml_agg->nullspace_vect, coarseLocalSize,
                                       ml_agg->nullspace_dim, ZcoarseTZcoarse, &MyComm);
#endif
          delete[] tmp;
        }
        options[AZ_max_iter] = coarseGlobalSize;
//        // Subtract off the rigid body modes
//        options[AZ_orth_kvecs] = AZ_TRUE;
//        options[AZ_keep_kvecs] = (Q) ? coarseGlobalSize - Q->NumVectors() : coarseGlobalSize;
//        options[AZ_apply_kvecs] = AZ_APPLY;
#ifdef MACOSX
        ML_Gen_SmootherAztec(ml_handle, AMG_NLevels-1, options, params, proc_config, status,
                             options[AZ_max_iter], ML_PRESMOOTHER, NULL);
#else
        ML_Gen_SmootherAztec(ml_handle, AMG_NLevels-1, options, params, proc_config, status,
                             options[AZ_max_iter], ML_PRESMOOTHER,
                             (Q) ?  singularCoarse::projection : NULL);
#endif
        break;
    }
  }

  // Type of cycle
  switch (cycle) {
    default:
    case 0:
      ML_Gen_Solver(ml_handle, ML_MGV, 0, AMG_NLevels-1);
      break;
    case 1:
      ML_Gen_Solver(ml_handle, ML_MGW, 0, AMG_NLevels-1);
      break;
  }

  Prec = new ML_Epetra::MultiLevelOperator(ml_handle, MyComm, K->OperatorDomainMap(),
                                           K->OperatorDomainMap());

  //------------------------------------
  // This definition of Prec calls the old class 'Epetra_ML_Operator'
  //Prec = new Epetra_ML_Operator(ml_handle, MyComm, K->OperatorDomainMap(),
  //                              K->OperatorDomainMap());
  //------------------------------------

}


AMGOperator::~AMGOperator() {

  if (Prec) {
    delete Prec;
    Prec = 0;
    ML_Destroy(&ml_handle);
    ML_Aggregate_Destroy(&ml_agg);
  }

  if (QtQ) {
    delete[] QtQ;
    QtQ = 0;
  }

  if (ZcoarseTZcoarse) {
    delete[] ZcoarseTZcoarse;
    ZcoarseTZcoarse = 0;
#ifndef MACOSX
    singularCoarse::setNullSpace(0, 0, 0, 0, 0);
#endif
  }

}


int AMGOperator::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

  return K->Apply(X, Y);

}


int AMGOperator::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

  int xcol = X.NumVectors();

  if (Y.NumVectors() < xcol)
    return -1;

  double *valX = (rightProjection == true) ? new double[xcol*X.MyLength()] : X.Values();

  Epetra_MultiVector X2(View, X.Map(), valX, X.MyLength(), xcol);

  if ((rightProjection == true) && (Q)) {
    int qcol = Q->NumVectors();
    double *work = new double[2*qcol*xcol];
    memcpy(X2.Values(), X.Values(), xcol*X.MyLength()*sizeof(double));
    callBLAS.GEMM('T', 'N', qcol, xcol, Q->MyLength(), 1.0, Q->Values(), Q->MyLength(),
                  X2.Values(), X2.MyLength(), 0.0, work + qcol*xcol, qcol);
    MyComm.SumAll(work + qcol*xcol, work, qcol*xcol);
    int info = 0;
    callLAPACK.POTRS('U', qcol, xcol, QtQ, qcol, work, qcol, &info);
    callBLAS.GEMM('N', 'N', X2.MyLength(), xcol, qcol, -1.0, Q->Values(), Q->MyLength(),
                  work, qcol, 1.0, X2.Values(), X2.MyLength());
    delete[] work;
  }

  Prec->ApplyInverse(X2, Y);

  if (rightProjection == true)
    delete[] valX;

  if ((leftProjection == true) && (Q)) {
    int qcol = Q->NumVectors();
    double *work = new double[2*qcol*xcol];
    callBLAS.GEMM('T', 'N', qcol, xcol, Q->MyLength(), 1.0, Q->Values(), Q->MyLength(),
                  Y.Values(), Y.MyLength(), 0.0, work + qcol*xcol, qcol);
    MyComm.SumAll(work + qcol*xcol, work, qcol*xcol);
    int info = 0;
    callLAPACK.POTRS('U', qcol, xcol, QtQ, qcol, work, qcol, &info);
    callBLAS.GEMM('N', 'N', Y.MyLength(), xcol, qcol, -1.0, Q->Values(), Q->MyLength(),
                  work, qcol, 1.0, Y.Values(), Y.MyLength());
    delete[] work;
  }

  return 0;

}


