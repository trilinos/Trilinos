// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//--------------------------------------------------------------------
// This file is a self-contained example of rebalancing the contents of
// a matrix and a vector among different MPI nodes. The matrix is formed,
// on purpose, in an unbalanced way.
//
// Originally, this example used EpetraExt + Isorropia for rebalance. It
// was a modified version of the matrix_1.cpp example in Isorropia.
//
// This example now uses Tpetra + Zoltan2, combining some aspects of the
// original example using Isorropia with aspects from example/graph/graph.cpp
// in Zoltan2.
//--------------------------------------------------------------------

#include <Tpetra_Rebalance_LinearProblem_def.hpp>

bool runExample(const Teuchos::RCP<const Teuchos::Comm<int> >& comm) {
  int numProcs( comm->getSize() );
  int localProc( comm->getRank() );

  using Map_t      = Tpetra::Map<>;
  using LocalId_t  = Map_t::local_ordinal_type;
  using GlobalId_t = Map_t::global_ordinal_type;
  using Node_t     = Map_t::node_type;
  using Scalar_t   = Tpetra::Details::DefaultTypes::scalar_type;
  using MultiV_t   = Tpetra::MultiVector<Scalar_t, LocalId_t, GlobalId_t, Node_t>;
  using Matrix_t   = Tpetra::CrsMatrix<Scalar_t, LocalId_t, GlobalId_t, Node_t>;
  using Problem_t  = Tpetra::LinearProblem<Scalar_t, LocalId_t, GlobalId_t, Node_t>;

  // ****************************************************************
  // Step 1/10: create a Tpetra::CrsMatrix with rows spread un-evenly
  //            over processors.
  // ****************************************************************
  LocalId_t  local_num_rows  = 200;
  LocalId_t  nnz_per_row     = local_num_rows/4+1;
  GlobalId_t global_num_rows = numProcs*local_num_rows;

  int mid_proc = numProcs/2;
  bool num_procs_even = numProcs%2==0 ? true : false;
  LocalId_t adjustment = local_num_rows/2;

  // Adjust local_num_rows so that it's not equal on all procs.
  if (localProc < mid_proc) {
    local_num_rows -= adjustment;
  }
  else {
    local_num_rows += adjustment;
  }

  // If numProcs is not an even number, undo the local_num_rows adjustment
  // on one proc so that the total will still be correct.
  if (localProc == numProcs-1) {
    if (num_procs_even == false) {
      local_num_rows -= adjustment;
    }
  }

  // Now we're ready to create a row-map.
  Teuchos::RCP< Map_t > rowmap = Teuchos::rcp( new Map_t(global_num_rows, local_num_rows, 0, comm) );
  Teuchos::RCP< Matrix_t > originalMatrix = Teuchos::rcp( new Matrix_t(rowmap, nnz_per_row) );

  std::vector<GlobalId_t> indices(nnz_per_row);
  std::vector<double> coefs(nnz_per_row);

  for (LocalId_t i(0); i < local_num_rows; ++i) {
    GlobalId_t global_row = rowmap->getGlobalElement(i);
    GlobalId_t first_col = global_row - nnz_per_row/2;

    if (first_col < 0) {
      first_col = 0;
    }
    else if (first_col > (global_num_rows - nnz_per_row)) {
      first_col = global_num_rows - nnz_per_row;
    }

    for (LocalId_t j(0); j < nnz_per_row; ++j) {
      indices[j] = first_col + j;
      coefs[j] = 1.0;
    }

    originalMatrix->insertGlobalValues(global_row,
                                       nnz_per_row,
                                       &coefs[0],
                                       &indices[0]);
  }

  originalMatrix->fillComplete();

  // ****************************************************************
  // Step 2/10: create the original linear problem
  // ****************************************************************
  size_t numVectors(2);

  Teuchos::RCP<MultiV_t> originalLhs( Teuchos::null );
  originalLhs  = Tpetra::createMultiVector<Scalar_t,LocalId_t,GlobalId_t>( originalMatrix->getDomainMap(), numVectors );
  originalLhs->randomize();

  Teuchos::RCP<MultiV_t> originalRhs( Teuchos::null );
  originalRhs  = Tpetra::createMultiVector<Scalar_t,LocalId_t,GlobalId_t>( originalMatrix->getRangeMap(), numVectors );
  originalMatrix->apply(*originalLhs, *originalRhs);

  Teuchos::RCP<Problem_t> originalLP = Teuchos::rcp<Problem_t>( new Problem_t(originalMatrix, originalLhs, originalRhs) );

  // ****************************************************************
  // Step 3/10: create the rebalance transform
  // ****************************************************************
  std::string partitioningMethod("block"); // "scotch" "parmetis"
  Teuchos::RCP<Teuchos::ParameterList> paramListForZoltan2PartitioningProblem = Teuchos::rcp<Teuchos::ParameterList>( new Teuchos::ParameterList() );
  paramListForZoltan2PartitioningProblem->set("partitioning_approach", "partition");
  paramListForZoltan2PartitioningProblem->set("algorithm", partitioningMethod);

  Tpetra::Rebalance_LinearProblem<Scalar_t, LocalId_t, GlobalId_t, Node_t> rebalanceTransform( paramListForZoltan2PartitioningProblem );

  // ****************************************************************
  // Step 4/10: rebalance the original linear problem (its matrix and vectors)
  // ****************************************************************
  Teuchos::RCP<Problem_t> transformedLP = rebalanceTransform( originalLP );
  rebalanceTransform.fwd();

  Teuchos::RCP<Matrix_t> rebalancedMatrix = Teuchos::rcp<Matrix_t>( dynamic_cast<Matrix_t *>(transformedLP->getMatrix().get()), false );
  Teuchos::RCP<MultiV_t> rebalancedLhs    = Teuchos::rcp<MultiV_t>( transformedLP->getLHS().get(), false );
  Teuchos::RCP<MultiV_t> rebalancedRhs    = Teuchos::rcp<MultiV_t>( transformedLP->getRHS().get(), false );

  // ****************************************************************
  // Step 5/10: compare the matrix balancing among MPI nodes, before
  //            and after the redistribution.
  // ****************************************************************
  bool allOk1(true);
  {
    LocalId_t originalLocalNumRows   = originalMatrix->getLocalNumRows();
    LocalId_t rebalancedLocalNumRows = rebalancedMatrix->getLocalNumRows();
    LocalId_t overallMin_originalLocalNumRows(0);
    LocalId_t overallMax_originalLocalNumRows(0);
    LocalId_t overallMin_rebalancedLocalNumRows(0);
    LocalId_t overallMax_rebalancedLocalNumRows(0);
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, 1, &originalLocalNumRows,   &overallMin_originalLocalNumRows  );
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, 1, &originalLocalNumRows,   &overallMax_originalLocalNumRows  );
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, 1, &rebalancedLocalNumRows, &overallMin_rebalancedLocalNumRows);
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, 1, &rebalancedLocalNumRows, &overallMax_rebalancedLocalNumRows);

    LocalId_t originalLocalNNZ   = originalMatrix->getLocalNumEntries();
    LocalId_t rebalancedLocalNNZ = rebalancedMatrix->getLocalNumEntries();
    LocalId_t overallMin_originalLocalNNZ(0);
    LocalId_t overallMax_originalLocalNNZ(0);
    LocalId_t overallMin_rebalancedLocalNNZ(0);
    LocalId_t overallMax_rebalancedLocalNNZ(0);
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, 1, &originalLocalNNZ,   &overallMin_originalLocalNNZ  );
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, 1, &originalLocalNNZ,   &overallMax_originalLocalNNZ  );
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, 1, &rebalancedLocalNNZ, &overallMin_rebalancedLocalNNZ);
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, 1, &rebalancedLocalNNZ, &overallMax_rebalancedLocalNNZ);

    comm->barrier();
    for (int p(0); p < numProcs; ++p) {
      comm->barrier();

      if (p != localProc) continue;

      std::cout << "proc " << p
                << ": original matrix local num rows = " << originalLocalNumRows
                << "( min = " << overallMin_originalLocalNumRows
                << ", max = " << overallMax_originalLocalNumRows
                << "); original local NNZ = " << originalLocalNNZ
                << "( min = " << overallMin_originalLocalNNZ
                << ", max = " << overallMax_originalLocalNNZ
                << ")"
                << std::endl;
      std::cout.flush();
      std::cout << "proc " << p
                << ": rebalanced matrix local num rows = " << rebalancedLocalNumRows
                << "( min = " << overallMin_rebalancedLocalNumRows
                << ", max = " << overallMax_rebalancedLocalNumRows
                << "); rebalanced local NNZ = " << rebalancedLocalNNZ
                << "( min = " << overallMin_rebalancedLocalNNZ
                << ", max = " << overallMax_rebalancedLocalNNZ
                << ")"
                << std::endl;
      std::cout.flush();
    }
    comm->barrier();

    if ((overallMin_originalLocalNumRows   <= overallMin_rebalancedLocalNumRows) &&
        (overallMax_rebalancedLocalNumRows <= overallMax_originalLocalNumRows  ) &&
        (overallMin_originalLocalNNZ       <= overallMin_rebalancedLocalNNZ    ) &&
        (overallMax_rebalancedLocalNNZ     <= overallMax_originalLocalNNZ      )) {
      // Ok
    }
    else {
      allOk1 = false;
    }
  }

  if (localProc == 0) {
    std::cout << "allOk1 = " << allOk1 << std::endl;
    std::cout.flush();
  }
  comm->barrier();
  
  // ****************************************************************
  // Step 6/10: compare matrix norms (before and after redistribution)
  // ****************************************************************
  bool allOk2(true);
  {
    Scalar_t originalFrobNorm  ( originalMatrix->getFrobeniusNorm()   );
    Scalar_t rebalancedFrobNorm( rebalancedMatrix->getFrobeniusNorm() );

    Scalar_t relativeDiff = (rebalancedFrobNorm - originalFrobNorm) / originalFrobNorm;
  
    if (localProc == 0) {
      std::cout << "||originalMat||_frob = " << originalFrobNorm
                << ", ||rebalancedMat||_2 = " << rebalancedFrobNorm
                << ", relativeDiff = " << relativeDiff
                << ", Scalar_t = " << typeid(Scalar_t).name()
                << ", Teuchos::ScalarTraits<Scalar_t>::eps() = " << Teuchos::ScalarTraits<Scalar_t>::eps()
                << std::endl;
      std::cout.flush();
    }
    allOk2 = ( std::fabs(relativeDiff) < 10. * Teuchos::ScalarTraits<Scalar_t>::eps() );
  }

  if (localProc == 0) {
    std::cout << "allOk2 = " << allOk2 << std::endl;
    std::cout.flush();
  }
  comm->barrier();

  // ****************************************************************
  // Step 7/10: compare lhs norms (before and after redistribution)
  // ****************************************************************
  bool allOk3(true);
  {
    std::vector<Scalar_t> originalLhsNorms2_vec( numVectors );
    Teuchos::ArrayView<Scalar_t> originalLhsNorms2_array( originalLhsNorms2_vec.data(), originalLhsNorms2_vec.size() );
    originalLhs->norm2(originalLhsNorms2_array);

    std::vector<Scalar_t> rebalancedLhsNorms2_vec( numVectors );
    Teuchos::ArrayView<Scalar_t> rebalancedLhsNorms2_array( rebalancedLhsNorms2_vec.data(), rebalancedLhsNorms2_vec.size() );
    rebalancedLhs->norm2(rebalancedLhsNorms2_array);

    for (size_t v(0); (v < numVectors) && allOk3; ++v) {
      Scalar_t originalNorm2   = originalLhsNorms2_array[v];
      Scalar_t rebalancedNorm2 = rebalancedLhsNorms2_array[v];

      Scalar_t relativeDiff = (rebalancedNorm2 - originalNorm2) / originalNorm2;
  
      if (localProc == 0) {
        std::cout << "||originalLhs[" << v << "]||_2 = " << originalNorm2
                  << ", ||rebalancedLhs[" << v << "]||_2 = " << rebalancedNorm2
                  << ", relativeDiff = " << relativeDiff
                  << ", Scalar_t = " << typeid(Scalar_t).name()
                  << ", Teuchos::ScalarTraits<Scalar_t>::eps() = " << Teuchos::ScalarTraits<Scalar_t>::eps()
                  << std::endl;
        std::cout.flush();
      }
      allOk3 = ( std::fabs(relativeDiff) < 10. * Teuchos::ScalarTraits<Scalar_t>::eps() );
    }
  }

  if (localProc == 0) {
    std::cout << "allOk3 = " << allOk3 << std::endl;
    std::cout.flush();
  }
  comm->barrier();

  // ****************************************************************
  // Step 8/10: compare rhs norms (before and after redistribution)
  // ****************************************************************
  bool allOk4(true);
  {
    std::vector<Scalar_t> originalRhsNorms2_vec( numVectors );
    Teuchos::ArrayView<Scalar_t> originalRhsNorms2_array( originalRhsNorms2_vec.data(), originalRhsNorms2_vec.size() );
    originalRhs->norm2(originalRhsNorms2_array);

    std::vector<Scalar_t> rebalancedRhsNorms2_vec( numVectors );
    Teuchos::ArrayView<Scalar_t> rebalancedRhsNorms2_array( rebalancedRhsNorms2_vec.data(), rebalancedRhsNorms2_vec.size() );
    rebalancedRhs->norm2(rebalancedRhsNorms2_array);

    for (size_t v(0); (v < numVectors) && allOk4; ++v) {
      Scalar_t originalNorm2   = originalRhsNorms2_array[v];
      Scalar_t rebalancedNorm2 = rebalancedRhsNorms2_array[v];

      Scalar_t relativeDiff = (rebalancedNorm2 - originalNorm2) / originalNorm2;
  
      if (localProc == 0) {
        std::cout << "||originalRhs[" << v << "]||_2 = " << originalNorm2
                  << ", ||rebalancedRhs[" << v << "]||_2 = " << rebalancedNorm2
                  << ", relativeDiff = " << relativeDiff
                  << ", Scalar_t = " << typeid(Scalar_t).name()
                  << ", Teuchos::ScalarTraits<Scalar_t>::eps() = " << Teuchos::ScalarTraits<Scalar_t>::eps()
                  << std::endl;
        std::cout.flush();
      }
      allOk4 = ( std::fabs(relativeDiff) < 10. * Teuchos::ScalarTraits<Scalar_t>::eps() );
    }
  }

  if (localProc == 0) {
    std::cout << "allOk4 = " << allOk4 << std::endl;
    std::cout.flush();
  }
  comm->barrier();

  // ****************************************************************
  // Step 9/10: compare mat * lhs norms (before and after redistribution)
  // ****************************************************************
  bool allOk5(true);
  {
    std::vector<Scalar_t> originalMatLhsNorms2_vec( numVectors );
    Teuchos::ArrayView<Scalar_t> originalMatLhsNorms2_array( originalMatLhsNorms2_vec.data(), originalMatLhsNorms2_vec.size() );

    Teuchos::RCP<MultiV_t> originalMatLhs = Tpetra::createMultiVector<Scalar_t,LocalId_t,GlobalId_t>( originalMatrix->getRangeMap(), numVectors );
    originalMatrix->apply(*originalLhs, *originalMatLhs);
    originalMatLhs->norm2(originalMatLhsNorms2_array);

    std::vector<Scalar_t> rebalancedMatLhsNorms2_vec( numVectors );
    Teuchos::ArrayView<Scalar_t> rebalancedMatLhsNorms2_array( rebalancedMatLhsNorms2_vec.data(), rebalancedMatLhsNorms2_vec.size() );

    Teuchos::RCP<MultiV_t> rebalancedMatLhs = Tpetra::createMultiVector<Scalar_t,LocalId_t,GlobalId_t>( rebalancedMatrix->getRangeMap(), numVectors );
    rebalancedMatrix->apply(*rebalancedLhs, *rebalancedMatLhs);
    rebalancedMatLhs->norm2(rebalancedMatLhsNorms2_array);

    for (size_t v(0); (v < numVectors) && allOk5; ++v) {
      Scalar_t originalNorm2   = originalMatLhsNorms2_array[v];
      Scalar_t rebalancedNorm2 = rebalancedMatLhsNorms2_array[v];

      Scalar_t relativeDiff = (rebalancedNorm2 - originalNorm2) / originalNorm2;
  
      if (localProc == 0) {
        std::cout << "||originalMatLhs[" << v << "]||_2 = " << originalNorm2
                  << ", ||rebalancedMatLhs[" << v << "]||_2 = " << rebalancedNorm2
                  << ", relativeDiff = " << relativeDiff
                  << ", Scalar_t = " << typeid(Scalar_t).name()
                  << ", Teuchos::ScalarTraits<Scalar_t>::eps() = " << Teuchos::ScalarTraits<Scalar_t>::eps()
                  << std::endl;
        std::cout.flush();
      }
      allOk5 = ( std::fabs(relativeDiff) < 10. * Teuchos::ScalarTraits<Scalar_t>::eps() );
    }
  }

  if (localProc == 0) {
    std::cout << "allOk5 = " << allOk5 << std::endl;
    std::cout.flush();
  }
  comm->barrier();

  // ****************************************************************
  // Step 10/10: check that rebalancedRhs == rebalanceMat * rebalancedLhs
  // ****************************************************************
  bool allOk6(true);
  {
    std::vector<Scalar_t> rebalancedRhsNorms2_vec( numVectors );
    Teuchos::ArrayView<Scalar_t> rebalancedRhsNorms2_array( rebalancedRhsNorms2_vec.data(), rebalancedRhsNorms2_vec.size() );
    rebalancedRhs->norm2(rebalancedRhsNorms2_array);

    Teuchos::RCP<MultiV_t> rebalancedMatLhs = Tpetra::createMultiVector<Scalar_t,LocalId_t,GlobalId_t>( rebalancedMatrix->getRangeMap(), numVectors );
    rebalancedMatrix->apply(*rebalancedLhs, *rebalancedMatLhs);
    MultiV_t diff(*rebalancedRhs);
    diff.update(-1., *rebalancedMatLhs, 1.); // diff = 1. * diff - 1. * rebalancedMatLhs

    std::vector<Scalar_t> diffNorms2_vec( numVectors );
    Teuchos::ArrayView<Scalar_t> diffNorms2_array( diffNorms2_vec.data(), diffNorms2_vec.size() );
    diff.norm2(diffNorms2_array);

    for (size_t v(0); (v < numVectors) && allOk6; ++v) {
      Scalar_t diffNorm2       = diffNorms2_array[v];
      Scalar_t rebalancedNorm2 = rebalancedRhsNorms2_array[v];

      Scalar_t ratio = diffNorm2 / rebalancedNorm2;

      if (localProc == 0) {
        std::cout << "||rabalancedRhs[" << v << "] - rebalancedMat * rebalancedLhs[" << v << "]||_2 = " << diffNorms2_array[v]
                  << ", ||rebalancedRhs[" << v << "]||_2 = " << rebalancedNorm2
                  << ", ratio = " << ratio
                  << ", Scalar_t = " << typeid(Scalar_t).name()
                  << ", eps = " << Teuchos::ScalarTraits<Scalar_t>::eps()
                  << std::endl;
        std::cout.flush();
      }
      allOk6 = ( std::fabs(ratio) < 10. * Teuchos::ScalarTraits<Scalar_t>::eps() );
    }
  }

  if (localProc == 0) {
    std::cout << "allOk6 = " << allOk6 << std::endl;
    std::cout.flush();
  }
  comm->barrier();

  return (allOk1 && allOk2 && allOk3 && allOk4 && allOk5 && allOk6);
}

int main(int argc, char** argv) {

  Tpetra::ScopeGuard tpetraScope(&argc, &argv);

  { // Avoiding objects to persist after MPI_Finalize or Kokkos::finalize
    Teuchos::RCP< const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
    if ( runExample(comm) ) {
      std::cout << "End Result: TEST PASSED" << std::endl;
    }
    else {
      std::cout << "End Result: TEST FAILED!" << std::endl;
    }
    
  }
  
  return 0;
}

