// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
// This driver constructs a simple tridiagonal matrix and constant RHS,
// and solves this system using the Belos Block GMRES method with a 'long double' 
// ScalarType.
//
// NOTE: No preconditioner is used in this case.
//
//
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include <BelosIteration.hpp>
#include "BelosLinearProblem.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosTypes.hpp"
#include "BelosConfigDefs.hpp"
#include "BelosTpetraAdapter.hpp"
#include "TpetraExt_MatrixMatrix.hpp"

#include <stdexcept>
#include <limits> 

#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Operator.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Details_DefaultTypes.hpp" 

#include "Teuchos_VerboseObject.hpp"

using namespace Teuchos;

int main(int argc, char *argv[]) {
  //
  typedef long double                      scalar_type;
  typedef int                              LO;
  typedef Tpetra::Map<LO, Tpetra::Details::DefaultTypes::global_ordinal_type>              Tpetra_Map;
  typedef Tpetra::CrsMatrix<scalar_type, LO, Tpetra::Details::DefaultTypes::global_ordinal_type>    Tpetra_CrsMatrix;
  typedef Tpetra::Vector<scalar_type, LO, Tpetra::Details::DefaultTypes::global_ordinal_type>       Tpetra_Vector;
  typedef Tpetra::Operator<scalar_type>             OP;
  typedef Tpetra::MultiVector<scalar_type>          MV;
  typedef Belos::LinearProblem<scalar_type, MV, OP> problem_type;

  Teuchos::GlobalMPISession session(&argc, &argv, NULL);

  bool success = false;
    
  auto comm = Tpetra::getDefaultComm ();

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  
  try {
    // The number of rows and columns in the matrix.
    const Tpetra::global_size_t numGblIndices = 50;

    // Construct a Map that puts approximately the same number of
    // equations on each processor.
    Teuchos::RCP<const Tpetra_Map> map = Teuchos::rcp (new Tpetra_Map (numGblIndices, 0, comm));
    auto numMyElements = map->getLocalNumElements ();

    // Create a Tpetra sparse matrix whose rows have distribution
    // given by the Map.  We expect at most three entries per row.
    Teuchos::RCP<Tpetra_CrsMatrix> A (new Tpetra_CrsMatrix (map, 3));
    // Fill the sparse matrix, one row at a time.
    const scalar_type two = static_cast<scalar_type> (2.0);
    const scalar_type negOne = static_cast<scalar_type> (-1.0);
    for (LO lclRow = 0; lclRow < static_cast<LO> (numMyElements); ++lclRow) {
      const Tpetra::Details::DefaultTypes::global_ordinal_type gblRow = map->getGlobalElement (lclRow);
      // A(0, 0:1) = [2, -1]
      if (gblRow == 0) {
        A->insertGlobalValues (gblRow, tuple<Tpetra::Details::DefaultTypes::global_ordinal_type> (gblRow, gblRow + 1),
             tuple<> (two, negOne));
      }
      // A(N-1, N-2:N-1) = [-1, 2]
      else if (static_cast<Tpetra::global_size_t> (gblRow) == numGblIndices - 1) {
        A->insertGlobalValues (gblRow, tuple<Tpetra::Details::DefaultTypes::global_ordinal_type> (gblRow - 1, gblRow),
             tuple<> (negOne, two));
      }
      // A(i, i-1:i+1) = [-1, 2, -1]
      else {
        A->insertGlobalValues (gblRow, tuple<Tpetra::Details::DefaultTypes::global_ordinal_type> (gblRow - 1, gblRow, gblRow + 1),
             tuple<> (negOne, two, negOne));
      }
    }
    // Tell the sparse matrix that we are done adding entries to it.
    A->fillComplete();

    //B = A*A^T
    Teuchos::RCP<Tpetra_CrsMatrix> B = Teuchos::rcp(new Tpetra_CrsMatrix(map, map->getGlobalNumElements()));
    Tpetra::MatrixMatrix::Multiply<scalar_type, LO, Tpetra::Details::DefaultTypes::global_ordinal_type>(*A, false, *A, true, *B);

    // Create vector of ones, b
    Teuchos::RCP<Tpetra_Vector> b (new Tpetra_Vector (map, true));
    b->putScalar(1.0); 

    //Create RHS: c = A*b
    Teuchos::RCP<Tpetra_Vector> c = Teuchos::rcp(new Tpetra_Vector(map, true));
    A->apply(*b, *c, Teuchos::NO_TRANS, Teuchos::ScalarTraits<scalar_type>::one (), Teuchos::ScalarTraits<scalar_type>::zero ());
    
    //Allocate solution vector x 
    Teuchos::RCP<Tpetra_Vector> x (new Tpetra_Vector (map, true));

    //Create linear problem: (A*A^T)*x = A*b.  It should have the same solution as A*x = b 
    //but executes additional capabilities in Tpetra to form.
    Teuchos::RCP<problem_type> my_problem (new problem_type (B, x, c));
    my_problem->setProblem();

    //Create BlockGmres solver
    Belos::BlockGmresSolMgr<scalar_type, MV, OP> my_solver;
    Teuchos::RCP<Teuchos::ParameterList> solverParams = Teuchos::parameterList ();
    my_solver.setParameters(solverParams);
    my_solver.setProblem (my_problem);

    //Perform solve
    my_solver.solve ();

    //Compute norm of solution vector
    scalar_type norm_x = x->norm2(); 
    *out << "mantissa length of long double ST = " << std::numeric_limits<scalar_type>::digits << "\n";
    typedef std::numeric_limits<scalar_type> ldbl; 
    out->precision(ldbl::max_digits10);
    *out << "cout precision = " << ldbl::max_digits10 << "\n"; 
    *out << "||x|| = " << norm_x << "\n";
    scalar_type norm_x_gold = 1695.64442027183031314;
    scalar_type diff = std::abs(norm_x-norm_x_gold); 
    *out << "diff = " << diff << "\n"; 
    if (diff < 1.0e-15) {
      success = true; 
    }  
  }

  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success)

  if (success) {
    *out << "End Result: TEST PASSED\n";
  }
  else {
    *out << "End Result: TEST FAILED\n";
  }

  return ( success ? 0 : 1 );

} // end test_bl_gmres_longdouble.cpp
