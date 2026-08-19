// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//
//  This test uses the MVOPTester.hpp functions to test the Belos adapters
//  to Tpetra and Thyra.
//

#include "BelosConfigDefs.hpp"
#include "BelosMVOPTester.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"

#include "BelosTypes.hpp"
#include "BelosThyraAdapter.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_TpetraLinearOp.hpp"

#include "Galeri_XpetraMaps.hpp"
#include "Galeri_XpetraProblemFactory.hpp"
#include "Galeri_XpetraParameters.hpp"


int main(int argc, char *argv[])
{
  Tpetra::initialize(&argc, &argv);
  int returnCode = 0;
  {
    using map_type = Tpetra::Map<>;
    using mv_type = Tpetra::MultiVector<>;
    using crs_type = Tpetra::CrsMatrix<>;
    using op_type = Tpetra::Operator<>;
    using Scalar = typename crs_type::scalar_type;
    using LocalOrdinal = typename crs_type::local_ordinal_type;
    using GlobalOrdinal = typename crs_type::global_ordinal_type;

    auto comm = Tpetra::getDefaultComm();

    // number of global elements
    int dim = 100;
    int blockSize = 3;

    // PID info
    int MyPID = comm->getRank();

    Teuchos::ParameterList galeriList;
    galeriList.set<GlobalOrdinal>("nx", dim);
    galeriList.set("mx", comm->getSize());

    // Construct a Map that puts approximately the same number of
    // equations on each processor.
    Teuchos::RCP<const map_type> Map = Galeri::Xpetra::CreateMap<LocalOrdinal, GlobalOrdinal, map_type>("Cartesian1D", comm, galeriList);

    // We are building a tridiagonal matrix where each row has (-1 2 -1)
    auto problem = Galeri::Xpetra::BuildProblem<Scalar, LocalOrdinal, GlobalOrdinal, map_type, crs_type, mv_type>("Laplace1D", Map, galeriList);
    auto A = problem->BuildMatrix();

    // Create an Tpetra_MultiVector for an initial std::vector to start the solver.
    // Note that this needs to have the same number of columns as the blocksize.
    auto ivec = Teuchos::rcp( new mv_type(Map, blockSize) );
    ivec->randomize();

    // Create an output manager to handle the I/O from the solver
    auto MyOM = Teuchos::rcp( new Belos::OutputManager<double>( MyPID ) );
    MyOM->setVerbosity( Belos::Errors + Belos::Warnings );

    typedef Thyra::MultiVectorBase<double> TMVB;
    typedef Thyra::LinearOpBase<double>    TLOB;
    // create thyra objects from the Tpetra objects

    // first, a Thyra::VectorSpaceBase
    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > tpetra_vs = Thyra::createVectorSpace<Scalar>(Map);

    // then, a MultiVectorBase (from the Tpetra_MultiVector)
    Teuchos::RCP<Thyra::MultiVectorBase<double> > thyra_ivec = Thyra::createMultiVector(ivec,tpetra_vs);

    // then, a LinearOpBase (from the Tpetra_CrsMatrix)
    Teuchos::RCP<Thyra::LinearOpBase<double> > thyra_op = Thyra::createLinearOp(Teuchos::rcp_dynamic_cast<op_type>(A, true), tpetra_vs, tpetra_vs);

    // test the Thyra adapter multivector
    bool full_success = true;
    bool success = Belos::TestMultiVecTraits<double,TMVB>(MyOM,thyra_ivec);
    full_success &= success;
    if (success) {
      if ( MyPID==0 ) {
        std::cout << "*** ThyraAdapter PASSED TestMultiVecTraits()" << std::endl;
      }
    } else {
      if ( MyPID==0 ) {
        std::cout << "*** ThyraAdapter FAILED TestMultiVecTraits() ***"
        << std::endl << std::endl;
      }
    }

    // test the Thyra adapter operator
    success = Belos::TestOperatorTraits<double,TMVB,TLOB>(MyOM,thyra_ivec,thyra_op);
    full_success &= success;
    if (success) {
      if ( MyPID==0 ) {
        std::cout << "*** ThyraAdapter PASSED TestOperatorTraits()" << std::endl;
      }
    } else {
      if ( MyPID==0 ) {
        std::cout << "*** ThyraAdapter FAILED TestOperatorTraits() ***"
        << std::endl << std::endl;
      }
    }

    if (!full_success) {
      if (MyPID==0)
        std::cout << "End Result: TEST FAILED" << std::endl;
      returnCode = -1;
    } else {
      //
      // Default return value
      //
      if (MyPID==0)
        std::cout << "End Result: TEST PASSED" << std::endl;
      returnCode = 0;
    }

  }
  Tpetra::finalize();
  return returnCode;
}
