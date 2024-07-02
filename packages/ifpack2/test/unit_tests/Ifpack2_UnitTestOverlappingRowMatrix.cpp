// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/*! \file Ifpack2_UnitTestRelaxation.cpp

\brief Ifpack2 Unit test for the OverlappingRowMatrix template.
*/


#include <Teuchos_ConfigDefs.hpp>
#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <iostream>

// Xpetra / Galeri
#ifdef HAVE_IFPACK2_XPETRA
#include <Xpetra_ConfigDefs.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include <Xpetra_Parameters.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_TpetraMap.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraMatrixTypes.hpp>
#endif

#include "Tpetra_Details_residual.hpp"
#include "KokkosSparse_spmv_impl.hpp"

#include <Ifpack2_UnitTestHelpers.hpp>
#include <Ifpack2_OverlappingRowMatrix.hpp>
#include <Ifpack2_CreateOverlapGraph.hpp>

namespace { // (anonymous)

using Teuchos::ArrayRCP;
using Teuchos::Comm;
using Teuchos::outArg;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::reduceAll;
using Teuchos::REDUCE_MIN;
using std::endl;
typedef tif_utest::Node Node;
typedef Tpetra::global_size_t GST;

// Macro used inside the unit test below.  It tests for global error,
// and if so, prints each process' error message and quits the test
// early.
//
// 'out' only prints on Process 0.  It's really not OK for other
// processes to print to stdout, but it usually works and we need to
// do it for debugging.
#define IFPACK2OVERLAPPINGROWMATRIX_REPORT_GLOBAL_ERR( WHAT_STRING ) do { \
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess)); \
  TEST_EQUALITY_CONST( gblSuccess, 1 ); \
  if (gblSuccess != 1) { \
    out << WHAT_STRING << " FAILED on one or more processes!" << endl; \
    for (int p = 0; p < numProcs; ++p) { \
      if (myRank == p && lclSuccess != 1) { \
        std::cout << errStrm.str () << std::flush; \
      } \
      comm->barrier (); \
      comm->barrier (); \
      comm->barrier (); \
    } \
    return; \
  } \
} while (false)




/***********************************************************************************/
template<class MatrixClass, class MultiVectorClass, class ConstMultiVectorClass>
void localReducedMatvec(const MatrixClass & A_lcl,
                        const ConstMultiVectorClass & X_lcl,
                        const int userNumRows,
                        const MultiVectorClass & Y_lcl) {
  using Teuchos::NO_TRANS;

  using execution_space = typename MatrixClass::execution_space;

  if (A_lcl.numRows() == 0 || userNumRows ==0 || userNumRows > A_lcl.numRows()) {
    return;
  }

  int team_size = -1;
  int vector_length = -1;
  int64_t rows_per_thread = -1;
  
  int64_t numLocalRows = userNumRows;
  int64_t myNnz = A_lcl.nnz();

  int64_t rows_per_team = KokkosSparse::Impl::spmv_launch_parameters<execution_space>(numLocalRows, myNnz, rows_per_thread, team_size, vector_length);
  int64_t worksets = (X_lcl.extent (0) + rows_per_team - 1) / rows_per_team;

  using policy_type = typename Kokkos::TeamPolicy<execution_space>;
  using team_member = typename policy_type::member_type;

  using residual_value_type = typename MultiVectorClass::non_const_value_type;
  using KAT = Kokkos::ArithTraits<residual_value_type>;
  using LO = int64_t;

  policy_type policy (1, 1);
  if (team_size < 0) {
    policy = policy_type (worksets, Kokkos::AUTO, vector_length);
  }
  else {
    policy = policy_type (worksets, team_size, vector_length);
  }

  bool is_vector = (X_lcl.extent(1) == 1);

  if(is_vector) {
    // Vector case
    // Kernel interior shamelessly horked from Ifpack2_Details_ScaledDampedResidual_def.hpp
    Kokkos::parallel_for("reduced-mv-vector",policy,KOKKOS_LAMBDA(const team_member& dev) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange (dev, 0, rows_per_team),[&] (const LO loop) {
            const LO lclRow = static_cast<LO> (dev.league_rank ()) * rows_per_team + loop;
            
            if (lclRow >= numLocalRows) {
              return;
            }

            const auto A_row = A_lcl.rowConst(lclRow);
            const LO row_length = A_row.length;
            residual_value_type A_x = KAT::zero ();          
            
            Kokkos::parallel_reduce(Kokkos::ThreadVectorRange (dev, row_length), [&] (const LO iEntry, residual_value_type& lsum) {
                const auto A_val = A_row.value(iEntry);
                lsum += A_val * X_lcl(A_row.colidx(iEntry),0);
              }, A_x);
            Y_lcl(lclRow,0) = A_x;
            
          });//end parallel_for TeamThreadRange
      });//end parallel_for "residual-vector"
  } else {
    // MultiVector case
    // Kernel interior shamelessly horked from Ifpack2_Details_ScaledDampedResidual_def.hpp
    Kokkos::parallel_for("reduced-mv-multivector",policy,KOKKOS_LAMBDA(const team_member& dev) {
        // NOTE: It looks like I should be able to get this data up above, but if I try to
        // we get internal compiler errors.  Who knew that gcc tried to "gimplify"?
        const LO numVectors = static_cast<LO>(X_lcl.extent(1));
        Kokkos::parallel_for(Kokkos::TeamThreadRange (dev, 0, rows_per_team),[=] (const LO loop) {
            const LO lclRow = static_cast<LO> (dev.league_rank ()) * rows_per_team + loop;
            
            if (lclRow >= numLocalRows) {
              return;
            }
            const auto A_row = A_lcl.rowConst(lclRow);
            const LO row_length = A_row.length;
            for(LO v=0; v<numVectors; v++) {
              residual_value_type A_x = KAT::zero ();          
              
              Kokkos::parallel_reduce(Kokkos::ThreadVectorRange (dev, row_length), [&] (const LO iEntry, residual_value_type& lsum) {
                  const auto A_val = A_row.value(iEntry);
                  lsum += A_val * X_lcl(A_row.colidx(iEntry),v);
                }, A_x);
              Y_lcl(lclRow,v) = A_x;
              
            }//end for numVectors
          });//end parallel_for TeamThreadRange
      });//end parallel_for "residual-multivector"
  }// end else                                     
}// end reducedMatvec


    
/***********************************************************************************/
template<class OverlappedMatrixClass, class MultiVectorClass>
void reducedMatvec(const OverlappedMatrixClass & A,
                   const MultiVectorClass & X,
                   const int overlapLevel,
                   MultiVectorClass & Y) {
  using crs_matrix_type = Tpetra::CrsMatrix<typename OverlappedMatrixClass::scalar_type,
    typename OverlappedMatrixClass::local_ordinal_type,
    typename OverlappedMatrixClass::global_ordinal_type,
    typename OverlappedMatrixClass::node_type>;

  // Assumes that X & Y are sufficiently overlapped for this to work
  RCP<const crs_matrix_type> undA = Teuchos::rcp_dynamic_cast<const crs_matrix_type>(A.getUnderlyingMatrix());
  auto hstarts = A.getExtHaloStartsHost();

  if(overlapLevel >= (int) hstarts.size()) 
    throw std::runtime_error("reducedMatvec: Exceeded available overlap");

  auto undA_lcl = undA->getLocalMatrixDevice ();
  auto extA_lcl = A.getExtMatrix()->getLocalMatrixDevice();
  auto X_lcl = X.getLocalViewDevice (Tpetra::Access::ReadOnly);
  auto Y_lcl = Y.getLocalViewDevice (Tpetra::Access::OverwriteAll);
  
  // Do the "Local part"
  auto numLocalRows = undA->getLocalNumRows();
  localReducedMatvec(undA_lcl,X_lcl,numLocalRows,Y_lcl);

  
  // Now, do the "overlapped part"
  if(overlapLevel > 0) {
    int yrange = hstarts[overlapLevel];
    auto Y_ext = Kokkos::subview(Y_lcl,std::make_pair(numLocalRows,numLocalRows+yrange),Kokkos::ALL());
    
    int xlimit = ( (overlapLevel == (int) hstarts.size()-1) ? X_lcl.extent(0) : numLocalRows+hstarts[overlapLevel+1] );
    auto X_ext = Kokkos::subview(X_lcl,std::make_pair(0,xlimit),Kokkos::ALL());
    
    localReducedMatvec(extA_lcl,X_ext,yrange,Y_ext);
  }

}







TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2OverlappingRowMatrix, Test0, Scalar, LO, GO)
{
  out << "Ifpack2::OverlappingRowMatrix unit test" << endl;
  Teuchos::OSTab tab0 (out);

#ifndef HAVE_IFPACK2_XPETRA
  out << "This test requires building with Xpetra enabled." << endl;
#else
  typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node>       CrsType;
  typedef Tpetra::RowMatrix<Scalar,LO,GO,Node>       row_matrix_type;
  typedef Xpetra::TpetraCrsMatrix<Scalar,LO,GO,Node> XCrsType;
  typedef Xpetra::Map<LO,GO,Node>                    XMapType;
  typedef Xpetra::MultiVector<Scalar,LO,GO,Node>     XMVectorType;
  typedef Tpetra::Vector<Scalar,LO,GO,Node>          VectorType;

  int lclSuccess = 1;
  int gblSuccess = 1;
  std::ostringstream errStrm; // for error collection

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();
  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();

  Teuchos::CommandLineProcessor clp;
  Xpetra::Parameters xpetraParameters (clp);
  Teuchos::ParameterList GaleriList;
  int nx = 447;  // ~200K unknowns per MPI rank
  size_t numElementsPerProc = nx*nx;
  GaleriList.set("nx", static_cast<GO> (nx));
  GaleriList.set("ny", static_cast<GO> (nx * numProcs));
  GaleriList.set("n", static_cast<GO> (numElementsPerProc*numProcs));

  // Short circuit --- this test should only be run in parallel.
  if (numProcs == 1) {
    out << "This test is only meaningful if run with multiple MPI processes."
        << endl;
    return;
  }

  const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
  RCP<XMapType> xmap;
  RCP<Galeri::Xpetra::Problem<XMapType, XCrsType, XMVectorType> > Pr;
  RCP<XCrsType> XA;
  RCP<CrsType> A;
  try {
    xmap = Xpetra::MapFactory<LO, GO>::Build (xpetraParameters.GetLib (), INVALID,
                                              numElementsPerProc, 0, comm);
    Pr = Galeri::Xpetra::BuildProblem<Scalar, LO, GO, XMapType, XCrsType, XMVectorType> (std::string ("Laplace2D"), xmap, GaleriList);
    XA = Pr->BuildMatrix ();
    A = XA->getTpetra_CrsMatrixNonConst ();
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Problem setup (Galeri or Xpetra) threw exception: " << e.what () << endl;
  }
  IFPACK2OVERLAPPINGROWMATRIX_REPORT_GLOBAL_ERR( "Problem setup" );

  if (A.is_null ()) {
    lclSuccess = 0;
  }
  IFPACK2OVERLAPPINGROWMATRIX_REPORT_GLOBAL_ERR( "A is null on at least one process." );

  VectorType X (A->getRowMap ());
  VectorType Y (A->getRowMap ());
  VectorType Z (A->getRowMap ());

  const int OverlapLevel = 5;

  // ======================================== //
  // Build the overlapping matrix using class //
  // Ifpack2::OverlappingRowMatrix.           //
  // ======================================== //
  RCP<Ifpack2::OverlappingRowMatrix<row_matrix_type> > B;


  try {
    B = rcp (new Ifpack2::OverlappingRowMatrix<row_matrix_type> (A, OverlapLevel));
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Ifpack2::OverlappingRowMatrix constructor threw exception: " << e.what () << endl;
  }
  IFPACK2OVERLAPPINGROWMATRIX_REPORT_GLOBAL_ERR( "Ifpack2::OverlappingRowMatrix constructor" );

  auto halo = B->getExtHaloStartsHost();
#if 0
  printf("Halo Starts:");
  for(size_t i=0; i< (size_t)halo.size(); i++)
    printf("%d ",(int) halo[i]);
  printf("\n");
#endif

  size_t NumGlobalRowsB = B->getGlobalNumRows ();
  size_t NumGlobalNonzerosB = B->getGlobalNumEntries ();

  {
  ArrayRCP<ArrayRCP<Scalar> > x_ptr = X.get2dViewNonConst ();
  for (LO i = 0 ; i < static_cast<LO> (A->getLocalNumRows ()); ++i) {
    x_ptr[0][i] = 1.0 * A->getRowMap ()->getGlobalElement (i);
  }
  }
  Y.putScalar (0.0);

  VectorType ExtX_B (B->getRowMap ());
  VectorType ExtY_B (B->getRowMap ());
  ExtY_B.putScalar (0.0);

  try {
    B->importMultiVector (X,ExtX_B);
    B->apply (ExtX_B,ExtY_B);
    B->exportMultiVector (ExtY_B, Y, Tpetra::ADD);
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Import, apply B, and Export: " << e.what () << endl;
  }
  IFPACK2OVERLAPPINGROWMATRIX_REPORT_GLOBAL_ERR( "Import, apply B, and Export" );

  // ================================================== //
  // Build the overlapping graph using                  //
  // CreateOverlappingMatrix.                            //
  // ================================================== //
  RCP<const CrsType> C;
  try {
    C = Ifpack2::createOverlapMatrix<CrsType> (A, OverlapLevel);
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Ifpack2::createOverlapMatrix threw an exception" << e.what () << endl;
  }
  IFPACK2OVERLAPPINGROWMATRIX_REPORT_GLOBAL_ERR( "Ifpack2::createOverlapMatrix" );

  // simple checks on global quantities
  const size_t NumGlobalRowsC = C->getGlobalNumRows ();
  const size_t NumGlobalNonzerosC = C->getGlobalNumEntries ();

  TEST_EQUALITY( NumGlobalRowsB, NumGlobalRowsC );
  TEST_EQUALITY( NumGlobalNonzerosB, NumGlobalNonzerosC );

  // Test fix to github issue #558. Check that all four maps report the same
  // number of local elements. This means that LocalFilter can filter based on
  // getDomainMap () and getRangeMap (), as desired, and see the overlap
  // pattern.
  {
    const auto n = B->getRowMap ()->getLocalNumElements ();
    TEST_EQUALITY( B->getColMap ()->getLocalNumElements (), n );
    TEST_EQUALITY( B->getRangeMap ()->getLocalNumElements (), n );
    TEST_EQUALITY( B->getDomainMap ()->getLocalNumElements (), n );
  }

  try {
    C->apply (X, Z);
  } catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Ifpack2::OverlappingRowMatrix::apply threw an exception" << e.what () << endl;
  }
  IFPACK2OVERLAPPINGROWMATRIX_REPORT_GLOBAL_ERR( "Ifpack2::OverlappingRowMatrix::apply" );

  TEST_COMPARE_FLOATING_ARRAYS( Y.get1dView (), Z.get1dView (), 1e4 * Teuchos::ScalarTraits<Scalar>::eps () );
#endif // HAVE_IFPACK2_XPETRA
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2OverlappingRowMatrix, getLocalDiag, Scalar, LO, GO)
{
  typedef Scalar scalar_type;
  typedef LO local_ordinal_type;
  typedef GO global_ordinal_type;
  typedef Tpetra::Map<>::node_type node_type;

  typedef Tpetra::CrsMatrix<scalar_type,
                            local_ordinal_type,
                            global_ordinal_type,
                            node_type> crs_matrix_type;
  typedef Tpetra::RowMatrix<scalar_type,
                            local_ordinal_type,
                            global_ordinal_type,
                            node_type> row_matrix_type;
  typedef Tpetra::Vector<scalar_type,
                         local_ordinal_type,
                         global_ordinal_type,
                         node_type> vec_type;
  typedef Tpetra::Map<local_ordinal_type,
                      global_ordinal_type,
                      node_type> map_type;
  using Teuchos::Array;
  using Teuchos::as;

  int OverlapLevel=2;
  out << "Testing that OverlappingRowMatrix's getLocalDiag method works properly." << std::endl;

  // This test assumes indexBase == 0.

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  // Short circuit --- this test should only be run in parallel.
  if (comm->getSize() == 1) {
    out << "This test is only meaningful if run with multiple MPI processes."
        << endl;
    return;
  }

  const size_t numLocalRows = 4;
  int myRank = comm->getRank();
  const Tpetra::global_size_t globalNumRows = comm->getSize() * numLocalRows;

  RCP<const map_type> rowMap(new map_type(globalNumRows, numLocalRows, 0, comm));
  RCP<const map_type> domMap = rowMap;
  RCP<const map_type> ranMap = rowMap;

  out << "Creating and filling matrix A" << endl;

  RCP<crs_matrix_type> A(new crs_matrix_type(rowMap, 3));

  Array<scalar_type> val(3);
  Array<global_ordinal_type> ind(3);
  scalar_type zero=0.0;
  scalar_type one=1.0;
  val[0] = static_cast<scalar_type>(-1.0);
  val[1] = static_cast<scalar_type>(2.0);
  val[2] = static_cast<scalar_type>(-1.0);

  GO gidOffset;
  GO nlr = static_cast<GO>(numLocalRows);
  Teuchos::scan (*comm, Teuchos::REDUCE_SUM, nlr, Teuchos::outArg (gidOffset));
  gidOffset -= numLocalRows;

  if (myRank == 0) {
    // Row GID 0: Insert [2 -1] into column GIDs [0 1]
    ind[1] = 0;
    ind[2] = 1;
    val[1] = zero;
    val[2] = -one;
    A->insertGlobalValues(0, ind.view(1, 2), val.view(1, 2));

    val[0] = val[2] = -one;
    for (LO i=1; i<static_cast<int>(numLocalRows); ++i) {
      ind[0]=gidOffset+i-1;
      ind[1]=gidOffset+i;
      ind[2]=gidOffset+i+1;
      val[1]=one*(gidOffset+i); //put gid on the diagonal
      A->insertGlobalValues(gidOffset+i, ind.view(0, 3), val.view(0, 3));
    }


  } else if (myRank == comm->getSize()-1) {

    val[0] = val[2] = -one;
    for (LO i=0; i<static_cast<int>(numLocalRows)-1; ++i) {
      ind[0]=gidOffset+i-1;
      ind[1]=gidOffset+i;
      ind[2]=gidOffset+i+1;
      val[1]=one*(gidOffset+i); //diagonal value is the row GID
      A->insertGlobalValues(gidOffset+i, ind.view(0, 3), val.view(0, 3));
    }
    //last global row
    ind[0] = gidOffset+numLocalRows-2;
    ind[1] = gidOffset+numLocalRows-1;
    val[0] = zero;
    val[1] = gidOffset+numLocalRows-1;
    A->insertGlobalValues(gidOffset+numLocalRows-1, ind.view(0, 2), val.view(0, 2));

  }
  else {

    val[0] = val[2] = -one;
    for (LO i=0; i<static_cast<int>(numLocalRows); ++i) {
      ind[0]=gidOffset+i-1;
      ind[1]=gidOffset+i;
      ind[2]=gidOffset+i+1;
      val[1]=one*(gidOffset+i); //diagonal value is the row GID
      A->insertGlobalValues(gidOffset+i, ind.view(0, 3), val.view(0, 3));
    }

  }

  A->fillComplete(domMap, ranMap);
  //out << "---- A matrix -----" << std::endl;
  //out.setOutputToRootOnly(-1);
  //A->describe(out,Teuchos::VERB_EXTREME);
  //out.setOutputToRootOnly(0);
  //out << "-------------------" << std::endl;

  RCP<Ifpack2::OverlappingRowMatrix<row_matrix_type> > ovA;

  out << "Building overlapping matrix" << std::endl;
  try {
    ovA = rcp (new Ifpack2::OverlappingRowMatrix<row_matrix_type> (A, OverlapLevel));
  } catch (std::exception& e) {
    out << "Ifpack2::OverlappingRowMatrix constructor threw exception: " << e.what () << endl;
  }
  RCP<const map_type> ovRowMap = ovA->getRowMap();
  //out << "-- row map --" << std::endl;
  //out.setOutputToRootOnly(-1);
  //ovRowMap->describe(out,Teuchos::VERB_EXTREME);
  //out.setOutputToRootOnly(0);
  //out << "-------------" << std::endl;
  RCP<vec_type> localDiag = rcp(new vec_type(ovRowMap));
  ovA->getLocalDiagCopy(*localDiag);
  //out << "-- diagonal --" << std::endl;
  //out.setOutputToRootOnly(-1);
  //localDiag->describe(out,Teuchos::VERB_EXTREME);
  //out.setOutputToRootOnly(0);
  //out << "--------------" << std::endl;

  Teuchos::Array<Scalar> ldGids(localDiag->getLocalLength());
  localDiag->get1dCopy(ldGids());
  auto ovrmGids = ovRowMap->getMyGlobalIndices();

  TEST_EQUALITY(static_cast<GO>(ldGids.size()),static_cast<GO>(ovrmGids.size()));
  for (size_t i=0; i<ovrmGids.size(); ++i)
    TEST_EQUALITY(ldGids[i],ovrmGids[i]);
}


TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2OverlappingRowMatrix, reducedMatvec, Scalar, LocalOrdinal, GlobalOrdinal)
{
  using SC = Scalar;
  using LO = LocalOrdinal;
  using GO = GlobalOrdinal;
  using NO = Tpetra::Map<>::node_type;
  using map_type = Tpetra::Map<LO,GO,NO>;
  using row_matrix_type =  Tpetra::RowMatrix<SC,LO,GO,NO>;
  using MV = Tpetra::MultiVector<SC,LO,GO,NO>;
  using Teuchos::RCP;
  Tpetra::global_size_t num_rows_per_proc = 5;
  
  const RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);
  // Only run on > 1 core
  if(rowmap->getComm()->getSize() == 1) return;  

  RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > A = tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);

  // This needs to be one less than the number of matvecs we test
  int overlapLevel = 2;
  Ifpack2::OverlappingRowMatrix<row_matrix_type> ovA(A, overlapLevel);

  //RCP<const row_matrix_type> ExtMatrix = ovA.getExtMatrix();
  auto ExtMatrix = ovA.getExtMatrix();
  SC one = Teuchos::ScalarTraits<SC>::one(), zero = Teuchos::ScalarTraits<SC>::zero();

  // Vectors in the non-overlapping space
  int numVecs = 2; 
  MV x(rowmap,numVecs), y_direct(rowmap,numVecs), y_overlap(rowmap,numVecs);
  x.putScalar(one); 

  {
    // Direct approach
    MV temp1(rowmap,numVecs), temp2(rowmap,numVecs);
    A->apply(x,temp1);
    A->apply(temp1,temp2);
    A->apply(temp2,y_direct);
  }

  {
    // Overlap approach
    RCP<const map_type> ovRowmap = ovA.getRowMap();
    RCP<const map_type> ovColmap = ovA.getColMap();
    MV ovX(ovRowmap,numVecs), ovY(ovRowmap,numVecs), temp1(ovRowmap,numVecs), temp2(ovRowmap,numVecs);
    ovX.putScalar(zero);
    auto hstarts = ovA.getExtHaloStartsHost();
    ovA.importMultiVector(x,ovX);
#if 0
    printf("Halo Starts:");
    for(size_t i=0; i< (size_t)hstarts.size(); i++)
      printf("%d ",(int) hstarts[i]);
    printf("\n");
#endif
    //    printf("Before matvec A is (locally)%dx%d x is of size %d, ovX is ov size %d\n",(int)A->getLocalNumRows(),(int)A->getLocalNumCols(),
    //           (int)x.getMap()->getLocalNumElements(),(int)ovX.getMap()->getLocalNumElements());
    //    printf("ovA->getUnderlyingMatrix() is (locally) %dx%d\n",(int)ovA.getUnderlyingMatrix()->getLocalNumRows(),(int)ovA.getUnderlyingMatrix()->getLocalNumCols());

    reducedMatvec(ovA,ovX,2,temp1);
    reducedMatvec(ovA,temp1,1,temp2);
    reducedMatvec(ovA,temp2,0,ovY);

    auto Y_lcl = y_overlap.getLocalViewDevice(Tpetra::Access::OverwriteAll);
    auto ovY_lcl = ovY.getLocalViewDevice(Tpetra::Access::ReadOnly);
    auto ovYsub = Kokkos::subview(ovY_lcl, std::make_pair<int, int>(0, Y_lcl.extent(0)), Kokkos::ALL);
    Kokkos::deep_copy(Y_lcl,ovYsub);

  }


  // Compare solutions
  TEST_COMPARE_FLOATING_ARRAYS( y_direct.get1dView (), y_overlap.get1dView (), 1e4 * Teuchos::ScalarTraits<Scalar>::eps () );
  
}


#define UNIT_TEST_GROUP_SCALAR_ORDINAL( Scalar, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2OverlappingRowMatrix, Test0, Scalar, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2OverlappingRowMatrix, getLocalDiag, Scalar, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2OverlappingRowMatrix, reducedMatvec, Scalar, LO, GO ) 

typedef Tpetra::MultiVector<>::scalar_type default_scalar_type;
typedef Tpetra::MultiVector<>::local_ordinal_type default_local_ordinal_type;
typedef Tpetra::MultiVector<>::global_ordinal_type default_global_ordinal_type;

UNIT_TEST_GROUP_SCALAR_ORDINAL(default_scalar_type, default_local_ordinal_type, default_global_ordinal_type)

} // namespace (anonymous)

