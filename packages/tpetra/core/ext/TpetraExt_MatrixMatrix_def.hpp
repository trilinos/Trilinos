// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
// @HEADER

#ifndef TPETRA_MATRIXMATRIX_DEF_HPP
#define TPETRA_MATRIXMATRIX_DEF_HPP

#include "TpetraExt_MatrixMatrix_decl.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Array.hpp"
#include "Tpetra_Util.hpp"
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "TpetraExt_MMHelpers_def.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_Import_Util.hpp"
#include "Tpetra_Import_Util2.hpp"
#include <algorithm>
#include "Teuchos_FancyOStream.hpp"

//#define COMPUTE_MMM_STATISTICS


/*! \file TpetraExt_MatrixMatrix_def.hpp

    The implementations for the members of class Tpetra::MatrixMatrixMultiply and related non-member constructors.
 */

namespace Tpetra {


namespace MatrixMatrix{


template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void Multiply(
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
  bool transposeA,
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
  bool transposeB,
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
  bool call_FillComplete_on_result,
  const std::string & label)
{

#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix = std::string("TpetraExt ")+ label + std::string(": ");
  using Teuchos::TimeMonitor;
  Teuchos::RCP<Teuchos::TimeMonitor> MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("MMM All Setup"))));
#endif

  //TEUCHOS_FUNC_TIME_MONITOR_DIFF("My Matrix Mult", mmm_multiply);
  typedef CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> Matrix_t;
  //
  //This method forms the matrix-matrix product C = op(A) * op(B), where
  //op(A) == A   if transposeA is false,
  //op(A) == A^T if transposeA is true,
  //and similarly for op(B).
  //

  //A and B should already be Filled.
  //(Should we go ahead and call FillComplete() on them if necessary?
  // or error out? For now, we choose to error out.)
  TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), std::runtime_error, "MatrixMatrix::Multiply(): Matrix A is not fill complete.");
  TEUCHOS_TEST_FOR_EXCEPTION(!B.isFillComplete(), std::runtime_error, "MatrixMatrix::Multiply(): Matrix B is not fill complete.");
  TEUCHOS_TEST_FOR_EXCEPTION(C.isLocallyIndexed() , std::runtime_error, "MatrixMatrix::Multiply(): Result matrix C must not be locally indexed.");

  //Convience typedefs
  typedef CrsMatrixStruct<
    Scalar,
    LocalOrdinal,
    GlobalOrdinal,
    Node> CrsMatrixStruct_t;
  typedef Map<LocalOrdinal, GlobalOrdinal, Node> Map_t;

  RCP<const Matrix_t > Aprime = null;
  RCP<const Matrix_t > Bprime = null;

  // Is this a "clean" matrix
  bool NewFlag=!C.getGraph()->isLocallyIndexed() && !C.getGraph()->isGloballyIndexed();

  bool use_optimized_ATB=false;
  if(transposeA && !transposeB && call_FillComplete_on_result && NewFlag) {
    use_optimized_ATB=true;
  }
#ifdef USE_OLD_TRANSPOSE // NOTE: For Grey Ballard's use.  Remove this later.
  use_optimized_ATB=false;
#endif


  if(!use_optimized_ATB && transposeA) {
    RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node> at (Teuchos::rcpFromRef (A));
    Aprime = at.createTranspose();
  }
  else{
    Aprime = rcpFromRef(A);
  }

  if(transposeB){
    RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node> bt (Teuchos::rcpFromRef (B));
    Bprime=bt.createTranspose();
  }
  else{
    Bprime = rcpFromRef(B);
  }


  //now check size compatibility
  global_size_t numACols = A.getDomainMap()->getGlobalNumElements();
  global_size_t numBCols = B.getDomainMap()->getGlobalNumElements();
  global_size_t Aouter = transposeA ? numACols : A.getGlobalNumRows();
  global_size_t Bouter = transposeB ? B.getGlobalNumRows() : numBCols;
  global_size_t Ainner = transposeA ? A.getGlobalNumRows() : numACols;
  global_size_t Binner = transposeB ? numBCols : B.getGlobalNumRows();
  TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), std::runtime_error,
    "MatrixMatrix::Multiply: ERROR, inner dimensions of op(A) and op(B) "
    "must match for matrix-matrix product. op(A) is "
    <<Aouter<<"x"<<Ainner << ", op(B) is "<<Binner<<"x"<<Bouter<<std::endl);

  //The result matrix C must at least have a row-map that reflects the
  //correct row-size. Don't check the number of columns because rectangular
  //matrices which were constructed with only one map can still end up
  //having the correct capacity and dimensions when filled.
  TEUCHOS_TEST_FOR_EXCEPTION(Aouter > C.getGlobalNumRows(), std::runtime_error,
    "MatrixMatrix::Multiply: ERROR, dimensions of result C must "
    "match dimensions of op(A) * op(B). C has "<<C.getGlobalNumRows()
     << " rows, should have at least "<<Aouter << std::endl);

  //It doesn't matter whether C is already Filled or not. If it is already
  //Filled, it must have space allocated for the positions that will be
  //referenced in forming C = op(A)*op(B). If it doesn't have enough space,
  //we'll error out later when trying to store result values.

  // CGB: However, matrix must be in active-fill
  TEUCHOS_TEST_FOR_EXCEPT( C.isFillActive() == false );

  //We're going to need to import remotely-owned sections of A and/or B
  //if more than 1 processor is performing this run, depending on the scenario.
  int numProcs = A.getComm()->getSize();

  //Declare a couple of structs that will be used to hold views of the data
  //of A and B, to be used for fast access during the matrix-multiplication.
  CrsMatrixStruct_t Aview;
  CrsMatrixStruct_t Bview;

  RCP<const Map_t > targetMap_A = Aprime->getRowMap();
  RCP<const Map_t > targetMap_B = Bprime->getRowMap();

#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("MMM All I&X"))));
#endif

  //Now import any needed remote rows and populate the Aview struct.
  // NOTE: We assert that an import isn't needed --- since we do the transpose above to handle that.
  if(!use_optimized_ATB) {
    RCP<const Import<LocalOrdinal,GlobalOrdinal, Node> > dummyImporter;
    MMdetails::import_and_extract_views(*Aprime, targetMap_A, Aview,dummyImporter,true,label);
  }


  //We will also need local access to all rows of B that correspond to the
  //column-map of op(A).
  if (numProcs > 1) {
    targetMap_B = Aprime->getColMap(); //colmap_op_A;
  }

  //Now import any needed remote rows and populate the Bview struct.
  if(!use_optimized_ATB)
    MMdetails::import_and_extract_views(*Bprime, targetMap_B, Bview, Aprime->getGraph()->getImporter(),false,label);

#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("MMM All Multiply"))));
#endif


  //Now call the appropriate method to perform the actual multiplication.
  CrsWrapper_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crsmat(C);

  if(use_optimized_ATB) {
    MMdetails::mult_AT_B_newmatrix(A, B, C,label);
  }
  else if(call_FillComplete_on_result && NewFlag ) {
    MMdetails::mult_A_B_newmatrix(Aview, Bview, C,label);
  }
  else {
    MMdetails::mult_A_B(Aview, Bview, crsmat,label);
#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("MMM All FillComplete"))));
#endif
    if (call_FillComplete_on_result) {
      //We'll call FillComplete on the C matrix before we exit, and give
      //it a domain-map and a range-map.
      //The domain-map will be the domain-map of B, unless
      //op(B)==transpose(B), in which case the range-map of B will be used.
      //The range-map will be the range-map of A, unless
      //op(A)==transpose(A), in which case the domain-map of A will be used.
      if (!C.isFillComplete()) {
        C.fillComplete(Bprime->getDomainMap(), Aprime->getRangeMap());
      }
    }
  }

}


template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void Jacobi(Scalar omega,
            const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> & Dinv,
            const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
            const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
            CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
            bool call_FillComplete_on_result,
	    const std::string & label)
{
#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix = std::string("TpetraExt ")+ label + std::string(": ");
  using Teuchos::TimeMonitor;
  Teuchos::RCP<Teuchos::TimeMonitor> MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix+std::string("Jacobi All Setup"))));
#endif
  typedef CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> Matrix_t;

  //A and B should already be Filled.
  //(Should we go ahead and call FillComplete() on them if necessary?
  // or error out? For now, we choose to error out.)
  TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), std::runtime_error, "MatrixMatrix::Multiply(): Matrix A is not fill complete.");
  TEUCHOS_TEST_FOR_EXCEPTION(!B.isFillComplete(), std::runtime_error, "MatrixMatrix::Multiply(): Matrix B is not fill complete.");
  TEUCHOS_TEST_FOR_EXCEPTION(C.isLocallyIndexed() , std::runtime_error, "MatrixMatrix::Multiply(): Result matrix C must not be locally indexed.");

  //Convience typedefs
  typedef CrsMatrixStruct<
    Scalar,
    LocalOrdinal,
    GlobalOrdinal,
    Node> CrsMatrixStruct_t;
  typedef Map<LocalOrdinal, GlobalOrdinal, Node> Map_t;

  RCP<const Matrix_t > Aprime = rcpFromRef(A);
  RCP<const Matrix_t > Bprime = rcpFromRef(B);

  //now check size compatibility
  global_size_t numACols = A.getDomainMap()->getGlobalNumElements();
  global_size_t numBCols = B.getDomainMap()->getGlobalNumElements();
  global_size_t Aouter = A.getGlobalNumRows();
  global_size_t Bouter = numBCols;
  global_size_t Ainner = numACols;
  global_size_t Binner = B.getGlobalNumRows();
  TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), std::runtime_error,
    "MatrixMatrix::Jacobi: ERROR, inner dimensions of op(A) and op(B) "
    "must match for matrix-matrix product. op(A) is "
    <<Aouter<<"x"<<Ainner << ", op(B) is "<<Binner<<"x"<<Bouter<<std::endl);

  //The result matrix C must at least have a row-map that reflects the
  //correct row-size. Don't check the number of columns because rectangular
  //matrices which were constructed with only one map can still end up
  //having the correct capacity and dimensions when filled.
  TEUCHOS_TEST_FOR_EXCEPTION(Aouter > C.getGlobalNumRows(), std::runtime_error,
    "MatrixMatrix::Multiply: ERROR, dimensions of result C must "
    "match dimensions of op(A) * op(B). C has "<<C.getGlobalNumRows()
     << " rows, should have at least "<<Aouter << std::endl);

  //It doesn't matter whether C is already Filled or not. If it is already
  //Filled, it must have space allocated for the positions that will be
  //referenced in forming C = op(A)*op(B). If it doesn't have enough space,
  //we'll error out later when trying to store result values.

  // CGB: However, matrix must be in active-fill
  TEUCHOS_TEST_FOR_EXCEPT( C.isFillActive() == false );

  //We're going to need to import remotely-owned sections of A and/or B
  //if more than 1 processor is performing this run, depending on the scenario.
  int numProcs = A.getComm()->getSize();

  //Declare a couple of structs that will be used to hold views of the data
  //of A and B, to be used for fast access during the matrix-multiplication.
  CrsMatrixStruct_t Aview;
  CrsMatrixStruct_t Bview;

  RCP<const Map_t > targetMap_A = Aprime->getRowMap();
  RCP<const Map_t > targetMap_B = Bprime->getRowMap();

#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix+std::string("Jacobi All I&X"))));
#endif

  //Now import any needed remote rows and populate the Aview struct.
  RCP<const Import<LocalOrdinal,GlobalOrdinal, Node> > dummyImporter;
  MMdetails::import_and_extract_views(*Aprime, targetMap_A, Aview,dummyImporter,false,label);

  //We will also need local access to all rows of B that correspond to the
  //column-map of op(A).
  if (numProcs > 1) {
    targetMap_B = Aprime->getColMap(); //colmap_op_A;
  }

  //Now import any needed remote rows and populate the Bview struct.
  MMdetails::import_and_extract_views(*Bprime, targetMap_B, Bview, Aprime->getGraph()->getImporter(),false,label);

#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix+std::string("Jacobi All Multiply"))));
#endif


  //Now call the appropriate method to perform the actual multiplication.
  CrsWrapper_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crsmat(C);

  // Is this a "clean" matrix
  bool NewFlag=!C.getGraph()->isLocallyIndexed() && !C.getGraph()->isGloballyIndexed();

  if(call_FillComplete_on_result && NewFlag ) {
    MMdetails::jacobi_A_B_newmatrix(omega,Dinv,Aview, Bview, C,label);
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::runtime_error,
      "jacobi_A_B_general not implemented");
    // FIXME (mfh 03 Apr 2014) This statement is unreachable, so I'm
    // commenting it out.
// #ifdef HAVE_TPETRA_MMM_TIMINGS
//     MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer("TpetraExt: Jacobi FillComplete")));
// #endif
    // FIXME (mfh 03 Apr 2014) This statement is unreachable, so I'm
    // commenting it out.
    // if (call_FillComplete_on_result) {
    //   //We'll call FillComplete on the C matrix before we exit, and give
    //   //it a domain-map and a range-map.
    //   //The domain-map will be the domain-map of B, unless
    //   //op(B)==transpose(B), in which case the range-map of B will be used.
    //   //The range-map will be the range-map of A, unless
    //   //op(A)==transpose(A), in which case the domain-map of A will be used.
    //   if (!C.isFillComplete()) {
    //     C.fillComplete(Bprime->getDomainMap(), Aprime->getRangeMap());
    //   }
    // }
  }
}



template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void Add(
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
  bool transposeA,
  Scalar scalarA,
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
  Scalar scalarB )
{
  TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), std::runtime_error,
    "MatrixMatrix::Add ERROR, input matrix A.isFillComplete() is false; it is required to be true. (Result matrix B is not required to be isFillComplete()).");
  TEUCHOS_TEST_FOR_EXCEPTION(B.isFillComplete() , std::runtime_error,
    "MatrixMatrix::Add ERROR, input matrix B must not be fill complete!");
  TEUCHOS_TEST_FOR_EXCEPTION(B.isStaticGraph() , std::runtime_error,
    "MatrixMatrix::Add ERROR, input matrix B must not have static graph!");
  TEUCHOS_TEST_FOR_EXCEPTION(B.isLocallyIndexed() , std::runtime_error,
    "MatrixMatrix::Add ERROR, input matrix B must not be locally indexed!");
  TEUCHOS_TEST_FOR_EXCEPTION(B.getProfileType()!=DynamicProfile, std::runtime_error,
    "MatrixMatrix::Add ERROR, input matrix B must have a dynamic profile!");
  //Convience typedef
  typedef CrsMatrix<
    Scalar,
    LocalOrdinal,
    GlobalOrdinal,
    Node> CrsMatrix_t;
  RCP<const CrsMatrix_t> Aprime = null;
  if( transposeA ){
          RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node> theTransposer(Teuchos::rcpFromRef (A));
    Aprime = theTransposer.createTranspose();
  }
  else{
    Aprime = rcpFromRef(A);
  }
  size_t a_numEntries;
  Array<GlobalOrdinal> a_inds(A.getNodeMaxNumRowEntries());
  Array<Scalar> a_vals(A.getNodeMaxNumRowEntries());
  GlobalOrdinal row;

  if(scalarB != ScalarTraits<Scalar>::one()){
    B.scale(scalarB);
  }

  bool bFilled = B.isFillComplete();
  size_t numMyRows = B.getNodeNumRows();
  if(scalarA != ScalarTraits<Scalar>::zero()){
    for(LocalOrdinal i = 0; (size_t)i < numMyRows; ++i){
      row = B.getRowMap()->getGlobalElement(i);
      Aprime->getGlobalRowCopy(row, a_inds(), a_vals(), a_numEntries);
      if(scalarA != ScalarTraits<Scalar>::one()){
        for(size_t j =0; j<a_numEntries; ++j){
          a_vals[j] *= scalarA;
        }
      }
      if(bFilled){
        B.sumIntoGlobalValues(row, a_inds(0,a_numEntries), a_vals(0,a_numEntries));
      }
      else{
        B.insertGlobalValues(row, a_inds(0,a_numEntries), a_vals(0,a_numEntries));
      }

    }
  }
}


template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
add (const Scalar& alpha,
     const bool transposeA,
     const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
     const Scalar& beta,
     const bool transposeB,
     const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
     const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& domainMap,
     const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& rangeMap,
     const Teuchos::RCP<Teuchos::ParameterList>& params)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::rcp_implicit_cast;
  typedef RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> row_matrix_type;
  typedef CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crs_matrix_type;
  typedef RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node> transposer_type;

  TEUCHOS_TEST_FOR_EXCEPTION(
    ! A.isFillComplete () || ! B.isFillComplete (), std::invalid_argument,
    "Tpetra::MatrixMatrix::add: A and B must both be fill complete.");

#ifdef HAVE_TPETRA_DEBUG
  // The matrices don't have domain or range Maps unless they are fill complete.
  if (A.isFillComplete () && B.isFillComplete ()) {
    const bool domainMapsSame =
      (! transposeA && ! transposeB && ! A.getDomainMap ()->isSameAs (* (B.getDomainMap ()))) ||
      (! transposeA && transposeB && ! A.getDomainMap ()->isSameAs (* (B.getRangeMap ()))) ||
      (transposeA && ! transposeB && ! A.getRangeMap ()->isSameAs (* (B.getDomainMap ())));
    TEUCHOS_TEST_FOR_EXCEPTION(
      domainMapsSame, std::invalid_argument,
      "Tpetra::MatrixMatrix::add: The domain Maps of Op(A) and Op(B) are not the same.");

    const bool rangeMapsSame =
      (! transposeA && ! transposeB && ! A.getRangeMap ()->isSameAs (* (B.getRangeMap ()))) ||
      (! transposeA && transposeB && ! A.getRangeMap ()->isSameAs (* (B.getDomainMap ()))) ||
      (transposeA && ! transposeB && ! A.getDomainMap ()->isSameAs (* (B.getRangeMap ())));
    TEUCHOS_TEST_FOR_EXCEPTION(
      rangeMapsSame, std::invalid_argument,
      "Tpetra::MatrixMatrix::add: The range Maps of Op(A) and Op(B) are not the same.");
  }
#endif // HAVE_TPETRA_DEBUG

  // Form the explicit transpose of A if necessary.
  RCP<const crs_matrix_type> Aprime;
  if (transposeA) {
    transposer_type theTransposer (rcpFromRef (A));
    Aprime = theTransposer.createTranspose ();
  } else {
    Aprime = rcpFromRef (A);
  }

#ifdef HAVE_TPETRA_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(Aprime.is_null (), std::logic_error,
    "Tpetra::MatrixMatrix::Add: Failed to compute Op(A).  "
    "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG

  // Form the explicit transpose of B if necessary.
  RCP<const crs_matrix_type> Bprime;
  if (transposeB) {
    transposer_type theTransposer (rcpFromRef (B));
    Bprime = theTransposer.createTranspose ();
  } else {
    Bprime = rcpFromRef (B);
  }

#ifdef HAVE_TPETRA_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(Bprime.is_null (), std::logic_error,
    "Tpetra::MatrixMatrix::Add: Failed to compute Op(B).  "
    "Please report this bug to the Tpetra developers.");

  TEUCHOS_TEST_FOR_EXCEPTION(
    ! Aprime->isFillComplete () || ! Bprime->isFillComplete (), std::invalid_argument,
    "Tpetra::MatrixMatrix::add: Aprime and Bprime must both be fill complete.  "
    "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG


  RCP<row_matrix_type> C =
    Bprime->add (alpha, *rcp_implicit_cast<const row_matrix_type> (Aprime),
                 beta, domainMap, rangeMap, params);
  return rcp_dynamic_cast<crs_matrix_type> (C);
}




template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void Add(
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
  bool transposeA,
  Scalar scalarA,
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
  bool transposeB,
  Scalar scalarB,
  RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > C)
{
  using Teuchos::as;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::rcpFromRef;
  using Teuchos::tuple;
  using std::endl;
  //  typedef typename ArrayView<const Scalar>::size_type size_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef Map<LocalOrdinal, GlobalOrdinal, Node> map_type;
  //  typedef Import<LocalOrdinal, GlobalOrdinal, Node> import_type;
  //  typedef RowGraph<LocalOrdinal, GlobalOrdinal, Node> row_graph_type;
  //  typedef CrsGraph<LocalOrdinal, GlobalOrdinal, Node> crs_graph_type;
  typedef CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crs_matrix_type;
  typedef RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node> transposer_type;

  TEUCHOS_TEST_FOR_EXCEPTION(C.is_null (), std::logic_error,
    "Tpetra::MatrixMatrix::Add: The case C == null does not actually work.  "
    "Fixing this will require an interface change.");

  TEUCHOS_TEST_FOR_EXCEPTION(
    ! A.isFillComplete () || ! B.isFillComplete (), std::invalid_argument,
    "Tpetra::MatrixMatrix::Add: Both input matrices must be fill complete "
    "before calling this function.");

#ifdef HAVE_TPETRA_DEBUG
  {
    const bool domainMapsSame =
      (! transposeA && ! transposeB && ! A.getDomainMap ()->isSameAs (* (B.getDomainMap ()))) ||
      (! transposeA && transposeB && ! A.getDomainMap ()->isSameAs (* (B.getRangeMap ()))) ||
      (transposeA && ! transposeB && ! A.getRangeMap ()->isSameAs (* (B.getDomainMap ())));
    TEUCHOS_TEST_FOR_EXCEPTION(
      domainMapsSame, std::invalid_argument,
      "Tpetra::MatrixMatrix::Add: The domain Maps of Op(A) and Op(B) are not the same.");

    const bool rangeMapsSame =
      (! transposeA && ! transposeB && ! A.getRangeMap ()->isSameAs (* (B.getRangeMap ()))) ||
      (! transposeA && transposeB && ! A.getRangeMap ()->isSameAs (* (B.getDomainMap ()))) ||
      (transposeA && ! transposeB && ! A.getDomainMap ()->isSameAs (* (B.getRangeMap ())));
    TEUCHOS_TEST_FOR_EXCEPTION(
      rangeMapsSame, std::invalid_argument,
      "Tpetra::MatrixMatrix::Add: The range Maps of Op(A) and Op(B) are not the same.");
  }
#endif // HAVE_TPETRA_DEBUG

  // Form the explicit transpose of A if necessary.
  RCP<const crs_matrix_type> Aprime;
  if (transposeA) {
    transposer_type theTransposer (rcpFromRef (A));
    Aprime = theTransposer.createTranspose ();
  } else {
    Aprime = rcpFromRef (A);
  }

#ifdef HAVE_TPETRA_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(Aprime.is_null (), std::logic_error,
    "Tpetra::MatrixMatrix::Add: Failed to compute Op(A).  "
    "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG

  // Form the explicit transpose of B if necessary.
  RCP<const crs_matrix_type> Bprime;
  if (transposeB) {
    transposer_type theTransposer (rcpFromRef (B));
    Bprime = theTransposer.createTranspose ();
  } else {
    Bprime = rcpFromRef (B);
  }

#ifdef HAVE_TPETRA_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(Bprime.is_null (), std::logic_error,
    "Tpetra::MatrixMatrix::Add: Failed to compute Op(B).  "
    "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG

  // Allocate or zero the entries of the result matrix.
  if (! C.is_null ()) {
    C->setAllToScalar (STS::zero ());
  } else {
#if 0
    // If Aprime and Bprime have the same row Map, and if C is null,
    // we can optimize construction and fillComplete of C.  For now,
    // we just check pointer equality, to avoid the all-reduce in
    // isSameAs.  It may be worth that all-reduce to check, however.
    //if (Aprime->getRowMap ().getRawPtr () == Bprime->getRowMap ().getRawPtr ()) {
    if (Aprime->getRowMap ()->isSameAs (* (Bprime->getRowMap ())) {
      RCP<const map_type> rowMap = Aprime->getRowMap ();

      RCP<const crs_graph_type> A_graph =
        rcp_dynamic_cast<const crs_graph_type> (Aprime->getGraph (), true);
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(A_graph.is_null (), std::logic_error,
        "Tpetra::MatrixMatrix::Add: Graph of Op(A) is null.  "
        "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
      RCP<const crs_graph_type> B_graph =
        rcp_dynamic_cast<const crs_graph_type> (Bprime->getGraph (), true);
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(B_graph.is_null (), std::logic_error,
        "Tpetra::MatrixMatrix::Add: Graph of Op(B) is null.  "
        "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
      RCP<const map_type> A_colMap = A_graph->getColMap ();
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(A_colMap.is_null (), std::logic_error,
        "Tpetra::MatrixMatrix::Add: Column Map of Op(A) is null.  "
        "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
      RCP<const map_type> B_colMap = B_graph->getColMap ();
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(B_colMap.is_null (), std::logic_error,
        "Tpetra::MatrixMatrix::Add: Column Map of Op(B) is null.  "
        "Please report this bug to the Tpetra developers.");
      TEUCHOS_TEST_FOR_EXCEPTION(A_graph->getImporter ().is_null (),
        std::logic_error,
        "Tpetra::MatrixMatrix::Add: Op(A)'s Import is null.  "
        "Please report this bug to the Tpetra developers.");
      TEUCHOS_TEST_FOR_EXCEPTION(B_graph->getImporter ().is_null (),
        std::logic_error,
        "Tpetra::MatrixMatrix::Add: Op(B)'s Import is null.  "
        "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG

      // Compute the (column Map and) Import of the matrix sum.
      RCP<const import_type> sumImport =
        A_graph->getImporter ()->setUnion (* (B_graph->getImporter ()));
      RCP<const map_type> C_colMap = sumImport->getTargetMap ();

      // First, count the number of entries in each row.  Then, go
      // back over the rows again, and compute the actual sum.
      // Remember that C may have a different column Map than Aprime
      // or Bprime, so its local indices may be different.  That's why
      // we have to convert from local to global indices.

      ArrayView<const LocalOrdinal> A_local_ind;
      Array<GlobalOrdinal> A_global_ind;
      ArrayView<const LocalOrdinal> B_local_ind;
      Array<GlobalOrdinal> B_global_ind;

      const size_t localNumRows = rowMap->getNodeNumElements ();
      ArrayRCP<size_t> numEntriesPerRow (localNumRows);
      // Compute the max number of entries in any row of A+B on this
      // process, so that we won't have to resize temporary arrays.
      size_t maxNumEntriesPerRow = 0;
      for (size_t localRow = 0; localRow < localNumRows; ++localRow) {
        // Get view of current row of A_graph, in its local indices.
        A_graph->getLocalRowView (as<LocalOrdinal> (localRow), A_local_ind);
        const size_type A_numEnt = A_local_ind.size ();
        if (A_numEnt > A_global_ind.size ()) {
          A_global_ind.resize (A_numEnt);
        }
        // Convert A's local indices to global indices.
        for (size_type k = 0; k < A_numEnt; ++k) {
          A_global_ind[k] = A_colMap->getGlobalElement (A_local_ind[k]);
        }

        // Get view of current row of B_graph, in its local indices.
        B_graph->getLocalRowView (as<LocalOrdinal> (localRow), B_local_ind);
        const size_type B_numEnt = B_local_ind.size ();
        if (B_numEnt > B_global_ind.size ()) {
          B_global_ind.resize (B_numEnt);
        }
        // Convert B's local indices to global indices.
        for (size_type k = 0; k < B_numEnt; ++k) {
          B_global_ind[k] = B_colMap->getGlobalElement (B_local_ind[k]);
        }

        // Count the number of entries in the merged row of A + B.
        const size_t curNumEntriesPerRow =
          keyMergeCount (A_global_ind.begin (), A_global_ind.end (),
                         B_global_ind.begin (), B_global_ind.end ());
        numEntriesPerRow[localRow] = curNumEntriesPerRow;
        maxNumEntriesPerRow = std::max (maxNumEntriesPerRow, curNumEntriesPerRow);
      }

      // Create C, using the sum column Map and number of entries per
      // row that we computed above.  Having the exact number of
      // entries per row lets us use static profile, making it valid
      // to call expertStaticFillComplete.
      C = rcp (new crs_matrix_type (rowMap, C_colMap, numEntriesPerRow, StaticProfile));

      // Go back through the rows and actually compute the sum.  We
      // don't ever have to resize A_global_ind or B_global_ind below,
      // since we've already done it above.
      ArrayView<const Scalar> A_val;
      ArrayView<const Scalar> B_val;

      Array<LocalOrdinal> AplusB_local_ind (maxNumEntriesPerRow);
      Array<GlobalOrdinal> AplusB_global_ind (maxNumEntriesPerRow);
      Array<Scalar> AplusB_val (maxNumEntriesPerRow);

      for (size_t localRow = 0; localRow < localNumRows; ++localRow) {
        // Get view of current row of A, in A's local indices.
        Aprime->getLocalRowView (as<LocalOrdinal> (localRow), A_local_ind, A_val);
        // Convert A's local indices to global indices.
        for (size_type k = 0; k < A_local_ind.size (); ++k) {
          A_global_ind[k] = A_colMap->getGlobalElement (A_local_ind[k]);
        }

        // Get view of current row of B, in B's local indices.
        Bprime->getLocalRowView (as<LocalOrdinal> (localRow), B_local_ind, B_val);
        // Convert B's local indices to global indices.
        for (size_type k = 0; k < B_local_ind.size (); ++k) {
          B_global_ind[k] = B_colMap->getGlobalElement (B_local_ind[k]);
        }

        const size_t curNumEntries = numEntriesPerRow[localRow];
        ArrayView<LocalOrdinal> C_local_ind = AplusB_local_ind (0, curNumEntries);
        ArrayView<GlobalOrdinal> C_global_ind = AplusB_global_ind (0, curNumEntries);
        ArrayView<Scalar> C_val = AplusB_val (0, curNumEntries);

        // Sum the entries in the current row of A plus B.
        keyValueMerge (A_global_ind.begin (), A_global_ind.end (),
                       A_val.begin (), A_val.end (),
                       B_global_ind.begin (), B_global_ind.end (),
                       B_val.begin (), B_val.end (),
                       C_global_ind.begin (), C_val.begin (),
                       std::plus<Scalar> ());
        // Convert the sum's global indices into C's local indices.
        for (size_type k = 0; k < as<size_type> (numEntriesPerRow[localRow]); ++k) {
          C_local_ind[k] = C_colMap->getLocalElement (C_global_ind[k]);
        }
        // Give the current row sum to C.
        C->replaceLocalValues (localRow, C_local_ind, C_val);
      }

      // Use "expert static fill complete" to bypass construction of
      // the Import and Export (if applicable) object(s).
      RCP<const map_type> domainMap = A_graph->getDomainMap ();
      RCP<const map_type> rangeMap = A_graph->getRangeMap ();
      C->expertStaticFillComplete (domainMap, rangeMap, sumImport, A_graph->getExporter ());

      return; // Now we're done!
    }
    else {
      // FIXME (mfh 08 May 2013) When I first looked at this method, I
      // noticed that C was being given the row Map of Aprime (the
      // possibly transposed version of A).  Is this what we want?
      C = rcp (new crs_matrix_type (Aprime->getRowMap (), null));
    }

#else
    // FIXME (mfh 08 May 2013) When I first looked at this method, I
    // noticed that C was being given the row Map of Aprime (the
    // possibly transposed version of A).  Is this what we want?
    C = rcp (new crs_matrix_type (Aprime->getRowMap (), 0));
#endif // 0
  }

#ifdef HAVE_TPETRA_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(Aprime.is_null (), std::logic_error,
    "Tpetra::MatrixMatrix::Add: At this point, Aprime is null.  "
    "Please report this bug to the Tpetra developers.");
  TEUCHOS_TEST_FOR_EXCEPTION(Bprime.is_null (), std::logic_error,
    "Tpetra::MatrixMatrix::Add: At this point, Bprime is null.  "
    "Please report this bug to the Tpetra developers.");
  TEUCHOS_TEST_FOR_EXCEPTION(C.is_null (), std::logic_error,
    "Tpetra::MatrixMatrix::Add: At this point, C is null.  "
    "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG

  Array<RCP<const crs_matrix_type> > Mat =
    tuple<RCP<const crs_matrix_type> > (Aprime, Bprime);
  Array<Scalar> scalar = tuple<Scalar> (scalarA, scalarB);

  // do a loop over each matrix to add: A reordering might be more efficient
  for (int k = 0; k < 2; ++k) {
    Array<GlobalOrdinal> Indices;
    Array<Scalar> Values;

    // Loop over each locally owned row of the current matrix (either
    // Aprime or Bprime), and sum its entries into the corresponding
    // row of C.  This works regardless of whether Aprime or Bprime
    // has the same row Map as C, because both sumIntoGlobalValues and
    // insertGlobalValues allow summing resp. inserting into nonowned
    // rows of C.
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(Mat[k].is_null (), std::logic_error,
      "Tpetra::MatrixMatrix::Add: At this point, curRowMap is null.  "
      "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
    RCP<const map_type> curRowMap = Mat[k]->getRowMap ();
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(curRowMap.is_null (), std::logic_error,
      "Tpetra::MatrixMatrix::Add: At this point, curRowMap is null.  "
      "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG

    const size_t localNumRows = Mat[k]->getNodeNumRows ();
    for (size_t i = 0; i < localNumRows; ++i) {
      const GlobalOrdinal globalRow = curRowMap->getGlobalElement (i);
      size_t numEntries = Mat[k]->getNumEntriesInGlobalRow (globalRow);
      if (numEntries > 0) {
        Indices.resize (numEntries);
        Values.resize (numEntries);
        Mat[k]->getGlobalRowCopy (globalRow, Indices (), Values (), numEntries);

        if (scalar[k] != STS::one ()) {
          for (size_t j = 0; j < numEntries; ++j) {
            Values[j] *= scalar[k];
          }
        }

        if (C->isFillComplete ()) {
          C->sumIntoGlobalValues (globalRow, Indices, Values);
        } else {
          C->insertGlobalValues (globalRow, Indices, Values);
        }
      }
    }
  }
}



} //End namespace MatrixMatrix

namespace MMdetails{

// Prints MMM-style statistics on communication done with an Import or Export object
template <class TransferType>
void printMultiplicationStatistics(Teuchos::RCP<TransferType > Transfer, const std::string &label){
  if(Transfer.is_null()) return;

  const Distributor & Distor                   = Transfer->getDistributor();
  Teuchos::RCP<const Teuchos::Comm<int> > Comm = Transfer->getSourceMap()->getComm();

  size_t rows_send   = Transfer->getNumExportIDs();
  size_t rows_recv   = Transfer->getNumRemoteIDs();

  size_t round1_send = Transfer->getNumExportIDs() * sizeof(size_t);
  size_t round1_recv = Transfer->getNumRemoteIDs() * sizeof(size_t);
  size_t num_send_neighbors = Distor.getNumSends();
  size_t num_recv_neighbors = Distor.getNumReceives();
  size_t round2_send, round2_recv;
  Distor.getLastDoStatistics(round2_send,round2_recv);

  int myPID    = Comm->getRank();
  int NumProcs = Comm->getSize();

  // Processor by processor statistics
  //    printf("[%d] %s Statistics: neigh[s/r]=%d/%d rows[s/r]=%d/%d r1bytes[s/r]=%d/%d r2bytes[s/r]=%d/%d\n",
  //    myPID,label.c_str(),num_send_neighbors,num_recv_neighbors,rows_send,rows_recv,round1_send,round1_recv,round2_send,round2_recv);

  // Global statistics
  size_t lstats[8] = {num_send_neighbors,num_recv_neighbors,rows_send,rows_recv,round1_send,round1_recv,round2_send,round2_recv};
  size_t gstats_min[8], gstats_max[8];

  double lstats_avg[8], gstats_avg[8];
  for(int i=0; i<8; i++)
    lstats_avg[i] = ((double)lstats[i])/NumProcs;

  Teuchos::reduceAll(*Comm(),Teuchos::REDUCE_MIN,8,lstats,gstats_min);
  Teuchos::reduceAll(*Comm(),Teuchos::REDUCE_MAX,8,lstats,gstats_max);
  Teuchos::reduceAll(*Comm(),Teuchos::REDUCE_SUM,8,lstats_avg,gstats_avg);

  if(!myPID) {
    printf("%s Send Statistics[min/avg/max]: neigh=%d/%4.1f/%d rows=%d/%4.1f/%d round1=%d/%4.1f/%d round2=%d/%4.1f/%d\n",label.c_str(),
           (int)gstats_min[0],gstats_avg[0],(int)gstats_max[0], (int)gstats_min[2],gstats_avg[2],(int)gstats_max[2],
           (int)gstats_min[4],gstats_avg[4],(int)gstats_max[4], (int)gstats_min[6],gstats_avg[6],(int)gstats_max[6]);
    printf("%s Recv Statistics[min/avg/max]: neigh=%d/%4.1f/%d rows=%d/%4.1f/%d round1=%d/%4.1f/%d round2=%d/%4.1f/%d\n",label.c_str(),
           (int)gstats_min[1],gstats_avg[1],(int)gstats_max[1], (int)gstats_min[3],gstats_avg[3],(int)gstats_max[3],
           (int)gstats_min[5],gstats_avg[5],(int)gstats_max[5], (int)gstats_min[7],gstats_avg[7],(int)gstats_max[7]);
  }
}


//kernel method for computing the local portion of C = A*B
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
void mult_AT_B_newmatrix(
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
  const std::string & label) {

  // Using &  Typedefs
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef CrsMatrixStruct<
    Scalar,
    LocalOrdinal,
    GlobalOrdinal,
    Node> CrsMatrixStruct_t;

#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix = std::string("TpetraExt ")+ label + std::string(": ");
  using Teuchos::TimeMonitor;
  Teuchos::RCP<Teuchos::TimeMonitor> MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix+std::string("MMM-T Transpose"))));
#endif

  /*************************************************************/
  /* 1) Local Transpose of A                                   */
  /*************************************************************/
  RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node> at (Teuchos::rcpFromRef (A));
  RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Atrans = at.createTransposeLocal();

  /*************************************************************/
  /* 2/3) Call mult_A_B_newmatrix w/ fillComplete              */
  /*************************************************************/
#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix+std::string("MMM-T I&X"))));
#endif

  // Get views, asserting that no import is required to speed up computation
  CrsMatrixStruct_t Aview;
  CrsMatrixStruct_t Bview;
  RCP<const Import<LocalOrdinal,GlobalOrdinal, Node> > dummyImporter;
  MMdetails::import_and_extract_views(*Atrans, Atrans->getRowMap(), Aview, dummyImporter,true,label);
  MMdetails::import_and_extract_views(B, B.getRowMap(), Bview, dummyImporter,true,label);

#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix+std::string("MMM-T AB-core"))));
#endif

  RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >Ctemp;

  // If Atrans has no Exporter, we can use C instead of having to create a temp matrix
  bool needs_final_export = !Atrans->getGraph()->getExporter().is_null();
  if(needs_final_export)
    Ctemp = rcp(new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(Atrans->getRowMap(),0));
  else
    Ctemp = rcp(&C,false);// don't allow deallocation

  // Multiply
  mult_A_B_newmatrix(Aview,Bview,*Ctemp,label);

  /*************************************************************/
  /* 4) exportAndFillComplete matrix                           */
  /*************************************************************/
#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix+std::string("MMM-T exportAndFillComplete"))));
#endif

  Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Crcp(&C,false);
  if(needs_final_export) {
    Teuchos::ParameterList labelList;
    labelList.set("Timer Label",label);
    Ctemp->exportAndFillComplete(Crcp,*Ctemp->getGraph()->getExporter(),
                                 B.getDomainMap(),A.getDomainMap(),rcp(&labelList,false));
  }
#ifdef COMPUTE_MMM_STATISTICS
  printMultiplicationStatistics(Ctemp->getGraph()->getExporter(),label+std::string(" AT_B MMM"));
#endif
}


//kernel method for computing the local portion of C = A*B
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
void mult_A_B(
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
  CrsWrapper<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
  const std::string & label)
{
  typedef Teuchos::ScalarTraits<Scalar> STS;
  //TEUCHOS_FUNC_TIME_MONITOR_DIFF("mult_A_B", mult_A_B);
  LocalOrdinal C_firstCol = Bview.colMap->getMinLocalIndex();
  LocalOrdinal C_lastCol = Bview.colMap->getMaxLocalIndex();

  LocalOrdinal C_firstCol_import = OrdinalTraits<LocalOrdinal>::zero();
  LocalOrdinal C_lastCol_import = OrdinalTraits<LocalOrdinal>::invalid();

  ArrayView<const GlobalOrdinal> bcols = Bview.colMap->getNodeElementList();
  ArrayView<const GlobalOrdinal> bcols_import = null;
  if (Bview.importColMap != null) {
    C_firstCol_import = Bview.importColMap->getMinLocalIndex();
    C_lastCol_import = Bview.importColMap->getMaxLocalIndex();

    bcols_import = Bview.importColMap->getNodeElementList();
  }

  size_t C_numCols = C_lastCol - C_firstCol +
                        OrdinalTraits<LocalOrdinal>::one();
  size_t C_numCols_import = C_lastCol_import - C_firstCol_import +
                                OrdinalTraits<LocalOrdinal>::one();

  if (C_numCols_import > C_numCols) C_numCols = C_numCols_import;

  Array<Scalar> dwork = Array<Scalar>(C_numCols);
  Array<GlobalOrdinal> iwork = Array<GlobalOrdinal>(C_numCols);
  Array<size_t> iwork2 = Array<size_t>(C_numCols);

  Array<Scalar> C_row_i = dwork;
  Array<GlobalOrdinal> C_cols = iwork;
  Array<size_t> c_index = iwork2;
  Array<GlobalOrdinal> combined_index = Array<GlobalOrdinal>(2*C_numCols);
  Array<Scalar> combined_values = Array<Scalar>(2*C_numCols);

  size_t C_row_i_length, j, k, last_index;

  // Run through all the hash table lookups once and for all
  LocalOrdinal LO_INVALID = OrdinalTraits<LocalOrdinal>::invalid();
  Array<LocalOrdinal> Acol2Brow(Aview.colMap->getNodeNumElements(),LO_INVALID);
  Array<LocalOrdinal> Acol2Irow(Aview.colMap->getNodeNumElements(),LO_INVALID);
  if(Aview.colMap->isSameAs(*Bview.origMatrix->getRowMap())){
    // Maps are the same: Use local IDs as the hash
    for(LocalOrdinal i=Aview.colMap->getMinLocalIndex(); i <=
            Aview.colMap->getMaxLocalIndex(); i++)
      Acol2Brow[i]=i;
  }
  else {
    // Maps are not the same:  Use the map's hash
    for(LocalOrdinal i=Aview.colMap->getMinLocalIndex(); i <=
          Aview.colMap->getMaxLocalIndex(); i++) {
      GlobalOrdinal GID = Aview.colMap->getGlobalElement(i);
      LocalOrdinal BLID = Bview.origMatrix->getRowMap()->getLocalElement(GID);
      if(BLID != LO_INVALID) Acol2Brow[i] = BLID;
      else Acol2Irow[i] = Bview.importMatrix->getRowMap()->getLocalElement(GID);
    }
  }

  //To form C = A*B we're going to execute this expression:
  //
  // C(i,j) = sum_k( A(i,k)*B(k,j) )
  //
  //Our goal, of course, is to navigate the data in A and B once, without
  //performing searches for column-indices, etc.
  ArrayRCP<const size_t> Arowptr_RCP, Browptr_RCP, Irowptr_RCP;
  ArrayRCP<const LocalOrdinal> Acolind_RCP, Bcolind_RCP, Icolind_RCP;
  ArrayRCP<const Scalar> Avals_RCP, Bvals_RCP, Ivals_RCP;
  ArrayView<const size_t> Arowptr, Browptr, Irowptr;
  ArrayView<const LocalOrdinal> Acolind, Bcolind, Icolind;
  ArrayView<const Scalar> Avals, Bvals, Ivals;
  Aview.origMatrix->getAllValues(Arowptr_RCP,Acolind_RCP,Avals_RCP);
  Bview.origMatrix->getAllValues(Browptr_RCP,Bcolind_RCP,Bvals_RCP);
  Arowptr = Arowptr_RCP();  Acolind = Acolind_RCP();  Avals = Avals_RCP();
  Browptr = Browptr_RCP();  Bcolind = Bcolind_RCP();  Bvals = Bvals_RCP();
  if(!Bview.importMatrix.is_null()) {
    Bview.importMatrix->getAllValues(Irowptr_RCP,Icolind_RCP,Ivals_RCP);
    Irowptr = Irowptr_RCP();  Icolind = Icolind_RCP();  Ivals = Ivals_RCP();
  }

  bool C_filled = C.isFillComplete();

  for (size_t i = 0; i < C_numCols; i++)
      c_index[i] = OrdinalTraits<size_t>::invalid();

  //loop over the rows of A.
  size_t Arows = Aview.rowMap->getNodeNumElements();
  for(size_t i=0; i<Arows; ++i) {

    //only navigate the local portion of Aview... which is, thankfully, all of A
    //since this routine doesn't do transpose modes
    GlobalOrdinal global_row = Aview.rowMap->getGlobalElement(i);

    //loop across the i-th row of A and for each corresponding row
    //in B, loop across colums and accumulate product
    //A(i,k)*B(k,j) into our partial sum quantities C_row_i. In other words,
    //as we stride across B(k,:) we're calculating updates for row i of the
    //result matrix C.


    C_row_i_length = OrdinalTraits<size_t>::zero();

    for(k = Arowptr[i]; k < Arowptr[i+1]; ++k) {
      LocalOrdinal Ak = Acol2Brow[Acolind[k]];
      Scalar Aval = Avals[k];
      if (Aval == STS::zero())
        continue;

      if (Ak==LO_INVALID) continue;

      for(j=Browptr[Ak]; j< Browptr[Ak+1]; ++j) {
          LocalOrdinal col = Bcolind[j];
          //assert(col >= 0 && col < C_numCols);
          if (c_index[col] == OrdinalTraits<size_t>::invalid()){
          //assert(C_row_i_length >= 0 && C_row_i_length < C_numCols);
            // This has to be a +=  so insertGlobalValue goes out
            C_row_i[C_row_i_length] = Aval*Bvals[j];
            C_cols[C_row_i_length] = col;
            c_index[col] = C_row_i_length;
            C_row_i_length++;
          }
          else {
            C_row_i[c_index[col]] += Aval*Bvals[j];
          }
        }
    }

    for (size_t ii = 0; ii < C_row_i_length; ii++) {
      c_index[C_cols[ii]] = OrdinalTraits<size_t>::invalid();
      C_cols[ii] = bcols[C_cols[ii]];
      combined_index[ii] = C_cols[ii];
      combined_values[ii] = C_row_i[ii];
    }
    last_index = C_row_i_length;

    //
    //Now put the C_row_i values into C.
    //
    // We might have to revamp this later.
    C_row_i_length = OrdinalTraits<size_t>::zero();

    for(k = Arowptr[i]; k < Arowptr[i+1]; ++k) {
      LocalOrdinal Ak = Acol2Brow[Acolind[k]];
      Scalar Aval = Avals[k];
      if (Aval == STS::zero())
        continue;

      if (Ak!=LO_INVALID) continue;

      Ak = Acol2Irow[Acolind[k]];
      for(j=Irowptr[Ak]; j< Irowptr[Ak+1]; ++j) {
          LocalOrdinal col = Icolind[j];
          //assert(col >= 0 && col < C_numCols);
          if (c_index[col] == OrdinalTraits<size_t>::invalid()){
          //assert(C_row_i_length >= 0 && C_row_i_length < C_numCols);
            // This has to be a +=  so insertGlobalValue goes out
            C_row_i[C_row_i_length] = Aval*Ivals[j];
            C_cols[C_row_i_length] = col;
            c_index[col] = C_row_i_length;
            C_row_i_length++;
            }
            else {
              // This has to be a +=  so insertGlobalValue goes out
              C_row_i[c_index[col]] += Aval*Ivals[j];
            }
        }
    }

    for (size_t ii = 0; ii < C_row_i_length; ii++) {
      c_index[C_cols[ii]] = OrdinalTraits<size_t>::invalid();
      C_cols[ii] = bcols_import[C_cols[ii]];
      combined_index[last_index] = C_cols[ii];
      combined_values[last_index] = C_row_i[ii];
      last_index++;
    }

      //
      //Now put the C_row_i values into C.
      //
      // We might have to revamp this later.
    C_filled ?
      C.sumIntoGlobalValues(
          global_row,
          combined_index.view(OrdinalTraits<size_t>::zero(), last_index),
          combined_values.view(OrdinalTraits<size_t>::zero(), last_index))
      :
      C.insertGlobalValues(
          global_row,
          combined_index.view(OrdinalTraits<size_t>::zero(), last_index),
          combined_values.view(OrdinalTraits<size_t>::zero(), last_index));

  }

}

template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
void setMaxNumEntriesPerRow(
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Mview)
{
  typedef typename Array<ArrayView<const LocalOrdinal> >::size_type  local_length_size;
  Mview.maxNumRowEntries = OrdinalTraits<local_length_size>::zero();
  if(Mview.indices.size() > OrdinalTraits<local_length_size>::zero() ){
    Mview.maxNumRowEntries = Mview.indices[0].size();
    for(local_length_size i = 1; i<Mview.indices.size(); ++i){
      if(Mview.indices[i].size() > Mview.maxNumRowEntries){
        Mview.maxNumRowEntries = Mview.indices[i].size();
      }
    }
  }
}


template<class CrsMatrixType>
size_t C_estimate_nnz(CrsMatrixType & A, CrsMatrixType &B){
  // Follows the NZ estimate in ML's ml_matmatmult.c
  size_t Aest = 100, Best=100;
  if(A.getNodeNumEntries() > 0)
    Aest = (A.getNodeNumRows()>0)? A.getNodeNumEntries()/A.getNodeNumEntries():100;
  if(B.getNodeNumEntries() > 0)
    Best=(B.getNodeNumRows()>0)? B.getNodeNumEntries()/B.getNodeNumEntries():100;

  size_t nnzperrow=(size_t)(sqrt((double)Aest) + sqrt((double)Best) - 1);
  nnzperrow*=nnzperrow;

  return (size_t)(A.getNodeNumRows()*nnzperrow*0.75 + 100);
}



//kernel method for computing the local portion of C = A*B
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
void mult_A_B_newmatrix(
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
  const std::string & label)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  typedef Import<LocalOrdinal, GlobalOrdinal, Node> import_type;
  typedef Map<LocalOrdinal, GlobalOrdinal, Node> map_type;

#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix = std::string("TpetraExt ")+ label + std::string(": ");
  using Teuchos::TimeMonitor;
  RCP<TimeMonitor> MM = rcp(new TimeMonitor(*(TimeMonitor::getNewTimer(prefix+std::string("MMM M5 Cmap")))));
#endif
  size_t ST_INVALID = Teuchos::OrdinalTraits<LocalOrdinal>::invalid();
  LocalOrdinal LO_INVALID = Teuchos::OrdinalTraits<LocalOrdinal>::invalid();


  // Build the final importer / column map, hash table lookups for C
  RCP<const import_type> Cimport;
  RCP<const map_type> Ccolmap;
  RCP<const import_type> Bimport = Bview.origMatrix->getGraph()->getImporter();
  RCP<const import_type> Iimport = Bview.importMatrix.is_null() ? Teuchos::null :  Bview.importMatrix->getGraph()->getImporter();
  Array<LocalOrdinal> Bcol2Ccol(Bview.colMap->getNodeNumElements()), Icol2Ccol;

  if(Bview.importMatrix.is_null()) {
    Cimport = Bimport;
    Ccolmap = Bview.colMap;
    // Bcol2Ccol is trivial
    for(size_t i=0; i<Bview.colMap->getNodeNumElements(); i++) {
      Bcol2Ccol[i] = Teuchos::as<LocalOrdinal>(i);
    }
  }
  else {
    // Choose the right variant of setUnion
    if(!Bimport.is_null() && !Iimport.is_null()){
      Cimport = Bimport->setUnion(*Iimport);
      Ccolmap = Cimport->getTargetMap();
    }
    else if(!Bimport.is_null() && Iimport.is_null()) {
      Cimport = Bimport->setUnion();
    }
    else if(Bimport.is_null() && !Iimport.is_null()) {
      Cimport = Iimport->setUnion();
    }
    else
      throw std::runtime_error("TpetraExt::MMM status of matrix importers is nonsensical");

    Ccolmap = Cimport->getTargetMap();

    if(!Cimport->getSourceMap()->isSameAs(*Bview.origMatrix->getDomainMap()))
      throw std::runtime_error("Tpetra::MMM: Import setUnion messed with the DomainMap in an unfortunate way");

    // NOTE: This is not efficient and should be folded into setUnion
    Icol2Ccol.resize(Bview.importMatrix->getColMap()->getNodeNumElements());
    ArrayView<const GlobalOrdinal> Bgid = Bview.origMatrix->getColMap()->getNodeElementList();
    ArrayView<const GlobalOrdinal> Igid = Bview.importMatrix->getColMap()->getNodeElementList();

    for(size_t i=0; i<Bview.origMatrix->getColMap()->getNodeNumElements(); i++)
      Bcol2Ccol[i] = Ccolmap->getLocalElement(Bgid[i]);
    for(size_t i=0; i<Bview.importMatrix->getColMap()->getNodeNumElements(); i++)
      Icol2Ccol[i] = Ccolmap->getLocalElement(Igid[i]);
  }

#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix+std::string("MMM Newmatrix SerialCore"))));
#endif

  // Sizes
  size_t m=Aview.origMatrix->getNodeNumRows();
  size_t n=Ccolmap->getNodeNumElements();

  // Get Data Pointers
  ArrayRCP<const size_t> Arowptr_RCP, Browptr_RCP, Irowptr_RCP;
  ArrayRCP<size_t> Crowptr_RCP;
  ArrayRCP<const LocalOrdinal> Acolind_RCP, Bcolind_RCP, Icolind_RCP;
  ArrayRCP<LocalOrdinal> Ccolind_RCP;
  ArrayRCP<const Scalar> Avals_RCP, Bvals_RCP, Ivals_RCP;
  ArrayRCP<Scalar> Cvals_RCP;

  Aview.origMatrix->getAllValues(Arowptr_RCP,Acolind_RCP,Avals_RCP);
  Bview.origMatrix->getAllValues(Browptr_RCP,Bcolind_RCP,Bvals_RCP);
  if(!Bview.importMatrix.is_null()) Bview.importMatrix->getAllValues(Irowptr_RCP,Icolind_RCP,Ivals_RCP);


  // For efficiency
  ArrayView<const size_t> Arowptr, Browptr, Irowptr;
  ArrayView<const LocalOrdinal> Acolind, Bcolind, Icolind;
  ArrayView<const Scalar> Avals, Bvals, Ivals;
  ArrayView<size_t> Crowptr;
  ArrayView<LocalOrdinal> Ccolind;
  ArrayView<Scalar> Cvals;
  Arowptr = Arowptr_RCP();  Acolind = Acolind_RCP();  Avals = Avals_RCP();
  Browptr = Browptr_RCP();  Bcolind = Bcolind_RCP();  Bvals = Bvals_RCP();
  if(!Bview.importMatrix.is_null()) {
    Irowptr = Irowptr_RCP();  Icolind = Icolind_RCP();  Ivals = Ivals_RCP();
  }

  // The status array will contain the index into colind where this entry was last deposited.
  // c_status[i] < CSR_ip - not in the row yet.
  // c_status[i] >= CSR_ip, this is the entry where you can find the data
  // We start with this filled with INVALID's indicating that there are no entries yet.
  // Sadly, this complicates the code due to the fact that size_t's are unsigned.
  size_t INVALID = Teuchos::OrdinalTraits<size_t>::invalid();
  Array<size_t> c_status(n, ST_INVALID);

  // Classic csr assembly (low memory edition)
  size_t CSR_alloc=std::max(C_estimate_nnz(*Aview.origMatrix,*Bview.origMatrix),n);
  size_t CSR_ip=0,OLD_ip=0;
  Crowptr_RCP.resize(m+1);       Crowptr = Crowptr_RCP();
  Ccolind_RCP.resize(CSR_alloc); Ccolind = Ccolind_RCP();
  Cvals_RCP.resize(CSR_alloc);   Cvals   = Cvals_RCP();

  // Run through all the hash table lookups once and for all
  Array<LocalOrdinal> targetMapToOrigRow(Aview.colMap->getNodeNumElements(),LO_INVALID);
  Array<LocalOrdinal> targetMapToImportRow(Aview.colMap->getNodeNumElements(),LO_INVALID);

  if(Aview.colMap->isSameAs(*Bview.rowMap)){
    // Maps are the same: Use local IDs as the hash
    for(LocalOrdinal i=Aview.colMap->getMinLocalIndex(); i <= Aview.colMap->getMaxLocalIndex(); i++) {
      LocalOrdinal B_LID = Bview.origMatrix->getRowMap()->getLocalElement(Aview.colMap->getGlobalElement(i));
      if(B_LID != LO_INVALID) targetMapToOrigRow[i] = B_LID;
      else {
        LocalOrdinal I_LID = Bview.importMatrix->getRowMap()->getLocalElement(Aview.colMap->getGlobalElement(i));
        targetMapToImportRow[i] = I_LID;
      }
    }
  }
  else {
    // Maps are not the same:  Use the map's hash
    for(LocalOrdinal i=Aview.colMap->getMinLocalIndex(); i <= Aview.colMap->getMaxLocalIndex(); i++) {
      LocalOrdinal B_LID = Bview.origMatrix->getRowMap()->getLocalElement(Aview.colMap->getGlobalElement(i));
      if(B_LID != LO_INVALID) targetMapToOrigRow[i] = B_LID;
      else {
        LocalOrdinal I_LID = Bview.importMatrix->getRowMap()->getLocalElement(Aview.colMap->getGlobalElement(i));
        targetMapToImportRow[i] = I_LID;
      }
    }
  }

  const Scalar SC_ZERO = Teuchos::ScalarTraits<Scalar>::zero();

  // For each row of A/C
  for(size_t i=0; i<m; i++){
    Crowptr[i]=CSR_ip;

    for(size_t k=Arowptr[i]; k<Arowptr[i+1]; k++){
      LocalOrdinal Ak      = Acolind[k];
      Scalar       Aval    = Avals[k];
      if(Aval==SC_ZERO) continue;

      if(targetMapToOrigRow[Ak] != LO_INVALID){
        // Local matrix
        size_t Bk = Teuchos::as<size_t>(targetMapToOrigRow[Ak]);

        for(size_t j=Browptr[Bk]; j<Browptr[Bk+1]; ++j) {
          LocalOrdinal Cj=Bcol2Ccol[Bcolind[j]];

          if(c_status[Cj]==INVALID || c_status[Cj]<OLD_ip){
            // New entry
            c_status[Cj]      = CSR_ip;
            Ccolind[CSR_ip]= Cj;
            Cvals[CSR_ip]  = Aval*Bvals[j];
            CSR_ip++;
          }
          else
            Cvals[c_status[Cj]]+=Aval*Bvals[j];
        }
      }
      else{
        // Remote matrix
        size_t Ik = Teuchos::as<size_t>(targetMapToImportRow[Ak]);
        for(size_t j=Irowptr[Ik]; j<Irowptr[Ik+1]; ++j) {
          LocalOrdinal Cj=Icol2Ccol[Icolind[j]];

          if(c_status[Cj]==INVALID || c_status[Cj]<OLD_ip){
            // New entry
            c_status[Cj]=CSR_ip;
            Ccolind[CSR_ip]=Cj;
            Cvals[CSR_ip]=Aval*Ivals[j];
            CSR_ip++;
          }
          else
            Cvals[c_status[Cj]]+=Aval*Ivals[j];
        }
      }
    }

    // Resize for next pass if needed
    if(CSR_ip + n > CSR_alloc){
      CSR_alloc*=2;
      Ccolind_RCP.resize(CSR_alloc); Ccolind = Ccolind_RCP();
      Cvals_RCP.resize(CSR_alloc);   Cvals   = Cvals_RCP();
    }
    OLD_ip=CSR_ip;
  }

  Crowptr[m]=CSR_ip;

  // Downward resize
  Cvals_RCP.resize(CSR_ip);
  Ccolind_RCP.resize(CSR_ip);


#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM = rcp (new TimeMonitor (* (TimeMonitor::getNewTimer(prefix+std::string("MMM Newmatrix Final Sort")))));
#endif

  // Replace the column map
  C.replaceColMap(Ccolmap);

  // Final sort & set of CRS arrays
  Import_Util::sortCrsEntries(Crowptr_RCP(),Ccolind_RCP(),Cvals_RCP());
  C.setAllValues(Crowptr_RCP,Ccolind_RCP,Cvals_RCP);

#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM = rcp (new TimeMonitor (* (TimeMonitor::getNewTimer(prefix+std::string("MMM Newmatrix ESFC")))));
#endif

  // Final FillComplete
  C.expertStaticFillComplete(Bview.origMatrix->getDomainMap(),Aview.origMatrix->getRangeMap(),Cimport);
}


//kernel method for computing the local portion of C = (I-omega D^{-1} A)*B
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
void jacobi_A_B_newmatrix(
  Scalar omega,
  const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> & Dinv,
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
  const std::string & label)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayView;
  typedef Import<LocalOrdinal, GlobalOrdinal, Node> import_type;
  typedef Map<LocalOrdinal, GlobalOrdinal, Node> map_type;

#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix = std::string("TpetraExt ")+ label + std::string(": ");
  using Teuchos::TimeMonitor;
  RCP<TimeMonitor> MM = rcp(new TimeMonitor(*(TimeMonitor::getNewTimer(prefix+std::string("Jacobi M5 Cmap")))));
#endif
  size_t ST_INVALID = Teuchos::OrdinalTraits<LocalOrdinal>::invalid();
  LocalOrdinal LO_INVALID = Teuchos::OrdinalTraits<LocalOrdinal>::invalid();


  // Build the final importer / column map, hash table lookups for C
  RCP<const import_type> Cimport;
  RCP<const map_type> Ccolmap;
  RCP<const import_type> Bimport = Bview.origMatrix->getGraph()->getImporter();
  RCP<const import_type> Iimport = Bview.importMatrix.is_null() ? Teuchos::null :  Bview.importMatrix->getGraph()->getImporter();
  Array<LocalOrdinal> Bcol2Ccol(Bview.colMap->getNodeNumElements()), Icol2Ccol;

  if(Bview.importMatrix.is_null()) {
    Cimport = Bimport;
    Ccolmap = Bview.colMap;
    // Bcol2Ccol is trivial
    for(size_t i=0; i<Bview.colMap->getNodeNumElements(); i++) {
      Bcol2Ccol[i] = Teuchos::as<LocalOrdinal>(i);
    }
  }
  else {
    // Choose the right variant of setUnion
    if(!Bimport.is_null() && !Iimport.is_null()){
      Cimport = Bimport->setUnion(*Iimport);
      Ccolmap = Cimport->getTargetMap();
    }
    else if(!Bimport.is_null() && Iimport.is_null()) {
      Cimport = Bimport->setUnion();
    }
    else if(Bimport.is_null() && !Iimport.is_null()) {
      Cimport = Iimport->setUnion();
    }
    else
      throw std::runtime_error("TpetraExt::Jacobi status of matrix importers is nonsensical");

    Ccolmap = Cimport->getTargetMap();

    if(!Cimport->getSourceMap()->isSameAs(*Bview.origMatrix->getDomainMap()))
      throw std::runtime_error("Tpetra:Jacobi Import setUnion messed with the DomainMap in an unfortunate way");

    // NOTE: This is not efficient and should be folded into setUnion
    Icol2Ccol.resize(Bview.importMatrix->getColMap()->getNodeNumElements());
    ArrayView<const GlobalOrdinal> Bgid = Bview.origMatrix->getColMap()->getNodeElementList();
    ArrayView<const GlobalOrdinal> Igid = Bview.importMatrix->getColMap()->getNodeElementList();

    for(size_t i=0; i<Bview.origMatrix->getColMap()->getNodeNumElements(); i++)
      Bcol2Ccol[i] = Ccolmap->getLocalElement(Bgid[i]);
    for(size_t i=0; i<Bview.importMatrix->getColMap()->getNodeNumElements(); i++)
      Icol2Ccol[i] = Ccolmap->getLocalElement(Igid[i]);
  }

#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix+std::string("Jacobi Newmatrix SerialCore"))));
#endif

  // Sizes
  size_t m=Aview.origMatrix->getNodeNumRows();
  size_t n=Ccolmap->getNodeNumElements();

  // Get Data Pointers
  ArrayRCP<const size_t> Arowptr_RCP, Browptr_RCP, Irowptr_RCP;
  ArrayRCP<size_t> Crowptr_RCP;
  ArrayRCP<const LocalOrdinal> Acolind_RCP, Bcolind_RCP, Icolind_RCP;
  ArrayRCP<LocalOrdinal> Ccolind_RCP;
  ArrayRCP<const Scalar> Avals_RCP, Bvals_RCP, Ivals_RCP;
  ArrayRCP<Scalar> Cvals_RCP;
  ArrayRCP<const Scalar> Dvals_RCP;

  Aview.origMatrix->getAllValues(Arowptr_RCP,Acolind_RCP,Avals_RCP);
  Bview.origMatrix->getAllValues(Browptr_RCP,Bcolind_RCP,Bvals_RCP);
  if(!Bview.importMatrix.is_null()) Bview.importMatrix->getAllValues(Irowptr_RCP,Icolind_RCP,Ivals_RCP);
  Dvals_RCP = Dinv.getData();

  // For efficiency
  ArrayView<const size_t> Arowptr, Browptr, Irowptr;
  ArrayView<const LocalOrdinal> Acolind, Bcolind, Icolind;
  ArrayView<const Scalar> Avals, Bvals, Ivals;
  ArrayView<size_t> Crowptr;
  ArrayView<LocalOrdinal> Ccolind;
  ArrayView<Scalar> Cvals;
  ArrayView<const Scalar> Dvals;
  Arowptr = Arowptr_RCP();  Acolind = Acolind_RCP();  Avals = Avals_RCP();
  Browptr = Browptr_RCP();  Bcolind = Bcolind_RCP();  Bvals = Bvals_RCP();
  if(!Bview.importMatrix.is_null()) {
    Irowptr = Irowptr_RCP();  Icolind = Icolind_RCP();  Ivals = Ivals_RCP();
  }
  Dvals = Dvals_RCP();

  // The status array will contain the index into colind where this entry was last deposited.
  // c_status[i] < CSR_ip - not in the row yet.
  // c_status[i] >= CSR_ip, this is the entry where you can find the data
  // We start with this filled with INVALID's indicating that there are no entries yet.
  // Sadly, this complicates the code due to the fact that size_t's are unsigned.
  size_t INVALID = Teuchos::OrdinalTraits<size_t>::invalid();
  Array<size_t> c_status(n, ST_INVALID);

  // Classic csr assembly (low memory edition)
  size_t CSR_alloc=std::max(C_estimate_nnz(*Aview.origMatrix,*Bview.origMatrix),n);
  size_t CSR_ip=0,OLD_ip=0;
  Crowptr_RCP.resize(m+1);       Crowptr = Crowptr_RCP();
  Ccolind_RCP.resize(CSR_alloc); Ccolind = Ccolind_RCP();
  Cvals_RCP.resize(CSR_alloc);   Cvals   = Cvals_RCP();

  // Run through all the hash table lookups once and for all
  Array<LocalOrdinal> targetMapToOrigRow(Aview.colMap->getNodeNumElements(),LO_INVALID);
  Array<LocalOrdinal> targetMapToImportRow(Aview.colMap->getNodeNumElements(),LO_INVALID);

  if(Aview.colMap->isSameAs(*Bview.rowMap)){
    // Maps are the same: Use local IDs as the hash
    for(LocalOrdinal i=Aview.colMap->getMinLocalIndex(); i <= Aview.colMap->getMaxLocalIndex(); i++) {
      LocalOrdinal B_LID = Bview.origMatrix->getRowMap()->getLocalElement(Aview.colMap->getGlobalElement(i));
      if(B_LID != LO_INVALID) targetMapToOrigRow[i] = B_LID;
      else {
        LocalOrdinal I_LID = Bview.importMatrix->getRowMap()->getLocalElement(Aview.colMap->getGlobalElement(i));
        targetMapToImportRow[i] = I_LID;
      }
    }
  }
  else {
    // Maps are not the same:  Use the map's hash
    for(LocalOrdinal i=Aview.colMap->getMinLocalIndex(); i <= Aview.colMap->getMaxLocalIndex(); i++) {
      LocalOrdinal B_LID = Bview.origMatrix->getRowMap()->getLocalElement(Aview.colMap->getGlobalElement(i));
      if(B_LID != LO_INVALID) targetMapToOrigRow[i] = B_LID;
      else {
        LocalOrdinal I_LID = Bview.importMatrix->getRowMap()->getLocalElement(Aview.colMap->getGlobalElement(i));
        targetMapToImportRow[i] = I_LID;
      }
    }
  }

  const Scalar SC_ZERO = Teuchos::ScalarTraits<Scalar>::zero();

  // For each row of A/C
  for(size_t i=0; i<m; i++){
    Crowptr[i]=CSR_ip;
    Scalar Dval = Dvals[i];

    // Entries of B
    for(size_t k=Browptr[i]; k<Browptr[i+1]; k++){
      Scalar Bval = Bvals[k];
      if(Bval==SC_ZERO) continue;
      LocalOrdinal Ck=Bcol2Ccol[Bcolind[k]];

      // Assume no repeated entries in B
      c_status[Ck]    = CSR_ip;
      Ccolind[CSR_ip] = Ck;
      Cvals[CSR_ip]   = Bvals[k];
      CSR_ip++;
    }


    // Entries of -omega * Dinv * A * B
    for(size_t k=Arowptr[i]; k<Arowptr[i+1]; k++){
      LocalOrdinal Ak      = Acolind[k];
      Scalar       Aval    = Avals[k];
      if(Aval==SC_ZERO) continue;

      if(targetMapToOrigRow[Ak] != LO_INVALID){
        // Local matrix
        size_t Bk = Teuchos::as<size_t>(targetMapToOrigRow[Ak]);

        for(size_t j=Browptr[Bk]; j<Browptr[Bk+1]; ++j) {
          LocalOrdinal Cj=Bcol2Ccol[Bcolind[j]];

          if(c_status[Cj]==INVALID || c_status[Cj]<OLD_ip){
            // New entry
            c_status[Cj]    = CSR_ip;
            Ccolind[CSR_ip] = Cj;
            Cvals[CSR_ip]   = - omega * Dval* Aval * Bvals[j];
            CSR_ip++;
          }
          else
            Cvals[c_status[Cj]] -= omega * Dval* Aval * Bvals[j];
        }
      }
      else{
        // Remote matrix
        size_t Ik = Teuchos::as<size_t>(targetMapToImportRow[Ak]);
        for(size_t j=Irowptr[Ik]; j<Irowptr[Ik+1]; ++j) {
          LocalOrdinal Cj=Icol2Ccol[Icolind[j]];

          if(c_status[Cj]==INVALID || c_status[Cj]<OLD_ip){
            // New entry
            c_status[Cj]    = CSR_ip;
            Ccolind[CSR_ip] = Cj;
            Cvals[CSR_ip]   = - omega * Dval* Aval * Ivals[j];
            CSR_ip++;
          }
          else
            Cvals[c_status[Cj]] -= omega * Dval* Aval * Ivals[j];
        }
      }
    }

    // Resize for next pass if needed
    if(CSR_ip + n > CSR_alloc){
      CSR_alloc*=2;
      Ccolind_RCP.resize(CSR_alloc); Ccolind = Ccolind_RCP();
      Cvals_RCP.resize(CSR_alloc);   Cvals   = Cvals_RCP();
    }
    OLD_ip=CSR_ip;
  }

  Crowptr[m]=CSR_ip;

  // Downward resize
  Cvals_RCP.resize(CSR_ip);
  Ccolind_RCP.resize(CSR_ip);


#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM = rcp (new TimeMonitor (* (TimeMonitor::getNewTimer(prefix+std::string("Jacobi Newmatrix Final Sort")))));
#endif

  // Replace the column map
  C.replaceColMap(Ccolmap);

  // Final sort & set of CRS arrays
  Import_Util::sortCrsEntries(Crowptr_RCP(),Ccolind_RCP(),Cvals_RCP());
  C.setAllValues(Crowptr_RCP,Ccolind_RCP,Cvals_RCP);

#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM = rcp (new TimeMonitor (* (TimeMonitor::getNewTimer(prefix+std::string("Jacobi Newmatrix ESFC")))));
#endif

  // Final FillComplete
  C.expertStaticFillComplete(Bview.origMatrix->getDomainMap(),Aview.origMatrix->getRangeMap(),Cimport);
}


template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
void import_and_extract_views(
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& M,
  RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > targetMap,
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Mview,
  RCP<const Import<LocalOrdinal, GlobalOrdinal, Node> > prototypeImporter,
  bool userAssertsThereAreNoRemotes,
  const std::string & label)
{
#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix = std::string("TpetraExt ")+ label + std::string(": ");
  using Teuchos::TimeMonitor;
  Teuchos::RCP<Teuchos::TimeMonitor> MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix+std::string("MMM I&X Alloc"))));
#endif

  //Convience typedef
  typedef Map<LocalOrdinal, GlobalOrdinal, Node> Map_t;
  typedef CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> CrsMatrix_t;
  // The goal of this method is to populate the 'Mview' struct with views of the
  // rows of M, including all rows that correspond to elements in 'targetMap'.
  //
  // If targetMap includes local elements that correspond to remotely-owned rows
  // of M, then those remotely-owned rows will be imported into
  // 'Mview.importMatrix', and views of them will be included in 'Mview'.
  Mview.deleteContents();

  RCP<const Map_t> Mrowmap = M.getRowMap();
  RCP<const Map_t> MremoteRowMap;
  const int numProcs = Mrowmap->getComm()->getSize();

  ArrayView<const GlobalOrdinal> Mrows = targetMap->getNodeElementList();

  size_t numRemote = 0;
  size_t numRows   = targetMap->getNodeNumElements();
  Mview.origMatrix = Teuchos::rcp(&M,false);
  Mview.origRowMap = M.getRowMap();
  Mview.rowMap = targetMap;
  Mview.colMap = M.getColMap();
  Mview.domainMap = M.getDomainMap();
  Mview.importColMap = null;

  // Short circuit if the user swears there are no remotes
  if(userAssertsThereAreNoRemotes) return;

#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix+std::string("MMM I&X RemoteMap"))));
#endif

  // mark each row in targetMap as local or remote, and go ahead and get a view for the local rows
  int mode = 0;
  if(!prototypeImporter.is_null() && prototypeImporter->getSourceMap()->isSameAs(*Mrowmap) && prototypeImporter->getTargetMap()->isSameAs(*targetMap)) {
    // We have a valid prototype importer --- ask it for the remotes
    numRemote = prototypeImporter->getNumRemoteIDs();
    Array<GlobalOrdinal> MremoteRows(numRemote);
    ArrayView<const LocalOrdinal> RemoteLIDs = prototypeImporter->getRemoteLIDs();
    for(size_t i=0; i<numRemote; i++) {
      MremoteRows[i] = targetMap->getGlobalElement(RemoteLIDs[i]);
    }

    MremoteRowMap=rcp(new Map_t(OrdinalTraits<global_size_t>::invalid(), MremoteRows(), Mrowmap->getIndexBase(), Mrowmap->getComm(), Mrowmap->getNode()));
    mode=1;
  }
  else if(prototypeImporter.is_null()) {
    // No prototype importer --- count the remotes the hard way
    Array<GlobalOrdinal> MremoteRows(numRows);
    for(size_t i=0; i < numRows; ++i) {
      const LocalOrdinal mlid = Mrowmap->getLocalElement(Mrows[i]);

      if (mlid == OrdinalTraits<LocalOrdinal>::invalid()) {
        MremoteRows[numRemote]=Mrows[i];
        ++numRemote;
      }
    }
    MremoteRows.resize(numRemote);
    MremoteRowMap=rcp(new Map_t(OrdinalTraits<global_size_t>::invalid(), MremoteRows(), Mrowmap->getIndexBase(), Mrowmap->getComm(), Mrowmap->getNode()));
    mode=2;
  }
  else {
    // prototypeImporter is bad.  But if we're in serial that's OK.
    mode=3;
  }

  if (numProcs < 2) {
    TEUCHOS_TEST_FOR_EXCEPTION(numRemote > 0, std::runtime_error,
      "MatrixMatrix::import_and_extract_views ERROR, numProcs < 2 but attempting to import remote matrix rows." <<std::endl);
    //If only one processor we don't need to import any remote rows, so return.
    return;
  }

  //
  // Now we will import the needed remote rows of M, if the global maximum
  // value of numRemote is greater than 0.
  //
#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix+std::string("MMM I&X Collective-0"))));
#endif

  global_size_t globalMaxNumRemote = 0;
  Teuchos::reduceAll(*(Mrowmap->getComm()) , Teuchos::REDUCE_MAX, (global_size_t)numRemote, Teuchos::outArg(globalMaxNumRemote) );


  if (globalMaxNumRemote > 0) {
#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix+std::string("MMM I&X Import-2"))));
#endif
    // Create an importer with target-map MremoteRowMap and source-map Mrowmap.
    RCP<const Import<LocalOrdinal, GlobalOrdinal, Node> > importer;

    if(mode==1)
      importer = prototypeImporter->createRemoteOnlyImport(MremoteRowMap);
    else if(mode==2)
      importer=rcp(new Import<LocalOrdinal, GlobalOrdinal, Node>(Mrowmap, MremoteRowMap));
    else
      throw std::runtime_error("prototypeImporter->SourceMap() does not match M.getRowMap()!");

#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix+std::string("MMM I&X Import-3"))));
#endif

    // Now create a new matrix into which we can import the remote rows of M that we need.
    Mview.importMatrix = Tpetra::importAndFillCompleteCrsMatrix<CrsMatrix_t>(Teuchos::rcp(&M,false),*importer,M.getDomainMap(),M.getRangeMap(),Teuchos::null);

#ifdef COMPUTE_MMM_STATISTICS
    printMultiplicationStatistics(importer,label+std::string(" I&X MMM"));
#endif


#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix+std::string("MMM I&X Import-4"))));
#endif

    // Save the column map of the imported matrix, so that we can convert indices back to global for arithmetic later
    Mview.importColMap = Mview.importMatrix->getColMap();
  }
}



} //End namepsace MMdetails

} //End namespace Tpetra
//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_MATRIXMATRIX_INSTANT(SCALAR,LO,GO,NODE) \
  \
  template \
  void MatrixMatrix::Multiply( \
    const CrsMatrix< SCALAR , LO , GO , NODE >& A, \
    bool transposeA, \
    const CrsMatrix< SCALAR , LO , GO , NODE >& B, \
    bool transposeB, \
    CrsMatrix< SCALAR , LO , GO , NODE >& C, \
    bool call_FillComplete_on_result, \
    const std::string & label); \
\
template \
  void MatrixMatrix::Jacobi( \
    SCALAR omega, \
    const Vector< SCALAR, LO, GO, NODE > & Dinv, \
    const CrsMatrix< SCALAR , LO , GO , NODE >& A, \
    const CrsMatrix< SCALAR , LO , GO , NODE >& B, \
    CrsMatrix< SCALAR , LO , GO , NODE >& C, \
    bool call_FillComplete_on_result, \
    const std::string & label); \
\
  template \
  void MatrixMatrix::Add( \
    const CrsMatrix< SCALAR , LO , GO , NODE >& A, \
    bool transposeA, \
    SCALAR scalarA, \
    const CrsMatrix< SCALAR , LO , GO , NODE >& B, \
    bool transposeB, \
    SCALAR scalarB, \
    RCP<CrsMatrix< SCALAR , LO , GO , NODE > > C); \
  \
  template \
  void MatrixMatrix::Add( \
    const CrsMatrix<SCALAR, LO, GO, NODE>& A, \
    bool transposeA, \
    SCALAR scalarA, \
    CrsMatrix<SCALAR, LO, GO, NODE>& B, \
    SCALAR scalarB ); \
  \
  template \
  Teuchos::RCP<CrsMatrix< SCALAR , LO , GO , NODE > > \
  MatrixMatrix::add (const SCALAR & alpha, \
                     const bool transposeA, \
                     const CrsMatrix< SCALAR , LO , GO , NODE >& A, \
                     const SCALAR & beta, \
                     const bool transposeB, \
                     const CrsMatrix< SCALAR , LO , GO , NODE >& B, \
                     const Teuchos::RCP<const Map< LO , GO , NODE > >& domainMap, \
                     const Teuchos::RCP<const Map< LO , GO , NODE > >& rangeMap, \
                     const Teuchos::RCP<Teuchos::ParameterList>& params); \
\


#endif // TPETRA_MATRIXMATRIX_DEF_HPP
