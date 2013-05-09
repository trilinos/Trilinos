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
#include <algorithm>
#include "Teuchos_FancyOStream.hpp"


/*! \file TpetraExt_MatrixMatrix_def.hpp

    The implementations for the members of class Tpetra::MatrixMatrixMultiply and related non-member constructors.
 */

namespace Tpetra {


namespace MatrixMatrix{

template <class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class Node,
           class SpMatOps >
void Multiply(
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& A,
  bool transposeA,
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& B,
  bool transposeB,
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& C,
  bool call_FillComplete_on_result)
{
  //TEUCHOS_FUNC_TIME_MONITOR_DIFF("My Matrix Mult", mmm_multiply);
  typedef CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps> Matrix_t;
  //
  //This method forms the matrix-matrix product C = op(A) * op(B), where
  //op(A) == A   if transposeA is false,
  //op(A) == A^T if transposeA is true,
  //and similarly for op(B).
  //

  //A and B should already be Filled.
  //(Should we go ahead and call FillComplete() on them if necessary?
  // or error out? For now, we choose to error out.)
  TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), std::runtime_error,
    "Uh oh. Looks like there's a bit of a problem here. No worries though. We'll help you figure it out. You're "
    "a fantastic programer and this just a minor bump in the road! Maybe the information below can help you out a bit."
    "\n\n MatrixMatrix::Multiply(): Matrix A is not fill complete.");
  TEUCHOS_TEST_FOR_EXCEPTION(!B.isFillComplete(), std::runtime_error,
    "Uh oh. Looks like there's a bit of a problem here. No worries though. We'll help you figure it out. You're "
    "a fantastic programer and this just a minor bump in the road! Maybe the information below can help you out a bit."
    "\n\n MatrixMatrix::Multiply(): Matrix B is not fill complete.");
  TEUCHOS_TEST_FOR_EXCEPTION(C.isLocallyIndexed() , std::runtime_error,
    "MatrixMatrix::Add ERROR, input matrix C must not be locally indexed!");

  //Convience typedefs
  typedef CrsMatrixStruct<
    Scalar,
    LocalOrdinal,
    GlobalOrdinal,
    Node,
    SpMatOps> CrsMatrixStruct_t;
  typedef Map<LocalOrdinal, GlobalOrdinal, Node> Map_t;

  RCP<const Matrix_t > Aprime = null;
  RCP<const Matrix_t > Bprime = null;
  if(transposeA){
    RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps> at (Teuchos::rcpFromRef (A));
    Aprime = at.createTranspose();
  }
  else{
    Aprime = rcpFromRef(A);
  }
  if(transposeB){
    RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps> bt (Teuchos::rcpFromRef (B));
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

  //Now import any needed remote rows and populate the Aview struct.
  MMdetails::import_and_extract_views(*Aprime, targetMap_A, Aview);


  //We will also need local access to all rows of B that correspond to the
  //column-map of op(A).
  if (numProcs > 1) {
    targetMap_B = Aprime->getColMap(); //colmap_op_A;
  }

  //Now import any needed remote rows and populate the Bview struct.
  MMdetails::import_and_extract_views(*Bprime, targetMap_B, Bview);


  //If the result matrix C is not already FillComplete'd, we will do a
  //preprocessing step to create the nonzero structure,
 /*
  if (!C.isFillComplete()) {
    CrsWrapper_GraphBuilder<Scalar, LocalOrdinal, GlobalOrdinal, Node> crsgraphbuilder(C.getRowMap());

    //pass the graph-builder object to the multiplication kernel to fill in all
    //the nonzero positions that will be used in the result matrix.
    MMdetails::mult_A_B(Aview, Bview, crsgraphbuilder, true);

    //now insert all of the nonzero positions into the result matrix.
    insert_matrix_locations(crsgraphbuilder, C);


    if (call_FillComplete_on_result) {
      C.fillComplete(Bprime->getDomainMap(), Aprime->getRangeMap());
      call_FillComplete_on_result = false;
    }
  }
  */

  //Now call the appropriate method to perform the actual multiplication.

  CrsWrapper_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crsmat(C);

  MMdetails::mult_A_B(Aview, Bview, crsmat);

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

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node,
          class SpMatOps >
void Add(
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& A,
  bool transposeA,
  Scalar scalarA,
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& B,
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
    Node,
    SpMatOps> CrsMatrix_t;
  RCP<const CrsMatrix_t> Aprime = null;
  if( transposeA ){
          RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps> theTransposer(Teuchos::rcpFromRef (A));
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
          class Node,
          class SpMatOps>
void Add(
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& A,
  bool transposeA,
  Scalar scalarA,
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& B,
  bool transposeB,
  Scalar scalarB,
  RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps> > C)
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
  typedef typename ArrayView<const Scalar>::size_type size_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef Map<LocalOrdinal, GlobalOrdinal, Node> map_type;
  typedef Import<LocalOrdinal, GlobalOrdinal, Node> import_type;
  typedef RowGraph<LocalOrdinal, GlobalOrdinal, Node> row_graph_type;
  typedef CrsGraph<LocalOrdinal, GlobalOrdinal, Node, SpMatOps> crs_graph_type;
  typedef CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps> crs_matrix_type;
  typedef RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps> transposer_type;

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


//kernel method for computing the local portion of C = A*B
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node,
         class SpMatOps>
void mult_A_B(
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& Aview,
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& Bview,
  CrsWrapper<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
  bool onlyCalculateStructure)
{
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

  Array<Scalar> dwork = onlyCalculateStructure ? Array<Scalar>() : Array<Scalar>(C_numCols);
  Array<GlobalOrdinal> iwork = Array<GlobalOrdinal>(C_numCols);
  Array<size_t> iwork2 = Array<size_t>(C_numCols);

  Array<Scalar> C_row_i = dwork;
  Array<GlobalOrdinal> C_cols = iwork;
  Array<size_t> c_index = iwork2;
  Array<GlobalOrdinal> combined_index = Array<GlobalOrdinal>(2*C_numCols);
  Array<Scalar> combined_values = onlyCalculateStructure ? Array<Scalar>() : Array<Scalar>(2*C_numCols);

  size_t C_row_i_length, j, k, lastj, last_index;

  // Run through all the hash table lookups once and for all
  Array<LocalOrdinal> Acol2Brow(Aview.colMap->getNodeNumElements());
  if(Aview.colMap->isSameAs(*Bview.rowMap)){
    // Maps are the same: Use local IDs as the hash
    for(LocalOrdinal i=Aview.colMap->getMinLocalIndex(); i <=
            Aview.colMap->getMaxLocalIndex(); i++)
      Acol2Brow[i]=i;
  }
  else {
    // Maps are not the same:  Use the map's hash
    for(LocalOrdinal i=Aview.colMap->getMinLocalIndex(); i <=
            Aview.colMap->getMaxLocalIndex(); i++)
      Acol2Brow[i]=Bview.rowMap->getLocalElement(
            Aview.colMap->getGlobalElement(i));
  }

  //To form C = A*B we're going to execute this expression:
  //
  // C(i,j) = sum_k( A(i,k)*B(k,j) )
  //
  //Our goal, of course, is to navigate the data in A and B once, without
  //performing searches for column-indices, etc.

  bool C_filled = C.isFillComplete();

  for (size_t i = 0; i < C_numCols; i++)
      c_index[i] = OrdinalTraits<size_t>::invalid();

  //loop over the rows of A.
  for(size_t i=0; i<Aview.numRows; ++i) {

    //only navigate the local portion of Aview... (It's probable that we
    //imported more of A than we need for A*B, because other cases like A^T*B
    //need the extra rows.)
    if (Aview.remote[i]) {
      continue;
    }


    ArrayView<const LocalOrdinal> Aindices_i = Aview.indices[i];
    ArrayView<const Scalar> Aval_i  = onlyCalculateStructure ? null : Aview.values[i];

    GlobalOrdinal global_row = Aview.rowMap->getGlobalElement(i);


    //loop across the i-th row of A and for each corresponding row
    //in B, loop across colums and accumulate product
    //A(i,k)*B(k,j) into our partial sum quantities C_row_i. In other words,
    //as we stride across B(k,:) we're calculating updates for row i of the
    //result matrix C.


    C_row_i_length = OrdinalTraits<size_t>::zero();

    for(k = OrdinalTraits<size_t>::zero(); k < Aview.numEntriesPerRow[i]; ++k) {
      LocalOrdinal Ak = Acol2Brow[Aindices_i[k]];
      Scalar Aval = onlyCalculateStructure ? Teuchos::as<Scalar>(0) : Aval_i[k];

      if (Bview.remote[Ak]) continue;

      const LocalOrdinal *Bcol_inds = Bview.indices[Ak].getRawPtr();
      const Scalar *Bvals_k = onlyCalculateStructure ? NULL :
                                     Bview.values[Ak].getRawPtr();

      lastj = Bview.numEntriesPerRow[Ak];

        for(j=OrdinalTraits<size_t>::zero(); j< lastj; ++j) {
          LocalOrdinal col = Bcol_inds[j];
          //assert(col >= 0 && col < C_numCols);
          if (c_index[col] == OrdinalTraits<size_t>::invalid()){
          //assert(C_row_i_length >= 0 && C_row_i_length < C_numCols);
            if(!onlyCalculateStructure){
                // This has to be a +=  so insertGlobalValue goes out
                C_row_i[C_row_i_length] = Aval*Bvals_k[j];
            }
            C_cols[C_row_i_length] = col;
            c_index[col] = C_row_i_length;
            C_row_i_length++;
          }
          else {
            if(!onlyCalculateStructure){
              C_row_i[c_index[col]] += Aval*Bvals_k[j];
            }
          }
        }
    }

    for (size_t ii = 0; ii < C_row_i_length; ii++) {
      c_index[C_cols[ii]] = OrdinalTraits<size_t>::invalid();
      C_cols[ii] = bcols[C_cols[ii]];
      combined_index[ii] = C_cols[ii];
      if (!onlyCalculateStructure)
          combined_values[ii] = C_row_i[ii];
    }
    last_index = C_row_i_length;

      //
      //Now put the C_row_i values into C.
      //
      // We might have to revamp this later.
    /*if (!onlyCalculateStructure)
    {
      C_filled ?
        C.sumIntoGlobalValues(
          global_row,
          C_cols.view(OrdinalTraits<size_t>::zero(), C_row_i_length),
          onlyCalculateStructure ? null :
          C_row_i.view(OrdinalTraits<size_t>::zero(), C_row_i_length))
        :
        C.insertGlobalValues(
          global_row,
          C_cols.view(OrdinalTraits<size_t>::zero(), C_row_i_length),
          onlyCalculateStructure ? null :
          C_row_i.view(OrdinalTraits<size_t>::zero(), C_row_i_length));
    }*/

    C_row_i_length = OrdinalTraits<size_t>::zero();

    for(k = OrdinalTraits<size_t>::zero(); k < Aview.numEntriesPerRow[i]; ++k) {
      LocalOrdinal Ak = Acol2Brow[Aindices_i[k]];
      Scalar Aval = onlyCalculateStructure ? Teuchos::as<Scalar>(0) : Aval_i[k];

      if (!Bview.remote[Ak]) continue;

      const LocalOrdinal *Bcol_inds = Bview.indices[Ak].getRawPtr();
      const Scalar *Bvals_k = onlyCalculateStructure ? NULL :
                                     Bview.values[Ak].getRawPtr();

      lastj = Bview.numEntriesPerRow[Ak];
        for(j=OrdinalTraits<size_t>::zero(); j< lastj; ++j) {
          LocalOrdinal col = Bcol_inds[j];
          //assert(col >= 0 && col < C_numCols);
          if (c_index[col] == OrdinalTraits<size_t>::invalid()){
          //assert(C_row_i_length >= 0 && C_row_i_length < C_numCols);
            if(!onlyCalculateStructure){
                // This has to be a +=  so insertGlobalValue goes out
                C_row_i[C_row_i_length] = Aval*Bvals_k[j];
            }
            C_cols[C_row_i_length] = col;
            c_index[col] = C_row_i_length;
            C_row_i_length++;
            }
            else {
              if(!onlyCalculateStructure){
                // This has to be a +=  so insertGlobalValue goes out
                C_row_i[c_index[col]] += Aval*Bvals_k[j];
              }
            }
        }
    }

    for (size_t ii = 0; ii < C_row_i_length; ii++) {
      c_index[C_cols[ii]] = OrdinalTraits<size_t>::invalid();
      C_cols[ii] = bcols_import[C_cols[ii]];
      combined_index[last_index] = C_cols[ii];
      if (!onlyCalculateStructure)
          combined_values[last_index] = C_row_i[ii];
      last_index++;
    }

      //
      //Now put the C_row_i values into C.
      //
      // We might have to revamp this later.
    if (!onlyCalculateStructure)
    {
      C_filled ?
        C.sumIntoGlobalValues(
          global_row,
          combined_index.view(OrdinalTraits<size_t>::zero(), last_index),
          onlyCalculateStructure ? null : combined_values.view(
          OrdinalTraits<size_t>::zero(), last_index))
        :
        C.insertGlobalValues(
          global_row,
          combined_index.view(OrdinalTraits<size_t>::zero(), last_index),
          onlyCalculateStructure ? null : combined_values.view(
          OrdinalTraits<size_t>::zero(), last_index));
    }

  }

}

template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node,
         class SpMatOps>
void setMaxNumEntriesPerRow(
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& Mview)
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

template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node,
         class SpMatOps>
void import_and_extract_views(
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& M,
  RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > targetMap,
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& Mview)
{
  //Convience typedef
  typedef Map<LocalOrdinal, GlobalOrdinal, Node> Map_t;
  // The goal of this method is to populate the 'Mview' struct with views of the
  // rows of M, including all rows that correspond to elements in 'targetMap'.
  //
  // If targetMap includes local elements that correspond to remotely-owned rows
  // of M, then those remotely-owned rows will be imported into
  // 'Mview.importMatrix', and views of them will be included in 'Mview'.
  Mview.deleteContents();

  RCP<const Map_t> Mrowmap = M.getRowMap();

  const int numProcs = Mrowmap->getComm()->getSize();

  ArrayView<const GlobalOrdinal> Mrows = targetMap->getNodeElementList();

  Mview.numRemote = 0;
  Mview.numRows = targetMap->getNodeNumElements();
  Mview.numEntriesPerRow.resize(Mview.numRows);
  Mview.indices.resize(         Mview.numRows);
  Mview.values.resize(          Mview.numRows);
  Mview.remote.resize(          Mview.numRows);


  Mview.origRowMap = M.getRowMap();
  Mview.rowMap = targetMap;
  Mview.colMap = M.getColMap();
  Mview.domainMap = M.getDomainMap();
  Mview.importColMap = null;


  // mark each row in targetMap as local or remote, and go ahead and get a view for the local rows

  for(size_t i=0; i < Mview.numRows; ++i)
  {
    const LocalOrdinal mlid = Mrowmap->getLocalElement(Mrows[i]);

    if (mlid == OrdinalTraits<LocalOrdinal>::invalid()) {
      Mview.remote[i] = true;
      ++Mview.numRemote;
    }
    else {
      Mview.remote[i] = false;
      M.getLocalRowView(mlid, Mview.indices[i], Mview.values[i]);
            Mview.numEntriesPerRow[i] = Mview.indices[i].size();
    }
  }

  if (numProcs < 2) {
    TEUCHOS_TEST_FOR_EXCEPTION(Mview.numRemote > 0, std::runtime_error,
      "MatrixMatrix::import_and_extract_views ERROR, numProcs < 2 but attempting to import remote matrix rows." <<std::endl);
    setMaxNumEntriesPerRow(Mview);
    //If only one processor we don't need to import any remote rows, so return.
    return;
  }

  //
  // Now we will import the needed remote rows of M, if the global maximum
  // value of numRemote is greater than 0.
  //

  global_size_t globalMaxNumRemote = 0;
  Teuchos::reduceAll(*(Mrowmap->getComm()) , Teuchos::REDUCE_MAX, Mview.numRemote, Teuchos::outArg(globalMaxNumRemote) );

  if (globalMaxNumRemote > 0) {
    // Create a map that describes the remote rows of M that we need.

    Array<GlobalOrdinal> MremoteRows(Mview.numRemote);


    global_size_t offset = 0;
    for(size_t i=0; i < Mview.numRows; ++i) {
      if (Mview.remote[i]) {
        MremoteRows[offset++] = Mrows[i];
      }
    }

    RCP<const Map_t> MremoteRowMap = rcp(new Map_t(
      OrdinalTraits<global_size_t>::invalid(),
      MremoteRows(),
      Mrowmap->getIndexBase(),
      Mrowmap->getComm(),
      Mrowmap->getNode()));

    // Create an importer with target-map MremoteRowMap and source-map Mrowmap.
    Import<LocalOrdinal, GlobalOrdinal, Node> importer(Mrowmap, MremoteRowMap);

    // Now create a new matrix into which we can import the remote rows of M that we need.
    Mview.importMatrix = rcp(new CrsMatrix<Scalar,LocalOrdinal, GlobalOrdinal, Node, SpMatOps>( MremoteRowMap, 1 ));
    Mview.importMatrix->doImport(M, importer, INSERT);
    Mview.importMatrix->fillComplete(M.getDomainMap(), M.getRangeMap());

    // Save the column map of the imported matrix, so that we can convert indices back to global for arithmetic later
    Mview.importColMap = Mview.importMatrix->getColMap();

    // Finally, use the freshly imported data to fill in the gaps in our views of rows of M.
    for(size_t i=0; i < Mview.numRows; ++i)
    {
      if (Mview.remote[i]) {
        const LocalOrdinal importLID = MremoteRowMap->getLocalElement(Mrows[i]);
        Mview.importMatrix->getLocalRowView(importLID,
                                             Mview.indices[i],
                                             Mview.values[i]);
        Mview.numEntriesPerRow[i] = Mview.indices[i].size();
      }
    }
  }
  setMaxNumEntriesPerRow(Mview);
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
    bool call_FillComplete_on_result); \
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
  template  \
  void MatrixMatrix::Add( \
    const CrsMatrix<SCALAR, LO, GO, NODE>& A, \
    bool transposeA, \
    SCALAR scalarA, \
    CrsMatrix<SCALAR, LO, GO, NODE>& B, \
    SCALAR scalarB ); \
  \


#endif // TPETRA_MATRIXMATRIX_DEF_HPP
