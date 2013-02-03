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
    RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps> theTransposer (Teuchos::rcpFromRef (A));
    Aprime = theTransposer.createTranspose(DoOptimizeStorage); 
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
  //
  //This method forms the matrix-matrix sum C = scalarA * op(A) + scalarB * op(B), where

  //Convience typedef
  typedef CrsMatrix<
    Scalar, 
    LocalOrdinal,
    GlobalOrdinal,
    Node,
    SpMatOps> CrsMatrix_t;

  //A and B should already be Filled. C should be an empty pointer.


  TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete() || !B.isFillComplete(), std::runtime_error,
    "EpetraExt::MatrixMatrix::Add ERROR, input matrix A.Filled() or B.Filled() is false,"
    "they are required to be true. (Result matrix C should be an empty pointer)" << std::endl);


  RCP<const CrsMatrix_t> Aprime = null;
  RCP<const CrsMatrix_t> Bprime = null;


  //explicit tranpose A formed as necessary
  if( transposeA ) {
    RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps> theTransposer (Teuchos::rcpFromRef (A));
    Aprime = theTransposer.createTranspose(DoOptimizeStorage);
  }
  else{
    Aprime = rcpFromRef(A);
  }

  //explicit tranpose B formed as necessary
  if( transposeB ) {
    RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps> theTransposer (Teuchos::rcpFromRef (B));
    Bprime = theTransposer.createTranspose(DoOptimizeStorage);
  }
  else{
    Bprime = rcpFromRef(B);
  }

  // allocate or zero the new matrix
  if(C != null)
     C->setAllToScalar(ScalarTraits<Scalar>::zero());
  else
    C = rcp(new CrsMatrix_t(Aprime->getRowMap(), null));

  Array<RCP<const CrsMatrix_t> > Mat = 
    Teuchos::tuple<RCP<const CrsMatrix_t> >(Aprime, Bprime);
  Array<Scalar> scalar = Teuchos::tuple<Scalar>(scalarA, scalarB);

  // do a loop over each matrix to add: A reordering might be more efficient
  for(int k=0;k<2;++k) {
    size_t NumEntries;
    Array<GlobalOrdinal> Indices;
    Array<Scalar> Values;
   
     size_t NumMyRows = Mat[k]->getNodeNumRows();
     GlobalOrdinal Row;
   
     //Loop over rows and sum into C
     for( size_t i = OrdinalTraits<size_t>::zero(); i < NumMyRows; ++i ) {
        Row = Mat[k]->getRowMap()->getGlobalElement(i);
        NumEntries = Mat[k]->getNumEntriesInGlobalRow(Row);
        if(NumEntries == OrdinalTraits<global_size_t>::zero()){
          continue;
        }

        Indices.resize(NumEntries);
        Values.resize(NumEntries);
		    Mat[k]->getGlobalRowCopy(Row, Indices(), Values(), NumEntries);
   
        if( scalar[k] != ScalarTraits<Scalar>::one() )
           for( size_t j = OrdinalTraits<size_t>::zero(); j < NumEntries; ++j ) Values[j] *= scalar[k];
   
        if(C->isFillComplete()) { // Sum in values
           C->sumIntoGlobalValues( Row, Indices, Values);
        } else { // just add it to the unfilled CRS Matrix
           C->insertGlobalValues( Row, Indices, Values);
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
  LocalOrdinal C_firstCol = Bview.colMap->getMinLocalIndex();
  LocalOrdinal C_lastCol = Bview.colMap->getMaxLocalIndex();

  LocalOrdinal C_firstCol_import = OrdinalTraits<LocalOrdinal>::zero();
  LocalOrdinal C_lastCol_import = OrdinalTraits<LocalOrdinal>::invalid();

  ArrayView<const GlobalOrdinal> bcols =Bview.colMap->getNodeElementList();
  ArrayView<const GlobalOrdinal> bcols_import = null;
  if (Bview.importColMap != null) {
    C_firstCol_import = Bview.importColMap->getMinLocalIndex();
    C_lastCol_import = Bview.importColMap->getMaxLocalIndex();

    bcols_import = Bview.importColMap->getNodeElementList();
  }

  size_t C_numCols = C_lastCol - C_firstCol + OrdinalTraits<LocalOrdinal>::one();
  size_t C_numCols_import = C_lastCol_import - C_firstCol_import + OrdinalTraits<LocalOrdinal>::one();

  if (C_numCols_import > C_numCols) C_numCols = C_numCols_import;

  Array<Scalar> dwork = onlyCalculateStructure ? Array<Scalar>() : Array<Scalar>(C_numCols);
  Array<GlobalOrdinal> iwork = Array<GlobalOrdinal>(C_numCols);

  Array<Scalar> C_row_i = dwork;
  Array<GlobalOrdinal> C_cols = iwork;

  size_t C_row_i_length, j, k;

  // Run through all the hash table lookups once and for all
  Array<LocalOrdinal> Acol2Brow(Aview.colMap->getNodeNumElements());
  if(Aview.colMap->isSameAs(*Bview.rowMap)){
    // Maps are the same: Use local IDs as the hash
    for(LocalOrdinal i=Aview.colMap->getMinLocalIndex();i<=Aview.colMap->getMaxLocalIndex();i++)
      Acol2Brow[i]=i;				
  }
  else {
    // Maps are not the same:  Use the map's hash
    for(LocalOrdinal i=Aview.colMap->getMinLocalIndex();i<=Aview.colMap->getMaxLocalIndex();i++)
      Acol2Brow[i]=Bview.rowMap->getLocalElement(Aview.colMap->getGlobalElement(i));
  }

  //To form C = A*B we're going to execute this expression:
  //
  // C(i,j) = sum_k( A(i,k)*B(k,j) )
  //
  //Our goal, of course, is to navigate the data in A and B once, without
  //performing searches for column-indices, etc.

  bool C_filled = C.isFillComplete();

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



    for(k=OrdinalTraits<size_t>::zero(); k<Aview.numEntriesPerRow[i]; ++k) {
      LocalOrdinal Ak=Acol2Brow[Aindices_i[k]];
      Scalar Aval = onlyCalculateStructure ? Teuchos::as<Scalar>(0) : Aval_i[k];

      ArrayView<const LocalOrdinal> Bcol_inds = Bview.indices[Ak];
      ArrayView<const Scalar> Bvals_k = onlyCalculateStructure ? null : Bview.values[Ak];

      C_row_i_length = OrdinalTraits<size_t>::zero();

      if (Bview.remote[Ak]) {
        for(j=OrdinalTraits<size_t>::zero(); j<Bview.numEntriesPerRow[Ak]; ++j) {
          if(!onlyCalculateStructure){
            C_row_i[C_row_i_length] = Aval*Bvals_k[j];
          }
          C_cols[C_row_i_length++] = bcols_import[Bcol_inds[j]];
        }
      }
      else {
        for(j=OrdinalTraits<size_t>::zero(); j<Bview.numEntriesPerRow[Ak]; ++j) {
          if(!onlyCalculateStructure){
            C_row_i[C_row_i_length] = Aval*Bvals_k[j];
          }
          C_cols[C_row_i_length++] = bcols[Bcol_inds[j]];
        }
      }

      //
      //Now put the C_row_i values into C.
      //

      C_filled ?
        C.sumIntoGlobalValues(
          global_row, 
          C_cols.view(OrdinalTraits<size_t>::zero(), C_row_i_length), 
          onlyCalculateStructure ? null : C_row_i.view(OrdinalTraits<size_t>::zero(), C_row_i_length))
        :
        C.insertGlobalValues(
          global_row, 
          C_cols.view(OrdinalTraits<size_t>::zero(), C_row_i_length), 
          onlyCalculateStructure ? null : C_row_i.view(OrdinalTraits<size_t>::zero(), C_row_i_length));

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
      OrdinalTraits<GlobalOrdinal>::invalid(), 
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
