//@HEADER
//************************************************************************
//
//              Isorropia: Partitioning and Load Balancing Package
//                Copyright (2006) Sandia Corporation
//
//Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
//license for use of this work by or on behalf of the U.S. Government.
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
//************************************************************************
//@HEADER

#ifndef _Isorropia_TpetraRedistributor_hpp_
#define _Isorropia_TpetraRedistributor_hpp_

#include <Isorropia_Redistributor.hpp>
//#include <Isorropia_ConfigDefs.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

//#include <Isorropia_Exception.hpp>
#include <EEP_Isorropia_Tpetra.hpp>
#include <EEP_Isorropia_TpetraPartitioner.hpp>

#ifdef USE_UTILS
#include <../../utils/ispatest_epetra_utils.hpp>
#include <../../utils/ispatest_lbeval_utils.hpp>
using namespace ispatest;
#endif

#include <Teuchos_RCP.hpp>

#include <Tpetra_Map_decl.hpp>
#include <Tpetra_Import_decl.hpp>
#include <Tpetra_Vector_decl.hpp>
#include <Tpetra_MultiVector_decl.hpp>
#include <Tpetra_CrsGraph_decl.hpp>
#include <Tpetra_CrsMatrix_decl.hpp>
#include <Tpetra_LinearProblem_decl.hpp>
//#include <Tpetra_Comm_decl.hpp>

namespace Isorropia {

namespace Tpetra {
  class Partitioner;

/** @ingroup partitioning_grp partitioning_rcp_grp partitioning_ptr_grp
     Class which is constructed with a Partitioner instance, and
     provides several methods for redistributing Epetra objects
     given the partitioning computed by the Partitioner object.
*/

class Redistributor : public Isorropia::Redistributor {
public:

  /** @ingroup partitioning_rcp_grp
      This constructor calls the Isorropia::Epetra::Partitioner::partition
      method on the @c partitioner if it has not already been called.
 
      \param partitioner (in) this input partitioner determines the new partitioning
            to be created when Isorropia::Epetra::Redistributor::redistribute is called
   */
  Redistributor(Teuchos::RCP<Isorropia::Epetra::Partitioner> partitioner);

  /** @ingroup partitioning_rcp_grp
      This constructor sets the target map for the redistribution.

      \param target_map (in) this input map determines the new matrices/vectors
      to be created when Isorropia::Epetra::Redistributor::redistribute is
      called
   */
  Redistributor(Teuchos::RCP<Epetra_Map> target_map);

  /** @ingroup partitioning_ptr_grp
      This constructor calls the Isorropia::Epetra::Partitioner::partition
      method on the @c partitioner if it has not already been called.
 
      \param partitioner (in) this input partitioner determines the new partitioning
            to be created when Isorropia::Epetra::Redistributor::redistribute is called
   */
  Redistributor(Isorropia::Epetra::Partitioner *partitioner);

  /** @ingroup partitioning_ptr_grp
      This constructor sets the target map for the redistribution.

      \param target_map (in) this input map determines the new matrices/vectors
      to be created when Isorropia::Epetra::Redistributor::redistribute is
      called
   */
  Redistributor(Epetra_Map *target_map);

  /** 
       Destructor
   */
  virtual ~Redistributor();

  /** @ingroup partitioning_grp
      Method to redistribute a Epetra_SrcDistObject into a
      Epetra_DistObject. The caller is required to have constructed
      the target object using the correct target map.
  */
  void redistribute(const Epetra_SrcDistObject& src,
		    Epetra_DistObject& target);

  /** @ingroup partitioning_rcp_grp
      Method to accept a Epetra_CrsGraph object, and
      return a redistributed Epetra_CrsGraph object.

      \param input_graph (in) the graph for which we want a new graph that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \param callFillComplete (in) The new graph is FillComplete'd if callFillComplete is @c true. 
      In that case, the range map is set to equal the row map. 
      The domain map will equal the range map, unless the
      input_graph has different domain and range maps, in which case
      the original domain map is preserved.  By default callFillComplete is @c true.

      \return a reference counted pointer to the new redistributed graph 
  */
  Teuchos::RCP<Epetra_CrsGraph>
     redistribute(const Epetra_CrsGraph& input_graph, bool callFillComplete= true);

  /** @ingroup partitioning_ptr_grp
      Method to accept a Epetra_CrsGraph object, and
      return a redistributed Epetra_CrsGraph object.

      \param input_graph (in) the graph for which we want a new graph that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \param outputGraphPtr (out) pointer to the new redistributed graph 

      \param callFillComplete (in) The new graph is FillComplete'd if callFillComplete is @c true. 
      In that case, the range map is set to equal the row map. 
      The domain map will equal the range map, unless the
      input_graph has different domain and range maps, in which case
      the original domain map is preserved.  By default callFillComplete is @c true.

  */
  void redistribute(const Epetra_CrsGraph& input_graph, Epetra_CrsGraph * &outputGraphPtr, bool callFillComplete= true);

  /** @ingroup partitioning_rcp_grp
      Method to accept a Epetra_CrsMatrix object, and
      return a redistributed Epetra_CrsMatrix object.

      \param input_matrix (in) the matrix for which we want a new matrix that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \param callFillComplete (in) The new matrix is FillComplete'd if callFillComplete is @c true. 
      In that case, the range map is set to equal the row map. 
      The domain map will equal the range map, unless the
      input_matrix has different domain and range maps, in which case
      the original domain map is preserved.  By default callFillComplete is @c true.

      \return a reference counted pointer to the new redistributed matrix
  */
  Teuchos::RCP<Epetra_CrsMatrix>
     redistribute(const Epetra_CrsMatrix& input_matrix, bool callFillComplete= true);

  /** @ingroup partitioning_ptr_grp
      Method to accept a Epetra_CrsMatrix object, and
      return a redistributed Epetra_CrsMatrix object.

      \param inputMatrix (in) the matrix for which we want a new matrix that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \param outputMatrix (out) pointer to the new redistributed matrix

      \param callFillComplete (in) The new matrix is FillComplete'd if callFillComplete is @c true. 
      In that case, the range map is set to equal the row map. 
      The domain map will equal the range map, unless the
      input_matrix has different domain and range maps, in which case
      the original domain map is preserved.  By default callFillComplete is @c true.
  */

  void redistribute(const Epetra_CrsMatrix& inputMatrix, Epetra_CrsMatrix * &outputMatrix, bool callFillComplete= true);

  /** @ingroup partitioning_rcp_grp
      Method to accept a Epetra_RowMatrix object, and
      return a redistributed Epetra_CrsMatrix object.

      \param input_matrix (in) the row matrix for which we want a new matrix that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \param callFillComplete (in) The new matrix is FillComplete'd if callFillComplete is @c true. 
      In that case, the range map is set to equal the row map. 
      The domain map will equal the range map, unless the
      input_matrix has different domain and range maps, in which case
      the original domain map is preserved.  By default callFillComplete is @c true.

      \return a reference counted pointer to the new redistributed matrix
  */
  Teuchos::RCP<Epetra_CrsMatrix>
     redistribute(const Epetra_RowMatrix& input_matrix, bool callFillComplete= true);


  /** @ingroup partitioning_ptr_grp
      Method to accept a Epetra_RowMatrix object, and
      return a redistributed Epetra_CrsMatrix object.

      \param inputMatrix (in) the row matrix for which we want a new matrix that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \param outputMatrix (out) pointer to the new redistributed matrix

      \param callFillComplete (in) The new matrix is FillComplete'd if callFillComplete is @c true. 
      In that case, the range map is set to equal the row map. 
      The domain map will equal the range map, unless the
      input_matrix has different domain and range maps, in which case
      the original domain map is preserved.  By default callFillComplete is @c true.
  */
     void 
     redistribute(const Epetra_RowMatrix& inputMatrix, Epetra_CrsMatrix * &outputMatrix, bool callFillComplete= true);

  /** @ingroup partitioning_rcp_grp
      Method to accept a Epetra_Vector object, and
      return a redistributed Epetra_Vector object.

      \param input_vector (in) the vector for which we want a new vector that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \return a reference counted pointer to the new redistributed vector
  */
  Teuchos::RCP<Epetra_Vector>
     redistribute(const Epetra_Vector& input_vector);

  /** @ingroup partitioning_ptr_grp
      Method to accept a Epetra_Vector object, and
      return a redistributed Epetra_Vector object.

      \param inputVector (in) the vector for which we want a new vector that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \param outputVector (out) pointer to the new redistributed vector
  */
  void
  redistribute(const Epetra_Vector& inputVector, Epetra_Vector * &outputVector);

  /** @ingroup partitioning_rcp_grp 
      Method to accept a Epetra_MultiVector object, and
      return a redistributed Epetra_MultiVector object.

      \param input_vector (in) the multi vector for which we want a new multi vector that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \return a reference counted pointer to the new redistributed multi vector

  */
  Teuchos::RCP<Epetra_MultiVector>  
     redistribute(const Epetra_MultiVector& input_vector);


  /** @ingroup partitioning_ptr_grp 
      Method to accept a Epetra_MultiVector object, and
      return a redistributed Epetra_MultiVector object.

      \param inputVector (in) the multi vector for which we want a new multi vector that is distributed
                      according to the partitioner with which this Redistributor was
                      created.

      \param outputVector (out) a reference counted pointer to the new redistributed multi vector

  */
  void 
  redistribute(const Epetra_MultiVector& inputVector, Epetra_MultiVector * &outputVector);

  /** @ingroup partitioning_grp 
      Reverse redistribute an Epetra_Vector.

      \param input_vector (in) a vector that is distributed according to the partitioner that was used to create this Redistributor

      \param output_vector (out) a copy of the @c input_vector which has been redistributed according
                    to the reverse of the partitioner that was used to create this Redistributor

  */
  void
     redistribute_reverse(const Epetra_Vector& input_vector, Epetra_Vector& output_vector);

  /** @ingroup partitioning_grp 
      Reverse redistribute an Epetra_MultiVector.

      \param input_vector (in) a multi vector that is distributed according to the partitioner that was used to create this Redistributor

      \param output_vector (out) a copy of the @c input_vector which has been redistributed according
                    to the reverse of the partitioner that was used to create this Redistributor
  */
  void
     redistribute_reverse(const Epetra_MultiVector& input_vector, Epetra_MultiVector& output_vector);

  Epetra_Import &get_importer() { return *importer_;}

private:
  /** @ingroup partitioning_grp
      Create an importer object to be used in the redistribution
      \param src_map (in) the map describing the pattern of the import operation
   */
  void create_importer(const Epetra_BlockMap& src_map);

  Teuchos::RCP<Isorropia::Epetra::Partitioner> partitioner_;
  Teuchos::RCP<Epetra_Import> importer_;
  Teuchos::RCP<Epetra_Map> target_map_;

}; //class Redistributor

Redistributor::Redistributor(Teuchos::RCP<Isorropia::Epetra::Partitioner> partitioner)
  : partitioner_(partitioner),
  importer_(),
  target_map_()
{
  if (!partitioner_->alreadyComputed()) {
    partitioner_->partition();
  }
}

Redistributor::Redistributor(Isorropia::Epetra::Partitioner *partitioner)
  : partitioner_(Teuchos::RCP<Isorropia::Epetra::Partitioner>(partitioner,false)),
  importer_(),
  target_map_()
{
  if (!partitioner_->alreadyComputed()) {
    partitioner_->partition();
  }
}

Redistributor::Redistributor(Teuchos::RCP<Epetra_Map> target_map)
  : partitioner_(),
  importer_(),
  target_map_(target_map)
{
    // Do not partition, target map is given.
}

Redistributor::Redistributor(Epetra_Map *target_map)
  : partitioner_(),
  importer_(),
  target_map_(Teuchos::RCP<Epetra_Map>(target_map, false))
{
    // Do not partition, target map is given.
}

Redistributor::~Redistributor()
{
}

void
Redistributor::redistribute(const Epetra_SrcDistObject& src,
				 Epetra_DistObject& target)
{
  create_importer(src.Map());

  target.Import(src, *importer_, Insert);
}

Teuchos::RCP<Epetra_CrsGraph>
Redistributor::redistribute(const Epetra_CrsGraph& input_graph, bool callFillComplete)
{
  Epetra_CrsGraph *outputGraphPtr=0;
  redistribute(input_graph, outputGraphPtr, callFillComplete);

  return Teuchos::RCP<Epetra_CrsGraph>(outputGraphPtr);
}

void 
Redistributor::redistribute(const Epetra_CrsGraph& input_graph, Epetra_CrsGraph * &outputGraphPtr, bool callFillComplete)
{
  create_importer(input_graph.RowMap());

  // First obtain the length of each of my new rows

  int myOldRows = input_graph.NumMyRows();
  int myNewRows = target_map_->NumMyElements();

  double *nnz = new double [myOldRows];
  for (int i=0; i < myOldRows; i++){
    nnz[i] = input_graph.NumMyIndices(i);
  }

  Epetra_Vector oldRowSizes(Copy, input_graph.RowMap(), nnz);

  if (myOldRows)
    delete [] nnz;

  Epetra_Vector newRowSizes(*target_map_);

  newRowSizes.Import(oldRowSizes, *importer_, Insert);

  int *rowSize = new int [myNewRows];
  for (int i=0; i< myNewRows; i++){
    rowSize[i] = static_cast<int>(newRowSizes[i]);
  }

  // Receive new rows, send old rows

  outputGraphPtr = new Epetra_CrsGraph(Copy, *target_map_, rowSize, true);

  if (myNewRows)
    delete [] rowSize;

  outputGraphPtr->Import(input_graph, *importer_, Insert);

  // Set the new domain map such that
  // (a) if old DomainMap == old RangeMap, preserve this property,
  // (b) otherwise, let the new DomainMap be the old DomainMap 
  const Epetra_BlockMap *newDomainMap;
  if (input_graph.DomainMap().SameAs(input_graph.RangeMap()))
     newDomainMap = &(outputGraphPtr->RangeMap());
  else
     newDomainMap = &(input_graph.DomainMap());

  if (callFillComplete && (!outputGraphPtr->Filled()))
    outputGraphPtr->FillComplete(*newDomainMap, *target_map_);

  return;
}

Teuchos::RCP<Epetra_CrsMatrix>
Redistributor::redistribute(const Epetra_CrsMatrix& inputMatrix, bool callFillComplete)
{
  Epetra_CrsMatrix *outputMatrix=0;
  redistribute(inputMatrix, outputMatrix, callFillComplete);

  return Teuchos::RCP<Epetra_CrsMatrix>(outputMatrix);
}

void 
Redistributor::redistribute(const Epetra_CrsMatrix& inputMatrix, Epetra_CrsMatrix * &outputMatrix, bool callFillComplete)
{
  create_importer(inputMatrix.RowMap());

  // First obtain the length of each of my new rows

  int myOldRows = inputMatrix.NumMyRows();
  int myNewRows = target_map_->NumMyElements();

  double *nnz = new double [myOldRows];
  for (int i=0; i < myOldRows; i++){
    nnz[i] = inputMatrix.NumMyEntries(i);
  }

  Epetra_Vector oldRowSizes(Copy, inputMatrix.RowMap(), nnz);

  if (myOldRows)
    delete [] nnz;

  Epetra_Vector newRowSizes(*target_map_);

  newRowSizes.Import(oldRowSizes, *importer_, Insert);

  int *rowSize=0;
  if(myNewRows){
    rowSize = new int [myNewRows];
    for (int i=0; i< myNewRows; i++){
      rowSize[i] = static_cast<int>(newRowSizes[i]);
    }
  }

  // Receive new rows, send old rows

  outputMatrix = new Epetra_CrsMatrix(Copy, *target_map_, rowSize, true);

  if (myNewRows)
    delete [] rowSize;

  outputMatrix->Import(inputMatrix, *importer_, Insert);

  // Set the new domain map such that
  // (a) if old DomainMap == old RangeMap, preserve this property,
  // (b) otherwise, let the new DomainMap be the old DomainMap 
  const Epetra_Map *newDomainMap;
  if (inputMatrix.DomainMap().SameAs(inputMatrix.RangeMap()))
     newDomainMap = &(outputMatrix->RangeMap());
  else
     newDomainMap = &(inputMatrix.DomainMap());

  if (callFillComplete && (!outputMatrix->Filled()))
    outputMatrix->FillComplete(*newDomainMap,  *target_map_);

#ifdef USE_UTILS
  // TODO: "PRINT ZOLTAN METRICS" should specify graph, hypergraph, or imbalance
  //     instead of just amount of output
  //     Metrics should be printed for each redistribute() method

  if (!Teuchos::is_null(partitioner_) &&
        partitioner_->printZoltanMetrics() > 0){
    std::vector<double> balance(2), cutWgt(2), cutn(2), cutl(2);
    std::vector<int> numCuts(2);
    CostDescriber *defaultCosts = NULL;

    Teuchos::RCP<CostDescriber> &cptr = partitioner_->getCosts();

    if (Teuchos::is_null(cptr)){
      defaultCosts = new CostDescriber;
      cptr = Teuchos::rcp(defaultCosts);
    }

    int fail = cptr->compareBeforeAndAfterGraph( inputMatrix, *outputMatrix, *importer_,
                      balance, numCuts, cutWgt, cutn, cutl);

    if (!fail){

      if (partitioner_->printZoltanMetrics() > 1){
        ispatest::printRowMatrix(inputMatrix, std::cout, "BEFORE:", true);
      }
      inputMatrix.Comm().Barrier();
      if (inputMatrix.Comm().MyPID() == 0){
        std::cout << "BEFORE: Imbalance " << balance[0] << ", cuts " << numCuts[0];
        std::cout << ", cut weight " << cutWgt[0];
        std::cout << ", CUTN " << cutn[0];
        std::cout << ", CUTL " << cutl[0] << std::endl;
      }

      if (partitioner_->printZoltanMetrics() > 1){
        ispatest::printRowMatrix(*outputMatrix, std::cout, "AFTER:", true);
      }
      inputMatrix.Comm().Barrier();
      if (inputMatrix.Comm().MyPID() == 0){
        std::cout << " AFTER: Imbalance " << balance[1] << ", cuts " << numCuts[1];
        std::cout << ", cut weight " << cutWgt[1];
        std::cout << ", CUTN " << cutn[1];
        std::cout << ", CUTL " << cutl[1] << std::endl;
      }
    }
    else{
      if (inputMatrix.Comm().MyPID() == 0){
        std::cout << " Error computing Zoltan quality metrics" << std::endl;
      }
    }
    inputMatrix.Comm().Barrier();
  }
#endif

  return;
}

Teuchos::RCP<Epetra_CrsMatrix>
Redistributor::redistribute(const Epetra_RowMatrix& inputMatrix, bool callFillComplete)
{
  Epetra_CrsMatrix *outputMatrix = 0;
  redistribute(inputMatrix,outputMatrix,callFillComplete);
  return Teuchos::RCP<Epetra_CrsMatrix>(outputMatrix);
}

void
Redistributor::redistribute(const Epetra_RowMatrix& inputMatrix, Epetra_CrsMatrix * &outputMatrix, bool callFillComplete)
{


  create_importer(inputMatrix.RowMatrixRowMap());
 // First obtain the length of each of my new rows

  int myOldRows = inputMatrix.NumMyRows();
  int myNewRows = target_map_->NumMyElements();


  double *nnz = new double [myOldRows];
  int val;
  for (int i=0; i < myOldRows; i++){
    inputMatrix.NumMyRowEntries(i, val);
    nnz[i] = static_cast<double>(val);
  }

  Epetra_Vector oldRowSizes(Copy, inputMatrix.RowMatrixRowMap(), nnz);

  if (myOldRows)
    delete [] nnz;

  Epetra_Vector newRowSizes(*target_map_);

  newRowSizes.Import(oldRowSizes, *importer_, Insert);

  int *rowSize = new int [myNewRows];
  for (int i=0; i< myNewRows; i++){
    rowSize[i] = static_cast<int>(newRowSizes[i]);
  }

  // Receive new rows, send old rows

  outputMatrix = new Epetra_CrsMatrix(Copy, *target_map_, rowSize, true);

  if (myNewRows)
    delete [] rowSize;

  outputMatrix->Import(inputMatrix, *importer_, Insert);

  // Set the new domain map such that
  // (a) if old DomainMap == old RangeMap, preserve this property,
  // (b) otherwise, use the original OperatorDomainMap
  //if (inputMatrix.NumGlobalRows() == inputMatrix.NumGlobalCols())
  if (inputMatrix.OperatorDomainMap().SameAs(inputMatrix.OperatorRangeMap())){
    if (callFillComplete && (!outputMatrix->Filled()))
      outputMatrix->FillComplete();
  }
  else {
    if (callFillComplete && (!outputMatrix->Filled()))
      outputMatrix->FillComplete(inputMatrix.OperatorDomainMap(), *target_map_);
  }

#ifdef USE_UTILS
  // TODO: "PRINT ZOLTAN METRICS" should specify graph, hypergraph, or imbalance
  //     instead of just amount of output
  //     Metrics should be printed for each redistribute() method

  if (!Teuchos::is_null(partitioner_) &&
        partitioner_->printZoltanMetrics() > 0){
    std::vector<double> balance(2), cutWgt(2), cutn(2), cutl(2);
    std::vector<int> numCuts(2);
    CostDescriber *defaultCosts = NULL;

    Teuchos::RCP<CostDescriber> &cptr = partitioner_->getCosts();

    if (Teuchos::is_null(cptr)){
      defaultCosts = new CostDescriber;
      cptr = Teuchos::rcp(defaultCosts);
    }

    int fail = cptr->compareBeforeAndAfterGraph( inputMatrix, *outputMatrix, *importer_,
                      balance, numCuts, cutWgt, cutn, cutl);
    if (!fail){

      if (partitioner_->printZoltanMetrics() > 1){
        ispatest::printRowMatrix(inputMatrix, std::cout, "BEFORE:", true);
      } 
      inputMatrix.Comm().Barrier();
      if (inputMatrix.Comm().MyPID() == 0){
        std::cout << "BEFORE: Imbalance " << balance[0] << ", cuts " << numCuts[0];
        std::cout << ", cut weight " << cutWgt[0];
        std::cout << ", CUTN " << cutn[0];
        std::cout << ", CUTL " << cutl[0] << std::endl;
      } 
  
      if (partitioner_->printZoltanMetrics() > 1){
        ispatest::printRowMatrix(*outputMatrix, std::cout, "AFTER:", true);
      } 
      inputMatrix.Comm().Barrier();
      if (inputMatrix.Comm().MyPID() == 0){
        std::cout << " AFTER: Imbalance " << balance[1] << ", cuts " << numCuts[1];
        std::cout << ", cut weight " << cutWgt[1];
        std::cout << ", CUTN " << cutn[1];
        std::cout << ", CUTL " << cutl[1] << std::endl;
      } 
    } 
    else{
      if (inputMatrix.Comm().MyPID() == 0){
        std::cout << " Error computing Zoltan quality metrics" << std::endl;
      }
    } 
    inputMatrix.Comm().Barrier();
  }
#endif

  return;
}

Teuchos::RCP<Epetra_Vector>
Redistributor::redistribute(const Epetra_Vector& input_vector)
{
  Epetra_Vector *outputVector = 0;
  redistribute(input_vector,outputVector);

  return Teuchos::RCP<Epetra_Vector>(outputVector);
}

void
Redistributor::redistribute(const Epetra_Vector& inputVector, Epetra_Vector * &outputVector)
{
  create_importer(inputVector.Map());

  outputVector = new Epetra_Vector(*target_map_);

  outputVector->Import(inputVector, *importer_, Insert);

  return;
}

Teuchos::RCP<Epetra_MultiVector>
Redistributor::redistribute(const Epetra_MultiVector& input_vector)
{
  Epetra_MultiVector *outputVector=0;
  redistribute(input_vector,outputVector);

  return Teuchos::RCP<Epetra_MultiVector>(outputVector);
}


void
Redistributor::redistribute(const Epetra_MultiVector& inputVector, Epetra_MultiVector * &outputVector)
{
  create_importer(inputVector.Map());

  outputVector = new Epetra_MultiVector(*target_map_, inputVector.NumVectors());

  outputVector->Import(inputVector, *importer_, Insert);

  return;
}


// Reverse redistribute methods (for vectors). 

void
Redistributor::redistribute_reverse(const Epetra_Vector& input_vector, Epetra_Vector& output_vector)
{
  create_importer(input_vector.Map());

  // Export using the importer
  output_vector.Export(input_vector, *importer_, Insert);

}

void
Redistributor::redistribute_reverse(const Epetra_MultiVector& input_vector, Epetra_MultiVector& output_vector)
{
  create_importer(input_vector.Map());

  // Export using the importer
  output_vector.Export(input_vector, *importer_, Insert);

}

void Redistributor::create_importer(const Epetra_BlockMap& src_map)
{

  if (!Teuchos::is_null(partitioner_) && partitioner_->numProperties() >
                                src_map.Comm().NumProc()) {
    throw Isorropia::Exception("Cannot redistribute: Too many parts for too few processors.");
  }

  if (Teuchos::is_null(target_map_) && !Teuchos::is_null(partitioner_))
      target_map_ = partitioner_->createNewMap();

  importer_ = Teuchos::rcp(new Epetra_Import(*target_map_, src_map));

}

}//namespace Tpetra

}//namespace Isorropia

#endif


#if defined(Isorropia_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Isorropia package is deprecated"
#endif
#endif

