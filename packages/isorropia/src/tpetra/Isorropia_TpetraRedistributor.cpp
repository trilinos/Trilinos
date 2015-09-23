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

#include <Isorropia_TpetraRedistributor.hpp>
#include <Isorropia_Exception.hpp>
#include <Isorropia_Tpetra.hpp>
#include <Isorropia_TpetraPartitioner.hpp>

#include <Teuchos_RCP.hpp>


namespace Isorropia {

namespace Tpetra {

#ifdef HAVE_ISORROPIA_TPETRA

template<typename Node>
Redistributor<Node>::Redistributor(Teuchos::RCP<Isorropia::Tpetra::Partitioner<Node> > partitioner)
  : partitioner_(partitioner),
  importer_(),
  target_map_(),
  created_importer_(false)
{
  if (!partitioner_->alreadyComputed()) {
    partitioner_->partition();
  }
}

template<class Node>
Redistributor<Node>::Redistributor(Isorropia::Tpetra::Partitioner<Node> *partitioner)
  : partitioner_(Teuchos::RCP<Isorropia::Tpetra::Partitioner<Node> >(partitioner,false)),
  importer_(),
  target_map_(),
  created_importer_(false)
{
  if (!partitioner_->alreadyComputed()) {
    partitioner_->partition();
  }
}

template<typename Node>
Redistributor<Node>::~Redistributor()
{
}

template<class Node>
Teuchos::RCP< ::Tpetra::CrsGraph<int,int,Node> >
Redistributor<Node>::redistribute(const ::Tpetra::CrsGraph<int,int,Node>& input_graph, bool callFillComplete)
{
  ::Tpetra::CrsGraph<int,int,Node> *outputGraphPtr=0;
  redistribute(input_graph, outputGraphPtr, callFillComplete);

  return Teuchos::RCP< ::Tpetra::CrsGraph<int,int,Node> >(outputGraphPtr);
}

template<class Node>
void 
Redistributor<Node>::redistribute(const ::Tpetra::CrsGraph<int,int,Node>& input_graph, ::Tpetra::CrsGraph<int,int,Node> * &outputGraphPtr, bool callFillComplete)
{
  if (!created_importer_) {
    create_importer(input_graph.RowMap());
  }

  // First obtain the length of each of my new rows

  int myOldRows = input_graph.NumMyRows();
  int myNewRows = target_map_->NumMyElements();

  double *nnz = new double [myOldRows];
  for (int i=0; i < myOldRows; i++){
    nnz[i] = input_graph.NumMyIndices(i);
  }

  Teuchos::ArrayView<double>::ArrayView   myNNZView(nnz, myOldRows);

  ::Tpetra::Vector<double,int,int,Node> oldRowSizes(input_graph.RowMap(), myNNZView);

  if (myOldRows)
    delete [] nnz;

  ::Tpetra::Vector<double,int,int,Node> newRowSizes(*target_map_);

  newRowSizes.Import(oldRowSizes, *importer_, ::Tpetra::INSERT);

  int *rowSize = new int [myNewRows];
  for (int i=0; i< myNewRows; i++){
    rowSize[i] = static_cast<int>(newRowSizes[i]);
  }

  Teuchos::ArrayView<int>::ArrayView rowSizeView(rowSize, myNewRows);

  // Receive new rows, send old rows

  outputGraphPtr = new ::Tpetra::CrsGraph<int,int,Node>(*target_map_, rowSizeView, true);

  if (myNewRows)
    delete [] rowSize;

  outputGraphPtr->Import(input_graph, *importer_, ::Tpetra::INSERT);

  // Set the new domain map such that
  // (a) if old DomainMap == old RangeMap, preserve this property,
  // (b) otherwise, let the new DomainMap be the old DomainMap 
  const ::Tpetra::Map<int,int,Node> *newDomainMap;
  if (input_graph.DomainMap().SameAs(input_graph.RangeMap()))
     newDomainMap = &(outputGraphPtr->RangeMap());
  else
     newDomainMap = &(input_graph.DomainMap());

  if (callFillComplete && (!outputGraphPtr->Filled()))
    outputGraphPtr->FillComplete(*newDomainMap, *target_map_);

  return;
}

template<class Node>
Teuchos::RCP< ::Tpetra::CrsMatrix<double,int,int,Node> >
Redistributor<Node>::redistribute(const ::Tpetra::CrsMatrix<double,int,int,Node>& inputMatrix, bool callFillComplete)
{
  ::Tpetra::CrsMatrix<double,int,int,Node> *outputMatrix=0;
  redistribute(inputMatrix, outputMatrix, callFillComplete);

  return Teuchos::RCP< ::Tpetra::CrsMatrix<double,int,int,Node> >(outputMatrix);
}

template<class Node>
void 
Redistributor<Node>::redistribute(const ::Tpetra::CrsMatrix<double,int,int,Node>& inputMatrix, 
                                  ::Tpetra::CrsMatrix<double,int,int,Node> * &outputMatrix, bool callFillComplete)
{
  if (!created_importer_) {
    create_importer(inputMatrix.RowMap());
  }

  // First obtain the length of each of my new rows

  int myOldRows = inputMatrix.NumMyRows();
  int myNewRows = target_map_->NumMyElements();

  double *nnz = new double [myOldRows];
  for (int i=0; i < myOldRows; i++){
    nnz[i] = inputMatrix.NumMyEntries(i);
  }

  Teuchos::ArrayView<double>::ArrayView   myNNZView(nnz, myOldRows);

  ::Tpetra::Vector<double,int,int,Node> oldRowSizes(inputMatrix.RowMap(), myNNZView);

  if (myOldRows)
    delete [] nnz;

  ::Tpetra::Vector<double,int,int,Node> newRowSizes(*target_map_);

  newRowSizes.Import(oldRowSizes, *importer_, ::Tpetra::INSERT);

  int *rowSize=0;
  if(myNewRows){
    rowSize = new int [myNewRows];
    for (int i=0; i< myNewRows; i++){
      rowSize[i] = static_cast<int>(newRowSizes[i]);
    }
  }

  Teuchos::ArrayView<int>::ArrayView rowSizeView(rowSize, myNewRows);

  // Receive new rows, send old rows

  outputMatrix = new ::Tpetra::CrsMatrix<double,int,int,Node> (*target_map_, rowSizeView, true);

  if (myNewRows)
    delete [] rowSize;

  outputMatrix->Import(inputMatrix, *importer_, ::Tpetra::INSERT);

  // Set the new domain map such that
  // (a) if old DomainMap == old RangeMap, preserve this property,
  // (b) otherwise, let the new DomainMap be the old DomainMap 
  const ::Tpetra::Map<int,int,Node> *newDomainMap;
  if (inputMatrix.DomainMap().SameAs(inputMatrix.RangeMap()))
     newDomainMap = &(outputMatrix->RangeMap());
  else
     newDomainMap = &(inputMatrix.DomainMap());

  if (callFillComplete && (!outputMatrix->Filled()))
    outputMatrix->FillComplete(*newDomainMap,  *target_map_);

  return;
}

template<class Node>
Teuchos::RCP< ::Tpetra::Vector<double,int,int,Node> >
Redistributor<Node>::redistribute(const ::Tpetra::Vector<double,int,int,Node> & input_vector)
{
  ::Tpetra::Vector<double,int,int,Node> *outputVector = 0;
  redistribute(input_vector,outputVector);

  return Teuchos::RCP< ::Tpetra::Vector<double,int,int,Node> >(outputVector);
}

template<class Node>
void
Redistributor<Node>::redistribute(const ::Tpetra::Vector<double,int,int,Node>& inputVector, 
                                  ::Tpetra::Vector<double,int,int,Node> * &outputVector)
{
  if (!created_importer_) {
    create_importer(inputVector.Map());
  }

  outputVector = new ::Tpetra::Vector<double,int,int,Node>(*target_map_);

  outputVector->Import(inputVector, *importer_, ::Tpetra::INSERT);

  return;
}

template<class Node>
Teuchos::RCP< ::Tpetra::MultiVector<double,int,int,Node> >
Redistributor<Node>::redistribute(const  ::Tpetra::MultiVector<double,int,int,Node> & input_vector)
{
   ::Tpetra::MultiVector<double,int,int,Node>  *outputVector=0;
  redistribute(input_vector,outputVector);

  return Teuchos::RCP< ::Tpetra::MultiVector<double,int,int,Node> >(outputVector);
}


template<class Node>
void
Redistributor<Node>::redistribute(const  ::Tpetra::MultiVector<double,int,int,Node> & inputVector,  ::Tpetra::MultiVector<double,int,int,Node>  * &outputVector)
{
  if (!created_importer_) {
    create_importer(inputVector.Map());
  }

  outputVector = new  ::Tpetra::MultiVector<double,int,int,Node> (*target_map_, inputVector.NumVectors());

  outputVector->Import(inputVector, *importer_, ::Tpetra::INSERT);

  return;
}




// Reverse redistribute methods (for vectors). 

template<class Node>
void
Redistributor<Node>::redistribute_reverse(const  ::Tpetra::Vector<double,int,int,Node> & input_vector,  ::Tpetra::Vector<double,int,int,Node> & output_vector)
{
  if (!created_importer_) {
    create_importer(input_vector.Map());
  }

  // Export using the importer
  output_vector.Export(input_vector, *importer_, ::Tpetra::INSERT);

}

template<class Node>
void
Redistributor<Node>::redistribute_reverse(const  ::Tpetra::MultiVector<double,int,int,Node> & input_vector,  ::Tpetra::MultiVector<double,int,int,Node> & output_vector)
{
  if (!created_importer_) {
    create_importer(input_vector.Map());
  }

  // Export using the importer
  output_vector.Export(input_vector, *importer_, ::Tpetra::INSERT);

}

template<class Node>
void Redistributor<Node>::create_importer(const ::Tpetra::Map<int,int,Node>& src_map)
{
  if (created_importer_) return;

  if (partitioner_->numProperties() > src_map.Comm().NumProc()) {
    throw Isorropia::Exception("Cannot redistribute: Too many parts for too few processors.");
  }

  target_map_ = partitioner_->createNewMap();

  importer_ = Teuchos::rcp(new ::Tpetra::Import<int,int,Node>(src_map, *target_map_));

  created_importer_ = true;
}

#endif //HAVE_ISORROPIA_TPETRA

}//namespace Tpetra

}//namespace Isorropia

