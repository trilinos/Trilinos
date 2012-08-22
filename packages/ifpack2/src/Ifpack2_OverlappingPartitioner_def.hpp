/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK2_OVERLAPPINGPARTITIONER_DEF_HPP
#define IFPACK2_OVERLAPPINGPARTITIONER_DEF_HPP
#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_OverlappingPartitioner_decl.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include <vector>
#include <string>

namespace Ifpack2 {
//==============================================================================
template<class GraphType>
OverlappingPartitioner<GraphType>::OverlappingPartitioner(const Teuchos::RCP<const Tpetra::RowGraph<LocalOrdinal,GlobalOrdinal,Node> >& Graph) :
  NumLocalParts_(1),
  Graph_(Graph),
  OverlappingLevel_(0),
  IsComputed_(false),
  verbose_(false)
{
}

//==============================================================================
template<class GraphType>
OverlappingPartitioner<GraphType>::~Ifpack2_OverlappingPartitioner()
{
}

//==============================================================================
template<class GraphType>
size_t OverlappingPartitioner<GraphType>::numLocalParts() const 
{
  return(NumLocalParts_);
}

//==============================================================================
template<class GraphType>
size_t OverlappingPartitioner<GraphType>::overlappingLevel() const
{
  return(OverlappingLevel_);
}

//==============================================================================
template<class GraphType>
typename GraphType::local_ordinal_type OverlappingPartitioner<GraphType>::operator() (LocalOrdinal MyRow) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(((MyRow < 0) || ((size_t)MyRow > Graph_->getNodeNumRows())), std::runtime_error, "Ifpack2::OverlappingPartitioner::operator() invalid row.");
   
  return(Partition_[MyRow]);
}


//==============================================================================
template<class GraphType>
typename GraphType::local_ordinal_type OverlappingPartitioner<GraphType>::operator() (LocalOrdinal i, LocalOrdinal j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION( ((i < 0) || ((i > NumLocalParts_) || (j < 0) || (j > Parts_[i].size()))),
			      std::runtime_error, "Ifpack2::OverlappingPartitioner::operator() invalid row or node.");
  return(Parts_[i][j]);
}

//==============================================================================
template<class GraphType>
size_t OverlappingPartitioner<GraphType>::numRowsInPart(const LocalOrdinal Part) const
{
  TEUCHOS_TEST_FOR_EXCEPTION( ((Part < 0) || (Part > NumLocalParts_)),
			      std::runtime_error, "Ifpack2::OverlappingPartitioner::numRowsInPart() invalid partition.");
  return(Parts_[Part].size());
}

//==============================================================================
template<class GraphType>
void OverlappingPartitioner<GraphType>::rowsInPart(const LocalOrdinal Part,  Teuchos::ArrayRCP<LocalOrdinal> & List) const
{
  // Let numRowsInPart do the sanity checking...
  size_t numrows= numRowsInPart(Part); 
  for (size_t i = 0 ; i < numrows; i++)
    List[i] = Parts_[Part][i];
}

//==============================================================================
template<class GraphType>
Teuchos::ArrayView<const typename GraphType::local_ordinal_type>  OverlappingPartitioner<GraphType>::nonOverlappingPartition() const
{
  return(Partition_.view(0, Graph_->getNodeNumRows()));
}

//==============================================================================
template<class GraphType>
void OverlappingPartitioner<GraphType>::setParameters(Teuchos::ParameterList& List)
{
  NumLocalParts_    = List.get("partitioner: local parts", NumLocalParts_);
  OverlappingLevel_ = List.get("partitioner: overlap", OverlappingLevel_);
  verbose_          = List.get("partitioner: print level", verbose_);

  if (NumLocalParts_ < 0)
    NumLocalParts_ = Graph_->getNodeNumRows() / (-NumLocalParts_);
  if (NumLocalParts_ == 0)
    NumLocalParts_ = 1;
  
  // Sanity checking
  TEUCHOS_TEST_FOR_EXCEPTION( ((NumLocalParts_ < 0) || (size_t)NumLocalParts_ > Graph_->getNodeNumRows()),
			      std::runtime_error, "Ifpack2::OverlappingPartitioner::setParameters() invalid NumLocalParts_");
  TEUCHOS_TEST_FOR_EXCEPTION( (OverlappingLevel_ < 0),
			      std::runtime_error, "Ifpack2::OverlappingPartitioner::setParameters() invalid OverlappingLevel_");

  setPartitionParameters(List);
}

//==============================================================================
template<class GraphType>
void OverlappingPartitioner<GraphType>::compute()
{
  TEUCHOS_TEST_FOR_EXCEPTION( ((NumLocalParts_ < 1) ||  (OverlappingLevel_ < 0)),
			      std::runtime_error, "Ifpack2::OverlappingPartitioner::compute() invalid NumLocalParts_ or OverlappingLevel_");

  std::string PrintMsg_("OverlappingPartitioner: ");

  // some output
  if (verbose_ && (Graph_->getComm()->getRank() == 0)) {
    std::cout << PrintMsg_ << "Number of local parts          = " << NumLocalParts_ << std::endl;
    std::cout << PrintMsg_ << "Approx. Number of global parts = " 
	      << NumLocalParts_ * Graph_->getComm()->getSize() << std::endl;
    std::cout << PrintMsg_ << "Amount of overlap              = " << OverlappingLevel_ << std::endl;
  }

  // 1.- allocate memory 
  Partition_.resize(Graph_->getNodeNumRows());
  Parts_.resize(NumLocalParts_);

  // 2.- sanity checks on input graph
  TEUCHOS_TEST_FOR_EXCEPTION( (!Graph_->isFillComplete() || Graph_->getGlobalNumRows() != Graph_->getGlobalNumCols()),
			      std::runtime_error, "Ifpack2::OverlappingPartitioner::compute() input graph error");
 
  // 3.- perform non-overlapping partition 
  computePartitions();

  // 4.- compute the partitions with overlapping
  computeOverlappingPartitions();

  // 5.- mark as computed
  IsComputed_ = true;
}

//==============================================================================
template<class GraphType>
void OverlappingPartitioner<GraphType>::computeOverlappingPartitions()
{
  using std::vector;

  // Old FIXME from Ifpack: the first part of this function should be elsewhere
 
  // start defining the subgraphs for no overlap

  vector<size_t> sizes;
  sizes.resize(NumLocalParts_);

  // 1.- compute how many rows are in each subgraph
  for (int i = 0 ; i < NumLocalParts_ ; ++i)
    sizes[i] = 0;

  for (size_t i = 0 ; i < Graph_->getNodeNumRows() ; ++i) {
    TEUCHOS_TEST_FOR_EXCEPTION( (Partition_[i] >= NumLocalParts_),
				std::runtime_error, "Ifpack2::OverlappingPartitioner::computeOverlappingPartitions() Partition_[i] > NumLocalParts_");

    // no singletons should be here, as the matrix is
    // supposed to be filtered through Ifpack2_SingletonFilter
    TEUCHOS_TEST_FOR_EXCEPTION((Partition_[i] == -1),
			       std::runtime_error, "Ifpack2::OverlappingPartitioner::computeOverlappingPartitions() singletons are forbidden");
    sizes[Partition_[i]]++;
  }

  // 2.- allocate space for each subgraph
  for (int i = 0 ; i < NumLocalParts_ ; ++i)
    Parts_[i].resize(sizes[i]);

  // 3.- cycle over all rows and populate the vectors
  for (int i = 0 ; i < NumLocalParts_ ; ++i)
    sizes[i] = 0;

  for (size_t i = 0 ; i < Graph_->getNodeNumRows() ; ++i) {
    LocalOrdinal part = Partition_[i];
    size_t count = sizes[part];
    Parts_[part][count] = i;
    sizes[part]++;
  }

  // If there is no overlap, we're done, so return
  if (OverlappingLevel_ == 0) return;

  // wider overlap requires further computations
  for (size_t level = 1 ; level <= OverlappingLevel_ ; ++level) {
    vector<vector<size_t> > tmp;
    tmp.resize(NumLocalParts_);

    // cycle over all rows in the local graph (that is the overlapping
    // graph). For each row, all columns will belong to the subgraph of
    // row `i'.

    int MaxNumEntries_tmp = Graph_->getNodeMaxNumRowEntries();
    Teuchos::Array<LocalOrdinal> Indices;
    Indices.resize(MaxNumEntries_tmp);

    for (int part = 0 ; part < NumLocalParts_ ; ++part) {
      for (size_t i = 0; i < (size_t)Parts_[part].size() ; ++i) {  
	LocalOrdinal LRID = Parts_[part][i];
	size_t NumIndices;
	Graph_->getLocalRowCopy(LRID,Indices(),NumIndices);

	for (size_t j = 0 ; j < NumIndices ; ++j) {
	  // use *local* indices only
	  LocalOrdinal col = Indices[j];
          if ((size_t)col >= Graph_->getNodeNumRows())
            continue;

	  // has this column already been inserted?
	  vector<size_t>::iterator
	    where = find(tmp[part].begin(), tmp[part].end(), col);

	  if (where == tmp[part].end()) {
	    tmp[part].push_back(col);
	  }
	}
      }
    }

    // now I convert the STL vectors into Teuchos Array RCP's
    for (int i = 0 ; i < NumLocalParts_ ; ++i) {
      Parts_[i].resize(tmp[i].size());
      for (size_t j = 0 ; j < tmp[i].size() ; ++j)
	Parts_[i][j] = tmp[i][j];
    }
  }
}

//==============================================================================
template<class GraphType>
bool OverlappingPartitioner<GraphType>::isComputed()
{
  return(IsComputed_);
}

//==============================================================================
template<class GraphType>
std::ostream& OverlappingPartitioner<GraphType>::print(std::ostream& os) const
{
  using std::endl;
  // NTS: Need to fix for FancyIO later

  if (Graph_->getComm()->getRank()) 
    return(os);

  os << "================================================================================" << endl;
  os << "Ifpack2::OverlappingPartitioner" << endl;
  os << "Number of local rows  = " << Graph_->getNodeNumRows() << endl;
  os << "Number of global rows = " << Graph_->getGlobalNumRows() << endl;
  os << "Number of local parts = " << NumLocalParts_ << endl;
  os << "Overlapping level     = " << OverlappingLevel_ << endl;
  os << "Is computed           = " << IsComputed_ << endl;
  os << "================================================================================" << endl;

  return(os);
}

//==============================================================================
template<class GraphType>
std::string OverlappingPartitioner<GraphType>::description() const
{
  std::ostringstream oss;
  oss << Teuchos::Describable::description();
  if (isComputed()) {
    oss << "{status = computed";
  }
  else {
    oss << "{status = is not computed";
  }
  oss <<"}";
  return oss.str();
}

//==============================================================================
template<class GraphType>
void  OverlappingPartitioner<GraphType>::describe(Teuchos::FancyOStream &os, const Teuchos::EVerbosityLevel verbLevel) const
{
  using std::endl;
  if(verbLevel==Teuchos::VERB_NONE) return;

  os << "================================================================================" << endl;
  os << "Ifpack2::OverlappingPartitioner" << endl;
  os << "Number of local rows  = " << Graph_->getNodeNumRows() << endl;
  os << "Number of global rows = " << Graph_->getGlobalNumRows() << endl;
  os << "Number of local parts = " << NumLocalParts_ << endl;
  os << "Overlapping level     = " << OverlappingLevel_ << endl;
  os << "Is computed           = " << IsComputed_ << endl;
  os << "================================================================================" << endl;
}


}// namespace Ifpack2

#endif // IFPACK2_OVERLAPPINGPARTITIONER_DEF_HPP
