/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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

template<class GraphType>
OverlappingPartitioner<GraphType>::
OverlappingPartitioner (const Teuchos::RCP<const row_graph_type>& graph) :
  NumLocalParts_ (1),
  Graph_ (graph),
  OverlappingLevel_ (0),
  IsComputed_ (false),
  verbose_ (false)
{}


template<class GraphType>
OverlappingPartitioner<GraphType>::~OverlappingPartitioner() {}


template<class GraphType>
int
OverlappingPartitioner<GraphType>::numLocalParts () const 
{
  return NumLocalParts_;
}


template<class GraphType>
int OverlappingPartitioner<GraphType>::overlappingLevel () const
{
  return OverlappingLevel_;
}


template<class GraphType>
typename GraphType::local_ordinal_type 
OverlappingPartitioner<GraphType>::
operator () (const local_ordinal_type MyRow) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    MyRow < 0 || Teuchos::as<size_t> (MyRow) > Graph_->getNodeNumRows (), 
    std::runtime_error, 
    "Ifpack2::OverlappingPartitioner::operator(): "
    "Invalid local row index " << MyRow << ".");
   
  return Partition_[MyRow];
}


//==============================================================================
template<class GraphType>
typename GraphType::local_ordinal_type 
OverlappingPartitioner<GraphType>::
operator() (const local_ordinal_type i, const local_ordinal_type j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    i < 0 || i > Teuchos::as<local_ordinal_type> (NumLocalParts_), 
    std::runtime_error, 
    "Ifpack2::OverlappingPartitioner::operator(): "
    "Invalid local row index i=" << i << ".");
  TEUCHOS_TEST_FOR_EXCEPTION(
    j < 0 || j > Teuchos::as<local_ordinal_type> (Parts_[i].size ()),
    std::runtime_error, 
    "Ifpack2::OverlappingPartitioner::operator(): "
    "Invalid node index j=" << j << ".");
  return Parts_[i][j];
}

//==============================================================================
template<class GraphType>
size_t
OverlappingPartitioner<GraphType>::
numRowsInPart (const local_ordinal_type Part) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    Part < 0 || Part > Teuchos::as<local_ordinal_type> (NumLocalParts_), 
    std::runtime_error, 
    "Ifpack2::OverlappingPartitioner::numRowsInPart: "
    "Invalid partition index Part=" << Part << ".");
  return Parts_[Part].size ();
}

//==============================================================================
template<class GraphType>
void
OverlappingPartitioner<GraphType>::
rowsInPart (const local_ordinal_type Part,
	    Teuchos::ArrayRCP<local_ordinal_type>& List) const
{
  // Let numRowsInPart do the sanity checking...
  const size_t numRows = numRowsInPart (Part); 
  for (size_t i = 0; i < numRows; ++i) {
    List[i] = Parts_[Part][i];
  }
}

//==============================================================================
template<class GraphType>
Teuchos::ArrayView<const typename GraphType::local_ordinal_type>
OverlappingPartitioner<GraphType>::nonOverlappingPartition () const
{
  return Partition_.view (0, Graph_->getNodeNumRows ());
}

//==============================================================================
template<class GraphType>
void
OverlappingPartitioner<GraphType>::
setParameters (Teuchos::ParameterList& List)
{
  NumLocalParts_    = List.get("partitioner: local parts", NumLocalParts_);
  OverlappingLevel_ = List.get("partitioner: overlap", OverlappingLevel_);
  verbose_          = List.get("partitioner: print level", verbose_);

  if (NumLocalParts_ < 0) {
    NumLocalParts_ = Graph_->getNodeNumRows() / (-NumLocalParts_);
  }
  if (NumLocalParts_ == 0) {
    NumLocalParts_ = 1;
  }
  
  // Sanity checking
  TEUCHOS_TEST_FOR_EXCEPTION(
    NumLocalParts_ < 0 || 
    Teuchos::as<size_t> (NumLocalParts_) > Graph_->getNodeNumRows(),
    std::runtime_error, 
    "Ifpack2::OverlappingPartitioner::setParameters: "
    "Invalid NumLocalParts_ = " << NumLocalParts_ << ".");
  TEUCHOS_TEST_FOR_EXCEPTION(
    OverlappingLevel_ < 0, std::runtime_error,
    "Ifpack2::OverlappingPartitioner::setParameters: "
    "Invalid OverlappingLevel_ = " << OverlappingLevel_ << ".");

  setPartitionParameters(List);
}

//==============================================================================
template<class GraphType>
void OverlappingPartitioner<GraphType>::compute()
{
  using std::cout;
  using std::endl;

  TEUCHOS_TEST_FOR_EXCEPTION(
    NumLocalParts_ < 1 || OverlappingLevel_ < 0,
    std::runtime_error, 
    "Ifpack2::OverlappingPartitioner::compute: "
    "Invalid NumLocalParts_ or OverlappingLevel_.");

  // std::string's constructor has some overhead, so it's better to
  // use const char[] for local constant strings.
  const char printMsg[] = "OverlappingPartitioner: ";

  if (verbose_ && (Graph_->getComm()->getRank() == 0)) {
    cout << printMsg << "Number of local parts          = " 
	 << NumLocalParts_ << endl;
    cout << printMsg << "Approx. Number of global parts = " 
	 << NumLocalParts_ * Graph_->getComm ()->getSize () << endl;
    cout << printMsg << "Amount of overlap              = " 
	 << OverlappingLevel_ << endl;
  }

  // 1.- allocate memory 
  Partition_.resize (Graph_->getNodeNumRows ());
  Parts_.resize (NumLocalParts_);

  // 2.- sanity checks on input graph
  TEUCHOS_TEST_FOR_EXCEPTION( 
    ! Graph_->isFillComplete (), std::runtime_error, 
    "Ifpack2::OverlappingPartitioner::compute: "
    "The input graph must be fill complete.");

  TEUCHOS_TEST_FOR_EXCEPTION( 
    Graph_->getGlobalNumRows () != Graph_->getGlobalNumCols (),
    std::runtime_error, 
    "Ifpack2::OverlappingPartitioner::compute: "
    "The input graph must be (globally) square.");
 
  // 3.- perform non-overlapping partition 
  computePartitions ();

  // 4.- compute the partitions with overlapping
  computeOverlappingPartitions ();

  // 5.- mark as computed
  IsComputed_ = true;
}

//==============================================================================
template<class GraphType>
void OverlappingPartitioner<GraphType>::computeOverlappingPartitions()
{
  const local_ordinal_type invalid = 
    Teuchos::OrdinalTraits<local_ordinal_type>::invalid();

  // Old FIXME from Ifpack: the first part of this function should be elsewhere
 
  // start defining the subgraphs for no overlap

  std::vector<size_t> sizes;
  sizes.resize (NumLocalParts_);

  // 1.- compute how many rows are in each subgraph
  for (int i = 0; i < NumLocalParts_; ++i) {
    sizes[i] = 0;
  }

  for (size_t i = 0; i < Graph_->getNodeNumRows (); ++i) {
    TEUCHOS_TEST_FOR_EXCEPTION(
      Partition_[i] >= NumLocalParts_, std::runtime_error, 
      "Ifpack2::OverlappingPartitioner::computeOverlappingPartitions: "
      "Partition_[i] > NumLocalParts_.");
    // invalid indicates that this unknown is not in a nonoverlapping
    // partition
    if (Partition_[i] != invalid) {
      sizes[Partition_[i]]++;
    }
  }

  // 2.- allocate space for each subgraph
  for (int i = 0; i < NumLocalParts_; ++i) {
    Parts_[i].resize (sizes[i]);
  }

  // 3.- cycle over all rows and populate the vectors
  for (int i = 0; i < NumLocalParts_; ++i) {
    sizes[i] = 0;
  }

  for (size_t i = 0; i < Graph_->getNodeNumRows (); ++i) {
    const local_ordinal_type part = Partition_[i];
    if (part != invalid) {
      const size_t count = sizes[part];
      Parts_[part][count] = i;
      sizes[part]++;
    }
  }

  // If there is no overlap, we're done, so return
  if (OverlappingLevel_ == 0) {
    return;
  }

  // wider overlap requires further computations
  for (int level = 1; level <= OverlappingLevel_; ++level) {
    std::vector<std::vector<size_t> > tmp;
    tmp.resize (NumLocalParts_);

    // cycle over all rows in the local graph (that is the overlapping
    // graph). For each row, all columns will belong to the subgraph
    // of row `i'.

    int MaxNumEntries_tmp = Graph_->getNodeMaxNumRowEntries();
    Teuchos::Array<local_ordinal_type> Indices;
    Indices.resize (MaxNumEntries_tmp);

    for (int part = 0; part < NumLocalParts_ ; ++part) {
      for (size_t i = 0; i < Teuchos::as<size_t> (Parts_[part].size ()); ++i) {  
	const local_ordinal_type LRID = Parts_[part][i];
        
	size_t NumIndices;
	Graph_->getLocalRowCopy (LRID, Indices (), NumIndices);

	for (size_t j = 0; j < NumIndices; ++j) {
	  // use *local* indices only
	  const local_ordinal_type col = Indices[j];
          if (Teuchos::as<size_t> (col) >= Graph_->getNodeNumRows ()) {
            continue;
	  }

	  // has this column already been inserted?
	  std::vector<size_t>::iterator where = 
	    std::find (tmp[part].begin (), tmp[part].end (), Teuchos::as<size_t> (col));

	  if (where == tmp[part].end()) {
	    tmp[part].push_back (col);
	  }
	}

        // has this column already been inserted?
	std::vector<size_t>::iterator where =
          std::find (tmp[part].begin (), tmp[part].end (), Teuchos::as<size_t> (LRID));
        
        // This happens here b/c Vanka on Stokes with Stabilized elements will have
        // a zero pivot entry if this gets pushed back first. So... Last.
        if (where == tmp[part].end ()) {
          tmp[part].push_back (LRID);
        }
      }
    }

    // now I convert the STL vectors into Teuchos Array RCP's
    //
    // FIXME (mfh 12 July 2013) You could have started with ArrayRCP
    // in the first place (which implements push_back and iterators)
    // and avoided the copy.
    for (int i = 0; i < NumLocalParts_; ++i) {
      Parts_[i].resize (tmp[i].size ());
      for (size_t j = 0; j < tmp[i].size (); ++j) {
	Parts_[i][j] = tmp[i][j];
      }
    }
  }
}

//==============================================================================
template<class GraphType>
bool OverlappingPartitioner<GraphType>::isComputed() const
{
  return IsComputed_;
}

//==============================================================================
template<class GraphType>
std::ostream& 
OverlappingPartitioner<GraphType>::print (std::ostream& os) const
{
  Teuchos::FancyOStream fos (Teuchos::rcpFromRef (os));
  fos.setOutputToRootOnly (0);
  describe (fos);
  return os;
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
  if (verbLevel == Teuchos::VERB_NONE) {
    return;
  }

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

#define IFPACK2_OVERLAPPINGPARTITIONER_INSTANT(LO,GO,N) \
  template class Ifpack2::OverlappingPartitioner<Tpetra::CrsGraph< LO, GO, N > >; \
  template class Ifpack2::OverlappingPartitioner<Tpetra::RowGraph< LO, GO, N > >;

#endif // IFPACK2_OVERLAPPINGPARTITIONER_DEF_HPP
