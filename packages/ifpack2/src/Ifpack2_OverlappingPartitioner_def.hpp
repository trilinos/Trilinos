// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_OVERLAPPINGPARTITIONER_DEF_HPP
#define IFPACK2_OVERLAPPINGPARTITIONER_DEF_HPP
#include <vector>
#include <string>
#include <algorithm>
#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_OverlappingPartitioner_decl.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"

namespace Ifpack2 {

template<class GraphType>
OverlappingPartitioner<GraphType>::
OverlappingPartitioner (const Teuchos::RCP<const row_graph_type>& graph) :
  NumLocalParts_ (1),
  Graph_ (graph),
  OverlappingLevel_ (0),
  IsComputed_ (false),
  verbose_ (false),
  maintainSparsity_(false)
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
    MyRow < 0 || Teuchos::as<size_t> (MyRow) > Graph_->getLocalNumRows (), 
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
  return Partition_.view (0, Graph_->getLocalNumRows ());
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
  maintainSparsity_ = List.get("partitioner: maintain sparsity", false);
  typedef Teuchos::RCP<   Tpetra::Map<typename GraphType::local_ordinal_type, typename GraphType::global_ordinal_type, typename GraphType::node_type> const >  map_type;
  typedef Teuchos::RCP<Tpetra::Import<typename GraphType::local_ordinal_type, typename GraphType::global_ordinal_type, typename GraphType::node_type> const >  import_type;

  // when using overlapping schwarz wth combineMode ADD, we need import and
  // overlap row map to adjust weights for ith dof. Specifically, we sum 
  // all blocks (including off processor ones) that contain i.
  import_type theImport = Graph_->getImporter();
  if (theImport != Teuchos::null)   List.set< import_type >("theImport",theImport);
  List.set< map_type >("OverlapRowMap",Graph_->getRowMap());

  if (NumLocalParts_ < 0) {
    NumLocalParts_ = Graph_->getLocalNumRows() / (-NumLocalParts_);
  }
  if (NumLocalParts_ == 0) {
    NumLocalParts_ = 1;
  }
  
  // Sanity checking
  TEUCHOS_TEST_FOR_EXCEPTION(
    NumLocalParts_ < 0 || 
    Teuchos::as<size_t> (NumLocalParts_) > Graph_->getLocalNumRows(),
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
  Partition_.resize (Graph_->getLocalNumRows ());
  //Parts_ is allocated in computeOverlappingPartitions_, where it is used

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
  //If user has explicitly specified parts, then Partition_ size has been set to 0.
  //In this case, there is no need to compute Parts_.
  if (Partition_.size() == 0)
    return;

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

  for (size_t i = 0; i < Graph_->getLocalNumRows (); ++i) {
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
  Parts_.resize (NumLocalParts_);
  for (int i = 0; i < NumLocalParts_; ++i) {
    Parts_[i].resize (sizes[i]);
  }

  // 3.- cycle over all rows and populate the vectors
  for (int i = 0; i < NumLocalParts_; ++i) {
    sizes[i] = 0;
  }

  for (size_t i = 0; i < Graph_->getLocalNumRows (); ++i) {
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

    int MaxNumEntries_tmp = Graph_->getLocalMaxNumRowEntries();
    nonconst_local_inds_host_view_type Indices("Indices",MaxNumEntries_tmp);
    nonconst_local_inds_host_view_type newIndices("newIndices",MaxNumEntries_tmp);
    
    if (!maintainSparsity_) {

      local_ordinal_type numLocalRows = Graph_->getLocalNumRows();
      for (int part = 0; part < NumLocalParts_ ; ++part) {
        for (size_t i = 0; i < Teuchos::as<size_t> (Parts_[part].size ()); ++i) {
          const local_ordinal_type LRID = Parts_[part][i];
          
          size_t numIndices;
          Graph_->getLocalRowCopy (LRID, Indices, numIndices);

          for (size_t j = 0; j < numIndices; ++j) {
            // use *local* indices only
            const local_ordinal_type col = Indices[j];
            if (col >= numLocalRows) {
              continue;
            }

            // has this column already been inserted?
            std::vector<size_t>::iterator where = 
              std::find (tmp[part].begin (), tmp[part].end (), Teuchos::as<size_t> (col));


            if (where == tmp[part].end()) {
              tmp[part].push_back (col);
            }
          }

          //10-Jan-2017 JHU : A possible optimization is to avoid the std::find's in the loop above,
          //but instead sort and make unique afterwards.  One would have to be careful to only sort/
          //make unique the insertions associated with the current row LRID, as for each row, LRID
          //may occur after all other col indices for row LRID (see comment about zero pivot below).

          // has this column already been inserted?
          std::vector<size_t>::iterator where =
            std::find (tmp[part].begin (), tmp[part].end (), Teuchos::as<size_t> (LRID));
          
          // This happens here b/c Vanka on Stokes with Stabilized elements will have
          // a zero pivot entry if this gets pushed back first. So... Last.
          if (where == tmp[part].end ()) {
            tmp[part].push_back (LRID);
          }
        }
      } //for (int part = 0; ...

    } else {
      //maintainSparsity_ == true

      for (int part = 0; part < NumLocalParts_ ; ++part) {
        for (size_t i = 0; i < Teuchos::as<size_t> (Parts_[part].size ()); ++i) {
          const local_ordinal_type LRID = Parts_[part][i];
          
          size_t numIndices;
          Graph_->getLocalRowCopy (LRID, Indices, numIndices);
          //JJH: the entries in Indices are already sorted.  However, the Tpetra documentation states
          //     that we can't count on this always being true, hence we sort.  Also note that there are
          //     unused entries at the end of Indices (it's sized to hold any row).  This means we can't
          //     just use Indices.end() in sorting and in std::includes
          Tpetra::sort(Indices,numIndices);

          for (size_t j = 0; j < numIndices; ++j) {
            // use *local* indices only
            const local_ordinal_type col = Indices[j];
            if (Teuchos::as<size_t> (col) >= Graph_->getLocalNumRows ()) {
              continue;
            }

            // has this column already been inserted?
            std::vector<size_t>::iterator where = 
              std::find (tmp[part].begin (), tmp[part].end (), Teuchos::as<size_t> (col));


            if (where == tmp[part].end()) {
              // Check if row associated with "col" increases connectivity already defined by row LRID's stencil.
              // If it does and maintainSparsity_ is true, do not add "col" to the current partition (block).
              size_t numNewIndices;
              Graph_->getLocalRowCopy(col, newIndices, numNewIndices);
              Tpetra::sort(newIndices,numNewIndices);
              auto Indices_rcp = Kokkos::Compat::persistingView<nonconst_local_inds_host_view_type>(Indices, 0, numIndices);
              auto newIndices_rcp = Kokkos::Compat::persistingView<nonconst_local_inds_host_view_type>(newIndices, 0, numNewIndices);
              bool isSubset = std::includes(Indices_rcp.begin(),Indices_rcp.begin()+numIndices,
                                   newIndices_rcp.begin(),newIndices_rcp.begin()+numNewIndices);
              if (isSubset) {
                tmp[part].push_back (col);
              }
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
      } //for (int part = 0;

    } //if (maintainSparsity_) ... else

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
  os << "Number of local rows  = " << Graph_->getLocalNumRows() << endl;
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
