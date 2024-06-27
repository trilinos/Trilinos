// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_LINE_PARTITIONER_DEF_HPP
#define IFPACK2_LINE_PARTITIONER_DEF_HPP

#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Util.hpp"

namespace Ifpack2 {

template<class T>
inline typename Teuchos::ScalarTraits<T>::magnitudeType square(T x) {
  return Teuchos::ScalarTraits<T>::magnitude(x) * Teuchos::ScalarTraits<T>::magnitude(x);
}

//==============================================================================
// Constructor
template<class GraphType,class Scalar>
LinePartitioner<GraphType,Scalar>::
LinePartitioner (const Teuchos::RCP<const row_graph_type>& graph) :
  OverlappingPartitioner<GraphType> (graph), NumEqns_(1), threshold_(0.0)
{}


template<class GraphType,class Scalar>
LinePartitioner<GraphType,Scalar>::~LinePartitioner() {}


template<class GraphType,class Scalar>
void
LinePartitioner<GraphType,Scalar>::
setPartitionParameters(Teuchos::ParameterList& List) {
  threshold_ = List.get("partitioner: line detection threshold",threshold_);
  TEUCHOS_TEST_FOR_EXCEPTION(threshold_ < 0.0 || threshold_ > 1.0,
                             std::runtime_error,"Ifpack2::LinePartitioner: threshold not valid");

  NumEqns_   = List.get("partitioner: PDE equations",NumEqns_);
  TEUCHOS_TEST_FOR_EXCEPTION(NumEqns_<1,std::runtime_error,"Ifpack2::LinePartitioner: NumEqns not valid");

  coord_   = List.get("partitioner: coordinates",coord_);
  TEUCHOS_TEST_FOR_EXCEPTION(coord_.is_null(),std::runtime_error,"Ifpack2::LinePartitioner: coordinates not defined");
}


template<class GraphType,class Scalar>
void LinePartitioner<GraphType,Scalar>::computePartitions() {
  const local_ordinal_type invalid = Teuchos::OrdinalTraits<local_ordinal_type>::invalid();

  // Sanity Checks
  TEUCHOS_TEST_FOR_EXCEPTION(coord_.is_null(),std::runtime_error,"Ifpack2::LinePartitioner: coordinates not defined");
  TEUCHOS_TEST_FOR_EXCEPTION((size_t)this->Partition_.size() != this->Graph_->getLocalNumRows(),std::runtime_error,"Ifpack2::LinePartitioner: partition size error");

  // Short circuit
  if(this->Partition_.size() == 0) {this->NumLocalParts_ = 0; return;}

  // Set partitions to invalid to initialize algorithm
  for(size_t i=0; i<this->Graph_->getLocalNumRows(); i++)
    this->Partition_[i] = invalid;

  // Use the auto partitioner
  this->NumLocalParts_ = this->Compute_Blocks_AutoLine(this->Partition_());

  // Resize Parts_
  this->Parts_.resize(this->NumLocalParts_);
}


// ============================================================================
template<class GraphType,class Scalar>
int LinePartitioner<GraphType,Scalar>::Compute_Blocks_AutoLine(Teuchos::ArrayView<local_ordinal_type> blockIndices) const {
  typedef local_ordinal_type LO;
  typedef magnitude_type MT;
  const LO invalid  = Teuchos::OrdinalTraits<LO>::invalid();
  const MT zero = Teuchos::ScalarTraits<MT>::zero();

  Teuchos::ArrayRCP<const MT>  xvalsRCP, yvalsRCP, zvalsRCP;
  Teuchos::ArrayView<const MT> xvals, yvals, zvals;
  xvalsRCP = coord_->getData(0); xvals = xvalsRCP();
  if(coord_->getNumVectors() > 1) { yvalsRCP = coord_->getData(1); yvals = yvalsRCP(); }
  if(coord_->getNumVectors() > 2) { zvalsRCP = coord_->getData(2); zvals = zvalsRCP(); }

  double tol             = threshold_;
  size_t N               = this->Graph_->getLocalNumRows();
  size_t allocated_space = this->Graph_->getLocalMaxNumRowEntries();

  nonconst_local_inds_host_view_type cols("cols",allocated_space);
  Teuchos::Array<LO>     indices(allocated_space);
  Teuchos::Array<MT> dist(allocated_space);

  Teuchos::Array<LO>     itemp(2*allocated_space);
  Teuchos::Array<MT> dtemp(allocated_space);

  LO num_lines = 0;

  for(LO i=0; i<(LO)N; i+=NumEqns_) {
    size_t nz=0;
    // Short circuit if I've already been blocked
    if(blockIndices[i] != invalid) continue;

    // Get neighbors and sort by distance
    this->Graph_->getLocalRowCopy(i,cols,nz);
    MT x0 = (!xvals.is_null()) ? xvals[i/NumEqns_] : zero;
    MT y0 = (!yvals.is_null()) ? yvals[i/NumEqns_] : zero;
    MT z0 = (!zvals.is_null()) ? zvals[i/NumEqns_] : zero;

    LO neighbor_len=0;
    for(size_t j=0; j<nz; j+=NumEqns_) {
      MT mydist = zero;
      LO nn = cols[j] / NumEqns_;
      if(cols[j] >=(LO)N) continue; // Check for off-proc entries
      if(!xvals.is_null()) mydist += square<MT>(x0 - xvals[nn]);
      if(!yvals.is_null()) mydist += square<MT>(y0 - yvals[nn]);
      if(!zvals.is_null()) mydist += square<MT>(z0 - zvals[nn]);
      dist[neighbor_len] = Teuchos::ScalarTraits<MT>::squareroot(mydist);
      indices[neighbor_len]=cols[j];
      neighbor_len++;
    }

    Teuchos::ArrayView<MT> dist_view = dist(0,neighbor_len);
    Tpetra::sort2(dist_view.begin(),dist_view.end(),indices.begin());

    // Number myself
    for(LO k=0; k<NumEqns_; k++)
      blockIndices[i + k] = num_lines;

    // Fire off a neighbor line search (nearest neighbor)
    if(neighbor_len > 2 && dist[1]/dist[neighbor_len-1] < tol) {
      local_automatic_line_search(NumEqns_,blockIndices,i,indices[1],num_lines,tol,itemp,dtemp);
    }
    // Fire off a neighbor line search (second nearest neighbor)
    if(neighbor_len > 3 && dist[2]/dist[neighbor_len-1] < tol) {
      local_automatic_line_search(NumEqns_,blockIndices,i,indices[2],num_lines,tol,itemp,dtemp);
    }

    num_lines++;
  }
  return num_lines;
}
// ============================================================================
template<class GraphType,class Scalar>
void LinePartitioner<GraphType,Scalar>::local_automatic_line_search(int NumEqns, Teuchos::ArrayView <local_ordinal_type> blockIndices, local_ordinal_type last, local_ordinal_type next,  local_ordinal_type LineID, double tol,  Teuchos::Array<local_ordinal_type> itemp, Teuchos::Array<magnitude_type> dtemp) const {
  typedef local_ordinal_type LO;
  typedef magnitude_type MT;
  const LO invalid  = Teuchos::OrdinalTraits<LO>::invalid();
  const MT zero = Teuchos::ScalarTraits<MT>::zero();

  Teuchos::ArrayRCP<const MT>  xvalsRCP, yvalsRCP, zvalsRCP;
  Teuchos::ArrayView<const MT> xvals, yvals, zvals;
  xvalsRCP = coord_->getData(0); xvals = xvalsRCP();
  if(coord_->getNumVectors() > 1) { yvalsRCP = coord_->getData(1); yvals = yvalsRCP(); }
  if(coord_->getNumVectors() > 2) { zvalsRCP = coord_->getData(2); zvals = zvalsRCP(); }

  size_t N               = this->Graph_->getLocalNumRows();
  size_t allocated_space = this->Graph_->getLocalMaxNumRowEntries();

  nonconst_local_inds_host_view_type cols(itemp.data(),allocated_space);
  Teuchos::ArrayView<LO>     indices = itemp.view(allocated_space,allocated_space);
  Teuchos::ArrayView<MT> dist= dtemp();

  while (blockIndices[next] == invalid) {
    // Get the next row
    size_t nz=0;
    LO neighbors_in_line=0;

    this->Graph_->getLocalRowCopy(next,cols,nz);
    MT x0 = (!xvals.is_null()) ? xvals[next/NumEqns_] : zero;
    MT y0 = (!yvals.is_null()) ? yvals[next/NumEqns_] : zero;
    MT z0 = (!zvals.is_null()) ? zvals[next/NumEqns_] : zero;

    // Calculate neighbor distances & sort
    LO neighbor_len=0;
    for(size_t i=0; i<nz; i+=NumEqns) {
      MT mydist = zero;
      if(cols[i] >=(LO)N) continue; // Check for off-proc entries
      LO nn = cols[i] / NumEqns;
      if(blockIndices[nn]==LineID) neighbors_in_line++;
      if(!xvals.is_null()) mydist += square<MT>(x0 - xvals[nn]);
      if(!yvals.is_null()) mydist += square<MT>(y0 - yvals[nn]);
      if(!zvals.is_null()) mydist += square<MT>(z0 - zvals[nn]);
      dist[neighbor_len] = Teuchos::ScalarTraits<MT>::squareroot(mydist);
      indices[neighbor_len]=cols[i];
      neighbor_len++;
    }
    // If more than one of my neighbors is already in this line.  I
    // can't be because I'd create a cycle
    if(neighbors_in_line > 1) break;

    // Otherwise add me to the line
    for(LO k=0; k<NumEqns; k++)
      blockIndices[next + k] = LineID;

    // Try to find the next guy in the line (only check the closest two that aren't element 0 (diagonal))
    Teuchos::ArrayView<MT> dist_view = dist(0,neighbor_len);
    Tpetra::sort2(dist_view.begin(),dist_view.end(),indices.begin());

    if(neighbor_len > 2 && indices[1] != last && blockIndices[indices[1]] == -1 && dist[1]/dist[neighbor_len-1] < tol) {
      last=next;
      next=indices[1];
    }
    else if(neighbor_len > 3 && indices[2] != last && blockIndices[indices[2]] == -1 && dist[2]/dist[neighbor_len-1] < tol) {
      last=next;
      next=indices[2];
    }
    else {
      // I have no further neighbors in this line
      break;
    }
  }
}



}// namespace Ifpack2

#define IFPACK2_LINEPARTITIONER_INSTANT(S,LO,GO,N) \
  template class Ifpack2::LinePartitioner<Tpetra::CrsGraph< LO, GO, N >,S >; \
  template class Ifpack2::LinePartitioner<Tpetra::RowGraph< LO, GO, N >,S >;

#endif // IFPACK2_LINEPARTITIONER_DEF_HPP
