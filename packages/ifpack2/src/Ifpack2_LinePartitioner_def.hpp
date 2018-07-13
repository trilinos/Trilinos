/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
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
  TEUCHOS_TEST_FOR_EXCEPTION((size_t)this->Partition_.size() != this->Graph_->getNodeNumRows(),std::runtime_error,"Ifpack2::LinePartitioner: partition size error");

  // Short circuit
  if(this->Partition_.size() == 0) {this->NumLocalParts_ = 0; return;}

  // Set partitions to invalid to initialize algorithm
  for(size_t i=0; i<this->Graph_->getNodeNumRows(); i++)
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
  const LO invalid  = Teuchos::OrdinalTraits<LO>::invalid();
  const double zero = Teuchos::ScalarTraits<double>::zero();

  Teuchos::ArrayRCP<const double>  xvalsRCP, yvalsRCP, zvalsRCP;
  Teuchos::ArrayView<const double> xvals, yvals, zvals;
  xvalsRCP = coord_->getData(0); xvals = xvalsRCP();
  if(coord_->getNumVectors() > 1) { yvalsRCP = coord_->getData(1); yvals = yvalsRCP(); }
  if(coord_->getNumVectors() > 2) { zvalsRCP = coord_->getData(2); zvals = zvalsRCP(); }

  double tol             = threshold_;
  size_t N               = this->Graph_->getNodeNumRows();
  size_t allocated_space = this->Graph_->getNodeMaxNumRowEntries();

  Teuchos::Array<LO>     cols(allocated_space);
  Teuchos::Array<LO>     indices(allocated_space);
  Teuchos::Array<double> dist(allocated_space);

  Teuchos::Array<LO>     itemp(2*allocated_space);
  Teuchos::Array<double> dtemp(allocated_space);

  LO num_lines = 0;

  for(LO i=0; i<(LO)N; i+=NumEqns_) {
    size_t nz=0;
    // Short circuit if I've already been blocked
    if(blockIndices[i] != invalid) continue;

    // Get neighbors and sort by distance
    this->Graph_->getLocalRowCopy(i,cols(),nz);
    double x0 = (!xvals.is_null()) ? xvals[i/NumEqns_] : zero;
    double y0 = (!yvals.is_null()) ? yvals[i/NumEqns_] : zero;
    double z0 = (!zvals.is_null()) ? zvals[i/NumEqns_] : zero;

    LO neighbor_len=0;
    for(size_t j=0; j<nz; j+=NumEqns_) {
      double mydist = zero;
      LO nn = cols[j] / NumEqns_;
      if(cols[j] >=(LO)N) continue; // Check for off-proc entries
      if(!xvals.is_null()) mydist += square<double>(x0 - xvals[nn]);
      if(!yvals.is_null()) mydist += square<double>(y0 - yvals[nn]);
      if(!zvals.is_null()) mydist += square<double>(z0 - zvals[nn]);
      dist[neighbor_len] = Teuchos::ScalarTraits<double>::squareroot(mydist);
      indices[neighbor_len]=cols[j];
      neighbor_len++;
    }

    Teuchos::ArrayView<double> dist_view = dist(0,neighbor_len);
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
void LinePartitioner<GraphType,Scalar>::local_automatic_line_search(int NumEqns, Teuchos::ArrayView <local_ordinal_type> blockIndices, local_ordinal_type last, local_ordinal_type next,  local_ordinal_type LineID, double tol,  Teuchos::Array<local_ordinal_type> itemp, Teuchos::Array<double> dtemp) const {
  typedef local_ordinal_type LO;
  const LO invalid  = Teuchos::OrdinalTraits<LO>::invalid();
  const double zero = Teuchos::ScalarTraits<double>::zero();

  Teuchos::ArrayRCP<const double>  xvalsRCP, yvalsRCP, zvalsRCP;
  Teuchos::ArrayView<const double> xvals, yvals, zvals;
  xvalsRCP = coord_->getData(0); xvals = xvalsRCP();
  if(coord_->getNumVectors() > 1) { yvalsRCP = coord_->getData(1); yvals = yvalsRCP(); }
  if(coord_->getNumVectors() > 2) { zvalsRCP = coord_->getData(2); zvals = zvalsRCP(); }

  size_t N               = this->Graph_->getNodeNumRows();
  size_t allocated_space = this->Graph_->getNodeMaxNumRowEntries();
  Teuchos::ArrayView<LO>     cols    = itemp();
  Teuchos::ArrayView<LO>     indices = itemp.view(allocated_space,allocated_space);
  Teuchos::ArrayView<double> dist= dtemp();

  while (blockIndices[next] == invalid) {
    // Get the next row
    size_t nz=0;
    LO neighbors_in_line=0;

    this->Graph_->getLocalRowCopy(next,cols(),nz);
    double x0 = (!xvals.is_null()) ? xvals[next/NumEqns_] : zero;
    double y0 = (!yvals.is_null()) ? yvals[next/NumEqns_] : zero;
    double z0 = (!zvals.is_null()) ? zvals[next/NumEqns_] : zero;

    // Calculate neighbor distances & sort
    LO neighbor_len=0;
    for(size_t i=0; i<nz; i+=NumEqns) {
      double mydist = zero;
      if(cols[i] >=(LO)N) continue; // Check for off-proc entries
      LO nn = cols[i] / NumEqns;
      if(blockIndices[nn]==LineID) neighbors_in_line++;
      if(!xvals.is_null()) mydist += square<double>(x0 - xvals[nn]);
      if(!yvals.is_null()) mydist += square<double>(y0 - yvals[nn]);
      if(!zvals.is_null()) mydist += square<double>(z0 - zvals[nn]);
      dist[neighbor_len] = Teuchos::ScalarTraits<double>::squareroot(mydist);
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
    Teuchos::ArrayView<double> dist_view = dist(0,neighbor_len);
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
