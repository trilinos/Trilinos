/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Partitioner.h"
#include "Ifpack_OverlappingPartitioner.h"
#include "Ifpack_LinePartitioner.h"
#include "Ifpack_Graph.h"
#include "Epetra_Util.h"

// ============================================================================
inline void Ifpack_LinePartitioner::local_automatic_line_search(int NumEqns, int * blockIndices, int last, int next, int LineID, double tol, int *itemp, double * dtemp) const {
  double *xvals=xcoord_, *yvals=ycoord_, *zvals=zcoord_;

  int N = NumMyRows();

  int allocated_space = MaxNumEntries();
  int * cols    = itemp;
  int * indices = &itemp[allocated_space];
  double * dist = dtemp;


  while (blockIndices[next] == -1) {
    // Get the next row
    int n=0;
    int neighbors_in_line=0;

    Graph_->ExtractMyRowCopy(next,allocated_space,n,cols);
    double x0 = (xvals) ? xvals[next/NumEqns] : 0.0;
    double y0 = (yvals) ? yvals[next/NumEqns] : 0.0;
    double z0 = (zvals) ? zvals[next/NumEqns] : 0.0;

    // Calculate neighbor distances & sort
    int neighbor_len=0;
    for(int i=0; i<n; i+=NumEqns) {
      double mydist = 0.0;
      if(cols[i] >=N) continue; // Check for off-proc entries
      int nn = cols[i] / NumEqns;
      if(blockIndices[nn]==LineID) neighbors_in_line++;
      if(xvals) mydist += (x0 - xvals[nn]) * (x0 - xvals[nn]);
      if(yvals) mydist += (y0 - yvals[nn]) * (y0 - yvals[nn]);
      if(zvals) mydist += (z0 - zvals[nn]) * (z0 - zvals[nn]);
      dist[neighbor_len] = sqrt(mydist);
      indices[neighbor_len]=cols[i];
      neighbor_len++;
    }
    // If more than one of my neighbors is already in this line.  I
    // can't be because I'd create a cycle
    if(neighbors_in_line > 1) break;

    // Otherwise add me to the line 
    for(int k=0; k<NumEqns; k++) 
      blockIndices[next + k] = LineID;
    
    // Try to find the next guy in the line (only check the closest two that aren't element 0 (diagonal))
    Epetra_Util::Sort(true,neighbor_len,dist,0,0,1,&indices,0,0);

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

// ============================================================================
int Ifpack_LinePartitioner::Compute_Blocks_AutoLine(int * blockIndices) const {
  double *xvals=xcoord_, *yvals=ycoord_, *zvals=zcoord_;
  int NumEqns = NumEqns_;
  double tol = threshold_;
  int N = NumMyRows();
  int allocated_space = MaxNumEntries();
    
  int * cols    = new int[2*allocated_space];
  int * indices = &cols[allocated_space];
  double * dist = new double[allocated_space];

  int * itemp   = new int[2*allocated_space];
  double *dtemp = new double[allocated_space];

  int num_lines = 0;

  for(int i=0; i<N; i+=NumEqns) {
    int nz=0;
    // Short circuit if I've already been blocked
    if(blockIndices[i] !=-1) continue;

    // Get neighbors and sort by distance
    Graph_->ExtractMyRowCopy(i,allocated_space,nz,cols);
    double x0 = (xvals) ? xvals[i/NumEqns] : 0.0;
    double y0 = (yvals) ? yvals[i/NumEqns] : 0.0;
    double z0 = (zvals) ? zvals[i/NumEqns] : 0.0;

    int neighbor_len=0;
    for(int j=0; j<nz; j+=NumEqns) {
      double mydist = 0.0;
      int nn = cols[j] / NumEqns;
      if(cols[j] >=N) continue; // Check for off-proc entries
      if(xvals) mydist += (x0 - xvals[nn]) * (x0 - xvals[nn]);
      if(yvals) mydist += (y0 - yvals[nn]) * (y0 - yvals[nn]);
      if(zvals) mydist += (z0 - zvals[nn]) * (z0 - zvals[nn]);
      dist[neighbor_len] = sqrt(mydist);
      indices[neighbor_len]=cols[j];
      neighbor_len++;
    }

    Epetra_Util::Sort(true,neighbor_len,dist,0,0,1,&indices,0,0);

    // Number myself
    for(int k=0; k<NumEqns; k++)
      blockIndices[i + k] = num_lines;

    // Fire off a neighbor line search (nearest neighbor)
    if(neighbor_len > 2 && dist[1]/dist[neighbor_len-1] < tol) {
      local_automatic_line_search(NumEqns,blockIndices,i,indices[1],num_lines,tol,itemp,dtemp);
    }
    // Fire off a neighbor line search (second nearest neighbor)
    if(neighbor_len > 3 && dist[2]/dist[neighbor_len-1] < tol) {
      local_automatic_line_search(NumEqns,blockIndices,i,indices[2],num_lines,tol,itemp,dtemp);
    }

    num_lines++;
  }
  
  // Cleanup
  delete [] cols;
  delete [] dist;
  delete [] itemp;
  delete [] dtemp;

  return num_lines;
}
//==============================================================================
int Ifpack_LinePartitioner::ComputePartitions()
{
  // Sanity Checks
  if(!xcoord_ && !ycoord_ && !zcoord_)  IFPACK_CHK_ERR(-1);
  if((int)Partition_.size() != NumMyRows())  IFPACK_CHK_ERR(-2);

  // Short circuit
  if(Partition_.size() == 0) {NumLocalParts_ = 0; return 0;}

  // Set partitions to -1 to initialize algorithm
  for(int i=0; i<NumMyRows(); i++)
    Partition_[i] = -1;

  // Use the auto partitioner 
  NumLocalParts_ = Compute_Blocks_AutoLine(&Partition_[0]);
  
  // Resize Parts_
  Parts_.resize(NumLocalParts_);
  return(0);
}
