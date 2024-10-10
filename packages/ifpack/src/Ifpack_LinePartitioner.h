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

#ifndef IFPACK_LINEPARTITIONER_H
#define IFPACK_LINEPARTITIONER_H

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Partitioner.h"
#include "Ifpack_OverlappingPartitioner.h"
#include "Ifpack_Graph_Epetra_RowMatrix.h"
#include "Teuchos_ParameterList.hpp"
class Epetra_Comm;
class Ifpack_Graph;
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_Import;
class Epetra_RowMatrix;

/* \brief Ifpack_LinePartitioner: A class to define partitions into a set of lines.

These "line" partitions could then be used in to do block Gauss-Seidel smoothing, for instance.

The current implementation uses a local (on processor) line detection in one of two forms.
Both are inspired by the work of Mavriplis for convection-diffusion (AIAA Journal, Vol 37(10), 1999).  


Algorithm 1: Matrix entries

Here we use the matrix entries, running lines across the "large" matrix entries, if they are sufficiently
large relative to the small ones.


Algorithms 2: Coordinates

Here we use coordinate information to pick "close" points if they are sufficiently far away
from the "far" points.  We also make sure the line can never double back on itself, so that
the associated sub-matrix could (in theory) be handed off to a fast triangular solver.  This 
implementation doesn't actual do that, however.

This implementation is deived from the related routine in ML.


Supported parameters:
  \c "partitioner: line type": "matrix entries" or "coordinates" (char *)
  \c "partitioner: line detection threshold": if ||x_j - x_i||^2 < thresh * max_k||x_k - x_i||^2, then the points are close enough to line smooth (double)
  \c "partitioner: x-coordinates": x coordinates of local nodes (double *)
  \c "partitioner: y-coordinates": y coordinates of local nodes (double *)
  \c "partitioner: z-coordinates": z coordinates of local nodes (double *)
  \c "partitioner: PDE equations": number of equations per node (integer)

*/

class Ifpack_LinePartitioner : public Ifpack_OverlappingPartitioner {

public:
  // Useful typedef
  typedef enum {COORDINATES=0, MATRIX_ENTRIES,} LINE_MODE;


  //! Constructor.
  Ifpack_LinePartitioner(const Ifpack_Graph* Graph) :
  Ifpack_OverlappingPartitioner(Graph),
    Matrix_(0),
    mode_(COORDINATES),
    NumEqns_(1),
    xcoord_(0),
    ycoord_(0),
    zcoord_(0),
    threshold_(0.0)    
  {
	
  }

 Ifpack_LinePartitioner(const Epetra_RowMatrix* Matrix) :
  Ifpack_OverlappingPartitioner(0),
    Matrix_(Matrix),
    mode_(COORDINATES),
    NumEqns_(1),
    xcoord_(0),
    ycoord_(0),
    zcoord_(0),
    threshold_(0.0)    
  {    
    GraphWrapper_ = Teuchos::rcp(new Ifpack_Graph_Epetra_RowMatrix(Teuchos::rcp(Matrix,false)));
    Graph_ = &*GraphWrapper_;
	    
  }


  //! Destructor.
  virtual ~Ifpack_LinePartitioner() {}

  //! Sets all the parameters for the partitioner
  int SetPartitionParameters(Teuchos::ParameterList& List)
  {
    std::string mymode;
    mode_=COORDINATES;
    mymode = List.get("partitioner: line mode",mymode);
    if(mymode=="coordinates")         mode_=COORDINATES;
    else if(mymode=="matrix entries") mode_=MATRIX_ENTRIES;

    threshold_ = List.get("partitioner: line detection threshold",threshold_);
    if(threshold_ < 0.0)  IFPACK_CHK_ERR(-1);
    if(mode_==COORDINATES && threshold_ > 1.0)  IFPACK_CHK_ERR(-1);


    NumEqns_   = List.get("partitioner: PDE equations",NumEqns_);
    if(NumEqns_ < 1 )  IFPACK_CHK_ERR(-2);
    
    xcoord_   = List.get("partitioner: x-coordinates",xcoord_);
    ycoord_   = List.get("partitioner: y-coordinates",ycoord_);
    zcoord_   = List.get("partitioner: z-coordinates",zcoord_);
    if(mode_==COORDINATES && !xcoord_ && !ycoord_ && !zcoord_) IFPACK_CHK_ERR(-3);

    return(0);
  }

  //! Computes the partitions. Returns 0 if successful.
  int ComputePartitions();

private:

  // Useful functions
  int Compute_Blocks_AutoLine(int * blockIndices) const;
  void local_automatic_line_search(int NumEqns, int * blockIndices, int last, int next, int LineID, double tol, int *itemp, double * dtemp) const;

  // Stuff I allocated
  Teuchos::RCP<const Ifpack_Graph> GraphWrapper_;

  // User data
  const Epetra_RowMatrix* Matrix_;
  LINE_MODE mode_;
  int NumEqns_;
  double * xcoord_;
  double * ycoord_;
  double * zcoord_;
  double threshold_;


  // State data
};

#endif // IFPACK_LINEPARTITIONER_H
