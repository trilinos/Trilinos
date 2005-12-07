// @HEADER
// ***********************************************************************
//
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
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
// @HEADER

#include "Amesos_Status.h"
void Amesos_Status::SetStatusParameters( const Teuchos::ParameterList &ParameterList) {

  // some verbose output:
  // 0 - no output at all
  // 1 - output as specified by other parameters
  // 2 - all possible output
  if( ParameterList.isParameter("OutputLevel") )
    verbose_ = ParameterList.get<int>("OutputLevel");

  // level of debug output:
  // 0 - no output at all
  // 1 - some debug output - set by some tests upon a test failure 
  // >1 - more debug output (unused at this point)
  if( ParameterList.isParameter("DebugLevel") )
    debug_ = ParameterList.get<int>("DebugLevel");

  // print some timing information (on process 0)
  if( ParameterList.isParameter("PrintTiming") )
    PrintTiming_ = ParameterList.get<bool>("PrintTiming");

  // print some statistics (on process 0). Do not include timing
  if( ParameterList.isParameter("PrintStatus") )
    PrintStatus_ = ParameterList.get<bool>("PrintStatus");

  // compute norms of some vectors
  if( ParameterList.isParameter("ComputeVectorNorms") )
    ComputeVectorNorms_ = ParameterList.get<bool>("ComputeVectorNorms");

  // compute the true residual Ax-b after solution
  if( ParameterList.isParameter("ComputeTrueResidual") )
    ComputeTrueResidual_ = ParameterList.get<bool>("ComputeTrueResidual");

}
