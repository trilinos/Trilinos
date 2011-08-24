//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2009) Sandia Corporation
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
// ************************************************************************
//@HEADER

/// \file Tsqr_ConfigDefs.hpp
/// \brief File to include in order to get TSQR's configure-time options.
///

#ifndef __Tsqr_ConfigDefs_hpp
#define __Tsqr_ConfigDefs_hpp

#include <Tsqr_Config.hpp>

/// \namespace TSQR
/// \brief Implementation of the Tall Skinny QR (TSQR) factorization.
///
/// This namespace contains a full hybrid-parallel (MPI + Kokkos)
/// implementation of the Tall Skinny QR (TSQR) factorization.  The
/// following paper describes the implementation:
///
/// Mark Hoemmen.  "A communication-avoiding, hybrid-parallel,
/// rank-revealing orthogonalization method."  IEEE International
/// Parallel and Distributed Processing Symposium (IPDPS), April 2011.
///
/// For further details, see the following:
///
/// Marghoob Mohiyuddin, Mark Hoemmen, James Demmel, and Kathy Yelick.
/// "Minimizing Communication in Sparse Matrix Solvers."  In
/// Proceedings of Supercomputing 2009, November 2009.
///
/// James Demmel, Laura Grigori, Mark Frederick Hoemmen, and Julien
/// Langou.  "Communication-optimal parallel and sequential QR and LU
/// factorizations."  Technical report, UCB/EECS-2008-89, August 2008.
///
namespace TSQR {
  //
  // We declare the TSQR namespace here so that Doxygen will find it
  // and pull in all its documentation.
  //


  /// \namespace Test
  /// \brief Accuracy and performance tests for TSQR.
  ///
  namespace Test {
    //
    // We declare the TSQR::Test namespace here so that Doxygen will
    // find it and pull in all its documentation.
    //
  } // namespace Test

} // namespace TSQR


#endif // __Tsqr_ConfigDefs_hpp
