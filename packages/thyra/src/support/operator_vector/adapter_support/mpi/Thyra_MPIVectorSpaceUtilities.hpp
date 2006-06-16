// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
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

#ifndef THYRA_MPI_VECTOR_SPACE_UTILITIES_HPP
#define THYRA_MPI_VECTOR_SPACE_UTILITIES_HPP

#include "Thyra_OperatorVectorTypes.hpp"
#ifdef HAVE_MPI
#  include "mpi.h"
#endif

namespace Thyra {

namespace MPIVectorSpaceUtilities {

/** \brief . */
Index computeMapCode( MPI_Comm mpiComm, const Index localSubDim );

/** \brief . */
Index computeLocalOffset( MPI_Comm mpiComm, const Index localSubDim );

/** \brief . */
Index computeGlobalDim( MPI_Comm mpiComm, const Index localSubDim );

/** \brief . */
void broadcast( MPI_Comm mpiComm, const int rootRank, Index* value );

} // namespace MPIVectorSpaceUtiltities

} // namespace Thyra

#endif // THYRA_MPI_VECTOR_SPACE_UTILITIES_HPP
