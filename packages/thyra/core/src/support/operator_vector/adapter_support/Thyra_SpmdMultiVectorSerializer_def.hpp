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

#ifndef THYRA_SPMD_MULTI_VECTOR_SERIALIZER_HPP
#define THYRA_SPMD_MULTI_VECTOR_SERIALIZER_HPP

#include "Thyra_SpmdMultiVectorSerializer_decl.hpp"
#include "Thyra_SpmdVectorSpaceDefaultBase.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"

namespace Thyra {

template<class Scalar>
SpmdMultiVectorSerializer<Scalar>::SpmdMultiVectorSerializer(
  const bool  my_binaryMode
  )
  :binaryMode_(my_binaryMode)
{}

template<class Scalar>
bool SpmdMultiVectorSerializer<Scalar>::isCompatible(
  const MultiVectorBase<Scalar> &mv
  ) const
{
  return 0!=dynamic_cast<const SpmdVectorSpaceBase<Scalar>*>(&*mv.range());
}

template<class Scalar>
void SpmdMultiVectorSerializer<Scalar>::serialize(
  const MultiVectorBase<Scalar>& mv, std::ostream& out
  ) const
{
  Teuchos::RCP<const SpmdVectorSpaceBase<Scalar> >
    mpi_vec_spc
    = Teuchos::rcp_dynamic_cast<const SpmdVectorSpaceBase<Scalar> >(mv.range());
  out.precision(std::numeric_limits<Scalar>::digits10+4);
  if( mpi_vec_spc.get() ) {
    // This is a mpi-based vector space so let's just write the local
    // multi-vector elements (row-by-row).
    const Ordinal
      localOffset = mpi_vec_spc->localOffset(),
      localSubDim = mpi_vec_spc->localSubDim();
    const Range1D localRng( localOffset, localOffset+localSubDim-1 ); 
    ConstDetachedMultiVectorView<Scalar> local_mv(mv,localRng,Range1D());
    out << localSubDim << " " << local_mv.numSubCols() << std::endl;
    if( binaryMode() ) {
      // Write column-wise for better cache performance
      for( Ordinal j = 0; j < local_mv.numSubCols(); ++j )
        out.write( reinterpret_cast<const char*>(&local_mv(0,j)), sizeof(Scalar)*localSubDim );
    }
    else {
      // Write row-wise for better readability
      for( Ordinal i = 0; i < localSubDim; ++i ) {
        out << " " << i;
        for( Ordinal j = 0; j < local_mv.numSubCols(); ++j ) {
          out << " " << local_mv(i,j);
        }
        out << std::endl;
      }
    }
  }
  else {
    //  This is a serial (or locally replicated) vector space so
    // just write all of the multi-vector elements here.
    TEST_FOR_EXCEPTION( true, std::logic_error, "Does not handle non-SPMD spaces yet" );
  }
}

template<class Scalar>
void SpmdMultiVectorSerializer<Scalar>::deserialize(
  std::istream& in, MultiVectorBase<Scalar>* mv
  ) const
{
  Teuchos::RCP<const SpmdVectorSpaceBase<Scalar> >
    mpi_vec_spc = Teuchos::rcp_dynamic_cast<const SpmdVectorSpaceBase<Scalar> >(mv->range());
  if( mpi_vec_spc.get() ) {
    // This is a mpi-based vector space so let's just read the local
    // multi-vector elements (row-by-row).
    const Ordinal
      localOffset = mpi_vec_spc->localOffset(),
      localSubDim = mpi_vec_spc->localSubDim();
    const Range1D localRng( localOffset, localOffset+localSubDim-1 ); 
    DetachedMultiVectorView<Scalar> local_mv(*mv,localRng,Range1D());
#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPTION(
      !in, std::logic_error
      ,"Error: The input stream given is empty before any reading has began!\n"
      "If this stream came from a file, then the file may not exist!"
      );
#endif
    Ordinal localSubDim_in;
    in >> localSubDim_in;
#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPTION(
      localSubDim != localSubDim_in, std::logic_error
      , "Error, localSubDim = "<<localSubDim<<" does not match the read in value of "
      "localSubDim_in = "<<localSubDim_in<<"!"
      );
#endif
    Ordinal numSubCols_in;
    in >> numSubCols_in;
#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPTION(
      local_mv.numSubCols() != numSubCols_in, std::logic_error
      , "Error, numSubCols = "<<local_mv.numSubCols()<<" does not match the read in value of "
      "numSubCols_in = "<<numSubCols_in<<"!"
      );
#endif
    // Get rid of extra newline after first line
    in >> std::ws;
    // Get the elements
    if( binaryMode() ) {
      // Column-wise
      for( Ordinal j = 0; j < local_mv.numSubCols(); ++j )
        in.read( reinterpret_cast<char*>(&local_mv(0,j)), sizeof(Scalar)*localSubDim );
    }
    else {
      // Row-wise
      for( Ordinal i = 0; i < localSubDim; ++i ) {
#ifdef TEUCHOS_DEBUG
        TEST_FOR_EXCEPTION( !in, std::logic_error, "Error, premature end of input!"	);
#endif
        Ordinal i_in;
        in >> i_in;
#ifdef TEUCHOS_DEBUG
        TEST_FOR_EXCEPTION(
          i != i_in, std::logic_error
          , "Error, i = "<<i<<" does not match the read in value of "
          "i_in = "<<i_in<<"!"
          );
#endif
        for( Ordinal j = 0; j < local_mv.numSubCols(); ++j ) {
#ifdef TEUCHOS_DEBUG
          TEST_FOR_EXCEPTION(
            !in, std::logic_error
            ,"Error: The input stream ran out at j="<<j<<" before"
            " reaching the promised " << local_mv.numSubCols()
            << " rows of the (multi)vector!"
            );
#endif
          in >> local_mv(i,j);
        }
      }
    }
  }
  else {
    //  This is a serial (or locally replicated) vector space so
    // just read all of the multi-vector elements here.
    TEST_FOR_EXCEPTION( true, std::logic_error, "Does not handle non-SPMD spaces yet" );
  }
}

} // end namespace Thyra

#endif // THYRA_SPMD_MULTI_VECTOR_SERIALIZER_HPP
