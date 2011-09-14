/** \HEADER
 *************************************************************************
 *
 *                            Kokkos
 *                 Copyright 2010 Sandia Corporation
 *
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the Corporation nor the names of the
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 *  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 *  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *************************************************************************
 */

template<class Scalar , class DeviceType >
struct Dirichlet;

template<class Scalar >
struct Dirichlet<Scalar , KOKKOS_MACRO_DEVICE >
{
  
  typedef KOKKOS_MACRO_DEVICE                               device_type;
  typedef device_type::size_type                            index_type;
  typedef Kokkos::MultiVectorView<Scalar, device_type>      scalar_vector;  
  typedef Kokkos::MultiVectorView<index_type, device_type>  index_vector;

  scalar_vector A ;
  index_vector  A_row ;
  index_vector  A_col ;
  scalar_vector b ;

  index_vector  row_flag ;
  scalar_vector value ;

  Dirichlet(
    const scalar_vector & arg_A, 
    const index_vector  & arg_A_row,
    const index_vector  & arg_A_col, 
    const scalar_vector & arg_b, 
    const index_vector  & arg_row_flag,
    const scalar_vector & arg_value )
  : A(     arg_A )
  , A_row( arg_A_row )
  , A_col( arg_A_col )
  , b(     arg_b )
  , row_flag( arg_row_flag )
  , value(    arg_value )
  { }
    

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( index_type irow ) const
  {
    //  Apply a dirichlet boundary condition to 'irow'
    //  to maintain the symmetry of the original 
    //  global stiffness matrix, zero out the columns
    //  that correspond to boundary conditions, and
    //  adjust the load vector accordingly

    const index_type iBeg = A_row(irow);
    const index_type iEnd = A_row(irow+1);

    if ( row_flag( irow ) ) {

      //  set the load vector equal to a specified value

      b(irow) = value( irow );

      //  zero each value on the row, and leave a one
      //  on the diagonal

      for(index_type i = iBeg ; i < iEnd ; i++){
        A(i) = irow == A_col(i) ? 1 : 0 ;
      }
    }
    else {
      //  Find any columns that are boundary conditions.
      //  Clear them and adjust the load vector

      for( index_type i = iBeg ; i < iEnd ; i++ ){
        const index_type j = A_col(i) ;
        if ( row_flag(j) ) {
          b( irow ) -= value(j) * A(i);
          A( i ) = 0;
        }
      }
    }
  }
};



