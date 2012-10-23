/*
//@HEADER
// ************************************************************************
// 
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

template< typename Scalar , class DeviceType >
struct initialize_node;

template<typename Scalar>
struct initialize_node<Scalar, KOKKOSARRAY_MACRO_DEVICE>
{
	typedef KOKKOSARRAY_MACRO_DEVICE     device_type ;
	typedef typename KokkosArray::MDArray<Scalar,device_type> array_type ;
  typedef typename KokkosArray::MDArray<int,device_type>    int_array_type ;

  enum { NumNodePerElement = 8 };

  typedef Region<Scalar,device_type> MyRegion;

  initialize_node( const MyRegion & arg_region )
    : region(arg_region)
  {}


  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( int inode )const {

    const int begin = region.node_elem_offset(inode);
    const int end = region.node_elem_offset(inode+1);

    Scalar node_mass = 0;
    for(int i = begin; i != end; ++i) {
      node_mass += region.elem_mass(region.node_elem_ids(i,0)) / NumNodePerElement;
    }

    region.nodal_mass(inode) = node_mass;

  }

  const MyRegion region;
};
