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
#include <KokkosArray_View.hpp>
namespace KokkosArray {
namespace Impl {



struct PhysicalLayout {
  enum LayoutType {Left,Right,Scalar,Vector};
  LayoutType layout_type;
  long long int stride[8]; //distance between two neighboring elements in a given dimension
  int rank;

  template<class DataType, class Device, class DataManagement, class Specialisation>
  PhysicalLayout(View<DataType,LayoutLeft,Device,DataManagement,Specialisation> view) {
	layout_type = Left;
	rank = view.Rank;
    for(int i=0;i<8;i++) stride[i] = 0;
	stride[0] = 1;
	stride[1] = view.m_stride;
    for(int i = 2;i<rank;i++)
    	stride[i] = view.dimension(i-1)*stride[i-1];
  }

  template<class DataType, class Device, class DataManagement, class Specialisation>
  PhysicalLayout(View<DataType,LayoutRight,Device,DataManagement,Specialisation> view) {
	layout_type = Right;
	rank = view.Rank;
    for(int i=0;i<8;i++) stride[i] = 0;
	stride[rank-1] = 1;
    for(int i = rank-2;i>=0;i++)
    	stride[i] = view.dimension(i+1)*stride[i+1];
	stride[0] = view.m_stride;
  }

  template<class DataType, class Device, class DataManagement, class Specialisation>
  PhysicalLayout(View<DataType,LayoutScalar,Device,DataManagement,Specialisation> view) {
	layout_type = Scalar;
	rank = 0;
    for(int i=0;i<8;i++) stride[i] = 0;
  }

  template<class DataType, class Device, class DataManagement, class Specialisation>
  PhysicalLayout(View<DataType,LayoutVector,Device,DataManagement,Specialisation> view) {
	layout_type = Vector;
	rank = 1;
    for(int i=0;i<8;i++) stride[i] = 0;
	stride[0] = view.dimension_0();
  }
};

}
}
