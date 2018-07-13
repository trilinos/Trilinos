// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef PANZER_COMMON_ARRAY_FACTORIES_IMPL_HPP
#define PANZER_COMMON_ARRAY_FACTORIES_IMPL_HPP

#include "Teuchos_RCP.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_KokkosViewFactory.hpp"

namespace panzer {

// Implementation for intrepid container factory
template <typename Scalar,typename T0>
Kokkos::DynRankView<Scalar,PHX::Device> Intrepid2FieldContainerFactory::
buildArray(const std::string & str,int d0) const
{
  static_assert(std::is_same<Scalar,double>::value,"ERROR: CommonArryFactory for DynRankView only supports double scalar type!");
  return Kokkos::DynRankView<Scalar,PHX::Device>(str,d0);
}

template <typename Scalar,typename T0,typename T1>
Kokkos::DynRankView<Scalar,PHX::Device> Intrepid2FieldContainerFactory::
buildArray(const std::string & str,int d0,int d1) const
{ 
  static_assert(std::is_same<Scalar,double>::value,"ERROR: CommonArryFactory for DynRankView only supports double scalar type!");
  return Kokkos::DynRankView<Scalar,PHX::Device>(str,d0,d1);
}

template <typename Scalar,typename T0,typename T1,typename T2>
Kokkos::DynRankView<Scalar,PHX::Device> Intrepid2FieldContainerFactory::
buildArray(const std::string & str,int d0,int d1,int d2) const
{
  static_assert(std::is_same<Scalar,double>::value,"ERROR: CommonArryFactory for DynRankView only supports double scalar type!");
  return Kokkos::DynRankView<Scalar,PHX::Device>(str,d0,d1,d2);
}

template <typename Scalar,typename T0,typename T1,typename T2,typename T3>
Kokkos::DynRankView<Scalar,PHX::Device> Intrepid2FieldContainerFactory::
buildArray(const std::string & str,int d0,int d1,int d2,int d3) const
{
  static_assert(std::is_same<Scalar,double>::value,"ERROR: CommonArryFactory for DynRankView only supports double scalar type!");
  return Kokkos::DynRankView<Scalar,PHX::Device>(str,d0,d1,d2,d3);
}

template <typename Scalar,typename T0,typename T1,typename T2,typename T3,typename T4>
Kokkos::DynRankView<Scalar,PHX::Device> Intrepid2FieldContainerFactory::
buildArray(const std::string & str,int d0,int d1,int d2,int d3,int d4) const
{
  static_assert(std::is_same<Scalar,double>::value,"ERROR: CommonArryFactory for DynRankView only supports double scalar type!");
  return Kokkos::DynRankView<Scalar,PHX::Device>(str,d0,d1,d2,d3,d4);
}

// Implementation for MDField array factory
template <typename Scalar,typename T0>
PHX::MDField<Scalar> MDFieldArrayFactory::
buildArray(const std::string & str,int d0) const
{ 
  typedef PHX::KokkosViewFactory<Scalar,PHX::Device> ViewFactory;

  PHX::MDField<Scalar> field = PHX::MDField<Scalar>(prefix_+str,Teuchos::rcp(new PHX::MDALayout<T0>(d0))); 

  if(allocArray_)
    field.setFieldData(ViewFactory::buildView(field.fieldTag(),ddims_));

  return field;
}

template <typename Scalar,typename T0,typename T1>
PHX::MDField<Scalar> MDFieldArrayFactory::
buildArray(const std::string & str,int d0,int d1) const
{ 
  typedef PHX::KokkosViewFactory<Scalar,PHX::Device> ViewFactory;

  PHX::MDField<Scalar> field = PHX::MDField<Scalar>(prefix_+str,Teuchos::rcp(new PHX::MDALayout<T0,T1 >(d0,d1))); 

  if(allocArray_)
    field.setFieldData(ViewFactory::buildView(field.fieldTag(),ddims_));

  return field;
}

template <typename Scalar,typename T0,typename T1,typename T2>
PHX::MDField<Scalar> MDFieldArrayFactory::
buildArray(const std::string & str,int d0,int d1,int d2) const
{ 
  typedef PHX::KokkosViewFactory<Scalar,PHX::Device> ViewFactory;

  PHX::MDField<Scalar> field = PHX::MDField<Scalar>(prefix_+str,Teuchos::rcp(new PHX::MDALayout<T0,T1,T2>(d0,d1,d2))); 

  if(allocArray_)
    field.setFieldData(ViewFactory::buildView(field.fieldTag(),ddims_));

  return field;
}

template <typename Scalar,typename T0,typename T1,typename T2,typename T3>
PHX::MDField<Scalar> MDFieldArrayFactory::
buildArray(const std::string & str,int d0,int d1,int d2,int d3) const
{ 
  typedef PHX::KokkosViewFactory<Scalar,PHX::Device> ViewFactory;

  PHX::MDField<Scalar> field = PHX::MDField<Scalar>(prefix_+str,Teuchos::rcp(new PHX::MDALayout<T0,T1,T2,T3>(d0,d1,d2,d3))); 

  if(allocArray_)
    field.setFieldData(ViewFactory::buildView(field.fieldTag(),ddims_));

  return field;
}

template <typename Scalar,typename T0,typename T1,typename T2,typename T3,typename T4>
PHX::MDField<Scalar> MDFieldArrayFactory::
buildArray(const std::string & str,int d0,int d1,int d2,int d3,int d4) const
{ 
  typedef PHX::KokkosViewFactory<Scalar,PHX::Device> ViewFactory;

  PHX::MDField<Scalar> field = PHX::MDField<Scalar>(prefix_+str,Teuchos::rcp(new PHX::MDALayout<T0,T1,T2,T3,T4>(d0,d1,d2,d3,d4))); 

  if(allocArray_)
    field.setFieldData(ViewFactory::buildView(field.fieldTag(),ddims_));

  return field;
}

// Implementation for MDField array factory
template <typename Scalar,typename T0>
PHX::MDField<Scalar,T0> MDFieldArrayFactory::
buildStaticArray(const std::string & str,int d0) const
{ 
  typedef PHX::KokkosViewFactory<Scalar,PHX::Device> ViewFactory;

  PHX::MDField<Scalar,T0> field = PHX::MDField<Scalar,T0>(prefix_+str,Teuchos::rcp(new PHX::MDALayout<T0>(d0))); 

  if(allocArray_)
    field.setFieldData(ViewFactory::buildView(field.fieldTag(),ddims_));

  return field;
}

template <typename Scalar,typename T0,typename T1>
PHX::MDField<Scalar,T0,T1> MDFieldArrayFactory::
buildStaticArray(const std::string & str,int d0,int d1) const
{ 
  typedef PHX::KokkosViewFactory<Scalar,PHX::Device> ViewFactory;

  PHX::MDField<Scalar,T0,T1> field = PHX::MDField<Scalar,T0,T1>(prefix_+str,Teuchos::rcp(new PHX::MDALayout<T0,T1>(d0,d1))); 

  if(allocArray_)
    field.setFieldData(ViewFactory::buildView(field.fieldTag(),ddims_));

  return field;
}

template <typename Scalar,typename T0,typename T1,typename T2>
PHX::MDField<Scalar,T0,T1,T2> MDFieldArrayFactory::
buildStaticArray(const std::string & str,int d0,int d1,int d2) const
{ 
  typedef PHX::KokkosViewFactory<Scalar,PHX::Device> ViewFactory;

  PHX::MDField<Scalar,T0,T1,T2> field = PHX::MDField<Scalar,T0,T1,T2>(prefix_+str,Teuchos::rcp(new PHX::MDALayout<T0,T1,T2>(d0,d1,d2))); 

  if(allocArray_)
    field.setFieldData(ViewFactory::buildView(field.fieldTag(),ddims_));

  return field;
}

template <typename Scalar,typename T0,typename T1,typename T2,typename T3>
PHX::MDField<Scalar,T0,T1,T2,T3> MDFieldArrayFactory::
buildStaticArray(const std::string & str,int d0,int d1,int d2,int d3) const
{ 
  typedef PHX::KokkosViewFactory<Scalar,PHX::Device> ViewFactory;

  PHX::MDField<Scalar,T0,T1,T2,T3> field = PHX::MDField<Scalar,T0,T1,T2,T3>(prefix_+str,Teuchos::rcp(new PHX::MDALayout<T0,T1,T2,T3>(d0,d1,d2,d3))); 

  if(allocArray_)
    field.setFieldData(ViewFactory::buildView(field.fieldTag(),ddims_));

  return field;
}

template <typename Scalar,typename T0,typename T1,typename T2,typename T3,typename T4>
PHX::MDField<Scalar,T0,T1,T2,T3,T4> MDFieldArrayFactory::
buildStaticArray(const std::string & str,int d0,int d1,int d2,int d3,int d4) const
{ 
  typedef PHX::KokkosViewFactory<Scalar,PHX::Device> ViewFactory;

  PHX::MDField<Scalar,T0,T1,T2,T3,T4> field = PHX::MDField<Scalar,T0,T1,T2,T3,T4>(prefix_+str,Teuchos::rcp(new PHX::MDALayout<T0,T1,T2,T3,T4>(d0,d1,d2,d3,d4))); 

  if(allocArray_)
    field.setFieldData(ViewFactory::buildView(field.fieldTag(),ddims_));

  return field;
}

} // end namespace panzer

#endif
