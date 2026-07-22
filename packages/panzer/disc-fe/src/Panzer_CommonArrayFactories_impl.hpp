// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
  typedef PHX::KokkosViewFactory<Scalar,typename PHX::DevLayout<Scalar>::type,PHX::Device> ViewFactory;

  PHX::MDField<Scalar> field = PHX::MDField<Scalar>(prefix_+str,Teuchos::rcp(new PHX::MDALayout<T0>(d0))); 

  if(allocArray_)
    field.setFieldData(ViewFactory::buildView(field.fieldTag(),ddims_));

  return field;
}

template <typename Scalar,typename T0,typename T1>
PHX::MDField<Scalar> MDFieldArrayFactory::
buildArray(const std::string & str,int d0,int d1) const
{ 
  typedef PHX::KokkosViewFactory<Scalar,typename PHX::DevLayout<Scalar>::type,PHX::Device> ViewFactory;

  PHX::MDField<Scalar> field = PHX::MDField<Scalar>(prefix_+str,Teuchos::rcp(new PHX::MDALayout<T0,T1 >(d0,d1))); 

  if(allocArray_)
    field.setFieldData(ViewFactory::buildView(field.fieldTag(),ddims_));

  return field;
}

template <typename Scalar,typename T0,typename T1,typename T2>
PHX::MDField<Scalar> MDFieldArrayFactory::
buildArray(const std::string & str,int d0,int d1,int d2) const
{ 
  typedef PHX::KokkosViewFactory<Scalar,typename PHX::DevLayout<Scalar>::type,PHX::Device> ViewFactory;

  PHX::MDField<Scalar> field = PHX::MDField<Scalar>(prefix_+str,Teuchos::rcp(new PHX::MDALayout<T0,T1,T2>(d0,d1,d2))); 

  if(allocArray_)
    field.setFieldData(ViewFactory::buildView(field.fieldTag(),ddims_));

  return field;
}

template <typename Scalar,typename T0,typename T1,typename T2,typename T3>
PHX::MDField<Scalar> MDFieldArrayFactory::
buildArray(const std::string & str,int d0,int d1,int d2,int d3) const
{ 
  typedef PHX::KokkosViewFactory<Scalar,typename PHX::DevLayout<Scalar>::type,PHX::Device> ViewFactory;

  PHX::MDField<Scalar> field = PHX::MDField<Scalar>(prefix_+str,Teuchos::rcp(new PHX::MDALayout<T0,T1,T2,T3>(d0,d1,d2,d3))); 

  if(allocArray_)
    field.setFieldData(ViewFactory::buildView(field.fieldTag(),ddims_));

  return field;
}

template <typename Scalar,typename T0,typename T1,typename T2,typename T3,typename T4>
PHX::MDField<Scalar> MDFieldArrayFactory::
buildArray(const std::string & str,int d0,int d1,int d2,int d3,int d4) const
{ 
  typedef PHX::KokkosViewFactory<Scalar,typename PHX::DevLayout<Scalar>::type,PHX::Device> ViewFactory;

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
  typedef PHX::KokkosViewFactory<Scalar,typename PHX::DevLayout<Scalar>::type,PHX::Device> ViewFactory;

  PHX::MDField<Scalar,T0> field = PHX::MDField<Scalar,T0>(prefix_+str,Teuchos::rcp(new PHX::MDALayout<T0>(d0))); 

  if(allocArray_)
    field.setFieldData(ViewFactory::buildView(field.fieldTag(),ddims_));

  return field;
}

template <typename Scalar,typename T0,typename T1>
PHX::MDField<Scalar,T0,T1> MDFieldArrayFactory::
buildStaticArray(const std::string & str,int d0,int d1) const
{ 
  typedef PHX::KokkosViewFactory<Scalar,typename PHX::DevLayout<Scalar>::type,PHX::Device> ViewFactory;

  PHX::MDField<Scalar,T0,T1> field = PHX::MDField<Scalar,T0,T1>(prefix_+str,Teuchos::rcp(new PHX::MDALayout<T0,T1>(d0,d1))); 

  if(allocArray_)
    field.setFieldData(ViewFactory::buildView(field.fieldTag(),ddims_));

  return field;
}

template <typename Scalar,typename T0,typename T1,typename T2>
PHX::MDField<Scalar,T0,T1,T2> MDFieldArrayFactory::
buildStaticArray(const std::string & str,int d0,int d1,int d2) const
{ 
  typedef PHX::KokkosViewFactory<Scalar,typename PHX::DevLayout<Scalar>::type,PHX::Device> ViewFactory;

  PHX::MDField<Scalar,T0,T1,T2> field = PHX::MDField<Scalar,T0,T1,T2>(prefix_+str,Teuchos::rcp(new PHX::MDALayout<T0,T1,T2>(d0,d1,d2))); 

  if(allocArray_)
    field.setFieldData(ViewFactory::buildView(field.fieldTag(),ddims_));

  return field;
}

template <typename Scalar,typename T0,typename T1,typename T2,typename T3>
PHX::MDField<Scalar,T0,T1,T2,T3> MDFieldArrayFactory::
buildStaticArray(const std::string & str,int d0,int d1,int d2,int d3) const
{ 
  typedef PHX::KokkosViewFactory<Scalar,typename PHX::DevLayout<Scalar>::type,PHX::Device> ViewFactory;

  PHX::MDField<Scalar,T0,T1,T2,T3> field = PHX::MDField<Scalar,T0,T1,T2,T3>(prefix_+str,Teuchos::rcp(new PHX::MDALayout<T0,T1,T2,T3>(d0,d1,d2,d3))); 

  if(allocArray_)
    field.setFieldData(ViewFactory::buildView(field.fieldTag(),ddims_));

  return field;
}

template <typename Scalar,typename T0,typename T1,typename T2,typename T3,typename T4>
PHX::MDField<Scalar,T0,T1,T2,T3,T4> MDFieldArrayFactory::
buildStaticArray(const std::string & str,int d0,int d1,int d2,int d3,int d4) const
{ 
  typedef PHX::KokkosViewFactory<Scalar,typename PHX::DevLayout<Scalar>::type,PHX::Device> ViewFactory;

  PHX::MDField<Scalar,T0,T1,T2,T3,T4> field = PHX::MDField<Scalar,T0,T1,T2,T3,T4>(prefix_+str,Teuchos::rcp(new PHX::MDALayout<T0,T1,T2,T3,T4>(d0,d1,d2,d3,d4))); 

  if(allocArray_)
    field.setFieldData(ViewFactory::buildView(field.fieldTag(),ddims_));

  return field;
}

} // end namespace panzer

#endif
