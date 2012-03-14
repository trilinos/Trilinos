#ifndef PANZER_COMMON_ARRAY_FACTORIES_IMPL_HPP
#define PANZER_COMMON_ARRAY_FACTORIES_IMPL_HPP

#include "Teuchos_RCP.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"

namespace panzer {

// Implementation for intrepid container factory
template <typename Scalar>
template <typename T0>
Intrepid::FieldContainer<Scalar> IntrepidFieldContainerFactory<Scalar>::
buildArray(const std::string & str,int d0) const
{ return Intrepid::FieldContainer<Scalar>(d0); }

template <typename Scalar>
template <typename T0,typename T1>
Intrepid::FieldContainer<Scalar> IntrepidFieldContainerFactory<Scalar>::
buildArray(const std::string & str,int d0,int d1) const
{ return Intrepid::FieldContainer<Scalar>(d0,d1); }

template <typename Scalar>
template <typename T0,typename T1,typename T2>
Intrepid::FieldContainer<Scalar> IntrepidFieldContainerFactory<Scalar>::
buildArray(const std::string & str,int d0,int d1,int d2) const
{ return Intrepid::FieldContainer<Scalar>(d0,d1,d2); }

template <typename Scalar>
template <typename T0,typename T1,typename T2,typename T3>
Intrepid::FieldContainer<Scalar> IntrepidFieldContainerFactory<Scalar>::
buildArray(const std::string & str,int d0,int d1,int d2,int d3) const
{ return Intrepid::FieldContainer<Scalar>(d0,d1,d2,d3); }

template <typename Scalar>
template <typename T0,typename T1,typename T2,typename T3,typename T4>
Intrepid::FieldContainer<Scalar> IntrepidFieldContainerFactory<Scalar>::
buildArray(const std::string & str,int d0,int d1,int d2,int d3,int d4) const
{ return Intrepid::FieldContainer<Scalar>(d0,d1,d2,d3,d4); }

// Implementation for MDField array factory
template <typename Scalar>
MDFieldArrayFactory<Scalar>::
MDFieldArrayFactory() : prefix_("") {}

template <typename Scalar>
MDFieldArrayFactory<Scalar>::
MDFieldArrayFactory(const std::string & prefix) : prefix_(prefix) {}

template <typename Scalar>
template <typename T0>
PHX::MDField<Scalar> MDFieldArrayFactory<Scalar>::
buildArray(const std::string & str,int d0) const
{ return PHX::MDField<Scalar>(prefix_+str,Teuchos::rcp(new PHX::MDALayout<T0>(d0))); }

template <typename Scalar>
template <typename T0,typename T1>
PHX::MDField<Scalar> MDFieldArrayFactory<Scalar>::
buildArray(const std::string & str,int d0,int d1) const
{ return PHX::MDField<Scalar>(prefix_+str,Teuchos::rcp(new PHX::MDALayout<T0,T1 >(d0,d1))); }

template <typename Scalar>
template <typename T0,typename T1,typename T2>
PHX::MDField<Scalar> MDFieldArrayFactory<Scalar>::
buildArray(const std::string & str,int d0,int d1,int d2) const
{ return PHX::MDField<Scalar>(prefix_+str,Teuchos::rcp(new PHX::MDALayout<T0,T1,T2>(d0,d1,d2))); }

template <typename Scalar>
template <typename T0,typename T1,typename T2,typename T3>
PHX::MDField<Scalar> MDFieldArrayFactory<Scalar>::
buildArray(const std::string & str,int d0,int d1,int d2,int d3) const
{ return PHX::MDField<Scalar>(prefix_+str,Teuchos::rcp(new PHX::MDALayout<T0,T1,T2,T3>(d0,d1,d2,d3))); }

template <typename Scalar>
template <typename T0,typename T1,typename T2,typename T3,typename T4>
PHX::MDField<Scalar> MDFieldArrayFactory<Scalar>::
buildArray(const std::string & str,int d0,int d1,int d2,int d3,int d4) const
{ return PHX::MDField<Scalar>(prefix_+str,Teuchos::rcp(new PHX::MDALayout<T0,T1,T2,T3,T4>(d0,d1,d2,d3,d4))); }

} // end namespace panzer

#endif
