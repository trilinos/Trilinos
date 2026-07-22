// $Id$ 
// $Source$ 
// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

template <typename TypeSeq, typename BaseT, typename ObjectT>
PHX::TemplateManager<TypeSeq,BaseT,ObjectT>::
TemplateManager()
{
  // Determine number of types
  int sz = Sacado::mpl::size<TypeSeq>::value;
  objects.resize(sz);
  disabled.resize(sz,false);
}

template <typename TypeSeq, typename BaseT, typename ObjectT>
PHX::TemplateManager<TypeSeq,BaseT,ObjectT>::
~TemplateManager()
{
}

template <typename TypeSeq, typename BaseT, typename ObjectT>
template <typename BuilderOpT>
void
PHX::TemplateManager<TypeSeq,BaseT,ObjectT>::
buildObjects(const BuilderOpT& builder)
{
  Sacado::mpl::for_each_no_kokkos<TypeSeq>(BuildObject<BuilderOpT>(objects,disabled,builder));
}

template <typename TypeSeq, typename BaseT, typename ObjectT>
void
PHX::TemplateManager<TypeSeq,BaseT,ObjectT>::
buildObjects()
{
  DefaultBuilderOp builder;
  (*this).template buildObjects<DefaultBuilderOp>(builder);
}

template <typename TypeSeq, typename BaseT, typename ObjectT>
template<typename ScalarT>
Teuchos::RCP<BaseT>
PHX::TemplateManager<TypeSeq,BaseT,ObjectT>::
getAsBase()
{
  int idx = Sacado::mpl::find<TypeSeq,ScalarT>::value;
  return objects[idx];
}

template <typename TypeSeq, typename BaseT, typename ObjectT>
template<typename ScalarT>
Teuchos::RCP<const BaseT>
PHX::TemplateManager<TypeSeq,BaseT,ObjectT>::getAsBase() const
{
  int idx = Sacado::mpl::find<TypeSeq,ScalarT>::value;
  return objects[idx];
}

template <typename TypeSeq, typename BaseT, typename ObjectT>
template<typename ScalarT>
Teuchos::RCP< typename Sacado::mpl::apply<ObjectT,ScalarT>::type >
PHX::TemplateManager<TypeSeq,BaseT,ObjectT>::
getAsObject()
{
  int idx = Sacado::mpl::find<TypeSeq,ScalarT>::value;
  return Teuchos::rcp_dynamic_cast< typename Sacado::mpl::apply<ObjectT,ScalarT>::type >(objects[idx], true);
}

template <typename TypeSeq, typename BaseT, typename ObjectT>
template<typename ScalarT>
Teuchos::RCP< const typename Sacado::mpl::apply<ObjectT,ScalarT>::type >
PHX::TemplateManager<TypeSeq,BaseT,ObjectT>::
getAsObject() const
{
  int idx = Sacado::mpl::find<TypeSeq,ScalarT>::value;
  return Teuchos::rcp_dynamic_cast< const typename Sacado::mpl::apply<ObjectT,ScalarT>::type >(objects[idx], 
							     true);
}

template <typename TypeSeq, typename BaseT, typename ObjectT>
typename PHX::TemplateManager<TypeSeq,BaseT,ObjectT>::iterator
PHX::TemplateManager<TypeSeq,BaseT,ObjectT>::
begin()
{
  return PHX::TemplateIterator<TypeSeq,BaseT,ObjectT>(*this,
							 objects.begin());
}

template <typename TypeSeq, typename BaseT, typename ObjectT>
typename PHX::TemplateManager<TypeSeq,BaseT,ObjectT>::const_iterator
PHX::TemplateManager<TypeSeq,BaseT,ObjectT>::
begin() const
{
  return PHX::ConstTemplateIterator<TypeSeq,BaseT,ObjectT>(*this,
							      objects.begin());
}

template <typename TypeSeq, typename BaseT, typename ObjectT>
typename PHX::TemplateManager<TypeSeq,BaseT,ObjectT>::iterator
PHX::TemplateManager<TypeSeq,BaseT,ObjectT>::
end()
{
  return PHX::TemplateIterator<TypeSeq,BaseT,ObjectT>(*this,
							 objects.end());
}

template <typename TypeSeq, typename BaseT, typename ObjectT>
typename PHX::TemplateManager<TypeSeq,BaseT,ObjectT>::const_iterator
PHX::TemplateManager<TypeSeq,BaseT,ObjectT>::
end() const
{
  return PHX::ConstTemplateIterator<TypeSeq,BaseT,ObjectT>(*this,
							      objects.end());
}

template <typename TypeSeq, typename BaseT, typename ObjectT>
template<typename ScalarT>
void
PHX::TemplateManager<TypeSeq,BaseT,ObjectT>::deleteType()
{
  int idx = Sacado::mpl::find<TypeSeq,ScalarT>::value;
  objects[idx] = Teuchos::null;
}

template <typename TypeSeq, typename BaseT, typename ObjectT>
template<typename ScalarT>
void
PHX::TemplateManager<TypeSeq,BaseT,ObjectT>::disableType()
{
  int idx = Sacado::mpl::find<TypeSeq,ScalarT>::value;
  disabled[idx] = true;
}
