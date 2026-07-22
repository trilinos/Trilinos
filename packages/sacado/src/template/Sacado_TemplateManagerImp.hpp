// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

template <typename TypeSeq, typename BaseT, typename ObjectT>
Sacado::TemplateManager<TypeSeq,BaseT,ObjectT>::
TemplateManager()
{
  // Determine number of types
  int sz = mpl::size<TypeSeq>::value;
  objects.resize(sz);
}

template <typename TypeSeq, typename BaseT, typename ObjectT>
Sacado::TemplateManager<TypeSeq,BaseT,ObjectT>::
~TemplateManager()
{
}

template <typename TypeSeq, typename BaseT, typename ObjectT>
template <typename BuilderOpT>
void
Sacado::TemplateManager<TypeSeq,BaseT,ObjectT>::
buildObjects(const BuilderOpT& builder)
{
  mpl::for_each<TypeSeq>(BuildObject<BuilderOpT>(objects,builder));
}

template <typename TypeSeq, typename BaseT, typename ObjectT>
void
Sacado::TemplateManager<TypeSeq,BaseT,ObjectT>::
buildObjects()
{
  DefaultBuilderOp builder;
  (*this).template buildObjects<DefaultBuilderOp>(builder);
}

template <typename TypeSeq, typename BaseT, typename ObjectT>
template<typename ScalarT>
Teuchos::RCP<BaseT>
Sacado::TemplateManager<TypeSeq,BaseT,ObjectT>::
getAsBase()
{
  int idx = mpl::find<TypeSeq,ScalarT>::value;
  return objects[idx];
}

template <typename TypeSeq, typename BaseT, typename ObjectT>
template<typename ScalarT>
Teuchos::RCP<const BaseT>
Sacado::TemplateManager<TypeSeq,BaseT,ObjectT>::getAsBase() const
{
  int idx = mpl::find<TypeSeq,ScalarT>::value;
  return objects[idx];
}

template <typename TypeSeq, typename BaseT, typename ObjectT>
template<typename ScalarT>
Teuchos::RCP< typename Sacado::mpl::apply<ObjectT,ScalarT>::type >
Sacado::TemplateManager<TypeSeq,BaseT,ObjectT>::
getAsObject()
{
  int idx = mpl::find<TypeSeq,ScalarT>::value;
  return Teuchos::rcp_dynamic_cast< typename Sacado::mpl::apply<ObjectT,ScalarT>::type >(objects[idx], true);
}

template <typename TypeSeq, typename BaseT, typename ObjectT>
template<typename ScalarT>
Teuchos::RCP< const typename Sacado::mpl::apply<ObjectT,ScalarT>::type >
Sacado::TemplateManager<TypeSeq,BaseT,ObjectT>::
getAsObject() const
{
  int idx = mpl::find<TypeSeq,ScalarT>::value;
  return Teuchos::rcp_dynamic_cast< const typename Sacado::mpl::apply<ObjectT,ScalarT>::type >(objects[idx], true);
}

template <typename TypeSeq, typename BaseT, typename ObjectT>
typename Sacado::TemplateManager<TypeSeq,BaseT,ObjectT>::iterator
Sacado::TemplateManager<TypeSeq,BaseT,ObjectT>::
begin()
{
  return Sacado::TemplateIterator<BaseT>(objects.begin());
}

template <typename TypeSeq, typename BaseT, typename ObjectT>
typename Sacado::TemplateManager<TypeSeq,BaseT,ObjectT>::const_iterator
Sacado::TemplateManager<TypeSeq,BaseT,ObjectT>::
begin() const
{
  return Sacado::ConstTemplateIterator<BaseT>(objects.begin());
}

template <typename TypeSeq, typename BaseT, typename ObjectT>
typename Sacado::TemplateManager<TypeSeq,BaseT,ObjectT>::iterator
Sacado::TemplateManager<TypeSeq,BaseT,ObjectT>::
end()
{
  return Sacado::TemplateIterator<BaseT>(objects.end());
}

template <typename TypeSeq, typename BaseT, typename ObjectT>
typename Sacado::TemplateManager<TypeSeq,BaseT,ObjectT>::const_iterator
Sacado::TemplateManager<TypeSeq,BaseT,ObjectT>::
end() const
{
  return Sacado::ConstTemplateIterator<BaseT>(objects.end());
}
