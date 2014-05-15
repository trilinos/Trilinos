// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
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
