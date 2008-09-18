// $Id$ 
// $Source$ 
// @HEADER
// ************************************************************************
// 
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                  Copyright (2008) Sandia Corporation
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

template <typename TypeSeq, typename BaseT, typename ObjectT>
PHX::TemplateManager<TypeSeq,BaseT,ObjectT>::
TemplateManager()
{
  // Determine number of types
  int sz = Sacado::mpl::size<TypeSeq>::value;
  objects.resize(sz);
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
  Sacado::mpl::for_each<TypeSeq>(BuildObject<BuilderOpT>(objects,builder));
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
Teuchos::RCP< typename boost::mpl::apply<ObjectT,ScalarT>::type >
PHX::TemplateManager<TypeSeq,BaseT,ObjectT>::
getAsObject()
{
  int idx = Sacado::mpl::find<TypeSeq,ScalarT>::value;
  return Teuchos::rcp_dynamic_cast< typename boost::mpl::apply<ObjectT,ScalarT>::type >(objects[idx], true);
}

template <typename TypeSeq, typename BaseT, typename ObjectT>
template<typename ScalarT>
Teuchos::RCP< const typename boost::mpl::apply<ObjectT,ScalarT>::type >
PHX::TemplateManager<TypeSeq,BaseT,ObjectT>::
getAsObject() const
{
  int idx = Sacado::mpl::find<TypeSeq,ScalarT>::value;
  return Teuchos::rcp_dynamic_cast< const typename boost::mpl::apply<ObjectT,ScalarT>::type >(objects[idx], 
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
