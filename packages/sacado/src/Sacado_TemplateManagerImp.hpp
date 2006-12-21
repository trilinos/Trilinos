// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

template <typename TypeSeq, typename BaseT, template<typename> class ObjectT>
Sacado::TemplateManager<TypeSeq,BaseT,ObjectT>::
TemplateManager()
{
  // Determine number of types
  int sz = mpl::size<TypeSeq>::value;
  objects.resize(sz);
}

template <typename TypeSeq, typename BaseT, template<typename> class ObjectT>
Sacado::TemplateManager<TypeSeq,BaseT,ObjectT>::
~TemplateManager()
{
}

template <typename TypeSeq, typename BaseT, template<typename> class ObjectT>
template <typename BuilderOpT>
void
Sacado::TemplateManager<TypeSeq,BaseT,ObjectT>::
buildObjects(const BuilderOpT& builder)
{
  mpl::for_each<TypeSeq>(BuildObject<BuilderOpT>(objects,builder));
}

template <typename TypeSeq, typename BaseT, template<typename> class ObjectT>
void
Sacado::TemplateManager<TypeSeq,BaseT,ObjectT>::
buildObjects()
{
  DefaultBuilderOp builder;
  (*this).template buildObjects<DefaultBuilderOp>(builder);
}

template <typename TypeSeq, typename BaseT, template<typename> class ObjectT>
template<typename ScalarT>
BaseT&
Sacado::TemplateManager<TypeSeq,BaseT,ObjectT>::
getAsBase()
{
  int idx = mpl::find<TypeSeq,ScalarT>::value;
  return *(objects[idx]);
}

template <typename TypeSeq, typename BaseT, template<typename> class ObjectT>
template<typename ScalarT>
const BaseT&
Sacado::TemplateManager<TypeSeq,BaseT,ObjectT>::getAsBase() const
{
  int idx = mpl::find<TypeSeq,ScalarT>::value;
  return *(objects[idx]);
}

template <typename TypeSeq, typename BaseT, template<typename> class ObjectT>
template<typename ScalarT>
ObjectT<ScalarT>&
Sacado::TemplateManager<TypeSeq,BaseT,ObjectT>::
getAsObject()
{
  return
    dynamic_cast< ObjectT<ScalarT>& >( (*this).template getAsBase<ScalarT>() );
}

template <typename TypeSeq, typename BaseT, template<typename> class ObjectT>
template<typename ScalarT>
const ObjectT<ScalarT>&
Sacado::TemplateManager<TypeSeq,BaseT,ObjectT>::
getAsObject() const
{
  return
    dynamic_cast< const ObjectT<ScalarT>& >( (*this).template getAsBase<ScalarT>() );
}

template <typename TypeSeq, typename BaseT, template<typename> class ObjectT>
typename Sacado::TemplateManager<TypeSeq,BaseT,ObjectT>::iterator
Sacado::TemplateManager<TypeSeq,BaseT,ObjectT>::
begin()
{
  return Sacado::TemplateIterator<TypeSeq,BaseT,ObjectT>(*this,
							 objects.begin());
}

template <typename TypeSeq, typename BaseT, template<typename> class ObjectT>
typename Sacado::TemplateManager<TypeSeq,BaseT,ObjectT>::iterator
Sacado::TemplateManager<TypeSeq,BaseT,ObjectT>::
end()
{
  return Sacado::TemplateIterator<TypeSeq,BaseT,ObjectT>(*this,
							 objects.end());
}
