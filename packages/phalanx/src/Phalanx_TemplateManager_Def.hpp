// $Id$ 
// $Source$ 
// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
