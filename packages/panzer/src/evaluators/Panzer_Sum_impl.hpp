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

#ifndef PANZER_SUM_IMPL_HPP
#define PANZER_SUM_IMPL_HPP

#include <cstddef>
#include <string>
#include <vector>

#include "Phalanx_MDField_Utilities.hpp"

#define PANZER_USE_FAST_SUM 1
// #define PANZER_USE_FAST_SUM 0

namespace panzer {

//**********************************************************************
PHX_EVALUATOR_CTOR(Sum,p)
{
  std::string sum_name = p.get<std::string>("Sum Name");
  Teuchos::RCP<std::vector<std::string> > value_names = 
    p.get<Teuchos::RCP<std::vector<std::string> > >("Values Names");
  Teuchos::RCP<PHX::DataLayout> data_layout = 
    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout");

  // check if the user wants to scale each term independently
  if(p.isType<Teuchos::RCP<const std::vector<double> > >("Scalars")) {
    scalars = *p.get<Teuchos::RCP<const std::vector<double> > >("Scalars");

    // safety/sanity check
    TEUCHOS_ASSERT(scalars.size()==value_names->size());
  }
  else {
    // otherwise use all ones (a simple sum)
    scalars = std::vector<double>(value_names->size(),1.0);
  }
  
  sum = PHX::MDField<ScalarT>(sum_name, data_layout);
  
  this->addEvaluatedField(sum);
 
  values.resize(value_names->size());
  for (std::size_t i=0; i < value_names->size(); ++i) {
    values[i] = PHX::MDField<const ScalarT>( (*value_names)[i], data_layout);
    this->addDependentField(values[i]);
  }
 
  std::string n = "Sum Evaluator";
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(Sum,worksets,fm)
{
  this->utils.setFieldData(sum,fm);
  for (std::size_t i=0; i < values.size(); ++i)
    this->utils.setFieldData(values[i],fm);

  cell_data_size = sum.size() / sum.fieldTag().dataLayout().dimension(0);
}

//**********************************************************************
PHX_EVALUATE_FIELDS(Sum,workset)
{ 
#if PANZER_USE_FAST_SUM 
  sum.deep_copy(ScalarT(0.0));
  for (std::size_t j = 0; j < values.size(); ++j) {
    
    PHX::MDFieldIterator<ScalarT> sum_it(sum);
    PHX::MDFieldIterator<const ScalarT> values_it(values[j]);
    // for (PHX::MDFieldIterator<ScalarT> sum_it(sum), values_it(values[j]);
    for ( ;
         ! (sum_it.done() || values_it.done());
         ++sum_it, ++values_it)
      *sum_it += scalars[j]*(*values_it);
  }
#else
  std::size_t length = workset.num_cells * cell_data_size;
  for (std::size_t i = 0; i < length; ++i) {
    sum[i] = 0.0;
    for (std::size_t j = 0; j < values.size(); ++j)
      sum[i] += scalars[j]*(values[j][i]);
  }
#endif
}

//**********************************************************************
//**********************************************************************

template<typename EvalT, typename TRAITS,typename Tag0>
SumStatic<EvalT,TRAITS,Tag0,void,void>::
SumStatic(const Teuchos::ParameterList& p)
{
  std::string sum_name = p.get<std::string>("Sum Name");
  Teuchos::RCP<std::vector<std::string> > value_names = 
    p.get<Teuchos::RCP<std::vector<std::string> > >("Values Names");
  Teuchos::RCP<PHX::DataLayout> data_layout = 
    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout");
  
  // sanity check
  TEUCHOS_ASSERT(data_layout->rank()==1);
  
  sum = PHX::MDField<ScalarT,Tag0>(sum_name, data_layout);
  
  this->addEvaluatedField(sum);
 
  values.resize(value_names->size());
  for (std::size_t i=0; i < value_names->size(); ++i) {
    values[i] = PHX::MDField<const ScalarT,Tag0>( (*value_names)[i], data_layout);
    this->addDependentField(values[i]);
  }
 
  std::string n = "SumStatic Rank 1 Evaluator";
  this->setName(n);
}

//**********************************************************************

template<typename EvalT, typename TRAITS,typename Tag0>
void SumStatic<EvalT,TRAITS,Tag0,void,void>::
postRegistrationSetup(typename TRAITS::SetupData d,
                      PHX::FieldManager<TRAITS>& fm)
{
  this->utils.setFieldData(sum,fm);
  for (std::size_t i=0; i < values.size(); ++i)
    this->utils.setFieldData(values[i],fm);
}

//**********************************************************************

template<typename EvalT, typename TRAITS,typename Tag0>
void SumStatic<EvalT,TRAITS,Tag0,void,void>::
evaluateFields(typename TRAITS::EvalData d)
{
  sum.deep_copy(ScalarT(0.0));
  
  for (std::size_t d = 0; d < values.size(); ++d)
    for (std::size_t i = 0; i < sum.dimension_0(); ++i)
      sum(i) += (values[d])(i);
}

//**********************************************************************
//**********************************************************************

template<typename EvalT, typename TRAITS,typename Tag0,typename Tag1>
SumStatic<EvalT,TRAITS,Tag0,Tag1,void>::
SumStatic(const Teuchos::ParameterList& p)
{
  std::string sum_name = p.get<std::string>("Sum Name");
  Teuchos::RCP<std::vector<std::string> > value_names = 
    p.get<Teuchos::RCP<std::vector<std::string> > >("Values Names");
  Teuchos::RCP<PHX::DataLayout> data_layout = 
    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout");
  
  // sanity check
  TEUCHOS_ASSERT(data_layout->rank()==2);
  
  sum = PHX::MDField<ScalarT,Tag0,Tag1>(sum_name, data_layout);
  
  this->addEvaluatedField(sum);
 
  values.resize(value_names->size());
  for (std::size_t i=0; i < value_names->size(); ++i) {
    values[i] = PHX::MDField<const ScalarT,Tag0,Tag1>( (*value_names)[i], data_layout);
    this->addDependentField(values[i]);
  }
 
  std::string n = "SumStatic Rank 2 Evaluator";
  this->setName(n);
}

//**********************************************************************

template<typename EvalT, typename TRAITS,typename Tag0,typename Tag1>
void SumStatic<EvalT,TRAITS,Tag0,Tag1,void>::
postRegistrationSetup(typename TRAITS::SetupData d,
                      PHX::FieldManager<TRAITS>& fm)
{
  this->utils.setFieldData(sum,fm);
  for (std::size_t i=0; i < values.size(); ++i)
    this->utils.setFieldData(values[i],fm);
}

//**********************************************************************

template<typename EvalT, typename TRAITS,typename Tag0,typename Tag1>
void SumStatic<EvalT,TRAITS,Tag0,Tag1,void>::
evaluateFields(typename TRAITS::EvalData d)
{
  sum.deep_copy(ScalarT(0.0));
  
  for (std::size_t d = 0; d < values.size(); ++d) {
    current_value = values[d];
    Kokkos::parallel_for(sum.dimension_0(), *this);
  }

/*
  for (std::size_t d = 0; d < values.size(); ++d)
    for (std::size_t i = 0; i < sum.dimension_0(); ++i)
      for (std::size_t j = 0; j < sum.dimension_1(); ++j)
        sum(i,j) += (values[d])(i,j);
*/
}

//**********************************************************************

template<typename EvalT, typename TRAITS,typename Tag0,typename Tag1>
KOKKOS_INLINE_FUNCTION
void SumStatic<EvalT,TRAITS,Tag0,Tag1,void>::
operator()( const unsigned c ) const
{
  for (std::size_t j = 0; j < sum.dimension_1(); ++j)
    sum(c,j) += current_value(c,j);
}


//**********************************************************************
//**********************************************************************

/*
template<typename EvalT, typename TRAITS,typename Tag0,typename Tag1,typename Tag2>
SumStatic<EvalT,TRAITS,Tag0,Tag1,Tag2>::
SumStatic(const Teuchos::ParameterList& p)
{
  std::string sum_name = p.get<std::string>("Sum Name");
  Teuchos::RCP<std::vector<std::string> > value_names = 
    p.get<Teuchos::RCP<std::vector<std::string> > >("Values Names");
  Teuchos::RCP<PHX::DataLayout> data_layout = 
    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout");
  
  // sanity check
  TEUCHOS_ASSERT(data_layout->rank()==3);
  
  sum = PHX::MDField<ScalarT,Tag0,Tag1,Tag2>(sum_name, data_layout);
  
  this->addEvaluatedField(sum);
 
  values.resize(value_names->size());
  for (std::size_t i=0; i < value_names->size(); ++i) {
    values[i] = PHX::MDField<ScalarT,Tag0,Tag1,Tag2>( (*value_names)[i], data_layout);
    this->addDependentField(values[i]);
  }
 
  std::string n = "Sum Evaluator";
  this->setName(n);
}
*/

//**********************************************************************
/*

template<typename EvalT, typename TRAITS,typename Tag0,typename Tag1,typename Tag2>
void SumStatic<EvalT,TRAITS,Tag0,Tag1,Tag2>::
postRegistrationSetup(typename TRAITS::SetupData d,
                      PHX::FieldManager<TRAITS>& fm)
{
  this->utils.setFieldData(sum,fm);
  for (std::size_t i=0; i < values.size(); ++i)
    this->utils.setFieldData(values[i],fm);
}
*/

//**********************************************************************

/*
template<typename EvalT, typename TRAITS,typename Tag0,typename Tag1,typename Tag2>
void SumStatic<EvalT,TRAITS,Tag0,Tag1,Tag2>::
evaluateFields(typename TRAITS::EvalData d)
{
  sum.deep_copy(ScalarT(0.0));
  
  for (std::size_t d = 0; d < values.size(); ++d)
    for (std::size_t i = 0; i < sum.dimension_0(); ++i)
      for (std::size_t j = 0; j < sum.dimension_1(); ++j)
        for (std::size_t k = 0; k < sum.dimension_2(); ++k)
          sum(i,j,k) += (values[d])(i);
}
*/

//**********************************************************************
//**********************************************************************

template<typename EvalT, typename TRAITS,typename Tag0,typename Tag1,typename Tag2>
Teuchos::RCP<PHX::Evaluator<TRAITS> > 
buildStaticSumEvaluator(const std::string & sum_name,
                        const std::vector<std::string> & value_names,
                        const Teuchos::RCP<PHX::DataLayout> & data_layout)
{
  Teuchos::ParameterList p;
  p.set<std::string>("Sum Name",sum_name);
  p.set<Teuchos::RCP<std::vector<std::string> > >("Values Names",Teuchos::rcp(new std::vector<std::string>(value_names)));
  p.set< Teuchos::RCP<PHX::DataLayout> >("Data Layout",data_layout);
   
  return Teuchos::rcp(new SumStatic<EvalT,TRAITS,Tag0,Tag1,Tag2>(p));
}

}

#endif
