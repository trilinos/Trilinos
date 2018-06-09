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

//#define PANZER_USE_FAST_SUM 1
#define PANZER_USE_FAST_SUM 0

namespace panzer {

//**********************************************************************
template<typename EvalT, typename Traits>
Sum<EvalT, Traits>::
Sum(
  const Teuchos::ParameterList& p)
{
  std::string sum_name = p.get<std::string>("Sum Name");
  Teuchos::RCP<std::vector<std::string> > value_names = 
    p.get<Teuchos::RCP<std::vector<std::string> > >("Values Names");
  Teuchos::RCP<PHX::DataLayout> data_layout = 
    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout");

  TEUCHOS_ASSERT(static_cast<int>(value_names->size()) < MAX_VALUES);

  // check if the user wants to scale each term independently
  auto local_scalars = Kokkos::View<double *,PHX::Device>("scalars",value_names->size());
  if(p.isType<Teuchos::RCP<const std::vector<double> > >("Scalars")) {
    auto scalars_v = *p.get<Teuchos::RCP<const std::vector<double> > >("Scalars");

    // safety/sanity check
    TEUCHOS_ASSERT(scalars_v.size()==value_names->size());

    for (std::size_t i=0; i < value_names->size(); ++i)
      local_scalars(i) = scalars_v[i];
  }
  else {
    for (std::size_t i=0; i < value_names->size(); ++i)
      local_scalars(i) = 1.0;
  }

  scalars = local_scalars; 
  
  sum = PHX::MDField<ScalarT>(sum_name, data_layout);
  
  this->addEvaluatedField(sum);
 
  for (std::size_t i=0; i < value_names->size(); ++i) {
    values[i] = PHX::MDField<const ScalarT>( (*value_names)[i], data_layout);
    this->addDependentField(values[i]);
  }
  /*
  values.resize(value_names->size());
  for (std::size_t i=0; i < value_names->size(); ++i) {
    values[i] = PHX::MDField<const ScalarT>( (*value_names)[i], data_layout);
    this->addDependentField(values[i]);
  }
  */
 
  std::string n = "Sum Evaluator";
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
Sum<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData  /* worksets */,
  PHX::FieldManager<Traits>&  fm)
{
  this->utils.setFieldData(sum,fm);
  for (std::size_t i=0; i < scalars.extent(0); ++i)
    this->utils.setFieldData(values[i],fm);

  cell_data_size = sum.size() / sum.fieldTag().dataLayout().extent(0);
}


//**********************************************************************
template<typename EvalT, typename TRAITS>
template<unsigned int RANK>
KOKKOS_INLINE_FUNCTION
void Sum<EvalT, TRAITS>::operator() (PanzerSumTag<RANK>, const int &i) const{
  auto num_vals = scalars.extent(0);


  if (RANK == 1 )
  {
    for (std::size_t iv = 0; iv < num_vals; ++iv)
      sum(i) += scalars(iv)*(values[iv](i));
  }
  else if (RANK == 2)
  {
    const size_t dim_1 = sum.extent(1);
    for (std::size_t j = 0; j < dim_1; ++j)
      for (std::size_t iv = 0; iv < num_vals; ++iv)
        sum(i,j) += scalars(iv)*(values[iv](i,j));
  }
  else if (RANK == 3)
  {
    const size_t dim_1 = sum.extent(1),dim_2 = sum.extent(2);
    for (std::size_t j = 0; j < dim_1; ++j)
      for (std::size_t k = 0; k < dim_2; ++k)
        for (std::size_t iv = 0; iv < num_vals; ++iv)
          sum(i,j,k) += scalars(iv)*(values[iv](i,j,k));
  }
  else if (RANK == 4)
  {
    const size_t dim_1 = sum.extent(1),dim_2 = sum.extent(2),dim_3 = sum.extent(3);
    for (std::size_t j = 0; j < dim_1; ++j)
      for (std::size_t k = 0; k < dim_2; ++k)
        for (std::size_t l = 0; l < dim_3; ++l)
          for (std::size_t iv = 0; iv < num_vals; ++iv)
            sum(i,j,k,l) += scalars(iv)*(values[iv](i,j,k,l));
  }
  else if (RANK == 5)
  {
    const size_t dim_1 = sum.extent(1),dim_2 = sum.extent(2),dim_3 = sum.extent(3),dim_4 = sum.extent(4);
    for (std::size_t j = 0; j < dim_1; ++j)
      for (std::size_t k = 0; k < dim_2; ++k)
        for (std::size_t l = 0; l < dim_3; ++l)
          for (std::size_t m = 0; m < dim_4; ++m)
            for (std::size_t iv = 0; iv < num_vals; ++iv)
              sum(i,j,k,l,m) += scalars(iv)*(values[iv](i,j,k,l,m));
  }
  else if (RANK == 6)
  {
    const size_t dim_1 = sum.extent(1),dim_2 = sum.extent(2),dim_3 = sum.extent(3),dim_4 = sum.extent(4),dim_5 = sum.extent(5);
    for (std::size_t j = 0; j < dim_1; ++j)
      for (std::size_t k = 0; k < dim_2; ++k)
        for (std::size_t l = 0; l < dim_3; ++l)
          for (std::size_t m = 0; m < dim_4; ++m)
            for (std::size_t n = 0; n < dim_5; ++n)
              for (std::size_t iv = 0; iv < num_vals; ++iv)
                sum(i,j,k,l,m,n) += scalars(iv)*(values[iv](i,j,k,l,m,n));
  }
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
Sum<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData  /* workset */)
{   

  sum.deep_copy(ScalarT(0.0));

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
  size_t rank = sum.rank();
  const size_t length = sum.extent(0);
  if (rank == 1 )
  {
    Kokkos::parallel_for(Kokkos::RangePolicy<PanzerSumTag<1> >(0, length), *this);
  }
  else if (rank == 2)
  {
    Kokkos::parallel_for(Kokkos::RangePolicy<PanzerSumTag<2> >(0, length), *this);
  }
  else if (rank == 3)
  {
    Kokkos::parallel_for(Kokkos::RangePolicy<PanzerSumTag<3> >(0, length), *this);
  }
  else if (rank == 4)
  {
    Kokkos::parallel_for(Kokkos::RangePolicy<PanzerSumTag<4> >(0, length), *this);
  }
  else if (rank == 5)
  {
    Kokkos::parallel_for(Kokkos::RangePolicy<PanzerSumTag<5> >(0, length), *this);
  }
  else if (rank == 6)
  {
    Kokkos::parallel_for(Kokkos::RangePolicy<PanzerSumTag<6> >(0, length), *this);
  }
  else
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "ERROR: rank of sum is higher than supported");
  }

#endif
}


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
postRegistrationSetup(typename TRAITS::SetupData /* d */,
                      PHX::FieldManager<TRAITS>& fm)
{
  this->utils.setFieldData(sum,fm);
  for (std::size_t i=0; i < values.size(); ++i)
    this->utils.setFieldData(values[i],fm);
}

//**********************************************************************

template<typename EvalT, typename TRAITS,typename Tag0>
void SumStatic<EvalT,TRAITS,Tag0,void,void>::
evaluateFields(typename TRAITS::EvalData /* d */)
{
  sum.deep_copy(ScalarT(0.0));
  for (std::size_t i = 0; i < sum.extent(0); ++i)
    for (std::size_t d = 0; d < values.size(); ++d)
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

  // check if the user wants to scale each term independently
  if(p.isType<Teuchos::RCP<const std::vector<double> > >("Scalars")) {
    Teuchos::RCP<const std::vector<double> > scalar_values
        = p.get<Teuchos::RCP<const std::vector<double> > >("Scalars");

    // safety/sanity check
    TEUCHOS_ASSERT(scalar_values->size()==value_names->size());
    useScalars = true;

    Kokkos::View<double*,PHX::Device> scalars_nc
        = Kokkos::View<double*,PHX::Device>("scalars",scalar_values->size());

    for(std::size_t i=0;i<scalar_values->size();i++)
      scalars_nc(i) = (*scalar_values)[i];

    scalars = scalars_nc;
  }
  else {
    useScalars = false;
  }

  // sanity check
  TEUCHOS_ASSERT(data_layout->rank()==2);
  
  sum = PHX::MDField<ScalarT,Tag0,Tag1>(sum_name, data_layout);
  
  this->addEvaluatedField(sum);
 
  values.resize(value_names->size());
  for (std::size_t i=0; i < value_names->size(); ++i) {
    values[i] = PHX::MDField<const ScalarT,Tag0,Tag1>( (*value_names)[i], data_layout);
    this->addDependentField(values[i]);
  }
  numValues = value_names->size();
  TEUCHOS_ASSERT(numValues<=MAX_VALUES);
 
  std::string n = "SumStatic Rank 2 Evaluator";
  this->setName(n);
}

//**********************************************************************

template<typename EvalT, typename TRAITS,typename Tag0,typename Tag1>
SumStatic<EvalT,TRAITS,Tag0,Tag1,void>::
SumStatic(const std::vector<PHX::Tag<typename EvalT::ScalarT>> & inputs,
          const std::vector<double> & scalar_values,
          const PHX::Tag<typename EvalT::ScalarT> & output)
{
  TEUCHOS_ASSERT(scalar_values.size()==inputs.size());

  // check if the user wants to scale each term independently
  if(scalars.size()==0) {
    useScalars = false;
  }
  else {
    useScalars = true;

    Kokkos::View<double*,PHX::Device> scalars_nc
        = Kokkos::View<double*,PHX::Device>("scalars",scalar_values.size());

    for(std::size_t i=0;i<scalar_values.size();i++)
      scalars_nc(i) = scalar_values[i];

    scalars = scalars_nc;
  }

  // sanity check
  TEUCHOS_ASSERT(inputs.size()<=MAX_VALUES);
  
  sum = output;
  this->addEvaluatedField(sum);
 
  values.resize(inputs.size());
  for (std::size_t i=0; i < inputs.size(); ++i) {
    values[i] = inputs[i];
    this->addDependentField(values[i]);
  }

  numValues = inputs.size();

  std::string n = "SumStatic Rank 2 Evaluator";
  this->setName(n);
}

//**********************************************************************

template<typename EvalT, typename TRAITS,typename Tag0,typename Tag1>
void SumStatic<EvalT,TRAITS,Tag0,Tag1,void>::
postRegistrationSetup(typename TRAITS::SetupData /* d */,
                      PHX::FieldManager<TRAITS>& fm)
{
  this->utils.setFieldData(sum,fm);
  for (std::size_t i=0; i < values.size(); ++i) {
    this->utils.setFieldData(values[i],fm);
    value_views[i] = values[i].get_static_view();
  }
}

//**********************************************************************

template<typename EvalT, typename TRAITS,typename Tag0,typename Tag1>
void SumStatic<EvalT,TRAITS,Tag0,Tag1,void>::
evaluateFields(typename TRAITS::EvalData /* d */)
{
  sum.deep_copy(ScalarT(0.0));

  // Kokkos::parallel_for(sum.extent(0), *this);
  if(useScalars) 
    Kokkos::parallel_for(Kokkos::RangePolicy<PHX::Device,ScalarsTag>(0,sum.extent(0)), *this);
  else
    Kokkos::parallel_for(Kokkos::RangePolicy<PHX::Device,NoScalarsTag>(0,sum.extent(0)), *this);
}

//**********************************************************************

template<typename EvalT, typename TRAITS,typename Tag0,typename Tag1>
KOKKOS_INLINE_FUNCTION
void SumStatic<EvalT,TRAITS,Tag0,Tag1,void>::
operator()(const ScalarsTag, const unsigned c ) const
{
  for (int i=0;i<numValues;i++) {
    for (int j = 0; j < sum.extent_int(1); ++j)
      sum(c,j) += scalars(i)*value_views[i](c,j);
  }
}

//**********************************************************************

template<typename EvalT, typename TRAITS,typename Tag0,typename Tag1>
KOKKOS_INLINE_FUNCTION
void SumStatic<EvalT,TRAITS,Tag0,Tag1,void>::
operator()(const NoScalarsTag, const unsigned c ) const
{
  for (int i=0;i<numValues;i++) {
    for (int j = 0; j < sum.extent_int(1); ++j)
      sum(c,j) += value_views[i](c,j);
  }
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
    for (std::size_t i = 0; i < sum.extent(0); ++i)
      for (std::size_t j = 0; j < sum.extent(1); ++j)
        for (std::size_t k = 0; k < sum.extent(2); ++k)
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
