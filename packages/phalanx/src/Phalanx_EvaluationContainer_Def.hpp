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


#ifndef PHX_SCALAR_CONTAINER_DEF_HPP
#define PHX_SCALAR_CONTAINER_DEF_HPP

#include "Teuchos_Assert.hpp"
#include "Phalanx_Evaluator.hpp"
#include "Phalanx_TypeStrings.hpp"
#include "Phalanx_KokkosViewFactoryFunctor.hpp"
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/at.hpp>
#include <sstream>

// *************************************************************************
template <typename EvalT, typename Traits>
PHX::EvaluationContainer<EvalT, Traits>::EvaluationContainer() :
  post_registration_setup_called_(false)
{
  this->vp_manager_.setEvaluationTypeName( PHX::typeAsString<EvalT>() );
}

// *************************************************************************
template <typename EvalT, typename Traits> 
PHX::EvaluationContainer<EvalT, Traits>::~EvaluationContainer()
{

}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
requireField(const PHX::FieldTag& f)
{
  this->vp_manager_.requireField(f);
}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
registerEvaluator(const Teuchos::RCP<PHX::Evaluator<Traits> >& p)
{
  this->vp_manager_.registerEvaluator(p);
}

// *************************************************************************
template <typename EvalT, typename Traits> 
void PHX::EvaluationContainer<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData d,
		      PHX::FieldManager<Traits>& fm)
{
  // Figure out all evaluator dependencies
  if ( !(this->vp_manager_.sortingCalled()) )
    this->vp_manager_.sortAndOrderEvaluators();
  
  const std::vector< Teuchos::RCP<PHX::FieldTag> >& var_list = 
    this->vp_manager_.getFieldTags();

  std::vector< Teuchos::RCP<PHX::FieldTag> >::const_iterator  var;

  for (var = var_list.begin(); var != var_list.end(); ++var) {
    typedef typename boost::mpl::at<typename Traits::EvalToDataMap,EvalT>::type EvalDataTypes;
    Sacado::mpl::for_each<EvalDataTypes>(PHX::KokkosViewFactoryFunctor<EvalT>(fields_,*(*var),kokkos_extended_data_type_dimensions_));

    TEUCHOS_TEST_FOR_EXCEPTION(fields_.find((*var)->identifier()) == fields_.end(),std::runtime_error,
			       "Error: PHX::EvaluationContainer::postRegistrationSetup(): could not build a Kokkos::View for field named \"" << (*var)->name() << "\" of type \"" << (*var)->dataTypeInfo().name() << "\" for the evaluation type \"" << PHX::typeAsString<EvalT>() << "\".");
  }

  // Allow fields in evaluators to grab pointers to relevant field data
  this->vp_manager_.postRegistrationSetup(d,fm);

  post_registration_setup_called_ = true;
}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
evaluateFields(typename Traits::EvalData d)
{
#ifdef PHX_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( !(this->setupCalled()) , std::logic_error,
		      "You must call post registration setup for each evaluation type before calling the evaluateFields() method for that type!");
#endif

  this->vp_manager_.evaluateFields(d);
}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
preEvaluate(typename Traits::PreEvalData d)
{
#ifdef PHX_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( !(this->setupCalled()) , std::logic_error,
		      "You must call post registration setup for each evaluation type before calling the preEvaluate() method for that type!");
#endif

  this->vp_manager_.preEvaluate(d);
}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
postEvaluate(typename Traits::PostEvalData d)
{
#ifdef PHX_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( !(this->setupCalled()) , std::logic_error,
		      "You must call post registration setup for each evaluation type before calling the postEvaluate() method for that type!");
#endif

  this->vp_manager_.postEvaluate(d);
}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
setKokkosExtendedDataTypeDimensions(const std::vector<PHX::index_size_type>& dims)
{
  kokkos_extended_data_type_dimensions_ = dims;
}

// *************************************************************************
template <typename EvalT, typename Traits>
const std::vector<PHX::index_size_type> & PHX::EvaluationContainer<EvalT, Traits>::
getKokkosExtendedDataTypeDimensions() const 
{
#ifdef PHX_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( !(this->setupCalled()) , std::logic_error,
		      "You must call post registration setup for each evaluation type before calling the postEvaluate() method for that type!");
#endif
  return kokkos_extended_data_type_dimensions_;
}

// *************************************************************************
template <typename EvalT, typename Traits>
boost::any
PHX::EvaluationContainer<EvalT, Traits>::getFieldData(const PHX::FieldTag& f)
{
  //return fields_[f.identifier()];
  boost::unordered_map<std::string,boost::any>::iterator a= fields_.find(f.identifier());
   if (a==fields_.end()){
    std::cout << " PHX::EvaluationContainer<EvalT, Traits>::getFieldData can't find an f.identifier() "<<  f.identifier() << std::endl;
   }
  return a->second;
}


// *************************************************************************
template <typename EvalT, typename Traits>
bool PHX::EvaluationContainer<EvalT, Traits>::setupCalled() const
{
  return post_registration_setup_called_;
}

// *************************************************************************
template <typename EvalT, typename Traits>
const std::string PHX::EvaluationContainer<EvalT, Traits>::
evaluationType() const
{
  return PHX::typeAsString<EvalT>();
}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::print(std::ostream& os) const
{
  std::string type = PHX::typeAsString<EvalT>();

  os << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  os << "Starting PHX::EvaluationContainer Output" << std::endl;
  os << "Evaluation Type = " << type << std::endl;
  os << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  os << this->vp_manager_ << std::endl;
  for (boost::unordered_map<std::string,boost::any>::const_iterator i = 
	 fields_.begin(); i != fields_.end(); ++i)
    os << i->first << std::endl;
  os << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  os << "Finished PHX::EvaluationContainer Output" << std::endl;
  os << "Evaluation Type = " << type << std::endl;
  os << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  os << std::endl;
}

// *************************************************************************

#endif 
