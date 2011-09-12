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


#ifndef PHX_FIELD_EVALUATOR_FACTORY_DEF_HPP
#define PHX_FIELD_EVALUATOR_FACTORY_DEF_HPP

#include <sstream>
#include "boost/mpl/size.hpp"
#include "boost/mpl/for_each.hpp"
#include "boost/mpl/range_c.hpp"
#include "Teuchos_TestForException.hpp"
#include "Phalanx_Evaluator_Factory_UFO.hpp"

//**********************************************************************
template<typename Traits, typename FactoryTraits>
Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator_TemplateManager<Traits> > > > 
PHX::EvaluatorFactory<Traits, FactoryTraits>::
buildEvaluators(const std::map<std::string, Teuchos::RCP<Teuchos::ParameterList> >& data)
{
  Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator_TemplateManager<Traits> > > > vector_tm = 
    Teuchos::rcp(new std::vector< Teuchos::RCP<Evaluator_TemplateManager<Traits> > >);

  std::map<std::string, Teuchos::RCP<Teuchos::ParameterList> >::const_iterator
    it = data.begin();
  
  for (; it != data.end(); ++it) {
    
    Teuchos::RCP< PHX::Evaluator_TemplateManager<Traits> > tm = 
      Teuchos::rcp(new PHX::Evaluator_TemplateManager<Traits>);
    
    Teuchos::RCP<Teuchos::ParameterList> p = it->second;

    bool found_type = false;
    int object_type = p->get<int>("Type");
    static const int size = boost::mpl::size<typename FactoryTraits::EvaluatorTypes>::value;
    boost::mpl::for_each< boost::mpl::range_c<int,0,size> >( UFO<Traits,FactoryTraits>(object_type, p, tm, found_type) );

    if (!found_type) {
      std::ostringstream msg;
      msg << "Unable to find model in EvaluatorFactory for "
	  << it->first
	  << ".  Please make sure you have set a valid integer "
	  << "for \"Type\" in the parameter list!";
      TEST_FOR_EXCEPTION(!found_type, std::logic_error, msg.str());
    }

    vector_tm->push_back(tm);

  }
  
  return vector_tm;

}
//**********************************************************************
template<typename Traits>
void PHX::registerEvaluators(const Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator_TemplateManager<Traits> > > >& providers, PHX::FieldManager<Traits>& fm)
{
  // Loop over each provider template manager
  typename std::vector< Teuchos::RCP<Evaluator_TemplateManager<Traits> > >::iterator 
    tm = providers->begin();
  for (; tm != providers->end(); ++tm) {
    
    // Loop over Evaluation Types
    typename PHX::FieldManager<Traits>::iterator vmit = fm.begin();
    typename Evaluator_TemplateManager<Traits>::iterator vpit = 
      (*tm)->begin();
    for (; vpit != (*tm)->end(); ++vpit) {
      Teuchos::RCP<PHX::Evaluator<Traits> > vp =
	Teuchos::rcp_dynamic_cast<PHX::Evaluator<Traits> >(vpit.rcp());
      fm.registerEvaluator(vmit, vp);
      ++vmit;
    } 
  }

}
//**********************************************************************

#endif
