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

#ifndef PANZER_GLOBAL_STATISTICS_IMPL_HPP
#define PANZER_GLOBAL_STATISTICS_IMPL_HPP

#include "Intrepid_FunctionSpaceTools.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_String_Utilities.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_GlobalData.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_CommHelpers.hpp"
#include <boost/io/ios_state.hpp>
#include <iomanip>

namespace panzer {

//**********************************************************************
PHX_EVALUATOR_CTOR(GlobalStatistics,p)
{
  comm = p.get< Teuchos::RCP<const Teuchos::Comm<int> > >("Comm");

  global_data = p.get<Teuchos::RCP<panzer::GlobalData> >("Global Data");

  // Expects a string that is a Colon separated list of field names to compute statistics on.
  // for example the string "UX:UY:UZ:PRESSURE" would be separated into a vector with
  // four fields, "UX", "UY", "UZ", and "PRESSURE".
  std::string names_string = p.get<std::string>("Names");
  std::vector<std::string> names;
  panzer::StringTokenizer(names, names_string);

  Teuchos::RCP<panzer::IntegrationRule> ir = p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR");

  field_values.clear();
  for (typename std::vector<std::string>::const_iterator name = names.begin(); name != names.end(); ++name)
    field_values.push_back(PHX::MDField<ScalarT,Cell,IP>(*name, ir->dl_scalar));

  Teuchos::RCP<PHX::MDALayout<Cell> > cell_dl = Teuchos::rcp(new PHX::MDALayout<Cell>(ir->dl_scalar->dimension(0)));
  volumes = PHX::MDField<ScalarT,Cell>("Cell Volumes",cell_dl);

  tmp = PHX::MDField<ScalarT,Cell>("GlobalStatistics:tmp:"+names_string,cell_dl);
  ones = PHX::MDField<ScalarT,Cell,IP>("GlobalStatistics:ones:"+names_string,ir->dl_scalar);

  this->addEvaluatedField(volumes);
  this->addEvaluatedField(tmp);
  this->addEvaluatedField(ones);
  for (typename std::vector<PHX::MDField<ScalarT,Cell,IP> >::const_iterator field = field_values.begin();
       field != field_values.end(); ++field) {
    this->addDependentField(*field);
  }

  averages.resize(field_values.size());
  maxs.resize(field_values.size());
  mins.resize(field_values.size());
  global_maxs.resize(field_values.size());
  global_mins.resize(field_values.size());
  global_averages.resize(field_values.size());

  ir_order = ir->cubature_degree;

  std::string n = "GlobalStatistics: " + names_string;
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(GlobalStatistics,sd,fm)
{
  this->utils.setFieldData(volumes,fm);
  this->utils.setFieldData(tmp,fm);
  this->utils.setFieldData(ones,fm);
  
  for (typename std::vector<PHX::MDField<ScalarT,Cell,IP> >::iterator field = field_values.begin();
       field != field_values.end(); ++field)
    this->utils.setFieldData(*field,fm);

  ir_index = panzer::getIntegrationRuleIndex(ir_order,(*sd.worksets_)[0]);

  for (typename PHX::MDField<ScalarT,Cell,IP>::size_type cell = 0; cell < ones.dimension(0); ++cell)
    for (typename PHX::MDField<ScalarT,Cell,IP>::size_type ip = 0; ip < ones.dimension(1); ++ip)
      ones(cell,ip) = 1.0;
}

//**********************************************************************
PHX_EVALUATE_FIELDS(GlobalStatistics,workset)
{ 
  if (workset.num_cells == 0)
    return;

  Intrepid::FunctionSpaceTools::integrate<ScalarT>(volumes, ones, 
						   (workset.int_rules[ir_index])->weighted_measure, 
						   Intrepid::COMP_BLAS);

  for (std::size_t cell = 0; cell < workset.num_cells; ++cell)
    total_volume += volumes(cell);

  typename std::vector<PHX::MDField<ScalarT,Cell,IP> >::size_type field_index = 0;
  for (typename std::vector<PHX::MDField<ScalarT,Cell,IP> >::iterator field = field_values.begin();
       field != field_values.end(); ++field,++field_index) {
    
    Intrepid::FunctionSpaceTools::integrate<ScalarT>(tmp, *field, 
						     (workset.int_rules[ir_index])->weighted_measure, 
						     Intrepid::COMP_BLAS);
    
    for (std::size_t cell = 0; cell < workset.num_cells; ++cell) {
      averages[field_index] += tmp(cell);

      for (typename PHX::MDField<ScalarT,Cell,IP>::size_type ip = 0; ip < (field->dimension(1)); ++ip) {
	maxs[field_index] = std::max( (*field)(cell,ip), maxs[field_index]);
	mins[field_index] = std::min( (*field)(cell,ip), mins[field_index]);
      }
    }
    
  }

}

//**********************************************************************
PHX_PRE_EVALUATE_FIELDS(GlobalStatistics,data)
{
  total_volume = Teuchos::ScalarTraits<ScalarT>::zero();

  for (typename std::vector<ScalarT>::iterator field = averages.begin(); field != averages.end(); ++field)
    *field = Teuchos::ScalarTraits<ScalarT>::zero();

  for (typename std::vector<ScalarT>::iterator field = maxs.begin(); field != maxs.end(); ++field)
    *field = Teuchos::ScalarTraits<ScalarT>::rmin();

  for (typename std::vector<ScalarT>::iterator field = mins.begin(); field != mins.end(); ++field)
    *field = Teuchos::ScalarTraits<ScalarT>::rmax();
}

//**********************************************************************
PHX_POST_EVALUATE_FIELDS(GlobalStatistics,data)
{
  this->postprocess(*(global_data->os));
}

//**********************************************************************
template<typename EvalT, typename Traits>
void GlobalStatistics<EvalT, Traits>::postprocess(std::ostream& os)
{
  // throw unless specialized for residual evaluations
  TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"SHOULD NEVER BE CALLED!");
}

//**********************************************************************
template<>
void GlobalStatistics<panzer::Traits::Residual, panzer::Traits>::postprocess(std::ostream& os)
{
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, static_cast<int>(1), &total_volume, &global_total_volume);
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, static_cast<int>(averages.size()), &averages[0], &global_averages[0]);
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, static_cast<int>(maxs.size()), &maxs[0], &global_maxs[0]);
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, static_cast<int>(mins.size()), &mins[0], &global_mins[0]);

  for (std::vector<ScalarT>::size_type i = 0; i < field_values.size(); ++i)
    global_averages[i] /= global_total_volume;

  if (comm->getRank() == 0) {
    
    boost::io::ios_all_saver saver(os);
    
    std::size_t precision = 8;
    os << std::scientific << std::showpoint << std::setprecision(precision) << std::left;
    
    std::size_t name_width = 0;
    for (std::vector<ScalarT>::size_type i = 0; i < field_values.size(); ++i)
      name_width = std::max(name_width,field_values[i].fieldTag().name().size());
    
    std::size_t value_width = precision + 7;
    
    os << std::setw(name_width) << "Field" 
       << " " << std::setw(value_width) << "Average" 
       << " " << std::setw(value_width) << "Maximum (@IP)" 
       << " " << std::setw(value_width) << "Minimum (@IP)" 
       << std::endl;
    
    for (std::vector<ScalarT>::size_type i = 0; i < field_values.size(); ++i) {
      os << std::setw(name_width) <<  field_values[i].fieldTag().name() 
	 << " " << std::setw(value_width) << global_averages[i]
	 << " " << std::setw(value_width) << global_maxs[i]
	 << " " << std::setw(value_width) << global_mins[i] << std::endl;
    }

  }

}

//**********************************************************************
template<typename EvalT, typename Traits>
const PHX::FieldTag& GlobalStatistics<EvalT, Traits>::getRequiredFieldTag()
{
  return tmp.fieldTag();
}


} // namespace panzer



#endif

