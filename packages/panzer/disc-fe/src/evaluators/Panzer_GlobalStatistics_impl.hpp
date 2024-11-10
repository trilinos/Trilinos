// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_GLOBAL_STATISTICS_IMPL_HPP
#define PANZER_GLOBAL_STATISTICS_IMPL_HPP

#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_String_Utilities.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_IosAllSaver.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_CommHelpers.hpp"
#include <iomanip>

namespace panzer {

//**********************************************************************
template<typename EvalT, typename Traits>
GlobalStatistics<EvalT, Traits>::
GlobalStatistics(
  const Teuchos::ParameterList& p)
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
    field_values.push_back(PHX::MDField<const ScalarT,Cell,IP>(*name, ir->dl_scalar));

  Teuchos::RCP<PHX::MDALayout<Cell> > cell_dl = Teuchos::rcp(new PHX::MDALayout<Cell>(ir->dl_scalar->extent(0)));
  volumes = PHX::MDField<ScalarT,Cell>("Cell Volumes",cell_dl);

  tmp = PHX::MDField<ScalarT,Cell>("GlobalStatistics:tmp:"+names_string,cell_dl);
  ones = PHX::MDField<ScalarT,Cell,IP>("GlobalStatistics:ones:"+names_string,ir->dl_scalar);

  this->addEvaluatedField(volumes);
  this->addEvaluatedField(tmp);
  this->addEvaluatedField(ones);
  for (typename std::vector<PHX::MDField<const ScalarT,Cell,IP> >::const_iterator field = field_values.begin();
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
template<typename EvalT, typename Traits>
void
GlobalStatistics<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  ir_index = panzer::getIntegrationRuleIndex(ir_order,(*sd.worksets_)[0], this->wda);
  auto l_ones = ones.get_static_view();
  Kokkos::parallel_for("GlobalStatistics", l_ones.extent(0), KOKKOS_LAMBDA(int cell) {
      for (std::size_t ip = 0; ip < l_ones.extent(1); ++ip)
	l_ones(cell,ip) = 1.0;
    });
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
GlobalStatistics<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{ 
  if (workset.num_cells == 0)
    return;

  Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::integrate(volumes.get_view(),
                                                                         ones.get_view(), 
                                                                         (this->wda(workset).int_rules[ir_index])->weighted_measure.get_view());
  auto volumes_h = Kokkos::create_mirror_view(as_view(volumes));
  Kokkos::deep_copy(volumes_h, as_view(volumes));

  for (index_t cell = 0; cell < workset.num_cells; ++cell)
    total_volume += volumes_h(cell);

  typename std::vector<PHX::MDField<ScalarT,Cell,IP> >::size_type field_index = 0;
  for (typename std::vector<PHX::MDField<const ScalarT,Cell,IP> >::iterator field = field_values.begin();
       field != field_values.end(); ++field,++field_index) {
    
    Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::integrate(tmp.get_view(),
                                                                           field->get_view(), 
                                                                           (this->wda(workset).int_rules[ir_index])->weighted_measure.get_view());
    auto tmp_h = Kokkos::create_mirror_view(tmp.get_static_view());
    auto field_h = Kokkos::create_mirror_view( field->get_static_view());
    Kokkos::deep_copy(tmp_h, tmp.get_static_view());
    Kokkos::deep_copy(field_h, field->get_static_view());


    for (index_t cell = 0; cell < workset.num_cells; ++cell) {
      averages[field_index] += tmp_h(cell);

      for (typename PHX::MDField<ScalarT,Cell,IP>::size_type ip = 0; ip < (field->extent(1)); ++ip) {
        maxs[field_index] = std::max( field_h(cell,ip), maxs[field_index]);
        mins[field_index] = std::min( field_h(cell,ip), mins[field_index]);
      }
    }
    
  }
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
GlobalStatistics<EvalT, Traits>::
preEvaluate(
  typename Traits::PreEvalData  /* data */)
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
template<typename EvalT, typename Traits>
void
GlobalStatistics<EvalT, Traits>::
postEvaluate(
  typename Traits::PostEvalData  /* data */)
{
  this->postprocess(*(global_data->os));
}

//**********************************************************************
template<typename EvalT, typename TRAITS>
void GlobalStatistics<EvalT, TRAITS>::postprocess(std::ostream& /* os */)
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

    panzer::ios_all_saver saver(os);
    
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
template<typename EvalT, typename TRAITS>
const PHX::FieldTag& GlobalStatistics<EvalT, TRAITS>::getRequiredFieldTag()
{
  return tmp.fieldTag();
}


} // namespace panzer



#endif

