#ifndef PANZER_GLOBAL_STATISTICS_HPP
#define PANZER_GLOBAL_STATISTICS_HPP

#include <iostream>
#include <string>
#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Teuchos_Comm.hpp"

namespace panzer {
    
PHX_EVALUATOR_CLASS_PP(GlobalStatistics)
  
  PHX::MDField<ScalarT,Cell> volumes;
    
  PHX::MDField<ScalarT,Cell> tmp;

  PHX::MDField<ScalarT,Cell,IP> ones;

  std::vector<PHX::MDField<ScalarT,Cell,IP> > field_values;

  ScalarT total_volume;
  std::vector<ScalarT> averages;
  std::vector<ScalarT> maxs;
  std::vector<ScalarT> mins;
  ScalarT global_total_volume;
  std::vector<ScalarT> global_averages;
  std::vector<ScalarT> global_maxs;
  std::vector<ScalarT> global_mins;

  int ir_order;
  std::size_t ir_index;

  Teuchos::RCP<const Teuchos::Comm<int> > comm;

  void postprocess(std::ostream& os);

public:
  const PHX::FieldTag& getRequiredFieldTag();

PHX_EVALUATOR_CLASS_END

}

#include "Panzer_GlobalStatisticsT.hpp"

#endif
