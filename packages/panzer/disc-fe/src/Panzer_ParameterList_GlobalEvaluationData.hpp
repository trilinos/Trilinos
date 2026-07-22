// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_ParameterList_GlobalEvaluationData_hpp__
#define __Panzer_ParameterList_GlobalEvaluationData_hpp__

#include <string>
#include <vector>

#include "Panzer_GlobalEvaluationData.hpp"

namespace panzer {

/** This class is used by the residual scatter
  * to determine the name and indices of the active parameters
  * for scattering to the residual vector.
  */
class ParameterList_GlobalEvaluationData : public GlobalEvaluationData {
public:
   ParameterList_GlobalEvaluationData(const std::vector<std::string> & activeParameters) 
     : activeParameters_(activeParameters) {}
   virtual ~ParameterList_GlobalEvaluationData() {} 

   virtual void ghostToGlobal(int /* mem */) {}
   virtual void globalToGhost(int /* mem */) {}

   virtual bool requiresDirichletAdjustment() const { return false; }

   virtual void initializeData() {}

   const std::vector<std::string> & getActiveParameters() const
   { return activeParameters_; }

private:
   std::vector<std::string> activeParameters_;
};

}

#endif
