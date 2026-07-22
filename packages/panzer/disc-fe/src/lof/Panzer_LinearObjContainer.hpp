// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_LinearObjContainer_hpp__
#define __Panzer_LinearObjContainer_hpp__

#include "PanzerDiscFE_config.hpp"

#include "Panzer_GlobalEvaluationData.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace panzer {

// **********************************************************************************************
// *********************** LINEAR OBJ CONTAINER *************************************************
// **********************************************************************************************

class LinearObjContainer : public GlobalEvaluationData_Default {
public:
   virtual ~LinearObjContainer() {}

   typedef enum { X=0x1, DxDt=0x2, F=0x4, Mat=0x8} Members;

   virtual void initialize() = 0;
};

}

#endif
