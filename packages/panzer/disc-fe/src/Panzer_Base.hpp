// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_BASE_H
#define PANZER_BASE_H

namespace panzer {

  //! Non-templated empty base class for template managers
  class Base {
    
  public:
    
    Base() {}
    
    virtual ~Base() {}
    
  };
  
}

#endif
