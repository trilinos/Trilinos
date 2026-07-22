// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_MATERIALMODEL_ENTRY_HPP
#define PANZER_MATERIALMODEL_ENTRY_HPP

#include <iostream>
#include <string>
#include "Teuchos_ParameterList.hpp"

namespace panzer {

  //! \brief Class the holds parsed input data for material models.
  class MaterialModelEntry
  {

  public:

    MaterialModelEntry();

    MaterialModelEntry(const std::string factory_name);

    MaterialModelEntry(const std::string factory_name, 
		       const Teuchos::ParameterList& p);
    
    std::string factoryName() const;

    const Teuchos::ParameterList& params() const;

    void operator=(const MaterialModelEntry& e);
    
    bool operator==(const MaterialModelEntry& e) const;

    bool operator!=(const MaterialModelEntry& e) const;

    void print(std::ostream& os) const;

  protected:

    std::string m_factory_name;

    Teuchos::ParameterList m_params;

  };

  std::ostream& operator<<(std::ostream& os, panzer::MaterialModelEntry& m);

}

#endif
