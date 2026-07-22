// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_PARAMETER_LIBRARY_ACCEPTOR_DEFAULT_IMPL_HPP
#define PANZER_PARAMETER_LIBRARY_ACCEPTOR_DEFAULT_IMPL_HPP

#include "Panzer_ParameterLibraryAcceptor.hpp"
#include "Teuchos_Assert.hpp"

namespace panzer {

  /** \brief Pure Virtual base class for accepting the parameter library

      This class is used to retrieve the parameter library from an
      object.
  */
  class ParameterLibraryAcceptor_DefaultImpl :
    public ParameterLibraryAcceptor {

  public:
    
    ParameterLibraryAcceptor_DefaultImpl() {}
    
    ParameterLibraryAcceptor_DefaultImpl(const Teuchos::RCP<panzer::ParamLib>& pl) : 
      m_param_lib(pl) {}
    
    virtual ~ParameterLibraryAcceptor_DefaultImpl() {}

    virtual void setParameterLibrary(const Teuchos::RCP<panzer::ParamLib>& pl)
    {
      TEUCHOS_TEST_FOR_EXCEPTION(nonnull(m_param_lib), std::logic_error,
				 "A parameter library has already been set on this object!");
      m_param_lib = pl;
    }

    virtual Teuchos::RCP<panzer::ParamLib> getParameterLibrary() const
    { return m_param_lib; }

  private:

    Teuchos::RCP<panzer::ParamLib> m_param_lib;

  };

}

#endif
