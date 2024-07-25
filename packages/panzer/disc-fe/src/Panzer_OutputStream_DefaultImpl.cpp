// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_OutputStream_DefaultImpl.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_Assert.hpp"

namespace panzer {
  
  OutputStreamDefaultImpl::OutputStreamDefaultImpl()
  { }
  
  OutputStreamDefaultImpl::~OutputStreamDefaultImpl()  
  { }  

  void OutputStreamDefaultImpl::
  setOStream(const Teuchos::RCP<Teuchos::FancyOStream>& os)
  {
    m_out = os;
    m_out->setOutputToRootOnly(0); 
    m_pout = Teuchos::rcp(new Teuchos::FancyOStream(os->getOStream()));
    m_pout->copyAllOutputOptions(*m_out);
    m_pout->setOutputToRootOnly(-1);
  }

  Teuchos::RCP<Teuchos::FancyOStream> 
  OutputStreamDefaultImpl::getOStream() const
  {
    return m_out;
  }

  Teuchos::FancyOStream& OutputStreamDefaultImpl::out() const
  {
    return *m_out;
  }

  Teuchos::FancyOStream& OutputStreamDefaultImpl::pout() const
  {
    return *m_pout;
  }

  void OutputStreamDefaultImpl::setVerbosityLevel(EVerbosityLevel vl)
  {
    m_level = vl;
  }

  EVerbosityLevel OutputStreamDefaultImpl::getVerbosityLevel() const
  {
    return m_level;
  }
  
  bool OutputStreamDefaultImpl::
  doOutput(EVerbosityLevel vl, bool only_for_exact_level) const
  {
    if ( !only_for_exact_level && 
	 (Teuchos::as<int>(vl) >= Teuchos::as<int>(m_level)) )
      return true;
    
    if (only_for_exact_level && (vl == m_level) )
      return true;
    
    return false;
  }

}
