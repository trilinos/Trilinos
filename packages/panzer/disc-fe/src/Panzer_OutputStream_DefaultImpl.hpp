// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_OUTPUT_STREAM_DEFAULT_IMPL_HPP
#define PANZER_OUTPUT_STREAM_DEFAULT_IMPL_HPP

#include "Teuchos_RCP.hpp"
#include "Panzer_OutputStream.hpp"

namespace panzer {
  
  /** \brief Default implementation 
      
     This class carries two ostreams used in SIMD applications.  The
     first, an ostream only prints from a single process that is
     designated as the print process. Teh second is an ostream that
     will print from all processes.  Classes can inherit off this base
     class to give all objects a common look and feel for output.
  */
  class OutputStreamDefaultImpl : public panzer::OutputStream {

  public:

    OutputStreamDefaultImpl();

    ~OutputStreamDefaultImpl();    

    void setOStream(const Teuchos::RCP<Teuchos::FancyOStream>& os);

    Teuchos::RCP<Teuchos::FancyOStream> getOStream() const;

    Teuchos::FancyOStream& out() const;

    Teuchos::FancyOStream& pout() const;

    void setVerbosityLevel(EVerbosityLevel vl);

    EVerbosityLevel getVerbosityLevel() const;

    bool doOutput(EVerbosityLevel vl, bool only_for_exact_level = false) const;

  private:

    Teuchos::RCP<Teuchos::FancyOStream> m_out;
    Teuchos::RCP<Teuchos::FancyOStream> m_pout;
    EVerbosityLevel m_level;

  };

}

#endif
