// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_OUTPUT_STREAM_HPP
#define PANZER_OUTPUT_STREAM_HPP

#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_RCP.hpp"

namespace panzer {
  
  enum EVerbosityLevel {
    VERB_NONE,
    VERB_LOW,
    VERB_MEDIUM,
    VERB_HIGH,
    VERB_EXTREME
  };

  /** \brief Interface for handling output in Panzer
      
     This class carries two ostreams used in SIMD applications.  The
     first (out), an ostream that only prints on a specific process
     that is designated as the print process. The second is an ostream
     (pout) that will print from all processes.  Classes can inherit
     off this base class to give all objects a common look and feel
     for output.
  */
  class OutputStream {

  public:

    virtual ~OutputStream() {}

    virtual void setOStream(const Teuchos::RCP<Teuchos::FancyOStream>& os) = 0;

    virtual Teuchos::RCP<Teuchos::FancyOStream> getOStream() const = 0;

    /** returns ostream that prints only to print process */
    virtual Teuchos::FancyOStream& out() const = 0;

    /** returns ostream that prints on all processes */
    virtual Teuchos::FancyOStream& pout() const = 0;

    virtual void setVerbosityLevel(EVerbosityLevel vl) = 0;

    virtual EVerbosityLevel getVerbosityLevel() const = 0;

    /** \brief Returns true if vl is equal to or greater than the object's verbosity level

    \param vl [in] Verbosity level for comparison
    \param only_for_exact_level [in] Forces the compaison to be the exact vebosity level instead of equal to or greater than
    */
    virtual bool doOutput(EVerbosityLevel vl, bool only_for_exact_level = false) const = 0;

  };

}

#endif
