// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
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
