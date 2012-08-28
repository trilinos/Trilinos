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
