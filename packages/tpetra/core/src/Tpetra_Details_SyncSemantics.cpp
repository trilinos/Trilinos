/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// ************************************************************************
// @HEADER
*/
#include "Tpetra_Details_SyncSemantics.hpp"
#include <stdexcept>

namespace Tpetra {
namespace Details {

  namespace SyncSematicsDetails {
    // Static variables
    bool useRelaxedSyncs = false;
  }// end SyncSemanticsDetails


  void enableRelaxedSyncs() {
    if(SyncSematicsDetails::useRelaxedSyncs)
      throw std::runtime_error("ERROR: Tpetra::Details::enableRelaxedSyncs(): Cannot call twice in a row");

    SyncSematicsDetails::useRelaxedSyncs = true;
    Kokkos::fence("Tpetra::enableRelaxedSyncs");
  }

  void disableRelaxedSyncs() {
    if(!SyncSematicsDetails::useRelaxedSyncs)
      throw std::runtime_error("ERROR: Tpetra::Details::disableRelaxedSyncs(): Cannot call twice in a row");

    SyncSematicsDetails::useRelaxedSyncs = false;
    // Kokkos::fence("Tpetra::disableRelaxedSyncs");
  }

  bool areRelaxedSyncsEnabled() {
    return SyncSematicsDetails::useRelaxedSyncs;
  }


  void fence() {
    if(areRelaxedSyncsEnabled()) {
      Kokkos::DefaultExecutionSpace().fence();
    }
    else {
      Kokkos::fence();
    }

  }

  void fence(const std::string& label) {
    if(areRelaxedSyncsEnabled()) {
      Kokkos::DefaultExecutionSpace().fence(label);
    }
    else {
      Kokkos::fence(label);
    }
  }

} // namespace Details


} // namespace Tpetra
