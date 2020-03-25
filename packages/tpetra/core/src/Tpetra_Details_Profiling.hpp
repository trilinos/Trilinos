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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_PROFILING_HPP
#define TPETRA_DETAILS_PROFILING_HPP

/// \file Tpetra_Details_Profiling.hpp
/// \brief Declaration of Tpetra::Details::Profiling, a scope guard
///   for Kokkos Profiling.

#include "TpetraCore_config.h"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_RCP.hpp"


namespace Tpetra {
namespace Details {

/// \brief Profile the given scope.
///
/// This class implements a "scope guard" for profiling a region of
/// code, using Kokkos Profiling.
///
/// If you disable Kokkos Profiling by setting some configuration
/// option such that KOKKOS_ENABLE_PROFILING is not defined, then this
/// class is harmless, but will do NOTHING.
///
/// Tpetra currently uses Kokkos::Profiling::pushRegion(name)
/// to start the region to profile, and Kokkos::Profiling::popRegion()
/// to close the region to profile.  The "scope guard" idiom means
/// that the destructor automatically closes the region, no matter how
/// the scope exits, even if an exception is thrown.
///
/// Here is an example of how to use this class:
/// \code
/// void someFunctionYouWrote () {
///   ProfilingRegion wholeFuncRegion ("someFunctionYouWrote");
///   doStuff ();
///   {
///     ProfilingRegion innerRegion1 ("inner-region-1");
///     doOtherStuff ();
///   }
///   doSomeMoreStuff ();
///   {
///     ProfilingRegion innerRegion2 ("yet-another-inner-region");
///     doYetMoreStuff ();
///   }
///   finallyDoEvenStillMoreStuff ();
/// }
/// \endcode
///
/// It's critical that you not forget to name the ProfilingRegion
/// variable.  If you forget, you will create a temporary
/// ProfilingRegion object that will first open, then immediately
/// close the "region" before anything else happens.  This is not what
/// you want.
///
/// You may name regions whatever you like.  They don't have to have
/// the same name as the region variable.  Kokkos Profiling has a
/// length limit of 512 characters for region names.
///
/// For details about Kokkos Profiling, see the kokkos-tools wiki:
///
/// https://github.com/kokkos/kokkos-tools/wiki
class ProfilingRegion {
public:
  //! Open region to profile; name the region \c name.
  ProfilingRegion (const char name[]);
  //! Close region to profile.
  ~ProfilingRegion ();

private:
  Teuchos::RCP<Teuchos::TimeMonitor> tm;

};

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_PROFILING_HPP
