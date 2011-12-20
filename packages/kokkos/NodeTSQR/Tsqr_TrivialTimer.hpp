/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
//@HEADER
*/

#ifndef __TSQR_TrivialTimer_hpp
#define __TSQR_TrivialTimer_hpp

#include <string>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  /// \class TrivialTimer
  /// \brief Satisfies TimerType concept trivially.
  ///
  /// This is a "prototype" for the TimerType concept; it satisfies
  /// the concept trivially.
  class TrivialTimer {
  public:
    /// Constructor.
    ///
    /// \param name [in] Timer label
    /// \param doStart [in] Whether to start timer on instantiation
    TrivialTimer (const std::string& theName, bool doStart = false);

    /// \brief Start the timer.
    ///
    /// This is a trivial timer, so this implementation does not
    /// actually return valid times.  However, it satisfies our
    /// TimerType concept.
    void start (bool reset = false);
    
    //! Stop the timer and return (fake) elapsed time.
    double stop ();

    //! Whether this timer is running
    bool isRunning () const { return isRunning_; }

    //! Name of this timer object
    const std::string& name() const { return name_; }

  private:
    //! Name of this timer
    std::string name_;

    /// \brief Counter, used to implement trivial timer feature.
    ///
    /// In order to compute the "resolution" of this fake timer, we
    /// add a counter which is incremented on each call to \c stop().
    /// The \c stop() method computes a fake timing result based on
    /// the counter value.
    size_t counter_;
    
    //! Whether this timer is running
    bool isRunning_;

    //! Verify the TimerType concept
    static void verifyConcept();
  };

} // namespace TSQR

#endif // __TSQR_TrivialTimer_hpp
