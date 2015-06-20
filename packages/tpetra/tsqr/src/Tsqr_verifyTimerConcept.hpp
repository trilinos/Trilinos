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

#ifndef __TSQR_Test_verifyTimerConcept_hpp
#define __TSQR_Test_verifyTimerConcept_hpp

#include <stdexcept>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Test {
    /// \function verifyTimerConcept
    /// \brief Ensure that TimerType has the required interface.
    ///
    /// Our TSQR benchmarks are templated on TimerType, in order to
    /// avoid depending on a particular timer implementation.
    /// TimerType should support the following interface (modeled after
    /// the \c Teuchos::Time class):
    ///
    /// - Construction using a const std::string& or something
    ///   convertible to that (to name the timer).  Semantically, the
    ///   constructor should not start the timer.
    /// - name(): Returns the name of the timer, as it was set in the
    ///   constructor (it does not change the string).
    /// - start(bool reset=false): Returns nothing, starts the timer,
    ///   resets it first if reset==true).
    /// - stop(): Returns a double-precision floating-point value,
    ///   which is the number of seconds elapsed since calling
    ///   start().
    /// - isRunning(): Returns a Boolean saying whether the timer is
    ///   currently running.  start() should make the timer start
    ///   running, and stop() shoudl make it stop running.
    ///
    /// TimerType need not be able to handle recursive calls, though
    /// this might be helpful.  The intended use case is that start()
    /// and stop() wrap some timing loop, and the loop does not
    /// reference the timer at all.  We include this concept check in
    /// all of our TSQR benchmark routines via
    ///
    /// \code
    /// verifyTimerConcept<TimerType>();
    /// \endcode
    ///
    /// If TimerType does not satisfy this interface, that line of
    /// code will fail to compile.  The compiler should give an
    /// informative error message about a missing method.
    /// verifyTimerConcept() also checks some semantic properties of
    /// TimerType.
    template<class TimerType>
    double
    verifyTimerConcept ()
    {
      TimerType timer (std::string("NameOfTimer"));

      std::string timerName = timer.name();
      if (timerName != "NameOfTimer")
	throw std::logic_error ("TimerType does not correctly store the timer name");

      // Test default argument of start()
      if (timer.isRunning())
	throw std::logic_error ("TimerType does not correctly initialize isRunning");
      timer.start ();
      if (! timer.isRunning())
	throw std::logic_error ("TimerType does not correctly set isRunning");
      double result1 = timer.stop();
      if (timer.isRunning())
	throw std::logic_error ("TimerType does not correctly reset isRunning");

      // Test nondefault argument of start()
      if (timer.isRunning())
	throw std::logic_error ("TimerType does not correctly initialize isRunning");
      timer.start (true);
      if (! timer.isRunning())
	throw std::logic_error ("TimerType does not correctly set isRunning");
      double result2 = timer.stop();
      if (timer.isRunning())
	throw std::logic_error ("TimerType does not correctly reset isRunning");

      return result1 + result2;
    }
  } // namespace Test
} // namespace TSQR

#endif // __TSQR_Test_verifyTimerConcept_hpp
