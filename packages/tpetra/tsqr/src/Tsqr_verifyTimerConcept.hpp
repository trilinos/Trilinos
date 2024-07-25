// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __TSQR_Test_verifyTimerConcept_hpp
#define __TSQR_Test_verifyTimerConcept_hpp

#include <stdexcept>

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
