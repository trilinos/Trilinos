//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef __TSQR_Test_verifyTimerConcept_hpp
#define __TSQR_Test_verifyTimerConcept_hpp

#include <stdexcept>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/// \function verifyTimerConcept
/// \brief Make sure TimerType has the required interface
///
/// Our TSQR benchmarks are templated on TimerType, in order to avoid
/// depending on a particular timer implementation.  TimerType should
/// support the following three methods (modeled after Trilinos'
/// Teuchos::Time class):
///
/// \li Construction using a const std::string& or something 
///     convertible to that (to name the timer).  Semantically,
///     the constructor should not start the timer.
/// \li start(bool reset=false) (returns nothing, starts the timer)
/// \li stop() (returns a double-precision floating-point value,
///     which is the number of seconds elapsed since calling start())
///
/// TimerType need not be able to handle recursive calls; the model
/// is that start() and stop() wrap some timing loop, and the loop 
/// does not reference the timer at all.  We include this concept 
/// check in all of our TSQR benchmark routines via 
///
/// \code
/// verifyTimerConcept< TimerType >();
/// \endcode
///
/// If TimerType does not satisfy this interface, that line of code
/// will fail to compile.  The compiler should give an informative
/// error message about a missing method.
namespace TSQR {
  namespace Test {
    template< class TimerType >
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
