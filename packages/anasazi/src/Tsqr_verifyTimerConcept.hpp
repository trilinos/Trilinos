#ifndef __TSQR_Test_verifyTimerConcept_hpp
#define __TSQR_Test_verifyTimerConcept_hpp

#include <string>

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
/// \li start() (returns nothing, starts the timer)
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
    void 
    verifyTimerConcept ()
    {
      TimerType timer (std::string("NameOfTimer"));
      timer.start ();
      double result = timer.stop();
    }
  } // namespace Test
} // namespace TSQR

#endif // __TSQR_Test_verifyTimerConcept_hpp
