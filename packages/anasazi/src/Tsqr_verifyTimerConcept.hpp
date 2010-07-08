#ifndef __TSQR_Test_verifyTimerConcept_hpp
#define __TSQR_Test_verifyTimerConcept_hpp

#include <ostream>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/// \function verifyTimerConcept
/// \brief Make sure TimerType has the required interface
///
/// Our TSQR benchmarks are templated on TimerType, in order to avoid
/// depending on a particular timer implementation.  TimerType should
/// support the following three methods:
///
/// \li Default construction
/// \li void start() 
/// \li double finish()
///
/// We include this concept check in all of our TSQR benchmark
/// routines via 
///
/// verifyTimerConcept< TimerType >();
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
      TimerType timer;
      timer.start ();
      double result = timer.finish();
    }
  } // namespace Test
} // namespace TSQR

#endif // __TSQR_Test_verifyTimerConcept_hpp
