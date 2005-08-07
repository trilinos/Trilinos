#ifndef UTILS_ENUMS_H
#define UTILS_ENUMS_H

#include "NOX_Utils.H"

// The NOX::Utils class utilizes nested classes, which cannot currently be
// wrapped by SWIG.  Therefore, we pull out the needed enumerations that are
// used by the NOX::Parameter::List.

class Utils {
public:
  enum MsgType { 
    Error                    = NOX::Utils::Error,
    Warning                  = NOX::Utils::Warning,
    OuterIteration           = NOX::Utils::OuterIteration,
    InnerIteration           = NOX::Utils::InnerIteration,
    Parameters               = NOX::Utils::Parameters,
    Details                  = NOX::Utils::Details,
    OuterIterationStatusTest = NOX::Utils::OuterIterationStatusTest,
    LinearSolverDetails      = NOX::Utils::LinearSolverDetails,
    TestDetails		     = NOX::Utils::TestDetails
  };
};

#endif
