#ifndef UTILS_ENUMS_H
#define UTILS_ENUMS_H

#include "LOCA_Utils.H"

class Utils {
 public:
  enum MsgType { 
    Error =                 LOCA::Utils::Error,
    Warning =               LOCA::Utils::Warning,
    StepperIteration =      LOCA::Utils::StepperIteration,
    StepperDetails =        LOCA::Utils::StepperDetails,
    Solver =                LOCA::Utils::Solver,
    SolverDetails =         LOCA::Utils::SolverDetails,
    Direction =             LOCA::Utils::Direction,
    Parameters =            LOCA::Utils::Parameters
  };
};

#endif
