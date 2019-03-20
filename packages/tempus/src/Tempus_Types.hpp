// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_Types_hpp
#define Tempus_Types_hpp

namespace Tempus {


/** \brief Status for the Integrator, the Stepper and the SolutionState */
enum Status {
  PASSED,
  FAILED,
  WORKING
};


/** \brief Convert Status to string. */
inline
const std::string toString(const Status status)
{
  std::string s = "Invalid Status!";
  switch(status) {
    case PASSED:  { s = "PASSED";          break; }
    case FAILED:  { s = "FAILED";          break; }
    case WORKING: { s = "WORKING";         break; }
    default:      { s = "Invalid Status!"; break; }
  }
  TEUCHOS_TEST_FOR_EXCEPTION(s == "Invalid Status!", std::logic_error,
   "Error - Invalid status = " << status << "\n");
  return s;
}


} // namespace Tempus
#endif // Tempus_Types_hpp
