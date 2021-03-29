#pragma once
#ifndef ROL2_HPP
#define ROL2_HPP

// C++ Standard Library Includes
//
//#include <algorithm>
//#include <array>
//#include <cmath>
//#include <cstdint>
//#include <exception>
//#include <initializer_list>
//#include <iostream>
//#include <limits>
//#include <map>
//#include <memory>
//#include <sstream>
//#include <stdexcept>
//#include <string>
//#include <tuple>
//#include <type_traits>
//#include <utility>
//#include <vector>


// Common ROL2 Components Includes

#include "ROL2_Utilities.hpp"
#include "ROL2_Elementwise.hpp"


// Forward Declarations
namespace ROL2 {
template<class> class StatusTest;
namespace TypeU {
template<class> class TrustRegionModel;
} // namespace TypeU
} // namespace ROL2

// ROL2 Class Declarations

#include "ROL2_UpdateType.hpp"
#include "ROL2_Vector_Decl.hpp"
#include "ROL2_StdVector_Decl.hpp"
#include "ROL2_Objective_Decl.hpp"

#include "ROL2_Algorithm_Decl.hpp"
#include "ROL2_StatusTest_Decl.hpp"
#include "ROL2_CombinedStatusTest_Decl.hpp"
#include "ROL2_TypeU_Algorithm_Decl.hpp"

#include "ROL2_LinearOperator.hpp"

#include "ROL2_Secant_Decl.hpp"
#include "ROL2_lBFGS_Decl.hpp"
#include "ROL2_lSR1_Decl.hpp"
#include "ROL2_lDFP_Decl.hpp"
#include "ROL2_BarzilaiBorwein_Decl.hpp"

#include "ROL2_Krylov_Decl.hpp"
#include "ROL2_ConjugateGradients_Decl.hpp"
#include "ROL2_ConjugateResiduals_Decl.hpp"

#include "ROL2_TypeU_TrustRegion_Decl.hpp"
#include "ROL2_TypeU_TrustRegionModel_Decl.hpp"
#include "ROL2_TypeU_TruncatedCG_Decl.hpp"
#include "ROL2_TypeU_CauchyPoint_Decl.hpp"
#include "ROL2_TypeU_DogLeg_Decl.hpp"
#include "ROL2_TypeU_DoubleDogLeg_Decl.hpp"
#include "ROL2_TypeU_SPGTrustRegion_Decl.hpp"

#include "ROL2_TypeU_TrustRegionAlgorithm_Decl.hpp"

#include "ROL2_TypeU_DescentDirection_Decl.hpp"
#include "ROL2_TypeU_Gradient_Decl.hpp"
#include "ROL2_TypeU_Newton_Decl.hpp"
#include "ROL2_TypeU_QuasiNewton_Decl.hpp"
#include "ROL2_TypeU_NewtonKrylov_Decl.hpp"
// #include "ROL2_TypeU_NonlinearCG_Decl.hpp"

#include "ROL2_TypeU_LineSearch_Decl.hpp"
#include "ROL2_TypeU_BackTracking_Decl.hpp"
#include "ROL2_TypeU_CubicInterp_Decl.hpp"
#include "ROL2_TypeU_IterationScaling_Decl.hpp"
#include "ROL2_TypeU_PathBasedTargetLevel_Decl.hpp"
//#include "ROL2_TypeU_ScalarMinimizationLineSearch_Decl.hpp"


// ROL2 Class Definitions

#include "ROL2_Vector_Def.hpp"
#include "ROL2_StdVector_Def.hpp"
#include "ROL2_Objective_Def.hpp"

#include "ROL2_Secant_Def.hpp"
#include "ROL2_lBFGS_Def.hpp"
#include "ROL2_lSR1_Def.hpp"
#include "ROL2_lDFP_Def.hpp"
#include "ROL2_BarzilaiBorwein_Def.hpp"

#include "ROL2_Krylov_Def.hpp"
#include "ROL2_ConjugateGradients_Def.hpp"
#include "ROL2_ConjugateResiduals_Def.hpp"

#include "ROL2_StatusTest_Def.hpp"
#include "ROL2_CombinedStatusTest_Def.hpp"

#include "ROL2_Algorithm_Def.hpp"
#include "ROL2_TypeU_Algorithm_Def.hpp"

#include "ROL2_TypeU_TrustRegionModel_Def.hpp"
#include "ROL2_TypeU_TrustRegionAlgorithm_Def.hpp"
#include "ROL2_TypeU_TruncatedCG_Def.hpp"
#include "ROL2_TypeU_CauchyPoint_Def.hpp"
#include "ROL2_TypeU_DogLeg_Def.hpp"
#include "ROL2_TypeU_DoubleDogLeg_Def.hpp"
#include "ROL2_TypeU_SPGTrustRegion_Def.hpp"
#include "ROL2_TypeU_TrustRegion_Def.hpp"

#include "ROL2_TypeU_DescentDirection_Def.hpp"
#include "ROL2_TypeU_Gradient_Def.hpp"
#include "ROL2_TypeU_Newton_Def.hpp"
#include "ROL2_TypeU_QuasiNewton_Def.hpp"
#include "ROL2_TypeU_NewtonKrylov_Def.hpp"

#include "ROL2_TypeU_DescentDirection_Def.hpp"
#include "ROL2_TypeU_Gradient_Def.hpp"
#include "ROL2_TypeU_Newton_Def.hpp"
#include "ROL2_TypeU_QuasiNewton_Def.hpp"
#include "ROL2_TypeU_NewtonKrylov_Def.hpp"
// #include "ROL2_TypeU_NonlinearCG_Def.hpp"

#include "ROL2_TypeU_LineSearch_Def.hpp"
#include "ROL2_TypeU_BackTracking_Def.hpp"
#include "ROL2_TypeU_CubicInterp_Def.hpp"
#include "ROL2_TypeU_IterationScaling_Def.hpp"
#include "ROL2_TypeU_PathBasedTargetLevel_Def.hpp"
//#include "ROL2_TypeU_ScalarMinimizationLineSearch_Def.hpp"



#endif // ROL2_HPP

