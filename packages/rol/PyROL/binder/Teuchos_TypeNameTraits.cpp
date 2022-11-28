#include <ROL_AffineTransformConstraint.hpp>
#include <ROL_AffineTransformObjective.hpp>
#include <ROL_AugmentedLagrangianObjective.hpp>
#include <ROL_BackTracking_U.hpp>
#include <ROL_BarzilaiBorwein.hpp>
#include <ROL_BisectionScalarMinimization.hpp>
#include <ROL_BoundConstraint.hpp>
#include <ROL_BoundConstraint_Partitioned.hpp>
#include <ROL_Bounds.hpp>
#include <ROL_Bracketing.hpp>
#include <ROL_BrentsProjection.hpp>
#include <ROL_BrentsScalarMinimization.hpp>
#include <ROL_BundleStatusTest.hpp>
#include <ROL_Bundle_U_AS.hpp>
#include <ROL_Bundle_U_TT.hpp>
#include <ROL_CauchyPoint_U.hpp>
#include <ROL_CombinedStatusTest.hpp>
#include <ROL_ConjugateGradients.hpp>
#include <ROL_ConjugateResiduals.hpp>
#include <ROL_ConstraintStatusTest.hpp>
#include <ROL_Constraint_Partitioned.hpp>
#include <ROL_CubicInterp_U.hpp>
#include <ROL_DaiFletcherProjection.hpp>
#include <ROL_DescentDirection_U.hpp>
#include <ROL_DogLeg_U.hpp>
#include <ROL_DoubleDogLeg_U.hpp>
#include <ROL_DouglasRachfordProjection.hpp>
#include <ROL_DykstraProjection.hpp>
#include <ROL_ElasticLinearConstraint.hpp>
#include <ROL_ElasticObjective.hpp>
#include <ROL_Elementwise_Function.hpp>
#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_FletcherObjectiveE.hpp>
#include <ROL_GMRES.hpp>
#include <ROL_GoldenSectionScalarMinimization.hpp>
#include <ROL_Gradient_U.hpp>
#include <ROL_IterationScaling_U.hpp>
#include <ROL_Krylov.hpp>
#include <ROL_LineSearch_U.hpp>
#include <ROL_LineSearch_U_Types.hpp>
#include <ROL_LinearConstraint.hpp>
#include <ROL_LinearOperator.hpp>
#include <ROL_MINRES.hpp>
#include <ROL_NewtonKrylov_U.hpp>
#include <ROL_Newton_U.hpp>
#include <ROL_NonlinearCG.hpp>
#include <ROL_NonlinearCG_U.hpp>
#include <ROL_NullSpaceOperator.hpp>
#include <ROL_PQNObjective.hpp>
#include <ROL_PartitionedVector.hpp>
#include <ROL_PathBasedTargetLevel_U.hpp>
#include <ROL_PolyhedralProjection.hpp>
#include <ROL_Problem.hpp>
#include <ROL_QuasiNewton_U.hpp>
#include <ROL_ReduceLinearConstraint.hpp>
#include <ROL_ReducedLinearConstraint.hpp>
#include <ROL_RiddersProjection.hpp>
#include <ROL_SPGTrustRegion_U.hpp>
#include <ROL_ScalarController.hpp>
#include <ROL_ScalarFunction.hpp>
#include <ROL_ScalarMinimization.hpp>
#include <ROL_ScalarMinimizationLineSearch_U.hpp>
#include <ROL_ScalarMinimizationStatusTest.hpp>
#include <ROL_Secant.hpp>
#include <ROL_SemismoothNewtonProjection.hpp>
#include <ROL_SingletonVector.hpp>
#include <ROL_SlacklessObjective.hpp>
#include <ROL_StatusTest.hpp>
#include <ROL_Stream.hpp>
#include <ROL_TruncatedCG_U.hpp>
#include <ROL_TrustRegionModel_U.hpp>
#include <ROL_TrustRegion_U_Types.hpp>
#include <ROL_TypeB_Algorithm.hpp>
#include <ROL_TypeB_ColemanLiAlgorithm.hpp>
#include <ROL_TypeB_GradientAlgorithm.hpp>
#include <ROL_TypeB_InteriorPointAlgorithm.hpp>
#include <ROL_TypeB_KelleySachsAlgorithm.hpp>
#include <ROL_TypeB_LSecantBAlgorithm.hpp>
#include <ROL_TypeB_LinMoreAlgorithm.hpp>
#include <ROL_TypeB_MoreauYosidaAlgorithm.hpp>
#include <ROL_TypeB_NewtonKrylovAlgorithm.hpp>
#include <ROL_TypeB_PrimalDualActiveSetAlgorithm.hpp>
#include <ROL_TypeB_QuasiNewtonAlgorithm.hpp>
#include <ROL_TypeB_SpectralGradientAlgorithm.hpp>
#include <ROL_TypeB_TrustRegionSPGAlgorithm.hpp>
#include <ROL_TypeE_Algorithm.hpp>
#include <ROL_TypeE_AugmentedLagrangianAlgorithm.hpp>
#include <ROL_TypeE_CompositeStepAlgorithm.hpp>
#include <ROL_TypeE_FletcherAlgorithm.hpp>
#include <ROL_TypeE_StabilizedLCLAlgorithm.hpp>
#include <ROL_TypeG_Algorithm.hpp>
#include <ROL_TypeG_AugmentedLagrangianAlgorithm.hpp>
#include <ROL_TypeG_InteriorPointAlgorithm.hpp>
#include <ROL_TypeG_MoreauYosidaAlgorithm.hpp>
#include <ROL_TypeG_StabilizedLCLAlgorithm.hpp>
#include <ROL_TypeU_Algorithm.hpp>
#include <ROL_TypeU_BundleAlgorithm.hpp>
#include <ROL_TypeU_LineSearchAlgorithm.hpp>
#include <ROL_TypeU_TrustRegionAlgorithm.hpp>
#include <ROL_Types.hpp>
#include <ROL_UpdateType.hpp>
#include <ROL_VectorController.hpp>
#include <ROL_lBFGS.hpp>
#include <ROL_lDFP.hpp>
#include <ROL_lSR1.hpp>
#include <Teuchos_BLAS_types.hpp>
#include <Teuchos_DataAccess.hpp>
#include <Teuchos_ENull.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_FilteredIterator.hpp>
#include <Teuchos_ParameterEntry.hpp>
#include <Teuchos_ParameterEntryValidator.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ParameterListModifier.hpp>
#include <Teuchos_RCPDecl.hpp>
#include <Teuchos_RCPNode.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_StringIndexedOrderedValueObjectContainer.hpp>
#include <Teuchos_TestForException.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_any.hpp>
#include <cwchar>
#include <deque>
#include <ios>
#include <iterator>
#include <locale>
#include <memory>
#include <ostream>
#include <sstream> // __str__
#include <streambuf>
#include <string>
#include <typeinfo>
#include <vector>

#include <functional>
#include <pybind11/pybind11.h>
#include <string>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <Teuchos_RCP.hpp>


#ifndef BINDER_PYBIND11_TYPE_CASTER
	#define BINDER_PYBIND11_TYPE_CASTER
	PYBIND11_DECLARE_HOLDER_TYPE(T, Teuchos::RCP<T>)
	PYBIND11_DECLARE_HOLDER_TYPE(T, T*)
	PYBIND11_MAKE_OPAQUE(Teuchos::RCP<void>)
#endif

void bind_Teuchos_TypeNameTraits(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	// Teuchos::demangleName(const std::string &) file:Teuchos_TypeNameTraits.hpp line:77
	M("Teuchos").def("demangleName", (std::string (*)(const std::string &)) &Teuchos::demangleName, "Demangle a C++ name if valid.\n\n The name must have come from typeid(...).name() in order to be\n valid name to pass to this function.\n\n \n\n \n\nC++: Teuchos::demangleName(const std::string &) --> std::string", pybind11::arg("mangledName"));

	// Teuchos::typeName(const class Teuchos::any::placeholder &) file:Teuchos_TypeNameTraits.hpp line:115
	M("Teuchos").def("typeName", (std::string (*)(const class Teuchos::any::placeholder &)) &Teuchos::typeName<Teuchos::any::placeholder>, "C++: Teuchos::typeName(const class Teuchos::any::placeholder &) --> std::string", pybind11::arg("t"));

	// Teuchos::TestForException_incrThrowNumber() file:Teuchos_TestForException.hpp line:61
	M("Teuchos").def("TestForException_incrThrowNumber", (void (*)()) &Teuchos::TestForException_incrThrowNumber, "Increment the throw number.  \n\nC++: Teuchos::TestForException_incrThrowNumber() --> void");

	// Teuchos::TestForException_getThrowNumber() file:Teuchos_TestForException.hpp line:64
	M("Teuchos").def("TestForException_getThrowNumber", (int (*)()) &Teuchos::TestForException_getThrowNumber, "Increment the throw number.  \n\nC++: Teuchos::TestForException_getThrowNumber() --> int");

	// Teuchos::TestForException_break(const std::string &) file:Teuchos_TestForException.hpp line:68
	M("Teuchos").def("TestForException_break", (void (*)(const std::string &)) &Teuchos::TestForException_break, "The only purpose for this function is to set a breakpoint.\n    \n\n\nC++: Teuchos::TestForException_break(const std::string &) --> void", pybind11::arg("msg"));

	// Teuchos::TestForException_setEnableStacktrace(bool) file:Teuchos_TestForException.hpp line:72
	M("Teuchos").def("TestForException_setEnableStacktrace", (void (*)(bool)) &Teuchos::TestForException_setEnableStacktrace, "Set at runtime if stacktracing functionality is enabled when *\n    exceptions are thrown.  \n\n\nC++: Teuchos::TestForException_setEnableStacktrace(bool) --> void", pybind11::arg("enableStrackTrace"));

	// Teuchos::TestForException_getEnableStacktrace() file:Teuchos_TestForException.hpp line:76
	M("Teuchos").def("TestForException_getEnableStacktrace", (bool (*)()) &Teuchos::TestForException_getEnableStacktrace, "Get at runtime if stacktracing functionality is enabled when\n exceptions are thrown. \n\nC++: Teuchos::TestForException_getEnableStacktrace() --> bool");

	// Teuchos::TestForTermination_terminate(const std::string &) file:Teuchos_TestForException.hpp line:79
	M("Teuchos").def("TestForTermination_terminate", (void (*)(const std::string &)) &Teuchos::TestForTermination_terminate, "Prints the message to std::cerr and calls std::terminate. \n\nC++: Teuchos::TestForTermination_terminate(const std::string &) --> void", pybind11::arg("msg"));

}
