#include <PyROL_Teuchos_Custom.hpp>
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
#include <ROL_CauchyPoint.hpp>
#include <ROL_CauchyPoint_U.hpp>
#include <ROL_ColemanLiModel.hpp>
#include <ROL_CombinedStatusTest.hpp>
#include <ROL_ConjugateGradients.hpp>
#include <ROL_ConjugateResiduals.hpp>
#include <ROL_Constraint.hpp>
#include <ROL_ConstraintStatusTest.hpp>
#include <ROL_Constraint_Partitioned.hpp>
#include <ROL_Constraint_SimOpt.hpp>
#include <ROL_CubicInterp_U.hpp>
#include <ROL_DaiFletcherProjection.hpp>
#include <ROL_DescentDirection_U.hpp>
#include <ROL_DogLeg.hpp>
#include <ROL_DogLeg_U.hpp>
#include <ROL_DoubleDogLeg.hpp>
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
#include <ROL_KelleySachsModel.hpp>
#include <ROL_Krylov.hpp>
#include <ROL_LinMore.hpp>
#include <ROL_LinMoreModel.hpp>
#include <ROL_LineSearch_U.hpp>
#include <ROL_LinearConstraint.hpp>
#include <ROL_LinearOperator.hpp>
#include <ROL_NewtonKrylov_U.hpp>
#include <ROL_Newton_U.hpp>
#include <ROL_NonlinearCG.hpp>
#include <ROL_NonlinearCG_U.hpp>
#include <ROL_NonlinearLeastSquaresObjective.hpp>
#include <ROL_NullSpaceOperator.hpp>
#include <ROL_Objective.hpp>
#include <ROL_Objective_FSsolver.hpp>
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
#include <ROL_SimConstraint.hpp>
#include <ROL_SingletonVector.hpp>
#include <ROL_SlacklessObjective.hpp>
#include <ROL_StatusTest.hpp>
#include <ROL_TruncatedCG.hpp>
#include <ROL_TruncatedCG_U.hpp>
#include <ROL_TrustRegionModel.hpp>
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
#include <ROL_Vector.hpp>
#include <ROL_VectorController.hpp>
#include <ROL_Vector_SimOpt.hpp>
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
#include <Teuchos_PtrDecl.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_RCPDecl.hpp>
#include <Teuchos_RCPNode.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_StringIndexedOrderedValueObjectContainer.hpp>
#include <Teuchos_any.hpp>
#include <cwchar>
#include <deque>
#include <ios>
#include <iterator>
#include <locale>
#include <memory>
#include <ostream>
#include <streambuf>
#include <string>
#include <typeinfo>
#include <vector>

#include <functional>
#include <pybind11/smart_holder.h>
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

void bind_Teuchos_RCP(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	// Teuchos::RCP_createNewRCPNodeRawPtrNonowned(std::ostream *) file:Teuchos_RCP.hpp line:75
	M("Teuchos").def("RCP_createNewRCPNodeRawPtrNonowned", (class Teuchos::RCPNode * (*)(std::ostream *)) &Teuchos::RCP_createNewRCPNodeRawPtrNonowned<std::ostream>, "C++: Teuchos::RCP_createNewRCPNodeRawPtrNonowned(std::ostream *) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"));

	// Teuchos::RCP_createNewRCPNodeRawPtrNonowned(class Teuchos::ParameterEntry *) file:Teuchos_RCP.hpp line:75
	M("Teuchos").def("RCP_createNewRCPNodeRawPtrNonowned", (class Teuchos::RCPNode * (*)(class Teuchos::ParameterEntry *)) &Teuchos::RCP_createNewRCPNodeRawPtrNonowned<Teuchos::ParameterEntry>, "C++: Teuchos::RCP_createNewRCPNodeRawPtrNonowned(class Teuchos::ParameterEntry *) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"));

	// Teuchos::RCP_createNewRCPNodeRawPtrNonowned(const class Teuchos::ParameterEntry *) file:Teuchos_RCP.hpp line:75
	M("Teuchos").def("RCP_createNewRCPNodeRawPtrNonowned", (class Teuchos::RCPNode * (*)(const class Teuchos::ParameterEntry *)) &Teuchos::RCP_createNewRCPNodeRawPtrNonowned<const Teuchos::ParameterEntry>, "C++: Teuchos::RCP_createNewRCPNodeRawPtrNonowned(const class Teuchos::ParameterEntry *) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"));

	// Teuchos::RCP_createNewRCPNodeRawPtrNonowned(const class ROL::Vector<double> *) file:Teuchos_RCP.hpp line:75
	M("Teuchos").def("RCP_createNewRCPNodeRawPtrNonowned", (class Teuchos::RCPNode * (*)(const class ROL::Vector<double> *)) &Teuchos::RCP_createNewRCPNodeRawPtrNonowned<const ROL::Vector<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtrNonowned(const class ROL::Vector<double> *) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"));

	// Teuchos::RCP_createNewRCPNodeRawPtrNonowned(class ROL::Constraint_SimOpt<double> *) file:Teuchos_RCP.hpp line:75
	M("Teuchos").def("RCP_createNewRCPNodeRawPtrNonowned", (class Teuchos::RCPNode * (*)(class ROL::Constraint_SimOpt<double> *)) &Teuchos::RCP_createNewRCPNodeRawPtrNonowned<ROL::Constraint_SimOpt<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtrNonowned(class ROL::Constraint_SimOpt<double> *) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"));

	// Teuchos::RCP_createNewRCPNodeRawPtrNonowned(class ROL::Objective<double> *) file:Teuchos_RCP.hpp line:75
	M("Teuchos").def("RCP_createNewRCPNodeRawPtrNonowned", (class Teuchos::RCPNode * (*)(class ROL::Objective<double> *)) &Teuchos::RCP_createNewRCPNodeRawPtrNonowned<ROL::Objective<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtrNonowned(class ROL::Objective<double> *) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"));

	// Teuchos::RCP_createNewRCPNodeRawPtrNonowned(class ROL::Constraint<double> *) file:Teuchos_RCP.hpp line:75
	M("Teuchos").def("RCP_createNewRCPNodeRawPtrNonowned", (class Teuchos::RCPNode * (*)(class ROL::Constraint<double> *)) &Teuchos::RCP_createNewRCPNodeRawPtrNonowned<ROL::Constraint<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtrNonowned(class ROL::Constraint<double> *) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"));

	// Teuchos::RCP_createNewRCPNodeRawPtrNonowned(class ROL::Vector<double> *) file:Teuchos_RCP.hpp line:75
	M("Teuchos").def("RCP_createNewRCPNodeRawPtrNonowned", (class Teuchos::RCPNode * (*)(class ROL::Vector<double> *)) &Teuchos::RCP_createNewRCPNodeRawPtrNonowned<ROL::Vector<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtrNonowned(class ROL::Vector<double> *) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"));

	// Teuchos::RCP_createNewRCPNodeRawPtrNonowned(class ROL::BoundConstraint<double> *) file:Teuchos_RCP.hpp line:75
	M("Teuchos").def("RCP_createNewRCPNodeRawPtrNonowned", (class Teuchos::RCPNode * (*)(class ROL::BoundConstraint<double> *)) &Teuchos::RCP_createNewRCPNodeRawPtrNonowned<ROL::BoundConstraint<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtrNonowned(class ROL::BoundConstraint<double> *) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"));

	// Teuchos::RCP_createNewRCPNodeRawPtrNonowned(class ROL::ElasticObjective<double> *) file:Teuchos_RCP.hpp line:75
	M("Teuchos").def("RCP_createNewRCPNodeRawPtrNonowned", (class Teuchos::RCPNode * (*)(class ROL::ElasticObjective<double> *)) &Teuchos::RCP_createNewRCPNodeRawPtrNonowned<ROL::ElasticObjective<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtrNonowned(class ROL::ElasticObjective<double> *) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(std::ostream *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(std::ostream *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<std::ostream>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(std::ostream *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class Teuchos::ParameterList *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class Teuchos::ParameterList *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<Teuchos::ParameterList>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class Teuchos::ParameterList *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::Vector_SimOpt<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::Vector_SimOpt<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::Vector_SimOpt<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::Vector_SimOpt<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::SimConstraint<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::SimConstraint<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::SimConstraint<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::SimConstraint<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::NonlinearLeastSquaresObjective<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::NonlinearLeastSquaresObjective<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::NonlinearLeastSquaresObjective<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::NonlinearLeastSquaresObjective<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::CombinedStatusTest<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::CombinedStatusTest<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::CombinedStatusTest<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::CombinedStatusTest<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(struct ROL::TypeU::AlgorithmState<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(struct ROL::TypeU::AlgorithmState<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::TypeU::AlgorithmState<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(struct ROL::TypeU::AlgorithmState<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::StatusTest<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::StatusTest<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::StatusTest<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::StatusTest<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::CauchyPoint_U<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::CauchyPoint_U<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::CauchyPoint_U<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::CauchyPoint_U<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::DogLeg_U<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::DogLeg_U<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::DogLeg_U<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::DogLeg_U<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::DoubleDogLeg_U<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::DoubleDogLeg_U<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::DoubleDogLeg_U<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::DoubleDogLeg_U<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TruncatedCG_U<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::TruncatedCG_U<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::TruncatedCG_U<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TruncatedCG_U<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::SPGTrustRegion_U<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::SPGTrustRegion_U<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::SPGTrustRegion_U<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::SPGTrustRegion_U<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TrustRegionModel_U<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::TrustRegionModel_U<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::TrustRegionModel_U<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TrustRegionModel_U<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::NullSpaceOperator<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::NullSpaceOperator<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::NullSpaceOperator<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::NullSpaceOperator<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeU::TrustRegionAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::TypeU::TrustRegionAlgorithm<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::TypeU::TrustRegionAlgorithm<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeU::TrustRegionAlgorithm<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::Objective_FSsolver<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::Objective_FSsolver<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::Objective_FSsolver<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::Objective_FSsolver<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(struct ROL::TypeE::AlgorithmState<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(struct ROL::TypeE::AlgorithmState<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::TypeE::AlgorithmState<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(struct ROL::TypeE::AlgorithmState<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::ConstraintStatusTest<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::ConstraintStatusTest<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::ConstraintStatusTest<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::ConstraintStatusTest<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::ScalarController<double, int> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::ScalarController<double, int> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::ScalarController<double, int>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::ScalarController<double, int> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::VectorController<double, int> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::VectorController<double, int> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::VectorController<double, int>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::VectorController<double, int> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::SingletonVector<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::SingletonVector<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::SingletonVector<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::SingletonVector<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeE::AugmentedLagrangianAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::TypeE::AugmentedLagrangianAlgorithm<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::TypeE::AugmentedLagrangianAlgorithm<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeE::AugmentedLagrangianAlgorithm<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::BundleStatusTest<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::BundleStatusTest<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::BundleStatusTest<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::BundleStatusTest<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::Bundle_U_TT<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::Bundle_U_TT<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::Bundle_U_TT<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::Bundle_U_TT<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::Bundle_U_AS<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::Bundle_U_AS<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::Bundle_U_AS<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::Bundle_U_AS<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::IterationScaling_U<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::IterationScaling_U<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::IterationScaling_U<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::IterationScaling_U<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::PathBasedTargetLevel_U<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::PathBasedTargetLevel_U<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::PathBasedTargetLevel_U<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::PathBasedTargetLevel_U<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::BackTracking_U<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::BackTracking_U<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::BackTracking_U<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::BackTracking_U<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::CubicInterp_U<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::CubicInterp_U<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::CubicInterp_U<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::CubicInterp_U<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::Bracketing<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::Bracketing<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::Bracketing<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::Bracketing<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::BrentsScalarMinimization<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::BrentsScalarMinimization<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::BrentsScalarMinimization<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::BrentsScalarMinimization<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::BisectionScalarMinimization<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::BisectionScalarMinimization<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::BisectionScalarMinimization<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::BisectionScalarMinimization<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::GoldenSectionScalarMinimization<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::GoldenSectionScalarMinimization<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::GoldenSectionScalarMinimization<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::GoldenSectionScalarMinimization<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::ScalarMinimizationLineSearch_U<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::ScalarMinimizationLineSearch_U<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::ScalarMinimizationLineSearch_U<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::ScalarMinimizationLineSearch_U<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeU::BundleAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::TypeU::BundleAlgorithm<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::TypeU::BundleAlgorithm<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeU::BundleAlgorithm<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::Gradient_U<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::Gradient_U<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::Gradient_U<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::Gradient_U<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(struct ROL::NonlinearCGState<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(struct ROL::NonlinearCGState<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::NonlinearCGState<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(struct ROL::NonlinearCGState<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::NonlinearCG<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::NonlinearCG<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::NonlinearCG<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::NonlinearCG<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::NonlinearCG_U<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::NonlinearCG_U<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::NonlinearCG_U<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::NonlinearCG_U<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::QuasiNewton_U<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::QuasiNewton_U<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::QuasiNewton_U<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::QuasiNewton_U<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::Newton_U<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::Newton_U<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::Newton_U<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::Newton_U<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::ConjugateResiduals<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::ConjugateResiduals<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::ConjugateResiduals<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::ConjugateResiduals<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::ConjugateGradients<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::ConjugateGradients<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::ConjugateGradients<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::ConjugateGradients<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class Teuchos::SerialDenseMatrix<int, double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class Teuchos::SerialDenseMatrix<int, double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<Teuchos::SerialDenseMatrix<int, double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class Teuchos::SerialDenseMatrix<int, double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class Teuchos::SerialDenseVector<int, double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class Teuchos::SerialDenseVector<int, double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<Teuchos::SerialDenseVector<int, double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class Teuchos::SerialDenseVector<int, double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::GMRES<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::GMRES<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::GMRES<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::GMRES<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::NewtonKrylov_U<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::NewtonKrylov_U<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::NewtonKrylov_U<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::NewtonKrylov_U<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeU::LineSearchAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::TypeU::LineSearchAlgorithm<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::TypeU::LineSearchAlgorithm<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeU::LineSearchAlgorithm<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(struct ROL::TypeB::AlgorithmState<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(struct ROL::TypeB::AlgorithmState<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::TypeB::AlgorithmState<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(struct ROL::TypeB::AlgorithmState<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::PolyhedralProjection<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::PolyhedralProjection<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::PolyhedralProjection<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::PolyhedralProjection<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeB::NewtonKrylovAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::TypeB::NewtonKrylovAlgorithm<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::TypeB::NewtonKrylovAlgorithm<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeB::NewtonKrylovAlgorithm<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::ReducedLinearConstraint<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::ReducedLinearConstraint<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::ReducedLinearConstraint<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::ReducedLinearConstraint<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeB::LSecantBAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::TypeB::LSecantBAlgorithm<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::TypeB::LSecantBAlgorithm<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeB::LSecantBAlgorithm<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::PQNObjective<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::PQNObjective<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::PQNObjective<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::PQNObjective<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::Problem<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::Problem<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::Problem<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::Problem<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeB::QuasiNewtonAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::TypeB::QuasiNewtonAlgorithm<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::TypeB::QuasiNewtonAlgorithm<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeB::QuasiNewtonAlgorithm<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeB::GradientAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::TypeB::GradientAlgorithm<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::TypeB::GradientAlgorithm<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeB::GradientAlgorithm<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeB::KelleySachsAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::TypeB::KelleySachsAlgorithm<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::TypeB::KelleySachsAlgorithm<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeB::KelleySachsAlgorithm<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeB::TrustRegionSPGAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::TypeB::TrustRegionSPGAlgorithm<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::TypeB::TrustRegionSPGAlgorithm<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeB::TrustRegionSPGAlgorithm<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeB::ColemanLiAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::TypeB::ColemanLiAlgorithm<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::TypeB::ColemanLiAlgorithm<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeB::ColemanLiAlgorithm<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeB::LinMoreAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::TypeB::LinMoreAlgorithm<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::TypeB::LinMoreAlgorithm<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeB::LinMoreAlgorithm<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeB::MoreauYosidaAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::TypeB::MoreauYosidaAlgorithm<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::TypeB::MoreauYosidaAlgorithm<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeB::MoreauYosidaAlgorithm<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::PartitionedVector<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::PartitionedVector<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::PartitionedVector<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::PartitionedVector<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeB::PrimalDualActiveSetAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::TypeB::PrimalDualActiveSetAlgorithm<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::TypeB::PrimalDualActiveSetAlgorithm<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeB::PrimalDualActiveSetAlgorithm<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeB::InteriorPointAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::TypeB::InteriorPointAlgorithm<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::TypeB::InteriorPointAlgorithm<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeB::InteriorPointAlgorithm<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeB::SpectralGradientAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::TypeB::SpectralGradientAlgorithm<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::TypeB::SpectralGradientAlgorithm<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeB::SpectralGradientAlgorithm<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeE::FletcherAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::TypeE::FletcherAlgorithm<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::TypeE::FletcherAlgorithm<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeE::FletcherAlgorithm<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeE::CompositeStepAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::TypeE::CompositeStepAlgorithm<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::TypeE::CompositeStepAlgorithm<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeE::CompositeStepAlgorithm<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::AugmentedLagrangianObjective<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::AugmentedLagrangianObjective<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::AugmentedLagrangianObjective<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::AugmentedLagrangianObjective<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::ElasticLinearConstraint<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::ElasticLinearConstraint<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::ElasticLinearConstraint<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::ElasticLinearConstraint<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::BoundConstraint<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::BoundConstraint<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::BoundConstraint<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::BoundConstraint<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::Bounds<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::Bounds<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::Bounds<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::Bounds<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::BoundConstraint_Partitioned<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::BoundConstraint_Partitioned<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::BoundConstraint_Partitioned<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::BoundConstraint_Partitioned<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeE::StabilizedLCLAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::TypeE::StabilizedLCLAlgorithm<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::TypeE::StabilizedLCLAlgorithm<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeE::StabilizedLCLAlgorithm<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(struct ROL::TypeG::AlgorithmState<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(struct ROL::TypeG::AlgorithmState<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::TypeG::AlgorithmState<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(struct ROL::TypeG::AlgorithmState<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeG::AugmentedLagrangianAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::TypeG::AugmentedLagrangianAlgorithm<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::TypeG::AugmentedLagrangianAlgorithm<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeG::AugmentedLagrangianAlgorithm<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeG::MoreauYosidaAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::TypeG::MoreauYosidaAlgorithm<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::TypeG::MoreauYosidaAlgorithm<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeG::MoreauYosidaAlgorithm<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeG::InteriorPointAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::TypeG::InteriorPointAlgorithm<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::TypeG::InteriorPointAlgorithm<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeG::InteriorPointAlgorithm<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeG::StabilizedLCLAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::TypeG::StabilizedLCLAlgorithm<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::TypeG::StabilizedLCLAlgorithm<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TypeG::StabilizedLCLAlgorithm<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(struct ROL::AlgorithmState<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(struct ROL::AlgorithmState<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::AlgorithmState<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(struct ROL::AlgorithmState<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::CauchyPoint<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::CauchyPoint<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::CauchyPoint<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::CauchyPoint<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::DogLeg<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::DogLeg<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::DogLeg<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::DogLeg<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::DoubleDogLeg<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::DoubleDogLeg<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::DoubleDogLeg<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::DoubleDogLeg<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TruncatedCG<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::TruncatedCG<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::TruncatedCG<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TruncatedCG<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::LinMore<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::LinMore<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::LinMore<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::LinMore<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(struct ROL::StepState<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(struct ROL::StepState<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::StepState<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(struct ROL::StepState<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(struct ROL::SecantState<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(struct ROL::SecantState<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::SecantState<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(struct ROL::SecantState<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::lBFGS<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::lBFGS<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::lBFGS<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::lBFGS<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::lDFP<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::lDFP<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::lDFP<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::lDFP<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::lSR1<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::lSR1<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::lSR1<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::lSR1<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::BarzilaiBorwein<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::BarzilaiBorwein<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::BarzilaiBorwein<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::BarzilaiBorwein<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::KelleySachsModel<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::KelleySachsModel<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::KelleySachsModel<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::KelleySachsModel<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::ColemanLiModel<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::ColemanLiModel<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::ColemanLiModel<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::ColemanLiModel<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::LinMoreModel<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::LinMoreModel<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::LinMoreModel<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::LinMoreModel<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TrustRegionModel<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::TrustRegionModel<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::TrustRegionModel<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::TrustRegionModel<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::Constraint_Partitioned<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::Constraint_Partitioned<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::Constraint_Partitioned<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::Constraint_Partitioned<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::SlacklessObjective<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::SlacklessObjective<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::SlacklessObjective<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::SlacklessObjective<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::ReduceLinearConstraint<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::ReduceLinearConstraint<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::ReduceLinearConstraint<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::ReduceLinearConstraint<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::LinearConstraint<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::LinearConstraint<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::LinearConstraint<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::LinearConstraint<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::AffineTransformObjective<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::AffineTransformObjective<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::AffineTransformObjective<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::AffineTransformObjective<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::AffineTransformConstraint<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::AffineTransformConstraint<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::AffineTransformConstraint<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::AffineTransformConstraint<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::DaiFletcherProjection<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::DaiFletcherProjection<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::DaiFletcherProjection<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::DaiFletcherProjection<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::DykstraProjection<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::DykstraProjection<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::DykstraProjection<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::DykstraProjection<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::DouglasRachfordProjection<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::DouglasRachfordProjection<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::DouglasRachfordProjection<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::DouglasRachfordProjection<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::SemismoothNewtonProjection<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::SemismoothNewtonProjection<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::SemismoothNewtonProjection<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::SemismoothNewtonProjection<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::RiddersProjection<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::RiddersProjection<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::RiddersProjection<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::RiddersProjection<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

	// Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::BrentsProjection<double> *, bool) file:Teuchos_RCP.hpp line:91
	M("Teuchos").def("RCP_createNewRCPNodeRawPtr", (class Teuchos::RCPNode * (*)(class ROL::BrentsProjection<double> *, bool)) &Teuchos::RCP_createNewRCPNodeRawPtr<ROL::BrentsProjection<double>>, "C++: Teuchos::RCP_createNewRCPNodeRawPtr(class ROL::BrentsProjection<double> *, bool) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("p"), pybind11::arg("has_ownership_in"));

}
