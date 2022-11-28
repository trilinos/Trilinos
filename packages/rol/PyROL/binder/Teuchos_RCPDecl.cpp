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
#include <Teuchos_PtrDecl.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_RCPDecl.hpp>
#include <Teuchos_RCPNode.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_StringIndexedOrderedValueObjectContainer.hpp>
#include <Teuchos_XMLObject.hpp>
#include <Teuchos_XMLObjectImplem.hpp>
#include <Teuchos_any.hpp>
#include <Teuchos_basic_oblackholestream.hpp>
#include <Teuchos_dyn_cast.hpp>
#include <cwchar>
#include <deque>
#include <functional>
#include <ios>
#include <iterator>
#include <locale>
#include <map>
#include <memory>
#include <ostream>
#include <sstream>
#include <sstream> // __str__
#include <streambuf>
#include <string>
#include <typeinfo>
#include <utility>
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

// Teuchos::m_bad_cast file:Teuchos_dyn_cast.hpp line:60
struct PyCallBack_Teuchos_m_bad_cast : public Teuchos::m_bad_cast {
	using Teuchos::m_bad_cast::m_bad_cast;

	const char * what() const throw() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::m_bad_cast *>(this), "what");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const char *>::value) {
				static pybind11::detail::override_caster_t<const char *> caster;
				return pybind11::detail::cast_ref<const char *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const char *>(std::move(o));
		}
		return m_bad_cast::what();
	}
};

void bind_Teuchos_RCPDecl(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<ROL::VectorController<double, int>>> cl(M("Teuchos"), "RCP_ROL_VectorController_double_int_t", "");
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<ROL::VectorController<double, int>>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("") );

		cl.def( pybind11::init( [](class ROL::VectorController<double, int> * a0){ return new Teuchos::RCP<ROL::VectorController<double, int>>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<class ROL::VectorController<double, int> *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership_in") );

		cl.def( pybind11::init( [](Teuchos::RCP<ROL::VectorController<double, int>> const &o){ return new Teuchos::RCP<ROL::VectorController<double, int>>(o); } ) );
		cl.def( pybind11::init<class ROL::VectorController<double, int> *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::VectorController<double, int> *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::VectorController<double, int> *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<class ROL::VectorController<double, int> > & (Teuchos::RCP<ROL::VectorController<double, int>>::*)(const class Teuchos::RCP<class ROL::VectorController<double, int> > &)) &Teuchos::RCP<ROL::VectorController<double, int> >::operator=, "C++: Teuchos::RCP<ROL::VectorController<double, int> >::operator=(const class Teuchos::RCP<class ROL::VectorController<double, int> > &) --> class Teuchos::RCP<class ROL::VectorController<double, int> > &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<class ROL::VectorController<double, int> > & (Teuchos::RCP<ROL::VectorController<double, int>>::*)(enum Teuchos::ENull)) &Teuchos::RCP<ROL::VectorController<double, int> >::operator=, "C++: Teuchos::RCP<ROL::VectorController<double, int> >::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<class ROL::VectorController<double, int> > &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<ROL::VectorController<double, int>>::*)(class Teuchos::RCP<class ROL::VectorController<double, int> > &)) &Teuchos::RCP<ROL::VectorController<double, int> >::swap, "C++: Teuchos::RCP<ROL::VectorController<double, int> >::swap(class Teuchos::RCP<class ROL::VectorController<double, int> > &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<ROL::VectorController<double, int>>::*)() const) &Teuchos::RCP<ROL::VectorController<double, int> >::is_null, "C++: Teuchos::RCP<ROL::VectorController<double, int> >::is_null() const --> bool");
		cl.def("arrow", (class ROL::VectorController<double, int> * (Teuchos::RCP<ROL::VectorController<double, int>>::*)() const) &Teuchos::RCP<ROL::VectorController<double, int> >::operator->, "C++: Teuchos::RCP<ROL::VectorController<double, int> >::operator->() const --> class ROL::VectorController<double, int> *", pybind11::return_value_policy::automatic);
		cl.def("__mul__", (class ROL::VectorController<double, int> & (Teuchos::RCP<ROL::VectorController<double, int>>::*)() const) &Teuchos::RCP<ROL::VectorController<double, int> >::operator*, "C++: Teuchos::RCP<ROL::VectorController<double, int> >::operator*() const --> class ROL::VectorController<double, int> &", pybind11::return_value_policy::automatic);
		cl.def("get", (class ROL::VectorController<double, int> * (Teuchos::RCP<ROL::VectorController<double, int>>::*)() const) &Teuchos::RCP<ROL::VectorController<double, int> >::get, "C++: Teuchos::RCP<ROL::VectorController<double, int> >::get() const --> class ROL::VectorController<double, int> *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class ROL::VectorController<double, int> * (Teuchos::RCP<ROL::VectorController<double, int>>::*)() const) &Teuchos::RCP<ROL::VectorController<double, int> >::getRawPtr, "C++: Teuchos::RCP<ROL::VectorController<double, int> >::getRawPtr() const --> class ROL::VectorController<double, int> *", pybind11::return_value_policy::automatic);
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<ROL::VectorController<double, int>>::*)() const) &Teuchos::RCP<ROL::VectorController<double, int> >::strength, "C++: Teuchos::RCP<ROL::VectorController<double, int> >::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<ROL::VectorController<double, int>>::*)() const) &Teuchos::RCP<ROL::VectorController<double, int> >::is_valid_ptr, "C++: Teuchos::RCP<ROL::VectorController<double, int> >::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<ROL::VectorController<double, int>>::*)() const) &Teuchos::RCP<ROL::VectorController<double, int> >::strong_count, "C++: Teuchos::RCP<ROL::VectorController<double, int> >::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<ROL::VectorController<double, int>>::*)() const) &Teuchos::RCP<ROL::VectorController<double, int> >::weak_count, "C++: Teuchos::RCP<ROL::VectorController<double, int> >::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<ROL::VectorController<double, int>>::*)() const) &Teuchos::RCP<ROL::VectorController<double, int> >::total_count, "C++: Teuchos::RCP<ROL::VectorController<double, int> >::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<ROL::VectorController<double, int>>::*)()) &Teuchos::RCP<ROL::VectorController<double, int> >::set_has_ownership, "C++: Teuchos::RCP<ROL::VectorController<double, int> >::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<ROL::VectorController<double, int>>::*)() const) &Teuchos::RCP<ROL::VectorController<double, int> >::has_ownership, "C++: Teuchos::RCP<ROL::VectorController<double, int> >::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<class ROL::VectorController<double, int> > (Teuchos::RCP<ROL::VectorController<double, int>>::*)() const) &Teuchos::RCP<ROL::VectorController<double, int> >::create_weak, "C++: Teuchos::RCP<ROL::VectorController<double, int> >::create_weak() const --> class Teuchos::RCP<class ROL::VectorController<double, int> >");
		cl.def("create_strong", (class Teuchos::RCP<class ROL::VectorController<double, int> > (Teuchos::RCP<ROL::VectorController<double, int>>::*)() const) &Teuchos::RCP<ROL::VectorController<double, int> >::create_strong, "C++: Teuchos::RCP<ROL::VectorController<double, int> >::create_strong() const --> class Teuchos::RCP<class ROL::VectorController<double, int> >");
		cl.def("assert_not_null", (const class Teuchos::RCP<class ROL::VectorController<double, int> > & (Teuchos::RCP<ROL::VectorController<double, int>>::*)() const) &Teuchos::RCP<ROL::VectorController<double, int> >::assert_not_null, "C++: Teuchos::RCP<ROL::VectorController<double, int> >::assert_not_null() const --> const class Teuchos::RCP<class ROL::VectorController<double, int> > &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<class ROL::VectorController<double, int> > & (Teuchos::RCP<ROL::VectorController<double, int>>::*)() const) &Teuchos::RCP<ROL::VectorController<double, int> >::assert_valid_ptr, "C++: Teuchos::RCP<ROL::VectorController<double, int> >::assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::VectorController<double, int> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<class ROL::VectorController<double, int> > & (Teuchos::RCP<ROL::VectorController<double, int>>::*)() const) &Teuchos::RCP<ROL::VectorController<double, int> >::debug_assert_not_null, "C++: Teuchos::RCP<ROL::VectorController<double, int> >::debug_assert_not_null() const --> const class Teuchos::RCP<class ROL::VectorController<double, int> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<class ROL::VectorController<double, int> > & (Teuchos::RCP<ROL::VectorController<double, int>>::*)() const) &Teuchos::RCP<ROL::VectorController<double, int> >::debug_assert_valid_ptr, "C++: Teuchos::RCP<ROL::VectorController<double, int> >::debug_assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::VectorController<double, int> > &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<ROL::VectorController<double, int>>::*)()) &Teuchos::RCP<ROL::VectorController<double, int> >::reset, "C++: Teuchos::RCP<ROL::VectorController<double, int> >::reset() --> void");
		cl.def("access_private_ptr", (class ROL::VectorController<double, int> * (Teuchos::RCP<ROL::VectorController<double, int>>::*)() const) &Teuchos::RCP<ROL::VectorController<double, int> >::access_private_ptr, "C++: Teuchos::RCP<ROL::VectorController<double, int> >::access_private_ptr() const --> class ROL::VectorController<double, int> *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::VectorController<double, int>>::*)()) &Teuchos::RCP<ROL::VectorController<double, int> >::nonconst_access_private_node, "C++: Teuchos::RCP<ROL::VectorController<double, int> >::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::VectorController<double, int>>::*)() const) &Teuchos::RCP<ROL::VectorController<double, int> >::access_private_node, "C++: Teuchos::RCP<ROL::VectorController<double, int> >::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<ROL::VectorController<double, int>> const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<ROL::SingletonVector<double>>> cl(M("Teuchos"), "RCP_ROL_SingletonVector_double_t", "");
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<ROL::SingletonVector<double>>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("null_arg") );

		cl.def( pybind11::init( [](class ROL::SingletonVector<double> * a0){ return new Teuchos::RCP<ROL::SingletonVector<double>>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<class ROL::SingletonVector<double> *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership_in") );

		cl.def( pybind11::init( [](Teuchos::RCP<ROL::SingletonVector<double>> const &o){ return new Teuchos::RCP<ROL::SingletonVector<double>>(o); } ) );
		cl.def( pybind11::init<class ROL::SingletonVector<double> *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::SingletonVector<double> *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::SingletonVector<double> *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<class ROL::SingletonVector<double> > & (Teuchos::RCP<ROL::SingletonVector<double>>::*)(const class Teuchos::RCP<class ROL::SingletonVector<double> > &)) &Teuchos::RCP<ROL::SingletonVector<double> >::operator=, "C++: Teuchos::RCP<ROL::SingletonVector<double> >::operator=(const class Teuchos::RCP<class ROL::SingletonVector<double> > &) --> class Teuchos::RCP<class ROL::SingletonVector<double> > &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<class ROL::SingletonVector<double> > & (Teuchos::RCP<ROL::SingletonVector<double>>::*)(enum Teuchos::ENull)) &Teuchos::RCP<ROL::SingletonVector<double> >::operator=, "C++: Teuchos::RCP<ROL::SingletonVector<double> >::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<class ROL::SingletonVector<double> > &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<ROL::SingletonVector<double>>::*)(class Teuchos::RCP<class ROL::SingletonVector<double> > &)) &Teuchos::RCP<ROL::SingletonVector<double> >::swap, "C++: Teuchos::RCP<ROL::SingletonVector<double> >::swap(class Teuchos::RCP<class ROL::SingletonVector<double> > &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<ROL::SingletonVector<double>>::*)() const) &Teuchos::RCP<ROL::SingletonVector<double> >::is_null, "C++: Teuchos::RCP<ROL::SingletonVector<double> >::is_null() const --> bool");
		cl.def("arrow", (class ROL::SingletonVector<double> * (Teuchos::RCP<ROL::SingletonVector<double>>::*)() const) &Teuchos::RCP<ROL::SingletonVector<double> >::operator->, "C++: Teuchos::RCP<ROL::SingletonVector<double> >::operator->() const --> class ROL::SingletonVector<double> *", pybind11::return_value_policy::automatic);
		cl.def("__mul__", (class ROL::SingletonVector<double> & (Teuchos::RCP<ROL::SingletonVector<double>>::*)() const) &Teuchos::RCP<ROL::SingletonVector<double> >::operator*, "C++: Teuchos::RCP<ROL::SingletonVector<double> >::operator*() const --> class ROL::SingletonVector<double> &", pybind11::return_value_policy::automatic);
		cl.def("get", (class ROL::SingletonVector<double> * (Teuchos::RCP<ROL::SingletonVector<double>>::*)() const) &Teuchos::RCP<ROL::SingletonVector<double> >::get, "C++: Teuchos::RCP<ROL::SingletonVector<double> >::get() const --> class ROL::SingletonVector<double> *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class ROL::SingletonVector<double> * (Teuchos::RCP<ROL::SingletonVector<double>>::*)() const) &Teuchos::RCP<ROL::SingletonVector<double> >::getRawPtr, "C++: Teuchos::RCP<ROL::SingletonVector<double> >::getRawPtr() const --> class ROL::SingletonVector<double> *", pybind11::return_value_policy::automatic);
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<ROL::SingletonVector<double>>::*)() const) &Teuchos::RCP<ROL::SingletonVector<double> >::strength, "C++: Teuchos::RCP<ROL::SingletonVector<double> >::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<ROL::SingletonVector<double>>::*)() const) &Teuchos::RCP<ROL::SingletonVector<double> >::is_valid_ptr, "C++: Teuchos::RCP<ROL::SingletonVector<double> >::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<ROL::SingletonVector<double>>::*)() const) &Teuchos::RCP<ROL::SingletonVector<double> >::strong_count, "C++: Teuchos::RCP<ROL::SingletonVector<double> >::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<ROL::SingletonVector<double>>::*)() const) &Teuchos::RCP<ROL::SingletonVector<double> >::weak_count, "C++: Teuchos::RCP<ROL::SingletonVector<double> >::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<ROL::SingletonVector<double>>::*)() const) &Teuchos::RCP<ROL::SingletonVector<double> >::total_count, "C++: Teuchos::RCP<ROL::SingletonVector<double> >::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<ROL::SingletonVector<double>>::*)()) &Teuchos::RCP<ROL::SingletonVector<double> >::set_has_ownership, "C++: Teuchos::RCP<ROL::SingletonVector<double> >::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<ROL::SingletonVector<double>>::*)() const) &Teuchos::RCP<ROL::SingletonVector<double> >::has_ownership, "C++: Teuchos::RCP<ROL::SingletonVector<double> >::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<class ROL::SingletonVector<double> > (Teuchos::RCP<ROL::SingletonVector<double>>::*)() const) &Teuchos::RCP<ROL::SingletonVector<double> >::create_weak, "C++: Teuchos::RCP<ROL::SingletonVector<double> >::create_weak() const --> class Teuchos::RCP<class ROL::SingletonVector<double> >");
		cl.def("create_strong", (class Teuchos::RCP<class ROL::SingletonVector<double> > (Teuchos::RCP<ROL::SingletonVector<double>>::*)() const) &Teuchos::RCP<ROL::SingletonVector<double> >::create_strong, "C++: Teuchos::RCP<ROL::SingletonVector<double> >::create_strong() const --> class Teuchos::RCP<class ROL::SingletonVector<double> >");
		cl.def("assert_not_null", (const class Teuchos::RCP<class ROL::SingletonVector<double> > & (Teuchos::RCP<ROL::SingletonVector<double>>::*)() const) &Teuchos::RCP<ROL::SingletonVector<double> >::assert_not_null, "C++: Teuchos::RCP<ROL::SingletonVector<double> >::assert_not_null() const --> const class Teuchos::RCP<class ROL::SingletonVector<double> > &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<class ROL::SingletonVector<double> > & (Teuchos::RCP<ROL::SingletonVector<double>>::*)() const) &Teuchos::RCP<ROL::SingletonVector<double> >::assert_valid_ptr, "C++: Teuchos::RCP<ROL::SingletonVector<double> >::assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::SingletonVector<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<class ROL::SingletonVector<double> > & (Teuchos::RCP<ROL::SingletonVector<double>>::*)() const) &Teuchos::RCP<ROL::SingletonVector<double> >::debug_assert_not_null, "C++: Teuchos::RCP<ROL::SingletonVector<double> >::debug_assert_not_null() const --> const class Teuchos::RCP<class ROL::SingletonVector<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<class ROL::SingletonVector<double> > & (Teuchos::RCP<ROL::SingletonVector<double>>::*)() const) &Teuchos::RCP<ROL::SingletonVector<double> >::debug_assert_valid_ptr, "C++: Teuchos::RCP<ROL::SingletonVector<double> >::debug_assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::SingletonVector<double> > &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<ROL::SingletonVector<double>>::*)()) &Teuchos::RCP<ROL::SingletonVector<double> >::reset, "C++: Teuchos::RCP<ROL::SingletonVector<double> >::reset() --> void");
		cl.def("access_private_ptr", (class ROL::SingletonVector<double> * (Teuchos::RCP<ROL::SingletonVector<double>>::*)() const) &Teuchos::RCP<ROL::SingletonVector<double> >::access_private_ptr, "C++: Teuchos::RCP<ROL::SingletonVector<double> >::access_private_ptr() const --> class ROL::SingletonVector<double> *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::SingletonVector<double>>::*)()) &Teuchos::RCP<ROL::SingletonVector<double> >::nonconst_access_private_node, "C++: Teuchos::RCP<ROL::SingletonVector<double> >::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::SingletonVector<double>>::*)() const) &Teuchos::RCP<ROL::SingletonVector<double> >::access_private_node, "C++: Teuchos::RCP<ROL::SingletonVector<double> >::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<ROL::SingletonVector<double>> const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double>>> cl(M("Teuchos"), "RCP_ROL_TypeE_AugmentedLagrangianAlgorithm_double_t", "");
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double>>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("null_arg") );

		cl.def( pybind11::init( [](class ROL::TypeE::AugmentedLagrangianAlgorithm<double> * a0){ return new Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double>>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<class ROL::TypeE::AugmentedLagrangianAlgorithm<double> *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership_in") );

		cl.def( pybind11::init( [](Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double>> const &o){ return new Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double>>(o); } ) );
		cl.def( pybind11::init<class ROL::TypeE::AugmentedLagrangianAlgorithm<double> *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::TypeE::AugmentedLagrangianAlgorithm<double> *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::TypeE::AugmentedLagrangianAlgorithm<double> *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<class ROL::TypeE::AugmentedLagrangianAlgorithm<double> > & (Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double>>::*)(const class Teuchos::RCP<class ROL::TypeE::AugmentedLagrangianAlgorithm<double> > &)) &Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::operator=, "C++: Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::operator=(const class Teuchos::RCP<class ROL::TypeE::AugmentedLagrangianAlgorithm<double> > &) --> class Teuchos::RCP<class ROL::TypeE::AugmentedLagrangianAlgorithm<double> > &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<class ROL::TypeE::AugmentedLagrangianAlgorithm<double> > & (Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double>>::*)(enum Teuchos::ENull)) &Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::operator=, "C++: Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<class ROL::TypeE::AugmentedLagrangianAlgorithm<double> > &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double>>::*)(class Teuchos::RCP<class ROL::TypeE::AugmentedLagrangianAlgorithm<double> > &)) &Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::swap, "C++: Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::swap(class Teuchos::RCP<class ROL::TypeE::AugmentedLagrangianAlgorithm<double> > &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::is_null, "C++: Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::is_null() const --> bool");
		cl.def("arrow", (class ROL::TypeE::AugmentedLagrangianAlgorithm<double> * (Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::operator->, "C++: Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::operator->() const --> class ROL::TypeE::AugmentedLagrangianAlgorithm<double> *", pybind11::return_value_policy::automatic);
		cl.def("__mul__", (class ROL::TypeE::AugmentedLagrangianAlgorithm<double> & (Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::operator*, "C++: Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::operator*() const --> class ROL::TypeE::AugmentedLagrangianAlgorithm<double> &", pybind11::return_value_policy::automatic);
		cl.def("get", (class ROL::TypeE::AugmentedLagrangianAlgorithm<double> * (Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::get, "C++: Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::get() const --> class ROL::TypeE::AugmentedLagrangianAlgorithm<double> *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class ROL::TypeE::AugmentedLagrangianAlgorithm<double> * (Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::getRawPtr, "C++: Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::getRawPtr() const --> class ROL::TypeE::AugmentedLagrangianAlgorithm<double> *", pybind11::return_value_policy::automatic);
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::strength, "C++: Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::is_valid_ptr, "C++: Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::strong_count, "C++: Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::weak_count, "C++: Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::total_count, "C++: Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double>>::*)()) &Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::set_has_ownership, "C++: Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::has_ownership, "C++: Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<class ROL::TypeE::AugmentedLagrangianAlgorithm<double> > (Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::create_weak, "C++: Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::create_weak() const --> class Teuchos::RCP<class ROL::TypeE::AugmentedLagrangianAlgorithm<double> >");
		cl.def("create_strong", (class Teuchos::RCP<class ROL::TypeE::AugmentedLagrangianAlgorithm<double> > (Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::create_strong, "C++: Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::create_strong() const --> class Teuchos::RCP<class ROL::TypeE::AugmentedLagrangianAlgorithm<double> >");
		cl.def("assert_not_null", (const class Teuchos::RCP<class ROL::TypeE::AugmentedLagrangianAlgorithm<double> > & (Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::assert_not_null, "C++: Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::assert_not_null() const --> const class Teuchos::RCP<class ROL::TypeE::AugmentedLagrangianAlgorithm<double> > &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<class ROL::TypeE::AugmentedLagrangianAlgorithm<double> > & (Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::assert_valid_ptr, "C++: Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::TypeE::AugmentedLagrangianAlgorithm<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<class ROL::TypeE::AugmentedLagrangianAlgorithm<double> > & (Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::debug_assert_not_null, "C++: Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::debug_assert_not_null() const --> const class Teuchos::RCP<class ROL::TypeE::AugmentedLagrangianAlgorithm<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<class ROL::TypeE::AugmentedLagrangianAlgorithm<double> > & (Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::debug_assert_valid_ptr, "C++: Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::debug_assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::TypeE::AugmentedLagrangianAlgorithm<double> > &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double>>::*)()) &Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::reset, "C++: Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::reset() --> void");
		cl.def("access_private_ptr", (class ROL::TypeE::AugmentedLagrangianAlgorithm<double> * (Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::access_private_ptr, "C++: Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::access_private_ptr() const --> class ROL::TypeE::AugmentedLagrangianAlgorithm<double> *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double>>::*)()) &Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::nonconst_access_private_node, "C++: Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::access_private_node, "C++: Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double> >::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<ROL::TypeE::AugmentedLagrangianAlgorithm<double>> const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double>>> cl(M("Teuchos"), "RCP_ROL_TypeE_FletcherAlgorithm_double_t", "");
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double>>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("null_arg") );

		cl.def( pybind11::init( [](class ROL::TypeE::FletcherAlgorithm<double> * a0){ return new Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double>>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<class ROL::TypeE::FletcherAlgorithm<double> *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership_in") );

		cl.def( pybind11::init( [](Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double>> const &o){ return new Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double>>(o); } ) );
		cl.def( pybind11::init<class ROL::TypeE::FletcherAlgorithm<double> *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::TypeE::FletcherAlgorithm<double> *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::TypeE::FletcherAlgorithm<double> *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<class ROL::TypeE::FletcherAlgorithm<double> > & (Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double>>::*)(const class Teuchos::RCP<class ROL::TypeE::FletcherAlgorithm<double> > &)) &Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::operator=, "C++: Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::operator=(const class Teuchos::RCP<class ROL::TypeE::FletcherAlgorithm<double> > &) --> class Teuchos::RCP<class ROL::TypeE::FletcherAlgorithm<double> > &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<class ROL::TypeE::FletcherAlgorithm<double> > & (Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double>>::*)(enum Teuchos::ENull)) &Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::operator=, "C++: Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<class ROL::TypeE::FletcherAlgorithm<double> > &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double>>::*)(class Teuchos::RCP<class ROL::TypeE::FletcherAlgorithm<double> > &)) &Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::swap, "C++: Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::swap(class Teuchos::RCP<class ROL::TypeE::FletcherAlgorithm<double> > &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::is_null, "C++: Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::is_null() const --> bool");
		cl.def("arrow", (class ROL::TypeE::FletcherAlgorithm<double> * (Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::operator->, "C++: Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::operator->() const --> class ROL::TypeE::FletcherAlgorithm<double> *", pybind11::return_value_policy::automatic);
		cl.def("__mul__", (class ROL::TypeE::FletcherAlgorithm<double> & (Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::operator*, "C++: Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::operator*() const --> class ROL::TypeE::FletcherAlgorithm<double> &", pybind11::return_value_policy::automatic);
		cl.def("get", (class ROL::TypeE::FletcherAlgorithm<double> * (Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::get, "C++: Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::get() const --> class ROL::TypeE::FletcherAlgorithm<double> *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class ROL::TypeE::FletcherAlgorithm<double> * (Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::getRawPtr, "C++: Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::getRawPtr() const --> class ROL::TypeE::FletcherAlgorithm<double> *", pybind11::return_value_policy::automatic);
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::strength, "C++: Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::is_valid_ptr, "C++: Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::strong_count, "C++: Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::weak_count, "C++: Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::total_count, "C++: Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double>>::*)()) &Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::set_has_ownership, "C++: Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::has_ownership, "C++: Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<class ROL::TypeE::FletcherAlgorithm<double> > (Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::create_weak, "C++: Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::create_weak() const --> class Teuchos::RCP<class ROL::TypeE::FletcherAlgorithm<double> >");
		cl.def("create_strong", (class Teuchos::RCP<class ROL::TypeE::FletcherAlgorithm<double> > (Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::create_strong, "C++: Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::create_strong() const --> class Teuchos::RCP<class ROL::TypeE::FletcherAlgorithm<double> >");
		cl.def("assert_not_null", (const class Teuchos::RCP<class ROL::TypeE::FletcherAlgorithm<double> > & (Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::assert_not_null, "C++: Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::assert_not_null() const --> const class Teuchos::RCP<class ROL::TypeE::FletcherAlgorithm<double> > &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<class ROL::TypeE::FletcherAlgorithm<double> > & (Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::assert_valid_ptr, "C++: Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::TypeE::FletcherAlgorithm<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<class ROL::TypeE::FletcherAlgorithm<double> > & (Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::debug_assert_not_null, "C++: Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::debug_assert_not_null() const --> const class Teuchos::RCP<class ROL::TypeE::FletcherAlgorithm<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<class ROL::TypeE::FletcherAlgorithm<double> > & (Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::debug_assert_valid_ptr, "C++: Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::debug_assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::TypeE::FletcherAlgorithm<double> > &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double>>::*)()) &Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::reset, "C++: Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::reset() --> void");
		cl.def("access_private_ptr", (class ROL::TypeE::FletcherAlgorithm<double> * (Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::access_private_ptr, "C++: Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::access_private_ptr() const --> class ROL::TypeE::FletcherAlgorithm<double> *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double>>::*)()) &Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::nonconst_access_private_node, "C++: Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::access_private_node, "C++: Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double> >::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<ROL::TypeE::FletcherAlgorithm<double>> const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double>>> cl(M("Teuchos"), "RCP_ROL_TypeE_CompositeStepAlgorithm_double_t", "");
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double>>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("null_arg") );

		cl.def( pybind11::init( [](class ROL::TypeE::CompositeStepAlgorithm<double> * a0){ return new Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double>>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<class ROL::TypeE::CompositeStepAlgorithm<double> *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership_in") );

		cl.def( pybind11::init( [](Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double>> const &o){ return new Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double>>(o); } ) );
		cl.def( pybind11::init<class ROL::TypeE::CompositeStepAlgorithm<double> *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::TypeE::CompositeStepAlgorithm<double> *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::TypeE::CompositeStepAlgorithm<double> *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<class ROL::TypeE::CompositeStepAlgorithm<double> > & (Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double>>::*)(const class Teuchos::RCP<class ROL::TypeE::CompositeStepAlgorithm<double> > &)) &Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::operator=, "C++: Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::operator=(const class Teuchos::RCP<class ROL::TypeE::CompositeStepAlgorithm<double> > &) --> class Teuchos::RCP<class ROL::TypeE::CompositeStepAlgorithm<double> > &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<class ROL::TypeE::CompositeStepAlgorithm<double> > & (Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double>>::*)(enum Teuchos::ENull)) &Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::operator=, "C++: Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<class ROL::TypeE::CompositeStepAlgorithm<double> > &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double>>::*)(class Teuchos::RCP<class ROL::TypeE::CompositeStepAlgorithm<double> > &)) &Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::swap, "C++: Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::swap(class Teuchos::RCP<class ROL::TypeE::CompositeStepAlgorithm<double> > &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::is_null, "C++: Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::is_null() const --> bool");
		cl.def("arrow", (class ROL::TypeE::CompositeStepAlgorithm<double> * (Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::operator->, "C++: Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::operator->() const --> class ROL::TypeE::CompositeStepAlgorithm<double> *", pybind11::return_value_policy::automatic);
		cl.def("__mul__", (class ROL::TypeE::CompositeStepAlgorithm<double> & (Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::operator*, "C++: Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::operator*() const --> class ROL::TypeE::CompositeStepAlgorithm<double> &", pybind11::return_value_policy::automatic);
		cl.def("get", (class ROL::TypeE::CompositeStepAlgorithm<double> * (Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::get, "C++: Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::get() const --> class ROL::TypeE::CompositeStepAlgorithm<double> *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class ROL::TypeE::CompositeStepAlgorithm<double> * (Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::getRawPtr, "C++: Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::getRawPtr() const --> class ROL::TypeE::CompositeStepAlgorithm<double> *", pybind11::return_value_policy::automatic);
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::strength, "C++: Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::is_valid_ptr, "C++: Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::strong_count, "C++: Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::weak_count, "C++: Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::total_count, "C++: Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double>>::*)()) &Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::set_has_ownership, "C++: Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::has_ownership, "C++: Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<class ROL::TypeE::CompositeStepAlgorithm<double> > (Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::create_weak, "C++: Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::create_weak() const --> class Teuchos::RCP<class ROL::TypeE::CompositeStepAlgorithm<double> >");
		cl.def("create_strong", (class Teuchos::RCP<class ROL::TypeE::CompositeStepAlgorithm<double> > (Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::create_strong, "C++: Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::create_strong() const --> class Teuchos::RCP<class ROL::TypeE::CompositeStepAlgorithm<double> >");
		cl.def("assert_not_null", (const class Teuchos::RCP<class ROL::TypeE::CompositeStepAlgorithm<double> > & (Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::assert_not_null, "C++: Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::assert_not_null() const --> const class Teuchos::RCP<class ROL::TypeE::CompositeStepAlgorithm<double> > &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<class ROL::TypeE::CompositeStepAlgorithm<double> > & (Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::assert_valid_ptr, "C++: Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::TypeE::CompositeStepAlgorithm<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<class ROL::TypeE::CompositeStepAlgorithm<double> > & (Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::debug_assert_not_null, "C++: Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::debug_assert_not_null() const --> const class Teuchos::RCP<class ROL::TypeE::CompositeStepAlgorithm<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<class ROL::TypeE::CompositeStepAlgorithm<double> > & (Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::debug_assert_valid_ptr, "C++: Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::debug_assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::TypeE::CompositeStepAlgorithm<double> > &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double>>::*)()) &Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::reset, "C++: Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::reset() --> void");
		cl.def("access_private_ptr", (class ROL::TypeE::CompositeStepAlgorithm<double> * (Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::access_private_ptr, "C++: Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::access_private_ptr() const --> class ROL::TypeE::CompositeStepAlgorithm<double> *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double>>::*)()) &Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::nonconst_access_private_node, "C++: Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::access_private_node, "C++: Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double> >::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<ROL::TypeE::CompositeStepAlgorithm<double>> const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double>>> cl(M("Teuchos"), "RCP_ROL_TypeE_StabilizedLCLAlgorithm_double_t", "");
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double>>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("null_arg") );

		cl.def( pybind11::init( [](class ROL::TypeE::StabilizedLCLAlgorithm<double> * a0){ return new Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double>>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<class ROL::TypeE::StabilizedLCLAlgorithm<double> *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership_in") );

		cl.def( pybind11::init( [](Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double>> const &o){ return new Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double>>(o); } ) );
		cl.def( pybind11::init<class ROL::TypeE::StabilizedLCLAlgorithm<double> *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::TypeE::StabilizedLCLAlgorithm<double> *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::TypeE::StabilizedLCLAlgorithm<double> *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<class ROL::TypeE::StabilizedLCLAlgorithm<double> > & (Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double>>::*)(const class Teuchos::RCP<class ROL::TypeE::StabilizedLCLAlgorithm<double> > &)) &Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::operator=, "C++: Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::operator=(const class Teuchos::RCP<class ROL::TypeE::StabilizedLCLAlgorithm<double> > &) --> class Teuchos::RCP<class ROL::TypeE::StabilizedLCLAlgorithm<double> > &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<class ROL::TypeE::StabilizedLCLAlgorithm<double> > & (Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double>>::*)(enum Teuchos::ENull)) &Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::operator=, "C++: Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<class ROL::TypeE::StabilizedLCLAlgorithm<double> > &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double>>::*)(class Teuchos::RCP<class ROL::TypeE::StabilizedLCLAlgorithm<double> > &)) &Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::swap, "C++: Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::swap(class Teuchos::RCP<class ROL::TypeE::StabilizedLCLAlgorithm<double> > &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::is_null, "C++: Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::is_null() const --> bool");
		cl.def("arrow", (class ROL::TypeE::StabilizedLCLAlgorithm<double> * (Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::operator->, "C++: Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::operator->() const --> class ROL::TypeE::StabilizedLCLAlgorithm<double> *", pybind11::return_value_policy::automatic);
		cl.def("__mul__", (class ROL::TypeE::StabilizedLCLAlgorithm<double> & (Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::operator*, "C++: Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::operator*() const --> class ROL::TypeE::StabilizedLCLAlgorithm<double> &", pybind11::return_value_policy::automatic);
		cl.def("get", (class ROL::TypeE::StabilizedLCLAlgorithm<double> * (Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::get, "C++: Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::get() const --> class ROL::TypeE::StabilizedLCLAlgorithm<double> *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class ROL::TypeE::StabilizedLCLAlgorithm<double> * (Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::getRawPtr, "C++: Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::getRawPtr() const --> class ROL::TypeE::StabilizedLCLAlgorithm<double> *", pybind11::return_value_policy::automatic);
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::strength, "C++: Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::is_valid_ptr, "C++: Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::strong_count, "C++: Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::weak_count, "C++: Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::total_count, "C++: Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double>>::*)()) &Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::set_has_ownership, "C++: Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::has_ownership, "C++: Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<class ROL::TypeE::StabilizedLCLAlgorithm<double> > (Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::create_weak, "C++: Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::create_weak() const --> class Teuchos::RCP<class ROL::TypeE::StabilizedLCLAlgorithm<double> >");
		cl.def("create_strong", (class Teuchos::RCP<class ROL::TypeE::StabilizedLCLAlgorithm<double> > (Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::create_strong, "C++: Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::create_strong() const --> class Teuchos::RCP<class ROL::TypeE::StabilizedLCLAlgorithm<double> >");
		cl.def("assert_not_null", (const class Teuchos::RCP<class ROL::TypeE::StabilizedLCLAlgorithm<double> > & (Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::assert_not_null, "C++: Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::assert_not_null() const --> const class Teuchos::RCP<class ROL::TypeE::StabilizedLCLAlgorithm<double> > &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<class ROL::TypeE::StabilizedLCLAlgorithm<double> > & (Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::assert_valid_ptr, "C++: Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::TypeE::StabilizedLCLAlgorithm<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<class ROL::TypeE::StabilizedLCLAlgorithm<double> > & (Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::debug_assert_not_null, "C++: Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::debug_assert_not_null() const --> const class Teuchos::RCP<class ROL::TypeE::StabilizedLCLAlgorithm<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<class ROL::TypeE::StabilizedLCLAlgorithm<double> > & (Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::debug_assert_valid_ptr, "C++: Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::debug_assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::TypeE::StabilizedLCLAlgorithm<double> > &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double>>::*)()) &Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::reset, "C++: Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::reset() --> void");
		cl.def("access_private_ptr", (class ROL::TypeE::StabilizedLCLAlgorithm<double> * (Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::access_private_ptr, "C++: Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::access_private_ptr() const --> class ROL::TypeE::StabilizedLCLAlgorithm<double> *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double>>::*)()) &Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::nonconst_access_private_node, "C++: Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::access_private_node, "C++: Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double> >::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<ROL::TypeE::StabilizedLCLAlgorithm<double>> const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<ROL::ConstraintStatusTest<double>>> cl(M("Teuchos"), "RCP_ROL_ConstraintStatusTest_double_t", "");
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<ROL::ConstraintStatusTest<double>>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("null_arg") );

		cl.def( pybind11::init( [](class ROL::ConstraintStatusTest<double> * a0){ return new Teuchos::RCP<ROL::ConstraintStatusTest<double>>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<class ROL::ConstraintStatusTest<double> *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership_in") );

		cl.def( pybind11::init( [](Teuchos::RCP<ROL::ConstraintStatusTest<double>> const &o){ return new Teuchos::RCP<ROL::ConstraintStatusTest<double>>(o); } ) );
		cl.def( pybind11::init<class ROL::ConstraintStatusTest<double> *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::ConstraintStatusTest<double> *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::ConstraintStatusTest<double> *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<class ROL::ConstraintStatusTest<double> > & (Teuchos::RCP<ROL::ConstraintStatusTest<double>>::*)(const class Teuchos::RCP<class ROL::ConstraintStatusTest<double> > &)) &Teuchos::RCP<ROL::ConstraintStatusTest<double> >::operator=, "C++: Teuchos::RCP<ROL::ConstraintStatusTest<double> >::operator=(const class Teuchos::RCP<class ROL::ConstraintStatusTest<double> > &) --> class Teuchos::RCP<class ROL::ConstraintStatusTest<double> > &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<class ROL::ConstraintStatusTest<double> > & (Teuchos::RCP<ROL::ConstraintStatusTest<double>>::*)(enum Teuchos::ENull)) &Teuchos::RCP<ROL::ConstraintStatusTest<double> >::operator=, "C++: Teuchos::RCP<ROL::ConstraintStatusTest<double> >::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<class ROL::ConstraintStatusTest<double> > &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<ROL::ConstraintStatusTest<double>>::*)(class Teuchos::RCP<class ROL::ConstraintStatusTest<double> > &)) &Teuchos::RCP<ROL::ConstraintStatusTest<double> >::swap, "C++: Teuchos::RCP<ROL::ConstraintStatusTest<double> >::swap(class Teuchos::RCP<class ROL::ConstraintStatusTest<double> > &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<ROL::ConstraintStatusTest<double>>::*)() const) &Teuchos::RCP<ROL::ConstraintStatusTest<double> >::is_null, "C++: Teuchos::RCP<ROL::ConstraintStatusTest<double> >::is_null() const --> bool");
		cl.def("arrow", (class ROL::ConstraintStatusTest<double> * (Teuchos::RCP<ROL::ConstraintStatusTest<double>>::*)() const) &Teuchos::RCP<ROL::ConstraintStatusTest<double> >::operator->, "C++: Teuchos::RCP<ROL::ConstraintStatusTest<double> >::operator->() const --> class ROL::ConstraintStatusTest<double> *", pybind11::return_value_policy::automatic);
		cl.def("__mul__", (class ROL::ConstraintStatusTest<double> & (Teuchos::RCP<ROL::ConstraintStatusTest<double>>::*)() const) &Teuchos::RCP<ROL::ConstraintStatusTest<double> >::operator*, "C++: Teuchos::RCP<ROL::ConstraintStatusTest<double> >::operator*() const --> class ROL::ConstraintStatusTest<double> &", pybind11::return_value_policy::automatic);
		cl.def("get", (class ROL::ConstraintStatusTest<double> * (Teuchos::RCP<ROL::ConstraintStatusTest<double>>::*)() const) &Teuchos::RCP<ROL::ConstraintStatusTest<double> >::get, "C++: Teuchos::RCP<ROL::ConstraintStatusTest<double> >::get() const --> class ROL::ConstraintStatusTest<double> *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class ROL::ConstraintStatusTest<double> * (Teuchos::RCP<ROL::ConstraintStatusTest<double>>::*)() const) &Teuchos::RCP<ROL::ConstraintStatusTest<double> >::getRawPtr, "C++: Teuchos::RCP<ROL::ConstraintStatusTest<double> >::getRawPtr() const --> class ROL::ConstraintStatusTest<double> *", pybind11::return_value_policy::automatic);
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<ROL::ConstraintStatusTest<double>>::*)() const) &Teuchos::RCP<ROL::ConstraintStatusTest<double> >::strength, "C++: Teuchos::RCP<ROL::ConstraintStatusTest<double> >::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<ROL::ConstraintStatusTest<double>>::*)() const) &Teuchos::RCP<ROL::ConstraintStatusTest<double> >::is_valid_ptr, "C++: Teuchos::RCP<ROL::ConstraintStatusTest<double> >::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<ROL::ConstraintStatusTest<double>>::*)() const) &Teuchos::RCP<ROL::ConstraintStatusTest<double> >::strong_count, "C++: Teuchos::RCP<ROL::ConstraintStatusTest<double> >::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<ROL::ConstraintStatusTest<double>>::*)() const) &Teuchos::RCP<ROL::ConstraintStatusTest<double> >::weak_count, "C++: Teuchos::RCP<ROL::ConstraintStatusTest<double> >::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<ROL::ConstraintStatusTest<double>>::*)() const) &Teuchos::RCP<ROL::ConstraintStatusTest<double> >::total_count, "C++: Teuchos::RCP<ROL::ConstraintStatusTest<double> >::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<ROL::ConstraintStatusTest<double>>::*)()) &Teuchos::RCP<ROL::ConstraintStatusTest<double> >::set_has_ownership, "C++: Teuchos::RCP<ROL::ConstraintStatusTest<double> >::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<ROL::ConstraintStatusTest<double>>::*)() const) &Teuchos::RCP<ROL::ConstraintStatusTest<double> >::has_ownership, "C++: Teuchos::RCP<ROL::ConstraintStatusTest<double> >::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<class ROL::ConstraintStatusTest<double> > (Teuchos::RCP<ROL::ConstraintStatusTest<double>>::*)() const) &Teuchos::RCP<ROL::ConstraintStatusTest<double> >::create_weak, "C++: Teuchos::RCP<ROL::ConstraintStatusTest<double> >::create_weak() const --> class Teuchos::RCP<class ROL::ConstraintStatusTest<double> >");
		cl.def("create_strong", (class Teuchos::RCP<class ROL::ConstraintStatusTest<double> > (Teuchos::RCP<ROL::ConstraintStatusTest<double>>::*)() const) &Teuchos::RCP<ROL::ConstraintStatusTest<double> >::create_strong, "C++: Teuchos::RCP<ROL::ConstraintStatusTest<double> >::create_strong() const --> class Teuchos::RCP<class ROL::ConstraintStatusTest<double> >");
		cl.def("assert_not_null", (const class Teuchos::RCP<class ROL::ConstraintStatusTest<double> > & (Teuchos::RCP<ROL::ConstraintStatusTest<double>>::*)() const) &Teuchos::RCP<ROL::ConstraintStatusTest<double> >::assert_not_null, "C++: Teuchos::RCP<ROL::ConstraintStatusTest<double> >::assert_not_null() const --> const class Teuchos::RCP<class ROL::ConstraintStatusTest<double> > &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<class ROL::ConstraintStatusTest<double> > & (Teuchos::RCP<ROL::ConstraintStatusTest<double>>::*)() const) &Teuchos::RCP<ROL::ConstraintStatusTest<double> >::assert_valid_ptr, "C++: Teuchos::RCP<ROL::ConstraintStatusTest<double> >::assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::ConstraintStatusTest<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<class ROL::ConstraintStatusTest<double> > & (Teuchos::RCP<ROL::ConstraintStatusTest<double>>::*)() const) &Teuchos::RCP<ROL::ConstraintStatusTest<double> >::debug_assert_not_null, "C++: Teuchos::RCP<ROL::ConstraintStatusTest<double> >::debug_assert_not_null() const --> const class Teuchos::RCP<class ROL::ConstraintStatusTest<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<class ROL::ConstraintStatusTest<double> > & (Teuchos::RCP<ROL::ConstraintStatusTest<double>>::*)() const) &Teuchos::RCP<ROL::ConstraintStatusTest<double> >::debug_assert_valid_ptr, "C++: Teuchos::RCP<ROL::ConstraintStatusTest<double> >::debug_assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::ConstraintStatusTest<double> > &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<ROL::ConstraintStatusTest<double>>::*)()) &Teuchos::RCP<ROL::ConstraintStatusTest<double> >::reset, "C++: Teuchos::RCP<ROL::ConstraintStatusTest<double> >::reset() --> void");
		cl.def("access_private_ptr", (class ROL::ConstraintStatusTest<double> * (Teuchos::RCP<ROL::ConstraintStatusTest<double>>::*)() const) &Teuchos::RCP<ROL::ConstraintStatusTest<double> >::access_private_ptr, "C++: Teuchos::RCP<ROL::ConstraintStatusTest<double> >::access_private_ptr() const --> class ROL::ConstraintStatusTest<double> *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::ConstraintStatusTest<double>>::*)()) &Teuchos::RCP<ROL::ConstraintStatusTest<double> >::nonconst_access_private_node, "C++: Teuchos::RCP<ROL::ConstraintStatusTest<double> >::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::ConstraintStatusTest<double>>::*)() const) &Teuchos::RCP<ROL::ConstraintStatusTest<double> >::access_private_node, "C++: Teuchos::RCP<ROL::ConstraintStatusTest<double> >::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<ROL::ConstraintStatusTest<double>> const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<ROL::AugmentedLagrangianObjective<double>>> cl(M("Teuchos"), "RCP_ROL_AugmentedLagrangianObjective_double_t", "");
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<ROL::AugmentedLagrangianObjective<double>>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("") );

		cl.def( pybind11::init( [](class ROL::AugmentedLagrangianObjective<double> * a0){ return new Teuchos::RCP<ROL::AugmentedLagrangianObjective<double>>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<class ROL::AugmentedLagrangianObjective<double> *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership_in") );

		cl.def( pybind11::init( [](Teuchos::RCP<ROL::AugmentedLagrangianObjective<double>> const &o){ return new Teuchos::RCP<ROL::AugmentedLagrangianObjective<double>>(o); } ) );
		cl.def( pybind11::init<class ROL::AugmentedLagrangianObjective<double> *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::AugmentedLagrangianObjective<double> *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::AugmentedLagrangianObjective<double> *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<class ROL::AugmentedLagrangianObjective<double> > & (Teuchos::RCP<ROL::AugmentedLagrangianObjective<double>>::*)(const class Teuchos::RCP<class ROL::AugmentedLagrangianObjective<double> > &)) &Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::operator=, "C++: Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::operator=(const class Teuchos::RCP<class ROL::AugmentedLagrangianObjective<double> > &) --> class Teuchos::RCP<class ROL::AugmentedLagrangianObjective<double> > &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<class ROL::AugmentedLagrangianObjective<double> > & (Teuchos::RCP<ROL::AugmentedLagrangianObjective<double>>::*)(enum Teuchos::ENull)) &Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::operator=, "C++: Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<class ROL::AugmentedLagrangianObjective<double> > &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<ROL::AugmentedLagrangianObjective<double>>::*)(class Teuchos::RCP<class ROL::AugmentedLagrangianObjective<double> > &)) &Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::swap, "C++: Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::swap(class Teuchos::RCP<class ROL::AugmentedLagrangianObjective<double> > &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<ROL::AugmentedLagrangianObjective<double>>::*)() const) &Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::is_null, "C++: Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::is_null() const --> bool");
		cl.def("arrow", (class ROL::AugmentedLagrangianObjective<double> * (Teuchos::RCP<ROL::AugmentedLagrangianObjective<double>>::*)() const) &Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::operator->, "C++: Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::operator->() const --> class ROL::AugmentedLagrangianObjective<double> *", pybind11::return_value_policy::automatic);
		cl.def("__mul__", (class ROL::AugmentedLagrangianObjective<double> & (Teuchos::RCP<ROL::AugmentedLagrangianObjective<double>>::*)() const) &Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::operator*, "C++: Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::operator*() const --> class ROL::AugmentedLagrangianObjective<double> &", pybind11::return_value_policy::automatic);
		cl.def("get", (class ROL::AugmentedLagrangianObjective<double> * (Teuchos::RCP<ROL::AugmentedLagrangianObjective<double>>::*)() const) &Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::get, "C++: Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::get() const --> class ROL::AugmentedLagrangianObjective<double> *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class ROL::AugmentedLagrangianObjective<double> * (Teuchos::RCP<ROL::AugmentedLagrangianObjective<double>>::*)() const) &Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::getRawPtr, "C++: Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::getRawPtr() const --> class ROL::AugmentedLagrangianObjective<double> *", pybind11::return_value_policy::automatic);
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<ROL::AugmentedLagrangianObjective<double>>::*)() const) &Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::strength, "C++: Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<ROL::AugmentedLagrangianObjective<double>>::*)() const) &Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::is_valid_ptr, "C++: Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<ROL::AugmentedLagrangianObjective<double>>::*)() const) &Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::strong_count, "C++: Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<ROL::AugmentedLagrangianObjective<double>>::*)() const) &Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::weak_count, "C++: Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<ROL::AugmentedLagrangianObjective<double>>::*)() const) &Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::total_count, "C++: Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<ROL::AugmentedLagrangianObjective<double>>::*)()) &Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::set_has_ownership, "C++: Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<ROL::AugmentedLagrangianObjective<double>>::*)() const) &Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::has_ownership, "C++: Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<class ROL::AugmentedLagrangianObjective<double> > (Teuchos::RCP<ROL::AugmentedLagrangianObjective<double>>::*)() const) &Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::create_weak, "C++: Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::create_weak() const --> class Teuchos::RCP<class ROL::AugmentedLagrangianObjective<double> >");
		cl.def("create_strong", (class Teuchos::RCP<class ROL::AugmentedLagrangianObjective<double> > (Teuchos::RCP<ROL::AugmentedLagrangianObjective<double>>::*)() const) &Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::create_strong, "C++: Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::create_strong() const --> class Teuchos::RCP<class ROL::AugmentedLagrangianObjective<double> >");
		cl.def("assert_not_null", (const class Teuchos::RCP<class ROL::AugmentedLagrangianObjective<double> > & (Teuchos::RCP<ROL::AugmentedLagrangianObjective<double>>::*)() const) &Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::assert_not_null, "C++: Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::assert_not_null() const --> const class Teuchos::RCP<class ROL::AugmentedLagrangianObjective<double> > &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<class ROL::AugmentedLagrangianObjective<double> > & (Teuchos::RCP<ROL::AugmentedLagrangianObjective<double>>::*)() const) &Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::assert_valid_ptr, "C++: Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::AugmentedLagrangianObjective<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<class ROL::AugmentedLagrangianObjective<double> > & (Teuchos::RCP<ROL::AugmentedLagrangianObjective<double>>::*)() const) &Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::debug_assert_not_null, "C++: Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::debug_assert_not_null() const --> const class Teuchos::RCP<class ROL::AugmentedLagrangianObjective<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<class ROL::AugmentedLagrangianObjective<double> > & (Teuchos::RCP<ROL::AugmentedLagrangianObjective<double>>::*)() const) &Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::debug_assert_valid_ptr, "C++: Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::debug_assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::AugmentedLagrangianObjective<double> > &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<ROL::AugmentedLagrangianObjective<double>>::*)()) &Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::reset, "C++: Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::reset() --> void");
		cl.def("access_private_ptr", (class ROL::AugmentedLagrangianObjective<double> * (Teuchos::RCP<ROL::AugmentedLagrangianObjective<double>>::*)() const) &Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::access_private_ptr, "C++: Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::access_private_ptr() const --> class ROL::AugmentedLagrangianObjective<double> *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::AugmentedLagrangianObjective<double>>::*)()) &Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::nonconst_access_private_node, "C++: Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::AugmentedLagrangianObjective<double>>::*)() const) &Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::access_private_node, "C++: Teuchos::RCP<ROL::AugmentedLagrangianObjective<double> >::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<ROL::AugmentedLagrangianObjective<double>> const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<ROL::ElasticLinearConstraint<double>>> cl(M("Teuchos"), "RCP_ROL_ElasticLinearConstraint_double_t", "");
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<ROL::ElasticLinearConstraint<double>>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("null_arg") );

		cl.def( pybind11::init( [](class ROL::ElasticLinearConstraint<double> * a0){ return new Teuchos::RCP<ROL::ElasticLinearConstraint<double>>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<class ROL::ElasticLinearConstraint<double> *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership_in") );

		cl.def( pybind11::init( [](Teuchos::RCP<ROL::ElasticLinearConstraint<double>> const &o){ return new Teuchos::RCP<ROL::ElasticLinearConstraint<double>>(o); } ) );
		cl.def( pybind11::init<class ROL::ElasticLinearConstraint<double> *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::ElasticLinearConstraint<double> *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::ElasticLinearConstraint<double> *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<class ROL::ElasticLinearConstraint<double> > & (Teuchos::RCP<ROL::ElasticLinearConstraint<double>>::*)(const class Teuchos::RCP<class ROL::ElasticLinearConstraint<double> > &)) &Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::operator=, "C++: Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::operator=(const class Teuchos::RCP<class ROL::ElasticLinearConstraint<double> > &) --> class Teuchos::RCP<class ROL::ElasticLinearConstraint<double> > &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<class ROL::ElasticLinearConstraint<double> > & (Teuchos::RCP<ROL::ElasticLinearConstraint<double>>::*)(enum Teuchos::ENull)) &Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::operator=, "C++: Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<class ROL::ElasticLinearConstraint<double> > &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<ROL::ElasticLinearConstraint<double>>::*)(class Teuchos::RCP<class ROL::ElasticLinearConstraint<double> > &)) &Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::swap, "C++: Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::swap(class Teuchos::RCP<class ROL::ElasticLinearConstraint<double> > &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<ROL::ElasticLinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::is_null, "C++: Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::is_null() const --> bool");
		cl.def("arrow", (class ROL::ElasticLinearConstraint<double> * (Teuchos::RCP<ROL::ElasticLinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::operator->, "C++: Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::operator->() const --> class ROL::ElasticLinearConstraint<double> *", pybind11::return_value_policy::automatic);
		cl.def("__mul__", (class ROL::ElasticLinearConstraint<double> & (Teuchos::RCP<ROL::ElasticLinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::operator*, "C++: Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::operator*() const --> class ROL::ElasticLinearConstraint<double> &", pybind11::return_value_policy::automatic);
		cl.def("get", (class ROL::ElasticLinearConstraint<double> * (Teuchos::RCP<ROL::ElasticLinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::get, "C++: Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::get() const --> class ROL::ElasticLinearConstraint<double> *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class ROL::ElasticLinearConstraint<double> * (Teuchos::RCP<ROL::ElasticLinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::getRawPtr, "C++: Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::getRawPtr() const --> class ROL::ElasticLinearConstraint<double> *", pybind11::return_value_policy::automatic);
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<ROL::ElasticLinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::strength, "C++: Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<ROL::ElasticLinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::is_valid_ptr, "C++: Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<ROL::ElasticLinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::strong_count, "C++: Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<ROL::ElasticLinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::weak_count, "C++: Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<ROL::ElasticLinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::total_count, "C++: Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<ROL::ElasticLinearConstraint<double>>::*)()) &Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::set_has_ownership, "C++: Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<ROL::ElasticLinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::has_ownership, "C++: Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<class ROL::ElasticLinearConstraint<double> > (Teuchos::RCP<ROL::ElasticLinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::create_weak, "C++: Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::create_weak() const --> class Teuchos::RCP<class ROL::ElasticLinearConstraint<double> >");
		cl.def("create_strong", (class Teuchos::RCP<class ROL::ElasticLinearConstraint<double> > (Teuchos::RCP<ROL::ElasticLinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::create_strong, "C++: Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::create_strong() const --> class Teuchos::RCP<class ROL::ElasticLinearConstraint<double> >");
		cl.def("assert_not_null", (const class Teuchos::RCP<class ROL::ElasticLinearConstraint<double> > & (Teuchos::RCP<ROL::ElasticLinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::assert_not_null, "C++: Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::assert_not_null() const --> const class Teuchos::RCP<class ROL::ElasticLinearConstraint<double> > &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<class ROL::ElasticLinearConstraint<double> > & (Teuchos::RCP<ROL::ElasticLinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::assert_valid_ptr, "C++: Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::ElasticLinearConstraint<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<class ROL::ElasticLinearConstraint<double> > & (Teuchos::RCP<ROL::ElasticLinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::debug_assert_not_null, "C++: Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::debug_assert_not_null() const --> const class Teuchos::RCP<class ROL::ElasticLinearConstraint<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<class ROL::ElasticLinearConstraint<double> > & (Teuchos::RCP<ROL::ElasticLinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::debug_assert_valid_ptr, "C++: Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::debug_assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::ElasticLinearConstraint<double> > &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<ROL::ElasticLinearConstraint<double>>::*)()) &Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::reset, "C++: Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::reset() --> void");
		cl.def("access_private_ptr", (class ROL::ElasticLinearConstraint<double> * (Teuchos::RCP<ROL::ElasticLinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::access_private_ptr, "C++: Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::access_private_ptr() const --> class ROL::ElasticLinearConstraint<double> *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::ElasticLinearConstraint<double>>::*)()) &Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::nonconst_access_private_node, "C++: Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::ElasticLinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::access_private_node, "C++: Teuchos::RCP<ROL::ElasticLinearConstraint<double> >::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<ROL::ElasticLinearConstraint<double>> const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<ROL::Bounds<double>>> cl(M("Teuchos"), "RCP_ROL_Bounds_double_t", "");
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<ROL::Bounds<double>>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("null_arg") );

		cl.def( pybind11::init( [](class ROL::Bounds<double> * a0){ return new Teuchos::RCP<ROL::Bounds<double>>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<class ROL::Bounds<double> *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership_in") );

		cl.def( pybind11::init( [](Teuchos::RCP<ROL::Bounds<double>> const &o){ return new Teuchos::RCP<ROL::Bounds<double>>(o); } ) );
		cl.def( pybind11::init<class ROL::Bounds<double> *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::Bounds<double> *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::Bounds<double> *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<class ROL::Bounds<double> > & (Teuchos::RCP<ROL::Bounds<double>>::*)(const class Teuchos::RCP<class ROL::Bounds<double> > &)) &Teuchos::RCP<ROL::Bounds<double> >::operator=, "C++: Teuchos::RCP<ROL::Bounds<double> >::operator=(const class Teuchos::RCP<class ROL::Bounds<double> > &) --> class Teuchos::RCP<class ROL::Bounds<double> > &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<class ROL::Bounds<double> > & (Teuchos::RCP<ROL::Bounds<double>>::*)(enum Teuchos::ENull)) &Teuchos::RCP<ROL::Bounds<double> >::operator=, "C++: Teuchos::RCP<ROL::Bounds<double> >::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<class ROL::Bounds<double> > &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<ROL::Bounds<double>>::*)(class Teuchos::RCP<class ROL::Bounds<double> > &)) &Teuchos::RCP<ROL::Bounds<double> >::swap, "C++: Teuchos::RCP<ROL::Bounds<double> >::swap(class Teuchos::RCP<class ROL::Bounds<double> > &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<ROL::Bounds<double>>::*)() const) &Teuchos::RCP<ROL::Bounds<double> >::is_null, "C++: Teuchos::RCP<ROL::Bounds<double> >::is_null() const --> bool");
		cl.def("arrow", (class ROL::Bounds<double> * (Teuchos::RCP<ROL::Bounds<double>>::*)() const) &Teuchos::RCP<ROL::Bounds<double> >::operator->, "C++: Teuchos::RCP<ROL::Bounds<double> >::operator->() const --> class ROL::Bounds<double> *", pybind11::return_value_policy::automatic);
		cl.def("__mul__", (class ROL::Bounds<double> & (Teuchos::RCP<ROL::Bounds<double>>::*)() const) &Teuchos::RCP<ROL::Bounds<double> >::operator*, "C++: Teuchos::RCP<ROL::Bounds<double> >::operator*() const --> class ROL::Bounds<double> &", pybind11::return_value_policy::automatic);
		cl.def("get", (class ROL::Bounds<double> * (Teuchos::RCP<ROL::Bounds<double>>::*)() const) &Teuchos::RCP<ROL::Bounds<double> >::get, "C++: Teuchos::RCP<ROL::Bounds<double> >::get() const --> class ROL::Bounds<double> *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class ROL::Bounds<double> * (Teuchos::RCP<ROL::Bounds<double>>::*)() const) &Teuchos::RCP<ROL::Bounds<double> >::getRawPtr, "C++: Teuchos::RCP<ROL::Bounds<double> >::getRawPtr() const --> class ROL::Bounds<double> *", pybind11::return_value_policy::automatic);
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<ROL::Bounds<double>>::*)() const) &Teuchos::RCP<ROL::Bounds<double> >::strength, "C++: Teuchos::RCP<ROL::Bounds<double> >::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<ROL::Bounds<double>>::*)() const) &Teuchos::RCP<ROL::Bounds<double> >::is_valid_ptr, "C++: Teuchos::RCP<ROL::Bounds<double> >::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<ROL::Bounds<double>>::*)() const) &Teuchos::RCP<ROL::Bounds<double> >::strong_count, "C++: Teuchos::RCP<ROL::Bounds<double> >::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<ROL::Bounds<double>>::*)() const) &Teuchos::RCP<ROL::Bounds<double> >::weak_count, "C++: Teuchos::RCP<ROL::Bounds<double> >::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<ROL::Bounds<double>>::*)() const) &Teuchos::RCP<ROL::Bounds<double> >::total_count, "C++: Teuchos::RCP<ROL::Bounds<double> >::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<ROL::Bounds<double>>::*)()) &Teuchos::RCP<ROL::Bounds<double> >::set_has_ownership, "C++: Teuchos::RCP<ROL::Bounds<double> >::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<ROL::Bounds<double>>::*)() const) &Teuchos::RCP<ROL::Bounds<double> >::has_ownership, "C++: Teuchos::RCP<ROL::Bounds<double> >::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<class ROL::Bounds<double> > (Teuchos::RCP<ROL::Bounds<double>>::*)() const) &Teuchos::RCP<ROL::Bounds<double> >::create_weak, "C++: Teuchos::RCP<ROL::Bounds<double> >::create_weak() const --> class Teuchos::RCP<class ROL::Bounds<double> >");
		cl.def("create_strong", (class Teuchos::RCP<class ROL::Bounds<double> > (Teuchos::RCP<ROL::Bounds<double>>::*)() const) &Teuchos::RCP<ROL::Bounds<double> >::create_strong, "C++: Teuchos::RCP<ROL::Bounds<double> >::create_strong() const --> class Teuchos::RCP<class ROL::Bounds<double> >");
		cl.def("assert_not_null", (const class Teuchos::RCP<class ROL::Bounds<double> > & (Teuchos::RCP<ROL::Bounds<double>>::*)() const) &Teuchos::RCP<ROL::Bounds<double> >::assert_not_null, "C++: Teuchos::RCP<ROL::Bounds<double> >::assert_not_null() const --> const class Teuchos::RCP<class ROL::Bounds<double> > &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<class ROL::Bounds<double> > & (Teuchos::RCP<ROL::Bounds<double>>::*)() const) &Teuchos::RCP<ROL::Bounds<double> >::assert_valid_ptr, "C++: Teuchos::RCP<ROL::Bounds<double> >::assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::Bounds<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<class ROL::Bounds<double> > & (Teuchos::RCP<ROL::Bounds<double>>::*)() const) &Teuchos::RCP<ROL::Bounds<double> >::debug_assert_not_null, "C++: Teuchos::RCP<ROL::Bounds<double> >::debug_assert_not_null() const --> const class Teuchos::RCP<class ROL::Bounds<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<class ROL::Bounds<double> > & (Teuchos::RCP<ROL::Bounds<double>>::*)() const) &Teuchos::RCP<ROL::Bounds<double> >::debug_assert_valid_ptr, "C++: Teuchos::RCP<ROL::Bounds<double> >::debug_assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::Bounds<double> > &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<ROL::Bounds<double>>::*)()) &Teuchos::RCP<ROL::Bounds<double> >::reset, "C++: Teuchos::RCP<ROL::Bounds<double> >::reset() --> void");
		cl.def("access_private_ptr", (class ROL::Bounds<double> * (Teuchos::RCP<ROL::Bounds<double>>::*)() const) &Teuchos::RCP<ROL::Bounds<double> >::access_private_ptr, "C++: Teuchos::RCP<ROL::Bounds<double> >::access_private_ptr() const --> class ROL::Bounds<double> *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::Bounds<double>>::*)()) &Teuchos::RCP<ROL::Bounds<double> >::nonconst_access_private_node, "C++: Teuchos::RCP<ROL::Bounds<double> >::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::Bounds<double>>::*)() const) &Teuchos::RCP<ROL::Bounds<double> >::access_private_node, "C++: Teuchos::RCP<ROL::Bounds<double> >::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<ROL::Bounds<double>> const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<ROL::BoundConstraint_Partitioned<double>>> cl(M("Teuchos"), "RCP_ROL_BoundConstraint_Partitioned_double_t", "");
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<ROL::BoundConstraint_Partitioned<double>>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("null_arg") );

		cl.def( pybind11::init( [](class ROL::BoundConstraint_Partitioned<double> * a0){ return new Teuchos::RCP<ROL::BoundConstraint_Partitioned<double>>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<class ROL::BoundConstraint_Partitioned<double> *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership_in") );

		cl.def( pybind11::init( [](Teuchos::RCP<ROL::BoundConstraint_Partitioned<double>> const &o){ return new Teuchos::RCP<ROL::BoundConstraint_Partitioned<double>>(o); } ) );
		cl.def( pybind11::init<class ROL::BoundConstraint_Partitioned<double> *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::BoundConstraint_Partitioned<double> *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::BoundConstraint_Partitioned<double> *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<class ROL::BoundConstraint_Partitioned<double> > & (Teuchos::RCP<ROL::BoundConstraint_Partitioned<double>>::*)(const class Teuchos::RCP<class ROL::BoundConstraint_Partitioned<double> > &)) &Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::operator=, "C++: Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::operator=(const class Teuchos::RCP<class ROL::BoundConstraint_Partitioned<double> > &) --> class Teuchos::RCP<class ROL::BoundConstraint_Partitioned<double> > &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<class ROL::BoundConstraint_Partitioned<double> > & (Teuchos::RCP<ROL::BoundConstraint_Partitioned<double>>::*)(enum Teuchos::ENull)) &Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::operator=, "C++: Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<class ROL::BoundConstraint_Partitioned<double> > &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<ROL::BoundConstraint_Partitioned<double>>::*)(class Teuchos::RCP<class ROL::BoundConstraint_Partitioned<double> > &)) &Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::swap, "C++: Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::swap(class Teuchos::RCP<class ROL::BoundConstraint_Partitioned<double> > &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<ROL::BoundConstraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::is_null, "C++: Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::is_null() const --> bool");
		cl.def("arrow", (class ROL::BoundConstraint_Partitioned<double> * (Teuchos::RCP<ROL::BoundConstraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::operator->, "C++: Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::operator->() const --> class ROL::BoundConstraint_Partitioned<double> *", pybind11::return_value_policy::automatic);
		cl.def("__mul__", (class ROL::BoundConstraint_Partitioned<double> & (Teuchos::RCP<ROL::BoundConstraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::operator*, "C++: Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::operator*() const --> class ROL::BoundConstraint_Partitioned<double> &", pybind11::return_value_policy::automatic);
		cl.def("get", (class ROL::BoundConstraint_Partitioned<double> * (Teuchos::RCP<ROL::BoundConstraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::get, "C++: Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::get() const --> class ROL::BoundConstraint_Partitioned<double> *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class ROL::BoundConstraint_Partitioned<double> * (Teuchos::RCP<ROL::BoundConstraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::getRawPtr, "C++: Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::getRawPtr() const --> class ROL::BoundConstraint_Partitioned<double> *", pybind11::return_value_policy::automatic);
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<ROL::BoundConstraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::strength, "C++: Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<ROL::BoundConstraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::is_valid_ptr, "C++: Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<ROL::BoundConstraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::strong_count, "C++: Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<ROL::BoundConstraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::weak_count, "C++: Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<ROL::BoundConstraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::total_count, "C++: Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<ROL::BoundConstraint_Partitioned<double>>::*)()) &Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::set_has_ownership, "C++: Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<ROL::BoundConstraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::has_ownership, "C++: Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<class ROL::BoundConstraint_Partitioned<double> > (Teuchos::RCP<ROL::BoundConstraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::create_weak, "C++: Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::create_weak() const --> class Teuchos::RCP<class ROL::BoundConstraint_Partitioned<double> >");
		cl.def("create_strong", (class Teuchos::RCP<class ROL::BoundConstraint_Partitioned<double> > (Teuchos::RCP<ROL::BoundConstraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::create_strong, "C++: Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::create_strong() const --> class Teuchos::RCP<class ROL::BoundConstraint_Partitioned<double> >");
		cl.def("assert_not_null", (const class Teuchos::RCP<class ROL::BoundConstraint_Partitioned<double> > & (Teuchos::RCP<ROL::BoundConstraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::assert_not_null, "C++: Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::assert_not_null() const --> const class Teuchos::RCP<class ROL::BoundConstraint_Partitioned<double> > &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<class ROL::BoundConstraint_Partitioned<double> > & (Teuchos::RCP<ROL::BoundConstraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::assert_valid_ptr, "C++: Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::BoundConstraint_Partitioned<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<class ROL::BoundConstraint_Partitioned<double> > & (Teuchos::RCP<ROL::BoundConstraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::debug_assert_not_null, "C++: Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::debug_assert_not_null() const --> const class Teuchos::RCP<class ROL::BoundConstraint_Partitioned<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<class ROL::BoundConstraint_Partitioned<double> > & (Teuchos::RCP<ROL::BoundConstraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::debug_assert_valid_ptr, "C++: Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::debug_assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::BoundConstraint_Partitioned<double> > &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<ROL::BoundConstraint_Partitioned<double>>::*)()) &Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::reset, "C++: Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::reset() --> void");
		cl.def("access_private_ptr", (class ROL::BoundConstraint_Partitioned<double> * (Teuchos::RCP<ROL::BoundConstraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::access_private_ptr, "C++: Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::access_private_ptr() const --> class ROL::BoundConstraint_Partitioned<double> *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::BoundConstraint_Partitioned<double>>::*)()) &Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::nonconst_access_private_node, "C++: Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::BoundConstraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::access_private_node, "C++: Teuchos::RCP<ROL::BoundConstraint_Partitioned<double> >::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<ROL::BoundConstraint_Partitioned<double>> const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<ROL::ElasticObjective<double>>> cl(M("Teuchos"), "RCP_ROL_ElasticObjective_double_t", "");
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<ROL::ElasticObjective<double>>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("null_arg") );

		cl.def( pybind11::init( [](class ROL::ElasticObjective<double> * a0){ return new Teuchos::RCP<ROL::ElasticObjective<double>>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<class ROL::ElasticObjective<double> *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership") );

		cl.def( pybind11::init( [](Teuchos::RCP<ROL::ElasticObjective<double>> const &o){ return new Teuchos::RCP<ROL::ElasticObjective<double>>(o); } ) );
		cl.def( pybind11::init<class ROL::ElasticObjective<double> *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::ElasticObjective<double> *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::ElasticObjective<double> *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<class ROL::ElasticObjective<double> > & (Teuchos::RCP<ROL::ElasticObjective<double>>::*)(const class Teuchos::RCP<class ROL::ElasticObjective<double> > &)) &Teuchos::RCP<ROL::ElasticObjective<double> >::operator=, "C++: Teuchos::RCP<ROL::ElasticObjective<double> >::operator=(const class Teuchos::RCP<class ROL::ElasticObjective<double> > &) --> class Teuchos::RCP<class ROL::ElasticObjective<double> > &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<class ROL::ElasticObjective<double> > & (Teuchos::RCP<ROL::ElasticObjective<double>>::*)(enum Teuchos::ENull)) &Teuchos::RCP<ROL::ElasticObjective<double> >::operator=, "C++: Teuchos::RCP<ROL::ElasticObjective<double> >::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<class ROL::ElasticObjective<double> > &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<ROL::ElasticObjective<double>>::*)(class Teuchos::RCP<class ROL::ElasticObjective<double> > &)) &Teuchos::RCP<ROL::ElasticObjective<double> >::swap, "C++: Teuchos::RCP<ROL::ElasticObjective<double> >::swap(class Teuchos::RCP<class ROL::ElasticObjective<double> > &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<ROL::ElasticObjective<double>>::*)() const) &Teuchos::RCP<ROL::ElasticObjective<double> >::is_null, "C++: Teuchos::RCP<ROL::ElasticObjective<double> >::is_null() const --> bool");
		cl.def("arrow", (class ROL::ElasticObjective<double> * (Teuchos::RCP<ROL::ElasticObjective<double>>::*)() const) &Teuchos::RCP<ROL::ElasticObjective<double> >::operator->, "C++: Teuchos::RCP<ROL::ElasticObjective<double> >::operator->() const --> class ROL::ElasticObjective<double> *", pybind11::return_value_policy::automatic);
		cl.def("__mul__", (class ROL::ElasticObjective<double> & (Teuchos::RCP<ROL::ElasticObjective<double>>::*)() const) &Teuchos::RCP<ROL::ElasticObjective<double> >::operator*, "C++: Teuchos::RCP<ROL::ElasticObjective<double> >::operator*() const --> class ROL::ElasticObjective<double> &", pybind11::return_value_policy::automatic);
		cl.def("get", (class ROL::ElasticObjective<double> * (Teuchos::RCP<ROL::ElasticObjective<double>>::*)() const) &Teuchos::RCP<ROL::ElasticObjective<double> >::get, "C++: Teuchos::RCP<ROL::ElasticObjective<double> >::get() const --> class ROL::ElasticObjective<double> *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class ROL::ElasticObjective<double> * (Teuchos::RCP<ROL::ElasticObjective<double>>::*)() const) &Teuchos::RCP<ROL::ElasticObjective<double> >::getRawPtr, "C++: Teuchos::RCP<ROL::ElasticObjective<double> >::getRawPtr() const --> class ROL::ElasticObjective<double> *", pybind11::return_value_policy::automatic);
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<ROL::ElasticObjective<double>>::*)() const) &Teuchos::RCP<ROL::ElasticObjective<double> >::strength, "C++: Teuchos::RCP<ROL::ElasticObjective<double> >::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<ROL::ElasticObjective<double>>::*)() const) &Teuchos::RCP<ROL::ElasticObjective<double> >::is_valid_ptr, "C++: Teuchos::RCP<ROL::ElasticObjective<double> >::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<ROL::ElasticObjective<double>>::*)() const) &Teuchos::RCP<ROL::ElasticObjective<double> >::strong_count, "C++: Teuchos::RCP<ROL::ElasticObjective<double> >::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<ROL::ElasticObjective<double>>::*)() const) &Teuchos::RCP<ROL::ElasticObjective<double> >::weak_count, "C++: Teuchos::RCP<ROL::ElasticObjective<double> >::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<ROL::ElasticObjective<double>>::*)() const) &Teuchos::RCP<ROL::ElasticObjective<double> >::total_count, "C++: Teuchos::RCP<ROL::ElasticObjective<double> >::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<ROL::ElasticObjective<double>>::*)()) &Teuchos::RCP<ROL::ElasticObjective<double> >::set_has_ownership, "C++: Teuchos::RCP<ROL::ElasticObjective<double> >::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<ROL::ElasticObjective<double>>::*)() const) &Teuchos::RCP<ROL::ElasticObjective<double> >::has_ownership, "C++: Teuchos::RCP<ROL::ElasticObjective<double> >::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<class ROL::ElasticObjective<double> > (Teuchos::RCP<ROL::ElasticObjective<double>>::*)() const) &Teuchos::RCP<ROL::ElasticObjective<double> >::create_weak, "C++: Teuchos::RCP<ROL::ElasticObjective<double> >::create_weak() const --> class Teuchos::RCP<class ROL::ElasticObjective<double> >");
		cl.def("create_strong", (class Teuchos::RCP<class ROL::ElasticObjective<double> > (Teuchos::RCP<ROL::ElasticObjective<double>>::*)() const) &Teuchos::RCP<ROL::ElasticObjective<double> >::create_strong, "C++: Teuchos::RCP<ROL::ElasticObjective<double> >::create_strong() const --> class Teuchos::RCP<class ROL::ElasticObjective<double> >");
		cl.def("assert_not_null", (const class Teuchos::RCP<class ROL::ElasticObjective<double> > & (Teuchos::RCP<ROL::ElasticObjective<double>>::*)() const) &Teuchos::RCP<ROL::ElasticObjective<double> >::assert_not_null, "C++: Teuchos::RCP<ROL::ElasticObjective<double> >::assert_not_null() const --> const class Teuchos::RCP<class ROL::ElasticObjective<double> > &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<class ROL::ElasticObjective<double> > & (Teuchos::RCP<ROL::ElasticObjective<double>>::*)() const) &Teuchos::RCP<ROL::ElasticObjective<double> >::assert_valid_ptr, "C++: Teuchos::RCP<ROL::ElasticObjective<double> >::assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::ElasticObjective<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<class ROL::ElasticObjective<double> > & (Teuchos::RCP<ROL::ElasticObjective<double>>::*)() const) &Teuchos::RCP<ROL::ElasticObjective<double> >::debug_assert_not_null, "C++: Teuchos::RCP<ROL::ElasticObjective<double> >::debug_assert_not_null() const --> const class Teuchos::RCP<class ROL::ElasticObjective<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<class ROL::ElasticObjective<double> > & (Teuchos::RCP<ROL::ElasticObjective<double>>::*)() const) &Teuchos::RCP<ROL::ElasticObjective<double> >::debug_assert_valid_ptr, "C++: Teuchos::RCP<ROL::ElasticObjective<double> >::debug_assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::ElasticObjective<double> > &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<ROL::ElasticObjective<double>>::*)()) &Teuchos::RCP<ROL::ElasticObjective<double> >::reset, "C++: Teuchos::RCP<ROL::ElasticObjective<double> >::reset() --> void");
		cl.def("access_private_ptr", (class ROL::ElasticObjective<double> * (Teuchos::RCP<ROL::ElasticObjective<double>>::*)() const) &Teuchos::RCP<ROL::ElasticObjective<double> >::access_private_ptr, "C++: Teuchos::RCP<ROL::ElasticObjective<double> >::access_private_ptr() const --> class ROL::ElasticObjective<double> *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::ElasticObjective<double>>::*)()) &Teuchos::RCP<ROL::ElasticObjective<double> >::nonconst_access_private_node, "C++: Teuchos::RCP<ROL::ElasticObjective<double> >::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::ElasticObjective<double>>::*)() const) &Teuchos::RCP<ROL::ElasticObjective<double> >::access_private_node, "C++: Teuchos::RCP<ROL::ElasticObjective<double> >::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<ROL::ElasticObjective<double>> const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double>>> cl(M("Teuchos"), "RCP_ROL_TypeG_AugmentedLagrangianAlgorithm_double_t", "");
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double>>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("null_arg") );

		cl.def( pybind11::init( [](class ROL::TypeG::AugmentedLagrangianAlgorithm<double> * a0){ return new Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double>>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<class ROL::TypeG::AugmentedLagrangianAlgorithm<double> *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership_in") );

		cl.def( pybind11::init( [](Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double>> const &o){ return new Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double>>(o); } ) );
		cl.def( pybind11::init<class ROL::TypeG::AugmentedLagrangianAlgorithm<double> *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::TypeG::AugmentedLagrangianAlgorithm<double> *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::TypeG::AugmentedLagrangianAlgorithm<double> *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<class ROL::TypeG::AugmentedLagrangianAlgorithm<double> > & (Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double>>::*)(const class Teuchos::RCP<class ROL::TypeG::AugmentedLagrangianAlgorithm<double> > &)) &Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::operator=, "C++: Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::operator=(const class Teuchos::RCP<class ROL::TypeG::AugmentedLagrangianAlgorithm<double> > &) --> class Teuchos::RCP<class ROL::TypeG::AugmentedLagrangianAlgorithm<double> > &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<class ROL::TypeG::AugmentedLagrangianAlgorithm<double> > & (Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double>>::*)(enum Teuchos::ENull)) &Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::operator=, "C++: Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<class ROL::TypeG::AugmentedLagrangianAlgorithm<double> > &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double>>::*)(class Teuchos::RCP<class ROL::TypeG::AugmentedLagrangianAlgorithm<double> > &)) &Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::swap, "C++: Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::swap(class Teuchos::RCP<class ROL::TypeG::AugmentedLagrangianAlgorithm<double> > &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::is_null, "C++: Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::is_null() const --> bool");
		cl.def("arrow", (class ROL::TypeG::AugmentedLagrangianAlgorithm<double> * (Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::operator->, "C++: Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::operator->() const --> class ROL::TypeG::AugmentedLagrangianAlgorithm<double> *", pybind11::return_value_policy::automatic);
		cl.def("__mul__", (class ROL::TypeG::AugmentedLagrangianAlgorithm<double> & (Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::operator*, "C++: Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::operator*() const --> class ROL::TypeG::AugmentedLagrangianAlgorithm<double> &", pybind11::return_value_policy::automatic);
		cl.def("get", (class ROL::TypeG::AugmentedLagrangianAlgorithm<double> * (Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::get, "C++: Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::get() const --> class ROL::TypeG::AugmentedLagrangianAlgorithm<double> *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class ROL::TypeG::AugmentedLagrangianAlgorithm<double> * (Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::getRawPtr, "C++: Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::getRawPtr() const --> class ROL::TypeG::AugmentedLagrangianAlgorithm<double> *", pybind11::return_value_policy::automatic);
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::strength, "C++: Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::is_valid_ptr, "C++: Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::strong_count, "C++: Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::weak_count, "C++: Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::total_count, "C++: Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double>>::*)()) &Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::set_has_ownership, "C++: Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::has_ownership, "C++: Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<class ROL::TypeG::AugmentedLagrangianAlgorithm<double> > (Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::create_weak, "C++: Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::create_weak() const --> class Teuchos::RCP<class ROL::TypeG::AugmentedLagrangianAlgorithm<double> >");
		cl.def("create_strong", (class Teuchos::RCP<class ROL::TypeG::AugmentedLagrangianAlgorithm<double> > (Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::create_strong, "C++: Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::create_strong() const --> class Teuchos::RCP<class ROL::TypeG::AugmentedLagrangianAlgorithm<double> >");
		cl.def("assert_not_null", (const class Teuchos::RCP<class ROL::TypeG::AugmentedLagrangianAlgorithm<double> > & (Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::assert_not_null, "C++: Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::assert_not_null() const --> const class Teuchos::RCP<class ROL::TypeG::AugmentedLagrangianAlgorithm<double> > &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<class ROL::TypeG::AugmentedLagrangianAlgorithm<double> > & (Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::assert_valid_ptr, "C++: Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::TypeG::AugmentedLagrangianAlgorithm<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<class ROL::TypeG::AugmentedLagrangianAlgorithm<double> > & (Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::debug_assert_not_null, "C++: Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::debug_assert_not_null() const --> const class Teuchos::RCP<class ROL::TypeG::AugmentedLagrangianAlgorithm<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<class ROL::TypeG::AugmentedLagrangianAlgorithm<double> > & (Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::debug_assert_valid_ptr, "C++: Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::debug_assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::TypeG::AugmentedLagrangianAlgorithm<double> > &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double>>::*)()) &Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::reset, "C++: Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::reset() --> void");
		cl.def("access_private_ptr", (class ROL::TypeG::AugmentedLagrangianAlgorithm<double> * (Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::access_private_ptr, "C++: Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::access_private_ptr() const --> class ROL::TypeG::AugmentedLagrangianAlgorithm<double> *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double>>::*)()) &Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::nonconst_access_private_node, "C++: Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::access_private_node, "C++: Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double> >::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double>> const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double>>> cl(M("Teuchos"), "RCP_ROL_TypeG_MoreauYosidaAlgorithm_double_t", "");
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double>>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("null_arg") );

		cl.def( pybind11::init( [](class ROL::TypeG::MoreauYosidaAlgorithm<double> * a0){ return new Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double>>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<class ROL::TypeG::MoreauYosidaAlgorithm<double> *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership_in") );

		cl.def( pybind11::init( [](Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double>> const &o){ return new Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double>>(o); } ) );
		cl.def( pybind11::init<class ROL::TypeG::MoreauYosidaAlgorithm<double> *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::TypeG::MoreauYosidaAlgorithm<double> *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::TypeG::MoreauYosidaAlgorithm<double> *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<class ROL::TypeG::MoreauYosidaAlgorithm<double> > & (Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double>>::*)(const class Teuchos::RCP<class ROL::TypeG::MoreauYosidaAlgorithm<double> > &)) &Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::operator=, "C++: Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::operator=(const class Teuchos::RCP<class ROL::TypeG::MoreauYosidaAlgorithm<double> > &) --> class Teuchos::RCP<class ROL::TypeG::MoreauYosidaAlgorithm<double> > &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<class ROL::TypeG::MoreauYosidaAlgorithm<double> > & (Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double>>::*)(enum Teuchos::ENull)) &Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::operator=, "C++: Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<class ROL::TypeG::MoreauYosidaAlgorithm<double> > &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double>>::*)(class Teuchos::RCP<class ROL::TypeG::MoreauYosidaAlgorithm<double> > &)) &Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::swap, "C++: Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::swap(class Teuchos::RCP<class ROL::TypeG::MoreauYosidaAlgorithm<double> > &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::is_null, "C++: Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::is_null() const --> bool");
		cl.def("arrow", (class ROL::TypeG::MoreauYosidaAlgorithm<double> * (Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::operator->, "C++: Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::operator->() const --> class ROL::TypeG::MoreauYosidaAlgorithm<double> *", pybind11::return_value_policy::automatic);
		cl.def("__mul__", (class ROL::TypeG::MoreauYosidaAlgorithm<double> & (Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::operator*, "C++: Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::operator*() const --> class ROL::TypeG::MoreauYosidaAlgorithm<double> &", pybind11::return_value_policy::automatic);
		cl.def("get", (class ROL::TypeG::MoreauYosidaAlgorithm<double> * (Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::get, "C++: Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::get() const --> class ROL::TypeG::MoreauYosidaAlgorithm<double> *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class ROL::TypeG::MoreauYosidaAlgorithm<double> * (Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::getRawPtr, "C++: Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::getRawPtr() const --> class ROL::TypeG::MoreauYosidaAlgorithm<double> *", pybind11::return_value_policy::automatic);
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::strength, "C++: Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::is_valid_ptr, "C++: Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::strong_count, "C++: Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::weak_count, "C++: Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::total_count, "C++: Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double>>::*)()) &Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::set_has_ownership, "C++: Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::has_ownership, "C++: Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<class ROL::TypeG::MoreauYosidaAlgorithm<double> > (Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::create_weak, "C++: Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::create_weak() const --> class Teuchos::RCP<class ROL::TypeG::MoreauYosidaAlgorithm<double> >");
		cl.def("create_strong", (class Teuchos::RCP<class ROL::TypeG::MoreauYosidaAlgorithm<double> > (Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::create_strong, "C++: Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::create_strong() const --> class Teuchos::RCP<class ROL::TypeG::MoreauYosidaAlgorithm<double> >");
		cl.def("assert_not_null", (const class Teuchos::RCP<class ROL::TypeG::MoreauYosidaAlgorithm<double> > & (Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::assert_not_null, "C++: Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::assert_not_null() const --> const class Teuchos::RCP<class ROL::TypeG::MoreauYosidaAlgorithm<double> > &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<class ROL::TypeG::MoreauYosidaAlgorithm<double> > & (Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::assert_valid_ptr, "C++: Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::TypeG::MoreauYosidaAlgorithm<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<class ROL::TypeG::MoreauYosidaAlgorithm<double> > & (Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::debug_assert_not_null, "C++: Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::debug_assert_not_null() const --> const class Teuchos::RCP<class ROL::TypeG::MoreauYosidaAlgorithm<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<class ROL::TypeG::MoreauYosidaAlgorithm<double> > & (Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::debug_assert_valid_ptr, "C++: Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::debug_assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::TypeG::MoreauYosidaAlgorithm<double> > &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double>>::*)()) &Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::reset, "C++: Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::reset() --> void");
		cl.def("access_private_ptr", (class ROL::TypeG::MoreauYosidaAlgorithm<double> * (Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::access_private_ptr, "C++: Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::access_private_ptr() const --> class ROL::TypeG::MoreauYosidaAlgorithm<double> *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double>>::*)()) &Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::nonconst_access_private_node, "C++: Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::access_private_node, "C++: Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double> >::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double>> const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double>>> cl(M("Teuchos"), "RCP_ROL_TypeG_InteriorPointAlgorithm_double_t", "");
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double>>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("null_arg") );

		cl.def( pybind11::init( [](class ROL::TypeG::InteriorPointAlgorithm<double> * a0){ return new Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double>>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<class ROL::TypeG::InteriorPointAlgorithm<double> *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership_in") );

		cl.def( pybind11::init( [](Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double>> const &o){ return new Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double>>(o); } ) );
		cl.def( pybind11::init<class ROL::TypeG::InteriorPointAlgorithm<double> *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::TypeG::InteriorPointAlgorithm<double> *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::TypeG::InteriorPointAlgorithm<double> *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<class ROL::TypeG::InteriorPointAlgorithm<double> > & (Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double>>::*)(const class Teuchos::RCP<class ROL::TypeG::InteriorPointAlgorithm<double> > &)) &Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::operator=, "C++: Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::operator=(const class Teuchos::RCP<class ROL::TypeG::InteriorPointAlgorithm<double> > &) --> class Teuchos::RCP<class ROL::TypeG::InteriorPointAlgorithm<double> > &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<class ROL::TypeG::InteriorPointAlgorithm<double> > & (Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double>>::*)(enum Teuchos::ENull)) &Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::operator=, "C++: Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<class ROL::TypeG::InteriorPointAlgorithm<double> > &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double>>::*)(class Teuchos::RCP<class ROL::TypeG::InteriorPointAlgorithm<double> > &)) &Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::swap, "C++: Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::swap(class Teuchos::RCP<class ROL::TypeG::InteriorPointAlgorithm<double> > &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::is_null, "C++: Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::is_null() const --> bool");
		cl.def("arrow", (class ROL::TypeG::InteriorPointAlgorithm<double> * (Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::operator->, "C++: Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::operator->() const --> class ROL::TypeG::InteriorPointAlgorithm<double> *", pybind11::return_value_policy::automatic);
		cl.def("__mul__", (class ROL::TypeG::InteriorPointAlgorithm<double> & (Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::operator*, "C++: Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::operator*() const --> class ROL::TypeG::InteriorPointAlgorithm<double> &", pybind11::return_value_policy::automatic);
		cl.def("get", (class ROL::TypeG::InteriorPointAlgorithm<double> * (Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::get, "C++: Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::get() const --> class ROL::TypeG::InteriorPointAlgorithm<double> *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class ROL::TypeG::InteriorPointAlgorithm<double> * (Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::getRawPtr, "C++: Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::getRawPtr() const --> class ROL::TypeG::InteriorPointAlgorithm<double> *", pybind11::return_value_policy::automatic);
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::strength, "C++: Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::is_valid_ptr, "C++: Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::strong_count, "C++: Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::weak_count, "C++: Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::total_count, "C++: Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double>>::*)()) &Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::set_has_ownership, "C++: Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::has_ownership, "C++: Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<class ROL::TypeG::InteriorPointAlgorithm<double> > (Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::create_weak, "C++: Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::create_weak() const --> class Teuchos::RCP<class ROL::TypeG::InteriorPointAlgorithm<double> >");
		cl.def("create_strong", (class Teuchos::RCP<class ROL::TypeG::InteriorPointAlgorithm<double> > (Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::create_strong, "C++: Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::create_strong() const --> class Teuchos::RCP<class ROL::TypeG::InteriorPointAlgorithm<double> >");
		cl.def("assert_not_null", (const class Teuchos::RCP<class ROL::TypeG::InteriorPointAlgorithm<double> > & (Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::assert_not_null, "C++: Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::assert_not_null() const --> const class Teuchos::RCP<class ROL::TypeG::InteriorPointAlgorithm<double> > &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<class ROL::TypeG::InteriorPointAlgorithm<double> > & (Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::assert_valid_ptr, "C++: Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::TypeG::InteriorPointAlgorithm<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<class ROL::TypeG::InteriorPointAlgorithm<double> > & (Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::debug_assert_not_null, "C++: Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::debug_assert_not_null() const --> const class Teuchos::RCP<class ROL::TypeG::InteriorPointAlgorithm<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<class ROL::TypeG::InteriorPointAlgorithm<double> > & (Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::debug_assert_valid_ptr, "C++: Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::debug_assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::TypeG::InteriorPointAlgorithm<double> > &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double>>::*)()) &Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::reset, "C++: Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::reset() --> void");
		cl.def("access_private_ptr", (class ROL::TypeG::InteriorPointAlgorithm<double> * (Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::access_private_ptr, "C++: Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::access_private_ptr() const --> class ROL::TypeG::InteriorPointAlgorithm<double> *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double>>::*)()) &Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::nonconst_access_private_node, "C++: Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::access_private_node, "C++: Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double> >::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double>> const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double>>> cl(M("Teuchos"), "RCP_ROL_TypeG_StabilizedLCLAlgorithm_double_t", "");
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double>>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("null_arg") );

		cl.def( pybind11::init( [](class ROL::TypeG::StabilizedLCLAlgorithm<double> * a0){ return new Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double>>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<class ROL::TypeG::StabilizedLCLAlgorithm<double> *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership_in") );

		cl.def( pybind11::init( [](Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double>> const &o){ return new Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double>>(o); } ) );
		cl.def( pybind11::init<class ROL::TypeG::StabilizedLCLAlgorithm<double> *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::TypeG::StabilizedLCLAlgorithm<double> *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::TypeG::StabilizedLCLAlgorithm<double> *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<class ROL::TypeG::StabilizedLCLAlgorithm<double> > & (Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double>>::*)(const class Teuchos::RCP<class ROL::TypeG::StabilizedLCLAlgorithm<double> > &)) &Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::operator=, "C++: Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::operator=(const class Teuchos::RCP<class ROL::TypeG::StabilizedLCLAlgorithm<double> > &) --> class Teuchos::RCP<class ROL::TypeG::StabilizedLCLAlgorithm<double> > &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<class ROL::TypeG::StabilizedLCLAlgorithm<double> > & (Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double>>::*)(enum Teuchos::ENull)) &Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::operator=, "C++: Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<class ROL::TypeG::StabilizedLCLAlgorithm<double> > &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double>>::*)(class Teuchos::RCP<class ROL::TypeG::StabilizedLCLAlgorithm<double> > &)) &Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::swap, "C++: Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::swap(class Teuchos::RCP<class ROL::TypeG::StabilizedLCLAlgorithm<double> > &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::is_null, "C++: Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::is_null() const --> bool");
		cl.def("arrow", (class ROL::TypeG::StabilizedLCLAlgorithm<double> * (Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::operator->, "C++: Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::operator->() const --> class ROL::TypeG::StabilizedLCLAlgorithm<double> *", pybind11::return_value_policy::automatic);
		cl.def("__mul__", (class ROL::TypeG::StabilizedLCLAlgorithm<double> & (Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::operator*, "C++: Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::operator*() const --> class ROL::TypeG::StabilizedLCLAlgorithm<double> &", pybind11::return_value_policy::automatic);
		cl.def("get", (class ROL::TypeG::StabilizedLCLAlgorithm<double> * (Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::get, "C++: Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::get() const --> class ROL::TypeG::StabilizedLCLAlgorithm<double> *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class ROL::TypeG::StabilizedLCLAlgorithm<double> * (Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::getRawPtr, "C++: Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::getRawPtr() const --> class ROL::TypeG::StabilizedLCLAlgorithm<double> *", pybind11::return_value_policy::automatic);
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::strength, "C++: Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::is_valid_ptr, "C++: Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::strong_count, "C++: Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::weak_count, "C++: Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::total_count, "C++: Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double>>::*)()) &Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::set_has_ownership, "C++: Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::has_ownership, "C++: Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<class ROL::TypeG::StabilizedLCLAlgorithm<double> > (Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::create_weak, "C++: Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::create_weak() const --> class Teuchos::RCP<class ROL::TypeG::StabilizedLCLAlgorithm<double> >");
		cl.def("create_strong", (class Teuchos::RCP<class ROL::TypeG::StabilizedLCLAlgorithm<double> > (Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::create_strong, "C++: Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::create_strong() const --> class Teuchos::RCP<class ROL::TypeG::StabilizedLCLAlgorithm<double> >");
		cl.def("assert_not_null", (const class Teuchos::RCP<class ROL::TypeG::StabilizedLCLAlgorithm<double> > & (Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::assert_not_null, "C++: Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::assert_not_null() const --> const class Teuchos::RCP<class ROL::TypeG::StabilizedLCLAlgorithm<double> > &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<class ROL::TypeG::StabilizedLCLAlgorithm<double> > & (Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::assert_valid_ptr, "C++: Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::TypeG::StabilizedLCLAlgorithm<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<class ROL::TypeG::StabilizedLCLAlgorithm<double> > & (Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::debug_assert_not_null, "C++: Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::debug_assert_not_null() const --> const class Teuchos::RCP<class ROL::TypeG::StabilizedLCLAlgorithm<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<class ROL::TypeG::StabilizedLCLAlgorithm<double> > & (Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::debug_assert_valid_ptr, "C++: Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::debug_assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::TypeG::StabilizedLCLAlgorithm<double> > &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double>>::*)()) &Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::reset, "C++: Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::reset() --> void");
		cl.def("access_private_ptr", (class ROL::TypeG::StabilizedLCLAlgorithm<double> * (Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::access_private_ptr, "C++: Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::access_private_ptr() const --> class ROL::TypeG::StabilizedLCLAlgorithm<double> *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double>>::*)()) &Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::nonconst_access_private_node, "C++: Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double>>::*)() const) &Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::access_private_node, "C++: Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double> >::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double>> const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<ROL::Constraint_Partitioned<double>>> cl(M("Teuchos"), "RCP_ROL_Constraint_Partitioned_double_t", "");
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<ROL::Constraint_Partitioned<double>>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("null_arg") );

		cl.def( pybind11::init( [](class ROL::Constraint_Partitioned<double> * a0){ return new Teuchos::RCP<ROL::Constraint_Partitioned<double>>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<class ROL::Constraint_Partitioned<double> *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership_in") );

		cl.def( pybind11::init( [](Teuchos::RCP<ROL::Constraint_Partitioned<double>> const &o){ return new Teuchos::RCP<ROL::Constraint_Partitioned<double>>(o); } ) );
		cl.def( pybind11::init<class ROL::Constraint_Partitioned<double> *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::Constraint_Partitioned<double> *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::Constraint_Partitioned<double> *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<class ROL::Constraint_Partitioned<double> > & (Teuchos::RCP<ROL::Constraint_Partitioned<double>>::*)(const class Teuchos::RCP<class ROL::Constraint_Partitioned<double> > &)) &Teuchos::RCP<ROL::Constraint_Partitioned<double> >::operator=, "C++: Teuchos::RCP<ROL::Constraint_Partitioned<double> >::operator=(const class Teuchos::RCP<class ROL::Constraint_Partitioned<double> > &) --> class Teuchos::RCP<class ROL::Constraint_Partitioned<double> > &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<class ROL::Constraint_Partitioned<double> > & (Teuchos::RCP<ROL::Constraint_Partitioned<double>>::*)(enum Teuchos::ENull)) &Teuchos::RCP<ROL::Constraint_Partitioned<double> >::operator=, "C++: Teuchos::RCP<ROL::Constraint_Partitioned<double> >::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<class ROL::Constraint_Partitioned<double> > &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<ROL::Constraint_Partitioned<double>>::*)(class Teuchos::RCP<class ROL::Constraint_Partitioned<double> > &)) &Teuchos::RCP<ROL::Constraint_Partitioned<double> >::swap, "C++: Teuchos::RCP<ROL::Constraint_Partitioned<double> >::swap(class Teuchos::RCP<class ROL::Constraint_Partitioned<double> > &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<ROL::Constraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::Constraint_Partitioned<double> >::is_null, "C++: Teuchos::RCP<ROL::Constraint_Partitioned<double> >::is_null() const --> bool");
		cl.def("arrow", (class ROL::Constraint_Partitioned<double> * (Teuchos::RCP<ROL::Constraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::Constraint_Partitioned<double> >::operator->, "C++: Teuchos::RCP<ROL::Constraint_Partitioned<double> >::operator->() const --> class ROL::Constraint_Partitioned<double> *", pybind11::return_value_policy::automatic);
		cl.def("__mul__", (class ROL::Constraint_Partitioned<double> & (Teuchos::RCP<ROL::Constraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::Constraint_Partitioned<double> >::operator*, "C++: Teuchos::RCP<ROL::Constraint_Partitioned<double> >::operator*() const --> class ROL::Constraint_Partitioned<double> &", pybind11::return_value_policy::automatic);
		cl.def("get", (class ROL::Constraint_Partitioned<double> * (Teuchos::RCP<ROL::Constraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::Constraint_Partitioned<double> >::get, "C++: Teuchos::RCP<ROL::Constraint_Partitioned<double> >::get() const --> class ROL::Constraint_Partitioned<double> *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class ROL::Constraint_Partitioned<double> * (Teuchos::RCP<ROL::Constraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::Constraint_Partitioned<double> >::getRawPtr, "C++: Teuchos::RCP<ROL::Constraint_Partitioned<double> >::getRawPtr() const --> class ROL::Constraint_Partitioned<double> *", pybind11::return_value_policy::automatic);
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<ROL::Constraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::Constraint_Partitioned<double> >::strength, "C++: Teuchos::RCP<ROL::Constraint_Partitioned<double> >::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<ROL::Constraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::Constraint_Partitioned<double> >::is_valid_ptr, "C++: Teuchos::RCP<ROL::Constraint_Partitioned<double> >::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<ROL::Constraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::Constraint_Partitioned<double> >::strong_count, "C++: Teuchos::RCP<ROL::Constraint_Partitioned<double> >::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<ROL::Constraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::Constraint_Partitioned<double> >::weak_count, "C++: Teuchos::RCP<ROL::Constraint_Partitioned<double> >::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<ROL::Constraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::Constraint_Partitioned<double> >::total_count, "C++: Teuchos::RCP<ROL::Constraint_Partitioned<double> >::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<ROL::Constraint_Partitioned<double>>::*)()) &Teuchos::RCP<ROL::Constraint_Partitioned<double> >::set_has_ownership, "C++: Teuchos::RCP<ROL::Constraint_Partitioned<double> >::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<ROL::Constraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::Constraint_Partitioned<double> >::has_ownership, "C++: Teuchos::RCP<ROL::Constraint_Partitioned<double> >::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<class ROL::Constraint_Partitioned<double> > (Teuchos::RCP<ROL::Constraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::Constraint_Partitioned<double> >::create_weak, "C++: Teuchos::RCP<ROL::Constraint_Partitioned<double> >::create_weak() const --> class Teuchos::RCP<class ROL::Constraint_Partitioned<double> >");
		cl.def("create_strong", (class Teuchos::RCP<class ROL::Constraint_Partitioned<double> > (Teuchos::RCP<ROL::Constraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::Constraint_Partitioned<double> >::create_strong, "C++: Teuchos::RCP<ROL::Constraint_Partitioned<double> >::create_strong() const --> class Teuchos::RCP<class ROL::Constraint_Partitioned<double> >");
		cl.def("assert_not_null", (const class Teuchos::RCP<class ROL::Constraint_Partitioned<double> > & (Teuchos::RCP<ROL::Constraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::Constraint_Partitioned<double> >::assert_not_null, "C++: Teuchos::RCP<ROL::Constraint_Partitioned<double> >::assert_not_null() const --> const class Teuchos::RCP<class ROL::Constraint_Partitioned<double> > &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<class ROL::Constraint_Partitioned<double> > & (Teuchos::RCP<ROL::Constraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::Constraint_Partitioned<double> >::assert_valid_ptr, "C++: Teuchos::RCP<ROL::Constraint_Partitioned<double> >::assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::Constraint_Partitioned<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<class ROL::Constraint_Partitioned<double> > & (Teuchos::RCP<ROL::Constraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::Constraint_Partitioned<double> >::debug_assert_not_null, "C++: Teuchos::RCP<ROL::Constraint_Partitioned<double> >::debug_assert_not_null() const --> const class Teuchos::RCP<class ROL::Constraint_Partitioned<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<class ROL::Constraint_Partitioned<double> > & (Teuchos::RCP<ROL::Constraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::Constraint_Partitioned<double> >::debug_assert_valid_ptr, "C++: Teuchos::RCP<ROL::Constraint_Partitioned<double> >::debug_assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::Constraint_Partitioned<double> > &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<ROL::Constraint_Partitioned<double>>::*)()) &Teuchos::RCP<ROL::Constraint_Partitioned<double> >::reset, "C++: Teuchos::RCP<ROL::Constraint_Partitioned<double> >::reset() --> void");
		cl.def("access_private_ptr", (class ROL::Constraint_Partitioned<double> * (Teuchos::RCP<ROL::Constraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::Constraint_Partitioned<double> >::access_private_ptr, "C++: Teuchos::RCP<ROL::Constraint_Partitioned<double> >::access_private_ptr() const --> class ROL::Constraint_Partitioned<double> *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::Constraint_Partitioned<double>>::*)()) &Teuchos::RCP<ROL::Constraint_Partitioned<double> >::nonconst_access_private_node, "C++: Teuchos::RCP<ROL::Constraint_Partitioned<double> >::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::Constraint_Partitioned<double>>::*)() const) &Teuchos::RCP<ROL::Constraint_Partitioned<double> >::access_private_node, "C++: Teuchos::RCP<ROL::Constraint_Partitioned<double> >::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<ROL::Constraint_Partitioned<double>> const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<ROL::AffineTransformObjective<double>>> cl(M("Teuchos"), "RCP_ROL_AffineTransformObjective_double_t", "");
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<ROL::AffineTransformObjective<double>>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("null_arg") );

		cl.def( pybind11::init( [](class ROL::AffineTransformObjective<double> * a0){ return new Teuchos::RCP<ROL::AffineTransformObjective<double>>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<class ROL::AffineTransformObjective<double> *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership_in") );

		cl.def( pybind11::init( [](Teuchos::RCP<ROL::AffineTransformObjective<double>> const &o){ return new Teuchos::RCP<ROL::AffineTransformObjective<double>>(o); } ) );
		cl.def( pybind11::init<class ROL::AffineTransformObjective<double> *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::AffineTransformObjective<double> *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::AffineTransformObjective<double> *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<class ROL::AffineTransformObjective<double> > & (Teuchos::RCP<ROL::AffineTransformObjective<double>>::*)(const class Teuchos::RCP<class ROL::AffineTransformObjective<double> > &)) &Teuchos::RCP<ROL::AffineTransformObjective<double> >::operator=, "C++: Teuchos::RCP<ROL::AffineTransformObjective<double> >::operator=(const class Teuchos::RCP<class ROL::AffineTransformObjective<double> > &) --> class Teuchos::RCP<class ROL::AffineTransformObjective<double> > &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<class ROL::AffineTransformObjective<double> > & (Teuchos::RCP<ROL::AffineTransformObjective<double>>::*)(enum Teuchos::ENull)) &Teuchos::RCP<ROL::AffineTransformObjective<double> >::operator=, "C++: Teuchos::RCP<ROL::AffineTransformObjective<double> >::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<class ROL::AffineTransformObjective<double> > &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<ROL::AffineTransformObjective<double>>::*)(class Teuchos::RCP<class ROL::AffineTransformObjective<double> > &)) &Teuchos::RCP<ROL::AffineTransformObjective<double> >::swap, "C++: Teuchos::RCP<ROL::AffineTransformObjective<double> >::swap(class Teuchos::RCP<class ROL::AffineTransformObjective<double> > &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<ROL::AffineTransformObjective<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformObjective<double> >::is_null, "C++: Teuchos::RCP<ROL::AffineTransformObjective<double> >::is_null() const --> bool");
		cl.def("arrow", (class ROL::AffineTransformObjective<double> * (Teuchos::RCP<ROL::AffineTransformObjective<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformObjective<double> >::operator->, "C++: Teuchos::RCP<ROL::AffineTransformObjective<double> >::operator->() const --> class ROL::AffineTransformObjective<double> *", pybind11::return_value_policy::automatic);
		cl.def("__mul__", (class ROL::AffineTransformObjective<double> & (Teuchos::RCP<ROL::AffineTransformObjective<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformObjective<double> >::operator*, "C++: Teuchos::RCP<ROL::AffineTransformObjective<double> >::operator*() const --> class ROL::AffineTransformObjective<double> &", pybind11::return_value_policy::automatic);
		cl.def("get", (class ROL::AffineTransformObjective<double> * (Teuchos::RCP<ROL::AffineTransformObjective<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformObjective<double> >::get, "C++: Teuchos::RCP<ROL::AffineTransformObjective<double> >::get() const --> class ROL::AffineTransformObjective<double> *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class ROL::AffineTransformObjective<double> * (Teuchos::RCP<ROL::AffineTransformObjective<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformObjective<double> >::getRawPtr, "C++: Teuchos::RCP<ROL::AffineTransformObjective<double> >::getRawPtr() const --> class ROL::AffineTransformObjective<double> *", pybind11::return_value_policy::automatic);
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<ROL::AffineTransformObjective<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformObjective<double> >::strength, "C++: Teuchos::RCP<ROL::AffineTransformObjective<double> >::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<ROL::AffineTransformObjective<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformObjective<double> >::is_valid_ptr, "C++: Teuchos::RCP<ROL::AffineTransformObjective<double> >::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<ROL::AffineTransformObjective<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformObjective<double> >::strong_count, "C++: Teuchos::RCP<ROL::AffineTransformObjective<double> >::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<ROL::AffineTransformObjective<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformObjective<double> >::weak_count, "C++: Teuchos::RCP<ROL::AffineTransformObjective<double> >::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<ROL::AffineTransformObjective<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformObjective<double> >::total_count, "C++: Teuchos::RCP<ROL::AffineTransformObjective<double> >::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<ROL::AffineTransformObjective<double>>::*)()) &Teuchos::RCP<ROL::AffineTransformObjective<double> >::set_has_ownership, "C++: Teuchos::RCP<ROL::AffineTransformObjective<double> >::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<ROL::AffineTransformObjective<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformObjective<double> >::has_ownership, "C++: Teuchos::RCP<ROL::AffineTransformObjective<double> >::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<class ROL::AffineTransformObjective<double> > (Teuchos::RCP<ROL::AffineTransformObjective<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformObjective<double> >::create_weak, "C++: Teuchos::RCP<ROL::AffineTransformObjective<double> >::create_weak() const --> class Teuchos::RCP<class ROL::AffineTransformObjective<double> >");
		cl.def("create_strong", (class Teuchos::RCP<class ROL::AffineTransformObjective<double> > (Teuchos::RCP<ROL::AffineTransformObjective<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformObjective<double> >::create_strong, "C++: Teuchos::RCP<ROL::AffineTransformObjective<double> >::create_strong() const --> class Teuchos::RCP<class ROL::AffineTransformObjective<double> >");
		cl.def("assert_not_null", (const class Teuchos::RCP<class ROL::AffineTransformObjective<double> > & (Teuchos::RCP<ROL::AffineTransformObjective<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformObjective<double> >::assert_not_null, "C++: Teuchos::RCP<ROL::AffineTransformObjective<double> >::assert_not_null() const --> const class Teuchos::RCP<class ROL::AffineTransformObjective<double> > &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<class ROL::AffineTransformObjective<double> > & (Teuchos::RCP<ROL::AffineTransformObjective<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformObjective<double> >::assert_valid_ptr, "C++: Teuchos::RCP<ROL::AffineTransformObjective<double> >::assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::AffineTransformObjective<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<class ROL::AffineTransformObjective<double> > & (Teuchos::RCP<ROL::AffineTransformObjective<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformObjective<double> >::debug_assert_not_null, "C++: Teuchos::RCP<ROL::AffineTransformObjective<double> >::debug_assert_not_null() const --> const class Teuchos::RCP<class ROL::AffineTransformObjective<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<class ROL::AffineTransformObjective<double> > & (Teuchos::RCP<ROL::AffineTransformObjective<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformObjective<double> >::debug_assert_valid_ptr, "C++: Teuchos::RCP<ROL::AffineTransformObjective<double> >::debug_assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::AffineTransformObjective<double> > &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<ROL::AffineTransformObjective<double>>::*)()) &Teuchos::RCP<ROL::AffineTransformObjective<double> >::reset, "C++: Teuchos::RCP<ROL::AffineTransformObjective<double> >::reset() --> void");
		cl.def("access_private_ptr", (class ROL::AffineTransformObjective<double> * (Teuchos::RCP<ROL::AffineTransformObjective<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformObjective<double> >::access_private_ptr, "C++: Teuchos::RCP<ROL::AffineTransformObjective<double> >::access_private_ptr() const --> class ROL::AffineTransformObjective<double> *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::AffineTransformObjective<double>>::*)()) &Teuchos::RCP<ROL::AffineTransformObjective<double> >::nonconst_access_private_node, "C++: Teuchos::RCP<ROL::AffineTransformObjective<double> >::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::AffineTransformObjective<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformObjective<double> >::access_private_node, "C++: Teuchos::RCP<ROL::AffineTransformObjective<double> >::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<ROL::AffineTransformObjective<double>> const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<ROL::LinearConstraint<double>>> cl(M("Teuchos"), "RCP_ROL_LinearConstraint_double_t", "");
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<ROL::LinearConstraint<double>>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("null_arg") );

		cl.def( pybind11::init( [](class ROL::LinearConstraint<double> * a0){ return new Teuchos::RCP<ROL::LinearConstraint<double>>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<class ROL::LinearConstraint<double> *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership_in") );

		cl.def( pybind11::init( [](Teuchos::RCP<ROL::LinearConstraint<double>> const &o){ return new Teuchos::RCP<ROL::LinearConstraint<double>>(o); } ) );
		cl.def( pybind11::init<class ROL::LinearConstraint<double> *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::LinearConstraint<double> *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::LinearConstraint<double> *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<class ROL::LinearConstraint<double> > & (Teuchos::RCP<ROL::LinearConstraint<double>>::*)(const class Teuchos::RCP<class ROL::LinearConstraint<double> > &)) &Teuchos::RCP<ROL::LinearConstraint<double> >::operator=, "C++: Teuchos::RCP<ROL::LinearConstraint<double> >::operator=(const class Teuchos::RCP<class ROL::LinearConstraint<double> > &) --> class Teuchos::RCP<class ROL::LinearConstraint<double> > &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<class ROL::LinearConstraint<double> > & (Teuchos::RCP<ROL::LinearConstraint<double>>::*)(enum Teuchos::ENull)) &Teuchos::RCP<ROL::LinearConstraint<double> >::operator=, "C++: Teuchos::RCP<ROL::LinearConstraint<double> >::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<class ROL::LinearConstraint<double> > &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<ROL::LinearConstraint<double>>::*)(class Teuchos::RCP<class ROL::LinearConstraint<double> > &)) &Teuchos::RCP<ROL::LinearConstraint<double> >::swap, "C++: Teuchos::RCP<ROL::LinearConstraint<double> >::swap(class Teuchos::RCP<class ROL::LinearConstraint<double> > &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<ROL::LinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::LinearConstraint<double> >::is_null, "C++: Teuchos::RCP<ROL::LinearConstraint<double> >::is_null() const --> bool");
		cl.def("arrow", (class ROL::LinearConstraint<double> * (Teuchos::RCP<ROL::LinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::LinearConstraint<double> >::operator->, "C++: Teuchos::RCP<ROL::LinearConstraint<double> >::operator->() const --> class ROL::LinearConstraint<double> *", pybind11::return_value_policy::automatic);
		cl.def("__mul__", (class ROL::LinearConstraint<double> & (Teuchos::RCP<ROL::LinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::LinearConstraint<double> >::operator*, "C++: Teuchos::RCP<ROL::LinearConstraint<double> >::operator*() const --> class ROL::LinearConstraint<double> &", pybind11::return_value_policy::automatic);
		cl.def("get", (class ROL::LinearConstraint<double> * (Teuchos::RCP<ROL::LinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::LinearConstraint<double> >::get, "C++: Teuchos::RCP<ROL::LinearConstraint<double> >::get() const --> class ROL::LinearConstraint<double> *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class ROL::LinearConstraint<double> * (Teuchos::RCP<ROL::LinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::LinearConstraint<double> >::getRawPtr, "C++: Teuchos::RCP<ROL::LinearConstraint<double> >::getRawPtr() const --> class ROL::LinearConstraint<double> *", pybind11::return_value_policy::automatic);
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<ROL::LinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::LinearConstraint<double> >::strength, "C++: Teuchos::RCP<ROL::LinearConstraint<double> >::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<ROL::LinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::LinearConstraint<double> >::is_valid_ptr, "C++: Teuchos::RCP<ROL::LinearConstraint<double> >::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<ROL::LinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::LinearConstraint<double> >::strong_count, "C++: Teuchos::RCP<ROL::LinearConstraint<double> >::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<ROL::LinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::LinearConstraint<double> >::weak_count, "C++: Teuchos::RCP<ROL::LinearConstraint<double> >::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<ROL::LinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::LinearConstraint<double> >::total_count, "C++: Teuchos::RCP<ROL::LinearConstraint<double> >::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<ROL::LinearConstraint<double>>::*)()) &Teuchos::RCP<ROL::LinearConstraint<double> >::set_has_ownership, "C++: Teuchos::RCP<ROL::LinearConstraint<double> >::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<ROL::LinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::LinearConstraint<double> >::has_ownership, "C++: Teuchos::RCP<ROL::LinearConstraint<double> >::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<class ROL::LinearConstraint<double> > (Teuchos::RCP<ROL::LinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::LinearConstraint<double> >::create_weak, "C++: Teuchos::RCP<ROL::LinearConstraint<double> >::create_weak() const --> class Teuchos::RCP<class ROL::LinearConstraint<double> >");
		cl.def("create_strong", (class Teuchos::RCP<class ROL::LinearConstraint<double> > (Teuchos::RCP<ROL::LinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::LinearConstraint<double> >::create_strong, "C++: Teuchos::RCP<ROL::LinearConstraint<double> >::create_strong() const --> class Teuchos::RCP<class ROL::LinearConstraint<double> >");
		cl.def("assert_not_null", (const class Teuchos::RCP<class ROL::LinearConstraint<double> > & (Teuchos::RCP<ROL::LinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::LinearConstraint<double> >::assert_not_null, "C++: Teuchos::RCP<ROL::LinearConstraint<double> >::assert_not_null() const --> const class Teuchos::RCP<class ROL::LinearConstraint<double> > &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<class ROL::LinearConstraint<double> > & (Teuchos::RCP<ROL::LinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::LinearConstraint<double> >::assert_valid_ptr, "C++: Teuchos::RCP<ROL::LinearConstraint<double> >::assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::LinearConstraint<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<class ROL::LinearConstraint<double> > & (Teuchos::RCP<ROL::LinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::LinearConstraint<double> >::debug_assert_not_null, "C++: Teuchos::RCP<ROL::LinearConstraint<double> >::debug_assert_not_null() const --> const class Teuchos::RCP<class ROL::LinearConstraint<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<class ROL::LinearConstraint<double> > & (Teuchos::RCP<ROL::LinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::LinearConstraint<double> >::debug_assert_valid_ptr, "C++: Teuchos::RCP<ROL::LinearConstraint<double> >::debug_assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::LinearConstraint<double> > &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<ROL::LinearConstraint<double>>::*)()) &Teuchos::RCP<ROL::LinearConstraint<double> >::reset, "C++: Teuchos::RCP<ROL::LinearConstraint<double> >::reset() --> void");
		cl.def("access_private_ptr", (class ROL::LinearConstraint<double> * (Teuchos::RCP<ROL::LinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::LinearConstraint<double> >::access_private_ptr, "C++: Teuchos::RCP<ROL::LinearConstraint<double> >::access_private_ptr() const --> class ROL::LinearConstraint<double> *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::LinearConstraint<double>>::*)()) &Teuchos::RCP<ROL::LinearConstraint<double> >::nonconst_access_private_node, "C++: Teuchos::RCP<ROL::LinearConstraint<double> >::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::LinearConstraint<double>>::*)() const) &Teuchos::RCP<ROL::LinearConstraint<double> >::access_private_node, "C++: Teuchos::RCP<ROL::LinearConstraint<double> >::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<ROL::LinearConstraint<double>> const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<ROL::AffineTransformConstraint<double>>> cl(M("Teuchos"), "RCP_ROL_AffineTransformConstraint_double_t", "");
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<ROL::AffineTransformConstraint<double>>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("null_arg") );

		cl.def( pybind11::init( [](class ROL::AffineTransformConstraint<double> * a0){ return new Teuchos::RCP<ROL::AffineTransformConstraint<double>>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<class ROL::AffineTransformConstraint<double> *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership_in") );

		cl.def( pybind11::init( [](Teuchos::RCP<ROL::AffineTransformConstraint<double>> const &o){ return new Teuchos::RCP<ROL::AffineTransformConstraint<double>>(o); } ) );
		cl.def( pybind11::init<class ROL::AffineTransformConstraint<double> *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::AffineTransformConstraint<double> *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::AffineTransformConstraint<double> *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<class ROL::AffineTransformConstraint<double> > & (Teuchos::RCP<ROL::AffineTransformConstraint<double>>::*)(const class Teuchos::RCP<class ROL::AffineTransformConstraint<double> > &)) &Teuchos::RCP<ROL::AffineTransformConstraint<double> >::operator=, "C++: Teuchos::RCP<ROL::AffineTransformConstraint<double> >::operator=(const class Teuchos::RCP<class ROL::AffineTransformConstraint<double> > &) --> class Teuchos::RCP<class ROL::AffineTransformConstraint<double> > &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<class ROL::AffineTransformConstraint<double> > & (Teuchos::RCP<ROL::AffineTransformConstraint<double>>::*)(enum Teuchos::ENull)) &Teuchos::RCP<ROL::AffineTransformConstraint<double> >::operator=, "C++: Teuchos::RCP<ROL::AffineTransformConstraint<double> >::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<class ROL::AffineTransformConstraint<double> > &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<ROL::AffineTransformConstraint<double>>::*)(class Teuchos::RCP<class ROL::AffineTransformConstraint<double> > &)) &Teuchos::RCP<ROL::AffineTransformConstraint<double> >::swap, "C++: Teuchos::RCP<ROL::AffineTransformConstraint<double> >::swap(class Teuchos::RCP<class ROL::AffineTransformConstraint<double> > &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<ROL::AffineTransformConstraint<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformConstraint<double> >::is_null, "C++: Teuchos::RCP<ROL::AffineTransformConstraint<double> >::is_null() const --> bool");
		cl.def("arrow", (class ROL::AffineTransformConstraint<double> * (Teuchos::RCP<ROL::AffineTransformConstraint<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformConstraint<double> >::operator->, "C++: Teuchos::RCP<ROL::AffineTransformConstraint<double> >::operator->() const --> class ROL::AffineTransformConstraint<double> *", pybind11::return_value_policy::automatic);
		cl.def("__mul__", (class ROL::AffineTransformConstraint<double> & (Teuchos::RCP<ROL::AffineTransformConstraint<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformConstraint<double> >::operator*, "C++: Teuchos::RCP<ROL::AffineTransformConstraint<double> >::operator*() const --> class ROL::AffineTransformConstraint<double> &", pybind11::return_value_policy::automatic);
		cl.def("get", (class ROL::AffineTransformConstraint<double> * (Teuchos::RCP<ROL::AffineTransformConstraint<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformConstraint<double> >::get, "C++: Teuchos::RCP<ROL::AffineTransformConstraint<double> >::get() const --> class ROL::AffineTransformConstraint<double> *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class ROL::AffineTransformConstraint<double> * (Teuchos::RCP<ROL::AffineTransformConstraint<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformConstraint<double> >::getRawPtr, "C++: Teuchos::RCP<ROL::AffineTransformConstraint<double> >::getRawPtr() const --> class ROL::AffineTransformConstraint<double> *", pybind11::return_value_policy::automatic);
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<ROL::AffineTransformConstraint<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformConstraint<double> >::strength, "C++: Teuchos::RCP<ROL::AffineTransformConstraint<double> >::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<ROL::AffineTransformConstraint<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformConstraint<double> >::is_valid_ptr, "C++: Teuchos::RCP<ROL::AffineTransformConstraint<double> >::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<ROL::AffineTransformConstraint<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformConstraint<double> >::strong_count, "C++: Teuchos::RCP<ROL::AffineTransformConstraint<double> >::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<ROL::AffineTransformConstraint<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformConstraint<double> >::weak_count, "C++: Teuchos::RCP<ROL::AffineTransformConstraint<double> >::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<ROL::AffineTransformConstraint<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformConstraint<double> >::total_count, "C++: Teuchos::RCP<ROL::AffineTransformConstraint<double> >::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<ROL::AffineTransformConstraint<double>>::*)()) &Teuchos::RCP<ROL::AffineTransformConstraint<double> >::set_has_ownership, "C++: Teuchos::RCP<ROL::AffineTransformConstraint<double> >::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<ROL::AffineTransformConstraint<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformConstraint<double> >::has_ownership, "C++: Teuchos::RCP<ROL::AffineTransformConstraint<double> >::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<class ROL::AffineTransformConstraint<double> > (Teuchos::RCP<ROL::AffineTransformConstraint<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformConstraint<double> >::create_weak, "C++: Teuchos::RCP<ROL::AffineTransformConstraint<double> >::create_weak() const --> class Teuchos::RCP<class ROL::AffineTransformConstraint<double> >");
		cl.def("create_strong", (class Teuchos::RCP<class ROL::AffineTransformConstraint<double> > (Teuchos::RCP<ROL::AffineTransformConstraint<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformConstraint<double> >::create_strong, "C++: Teuchos::RCP<ROL::AffineTransformConstraint<double> >::create_strong() const --> class Teuchos::RCP<class ROL::AffineTransformConstraint<double> >");
		cl.def("assert_not_null", (const class Teuchos::RCP<class ROL::AffineTransformConstraint<double> > & (Teuchos::RCP<ROL::AffineTransformConstraint<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformConstraint<double> >::assert_not_null, "C++: Teuchos::RCP<ROL::AffineTransformConstraint<double> >::assert_not_null() const --> const class Teuchos::RCP<class ROL::AffineTransformConstraint<double> > &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<class ROL::AffineTransformConstraint<double> > & (Teuchos::RCP<ROL::AffineTransformConstraint<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformConstraint<double> >::assert_valid_ptr, "C++: Teuchos::RCP<ROL::AffineTransformConstraint<double> >::assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::AffineTransformConstraint<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<class ROL::AffineTransformConstraint<double> > & (Teuchos::RCP<ROL::AffineTransformConstraint<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformConstraint<double> >::debug_assert_not_null, "C++: Teuchos::RCP<ROL::AffineTransformConstraint<double> >::debug_assert_not_null() const --> const class Teuchos::RCP<class ROL::AffineTransformConstraint<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<class ROL::AffineTransformConstraint<double> > & (Teuchos::RCP<ROL::AffineTransformConstraint<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformConstraint<double> >::debug_assert_valid_ptr, "C++: Teuchos::RCP<ROL::AffineTransformConstraint<double> >::debug_assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::AffineTransformConstraint<double> > &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<ROL::AffineTransformConstraint<double>>::*)()) &Teuchos::RCP<ROL::AffineTransformConstraint<double> >::reset, "C++: Teuchos::RCP<ROL::AffineTransformConstraint<double> >::reset() --> void");
		cl.def("access_private_ptr", (class ROL::AffineTransformConstraint<double> * (Teuchos::RCP<ROL::AffineTransformConstraint<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformConstraint<double> >::access_private_ptr, "C++: Teuchos::RCP<ROL::AffineTransformConstraint<double> >::access_private_ptr() const --> class ROL::AffineTransformConstraint<double> *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::AffineTransformConstraint<double>>::*)()) &Teuchos::RCP<ROL::AffineTransformConstraint<double> >::nonconst_access_private_node, "C++: Teuchos::RCP<ROL::AffineTransformConstraint<double> >::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::AffineTransformConstraint<double>>::*)() const) &Teuchos::RCP<ROL::AffineTransformConstraint<double> >::access_private_node, "C++: Teuchos::RCP<ROL::AffineTransformConstraint<double> >::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<ROL::AffineTransformConstraint<double>> const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<ROL::DaiFletcherProjection<double>>> cl(M("Teuchos"), "RCP_ROL_DaiFletcherProjection_double_t", "");
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<ROL::DaiFletcherProjection<double>>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("null_arg") );

		cl.def( pybind11::init( [](class ROL::DaiFletcherProjection<double> * a0){ return new Teuchos::RCP<ROL::DaiFletcherProjection<double>>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<class ROL::DaiFletcherProjection<double> *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership_in") );

		cl.def( pybind11::init( [](Teuchos::RCP<ROL::DaiFletcherProjection<double>> const &o){ return new Teuchos::RCP<ROL::DaiFletcherProjection<double>>(o); } ) );
		cl.def( pybind11::init<class ROL::DaiFletcherProjection<double> *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::DaiFletcherProjection<double> *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::DaiFletcherProjection<double> *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<class ROL::DaiFletcherProjection<double> > & (Teuchos::RCP<ROL::DaiFletcherProjection<double>>::*)(const class Teuchos::RCP<class ROL::DaiFletcherProjection<double> > &)) &Teuchos::RCP<ROL::DaiFletcherProjection<double> >::operator=, "C++: Teuchos::RCP<ROL::DaiFletcherProjection<double> >::operator=(const class Teuchos::RCP<class ROL::DaiFletcherProjection<double> > &) --> class Teuchos::RCP<class ROL::DaiFletcherProjection<double> > &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<class ROL::DaiFletcherProjection<double> > & (Teuchos::RCP<ROL::DaiFletcherProjection<double>>::*)(enum Teuchos::ENull)) &Teuchos::RCP<ROL::DaiFletcherProjection<double> >::operator=, "C++: Teuchos::RCP<ROL::DaiFletcherProjection<double> >::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<class ROL::DaiFletcherProjection<double> > &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<ROL::DaiFletcherProjection<double>>::*)(class Teuchos::RCP<class ROL::DaiFletcherProjection<double> > &)) &Teuchos::RCP<ROL::DaiFletcherProjection<double> >::swap, "C++: Teuchos::RCP<ROL::DaiFletcherProjection<double> >::swap(class Teuchos::RCP<class ROL::DaiFletcherProjection<double> > &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<ROL::DaiFletcherProjection<double>>::*)() const) &Teuchos::RCP<ROL::DaiFletcherProjection<double> >::is_null, "C++: Teuchos::RCP<ROL::DaiFletcherProjection<double> >::is_null() const --> bool");
		cl.def("arrow", (class ROL::DaiFletcherProjection<double> * (Teuchos::RCP<ROL::DaiFletcherProjection<double>>::*)() const) &Teuchos::RCP<ROL::DaiFletcherProjection<double> >::operator->, "C++: Teuchos::RCP<ROL::DaiFletcherProjection<double> >::operator->() const --> class ROL::DaiFletcherProjection<double> *", pybind11::return_value_policy::automatic);
		cl.def("__mul__", (class ROL::DaiFletcherProjection<double> & (Teuchos::RCP<ROL::DaiFletcherProjection<double>>::*)() const) &Teuchos::RCP<ROL::DaiFletcherProjection<double> >::operator*, "C++: Teuchos::RCP<ROL::DaiFletcherProjection<double> >::operator*() const --> class ROL::DaiFletcherProjection<double> &", pybind11::return_value_policy::automatic);
		cl.def("get", (class ROL::DaiFletcherProjection<double> * (Teuchos::RCP<ROL::DaiFletcherProjection<double>>::*)() const) &Teuchos::RCP<ROL::DaiFletcherProjection<double> >::get, "C++: Teuchos::RCP<ROL::DaiFletcherProjection<double> >::get() const --> class ROL::DaiFletcherProjection<double> *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class ROL::DaiFletcherProjection<double> * (Teuchos::RCP<ROL::DaiFletcherProjection<double>>::*)() const) &Teuchos::RCP<ROL::DaiFletcherProjection<double> >::getRawPtr, "C++: Teuchos::RCP<ROL::DaiFletcherProjection<double> >::getRawPtr() const --> class ROL::DaiFletcherProjection<double> *", pybind11::return_value_policy::automatic);
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<ROL::DaiFletcherProjection<double>>::*)() const) &Teuchos::RCP<ROL::DaiFletcherProjection<double> >::strength, "C++: Teuchos::RCP<ROL::DaiFletcherProjection<double> >::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<ROL::DaiFletcherProjection<double>>::*)() const) &Teuchos::RCP<ROL::DaiFletcherProjection<double> >::is_valid_ptr, "C++: Teuchos::RCP<ROL::DaiFletcherProjection<double> >::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<ROL::DaiFletcherProjection<double>>::*)() const) &Teuchos::RCP<ROL::DaiFletcherProjection<double> >::strong_count, "C++: Teuchos::RCP<ROL::DaiFletcherProjection<double> >::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<ROL::DaiFletcherProjection<double>>::*)() const) &Teuchos::RCP<ROL::DaiFletcherProjection<double> >::weak_count, "C++: Teuchos::RCP<ROL::DaiFletcherProjection<double> >::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<ROL::DaiFletcherProjection<double>>::*)() const) &Teuchos::RCP<ROL::DaiFletcherProjection<double> >::total_count, "C++: Teuchos::RCP<ROL::DaiFletcherProjection<double> >::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<ROL::DaiFletcherProjection<double>>::*)()) &Teuchos::RCP<ROL::DaiFletcherProjection<double> >::set_has_ownership, "C++: Teuchos::RCP<ROL::DaiFletcherProjection<double> >::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<ROL::DaiFletcherProjection<double>>::*)() const) &Teuchos::RCP<ROL::DaiFletcherProjection<double> >::has_ownership, "C++: Teuchos::RCP<ROL::DaiFletcherProjection<double> >::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<class ROL::DaiFletcherProjection<double> > (Teuchos::RCP<ROL::DaiFletcherProjection<double>>::*)() const) &Teuchos::RCP<ROL::DaiFletcherProjection<double> >::create_weak, "C++: Teuchos::RCP<ROL::DaiFletcherProjection<double> >::create_weak() const --> class Teuchos::RCP<class ROL::DaiFletcherProjection<double> >");
		cl.def("create_strong", (class Teuchos::RCP<class ROL::DaiFletcherProjection<double> > (Teuchos::RCP<ROL::DaiFletcherProjection<double>>::*)() const) &Teuchos::RCP<ROL::DaiFletcherProjection<double> >::create_strong, "C++: Teuchos::RCP<ROL::DaiFletcherProjection<double> >::create_strong() const --> class Teuchos::RCP<class ROL::DaiFletcherProjection<double> >");
		cl.def("assert_not_null", (const class Teuchos::RCP<class ROL::DaiFletcherProjection<double> > & (Teuchos::RCP<ROL::DaiFletcherProjection<double>>::*)() const) &Teuchos::RCP<ROL::DaiFletcherProjection<double> >::assert_not_null, "C++: Teuchos::RCP<ROL::DaiFletcherProjection<double> >::assert_not_null() const --> const class Teuchos::RCP<class ROL::DaiFletcherProjection<double> > &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<class ROL::DaiFletcherProjection<double> > & (Teuchos::RCP<ROL::DaiFletcherProjection<double>>::*)() const) &Teuchos::RCP<ROL::DaiFletcherProjection<double> >::assert_valid_ptr, "C++: Teuchos::RCP<ROL::DaiFletcherProjection<double> >::assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::DaiFletcherProjection<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<class ROL::DaiFletcherProjection<double> > & (Teuchos::RCP<ROL::DaiFletcherProjection<double>>::*)() const) &Teuchos::RCP<ROL::DaiFletcherProjection<double> >::debug_assert_not_null, "C++: Teuchos::RCP<ROL::DaiFletcherProjection<double> >::debug_assert_not_null() const --> const class Teuchos::RCP<class ROL::DaiFletcherProjection<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<class ROL::DaiFletcherProjection<double> > & (Teuchos::RCP<ROL::DaiFletcherProjection<double>>::*)() const) &Teuchos::RCP<ROL::DaiFletcherProjection<double> >::debug_assert_valid_ptr, "C++: Teuchos::RCP<ROL::DaiFletcherProjection<double> >::debug_assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::DaiFletcherProjection<double> > &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<ROL::DaiFletcherProjection<double>>::*)()) &Teuchos::RCP<ROL::DaiFletcherProjection<double> >::reset, "C++: Teuchos::RCP<ROL::DaiFletcherProjection<double> >::reset() --> void");
		cl.def("access_private_ptr", (class ROL::DaiFletcherProjection<double> * (Teuchos::RCP<ROL::DaiFletcherProjection<double>>::*)() const) &Teuchos::RCP<ROL::DaiFletcherProjection<double> >::access_private_ptr, "C++: Teuchos::RCP<ROL::DaiFletcherProjection<double> >::access_private_ptr() const --> class ROL::DaiFletcherProjection<double> *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::DaiFletcherProjection<double>>::*)()) &Teuchos::RCP<ROL::DaiFletcherProjection<double> >::nonconst_access_private_node, "C++: Teuchos::RCP<ROL::DaiFletcherProjection<double> >::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::DaiFletcherProjection<double>>::*)() const) &Teuchos::RCP<ROL::DaiFletcherProjection<double> >::access_private_node, "C++: Teuchos::RCP<ROL::DaiFletcherProjection<double> >::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<ROL::DaiFletcherProjection<double>> const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<ROL::DykstraProjection<double>>> cl(M("Teuchos"), "RCP_ROL_DykstraProjection_double_t", "");
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<ROL::DykstraProjection<double>>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("null_arg") );

		cl.def( pybind11::init( [](class ROL::DykstraProjection<double> * a0){ return new Teuchos::RCP<ROL::DykstraProjection<double>>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<class ROL::DykstraProjection<double> *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership_in") );

		cl.def( pybind11::init( [](Teuchos::RCP<ROL::DykstraProjection<double>> const &o){ return new Teuchos::RCP<ROL::DykstraProjection<double>>(o); } ) );
		cl.def( pybind11::init<class ROL::DykstraProjection<double> *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::DykstraProjection<double> *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::DykstraProjection<double> *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<class ROL::DykstraProjection<double> > & (Teuchos::RCP<ROL::DykstraProjection<double>>::*)(const class Teuchos::RCP<class ROL::DykstraProjection<double> > &)) &Teuchos::RCP<ROL::DykstraProjection<double> >::operator=, "C++: Teuchos::RCP<ROL::DykstraProjection<double> >::operator=(const class Teuchos::RCP<class ROL::DykstraProjection<double> > &) --> class Teuchos::RCP<class ROL::DykstraProjection<double> > &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<class ROL::DykstraProjection<double> > & (Teuchos::RCP<ROL::DykstraProjection<double>>::*)(enum Teuchos::ENull)) &Teuchos::RCP<ROL::DykstraProjection<double> >::operator=, "C++: Teuchos::RCP<ROL::DykstraProjection<double> >::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<class ROL::DykstraProjection<double> > &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<ROL::DykstraProjection<double>>::*)(class Teuchos::RCP<class ROL::DykstraProjection<double> > &)) &Teuchos::RCP<ROL::DykstraProjection<double> >::swap, "C++: Teuchos::RCP<ROL::DykstraProjection<double> >::swap(class Teuchos::RCP<class ROL::DykstraProjection<double> > &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<ROL::DykstraProjection<double>>::*)() const) &Teuchos::RCP<ROL::DykstraProjection<double> >::is_null, "C++: Teuchos::RCP<ROL::DykstraProjection<double> >::is_null() const --> bool");
		cl.def("arrow", (class ROL::DykstraProjection<double> * (Teuchos::RCP<ROL::DykstraProjection<double>>::*)() const) &Teuchos::RCP<ROL::DykstraProjection<double> >::operator->, "C++: Teuchos::RCP<ROL::DykstraProjection<double> >::operator->() const --> class ROL::DykstraProjection<double> *", pybind11::return_value_policy::automatic);
		cl.def("__mul__", (class ROL::DykstraProjection<double> & (Teuchos::RCP<ROL::DykstraProjection<double>>::*)() const) &Teuchos::RCP<ROL::DykstraProjection<double> >::operator*, "C++: Teuchos::RCP<ROL::DykstraProjection<double> >::operator*() const --> class ROL::DykstraProjection<double> &", pybind11::return_value_policy::automatic);
		cl.def("get", (class ROL::DykstraProjection<double> * (Teuchos::RCP<ROL::DykstraProjection<double>>::*)() const) &Teuchos::RCP<ROL::DykstraProjection<double> >::get, "C++: Teuchos::RCP<ROL::DykstraProjection<double> >::get() const --> class ROL::DykstraProjection<double> *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class ROL::DykstraProjection<double> * (Teuchos::RCP<ROL::DykstraProjection<double>>::*)() const) &Teuchos::RCP<ROL::DykstraProjection<double> >::getRawPtr, "C++: Teuchos::RCP<ROL::DykstraProjection<double> >::getRawPtr() const --> class ROL::DykstraProjection<double> *", pybind11::return_value_policy::automatic);
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<ROL::DykstraProjection<double>>::*)() const) &Teuchos::RCP<ROL::DykstraProjection<double> >::strength, "C++: Teuchos::RCP<ROL::DykstraProjection<double> >::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<ROL::DykstraProjection<double>>::*)() const) &Teuchos::RCP<ROL::DykstraProjection<double> >::is_valid_ptr, "C++: Teuchos::RCP<ROL::DykstraProjection<double> >::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<ROL::DykstraProjection<double>>::*)() const) &Teuchos::RCP<ROL::DykstraProjection<double> >::strong_count, "C++: Teuchos::RCP<ROL::DykstraProjection<double> >::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<ROL::DykstraProjection<double>>::*)() const) &Teuchos::RCP<ROL::DykstraProjection<double> >::weak_count, "C++: Teuchos::RCP<ROL::DykstraProjection<double> >::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<ROL::DykstraProjection<double>>::*)() const) &Teuchos::RCP<ROL::DykstraProjection<double> >::total_count, "C++: Teuchos::RCP<ROL::DykstraProjection<double> >::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<ROL::DykstraProjection<double>>::*)()) &Teuchos::RCP<ROL::DykstraProjection<double> >::set_has_ownership, "C++: Teuchos::RCP<ROL::DykstraProjection<double> >::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<ROL::DykstraProjection<double>>::*)() const) &Teuchos::RCP<ROL::DykstraProjection<double> >::has_ownership, "C++: Teuchos::RCP<ROL::DykstraProjection<double> >::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<class ROL::DykstraProjection<double> > (Teuchos::RCP<ROL::DykstraProjection<double>>::*)() const) &Teuchos::RCP<ROL::DykstraProjection<double> >::create_weak, "C++: Teuchos::RCP<ROL::DykstraProjection<double> >::create_weak() const --> class Teuchos::RCP<class ROL::DykstraProjection<double> >");
		cl.def("create_strong", (class Teuchos::RCP<class ROL::DykstraProjection<double> > (Teuchos::RCP<ROL::DykstraProjection<double>>::*)() const) &Teuchos::RCP<ROL::DykstraProjection<double> >::create_strong, "C++: Teuchos::RCP<ROL::DykstraProjection<double> >::create_strong() const --> class Teuchos::RCP<class ROL::DykstraProjection<double> >");
		cl.def("assert_not_null", (const class Teuchos::RCP<class ROL::DykstraProjection<double> > & (Teuchos::RCP<ROL::DykstraProjection<double>>::*)() const) &Teuchos::RCP<ROL::DykstraProjection<double> >::assert_not_null, "C++: Teuchos::RCP<ROL::DykstraProjection<double> >::assert_not_null() const --> const class Teuchos::RCP<class ROL::DykstraProjection<double> > &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<class ROL::DykstraProjection<double> > & (Teuchos::RCP<ROL::DykstraProjection<double>>::*)() const) &Teuchos::RCP<ROL::DykstraProjection<double> >::assert_valid_ptr, "C++: Teuchos::RCP<ROL::DykstraProjection<double> >::assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::DykstraProjection<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<class ROL::DykstraProjection<double> > & (Teuchos::RCP<ROL::DykstraProjection<double>>::*)() const) &Teuchos::RCP<ROL::DykstraProjection<double> >::debug_assert_not_null, "C++: Teuchos::RCP<ROL::DykstraProjection<double> >::debug_assert_not_null() const --> const class Teuchos::RCP<class ROL::DykstraProjection<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<class ROL::DykstraProjection<double> > & (Teuchos::RCP<ROL::DykstraProjection<double>>::*)() const) &Teuchos::RCP<ROL::DykstraProjection<double> >::debug_assert_valid_ptr, "C++: Teuchos::RCP<ROL::DykstraProjection<double> >::debug_assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::DykstraProjection<double> > &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<ROL::DykstraProjection<double>>::*)()) &Teuchos::RCP<ROL::DykstraProjection<double> >::reset, "C++: Teuchos::RCP<ROL::DykstraProjection<double> >::reset() --> void");
		cl.def("access_private_ptr", (class ROL::DykstraProjection<double> * (Teuchos::RCP<ROL::DykstraProjection<double>>::*)() const) &Teuchos::RCP<ROL::DykstraProjection<double> >::access_private_ptr, "C++: Teuchos::RCP<ROL::DykstraProjection<double> >::access_private_ptr() const --> class ROL::DykstraProjection<double> *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::DykstraProjection<double>>::*)()) &Teuchos::RCP<ROL::DykstraProjection<double> >::nonconst_access_private_node, "C++: Teuchos::RCP<ROL::DykstraProjection<double> >::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::DykstraProjection<double>>::*)() const) &Teuchos::RCP<ROL::DykstraProjection<double> >::access_private_node, "C++: Teuchos::RCP<ROL::DykstraProjection<double> >::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<ROL::DykstraProjection<double>> const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<ROL::DouglasRachfordProjection<double>>> cl(M("Teuchos"), "RCP_ROL_DouglasRachfordProjection_double_t", "");
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<ROL::DouglasRachfordProjection<double>>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("null_arg") );

		cl.def( pybind11::init( [](class ROL::DouglasRachfordProjection<double> * a0){ return new Teuchos::RCP<ROL::DouglasRachfordProjection<double>>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<class ROL::DouglasRachfordProjection<double> *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership_in") );

		cl.def( pybind11::init( [](Teuchos::RCP<ROL::DouglasRachfordProjection<double>> const &o){ return new Teuchos::RCP<ROL::DouglasRachfordProjection<double>>(o); } ) );
		cl.def( pybind11::init<class ROL::DouglasRachfordProjection<double> *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::DouglasRachfordProjection<double> *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::DouglasRachfordProjection<double> *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<class ROL::DouglasRachfordProjection<double> > & (Teuchos::RCP<ROL::DouglasRachfordProjection<double>>::*)(const class Teuchos::RCP<class ROL::DouglasRachfordProjection<double> > &)) &Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::operator=, "C++: Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::operator=(const class Teuchos::RCP<class ROL::DouglasRachfordProjection<double> > &) --> class Teuchos::RCP<class ROL::DouglasRachfordProjection<double> > &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<class ROL::DouglasRachfordProjection<double> > & (Teuchos::RCP<ROL::DouglasRachfordProjection<double>>::*)(enum Teuchos::ENull)) &Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::operator=, "C++: Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<class ROL::DouglasRachfordProjection<double> > &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<ROL::DouglasRachfordProjection<double>>::*)(class Teuchos::RCP<class ROL::DouglasRachfordProjection<double> > &)) &Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::swap, "C++: Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::swap(class Teuchos::RCP<class ROL::DouglasRachfordProjection<double> > &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<ROL::DouglasRachfordProjection<double>>::*)() const) &Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::is_null, "C++: Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::is_null() const --> bool");
		cl.def("arrow", (class ROL::DouglasRachfordProjection<double> * (Teuchos::RCP<ROL::DouglasRachfordProjection<double>>::*)() const) &Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::operator->, "C++: Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::operator->() const --> class ROL::DouglasRachfordProjection<double> *", pybind11::return_value_policy::automatic);
		cl.def("__mul__", (class ROL::DouglasRachfordProjection<double> & (Teuchos::RCP<ROL::DouglasRachfordProjection<double>>::*)() const) &Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::operator*, "C++: Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::operator*() const --> class ROL::DouglasRachfordProjection<double> &", pybind11::return_value_policy::automatic);
		cl.def("get", (class ROL::DouglasRachfordProjection<double> * (Teuchos::RCP<ROL::DouglasRachfordProjection<double>>::*)() const) &Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::get, "C++: Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::get() const --> class ROL::DouglasRachfordProjection<double> *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class ROL::DouglasRachfordProjection<double> * (Teuchos::RCP<ROL::DouglasRachfordProjection<double>>::*)() const) &Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::getRawPtr, "C++: Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::getRawPtr() const --> class ROL::DouglasRachfordProjection<double> *", pybind11::return_value_policy::automatic);
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<ROL::DouglasRachfordProjection<double>>::*)() const) &Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::strength, "C++: Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<ROL::DouglasRachfordProjection<double>>::*)() const) &Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::is_valid_ptr, "C++: Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<ROL::DouglasRachfordProjection<double>>::*)() const) &Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::strong_count, "C++: Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<ROL::DouglasRachfordProjection<double>>::*)() const) &Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::weak_count, "C++: Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<ROL::DouglasRachfordProjection<double>>::*)() const) &Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::total_count, "C++: Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<ROL::DouglasRachfordProjection<double>>::*)()) &Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::set_has_ownership, "C++: Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<ROL::DouglasRachfordProjection<double>>::*)() const) &Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::has_ownership, "C++: Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<class ROL::DouglasRachfordProjection<double> > (Teuchos::RCP<ROL::DouglasRachfordProjection<double>>::*)() const) &Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::create_weak, "C++: Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::create_weak() const --> class Teuchos::RCP<class ROL::DouglasRachfordProjection<double> >");
		cl.def("create_strong", (class Teuchos::RCP<class ROL::DouglasRachfordProjection<double> > (Teuchos::RCP<ROL::DouglasRachfordProjection<double>>::*)() const) &Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::create_strong, "C++: Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::create_strong() const --> class Teuchos::RCP<class ROL::DouglasRachfordProjection<double> >");
		cl.def("assert_not_null", (const class Teuchos::RCP<class ROL::DouglasRachfordProjection<double> > & (Teuchos::RCP<ROL::DouglasRachfordProjection<double>>::*)() const) &Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::assert_not_null, "C++: Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::assert_not_null() const --> const class Teuchos::RCP<class ROL::DouglasRachfordProjection<double> > &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<class ROL::DouglasRachfordProjection<double> > & (Teuchos::RCP<ROL::DouglasRachfordProjection<double>>::*)() const) &Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::assert_valid_ptr, "C++: Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::DouglasRachfordProjection<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<class ROL::DouglasRachfordProjection<double> > & (Teuchos::RCP<ROL::DouglasRachfordProjection<double>>::*)() const) &Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::debug_assert_not_null, "C++: Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::debug_assert_not_null() const --> const class Teuchos::RCP<class ROL::DouglasRachfordProjection<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<class ROL::DouglasRachfordProjection<double> > & (Teuchos::RCP<ROL::DouglasRachfordProjection<double>>::*)() const) &Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::debug_assert_valid_ptr, "C++: Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::debug_assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::DouglasRachfordProjection<double> > &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<ROL::DouglasRachfordProjection<double>>::*)()) &Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::reset, "C++: Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::reset() --> void");
		cl.def("access_private_ptr", (class ROL::DouglasRachfordProjection<double> * (Teuchos::RCP<ROL::DouglasRachfordProjection<double>>::*)() const) &Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::access_private_ptr, "C++: Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::access_private_ptr() const --> class ROL::DouglasRachfordProjection<double> *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::DouglasRachfordProjection<double>>::*)()) &Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::nonconst_access_private_node, "C++: Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::DouglasRachfordProjection<double>>::*)() const) &Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::access_private_node, "C++: Teuchos::RCP<ROL::DouglasRachfordProjection<double> >::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<ROL::DouglasRachfordProjection<double>> const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<ROL::SemismoothNewtonProjection<double>>> cl(M("Teuchos"), "RCP_ROL_SemismoothNewtonProjection_double_t", "");
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<ROL::SemismoothNewtonProjection<double>>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("null_arg") );

		cl.def( pybind11::init( [](class ROL::SemismoothNewtonProjection<double> * a0){ return new Teuchos::RCP<ROL::SemismoothNewtonProjection<double>>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<class ROL::SemismoothNewtonProjection<double> *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership_in") );

		cl.def( pybind11::init( [](Teuchos::RCP<ROL::SemismoothNewtonProjection<double>> const &o){ return new Teuchos::RCP<ROL::SemismoothNewtonProjection<double>>(o); } ) );
		cl.def( pybind11::init<class ROL::SemismoothNewtonProjection<double> *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::SemismoothNewtonProjection<double> *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::SemismoothNewtonProjection<double> *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<class ROL::SemismoothNewtonProjection<double> > & (Teuchos::RCP<ROL::SemismoothNewtonProjection<double>>::*)(const class Teuchos::RCP<class ROL::SemismoothNewtonProjection<double> > &)) &Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::operator=, "C++: Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::operator=(const class Teuchos::RCP<class ROL::SemismoothNewtonProjection<double> > &) --> class Teuchos::RCP<class ROL::SemismoothNewtonProjection<double> > &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<class ROL::SemismoothNewtonProjection<double> > & (Teuchos::RCP<ROL::SemismoothNewtonProjection<double>>::*)(enum Teuchos::ENull)) &Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::operator=, "C++: Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<class ROL::SemismoothNewtonProjection<double> > &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<ROL::SemismoothNewtonProjection<double>>::*)(class Teuchos::RCP<class ROL::SemismoothNewtonProjection<double> > &)) &Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::swap, "C++: Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::swap(class Teuchos::RCP<class ROL::SemismoothNewtonProjection<double> > &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<ROL::SemismoothNewtonProjection<double>>::*)() const) &Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::is_null, "C++: Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::is_null() const --> bool");
		cl.def("arrow", (class ROL::SemismoothNewtonProjection<double> * (Teuchos::RCP<ROL::SemismoothNewtonProjection<double>>::*)() const) &Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::operator->, "C++: Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::operator->() const --> class ROL::SemismoothNewtonProjection<double> *", pybind11::return_value_policy::automatic);
		cl.def("__mul__", (class ROL::SemismoothNewtonProjection<double> & (Teuchos::RCP<ROL::SemismoothNewtonProjection<double>>::*)() const) &Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::operator*, "C++: Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::operator*() const --> class ROL::SemismoothNewtonProjection<double> &", pybind11::return_value_policy::automatic);
		cl.def("get", (class ROL::SemismoothNewtonProjection<double> * (Teuchos::RCP<ROL::SemismoothNewtonProjection<double>>::*)() const) &Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::get, "C++: Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::get() const --> class ROL::SemismoothNewtonProjection<double> *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class ROL::SemismoothNewtonProjection<double> * (Teuchos::RCP<ROL::SemismoothNewtonProjection<double>>::*)() const) &Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::getRawPtr, "C++: Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::getRawPtr() const --> class ROL::SemismoothNewtonProjection<double> *", pybind11::return_value_policy::automatic);
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<ROL::SemismoothNewtonProjection<double>>::*)() const) &Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::strength, "C++: Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<ROL::SemismoothNewtonProjection<double>>::*)() const) &Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::is_valid_ptr, "C++: Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<ROL::SemismoothNewtonProjection<double>>::*)() const) &Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::strong_count, "C++: Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<ROL::SemismoothNewtonProjection<double>>::*)() const) &Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::weak_count, "C++: Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<ROL::SemismoothNewtonProjection<double>>::*)() const) &Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::total_count, "C++: Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<ROL::SemismoothNewtonProjection<double>>::*)()) &Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::set_has_ownership, "C++: Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<ROL::SemismoothNewtonProjection<double>>::*)() const) &Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::has_ownership, "C++: Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<class ROL::SemismoothNewtonProjection<double> > (Teuchos::RCP<ROL::SemismoothNewtonProjection<double>>::*)() const) &Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::create_weak, "C++: Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::create_weak() const --> class Teuchos::RCP<class ROL::SemismoothNewtonProjection<double> >");
		cl.def("create_strong", (class Teuchos::RCP<class ROL::SemismoothNewtonProjection<double> > (Teuchos::RCP<ROL::SemismoothNewtonProjection<double>>::*)() const) &Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::create_strong, "C++: Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::create_strong() const --> class Teuchos::RCP<class ROL::SemismoothNewtonProjection<double> >");
		cl.def("assert_not_null", (const class Teuchos::RCP<class ROL::SemismoothNewtonProjection<double> > & (Teuchos::RCP<ROL::SemismoothNewtonProjection<double>>::*)() const) &Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::assert_not_null, "C++: Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::assert_not_null() const --> const class Teuchos::RCP<class ROL::SemismoothNewtonProjection<double> > &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<class ROL::SemismoothNewtonProjection<double> > & (Teuchos::RCP<ROL::SemismoothNewtonProjection<double>>::*)() const) &Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::assert_valid_ptr, "C++: Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::SemismoothNewtonProjection<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<class ROL::SemismoothNewtonProjection<double> > & (Teuchos::RCP<ROL::SemismoothNewtonProjection<double>>::*)() const) &Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::debug_assert_not_null, "C++: Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::debug_assert_not_null() const --> const class Teuchos::RCP<class ROL::SemismoothNewtonProjection<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<class ROL::SemismoothNewtonProjection<double> > & (Teuchos::RCP<ROL::SemismoothNewtonProjection<double>>::*)() const) &Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::debug_assert_valid_ptr, "C++: Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::debug_assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::SemismoothNewtonProjection<double> > &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<ROL::SemismoothNewtonProjection<double>>::*)()) &Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::reset, "C++: Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::reset() --> void");
		cl.def("access_private_ptr", (class ROL::SemismoothNewtonProjection<double> * (Teuchos::RCP<ROL::SemismoothNewtonProjection<double>>::*)() const) &Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::access_private_ptr, "C++: Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::access_private_ptr() const --> class ROL::SemismoothNewtonProjection<double> *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::SemismoothNewtonProjection<double>>::*)()) &Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::nonconst_access_private_node, "C++: Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::SemismoothNewtonProjection<double>>::*)() const) &Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::access_private_node, "C++: Teuchos::RCP<ROL::SemismoothNewtonProjection<double> >::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<ROL::SemismoothNewtonProjection<double>> const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<ROL::RiddersProjection<double>>> cl(M("Teuchos"), "RCP_ROL_RiddersProjection_double_t", "");
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<ROL::RiddersProjection<double>>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("null_arg") );

		cl.def( pybind11::init( [](class ROL::RiddersProjection<double> * a0){ return new Teuchos::RCP<ROL::RiddersProjection<double>>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<class ROL::RiddersProjection<double> *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership_in") );

		cl.def( pybind11::init( [](Teuchos::RCP<ROL::RiddersProjection<double>> const &o){ return new Teuchos::RCP<ROL::RiddersProjection<double>>(o); } ) );
		cl.def( pybind11::init<class ROL::RiddersProjection<double> *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::RiddersProjection<double> *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::RiddersProjection<double> *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<class ROL::RiddersProjection<double> > & (Teuchos::RCP<ROL::RiddersProjection<double>>::*)(const class Teuchos::RCP<class ROL::RiddersProjection<double> > &)) &Teuchos::RCP<ROL::RiddersProjection<double> >::operator=, "C++: Teuchos::RCP<ROL::RiddersProjection<double> >::operator=(const class Teuchos::RCP<class ROL::RiddersProjection<double> > &) --> class Teuchos::RCP<class ROL::RiddersProjection<double> > &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<class ROL::RiddersProjection<double> > & (Teuchos::RCP<ROL::RiddersProjection<double>>::*)(enum Teuchos::ENull)) &Teuchos::RCP<ROL::RiddersProjection<double> >::operator=, "C++: Teuchos::RCP<ROL::RiddersProjection<double> >::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<class ROL::RiddersProjection<double> > &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<ROL::RiddersProjection<double>>::*)(class Teuchos::RCP<class ROL::RiddersProjection<double> > &)) &Teuchos::RCP<ROL::RiddersProjection<double> >::swap, "C++: Teuchos::RCP<ROL::RiddersProjection<double> >::swap(class Teuchos::RCP<class ROL::RiddersProjection<double> > &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<ROL::RiddersProjection<double>>::*)() const) &Teuchos::RCP<ROL::RiddersProjection<double> >::is_null, "C++: Teuchos::RCP<ROL::RiddersProjection<double> >::is_null() const --> bool");
		cl.def("arrow", (class ROL::RiddersProjection<double> * (Teuchos::RCP<ROL::RiddersProjection<double>>::*)() const) &Teuchos::RCP<ROL::RiddersProjection<double> >::operator->, "C++: Teuchos::RCP<ROL::RiddersProjection<double> >::operator->() const --> class ROL::RiddersProjection<double> *", pybind11::return_value_policy::automatic);
		cl.def("__mul__", (class ROL::RiddersProjection<double> & (Teuchos::RCP<ROL::RiddersProjection<double>>::*)() const) &Teuchos::RCP<ROL::RiddersProjection<double> >::operator*, "C++: Teuchos::RCP<ROL::RiddersProjection<double> >::operator*() const --> class ROL::RiddersProjection<double> &", pybind11::return_value_policy::automatic);
		cl.def("get", (class ROL::RiddersProjection<double> * (Teuchos::RCP<ROL::RiddersProjection<double>>::*)() const) &Teuchos::RCP<ROL::RiddersProjection<double> >::get, "C++: Teuchos::RCP<ROL::RiddersProjection<double> >::get() const --> class ROL::RiddersProjection<double> *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class ROL::RiddersProjection<double> * (Teuchos::RCP<ROL::RiddersProjection<double>>::*)() const) &Teuchos::RCP<ROL::RiddersProjection<double> >::getRawPtr, "C++: Teuchos::RCP<ROL::RiddersProjection<double> >::getRawPtr() const --> class ROL::RiddersProjection<double> *", pybind11::return_value_policy::automatic);
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<ROL::RiddersProjection<double>>::*)() const) &Teuchos::RCP<ROL::RiddersProjection<double> >::strength, "C++: Teuchos::RCP<ROL::RiddersProjection<double> >::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<ROL::RiddersProjection<double>>::*)() const) &Teuchos::RCP<ROL::RiddersProjection<double> >::is_valid_ptr, "C++: Teuchos::RCP<ROL::RiddersProjection<double> >::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<ROL::RiddersProjection<double>>::*)() const) &Teuchos::RCP<ROL::RiddersProjection<double> >::strong_count, "C++: Teuchos::RCP<ROL::RiddersProjection<double> >::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<ROL::RiddersProjection<double>>::*)() const) &Teuchos::RCP<ROL::RiddersProjection<double> >::weak_count, "C++: Teuchos::RCP<ROL::RiddersProjection<double> >::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<ROL::RiddersProjection<double>>::*)() const) &Teuchos::RCP<ROL::RiddersProjection<double> >::total_count, "C++: Teuchos::RCP<ROL::RiddersProjection<double> >::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<ROL::RiddersProjection<double>>::*)()) &Teuchos::RCP<ROL::RiddersProjection<double> >::set_has_ownership, "C++: Teuchos::RCP<ROL::RiddersProjection<double> >::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<ROL::RiddersProjection<double>>::*)() const) &Teuchos::RCP<ROL::RiddersProjection<double> >::has_ownership, "C++: Teuchos::RCP<ROL::RiddersProjection<double> >::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<class ROL::RiddersProjection<double> > (Teuchos::RCP<ROL::RiddersProjection<double>>::*)() const) &Teuchos::RCP<ROL::RiddersProjection<double> >::create_weak, "C++: Teuchos::RCP<ROL::RiddersProjection<double> >::create_weak() const --> class Teuchos::RCP<class ROL::RiddersProjection<double> >");
		cl.def("create_strong", (class Teuchos::RCP<class ROL::RiddersProjection<double> > (Teuchos::RCP<ROL::RiddersProjection<double>>::*)() const) &Teuchos::RCP<ROL::RiddersProjection<double> >::create_strong, "C++: Teuchos::RCP<ROL::RiddersProjection<double> >::create_strong() const --> class Teuchos::RCP<class ROL::RiddersProjection<double> >");
		cl.def("assert_not_null", (const class Teuchos::RCP<class ROL::RiddersProjection<double> > & (Teuchos::RCP<ROL::RiddersProjection<double>>::*)() const) &Teuchos::RCP<ROL::RiddersProjection<double> >::assert_not_null, "C++: Teuchos::RCP<ROL::RiddersProjection<double> >::assert_not_null() const --> const class Teuchos::RCP<class ROL::RiddersProjection<double> > &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<class ROL::RiddersProjection<double> > & (Teuchos::RCP<ROL::RiddersProjection<double>>::*)() const) &Teuchos::RCP<ROL::RiddersProjection<double> >::assert_valid_ptr, "C++: Teuchos::RCP<ROL::RiddersProjection<double> >::assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::RiddersProjection<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<class ROL::RiddersProjection<double> > & (Teuchos::RCP<ROL::RiddersProjection<double>>::*)() const) &Teuchos::RCP<ROL::RiddersProjection<double> >::debug_assert_not_null, "C++: Teuchos::RCP<ROL::RiddersProjection<double> >::debug_assert_not_null() const --> const class Teuchos::RCP<class ROL::RiddersProjection<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<class ROL::RiddersProjection<double> > & (Teuchos::RCP<ROL::RiddersProjection<double>>::*)() const) &Teuchos::RCP<ROL::RiddersProjection<double> >::debug_assert_valid_ptr, "C++: Teuchos::RCP<ROL::RiddersProjection<double> >::debug_assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::RiddersProjection<double> > &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<ROL::RiddersProjection<double>>::*)()) &Teuchos::RCP<ROL::RiddersProjection<double> >::reset, "C++: Teuchos::RCP<ROL::RiddersProjection<double> >::reset() --> void");
		cl.def("access_private_ptr", (class ROL::RiddersProjection<double> * (Teuchos::RCP<ROL::RiddersProjection<double>>::*)() const) &Teuchos::RCP<ROL::RiddersProjection<double> >::access_private_ptr, "C++: Teuchos::RCP<ROL::RiddersProjection<double> >::access_private_ptr() const --> class ROL::RiddersProjection<double> *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::RiddersProjection<double>>::*)()) &Teuchos::RCP<ROL::RiddersProjection<double> >::nonconst_access_private_node, "C++: Teuchos::RCP<ROL::RiddersProjection<double> >::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::RiddersProjection<double>>::*)() const) &Teuchos::RCP<ROL::RiddersProjection<double> >::access_private_node, "C++: Teuchos::RCP<ROL::RiddersProjection<double> >::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<ROL::RiddersProjection<double>> const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<ROL::BrentsProjection<double>>> cl(M("Teuchos"), "RCP_ROL_BrentsProjection_double_t", "");
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<ROL::BrentsProjection<double>>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("null_arg") );

		cl.def( pybind11::init( [](class ROL::BrentsProjection<double> * a0){ return new Teuchos::RCP<ROL::BrentsProjection<double>>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<class ROL::BrentsProjection<double> *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership_in") );

		cl.def( pybind11::init( [](Teuchos::RCP<ROL::BrentsProjection<double>> const &o){ return new Teuchos::RCP<ROL::BrentsProjection<double>>(o); } ) );
		cl.def( pybind11::init<class ROL::BrentsProjection<double> *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::BrentsProjection<double> *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class ROL::BrentsProjection<double> *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<class ROL::BrentsProjection<double> > & (Teuchos::RCP<ROL::BrentsProjection<double>>::*)(const class Teuchos::RCP<class ROL::BrentsProjection<double> > &)) &Teuchos::RCP<ROL::BrentsProjection<double> >::operator=, "C++: Teuchos::RCP<ROL::BrentsProjection<double> >::operator=(const class Teuchos::RCP<class ROL::BrentsProjection<double> > &) --> class Teuchos::RCP<class ROL::BrentsProjection<double> > &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<class ROL::BrentsProjection<double> > & (Teuchos::RCP<ROL::BrentsProjection<double>>::*)(enum Teuchos::ENull)) &Teuchos::RCP<ROL::BrentsProjection<double> >::operator=, "C++: Teuchos::RCP<ROL::BrentsProjection<double> >::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<class ROL::BrentsProjection<double> > &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<ROL::BrentsProjection<double>>::*)(class Teuchos::RCP<class ROL::BrentsProjection<double> > &)) &Teuchos::RCP<ROL::BrentsProjection<double> >::swap, "C++: Teuchos::RCP<ROL::BrentsProjection<double> >::swap(class Teuchos::RCP<class ROL::BrentsProjection<double> > &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<ROL::BrentsProjection<double>>::*)() const) &Teuchos::RCP<ROL::BrentsProjection<double> >::is_null, "C++: Teuchos::RCP<ROL::BrentsProjection<double> >::is_null() const --> bool");
		cl.def("arrow", (class ROL::BrentsProjection<double> * (Teuchos::RCP<ROL::BrentsProjection<double>>::*)() const) &Teuchos::RCP<ROL::BrentsProjection<double> >::operator->, "C++: Teuchos::RCP<ROL::BrentsProjection<double> >::operator->() const --> class ROL::BrentsProjection<double> *", pybind11::return_value_policy::automatic);
		cl.def("__mul__", (class ROL::BrentsProjection<double> & (Teuchos::RCP<ROL::BrentsProjection<double>>::*)() const) &Teuchos::RCP<ROL::BrentsProjection<double> >::operator*, "C++: Teuchos::RCP<ROL::BrentsProjection<double> >::operator*() const --> class ROL::BrentsProjection<double> &", pybind11::return_value_policy::automatic);
		cl.def("get", (class ROL::BrentsProjection<double> * (Teuchos::RCP<ROL::BrentsProjection<double>>::*)() const) &Teuchos::RCP<ROL::BrentsProjection<double> >::get, "C++: Teuchos::RCP<ROL::BrentsProjection<double> >::get() const --> class ROL::BrentsProjection<double> *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class ROL::BrentsProjection<double> * (Teuchos::RCP<ROL::BrentsProjection<double>>::*)() const) &Teuchos::RCP<ROL::BrentsProjection<double> >::getRawPtr, "C++: Teuchos::RCP<ROL::BrentsProjection<double> >::getRawPtr() const --> class ROL::BrentsProjection<double> *", pybind11::return_value_policy::automatic);
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<ROL::BrentsProjection<double>>::*)() const) &Teuchos::RCP<ROL::BrentsProjection<double> >::strength, "C++: Teuchos::RCP<ROL::BrentsProjection<double> >::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<ROL::BrentsProjection<double>>::*)() const) &Teuchos::RCP<ROL::BrentsProjection<double> >::is_valid_ptr, "C++: Teuchos::RCP<ROL::BrentsProjection<double> >::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<ROL::BrentsProjection<double>>::*)() const) &Teuchos::RCP<ROL::BrentsProjection<double> >::strong_count, "C++: Teuchos::RCP<ROL::BrentsProjection<double> >::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<ROL::BrentsProjection<double>>::*)() const) &Teuchos::RCP<ROL::BrentsProjection<double> >::weak_count, "C++: Teuchos::RCP<ROL::BrentsProjection<double> >::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<ROL::BrentsProjection<double>>::*)() const) &Teuchos::RCP<ROL::BrentsProjection<double> >::total_count, "C++: Teuchos::RCP<ROL::BrentsProjection<double> >::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<ROL::BrentsProjection<double>>::*)()) &Teuchos::RCP<ROL::BrentsProjection<double> >::set_has_ownership, "C++: Teuchos::RCP<ROL::BrentsProjection<double> >::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<ROL::BrentsProjection<double>>::*)() const) &Teuchos::RCP<ROL::BrentsProjection<double> >::has_ownership, "C++: Teuchos::RCP<ROL::BrentsProjection<double> >::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<class ROL::BrentsProjection<double> > (Teuchos::RCP<ROL::BrentsProjection<double>>::*)() const) &Teuchos::RCP<ROL::BrentsProjection<double> >::create_weak, "C++: Teuchos::RCP<ROL::BrentsProjection<double> >::create_weak() const --> class Teuchos::RCP<class ROL::BrentsProjection<double> >");
		cl.def("create_strong", (class Teuchos::RCP<class ROL::BrentsProjection<double> > (Teuchos::RCP<ROL::BrentsProjection<double>>::*)() const) &Teuchos::RCP<ROL::BrentsProjection<double> >::create_strong, "C++: Teuchos::RCP<ROL::BrentsProjection<double> >::create_strong() const --> class Teuchos::RCP<class ROL::BrentsProjection<double> >");
		cl.def("assert_not_null", (const class Teuchos::RCP<class ROL::BrentsProjection<double> > & (Teuchos::RCP<ROL::BrentsProjection<double>>::*)() const) &Teuchos::RCP<ROL::BrentsProjection<double> >::assert_not_null, "C++: Teuchos::RCP<ROL::BrentsProjection<double> >::assert_not_null() const --> const class Teuchos::RCP<class ROL::BrentsProjection<double> > &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<class ROL::BrentsProjection<double> > & (Teuchos::RCP<ROL::BrentsProjection<double>>::*)() const) &Teuchos::RCP<ROL::BrentsProjection<double> >::assert_valid_ptr, "C++: Teuchos::RCP<ROL::BrentsProjection<double> >::assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::BrentsProjection<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<class ROL::BrentsProjection<double> > & (Teuchos::RCP<ROL::BrentsProjection<double>>::*)() const) &Teuchos::RCP<ROL::BrentsProjection<double> >::debug_assert_not_null, "C++: Teuchos::RCP<ROL::BrentsProjection<double> >::debug_assert_not_null() const --> const class Teuchos::RCP<class ROL::BrentsProjection<double> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<class ROL::BrentsProjection<double> > & (Teuchos::RCP<ROL::BrentsProjection<double>>::*)() const) &Teuchos::RCP<ROL::BrentsProjection<double> >::debug_assert_valid_ptr, "C++: Teuchos::RCP<ROL::BrentsProjection<double> >::debug_assert_valid_ptr() const --> const class Teuchos::RCP<class ROL::BrentsProjection<double> > &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<ROL::BrentsProjection<double>>::*)()) &Teuchos::RCP<ROL::BrentsProjection<double> >::reset, "C++: Teuchos::RCP<ROL::BrentsProjection<double> >::reset() --> void");
		cl.def("access_private_ptr", (class ROL::BrentsProjection<double> * (Teuchos::RCP<ROL::BrentsProjection<double>>::*)() const) &Teuchos::RCP<ROL::BrentsProjection<double> >::access_private_ptr, "C++: Teuchos::RCP<ROL::BrentsProjection<double> >::access_private_ptr() const --> class ROL::BrentsProjection<double> *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::BrentsProjection<double>>::*)()) &Teuchos::RCP<ROL::BrentsProjection<double> >::nonconst_access_private_node, "C++: Teuchos::RCP<ROL::BrentsProjection<double> >::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<ROL::BrentsProjection<double>>::*)() const) &Teuchos::RCP<ROL::BrentsProjection<double> >::access_private_node, "C++: Teuchos::RCP<ROL::BrentsProjection<double> >::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<ROL::BrentsProjection<double>> const &o) -> std::string { std::ostringstream s; s << o; return s.str(); } );
	}
	{ // Teuchos::RCPComp file:Teuchos_RCPDecl.hpp line:955
		pybind11::class_<Teuchos::RCPComp> cl(M("Teuchos"), "RCPComp", "Struct for comparing two RCPs. Simply compares\n the raw pointers contained within the RCPs");
		cl.def( pybind11::init( [](Teuchos::RCPComp const &o){ return new Teuchos::RCPComp(o); } ) );
		cl.def( pybind11::init( [](){ return new Teuchos::RCPComp(); } ) );
	}
	{ // Teuchos::RCPConstComp file:Teuchos_RCPDecl.hpp line:965
		pybind11::class_<Teuchos::RCPConstComp> cl(M("Teuchos"), "RCPConstComp", "Struct for comparing two RCPs. Simply compares\n the raw pointers contained within the RCPs");
		cl.def( pybind11::init( [](Teuchos::RCPConstComp const &o){ return new Teuchos::RCPConstComp(o); } ) );
		cl.def( pybind11::init( [](){ return new Teuchos::RCPConstComp(); } ) );
		cl.def("__call__", (bool (Teuchos::RCPConstComp::*)(const class Teuchos::RCP<const class Teuchos::ParameterEntry>, const class Teuchos::RCP<const class Teuchos::ParameterEntry>) const) &Teuchos::RCPConstComp::operator()<Teuchos::ParameterEntry,Teuchos::ParameterEntry>, "C++: Teuchos::RCPConstComp::operator()(const class Teuchos::RCP<const class Teuchos::ParameterEntry>, const class Teuchos::RCP<const class Teuchos::ParameterEntry>) const --> bool", pybind11::arg("p1"), pybind11::arg("p2"));
	}
	// Teuchos::rcp(class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > * a0) -> Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > > (*)(class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > *, bool)) &Teuchos::rcp<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>, "C++: Teuchos::rcp(class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > *, bool) --> class Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > * a0) -> Teuchos::RCP<class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > > (*)(class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > *, bool)) &Teuchos::rcp<Teuchos::basic_oblackholestream<char, std::char_traits<char> >>, "C++: Teuchos::rcp(class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > *, bool) --> class Teuchos::RCP<class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class std::basic_ostringstream<char> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class std::basic_ostringstream<char> * a0) -> Teuchos::RCP<class std::basic_ostringstream<char> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class std::basic_ostringstream<char> > (*)(class std::basic_ostringstream<char> *, bool)) &Teuchos::rcp<std::basic_ostringstream<char>>, "C++: Teuchos::rcp(class std::basic_ostringstream<char> *, bool) --> class Teuchos::RCP<class std::basic_ostringstream<char> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class Teuchos::ParameterList *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class Teuchos::ParameterList * a0) -> Teuchos::RCP<class Teuchos::ParameterList> { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class Teuchos::ParameterList> (*)(class Teuchos::ParameterList *, bool)) &Teuchos::rcp<Teuchos::ParameterList>, "C++: Teuchos::rcp(class Teuchos::ParameterList *, bool) --> class Teuchos::RCP<class Teuchos::ParameterList>", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::TypeU::BundleAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::TypeU::BundleAlgorithm<double> * a0) -> Teuchos::RCP<class ROL::TypeU::BundleAlgorithm<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::TypeU::BundleAlgorithm<double> > (*)(class ROL::TypeU::BundleAlgorithm<double> *, bool)) &Teuchos::rcp<ROL::TypeU::BundleAlgorithm<double>>, "C++: Teuchos::rcp(class ROL::TypeU::BundleAlgorithm<double> *, bool) --> class Teuchos::RCP<class ROL::TypeU::BundleAlgorithm<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::CombinedStatusTest<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::CombinedStatusTest<double> * a0) -> Teuchos::RCP<class ROL::CombinedStatusTest<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::CombinedStatusTest<double> > (*)(class ROL::CombinedStatusTest<double> *, bool)) &Teuchos::rcp<ROL::CombinedStatusTest<double>>, "C++: Teuchos::rcp(class ROL::CombinedStatusTest<double> *, bool) --> class Teuchos::RCP<class ROL::CombinedStatusTest<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(struct ROL::TypeU::AlgorithmState<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](struct ROL::TypeU::AlgorithmState<double> * a0) -> Teuchos::RCP<struct ROL::TypeU::AlgorithmState<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<struct ROL::TypeU::AlgorithmState<double> > (*)(struct ROL::TypeU::AlgorithmState<double> *, bool)) &Teuchos::rcp<ROL::TypeU::AlgorithmState<double>>, "C++: Teuchos::rcp(struct ROL::TypeU::AlgorithmState<double> *, bool) --> class Teuchos::RCP<struct ROL::TypeU::AlgorithmState<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::StatusTest<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::StatusTest<double> * a0) -> Teuchos::RCP<class ROL::StatusTest<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::StatusTest<double> > (*)(class ROL::StatusTest<double> *, bool)) &Teuchos::rcp<ROL::StatusTest<double>>, "C++: Teuchos::rcp(class ROL::StatusTest<double> *, bool) --> class Teuchos::RCP<class ROL::StatusTest<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::BundleStatusTest<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::BundleStatusTest<double> * a0) -> Teuchos::RCP<class ROL::BundleStatusTest<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::BundleStatusTest<double> > (*)(class ROL::BundleStatusTest<double> *, bool)) &Teuchos::rcp<ROL::BundleStatusTest<double>>, "C++: Teuchos::rcp(class ROL::BundleStatusTest<double> *, bool) --> class Teuchos::RCP<class ROL::BundleStatusTest<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::Bundle_U_TT<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::Bundle_U_TT<double> * a0) -> Teuchos::RCP<class ROL::Bundle_U_TT<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::Bundle_U_TT<double> > (*)(class ROL::Bundle_U_TT<double> *, bool)) &Teuchos::rcp<ROL::Bundle_U_TT<double>>, "C++: Teuchos::rcp(class ROL::Bundle_U_TT<double> *, bool) --> class Teuchos::RCP<class ROL::Bundle_U_TT<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::Bundle_U_AS<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::Bundle_U_AS<double> * a0) -> Teuchos::RCP<class ROL::Bundle_U_AS<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::Bundle_U_AS<double> > (*)(class ROL::Bundle_U_AS<double> *, bool)) &Teuchos::rcp<ROL::Bundle_U_AS<double>>, "C++: Teuchos::rcp(class ROL::Bundle_U_AS<double> *, bool) --> class Teuchos::RCP<class ROL::Bundle_U_AS<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::IterationScaling_U<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::IterationScaling_U<double> * a0) -> Teuchos::RCP<class ROL::IterationScaling_U<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::IterationScaling_U<double> > (*)(class ROL::IterationScaling_U<double> *, bool)) &Teuchos::rcp<ROL::IterationScaling_U<double>>, "C++: Teuchos::rcp(class ROL::IterationScaling_U<double> *, bool) --> class Teuchos::RCP<class ROL::IterationScaling_U<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::PathBasedTargetLevel_U<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::PathBasedTargetLevel_U<double> * a0) -> Teuchos::RCP<class ROL::PathBasedTargetLevel_U<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::PathBasedTargetLevel_U<double> > (*)(class ROL::PathBasedTargetLevel_U<double> *, bool)) &Teuchos::rcp<ROL::PathBasedTargetLevel_U<double>>, "C++: Teuchos::rcp(class ROL::PathBasedTargetLevel_U<double> *, bool) --> class Teuchos::RCP<class ROL::PathBasedTargetLevel_U<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::BackTracking_U<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::BackTracking_U<double> * a0) -> Teuchos::RCP<class ROL::BackTracking_U<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::BackTracking_U<double> > (*)(class ROL::BackTracking_U<double> *, bool)) &Teuchos::rcp<ROL::BackTracking_U<double>>, "C++: Teuchos::rcp(class ROL::BackTracking_U<double> *, bool) --> class Teuchos::RCP<class ROL::BackTracking_U<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::CubicInterp_U<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::CubicInterp_U<double> * a0) -> Teuchos::RCP<class ROL::CubicInterp_U<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::CubicInterp_U<double> > (*)(class ROL::CubicInterp_U<double> *, bool)) &Teuchos::rcp<ROL::CubicInterp_U<double>>, "C++: Teuchos::rcp(class ROL::CubicInterp_U<double> *, bool) --> class Teuchos::RCP<class ROL::CubicInterp_U<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::ScalarMinimizationLineSearch_U<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::ScalarMinimizationLineSearch_U<double> * a0) -> Teuchos::RCP<class ROL::ScalarMinimizationLineSearch_U<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::ScalarMinimizationLineSearch_U<double> > (*)(class ROL::ScalarMinimizationLineSearch_U<double> *, bool)) &Teuchos::rcp<ROL::ScalarMinimizationLineSearch_U<double>>, "C++: Teuchos::rcp(class ROL::ScalarMinimizationLineSearch_U<double> *, bool) --> class Teuchos::RCP<class ROL::ScalarMinimizationLineSearch_U<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::Bracketing<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::Bracketing<double> * a0) -> Teuchos::RCP<class ROL::Bracketing<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::Bracketing<double> > (*)(class ROL::Bracketing<double> *, bool)) &Teuchos::rcp<ROL::Bracketing<double>>, "C++: Teuchos::rcp(class ROL::Bracketing<double> *, bool) --> class Teuchos::RCP<class ROL::Bracketing<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::BrentsScalarMinimization<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::BrentsScalarMinimization<double> * a0) -> Teuchos::RCP<class ROL::BrentsScalarMinimization<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::BrentsScalarMinimization<double> > (*)(class ROL::BrentsScalarMinimization<double> *, bool)) &Teuchos::rcp<ROL::BrentsScalarMinimization<double>>, "C++: Teuchos::rcp(class ROL::BrentsScalarMinimization<double> *, bool) --> class Teuchos::RCP<class ROL::BrentsScalarMinimization<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::BisectionScalarMinimization<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::BisectionScalarMinimization<double> * a0) -> Teuchos::RCP<class ROL::BisectionScalarMinimization<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::BisectionScalarMinimization<double> > (*)(class ROL::BisectionScalarMinimization<double> *, bool)) &Teuchos::rcp<ROL::BisectionScalarMinimization<double>>, "C++: Teuchos::rcp(class ROL::BisectionScalarMinimization<double> *, bool) --> class Teuchos::RCP<class ROL::BisectionScalarMinimization<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::GoldenSectionScalarMinimization<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::GoldenSectionScalarMinimization<double> * a0) -> Teuchos::RCP<class ROL::GoldenSectionScalarMinimization<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::GoldenSectionScalarMinimization<double> > (*)(class ROL::GoldenSectionScalarMinimization<double> *, bool)) &Teuchos::rcp<ROL::GoldenSectionScalarMinimization<double>>, "C++: Teuchos::rcp(class ROL::GoldenSectionScalarMinimization<double> *, bool) --> class Teuchos::RCP<class ROL::GoldenSectionScalarMinimization<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::NullSpaceOperator<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::NullSpaceOperator<double> * a0) -> Teuchos::RCP<class ROL::NullSpaceOperator<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::NullSpaceOperator<double> > (*)(class ROL::NullSpaceOperator<double> *, bool)) &Teuchos::rcp<ROL::NullSpaceOperator<double>>, "C++: Teuchos::rcp(class ROL::NullSpaceOperator<double> *, bool) --> class Teuchos::RCP<class ROL::NullSpaceOperator<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::TypeU::LineSearchAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::TypeU::LineSearchAlgorithm<double> * a0) -> Teuchos::RCP<class ROL::TypeU::LineSearchAlgorithm<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::TypeU::LineSearchAlgorithm<double> > (*)(class ROL::TypeU::LineSearchAlgorithm<double> *, bool)) &Teuchos::rcp<ROL::TypeU::LineSearchAlgorithm<double>>, "C++: Teuchos::rcp(class ROL::TypeU::LineSearchAlgorithm<double> *, bool) --> class Teuchos::RCP<class ROL::TypeU::LineSearchAlgorithm<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::Gradient_U<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::Gradient_U<double> * a0) -> Teuchos::RCP<class ROL::Gradient_U<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::Gradient_U<double> > (*)(class ROL::Gradient_U<double> *, bool)) &Teuchos::rcp<ROL::Gradient_U<double>>, "C++: Teuchos::rcp(class ROL::Gradient_U<double> *, bool) --> class Teuchos::RCP<class ROL::Gradient_U<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::NonlinearCG_U<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::NonlinearCG_U<double> * a0) -> Teuchos::RCP<class ROL::NonlinearCG_U<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::NonlinearCG_U<double> > (*)(class ROL::NonlinearCG_U<double> *, bool)) &Teuchos::rcp<ROL::NonlinearCG_U<double>>, "C++: Teuchos::rcp(class ROL::NonlinearCG_U<double> *, bool) --> class Teuchos::RCP<class ROL::NonlinearCG_U<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::NonlinearCG<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::NonlinearCG<double> * a0) -> Teuchos::RCP<class ROL::NonlinearCG<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::NonlinearCG<double> > (*)(class ROL::NonlinearCG<double> *, bool)) &Teuchos::rcp<ROL::NonlinearCG<double>>, "C++: Teuchos::rcp(class ROL::NonlinearCG<double> *, bool) --> class Teuchos::RCP<class ROL::NonlinearCG<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(struct ROL::NonlinearCGState<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](struct ROL::NonlinearCGState<double> * a0) -> Teuchos::RCP<struct ROL::NonlinearCGState<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<struct ROL::NonlinearCGState<double> > (*)(struct ROL::NonlinearCGState<double> *, bool)) &Teuchos::rcp<ROL::NonlinearCGState<double>>, "C++: Teuchos::rcp(struct ROL::NonlinearCGState<double> *, bool) --> class Teuchos::RCP<struct ROL::NonlinearCGState<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::QuasiNewton_U<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::QuasiNewton_U<double> * a0) -> Teuchos::RCP<class ROL::QuasiNewton_U<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::QuasiNewton_U<double> > (*)(class ROL::QuasiNewton_U<double> *, bool)) &Teuchos::rcp<ROL::QuasiNewton_U<double>>, "C++: Teuchos::rcp(class ROL::QuasiNewton_U<double> *, bool) --> class Teuchos::RCP<class ROL::QuasiNewton_U<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::lBFGS<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::lBFGS<double> * a0) -> Teuchos::RCP<class ROL::lBFGS<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::lBFGS<double> > (*)(class ROL::lBFGS<double> *, bool)) &Teuchos::rcp<ROL::lBFGS<double>>, "C++: Teuchos::rcp(class ROL::lBFGS<double> *, bool) --> class Teuchos::RCP<class ROL::lBFGS<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(struct ROL::SecantState<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](struct ROL::SecantState<double> * a0) -> Teuchos::RCP<struct ROL::SecantState<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<struct ROL::SecantState<double> > (*)(struct ROL::SecantState<double> *, bool)) &Teuchos::rcp<ROL::SecantState<double>>, "C++: Teuchos::rcp(struct ROL::SecantState<double> *, bool) --> class Teuchos::RCP<struct ROL::SecantState<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::lDFP<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::lDFP<double> * a0) -> Teuchos::RCP<class ROL::lDFP<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::lDFP<double> > (*)(class ROL::lDFP<double> *, bool)) &Teuchos::rcp<ROL::lDFP<double>>, "C++: Teuchos::rcp(class ROL::lDFP<double> *, bool) --> class Teuchos::RCP<class ROL::lDFP<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::lSR1<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::lSR1<double> * a0) -> Teuchos::RCP<class ROL::lSR1<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::lSR1<double> > (*)(class ROL::lSR1<double> *, bool)) &Teuchos::rcp<ROL::lSR1<double>>, "C++: Teuchos::rcp(class ROL::lSR1<double> *, bool) --> class Teuchos::RCP<class ROL::lSR1<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::BarzilaiBorwein<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::BarzilaiBorwein<double> * a0) -> Teuchos::RCP<class ROL::BarzilaiBorwein<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::BarzilaiBorwein<double> > (*)(class ROL::BarzilaiBorwein<double> *, bool)) &Teuchos::rcp<ROL::BarzilaiBorwein<double>>, "C++: Teuchos::rcp(class ROL::BarzilaiBorwein<double> *, bool) --> class Teuchos::RCP<class ROL::BarzilaiBorwein<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::Newton_U<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::Newton_U<double> * a0) -> Teuchos::RCP<class ROL::Newton_U<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::Newton_U<double> > (*)(class ROL::Newton_U<double> *, bool)) &Teuchos::rcp<ROL::Newton_U<double>>, "C++: Teuchos::rcp(class ROL::Newton_U<double> *, bool) --> class Teuchos::RCP<class ROL::Newton_U<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::NewtonKrylov_U<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::NewtonKrylov_U<double> * a0) -> Teuchos::RCP<class ROL::NewtonKrylov_U<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::NewtonKrylov_U<double> > (*)(class ROL::NewtonKrylov_U<double> *, bool)) &Teuchos::rcp<ROL::NewtonKrylov_U<double>>, "C++: Teuchos::rcp(class ROL::NewtonKrylov_U<double> *, bool) --> class Teuchos::RCP<class ROL::NewtonKrylov_U<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::ConjugateResiduals<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::ConjugateResiduals<double> * a0) -> Teuchos::RCP<class ROL::ConjugateResiduals<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::ConjugateResiduals<double> > (*)(class ROL::ConjugateResiduals<double> *, bool)) &Teuchos::rcp<ROL::ConjugateResiduals<double>>, "C++: Teuchos::rcp(class ROL::ConjugateResiduals<double> *, bool) --> class Teuchos::RCP<class ROL::ConjugateResiduals<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::ConjugateGradients<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::ConjugateGradients<double> * a0) -> Teuchos::RCP<class ROL::ConjugateGradients<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::ConjugateGradients<double> > (*)(class ROL::ConjugateGradients<double> *, bool)) &Teuchos::rcp<ROL::ConjugateGradients<double>>, "C++: Teuchos::rcp(class ROL::ConjugateGradients<double> *, bool) --> class Teuchos::RCP<class ROL::ConjugateGradients<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::GMRES<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::GMRES<double> * a0) -> Teuchos::RCP<class ROL::GMRES<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::GMRES<double> > (*)(class ROL::GMRES<double> *, bool)) &Teuchos::rcp<ROL::GMRES<double>>, "C++: Teuchos::rcp(class ROL::GMRES<double> *, bool) --> class Teuchos::RCP<class ROL::GMRES<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class Teuchos::SerialDenseMatrix<int, double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class Teuchos::SerialDenseMatrix<int, double> * a0) -> Teuchos::RCP<class Teuchos::SerialDenseMatrix<int, double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class Teuchos::SerialDenseMatrix<int, double> > (*)(class Teuchos::SerialDenseMatrix<int, double> *, bool)) &Teuchos::rcp<Teuchos::SerialDenseMatrix<int, double>>, "C++: Teuchos::rcp(class Teuchos::SerialDenseMatrix<int, double> *, bool) --> class Teuchos::RCP<class Teuchos::SerialDenseMatrix<int, double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class Teuchos::SerialDenseVector<int, double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class Teuchos::SerialDenseVector<int, double> * a0) -> Teuchos::RCP<class Teuchos::SerialDenseVector<int, double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class Teuchos::SerialDenseVector<int, double> > (*)(class Teuchos::SerialDenseVector<int, double> *, bool)) &Teuchos::rcp<Teuchos::SerialDenseVector<int, double>>, "C++: Teuchos::rcp(class Teuchos::SerialDenseVector<int, double> *, bool) --> class Teuchos::RCP<class Teuchos::SerialDenseVector<int, double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::TypeU::TrustRegionAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::TypeU::TrustRegionAlgorithm<double> * a0) -> Teuchos::RCP<class ROL::TypeU::TrustRegionAlgorithm<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::TypeU::TrustRegionAlgorithm<double> > (*)(class ROL::TypeU::TrustRegionAlgorithm<double> *, bool)) &Teuchos::rcp<ROL::TypeU::TrustRegionAlgorithm<double>>, "C++: Teuchos::rcp(class ROL::TypeU::TrustRegionAlgorithm<double> *, bool) --> class Teuchos::RCP<class ROL::TypeU::TrustRegionAlgorithm<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::CauchyPoint_U<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::CauchyPoint_U<double> * a0) -> Teuchos::RCP<class ROL::CauchyPoint_U<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::CauchyPoint_U<double> > (*)(class ROL::CauchyPoint_U<double> *, bool)) &Teuchos::rcp<ROL::CauchyPoint_U<double>>, "C++: Teuchos::rcp(class ROL::CauchyPoint_U<double> *, bool) --> class Teuchos::RCP<class ROL::CauchyPoint_U<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::DogLeg_U<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::DogLeg_U<double> * a0) -> Teuchos::RCP<class ROL::DogLeg_U<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::DogLeg_U<double> > (*)(class ROL::DogLeg_U<double> *, bool)) &Teuchos::rcp<ROL::DogLeg_U<double>>, "C++: Teuchos::rcp(class ROL::DogLeg_U<double> *, bool) --> class Teuchos::RCP<class ROL::DogLeg_U<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::DoubleDogLeg_U<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::DoubleDogLeg_U<double> * a0) -> Teuchos::RCP<class ROL::DoubleDogLeg_U<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::DoubleDogLeg_U<double> > (*)(class ROL::DoubleDogLeg_U<double> *, bool)) &Teuchos::rcp<ROL::DoubleDogLeg_U<double>>, "C++: Teuchos::rcp(class ROL::DoubleDogLeg_U<double> *, bool) --> class Teuchos::RCP<class ROL::DoubleDogLeg_U<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::TruncatedCG_U<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::TruncatedCG_U<double> * a0) -> Teuchos::RCP<class ROL::TruncatedCG_U<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::TruncatedCG_U<double> > (*)(class ROL::TruncatedCG_U<double> *, bool)) &Teuchos::rcp<ROL::TruncatedCG_U<double>>, "C++: Teuchos::rcp(class ROL::TruncatedCG_U<double> *, bool) --> class Teuchos::RCP<class ROL::TruncatedCG_U<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::SPGTrustRegion_U<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::SPGTrustRegion_U<double> * a0) -> Teuchos::RCP<class ROL::SPGTrustRegion_U<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::SPGTrustRegion_U<double> > (*)(class ROL::SPGTrustRegion_U<double> *, bool)) &Teuchos::rcp<ROL::SPGTrustRegion_U<double>>, "C++: Teuchos::rcp(class ROL::SPGTrustRegion_U<double> *, bool) --> class Teuchos::RCP<class ROL::SPGTrustRegion_U<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::TrustRegionModel_U<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::TrustRegionModel_U<double> * a0) -> Teuchos::RCP<class ROL::TrustRegionModel_U<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::TrustRegionModel_U<double> > (*)(class ROL::TrustRegionModel_U<double> *, bool)) &Teuchos::rcp<ROL::TrustRegionModel_U<double>>, "C++: Teuchos::rcp(class ROL::TrustRegionModel_U<double> *, bool) --> class Teuchos::RCP<class ROL::TrustRegionModel_U<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::TypeB::NewtonKrylovAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::TypeB::NewtonKrylovAlgorithm<double> * a0) -> Teuchos::RCP<class ROL::TypeB::NewtonKrylovAlgorithm<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::TypeB::NewtonKrylovAlgorithm<double> > (*)(class ROL::TypeB::NewtonKrylovAlgorithm<double> *, bool)) &Teuchos::rcp<ROL::TypeB::NewtonKrylovAlgorithm<double>>, "C++: Teuchos::rcp(class ROL::TypeB::NewtonKrylovAlgorithm<double> *, bool) --> class Teuchos::RCP<class ROL::TypeB::NewtonKrylovAlgorithm<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(struct ROL::TypeB::AlgorithmState<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](struct ROL::TypeB::AlgorithmState<double> * a0) -> Teuchos::RCP<struct ROL::TypeB::AlgorithmState<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<struct ROL::TypeB::AlgorithmState<double> > (*)(struct ROL::TypeB::AlgorithmState<double> *, bool)) &Teuchos::rcp<ROL::TypeB::AlgorithmState<double>>, "C++: Teuchos::rcp(struct ROL::TypeB::AlgorithmState<double> *, bool) --> class Teuchos::RCP<struct ROL::TypeB::AlgorithmState<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::PolyhedralProjection<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::PolyhedralProjection<double> * a0) -> Teuchos::RCP<class ROL::PolyhedralProjection<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::PolyhedralProjection<double> > (*)(class ROL::PolyhedralProjection<double> *, bool)) &Teuchos::rcp<ROL::PolyhedralProjection<double>>, "C++: Teuchos::rcp(class ROL::PolyhedralProjection<double> *, bool) --> class Teuchos::RCP<class ROL::PolyhedralProjection<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::TypeB::LSecantBAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::TypeB::LSecantBAlgorithm<double> * a0) -> Teuchos::RCP<class ROL::TypeB::LSecantBAlgorithm<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::TypeB::LSecantBAlgorithm<double> > (*)(class ROL::TypeB::LSecantBAlgorithm<double> *, bool)) &Teuchos::rcp<ROL::TypeB::LSecantBAlgorithm<double>>, "C++: Teuchos::rcp(class ROL::TypeB::LSecantBAlgorithm<double> *, bool) --> class Teuchos::RCP<class ROL::TypeB::LSecantBAlgorithm<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::ReducedLinearConstraint<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::ReducedLinearConstraint<double> * a0) -> Teuchos::RCP<class ROL::ReducedLinearConstraint<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::ReducedLinearConstraint<double> > (*)(class ROL::ReducedLinearConstraint<double> *, bool)) &Teuchos::rcp<ROL::ReducedLinearConstraint<double>>, "C++: Teuchos::rcp(class ROL::ReducedLinearConstraint<double> *, bool) --> class Teuchos::RCP<class ROL::ReducedLinearConstraint<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::TypeB::QuasiNewtonAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::TypeB::QuasiNewtonAlgorithm<double> * a0) -> Teuchos::RCP<class ROL::TypeB::QuasiNewtonAlgorithm<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::TypeB::QuasiNewtonAlgorithm<double> > (*)(class ROL::TypeB::QuasiNewtonAlgorithm<double> *, bool)) &Teuchos::rcp<ROL::TypeB::QuasiNewtonAlgorithm<double>>, "C++: Teuchos::rcp(class ROL::TypeB::QuasiNewtonAlgorithm<double> *, bool) --> class Teuchos::RCP<class ROL::TypeB::QuasiNewtonAlgorithm<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::PQNObjective<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::PQNObjective<double> * a0) -> Teuchos::RCP<class ROL::PQNObjective<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::PQNObjective<double> > (*)(class ROL::PQNObjective<double> *, bool)) &Teuchos::rcp<ROL::PQNObjective<double>>, "C++: Teuchos::rcp(class ROL::PQNObjective<double> *, bool) --> class Teuchos::RCP<class ROL::PQNObjective<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::Problem<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::Problem<double> * a0) -> Teuchos::RCP<class ROL::Problem<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::Problem<double> > (*)(class ROL::Problem<double> *, bool)) &Teuchos::rcp<ROL::Problem<double>>, "C++: Teuchos::rcp(class ROL::Problem<double> *, bool) --> class Teuchos::RCP<class ROL::Problem<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::TypeB::GradientAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::TypeB::GradientAlgorithm<double> * a0) -> Teuchos::RCP<class ROL::TypeB::GradientAlgorithm<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::TypeB::GradientAlgorithm<double> > (*)(class ROL::TypeB::GradientAlgorithm<double> *, bool)) &Teuchos::rcp<ROL::TypeB::GradientAlgorithm<double>>, "C++: Teuchos::rcp(class ROL::TypeB::GradientAlgorithm<double> *, bool) --> class Teuchos::RCP<class ROL::TypeB::GradientAlgorithm<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::TypeB::KelleySachsAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::TypeB::KelleySachsAlgorithm<double> * a0) -> Teuchos::RCP<class ROL::TypeB::KelleySachsAlgorithm<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::TypeB::KelleySachsAlgorithm<double> > (*)(class ROL::TypeB::KelleySachsAlgorithm<double> *, bool)) &Teuchos::rcp<ROL::TypeB::KelleySachsAlgorithm<double>>, "C++: Teuchos::rcp(class ROL::TypeB::KelleySachsAlgorithm<double> *, bool) --> class Teuchos::RCP<class ROL::TypeB::KelleySachsAlgorithm<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::TypeB::TrustRegionSPGAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::TypeB::TrustRegionSPGAlgorithm<double> * a0) -> Teuchos::RCP<class ROL::TypeB::TrustRegionSPGAlgorithm<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::TypeB::TrustRegionSPGAlgorithm<double> > (*)(class ROL::TypeB::TrustRegionSPGAlgorithm<double> *, bool)) &Teuchos::rcp<ROL::TypeB::TrustRegionSPGAlgorithm<double>>, "C++: Teuchos::rcp(class ROL::TypeB::TrustRegionSPGAlgorithm<double> *, bool) --> class Teuchos::RCP<class ROL::TypeB::TrustRegionSPGAlgorithm<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::TypeB::ColemanLiAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::TypeB::ColemanLiAlgorithm<double> * a0) -> Teuchos::RCP<class ROL::TypeB::ColemanLiAlgorithm<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::TypeB::ColemanLiAlgorithm<double> > (*)(class ROL::TypeB::ColemanLiAlgorithm<double> *, bool)) &Teuchos::rcp<ROL::TypeB::ColemanLiAlgorithm<double>>, "C++: Teuchos::rcp(class ROL::TypeB::ColemanLiAlgorithm<double> *, bool) --> class Teuchos::RCP<class ROL::TypeB::ColemanLiAlgorithm<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::TypeB::LinMoreAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::TypeB::LinMoreAlgorithm<double> * a0) -> Teuchos::RCP<class ROL::TypeB::LinMoreAlgorithm<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::TypeB::LinMoreAlgorithm<double> > (*)(class ROL::TypeB::LinMoreAlgorithm<double> *, bool)) &Teuchos::rcp<ROL::TypeB::LinMoreAlgorithm<double>>, "C++: Teuchos::rcp(class ROL::TypeB::LinMoreAlgorithm<double> *, bool) --> class Teuchos::RCP<class ROL::TypeB::LinMoreAlgorithm<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::TypeB::MoreauYosidaAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::TypeB::MoreauYosidaAlgorithm<double> * a0) -> Teuchos::RCP<class ROL::TypeB::MoreauYosidaAlgorithm<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::TypeB::MoreauYosidaAlgorithm<double> > (*)(class ROL::TypeB::MoreauYosidaAlgorithm<double> *, bool)) &Teuchos::rcp<ROL::TypeB::MoreauYosidaAlgorithm<double>>, "C++: Teuchos::rcp(class ROL::TypeB::MoreauYosidaAlgorithm<double> *, bool) --> class Teuchos::RCP<class ROL::TypeB::MoreauYosidaAlgorithm<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::ScalarController<double, int> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::ScalarController<double, int> * a0) -> Teuchos::RCP<class ROL::ScalarController<double, int> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::ScalarController<double, int> > (*)(class ROL::ScalarController<double, int> *, bool)) &Teuchos::rcp<ROL::ScalarController<double, int>>, "C++: Teuchos::rcp(class ROL::ScalarController<double, int> *, bool) --> class Teuchos::RCP<class ROL::ScalarController<double, int> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::VectorController<double, int> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::VectorController<double, int> * a0) -> Teuchos::RCP<class ROL::VectorController<double, int> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::VectorController<double, int> > (*)(class ROL::VectorController<double, int> *, bool)) &Teuchos::rcp<ROL::VectorController<double, int>>, "C++: Teuchos::rcp(class ROL::VectorController<double, int> *, bool) --> class Teuchos::RCP<class ROL::VectorController<double, int> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::SingletonVector<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::SingletonVector<double> * a0) -> Teuchos::RCP<class ROL::SingletonVector<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::SingletonVector<double> > (*)(class ROL::SingletonVector<double> *, bool)) &Teuchos::rcp<ROL::SingletonVector<double>>, "C++: Teuchos::rcp(class ROL::SingletonVector<double> *, bool) --> class Teuchos::RCP<class ROL::SingletonVector<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::TypeB::PrimalDualActiveSetAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::TypeB::PrimalDualActiveSetAlgorithm<double> * a0) -> Teuchos::RCP<class ROL::TypeB::PrimalDualActiveSetAlgorithm<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::TypeB::PrimalDualActiveSetAlgorithm<double> > (*)(class ROL::TypeB::PrimalDualActiveSetAlgorithm<double> *, bool)) &Teuchos::rcp<ROL::TypeB::PrimalDualActiveSetAlgorithm<double>>, "C++: Teuchos::rcp(class ROL::TypeB::PrimalDualActiveSetAlgorithm<double> *, bool) --> class Teuchos::RCP<class ROL::TypeB::PrimalDualActiveSetAlgorithm<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::PartitionedVector<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::PartitionedVector<double> * a0) -> Teuchos::RCP<class ROL::PartitionedVector<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::PartitionedVector<double> > (*)(class ROL::PartitionedVector<double> *, bool)) &Teuchos::rcp<ROL::PartitionedVector<double>>, "C++: Teuchos::rcp(class ROL::PartitionedVector<double> *, bool) --> class Teuchos::RCP<class ROL::PartitionedVector<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::TypeB::InteriorPointAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::TypeB::InteriorPointAlgorithm<double> * a0) -> Teuchos::RCP<class ROL::TypeB::InteriorPointAlgorithm<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::TypeB::InteriorPointAlgorithm<double> > (*)(class ROL::TypeB::InteriorPointAlgorithm<double> *, bool)) &Teuchos::rcp<ROL::TypeB::InteriorPointAlgorithm<double>>, "C++: Teuchos::rcp(class ROL::TypeB::InteriorPointAlgorithm<double> *, bool) --> class Teuchos::RCP<class ROL::TypeB::InteriorPointAlgorithm<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::TypeB::SpectralGradientAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::TypeB::SpectralGradientAlgorithm<double> * a0) -> Teuchos::RCP<class ROL::TypeB::SpectralGradientAlgorithm<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::TypeB::SpectralGradientAlgorithm<double> > (*)(class ROL::TypeB::SpectralGradientAlgorithm<double> *, bool)) &Teuchos::rcp<ROL::TypeB::SpectralGradientAlgorithm<double>>, "C++: Teuchos::rcp(class ROL::TypeB::SpectralGradientAlgorithm<double> *, bool) --> class Teuchos::RCP<class ROL::TypeB::SpectralGradientAlgorithm<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::TypeE::AugmentedLagrangianAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::TypeE::AugmentedLagrangianAlgorithm<double> * a0) -> Teuchos::RCP<class ROL::TypeE::AugmentedLagrangianAlgorithm<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::TypeE::AugmentedLagrangianAlgorithm<double> > (*)(class ROL::TypeE::AugmentedLagrangianAlgorithm<double> *, bool)) &Teuchos::rcp<ROL::TypeE::AugmentedLagrangianAlgorithm<double>>, "C++: Teuchos::rcp(class ROL::TypeE::AugmentedLagrangianAlgorithm<double> *, bool) --> class Teuchos::RCP<class ROL::TypeE::AugmentedLagrangianAlgorithm<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(struct ROL::TypeE::AlgorithmState<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](struct ROL::TypeE::AlgorithmState<double> * a0) -> Teuchos::RCP<struct ROL::TypeE::AlgorithmState<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<struct ROL::TypeE::AlgorithmState<double> > (*)(struct ROL::TypeE::AlgorithmState<double> *, bool)) &Teuchos::rcp<ROL::TypeE::AlgorithmState<double>>, "C++: Teuchos::rcp(struct ROL::TypeE::AlgorithmState<double> *, bool) --> class Teuchos::RCP<struct ROL::TypeE::AlgorithmState<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::ConstraintStatusTest<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::ConstraintStatusTest<double> * a0) -> Teuchos::RCP<class ROL::ConstraintStatusTest<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::ConstraintStatusTest<double> > (*)(class ROL::ConstraintStatusTest<double> *, bool)) &Teuchos::rcp<ROL::ConstraintStatusTest<double>>, "C++: Teuchos::rcp(class ROL::ConstraintStatusTest<double> *, bool) --> class Teuchos::RCP<class ROL::ConstraintStatusTest<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::TypeE::FletcherAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::TypeE::FletcherAlgorithm<double> * a0) -> Teuchos::RCP<class ROL::TypeE::FletcherAlgorithm<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::TypeE::FletcherAlgorithm<double> > (*)(class ROL::TypeE::FletcherAlgorithm<double> *, bool)) &Teuchos::rcp<ROL::TypeE::FletcherAlgorithm<double>>, "C++: Teuchos::rcp(class ROL::TypeE::FletcherAlgorithm<double> *, bool) --> class Teuchos::RCP<class ROL::TypeE::FletcherAlgorithm<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::TypeE::CompositeStepAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::TypeE::CompositeStepAlgorithm<double> * a0) -> Teuchos::RCP<class ROL::TypeE::CompositeStepAlgorithm<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::TypeE::CompositeStepAlgorithm<double> > (*)(class ROL::TypeE::CompositeStepAlgorithm<double> *, bool)) &Teuchos::rcp<ROL::TypeE::CompositeStepAlgorithm<double>>, "C++: Teuchos::rcp(class ROL::TypeE::CompositeStepAlgorithm<double> *, bool) --> class Teuchos::RCP<class ROL::TypeE::CompositeStepAlgorithm<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::TypeE::StabilizedLCLAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::TypeE::StabilizedLCLAlgorithm<double> * a0) -> Teuchos::RCP<class ROL::TypeE::StabilizedLCLAlgorithm<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::TypeE::StabilizedLCLAlgorithm<double> > (*)(class ROL::TypeE::StabilizedLCLAlgorithm<double> *, bool)) &Teuchos::rcp<ROL::TypeE::StabilizedLCLAlgorithm<double>>, "C++: Teuchos::rcp(class ROL::TypeE::StabilizedLCLAlgorithm<double> *, bool) --> class Teuchos::RCP<class ROL::TypeE::StabilizedLCLAlgorithm<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::AugmentedLagrangianObjective<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::AugmentedLagrangianObjective<double> * a0) -> Teuchos::RCP<class ROL::AugmentedLagrangianObjective<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::AugmentedLagrangianObjective<double> > (*)(class ROL::AugmentedLagrangianObjective<double> *, bool)) &Teuchos::rcp<ROL::AugmentedLagrangianObjective<double>>, "C++: Teuchos::rcp(class ROL::AugmentedLagrangianObjective<double> *, bool) --> class Teuchos::RCP<class ROL::AugmentedLagrangianObjective<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::ElasticLinearConstraint<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::ElasticLinearConstraint<double> * a0) -> Teuchos::RCP<class ROL::ElasticLinearConstraint<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::ElasticLinearConstraint<double> > (*)(class ROL::ElasticLinearConstraint<double> *, bool)) &Teuchos::rcp<ROL::ElasticLinearConstraint<double>>, "C++: Teuchos::rcp(class ROL::ElasticLinearConstraint<double> *, bool) --> class Teuchos::RCP<class ROL::ElasticLinearConstraint<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::BoundConstraint<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::BoundConstraint<double> * a0) -> Teuchos::RCP<class ROL::BoundConstraint<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::BoundConstraint<double> > (*)(class ROL::BoundConstraint<double> *, bool)) &Teuchos::rcp<ROL::BoundConstraint<double>>, "C++: Teuchos::rcp(class ROL::BoundConstraint<double> *, bool) --> class Teuchos::RCP<class ROL::BoundConstraint<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::Bounds<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::Bounds<double> * a0) -> Teuchos::RCP<class ROL::Bounds<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::Bounds<double> > (*)(class ROL::Bounds<double> *, bool)) &Teuchos::rcp<ROL::Bounds<double>>, "C++: Teuchos::rcp(class ROL::Bounds<double> *, bool) --> class Teuchos::RCP<class ROL::Bounds<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::BoundConstraint_Partitioned<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::BoundConstraint_Partitioned<double> * a0) -> Teuchos::RCP<class ROL::BoundConstraint_Partitioned<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::BoundConstraint_Partitioned<double> > (*)(class ROL::BoundConstraint_Partitioned<double> *, bool)) &Teuchos::rcp<ROL::BoundConstraint_Partitioned<double>>, "C++: Teuchos::rcp(class ROL::BoundConstraint_Partitioned<double> *, bool) --> class Teuchos::RCP<class ROL::BoundConstraint_Partitioned<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::TypeG::AugmentedLagrangianAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::TypeG::AugmentedLagrangianAlgorithm<double> * a0) -> Teuchos::RCP<class ROL::TypeG::AugmentedLagrangianAlgorithm<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::TypeG::AugmentedLagrangianAlgorithm<double> > (*)(class ROL::TypeG::AugmentedLagrangianAlgorithm<double> *, bool)) &Teuchos::rcp<ROL::TypeG::AugmentedLagrangianAlgorithm<double>>, "C++: Teuchos::rcp(class ROL::TypeG::AugmentedLagrangianAlgorithm<double> *, bool) --> class Teuchos::RCP<class ROL::TypeG::AugmentedLagrangianAlgorithm<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(struct ROL::TypeG::AlgorithmState<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](struct ROL::TypeG::AlgorithmState<double> * a0) -> Teuchos::RCP<struct ROL::TypeG::AlgorithmState<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<struct ROL::TypeG::AlgorithmState<double> > (*)(struct ROL::TypeG::AlgorithmState<double> *, bool)) &Teuchos::rcp<ROL::TypeG::AlgorithmState<double>>, "C++: Teuchos::rcp(struct ROL::TypeG::AlgorithmState<double> *, bool) --> class Teuchos::RCP<struct ROL::TypeG::AlgorithmState<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::TypeG::MoreauYosidaAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::TypeG::MoreauYosidaAlgorithm<double> * a0) -> Teuchos::RCP<class ROL::TypeG::MoreauYosidaAlgorithm<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::TypeG::MoreauYosidaAlgorithm<double> > (*)(class ROL::TypeG::MoreauYosidaAlgorithm<double> *, bool)) &Teuchos::rcp<ROL::TypeG::MoreauYosidaAlgorithm<double>>, "C++: Teuchos::rcp(class ROL::TypeG::MoreauYosidaAlgorithm<double> *, bool) --> class Teuchos::RCP<class ROL::TypeG::MoreauYosidaAlgorithm<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::TypeG::InteriorPointAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::TypeG::InteriorPointAlgorithm<double> * a0) -> Teuchos::RCP<class ROL::TypeG::InteriorPointAlgorithm<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::TypeG::InteriorPointAlgorithm<double> > (*)(class ROL::TypeG::InteriorPointAlgorithm<double> *, bool)) &Teuchos::rcp<ROL::TypeG::InteriorPointAlgorithm<double>>, "C++: Teuchos::rcp(class ROL::TypeG::InteriorPointAlgorithm<double> *, bool) --> class Teuchos::RCP<class ROL::TypeG::InteriorPointAlgorithm<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::TypeG::StabilizedLCLAlgorithm<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::TypeG::StabilizedLCLAlgorithm<double> * a0) -> Teuchos::RCP<class ROL::TypeG::StabilizedLCLAlgorithm<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::TypeG::StabilizedLCLAlgorithm<double> > (*)(class ROL::TypeG::StabilizedLCLAlgorithm<double> *, bool)) &Teuchos::rcp<ROL::TypeG::StabilizedLCLAlgorithm<double>>, "C++: Teuchos::rcp(class ROL::TypeG::StabilizedLCLAlgorithm<double> *, bool) --> class Teuchos::RCP<class ROL::TypeG::StabilizedLCLAlgorithm<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::Constraint_Partitioned<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::Constraint_Partitioned<double> * a0) -> Teuchos::RCP<class ROL::Constraint_Partitioned<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::Constraint_Partitioned<double> > (*)(class ROL::Constraint_Partitioned<double> *, bool)) &Teuchos::rcp<ROL::Constraint_Partitioned<double>>, "C++: Teuchos::rcp(class ROL::Constraint_Partitioned<double> *, bool) --> class Teuchos::RCP<class ROL::Constraint_Partitioned<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::SlacklessObjective<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::SlacklessObjective<double> * a0) -> Teuchos::RCP<class ROL::SlacklessObjective<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::SlacklessObjective<double> > (*)(class ROL::SlacklessObjective<double> *, bool)) &Teuchos::rcp<ROL::SlacklessObjective<double>>, "C++: Teuchos::rcp(class ROL::SlacklessObjective<double> *, bool) --> class Teuchos::RCP<class ROL::SlacklessObjective<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::ReduceLinearConstraint<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::ReduceLinearConstraint<double> * a0) -> Teuchos::RCP<class ROL::ReduceLinearConstraint<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::ReduceLinearConstraint<double> > (*)(class ROL::ReduceLinearConstraint<double> *, bool)) &Teuchos::rcp<ROL::ReduceLinearConstraint<double>>, "C++: Teuchos::rcp(class ROL::ReduceLinearConstraint<double> *, bool) --> class Teuchos::RCP<class ROL::ReduceLinearConstraint<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::AffineTransformObjective<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::AffineTransformObjective<double> * a0) -> Teuchos::RCP<class ROL::AffineTransformObjective<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::AffineTransformObjective<double> > (*)(class ROL::AffineTransformObjective<double> *, bool)) &Teuchos::rcp<ROL::AffineTransformObjective<double>>, "C++: Teuchos::rcp(class ROL::AffineTransformObjective<double> *, bool) --> class Teuchos::RCP<class ROL::AffineTransformObjective<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::LinearConstraint<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::LinearConstraint<double> * a0) -> Teuchos::RCP<class ROL::LinearConstraint<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::LinearConstraint<double> > (*)(class ROL::LinearConstraint<double> *, bool)) &Teuchos::rcp<ROL::LinearConstraint<double>>, "C++: Teuchos::rcp(class ROL::LinearConstraint<double> *, bool) --> class Teuchos::RCP<class ROL::LinearConstraint<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::AffineTransformConstraint<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::AffineTransformConstraint<double> * a0) -> Teuchos::RCP<class ROL::AffineTransformConstraint<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::AffineTransformConstraint<double> > (*)(class ROL::AffineTransformConstraint<double> *, bool)) &Teuchos::rcp<ROL::AffineTransformConstraint<double>>, "C++: Teuchos::rcp(class ROL::AffineTransformConstraint<double> *, bool) --> class Teuchos::RCP<class ROL::AffineTransformConstraint<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::DaiFletcherProjection<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::DaiFletcherProjection<double> * a0) -> Teuchos::RCP<class ROL::DaiFletcherProjection<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::DaiFletcherProjection<double> > (*)(class ROL::DaiFletcherProjection<double> *, bool)) &Teuchos::rcp<ROL::DaiFletcherProjection<double>>, "C++: Teuchos::rcp(class ROL::DaiFletcherProjection<double> *, bool) --> class Teuchos::RCP<class ROL::DaiFletcherProjection<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::DykstraProjection<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::DykstraProjection<double> * a0) -> Teuchos::RCP<class ROL::DykstraProjection<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::DykstraProjection<double> > (*)(class ROL::DykstraProjection<double> *, bool)) &Teuchos::rcp<ROL::DykstraProjection<double>>, "C++: Teuchos::rcp(class ROL::DykstraProjection<double> *, bool) --> class Teuchos::RCP<class ROL::DykstraProjection<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::DouglasRachfordProjection<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::DouglasRachfordProjection<double> * a0) -> Teuchos::RCP<class ROL::DouglasRachfordProjection<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::DouglasRachfordProjection<double> > (*)(class ROL::DouglasRachfordProjection<double> *, bool)) &Teuchos::rcp<ROL::DouglasRachfordProjection<double>>, "C++: Teuchos::rcp(class ROL::DouglasRachfordProjection<double> *, bool) --> class Teuchos::RCP<class ROL::DouglasRachfordProjection<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::SemismoothNewtonProjection<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::SemismoothNewtonProjection<double> * a0) -> Teuchos::RCP<class ROL::SemismoothNewtonProjection<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::SemismoothNewtonProjection<double> > (*)(class ROL::SemismoothNewtonProjection<double> *, bool)) &Teuchos::rcp<ROL::SemismoothNewtonProjection<double>>, "C++: Teuchos::rcp(class ROL::SemismoothNewtonProjection<double> *, bool) --> class Teuchos::RCP<class ROL::SemismoothNewtonProjection<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::RiddersProjection<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::RiddersProjection<double> * a0) -> Teuchos::RCP<class ROL::RiddersProjection<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::RiddersProjection<double> > (*)(class ROL::RiddersProjection<double> *, bool)) &Teuchos::rcp<ROL::RiddersProjection<double>>, "C++: Teuchos::rcp(class ROL::RiddersProjection<double> *, bool) --> class Teuchos::RCP<class ROL::RiddersProjection<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::BrentsProjection<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::BrentsProjection<double> * a0) -> Teuchos::RCP<class ROL::BrentsProjection<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::BrentsProjection<double> > (*)(class ROL::BrentsProjection<double> *, bool)) &Teuchos::rcp<ROL::BrentsProjection<double>>, "C++: Teuchos::rcp(class ROL::BrentsProjection<double> *, bool) --> class Teuchos::RCP<class ROL::BrentsProjection<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcpFromRef(class Teuchos::ParameterEntry &) file:Teuchos_RCPDecl.hpp line:1297
	M("Teuchos").def("rcpFromRef", (class Teuchos::RCP<class Teuchos::ParameterEntry> (*)(class Teuchos::ParameterEntry &)) &Teuchos::rcpFromRef<Teuchos::ParameterEntry>, "C++: Teuchos::rcpFromRef(class Teuchos::ParameterEntry &) --> class Teuchos::RCP<class Teuchos::ParameterEntry>", pybind11::arg("r"));

	// Teuchos::rcpFromRef(const class Teuchos::ParameterEntry &) file:Teuchos_RCPDecl.hpp line:1297
	M("Teuchos").def("rcpFromRef", (class Teuchos::RCP<const class Teuchos::ParameterEntry> (*)(const class Teuchos::ParameterEntry &)) &Teuchos::rcpFromRef<const Teuchos::ParameterEntry>, "C++: Teuchos::rcpFromRef(const class Teuchos::ParameterEntry &) --> class Teuchos::RCP<const class Teuchos::ParameterEntry>", pybind11::arg("r"));

	// Teuchos::rcpFromRef(const class ROL::Vector<double> &) file:Teuchos_RCP.hpp line:648
	M("Teuchos").def("rcpFromRef", (class Teuchos::RCP<const class ROL::Vector<double> > (*)(const class ROL::Vector<double> &)) &Teuchos::rcpFromRef<const ROL::Vector<double>>, "C++: Teuchos::rcpFromRef(const class ROL::Vector<double> &) --> class Teuchos::RCP<const class ROL::Vector<double> >", pybind11::arg("r"));

	// Teuchos::rcpFromRef(class ROL::Objective<double> &) file:Teuchos_RCP.hpp line:648
	M("Teuchos").def("rcpFromRef", (class Teuchos::RCP<class ROL::Objective<double> > (*)(class ROL::Objective<double> &)) &Teuchos::rcpFromRef<ROL::Objective<double>>, "C++: Teuchos::rcpFromRef(class ROL::Objective<double> &) --> class Teuchos::RCP<class ROL::Objective<double> >", pybind11::arg("r"));

	// Teuchos::rcpFromRef(class ROL::Constraint<double> &) file:Teuchos_RCP.hpp line:648
	M("Teuchos").def("rcpFromRef", (class Teuchos::RCP<class ROL::Constraint<double> > (*)(class ROL::Constraint<double> &)) &Teuchos::rcpFromRef<ROL::Constraint<double>>, "C++: Teuchos::rcpFromRef(class ROL::Constraint<double> &) --> class Teuchos::RCP<class ROL::Constraint<double> >", pybind11::arg("r"));

	// Teuchos::rcpFromRef(class ROL::BoundConstraint<double> &) file:Teuchos_RCP.hpp line:648
	M("Teuchos").def("rcpFromRef", (class Teuchos::RCP<class ROL::BoundConstraint<double> > (*)(class ROL::BoundConstraint<double> &)) &Teuchos::rcpFromRef<ROL::BoundConstraint<double>>, "C++: Teuchos::rcpFromRef(class ROL::BoundConstraint<double> &) --> class Teuchos::RCP<class ROL::BoundConstraint<double> >", pybind11::arg("r"));

	// Teuchos::rcpFromRef(class ROL::Vector<double> &) file:Teuchos_RCP.hpp line:648
	M("Teuchos").def("rcpFromRef", (class Teuchos::RCP<class ROL::Vector<double> > (*)(class ROL::Vector<double> &)) &Teuchos::rcpFromRef<ROL::Vector<double>>, "C++: Teuchos::rcpFromRef(class ROL::Vector<double> &) --> class Teuchos::RCP<class ROL::Vector<double> >", pybind11::arg("r"));

	// Teuchos::rcpFromRef(class ROL::ElasticObjective<double> &) file:Teuchos_RCP.hpp line:648
	M("Teuchos").def("rcpFromRef", (class Teuchos::RCP<class ROL::ElasticObjective<double> > (*)(class ROL::ElasticObjective<double> &)) &Teuchos::rcpFromRef<ROL::ElasticObjective<double>>, "C++: Teuchos::rcpFromRef(class ROL::ElasticObjective<double> &) --> class Teuchos::RCP<class ROL::ElasticObjective<double> >", pybind11::arg("r"));

	// Teuchos::rcpWithEmbeddedObjPostDestroy(class Teuchos::ParameterList *, const class Teuchos::RCP<class Teuchos::ParameterList> &, bool) file:Teuchos_RCP.hpp line:676
	M("Teuchos").def("rcpWithEmbeddedObjPostDestroy", [](class Teuchos::ParameterList * a0, const class Teuchos::RCP<class Teuchos::ParameterList> & a1) -> Teuchos::RCP<class Teuchos::ParameterList> { return Teuchos::rcpWithEmbeddedObjPostDestroy(a0, a1); }, "", pybind11::arg("p"), pybind11::arg("embedded"));
	M("Teuchos").def("rcpWithEmbeddedObjPostDestroy", (class Teuchos::RCP<class Teuchos::ParameterList> (*)(class Teuchos::ParameterList *, const class Teuchos::RCP<class Teuchos::ParameterList> &, bool)) &Teuchos::rcpWithEmbeddedObjPostDestroy<Teuchos::ParameterList,Teuchos::RCP<Teuchos::ParameterList>>, "C++: Teuchos::rcpWithEmbeddedObjPostDestroy(class Teuchos::ParameterList *, const class Teuchos::RCP<class Teuchos::ParameterList> &, bool) --> class Teuchos::RCP<class Teuchos::ParameterList>", pybind11::arg("p"), pybind11::arg("embedded"), pybind11::arg("owns_mem"));

	// Teuchos::rcpWithEmbeddedObjPostDestroy(const class Teuchos::ParameterList *, const class Teuchos::RCP<const class Teuchos::ParameterList> &, bool) file:Teuchos_RCP.hpp line:676
	M("Teuchos").def("rcpWithEmbeddedObjPostDestroy", [](const class Teuchos::ParameterList * a0, const class Teuchos::RCP<const class Teuchos::ParameterList> & a1) -> Teuchos::RCP<const class Teuchos::ParameterList> { return Teuchos::rcpWithEmbeddedObjPostDestroy(a0, a1); }, "", pybind11::arg("p"), pybind11::arg("embedded"));
	M("Teuchos").def("rcpWithEmbeddedObjPostDestroy", (class Teuchos::RCP<const class Teuchos::ParameterList> (*)(const class Teuchos::ParameterList *, const class Teuchos::RCP<const class Teuchos::ParameterList> &, bool)) &Teuchos::rcpWithEmbeddedObjPostDestroy<const Teuchos::ParameterList,Teuchos::RCP<const Teuchos::ParameterList>>, "C++: Teuchos::rcpWithEmbeddedObjPostDestroy(const class Teuchos::ParameterList *, const class Teuchos::RCP<const class Teuchos::ParameterList> &, bool) --> class Teuchos::RCP<const class Teuchos::ParameterList>", pybind11::arg("p"), pybind11::arg("embedded"), pybind11::arg("owns_mem"));

	// Teuchos::is_null(const class Teuchos::RCP<class Teuchos::XMLObjectImplem> &) file:Teuchos_RCP.hpp line:715
	M("Teuchos").def("is_null", (bool (*)(const class Teuchos::RCP<class Teuchos::XMLObjectImplem> &)) &Teuchos::is_null<Teuchos::XMLObjectImplem>, "C++: Teuchos::is_null(const class Teuchos::RCP<class Teuchos::XMLObjectImplem> &) --> bool", pybind11::arg("p"));

	// Teuchos::nonnull(const class Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > > &) file:Teuchos_RCP.hpp line:723
	M("Teuchos").def("nonnull", (bool (*)(const class Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > > &)) &Teuchos::nonnull<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>, "C++: Teuchos::nonnull(const class Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > > &) --> bool", pybind11::arg("p"));

	// Teuchos::nonnull(const class Teuchos::RCP<const class Teuchos::ParameterEntryValidator> &) file:Teuchos_RCP.hpp line:723
	M("Teuchos").def("nonnull", (bool (*)(const class Teuchos::RCP<const class Teuchos::ParameterEntryValidator> &)) &Teuchos::nonnull<const Teuchos::ParameterEntryValidator>, "C++: Teuchos::nonnull(const class Teuchos::RCP<const class Teuchos::ParameterEntryValidator> &) --> bool", pybind11::arg("p"));

	// Teuchos::rcp_static_cast(const class Teuchos::RCP<const struct ROL::TypeB::AlgorithmState<double> > &) file:Teuchos_RCP.hpp line:775
	M("Teuchos").def("rcp_static_cast", (class Teuchos::RCP<const struct ROL::TypeB::AlgorithmState<double> > (*)(const class Teuchos::RCP<const struct ROL::TypeB::AlgorithmState<double> > &)) &Teuchos::rcp_static_cast<const ROL::TypeB::AlgorithmState<double>,const ROL::TypeB::AlgorithmState<double>>, "C++: Teuchos::rcp_static_cast(const class Teuchos::RCP<const struct ROL::TypeB::AlgorithmState<double> > &) --> class Teuchos::RCP<const struct ROL::TypeB::AlgorithmState<double> >", pybind11::arg("p1"));

	// Teuchos::rcp_const_cast(const class Teuchos::RCP<const class ROL::Vector<double> > &) file:Teuchos_RCP.hpp line:786
	M("Teuchos").def("rcp_const_cast", (class Teuchos::RCP<class ROL::Vector<double> > (*)(const class Teuchos::RCP<const class ROL::Vector<double> > &)) &Teuchos::rcp_const_cast<ROL::Vector<double>,const ROL::Vector<double>>, "C++: Teuchos::rcp_const_cast(const class Teuchos::RCP<const class ROL::Vector<double> > &) --> class Teuchos::RCP<class ROL::Vector<double> >", pybind11::arg("p1"));

	{ // Teuchos::m_bad_cast file:Teuchos_dyn_cast.hpp line:60
		pybind11::class_<Teuchos::m_bad_cast, Teuchos::RCP<Teuchos::m_bad_cast>, PyCallBack_Teuchos_m_bad_cast, std::bad_cast> cl(M("Teuchos"), "m_bad_cast", "Exception class for bad cast.\n\n\nWe create this class so that we may throw a bad_cast when appropriate and\nstill use the TEUCHOS_TEST_FOR_EXCEPTION macro.  We recommend users try to catch a\nbad_cast.");
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("what_arg") );

		cl.def( pybind11::init( [](PyCallBack_Teuchos_m_bad_cast const &o){ return new PyCallBack_Teuchos_m_bad_cast(o); } ) );
		cl.def( pybind11::init( [](Teuchos::m_bad_cast const &o){ return new Teuchos::m_bad_cast(o); } ) );
		cl.def("what", (const char * (Teuchos::m_bad_cast::*)() const) &Teuchos::m_bad_cast::what, "C++: Teuchos::m_bad_cast::what() const --> const char *", pybind11::return_value_policy::automatic);
		cl.def("assign", (class Teuchos::m_bad_cast & (Teuchos::m_bad_cast::*)(const class Teuchos::m_bad_cast &)) &Teuchos::m_bad_cast::operator=, "C++: Teuchos::m_bad_cast::operator=(const class Teuchos::m_bad_cast &) --> class Teuchos::m_bad_cast &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	// Teuchos::dyn_cast_throw_exception(const std::string &, const std::string &, const std::string &) file:Teuchos_dyn_cast.hpp line:70
	M("Teuchos").def("dyn_cast_throw_exception", (void (*)(const std::string &, const std::string &, const std::string &)) &Teuchos::dyn_cast_throw_exception, "C++: Teuchos::dyn_cast_throw_exception(const std::string &, const std::string &, const std::string &) --> void", pybind11::arg("T_from"), pybind11::arg("T_from_concr"), pybind11::arg("T_to"));

	// Teuchos::ptrFromRef(class Teuchos::ParameterEntry &) file:Teuchos_PtrDecl.hpp line:298
	M("Teuchos").def("ptrFromRef", (class Teuchos::Ptr<class Teuchos::ParameterEntry> (*)(class Teuchos::ParameterEntry &)) &Teuchos::ptrFromRef<Teuchos::ParameterEntry>, "C++: Teuchos::ptrFromRef(class Teuchos::ParameterEntry &) --> class Teuchos::Ptr<class Teuchos::ParameterEntry>", pybind11::arg("arg"));

	// Teuchos::ptrFromRef(const class Teuchos::ParameterEntry &) file:Teuchos_PtrDecl.hpp line:298
	M("Teuchos").def("ptrFromRef", (class Teuchos::Ptr<const class Teuchos::ParameterEntry> (*)(const class Teuchos::ParameterEntry &)) &Teuchos::ptrFromRef<const Teuchos::ParameterEntry>, "C++: Teuchos::ptrFromRef(const class Teuchos::ParameterEntry &) --> class Teuchos::Ptr<const class Teuchos::ParameterEntry>", pybind11::arg("arg"));

	// Teuchos::rcpFromPtr(const class Teuchos::Ptr<class Teuchos::ParameterEntry> &) file:Teuchos_PtrDecl.hpp line:309
	M("Teuchos").def("rcpFromPtr", (class Teuchos::RCP<class Teuchos::ParameterEntry> (*)(const class Teuchos::Ptr<class Teuchos::ParameterEntry> &)) &Teuchos::rcpFromPtr<Teuchos::ParameterEntry>, "C++: Teuchos::rcpFromPtr(const class Teuchos::Ptr<class Teuchos::ParameterEntry> &) --> class Teuchos::RCP<class Teuchos::ParameterEntry>", pybind11::arg("ptr"));

	// Teuchos::rcpFromPtr(const class Teuchos::Ptr<const class Teuchos::ParameterEntry> &) file:Teuchos_PtrDecl.hpp line:309
	M("Teuchos").def("rcpFromPtr", (class Teuchos::RCP<const class Teuchos::ParameterEntry> (*)(const class Teuchos::Ptr<const class Teuchos::ParameterEntry> &)) &Teuchos::rcpFromPtr<const Teuchos::ParameterEntry>, "C++: Teuchos::rcpFromPtr(const class Teuchos::Ptr<const class Teuchos::ParameterEntry> &) --> class Teuchos::RCP<const class Teuchos::ParameterEntry>", pybind11::arg("ptr"));

	// Teuchos::is_null(const class Teuchos::Ptr<class Teuchos::ParameterEntry> &) file:Teuchos_PtrDecl.hpp line:355
	M("Teuchos").def("is_null", (bool (*)(const class Teuchos::Ptr<class Teuchos::ParameterEntry> &)) &Teuchos::is_null<Teuchos::ParameterEntry>, "C++: Teuchos::is_null(const class Teuchos::Ptr<class Teuchos::ParameterEntry> &) --> bool", pybind11::arg("p"));

	// Teuchos::is_null(const class Teuchos::Ptr<const class Teuchos::ParameterEntry> &) file:Teuchos_PtrDecl.hpp line:355
	M("Teuchos").def("is_null", (bool (*)(const class Teuchos::Ptr<const class Teuchos::ParameterEntry> &)) &Teuchos::is_null<const Teuchos::ParameterEntry>, "C++: Teuchos::is_null(const class Teuchos::Ptr<const class Teuchos::ParameterEntry> &) --> bool", pybind11::arg("p"));

}
