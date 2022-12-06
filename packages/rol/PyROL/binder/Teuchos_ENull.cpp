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
#include <ROL_Bundle_U.hpp>
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
#include <ROL_LineSearch_U_Types.hpp>
#include <ROL_LinearConstraint.hpp>
#include <ROL_LinearOperator.hpp>
#include <ROL_MINRES.hpp>
#include <ROL_NewtonKrylov_U.hpp>
#include <ROL_Newton_U.hpp>
#include <ROL_NonlinearCG.hpp>
#include <ROL_NonlinearCG_U.hpp>
#include <ROL_NullSpaceOperator.hpp>
#include <ROL_Objective.hpp>
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
#include <ROL_Step.hpp>
#include <ROL_Stream.hpp>
#include <ROL_TruncatedCG.hpp>
#include <ROL_TruncatedCG_U.hpp>
#include <ROL_TrustRegion.hpp>
#include <ROL_TrustRegionModel.hpp>
#include <ROL_TrustRegionModel_U.hpp>
#include <ROL_TrustRegionTypes.hpp>
#include <ROL_TrustRegion_U.hpp>
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
#include <ROL_lBFGS.hpp>
#include <ROL_lDFP.hpp>
#include <ROL_lSR1.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_BLAS_types.hpp>
#include <Teuchos_DataAccess.hpp>
#include <Teuchos_Dependency.hpp>
#include <Teuchos_DependencySheet.hpp>
#include <Teuchos_ENull.hpp>
#include <Teuchos_Exceptions.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_FilteredIterator.hpp>
#include <Teuchos_ParameterEntry.hpp>
#include <Teuchos_ParameterEntryValidator.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ParameterListModifier.hpp>
#include <Teuchos_Ptr.hpp>
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
#include <Teuchos_toString.hpp>
#include <cwchar>
#include <deque>
#include <functional>
#include <ios>
#include <iterator>
#include <locale>
#include <map>
#include <memory>
#include <ostream>
#include <random>
#include <set>
#include <sstream>
#include <sstream> // __str__
#include <stdexcept>
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

// Teuchos::ExceptionBase file:Teuchos_Exceptions.hpp line:57
struct PyCallBack_Teuchos_ExceptionBase : public Teuchos::ExceptionBase {
	using Teuchos::ExceptionBase::ExceptionBase;

	const char * what() const noexcept override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::ExceptionBase *>(this), "what");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const char *>::value) {
				static pybind11::detail::override_caster_t<const char *> caster;
				return pybind11::detail::cast_ref<const char *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const char *>(std::move(o));
		}
		return logic_error::what();
	}
};

// Teuchos::DuplicateOwningRCPError file:Teuchos_Exceptions.hpp line:69
struct PyCallBack_Teuchos_DuplicateOwningRCPError : public Teuchos::DuplicateOwningRCPError {
	using Teuchos::DuplicateOwningRCPError::DuplicateOwningRCPError;

	const char * what() const noexcept override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::DuplicateOwningRCPError *>(this), "what");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const char *>::value) {
				static pybind11::detail::override_caster_t<const char *> caster;
				return pybind11::detail::cast_ref<const char *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const char *>(std::move(o));
		}
		return logic_error::what();
	}
};

// Teuchos::NullReferenceError file:Teuchos_Exceptions.hpp line:77
struct PyCallBack_Teuchos_NullReferenceError : public Teuchos::NullReferenceError {
	using Teuchos::NullReferenceError::NullReferenceError;

	const char * what() const noexcept override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::NullReferenceError *>(this), "what");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const char *>::value) {
				static pybind11::detail::override_caster_t<const char *> caster;
				return pybind11::detail::cast_ref<const char *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const char *>(std::move(o));
		}
		return logic_error::what();
	}
};

// Teuchos::NonconstAccessError file:Teuchos_Exceptions.hpp line:85
struct PyCallBack_Teuchos_NonconstAccessError : public Teuchos::NonconstAccessError {
	using Teuchos::NonconstAccessError::NonconstAccessError;

	const char * what() const noexcept override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::NonconstAccessError *>(this), "what");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const char *>::value) {
				static pybind11::detail::override_caster_t<const char *> caster;
				return pybind11::detail::cast_ref<const char *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const char *>(std::move(o));
		}
		return logic_error::what();
	}
};

// Teuchos::RangeError file:Teuchos_Exceptions.hpp line:93
struct PyCallBack_Teuchos_RangeError : public Teuchos::RangeError {
	using Teuchos::RangeError::RangeError;

	const char * what() const noexcept override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::RangeError *>(this), "what");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const char *>::value) {
				static pybind11::detail::override_caster_t<const char *> caster;
				return pybind11::detail::cast_ref<const char *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const char *>(std::move(o));
		}
		return logic_error::what();
	}
};

// Teuchos::DanglingReferenceError file:Teuchos_Exceptions.hpp line:101
struct PyCallBack_Teuchos_DanglingReferenceError : public Teuchos::DanglingReferenceError {
	using Teuchos::DanglingReferenceError::DanglingReferenceError;

	const char * what() const noexcept override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::DanglingReferenceError *>(this), "what");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const char *>::value) {
				static pybind11::detail::override_caster_t<const char *> caster;
				return pybind11::detail::cast_ref<const char *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const char *>(std::move(o));
		}
		return logic_error::what();
	}
};

// Teuchos::IncompatibleIteratorsError file:Teuchos_Exceptions.hpp line:109
struct PyCallBack_Teuchos_IncompatibleIteratorsError : public Teuchos::IncompatibleIteratorsError {
	using Teuchos::IncompatibleIteratorsError::IncompatibleIteratorsError;

	const char * what() const noexcept override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::IncompatibleIteratorsError *>(this), "what");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const char *>::value) {
				static pybind11::detail::override_caster_t<const char *> caster;
				return pybind11::detail::cast_ref<const char *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const char *>(std::move(o));
		}
		return logic_error::what();
	}
};

// Teuchos::DuplicateParameterSublist file:Teuchos_Exceptions.hpp line:118
struct PyCallBack_Teuchos_DuplicateParameterSublist : public Teuchos::DuplicateParameterSublist {
	using Teuchos::DuplicateParameterSublist::DuplicateParameterSublist;

	const char * what() const noexcept override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::DuplicateParameterSublist *>(this), "what");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const char *>::value) {
				static pybind11::detail::override_caster_t<const char *> caster;
				return pybind11::detail::cast_ref<const char *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const char *>(std::move(o));
		}
		return logic_error::what();
	}
};

// Teuchos::DuplicateParameterEntryException file:Teuchos_Exceptions.hpp line:132
struct PyCallBack_Teuchos_DuplicateParameterEntryException : public Teuchos::DuplicateParameterEntryException {
	using Teuchos::DuplicateParameterEntryException::DuplicateParameterEntryException;

	const char * what() const noexcept override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::DuplicateParameterEntryException *>(this), "what");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const char *>::value) {
				static pybind11::detail::override_caster_t<const char *> caster;
				return pybind11::detail::cast_ref<const char *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const char *>(std::move(o));
		}
		return logic_error::what();
	}
};

// Teuchos::DuplicateParameterEntryIDException file:Teuchos_Exceptions.hpp line:145
struct PyCallBack_Teuchos_DuplicateParameterEntryIDException : public Teuchos::DuplicateParameterEntryIDException {
	using Teuchos::DuplicateParameterEntryIDException::DuplicateParameterEntryIDException;

	const char * what() const noexcept override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::DuplicateParameterEntryIDException *>(this), "what");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const char *>::value) {
				static pybind11::detail::override_caster_t<const char *> caster;
				return pybind11::detail::cast_ref<const char *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const char *>(std::move(o));
		}
		return logic_error::what();
	}
};

// Teuchos::DuplicateValidatorIDException file:Teuchos_Exceptions.hpp line:158
struct PyCallBack_Teuchos_DuplicateValidatorIDException : public Teuchos::DuplicateValidatorIDException {
	using Teuchos::DuplicateValidatorIDException::DuplicateValidatorIDException;

	const char * what() const noexcept override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::DuplicateValidatorIDException *>(this), "what");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const char *>::value) {
				static pybind11::detail::override_caster_t<const char *> caster;
				return pybind11::detail::cast_ref<const char *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const char *>(std::move(o));
		}
		return logic_error::what();
	}
};

// Teuchos::RCPNode file:Teuchos_RCPNode.hpp line:153
struct PyCallBack_Teuchos_RCPNode : public Teuchos::RCPNode {
	using Teuchos::RCPNode::RCPNode;

	bool is_valid_ptr() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::RCPNode *>(this), "is_valid_ptr");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"RCPNode::is_valid_ptr\"");
	}
	void delete_obj() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::RCPNode *>(this), "delete_obj");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"RCPNode::delete_obj\"");
	}
	void throw_invalid_obj_exception(const std::string & a0, const void * a1, const class Teuchos::RCPNode * a2, const void * a3) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::RCPNode *>(this), "throw_invalid_obj_exception");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"RCPNode::throw_invalid_obj_exception\"");
	}
	const std::string get_base_obj_type_name() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::RCPNode *>(this), "get_base_obj_type_name");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const std::string>::value) {
				static pybind11::detail::override_caster_t<const std::string> caster;
				return pybind11::detail::cast_ref<const std::string>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const std::string>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"RCPNode::get_base_obj_type_name\"");
	}
};

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

void bind_Teuchos_ENull(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	// Teuchos::ENull file:Teuchos_ENull.hpp line:54
	pybind11::enum_<Teuchos::ENull>(M("Teuchos"), "ENull", pybind11::arithmetic(), "Used to initialize a RCP object to NULL using an\n implicit conversion!\n\n \n\n ", pybind11::module_local())
		.value("null", Teuchos::null)
		.export_values();

;

	{ // Teuchos::ExceptionBase file:Teuchos_Exceptions.hpp line:57
		pybind11::class_<Teuchos::ExceptionBase, Teuchos::RCP<Teuchos::ExceptionBase>, PyCallBack_Teuchos_ExceptionBase, std::logic_error> cl(M("Teuchos"), "ExceptionBase", "Base exception class for Teuchos\n\n \n\n ", pybind11::module_local());
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("what_arg") );

		cl.def( pybind11::init( [](PyCallBack_Teuchos_ExceptionBase const &o){ return new PyCallBack_Teuchos_ExceptionBase(o); } ) );
		cl.def( pybind11::init( [](Teuchos::ExceptionBase const &o){ return new Teuchos::ExceptionBase(o); } ) );
		cl.def("assign", (class Teuchos::ExceptionBase & (Teuchos::ExceptionBase::*)(const class Teuchos::ExceptionBase &)) &Teuchos::ExceptionBase::operator=, "C++: Teuchos::ExceptionBase::operator=(const class Teuchos::ExceptionBase &) --> class Teuchos::ExceptionBase &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // Teuchos::DuplicateOwningRCPError file:Teuchos_Exceptions.hpp line:69
		pybind11::class_<Teuchos::DuplicateOwningRCPError, Teuchos::RCP<Teuchos::DuplicateOwningRCPError>, PyCallBack_Teuchos_DuplicateOwningRCPError, Teuchos::ExceptionBase> cl(M("Teuchos"), "DuplicateOwningRCPError", "Thrown if a duplicate owning RCP is creatd the the same object.\n\n \n\n ", pybind11::module_local());
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("what_arg") );

		cl.def("assign", (class Teuchos::DuplicateOwningRCPError & (Teuchos::DuplicateOwningRCPError::*)(const class Teuchos::DuplicateOwningRCPError &)) &Teuchos::DuplicateOwningRCPError::operator=, "C++: Teuchos::DuplicateOwningRCPError::operator=(const class Teuchos::DuplicateOwningRCPError &) --> class Teuchos::DuplicateOwningRCPError &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // Teuchos::NullReferenceError file:Teuchos_Exceptions.hpp line:77
		pybind11::class_<Teuchos::NullReferenceError, Teuchos::RCP<Teuchos::NullReferenceError>, PyCallBack_Teuchos_NullReferenceError, Teuchos::ExceptionBase> cl(M("Teuchos"), "NullReferenceError", "Null reference error exception class.\n\n \n\n ", pybind11::module_local());
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("what_arg") );

		cl.def( pybind11::init( [](PyCallBack_Teuchos_NullReferenceError const &o){ return new PyCallBack_Teuchos_NullReferenceError(o); } ) );
		cl.def( pybind11::init( [](Teuchos::NullReferenceError const &o){ return new Teuchos::NullReferenceError(o); } ) );
		cl.def("assign", (class Teuchos::NullReferenceError & (Teuchos::NullReferenceError::*)(const class Teuchos::NullReferenceError &)) &Teuchos::NullReferenceError::operator=, "C++: Teuchos::NullReferenceError::operator=(const class Teuchos::NullReferenceError &) --> class Teuchos::NullReferenceError &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // Teuchos::NonconstAccessError file:Teuchos_Exceptions.hpp line:85
		pybind11::class_<Teuchos::NonconstAccessError, Teuchos::RCP<Teuchos::NonconstAccessError>, PyCallBack_Teuchos_NonconstAccessError, Teuchos::ExceptionBase> cl(M("Teuchos"), "NonconstAccessError", "Null reference error exception class.\n\n \n\n ", pybind11::module_local());
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("what_arg") );

		cl.def("assign", (class Teuchos::NonconstAccessError & (Teuchos::NonconstAccessError::*)(const class Teuchos::NonconstAccessError &)) &Teuchos::NonconstAccessError::operator=, "C++: Teuchos::NonconstAccessError::operator=(const class Teuchos::NonconstAccessError &) --> class Teuchos::NonconstAccessError &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // Teuchos::RangeError file:Teuchos_Exceptions.hpp line:93
		pybind11::class_<Teuchos::RangeError, Teuchos::RCP<Teuchos::RangeError>, PyCallBack_Teuchos_RangeError, Teuchos::ExceptionBase> cl(M("Teuchos"), "RangeError", "Range error exception class.\n\n \n\n ", pybind11::module_local());
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("what_arg") );

		cl.def( pybind11::init( [](PyCallBack_Teuchos_RangeError const &o){ return new PyCallBack_Teuchos_RangeError(o); } ) );
		cl.def( pybind11::init( [](Teuchos::RangeError const &o){ return new Teuchos::RangeError(o); } ) );
		cl.def("assign", (class Teuchos::RangeError & (Teuchos::RangeError::*)(const class Teuchos::RangeError &)) &Teuchos::RangeError::operator=, "C++: Teuchos::RangeError::operator=(const class Teuchos::RangeError &) --> class Teuchos::RangeError &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // Teuchos::DanglingReferenceError file:Teuchos_Exceptions.hpp line:101
		pybind11::class_<Teuchos::DanglingReferenceError, Teuchos::RCP<Teuchos::DanglingReferenceError>, PyCallBack_Teuchos_DanglingReferenceError, Teuchos::ExceptionBase> cl(M("Teuchos"), "DanglingReferenceError", "Dangling reference error exception class.\n\n \n\n ", pybind11::module_local());
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("what_arg") );

		cl.def( pybind11::init( [](PyCallBack_Teuchos_DanglingReferenceError const &o){ return new PyCallBack_Teuchos_DanglingReferenceError(o); } ) );
		cl.def( pybind11::init( [](Teuchos::DanglingReferenceError const &o){ return new Teuchos::DanglingReferenceError(o); } ) );
		cl.def("assign", (class Teuchos::DanglingReferenceError & (Teuchos::DanglingReferenceError::*)(const class Teuchos::DanglingReferenceError &)) &Teuchos::DanglingReferenceError::operator=, "C++: Teuchos::DanglingReferenceError::operator=(const class Teuchos::DanglingReferenceError &) --> class Teuchos::DanglingReferenceError &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // Teuchos::IncompatibleIteratorsError file:Teuchos_Exceptions.hpp line:109
		pybind11::class_<Teuchos::IncompatibleIteratorsError, Teuchos::RCP<Teuchos::IncompatibleIteratorsError>, PyCallBack_Teuchos_IncompatibleIteratorsError, Teuchos::ExceptionBase> cl(M("Teuchos"), "IncompatibleIteratorsError", "Incompatiable iterators error exception class.\n\n \n\n ", pybind11::module_local());
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("what_arg") );

		cl.def("assign", (class Teuchos::IncompatibleIteratorsError & (Teuchos::IncompatibleIteratorsError::*)(const class Teuchos::IncompatibleIteratorsError &)) &Teuchos::IncompatibleIteratorsError::operator=, "C++: Teuchos::IncompatibleIteratorsError::operator=(const class Teuchos::IncompatibleIteratorsError &) --> class Teuchos::IncompatibleIteratorsError &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // Teuchos::DuplicateParameterSublist file:Teuchos_Exceptions.hpp line:118
		pybind11::class_<Teuchos::DuplicateParameterSublist, Teuchos::RCP<Teuchos::DuplicateParameterSublist>, PyCallBack_Teuchos_DuplicateParameterSublist, Teuchos::ExceptionBase> cl(M("Teuchos"), "DuplicateParameterSublist", "Optionally thrown when a sublist is set twice by either\n updateParametersFromXmlFile(), updateParametersFromXmlFileAndUpdate() or\n updateParametersFromXmlString()\n\n \n\n ", pybind11::module_local());
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("what_arg") );

		cl.def("assign", (class Teuchos::DuplicateParameterSublist & (Teuchos::DuplicateParameterSublist::*)(const class Teuchos::DuplicateParameterSublist &)) &Teuchos::DuplicateParameterSublist::operator=, "C++: Teuchos::DuplicateParameterSublist::operator=(const class Teuchos::DuplicateParameterSublist &) --> class Teuchos::DuplicateParameterSublist &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // Teuchos::DuplicateParameterEntryException file:Teuchos_Exceptions.hpp line:132
		pybind11::class_<Teuchos::DuplicateParameterEntryException, Teuchos::RCP<Teuchos::DuplicateParameterEntryException>, PyCallBack_Teuchos_DuplicateParameterEntryException, Teuchos::ExceptionBase> cl(M("Teuchos"), "DuplicateParameterEntryException", "Thrown when a Parameter Entry that is already being tracked\n is attempted to be inserted again into the masterParameterEntryMap\n and masterIDMap\n\n \n\n ", pybind11::module_local());
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("what_arg") );

		cl.def("assign", (class Teuchos::DuplicateParameterEntryException & (Teuchos::DuplicateParameterEntryException::*)(const class Teuchos::DuplicateParameterEntryException &)) &Teuchos::DuplicateParameterEntryException::operator=, "C++: Teuchos::DuplicateParameterEntryException::operator=(const class Teuchos::DuplicateParameterEntryException &) --> class Teuchos::DuplicateParameterEntryException &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // Teuchos::DuplicateParameterEntryIDException file:Teuchos_Exceptions.hpp line:145
		pybind11::class_<Teuchos::DuplicateParameterEntryIDException, Teuchos::RCP<Teuchos::DuplicateParameterEntryIDException>, PyCallBack_Teuchos_DuplicateParameterEntryIDException, Teuchos::ExceptionBase> cl(M("Teuchos"), "DuplicateParameterEntryIDException", "Thrown when a Parameter Entry ID that is already being used\n is attempted to be reused again.\n\n \n\n ", pybind11::module_local());
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("what_arg") );

		cl.def("assign", (class Teuchos::DuplicateParameterEntryIDException & (Teuchos::DuplicateParameterEntryIDException::*)(const class Teuchos::DuplicateParameterEntryIDException &)) &Teuchos::DuplicateParameterEntryIDException::operator=, "C++: Teuchos::DuplicateParameterEntryIDException::operator=(const class Teuchos::DuplicateParameterEntryIDException &) --> class Teuchos::DuplicateParameterEntryIDException &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // Teuchos::DuplicateValidatorIDException file:Teuchos_Exceptions.hpp line:158
		pybind11::class_<Teuchos::DuplicateValidatorIDException, Teuchos::RCP<Teuchos::DuplicateValidatorIDException>, PyCallBack_Teuchos_DuplicateValidatorIDException, Teuchos::ExceptionBase> cl(M("Teuchos"), "DuplicateValidatorIDException", "Thrown when a ParameterEntryValidatorID that\n is already being used is attempted to be reused again.\n\n \n\n ", pybind11::module_local());
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("what_arg") );

		cl.def("assign", (class Teuchos::DuplicateValidatorIDException & (Teuchos::DuplicateValidatorIDException::*)(const class Teuchos::DuplicateValidatorIDException &)) &Teuchos::DuplicateValidatorIDException::operator=, "C++: Teuchos::DuplicateValidatorIDException::operator=(const class Teuchos::DuplicateValidatorIDException &) --> class Teuchos::DuplicateValidatorIDException &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // Teuchos::ToStringTraits file:Teuchos_toString.hpp line:60
		pybind11::class_<Teuchos::ToStringTraits<int>, Teuchos::RCP<Teuchos::ToStringTraits<int>>> cl(M("Teuchos"), "ToStringTraits_int_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::ToStringTraits<int>(); } ) );
		cl.def_static("toString", (std::string (*)(const int &)) &Teuchos::ToStringTraits<int>::toString, "C++: Teuchos::ToStringTraits<int>::toString(const int &) --> std::string", pybind11::arg("t"));
	}
	// Teuchos::toString(const double &) file:Teuchos_toString.hpp line:82
	M("Teuchos").def("toString", (std::string (*)(const double &)) &Teuchos::toString<double>, "C++: Teuchos::toString(const double &) --> std::string", pybind11::arg("t"));

	// Teuchos::toString(const int &) file:Teuchos_toString.hpp line:82
	M("Teuchos").def("toString", (std::string (*)(const int &)) &Teuchos::toString<int>, "C++: Teuchos::toString(const int &) --> std::string", pybind11::arg("t"));

	// Teuchos::toString(const bool &) file:Teuchos_toString.hpp line:82
	M("Teuchos").def("toString", (std::string (*)(const bool &)) &Teuchos::toString<bool>, "C++: Teuchos::toString(const bool &) --> std::string", pybind11::arg("t"));

	// Teuchos::toString(const std::string &) file:Teuchos_toString.hpp line:82
	M("Teuchos").def("toString", (std::string (*)(const std::string &)) &Teuchos::toString<std::string>, "C++: Teuchos::toString(const std::string &) --> std::string", pybind11::arg("t"));

	{ // Teuchos::ToStringTraits file:Teuchos_toString.hpp line:90
		pybind11::class_<Teuchos::ToStringTraits<bool>, Teuchos::RCP<Teuchos::ToStringTraits<bool>>> cl(M("Teuchos"), "ToStringTraits_bool_t", "Specialization for bool. ", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::ToStringTraits<bool>(); } ) );
		cl.def_static("toString", (std::string (*)(const bool &)) &Teuchos::ToStringTraits<bool>::toString, "C++: Teuchos::ToStringTraits<bool>::toString(const bool &) --> std::string", pybind11::arg("t"));
	}
	{ // Teuchos::ToStringTraits file:Teuchos_toString.hpp line:103
		pybind11::class_<Teuchos::ToStringTraits<std::string>, Teuchos::RCP<Teuchos::ToStringTraits<std::string>>> cl(M("Teuchos"), "ToStringTraits_std_string_t", "Specialization for std::string. ", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::ToStringTraits<std::string>(); } ) );
		cl.def_static("toString", (std::string (*)(const std::string &)) &Teuchos::ToStringTraits<std::string >::toString, "C++: Teuchos::ToStringTraits<std::string >::toString(const std::string &) --> std::string", pybind11::arg("t"));
	}
	{ // Teuchos::ToStringTraits file:Teuchos_toString.hpp line:113
		pybind11::class_<Teuchos::ToStringTraits<double>, Teuchos::RCP<Teuchos::ToStringTraits<double>>> cl(M("Teuchos"), "ToStringTraits_double_t", "Specialization for double. ", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::ToStringTraits<double>(); } ) );
		cl.def_static("toString", (std::string (*)(const double &)) &Teuchos::ToStringTraits<double>::toString, "C++: Teuchos::ToStringTraits<double>::toString(const double &) --> std::string", pybind11::arg("t"));
	}
	{ // Teuchos::ToStringTraits file:Teuchos_toString.hpp line:149
		pybind11::class_<Teuchos::ToStringTraits<float>, Teuchos::RCP<Teuchos::ToStringTraits<float>>> cl(M("Teuchos"), "ToStringTraits_float_t", "Specialization for float. ", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::ToStringTraits<float>(); } ) );
		cl.def_static("toString", (std::string (*)(const float &)) &Teuchos::ToStringTraits<float>::toString, "C++: Teuchos::ToStringTraits<float>::toString(const float &) --> std::string", pybind11::arg("t"));
	}
	// Teuchos::EPrePostDestruction file:Teuchos_RCPNode.hpp line:79
	pybind11::enum_<Teuchos::EPrePostDestruction>(M("Teuchos"), "EPrePostDestruction", pybind11::arithmetic(), "Used to specify a pre or post destruction of extra data\n\n \n\n ", pybind11::module_local())
		.value("PRE_DESTROY", Teuchos::PRE_DESTROY)
		.value("POST_DESTROY", Teuchos::POST_DESTROY)
		.export_values();

;

	// Teuchos::ERCPStrength file:Teuchos_RCPNode.hpp line:85
	pybind11::enum_<Teuchos::ERCPStrength>(M("Teuchos"), "ERCPStrength", pybind11::arithmetic(), "Used to specify if the pointer is weak or strong.\n\n \n\n ", pybind11::module_local())
		.value("RCP_STRONG", Teuchos::RCP_STRONG)
		.value("RCP_WEAK", Teuchos::RCP_WEAK)
		.export_values();

;

	// Teuchos::ERCPNodeLookup file:Teuchos_RCPNode.hpp line:91
	pybind11::enum_<Teuchos::ERCPNodeLookup>(M("Teuchos"), "ERCPNodeLookup", pybind11::arithmetic(), "Used to determine if RCPNode lookup is performed or not.\n\n \n\n ", pybind11::module_local())
		.value("RCP_ENABLE_NODE_LOOKUP", Teuchos::RCP_ENABLE_NODE_LOOKUP)
		.value("RCP_DISABLE_NODE_LOOKUP", Teuchos::RCP_DISABLE_NODE_LOOKUP)
		.export_values();

;

	// Teuchos::debugAssertStrength(enum Teuchos::ERCPStrength) file:Teuchos_RCPNode.hpp line:94
	M("Teuchos").def("debugAssertStrength", (void (*)(enum Teuchos::ERCPStrength)) &Teuchos::debugAssertStrength, ". \n\nC++: Teuchos::debugAssertStrength(enum Teuchos::ERCPStrength) --> void", pybind11::arg("strength"));

	{ // Teuchos::ToStringTraits file:Teuchos_RCPNode.hpp line:119
		pybind11::class_<Teuchos::ToStringTraits<Teuchos::ERCPStrength>, Teuchos::RCP<Teuchos::ToStringTraits<Teuchos::ERCPStrength>>> cl(M("Teuchos"), "ToStringTraits_Teuchos_ERCPStrength_t", "Traits class specialization for toString(...) function for\n converting from ERCPStrength to std::string.\n\n \n\n ", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::ToStringTraits<Teuchos::ERCPStrength>(); } ) );
		cl.def_static("toString", (std::string (*)(const enum Teuchos::ERCPStrength &)) &Teuchos::ToStringTraits<Teuchos::ERCPStrength>::toString, "C++: Teuchos::ToStringTraits<Teuchos::ERCPStrength>::toString(const enum Teuchos::ERCPStrength &) --> std::string", pybind11::arg("t"));
	}
	{ // Teuchos::RCPNode file:Teuchos_RCPNode.hpp line:153
		pybind11::class_<Teuchos::RCPNode, Teuchos::RCP<Teuchos::RCPNode>, PyCallBack_Teuchos_RCPNode> cl(M("Teuchos"), "RCPNode", "Node class to keep track of address and the reference count for a\n reference-counted utility class and delete the object.\n\n This is not a general user-level class.  This is used in the implementation\n of all of the reference-counting utility classes.\n\n NOTE: The reference counts all start a 0 so the client (i.e. RCPNodeHandle)\n must increment them from 0 after creation.\n\n \n\n ", pybind11::module_local());
		cl.def( pybind11::init<bool>(), pybind11::arg("has_ownership_in") );

		cl.def("attemptIncrementStrongCountFromNonZeroValue", (bool (Teuchos::RCPNode::*)()) &Teuchos::RCPNode::attemptIncrementStrongCountFromNonZeroValue, "attemptIncrementStrongCountFromNonZeroValue() supports weak\n to strong conversion but this is forward looking code.\n\nC++: Teuchos::RCPNode::attemptIncrementStrongCountFromNonZeroValue() --> bool");
		cl.def("strong_count", (int (Teuchos::RCPNode::*)() const) &Teuchos::RCPNode::strong_count, ". \n\nC++: Teuchos::RCPNode::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCPNode::*)() const) &Teuchos::RCPNode::weak_count, ". \n\nC++: Teuchos::RCPNode::weak_count() const --> int");
		cl.def("incr_count", (void (Teuchos::RCPNode::*)(const enum Teuchos::ERCPStrength)) &Teuchos::RCPNode::incr_count, ". \n\nC++: Teuchos::RCPNode::incr_count(const enum Teuchos::ERCPStrength) --> void", pybind11::arg("strength"));
		cl.def("deincr_count", (int (Teuchos::RCPNode::*)(const enum Teuchos::ERCPStrength)) &Teuchos::RCPNode::deincr_count, ". \n\nC++: Teuchos::RCPNode::deincr_count(const enum Teuchos::ERCPStrength) --> int", pybind11::arg("strength"));
		cl.def("has_ownership", (void (Teuchos::RCPNode::*)(bool)) &Teuchos::RCPNode::has_ownership, ". \n\nC++: Teuchos::RCPNode::has_ownership(bool) --> void", pybind11::arg("has_ownership_in"));
		cl.def("has_ownership", (bool (Teuchos::RCPNode::*)() const) &Teuchos::RCPNode::has_ownership, ". \n\nC++: Teuchos::RCPNode::has_ownership() const --> bool");
		cl.def("set_extra_data", (void (Teuchos::RCPNode::*)(const class Teuchos::any &, const std::string &, enum Teuchos::EPrePostDestruction, bool)) &Teuchos::RCPNode::set_extra_data, ". \n\nC++: Teuchos::RCPNode::set_extra_data(const class Teuchos::any &, const std::string &, enum Teuchos::EPrePostDestruction, bool) --> void", pybind11::arg("extra_data"), pybind11::arg("name"), pybind11::arg("destroy_when"), pybind11::arg("force_unique"));
		cl.def("get_extra_data", (class Teuchos::any & (Teuchos::RCPNode::*)(const std::string &, const std::string &)) &Teuchos::RCPNode::get_extra_data, ". \n\nC++: Teuchos::RCPNode::get_extra_data(const std::string &, const std::string &) --> class Teuchos::any &", pybind11::return_value_policy::automatic, pybind11::arg("type_name"), pybind11::arg("name"));
		cl.def("get_optional_extra_data", (class Teuchos::any * (Teuchos::RCPNode::*)(const std::string &, const std::string &)) &Teuchos::RCPNode::get_optional_extra_data, ". \n\nC++: Teuchos::RCPNode::get_optional_extra_data(const std::string &, const std::string &) --> class Teuchos::any *", pybind11::return_value_policy::automatic, pybind11::arg("type_name"), pybind11::arg("name"));
		cl.def("is_valid_ptr", (bool (Teuchos::RCPNode::*)() const) &Teuchos::RCPNode::is_valid_ptr, ". \n\nC++: Teuchos::RCPNode::is_valid_ptr() const --> bool");
		cl.def("delete_obj", (void (Teuchos::RCPNode::*)()) &Teuchos::RCPNode::delete_obj, ". \n\nC++: Teuchos::RCPNode::delete_obj() --> void");
		cl.def("throw_invalid_obj_exception", (void (Teuchos::RCPNode::*)(const std::string &, const void *, const class Teuchos::RCPNode *, const void *) const) &Teuchos::RCPNode::throw_invalid_obj_exception, ". \n\nC++: Teuchos::RCPNode::throw_invalid_obj_exception(const std::string &, const void *, const class Teuchos::RCPNode *, const void *) const --> void", pybind11::arg("rcp_type_name"), pybind11::arg("rcp_ptr"), pybind11::arg("rcp_node_ptr"), pybind11::arg("rcp_obj_ptr"));
		cl.def("get_base_obj_type_name", (const std::string (Teuchos::RCPNode::*)() const) &Teuchos::RCPNode::get_base_obj_type_name, ". \n\nC++: Teuchos::RCPNode::get_base_obj_type_name() const --> const std::string");
	}
	// Teuchos::throw_null_ptr_error(const std::string &) file:Teuchos_RCPNode.hpp line:334
	M("Teuchos").def("throw_null_ptr_error", (void (*)(const std::string &)) &Teuchos::throw_null_ptr_error, "Throw that a pointer passed into an RCP object is null.\n\n \n\n \n\nC++: Teuchos::throw_null_ptr_error(const std::string &) --> void", pybind11::arg("type_name"));

	{ // Teuchos::RCPNodeTracer file:Teuchos_RCPNode.hpp line:368
		pybind11::class_<Teuchos::RCPNodeTracer, Teuchos::RCP<Teuchos::RCPNodeTracer>> cl(M("Teuchos"), "RCPNodeTracer", "Debug-mode RCPNode tracing class.\n\n This is a static class that is used to trace all RCP nodes that are created\n and destroyed and to look-up RCPNodes given an an object's address.  This\n database is used for several different types of debug-mode runtime checking\n including a) the detection of cicular references, b) detecting the creation\n of duplicate owning RCPNode objects for the same reference-counted object,\n and c) to create weak RCP objects for existing RCPNode objects.\n\n This is primarily an internal implementation class but there are a few\n functions (maked as such below) that can be called by general users to turn\n on and off node tracing and to print the active RCPNode objects at any\n time.\n\n \n\n ", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::RCPNodeTracer(); } ) );
		cl.def_static("isTracingActiveRCPNodes", (bool (*)()) &Teuchos::RCPNodeTracer::isTracingActiveRCPNodes, "Return if we are tracing active nodes or not.\n\n NOTE: This will always return false when TEUCHOS_DEBUG is\n not defined.\n\nC++: Teuchos::RCPNodeTracer::isTracingActiveRCPNodes() --> bool");
		cl.def_static("numActiveRCPNodes", (int (*)()) &Teuchos::RCPNodeTracer::numActiveRCPNodes, "Print the number of active RCPNode objects currently being\n tracked.\n\nC++: Teuchos::RCPNodeTracer::numActiveRCPNodes() --> int");
		cl.def_static("getRCPNodeStatistics", (struct Teuchos::RCPNodeTracer::RCPNodeStatistics (*)()) &Teuchos::RCPNodeTracer::getRCPNodeStatistics, "Return the statistics on RCPNode allocations. \n\nC++: Teuchos::RCPNodeTracer::getRCPNodeStatistics() --> struct Teuchos::RCPNodeTracer::RCPNodeStatistics");
		cl.def_static("setPrintRCPNodeStatisticsOnExit", (void (*)(bool)) &Teuchos::RCPNodeTracer::setPrintRCPNodeStatisticsOnExit, "Set if RCPNode usage statistics will be printed when the program\n ends or not.\n\nC++: Teuchos::RCPNodeTracer::setPrintRCPNodeStatisticsOnExit(bool) --> void", pybind11::arg("printRCPNodeStatisticsOnExit"));
		cl.def_static("getPrintRCPNodeStatisticsOnExit", (bool (*)()) &Teuchos::RCPNodeTracer::getPrintRCPNodeStatisticsOnExit, "Return if RCPNode usage statistics will be printed when the\n program ends or not.\n\nC++: Teuchos::RCPNodeTracer::getPrintRCPNodeStatisticsOnExit() --> bool");
		cl.def_static("setPrintActiveRcpNodesOnExit", (void (*)(bool)) &Teuchos::RCPNodeTracer::setPrintActiveRcpNodesOnExit, "Set if printActiveRCPNodes() is called on exit from the\n program.\n\nC++: Teuchos::RCPNodeTracer::setPrintActiveRcpNodesOnExit(bool) --> void", pybind11::arg("printActiveRcpNodesOnExit"));
		cl.def_static("getPrintActiveRcpNodesOnExit", (bool (*)()) &Teuchos::RCPNodeTracer::getPrintActiveRcpNodesOnExit, "Return if printActiveRCPNodes() is called on exit from the\n program.\n\nC++: Teuchos::RCPNodeTracer::getPrintActiveRcpNodesOnExit() --> bool");
		cl.def_static("addNewRCPNode", (void (*)(class Teuchos::RCPNode *, const std::string &)) &Teuchos::RCPNodeTracer::addNewRCPNode, "Add new RCPNode to the global list.\n\n Only gets called when RCPNode tracing has been activated.\n\nC++: Teuchos::RCPNodeTracer::addNewRCPNode(class Teuchos::RCPNode *, const std::string &) --> void", pybind11::arg("rcp_node"), pybind11::arg("info"));
		cl.def_static("removeRCPNode", (void (*)(class Teuchos::RCPNode *)) &Teuchos::RCPNodeTracer::removeRCPNode, "Remove an RCPNode from global list.\n\n Always gets called in a debug build (TEUCHOS_DEBUG defined) when\n node tracing is enabled.\n\nC++: Teuchos::RCPNodeTracer::removeRCPNode(class Teuchos::RCPNode *) --> void", pybind11::arg("rcp_node"));
		cl.def_static("getExistingRCPNodeGivenLookupKey", (class Teuchos::RCPNode * (*)(const void *)) &Teuchos::RCPNodeTracer::getExistingRCPNodeGivenLookupKey, "Return a raw pointer to an existing owning RCPNode given its\n lookup key.\n\n \n returnVal != 0 if an owning RCPNode exists, 0\n otherwsise.\n\nC++: Teuchos::RCPNodeTracer::getExistingRCPNodeGivenLookupKey(const void *) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("lookupKey"));
		cl.def_static("getActiveRCPNodeHeaderString", (std::string (*)()) &Teuchos::RCPNodeTracer::getActiveRCPNodeHeaderString, "Header string used in printActiveRCPNodes(). \n\nC++: Teuchos::RCPNodeTracer::getActiveRCPNodeHeaderString() --> std::string");
		cl.def_static("getCommonDebugNotesString", (std::string (*)()) &Teuchos::RCPNodeTracer::getCommonDebugNotesString, "Common error message string on how to debug RCPNode problems. \n\nC++: Teuchos::RCPNodeTracer::getCommonDebugNotesString() --> std::string");

		{ // Teuchos::RCPNodeTracer::RCPNodeStatistics file:Teuchos_RCPNode.hpp line:375
			auto & enclosing_class = cl;
			pybind11::class_<Teuchos::RCPNodeTracer::RCPNodeStatistics, Teuchos::RCP<Teuchos::RCPNodeTracer::RCPNodeStatistics>> cl(enclosing_class, "RCPNodeStatistics", "RCP statistics struct. ", pybind11::module_local());
			cl.def( pybind11::init( [](){ return new Teuchos::RCPNodeTracer::RCPNodeStatistics(); } ) );
			cl.def_readwrite("maxNumRCPNodes", &Teuchos::RCPNodeTracer::RCPNodeStatistics::maxNumRCPNodes);
			cl.def_readwrite("totalNumRCPNodeAllocations", &Teuchos::RCPNodeTracer::RCPNodeStatistics::totalNumRCPNodeAllocations);
			cl.def_readwrite("totalNumRCPNodeDeletions", &Teuchos::RCPNodeTracer::RCPNodeStatistics::totalNumRCPNodeDeletions);
		}

	}
	{ // Teuchos::ActiveRCPNodesSetup file:Teuchos_RCPNode.hpp line:703
		pybind11::class_<Teuchos::ActiveRCPNodesSetup, Teuchos::RCP<Teuchos::ActiveRCPNodesSetup>> cl(M("Teuchos"), "ActiveRCPNodesSetup", "Sets up node tracing and prints remaining RCPNodes on destruction.\n\n This class is used by automataic code that sets up support for RCPNode\n tracing and for printing of remaining nodes on destruction.\n\n \n\n ", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::ActiveRCPNodesSetup(); } ) );
		cl.def( pybind11::init( [](Teuchos::ActiveRCPNodesSetup const &o){ return new Teuchos::ActiveRCPNodesSetup(o); } ) );
		cl.def("foo", (void (Teuchos::ActiveRCPNodesSetup::*)()) &Teuchos::ActiveRCPNodesSetup::foo, ". \n\nC++: Teuchos::ActiveRCPNodesSetup::foo() --> void");
	}
	{ // Teuchos::RCPNodeHandle file:Teuchos_RCPNode.hpp line:748
		pybind11::class_<Teuchos::RCPNodeHandle, Teuchos::RCP<Teuchos::RCPNodeHandle>> cl(M("Teuchos"), "RCPNodeHandle", "Handle class that manages the RCPNode's reference counting.\n\n \n This class is not intended for Teuchos users.  It\n   is an implementation detail of Teuchos' reference-counting\n   \"smart\" pointer (RCP) and array (ArrayRCP) classes.\n\n NOTE: I (Ross Bartlett) am not generally a big fan of handle classes and\n greatly prefer smart pointers.  However, this is one case where a handle\n class makes sense.  First, I want special behavior in some functions when\n the wrapped RCPNode pointer is null.  Second, I can't use one of the\n smart-pointer classes because this class is used to implement all of those\n smart-pointer classes!\n\n \n\n ", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::RCPNodeHandle(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("null_arg") );

		cl.def( pybind11::init( [](class Teuchos::RCPNode * a0){ return new Teuchos::RCPNodeHandle(a0); } ), "doc" , pybind11::arg("node"));
		cl.def( pybind11::init( [](class Teuchos::RCPNode * a0, enum Teuchos::ERCPStrength const & a1){ return new Teuchos::RCPNodeHandle(a0, a1); } ), "doc" , pybind11::arg("node"), pybind11::arg("strength_in"));
		cl.def( pybind11::init<class Teuchos::RCPNode *, enum Teuchos::ERCPStrength, bool>(), pybind11::arg("node"), pybind11::arg("strength_in"), pybind11::arg("newNode") );

		cl.def( pybind11::init( [](Teuchos::RCPNodeHandle const &o){ return new Teuchos::RCPNodeHandle(o); } ) );
		cl.def("swap", (void (Teuchos::RCPNodeHandle::*)(class Teuchos::RCPNodeHandle &)) &Teuchos::RCPNodeHandle::swap, "Swap the contents of  with \n\nC++: Teuchos::RCPNodeHandle::swap(class Teuchos::RCPNodeHandle &) --> void", pybind11::arg("node_ref"));
		cl.def("assign", (class Teuchos::RCPNodeHandle & (Teuchos::RCPNodeHandle::*)(enum Teuchos::ENull)) &Teuchos::RCPNodeHandle::operator=, "Null assignment.\n\n This method satisfies the strong exception guarantee: It either\n returns successfully, or throws an exception without modifying\n any user-visible state.\n\nC++: Teuchos::RCPNodeHandle::operator=(enum Teuchos::ENull) --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("assign", (class Teuchos::RCPNodeHandle & (Teuchos::RCPNodeHandle::*)(const class Teuchos::RCPNodeHandle &)) &Teuchos::RCPNodeHandle::operator=, "Copy assignment operator.\n\n This method satisfies the strong exception guarantee: It either\n returns successfully, or throws an exception without modifying\n any user-visible state.\n\nC++: Teuchos::RCPNodeHandle::operator=(const class Teuchos::RCPNodeHandle &) --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic, pybind11::arg("node_ref"));
		cl.def("create_strong_lock", (class Teuchos::RCPNodeHandle (Teuchos::RCPNodeHandle::*)() const) &Teuchos::RCPNodeHandle::create_strong_lock, "Return a strong handle if possible using thread safe atomics\n\nC++: Teuchos::RCPNodeHandle::create_strong_lock() const --> class Teuchos::RCPNodeHandle");
		cl.def("create_weak", (class Teuchos::RCPNodeHandle (Teuchos::RCPNodeHandle::*)() const) &Teuchos::RCPNodeHandle::create_weak, "Return a weak handle.\n\nC++: Teuchos::RCPNodeHandle::create_weak() const --> class Teuchos::RCPNodeHandle");
		cl.def("create_strong", (class Teuchos::RCPNodeHandle (Teuchos::RCPNodeHandle::*)() const) &Teuchos::RCPNodeHandle::create_strong, "Return a strong handle.\n\nC++: Teuchos::RCPNodeHandle::create_strong() const --> class Teuchos::RCPNodeHandle");
		cl.def("node_ptr", (class Teuchos::RCPNode * (Teuchos::RCPNodeHandle::*)() const) &Teuchos::RCPNodeHandle::node_ptr, "Return a pointer to the underlying RCPNode.\n\nC++: Teuchos::RCPNodeHandle::node_ptr() const --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic);
		cl.def("is_node_null", (bool (Teuchos::RCPNodeHandle::*)() const) &Teuchos::RCPNodeHandle::is_node_null, "Whether the underlying RCPNode is NULL.\n\nC++: Teuchos::RCPNodeHandle::is_node_null() const --> bool");
		cl.def("is_valid_ptr", (bool (Teuchos::RCPNodeHandle::*)() const) &Teuchos::RCPNodeHandle::is_valid_ptr, "Whether the underlying pointer is valid.\n\n \n NULL is a valid pointer; this method returns true in that case.\n\nC++: Teuchos::RCPNodeHandle::is_valid_ptr() const --> bool");
		cl.def("same_node", (bool (Teuchos::RCPNodeHandle::*)(const class Teuchos::RCPNodeHandle &) const) &Teuchos::RCPNodeHandle::same_node, "Whether the RCPNode for which  is a handle is the\n   same RCPNode as this object's RCPNode.\n\nC++: Teuchos::RCPNodeHandle::same_node(const class Teuchos::RCPNodeHandle &) const --> bool", pybind11::arg("node2"));
		cl.def("strong_count", (int (Teuchos::RCPNodeHandle::*)() const) &Teuchos::RCPNodeHandle::strong_count, "The strong count for this RCPNode, or 0 if the node is NULL.\n\nC++: Teuchos::RCPNodeHandle::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCPNodeHandle::*)() const) &Teuchos::RCPNodeHandle::weak_count, "The weak count for this RCPNode, or 0 if the node is NULL.\n\nC++: Teuchos::RCPNodeHandle::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCPNodeHandle::*)() const) &Teuchos::RCPNodeHandle::total_count, "The sum of the weak and string counts.\n\nC++: Teuchos::RCPNodeHandle::total_count() const --> int");
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCPNodeHandle::*)() const) &Teuchos::RCPNodeHandle::strength, "The strength of this handle.\n\nC++: Teuchos::RCPNodeHandle::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("has_ownership", (void (Teuchos::RCPNodeHandle::*)(bool)) &Teuchos::RCPNodeHandle::has_ownership, ". \n\nC++: Teuchos::RCPNodeHandle::has_ownership(bool) --> void", pybind11::arg("has_ownership_in"));
		cl.def("has_ownership", (bool (Teuchos::RCPNodeHandle::*)() const) &Teuchos::RCPNodeHandle::has_ownership, ". \n\nC++: Teuchos::RCPNodeHandle::has_ownership() const --> bool");
		cl.def("set_extra_data", (void (Teuchos::RCPNodeHandle::*)(const class Teuchos::any &, const std::string &, enum Teuchos::EPrePostDestruction, bool)) &Teuchos::RCPNodeHandle::set_extra_data, ". \n\nC++: Teuchos::RCPNodeHandle::set_extra_data(const class Teuchos::any &, const std::string &, enum Teuchos::EPrePostDestruction, bool) --> void", pybind11::arg("extra_data"), pybind11::arg("name"), pybind11::arg("destroy_when"), pybind11::arg("force_unique"));
		cl.def("get_extra_data", (class Teuchos::any & (Teuchos::RCPNodeHandle::*)(const std::string &, const std::string &)) &Teuchos::RCPNodeHandle::get_extra_data, ". \n\nC++: Teuchos::RCPNodeHandle::get_extra_data(const std::string &, const std::string &) --> class Teuchos::any &", pybind11::return_value_policy::automatic, pybind11::arg("type_name"), pybind11::arg("name"));
		cl.def("get_optional_extra_data", (class Teuchos::any * (Teuchos::RCPNodeHandle::*)(const std::string &, const std::string &)) &Teuchos::RCPNodeHandle::get_optional_extra_data, ". \n\nC++: Teuchos::RCPNodeHandle::get_optional_extra_data(const std::string &, const std::string &) --> class Teuchos::any *", pybind11::return_value_policy::automatic, pybind11::arg("type_name"), pybind11::arg("name"));
		cl.def("debug_assert_not_null", (void (Teuchos::RCPNodeHandle::*)() const) &Teuchos::RCPNodeHandle::debug_assert_not_null, ". \n\nC++: Teuchos::RCPNodeHandle::debug_assert_not_null() const --> void");

		cl.def("__str__", [](Teuchos::RCPNodeHandle const &o) -> std::string { std::ostringstream s; Teuchos::operator<<(s, o); return s.str(); } );
	}
	{ // Teuchos::RCPNodeThrowDeleter file:Teuchos_RCPNode.hpp line:1110
		pybind11::class_<Teuchos::RCPNodeThrowDeleter, Teuchos::RCP<Teuchos::RCPNodeThrowDeleter>> cl(M("Teuchos"), "RCPNodeThrowDeleter", "Deletes a (non-owning) RCPNode but not it's underlying object in\n case of a throw.\n\n This class is used in contexts where RCPNodeTracer::addNewRCPNode(...)\n might thrown an exception for a duplicate node being added.  The assumption\n is that there must already be an owning (or non-owning) RCP object that\n will delete the underlying object and therefore this class should *not*\n call delete_obj()!", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::RCPNode *>(), pybind11::arg("node") );

		cl.def("get", (class Teuchos::RCPNode * (Teuchos::RCPNodeThrowDeleter::*)() const) &Teuchos::RCPNodeThrowDeleter::get, ". \n\nC++: Teuchos::RCPNodeThrowDeleter::get() const --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic);
		cl.def("release", (void (Teuchos::RCPNodeThrowDeleter::*)()) &Teuchos::RCPNodeThrowDeleter::release, "Releaes the RCPNode pointer before the destructor is called. \n\nC++: Teuchos::RCPNodeThrowDeleter::release() --> void");
	}
	{ // Teuchos::Ptr file:Teuchos_PtrDecl.hpp line:104
		pybind11::class_<Teuchos::Ptr<const Teuchos::ParameterEntry>, Teuchos::RCP<Teuchos::Ptr<const Teuchos::ParameterEntry>>> cl(M("Teuchos"), "Ptr_const_Teuchos_ParameterEntry_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::Ptr<const Teuchos::ParameterEntry>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("null_in") );

		cl.def( pybind11::init<const class Teuchos::ParameterEntry *>(), pybind11::arg("ptr_in") );

		cl.def( pybind11::init( [](Teuchos::Ptr<const Teuchos::ParameterEntry> const &o){ return new Teuchos::Ptr<const Teuchos::ParameterEntry>(o); } ) );
		cl.def("assign", (class Teuchos::Ptr<const class Teuchos::ParameterEntry> & (Teuchos::Ptr<const Teuchos::ParameterEntry>::*)(const class Teuchos::Ptr<const class Teuchos::ParameterEntry> &)) &Teuchos::Ptr<const Teuchos::ParameterEntry>::operator=, "C++: Teuchos::Ptr<const Teuchos::ParameterEntry>::operator=(const class Teuchos::Ptr<const class Teuchos::ParameterEntry> &) --> class Teuchos::Ptr<const class Teuchos::ParameterEntry> &", pybind11::return_value_policy::automatic, pybind11::arg("ptr"));
		cl.def("arrow", (const class Teuchos::ParameterEntry * (Teuchos::Ptr<const Teuchos::ParameterEntry>::*)() const) &Teuchos::Ptr<const Teuchos::ParameterEntry>::operator->, "C++: Teuchos::Ptr<const Teuchos::ParameterEntry>::operator->() const --> const class Teuchos::ParameterEntry *", pybind11::return_value_policy::automatic);
		cl.def("dereference", (const class Teuchos::ParameterEntry & (Teuchos::Ptr<const Teuchos::ParameterEntry>::*)() const) &Teuchos::Ptr<const Teuchos::ParameterEntry>::operator*, "C++: Teuchos::Ptr<const Teuchos::ParameterEntry>::operator*() const --> const class Teuchos::ParameterEntry &", pybind11::return_value_policy::automatic);
		cl.def("get", (const class Teuchos::ParameterEntry * (Teuchos::Ptr<const Teuchos::ParameterEntry>::*)() const) &Teuchos::Ptr<const Teuchos::ParameterEntry>::get, "C++: Teuchos::Ptr<const Teuchos::ParameterEntry>::get() const --> const class Teuchos::ParameterEntry *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (const class Teuchos::ParameterEntry * (Teuchos::Ptr<const Teuchos::ParameterEntry>::*)() const) &Teuchos::Ptr<const Teuchos::ParameterEntry>::getRawPtr, "C++: Teuchos::Ptr<const Teuchos::ParameterEntry>::getRawPtr() const --> const class Teuchos::ParameterEntry *", pybind11::return_value_policy::automatic);
		cl.def("is_null", (bool (Teuchos::Ptr<const Teuchos::ParameterEntry>::*)() const) &Teuchos::Ptr<const Teuchos::ParameterEntry>::is_null, "C++: Teuchos::Ptr<const Teuchos::ParameterEntry>::is_null() const --> bool");
		cl.def("assert_not_null", (const class Teuchos::Ptr<const class Teuchos::ParameterEntry> & (Teuchos::Ptr<const Teuchos::ParameterEntry>::*)() const) &Teuchos::Ptr<const Teuchos::ParameterEntry>::assert_not_null, "C++: Teuchos::Ptr<const Teuchos::ParameterEntry>::assert_not_null() const --> const class Teuchos::Ptr<const class Teuchos::ParameterEntry> &", pybind11::return_value_policy::automatic);
		cl.def("ptr", (const class Teuchos::Ptr<const class Teuchos::ParameterEntry> (Teuchos::Ptr<const Teuchos::ParameterEntry>::*)() const) &Teuchos::Ptr<const Teuchos::ParameterEntry>::ptr, "C++: Teuchos::Ptr<const Teuchos::ParameterEntry>::ptr() const --> const class Teuchos::Ptr<const class Teuchos::ParameterEntry>");
		cl.def("getConst", (class Teuchos::Ptr<const class Teuchos::ParameterEntry> (Teuchos::Ptr<const Teuchos::ParameterEntry>::*)() const) &Teuchos::Ptr<const Teuchos::ParameterEntry>::getConst, "C++: Teuchos::Ptr<const Teuchos::ParameterEntry>::getConst() const --> class Teuchos::Ptr<const class Teuchos::ParameterEntry>");

		cl.def("__str__", [](Teuchos::Ptr<const Teuchos::ParameterEntry> const &o) -> std::string { std::ostringstream s; Teuchos::operator<<(s, o); return s.str(); } );
	}
	{ // Teuchos::Ptr file:Teuchos_PtrDecl.hpp line:104
		pybind11::class_<Teuchos::Ptr<Teuchos::ParameterEntry>, Teuchos::RCP<Teuchos::Ptr<Teuchos::ParameterEntry>>> cl(M("Teuchos"), "Ptr_Teuchos_ParameterEntry_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::Ptr<Teuchos::ParameterEntry>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("null_in") );

		cl.def( pybind11::init<class Teuchos::ParameterEntry *>(), pybind11::arg("ptr_in") );

		cl.def( pybind11::init( [](Teuchos::Ptr<Teuchos::ParameterEntry> const &o){ return new Teuchos::Ptr<Teuchos::ParameterEntry>(o); } ) );
		cl.def("assign", (class Teuchos::Ptr<class Teuchos::ParameterEntry> & (Teuchos::Ptr<Teuchos::ParameterEntry>::*)(const class Teuchos::Ptr<class Teuchos::ParameterEntry> &)) &Teuchos::Ptr<Teuchos::ParameterEntry>::operator=, "C++: Teuchos::Ptr<Teuchos::ParameterEntry>::operator=(const class Teuchos::Ptr<class Teuchos::ParameterEntry> &) --> class Teuchos::Ptr<class Teuchos::ParameterEntry> &", pybind11::return_value_policy::automatic, pybind11::arg("ptr"));
		cl.def("arrow", (class Teuchos::ParameterEntry * (Teuchos::Ptr<Teuchos::ParameterEntry>::*)() const) &Teuchos::Ptr<Teuchos::ParameterEntry>::operator->, "C++: Teuchos::Ptr<Teuchos::ParameterEntry>::operator->() const --> class Teuchos::ParameterEntry *", pybind11::return_value_policy::automatic);
		cl.def("dereference", (class Teuchos::ParameterEntry & (Teuchos::Ptr<Teuchos::ParameterEntry>::*)() const) &Teuchos::Ptr<Teuchos::ParameterEntry>::operator*, "C++: Teuchos::Ptr<Teuchos::ParameterEntry>::operator*() const --> class Teuchos::ParameterEntry &", pybind11::return_value_policy::automatic);
		cl.def("get", (class Teuchos::ParameterEntry * (Teuchos::Ptr<Teuchos::ParameterEntry>::*)() const) &Teuchos::Ptr<Teuchos::ParameterEntry>::get, "C++: Teuchos::Ptr<Teuchos::ParameterEntry>::get() const --> class Teuchos::ParameterEntry *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class Teuchos::ParameterEntry * (Teuchos::Ptr<Teuchos::ParameterEntry>::*)() const) &Teuchos::Ptr<Teuchos::ParameterEntry>::getRawPtr, "C++: Teuchos::Ptr<Teuchos::ParameterEntry>::getRawPtr() const --> class Teuchos::ParameterEntry *", pybind11::return_value_policy::automatic);
		cl.def("is_null", (bool (Teuchos::Ptr<Teuchos::ParameterEntry>::*)() const) &Teuchos::Ptr<Teuchos::ParameterEntry>::is_null, "C++: Teuchos::Ptr<Teuchos::ParameterEntry>::is_null() const --> bool");
		cl.def("assert_not_null", (const class Teuchos::Ptr<class Teuchos::ParameterEntry> & (Teuchos::Ptr<Teuchos::ParameterEntry>::*)() const) &Teuchos::Ptr<Teuchos::ParameterEntry>::assert_not_null, "C++: Teuchos::Ptr<Teuchos::ParameterEntry>::assert_not_null() const --> const class Teuchos::Ptr<class Teuchos::ParameterEntry> &", pybind11::return_value_policy::automatic);
		cl.def("ptr", (const class Teuchos::Ptr<class Teuchos::ParameterEntry> (Teuchos::Ptr<Teuchos::ParameterEntry>::*)() const) &Teuchos::Ptr<Teuchos::ParameterEntry>::ptr, "C++: Teuchos::Ptr<Teuchos::ParameterEntry>::ptr() const --> const class Teuchos::Ptr<class Teuchos::ParameterEntry>");
		cl.def("getConst", (class Teuchos::Ptr<const class Teuchos::ParameterEntry> (Teuchos::Ptr<Teuchos::ParameterEntry>::*)() const) &Teuchos::Ptr<Teuchos::ParameterEntry>::getConst, "C++: Teuchos::Ptr<Teuchos::ParameterEntry>::getConst() const --> class Teuchos::Ptr<const class Teuchos::ParameterEntry>");

		cl.def("__str__", [](Teuchos::Ptr<Teuchos::ParameterEntry> const &o) -> std::string { std::ostringstream s; Teuchos::operator<<(s, o); return s.str(); } );
	}
	{ // Teuchos::Ptr file:Teuchos_PtrDecl.hpp line:104
		pybind11::class_<Teuchos::Ptr<Teuchos::ParameterList>, Teuchos::RCP<Teuchos::Ptr<Teuchos::ParameterList>>> cl(M("Teuchos"), "Ptr_Teuchos_ParameterList_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::Ptr<Teuchos::ParameterList>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("null_in") );

		cl.def( pybind11::init<class Teuchos::ParameterList *>(), pybind11::arg("ptr_in") );

		cl.def( pybind11::init( [](Teuchos::Ptr<Teuchos::ParameterList> const &o){ return new Teuchos::Ptr<Teuchos::ParameterList>(o); } ) );
		cl.def("assign", (class Teuchos::Ptr<class Teuchos::ParameterList> & (Teuchos::Ptr<Teuchos::ParameterList>::*)(const class Teuchos::Ptr<class Teuchos::ParameterList> &)) &Teuchos::Ptr<Teuchos::ParameterList>::operator=, "C++: Teuchos::Ptr<Teuchos::ParameterList>::operator=(const class Teuchos::Ptr<class Teuchos::ParameterList> &) --> class Teuchos::Ptr<class Teuchos::ParameterList> &", pybind11::return_value_policy::automatic, pybind11::arg("ptr"));
		cl.def("arrow", (class Teuchos::ParameterList * (Teuchos::Ptr<Teuchos::ParameterList>::*)() const) &Teuchos::Ptr<Teuchos::ParameterList>::operator->, "C++: Teuchos::Ptr<Teuchos::ParameterList>::operator->() const --> class Teuchos::ParameterList *", pybind11::return_value_policy::automatic);
		cl.def("dereference", (class Teuchos::ParameterList & (Teuchos::Ptr<Teuchos::ParameterList>::*)() const) &Teuchos::Ptr<Teuchos::ParameterList>::operator*, "C++: Teuchos::Ptr<Teuchos::ParameterList>::operator*() const --> class Teuchos::ParameterList &", pybind11::return_value_policy::automatic);
		cl.def("get", (class Teuchos::ParameterList * (Teuchos::Ptr<Teuchos::ParameterList>::*)() const) &Teuchos::Ptr<Teuchos::ParameterList>::get, "C++: Teuchos::Ptr<Teuchos::ParameterList>::get() const --> class Teuchos::ParameterList *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class Teuchos::ParameterList * (Teuchos::Ptr<Teuchos::ParameterList>::*)() const) &Teuchos::Ptr<Teuchos::ParameterList>::getRawPtr, "C++: Teuchos::Ptr<Teuchos::ParameterList>::getRawPtr() const --> class Teuchos::ParameterList *", pybind11::return_value_policy::automatic);
		cl.def("is_null", (bool (Teuchos::Ptr<Teuchos::ParameterList>::*)() const) &Teuchos::Ptr<Teuchos::ParameterList>::is_null, "C++: Teuchos::Ptr<Teuchos::ParameterList>::is_null() const --> bool");
		cl.def("assert_not_null", (const class Teuchos::Ptr<class Teuchos::ParameterList> & (Teuchos::Ptr<Teuchos::ParameterList>::*)() const) &Teuchos::Ptr<Teuchos::ParameterList>::assert_not_null, "C++: Teuchos::Ptr<Teuchos::ParameterList>::assert_not_null() const --> const class Teuchos::Ptr<class Teuchos::ParameterList> &", pybind11::return_value_policy::automatic);
		cl.def("ptr", (const class Teuchos::Ptr<class Teuchos::ParameterList> (Teuchos::Ptr<Teuchos::ParameterList>::*)() const) &Teuchos::Ptr<Teuchos::ParameterList>::ptr, "C++: Teuchos::Ptr<Teuchos::ParameterList>::ptr() const --> const class Teuchos::Ptr<class Teuchos::ParameterList>");

		cl.def("__str__", [](Teuchos::Ptr<Teuchos::ParameterList> const &o) -> std::string { std::ostringstream s; Teuchos::operator<<(s, o); return s.str(); } );
	}
	// Teuchos::ERCPWeakNoDealloc file:Teuchos_RCPDecl.hpp line:75
	pybind11::enum_<Teuchos::ERCPWeakNoDealloc>(M("Teuchos"), "ERCPWeakNoDealloc", pybind11::arithmetic(), "", pybind11::module_local())
		.value("RCP_WEAK_NO_DEALLOC", Teuchos::RCP_WEAK_NO_DEALLOC)
		.export_values();

;

	// Teuchos::ERCPUndefinedWeakNoDealloc file:Teuchos_RCPDecl.hpp line:76
	pybind11::enum_<Teuchos::ERCPUndefinedWeakNoDealloc>(M("Teuchos"), "ERCPUndefinedWeakNoDealloc", pybind11::arithmetic(), "", pybind11::module_local())
		.value("RCP_UNDEFINED_WEAK_NO_DEALLOC", Teuchos::RCP_UNDEFINED_WEAK_NO_DEALLOC)
		.export_values();

;

	// Teuchos::ERCPUndefinedWithDealloc file:Teuchos_RCPDecl.hpp line:77
	pybind11::enum_<Teuchos::ERCPUndefinedWithDealloc>(M("Teuchos"), "ERCPUndefinedWithDealloc", pybind11::arithmetic(), "", pybind11::module_local())
		.value("RCP_UNDEFINED_WITH_DEALLOC", Teuchos::RCP_UNDEFINED_WITH_DEALLOC)
		.export_values();

;

	{ // Teuchos::RCPComp file:Teuchos_RCPDecl.hpp line:955
		pybind11::class_<Teuchos::RCPComp, Teuchos::RCP<Teuchos::RCPComp>> cl(M("Teuchos"), "RCPComp", "Struct for comparing two RCPs. Simply compares\n the raw pointers contained within the RCPs", pybind11::module_local());
		cl.def( pybind11::init( [](Teuchos::RCPComp const &o){ return new Teuchos::RCPComp(o); } ) );
		cl.def( pybind11::init( [](){ return new Teuchos::RCPComp(); } ) );
	}
	{ // Teuchos::RCPConstComp file:Teuchos_RCPDecl.hpp line:965
		pybind11::class_<Teuchos::RCPConstComp, Teuchos::RCP<Teuchos::RCPConstComp>> cl(M("Teuchos"), "RCPConstComp", "Struct for comparing two RCPs. Simply compares\n the raw pointers contained within the RCPs", pybind11::module_local());
		cl.def( pybind11::init( [](Teuchos::RCPConstComp const &o){ return new Teuchos::RCPConstComp(o); } ) );
		cl.def( pybind11::init( [](){ return new Teuchos::RCPConstComp(); } ) );
		cl.def("__call__", (bool (Teuchos::RCPConstComp::*)(const class Teuchos::RCP<const class Teuchos::ParameterEntry>, const class Teuchos::RCP<const class Teuchos::ParameterEntry>) const) &Teuchos::RCPConstComp::operator()<Teuchos::ParameterEntry,Teuchos::ParameterEntry>, "C++: Teuchos::RCPConstComp::operator()(const class Teuchos::RCP<const class Teuchos::ParameterEntry>, const class Teuchos::RCP<const class Teuchos::ParameterEntry>) const --> bool", pybind11::arg("p1"), pybind11::arg("p2"));
	}
	{ // Teuchos::DeallocNull file:Teuchos_RCPDecl.hpp line:996
		pybind11::class_<Teuchos::DeallocNull<Teuchos::ParameterEntry>, Teuchos::RCP<Teuchos::DeallocNull<Teuchos::ParameterEntry>>> cl(M("Teuchos"), "DeallocNull_Teuchos_ParameterEntry_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::DeallocNull<Teuchos::ParameterEntry>(); } ) );
		cl.def( pybind11::init( [](Teuchos::DeallocNull<Teuchos::ParameterEntry> const &o){ return new Teuchos::DeallocNull<Teuchos::ParameterEntry>(o); } ) );
		cl.def("free", (void (Teuchos::DeallocNull<Teuchos::ParameterEntry>::*)(class Teuchos::ParameterEntry *)) &Teuchos::DeallocNull<Teuchos::ParameterEntry>::free, "C++: Teuchos::DeallocNull<Teuchos::ParameterEntry>::free(class Teuchos::ParameterEntry *) --> void", pybind11::arg("ptr"));
	}
	{ // Teuchos::DeallocNull file:Teuchos_RCPDecl.hpp line:996
		pybind11::class_<Teuchos::DeallocNull<const Teuchos::ParameterEntry>, Teuchos::RCP<Teuchos::DeallocNull<const Teuchos::ParameterEntry>>> cl(M("Teuchos"), "DeallocNull_const_Teuchos_ParameterEntry_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::DeallocNull<const Teuchos::ParameterEntry>(); } ) );
		cl.def( pybind11::init( [](Teuchos::DeallocNull<const Teuchos::ParameterEntry> const &o){ return new Teuchos::DeallocNull<const Teuchos::ParameterEntry>(o); } ) );
		cl.def("free", (void (Teuchos::DeallocNull<const Teuchos::ParameterEntry>::*)(const class Teuchos::ParameterEntry *)) &Teuchos::DeallocNull<const Teuchos::ParameterEntry>::free, "C++: Teuchos::DeallocNull<const Teuchos::ParameterEntry>::free(const class Teuchos::ParameterEntry *) --> void", pybind11::arg("ptr"));
	}
	{ // Teuchos::DeallocNull file:Teuchos_RCPDecl.hpp line:996
		pybind11::class_<Teuchos::DeallocNull<const ROL::Vector<double>>, Teuchos::RCP<Teuchos::DeallocNull<const ROL::Vector<double>>>> cl(M("Teuchos"), "DeallocNull_const_ROL_Vector_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::DeallocNull<const ROL::Vector<double>>(); } ) );
		cl.def( pybind11::init( [](Teuchos::DeallocNull<const ROL::Vector<double>> const &o){ return new Teuchos::DeallocNull<const ROL::Vector<double>>(o); } ) );
		cl.def("free", (void (Teuchos::DeallocNull<const ROL::Vector<double>>::*)(const class ROL::Vector<double> *)) &Teuchos::DeallocNull<const ROL::Vector<double> >::free, "C++: Teuchos::DeallocNull<const ROL::Vector<double> >::free(const class ROL::Vector<double> *) --> void", pybind11::arg("ptr"));
	}
	{ // Teuchos::DeallocNull file:Teuchos_RCPDecl.hpp line:996
		pybind11::class_<Teuchos::DeallocNull<ROL::Objective<double>>, Teuchos::RCP<Teuchos::DeallocNull<ROL::Objective<double>>>> cl(M("Teuchos"), "DeallocNull_ROL_Objective_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::DeallocNull<ROL::Objective<double>>(); } ) );
		cl.def( pybind11::init( [](Teuchos::DeallocNull<ROL::Objective<double>> const &o){ return new Teuchos::DeallocNull<ROL::Objective<double>>(o); } ) );
		cl.def("free", (void (Teuchos::DeallocNull<ROL::Objective<double>>::*)(class ROL::Objective<double> *)) &Teuchos::DeallocNull<ROL::Objective<double> >::free, "C++: Teuchos::DeallocNull<ROL::Objective<double> >::free(class ROL::Objective<double> *) --> void", pybind11::arg("ptr"));
	}
	{ // Teuchos::DeallocNull file:Teuchos_RCPDecl.hpp line:996
		pybind11::class_<Teuchos::DeallocNull<ROL::Constraint<double>>, Teuchos::RCP<Teuchos::DeallocNull<ROL::Constraint<double>>>> cl(M("Teuchos"), "DeallocNull_ROL_Constraint_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::DeallocNull<ROL::Constraint<double>>(); } ) );
		cl.def( pybind11::init( [](Teuchos::DeallocNull<ROL::Constraint<double>> const &o){ return new Teuchos::DeallocNull<ROL::Constraint<double>>(o); } ) );
		cl.def("free", (void (Teuchos::DeallocNull<ROL::Constraint<double>>::*)(class ROL::Constraint<double> *)) &Teuchos::DeallocNull<ROL::Constraint<double> >::free, "C++: Teuchos::DeallocNull<ROL::Constraint<double> >::free(class ROL::Constraint<double> *) --> void", pybind11::arg("ptr"));
	}
	{ // Teuchos::DeallocNull file:Teuchos_RCPDecl.hpp line:996
		pybind11::class_<Teuchos::DeallocNull<ROL::BoundConstraint<double>>, Teuchos::RCP<Teuchos::DeallocNull<ROL::BoundConstraint<double>>>> cl(M("Teuchos"), "DeallocNull_ROL_BoundConstraint_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::DeallocNull<ROL::BoundConstraint<double>>(); } ) );
		cl.def( pybind11::init( [](Teuchos::DeallocNull<ROL::BoundConstraint<double>> const &o){ return new Teuchos::DeallocNull<ROL::BoundConstraint<double>>(o); } ) );
		cl.def("free", (void (Teuchos::DeallocNull<ROL::BoundConstraint<double>>::*)(class ROL::BoundConstraint<double> *)) &Teuchos::DeallocNull<ROL::BoundConstraint<double> >::free, "C++: Teuchos::DeallocNull<ROL::BoundConstraint<double> >::free(class ROL::BoundConstraint<double> *) --> void", pybind11::arg("ptr"));
	}
	{ // Teuchos::DeallocNull file:Teuchos_RCPDecl.hpp line:996
		pybind11::class_<Teuchos::DeallocNull<ROL::Vector<double>>, Teuchos::RCP<Teuchos::DeallocNull<ROL::Vector<double>>>> cl(M("Teuchos"), "DeallocNull_ROL_Vector_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::DeallocNull<ROL::Vector<double>>(); } ) );
		cl.def( pybind11::init( [](Teuchos::DeallocNull<ROL::Vector<double>> const &o){ return new Teuchos::DeallocNull<ROL::Vector<double>>(o); } ) );
		cl.def("free", (void (Teuchos::DeallocNull<ROL::Vector<double>>::*)(class ROL::Vector<double> *)) &Teuchos::DeallocNull<ROL::Vector<double> >::free, "C++: Teuchos::DeallocNull<ROL::Vector<double> >::free(class ROL::Vector<double> *) --> void", pybind11::arg("ptr"));
	}
	{ // Teuchos::DeallocNull file:Teuchos_RCPDecl.hpp line:996
		pybind11::class_<Teuchos::DeallocNull<ROL::ElasticObjective<double>>, Teuchos::RCP<Teuchos::DeallocNull<ROL::ElasticObjective<double>>>> cl(M("Teuchos"), "DeallocNull_ROL_ElasticObjective_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::DeallocNull<ROL::ElasticObjective<double>>(); } ) );
		cl.def( pybind11::init( [](Teuchos::DeallocNull<ROL::ElasticObjective<double>> const &o){ return new Teuchos::DeallocNull<ROL::ElasticObjective<double>>(o); } ) );
		cl.def("free", (void (Teuchos::DeallocNull<ROL::ElasticObjective<double>>::*)(class ROL::ElasticObjective<double> *)) &Teuchos::DeallocNull<ROL::ElasticObjective<double> >::free, "C++: Teuchos::DeallocNull<ROL::ElasticObjective<double> >::free(class ROL::ElasticObjective<double> *) --> void", pybind11::arg("ptr"));
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

	// Teuchos::rcp(struct ROL::AlgorithmState<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](struct ROL::AlgorithmState<double> * a0) -> Teuchos::RCP<struct ROL::AlgorithmState<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<struct ROL::AlgorithmState<double> > (*)(struct ROL::AlgorithmState<double> *, bool)) &Teuchos::rcp<ROL::AlgorithmState<double>>, "C++: Teuchos::rcp(struct ROL::AlgorithmState<double> *, bool) --> class Teuchos::RCP<struct ROL::AlgorithmState<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::CauchyPoint<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::CauchyPoint<double> * a0) -> Teuchos::RCP<class ROL::CauchyPoint<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::CauchyPoint<double> > (*)(class ROL::CauchyPoint<double> *, bool)) &Teuchos::rcp<ROL::CauchyPoint<double>>, "C++: Teuchos::rcp(class ROL::CauchyPoint<double> *, bool) --> class Teuchos::RCP<class ROL::CauchyPoint<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::DogLeg<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::DogLeg<double> * a0) -> Teuchos::RCP<class ROL::DogLeg<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::DogLeg<double> > (*)(class ROL::DogLeg<double> *, bool)) &Teuchos::rcp<ROL::DogLeg<double>>, "C++: Teuchos::rcp(class ROL::DogLeg<double> *, bool) --> class Teuchos::RCP<class ROL::DogLeg<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::DoubleDogLeg<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::DoubleDogLeg<double> * a0) -> Teuchos::RCP<class ROL::DoubleDogLeg<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::DoubleDogLeg<double> > (*)(class ROL::DoubleDogLeg<double> *, bool)) &Teuchos::rcp<ROL::DoubleDogLeg<double>>, "C++: Teuchos::rcp(class ROL::DoubleDogLeg<double> *, bool) --> class Teuchos::RCP<class ROL::DoubleDogLeg<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::TruncatedCG<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::TruncatedCG<double> * a0) -> Teuchos::RCP<class ROL::TruncatedCG<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::TruncatedCG<double> > (*)(class ROL::TruncatedCG<double> *, bool)) &Teuchos::rcp<ROL::TruncatedCG<double>>, "C++: Teuchos::rcp(class ROL::TruncatedCG<double> *, bool) --> class Teuchos::RCP<class ROL::TruncatedCG<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::LinMore<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::LinMore<double> * a0) -> Teuchos::RCP<class ROL::LinMore<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::LinMore<double> > (*)(class ROL::LinMore<double> *, bool)) &Teuchos::rcp<ROL::LinMore<double>>, "C++: Teuchos::rcp(class ROL::LinMore<double> *, bool) --> class Teuchos::RCP<class ROL::LinMore<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(struct ROL::StepState<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](struct ROL::StepState<double> * a0) -> Teuchos::RCP<struct ROL::StepState<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<struct ROL::StepState<double> > (*)(struct ROL::StepState<double> *, bool)) &Teuchos::rcp<ROL::StepState<double>>, "C++: Teuchos::rcp(struct ROL::StepState<double> *, bool) --> class Teuchos::RCP<struct ROL::StepState<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

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

	// Teuchos::rcp(class ROL::KelleySachsModel<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::KelleySachsModel<double> * a0) -> Teuchos::RCP<class ROL::KelleySachsModel<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::KelleySachsModel<double> > (*)(class ROL::KelleySachsModel<double> *, bool)) &Teuchos::rcp<ROL::KelleySachsModel<double>>, "C++: Teuchos::rcp(class ROL::KelleySachsModel<double> *, bool) --> class Teuchos::RCP<class ROL::KelleySachsModel<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::ColemanLiModel<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::ColemanLiModel<double> * a0) -> Teuchos::RCP<class ROL::ColemanLiModel<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::ColemanLiModel<double> > (*)(class ROL::ColemanLiModel<double> *, bool)) &Teuchos::rcp<ROL::ColemanLiModel<double>>, "C++: Teuchos::rcp(class ROL::ColemanLiModel<double> *, bool) --> class Teuchos::RCP<class ROL::ColemanLiModel<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::LinMoreModel<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::LinMoreModel<double> * a0) -> Teuchos::RCP<class ROL::LinMoreModel<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::LinMoreModel<double> > (*)(class ROL::LinMoreModel<double> *, bool)) &Teuchos::rcp<ROL::LinMoreModel<double>>, "C++: Teuchos::rcp(class ROL::LinMoreModel<double> *, bool) --> class Teuchos::RCP<class ROL::LinMoreModel<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class ROL::TrustRegionModel<double> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class ROL::TrustRegionModel<double> * a0) -> Teuchos::RCP<class ROL::TrustRegionModel<double> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class ROL::TrustRegionModel<double> > (*)(class ROL::TrustRegionModel<double> *, bool)) &Teuchos::rcp<ROL::TrustRegionModel<double>>, "C++: Teuchos::rcp(class ROL::TrustRegionModel<double> *, bool) --> class Teuchos::RCP<class ROL::TrustRegionModel<double> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

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

	// Teuchos::is_null(const class Teuchos::RCP<class ROL::TrustRegionModel<double> > &) file:Teuchos_RCP.hpp line:715
	M("Teuchos").def("is_null", (bool (*)(const class Teuchos::RCP<class ROL::TrustRegionModel<double> > &)) &Teuchos::is_null<ROL::TrustRegionModel<double>>, "C++: Teuchos::is_null(const class Teuchos::RCP<class ROL::TrustRegionModel<double> > &) --> bool", pybind11::arg("p"));

	// Teuchos::nonnull(const class Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > > &) file:Teuchos_RCP.hpp line:723
	M("Teuchos").def("nonnull", (bool (*)(const class Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > > &)) &Teuchos::nonnull<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>, "C++: Teuchos::nonnull(const class Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > > &) --> bool", pybind11::arg("p"));

	// Teuchos::nonnull(const class Teuchos::RCP<const class Teuchos::ParameterEntryValidator> &) file:Teuchos_RCP.hpp line:723
	M("Teuchos").def("nonnull", (bool (*)(const class Teuchos::RCP<const class Teuchos::ParameterEntryValidator> &)) &Teuchos::nonnull<const Teuchos::ParameterEntryValidator>, "C++: Teuchos::nonnull(const class Teuchos::RCP<const class Teuchos::ParameterEntryValidator> &) --> bool", pybind11::arg("p"));

	// Teuchos::rcp_static_cast(const class Teuchos::RCP<const struct ROL::TypeB::AlgorithmState<double> > &) file:Teuchos_RCP.hpp line:775
	M("Teuchos").def("rcp_static_cast", (class Teuchos::RCP<const struct ROL::TypeB::AlgorithmState<double> > (*)(const class Teuchos::RCP<const struct ROL::TypeB::AlgorithmState<double> > &)) &Teuchos::rcp_static_cast<const ROL::TypeB::AlgorithmState<double>,const ROL::TypeB::AlgorithmState<double>>, "C++: Teuchos::rcp_static_cast(const class Teuchos::RCP<const struct ROL::TypeB::AlgorithmState<double> > &) --> class Teuchos::RCP<const struct ROL::TypeB::AlgorithmState<double> >", pybind11::arg("p1"));

	// Teuchos::rcp_const_cast(const class Teuchos::RCP<const class ROL::Vector<double> > &) file:Teuchos_RCP.hpp line:786
	M("Teuchos").def("rcp_const_cast", (class Teuchos::RCP<class ROL::Vector<double> > (*)(const class Teuchos::RCP<const class ROL::Vector<double> > &)) &Teuchos::rcp_const_cast<ROL::Vector<double>,const ROL::Vector<double>>, "C++: Teuchos::rcp_const_cast(const class Teuchos::RCP<const class ROL::Vector<double> > &) --> class Teuchos::RCP<class ROL::Vector<double> >", pybind11::arg("p1"));

	{ // Teuchos::m_bad_cast file:Teuchos_dyn_cast.hpp line:60
		pybind11::class_<Teuchos::m_bad_cast, Teuchos::RCP<Teuchos::m_bad_cast>, PyCallBack_Teuchos_m_bad_cast, std::bad_cast> cl(M("Teuchos"), "m_bad_cast", "Exception class for bad cast.\n\n\nWe create this class so that we may throw a bad_cast when appropriate and\nstill use the TEUCHOS_TEST_FOR_EXCEPTION macro.  We recommend users try to catch a\nbad_cast.", pybind11::module_local());
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("what_arg") );

		cl.def( pybind11::init( [](PyCallBack_Teuchos_m_bad_cast const &o){ return new PyCallBack_Teuchos_m_bad_cast(o); } ) );
		cl.def( pybind11::init( [](Teuchos::m_bad_cast const &o){ return new Teuchos::m_bad_cast(o); } ) );
		cl.def("what", (const char * (Teuchos::m_bad_cast::*)() const) &Teuchos::m_bad_cast::what, "C++: Teuchos::m_bad_cast::what() const --> const char *", pybind11::return_value_policy::automatic);
		cl.def("assign", (class Teuchos::m_bad_cast & (Teuchos::m_bad_cast::*)(const class Teuchos::m_bad_cast &)) &Teuchos::m_bad_cast::operator=, "C++: Teuchos::m_bad_cast::operator=(const class Teuchos::m_bad_cast &) --> class Teuchos::m_bad_cast &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	// Teuchos::dyn_cast_throw_exception(const std::string &, const std::string &, const std::string &) file:Teuchos_dyn_cast.hpp line:70
	M("Teuchos").def("dyn_cast_throw_exception", (void (*)(const std::string &, const std::string &, const std::string &)) &Teuchos::dyn_cast_throw_exception, "C++: Teuchos::dyn_cast_throw_exception(const std::string &, const std::string &, const std::string &) --> void", pybind11::arg("T_from"), pybind11::arg("T_from_concr"), pybind11::arg("T_to"));

	// Teuchos::dyn_cast(class ROL::TrustRegionModel<double> &) file:Teuchos_dyn_cast.hpp line:173
	M("Teuchos").def("dyn_cast", (class ROL::KelleySachsModel<double> & (*)(class ROL::TrustRegionModel<double> &)) &Teuchos::dyn_cast<ROL::KelleySachsModel<double>,ROL::TrustRegionModel<double>>, "C++: Teuchos::dyn_cast(class ROL::TrustRegionModel<double> &) --> class ROL::KelleySachsModel<double> &", pybind11::return_value_policy::automatic, pybind11::arg("from"));

	// Teuchos::dyn_cast(class ROL::TrustRegionModel<double> &) file:Teuchos_dyn_cast.hpp line:173
	M("Teuchos").def("dyn_cast", (class ROL::ColemanLiModel<double> & (*)(class ROL::TrustRegionModel<double> &)) &Teuchos::dyn_cast<ROL::ColemanLiModel<double>,ROL::TrustRegionModel<double>>, "C++: Teuchos::dyn_cast(class ROL::TrustRegionModel<double> &) --> class ROL::ColemanLiModel<double> &", pybind11::return_value_policy::automatic, pybind11::arg("from"));

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
