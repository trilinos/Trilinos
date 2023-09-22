#include <ROL_BoundConstraint.hpp>
#include <ROL_Constraint.hpp>
#include <ROL_Elementwise_Function.hpp>
#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_Objective.hpp>
#include <ROL_PolyhedralProjection.hpp>
#include <ROL_Problem.hpp>
#include <ROL_StatusTest.hpp>
#include <ROL_TypeG_Algorithm.hpp>
#include <ROL_TypeG_AlgorithmFactory.hpp>
#include <ROL_TypeG_AugmentedLagrangianAlgorithm.hpp>
#include <ROL_TypeG_InteriorPointAlgorithm.hpp>
#include <ROL_TypeG_MoreauYosidaAlgorithm.hpp>
#include <ROL_TypeG_StabilizedLCLAlgorithm.hpp>
#include <ROL_Types.hpp>
#include <ROL_UpdateType.hpp>
#include <ROL_Vector.hpp>
#include <Teuchos_ENull.hpp>
#include <Teuchos_FilteredIterator.hpp>
#include <Teuchos_ParameterEntry.hpp>
#include <Teuchos_ParameterEntryValidator.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ParameterListModifier.hpp>
#include <Teuchos_PtrDecl.hpp>
#include <Teuchos_RCPDecl.hpp>
#include <Teuchos_RCPNode.hpp>
#include <Teuchos_StringIndexedOrderedValueObjectContainer.hpp>
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
#include <vector>

#include <functional>
#include "PyROL_Smart_Holder.hpp"
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

// ROL::TypeG::Algorithm file:ROL_TypeG_Algorithm.hpp line:90
struct PyCallBack_ROL_TypeG_Algorithm_double_t : public ROL::TypeG::Algorithm<double> {
	using ROL::TypeG::Algorithm<double>::Algorithm;

	void run(class ROL::Problem<double> & a0, std::ostream & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::Algorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, std::ostream & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::Algorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4, std::ostream & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::Algorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, std::ostream & a6) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::Algorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, std::ostream & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::Algorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7, std::ostream & a8) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::Algorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, std::ostream & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::Algorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"Algorithm::run\"");
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6, std::ostream & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::Algorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7, std::ostream & a8) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::Algorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9, std::ostream & a10) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::Algorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10, std::ostream & a11) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::Algorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, std::ostream & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::Algorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, std::ostream & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::Algorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, std::ostream & a8) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::Algorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, std::ostream & a9) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::Algorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9, std::ostream & a10) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::Algorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9, std::ostream & a10) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::Algorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9, std::ostream & a10) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::Algorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9, const class ROL::Vector<double> & a10, std::ostream & a11) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::Algorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9, class ROL::Constraint<double> & a10, class ROL::Vector<double> & a11, const class ROL::Vector<double> & a12, std::ostream & a13) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::Algorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10, class ROL::Constraint<double> & a11, class ROL::Vector<double> & a12, const class ROL::Vector<double> & a13, std::ostream & a14) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::Algorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14);
	}
	void writeHeader(std::ostream & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::Algorithm<double> *>(this), "writeHeader");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::writeHeader(a0);
	}
	void writeName(std::ostream & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::Algorithm<double> *>(this), "writeName");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::writeName(a0);
	}
	void writeOutput(std::ostream & a0, const bool a1) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::Algorithm<double> *>(this), "writeOutput");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::writeOutput(a0, a1);
	}
	void writeExitStatus(std::ostream & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::Algorithm<double> *>(this), "writeExitStatus");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::writeExitStatus(a0);
	}
};

// ROL::TypeG::AugmentedLagrangianAlgorithm file:ROL_TypeG_AugmentedLagrangianAlgorithm.hpp line:60
struct PyCallBack_ROL_TypeG_AugmentedLagrangianAlgorithm_double_t : public ROL::TypeG::AugmentedLagrangianAlgorithm<double> {
	using ROL::TypeG::AugmentedLagrangianAlgorithm<double>::AugmentedLagrangianAlgorithm;

	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, std::ostream & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::AugmentedLagrangianAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return AugmentedLagrangianAlgorithm::run(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	void writeHeader(std::ostream & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::AugmentedLagrangianAlgorithm<double> *>(this), "writeHeader");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return AugmentedLagrangianAlgorithm::writeHeader(a0);
	}
	void writeName(std::ostream & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::AugmentedLagrangianAlgorithm<double> *>(this), "writeName");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return AugmentedLagrangianAlgorithm::writeName(a0);
	}
	void writeOutput(std::ostream & a0, const bool a1) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::AugmentedLagrangianAlgorithm<double> *>(this), "writeOutput");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return AugmentedLagrangianAlgorithm::writeOutput(a0, a1);
	}
	void run(class ROL::Problem<double> & a0, std::ostream & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::AugmentedLagrangianAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, std::ostream & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::AugmentedLagrangianAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4, std::ostream & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::AugmentedLagrangianAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, std::ostream & a6) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::AugmentedLagrangianAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, std::ostream & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::AugmentedLagrangianAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7, std::ostream & a8) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::AugmentedLagrangianAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6, std::ostream & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::AugmentedLagrangianAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7, std::ostream & a8) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::AugmentedLagrangianAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9, std::ostream & a10) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::AugmentedLagrangianAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10, std::ostream & a11) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::AugmentedLagrangianAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, std::ostream & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::AugmentedLagrangianAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, std::ostream & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::AugmentedLagrangianAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, std::ostream & a8) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::AugmentedLagrangianAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, std::ostream & a9) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::AugmentedLagrangianAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9, std::ostream & a10) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::AugmentedLagrangianAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9, std::ostream & a10) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::AugmentedLagrangianAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9, std::ostream & a10) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::AugmentedLagrangianAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9, const class ROL::Vector<double> & a10, std::ostream & a11) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::AugmentedLagrangianAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9, class ROL::Constraint<double> & a10, class ROL::Vector<double> & a11, const class ROL::Vector<double> & a12, std::ostream & a13) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::AugmentedLagrangianAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10, class ROL::Constraint<double> & a11, class ROL::Vector<double> & a12, const class ROL::Vector<double> & a13, std::ostream & a14) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::AugmentedLagrangianAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14);
	}
	void writeExitStatus(std::ostream & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::AugmentedLagrangianAlgorithm<double> *>(this), "writeExitStatus");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::writeExitStatus(a0);
	}
};

// ROL::TypeG::MoreauYosidaAlgorithm file:ROL_TypeG_MoreauYosidaAlgorithm.hpp line:58
struct PyCallBack_ROL_TypeG_MoreauYosidaAlgorithm_double_t : public ROL::TypeG::MoreauYosidaAlgorithm<double> {
	using ROL::TypeG::MoreauYosidaAlgorithm<double>::MoreauYosidaAlgorithm;

	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, std::ostream & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::MoreauYosidaAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return MoreauYosidaAlgorithm::run(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	void writeHeader(std::ostream & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::MoreauYosidaAlgorithm<double> *>(this), "writeHeader");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return MoreauYosidaAlgorithm::writeHeader(a0);
	}
	void writeName(std::ostream & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::MoreauYosidaAlgorithm<double> *>(this), "writeName");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return MoreauYosidaAlgorithm::writeName(a0);
	}
	void writeOutput(std::ostream & a0, const bool a1) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::MoreauYosidaAlgorithm<double> *>(this), "writeOutput");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return MoreauYosidaAlgorithm::writeOutput(a0, a1);
	}
	void run(class ROL::Problem<double> & a0, std::ostream & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::MoreauYosidaAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, std::ostream & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::MoreauYosidaAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4, std::ostream & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::MoreauYosidaAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, std::ostream & a6) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::MoreauYosidaAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, std::ostream & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::MoreauYosidaAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7, std::ostream & a8) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::MoreauYosidaAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6, std::ostream & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::MoreauYosidaAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7, std::ostream & a8) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::MoreauYosidaAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9, std::ostream & a10) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::MoreauYosidaAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10, std::ostream & a11) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::MoreauYosidaAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, std::ostream & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::MoreauYosidaAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, std::ostream & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::MoreauYosidaAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, std::ostream & a8) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::MoreauYosidaAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, std::ostream & a9) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::MoreauYosidaAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9, std::ostream & a10) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::MoreauYosidaAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9, std::ostream & a10) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::MoreauYosidaAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9, std::ostream & a10) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::MoreauYosidaAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9, const class ROL::Vector<double> & a10, std::ostream & a11) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::MoreauYosidaAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9, class ROL::Constraint<double> & a10, class ROL::Vector<double> & a11, const class ROL::Vector<double> & a12, std::ostream & a13) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::MoreauYosidaAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10, class ROL::Constraint<double> & a11, class ROL::Vector<double> & a12, const class ROL::Vector<double> & a13, std::ostream & a14) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::MoreauYosidaAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14);
	}
	void writeExitStatus(std::ostream & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::MoreauYosidaAlgorithm<double> *>(this), "writeExitStatus");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::writeExitStatus(a0);
	}
};

// ROL::TypeG::InteriorPointAlgorithm file:ROL_TypeG_InteriorPointAlgorithm.hpp line:58
struct PyCallBack_ROL_TypeG_InteriorPointAlgorithm_double_t : public ROL::TypeG::InteriorPointAlgorithm<double> {
	using ROL::TypeG::InteriorPointAlgorithm<double>::InteriorPointAlgorithm;

	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, std::ostream & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::InteriorPointAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return InteriorPointAlgorithm::run(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	void writeHeader(std::ostream & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::InteriorPointAlgorithm<double> *>(this), "writeHeader");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return InteriorPointAlgorithm::writeHeader(a0);
	}
	void writeName(std::ostream & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::InteriorPointAlgorithm<double> *>(this), "writeName");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return InteriorPointAlgorithm::writeName(a0);
	}
	void writeOutput(std::ostream & a0, const bool a1) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::InteriorPointAlgorithm<double> *>(this), "writeOutput");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return InteriorPointAlgorithm::writeOutput(a0, a1);
	}
	void run(class ROL::Problem<double> & a0, std::ostream & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::InteriorPointAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, std::ostream & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::InteriorPointAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4, std::ostream & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::InteriorPointAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, std::ostream & a6) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::InteriorPointAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, std::ostream & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::InteriorPointAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7, std::ostream & a8) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::InteriorPointAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6, std::ostream & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::InteriorPointAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7, std::ostream & a8) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::InteriorPointAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9, std::ostream & a10) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::InteriorPointAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10, std::ostream & a11) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::InteriorPointAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, std::ostream & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::InteriorPointAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, std::ostream & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::InteriorPointAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, std::ostream & a8) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::InteriorPointAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, std::ostream & a9) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::InteriorPointAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9, std::ostream & a10) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::InteriorPointAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9, std::ostream & a10) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::InteriorPointAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9, std::ostream & a10) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::InteriorPointAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9, const class ROL::Vector<double> & a10, std::ostream & a11) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::InteriorPointAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9, class ROL::Constraint<double> & a10, class ROL::Vector<double> & a11, const class ROL::Vector<double> & a12, std::ostream & a13) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::InteriorPointAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10, class ROL::Constraint<double> & a11, class ROL::Vector<double> & a12, const class ROL::Vector<double> & a13, std::ostream & a14) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::InteriorPointAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14);
	}
	void writeExitStatus(std::ostream & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::InteriorPointAlgorithm<double> *>(this), "writeExitStatus");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::writeExitStatus(a0);
	}
};

// ROL::TypeG::StabilizedLCLAlgorithm file:ROL_TypeG_StabilizedLCLAlgorithm.hpp line:61
struct PyCallBack_ROL_TypeG_StabilizedLCLAlgorithm_double_t : public ROL::TypeG::StabilizedLCLAlgorithm<double> {
	using ROL::TypeG::StabilizedLCLAlgorithm<double>::StabilizedLCLAlgorithm;

	void run(class ROL::Problem<double> & a0, std::ostream & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::StabilizedLCLAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return StabilizedLCLAlgorithm::run(a0, a1);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, std::ostream & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::StabilizedLCLAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return StabilizedLCLAlgorithm::run(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	void writeHeader(std::ostream & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::StabilizedLCLAlgorithm<double> *>(this), "writeHeader");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return StabilizedLCLAlgorithm::writeHeader(a0);
	}
	void writeName(std::ostream & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::StabilizedLCLAlgorithm<double> *>(this), "writeName");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return StabilizedLCLAlgorithm::writeName(a0);
	}
	void writeOutput(std::ostream & a0, const bool a1) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::StabilizedLCLAlgorithm<double> *>(this), "writeOutput");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return StabilizedLCLAlgorithm::writeOutput(a0, a1);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, std::ostream & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::StabilizedLCLAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4, std::ostream & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::StabilizedLCLAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, std::ostream & a6) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::StabilizedLCLAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, std::ostream & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::StabilizedLCLAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7, std::ostream & a8) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::StabilizedLCLAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6, std::ostream & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::StabilizedLCLAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7, std::ostream & a8) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::StabilizedLCLAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9, std::ostream & a10) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::StabilizedLCLAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10, std::ostream & a11) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::StabilizedLCLAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, std::ostream & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::StabilizedLCLAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, std::ostream & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::StabilizedLCLAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, std::ostream & a8) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::StabilizedLCLAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, std::ostream & a9) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::StabilizedLCLAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9, std::ostream & a10) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::StabilizedLCLAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9, std::ostream & a10) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::StabilizedLCLAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9, std::ostream & a10) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::StabilizedLCLAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9, const class ROL::Vector<double> & a10, std::ostream & a11) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::StabilizedLCLAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9, class ROL::Constraint<double> & a10, class ROL::Vector<double> & a11, const class ROL::Vector<double> & a12, std::ostream & a13) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::StabilizedLCLAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10, class ROL::Constraint<double> & a11, class ROL::Vector<double> & a12, const class ROL::Vector<double> & a13, std::ostream & a14) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::StabilizedLCLAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14);
	}
	void writeExitStatus(std::ostream & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeG::StabilizedLCLAlgorithm<double> *>(this), "writeExitStatus");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::writeExitStatus(a0);
	}
};

void bind_ROL_TypeG_Algorithm(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // ROL::TypeG::AlgorithmState file:ROL_TypeG_Algorithm.hpp line:62
		pybind11::class_<ROL::TypeG::AlgorithmState<double>, Teuchos::RCP<ROL::TypeG::AlgorithmState<double>>, ROL::AlgorithmState<double>> cl(M("ROL::TypeG"), "AlgorithmState_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new ROL::TypeG::AlgorithmState<double>(); } ) );
		cl.def( pybind11::init( [](ROL::TypeG::AlgorithmState<double> const &o){ return new ROL::TypeG::AlgorithmState<double>(o); } ) );
		cl.def_readwrite("searchSize", &ROL::TypeG::AlgorithmState<double>::searchSize);
		cl.def_readwrite("stepVec", &ROL::TypeG::AlgorithmState<double>::stepVec);
		cl.def_readwrite("gradientVec", &ROL::TypeG::AlgorithmState<double>::gradientVec);
		cl.def_readwrite("constraintVec", &ROL::TypeG::AlgorithmState<double>::constraintVec);
		cl.def("reset", (void (ROL::TypeG::AlgorithmState<double>::*)()) &ROL::TypeG::AlgorithmState<double>::reset, "C++: ROL::TypeG::AlgorithmState<double>::reset() --> void");
		cl.def("assign", (struct ROL::TypeG::AlgorithmState<double> & (ROL::TypeG::AlgorithmState<double>::*)(const struct ROL::TypeG::AlgorithmState<double> &)) &ROL::TypeG::AlgorithmState<double>::operator=, "C++: ROL::TypeG::AlgorithmState<double>::operator=(const struct ROL::TypeG::AlgorithmState<double> &) --> struct ROL::TypeG::AlgorithmState<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def_readwrite("iter", &ROL::AlgorithmState<double>::iter);
		cl.def_readwrite("minIter", &ROL::AlgorithmState<double>::minIter);
		cl.def_readwrite("nfval", &ROL::AlgorithmState<double>::nfval);
		cl.def_readwrite("ncval", &ROL::AlgorithmState<double>::ncval);
		cl.def_readwrite("ngrad", &ROL::AlgorithmState<double>::ngrad);
		cl.def_readwrite("value", &ROL::AlgorithmState<double>::value);
		cl.def_readwrite("minValue", &ROL::AlgorithmState<double>::minValue);
		cl.def_readwrite("gnorm", &ROL::AlgorithmState<double>::gnorm);
		cl.def_readwrite("cnorm", &ROL::AlgorithmState<double>::cnorm);
		cl.def_readwrite("snorm", &ROL::AlgorithmState<double>::snorm);
		cl.def_readwrite("aggregateGradientNorm", &ROL::AlgorithmState<double>::aggregateGradientNorm);
		cl.def_readwrite("aggregateModelError", &ROL::AlgorithmState<double>::aggregateModelError);
		cl.def_readwrite("flag", &ROL::AlgorithmState<double>::flag);
		cl.def_readwrite("iterateVec", &ROL::AlgorithmState<double>::iterateVec);
		cl.def_readwrite("lagmultVec", &ROL::AlgorithmState<double>::lagmultVec);
		cl.def_readwrite("minIterVec", &ROL::AlgorithmState<double>::minIterVec);
		cl.def_readwrite("statusFlag", &ROL::AlgorithmState<double>::statusFlag);
		cl.def("reset", (void (ROL::AlgorithmState<double>::*)()) &ROL::AlgorithmState<double>::reset, "C++: ROL::AlgorithmState<double>::reset() --> void");
		cl.def("assign", (struct ROL::AlgorithmState<double> & (ROL::AlgorithmState<double>::*)(const struct ROL::AlgorithmState<double> &)) &ROL::AlgorithmState<double>::operator=, "C++: ROL::AlgorithmState<double>::operator=(const struct ROL::AlgorithmState<double> &) --> struct ROL::AlgorithmState<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::TypeG::Algorithm file:ROL_TypeG_Algorithm.hpp line:90
		pybind11::class_<ROL::TypeG::Algorithm<double>, Teuchos::RCP<ROL::TypeG::Algorithm<double>>, PyCallBack_ROL_TypeG_Algorithm_double_t> cl(M("ROL::TypeG"), "Algorithm_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new PyCallBack_ROL_TypeG_Algorithm_double_t(); } ) );
		cl.def(pybind11::init<PyCallBack_ROL_TypeG_Algorithm_double_t const &>());
		cl.def("setStatusTest", [](ROL::TypeG::Algorithm<double> &o, const class Teuchos::RCP<class ROL::StatusTest<double> > & a0) -> void { return o.setStatusTest(a0); }, "", pybind11::arg("status"));
		cl.def("setStatusTest", (void (ROL::TypeG::Algorithm<double>::*)(const class Teuchos::RCP<class ROL::StatusTest<double> > &, bool)) &ROL::TypeG::Algorithm<double>::setStatusTest, "C++: ROL::TypeG::Algorithm<double>::setStatusTest(const class Teuchos::RCP<class ROL::StatusTest<double> > &, bool) --> void", pybind11::arg("status"), pybind11::arg("combineStatus"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Problem<double> & a0) -> void { return o.run(a0); }, "", pybind11::arg("problem"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Problem<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Problem<double> &, std::ostream &) --> void", pybind11::arg("problem"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4) -> void { return o.run(a0, a1, a2, a3, a4); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4) -> void { return o.run(a0, a1, a2, a3, a4); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5) -> void { return o.run(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9, const class ROL::Vector<double> & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9, class ROL::Constraint<double> & a10, class ROL::Vector<double> & a11, const class ROL::Vector<double> & a12) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10, class ROL::Constraint<double> & a11, class ROL::Vector<double> & a12, const class ROL::Vector<double> & a13) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("writeHeader", (void (ROL::TypeG::Algorithm<double>::*)(std::ostream &) const) &ROL::TypeG::Algorithm<double>::writeHeader, "C++: ROL::TypeG::Algorithm<double>::writeHeader(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("writeName", (void (ROL::TypeG::Algorithm<double>::*)(std::ostream &) const) &ROL::TypeG::Algorithm<double>::writeName, "C++: ROL::TypeG::Algorithm<double>::writeName(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("writeOutput", [](ROL::TypeG::Algorithm<double> const &o, std::ostream & a0) -> void { return o.writeOutput(a0); }, "", pybind11::arg("os"));
		cl.def("writeOutput", (void (ROL::TypeG::Algorithm<double>::*)(std::ostream &, const bool) const) &ROL::TypeG::Algorithm<double>::writeOutput, "C++: ROL::TypeG::Algorithm<double>::writeOutput(std::ostream &, const bool) const --> void", pybind11::arg("os"), pybind11::arg("write_header"));
		cl.def("writeExitStatus", (void (ROL::TypeG::Algorithm<double>::*)(std::ostream &) const) &ROL::TypeG::Algorithm<double>::writeExitStatus, "C++: ROL::TypeG::Algorithm<double>::writeExitStatus(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("getState", (class Teuchos::RCP<const struct ROL::TypeG::AlgorithmState<double> > (ROL::TypeG::Algorithm<double>::*)() const) &ROL::TypeG::Algorithm<double>::getState, "C++: ROL::TypeG::Algorithm<double>::getState() const --> class Teuchos::RCP<const struct ROL::TypeG::AlgorithmState<double> >");
		cl.def("reset", (void (ROL::TypeG::Algorithm<double>::*)()) &ROL::TypeG::Algorithm<double>::reset, "C++: ROL::TypeG::Algorithm<double>::reset() --> void");
	}
	{ // ROL::TypeG::AugmentedLagrangianAlgorithm file:ROL_TypeG_AugmentedLagrangianAlgorithm.hpp line:60
		pybind11::class_<ROL::TypeG::AugmentedLagrangianAlgorithm<double>, Teuchos::RCP<ROL::TypeG::AugmentedLagrangianAlgorithm<double>>, PyCallBack_ROL_TypeG_AugmentedLagrangianAlgorithm_double_t, ROL::TypeG::Algorithm<double>> cl(M("ROL::TypeG"), "AugmentedLagrangianAlgorithm_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("list") );

		cl.def( pybind11::init( [](PyCallBack_ROL_TypeG_AugmentedLagrangianAlgorithm_double_t const &o){ return new PyCallBack_ROL_TypeG_AugmentedLagrangianAlgorithm_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> const &o){ return new ROL::TypeG::AugmentedLagrangianAlgorithm<double>(o); } ) );
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10, class ROL::Constraint<double> & a11, class ROL::Vector<double> & a12, const class ROL::Vector<double> & a13) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10, class ROL::Constraint<double> & a11, class ROL::Vector<double> & a12, const class ROL::Vector<double> & a13, std::ostream & a14) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9, class ROL::Constraint<double> & a10, class ROL::Vector<double> & a11, const class ROL::Vector<double> & a12) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9, class ROL::Constraint<double> & a10, class ROL::Vector<double> & a11, const class ROL::Vector<double> & a12, std::ostream & a13) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9, const class ROL::Vector<double> & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9, const class ROL::Vector<double> & a10, std::ostream & a11) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9, std::ostream & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9, std::ostream & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9, std::ostream & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, std::ostream & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, std::ostream & a8) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, std::ostream & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, std::ostream & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10, std::ostream & a11) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9, std::ostream & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7, std::ostream & a8) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6, std::ostream & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7, std::ostream & a8) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, std::ostream & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5) -> void { return o.run(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, std::ostream & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4) -> void { return o.run(a0, a1, a2, a3, a4); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4, std::ostream & a5) -> void { return o.run(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4) -> void { return o.run(a0, a1, a2, a3, a4); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, std::ostream & a5) -> void { return o.run(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Problem<double> & a0) -> void { return o.run(a0); }, "", pybind11::arg("problem"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Problem<double> & a0, std::ostream & a1) -> void { return o.run(a0, a1); }, "", pybind11::arg("problem"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"));
		cl.def("run", (void (ROL::TypeG::AugmentedLagrangianAlgorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::AugmentedLagrangianAlgorithm<double>::run, "C++: ROL::TypeG::AugmentedLagrangianAlgorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("outStream"));
		cl.def("writeHeader", (void (ROL::TypeG::AugmentedLagrangianAlgorithm<double>::*)(std::ostream &) const) &ROL::TypeG::AugmentedLagrangianAlgorithm<double>::writeHeader, "C++: ROL::TypeG::AugmentedLagrangianAlgorithm<double>::writeHeader(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("writeName", (void (ROL::TypeG::AugmentedLagrangianAlgorithm<double>::*)(std::ostream &) const) &ROL::TypeG::AugmentedLagrangianAlgorithm<double>::writeName, "C++: ROL::TypeG::AugmentedLagrangianAlgorithm<double>::writeName(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("writeOutput", [](ROL::TypeG::AugmentedLagrangianAlgorithm<double> const &o, std::ostream & a0) -> void { return o.writeOutput(a0); }, "", pybind11::arg("os"));
		cl.def("writeOutput", (void (ROL::TypeG::AugmentedLagrangianAlgorithm<double>::*)(std::ostream &, const bool) const) &ROL::TypeG::AugmentedLagrangianAlgorithm<double>::writeOutput, "C++: ROL::TypeG::AugmentedLagrangianAlgorithm<double>::writeOutput(std::ostream &, const bool) const --> void", pybind11::arg("os"), pybind11::arg("print_header"));
		cl.def("setStatusTest", [](ROL::TypeG::Algorithm<double> &o, const class Teuchos::RCP<class ROL::StatusTest<double> > & a0) -> void { return o.setStatusTest(a0); }, "", pybind11::arg("status"));
		cl.def("setStatusTest", (void (ROL::TypeG::Algorithm<double>::*)(const class Teuchos::RCP<class ROL::StatusTest<double> > &, bool)) &ROL::TypeG::Algorithm<double>::setStatusTest, "C++: ROL::TypeG::Algorithm<double>::setStatusTest(const class Teuchos::RCP<class ROL::StatusTest<double> > &, bool) --> void", pybind11::arg("status"), pybind11::arg("combineStatus"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Problem<double> & a0) -> void { return o.run(a0); }, "", pybind11::arg("problem"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Problem<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Problem<double> &, std::ostream &) --> void", pybind11::arg("problem"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4) -> void { return o.run(a0, a1, a2, a3, a4); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4) -> void { return o.run(a0, a1, a2, a3, a4); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5) -> void { return o.run(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9, const class ROL::Vector<double> & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9, class ROL::Constraint<double> & a10, class ROL::Vector<double> & a11, const class ROL::Vector<double> & a12) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10, class ROL::Constraint<double> & a11, class ROL::Vector<double> & a12, const class ROL::Vector<double> & a13) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("writeHeader", (void (ROL::TypeG::Algorithm<double>::*)(std::ostream &) const) &ROL::TypeG::Algorithm<double>::writeHeader, "C++: ROL::TypeG::Algorithm<double>::writeHeader(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("writeName", (void (ROL::TypeG::Algorithm<double>::*)(std::ostream &) const) &ROL::TypeG::Algorithm<double>::writeName, "C++: ROL::TypeG::Algorithm<double>::writeName(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("writeOutput", [](ROL::TypeG::Algorithm<double> const &o, std::ostream & a0) -> void { return o.writeOutput(a0); }, "", pybind11::arg("os"));
		cl.def("writeOutput", (void (ROL::TypeG::Algorithm<double>::*)(std::ostream &, const bool) const) &ROL::TypeG::Algorithm<double>::writeOutput, "C++: ROL::TypeG::Algorithm<double>::writeOutput(std::ostream &, const bool) const --> void", pybind11::arg("os"), pybind11::arg("write_header"));
		cl.def("writeExitStatus", (void (ROL::TypeG::Algorithm<double>::*)(std::ostream &) const) &ROL::TypeG::Algorithm<double>::writeExitStatus, "C++: ROL::TypeG::Algorithm<double>::writeExitStatus(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("getState", (class Teuchos::RCP<const struct ROL::TypeG::AlgorithmState<double> > (ROL::TypeG::Algorithm<double>::*)() const) &ROL::TypeG::Algorithm<double>::getState, "C++: ROL::TypeG::Algorithm<double>::getState() const --> class Teuchos::RCP<const struct ROL::TypeG::AlgorithmState<double> >");
		cl.def("reset", (void (ROL::TypeG::Algorithm<double>::*)()) &ROL::TypeG::Algorithm<double>::reset, "C++: ROL::TypeG::Algorithm<double>::reset() --> void");
	}
	{ // ROL::TypeG::MoreauYosidaAlgorithm file:ROL_TypeG_MoreauYosidaAlgorithm.hpp line:58
		pybind11::class_<ROL::TypeG::MoreauYosidaAlgorithm<double>, Teuchos::RCP<ROL::TypeG::MoreauYosidaAlgorithm<double>>, PyCallBack_ROL_TypeG_MoreauYosidaAlgorithm_double_t, ROL::TypeG::Algorithm<double>> cl(M("ROL::TypeG"), "MoreauYosidaAlgorithm_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("list") );

		cl.def( pybind11::init( [](PyCallBack_ROL_TypeG_MoreauYosidaAlgorithm_double_t const &o){ return new PyCallBack_ROL_TypeG_MoreauYosidaAlgorithm_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::TypeG::MoreauYosidaAlgorithm<double> const &o){ return new ROL::TypeG::MoreauYosidaAlgorithm<double>(o); } ) );
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10, class ROL::Constraint<double> & a11, class ROL::Vector<double> & a12, const class ROL::Vector<double> & a13) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10, class ROL::Constraint<double> & a11, class ROL::Vector<double> & a12, const class ROL::Vector<double> & a13, std::ostream & a14) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9, class ROL::Constraint<double> & a10, class ROL::Vector<double> & a11, const class ROL::Vector<double> & a12) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9, class ROL::Constraint<double> & a10, class ROL::Vector<double> & a11, const class ROL::Vector<double> & a12, std::ostream & a13) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9, const class ROL::Vector<double> & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9, const class ROL::Vector<double> & a10, std::ostream & a11) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9, std::ostream & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9, std::ostream & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9, std::ostream & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, std::ostream & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, std::ostream & a8) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, std::ostream & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, std::ostream & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10, std::ostream & a11) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9, std::ostream & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7, std::ostream & a8) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6, std::ostream & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7, std::ostream & a8) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, std::ostream & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5) -> void { return o.run(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, std::ostream & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4) -> void { return o.run(a0, a1, a2, a3, a4); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4, std::ostream & a5) -> void { return o.run(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4) -> void { return o.run(a0, a1, a2, a3, a4); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, std::ostream & a5) -> void { return o.run(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Problem<double> & a0) -> void { return o.run(a0); }, "", pybind11::arg("problem"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Problem<double> & a0, std::ostream & a1) -> void { return o.run(a0, a1); }, "", pybind11::arg("problem"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::MoreauYosidaAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"));
		cl.def("run", (void (ROL::TypeG::MoreauYosidaAlgorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::MoreauYosidaAlgorithm<double>::run, "C++: ROL::TypeG::MoreauYosidaAlgorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("outStream"));
		cl.def("writeHeader", (void (ROL::TypeG::MoreauYosidaAlgorithm<double>::*)(std::ostream &) const) &ROL::TypeG::MoreauYosidaAlgorithm<double>::writeHeader, "C++: ROL::TypeG::MoreauYosidaAlgorithm<double>::writeHeader(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("writeName", (void (ROL::TypeG::MoreauYosidaAlgorithm<double>::*)(std::ostream &) const) &ROL::TypeG::MoreauYosidaAlgorithm<double>::writeName, "C++: ROL::TypeG::MoreauYosidaAlgorithm<double>::writeName(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("writeOutput", [](ROL::TypeG::MoreauYosidaAlgorithm<double> const &o, std::ostream & a0) -> void { return o.writeOutput(a0); }, "", pybind11::arg("os"));
		cl.def("writeOutput", (void (ROL::TypeG::MoreauYosidaAlgorithm<double>::*)(std::ostream &, const bool) const) &ROL::TypeG::MoreauYosidaAlgorithm<double>::writeOutput, "C++: ROL::TypeG::MoreauYosidaAlgorithm<double>::writeOutput(std::ostream &, const bool) const --> void", pybind11::arg("os"), pybind11::arg("print_header"));
		cl.def("setStatusTest", [](ROL::TypeG::Algorithm<double> &o, const class Teuchos::RCP<class ROL::StatusTest<double> > & a0) -> void { return o.setStatusTest(a0); }, "", pybind11::arg("status"));
		cl.def("setStatusTest", (void (ROL::TypeG::Algorithm<double>::*)(const class Teuchos::RCP<class ROL::StatusTest<double> > &, bool)) &ROL::TypeG::Algorithm<double>::setStatusTest, "C++: ROL::TypeG::Algorithm<double>::setStatusTest(const class Teuchos::RCP<class ROL::StatusTest<double> > &, bool) --> void", pybind11::arg("status"), pybind11::arg("combineStatus"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Problem<double> & a0) -> void { return o.run(a0); }, "", pybind11::arg("problem"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Problem<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Problem<double> &, std::ostream &) --> void", pybind11::arg("problem"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4) -> void { return o.run(a0, a1, a2, a3, a4); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4) -> void { return o.run(a0, a1, a2, a3, a4); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5) -> void { return o.run(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9, const class ROL::Vector<double> & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9, class ROL::Constraint<double> & a10, class ROL::Vector<double> & a11, const class ROL::Vector<double> & a12) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10, class ROL::Constraint<double> & a11, class ROL::Vector<double> & a12, const class ROL::Vector<double> & a13) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("writeHeader", (void (ROL::TypeG::Algorithm<double>::*)(std::ostream &) const) &ROL::TypeG::Algorithm<double>::writeHeader, "C++: ROL::TypeG::Algorithm<double>::writeHeader(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("writeName", (void (ROL::TypeG::Algorithm<double>::*)(std::ostream &) const) &ROL::TypeG::Algorithm<double>::writeName, "C++: ROL::TypeG::Algorithm<double>::writeName(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("writeOutput", [](ROL::TypeG::Algorithm<double> const &o, std::ostream & a0) -> void { return o.writeOutput(a0); }, "", pybind11::arg("os"));
		cl.def("writeOutput", (void (ROL::TypeG::Algorithm<double>::*)(std::ostream &, const bool) const) &ROL::TypeG::Algorithm<double>::writeOutput, "C++: ROL::TypeG::Algorithm<double>::writeOutput(std::ostream &, const bool) const --> void", pybind11::arg("os"), pybind11::arg("write_header"));
		cl.def("writeExitStatus", (void (ROL::TypeG::Algorithm<double>::*)(std::ostream &) const) &ROL::TypeG::Algorithm<double>::writeExitStatus, "C++: ROL::TypeG::Algorithm<double>::writeExitStatus(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("getState", (class Teuchos::RCP<const struct ROL::TypeG::AlgorithmState<double> > (ROL::TypeG::Algorithm<double>::*)() const) &ROL::TypeG::Algorithm<double>::getState, "C++: ROL::TypeG::Algorithm<double>::getState() const --> class Teuchos::RCP<const struct ROL::TypeG::AlgorithmState<double> >");
		cl.def("reset", (void (ROL::TypeG::Algorithm<double>::*)()) &ROL::TypeG::Algorithm<double>::reset, "C++: ROL::TypeG::Algorithm<double>::reset() --> void");
	}
	{ // ROL::TypeG::InteriorPointAlgorithm file:ROL_TypeG_InteriorPointAlgorithm.hpp line:58
		pybind11::class_<ROL::TypeG::InteriorPointAlgorithm<double>, Teuchos::RCP<ROL::TypeG::InteriorPointAlgorithm<double>>, PyCallBack_ROL_TypeG_InteriorPointAlgorithm_double_t, ROL::TypeG::Algorithm<double>> cl(M("ROL::TypeG"), "InteriorPointAlgorithm_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("list") );

		cl.def( pybind11::init( [](PyCallBack_ROL_TypeG_InteriorPointAlgorithm_double_t const &o){ return new PyCallBack_ROL_TypeG_InteriorPointAlgorithm_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::TypeG::InteriorPointAlgorithm<double> const &o){ return new ROL::TypeG::InteriorPointAlgorithm<double>(o); } ) );
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10, class ROL::Constraint<double> & a11, class ROL::Vector<double> & a12, const class ROL::Vector<double> & a13) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10, class ROL::Constraint<double> & a11, class ROL::Vector<double> & a12, const class ROL::Vector<double> & a13, std::ostream & a14) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9, class ROL::Constraint<double> & a10, class ROL::Vector<double> & a11, const class ROL::Vector<double> & a12) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9, class ROL::Constraint<double> & a10, class ROL::Vector<double> & a11, const class ROL::Vector<double> & a12, std::ostream & a13) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9, const class ROL::Vector<double> & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9, const class ROL::Vector<double> & a10, std::ostream & a11) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9, std::ostream & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9, std::ostream & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9, std::ostream & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, std::ostream & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, std::ostream & a8) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, std::ostream & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, std::ostream & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10, std::ostream & a11) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9, std::ostream & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7, std::ostream & a8) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6, std::ostream & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7, std::ostream & a8) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, std::ostream & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5) -> void { return o.run(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, std::ostream & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4) -> void { return o.run(a0, a1, a2, a3, a4); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4, std::ostream & a5) -> void { return o.run(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4) -> void { return o.run(a0, a1, a2, a3, a4); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, std::ostream & a5) -> void { return o.run(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Problem<double> & a0) -> void { return o.run(a0); }, "", pybind11::arg("problem"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Problem<double> & a0, std::ostream & a1) -> void { return o.run(a0, a1); }, "", pybind11::arg("problem"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::InteriorPointAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"));
		cl.def("run", (void (ROL::TypeG::InteriorPointAlgorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::InteriorPointAlgorithm<double>::run, "C++: ROL::TypeG::InteriorPointAlgorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("outStream"));
		cl.def("writeHeader", (void (ROL::TypeG::InteriorPointAlgorithm<double>::*)(std::ostream &) const) &ROL::TypeG::InteriorPointAlgorithm<double>::writeHeader, "C++: ROL::TypeG::InteriorPointAlgorithm<double>::writeHeader(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("writeName", (void (ROL::TypeG::InteriorPointAlgorithm<double>::*)(std::ostream &) const) &ROL::TypeG::InteriorPointAlgorithm<double>::writeName, "C++: ROL::TypeG::InteriorPointAlgorithm<double>::writeName(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("writeOutput", [](ROL::TypeG::InteriorPointAlgorithm<double> const &o, std::ostream & a0) -> void { return o.writeOutput(a0); }, "", pybind11::arg("os"));
		cl.def("writeOutput", (void (ROL::TypeG::InteriorPointAlgorithm<double>::*)(std::ostream &, const bool) const) &ROL::TypeG::InteriorPointAlgorithm<double>::writeOutput, "C++: ROL::TypeG::InteriorPointAlgorithm<double>::writeOutput(std::ostream &, const bool) const --> void", pybind11::arg("os"), pybind11::arg("print_header"));
		cl.def("setStatusTest", [](ROL::TypeG::Algorithm<double> &o, const class Teuchos::RCP<class ROL::StatusTest<double> > & a0) -> void { return o.setStatusTest(a0); }, "", pybind11::arg("status"));
		cl.def("setStatusTest", (void (ROL::TypeG::Algorithm<double>::*)(const class Teuchos::RCP<class ROL::StatusTest<double> > &, bool)) &ROL::TypeG::Algorithm<double>::setStatusTest, "C++: ROL::TypeG::Algorithm<double>::setStatusTest(const class Teuchos::RCP<class ROL::StatusTest<double> > &, bool) --> void", pybind11::arg("status"), pybind11::arg("combineStatus"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Problem<double> & a0) -> void { return o.run(a0); }, "", pybind11::arg("problem"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Problem<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Problem<double> &, std::ostream &) --> void", pybind11::arg("problem"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4) -> void { return o.run(a0, a1, a2, a3, a4); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4) -> void { return o.run(a0, a1, a2, a3, a4); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5) -> void { return o.run(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9, const class ROL::Vector<double> & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9, class ROL::Constraint<double> & a10, class ROL::Vector<double> & a11, const class ROL::Vector<double> & a12) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10, class ROL::Constraint<double> & a11, class ROL::Vector<double> & a12, const class ROL::Vector<double> & a13) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("writeHeader", (void (ROL::TypeG::Algorithm<double>::*)(std::ostream &) const) &ROL::TypeG::Algorithm<double>::writeHeader, "C++: ROL::TypeG::Algorithm<double>::writeHeader(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("writeName", (void (ROL::TypeG::Algorithm<double>::*)(std::ostream &) const) &ROL::TypeG::Algorithm<double>::writeName, "C++: ROL::TypeG::Algorithm<double>::writeName(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("writeOutput", [](ROL::TypeG::Algorithm<double> const &o, std::ostream & a0) -> void { return o.writeOutput(a0); }, "", pybind11::arg("os"));
		cl.def("writeOutput", (void (ROL::TypeG::Algorithm<double>::*)(std::ostream &, const bool) const) &ROL::TypeG::Algorithm<double>::writeOutput, "C++: ROL::TypeG::Algorithm<double>::writeOutput(std::ostream &, const bool) const --> void", pybind11::arg("os"), pybind11::arg("write_header"));
		cl.def("writeExitStatus", (void (ROL::TypeG::Algorithm<double>::*)(std::ostream &) const) &ROL::TypeG::Algorithm<double>::writeExitStatus, "C++: ROL::TypeG::Algorithm<double>::writeExitStatus(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("getState", (class Teuchos::RCP<const struct ROL::TypeG::AlgorithmState<double> > (ROL::TypeG::Algorithm<double>::*)() const) &ROL::TypeG::Algorithm<double>::getState, "C++: ROL::TypeG::Algorithm<double>::getState() const --> class Teuchos::RCP<const struct ROL::TypeG::AlgorithmState<double> >");
		cl.def("reset", (void (ROL::TypeG::Algorithm<double>::*)()) &ROL::TypeG::Algorithm<double>::reset, "C++: ROL::TypeG::Algorithm<double>::reset() --> void");
	}
	{ // ROL::TypeG::StabilizedLCLAlgorithm file:ROL_TypeG_StabilizedLCLAlgorithm.hpp line:61
		pybind11::class_<ROL::TypeG::StabilizedLCLAlgorithm<double>, Teuchos::RCP<ROL::TypeG::StabilizedLCLAlgorithm<double>>, PyCallBack_ROL_TypeG_StabilizedLCLAlgorithm_double_t, ROL::TypeG::Algorithm<double>> cl(M("ROL::TypeG"), "StabilizedLCLAlgorithm_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("list") );

		cl.def( pybind11::init( [](PyCallBack_ROL_TypeG_StabilizedLCLAlgorithm_double_t const &o){ return new PyCallBack_ROL_TypeG_StabilizedLCLAlgorithm_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::TypeG::StabilizedLCLAlgorithm<double> const &o){ return new ROL::TypeG::StabilizedLCLAlgorithm<double>(o); } ) );
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10, class ROL::Constraint<double> & a11, class ROL::Vector<double> & a12, const class ROL::Vector<double> & a13) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10, class ROL::Constraint<double> & a11, class ROL::Vector<double> & a12, const class ROL::Vector<double> & a13, std::ostream & a14) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9, class ROL::Constraint<double> & a10, class ROL::Vector<double> & a11, const class ROL::Vector<double> & a12) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9, class ROL::Constraint<double> & a10, class ROL::Vector<double> & a11, const class ROL::Vector<double> & a12, std::ostream & a13) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9, const class ROL::Vector<double> & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9, const class ROL::Vector<double> & a10, std::ostream & a11) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9, std::ostream & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9, std::ostream & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9, std::ostream & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, std::ostream & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, std::ostream & a8) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, std::ostream & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, std::ostream & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10, std::ostream & a11) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9, std::ostream & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7, std::ostream & a8) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6, std::ostream & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7, std::ostream & a8) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, std::ostream & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5) -> void { return o.run(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, std::ostream & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4) -> void { return o.run(a0, a1, a2, a3, a4); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4, std::ostream & a5) -> void { return o.run(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4) -> void { return o.run(a0, a1, a2, a3, a4); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, std::ostream & a5) -> void { return o.run(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Problem<double> & a0) -> void { return o.run(a0); }, "", pybind11::arg("problem"));
		cl.def("run", (void (ROL::TypeG::StabilizedLCLAlgorithm<double>::*)(class ROL::Problem<double> &, std::ostream &)) &ROL::TypeG::StabilizedLCLAlgorithm<double>::run, "C++: ROL::TypeG::StabilizedLCLAlgorithm<double>::run(class ROL::Problem<double> &, std::ostream &) --> void", pybind11::arg("problem"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::StabilizedLCLAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"));
		cl.def("run", (void (ROL::TypeG::StabilizedLCLAlgorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::StabilizedLCLAlgorithm<double>::run, "C++: ROL::TypeG::StabilizedLCLAlgorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("outStream"));
		cl.def("writeHeader", (void (ROL::TypeG::StabilizedLCLAlgorithm<double>::*)(std::ostream &) const) &ROL::TypeG::StabilizedLCLAlgorithm<double>::writeHeader, "C++: ROL::TypeG::StabilizedLCLAlgorithm<double>::writeHeader(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("writeName", (void (ROL::TypeG::StabilizedLCLAlgorithm<double>::*)(std::ostream &) const) &ROL::TypeG::StabilizedLCLAlgorithm<double>::writeName, "C++: ROL::TypeG::StabilizedLCLAlgorithm<double>::writeName(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("writeOutput", [](ROL::TypeG::StabilizedLCLAlgorithm<double> const &o, std::ostream & a0) -> void { return o.writeOutput(a0); }, "", pybind11::arg("os"));
		cl.def("writeOutput", (void (ROL::TypeG::StabilizedLCLAlgorithm<double>::*)(std::ostream &, const bool) const) &ROL::TypeG::StabilizedLCLAlgorithm<double>::writeOutput, "C++: ROL::TypeG::StabilizedLCLAlgorithm<double>::writeOutput(std::ostream &, const bool) const --> void", pybind11::arg("os"), pybind11::arg("print_header"));
		cl.def("setStatusTest", [](ROL::TypeG::Algorithm<double> &o, const class Teuchos::RCP<class ROL::StatusTest<double> > & a0) -> void { return o.setStatusTest(a0); }, "", pybind11::arg("status"));
		cl.def("setStatusTest", (void (ROL::TypeG::Algorithm<double>::*)(const class Teuchos::RCP<class ROL::StatusTest<double> > &, bool)) &ROL::TypeG::Algorithm<double>::setStatusTest, "C++: ROL::TypeG::Algorithm<double>::setStatusTest(const class Teuchos::RCP<class ROL::StatusTest<double> > &, bool) --> void", pybind11::arg("status"), pybind11::arg("combineStatus"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Problem<double> & a0) -> void { return o.run(a0); }, "", pybind11::arg("problem"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Problem<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Problem<double> &, std::ostream &) --> void", pybind11::arg("problem"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4) -> void { return o.run(a0, a1, a2, a3, a4); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4) -> void { return o.run(a0, a1, a2, a3, a4); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5) -> void { return o.run(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::BoundConstraint<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::Constraint<double> & a2, class ROL::Vector<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, const class ROL::Vector<double> & a9) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7, class ROL::Constraint<double> & a8, class ROL::Vector<double> & a9, const class ROL::Vector<double> & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Constraint<double> & a6, class ROL::Vector<double> & a7, class ROL::BoundConstraint<double> & a8, const class ROL::Vector<double> & a9, class ROL::Constraint<double> & a10, class ROL::Vector<double> & a11, const class ROL::Vector<double> & a12) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeG::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10, class ROL::Constraint<double> & a11, class ROL::Vector<double> & a12, const class ROL::Vector<double> & a13) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", (void (ROL::TypeG::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeG::Algorithm<double>::run, "C++: ROL::TypeG::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("writeHeader", (void (ROL::TypeG::Algorithm<double>::*)(std::ostream &) const) &ROL::TypeG::Algorithm<double>::writeHeader, "C++: ROL::TypeG::Algorithm<double>::writeHeader(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("writeName", (void (ROL::TypeG::Algorithm<double>::*)(std::ostream &) const) &ROL::TypeG::Algorithm<double>::writeName, "C++: ROL::TypeG::Algorithm<double>::writeName(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("writeOutput", [](ROL::TypeG::Algorithm<double> const &o, std::ostream & a0) -> void { return o.writeOutput(a0); }, "", pybind11::arg("os"));
		cl.def("writeOutput", (void (ROL::TypeG::Algorithm<double>::*)(std::ostream &, const bool) const) &ROL::TypeG::Algorithm<double>::writeOutput, "C++: ROL::TypeG::Algorithm<double>::writeOutput(std::ostream &, const bool) const --> void", pybind11::arg("os"), pybind11::arg("write_header"));
		cl.def("writeExitStatus", (void (ROL::TypeG::Algorithm<double>::*)(std::ostream &) const) &ROL::TypeG::Algorithm<double>::writeExitStatus, "C++: ROL::TypeG::Algorithm<double>::writeExitStatus(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("getState", (class Teuchos::RCP<const struct ROL::TypeG::AlgorithmState<double> > (ROL::TypeG::Algorithm<double>::*)() const) &ROL::TypeG::Algorithm<double>::getState, "C++: ROL::TypeG::Algorithm<double>::getState() const --> class Teuchos::RCP<const struct ROL::TypeG::AlgorithmState<double> >");
		cl.def("reset", (void (ROL::TypeG::Algorithm<double>::*)()) &ROL::TypeG::Algorithm<double>::reset, "C++: ROL::TypeG::Algorithm<double>::reset() --> void");
	}
	// ROL::TypeG::EAlgorithmG file:ROL_TypeG_AlgorithmFactory.hpp line:64
	pybind11::enum_<ROL::TypeG::EAlgorithmG>(M("ROL::TypeG"), "EAlgorithmG", pybind11::arithmetic(), "Enumeration of generally constrained algorithm types.\n\n    \n    ALGORITHM_G_AUGMENTEDLAGRANGIAN describe\n    \n\n    ALGORITHM_G_MOREAUYOSIDA        describe\n    \n\n    ALGORITHM_G_INTERIORPOINT       describe\n    \n\n    ALGORITHM_G_STABILIZEDLCL       describe", pybind11::module_local())
		.value("ALGORITHM_G_AUGMENTEDLAGRANGIAN", ROL::TypeG::ALGORITHM_G_AUGMENTEDLAGRANGIAN)
		.value("ALGORITHM_G_MOREAUYOSIDA", ROL::TypeG::ALGORITHM_G_MOREAUYOSIDA)
		.value("ALGORITHM_G_INTERIORPOINT", ROL::TypeG::ALGORITHM_G_INTERIORPOINT)
		.value("ALGORITHM_G_STABILIZEDLCL", ROL::TypeG::ALGORITHM_G_STABILIZEDLCL)
		.value("ALGORITHM_G_LAST", ROL::TypeG::ALGORITHM_G_LAST)
		.export_values();

;

	// ROL::TypeG::EAlgorithmGToString(enum ROL::TypeG::EAlgorithmG) file:ROL_TypeG_AlgorithmFactory.hpp line:72
	M("ROL::TypeG").def("EAlgorithmGToString", (std::string (*)(enum ROL::TypeG::EAlgorithmG)) &ROL::TypeG::EAlgorithmGToString, "C++: ROL::TypeG::EAlgorithmGToString(enum ROL::TypeG::EAlgorithmG) --> std::string", pybind11::arg("alg"));

	// ROL::TypeG::isValidAlgorithmG(enum ROL::TypeG::EAlgorithmG) file:ROL_TypeG_AlgorithmFactory.hpp line:90
	M("ROL::TypeG").def("isValidAlgorithmG", (int (*)(enum ROL::TypeG::EAlgorithmG)) &ROL::TypeG::isValidAlgorithmG, "Verifies validity of a AlgorithmG enum.\n\n    \n  [in]  - enum of the AlgorithmG\n    \n\n 1 if the argument is a valid AlgorithmG; 0 otherwise.\n\nC++: ROL::TypeG::isValidAlgorithmG(enum ROL::TypeG::EAlgorithmG) --> int", pybind11::arg("alg"));

	// ROL::TypeG::StringToEAlgorithmG(std::string) file:ROL_TypeG_AlgorithmFactory.hpp line:119
	M("ROL::TypeG").def("StringToEAlgorithmG", (enum ROL::TypeG::EAlgorithmG (*)(std::string)) &ROL::TypeG::StringToEAlgorithmG, "C++: ROL::TypeG::StringToEAlgorithmG(std::string) --> enum ROL::TypeG::EAlgorithmG", pybind11::arg("s"));

	// ROL::TypeG::AlgorithmFactory(class Teuchos::ParameterList &) file:ROL_TypeG_AlgorithmFactory.hpp line:130
	M("ROL::TypeG").def("AlgorithmFactory", (class Teuchos::RCP<class ROL::TypeG::Algorithm<double> > (*)(class Teuchos::ParameterList &)) &ROL::TypeG::AlgorithmFactory<double>, "C++: ROL::TypeG::AlgorithmFactory(class Teuchos::ParameterList &) --> class Teuchos::RCP<class ROL::TypeG::Algorithm<double> >", pybind11::arg("parlist"));

}
