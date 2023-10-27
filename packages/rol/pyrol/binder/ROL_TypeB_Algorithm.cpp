#include <ROL_BoundConstraint.hpp>
#include <ROL_Constraint.hpp>
#include <ROL_Elementwise_Function.hpp>
#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_Krylov.hpp>
#include <ROL_Objective.hpp>
#include <ROL_PolyhedralProjection.hpp>
#include <ROL_Problem.hpp>
#include <ROL_Secant.hpp>
#include <ROL_StatusTest.hpp>
#include <ROL_TypeB_Algorithm.hpp>
#include <ROL_TypeB_GradientAlgorithm.hpp>
#include <ROL_TypeB_NewtonKrylovAlgorithm.hpp>
#include <ROL_Types.hpp>
#include <ROL_UpdateType.hpp>
#include <ROL_Vector.hpp>
#include <Teuchos_ENull.hpp>
#include <Teuchos_FilteredIterator.hpp>
#include <Teuchos_ParameterEntry.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ParameterListModifier.hpp>
#include <Teuchos_RCPDecl.hpp>
#include <Teuchos_RCPNode.hpp>
#include <Teuchos_StringIndexedOrderedValueObjectContainer.hpp>
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

// ROL::TypeB::Algorithm file:ROL_TypeB_Algorithm.hpp line:84
struct PyCallBack_ROL_TypeB_Algorithm_double_t : public ROL::TypeB::Algorithm<double> {
	using ROL::TypeB::Algorithm<double>::Algorithm;

	void run(class ROL::Problem<double> & a0, std::ostream & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::Algorithm<double> *>(this), "run");
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
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, std::ostream & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::Algorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, std::ostream & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::Algorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"Algorithm::run\"");
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, std::ostream & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::Algorithm<double> *>(this), "run");
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
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, std::ostream & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::Algorithm<double> *>(this), "run");
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
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, std::ostream & a6) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::Algorithm<double> *>(this), "run");
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
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7, std::ostream & a8) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::Algorithm<double> *>(this), "run");
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
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7, std::ostream & a8) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::Algorithm<double> *>(this), "run");
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
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10, std::ostream & a11) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::Algorithm<double> *>(this), "run");
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
	void writeHeader(std::ostream & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::Algorithm<double> *>(this), "writeHeader");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::Algorithm<double> *>(this), "writeName");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::Algorithm<double> *>(this), "writeOutput");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::Algorithm<double> *>(this), "writeExitStatus");
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

// ROL::TypeB::GradientAlgorithm file:ROL_TypeB_GradientAlgorithm.hpp line:57
struct PyCallBack_ROL_TypeB_GradientAlgorithm_double_t : public ROL::TypeB::GradientAlgorithm<double> {
	using ROL::TypeB::GradientAlgorithm<double>::GradientAlgorithm;

	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, std::ostream & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::GradientAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return GradientAlgorithm::run(a0, a1, a2, a3, a4);
	}
	void writeHeader(std::ostream & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::GradientAlgorithm<double> *>(this), "writeHeader");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return GradientAlgorithm::writeHeader(a0);
	}
	void writeName(std::ostream & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::GradientAlgorithm<double> *>(this), "writeName");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return GradientAlgorithm::writeName(a0);
	}
	void writeOutput(std::ostream & a0, const bool a1) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::GradientAlgorithm<double> *>(this), "writeOutput");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return GradientAlgorithm::writeOutput(a0, a1);
	}
	void run(class ROL::Problem<double> & a0, std::ostream & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::GradientAlgorithm<double> *>(this), "run");
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
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, std::ostream & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::GradientAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, std::ostream & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::GradientAlgorithm<double> *>(this), "run");
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
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, std::ostream & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::GradientAlgorithm<double> *>(this), "run");
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
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, std::ostream & a6) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::GradientAlgorithm<double> *>(this), "run");
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
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7, std::ostream & a8) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::GradientAlgorithm<double> *>(this), "run");
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
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7, std::ostream & a8) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::GradientAlgorithm<double> *>(this), "run");
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
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10, std::ostream & a11) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::GradientAlgorithm<double> *>(this), "run");
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
	void writeExitStatus(std::ostream & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::GradientAlgorithm<double> *>(this), "writeExitStatus");
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

// ROL::TypeB::NewtonKrylovAlgorithm file:ROL_TypeB_NewtonKrylovAlgorithm.hpp line:59
struct PyCallBack_ROL_TypeB_NewtonKrylovAlgorithm_double_t : public ROL::TypeB::NewtonKrylovAlgorithm<double> {
	using ROL::TypeB::NewtonKrylovAlgorithm<double>::NewtonKrylovAlgorithm;

	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, std::ostream & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::NewtonKrylovAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NewtonKrylovAlgorithm::run(a0, a1, a2, a3, a4);
	}
	void run(class ROL::Problem<double> & a0, std::ostream & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::NewtonKrylovAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NewtonKrylovAlgorithm::run(a0, a1);
	}
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, std::ostream & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::NewtonKrylovAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NewtonKrylovAlgorithm::run(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	void writeHeader(std::ostream & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::NewtonKrylovAlgorithm<double> *>(this), "writeHeader");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NewtonKrylovAlgorithm::writeHeader(a0);
	}
	void writeName(std::ostream & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::NewtonKrylovAlgorithm<double> *>(this), "writeName");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NewtonKrylovAlgorithm::writeName(a0);
	}
	void writeOutput(std::ostream & a0, const bool a1) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::NewtonKrylovAlgorithm<double> *>(this), "writeOutput");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NewtonKrylovAlgorithm::writeOutput(a0, a1);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, std::ostream & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::NewtonKrylovAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run(a0, a1, a2, a3);
	}
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, std::ostream & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::NewtonKrylovAlgorithm<double> *>(this), "run");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::NewtonKrylovAlgorithm<double> *>(this), "run");
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
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7, std::ostream & a8) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::NewtonKrylovAlgorithm<double> *>(this), "run");
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
	void run(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7, std::ostream & a8) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::NewtonKrylovAlgorithm<double> *>(this), "run");
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
	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10, std::ostream & a11) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::NewtonKrylovAlgorithm<double> *>(this), "run");
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
	void writeExitStatus(std::ostream & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::NewtonKrylovAlgorithm<double> *>(this), "writeExitStatus");
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

void bind_ROL_TypeB_Algorithm(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // ROL::TypeB::AlgorithmState file:ROL_TypeB_Algorithm.hpp line:62
		pybind11::class_<ROL::TypeB::AlgorithmState<double>, Teuchos::RCP<ROL::TypeB::AlgorithmState<double>>, ROL::AlgorithmState<double>> cl(M("ROL::TypeB"), "AlgorithmState_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new ROL::TypeB::AlgorithmState<double>(); } ) );
		cl.def( pybind11::init( [](ROL::TypeB::AlgorithmState<double> const &o){ return new ROL::TypeB::AlgorithmState<double>(o); } ) );
		cl.def_readwrite("searchSize", &ROL::TypeB::AlgorithmState<double>::searchSize);
		cl.def_readwrite("stepVec", &ROL::TypeB::AlgorithmState<double>::stepVec);
		cl.def_readwrite("gradientVec", &ROL::TypeB::AlgorithmState<double>::gradientVec);
		cl.def_readwrite("nproj", &ROL::TypeB::AlgorithmState<double>::nproj);
		cl.def("reset", (void (ROL::TypeB::AlgorithmState<double>::*)()) &ROL::TypeB::AlgorithmState<double>::reset, "C++: ROL::TypeB::AlgorithmState<double>::reset() --> void");
		cl.def("assign", (struct ROL::TypeB::AlgorithmState<double> & (ROL::TypeB::AlgorithmState<double>::*)(const struct ROL::TypeB::AlgorithmState<double> &)) &ROL::TypeB::AlgorithmState<double>::operator=, "C++: ROL::TypeB::AlgorithmState<double>::operator=(const struct ROL::TypeB::AlgorithmState<double> &) --> struct ROL::TypeB::AlgorithmState<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
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
	{ // ROL::TypeB::Algorithm file:ROL_TypeB_Algorithm.hpp line:84
		pybind11::class_<ROL::TypeB::Algorithm<double>, Teuchos::RCP<ROL::TypeB::Algorithm<double>>, PyCallBack_ROL_TypeB_Algorithm_double_t> cl(M("ROL::TypeB"), "Algorithm_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new PyCallBack_ROL_TypeB_Algorithm_double_t(); } ) );
		cl.def(pybind11::init<PyCallBack_ROL_TypeB_Algorithm_double_t const &>());
		cl.def("setStatusTest", [](ROL::TypeB::Algorithm<double> &o, const class Teuchos::RCP<class ROL::StatusTest<double> > & a0) -> void { return o.setStatusTest(a0); }, "", pybind11::arg("status"));
		cl.def("setStatusTest", (void (ROL::TypeB::Algorithm<double>::*)(const class Teuchos::RCP<class ROL::StatusTest<double> > &, const bool)) &ROL::TypeB::Algorithm<double>::setStatusTest, "C++: ROL::TypeB::Algorithm<double>::setStatusTest(const class Teuchos::RCP<class ROL::StatusTest<double> > &, const bool) --> void", pybind11::arg("status"), pybind11::arg("combineStatus"));
		cl.def("run", [](ROL::TypeB::Algorithm<double> &o, class ROL::Problem<double> & a0) -> void { return o.run(a0); }, "", pybind11::arg("problem"));
		cl.def("run", (void (ROL::TypeB::Algorithm<double>::*)(class ROL::Problem<double> &, std::ostream &)) &ROL::TypeB::Algorithm<double>::run, "C++: ROL::TypeB::Algorithm<double>::run(class ROL::Problem<double> &, std::ostream &) --> void", pybind11::arg("problem"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2) -> void { return o.run(a0, a1, a2); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"));
		cl.def("run", (void (ROL::TypeB::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeB::Algorithm<double>::run, "C++: ROL::TypeB::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3) -> void { return o.run(a0, a1, a2, a3); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"));
		cl.def("run", (void (ROL::TypeB::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeB::Algorithm<double>::run, "C++: ROL::TypeB::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4) -> void { return o.run(a0, a1, a2, a3, a4); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", (void (ROL::TypeB::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeB::Algorithm<double>::run, "C++: ROL::TypeB::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", (void (ROL::TypeB::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeB::Algorithm<double>::run, "C++: ROL::TypeB::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5) -> void { return o.run(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"));
		cl.def("run", (void (ROL::TypeB::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeB::Algorithm<double>::run, "C++: ROL::TypeB::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("linear_ires"));
		cl.def("run", (void (ROL::TypeB::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeB::Algorithm<double>::run, "C++: ROL::TypeB::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("linear_ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"));
		cl.def("run", (void (ROL::TypeB::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeB::Algorithm<double>::run, "C++: ROL::TypeB::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("linear_ires"));
		cl.def("run", (void (ROL::TypeB::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeB::Algorithm<double>::run, "C++: ROL::TypeB::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("linear_ires"), pybind11::arg("outStream"));
		cl.def("writeHeader", (void (ROL::TypeB::Algorithm<double>::*)(std::ostream &) const) &ROL::TypeB::Algorithm<double>::writeHeader, "C++: ROL::TypeB::Algorithm<double>::writeHeader(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("writeName", (void (ROL::TypeB::Algorithm<double>::*)(std::ostream &) const) &ROL::TypeB::Algorithm<double>::writeName, "C++: ROL::TypeB::Algorithm<double>::writeName(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("writeOutput", [](ROL::TypeB::Algorithm<double> const &o, std::ostream & a0) -> void { return o.writeOutput(a0); }, "", pybind11::arg("os"));
		cl.def("writeOutput", (void (ROL::TypeB::Algorithm<double>::*)(std::ostream &, const bool) const) &ROL::TypeB::Algorithm<double>::writeOutput, "C++: ROL::TypeB::Algorithm<double>::writeOutput(std::ostream &, const bool) const --> void", pybind11::arg("os"), pybind11::arg("write_header"));
		cl.def("writeExitStatus", (void (ROL::TypeB::Algorithm<double>::*)(std::ostream &) const) &ROL::TypeB::Algorithm<double>::writeExitStatus, "C++: ROL::TypeB::Algorithm<double>::writeExitStatus(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("getState", (class Teuchos::RCP<const struct ROL::TypeB::AlgorithmState<double> > (ROL::TypeB::Algorithm<double>::*)() const) &ROL::TypeB::Algorithm<double>::getState, "C++: ROL::TypeB::Algorithm<double>::getState() const --> class Teuchos::RCP<const struct ROL::TypeB::AlgorithmState<double> >");
		cl.def("reset", (void (ROL::TypeB::Algorithm<double>::*)()) &ROL::TypeB::Algorithm<double>::reset, "C++: ROL::TypeB::Algorithm<double>::reset() --> void");
	}
	{ // ROL::TypeB::GradientAlgorithm file:ROL_TypeB_GradientAlgorithm.hpp line:57
		pybind11::class_<ROL::TypeB::GradientAlgorithm<double>, Teuchos::RCP<ROL::TypeB::GradientAlgorithm<double>>, PyCallBack_ROL_TypeB_GradientAlgorithm_double_t, ROL::TypeB::Algorithm<double>> cl(M("ROL::TypeB"), "GradientAlgorithm_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("list") );

		cl.def( pybind11::init( [](PyCallBack_ROL_TypeB_GradientAlgorithm_double_t const &o){ return new PyCallBack_ROL_TypeB_GradientAlgorithm_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::TypeB::GradientAlgorithm<double> const &o){ return new ROL::TypeB::GradientAlgorithm<double>(o); } ) );
		cl.def("run", [](ROL::TypeB::GradientAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("linear_ires"));
		cl.def("run", [](ROL::TypeB::GradientAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10, std::ostream & a11) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("linear_ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::GradientAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"));
		cl.def("run", [](ROL::TypeB::GradientAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7, std::ostream & a8) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::GradientAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("linear_ires"));
		cl.def("run", [](ROL::TypeB::GradientAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7, std::ostream & a8) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("linear_ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::GradientAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5) -> void { return o.run(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"));
		cl.def("run", [](ROL::TypeB::GradientAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, std::ostream & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::GradientAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", [](ROL::TypeB::GradientAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, std::ostream & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::GradientAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4) -> void { return o.run(a0, a1, a2, a3, a4); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", [](ROL::TypeB::GradientAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, std::ostream & a5) -> void { return o.run(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::GradientAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2) -> void { return o.run(a0, a1, a2); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"));
		cl.def("run", [](ROL::TypeB::GradientAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, std::ostream & a3) -> void { return o.run(a0, a1, a2, a3); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::GradientAlgorithm<double> &o, class ROL::Problem<double> & a0) -> void { return o.run(a0); }, "", pybind11::arg("problem"));
		cl.def("run", [](ROL::TypeB::GradientAlgorithm<double> &o, class ROL::Problem<double> & a0, std::ostream & a1) -> void { return o.run(a0, a1); }, "", pybind11::arg("problem"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::GradientAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3) -> void { return o.run(a0, a1, a2, a3); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"));
		cl.def("run", (void (ROL::TypeB::GradientAlgorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeB::GradientAlgorithm<double>::run, "C++: ROL::TypeB::GradientAlgorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("outStream"));
		cl.def("writeHeader", (void (ROL::TypeB::GradientAlgorithm<double>::*)(std::ostream &) const) &ROL::TypeB::GradientAlgorithm<double>::writeHeader, "C++: ROL::TypeB::GradientAlgorithm<double>::writeHeader(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("writeName", (void (ROL::TypeB::GradientAlgorithm<double>::*)(std::ostream &) const) &ROL::TypeB::GradientAlgorithm<double>::writeName, "C++: ROL::TypeB::GradientAlgorithm<double>::writeName(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("writeOutput", [](ROL::TypeB::GradientAlgorithm<double> const &o, std::ostream & a0) -> void { return o.writeOutput(a0); }, "", pybind11::arg("os"));
		cl.def("writeOutput", (void (ROL::TypeB::GradientAlgorithm<double>::*)(std::ostream &, const bool) const) &ROL::TypeB::GradientAlgorithm<double>::writeOutput, "C++: ROL::TypeB::GradientAlgorithm<double>::writeOutput(std::ostream &, const bool) const --> void", pybind11::arg("os"), pybind11::arg("write_header"));
		cl.def("setStatusTest", [](ROL::TypeB::Algorithm<double> &o, const class Teuchos::RCP<class ROL::StatusTest<double> > & a0) -> void { return o.setStatusTest(a0); }, "", pybind11::arg("status"));
		cl.def("setStatusTest", (void (ROL::TypeB::Algorithm<double>::*)(const class Teuchos::RCP<class ROL::StatusTest<double> > &, const bool)) &ROL::TypeB::Algorithm<double>::setStatusTest, "C++: ROL::TypeB::Algorithm<double>::setStatusTest(const class Teuchos::RCP<class ROL::StatusTest<double> > &, const bool) --> void", pybind11::arg("status"), pybind11::arg("combineStatus"));
		cl.def("run", [](ROL::TypeB::Algorithm<double> &o, class ROL::Problem<double> & a0) -> void { return o.run(a0); }, "", pybind11::arg("problem"));
		cl.def("run", (void (ROL::TypeB::Algorithm<double>::*)(class ROL::Problem<double> &, std::ostream &)) &ROL::TypeB::Algorithm<double>::run, "C++: ROL::TypeB::Algorithm<double>::run(class ROL::Problem<double> &, std::ostream &) --> void", pybind11::arg("problem"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2) -> void { return o.run(a0, a1, a2); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"));
		cl.def("run", (void (ROL::TypeB::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeB::Algorithm<double>::run, "C++: ROL::TypeB::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3) -> void { return o.run(a0, a1, a2, a3); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"));
		cl.def("run", (void (ROL::TypeB::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeB::Algorithm<double>::run, "C++: ROL::TypeB::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4) -> void { return o.run(a0, a1, a2, a3, a4); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", (void (ROL::TypeB::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeB::Algorithm<double>::run, "C++: ROL::TypeB::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", (void (ROL::TypeB::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeB::Algorithm<double>::run, "C++: ROL::TypeB::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5) -> void { return o.run(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"));
		cl.def("run", (void (ROL::TypeB::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeB::Algorithm<double>::run, "C++: ROL::TypeB::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("linear_ires"));
		cl.def("run", (void (ROL::TypeB::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeB::Algorithm<double>::run, "C++: ROL::TypeB::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("linear_ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"));
		cl.def("run", (void (ROL::TypeB::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeB::Algorithm<double>::run, "C++: ROL::TypeB::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("linear_ires"));
		cl.def("run", (void (ROL::TypeB::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeB::Algorithm<double>::run, "C++: ROL::TypeB::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("linear_ires"), pybind11::arg("outStream"));
		cl.def("writeHeader", (void (ROL::TypeB::Algorithm<double>::*)(std::ostream &) const) &ROL::TypeB::Algorithm<double>::writeHeader, "C++: ROL::TypeB::Algorithm<double>::writeHeader(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("writeName", (void (ROL::TypeB::Algorithm<double>::*)(std::ostream &) const) &ROL::TypeB::Algorithm<double>::writeName, "C++: ROL::TypeB::Algorithm<double>::writeName(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("writeOutput", [](ROL::TypeB::Algorithm<double> const &o, std::ostream & a0) -> void { return o.writeOutput(a0); }, "", pybind11::arg("os"));
		cl.def("writeOutput", (void (ROL::TypeB::Algorithm<double>::*)(std::ostream &, const bool) const) &ROL::TypeB::Algorithm<double>::writeOutput, "C++: ROL::TypeB::Algorithm<double>::writeOutput(std::ostream &, const bool) const --> void", pybind11::arg("os"), pybind11::arg("write_header"));
		cl.def("writeExitStatus", (void (ROL::TypeB::Algorithm<double>::*)(std::ostream &) const) &ROL::TypeB::Algorithm<double>::writeExitStatus, "C++: ROL::TypeB::Algorithm<double>::writeExitStatus(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("getState", (class Teuchos::RCP<const struct ROL::TypeB::AlgorithmState<double> > (ROL::TypeB::Algorithm<double>::*)() const) &ROL::TypeB::Algorithm<double>::getState, "C++: ROL::TypeB::Algorithm<double>::getState() const --> class Teuchos::RCP<const struct ROL::TypeB::AlgorithmState<double> >");
		cl.def("reset", (void (ROL::TypeB::Algorithm<double>::*)()) &ROL::TypeB::Algorithm<double>::reset, "C++: ROL::TypeB::Algorithm<double>::reset() --> void");
	}
	{ // ROL::TypeB::NewtonKrylovAlgorithm file:ROL_TypeB_NewtonKrylovAlgorithm.hpp line:59
		pybind11::class_<ROL::TypeB::NewtonKrylovAlgorithm<double>, Teuchos::RCP<ROL::TypeB::NewtonKrylovAlgorithm<double>>, PyCallBack_ROL_TypeB_NewtonKrylovAlgorithm_double_t, ROL::TypeB::Algorithm<double>> cl(M("ROL::TypeB"), "NewtonKrylovAlgorithm_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](class Teuchos::ParameterList & a0){ return new ROL::TypeB::NewtonKrylovAlgorithm<double>(a0); }, [](class Teuchos::ParameterList & a0){ return new PyCallBack_ROL_TypeB_NewtonKrylovAlgorithm_double_t(a0); } ), "doc");
		cl.def( pybind11::init<class Teuchos::ParameterList &, const class Teuchos::RCP<class ROL::Secant<double> > &>(), pybind11::arg("list"), pybind11::arg("secant") );

		cl.def( pybind11::init( [](class Teuchos::ParameterList & a0, const class Teuchos::RCP<class ROL::Krylov<double> > & a1){ return new ROL::TypeB::NewtonKrylovAlgorithm<double>(a0, a1); }, [](class Teuchos::ParameterList & a0, const class Teuchos::RCP<class ROL::Krylov<double> > & a1){ return new PyCallBack_ROL_TypeB_NewtonKrylovAlgorithm_double_t(a0, a1); } ), "doc");
		cl.def( pybind11::init<class Teuchos::ParameterList &, const class Teuchos::RCP<class ROL::Krylov<double> > &, const class Teuchos::RCP<class ROL::Secant<double> > &>(), pybind11::arg("list"), pybind11::arg("krylov"), pybind11::arg("secant") );

		cl.def( pybind11::init( [](PyCallBack_ROL_TypeB_NewtonKrylovAlgorithm_double_t const &o){ return new PyCallBack_ROL_TypeB_NewtonKrylovAlgorithm_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::TypeB::NewtonKrylovAlgorithm<double> const &o){ return new ROL::TypeB::NewtonKrylovAlgorithm<double>(o); } ) );
		cl.def("run", [](ROL::TypeB::NewtonKrylovAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("linear_ires"));
		cl.def("run", [](ROL::TypeB::NewtonKrylovAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10, std::ostream & a11) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("linear_ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::NewtonKrylovAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"));
		cl.def("run", [](ROL::TypeB::NewtonKrylovAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7, std::ostream & a8) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::NewtonKrylovAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("linear_ires"));
		cl.def("run", [](ROL::TypeB::NewtonKrylovAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7, std::ostream & a8) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("linear_ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::NewtonKrylovAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5) -> void { return o.run(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"));
		cl.def("run", [](ROL::TypeB::NewtonKrylovAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, std::ostream & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::NewtonKrylovAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4) -> void { return o.run(a0, a1, a2, a3, a4); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", [](ROL::TypeB::NewtonKrylovAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, std::ostream & a5) -> void { return o.run(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::NewtonKrylovAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2) -> void { return o.run(a0, a1, a2); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"));
		cl.def("run", [](ROL::TypeB::NewtonKrylovAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, std::ostream & a3) -> void { return o.run(a0, a1, a2, a3); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::NewtonKrylovAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3) -> void { return o.run(a0, a1, a2, a3); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"));
		cl.def("run", (void (ROL::TypeB::NewtonKrylovAlgorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeB::NewtonKrylovAlgorithm<double>::run, "C++: ROL::TypeB::NewtonKrylovAlgorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::NewtonKrylovAlgorithm<double> &o, class ROL::Problem<double> & a0) -> void { return o.run(a0); }, "", pybind11::arg("problem"));
		cl.def("run", (void (ROL::TypeB::NewtonKrylovAlgorithm<double>::*)(class ROL::Problem<double> &, std::ostream &)) &ROL::TypeB::NewtonKrylovAlgorithm<double>::run, "C++: ROL::TypeB::NewtonKrylovAlgorithm<double>::run(class ROL::Problem<double> &, std::ostream &) --> void", pybind11::arg("problem"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::NewtonKrylovAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", (void (ROL::TypeB::NewtonKrylovAlgorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeB::NewtonKrylovAlgorithm<double>::run, "C++: ROL::TypeB::NewtonKrylovAlgorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("writeHeader", (void (ROL::TypeB::NewtonKrylovAlgorithm<double>::*)(std::ostream &) const) &ROL::TypeB::NewtonKrylovAlgorithm<double>::writeHeader, "C++: ROL::TypeB::NewtonKrylovAlgorithm<double>::writeHeader(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("writeName", (void (ROL::TypeB::NewtonKrylovAlgorithm<double>::*)(std::ostream &) const) &ROL::TypeB::NewtonKrylovAlgorithm<double>::writeName, "C++: ROL::TypeB::NewtonKrylovAlgorithm<double>::writeName(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("writeOutput", [](ROL::TypeB::NewtonKrylovAlgorithm<double> const &o, std::ostream & a0) -> void { return o.writeOutput(a0); }, "", pybind11::arg("os"));
		cl.def("writeOutput", (void (ROL::TypeB::NewtonKrylovAlgorithm<double>::*)(std::ostream &, const bool) const) &ROL::TypeB::NewtonKrylovAlgorithm<double>::writeOutput, "C++: ROL::TypeB::NewtonKrylovAlgorithm<double>::writeOutput(std::ostream &, const bool) const --> void", pybind11::arg("os"), pybind11::arg("write_header"));
		cl.def("setStatusTest", [](ROL::TypeB::Algorithm<double> &o, const class Teuchos::RCP<class ROL::StatusTest<double> > & a0) -> void { return o.setStatusTest(a0); }, "", pybind11::arg("status"));
		cl.def("setStatusTest", (void (ROL::TypeB::Algorithm<double>::*)(const class Teuchos::RCP<class ROL::StatusTest<double> > &, const bool)) &ROL::TypeB::Algorithm<double>::setStatusTest, "C++: ROL::TypeB::Algorithm<double>::setStatusTest(const class Teuchos::RCP<class ROL::StatusTest<double> > &, const bool) --> void", pybind11::arg("status"), pybind11::arg("combineStatus"));
		cl.def("run", [](ROL::TypeB::Algorithm<double> &o, class ROL::Problem<double> & a0) -> void { return o.run(a0); }, "", pybind11::arg("problem"));
		cl.def("run", (void (ROL::TypeB::Algorithm<double>::*)(class ROL::Problem<double> &, std::ostream &)) &ROL::TypeB::Algorithm<double>::run, "C++: ROL::TypeB::Algorithm<double>::run(class ROL::Problem<double> &, std::ostream &) --> void", pybind11::arg("problem"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2) -> void { return o.run(a0, a1, a2); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"));
		cl.def("run", (void (ROL::TypeB::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeB::Algorithm<double>::run, "C++: ROL::TypeB::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3) -> void { return o.run(a0, a1, a2, a3); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"));
		cl.def("run", (void (ROL::TypeB::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeB::Algorithm<double>::run, "C++: ROL::TypeB::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4) -> void { return o.run(a0, a1, a2, a3, a4); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", (void (ROL::TypeB::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &)) &ROL::TypeB::Algorithm<double>::run, "C++: ROL::TypeB::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", (void (ROL::TypeB::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeB::Algorithm<double>::run, "C++: ROL::TypeB::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5) -> void { return o.run(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"));
		cl.def("run", (void (ROL::TypeB::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeB::Algorithm<double>::run, "C++: ROL::TypeB::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("linear_ires"));
		cl.def("run", (void (ROL::TypeB::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeB::Algorithm<double>::run, "C++: ROL::TypeB::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("linear_ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::Algorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"));
		cl.def("run", (void (ROL::TypeB::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeB::Algorithm<double>::run, "C++: ROL::TypeB::Algorithm<double>::run(class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::Algorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("linear_ires"));
		cl.def("run", (void (ROL::TypeB::Algorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &)) &ROL::TypeB::Algorithm<double>::run, "C++: ROL::TypeB::Algorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Constraint<double> &, class ROL::Vector<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("linear_ires"), pybind11::arg("outStream"));
		cl.def("writeHeader", (void (ROL::TypeB::Algorithm<double>::*)(std::ostream &) const) &ROL::TypeB::Algorithm<double>::writeHeader, "C++: ROL::TypeB::Algorithm<double>::writeHeader(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("writeName", (void (ROL::TypeB::Algorithm<double>::*)(std::ostream &) const) &ROL::TypeB::Algorithm<double>::writeName, "C++: ROL::TypeB::Algorithm<double>::writeName(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("writeOutput", [](ROL::TypeB::Algorithm<double> const &o, std::ostream & a0) -> void { return o.writeOutput(a0); }, "", pybind11::arg("os"));
		cl.def("writeOutput", (void (ROL::TypeB::Algorithm<double>::*)(std::ostream &, const bool) const) &ROL::TypeB::Algorithm<double>::writeOutput, "C++: ROL::TypeB::Algorithm<double>::writeOutput(std::ostream &, const bool) const --> void", pybind11::arg("os"), pybind11::arg("write_header"));
		cl.def("writeExitStatus", (void (ROL::TypeB::Algorithm<double>::*)(std::ostream &) const) &ROL::TypeB::Algorithm<double>::writeExitStatus, "C++: ROL::TypeB::Algorithm<double>::writeExitStatus(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("getState", (class Teuchos::RCP<const struct ROL::TypeB::AlgorithmState<double> > (ROL::TypeB::Algorithm<double>::*)() const) &ROL::TypeB::Algorithm<double>::getState, "C++: ROL::TypeB::Algorithm<double>::getState() const --> class Teuchos::RCP<const struct ROL::TypeB::AlgorithmState<double> >");
		cl.def("reset", (void (ROL::TypeB::Algorithm<double>::*)()) &ROL::TypeB::Algorithm<double>::reset, "C++: ROL::TypeB::Algorithm<double>::reset() --> void");
	}
}
