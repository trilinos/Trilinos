#include <ROL_BoundConstraint.hpp>
#include <ROL_Constraint.hpp>
#include <ROL_Elementwise_Function.hpp>
#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_Objective.hpp>
#include <ROL_PolyhedralProjection.hpp>
#include <ROL_Problem.hpp>
#include <ROL_Secant.hpp>
#include <ROL_StatusTest.hpp>
#include <ROL_TypeB_Algorithm.hpp>
#include <ROL_TypeB_AlgorithmFactory.hpp>
#include <ROL_TypeB_ColemanLiAlgorithm.hpp>
#include <ROL_TypeB_TrustRegionSPGAlgorithm.hpp>
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

// ROL::TypeB::TrustRegionSPGAlgorithm file:ROL_TypeB_TrustRegionSPGAlgorithm.hpp line:59
struct PyCallBack_ROL_TypeB_TrustRegionSPGAlgorithm_double_t : public ROL::TypeB::TrustRegionSPGAlgorithm<double> {
	using ROL::TypeB::TrustRegionSPGAlgorithm<double>::TrustRegionSPGAlgorithm;

	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, std::ostream & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::TrustRegionSPGAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TrustRegionSPGAlgorithm::run(a0, a1, a2, a3, a4);
	}
	void writeHeader(std::ostream & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::TrustRegionSPGAlgorithm<double> *>(this), "writeHeader");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TrustRegionSPGAlgorithm::writeHeader(a0);
	}
	void writeName(std::ostream & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::TrustRegionSPGAlgorithm<double> *>(this), "writeName");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TrustRegionSPGAlgorithm::writeName(a0);
	}
	void writeOutput(std::ostream & a0, const bool a1) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::TrustRegionSPGAlgorithm<double> *>(this), "writeOutput");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TrustRegionSPGAlgorithm::writeOutput(a0, a1);
	}
	void run(class ROL::Problem<double> & a0, std::ostream & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::TrustRegionSPGAlgorithm<double> *>(this), "run");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::TrustRegionSPGAlgorithm<double> *>(this), "run");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::TrustRegionSPGAlgorithm<double> *>(this), "run");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::TrustRegionSPGAlgorithm<double> *>(this), "run");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::TrustRegionSPGAlgorithm<double> *>(this), "run");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::TrustRegionSPGAlgorithm<double> *>(this), "run");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::TrustRegionSPGAlgorithm<double> *>(this), "run");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::TrustRegionSPGAlgorithm<double> *>(this), "run");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::TrustRegionSPGAlgorithm<double> *>(this), "writeExitStatus");
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

// ROL::TypeB::ColemanLiAlgorithm file:ROL_TypeB_ColemanLiAlgorithm.hpp line:61
struct PyCallBack_ROL_TypeB_ColemanLiAlgorithm_double_t : public ROL::TypeB::ColemanLiAlgorithm<double> {
	using ROL::TypeB::ColemanLiAlgorithm<double>::ColemanLiAlgorithm;

	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, std::ostream & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::ColemanLiAlgorithm<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ColemanLiAlgorithm::run(a0, a1, a2, a3, a4);
	}
	void writeHeader(std::ostream & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::ColemanLiAlgorithm<double> *>(this), "writeHeader");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ColemanLiAlgorithm::writeHeader(a0);
	}
	void writeName(std::ostream & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::ColemanLiAlgorithm<double> *>(this), "writeName");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ColemanLiAlgorithm::writeName(a0);
	}
	void writeOutput(std::ostream & a0, const bool a1) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::ColemanLiAlgorithm<double> *>(this), "writeOutput");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ColemanLiAlgorithm::writeOutput(a0, a1);
	}
	void run(class ROL::Problem<double> & a0, std::ostream & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::ColemanLiAlgorithm<double> *>(this), "run");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::ColemanLiAlgorithm<double> *>(this), "run");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::ColemanLiAlgorithm<double> *>(this), "run");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::ColemanLiAlgorithm<double> *>(this), "run");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::ColemanLiAlgorithm<double> *>(this), "run");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::ColemanLiAlgorithm<double> *>(this), "run");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::ColemanLiAlgorithm<double> *>(this), "run");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::ColemanLiAlgorithm<double> *>(this), "run");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TypeB::ColemanLiAlgorithm<double> *>(this), "writeExitStatus");
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

void bind_ROL_TypeB_TrustRegionSPGAlgorithm(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // ROL::TypeB::TrustRegionSPGAlgorithm file:ROL_TypeB_TrustRegionSPGAlgorithm.hpp line:59
		pybind11::class_<ROL::TypeB::TrustRegionSPGAlgorithm<double>, Teuchos::RCP<ROL::TypeB::TrustRegionSPGAlgorithm<double>>, PyCallBack_ROL_TypeB_TrustRegionSPGAlgorithm_double_t, ROL::TypeB::Algorithm<double>> cl(M("ROL::TypeB"), "TrustRegionSPGAlgorithm_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](class Teuchos::ParameterList & a0){ return new ROL::TypeB::TrustRegionSPGAlgorithm<double>(a0); }, [](class Teuchos::ParameterList & a0){ return new PyCallBack_ROL_TypeB_TrustRegionSPGAlgorithm_double_t(a0); } ), "doc");
		cl.def( pybind11::init<class Teuchos::ParameterList &, const class Teuchos::RCP<class ROL::Secant<double> > &>(), pybind11::arg("list"), pybind11::arg("secant") );

		cl.def( pybind11::init( [](PyCallBack_ROL_TypeB_TrustRegionSPGAlgorithm_double_t const &o){ return new PyCallBack_ROL_TypeB_TrustRegionSPGAlgorithm_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::TypeB::TrustRegionSPGAlgorithm<double> const &o){ return new ROL::TypeB::TrustRegionSPGAlgorithm<double>(o); } ) );
		cl.def("run", [](ROL::TypeB::TrustRegionSPGAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("linear_ires"));
		cl.def("run", [](ROL::TypeB::TrustRegionSPGAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10, std::ostream & a11) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("linear_ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::TrustRegionSPGAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"));
		cl.def("run", [](ROL::TypeB::TrustRegionSPGAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7, std::ostream & a8) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::TrustRegionSPGAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("linear_ires"));
		cl.def("run", [](ROL::TypeB::TrustRegionSPGAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7, std::ostream & a8) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("linear_ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::TrustRegionSPGAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5) -> void { return o.run(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"));
		cl.def("run", [](ROL::TypeB::TrustRegionSPGAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, std::ostream & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::TrustRegionSPGAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", [](ROL::TypeB::TrustRegionSPGAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, std::ostream & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::TrustRegionSPGAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4) -> void { return o.run(a0, a1, a2, a3, a4); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", [](ROL::TypeB::TrustRegionSPGAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, std::ostream & a5) -> void { return o.run(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::TrustRegionSPGAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2) -> void { return o.run(a0, a1, a2); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"));
		cl.def("run", [](ROL::TypeB::TrustRegionSPGAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, std::ostream & a3) -> void { return o.run(a0, a1, a2, a3); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::TrustRegionSPGAlgorithm<double> &o, class ROL::Problem<double> & a0) -> void { return o.run(a0); }, "", pybind11::arg("problem"));
		cl.def("run", [](ROL::TypeB::TrustRegionSPGAlgorithm<double> &o, class ROL::Problem<double> & a0, std::ostream & a1) -> void { return o.run(a0, a1); }, "", pybind11::arg("problem"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::TrustRegionSPGAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3) -> void { return o.run(a0, a1, a2, a3); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"));
		cl.def("run", (void (ROL::TypeB::TrustRegionSPGAlgorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeB::TrustRegionSPGAlgorithm<double>::run, "C++: ROL::TypeB::TrustRegionSPGAlgorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("outStream"));
		cl.def("writeHeader", (void (ROL::TypeB::TrustRegionSPGAlgorithm<double>::*)(std::ostream &) const) &ROL::TypeB::TrustRegionSPGAlgorithm<double>::writeHeader, "C++: ROL::TypeB::TrustRegionSPGAlgorithm<double>::writeHeader(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("writeName", (void (ROL::TypeB::TrustRegionSPGAlgorithm<double>::*)(std::ostream &) const) &ROL::TypeB::TrustRegionSPGAlgorithm<double>::writeName, "C++: ROL::TypeB::TrustRegionSPGAlgorithm<double>::writeName(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("writeOutput", [](ROL::TypeB::TrustRegionSPGAlgorithm<double> const &o, std::ostream & a0) -> void { return o.writeOutput(a0); }, "", pybind11::arg("os"));
		cl.def("writeOutput", (void (ROL::TypeB::TrustRegionSPGAlgorithm<double>::*)(std::ostream &, const bool) const) &ROL::TypeB::TrustRegionSPGAlgorithm<double>::writeOutput, "C++: ROL::TypeB::TrustRegionSPGAlgorithm<double>::writeOutput(std::ostream &, const bool) const --> void", pybind11::arg("os"), pybind11::arg("write_header"));
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
	{ // ROL::TypeB::ColemanLiAlgorithm file:ROL_TypeB_ColemanLiAlgorithm.hpp line:61
		pybind11::class_<ROL::TypeB::ColemanLiAlgorithm<double>, Teuchos::RCP<ROL::TypeB::ColemanLiAlgorithm<double>>, PyCallBack_ROL_TypeB_ColemanLiAlgorithm_double_t, ROL::TypeB::Algorithm<double>> cl(M("ROL::TypeB"), "ColemanLiAlgorithm_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](class Teuchos::ParameterList & a0){ return new ROL::TypeB::ColemanLiAlgorithm<double>(a0); }, [](class Teuchos::ParameterList & a0){ return new PyCallBack_ROL_TypeB_ColemanLiAlgorithm_double_t(a0); } ), "doc");
		cl.def( pybind11::init<class Teuchos::ParameterList &, const class Teuchos::RCP<class ROL::Secant<double> > &>(), pybind11::arg("list"), pybind11::arg("secant") );

		cl.def( pybind11::init( [](PyCallBack_ROL_TypeB_ColemanLiAlgorithm_double_t const &o){ return new PyCallBack_ROL_TypeB_ColemanLiAlgorithm_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::TypeB::ColemanLiAlgorithm<double> const &o){ return new ROL::TypeB::ColemanLiAlgorithm<double>(o); } ) );
		cl.def("run", [](ROL::TypeB::ColemanLiAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("linear_ires"));
		cl.def("run", [](ROL::TypeB::ColemanLiAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Constraint<double> & a7, class ROL::Vector<double> & a8, class ROL::BoundConstraint<double> & a9, const class ROL::Vector<double> & a10, std::ostream & a11) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("linear_ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::ColemanLiAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"));
		cl.def("run", [](ROL::TypeB::ColemanLiAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::Constraint<double> & a5, class ROL::Vector<double> & a6, class ROL::BoundConstraint<double> & a7, std::ostream & a8) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::ColemanLiAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("linear_ires"));
		cl.def("run", [](ROL::TypeB::ColemanLiAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, class ROL::BoundConstraint<double> & a6, const class ROL::Vector<double> & a7, std::ostream & a8) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("linear_ires"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::ColemanLiAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5) -> void { return o.run(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"));
		cl.def("run", [](ROL::TypeB::ColemanLiAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, class ROL::BoundConstraint<double> & a5, std::ostream & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::ColemanLiAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("run", [](ROL::TypeB::ColemanLiAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, class ROL::Constraint<double> & a4, class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, std::ostream & a7) -> void { return o.run(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::ColemanLiAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4) -> void { return o.run(a0, a1, a2, a3, a4); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("run", [](ROL::TypeB::ColemanLiAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, class ROL::Constraint<double> & a3, class ROL::Vector<double> & a4, std::ostream & a5) -> void { return o.run(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::ColemanLiAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2) -> void { return o.run(a0, a1, a2); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"));
		cl.def("run", [](ROL::TypeB::ColemanLiAlgorithm<double> &o, class ROL::Vector<double> & a0, class ROL::Objective<double> & a1, class ROL::BoundConstraint<double> & a2, std::ostream & a3) -> void { return o.run(a0, a1, a2, a3); }, "", pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::ColemanLiAlgorithm<double> &o, class ROL::Problem<double> & a0) -> void { return o.run(a0); }, "", pybind11::arg("problem"));
		cl.def("run", [](ROL::TypeB::ColemanLiAlgorithm<double> &o, class ROL::Problem<double> & a0, std::ostream & a1) -> void { return o.run(a0, a1); }, "", pybind11::arg("problem"), pybind11::arg("outStream"));
		cl.def("run", [](ROL::TypeB::ColemanLiAlgorithm<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3) -> void { return o.run(a0, a1, a2, a3); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"));
		cl.def("run", (void (ROL::TypeB::ColemanLiAlgorithm<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, std::ostream &)) &ROL::TypeB::ColemanLiAlgorithm<double>::run, "C++: ROL::TypeB::ColemanLiAlgorithm<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, std::ostream &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("outStream"));
		cl.def("writeHeader", (void (ROL::TypeB::ColemanLiAlgorithm<double>::*)(std::ostream &) const) &ROL::TypeB::ColemanLiAlgorithm<double>::writeHeader, "C++: ROL::TypeB::ColemanLiAlgorithm<double>::writeHeader(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("writeName", (void (ROL::TypeB::ColemanLiAlgorithm<double>::*)(std::ostream &) const) &ROL::TypeB::ColemanLiAlgorithm<double>::writeName, "C++: ROL::TypeB::ColemanLiAlgorithm<double>::writeName(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("writeOutput", [](ROL::TypeB::ColemanLiAlgorithm<double> const &o, std::ostream & a0) -> void { return o.writeOutput(a0); }, "", pybind11::arg("os"));
		cl.def("writeOutput", (void (ROL::TypeB::ColemanLiAlgorithm<double>::*)(std::ostream &, const bool) const) &ROL::TypeB::ColemanLiAlgorithm<double>::writeOutput, "C++: ROL::TypeB::ColemanLiAlgorithm<double>::writeOutput(std::ostream &, const bool) const --> void", pybind11::arg("os"), pybind11::arg("write_header"));
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
	// ROL::TypeB::EAlgorithmB file:ROL_TypeB_AlgorithmFactory.hpp line:74
	pybind11::enum_<ROL::TypeB::EAlgorithmB>(M("ROL::TypeB"), "EAlgorithmB", pybind11::arithmetic(), "Enumeration of bound constrained algorithm types.\n\n    \n    ALGORITHM_B_LINESEARCH          describe\n    \n\n    ALGORITHM_B_TRUSTREGION         describe\n    \n\n    ALGORITHM_B_MOREAUYOSIDA        describe\n    \n\n    ALGORITHM_B_PRIMALDUALACTIVESET describe\n    \n\n    ALGORITHM_B_INTERIORPOINT       describe\n    \n\n    ALGORITHM_B_SPECTRALGRADIENT    describe", pybind11::module_local())
		.value("ALGORITHM_B_LINESEARCH", ROL::TypeB::ALGORITHM_B_LINESEARCH)
		.value("ALGORITHM_B_TRUSTREGION", ROL::TypeB::ALGORITHM_B_TRUSTREGION)
		.value("ALGORITHM_B_MOREAUYOSIDA", ROL::TypeB::ALGORITHM_B_MOREAUYOSIDA)
		.value("ALGORITHM_B_PRIMALDUALACTIVESET", ROL::TypeB::ALGORITHM_B_PRIMALDUALACTIVESET)
		.value("ALGORITHM_B_INTERIORPOINT", ROL::TypeB::ALGORITHM_B_INTERIORPOINT)
		.value("ALGORITHM_B_SPECTRALGRADIENT", ROL::TypeB::ALGORITHM_B_SPECTRALGRADIENT)
		.value("ALGORITHM_B_LAST", ROL::TypeB::ALGORITHM_B_LAST)
		.export_values();

;

	// ROL::TypeB::EAlgorithmBToString(enum ROL::TypeB::EAlgorithmB) file:ROL_TypeB_AlgorithmFactory.hpp line:84
	M("ROL::TypeB").def("EAlgorithmBToString", (std::string (*)(enum ROL::TypeB::EAlgorithmB)) &ROL::TypeB::EAlgorithmBToString, "C++: ROL::TypeB::EAlgorithmBToString(enum ROL::TypeB::EAlgorithmB) --> std::string", pybind11::arg("alg"));

	// ROL::TypeB::isValidAlgorithmB(enum ROL::TypeB::EAlgorithmB) file:ROL_TypeB_AlgorithmFactory.hpp line:104
	M("ROL::TypeB").def("isValidAlgorithmB", (int (*)(enum ROL::TypeB::EAlgorithmB)) &ROL::TypeB::isValidAlgorithmB, "Verifies validity of a AlgorithmB enum.\n\n    \n  [in]  - enum of the AlgorithmB\n    \n\n 1 if the argument is a valid AlgorithmB; 0 otherwise.\n\nC++: ROL::TypeB::isValidAlgorithmB(enum ROL::TypeB::EAlgorithmB) --> int", pybind11::arg("alg"));

	// ROL::TypeB::StringToEAlgorithmB(std::string) file:ROL_TypeB_AlgorithmFactory.hpp line:135
	M("ROL::TypeB").def("StringToEAlgorithmB", (enum ROL::TypeB::EAlgorithmB (*)(std::string)) &ROL::TypeB::StringToEAlgorithmB, "C++: ROL::TypeB::StringToEAlgorithmB(std::string) --> enum ROL::TypeB::EAlgorithmB", pybind11::arg("s"));

	// ROL::TypeB::AlgorithmFactory(class Teuchos::ParameterList &) file:ROL_TypeB_AlgorithmFactory.hpp line:146
	M("ROL::TypeB").def("AlgorithmFactory", (class Teuchos::RCP<class ROL::TypeB::Algorithm<double> > (*)(class Teuchos::ParameterList &)) &ROL::TypeB::AlgorithmFactory<double>, "C++: ROL::TypeB::AlgorithmFactory(class Teuchos::ParameterList &) --> class Teuchos::RCP<class ROL::TypeB::Algorithm<double> >", pybind11::arg("parlist"));

}
