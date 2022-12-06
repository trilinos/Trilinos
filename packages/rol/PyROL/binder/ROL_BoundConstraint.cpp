#include <ROL_BoundConstraint.hpp>
#include <ROL_BoundConstraint_Partitioned.hpp>
#include <ROL_Bounds.hpp>
#include <ROL_Constraint.hpp>
#include <ROL_ConstraintAssembler.hpp>
#include <ROL_Elementwise_Function.hpp>
#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_LinearConstraint.hpp>
#include <ROL_LinearOperator.hpp>
#include <ROL_SlacklessObjective.hpp>
#include <ROL_UpdateType.hpp>
#include <ROL_Vector.hpp>
#include <ROL_VectorController.hpp>
#include <Teuchos_ENull.hpp>
#include <Teuchos_RCPDecl.hpp>
#include <Teuchos_RCPNode.hpp>
#include <Teuchos_any.hpp>
#include <cwchar>
#include <functional>
#include <ios>
#include <iterator>
#include <memory>
#include <ostream>
#include <sstream> // __str__
#include <streambuf>
#include <string>
#include <unordered_map>
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

// ROL::BoundConstraint file:ROL_BoundConstraint.hpp line:73
struct PyCallBack_ROL_BoundConstraint_double_t : public ROL::BoundConstraint<double> {
	using ROL::BoundConstraint<double>::BoundConstraint;

	void project(class ROL::Vector<double> & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BoundConstraint<double> *>(this), "project");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BoundConstraint::project(a0);
	}
	void projectInterior(class ROL::Vector<double> & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BoundConstraint<double> *>(this), "projectInterior");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BoundConstraint::projectInterior(a0);
	}
	void pruneUpperActive(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BoundConstraint<double> *>(this), "pruneUpperActive");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BoundConstraint::pruneUpperActive(a0, a1, a2);
	}
	void pruneUpperActive(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double a3, double a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BoundConstraint<double> *>(this), "pruneUpperActive");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BoundConstraint::pruneUpperActive(a0, a1, a2, a3, a4);
	}
	void pruneLowerActive(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BoundConstraint<double> *>(this), "pruneLowerActive");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BoundConstraint::pruneLowerActive(a0, a1, a2);
	}
	void pruneLowerActive(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double a3, double a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BoundConstraint<double> *>(this), "pruneLowerActive");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BoundConstraint::pruneLowerActive(a0, a1, a2, a3, a4);
	}
	const class Teuchos::RCP<const class ROL::Vector<double> > getLowerBound() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BoundConstraint<double> *>(this), "getLowerBound");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const class Teuchos::RCP<const class ROL::Vector<double> >>::value) {
				static pybind11::detail::override_caster_t<const class Teuchos::RCP<const class ROL::Vector<double> >> caster;
				return pybind11::detail::cast_ref<const class Teuchos::RCP<const class ROL::Vector<double> >>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const class Teuchos::RCP<const class ROL::Vector<double> >>(std::move(o));
		}
		return BoundConstraint::getLowerBound();
	}
	const class Teuchos::RCP<const class ROL::Vector<double> > getUpperBound() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BoundConstraint<double> *>(this), "getUpperBound");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const class Teuchos::RCP<const class ROL::Vector<double> >>::value) {
				static pybind11::detail::override_caster_t<const class Teuchos::RCP<const class ROL::Vector<double> >> caster;
				return pybind11::detail::cast_ref<const class Teuchos::RCP<const class ROL::Vector<double> >>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const class Teuchos::RCP<const class ROL::Vector<double> >>(std::move(o));
		}
		return BoundConstraint::getUpperBound();
	}
	bool isFeasible(const class ROL::Vector<double> & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BoundConstraint<double> *>(this), "isFeasible");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return BoundConstraint::isFeasible(a0);
	}
	void applyInverseScalingFunction(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BoundConstraint<double> *>(this), "applyInverseScalingFunction");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BoundConstraint::applyInverseScalingFunction(a0, a1, a2, a3);
	}
	void applyScalingFunctionJacobian(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BoundConstraint<double> *>(this), "applyScalingFunctionJacobian");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BoundConstraint::applyScalingFunctionJacobian(a0, a1, a2, a3);
	}
};

// ROL::BoundConstraint_Partitioned file:ROL_BoundConstraint_Partitioned.hpp line:61
struct PyCallBack_ROL_BoundConstraint_Partitioned_double_t : public ROL::BoundConstraint_Partitioned<double> {
	using ROL::BoundConstraint_Partitioned<double>::BoundConstraint_Partitioned;

	void project(class ROL::Vector<double> & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BoundConstraint_Partitioned<double> *>(this), "project");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BoundConstraint_Partitioned::project(a0);
	}
	void projectInterior(class ROL::Vector<double> & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BoundConstraint_Partitioned<double> *>(this), "projectInterior");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BoundConstraint_Partitioned::projectInterior(a0);
	}
	void pruneUpperActive(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BoundConstraint_Partitioned<double> *>(this), "pruneUpperActive");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BoundConstraint_Partitioned::pruneUpperActive(a0, a1, a2);
	}
	void pruneUpperActive(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double a3, double a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BoundConstraint_Partitioned<double> *>(this), "pruneUpperActive");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BoundConstraint_Partitioned::pruneUpperActive(a0, a1, a2, a3, a4);
	}
	void pruneLowerActive(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BoundConstraint_Partitioned<double> *>(this), "pruneLowerActive");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BoundConstraint_Partitioned::pruneLowerActive(a0, a1, a2);
	}
	void pruneLowerActive(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double a3, double a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BoundConstraint_Partitioned<double> *>(this), "pruneLowerActive");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BoundConstraint_Partitioned::pruneLowerActive(a0, a1, a2, a3, a4);
	}
	bool isFeasible(const class ROL::Vector<double> & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BoundConstraint_Partitioned<double> *>(this), "isFeasible");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return BoundConstraint_Partitioned::isFeasible(a0);
	}
	void applyInverseScalingFunction(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BoundConstraint_Partitioned<double> *>(this), "applyInverseScalingFunction");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BoundConstraint_Partitioned::applyInverseScalingFunction(a0, a1, a2, a3);
	}
	void applyScalingFunctionJacobian(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BoundConstraint_Partitioned<double> *>(this), "applyScalingFunctionJacobian");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BoundConstraint_Partitioned::applyScalingFunctionJacobian(a0, a1, a2, a3);
	}
	const class Teuchos::RCP<const class ROL::Vector<double> > getLowerBound() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BoundConstraint_Partitioned<double> *>(this), "getLowerBound");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const class Teuchos::RCP<const class ROL::Vector<double> >>::value) {
				static pybind11::detail::override_caster_t<const class Teuchos::RCP<const class ROL::Vector<double> >> caster;
				return pybind11::detail::cast_ref<const class Teuchos::RCP<const class ROL::Vector<double> >>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const class Teuchos::RCP<const class ROL::Vector<double> >>(std::move(o));
		}
		return BoundConstraint::getLowerBound();
	}
	const class Teuchos::RCP<const class ROL::Vector<double> > getUpperBound() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BoundConstraint_Partitioned<double> *>(this), "getUpperBound");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const class Teuchos::RCP<const class ROL::Vector<double> >>::value) {
				static pybind11::detail::override_caster_t<const class Teuchos::RCP<const class ROL::Vector<double> >> caster;
				return pybind11::detail::cast_ref<const class Teuchos::RCP<const class ROL::Vector<double> >>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const class Teuchos::RCP<const class ROL::Vector<double> >>(std::move(o));
		}
		return BoundConstraint::getUpperBound();
	}
};

// ROL::SlacklessObjective file:ROL_SlacklessObjective.hpp line:59
struct PyCallBack_ROL_SlacklessObjective_double_t : public ROL::SlacklessObjective<double> {
	using ROL::SlacklessObjective<double>::SlacklessObjective;

	void update(const class ROL::Vector<double> & a0, enum ROL::UpdateType a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::SlacklessObjective<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SlacklessObjective::update(a0, a1, a2);
	}
	void update(const class ROL::Vector<double> & a0, bool a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::SlacklessObjective<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SlacklessObjective::update(a0, a1, a2);
	}
	double value(const class ROL::Vector<double> & a0, double & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::SlacklessObjective<double> *>(this), "value");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SlacklessObjective::value(a0, a1);
	}
	double dirDeriv(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double & a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::SlacklessObjective<double> *>(this), "dirDeriv");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SlacklessObjective::dirDeriv(a0, a1, a2);
	}
	void gradient(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double & a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::SlacklessObjective<double> *>(this), "gradient");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SlacklessObjective::gradient(a0, a1, a2);
	}
	void hessVec(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::SlacklessObjective<double> *>(this), "hessVec");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SlacklessObjective::hessVec(a0, a1, a2, a3);
	}
	void invHessVec(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::SlacklessObjective<double> *>(this), "invHessVec");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SlacklessObjective::invHessVec(a0, a1, a2, a3);
	}
	void precond(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::SlacklessObjective<double> *>(this), "precond");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SlacklessObjective::precond(a0, a1, a2, a3);
	}
};

// ROL::LinearConstraint file:ROL_LinearConstraint.hpp line:60
struct PyCallBack_ROL_LinearConstraint_double_t : public ROL::LinearConstraint<double> {
	using ROL::LinearConstraint<double>::LinearConstraint;

	void update(const class ROL::Vector<double> & a0, enum ROL::UpdateType a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LinearConstraint<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return LinearConstraint::update(a0, a1, a2);
	}
	void update(const class ROL::Vector<double> & a0, bool a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LinearConstraint<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return LinearConstraint::update(a0, a1, a2);
	}
	void value(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double & a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LinearConstraint<double> *>(this), "value");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return LinearConstraint::value(a0, a1, a2);
	}
	void applyJacobian(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LinearConstraint<double> *>(this), "applyJacobian");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return LinearConstraint::applyJacobian(a0, a1, a2, a3);
	}
	void applyAdjointJacobian(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LinearConstraint<double> *>(this), "applyAdjointJacobian");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return LinearConstraint::applyAdjointJacobian(a0, a1, a2, a3);
	}
	void applyAdjointJacobian(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, double & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LinearConstraint<double> *>(this), "applyAdjointJacobian");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return LinearConstraint::applyAdjointJacobian(a0, a1, a2, a3, a4);
	}
	void applyAdjointHessian(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, double & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LinearConstraint<double> *>(this), "applyAdjointHessian");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return LinearConstraint::applyAdjointHessian(a0, a1, a2, a3, a4);
	}
	void applyPreconditioner(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, double & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LinearConstraint<double> *>(this), "applyPreconditioner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint::applyPreconditioner(a0, a1, a2, a3, a4);
	}
};

// ROL::Bounds file:ROL_Bounds.hpp line:59
struct PyCallBack_ROL_Bounds_double_t : public ROL::Bounds<double> {
	using ROL::Bounds<double>::Bounds;

	void project(class ROL::Vector<double> & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Bounds<double> *>(this), "project");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Bounds::project(a0);
	}
	void projectInterior(class ROL::Vector<double> & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Bounds<double> *>(this), "projectInterior");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Bounds::projectInterior(a0);
	}
	void pruneUpperActive(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Bounds<double> *>(this), "pruneUpperActive");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Bounds::pruneUpperActive(a0, a1, a2);
	}
	void pruneUpperActive(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double a3, double a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Bounds<double> *>(this), "pruneUpperActive");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Bounds::pruneUpperActive(a0, a1, a2, a3, a4);
	}
	void pruneLowerActive(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Bounds<double> *>(this), "pruneLowerActive");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Bounds::pruneLowerActive(a0, a1, a2);
	}
	void pruneLowerActive(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double a3, double a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Bounds<double> *>(this), "pruneLowerActive");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Bounds::pruneLowerActive(a0, a1, a2, a3, a4);
	}
	bool isFeasible(const class ROL::Vector<double> & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Bounds<double> *>(this), "isFeasible");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return Bounds::isFeasible(a0);
	}
	void applyInverseScalingFunction(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Bounds<double> *>(this), "applyInverseScalingFunction");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Bounds::applyInverseScalingFunction(a0, a1, a2, a3);
	}
	void applyScalingFunctionJacobian(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Bounds<double> *>(this), "applyScalingFunctionJacobian");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Bounds::applyScalingFunctionJacobian(a0, a1, a2, a3);
	}
	const class Teuchos::RCP<const class ROL::Vector<double> > getLowerBound() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Bounds<double> *>(this), "getLowerBound");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const class Teuchos::RCP<const class ROL::Vector<double> >>::value) {
				static pybind11::detail::override_caster_t<const class Teuchos::RCP<const class ROL::Vector<double> >> caster;
				return pybind11::detail::cast_ref<const class Teuchos::RCP<const class ROL::Vector<double> >>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const class Teuchos::RCP<const class ROL::Vector<double> >>(std::move(o));
		}
		return BoundConstraint::getLowerBound();
	}
	const class Teuchos::RCP<const class ROL::Vector<double> > getUpperBound() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Bounds<double> *>(this), "getUpperBound");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const class Teuchos::RCP<const class ROL::Vector<double> >>::value) {
				static pybind11::detail::override_caster_t<const class Teuchos::RCP<const class ROL::Vector<double> >> caster;
				return pybind11::detail::cast_ref<const class Teuchos::RCP<const class ROL::Vector<double> >>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const class Teuchos::RCP<const class ROL::Vector<double> >>(std::move(o));
		}
		return BoundConstraint::getUpperBound();
	}
};

void bind_ROL_BoundConstraint(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // ROL::BoundConstraint file:ROL_BoundConstraint.hpp line:73
		pybind11::class_<ROL::BoundConstraint<double>, Teuchos::RCP<ROL::BoundConstraint<double>>, PyCallBack_ROL_BoundConstraint_double_t> cl(M("ROL"), "BoundConstraint_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new ROL::BoundConstraint<double>(); }, [](){ return new PyCallBack_ROL_BoundConstraint_double_t(); } ) );
		cl.def( pybind11::init<const class ROL::Vector<double> &>(), pybind11::arg("x") );

		cl.def( pybind11::init( [](PyCallBack_ROL_BoundConstraint_double_t const &o){ return new PyCallBack_ROL_BoundConstraint_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::BoundConstraint<double> const &o){ return new ROL::BoundConstraint<double>(o); } ) );
		cl.def("project", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &)) &ROL::BoundConstraint<double>::project, "C++: ROL::BoundConstraint<double>::project(class ROL::Vector<double> &) --> void", pybind11::arg("x"));
		cl.def("projectInterior", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &)) &ROL::BoundConstraint<double>::projectInterior, "C++: ROL::BoundConstraint<double>::projectInterior(class ROL::Vector<double> &) --> void", pybind11::arg("x"));
		cl.def("pruneUpperActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) -> void { return o.pruneUpperActive(a0, a1); }, "", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("pruneUpperActive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double)) &ROL::BoundConstraint<double>::pruneUpperActive, "C++: ROL::BoundConstraint<double>::pruneUpperActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, double) --> void", pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("eps"));
		cl.def("pruneUpperActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> void { return o.pruneUpperActive(a0, a1, a2); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"));
		cl.def("pruneUpperActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double const & a3) -> void { return o.pruneUpperActive(a0, a1, a2, a3); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"));
		cl.def("pruneUpperActive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double)) &ROL::BoundConstraint<double>::pruneUpperActive, "C++: ROL::BoundConstraint<double>::pruneUpperActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double) --> void", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"), pybind11::arg("geps"));
		cl.def("pruneLowerActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) -> void { return o.pruneLowerActive(a0, a1); }, "", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("pruneLowerActive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double)) &ROL::BoundConstraint<double>::pruneLowerActive, "C++: ROL::BoundConstraint<double>::pruneLowerActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, double) --> void", pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("eps"));
		cl.def("pruneLowerActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> void { return o.pruneLowerActive(a0, a1, a2); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"));
		cl.def("pruneLowerActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double const & a3) -> void { return o.pruneLowerActive(a0, a1, a2, a3); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"));
		cl.def("pruneLowerActive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double)) &ROL::BoundConstraint<double>::pruneLowerActive, "C++: ROL::BoundConstraint<double>::pruneLowerActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double) --> void", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"), pybind11::arg("geps"));
		cl.def("getLowerBound", (const class Teuchos::RCP<const class ROL::Vector<double> > (ROL::BoundConstraint<double>::*)() const) &ROL::BoundConstraint<double>::getLowerBound, "C++: ROL::BoundConstraint<double>::getLowerBound() const --> const class Teuchos::RCP<const class ROL::Vector<double> >");
		cl.def("getUpperBound", (const class Teuchos::RCP<const class ROL::Vector<double> > (ROL::BoundConstraint<double>::*)() const) &ROL::BoundConstraint<double>::getUpperBound, "C++: ROL::BoundConstraint<double>::getUpperBound() const --> const class Teuchos::RCP<const class ROL::Vector<double> >");
		cl.def("isFeasible", (bool (ROL::BoundConstraint<double>::*)(const class ROL::Vector<double> &)) &ROL::BoundConstraint<double>::isFeasible, "C++: ROL::BoundConstraint<double>::isFeasible(const class ROL::Vector<double> &) --> bool", pybind11::arg("v"));
		cl.def("applyInverseScalingFunction", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) const) &ROL::BoundConstraint<double>::applyInverseScalingFunction, "C++: ROL::BoundConstraint<double>::applyInverseScalingFunction(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) const --> void", pybind11::arg("dv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("g"));
		cl.def("applyScalingFunctionJacobian", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) const) &ROL::BoundConstraint<double>::applyScalingFunctionJacobian, "C++: ROL::BoundConstraint<double>::applyScalingFunctionJacobian(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) const --> void", pybind11::arg("dv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("g"));
		cl.def("activateLower", (void (ROL::BoundConstraint<double>::*)()) &ROL::BoundConstraint<double>::activateLower, "C++: ROL::BoundConstraint<double>::activateLower() --> void");
		cl.def("activateUpper", (void (ROL::BoundConstraint<double>::*)()) &ROL::BoundConstraint<double>::activateUpper, "C++: ROL::BoundConstraint<double>::activateUpper() --> void");
		cl.def("activate", (void (ROL::BoundConstraint<double>::*)()) &ROL::BoundConstraint<double>::activate, "C++: ROL::BoundConstraint<double>::activate() --> void");
		cl.def("deactivateLower", (void (ROL::BoundConstraint<double>::*)()) &ROL::BoundConstraint<double>::deactivateLower, "C++: ROL::BoundConstraint<double>::deactivateLower() --> void");
		cl.def("deactivateUpper", (void (ROL::BoundConstraint<double>::*)()) &ROL::BoundConstraint<double>::deactivateUpper, "C++: ROL::BoundConstraint<double>::deactivateUpper() --> void");
		cl.def("deactivate", (void (ROL::BoundConstraint<double>::*)()) &ROL::BoundConstraint<double>::deactivate, "C++: ROL::BoundConstraint<double>::deactivate() --> void");
		cl.def("isLowerActivated", (bool (ROL::BoundConstraint<double>::*)() const) &ROL::BoundConstraint<double>::isLowerActivated, "C++: ROL::BoundConstraint<double>::isLowerActivated() const --> bool");
		cl.def("isUpperActivated", (bool (ROL::BoundConstraint<double>::*)() const) &ROL::BoundConstraint<double>::isUpperActivated, "C++: ROL::BoundConstraint<double>::isUpperActivated() const --> bool");
		cl.def("isActivated", (bool (ROL::BoundConstraint<double>::*)() const) &ROL::BoundConstraint<double>::isActivated, "C++: ROL::BoundConstraint<double>::isActivated() const --> bool");
		cl.def("pruneActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) -> void { return o.pruneActive(a0, a1); }, "", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("pruneActive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double)) &ROL::BoundConstraint<double>::pruneActive, "C++: ROL::BoundConstraint<double>::pruneActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, double) --> void", pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("eps"));
		cl.def("pruneActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> void { return o.pruneActive(a0, a1, a2); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"));
		cl.def("pruneActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double const & a3) -> void { return o.pruneActive(a0, a1, a2, a3); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"));
		cl.def("pruneActive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double)) &ROL::BoundConstraint<double>::pruneActive, "C++: ROL::BoundConstraint<double>::pruneActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double) --> void", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"), pybind11::arg("geps"));
		cl.def("pruneLowerInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) -> void { return o.pruneLowerInactive(a0, a1); }, "", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("pruneLowerInactive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double)) &ROL::BoundConstraint<double>::pruneLowerInactive, "C++: ROL::BoundConstraint<double>::pruneLowerInactive(class ROL::Vector<double> &, const class ROL::Vector<double> &, double) --> void", pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("eps"));
		cl.def("pruneUpperInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) -> void { return o.pruneUpperInactive(a0, a1); }, "", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("pruneUpperInactive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double)) &ROL::BoundConstraint<double>::pruneUpperInactive, "C++: ROL::BoundConstraint<double>::pruneUpperInactive(class ROL::Vector<double> &, const class ROL::Vector<double> &, double) --> void", pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("eps"));
		cl.def("pruneLowerInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> void { return o.pruneLowerInactive(a0, a1, a2); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"));
		cl.def("pruneLowerInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double const & a3) -> void { return o.pruneLowerInactive(a0, a1, a2, a3); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"));
		cl.def("pruneLowerInactive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double)) &ROL::BoundConstraint<double>::pruneLowerInactive, "C++: ROL::BoundConstraint<double>::pruneLowerInactive(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double) --> void", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"), pybind11::arg("geps"));
		cl.def("pruneUpperInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> void { return o.pruneUpperInactive(a0, a1, a2); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"));
		cl.def("pruneUpperInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double const & a3) -> void { return o.pruneUpperInactive(a0, a1, a2, a3); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"));
		cl.def("pruneUpperInactive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double)) &ROL::BoundConstraint<double>::pruneUpperInactive, "C++: ROL::BoundConstraint<double>::pruneUpperInactive(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double) --> void", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"), pybind11::arg("geps"));
		cl.def("pruneInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) -> void { return o.pruneInactive(a0, a1); }, "", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("pruneInactive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double)) &ROL::BoundConstraint<double>::pruneInactive, "C++: ROL::BoundConstraint<double>::pruneInactive(class ROL::Vector<double> &, const class ROL::Vector<double> &, double) --> void", pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("eps"));
		cl.def("pruneInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> void { return o.pruneInactive(a0, a1, a2); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"));
		cl.def("pruneInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double const & a3) -> void { return o.pruneInactive(a0, a1, a2, a3); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"));
		cl.def("pruneInactive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double)) &ROL::BoundConstraint<double>::pruneInactive, "C++: ROL::BoundConstraint<double>::pruneInactive(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double) --> void", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"), pybind11::arg("geps"));
		cl.def("computeProjectedGradient", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::BoundConstraint<double>::computeProjectedGradient, "C++: ROL::BoundConstraint<double>::computeProjectedGradient(class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("g"), pybind11::arg("x"));
		cl.def("computeProjectedStep", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::BoundConstraint<double>::computeProjectedStep, "C++: ROL::BoundConstraint<double>::computeProjectedStep(class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("assign", (class ROL::BoundConstraint<double> & (ROL::BoundConstraint<double>::*)(const class ROL::BoundConstraint<double> &)) &ROL::BoundConstraint<double>::operator=, "C++: ROL::BoundConstraint<double>::operator=(const class ROL::BoundConstraint<double> &) --> class ROL::BoundConstraint<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::BoundConstraint_Partitioned file:ROL_BoundConstraint_Partitioned.hpp line:61
		pybind11::class_<ROL::BoundConstraint_Partitioned<double>, Teuchos::RCP<ROL::BoundConstraint_Partitioned<double>>, PyCallBack_ROL_BoundConstraint_Partitioned_double_t, ROL::BoundConstraint<double>> cl(M("ROL"), "BoundConstraint_Partitioned_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](PyCallBack_ROL_BoundConstraint_Partitioned_double_t const &o){ return new PyCallBack_ROL_BoundConstraint_Partitioned_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::BoundConstraint_Partitioned<double> const &o){ return new ROL::BoundConstraint_Partitioned<double>(o); } ) );
		cl.def("update", [](ROL::BoundConstraint_Partitioned<double> &o, const class ROL::Vector<double> & a0) -> void { return o.update(a0); }, "", pybind11::arg("x"));
		cl.def("update", [](ROL::BoundConstraint_Partitioned<double> &o, const class ROL::Vector<double> & a0, bool const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("flag"));
		cl.def("update", (void (ROL::BoundConstraint_Partitioned<double>::*)(const class ROL::Vector<double> &, bool, int)) &ROL::BoundConstraint_Partitioned<double>::update, "C++: ROL::BoundConstraint_Partitioned<double>::update(const class ROL::Vector<double> &, bool, int) --> void", pybind11::arg("x"), pybind11::arg("flag"), pybind11::arg("iter"));
		cl.def("project", (void (ROL::BoundConstraint_Partitioned<double>::*)(class ROL::Vector<double> &)) &ROL::BoundConstraint_Partitioned<double>::project, "C++: ROL::BoundConstraint_Partitioned<double>::project(class ROL::Vector<double> &) --> void", pybind11::arg("x"));
		cl.def("projectInterior", (void (ROL::BoundConstraint_Partitioned<double>::*)(class ROL::Vector<double> &)) &ROL::BoundConstraint_Partitioned<double>::projectInterior, "C++: ROL::BoundConstraint_Partitioned<double>::projectInterior(class ROL::Vector<double> &) --> void", pybind11::arg("x"));
		cl.def("pruneUpperActive", [](ROL::BoundConstraint_Partitioned<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) -> void { return o.pruneUpperActive(a0, a1); }, "", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("pruneUpperActive", (void (ROL::BoundConstraint_Partitioned<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double)) &ROL::BoundConstraint_Partitioned<double>::pruneUpperActive, "C++: ROL::BoundConstraint_Partitioned<double>::pruneUpperActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, double) --> void", pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("eps"));
		cl.def("pruneUpperActive", [](ROL::BoundConstraint_Partitioned<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> void { return o.pruneUpperActive(a0, a1, a2); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"));
		cl.def("pruneUpperActive", [](ROL::BoundConstraint_Partitioned<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double const & a3) -> void { return o.pruneUpperActive(a0, a1, a2, a3); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"));
		cl.def("pruneUpperActive", (void (ROL::BoundConstraint_Partitioned<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double)) &ROL::BoundConstraint_Partitioned<double>::pruneUpperActive, "C++: ROL::BoundConstraint_Partitioned<double>::pruneUpperActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double) --> void", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"), pybind11::arg("geps"));
		cl.def("pruneLowerActive", [](ROL::BoundConstraint_Partitioned<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) -> void { return o.pruneLowerActive(a0, a1); }, "", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("pruneLowerActive", (void (ROL::BoundConstraint_Partitioned<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double)) &ROL::BoundConstraint_Partitioned<double>::pruneLowerActive, "C++: ROL::BoundConstraint_Partitioned<double>::pruneLowerActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, double) --> void", pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("eps"));
		cl.def("pruneLowerActive", [](ROL::BoundConstraint_Partitioned<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> void { return o.pruneLowerActive(a0, a1, a2); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"));
		cl.def("pruneLowerActive", [](ROL::BoundConstraint_Partitioned<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double const & a3) -> void { return o.pruneLowerActive(a0, a1, a2, a3); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"));
		cl.def("pruneLowerActive", (void (ROL::BoundConstraint_Partitioned<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double)) &ROL::BoundConstraint_Partitioned<double>::pruneLowerActive, "C++: ROL::BoundConstraint_Partitioned<double>::pruneLowerActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double) --> void", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"), pybind11::arg("geps"));
		cl.def("isFeasible", (bool (ROL::BoundConstraint_Partitioned<double>::*)(const class ROL::Vector<double> &)) &ROL::BoundConstraint_Partitioned<double>::isFeasible, "C++: ROL::BoundConstraint_Partitioned<double>::isFeasible(const class ROL::Vector<double> &) --> bool", pybind11::arg("v"));
		cl.def("applyInverseScalingFunction", (void (ROL::BoundConstraint_Partitioned<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) const) &ROL::BoundConstraint_Partitioned<double>::applyInverseScalingFunction, "C++: ROL::BoundConstraint_Partitioned<double>::applyInverseScalingFunction(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) const --> void", pybind11::arg("dv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("g"));
		cl.def("applyScalingFunctionJacobian", (void (ROL::BoundConstraint_Partitioned<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) const) &ROL::BoundConstraint_Partitioned<double>::applyScalingFunctionJacobian, "C++: ROL::BoundConstraint_Partitioned<double>::applyScalingFunctionJacobian(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) const --> void", pybind11::arg("dv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("g"));
		cl.def("assign", (class ROL::BoundConstraint_Partitioned<double> & (ROL::BoundConstraint_Partitioned<double>::*)(const class ROL::BoundConstraint_Partitioned<double> &)) &ROL::BoundConstraint_Partitioned<double>::operator=, "C++: ROL::BoundConstraint_Partitioned<double>::operator=(const class ROL::BoundConstraint_Partitioned<double> &) --> class ROL::BoundConstraint_Partitioned<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("project", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &)) &ROL::BoundConstraint<double>::project, "C++: ROL::BoundConstraint<double>::project(class ROL::Vector<double> &) --> void", pybind11::arg("x"));
		cl.def("projectInterior", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &)) &ROL::BoundConstraint<double>::projectInterior, "C++: ROL::BoundConstraint<double>::projectInterior(class ROL::Vector<double> &) --> void", pybind11::arg("x"));
		cl.def("pruneUpperActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) -> void { return o.pruneUpperActive(a0, a1); }, "", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("pruneUpperActive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double)) &ROL::BoundConstraint<double>::pruneUpperActive, "C++: ROL::BoundConstraint<double>::pruneUpperActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, double) --> void", pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("eps"));
		cl.def("pruneUpperActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> void { return o.pruneUpperActive(a0, a1, a2); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"));
		cl.def("pruneUpperActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double const & a3) -> void { return o.pruneUpperActive(a0, a1, a2, a3); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"));
		cl.def("pruneUpperActive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double)) &ROL::BoundConstraint<double>::pruneUpperActive, "C++: ROL::BoundConstraint<double>::pruneUpperActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double) --> void", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"), pybind11::arg("geps"));
		cl.def("pruneLowerActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) -> void { return o.pruneLowerActive(a0, a1); }, "", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("pruneLowerActive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double)) &ROL::BoundConstraint<double>::pruneLowerActive, "C++: ROL::BoundConstraint<double>::pruneLowerActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, double) --> void", pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("eps"));
		cl.def("pruneLowerActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> void { return o.pruneLowerActive(a0, a1, a2); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"));
		cl.def("pruneLowerActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double const & a3) -> void { return o.pruneLowerActive(a0, a1, a2, a3); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"));
		cl.def("pruneLowerActive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double)) &ROL::BoundConstraint<double>::pruneLowerActive, "C++: ROL::BoundConstraint<double>::pruneLowerActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double) --> void", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"), pybind11::arg("geps"));
		cl.def("getLowerBound", (const class Teuchos::RCP<const class ROL::Vector<double> > (ROL::BoundConstraint<double>::*)() const) &ROL::BoundConstraint<double>::getLowerBound, "C++: ROL::BoundConstraint<double>::getLowerBound() const --> const class Teuchos::RCP<const class ROL::Vector<double> >");
		cl.def("getUpperBound", (const class Teuchos::RCP<const class ROL::Vector<double> > (ROL::BoundConstraint<double>::*)() const) &ROL::BoundConstraint<double>::getUpperBound, "C++: ROL::BoundConstraint<double>::getUpperBound() const --> const class Teuchos::RCP<const class ROL::Vector<double> >");
		cl.def("isFeasible", (bool (ROL::BoundConstraint<double>::*)(const class ROL::Vector<double> &)) &ROL::BoundConstraint<double>::isFeasible, "C++: ROL::BoundConstraint<double>::isFeasible(const class ROL::Vector<double> &) --> bool", pybind11::arg("v"));
		cl.def("applyInverseScalingFunction", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) const) &ROL::BoundConstraint<double>::applyInverseScalingFunction, "C++: ROL::BoundConstraint<double>::applyInverseScalingFunction(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) const --> void", pybind11::arg("dv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("g"));
		cl.def("applyScalingFunctionJacobian", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) const) &ROL::BoundConstraint<double>::applyScalingFunctionJacobian, "C++: ROL::BoundConstraint<double>::applyScalingFunctionJacobian(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) const --> void", pybind11::arg("dv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("g"));
		cl.def("activateLower", (void (ROL::BoundConstraint<double>::*)()) &ROL::BoundConstraint<double>::activateLower, "C++: ROL::BoundConstraint<double>::activateLower() --> void");
		cl.def("activateUpper", (void (ROL::BoundConstraint<double>::*)()) &ROL::BoundConstraint<double>::activateUpper, "C++: ROL::BoundConstraint<double>::activateUpper() --> void");
		cl.def("activate", (void (ROL::BoundConstraint<double>::*)()) &ROL::BoundConstraint<double>::activate, "C++: ROL::BoundConstraint<double>::activate() --> void");
		cl.def("deactivateLower", (void (ROL::BoundConstraint<double>::*)()) &ROL::BoundConstraint<double>::deactivateLower, "C++: ROL::BoundConstraint<double>::deactivateLower() --> void");
		cl.def("deactivateUpper", (void (ROL::BoundConstraint<double>::*)()) &ROL::BoundConstraint<double>::deactivateUpper, "C++: ROL::BoundConstraint<double>::deactivateUpper() --> void");
		cl.def("deactivate", (void (ROL::BoundConstraint<double>::*)()) &ROL::BoundConstraint<double>::deactivate, "C++: ROL::BoundConstraint<double>::deactivate() --> void");
		cl.def("isLowerActivated", (bool (ROL::BoundConstraint<double>::*)() const) &ROL::BoundConstraint<double>::isLowerActivated, "C++: ROL::BoundConstraint<double>::isLowerActivated() const --> bool");
		cl.def("isUpperActivated", (bool (ROL::BoundConstraint<double>::*)() const) &ROL::BoundConstraint<double>::isUpperActivated, "C++: ROL::BoundConstraint<double>::isUpperActivated() const --> bool");
		cl.def("isActivated", (bool (ROL::BoundConstraint<double>::*)() const) &ROL::BoundConstraint<double>::isActivated, "C++: ROL::BoundConstraint<double>::isActivated() const --> bool");
		cl.def("pruneActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) -> void { return o.pruneActive(a0, a1); }, "", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("pruneActive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double)) &ROL::BoundConstraint<double>::pruneActive, "C++: ROL::BoundConstraint<double>::pruneActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, double) --> void", pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("eps"));
		cl.def("pruneActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> void { return o.pruneActive(a0, a1, a2); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"));
		cl.def("pruneActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double const & a3) -> void { return o.pruneActive(a0, a1, a2, a3); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"));
		cl.def("pruneActive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double)) &ROL::BoundConstraint<double>::pruneActive, "C++: ROL::BoundConstraint<double>::pruneActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double) --> void", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"), pybind11::arg("geps"));
		cl.def("pruneLowerInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) -> void { return o.pruneLowerInactive(a0, a1); }, "", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("pruneLowerInactive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double)) &ROL::BoundConstraint<double>::pruneLowerInactive, "C++: ROL::BoundConstraint<double>::pruneLowerInactive(class ROL::Vector<double> &, const class ROL::Vector<double> &, double) --> void", pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("eps"));
		cl.def("pruneUpperInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) -> void { return o.pruneUpperInactive(a0, a1); }, "", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("pruneUpperInactive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double)) &ROL::BoundConstraint<double>::pruneUpperInactive, "C++: ROL::BoundConstraint<double>::pruneUpperInactive(class ROL::Vector<double> &, const class ROL::Vector<double> &, double) --> void", pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("eps"));
		cl.def("pruneLowerInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> void { return o.pruneLowerInactive(a0, a1, a2); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"));
		cl.def("pruneLowerInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double const & a3) -> void { return o.pruneLowerInactive(a0, a1, a2, a3); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"));
		cl.def("pruneLowerInactive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double)) &ROL::BoundConstraint<double>::pruneLowerInactive, "C++: ROL::BoundConstraint<double>::pruneLowerInactive(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double) --> void", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"), pybind11::arg("geps"));
		cl.def("pruneUpperInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> void { return o.pruneUpperInactive(a0, a1, a2); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"));
		cl.def("pruneUpperInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double const & a3) -> void { return o.pruneUpperInactive(a0, a1, a2, a3); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"));
		cl.def("pruneUpperInactive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double)) &ROL::BoundConstraint<double>::pruneUpperInactive, "C++: ROL::BoundConstraint<double>::pruneUpperInactive(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double) --> void", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"), pybind11::arg("geps"));
		cl.def("pruneInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) -> void { return o.pruneInactive(a0, a1); }, "", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("pruneInactive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double)) &ROL::BoundConstraint<double>::pruneInactive, "C++: ROL::BoundConstraint<double>::pruneInactive(class ROL::Vector<double> &, const class ROL::Vector<double> &, double) --> void", pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("eps"));
		cl.def("pruneInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> void { return o.pruneInactive(a0, a1, a2); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"));
		cl.def("pruneInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double const & a3) -> void { return o.pruneInactive(a0, a1, a2, a3); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"));
		cl.def("pruneInactive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double)) &ROL::BoundConstraint<double>::pruneInactive, "C++: ROL::BoundConstraint<double>::pruneInactive(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double) --> void", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"), pybind11::arg("geps"));
		cl.def("computeProjectedGradient", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::BoundConstraint<double>::computeProjectedGradient, "C++: ROL::BoundConstraint<double>::computeProjectedGradient(class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("g"), pybind11::arg("x"));
		cl.def("computeProjectedStep", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::BoundConstraint<double>::computeProjectedStep, "C++: ROL::BoundConstraint<double>::computeProjectedStep(class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("assign", (class ROL::BoundConstraint<double> & (ROL::BoundConstraint<double>::*)(const class ROL::BoundConstraint<double> &)) &ROL::BoundConstraint<double>::operator=, "C++: ROL::BoundConstraint<double>::operator=(const class ROL::BoundConstraint<double> &) --> class ROL::BoundConstraint<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::ConstraintData file:ROL_ConstraintAssembler.hpp line:61
		pybind11::class_<ROL::ConstraintData<double>, Teuchos::RCP<ROL::ConstraintData<double>>> cl(M("ROL"), "ConstraintData_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](const class Teuchos::RCP<class ROL::Constraint<double> > & a0, const class Teuchos::RCP<class ROL::Vector<double> > & a1){ return new ROL::ConstraintData<double>(a0, a1); } ), "doc" , pybind11::arg("con"), pybind11::arg("mul"));
		cl.def( pybind11::init( [](const class Teuchos::RCP<class ROL::Constraint<double> > & a0, const class Teuchos::RCP<class ROL::Vector<double> > & a1, const class Teuchos::RCP<class ROL::Vector<double> > & a2){ return new ROL::ConstraintData<double>(a0, a1, a2); } ), "doc" , pybind11::arg("con"), pybind11::arg("mul"), pybind11::arg("res"));
		cl.def( pybind11::init<const class Teuchos::RCP<class ROL::Constraint<double> > &, const class Teuchos::RCP<class ROL::Vector<double> > &, const class Teuchos::RCP<class ROL::Vector<double> > &, const class Teuchos::RCP<class ROL::BoundConstraint<double> > &>(), pybind11::arg("con"), pybind11::arg("mul"), pybind11::arg("res"), pybind11::arg("bnd") );

		cl.def( pybind11::init( [](ROL::ConstraintData<double> const &o){ return new ROL::ConstraintData<double>(o); } ) );
		cl.def_readonly("constraint", &ROL::ConstraintData<double>::constraint);
		cl.def_readonly("multiplier", &ROL::ConstraintData<double>::multiplier);
		cl.def_readonly("residual", &ROL::ConstraintData<double>::residual);
		cl.def_readonly("bounds", &ROL::ConstraintData<double>::bounds);
	}
	{ // ROL::SlacklessObjective file:ROL_SlacklessObjective.hpp line:59
		pybind11::class_<ROL::SlacklessObjective<double>, Teuchos::RCP<ROL::SlacklessObjective<double>>, PyCallBack_ROL_SlacklessObjective_double_t, ROL::Objective<double>> cl(M("ROL"), "SlacklessObjective_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<const class Teuchos::RCP<class ROL::Objective<double> > &>(), pybind11::arg("obj") );

		cl.def( pybind11::init( [](PyCallBack_ROL_SlacklessObjective_double_t const &o){ return new PyCallBack_ROL_SlacklessObjective_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::SlacklessObjective<double> const &o){ return new ROL::SlacklessObjective<double>(o); } ) );
		cl.def("getObjective", (class Teuchos::RCP<class ROL::Objective<double> > (ROL::SlacklessObjective<double>::*)() const) &ROL::SlacklessObjective<double>::getObjective, "C++: ROL::SlacklessObjective<double>::getObjective() const --> class Teuchos::RCP<class ROL::Objective<double> >");
		cl.def("update", [](ROL::SlacklessObjective<double> &o, const class ROL::Vector<double> & a0, enum ROL::UpdateType const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("type"));
		cl.def("update", (void (ROL::SlacklessObjective<double>::*)(const class ROL::Vector<double> &, enum ROL::UpdateType, int)) &ROL::SlacklessObjective<double>::update, "C++: ROL::SlacklessObjective<double>::update(const class ROL::Vector<double> &, enum ROL::UpdateType, int) --> void", pybind11::arg("x"), pybind11::arg("type"), pybind11::arg("iter"));
		cl.def("update", [](ROL::SlacklessObjective<double> &o, const class ROL::Vector<double> & a0) -> void { return o.update(a0); }, "", pybind11::arg("x"));
		cl.def("update", [](ROL::SlacklessObjective<double> &o, const class ROL::Vector<double> & a0, bool const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("flag"));
		cl.def("update", (void (ROL::SlacklessObjective<double>::*)(const class ROL::Vector<double> &, bool, int)) &ROL::SlacklessObjective<double>::update, "C++: ROL::SlacklessObjective<double>::update(const class ROL::Vector<double> &, bool, int) --> void", pybind11::arg("x"), pybind11::arg("flag"), pybind11::arg("iter"));
		cl.def("value", (double (ROL::SlacklessObjective<double>::*)(const class ROL::Vector<double> &, double &)) &ROL::SlacklessObjective<double>::value, "C++: ROL::SlacklessObjective<double>::value(const class ROL::Vector<double> &, double &) --> double", pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("dirDeriv", (double (ROL::SlacklessObjective<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::SlacklessObjective<double>::dirDeriv, "C++: ROL::SlacklessObjective<double>::dirDeriv(const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> double", pybind11::arg("x"), pybind11::arg("d"), pybind11::arg("tol"));
		cl.def("gradient", (void (ROL::SlacklessObjective<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::SlacklessObjective<double>::gradient, "C++: ROL::SlacklessObjective<double>::gradient(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("hessVec", (void (ROL::SlacklessObjective<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::SlacklessObjective<double>::hessVec, "C++: ROL::SlacklessObjective<double>::hessVec(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("hv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("invHessVec", (void (ROL::SlacklessObjective<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::SlacklessObjective<double>::invHessVec, "C++: ROL::SlacklessObjective<double>::invHessVec(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("ihv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("precond", (void (ROL::SlacklessObjective<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::SlacklessObjective<double>::precond, "C++: ROL::SlacklessObjective<double>::precond(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("Pv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("update", [](ROL::Objective<double> &o, const class ROL::Vector<double> & a0, enum ROL::UpdateType const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("type"));
		cl.def("update", (void (ROL::Objective<double>::*)(const class ROL::Vector<double> &, enum ROL::UpdateType, int)) &ROL::Objective<double>::update, "Update objective function. \n\n      This function updates the objective function at new iterations. \n      \n\n      is the new iterate. \n      \n\n   is the type of update requested.\n      \n\n   is the outer algorithm iterations count.\n\nC++: ROL::Objective<double>::update(const class ROL::Vector<double> &, enum ROL::UpdateType, int) --> void", pybind11::arg("x"), pybind11::arg("type"), pybind11::arg("iter"));
		cl.def("update", [](ROL::Objective<double> &o, const class ROL::Vector<double> & a0) -> void { return o.update(a0); }, "", pybind11::arg("x"));
		cl.def("update", [](ROL::Objective<double> &o, const class ROL::Vector<double> & a0, bool const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("flag"));
		cl.def("update", (void (ROL::Objective<double>::*)(const class ROL::Vector<double> &, bool, int)) &ROL::Objective<double>::update, "Update objective function. \n\n      This function updates the objective function at new iterations. \n      \n\n      is the new iterate. \n      \n\n   is true if the iterate has changed.\n      \n\n   is the outer algorithm iterations count.\n\nC++: ROL::Objective<double>::update(const class ROL::Vector<double> &, bool, int) --> void", pybind11::arg("x"), pybind11::arg("flag"), pybind11::arg("iter"));
		cl.def("value", (double (ROL::Objective<double>::*)(const class ROL::Vector<double> &, double &)) &ROL::Objective<double>::value, "C++: ROL::Objective<double>::value(const class ROL::Vector<double> &, double &) --> double", pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("gradient", (void (ROL::Objective<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Objective<double>::gradient, "Compute gradient.\n\n      This function returns the objective function gradient.\n      \n\n   is the gradient.\n      \n\n   is the current iterate.\n      \n\n is a tolerance for inexact objective function computation.\n\n      The default implementation is a finite-difference approximation based on the function value.\n      This requires the definition of a basis \n\n for the optimization vectors x and\n      the definition of a basis \n\n for the dual optimization vectors (gradient vectors g).\n      The bases must be related through the Riesz map, i.e., \n\n,\n      and this must be reflected in the implementation of the ROL::Vector::dual() method.\n\nC++: ROL::Objective<double>::gradient(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("dirDeriv", (double (ROL::Objective<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Objective<double>::dirDeriv, "Compute directional derivative.\n\n      This function returns the directional derivative of the objective function in the \n direction.\n      \n\n   is the current iterate.\n      \n\n   is the direction.\n      \n\n is a tolerance for inexact objective function computation.\n\nC++: ROL::Objective<double>::dirDeriv(const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> double", pybind11::arg("x"), pybind11::arg("d"), pybind11::arg("tol"));
		cl.def("hessVec", (void (ROL::Objective<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Objective<double>::hessVec, "Apply Hessian approximation to vector.\n\n      This function applies the Hessian of the objective function to the vector \n.\n      \n\n  is the the action of the Hessian on \n.\n      \n\n   is the direction vector.\n      \n\n   is the current iterate.\n      \n\n is a tolerance for inexact objective function computation.\n\nC++: ROL::Objective<double>::hessVec(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("hv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("invHessVec", (void (ROL::Objective<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Objective<double>::invHessVec, "Apply inverse Hessian approximation to vector.\n\n      This function applies the inverse Hessian of the objective function to the vector \n.\n      \n\n  is the action of the inverse Hessian on \n.\n      \n\n   is the direction vector.\n      \n\n   is the current iterate.\n      \n\n is a tolerance for inexact objective function computation.\n\nC++: ROL::Objective<double>::invHessVec(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("hv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("precond", (void (ROL::Objective<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Objective<double>::precond, "Apply preconditioner to vector.\n\n      This function applies a preconditioner for the Hessian of the objective function to the vector \n.\n      \n\n  is the action of the Hessian preconditioner on \n.\n      \n\n   is the direction vector.\n      \n\n   is the current iterate.\n      \n\n is a tolerance for inexact objective function computation.\n\nC++: ROL::Objective<double>::precond(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("Pv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("assign", (class ROL::Objective<double> & (ROL::Objective<double>::*)(const class ROL::Objective<double> &)) &ROL::Objective<double>::operator=, "C++: ROL::Objective<double>::operator=(const class ROL::Objective<double> &) --> class ROL::Objective<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::LinearConstraint file:ROL_LinearConstraint.hpp line:60
		pybind11::class_<ROL::LinearConstraint<double>, Teuchos::RCP<ROL::LinearConstraint<double>>, PyCallBack_ROL_LinearConstraint_double_t, ROL::Constraint<double>> cl(M("ROL"), "LinearConstraint_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<const class Teuchos::RCP<const class ROL::LinearOperator<double> > &, const class Teuchos::RCP<const class ROL::Vector<double> > &>(), pybind11::arg("A"), pybind11::arg("b") );

		cl.def( pybind11::init( [](PyCallBack_ROL_LinearConstraint_double_t const &o){ return new PyCallBack_ROL_LinearConstraint_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::LinearConstraint<double> const &o){ return new ROL::LinearConstraint<double>(o); } ) );
		cl.def("update", [](ROL::LinearConstraint<double> &o, const class ROL::Vector<double> & a0, enum ROL::UpdateType const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("type"));
		cl.def("update", (void (ROL::LinearConstraint<double>::*)(const class ROL::Vector<double> &, enum ROL::UpdateType, int)) &ROL::LinearConstraint<double>::update, "C++: ROL::LinearConstraint<double>::update(const class ROL::Vector<double> &, enum ROL::UpdateType, int) --> void", pybind11::arg("x"), pybind11::arg("type"), pybind11::arg("iter"));
		cl.def("update", [](ROL::LinearConstraint<double> &o, const class ROL::Vector<double> & a0) -> void { return o.update(a0); }, "", pybind11::arg("x"));
		cl.def("update", [](ROL::LinearConstraint<double> &o, const class ROL::Vector<double> & a0, bool const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("flag"));
		cl.def("update", (void (ROL::LinearConstraint<double>::*)(const class ROL::Vector<double> &, bool, int)) &ROL::LinearConstraint<double>::update, "C++: ROL::LinearConstraint<double>::update(const class ROL::Vector<double> &, bool, int) --> void", pybind11::arg("x"), pybind11::arg("flag"), pybind11::arg("iter"));
		cl.def("value", (void (ROL::LinearConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::LinearConstraint<double>::value, "C++: ROL::LinearConstraint<double>::value(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("c"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("applyJacobian", (void (ROL::LinearConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::LinearConstraint<double>::applyJacobian, "C++: ROL::LinearConstraint<double>::applyJacobian(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("jv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("applyAdjointJacobian", (void (ROL::LinearConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::LinearConstraint<double>::applyAdjointJacobian, "C++: ROL::LinearConstraint<double>::applyAdjointJacobian(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("ajv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("applyAdjointJacobian", (void (ROL::LinearConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::LinearConstraint<double>::applyAdjointJacobian, "C++: ROL::LinearConstraint<double>::applyAdjointJacobian(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("ajv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("dualv"), pybind11::arg("tol"));
		cl.def("applyAdjointHessian", (void (ROL::LinearConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::LinearConstraint<double>::applyAdjointHessian, "C++: ROL::LinearConstraint<double>::applyAdjointHessian(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("ahuv"), pybind11::arg("u"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("createRangeSpaceVector", (class Teuchos::RCP<class ROL::Vector<double> > (ROL::LinearConstraint<double>::*)() const) &ROL::LinearConstraint<double>::createRangeSpaceVector, "C++: ROL::LinearConstraint<double>::createRangeSpaceVector() const --> class Teuchos::RCP<class ROL::Vector<double> >");
		cl.def("update", [](ROL::Constraint<double> &o, const class ROL::Vector<double> & a0, enum ROL::UpdateType const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("type"));
		cl.def("update", (void (ROL::Constraint<double>::*)(const class ROL::Vector<double> &, enum ROL::UpdateType, int)) &ROL::Constraint<double>::update, "Update constraint function. \n\n      This function updates the constraint function at new iterations. \n      \n\n      is the new iterate. \n      \n\n   is the type of update requested.\n      \n\n   is the outer algorithm iterations count.\n\nC++: ROL::Constraint<double>::update(const class ROL::Vector<double> &, enum ROL::UpdateType, int) --> void", pybind11::arg("x"), pybind11::arg("type"), pybind11::arg("iter"));
		cl.def("update", [](ROL::Constraint<double> &o, const class ROL::Vector<double> & a0) -> void { return o.update(a0); }, "", pybind11::arg("x"));
		cl.def("update", [](ROL::Constraint<double> &o, const class ROL::Vector<double> & a0, bool const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("flag"));
		cl.def("update", (void (ROL::Constraint<double>::*)(const class ROL::Vector<double> &, bool, int)) &ROL::Constraint<double>::update, "Update constraint functions.  \n                x is the optimization variable, \n                flag = true if optimization variable is changed,\n                iter is the outer algorithm iterations count.\n\nC++: ROL::Constraint<double>::update(const class ROL::Vector<double> &, bool, int) --> void", pybind11::arg("x"), pybind11::arg("flag"), pybind11::arg("iter"));
		cl.def("value", (void (ROL::Constraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Constraint<double>::value, "C++: ROL::Constraint<double>::value(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("c"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("applyJacobian", (void (ROL::Constraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Constraint<double>::applyJacobian, "Apply the constraint Jacobian at \n, \n,\n             to vector \n\n.\n\n             \n  is the result of applying the constraint Jacobian to  at  a constraint-space vector\n             \n\n   is an optimization-space vector\n             \n\n   is the constraint argument; an optimization-space vector\n             \n\n is a tolerance for inexact evaluations; currently unused\n\n             On return, \n, where\n             \n\n, \n. \n             The default implementation is a finite-difference approximation.\n\n             ---\n\nC++: ROL::Constraint<double>::applyJacobian(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("jv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("applyAdjointJacobian", (void (ROL::Constraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Constraint<double>::applyAdjointJacobian, "Apply the adjoint of the the constraint Jacobian at \n, \n,\n             to vector \n\n.\n\n             \n is the result of applying the adjoint of the constraint Jacobian to  at  a dual optimization-space vector\n             \n\n   is a dual constraint-space vector\n             \n\n   is the constraint argument; an optimization-space vector\n             \n\n is a tolerance for inexact evaluations; currently unused\n\n             On return, \n, where\n             \n\n, \n. \n             The default implementation is a finite-difference approximation.\n\n             ---\n\nC++: ROL::Constraint<double>::applyAdjointJacobian(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("ajv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("applyAdjointJacobian", (void (ROL::Constraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Constraint<double>::applyAdjointJacobian, "Apply the adjoint of the the constraint Jacobian at \n, \n,\n             to vector \n\n.\n\n             \n is the result of applying the adjoint of the constraint Jacobian to  at  a dual optimization-space vector\n             \n\n   is a dual constraint-space vector\n             \n\n   is the constraint argument; an optimization-space vector\n             \n\n  is a vector used for temporary variables; a constraint-space vector\n             \n\n is a tolerance for inexact evaluations; currently unused\n\n             On return, \n, where\n             \n\n, \n. \n             The default implementation is a finite-difference approximation.\n\n             ---\n\nC++: ROL::Constraint<double>::applyAdjointJacobian(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("ajv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("dualv"), pybind11::arg("tol"));
		cl.def("applyAdjointHessian", (void (ROL::Constraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Constraint<double>::applyAdjointHessian, "Apply the derivative of the adjoint of the constraint Jacobian at \n\n             to vector \n in direction \n,\n             according to \n\n.\n\n             \n is the result of applying the derivative of the adjoint of the constraint Jacobian at  to vector  in direction  a dual optimization-space vector\n             \n\n    is the direction vector; a dual constraint-space vector\n             \n\n    is an optimization-space vector\n             \n\n    is the constraint argument; an optimization-space vector\n             \n\n  is a tolerance for inexact evaluations; currently unused\n\n             On return, \n, where\n             \n\n, \n, and \n. \n             The default implementation is a finite-difference approximation based on the adjoint Jacobian.\n\n             ---\n\nC++: ROL::Constraint<double>::applyAdjointHessian(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("huv"), pybind11::arg("u"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("applyPreconditioner", (void (ROL::Constraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Constraint<double>::applyPreconditioner, "Apply a constraint preconditioner at \n, \n,\n             to vector \n\n.  Ideally, this preconditioner satisfies the following relationship:\n             \n\n\n\n             where R is the appropriate Riesz map in \n.  It is used by the #solveAugmentedSystem method.\n\n             \n  is the result of applying the constraint preconditioner to  at  a dual constraint-space vector\n             \n\n   is a constraint-space vector\n             \n\n   is the preconditioner argument; an optimization-space vector\n             \n\n   is the preconditioner argument; a dual optimization-space vector, unused\n             \n\n is a tolerance for inexact evaluations\n\n             On return, \n, where\n             \n\n, \n. \n             The default implementation is the Riesz map in \n\n.\n\n             ---\n\nC++: ROL::Constraint<double>::applyPreconditioner(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("pv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("tol"));
		cl.def("activate", (void (ROL::Constraint<double>::*)()) &ROL::Constraint<double>::activate, "Turn on constraints \n\nC++: ROL::Constraint<double>::activate() --> void");
		cl.def("deactivate", (void (ROL::Constraint<double>::*)()) &ROL::Constraint<double>::deactivate, "Turn off constraints\n\nC++: ROL::Constraint<double>::deactivate() --> void");
		cl.def("isActivated", (bool (ROL::Constraint<double>::*)()) &ROL::Constraint<double>::isActivated, "Check if constraints are on\n\nC++: ROL::Constraint<double>::isActivated() --> bool");
		cl.def("assign", (class ROL::Constraint<double> & (ROL::Constraint<double>::*)(const class ROL::Constraint<double> &)) &ROL::Constraint<double>::operator=, "C++: ROL::Constraint<double>::operator=(const class ROL::Constraint<double> &) --> class ROL::Constraint<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::Bounds file:ROL_Bounds.hpp line:59
		pybind11::class_<ROL::Bounds<double>, Teuchos::RCP<ROL::Bounds<double>>, PyCallBack_ROL_Bounds_double_t, ROL::BoundConstraint<double>> cl(M("ROL"), "Bounds_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](const class ROL::Vector<double> & a0){ return new ROL::Bounds<double>(a0); }, [](const class ROL::Vector<double> & a0){ return new PyCallBack_ROL_Bounds_double_t(a0); } ), "doc");
		cl.def( pybind11::init( [](const class ROL::Vector<double> & a0, bool const & a1){ return new ROL::Bounds<double>(a0, a1); }, [](const class ROL::Vector<double> & a0, bool const & a1){ return new PyCallBack_ROL_Bounds_double_t(a0, a1); } ), "doc");
		cl.def( pybind11::init( [](const class ROL::Vector<double> & a0, bool const & a1, double const & a2){ return new ROL::Bounds<double>(a0, a1, a2); }, [](const class ROL::Vector<double> & a0, bool const & a1, double const & a2){ return new PyCallBack_ROL_Bounds_double_t(a0, a1, a2); } ), "doc");
		cl.def( pybind11::init<const class ROL::Vector<double> &, bool, double, double>(), pybind11::arg("x"), pybind11::arg("isLower"), pybind11::arg("scale"), pybind11::arg("feasTol") );

		cl.def( pybind11::init( [](const class Teuchos::RCP<class ROL::Vector<double> > & a0, const class Teuchos::RCP<class ROL::Vector<double> > & a1){ return new ROL::Bounds<double>(a0, a1); }, [](const class Teuchos::RCP<class ROL::Vector<double> > & a0, const class Teuchos::RCP<class ROL::Vector<double> > & a1){ return new PyCallBack_ROL_Bounds_double_t(a0, a1); } ), "doc");
		cl.def( pybind11::init( [](const class Teuchos::RCP<class ROL::Vector<double> > & a0, const class Teuchos::RCP<class ROL::Vector<double> > & a1, const double & a2){ return new ROL::Bounds<double>(a0, a1, a2); }, [](const class Teuchos::RCP<class ROL::Vector<double> > & a0, const class Teuchos::RCP<class ROL::Vector<double> > & a1, const double & a2){ return new PyCallBack_ROL_Bounds_double_t(a0, a1, a2); } ), "doc");
		cl.def( pybind11::init<const class Teuchos::RCP<class ROL::Vector<double> > &, const class Teuchos::RCP<class ROL::Vector<double> > &, const double, const double>(), pybind11::arg("x_lo"), pybind11::arg("x_up"), pybind11::arg("scale"), pybind11::arg("feasTol") );

		cl.def( pybind11::init( [](PyCallBack_ROL_Bounds_double_t const &o){ return new PyCallBack_ROL_Bounds_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::Bounds<double> const &o){ return new ROL::Bounds<double>(o); } ) );
		cl.def("project", (void (ROL::Bounds<double>::*)(class ROL::Vector<double> &)) &ROL::Bounds<double>::project, "C++: ROL::Bounds<double>::project(class ROL::Vector<double> &) --> void", pybind11::arg("x"));
		cl.def("projectInterior", (void (ROL::Bounds<double>::*)(class ROL::Vector<double> &)) &ROL::Bounds<double>::projectInterior, "C++: ROL::Bounds<double>::projectInterior(class ROL::Vector<double> &) --> void", pybind11::arg("x"));
		cl.def("pruneUpperActive", [](ROL::Bounds<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) -> void { return o.pruneUpperActive(a0, a1); }, "", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("pruneUpperActive", (void (ROL::Bounds<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double)) &ROL::Bounds<double>::pruneUpperActive, "C++: ROL::Bounds<double>::pruneUpperActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, double) --> void", pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("eps"));
		cl.def("pruneUpperActive", [](ROL::Bounds<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> void { return o.pruneUpperActive(a0, a1, a2); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"));
		cl.def("pruneUpperActive", [](ROL::Bounds<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double const & a3) -> void { return o.pruneUpperActive(a0, a1, a2, a3); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"));
		cl.def("pruneUpperActive", (void (ROL::Bounds<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double)) &ROL::Bounds<double>::pruneUpperActive, "C++: ROL::Bounds<double>::pruneUpperActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double) --> void", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"), pybind11::arg("geps"));
		cl.def("pruneLowerActive", [](ROL::Bounds<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) -> void { return o.pruneLowerActive(a0, a1); }, "", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("pruneLowerActive", (void (ROL::Bounds<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double)) &ROL::Bounds<double>::pruneLowerActive, "C++: ROL::Bounds<double>::pruneLowerActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, double) --> void", pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("eps"));
		cl.def("pruneLowerActive", [](ROL::Bounds<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> void { return o.pruneLowerActive(a0, a1, a2); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"));
		cl.def("pruneLowerActive", [](ROL::Bounds<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double const & a3) -> void { return o.pruneLowerActive(a0, a1, a2, a3); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"));
		cl.def("pruneLowerActive", (void (ROL::Bounds<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double)) &ROL::Bounds<double>::pruneLowerActive, "C++: ROL::Bounds<double>::pruneLowerActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double) --> void", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"), pybind11::arg("geps"));
		cl.def("isFeasible", (bool (ROL::Bounds<double>::*)(const class ROL::Vector<double> &)) &ROL::Bounds<double>::isFeasible, "C++: ROL::Bounds<double>::isFeasible(const class ROL::Vector<double> &) --> bool", pybind11::arg("v"));
		cl.def("applyInverseScalingFunction", (void (ROL::Bounds<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) const) &ROL::Bounds<double>::applyInverseScalingFunction, "C++: ROL::Bounds<double>::applyInverseScalingFunction(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) const --> void", pybind11::arg("dv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("g"));
		cl.def("applyScalingFunctionJacobian", (void (ROL::Bounds<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) const) &ROL::Bounds<double>::applyScalingFunctionJacobian, "C++: ROL::Bounds<double>::applyScalingFunctionJacobian(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) const --> void", pybind11::arg("dv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("g"));
		cl.def("project", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &)) &ROL::BoundConstraint<double>::project, "C++: ROL::BoundConstraint<double>::project(class ROL::Vector<double> &) --> void", pybind11::arg("x"));
		cl.def("projectInterior", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &)) &ROL::BoundConstraint<double>::projectInterior, "C++: ROL::BoundConstraint<double>::projectInterior(class ROL::Vector<double> &) --> void", pybind11::arg("x"));
		cl.def("pruneUpperActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) -> void { return o.pruneUpperActive(a0, a1); }, "", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("pruneUpperActive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double)) &ROL::BoundConstraint<double>::pruneUpperActive, "C++: ROL::BoundConstraint<double>::pruneUpperActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, double) --> void", pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("eps"));
		cl.def("pruneUpperActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> void { return o.pruneUpperActive(a0, a1, a2); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"));
		cl.def("pruneUpperActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double const & a3) -> void { return o.pruneUpperActive(a0, a1, a2, a3); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"));
		cl.def("pruneUpperActive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double)) &ROL::BoundConstraint<double>::pruneUpperActive, "C++: ROL::BoundConstraint<double>::pruneUpperActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double) --> void", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"), pybind11::arg("geps"));
		cl.def("pruneLowerActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) -> void { return o.pruneLowerActive(a0, a1); }, "", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("pruneLowerActive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double)) &ROL::BoundConstraint<double>::pruneLowerActive, "C++: ROL::BoundConstraint<double>::pruneLowerActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, double) --> void", pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("eps"));
		cl.def("pruneLowerActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> void { return o.pruneLowerActive(a0, a1, a2); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"));
		cl.def("pruneLowerActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double const & a3) -> void { return o.pruneLowerActive(a0, a1, a2, a3); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"));
		cl.def("pruneLowerActive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double)) &ROL::BoundConstraint<double>::pruneLowerActive, "C++: ROL::BoundConstraint<double>::pruneLowerActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double) --> void", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"), pybind11::arg("geps"));
		cl.def("getLowerBound", (const class Teuchos::RCP<const class ROL::Vector<double> > (ROL::BoundConstraint<double>::*)() const) &ROL::BoundConstraint<double>::getLowerBound, "C++: ROL::BoundConstraint<double>::getLowerBound() const --> const class Teuchos::RCP<const class ROL::Vector<double> >");
		cl.def("getUpperBound", (const class Teuchos::RCP<const class ROL::Vector<double> > (ROL::BoundConstraint<double>::*)() const) &ROL::BoundConstraint<double>::getUpperBound, "C++: ROL::BoundConstraint<double>::getUpperBound() const --> const class Teuchos::RCP<const class ROL::Vector<double> >");
		cl.def("isFeasible", (bool (ROL::BoundConstraint<double>::*)(const class ROL::Vector<double> &)) &ROL::BoundConstraint<double>::isFeasible, "C++: ROL::BoundConstraint<double>::isFeasible(const class ROL::Vector<double> &) --> bool", pybind11::arg("v"));
		cl.def("applyInverseScalingFunction", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) const) &ROL::BoundConstraint<double>::applyInverseScalingFunction, "C++: ROL::BoundConstraint<double>::applyInverseScalingFunction(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) const --> void", pybind11::arg("dv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("g"));
		cl.def("applyScalingFunctionJacobian", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) const) &ROL::BoundConstraint<double>::applyScalingFunctionJacobian, "C++: ROL::BoundConstraint<double>::applyScalingFunctionJacobian(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) const --> void", pybind11::arg("dv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("g"));
		cl.def("activateLower", (void (ROL::BoundConstraint<double>::*)()) &ROL::BoundConstraint<double>::activateLower, "C++: ROL::BoundConstraint<double>::activateLower() --> void");
		cl.def("activateUpper", (void (ROL::BoundConstraint<double>::*)()) &ROL::BoundConstraint<double>::activateUpper, "C++: ROL::BoundConstraint<double>::activateUpper() --> void");
		cl.def("activate", (void (ROL::BoundConstraint<double>::*)()) &ROL::BoundConstraint<double>::activate, "C++: ROL::BoundConstraint<double>::activate() --> void");
		cl.def("deactivateLower", (void (ROL::BoundConstraint<double>::*)()) &ROL::BoundConstraint<double>::deactivateLower, "C++: ROL::BoundConstraint<double>::deactivateLower() --> void");
		cl.def("deactivateUpper", (void (ROL::BoundConstraint<double>::*)()) &ROL::BoundConstraint<double>::deactivateUpper, "C++: ROL::BoundConstraint<double>::deactivateUpper() --> void");
		cl.def("deactivate", (void (ROL::BoundConstraint<double>::*)()) &ROL::BoundConstraint<double>::deactivate, "C++: ROL::BoundConstraint<double>::deactivate() --> void");
		cl.def("isLowerActivated", (bool (ROL::BoundConstraint<double>::*)() const) &ROL::BoundConstraint<double>::isLowerActivated, "C++: ROL::BoundConstraint<double>::isLowerActivated() const --> bool");
		cl.def("isUpperActivated", (bool (ROL::BoundConstraint<double>::*)() const) &ROL::BoundConstraint<double>::isUpperActivated, "C++: ROL::BoundConstraint<double>::isUpperActivated() const --> bool");
		cl.def("isActivated", (bool (ROL::BoundConstraint<double>::*)() const) &ROL::BoundConstraint<double>::isActivated, "C++: ROL::BoundConstraint<double>::isActivated() const --> bool");
		cl.def("pruneActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) -> void { return o.pruneActive(a0, a1); }, "", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("pruneActive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double)) &ROL::BoundConstraint<double>::pruneActive, "C++: ROL::BoundConstraint<double>::pruneActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, double) --> void", pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("eps"));
		cl.def("pruneActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> void { return o.pruneActive(a0, a1, a2); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"));
		cl.def("pruneActive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double const & a3) -> void { return o.pruneActive(a0, a1, a2, a3); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"));
		cl.def("pruneActive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double)) &ROL::BoundConstraint<double>::pruneActive, "C++: ROL::BoundConstraint<double>::pruneActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double) --> void", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"), pybind11::arg("geps"));
		cl.def("pruneLowerInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) -> void { return o.pruneLowerInactive(a0, a1); }, "", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("pruneLowerInactive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double)) &ROL::BoundConstraint<double>::pruneLowerInactive, "C++: ROL::BoundConstraint<double>::pruneLowerInactive(class ROL::Vector<double> &, const class ROL::Vector<double> &, double) --> void", pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("eps"));
		cl.def("pruneUpperInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) -> void { return o.pruneUpperInactive(a0, a1); }, "", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("pruneUpperInactive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double)) &ROL::BoundConstraint<double>::pruneUpperInactive, "C++: ROL::BoundConstraint<double>::pruneUpperInactive(class ROL::Vector<double> &, const class ROL::Vector<double> &, double) --> void", pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("eps"));
		cl.def("pruneLowerInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> void { return o.pruneLowerInactive(a0, a1, a2); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"));
		cl.def("pruneLowerInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double const & a3) -> void { return o.pruneLowerInactive(a0, a1, a2, a3); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"));
		cl.def("pruneLowerInactive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double)) &ROL::BoundConstraint<double>::pruneLowerInactive, "C++: ROL::BoundConstraint<double>::pruneLowerInactive(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double) --> void", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"), pybind11::arg("geps"));
		cl.def("pruneUpperInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> void { return o.pruneUpperInactive(a0, a1, a2); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"));
		cl.def("pruneUpperInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double const & a3) -> void { return o.pruneUpperInactive(a0, a1, a2, a3); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"));
		cl.def("pruneUpperInactive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double)) &ROL::BoundConstraint<double>::pruneUpperInactive, "C++: ROL::BoundConstraint<double>::pruneUpperInactive(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double) --> void", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"), pybind11::arg("geps"));
		cl.def("pruneInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) -> void { return o.pruneInactive(a0, a1); }, "", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("pruneInactive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double)) &ROL::BoundConstraint<double>::pruneInactive, "C++: ROL::BoundConstraint<double>::pruneInactive(class ROL::Vector<double> &, const class ROL::Vector<double> &, double) --> void", pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("eps"));
		cl.def("pruneInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> void { return o.pruneInactive(a0, a1, a2); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"));
		cl.def("pruneInactive", [](ROL::BoundConstraint<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double const & a3) -> void { return o.pruneInactive(a0, a1, a2, a3); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"));
		cl.def("pruneInactive", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double)) &ROL::BoundConstraint<double>::pruneInactive, "C++: ROL::BoundConstraint<double>::pruneInactive(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double) --> void", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"), pybind11::arg("geps"));
		cl.def("computeProjectedGradient", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::BoundConstraint<double>::computeProjectedGradient, "C++: ROL::BoundConstraint<double>::computeProjectedGradient(class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("g"), pybind11::arg("x"));
		cl.def("computeProjectedStep", (void (ROL::BoundConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::BoundConstraint<double>::computeProjectedStep, "C++: ROL::BoundConstraint<double>::computeProjectedStep(class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("assign", (class ROL::BoundConstraint<double> & (ROL::BoundConstraint<double>::*)(const class ROL::BoundConstraint<double> &)) &ROL::BoundConstraint<double>::operator=, "C++: ROL::BoundConstraint<double>::operator=(const class ROL::BoundConstraint<double> &) --> class ROL::BoundConstraint<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
}
