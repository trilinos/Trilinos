#include <ROL_BoundConstraint.hpp>
#include <ROL_Constraint.hpp>
#include <ROL_DescentDirection_U.hpp>
#include <ROL_Elementwise_Function.hpp>
#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_Gradient_U.hpp>
#include <ROL_Krylov.hpp>
#include <ROL_NewtonKrylov_U.hpp>
#include <ROL_Newton_U.hpp>
#include <ROL_NonlinearCG.hpp>
#include <ROL_NonlinearCG_U.hpp>
#include <ROL_Objective.hpp>
#include <ROL_QuasiNewton_U.hpp>
#include <ROL_Secant.hpp>
#include <ROL_Step.hpp>
#include <ROL_TrustRegionModel_U.hpp>
#include <ROL_TrustRegion_U.hpp>
#include <ROL_TrustRegion_U_Types.hpp>
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
#include <Teuchos_any.hpp>
#include <cwchar>
#include <deque>
#include <ios>
#include <iterator>
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

// ROL::Step file:ROL_Step.hpp line:68
struct PyCallBack_ROL_Step_double_t : public ROL::Step<double> {
	using ROL::Step<double>::Step;

	void initialize(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, struct ROL::AlgorithmState<double> & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Step<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Step::initialize(a0, a1, a2, a3, a4);
	}
	void initialize(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, class ROL::Objective<double> & a3, class ROL::BoundConstraint<double> & a4, struct ROL::AlgorithmState<double> & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Step<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Step::initialize(a0, a1, a2, a3, a4, a5);
	}
	void initialize(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, class ROL::Objective<double> & a4, class ROL::Constraint<double> & a5, struct ROL::AlgorithmState<double> & a6) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Step<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Step::initialize(a0, a1, a2, a3, a4, a5, a6);
	}
	void initialize(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, class ROL::Objective<double> & a4, class ROL::Constraint<double> & a5, class ROL::BoundConstraint<double> & a6, struct ROL::AlgorithmState<double> & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Step<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Step::initialize(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	void compute(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, struct ROL::AlgorithmState<double> & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Step<double> *>(this), "compute");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Step::compute(a0, a1, a2, a3, a4);
	}
	void update(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, struct ROL::AlgorithmState<double> & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Step<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Step::update(a0, a1, a2, a3, a4);
	}
	void compute(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, class ROL::Objective<double> & a3, class ROL::Constraint<double> & a4, struct ROL::AlgorithmState<double> & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Step<double> *>(this), "compute");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Step::compute(a0, a1, a2, a3, a4, a5);
	}
	void update(class ROL::Vector<double> & a0, class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, class ROL::Objective<double> & a3, class ROL::Constraint<double> & a4, struct ROL::AlgorithmState<double> & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Step<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Step::update(a0, a1, a2, a3, a4, a5);
	}
	void compute(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, class ROL::Objective<double> & a3, class ROL::Constraint<double> & a4, class ROL::BoundConstraint<double> & a5, struct ROL::AlgorithmState<double> & a6) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Step<double> *>(this), "compute");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Step::compute(a0, a1, a2, a3, a4, a5, a6);
	}
	void update(class ROL::Vector<double> & a0, class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, class ROL::Objective<double> & a3, class ROL::Constraint<double> & a4, class ROL::BoundConstraint<double> & a5, struct ROL::AlgorithmState<double> & a6) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Step<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Step::update(a0, a1, a2, a3, a4, a5, a6);
	}
	std::string printHeader() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Step<double> *>(this), "printHeader");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<std::string>::value) {
				static pybind11::detail::override_caster_t<std::string> caster;
				return pybind11::detail::cast_ref<std::string>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<std::string>(std::move(o));
		}
		return Step::printHeader();
	}
	std::string printName() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Step<double> *>(this), "printName");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<std::string>::value) {
				static pybind11::detail::override_caster_t<std::string> caster;
				return pybind11::detail::cast_ref<std::string>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<std::string>(std::move(o));
		}
		return Step::printName();
	}
	std::string print(struct ROL::AlgorithmState<double> & a0, bool a1) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Step<double> *>(this), "print");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<std::string>::value) {
				static pybind11::detail::override_caster_t<std::string> caster;
				return pybind11::detail::cast_ref<std::string>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<std::string>(std::move(o));
		}
		return Step::print(a0, a1);
	}
};

// ROL::Gradient_U file:ROL_Gradient_U.hpp line:60
struct PyCallBack_ROL_Gradient_U_double_t : public ROL::Gradient_U<double> {
	using ROL::Gradient_U<double>::Gradient_U;

	void compute(class ROL::Vector<double> & a0, double & a1, double & a2, int & a3, int & a4, const class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Objective<double> & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Gradient_U<double> *>(this), "compute");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Gradient_U::compute(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	std::string printName() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Gradient_U<double> *>(this), "printName");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<std::string>::value) {
				static pybind11::detail::override_caster_t<std::string> caster;
				return pybind11::detail::cast_ref<std::string>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<std::string>(std::move(o));
		}
		return Gradient_U::printName();
	}
	void initialize(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Gradient_U<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return DescentDirection_U::initialize(a0, a1);
	}
	void update(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const double a4, const int a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Gradient_U<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return DescentDirection_U::update(a0, a1, a2, a3, a4, a5);
	}
};

// ROL::QuasiNewton_U file:ROL_QuasiNewton_U.hpp line:61
struct PyCallBack_ROL_QuasiNewton_U_double_t : public ROL::QuasiNewton_U<double> {
	using ROL::QuasiNewton_U<double>::QuasiNewton_U;

	void compute(class ROL::Vector<double> & a0, double & a1, double & a2, int & a3, int & a4, const class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Objective<double> & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::QuasiNewton_U<double> *>(this), "compute");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return QuasiNewton_U::compute(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	void update(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const double a4, const int a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::QuasiNewton_U<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return QuasiNewton_U::update(a0, a1, a2, a3, a4, a5);
	}
	std::string printName() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::QuasiNewton_U<double> *>(this), "printName");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<std::string>::value) {
				static pybind11::detail::override_caster_t<std::string> caster;
				return pybind11::detail::cast_ref<std::string>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<std::string>(std::move(o));
		}
		return QuasiNewton_U::printName();
	}
	void initialize(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::QuasiNewton_U<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return DescentDirection_U::initialize(a0, a1);
	}
};

// ROL::NonlinearCG file:ROL_NonlinearCG.hpp line:86
struct PyCallBack_ROL_NonlinearCG_double_t : public ROL::NonlinearCG<double> {
	using ROL::NonlinearCG<double>::NonlinearCG;

	void run(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, class ROL::Objective<double> & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::NonlinearCG<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NonlinearCG::run(a0, a1, a2, a3);
	}
};

// ROL::NonlinearCG_U file:ROL_NonlinearCG_U.hpp line:62
struct PyCallBack_ROL_NonlinearCG_U_double_t : public ROL::NonlinearCG_U<double> {
	using ROL::NonlinearCG_U<double>::NonlinearCG_U;

	void compute(class ROL::Vector<double> & a0, double & a1, double & a2, int & a3, int & a4, const class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Objective<double> & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::NonlinearCG_U<double> *>(this), "compute");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NonlinearCG_U::compute(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	std::string printName() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::NonlinearCG_U<double> *>(this), "printName");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<std::string>::value) {
				static pybind11::detail::override_caster_t<std::string> caster;
				return pybind11::detail::cast_ref<std::string>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<std::string>(std::move(o));
		}
		return NonlinearCG_U::printName();
	}
	void initialize(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::NonlinearCG_U<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return DescentDirection_U::initialize(a0, a1);
	}
	void update(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const double a4, const int a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::NonlinearCG_U<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return DescentDirection_U::update(a0, a1, a2, a3, a4, a5);
	}
};

// ROL::Newton_U file:ROL_Newton_U.hpp line:59
struct PyCallBack_ROL_Newton_U_double_t : public ROL::Newton_U<double> {
	using ROL::Newton_U<double>::Newton_U;

	void compute(class ROL::Vector<double> & a0, double & a1, double & a2, int & a3, int & a4, const class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Objective<double> & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Newton_U<double> *>(this), "compute");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Newton_U::compute(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	std::string printName() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Newton_U<double> *>(this), "printName");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<std::string>::value) {
				static pybind11::detail::override_caster_t<std::string> caster;
				return pybind11::detail::cast_ref<std::string>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<std::string>(std::move(o));
		}
		return Newton_U::printName();
	}
	void initialize(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Newton_U<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return DescentDirection_U::initialize(a0, a1);
	}
	void update(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const double a4, const int a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Newton_U<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return DescentDirection_U::update(a0, a1, a2, a3, a4, a5);
	}
};

// ROL::NewtonKrylov_U file:ROL_NewtonKrylov_U.hpp line:63
struct PyCallBack_ROL_NewtonKrylov_U_double_t : public ROL::NewtonKrylov_U<double> {
	using ROL::NewtonKrylov_U<double>::NewtonKrylov_U;

	void compute(class ROL::Vector<double> & a0, double & a1, double & a2, int & a3, int & a4, const class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Objective<double> & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::NewtonKrylov_U<double> *>(this), "compute");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NewtonKrylov_U::compute(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	void update(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const double a4, const int a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::NewtonKrylov_U<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NewtonKrylov_U::update(a0, a1, a2, a3, a4, a5);
	}
	std::string printName() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::NewtonKrylov_U<double> *>(this), "printName");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<std::string>::value) {
				static pybind11::detail::override_caster_t<std::string> caster;
				return pybind11::detail::cast_ref<std::string>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<std::string>(std::move(o));
		}
		return NewtonKrylov_U::printName();
	}
	void initialize(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::NewtonKrylov_U<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return DescentDirection_U::initialize(a0, a1);
	}
};

// ROL::TrustRegionModel_U file:ROL_TrustRegionModel_U.hpp line:66
struct PyCallBack_ROL_TrustRegionModel_U_double_t : public ROL::TrustRegionModel_U<double> {
	using ROL::TrustRegionModel_U<double>::TrustRegionModel_U;

	void setData(class ROL::Objective<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionModel_U<double> *>(this), "setData");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TrustRegionModel_U::setData(a0, a1, a2);
	}
	double value(const class ROL::Vector<double> & a0, double & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionModel_U<double> *>(this), "value");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return TrustRegionModel_U::value(a0, a1);
	}
	void gradient(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double & a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionModel_U<double> *>(this), "gradient");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TrustRegionModel_U::gradient(a0, a1, a2);
	}
	void hessVec(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionModel_U<double> *>(this), "hessVec");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TrustRegionModel_U::hessVec(a0, a1, a2, a3);
	}
	void invHessVec(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionModel_U<double> *>(this), "invHessVec");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TrustRegionModel_U::invHessVec(a0, a1, a2, a3);
	}
	void precond(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionModel_U<double> *>(this), "precond");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TrustRegionModel_U::precond(a0, a1, a2, a3);
	}
	const class Teuchos::RCP<const class ROL::Vector<double> > getGradient() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionModel_U<double> *>(this), "getGradient");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const class Teuchos::RCP<const class ROL::Vector<double> >>::value) {
				static pybind11::detail::override_caster_t<const class Teuchos::RCP<const class ROL::Vector<double> >> caster;
				return pybind11::detail::cast_ref<const class Teuchos::RCP<const class ROL::Vector<double> >>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const class Teuchos::RCP<const class ROL::Vector<double> >>(std::move(o));
		}
		return TrustRegionModel_U::getGradient();
	}
	const class Teuchos::RCP<const class ROL::Vector<double> > getIterate() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionModel_U<double> *>(this), "getIterate");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const class Teuchos::RCP<const class ROL::Vector<double> >>::value) {
				static pybind11::detail::override_caster_t<const class Teuchos::RCP<const class ROL::Vector<double> >> caster;
				return pybind11::detail::cast_ref<const class Teuchos::RCP<const class ROL::Vector<double> >>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const class Teuchos::RCP<const class ROL::Vector<double> >>(std::move(o));
		}
		return TrustRegionModel_U::getIterate();
	}
	const class Teuchos::RCP<class ROL::Objective<double> > getObjective() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionModel_U<double> *>(this), "getObjective");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const class Teuchos::RCP<class ROL::Objective<double> >>::value) {
				static pybind11::detail::override_caster_t<const class Teuchos::RCP<class ROL::Objective<double> >> caster;
				return pybind11::detail::cast_ref<const class Teuchos::RCP<class ROL::Objective<double> >>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const class Teuchos::RCP<class ROL::Objective<double> >>(std::move(o));
		}
		return TrustRegionModel_U::getObjective();
	}
	void update(const class ROL::Vector<double> & a0, enum ROL::UpdateType a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionModel_U<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Objective::update(a0, a1, a2);
	}
	void update(const class ROL::Vector<double> & a0, bool a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionModel_U<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Objective::update(a0, a1, a2);
	}
	double dirDeriv(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double & a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionModel_U<double> *>(this), "dirDeriv");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return Objective::dirDeriv(a0, a1, a2);
	}
	void prox(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionModel_U<double> *>(this), "prox");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Objective::prox(a0, a1, a2, a3);
	}
};

// ROL::TrustRegion_U file:ROL_TrustRegion_U.hpp line:57
struct PyCallBack_ROL_TrustRegion_U_double_t : public ROL::TrustRegion_U<double> {
	using ROL::TrustRegion_U<double>::TrustRegion_U;

	void initialize(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegion_U<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TrustRegion_U::initialize(a0, a1);
	}
	void solve(class ROL::Vector<double> & a0, double & a1, double & a2, int & a3, int & a4, const double a5, class ROL::TrustRegionModel_U<double> & a6) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegion_U<double> *>(this), "solve");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"TrustRegion_U::solve\"");
	}
};

void bind_ROL_Step(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // ROL::Step file:ROL_Step.hpp line:68
		pybind11::class_<ROL::Step<double>, Teuchos::RCP<ROL::Step<double>>, PyCallBack_ROL_Step_double_t> cl(M("ROL"), "Step_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new ROL::Step<double>(); }, [](){ return new PyCallBack_ROL_Step_double_t(); } ) );
		cl.def( pybind11::init( [](PyCallBack_ROL_Step_double_t const &o){ return new PyCallBack_ROL_Step_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::Step<double> const &o){ return new ROL::Step<double>(o); } ) );
		cl.def("initialize", (void (ROL::Step<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &)) &ROL::Step<double>::initialize, "C++: ROL::Step<double>::initialize(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("con"), pybind11::arg("algo_state"));
		cl.def("initialize", (void (ROL::Step<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &)) &ROL::Step<double>::initialize, "C++: ROL::Step<double>::initialize(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &) --> void", pybind11::arg("x"), pybind11::arg("s"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("con"), pybind11::arg("algo_state"));
		cl.def("initialize", (void (ROL::Step<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, struct ROL::AlgorithmState<double> &)) &ROL::Step<double>::initialize, "C++: ROL::Step<double>::initialize(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, struct ROL::AlgorithmState<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("l"), pybind11::arg("c"), pybind11::arg("obj"), pybind11::arg("con"), pybind11::arg("algo_state"));
		cl.def("initialize", (void (ROL::Step<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &)) &ROL::Step<double>::initialize, "C++: ROL::Step<double>::initialize(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("l"), pybind11::arg("c"), pybind11::arg("obj"), pybind11::arg("con"), pybind11::arg("bnd"), pybind11::arg("algo_state"));
		cl.def("compute", (void (ROL::Step<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &)) &ROL::Step<double>::compute, "C++: ROL::Step<double>::compute(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &) --> void", pybind11::arg("s"), pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("algo_state"));
		cl.def("update", (void (ROL::Step<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &)) &ROL::Step<double>::update, "C++: ROL::Step<double>::update(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &) --> void", pybind11::arg("x"), pybind11::arg("s"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("algo_state"));
		cl.def("compute", (void (ROL::Step<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, struct ROL::AlgorithmState<double> &)) &ROL::Step<double>::compute, "C++: ROL::Step<double>::compute(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, struct ROL::AlgorithmState<double> &) --> void", pybind11::arg("s"), pybind11::arg("x"), pybind11::arg("l"), pybind11::arg("obj"), pybind11::arg("con"), pybind11::arg("algo_state"));
		cl.def("update", (void (ROL::Step<double>::*)(class ROL::Vector<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, struct ROL::AlgorithmState<double> &)) &ROL::Step<double>::update, "C++: ROL::Step<double>::update(class ROL::Vector<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, struct ROL::AlgorithmState<double> &) --> void", pybind11::arg("x"), pybind11::arg("l"), pybind11::arg("s"), pybind11::arg("obj"), pybind11::arg("con"), pybind11::arg("algo_state"));
		cl.def("compute", (void (ROL::Step<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &)) &ROL::Step<double>::compute, "C++: ROL::Step<double>::compute(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &) --> void", pybind11::arg("s"), pybind11::arg("x"), pybind11::arg("l"), pybind11::arg("obj"), pybind11::arg("con"), pybind11::arg("bnd"), pybind11::arg("algo_state"));
		cl.def("update", (void (ROL::Step<double>::*)(class ROL::Vector<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &)) &ROL::Step<double>::update, "C++: ROL::Step<double>::update(class ROL::Vector<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::Constraint<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &) --> void", pybind11::arg("x"), pybind11::arg("l"), pybind11::arg("s"), pybind11::arg("obj"), pybind11::arg("con"), pybind11::arg("bnd"), pybind11::arg("algo_state"));
		cl.def("printHeader", (std::string (ROL::Step<double>::*)() const) &ROL::Step<double>::printHeader, "C++: ROL::Step<double>::printHeader() const --> std::string");
		cl.def("printName", (std::string (ROL::Step<double>::*)() const) &ROL::Step<double>::printName, "C++: ROL::Step<double>::printName() const --> std::string");
		cl.def("print", [](ROL::Step<double> const &o, struct ROL::AlgorithmState<double> & a0) -> std::string { return o.print(a0); }, "", pybind11::arg("algo_state"));
		cl.def("print", (std::string (ROL::Step<double>::*)(struct ROL::AlgorithmState<double> &, bool) const) &ROL::Step<double>::print, "C++: ROL::Step<double>::print(struct ROL::AlgorithmState<double> &, bool) const --> std::string", pybind11::arg("algo_state"), pybind11::arg("printHeader"));
		cl.def("getStepState", (const class Teuchos::RCP<const struct ROL::StepState<double> > (ROL::Step<double>::*)() const) &ROL::Step<double>::getStepState, "C++: ROL::Step<double>::getStepState() const --> const class Teuchos::RCP<const struct ROL::StepState<double> >");
		cl.def("reset", [](ROL::Step<double> &o) -> void { return o.reset(); }, "");
		cl.def("reset", (void (ROL::Step<double>::*)(const double)) &ROL::Step<double>::reset, "C++: ROL::Step<double>::reset(const double) --> void", pybind11::arg("searchSize"));
		cl.def("assign", (class ROL::Step<double> & (ROL::Step<double>::*)(const class ROL::Step<double> &)) &ROL::Step<double>::operator=, "C++: ROL::Step<double>::operator=(const class ROL::Step<double> &) --> class ROL::Step<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::Gradient_U file:ROL_Gradient_U.hpp line:60
		pybind11::class_<ROL::Gradient_U<double>, Teuchos::RCP<ROL::Gradient_U<double>>, PyCallBack_ROL_Gradient_U_double_t, ROL::DescentDirection_U<double>> cl(M("ROL"), "Gradient_U_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new ROL::Gradient_U<double>(); }, [](){ return new PyCallBack_ROL_Gradient_U_double_t(); } ) );
		cl.def( pybind11::init( [](PyCallBack_ROL_Gradient_U_double_t const &o){ return new PyCallBack_ROL_Gradient_U_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::Gradient_U<double> const &o){ return new ROL::Gradient_U<double>(o); } ) );
		cl.def("compute", (void (ROL::Gradient_U<double>::*)(class ROL::Vector<double> &, double &, double &, int &, int &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &)) &ROL::Gradient_U<double>::compute, "C++: ROL::Gradient_U<double>::compute(class ROL::Vector<double> &, double &, double &, int &, int &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &) --> void", pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("sdotg"), pybind11::arg("iter"), pybind11::arg("flag"), pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"));
		cl.def("printName", (std::string (ROL::Gradient_U<double>::*)() const) &ROL::Gradient_U<double>::printName, "C++: ROL::Gradient_U<double>::printName() const --> std::string");
		cl.def("assign", (class ROL::Gradient_U<double> & (ROL::Gradient_U<double>::*)(const class ROL::Gradient_U<double> &)) &ROL::Gradient_U<double>::operator=, "C++: ROL::Gradient_U<double>::operator=(const class ROL::Gradient_U<double> &) --> class ROL::Gradient_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("initialize", (void (ROL::DescentDirection_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::DescentDirection_U<double>::initialize, "C++: ROL::DescentDirection_U<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"));
		cl.def("compute", (void (ROL::DescentDirection_U<double>::*)(class ROL::Vector<double> &, double &, double &, int &, int &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &)) &ROL::DescentDirection_U<double>::compute, "C++: ROL::DescentDirection_U<double>::compute(class ROL::Vector<double> &, double &, double &, int &, int &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &) --> void", pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("sdotg"), pybind11::arg("iter"), pybind11::arg("flag"), pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"));
		cl.def("update", (void (ROL::DescentDirection_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const double, const int)) &ROL::DescentDirection_U<double>::update, "C++: ROL::DescentDirection_U<double>::update(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const double, const int) --> void", pybind11::arg("x"), pybind11::arg("s"), pybind11::arg("gold"), pybind11::arg("gnew"), pybind11::arg("snorm"), pybind11::arg("iter"));
		cl.def("printName", (std::string (ROL::DescentDirection_U<double>::*)() const) &ROL::DescentDirection_U<double>::printName, "C++: ROL::DescentDirection_U<double>::printName() const --> std::string");
		cl.def("assign", (class ROL::DescentDirection_U<double> & (ROL::DescentDirection_U<double>::*)(const class ROL::DescentDirection_U<double> &)) &ROL::DescentDirection_U<double>::operator=, "C++: ROL::DescentDirection_U<double>::operator=(const class ROL::DescentDirection_U<double> &) --> class ROL::DescentDirection_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::QuasiNewton_U file:ROL_QuasiNewton_U.hpp line:61
		pybind11::class_<ROL::QuasiNewton_U<double>, Teuchos::RCP<ROL::QuasiNewton_U<double>>, PyCallBack_ROL_QuasiNewton_U_double_t, ROL::DescentDirection_U<double>> cl(M("ROL"), "QuasiNewton_U_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](class Teuchos::ParameterList & a0){ return new ROL::QuasiNewton_U<double>(a0); }, [](class Teuchos::ParameterList & a0){ return new PyCallBack_ROL_QuasiNewton_U_double_t(a0); } ), "doc");
		cl.def( pybind11::init<class Teuchos::ParameterList &, const class Teuchos::RCP<class ROL::Secant<double> > &>(), pybind11::arg("parlist"), pybind11::arg("secant") );

		cl.def( pybind11::init( [](PyCallBack_ROL_QuasiNewton_U_double_t const &o){ return new PyCallBack_ROL_QuasiNewton_U_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::QuasiNewton_U<double> const &o){ return new ROL::QuasiNewton_U<double>(o); } ) );
		cl.def("compute", (void (ROL::QuasiNewton_U<double>::*)(class ROL::Vector<double> &, double &, double &, int &, int &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &)) &ROL::QuasiNewton_U<double>::compute, "C++: ROL::QuasiNewton_U<double>::compute(class ROL::Vector<double> &, double &, double &, int &, int &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &) --> void", pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("sdotg"), pybind11::arg("iter"), pybind11::arg("flag"), pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"));
		cl.def("update", (void (ROL::QuasiNewton_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const double, const int)) &ROL::QuasiNewton_U<double>::update, "C++: ROL::QuasiNewton_U<double>::update(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const double, const int) --> void", pybind11::arg("x"), pybind11::arg("s"), pybind11::arg("gold"), pybind11::arg("gnew"), pybind11::arg("snorm"), pybind11::arg("iter"));
		cl.def("printName", (std::string (ROL::QuasiNewton_U<double>::*)() const) &ROL::QuasiNewton_U<double>::printName, "C++: ROL::QuasiNewton_U<double>::printName() const --> std::string");
		cl.def("assign", (class ROL::QuasiNewton_U<double> & (ROL::QuasiNewton_U<double>::*)(const class ROL::QuasiNewton_U<double> &)) &ROL::QuasiNewton_U<double>::operator=, "C++: ROL::QuasiNewton_U<double>::operator=(const class ROL::QuasiNewton_U<double> &) --> class ROL::QuasiNewton_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("initialize", (void (ROL::DescentDirection_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::DescentDirection_U<double>::initialize, "C++: ROL::DescentDirection_U<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"));
		cl.def("compute", (void (ROL::DescentDirection_U<double>::*)(class ROL::Vector<double> &, double &, double &, int &, int &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &)) &ROL::DescentDirection_U<double>::compute, "C++: ROL::DescentDirection_U<double>::compute(class ROL::Vector<double> &, double &, double &, int &, int &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &) --> void", pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("sdotg"), pybind11::arg("iter"), pybind11::arg("flag"), pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"));
		cl.def("update", (void (ROL::DescentDirection_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const double, const int)) &ROL::DescentDirection_U<double>::update, "C++: ROL::DescentDirection_U<double>::update(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const double, const int) --> void", pybind11::arg("x"), pybind11::arg("s"), pybind11::arg("gold"), pybind11::arg("gnew"), pybind11::arg("snorm"), pybind11::arg("iter"));
		cl.def("printName", (std::string (ROL::DescentDirection_U<double>::*)() const) &ROL::DescentDirection_U<double>::printName, "C++: ROL::DescentDirection_U<double>::printName() const --> std::string");
		cl.def("assign", (class ROL::DescentDirection_U<double> & (ROL::DescentDirection_U<double>::*)(const class ROL::DescentDirection_U<double> &)) &ROL::DescentDirection_U<double>::operator=, "C++: ROL::DescentDirection_U<double>::operator=(const class ROL::DescentDirection_U<double> &) --> class ROL::DescentDirection_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::NonlinearCGState file:ROL_NonlinearCG.hpp line:77
		pybind11::class_<ROL::NonlinearCGState<double>, Teuchos::RCP<ROL::NonlinearCGState<double>>> cl(M("ROL"), "NonlinearCGState_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](ROL::NonlinearCGState<double> const &o){ return new ROL::NonlinearCGState<double>(o); } ) );
		cl.def( pybind11::init( [](){ return new ROL::NonlinearCGState<double>(); } ) );
		cl.def_readwrite("grad", &ROL::NonlinearCGState<double>::grad);
		cl.def_readwrite("pstep", &ROL::NonlinearCGState<double>::pstep);
		cl.def_readwrite("iter", &ROL::NonlinearCGState<double>::iter);
		cl.def_readwrite("restart", &ROL::NonlinearCGState<double>::restart);
		cl.def_readwrite("nlcg_type", &ROL::NonlinearCGState<double>::nlcg_type);
	}
	{ // ROL::NonlinearCG file:ROL_NonlinearCG.hpp line:86
		pybind11::class_<ROL::NonlinearCG<double>, Teuchos::RCP<ROL::NonlinearCG<double>>, PyCallBack_ROL_NonlinearCG_double_t> cl(M("ROL"), "NonlinearCG_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](enum ROL::ENonlinearCG const & a0){ return new ROL::NonlinearCG<double>(a0); }, [](enum ROL::ENonlinearCG const & a0){ return new PyCallBack_ROL_NonlinearCG_double_t(a0); } ), "doc");
		cl.def( pybind11::init<enum ROL::ENonlinearCG, int>(), pybind11::arg("type"), pybind11::arg("restart") );

		cl.def( pybind11::init( [](PyCallBack_ROL_NonlinearCG_double_t const &o){ return new PyCallBack_ROL_NonlinearCG_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::NonlinearCG<double> const &o){ return new ROL::NonlinearCG<double>(o); } ) );
		cl.def("run", (void (ROL::NonlinearCG<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &)) &ROL::NonlinearCG<double>::run, "C++: ROL::NonlinearCG<double>::run(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &) --> void", pybind11::arg("s"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("obj"));
		cl.def("assign", (class ROL::NonlinearCG<double> & (ROL::NonlinearCG<double>::*)(const class ROL::NonlinearCG<double> &)) &ROL::NonlinearCG<double>::operator=, "C++: ROL::NonlinearCG<double>::operator=(const class ROL::NonlinearCG<double> &) --> class ROL::NonlinearCG<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::NonlinearCG_U file:ROL_NonlinearCG_U.hpp line:62
		pybind11::class_<ROL::NonlinearCG_U<double>, Teuchos::RCP<ROL::NonlinearCG_U<double>>, PyCallBack_ROL_NonlinearCG_U_double_t, ROL::DescentDirection_U<double>> cl(M("ROL"), "NonlinearCG_U_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](class Teuchos::ParameterList & a0){ return new ROL::NonlinearCG_U<double>(a0); }, [](class Teuchos::ParameterList & a0){ return new PyCallBack_ROL_NonlinearCG_U_double_t(a0); } ), "doc");
		cl.def( pybind11::init<class Teuchos::ParameterList &, const class Teuchos::RCP<class ROL::NonlinearCG<double> > &>(), pybind11::arg("parlist"), pybind11::arg("nlcg") );

		cl.def( pybind11::init( [](PyCallBack_ROL_NonlinearCG_U_double_t const &o){ return new PyCallBack_ROL_NonlinearCG_U_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::NonlinearCG_U<double> const &o){ return new ROL::NonlinearCG_U<double>(o); } ) );
		cl.def("compute", (void (ROL::NonlinearCG_U<double>::*)(class ROL::Vector<double> &, double &, double &, int &, int &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &)) &ROL::NonlinearCG_U<double>::compute, "C++: ROL::NonlinearCG_U<double>::compute(class ROL::Vector<double> &, double &, double &, int &, int &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &) --> void", pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("sdotg"), pybind11::arg("iter"), pybind11::arg("flag"), pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"));
		cl.def("printName", (std::string (ROL::NonlinearCG_U<double>::*)() const) &ROL::NonlinearCG_U<double>::printName, "C++: ROL::NonlinearCG_U<double>::printName() const --> std::string");
		cl.def("assign", (class ROL::NonlinearCG_U<double> & (ROL::NonlinearCG_U<double>::*)(const class ROL::NonlinearCG_U<double> &)) &ROL::NonlinearCG_U<double>::operator=, "C++: ROL::NonlinearCG_U<double>::operator=(const class ROL::NonlinearCG_U<double> &) --> class ROL::NonlinearCG_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("initialize", (void (ROL::DescentDirection_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::DescentDirection_U<double>::initialize, "C++: ROL::DescentDirection_U<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"));
		cl.def("compute", (void (ROL::DescentDirection_U<double>::*)(class ROL::Vector<double> &, double &, double &, int &, int &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &)) &ROL::DescentDirection_U<double>::compute, "C++: ROL::DescentDirection_U<double>::compute(class ROL::Vector<double> &, double &, double &, int &, int &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &) --> void", pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("sdotg"), pybind11::arg("iter"), pybind11::arg("flag"), pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"));
		cl.def("update", (void (ROL::DescentDirection_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const double, const int)) &ROL::DescentDirection_U<double>::update, "C++: ROL::DescentDirection_U<double>::update(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const double, const int) --> void", pybind11::arg("x"), pybind11::arg("s"), pybind11::arg("gold"), pybind11::arg("gnew"), pybind11::arg("snorm"), pybind11::arg("iter"));
		cl.def("printName", (std::string (ROL::DescentDirection_U<double>::*)() const) &ROL::DescentDirection_U<double>::printName, "C++: ROL::DescentDirection_U<double>::printName() const --> std::string");
		cl.def("assign", (class ROL::DescentDirection_U<double> & (ROL::DescentDirection_U<double>::*)(const class ROL::DescentDirection_U<double> &)) &ROL::DescentDirection_U<double>::operator=, "C++: ROL::DescentDirection_U<double>::operator=(const class ROL::DescentDirection_U<double> &) --> class ROL::DescentDirection_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::Newton_U file:ROL_Newton_U.hpp line:59
		pybind11::class_<ROL::Newton_U<double>, Teuchos::RCP<ROL::Newton_U<double>>, PyCallBack_ROL_Newton_U_double_t, ROL::DescentDirection_U<double>> cl(M("ROL"), "Newton_U_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new ROL::Newton_U<double>(); }, [](){ return new PyCallBack_ROL_Newton_U_double_t(); } ) );
		cl.def( pybind11::init( [](PyCallBack_ROL_Newton_U_double_t const &o){ return new PyCallBack_ROL_Newton_U_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::Newton_U<double> const &o){ return new ROL::Newton_U<double>(o); } ) );
		cl.def("compute", (void (ROL::Newton_U<double>::*)(class ROL::Vector<double> &, double &, double &, int &, int &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &)) &ROL::Newton_U<double>::compute, "C++: ROL::Newton_U<double>::compute(class ROL::Vector<double> &, double &, double &, int &, int &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &) --> void", pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("sdotg"), pybind11::arg("iter"), pybind11::arg("flag"), pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"));
		cl.def("printName", (std::string (ROL::Newton_U<double>::*)() const) &ROL::Newton_U<double>::printName, "C++: ROL::Newton_U<double>::printName() const --> std::string");
		cl.def("assign", (class ROL::Newton_U<double> & (ROL::Newton_U<double>::*)(const class ROL::Newton_U<double> &)) &ROL::Newton_U<double>::operator=, "C++: ROL::Newton_U<double>::operator=(const class ROL::Newton_U<double> &) --> class ROL::Newton_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("initialize", (void (ROL::DescentDirection_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::DescentDirection_U<double>::initialize, "C++: ROL::DescentDirection_U<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"));
		cl.def("compute", (void (ROL::DescentDirection_U<double>::*)(class ROL::Vector<double> &, double &, double &, int &, int &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &)) &ROL::DescentDirection_U<double>::compute, "C++: ROL::DescentDirection_U<double>::compute(class ROL::Vector<double> &, double &, double &, int &, int &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &) --> void", pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("sdotg"), pybind11::arg("iter"), pybind11::arg("flag"), pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"));
		cl.def("update", (void (ROL::DescentDirection_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const double, const int)) &ROL::DescentDirection_U<double>::update, "C++: ROL::DescentDirection_U<double>::update(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const double, const int) --> void", pybind11::arg("x"), pybind11::arg("s"), pybind11::arg("gold"), pybind11::arg("gnew"), pybind11::arg("snorm"), pybind11::arg("iter"));
		cl.def("printName", (std::string (ROL::DescentDirection_U<double>::*)() const) &ROL::DescentDirection_U<double>::printName, "C++: ROL::DescentDirection_U<double>::printName() const --> std::string");
		cl.def("assign", (class ROL::DescentDirection_U<double> & (ROL::DescentDirection_U<double>::*)(const class ROL::DescentDirection_U<double> &)) &ROL::DescentDirection_U<double>::operator=, "C++: ROL::DescentDirection_U<double>::operator=(const class ROL::DescentDirection_U<double> &) --> class ROL::DescentDirection_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::NewtonKrylov_U file:ROL_NewtonKrylov_U.hpp line:63
		pybind11::class_<ROL::NewtonKrylov_U<double>, Teuchos::RCP<ROL::NewtonKrylov_U<double>>, PyCallBack_ROL_NewtonKrylov_U_double_t, ROL::DescentDirection_U<double>> cl(M("ROL"), "NewtonKrylov_U_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("parlist") );

		cl.def( pybind11::init( [](class Teuchos::ParameterList & a0, const class Teuchos::RCP<class ROL::Krylov<double> > & a1, const class Teuchos::RCP<class ROL::Secant<double> > & a2){ return new ROL::NewtonKrylov_U<double>(a0, a1, a2); }, [](class Teuchos::ParameterList & a0, const class Teuchos::RCP<class ROL::Krylov<double> > & a1, const class Teuchos::RCP<class ROL::Secant<double> > & a2){ return new PyCallBack_ROL_NewtonKrylov_U_double_t(a0, a1, a2); } ), "doc");
		cl.def( pybind11::init<class Teuchos::ParameterList &, const class Teuchos::RCP<class ROL::Krylov<double> > &, const class Teuchos::RCP<class ROL::Secant<double> > &, const bool>(), pybind11::arg("parlist"), pybind11::arg("krylov"), pybind11::arg("secant"), pybind11::arg("computeObj") );

		cl.def( pybind11::init( [](PyCallBack_ROL_NewtonKrylov_U_double_t const &o){ return new PyCallBack_ROL_NewtonKrylov_U_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::NewtonKrylov_U<double> const &o){ return new ROL::NewtonKrylov_U<double>(o); } ) );
		cl.def("compute", (void (ROL::NewtonKrylov_U<double>::*)(class ROL::Vector<double> &, double &, double &, int &, int &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &)) &ROL::NewtonKrylov_U<double>::compute, "C++: ROL::NewtonKrylov_U<double>::compute(class ROL::Vector<double> &, double &, double &, int &, int &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &) --> void", pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("sdotg"), pybind11::arg("iter"), pybind11::arg("flag"), pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"));
		cl.def("update", (void (ROL::NewtonKrylov_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const double, const int)) &ROL::NewtonKrylov_U<double>::update, "C++: ROL::NewtonKrylov_U<double>::update(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const double, const int) --> void", pybind11::arg("x"), pybind11::arg("s"), pybind11::arg("gold"), pybind11::arg("gnew"), pybind11::arg("snorm"), pybind11::arg("iter"));
		cl.def("printName", (std::string (ROL::NewtonKrylov_U<double>::*)() const) &ROL::NewtonKrylov_U<double>::printName, "C++: ROL::NewtonKrylov_U<double>::printName() const --> std::string");
		cl.def("assign", (class ROL::NewtonKrylov_U<double> & (ROL::NewtonKrylov_U<double>::*)(const class ROL::NewtonKrylov_U<double> &)) &ROL::NewtonKrylov_U<double>::operator=, "C++: ROL::NewtonKrylov_U<double>::operator=(const class ROL::NewtonKrylov_U<double> &) --> class ROL::NewtonKrylov_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("initialize", (void (ROL::DescentDirection_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::DescentDirection_U<double>::initialize, "C++: ROL::DescentDirection_U<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"));
		cl.def("compute", (void (ROL::DescentDirection_U<double>::*)(class ROL::Vector<double> &, double &, double &, int &, int &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &)) &ROL::DescentDirection_U<double>::compute, "C++: ROL::DescentDirection_U<double>::compute(class ROL::Vector<double> &, double &, double &, int &, int &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &) --> void", pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("sdotg"), pybind11::arg("iter"), pybind11::arg("flag"), pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"));
		cl.def("update", (void (ROL::DescentDirection_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const double, const int)) &ROL::DescentDirection_U<double>::update, "C++: ROL::DescentDirection_U<double>::update(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const double, const int) --> void", pybind11::arg("x"), pybind11::arg("s"), pybind11::arg("gold"), pybind11::arg("gnew"), pybind11::arg("snorm"), pybind11::arg("iter"));
		cl.def("printName", (std::string (ROL::DescentDirection_U<double>::*)() const) &ROL::DescentDirection_U<double>::printName, "C++: ROL::DescentDirection_U<double>::printName() const --> std::string");
		cl.def("assign", (class ROL::DescentDirection_U<double> & (ROL::DescentDirection_U<double>::*)(const class ROL::DescentDirection_U<double> &)) &ROL::DescentDirection_U<double>::operator=, "C++: ROL::DescentDirection_U<double>::operator=(const class ROL::DescentDirection_U<double> &) --> class ROL::DescentDirection_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	// ROL::ETrustRegionU file:ROL_TrustRegion_U_Types.hpp line:63
	pybind11::enum_<ROL::ETrustRegionU>(M("ROL"), "ETrustRegionU", pybind11::arithmetic(), "Enumeration of trust-region solver types.\n\n      \n    TRUSTREGION_U_CAUCHYPOINT     describe\n      \n\n    TRUSTREGION_U_TRUNCATEDCG     describe\n      \n\n    TRUSTREGION_U_SPG             describe\n      \n\n    TRUSTREGION_U_DOGLEG          describe\n      \n\n    TRUSTREGION_U_DOUBLEDOGLEG    describe", pybind11::module_local())
		.value("TRUSTREGION_U_CAUCHYPOINT", ROL::TRUSTREGION_U_CAUCHYPOINT)
		.value("TRUSTREGION_U_TRUNCATEDCG", ROL::TRUSTREGION_U_TRUNCATEDCG)
		.value("TRUSTREGION_U_SPG", ROL::TRUSTREGION_U_SPG)
		.value("TRUSTREGION_U_DOGLEG", ROL::TRUSTREGION_U_DOGLEG)
		.value("TRUSTREGION_U_DOUBLEDOGLEG", ROL::TRUSTREGION_U_DOUBLEDOGLEG)
		.value("TRUSTREGION_U_LAST", ROL::TRUSTREGION_U_LAST)
		.export_values();

;

	// ROL::ETrustRegionUToString(enum ROL::ETrustRegionU) file:ROL_TrustRegion_U_Types.hpp line:72
	M("ROL").def("ETrustRegionUToString", (std::string (*)(enum ROL::ETrustRegionU)) &ROL::ETrustRegionUToString, "C++: ROL::ETrustRegionUToString(enum ROL::ETrustRegionU) --> std::string", pybind11::arg("tr"));

	// ROL::isValidTrustRegionU(enum ROL::ETrustRegionU) file:ROL_TrustRegion_U_Types.hpp line:91
	M("ROL").def("isValidTrustRegionU", (int (*)(enum ROL::ETrustRegionU)) &ROL::isValidTrustRegionU, "Verifies validity of a TrustRegionU enum.\n\n      \n  [in]  - enum of the TrustRegionU\n      \n\n 1 if the argument is a valid TrustRegionU; 0 otherwise.\n\nC++: ROL::isValidTrustRegionU(enum ROL::ETrustRegionU) --> int", pybind11::arg("ls"));

	// ROL::StringToETrustRegionU(std::string) file:ROL_TrustRegion_U_Types.hpp line:120
	M("ROL").def("StringToETrustRegionU", (enum ROL::ETrustRegionU (*)(std::string)) &ROL::StringToETrustRegionU, "C++: ROL::StringToETrustRegionU(std::string) --> enum ROL::ETrustRegionU", pybind11::arg("s"));

	{ // ROL::TrustRegionModel_U file:ROL_TrustRegionModel_U.hpp line:66
		pybind11::class_<ROL::TrustRegionModel_U<double>, Teuchos::RCP<ROL::TrustRegionModel_U<double>>, PyCallBack_ROL_TrustRegionModel_U_double_t, ROL::Objective<double>> cl(M("ROL"), "TrustRegionModel_U_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](class Teuchos::ParameterList & a0){ return new ROL::TrustRegionModel_U<double>(a0); }, [](class Teuchos::ParameterList & a0){ return new PyCallBack_ROL_TrustRegionModel_U_double_t(a0); } ), "doc");
		cl.def( pybind11::init( [](class Teuchos::ParameterList & a0, const class Teuchos::RCP<class ROL::Secant<double> > & a1){ return new ROL::TrustRegionModel_U<double>(a0, a1); }, [](class Teuchos::ParameterList & a0, const class Teuchos::RCP<class ROL::Secant<double> > & a1){ return new PyCallBack_ROL_TrustRegionModel_U_double_t(a0, a1); } ), "doc");
		cl.def( pybind11::init<class Teuchos::ParameterList &, const class Teuchos::RCP<class ROL::Secant<double> > &, enum ROL::ESecantMode>(), pybind11::arg("list"), pybind11::arg("secant"), pybind11::arg("mode") );

		cl.def( pybind11::init( [](PyCallBack_ROL_TrustRegionModel_U_double_t const &o){ return new PyCallBack_ROL_TrustRegionModel_U_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::TrustRegionModel_U<double> const &o){ return new ROL::TrustRegionModel_U<double>(o); } ) );
		cl.def("update", [](ROL::TrustRegionModel_U<double> &o, const class ROL::Vector<double> & a0) -> void { return o.update(a0); }, "", pybind11::arg("x"));
		cl.def("update", [](ROL::TrustRegionModel_U<double> &o, const class ROL::Vector<double> & a0, bool const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("flag"));
		cl.def("update", [](ROL::TrustRegionModel_U<double> &o, const class ROL::Vector<double> & a0, bool const & a1, int const & a2) -> void { return o.update(a0, a1, a2); }, "", pybind11::arg("x"), pybind11::arg("flag"), pybind11::arg("iter"));
		cl.def("update", [](ROL::TrustRegionModel_U<double> &o, const class ROL::Vector<double> & a0, enum ROL::UpdateType const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("type"));
		cl.def("update", [](ROL::TrustRegionModel_U<double> &o, const class ROL::Vector<double> & a0, enum ROL::UpdateType const & a1, int const & a2) -> void { return o.update(a0, a1, a2); }, "", pybind11::arg("x"), pybind11::arg("type"), pybind11::arg("iter"));
		cl.def("initialize", (void (ROL::TrustRegionModel_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::TrustRegionModel_U<double>::initialize, "C++: ROL::TrustRegionModel_U<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"));
		cl.def("validate", (void (ROL::TrustRegionModel_U<double>::*)(class ROL::Objective<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, enum ROL::ETrustRegionU)) &ROL::TrustRegionModel_U<double>::validate, "C++: ROL::TrustRegionModel_U<double>::validate(class ROL::Objective<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, enum ROL::ETrustRegionU) --> void", pybind11::arg("obj"), pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("etr"));
		cl.def("setData", (void (ROL::TrustRegionModel_U<double>::*)(class ROL::Objective<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::TrustRegionModel_U<double>::setData, "C++: ROL::TrustRegionModel_U<double>::setData(class ROL::Objective<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("obj"), pybind11::arg("x"), pybind11::arg("g"));
		cl.def("update", (void (ROL::TrustRegionModel_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const double, const int)) &ROL::TrustRegionModel_U<double>::update, "C++: ROL::TrustRegionModel_U<double>::update(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const double, const int) --> void", pybind11::arg("x"), pybind11::arg("s"), pybind11::arg("gold"), pybind11::arg("gnew"), pybind11::arg("snorm"), pybind11::arg("iter"));
		cl.def("value", (double (ROL::TrustRegionModel_U<double>::*)(const class ROL::Vector<double> &, double &)) &ROL::TrustRegionModel_U<double>::value, "C++: ROL::TrustRegionModel_U<double>::value(const class ROL::Vector<double> &, double &) --> double", pybind11::arg("s"), pybind11::arg("tol"));
		cl.def("gradient", (void (ROL::TrustRegionModel_U<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::TrustRegionModel_U<double>::gradient, "C++: ROL::TrustRegionModel_U<double>::gradient(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("g"), pybind11::arg("s"), pybind11::arg("tol"));
		cl.def("hessVec", (void (ROL::TrustRegionModel_U<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::TrustRegionModel_U<double>::hessVec, "C++: ROL::TrustRegionModel_U<double>::hessVec(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("hv"), pybind11::arg("v"), pybind11::arg("s"), pybind11::arg("tol"));
		cl.def("invHessVec", (void (ROL::TrustRegionModel_U<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::TrustRegionModel_U<double>::invHessVec, "C++: ROL::TrustRegionModel_U<double>::invHessVec(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("hv"), pybind11::arg("v"), pybind11::arg("s"), pybind11::arg("tol"));
		cl.def("precond", (void (ROL::TrustRegionModel_U<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::TrustRegionModel_U<double>::precond, "C++: ROL::TrustRegionModel_U<double>::precond(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("Pv"), pybind11::arg("v"), pybind11::arg("s"), pybind11::arg("tol"));
		cl.def("getGradient", (const class Teuchos::RCP<const class ROL::Vector<double> > (ROL::TrustRegionModel_U<double>::*)() const) &ROL::TrustRegionModel_U<double>::getGradient, "C++: ROL::TrustRegionModel_U<double>::getGradient() const --> const class Teuchos::RCP<const class ROL::Vector<double> >");
		cl.def("getIterate", (const class Teuchos::RCP<const class ROL::Vector<double> > (ROL::TrustRegionModel_U<double>::*)() const) &ROL::TrustRegionModel_U<double>::getIterate, "C++: ROL::TrustRegionModel_U<double>::getIterate() const --> const class Teuchos::RCP<const class ROL::Vector<double> >");
		cl.def("getObjective", (const class Teuchos::RCP<class ROL::Objective<double> > (ROL::TrustRegionModel_U<double>::*)() const) &ROL::TrustRegionModel_U<double>::getObjective, "C++: ROL::TrustRegionModel_U<double>::getObjective() const --> const class Teuchos::RCP<class ROL::Objective<double> >");
		cl.def("assign", (class ROL::TrustRegionModel_U<double> & (ROL::TrustRegionModel_U<double>::*)(const class ROL::TrustRegionModel_U<double> &)) &ROL::TrustRegionModel_U<double>::operator=, "C++: ROL::TrustRegionModel_U<double>::operator=(const class ROL::TrustRegionModel_U<double> &) --> class ROL::TrustRegionModel_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
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
		cl.def("prox", (void (ROL::Objective<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double &)) &ROL::Objective<double>::prox, "C++: ROL::Objective<double>::prox(class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double &) --> void", pybind11::arg("Pv"), pybind11::arg("v"), pybind11::arg("t"), pybind11::arg("tol"));
		cl.def("assign", (class ROL::Objective<double> & (ROL::Objective<double>::*)(const class ROL::Objective<double> &)) &ROL::Objective<double>::operator=, "C++: ROL::Objective<double>::operator=(const class ROL::Objective<double> &) --> class ROL::Objective<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::TrustRegion_U file:ROL_TrustRegion_U.hpp line:57
		pybind11::class_<ROL::TrustRegion_U<double>, Teuchos::RCP<ROL::TrustRegion_U<double>>, PyCallBack_ROL_TrustRegion_U_double_t> cl(M("ROL"), "TrustRegion_U_double_t", "", pybind11::module_local());
		cl.def(pybind11::init<PyCallBack_ROL_TrustRegion_U_double_t const &>());
		cl.def( pybind11::init( [](){ return new PyCallBack_ROL_TrustRegion_U_double_t(); } ) );
		cl.def("initialize", (void (ROL::TrustRegion_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::TrustRegion_U<double>::initialize, "C++: ROL::TrustRegion_U<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"));
		cl.def("solve", (void (ROL::TrustRegion_U<double>::*)(class ROL::Vector<double> &, double &, double &, int &, int &, const double, class ROL::TrustRegionModel_U<double> &)) &ROL::TrustRegion_U<double>::solve, "C++: ROL::TrustRegion_U<double>::solve(class ROL::Vector<double> &, double &, double &, int &, int &, const double, class ROL::TrustRegionModel_U<double> &) --> void", pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("pRed"), pybind11::arg("iflag"), pybind11::arg("iter"), pybind11::arg("del"), pybind11::arg("model"));
		cl.def("assign", (class ROL::TrustRegion_U<double> & (ROL::TrustRegion_U<double>::*)(const class ROL::TrustRegion_U<double> &)) &ROL::TrustRegion_U<double>::operator=, "C++: ROL::TrustRegion_U<double>::operator=(const class ROL::TrustRegion_U<double> &) --> class ROL::TrustRegion_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
}
