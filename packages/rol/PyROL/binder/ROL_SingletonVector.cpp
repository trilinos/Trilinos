#include <ROL_BoundConstraint.hpp>
#include <ROL_Elementwise_Function.hpp>
#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_MoreauYosidaObjective.hpp>
#include <ROL_Objective.hpp>
#include <ROL_ScalarController.hpp>
#include <ROL_SingletonVector.hpp>
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

// ROL::SingletonVector file:ROL_SingletonVector.hpp line:59
struct PyCallBack_ROL_SingletonVector_double_t : public ROL::SingletonVector<double> {
	using ROL::SingletonVector<double>::SingletonVector;

	void set(const class ROL::Vector<double> & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::SingletonVector<double> *>(this), "set");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SingletonVector::set(a0);
	}
	void plus(const class ROL::Vector<double> & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::SingletonVector<double> *>(this), "plus");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SingletonVector::plus(a0);
	}
	void axpy(const double a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::SingletonVector<double> *>(this), "axpy");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SingletonVector::axpy(a0, a1);
	}
	void scale(const double a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::SingletonVector<double> *>(this), "scale");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SingletonVector::scale(a0);
	}
	double dot(const class ROL::Vector<double> & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::SingletonVector<double> *>(this), "dot");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SingletonVector::dot(a0);
	}
	double norm() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::SingletonVector<double> *>(this), "norm");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SingletonVector::norm();
	}
	class Teuchos::RCP<class ROL::Vector<double> > clone() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::SingletonVector<double> *>(this), "clone");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<class Teuchos::RCP<class ROL::Vector<double> >>::value) {
				static pybind11::detail::override_caster_t<class Teuchos::RCP<class ROL::Vector<double> >> caster;
				return pybind11::detail::cast_ref<class Teuchos::RCP<class ROL::Vector<double> >>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Teuchos::RCP<class ROL::Vector<double> >>(std::move(o));
		}
		return SingletonVector::clone();
	}
	class Teuchos::RCP<class ROL::Vector<double> > basis(const int a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::SingletonVector<double> *>(this), "basis");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<class Teuchos::RCP<class ROL::Vector<double> >>::value) {
				static pybind11::detail::override_caster_t<class Teuchos::RCP<class ROL::Vector<double> >> caster;
				return pybind11::detail::cast_ref<class Teuchos::RCP<class ROL::Vector<double> >>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Teuchos::RCP<class ROL::Vector<double> >>(std::move(o));
		}
		return SingletonVector::basis(a0);
	}
	int dimension() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::SingletonVector<double> *>(this), "dimension");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return SingletonVector::dimension();
	}
	void applyUnary(const class ROL::Elementwise::UnaryFunction<double> & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::SingletonVector<double> *>(this), "applyUnary");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SingletonVector::applyUnary(a0);
	}
	void applyBinary(const class ROL::Elementwise::BinaryFunction<double> & a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::SingletonVector<double> *>(this), "applyBinary");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SingletonVector::applyBinary(a0, a1);
	}
	double reduce(const class ROL::Elementwise::ReductionOp<double> & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::SingletonVector<double> *>(this), "reduce");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SingletonVector::reduce(a0);
	}
	void setScalar(const double a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::SingletonVector<double> *>(this), "setScalar");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SingletonVector::setScalar(a0);
	}
	void randomize(const double a0, const double a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::SingletonVector<double> *>(this), "randomize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SingletonVector::randomize(a0, a1);
	}
	void print(std::ostream & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::SingletonVector<double> *>(this), "print");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SingletonVector::print(a0);
	}
	void zero() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::SingletonVector<double> *>(this), "zero");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Vector::zero();
	}
	const class ROL::Vector<double> & dual() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::SingletonVector<double> *>(this), "dual");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const class ROL::Vector<double> &>::value) {
				static pybind11::detail::override_caster_t<const class ROL::Vector<double> &> caster;
				return pybind11::detail::cast_ref<const class ROL::Vector<double> &>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const class ROL::Vector<double> &>(std::move(o));
		}
		return Vector::dual();
	}
	double apply(const class ROL::Vector<double> & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::SingletonVector<double> *>(this), "apply");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return Vector::apply(a0);
	}
};

// ROL::MoreauYosidaObjective file:ROL_MoreauYosidaObjective.hpp line:67
struct PyCallBack_ROL_MoreauYosidaObjective_double_t : public ROL::MoreauYosidaObjective<double> {
	using ROL::MoreauYosidaObjective<double>::MoreauYosidaObjective;

	void update(const class ROL::Vector<double> & a0, enum ROL::UpdateType a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::MoreauYosidaObjective<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return MoreauYosidaObjective::update(a0, a1, a2);
	}
	double value(const class ROL::Vector<double> & a0, double & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::MoreauYosidaObjective<double> *>(this), "value");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return MoreauYosidaObjective::value(a0, a1);
	}
	void gradient(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double & a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::MoreauYosidaObjective<double> *>(this), "gradient");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return MoreauYosidaObjective::gradient(a0, a1, a2);
	}
	void hessVec(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::MoreauYosidaObjective<double> *>(this), "hessVec");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return MoreauYosidaObjective::hessVec(a0, a1, a2, a3);
	}
	void update(const class ROL::Vector<double> & a0, bool a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::MoreauYosidaObjective<double> *>(this), "update");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::MoreauYosidaObjective<double> *>(this), "dirDeriv");
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
	void invHessVec(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::MoreauYosidaObjective<double> *>(this), "invHessVec");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Objective::invHessVec(a0, a1, a2, a3);
	}
	void precond(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::MoreauYosidaObjective<double> *>(this), "precond");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Objective::precond(a0, a1, a2, a3);
	}
};

void bind_ROL_SingletonVector(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // ROL::SingletonVector file:ROL_SingletonVector.hpp line:59
		pybind11::class_<ROL::SingletonVector<double>, Teuchos::RCP<ROL::SingletonVector<double>>, PyCallBack_ROL_SingletonVector_double_t, ROL::Vector<double>> cl(M("ROL"), "SingletonVector_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new ROL::SingletonVector<double>(); }, [](){ return new PyCallBack_ROL_SingletonVector_double_t(); } ), "doc");
		cl.def( pybind11::init<double>(), pybind11::arg("value") );

		cl.def( pybind11::init( [](PyCallBack_ROL_SingletonVector_double_t const &o){ return new PyCallBack_ROL_SingletonVector_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::SingletonVector<double> const &o){ return new ROL::SingletonVector<double>(o); } ) );
		cl.def("getValue", (double (ROL::SingletonVector<double>::*)() const) &ROL::SingletonVector<double>::getValue, "C++: ROL::SingletonVector<double>::getValue() const --> double");
		cl.def("setValue", (void (ROL::SingletonVector<double>::*)(double)) &ROL::SingletonVector<double>::setValue, "C++: ROL::SingletonVector<double>::setValue(double) --> void", pybind11::arg("v"));
		cl.def("set", (void (ROL::SingletonVector<double>::*)(const class ROL::Vector<double> &)) &ROL::SingletonVector<double>::set, "C++: ROL::SingletonVector<double>::set(const class ROL::Vector<double> &) --> void", pybind11::arg("x"));
		cl.def("plus", (void (ROL::SingletonVector<double>::*)(const class ROL::Vector<double> &)) &ROL::SingletonVector<double>::plus, "C++: ROL::SingletonVector<double>::plus(const class ROL::Vector<double> &) --> void", pybind11::arg("x"));
		cl.def("axpy", (void (ROL::SingletonVector<double>::*)(const double, const class ROL::Vector<double> &)) &ROL::SingletonVector<double>::axpy, "C++: ROL::SingletonVector<double>::axpy(const double, const class ROL::Vector<double> &) --> void", pybind11::arg("alpha"), pybind11::arg("x"));
		cl.def("scale", (void (ROL::SingletonVector<double>::*)(const double)) &ROL::SingletonVector<double>::scale, "C++: ROL::SingletonVector<double>::scale(const double) --> void", pybind11::arg("alpha"));
		cl.def("dot", (double (ROL::SingletonVector<double>::*)(const class ROL::Vector<double> &) const) &ROL::SingletonVector<double>::dot, "C++: ROL::SingletonVector<double>::dot(const class ROL::Vector<double> &) const --> double", pybind11::arg("x"));
		cl.def("norm", (double (ROL::SingletonVector<double>::*)() const) &ROL::SingletonVector<double>::norm, "C++: ROL::SingletonVector<double>::norm() const --> double");
		cl.def("clone", (class Teuchos::RCP<class ROL::Vector<double> > (ROL::SingletonVector<double>::*)() const) &ROL::SingletonVector<double>::clone, "C++: ROL::SingletonVector<double>::clone() const --> class Teuchos::RCP<class ROL::Vector<double> >");
		cl.def("basis", (class Teuchos::RCP<class ROL::Vector<double> > (ROL::SingletonVector<double>::*)(const int) const) &ROL::SingletonVector<double>::basis, "C++: ROL::SingletonVector<double>::basis(const int) const --> class Teuchos::RCP<class ROL::Vector<double> >", pybind11::arg("i"));
		cl.def("dimension", (int (ROL::SingletonVector<double>::*)() const) &ROL::SingletonVector<double>::dimension, "C++: ROL::SingletonVector<double>::dimension() const --> int");
		cl.def("applyUnary", (void (ROL::SingletonVector<double>::*)(const class ROL::Elementwise::UnaryFunction<double> &)) &ROL::SingletonVector<double>::applyUnary, "C++: ROL::SingletonVector<double>::applyUnary(const class ROL::Elementwise::UnaryFunction<double> &) --> void", pybind11::arg("f"));
		cl.def("applyBinary", (void (ROL::SingletonVector<double>::*)(const class ROL::Elementwise::BinaryFunction<double> &, const class ROL::Vector<double> &)) &ROL::SingletonVector<double>::applyBinary, "C++: ROL::SingletonVector<double>::applyBinary(const class ROL::Elementwise::BinaryFunction<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("f"), pybind11::arg("x"));
		cl.def("reduce", (double (ROL::SingletonVector<double>::*)(const class ROL::Elementwise::ReductionOp<double> &) const) &ROL::SingletonVector<double>::reduce, "C++: ROL::SingletonVector<double>::reduce(const class ROL::Elementwise::ReductionOp<double> &) const --> double", pybind11::arg("r"));
		cl.def("setScalar", (void (ROL::SingletonVector<double>::*)(const double)) &ROL::SingletonVector<double>::setScalar, "C++: ROL::SingletonVector<double>::setScalar(const double) --> void", pybind11::arg("C"));
		cl.def("randomize", [](ROL::SingletonVector<double> &o) -> void { return o.randomize(); }, "");
		cl.def("randomize", [](ROL::SingletonVector<double> &o, const double & a0) -> void { return o.randomize(a0); }, "", pybind11::arg("l"));
		cl.def("randomize", (void (ROL::SingletonVector<double>::*)(const double, const double)) &ROL::SingletonVector<double>::randomize, "C++: ROL::SingletonVector<double>::randomize(const double, const double) --> void", pybind11::arg("l"), pybind11::arg("u"));
		cl.def("print", (void (ROL::SingletonVector<double>::*)(std::ostream &) const) &ROL::SingletonVector<double>::print, "C++: ROL::SingletonVector<double>::print(std::ostream &) const --> void", pybind11::arg("os"));
		cl.def("assign", (class ROL::SingletonVector<double> & (ROL::SingletonVector<double>::*)(const class ROL::SingletonVector<double> &)) &ROL::SingletonVector<double>::operator=, "C++: ROL::SingletonVector<double>::operator=(const class ROL::SingletonVector<double> &) --> class ROL::SingletonVector<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("plus", (void (ROL::Vector<double>::*)(const class ROL::Vector<double> &)) &ROL::Vector<double>::plus, "C++: ROL::Vector<double>::plus(const class ROL::Vector<double> &) --> void", pybind11::arg("x"));
		cl.def("scale", (void (ROL::Vector<double>::*)(const double)) &ROL::Vector<double>::scale, "C++: ROL::Vector<double>::scale(const double) --> void", pybind11::arg("alpha"));
		cl.def("dot", (double (ROL::Vector<double>::*)(const class ROL::Vector<double> &) const) &ROL::Vector<double>::dot, "C++: ROL::Vector<double>::dot(const class ROL::Vector<double> &) const --> double", pybind11::arg("x"));
		cl.def("norm", (double (ROL::Vector<double>::*)() const) &ROL::Vector<double>::norm, "C++: ROL::Vector<double>::norm() const --> double");
		cl.def("clone", (class Teuchos::RCP<class ROL::Vector<double> > (ROL::Vector<double>::*)() const) &ROL::Vector<double>::clone, "C++: ROL::Vector<double>::clone() const --> class Teuchos::RCP<class ROL::Vector<double> >");
		cl.def("axpy", (void (ROL::Vector<double>::*)(const double, const class ROL::Vector<double> &)) &ROL::Vector<double>::axpy, "Compute \n where \n.\n\n             \n is the scaling of \n             \n\n     is a vector.\n\n             On return \n.\n             Uses #clone, #set, #scale and #plus for the computation.\n             Please overload if a more efficient implementation is needed.\n\n             ---\n\nC++: ROL::Vector<double>::axpy(const double, const class ROL::Vector<double> &) --> void", pybind11::arg("alpha"), pybind11::arg("x"));
		cl.def("zero", (void (ROL::Vector<double>::*)()) &ROL::Vector<double>::zero, "Set to zero vector.\n\n             Uses #scale by zero for the computation.\n             Please overload if a more efficient implementation is needed.\n\n             ---\n\nC++: ROL::Vector<double>::zero() --> void");
		cl.def("basis", (class Teuchos::RCP<class ROL::Vector<double> > (ROL::Vector<double>::*)(const int) const) &ROL::Vector<double>::basis, "Return i-th basis vector.\n\n             \n is the index of the basis function.\n             \n\n A reference-counted pointer to the basis vector with index \n\n             Overloading the basis is only required if the default gradient implementation\n             is used, which computes a finite-difference approximation.\n\n             ---\n\nC++: ROL::Vector<double>::basis(const int) const --> class Teuchos::RCP<class ROL::Vector<double> >", pybind11::arg("i"));
		cl.def("dimension", (int (ROL::Vector<double>::*)() const) &ROL::Vector<double>::dimension, "Return dimension of the vector space.\n\n             \n The dimension of the vector space, i.e., the total number of basis vectors.\n\n             Overload if the basis is overloaded.\n\n             ---\n\nC++: ROL::Vector<double>::dimension() const --> int");
		cl.def("set", (void (ROL::Vector<double>::*)(const class ROL::Vector<double> &)) &ROL::Vector<double>::set, "Set \n where \n.\n\n             \n     is a vector.\n\n             On return \n.\n             Uses #zero and #plus methods for the computation.\n             Please overload if a more efficient implementation is needed.\n\n             ---\n\nC++: ROL::Vector<double>::set(const class ROL::Vector<double> &) --> void", pybind11::arg("x"));
		cl.def("dual", (const class ROL::Vector<double> & (ROL::Vector<double>::*)() const) &ROL::Vector<double>::dual, "Return dual representation of \n, for example,\n             the result of applying a Riesz map, or change of basis, or\n             change of memory layout.\n\n             \n         A const reference to dual representation.\n\n             By default, returns the current object.\n             Please overload if you need a dual representation.\n\n             ---\n\nC++: ROL::Vector<double>::dual() const --> const class ROL::Vector<double> &", pybind11::return_value_policy::automatic);
		cl.def("apply", (double (ROL::Vector<double>::*)(const class ROL::Vector<double> &) const) &ROL::Vector<double>::apply, "Apply \n to a dual vector.  This is equivalent\n             to the call \n\n.\n\n             \n      is a vector\n             \n\n         The number equal to \n.\n\n             ---\n\nC++: ROL::Vector<double>::apply(const class ROL::Vector<double> &) const --> double", pybind11::arg("x"));
		cl.def("applyUnary", (void (ROL::Vector<double>::*)(const class ROL::Elementwise::UnaryFunction<double> &)) &ROL::Vector<double>::applyUnary, "C++: ROL::Vector<double>::applyUnary(const class ROL::Elementwise::UnaryFunction<double> &) --> void", pybind11::arg("f"));
		cl.def("applyBinary", (void (ROL::Vector<double>::*)(const class ROL::Elementwise::BinaryFunction<double> &, const class ROL::Vector<double> &)) &ROL::Vector<double>::applyBinary, "C++: ROL::Vector<double>::applyBinary(const class ROL::Elementwise::BinaryFunction<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("f"), pybind11::arg("x"));
		cl.def("reduce", (double (ROL::Vector<double>::*)(const class ROL::Elementwise::ReductionOp<double> &) const) &ROL::Vector<double>::reduce, "C++: ROL::Vector<double>::reduce(const class ROL::Elementwise::ReductionOp<double> &) const --> double", pybind11::arg("r"));
		cl.def("print", (void (ROL::Vector<double>::*)(std::ostream &) const) &ROL::Vector<double>::print, "C++: ROL::Vector<double>::print(std::ostream &) const --> void", pybind11::arg("outStream"));
		cl.def("setScalar", (void (ROL::Vector<double>::*)(const double)) &ROL::Vector<double>::setScalar, "Set \n where \n.\n\n             \n     is a scalar.\n\n             On return \n.\n             Uses #applyUnary methods for the computation.\n             Please overload if a more efficient implementation is needed.\n\n             ---\n\nC++: ROL::Vector<double>::setScalar(const double) --> void", pybind11::arg("C"));
		cl.def("randomize", [](ROL::Vector<double> &o) -> void { return o.randomize(); }, "");
		cl.def("randomize", [](ROL::Vector<double> &o, const double & a0) -> void { return o.randomize(a0); }, "", pybind11::arg("l"));
		cl.def("randomize", (void (ROL::Vector<double>::*)(const double, const double)) &ROL::Vector<double>::randomize, "Set vector to be uniform random between [l,u].\n\n             \n     is a the lower bound.\n             \n\n     is a the upper bound.\n\n             On return the components of \n are uniform\n             random numbers on the interval \n\n.\n       	     The default implementation uses #applyUnary methods for the\n       	     computation. Please overload if a more efficient implementation is\n             needed.\n\n             ---\n\nC++: ROL::Vector<double>::randomize(const double, const double) --> void", pybind11::arg("l"), pybind11::arg("u"));
		cl.def("assign", (class ROL::Vector<double> & (ROL::Vector<double>::*)(const class ROL::Vector<double> &)) &ROL::Vector<double>::operator=, "C++: ROL::Vector<double>::operator=(const class ROL::Vector<double> &) --> class ROL::Vector<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::ScalarController file:ROL_ScalarController.hpp line:54
		pybind11::class_<ROL::ScalarController<double,int>, Teuchos::RCP<ROL::ScalarController<double,int>>, ROL::VectorController<double,int>> cl(M("ROL"), "ScalarController_double_int_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new ROL::ScalarController<double,int>(); } ) );
		cl.def( pybind11::init( [](ROL::ScalarController<double,int> const &o){ return new ROL::ScalarController<double,int>(o); } ) );
		cl.def("get", (bool (ROL::ScalarController<double,int>::*)(double &, const int &)) &ROL::ScalarController<double, int>::get, "C++: ROL::ScalarController<double, int>::get(double &, const int &) --> bool", pybind11::arg("x"), pybind11::arg("param"));
		cl.def("set", (void (ROL::ScalarController<double,int>::*)(double, const int &)) &ROL::ScalarController<double, int>::set, "C++: ROL::ScalarController<double, int>::set(double, const int &) --> void", pybind11::arg("x"), pybind11::arg("param"));
		cl.def("reset", [](ROL::VectorController<double,int> &o) -> void { return o.reset(); }, "");
		cl.def("reset", (void (ROL::VectorController<double,int>::*)(bool)) &ROL::VectorController<double, int>::reset, "C++: ROL::VectorController<double, int>::reset(bool) --> void", pybind11::arg("flag"));
		cl.def("objectiveUpdate", [](ROL::VectorController<double,int> &o) -> void { return o.objectiveUpdate(); }, "");
		cl.def("objectiveUpdate", (void (ROL::VectorController<double,int>::*)(bool)) &ROL::VectorController<double, int>::objectiveUpdate, "C++: ROL::VectorController<double, int>::objectiveUpdate(bool) --> void", pybind11::arg("flag"));
		cl.def("constraintUpdate", [](ROL::VectorController<double,int> &o) -> void { return o.constraintUpdate(); }, "");
		cl.def("constraintUpdate", (void (ROL::VectorController<double,int>::*)(bool)) &ROL::VectorController<double, int>::constraintUpdate, "C++: ROL::VectorController<double, int>::constraintUpdate(bool) --> void", pybind11::arg("flag"));
		cl.def("objectiveUpdate", (void (ROL::VectorController<double,int>::*)(enum ROL::UpdateType)) &ROL::VectorController<double, int>::objectiveUpdate, "C++: ROL::VectorController<double, int>::objectiveUpdate(enum ROL::UpdateType) --> void", pybind11::arg("type"));
		cl.def("constraintUpdate", (void (ROL::VectorController<double,int>::*)(enum ROL::UpdateType)) &ROL::VectorController<double, int>::constraintUpdate, "C++: ROL::VectorController<double, int>::constraintUpdate(enum ROL::UpdateType) --> void", pybind11::arg("type"));
		cl.def("isNull", (bool (ROL::VectorController<double,int>::*)(const int &) const) &ROL::VectorController<double, int>::isNull, "C++: ROL::VectorController<double, int>::isNull(const int &) const --> bool", pybind11::arg("param"));
		cl.def("isComputed", (bool (ROL::VectorController<double,int>::*)(const int &) const) &ROL::VectorController<double, int>::isComputed, "C++: ROL::VectorController<double, int>::isComputed(const int &) const --> bool", pybind11::arg("param"));
		cl.def("allocate", (void (ROL::VectorController<double,int>::*)(const class ROL::Vector<double> &, const int &)) &ROL::VectorController<double, int>::allocate, "C++: ROL::VectorController<double, int>::allocate(const class ROL::Vector<double> &, const int &) --> void", pybind11::arg("x"), pybind11::arg("param"));
		cl.def("set", (const class Teuchos::RCP<class ROL::Vector<double> > (ROL::VectorController<double,int>::*)(const int &)) &ROL::VectorController<double, int>::set, "C++: ROL::VectorController<double, int>::set(const int &) --> const class Teuchos::RCP<class ROL::Vector<double> >", pybind11::arg("param"));
		cl.def("get", (const class Teuchos::RCP<const class ROL::Vector<double> > (ROL::VectorController<double,int>::*)(const int &) const) &ROL::VectorController<double, int>::get, "C++: ROL::VectorController<double, int>::get(const int &) const --> const class Teuchos::RCP<const class ROL::Vector<double> >", pybind11::arg("param"));
		cl.def("get", (bool (ROL::VectorController<double,int>::*)(class ROL::Vector<double> &, const int &)) &ROL::VectorController<double, int>::get, "C++: ROL::VectorController<double, int>::get(class ROL::Vector<double> &, const int &) --> bool", pybind11::arg("x"), pybind11::arg("param"));
		cl.def("set", (void (ROL::VectorController<double,int>::*)(const class ROL::Vector<double> &, const int &)) &ROL::VectorController<double, int>::set, "C++: ROL::VectorController<double, int>::set(const class ROL::Vector<double> &, const int &) --> void", pybind11::arg("x"), pybind11::arg("param"));
		cl.def("push", (void (ROL::VectorController<double,int>::*)(class ROL::VectorController<double, int> &) const) &ROL::VectorController<double, int>::push, "C++: ROL::VectorController<double, int>::push(class ROL::VectorController<double, int> &) const --> void", pybind11::arg("to"));
	}
	{ // ROL::MoreauYosidaObjective file:ROL_MoreauYosidaObjective.hpp line:67
		pybind11::class_<ROL::MoreauYosidaObjective<double>, Teuchos::RCP<ROL::MoreauYosidaObjective<double>>, PyCallBack_ROL_MoreauYosidaObjective_double_t, ROL::Objective<double>> cl(M("ROL"), "MoreauYosidaObjective_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](const class Teuchos::RCP<class ROL::Objective<double> > & a0, const class Teuchos::RCP<class ROL::BoundConstraint<double> > & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3){ return new ROL::MoreauYosidaObjective<double>(a0, a1, a2, a3); }, [](const class Teuchos::RCP<class ROL::Objective<double> > & a0, const class Teuchos::RCP<class ROL::BoundConstraint<double> > & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3){ return new PyCallBack_ROL_MoreauYosidaObjective_double_t(a0, a1, a2, a3); } ), "doc");
		cl.def( pybind11::init( [](const class Teuchos::RCP<class ROL::Objective<double> > & a0, const class Teuchos::RCP<class ROL::BoundConstraint<double> > & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const double & a4){ return new ROL::MoreauYosidaObjective<double>(a0, a1, a2, a3, a4); }, [](const class Teuchos::RCP<class ROL::Objective<double> > & a0, const class Teuchos::RCP<class ROL::BoundConstraint<double> > & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const double & a4){ return new PyCallBack_ROL_MoreauYosidaObjective_double_t(a0, a1, a2, a3, a4); } ), "doc");
		cl.def( pybind11::init( [](const class Teuchos::RCP<class ROL::Objective<double> > & a0, const class Teuchos::RCP<class ROL::BoundConstraint<double> > & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const double & a4, const bool & a5){ return new ROL::MoreauYosidaObjective<double>(a0, a1, a2, a3, a4, a5); }, [](const class Teuchos::RCP<class ROL::Objective<double> > & a0, const class Teuchos::RCP<class ROL::BoundConstraint<double> > & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const double & a4, const bool & a5){ return new PyCallBack_ROL_MoreauYosidaObjective_double_t(a0, a1, a2, a3, a4, a5); } ), "doc");
		cl.def( pybind11::init<const class Teuchos::RCP<class ROL::Objective<double> > &, const class Teuchos::RCP<class ROL::BoundConstraint<double> > &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const double, const bool, const bool>(), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("mu"), pybind11::arg("updateMultiplier"), pybind11::arg("updatePenalty") );

		cl.def( pybind11::init( [](const class Teuchos::RCP<class ROL::Objective<double> > & a0, const class Teuchos::RCP<class ROL::BoundConstraint<double> > & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class ROL::Vector<double> & a4){ return new ROL::MoreauYosidaObjective<double>(a0, a1, a2, a3, a4); }, [](const class Teuchos::RCP<class ROL::Objective<double> > & a0, const class Teuchos::RCP<class ROL::BoundConstraint<double> > & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class ROL::Vector<double> & a4){ return new PyCallBack_ROL_MoreauYosidaObjective_double_t(a0, a1, a2, a3, a4); } ), "doc");
		cl.def( pybind11::init( [](const class Teuchos::RCP<class ROL::Objective<double> > & a0, const class Teuchos::RCP<class ROL::BoundConstraint<double> > & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class ROL::Vector<double> & a4, const double & a5){ return new ROL::MoreauYosidaObjective<double>(a0, a1, a2, a3, a4, a5); }, [](const class Teuchos::RCP<class ROL::Objective<double> > & a0, const class Teuchos::RCP<class ROL::BoundConstraint<double> > & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class ROL::Vector<double> & a4, const double & a5){ return new PyCallBack_ROL_MoreauYosidaObjective_double_t(a0, a1, a2, a3, a4, a5); } ), "doc");
		cl.def( pybind11::init( [](const class Teuchos::RCP<class ROL::Objective<double> > & a0, const class Teuchos::RCP<class ROL::BoundConstraint<double> > & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class ROL::Vector<double> & a4, const double & a5, const bool & a6){ return new ROL::MoreauYosidaObjective<double>(a0, a1, a2, a3, a4, a5, a6); }, [](const class Teuchos::RCP<class ROL::Objective<double> > & a0, const class Teuchos::RCP<class ROL::BoundConstraint<double> > & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class ROL::Vector<double> & a4, const double & a5, const bool & a6){ return new PyCallBack_ROL_MoreauYosidaObjective_double_t(a0, a1, a2, a3, a4, a5, a6); } ), "doc");
		cl.def( pybind11::init<const class Teuchos::RCP<class ROL::Objective<double> > &, const class Teuchos::RCP<class ROL::BoundConstraint<double> > &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const double, const bool, const bool>(), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("lam"), pybind11::arg("mu"), pybind11::arg("updateMultiplier"), pybind11::arg("updatePenalty") );

		cl.def( pybind11::init<const class Teuchos::RCP<class ROL::Objective<double> > &, const class Teuchos::RCP<class ROL::BoundConstraint<double> > &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class Teuchos::ParameterList &>(), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("parlist") );

		cl.def( pybind11::init<const class Teuchos::RCP<class ROL::Objective<double> > &, const class Teuchos::RCP<class ROL::BoundConstraint<double> > &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class Teuchos::ParameterList &>(), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("lam"), pybind11::arg("parlist") );

		cl.def( pybind11::init( [](PyCallBack_ROL_MoreauYosidaObjective_double_t const &o){ return new PyCallBack_ROL_MoreauYosidaObjective_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::MoreauYosidaObjective<double> const &o){ return new ROL::MoreauYosidaObjective<double>(o); } ) );
		cl.def("updateMultipliers", (void (ROL::MoreauYosidaObjective<double>::*)(double, const class ROL::Vector<double> &)) &ROL::MoreauYosidaObjective<double>::updateMultipliers, "C++: ROL::MoreauYosidaObjective<double>::updateMultipliers(double, const class ROL::Vector<double> &) --> void", pybind11::arg("mu"), pybind11::arg("x"));
		cl.def("reset", (void (ROL::MoreauYosidaObjective<double>::*)(const double)) &ROL::MoreauYosidaObjective<double>::reset, "C++: ROL::MoreauYosidaObjective<double>::reset(const double) --> void", pybind11::arg("mu"));
		cl.def("testComplementarity", (double (ROL::MoreauYosidaObjective<double>::*)(const class ROL::Vector<double> &)) &ROL::MoreauYosidaObjective<double>::testComplementarity, "C++: ROL::MoreauYosidaObjective<double>::testComplementarity(const class ROL::Vector<double> &) --> double", pybind11::arg("x"));
		cl.def("getObjectiveValue", (double (ROL::MoreauYosidaObjective<double>::*)(const class ROL::Vector<double> &, double &)) &ROL::MoreauYosidaObjective<double>::getObjectiveValue, "C++: ROL::MoreauYosidaObjective<double>::getObjectiveValue(const class ROL::Vector<double> &, double &) --> double", pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("getObjectiveGradient", (void (ROL::MoreauYosidaObjective<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::MoreauYosidaObjective<double>::getObjectiveGradient, "C++: ROL::MoreauYosidaObjective<double>::getObjectiveGradient(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("getNumberFunctionEvaluations", (int (ROL::MoreauYosidaObjective<double>::*)()) &ROL::MoreauYosidaObjective<double>::getNumberFunctionEvaluations, "C++: ROL::MoreauYosidaObjective<double>::getNumberFunctionEvaluations() --> int");
		cl.def("getNumberGradientEvaluations", (int (ROL::MoreauYosidaObjective<double>::*)()) &ROL::MoreauYosidaObjective<double>::getNumberGradientEvaluations, "C++: ROL::MoreauYosidaObjective<double>::getNumberGradientEvaluations() --> int");
		cl.def("update", [](ROL::MoreauYosidaObjective<double> &o, const class ROL::Vector<double> & a0, enum ROL::UpdateType const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("type"));
		cl.def("update", (void (ROL::MoreauYosidaObjective<double>::*)(const class ROL::Vector<double> &, enum ROL::UpdateType, int)) &ROL::MoreauYosidaObjective<double>::update, "C++: ROL::MoreauYosidaObjective<double>::update(const class ROL::Vector<double> &, enum ROL::UpdateType, int) --> void", pybind11::arg("x"), pybind11::arg("type"), pybind11::arg("iter"));
		cl.def("value", (double (ROL::MoreauYosidaObjective<double>::*)(const class ROL::Vector<double> &, double &)) &ROL::MoreauYosidaObjective<double>::value, "C++: ROL::MoreauYosidaObjective<double>::value(const class ROL::Vector<double> &, double &) --> double", pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("gradient", (void (ROL::MoreauYosidaObjective<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::MoreauYosidaObjective<double>::gradient, "C++: ROL::MoreauYosidaObjective<double>::gradient(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("hessVec", (void (ROL::MoreauYosidaObjective<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::MoreauYosidaObjective<double>::hessVec, "C++: ROL::MoreauYosidaObjective<double>::hessVec(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("hv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("tol"));
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
}
