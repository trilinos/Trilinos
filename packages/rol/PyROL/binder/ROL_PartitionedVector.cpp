#include <ROL_Constraint.hpp>
#include <ROL_Constraint_Partitioned.hpp>
#include <ROL_Elementwise_Function.hpp>
#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_PartitionedVector.hpp>
#include <ROL_UpdateType.hpp>
#include <ROL_Vector.hpp>
#include <Teuchos_ENull.hpp>
#include <Teuchos_RCPDecl.hpp>
#include <Teuchos_RCPNode.hpp>
#include <Teuchos_any.hpp>
#include <cwchar>
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

// ROL::PartitionedVector file:ROL_PartitionedVector.hpp line:60
struct PyCallBack_ROL_PartitionedVector_double_t : public ROL::PartitionedVector<double> {
	using ROL::PartitionedVector<double>::PartitionedVector;

	void set(const class ROL::Vector<double> & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::PartitionedVector<double> *>(this), "set");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PartitionedVector::set(a0);
	}
	void plus(const class ROL::Vector<double> & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::PartitionedVector<double> *>(this), "plus");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PartitionedVector::plus(a0);
	}
	void scale(const double a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::PartitionedVector<double> *>(this), "scale");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PartitionedVector::scale(a0);
	}
	void axpy(const double a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::PartitionedVector<double> *>(this), "axpy");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PartitionedVector::axpy(a0, a1);
	}
	double dot(const class ROL::Vector<double> & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::PartitionedVector<double> *>(this), "dot");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PartitionedVector::dot(a0);
	}
	double norm() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::PartitionedVector<double> *>(this), "norm");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PartitionedVector::norm();
	}
	class Teuchos::RCP<class ROL::Vector<double> > clone() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::PartitionedVector<double> *>(this), "clone");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<class Teuchos::RCP<class ROL::Vector<double> >>::value) {
				static pybind11::detail::override_caster_t<class Teuchos::RCP<class ROL::Vector<double> >> caster;
				return pybind11::detail::cast_ref<class Teuchos::RCP<class ROL::Vector<double> >>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Teuchos::RCP<class ROL::Vector<double> >>(std::move(o));
		}
		return PartitionedVector::clone();
	}
	const class ROL::Vector<double> & dual() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::PartitionedVector<double> *>(this), "dual");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const class ROL::Vector<double> &>::value) {
				static pybind11::detail::override_caster_t<const class ROL::Vector<double> &> caster;
				return pybind11::detail::cast_ref<const class ROL::Vector<double> &>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const class ROL::Vector<double> &>(std::move(o));
		}
		return PartitionedVector::dual();
	}
	double apply(const class ROL::Vector<double> & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::PartitionedVector<double> *>(this), "apply");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PartitionedVector::apply(a0);
	}
	class Teuchos::RCP<class ROL::Vector<double> > basis(const int a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::PartitionedVector<double> *>(this), "basis");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<class Teuchos::RCP<class ROL::Vector<double> >>::value) {
				static pybind11::detail::override_caster_t<class Teuchos::RCP<class ROL::Vector<double> >> caster;
				return pybind11::detail::cast_ref<class Teuchos::RCP<class ROL::Vector<double> >>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Teuchos::RCP<class ROL::Vector<double> >>(std::move(o));
		}
		return PartitionedVector::basis(a0);
	}
	int dimension() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::PartitionedVector<double> *>(this), "dimension");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return PartitionedVector::dimension();
	}
	void zero() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::PartitionedVector<double> *>(this), "zero");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PartitionedVector::zero();
	}
	void applyUnary(const class ROL::Elementwise::UnaryFunction<double> & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::PartitionedVector<double> *>(this), "applyUnary");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PartitionedVector::applyUnary(a0);
	}
	void applyBinary(const class ROL::Elementwise::BinaryFunction<double> & a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::PartitionedVector<double> *>(this), "applyBinary");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PartitionedVector::applyBinary(a0, a1);
	}
	double reduce(const class ROL::Elementwise::ReductionOp<double> & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::PartitionedVector<double> *>(this), "reduce");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return PartitionedVector::reduce(a0);
	}
	void setScalar(const double a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::PartitionedVector<double> *>(this), "setScalar");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PartitionedVector::setScalar(a0);
	}
	void randomize(const double a0, const double a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::PartitionedVector<double> *>(this), "randomize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PartitionedVector::randomize(a0, a1);
	}
};

// ROL::Constraint_Partitioned file:ROL_Constraint_Partitioned.hpp line:56
struct PyCallBack_ROL_Constraint_Partitioned_double_t : public ROL::Constraint_Partitioned<double> {
	using ROL::Constraint_Partitioned<double>::Constraint_Partitioned;

	void update(const class ROL::Vector<double> & a0, enum ROL::UpdateType a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_Partitioned<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_Partitioned::update(a0, a1, a2);
	}
	void update(const class ROL::Vector<double> & a0, bool a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_Partitioned<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_Partitioned::update(a0, a1, a2);
	}
	void value(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double & a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_Partitioned<double> *>(this), "value");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_Partitioned::value(a0, a1, a2);
	}
	void applyJacobian(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_Partitioned<double> *>(this), "applyJacobian");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_Partitioned::applyJacobian(a0, a1, a2, a3);
	}
	void applyAdjointJacobian(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_Partitioned<double> *>(this), "applyAdjointJacobian");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_Partitioned::applyAdjointJacobian(a0, a1, a2, a3);
	}
	void applyAdjointHessian(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, double & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_Partitioned<double> *>(this), "applyAdjointHessian");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_Partitioned::applyAdjointHessian(a0, a1, a2, a3, a4);
	}
	void applyPreconditioner(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, double & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_Partitioned<double> *>(this), "applyPreconditioner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_Partitioned::applyPreconditioner(a0, a1, a2, a3, a4);
	}
	void applyAdjointJacobian(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, double & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_Partitioned<double> *>(this), "applyAdjointJacobian");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint::applyAdjointJacobian(a0, a1, a2, a3, a4);
	}
};

void bind_ROL_PartitionedVector(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // ROL::PartitionedVector file:ROL_PartitionedVector.hpp line:60
		pybind11::class_<ROL::PartitionedVector<double>, Teuchos::RCP<ROL::PartitionedVector<double>>, PyCallBack_ROL_PartitionedVector_double_t, ROL::Vector<double>> cl(M("ROL"), "PartitionedVector_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](PyCallBack_ROL_PartitionedVector_double_t const &o){ return new PyCallBack_ROL_PartitionedVector_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::PartitionedVector<double> const &o){ return new ROL::PartitionedVector<double>(o); } ) );
		cl.def("set", (void (ROL::PartitionedVector<double>::*)(const class ROL::Vector<double> &)) &ROL::PartitionedVector<double>::set, "C++: ROL::PartitionedVector<double>::set(const class ROL::Vector<double> &) --> void", pybind11::arg("x"));
		cl.def("plus", (void (ROL::PartitionedVector<double>::*)(const class ROL::Vector<double> &)) &ROL::PartitionedVector<double>::plus, "C++: ROL::PartitionedVector<double>::plus(const class ROL::Vector<double> &) --> void", pybind11::arg("x"));
		cl.def("scale", (void (ROL::PartitionedVector<double>::*)(const double)) &ROL::PartitionedVector<double>::scale, "C++: ROL::PartitionedVector<double>::scale(const double) --> void", pybind11::arg("alpha"));
		cl.def("axpy", (void (ROL::PartitionedVector<double>::*)(const double, const class ROL::Vector<double> &)) &ROL::PartitionedVector<double>::axpy, "C++: ROL::PartitionedVector<double>::axpy(const double, const class ROL::Vector<double> &) --> void", pybind11::arg("alpha"), pybind11::arg("x"));
		cl.def("dot", (double (ROL::PartitionedVector<double>::*)(const class ROL::Vector<double> &) const) &ROL::PartitionedVector<double>::dot, "C++: ROL::PartitionedVector<double>::dot(const class ROL::Vector<double> &) const --> double", pybind11::arg("x"));
		cl.def("norm", (double (ROL::PartitionedVector<double>::*)() const) &ROL::PartitionedVector<double>::norm, "C++: ROL::PartitionedVector<double>::norm() const --> double");
		cl.def("clone", (class Teuchos::RCP<class ROL::Vector<double> > (ROL::PartitionedVector<double>::*)() const) &ROL::PartitionedVector<double>::clone, "C++: ROL::PartitionedVector<double>::clone() const --> class Teuchos::RCP<class ROL::Vector<double> >");
		cl.def("dual", (const class ROL::Vector<double> & (ROL::PartitionedVector<double>::*)() const) &ROL::PartitionedVector<double>::dual, "C++: ROL::PartitionedVector<double>::dual() const --> const class ROL::Vector<double> &", pybind11::return_value_policy::automatic);
		cl.def("apply", (double (ROL::PartitionedVector<double>::*)(const class ROL::Vector<double> &) const) &ROL::PartitionedVector<double>::apply, "C++: ROL::PartitionedVector<double>::apply(const class ROL::Vector<double> &) const --> double", pybind11::arg("x"));
		cl.def("basis", (class Teuchos::RCP<class ROL::Vector<double> > (ROL::PartitionedVector<double>::*)(const int) const) &ROL::PartitionedVector<double>::basis, "C++: ROL::PartitionedVector<double>::basis(const int) const --> class Teuchos::RCP<class ROL::Vector<double> >", pybind11::arg("i"));
		cl.def("dimension", (int (ROL::PartitionedVector<double>::*)() const) &ROL::PartitionedVector<double>::dimension, "C++: ROL::PartitionedVector<double>::dimension() const --> int");
		cl.def("zero", (void (ROL::PartitionedVector<double>::*)()) &ROL::PartitionedVector<double>::zero, "C++: ROL::PartitionedVector<double>::zero() --> void");
		cl.def("applyUnary", (void (ROL::PartitionedVector<double>::*)(const class ROL::Elementwise::UnaryFunction<double> &)) &ROL::PartitionedVector<double>::applyUnary, "C++: ROL::PartitionedVector<double>::applyUnary(const class ROL::Elementwise::UnaryFunction<double> &) --> void", pybind11::arg("f"));
		cl.def("applyBinary", (void (ROL::PartitionedVector<double>::*)(const class ROL::Elementwise::BinaryFunction<double> &, const class ROL::Vector<double> &)) &ROL::PartitionedVector<double>::applyBinary, "C++: ROL::PartitionedVector<double>::applyBinary(const class ROL::Elementwise::BinaryFunction<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("f"), pybind11::arg("x"));
		cl.def("reduce", (double (ROL::PartitionedVector<double>::*)(const class ROL::Elementwise::ReductionOp<double> &) const) &ROL::PartitionedVector<double>::reduce, "C++: ROL::PartitionedVector<double>::reduce(const class ROL::Elementwise::ReductionOp<double> &) const --> double", pybind11::arg("r"));
		cl.def("setScalar", (void (ROL::PartitionedVector<double>::*)(const double)) &ROL::PartitionedVector<double>::setScalar, "C++: ROL::PartitionedVector<double>::setScalar(const double) --> void", pybind11::arg("C"));
		cl.def("randomize", [](ROL::PartitionedVector<double> &o) -> void { return o.randomize(); }, "");
		cl.def("randomize", [](ROL::PartitionedVector<double> &o, const double & a0) -> void { return o.randomize(a0); }, "", pybind11::arg("l"));
		cl.def("randomize", (void (ROL::PartitionedVector<double>::*)(const double, const double)) &ROL::PartitionedVector<double>::randomize, "C++: ROL::PartitionedVector<double>::randomize(const double, const double) --> void", pybind11::arg("l"), pybind11::arg("u"));
		cl.def("__getitem__", (class ROL::Vector<double> & (ROL::PartitionedVector<double>::*)(unsigned long)) &ROL::PartitionedVector<double>::operator[], "C++: ROL::PartitionedVector<double>::operator[](unsigned long) --> class ROL::Vector<double> &", pybind11::return_value_policy::automatic, pybind11::arg("i"));
		cl.def("get", (class Teuchos::RCP<class ROL::Vector<double> > (ROL::PartitionedVector<double>::*)(unsigned long)) &ROL::PartitionedVector<double>::get, "C++: ROL::PartitionedVector<double>::get(unsigned long) --> class Teuchos::RCP<class ROL::Vector<double> >", pybind11::arg("i"));
		cl.def("set", (void (ROL::PartitionedVector<double>::*)(unsigned long, const class ROL::Vector<double> &)) &ROL::PartitionedVector<double>::set, "C++: ROL::PartitionedVector<double>::set(unsigned long, const class ROL::Vector<double> &) --> void", pybind11::arg("i"), pybind11::arg("x"));
		cl.def("zero", (void (ROL::PartitionedVector<double>::*)(unsigned long)) &ROL::PartitionedVector<double>::zero, "C++: ROL::PartitionedVector<double>::zero(unsigned long) --> void", pybind11::arg("i"));
		cl.def("numVectors", (unsigned long (ROL::PartitionedVector<double>::*)() const) &ROL::PartitionedVector<double>::numVectors, "C++: ROL::PartitionedVector<double>::numVectors() const --> unsigned long");
		cl.def_static("create", (class Teuchos::RCP<class ROL::PartitionedVector<double> > (*)(const class ROL::Vector<double> &, unsigned long)) &ROL::PartitionedVector<double>::create, "C++: ROL::PartitionedVector<double>::create(const class ROL::Vector<double> &, unsigned long) --> class Teuchos::RCP<class ROL::PartitionedVector<double> >", pybind11::arg("x"), pybind11::arg("N"));
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
		cl.def("setScalar", (void (ROL::Vector<double>::*)(const double)) &ROL::Vector<double>::setScalar, "Set \n where \n.\n\n             \n     is a scalar.\n\n             On return \n.\n             Uses #applyUnary methods for the computation.\n             Please overload if a more efficient implementation is needed.\n\n             ---\n\nC++: ROL::Vector<double>::setScalar(const double) --> void", pybind11::arg("C"));
		cl.def("randomize", [](ROL::Vector<double> &o) -> void { return o.randomize(); }, "");
		cl.def("randomize", [](ROL::Vector<double> &o, const double & a0) -> void { return o.randomize(a0); }, "", pybind11::arg("l"));
		cl.def("randomize", (void (ROL::Vector<double>::*)(const double, const double)) &ROL::Vector<double>::randomize, "Set vector to be uniform random between [l,u].\n\n             \n     is a the lower bound.\n             \n\n     is a the upper bound.\n\n             On return the components of \n are uniform\n             random numbers on the interval \n\n.\n       	     The default implementation uses #applyUnary methods for the\n       	     computation. Please overload if a more efficient implementation is\n             needed.\n\n             ---\n\nC++: ROL::Vector<double>::randomize(const double, const double) --> void", pybind11::arg("l"), pybind11::arg("u"));
		cl.def("assign", (class ROL::Vector<double> & (ROL::Vector<double>::*)(const class ROL::Vector<double> &)) &ROL::Vector<double>::operator=, "C++: ROL::Vector<double>::operator=(const class ROL::Vector<double> &) --> class ROL::Vector<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::Constraint_Partitioned file:ROL_Constraint_Partitioned.hpp line:56
		pybind11::class_<ROL::Constraint_Partitioned<double>, Teuchos::RCP<ROL::Constraint_Partitioned<double>>, PyCallBack_ROL_Constraint_Partitioned_double_t, ROL::Constraint<double>> cl(M("ROL"), "Constraint_Partitioned_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](PyCallBack_ROL_Constraint_Partitioned_double_t const &o){ return new PyCallBack_ROL_Constraint_Partitioned_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::Constraint_Partitioned<double> const &o){ return new ROL::Constraint_Partitioned<double>(o); } ) );
		cl.def("applyAdjointJacobian", [](ROL::Constraint_Partitioned<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, double & a4) -> void { return o.applyAdjointJacobian(a0, a1, a2, a3, a4); }, "", pybind11::arg("ajv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("dualv"), pybind11::arg("tol"));
		cl.def("getNumberConstraintEvaluations", (int (ROL::Constraint_Partitioned<double>::*)() const) &ROL::Constraint_Partitioned<double>::getNumberConstraintEvaluations, "C++: ROL::Constraint_Partitioned<double>::getNumberConstraintEvaluations() const --> int");
		cl.def("get", [](ROL::Constraint_Partitioned<double> const &o) -> Teuchos::RCP<class ROL::Constraint<double> > { return o.get(); }, "");
		cl.def("get", (class Teuchos::RCP<class ROL::Constraint<double> > (ROL::Constraint_Partitioned<double>::*)(int) const) &ROL::Constraint_Partitioned<double>::get, "C++: ROL::Constraint_Partitioned<double>::get(int) const --> class Teuchos::RCP<class ROL::Constraint<double> >", pybind11::arg("ind"));
		cl.def("update", [](ROL::Constraint_Partitioned<double> &o, const class ROL::Vector<double> & a0, enum ROL::UpdateType const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("type"));
		cl.def("update", (void (ROL::Constraint_Partitioned<double>::*)(const class ROL::Vector<double> &, enum ROL::UpdateType, int)) &ROL::Constraint_Partitioned<double>::update, "C++: ROL::Constraint_Partitioned<double>::update(const class ROL::Vector<double> &, enum ROL::UpdateType, int) --> void", pybind11::arg("x"), pybind11::arg("type"), pybind11::arg("iter"));
		cl.def("update", [](ROL::Constraint_Partitioned<double> &o, const class ROL::Vector<double> & a0) -> void { return o.update(a0); }, "", pybind11::arg("x"));
		cl.def("update", [](ROL::Constraint_Partitioned<double> &o, const class ROL::Vector<double> & a0, bool const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("flag"));
		cl.def("update", (void (ROL::Constraint_Partitioned<double>::*)(const class ROL::Vector<double> &, bool, int)) &ROL::Constraint_Partitioned<double>::update, "C++: ROL::Constraint_Partitioned<double>::update(const class ROL::Vector<double> &, bool, int) --> void", pybind11::arg("x"), pybind11::arg("flag"), pybind11::arg("iter"));
		cl.def("value", (void (ROL::Constraint_Partitioned<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Constraint_Partitioned<double>::value, "C++: ROL::Constraint_Partitioned<double>::value(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("c"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("applyJacobian", (void (ROL::Constraint_Partitioned<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Constraint_Partitioned<double>::applyJacobian, "C++: ROL::Constraint_Partitioned<double>::applyJacobian(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("jv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("applyAdjointJacobian", (void (ROL::Constraint_Partitioned<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Constraint_Partitioned<double>::applyAdjointJacobian, "C++: ROL::Constraint_Partitioned<double>::applyAdjointJacobian(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("ajv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("applyAdjointHessian", (void (ROL::Constraint_Partitioned<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Constraint_Partitioned<double>::applyAdjointHessian, "C++: ROL::Constraint_Partitioned<double>::applyAdjointHessian(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("ahuv"), pybind11::arg("u"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("applyPreconditioner", (void (ROL::Constraint_Partitioned<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Constraint_Partitioned<double>::applyPreconditioner, "C++: ROL::Constraint_Partitioned<double>::applyPreconditioner(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("pv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("tol"));
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
}
