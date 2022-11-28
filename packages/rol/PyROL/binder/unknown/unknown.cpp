#include <ROL_AffineTransformConstraint.hpp>
#include <ROL_AffineTransformObjective.hpp>
#include <ROL_BoundConstraint.hpp>
#include <ROL_BoundConstraint_Partitioned.hpp>
#include <ROL_ConjugateGradients.hpp>
#include <ROL_ConjugateResiduals.hpp>
#include <ROL_ConstraintAssembler.hpp>
#include <ROL_Constraint_Partitioned.hpp>
#include <ROL_DaiFletcherProjection.hpp>
#include <ROL_DouglasRachfordProjection.hpp>
#include <ROL_DykstraProjection.hpp>
#include <ROL_Elementwise_Function.hpp>
#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_GMRES.hpp>
#include <ROL_Krylov.hpp>
#include <ROL_LinearConstraint.hpp>
#include <ROL_LinearOperator.hpp>
#include <ROL_NullSpaceOperator.hpp>
#include <ROL_PartitionedVector.hpp>
#include <ROL_PolyhedralProjection.hpp>
#include <ROL_Problem.hpp>
#include <ROL_RangeSpaceOperator.hpp>
#include <ROL_ReduceLinearConstraint.hpp>
#include <ROL_SlacklessObjective.hpp>
#include <ROL_StatusTest.hpp>
#include <ROL_Types.hpp>
#include <ROL_UpdateType.hpp>
#include <ROL_VectorController.hpp>
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

// ROL::Vector file: line:14
struct PyCallBack_ROL_Vector_double_t : public ROL::Vector<double> {
	using ROL::Vector<double>::Vector;

	void plus(const class ROL::Vector<double> & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Vector<double> *>(this), "plus");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"Vector::plus\"");
	}
	void scale(const double a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Vector<double> *>(this), "scale");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"Vector::scale\"");
	}
	double dot(const class ROL::Vector<double> & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Vector<double> *>(this), "dot");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"Vector::dot\"");
	}
	double norm() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Vector<double> *>(this), "norm");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"Vector::norm\"");
	}
	class Teuchos::RCP<class ROL::Vector<double> > clone() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Vector<double> *>(this), "clone");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<class Teuchos::RCP<class ROL::Vector<double> >>::value) {
				static pybind11::detail::override_caster_t<class Teuchos::RCP<class ROL::Vector<double> >> caster;
				return pybind11::detail::cast_ref<class Teuchos::RCP<class ROL::Vector<double> >>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Teuchos::RCP<class ROL::Vector<double> >>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"Vector::clone\"");
	}
	void axpy(const double a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Vector<double> *>(this), "axpy");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Vector::axpy(a0, a1);
	}
	void zero() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Vector<double> *>(this), "zero");
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
	class Teuchos::RCP<class ROL::Vector<double> > basis(const int a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Vector<double> *>(this), "basis");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<class Teuchos::RCP<class ROL::Vector<double> >>::value) {
				static pybind11::detail::override_caster_t<class Teuchos::RCP<class ROL::Vector<double> >> caster;
				return pybind11::detail::cast_ref<class Teuchos::RCP<class ROL::Vector<double> >>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Teuchos::RCP<class ROL::Vector<double> >>(std::move(o));
		}
		return Vector::basis(a0);
	}
	int dimension() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Vector<double> *>(this), "dimension");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return Vector::dimension();
	}
	void set(const class ROL::Vector<double> & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Vector<double> *>(this), "set");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Vector::set(a0);
	}
	const class ROL::Vector<double> & dual() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Vector<double> *>(this), "dual");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Vector<double> *>(this), "apply");
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
	void applyUnary(const class ROL::Elementwise::UnaryFunction<double> & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Vector<double> *>(this), "applyUnary");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Vector::applyUnary(a0);
	}
	void applyBinary(const class ROL::Elementwise::BinaryFunction<double> & a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Vector<double> *>(this), "applyBinary");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Vector::applyBinary(a0, a1);
	}
	double reduce(const class ROL::Elementwise::ReductionOp<double> & a0) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Vector<double> *>(this), "reduce");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return Vector::reduce(a0);
	}
	void setScalar(const double a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Vector<double> *>(this), "setScalar");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Vector::setScalar(a0);
	}
	void randomize(const double a0, const double a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Vector<double> *>(this), "randomize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Vector::randomize(a0, a1);
	}
};

// ROL::Objective file: line:17
struct PyCallBack_ROL_Objective_double_t : public ROL::Objective<double> {
	using ROL::Objective<double>::Objective;

	void update(const class ROL::Vector<double> & a0, enum ROL::UpdateType a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Objective<double> *>(this), "update");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Objective<double> *>(this), "update");
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
	double value(const class ROL::Vector<double> & a0, double & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Objective<double> *>(this), "value");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"Objective::value\"");
	}
	void gradient(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double & a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Objective<double> *>(this), "gradient");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Objective::gradient(a0, a1, a2);
	}
	double dirDeriv(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double & a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Objective<double> *>(this), "dirDeriv");
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
	void hessVec(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Objective<double> *>(this), "hessVec");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Objective::hessVec(a0, a1, a2, a3);
	}
	void invHessVec(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Objective<double> *>(this), "invHessVec");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Objective<double> *>(this), "precond");
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

// ROL::Constraint file: line:21
struct PyCallBack_ROL_Constraint_double_t : public ROL::Constraint<double> {
	using ROL::Constraint<double>::Constraint;

	void update(const class ROL::Vector<double> & a0, enum ROL::UpdateType a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint::update(a0, a1, a2);
	}
	void update(const class ROL::Vector<double> & a0, bool a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint::update(a0, a1, a2);
	}
	void value(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double & a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint<double> *>(this), "value");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"Constraint::value\"");
	}
	void applyJacobian(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint<double> *>(this), "applyJacobian");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint::applyJacobian(a0, a1, a2, a3);
	}
	void applyAdjointJacobian(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint<double> *>(this), "applyAdjointJacobian");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint::applyAdjointJacobian(a0, a1, a2, a3);
	}
	void applyAdjointJacobian(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, double & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint<double> *>(this), "applyAdjointJacobian");
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
	void applyAdjointHessian(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, double & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint<double> *>(this), "applyAdjointHessian");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint::applyAdjointHessian(a0, a1, a2, a3, a4);
	}
	void applyPreconditioner(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, double & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint<double> *>(this), "applyPreconditioner");
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

// ROL::Problem file: line:29
struct PyCallBack_ROL_Problem_double_t : public ROL::Problem<double> {
	using ROL::Problem<double>::Problem;

	void edit() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Problem<double> *>(this), "edit");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Problem::edit();
	}
};

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

// ROL::AffineTransformObjective file:ROL_AffineTransformObjective.hpp line:62
struct PyCallBack_ROL_AffineTransformObjective_double_t : public ROL::AffineTransformObjective<double> {
	using ROL::AffineTransformObjective<double>::AffineTransformObjective;

	void update(const class ROL::Vector<double> & a0, enum ROL::UpdateType a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::AffineTransformObjective<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return AffineTransformObjective::update(a0, a1, a2);
	}
	void update(const class ROL::Vector<double> & a0, bool a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::AffineTransformObjective<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return AffineTransformObjective::update(a0, a1, a2);
	}
	double value(const class ROL::Vector<double> & a0, double & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::AffineTransformObjective<double> *>(this), "value");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return AffineTransformObjective::value(a0, a1);
	}
	void gradient(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double & a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::AffineTransformObjective<double> *>(this), "gradient");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return AffineTransformObjective::gradient(a0, a1, a2);
	}
	void hessVec(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::AffineTransformObjective<double> *>(this), "hessVec");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return AffineTransformObjective::hessVec(a0, a1, a2, a3);
	}
	double dirDeriv(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double & a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::AffineTransformObjective<double> *>(this), "dirDeriv");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::AffineTransformObjective<double> *>(this), "invHessVec");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::AffineTransformObjective<double> *>(this), "precond");
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

// ROL::AffineTransformConstraint file:ROL_AffineTransformConstraint.hpp line:62
struct PyCallBack_ROL_AffineTransformConstraint_double_t : public ROL::AffineTransformConstraint<double> {
	using ROL::AffineTransformConstraint<double>::AffineTransformConstraint;

	void update(const class ROL::Vector<double> & a0, enum ROL::UpdateType a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::AffineTransformConstraint<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return AffineTransformConstraint::update(a0, a1, a2);
	}
	void update(const class ROL::Vector<double> & a0, bool a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::AffineTransformConstraint<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return AffineTransformConstraint::update(a0, a1, a2);
	}
	void value(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double & a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::AffineTransformConstraint<double> *>(this), "value");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return AffineTransformConstraint::value(a0, a1, a2);
	}
	void applyJacobian(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::AffineTransformConstraint<double> *>(this), "applyJacobian");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return AffineTransformConstraint::applyJacobian(a0, a1, a2, a3);
	}
	void applyAdjointJacobian(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::AffineTransformConstraint<double> *>(this), "applyAdjointJacobian");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return AffineTransformConstraint::applyAdjointJacobian(a0, a1, a2, a3);
	}
	void applyAdjointHessian(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, double & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::AffineTransformConstraint<double> *>(this), "applyAdjointHessian");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return AffineTransformConstraint::applyAdjointHessian(a0, a1, a2, a3, a4);
	}
	void applyAdjointJacobian(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, double & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::AffineTransformConstraint<double> *>(this), "applyAdjointJacobian");
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
	void applyPreconditioner(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, double & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::AffineTransformConstraint<double> *>(this), "applyPreconditioner");
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

// ROL::NullSpaceOperator file:ROL_NullSpaceOperator.hpp line:60
struct PyCallBack_ROL_NullSpaceOperator_double_t : public ROL::NullSpaceOperator<double> {
	using ROL::NullSpaceOperator<double>::NullSpaceOperator;

	void update(const class ROL::Vector<double> & a0, bool a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::NullSpaceOperator<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NullSpaceOperator::update(a0, a1, a2);
	}
	void apply(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double & a2) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::NullSpaceOperator<double> *>(this), "apply");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NullSpaceOperator::apply(a0, a1, a2);
	}
	void applyAdjoint(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double & a2) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::NullSpaceOperator<double> *>(this), "applyAdjoint");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NullSpaceOperator::applyAdjoint(a0, a1, a2);
	}
	void applyInverse(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double & a2) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::NullSpaceOperator<double> *>(this), "applyInverse");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NullSpaceOperator::applyInverse(a0, a1, a2);
	}
	void applyAdjointInverse(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double & a2) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::NullSpaceOperator<double> *>(this), "applyAdjointInverse");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NullSpaceOperator::applyAdjointInverse(a0, a1, a2);
	}
};

// ROL::PolyhedralProjection file:ROL_PolyhedralProjection.hpp line:55
struct PyCallBack_ROL_PolyhedralProjection_double_t : public ROL::PolyhedralProjection<double> {
	using ROL::PolyhedralProjection<double>::PolyhedralProjection;

};

// ROL::DaiFletcherProjection file:ROL_DaiFletcherProjection.hpp line:54
struct PyCallBack_ROL_DaiFletcherProjection_double_t : public ROL::DaiFletcherProjection<double> {
	using ROL::DaiFletcherProjection<double>::DaiFletcherProjection;

};

// ROL::DykstraProjection file:ROL_DykstraProjection.hpp line:54
struct PyCallBack_ROL_DykstraProjection_double_t : public ROL::DykstraProjection<double> {
	using ROL::DykstraProjection<double>::DykstraProjection;

};

// ROL::DouglasRachfordProjection file:ROL_DouglasRachfordProjection.hpp line:54
struct PyCallBack_ROL_DouglasRachfordProjection_double_t : public ROL::DouglasRachfordProjection<double> {
	using ROL::DouglasRachfordProjection<double>::DouglasRachfordProjection;

};

// ROL::Krylov file:ROL_Krylov.hpp line:58
struct PyCallBack_ROL_Krylov_double_t : public ROL::Krylov<double> {
	using ROL::Krylov<double>::Krylov;

	double run(class ROL::Vector<double> & a0, class ROL::LinearOperator<double> & a1, const class ROL::Vector<double> & a2, class ROL::LinearOperator<double> & a3, int & a4, int & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Krylov<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"Krylov::run\"");
	}
};

// ROL::ConjugateGradients file:ROL_ConjugateGradients.hpp line:57
struct PyCallBack_ROL_ConjugateGradients_double_t : public ROL::ConjugateGradients<double> {
	using ROL::ConjugateGradients<double>::ConjugateGradients;

	double run(class ROL::Vector<double> & a0, class ROL::LinearOperator<double> & a1, const class ROL::Vector<double> & a2, class ROL::LinearOperator<double> & a3, int & a4, int & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::ConjugateGradients<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return ConjugateGradients::run(a0, a1, a2, a3, a4, a5);
	}
};

// ROL::ConjugateResiduals file:ROL_ConjugateResiduals.hpp line:57
struct PyCallBack_ROL_ConjugateResiduals_double_t : public ROL::ConjugateResiduals<double> {
	using ROL::ConjugateResiduals<double>::ConjugateResiduals;

	double run(class ROL::Vector<double> & a0, class ROL::LinearOperator<double> & a1, const class ROL::Vector<double> & a2, class ROL::LinearOperator<double> & a3, int & a4, int & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::ConjugateResiduals<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return ConjugateResiduals::run(a0, a1, a2, a3, a4, a5);
	}
};

// ROL::GMRES file:ROL_GMRES.hpp line:60
struct PyCallBack_ROL_GMRES_double_t : public ROL::GMRES<double> {
	using ROL::GMRES<double>::GMRES;

	double run(class ROL::Vector<double> & a0, class ROL::LinearOperator<double> & a1, const class ROL::Vector<double> & a2, class ROL::LinearOperator<double> & a3, int & a4, int & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::GMRES<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return GMRES::run(a0, a1, a2, a3, a4, a5);
	}
};

void bind_unknown_unknown(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // ROL::Vector file: line:14
		pybind11::class_<ROL::Vector<double>, Teuchos::RCP<ROL::Vector<double>>, PyCallBack_ROL_Vector_double_t> cl(M("ROL"), "Vector_double_t", "");
		cl.def( pybind11::init( [](){ return new PyCallBack_ROL_Vector_double_t(); } ) );
		cl.def(pybind11::init<PyCallBack_ROL_Vector_double_t const &>());
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
	{ // ROL::Objective file: line:17
		pybind11::class_<ROL::Objective<double>, Teuchos::RCP<ROL::Objective<double>>, PyCallBack_ROL_Objective_double_t> cl(M("ROL"), "Objective_double_t", "");
		cl.def( pybind11::init( [](){ return new PyCallBack_ROL_Objective_double_t(); } ) );
		cl.def(pybind11::init<PyCallBack_ROL_Objective_double_t const &>());
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
	{ // ROL::Constraint file: line:21
		pybind11::class_<ROL::Constraint<double>, Teuchos::RCP<ROL::Constraint<double>>, PyCallBack_ROL_Constraint_double_t> cl(M("ROL"), "Constraint_double_t", "");
		cl.def( pybind11::init( [](){ return new PyCallBack_ROL_Constraint_double_t(); } ) );
		cl.def(pybind11::init<PyCallBack_ROL_Constraint_double_t const &>());
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
	{ // ROL::Problem file: line:29
		pybind11::class_<ROL::Problem<double>, Teuchos::RCP<ROL::Problem<double>>, PyCallBack_ROL_Problem_double_t> cl(M("ROL"), "Problem_double_t", "");
		cl.def( pybind11::init( [](const class Teuchos::RCP<class ROL::Objective<double> > & a0, const class Teuchos::RCP<class ROL::Vector<double> > & a1){ return new ROL::Problem<double>(a0, a1); }, [](const class Teuchos::RCP<class ROL::Objective<double> > & a0, const class Teuchos::RCP<class ROL::Vector<double> > & a1){ return new PyCallBack_ROL_Problem_double_t(a0, a1); } ), "doc");
		cl.def( pybind11::init<const class Teuchos::RCP<class ROL::Objective<double> > &, const class Teuchos::RCP<class ROL::Vector<double> > &, const class Teuchos::RCP<class ROL::Vector<double> > &>(), pybind11::arg("obj"), pybind11::arg("x"), pybind11::arg("g") );

		cl.def( pybind11::init( [](PyCallBack_ROL_Problem_double_t const &o){ return new PyCallBack_ROL_Problem_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::Problem<double> const &o){ return new ROL::Problem<double>(o); } ) );
		cl.def("addBoundConstraint", (void (ROL::Problem<double>::*)(const class Teuchos::RCP<class ROL::BoundConstraint<double> > &)) &ROL::Problem<double>::addBoundConstraint, "Add a bound constraint.\n\n      \n  bound constraint object\n\nC++: ROL::Problem<double>::addBoundConstraint(const class Teuchos::RCP<class ROL::BoundConstraint<double> > &) --> void", pybind11::arg("bnd"));
		cl.def("removeBoundConstraint", (void (ROL::Problem<double>::*)()) &ROL::Problem<double>::removeBoundConstraint, "Remove an existing bound constraint.\n\nC++: ROL::Problem<double>::removeBoundConstraint() --> void");
		cl.def("addConstraint", [](ROL::Problem<double> &o, std::string const & a0, const class Teuchos::RCP<class ROL::Constraint<double> > & a1, const class Teuchos::RCP<class ROL::Vector<double> > & a2) -> void { return o.addConstraint(a0, a1, a2); }, "", pybind11::arg("name"), pybind11::arg("econ"), pybind11::arg("emul"));
		cl.def("addConstraint", [](ROL::Problem<double> &o, std::string const & a0, const class Teuchos::RCP<class ROL::Constraint<double> > & a1, const class Teuchos::RCP<class ROL::Vector<double> > & a2, const class Teuchos::RCP<class ROL::Vector<double> > & a3) -> void { return o.addConstraint(a0, a1, a2, a3); }, "", pybind11::arg("name"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"));
		cl.def("addConstraint", (void (ROL::Problem<double>::*)(std::string, const class Teuchos::RCP<class ROL::Constraint<double> > &, const class Teuchos::RCP<class ROL::Vector<double> > &, const class Teuchos::RCP<class ROL::Vector<double> > &, bool)) &ROL::Problem<double>::addConstraint, "Add an equality constraint.\n\n      \n   the unique constraint identifier\n      \n\n   constraint object\n      \n\n   dual constraint space vector\n      \n\n   primal constraint space vector\n      \n\n  whether or not to clear constraint container\n\nC++: ROL::Problem<double>::addConstraint(std::string, const class Teuchos::RCP<class ROL::Constraint<double> > &, const class Teuchos::RCP<class ROL::Vector<double> > &, const class Teuchos::RCP<class ROL::Vector<double> > &, bool) --> void", pybind11::arg("name"), pybind11::arg("econ"), pybind11::arg("emul"), pybind11::arg("eres"), pybind11::arg("reset"));
		cl.def("addConstraint", [](ROL::Problem<double> &o, std::string const & a0, const class Teuchos::RCP<class ROL::Constraint<double> > & a1, const class Teuchos::RCP<class ROL::Vector<double> > & a2, const class Teuchos::RCP<class ROL::BoundConstraint<double> > & a3) -> void { return o.addConstraint(a0, a1, a2, a3); }, "", pybind11::arg("name"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"));
		cl.def("addConstraint", [](ROL::Problem<double> &o, std::string const & a0, const class Teuchos::RCP<class ROL::Constraint<double> > & a1, const class Teuchos::RCP<class ROL::Vector<double> > & a2, const class Teuchos::RCP<class ROL::BoundConstraint<double> > & a3, const class Teuchos::RCP<class ROL::Vector<double> > & a4) -> void { return o.addConstraint(a0, a1, a2, a3, a4); }, "", pybind11::arg("name"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"));
		cl.def("addConstraint", (void (ROL::Problem<double>::*)(std::string, const class Teuchos::RCP<class ROL::Constraint<double> > &, const class Teuchos::RCP<class ROL::Vector<double> > &, const class Teuchos::RCP<class ROL::BoundConstraint<double> > &, const class Teuchos::RCP<class ROL::Vector<double> > &, bool)) &ROL::Problem<double>::addConstraint, "Add an inequality constraint.\n\n      \n   the unique constraint identifier\n      \n\n   constraint object\n      \n\n   dual constraint space vector\n      \n\n   bound constraint\n      \n\n   primal constraint space vector\n      \n\n  whether or not to clear constraint container\n\nC++: ROL::Problem<double>::addConstraint(std::string, const class Teuchos::RCP<class ROL::Constraint<double> > &, const class Teuchos::RCP<class ROL::Vector<double> > &, const class Teuchos::RCP<class ROL::BoundConstraint<double> > &, const class Teuchos::RCP<class ROL::Vector<double> > &, bool) --> void", pybind11::arg("name"), pybind11::arg("icon"), pybind11::arg("imul"), pybind11::arg("ibnd"), pybind11::arg("ires"), pybind11::arg("reset"));
		cl.def("removeConstraint", (void (ROL::Problem<double>::*)(std::string)) &ROL::Problem<double>::removeConstraint, "Remove an existing constraint.\n\n      \n  the unique constraint identifier\n\nC++: ROL::Problem<double>::removeConstraint(std::string) --> void", pybind11::arg("name"));
		cl.def("addLinearConstraint", [](ROL::Problem<double> &o, std::string const & a0, const class Teuchos::RCP<class ROL::Constraint<double> > & a1, const class Teuchos::RCP<class ROL::Vector<double> > & a2) -> void { return o.addLinearConstraint(a0, a1, a2); }, "", pybind11::arg("name"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"));
		cl.def("addLinearConstraint", [](ROL::Problem<double> &o, std::string const & a0, const class Teuchos::RCP<class ROL::Constraint<double> > & a1, const class Teuchos::RCP<class ROL::Vector<double> > & a2, const class Teuchos::RCP<class ROL::Vector<double> > & a3) -> void { return o.addLinearConstraint(a0, a1, a2, a3); }, "", pybind11::arg("name"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"));
		cl.def("addLinearConstraint", (void (ROL::Problem<double>::*)(std::string, const class Teuchos::RCP<class ROL::Constraint<double> > &, const class Teuchos::RCP<class ROL::Vector<double> > &, const class Teuchos::RCP<class ROL::Vector<double> > &, bool)) &ROL::Problem<double>::addLinearConstraint, "Add a linear equality constraint.\n\n      \n         the unique constraint identifier\n      \n\n  constraint object\n      \n\n  dual constraint space vector\n      \n\n  primal constraint space vector\n      \n\n        whether or not to clear linear constraint container\n\nC++: ROL::Problem<double>::addLinearConstraint(std::string, const class Teuchos::RCP<class ROL::Constraint<double> > &, const class Teuchos::RCP<class ROL::Vector<double> > &, const class Teuchos::RCP<class ROL::Vector<double> > &, bool) --> void", pybind11::arg("name"), pybind11::arg("linear_econ"), pybind11::arg("linear_emul"), pybind11::arg("linear_eres"), pybind11::arg("reset"));
		cl.def("addLinearConstraint", [](ROL::Problem<double> &o, std::string const & a0, const class Teuchos::RCP<class ROL::Constraint<double> > & a1, const class Teuchos::RCP<class ROL::Vector<double> > & a2, const class Teuchos::RCP<class ROL::BoundConstraint<double> > & a3) -> void { return o.addLinearConstraint(a0, a1, a2, a3); }, "", pybind11::arg("name"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"));
		cl.def("addLinearConstraint", [](ROL::Problem<double> &o, std::string const & a0, const class Teuchos::RCP<class ROL::Constraint<double> > & a1, const class Teuchos::RCP<class ROL::Vector<double> > & a2, const class Teuchos::RCP<class ROL::BoundConstraint<double> > & a3, const class Teuchos::RCP<class ROL::Vector<double> > & a4) -> void { return o.addLinearConstraint(a0, a1, a2, a3, a4); }, "", pybind11::arg("name"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("linear_ires"));
		cl.def("addLinearConstraint", (void (ROL::Problem<double>::*)(std::string, const class Teuchos::RCP<class ROL::Constraint<double> > &, const class Teuchos::RCP<class ROL::Vector<double> > &, const class Teuchos::RCP<class ROL::BoundConstraint<double> > &, const class Teuchos::RCP<class ROL::Vector<double> > &, bool)) &ROL::Problem<double>::addLinearConstraint, "Add a linear inequality constraint.\n\n      \n         the unique constraint identifier\n      \n\n  constraint object\n      \n\n  dual constraint space vector\n      \n\n  bound constraint\n      \n\n  primal constraint space vector\n      \n\n        whether or not to clear linear constraint container\n\nC++: ROL::Problem<double>::addLinearConstraint(std::string, const class Teuchos::RCP<class ROL::Constraint<double> > &, const class Teuchos::RCP<class ROL::Vector<double> > &, const class Teuchos::RCP<class ROL::BoundConstraint<double> > &, const class Teuchos::RCP<class ROL::Vector<double> > &, bool) --> void", pybind11::arg("name"), pybind11::arg("linear_icon"), pybind11::arg("linear_imul"), pybind11::arg("linear_ibnd"), pybind11::arg("linear_ires"), pybind11::arg("reset"));
		cl.def("removeLinearConstraint", (void (ROL::Problem<double>::*)(std::string)) &ROL::Problem<double>::removeLinearConstraint, "Remove an existing linear constraint.\n\n      \n  the unique constraint identifier\n\nC++: ROL::Problem<double>::removeLinearConstraint(std::string) --> void", pybind11::arg("name"));
		cl.def("setProjectionAlgorithm", (void (ROL::Problem<double>::*)(class Teuchos::ParameterList &)) &ROL::Problem<double>::setProjectionAlgorithm, "Set polyhedral projection algorithm.\n\n      \n  polyhedral projection algorithm\n\nC++: ROL::Problem<double>::setProjectionAlgorithm(class Teuchos::ParameterList &) --> void", pybind11::arg("list"));
		cl.def("getObjective", (const class Teuchos::RCP<class ROL::Objective<double> > & (ROL::Problem<double>::*)()) &ROL::Problem<double>::getObjective, "Get the objective function.\n\nC++: ROL::Problem<double>::getObjective() --> const class Teuchos::RCP<class ROL::Objective<double> > &", pybind11::return_value_policy::automatic);
		cl.def("getPrimalOptimizationVector", (const class Teuchos::RCP<class ROL::Vector<double> > & (ROL::Problem<double>::*)()) &ROL::Problem<double>::getPrimalOptimizationVector, "Get the primal optimization space vector.\n\nC++: ROL::Problem<double>::getPrimalOptimizationVector() --> const class Teuchos::RCP<class ROL::Vector<double> > &", pybind11::return_value_policy::automatic);
		cl.def("getDualOptimizationVector", (const class Teuchos::RCP<class ROL::Vector<double> > & (ROL::Problem<double>::*)()) &ROL::Problem<double>::getDualOptimizationVector, "Get the dual optimization space vector.\n\nC++: ROL::Problem<double>::getDualOptimizationVector() --> const class Teuchos::RCP<class ROL::Vector<double> > &", pybind11::return_value_policy::automatic);
		cl.def("getBoundConstraint", (const class Teuchos::RCP<class ROL::BoundConstraint<double> > & (ROL::Problem<double>::*)()) &ROL::Problem<double>::getBoundConstraint, "Get the bound constraint.\n\nC++: ROL::Problem<double>::getBoundConstraint() --> const class Teuchos::RCP<class ROL::BoundConstraint<double> > &", pybind11::return_value_policy::automatic);
		cl.def("getConstraint", (const class Teuchos::RCP<class ROL::Constraint<double> > & (ROL::Problem<double>::*)()) &ROL::Problem<double>::getConstraint, "Get the equality constraint.\n\nC++: ROL::Problem<double>::getConstraint() --> const class Teuchos::RCP<class ROL::Constraint<double> > &", pybind11::return_value_policy::automatic);
		cl.def("getMultiplierVector", (const class Teuchos::RCP<class ROL::Vector<double> > & (ROL::Problem<double>::*)()) &ROL::Problem<double>::getMultiplierVector, "Get the dual constraint space vector.\n\nC++: ROL::Problem<double>::getMultiplierVector() --> const class Teuchos::RCP<class ROL::Vector<double> > &", pybind11::return_value_policy::automatic);
		cl.def("getResidualVector", (const class Teuchos::RCP<class ROL::Vector<double> > & (ROL::Problem<double>::*)()) &ROL::Problem<double>::getResidualVector, "Get the primal constraint space vector.\n\nC++: ROL::Problem<double>::getResidualVector() --> const class Teuchos::RCP<class ROL::Vector<double> > &", pybind11::return_value_policy::automatic);
		cl.def("getPolyhedralProjection", (const class Teuchos::RCP<class ROL::PolyhedralProjection<double> > & (ROL::Problem<double>::*)()) &ROL::Problem<double>::getPolyhedralProjection, "Get the polyhedral projection object.  This is a null pointer if\n             no linear constraints and/or bounds are present.\n\nC++: ROL::Problem<double>::getPolyhedralProjection() --> const class Teuchos::RCP<class ROL::PolyhedralProjection<double> > &", pybind11::return_value_policy::automatic);
		cl.def("getProblemType", (enum ROL::EProblem (ROL::Problem<double>::*)()) &ROL::Problem<double>::getProblemType, "Get the optimization problem type (U, B, E, or G).\n\nC++: ROL::Problem<double>::getProblemType() --> enum ROL::EProblem");
		cl.def("isFinalized", (bool (ROL::Problem<double>::*)() const) &ROL::Problem<double>::isFinalized, "Indicate whether or no finalize has been called.\n\nC++: ROL::Problem<double>::isFinalized() const --> bool");
		cl.def("edit", (void (ROL::Problem<double>::*)()) &ROL::Problem<double>::edit, "Resume editting optimization problem after finalize has been called.\n\nC++: ROL::Problem<double>::edit() --> void");
		cl.def("finalizeIteration", (void (ROL::Problem<double>::*)()) &ROL::Problem<double>::finalizeIteration, "Transform the optimization variables to the native\n             parameterization after an optimization algorithm has finished.\n\nC++: ROL::Problem<double>::finalizeIteration() --> void");
		cl.def("assign", (class ROL::Problem<double> & (ROL::Problem<double>::*)(const class ROL::Problem<double> &)) &ROL::Problem<double>::operator=, "C++: ROL::Problem<double>::operator=(const class ROL::Problem<double> &) --> class ROL::Problem<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::PartitionedVector file:ROL_PartitionedVector.hpp line:60
		pybind11::class_<ROL::PartitionedVector<double>, Teuchos::RCP<ROL::PartitionedVector<double>>, PyCallBack_ROL_PartitionedVector_double_t, ROL::Vector<double>> cl(M("ROL"), "PartitionedVector_double_t", "");
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
		pybind11::class_<ROL::Constraint_Partitioned<double>, Teuchos::RCP<ROL::Constraint_Partitioned<double>>, PyCallBack_ROL_Constraint_Partitioned_double_t, ROL::Constraint<double>> cl(M("ROL"), "Constraint_Partitioned_double_t", "");
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
	{ // ROL::BoundConstraint file:ROL_BoundConstraint.hpp line:73
		pybind11::class_<ROL::BoundConstraint<double>, Teuchos::RCP<ROL::BoundConstraint<double>>, PyCallBack_ROL_BoundConstraint_double_t> cl(M("ROL"), "BoundConstraint_double_t", "");
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
		pybind11::class_<ROL::BoundConstraint_Partitioned<double>, Teuchos::RCP<ROL::BoundConstraint_Partitioned<double>>, PyCallBack_ROL_BoundConstraint_Partitioned_double_t, ROL::BoundConstraint<double>> cl(M("ROL"), "BoundConstraint_Partitioned_double_t", "");
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
	{ // ROL::SlacklessObjective file:ROL_SlacklessObjective.hpp line:59
		pybind11::class_<ROL::SlacklessObjective<double>, Teuchos::RCP<ROL::SlacklessObjective<double>>, PyCallBack_ROL_SlacklessObjective_double_t, ROL::Objective<double>> cl(M("ROL"), "SlacklessObjective_double_t", "");
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
		pybind11::class_<ROL::LinearConstraint<double>, Teuchos::RCP<ROL::LinearConstraint<double>>, PyCallBack_ROL_LinearConstraint_double_t, ROL::Constraint<double>> cl(M("ROL"), "LinearConstraint_double_t", "");
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
	{ // ROL::VectorController file:ROL_VectorController.hpp line:54
		pybind11::class_<ROL::VectorController<double,int>, Teuchos::RCP<ROL::VectorController<double,int>>> cl(M("ROL"), "VectorController_double_int_t", "");
		cl.def( pybind11::init( [](){ return new ROL::VectorController<double,int>(); } ) );
		cl.def( pybind11::init( [](ROL::VectorController<double,int> const &o){ return new ROL::VectorController<double,int>(o); } ) );
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
	{ // ROL::AffineTransformObjective file:ROL_AffineTransformObjective.hpp line:62
		pybind11::class_<ROL::AffineTransformObjective<double>, Teuchos::RCP<ROL::AffineTransformObjective<double>>, PyCallBack_ROL_AffineTransformObjective_double_t, ROL::Objective<double>> cl(M("ROL"), "AffineTransformObjective_double_t", "");
		cl.def( pybind11::init( [](PyCallBack_ROL_AffineTransformObjective_double_t const &o){ return new PyCallBack_ROL_AffineTransformObjective_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::AffineTransformObjective<double> const &o){ return new ROL::AffineTransformObjective<double>(o); } ) );
		cl.def("update", [](ROL::AffineTransformObjective<double> &o, const class ROL::Vector<double> & a0, enum ROL::UpdateType const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("type"));
		cl.def("update", (void (ROL::AffineTransformObjective<double>::*)(const class ROL::Vector<double> &, enum ROL::UpdateType, int)) &ROL::AffineTransformObjective<double>::update, "C++: ROL::AffineTransformObjective<double>::update(const class ROL::Vector<double> &, enum ROL::UpdateType, int) --> void", pybind11::arg("x"), pybind11::arg("type"), pybind11::arg("iter"));
		cl.def("update", [](ROL::AffineTransformObjective<double> &o, const class ROL::Vector<double> & a0) -> void { return o.update(a0); }, "", pybind11::arg("x"));
		cl.def("update", [](ROL::AffineTransformObjective<double> &o, const class ROL::Vector<double> & a0, bool const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("flag"));
		cl.def("update", (void (ROL::AffineTransformObjective<double>::*)(const class ROL::Vector<double> &, bool, int)) &ROL::AffineTransformObjective<double>::update, "C++: ROL::AffineTransformObjective<double>::update(const class ROL::Vector<double> &, bool, int) --> void", pybind11::arg("x"), pybind11::arg("flag"), pybind11::arg("iter"));
		cl.def("value", (double (ROL::AffineTransformObjective<double>::*)(const class ROL::Vector<double> &, double &)) &ROL::AffineTransformObjective<double>::value, "C++: ROL::AffineTransformObjective<double>::value(const class ROL::Vector<double> &, double &) --> double", pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("gradient", (void (ROL::AffineTransformObjective<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::AffineTransformObjective<double>::gradient, "C++: ROL::AffineTransformObjective<double>::gradient(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("hessVec", (void (ROL::AffineTransformObjective<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::AffineTransformObjective<double>::hessVec, "C++: ROL::AffineTransformObjective<double>::hessVec(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("hv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("tol"));
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
	{ // ROL::AffineTransformConstraint file:ROL_AffineTransformConstraint.hpp line:62
		pybind11::class_<ROL::AffineTransformConstraint<double>, Teuchos::RCP<ROL::AffineTransformConstraint<double>>, PyCallBack_ROL_AffineTransformConstraint_double_t, ROL::Constraint<double>> cl(M("ROL"), "AffineTransformConstraint_double_t", "");
		cl.def( pybind11::init( [](PyCallBack_ROL_AffineTransformConstraint_double_t const &o){ return new PyCallBack_ROL_AffineTransformConstraint_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::AffineTransformConstraint<double> const &o){ return new ROL::AffineTransformConstraint<double>(o); } ) );
		cl.def("update", [](ROL::AffineTransformConstraint<double> &o, const class ROL::Vector<double> & a0, enum ROL::UpdateType const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("type"));
		cl.def("update", (void (ROL::AffineTransformConstraint<double>::*)(const class ROL::Vector<double> &, enum ROL::UpdateType, int)) &ROL::AffineTransformConstraint<double>::update, "C++: ROL::AffineTransformConstraint<double>::update(const class ROL::Vector<double> &, enum ROL::UpdateType, int) --> void", pybind11::arg("x"), pybind11::arg("type"), pybind11::arg("iter"));
		cl.def("update", [](ROL::AffineTransformConstraint<double> &o, const class ROL::Vector<double> & a0) -> void { return o.update(a0); }, "", pybind11::arg("x"));
		cl.def("update", [](ROL::AffineTransformConstraint<double> &o, const class ROL::Vector<double> & a0, bool const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("flag"));
		cl.def("update", (void (ROL::AffineTransformConstraint<double>::*)(const class ROL::Vector<double> &, bool, int)) &ROL::AffineTransformConstraint<double>::update, "C++: ROL::AffineTransformConstraint<double>::update(const class ROL::Vector<double> &, bool, int) --> void", pybind11::arg("x"), pybind11::arg("flag"), pybind11::arg("iter"));
		cl.def("value", (void (ROL::AffineTransformConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::AffineTransformConstraint<double>::value, "C++: ROL::AffineTransformConstraint<double>::value(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("c"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("applyJacobian", (void (ROL::AffineTransformConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::AffineTransformConstraint<double>::applyJacobian, "C++: ROL::AffineTransformConstraint<double>::applyJacobian(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("jv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("applyAdjointJacobian", (void (ROL::AffineTransformConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::AffineTransformConstraint<double>::applyAdjointJacobian, "C++: ROL::AffineTransformConstraint<double>::applyAdjointJacobian(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("ajv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("applyAdjointHessian", (void (ROL::AffineTransformConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::AffineTransformConstraint<double>::applyAdjointHessian, "C++: ROL::AffineTransformConstraint<double>::applyAdjointHessian(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("ahuv"), pybind11::arg("u"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("tol"));
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
	{ // ROL::NullSpaceOperator file:ROL_NullSpaceOperator.hpp line:60
		pybind11::class_<ROL::NullSpaceOperator<double>, Teuchos::RCP<ROL::NullSpaceOperator<double>>, PyCallBack_ROL_NullSpaceOperator_double_t, ROL::LinearOperator<double>> cl(M("ROL"), "NullSpaceOperator_double_t", "");
		cl.def( pybind11::init( [](const class Teuchos::RCP<class ROL::Constraint<double> > & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2){ return new ROL::NullSpaceOperator<double>(a0, a1, a2); }, [](const class Teuchos::RCP<class ROL::Constraint<double> > & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2){ return new PyCallBack_ROL_NullSpaceOperator_double_t(a0, a1, a2); } ), "doc");
		cl.def( pybind11::init<const class Teuchos::RCP<class ROL::Constraint<double> > &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const bool>(), pybind11::arg("con"), pybind11::arg("dom"), pybind11::arg("ran"), pybind11::arg("useAugSys") );

		cl.def( pybind11::init<const class Teuchos::RCP<class ROL::Constraint<double> > &, const class Teuchos::RCP<const class ROL::Vector<double> > &, const class Teuchos::RCP<const class ROL::Vector<double> > &>(), pybind11::arg("con"), pybind11::arg("dom"), pybind11::arg("ran") );

		cl.def( pybind11::init( [](PyCallBack_ROL_NullSpaceOperator_double_t const &o){ return new PyCallBack_ROL_NullSpaceOperator_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::NullSpaceOperator<double> const &o){ return new ROL::NullSpaceOperator<double>(o); } ) );
		cl.def("update", [](ROL::NullSpaceOperator<double> &o, const class ROL::Vector<double> & a0) -> void { return o.update(a0); }, "", pybind11::arg("x"));
		cl.def("update", [](ROL::NullSpaceOperator<double> &o, const class ROL::Vector<double> & a0, bool const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("flag"));
		cl.def("update", (void (ROL::NullSpaceOperator<double>::*)(const class ROL::Vector<double> &, bool, int)) &ROL::NullSpaceOperator<double>::update, "C++: ROL::NullSpaceOperator<double>::update(const class ROL::Vector<double> &, bool, int) --> void", pybind11::arg("x"), pybind11::arg("flag"), pybind11::arg("iter"));
		cl.def("apply", (void (ROL::NullSpaceOperator<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) const) &ROL::NullSpaceOperator<double>::apply, "C++: ROL::NullSpaceOperator<double>::apply(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) const --> void", pybind11::arg("Hv"), pybind11::arg("v"), pybind11::arg("tol"));
		cl.def("applyAdjoint", (void (ROL::NullSpaceOperator<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) const) &ROL::NullSpaceOperator<double>::applyAdjoint, "C++: ROL::NullSpaceOperator<double>::applyAdjoint(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) const --> void", pybind11::arg("Hv"), pybind11::arg("v"), pybind11::arg("tol"));
		cl.def("applyInverse", (void (ROL::NullSpaceOperator<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) const) &ROL::NullSpaceOperator<double>::applyInverse, "C++: ROL::NullSpaceOperator<double>::applyInverse(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) const --> void", pybind11::arg("Hv"), pybind11::arg("v"), pybind11::arg("tol"));
		cl.def("applyAdjointInverse", (void (ROL::NullSpaceOperator<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) const) &ROL::NullSpaceOperator<double>::applyAdjointInverse, "C++: ROL::NullSpaceOperator<double>::applyAdjointInverse(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) const --> void", pybind11::arg("Hv"), pybind11::arg("v"), pybind11::arg("tol"));
		cl.def("update", [](ROL::LinearOperator<double> &o, const class ROL::Vector<double> & a0) -> void { return o.update(a0); }, "", pybind11::arg("x"));
		cl.def("update", [](ROL::LinearOperator<double> &o, const class ROL::Vector<double> & a0, bool const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("flag"));
		cl.def("update", (void (ROL::LinearOperator<double>::*)(const class ROL::Vector<double> &, bool, int)) &ROL::LinearOperator<double>::update, "C++: ROL::LinearOperator<double>::update(const class ROL::Vector<double> &, bool, int) --> void", pybind11::arg("x"), pybind11::arg("flag"), pybind11::arg("iter"));
		cl.def("apply", (void (ROL::LinearOperator<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) const) &ROL::LinearOperator<double>::apply, "C++: ROL::LinearOperator<double>::apply(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) const --> void", pybind11::arg("Hv"), pybind11::arg("v"), pybind11::arg("tol"));
		cl.def("applyInverse", (void (ROL::LinearOperator<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) const) &ROL::LinearOperator<double>::applyInverse, "C++: ROL::LinearOperator<double>::applyInverse(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) const --> void", pybind11::arg("Hv"), pybind11::arg("v"), pybind11::arg("tol"));
		cl.def("applyAdjoint", (void (ROL::LinearOperator<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) const) &ROL::LinearOperator<double>::applyAdjoint, "C++: ROL::LinearOperator<double>::applyAdjoint(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) const --> void", pybind11::arg("Hv"), pybind11::arg("v"), pybind11::arg("tol"));
		cl.def("applyAdjointInverse", (void (ROL::LinearOperator<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) const) &ROL::LinearOperator<double>::applyAdjointInverse, "C++: ROL::LinearOperator<double>::applyAdjointInverse(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) const --> void", pybind11::arg("Hv"), pybind11::arg("v"), pybind11::arg("tol"));
		cl.def("assign", (class ROL::LinearOperator<double> & (ROL::LinearOperator<double>::*)(const class ROL::LinearOperator<double> &)) &ROL::LinearOperator<double>::operator=, "C++: ROL::LinearOperator<double>::operator=(const class ROL::LinearOperator<double> &) --> class ROL::LinearOperator<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::ReduceLinearConstraint file:ROL_ReduceLinearConstraint.hpp line:63
		pybind11::class_<ROL::ReduceLinearConstraint<double>, Teuchos::RCP<ROL::ReduceLinearConstraint<double>>> cl(M("ROL"), "ReduceLinearConstraint_double_t", "");
		cl.def( pybind11::init<const class Teuchos::RCP<class ROL::Constraint<double> > &, const class Teuchos::RCP<class ROL::Vector<double> > &, const class Teuchos::RCP<const class ROL::Vector<double> > &>(), pybind11::arg("lcon"), pybind11::arg("x"), pybind11::arg("c") );

		cl.def( pybind11::init( [](ROL::ReduceLinearConstraint<double> const &o){ return new ROL::ReduceLinearConstraint<double>(o); } ) );
		cl.def("transform", (class Teuchos::RCP<class ROL::Objective<double> > (ROL::ReduceLinearConstraint<double>::*)(const class Teuchos::RCP<class ROL::Objective<double> > &) const) &ROL::ReduceLinearConstraint<double>::transform, "C++: ROL::ReduceLinearConstraint<double>::transform(const class Teuchos::RCP<class ROL::Objective<double> > &) const --> class Teuchos::RCP<class ROL::Objective<double> >", pybind11::arg("obj"));
		cl.def("transform", (class Teuchos::RCP<class ROL::Constraint<double> > (ROL::ReduceLinearConstraint<double>::*)(const class Teuchos::RCP<class ROL::Constraint<double> > &) const) &ROL::ReduceLinearConstraint<double>::transform, "C++: ROL::ReduceLinearConstraint<double>::transform(const class Teuchos::RCP<class ROL::Constraint<double> > &) const --> class Teuchos::RCP<class ROL::Constraint<double> >", pybind11::arg("con"));
		cl.def("getLinearConstraint", (class Teuchos::RCP<class ROL::Constraint<double> > (ROL::ReduceLinearConstraint<double>::*)() const) &ROL::ReduceLinearConstraint<double>::getLinearConstraint, "C++: ROL::ReduceLinearConstraint<double>::getLinearConstraint() const --> class Teuchos::RCP<class ROL::Constraint<double> >");
		cl.def("getFeasibleVector", (class Teuchos::RCP<const class ROL::Vector<double> > (ROL::ReduceLinearConstraint<double>::*)() const) &ROL::ReduceLinearConstraint<double>::getFeasibleVector, "C++: ROL::ReduceLinearConstraint<double>::getFeasibleVector() const --> class Teuchos::RCP<const class ROL::Vector<double> >");
		cl.def("project", (void (ROL::ReduceLinearConstraint<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &) const) &ROL::ReduceLinearConstraint<double>::project, "C++: ROL::ReduceLinearConstraint<double>::project(class ROL::Vector<double> &, const class ROL::Vector<double> &) const --> void", pybind11::arg("x"), pybind11::arg("y"));
		cl.def("project", (void (ROL::ReduceLinearConstraint<double>::*)(const class Teuchos::RCP<class ROL::Vector<double> > &, const class Teuchos::RCP<const class ROL::Vector<double> > &) const) &ROL::ReduceLinearConstraint<double>::project, "C++: ROL::ReduceLinearConstraint<double>::project(const class Teuchos::RCP<class ROL::Vector<double> > &, const class Teuchos::RCP<const class ROL::Vector<double> > &) const --> void", pybind11::arg("x"), pybind11::arg("y"));
	}
	{ // ROL::PolyhedralProjection file:ROL_PolyhedralProjection.hpp line:55
		pybind11::class_<ROL::PolyhedralProjection<double>, Teuchos::RCP<ROL::PolyhedralProjection<double>>, PyCallBack_ROL_PolyhedralProjection_double_t> cl(M("ROL"), "PolyhedralProjection_double_t", "");
		cl.def( pybind11::init<const class Teuchos::RCP<class ROL::BoundConstraint<double> > &>(), pybind11::arg("bnd") );

		cl.def( pybind11::init<const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class Teuchos::RCP<class ROL::BoundConstraint<double> > &, const class Teuchos::RCP<class ROL::Constraint<double> > &, const class ROL::Vector<double> &, const class ROL::Vector<double> &>(), pybind11::arg("xprim"), pybind11::arg("xdual"), pybind11::arg("bnd"), pybind11::arg("con"), pybind11::arg("mul"), pybind11::arg("res") );

		cl.def( pybind11::init( [](PyCallBack_ROL_PolyhedralProjection_double_t const &o){ return new PyCallBack_ROL_PolyhedralProjection_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::PolyhedralProjection<double> const &o){ return new ROL::PolyhedralProjection<double>(o); } ) );
		cl.def("getLinearConstraint", (const class Teuchos::RCP<class ROL::Constraint<double> > (ROL::PolyhedralProjection<double>::*)() const) &ROL::PolyhedralProjection<double>::getLinearConstraint, "C++: ROL::PolyhedralProjection<double>::getLinearConstraint() const --> const class Teuchos::RCP<class ROL::Constraint<double> >");
		cl.def("getBoundConstraint", (const class Teuchos::RCP<class ROL::BoundConstraint<double> > (ROL::PolyhedralProjection<double>::*)() const) &ROL::PolyhedralProjection<double>::getBoundConstraint, "C++: ROL::PolyhedralProjection<double>::getBoundConstraint() const --> const class Teuchos::RCP<class ROL::BoundConstraint<double> >");
		cl.def("getMultiplier", (const class Teuchos::RCP<class ROL::Vector<double> > (ROL::PolyhedralProjection<double>::*)() const) &ROL::PolyhedralProjection<double>::getMultiplier, "C++: ROL::PolyhedralProjection<double>::getMultiplier() const --> const class Teuchos::RCP<class ROL::Vector<double> >");
		cl.def("getResidual", (const class Teuchos::RCP<class ROL::Vector<double> > (ROL::PolyhedralProjection<double>::*)() const) &ROL::PolyhedralProjection<double>::getResidual, "C++: ROL::PolyhedralProjection<double>::getResidual() const --> const class Teuchos::RCP<class ROL::Vector<double> >");
	}
	{ // ROL::DaiFletcherProjection file:ROL_DaiFletcherProjection.hpp line:54
		pybind11::class_<ROL::DaiFletcherProjection<double>, Teuchos::RCP<ROL::DaiFletcherProjection<double>>, PyCallBack_ROL_DaiFletcherProjection_double_t, ROL::PolyhedralProjection<double>> cl(M("ROL"), "DaiFletcherProjection_double_t", "");
		cl.def( pybind11::init<const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class Teuchos::RCP<class ROL::BoundConstraint<double> > &, const class Teuchos::RCP<class ROL::Constraint<double> > &, const class ROL::Vector<double> &, const class ROL::Vector<double> &>(), pybind11::arg("xprim"), pybind11::arg("xdual"), pybind11::arg("bnd"), pybind11::arg("con"), pybind11::arg("mul"), pybind11::arg("res") );

		cl.def( pybind11::init<const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class Teuchos::RCP<class ROL::BoundConstraint<double> > &, const class Teuchos::RCP<class ROL::Constraint<double> > &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class Teuchos::ParameterList &>(), pybind11::arg("xprim"), pybind11::arg("xdual"), pybind11::arg("bnd"), pybind11::arg("con"), pybind11::arg("mul"), pybind11::arg("res"), pybind11::arg("list") );

		cl.def( pybind11::init( [](PyCallBack_ROL_DaiFletcherProjection_double_t const &o){ return new PyCallBack_ROL_DaiFletcherProjection_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::DaiFletcherProjection<double> const &o){ return new ROL::DaiFletcherProjection<double>(o); } ) );
		cl.def("getLinearConstraint", (const class Teuchos::RCP<class ROL::Constraint<double> > (ROL::PolyhedralProjection<double>::*)() const) &ROL::PolyhedralProjection<double>::getLinearConstraint, "C++: ROL::PolyhedralProjection<double>::getLinearConstraint() const --> const class Teuchos::RCP<class ROL::Constraint<double> >");
		cl.def("getBoundConstraint", (const class Teuchos::RCP<class ROL::BoundConstraint<double> > (ROL::PolyhedralProjection<double>::*)() const) &ROL::PolyhedralProjection<double>::getBoundConstraint, "C++: ROL::PolyhedralProjection<double>::getBoundConstraint() const --> const class Teuchos::RCP<class ROL::BoundConstraint<double> >");
		cl.def("getMultiplier", (const class Teuchos::RCP<class ROL::Vector<double> > (ROL::PolyhedralProjection<double>::*)() const) &ROL::PolyhedralProjection<double>::getMultiplier, "C++: ROL::PolyhedralProjection<double>::getMultiplier() const --> const class Teuchos::RCP<class ROL::Vector<double> >");
		cl.def("getResidual", (const class Teuchos::RCP<class ROL::Vector<double> > (ROL::PolyhedralProjection<double>::*)() const) &ROL::PolyhedralProjection<double>::getResidual, "C++: ROL::PolyhedralProjection<double>::getResidual() const --> const class Teuchos::RCP<class ROL::Vector<double> >");
	}
	{ // ROL::DykstraProjection file:ROL_DykstraProjection.hpp line:54
		pybind11::class_<ROL::DykstraProjection<double>, Teuchos::RCP<ROL::DykstraProjection<double>>, PyCallBack_ROL_DykstraProjection_double_t, ROL::PolyhedralProjection<double>> cl(M("ROL"), "DykstraProjection_double_t", "");
		cl.def( pybind11::init<const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class Teuchos::RCP<class ROL::BoundConstraint<double> > &, const class Teuchos::RCP<class ROL::Constraint<double> > &, const class ROL::Vector<double> &, const class ROL::Vector<double> &>(), pybind11::arg("xprim"), pybind11::arg("xdual"), pybind11::arg("bnd"), pybind11::arg("con"), pybind11::arg("mul"), pybind11::arg("res") );

		cl.def( pybind11::init<const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class Teuchos::RCP<class ROL::BoundConstraint<double> > &, const class Teuchos::RCP<class ROL::Constraint<double> > &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class Teuchos::ParameterList &>(), pybind11::arg("xprim"), pybind11::arg("xdual"), pybind11::arg("bnd"), pybind11::arg("con"), pybind11::arg("mul"), pybind11::arg("res"), pybind11::arg("list") );

		cl.def( pybind11::init( [](PyCallBack_ROL_DykstraProjection_double_t const &o){ return new PyCallBack_ROL_DykstraProjection_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::DykstraProjection<double> const &o){ return new ROL::DykstraProjection<double>(o); } ) );
		cl.def("getLinearConstraint", (const class Teuchos::RCP<class ROL::Constraint<double> > (ROL::PolyhedralProjection<double>::*)() const) &ROL::PolyhedralProjection<double>::getLinearConstraint, "C++: ROL::PolyhedralProjection<double>::getLinearConstraint() const --> const class Teuchos::RCP<class ROL::Constraint<double> >");
		cl.def("getBoundConstraint", (const class Teuchos::RCP<class ROL::BoundConstraint<double> > (ROL::PolyhedralProjection<double>::*)() const) &ROL::PolyhedralProjection<double>::getBoundConstraint, "C++: ROL::PolyhedralProjection<double>::getBoundConstraint() const --> const class Teuchos::RCP<class ROL::BoundConstraint<double> >");
		cl.def("getMultiplier", (const class Teuchos::RCP<class ROL::Vector<double> > (ROL::PolyhedralProjection<double>::*)() const) &ROL::PolyhedralProjection<double>::getMultiplier, "C++: ROL::PolyhedralProjection<double>::getMultiplier() const --> const class Teuchos::RCP<class ROL::Vector<double> >");
		cl.def("getResidual", (const class Teuchos::RCP<class ROL::Vector<double> > (ROL::PolyhedralProjection<double>::*)() const) &ROL::PolyhedralProjection<double>::getResidual, "C++: ROL::PolyhedralProjection<double>::getResidual() const --> const class Teuchos::RCP<class ROL::Vector<double> >");
	}
	{ // ROL::DouglasRachfordProjection file:ROL_DouglasRachfordProjection.hpp line:54
		pybind11::class_<ROL::DouglasRachfordProjection<double>, Teuchos::RCP<ROL::DouglasRachfordProjection<double>>, PyCallBack_ROL_DouglasRachfordProjection_double_t, ROL::PolyhedralProjection<double>> cl(M("ROL"), "DouglasRachfordProjection_double_t", "");
		cl.def( pybind11::init<const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class Teuchos::RCP<class ROL::BoundConstraint<double> > &, const class Teuchos::RCP<class ROL::Constraint<double> > &, const class ROL::Vector<double> &, const class ROL::Vector<double> &>(), pybind11::arg("xprim"), pybind11::arg("xdual"), pybind11::arg("bnd"), pybind11::arg("con"), pybind11::arg("mul"), pybind11::arg("res") );

		cl.def( pybind11::init<const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class Teuchos::RCP<class ROL::BoundConstraint<double> > &, const class Teuchos::RCP<class ROL::Constraint<double> > &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class Teuchos::ParameterList &>(), pybind11::arg("xprim"), pybind11::arg("xdual"), pybind11::arg("bnd"), pybind11::arg("con"), pybind11::arg("mul"), pybind11::arg("res"), pybind11::arg("list") );

		cl.def( pybind11::init( [](PyCallBack_ROL_DouglasRachfordProjection_double_t const &o){ return new PyCallBack_ROL_DouglasRachfordProjection_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::DouglasRachfordProjection<double> const &o){ return new ROL::DouglasRachfordProjection<double>(o); } ) );
		cl.def("getLinearConstraint", (const class Teuchos::RCP<class ROL::Constraint<double> > (ROL::PolyhedralProjection<double>::*)() const) &ROL::PolyhedralProjection<double>::getLinearConstraint, "C++: ROL::PolyhedralProjection<double>::getLinearConstraint() const --> const class Teuchos::RCP<class ROL::Constraint<double> >");
		cl.def("getBoundConstraint", (const class Teuchos::RCP<class ROL::BoundConstraint<double> > (ROL::PolyhedralProjection<double>::*)() const) &ROL::PolyhedralProjection<double>::getBoundConstraint, "C++: ROL::PolyhedralProjection<double>::getBoundConstraint() const --> const class Teuchos::RCP<class ROL::BoundConstraint<double> >");
		cl.def("getMultiplier", (const class Teuchos::RCP<class ROL::Vector<double> > (ROL::PolyhedralProjection<double>::*)() const) &ROL::PolyhedralProjection<double>::getMultiplier, "C++: ROL::PolyhedralProjection<double>::getMultiplier() const --> const class Teuchos::RCP<class ROL::Vector<double> >");
		cl.def("getResidual", (const class Teuchos::RCP<class ROL::Vector<double> > (ROL::PolyhedralProjection<double>::*)() const) &ROL::PolyhedralProjection<double>::getResidual, "C++: ROL::PolyhedralProjection<double>::getResidual() const --> const class Teuchos::RCP<class ROL::Vector<double> >");
	}
	{ // ROL::Krylov file:ROL_Krylov.hpp line:58
		pybind11::class_<ROL::Krylov<double>, Teuchos::RCP<ROL::Krylov<double>>, PyCallBack_ROL_Krylov_double_t> cl(M("ROL"), "Krylov_double_t", "");
		cl.def( pybind11::init( [](){ return new PyCallBack_ROL_Krylov_double_t(); } ), "doc");
		cl.def( pybind11::init( [](double const & a0){ return new PyCallBack_ROL_Krylov_double_t(a0); } ), "doc");
		cl.def( pybind11::init( [](double const & a0, double const & a1){ return new PyCallBack_ROL_Krylov_double_t(a0, a1); } ), "doc");
		cl.def( pybind11::init<double, double, unsigned int>(), pybind11::arg("absTol"), pybind11::arg("relTol"), pybind11::arg("maxit") );

		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("parlist") );

		cl.def(pybind11::init<PyCallBack_ROL_Krylov_double_t const &>());
		cl.def("run", (double (ROL::Krylov<double>::*)(class ROL::Vector<double> &, class ROL::LinearOperator<double> &, const class ROL::Vector<double> &, class ROL::LinearOperator<double> &, int &, int &)) &ROL::Krylov<double>::run, "C++: ROL::Krylov<double>::run(class ROL::Vector<double> &, class ROL::LinearOperator<double> &, const class ROL::Vector<double> &, class ROL::LinearOperator<double> &, int &, int &) --> double", pybind11::arg("x"), pybind11::arg("A"), pybind11::arg("b"), pybind11::arg("M"), pybind11::arg("iter"), pybind11::arg("flag"));
		cl.def("resetAbsoluteTolerance", (void (ROL::Krylov<double>::*)(const double)) &ROL::Krylov<double>::resetAbsoluteTolerance, "C++: ROL::Krylov<double>::resetAbsoluteTolerance(const double) --> void", pybind11::arg("absTol"));
		cl.def("resetRelativeTolerance", (void (ROL::Krylov<double>::*)(const double)) &ROL::Krylov<double>::resetRelativeTolerance, "C++: ROL::Krylov<double>::resetRelativeTolerance(const double) --> void", pybind11::arg("relTol"));
		cl.def("resetMaximumIteration", (void (ROL::Krylov<double>::*)(const unsigned int)) &ROL::Krylov<double>::resetMaximumIteration, "C++: ROL::Krylov<double>::resetMaximumIteration(const unsigned int) --> void", pybind11::arg("maxit"));
		cl.def("getAbsoluteTolerance", (double (ROL::Krylov<double>::*)() const) &ROL::Krylov<double>::getAbsoluteTolerance, "C++: ROL::Krylov<double>::getAbsoluteTolerance() const --> double");
		cl.def("getRelativeTolerance", (double (ROL::Krylov<double>::*)() const) &ROL::Krylov<double>::getRelativeTolerance, "C++: ROL::Krylov<double>::getRelativeTolerance() const --> double");
		cl.def("getMaximumIteration", (unsigned int (ROL::Krylov<double>::*)() const) &ROL::Krylov<double>::getMaximumIteration, "C++: ROL::Krylov<double>::getMaximumIteration() const --> unsigned int");
		cl.def("assign", (class ROL::Krylov<double> & (ROL::Krylov<double>::*)(const class ROL::Krylov<double> &)) &ROL::Krylov<double>::operator=, "C++: ROL::Krylov<double>::operator=(const class ROL::Krylov<double> &) --> class ROL::Krylov<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::ConjugateGradients file:ROL_ConjugateGradients.hpp line:57
		pybind11::class_<ROL::ConjugateGradients<double>, Teuchos::RCP<ROL::ConjugateGradients<double>>, PyCallBack_ROL_ConjugateGradients_double_t, ROL::Krylov<double>> cl(M("ROL"), "ConjugateGradients_double_t", "");
		cl.def( pybind11::init( [](){ return new ROL::ConjugateGradients<double>(); }, [](){ return new PyCallBack_ROL_ConjugateGradients_double_t(); } ), "doc");
		cl.def( pybind11::init( [](double const & a0){ return new ROL::ConjugateGradients<double>(a0); }, [](double const & a0){ return new PyCallBack_ROL_ConjugateGradients_double_t(a0); } ), "doc");
		cl.def( pybind11::init( [](double const & a0, double const & a1){ return new ROL::ConjugateGradients<double>(a0, a1); }, [](double const & a0, double const & a1){ return new PyCallBack_ROL_ConjugateGradients_double_t(a0, a1); } ), "doc");
		cl.def( pybind11::init( [](double const & a0, double const & a1, unsigned int const & a2){ return new ROL::ConjugateGradients<double>(a0, a1, a2); }, [](double const & a0, double const & a1, unsigned int const & a2){ return new PyCallBack_ROL_ConjugateGradients_double_t(a0, a1, a2); } ), "doc");
		cl.def( pybind11::init<double, double, unsigned int, bool>(), pybind11::arg("absTol"), pybind11::arg("relTol"), pybind11::arg("maxit"), pybind11::arg("useInexact") );

		cl.def( pybind11::init( [](PyCallBack_ROL_ConjugateGradients_double_t const &o){ return new PyCallBack_ROL_ConjugateGradients_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::ConjugateGradients<double> const &o){ return new ROL::ConjugateGradients<double>(o); } ) );
		cl.def("run", (double (ROL::ConjugateGradients<double>::*)(class ROL::Vector<double> &, class ROL::LinearOperator<double> &, const class ROL::Vector<double> &, class ROL::LinearOperator<double> &, int &, int &)) &ROL::ConjugateGradients<double>::run, "C++: ROL::ConjugateGradients<double>::run(class ROL::Vector<double> &, class ROL::LinearOperator<double> &, const class ROL::Vector<double> &, class ROL::LinearOperator<double> &, int &, int &) --> double", pybind11::arg("x"), pybind11::arg("A"), pybind11::arg("b"), pybind11::arg("M"), pybind11::arg("iter"), pybind11::arg("flag"));
		cl.def("assign", (class ROL::ConjugateGradients<double> & (ROL::ConjugateGradients<double>::*)(const class ROL::ConjugateGradients<double> &)) &ROL::ConjugateGradients<double>::operator=, "C++: ROL::ConjugateGradients<double>::operator=(const class ROL::ConjugateGradients<double> &) --> class ROL::ConjugateGradients<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("run", (double (ROL::Krylov<double>::*)(class ROL::Vector<double> &, class ROL::LinearOperator<double> &, const class ROL::Vector<double> &, class ROL::LinearOperator<double> &, int &, int &)) &ROL::Krylov<double>::run, "C++: ROL::Krylov<double>::run(class ROL::Vector<double> &, class ROL::LinearOperator<double> &, const class ROL::Vector<double> &, class ROL::LinearOperator<double> &, int &, int &) --> double", pybind11::arg("x"), pybind11::arg("A"), pybind11::arg("b"), pybind11::arg("M"), pybind11::arg("iter"), pybind11::arg("flag"));
		cl.def("resetAbsoluteTolerance", (void (ROL::Krylov<double>::*)(const double)) &ROL::Krylov<double>::resetAbsoluteTolerance, "C++: ROL::Krylov<double>::resetAbsoluteTolerance(const double) --> void", pybind11::arg("absTol"));
		cl.def("resetRelativeTolerance", (void (ROL::Krylov<double>::*)(const double)) &ROL::Krylov<double>::resetRelativeTolerance, "C++: ROL::Krylov<double>::resetRelativeTolerance(const double) --> void", pybind11::arg("relTol"));
		cl.def("resetMaximumIteration", (void (ROL::Krylov<double>::*)(const unsigned int)) &ROL::Krylov<double>::resetMaximumIteration, "C++: ROL::Krylov<double>::resetMaximumIteration(const unsigned int) --> void", pybind11::arg("maxit"));
		cl.def("getAbsoluteTolerance", (double (ROL::Krylov<double>::*)() const) &ROL::Krylov<double>::getAbsoluteTolerance, "C++: ROL::Krylov<double>::getAbsoluteTolerance() const --> double");
		cl.def("getRelativeTolerance", (double (ROL::Krylov<double>::*)() const) &ROL::Krylov<double>::getRelativeTolerance, "C++: ROL::Krylov<double>::getRelativeTolerance() const --> double");
		cl.def("getMaximumIteration", (unsigned int (ROL::Krylov<double>::*)() const) &ROL::Krylov<double>::getMaximumIteration, "C++: ROL::Krylov<double>::getMaximumIteration() const --> unsigned int");
		cl.def("assign", (class ROL::Krylov<double> & (ROL::Krylov<double>::*)(const class ROL::Krylov<double> &)) &ROL::Krylov<double>::operator=, "C++: ROL::Krylov<double>::operator=(const class ROL::Krylov<double> &) --> class ROL::Krylov<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::ConjugateResiduals file:ROL_ConjugateResiduals.hpp line:57
		pybind11::class_<ROL::ConjugateResiduals<double>, Teuchos::RCP<ROL::ConjugateResiduals<double>>, PyCallBack_ROL_ConjugateResiduals_double_t, ROL::Krylov<double>> cl(M("ROL"), "ConjugateResiduals_double_t", "");
		cl.def( pybind11::init( [](){ return new ROL::ConjugateResiduals<double>(); }, [](){ return new PyCallBack_ROL_ConjugateResiduals_double_t(); } ), "doc");
		cl.def( pybind11::init( [](double const & a0){ return new ROL::ConjugateResiduals<double>(a0); }, [](double const & a0){ return new PyCallBack_ROL_ConjugateResiduals_double_t(a0); } ), "doc");
		cl.def( pybind11::init( [](double const & a0, double const & a1){ return new ROL::ConjugateResiduals<double>(a0, a1); }, [](double const & a0, double const & a1){ return new PyCallBack_ROL_ConjugateResiduals_double_t(a0, a1); } ), "doc");
		cl.def( pybind11::init( [](double const & a0, double const & a1, int const & a2){ return new ROL::ConjugateResiduals<double>(a0, a1, a2); }, [](double const & a0, double const & a1, int const & a2){ return new PyCallBack_ROL_ConjugateResiduals_double_t(a0, a1, a2); } ), "doc");
		cl.def( pybind11::init<double, double, int, bool>(), pybind11::arg("absTol"), pybind11::arg("relTol"), pybind11::arg("maxit"), pybind11::arg("useInexact") );

		cl.def( pybind11::init( [](PyCallBack_ROL_ConjugateResiduals_double_t const &o){ return new PyCallBack_ROL_ConjugateResiduals_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::ConjugateResiduals<double> const &o){ return new ROL::ConjugateResiduals<double>(o); } ) );
		cl.def("run", (double (ROL::ConjugateResiduals<double>::*)(class ROL::Vector<double> &, class ROL::LinearOperator<double> &, const class ROL::Vector<double> &, class ROL::LinearOperator<double> &, int &, int &)) &ROL::ConjugateResiduals<double>::run, "C++: ROL::ConjugateResiduals<double>::run(class ROL::Vector<double> &, class ROL::LinearOperator<double> &, const class ROL::Vector<double> &, class ROL::LinearOperator<double> &, int &, int &) --> double", pybind11::arg("x"), pybind11::arg("A"), pybind11::arg("b"), pybind11::arg("M"), pybind11::arg("iter"), pybind11::arg("flag"));
		cl.def("assign", (class ROL::ConjugateResiduals<double> & (ROL::ConjugateResiduals<double>::*)(const class ROL::ConjugateResiduals<double> &)) &ROL::ConjugateResiduals<double>::operator=, "C++: ROL::ConjugateResiduals<double>::operator=(const class ROL::ConjugateResiduals<double> &) --> class ROL::ConjugateResiduals<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("run", (double (ROL::Krylov<double>::*)(class ROL::Vector<double> &, class ROL::LinearOperator<double> &, const class ROL::Vector<double> &, class ROL::LinearOperator<double> &, int &, int &)) &ROL::Krylov<double>::run, "C++: ROL::Krylov<double>::run(class ROL::Vector<double> &, class ROL::LinearOperator<double> &, const class ROL::Vector<double> &, class ROL::LinearOperator<double> &, int &, int &) --> double", pybind11::arg("x"), pybind11::arg("A"), pybind11::arg("b"), pybind11::arg("M"), pybind11::arg("iter"), pybind11::arg("flag"));
		cl.def("resetAbsoluteTolerance", (void (ROL::Krylov<double>::*)(const double)) &ROL::Krylov<double>::resetAbsoluteTolerance, "C++: ROL::Krylov<double>::resetAbsoluteTolerance(const double) --> void", pybind11::arg("absTol"));
		cl.def("resetRelativeTolerance", (void (ROL::Krylov<double>::*)(const double)) &ROL::Krylov<double>::resetRelativeTolerance, "C++: ROL::Krylov<double>::resetRelativeTolerance(const double) --> void", pybind11::arg("relTol"));
		cl.def("resetMaximumIteration", (void (ROL::Krylov<double>::*)(const unsigned int)) &ROL::Krylov<double>::resetMaximumIteration, "C++: ROL::Krylov<double>::resetMaximumIteration(const unsigned int) --> void", pybind11::arg("maxit"));
		cl.def("getAbsoluteTolerance", (double (ROL::Krylov<double>::*)() const) &ROL::Krylov<double>::getAbsoluteTolerance, "C++: ROL::Krylov<double>::getAbsoluteTolerance() const --> double");
		cl.def("getRelativeTolerance", (double (ROL::Krylov<double>::*)() const) &ROL::Krylov<double>::getRelativeTolerance, "C++: ROL::Krylov<double>::getRelativeTolerance() const --> double");
		cl.def("getMaximumIteration", (unsigned int (ROL::Krylov<double>::*)() const) &ROL::Krylov<double>::getMaximumIteration, "C++: ROL::Krylov<double>::getMaximumIteration() const --> unsigned int");
		cl.def("assign", (class ROL::Krylov<double> & (ROL::Krylov<double>::*)(const class ROL::Krylov<double> &)) &ROL::Krylov<double>::operator=, "C++: ROL::Krylov<double>::operator=(const class ROL::Krylov<double> &) --> class ROL::Krylov<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::GMRES file:ROL_GMRES.hpp line:60
		pybind11::class_<ROL::GMRES<double>, Teuchos::RCP<ROL::GMRES<double>>, PyCallBack_ROL_GMRES_double_t, ROL::Krylov<double>> cl(M("ROL"), "GMRES_double_t", "");
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("parlist") );

		cl.def( pybind11::init( [](PyCallBack_ROL_GMRES_double_t const &o){ return new PyCallBack_ROL_GMRES_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::GMRES<double> const &o){ return new ROL::GMRES<double>(o); } ) );
		cl.def("run", (double (ROL::GMRES<double>::*)(class ROL::Vector<double> &, class ROL::LinearOperator<double> &, const class ROL::Vector<double> &, class ROL::LinearOperator<double> &, int &, int &)) &ROL::GMRES<double>::run, "C++: ROL::GMRES<double>::run(class ROL::Vector<double> &, class ROL::LinearOperator<double> &, const class ROL::Vector<double> &, class ROL::LinearOperator<double> &, int &, int &) --> double", pybind11::arg("x"), pybind11::arg("A"), pybind11::arg("b"), pybind11::arg("M"), pybind11::arg("iter"), pybind11::arg("flag"));
		cl.def("disableOutput", (void (ROL::GMRES<double>::*)()) &ROL::GMRES<double>::disableOutput, "C++: ROL::GMRES<double>::disableOutput() --> void");
		cl.def("assign", (class ROL::GMRES<double> & (ROL::GMRES<double>::*)(const class ROL::GMRES<double> &)) &ROL::GMRES<double>::operator=, "C++: ROL::GMRES<double>::operator=(const class ROL::GMRES<double> &) --> class ROL::GMRES<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("run", (double (ROL::Krylov<double>::*)(class ROL::Vector<double> &, class ROL::LinearOperator<double> &, const class ROL::Vector<double> &, class ROL::LinearOperator<double> &, int &, int &)) &ROL::Krylov<double>::run, "C++: ROL::Krylov<double>::run(class ROL::Vector<double> &, class ROL::LinearOperator<double> &, const class ROL::Vector<double> &, class ROL::LinearOperator<double> &, int &, int &) --> double", pybind11::arg("x"), pybind11::arg("A"), pybind11::arg("b"), pybind11::arg("M"), pybind11::arg("iter"), pybind11::arg("flag"));
		cl.def("resetAbsoluteTolerance", (void (ROL::Krylov<double>::*)(const double)) &ROL::Krylov<double>::resetAbsoluteTolerance, "C++: ROL::Krylov<double>::resetAbsoluteTolerance(const double) --> void", pybind11::arg("absTol"));
		cl.def("resetRelativeTolerance", (void (ROL::Krylov<double>::*)(const double)) &ROL::Krylov<double>::resetRelativeTolerance, "C++: ROL::Krylov<double>::resetRelativeTolerance(const double) --> void", pybind11::arg("relTol"));
		cl.def("resetMaximumIteration", (void (ROL::Krylov<double>::*)(const unsigned int)) &ROL::Krylov<double>::resetMaximumIteration, "C++: ROL::Krylov<double>::resetMaximumIteration(const unsigned int) --> void", pybind11::arg("maxit"));
		cl.def("getAbsoluteTolerance", (double (ROL::Krylov<double>::*)() const) &ROL::Krylov<double>::getAbsoluteTolerance, "C++: ROL::Krylov<double>::getAbsoluteTolerance() const --> double");
		cl.def("getRelativeTolerance", (double (ROL::Krylov<double>::*)() const) &ROL::Krylov<double>::getRelativeTolerance, "C++: ROL::Krylov<double>::getRelativeTolerance() const --> double");
		cl.def("getMaximumIteration", (unsigned int (ROL::Krylov<double>::*)() const) &ROL::Krylov<double>::getMaximumIteration, "C++: ROL::Krylov<double>::getMaximumIteration() const --> unsigned int");
		cl.def("assign", (class ROL::Krylov<double> & (ROL::Krylov<double>::*)(const class ROL::Krylov<double> &)) &ROL::Krylov<double>::operator=, "C++: ROL::Krylov<double>::operator=(const class ROL::Krylov<double> &) --> class ROL::Krylov<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
}
