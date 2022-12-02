#include <ROL_BoundConstraint.hpp>
#include <ROL_ColemanLiModel.hpp>
#include <ROL_Elementwise_Function.hpp>
#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_Objective.hpp>
#include <ROL_Secant.hpp>
#include <ROL_TrustRegionModel.hpp>
#include <ROL_TrustRegionTypes.hpp>
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

// ROL::TrustRegionModel file:ROL_TrustRegionModel.hpp line:67
struct PyCallBack_ROL_TrustRegionModel_double_t : public ROL::TrustRegionModel<double> {
	using ROL::TrustRegionModel<double>::TrustRegionModel;

	void update(class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class Teuchos::RCP<class ROL::Secant<double> > & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionModel<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TrustRegionModel::update(a0, a1, a2, a3, a4);
	}
	double value(const class ROL::Vector<double> & a0, double & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionModel<double> *>(this), "value");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return TrustRegionModel::value(a0, a1);
	}
	void gradient(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double & a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionModel<double> *>(this), "gradient");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TrustRegionModel::gradient(a0, a1, a2);
	}
	void hessVec(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionModel<double> *>(this), "hessVec");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TrustRegionModel::hessVec(a0, a1, a2, a3);
	}
	void invHessVec(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionModel<double> *>(this), "invHessVec");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TrustRegionModel::invHessVec(a0, a1, a2, a3);
	}
	void precond(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionModel<double> *>(this), "precond");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TrustRegionModel::precond(a0, a1, a2, a3);
	}
	const class Teuchos::RCP<const class ROL::Vector<double> > getGradient() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionModel<double> *>(this), "getGradient");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const class Teuchos::RCP<const class ROL::Vector<double> >>::value) {
				static pybind11::detail::override_caster_t<const class Teuchos::RCP<const class ROL::Vector<double> >> caster;
				return pybind11::detail::cast_ref<const class Teuchos::RCP<const class ROL::Vector<double> >>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const class Teuchos::RCP<const class ROL::Vector<double> >>(std::move(o));
		}
		return TrustRegionModel::getGradient();
	}
	const class Teuchos::RCP<const class ROL::Vector<double> > getIterate() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionModel<double> *>(this), "getIterate");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const class Teuchos::RCP<const class ROL::Vector<double> >>::value) {
				static pybind11::detail::override_caster_t<const class Teuchos::RCP<const class ROL::Vector<double> >> caster;
				return pybind11::detail::cast_ref<const class Teuchos::RCP<const class ROL::Vector<double> >>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const class Teuchos::RCP<const class ROL::Vector<double> >>(std::move(o));
		}
		return TrustRegionModel::getIterate();
	}
	const class Teuchos::RCP<class ROL::Objective<double> > getObjective() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionModel<double> *>(this), "getObjective");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const class Teuchos::RCP<class ROL::Objective<double> >>::value) {
				static pybind11::detail::override_caster_t<const class Teuchos::RCP<class ROL::Objective<double> >> caster;
				return pybind11::detail::cast_ref<const class Teuchos::RCP<class ROL::Objective<double> >>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const class Teuchos::RCP<class ROL::Objective<double> >>(std::move(o));
		}
		return TrustRegionModel::getObjective();
	}
	const class Teuchos::RCP<class ROL::BoundConstraint<double> > getBoundConstraint() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionModel<double> *>(this), "getBoundConstraint");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const class Teuchos::RCP<class ROL::BoundConstraint<double> >>::value) {
				static pybind11::detail::override_caster_t<const class Teuchos::RCP<class ROL::BoundConstraint<double> >> caster;
				return pybind11::detail::cast_ref<const class Teuchos::RCP<class ROL::BoundConstraint<double> >>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const class Teuchos::RCP<class ROL::BoundConstraint<double> >>(std::move(o));
		}
		return TrustRegionModel::getBoundConstraint();
	}
	void dualTransform(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionModel<double> *>(this), "dualTransform");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TrustRegionModel::dualTransform(a0, a1);
	}
	void primalTransform(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionModel<double> *>(this), "primalTransform");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TrustRegionModel::primalTransform(a0, a1);
	}
	void updatePredictedReduction(double & a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionModel<double> *>(this), "updatePredictedReduction");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TrustRegionModel::updatePredictedReduction(a0, a1);
	}
	void updateActualReduction(double & a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionModel<double> *>(this), "updateActualReduction");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TrustRegionModel::updateActualReduction(a0, a1);
	}
	void update(const class ROL::Vector<double> & a0, enum ROL::UpdateType a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionModel<double> *>(this), "update");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionModel<double> *>(this), "update");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionModel<double> *>(this), "dirDeriv");
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
};

// ROL::ColemanLiModel file:ROL_ColemanLiModel.hpp line:60
struct PyCallBack_ROL_ColemanLiModel_double_t : public ROL::ColemanLiModel<double> {
	using ROL::ColemanLiModel<double>::ColemanLiModel;

	void update(class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class Teuchos::RCP<class ROL::Secant<double> > & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::ColemanLiModel<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ColemanLiModel::update(a0, a1, a2, a3, a4);
	}
	double value(const class ROL::Vector<double> & a0, double & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::ColemanLiModel<double> *>(this), "value");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return ColemanLiModel::value(a0, a1);
	}
	void gradient(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double & a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::ColemanLiModel<double> *>(this), "gradient");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ColemanLiModel::gradient(a0, a1, a2);
	}
	void hessVec(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::ColemanLiModel<double> *>(this), "hessVec");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ColemanLiModel::hessVec(a0, a1, a2, a3);
	}
	void dualTransform(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::ColemanLiModel<double> *>(this), "dualTransform");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ColemanLiModel::dualTransform(a0, a1);
	}
	void primalTransform(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::ColemanLiModel<double> *>(this), "primalTransform");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ColemanLiModel::primalTransform(a0, a1);
	}
	void updatePredictedReduction(double & a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::ColemanLiModel<double> *>(this), "updatePredictedReduction");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ColemanLiModel::updatePredictedReduction(a0, a1);
	}
	void updateActualReduction(double & a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::ColemanLiModel<double> *>(this), "updateActualReduction");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ColemanLiModel::updateActualReduction(a0, a1);
	}
	void invHessVec(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::ColemanLiModel<double> *>(this), "invHessVec");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TrustRegionModel::invHessVec(a0, a1, a2, a3);
	}
	void precond(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::ColemanLiModel<double> *>(this), "precond");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TrustRegionModel::precond(a0, a1, a2, a3);
	}
	const class Teuchos::RCP<const class ROL::Vector<double> > getGradient() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::ColemanLiModel<double> *>(this), "getGradient");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const class Teuchos::RCP<const class ROL::Vector<double> >>::value) {
				static pybind11::detail::override_caster_t<const class Teuchos::RCP<const class ROL::Vector<double> >> caster;
				return pybind11::detail::cast_ref<const class Teuchos::RCP<const class ROL::Vector<double> >>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const class Teuchos::RCP<const class ROL::Vector<double> >>(std::move(o));
		}
		return TrustRegionModel::getGradient();
	}
	const class Teuchos::RCP<const class ROL::Vector<double> > getIterate() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::ColemanLiModel<double> *>(this), "getIterate");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const class Teuchos::RCP<const class ROL::Vector<double> >>::value) {
				static pybind11::detail::override_caster_t<const class Teuchos::RCP<const class ROL::Vector<double> >> caster;
				return pybind11::detail::cast_ref<const class Teuchos::RCP<const class ROL::Vector<double> >>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const class Teuchos::RCP<const class ROL::Vector<double> >>(std::move(o));
		}
		return TrustRegionModel::getIterate();
	}
	const class Teuchos::RCP<class ROL::Objective<double> > getObjective() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::ColemanLiModel<double> *>(this), "getObjective");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const class Teuchos::RCP<class ROL::Objective<double> >>::value) {
				static pybind11::detail::override_caster_t<const class Teuchos::RCP<class ROL::Objective<double> >> caster;
				return pybind11::detail::cast_ref<const class Teuchos::RCP<class ROL::Objective<double> >>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const class Teuchos::RCP<class ROL::Objective<double> >>(std::move(o));
		}
		return TrustRegionModel::getObjective();
	}
	const class Teuchos::RCP<class ROL::BoundConstraint<double> > getBoundConstraint() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::ColemanLiModel<double> *>(this), "getBoundConstraint");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const class Teuchos::RCP<class ROL::BoundConstraint<double> >>::value) {
				static pybind11::detail::override_caster_t<const class Teuchos::RCP<class ROL::BoundConstraint<double> >> caster;
				return pybind11::detail::cast_ref<const class Teuchos::RCP<class ROL::BoundConstraint<double> >>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const class Teuchos::RCP<class ROL::BoundConstraint<double> >>(std::move(o));
		}
		return TrustRegionModel::getBoundConstraint();
	}
	void update(const class ROL::Vector<double> & a0, enum ROL::UpdateType a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::ColemanLiModel<double> *>(this), "update");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::ColemanLiModel<double> *>(this), "update");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::ColemanLiModel<double> *>(this), "dirDeriv");
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
};

void bind_ROL_TrustRegionTypes(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	// ROL::ETrustRegion file:ROL_TrustRegionTypes.hpp line:64
	pybind11::enum_<ROL::ETrustRegion>(M("ROL"), "ETrustRegion", pybind11::arithmetic(), "Enumeration of trust-region solver types.\n\n      \n    CAUCHYPOINT     describe\n      \n\n    TRUNCATEDCG     describe\n      \n\n    DOGLEG          describe\n      \n\n    DOUBLEDOGLEG    describe")
		.value("TRUSTREGION_CAUCHYPOINT", ROL::TRUSTREGION_CAUCHYPOINT)
		.value("TRUSTREGION_TRUNCATEDCG", ROL::TRUSTREGION_TRUNCATEDCG)
		.value("TRUSTREGION_DOGLEG", ROL::TRUSTREGION_DOGLEG)
		.value("TRUSTREGION_DOUBLEDOGLEG", ROL::TRUSTREGION_DOUBLEDOGLEG)
		.value("TRUSTREGION_LINMORE", ROL::TRUSTREGION_LINMORE)
		.value("TRUSTREGION_LAST", ROL::TRUSTREGION_LAST)
		.export_values();

;

	// ROL::ETrustRegionToString(enum ROL::ETrustRegion) file:ROL_TrustRegionTypes.hpp line:73
	M("ROL").def("ETrustRegionToString", (std::string (*)(enum ROL::ETrustRegion)) &ROL::ETrustRegionToString, "C++: ROL::ETrustRegionToString(enum ROL::ETrustRegion) --> std::string", pybind11::arg("tr"));

	// ROL::isValidTrustRegion(enum ROL::ETrustRegion) file:ROL_TrustRegionTypes.hpp line:92
	M("ROL").def("isValidTrustRegion", (int (*)(enum ROL::ETrustRegion)) &ROL::isValidTrustRegion, "Verifies validity of a TrustRegion enum.\n\n      \n  [in]  - enum of the TrustRegion\n      \n\n 1 if the argument is a valid TrustRegion; 0 otherwise.\n\nC++: ROL::isValidTrustRegion(enum ROL::ETrustRegion) --> int", pybind11::arg("ls"));

	// ROL::StringToETrustRegion(std::string) file:ROL_TrustRegionTypes.hpp line:121
	M("ROL").def("StringToETrustRegion", (enum ROL::ETrustRegion (*)(std::string)) &ROL::StringToETrustRegion, "C++: ROL::StringToETrustRegion(std::string) --> enum ROL::ETrustRegion", pybind11::arg("s"));

	// ROL::ETrustRegionModel file:ROL_TrustRegionTypes.hpp line:137
	pybind11::enum_<ROL::ETrustRegionModel>(M("ROL"), "ETrustRegionModel", pybind11::arithmetic(), "Enumeration of trust-region model types.\n\n      \n    COLEMANLI   describe\n      \n\n    KELLEYSACHS describe")
		.value("TRUSTREGION_MODEL_COLEMANLI", ROL::TRUSTREGION_MODEL_COLEMANLI)
		.value("TRUSTREGION_MODEL_KELLEYSACHS", ROL::TRUSTREGION_MODEL_KELLEYSACHS)
		.value("TRUSTREGION_MODEL_LINMORE", ROL::TRUSTREGION_MODEL_LINMORE)
		.value("TRUSTREGION_MODEL_LAST", ROL::TRUSTREGION_MODEL_LAST)
		.export_values();

;

	// ROL::ETrustRegionModelToString(enum ROL::ETrustRegionModel) file:ROL_TrustRegionTypes.hpp line:144
	M("ROL").def("ETrustRegionModelToString", (std::string (*)(enum ROL::ETrustRegionModel)) &ROL::ETrustRegionModelToString, "C++: ROL::ETrustRegionModelToString(enum ROL::ETrustRegionModel) --> std::string", pybind11::arg("tr"));

	// ROL::isValidTrustRegionModel(enum ROL::ETrustRegionModel) file:ROL_TrustRegionTypes.hpp line:161
	M("ROL").def("isValidTrustRegionModel", (int (*)(enum ROL::ETrustRegionModel)) &ROL::isValidTrustRegionModel, "Verifies validity of a TrustRegionModel enum.\n\n      \n  [in]  - enum of the TrustRegionModel\n      \n\n 1 if the argument is a valid TrustRegionModel; 0 otherwise.\n\nC++: ROL::isValidTrustRegionModel(enum ROL::ETrustRegionModel) --> int", pybind11::arg("ls"));

	// ROL::StringToETrustRegionModel(std::string) file:ROL_TrustRegionTypes.hpp line:188
	M("ROL").def("StringToETrustRegionModel", (enum ROL::ETrustRegionModel (*)(std::string)) &ROL::StringToETrustRegionModel, "C++: ROL::StringToETrustRegionModel(std::string) --> enum ROL::ETrustRegionModel", pybind11::arg("s"));

	// ROL::isValidTrustRegionSubproblem(enum ROL::ETrustRegion, enum ROL::ETrustRegionModel, bool) file:ROL_TrustRegionTypes.hpp line:198
	M("ROL").def("isValidTrustRegionSubproblem", (bool (*)(enum ROL::ETrustRegion, enum ROL::ETrustRegionModel, bool)) &ROL::isValidTrustRegionSubproblem, "C++: ROL::isValidTrustRegionSubproblem(enum ROL::ETrustRegion, enum ROL::ETrustRegionModel, bool) --> bool", pybind11::arg("etr"), pybind11::arg("etrm"), pybind11::arg("isBnd"));

	// ROL::ETrustRegionFlag file:ROL_TrustRegionTypes.hpp line:232
	pybind11::enum_<ROL::ETrustRegionFlag>(M("ROL"), "ETrustRegionFlag", pybind11::arithmetic(), "Enumation of flags used by trust-region solvers.\n\n      \n TRUSTREGION_FLAG_SUCCESS        Actual and predicted reductions are positive \n      \n\n TRUSTREGION_FLAG_POSPREDNEG     Reduction is positive, predicted negative (impossible)\n      \n\n TRUSTREGION_FLAG_NPOSPREDPOS    Reduction is nonpositive, predicted positive\n      \n\n TRUSTREGION_FLAG_NPOSPREDNEG    Reduction is nonpositive, predicted negative (impossible)\n      \n\n TRUSTREGION_FLAG_QMINSUFDEC     Insufficient decrease of the quadratic model (bound constraint only)\n      \n\n TRUSTREGION_FLAG_NAN            Actual and/or predicted reduction is NaN\n\n  ")
		.value("TRUSTREGION_FLAG_SUCCESS", ROL::TRUSTREGION_FLAG_SUCCESS)
		.value("TRUSTREGION_FLAG_POSPREDNEG", ROL::TRUSTREGION_FLAG_POSPREDNEG)
		.value("TRUSTREGION_FLAG_NPOSPREDPOS", ROL::TRUSTREGION_FLAG_NPOSPREDPOS)
		.value("TRUSTREGION_FLAG_NPOSPREDNEG", ROL::TRUSTREGION_FLAG_NPOSPREDNEG)
		.value("TRUSTREGION_FLAG_QMINSUFDEC", ROL::TRUSTREGION_FLAG_QMINSUFDEC)
		.value("TRUSTREGION_FLAG_NAN", ROL::TRUSTREGION_FLAG_NAN)
		.value("TRUSTREGION_FLAG_UNDEFINED", ROL::TRUSTREGION_FLAG_UNDEFINED)
		.export_values();

;

	// ROL::ETrustRegionFlagToString(enum ROL::ETrustRegionFlag) file:ROL_TrustRegionTypes.hpp line:243
	M("ROL").def("ETrustRegionFlagToString", (std::string (*)(enum ROL::ETrustRegionFlag)) &ROL::ETrustRegionFlagToString, "C++: ROL::ETrustRegionFlagToString(enum ROL::ETrustRegionFlag) --> std::string", pybind11::arg("trf"));

	{ // ROL::TrustRegionModel file:ROL_TrustRegionModel.hpp line:67
		pybind11::class_<ROL::TrustRegionModel<double>, Teuchos::RCP<ROL::TrustRegionModel<double>>, PyCallBack_ROL_TrustRegionModel_double_t, ROL::Objective<double>> cl(M("ROL"), "TrustRegionModel_double_t", "");
		cl.def( pybind11::init( [](class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3){ return new ROL::TrustRegionModel<double>(a0, a1, a2, a3); }, [](class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3){ return new PyCallBack_ROL_TrustRegionModel_double_t(a0, a1, a2, a3); } ), "doc");
		cl.def( pybind11::init( [](class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class Teuchos::RCP<class ROL::Secant<double> > & a4){ return new ROL::TrustRegionModel<double>(a0, a1, a2, a3, a4); }, [](class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class Teuchos::RCP<class ROL::Secant<double> > & a4){ return new PyCallBack_ROL_TrustRegionModel_double_t(a0, a1, a2, a3, a4); } ), "doc");
		cl.def( pybind11::init( [](class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class Teuchos::RCP<class ROL::Secant<double> > & a4, const bool & a5){ return new ROL::TrustRegionModel<double>(a0, a1, a2, a3, a4, a5); }, [](class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class Teuchos::RCP<class ROL::Secant<double> > & a4, const bool & a5){ return new PyCallBack_ROL_TrustRegionModel_double_t(a0, a1, a2, a3, a4, a5); } ), "doc");
		cl.def( pybind11::init<class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class Teuchos::RCP<class ROL::Secant<double> > &, const bool, const bool>(), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("secant"), pybind11::arg("useSecantPrecond"), pybind11::arg("useSecantHessVec") );

		cl.def( pybind11::init( [](PyCallBack_ROL_TrustRegionModel_double_t const &o){ return new PyCallBack_ROL_TrustRegionModel_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::TrustRegionModel<double> const &o){ return new ROL::TrustRegionModel<double>(o); } ) );
		cl.def("update", [](ROL::TrustRegionModel<double> &o, const class ROL::Vector<double> & a0) -> void { return o.update(a0); }, "", pybind11::arg("x"));
		cl.def("update", [](ROL::TrustRegionModel<double> &o, const class ROL::Vector<double> & a0, bool const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("flag"));
		cl.def("update", [](ROL::TrustRegionModel<double> &o, const class ROL::Vector<double> & a0, bool const & a1, int const & a2) -> void { return o.update(a0, a1, a2); }, "", pybind11::arg("x"), pybind11::arg("flag"), pybind11::arg("iter"));
		cl.def("update", [](ROL::TrustRegionModel<double> &o, const class ROL::Vector<double> & a0, enum ROL::UpdateType const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("type"));
		cl.def("update", [](ROL::TrustRegionModel<double> &o, const class ROL::Vector<double> & a0, enum ROL::UpdateType const & a1, int const & a2) -> void { return o.update(a0, a1, a2); }, "", pybind11::arg("x"), pybind11::arg("type"), pybind11::arg("iter"));
		cl.def("update", [](ROL::TrustRegionModel<double> &o, class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3) -> void { return o.update(a0, a1, a2, a3); }, "", pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("x"), pybind11::arg("g"));
		cl.def("update", (void (ROL::TrustRegionModel<double>::*)(class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class Teuchos::RCP<class ROL::Secant<double> > &)) &ROL::TrustRegionModel<double>::update, "C++: ROL::TrustRegionModel<double>::update(class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class Teuchos::RCP<class ROL::Secant<double> > &) --> void", pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("secant"));
		cl.def("value", (double (ROL::TrustRegionModel<double>::*)(const class ROL::Vector<double> &, double &)) &ROL::TrustRegionModel<double>::value, "C++: ROL::TrustRegionModel<double>::value(const class ROL::Vector<double> &, double &) --> double", pybind11::arg("s"), pybind11::arg("tol"));
		cl.def("gradient", (void (ROL::TrustRegionModel<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::TrustRegionModel<double>::gradient, "C++: ROL::TrustRegionModel<double>::gradient(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("g"), pybind11::arg("s"), pybind11::arg("tol"));
		cl.def("hessVec", (void (ROL::TrustRegionModel<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::TrustRegionModel<double>::hessVec, "C++: ROL::TrustRegionModel<double>::hessVec(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("hv"), pybind11::arg("v"), pybind11::arg("s"), pybind11::arg("tol"));
		cl.def("invHessVec", (void (ROL::TrustRegionModel<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::TrustRegionModel<double>::invHessVec, "C++: ROL::TrustRegionModel<double>::invHessVec(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("hv"), pybind11::arg("v"), pybind11::arg("s"), pybind11::arg("tol"));
		cl.def("precond", (void (ROL::TrustRegionModel<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::TrustRegionModel<double>::precond, "C++: ROL::TrustRegionModel<double>::precond(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("Pv"), pybind11::arg("v"), pybind11::arg("s"), pybind11::arg("tol"));
		cl.def("getGradient", (const class Teuchos::RCP<const class ROL::Vector<double> > (ROL::TrustRegionModel<double>::*)() const) &ROL::TrustRegionModel<double>::getGradient, "C++: ROL::TrustRegionModel<double>::getGradient() const --> const class Teuchos::RCP<const class ROL::Vector<double> >");
		cl.def("getIterate", (const class Teuchos::RCP<const class ROL::Vector<double> > (ROL::TrustRegionModel<double>::*)() const) &ROL::TrustRegionModel<double>::getIterate, "C++: ROL::TrustRegionModel<double>::getIterate() const --> const class Teuchos::RCP<const class ROL::Vector<double> >");
		cl.def("getObjective", (const class Teuchos::RCP<class ROL::Objective<double> > (ROL::TrustRegionModel<double>::*)() const) &ROL::TrustRegionModel<double>::getObjective, "C++: ROL::TrustRegionModel<double>::getObjective() const --> const class Teuchos::RCP<class ROL::Objective<double> >");
		cl.def("getBoundConstraint", (const class Teuchos::RCP<class ROL::BoundConstraint<double> > (ROL::TrustRegionModel<double>::*)() const) &ROL::TrustRegionModel<double>::getBoundConstraint, "C++: ROL::TrustRegionModel<double>::getBoundConstraint() const --> const class Teuchos::RCP<class ROL::BoundConstraint<double> >");
		cl.def("dualTransform", (void (ROL::TrustRegionModel<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::TrustRegionModel<double>::dualTransform, "C++: ROL::TrustRegionModel<double>::dualTransform(class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("tv"), pybind11::arg("v"));
		cl.def("primalTransform", (void (ROL::TrustRegionModel<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::TrustRegionModel<double>::primalTransform, "C++: ROL::TrustRegionModel<double>::primalTransform(class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("tv"), pybind11::arg("v"));
		cl.def("updatePredictedReduction", (void (ROL::TrustRegionModel<double>::*)(double &, const class ROL::Vector<double> &)) &ROL::TrustRegionModel<double>::updatePredictedReduction, "C++: ROL::TrustRegionModel<double>::updatePredictedReduction(double &, const class ROL::Vector<double> &) --> void", pybind11::arg("pred"), pybind11::arg("s"));
		cl.def("updateActualReduction", (void (ROL::TrustRegionModel<double>::*)(double &, const class ROL::Vector<double> &)) &ROL::TrustRegionModel<double>::updateActualReduction, "C++: ROL::TrustRegionModel<double>::updateActualReduction(double &, const class ROL::Vector<double> &) --> void", pybind11::arg("ared"), pybind11::arg("s"));
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
	{ // ROL::ColemanLiModel file:ROL_ColemanLiModel.hpp line:60
		pybind11::class_<ROL::ColemanLiModel<double>, Teuchos::RCP<ROL::ColemanLiModel<double>>, PyCallBack_ROL_ColemanLiModel_double_t, ROL::TrustRegionModel<double>> cl(M("ROL"), "ColemanLiModel_double_t", "");
		cl.def( pybind11::init( [](class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3){ return new ROL::ColemanLiModel<double>(a0, a1, a2, a3); }, [](class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3){ return new PyCallBack_ROL_ColemanLiModel_double_t(a0, a1, a2, a3); } ), "doc");
		cl.def( pybind11::init( [](class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const double & a4){ return new ROL::ColemanLiModel<double>(a0, a1, a2, a3, a4); }, [](class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const double & a4){ return new PyCallBack_ROL_ColemanLiModel_double_t(a0, a1, a2, a3, a4); } ), "doc");
		cl.def( pybind11::init( [](class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const double & a4, const double & a5){ return new ROL::ColemanLiModel<double>(a0, a1, a2, a3, a4, a5); }, [](class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const double & a4, const double & a5){ return new PyCallBack_ROL_ColemanLiModel_double_t(a0, a1, a2, a3, a4, a5); } ), "doc");
		cl.def( pybind11::init( [](class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const double & a4, const double & a5, const bool & a6){ return new ROL::ColemanLiModel<double>(a0, a1, a2, a3, a4, a5, a6); }, [](class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const double & a4, const double & a5, const bool & a6){ return new PyCallBack_ROL_ColemanLiModel_double_t(a0, a1, a2, a3, a4, a5, a6); } ), "doc");
		cl.def( pybind11::init( [](class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const double & a4, const double & a5, const bool & a6, const class Teuchos::RCP<class ROL::Secant<double> > & a7){ return new ROL::ColemanLiModel<double>(a0, a1, a2, a3, a4, a5, a6, a7); }, [](class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const double & a4, const double & a5, const bool & a6, const class Teuchos::RCP<class ROL::Secant<double> > & a7){ return new PyCallBack_ROL_ColemanLiModel_double_t(a0, a1, a2, a3, a4, a5, a6, a7); } ), "doc");
		cl.def( pybind11::init( [](class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const double & a4, const double & a5, const bool & a6, const class Teuchos::RCP<class ROL::Secant<double> > & a7, const bool & a8){ return new ROL::ColemanLiModel<double>(a0, a1, a2, a3, a4, a5, a6, a7, a8); }, [](class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const double & a4, const double & a5, const bool & a6, const class Teuchos::RCP<class ROL::Secant<double> > & a7, const bool & a8){ return new PyCallBack_ROL_ColemanLiModel_double_t(a0, a1, a2, a3, a4, a5, a6, a7, a8); } ), "doc");
		cl.def( pybind11::init<class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const double, const double, const bool, const class Teuchos::RCP<class ROL::Secant<double> > &, const bool, const bool>(), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("stepBackMax"), pybind11::arg("stepBackScale"), pybind11::arg("singleReflect"), pybind11::arg("secant"), pybind11::arg("useSecantPrecond"), pybind11::arg("useSecantHessVec") );

		cl.def( pybind11::init( [](PyCallBack_ROL_ColemanLiModel_double_t const &o){ return new PyCallBack_ROL_ColemanLiModel_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::ColemanLiModel<double> const &o){ return new ROL::ColemanLiModel<double>(o); } ) );
		cl.def("update", [](ROL::ColemanLiModel<double> &o, const class ROL::Vector<double> & a0) -> void { return o.update(a0); }, "", pybind11::arg("x"));
		cl.def("update", [](ROL::ColemanLiModel<double> &o, const class ROL::Vector<double> & a0, bool const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("flag"));
		cl.def("update", [](ROL::ColemanLiModel<double> &o, const class ROL::Vector<double> & a0, bool const & a1, int const & a2) -> void { return o.update(a0, a1, a2); }, "", pybind11::arg("x"), pybind11::arg("flag"), pybind11::arg("iter"));
		cl.def("update", [](ROL::ColemanLiModel<double> &o, const class ROL::Vector<double> & a0, enum ROL::UpdateType const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("type"));
		cl.def("update", [](ROL::ColemanLiModel<double> &o, const class ROL::Vector<double> & a0, enum ROL::UpdateType const & a1, int const & a2) -> void { return o.update(a0, a1, a2); }, "", pybind11::arg("x"), pybind11::arg("type"), pybind11::arg("iter"));
		cl.def("update", [](ROL::ColemanLiModel<double> &o, class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3) -> void { return o.update(a0, a1, a2, a3); }, "", pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("x"), pybind11::arg("g"));
		cl.def("update", (void (ROL::ColemanLiModel<double>::*)(class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class Teuchos::RCP<class ROL::Secant<double> > &)) &ROL::ColemanLiModel<double>::update, "C++: ROL::ColemanLiModel<double>::update(class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class Teuchos::RCP<class ROL::Secant<double> > &) --> void", pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("secant"));
		cl.def("setRadius", (void (ROL::ColemanLiModel<double>::*)(const double)) &ROL::ColemanLiModel<double>::setRadius, "C++: ROL::ColemanLiModel<double>::setRadius(const double) --> void", pybind11::arg("del"));
		cl.def("value", (double (ROL::ColemanLiModel<double>::*)(const class ROL::Vector<double> &, double &)) &ROL::ColemanLiModel<double>::value, "C++: ROL::ColemanLiModel<double>::value(const class ROL::Vector<double> &, double &) --> double", pybind11::arg("s"), pybind11::arg("tol"));
		cl.def("gradient", (void (ROL::ColemanLiModel<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::ColemanLiModel<double>::gradient, "C++: ROL::ColemanLiModel<double>::gradient(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("g"), pybind11::arg("s"), pybind11::arg("tol"));
		cl.def("hessVec", (void (ROL::ColemanLiModel<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::ColemanLiModel<double>::hessVec, "C++: ROL::ColemanLiModel<double>::hessVec(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("hv"), pybind11::arg("v"), pybind11::arg("s"), pybind11::arg("tol"));
		cl.def("dualTransform", (void (ROL::ColemanLiModel<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::ColemanLiModel<double>::dualTransform, "C++: ROL::ColemanLiModel<double>::dualTransform(class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("tv"), pybind11::arg("v"));
		cl.def("primalTransform", (void (ROL::ColemanLiModel<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::ColemanLiModel<double>::primalTransform, "C++: ROL::ColemanLiModel<double>::primalTransform(class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("tiv"), pybind11::arg("v"));
		cl.def("updatePredictedReduction", (void (ROL::ColemanLiModel<double>::*)(double &, const class ROL::Vector<double> &)) &ROL::ColemanLiModel<double>::updatePredictedReduction, "C++: ROL::ColemanLiModel<double>::updatePredictedReduction(double &, const class ROL::Vector<double> &) --> void", pybind11::arg("pred"), pybind11::arg("s"));
		cl.def("updateActualReduction", (void (ROL::ColemanLiModel<double>::*)(double &, const class ROL::Vector<double> &)) &ROL::ColemanLiModel<double>::updateActualReduction, "C++: ROL::ColemanLiModel<double>::updateActualReduction(double &, const class ROL::Vector<double> &) --> void", pybind11::arg("ared"), pybind11::arg("s"));
		cl.def("update", [](ROL::TrustRegionModel<double> &o, const class ROL::Vector<double> & a0) -> void { return o.update(a0); }, "", pybind11::arg("x"));
		cl.def("update", [](ROL::TrustRegionModel<double> &o, const class ROL::Vector<double> & a0, bool const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("flag"));
		cl.def("update", [](ROL::TrustRegionModel<double> &o, const class ROL::Vector<double> & a0, bool const & a1, int const & a2) -> void { return o.update(a0, a1, a2); }, "", pybind11::arg("x"), pybind11::arg("flag"), pybind11::arg("iter"));
		cl.def("update", [](ROL::TrustRegionModel<double> &o, const class ROL::Vector<double> & a0, enum ROL::UpdateType const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("type"));
		cl.def("update", [](ROL::TrustRegionModel<double> &o, const class ROL::Vector<double> & a0, enum ROL::UpdateType const & a1, int const & a2) -> void { return o.update(a0, a1, a2); }, "", pybind11::arg("x"), pybind11::arg("type"), pybind11::arg("iter"));
		cl.def("update", [](ROL::TrustRegionModel<double> &o, class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3) -> void { return o.update(a0, a1, a2, a3); }, "", pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("x"), pybind11::arg("g"));
		cl.def("update", (void (ROL::TrustRegionModel<double>::*)(class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class Teuchos::RCP<class ROL::Secant<double> > &)) &ROL::TrustRegionModel<double>::update, "C++: ROL::TrustRegionModel<double>::update(class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class Teuchos::RCP<class ROL::Secant<double> > &) --> void", pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("secant"));
		cl.def("value", (double (ROL::TrustRegionModel<double>::*)(const class ROL::Vector<double> &, double &)) &ROL::TrustRegionModel<double>::value, "C++: ROL::TrustRegionModel<double>::value(const class ROL::Vector<double> &, double &) --> double", pybind11::arg("s"), pybind11::arg("tol"));
		cl.def("gradient", (void (ROL::TrustRegionModel<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::TrustRegionModel<double>::gradient, "C++: ROL::TrustRegionModel<double>::gradient(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("g"), pybind11::arg("s"), pybind11::arg("tol"));
		cl.def("hessVec", (void (ROL::TrustRegionModel<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::TrustRegionModel<double>::hessVec, "C++: ROL::TrustRegionModel<double>::hessVec(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("hv"), pybind11::arg("v"), pybind11::arg("s"), pybind11::arg("tol"));
		cl.def("invHessVec", (void (ROL::TrustRegionModel<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::TrustRegionModel<double>::invHessVec, "C++: ROL::TrustRegionModel<double>::invHessVec(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("hv"), pybind11::arg("v"), pybind11::arg("s"), pybind11::arg("tol"));
		cl.def("precond", (void (ROL::TrustRegionModel<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::TrustRegionModel<double>::precond, "C++: ROL::TrustRegionModel<double>::precond(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("Pv"), pybind11::arg("v"), pybind11::arg("s"), pybind11::arg("tol"));
		cl.def("getGradient", (const class Teuchos::RCP<const class ROL::Vector<double> > (ROL::TrustRegionModel<double>::*)() const) &ROL::TrustRegionModel<double>::getGradient, "C++: ROL::TrustRegionModel<double>::getGradient() const --> const class Teuchos::RCP<const class ROL::Vector<double> >");
		cl.def("getIterate", (const class Teuchos::RCP<const class ROL::Vector<double> > (ROL::TrustRegionModel<double>::*)() const) &ROL::TrustRegionModel<double>::getIterate, "C++: ROL::TrustRegionModel<double>::getIterate() const --> const class Teuchos::RCP<const class ROL::Vector<double> >");
		cl.def("getObjective", (const class Teuchos::RCP<class ROL::Objective<double> > (ROL::TrustRegionModel<double>::*)() const) &ROL::TrustRegionModel<double>::getObjective, "C++: ROL::TrustRegionModel<double>::getObjective() const --> const class Teuchos::RCP<class ROL::Objective<double> >");
		cl.def("getBoundConstraint", (const class Teuchos::RCP<class ROL::BoundConstraint<double> > (ROL::TrustRegionModel<double>::*)() const) &ROL::TrustRegionModel<double>::getBoundConstraint, "C++: ROL::TrustRegionModel<double>::getBoundConstraint() const --> const class Teuchos::RCP<class ROL::BoundConstraint<double> >");
		cl.def("dualTransform", (void (ROL::TrustRegionModel<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::TrustRegionModel<double>::dualTransform, "C++: ROL::TrustRegionModel<double>::dualTransform(class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("tv"), pybind11::arg("v"));
		cl.def("primalTransform", (void (ROL::TrustRegionModel<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::TrustRegionModel<double>::primalTransform, "C++: ROL::TrustRegionModel<double>::primalTransform(class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("tv"), pybind11::arg("v"));
		cl.def("updatePredictedReduction", (void (ROL::TrustRegionModel<double>::*)(double &, const class ROL::Vector<double> &)) &ROL::TrustRegionModel<double>::updatePredictedReduction, "C++: ROL::TrustRegionModel<double>::updatePredictedReduction(double &, const class ROL::Vector<double> &) --> void", pybind11::arg("pred"), pybind11::arg("s"));
		cl.def("updateActualReduction", (void (ROL::TrustRegionModel<double>::*)(double &, const class ROL::Vector<double> &)) &ROL::TrustRegionModel<double>::updateActualReduction, "C++: ROL::TrustRegionModel<double>::updateActualReduction(double &, const class ROL::Vector<double> &) --> void", pybind11::arg("ared"), pybind11::arg("s"));
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
