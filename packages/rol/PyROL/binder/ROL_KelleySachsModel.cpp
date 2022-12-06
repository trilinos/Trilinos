#include <ROL_BoundConstraint.hpp>
#include <ROL_CauchyPoint.hpp>
#include <ROL_Constraint.hpp>
#include <ROL_DogLeg.hpp>
#include <ROL_DoubleDogLeg.hpp>
#include <ROL_Elementwise_Function.hpp>
#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_KelleySachsModel.hpp>
#include <ROL_LinMore.hpp>
#include <ROL_LinMoreModel.hpp>
#include <ROL_Objective.hpp>
#include <ROL_Secant.hpp>
#include <ROL_Step.hpp>
#include <ROL_TruncatedCG.hpp>
#include <ROL_TrustRegion.hpp>
#include <ROL_TrustRegionFactory.hpp>
#include <ROL_TrustRegionModel.hpp>
#include <ROL_TrustRegionStep.hpp>
#include <ROL_TrustRegionTypes.hpp>
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

// ROL::KelleySachsModel file:ROL_KelleySachsModel.hpp line:61
struct PyCallBack_ROL_KelleySachsModel_double_t : public ROL::KelleySachsModel<double> {
	using ROL::KelleySachsModel<double>::KelleySachsModel;

	double value(const class ROL::Vector<double> & a0, double & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::KelleySachsModel<double> *>(this), "value");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return KelleySachsModel::value(a0, a1);
	}
	void gradient(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double & a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::KelleySachsModel<double> *>(this), "gradient");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return KelleySachsModel::gradient(a0, a1, a2);
	}
	void hessVec(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::KelleySachsModel<double> *>(this), "hessVec");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return KelleySachsModel::hessVec(a0, a1, a2, a3);
	}
	void invHessVec(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::KelleySachsModel<double> *>(this), "invHessVec");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return KelleySachsModel::invHessVec(a0, a1, a2, a3);
	}
	void precond(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::KelleySachsModel<double> *>(this), "precond");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return KelleySachsModel::precond(a0, a1, a2, a3);
	}
	void dualTransform(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::KelleySachsModel<double> *>(this), "dualTransform");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return KelleySachsModel::dualTransform(a0, a1);
	}
	void primalTransform(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::KelleySachsModel<double> *>(this), "primalTransform");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return KelleySachsModel::primalTransform(a0, a1);
	}
	void update(class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class Teuchos::RCP<class ROL::Secant<double> > & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::KelleySachsModel<double> *>(this), "update");
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
	const class Teuchos::RCP<const class ROL::Vector<double> > getGradient() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::KelleySachsModel<double> *>(this), "getGradient");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::KelleySachsModel<double> *>(this), "getIterate");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::KelleySachsModel<double> *>(this), "getObjective");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::KelleySachsModel<double> *>(this), "getBoundConstraint");
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
	void updatePredictedReduction(double & a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::KelleySachsModel<double> *>(this), "updatePredictedReduction");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::KelleySachsModel<double> *>(this), "updateActualReduction");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::KelleySachsModel<double> *>(this), "update");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::KelleySachsModel<double> *>(this), "update");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::KelleySachsModel<double> *>(this), "dirDeriv");
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

// ROL::TrustRegion file:ROL_TrustRegion.hpp line:60
struct PyCallBack_ROL_TrustRegion_double_t : public ROL::TrustRegion<double> {
	using ROL::TrustRegion<double>::TrustRegion;

	void initialize(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegion<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TrustRegion::initialize(a0, a1, a2);
	}
	void update(class ROL::Vector<double> & a0, double & a1, double & a2, int & a3, int & a4, enum ROL::ETrustRegionFlag & a5, const class ROL::Vector<double> & a6, const double a7, const double a8, const class ROL::Vector<double> & a9, int a10, class ROL::Objective<double> & a11, class ROL::BoundConstraint<double> & a12, class ROL::TrustRegionModel<double> & a13) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegion<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TrustRegion::update(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13);
	}
	void run(class ROL::Vector<double> & a0, double & a1, int & a2, int & a3, const double a4, class ROL::TrustRegionModel<double> & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegion<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"TrustRegion::run\"");
	}
};

// ROL::CauchyPoint file:ROL_CauchyPoint.hpp line:65
struct PyCallBack_ROL_CauchyPoint_double_t : public ROL::CauchyPoint<double> {
	using ROL::CauchyPoint<double>::CauchyPoint;

	void initialize(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::CauchyPoint<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return CauchyPoint::initialize(a0, a1, a2);
	}
	void run(class ROL::Vector<double> & a0, double & a1, int & a2, int & a3, const double a4, class ROL::TrustRegionModel<double> & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::CauchyPoint<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return CauchyPoint::run(a0, a1, a2, a3, a4, a5);
	}
	void update(class ROL::Vector<double> & a0, double & a1, double & a2, int & a3, int & a4, enum ROL::ETrustRegionFlag & a5, const class ROL::Vector<double> & a6, const double a7, const double a8, const class ROL::Vector<double> & a9, int a10, class ROL::Objective<double> & a11, class ROL::BoundConstraint<double> & a12, class ROL::TrustRegionModel<double> & a13) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::CauchyPoint<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TrustRegion::update(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13);
	}
};

// ROL::DogLeg file:ROL_DogLeg.hpp line:57
struct PyCallBack_ROL_DogLeg_double_t : public ROL::DogLeg<double> {
	using ROL::DogLeg<double>::DogLeg;

	void initialize(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::DogLeg<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return DogLeg::initialize(a0, a1, a2);
	}
	void run(class ROL::Vector<double> & a0, double & a1, int & a2, int & a3, const double a4, class ROL::TrustRegionModel<double> & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::DogLeg<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return DogLeg::run(a0, a1, a2, a3, a4, a5);
	}
	void update(class ROL::Vector<double> & a0, double & a1, double & a2, int & a3, int & a4, enum ROL::ETrustRegionFlag & a5, const class ROL::Vector<double> & a6, const double a7, const double a8, const class ROL::Vector<double> & a9, int a10, class ROL::Objective<double> & a11, class ROL::BoundConstraint<double> & a12, class ROL::TrustRegionModel<double> & a13) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::DogLeg<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TrustRegion::update(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13);
	}
};

// ROL::DoubleDogLeg file:ROL_DoubleDogLeg.hpp line:57
struct PyCallBack_ROL_DoubleDogLeg_double_t : public ROL::DoubleDogLeg<double> {
	using ROL::DoubleDogLeg<double>::DoubleDogLeg;

	void initialize(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::DoubleDogLeg<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return DoubleDogLeg::initialize(a0, a1, a2);
	}
	void run(class ROL::Vector<double> & a0, double & a1, int & a2, int & a3, const double a4, class ROL::TrustRegionModel<double> & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::DoubleDogLeg<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return DoubleDogLeg::run(a0, a1, a2, a3, a4, a5);
	}
	void update(class ROL::Vector<double> & a0, double & a1, double & a2, int & a3, int & a4, enum ROL::ETrustRegionFlag & a5, const class ROL::Vector<double> & a6, const double a7, const double a8, const class ROL::Vector<double> & a9, int a10, class ROL::Objective<double> & a11, class ROL::BoundConstraint<double> & a12, class ROL::TrustRegionModel<double> & a13) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::DoubleDogLeg<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TrustRegion::update(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13);
	}
};

// ROL::TruncatedCG file:ROL_TruncatedCG.hpp line:57
struct PyCallBack_ROL_TruncatedCG_double_t : public ROL::TruncatedCG<double> {
	using ROL::TruncatedCG<double>::TruncatedCG;

	void initialize(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TruncatedCG<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TruncatedCG::initialize(a0, a1, a2);
	}
	void run(class ROL::Vector<double> & a0, double & a1, int & a2, int & a3, const double a4, class ROL::TrustRegionModel<double> & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TruncatedCG<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TruncatedCG::run(a0, a1, a2, a3, a4, a5);
	}
	void update(class ROL::Vector<double> & a0, double & a1, double & a2, int & a3, int & a4, enum ROL::ETrustRegionFlag & a5, const class ROL::Vector<double> & a6, const double a7, const double a8, const class ROL::Vector<double> & a9, int a10, class ROL::Objective<double> & a11, class ROL::BoundConstraint<double> & a12, class ROL::TrustRegionModel<double> & a13) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TruncatedCG<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TrustRegion::update(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13);
	}
};

// ROL::LinMoreModel file:ROL_LinMoreModel.hpp line:61
struct PyCallBack_ROL_LinMoreModel_double_t : public ROL::LinMoreModel<double> {
	using ROL::LinMoreModel<double>::LinMoreModel;

	void update(class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class Teuchos::RCP<class ROL::Secant<double> > & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LinMoreModel<double> *>(this), "update");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LinMoreModel<double> *>(this), "value");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LinMoreModel<double> *>(this), "gradient");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LinMoreModel<double> *>(this), "hessVec");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LinMoreModel<double> *>(this), "invHessVec");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LinMoreModel<double> *>(this), "precond");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LinMoreModel<double> *>(this), "getGradient");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LinMoreModel<double> *>(this), "getIterate");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LinMoreModel<double> *>(this), "getObjective");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LinMoreModel<double> *>(this), "getBoundConstraint");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LinMoreModel<double> *>(this), "dualTransform");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LinMoreModel<double> *>(this), "primalTransform");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LinMoreModel<double> *>(this), "updatePredictedReduction");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LinMoreModel<double> *>(this), "updateActualReduction");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LinMoreModel<double> *>(this), "update");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LinMoreModel<double> *>(this), "update");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LinMoreModel<double> *>(this), "dirDeriv");
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

// ROL::LinMore file:ROL_LinMore.hpp line:59
struct PyCallBack_ROL_LinMore_double_t : public ROL::LinMore<double> {
	using ROL::LinMore<double>::LinMore;

	void initialize(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LinMore<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return LinMore::initialize(a0, a1, a2);
	}
	void run(class ROL::Vector<double> & a0, double & a1, int & a2, int & a3, const double a4, class ROL::TrustRegionModel<double> & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LinMore<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return LinMore::run(a0, a1, a2, a3, a4, a5);
	}
	void update(class ROL::Vector<double> & a0, double & a1, double & a2, int & a3, int & a4, enum ROL::ETrustRegionFlag & a5, const class ROL::Vector<double> & a6, const double a7, const double a8, const class ROL::Vector<double> & a9, int a10, class ROL::Objective<double> & a11, class ROL::BoundConstraint<double> & a12, class ROL::TrustRegionModel<double> & a13) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::LinMore<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TrustRegion::update(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13);
	}
};

// ROL::TrustRegionStep file: line:35
struct PyCallBack_ROL_TrustRegionStep_double_t : public ROL::TrustRegionStep<double> {
	using ROL::TrustRegionStep<double>::TrustRegionStep;

	void initialize(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, class ROL::Objective<double> & a3, class ROL::BoundConstraint<double> & a4, struct ROL::AlgorithmState<double> & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionStep<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TrustRegionStep::initialize(a0, a1, a2, a3, a4, a5);
	}
	void compute(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, struct ROL::AlgorithmState<double> & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionStep<double> *>(this), "compute");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TrustRegionStep::compute(a0, a1, a2, a3, a4);
	}
	void update(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, struct ROL::AlgorithmState<double> & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionStep<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TrustRegionStep::update(a0, a1, a2, a3, a4);
	}
	std::string printHeader() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionStep<double> *>(this), "printHeader");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<std::string>::value) {
				static pybind11::detail::override_caster_t<std::string> caster;
				return pybind11::detail::cast_ref<std::string>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<std::string>(std::move(o));
		}
		return TrustRegionStep::printHeader();
	}
	std::string printName() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionStep<double> *>(this), "printName");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<std::string>::value) {
				static pybind11::detail::override_caster_t<std::string> caster;
				return pybind11::detail::cast_ref<std::string>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<std::string>(std::move(o));
		}
		return TrustRegionStep::printName();
	}
	std::string print(struct ROL::AlgorithmState<double> & a0, bool a1) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionStep<double> *>(this), "print");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<std::string>::value) {
				static pybind11::detail::override_caster_t<std::string> caster;
				return pybind11::detail::cast_ref<std::string>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<std::string>(std::move(o));
		}
		return TrustRegionStep::print(a0, a1);
	}
	void initialize(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, struct ROL::AlgorithmState<double> & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionStep<double> *>(this), "initialize");
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
	void initialize(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, class ROL::Objective<double> & a4, class ROL::Constraint<double> & a5, struct ROL::AlgorithmState<double> & a6) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionStep<double> *>(this), "initialize");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionStep<double> *>(this), "initialize");
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
	void compute(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, class ROL::Objective<double> & a3, class ROL::Constraint<double> & a4, struct ROL::AlgorithmState<double> & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionStep<double> *>(this), "compute");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionStep<double> *>(this), "update");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionStep<double> *>(this), "compute");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TrustRegionStep<double> *>(this), "update");
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
};

void bind_ROL_KelleySachsModel(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // ROL::KelleySachsModel file:ROL_KelleySachsModel.hpp line:61
		pybind11::class_<ROL::KelleySachsModel<double>, Teuchos::RCP<ROL::KelleySachsModel<double>>, PyCallBack_ROL_KelleySachsModel_double_t, ROL::TrustRegionModel<double>> cl(M("ROL"), "KelleySachsModel_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3){ return new ROL::KelleySachsModel<double>(a0, a1, a2, a3); }, [](class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3){ return new PyCallBack_ROL_KelleySachsModel_double_t(a0, a1, a2, a3); } ), "doc");
		cl.def( pybind11::init( [](class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class Teuchos::RCP<class ROL::Secant<double> > & a4){ return new ROL::KelleySachsModel<double>(a0, a1, a2, a3, a4); }, [](class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class Teuchos::RCP<class ROL::Secant<double> > & a4){ return new PyCallBack_ROL_KelleySachsModel_double_t(a0, a1, a2, a3, a4); } ), "doc");
		cl.def( pybind11::init( [](class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class Teuchos::RCP<class ROL::Secant<double> > & a4, const bool & a5){ return new ROL::KelleySachsModel<double>(a0, a1, a2, a3, a4, a5); }, [](class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class Teuchos::RCP<class ROL::Secant<double> > & a4, const bool & a5){ return new PyCallBack_ROL_KelleySachsModel_double_t(a0, a1, a2, a3, a4, a5); } ), "doc");
		cl.def( pybind11::init<class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class Teuchos::RCP<class ROL::Secant<double> > &, const bool, const bool>(), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("secant"), pybind11::arg("useSecantPrecond"), pybind11::arg("useSecantHessVec") );

		cl.def( pybind11::init( [](PyCallBack_ROL_KelleySachsModel_double_t const &o){ return new PyCallBack_ROL_KelleySachsModel_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::KelleySachsModel<double> const &o){ return new ROL::KelleySachsModel<double>(o); } ) );
		cl.def("setEpsilon", (void (ROL::KelleySachsModel<double>::*)(const double)) &ROL::KelleySachsModel<double>::setEpsilon, "C++: ROL::KelleySachsModel<double>::setEpsilon(const double) --> void", pybind11::arg("eps"));
		cl.def("value", (double (ROL::KelleySachsModel<double>::*)(const class ROL::Vector<double> &, double &)) &ROL::KelleySachsModel<double>::value, "C++: ROL::KelleySachsModel<double>::value(const class ROL::Vector<double> &, double &) --> double", pybind11::arg("s"), pybind11::arg("tol"));
		cl.def("gradient", (void (ROL::KelleySachsModel<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::KelleySachsModel<double>::gradient, "C++: ROL::KelleySachsModel<double>::gradient(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("g"), pybind11::arg("s"), pybind11::arg("tol"));
		cl.def("hessVec", (void (ROL::KelleySachsModel<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::KelleySachsModel<double>::hessVec, "C++: ROL::KelleySachsModel<double>::hessVec(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("Hv"), pybind11::arg("v"), pybind11::arg("s"), pybind11::arg("tol"));
		cl.def("invHessVec", (void (ROL::KelleySachsModel<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::KelleySachsModel<double>::invHessVec, "C++: ROL::KelleySachsModel<double>::invHessVec(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("Hv"), pybind11::arg("v"), pybind11::arg("s"), pybind11::arg("tol"));
		cl.def("precond", (void (ROL::KelleySachsModel<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::KelleySachsModel<double>::precond, "C++: ROL::KelleySachsModel<double>::precond(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("Mv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("dualTransform", (void (ROL::KelleySachsModel<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::KelleySachsModel<double>::dualTransform, "C++: ROL::KelleySachsModel<double>::dualTransform(class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("tv"), pybind11::arg("v"));
		cl.def("primalTransform", (void (ROL::KelleySachsModel<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::KelleySachsModel<double>::primalTransform, "C++: ROL::KelleySachsModel<double>::primalTransform(class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("tv"), pybind11::arg("v"));
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
	{ // ROL::TrustRegion file:ROL_TrustRegion.hpp line:60
		pybind11::class_<ROL::TrustRegion<double>, Teuchos::RCP<ROL::TrustRegion<double>>, PyCallBack_ROL_TrustRegion_double_t> cl(M("ROL"), "TrustRegion_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("parlist") );

		cl.def(pybind11::init<PyCallBack_ROL_TrustRegion_double_t const &>());
		cl.def("initialize", (void (ROL::TrustRegion<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::TrustRegion<double>::initialize, "C++: ROL::TrustRegion<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("s"), pybind11::arg("g"));
		cl.def("update", (void (ROL::TrustRegion<double>::*)(class ROL::Vector<double> &, double &, double &, int &, int &, enum ROL::ETrustRegionFlag &, const class ROL::Vector<double> &, const double, const double, const class ROL::Vector<double> &, int, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::TrustRegionModel<double> &)) &ROL::TrustRegion<double>::update, "C++: ROL::TrustRegion<double>::update(class ROL::Vector<double> &, double &, double &, int &, int &, enum ROL::ETrustRegionFlag &, const class ROL::Vector<double> &, const double, const double, const class ROL::Vector<double> &, int, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::TrustRegionModel<double> &) --> void", pybind11::arg("x"), pybind11::arg("fnew"), pybind11::arg("del"), pybind11::arg("nfval"), pybind11::arg("ngrad"), pybind11::arg("flagTR"), pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("fold"), pybind11::arg("g"), pybind11::arg("iter"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("model"));
		cl.def("run", (void (ROL::TrustRegion<double>::*)(class ROL::Vector<double> &, double &, int &, int &, const double, class ROL::TrustRegionModel<double> &)) &ROL::TrustRegion<double>::run, "C++: ROL::TrustRegion<double>::run(class ROL::Vector<double> &, double &, int &, int &, const double, class ROL::TrustRegionModel<double> &) --> void", pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("iflag"), pybind11::arg("iter"), pybind11::arg("del"), pybind11::arg("model"));
		cl.def("setPredictedReduction", (void (ROL::TrustRegion<double>::*)(const double)) &ROL::TrustRegion<double>::setPredictedReduction, "C++: ROL::TrustRegion<double>::setPredictedReduction(const double) --> void", pybind11::arg("pRed"));
		cl.def("getPredictedReduction", (double (ROL::TrustRegion<double>::*)() const) &ROL::TrustRegion<double>::getPredictedReduction, "C++: ROL::TrustRegion<double>::getPredictedReduction() const --> double");
		cl.def("assign", (class ROL::TrustRegion<double> & (ROL::TrustRegion<double>::*)(const class ROL::TrustRegion<double> &)) &ROL::TrustRegion<double>::operator=, "C++: ROL::TrustRegion<double>::operator=(const class ROL::TrustRegion<double> &) --> class ROL::TrustRegion<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::CauchyPoint file:ROL_CauchyPoint.hpp line:65
		pybind11::class_<ROL::CauchyPoint<double>, Teuchos::RCP<ROL::CauchyPoint<double>>, PyCallBack_ROL_CauchyPoint_double_t, ROL::TrustRegion<double>> cl(M("ROL"), "CauchyPoint_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("parlist") );

		cl.def( pybind11::init( [](PyCallBack_ROL_CauchyPoint_double_t const &o){ return new PyCallBack_ROL_CauchyPoint_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::CauchyPoint<double> const &o){ return new ROL::CauchyPoint<double>(o); } ) );
		cl.def("initialize", (void (ROL::CauchyPoint<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::CauchyPoint<double>::initialize, "C++: ROL::CauchyPoint<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("s"), pybind11::arg("g"));
		cl.def("run", (void (ROL::CauchyPoint<double>::*)(class ROL::Vector<double> &, double &, int &, int &, const double, class ROL::TrustRegionModel<double> &)) &ROL::CauchyPoint<double>::run, "C++: ROL::CauchyPoint<double>::run(class ROL::Vector<double> &, double &, int &, int &, const double, class ROL::TrustRegionModel<double> &) --> void", pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("iflag"), pybind11::arg("iter"), pybind11::arg("del"), pybind11::arg("model"));
		cl.def("assign", (class ROL::CauchyPoint<double> & (ROL::CauchyPoint<double>::*)(const class ROL::CauchyPoint<double> &)) &ROL::CauchyPoint<double>::operator=, "C++: ROL::CauchyPoint<double>::operator=(const class ROL::CauchyPoint<double> &) --> class ROL::CauchyPoint<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("initialize", (void (ROL::TrustRegion<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::TrustRegion<double>::initialize, "C++: ROL::TrustRegion<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("s"), pybind11::arg("g"));
		cl.def("update", (void (ROL::TrustRegion<double>::*)(class ROL::Vector<double> &, double &, double &, int &, int &, enum ROL::ETrustRegionFlag &, const class ROL::Vector<double> &, const double, const double, const class ROL::Vector<double> &, int, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::TrustRegionModel<double> &)) &ROL::TrustRegion<double>::update, "C++: ROL::TrustRegion<double>::update(class ROL::Vector<double> &, double &, double &, int &, int &, enum ROL::ETrustRegionFlag &, const class ROL::Vector<double> &, const double, const double, const class ROL::Vector<double> &, int, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::TrustRegionModel<double> &) --> void", pybind11::arg("x"), pybind11::arg("fnew"), pybind11::arg("del"), pybind11::arg("nfval"), pybind11::arg("ngrad"), pybind11::arg("flagTR"), pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("fold"), pybind11::arg("g"), pybind11::arg("iter"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("model"));
		cl.def("run", (void (ROL::TrustRegion<double>::*)(class ROL::Vector<double> &, double &, int &, int &, const double, class ROL::TrustRegionModel<double> &)) &ROL::TrustRegion<double>::run, "C++: ROL::TrustRegion<double>::run(class ROL::Vector<double> &, double &, int &, int &, const double, class ROL::TrustRegionModel<double> &) --> void", pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("iflag"), pybind11::arg("iter"), pybind11::arg("del"), pybind11::arg("model"));
		cl.def("setPredictedReduction", (void (ROL::TrustRegion<double>::*)(const double)) &ROL::TrustRegion<double>::setPredictedReduction, "C++: ROL::TrustRegion<double>::setPredictedReduction(const double) --> void", pybind11::arg("pRed"));
		cl.def("getPredictedReduction", (double (ROL::TrustRegion<double>::*)() const) &ROL::TrustRegion<double>::getPredictedReduction, "C++: ROL::TrustRegion<double>::getPredictedReduction() const --> double");
		cl.def("assign", (class ROL::TrustRegion<double> & (ROL::TrustRegion<double>::*)(const class ROL::TrustRegion<double> &)) &ROL::TrustRegion<double>::operator=, "C++: ROL::TrustRegion<double>::operator=(const class ROL::TrustRegion<double> &) --> class ROL::TrustRegion<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::DogLeg file:ROL_DogLeg.hpp line:57
		pybind11::class_<ROL::DogLeg<double>, Teuchos::RCP<ROL::DogLeg<double>>, PyCallBack_ROL_DogLeg_double_t, ROL::TrustRegion<double>> cl(M("ROL"), "DogLeg_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("parlist") );

		cl.def( pybind11::init( [](PyCallBack_ROL_DogLeg_double_t const &o){ return new PyCallBack_ROL_DogLeg_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::DogLeg<double> const &o){ return new ROL::DogLeg<double>(o); } ) );
		cl.def("initialize", (void (ROL::DogLeg<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::DogLeg<double>::initialize, "C++: ROL::DogLeg<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("s"), pybind11::arg("g"));
		cl.def("run", (void (ROL::DogLeg<double>::*)(class ROL::Vector<double> &, double &, int &, int &, const double, class ROL::TrustRegionModel<double> &)) &ROL::DogLeg<double>::run, "C++: ROL::DogLeg<double>::run(class ROL::Vector<double> &, double &, int &, int &, const double, class ROL::TrustRegionModel<double> &) --> void", pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("iflag"), pybind11::arg("iter"), pybind11::arg("del"), pybind11::arg("model"));
		cl.def("assign", (class ROL::DogLeg<double> & (ROL::DogLeg<double>::*)(const class ROL::DogLeg<double> &)) &ROL::DogLeg<double>::operator=, "C++: ROL::DogLeg<double>::operator=(const class ROL::DogLeg<double> &) --> class ROL::DogLeg<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("initialize", (void (ROL::TrustRegion<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::TrustRegion<double>::initialize, "C++: ROL::TrustRegion<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("s"), pybind11::arg("g"));
		cl.def("update", (void (ROL::TrustRegion<double>::*)(class ROL::Vector<double> &, double &, double &, int &, int &, enum ROL::ETrustRegionFlag &, const class ROL::Vector<double> &, const double, const double, const class ROL::Vector<double> &, int, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::TrustRegionModel<double> &)) &ROL::TrustRegion<double>::update, "C++: ROL::TrustRegion<double>::update(class ROL::Vector<double> &, double &, double &, int &, int &, enum ROL::ETrustRegionFlag &, const class ROL::Vector<double> &, const double, const double, const class ROL::Vector<double> &, int, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::TrustRegionModel<double> &) --> void", pybind11::arg("x"), pybind11::arg("fnew"), pybind11::arg("del"), pybind11::arg("nfval"), pybind11::arg("ngrad"), pybind11::arg("flagTR"), pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("fold"), pybind11::arg("g"), pybind11::arg("iter"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("model"));
		cl.def("run", (void (ROL::TrustRegion<double>::*)(class ROL::Vector<double> &, double &, int &, int &, const double, class ROL::TrustRegionModel<double> &)) &ROL::TrustRegion<double>::run, "C++: ROL::TrustRegion<double>::run(class ROL::Vector<double> &, double &, int &, int &, const double, class ROL::TrustRegionModel<double> &) --> void", pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("iflag"), pybind11::arg("iter"), pybind11::arg("del"), pybind11::arg("model"));
		cl.def("setPredictedReduction", (void (ROL::TrustRegion<double>::*)(const double)) &ROL::TrustRegion<double>::setPredictedReduction, "C++: ROL::TrustRegion<double>::setPredictedReduction(const double) --> void", pybind11::arg("pRed"));
		cl.def("getPredictedReduction", (double (ROL::TrustRegion<double>::*)() const) &ROL::TrustRegion<double>::getPredictedReduction, "C++: ROL::TrustRegion<double>::getPredictedReduction() const --> double");
		cl.def("assign", (class ROL::TrustRegion<double> & (ROL::TrustRegion<double>::*)(const class ROL::TrustRegion<double> &)) &ROL::TrustRegion<double>::operator=, "C++: ROL::TrustRegion<double>::operator=(const class ROL::TrustRegion<double> &) --> class ROL::TrustRegion<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::DoubleDogLeg file:ROL_DoubleDogLeg.hpp line:57
		pybind11::class_<ROL::DoubleDogLeg<double>, Teuchos::RCP<ROL::DoubleDogLeg<double>>, PyCallBack_ROL_DoubleDogLeg_double_t, ROL::TrustRegion<double>> cl(M("ROL"), "DoubleDogLeg_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("parlist") );

		cl.def( pybind11::init( [](PyCallBack_ROL_DoubleDogLeg_double_t const &o){ return new PyCallBack_ROL_DoubleDogLeg_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::DoubleDogLeg<double> const &o){ return new ROL::DoubleDogLeg<double>(o); } ) );
		cl.def("initialize", (void (ROL::DoubleDogLeg<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::DoubleDogLeg<double>::initialize, "C++: ROL::DoubleDogLeg<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("s"), pybind11::arg("g"));
		cl.def("run", (void (ROL::DoubleDogLeg<double>::*)(class ROL::Vector<double> &, double &, int &, int &, const double, class ROL::TrustRegionModel<double> &)) &ROL::DoubleDogLeg<double>::run, "C++: ROL::DoubleDogLeg<double>::run(class ROL::Vector<double> &, double &, int &, int &, const double, class ROL::TrustRegionModel<double> &) --> void", pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("iflag"), pybind11::arg("iter"), pybind11::arg("del"), pybind11::arg("model"));
		cl.def("assign", (class ROL::DoubleDogLeg<double> & (ROL::DoubleDogLeg<double>::*)(const class ROL::DoubleDogLeg<double> &)) &ROL::DoubleDogLeg<double>::operator=, "C++: ROL::DoubleDogLeg<double>::operator=(const class ROL::DoubleDogLeg<double> &) --> class ROL::DoubleDogLeg<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("initialize", (void (ROL::TrustRegion<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::TrustRegion<double>::initialize, "C++: ROL::TrustRegion<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("s"), pybind11::arg("g"));
		cl.def("update", (void (ROL::TrustRegion<double>::*)(class ROL::Vector<double> &, double &, double &, int &, int &, enum ROL::ETrustRegionFlag &, const class ROL::Vector<double> &, const double, const double, const class ROL::Vector<double> &, int, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::TrustRegionModel<double> &)) &ROL::TrustRegion<double>::update, "C++: ROL::TrustRegion<double>::update(class ROL::Vector<double> &, double &, double &, int &, int &, enum ROL::ETrustRegionFlag &, const class ROL::Vector<double> &, const double, const double, const class ROL::Vector<double> &, int, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::TrustRegionModel<double> &) --> void", pybind11::arg("x"), pybind11::arg("fnew"), pybind11::arg("del"), pybind11::arg("nfval"), pybind11::arg("ngrad"), pybind11::arg("flagTR"), pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("fold"), pybind11::arg("g"), pybind11::arg("iter"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("model"));
		cl.def("run", (void (ROL::TrustRegion<double>::*)(class ROL::Vector<double> &, double &, int &, int &, const double, class ROL::TrustRegionModel<double> &)) &ROL::TrustRegion<double>::run, "C++: ROL::TrustRegion<double>::run(class ROL::Vector<double> &, double &, int &, int &, const double, class ROL::TrustRegionModel<double> &) --> void", pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("iflag"), pybind11::arg("iter"), pybind11::arg("del"), pybind11::arg("model"));
		cl.def("setPredictedReduction", (void (ROL::TrustRegion<double>::*)(const double)) &ROL::TrustRegion<double>::setPredictedReduction, "C++: ROL::TrustRegion<double>::setPredictedReduction(const double) --> void", pybind11::arg("pRed"));
		cl.def("getPredictedReduction", (double (ROL::TrustRegion<double>::*)() const) &ROL::TrustRegion<double>::getPredictedReduction, "C++: ROL::TrustRegion<double>::getPredictedReduction() const --> double");
		cl.def("assign", (class ROL::TrustRegion<double> & (ROL::TrustRegion<double>::*)(const class ROL::TrustRegion<double> &)) &ROL::TrustRegion<double>::operator=, "C++: ROL::TrustRegion<double>::operator=(const class ROL::TrustRegion<double> &) --> class ROL::TrustRegion<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::TruncatedCG file:ROL_TruncatedCG.hpp line:57
		pybind11::class_<ROL::TruncatedCG<double>, Teuchos::RCP<ROL::TruncatedCG<double>>, PyCallBack_ROL_TruncatedCG_double_t, ROL::TrustRegion<double>> cl(M("ROL"), "TruncatedCG_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("parlist") );

		cl.def( pybind11::init( [](PyCallBack_ROL_TruncatedCG_double_t const &o){ return new PyCallBack_ROL_TruncatedCG_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::TruncatedCG<double> const &o){ return new ROL::TruncatedCG<double>(o); } ) );
		cl.def("initialize", (void (ROL::TruncatedCG<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::TruncatedCG<double>::initialize, "C++: ROL::TruncatedCG<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("s"), pybind11::arg("g"));
		cl.def("run", (void (ROL::TruncatedCG<double>::*)(class ROL::Vector<double> &, double &, int &, int &, const double, class ROL::TrustRegionModel<double> &)) &ROL::TruncatedCG<double>::run, "C++: ROL::TruncatedCG<double>::run(class ROL::Vector<double> &, double &, int &, int &, const double, class ROL::TrustRegionModel<double> &) --> void", pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("iflag"), pybind11::arg("iter"), pybind11::arg("del"), pybind11::arg("model"));
		cl.def("assign", (class ROL::TruncatedCG<double> & (ROL::TruncatedCG<double>::*)(const class ROL::TruncatedCG<double> &)) &ROL::TruncatedCG<double>::operator=, "C++: ROL::TruncatedCG<double>::operator=(const class ROL::TruncatedCG<double> &) --> class ROL::TruncatedCG<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("initialize", (void (ROL::TrustRegion<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::TrustRegion<double>::initialize, "C++: ROL::TrustRegion<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("s"), pybind11::arg("g"));
		cl.def("update", (void (ROL::TrustRegion<double>::*)(class ROL::Vector<double> &, double &, double &, int &, int &, enum ROL::ETrustRegionFlag &, const class ROL::Vector<double> &, const double, const double, const class ROL::Vector<double> &, int, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::TrustRegionModel<double> &)) &ROL::TrustRegion<double>::update, "C++: ROL::TrustRegion<double>::update(class ROL::Vector<double> &, double &, double &, int &, int &, enum ROL::ETrustRegionFlag &, const class ROL::Vector<double> &, const double, const double, const class ROL::Vector<double> &, int, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::TrustRegionModel<double> &) --> void", pybind11::arg("x"), pybind11::arg("fnew"), pybind11::arg("del"), pybind11::arg("nfval"), pybind11::arg("ngrad"), pybind11::arg("flagTR"), pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("fold"), pybind11::arg("g"), pybind11::arg("iter"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("model"));
		cl.def("run", (void (ROL::TrustRegion<double>::*)(class ROL::Vector<double> &, double &, int &, int &, const double, class ROL::TrustRegionModel<double> &)) &ROL::TrustRegion<double>::run, "C++: ROL::TrustRegion<double>::run(class ROL::Vector<double> &, double &, int &, int &, const double, class ROL::TrustRegionModel<double> &) --> void", pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("iflag"), pybind11::arg("iter"), pybind11::arg("del"), pybind11::arg("model"));
		cl.def("setPredictedReduction", (void (ROL::TrustRegion<double>::*)(const double)) &ROL::TrustRegion<double>::setPredictedReduction, "C++: ROL::TrustRegion<double>::setPredictedReduction(const double) --> void", pybind11::arg("pRed"));
		cl.def("getPredictedReduction", (double (ROL::TrustRegion<double>::*)() const) &ROL::TrustRegion<double>::getPredictedReduction, "C++: ROL::TrustRegion<double>::getPredictedReduction() const --> double");
		cl.def("assign", (class ROL::TrustRegion<double> & (ROL::TrustRegion<double>::*)(const class ROL::TrustRegion<double> &)) &ROL::TrustRegion<double>::operator=, "C++: ROL::TrustRegion<double>::operator=(const class ROL::TrustRegion<double> &) --> class ROL::TrustRegion<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::LinMoreModel file:ROL_LinMoreModel.hpp line:61
		pybind11::class_<ROL::LinMoreModel<double>, Teuchos::RCP<ROL::LinMoreModel<double>>, PyCallBack_ROL_LinMoreModel_double_t, ROL::TrustRegionModel<double>> cl(M("ROL"), "LinMoreModel_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3){ return new ROL::LinMoreModel<double>(a0, a1, a2, a3); }, [](class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3){ return new PyCallBack_ROL_LinMoreModel_double_t(a0, a1, a2, a3); } ), "doc");
		cl.def( pybind11::init( [](class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class Teuchos::RCP<class ROL::Secant<double> > & a4){ return new ROL::LinMoreModel<double>(a0, a1, a2, a3, a4); }, [](class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class Teuchos::RCP<class ROL::Secant<double> > & a4){ return new PyCallBack_ROL_LinMoreModel_double_t(a0, a1, a2, a3, a4); } ), "doc");
		cl.def( pybind11::init( [](class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class Teuchos::RCP<class ROL::Secant<double> > & a4, const bool & a5){ return new ROL::LinMoreModel<double>(a0, a1, a2, a3, a4, a5); }, [](class ROL::Objective<double> & a0, class ROL::BoundConstraint<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class Teuchos::RCP<class ROL::Secant<double> > & a4, const bool & a5){ return new PyCallBack_ROL_LinMoreModel_double_t(a0, a1, a2, a3, a4, a5); } ), "doc");
		cl.def( pybind11::init<class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class Teuchos::RCP<class ROL::Secant<double> > &, const bool, const bool>(), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("secant"), pybind11::arg("useSecantPrecond"), pybind11::arg("useSecantHessVec") );

		cl.def( pybind11::init( [](PyCallBack_ROL_LinMoreModel_double_t const &o){ return new PyCallBack_ROL_LinMoreModel_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::LinMoreModel<double> const &o){ return new ROL::LinMoreModel<double>(o); } ) );
		cl.def("applyFullHessian", (void (ROL::LinMoreModel<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::LinMoreModel<double>::applyFullHessian, "C++: ROL::LinMoreModel<double>::applyFullHessian(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("hv"), pybind11::arg("v"), pybind11::arg("tol"));
		cl.def("applyFreeHessian", (void (ROL::LinMoreModel<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::LinMoreModel<double>::applyFreeHessian, "C++: ROL::LinMoreModel<double>::applyFreeHessian(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("hv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("applyFullPrecond", (void (ROL::LinMoreModel<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::LinMoreModel<double>::applyFullPrecond, "C++: ROL::LinMoreModel<double>::applyFullPrecond(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("pv"), pybind11::arg("v"), pybind11::arg("tol"));
		cl.def("applyFreePrecond", (void (ROL::LinMoreModel<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::LinMoreModel<double>::applyFreePrecond, "C++: ROL::LinMoreModel<double>::applyFreePrecond(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("pv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("tol"));
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
	{ // ROL::LinMore file:ROL_LinMore.hpp line:59
		pybind11::class_<ROL::LinMore<double>, Teuchos::RCP<ROL::LinMore<double>>, PyCallBack_ROL_LinMore_double_t, ROL::TrustRegion<double>> cl(M("ROL"), "LinMore_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("parlist") );

		cl.def( pybind11::init( [](PyCallBack_ROL_LinMore_double_t const &o){ return new PyCallBack_ROL_LinMore_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::LinMore<double> const &o){ return new ROL::LinMore<double>(o); } ) );
		cl.def("initialize", (void (ROL::LinMore<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::LinMore<double>::initialize, "C++: ROL::LinMore<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("s"), pybind11::arg("g"));
		cl.def("run", (void (ROL::LinMore<double>::*)(class ROL::Vector<double> &, double &, int &, int &, const double, class ROL::TrustRegionModel<double> &)) &ROL::LinMore<double>::run, "C++: ROL::LinMore<double>::run(class ROL::Vector<double> &, double &, int &, int &, const double, class ROL::TrustRegionModel<double> &) --> void", pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("iflag"), pybind11::arg("iter"), pybind11::arg("del"), pybind11::arg("model"));
		cl.def("assign", (class ROL::LinMore<double> & (ROL::LinMore<double>::*)(const class ROL::LinMore<double> &)) &ROL::LinMore<double>::operator=, "C++: ROL::LinMore<double>::operator=(const class ROL::LinMore<double> &) --> class ROL::LinMore<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("initialize", (void (ROL::TrustRegion<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::TrustRegion<double>::initialize, "C++: ROL::TrustRegion<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("s"), pybind11::arg("g"));
		cl.def("update", (void (ROL::TrustRegion<double>::*)(class ROL::Vector<double> &, double &, double &, int &, int &, enum ROL::ETrustRegionFlag &, const class ROL::Vector<double> &, const double, const double, const class ROL::Vector<double> &, int, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::TrustRegionModel<double> &)) &ROL::TrustRegion<double>::update, "C++: ROL::TrustRegion<double>::update(class ROL::Vector<double> &, double &, double &, int &, int &, enum ROL::ETrustRegionFlag &, const class ROL::Vector<double> &, const double, const double, const class ROL::Vector<double> &, int, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, class ROL::TrustRegionModel<double> &) --> void", pybind11::arg("x"), pybind11::arg("fnew"), pybind11::arg("del"), pybind11::arg("nfval"), pybind11::arg("ngrad"), pybind11::arg("flagTR"), pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("fold"), pybind11::arg("g"), pybind11::arg("iter"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("model"));
		cl.def("run", (void (ROL::TrustRegion<double>::*)(class ROL::Vector<double> &, double &, int &, int &, const double, class ROL::TrustRegionModel<double> &)) &ROL::TrustRegion<double>::run, "C++: ROL::TrustRegion<double>::run(class ROL::Vector<double> &, double &, int &, int &, const double, class ROL::TrustRegionModel<double> &) --> void", pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("iflag"), pybind11::arg("iter"), pybind11::arg("del"), pybind11::arg("model"));
		cl.def("setPredictedReduction", (void (ROL::TrustRegion<double>::*)(const double)) &ROL::TrustRegion<double>::setPredictedReduction, "C++: ROL::TrustRegion<double>::setPredictedReduction(const double) --> void", pybind11::arg("pRed"));
		cl.def("getPredictedReduction", (double (ROL::TrustRegion<double>::*)() const) &ROL::TrustRegion<double>::getPredictedReduction, "C++: ROL::TrustRegion<double>::getPredictedReduction() const --> double");
		cl.def("assign", (class ROL::TrustRegion<double> & (ROL::TrustRegion<double>::*)(const class ROL::TrustRegion<double> &)) &ROL::TrustRegion<double>::operator=, "C++: ROL::TrustRegion<double>::operator=(const class ROL::TrustRegion<double> &) --> class ROL::TrustRegion<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	// ROL::TrustRegionFactory(class Teuchos::ParameterList &) file:ROL_TrustRegionFactory.hpp line:61
	M("ROL").def("TrustRegionFactory", (class Teuchos::RCP<class ROL::TrustRegion<double> > (*)(class Teuchos::ParameterList &)) &ROL::TrustRegionFactory<double>, "C++: ROL::TrustRegionFactory(class Teuchos::ParameterList &) --> class Teuchos::RCP<class ROL::TrustRegion<double> >", pybind11::arg("parlist"));

	{ // ROL::TrustRegionStep file: line:35
		pybind11::class_<ROL::TrustRegionStep<double>, Teuchos::RCP<ROL::TrustRegionStep<double>>, PyCallBack_ROL_TrustRegionStep_double_t, ROL::Step<double>> cl(M("ROL"), "TrustRegionStep_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("parlist") );

		cl.def( pybind11::init<class Teuchos::RCP<class ROL::Secant<double> > &, class Teuchos::ParameterList &>(), pybind11::arg("secant"), pybind11::arg("parlist") );

		cl.def( pybind11::init( [](PyCallBack_ROL_TrustRegionStep_double_t const &o){ return new PyCallBack_ROL_TrustRegionStep_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::TrustRegionStep<double> const &o){ return new ROL::TrustRegionStep<double>(o); } ) );
		cl.def("initialize", [](ROL::TrustRegionStep<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, class ROL::Objective<double> & a4, class ROL::Constraint<double> & a5, class ROL::BoundConstraint<double> & a6, struct ROL::AlgorithmState<double> & a7) -> void { return o.initialize(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("l"), pybind11::arg("c"), pybind11::arg("obj"), pybind11::arg("con"), pybind11::arg("bnd"), pybind11::arg("algo_state"));
		cl.def("initialize", [](ROL::TrustRegionStep<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, class ROL::Objective<double> & a4, class ROL::Constraint<double> & a5, struct ROL::AlgorithmState<double> & a6) -> void { return o.initialize(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("l"), pybind11::arg("c"), pybind11::arg("obj"), pybind11::arg("con"), pybind11::arg("algo_state"));
		cl.def("initialize", [](ROL::TrustRegionStep<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, class ROL::Objective<double> & a2, class ROL::BoundConstraint<double> & a3, struct ROL::AlgorithmState<double> & a4) -> void { return o.initialize(a0, a1, a2, a3, a4); }, "", pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("con"), pybind11::arg("algo_state"));
		cl.def("compute", [](ROL::TrustRegionStep<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, class ROL::Objective<double> & a3, class ROL::Constraint<double> & a4, class ROL::BoundConstraint<double> & a5, struct ROL::AlgorithmState<double> & a6) -> void { return o.compute(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("s"), pybind11::arg("x"), pybind11::arg("l"), pybind11::arg("obj"), pybind11::arg("con"), pybind11::arg("bnd"), pybind11::arg("algo_state"));
		cl.def("compute", [](ROL::TrustRegionStep<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, class ROL::Objective<double> & a3, class ROL::Constraint<double> & a4, struct ROL::AlgorithmState<double> & a5) -> void { return o.compute(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("s"), pybind11::arg("x"), pybind11::arg("l"), pybind11::arg("obj"), pybind11::arg("con"), pybind11::arg("algo_state"));
		cl.def("update", [](ROL::TrustRegionStep<double> &o, class ROL::Vector<double> & a0, class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, class ROL::Objective<double> & a3, class ROL::Constraint<double> & a4, class ROL::BoundConstraint<double> & a5, struct ROL::AlgorithmState<double> & a6) -> void { return o.update(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("x"), pybind11::arg("l"), pybind11::arg("s"), pybind11::arg("obj"), pybind11::arg("con"), pybind11::arg("bnd"), pybind11::arg("algo_state"));
		cl.def("update", [](ROL::TrustRegionStep<double> &o, class ROL::Vector<double> & a0, class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, class ROL::Objective<double> & a3, class ROL::Constraint<double> & a4, struct ROL::AlgorithmState<double> & a5) -> void { return o.update(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("x"), pybind11::arg("l"), pybind11::arg("s"), pybind11::arg("obj"), pybind11::arg("con"), pybind11::arg("algo_state"));
		cl.def("initialize", (void (ROL::TrustRegionStep<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &)) &ROL::TrustRegionStep<double>::initialize, "Initialize step.\n\n      This function initializes the information necessary to run the trust-region algorithm.\n      \n\n           is the initial guess for the optimization vector.\n      \n\n         is the objective function.\n      \n\n         is the bound constraint.\n      \n\n  is the algorithm state.\n\nC++: ROL::TrustRegionStep<double>::initialize(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &) --> void", pybind11::arg("x"), pybind11::arg("s"), pybind11::arg("g"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("algo_state"));
		cl.def("compute", (void (ROL::TrustRegionStep<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &)) &ROL::TrustRegionStep<double>::compute, "Compute step.\n\n      Computes a trial step, \n by solving the trust-region subproblem.  \n      The trust-region subproblem solver is defined by the enum ETrustRegion.  \n      \n\n          is the computed trial step\n      \n\n          is the current iterate\n      \n\n        is the objective function\n      \n\n        are the bound constraints\n      \n\n contains the current state of the algorithm\n\nC++: ROL::TrustRegionStep<double>::compute(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &) --> void", pybind11::arg("s"), pybind11::arg("x"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("algo_state"));
		cl.def("update", (void (ROL::TrustRegionStep<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &)) &ROL::TrustRegionStep<double>::update, "Update step, if successful.\n\n      Given a trial step, \n, this function updates \n. \n      This function also updates the secant approximation.\n\n      \n          is the updated iterate\n      \n\n          is the computed trial step\n      \n\n        is the objective function\n      \n\n        are the bound constraints\n      \n\n contains the current state of the algorithm\n\nC++: ROL::TrustRegionStep<double>::update(class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &, class ROL::BoundConstraint<double> &, struct ROL::AlgorithmState<double> &) --> void", pybind11::arg("x"), pybind11::arg("s"), pybind11::arg("obj"), pybind11::arg("bnd"), pybind11::arg("algo_state"));
		cl.def("printHeader", (std::string (ROL::TrustRegionStep<double>::*)() const) &ROL::TrustRegionStep<double>::printHeader, "Print iterate header.\n\n      This function produces a string containing header information.\n\nC++: ROL::TrustRegionStep<double>::printHeader() const --> std::string");
		cl.def("printName", (std::string (ROL::TrustRegionStep<double>::*)() const) &ROL::TrustRegionStep<double>::printName, "Print step name.\n\n      This function produces a string containing the algorithmic step information.\n\nC++: ROL::TrustRegionStep<double>::printName() const --> std::string");
		cl.def("print", [](ROL::TrustRegionStep<double> const &o, struct ROL::AlgorithmState<double> & a0) -> std::string { return o.print(a0); }, "", pybind11::arg("algo_state"));
		cl.def("print", (std::string (ROL::TrustRegionStep<double>::*)(struct ROL::AlgorithmState<double> &, bool) const) &ROL::TrustRegionStep<double>::print, "Print iterate status.\n\n      This function prints the iteration status.\n\n      \n    is the current state of the algorithm\n      \n\n   if ste to true will print the header at each iteration\n\nC++: ROL::TrustRegionStep<double>::print(struct ROL::AlgorithmState<double> &, bool) const --> std::string", pybind11::arg("algo_state"), pybind11::arg("print_header"));
		cl.def("assign", (class ROL::TrustRegionStep<double> & (ROL::TrustRegionStep<double>::*)(const class ROL::TrustRegionStep<double> &)) &ROL::TrustRegionStep<double>::operator=, "C++: ROL::TrustRegionStep<double>::operator=(const class ROL::TrustRegionStep<double> &) --> class ROL::TrustRegionStep<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
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
}
