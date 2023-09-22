#include <ROL_Algorithm.hpp>
#include <ROL_BatchManager.hpp>
#include <ROL_BoundConstraint.hpp>
#include <ROL_BoundConstraint_SimOpt.hpp>
#include <ROL_CauchyPoint.hpp>
#include <ROL_Constraint.hpp>
#include <ROL_Constraint_SimOpt.hpp>
#include <ROL_DogLeg.hpp>
#include <ROL_DoubleDogLeg.hpp>
#include <ROL_Elementwise_Function.hpp>
#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_KelleySachsModel.hpp>
#include <ROL_LinMore.hpp>
#include <ROL_LinMoreModel.hpp>
#include <ROL_NonlinearLeastSquaresObjective.hpp>
#include <ROL_Objective.hpp>
#include <ROL_Objective_FSsolver.hpp>
#include <ROL_Secant.hpp>
#include <ROL_Solver.hpp>
#include <ROL_StatusTest.hpp>
#include <ROL_Step.hpp>
#include <ROL_TruncatedCG.hpp>
#include <ROL_TrustRegion.hpp>
#include <ROL_TrustRegionFactory.hpp>
#include <ROL_TrustRegionModel.hpp>
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

PYBIND11_TYPE_CASTER_BASE_HOLDER(ROL::Constraint_SimOpt<double>, Teuchos::RCP<ROL::Constraint_SimOpt<double>>)

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

// ROL::Constraint_SimOpt file:ROL_Constraint_SimOpt.hpp line:108
struct PyCallBack_ROL_Constraint_SimOpt_double_t : public ROL::Constraint_SimOpt<double> {
	using ROL::Constraint_SimOpt<double>::Constraint_SimOpt;

	void update(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, bool a2, int a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_SimOpt::update(a0, a1, a2, a3);
	}
	void update(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, enum ROL::UpdateType a2, int a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_SimOpt::update(a0, a1, a2, a3);
	}
	void update_1(const class ROL::Vector<double> & a0, bool a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "update_1");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_SimOpt::update_1(a0, a1, a2);
	}
	void update_1(const class ROL::Vector<double> & a0, enum ROL::UpdateType a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "update_1");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_SimOpt::update_1(a0, a1, a2);
	}
	void update_2(const class ROL::Vector<double> & a0, bool a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "update_2");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_SimOpt::update_2(a0, a1, a2);
	}
	void update_2(const class ROL::Vector<double> & a0, enum ROL::UpdateType a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "update_2");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_SimOpt::update_2(a0, a1, a2);
	}
	void solve_update(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, enum ROL::UpdateType a2, int a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "solve_update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_SimOpt::solve_update(a0, a1, a2, a3);
	}
	void value(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "value");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"Constraint_SimOpt::value\"");
	}
	void solve(class ROL::Vector<double> & a0, class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "solve");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_SimOpt::solve(a0, a1, a2, a3);
	}
	void setSolveParameters(class Teuchos::ParameterList & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "setSolveParameters");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_SimOpt::setSolveParameters(a0);
	}
	void applyJacobian_1(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, double & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "applyJacobian_1");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_SimOpt::applyJacobian_1(a0, a1, a2, a3, a4);
	}
	void applyJacobian_2(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, double & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "applyJacobian_2");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_SimOpt::applyJacobian_2(a0, a1, a2, a3, a4);
	}
	void applyInverseJacobian_1(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, double & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "applyInverseJacobian_1");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_SimOpt::applyInverseJacobian_1(a0, a1, a2, a3, a4);
	}
	void applyAdjointJacobian_1(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, double & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "applyAdjointJacobian_1");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_SimOpt::applyAdjointJacobian_1(a0, a1, a2, a3, a4);
	}
	void applyAdjointJacobian_1(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class ROL::Vector<double> & a4, double & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "applyAdjointJacobian_1");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_SimOpt::applyAdjointJacobian_1(a0, a1, a2, a3, a4, a5);
	}
	void applyAdjointJacobian_2(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, double & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "applyAdjointJacobian_2");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_SimOpt::applyAdjointJacobian_2(a0, a1, a2, a3, a4);
	}
	void applyAdjointJacobian_2(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class ROL::Vector<double> & a4, double & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "applyAdjointJacobian_2");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_SimOpt::applyAdjointJacobian_2(a0, a1, a2, a3, a4, a5);
	}
	void applyInverseAdjointJacobian_1(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, double & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "applyInverseAdjointJacobian_1");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_SimOpt::applyInverseAdjointJacobian_1(a0, a1, a2, a3, a4);
	}
	void applyAdjointHessian_11(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class ROL::Vector<double> & a4, double & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "applyAdjointHessian_11");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_SimOpt::applyAdjointHessian_11(a0, a1, a2, a3, a4, a5);
	}
	void applyAdjointHessian_12(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class ROL::Vector<double> & a4, double & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "applyAdjointHessian_12");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_SimOpt::applyAdjointHessian_12(a0, a1, a2, a3, a4, a5);
	}
	void applyAdjointHessian_21(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class ROL::Vector<double> & a4, double & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "applyAdjointHessian_21");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_SimOpt::applyAdjointHessian_21(a0, a1, a2, a3, a4, a5);
	}
	void applyAdjointHessian_22(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class ROL::Vector<double> & a4, double & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "applyAdjointHessian_22");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_SimOpt::applyAdjointHessian_22(a0, a1, a2, a3, a4, a5);
	}
	void applyPreconditioner(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, double & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "applyPreconditioner");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_SimOpt::applyPreconditioner(a0, a1, a2, a3, a4);
	}
	void update(const class ROL::Vector<double> & a0, bool a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_SimOpt::update(a0, a1, a2);
	}
	void update(const class ROL::Vector<double> & a0, enum ROL::UpdateType a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_SimOpt::update(a0, a1, a2);
	}
	void value(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double & a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "value");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_SimOpt::value(a0, a1, a2);
	}
	void applyJacobian(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "applyJacobian");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_SimOpt::applyJacobian(a0, a1, a2, a3);
	}
	void applyAdjointJacobian(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "applyAdjointJacobian");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_SimOpt::applyAdjointJacobian(a0, a1, a2, a3);
	}
	void applyAdjointHessian(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, double & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "applyAdjointHessian");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Constraint_SimOpt::applyAdjointHessian(a0, a1, a2, a3, a4);
	}
	double checkSolve(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const bool a3, std::ostream & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "checkSolve");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return Constraint_SimOpt::checkSolve(a0, a1, a2, a3, a4);
	}
	double checkAdjointConsistencyJacobian_1(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const bool a4, std::ostream & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "checkAdjointConsistencyJacobian_1");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return Constraint_SimOpt::checkAdjointConsistencyJacobian_1(a0, a1, a2, a3, a4, a5);
	}
	double checkAdjointConsistencyJacobian_1(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, const bool a6, std::ostream & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "checkAdjointConsistencyJacobian_1");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return Constraint_SimOpt::checkAdjointConsistencyJacobian_1(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	double checkAdjointConsistencyJacobian_2(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const bool a4, std::ostream & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "checkAdjointConsistencyJacobian_2");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return Constraint_SimOpt::checkAdjointConsistencyJacobian_2(a0, a1, a2, a3, a4, a5);
	}
	double checkAdjointConsistencyJacobian_2(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, const bool a6, std::ostream & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "checkAdjointConsistencyJacobian_2");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return Constraint_SimOpt::checkAdjointConsistencyJacobian_2(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	double checkInverseJacobian_1(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const bool a4, std::ostream & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "checkInverseJacobian_1");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return Constraint_SimOpt::checkInverseJacobian_1(a0, a1, a2, a3, a4, a5);
	}
	double checkInverseAdjointJacobian_1(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const bool a4, std::ostream & a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "checkInverseAdjointJacobian_1");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return Constraint_SimOpt::checkInverseAdjointJacobian_1(a0, a1, a2, a3, a4, a5);
	}
	void applyAdjointJacobian(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, double & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "applyAdjointJacobian");
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
	double checkAdjointConsistencyJacobian(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const bool a3, std::ostream & a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "checkAdjointConsistencyJacobian");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return Constraint::checkAdjointConsistencyJacobian(a0, a1, a2, a3, a4);
	}
	double checkAdjointConsistencyJacobian(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class ROL::Vector<double> & a4, const bool a5, std::ostream & a6) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Constraint_SimOpt<double> *>(this), "checkAdjointConsistencyJacobian");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return Constraint::checkAdjointConsistencyJacobian(a0, a1, a2, a3, a4, a5, a6);
	}
};

// ROL::NonlinearLeastSquaresObjective file:ROL_NonlinearLeastSquaresObjective.hpp line:73
struct PyCallBack_ROL_NonlinearLeastSquaresObjective_double_t : public ROL::NonlinearLeastSquaresObjective<double> {
	using ROL::NonlinearLeastSquaresObjective<double>::NonlinearLeastSquaresObjective;

	void update(const class ROL::Vector<double> & a0, enum ROL::UpdateType a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::NonlinearLeastSquaresObjective<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NonlinearLeastSquaresObjective::update(a0, a1, a2);
	}
	void update(const class ROL::Vector<double> & a0, bool a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::NonlinearLeastSquaresObjective<double> *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NonlinearLeastSquaresObjective::update(a0, a1, a2);
	}
	double value(const class ROL::Vector<double> & a0, double & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::NonlinearLeastSquaresObjective<double> *>(this), "value");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return NonlinearLeastSquaresObjective::value(a0, a1);
	}
	void gradient(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double & a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::NonlinearLeastSquaresObjective<double> *>(this), "gradient");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NonlinearLeastSquaresObjective::gradient(a0, a1, a2);
	}
	void hessVec(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::NonlinearLeastSquaresObjective<double> *>(this), "hessVec");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NonlinearLeastSquaresObjective::hessVec(a0, a1, a2, a3);
	}
	void precond(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::NonlinearLeastSquaresObjective<double> *>(this), "precond");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NonlinearLeastSquaresObjective::precond(a0, a1, a2, a3);
	}
	double dirDeriv(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double & a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::NonlinearLeastSquaresObjective<double> *>(this), "dirDeriv");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::NonlinearLeastSquaresObjective<double> *>(this), "invHessVec");
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
};

// ROL::Objective_FSsolver file:ROL_Objective_FSsolver.hpp line:52
struct PyCallBack_ROL_Objective_FSsolver_double_t : public ROL::Objective_FSsolver<double> {
	using ROL::Objective_FSsolver<double>::Objective_FSsolver;

	double value(const class ROL::Vector<double> & a0, double & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Objective_FSsolver<double> *>(this), "value");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return Objective_FSsolver::value(a0, a1);
	}
	void gradient(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double & a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Objective_FSsolver<double> *>(this), "gradient");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Objective_FSsolver::gradient(a0, a1, a2);
	}
	void hessVec(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Objective_FSsolver<double> *>(this), "hessVec");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Objective_FSsolver::hessVec(a0, a1, a2, a3);
	}
	void update(const class ROL::Vector<double> & a0, enum ROL::UpdateType a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Objective_FSsolver<double> *>(this), "update");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Objective_FSsolver<double> *>(this), "update");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Objective_FSsolver<double> *>(this), "dirDeriv");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Objective_FSsolver<double> *>(this), "invHessVec");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Objective_FSsolver<double> *>(this), "precond");
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

// ROL::BatchManager file:ROL_BatchManager.hpp line:53
struct PyCallBack_ROL_BatchManager_double_t : public ROL::BatchManager<double> {
	using ROL::BatchManager<double>::BatchManager;

	int batchID() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BatchManager<double> *>(this), "batchID");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return BatchManager::batchID();
	}
	int numBatches() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BatchManager<double> *>(this), "numBatches");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return BatchManager::numBatches();
	}
	void sumAll(double * a0, double * a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BatchManager<double> *>(this), "sumAll");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BatchManager::sumAll(a0, a1, a2);
	}
	void sumAll(class ROL::Vector<double> & a0, class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BatchManager<double> *>(this), "sumAll");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BatchManager::sumAll(a0, a1);
	}
	void reduceAll(double * a0, double * a1, int a2, const class ROL::Elementwise::ReductionOp<double> & a3) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BatchManager<double> *>(this), "reduceAll");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BatchManager::reduceAll(a0, a1, a2, a3);
	}
	void gatherAll(const double * a0, const int a1, double * a2, const int a3) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BatchManager<double> *>(this), "gatherAll");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BatchManager::gatherAll(a0, a1, a2, a3);
	}
	void broadcast(double * a0, int a1, int a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BatchManager<double> *>(this), "broadcast");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BatchManager::broadcast(a0, a1, a2);
	}
	void barrier() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BatchManager<double> *>(this), "barrier");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BatchManager::barrier();
	}
};

// ROL::BoundConstraint_SimOpt file: line:41
struct PyCallBack_ROL_BoundConstraint_SimOpt_double_t : public ROL::BoundConstraint_SimOpt<double> {
	using ROL::BoundConstraint_SimOpt<double>::BoundConstraint_SimOpt;

	void project(class ROL::Vector<double> & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BoundConstraint_SimOpt<double> *>(this), "project");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BoundConstraint_SimOpt::project(a0);
	}
	void projectInterior(class ROL::Vector<double> & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BoundConstraint_SimOpt<double> *>(this), "projectInterior");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BoundConstraint_SimOpt::projectInterior(a0);
	}
	void pruneUpperActive(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BoundConstraint_SimOpt<double> *>(this), "pruneUpperActive");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BoundConstraint_SimOpt::pruneUpperActive(a0, a1, a2);
	}
	void pruneUpperActive(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double a3, double a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BoundConstraint_SimOpt<double> *>(this), "pruneUpperActive");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BoundConstraint_SimOpt::pruneUpperActive(a0, a1, a2, a3, a4);
	}
	void pruneLowerActive(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, double a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BoundConstraint_SimOpt<double> *>(this), "pruneLowerActive");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BoundConstraint_SimOpt::pruneLowerActive(a0, a1, a2);
	}
	void pruneLowerActive(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double a3, double a4) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BoundConstraint_SimOpt<double> *>(this), "pruneLowerActive");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BoundConstraint_SimOpt::pruneLowerActive(a0, a1, a2, a3, a4);
	}
	const class Teuchos::RCP<const class ROL::Vector<double> > getLowerBound() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BoundConstraint_SimOpt<double> *>(this), "getLowerBound");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const class Teuchos::RCP<const class ROL::Vector<double> >>::value) {
				static pybind11::detail::override_caster_t<const class Teuchos::RCP<const class ROL::Vector<double> >> caster;
				return pybind11::detail::cast_ref<const class Teuchos::RCP<const class ROL::Vector<double> >>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const class Teuchos::RCP<const class ROL::Vector<double> >>(std::move(o));
		}
		return BoundConstraint_SimOpt::getLowerBound();
	}
	const class Teuchos::RCP<const class ROL::Vector<double> > getUpperBound() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BoundConstraint_SimOpt<double> *>(this), "getUpperBound");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const class Teuchos::RCP<const class ROL::Vector<double> >>::value) {
				static pybind11::detail::override_caster_t<const class Teuchos::RCP<const class ROL::Vector<double> >> caster;
				return pybind11::detail::cast_ref<const class Teuchos::RCP<const class ROL::Vector<double> >>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const class Teuchos::RCP<const class ROL::Vector<double> >>(std::move(o));
		}
		return BoundConstraint_SimOpt::getUpperBound();
	}
	bool isFeasible(const class ROL::Vector<double> & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BoundConstraint_SimOpt<double> *>(this), "isFeasible");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return BoundConstraint_SimOpt::isFeasible(a0);
	}
	void applyInverseScalingFunction(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BoundConstraint_SimOpt<double> *>(this), "applyInverseScalingFunction");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BoundConstraint_SimOpt::applyInverseScalingFunction(a0, a1, a2, a3);
	}
	void applyScalingFunctionJacobian(class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BoundConstraint_SimOpt<double> *>(this), "applyScalingFunctionJacobian");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BoundConstraint_SimOpt::applyScalingFunctionJacobian(a0, a1, a2, a3);
	}
};

// ROL::Algorithm file: line:45
struct PyCallBack_ROL_Algorithm_double_t : public ROL::Algorithm<double> {
	using ROL::Algorithm<double>::Algorithm;

	void run_void(class ROL::Vector<double> & a0, class ROL::Objective<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Algorithm<double> *>(this), "run_void");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Algorithm::run_void(a0, a1);
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

	{ // ROL::Constraint_SimOpt file:ROL_Constraint_SimOpt.hpp line:108
		pybind11::class_<ROL::Constraint_SimOpt<double>, Teuchos::RCP<ROL::Constraint_SimOpt<double>>, PyCallBack_ROL_Constraint_SimOpt_double_t, ROL::Constraint<double>> cl(M("ROL"), "Constraint_SimOpt_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new PyCallBack_ROL_Constraint_SimOpt_double_t(); } ) );
		cl.def(pybind11::init<PyCallBack_ROL_Constraint_SimOpt_double_t const &>());
		cl.def("applyAdjointJacobian", [](ROL::Constraint_SimOpt<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, double & a4) -> void { return o.applyAdjointJacobian(a0, a1, a2, a3, a4); }, "", pybind11::arg("ajv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("dualv"), pybind11::arg("tol"));
		cl.def("update", [](ROL::Constraint_SimOpt<double> &o, const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("u"), pybind11::arg("z"));
		cl.def("update", [](ROL::Constraint_SimOpt<double> &o, const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, bool const & a2) -> void { return o.update(a0, a1, a2); }, "", pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("flag"));
		cl.def("update", (void (ROL::Constraint_SimOpt<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, bool, int)) &ROL::Constraint_SimOpt<double>::update, "C++: ROL::Constraint_SimOpt<double>::update(const class ROL::Vector<double> &, const class ROL::Vector<double> &, bool, int) --> void", pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("flag"), pybind11::arg("iter"));
		cl.def("update", [](ROL::Constraint_SimOpt<double> &o, const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, enum ROL::UpdateType const & a2) -> void { return o.update(a0, a1, a2); }, "", pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("type"));
		cl.def("update", (void (ROL::Constraint_SimOpt<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, enum ROL::UpdateType, int)) &ROL::Constraint_SimOpt<double>::update, "C++: ROL::Constraint_SimOpt<double>::update(const class ROL::Vector<double> &, const class ROL::Vector<double> &, enum ROL::UpdateType, int) --> void", pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("type"), pybind11::arg("iter"));
		cl.def("update_1", [](ROL::Constraint_SimOpt<double> &o, const class ROL::Vector<double> & a0) -> void { return o.update_1(a0); }, "", pybind11::arg("u"));
		cl.def("update_1", [](ROL::Constraint_SimOpt<double> &o, const class ROL::Vector<double> & a0, bool const & a1) -> void { return o.update_1(a0, a1); }, "", pybind11::arg("u"), pybind11::arg("flag"));
		cl.def("update_1", (void (ROL::Constraint_SimOpt<double>::*)(const class ROL::Vector<double> &, bool, int)) &ROL::Constraint_SimOpt<double>::update_1, "C++: ROL::Constraint_SimOpt<double>::update_1(const class ROL::Vector<double> &, bool, int) --> void", pybind11::arg("u"), pybind11::arg("flag"), pybind11::arg("iter"));
		cl.def("update_1", [](ROL::Constraint_SimOpt<double> &o, const class ROL::Vector<double> & a0, enum ROL::UpdateType const & a1) -> void { return o.update_1(a0, a1); }, "", pybind11::arg("u"), pybind11::arg("type"));
		cl.def("update_1", (void (ROL::Constraint_SimOpt<double>::*)(const class ROL::Vector<double> &, enum ROL::UpdateType, int)) &ROL::Constraint_SimOpt<double>::update_1, "C++: ROL::Constraint_SimOpt<double>::update_1(const class ROL::Vector<double> &, enum ROL::UpdateType, int) --> void", pybind11::arg("u"), pybind11::arg("type"), pybind11::arg("iter"));
		cl.def("update_2", [](ROL::Constraint_SimOpt<double> &o, const class ROL::Vector<double> & a0) -> void { return o.update_2(a0); }, "", pybind11::arg("z"));
		cl.def("update_2", [](ROL::Constraint_SimOpt<double> &o, const class ROL::Vector<double> & a0, bool const & a1) -> void { return o.update_2(a0, a1); }, "", pybind11::arg("z"), pybind11::arg("flag"));
		cl.def("update_2", (void (ROL::Constraint_SimOpt<double>::*)(const class ROL::Vector<double> &, bool, int)) &ROL::Constraint_SimOpt<double>::update_2, "C++: ROL::Constraint_SimOpt<double>::update_2(const class ROL::Vector<double> &, bool, int) --> void", pybind11::arg("z"), pybind11::arg("flag"), pybind11::arg("iter"));
		cl.def("update_2", [](ROL::Constraint_SimOpt<double> &o, const class ROL::Vector<double> & a0, enum ROL::UpdateType const & a1) -> void { return o.update_2(a0, a1); }, "", pybind11::arg("z"), pybind11::arg("type"));
		cl.def("update_2", (void (ROL::Constraint_SimOpt<double>::*)(const class ROL::Vector<double> &, enum ROL::UpdateType, int)) &ROL::Constraint_SimOpt<double>::update_2, "C++: ROL::Constraint_SimOpt<double>::update_2(const class ROL::Vector<double> &, enum ROL::UpdateType, int) --> void", pybind11::arg("z"), pybind11::arg("type"), pybind11::arg("iter"));
		cl.def("solve_update", [](ROL::Constraint_SimOpt<double> &o, const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, enum ROL::UpdateType const & a2) -> void { return o.solve_update(a0, a1, a2); }, "", pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("type"));
		cl.def("solve_update", (void (ROL::Constraint_SimOpt<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, enum ROL::UpdateType, int)) &ROL::Constraint_SimOpt<double>::solve_update, "C++: ROL::Constraint_SimOpt<double>::solve_update(const class ROL::Vector<double> &, const class ROL::Vector<double> &, enum ROL::UpdateType, int) --> void", pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("type"), pybind11::arg("iter"));
		cl.def("value", (void (ROL::Constraint_SimOpt<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Constraint_SimOpt<double>::value, "C++: ROL::Constraint_SimOpt<double>::value(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("c"), pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("tol"));
		cl.def("solve", (void (ROL::Constraint_SimOpt<double>::*)(class ROL::Vector<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Constraint_SimOpt<double>::solve, "C++: ROL::Constraint_SimOpt<double>::solve(class ROL::Vector<double> &, class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("c"), pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("tol"));
		cl.def("setSolveParameters", (void (ROL::Constraint_SimOpt<double>::*)(class Teuchos::ParameterList &)) &ROL::Constraint_SimOpt<double>::setSolveParameters, "C++: ROL::Constraint_SimOpt<double>::setSolveParameters(class Teuchos::ParameterList &) --> void", pybind11::arg("parlist"));
		cl.def("applyJacobian_1", (void (ROL::Constraint_SimOpt<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Constraint_SimOpt<double>::applyJacobian_1, "C++: ROL::Constraint_SimOpt<double>::applyJacobian_1(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("jv"), pybind11::arg("v"), pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("tol"));
		cl.def("applyJacobian_2", (void (ROL::Constraint_SimOpt<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Constraint_SimOpt<double>::applyJacobian_2, "C++: ROL::Constraint_SimOpt<double>::applyJacobian_2(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("jv"), pybind11::arg("v"), pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("tol"));
		cl.def("applyInverseJacobian_1", (void (ROL::Constraint_SimOpt<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Constraint_SimOpt<double>::applyInverseJacobian_1, "C++: ROL::Constraint_SimOpt<double>::applyInverseJacobian_1(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("ijv"), pybind11::arg("v"), pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("tol"));
		cl.def("applyAdjointJacobian_1", (void (ROL::Constraint_SimOpt<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Constraint_SimOpt<double>::applyAdjointJacobian_1, "C++: ROL::Constraint_SimOpt<double>::applyAdjointJacobian_1(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("ajv"), pybind11::arg("v"), pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("tol"));
		cl.def("applyAdjointJacobian_1", (void (ROL::Constraint_SimOpt<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Constraint_SimOpt<double>::applyAdjointJacobian_1, "C++: ROL::Constraint_SimOpt<double>::applyAdjointJacobian_1(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("ajv"), pybind11::arg("v"), pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("dualv"), pybind11::arg("tol"));
		cl.def("applyAdjointJacobian_2", (void (ROL::Constraint_SimOpt<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Constraint_SimOpt<double>::applyAdjointJacobian_2, "C++: ROL::Constraint_SimOpt<double>::applyAdjointJacobian_2(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("ajv"), pybind11::arg("v"), pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("tol"));
		cl.def("applyAdjointJacobian_2", (void (ROL::Constraint_SimOpt<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Constraint_SimOpt<double>::applyAdjointJacobian_2, "C++: ROL::Constraint_SimOpt<double>::applyAdjointJacobian_2(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("ajv"), pybind11::arg("v"), pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("dualv"), pybind11::arg("tol"));
		cl.def("applyInverseAdjointJacobian_1", (void (ROL::Constraint_SimOpt<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Constraint_SimOpt<double>::applyInverseAdjointJacobian_1, "C++: ROL::Constraint_SimOpt<double>::applyInverseAdjointJacobian_1(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("iajv"), pybind11::arg("v"), pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("tol"));
		cl.def("applyAdjointHessian_11", (void (ROL::Constraint_SimOpt<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Constraint_SimOpt<double>::applyAdjointHessian_11, "C++: ROL::Constraint_SimOpt<double>::applyAdjointHessian_11(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("ahwv"), pybind11::arg("w"), pybind11::arg("v"), pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("tol"));
		cl.def("applyAdjointHessian_12", (void (ROL::Constraint_SimOpt<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Constraint_SimOpt<double>::applyAdjointHessian_12, "C++: ROL::Constraint_SimOpt<double>::applyAdjointHessian_12(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("ahwv"), pybind11::arg("w"), pybind11::arg("v"), pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("tol"));
		cl.def("applyAdjointHessian_21", (void (ROL::Constraint_SimOpt<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Constraint_SimOpt<double>::applyAdjointHessian_21, "C++: ROL::Constraint_SimOpt<double>::applyAdjointHessian_21(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("ahwv"), pybind11::arg("w"), pybind11::arg("v"), pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("tol"));
		cl.def("applyAdjointHessian_22", (void (ROL::Constraint_SimOpt<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Constraint_SimOpt<double>::applyAdjointHessian_22, "C++: ROL::Constraint_SimOpt<double>::applyAdjointHessian_22(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("ahwv"), pybind11::arg("w"), pybind11::arg("v"), pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("tol"));
		cl.def("applyPreconditioner", (void (ROL::Constraint_SimOpt<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Constraint_SimOpt<double>::applyPreconditioner, "C++: ROL::Constraint_SimOpt<double>::applyPreconditioner(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("pv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("tol"));
		cl.def("update", [](ROL::Constraint_SimOpt<double> &o, const class ROL::Vector<double> & a0) -> void { return o.update(a0); }, "", pybind11::arg("x"));
		cl.def("update", [](ROL::Constraint_SimOpt<double> &o, const class ROL::Vector<double> & a0, bool const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("flag"));
		cl.def("update", (void (ROL::Constraint_SimOpt<double>::*)(const class ROL::Vector<double> &, bool, int)) &ROL::Constraint_SimOpt<double>::update, "C++: ROL::Constraint_SimOpt<double>::update(const class ROL::Vector<double> &, bool, int) --> void", pybind11::arg("x"), pybind11::arg("flag"), pybind11::arg("iter"));
		cl.def("update", [](ROL::Constraint_SimOpt<double> &o, const class ROL::Vector<double> & a0, enum ROL::UpdateType const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("type"));
		cl.def("update", (void (ROL::Constraint_SimOpt<double>::*)(const class ROL::Vector<double> &, enum ROL::UpdateType, int)) &ROL::Constraint_SimOpt<double>::update, "C++: ROL::Constraint_SimOpt<double>::update(const class ROL::Vector<double> &, enum ROL::UpdateType, int) --> void", pybind11::arg("x"), pybind11::arg("type"), pybind11::arg("iter"));
		cl.def("value", (void (ROL::Constraint_SimOpt<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Constraint_SimOpt<double>::value, "C++: ROL::Constraint_SimOpt<double>::value(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("c"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("applyJacobian", (void (ROL::Constraint_SimOpt<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Constraint_SimOpt<double>::applyJacobian, "C++: ROL::Constraint_SimOpt<double>::applyJacobian(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("jv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("applyAdjointJacobian", (void (ROL::Constraint_SimOpt<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Constraint_SimOpt<double>::applyAdjointJacobian, "C++: ROL::Constraint_SimOpt<double>::applyAdjointJacobian(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("ajv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("applyAdjointHessian", (void (ROL::Constraint_SimOpt<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Constraint_SimOpt<double>::applyAdjointHessian, "C++: ROL::Constraint_SimOpt<double>::applyAdjointHessian(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("ahwv"), pybind11::arg("w"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("checkSolve", [](ROL::Constraint_SimOpt<double> &o, const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> double { return o.checkSolve(a0, a1, a2); }, "", pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("c"));
		cl.def("checkSolve", [](ROL::Constraint_SimOpt<double> &o, const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const bool & a3) -> double { return o.checkSolve(a0, a1, a2, a3); }, "", pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("c"), pybind11::arg("printToStream"));
		cl.def("checkSolve", (double (ROL::Constraint_SimOpt<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const bool, std::ostream &)) &ROL::Constraint_SimOpt<double>::checkSolve, "C++: ROL::Constraint_SimOpt<double>::checkSolve(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const bool, std::ostream &) --> double", pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("c"), pybind11::arg("printToStream"), pybind11::arg("outStream"));
		cl.def("checkAdjointConsistencyJacobian_1", [](ROL::Constraint_SimOpt<double> &o, const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3) -> double { return o.checkAdjointConsistencyJacobian_1(a0, a1, a2, a3); }, "", pybind11::arg("w"), pybind11::arg("v"), pybind11::arg("u"), pybind11::arg("z"));
		cl.def("checkAdjointConsistencyJacobian_1", [](ROL::Constraint_SimOpt<double> &o, const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const bool & a4) -> double { return o.checkAdjointConsistencyJacobian_1(a0, a1, a2, a3, a4); }, "", pybind11::arg("w"), pybind11::arg("v"), pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("printToStream"));
		cl.def("checkAdjointConsistencyJacobian_1", (double (ROL::Constraint_SimOpt<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const bool, std::ostream &)) &ROL::Constraint_SimOpt<double>::checkAdjointConsistencyJacobian_1, "C++: ROL::Constraint_SimOpt<double>::checkAdjointConsistencyJacobian_1(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const bool, std::ostream &) --> double", pybind11::arg("w"), pybind11::arg("v"), pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("printToStream"), pybind11::arg("outStream"));
		cl.def("checkAdjointConsistencyJacobian_1", [](ROL::Constraint_SimOpt<double> &o, const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5) -> double { return o.checkAdjointConsistencyJacobian_1(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("w"), pybind11::arg("v"), pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("dualw"), pybind11::arg("dualv"));
		cl.def("checkAdjointConsistencyJacobian_1", [](ROL::Constraint_SimOpt<double> &o, const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, const bool & a6) -> double { return o.checkAdjointConsistencyJacobian_1(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("w"), pybind11::arg("v"), pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("dualw"), pybind11::arg("dualv"), pybind11::arg("printToStream"));
		cl.def("checkAdjointConsistencyJacobian_1", (double (ROL::Constraint_SimOpt<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const bool, std::ostream &)) &ROL::Constraint_SimOpt<double>::checkAdjointConsistencyJacobian_1, "C++: ROL::Constraint_SimOpt<double>::checkAdjointConsistencyJacobian_1(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const bool, std::ostream &) --> double", pybind11::arg("w"), pybind11::arg("v"), pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("dualw"), pybind11::arg("dualv"), pybind11::arg("printToStream"), pybind11::arg("outStream"));
		cl.def("checkAdjointConsistencyJacobian_2", [](ROL::Constraint_SimOpt<double> &o, const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3) -> double { return o.checkAdjointConsistencyJacobian_2(a0, a1, a2, a3); }, "", pybind11::arg("w"), pybind11::arg("v"), pybind11::arg("u"), pybind11::arg("z"));
		cl.def("checkAdjointConsistencyJacobian_2", [](ROL::Constraint_SimOpt<double> &o, const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const bool & a4) -> double { return o.checkAdjointConsistencyJacobian_2(a0, a1, a2, a3, a4); }, "", pybind11::arg("w"), pybind11::arg("v"), pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("printToStream"));
		cl.def("checkAdjointConsistencyJacobian_2", (double (ROL::Constraint_SimOpt<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const bool, std::ostream &)) &ROL::Constraint_SimOpt<double>::checkAdjointConsistencyJacobian_2, "C++: ROL::Constraint_SimOpt<double>::checkAdjointConsistencyJacobian_2(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const bool, std::ostream &) --> double", pybind11::arg("w"), pybind11::arg("v"), pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("printToStream"), pybind11::arg("outStream"));
		cl.def("checkAdjointConsistencyJacobian_2", [](ROL::Constraint_SimOpt<double> &o, const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5) -> double { return o.checkAdjointConsistencyJacobian_2(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("w"), pybind11::arg("v"), pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("dualw"), pybind11::arg("dualv"));
		cl.def("checkAdjointConsistencyJacobian_2", [](ROL::Constraint_SimOpt<double> &o, const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, const bool & a6) -> double { return o.checkAdjointConsistencyJacobian_2(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("w"), pybind11::arg("v"), pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("dualw"), pybind11::arg("dualv"), pybind11::arg("printToStream"));
		cl.def("checkAdjointConsistencyJacobian_2", (double (ROL::Constraint_SimOpt<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const bool, std::ostream &)) &ROL::Constraint_SimOpt<double>::checkAdjointConsistencyJacobian_2, "C++: ROL::Constraint_SimOpt<double>::checkAdjointConsistencyJacobian_2(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const bool, std::ostream &) --> double", pybind11::arg("w"), pybind11::arg("v"), pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("dualw"), pybind11::arg("dualv"), pybind11::arg("printToStream"), pybind11::arg("outStream"));
		cl.def("checkInverseJacobian_1", [](ROL::Constraint_SimOpt<double> &o, const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3) -> double { return o.checkInverseJacobian_1(a0, a1, a2, a3); }, "", pybind11::arg("jv"), pybind11::arg("v"), pybind11::arg("u"), pybind11::arg("z"));
		cl.def("checkInverseJacobian_1", [](ROL::Constraint_SimOpt<double> &o, const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const bool & a4) -> double { return o.checkInverseJacobian_1(a0, a1, a2, a3, a4); }, "", pybind11::arg("jv"), pybind11::arg("v"), pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("printToStream"));
		cl.def("checkInverseJacobian_1", (double (ROL::Constraint_SimOpt<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const bool, std::ostream &)) &ROL::Constraint_SimOpt<double>::checkInverseJacobian_1, "C++: ROL::Constraint_SimOpt<double>::checkInverseJacobian_1(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const bool, std::ostream &) --> double", pybind11::arg("jv"), pybind11::arg("v"), pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("printToStream"), pybind11::arg("outStream"));
		cl.def("checkInverseAdjointJacobian_1", [](ROL::Constraint_SimOpt<double> &o, const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3) -> double { return o.checkInverseAdjointJacobian_1(a0, a1, a2, a3); }, "", pybind11::arg("jv"), pybind11::arg("v"), pybind11::arg("u"), pybind11::arg("z"));
		cl.def("checkInverseAdjointJacobian_1", [](ROL::Constraint_SimOpt<double> &o, const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const bool & a4) -> double { return o.checkInverseAdjointJacobian_1(a0, a1, a2, a3, a4); }, "", pybind11::arg("jv"), pybind11::arg("v"), pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("printToStream"));
		cl.def("checkInverseAdjointJacobian_1", (double (ROL::Constraint_SimOpt<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const bool, std::ostream &)) &ROL::Constraint_SimOpt<double>::checkInverseAdjointJacobian_1, "C++: ROL::Constraint_SimOpt<double>::checkInverseAdjointJacobian_1(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const bool, std::ostream &) --> double", pybind11::arg("jv"), pybind11::arg("v"), pybind11::arg("u"), pybind11::arg("z"), pybind11::arg("printToStream"), pybind11::arg("outStream"));
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
		cl.def("checkAdjointConsistencyJacobian", [](ROL::Constraint<double> &o, const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> double { return o.checkAdjointConsistencyJacobian(a0, a1, a2); }, "", pybind11::arg("w"), pybind11::arg("v"), pybind11::arg("x"));
		cl.def("checkAdjointConsistencyJacobian", [](ROL::Constraint<double> &o, const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const bool & a3) -> double { return o.checkAdjointConsistencyJacobian(a0, a1, a2, a3); }, "", pybind11::arg("w"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("printToStream"));
		cl.def("checkAdjointConsistencyJacobian", (double (ROL::Constraint<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const bool, std::ostream &)) &ROL::Constraint<double>::checkAdjointConsistencyJacobian, "C++: ROL::Constraint<double>::checkAdjointConsistencyJacobian(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const bool, std::ostream &) --> double", pybind11::arg("w"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("printToStream"), pybind11::arg("outStream"));
		cl.def("checkAdjointConsistencyJacobian", [](ROL::Constraint<double> &o, const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class ROL::Vector<double> & a4) -> double { return o.checkAdjointConsistencyJacobian(a0, a1, a2, a3, a4); }, "", pybind11::arg("w"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("dualw"), pybind11::arg("dualv"));
		cl.def("checkAdjointConsistencyJacobian", [](ROL::Constraint<double> &o, const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const class ROL::Vector<double> & a4, const bool & a5) -> double { return o.checkAdjointConsistencyJacobian(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("w"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("dualw"), pybind11::arg("dualv"), pybind11::arg("printToStream"));
		cl.def("checkAdjointConsistencyJacobian", (double (ROL::Constraint<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const bool, std::ostream &)) &ROL::Constraint<double>::checkAdjointConsistencyJacobian, "C++: ROL::Constraint<double>::checkAdjointConsistencyJacobian(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const bool, std::ostream &) --> double", pybind11::arg("w"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("dualw"), pybind11::arg("dualv"), pybind11::arg("printToStream"), pybind11::arg("outStream"));
		cl.def("assign", (class ROL::Constraint<double> & (ROL::Constraint<double>::*)(const class ROL::Constraint<double> &)) &ROL::Constraint<double>::operator=, "C++: ROL::Constraint<double>::operator=(const class ROL::Constraint<double> &) --> class ROL::Constraint<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::NonlinearLeastSquaresObjective file:ROL_NonlinearLeastSquaresObjective.hpp line:73
		pybind11::class_<ROL::NonlinearLeastSquaresObjective<double>, Teuchos::RCP<ROL::NonlinearLeastSquaresObjective<double>>, PyCallBack_ROL_NonlinearLeastSquaresObjective_double_t, ROL::Objective<double>> cl(M("ROL"), "NonlinearLeastSquaresObjective_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](const class Teuchos::RCP<class ROL::Constraint<double> > & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2){ return new ROL::NonlinearLeastSquaresObjective<double>(a0, a1, a2); }, [](const class Teuchos::RCP<class ROL::Constraint<double> > & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2){ return new PyCallBack_ROL_NonlinearLeastSquaresObjective_double_t(a0, a1, a2); } ), "doc");
		cl.def( pybind11::init<const class Teuchos::RCP<class ROL::Constraint<double> > &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const bool>(), pybind11::arg("con"), pybind11::arg("optvec"), pybind11::arg("convec"), pybind11::arg("GNH") );

		cl.def( pybind11::init( [](PyCallBack_ROL_NonlinearLeastSquaresObjective_double_t const &o){ return new PyCallBack_ROL_NonlinearLeastSquaresObjective_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::NonlinearLeastSquaresObjective<double> const &o){ return new ROL::NonlinearLeastSquaresObjective<double>(o); } ) );
		cl.def("update", [](ROL::NonlinearLeastSquaresObjective<double> &o, const class ROL::Vector<double> & a0, enum ROL::UpdateType const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("type"));
		cl.def("update", (void (ROL::NonlinearLeastSquaresObjective<double>::*)(const class ROL::Vector<double> &, enum ROL::UpdateType, int)) &ROL::NonlinearLeastSquaresObjective<double>::update, "C++: ROL::NonlinearLeastSquaresObjective<double>::update(const class ROL::Vector<double> &, enum ROL::UpdateType, int) --> void", pybind11::arg("x"), pybind11::arg("type"), pybind11::arg("iter"));
		cl.def("update", [](ROL::NonlinearLeastSquaresObjective<double> &o, const class ROL::Vector<double> & a0) -> void { return o.update(a0); }, "", pybind11::arg("x"));
		cl.def("update", [](ROL::NonlinearLeastSquaresObjective<double> &o, const class ROL::Vector<double> & a0, bool const & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("x"), pybind11::arg("flag"));
		cl.def("update", (void (ROL::NonlinearLeastSquaresObjective<double>::*)(const class ROL::Vector<double> &, bool, int)) &ROL::NonlinearLeastSquaresObjective<double>::update, "C++: ROL::NonlinearLeastSquaresObjective<double>::update(const class ROL::Vector<double> &, bool, int) --> void", pybind11::arg("x"), pybind11::arg("flag"), pybind11::arg("iter"));
		cl.def("value", (double (ROL::NonlinearLeastSquaresObjective<double>::*)(const class ROL::Vector<double> &, double &)) &ROL::NonlinearLeastSquaresObjective<double>::value, "C++: ROL::NonlinearLeastSquaresObjective<double>::value(const class ROL::Vector<double> &, double &) --> double", pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("gradient", (void (ROL::NonlinearLeastSquaresObjective<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::NonlinearLeastSquaresObjective<double>::gradient, "C++: ROL::NonlinearLeastSquaresObjective<double>::gradient(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("hessVec", (void (ROL::NonlinearLeastSquaresObjective<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::NonlinearLeastSquaresObjective<double>::hessVec, "C++: ROL::NonlinearLeastSquaresObjective<double>::hessVec(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("hv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("tol"));
		cl.def("precond", (void (ROL::NonlinearLeastSquaresObjective<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::NonlinearLeastSquaresObjective<double>::precond, "C++: ROL::NonlinearLeastSquaresObjective<double>::precond(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("Pv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("tol"));
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
	{ // ROL::Objective_FSsolver file:ROL_Objective_FSsolver.hpp line:52
		pybind11::class_<ROL::Objective_FSsolver<double>, Teuchos::RCP<ROL::Objective_FSsolver<double>>, PyCallBack_ROL_Objective_FSsolver_double_t, ROL::Objective<double>> cl(M("ROL"), "Objective_FSsolver_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](PyCallBack_ROL_Objective_FSsolver_double_t const &o){ return new PyCallBack_ROL_Objective_FSsolver_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::Objective_FSsolver<double> const &o){ return new ROL::Objective_FSsolver<double>(o); } ) );
		cl.def( pybind11::init( [](){ return new ROL::Objective_FSsolver<double>(); }, [](){ return new PyCallBack_ROL_Objective_FSsolver_double_t(); } ) );
		cl.def("value", (double (ROL::Objective_FSsolver<double>::*)(const class ROL::Vector<double> &, double &)) &ROL::Objective_FSsolver<double>::value, "C++: ROL::Objective_FSsolver<double>::value(const class ROL::Vector<double> &, double &) --> double", pybind11::arg("u"), pybind11::arg("tol"));
		cl.def("gradient", (void (ROL::Objective_FSsolver<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Objective_FSsolver<double>::gradient, "C++: ROL::Objective_FSsolver<double>::gradient(class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("g"), pybind11::arg("u"), pybind11::arg("tol"));
		cl.def("hessVec", (void (ROL::Objective_FSsolver<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &)) &ROL::Objective_FSsolver<double>::hessVec, "C++: ROL::Objective_FSsolver<double>::hessVec(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double &) --> void", pybind11::arg("hv"), pybind11::arg("v"), pybind11::arg("u"), pybind11::arg("tol"));
		cl.def("assign", (class ROL::Objective_FSsolver<double> & (ROL::Objective_FSsolver<double>::*)(const class ROL::Objective_FSsolver<double> &)) &ROL::Objective_FSsolver<double>::operator=, "C++: ROL::Objective_FSsolver<double>::operator=(const class ROL::Objective_FSsolver<double> &) --> class ROL::Objective_FSsolver<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
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
	{ // ROL::BatchManager file:ROL_BatchManager.hpp line:53
		pybind11::class_<ROL::BatchManager<double>, Teuchos::RCP<ROL::BatchManager<double>>, PyCallBack_ROL_BatchManager_double_t> cl(M("ROL"), "BatchManager_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new ROL::BatchManager<double>(); }, [](){ return new PyCallBack_ROL_BatchManager_double_t(); } ) );
		cl.def("batchID", (int (ROL::BatchManager<double>::*)()) &ROL::BatchManager<double>::batchID, "C++: ROL::BatchManager<double>::batchID() --> int");
		cl.def("numBatches", (int (ROL::BatchManager<double>::*)()) &ROL::BatchManager<double>::numBatches, "C++: ROL::BatchManager<double>::numBatches() --> int");
		cl.def("sumAll", (void (ROL::BatchManager<double>::*)(double *, double *, int)) &ROL::BatchManager<double>::sumAll, "C++: ROL::BatchManager<double>::sumAll(double *, double *, int) --> void", pybind11::arg("input"), pybind11::arg("output"), pybind11::arg("dim"));
		cl.def("sumAll", (void (ROL::BatchManager<double>::*)(class ROL::Vector<double> &, class ROL::Vector<double> &)) &ROL::BatchManager<double>::sumAll, "C++: ROL::BatchManager<double>::sumAll(class ROL::Vector<double> &, class ROL::Vector<double> &) --> void", pybind11::arg("input"), pybind11::arg("output"));
		cl.def("reduceAll", (void (ROL::BatchManager<double>::*)(double *, double *, int, const class ROL::Elementwise::ReductionOp<double> &)) &ROL::BatchManager<double>::reduceAll, "C++: ROL::BatchManager<double>::reduceAll(double *, double *, int, const class ROL::Elementwise::ReductionOp<double> &) --> void", pybind11::arg("input"), pybind11::arg("output"), pybind11::arg("dim"), pybind11::arg("r"));
		cl.def("gatherAll", (void (ROL::BatchManager<double>::*)(const double *, const int, double *, const int) const) &ROL::BatchManager<double>::gatherAll, "C++: ROL::BatchManager<double>::gatherAll(const double *, const int, double *, const int) const --> void", pybind11::arg("send"), pybind11::arg("ssize"), pybind11::arg("receive"), pybind11::arg("rsize"));
		cl.def("broadcast", (void (ROL::BatchManager<double>::*)(double *, int, int)) &ROL::BatchManager<double>::broadcast, "C++: ROL::BatchManager<double>::broadcast(double *, int, int) --> void", pybind11::arg("input"), pybind11::arg("cnt"), pybind11::arg("root"));
		cl.def("barrier", (void (ROL::BatchManager<double>::*)()) &ROL::BatchManager<double>::barrier, "C++: ROL::BatchManager<double>::barrier() --> void");
		cl.def("assign", (class ROL::BatchManager<double> & (ROL::BatchManager<double>::*)(const class ROL::BatchManager<double> &)) &ROL::BatchManager<double>::operator=, "C++: ROL::BatchManager<double>::operator=(const class ROL::BatchManager<double> &) --> class ROL::BatchManager<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::BoundConstraint_SimOpt file: line:41
		pybind11::class_<ROL::BoundConstraint_SimOpt<double>, Teuchos::RCP<ROL::BoundConstraint_SimOpt<double>>, PyCallBack_ROL_BoundConstraint_SimOpt_double_t, ROL::BoundConstraint<double>> cl(M("ROL"), "BoundConstraint_SimOpt_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<const class Teuchos::RCP<class ROL::BoundConstraint<double> > &, const class Teuchos::RCP<class ROL::BoundConstraint<double> > &>(), pybind11::arg("bnd1"), pybind11::arg("bnd2") );

		cl.def( pybind11::init( [](PyCallBack_ROL_BoundConstraint_SimOpt_double_t const &o){ return new PyCallBack_ROL_BoundConstraint_SimOpt_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::BoundConstraint_SimOpt<double> const &o){ return new ROL::BoundConstraint_SimOpt<double>(o); } ) );
		cl.def("project", (void (ROL::BoundConstraint_SimOpt<double>::*)(class ROL::Vector<double> &)) &ROL::BoundConstraint_SimOpt<double>::project, "Project optimization variables onto the bounds.\n\n      This function implements the projection of \n onto the bounds, i.e.,\n      \n\n\n\n       \n is the optimization variable.\n\nC++: ROL::BoundConstraint_SimOpt<double>::project(class ROL::Vector<double> &) --> void", pybind11::arg("x"));
		cl.def("projectInterior", (void (ROL::BoundConstraint_SimOpt<double>::*)(class ROL::Vector<double> &)) &ROL::BoundConstraint_SimOpt<double>::projectInterior, "Project optimization variables into the interior of the feasible set.\n\n      This function implements the projection of \n into the interior of the\n      feasible set, i.e.,\n      \n\n\n\n\n       \n is the optimization variable.\n\nC++: ROL::BoundConstraint_SimOpt<double>::projectInterior(class ROL::Vector<double> &) --> void", pybind11::arg("x"));
		cl.def("pruneUpperActive", [](ROL::BoundConstraint_SimOpt<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) -> void { return o.pruneUpperActive(a0, a1); }, "", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("pruneUpperActive", (void (ROL::BoundConstraint_SimOpt<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double)) &ROL::BoundConstraint_SimOpt<double>::pruneUpperActive, "Set variables to zero if they correspond to the upper \n-active set.\n\n      This function sets \n if \n.  Here,\n      the upper \n\n-active set is defined as\n      \n\n\n\n      \n   is the variable to be pruned.\n      \n\n   is the current optimization variable.\n      \n\n is the active-set tolerance \n.\n\nC++: ROL::BoundConstraint_SimOpt<double>::pruneUpperActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, double) --> void", pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("eps"));
		cl.def("pruneUpperActive", [](ROL::BoundConstraint_SimOpt<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> void { return o.pruneUpperActive(a0, a1, a2); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"));
		cl.def("pruneUpperActive", [](ROL::BoundConstraint_SimOpt<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double const & a3) -> void { return o.pruneUpperActive(a0, a1, a2, a3); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"));
		cl.def("pruneUpperActive", (void (ROL::BoundConstraint_SimOpt<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double)) &ROL::BoundConstraint_SimOpt<double>::pruneUpperActive, "Set variables to zero if they correspond to the upper \n-binding set.\n\n      This function sets \n if \n.  Here,\n      the upper \n\n-binding set is defined as\n      \n\n\n\n\n      \n   is the variable to be pruned.\n      \n\n   is the current optimization variable.\n      \n\n   is the negative search direction.\n      \n\n is the active-set tolerance \n.\n\nC++: ROL::BoundConstraint_SimOpt<double>::pruneUpperActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double) --> void", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"), pybind11::arg("geps"));
		cl.def("pruneLowerActive", [](ROL::BoundConstraint_SimOpt<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) -> void { return o.pruneLowerActive(a0, a1); }, "", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("pruneLowerActive", (void (ROL::BoundConstraint_SimOpt<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double)) &ROL::BoundConstraint_SimOpt<double>::pruneLowerActive, "Set variables to zero if they correspond to the lower \n-active set.\n\n      This function sets \n if \n.  Here,\n      the lower \n\n-active set is defined as\n      \n\n\n\n      \n   is the variable to be pruned.\n      \n\n   is the current optimization variable.\n      \n\n is the active-set tolerance \n.\n\nC++: ROL::BoundConstraint_SimOpt<double>::pruneLowerActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, double) --> void", pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("eps"));
		cl.def("pruneLowerActive", [](ROL::BoundConstraint_SimOpt<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> void { return o.pruneLowerActive(a0, a1, a2); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"));
		cl.def("pruneLowerActive", [](ROL::BoundConstraint_SimOpt<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double const & a3) -> void { return o.pruneLowerActive(a0, a1, a2, a3); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"));
		cl.def("pruneLowerActive", (void (ROL::BoundConstraint_SimOpt<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double)) &ROL::BoundConstraint_SimOpt<double>::pruneLowerActive, "Set variables to zero if they correspond to the lower \n-binding set.\n\n      This function sets \n if \n.  Here,\n      the lower \n\n-binding set is defined as\n      \n\n\n\n\n      \n   is the variable to be pruned.\n      \n\n   is the current optimization variable.\n      \n\n   is the negative search direction.\n      \n\n is the active-set tolerance \n.\n\nC++: ROL::BoundConstraint_SimOpt<double>::pruneLowerActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double) --> void", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"), pybind11::arg("geps"));
		cl.def("getLowerBound", (const class Teuchos::RCP<const class ROL::Vector<double> > (ROL::BoundConstraint_SimOpt<double>::*)() const) &ROL::BoundConstraint_SimOpt<double>::getLowerBound, "C++: ROL::BoundConstraint_SimOpt<double>::getLowerBound() const --> const class Teuchos::RCP<const class ROL::Vector<double> >");
		cl.def("getUpperBound", (const class Teuchos::RCP<const class ROL::Vector<double> > (ROL::BoundConstraint_SimOpt<double>::*)() const) &ROL::BoundConstraint_SimOpt<double>::getUpperBound, "C++: ROL::BoundConstraint_SimOpt<double>::getUpperBound() const --> const class Teuchos::RCP<const class ROL::Vector<double> >");
		cl.def("pruneActive", [](ROL::BoundConstraint_SimOpt<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) -> void { return o.pruneActive(a0, a1); }, "", pybind11::arg("v"), pybind11::arg("x"));
		cl.def("pruneActive", (void (ROL::BoundConstraint_SimOpt<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, double)) &ROL::BoundConstraint_SimOpt<double>::pruneActive, "Set variables to zero if they correspond to the \n-active set.\n\n      This function sets \n if \n.  Here,\n      the \n\n-active set is defined as\n      \n\n\n\n      \n   is the variable to be pruned.\n      \n\n   is the current optimization variable.\n      \n\n is the active-set tolerance \n.\n\nC++: ROL::BoundConstraint_SimOpt<double>::pruneActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, double) --> void", pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("eps"));
		cl.def("pruneActive", [](ROL::BoundConstraint_SimOpt<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2) -> void { return o.pruneActive(a0, a1, a2); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"));
		cl.def("pruneActive", [](ROL::BoundConstraint_SimOpt<double> &o, class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, double const & a3) -> void { return o.pruneActive(a0, a1, a2, a3); }, "", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"));
		cl.def("pruneActive", (void (ROL::BoundConstraint_SimOpt<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double)) &ROL::BoundConstraint_SimOpt<double>::pruneActive, "Set variables to zero if they correspond to the \n-binding set.\n\n      This function sets \n if \n.  Here,\n      the \n\n-binding set is defined as\n      \n\n\n\n      \n   is the variable to be pruned.\n      \n\n   is the current optimization variable.\n      \n\n   is the negative search direction.\n      \n\n is the active-set tolerance \n.\n\nC++: ROL::BoundConstraint_SimOpt<double>::pruneActive(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, double, double) --> void", pybind11::arg("v"), pybind11::arg("g"), pybind11::arg("x"), pybind11::arg("xeps"), pybind11::arg("geps"));
		cl.def("isFeasible", (bool (ROL::BoundConstraint_SimOpt<double>::*)(const class ROL::Vector<double> &)) &ROL::BoundConstraint_SimOpt<double>::isFeasible, "Check if the vector, v, is feasible.\n\n      This function returns true if \n.\n      \n\n   is the vector to be checked.\n\nC++: ROL::BoundConstraint_SimOpt<double>::isFeasible(const class ROL::Vector<double> &) --> bool", pybind11::arg("v"));
		cl.def("applyInverseScalingFunction", (void (ROL::BoundConstraint_SimOpt<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) const) &ROL::BoundConstraint_SimOpt<double>::applyInverseScalingFunction, "Apply inverse scaling function.\n\n      This function applies the inverse scaling function \n to\n      a vector \n\n, i.e., the output is \n.\n      The scaling function must satisfy:\n      (i) \n\n if \n and \n;\n      (ii) \n\n if \n and \n; and\n      (iii) \n\n otherwise.\n      \n\n   is the inverse scaling function applied to v.\n      \n\n   is the vector being scaled.\n      \n\n   is the primal vector at which the scaling function is evaluated.\n      \n\n   is the dual vector at which the scaling function is evaluated.\n\nC++: ROL::BoundConstraint_SimOpt<double>::applyInverseScalingFunction(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) const --> void", pybind11::arg("dv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("g"));
		cl.def("applyScalingFunctionJacobian", (void (ROL::BoundConstraint_SimOpt<double>::*)(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) const) &ROL::BoundConstraint_SimOpt<double>::applyScalingFunctionJacobian, "Apply scaling function Jacobian.\n\n      This function applies the Jacobian of the scaling function \n to\n      a vector \n\n.  The output is \n.  The\n      scaling function must satisfy:\n      (i) \n\n if \n and \n;\n      (ii) \n\n if \n and \n; and\n      (iii) \n\n otherwise.\n      \n\n   is the scaling function Jacobian applied to v.\n      \n\n   is the vector being scaled.\n      \n\n   is the primal vector at which the scaling function is evaluated.\n      \n\n   is the dual vector at which the scaling function is evaluated.\n\nC++: ROL::BoundConstraint_SimOpt<double>::applyScalingFunctionJacobian(class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &) const --> void", pybind11::arg("dv"), pybind11::arg("v"), pybind11::arg("x"), pybind11::arg("g"));
		cl.def("assign", (class ROL::BoundConstraint_SimOpt<double> & (ROL::BoundConstraint_SimOpt<double>::*)(const class ROL::BoundConstraint_SimOpt<double> &)) &ROL::BoundConstraint_SimOpt<double>::operator=, "C++: ROL::BoundConstraint_SimOpt<double>::operator=(const class ROL::BoundConstraint_SimOpt<double> &) --> class ROL::BoundConstraint_SimOpt<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
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
	{ // ROL::Solver file: line:44
		pybind11::class_<ROL::Solver<double>, Teuchos::RCP<ROL::Solver<double>>> cl(M("ROL"), "Solver_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<const class Teuchos::RCP<class ROL::Problem<double> > &, class Teuchos::ParameterList &>(), pybind11::arg("opt"), pybind11::arg("parlist") );

		cl.def( pybind11::init( [](ROL::Solver<double> const &o){ return new ROL::Solver<double>(o); } ) );
		cl.def("solve", [](ROL::Solver<double> &o) -> int { return o.solve(); }, "");
		cl.def("solve", [](ROL::Solver<double> &o, const class Teuchos::RCP<class ROL::StatusTest<double> > & a0) -> int { return o.solve(a0); }, "", pybind11::arg("status"));
		cl.def("solve", (int (ROL::Solver<double>::*)(const class Teuchos::RCP<class ROL::StatusTest<double> > &, bool)) &ROL::Solver<double>::solve, "Solve optimization problem with no iteration output.\n\n      \n          is a user-defined StatusTest\n      \n\n   if true, the user-defined StatusTest will be combined with the default StatusTest\n\n      ---\n\nC++: ROL::Solver<double>::solve(const class Teuchos::RCP<class ROL::StatusTest<double> > &, bool) --> int", pybind11::arg("status"), pybind11::arg("combineStatus"));
		cl.def("solve", [](ROL::Solver<double> &o, std::ostream & a0) -> int { return o.solve(a0); }, "", pybind11::arg("outStream"));
		cl.def("solve", [](ROL::Solver<double> &o, std::ostream & a0, const class Teuchos::RCP<class ROL::StatusTest<double> > & a1) -> int { return o.solve(a0, a1); }, "", pybind11::arg("outStream"), pybind11::arg("status"));
		cl.def("solve", (int (ROL::Solver<double>::*)(std::ostream &, const class Teuchos::RCP<class ROL::StatusTest<double> > &, bool)) &ROL::Solver<double>::solve, "Solve optimization problem.\n\n      \n       is the output stream to collect iteration history\n      \n\n          is a user-defined StatusTest\n      \n\n   if true, the user-defined StatusTest will be combined with the default StatusTest\n\n      ---\n\nC++: ROL::Solver<double>::solve(std::ostream &, const class Teuchos::RCP<class ROL::StatusTest<double> > &, bool) --> int", pybind11::arg("outStream"), pybind11::arg("status"), pybind11::arg("combineStatus"));
		cl.def("getAlgorithmState", (class Teuchos::RCP<const struct ROL::AlgorithmState<double> > (ROL::Solver<double>::*)() const) &ROL::Solver<double>::getAlgorithmState, "C++: ROL::Solver<double>::getAlgorithmState() const --> class Teuchos::RCP<const struct ROL::AlgorithmState<double> >");
		cl.def("reset", (void (ROL::Solver<double>::*)()) &ROL::Solver<double>::reset, "Reset both Algorithm and Step.\n\n      This function will reset the AlgorithmState and reinitialize the\n      Step.  This function does not permit changing the Step specified\n      upon construction.  To change the Step, reinitialize the\n      OptimizationSolver.\n\n      ---\n\nC++: ROL::Solver<double>::reset() --> void");
	}
	{ // ROL::Algorithm file: line:45
		pybind11::class_<ROL::Algorithm<double>, Teuchos::RCP<ROL::Algorithm<double>>, PyCallBack_ROL_Algorithm_double_t> cl(M("ROL"), "Algorithm_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](const class Teuchos::RCP<class ROL::Step<double> > & a0, const class Teuchos::RCP<class ROL::StatusTest<double> > & a1){ return new ROL::Algorithm<double>(a0, a1); }, [](const class Teuchos::RCP<class ROL::Step<double> > & a0, const class Teuchos::RCP<class ROL::StatusTest<double> > & a1){ return new PyCallBack_ROL_Algorithm_double_t(a0, a1); } ), "doc");
		cl.def( pybind11::init<const class Teuchos::RCP<class ROL::Step<double> > &, const class Teuchos::RCP<class ROL::StatusTest<double> > &, bool>(), pybind11::arg("step"), pybind11::arg("status"), pybind11::arg("printHeader") );

		cl.def( pybind11::init( [](const class Teuchos::RCP<class ROL::Step<double> > & a0, const class Teuchos::RCP<class ROL::StatusTest<double> > & a1, const class Teuchos::RCP<struct ROL::AlgorithmState<double> > & a2){ return new ROL::Algorithm<double>(a0, a1, a2); }, [](const class Teuchos::RCP<class ROL::Step<double> > & a0, const class Teuchos::RCP<class ROL::StatusTest<double> > & a1, const class Teuchos::RCP<struct ROL::AlgorithmState<double> > & a2){ return new PyCallBack_ROL_Algorithm_double_t(a0, a1, a2); } ), "doc");
		cl.def( pybind11::init<const class Teuchos::RCP<class ROL::Step<double> > &, const class Teuchos::RCP<class ROL::StatusTest<double> > &, const class Teuchos::RCP<struct ROL::AlgorithmState<double> > &, bool>(), pybind11::arg("step"), pybind11::arg("status"), pybind11::arg("state"), pybind11::arg("printHeader") );

		cl.def( pybind11::init( [](PyCallBack_ROL_Algorithm_double_t const &o){ return new PyCallBack_ROL_Algorithm_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::Algorithm<double> const &o){ return new ROL::Algorithm<double>(o); } ) );
		cl.def("run_void", (void (ROL::Algorithm<double>::*)(class ROL::Vector<double> &, class ROL::Objective<double> &)) &ROL::Algorithm<double>::run_void, "Run algorithm on unconstrained problems (Type-U).\n             This is the primary Type-U interface.\n\nC++: ROL::Algorithm<double>::run_void(class ROL::Vector<double> &, class ROL::Objective<double> &) --> void", pybind11::arg("x"), pybind11::arg("obj"));
		cl.def("getIterHeader", (std::string (ROL::Algorithm<double>::*)()) &ROL::Algorithm<double>::getIterHeader, "C++: ROL::Algorithm<double>::getIterHeader() --> std::string");
		cl.def("getIterInfo", [](ROL::Algorithm<double> &o) -> std::string { return o.getIterInfo(); }, "");
		cl.def("getIterInfo", (std::string (ROL::Algorithm<double>::*)(bool)) &ROL::Algorithm<double>::getIterInfo, "C++: ROL::Algorithm<double>::getIterInfo(bool) --> std::string", pybind11::arg("withHeader"));
		cl.def("getState", (class Teuchos::RCP<const struct ROL::AlgorithmState<double> > (ROL::Algorithm<double>::*)() const) &ROL::Algorithm<double>::getState, "C++: ROL::Algorithm<double>::getState() const --> class Teuchos::RCP<const struct ROL::AlgorithmState<double> >");
		cl.def("reset", (void (ROL::Algorithm<double>::*)()) &ROL::Algorithm<double>::reset, "C++: ROL::Algorithm<double>::reset() --> void");
		cl.def("assign", (class ROL::Algorithm<double> & (ROL::Algorithm<double>::*)(const class ROL::Algorithm<double> &)) &ROL::Algorithm<double>::operator=, "C++: ROL::Algorithm<double>::operator=(const class ROL::Algorithm<double> &) --> class ROL::Algorithm<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
}
