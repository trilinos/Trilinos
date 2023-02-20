#include <ROL_BackTracking_U.hpp>
#include <ROL_BisectionScalarMinimization.hpp>
#include <ROL_Bracketing.hpp>
#include <ROL_BrentsScalarMinimization.hpp>
#include <ROL_BundleStatusTest.hpp>
#include <ROL_Bundle_U.hpp>
#include <ROL_Bundle_U_AS.hpp>
#include <ROL_Bundle_U_TT.hpp>
#include <ROL_CubicInterp_U.hpp>
#include <ROL_DescentDirection_U.hpp>
#include <ROL_Elementwise_Function.hpp>
#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_GoldenSectionScalarMinimization.hpp>
#include <ROL_IterationScaling_U.hpp>
#include <ROL_LineSearch_U.hpp>
#include <ROL_LineSearch_U_Factory.hpp>
#include <ROL_LineSearch_U_Types.hpp>
#include <ROL_Objective.hpp>
#include <ROL_PartitionedVector.hpp>
#include <ROL_PathBasedTargetLevel_U.hpp>
#include <ROL_ScalarFunction.hpp>
#include <ROL_ScalarMinimization.hpp>
#include <ROL_ScalarMinimizationLineSearch_U.hpp>
#include <ROL_ScalarMinimizationStatusTest.hpp>
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
#include "PyROL_Smart_Holder.hpp"
#include <string>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>


#ifndef BINDER_PYBIND11_TYPE_CASTER
	#define BINDER_PYBIND11_TYPE_CASTER
	PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)
	PYBIND11_DECLARE_HOLDER_TYPE(T, T*)
	PYBIND11_MAKE_OPAQUE(std::shared_ptr<void>)
#endif

// ROL::BundleStatusTest file:ROL_BundleStatusTest.hpp line:53
struct PyCallBack_ROL_BundleStatusTest_double_t : public ROL::BundleStatusTest<double> {
	using ROL::BundleStatusTest<double>::BundleStatusTest;

	bool check(struct ROL::AlgorithmState<double> & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BundleStatusTest<double> *>(this), "check");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return BundleStatusTest::check(a0);
	}
};

// ROL::Bundle_U_AS file:ROL_Bundle_U_AS.hpp line:56
struct PyCallBack_ROL_Bundle_U_AS_double_t : public ROL::Bundle_U_AS<double> {
	using ROL::Bundle_U_AS<double>::Bundle_U_AS;

	void initialize(const class ROL::Vector<double> & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Bundle_U_AS<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Bundle_U_AS::initialize(a0);
	}
	unsigned int solveDual(const double a0, const unsigned int a1, const double a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Bundle_U_AS<double> *>(this), "solveDual");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<unsigned int>::value) {
				static pybind11::detail::override_caster_t<unsigned int> caster;
				return pybind11::detail::cast_ref<unsigned int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<unsigned int>(std::move(o));
		}
		return Bundle_U_AS::solveDual(a0, a1, a2);
	}
};

// ROL::Bundle_U_TT file:ROL_Bundle_U_TT.hpp line:61
struct PyCallBack_ROL_Bundle_U_TT_double_t : public ROL::Bundle_U_TT<double> {
	using ROL::Bundle_U_TT<double>::Bundle_U_TT;

	unsigned int solveDual(const double a0, const unsigned int a1, const double a2) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Bundle_U_TT<double> *>(this), "solveDual");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<unsigned int>::value) {
				static pybind11::detail::override_caster_t<unsigned int> caster;
				return pybind11::detail::cast_ref<unsigned int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<unsigned int>(std::move(o));
		}
		return Bundle_U_TT::solveDual(a0, a1, a2);
	}
	void initialize(const class ROL::Vector<double> & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Bundle_U_TT<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Bundle_U::initialize(a0);
	}
};

// ROL::IterationScaling_U file:ROL_IterationScaling_U.hpp line:56
struct PyCallBack_ROL_IterationScaling_U_double_t : public ROL::IterationScaling_U<double> {
	using ROL::IterationScaling_U<double>::IterationScaling_U;

	void initialize(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::IterationScaling_U<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return IterationScaling_U::initialize(a0, a1);
	}
	void run(double & a0, double & a1, int & a2, int & a3, const double & a4, const class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Objective<double> & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::IterationScaling_U<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return IterationScaling_U::run(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	bool status(const enum ROL::ELineSearchU a0, int & a1, int & a2, const double a3, const double a4, const double a5, const double a6, const class ROL::Vector<double> & a7, const class ROL::Vector<double> & a8, class ROL::Objective<double> & a9) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::IterationScaling_U<double> *>(this), "status");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LineSearch_U::status(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9);
	}
	double getInitialAlpha(int & a0, int & a1, const double a2, const double a3, const class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Objective<double> & a6) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::IterationScaling_U<double> *>(this), "getInitialAlpha");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return LineSearch_U::getInitialAlpha(a0, a1, a2, a3, a4, a5, a6);
	}
};

// ROL::PathBasedTargetLevel_U file:ROL_PathBasedTargetLevel_U.hpp line:56
struct PyCallBack_ROL_PathBasedTargetLevel_U_double_t : public ROL::PathBasedTargetLevel_U<double> {
	using ROL::PathBasedTargetLevel_U<double>::PathBasedTargetLevel_U;

	void initialize(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::PathBasedTargetLevel_U<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PathBasedTargetLevel_U::initialize(a0, a1);
	}
	void run(double & a0, double & a1, int & a2, int & a3, const double & a4, const class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Objective<double> & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::PathBasedTargetLevel_U<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PathBasedTargetLevel_U::run(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	bool status(const enum ROL::ELineSearchU a0, int & a1, int & a2, const double a3, const double a4, const double a5, const double a6, const class ROL::Vector<double> & a7, const class ROL::Vector<double> & a8, class ROL::Objective<double> & a9) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::PathBasedTargetLevel_U<double> *>(this), "status");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LineSearch_U::status(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9);
	}
	double getInitialAlpha(int & a0, int & a1, const double a2, const double a3, const class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Objective<double> & a6) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::PathBasedTargetLevel_U<double> *>(this), "getInitialAlpha");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return LineSearch_U::getInitialAlpha(a0, a1, a2, a3, a4, a5, a6);
	}
};

// ROL::BackTracking_U file:ROL_BackTracking_U.hpp line:56
struct PyCallBack_ROL_BackTracking_U_double_t : public ROL::BackTracking_U<double> {
	using ROL::BackTracking_U<double>::BackTracking_U;

	void initialize(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BackTracking_U<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BackTracking_U::initialize(a0, a1);
	}
	void run(double & a0, double & a1, int & a2, int & a3, const double & a4, const class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Objective<double> & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BackTracking_U<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BackTracking_U::run(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	bool status(const enum ROL::ELineSearchU a0, int & a1, int & a2, const double a3, const double a4, const double a5, const double a6, const class ROL::Vector<double> & a7, const class ROL::Vector<double> & a8, class ROL::Objective<double> & a9) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BackTracking_U<double> *>(this), "status");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LineSearch_U::status(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9);
	}
	double getInitialAlpha(int & a0, int & a1, const double a2, const double a3, const class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Objective<double> & a6) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BackTracking_U<double> *>(this), "getInitialAlpha");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return LineSearch_U::getInitialAlpha(a0, a1, a2, a3, a4, a5, a6);
	}
};

// ROL::CubicInterp_U file:ROL_CubicInterp_U.hpp line:56
struct PyCallBack_ROL_CubicInterp_U_double_t : public ROL::CubicInterp_U<double> {
	using ROL::CubicInterp_U<double>::CubicInterp_U;

	void initialize(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::CubicInterp_U<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return CubicInterp_U::initialize(a0, a1);
	}
	void run(double & a0, double & a1, int & a2, int & a3, const double & a4, const class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Objective<double> & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::CubicInterp_U<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return CubicInterp_U::run(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	bool status(const enum ROL::ELineSearchU a0, int & a1, int & a2, const double a3, const double a4, const double a5, const double a6, const class ROL::Vector<double> & a7, const class ROL::Vector<double> & a8, class ROL::Objective<double> & a9) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::CubicInterp_U<double> *>(this), "status");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LineSearch_U::status(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9);
	}
	double getInitialAlpha(int & a0, int & a1, const double a2, const double a3, const class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Objective<double> & a6) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::CubicInterp_U<double> *>(this), "getInitialAlpha");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return LineSearch_U::getInitialAlpha(a0, a1, a2, a3, a4, a5, a6);
	}
};

// ROL::ScalarMinimizationStatusTest file:ROL_ScalarMinimizationStatusTest.hpp line:54
struct PyCallBack_ROL_ScalarMinimizationStatusTest_double_t : public ROL::ScalarMinimizationStatusTest<double> {
	using ROL::ScalarMinimizationStatusTest<double>::ScalarMinimizationStatusTest;

	bool check(double & a0, double & a1, double & a2, int & a3, int & a4, const bool a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::ScalarMinimizationStatusTest<double> *>(this), "check");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ScalarMinimizationStatusTest::check(a0, a1, a2, a3, a4, a5);
	}
};

// ROL::ScalarMinimization file:ROL_ScalarMinimization.hpp line:59
struct PyCallBack_ROL_ScalarMinimization_double_t : public ROL::ScalarMinimization<double> {
	using ROL::ScalarMinimization<double>::ScalarMinimization;

	void run(double & a0, double & a1, int & a2, int & a3, class ROL::ScalarFunction<double> & a4, const double a5, const double a6, class ROL::ScalarMinimizationStatusTest<double> & a7) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::ScalarMinimization<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"ScalarMinimization::run\"");
	}
};

// ROL::BrentsScalarMinimization file:ROL_BrentsScalarMinimization.hpp line:60
struct PyCallBack_ROL_BrentsScalarMinimization_double_t : public ROL::BrentsScalarMinimization<double> {
	using ROL::BrentsScalarMinimization<double>::BrentsScalarMinimization;

	void run(double & a0, double & a1, int & a2, int & a3, class ROL::ScalarFunction<double> & a4, const double a5, const double a6, class ROL::ScalarMinimizationStatusTest<double> & a7) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BrentsScalarMinimization<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BrentsScalarMinimization::run(a0, a1, a2, a3, a4, a5, a6, a7);
	}
};

// ROL::BisectionScalarMinimization file:ROL_BisectionScalarMinimization.hpp line:60
struct PyCallBack_ROL_BisectionScalarMinimization_double_t : public ROL::BisectionScalarMinimization<double> {
	using ROL::BisectionScalarMinimization<double>::BisectionScalarMinimization;

	void run(double & a0, double & a1, int & a2, int & a3, class ROL::ScalarFunction<double> & a4, const double a5, const double a6, class ROL::ScalarMinimizationStatusTest<double> & a7) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::BisectionScalarMinimization<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return BisectionScalarMinimization::run(a0, a1, a2, a3, a4, a5, a6, a7);
	}
};

// ROL::GoldenSectionScalarMinimization file:ROL_GoldenSectionScalarMinimization.hpp line:59
struct PyCallBack_ROL_GoldenSectionScalarMinimization_double_t : public ROL::GoldenSectionScalarMinimization<double> {
	using ROL::GoldenSectionScalarMinimization<double>::GoldenSectionScalarMinimization;

	void run(double & a0, double & a1, int & a2, int & a3, class ROL::ScalarFunction<double> & a4, const double a5, const double a6, class ROL::ScalarMinimizationStatusTest<double> & a7) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::GoldenSectionScalarMinimization<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return GoldenSectionScalarMinimization::run(a0, a1, a2, a3, a4, a5, a6, a7);
	}
};

// ROL::Bracketing file:ROL_Bracketing.hpp line:57
struct PyCallBack_ROL_Bracketing_double_t : public ROL::Bracketing<double> {
	using ROL::Bracketing<double>::Bracketing;

	void run(double & a0, double & a1, double & a2, double & a3, double & a4, double & a5, int & a6, int & a7, class ROL::ScalarFunction<double> & a8, class ROL::ScalarMinimizationStatusTest<double> & a9) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::Bracketing<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return Bracketing::run(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9);
	}
};

// ROL::ScalarMinimizationLineSearch_U file:ROL_ScalarMinimizationLineSearch_U.hpp line:62
struct PyCallBack_ROL_ScalarMinimizationLineSearch_U_double_t : public ROL::ScalarMinimizationLineSearch_U<double> {
	using ROL::ScalarMinimizationLineSearch_U<double>::ScalarMinimizationLineSearch_U;

	void initialize(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::ScalarMinimizationLineSearch_U<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ScalarMinimizationLineSearch_U::initialize(a0, a1);
	}
	void run(double & a0, double & a1, int & a2, int & a3, const double & a4, const class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Objective<double> & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::ScalarMinimizationLineSearch_U<double> *>(this), "run");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ScalarMinimizationLineSearch_U::run(a0, a1, a2, a3, a4, a5, a6, a7);
	}
	bool status(const enum ROL::ELineSearchU a0, int & a1, int & a2, const double a3, const double a4, const double a5, const double a6, const class ROL::Vector<double> & a7, const class ROL::Vector<double> & a8, class ROL::Objective<double> & a9) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::ScalarMinimizationLineSearch_U<double> *>(this), "status");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return LineSearch_U::status(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9);
	}
	double getInitialAlpha(int & a0, int & a1, const double a2, const double a3, const class ROL::Vector<double> & a4, const class ROL::Vector<double> & a5, class ROL::Objective<double> & a6) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::ScalarMinimizationLineSearch_U<double> *>(this), "getInitialAlpha");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return LineSearch_U::getInitialAlpha(a0, a1, a2, a3, a4, a5, a6);
	}
};

// ROL::DescentDirection_U file:ROL_DescentDirection_U.hpp line:58
struct PyCallBack_ROL_DescentDirection_U_double_t : public ROL::DescentDirection_U<double> {
	using ROL::DescentDirection_U<double>::DescentDirection_U;

	void initialize(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::DescentDirection_U<double> *>(this), "initialize");
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
	void compute(class ROL::Vector<double> & a0, double & a1, double & a2, int & a3, int & a4, const class ROL::Vector<double> & a5, const class ROL::Vector<double> & a6, class ROL::Objective<double> & a7) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::DescentDirection_U<double> *>(this), "compute");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6, a7);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"DescentDirection_U::compute\"");
	}
	void update(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1, const class ROL::Vector<double> & a2, const class ROL::Vector<double> & a3, const double a4, const int a5) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::DescentDirection_U<double> *>(this), "update");
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
	std::string printName() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::DescentDirection_U<double> *>(this), "printName");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<std::string>::value) {
				static pybind11::detail::override_caster_t<std::string> caster;
				return pybind11::detail::cast_ref<std::string>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<std::string>(std::move(o));
		}
		return DescentDirection_U::printName();
	}
};

void bind_ROL_BundleStatusTest(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // ROL::BundleStatusTest file:ROL_BundleStatusTest.hpp line:53
		pybind11::class_<ROL::BundleStatusTest<double>, std::shared_ptr<ROL::BundleStatusTest<double>>, PyCallBack_ROL_BundleStatusTest_double_t, ROL::StatusTest<double>> cl(M("ROL"), "BundleStatusTest_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("parlist") );

		cl.def( pybind11::init( [](){ return new ROL::BundleStatusTest<double>(); }, [](){ return new PyCallBack_ROL_BundleStatusTest_double_t(); } ), "doc");
		cl.def( pybind11::init( [](double const & a0){ return new ROL::BundleStatusTest<double>(a0); }, [](double const & a0){ return new PyCallBack_ROL_BundleStatusTest_double_t(a0); } ), "doc");
		cl.def( pybind11::init<double, int>(), pybind11::arg("tol"), pybind11::arg("max_iter") );

		cl.def( pybind11::init( [](PyCallBack_ROL_BundleStatusTest_double_t const &o){ return new PyCallBack_ROL_BundleStatusTest_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::BundleStatusTest<double> const &o){ return new ROL::BundleStatusTest<double>(o); } ) );
		cl.def("check", (bool (ROL::BundleStatusTest<double>::*)(struct ROL::AlgorithmState<double> &)) &ROL::BundleStatusTest<double>::check, "C++: ROL::BundleStatusTest<double>::check(struct ROL::AlgorithmState<double> &) --> bool", pybind11::arg("state"));
		cl.def("assign", (class ROL::BundleStatusTest<double> & (ROL::BundleStatusTest<double>::*)(const class ROL::BundleStatusTest<double> &)) &ROL::BundleStatusTest<double>::operator=, "C++: ROL::BundleStatusTest<double>::operator=(const class ROL::BundleStatusTest<double> &) --> class ROL::BundleStatusTest<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("check", (bool (ROL::StatusTest<double>::*)(struct ROL::AlgorithmState<double> &)) &ROL::StatusTest<double>::check, "C++: ROL::StatusTest<double>::check(struct ROL::AlgorithmState<double> &) --> bool", pybind11::arg("state"));
		cl.def("assign", (class ROL::StatusTest<double> & (ROL::StatusTest<double>::*)(const class ROL::StatusTest<double> &)) &ROL::StatusTest<double>::operator=, "C++: ROL::StatusTest<double>::operator=(const class ROL::StatusTest<double> &) --> class ROL::StatusTest<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::Bundle_U_AS file:ROL_Bundle_U_AS.hpp line:56
		pybind11::class_<ROL::Bundle_U_AS<double>, std::shared_ptr<ROL::Bundle_U_AS<double>>, PyCallBack_ROL_Bundle_U_AS_double_t, ROL::Bundle_U<double>> cl(M("ROL"), "Bundle_U_AS_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new ROL::Bundle_U_AS<double>(); }, [](){ return new PyCallBack_ROL_Bundle_U_AS_double_t(); } ), "doc");
		cl.def( pybind11::init( [](const unsigned int & a0){ return new ROL::Bundle_U_AS<double>(a0); }, [](const unsigned int & a0){ return new PyCallBack_ROL_Bundle_U_AS_double_t(a0); } ), "doc");
		cl.def( pybind11::init( [](const unsigned int & a0, const double & a1){ return new ROL::Bundle_U_AS<double>(a0, a1); }, [](const unsigned int & a0, const double & a1){ return new PyCallBack_ROL_Bundle_U_AS_double_t(a0, a1); } ), "doc");
		cl.def( pybind11::init( [](const unsigned int & a0, const double & a1, const double & a2){ return new ROL::Bundle_U_AS<double>(a0, a1, a2); }, [](const unsigned int & a0, const double & a1, const double & a2){ return new PyCallBack_ROL_Bundle_U_AS_double_t(a0, a1, a2); } ), "doc");
		cl.def( pybind11::init<const unsigned int, const double, const double, const unsigned int>(), pybind11::arg("maxSize"), pybind11::arg("coeff"), pybind11::arg("omega"), pybind11::arg("remSize") );

		cl.def( pybind11::init( [](PyCallBack_ROL_Bundle_U_AS_double_t const &o){ return new PyCallBack_ROL_Bundle_U_AS_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::Bundle_U_AS<double> const &o){ return new ROL::Bundle_U_AS<double>(o); } ) );
		cl.def("initialize", (void (ROL::Bundle_U_AS<double>::*)(const class ROL::Vector<double> &)) &ROL::Bundle_U_AS<double>::initialize, "C++: ROL::Bundle_U_AS<double>::initialize(const class ROL::Vector<double> &) --> void", pybind11::arg("g"));
		cl.def("solveDual", [](ROL::Bundle_U_AS<double> &o, const double & a0) -> unsigned int { return o.solveDual(a0); }, "", pybind11::arg("t"));
		cl.def("solveDual", [](ROL::Bundle_U_AS<double> &o, const double & a0, const unsigned int & a1) -> unsigned int { return o.solveDual(a0, a1); }, "", pybind11::arg("t"), pybind11::arg("maxit"));
		cl.def("solveDual", (unsigned int (ROL::Bundle_U_AS<double>::*)(const double, const unsigned int, const double)) &ROL::Bundle_U_AS<double>::solveDual, "C++: ROL::Bundle_U_AS<double>::solveDual(const double, const unsigned int, const double) --> unsigned int", pybind11::arg("t"), pybind11::arg("maxit"), pybind11::arg("tol"));
		cl.def("assign", (class ROL::Bundle_U_AS<double> & (ROL::Bundle_U_AS<double>::*)(const class ROL::Bundle_U_AS<double> &)) &ROL::Bundle_U_AS<double>::operator=, "C++: ROL::Bundle_U_AS<double>::operator=(const class ROL::Bundle_U_AS<double> &) --> class ROL::Bundle_U_AS<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("initialize", (void (ROL::Bundle_U<double>::*)(const class ROL::Vector<double> &)) &ROL::Bundle_U<double>::initialize, "C++: ROL::Bundle_U<double>::initialize(const class ROL::Vector<double> &) --> void", pybind11::arg("g"));
		cl.def("solveDual", [](ROL::Bundle_U<double> &o, const double & a0) -> unsigned int { return o.solveDual(a0); }, "", pybind11::arg("t"));
		cl.def("solveDual", [](ROL::Bundle_U<double> &o, const double & a0, const unsigned int & a1) -> unsigned int { return o.solveDual(a0, a1); }, "", pybind11::arg("t"), pybind11::arg("maxit"));
		cl.def("solveDual", (unsigned int (ROL::Bundle_U<double>::*)(const double, const unsigned int, const double)) &ROL::Bundle_U<double>::solveDual, "C++: ROL::Bundle_U<double>::solveDual(const double, const unsigned int, const double) --> unsigned int", pybind11::arg("t"), pybind11::arg("maxit"), pybind11::arg("tol"));
		cl.def("linearizationError", (const double (ROL::Bundle_U<double>::*)(const unsigned int) const) &ROL::Bundle_U<double>::linearizationError, "C++: ROL::Bundle_U<double>::linearizationError(const unsigned int) const --> const double", pybind11::arg("i"));
		cl.def("distanceMeasure", (const double (ROL::Bundle_U<double>::*)(const unsigned int) const) &ROL::Bundle_U<double>::distanceMeasure, "C++: ROL::Bundle_U<double>::distanceMeasure(const unsigned int) const --> const double", pybind11::arg("i"));
		cl.def("subgradient", (const class ROL::Vector<double> & (ROL::Bundle_U<double>::*)(const unsigned int) const) &ROL::Bundle_U<double>::subgradient, "C++: ROL::Bundle_U<double>::subgradient(const unsigned int) const --> const class ROL::Vector<double> &", pybind11::return_value_policy::automatic, pybind11::arg("i"));
		cl.def("getDualVariable", (const double (ROL::Bundle_U<double>::*)(const unsigned int) const) &ROL::Bundle_U<double>::getDualVariable, "C++: ROL::Bundle_U<double>::getDualVariable(const unsigned int) const --> const double", pybind11::arg("i"));
		cl.def("setDualVariable", (void (ROL::Bundle_U<double>::*)(const unsigned int, const double)) &ROL::Bundle_U<double>::setDualVariable, "C++: ROL::Bundle_U<double>::setDualVariable(const unsigned int, const double) --> void", pybind11::arg("i"), pybind11::arg("val"));
		cl.def("resetDualVariables", (void (ROL::Bundle_U<double>::*)()) &ROL::Bundle_U<double>::resetDualVariables, "C++: ROL::Bundle_U<double>::resetDualVariables() --> void");
		cl.def("computeAlpha", (const double (ROL::Bundle_U<double>::*)(const double, const double) const) &ROL::Bundle_U<double>::computeAlpha, "C++: ROL::Bundle_U<double>::computeAlpha(const double, const double) const --> const double", pybind11::arg("dm"), pybind11::arg("le"));
		cl.def("alpha", (const double (ROL::Bundle_U<double>::*)(const unsigned int) const) &ROL::Bundle_U<double>::alpha, "C++: ROL::Bundle_U<double>::alpha(const unsigned int) const --> const double", pybind11::arg("i"));
		cl.def("size", (unsigned int (ROL::Bundle_U<double>::*)() const) &ROL::Bundle_U<double>::size, "C++: ROL::Bundle_U<double>::size() const --> unsigned int");
		cl.def("aggregate", (void (ROL::Bundle_U<double>::*)(class ROL::Vector<double> &, double &, double &) const) &ROL::Bundle_U<double>::aggregate, "C++: ROL::Bundle_U<double>::aggregate(class ROL::Vector<double> &, double &, double &) const --> void", pybind11::arg("aggSubGrad"), pybind11::arg("aggLinErr"), pybind11::arg("aggDistMeas"));
		cl.def("reset", (void (ROL::Bundle_U<double>::*)(const class ROL::Vector<double> &, const double, const double)) &ROL::Bundle_U<double>::reset, "C++: ROL::Bundle_U<double>::reset(const class ROL::Vector<double> &, const double, const double) --> void", pybind11::arg("g"), pybind11::arg("le"), pybind11::arg("dm"));
		cl.def("update", (void (ROL::Bundle_U<double>::*)(const bool, const double, const double, const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::Bundle_U<double>::update, "C++: ROL::Bundle_U<double>::update(const bool, const double, const double, const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("flag"), pybind11::arg("linErr"), pybind11::arg("distMeas"), pybind11::arg("g"), pybind11::arg("s"));
		cl.def("assign", (class ROL::Bundle_U<double> & (ROL::Bundle_U<double>::*)(const class ROL::Bundle_U<double> &)) &ROL::Bundle_U<double>::operator=, "C++: ROL::Bundle_U<double>::operator=(const class ROL::Bundle_U<double> &) --> class ROL::Bundle_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::Bundle_U_TT file:ROL_Bundle_U_TT.hpp line:61
		pybind11::class_<ROL::Bundle_U_TT<double>, std::shared_ptr<ROL::Bundle_U_TT<double>>, PyCallBack_ROL_Bundle_U_TT_double_t, ROL::Bundle_U<double>> cl(M("ROL"), "Bundle_U_TT_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new ROL::Bundle_U_TT<double>(); }, [](){ return new PyCallBack_ROL_Bundle_U_TT_double_t(); } ), "doc");
		cl.def( pybind11::init( [](const unsigned int & a0){ return new ROL::Bundle_U_TT<double>(a0); }, [](const unsigned int & a0){ return new PyCallBack_ROL_Bundle_U_TT_double_t(a0); } ), "doc");
		cl.def( pybind11::init( [](const unsigned int & a0, const double & a1){ return new ROL::Bundle_U_TT<double>(a0, a1); }, [](const unsigned int & a0, const double & a1){ return new PyCallBack_ROL_Bundle_U_TT_double_t(a0, a1); } ), "doc");
		cl.def( pybind11::init( [](const unsigned int & a0, const double & a1, const double & a2){ return new ROL::Bundle_U_TT<double>(a0, a1, a2); }, [](const unsigned int & a0, const double & a1, const double & a2){ return new PyCallBack_ROL_Bundle_U_TT_double_t(a0, a1, a2); } ), "doc");
		cl.def( pybind11::init<const unsigned int, const double, const double, const unsigned int>(), pybind11::arg("maxSize"), pybind11::arg("coeff"), pybind11::arg("omega"), pybind11::arg("remSize") );

		cl.def( pybind11::init( [](PyCallBack_ROL_Bundle_U_TT_double_t const &o){ return new PyCallBack_ROL_Bundle_U_TT_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::Bundle_U_TT<double> const &o){ return new ROL::Bundle_U_TT<double>(o); } ) );
		cl.def("solveDual", [](ROL::Bundle_U_TT<double> &o, const double & a0) -> unsigned int { return o.solveDual(a0); }, "", pybind11::arg("t"));
		cl.def("solveDual", [](ROL::Bundle_U_TT<double> &o, const double & a0, const unsigned int & a1) -> unsigned int { return o.solveDual(a0, a1); }, "", pybind11::arg("t"), pybind11::arg("maxit"));
		cl.def("solveDual", (unsigned int (ROL::Bundle_U_TT<double>::*)(const double, const unsigned int, const double)) &ROL::Bundle_U_TT<double>::solveDual, "C++: ROL::Bundle_U_TT<double>::solveDual(const double, const unsigned int, const double) --> unsigned int", pybind11::arg("t"), pybind11::arg("maxit"), pybind11::arg("tol"));
		cl.def("assign", (class ROL::Bundle_U_TT<double> & (ROL::Bundle_U_TT<double>::*)(const class ROL::Bundle_U_TT<double> &)) &ROL::Bundle_U_TT<double>::operator=, "C++: ROL::Bundle_U_TT<double>::operator=(const class ROL::Bundle_U_TT<double> &) --> class ROL::Bundle_U_TT<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("initialize", (void (ROL::Bundle_U<double>::*)(const class ROL::Vector<double> &)) &ROL::Bundle_U<double>::initialize, "C++: ROL::Bundle_U<double>::initialize(const class ROL::Vector<double> &) --> void", pybind11::arg("g"));
		cl.def("solveDual", [](ROL::Bundle_U<double> &o, const double & a0) -> unsigned int { return o.solveDual(a0); }, "", pybind11::arg("t"));
		cl.def("solveDual", [](ROL::Bundle_U<double> &o, const double & a0, const unsigned int & a1) -> unsigned int { return o.solveDual(a0, a1); }, "", pybind11::arg("t"), pybind11::arg("maxit"));
		cl.def("solveDual", (unsigned int (ROL::Bundle_U<double>::*)(const double, const unsigned int, const double)) &ROL::Bundle_U<double>::solveDual, "C++: ROL::Bundle_U<double>::solveDual(const double, const unsigned int, const double) --> unsigned int", pybind11::arg("t"), pybind11::arg("maxit"), pybind11::arg("tol"));
		cl.def("linearizationError", (const double (ROL::Bundle_U<double>::*)(const unsigned int) const) &ROL::Bundle_U<double>::linearizationError, "C++: ROL::Bundle_U<double>::linearizationError(const unsigned int) const --> const double", pybind11::arg("i"));
		cl.def("distanceMeasure", (const double (ROL::Bundle_U<double>::*)(const unsigned int) const) &ROL::Bundle_U<double>::distanceMeasure, "C++: ROL::Bundle_U<double>::distanceMeasure(const unsigned int) const --> const double", pybind11::arg("i"));
		cl.def("subgradient", (const class ROL::Vector<double> & (ROL::Bundle_U<double>::*)(const unsigned int) const) &ROL::Bundle_U<double>::subgradient, "C++: ROL::Bundle_U<double>::subgradient(const unsigned int) const --> const class ROL::Vector<double> &", pybind11::return_value_policy::automatic, pybind11::arg("i"));
		cl.def("getDualVariable", (const double (ROL::Bundle_U<double>::*)(const unsigned int) const) &ROL::Bundle_U<double>::getDualVariable, "C++: ROL::Bundle_U<double>::getDualVariable(const unsigned int) const --> const double", pybind11::arg("i"));
		cl.def("setDualVariable", (void (ROL::Bundle_U<double>::*)(const unsigned int, const double)) &ROL::Bundle_U<double>::setDualVariable, "C++: ROL::Bundle_U<double>::setDualVariable(const unsigned int, const double) --> void", pybind11::arg("i"), pybind11::arg("val"));
		cl.def("resetDualVariables", (void (ROL::Bundle_U<double>::*)()) &ROL::Bundle_U<double>::resetDualVariables, "C++: ROL::Bundle_U<double>::resetDualVariables() --> void");
		cl.def("computeAlpha", (const double (ROL::Bundle_U<double>::*)(const double, const double) const) &ROL::Bundle_U<double>::computeAlpha, "C++: ROL::Bundle_U<double>::computeAlpha(const double, const double) const --> const double", pybind11::arg("dm"), pybind11::arg("le"));
		cl.def("alpha", (const double (ROL::Bundle_U<double>::*)(const unsigned int) const) &ROL::Bundle_U<double>::alpha, "C++: ROL::Bundle_U<double>::alpha(const unsigned int) const --> const double", pybind11::arg("i"));
		cl.def("size", (unsigned int (ROL::Bundle_U<double>::*)() const) &ROL::Bundle_U<double>::size, "C++: ROL::Bundle_U<double>::size() const --> unsigned int");
		cl.def("aggregate", (void (ROL::Bundle_U<double>::*)(class ROL::Vector<double> &, double &, double &) const) &ROL::Bundle_U<double>::aggregate, "C++: ROL::Bundle_U<double>::aggregate(class ROL::Vector<double> &, double &, double &) const --> void", pybind11::arg("aggSubGrad"), pybind11::arg("aggLinErr"), pybind11::arg("aggDistMeas"));
		cl.def("reset", (void (ROL::Bundle_U<double>::*)(const class ROL::Vector<double> &, const double, const double)) &ROL::Bundle_U<double>::reset, "C++: ROL::Bundle_U<double>::reset(const class ROL::Vector<double> &, const double, const double) --> void", pybind11::arg("g"), pybind11::arg("le"), pybind11::arg("dm"));
		cl.def("update", (void (ROL::Bundle_U<double>::*)(const bool, const double, const double, const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::Bundle_U<double>::update, "C++: ROL::Bundle_U<double>::update(const bool, const double, const double, const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("flag"), pybind11::arg("linErr"), pybind11::arg("distMeas"), pybind11::arg("g"), pybind11::arg("s"));
		cl.def("assign", (class ROL::Bundle_U<double> & (ROL::Bundle_U<double>::*)(const class ROL::Bundle_U<double> &)) &ROL::Bundle_U<double>::operator=, "C++: ROL::Bundle_U<double>::operator=(const class ROL::Bundle_U<double> &) --> class ROL::Bundle_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::IterationScaling_U file:ROL_IterationScaling_U.hpp line:56
		pybind11::class_<ROL::IterationScaling_U<double>, std::shared_ptr<ROL::IterationScaling_U<double>>, PyCallBack_ROL_IterationScaling_U_double_t, ROL::LineSearch_U<double>> cl(M("ROL"), "IterationScaling_U_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("parlist") );

		cl.def( pybind11::init( [](PyCallBack_ROL_IterationScaling_U_double_t const &o){ return new PyCallBack_ROL_IterationScaling_U_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::IterationScaling_U<double> const &o){ return new ROL::IterationScaling_U<double>(o); } ) );
		cl.def("initialize", (void (ROL::IterationScaling_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::IterationScaling_U<double>::initialize, "C++: ROL::IterationScaling_U<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"));
		cl.def("run", (void (ROL::IterationScaling_U<double>::*)(double &, double &, int &, int &, const double &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &)) &ROL::IterationScaling_U<double>::run, "C++: ROL::IterationScaling_U<double>::run(double &, double &, int &, int &, const double &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &) --> void", pybind11::arg("alpha"), pybind11::arg("fval"), pybind11::arg("ls_neval"), pybind11::arg("ls_ngrad"), pybind11::arg("gs"), pybind11::arg("s"), pybind11::arg("x"), pybind11::arg("obj"));
		cl.def("assign", (class ROL::IterationScaling_U<double> & (ROL::IterationScaling_U<double>::*)(const class ROL::IterationScaling_U<double> &)) &ROL::IterationScaling_U<double>::operator=, "C++: ROL::IterationScaling_U<double>::operator=(const class ROL::IterationScaling_U<double> &) --> class ROL::IterationScaling_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("initialize", (void (ROL::LineSearch_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::LineSearch_U<double>::initialize, "C++: ROL::LineSearch_U<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"));
		cl.def("run", (void (ROL::LineSearch_U<double>::*)(double &, double &, int &, int &, const double &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &)) &ROL::LineSearch_U<double>::run, "C++: ROL::LineSearch_U<double>::run(double &, double &, int &, int &, const double &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &) --> void", pybind11::arg("alpha"), pybind11::arg("fval"), pybind11::arg("ls_neval"), pybind11::arg("ls_ngrad"), pybind11::arg("gs"), pybind11::arg("s"), pybind11::arg("x"), pybind11::arg("obj"));
		cl.def("setMaxitUpdate", (void (ROL::LineSearch_U<double>::*)(double &, double &, const double &)) &ROL::LineSearch_U<double>::setMaxitUpdate, "C++: ROL::LineSearch_U<double>::setMaxitUpdate(double &, double &, const double &) --> void", pybind11::arg("alpha"), pybind11::arg("fnew"), pybind11::arg("fold"));
		cl.def("assign", (class ROL::LineSearch_U<double> & (ROL::LineSearch_U<double>::*)(const class ROL::LineSearch_U<double> &)) &ROL::LineSearch_U<double>::operator=, "C++: ROL::LineSearch_U<double>::operator=(const class ROL::LineSearch_U<double> &) --> class ROL::LineSearch_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::PathBasedTargetLevel_U file:ROL_PathBasedTargetLevel_U.hpp line:56
		pybind11::class_<ROL::PathBasedTargetLevel_U<double>, std::shared_ptr<ROL::PathBasedTargetLevel_U<double>>, PyCallBack_ROL_PathBasedTargetLevel_U_double_t, ROL::LineSearch_U<double>> cl(M("ROL"), "PathBasedTargetLevel_U_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("parlist") );

		cl.def( pybind11::init( [](PyCallBack_ROL_PathBasedTargetLevel_U_double_t const &o){ return new PyCallBack_ROL_PathBasedTargetLevel_U_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::PathBasedTargetLevel_U<double> const &o){ return new ROL::PathBasedTargetLevel_U<double>(o); } ) );
		cl.def("initialize", (void (ROL::PathBasedTargetLevel_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::PathBasedTargetLevel_U<double>::initialize, "C++: ROL::PathBasedTargetLevel_U<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"));
		cl.def("run", (void (ROL::PathBasedTargetLevel_U<double>::*)(double &, double &, int &, int &, const double &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &)) &ROL::PathBasedTargetLevel_U<double>::run, "C++: ROL::PathBasedTargetLevel_U<double>::run(double &, double &, int &, int &, const double &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &) --> void", pybind11::arg("alpha"), pybind11::arg("fval"), pybind11::arg("ls_neval"), pybind11::arg("ls_ngrad"), pybind11::arg("gs"), pybind11::arg("s"), pybind11::arg("x"), pybind11::arg("obj"));
		cl.def("assign", (class ROL::PathBasedTargetLevel_U<double> & (ROL::PathBasedTargetLevel_U<double>::*)(const class ROL::PathBasedTargetLevel_U<double> &)) &ROL::PathBasedTargetLevel_U<double>::operator=, "C++: ROL::PathBasedTargetLevel_U<double>::operator=(const class ROL::PathBasedTargetLevel_U<double> &) --> class ROL::PathBasedTargetLevel_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("initialize", (void (ROL::LineSearch_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::LineSearch_U<double>::initialize, "C++: ROL::LineSearch_U<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"));
		cl.def("run", (void (ROL::LineSearch_U<double>::*)(double &, double &, int &, int &, const double &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &)) &ROL::LineSearch_U<double>::run, "C++: ROL::LineSearch_U<double>::run(double &, double &, int &, int &, const double &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &) --> void", pybind11::arg("alpha"), pybind11::arg("fval"), pybind11::arg("ls_neval"), pybind11::arg("ls_ngrad"), pybind11::arg("gs"), pybind11::arg("s"), pybind11::arg("x"), pybind11::arg("obj"));
		cl.def("setMaxitUpdate", (void (ROL::LineSearch_U<double>::*)(double &, double &, const double &)) &ROL::LineSearch_U<double>::setMaxitUpdate, "C++: ROL::LineSearch_U<double>::setMaxitUpdate(double &, double &, const double &) --> void", pybind11::arg("alpha"), pybind11::arg("fnew"), pybind11::arg("fold"));
		cl.def("assign", (class ROL::LineSearch_U<double> & (ROL::LineSearch_U<double>::*)(const class ROL::LineSearch_U<double> &)) &ROL::LineSearch_U<double>::operator=, "C++: ROL::LineSearch_U<double>::operator=(const class ROL::LineSearch_U<double> &) --> class ROL::LineSearch_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::BackTracking_U file:ROL_BackTracking_U.hpp line:56
		pybind11::class_<ROL::BackTracking_U<double>, std::shared_ptr<ROL::BackTracking_U<double>>, PyCallBack_ROL_BackTracking_U_double_t, ROL::LineSearch_U<double>> cl(M("ROL"), "BackTracking_U_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("parlist") );

		cl.def( pybind11::init( [](PyCallBack_ROL_BackTracking_U_double_t const &o){ return new PyCallBack_ROL_BackTracking_U_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::BackTracking_U<double> const &o){ return new ROL::BackTracking_U<double>(o); } ) );
		cl.def("initialize", (void (ROL::BackTracking_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::BackTracking_U<double>::initialize, "C++: ROL::BackTracking_U<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"));
		cl.def("run", (void (ROL::BackTracking_U<double>::*)(double &, double &, int &, int &, const double &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &)) &ROL::BackTracking_U<double>::run, "C++: ROL::BackTracking_U<double>::run(double &, double &, int &, int &, const double &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &) --> void", pybind11::arg("alpha"), pybind11::arg("fval"), pybind11::arg("ls_neval"), pybind11::arg("ls_ngrad"), pybind11::arg("gs"), pybind11::arg("s"), pybind11::arg("x"), pybind11::arg("obj"));
		cl.def("assign", (class ROL::BackTracking_U<double> & (ROL::BackTracking_U<double>::*)(const class ROL::BackTracking_U<double> &)) &ROL::BackTracking_U<double>::operator=, "C++: ROL::BackTracking_U<double>::operator=(const class ROL::BackTracking_U<double> &) --> class ROL::BackTracking_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("initialize", (void (ROL::LineSearch_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::LineSearch_U<double>::initialize, "C++: ROL::LineSearch_U<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"));
		cl.def("run", (void (ROL::LineSearch_U<double>::*)(double &, double &, int &, int &, const double &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &)) &ROL::LineSearch_U<double>::run, "C++: ROL::LineSearch_U<double>::run(double &, double &, int &, int &, const double &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &) --> void", pybind11::arg("alpha"), pybind11::arg("fval"), pybind11::arg("ls_neval"), pybind11::arg("ls_ngrad"), pybind11::arg("gs"), pybind11::arg("s"), pybind11::arg("x"), pybind11::arg("obj"));
		cl.def("setMaxitUpdate", (void (ROL::LineSearch_U<double>::*)(double &, double &, const double &)) &ROL::LineSearch_U<double>::setMaxitUpdate, "C++: ROL::LineSearch_U<double>::setMaxitUpdate(double &, double &, const double &) --> void", pybind11::arg("alpha"), pybind11::arg("fnew"), pybind11::arg("fold"));
		cl.def("assign", (class ROL::LineSearch_U<double> & (ROL::LineSearch_U<double>::*)(const class ROL::LineSearch_U<double> &)) &ROL::LineSearch_U<double>::operator=, "C++: ROL::LineSearch_U<double>::operator=(const class ROL::LineSearch_U<double> &) --> class ROL::LineSearch_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::CubicInterp_U file:ROL_CubicInterp_U.hpp line:56
		pybind11::class_<ROL::CubicInterp_U<double>, std::shared_ptr<ROL::CubicInterp_U<double>>, PyCallBack_ROL_CubicInterp_U_double_t, ROL::LineSearch_U<double>> cl(M("ROL"), "CubicInterp_U_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("parlist") );

		cl.def( pybind11::init( [](PyCallBack_ROL_CubicInterp_U_double_t const &o){ return new PyCallBack_ROL_CubicInterp_U_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::CubicInterp_U<double> const &o){ return new ROL::CubicInterp_U<double>(o); } ) );
		cl.def("initialize", (void (ROL::CubicInterp_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::CubicInterp_U<double>::initialize, "C++: ROL::CubicInterp_U<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"));
		cl.def("run", (void (ROL::CubicInterp_U<double>::*)(double &, double &, int &, int &, const double &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &)) &ROL::CubicInterp_U<double>::run, "C++: ROL::CubicInterp_U<double>::run(double &, double &, int &, int &, const double &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &) --> void", pybind11::arg("alpha"), pybind11::arg("fval"), pybind11::arg("ls_neval"), pybind11::arg("ls_ngrad"), pybind11::arg("gs"), pybind11::arg("s"), pybind11::arg("x"), pybind11::arg("obj"));
		cl.def("assign", (class ROL::CubicInterp_U<double> & (ROL::CubicInterp_U<double>::*)(const class ROL::CubicInterp_U<double> &)) &ROL::CubicInterp_U<double>::operator=, "C++: ROL::CubicInterp_U<double>::operator=(const class ROL::CubicInterp_U<double> &) --> class ROL::CubicInterp_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("initialize", (void (ROL::LineSearch_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::LineSearch_U<double>::initialize, "C++: ROL::LineSearch_U<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"));
		cl.def("run", (void (ROL::LineSearch_U<double>::*)(double &, double &, int &, int &, const double &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &)) &ROL::LineSearch_U<double>::run, "C++: ROL::LineSearch_U<double>::run(double &, double &, int &, int &, const double &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &) --> void", pybind11::arg("alpha"), pybind11::arg("fval"), pybind11::arg("ls_neval"), pybind11::arg("ls_ngrad"), pybind11::arg("gs"), pybind11::arg("s"), pybind11::arg("x"), pybind11::arg("obj"));
		cl.def("setMaxitUpdate", (void (ROL::LineSearch_U<double>::*)(double &, double &, const double &)) &ROL::LineSearch_U<double>::setMaxitUpdate, "C++: ROL::LineSearch_U<double>::setMaxitUpdate(double &, double &, const double &) --> void", pybind11::arg("alpha"), pybind11::arg("fnew"), pybind11::arg("fold"));
		cl.def("assign", (class ROL::LineSearch_U<double> & (ROL::LineSearch_U<double>::*)(const class ROL::LineSearch_U<double> &)) &ROL::LineSearch_U<double>::operator=, "C++: ROL::LineSearch_U<double>::operator=(const class ROL::LineSearch_U<double> &) --> class ROL::LineSearch_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::ScalarMinimizationStatusTest file:ROL_ScalarMinimizationStatusTest.hpp line:54
		pybind11::class_<ROL::ScalarMinimizationStatusTest<double>, std::shared_ptr<ROL::ScalarMinimizationStatusTest<double>>, PyCallBack_ROL_ScalarMinimizationStatusTest_double_t> cl(M("ROL"), "ScalarMinimizationStatusTest_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](PyCallBack_ROL_ScalarMinimizationStatusTest_double_t const &o){ return new PyCallBack_ROL_ScalarMinimizationStatusTest_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::ScalarMinimizationStatusTest<double> const &o){ return new ROL::ScalarMinimizationStatusTest<double>(o); } ) );
		cl.def( pybind11::init( [](){ return new ROL::ScalarMinimizationStatusTest<double>(); }, [](){ return new PyCallBack_ROL_ScalarMinimizationStatusTest_double_t(); } ) );
		cl.def("check", [](ROL::ScalarMinimizationStatusTest<double> &o, double & a0, double & a1, double & a2, int & a3, int & a4) -> bool { return o.check(a0, a1, a2, a3, a4); }, "", pybind11::arg("x"), pybind11::arg("fx"), pybind11::arg("gx"), pybind11::arg("nfval"), pybind11::arg("ngrad"));
		cl.def("check", (bool (ROL::ScalarMinimizationStatusTest<double>::*)(double &, double &, double &, int &, int &, const bool)) &ROL::ScalarMinimizationStatusTest<double>::check, "C++: ROL::ScalarMinimizationStatusTest<double>::check(double &, double &, double &, int &, int &, const bool) --> bool", pybind11::arg("x"), pybind11::arg("fx"), pybind11::arg("gx"), pybind11::arg("nfval"), pybind11::arg("ngrad"), pybind11::arg("deriv"));
		cl.def("assign", (class ROL::ScalarMinimizationStatusTest<double> & (ROL::ScalarMinimizationStatusTest<double>::*)(const class ROL::ScalarMinimizationStatusTest<double> &)) &ROL::ScalarMinimizationStatusTest<double>::operator=, "C++: ROL::ScalarMinimizationStatusTest<double>::operator=(const class ROL::ScalarMinimizationStatusTest<double> &) --> class ROL::ScalarMinimizationStatusTest<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::ScalarMinimization file:ROL_ScalarMinimization.hpp line:59
		pybind11::class_<ROL::ScalarMinimization<double>, std::shared_ptr<ROL::ScalarMinimization<double>>, PyCallBack_ROL_ScalarMinimization_double_t> cl(M("ROL"), "ScalarMinimization_double_t", "", pybind11::module_local());
		cl.def(pybind11::init<PyCallBack_ROL_ScalarMinimization_double_t const &>());
		cl.def( pybind11::init( [](){ return new PyCallBack_ROL_ScalarMinimization_double_t(); } ) );
		cl.def("run", (void (ROL::ScalarMinimization<double>::*)(double &, double &, int &, int &, class ROL::ScalarFunction<double> &, const double, const double) const) &ROL::ScalarMinimization<double>::run, "C++: ROL::ScalarMinimization<double>::run(double &, double &, int &, int &, class ROL::ScalarFunction<double> &, const double, const double) const --> void", pybind11::arg("fx"), pybind11::arg("x"), pybind11::arg("nfval"), pybind11::arg("ngrad"), pybind11::arg("f"), pybind11::arg("A"), pybind11::arg("B"));
		cl.def("run", (void (ROL::ScalarMinimization<double>::*)(double &, double &, int &, int &, class ROL::ScalarFunction<double> &, const double, const double, class ROL::ScalarMinimizationStatusTest<double> &) const) &ROL::ScalarMinimization<double>::run, "C++: ROL::ScalarMinimization<double>::run(double &, double &, int &, int &, class ROL::ScalarFunction<double> &, const double, const double, class ROL::ScalarMinimizationStatusTest<double> &) const --> void", pybind11::arg("fx"), pybind11::arg("x"), pybind11::arg("nfval"), pybind11::arg("ngrad"), pybind11::arg("f"), pybind11::arg("A"), pybind11::arg("B"), pybind11::arg("test"));
		cl.def("assign", (class ROL::ScalarMinimization<double> & (ROL::ScalarMinimization<double>::*)(const class ROL::ScalarMinimization<double> &)) &ROL::ScalarMinimization<double>::operator=, "C++: ROL::ScalarMinimization<double>::operator=(const class ROL::ScalarMinimization<double> &) --> class ROL::ScalarMinimization<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::BrentsScalarMinimization file:ROL_BrentsScalarMinimization.hpp line:60
		pybind11::class_<ROL::BrentsScalarMinimization<double>, std::shared_ptr<ROL::BrentsScalarMinimization<double>>, PyCallBack_ROL_BrentsScalarMinimization_double_t, ROL::ScalarMinimization<double>> cl(M("ROL"), "BrentsScalarMinimization_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("parlist") );

		cl.def( pybind11::init( [](PyCallBack_ROL_BrentsScalarMinimization_double_t const &o){ return new PyCallBack_ROL_BrentsScalarMinimization_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::BrentsScalarMinimization<double> const &o){ return new ROL::BrentsScalarMinimization<double>(o); } ) );
		cl.def("run", (void (ROL::BrentsScalarMinimization<double>::*)(double &, double &, int &, int &, class ROL::ScalarFunction<double> &, const double, const double, class ROL::ScalarMinimizationStatusTest<double> &) const) &ROL::BrentsScalarMinimization<double>::run, "C++: ROL::BrentsScalarMinimization<double>::run(double &, double &, int &, int &, class ROL::ScalarFunction<double> &, const double, const double, class ROL::ScalarMinimizationStatusTest<double> &) const --> void", pybind11::arg("fx"), pybind11::arg("x"), pybind11::arg("nfval"), pybind11::arg("ngrad"), pybind11::arg("f"), pybind11::arg("A"), pybind11::arg("B"), pybind11::arg("test"));
		cl.def("assign", (class ROL::BrentsScalarMinimization<double> & (ROL::BrentsScalarMinimization<double>::*)(const class ROL::BrentsScalarMinimization<double> &)) &ROL::BrentsScalarMinimization<double>::operator=, "C++: ROL::BrentsScalarMinimization<double>::operator=(const class ROL::BrentsScalarMinimization<double> &) --> class ROL::BrentsScalarMinimization<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("run", (void (ROL::ScalarMinimization<double>::*)(double &, double &, int &, int &, class ROL::ScalarFunction<double> &, const double, const double) const) &ROL::ScalarMinimization<double>::run, "C++: ROL::ScalarMinimization<double>::run(double &, double &, int &, int &, class ROL::ScalarFunction<double> &, const double, const double) const --> void", pybind11::arg("fx"), pybind11::arg("x"), pybind11::arg("nfval"), pybind11::arg("ngrad"), pybind11::arg("f"), pybind11::arg("A"), pybind11::arg("B"));
		cl.def("run", (void (ROL::ScalarMinimization<double>::*)(double &, double &, int &, int &, class ROL::ScalarFunction<double> &, const double, const double, class ROL::ScalarMinimizationStatusTest<double> &) const) &ROL::ScalarMinimization<double>::run, "C++: ROL::ScalarMinimization<double>::run(double &, double &, int &, int &, class ROL::ScalarFunction<double> &, const double, const double, class ROL::ScalarMinimizationStatusTest<double> &) const --> void", pybind11::arg("fx"), pybind11::arg("x"), pybind11::arg("nfval"), pybind11::arg("ngrad"), pybind11::arg("f"), pybind11::arg("A"), pybind11::arg("B"), pybind11::arg("test"));
		cl.def("assign", (class ROL::ScalarMinimization<double> & (ROL::ScalarMinimization<double>::*)(const class ROL::ScalarMinimization<double> &)) &ROL::ScalarMinimization<double>::operator=, "C++: ROL::ScalarMinimization<double>::operator=(const class ROL::ScalarMinimization<double> &) --> class ROL::ScalarMinimization<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::BisectionScalarMinimization file:ROL_BisectionScalarMinimization.hpp line:60
		pybind11::class_<ROL::BisectionScalarMinimization<double>, std::shared_ptr<ROL::BisectionScalarMinimization<double>>, PyCallBack_ROL_BisectionScalarMinimization_double_t, ROL::ScalarMinimization<double>> cl(M("ROL"), "BisectionScalarMinimization_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("parlist") );

		cl.def( pybind11::init( [](PyCallBack_ROL_BisectionScalarMinimization_double_t const &o){ return new PyCallBack_ROL_BisectionScalarMinimization_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::BisectionScalarMinimization<double> const &o){ return new ROL::BisectionScalarMinimization<double>(o); } ) );
		cl.def("run", (void (ROL::BisectionScalarMinimization<double>::*)(double &, double &, int &, int &, class ROL::ScalarFunction<double> &, const double, const double, class ROL::ScalarMinimizationStatusTest<double> &) const) &ROL::BisectionScalarMinimization<double>::run, "C++: ROL::BisectionScalarMinimization<double>::run(double &, double &, int &, int &, class ROL::ScalarFunction<double> &, const double, const double, class ROL::ScalarMinimizationStatusTest<double> &) const --> void", pybind11::arg("fx"), pybind11::arg("x"), pybind11::arg("nfval"), pybind11::arg("ngrad"), pybind11::arg("f"), pybind11::arg("A"), pybind11::arg("B"), pybind11::arg("test"));
		cl.def("assign", (class ROL::BisectionScalarMinimization<double> & (ROL::BisectionScalarMinimization<double>::*)(const class ROL::BisectionScalarMinimization<double> &)) &ROL::BisectionScalarMinimization<double>::operator=, "C++: ROL::BisectionScalarMinimization<double>::operator=(const class ROL::BisectionScalarMinimization<double> &) --> class ROL::BisectionScalarMinimization<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("run", (void (ROL::ScalarMinimization<double>::*)(double &, double &, int &, int &, class ROL::ScalarFunction<double> &, const double, const double) const) &ROL::ScalarMinimization<double>::run, "C++: ROL::ScalarMinimization<double>::run(double &, double &, int &, int &, class ROL::ScalarFunction<double> &, const double, const double) const --> void", pybind11::arg("fx"), pybind11::arg("x"), pybind11::arg("nfval"), pybind11::arg("ngrad"), pybind11::arg("f"), pybind11::arg("A"), pybind11::arg("B"));
		cl.def("run", (void (ROL::ScalarMinimization<double>::*)(double &, double &, int &, int &, class ROL::ScalarFunction<double> &, const double, const double, class ROL::ScalarMinimizationStatusTest<double> &) const) &ROL::ScalarMinimization<double>::run, "C++: ROL::ScalarMinimization<double>::run(double &, double &, int &, int &, class ROL::ScalarFunction<double> &, const double, const double, class ROL::ScalarMinimizationStatusTest<double> &) const --> void", pybind11::arg("fx"), pybind11::arg("x"), pybind11::arg("nfval"), pybind11::arg("ngrad"), pybind11::arg("f"), pybind11::arg("A"), pybind11::arg("B"), pybind11::arg("test"));
		cl.def("assign", (class ROL::ScalarMinimization<double> & (ROL::ScalarMinimization<double>::*)(const class ROL::ScalarMinimization<double> &)) &ROL::ScalarMinimization<double>::operator=, "C++: ROL::ScalarMinimization<double>::operator=(const class ROL::ScalarMinimization<double> &) --> class ROL::ScalarMinimization<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::GoldenSectionScalarMinimization file:ROL_GoldenSectionScalarMinimization.hpp line:59
		pybind11::class_<ROL::GoldenSectionScalarMinimization<double>, std::shared_ptr<ROL::GoldenSectionScalarMinimization<double>>, PyCallBack_ROL_GoldenSectionScalarMinimization_double_t, ROL::ScalarMinimization<double>> cl(M("ROL"), "GoldenSectionScalarMinimization_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("parlist") );

		cl.def( pybind11::init( [](PyCallBack_ROL_GoldenSectionScalarMinimization_double_t const &o){ return new PyCallBack_ROL_GoldenSectionScalarMinimization_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::GoldenSectionScalarMinimization<double> const &o){ return new ROL::GoldenSectionScalarMinimization<double>(o); } ) );
		cl.def("run", (void (ROL::GoldenSectionScalarMinimization<double>::*)(double &, double &, int &, int &, class ROL::ScalarFunction<double> &, const double, const double, class ROL::ScalarMinimizationStatusTest<double> &) const) &ROL::GoldenSectionScalarMinimization<double>::run, "C++: ROL::GoldenSectionScalarMinimization<double>::run(double &, double &, int &, int &, class ROL::ScalarFunction<double> &, const double, const double, class ROL::ScalarMinimizationStatusTest<double> &) const --> void", pybind11::arg("fx"), pybind11::arg("x"), pybind11::arg("nfval"), pybind11::arg("ngrad"), pybind11::arg("f"), pybind11::arg("A"), pybind11::arg("B"), pybind11::arg("test"));
		cl.def("assign", (class ROL::GoldenSectionScalarMinimization<double> & (ROL::GoldenSectionScalarMinimization<double>::*)(const class ROL::GoldenSectionScalarMinimization<double> &)) &ROL::GoldenSectionScalarMinimization<double>::operator=, "C++: ROL::GoldenSectionScalarMinimization<double>::operator=(const class ROL::GoldenSectionScalarMinimization<double> &) --> class ROL::GoldenSectionScalarMinimization<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("run", (void (ROL::ScalarMinimization<double>::*)(double &, double &, int &, int &, class ROL::ScalarFunction<double> &, const double, const double) const) &ROL::ScalarMinimization<double>::run, "C++: ROL::ScalarMinimization<double>::run(double &, double &, int &, int &, class ROL::ScalarFunction<double> &, const double, const double) const --> void", pybind11::arg("fx"), pybind11::arg("x"), pybind11::arg("nfval"), pybind11::arg("ngrad"), pybind11::arg("f"), pybind11::arg("A"), pybind11::arg("B"));
		cl.def("run", (void (ROL::ScalarMinimization<double>::*)(double &, double &, int &, int &, class ROL::ScalarFunction<double> &, const double, const double, class ROL::ScalarMinimizationStatusTest<double> &) const) &ROL::ScalarMinimization<double>::run, "C++: ROL::ScalarMinimization<double>::run(double &, double &, int &, int &, class ROL::ScalarFunction<double> &, const double, const double, class ROL::ScalarMinimizationStatusTest<double> &) const --> void", pybind11::arg("fx"), pybind11::arg("x"), pybind11::arg("nfval"), pybind11::arg("ngrad"), pybind11::arg("f"), pybind11::arg("A"), pybind11::arg("B"), pybind11::arg("test"));
		cl.def("assign", (class ROL::ScalarMinimization<double> & (ROL::ScalarMinimization<double>::*)(const class ROL::ScalarMinimization<double> &)) &ROL::ScalarMinimization<double>::operator=, "C++: ROL::ScalarMinimization<double>::operator=(const class ROL::ScalarMinimization<double> &) --> class ROL::ScalarMinimization<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::Bracketing file:ROL_Bracketing.hpp line:57
		pybind11::class_<ROL::Bracketing<double>, std::shared_ptr<ROL::Bracketing<double>>, PyCallBack_ROL_Bracketing_double_t> cl(M("ROL"), "Bracketing_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new ROL::Bracketing<double>(); }, [](){ return new PyCallBack_ROL_Bracketing_double_t(); } ) );
		cl.def( pybind11::init( [](PyCallBack_ROL_Bracketing_double_t const &o){ return new PyCallBack_ROL_Bracketing_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::Bracketing<double> const &o){ return new ROL::Bracketing<double>(o); } ) );
		cl.def("run", (void (ROL::Bracketing<double>::*)(double &, double &, double &, double &, double &, double &, int &, int &, class ROL::ScalarFunction<double> &) const) &ROL::Bracketing<double>::run, "C++: ROL::Bracketing<double>::run(double &, double &, double &, double &, double &, double &, int &, int &, class ROL::ScalarFunction<double> &) const --> void", pybind11::arg("x"), pybind11::arg("fx"), pybind11::arg("a"), pybind11::arg("fa"), pybind11::arg("b"), pybind11::arg("fb"), pybind11::arg("nfval"), pybind11::arg("ngrad"), pybind11::arg("f"));
		cl.def("run", (void (ROL::Bracketing<double>::*)(double &, double &, double &, double &, double &, double &, int &, int &, class ROL::ScalarFunction<double> &, class ROL::ScalarMinimizationStatusTest<double> &) const) &ROL::Bracketing<double>::run, "C++: ROL::Bracketing<double>::run(double &, double &, double &, double &, double &, double &, int &, int &, class ROL::ScalarFunction<double> &, class ROL::ScalarMinimizationStatusTest<double> &) const --> void", pybind11::arg("x"), pybind11::arg("fx"), pybind11::arg("a"), pybind11::arg("fa"), pybind11::arg("b"), pybind11::arg("fb"), pybind11::arg("nfval"), pybind11::arg("ngrad"), pybind11::arg("f"), pybind11::arg("test"));
		cl.def("assign", (class ROL::Bracketing<double> & (ROL::Bracketing<double>::*)(const class ROL::Bracketing<double> &)) &ROL::Bracketing<double>::operator=, "C++: ROL::Bracketing<double>::operator=(const class ROL::Bracketing<double> &) --> class ROL::Bracketing<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::ScalarMinimizationLineSearch_U file:ROL_ScalarMinimizationLineSearch_U.hpp line:62
		pybind11::class_<ROL::ScalarMinimizationLineSearch_U<double>, std::shared_ptr<ROL::ScalarMinimizationLineSearch_U<double>>, PyCallBack_ROL_ScalarMinimizationLineSearch_U_double_t, ROL::LineSearch_U<double>> cl(M("ROL"), "ScalarMinimizationLineSearch_U_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](class Teuchos::ParameterList & a0){ return new ROL::ScalarMinimizationLineSearch_U<double>(a0); }, [](class Teuchos::ParameterList & a0){ return new PyCallBack_ROL_ScalarMinimizationLineSearch_U_double_t(a0); } ), "doc");
		cl.def( pybind11::init( [](class Teuchos::ParameterList & a0, const class std::shared_ptr<class ROL::ScalarMinimization<double> > & a1){ return new ROL::ScalarMinimizationLineSearch_U<double>(a0, a1); }, [](class Teuchos::ParameterList & a0, const class std::shared_ptr<class ROL::ScalarMinimization<double> > & a1){ return new PyCallBack_ROL_ScalarMinimizationLineSearch_U_double_t(a0, a1); } ), "doc");
		cl.def( pybind11::init( [](class Teuchos::ParameterList & a0, const class std::shared_ptr<class ROL::ScalarMinimization<double> > & a1, const class std::shared_ptr<class ROL::Bracketing<double> > & a2){ return new ROL::ScalarMinimizationLineSearch_U<double>(a0, a1, a2); }, [](class Teuchos::ParameterList & a0, const class std::shared_ptr<class ROL::ScalarMinimization<double> > & a1, const class std::shared_ptr<class ROL::Bracketing<double> > & a2){ return new PyCallBack_ROL_ScalarMinimizationLineSearch_U_double_t(a0, a1, a2); } ), "doc");
		cl.def( pybind11::init<class Teuchos::ParameterList &, const class std::shared_ptr<class ROL::ScalarMinimization<double> > &, const class std::shared_ptr<class ROL::Bracketing<double> > &, const class std::shared_ptr<class ROL::ScalarFunction<double> > &>(), pybind11::arg("parlist"), pybind11::arg("sm"), pybind11::arg("br"), pybind11::arg("sf") );

		cl.def( pybind11::init( [](PyCallBack_ROL_ScalarMinimizationLineSearch_U_double_t const &o){ return new PyCallBack_ROL_ScalarMinimizationLineSearch_U_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::ScalarMinimizationLineSearch_U<double> const &o){ return new ROL::ScalarMinimizationLineSearch_U<double>(o); } ) );
		cl.def("initialize", (void (ROL::ScalarMinimizationLineSearch_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::ScalarMinimizationLineSearch_U<double>::initialize, "C++: ROL::ScalarMinimizationLineSearch_U<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"));
		cl.def("run", (void (ROL::ScalarMinimizationLineSearch_U<double>::*)(double &, double &, int &, int &, const double &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &)) &ROL::ScalarMinimizationLineSearch_U<double>::run, "C++: ROL::ScalarMinimizationLineSearch_U<double>::run(double &, double &, int &, int &, const double &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &) --> void", pybind11::arg("alpha"), pybind11::arg("fval"), pybind11::arg("ls_neval"), pybind11::arg("ls_ngrad"), pybind11::arg("gs"), pybind11::arg("s"), pybind11::arg("x"), pybind11::arg("obj"));
		cl.def("assign", (class ROL::ScalarMinimizationLineSearch_U<double> & (ROL::ScalarMinimizationLineSearch_U<double>::*)(const class ROL::ScalarMinimizationLineSearch_U<double> &)) &ROL::ScalarMinimizationLineSearch_U<double>::operator=, "C++: ROL::ScalarMinimizationLineSearch_U<double>::operator=(const class ROL::ScalarMinimizationLineSearch_U<double> &) --> class ROL::ScalarMinimizationLineSearch_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("initialize", (void (ROL::LineSearch_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::LineSearch_U<double>::initialize, "C++: ROL::LineSearch_U<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"));
		cl.def("run", (void (ROL::LineSearch_U<double>::*)(double &, double &, int &, int &, const double &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &)) &ROL::LineSearch_U<double>::run, "C++: ROL::LineSearch_U<double>::run(double &, double &, int &, int &, const double &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &) --> void", pybind11::arg("alpha"), pybind11::arg("fval"), pybind11::arg("ls_neval"), pybind11::arg("ls_ngrad"), pybind11::arg("gs"), pybind11::arg("s"), pybind11::arg("x"), pybind11::arg("obj"));
		cl.def("setMaxitUpdate", (void (ROL::LineSearch_U<double>::*)(double &, double &, const double &)) &ROL::LineSearch_U<double>::setMaxitUpdate, "C++: ROL::LineSearch_U<double>::setMaxitUpdate(double &, double &, const double &) --> void", pybind11::arg("alpha"), pybind11::arg("fnew"), pybind11::arg("fold"));
		cl.def("assign", (class ROL::LineSearch_U<double> & (ROL::LineSearch_U<double>::*)(const class ROL::LineSearch_U<double> &)) &ROL::LineSearch_U<double>::operator=, "C++: ROL::LineSearch_U<double>::operator=(const class ROL::LineSearch_U<double> &) --> class ROL::LineSearch_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	// ROL::LineSearchUFactory(class Teuchos::ParameterList &) file:ROL_LineSearch_U_Factory.hpp line:55
	M("ROL").def("LineSearchUFactory", (class std::shared_ptr<class ROL::LineSearch_U<double> > (*)(class Teuchos::ParameterList &)) &ROL::LineSearchUFactory<double>, "C++: ROL::LineSearchUFactory(class Teuchos::ParameterList &) --> class std::shared_ptr<class ROL::LineSearch_U<double> >", pybind11::arg("parlist"));

	{ // ROL::DescentDirection_U file:ROL_DescentDirection_U.hpp line:58
		pybind11::class_<ROL::DescentDirection_U<double>, std::shared_ptr<ROL::DescentDirection_U<double>>, PyCallBack_ROL_DescentDirection_U_double_t> cl(M("ROL"), "DescentDirection_U_double_t", "", pybind11::module_local());
		cl.def(pybind11::init<PyCallBack_ROL_DescentDirection_U_double_t const &>());
		cl.def( pybind11::init( [](){ return new PyCallBack_ROL_DescentDirection_U_double_t(); } ) );
		cl.def("initialize", (void (ROL::DescentDirection_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::DescentDirection_U<double>::initialize, "C++: ROL::DescentDirection_U<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"));
		cl.def("compute", (void (ROL::DescentDirection_U<double>::*)(class ROL::Vector<double> &, double &, double &, int &, int &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &)) &ROL::DescentDirection_U<double>::compute, "C++: ROL::DescentDirection_U<double>::compute(class ROL::Vector<double> &, double &, double &, int &, int &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, class ROL::Objective<double> &) --> void", pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("sdotg"), pybind11::arg("iter"), pybind11::arg("flag"), pybind11::arg("x"), pybind11::arg("g"), pybind11::arg("obj"));
		cl.def("update", (void (ROL::DescentDirection_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const double, const int)) &ROL::DescentDirection_U<double>::update, "C++: ROL::DescentDirection_U<double>::update(const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const class ROL::Vector<double> &, const double, const int) --> void", pybind11::arg("x"), pybind11::arg("s"), pybind11::arg("gold"), pybind11::arg("gnew"), pybind11::arg("snorm"), pybind11::arg("iter"));
		cl.def("printName", (std::string (ROL::DescentDirection_U<double>::*)() const) &ROL::DescentDirection_U<double>::printName, "C++: ROL::DescentDirection_U<double>::printName() const --> std::string");
		cl.def("assign", (class ROL::DescentDirection_U<double> & (ROL::DescentDirection_U<double>::*)(const class ROL::DescentDirection_U<double> &)) &ROL::DescentDirection_U<double>::operator=, "C++: ROL::DescentDirection_U<double>::operator=(const class ROL::DescentDirection_U<double> &) --> class ROL::DescentDirection_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
}
