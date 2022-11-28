#include <ROL_CauchyPoint_U.hpp>
#include <ROL_DogLeg_U.hpp>
#include <ROL_DoubleDogLeg_U.hpp>
#include <ROL_Elementwise_Function.hpp>
#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_SPGTrustRegion_U.hpp>
#include <ROL_Secant.hpp>
#include <ROL_TruncatedCG_U.hpp>
#include <ROL_TrustRegionModel_U.hpp>
#include <ROL_TrustRegion_U.hpp>
#include <ROL_TrustRegion_U_Factory.hpp>
#include <ROL_TrustRegion_U_Types.hpp>
#include <ROL_UpdateType.hpp>
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

// ROL::CauchyPoint_U file:ROL_CauchyPoint_U.hpp line:57
struct PyCallBack_ROL_CauchyPoint_U_double_t : public ROL::CauchyPoint_U<double> {
	using ROL::CauchyPoint_U<double>::CauchyPoint_U;

	void initialize(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::CauchyPoint_U<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return CauchyPoint_U::initialize(a0, a1);
	}
	void solve(class ROL::Vector<double> & a0, double & a1, double & a2, int & a3, int & a4, const double a5, class ROL::TrustRegionModel_U<double> & a6) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::CauchyPoint_U<double> *>(this), "solve");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return CauchyPoint_U::solve(a0, a1, a2, a3, a4, a5, a6);
	}
};

// ROL::DogLeg_U file:ROL_DogLeg_U.hpp line:57
struct PyCallBack_ROL_DogLeg_U_double_t : public ROL::DogLeg_U<double> {
	using ROL::DogLeg_U<double>::DogLeg_U;

	void initialize(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::DogLeg_U<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return DogLeg_U::initialize(a0, a1);
	}
	void solve(class ROL::Vector<double> & a0, double & a1, double & a2, int & a3, int & a4, const double a5, class ROL::TrustRegionModel_U<double> & a6) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::DogLeg_U<double> *>(this), "solve");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return DogLeg_U::solve(a0, a1, a2, a3, a4, a5, a6);
	}
};

// ROL::DoubleDogLeg_U file:ROL_DoubleDogLeg_U.hpp line:57
struct PyCallBack_ROL_DoubleDogLeg_U_double_t : public ROL::DoubleDogLeg_U<double> {
	using ROL::DoubleDogLeg_U<double>::DoubleDogLeg_U;

	void initialize(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::DoubleDogLeg_U<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return DoubleDogLeg_U::initialize(a0, a1);
	}
	void solve(class ROL::Vector<double> & a0, double & a1, double & a2, int & a3, int & a4, const double a5, class ROL::TrustRegionModel_U<double> & a6) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::DoubleDogLeg_U<double> *>(this), "solve");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return DoubleDogLeg_U::solve(a0, a1, a2, a3, a4, a5, a6);
	}
};

// ROL::TruncatedCG_U file:ROL_TruncatedCG_U.hpp line:57
struct PyCallBack_ROL_TruncatedCG_U_double_t : public ROL::TruncatedCG_U<double> {
	using ROL::TruncatedCG_U<double>::TruncatedCG_U;

	void initialize(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TruncatedCG_U<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TruncatedCG_U::initialize(a0, a1);
	}
	void solve(class ROL::Vector<double> & a0, double & a1, double & a2, int & a3, int & a4, const double a5, class ROL::TrustRegionModel_U<double> & a6) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::TruncatedCG_U<double> *>(this), "solve");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return TruncatedCG_U::solve(a0, a1, a2, a3, a4, a5, a6);
	}
};

// ROL::SPGTrustRegion_U file:ROL_SPGTrustRegion_U.hpp line:59
struct PyCallBack_ROL_SPGTrustRegion_U_double_t : public ROL::SPGTrustRegion_U<double> {
	using ROL::SPGTrustRegion_U<double>::SPGTrustRegion_U;

	void initialize(const class ROL::Vector<double> & a0, const class ROL::Vector<double> & a1) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::SPGTrustRegion_U<double> *>(this), "initialize");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SPGTrustRegion_U::initialize(a0, a1);
	}
	void solve(class ROL::Vector<double> & a0, double & a1, double & a2, int & a3, int & a4, const double a5, class ROL::TrustRegionModel_U<double> & a6) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::SPGTrustRegion_U<double> *>(this), "solve");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SPGTrustRegion_U::solve(a0, a1, a2, a3, a4, a5, a6);
	}
};

void bind_ROL_CauchyPoint_U(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // ROL::CauchyPoint_U file:ROL_CauchyPoint_U.hpp line:57
		pybind11::class_<ROL::CauchyPoint_U<double>, Teuchos::RCP<ROL::CauchyPoint_U<double>>, PyCallBack_ROL_CauchyPoint_U_double_t, ROL::TrustRegion_U<double>> cl(M("ROL"), "CauchyPoint_U_double_t", "");
		cl.def( pybind11::init( [](){ return new ROL::CauchyPoint_U<double>(); }, [](){ return new PyCallBack_ROL_CauchyPoint_U_double_t(); } ) );
		cl.def( pybind11::init( [](PyCallBack_ROL_CauchyPoint_U_double_t const &o){ return new PyCallBack_ROL_CauchyPoint_U_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::CauchyPoint_U<double> const &o){ return new ROL::CauchyPoint_U<double>(o); } ) );
		cl.def("initialize", (void (ROL::CauchyPoint_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::CauchyPoint_U<double>::initialize, "C++: ROL::CauchyPoint_U<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"));
		cl.def("solve", (void (ROL::CauchyPoint_U<double>::*)(class ROL::Vector<double> &, double &, double &, int &, int &, const double, class ROL::TrustRegionModel_U<double> &)) &ROL::CauchyPoint_U<double>::solve, "C++: ROL::CauchyPoint_U<double>::solve(class ROL::Vector<double> &, double &, double &, int &, int &, const double, class ROL::TrustRegionModel_U<double> &) --> void", pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("pRed"), pybind11::arg("iflag"), pybind11::arg("iter"), pybind11::arg("del"), pybind11::arg("model"));
		cl.def("assign", (class ROL::CauchyPoint_U<double> & (ROL::CauchyPoint_U<double>::*)(const class ROL::CauchyPoint_U<double> &)) &ROL::CauchyPoint_U<double>::operator=, "C++: ROL::CauchyPoint_U<double>::operator=(const class ROL::CauchyPoint_U<double> &) --> class ROL::CauchyPoint_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("initialize", (void (ROL::TrustRegion_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::TrustRegion_U<double>::initialize, "C++: ROL::TrustRegion_U<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"));
		cl.def("solve", (void (ROL::TrustRegion_U<double>::*)(class ROL::Vector<double> &, double &, double &, int &, int &, const double, class ROL::TrustRegionModel_U<double> &)) &ROL::TrustRegion_U<double>::solve, "C++: ROL::TrustRegion_U<double>::solve(class ROL::Vector<double> &, double &, double &, int &, int &, const double, class ROL::TrustRegionModel_U<double> &) --> void", pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("pRed"), pybind11::arg("iflag"), pybind11::arg("iter"), pybind11::arg("del"), pybind11::arg("model"));
		cl.def("assign", (class ROL::TrustRegion_U<double> & (ROL::TrustRegion_U<double>::*)(const class ROL::TrustRegion_U<double> &)) &ROL::TrustRegion_U<double>::operator=, "C++: ROL::TrustRegion_U<double>::operator=(const class ROL::TrustRegion_U<double> &) --> class ROL::TrustRegion_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::DogLeg_U file:ROL_DogLeg_U.hpp line:57
		pybind11::class_<ROL::DogLeg_U<double>, Teuchos::RCP<ROL::DogLeg_U<double>>, PyCallBack_ROL_DogLeg_U_double_t, ROL::TrustRegion_U<double>> cl(M("ROL"), "DogLeg_U_double_t", "");
		cl.def( pybind11::init( [](){ return new ROL::DogLeg_U<double>(); }, [](){ return new PyCallBack_ROL_DogLeg_U_double_t(); } ) );
		cl.def( pybind11::init( [](PyCallBack_ROL_DogLeg_U_double_t const &o){ return new PyCallBack_ROL_DogLeg_U_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::DogLeg_U<double> const &o){ return new ROL::DogLeg_U<double>(o); } ) );
		cl.def("initialize", (void (ROL::DogLeg_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::DogLeg_U<double>::initialize, "C++: ROL::DogLeg_U<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"));
		cl.def("solve", (void (ROL::DogLeg_U<double>::*)(class ROL::Vector<double> &, double &, double &, int &, int &, const double, class ROL::TrustRegionModel_U<double> &)) &ROL::DogLeg_U<double>::solve, "C++: ROL::DogLeg_U<double>::solve(class ROL::Vector<double> &, double &, double &, int &, int &, const double, class ROL::TrustRegionModel_U<double> &) --> void", pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("pRed"), pybind11::arg("iflag"), pybind11::arg("iter"), pybind11::arg("del"), pybind11::arg("model"));
		cl.def("assign", (class ROL::DogLeg_U<double> & (ROL::DogLeg_U<double>::*)(const class ROL::DogLeg_U<double> &)) &ROL::DogLeg_U<double>::operator=, "C++: ROL::DogLeg_U<double>::operator=(const class ROL::DogLeg_U<double> &) --> class ROL::DogLeg_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("initialize", (void (ROL::TrustRegion_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::TrustRegion_U<double>::initialize, "C++: ROL::TrustRegion_U<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"));
		cl.def("solve", (void (ROL::TrustRegion_U<double>::*)(class ROL::Vector<double> &, double &, double &, int &, int &, const double, class ROL::TrustRegionModel_U<double> &)) &ROL::TrustRegion_U<double>::solve, "C++: ROL::TrustRegion_U<double>::solve(class ROL::Vector<double> &, double &, double &, int &, int &, const double, class ROL::TrustRegionModel_U<double> &) --> void", pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("pRed"), pybind11::arg("iflag"), pybind11::arg("iter"), pybind11::arg("del"), pybind11::arg("model"));
		cl.def("assign", (class ROL::TrustRegion_U<double> & (ROL::TrustRegion_U<double>::*)(const class ROL::TrustRegion_U<double> &)) &ROL::TrustRegion_U<double>::operator=, "C++: ROL::TrustRegion_U<double>::operator=(const class ROL::TrustRegion_U<double> &) --> class ROL::TrustRegion_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::DoubleDogLeg_U file:ROL_DoubleDogLeg_U.hpp line:57
		pybind11::class_<ROL::DoubleDogLeg_U<double>, Teuchos::RCP<ROL::DoubleDogLeg_U<double>>, PyCallBack_ROL_DoubleDogLeg_U_double_t, ROL::TrustRegion_U<double>> cl(M("ROL"), "DoubleDogLeg_U_double_t", "");
		cl.def( pybind11::init( [](){ return new ROL::DoubleDogLeg_U<double>(); }, [](){ return new PyCallBack_ROL_DoubleDogLeg_U_double_t(); } ) );
		cl.def( pybind11::init( [](PyCallBack_ROL_DoubleDogLeg_U_double_t const &o){ return new PyCallBack_ROL_DoubleDogLeg_U_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::DoubleDogLeg_U<double> const &o){ return new ROL::DoubleDogLeg_U<double>(o); } ) );
		cl.def("initialize", (void (ROL::DoubleDogLeg_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::DoubleDogLeg_U<double>::initialize, "C++: ROL::DoubleDogLeg_U<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"));
		cl.def("solve", (void (ROL::DoubleDogLeg_U<double>::*)(class ROL::Vector<double> &, double &, double &, int &, int &, const double, class ROL::TrustRegionModel_U<double> &)) &ROL::DoubleDogLeg_U<double>::solve, "C++: ROL::DoubleDogLeg_U<double>::solve(class ROL::Vector<double> &, double &, double &, int &, int &, const double, class ROL::TrustRegionModel_U<double> &) --> void", pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("pRed"), pybind11::arg("iflag"), pybind11::arg("iter"), pybind11::arg("del"), pybind11::arg("model"));
		cl.def("assign", (class ROL::DoubleDogLeg_U<double> & (ROL::DoubleDogLeg_U<double>::*)(const class ROL::DoubleDogLeg_U<double> &)) &ROL::DoubleDogLeg_U<double>::operator=, "C++: ROL::DoubleDogLeg_U<double>::operator=(const class ROL::DoubleDogLeg_U<double> &) --> class ROL::DoubleDogLeg_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("initialize", (void (ROL::TrustRegion_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::TrustRegion_U<double>::initialize, "C++: ROL::TrustRegion_U<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"));
		cl.def("solve", (void (ROL::TrustRegion_U<double>::*)(class ROL::Vector<double> &, double &, double &, int &, int &, const double, class ROL::TrustRegionModel_U<double> &)) &ROL::TrustRegion_U<double>::solve, "C++: ROL::TrustRegion_U<double>::solve(class ROL::Vector<double> &, double &, double &, int &, int &, const double, class ROL::TrustRegionModel_U<double> &) --> void", pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("pRed"), pybind11::arg("iflag"), pybind11::arg("iter"), pybind11::arg("del"), pybind11::arg("model"));
		cl.def("assign", (class ROL::TrustRegion_U<double> & (ROL::TrustRegion_U<double>::*)(const class ROL::TrustRegion_U<double> &)) &ROL::TrustRegion_U<double>::operator=, "C++: ROL::TrustRegion_U<double>::operator=(const class ROL::TrustRegion_U<double> &) --> class ROL::TrustRegion_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::TruncatedCG_U file:ROL_TruncatedCG_U.hpp line:57
		pybind11::class_<ROL::TruncatedCG_U<double>, Teuchos::RCP<ROL::TruncatedCG_U<double>>, PyCallBack_ROL_TruncatedCG_U_double_t, ROL::TrustRegion_U<double>> cl(M("ROL"), "TruncatedCG_U_double_t", "");
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("parlist") );

		cl.def( pybind11::init( [](PyCallBack_ROL_TruncatedCG_U_double_t const &o){ return new PyCallBack_ROL_TruncatedCG_U_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::TruncatedCG_U<double> const &o){ return new ROL::TruncatedCG_U<double>(o); } ) );
		cl.def("initialize", (void (ROL::TruncatedCG_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::TruncatedCG_U<double>::initialize, "C++: ROL::TruncatedCG_U<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"));
		cl.def("solve", (void (ROL::TruncatedCG_U<double>::*)(class ROL::Vector<double> &, double &, double &, int &, int &, const double, class ROL::TrustRegionModel_U<double> &)) &ROL::TruncatedCG_U<double>::solve, "C++: ROL::TruncatedCG_U<double>::solve(class ROL::Vector<double> &, double &, double &, int &, int &, const double, class ROL::TrustRegionModel_U<double> &) --> void", pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("pRed"), pybind11::arg("iflag"), pybind11::arg("iter"), pybind11::arg("del"), pybind11::arg("model"));
		cl.def("assign", (class ROL::TruncatedCG_U<double> & (ROL::TruncatedCG_U<double>::*)(const class ROL::TruncatedCG_U<double> &)) &ROL::TruncatedCG_U<double>::operator=, "C++: ROL::TruncatedCG_U<double>::operator=(const class ROL::TruncatedCG_U<double> &) --> class ROL::TruncatedCG_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("initialize", (void (ROL::TrustRegion_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::TrustRegion_U<double>::initialize, "C++: ROL::TrustRegion_U<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"));
		cl.def("solve", (void (ROL::TrustRegion_U<double>::*)(class ROL::Vector<double> &, double &, double &, int &, int &, const double, class ROL::TrustRegionModel_U<double> &)) &ROL::TrustRegion_U<double>::solve, "C++: ROL::TrustRegion_U<double>::solve(class ROL::Vector<double> &, double &, double &, int &, int &, const double, class ROL::TrustRegionModel_U<double> &) --> void", pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("pRed"), pybind11::arg("iflag"), pybind11::arg("iter"), pybind11::arg("del"), pybind11::arg("model"));
		cl.def("assign", (class ROL::TrustRegion_U<double> & (ROL::TrustRegion_U<double>::*)(const class ROL::TrustRegion_U<double> &)) &ROL::TrustRegion_U<double>::operator=, "C++: ROL::TrustRegion_U<double>::operator=(const class ROL::TrustRegion_U<double> &) --> class ROL::TrustRegion_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::SPGTrustRegion_U file:ROL_SPGTrustRegion_U.hpp line:59
		pybind11::class_<ROL::SPGTrustRegion_U<double>, Teuchos::RCP<ROL::SPGTrustRegion_U<double>>, PyCallBack_ROL_SPGTrustRegion_U_double_t, ROL::TrustRegion_U<double>> cl(M("ROL"), "SPGTrustRegion_U_double_t", "");
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("parlist") );

		cl.def( pybind11::init( [](PyCallBack_ROL_SPGTrustRegion_U_double_t const &o){ return new PyCallBack_ROL_SPGTrustRegion_U_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::SPGTrustRegion_U<double> const &o){ return new ROL::SPGTrustRegion_U<double>(o); } ) );
		cl.def("initialize", (void (ROL::SPGTrustRegion_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::SPGTrustRegion_U<double>::initialize, "C++: ROL::SPGTrustRegion_U<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"));
		cl.def("solve", (void (ROL::SPGTrustRegion_U<double>::*)(class ROL::Vector<double> &, double &, double &, int &, int &, const double, class ROL::TrustRegionModel_U<double> &)) &ROL::SPGTrustRegion_U<double>::solve, "C++: ROL::SPGTrustRegion_U<double>::solve(class ROL::Vector<double> &, double &, double &, int &, int &, const double, class ROL::TrustRegionModel_U<double> &) --> void", pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("pRed"), pybind11::arg("iflag"), pybind11::arg("iter"), pybind11::arg("del"), pybind11::arg("model"));
		cl.def("assign", (class ROL::SPGTrustRegion_U<double> & (ROL::SPGTrustRegion_U<double>::*)(const class ROL::SPGTrustRegion_U<double> &)) &ROL::SPGTrustRegion_U<double>::operator=, "C++: ROL::SPGTrustRegion_U<double>::operator=(const class ROL::SPGTrustRegion_U<double> &) --> class ROL::SPGTrustRegion_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("initialize", (void (ROL::TrustRegion_U<double>::*)(const class ROL::Vector<double> &, const class ROL::Vector<double> &)) &ROL::TrustRegion_U<double>::initialize, "C++: ROL::TrustRegion_U<double>::initialize(const class ROL::Vector<double> &, const class ROL::Vector<double> &) --> void", pybind11::arg("x"), pybind11::arg("g"));
		cl.def("solve", (void (ROL::TrustRegion_U<double>::*)(class ROL::Vector<double> &, double &, double &, int &, int &, const double, class ROL::TrustRegionModel_U<double> &)) &ROL::TrustRegion_U<double>::solve, "C++: ROL::TrustRegion_U<double>::solve(class ROL::Vector<double> &, double &, double &, int &, int &, const double, class ROL::TrustRegionModel_U<double> &) --> void", pybind11::arg("s"), pybind11::arg("snorm"), pybind11::arg("pRed"), pybind11::arg("iflag"), pybind11::arg("iter"), pybind11::arg("del"), pybind11::arg("model"));
		cl.def("assign", (class ROL::TrustRegion_U<double> & (ROL::TrustRegion_U<double>::*)(const class ROL::TrustRegion_U<double> &)) &ROL::TrustRegion_U<double>::operator=, "C++: ROL::TrustRegion_U<double>::operator=(const class ROL::TrustRegion_U<double> &) --> class ROL::TrustRegion_U<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	// ROL::TrustRegionUFactory(class Teuchos::ParameterList &) file:ROL_TrustRegion_U_Factory.hpp line:61
	M("ROL").def("TrustRegionUFactory", (class Teuchos::RCP<class ROL::TrustRegion_U<double> > (*)(class Teuchos::ParameterList &)) &ROL::TrustRegionUFactory<double>, "C++: ROL::TrustRegionUFactory(class Teuchos::ParameterList &) --> class Teuchos::RCP<class ROL::TrustRegion_U<double> >", pybind11::arg("list"));

}
