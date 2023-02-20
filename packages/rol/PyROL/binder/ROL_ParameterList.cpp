#include <ROL_CombinedStatusTest.hpp>
#include <ROL_ParameterList.hpp>
#include <ROL_StatusTest.hpp>
#include <ROL_Types.hpp>
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

#include <functional>
#include <pybind11/pybind11.h>
#include <string>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>


#ifndef BINDER_PYBIND11_TYPE_CASTER
	#define BINDER_PYBIND11_TYPE_CASTER
	PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)
	PYBIND11_DECLARE_HOLDER_TYPE(T, T*)
	PYBIND11_MAKE_OPAQUE(std::shared_ptr<void>)
#endif

// ROL::StatusTest file:ROL_StatusTest.hpp line:58
struct PyCallBack_ROL_StatusTest_double_t : public ROL::StatusTest<double> {
	using ROL::StatusTest<double>::StatusTest;

	bool check(struct ROL::AlgorithmState<double> & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::StatusTest<double> *>(this), "check");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return StatusTest::check(a0);
	}
};

// ROL::CombinedStatusTest file:ROL_CombinedStatusTest.hpp line:57
struct PyCallBack_ROL_CombinedStatusTest_double_t : public ROL::CombinedStatusTest<double> {
	using ROL::CombinedStatusTest<double>::CombinedStatusTest;

	bool check(struct ROL::AlgorithmState<double> & a0) override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const ROL::CombinedStatusTest<double> *>(this), "check");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return CombinedStatusTest::check(a0);
	}
};

void bind_ROL_ParameterList(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	// ROL::readParametersFromXml(const std::string &, class Teuchos::ParameterList &) file:ROL_ParameterList.hpp line:53
	M("ROL").def("readParametersFromXml", (void (*)(const std::string &, class Teuchos::ParameterList &)) &ROL::readParametersFromXml, "C++: ROL::readParametersFromXml(const std::string &, class Teuchos::ParameterList &) --> void", pybind11::arg("filename"), pybind11::arg("parlist"));

	// ROL::writeParameterListToXmlFile(class Teuchos::ParameterList &, const std::string &) file:ROL_ParameterList.hpp line:59
	M("ROL").def("writeParameterListToXmlFile", (void (*)(class Teuchos::ParameterList &, const std::string &)) &ROL::writeParameterListToXmlFile, "C++: ROL::writeParameterListToXmlFile(class Teuchos::ParameterList &, const std::string &) --> void", pybind11::arg("parlist"), pybind11::arg("filename"));

	// ROL::updateParametersFromXmlFile(const std::string &, class Teuchos::ParameterList &) file:ROL_ParameterList.hpp line:64
	M("ROL").def("updateParametersFromXmlFile", (void (*)(const std::string &, class Teuchos::ParameterList &)) &ROL::updateParametersFromXmlFile, "C++: ROL::updateParametersFromXmlFile(const std::string &, class Teuchos::ParameterList &) --> void", pybind11::arg("filename"), pybind11::arg("parlist"));

	// ROL::getParametersFromXmlFile(const std::string &) file:ROL_ParameterList.hpp line:73
	M("ROL").def("getParametersFromXmlFile", (class std::shared_ptr<class Teuchos::ParameterList> (*)(const std::string &)) &ROL::getParametersFromXmlFile, "C++: ROL::getParametersFromXmlFile(const std::string &) --> class std::shared_ptr<class Teuchos::ParameterList>", pybind11::arg("filename"));

	{ // ROL::StatusTest file:ROL_StatusTest.hpp line:58
		pybind11::class_<ROL::StatusTest<double>, std::shared_ptr<ROL::StatusTest<double>>, PyCallBack_ROL_StatusTest_double_t> cl(M("ROL"), "StatusTest_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("parlist") );

		cl.def( pybind11::init( [](){ return new ROL::StatusTest<double>(); }, [](){ return new PyCallBack_ROL_StatusTest_double_t(); } ), "doc");
		cl.def( pybind11::init( [](double const & a0){ return new ROL::StatusTest<double>(a0); }, [](double const & a0){ return new PyCallBack_ROL_StatusTest_double_t(a0); } ), "doc");
		cl.def( pybind11::init( [](double const & a0, double const & a1){ return new ROL::StatusTest<double>(a0, a1); }, [](double const & a0, double const & a1){ return new PyCallBack_ROL_StatusTest_double_t(a0, a1); } ), "doc");
		cl.def( pybind11::init( [](double const & a0, double const & a1, int const & a2){ return new ROL::StatusTest<double>(a0, a1, a2); }, [](double const & a0, double const & a1, int const & a2){ return new PyCallBack_ROL_StatusTest_double_t(a0, a1, a2); } ), "doc");
		cl.def( pybind11::init<double, double, int, bool>(), pybind11::arg("gtol"), pybind11::arg("stol"), pybind11::arg("max_iter"), pybind11::arg("use_rel") );

		cl.def( pybind11::init( [](PyCallBack_ROL_StatusTest_double_t const &o){ return new PyCallBack_ROL_StatusTest_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::StatusTest<double> const &o){ return new ROL::StatusTest<double>(o); } ) );
		cl.def("check", (bool (ROL::StatusTest<double>::*)(struct ROL::AlgorithmState<double> &)) &ROL::StatusTest<double>::check, "C++: ROL::StatusTest<double>::check(struct ROL::AlgorithmState<double> &) --> bool", pybind11::arg("state"));
		cl.def("assign", (class ROL::StatusTest<double> & (ROL::StatusTest<double>::*)(const class ROL::StatusTest<double> &)) &ROL::StatusTest<double>::operator=, "C++: ROL::StatusTest<double>::operator=(const class ROL::StatusTest<double> &) --> class ROL::StatusTest<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // ROL::CombinedStatusTest file:ROL_CombinedStatusTest.hpp line:57
		pybind11::class_<ROL::CombinedStatusTest<double>, std::shared_ptr<ROL::CombinedStatusTest<double>>, PyCallBack_ROL_CombinedStatusTest_double_t, ROL::StatusTest<double>> cl(M("ROL"), "CombinedStatusTest_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new ROL::CombinedStatusTest<double>(); }, [](){ return new PyCallBack_ROL_CombinedStatusTest_double_t(); } ) );
		cl.def( pybind11::init( [](PyCallBack_ROL_CombinedStatusTest_double_t const &o){ return new PyCallBack_ROL_CombinedStatusTest_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::CombinedStatusTest<double> const &o){ return new ROL::CombinedStatusTest<double>(o); } ) );
		cl.def("reset", (void (ROL::CombinedStatusTest<double>::*)()) &ROL::CombinedStatusTest<double>::reset, "C++: ROL::CombinedStatusTest<double>::reset() --> void");
		cl.def("add", (void (ROL::CombinedStatusTest<double>::*)(const class std::shared_ptr<class ROL::StatusTest<double> > &)) &ROL::CombinedStatusTest<double>::add, "C++: ROL::CombinedStatusTest<double>::add(const class std::shared_ptr<class ROL::StatusTest<double> > &) --> void", pybind11::arg("status"));
		cl.def("check", (bool (ROL::CombinedStatusTest<double>::*)(struct ROL::AlgorithmState<double> &)) &ROL::CombinedStatusTest<double>::check, "C++: ROL::CombinedStatusTest<double>::check(struct ROL::AlgorithmState<double> &) --> bool", pybind11::arg("state"));
		cl.def("assign", (class ROL::CombinedStatusTest<double> & (ROL::CombinedStatusTest<double>::*)(const class ROL::CombinedStatusTest<double> &)) &ROL::CombinedStatusTest<double>::operator=, "C++: ROL::CombinedStatusTest<double>::operator=(const class ROL::CombinedStatusTest<double> &) --> class ROL::CombinedStatusTest<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("check", (bool (ROL::StatusTest<double>::*)(struct ROL::AlgorithmState<double> &)) &ROL::StatusTest<double>::check, "C++: ROL::StatusTest<double>::check(struct ROL::AlgorithmState<double> &) --> bool", pybind11::arg("state"));
		cl.def("assign", (class ROL::StatusTest<double> & (ROL::StatusTest<double>::*)(const class ROL::StatusTest<double> &)) &ROL::StatusTest<double>::operator=, "C++: ROL::StatusTest<double>::operator=(const class ROL::StatusTest<double> &) --> class ROL::StatusTest<double> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
}
