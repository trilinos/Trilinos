#include <ROL_DescentDirection_U.hpp>
#include <ROL_Elementwise_Function.hpp>
#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_LineSearch_U.hpp>
#include <ROL_Objective.hpp>
#include <ROL_TypeU_LineSearchAlgorithm.hpp>
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
#include <cwchar>
#include <deque>
#include <ios>
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

// ROL::TypeU::LineSearchAlgorithm file:ROL_TypeU_LineSearchAlgorithm.hpp line:100
struct PyCallBack_ROL_TypeU_LineSearchAlgorithm_double_t : public ROL::TypeU::LineSearchAlgorithm<double> {
	using ROL::TypeU::LineSearchAlgorithm<double>::LineSearchAlgorithm;

};

void bind_ROL_TypeU_LineSearchAlgorithm(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // ROL::TypeU::LineSearchAlgorithm file:ROL_TypeU_LineSearchAlgorithm.hpp line:100
		pybind11::class_<ROL::TypeU::LineSearchAlgorithm<double>, Teuchos::RCP<ROL::TypeU::LineSearchAlgorithm<double>>, PyCallBack_ROL_TypeU_LineSearchAlgorithm_double_t, ROL::TypeU::Algorithm<double>> cl(M("ROL::TypeU"), "LineSearchAlgorithm_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](class Teuchos::ParameterList & a0){ return new ROL::TypeU::LineSearchAlgorithm<double>(a0); }, [](class Teuchos::ParameterList & a0){ return new PyCallBack_ROL_TypeU_LineSearchAlgorithm_double_t(a0); } ), "doc");
		cl.def( pybind11::init( [](class Teuchos::ParameterList & a0, const class Teuchos::RCP<class ROL::DescentDirection_U<double> > & a1){ return new ROL::TypeU::LineSearchAlgorithm<double>(a0, a1); }, [](class Teuchos::ParameterList & a0, const class Teuchos::RCP<class ROL::DescentDirection_U<double> > & a1){ return new PyCallBack_ROL_TypeU_LineSearchAlgorithm_double_t(a0, a1); } ), "doc");
		cl.def( pybind11::init<class Teuchos::ParameterList &, const class Teuchos::RCP<class ROL::DescentDirection_U<double> > &, const class Teuchos::RCP<class ROL::LineSearch_U<double> > &>(), pybind11::arg("parlist"), pybind11::arg("descent"), pybind11::arg("lineSearch") );

		cl.def( pybind11::init( [](PyCallBack_ROL_TypeU_LineSearchAlgorithm_double_t const &o){ return new PyCallBack_ROL_TypeU_LineSearchAlgorithm_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::TypeU::LineSearchAlgorithm<double> const &o){ return new ROL::TypeU::LineSearchAlgorithm<double>(o); } ) );
		cl.def("setStatusTest", [](ROL::TypeU::Algorithm<double> &o, const class Teuchos::RCP<class ROL::StatusTest<double> > & a0) -> void { return o.setStatusTest(a0); }, "", pybind11::arg("status"));
		cl.def("setStatusTest", (void (ROL::TypeU::Algorithm<double>::*)(const class Teuchos::RCP<class ROL::StatusTest<double> > &, bool)) &ROL::TypeU::Algorithm<double>::setStatusTest, "C++: ROL::TypeU::Algorithm<double>::setStatusTest(const class Teuchos::RCP<class ROL::StatusTest<double> > &, bool) --> void", pybind11::arg("status"), pybind11::arg("combineStatus"));
		cl.def("getState", (class Teuchos::RCP<const struct ROL::TypeU::AlgorithmState<double> > (ROL::TypeU::Algorithm<double>::*)() const) &ROL::TypeU::Algorithm<double>::getState, "C++: ROL::TypeU::Algorithm<double>::getState() const --> class Teuchos::RCP<const struct ROL::TypeU::AlgorithmState<double> >");
		cl.def("reset", (void (ROL::TypeU::Algorithm<double>::*)()) &ROL::TypeU::Algorithm<double>::reset, "C++: ROL::TypeU::Algorithm<double>::reset() --> void");
	}
}
