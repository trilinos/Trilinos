#include <ROL_BoundConstraint.hpp>
#include <ROL_Elementwise_Function.hpp>
#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_Objective.hpp>
#include <ROL_PartitionedVector.hpp>
#include <ROL_Secant.hpp>
#include <ROL_TypeB_QuasiNewtonAlgorithm.hpp>
#include <ROL_UpdateType.hpp>
#include <ROL_Vector.hpp>
#include <Teuchos_FilteredIterator.hpp>
#include <Teuchos_ParameterEntry.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ParameterListModifier.hpp>
#include <Teuchos_RCPDecl.hpp>
#include <Teuchos_StringIndexedOrderedValueObjectContainer.hpp>
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
#include <pybind11/smart_holder.h>
#include <string>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>


#ifndef BINDER_PYBIND11_TYPE_CASTER
	#define BINDER_PYBIND11_TYPE_CASTER
	PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)
	PYBIND11_DECLARE_HOLDER_TYPE(T, T*)
	PYBIND11_MAKE_OPAQUE(std::shared_ptr<void>)
#endif

// ROL::TypeB::QuasiNewtonAlgorithm file:ROL_TypeB_QuasiNewtonAlgorithm.hpp line:58
struct PyCallBack_ROL_TypeB_QuasiNewtonAlgorithm_double_t : public ROL::TypeB::QuasiNewtonAlgorithm<double> {
	using ROL::TypeB::QuasiNewtonAlgorithm<double>::QuasiNewtonAlgorithm;

};

void bind_ROL_TypeB_QuasiNewtonAlgorithm(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // ROL::TypeB::QuasiNewtonAlgorithm file:ROL_TypeB_QuasiNewtonAlgorithm.hpp line:58
		pybind11::class_<ROL::TypeB::QuasiNewtonAlgorithm<double>, std::shared_ptr<ROL::TypeB::QuasiNewtonAlgorithm<double>>, PyCallBack_ROL_TypeB_QuasiNewtonAlgorithm_double_t, ROL::TypeB::Algorithm<double>> cl(M("ROL::TypeB"), "QuasiNewtonAlgorithm_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](class Teuchos::ParameterList & a0){ return new ROL::TypeB::QuasiNewtonAlgorithm<double>(a0); }, [](class Teuchos::ParameterList & a0){ return new PyCallBack_ROL_TypeB_QuasiNewtonAlgorithm_double_t(a0); } ), "doc");
		cl.def( pybind11::init<class Teuchos::ParameterList &, const class std::shared_ptr<class ROL::Secant<double> > &>(), pybind11::arg("list"), pybind11::arg("secant") );

		cl.def( pybind11::init( [](PyCallBack_ROL_TypeB_QuasiNewtonAlgorithm_double_t const &o){ return new PyCallBack_ROL_TypeB_QuasiNewtonAlgorithm_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::TypeB::QuasiNewtonAlgorithm<double> const &o){ return new ROL::TypeB::QuasiNewtonAlgorithm<double>(o); } ) );
		cl.def("setStatusTest", [](ROL::TypeB::Algorithm<double> &o, const class std::shared_ptr<class ROL::StatusTest<double> > & a0) -> void { return o.setStatusTest(a0); }, "", pybind11::arg("status"));
		cl.def("setStatusTest", (void (ROL::TypeB::Algorithm<double>::*)(const class std::shared_ptr<class ROL::StatusTest<double> > &, const bool)) &ROL::TypeB::Algorithm<double>::setStatusTest, "C++: ROL::TypeB::Algorithm<double>::setStatusTest(const class std::shared_ptr<class ROL::StatusTest<double> > &, const bool) --> void", pybind11::arg("status"), pybind11::arg("combineStatus"));
		cl.def("getState", (class std::shared_ptr<const struct ROL::TypeB::AlgorithmState<double> > (ROL::TypeB::Algorithm<double>::*)() const) &ROL::TypeB::Algorithm<double>::getState, "C++: ROL::TypeB::Algorithm<double>::getState() const --> class std::shared_ptr<const struct ROL::TypeB::AlgorithmState<double> >");
		cl.def("reset", (void (ROL::TypeB::Algorithm<double>::*)()) &ROL::TypeB::Algorithm<double>::reset, "C++: ROL::TypeB::Algorithm<double>::reset() --> void");
	}
}
