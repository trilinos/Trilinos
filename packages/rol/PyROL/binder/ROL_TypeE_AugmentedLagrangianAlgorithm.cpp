#include <ROL_Constraint.hpp>
#include <ROL_Elementwise_Function.hpp>
#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_Objective.hpp>
#include <ROL_PartitionedVector.hpp>
#include <ROL_TypeE_AugmentedLagrangianAlgorithm.hpp>
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

// ROL::TypeE::AugmentedLagrangianAlgorithm file:ROL_TypeE_AugmentedLagrangianAlgorithm.hpp line:59
struct PyCallBack_ROL_TypeE_AugmentedLagrangianAlgorithm_double_t : public ROL::TypeE::AugmentedLagrangianAlgorithm<double> {
	using ROL::TypeE::AugmentedLagrangianAlgorithm<double>::AugmentedLagrangianAlgorithm;

};

void bind_ROL_TypeE_AugmentedLagrangianAlgorithm(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // ROL::TypeE::AugmentedLagrangianAlgorithm file:ROL_TypeE_AugmentedLagrangianAlgorithm.hpp line:59
		pybind11::class_<ROL::TypeE::AugmentedLagrangianAlgorithm<double>, std::shared_ptr<ROL::TypeE::AugmentedLagrangianAlgorithm<double>>, PyCallBack_ROL_TypeE_AugmentedLagrangianAlgorithm_double_t, ROL::TypeE::Algorithm<double>> cl(M("ROL::TypeE"), "AugmentedLagrangianAlgorithm_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("list") );

		cl.def( pybind11::init( [](PyCallBack_ROL_TypeE_AugmentedLagrangianAlgorithm_double_t const &o){ return new PyCallBack_ROL_TypeE_AugmentedLagrangianAlgorithm_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::TypeE::AugmentedLagrangianAlgorithm<double> const &o){ return new ROL::TypeE::AugmentedLagrangianAlgorithm<double>(o); } ) );
		cl.def("setStatusTest", [](ROL::TypeE::Algorithm<double> &o, const class std::shared_ptr<class ROL::StatusTest<double> > & a0) -> void { return o.setStatusTest(a0); }, "", pybind11::arg("status"));
		cl.def("setStatusTest", (void (ROL::TypeE::Algorithm<double>::*)(const class std::shared_ptr<class ROL::StatusTest<double> > &, bool)) &ROL::TypeE::Algorithm<double>::setStatusTest, "C++: ROL::TypeE::Algorithm<double>::setStatusTest(const class std::shared_ptr<class ROL::StatusTest<double> > &, bool) --> void", pybind11::arg("status"), pybind11::arg("combineStatus"));
		cl.def("getState", (class std::shared_ptr<const struct ROL::TypeE::AlgorithmState<double> > (ROL::TypeE::Algorithm<double>::*)() const) &ROL::TypeE::Algorithm<double>::getState, "C++: ROL::TypeE::Algorithm<double>::getState() const --> class std::shared_ptr<const struct ROL::TypeE::AlgorithmState<double> >");
		cl.def("reset", (void (ROL::TypeE::Algorithm<double>::*)()) &ROL::TypeE::Algorithm<double>::reset, "C++: ROL::TypeE::Algorithm<double>::reset() --> void");
	}
}
