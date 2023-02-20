#include <ROL_BoundConstraint.hpp>
#include <ROL_BoundConstraint_Partitioned.hpp>
#include <ROL_Bounds.hpp>
#include <ROL_Constraint.hpp>
#include <ROL_Constraint_Partitioned.hpp>
#include <ROL_Elementwise_Function.hpp>
#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_Objective.hpp>
#include <ROL_PartitionedVector.hpp>
#include <ROL_PolyhedralProjection.hpp>
#include <ROL_Problem.hpp>
#include <ROL_Secant.hpp>
#include <ROL_SlacklessObjective.hpp>
#include <ROL_TypeB_InteriorPointAlgorithm.hpp>
#include <ROL_TypeB_KelleySachsAlgorithm.hpp>
#include <ROL_TypeB_LSecantBAlgorithm.hpp>
#include <ROL_TypeB_PrimalDualActiveSetAlgorithm.hpp>
#include <ROL_TypeB_SpectralGradientAlgorithm.hpp>
#include <ROL_Types.hpp>
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

// ROL::TypeB::InteriorPointAlgorithm file:ROL_TypeB_InteriorPointAlgorithm.hpp line:58
struct PyCallBack_ROL_TypeB_InteriorPointAlgorithm_double_t : public ROL::TypeB::InteriorPointAlgorithm<double> {
	using ROL::TypeB::InteriorPointAlgorithm<double>::InteriorPointAlgorithm;

};

// ROL::TypeB::PrimalDualActiveSetAlgorithm file:ROL_TypeB_PrimalDualActiveSetAlgorithm.hpp line:59
struct PyCallBack_ROL_TypeB_PrimalDualActiveSetAlgorithm_double_t : public ROL::TypeB::PrimalDualActiveSetAlgorithm<double> {
	using ROL::TypeB::PrimalDualActiveSetAlgorithm<double>::PrimalDualActiveSetAlgorithm;

};

// ROL::TypeB::KelleySachsAlgorithm file:ROL_TypeB_KelleySachsAlgorithm.hpp line:61
struct PyCallBack_ROL_TypeB_KelleySachsAlgorithm_double_t : public ROL::TypeB::KelleySachsAlgorithm<double> {
	using ROL::TypeB::KelleySachsAlgorithm<double>::KelleySachsAlgorithm;

};

// ROL::TypeB::SpectralGradientAlgorithm file:ROL_TypeB_SpectralGradientAlgorithm.hpp line:57
struct PyCallBack_ROL_TypeB_SpectralGradientAlgorithm_double_t : public ROL::TypeB::SpectralGradientAlgorithm<double> {
	using ROL::TypeB::SpectralGradientAlgorithm<double>::SpectralGradientAlgorithm;

};

// ROL::TypeB::LSecantBAlgorithm file:ROL_TypeB_LSecantBAlgorithm.hpp line:61
struct PyCallBack_ROL_TypeB_LSecantBAlgorithm_double_t : public ROL::TypeB::LSecantBAlgorithm<double> {
	using ROL::TypeB::LSecantBAlgorithm<double>::LSecantBAlgorithm;

};

void bind_ROL_TypeB_InteriorPointAlgorithm(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // ROL::TypeB::InteriorPointAlgorithm file:ROL_TypeB_InteriorPointAlgorithm.hpp line:58
		pybind11::class_<ROL::TypeB::InteriorPointAlgorithm<double>, std::shared_ptr<ROL::TypeB::InteriorPointAlgorithm<double>>, PyCallBack_ROL_TypeB_InteriorPointAlgorithm_double_t, ROL::TypeB::Algorithm<double>> cl(M("ROL::TypeB"), "InteriorPointAlgorithm_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("list") );

		cl.def( pybind11::init( [](PyCallBack_ROL_TypeB_InteriorPointAlgorithm_double_t const &o){ return new PyCallBack_ROL_TypeB_InteriorPointAlgorithm_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::TypeB::InteriorPointAlgorithm<double> const &o){ return new ROL::TypeB::InteriorPointAlgorithm<double>(o); } ) );
		cl.def("setStatusTest", [](ROL::TypeB::Algorithm<double> &o, const class std::shared_ptr<class ROL::StatusTest<double> > & a0) -> void { return o.setStatusTest(a0); }, "", pybind11::arg("status"));
		cl.def("setStatusTest", (void (ROL::TypeB::Algorithm<double>::*)(const class std::shared_ptr<class ROL::StatusTest<double> > &, const bool)) &ROL::TypeB::Algorithm<double>::setStatusTest, "C++: ROL::TypeB::Algorithm<double>::setStatusTest(const class std::shared_ptr<class ROL::StatusTest<double> > &, const bool) --> void", pybind11::arg("status"), pybind11::arg("combineStatus"));
		cl.def("getState", (class std::shared_ptr<const struct ROL::TypeB::AlgorithmState<double> > (ROL::TypeB::Algorithm<double>::*)() const) &ROL::TypeB::Algorithm<double>::getState, "C++: ROL::TypeB::Algorithm<double>::getState() const --> class std::shared_ptr<const struct ROL::TypeB::AlgorithmState<double> >");
		cl.def("reset", (void (ROL::TypeB::Algorithm<double>::*)()) &ROL::TypeB::Algorithm<double>::reset, "C++: ROL::TypeB::Algorithm<double>::reset() --> void");
	}
	{ // ROL::TypeB::PrimalDualActiveSetAlgorithm file:ROL_TypeB_PrimalDualActiveSetAlgorithm.hpp line:59
		pybind11::class_<ROL::TypeB::PrimalDualActiveSetAlgorithm<double>, std::shared_ptr<ROL::TypeB::PrimalDualActiveSetAlgorithm<double>>, PyCallBack_ROL_TypeB_PrimalDualActiveSetAlgorithm_double_t, ROL::TypeB::Algorithm<double>> cl(M("ROL::TypeB"), "PrimalDualActiveSetAlgorithm_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](class Teuchos::ParameterList & a0){ return new ROL::TypeB::PrimalDualActiveSetAlgorithm<double>(a0); }, [](class Teuchos::ParameterList & a0){ return new PyCallBack_ROL_TypeB_PrimalDualActiveSetAlgorithm_double_t(a0); } ), "doc");
		cl.def( pybind11::init<class Teuchos::ParameterList &, const class std::shared_ptr<class ROL::Secant<double> > &>(), pybind11::arg("list"), pybind11::arg("secant") );

		cl.def( pybind11::init( [](PyCallBack_ROL_TypeB_PrimalDualActiveSetAlgorithm_double_t const &o){ return new PyCallBack_ROL_TypeB_PrimalDualActiveSetAlgorithm_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::TypeB::PrimalDualActiveSetAlgorithm<double> const &o){ return new ROL::TypeB::PrimalDualActiveSetAlgorithm<double>(o); } ) );
		cl.def("setStatusTest", [](ROL::TypeB::Algorithm<double> &o, const class std::shared_ptr<class ROL::StatusTest<double> > & a0) -> void { return o.setStatusTest(a0); }, "", pybind11::arg("status"));
		cl.def("setStatusTest", (void (ROL::TypeB::Algorithm<double>::*)(const class std::shared_ptr<class ROL::StatusTest<double> > &, const bool)) &ROL::TypeB::Algorithm<double>::setStatusTest, "C++: ROL::TypeB::Algorithm<double>::setStatusTest(const class std::shared_ptr<class ROL::StatusTest<double> > &, const bool) --> void", pybind11::arg("status"), pybind11::arg("combineStatus"));
		cl.def("getState", (class std::shared_ptr<const struct ROL::TypeB::AlgorithmState<double> > (ROL::TypeB::Algorithm<double>::*)() const) &ROL::TypeB::Algorithm<double>::getState, "C++: ROL::TypeB::Algorithm<double>::getState() const --> class std::shared_ptr<const struct ROL::TypeB::AlgorithmState<double> >");
		cl.def("reset", (void (ROL::TypeB::Algorithm<double>::*)()) &ROL::TypeB::Algorithm<double>::reset, "C++: ROL::TypeB::Algorithm<double>::reset() --> void");
	}
	{ // ROL::TypeB::KelleySachsAlgorithm file:ROL_TypeB_KelleySachsAlgorithm.hpp line:61
		pybind11::class_<ROL::TypeB::KelleySachsAlgorithm<double>, std::shared_ptr<ROL::TypeB::KelleySachsAlgorithm<double>>, PyCallBack_ROL_TypeB_KelleySachsAlgorithm_double_t, ROL::TypeB::Algorithm<double>> cl(M("ROL::TypeB"), "KelleySachsAlgorithm_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](class Teuchos::ParameterList & a0){ return new ROL::TypeB::KelleySachsAlgorithm<double>(a0); }, [](class Teuchos::ParameterList & a0){ return new PyCallBack_ROL_TypeB_KelleySachsAlgorithm_double_t(a0); } ), "doc");
		cl.def( pybind11::init<class Teuchos::ParameterList &, const class std::shared_ptr<class ROL::Secant<double> > &>(), pybind11::arg("list"), pybind11::arg("secant") );

		cl.def( pybind11::init( [](PyCallBack_ROL_TypeB_KelleySachsAlgorithm_double_t const &o){ return new PyCallBack_ROL_TypeB_KelleySachsAlgorithm_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::TypeB::KelleySachsAlgorithm<double> const &o){ return new ROL::TypeB::KelleySachsAlgorithm<double>(o); } ) );
		cl.def("setStatusTest", [](ROL::TypeB::Algorithm<double> &o, const class std::shared_ptr<class ROL::StatusTest<double> > & a0) -> void { return o.setStatusTest(a0); }, "", pybind11::arg("status"));
		cl.def("setStatusTest", (void (ROL::TypeB::Algorithm<double>::*)(const class std::shared_ptr<class ROL::StatusTest<double> > &, const bool)) &ROL::TypeB::Algorithm<double>::setStatusTest, "C++: ROL::TypeB::Algorithm<double>::setStatusTest(const class std::shared_ptr<class ROL::StatusTest<double> > &, const bool) --> void", pybind11::arg("status"), pybind11::arg("combineStatus"));
		cl.def("getState", (class std::shared_ptr<const struct ROL::TypeB::AlgorithmState<double> > (ROL::TypeB::Algorithm<double>::*)() const) &ROL::TypeB::Algorithm<double>::getState, "C++: ROL::TypeB::Algorithm<double>::getState() const --> class std::shared_ptr<const struct ROL::TypeB::AlgorithmState<double> >");
		cl.def("reset", (void (ROL::TypeB::Algorithm<double>::*)()) &ROL::TypeB::Algorithm<double>::reset, "C++: ROL::TypeB::Algorithm<double>::reset() --> void");
	}
	{ // ROL::TypeB::SpectralGradientAlgorithm file:ROL_TypeB_SpectralGradientAlgorithm.hpp line:57
		pybind11::class_<ROL::TypeB::SpectralGradientAlgorithm<double>, std::shared_ptr<ROL::TypeB::SpectralGradientAlgorithm<double>>, PyCallBack_ROL_TypeB_SpectralGradientAlgorithm_double_t, ROL::TypeB::Algorithm<double>> cl(M("ROL::TypeB"), "SpectralGradientAlgorithm_double_t", "", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::ParameterList &>(), pybind11::arg("list") );

		cl.def( pybind11::init( [](PyCallBack_ROL_TypeB_SpectralGradientAlgorithm_double_t const &o){ return new PyCallBack_ROL_TypeB_SpectralGradientAlgorithm_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::TypeB::SpectralGradientAlgorithm<double> const &o){ return new ROL::TypeB::SpectralGradientAlgorithm<double>(o); } ) );
		cl.def("setStatusTest", [](ROL::TypeB::Algorithm<double> &o, const class std::shared_ptr<class ROL::StatusTest<double> > & a0) -> void { return o.setStatusTest(a0); }, "", pybind11::arg("status"));
		cl.def("setStatusTest", (void (ROL::TypeB::Algorithm<double>::*)(const class std::shared_ptr<class ROL::StatusTest<double> > &, const bool)) &ROL::TypeB::Algorithm<double>::setStatusTest, "C++: ROL::TypeB::Algorithm<double>::setStatusTest(const class std::shared_ptr<class ROL::StatusTest<double> > &, const bool) --> void", pybind11::arg("status"), pybind11::arg("combineStatus"));
		cl.def("getState", (class std::shared_ptr<const struct ROL::TypeB::AlgorithmState<double> > (ROL::TypeB::Algorithm<double>::*)() const) &ROL::TypeB::Algorithm<double>::getState, "C++: ROL::TypeB::Algorithm<double>::getState() const --> class std::shared_ptr<const struct ROL::TypeB::AlgorithmState<double> >");
		cl.def("reset", (void (ROL::TypeB::Algorithm<double>::*)()) &ROL::TypeB::Algorithm<double>::reset, "C++: ROL::TypeB::Algorithm<double>::reset() --> void");
	}
	{ // ROL::TypeB::LSecantBAlgorithm file:ROL_TypeB_LSecantBAlgorithm.hpp line:61
		pybind11::class_<ROL::TypeB::LSecantBAlgorithm<double>, std::shared_ptr<ROL::TypeB::LSecantBAlgorithm<double>>, PyCallBack_ROL_TypeB_LSecantBAlgorithm_double_t, ROL::TypeB::Algorithm<double>> cl(M("ROL::TypeB"), "LSecantBAlgorithm_double_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](class Teuchos::ParameterList & a0){ return new ROL::TypeB::LSecantBAlgorithm<double>(a0); }, [](class Teuchos::ParameterList & a0){ return new PyCallBack_ROL_TypeB_LSecantBAlgorithm_double_t(a0); } ), "doc");
		cl.def( pybind11::init<class Teuchos::ParameterList &, const class std::shared_ptr<class ROL::Secant<double> > &>(), pybind11::arg("list"), pybind11::arg("secant") );

		cl.def( pybind11::init( [](PyCallBack_ROL_TypeB_LSecantBAlgorithm_double_t const &o){ return new PyCallBack_ROL_TypeB_LSecantBAlgorithm_double_t(o); } ) );
		cl.def( pybind11::init( [](ROL::TypeB::LSecantBAlgorithm<double> const &o){ return new ROL::TypeB::LSecantBAlgorithm<double>(o); } ) );
		cl.def("setStatusTest", [](ROL::TypeB::Algorithm<double> &o, const class std::shared_ptr<class ROL::StatusTest<double> > & a0) -> void { return o.setStatusTest(a0); }, "", pybind11::arg("status"));
		cl.def("setStatusTest", (void (ROL::TypeB::Algorithm<double>::*)(const class std::shared_ptr<class ROL::StatusTest<double> > &, const bool)) &ROL::TypeB::Algorithm<double>::setStatusTest, "C++: ROL::TypeB::Algorithm<double>::setStatusTest(const class std::shared_ptr<class ROL::StatusTest<double> > &, const bool) --> void", pybind11::arg("status"), pybind11::arg("combineStatus"));
		cl.def("getState", (class std::shared_ptr<const struct ROL::TypeB::AlgorithmState<double> > (ROL::TypeB::Algorithm<double>::*)() const) &ROL::TypeB::Algorithm<double>::getState, "C++: ROL::TypeB::Algorithm<double>::getState() const --> class std::shared_ptr<const struct ROL::TypeB::AlgorithmState<double> >");
		cl.def("reset", (void (ROL::TypeB::Algorithm<double>::*)()) &ROL::TypeB::Algorithm<double>::reset, "C++: ROL::TypeB::Algorithm<double>::reset() --> void");
	}
}
