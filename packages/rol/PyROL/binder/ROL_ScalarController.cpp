#include <ROL_BoundConstraint.hpp>
#include <ROL_Elementwise_Function.hpp>
#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_MoreauYosidaObjective.hpp>
#include <ROL_ScalarController.hpp>
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
#include <deque>
#include <memory>
#include <ostream>
#include <sstream> // __str__
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

void bind_ROL_ScalarController(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // ROL::ScalarController file:ROL_ScalarController.hpp line:54
		pybind11::class_<ROL::ScalarController<double,int>, Teuchos::RCP<ROL::ScalarController<double,int>>, ROL::VectorController<double,int>> cl(M("ROL"), "ScalarController_double_int_t", "");
		cl.def( pybind11::init( [](){ return new ROL::ScalarController<double,int>(); } ) );
		cl.def( pybind11::init( [](ROL::ScalarController<double,int> const &o){ return new ROL::ScalarController<double,int>(o); } ) );
		cl.def("get", (bool (ROL::ScalarController<double,int>::*)(double &, const int &)) &ROL::ScalarController<double, int>::get, "C++: ROL::ScalarController<double, int>::get(double &, const int &) --> bool", pybind11::arg("x"), pybind11::arg("param"));
		cl.def("set", (void (ROL::ScalarController<double,int>::*)(double, const int &)) &ROL::ScalarController<double, int>::set, "C++: ROL::ScalarController<double, int>::set(double, const int &) --> void", pybind11::arg("x"), pybind11::arg("param"));
		cl.def("reset", [](ROL::VectorController<double,int> &o) -> void { return o.reset(); }, "");
		cl.def("reset", (void (ROL::VectorController<double,int>::*)(bool)) &ROL::VectorController<double, int>::reset, "C++: ROL::VectorController<double, int>::reset(bool) --> void", pybind11::arg("flag"));
		cl.def("objectiveUpdate", [](ROL::VectorController<double,int> &o) -> void { return o.objectiveUpdate(); }, "");
		cl.def("objectiveUpdate", (void (ROL::VectorController<double,int>::*)(bool)) &ROL::VectorController<double, int>::objectiveUpdate, "C++: ROL::VectorController<double, int>::objectiveUpdate(bool) --> void", pybind11::arg("flag"));
		cl.def("constraintUpdate", [](ROL::VectorController<double,int> &o) -> void { return o.constraintUpdate(); }, "");
		cl.def("constraintUpdate", (void (ROL::VectorController<double,int>::*)(bool)) &ROL::VectorController<double, int>::constraintUpdate, "C++: ROL::VectorController<double, int>::constraintUpdate(bool) --> void", pybind11::arg("flag"));
		cl.def("objectiveUpdate", (void (ROL::VectorController<double,int>::*)(enum ROL::UpdateType)) &ROL::VectorController<double, int>::objectiveUpdate, "C++: ROL::VectorController<double, int>::objectiveUpdate(enum ROL::UpdateType) --> void", pybind11::arg("type"));
		cl.def("constraintUpdate", (void (ROL::VectorController<double,int>::*)(enum ROL::UpdateType)) &ROL::VectorController<double, int>::constraintUpdate, "C++: ROL::VectorController<double, int>::constraintUpdate(enum ROL::UpdateType) --> void", pybind11::arg("type"));
		cl.def("isNull", (bool (ROL::VectorController<double,int>::*)(const int &) const) &ROL::VectorController<double, int>::isNull, "C++: ROL::VectorController<double, int>::isNull(const int &) const --> bool", pybind11::arg("param"));
		cl.def("isComputed", (bool (ROL::VectorController<double,int>::*)(const int &) const) &ROL::VectorController<double, int>::isComputed, "C++: ROL::VectorController<double, int>::isComputed(const int &) const --> bool", pybind11::arg("param"));
		cl.def("allocate", (void (ROL::VectorController<double,int>::*)(const class ROL::Vector<double> &, const int &)) &ROL::VectorController<double, int>::allocate, "C++: ROL::VectorController<double, int>::allocate(const class ROL::Vector<double> &, const int &) --> void", pybind11::arg("x"), pybind11::arg("param"));
		cl.def("set", (const class Teuchos::RCP<class ROL::Vector<double> > (ROL::VectorController<double,int>::*)(const int &)) &ROL::VectorController<double, int>::set, "C++: ROL::VectorController<double, int>::set(const int &) --> const class Teuchos::RCP<class ROL::Vector<double> >", pybind11::arg("param"));
		cl.def("get", (const class Teuchos::RCP<const class ROL::Vector<double> > (ROL::VectorController<double,int>::*)(const int &) const) &ROL::VectorController<double, int>::get, "C++: ROL::VectorController<double, int>::get(const int &) const --> const class Teuchos::RCP<const class ROL::Vector<double> >", pybind11::arg("param"));
		cl.def("get", (bool (ROL::VectorController<double,int>::*)(class ROL::Vector<double> &, const int &)) &ROL::VectorController<double, int>::get, "C++: ROL::VectorController<double, int>::get(class ROL::Vector<double> &, const int &) --> bool", pybind11::arg("x"), pybind11::arg("param"));
		cl.def("set", (void (ROL::VectorController<double,int>::*)(const class ROL::Vector<double> &, const int &)) &ROL::VectorController<double, int>::set, "C++: ROL::VectorController<double, int>::set(const class ROL::Vector<double> &, const int &) --> void", pybind11::arg("x"), pybind11::arg("param"));
		cl.def("push", (void (ROL::VectorController<double,int>::*)(class ROL::VectorController<double, int> &) const) &ROL::VectorController<double, int>::push, "C++: ROL::VectorController<double, int>::push(class ROL::VectorController<double, int> &) const --> void", pybind11::arg("to"));
	}
}
