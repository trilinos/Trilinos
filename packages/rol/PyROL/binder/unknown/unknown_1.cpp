#include <PyROL_ETI_helper.hpp>
#include <ROL_Elementwise_Function.hpp>
#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_LinearOperator.hpp>
#include <ROL_QuadraticObjective.hpp>
#include <ROL_Vector.hpp>
#include <Teuchos_ENull.hpp>
#include <Teuchos_RCPDecl.hpp>
#include <Teuchos_RCPNode.hpp>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include <functional>
#include <pybind11/smart_holder.h>
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

void bind_unknown_unknown_1(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	// ROL::PyROL::foo(class ROL::QuadraticObjective<double>) file: line:26
	M("ROL::PyROL").def("foo", (void (*)(class ROL::QuadraticObjective<double>)) &ROL::PyROL::foo<ROL::QuadraticObjective<double>>, "C++: ROL::PyROL::foo(class ROL::QuadraticObjective<double>) --> void", pybind11::arg("a"));

}
