#include <Teuchos_Ptr.hpp>
#include <iterator>
#include <memory>
#include <string>

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

void bind_Teuchos_Ptr(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	// Teuchos::PtrPrivateUtilityPack::throw_null(const std::string &) file:Teuchos_Ptr.hpp line:55
	M("Teuchos::PtrPrivateUtilityPack").def("throw_null", (void (*)(const std::string &)) &Teuchos::PtrPrivateUtilityPack::throw_null, "C++: Teuchos::PtrPrivateUtilityPack::throw_null(const std::string &) --> void", pybind11::arg("type_name"));

}
