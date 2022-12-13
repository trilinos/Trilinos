#include "Teuchos_RCP.hpp"
#include "pybind11/pybind11.h"

template <typename A, typename trampoline_A>
Teuchos::RCP<A> customClone(const trampoline_A *ptr_to_this, std::string const &function_name) {
    pybind11::gil_scoped_acquire gil;
    pybind11::function overload = pybind11::get_overload(static_cast<const A *>(ptr_to_this), function_name.c_str());
    if (overload) {
        auto o = overload.operator()<pybind11::return_value_policy::reference>();
        auto self = pybind11::cast(ptr_to_this);
        auto cloned = self.attr(function_name.c_str())();

        auto keep_python_state_alive = Teuchos::rcp<pybind11::object>(new pybind11::object(cloned));
        auto ptr = pybind11::cast<trampoline_A *>(cloned);

        // aliasing shared_ptr: points to `A_trampoline* ptr` but refcounts the Python object
        return Teuchos::RCP<A>(keep_python_state_alive, ptr);		
    }
    pybind11::pybind11_fail("Tried to call pure virtual function \"customClone\"");
}