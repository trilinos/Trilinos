#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_VerbosityLevel.hpp"
#include "Tpetra_Access.hpp"
#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_MultiVector.hpp>

#include <pybind11/embed.h>

#include <torch/torch.h>
#include <torch/extension.h>

namespace py = pybind11;
using namespace py::literals;

using map_type = Tpetra::Map<>;
using vec_type = Tpetra::Vector<>;
using mv_type = Tpetra::MultiVector<>;
using scalar_type = vec_type::scalar_type;
using local_ordinal_type = vec_type::local_ordinal_type;
using global_ordinal_type = vec_type::global_ordinal_type;


template <class T>
torch::Tensor convert_kokkos_to_torch(T kokkos_array) {
  torch::Tensor tensor;
  if constexpr (T::rank == 1)
    tensor = torch::from_blob(kokkos_array.data(), {kokkos_array.extent(0)}, torch::CppTypeToScalarType<typename T::value_type>::value);
  else if constexpr (T::rank == 2)
    tensor = torch::from_blob(kokkos_array.data(), {kokkos_array.extent(0), kokkos_array.extent(1)}, torch::CppTypeToScalarType<typename T::value_type>::value);
  return tensor;
}

int main(int argc, char *argv[]) {
  Tpetra::ScopeGuard tpetraScope(&argc, &argv);
  {
    auto comm = Tpetra::getDefaultComm();
    const Tpetra::global_size_t numGblIndices = 50;

    auto out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

    const global_ordinal_type indexBase = 0;
    auto map = rcp(new map_type(numGblIndices, indexBase, comm));

    auto vec = rcp(new vec_type(map));
    vec->putScalar(1.);
    vec->describe(*out, Teuchos::VERB_EXTREME);

    {
      // init python interpreter
      py::scoped_interpreter guard{};
      py::module_ PyTrilinos2 = py::module_::import("PyTrilinos2");
      py::module_ Torch = py::module_::import("torch");

      /////////////////////////////////////////////////////////////////
      // Vector
      {
        auto locals = py::dict("vec"_a = *vec);
        py::exec(R"(
        vec.putScalar(2.)
        data = vec.numpy_view_host()
        data[:] = 3.
                  )", py::globals(), locals);
      }

      vec->describe(*out, Teuchos::VERB_EXTREME);

      /////////////////////////////////////////////////////////////////
      // MultiVector

      auto mv = rcp(new mv_type(map, 2));
      mv->putScalar(1.);
      mv->describe(*out, Teuchos::VERB_EXTREME);

      {
        auto locals = py::dict("mv"_a = *mv);
        py::exec(R"(
        mv.putScalar(2.)
        data = mv.numpy_view_host()
        data[:, :] = 3.
                  )", py::globals(), locals);
      }
      mv->describe(*out, Teuchos::VERB_EXTREME);

      {
        auto data = vec->getLocalViewHost(Tpetra::Access::ReadWrite);
        std::cout << convert_kokkos_to_torch(data) << std::endl;

        auto data1d = Kokkos::subview(data, Kokkos::ALL(), 0);
        auto torch_tensor1d = convert_kokkos_to_torch(data1d);
        std::cout << torch_tensor1d << std::endl;

        {
          auto locals = py::dict("tens"_a = torch_tensor1d);
          py::exec(R"(
        print(tens)
        tens[0]=4.
        print(tens)
                    )", py::globals(), locals);
        }
      }

      vec->describe(*out, Teuchos::VERB_EXTREME);

    } // Python interpreter goes away

  }
  return 0;
}
