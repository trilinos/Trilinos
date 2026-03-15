#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_VerbosityLevel.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"

#include "pybind11/embed.h"

class scoped_interpreter : public pybind11::scoped_interpreter {

public:
  scoped_interpreter() {
    pybind11::module_::import("PyTrilinos2");
    pybind11::module_::import("torch");
  }
};

using mv_type = Tpetra::MultiVector<>;
using scalar_type = mv_type::scalar_type;
using local_ordinal_type = mv_type::local_ordinal_type;
using global_ordinal_type = mv_type::global_ordinal_type;
using node_type = mv_type::node_type;
using map_type =
    Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type>;

int main(int argc, char *argv[]) {
  Tpetra::ScopeGuard tpetraScope(&argc, &argv);
  {
    auto comm = Tpetra::getDefaultComm();
    auto out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

    const Tpetra::global_size_t numGblIndices = 10;
    const global_ordinal_type indexBase = 0;
    const int numVectors = 2;
    auto map = rcp(new map_type(numGblIndices, indexBase, comm));

    auto vec = rcp(new mv_type(map, numVectors));
    vec->putScalar(1.);
    vec->describe(*out, Teuchos::VERB_EXTREME);

    {
      // init python interpreter
      scoped_interpreter guard{};

      /////////////////////////////////////////////////////////////////
      // Vector
      {
        auto locals = pybind11::dict(pybind11::arg("vec") = *vec);
        pybind11::exec(R"(
        import PyTrilinos2
        import numpy as np

        vec.putScalar(2.)

        dataKokkos = vec.getLocalViewHost(PyTrilinos2.Tpetra.Access.ReadWriteStruct())
        print(dataKokkos, dataKokkos.extent(0))

        data_numpy = dataKokkos.numpy()
        print(type(data_numpy), np.info(data_numpy))
        data_numpy[1, 0] = 3.

        data_torch = dataKokkos.torch()
        print(type(data_torch))
        data_torch[2, 0] = 15.
                  )",
                       pybind11::globals(), locals);
      }

      vec->describe(*out, Teuchos::VERB_EXTREME);

    } // Python interpreter goes away
  }
  return 0;
}
