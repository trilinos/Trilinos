#ifndef MPIDATATYPEGENERATOR_H
#define MPIDATATYPEGENERATOR_H

#include <stddef.h>
#include <map>
#include <vector>

#include "Parallel.hpp"
#include "MPIFinalizationCallback.hpp"

namespace stk {
namespace impl {

class MPIDatatypeGenerator
{
  public:
    MPI_Datatype get_datatype(size_t size);

  private:
    MPIDatatypeGenerator();

    ~MPIDatatypeGenerator();

    void generate_datatype(size_t size);

    MPI_Datatype generate_datatype_contiguous(size_t size);

    MPI_Datatype generate_datatype_struct(size_t size);

    void destructor();

    std::map<size_t, MPI_Datatype> m_datatypes;

    const static std::vector<MPI_Datatype> m_candidateTypes;
    const static std::vector<size_t> m_candidateSizes;
    const static size_t m_numTypesForContiguous;
    MPIFinalizationCallback m_destructor;

    friend MPIDatatypeGenerator& get_datatype_generator();
};


MPIDatatypeGenerator& get_datatype_generator();

}

template <typename T>
MPI_Datatype generate_mpi_datatype()
{
  return impl::get_datatype_generator().get_datatype(sizeof(T));
}

}

#endif
