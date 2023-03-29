#include "MPIDatatypeGenerator.hpp"
#include "stk_util/util/ReportHandler.hpp"


namespace stk {
namespace impl {

const std::vector<MPI_Datatype> MPIDatatypeGenerator::m_candidateTypes = {MPI_LONG_LONG, MPI_LONG, MPI_INT, MPI_SHORT, MPI_BYTE};
const std::vector<size_t> MPIDatatypeGenerator::m_candidateSizes       = {sizeof(long long), sizeof(long), sizeof(int), sizeof(short), 1};
const size_t MPIDatatypeGenerator::m_numTypesForContiguous = 3;

MPIDatatypeGenerator::MPIDatatypeGenerator() : m_destructor([this]() { destructor(); }) {}

MPIDatatypeGenerator::~MPIDatatypeGenerator() 
{ 
  m_destructor.destructor();
}


MPI_Datatype MPIDatatypeGenerator::get_datatype(size_t size)
{
  STK_ThrowRequireMsg(size > 0, "The sizeof(type) cannot return 0.  Where did you get this type from?");
  if (m_datatypes.find(size) == m_datatypes.end())
  {
    generate_datatype(size);
  }

  return m_datatypes[size];
}  


void MPIDatatypeGenerator::generate_datatype(size_t size)
{
  STK_ThrowAssertMsg(m_datatypes.find(size) == m_datatypes.end(), "Cannot create new datatype for size that already exists");

  MPI_Datatype newDatatype = generate_datatype_contiguous(size);
  if (newDatatype == MPI_DATATYPE_NULL)
  {
    newDatatype = generate_datatype_struct(size);
  }

  MPI_Datatype newDatatypeWithoutPadding;

  MPI_Type_create_resized(newDatatype, 0, size, &newDatatypeWithoutPadding);
  MPI_Type_commit(&newDatatypeWithoutPadding);

  MPI_Type_free(&newDatatype);

  m_datatypes[size] = newDatatypeWithoutPadding;
}


MPI_Datatype MPIDatatypeGenerator::generate_datatype_contiguous(size_t size)
{
  for (size_t i=0; i <= m_numTypesForContiguous; ++i)
  {
    if (size % m_candidateSizes[i] == 0)
    {
      MPI_Datatype newDatatype;
      MPI_Type_contiguous(size / m_candidateSizes[i], m_candidateTypes[i], &newDatatype);

      return newDatatype;
    }
  }

  return MPI_DATATYPE_NULL;
}

MPI_Datatype MPIDatatypeGenerator::generate_datatype_struct(size_t size)
{
  std::vector<int> blocklengths;
  std::vector<MPI_Aint> displacements;
  std::vector<MPI_Datatype> types;

  size_t size_remaining = size;
  for (size_t i=0; i < m_candidateTypes.size(); ++i)
  {
    size_t num_entries = size_remaining/m_candidateSizes[i];
    if (num_entries > 0)
    {
      types.push_back(m_candidateTypes[i]);
      blocklengths.push_back(num_entries);
      displacements.push_back(size - size_remaining);

      size_t block_size = num_entries * m_candidateSizes[i];
      STK_ThrowRequireMsg(block_size <= size_remaining, "remaining size cannot be less than block size");
      size_remaining -= block_size;
    }
  }

  MPI_Datatype newDatatype;
  MPI_Type_create_struct(types.size(), blocklengths.data(), displacements.data(), types.data(), &newDatatype);

  return newDatatype;
}

void MPIDatatypeGenerator::destructor()
{
  for (auto& sizeDatatypePair : m_datatypes)
  {
    MPI_Type_free(&sizeDatatypePair.second);
  }
}


MPIDatatypeGenerator& get_datatype_generator()
{
  static MPIDatatypeGenerator generator;
  return generator;
}

}

}