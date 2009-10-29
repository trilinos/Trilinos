
#include <stk_linsys/LinearSystem.hpp>
#include <stk_linsys/ImplDetails.hpp>
#include <stk_mesh/base/GetBuckets.hpp>

#include <stk_linsys/LinsysFunctions.hpp>

namespace stk {
namespace linsys {

LinearSystem::LinearSystem(MPI_Comm comm, fei::SharedPtr<fei::Factory> factory)
 : m_fei_factory(factory),
   m_dof_mapper(comm),
   m_fei_mgraph(new fei::MatrixGraph_Impl2(m_dof_mapper.get_fei_VectorSpace(), fei::SharedPtr<fei::VectorSpace>())),
   m_fei_linearsystem()
{
}

LinearSystem::~LinearSystem()
{
}

void
LinearSystem::synchronize_mappings_and_structure()
{
  m_fei_mgraph->initComplete();
  m_dof_mapper.finalize();
}

void
LinearSystem::create_fei_LinearSystem()
{
  m_fei_linearsystem = m_fei_factory->createLinearSystem(m_fei_mgraph);

  fei::SharedPtr<fei::Matrix> mtx = m_fei_factory->createMatrix(m_fei_mgraph);
  m_fei_linearsystem->setMatrix(mtx);
  fei::SharedPtr<fei::Vector> rhs = m_fei_factory->createVector(m_fei_mgraph);
  m_fei_linearsystem->setRHS(rhs);
  fei::SharedPtr<fei::Vector> soln = m_fei_factory->createVector(m_fei_mgraph,true);
  m_fei_linearsystem->setSolutionVector(soln);
}

void
LinearSystem::finalize_assembly()
{
  m_fei_linearsystem->loadComplete();
}

const DofMapper&
LinearSystem::get_DofMapper() const
{
  return m_dof_mapper;
}

DofMapper&
LinearSystem::get_DofMapper()
{
  return m_dof_mapper;
}

const fei::SharedPtr<fei::MatrixGraph>
LinearSystem::get_fei_MatrixGraph() const
{
  return m_fei_mgraph;
}

fei::SharedPtr<fei::MatrixGraph>
LinearSystem::get_fei_MatrixGraph()
{
  return m_fei_mgraph;
}

const fei::SharedPtr<fei::LinearSystem>
LinearSystem::get_fei_LinearSystem() const
{
  return m_fei_linearsystem;
}

fei::SharedPtr<fei::LinearSystem>
LinearSystem::get_fei_LinearSystem()
{
  return m_fei_linearsystem;
}

int
LinearSystem::solve(int &status, const Teuchos::ParameterList & params )
{
  return fei_solve(status, *m_fei_linearsystem, params);
}

}//namespace linsys
}//namespace stk

