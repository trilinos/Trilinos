/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stk_linsys/LinearSystem.hpp>
#include <stk_linsys/ImplDetails.hpp>
#include <stk_mesh/base/GetBuckets.hpp>

#include <stk_linsys/LinsysFunctions.hpp>

#include <fei_Trilinos_Helpers.hpp>

namespace stk_classic {
namespace linsys {

LinearSystem::LinearSystem(MPI_Comm comm, fei::SharedPtr<fei::Factory> factory)
 : m_fei_factory(factory),
   m_dof_mapper(comm),
   m_fei_mgraph(new fei::MatrixGraph_Impl2(m_dof_mapper.get_fei_VectorSpace(), fei::SharedPtr<fei::VectorSpace>())),
   m_fei_linearsystem(),
   m_param_set()
{
}

LinearSystem::~LinearSystem()
{
}

void
LinearSystem::set_parameters(Teuchos::ParameterList& paramlist)
{
  Trilinos_Helpers::copy_parameterlist(paramlist, m_param_set);
  m_dof_mapper.get_fei_VectorSpace()->setParameters(m_param_set);
  if (m_fei_factory.get() != NULL) {
    m_fei_factory->parameters(m_param_set);
  }
  if (m_fei_mgraph.get() != NULL) {
    m_fei_mgraph->setParameters(m_param_set);
  }
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
  m_fei_linearsystem->parameters(m_param_set);

  fei::SharedPtr<fei::Matrix> mtx = m_fei_factory->createMatrix(m_fei_mgraph);
  mtx->parameters(m_param_set);
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

void
LinearSystem::reset_to_zero()
{
  fei::SharedPtr<fei::Matrix> mtx = m_fei_linearsystem->getMatrix();
  if (mtx.get() != NULL) {
    mtx->putScalar(0);
  }

  fei::SharedPtr<fei::Vector> rhs = m_fei_linearsystem->getRHS();
  if (rhs.get() != NULL) {
    rhs->putScalar(0);
  }
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

void
LinearSystem::write_files(const std::string& base_name) const
{
  fei::SharedPtr<fei::Matrix> A = m_fei_linearsystem->getMatrix();
  if (A.get() != NULL) {
    std::ostringstream ossA;
    ossA << "A_" << base_name << ".mtx";
    std::string Aname = ossA.str();
    A->writeToFile(Aname.c_str());
  }
  fei::SharedPtr<fei::Vector> b = m_fei_linearsystem->getRHS();
  if (b.get() != NULL) {
    std::ostringstream ossb;
    ossb << "b_" << base_name << ".vec";
    std::string bname = ossb.str();
    b->writeToFile(bname.c_str());
  }
}

int
LinearSystem::solve(int &status, const Teuchos::ParameterList & params )
{
  return fei_solve(status, *m_fei_linearsystem, params);
}

}//namespace linsys
}//namespace stk_classic

