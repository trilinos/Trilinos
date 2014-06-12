/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stk_linsys/AggregateLinearSystem.hpp>
#include <stk_linsys/LinearSystem.hpp>
#include <stk_linsys/ImplDetails.hpp>
#include <stk_mesh/base/GetBuckets.hpp>

#include <stk_linsys/LinsysFunctions.hpp>

namespace stk_classic {
namespace linsys {

AggregateLinearSystem::AggregateLinearSystem(MPI_Comm comm, fei::SharedPtr<fei::Factory> factory, size_t num_matrices, size_t num_rhsvecs)
 : m_fei_factory(factory),
   m_linear_system(comm, factory),
   m_matrices(num_matrices),
   m_rhsvecs(num_rhsvecs)
{
}

AggregateLinearSystem::~AggregateLinearSystem()
{
}

void
AggregateLinearSystem::set_parameters(Teuchos::ParameterList& paramlist)
{
  m_linear_system.set_parameters(paramlist);
}

void
AggregateLinearSystem::set_num_matrices_rhsvecs(size_t num_matrices, size_t num_rhsvecs)
{
  m_matrices.resize(num_matrices);
  m_rhsvecs.resize(num_rhsvecs);
}

void
AggregateLinearSystem::synchronize_mappings_and_structure()
{
  m_linear_system.synchronize_mappings_and_structure();
}

void
AggregateLinearSystem::create_fei_LinearSystem()
{
  m_linear_system.create_fei_LinearSystem();

  fei::SharedPtr<fei::MatrixGraph> mgraph = m_linear_system.get_fei_MatrixGraph();

  for(size_t i=0; i<m_matrices.size(); ++i) {
    m_matrices[i] = m_fei_factory->createMatrix(mgraph);
  }

  bool is_soln_vec = false;
  for(size_t i=0; i<m_rhsvecs.size(); ++i) {
    m_rhsvecs[i] = m_fei_factory->createVector(mgraph,is_soln_vec);
  }
}

fei::SharedPtr<fei::Matrix>
AggregateLinearSystem::get_matrix(size_t index)
{
  if (index >= m_matrices.size()) {
    throw std::runtime_error("stk_classic::linsys::AggregateLinearSystem::get_matrix ERROR, index out of range.");
  }

  return m_matrices[index];
}

fei::SharedPtr<fei::Vector>
AggregateLinearSystem::get_rhsvec(size_t index)
{
  if (index >= m_rhsvecs.size()) {
    throw std::runtime_error("stk_classic::linsys::AggregateLinearSystem::get_rhsvec ERROR, index out of range.");
  }

  return m_rhsvecs[index];
}

void
AggregateLinearSystem::aggregate_system(const std::vector<double>& mat_scalars,
                                        const std::vector<double>& rhs_scalars)
{
  if (mat_scalars.size() != m_matrices.size()) {
    throw std::runtime_error("stk_classic::linsys::AggregateLinearSystem::aggregate_system ERROR, mat_scalars.size() != m_matrices.size().");
  }

  if (rhs_scalars.size() != m_rhsvecs.size()) {
    throw std::runtime_error("stk_classic::linsys::AggregateLinearSystem::aggregate_system ERROR, rhs_scalars.size() != m_rhsvecs.size().");
  }

  fei::SharedPtr<fei::LinearSystem> fei_linsys = m_linear_system.get_fei_LinearSystem();
  fei::SharedPtr<fei::Matrix> matrix = fei_linsys->getMatrix();
  fei::SharedPtr<fei::Vector> rhsvec = fei_linsys->getRHS();

  matrix->gatherFromOverlap();
  matrix->putScalar(0.0);

  for(size_t i=0; i<m_matrices.size(); ++i) {
    stk_classic::linsys::add_matrix_to_matrix(mat_scalars[i], *m_matrices[i], *matrix);
  }

  rhsvec->gatherFromOverlap();
  rhsvec->putScalar(0.0);

  for(size_t i=0; i<m_rhsvecs.size(); ++i) {
    stk_classic::linsys::add_vector_to_vector(rhs_scalars[i], *m_rhsvecs[i], *rhsvec);
  }
}

void
AggregateLinearSystem::finalize_assembly()
{
  m_linear_system.finalize_assembly();
}

const DofMapper&
AggregateLinearSystem::get_DofMapper() const
{
  return m_linear_system.get_DofMapper();
}

DofMapper&
AggregateLinearSystem::get_DofMapper()
{
  return m_linear_system.get_DofMapper();
}

void
AggregateLinearSystem::reset_to_zero()
{
  for(size_t i=0; i<m_matrices.size(); ++i) {
    if (m_matrices[i].get() != NULL) m_matrices[i]->putScalar(0);
  }
  for(size_t i=0; i<m_rhsvecs.size(); ++i) {
    if (m_rhsvecs[i].get() != NULL) m_rhsvecs[i]->putScalar(0);
  }
}

const fei::SharedPtr<fei::MatrixGraph>
AggregateLinearSystem::get_fei_MatrixGraph() const
{
  return m_linear_system.get_fei_MatrixGraph();
}

fei::SharedPtr<fei::MatrixGraph>
AggregateLinearSystem::get_fei_MatrixGraph()
{
  return m_linear_system.get_fei_MatrixGraph();
}

const fei::SharedPtr<fei::LinearSystem>
AggregateLinearSystem::get_fei_LinearSystem() const
{
  return m_linear_system.get_fei_LinearSystem();
}

fei::SharedPtr<fei::LinearSystem>
AggregateLinearSystem::get_fei_LinearSystem()
{
  return m_linear_system.get_fei_LinearSystem();
}

void
AggregateLinearSystem::write_files(const std::string& base_name) const
{
  for(size_t i=0; i<m_matrices.size(); ++i) {
    if (m_matrices[i].get() == NULL) continue;
    std::ostringstream ossA;
    ossA << "A_" << base_name << ".mat"<<i<<".mtx";
    std::string Aname = ossA.str();
    m_matrices[i]->writeToFile(Aname.c_str());
  }
  for(size_t i=0; i<m_rhsvecs.size(); ++i) {
    if (m_rhsvecs[i].get() == NULL) continue;
    std::ostringstream ossb;
    ossb << "b_" << base_name << ".vec"<<i<<".vec";
    std::string bname = ossb.str();
    m_rhsvecs[i]->writeToFile(bname.c_str());
  }
}

int
AggregateLinearSystem::solve(int &status, const Teuchos::ParameterList & params )
{
  return m_linear_system.solve(status, params);
}

}//namespace linsys
}//namespace stk_classic

