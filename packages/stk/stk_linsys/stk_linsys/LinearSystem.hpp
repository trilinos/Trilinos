#ifndef stk_linsys_LinearSystem_hpp
#define stk_linsys_LinearSystem_hpp

#include <stk_linsys/FeiBaseIncludes.hpp>
#include <stk_linsys/DofMapper.hpp>

#include <Teuchos_ParameterList.hpp>

namespace stk {
namespace linsys {

/** Container of linear-system (matrix, vectors) and mapping objects.
 *
 */
class LinearSystem {
 public:
  /** Constructor */
  LinearSystem(MPI_Comm comm, fei::SharedPtr<fei::Factory> factory);

  /** Destructor */
  virtual ~LinearSystem();

  /** This is a collective call -- will hang if only a subset of processors
   * call it.
   * Internally calls fei::MatrixGraph::initComplete() and
   * DofMapper::finalize().
   */
  void synchronize_mappings_and_structure();

  /** Uses the fei::Factory (that was passed as a constructor argument) to
   *  create a fei::LinearSystem and populate it with a fei::Matrix and
   *  fei::Vectors.
   */
  void create_fei_LinearSystem();

  /** This is a collective call -- will hang if only a subset of processors
   *  call it.
   *  Internally calls fei::LinearSystem::loadComplete(), which in turn calls
   *  fei::Matrix::globalAssemble() and fei::Vector::gatherFromOverlap().
   * These operations perform communication to move shared contributions
   * to owning processors, etc.
   */
  void finalize_assembly();

  /** Return DOF-mapping object */
  const DofMapper& get_DofMapper() const;

  /** Return DOF-mapping object */
  DofMapper& get_DofMapper();

  /** Return fei::MatrixGraph object */
  const fei::SharedPtr<fei::MatrixGraph> get_fei_MatrixGraph() const;

  /** Return fei::MatrixGraph object */
  fei::SharedPtr<fei::MatrixGraph> get_fei_MatrixGraph();

  /** Return fei::LinearSystem object */
  const fei::SharedPtr<fei::LinearSystem> get_fei_LinearSystem() const;

  /** Return fei::LinearSystem object */
  fei::SharedPtr<fei::LinearSystem> get_fei_LinearSystem();

/** Solve the linear system
 *
 * @param status Output flag indicating the termination condition of the
 *  underlying linear-solver. Values are solver-specific. In general, 0
 *  indicates that the solver achieved a solution that satisfied the
 *  stopping test, within the iteration limit, etc. If an iterative solver
 *  fails to converge, this status value will generally be non-zero, but
 *  the actual value can vary by solver-library.
 *
 * @param params Teuchos::ParameterList for the solver
 *
 * @return error-code 0 if successful. Note that a 0 error-return does not
 *  necessarily mean that the underlying solver achieved a solution. It
 *  simply means that no fatal errors were encountered, such as allocation
 *  failures, etc.
 */
  int solve(int & status, const Teuchos::ParameterList & params);

 private:

  fei::SharedPtr<fei::Factory> m_fei_factory;
  DofMapper m_dof_mapper;
  fei::SharedPtr<fei::MatrixGraph> m_fei_mgraph;

  fei::SharedPtr<fei::LinearSystem> m_fei_linearsystem;
};//struct LinearSystem

}//namespace linsys
}//namespace stk

#endif

