/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_linsys_AggregateLinearSystem_hpp
#define stk_linsys_AggregateLinearSystem_hpp

#include <stk_linsys/FeiBaseIncludes.hpp>
#include <stk_linsys/DofMapper.hpp>
#include <stk_linsys/LinearSystemInterface.hpp>
#include <stk_linsys/LinearSystem.hpp>

#include <Teuchos_ParameterList.hpp>

namespace stk_classic {
namespace linsys {

/** Container for holding and manipulating collections of matrices and vectors.
 *
 * This class provides the ability to form a linear system in which the matrix is
 * a linear combination of other matrices, and the right-hand-side is a linear combination
 * of other vectors.
 * If n matrices A[0..n-1] and m vectors b[0..n-1] are each assembled, then the
 * 'aggregate' linear-system is formed with A = sum(alpha[i]*A[i]) and
 * b = sum(beta[i]*b[i]) where alpha and beta are arrays of scalars provided by
 * the calling code.
 */
class AggregateLinearSystem : public LinearSystemInterface {
 public:
  /** Constructor */
  AggregateLinearSystem(MPI_Comm comm, fei::SharedPtr<fei::Factory> factory, size_t num_matrices=1, size_t num_rhsvecs=1);

  /** Destructor */
  virtual ~AggregateLinearSystem();

  void set_parameters(Teuchos::ParameterList& paramlist);

  /** set the number of matrices and right-hand-sides */
  void set_num_matrices_rhsvecs(size_t num_matrices, size_t num_rhsvecs);

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

  /** Return the matrix at offset 'index' in the internally-stored array of matrices.
  */
  fei::SharedPtr<fei::Matrix> get_matrix(size_t index);

  /** Return the rhs-vec at offset 'index' in the internally-stored array of rhs-vectors.
  */
  fei::SharedPtr<fei::Vector> get_rhsvec(size_t index);

  /** Given arrays of scalars (which must have the same lengths as specified when
   * this class was constructed), form an aggregate linear system as described in
   * the class-description comments above.
   */
  void aggregate_system(const std::vector<double>& mat_scalars,
                        const std::vector<double>& rhs_scalars);

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

  void reset_to_zero();

  /** Return fei::MatrixGraph object */
  const fei::SharedPtr<fei::MatrixGraph> get_fei_MatrixGraph() const;

  /** Return fei::MatrixGraph object */
  fei::SharedPtr<fei::MatrixGraph> get_fei_MatrixGraph();

  /** Return fei::LinearSystem object */
  const fei::SharedPtr<fei::LinearSystem> get_fei_LinearSystem() const;

  /** Return fei::LinearSystem object */
  fei::SharedPtr<fei::LinearSystem> get_fei_LinearSystem();

  void write_files(const std::string& base_name) const;

  /** Solve the linear system
   * Note that the caller is expected to have already called the method
   * 'aggregate_system' if multiple matrices/rhs-vectors are being used.
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
  stk_classic::linsys::LinearSystem m_linear_system;

  std::vector<fei::SharedPtr<fei::Matrix> > m_matrices;
  std::vector<fei::SharedPtr<fei::Vector> > m_rhsvecs;
};//class AggregateLinearSystem

}//namespace linsys
}//namespace stk_classic

#endif

