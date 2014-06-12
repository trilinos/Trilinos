/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_linsys_LinearSystemInterface_hpp
#define stk_linsys_LinearSystemInterface_hpp

#include <stk_linsys/FeiBaseIncludes.hpp>
#include <stk_linsys/DofMapper.hpp>

#include <Teuchos_ParameterList.hpp>

namespace stk_classic {
namespace linsys {

class LinearSystemInterface {
 public:
  virtual ~LinearSystemInterface() {}

  virtual void set_parameters(Teuchos::ParameterList& paramlist) = 0;

  virtual void synchronize_mappings_and_structure() = 0;
  virtual void create_fei_LinearSystem() = 0;
  virtual void finalize_assembly() = 0;

  /** Return DOF-mapping object */
  virtual const DofMapper& get_DofMapper() const = 0;

  /** Return DOF-mapping object */
  virtual DofMapper& get_DofMapper() = 0;

  /** set all matrix/vector coefficients to zero (if allocated) */
  virtual void reset_to_zero() = 0;

  /** Return fei::MatrixGraph object */
  virtual const fei::SharedPtr<fei::MatrixGraph> get_fei_MatrixGraph() const = 0;

  /** Return fei::MatrixGraph object */
  virtual fei::SharedPtr<fei::MatrixGraph> get_fei_MatrixGraph() = 0;

  /** Return fei::LinearSystem object */
  virtual const fei::SharedPtr<fei::LinearSystem> get_fei_LinearSystem() const = 0;

  /** Return fei::LinearSystem object */
  virtual fei::SharedPtr<fei::LinearSystem> get_fei_LinearSystem() = 0;

  /** Write matrix and vector objects out to file.
   * Use 'base_name' in the file-names.
   */
  virtual void write_files(const std::string& base_name) const = 0;

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
  virtual int solve(int & status, const Teuchos::ParameterList & params) = 0;

};//class LinearSystemInterface

}//namespace linsys
}//namespace stk_classic

#endif

