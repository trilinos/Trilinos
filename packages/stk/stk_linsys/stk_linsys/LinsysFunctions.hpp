/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_linsys_LinsysFunctions_hpp
#define stk_linsys_LinsysFunctions_hpp

#include <stk_linsys/LinearSystemInterface.hpp>

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/BulkData.hpp>

#include <Teuchos_ParameterList.hpp>


namespace stk {
namespace linsys {

/** Add connectivities (matrix-graph sparsity contributions) to the
    fei::MatrixGraph object in the specified LinearSystemInterface object.

  Connectivities are connections between two types of entities that
  are related to each other in the mesh. The most common example of
  this is element-to-node connectivities. So when from_type is element
  and to_connected_type is node, the matrix-graph will be populated with
  connectivities for nodes connected to each element selected by the
  given selector.

*/
void add_connectivities(stk::linsys::LinearSystemInterface& ls,
                        stk::mesh::EntityRank from_type,
                        stk::mesh::EntityRank to_connected_type,
                        const stk::mesh::FieldBase& field,
                        const stk::mesh::Selector& selector,
                        const stk::mesh::BulkData& mesh_bulk);

/** Apply a Dirichlet boundary-condition for the specified field, on
 * entities of the specified entity-type which are members of the
 * specified Part.
 * Obtains the list of entity-ids and passes that along with
 * the prescribed value, etc., to the fei::LinearSystem object.
 * The corresponding modifications to the matrix and vector will be made
 * when LinearSystemInterface::finalize_assembly() is called.
 */
void dirichlet_bc(stk::linsys::LinearSystemInterface& ls,
                  const stk::mesh::BulkData& mesh,
                  const stk::mesh::Part& bcpart,
                  stk::mesh::EntityRank entity_rank,
                  const stk::mesh::FieldBase& field,
                  unsigned field_component,
                  double prescribed_value);


/** Create an fei::Solver instance and perform a linear-solve on the
 * matrix and vectors contained in the given fei::LinearSystem instance.
 *
 * @param status Output flag indicating the termination condition of the
 *  underlying linear-solver. Values are solver-specific. In general, 0
 *  indicates that the solver achieved a solution that satisfied the
 *  stopping test, within the iteration limit, etc. If an iterative solver
 *  fails to converge, this status value will generally be non-zero, but
 *  the actual value can vary by solver-library.
 *
 * @return error-code 0 if successful. Note that a 0 error-return does not
 *  necessarily mean that the underlying solver achieved a solution. It
 *  simply means that no fatal errors were encountered, such as allocation
 *  failures, etc.
 */
int fei_solve(int & status, fei::LinearSystem &fei_ls, const Teuchos::ParameterList & params);

/** Return a 2-norm of |b - A*x|
 */
double compute_residual_norm2(fei::LinearSystem& fei_ls, fei::Vector& r);

/** Copy the contents of an fei::Vector to the corresponding field-data
 * locations in a stk::mesh::BulkData instance.
 * This function first calls vec.scatterToOverlap() to ensure that coefficients
 * for shared-entities are available on all sharing processors.
 */
void copy_vector_to_mesh( fei::Vector & vec,
                          const DofMapper & dof,
                          stk::mesh::BulkData & mesh_bulk_data
                        );

/** Scale matrix by a scalar: matrix = scalar*matrix
 */
void scale_matrix(double scalar, fei::Matrix& matrix);

/** Add a scaled matrix to another: dest += scalar*src
 */
void add_matrix_to_matrix(double scalar,
                          const fei::Matrix& src_matrix,
                          fei::Matrix& dest_matrix);

/** Scale vector by a scalar: vec = scalar*vec
 */
void scale_vector(double scalar,
                  fei::Vector& vec);

/** Add a scaled vector to another: dest += scalar*src
 */
void add_vector_to_vector(double scalar,
                          const fei::Vector& src_vector,
                          fei::Vector& dest_vector);

}//namespace linsys
}//namespace stk

#endif

