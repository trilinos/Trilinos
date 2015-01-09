/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_linsys_unit_tests_UnitTest_helpers_hpp
#define stk_linsys_unit_tests_UnitTest_helpers_hpp

#include <stk_linsys/LinearSystemInterface.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>

typedef stk_classic::mesh::Field<double>                          ScalarField ;
typedef stk_classic::mesh::Field<double, stk_classic::mesh::Cartesian>    VectorField ;

void fill_utest_mesh_meta_data(stk_classic::mesh::fem::FEMMetaData& fem_meta, bool use_temperature=true);
void fill_utest_mesh_bulk_data(stk_classic::mesh::BulkData& bulk_data);

void assemble_elem_matrices_and_vectors(stk_classic::mesh::BulkData& mesh, ScalarField& field, stk_classic::linsys::LinearSystemInterface& ls);

void assemble_elem_matrices_and_vectors(stk_classic::mesh::BulkData& mesh, ScalarField& field, stk_classic::linsys::DofMapper& dof_mapper, fei::Matrix& matrix, fei::Vector& rhs);

#endif

