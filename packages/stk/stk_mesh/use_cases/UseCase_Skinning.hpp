/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_UseCase_Skinning_hpp
#define stk_mesh_UseCase_Skinning_hpp

#include <vector>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>


bool skinning_use_case_1(stk::ParallelMachine pm);
bool skinning_use_case_1b(stk::ParallelMachine pm);
bool skinning_use_case_2(stk::ParallelMachine pm);


void separate_and_skin_mesh(
    stk::mesh::fem::FEMMetaData & fem_meta,
    stk::mesh::BulkData & mesh,
    stk::mesh::Part     & skin_part,
    std::vector< stk::mesh::EntityId > elements_to_separate,
    const stk::mesh::EntityRank rank_of_element
    );

#endif
