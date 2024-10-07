# CHANGELOG

5.21.5-2 (STK_VERSION 5210502) 10/07/2024
  stk_search: Fixed HIP sort error.
  stk_mesh: add multi-field NGP-FieldBLAS field_fill
5.21.5-1 (STK_VERSION 5210501) 9/27/2024
  stk_mesh: deprecate BulkData::relation_exist

5.21.5 (STK_VERSION 5210500) 9/25/2024
   general: Fixed MI300A unified memory build errors.
   stk_search: Turned off sorted results by default.

5.21.4-1 (STK_VERSION 5210401) 9/04/2024
   Fix cmake configuration errors that occurred on AMD MI300A platform

5.21.4 (STK_VERSION 5210400) 8/29/2024
   minor fixes, no signficant API changes

5.21.3-1 (STK_VERSION 5210301) 8/19/2024
   stk_mesh: fix ~65K limitation on per-bucket size of upward-connectivity tables
    - This is an implementation detail that only arises if user sets large bucket
      capacities such that 'entities-per-bucket'*'avg-upward-connections-per-entity' > 65K.
      By default buckets are capped at 512 and this limit won't be reached in practice.
      - Was using 16bit index type, now defaults to 32bit index type, but is settable
        via cmake-settable option 'STK_ENABLE_16BIT_UPWARDCONN_INDEX_TYPE=ON/OFF'

5.21.3 (STK_VERSION 5210300) 8/12/2024
   general: compile-warnings/errors fixed for gcc 12 and arm
   stk_mesh: BulkData::change_entity_owner now returns a bool
    - indicates whether any entities actually changed owner.
   stk_mesh: continue "simple field" transition
    - remove some previously-deprecated functions/behaviors
    - 'use_simple_field()' calls now unnecessary, will be deprecated in future
   stk_mesh: deleted previously-deprecated types/methods
    - DeviceMeshIndex and usage in DeviceMesh/HostMesh
    - NgpFieldBase::rotate_multistate_data()
    - DeviceMesh::host_get_entity(..)
   stk_search: changing template parameters for input Views
    - now less restrictive on how users' views were declared

5.21.1 (STK_VERSION 5210100) 6/17/2024
   stk_mesh: Deprecated BulkData::find_permutation and BulkData::check_permutation
             (replaced by free-functions in FindPermutation.hpp)
   stk_mesh: Added destroy_relations free-function in DestroyRelations.hpp
   stk_mesh: Added overloads of field_fill, field_copy, field_axpbyz in NgpFieldBLAS.hpp

5.19.4 (STK_VERSION 5190401) 5/29/2024
   stk_search: fixed bug in morton: (accessing device view on host)
   stk_search: fixed implementations to respect execution-space
   stk_mesh: change default bucket capacity from 512 to 16 (now grows as needed)

5.19.4 (STK_VERSION 5190400) 5/18/2024
   Added field_fill with Selector in stk_mesh/base/NgpFieldBLAS.hpp
   Added field_copy with and without Selector in stk_mesh/base/NgpFieldBLAS.hpp
   More work on stk-search coarse_search APIs and implementation correctness.
   Minor cleanups in cmake support.
   Fixed stk_search KDTREE OpenMP build

5.19.2 (STK_VERSION 5190201) 4/22/2024
   Fixed broken stk_search install targets

5.19.2 (STK_VERSION 5190200) 4/17/2024
   Added bool argument to stk::mesh::BulkData::update_field_data_states()
   Initial creation of STK_VERSION macro

