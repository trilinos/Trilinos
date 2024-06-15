# CHANGELOG

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

