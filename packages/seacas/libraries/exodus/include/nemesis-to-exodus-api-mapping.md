\page nemesis-mapping Mapping of nemesis API functions to Exodus API functions

\section initial Initial Information Routines

Nemesis API	      |      Exodus API
----------------------|------------------------------
ne_get_init_info      |      ex_get_init_info()
ne_put_init_info      |      ex_put_init_info()
ne_get_init_global    |      ex_get_init_global()
ne_put_init_global    |      ex_put_init_global()
ne_put_version        |      ex_put_nemesis_version()


\section lb Loadbalance Parameter Routines

Nemesis API	      |      Exodus API
----------------------|------------------------------
ne_get_loadbal_param  |       ex_get_loadbal_param()
ne_put_loadbal_param  |       ex_put_loadbal_param()
ne_put_loadbal_param_cc |     ex_put_loadbal_param_cc()

\section param Nodeset, Sideset & Element Block Global Parameter Routines

Nemesis API	        |      Exodus API
------------------------|------------------------------
ne_get_ns_param_global  |     ex_get_ns_param_global()
ne_put_ns_param_global  |     ex_put_ns_param_global()
ne_get_ss_param_global  |     ex_get_ss_param_global()
ne_put_ss_param_global  |     ex_put_ss_param_global()
ne_get_eb_info_global   |     ex_get_eb_info_global()
ne_put_eb_info_global   |     ex_put_eb_info_global()

\section subset Nodeset, Sideset & Element Block Subset Routines

Nemesis API	        |      Exodus API
------------------------|------------------------------
ne_get_n_side_set       |       ex_get_partial_set()
ne_put_n_side_set       |     ex_put_partial_set()
ne_get_n_side_set_df    |     ex_get_partial_set_dist_fact()
ne_put_n_side_set_df    |     ex_put_partial_set_dist_fact()
ne_get_n_node_set       |     ex_get_partial_set()
ne_put_n_node_set       |     ex_put_partial_set()
ne_get_n_node_set_df    |     ex_get_partial_set_dist_fact()
ne_put_n_node_set_df    |     ex_put_partial_set_dist_fact()
ne_get_n_coord          |     ex_get_partial_coord()
ne_put_n_coord          |     ex_put_partial_coord()
ne_get_n_elem_conn      |     ex_get_partial_conn()
ne_put_n_elem_conn      |     ex_put_partial_conn()
ne_get_n_elem_attr      |     ex_get_partial_attr()
ne_put_n_elem_attr      |     ex_put_partial_attr()
ne_get_elem_type        |     ex_get_elem_type()

\section variable Variable Routines

Nemesis API	        |      Exodus API
------------------------|------------------------------
ne_get_n_elem_var       |     ex_get_partial_var()
ne_put_elem_var_slab    |     ex_put_partial_var()
ne_get_n_nodal_var      |     ex_get_partial_var()
ne_put_nodal_var_slab   |     ex_put_partial_var()

\section map Number Map Routines

Nemesis API	        |      Exodus API
------------------------|------------------------------
ne_get_n_elem_num_map   |     ex_get_partial_id_map()
ne_put_n_elem_num_map   |     ex_put_partial_id_map()
ne_get_n_node_num_map   |     ex_get_partial_id_map()
ne_put_n_node_num_map   |     ex_put_partial_id_map()
ne_get_node_map         |     ex_get_processor_node_maps()
ne_put_node_map         |     ex_put_processor_node_maps()
ne_get_elem_map         |     ex_get_processor_elem_maps()
ne_put_elem_map         |     ex_put_processor_elem_maps()

\section comm  Communications Maps Routines

Nemesis API	      |      Exodus API
----------------------|------------------------------
ne_get_cmap_params     |      ex_get_cmap_params()
ne_put_cmap_params      |     ex_put_cmap_params()
ne_put_cmap_params_cc   |     ex_put_cmap_params_cc()
ne_get_node_cmap        |     ex_get_node_cmap()
ne_put_node_cmap        |     ex_put_node_cmap()
ne_get_elem_cmap        |     ex_get_elem_cmap()
ne_put_elem_cmap        |     ex_put_elem_cmap()
ne_get_idx              |     ex_get_idx()

