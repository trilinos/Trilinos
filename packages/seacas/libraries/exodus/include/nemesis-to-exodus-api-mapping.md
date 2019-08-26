\page nemesis-mapping Mapping of nemesis API functions to Exodus API functions
The nemesis library was originally an extension to the ExodusII
library which provided routines required to support use of Exodus
databases in a parallel setting; typically with a file-per-processor
usage.

Since the use of Exodus in parallel executions is now very common, the
Nemesis library routines have been integrated into the Exodus library
API. In most cases, the exodus API function corresponding to a nemesis
API function is obtained by replacing the `ne_` prefix with an `ex_`
prefix.  There are a few routines where this results in a name
collision or confusion (e.g. ne_put_version() is ex_put_nemesis_version() since
it would be confusing to call it ex_put_version()).  The partial read/write
functions which in nemesis are indicated by a `_n_` in the function
name have been in replaced by `_partial_` (although the corresponding
`ex_*_n_*` function does exist in the deprecated functions).

The tables below list all Nemesis API functions and the corresponding
Exodus API function. In many cases, the only change needed is
replacing `ne_` by `ex_`, but the routines which were made more
"generic" (e.g. ne_get_n_side_set() and ne_get_n_node_set() directly
map to ex_get_n_side_set() and ex_get_n_node_set() which are
deprecated, so the table below shows the recommended
ex_get_partial_set()) additional arguments are required.

The nemesis library can still be used since it is still built upon
request and its implementation is simply wrapper routines which
forward all existing nemesis function calls to the appropriate Exodus
function call (with needed argument changes).

\section initial Initial Information Routines

Nemesis API             |      Exodus API
------------------------|------------------------------
ne_get_init_info        |     ex_get_init_info()
ne_put_init_info        |     ex_put_init_info()
ne_get_init_global      |     ex_get_init_global()
ne_put_init_global      |     ex_put_init_global()
ne_put_version          |     ex_put_nemesis_version()

\section lb Loadbalance Parameter Routines

Nemesis API             |      Exodus API
------------------------|------------------------------
ne_get_loadbal_param    |     ex_get_loadbal_param()
ne_put_loadbal_param    |     ex_put_loadbal_param()
ne_put_loadbal_param_cc |     ex_put_loadbal_param_cc()

\section param Nodeset, Sideset & Element Block Global Parameter Routines

Nemesis API             |      Exodus API
------------------------|------------------------------
ne_get_ns_param_global  |     ex_get_ns_param_global()
ne_put_ns_param_global  |     ex_put_ns_param_global()
ne_get_ss_param_global  |     ex_get_ss_param_global()
ne_put_ss_param_global  |     ex_put_ss_param_global()
ne_get_eb_info_global   |     ex_get_eb_info_global()
ne_put_eb_info_global   |     ex_put_eb_info_global()

\section subset Nodeset, Sideset & Element Block Subset Routines

Nemesis API             |      Exodus API
------------------------|------------------------------
ne_get_n_side_set       |     ex_get_partial_set()
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

Nemesis API             |      Exodus API
------------------------|------------------------------
ne_get_n_elem_var       |     ex_get_partial_var()
ne_put_elem_var_slab    |     ex_put_partial_var()
ne_get_n_nodal_var      |     ex_get_partial_var()
ne_put_nodal_var_slab   |     ex_put_partial_var()

\section map Number Map Routines

Nemesis API             |      Exodus API
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

Nemesis API             |      Exodus API
------------------------|------------------------------
ne_get_cmap_params      |     ex_get_cmap_params()
ne_put_cmap_params      |     ex_put_cmap_params()
ne_put_cmap_params_cc   |     ex_put_cmap_params_cc()
ne_get_node_cmap        |     ex_get_node_cmap()
ne_put_node_cmap        |     ex_put_node_cmap()
ne_get_elem_cmap        |     ex_get_elem_cmap()
ne_put_elem_cmap        |     ex_put_elem_cmap()
ne_get_idx              |     ex_get_idx()
