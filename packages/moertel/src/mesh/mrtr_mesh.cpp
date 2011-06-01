#include "mrtr_mesh.H"
#include <iostream>
#include <algorithm>
#include <sstream>
#include <string>


#ifdef HAVE_MOERTEL_NEMESIS
#include "ne_nemesisI.h"
#endif

/*
Public interface to the Mesh class
*/

#ifdef HAVE_MOERTEL_NEMESIS
MOERTEL::Mesh::Mesh( const int pid, const int np, const bool v ):
    proc_id(pid), nprocs(np), verbose(v) {

    nodal_names_written = 0;
    elem_names_written = 0;
}
// Invalid constructor when running in parallel
//
MOERTEL::Mesh::Mesh(){

	std::stringstream outp;
			outp << "Must specify processor ID and number of processors in MOERTEL::Mesh constructor" 
				 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	throw ReportError(outp);

}

#else

MOERTEL::Mesh::Mesh( const int pid, const bool v ):
    verbose(v) {

    nodal_names_written = 0;
    elem_names_written = 0;
}
#endif

MOERTEL::Mesh::~Mesh() {}

int MOERTEL::Mesh::read_exodus(const char * filename) {

    int comp_ws = sizeof(double);//cn send this to exodus to tell it we are using doubles
    int io_ws = 0;
    std::vector<int>::iterator a;

    int ex_id = ex_open(filename,//cn open file
                        EX_READ,
                        &comp_ws,
                        &io_ws,
                        &exodus_version);

    if(ex_id < 0) {

		std::stringstream outp;
			outp << "Cannot open file" << filename
				 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
		throw ReportError(outp);

    }

    char _title[MAX_LINE_LENGTH];

    int ex_err = 0;

    ex_err = ex_get_init(ex_id,//cn read header
                         _title,
                         &num_dim,
                         &num_nodes,
                         &num_elem,
                         &num_elem_blk,
                         &num_node_sets,
                         &num_side_sets);

    check_exodus_error(ex_err,"MOERTEL::Mesh::read_exodus ex_get_init");

    title = _title;

    if(verbose)

        std::cout<<"=== ExodusII Read Info ==="<<std::endl
                 <<" ExodusII "<<exodus_version<<std::endl
                 <<" File "<<filename<<std::endl
                 <<" Exodus ID "<<ex_id<<std::endl
                 <<" comp_ws "<<comp_ws<<std::endl
                 <<" io_ws "<<io_ws<<std::endl
                 <<" Title "<<title<<std::endl
                 <<" num_dim "<<num_dim<<std::endl
                 <<" num_nodes "<<num_nodes<<std::endl
                 <<" num_elem "<<num_elem<<std::endl
                 <<" num_elem_blk "<<num_elem_blk<<std::endl
                 <<" num_node_sets "<<num_node_sets<<std::endl
                 <<" num_side_sets "<<num_side_sets<<std::endl<<std::endl;

    x.resize(num_nodes);
    y.resize(num_nodes);
    z.resize(num_nodes);

#ifdef HAVE_MOERTEL_NEMESIS

    if ( num_dim == 2 ) {

        ne_get_n_coord(ex_id, 1, num_nodes,
                       &x[0],
                       &y[0],
                       NULL);

        std::fill(z.begin(), z.end(), 0);
    } else {

        ne_get_n_coord(ex_id, 1, num_nodes,
                       &x[0],
                       &y[0],
                       &z[0]);
    }

    int num_proc;


    ne_get_init_info(ex_id, &num_proc, &nprocs_infile, &filetype);

    if(num_proc != nprocs) {

		std::stringstream outp;
			outp << "ERROR in file read: number of processors does not match number of input files"
				 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
		throw ReportError(outp);

    }


#else

    if ( num_dim == 2 ) {

        ex_err = ex_get_coord(ex_id,
                              &x[0],
                              &y[0],
                              0);

        std::fill(z.begin(), z.end(), 0);
    } else {

        ex_err = ex_get_coord(ex_id,
                              &x[0],
                              &y[0],
                              &z[0]);
    }

    check_exodus_error(ex_err,"MOERTEL::Mesh::read_exodus ex_get_coord");

#endif

    blk_ids.resize(num_elem_blk);
    num_elem_in_blk.resize(num_elem_blk);
    num_node_per_elem_in_blk.resize(num_elem_blk);
    connect.resize(num_elem_blk);

    ex_err = ex_get_elem_blk_ids(ex_id,//cn read block ids
                                 &blk_ids[0]);

    check_exodus_error(ex_err,"MOERTEL::Mesh::read_exodus ex_get_elem_blk_ids");

    for (int i = 0; i < num_elem_blk; i++) { //cn loop over blocks   we should have block ids start at 0

        char elem_type[MAX_STR_LENGTH];
        int num_attr = 0;

        ex_err = ex_get_elem_block(ex_id,//cn read elems for this blk[i]
                                   blk_ids[i],
                                   elem_type,
                                   &num_elem_in_blk[i],
                                   &num_node_per_elem_in_blk[i],
                                   &num_attr);

        check_exodus_error(ex_err,"MOERTEL::Mesh::read_exodus ex_get_elem_block");

        blk_elem_type.push_back(elem_type);

        if(verbose)

            std::cout<<" +++ Element Block Info +++"<<std::endl
                     <<"  i "<<i<<std::endl
                     <<"  block_ids[i] "<<blk_ids[i]<<std::endl
                     <<"  elem_type "<<blk_elem_type[i]<<std::endl
                     <<"  num_elem_this_blk "<<num_elem_in_blk[i]<<std::endl
                     <<"  num_node_per_elem "<<num_node_per_elem_in_blk[i]<<std::endl
                     <<"  num_attr "<<num_attr<<std::endl<<std::endl;

        //cn get connectivity info here

        connect[i].resize(num_elem_in_blk[i]*num_node_per_elem_in_blk[i]);

#ifdef HAVE_MOERTEL_NEMESIS

        ne_get_n_elem_conn(ex_id, blk_ids[i], 1, num_elem_in_blk[i], &connect[i][0]);

#else

        ex_err = ex_get_elem_conn(ex_id,
                                  blk_ids[i],
                                  &connect[i][0]);

        check_exodus_error(ex_err,"MOERTEL::Mesh::read_exodus ex_get_elem_conn");

#endif

        for(a = connect[i].begin(); a != connect[i].end(); a++)  // fix FORTRAN indexing

            (*a)--;

    }

    if(num_side_sets > 0) {

        ss_ids.resize(num_side_sets);
        num_sides_per_ss.resize(num_side_sets);
        num_df_per_ss.resize(num_side_sets);
        ss_ctr_list.resize(num_side_sets);
        ss_node_list.resize(num_side_sets);
        ss_elem_list.resize(num_side_sets);
        ss_side_list.resize(num_side_sets);

        side_set_node_map.resize(num_nodes);

        std::fill(side_set_node_map.begin(), side_set_node_map.end(), -1); // initialize it

        ex_err = ex_get_side_set_ids(ex_id,//cn side set ids
                                     &ss_ids[0]);

        check_exodus_error(ex_err,"MOERTEL::Mesh::read_exodus ex_get_side_set_ids");

        for (int i = 0; i < num_side_sets; i++) { //cn loop over sidesets

            ex_err = ex_get_side_set_param(ex_id,
                                           ss_ids[i],
                                           &num_sides_per_ss[i],
                                           &num_df_per_ss[i]);

            check_exodus_error(ex_err,"MOERTEL::Mesh::read_exodus ex_get_side_set_param");

            ss_ctr_list[i].resize(num_sides_per_ss[i]);
            ss_node_list[i].resize(num_df_per_ss[i]);

            ss_elem_list[i].resize(num_sides_per_ss[i]);
            ss_side_list[i].resize(num_sides_per_ss[i]);

            if(verbose)

                std::cout<<" +++ Sideset Info +++"<<std::endl
                         <<"  i "<<i<<std::endl
                         <<"  ss_ids[i] "<<ss_ids[i]<<std::endl
                         <<"  num_sides_per_set[i] "<<num_sides_per_ss[i]<<std::endl
                         <<"  num_df_per_sideset[i] "<<num_df_per_ss[i]<<std::endl<<std::endl;

            ex_err = ex_get_side_set(ex_id, ss_ids[i], &ss_elem_list[i][0], &ss_side_list[i][0]);

            ex_err = ex_get_side_set_node_list(ex_id, ss_ids[i], &ss_ctr_list[i][0], &ss_node_list[i][0]);

            for(a = ss_node_list[i].begin(); a != ss_node_list[i].end(); a++) { // fix FORTRAN indexing

                (*a)--;

                side_set_node_map[*a] = ss_ids[i];

            }

        } // end loop over side sets


    } // end if sidesets > 0

    if(num_node_sets > 0) {

        ns_ids.resize(num_node_sets);
        num_nodes_per_ns.resize(num_node_sets);
        num_df_per_ns.resize(num_node_sets);
        ns_node_list.resize(num_node_sets);

        node_set_map.resize(num_nodes);

        std::fill(node_set_map.begin(), node_set_map.end(), -1); // initialize it

        ex_err = ex_get_node_set_ids(ex_id,
                                     &ns_ids[0]);

        check_exodus_error(ex_err,"MOERTEL::Mesh::read_exodus ex_get_node_set_ids");

        for (int i = 0; i < num_node_sets; i++) { //cn loop over sidesets

            ex_err = ex_get_node_set_param(ex_id,
                                           ns_ids[i],
                                           &num_nodes_per_ns[i],
                                           &num_df_per_ns[i]);

            check_exodus_error(ex_err,"MOERTEL::Mesh::read_exodus ex_get_node_set_param");

            ns_node_list[i].resize(num_nodes_per_ns[i]);

            if(verbose)

                std::cout<<" +++ Nodeset Info +++"<<std::endl
                         <<"  i "<<i<<std::endl
                         <<"  ns_ids[i] "<<ns_ids[i]<<std::endl
                         <<"  num_nodes_per_set[i] "<<num_nodes_per_ns[i]<<std::endl
                         <<"  num_df_per_sideset[i] "<<num_df_per_ns[i]<<std::endl<<std::endl;


//#ifdef HAVE_MOERTEL_NEMESIS

//        ex_err = ne_get_n_node_set(ex_id, ns_ids[i], 1, num_nodes_per_ns[i], &ns_node_list[i][0]);

//#else

            ex_err = ex_get_node_set(ex_id, ns_ids[i], &ns_node_list[i][0]);

//#endif

            for(a = ns_node_list[i].begin(); a != ns_node_list[i].end(); a++) {

                // fix FORTRAN indexing

                (*a)--;

                // fill the hash table

                node_set_map[*a] = ns_ids[i];

            }

        } // end loop over node sets

    } // end if nodesets > 0

    node_num_map.resize(num_nodes);
    elem_num_map.resize(num_elem);

#ifdef HAVE_MOERTEL_NEMESIS

    ne_get_init_global(ex_id, &ne_num_global_nodes, &ne_num_global_elems, &ne_num_global_elem_blks,
                       &ne_num_global_node_sets, &ne_num_global_side_sets);

    ne_get_loadbal_param(ex_id, &num_internal_nodes, &num_border_nodes, &num_external_nodes,
                         &num_internal_elems, &num_border_elems, &num_node_cmaps, &num_elem_cmaps, proc_id);

    ne_get_n_node_num_map(ex_id, 1, num_nodes, &node_num_map[0]);
    for(a = node_num_map.begin(); a != node_num_map.end(); a++) (*a)--;
    ne_get_n_elem_num_map(ex_id, 1, num_elem, &elem_num_map[0]);
    for(a = elem_num_map.begin(); a != elem_num_map.end(); a++) (*a)--;

    elem_mapi.resize(num_internal_elems);
    elem_mapb.resize(num_border_elems);

    ne_get_elem_map(ex_id, &elem_mapi[0], &elem_mapb[0], proc_id);
    for(a = elem_mapi.begin(); a != elem_mapi.end(); a++) (*a)--;
    for(a = elem_mapb.begin(); a != elem_mapb.end(); a++) (*a)--;

    node_mapi.resize(num_internal_nodes);
    node_mapb.resize(num_border_nodes);
    node_mape.resize(num_external_nodes);

    ne_get_node_map(ex_id, &node_mapi[0], &node_mapb[0], &node_mape[0], proc_id);
    for(a = node_mapi.begin(); a != node_mapi.end(); a++) (*a)--;
    for(a = node_mapb.begin(); a != node_mapb.end(); a++) (*a)--;
    for(a = node_mape.begin(); a != node_mape.end(); a++) (*a)--;

    my_node_num_map = node_mapi;  // nodes this proc is responsible for
    my_node_num_map.insert(my_node_num_map.end(), node_mapb.begin(), node_mapb.end());

    if(ne_num_global_node_sets > 0) {

        global_ns_ids.resize(ne_num_global_node_sets);
        num_global_node_counts.resize(ne_num_global_node_sets);
        num_global_node_df_counts.resize(ne_num_global_node_sets);

        ne_get_ns_param_global(ex_id, &global_ns_ids[0], &num_global_node_counts[0],
                               &num_global_node_df_counts[0]);

    }

    if(ne_num_global_side_sets > 0) {

        global_ss_ids.resize(ne_num_global_side_sets);
        num_global_side_counts.resize(ne_num_global_side_sets);
        num_global_side_df_counts.resize(ne_num_global_side_sets);

        ne_get_ss_param_global(ex_id, &global_ss_ids[0], &num_global_side_counts[0],
                               &num_global_side_df_counts[0]);

    }

    global_elem_blk_ids.resize(ne_num_global_elem_blks);
    global_elem_blk_cnts.resize(ne_num_global_elem_blks);

    ne_get_eb_info_global(ex_id, &global_elem_blk_ids[0], &global_elem_blk_cnts[0]);

    node_cmap_ids.resize(num_node_cmaps);
    node_cmap_node_cnts.resize(num_node_cmaps);

    elem_cmap_ids.resize(num_elem_cmaps);
    elem_cmap_elem_cnts.resize(num_elem_cmaps);

    ne_get_cmap_params(ex_id, &node_cmap_ids[0], &node_cmap_node_cnts[0],
                       &elem_cmap_ids[0], &elem_cmap_elem_cnts[0],
                       proc_id);

    node_ids_in_cmap.resize(num_node_cmaps);
    n_proc_ids_in_cmap.resize(num_node_cmaps);

    for(int i = 0; i < num_node_cmaps; i++) {

        node_ids_in_cmap[i].resize(node_cmap_node_cnts[i]);
        n_proc_ids_in_cmap[i].resize(node_cmap_node_cnts[i]);

        ne_get_node_cmap(ex_id, node_cmap_ids[i], &node_ids_in_cmap[i][0], &n_proc_ids_in_cmap[i][0], proc_id);

    }

    elem_ids_in_cmap.resize(num_elem_cmaps);
    e_side_ids_in_cmap.resize(num_elem_cmaps);
    e_proc_ids_in_cmap.resize(num_elem_cmaps);

    for(int i = 0; i < num_elem_cmaps; i++) {

        elem_ids_in_cmap[i].resize(elem_cmap_elem_cnts[i]);
        e_side_ids_in_cmap[i].resize(elem_cmap_elem_cnts[i]);
        e_proc_ids_in_cmap[i].resize(elem_cmap_elem_cnts[i]);

        ne_get_elem_cmap(ex_id, elem_cmap_ids[i], &elem_ids_in_cmap[i][0],
                         &e_side_ids_in_cmap[i][0], &e_proc_ids_in_cmap[i][0], proc_id);

    }

#if 0

    my_node_num_map = node_mapi;  // start with the nodes internal to this processor

    for(int i = 0; i < num_node_cmaps; i++) {

        std::vector<int> cmap_node_ids(node_cmap_node_cnts[i]);
        std::vector<int> cmap_node_procids(node_cmap_node_cnts[i]);
        ne_get_node_cmap(ex_id, node_cmap_ids[i], &cmap_node_ids[0], &cmap_node_procids[0], proc_id);

        // build the node array for the nodes local to this processor
        //
        // HACK ALERT!!!!
        //
        // Put nodes on this processor that it shares with a lower numbered processor

        for(int j = 0; j < node_cmap_node_cnts[i]; j++) {

            if(cmap_node_procids[j] < proc_id) // add the node to my_node_map

                my_node_num_map.push_back(cmap_node_ids[j]);
        }

    }
#endif


#else

    ex_err = ex_get_node_num_map(ex_id, &node_num_map[0]);
    for(a = node_num_map.begin(); a != node_num_map.end(); a++) (*a)--;
    ex_err = ex_get_map(ex_id, &elem_num_map[0]);
    for(a = elem_num_map.begin(); a != elem_num_map.end(); a++) (*a)--;

    my_node_num_map = node_num_map;  // same in serial

#endif

    ex_err = close_exodus(ex_id);//cn close file

    if(verbose) {

        std::cout << "There are " << num_elem << " elements in this mesh." << std::endl;
        std::cout<<"=== End ExodusII Read Info ==="<<std::endl<<std::endl;

    }


    return 0;

}

/*
Compute nodal adjacency for a standard serial matrix graph
note that we do not store the diagonal--maybe we should
*/

void MOERTEL::Mesh::compute_nodal_adj() {

    nodal_adj.resize(num_nodes, std::vector<int>(0));
    nodal_adj_idx.resize(num_nodes + 1);
    nodal_adj_idx[0] = 0;   //probably wont work in parallel, or need to start somewhere else

    if(verbose)

        std::cout<<"=== void mesh::compute_nodal_adj ==="<<std::endl<<" nodal_adj"<<std::endl<<"  ";

    for (int blk = 0; blk < num_elem_blk; blk++) {

        std::vector<int> temp(num_node_per_elem_in_blk[blk]);

        for(int i = 0; i < num_elem_in_blk[blk]; i++) {

            for (int j = 0; j < num_node_per_elem_in_blk[blk]; j++) {

                temp[j] = connect[blk][i * num_node_per_elem_in_blk[blk] + j]; //load up nodes on each element
                //std::cout<<temp[j]<<std::endl;

            }

            for (int j = 0; j < num_node_per_elem_in_blk[blk]; j++) {

                for(int k = 0; k < num_node_per_elem_in_blk[blk]; k++) {

                    if(temp[j] != temp[k] ) { //cn skip the diagonal and load up nodes

                        nodal_adj[temp[j]].push_back(temp[k]);

                        //std::cout<<temp[j]<<","<<temp[k]<<std::endl;
                    }

                }
            }
        }
    }

    for(int i = 0; i < num_nodes; i++) {

        //cn sort and remove duplicates

        std::sort(nodal_adj[i].begin(), nodal_adj[i].end());

        std::vector<int>::iterator unique_end =

            std::unique(nodal_adj[i].begin(), nodal_adj[i].end());

        nodal_adj[i].erase(unique_end, nodal_adj[i].end());

        nodal_adj_idx[i + 1] = nodal_adj_idx[i] + nodal_adj[i].size();

        for( int j = nodal_adj_idx[i]; j < nodal_adj_idx[i + 1]; j++)

            nodal_adj_array.push_back(nodal_adj[i][j - nodal_adj_idx[i]]);

        if(verbose) {

            for( unsigned int j = 0;  j < nodal_adj[i].size(); j++)

                std::cout<<nodal_adj[i][j]<<" ";

            std::cout<<std::endl<<"  ";

        }
    }

    if(verbose) {

        std::cout<<std::endl<<" nodal_adj_idx"<<std::endl;

        for( int i = 0; i < num_nodes + 1; i++)

            std::cout<<"  "<<nodal_adj_idx[i]<<std::endl;

        std::cout<<std::endl<<" nodal_adj_array"<<std::endl;

        for( int i = 0; i < nodal_adj_idx[num_nodes]; i++)

            std::cout<<"  "<<nodal_adj_array[i]<<std::endl;

        std::cout<<"=== End void mesh::compute_nodal_adj ==="<<std::endl<<std::endl;
    }

}

int MOERTEL::Mesh::get_boundary_status(int blk, int elem) {

    int status;

    for(int i = 0; i < num_node_per_elem_in_blk[blk]; i++)

        if((status = node_set_map[connect[blk][elem * num_node_per_elem_in_blk[blk] + i]]) >= 0)

            return status;

    return -1;

}

int MOERTEL::Mesh::get_node_boundary_status(int nodeid) {

    int status;

    if((status = node_set_map[nodeid]) >= 0)

        return status;

    return -1;

}

/*
Private interface to MOERTEL::Mesh class
*/

int MOERTEL::Mesh::close_exodus(int ex_id) {

    int ex_err = ex_update (ex_id);

    check_exodus_error(ex_err,"MOERTEL::Mesh::close_exodus ex_close");

    ex_err = ex_close(ex_id);

    check_exodus_error(ex_err,"MOERTEL::Mesh::close_exodus ex_close");

    if(verbose)

        std::cout<<"=== ExodusII Close Info ==="<<std::endl
                 <<" Exodus ID "<<ex_id<<std::endl
                 <<"=== End ExodusII Close Info ==="<<std::endl<<std::endl;

    return ex_err;

}

void MOERTEL::Mesh::check_exodus_error(const int ex_err, const std::string msg) {

    if (ex_err < 0)

        std::cout<<"ExodusII error:  "<<msg<<std::endl<<std::endl;

    return;

}

int MOERTEL::Mesh::write_exodus(const char * filename) {

    int ex_id = create_exodus(filename);

    write_nodal_coordinates_exodus(ex_id);
    write_element_blocks_exodus(ex_id);
    write_nodal_data_exodus(ex_id);

    return 0;

}

int MOERTEL::Mesh::open_exodus_write_mesh(const char * filename) {

    staged_ex_id = create_exodus(filename);

    // Write the nodal coordinates of the mesh
    write_nodal_coordinates_exodus(staged_ex_id);
    // Write the mesh elements
    write_element_blocks_exodus(staged_ex_id);

    return 0;

}

void MOERTEL::Mesh::set_time_step(int its, double itime) {

    ts = its;
    time = itime;
}

int MOERTEL::Mesh::update_file() {

    int ex_err = 0;
    char **var_names;

    if(verbose){

        std::cout<<"=== Write Nodal Data Exodus ==="<<std::endl
                 <<" num_nodal_fields "<<nodal_fields.size()<<std::endl;

        std::cout<<"=== Write Element Data Exodus ==="<<std::endl
                 <<" num_elem_fields "<<elem_fields.size()<<std::endl;

	}

    if(nodal_fields.size() == 0 && elem_fields.size() == 0) return 0;  // Bail if there is nothing to write

    ex_err = ex_put_time (staged_ex_id, ts + 1, &time);


	if(nodal_names_written == 0) {

        ex_err = ex_put_var_param (staged_ex_id, "N", nodal_fields.size());

        var_names = new char*[nodal_fields.size()];

        for(unsigned int i = 0; i < nodal_fields.size(); i++) {

            var_names[i] = (char *)&nodal_field_names[i][0];

            if(verbose)

				std::cout<<" nodal var name  "<<var_names[i]<<std::endl<<std::endl;

        }

        ex_err = ex_put_var_names (staged_ex_id, "N", nodal_fields.size(), var_names);

        delete [] var_names;

		nodal_names_written = 1;

    }

	if(elem_names_written == 0) {

		ex_err = ex_put_var_param (staged_ex_id, "E", elem_fields.size());

		var_names = new char*[elem_fields.size()];

		for(unsigned int i = 0; i < elem_fields.size(); i++) {

			var_names[i] = (char *)&elem_field_names[i][0];

			if(verbose)

				std::cout<<" elem var name  "<<var_names[i]<<std::endl<<std::endl;

		}

		ex_err = ex_put_var_names (staged_ex_id, "E", elem_fields.size(), var_names);

		delete [] var_names;

		elem_names_written = 1;

	}

	for(unsigned int i = 0; i < nodal_fields.size(); i++){

        ex_err = ex_put_nodal_var (staged_ex_id, ts + 1, i + 1, num_nodes, &nodal_fields[i][0]);

	}

	int cnt = 0;

	for(unsigned int i = 0; i < elem_fields.size(); i++){
		for(int j = 0; j < num_elem_blk; j++){

			ex_err = ex_put_elem_var (staged_ex_id, ts + 1, i + 1, blk_ids[j], num_elem_in_blk[j], 
					&elem_fields[i][cnt]);

			cnt += num_elem_in_blk[j];
		}

	}


    ex_err = ex_update(staged_ex_id);

    nodal_fields.clear();
    nodal_field_names.clear();
	elem_fields.clear();
	elem_field_names.clear();

    return ex_err;

}

int MOERTEL::Mesh::write_nodal_coordinates_exodus(int ex_id) {

    char ** var_names;
    int ex_err;
    std::vector<int> tmpvec, tmpvec1, tmpvec2;
    std::vector<int>::iterator a;

#ifdef HAVE_MOERTEL_NEMESIS

    ne_put_init_global(ex_id, ne_num_global_nodes, ne_num_global_elems, ne_num_global_elem_blks,
                       ne_num_global_node_sets, ne_num_global_side_sets);

    ne_put_init_info(ex_id, nprocs, nprocs_infile, &filetype);

    ne_put_loadbal_param(ex_id, num_internal_nodes, num_border_nodes, num_external_nodes,
                         num_internal_elems, num_border_elems, num_node_cmaps, num_elem_cmaps, proc_id);

    tmpvec = node_num_map;
    for(a = tmpvec.begin(); a != tmpvec.end(); a++) (*a)++;
    ne_put_n_node_num_map(ex_id, 1, num_nodes, &tmpvec[0]);

    tmpvec = elem_num_map;
    for(a = tmpvec.begin(); a != tmpvec.end(); a++) (*a)++;
    ne_put_n_elem_num_map(ex_id, 1, num_elem, &tmpvec[0]);

    tmpvec = elem_mapi;
    for(a = tmpvec.begin(); a != tmpvec.end(); a++) (*a)++;
    tmpvec1 = elem_mapb;
    for(a = tmpvec1.begin(); a != tmpvec1.end(); a++) (*a)++;
    ne_put_elem_map(ex_id, &tmpvec[0], &tmpvec1[0], proc_id);

    tmpvec = node_mapi;
    for(a = tmpvec.begin(); a != tmpvec.end(); a++) (*a)++;
    tmpvec1 = node_mapb;
    for(a = tmpvec1.begin(); a != tmpvec1.end(); a++) (*a)++;
    tmpvec2 = node_mape;
    for(a = tmpvec2.begin(); a != tmpvec2.end(); a++) (*a)++;
    ne_put_node_map(ex_id, &tmpvec[0], &tmpvec1[0], &tmpvec2[0], proc_id);

    ne_put_n_coord(ex_id, 1, num_nodes, &x[0], &y[0], &z[0]);

    if(ne_num_global_node_sets > 0)

        ne_put_ns_param_global(ex_id, &global_ns_ids[0], &num_global_node_counts[0],
                               &num_global_node_df_counts[0]);

    if(ne_num_global_side_sets > 0)

        ne_put_ss_param_global(ex_id, &global_ss_ids[0], &num_global_side_counts[0],
                               &num_global_side_df_counts[0]);

    ne_put_eb_info_global(ex_id, &global_elem_blk_ids[0], &global_elem_blk_cnts[0]);

    ne_put_cmap_params(ex_id, &node_cmap_ids[0], &node_cmap_node_cnts[0],
                       &elem_cmap_ids[0], &elem_cmap_elem_cnts[0],
                       proc_id);

    for(int i = 0; i < num_node_cmaps; i++) {

        ne_put_node_cmap(ex_id, node_cmap_ids[i], &node_ids_in_cmap[i][0], &n_proc_ids_in_cmap[i][0], proc_id);

    }

    for(int i = 0; i < num_elem_cmaps; i++) {

        ne_put_elem_cmap(ex_id, elem_cmap_ids[i], &elem_ids_in_cmap[i][0],
                         &e_side_ids_in_cmap[i][0], &e_proc_ids_in_cmap[i][0], proc_id);

    }


#else

    ex_err = ex_put_coord(ex_id, &x[0], &y[0], &z[0]);

#endif

    if(num_node_sets > 0) {

//    ex_err = ex_put_node_set_ids(ex_id,
//  				 &ns_ids[0]);

        for (int i = 0; i < num_node_sets; i++) { //cn loop over sidesets

            ex_err = ex_put_node_set_param(ex_id,
                                           ns_ids[i],
                                           num_nodes_per_ns[i],
                                           num_df_per_ns[i]);

            tmpvec = ns_node_list[i];
            for(a = tmpvec.begin(); a != tmpvec.end(); a++) (*a)++;

#ifdef HAVE_MOERTEL_NEMESIS

            ex_err = ne_put_n_node_set(ex_id, ns_ids[i], 1, num_nodes_per_ns[i], &tmpvec[0]);

#else

            ex_err = ex_put_node_set(ex_id, ns_ids[i], &tmpvec[0]);

#endif

        }

    } // end if nodesets > 0

    if(num_side_sets > 0) {

//    ex_err = ex_get_side_set_ids(ex_id,//cn side set ids
//				 &ss_ids[0]);


        for (int i = 0; i < num_side_sets; i++) { //cn loop over sidesets

            ex_err = ex_put_side_set_param(ex_id,
                                           ss_ids[i],
                                           num_sides_per_ss[i],
                                           num_df_per_ss[i]);

            ex_err = ex_put_side_set(ex_id, ss_ids[i], &ss_elem_list[i][0], &ss_side_list[i][0]);

//	tmpvec = ss_node_list[i];
//        for(a = tmpvec.begin(); a != tmpvec.end(); a++) (*a)++;

//        ex_err = ex_put_side_set_node_list(ex_id, ss_ids[i], &ss_ctr_list[i][0], &tmpvec[0]);

        } // end loop over side sets

    } // end if sidesets > 0

    var_names = new char*[3];
    var_names[0] = new char[4];
    var_names[1] = new char[4];
    var_names[2] = new char[4];

    strcpy(var_names[0], "\"x\"");
    strcpy(var_names[1], "\"y\"");
    strcpy(var_names[2], "\"z\"");

    ex_err = ex_put_coord_names(ex_id, var_names);

    delete [] var_names[0];
    delete [] var_names[1];
    delete [] var_names[2];
    delete [] var_names;


    return ex_err;

}

int MOERTEL::Mesh::write_element_blocks_exodus(int ex_id) {

    int ex_err = 0;

    for(int i = 0; i < num_elem_blk; i++) {

        ex_err = ex_put_elem_block(ex_id, blk_ids[i], &blk_elem_type[i][0], num_elem_in_blk[i],
                                   num_node_per_elem_in_blk[i], 0);

        std::vector<int> connect_tmp(num_node_per_elem_in_blk[i] * num_elem_in_blk[i]);

        for ( int j = 0; j < num_node_per_elem_in_blk[i] * num_elem_in_blk[i]; j++ )

            connect_tmp[j] = connect[i][j] + 1;


#ifdef HAVE_MOERTEL_NEMESIS

        ne_put_n_elem_conn(ex_id, blk_ids[i], 1, num_elem_in_blk[i], &connect_tmp[0]);

#else

        ex_err = ex_put_elem_conn(ex_id, blk_ids[i], &connect_tmp[0]);

#endif
    }


    return ex_err;

}

int MOERTEL::Mesh::write_nodal_data_exodus(int ex_id) {

    int ex_err;
    char **var_names;


    if(verbose)

        std::cout<<"=== Write Nodal Data Exodus ==="<<std::endl
                 <<" num_nodal_fields "<<nodal_fields.size()<<std::endl;

    if(nodal_fields.size() == 0) return 0;

    ex_err = ex_put_var_param (ex_id, "N", nodal_fields.size());

    var_names = new char*[nodal_fields.size()];

    for(unsigned int i = 0; i < nodal_fields.size(); i++) {

        var_names[i] = (char *)&nodal_field_names[i][0];

        if(verbose)

            std::cout<<" name  "<<var_names[i]<<std::endl<<std::endl;

    }

    ex_err = ex_put_var_names (ex_id, "N", nodal_fields.size(), var_names);

    for(unsigned int i = 0; i < nodal_fields.size(); i++)

        ex_err = ex_put_nodal_var (ex_id, 1, i + 1, num_nodes, &nodal_fields[i][0]);


    delete [] var_names;

    return ex_err;

}

int MOERTEL::Mesh::write_elem_data_exodus(int ex_id) {

    int ex_err;
    char **var_names;


    if(verbose)

        std::cout<<"=== Write Elem Data Exodus ==="<<std::endl
                 <<" num_elem_fields "<<elem_fields.size()<<std::endl;

    if(elem_fields.size() == 0) return 0;

    ex_err = ex_put_var_param (ex_id, "E", elem_fields.size());

    var_names = new char*[elem_fields.size()];

    for(unsigned int i = 0; i < elem_fields.size(); i++) {

        var_names[i] = (char *)&elem_field_names[i][0];

        if(verbose)

            std::cout<<" elem field name  "<<var_names[i]<<std::endl<<std::endl;

    }

    ex_err = ex_put_var_names (ex_id, "E", elem_fields.size(), var_names);

	int cnt = 0;

	for(unsigned int i = 0; i < elem_fields.size(); i++){
		for(int j = 0; j < num_elem_blk; j++){

			ex_err = ex_put_elem_var (staged_ex_id, 1, i + 1, blk_ids[j], num_elem_in_blk[j], 
					&elem_fields[i][cnt]);

			cnt += num_elem_in_blk[j];
		}

	}


    delete [] var_names;

    return ex_err;

}

/*
int Mesh::add_nodal_data(std::string &name, std::vector<double> &data){

   if(num_nodes != data.size()){

		std::stringstream oss;
      oss	<<"ERROR in add_node_data: node data field differs in length from the mesh";
	  throw ReportError(oss);

   }

   num_nodal_fields++;

//   nodal_field_names.push_back(name);
//   nodal_fields.push_back(data);

   return 1;

}
*/
//int Mesh::add_nodal_data(std::basic_string<char, std::char_traits<char> > name, double *data){return 1;}
int MOERTEL::Mesh::add_nodal_data(std::string name, std::vector<double> &data) {


    nodal_field_names.push_back(name);
    nodal_fields.push_back(data);

    if(verbose)

        std::cout<<"=== Add nodal fields ==="<<std::endl
                 <<" num_nodal_fields "<<nodal_fields.size()<<std::endl
                 <<" sizeof nodal_field_names "<<nodal_field_names.size()<<std::endl
                 <<" sizeof nodal_fields "<<nodal_fields.size()<<std::endl<<std::endl;

    return 1;

}

int MOERTEL::Mesh::add_nodal_data(std::string name, double *data) {


    std::vector<double> a(data, data + num_nodes);;

    nodal_field_names.push_back(name);
    nodal_fields.push_back(a);

    if(verbose)

        std::cout<<"=== Add nodal fields ==="<<std::endl
                 <<" num_nodal_fields "<<nodal_fields.size()<<std::endl
                 <<" sizeof nodal_field_names "<<nodal_field_names.size()<<std::endl
                 <<" sizeof nodal_fields "<<nodal_fields.size()<<std::endl<<std::endl;

    return 1;

}

int MOERTEL::Mesh::add_elem_data(std::string name, std::vector<double> &data) {


    elem_field_names.push_back(name);
    elem_fields.push_back(data);

    if(verbose)

        std::cout<<"=== Add elem fields ==="<<std::endl
                 <<" num_elem_fields "<<elem_fields.size()<<std::endl
                 <<" sizeof elem_field_names "<<elem_field_names.size()<<std::endl
                 <<" sizeof elem_fields "<<elem_fields.size()<<std::endl<<std::endl;

    return 1;

}

int MOERTEL::Mesh::add_elem_data(std::string name, double *data) {


    std::vector<double> a(data, data + num_nodes);;

    elem_field_names.push_back(name);
    elem_fields.push_back(a);

    if(verbose)

        std::cout<<"=== Add elem fields ==="<<std::endl
                 <<" num_elem_fields "<<elem_fields.size()<<std::endl
                 <<" sizeof elem_field_names "<<elem_field_names.size()<<std::endl
                 <<" sizeof elem_fields "<<elem_fields.size()<<std::endl<<std::endl;

    return 1;

}


int MOERTEL::Mesh::create_exodus(const char * filename) {

    //Store things as doubles
    int comp_ws = sizeof(double);// = 8
    int io_ws = sizeof(double);// = 8

    int ex_id = ex_create(filename, EX_CLOBBER, &comp_ws, &io_ws);

    ex_id = ex_open(filename,
                    EX_WRITE,
                    &comp_ws,
                    &io_ws,
                    &exodus_version);

    if(verbose)

        std::cout<<"=== ExodusII Create Info ==="<<std::endl
                 <<" ExodusII "<<exodus_version<<std::endl
                 <<" File "<<filename<<std::endl
                 <<" Exodus ID "<<ex_id<<std::endl
                 <<" comp_ws "<<comp_ws<<std::endl
                 <<" io_ws "<<io_ws<<std::endl;

    char * title = new char[16];

    strcpy(title, "\"Exodus output\"");

    int ex_err = ex_put_init(ex_id, title, num_dim,
                             num_nodes, num_elem, num_elem_blk,
                             num_node_sets, num_side_sets);

    check_exodus_error(ex_err,"MOERTEL::Mesh::create_exodus ex_put_init");

    if(verbose)

        std::cout<<" Title "<<title<<std::endl
                 <<" num_dim "<<num_dim<<std::endl
                 <<" num_nodes "<<num_nodes<<std::endl
                 <<" num_elem "<<num_elem<<std::endl
                 <<" num_elem_blk "<<num_elem_blk<<std::endl
                 <<" num_node_sets "<<num_node_sets<<std::endl
                 <<" num_side_sets "<<num_side_sets<<std::endl;

    if(verbose)

        std::cout<<"=== End ExodusII Create Info ==="<<std::endl<<std::endl;

    delete [] title;

    return ex_id;

}

void MOERTEL::Mesh::set_vertex_map() {

    num_vertices = 0;

    for (int blk = 0; blk < num_elem_blk; blk++) {

        std::vector<int> temp(num_node_per_elem_in_blk[blk]);
        int num_vertices_in_elem = 3;
        if(4 == num_elem_in_blk[blk] || 9 == num_elem_in_blk[blk])
            num_vertices_in_elem = 4;

        for(int i = 0; i < num_elem_in_blk[blk]; i++) {
            for(int j = 0; j < num_vertices_in_elem; j++) {
                int nodeid =  get_node_id(blk, i, j);
                if( vertex_map.find( nodeid ) == vertex_map.end() ) {
                    vertex_map.insert(std::pair<int, int>(nodeid, num_vertices));
                    //std::cout<<get_node_id(blk, i, j)<<"   "<<num_vertices<<std::endl;
                    num_vertices++;
                }
            }
        }
    }
    return;

}

/*----------------------------------------------------------------------*
 | Report errors to std::cerr											|
 *----------------------------------------------------------------------*/

int MOERTEL::Mesh::ReportError(const std::stringstream &Message) {
	std::cerr << std::endl << Message.str() << std::endl;
	return(-1);
}
