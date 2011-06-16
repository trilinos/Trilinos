//@HEADER
/* 
 * 
 * Copyright (2010) Glen Hansen
 * Idaho National Laboratory
 *
 * This library is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of the
 * License, or (at your option) any later version.
 *  
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *                                                                           
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 * USA                                                                       
 * 
 * Questions? Contact Glen Hansen (Glen.Hansen@inl.gov), Idaho National Laboratory
 * 
 */
//@HEADER

#ifndef MORTEL_MESH_H
#define MORTEL_MESH_H

#include "exodusII.h"

#include <vector>
#include <string>
#include <map>

namespace MOERTEL {

class Mesh {
public:

#ifdef HAVE_MOERTEL_NEMESIS  // Use NEMESIS for parallel IO
    Mesh(const int proc_id, const int nprocs, const bool v = false);
#else
    Mesh(const int proc_id = 0, const bool v = false);
#endif

    ~Mesh();

// We need to read and write Exodus files

    int read_exodus(const char * filename);

// This class supports two types of file output. The first is simple, where the mesh data and results are
// stored in the class and written at the very end of execution, using the below.

    int write_exodus(const char * filename);

// The second method is for larger files and transient simulations, where the file is written as the simulation
// proceeds. In this approach, the following opens the output file and writes the mesh to it.

    int open_exodus_write_mesh(const char * filename);
    void set_time_step(int its, double itime);
    int update_file();

// We need a set of convenient functions to retrieve data from this object, and write data to it

    void compute_nodal_adj();

    int add_nodal_data(std::string name, std::vector<double> &data);
    int add_nodal_data(std::string name, double *data);

    int add_elem_data(std::string name, std::vector<double> &data);
    int add_elem_data(std::string name, double *data);

    void set_verbose(const bool v = true);

#ifdef HAVE_MOERTEL_NEMESIS
    int get_num_global_nodes() {
        return ne_num_global_nodes;
    }
#else
    int get_num_global_nodes() {
        return num_nodes;
    }
#endif
    int get_num_nodes() {
        return num_nodes;
    }
    int get_num_my_nodes() {
        return my_node_num_map.size();
    }
    int get_num_dim() {
        return num_dim;
    }
    int get_num_elem() {
        return num_elem;
    }
    int get_num_elem_blks() {
        return num_elem_blk;
    }
    int get_num_node_sets() {
        return num_node_sets;
    }
    int get_num_side_sets() {
        return num_side_sets;
    }

    int get_num_elem_in_blk(int blk) {
        return num_elem_in_blk[blk];
    }
    int get_num_nodes_per_elem_in_blk(int blk) {
        return num_node_per_elem_in_blk[blk];
    }
	std::string get_elem_type_in_blk(int blk) {
        return blk_elem_type[blk];
    }
    int *get_elem_nodes(int blk, int elem) {
        return &connect[blk][elem * num_node_per_elem_in_blk[blk]];
    }
    int& get_node_id(int blk, int elem, int offset) {
        return connect[blk][elem * num_node_per_elem_in_blk[blk] + offset];
    }
    int get_boundary_status(int blk, int elem);
    int get_node_boundary_status(int nodeid);

    std::vector<int> *get_node_num_map() {
        return &node_num_map;
    }
    std::vector<int> *get_elem_num_map() {
        return &elem_num_map;
    }
    std::vector<int> *get_my_node_num_map() {
        return &my_node_num_map;
    }

    double get_x(int i) {
        return x[i];
    }
    double get_y(int i) {
        return y[i];
    }
    double get_z(int i) {
        return z[i];
    }

    std::vector<double> *get_x_vector() {
        return &x;
    }
    std::vector<double> *get_y_vector() {
        return &y;
    }
    std::vector<double> *get_z_vector() {
        return &z;
    }

    std::vector<int>& get_nodal_adj(int i) {
        return nodal_adj[i];
    }
    std::vector<int>& get_node_set(int i) {
        return ns_node_list[i];
    }
    std::vector<int>& get_side_set(int i) {
        return ss_node_list[i];
    }
    int get_node_set_entry(int i, int j) {
        return ns_node_list[i][j];
    }
	int get_num_sides_in_set(int i){
			return num_sides_per_ss[i];
	}
	int get_num_nodes_on_side(int i, int j){
			return ss_ctr_list[i][j];
	}
    int get_side_set_node_entry(int i, int j) {
        return ss_node_list[i][j];
    }
    int get_node_set_value(int i) {
        return node_set_map[i];
    }
    int get_side_set_node_value(int i) {
        return side_set_node_map[i];
    }
    std::map<int,int>& get_vertex_map() {
        return vertex_map;
    }
    void set_vertex_map();
    int get_num_vertices() {
        return num_vertices;
    }

private:

    bool verbose;
    std::string title;
    float exodus_version;

    int num_dim;
    int num_nodes;
    int num_elem;
    int num_elem_blk;
    int num_node_sets;
    int num_side_sets;
    int num_vertices;


    std::vector<double> x;               // x locations of node points
    std::vector<double> y;               // y locations of node points
    std::vector<double> z;               // z locations of node points

    std::vector<int> blk_ids;
    std::vector<int> num_elem_in_blk;
    std::vector<std::string> blk_elem_type;
    std::vector<int> num_node_per_elem_in_blk;
    std::vector<std::vector<int> > connect;

    std::vector<int> ss_ids;
    std::vector<int> num_sides_per_ss;
    std::vector<int> num_df_per_ss;
    std::vector<std::vector<int> > ss_elem_list;
    std::vector<std::vector<int> > ss_side_list;
    std::vector<std::vector<int> > ss_node_list;
    std::vector<std::vector<int> > ss_ctr_list;

    std::vector<int> ns_ids;
    std::vector<int> num_nodes_per_ns;
    std::vector<int> num_df_per_ns;
    std::vector<std::vector<int> > ns_node_list;
    std::vector<std::vector<int> > ns_ctr_list;

    std::vector<std::vector<int> > nodal_adj; //cn we may only need this for epetra
    std::vector<int> nodal_adj_idx;
    std::vector<int> nodal_adj_array;

    std::vector<std::string> nodal_field_names;
    std::vector<std::vector<double> > nodal_fields;

    std::vector<std::string> elem_field_names;
    std::vector<std::vector<double> > elem_fields;

    std::map<int,int> vertex_map;

    int create_exodus(const char * filename);
    int write_nodal_coordinates_exodus(int ex_id);
    int write_element_blocks_exodus(int ex_id);
    int write_nodal_data_exodus(int ex_id);
    int write_elem_data_exodus(int ex_id);
    int close_exodus(int ex_id);
    void check_exodus_error(const int ex_err,const std::string msg);

    std::vector<int> node_num_map;
    std::vector<int> elem_num_map;

    std::vector<int> node_set_map;
    std::vector<int> side_set_node_map;
    std::vector<int> my_node_num_map;

	// Error functions
	
	int ReportError(const std::stringstream &Message);

    // Variables for staged output of exodus files

    int staged_ex_id;
    int ts;
    double time;
    int nodal_names_written, elem_names_written;


#ifdef HAVE_MOERTEL_NEMESIS
    int ne_num_global_nodes, ne_num_global_elems, ne_num_global_elem_blks,
        ne_num_global_node_sets, ne_num_global_side_sets;

    int num_internal_nodes, num_border_nodes, num_external_nodes,
        num_internal_elems, num_border_elems, num_node_cmaps, num_elem_cmaps;

    std::vector<int> node_cmap_ids, node_cmap_node_cnts, elem_cmap_ids, elem_cmap_elem_cnts;

    std::vector<int> elem_mapi, elem_mapb, node_mapi, node_mapb, node_mape;

    std::vector<int> global_ns_ids, num_global_node_counts, num_global_node_df_counts;
    std::vector<int> global_ss_ids, num_global_side_counts, num_global_side_df_counts;
    std::vector<int> global_elem_blk_ids, global_elem_blk_cnts;

    std::vector<std::vector<int> > node_ids_in_cmap, n_proc_ids_in_cmap;
    std::vector<std::vector<int> > elem_ids_in_cmap, e_side_ids_in_cmap, e_proc_ids_in_cmap;

    int proc_id, nprocs, nprocs_infile;

    char filetype;

#endif

};

}
#endif
