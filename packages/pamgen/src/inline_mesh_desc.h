// $Id$
#ifndef inline_mesh_descH
#define inline_mesh_descH

#include "inline_geometries.h"
#include "topology_enum.h"
#include "bc_specification.h"
#include "element_density_function.h"
#include "geometry_transform.h"
#include "uns_inline_decomp.h"
#include <sstream>
#include <string>

class Inline_Mesh_Desc;

namespace PAMGEN_NEVADA {





class Inline_Mesh_Desc
{
public:

  Inline_Mesh_Desc(){
    ZeroSet();
  };
  
  virtual ~Inline_Mesh_Desc();
  
  long long reportSize(const long long &, const long long &, const long long &,std::stringstream & ss,long long max_int);

  virtual long long numBlocks(){return inline_bx*inline_by*inline_bz;};
  virtual long long blockKstride(){return inline_bx*inline_by;};
  virtual long long GlobalNumElements(){return nelx_tot*nely_tot*nelz_tot;};

  virtual void Display_Class(std::ostream&, const std::string &indent); //called by Class_Display

  virtual long long Set_Up() = 0;
  long long Check_Spans();

  virtual void calculateSize(long long & total_el_count, 
			     long long & total_node_count, 
			     long long & total_edge_count){};

  bool Debug_Location() {return debug_mode;};
  void Debug_Location(bool dm) {debug_mode = dm;};

  std::string getErrorString(){return error_stream.str();}
  std::string getInfoString(){return info_stream.str();}
  std::string getWarningString(){return warning_stream.str();}
 
  long long Check_Block_BC_Sets();
 
  virtual long long Rename_Block_BC_Sets(){return 0;};

  virtual void Calc_Intervals(){};

  virtual void setStrides();

  virtual void Calc_Serial_Component(Partition * my_part,
			     std::vector <long long> & element_vector,
			     std::list <long long> & global_node_list,
			     std::vector<long long> & global_node_vector,
			     std::map <long long, long long> & global_node_map,
			     std::map <long long, long long> & global_element_map);

virtual  void Calc_Parallel_Info(
			 std::vector <long long> & element_vector,
			 std::vector<long long> & global_node_vector,
			 std::map <long long, long long> & global_node_map,                             
			 std::list <long long> & internal_node_list,
			 std::list <long long> & border_nodes_list,
			 std::list <long long> & internal_element_list,
			 std::list <long long> & border_elements_list,
			 std::list <long long> & node_proc_id_list,
			 std::list <long long> & element_proc_id_list,
			 std::vector <long long> & node_neighbor_vector,
			 std::list <long long>  * & boundary_node_list,
			 std::vector <long long> & element_neighbor_vector,
			 std::list <std::pair <long long ,Topo_Loc > > * & boundary_element_list);

 virtual void getGlobal_Element_Block_Totals(long long *);

  virtual long long Populate_Sideset_Info(std::map <long long, long long> & global_element_map,
			     std::map <long long, long long> & global_node_map,
			     long long * const * side_set_elements,
			     long long * const * side_set_faces,
			     long long * const * side_set_nodes,
			     long long * const * side_set_node_counter);

  void Populate_Nodeset_Info(long long * const * node_set_nodes,
			     std::map <long long, long long> & global_node_map);

  virtual void Populate_Map_and_Global_Element_List(long long * map, long long * gel);
  virtual void Populate_Connectivity(long long * const * conn_array, 
			     std::map <long long, long long> & global_node_map);

  void Populate_Border_Nodes_Elements( long long * internal_elements,
				       long long * internal_nodes,
				       long long * border_elements,
				       long long * border_nodes,
				       std::list <long long> & internal_node_list,	
				       std::list <long long> & border_nodes_list,
				       std::list <long long> & internal_element_list,
				       std::list <long long> & border_elements_list,
				       std::map <long long, long long> & global_node_map,
				       std::map <long long, long long> & global_element_map);

  void Populate_Parallel_Info( long long* const * comm_node_ids ,
			       long long* const * comm_node_proc_ids,
			       long long* const * comm_elem_ids,
			       long long* const * comm_side_ids,
			       long long* const * comm_elem_proc_ids,
			       std::vector <long long> & node_neighbor_vector,
			       std::vector <long long> & element_neighbor_vector,
			       std::list <long long>  * & boundary_node_list,                   
			       std::map <long long, long long> & global_node_map,
			       std::list <std::pair <long long ,Topo_Loc > > * & boundary_element_list,
			       std::map <long long, long long> & global_element_map);


  void Populate_Cmap( long long * node_cmap_node_cnts,
		      long long * node_cmap_ids,
		      long long * elem_cmap_elem_cnts,
		      long long * elem_cmap_ids,
		      std::vector <long long> & node_neighbor_vector,
		      std::vector <long long> & element_neighbor_vector,
		      std::list <long long>  * & boundary_node_list,                   
		      std::list <std::pair <long long ,Topo_Loc > > * & boundary_element_list);

  void Size_BC_Sets(long long nnx, 
		    long long nny, 
		    long long nnz);

  virtual void get_l_i_j_k_from_element_number(long long el,
				       long long & l,
				       long long & i,
				       long long & j,
				       long long & k);

  virtual void get_l_i_j_k_from_node_number(long long nn,
			     long long & l,
			     long long & i,
			     long long & j,
			     long long & k);

  virtual long long get_element_number_from_l_i_j_k(long long l, long long i, long long j, long long k);
  virtual long long get_node_number_from_l_i_j_k(long long l, long long i, long long j, long long k);

  long long get_neighbor(Topo_Loc tl,
		   long long ll, 
		   long long li, 
		   long long lj, 
		   long long lk);

  void get_face_nodes(Topo_Loc tl,long long global_element_id,long long  the_nodes[4]);

  void ZeroSet();

  /* DATA */
  bool debug_mode;
  std::stringstream error_stream;
  std::stringstream info_stream;
  std::stringstream warning_stream;

  InlineGeometryType inline_geometry_type;
  InlineDecompositionType inline_decomposition_type;

  long long dimension;
  long long trisection_blocks;
  long long inline_bx;
  long long inline_by;
  long long inline_bz;
  long long inline_nx;
  long long inline_ny;
  long long inline_nz;
  long long * a_inline_nx;//individual block values
  long long * a_inline_ny;
  long long * a_inline_nz;
  long long * c_inline_nx;//cumulative totals
  long long * c_inline_ny;
  long long * c_inline_nz;
  long long * cum_block_totals;
  long long * els_in_block;
  long long nelx_tot;
  long long nely_tot;
  long long nelz_tot;
  long long inc_nels[3];
  bool inc_nocuts[3];
  long long inline_nprocs[3];
  double ** block_dist;
  double ** c_block_dist;
  double ** first_size;
  double ** last_size;
  long long ** interval;
  double inline_gminx;
  double inline_gminy;
  double inline_gminz;
  double inline_gmaxx;
  double inline_gmaxy;
  double inline_gmaxz;

  double transition_radius;

  std::list < PG_BC_Specification * > nodeset_list; 
  std::list < PG_BC_Specification * > sideset_list; 

  PG_BC_Specification * getSideset_by_id(long long test_id){
    std::list < PG_BC_Specification * > :: iterator it;
    for(it = sideset_list.begin(); it != sideset_list.end(); it ++){
      if((*it)->id == test_id)return (*it);
    }
    return NULL;
  };

  PG_BC_Specification * getNodeset_by_id(long long test_id){
    std::list < PG_BC_Specification * > :: iterator it;
    for(it = nodeset_list.begin(); it != nodeset_list.end(); it ++){
      if((*it)->id == test_id)return (*it);
    }
    return NULL;
  };


  bool periodic_i;
  bool periodic_j;
  bool periodic_k;
  Element_Density_Function * Element_Density_Functions[3];// array of pointers to density functions
  Geometry_Transform *Geometry_Transform_Function;
  bool try_squared;
  bool enforce_periodic;
  /* DATA */


  long long instride;
  long long jnstride;
  long long knstride;

  long long iestride;
  long long jestride;
  long long kestride;
// protected:

// private:
  virtual long long Calc_Coord_Vectors(){return 0;}
  virtual void Populate_Coords(double * coords,   
		       std::vector<long long> & global_node_vector, 
		       std::map <long long, long long> & global_node_map,
			       long long num_nodes){};

  void Customize_Coords(double * coords,long long num_nodes,long long dim);

  double * Icoors;
  double * Jcoors;
  double * Kcoors;
  long long topo_loc_to_exo_face[6];

  LoopLimits getLimits( Topo_Loc the_set_location,
			long long sx, long long nx, 
			long long sy, long long ny, 
			long long sz, long long nz,
			long long irange, long long jrange);
  
  std::vector <long long> * element_block_lists;

  long long get_map_entry(std::map < long long, long long > & the_map, const long long & key);


  Partition * base_partition;

  virtual void Build_Global_Lists(std::list <long long> & element_list, 
			  std::vector <long long> & element_vector,
			  std::list <long long> & global_node_list,
			  std::vector <long long> & global_node_vector,
			  std::map <long long, long long> & global_node_map,
			  std::map <long long, long long> & global_element_map);
  virtual long long Element_Proc(long long);
  virtual Partition * Decompose(std::list <long long> & global_el_ids,long long & err_code);
  long long get_block_index(long long ordinal_val, long long count, long long * cumulative);
  
  long long my_rank;
  long long num_processors;

  std::vector < std::pair < long long, Topo_Loc > > *  sideset_vectors;
  std::vector < long long > *  nodeset_vectors;
  static Inline_Mesh_Desc * static_storage;
  static std::stringstream echo_stream;
  


};

Inline_Mesh_Desc* Parse_Inline_Mesh(std::string & file_name,  
				    std::stringstream & input_stream,
				    long long & sparse_error_count,
				    long long dim,
				    long long max_int);


}//end of namespace PAMGEN_NEVADA
#endif
