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
  
  int reportSize(const long long &, const long long &, const long long &,std::stringstream & ss);

  virtual int numBlocks(){return inline_bx*inline_by*inline_bz;};
  virtual int blockKstride(){return inline_bx*inline_by;};
  virtual int GlobalNumElements(){return nelx_tot*nely_tot*nelz_tot;};

  virtual void Display_Class(std::ostream&, const std::string &indent); //called by Class_Display

  virtual int Set_Up() = 0;
  int Check_Spans();

  virtual void calculateSize(long long & total_el_count, 
			     long long & total_node_count, 
			     long long & total_edge_count){};

  bool Debug_Location() {return debug_mode;};
  void Debug_Location(bool dm) {debug_mode = dm;};

  std::string getErrorString(){return error_stream.str();}
  std::string getInfoString(){return info_stream.str();}
  std::string getWarningString(){return warning_stream.str();}
 
  int Check_Block_BC_Sets();
 
  virtual int Rename_Block_BC_Sets(){return 0;};

  virtual void Calc_Intervals(){};

  virtual void setStrides();

  virtual void Calc_Serial_Component(Partition * my_part,
			     std::vector <int> & element_vector,
			     std::list <int> & global_node_list,
			     std::vector<int> & global_node_vector,
			     std::map <int, int> & global_node_map,
			     std::map <int, int> & global_element_map);

virtual  void Calc_Parallel_Info(
			 std::vector <int> & element_vector,
			 std::vector<int> & global_node_vector,
			 std::map <int, int> & global_node_map,                             
			 std::list <int> & internal_node_list,
			 std::list <int> & border_nodes_list,
			 std::list <int> & internal_element_list,
			 std::list <int> & border_elements_list,
			 std::list <int> & node_proc_id_list,
			 std::list <int> & element_proc_id_list,
			 std::vector <int> & node_neighbor_vector,
			 std::list <int>  * & boundary_node_list,
			 std::vector <int> & element_neighbor_vector,
			 std::list <std::pair <int ,Topo_Loc > > * & boundary_element_list);

 virtual void getGlobal_Element_Block_Totals(int *);

  virtual int Populate_Sideset_Info(std::map <int, int> & global_element_map,
			     std::map <int, int> & global_node_map,
			     int * const * side_set_elements,
			     int * const * side_set_faces,
			     int * const * side_set_nodes,
			     int * const * side_set_node_counter);

  void Populate_Nodeset_Info(int * const * node_set_nodes,
			     std::map <int, int> & global_node_map);

  virtual void Populate_Map_and_Global_Element_List(int * map, int * gel);
  virtual void Populate_Connectivity(int * const * conn_array, 
			     std::map <int, int> & global_node_map);

  void Populate_Border_Nodes_Elements( int * internal_elements,
				       int * internal_nodes,
				       int * border_elements,
				       int * border_nodes,
				       std::list <int> & internal_node_list,	
				       std::list <int> & border_nodes_list,
				       std::list <int> & internal_element_list,
				       std::list <int> & border_elements_list,
				       std::map <int, int> & global_node_map,
				       std::map <int, int> & global_element_map);

  void Populate_Parallel_Info( int* const * comm_node_ids ,
			       int* const * comm_node_proc_ids,
			       int* const * comm_elem_ids,
			       int* const * comm_side_ids,
			       int* const * comm_elem_proc_ids,
			       std::vector <int> & node_neighbor_vector,
			       std::vector <int> & element_neighbor_vector,
			       std::list <int>  * & boundary_node_list,                   
			       std::map <int, int> & global_node_map,
			       std::list <std::pair <int ,Topo_Loc > > * & boundary_element_list,
			       std::map <int, int> & global_element_map);


  void Populate_Cmap( int * node_cmap_node_cnts,
		      int * node_cmap_ids,
		      int * elem_cmap_elem_cnts,
		      int * elem_cmap_ids,
		      std::vector <int> & node_neighbor_vector,
		      std::vector <int> & element_neighbor_vector,
		      std::list <int>  * & boundary_node_list,                   
		      std::list <std::pair <int ,Topo_Loc > > * & boundary_element_list);

  void Size_BC_Sets(int nnx, 
		    int nny, 
		    int nnz);

  virtual void get_l_i_j_k_from_element_number(int el,
				       int & l,
				       int & i,
				       int & j,
				       int & k);

  virtual void get_l_i_j_k_from_node_number(int nn,
			     int & l,
			     int & i,
			     int & j,
			     int & k);

  virtual int get_element_number_from_l_i_j_k(int l, int i, int j, int k);
  virtual int get_node_number_from_l_i_j_k(int l, int i, int j, int k);

  int get_neighbor(Topo_Loc tl,
		   int ll, 
		   int li, 
		   int lj, 
		   int lk);

  void get_face_nodes(Topo_Loc tl,int global_element_id,int  the_nodes[4]);

  void ZeroSet();

  /* DATA */
  bool debug_mode;
  std::stringstream error_stream;
  std::stringstream info_stream;
  std::stringstream warning_stream;

  InlineGeometryType inline_geometry_type;
  InlineDecompositionType inline_decomposition_type;

  int dimension;
  int trisection_blocks;
  int inline_bx;
  int inline_by;
  int inline_bz;
  int inline_nx;
  int inline_ny;
  int inline_nz;
  int * a_inline_nx;//individual block values
  int * a_inline_ny;
  int * a_inline_nz;
  int * c_inline_nx;//cumulative totals
  int * c_inline_ny;
  int * c_inline_nz;
  int * cum_block_totals;
  int * els_in_block;
  int nelx_tot;
  int nely_tot;
  int nelz_tot;
  int inc_nels[3];
  bool inc_nocuts[3];
  int inline_nprocs[3];
  double ** block_dist;
  double ** c_block_dist;
  double ** first_size;
  double ** last_size;
  int ** interval;
  double inline_gminx;
  double inline_gminy;
  double inline_gminz;
  double inline_gmaxx;
  double inline_gmaxy;
  double inline_gmaxz;

  double transition_radius;

  std::list < PG_BC_Specification * > nodeset_list; 
  std::list < PG_BC_Specification * > sideset_list; 

  PG_BC_Specification * getSideset_by_id(int test_id){
    std::list < PG_BC_Specification * > :: iterator it;
    for(it = sideset_list.begin(); it != sideset_list.end(); it ++){
      if((*it)->id == test_id)return (*it);
    }
    return NULL;
  };

  PG_BC_Specification * getNodeset_by_id(int test_id){
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


  int instride;
  int jnstride;
  int knstride;

  int iestride;
  int jestride;
  int kestride;
// protected:

// private:
  virtual int Calc_Coord_Vectors(){return 0;}
  virtual void Populate_Coords(double * coords,   
		       std::vector<int> & global_node_vector, 
		       std::map <int, int> & global_node_map,
			       int num_nodes){};

  void Customize_Coords(double * coords,int num_nodes,int dim);

  double * Icoors;
  double * Jcoors;
  double * Kcoors;
  int topo_loc_to_exo_face[6];

  LoopLimits getLimits( Topo_Loc the_set_location,
			int sx, int nx, 
			int sy, int ny, 
			int sz, int nz,
			int irange, int jrange);
  
  std::vector <int> * element_block_lists;

  int get_map_entry(std::map < int, int > & the_map, const int & key);


  Partition * base_partition;

  virtual void Build_Global_Lists(std::list <int> & element_list, 
			  std::vector <int> & element_vector,
			  std::list <int> & global_node_list,
			  std::vector <int> & global_node_vector,
			  std::map <int, int> & global_node_map,
			  std::map <int, int> & global_element_map);
  virtual int Element_Proc(int);
  virtual Partition * Decompose(std::list <int> & global_el_ids,int & err_code);
  int get_block_index(int ordinal_val, int count, int * cumulative);
  
  unsigned my_rank;
  unsigned num_processors;

  std::vector < std::pair < int, Topo_Loc > > *  sideset_vectors;
  std::vector < int > *  nodeset_vectors;
  static Inline_Mesh_Desc * static_storage;
  static std::stringstream echo_stream;
  


};

Inline_Mesh_Desc* Parse_Inline_Mesh(std::string & file_name,  
				     std::stringstream & input_stream,
				     int & sparse_error_count,
				     int dim);


}//end of namespace PAMGEN_NEVADA
#endif
