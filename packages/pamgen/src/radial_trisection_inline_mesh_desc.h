// $Id$
#ifndef radial_trisection_inline_mesh_descH
#define radial_trisection_inline_mesh_descH
#include "radial_inline_mesh_desc.h"
#include "Vector.h"
class LoopLimits;
class Partition;
namespace PAMGEN_NEVADA {

class Quad_Patch
{
public:
  Quad_Patch(Vector v0,Vector v1, Vector v2, Vector v3){
    n0 = v0;
    n1 = v1;
    n2 = v2;
    n3 = v3;
  };
  virtual ~Quad_Patch(){};
  Vector n0;
  Vector n1;
  Vector n2;
  Vector n3;
  
  Vector Interpolate(Vector v){
    v *=2.;
    v -= Vector(1,1,1);
    
    Real chi0 = (1.0 - v.x) * (1.0 - v.y);
    Real chi1 = (1.0 + v.x) * (1.0 - v.y);
    Real chi2 = (1.0 + v.x) * (1.0 + v.y);
    Real chi3 = (1.0 - v.x) * (1.0 + v.y);
    
    Vector result  = chi0 * n0;
    result += chi1 * n1;
    result += chi2 * n2;
    result += chi3 * n3;
    
    result /= 4.;
    
    return result;
  }
  void print(){
    std::cout << "Quad " << std::endl;
    std::cout << "n = " << 0 << " " << n0 << std::endl;
    std::cout << "n = " << 1 << " " << n1 << std::endl;
    std::cout << "n = " << 2 << " " << n2 << std::endl;
    std::cout << "n = " << 3 << " " << n3 << std::endl << std::endl;
  }
};


class Radial_Trisection_Inline_Mesh_Desc : public Radial_Inline_Mesh_Desc
{
public:
  
  Radial_Trisection_Inline_Mesh_Desc(int dim){dimension = dim;};
  
  virtual ~Radial_Trisection_Inline_Mesh_Desc(){
    if(tri_block_cum_nn)delete []  tri_block_cum_nn;
  };
  
  virtual int Set_Up();
  virtual void calculateSize(long long & total_el_count, 
			     long long & total_node_count, 
			     long long & total_edge_count);
  virtual void setStrides();
  virtual void Calc_Intervals();
  virtual int Calc_Coord_Vectors();

  virtual int numBlocks(){return((inline_bx-1)*inline_by+1)*inline_bz;}
  virtual int blockKstride(){return(inline_bx-1)*inline_by+1;}

  virtual int GlobalNumElements();


  Vector calc_coords_periodic_trisect_blocks(double total_theta,
					     int nl,
					     int ni, 
					     int nj, 
					     int nk,
					     Quad_Patch ** quads);

  virtual void Populate_Coords(Real * coords,   
		       std::vector<int> & global_node_vector, 
		       std::map <int, int> & global_node_map,
		       int num_nodes);

  LoopLimits get_tri_block_limits(int l, Topo_Loc tl, Topo_Loc & ntl, int nodeset_plus_1);

  int GetBlockBasedGlobalID(int the_el,int bct);

  virtual void getGlobal_Element_Block_Totals(int *);

  virtual int Populate_Sideset_Info(std::map <int, int> & global_element_map,
			     std::map <int, int> & global_node_map,
			     int * const * side_set_elements,
			     int * const * side_set_faces,
			     int * const * side_set_nodes,
			     int * const * side_set_node_counter);

  virtual Partition * Decompose(std::list <int> & global_el_ids,int & err_codes);

  virtual void Build_Global_Lists(std::list <int> & element_list, 
			  std::vector <int> & element_vector,
			  std::list <int> & global_node_list,
			  std::vector <int> & global_node_vector,
			  std::map <int, int> & global_node_map,
			  std::map <int, int> & global_element_map);

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

  virtual void Populate_Connectivity(int * const * conn_array, 
				     std::map <int, int> & global_node_map);

  virtual void Populate_Map_and_Global_Element_List(int * map, int * gel);

  virtual int Rename_Block_BC_Sets();

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
  
  void get_l_and_remainder_from_elno(int el, int & l, int & remainder);
  void get_l_and_remainder_from_node_number(int el, int & l, int & remainder);

  virtual int Element_Proc(int);
    
  int * tri_block_cum_nn;
  int nn_center;
  int div;
  int mod;
};

}//end namespace
#endif
