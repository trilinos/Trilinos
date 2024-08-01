// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

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

      Radial_Trisection_Inline_Mesh_Desc(long long dim):
        tri_block_cum_nn(NULL),
        nn_center(0),
        div(0),
        mod(0)
      {
        dimension = dim;
      };

      virtual ~Radial_Trisection_Inline_Mesh_Desc(){
        if(tri_block_cum_nn)delete []  tri_block_cum_nn;
      };

      virtual long long Set_Up();
      virtual void calculateSize(long long & total_el_count,
          long long & total_node_count,
          long long & total_edge_count);
      virtual void setStrides();
      virtual std::string Calc_Intervals();
      virtual long long Calc_Coord_Vectors();

      virtual long long numBlocks(){return((inline_b[0]-1)*inline_b[1]+1)*inline_b[2];}
      virtual long long blockKstride(){return(inline_b[0]-1)*inline_b[1]+1;}

      virtual long long GlobalNumElements();


      Vector calc_coords_periodic_trisect_blocks(double total_theta,
          long long nl,
          long long ni,
          long long nj,
          long long nk,
          Quad_Patch ** quads);

      virtual void Populate_Coords(Real * coords,
          std::vector<long long> & global_node_vector,
          std::map <long long, long long> & global_node_map,
          long long num_nodes);

      LoopLimits get_tri_block_limits(bool is_block_bc, long long bid, long long l, Topo_Loc tl, Topo_Loc & ntl, long long nodeset_plus_1);

      long long GetBlockBasedGlobalID(long long the_el,long long bct);

      virtual void getGlobal_Element_Block_Totals(long long *);

      virtual long long Populate_Sideset_Info(std::map <long long, long long> & global_element_map,
          std::map <long long, long long> & global_node_map,
          long long * const * side_set_elements,
          long long * const * side_set_faces,
          long long * const * side_set_nodes,
          long long * const * side_set_node_counter);

      virtual long long  Decompose(std::set <long long> & global_el_ids);

      virtual void Build_Global_Lists(const std::set <long long> & element_list,
          std::vector <long long> & element_vector,
          std::list <long long> & global_node_list,
          std::vector <long long> & global_node_vector,
          std::map <long long, long long> & global_node_map,
          std::map <long long, long long> & global_element_map);

      virtual void Calc_Serial_Component(const std::set <long long> & gloabl_element_ids,
          const std::vector<long long> & global_node_vector);

      virtual  void Calc_Parallel_Info(
          const std::vector <long long> & element_vector,
          const std::vector<long long> & global_node_vector,
          const std::map <long long, long long> & global_node_map,
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

      virtual void Populate_Map_and_Global_Element_List(long long * map, long long * gel);

      virtual long long Rename_Block_BC_Sets();

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

      void get_l_and_remainder_from_elno(long long el, long long & l, long long & remainder);
      void get_l_and_remainder_from_node_number(long long el, long long & l, long long & remainder);

      virtual long long getBlockFromElementNumber(long long the_element);

      virtual long long Element_Proc(long long);

      long long * tri_block_cum_nn;
      long long nn_center;
      long long div;
      long long mod;
  };

}//end namespace
#endif
