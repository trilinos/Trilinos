// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef inline_mesh_descH
#define inline_mesh_descH

#ifdef _MSC_VER
# ifndef M_PI
# define M_PI        3.14159265358979323846
# endif
#endif

#include "inline_geometries.h"
#include "topology_enum.h"
#include "pamgen_bc_specification.h"
#include "element_density_function.h"
#include "geometry_transform.h"
#include "uns_inline_decomp.h"
#include <sstream>
#include <string>
#include <set>

class Inline_Mesh_Desc;

namespace PAMGEN_NEVADA {

  class Inline_Mesh_Desc
  {
    public:

      Inline_Mesh_Desc():
        dimension(0),
        instride(0),
        jnstride(0),
        knstride(0),
        iestride(0),
        jestride(0),
        kestride(0),
        element_block_lists(NULL)
      {
        ZeroSet();
      };

      virtual ~Inline_Mesh_Desc();

      long long reportSize(
          const long long &,
          const long long &,
          const long long &,
          std::stringstream & ss,
          long long max_int
          );

      virtual long long numBlocks()
      {
        return inline_b[0]*inline_b[1]*inline_b[2];
      };

      virtual long long blockKstride()
      {
        return inline_b[0]*inline_b[1];
      };

      virtual long long GlobalNumElements()
      {
        return nel_tot[0]*nel_tot[1]*nel_tot[2];
      };

      virtual void Display_Class(std::ostream&, const std::string &indent); //called by Class_Display

      virtual long long Set_Up() = 0;
      long long Check_Spans();

      virtual void calculateSize(
          long long & /* total_el_count */,
          long long & /* total_node_count */,
          long long & /* total_edge_count */
          )
      {
      };

      bool Debug_Location() {return debug_mode;};
      void Debug_Location(bool dm) {debug_mode = dm;};

      std::string getErrorString(){
        std::stringstream tstring;
        tstring << error_stream.str();
        Inline_Mesh_Desc * imd = next;
        while (imd){
          if(imd->error_stream.str().size()>0){
            tstring << "\n";
            tstring << imd->error_stream.str();
          }
          imd = imd->next;
        }
        return tstring.str();
      }

      std::string getInfoString(){
        std::stringstream tstring;
        tstring << info_stream.str();
        Inline_Mesh_Desc * imd = next;
        while (imd){
          if(imd->info_stream.str().size()>0){
            tstring << "\n";
            tstring << imd->info_stream.str();
          }
          imd = imd->next;
        }
        return tstring.str();
      }

      std::string getWarningString(){
        std::stringstream tstring;
        tstring << warning_stream.str();
        Inline_Mesh_Desc * imd = next;
        while (imd){
          if(imd->warning_stream.str().size()>0){
            tstring << "\n";
            tstring << imd->warning_stream.str();
          }
          imd = imd->next;
        }
        return tstring.str();
      }

      long long Check_Block_BC_Sets();
      long long Check_Blocks();

      void brokerDecompositionStrategy(){
        /* if any blocks are suppressed set the decomposition to sequential */
        if(!suppressed_blocks.empty()) inline_decomposition_type = SEQUENTIAL;
      };

      virtual long long Rename_Block_BC_Sets(){return 0;};

      virtual std::string Calc_Intervals(){
        std::string errorString;
        return errorString;
      };

      virtual void setStrides();

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
      long long inline_b[3];
      long long inline_n[3];
      long long * a_inline_n[3];//individual block values
      long long * c_inline_n[3];//cumulative totals
      long long * cum_block_totals;
      long long * els_in_block;
      long long nel_tot[3];
      long long inc_nels[3];
      bool inc_nocuts[3];
      long long inline_nprocs[3];
      double ** block_dist;
      double ** c_block_dist;
      double ** first_size;
      double ** last_size;
      long long ** interval;
      long long inline_block_start;
      double inline_offset[3];
      double inline_gmin[3];
      double inline_gmax[3];

      double transition_radius;

      long long total_unsupressed_elements;
      std::vector < long long > sequential_decomp_limits;

      std::list < PG_BC_Specification * > nodeset_list;
      std::list < PG_BC_Specification * > sideset_list;

      std::set < long long > suppressed_blocks;

      void addSuppressedBlock(int ablock){
        suppressed_blocks.insert(ablock);
      }

      bool isBlockSuppressed(int ablock){
        if(suppressed_blocks.find(ablock) != suppressed_blocks.end())return true;
        return false;
      }

      bool isElementSuppressed(long long elno){
        long long bid = getBlockFromElementNumber(elno);
        bool is_sup = isBlockSuppressed(bid+1);
        return is_sup;
      }

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
      virtual long long Calc_Coord_Vectors();
      virtual void Populate_Coords(
          double * /* coords */,
          std::vector<long long> & /* global_node_vector */,
          std::map <long long, long long> & /* global_node_map */,
          long long /* num_nodes */
          )
      {};

      void Offset_Coords(double * coords,long long num_nodes,long long dim);

      void Customize_Coords(double * coords,long long num_nodes,long long dim);

      double * IJKcoors[3];

      long long topo_loc_to_exo_face[6];

      LoopLimits getLimits(
          Topo_Loc the_set_location,
          long long sx, long long nx,
          long long sy, long long ny,
          long long sz, long long nz,
          long long irange, long long jrange
          );

      std::vector <long long> * element_block_lists;

      long long get_map_entry(const std::map < long long, long long > & the_map, const long long & key);


      Partition * base_partition;

      virtual void Build_Global_Lists(
          const std::set <long long> & element_list,
          std::vector <long long> & element_vector,
          std::list <long long> & global_node_list,
          std::vector <long long> & global_node_vector,
          std::map <long long, long long> & global_node_map,
          std::map <long long, long long> & global_element_map
          );
      virtual long long Element_Proc(long long);
      virtual long long Decompose(std::set <long long> & global_el_ids);
      virtual long long Decompose(std::set <long long> & global_el_ids, long long* lNPD,
                                  long long* gNPD);
      long long DecomposeSequential(std::set <long long> & global_el_ids);
      long long get_block_index(long long ordinal_val, long long count, long long * cumulative);

      virtual long long getBlockFromElementNumber(long long the_element);

      unsigned my_rank;
      unsigned num_processors;

      std::vector < std::pair < long long, Topo_Loc > > *  sideset_vectors;
      long long  *  sideset_global_count;
      long long  *  nodeset_global_count;
      std::vector < long long > *  nodeset_vectors;
      static Inline_Mesh_Desc * im_static_storage;
      static Inline_Mesh_Desc * first_im_static_storage;

      static void addDisc(Inline_Mesh_Desc * imd){
        /*set the first pointer if unset*/
        /*add the new entry to next if there is an im_static_storage*/
        /*set im_static_storage*/
        if(!first_im_static_storage)first_im_static_storage = imd;
        if(im_static_storage)im_static_storage->next = imd;
        im_static_storage = imd;
      }

      Inline_Mesh_Desc * next;
      //   static std::vector <Inline_Mesh_Desc *> InlineMeshDescVector;
      //   static int curr_inline_mesh;
      static std::stringstream echo_stream;

    private:
	std::vector <Partition *> sorted_partition_list;
    public:
	const std::vector<Partition *> & get_sorted_partition_list() const {return sorted_partition_list;}

  };

  Inline_Mesh_Desc* Parse_Inline_Mesh(std::string & file_name,
      std::stringstream & input_stream,
      long long & sparse_error_count,
      long long dim,
      long long max_int);


}//end of namespace PAMGEN_NEVADA
#endif
