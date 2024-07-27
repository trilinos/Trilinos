// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "calc_decomp_cuts.h"
#include <math.h>
#include <stdlib.h>
#include "uns_inline_decomp.h"
#include "inline_mesh_desc.h"
#include <sstream>
#include <set>
#include "Random.h"
#include <iostream>
#include <algorithm>
#include <limits.h>
//#pragma warning(disable:981)
//#pragma warning(disable:383)
namespace PAMGEN_NEVADA {

  Inline_Mesh_Desc * Inline_Mesh_Desc::im_static_storage = NULL;
  Inline_Mesh_Desc * Inline_Mesh_Desc::first_im_static_storage = NULL;
  std::stringstream Inline_Mesh_Desc::echo_stream;



bool part_compare_size(const Partition *a, const Partition *b) {
  if(a->numels < b->numels)
    return true;
  if(a->numels == b->numels)
    return (a->unique_id < b->unique_id);
  return false;
}

bool part_compare_centroid(const Partition *a, const Partition *b) {
  if(a->centroid < b->centroid)
    return true;
  if(a->centroid == b->centroid)
    return (a->unique_id < b->unique_id);
  return false;
}




/*****************************************************************************/
Inline_Mesh_Desc::~Inline_Mesh_Desc()
/*****************************************************************************/
{
  if(next)delete next;
  next = NULL;
  for(long long i = 0; i < 3; i ++){
    if(block_dist[i])delete []  block_dist[i];
    if(c_block_dist[i])delete []  c_block_dist[i];
    if(first_size[i])delete  [] first_size[i];
    if(last_size[i])delete  [] last_size[i];
    if(interval[i])delete  [] interval[i];
  }
  delete [] block_dist;
  delete [] c_block_dist;
  delete [] first_size;
  delete [] last_size;
  delete [] interval;

  delete [] a_inline_n[0];
  delete [] a_inline_n[1];
  delete [] a_inline_n[2];

  delete [] c_inline_n[0];
  delete [] c_inline_n[1];
  delete [] c_inline_n[2];

  if(cum_block_totals)delete [] cum_block_totals;
  if(els_in_block) delete [] els_in_block;

  if(Element_Density_Functions[0])delete Element_Density_Functions[0];
  if(Element_Density_Functions[1])delete Element_Density_Functions[1];
  if(Element_Density_Functions[2])delete Element_Density_Functions[2];
  if(Geometry_Transform_Function)delete Geometry_Transform_Function;
  if(base_partition){
    base_partition->empty();
    delete base_partition;
    base_partition = NULL;
  }
  if (IJKcoors[0])delete [] IJKcoors[0];
  if (IJKcoors[1])delete [] IJKcoors[1];
  if (IJKcoors[2])delete [] IJKcoors[2];
  if(element_block_lists)delete [] element_block_lists;
  if(sideset_list.size())delete [] sideset_vectors;
  if(sideset_list.size())delete [] sideset_global_count;
  if(nodeset_list.size())delete [] nodeset_vectors;

  std::list < PG_BC_Specification * > :: iterator it;
  for(it = nodeset_list.begin(); it != nodeset_list.end(); it ++){
    PG_BC_Specification * bcs = *it;
    delete bcs;
  }
  for(it = sideset_list.begin(); it != sideset_list.end(); it ++){
    PG_BC_Specification * bcs = *it;
    delete bcs;
  }
}

/****************************************************************************/
long long Inline_Mesh_Desc::reportSize(const long long & total_el_count, 
				       const long long & total_node_count, 
				       const long long & total_edge_count,
				       std::stringstream & es,
				       long long max_int)
/****************************************************************************/
{
  long long status = 0;
  if(total_el_count > max_int){
    es << "Terminating from Inline_Mesh_Desc, ";
    es << total_el_count ;
    es << " elements requested, max ";
    es << max_int;
    es << " elements permitted.";
    status = 1;
  }

  if(total_node_count > max_int){
    es << "Terminating from Inline_Mesh_Desc,";
    es << total_node_count ;
    es << " nodes requested, max ";
    es << max_int;
    es << " nodes permitted.";
    status = 1;
  }

  if(total_edge_count > max_int){
    es << "Terminating from Inline_Mesh_Desc,";
    es << total_edge_count ;
    es << " edges requested, max ";
    es << max_int;
    es << " edges permitted.";
    status = 1;
  }
  return status;
}

//! This is an access function to protect agains the stl maps
//! behaviour of always returning an entry for every [] operator
//! query.
/****************************************************************************/
long long Inline_Mesh_Desc::get_map_entry(const std::map < long long, long long > & the_map, const long long & key)
/****************************************************************************/
{
  std::map <long long, long long > ::const_iterator foo;
  foo = the_map.find(key);
  if(foo == the_map.end()){
    error_stream << "Looking for but not finding key entry " << key << "\n";
  }
  return (*foo).second;
}

/****************************************************************************/
long long Inline_Mesh_Desc::get_block_index(long long ordinal_val, 
				   long long count ,
				   long long * cumulative)//c_inline_n[2]);
/****************************************************************************/

{
  long long i = 0;
  while(i < count-1){
    if(ordinal_val >= cumulative[i] && ordinal_val < cumulative[i+1])return i;
    i ++;
  }
  return i;
}


//! Queries which processor an element lies on.
//! Calls the recursive Partition::Element_Proc function.
/****************************************************************************/
long long Inline_Mesh_Desc::Element_Proc(long long global_element_id)
/****************************************************************************/
{
  long long proc = 0;
  if(inline_decomposition_type == SEQUENTIAL){
    for(unsigned ict = 0; ict < sequential_decomp_limits.size();ict ++){
      if(sequential_decomp_limits[ict] >= global_element_id)return ict;
    }
  }
  else if(inline_decomposition_type == RANDOM){
    SRANDOM(global_element_id);
    long long rand_num = RANDOM();
    proc = rand_num%num_processors;
  }
  else if((inline_decomposition_type == BISECTION) || (inline_decomposition_type == PROCESSOR_LAYOUT)){
    long long l,i,j,k;
    get_l_i_j_k_from_element_number(global_element_id,l,i,j,k);
    long long ginds[4];
    ginds[0] = i;
    ginds[1] = j;
    ginds[2] = k;
    ginds[3] = l;
    return base_partition->Element_Proc(ginds);
  }

  return proc;
}

/****************************************************************************/
long long Inline_Mesh_Desc::DecomposeSequential(std::set <long long> & global_el_ids)
/****************************************************************************/
{
  sequential_decomp_limits.resize(num_processors);
  long long ltotal =  kestride * nel_tot[2];
  long long total = total_unsupressed_elements;
  long long num_per_proc = total/num_processors;
  long long remainder = total - num_per_proc*num_processors;
  long long my_start = my_rank * num_per_proc;
  long long my_end = my_start + num_per_proc;
  if(my_rank == num_processors-1)my_end +=remainder;
  long long mtotal = 0;
  for(long long acount = 0; acount < ltotal; acount ++){
    if(!isElementSuppressed(acount)){
      if(mtotal >= my_start && mtotal < my_end){
	global_el_ids.insert(acount);
      }
      long long the_proc = mtotal/num_per_proc;
      if(the_proc > (num_processors -1))the_proc = num_processors - 1;
      sequential_decomp_limits[the_proc] = acount;
      mtotal ++;
    }
  }
  sequential_decomp_limits[num_processors-1] = ltotal;
  return 0;
}



long long Inline_Mesh_Desc::Decompose(std::set <long long> & global_el_ids){
  long long local_ijk[3];
  long long global_ijk[6];
  long long return_value;
  return_value = Inline_Mesh_Desc::Decompose(global_el_ids,local_ijk,global_ijk);
  return return_value;
}

//! Partitions all elements by a recursive bisection into 
//! rectangular chunks. Alternative decompositions are possible.
//! This one is simple to implement and easy to run on all processors.
//! Zoltan could be used for this. The elements on the processor are
//! packed into the stl list.
/****************************************************************************/
  long long Inline_Mesh_Desc::Decompose(std::set <long long> & global_el_ids, long long* local_ijk,
                                        long long* global_ijk)
/****************************************************************************/
{
  //Recursive Bisection decomposition
  //Create a partition object having the entire domain
  //Place it in a list sorted by number of elements


  //While sorted_list.size < num_procs
  //  get first entry in list and the partition object it points to
  //  bisect the partition object
  //  remove the entry from the sorted_list
  //  add the children of the partition object to the sorted list
  //  re-sort the list

  //Then put the partition objects into a list sorted by centroid distance from bottom left corner
  //Assign proc id numbers to the centroid sorted list entries

  //Then loop over the elements for the current processor and push them into
  //the global el_ids list

  /**********************************************/
  //Recursive Bisection decomposition
  //Create a partition object having the entire domain

  if( Debug_Location())
    std::cout << "Inline_Mesh_Desc::Decompose()" << std::endl;
  

if(inline_decomposition_type == PROCESSOR_LAYOUT){

  long long remaining_cuts[3];
  remaining_cuts[0] = inline_nprocs[0];
  remaining_cuts[1] = inline_nprocs[1];
  remaining_cuts[2] = inline_nprocs[2];

  base_partition  = new Partition(0,0,0,0,1,nel_tot[0],nel_tot[1],nel_tot[2],inline_decomposition_type,remaining_cuts);
  //Place it in a list sorted by number of elements
  sorted_partition_list.push_back(base_partition);
  Partition* biggest;

  if(num_processors != (long long)inline_nprocs[0]*inline_nprocs[1]*inline_nprocs[2]){
    error_stream << "Inline_Mesh_Desc::Decompose "
		 << "The specified inline processor layout " 
		 << inline_nprocs[0] << " X " 
		 << inline_nprocs[1] << " X " 
		 << inline_nprocs[2] << " does not correspond to the number of processors " << num_processors << "." ;
    return 1;
   }
    inc_nels[0] = nel_tot[0]/inline_nprocs[0];
    inc_nels[1] = nel_tot[1]/inline_nprocs[1];
    inc_nels[2] = nel_tot[2]/inline_nprocs[2];


    if(inc_nels[0] == 0 ){
      error_stream << "Inline_Mesh_Desc::Decompose"
		   << " Value for numprocs specified in I direction " << inline_nprocs[0] 
		   << " is greater than the number of elements in I " << nel_tot[0] << ".";
      return 1;
    }

    if(inc_nels[1] == 0 ){
      error_stream << "Inline_Mesh_Desc::Decompose"
		   << " Value for numprocs specified in J direction " << inline_nprocs[1] 
		   << " is greater than the number of elements in J " << nel_tot[1] << ".";
      return 1;
    }

    if(inc_nels[2] == 0 ){
      error_stream << "Inline_Mesh_Desc::Decompose"
		   << " Value for numprocs specified in K direction " << inline_nprocs[2] 
		   << " is greater than the number of elements in K " << nel_tot[2] << ".";
      return 1;
    }
    
    info_stream << "Using PROCESSOR LAYOUT decomposition.\n";
    info_stream << "Number of elements/segment in directions I/X/R \t\t" << inc_nels[0] << "\n";
    info_stream << "Number of elements/segment in directions J/Y/THETA \t" << inc_nels[1] << "\n";
    info_stream << "Number of elements/segment in directions K/Z/PHI \t" << inc_nels[2] << "\n";
    info_stream << "Number of mesh segments in directions I/X/R \t\t" << remaining_cuts[0] << "\n";
    info_stream << "Number of mesh segments in directions J/Y/THETA \t" << remaining_cuts[1] << "\n";
    info_stream << "Number of mesh segments in directions K/Z/PHI \t" << remaining_cuts[2] << "\n";

    while(sorted_partition_list.size() < num_processors){
      //  get first entry in list and the partition object it points to
      biggest = sorted_partition_list.back();
      //  remove the entry from the sorted_list
      sorted_partition_list.pop_back();
      
      //  bisect the partition object and
      //  add the children of the partition object to the sorted list
      biggest->Processor_Partition(sorted_partition_list,inc_nels);
      //  re-sort the list
      std::sort(sorted_partition_list.begin(),sorted_partition_list.end(),part_compare_size);
    }
  }
  else if (inline_decomposition_type == BISECTION){

    //SETUP for bisection
  long long remaining_cuts[3];
  long long decomp_result = 0;
  if(dimension == 3){
    decomp_result = dom_decomp_3d(nel_tot[0],nel_tot[1],nel_tot[2],num_processors,&(inc_nels[0]),&(inc_nels[1]),&(inc_nels[2]));
  }
  else{
    decomp_result = dom_decomp_2d(nel_tot[0],nel_tot[1],num_processors,&(inc_nels[0]),&(inc_nels[1]));
  }

    if(decomp_result != 0){
      error_stream << "Terminating from Inline_Mesh_Desc::Decompose, ";
      error_stream << "non-zero return value from dom_decomp_2/3d ";
      error_stream << decomp_result;
      return 1;
    }

  if(dimension == 3){
    remaining_cuts[0] = inc_nels[0];
    remaining_cuts[1] = inc_nels[1];
    remaining_cuts[2] = inc_nels[2];
    inc_nels[0] = nel_tot[0]/inc_nels[0];
    inc_nels[1] = nel_tot[1]/inc_nels[1];
    inc_nels[2] = nel_tot[2]/inc_nels[2];
  }
  else{
    remaining_cuts[0] = inc_nels[0];
    remaining_cuts[1] = inc_nels[1];
    remaining_cuts[2] = 1;
    inc_nels[0] = nel_tot[0]/inc_nels[0];
    inc_nels[1] = nel_tot[1]/inc_nels[1];
    inc_nels[2] = 1;
  }

    {
      info_stream << "Using BISECTION LAYOUT decomposition.\n";
      info_stream << "Number of elements/segment in directions I/X/R \t\t" << inc_nels[0] << "\n";
      info_stream << "Number of elements/segment in directions J/Y/THETA \t" << inc_nels[1] << "\n";
      info_stream << "Number of elements/segment in directions K/Z/PHI \t" << inc_nels[2] << "\n";
      info_stream << "Number of mesh segments in directions I/X/R \t\t" << remaining_cuts[0] << "\n";
      info_stream << "Number of mesh segments in directions J/Y/THETA \t" << remaining_cuts[1] << "\n";
      info_stream << "Number of mesh segments in directions K/Z/PHI \t" << remaining_cuts[2] << "\n";
    }

    base_partition  = new Partition(0,0,0,0,1,nel_tot[0],nel_tot[1],nel_tot[2],inline_decomposition_type,remaining_cuts);
  //Place it in a list sorted by number of elements
  sorted_partition_list.push_back(base_partition);
  Partition* biggest;

    while(sorted_partition_list.size() < num_processors){
      //  get first entry in list and the partition object it points to
      biggest = sorted_partition_list.back();
      //  remove the entry from the sorted_list
      sorted_partition_list.pop_back();

      {
        biggest->Processor_Partition(sorted_partition_list, inc_nels);
      }
      //  re-sort the list
      std::sort(sorted_partition_list.begin(),sorted_partition_list.end(),part_compare_size);//sorted_partition_list.sort();
    }
  }
  else if(inline_decomposition_type == RANDOM){
    long long total = kestride * nel_tot[2];
    for(long long i = 0; i < total; i ++){
      SRANDOM(i);
      long long rand_num = RANDOM();
      unsigned proc = rand_num%num_processors;
      if(proc == my_rank)global_el_ids.insert(i);
    }
    return 0;
  }
  else if(inline_decomposition_type == SEQUENTIAL){
    return DecomposeSequential(global_el_ids);
  }  

  std::sort(sorted_partition_list.begin(),sorted_partition_list.end(),part_compare_centroid);//sorted_partition_list.sort();


  //Assign proc id numbers to the centroid sorted list entries
  std::vector < Partition * > :: iterator citer;
  long long proc_cnt = 0;
  Partition *my_part = NULL;
  for(citer = sorted_partition_list.begin();citer != sorted_partition_list.end();citer ++,proc_cnt++){
    (*citer)->proc_id = proc_cnt;
    if(proc_cnt == my_rank)my_part = (*citer);
  }


  // Get local ijk
  local_ijk[0] = my_part->highs[0] - my_part->lows[0] + 1;
  local_ijk[1] = my_part->highs[1] - my_part->lows[1] + 1;
  local_ijk[2] = my_part->highs[2] - my_part->lows[2] + 1;
  global_ijk[0] = my_part->lows[0];
  global_ijk[1] = my_part->highs[0];
  global_ijk[2] = my_part->lows[1];
  global_ijk[3] = my_part->highs[1];
  global_ijk[4] = my_part->lows[2];
  global_ijk[5] = my_part->highs[2];
  //
  
  //Then loop over the elements for the current processor and push them into
  //the global el_ids list
  for(long long k = my_part->lows[2]; k < my_part->highs[2]; k++){
    for(long long j = my_part->lows[1]; j < my_part->highs[1]; j++){
      for(long long i = my_part->lows[0]; i < my_part->highs[0]; i++){
        long long elnum = iestride*i+jestride*j+kestride*k;
        global_el_ids.insert(elnum);
      }
    }
  }

  if( Debug_Location())
    std::cout << "Inline_Mesh_Desc::Leaving-Decompose()" << std::endl;

  return 0;
}

/****************************************************************************/
long long Inline_Mesh_Desc::getBlockFromElementNumber(long long the_element)
/****************************************************************************/
{
  long long global_k = the_element/(kestride);
  long long global_j = (the_element - global_k*(kestride))/(jestride);
  long long global_i = the_element - global_k*(kestride)-global_j*(jestride);
  
  // these are the indices of the block in which the element resides
  //       long long block_k = global_k/(inline_nz);
  //       long long block_j = global_j/(inline_ny);
  //       long long block_i = global_i/(inline_nx);
  long long block_k = get_block_index(global_k,inline_b[2],c_inline_n[2]);
  long long block_j = get_block_index(global_j,inline_b[1],c_inline_n[1]);
  long long block_i = get_block_index(global_i,inline_b[0],c_inline_n[0]);
  
  // This is the ordinal number of the block the element resides in
  long long local_block = block_i + block_j*(inline_b[0])+ block_k*(blockKstride());
  return local_block;
}

//! A utility function to build up required bookkeeping objects.
/****************************************************************************/
void Inline_Mesh_Desc::Build_Global_Lists(const std::set <long long> & global_element_ids,
                                       std::vector <long long> & element_vector,
                                       std::list <long long> & global_node_list,
                                       std::vector <long long> & global_node_vector,
                                       std::map <long long, long long> & global_node_map,
                                       std::map <long long, long long> & global_element_map)
/****************************************************************************/
{
  if( Debug_Location()) 
    std::cout << "Inline_Mesh_Desc::Build_Global_Lists()" << std::endl;

  for(std::set <long long>::iterator the_it = global_element_ids.begin();the_it != global_element_ids.end();the_it++){
    element_vector.push_back(*the_it);
  }
  element_block_lists = new std::vector <long long> [numBlocks()];

  std::set <long long> ::iterator lit;
  for(lit = global_element_ids.begin();lit != global_element_ids.end();lit ++){
    long long the_element = *lit;

    long long global_k = the_element/(kestride);
    long long global_j = (the_element - global_k*(kestride))/(jestride);
    long long global_i = the_element - global_k*(kestride)-global_j*(jestride);
    
    // This is the ordinal number of the block the element resides in
    long long local_block = getBlockFromElementNumber(the_element);
    element_block_lists[local_block].push_back(the_element);
    long long nn;
    long long block_j = get_block_index(global_j,inline_b[1],c_inline_n[1]);
    if(periodic_j && (block_j == (inline_b[1]-1)) && (global_j == (nel_tot[1]-1))){
      if(dimension == 2){
      nn = (global_i+0)*instride + (global_j+0)*jnstride;                         global_node_list.push_back(nn);
      nn = (global_i+1)*instride + (global_j+0)*jnstride;                         global_node_list.push_back(nn);
      nn = (global_i+1)*instride + (0+0)*jnstride+(global_k+0)*knstride;      global_node_list.push_back(nn);
      nn = (global_i+0)*instride + (0+0)*jnstride+(global_k+0)*knstride;      global_node_list.push_back(nn);
      }
      else{
      nn = (global_i+0)*instride + (global_j+0)*jnstride + (global_k+0)*knstride; global_node_list.push_back(nn);
      nn = (global_i+1)*instride + (global_j+0)*jnstride + (global_k+0)*knstride; global_node_list.push_back(nn);
      nn = (global_i+1)*instride + (0+0)*jnstride+(global_k+0)*knstride;      global_node_list.push_back(nn);
      nn = (global_i+0)*instride + (0+0)*jnstride+(global_k+0)*knstride;      global_node_list.push_back(nn);

      nn = (global_i+0)*instride + (global_j+0)*jnstride + (global_k+1)*knstride; global_node_list.push_back(nn);
      nn = (global_i+1)*instride + (global_j+0)*jnstride + (global_k+1)*knstride; global_node_list.push_back(nn);
      nn = (global_i+1)*instride + (0+0)*jnstride+(global_k+1)*knstride;      global_node_list.push_back(nn);
      nn = (global_i+0)*instride + (0+0)*jnstride+(global_k+1)*knstride;      global_node_list.push_back(nn);
      }
    }
    else{
      if(dimension == 2){
      nn = (global_i+0)*instride + (global_j+0)*jnstride; global_node_list.push_back(nn);
      nn = (global_i+1)*instride + (global_j+0)*jnstride; global_node_list.push_back(nn);
      nn = (global_i+1)*instride + (global_j+1)*jnstride; global_node_list.push_back(nn);
      nn = (global_i+0)*instride + (global_j+1)*jnstride; global_node_list.push_back(nn);
      }
      else{
      nn = (global_i+0)*instride + (global_j+0)*jnstride + (global_k+0)*knstride; global_node_list.push_back(nn);
      nn = (global_i+1)*instride + (global_j+0)*jnstride + (global_k+0)*knstride; global_node_list.push_back(nn);
      nn = (global_i+1)*instride + (global_j+1)*jnstride + (global_k+0)*knstride; global_node_list.push_back(nn);
      nn = (global_i+0)*instride + (global_j+1)*jnstride + (global_k+0)*knstride; global_node_list.push_back(nn);

      nn = (global_i+0)*instride + (global_j+0)*jnstride + (global_k+1)*knstride; global_node_list.push_back(nn);
      nn = (global_i+1)*instride + (global_j+0)*jnstride + (global_k+1)*knstride; global_node_list.push_back(nn);
      nn = (global_i+1)*instride + (global_j+1)*jnstride + (global_k+1)*knstride; global_node_list.push_back(nn);
      nn = (global_i+0)*instride + (global_j+1)*jnstride + (global_k+1)*knstride; global_node_list.push_back(nn);
      }
    }
  }
  global_node_list.sort();
  global_node_list.unique();


  // Create the global_node_map
  std::list <long long> ::iterator nit;
  long long total = 0;
  for(nit = global_node_list.begin();nit != global_node_list.end();nit ++,total++){
    global_node_vector.push_back(*nit);
    global_node_map[*nit] = total;
  }


  // Create the global_element_map
  long long total_element_count = 0;
  for(unsigned bct = 0; bct < numBlocks();bct ++ ){
    for(unsigned elct = 0;elct < element_block_lists[bct].size();elct++,total_element_count++){
      long long the_el = element_block_lists[bct][elct];
      global_element_map[the_el] = total_element_count;
    }
  }

  if( Debug_Location())
    std::cout << "Inline_Mesh_Desc::Leaving - Build_Global_Lists()" << std::endl;
}

/****************************************************************************/
void Inline_Mesh_Desc::Populate_Nodeset_Info(long long * const * node_set_nodes,
					     std::map <long long, long long> & global_node_map)
/****************************************************************************/
{
  std::list < PG_BC_Specification *> ::iterator setit;
  long long nsct = 0;
  for(setit = nodeset_list.begin(); setit != nodeset_list.end();setit++,nsct ++){
    long long * the_nodes = node_set_nodes[nsct];
    for(unsigned elct = 0; elct < nodeset_vectors[nsct].size();elct++){
      the_nodes[elct] = get_map_entry(global_node_map,nodeset_vectors[nsct][elct])+1;// add one for exodus convention
    }
  }
}

/****************************************************************************/
long long Inline_Mesh_Desc::Populate_Sideset_Info(std::map <long long, long long> & global_element_map,
					  std::map <long long, long long> & global_node_map,
					  long long * const * side_set_elements,
					  long long * const * side_set_faces,
					  long long * const * side_set_nodes,
					  long long * const * side_set_node_counter)
/****************************************************************************/
{
  long long num_nodes_per_face = 2;
  if(dimension == 3){
    num_nodes_per_face = 4;
  }

  long long nsct = 0;
   std::list < PG_BC_Specification *> ::iterator setit;

  for(setit = sideset_list.begin(); setit != sideset_list.end();setit++,nsct ++){

    long long * the_elements = side_set_elements[nsct];
    long long * the_faces = side_set_faces[nsct];
    long long * the_nodes = side_set_nodes[nsct];
    long long * the_node_counter = side_set_node_counter[nsct];
    

    for(unsigned elct = 0; elct < sideset_vectors[nsct].size();elct ++){
      long long the_element = sideset_vectors[nsct][elct].first;
      Topo_Loc the_location =  sideset_vectors[nsct][elct].second;

      // These are the indices of the element in the entire domain
      long long gl_k = the_element/(kestride);
      long long gl_j = (the_element-gl_k*(kestride))/(jestride);
      long long gl_i = the_element - gl_k*(kestride)-gl_j*(jestride);

      //The result of the calculation is the id of the element in the block as found in the connectivity array
      the_elements[elct] = get_map_entry(global_element_map,the_element)+1;
      if(!getErrorString().empty()){return 1;}

      the_faces[elct] = topo_loc_to_exo_face[the_location];
      the_node_counter[elct] = num_nodes_per_face;
      //Since the nodes are numbered across blocks, the global indices are the actual indices of the node
      //It is only required to permute the indices to collect the appropriate nodes for the indicated face and 
      //calculated the node index using the nodal stride values for the entire mesh and the indices

      // adjust for periodicity in j 
      long long glj_plus1 = gl_j + 1;
      if(periodic_j){
        if(glj_plus1 == nel_tot[1]){
          glj_plus1 = 0;
        } 
      }

      switch(the_location) {

        case MINUS_I:{
      if(dimension == 2){
	  the_nodes[elct*num_nodes_per_face + 0] =1+ get_map_entry(global_node_map,(gl_i+0)*instride + (glj_plus1)*jnstride + (gl_k + 0)*knstride);
	  the_nodes[elct*num_nodes_per_face + 1] =1+ get_map_entry(global_node_map,(gl_i+0)*instride + (gl_j + 0)*jnstride + (gl_k + 0)*knstride);
      }
      else{
	  the_nodes[elct*num_nodes_per_face + 0] =1+ get_map_entry(global_node_map,(gl_i+0)*instride + (gl_j + 0)*jnstride + (gl_k + 0)*knstride);
	  the_nodes[elct*num_nodes_per_face + 1] =1+ get_map_entry(global_node_map,(gl_i+0)*instride + (gl_j + 0)*jnstride + (gl_k + 1)*knstride);
	  the_nodes[elct*num_nodes_per_face + 2] =1+ get_map_entry(global_node_map,(gl_i+0)*instride + (glj_plus1)*jnstride + (gl_k + 1)*knstride);
	  the_nodes[elct*num_nodes_per_face + 3] =1+ get_map_entry(global_node_map,(gl_i+0)*instride + (glj_plus1)*jnstride + (gl_k + 0)*knstride);
      }
	  break;
	}
      case PLUS_I:{
      if(dimension == 2){
	the_nodes[elct*num_nodes_per_face + 0] =1+ get_map_entry(global_node_map,(gl_i+1)*instride + (gl_j + 0)*jnstride + (gl_k + 0)*knstride);
	the_nodes[elct*num_nodes_per_face + 1] =1+ get_map_entry(global_node_map,(gl_i+1)*instride + (glj_plus1)*jnstride + (gl_k + 0)*knstride);
      }
      else{
	the_nodes[elct*num_nodes_per_face + 0] =1+ get_map_entry(global_node_map,(gl_i+1)*instride + (gl_j + 0)*jnstride + (gl_k + 0)*knstride);
	the_nodes[elct*num_nodes_per_face + 1] =1+ get_map_entry(global_node_map,(gl_i+1)*instride + (glj_plus1)*jnstride + (gl_k + 0)*knstride);
	the_nodes[elct*num_nodes_per_face + 2] =1+ get_map_entry(global_node_map,(gl_i+1)*instride + (glj_plus1)*jnstride + (gl_k + 1)*knstride);
	the_nodes[elct*num_nodes_per_face + 3] =1+ get_map_entry(global_node_map,(gl_i+1)*instride + (gl_j + 0)*jnstride + (gl_k + 1)*knstride);
      }
	break;
      }
      case MINUS_J:{
      if(dimension == 2){
	the_nodes[elct*num_nodes_per_face + 0] =1+ get_map_entry(global_node_map,(gl_i+0)*instride + (gl_j + 0)*jnstride + (gl_k + 0)*knstride);
	the_nodes[elct*num_nodes_per_face + 1] =1+ get_map_entry(global_node_map,(gl_i+1)*instride + (gl_j + 0)*jnstride + (gl_k + 0)*knstride);
      }
      else{
	the_nodes[elct*num_nodes_per_face + 0] =1+ get_map_entry(global_node_map,(gl_i+0)*instride + (gl_j + 0)*jnstride + (gl_k + 0)*knstride);
	the_nodes[elct*num_nodes_per_face + 1] =1+ get_map_entry(global_node_map,(gl_i+1)*instride + (gl_j + 0)*jnstride + (gl_k + 0)*knstride);
	the_nodes[elct*num_nodes_per_face + 2] =1+ get_map_entry(global_node_map,(gl_i+1)*instride + (gl_j + 0)*jnstride + (gl_k + 1)*knstride);
	the_nodes[elct*num_nodes_per_face + 3] =1+ get_map_entry(global_node_map,(gl_i+0)*instride + (gl_j + 0)*jnstride + (gl_k + 1)*knstride);
      }
	break;
      }
      case PLUS_J:{
      if(dimension == 2){
	the_nodes[elct*num_nodes_per_face + 0] =1+ get_map_entry(global_node_map,(gl_i+1)*instride + (glj_plus1)*jnstride + (gl_k + 0)*knstride);
	the_nodes[elct*num_nodes_per_face + 1] =1+ get_map_entry(global_node_map,(gl_i+0)*instride + (glj_plus1)*jnstride + (gl_k + 0)*knstride);
      }
      else{
	the_nodes[elct*num_nodes_per_face + 0] =1+ get_map_entry(global_node_map,(gl_i+1)*instride + (glj_plus1)*jnstride + (gl_k + 0)*knstride);
	the_nodes[elct*num_nodes_per_face + 1] =1+ get_map_entry(global_node_map,(gl_i+0)*instride + (glj_plus1)*jnstride + (gl_k + 0)*knstride);
	the_nodes[elct*num_nodes_per_face + 2] =1+ get_map_entry(global_node_map,(gl_i+0)*instride + (glj_plus1)*jnstride + (gl_k + 1)*knstride);
	the_nodes[elct*num_nodes_per_face + 3] =1+ get_map_entry(global_node_map,(gl_i+1)*instride + (glj_plus1)*jnstride + (gl_k + 1)*knstride);
      }
	break;
      }
      case MINUS_K:{
      if(dimension == 2){
      }
      else{
	the_nodes[elct*num_nodes_per_face + 0] =1+ get_map_entry(global_node_map,(gl_i+0)*instride + (gl_j + 0)*jnstride + (gl_k + 0)*knstride);
	the_nodes[elct*num_nodes_per_face + 1] =1+ get_map_entry(global_node_map,(gl_i+0)*instride + (glj_plus1)*jnstride + (gl_k + 0)*knstride);
	the_nodes[elct*num_nodes_per_face + 2] =1+ get_map_entry(global_node_map,(gl_i+1)*instride + (glj_plus1)*jnstride + (gl_k + 0)*knstride);
	the_nodes[elct*num_nodes_per_face + 3] =1+ get_map_entry(global_node_map,(gl_i+1)*instride + (gl_j + 0)*jnstride + (gl_k + 0)*knstride);
      }
	break;
      }
      case PLUS_K:{
      if(dimension == 2){
      }
      else{
	the_nodes[elct*num_nodes_per_face + 0] =1+ get_map_entry(global_node_map,(gl_i+0)*instride + (gl_j + 0)*jnstride + (gl_k + 1)*knstride);
	the_nodes[elct*num_nodes_per_face + 1] =1+ get_map_entry(global_node_map,(gl_i+1)*instride + (gl_j + 0)*jnstride + (gl_k + 1)*knstride);
	the_nodes[elct*num_nodes_per_face + 2] =1+ get_map_entry(global_node_map,(gl_i+1)*instride + (glj_plus1)*jnstride + (gl_k + 1)*knstride);
	the_nodes[elct*num_nodes_per_face + 3] =1+ get_map_entry(global_node_map,(gl_i+0)*instride + (glj_plus1)*jnstride + (gl_k + 1)*knstride);
      }
	break;		     
      }
      default:
	error_stream << "Inline_Mesh_Desc::Read_mesh(): "
		     << "Sideset applied to unknown Topo Location.";
	return 1;
	break;
      }   
    }
  }
  //END of SIDESETS
  return 0;
}




/****************************************************************************/
void Inline_Mesh_Desc::Populate_Connectivity(long long * const * conn_array,                                         
					  std::map <long long, long long> & global_node_map
)
/****************************************************************************/
{
  long long num_nodes_per_element = 4;
  if(dimension == 3){
    num_nodes_per_element = 8;
  }
  //Nodes are ordered across the entire block domain.
  //To calculate connectivity for a given element in a given block.
  //It is necessary to calculate the global element indices in i,j,k
  //that identify that element in the global space.
  long long total_element_count = 0;
  for(long long bct = 0; bct < numBlocks(); bct ++ ){
    long long * conn = conn_array[bct];
    //build connectivity for each element
    //nodes are numbered 1-tot_num_nodes+1 across all blocks
    //incrementing i fastest
    for(unsigned elct = 0;elct < element_block_lists[bct].size();elct++,total_element_count++){
      long long the_el = element_block_lists[bct][elct];
      long long Kg;
      long long Jg;
      long long Ig;
      long long l;
      get_l_i_j_k_from_element_number(the_el,l,Ig,Jg,Kg);


      conn[elct*num_nodes_per_element + 0 + 0] = get_map_entry(global_node_map,get_node_number_from_l_i_j_k(l, Ig + 0, Jg + 0, Kg + 0))+1;
      conn[elct*num_nodes_per_element + 1 + 0] = get_map_entry(global_node_map,get_node_number_from_l_i_j_k(l, Ig + 1, Jg + 0, Kg + 0))+1;
      conn[elct*num_nodes_per_element + 2 + 0] = get_map_entry(global_node_map,get_node_number_from_l_i_j_k(l, Ig + 1, Jg + 1, Kg + 0))+1;
      conn[elct*num_nodes_per_element + 3 + 0] = get_map_entry(global_node_map,get_node_number_from_l_i_j_k(l, Ig + 0, Jg + 1, Kg + 0))+1;
      if(dimension == 3){
	conn[elct*num_nodes_per_element + 0 + 4] = get_map_entry(global_node_map,get_node_number_from_l_i_j_k(l, Ig + 0, Jg + 0, Kg + 1))+1;
	conn[elct*num_nodes_per_element + 1 + 4] = get_map_entry(global_node_map,get_node_number_from_l_i_j_k(l, Ig + 1, Jg + 0, Kg + 1))+1;
	conn[elct*num_nodes_per_element + 2 + 4] = get_map_entry(global_node_map,get_node_number_from_l_i_j_k(l, Ig + 1, Jg + 1, Kg + 1))+1;
	conn[elct*num_nodes_per_element + 3 + 4] = get_map_entry(global_node_map,get_node_number_from_l_i_j_k(l, Ig + 0, Jg + 1, Kg + 1))+1;
      }
    }
  }
}

/****************************************************************************/
void Inline_Mesh_Desc::Populate_Map_and_Global_Element_List(long long * the_map, 
							 long long * global_element_numbers)
/****************************************************************************/
{
  long long total_count = 0;
  for(long long bct = 0; bct < numBlocks();bct ++ ){
    for(unsigned ect = 0; ect < element_block_lists[bct].size();ect ++){
      long long the_el = element_block_lists[bct][ect];
      //global element indices
      long long Kg = the_el/kestride;
      long long Jg = (the_el - Kg*kestride)/(jestride);
      long long Ig = the_el  - Kg*kestride - Jg*jestride;

      long long Kbg = get_block_index(Kg,inline_b[2],c_inline_n[2]);
      long long Jbg = get_block_index(Jg,inline_b[1],c_inline_n[1]);
      long long Ibg = get_block_index(Ig,inline_b[0],c_inline_n[0]);


      //ordinal of the block
      long long the_block = Ibg + Jbg*inline_b[0] + Kbg*blockKstride();

      //indices inside the block
      long long Kblock = Kg-c_inline_n[2][Kbg];
      long long Jblock = Jg-c_inline_n[1][Jbg];
      long long Iblock = Ig-c_inline_n[0][Ibg];


      //product
      //       long long elid = the_block*inline_nx*inline_ny*inline_nz + Iblock + Jblock*inline_nx + Kblock*inline_nx*inline_ny;
      long long elid = cum_block_totals[the_block] +
        Iblock + Jblock*a_inline_n[0][Ibg] + Kblock*a_inline_n[0][Ibg]*a_inline_n[1][Jbg];

      the_map[total_count] = elid + 1;
      global_element_numbers[total_count] = elid + 1;
      total_count ++;
    }
  }
}

/****************************************************************************/
long long Inline_Mesh_Desc::Check_Spans()
/****************************************************************************/
{
  if((inline_gmax[0] - inline_gmin[0]) <= 0.){
    error_stream << "Invalid span for 'X' range of inline mesh specification. The span must be positive";
    return 1;
  }
  
  if((inline_gmax[1] - inline_gmin[1]) <= 0.){
    error_stream << "Invalid span for 'Y' range of inline mesh specification. The span must be positive";
    return 1;
  }
  if(dimension == 3){
    if((inline_gmax[2]-inline_gmin[2]) <= 0.){
      error_stream << "Invalid span for 'Z' range of inline mesh specification. The span must be positive";
      return 1;
    }
  }
  return 0;
}

/****************************************************************************/
long long Inline_Mesh_Desc::Check_Block_BC_Sets()
/****************************************************************************/
{
  std::list < PG_BC_Specification * > ::iterator setit;
  for(setit = nodeset_list.begin(); setit != nodeset_list.end();setit++){
    for(unsigned ict = 0;ict < (*setit)->the_locs.size();ict ++){
      if((*setit)->the_locs[ict].block_boundary_set){
	long long bid = (*setit)->the_locs[ict].block_id;
	long long bmax = numBlocks();
	if (bid < 1 || bid > bmax){
	  error_stream << "Terminating from Inline_Mesh_Desc::Check_Block_BC_Sets,block index ";
	  error_stream << bid ;
	  error_stream << " is outside the range of blocks present in the mesh  1 to ";
	  error_stream << bmax;
	  error_stream << ".";
	  return 1;
	}
	/*check if block suppressed*/
	if(isBlockSuppressed(bid)){
	  error_stream << "Terminating from Inline_Mesh_Desc::Check_Block_BC_Sets,block index ";
	  error_stream << bid ;
	  error_stream << " is suppressed and may not accept nodesets.";
	  return 1;
	}
      }
    }
  }

  for(setit = sideset_list.begin(); setit != sideset_list.end();setit++){
    for(unsigned ict = 0;ict < (*setit)->the_locs.size();ict ++){
      if((*setit)->the_locs[ict].block_boundary_set){
	long long bid = (*setit)->the_locs[ict].block_id;
	long long bmax = numBlocks();
	if (bid < 1 || bid > bmax){
	  error_stream << "Terminating from Inline_Mesh_Desc::Check_Block_BC_Sets,block index ";
	  error_stream << bid ;
	  error_stream << " is outside the range of blocks present in the mesh  1 to ";
	  error_stream << bmax;
	  error_stream << ".";
	  return 1;
	}
	/*check if block suppressed*/
	if(isBlockSuppressed(bid)){
	  error_stream << "Terminating from Inline_Mesh_Desc::Check_Block_BC_Sets,block index ";
	  error_stream << bid ;
	  error_stream << " is suppressed and may not accept sidesets.";
	  return 1;
	}
      }
    }
  }

  return Rename_Block_BC_Sets();
}

/****************************************************************************/
long long Inline_Mesh_Desc::Check_Blocks()
/****************************************************************************/
{
  if(inline_b[0] == 0 || inline_b[1] == 0 || inline_b[2] == 0){
    error_stream << "Terminating from Inline_Mesh_Desc::Check_Blocks, zero value found,";
    error_stream << " inline_b[0] " << inline_b[0];
    error_stream << " inline_b[1] " << inline_b[1];
    error_stream << " inline_b[2] " << inline_b[2];
    return 1;
  }
  /*check suppressed blocks are properly clled out*/
  long long bmax = numBlocks();
  std::set<long long >::iterator sit;
  for(sit = suppressed_blocks.begin();sit != suppressed_blocks.end();sit ++){
    if(*sit < 0 || (*sit)>bmax){
      error_stream << "Terminating from Inline_Mesh_Desc::Check_Blocks block ";
      error_stream << *sit ;
      error_stream << " may not be suppressed as it does not exist.";
      return 1;
    }
  }

  return 0;
}

/****************************************************************************/
void Inline_Mesh_Desc::Size_BC_Sets(long long nnx, 
				 long long nny, 
				 long long nnz)
/****************************************************************************/
{
  //Nodesets
  std::list < PG_BC_Specification * > ::iterator setit;
  for(setit = nodeset_list.begin(); setit != nodeset_list.end();setit++){
    for(unsigned ict = 0;ict < (*setit)->the_locs.size();ict ++){
      Topo_Loc the_location = (*setit)->the_locs[ict].location;
      if((*setit)->the_locs[ict].block_boundary_set){
	long long bid = (*setit)->the_locs[ict].block_id-1;
	long long kind = bid/(inline_b[0]*inline_b[1]);
	long long jind = (bid-kind*(inline_b[0] * inline_b[1]))/inline_b[0];
	long long iind = bid - jind *(inline_b[0]) - kind*(inline_b[0] * inline_b[1]);
	
	(*setit)->the_locs[ict].limits = getLimits(the_location,
						   c_inline_n[0][iind],c_inline_n[0][iind+1]+1,
						   c_inline_n[1][jind],c_inline_n[1][jind+1]+1,
						   c_inline_n[2][kind],c_inline_n[2][kind+1]+1,
						   nnx, 
						   nny);
      }
      else{
	(*setit)->the_locs[ict].limits = getLimits(the_location,
						   0,nnx,
						   0,nny,
						   0,nnz,
						   nnx,
						   nny);
      }
    }
  }
  //Sidesets
  for(setit = sideset_list.begin(); setit != sideset_list.end();setit++){
    for(unsigned ict = 0;ict < (*setit)->the_locs.size();ict ++){
      Topo_Loc the_location = (*setit)->the_locs[ict].location;
      if((*setit)->the_locs[ict].block_boundary_set){
	long long bid = (*setit)->the_locs[ict].block_id-1;
	long long kind = bid/(inline_b[0]*inline_b[1]);
	long long jind = (bid-kind*(inline_b[0] * inline_b[1]))/inline_b[0];
	long long iind = bid - jind *(inline_b[0]) - kind*(inline_b[0] * inline_b[1]);
	(*setit)->the_locs[ict].limits = getLimits(the_location,
						   c_inline_n[0][iind],c_inline_n[0][iind+1],
						   c_inline_n[1][jind],c_inline_n[1][jind+1],
						   c_inline_n[2][kind],c_inline_n[2][kind+1],
						   nel_tot[0], 
						   nel_tot[1]);
      }
      else{
	(*setit)->the_locs[ict].limits = getLimits(the_location,
						   0,nel_tot[0],
						   0,nel_tot[1],
						   0,nel_tot[2],
						   nel_tot[0],
						   nel_tot[1]);
      } 
    }
  }
}

/******************************************************************************/
LoopLimits Inline_Mesh_Desc::getLimits( Topo_Loc the_set_location,
				     long long sx, long long nx, 
				     long long sy, long long ny, 
				     long long sz, long long nz, 
				     long long irange, long long jrange)
/******************************************************************************/
{
  LoopLimits ll;
  long long istart = sx;
  long long iend = nx;
  long long jstart = sy;
  long long jend = ny;
  long long kstart = sz;
  long long kend = nz;

  switch(the_set_location) {
  case ALL_NODES:{
    break;
  }
  case MINUS_I:{
    iend = istart + 1;
    break;
  }
  case PLUS_I:{
    istart = iend - 1;
    break;
  }
  case MINUS_J:{
    jend = jstart + 1;
    break;
  }
  case PLUS_J:{
    jstart = jend -1;
    break;
  }
  case MINUS_K:{
    kend = kstart + 1;
    break;
  }
  case PLUS_K:{
    kstart = kend -1;
    break;		     
  }
  case EDGE0:{
    jend = jstart + 1;
    kend = kstart + 1;
    break;		     
  }
  case EDGE1:{
    istart = iend - 1;
    kend = kstart + 1;
    break;		     
  }
  case EDGE2:{
    jstart = jend - 1;
    kend = kstart + 1;
    break;		     
  }
  case EDGE3:{
    iend = istart + 1;
    kend = kstart + 1;
    break;		     
  }
  case EDGE4:{
    iend = istart + 1;
    jend = jstart + 1;
    break;		     
  }
  case EDGE5:{
    istart = iend - 1;
    jend = jstart + 1;
    break;		     
  }
  case EDGE6:{
    istart = iend - 1;
    jstart = jend - 1;
    break;		     
  }
  case EDGE7:{
    iend = istart + 1;
    jstart = jend - 1;
    break;		     
  }
  case EDGE8:{
    jend = jstart + 1;
    kstart = kend - 1;
    break;		     
  }
  case EDGE9:{
    kstart = kend - 1;
    istart = iend - 1;
    break;		     
  }
  case EDGE10:{
    jstart = jend - 1;
    kstart = kend - 1;
    break;		     
  }
  case EDGE11:{
    iend = istart + 1;
    kstart = kend - 1;
    break;		     
  }
  case VERTEX0:{
    iend = istart + 1;
    jend = jstart + 1;
    kend = kstart + 1;
    break;		     
  }
  case VERTEX1:{
    istart = iend - 1;
    jend = jstart + 1;
    kend = kstart + 1;
    break;		     
  }
  case VERTEX2:{
    istart = iend - 1;
    jstart = jend - 1;
    kend = kstart + 1;
    break;		     
  }
  case VERTEX3:{
    iend = istart + 1;
    jstart = jend - 1;
    kend = kstart + 1;
    break;		     
  }
  case VERTEX4:{
    iend = istart + 1;
    jend = jstart + 1;
    kstart = kend - 1;
    break;		     
  }
  case VERTEX5:{
    istart = iend - 1;
    jend = jstart + 1;
    kstart = kend - 1;
    break;		     
  }
  case VERTEX6:{
    istart = iend - 1;
    jstart = jend - 1;
    kstart = kend - 1;
    break;		     
  }
  case VERTEX7:{
    iend = istart + 1;
    jstart = jend - 1;
    kstart = kend - 1;
    break;		     
  }
  case NUM_TOPO_CONNECTIONS:{
    break;		     
  }
  case PROCESSOR_NODE:{
    break;		     
  }
  default:
    iend = istart;
    jend = jstart;
    kend = kstart;
    break;
  }
  ll.is = istart;
  ll.ie = iend;
  ll.js = jstart;
  ll.je = jend;
  ll.jstride = irange;
  ll.total = (iend-istart)*(jend-jstart);
  if(dimension == 3){
    ll.ks = kstart;
    ll.ke = kend;
    ll.kstride = irange*jrange;
    ll.total = (iend-istart)*(jend-jstart)*(kend-kstart);
  }
  return ll;
}




/*****************************************************************************/
void Inline_Mesh_Desc::ZeroSet()
/*****************************************************************************/
{ 


  inline_geometry_type = UNKNOWN;
  inline_decomposition_type = BISECTION;
  trisection_blocks = 0;
  inline_b[0] = 0;
  inline_b[1] = 0;
  inline_b[2] = 1;
  inline_n[0] = 1;
  inline_n[1] = 1;
  inline_n[2] = 1;
  inline_block_start = 1;
  inline_offset[0] = 0.;
  inline_offset[1] = 0.;
  inline_offset[2] = 0.;
  inline_gmin[0] = 0.;
  inline_gmin[1] = 0.;
  inline_gmin[2] = 0.;
  inline_gmax[0] = 0.;
  inline_gmax[1] = 0.;
  inline_gmax[2] = 0.;
  inc_nels[0] = 0;
  inc_nels[1] = 0;
  inc_nels[2] = 0;
  inc_nocuts[0] = false;
  inc_nocuts[1] = false;
  inc_nocuts[2] = false;
  inline_nprocs[0] = 1;
  inline_nprocs[1] = 1;
  inline_nprocs[2] = 1;
  periodic_i = false;
  periodic_j = false;
  periodic_k = false;
  try_squared = false;
  enforce_periodic = false;
  Element_Density_Functions[0] = NULL;
  Element_Density_Functions[1] = NULL;
  Element_Density_Functions[2] = NULL;
  Geometry_Transform_Function = NULL;
  a_inline_n[0] = NULL;
  a_inline_n[1] = NULL;
  a_inline_n[2] = NULL;
  c_inline_n[0] = NULL;
  c_inline_n[1] = NULL;
  c_inline_n[2] = NULL;
  cum_block_totals = NULL;
  els_in_block = NULL;
  nel_tot[0] = 0;
  nel_tot[1] = 0;
  nel_tot[2] = 0;
  block_dist = new double * [3];
  c_block_dist = new double * [3];
  first_size = new double * [3];
  last_size = new double * [3];
  interval = new long long * [3];
  for(long long i = 0; i < 3; i ++){
    block_dist[i] = NULL;
    c_block_dist[i] = NULL;
    first_size[i] = NULL;
    last_size[i] = NULL;
    interval[i] = NULL;
  }

    IJKcoors[0] = NULL;
    IJKcoors[1] = NULL;
    IJKcoors[2] = NULL;
    base_partition = NULL;
    
    transition_radius = -1;

    total_unsupressed_elements = 0;

  my_rank = 0;
  num_processors = 1;

  topo_loc_to_exo_face[MINUS_I] = 4;
  topo_loc_to_exo_face[PLUS_I] = 2;
  topo_loc_to_exo_face[MINUS_J] = 1;
  topo_loc_to_exo_face[PLUS_J] = 3;
  topo_loc_to_exo_face[MINUS_K] = 5;
  topo_loc_to_exo_face[PLUS_K] = 6;

  sideset_vectors = NULL;
  sideset_global_count = NULL;
  nodeset_vectors = NULL;

  debug_mode = false;
  next = NULL;

}


/*****************************************************************************/
void Inline_Mesh_Desc::Display_Class(std::ostream& s, const std::string &indent)
/*****************************************************************************/
{
  for(long long i = 0; i < dimension; i++){
    if(Element_Density_Functions[i])Element_Density_Functions[i]->Display_Class(s,indent);
  }
}

/*****************************************************************************/
void Inline_Mesh_Desc::Populate_Border_Nodes_Elements( long long * internal_elements,
						       long long * internal_nodes,
						       long long * border_elements,
						       long long * border_nodes,
						       std::list   <long long> & internal_node_list,	
						       std::list   <long long> & border_nodes_list,
						       std::list   <long long> & internal_element_list,
						       std::list   <long long> & border_elements_list,
						       std::map <long long, long long> & global_node_map,
						       std::map <long long, long long> & global_element_map)
/*****************************************************************************/
{
  std::list <long long> :: iterator eit;
  long long the_count = 0;
  for(eit = internal_element_list.begin();eit != internal_element_list.end();eit ++,the_count++)
    internal_elements[the_count]= get_map_entry(global_element_map,(*eit))+1;
  the_count = 0;
  for(eit = internal_node_list.begin();eit != internal_node_list.end();eit ++,the_count++)
    internal_nodes[the_count]= get_map_entry(global_node_map,(*eit))+1;
  the_count = 0;
  for(eit = border_elements_list.begin();eit != border_elements_list.end();eit ++,the_count++)
    border_elements[the_count]= get_map_entry(global_element_map,(*eit))+1;
  the_count = 0;
  for(eit = border_nodes_list.begin();eit != border_nodes_list.end();eit ++,the_count++)
    border_nodes[the_count]= get_map_entry(global_node_map,(*eit))+1;
}


/*****************************************************************************/
void Inline_Mesh_Desc::Populate_Cmap( long long * node_cmap_node_cnts,
				      long long * node_cmap_ids,
				      long long * elem_cmap_elem_cnts,
				      long long * elem_cmap_ids,
				      std::vector <long long> & node_neighbor_vector,
				      std::vector <long long> & element_neighbor_vector,
				      std::list <long long>  * & boundary_node_list,                   
				      std::list < std::pair <long long ,Topo_Loc > > * & boundary_element_list)
/*****************************************************************************/
{
  for(unsigned i = 0; i < node_neighbor_vector.size();i++){
    node_cmap_node_cnts[i] = boundary_node_list[i].size();
    node_cmap_ids[i] = node_neighbor_vector[i];
  }
  for(unsigned i = 0; i < element_neighbor_vector.size();i++){
    elem_cmap_elem_cnts[i] = boundary_element_list[i].size();
    elem_cmap_ids[i] = element_neighbor_vector[i];
  }
}

/*****************************************************************************/
void Inline_Mesh_Desc::Populate_Parallel_Info( long long * const * comm_node_ids ,
					       long long * const * comm_node_proc_ids,
					       long long * const * comm_elem_ids,
					       long long * const * comm_side_ids,
					       long long * const * comm_elem_proc_ids,
					       std::vector <long long> & node_neighbor_vector,
					       std::vector <long long> & element_neighbor_vector,		
					       std::list <long long>  * & boundary_node_list,
					       std::map <long long, long long> & global_node_map,
					       std::list <std::pair <long long ,Topo_Loc > > * & boundary_element_list,
					       std::map <long long, long long> & global_element_map)
/*****************************************************************************/
{
  for(unsigned i = 0; i < node_neighbor_vector.size();i++){
    long long * comm_nodes = comm_node_ids[i];
    long long * comm_node_procs = comm_node_proc_ids[i];
    std::list <long long> :: iterator nlit;
    long long nct = 0;
    for(nlit = boundary_node_list[i].begin();nlit != boundary_node_list[i].end();nlit++,nct ++){
      comm_nodes[nct] = get_map_entry(global_node_map,(*nlit))+1;
      comm_node_procs[nct] = node_neighbor_vector[i];// is this right?
    }
  }
  for(unsigned i = 0; i < element_neighbor_vector.size();i++){
    long long * comm_elements = comm_elem_ids[i];
    long long * comm_sides = comm_side_ids[i];
    long long * comm_elem_procs = comm_elem_proc_ids[i];
    std::list < std::pair <long long ,Topo_Loc > > :: iterator nlit;
    long long nct = 0;
    for(nlit = boundary_element_list[i].begin();nlit != boundary_element_list[i].end();nlit++,nct ++){
      comm_elements[nct] = get_map_entry(global_element_map,(*nlit).first)+1;//foo bar
      comm_sides[nct] = topo_loc_to_exo_face[(*nlit).second];
      comm_elem_procs[nct] = element_neighbor_vector[i];// is this right?
    }
  }
}

/*****************************************************************************/
void Inline_Mesh_Desc::setStrides()
/*****************************************************************************/
{
  instride = 1;
  jnstride = nel_tot[0]+1;
  knstride = (nel_tot[0]+1)*(nel_tot[1]+1);
  
  if(inline_geometry_type == RADIAL && periodic_j){
    instride = 1;
    jnstride = nel_tot[0]+1;
    knstride = (nel_tot[0]+1)*(nel_tot[1]);
  }
  
  iestride = 1;
  jestride = nel_tot[0];
  kestride = (nel_tot[0])*(nel_tot[1]);
  
  for(long long i = 0; i < kestride*nel_tot[2];i++){
    if(!isElementSuppressed(i))total_unsupressed_elements ++;
  }
}

//! A utility function to build up required bookkeeping objects.
/****************************************************************************/
void Inline_Mesh_Desc::get_l_i_j_k_from_element_number(long long el,
								     long long & l,
								     long long & i,
								     long long & j,
								     long long & k)
/****************************************************************************/
{
  l = 0;
  k = el/kestride;
  long long remainder = el-k*kestride;
  
  j = (remainder)/(jestride);
  i = remainder - j*(jestride);

}

//! A utility function to build up required bookkeeping objects.
/****************************************************************************/
void Inline_Mesh_Desc::get_l_i_j_k_from_node_number(long long nn,
								 long long & /* l */,
								 long long & i,
								 long long & j,
								 long long & k)
/****************************************************************************/
{
  k = nn/knstride;
  long long remainder = nn - k*knstride;  
  j = remainder / (nel_tot[0]+1);
  i = remainder - j*(nel_tot[0]+1);
}

/****************************************************************************/
long long Inline_Mesh_Desc::get_element_number_from_l_i_j_k( long long  /* l */,
						       long long  i,
						       long long  j,
						       long long  k)
/****************************************************************************/
{
  //adding bounds checking to deal with neighbor calcultion
  // return -1 if off mesh
  // adjust and recurse if i or j throw element onto adjacent block
  if(k < 0) return -1;
  if(k >= nel_tot[2]) return -1;
  if(periodic_j){
    if(j < 0) j = nel_tot[1]-1;
    if(j >= nel_tot[1])j = 0;
  }
  else{
  if(j < 0) return -1;
  if(j >= nel_tot[1])return -1;
  }

  if(i<0)return -1;
  if(i >= nel_tot[0])return -1;

  long long elno;
  elno = k*kestride;
  elno += j*jestride;
  elno += i;

  return elno;
}

/****************************************************************************/
long long Inline_Mesh_Desc::get_node_number_from_l_i_j_k( long long  /* l */,
							  long long  i,
							  long long  j,
							  long long  k)
  /****************************************************************************/
{
  if(periodic_j){
    if( j == nel_tot[1])j = 0;
  }
  long long nno = k*knstride;
  nno += j*jnstride;
  nno += i;
  return nno;
}


/****************************************************************************/
long long Inline_Mesh_Desc::get_neighbor(Topo_Loc tl,
				   long long ll, 
				   long long li, 
				   long long lj, 
				   long long lk)
/****************************************************************************/
{
  long long result = -1;
  switch(tl) {
  
  case MINUS_I:{
    result =  get_element_number_from_l_i_j_k(ll,li-1,lj,lk);
    break;
  }
  case PLUS_I:{
    result =  get_element_number_from_l_i_j_k(ll,li+1,lj,lk);
    break;
  }
  case MINUS_J:{
    result =  get_element_number_from_l_i_j_k(ll,li,lj-1,lk);
    break;
  }
  case PLUS_J:{
    result =  get_element_number_from_l_i_j_k(ll,li,lj+1,lk);
    break;
  }
  case MINUS_K:{
    result =  get_element_number_from_l_i_j_k(ll,li,lj,lk-1);
    break;
  }
  case PLUS_K:{
    result =  get_element_number_from_l_i_j_k(ll,li,lj,lk+1);
    break;
  }
  case EDGE0:{
    result =  get_element_number_from_l_i_j_k(ll,li,lj-1,lk-1);
    break;
  }
  case EDGE1:{
    result =  get_element_number_from_l_i_j_k(ll,li+1,lj,lk-1);
    break;
  }
  case EDGE2:{
    result =  get_element_number_from_l_i_j_k(ll,li,lj+1,lk-1);
    break;
  }
  case EDGE3:{
    result =  get_element_number_from_l_i_j_k(ll,li-1,lj,lk-1);
    break;
  }
  case EDGE4:{
    result =  get_element_number_from_l_i_j_k(ll,li-1,lj-1,lk);
    break;
  }
  case EDGE5:{
    result =  get_element_number_from_l_i_j_k(ll,li+1,lj-1,lk);
    break;
  }
  case EDGE6:{
    result =  get_element_number_from_l_i_j_k(ll,li+1,lj+1,lk);
    break;
  }
  case EDGE7:{
    result =  get_element_number_from_l_i_j_k(ll,li-1,lj+1,lk);
    break;
  }
  case EDGE8:{
    result =  get_element_number_from_l_i_j_k(ll,li,lj-1,lk+1);
    break;
  }
  case EDGE9:{
    result =  get_element_number_from_l_i_j_k(ll,li+1,lj,lk+1);
    break;
  }
  case EDGE10:{
    result =  get_element_number_from_l_i_j_k(ll,li,lj+1,lk+1);
    break;
  }
  case EDGE11:{
    result =  get_element_number_from_l_i_j_k(ll,li-1,lj,lk+1);
    break;
  }
  case VERTEX0:{
    result =  get_element_number_from_l_i_j_k(ll,li-1,lj-1,lk-1);
    break;
  }
  case VERTEX1:{
    result =  get_element_number_from_l_i_j_k(ll,li+1,lj-1,lk-1);
    break;
  }
  case VERTEX2:{
    result =  get_element_number_from_l_i_j_k(ll,li+1,lj+1,lk-1);
    break;
  }
  case VERTEX3:{
    result =  get_element_number_from_l_i_j_k(ll,li-1,lj+1,lk-1);
    break;
  }
  case VERTEX4:{
    result =  get_element_number_from_l_i_j_k(ll,li-1,lj-1,lk+1);
    break;
  }
  case VERTEX5:{
    result =  get_element_number_from_l_i_j_k(ll,li+1,lj-1,lk+1);
    break;
  }
  case VERTEX6:{
    result =  get_element_number_from_l_i_j_k(ll,li+1,lj+1,lk+1);
    break;
  }
  case VERTEX7:{
    result =  get_element_number_from_l_i_j_k(ll,li-1,lj+1,lk+1);
    break;
  }
  default:
    return -1;
    break;
  }
  if(result != -1){
    if(isElementSuppressed(result)){
      result = -1;
    }
  }
  return result;
}

/****************************************************************************/
void Inline_Mesh_Desc::get_face_nodes(Topo_Loc tl,
				      long long global_element_id,
				      long long the_nodes[])
/****************************************************************************/
{
  long long gl_l,gl_i,gl_j,gl_k;
  get_l_i_j_k_from_element_number(global_element_id,gl_l,gl_i,gl_j,gl_k);
  switch(tl) {
    
  case MINUS_I:{
      if(dimension == 2){
	the_nodes[0] = get_node_number_from_l_i_j_k(gl_l, gl_i + 0, gl_j + 1, gl_k + 0);
	the_nodes[1] = get_node_number_from_l_i_j_k(gl_l, gl_i + 0, gl_j + 0, gl_k + 0);
      }
      else{
	the_nodes[0] = get_node_number_from_l_i_j_k(gl_l, gl_i + 0, gl_j + 0, gl_k + 0);
	the_nodes[1] = get_node_number_from_l_i_j_k(gl_l, gl_i + 0, gl_j + 0, gl_k + 1);
	the_nodes[2] = get_node_number_from_l_i_j_k(gl_l, gl_i + 0, gl_j + 1, gl_k + 1);
	the_nodes[3] = get_node_number_from_l_i_j_k(gl_l, gl_i + 0, gl_j + 1, gl_k + 0);
      }
    break;
  }
  case PLUS_I:{
      if(dimension == 2){
	the_nodes[0] = get_node_number_from_l_i_j_k(gl_l, gl_i + 1, gl_j + 0, gl_k + 0);
	the_nodes[1] = get_node_number_from_l_i_j_k(gl_l, gl_i + 1, gl_j + 1, gl_k + 0);
      }    
      else{
	the_nodes[0] = get_node_number_from_l_i_j_k(gl_l, gl_i + 1, gl_j + 0, gl_k + 0);
	the_nodes[1] = get_node_number_from_l_i_j_k(gl_l, gl_i + 1, gl_j + 1, gl_k + 0);
	the_nodes[2] = get_node_number_from_l_i_j_k(gl_l, gl_i + 1, gl_j + 1, gl_k + 1);
	the_nodes[3] = get_node_number_from_l_i_j_k(gl_l, gl_i + 1, gl_j + 0, gl_k + 1);
      }
    break;
  }
  case MINUS_J:{
      if(dimension == 2){
	the_nodes[0] = get_node_number_from_l_i_j_k(gl_l, gl_i + 0, gl_j + 0, gl_k + 0);
	the_nodes[1] = get_node_number_from_l_i_j_k(gl_l, gl_i + 1, gl_j + 0, gl_k + 0);
      }
      else{
	the_nodes[0] = get_node_number_from_l_i_j_k(gl_l, gl_i + 0, gl_j + 0, gl_k + 0);
	the_nodes[1] = get_node_number_from_l_i_j_k(gl_l, gl_i + 1, gl_j + 0, gl_k + 0);
	the_nodes[2] = get_node_number_from_l_i_j_k(gl_l, gl_i + 1, gl_j + 0, gl_k + 1);
	the_nodes[3] = get_node_number_from_l_i_j_k(gl_l, gl_i + 0, gl_j + 0, gl_k + 1);
      }
    break;
  }
  case PLUS_J:{
      if(dimension == 2){
	the_nodes[0] = get_node_number_from_l_i_j_k(gl_l, gl_i + 1, gl_j + 1, gl_k + 0);
	the_nodes[1] = get_node_number_from_l_i_j_k(gl_l, gl_i + 0, gl_j + 1, gl_k + 0);
      }
      else{
	the_nodes[0] = get_node_number_from_l_i_j_k(gl_l, gl_i + 1, gl_j + 1, gl_k + 0);
	the_nodes[1] = get_node_number_from_l_i_j_k(gl_l, gl_i + 0, gl_j + 1, gl_k + 0);
	the_nodes[2] = get_node_number_from_l_i_j_k(gl_l, gl_i + 0, gl_j + 1, gl_k + 1);
	the_nodes[3] = get_node_number_from_l_i_j_k(gl_l, gl_i + 1, gl_j + 1, gl_k + 1);
      }
    break;
  }
  case MINUS_K:{
      if(dimension == 2){
      }
      else{
    the_nodes[0] = get_node_number_from_l_i_j_k(gl_l, gl_i + 0, gl_j + 0, gl_k + 0);
    the_nodes[1] = get_node_number_from_l_i_j_k(gl_l, gl_i + 0, gl_j + 1, gl_k + 0);
    the_nodes[2] = get_node_number_from_l_i_j_k(gl_l, gl_i + 1, gl_j + 1, gl_k + 0);
    the_nodes[3] = get_node_number_from_l_i_j_k(gl_l, gl_i + 1, gl_j + 0, gl_k + 0);
      }
    break;
  }
  case PLUS_K:{
      if(dimension == 2){
      }
      else{
    the_nodes[0] = get_node_number_from_l_i_j_k(gl_l, gl_i + 0, gl_j + 0, gl_k + 1);
    the_nodes[1] = get_node_number_from_l_i_j_k(gl_l, gl_i + 1, gl_j + 0, gl_k + 1);
    the_nodes[2] = get_node_number_from_l_i_j_k(gl_l, gl_i + 1, gl_j + 1, gl_k + 1);
    the_nodes[3] = get_node_number_from_l_i_j_k(gl_l, gl_i + 0, gl_j + 1, gl_k + 1);
      }
    break;		     
  }
  case EDGE0:{
    the_nodes[0] = get_node_number_from_l_i_j_k(gl_l,gl_i+0,gl_j+0,gl_k+0);
    the_nodes[1] = get_node_number_from_l_i_j_k(gl_l,gl_i+1,gl_j+0,gl_k+0);
    break;		     
  }
  case EDGE1:{
    the_nodes[0] = get_node_number_from_l_i_j_k(gl_l,gl_i+1,gl_j+0,gl_k+0);
    the_nodes[1] = get_node_number_from_l_i_j_k(gl_l,gl_i+1,gl_j+1,gl_k+0);
    break;		     
  }
  case EDGE2:{
    the_nodes[0] = get_node_number_from_l_i_j_k(gl_l,gl_i+0,gl_j+1,gl_k+0);
    the_nodes[1] = get_node_number_from_l_i_j_k(gl_l,gl_i+1,gl_j+1,gl_k+0);
    break;		     
  }
  case EDGE3:{
    the_nodes[0] = get_node_number_from_l_i_j_k(gl_l,gl_i+0,gl_j+0,gl_k+0);
    the_nodes[1] = get_node_number_from_l_i_j_k(gl_l,gl_i+0,gl_j+1,gl_k+0);
    break;		     
  }
  case EDGE4:{
    the_nodes[0] = get_node_number_from_l_i_j_k(gl_l,gl_i+0,gl_j+0,gl_k+0);
    if(dimension == 3){
      the_nodes[1] = get_node_number_from_l_i_j_k(gl_l,gl_i+0,gl_j+0,gl_k+1);
    }
    break;		     
  }
  case EDGE5:{
    the_nodes[0] = get_node_number_from_l_i_j_k(gl_l,gl_i+1,gl_j+0,gl_k+0);
    if(dimension == 3){
      the_nodes[1] = get_node_number_from_l_i_j_k(gl_l,gl_i+1,gl_j+0,gl_k+1);
    }
    break;		     
  }
  case EDGE6:{
    the_nodes[0] = get_node_number_from_l_i_j_k(gl_l,gl_i+1,gl_j+1,gl_k+0);
    if(dimension == 3){
      the_nodes[1] = get_node_number_from_l_i_j_k(gl_l,gl_i+1,gl_j+1,gl_k+1);
    }
    break;		     
  }
  case EDGE7:{
    the_nodes[0] = get_node_number_from_l_i_j_k(gl_l,gl_i+0,gl_j+1,gl_k+0);
    if(dimension == 3){
      the_nodes[1] = get_node_number_from_l_i_j_k(gl_l,gl_i+0,gl_j+1,gl_k+1);
    }
    break;		     
  }
  case EDGE8:{
    the_nodes[0] = get_node_number_from_l_i_j_k(gl_l,gl_i+0,gl_j+0,gl_k+1);
    the_nodes[1] = get_node_number_from_l_i_j_k(gl_l,gl_i+1,gl_j+0,gl_k+1);
    break;		     
  }
  case EDGE9:{
    the_nodes[0] = get_node_number_from_l_i_j_k(gl_l,gl_i+1,gl_j+0,gl_k+1);
    the_nodes[1] = get_node_number_from_l_i_j_k(gl_l,gl_i+1,gl_j+1,gl_k+1);
    break;		     
  }
  case EDGE10:{
    the_nodes[0] = get_node_number_from_l_i_j_k(gl_l,gl_i+1,gl_j+1,gl_k+1);
    the_nodes[1] = get_node_number_from_l_i_j_k(gl_l,gl_i+0,gl_j+1,gl_k+1);
    break;		     
  }
  case EDGE11:{
    the_nodes[0] = get_node_number_from_l_i_j_k(gl_l,gl_i+0,gl_j+0,gl_k+1);
    the_nodes[1] = get_node_number_from_l_i_j_k(gl_l,gl_i+0,gl_j+1,gl_k+1);
    break;		     
  }
  case VERTEX0:{
    the_nodes[0] = get_node_number_from_l_i_j_k(gl_l,gl_i+0,gl_j+0,gl_k+0);
    break;		     
  }
  case VERTEX1:{
    the_nodes[0] = get_node_number_from_l_i_j_k(gl_l,gl_i+1,gl_j+0,gl_k+0);
    break;		     
  }
  case VERTEX2:{
    the_nodes[0] = get_node_number_from_l_i_j_k(gl_l,gl_i+1,gl_j+1,gl_k+0);
    break;		     
  }
  case VERTEX3:{
    the_nodes[0] = get_node_number_from_l_i_j_k(gl_l,gl_i+0,gl_j+1,gl_k+0);
    break;		     
  }
  case VERTEX4:{
    the_nodes[0] = get_node_number_from_l_i_j_k(gl_l,gl_i+0,gl_j+0,gl_k+1);
    break;		     
  }
  case VERTEX5:{
    the_nodes[0] = get_node_number_from_l_i_j_k(gl_l,gl_i+1,gl_j+0,gl_k+1);
    break;		     
  }
  case VERTEX6:{
    the_nodes[0] = get_node_number_from_l_i_j_k(gl_l,gl_i+1,gl_j+1,gl_k+1);
    break;		     
  }
  case VERTEX7:{
    the_nodes[0] = get_node_number_from_l_i_j_k(gl_l,gl_i+0,gl_j+1,gl_k+1);
    break;		     
  }
  default:
    
    break;
  }   
}


/****************************************************************************/
void Inline_Mesh_Desc::Calc_Parallel_Info(                                          
					  const std::vector <long long> & element_vector,
					  const std::vector<long long> & global_node_vector,   
					  const std::map <long long, long long> & global_node_map,                          
					  std::list   <long long> & internal_node_list,	
					  std::list   <long long> & border_nodes_list,
					  std::list   <long long> & internal_element_list,
					  std::list   <long long> & border_elements_list,
					  std::list   <long long> & node_proc_id_list,
					  std::list   <long long> & element_proc_id_list,
					  std::vector <long long> & node_neighbor_vector,
					  std::list <long long> * &  boundary_node_list,
					  std::vector <long long> & element_neighbor_vector,
					  std::list <std::pair <long long ,Topo_Loc > > * & boundary_element_list)
  /****************************************************************************/
{

  Tel *el_array = NULL;
  if(element_vector.size()>0){
    el_array = new Tel[element_vector.size()];
    for(unsigned gev = 0; gev < element_vector.size(); gev ++){
      el_array[gev].global_id = element_vector[gev];
      el_array[gev].real_element = true;
    }
  }
  Tel * node_array = NULL;
  if(global_node_vector.size()>0){
    node_array = new Tel[global_node_vector.size()];
    
    for(unsigned gnv = 0;gnv < global_node_vector.size();gnv ++){
      node_array[gnv].global_id = global_node_vector[gnv];
    }
  }

  // walk all local elements
  // peer through their neighboring faces and ask if that neighbor is on my processor
  // if it is do nothing
  // if it is not
  //   I am a border
  //   all nodes on that face are border nodes 
  //   expand node_array and el_array objects with neighbor proc ids and topology directions

  // walk all local elements

  long long nfn = 2;
  long long nfaces = 4;
  long long nen = 1;
  
  if(dimension == 3){
    nfn = 4;
    nfaces = 6;
    nen = 2;
  }




  for(unsigned gev = 0; gev < element_vector.size(); gev ++){
    long long face_nodes_array[4];
    long long my_id = el_array[gev].global_id;
    long long ll,li,lj,lk;
    get_l_i_j_k_from_element_number(my_id,ll,li,lj,lk);
    el_array[gev].visits ++;
    
    for(long long face_count = 0; face_count < nfaces; face_count ++){
      Topo_Loc tl = (Topo_Loc)face_count;
      long long neighbor = get_neighbor(tl,ll,li,lj,lk);
      if(neighbor >= 0){
	long long neighbor_proc_id = Element_Proc(neighbor);
	if(neighbor_proc_id != my_rank){
	  std::pair < long long,Topo_Loc> conn_pair(neighbor_proc_id,tl);
	  el_array[gev].conn_connections.push_back(conn_pair);

	  el_array[gev].visits ++;
	  get_face_nodes(tl,my_id,face_nodes_array);
     
	  for(long long fnc = 0; fnc < nfn; fnc ++){
	    node_array[get_map_entry(global_node_map,face_nodes_array[fnc])].visits ++;
	    node_array[get_map_entry(global_node_map,face_nodes_array[fnc])].proc_neighbors.push_back(neighbor_proc_id);
	  }
	}
      }
    }


    // need to do edges and vertex neighbors for nodes
    long long edge_start = EDGE4;
    long long edge_end = EDGE8;
    if(dimension == 3){
      edge_start = EDGE0;
      edge_end = VERTEX0;
    }
    for(long long edge_count = edge_start; edge_count < edge_end; edge_count ++){
      Topo_Loc tl = (Topo_Loc)edge_count;
      long long neighbor = get_neighbor(tl,ll,li,lj,lk);
      if(neighbor >= 0){
	long long neighbor_proc_id = Element_Proc(neighbor);
	if(neighbor_proc_id != my_rank){
	  get_face_nodes(tl,my_id,face_nodes_array);
	  
	  for(long long fnc = 0; fnc < nen; fnc ++){
	    node_array[get_map_entry(global_node_map,face_nodes_array[fnc])].visits ++;
	    node_array[get_map_entry(global_node_map,face_nodes_array[fnc])].proc_neighbors.push_back(neighbor_proc_id);
	  }
	}
      }
    }
    if(dimension == 3){
      // need to do vertices and vertex neighbors for nodes
      for(long long vertex_count = EDGE11; vertex_count < NUM_TOPO_CONNECTIONS; vertex_count ++){
	Topo_Loc tl = (Topo_Loc)vertex_count;
	long long neighbor = get_neighbor(tl,ll,li,lj,lk);
	if(neighbor >= 0){
	  long long neighbor_proc_id = Element_Proc(neighbor);
	  if(neighbor_proc_id != my_rank){
	    get_face_nodes(tl,my_id,face_nodes_array);
	    for(long long fnc = 0; fnc < 1; fnc ++){
	      node_array[get_map_entry(global_node_map,face_nodes_array[fnc])].visits ++;
	      node_array[get_map_entry(global_node_map,face_nodes_array[fnc])].proc_neighbors.push_back(neighbor_proc_id);
	    }
	  }
	}
      }
    }
  }

  for(unsigned i = 0; i < element_vector.size();i ++){
    if(el_array[i].visits > 1){
      // loop over all conn_connections
      std::list < std::pair < long long , Topo_Loc > > ::iterator conit;
      for(conit  = el_array[i].conn_connections.begin();
          conit != el_array[i].conn_connections.end();
          conit ++){
        element_proc_id_list.push_back((*conit).first);
      }
    }
  }
  // sort and uniq element_proc_id_list
  element_proc_id_list.sort();
  element_proc_id_list.unique();

  if(element_proc_id_list.size()){
    boundary_element_list = new std::list < std::pair <long long ,Topo_Loc > > [element_proc_id_list.size()];
  }

  std::map <long long,long long> element_neighbor_proc_map; //key is proc_id value is ordinal
  std::list <long long> ::iterator listit;
  long long the_count = 0;
  for(listit = element_proc_id_list.begin(); listit != element_proc_id_list.end(); listit++,the_count ++){
    element_neighbor_proc_map[*listit] = the_count;
    element_neighbor_vector.push_back(*listit);
  }

  // now populate the maps

  for(unsigned i = 0; i < element_vector.size();i ++){
    long long the_element = element_vector[i];
    if(el_array[i].visits == 1){
      internal_element_list.push_back(the_element);
    }
    if(el_array[i].visits > 1){
      // loop over all conn_connections
      
      border_elements_list.push_back(the_element);
      std::list < std::pair < long long , Topo_Loc > > ::iterator conit;
      for(conit  = el_array[i].conn_connections.begin();
          conit != el_array[i].conn_connections.end();
          conit ++){
	
        long long index = get_map_entry(element_neighbor_proc_map,(*conit).first);
        boundary_element_list[index].push_back( std::pair <long long ,Topo_Loc >(the_element,(*conit).second));
      }
    }
  }

  border_elements_list.sort();
  border_elements_list.unique();


  for(unsigned gnv = 0;gnv < global_node_vector.size();gnv ++){
    if(node_array[gnv].visits > 0){
      // loop over all conn_connections
      std::list < long long > ::iterator conit;
      for(conit  = node_array[gnv].proc_neighbors.begin();
	  conit != node_array[gnv].proc_neighbors.end();
	  conit ++){
	node_proc_id_list.push_back((*conit));
      }
    }
  }
    
  node_proc_id_list.sort();
  node_proc_id_list.unique();

  std::map <long long,long long> node_neighbor_proc_map; //key is proc_id value is ordinal
  std::list <long long> ::iterator nlistit;
  the_count = 0;
  for(nlistit = node_proc_id_list.begin(); nlistit != node_proc_id_list.end(); nlistit++,the_count ++){
    node_neighbor_proc_map[*nlistit] = the_count;
    node_neighbor_vector.push_back(*nlistit);
  }

  if(node_proc_id_list.size()){
    boundary_node_list = new std::list <long long> [node_proc_id_list.size()];
  }



  //node array needs global_id!!!!
  for(unsigned i = 0;i < global_node_vector.size();i ++){
    if(node_array[i].visits == 0){
      long long the_node = node_array[i].global_id;
      internal_node_list.push_back(the_node);
    }
    else if(node_array[i].visits > 0){
      long long the_node = node_array[i].global_id;
      // loop over all conn_connections
      std::list < long long > ::iterator conit;
      for(conit  = node_array[i].proc_neighbors.begin();
	  conit != node_array[i].proc_neighbors.end();
	  conit ++){
	long long index = get_map_entry(node_neighbor_proc_map,(*conit));
	boundary_node_list[index].push_back(the_node);
	border_nodes_list.push_back(the_node);
      }
    }
  }
  // sort the boundary_node_list
  for(unsigned i = 0; i < node_proc_id_list.size();i++){
    boundary_node_list[i].sort();
    boundary_node_list[i].unique();    
  }
  border_nodes_list.sort();
  border_nodes_list.unique();

  
  
  
  delete [] el_array;
  delete[] node_array;

  //check number node comm_maps
}

//! Reads in/creates the serial component of the unstructured mesh in parallel
/****************************************************************************/
void Inline_Mesh_Desc::Calc_Serial_Component(const std::set <long long> & global_element_ids,
					     const std::vector<long long> & global_node_vector)
/****************************************************************************/
{
  //NODESETS
  // Nodesets are simple because we are numbering nodes across the entire domain
  // sequentially in i,j,k
  // The loop limits are the index bounds that allow traversal of the nodes of interest.
  //DMHMOD

    // SET_INTERSECT
  std::list < long long > tnodes_vector;
  std::list < long long > tglobal_nodes_vector;
  
  std::set < long long > global_nodes_set;
  
  if(nodeset_list.size() > 0){
    nodeset_vectors = new std::vector <long long> [nodeset_list.size()];

    for(unsigned the_nct = 0; the_nct < global_node_vector.size();the_nct++){
      tglobal_nodes_vector.push_back(global_node_vector[the_nct]);
    }
    
    tglobal_nodes_vector.sort();
    tglobal_nodes_vector.unique();
    global_nodes_set.insert(tglobal_nodes_vector.begin(),tglobal_nodes_vector.end());
  }

  std::list < PG_BC_Specification *> ::iterator setit;
  long long nsct = 0;
  for(setit = nodeset_list.begin(); setit != nodeset_list.end();setit++,nsct ++){
    for(unsigned ict = 0;ict < (*setit)->the_locs.size();ict ++){
      LoopLimits ll = (*setit)->the_locs[ict].limits;
      if(dimension == 3){
	for ( long long _nk_ = ll.ks; _nk_ < ll.ke; _nk_ ++){ 
	  for ( long long _nj_ = ll.js; _nj_ < ll.je; _nj_ ++){ 
	    for ( long long _ni_ = ll.is; _ni_ < ll.ie; _ni_ ++){ 
	      long long global_node_id = get_node_number_from_l_i_j_k(trisection_blocks,_ni_,_nj_,_nk_);
	      if(global_nodes_set.find(global_node_id) != global_nodes_set.end())
		nodeset_vectors[nsct].push_back(global_node_id);
	    }
	  }
	}
      }
      else{
	long long _nk_ = 0;
	for ( long long _nj_ = ll.js; _nj_ < ll.je; _nj_ ++){ 
	  for ( long long _ni_ = ll.is; _ni_ < ll.ie; _ni_ ++){ 
	    long long global_node_id = get_node_number_from_l_i_j_k(trisection_blocks,_ni_,_nj_,_nk_);
	    if(global_nodes_set.find(global_node_id) != global_nodes_set.end())
	      nodeset_vectors[nsct].push_back(global_node_id);
	  }
	}
      }
    }
    std::vector < long long > :: iterator vit;
    std::sort  (nodeset_vectors[nsct].begin(),nodeset_vectors[nsct].end());
    vit = std::unique(nodeset_vectors[nsct].begin(),nodeset_vectors[nsct].end());
    nodeset_vectors[nsct].resize(vit-nodeset_vectors[nsct].begin());
  }
  // END OF NODESETS

  //SIDESETS
  // Sidesets are defined by the element id and the id of the face to which the set
  // applies.

  // For the purposes of calculation the sideset definition process starts by 
  // considering all elements to be numbered sequentially in i then j then k
  // across all blocks. This is not actually the case since exodus requires elements
  // to be numbered sequentially within a block.

  // If we consider the elements to be numbered sequentially across blocks in i,j, and k
  // then we can use a loop limits type call to get the i,j,and k index limits of the
  // elements in the sideset. We can then loop over these limits to calculate the 
  // index of the element if it were numbered across blocks. It is then left to 
  // calculate the index of this element if it were numbered sequentially within blocks.
  // This is done by calculating the global i,j,and k indices of the element within the 
  // entire domain, calculating the indices of the block within which the element resides,
  // calculating the ordinal number of the block in which the element resides, and the 
  // local i,j,and k indices of the element within the block it resides.

  //  These values are combined to calculate the index of the element that corresponds
  // to numbering the elements sequentially within blocks. 
  nsct = 0;
  if(sideset_list.size() > 0){
    sideset_vectors = new std::vector < std::pair <long long ,Topo_Loc > > [sideset_list.size()];  
    sideset_global_count = new long long [sideset_list.size()];  
    for(unsigned ict = 0; ict < sideset_list.size();ict++)sideset_global_count[ict] = 0.;
  }
  for(setit = sideset_list.begin(); setit != sideset_list.end();setit++,nsct ++){
    for(unsigned ict = 0;ict < (*setit)->the_locs.size();ict ++){
      
      LoopLimits ll = (*setit)->the_locs[ict].limits;
      Topo_Loc the_location =  (*setit)->the_locs[ict].location;
      //Sidesets allowed only on faces of block, not on edges or corners
            
      if(dimension == 3){
	for ( long long _nk_ = ll.ks; _nk_ < ll.ke; _nk_ ++){ 
	  for ( long long _nj_ = ll.js; _nj_ < ll.je; _nj_ ++){ 
	    for ( long long _ni_ = ll.is; _ni_ < ll.ie; _ni_ ++) {
	      long long elnumber = get_element_number_from_l_i_j_k(trisection_blocks,_ni_,_nj_,_nk_);
	      if(!isElementSuppressed(elnumber))sideset_global_count[nsct]++;
	      if(global_element_ids.find(elnumber)!=global_element_ids.end()){
		std::pair <long long ,Topo_Loc > el_loc_pair(elnumber,the_location);
		sideset_vectors[nsct].push_back(el_loc_pair);
	      }
	    }
	  }
	}
      }
      else{
	long long _nk_ = 0;
	for ( long long _nj_ = ll.js; _nj_ < ll.je; _nj_ ++){ 
	  for ( long long _ni_ = ll.is; _ni_ < ll.ie; _ni_ ++) {
	    long long elnumber = get_element_number_from_l_i_j_k(trisection_blocks,_ni_,_nj_,_nk_);
	    if(!isElementSuppressed(elnumber))sideset_global_count[nsct]++;
	    if(global_element_ids.find(elnumber)!=global_element_ids.end()){	  
	      std::pair <long long ,Topo_Loc > el_loc_pair(elnumber,the_location);
	      sideset_vectors[nsct].push_back(el_loc_pair);
	    }
	  }
	}
      }
    }
  }
  //END of SIDESETS
}

/****************************************************************************/
void Inline_Mesh_Desc::getGlobal_Element_Block_Totals(long long * totals_array)
/****************************************************************************/
{
  for(long long bct = 0; bct < numBlocks();bct ++ ){
    long long kind = bct/blockKstride();
    long long jind = (bct - kind * blockKstride())/inline_b[0];
    long long iind = bct - jind * inline_b[0] - kind * blockKstride();
    totals_array[bct] = a_inline_n[0][iind]*a_inline_n[1][jind]*a_inline_n[2][kind];
  } 
}


//! Queries which processor an element lies on.
//! Calls the recursive Partition::Element_Proc function.
/****************************************************************************/
void Inline_Mesh_Desc::Customize_Coords(double * coords, long long num_nodes,long long dim)
/****************************************************************************/
{
  if(!Geometry_Transform_Function)return;
  Geometry_Transform_Function->Operate(coords,num_nodes,dim);
}

/****************************************************************************/
void Inline_Mesh_Desc::Offset_Coords(double * coords, long long num_nodes,long long dim)
/****************************************************************************/
{
  for(long long ict = 0; ict < num_nodes; ict ++){
    for(long long idim = 0; idim < dim; idim ++){
      coords[idim*num_nodes + ict] = coords[idim*num_nodes + ict] + inline_offset[idim];
    }
  }
}

/****************************************************************************/
long long Inline_Mesh_Desc::Calc_Coord_Vectors()
/****************************************************************************/
{

  for(long long axis = 0; axis < dimension; axis ++){
    IJKcoors[axis] = new double[nel_tot[axis]+1];

    long long nct = 0;
    for(long long i = 0; i < inline_b[axis]; i ++){
      double sum = 0.;
      for(long long j = 0; j < a_inline_n[axis][i]; j ++){
	if((first_size[axis][i] > 0.) && (last_size[axis][i] > 0.)){
	  IJKcoors[axis][nct] = c_block_dist[axis][i]+sum;
	  sum += first_size[axis][i];
	  if(interval[axis][i]-1) sum += (double)j*(last_size[axis][i]-first_size[axis][i])/((double)interval[axis][i]-1);
	  IJKcoors[axis][nct+1] = c_block_dist[axis][i+1];
	}
	else{
	  IJKcoors[axis][nct] = c_block_dist[axis][i]+j*block_dist[axis][i]/(double)a_inline_n[axis][i];
	  IJKcoors[axis][nct+1] = c_block_dist[axis][i]+(j+1)*block_dist[axis][i]/(double)a_inline_n[axis][i];
	}
	nct ++;
      }
    }
    if(Element_Density_Functions[axis]){
      Element_Density_Functions[axis]->Integrate(inline_gmin[axis],inline_gmax[axis], error_stream);
      if(!error_stream.str().empty()){return 1;}
      double delta = inline_gmax[axis]-inline_gmin[axis];
      for(long long ict = 0; ict < (nel_tot[axis]+1); ict ++){
	double factor = (IJKcoors[axis][ict]-inline_gmin[axis])/delta;
	double interpolant =  Element_Density_Functions[axis]->Interpolate(factor, error_stream);if(!error_stream.str().empty())return 1;
	double new_coord = inline_gmin[axis]+interpolant*delta;
	IJKcoors[axis][ict] = new_coord;
      }
    }
  }

  return 0;
}



}//end namespace PAMGEN_NEVADA
