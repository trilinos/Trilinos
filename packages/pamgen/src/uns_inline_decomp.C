// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "uns_inline_decomp.h"
#include <assert.h>
#include <math.h>
#include "uns_inline_decomp.h"
#include "inline_mesh_desc.h"
#include <ostream>
#include <sstream>

namespace PAMGEN_NEVADA {

//! Recursive function to find the id of the processor on which 
//! the element in question resides.
/****************************************************************************/
  long long Partition::Element_Proc(long long Ginds[])
  /****************************************************************************/
  {
    if(high == NULL)return proc_id;
    if(Ginds[split_direction]>=split_value){
      return high->Element_Proc(Ginds);
    }
    else{
      return low->Element_Proc(Ginds);
    }
  }
  
  /****************************************************************************/
  void Partition::print()
  /****************************************************************************/
  {
  }
}
long long PAMGEN_NEVADA::Partition::partition_count = 0;

namespace PAMGEN_NEVADA {
  
  //! Bisects the a partition and adds the results to the_list.
  /****************************************************************************/
  void   Partition::Processor_Partition(std::vector < Partition * > & the_list,
					long long inc_nels[])
  /****************************************************************************/
  {
    // the following should never happen
    assert(high == NULL);
    assert(low  == NULL);
    long long l_lows[3];
    long long h_lows[3];
    long long h_highs[3];
    long long l_highs[3];
    long long dels[3];
    //assign same highs and lows initially and calculate ranges
    for(long long i = 0 ; i < 3; i ++){
      dels[i] = highs[i]-lows[i];
      l_lows[i] = lows[i];
      h_lows[i] = lows[i];
      h_highs[i] = highs[i];
      l_highs[i] = highs[i];
    }
    
    
    // find largest range
    split_direction = 0;
    if((dels[0] / inc_nels[0]) > 1 )split_direction = 0;
    if((dels[1] / inc_nels[1]) > 1 )split_direction = 1;
    if((dels[2] / inc_nels[2]) > 1 )split_direction = 2;
    
    
    assert(dels[split_direction]>1);
    //The following prevents remainder building up on last cut.
    long long local_split_size = inc_nels[split_direction];
    
    if((dels[split_direction]%remaining_cuts[split_direction]) > 1)local_split_size ++;
    
    
    split_value = lows[split_direction] + local_split_size;
    h_lows[split_direction] = split_value;
    l_highs[split_direction] = split_value;
    remaining_cuts[split_direction] =  remaining_cuts[split_direction]-1;
    
    high = new Partition(lows[3],h_lows[0],
			 h_lows[1],h_lows[2],
			 highs[3],h_highs[0],
			 h_highs[1],h_highs[2],
			 inline_decomposition_type,remaining_cuts);
    low  = new Partition(lows[3],l_lows[0],
			 l_lows[1],l_lows[2],
			 highs[3],l_highs[0],
			 l_highs[1],l_highs[2],
			 inline_decomposition_type,remaining_cuts);
    
    the_list.push_back(low);
    the_list.push_back(high);
  }
  
}// end namespace
