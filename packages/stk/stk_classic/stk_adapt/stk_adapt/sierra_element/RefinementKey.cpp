/*--------------------------------------------------------------------*/
/*    Copyright 2003 - 2009 Sandia Corporation.                       */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

// Copyright 2003 Sandia Corporation, Albuquerque, NM.

#include <stk_adapt/sierra_element/stk_percept_code_types.hpp>

#include <stk_adapt/sierra_element/RefinementKey.hpp>

namespace stk_classic { namespace adapt {
namespace Elem {

/*
  Derek Gaston: I've pretty much disabled this class so it always returns 0x0.
  This is because I've removed heterogeneous refinement from the framework.
  Hopefully this class will get completely removed one day...
*/

RefinementKey::RefinementKey()
{
  //Default constructor will assign 0x0 to value_
  //This will correspond to normal homogeneous refinement patterns
  value_=0x0;
}
RefinementKey::RefinementKey(UInt valueArgue)
{
  //constructor where value is known
//  value_=valueArgue;
  value_=0x0;

}

RefinementKey::~RefinementKey()
{
//Destructor
//Do Nothing
}

UInt Elem::RefinementKey::value()
{
  return value_;
}

void Elem::RefinementKey::assign_value(UInt valueArgue)
{
  //Assigns value_ based on a given id #
  //Usefull when changing value to default during unrefinement
  //as well as when assigning the unrefinement template id #

//  value_ = valueArgue;

}

void Elem::RefinementKey::assign_value(std::vector<UInt> & edge_order)
{
  //Assigns value_ based on a given edge order cut vector
  //Used durring heterogeneous marker resolution

/*  value_ = 0x0;

  for (UInt position = 0; position < edge_order.size(); ++position) {
    value_ |= (edge_order[position]+1)<<(position*4);
  }
*/
}

std::vector<UInt> Elem::RefinementKey::ordered_cut_edges( UInt numEdges ) const
{
  //Given the number of edges of a mesh object populates a vector
  //of edges to be refined in the correct order based on value_

  std::vector<UInt> edge_order;
  edge_order.clear();

  for (UInt iedge = 0; iedge<numEdges; ++iedge){
    UInt edge = value_>>(iedge*4) & 0xF;
    if(edge != 0x0) edge_order.push_back(edge-1);
  }
  return edge_order;

}


bool Elem::RefinementKey::full_refinement( UInt numEdges )
{
  //if number of edges on refinement key equal numEdges than full refinement
  std::vector<UInt> edgeOrder = ordered_cut_edges(numEdges);
  UInt edgeOrderSize = edgeOrder.size();
  return edgeOrderSize == numEdges;
}

} // namespace Elem
} // namespace adapt
} // namespace stk_classic
