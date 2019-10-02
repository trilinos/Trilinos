// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

/**
 * @file
 * @specifications for heterogeneous key id to give refinement topologies
 * @author Warren B. Tauber
 *
 * @par Description:
 *
 * In order to support removal of hanging nodes elements can be refined
 * according to topologies that do not refine every edge and or face.  To
 * account for the different possibilities a key id is given to each mesh
 * object.  This key id will contain the edges to be cut information .
 * Given the key id number member functions of this class will be able to
 * give a vector of ordered edges to be cut.  Also the reverse will be
 * available as well.  These functions will be different for each base
 * topology class.
 *
 * value_ will be such that for 2d solid triangular elements
 * a value of 0x021 would be edge 0 cut first then edge 1
 * a value 0f 0x321 would be edge 0 then edge 1 then edge 2
 * value of digits in key_id corresponds to (edge# +1 ) in left to right
 * ordering
 */

/*--------------------------------------------------------------------*/

#ifndef adapt_sierra_element_Elem_RefinementKey_hpp
#define adapt_sierra_element_Elem_RefinementKey_hpp

#include <vector>
#include <adapt/sierra_element/percept_code_types.hpp>

namespace percept {
  namespace Elem {

    class RefinementKey {
    private:
      /** Unsigned integer that will store the correct template **/
      UInt value_;

    public:

      /** Default Constructer that assigns the value_ to 0x0 **/
      RefinementKey();

      /** Constructer where value is known **/
      explicit RefinementKey(UInt valueArg);

      /** Destructor **/
      ~RefinementKey();

      /** Returns the value_ **/
      UInt value();

      /** Assigns value_**/
      void assign_value(UInt valueArg) ;

      /** function that assigns the value_ based on edges to be cut **/
      void assign_value( std::vector<UInt> & edge_order ) ;

      /** function that returns ordered edges to be cut **/
      /** for this method need the objTopologyes number of edges **/
      std::vector<UInt> ordered_cut_edges( UInt numEdges ) const ;

      /** Function returns boolean for full refinement **/
      bool full_refinement( UInt numEdges );
    };

  } // namespace Elem
} // namespace percept

#endif // adapt_sierra_element_Elem_RefinementKey_hpp
