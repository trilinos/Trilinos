/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef Stk_Rebalance_Use_Cases_UseCase_1_hpp
#define Stk_Rebalance_Use_Cases_UseCase_1_hpp

#include <stk_mesh/fixtures/HexFixture.hpp>

namespace stk {
namespace rebalance {
namespace use_cases {

  bool test_unequal_weights( stk::ParallelMachine comm );

} //namespace use_cases
} //namespace rebalance
} //namespace stk

#endif // Stk_Rebalance_Use_Cases_UseCase_1_hpp

/// \page stk_rebalance_use_case_1
///  \ingroup stk_rebalance_use_case_module
///
/// \section stk_rebalance_use_case_1_description Use Case 1: Unequal element weights
///
/// This use case demonstrates unequal element weights. A
/// single-element-thick column of hex elements is constructed in the
/// x-direction. Weights are assigned to elements such that a perfect
/// rebalance is possible with Proc_Rank+1 elements placed on each processor.
/// This is achieved when running on #procs < 4, but Zoltan does not
/// produce the optimal new distribution for 4 or more procs.  
/// 
///
/// The following hex8 mesh is used in this use case:
///
///       Global node and element numbering
///      <pre>
///           3       7      11      15      19                 
///           +-------+-------+-------+-------+                +-------+
///          /       /       /       /       /|               /       /|
///        4/      8/     12/     16/     20/ |              /       / |    Z  Y
///        +-------+-------+-------+-------+  |     ......  +-------+  |    | / 
///        |       |       |       |       |  +18           |       |  /    |/  
///        |  e1   |  e2   |  e3   |  e4   | /              |  eN   | /     *--X       
///        |       |       |       |       |/               |       |/         
///        +-------+-------+-------+-------+                +-------+
///        1       5      9       13      17                  
///      </pre>
///       where N = #elements = (#procs)(#procs+1)/2
///      
///       Local node numbering
///      <pre>
///           8       7
///           +-------+
///          /       /|
///        5/      6/ |
///        +-------+  |
///        |       |  +3
///        |  e1   | /
///        |       |/
///        +-------+
///        1       2
///      </pre>
///
/// The mesh is constructed on proc 0 using the HexFixture class
/// and unequal element weights are assigned as follows:
/// \dontinclude UseCase_Rebal_1.cpp
/// \skip Assign weights
/// \until end assign weights
/// 
/// See \ref UseCase_Rebal_1.cpp for the complete source listing.
/// 
/// 
/// 
/// 
