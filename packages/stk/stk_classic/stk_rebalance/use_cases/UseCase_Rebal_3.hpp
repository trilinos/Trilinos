/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef Stk_Rebalance_Use_Cases_UseCase_3_hpp
#define Stk_Rebalance_Use_Cases_UseCase_3_hpp

#include <stk_mesh/fixtures/HexFixture.hpp>

namespace stk_classic {
namespace rebalance {
namespace use_cases {

  bool test_contact_surfaces( stk_classic::ParallelMachine comm );

} //namespace use_cases
} //namespace rebalance
} //namespace stk_classic

#endif // Stk_Rebalance_Use_Cases_UseCase_3_hpp

/// \page stk_rebalance_use_case_3
///  \ingroup stk_rebalance_use_case_module
///
/// \section stk_rebalance_use_case_3 Use Case 3: Periodic Boundary via Constraint Relations
///
/// This use case sets up a 2D mesh of quad4 elements and then establishes a constraint
/// realtion between the top and bottom of the mesh as would be needed to enforce
/// periodic boundary conditions.
///
/// The following quad4 mesh is manually constructed on proc 0:
///
///       Global node and element numbering
///      <pre>
///     
///        21      22      23      24      25
///        +-------+-------+-------+-------+   y = top
///        |       |       |       |       |
///        |  e13  |  e14  |  e15  |  e16  |
///        |       |       |       |5      |
///     16 +-------+-------+-------+-------+ 20
///        |       |       |       |       |      
///        |  e9   |  e10  |  e11  |  e12  |      
///        |       |       |       |       |      
///     11 +-------+-------+-------+-------+ 15
///        |       |       |       |       |
///        |  e5   |  e6   |  e7   |  e8   |
///        |       |5      |       |       |  
///      6 +-------+-------+-------+-------+ 10
///        |       |       |       |       |  
///        |  e1   |  e2   |  e3   |  e4   |    
///        |       |       |       |       | 
///        +-------+-------+-------+-------+   y = bottom
///        1       2       3       4       5
///      </pre>
///
///  Local node numbering:
///
///      <pre>
///        3       4
///        +-------+     Y
///        |       |     |
///        |  e1   |     |
///        |       |     *--> X
///        +-------+
///        1       2
///      </pre>
///
/// and the two sets of nodes at y=bottom and y=top are
/// related through constraint relations as follows:
/// \dontinclude UseCase_Rebal_3.cpp
/// \skip Assign constraint relations between nodes
/// \until end snippet
/// 
/// The use case passes if the load imbalance of the new partition
/// is below the nominal value of 1.5.
///
/// See \ref UseCase_Rebal_3.cpp for the complete source listing.
///
