/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef Stk_Rebalance_Use_Cases_UseCase_2_hpp
#define Stk_Rebalance_Use_Cases_UseCase_2_hpp

#include <stk_mesh/fixtures/HexFixture.hpp>

namespace stk {
namespace rebalance {
namespace use_cases {

  bool test_heavy_nodes( stk::ParallelMachine comm );

} //namespace use_cases
} //namespace rebalance
} //namespace stk

#endif // Stk_Rebalance_Use_Cases_UseCase_2_hpp

/// \page stk_rebalance_use_case_2
///  \ingroup stk_rebalance_use_case_module
///
/// \section stk_rebalance_use_case_2_description Use Case 2: Node, Edge, Face and Element weights on subset
///
/// This use case demonstrates element weights comprised of contributions
/// from nodes, edges, faces and elements over a subset of the mesh.
/// A 3x3x3 cube of hex8 elements is constructed on proc 0
/// and weights are assigned to a subset of mesh entities.  The mesh and weights
/// are assigned as follows:
///
///       Global node and element numbering
///      <pre>
///
///                +-------+-------+-------+               +                   +                 +
///               /       /       /       /|              /|                  /|         
///              /       /       /       / |             / |                 / |        
///             +-------+-------+-------+  |            +  |                +  |    
///            /       /       /       /|  +           /|  +               /   +   
///           /       /       /       / | /|          / | /|              /    |  
///          +-------+-------+-------+  |/ |         +  |/ |             +     | 
///         /       /       /       /|  +  |        /|  +  |            /      |
///        /       /       /       / | /|  +       / | /|  +           /       +
///       +-------+-------+-------+  |/ | /|      +  |/ | /|          +        |        +
///       |       |       |       |  +  |/ |      |  +  |/ |          |        |    
///       |  e1   |  e2   |  e3   | /|  +  |      | /|  +  |          |        |        
///       |       |       |       |/ | /|  +      |/ | /|  +          |        +                 +
///       +-------+-------+-------+  |/ | /       +  |/ | /           +       /         
///       |       |       |       |  +  |/        |  +  |/            |      /      
///       |  e1   |  e2   |  e3   | /|  +         | /|  +             |     +      
///       |       |       |       |/ | /          |/ | /              |    /      
///       +-------+-------+-------+  |/           +  |/               +   /      
///       |       |       |       |  +            |  +                |  +      
///       |  e1   |  e2   |  e3   | /             | /                 | /      
///       |       |       |       |/              |/                  |/      
///       +-------+-------+-------+               +                   +                 +
///     x = 0
///      </pre>
///
///   <pre>
///   Weight_elems = 1.0                 Z  Y      Local node numbering
///   Weight_faces = 10.0                | /      
///   Weight_edges = 100.0               |/            8       7       
///   Weight_nodes = 1000.0              *--X          +-------+              
///                                                   /       /|       
///                                                 5/      6/ |       
///                                                 +-------+  |   
///                                                 |       |  +3   
///                                                 |  e1   | /
///                                                 |       |/
///                                                 +-------+
///                                                 1       2
///   </pre>
///
/// where all 27 elements are assigned weights along with the 9 faces, 12 edges and 4 nodes
/// on the plane at x = 0.
/// \dontinclude UseCase_Rebal_2.cpp
/// \skip bulk.modification_begin
/// \until bulk.modification_end
///
/// The use case passes if the amount of imbalance following a rebalance is
/// below 1.45 for 3 procs and below 1.1 for 2 or 4 procs.
/// 
/// See \ref UseCase_Rebal_2.cpp for the complete source listing.
