// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_CellToolsDocumentation.hpp
    \brief  Header file with additional documentation for the Intrepid2::CellTools class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_CELLTOOLS_DOCUMENTATION_HPP__
#define __INTREPID2_CELLTOOLS_DOCUMENTATION_HPP__

/***************************************************************************************************
 **                                                                                               **
 **                           D O C U M E N T A T I O N   P A G E S                               **
 **                                                                                               **
 **************************************************************************************************/

/**  
\page cell_tools_page                 Cell tools

<b>Table of contents </b>

-   \ref cell_topology_sec
-   \ref cell_topology_ref_cells
  - \ref sec_cell_topology_ref_map
  - \ref sec_cell_topology_ref_map_DF
  - \ref sec_cell_topology_subcell_map
  - \ref sec_cell_topology_subcell_wset

\section cell_topology_sec            Cell topologies

The range of admissible cell shapes in Intrepid is restricted to <var>d</var>-dimensional 
<strong>polytopes</strong>, <var>d=1,2,3</var>. A polytope is defined by a set of its vertices 
\f$ \{ {\bf v}_0,\ldots {\bf v}_V\} \f$ and a <strong>base topology</strong>  <var>BT</var> that 
defines how these verices are connected into <var>k</var>-dimensional, <var>k < d</var> facets 
(k-subcells) of that polytope.  

The base topology of any polytope can be extended by augmenting the set of its vertices by 
an additional set of points \f$\{ {\bf p}_0,\ldots {\bf p}_P\} \f$. The <strong>extended topology</strong>
<var>ET</var> is defined by specifying the connectivity of the set
\f$ \{ {\bf v}_0,\ldots {\bf v}_V\}\cup \{ {\bf p}_0,\ldots {\bf p}_P\} \f$
relative to the subcells specified by its base topology <var>BT</var>.
The vertices and the extra points are collectively referred to as <strong>nodes</strong>. Thus,
a polytope with <strong>extended topology</strong> <var>ET</var> is defined by a set of nodes 
\f$\{{\bf q}_0,\ldots,{\bf q}_N\}\f$, where \f$N = V + P\f$, and a connectivity rule for these nodes.

Intrepid requires any cell to have a valid base topology. The nodes of the cell should always be 
ordered by listing its vertices <strong>first</strong>, i.e.,
\f[ 
  \{{\bf q}_0,\ldots,{\bf q}_N\} = \{ {\bf v}_0,\ldots {\bf v}_V,{\bf p}_0,\ldots {\bf p}_P\} 
  \f]
To manage cell topologies Intrepid uses the Shards package http://trilinos.org/packages/shards .
Shards provides definitions for a standard set of base and extended cell topologies plus tools to
construct custom, user defined cell topologies, such as arbitrary polyhedral cells. For further
details see Shards documentation. 



\section cell_topology_ref_cells      Reference cells

For some cell topologies there exist simple, e.g., polynomial, mappings that allow to obtain any 
cell having that topology as an image of a single "standard" cell. We refer to such standard cells
as <strong>reference</strong> cells. 

Just like in the general case, a reference cell with a base topology <var>BT</var> is defined by a 
set of vertices, and a reference cell with extended topology <var>ET</var> is defined by a set of 
nodes that include the original vertices and some additional points. 

The actual vertex and node coordinates for the reference cells can be chosen arbitrarily; however, 
once selected they should not be changed because in many cases, e.g., in finite element reconstructions,
all calculations are done on a reference cell and then transformed to physical cells by an appropriate
pullback (see Section \ref sec_pullbacks).

In Intrepid base and extended reference cell topologies are defined using the following selections
of vertex and node coordinates:

\verbatim
|=======================|==============================================================================|
| Topology family       |    reference cell vertices/additional nodes defining extended topology       |
|=======================|==============================================================================|
| Line<2>               |                                                                              |
| Beam<2>               | {(-1,0,0),(1,0,0)}                                                           |
| ShellLine<2>          |                                                                              |
|-----------------------|------------------------------------------------------------------------------|
| Line<3>               |                                                                              |
| Beam<3>               | {(0,0,0)}                                                                    |
| ShellLine<3>          |                                                                              |
|=======================|==============================================================================|
| Triangle<3>           | {(0,0,0),(1,0,0),(0,1,0)}                                                    |
| ShellTriangle<3>      |                                                                              |
|-----------------------|------------------------------------------------------------------------------|
| Triangle<4>           | {(1/3,1/3,0)}                                                                |
|.......................|..............................................................................|
| Triangle<6>           | {(1/2,0,0),(1/2,1/2,0),(0,1/2,0)}                                            |
| ShellTriangle<6>      |                                                                              |
|=======================|==============================================================================|
| Quadrilateral<4>      | {(-1,-1,0),(1,-1,0), (1,1,0),(-1,1,0)}                                       |
| ShellQuadrilateral<4> |                                                                              |
|-----------------------|------------------------------------------------------------------------------|
| Quadrilateral<8>      | {(0,-1,0),(1,0,0),(0,1,0),(-1,0,0)}                                          |
| ShellQuadrilateral<8> |                                                                              |
|.......................|..............................................................................|
| Quadrilateral<9>      | {(0,-1,0),(1,0,0),(0,1,0),(-1,0,0),(0,0,0)}                                  |
| ShellQuadrilateral<9> |                                                                              |
|=======================|==============================================================================|
| Tetrahedron<4>        | {(0,0,0),(1,0,0),(0,1,0),(0,0,1)}                                            |
|-----------------------|------------------------------------------------------------------------------|
| Tetrahedron<8>        | {(1/2,0,0),(1/2,1/2,0),(0,1/2,0),(1/3,1/3,1/3)}                              |
| Tetrahedron<10>       | {(1/2,0,0),(1/2,1/2,0),(0,1/2,0),(0,0,1/2),(1/2,0,1/2),(0,1/2,1/2)}          |
|=======================|==============================================================================|
| Pyramid<5>            | {(-1,-1,0),(1,-1,0),(1,1,0),(-1,1,0),(0,0,1)}                                |
|-----------------------|------------------------------------------------------------------------------|
| Pyramid<13>           | {(0,-1,0),(1,0,0),(0,1,0),(-1,0,0), 1/2((-1,-1,1),(1,-1,1),(1,1,1),(-1,1,1))}|
| Pyramid<14>           | all of the above and (0,0,0)                                                 |
|=======================|==============================================================================|
| Wedge<6>              | {(0,0,-1),(1,0,-1),(0,1,-1),(0,0,1),(1,0,1),(0,1,1)}                         |
|-----------------------|------------------------------------------------------------------------------|
| Wedge<15>             | {(1/2,0,-1),(1/2,1/2,-1),(0,1/2,-1), (0,0,0),(1,0,0),(0,1,0),                |  
|                       |  (1/2,0, 1),(1/2,1/2, 1),(0,1/2, 1)                                          |
|.......................|..............................................................................|
| Wedge<18>             | All of the above plus {(1/2,0,0),(1/2,1/2,0),(0,1/2,0)}                      |
|=======================|==============================================================================|
| Hexahedron<8>         | {(-1,-1,-1),(1,-1,-1),(1,1,-1),(-1,1,-1),(-1,-1,1),(1,-1,1),(1,1,1),(-1,1,1)}|
|-----------------------|------------------------------------------------------------------------------|
| Hexahedron<20>        | {(0,-1,-1),(1,0,-1),(0,1,-1),(-1,0,-1), (0,-1,0),(1,0,0),(0,1,0),(-1,0,0),   |
|                       |  (0,-1, 1),(1,0, 1),(0,1, 1),(-1,0, 1) }                                     |
|.......................|..............................................................................|
| Hexahedron<27>        | All of the above plus center point and face midpoints:                       |
|                       | {(0,0,0), (0,0,-1),(0,0,1), (-1,0,0),(1,0,0), (0,-1,0),(0,1,0)}              |
|=======================|==============================================================================|
\endverbatim

Finite element reconstruction methods based on pullbacks (see Section \ref sec_pullbacks) are 
restricted to the above cell topologies.  


  
\subsection sec_cell_topology_ref_map      Reference-to-physical cell mapping

The mapping that takes a given reference cell to a physical cell with the same topology is defined 
using a nodal Lagrangian basis corresponding to the nodes of the reference cell. In other words, the
mapping is constructed using basis functions that are dual to the nodes of the reference cell.
Implementation details are as follows.

Assume that \f$ \hat{\mathcal{C}} \f$ is a reference cell with topology <var>T</var> and nodes
\f$\{\hat{{\bf q}}_0,\ldots,\hat{{\bf q}}_{N}\}\f$, and that \f$ \{\hat{\phi}_i\}_{i=0}^{N} \f$ is
the Lagrangian basis dual to these nodes, i.e., \f$ \hat{\phi}_i( \hat{{\bf q}}_j)=\delta_{ij} \f$.
A physical cell \f$\mathcal{C}\f$ with the same topology <var>T</var> as \f$\hat{\mathcal{C}}\f$ is  
then defined as the image of \f$ \hat{\mathcal{C}} \f$ under the mapping
\f[               
  F_\mathcal{C}(\hat{\bf x}) = \sum_{m=0}^{N}  {\bf q}_m(\mathcal{C})  \hat{\phi}_m(\hat{\bf x})          
\f]
where \f$\{{\bf q}_0(\mathcal{C}),\ldots,{\bf q}_N(\mathcal{C})\}\f$ is the set of <strong>physical nodes</strong>
that defines \f$\mathcal{C}\f$. The number of physical nodes is required to match the number of reference  
nodes in the specified cell topology <var>T</var>. The <var>i</var>-th coordinate function of the 
reference-to-physical mapping is given by
\f[        
  \big(F_{\mathcal{C}}(\hat{\bf x})\big)_i = 
              \sum_{m=0}^{N}  \big({\bf q}_m(\mathcal{C})\big)_i\hat{\phi}_m(\hat{\bf x})           
\f]
where \f$ \big({\bf q}_m(\mathcal{C})\big)_i \f$ is the <var>i</var>-th spatial coordinate of the <var>m</var>-th node.

For simplicity, unless there's a chance for confusion, the cell symbol will be ommitted from the
designations of physical points and reference-to-physical maps, i.e., we shall simply write
\f$F(\hat{\bf x})\ \mbox{and}\ {\bf q}_m\f$.

\par Summary
  
\li      \f$F(\hat{\bf x})\f$: implemented in Intrepid2::CellTools::mapToPhysicalFrame 
\li      \f$F^{-1}({\bf x})\f$: implemented in Intrepid2::CellTools::mapToReferenceFrame 
  
\warning Intrepid2::CellTools does not check for non-degeneracy of the physical cell obtained from a
given set of physical nodes. As a result, <var>F</var> is not guaranteed to be a diffeomorphism,
i.e., it may not have a continuously differentiable inverse. In this case some 
Intrepid2::CellTools methods, such as Intrepid2::CellTools::setJacobianInv, 
and Intrepid2::CellTools::mapToReferenceFrame will fail.
    

  
\subsection sec_cell_topology_ref_map_DF   Jacobian of the reference-to-physical cell mapping

Intrepid follows the convention that the rows of the Jacobian are the transposed gradients of the
coordinate functions of the mapping, i.e.,
\f[                    
      DF_{ij}(\hat{{\bf x}}) = \frac{\partial F_i(\hat{{\bf x}})}{\partial\hat{{\bf x}}_j}                    
\f]
In light of the definition of <var>F</var> in Section \ref sec_cell_topology_ref_map, it follows that
\f[  
      DF_{ij}(\hat{{\bf x}}) = \sum_{m=0}^{N}
                ({\bf q}_m)_i\frac{\partial\hat{\phi}_m(\hat{{\bf x}})}{\partial\hat{{\bf x}}_j} \,.       
\f]

\par Summary
  
\li     \f$DF_{ij}(\hat{{\bf x}}) \f$: implemented in Intrepid2::CellTools::setJacobian
\li     \f$DF_{ij}^{-1}(\hat{{\bf x}}) \f$: implemented in Intrepid2::CellTools::setJacobianInv
\li     \f$\mbox{det} DF_{ij}(\hat{{\bf x}}) \f$: implemented in  Intrepid2::CellTools::setJacobianDet                      


  
\subsection sec_cell_topology_subcell_map  Parametrization of physical 1- and 2-subcells
  
Parametrization of a given physical k-subcell \f$\mathcal{S}_i\f$, k=1,2, is a map from a 
k-dimensional parametrization domain \e R to that subcell:
\f[
    \Phi : R \mapsto \mathcal{S} \,.
\f]
Parametrization domains play role similar to that of reference cells in the sense that they allow 
computation of line and surface integrals on 1- and 2-subcells (edges and faces) to be reduced to   
computation of integrals on \e R . 

Parametrization maps are supported for 1- and 2-subcells (edges and faces) that belong to physical 
cells with reference cells. The reason is that these maps are defined by the composition of  the
parametrization maps for reference edges and faces with the mapping \e F defined in \ref sec_cell_topology_ref_map.
As a result, parametrization of a given physical k-subcell requires selection of a 
<strong>parent cell</strong> that contains the subcell. 

\remark  
Because a given k-subcell may belong to more than one physical cell, its parent cell is not unique. For a single
k-subcell the choice of a parent cell is not important, however, when dealing with subcell worksets
parent cells must all have the same topology (see \ref sec_cell_topology_subcell_wset for details about
subcell worksets).
  
  
Implementation of subcell parametrization is as follows. Assume that \f$\mathcal{S}_i\f$ is a k-subcell   
with parent cell \f$\mathcal{C}\f$; \f$\mathcal{\hat{C}}\f$ is the associated reference cell and 
\e i is the local ordinal of the subcell relative to the reference cell. To this physical subcell
corresponds a reference subcell \f$\hat{\mathcal{S}}_i\f$ having the same local ordinal. Parametrization of
the reference k-subcell is a map from the k-dimensional parametrization domain \e R to that subcell:
\f[
    \hat{\Phi} : R \mapsto \hat{\mathcal{S}}_i \,.
\f]
Parametrization of \f$\mathcal{S}_i\f$ is then defined as
\f[
    \Phi = F\circ \hat{\Phi} \,,
\f]
where \e F is the reference-to-physical mapping between the parent cell and its reference cell. 
  

A 1-subcell (edge) always has \c Line<N> topology and so, the parametrization domain for edges is the standard 1-cube:
\f[
      R = [-1,1]\,.
\f]
On the other hand, faces of reference cells can have \c Triangle<N> and/or \c Quadrilateral<N> topologies. 
Thus, the parametrization domain for a 2-subcell depends on its topology and is either the standard 
2-simplex or the standard 2-cube:
\f[
      R = \left\{\begin{array}{rl} 
          \{(0,0),(1,0),(0,1)\} & \mbox{if the face is Triangle} \\[1ex]
            [-1,1]\times [-1,1] & \mbox{if the face is Quadrilateral}
          \end{array}\right.
\f]
\par Summary

-    \f$ \Phi : R \mapsto {\mathcal{S}}_i \f$ requires two steps: 
  -# Intrepid2::CellTools::mapToReferenceSubcell to apply \f$\hat{\Phi}: R \mapsto \hat{\mathcal{S}}_i\f$;
  -# Intrepid2::CellTools::mapToPhysicalFrame to apply \e F



\subsection sec_cell_topology_subcell_wset Subcell worksets

A subcell workset comprises of 1- or 2-subcells and associated parent cells that satisfy the 
following conditions

\li     all subcells have the same cell topology;
\li     all parent cells have the same cell topology;
\li     The parent cell topology has a reference cell;
\li     relative to that reference cell, all subcells in the workset have the same local ordinal

Therefore, a subcell workset is defined by 

-#      collecting a set of 1- or 2-subcells having the same topology
-#      selecting a parent cell for every subcell in such a way that 
  -#      all parent cells have the same cell topology
  -#      all subcells in the workset have the same local ordinal relative to the parent cell topology

  Obviously, a subcell can have multiple parent cells. For example, in a mesh consisiting of Triangle<3> cells, every
  edge is shared by 2 triangle cells. To define an edge workset we can use either one of the two traingles
  sharing the cell.

  Suppose now that the mesh comprises of Triangle<3> and Quadrilateral<4> cells and we want to define an
  edge workset. Let's say the first few edges in our workset happen to be shared by 2 triangles and so 
  we select one of them as the parent cell. Now suppose the next edge is shared by a traingle and a 
  quadrilateral. Because all parent cells in the workset must have the same cell topology we cannot use
  the quadrilateral as a parent cell and so we choose the triangle. Finally suppose that one of the candidate 
  edges for our workset is shared by 2 quadrilaterals. Because of the requirement that all parent cells
  have the same topology, we will have to reject this edge because it does not posses a potential parent cell
  with the same topology as the rest of the edges in our workset.

  A subcell workset is denoted by \f$ \{\mathcal{S}_{c,i}\}_{c=0}^N \f$, where

  \li <var>c</var> is parent cell ordinal;
  \li <var>i</var> is the local subcell ordinal (relative to the topology of the parent cell) shared
  by all subcells in the workset. 

  */
#endif
