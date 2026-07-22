# Shards : Shared Discretization Tools

The purpose of Shards is to provide a suite of common tools for numerical and topological data that facilitate interoperability between typical software modules used to solve Partial Differential Equations (PDEs) by finite element, finite volume and finite difference methods. Shards comprises of two categories of tools: methods to manage and access information about cell topologies used in mesh-based methods for PDEs, and methods to work with multi-dimensional arrays used to store numerical data in corresponding computer codes. The basic cell topology functionality of Shards includes methods to query adjacencies of subcells, find subcell permutation with respect to a global cell and create user-defined custom cell topologies. Multi-dimensional array part of the package provides specialized compile-time dimension tags, multi-index access methods, rank and dimension queries.

Shards design is based on a domain model for cell topology data motivated by algebraic topology. In this approach the mesh is viewed as a chain complex consisting of 0,1,2 and 3-dimensional cells representing the nodes, edges, faces and elements in the mesh. Cell topologies are explicitly defined by a composition of subcells of dimension less or equal to that of the parent cell.

## Overview

<strong>Cells and cell topology</strong>

In Shards <strong>cell</strong> refers to a _d_-dimensional polytope, _d=1,2,3_. A polytope is defined by a set of its vertices_{V0,V1,…,Vv}_ and a <strong>base topology _BT_</strong> that defines how these verices are connected into _k_-dimensional, _k_ < d facets (_k_-subcells) of that polytope. The base topology of any polytope can be extended by augmenting the set of its vertices by an additional set of points _{P0,P1,…,Pp}_. The <strong>extended topology</strong> _<strong>ET</strong>_ is defined by specifying the connectivity of the set _{V0,V1,…,Vv}_ + _{P0,P1,…,Pp}_ relative to the subcells specified by its base topology _BT_. The vertices and the extra points are collectively referred to as <strong>nodes</strong>. Thus, a polytope with extended topology _ET_ is defined by a set of nodes _{N0,N1,…,Nn}_ where _n = v + p_, and a connectivity rule for these nodes.

Shards provides definitions for a standard set of base and extended cell topologies plus tools to construct custom, user defined cell topologies, such as arbitrary polyhedral cells. The nodes in all Shards topologies are ordered by listing the cell vertices first. For example, the nodes of Triangle<6>are ordered as _{N0,N1,N2,N3,N4,N5}_ where _{N0,N1,N2}_ are the three vertices and _{N3,N4,N5}_ are the three edge midpoints. Every cell with an extended topology also has a base topology. The base topology of Triangle<6> is Triangle<3>, defined by its 3 vertices _{V0,V1,V2}._

<strong>Remark.</strong> Shards provides only cell topologies, i.e., the connectivity rules for a cell and its lower-dimensional subcells. Shards does not specify what are the coordinates of the nodes forming a cell.

A list of all Shards cell topologies follows. Bold face numerals indicate the vertices of each cell topology.

<table width="700" border="1"><caption>0- and 1-dimensional cell topologies</caption>

<tbody>

<tr>

<th scope="col">Name</th>

<th scope="col">nodes</th>

<th scope="col">Description</th>

</tr>

<tr>

<th scope="row"><span style="font-family: Courier New,Courier,monospace;">Node<></span></th>

<td>{<strong>0</strong>}</td>

<td>A single node topology</td>

</tr>

<tr>

<th scope="row"><span style="font-family: Courier New,Courier,monospace;">Line</span><span style="font-family: Courier New,Courier,monospace;"><2></span></th>

<td>{<strong>0,1</strong>}</td>

<td>Base line topology, equivalent to <span style="font-family: Courier New,Courier,monospace;">Line<></span></td>

</tr>

<tr>

<th scope="row"><span style="font-family: Courier New,Courier,monospace;">Line</span><span style="font-family: Courier New,Courier,monospace;"><3></span></th>

<td>{<strong>0,1</strong>,2}</td>

<td>Extended line topology with 1 edge node.</td>

</tr>

</tbody>

</table>

<table width="700" border="1"><caption>2-dimensional cell topologies</caption>

<tbody>

<tr>

<th scope="col">Name</th>

<th scope="col">nodes</th>

<th scope="col">Description</th>

</tr>

<tr>

<th scope="row"><span style="font-family: Courier New,Courier,monospace;">Triangle</span><span style="font-family: Courier New,Courier,monospace;"><3></span></th>

<td>{<strong>0,1,2</strong>}</td>

<td>Base triangle (2-simplex), same as <span style="font-family: Courier New,Courier,monospace;">Triangle</span><span style="font-family: Courier New,Courier,monospace;"><></span></td>

</tr>

<tr>

<th scope="row"><span style="font-family: Courier New,Courier,monospace;">Triangle<4></span></th>

<td>{<strong>0,1,2</strong>,3}</td>

<td>Extended triangle with 1 interior node</td>

</tr>

<tr>

<th scope="row"><span style="font-family: Courier New,Courier,monospace;">Triangle<6></span></th>

<td>{<strong>0,1,2</strong>,3,4,5}</td>

<td>Extended triangle with 3 edge nodes</td>

</tr>

<tr>

<th scope="row"><span style="font-family: Courier New,Courier,monospace;">Quadrilateral<4></span></th>

<td>{<strong>0,1,2,3</strong>}</td>

<td>Base quadrilateral (2-cube), same as <span style="font-family: Courier New,Courier,monospace;">Quadrilateral<></span></td>

</tr>

<tr>

<th scope="row"><span style="font-family: Courier New,Courier,monospace;">Quadrilateral<8></span></th>

<td>{<strong>0,1,2,3</strong>,4,5,6,7}</td>

<td>Extended quadrilateral with 4 edge nodes</td>

</tr>

<tr>

<th scope="row"><span style="font-family: Courier New,Courier,monospace;">Quadrilateral<9></span></th>

<td>{<strong>0,1,2,3</strong>,4,5,6,7,8}</td>

<td>Extended quadrilateral with 4 edge and 1 interior nodes</td>

</tr>

<tr>

<th scope="row"><span style="font-family: Courier New,Courier,monospace;">ShellLine<2></span></th>

<td>{<strong>0,1</strong>}</td>

<td>Base shell line, same as <span style="font-family: Courier New,Courier,monospace;">ShellLine<></span></td>

</tr>

<tr>

<th scope="row"><span style="font-family: Courier New,Courier,monospace;">ShellLine<3></span></th>

<td>{<strong>0,1</strong>,2}</td>

<td>Extended shell line with 1 edge node</td>

</tr>

<tr>

<th scope="row"><span style="font-family: Courier New,Courier,monospace;">Beam</span><span style="font-family: Courier New,Courier,monospace;"><2></span></th>

<td>{<strong>0,1</strong>}</td>

<td>Base beam, same as <span style="font-family: Courier New,Courier,monospace;">Beam<></span></td>

</tr>

<tr>

<th scope="row"><span style="font-family: Courier New,Courier,monospace;">Beam<3></span></th>

<td>{<strong>0,1</strong>,2}</td>

<td>Extended beam with 1 edge node</td>

</tr>

<tr>

<th scope="row"><span style="font-family: Courier New,Courier,monospace;">Pentagon</span><span style="font-family: Courier New,Courier,monospace;"><5></span></th>

<td>{<strong>0,1,2,3,4</strong>}</td>

<td>Base pentagon, same as <span style="font-family: Courier New,Courier,monospace;">Pentagon<></span></td>

</tr>

<tr>

<th scope="row"><span style="font-family: Courier New,Courier,monospace;">Hexagon</span><span style="font-family: Courier New,Courier,monospace;"><6></span></th>

<td>{<strong>0,1,2,3,4,5</strong>}</td>

<td>Base hexagon, same as <span style="font-family: Courier New,Courier,monospace;">Hexagon<></span></td>

</tr>

</tbody>

</table>

<table width="700" border="1"><caption>3-dimensional cell topologies</caption>

<tbody>

<tr>

<th scope="col">Name</th>

<th scope="col">nodes</th>

<th scope="col">Description</th>

</tr>

<tr>

<th scope="row"><span style="font-family: Courier New,Courier,monospace;">Tetrahedron<4></span></th>

<td>{<strong>0,1,2,3</strong>}</td>

<td>Base tetrahedron (3-simplex), same as <span style="font-family: Courier New,Courier,monospace;">Tetrahedron<></span></td>

</tr>

<tr>

<th scope="row"><span style="font-family: Courier New,Courier,monospace;">Tetrahedron<8></span></th>

<td>{<strong>0,1,2,3</strong>,4,5,6,7}</td>

<td>Extended tetrahedron with 4 face nodes</td>

</tr>

<tr>

<th scope="row"><span style="font-family: Courier New,Courier,monospace;">Tetrahedron<10></span></th>

<td>{<strong>0,1,2,3</strong>,4,5,6,7,8,9}</td>

<td>Extended tetrahedron with 6 edge nodes</td>

</tr>

<tr>

<th scope="row"><span style="font-family: Courier New,Courier,monospace;">Hexahedron<8></span></th>

<td>{<strong>0,1,2,3,4,5,6,7</strong>}</td>

<td>Base hexahedron, same as <span style="font-family: Courier New,Courier,monospace;">Hexahedron</span><span style="font-family: Courier New,Courier,monospace;"><></span></td>

</tr>

<tr>

<th scope="row"><span style="font-family: Courier New,Courier,monospace;">Hexahedron<20></span></th>

<td>{<strong>0,1,2,3,4,5,6,7</strong>,…,20}</td>

<td>Extended hexahedron with 12 edge nodes</td>

</tr>

<tr>

<th scope="row"><span style="font-family: Courier New,Courier,monospace;">Hexahedron<27></span></th>

<td>{<strong>0,1,2,3,4,5,6,7</strong>,…,27}</td>

<td>Extended hexahedron with 12 edge, 6 face & 1 interior nodes</td>

</tr>

<tr>

<th scope="row"><span style="font-family: Courier New,Courier,monospace;">Pyramid<5></span></th>

<td>{<strong>0,1,2,3,4</strong>}</td>

<td>Base pyramid, same as <span style="font-family: Courier New,Courier,monospace;">Pyramid</span><span style="font-family: Courier New,Courier,monospace;"><></span></td>

</tr>

<tr>

<th scope="row"><span style="font-family: Courier New,Courier,monospace;">Pyramid<13></span></th>

<td>{<strong>0,1,2,3,4</strong>,5,…,12}</td>

<td>Extended pyramid with 8 edge nodes</td>

</tr>

<tr>

<th scope="row"><span style="font-family: Courier New,Courier,monospace;">Pyramid<14></span></th>

<td>{<strong>0,1,2,3,4</strong>,5,…,13}</td>

<td>Extended pyramid with 8 edge and 1 quad face nodes<span style="font-family: Courier New,Courier,monospace;"></span></td>

</tr>

<tr>

<th scope="row"><span style="font-family: Courier New,Courier,monospace;">Wedge<6></span></th>

<td>{<strong>0,1,2,3,4,5</strong>}</td>

<td>Base wedge, same as <span style="font-family: Courier New,Courier,monospace;">Wedge</span><span style="font-family: Courier New,Courier,monospace;"><></span></td>

</tr>

<tr>

<th scope="row"><span style="font-family: Courier New,Courier,monospace;">Wedge<15></span></th>

<td>{<strong>0,1,2,3,4,5</strong>,6,…,14}</td>

<td>Extended wedge with 9 edge nodes<span style="font-family: Courier New,Courier,monospace;"></span></td>

</tr>

<tr>

<th scope="row"><span style="font-family: Courier New,Courier,monospace;">Wedge<18></span></th>

<td>{<strong>0,1,2,3,4,5</strong>,6,…,17}</td>

<td>Extended wedge with 9 edge and 3 quad face nodes<span style="font-family: Courier New,Courier,monospace;"></span></td>

</tr>

</tbody>

</table>

## Documentation

Shards is part of the [Trilinos Project](https://trilinos.github.io), and additional information (e.g., examples, tutorials, and source code documentation) is available through [Shards's Doxygen webpages](https://trilinos.github.io/docs/shards/index.html).

## Questions?

Contact the lead developers:

- **Shards team**:    GitHub handle: @trilinos/Shards
- **Mauro Perego**:   GitHub handle: [mperego](https://github.com/mperego) or mperego@sandia.gov


## Copyright and License

For general copyright and license information, refer to the Trilinos [License and Copyright](https://trilinos.github.io/about.html#license-and-copyright) page.

For Shards-specific copyright and license details, refer to the [shards/COPYRIGHT](COPYRIGHT) and [shards/LICENSE](LICENSE) files located in the `shards` directory. Additional copyright information may also be found in the headers of individual source files.

For developers, general guidance on documenting copyrights and licenses can be found in the Trilinos [Guidance on Copyrights and Licenses](https://github.com/trilinos/Trilinos/wiki/Guidance-on-Copyrights-and-Licenses) document.
