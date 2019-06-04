.. _fem_assembly_index:

Finite Element Assembly
#######################

Overview
========

Finite element assembly is the process by which a global system of equations over a finite element mesh is generated.  The global system of equations typically represent the linearized discrete forms of one, or more, conservation equations expressing some physics being solved by the finite element model.  Typically, the global equations are formed by looping over the elements in a mesh, calculating the element contribution, and merging the contribution in to global system.  In many production finite element programs, assembly is the most complicated phase of the solution procedure due to the complexity of supporting a large library of element formulations.  In the following examples, the assembly process is simplified significantly by restricting attention to the two-dimensional mesh of four-node elements having a single degree of freedom (DOF) at each node, shown in Figure :numref:`FEMA_Mesh_Example`.  Simplifying the element equation construction allows focusing, instead, on the creation of the Tpetra objects used to represent the global system of equations and merging element equations in to it.

.. _FEMA_Mesh_Example:
.. figure:: /Images/Tpetra-FEM-Assembly-Example-Mesh.png
    :width: 400px
    :align: center
    :alt: Simple 3x3 mesh with nodes.

    Two-dimensional default mesh used in each assembly example.

Phases of Finite Element Assembly
=================================

Construction of the global system of equations typically follows two distinct phases:

- construction of a |tpetra_crsgraph|_ representing the mesh's nodal connectivity; and
- construction of a |tpetra_crsmatrix|_ from the previously constructed graph containing the numeric values of the global system of equations.

Graph Construction
------------------

A |tpetra_crsgraph|_ represents the structure of the matrix associated with the global system of equations.  The graph is used to compute the communication patterns for later algorithms on the data.  Nodes of the graph are related to, but not the same, as nodes of the finite element mesh.  The nodes of the graph are the rows in its (the graph's) matrix representation and represent a single equation, or degree of freedom, in the global system.  The nodes of the finite element mesh, on the other hand, represent the mesh's vertices and can be associated with no or multiple degrees of freedom.  For the simple mesh used in the following examples, the terms node, degree of freedom, or vertex can be used interchangeably to refer to either the graph or mesh entity.  However, for more complex systems containing elements of different types, multifreedom constraints, multiple degrees of freedom per node, etc., one must be careful to distinguish between node and mesh entities.

Element connectivity is represented by edges of the graph.  An edge ``(i,j)`` represents connectivity between degrees of freedom ``i`` and ``j``.  Thus, an element's connectivity is represented by the clique of vertices defining it. For example, vertices 0, 1, 4, and 5 are the clique of vertices defining element 0, as shown in Figure :numref:`FEMA_E0_Clique`.

.. _FEMA_E0_Clique:
.. figure:: /Images/Tpetra-FEM-Assembly-Example-E0-Clique.png
    :width: 400px
    :align: center
    :alt: Element 0 Clique

    Element 0 Clique

Each node contains edges to its adjacent nodes.  Nodes are adjacent if their respective mesh nodes belong to the same element.  For example, node 0 contains a self-edge as well as edges to nodes 1, 4, and 5.  Nodes belonging to multiple elements will contain edges to all of its neighboring nodes, as shown in :numref:`FEMA_Node_Adjacency`.

.. _FEMA_Node_Adjacency:
.. figure:: /Images/Tpetra-FEM-Assembly-Example-Node-Connectivity.png
    :width: 400px
    :align: center
    :alt: Node Connectivity

    Node Connectivity

    Examples showing how the node connectivity of the ultimate graph.  Nodes are 'adjacent' if they belong to the same element.

The assembly process results in the CRS graph shown whose matrix is shown in :numref:`FEMA_CrsGraph_Adjacency_Matrix`:

.. _FEMA_CrsGraph_Adjacency_Matrix:
.. figure:: /Images/Tpetra-FEM-Assembly-Example-CrsGraph.png
    :width: 400px
    :align: center
    :alt: CrsGraph Adjacency

    Matrix representation of the CRS graph for the example problems' finite
    element mesh.

Matrix Construction
-------------------

The |tpetra_crsmatrix|_ is constructed from the mesh's graph representation and contains the global system's numeric values.  Once constructed, the values in the matrix may be changed, but the structure (graph) of the matrix should not change (i.e., zeros cannot become nonzeros) without constructing a new graph.

Assembly Types
==============

Examples of three assembly "types" - Types 1-3, respectively - are presented:

* :ref:`Type 1 <fem_assembly_type1>`
* :ref:`Type 2 <fem_assembly_type2>`
* :ref:`Type 3 <fem_assembly_type3>`

.. toctree::
   :maxdepth: 1
   :hidden:

   Type-1
   Type-2
   Type-3

Source Code
===========


The source code for each example can be found in the Trilinos repository: |tpetra_example_fem_assembly|_


.. include:: /links.txt


