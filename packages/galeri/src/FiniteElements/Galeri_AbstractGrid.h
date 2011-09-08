// @HEADER
// ************************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions about Galeri? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
//
// ************************************************************************
// @HEADER

#ifndef GALERI_ABSTRACTGRID_H
#define GALERI_ABSTRACTGRID_H

/* 
 * \file Galeri_AbstractGrid.h
 */

class Epetra_Map;
class Epetra_DistObject;
class Epetra_Comm;

namespace Galeri {
namespace FiniteElements {

/*!
 * \class AbstractGrid
 *
 * \brief Abstract interface to access finite element grids.
 *
AbstractGrid is a pure virtual function, that specifies the interface
methods that a grid class must implement. Following an approach similar to
getrow() for matrices, there is no grid format; instead, the user 
must implement a set of getelement() and similar methods. Therefore, it is
possible to define grids that are built on-the-fly, a feature that is
particularly convenient in testing phase.

This format is based on the following assumptions:
- All grid elements are of the same type (for examples, all triangles).
  It is not possible to have mixed grids (i.e., with some triangles and
  some quadrilateral).
- All grid elements are 3D. For 2D problems, the user must specify a
  z-coordinate, for example 0.0.
- Elements, vertices and faces are numbered locally. The local-to-global
  mapping is required for vertices only, and must be defined using 
  Epetra_Map's.

Two Epetra_Map's must be created:
- The VertexMap() locally contains all the vertices that belong to the
  local finite elements. A global vertex can be replicated over more 
  than one process. 
- The RowMap() globally contains all the global vertices. A global vertex
  is assigned to exactly one process.

We require two maps for the following reason:
In parallel runs, vertices corresponding to internal boundary
faces can are replicated over processors. This makes the contruction
of the finite element matrix easier, since it is possible to
work on local quantities only. However, such a distribution
is not compatible with AztecOO and ML (and several other Trilinos
packages), since these libraries require each row to belong to 
exactly one processor. Methods ExportToVertexMap() and ExportToRowMap()
export objects from one map to the other. Typically, ExportToRowMap() is
used in the assembly phase, whese local stiffness matrix and right-hand
side are built using VertexMap(), then exported to RowMap().
ExportToVertexMap(), instead, is used after the solution of the
linear system is computed, to update the values of the solution on
the local vertices, so that norms can be computed, and the solution
visualized. These methods are trivial in the serial case.
 *
 * \author Marzio Sala, SNL 9214
 *
 * \date Last updated on 12-Apr-05.
 */

class AbstractGrid {

public:

  // @{ Constructors and destructors
  //
  //! Destructor.
  virtual ~AbstractGrid() {}

  // @}
  // @{ Query methods
  
  //! Returns the number of dimensions of the grid.
  virtual int NumDimensions() const = 0;

  //! Returns the number of vertices contained in each element.
  virtual int NumVerticesPerElement() const = 0;

  //! Returns the number of faces contained in each element.
  virtual int NumFacesPerElement() const = 0;

  //! Returns the number of vertices contained in each face.
  virtual int NumVerticesPerFace() const = 0;

  //! Returns a string containing the element type.
  /*! 
   *  Returns a string containing the type of element. This
   *  string is used in the quadrature class.
   *  Currently supported options are:
   *  - "GALERI_TRIANGLE"
   *  - "GALERI_QUAD"
   *  - "GALERI_HEX"
   *  - "GALERI_TET"
   */
  virtual string ElementType() const = 0;

  //! Returns the number of neighboring elements.
  virtual int NumNeighborsPerElement() const = 0;

  //! Returns the number of finite elements on the calling process.
  virtual int NumMyElements() const = 0;

  //! Returns the global number of finite elements.
  virtual int NumGlobalElements() const = 0;

  //! Returns the number of vertices on the calling process.
  virtual int NumMyVertices() const = 0;

  //! Returns the global number of vertices.
  virtual int NumGlobalVertices() const = 0;

  //! Returns the number of boundary faces on the calling process.
  virtual int NumMyBoundaryFaces() const = 0;

  //! Returns the global number of boundary faces.
  virtual int NumGlobalBoundaryFaces() const = 0;

  //! Returns the volume of all local elements.
  virtual double MyVolume() const = 0;

  //! Returns the global volume of the grid.
  virtual double GlobalVolume() const = 0;

  //! Returns the coordinates of local vertex \c LocalVertex in vector \c coord.
  /*!
   * \param LocalVertex - (In) Local ID of the vertex for whic
   *                      coordinates are required. Must be 
   *                      contained in the interval [0, NumMyVertices())
   * 
   * \param coord - (Out) double array of size 3. In output, contains
   *                the x-, y- and z-coordinate of the specified vertex.
   *
   * \note Parameter \c coord must be allocated of size 3 for both
   *       2D and 3D problems.
   */
  virtual void VertexCoord(const int LocalVertex, double* coord) const = 0;

  //! Returns the coordinates of specified local vertices.
  /*!
   *  \param Length - (In) Length of array \c IDs.
   *
   *  \param IDs - (In) Contains the list of vertices of which coordinates
   *                    are required.
   *
   *  \param x - (Out) double array of size \c Length. In output, contains
   *                   the x-coordinates of the specified vertices.
   *
   *  \param y - (Out) double array of size \c Length. In output, contains
   *                   the y-coordinates of the specified vertices.
   *
   *  \param z - (Out) double array of size \c Length. In output, contains
   *                   the z-coordinates of the specified vertices.
   *
   *  \note The \c z array must be allocated for both 2D and 3D problems.
   */
  virtual void VertexCoord(const int Length, const int* IDs, double* x, 
                           double* y, double* z) const = 0;

  //! Returns the local vertex IDs of the specified local finite element.
  /*!
   *  \param LocalElement - (In) ID of the required local element.
   *
   *  \param elements - (Out) array of length NumElementVertices(), in
   *                          output will contain the local ID of the
   *                          vertices of the specified element.
   */
  virtual void ElementVertices(const int LocalElement, int* elements) const = 0;

  //! Returns the local IDs of neighboring elements.
  virtual void ElementNeighbors(const int LocalElement, int* elements) const = 0;

  //! Returns the local vertex IDs of vertices contained in the specified boundary face.
  virtual void FaceVertices(const int LocalFace, int& tag, int* IDs) const = 0;

  //! Returns the patch ID of the specified face.
  /*! Returns an integer ID that identifies the given boundary face
   *  as belonging to a given part of the domain. It can be used by the
   *  user to specify the value and the type of the boundary condition.
   */
  virtual int FacePatch(const int LocalFace) const = 0;

  //! Returns the volume of the specified local finite element.
  virtual double ElementMinLength(const int LocalElement) const = 0;

  //! Returns the volume of the specified local finite element.
  virtual double ElementMaxLength(const int LocalElement) const = 0;

  //! Returns the volume of the specified local finite element.
  /*! Returns the area (in 2D) or the volume (in 3D) of the specified
   *  local element
   */
  virtual double ElementVolume(const int LocalElement) const = 0;

  //! Returns the area of the specified local face.
  /*! Returns the length (in 2D) or the area (in 3D) of the
   *  specified boundary face
   */
  virtual double FaceArea(const int LocalFace) const = 0;

  // @}
  // @{ Maps and import/export

  //! Returns a reference to the map representing the vertex distribution.
  virtual const Epetra_Map& VertexMap() const = 0;

  //! Returns a reference to the map representing the distribution of rows.
  virtual const Epetra_Map& RowMap() const = 0;

  //! Exports distributed object from RowMap() to VertexMap().
  virtual void ExportToVertexMap(const Epetra_DistObject& RowObject,
                                 Epetra_DistObject& VertexObject) const = 0;

  //! Exports distributed object from VertexMap() to RowMap().
  virtual void ExportToRowMap(const Epetra_DistObject& RowObject,
                              Epetra_DistObject& VertexObject) const = 0;

  //! Returns a reference to the communicator object.
  virtual const Epetra_Comm& Comm() const = 0;

  // @}
};

} // namespace FiniteElements
} // namespace Galeri
#endif
