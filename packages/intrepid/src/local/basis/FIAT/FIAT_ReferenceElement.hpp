// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Robert Kirby (robert.c.kirby@ttu.edu) 
// ************************************************************************
// @HEADER

/** \file FIAT_ReferenceElement.hpp
    \brief  Contains reference elements used by FIAT
    \author Created by R. Kirby.
*/

#ifndef FIAT_REFERENCE_ELEMENT_HPP
#define FIAT_REFERENCE_ELEMENT_HPP

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Array.hpp"

using Teuchos::SerialDenseMatrix;
using Teuchos::RefCountPtr;
using Teuchos::Array;
using Teuchos::rcp;

namespace FIAT
{

  /** \enum FIAT::ElementShape
      \brief Indicates the element shapes FIAT understands
  */
  enum ElementShape { 
    FIAT_LINE = 0 
  };


  /** \class FIAT::ReferenceElement
      \brief base class for all reference elements used in FIAT
  */
  template<class Scalar>
  class ReferenceElement {
  public:
    /** \brief constructor */
    ReferenceElement( const int spatialDim , 
		      const ElementShape elementShape ,
		      RefCountPtr<SerialDenseMatrix<int,Scalar> > verts,
		      RefCountPtr<Array<Array<Array<int> > > > topology,
		      RefCountPtr<SerialDenseMatrix<int,Scalar> > normals,
		      RefCountPtr<SerialDenseMatrix<int,Scalar> > scaledNormals,
		      RefCountPtr<SerialDenseMatrix<int,Scalar> > tangents,
		      RefCountPtr<SerialDenseMatrix<int,Scalar> > scaledTangents,
		      RefCountPtr<SerialDenseMatrix<int,Scalar> > faceTangents ,
		      RefCountPtr<Array<Array<Scalar> > > measure ) :
      spatialDim_( spatialDim ) , elementShape_( elementShape ) ,
      verts_( verts ) , topology_( topology ) , normals_( normals ) , scaledNormals_( scaledNormals ) ,
      tangents_( tangents ) , scaledTangents_( scaledTangents ) , faceTangents_( faceTangents ) ,
      measure_( measure ) { }

    /** \brief destructor */
    virtual ~ReferenceElement() {}

    /** \brief Returns the ElementShape of the element */
    inline ElementShape getShape() const 
    { 
      return elementShape_; 
    }

    /** \brief Returns the spatial dimension in which the 
	reference element lives */
    inline int getSpatialDimension() const 
    { 
      return spatialDim_; 
    }

    /** \brief returns a (refcountptr to) the array of vertices, stored as Scalar. 
	the first index is for the points, the second for the coordinates of the point
     */
    inline RefCountPtr<SerialDenseMatrix<int,Scalar> > getVertices( ) const 
    {
      return verts_;
    }

    /** \brief Returns the connectivity of the reference element.
	On output, (*top)[i][j] is the array of vertex ids comprising
	facet j of topological dimension i.
    */
    inline RefCountPtr<Array<Array<Array<int> > > > getTopology( ) const
    {
      return topology_;
    }

    /** \brief returns the number of facets of the given spatial
     dimension, raises an exception if an illegal dimension is specified */
    int getNumFacets( const int topdim ) const;

    /** \brief Returns an array containing the coordinates of vertices
	of the given facet of the given topological dimension, raising an
	exception if these are out of bounds.
     */
    RefCountPtr<SerialDenseMatrix<int,Scalar> > 
    getVerticesOfSubcomplex( const int topdim ,
			     const int facet ) const;

    /** \brief returns a RefCountPtr to an array of of the unit normals to facets
	of codimension 1.
	The first index of the array runs over the facets and the second index
	runs over the components of the normal vector.
     */
    inline RefCountPtr<SerialDenseMatrix<int,Scalar> > getNormals( ) const
    {
      return normals_;
    }

    /** \brief returns a RefCountPtr to an array of of the unit normals to facets
	of codimension 1, each scaled by the size of that facet.
	The first index of the array runs over the facets and the second index
	runs over the components of the normal vector.
     */
    inline RefCountPtr<SerialDenseMatrix<int,Scalar> > getScaledNormals( ) const
    {
      return scaledNormals_;
    }

    /** \brief returns a RefCountPtr to an array of of the unit tangents to facets
	of dimension 1 (edges)
	The first index of the array runs over the edges and the second index
	runs over the components of the tangent vector.
     */
    inline RefCountPtr<SerialDenseMatrix<int,Scalar> > getTangents( ) const
    {
      return tangents_;
    }

    /** \brief returns a RefCountPtr to an array of of the unit tangents to facets
	of dimension 1 (edges), each scaled by the size of that facet.
	The first index of the array runs over the facets and the second index
	runs over the components of the tangent vector.
     */
    inline RefCountPtr<SerialDenseMatrix<int,Scalar> > getScaledTangents( ) const
    {
      return scaledTangents_;
    }
    
    /** \brief returns a RefCountPtr to an array of of the scaled tangents to facets
	of codimension 1.
	The first index of the array runs over the facets and the second index
	runs over the components of the tangent vector.
	Raises an exception if the spatial dimension is not 3.
    */
    RefCountPtr<SerialDenseMatrix<int,Scalar> > getFaceTangents( ) const;
    
    
    /** \brief	Constructs the coordinates of a lattice of points on the reference element.
	n is the order of the lattice, so that if n == 1, the lattice includes the vertices.
	If n == 2, the lattice includes the vertices and edge midpoints, and so on.
	If the interior flag is other than zero, then that many points from the boundary in
	are skipped.  For example, if n = 3 and interior = 1 on a triangle, only the triangle's
	midpoint is included.
	The return array stores the points id in the first index and
	the spatial coordinate in the second.  
    */
    virtual RefCountPtr<SerialDenseMatrix<int,Scalar> > 
    makeLattice( const int order ,
		 const int interior = 0 ) const = 0;

    /** \brief Constructs the lattice of points on the interior of the given facet
	of the given topological dimension.  
	The resize method of pts may be invoked.
	An exception is raised on illegal facet number.
     */
    // this is not virtual, but is implemented in terms of makeLattice
    // and an AffineMap
    RefCountPtr<SerialDenseMatrix<int,Scalar> >
    makePoints( const int topdim , const int facet , const int n ) const;

    /** \brief Returns the size of the reference element (length in 1d, area in 2d, etc)
     */
    inline Scalar getMeasure( ) const
    {
      return (*measure_)[ spatialDim_ ][ 0 ];
    }

    /** \brief Returns the size of the given facet, raising an exception on illegal input */
    Scalar getMeasure( const int topdim , const int facet ) const;

  private:
    const int spatialDim_;
    const ElementShape elementShape_;
    RefCountPtr<SerialDenseMatrix<int,Scalar> > verts_;
    RefCountPtr<Array<Array<Array<int> > > > topology_;
    RefCountPtr<SerialDenseMatrix<int,Scalar> > normals_;
    RefCountPtr<SerialDenseMatrix<int,Scalar> > scaledNormals_;
    RefCountPtr<SerialDenseMatrix<int,Scalar> > tangents_;
    RefCountPtr<SerialDenseMatrix<int,Scalar> > scaledTangents_;
    RefCountPtr<SerialDenseMatrix<int,Scalar> > faceTangents_;
    RefCountPtr<Array<Array<Scalar> > > measure_;
  };

}

#include "FIAT_ReferenceElementDef.hpp"
 

#endif
