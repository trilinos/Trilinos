// @HEADER
// ************************************************************************
//
//            Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
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
//
// Questions about Galeri? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
//
// ************************************************************************
// @HEADER

/*! \file Galeri_grid_Loadable.h
 *
 * \brief A flexible grid data structure container for distributed problems.
 *
 * \author Marzio Sala
 *
 * \date Last modified on Aug-06
 */

#ifndef GALERI_LOADABLE_GRID_H
#define GALERI_LOADABLE_GRID_H

#include "../src/Galeri_ConfigDefs.h"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Hashtable.hpp"

#include "Galeri_core_Object.h"
#include "Galeri_core_Workspace.h"
#include "Galeri_grid_Element.h"
#include "Galeri_grid_Point.h"
#include "Galeri_grid_Segment.h"
#include "Galeri_grid_Triangle.h"
#include "Galeri_grid_Quad.h"
#include "Galeri_grid_Tet.h"
#include "Galeri_grid_Hex.h"

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Export.h"

#include "EpetraExt_DistArray.h"

using namespace Teuchos;

namespace Galeri {

namespace grid {

/*! \class Loadable
 *
 * \brief A flexible grid data structure container for distributed problems.
 *
 * Galeri::grid::Loadable is a loadable container for grid data structures. It
 * allows the allocation, setting and getting of all grid objects.
 *
 * A Galeri::grid::Loadable object is defined by the following entities:
 * - an Epetra_Comm;
 * - the number of global elements;
 * - the number of local elements (typically different on each processor).
 * - the element type;
 * - the number of global vertices;
 * - the number of local vertices.
 *  
 * Optionally, it is possible to assigned to each grid element and grid
 * vertex an arbitrary number of (double type) data.
 *
 * \warning Only one element type is allowed. If you need more than one
 * element type in your problem, you can simply create more than one grid
 * object.
 *
 * \warning There is no concept of <i>boundaries</i> in galeri/pfem.
 * Boundaries are indeed defined by independent grid objects.
 */
class Loadable : public core::Object
{
  public:
    // @{ \name Constructors and destructors.
    //! Empty constructor.
    Loadable()
    {}

    //! Constructor with specified Epetra_Comm, number of global elements, etc.
    /*! @param comm [In] communicator object
     *
     *  @param numGlobalElements [In] number of global elements in \c this
     *                              grid object.
     *
     *  @param numMyElements [In] number of local elements in \c this
     *                          grid object.
     *
     *  @param elementType  [In] a string value which defines the element
     *  type. Valid values are: \c Point, \c Segment, \c Triangle, \c Quad,
     *  \c Tet and \c Hex.
     *
     *  @param numElementData  [In] number of additional double-typed data
     *  to be stored on each element.
     *
     *  @param numVertexData  [In] number of additional double-typed data
     *  to be stored on each vertex.
     *
     *  This is just a shortcut for method initialize().
     */
    Loadable(const Epetra_Comm& comm,
             const int numGlobalElements,
             const int numMyElements,
             const string& elementType,
             const int numElementData = 0,
             const int numVertexData = 0)
    { 
      initialize(comm, numGlobalElements, numMyElements, elementType, 0,
                 numElementData, numVertexData);
    }

#if 0
    Loadable(const Epetra_Map& ElementMap, 
             const grid::Element& Element) :
      ElementMap_(new Epetra_Map(ElementMap)),
      GridElement_(Element)
    { 
      ADJ_ = rcp(new EpetraExt::DistArray<int>(*ElementMap_, GridElement_->getNumVertices()));
    }
#endif

    //! Destructor.
    ~Loadable() {}

    //! Initialization method.
    /*! @param comm [In] communicator object
     *
     *  @param numGlobalElements [In] number of global elements in \c this
     *                              grid object.
     *
     *  @param numMyElements [In] number of local elements in \c this
     *                          grid object.
     *
     *  @param elementType  [In] a string value which defines the element
     *  type. Valid values are: \c Point, \c Segment, \c Triangle, \c Quad,
     *  \c Tet and \c Hex.
     *
     *  @param myGlobalElements [In] array of integers, of size \c
     *  numMyElements, which contains the global ID of all local elements.
     *  By using this array, one can introduce global numbering to grid
     *  elements.
     *
     *  @param numElementData  [In] number of additional double-typed data
     *  to be stored on each element.
     *
     *  @param numVertexData  [In] number of additional double-typed data
     *  to be stored on each vertex.
     */
    void initialize(const Epetra_Comm& comm,
                    const int numGlobalElements,
                    const int numMyElements,
                    const string& elementType,
                    const int* myGlobalElements = 0,
                    const int numElementData = 0,
                    const int numVertexData = 0)
    { 
      if (myGlobalElements != 0)
        ElementMap_ = rcp(new Epetra_Map(numGlobalElements, numMyElements, 
                                         myGlobalElements, 0, comm));
      else if (numMyElements != -1)
        ElementMap_ = rcp(new Epetra_Map(numGlobalElements, numMyElements, 0, comm));
      else
        ElementMap_ = rcp(new Epetra_Map(numGlobalElements, 0, comm));

      elementType_ = elementType;

      if (elementType == "Point")
        GridElement_ = rcp(new Galeri::grid::Point);
      else if (elementType == "Segment")
        GridElement_ = rcp(new Galeri::grid::Segment);
      else if (elementType == "Triangle")
        GridElement_ = rcp(new Galeri::grid::Triangle);
      else if (elementType == "Quad")
        GridElement_ = rcp(new Galeri::grid::Quad);
      else if (elementType == "Tet")
        GridElement_ = rcp(new Galeri::grid::Tet);
      else if (elementType == "Hex")
        GridElement_ = rcp(new Galeri::grid::Hex);
      else
        TEST_FOR_EXCEPTION(true, std::logic_error,
                           "input elementType not recognized, " << elementType);

      ADJ_ = rcp(new EpetraExt::DistArray<int>(*ElementMap_, GridElement_->getNumVertices()));

      numElementData_ = numElementData;
      numVertexData_ = numVertexData;

      if (numElementData_ > 0)
        elementData_ = rcp(new Epetra_MultiVector(*ElementMap_, numElementData_));
    }

    // @}
    // @{ \name Get methods.
    
    //! Returns the global number of grid elements in \c this object.
    inline int getNumGlobalElements() const 
    {
      return(ElementMap_->NumGlobalElements());
    }

    //! Returns the local number of grid elements in \c this object.
    inline int getNumMyElements() const 
    {
      return(ElementMap_->NumMyElements());
    }

    //! Returns the global number of grid vertices in \c this object.
    inline int getNumGlobalVertices() const 
    {
      // NOTE: this requires IndexBase == 0
      return(VertexMap_->MaxAllGID() + 1);
    }

    //! Returns the local number of grid vertices in \c this object.
    inline int getNumMyVertices() const 
    {
      return(VertexMap_->NumMyElements());
    }

    //! Returns the number of vertices per element (a constant value).
    inline int getNumVerticesPerElement() const
    {
      return(GridElement_->getNumVertices());
    }

    //! Returns the Epetra_Map associated to the element distribution.
    inline const Epetra_Map getElementMap() const
    {
      return(*ElementMap_);
    }

    //! Returns the Epetra_Map associated to the vertex distribution.
    inline const Epetra_Map getVertexMap() const
    {
      return(*VertexMap_);
    }

    //! Returns the global grid element ID for the specified local (and locally owned) grid element ID.
    inline int getGEID(const int LEID) const
    {
      return(ElementMap_->GID(LEID));
    }

    //! Returns the global grid vertex ID for the specified (and locally owned) local grid vertex ID.
    inline int getGVID(const int LVID) const
    {
      return(VertexMap_->GID(LVID));
    }

    //! Returns the local grid element ID for the specified (and locally owned) global grid element ID.
    inline int getLEID(const int GEID) const
    {
      return(ElementMap_->LID(GEID));
    }

    //! Returns the local grid vertex ID for the specified (and locally owned) global grid vertex ID.
    inline int getLVID(const int GVID) const
    {
      return(VertexMap_->LID(GVID));
    }

    //! Returns the Galeri::grid::Element object of \c this object.
    const grid::Element getElement() const
    {
      return(*GridElement_);
    }

    //! Returns the number of optional double-typed data associated to each grid vertex.
    inline int getNumVertexData() const
    {
      return(numVertexData_);
    }

    //! Returns the number of optional double-typed data associated to each grid element.
    inline int getNumElementData() const
    {
      return(numElementData_);
    }

    // @}
    // @{ Optional element and vertex data

    //! Returns the optional data associated to the specified (and locally owned) global grid element ID, stored in position \c which in the data array.
    inline double getElementData(const int GEID, const int which) const
    {
      int LEID = getLEID(GEID);
      return((*elementData_)[which][LEID]);
    }

    //! Sets the optional data associated to the specified (and locally owned) global grid element ID, are stores it in position \c which in the data array.
    inline void setElementData(const int GEID, const int which, const double val)
    {
      int LEID = getLEID(GEID);
      (*elementData_)[which][LEID] = val;
    }
    //! Returns the optional data associated to the specified (and locally owned) global grid vertex ID, stored in position \c which in the data array.
    inline double getVertexData(const int GVID, const int which) const
    {
      int LVID = getLVID(GVID);
      return((*vertexData_)[which][LVID]);
    }

    //! Sets the optional data associated to the specified (and locally owned) global grid vertex ID, are stores it in position \c which in the data array.
    inline void setVertexData(const int GVID, const int which, const double val)
    {
      int LVID = getLVID(GVID);
      (*vertexData_)[which][LVID] = val;
    }

    // @}
    // @{ \name Data access methods
    
    //! Sets the \c index coordinate of the specified (and locally owned) global grid vertex ID to \c value.
    inline void setGlobalCoordinates(const int GID, const int index, const double value)
    {
      COO_->ReplaceGlobalValue(GID, index, value);
    }

    //! Sets the coordinates of the specified (and locally owned) global grid vertex ID to \c value.
    inline double& getGlobalCoordinates(const int GID, const int index)
    {
      int LID = VertexMap_->LID(GID);
      return((*COO_)[index][LID]);
    }

    //! Sets the \c index coordinate of the specified (and locally owned) local grid vertex ID to \c value.
    inline double& getMyCoordinates(const int LID, const int index)
    {
      return((*COO_)[index][LID]);
    }

    //! Sets the \c index coordinate of the specified (and locally owned) local grid vertex ID to \c value.
    inline void setGlobalConnectivity(const int GID, const int index, const int what)
    {
      int LID = ElementMap_->LID(GID);
      assert (LID != -1);
      (*ADJ_)(LID, index) = what;
    }

    //! Sets the \c index-th component of the specified (and locally owned) global grid element ID to \c value.
    inline int& getGlobalConnectivity(const int GID, const int index)
    {
      int LID = ElementMap_->LID(GID);
      return ((*ADJ_)(LID, index));
    }

    // FIXME???
    //! Sets the \c index-th component of the specified (and locally owned) local grid element ID to \c value.
    inline int& getMyConnectivity(const int LID, const int index)
    {
      return ((*ADJ_)(LID, index));
    }

    // FIXME??
    inline int& ADJ(const int LVID, const int index)
    {
      return((*ADJ_)(LVID, index));
    }

    // @} 
    // @{ \name 
    
    //! Freezes the grid connectivity, which cannot be modified any longer.
    virtual void freezeConnectivity()
    {
      const Epetra_Comm& Comm = ADJ_->Comm();

      // at this point all the elements have been inserted; we look
      // for the vertex map (with overlap). Note that the vertices
      // of all elements are in global numbering.

      int MaxSize = getElementMap().NumMyElements() * getElement().getNumVertices();

      vector<int> MyGlobalElements(MaxSize);
      int count = 0;

      // insert all elements in a hash table
      Teuchos::Hashtable<int, short int> hash(MaxSize * 2);

      for (int i = 0; i < getElementMap().NumMyElements(); ++i)
      {
        for (int j = 0; j < getElement().getNumVertices(); ++j)
        {
          const int& GVID = ADJ(i, j);
          if (hash.containsKey(GVID) == false)
          {
            MyGlobalElements[count++] = GVID;
            hash.put(GVID, 1);
          }
        }
      }

      VertexMap_ = rcp(new Epetra_Map(-1, count, &MyGlobalElements[0], 0,  Comm));

      COO_ = rcp(new Epetra_MultiVector(*VertexMap_, 
                                        Galeri::core::Workspace::getNumDimensions()));

      if (numVertexData_ > 0)
        vertexData_ = rcp(new Epetra_MultiVector(*VertexMap_, numVertexData_));
    }

    //! Freezes the grid coordinates, which cannot be modified any longer.
    void freezeCoordinates()
    {
      // do-nothing at this point
    }

    //! Prints the grid on \c os.
    virtual void print(ostream & os) const
    {
      // FIXME: add label, see parallel, add Barrier()??
      cout << *ElementMap_;

      cout << *VertexMap_;

      cout << *ADJ_;

      cout << *COO_;
    }

    //! Returns the Epetra_Map associated with grid vertices.
    Epetra_Map getLinearVertexMap()
    {
      if (linearVertexMap_ == Teuchos::null)
      {
        linearVertexMap_ = rcp(new Epetra_Map(getNumGlobalVertices(), 0, ElementMap_->Comm()));
      }
      return(*linearVertexMap_);
    }

    //! Returns the Epetra_MultiVector containing the grid coordinates.
    const Epetra_MultiVector& getLinearCoordinates()
    {
      if (linearCOO_ == Teuchos::null)
      {
        Epetra_Map linearVertexMap = getLinearVertexMap();
        linearCOO_ = rcp(new Epetra_MultiVector(linearVertexMap, 
                                                Galeri::core::Workspace::getNumDimensions()));
        linearExporter_ = rcp(new Epetra_Export(getVertexMap(), linearVertexMap));
        linearCOO_->Export(*COO_, *linearExporter_, Insert);
      }
      return(*linearCOO_);
    }

    // FIXME: delete this?
    string getElementType() const
    {
      return(elementType_);
    }

  private:

    // @}
    // @{ \name private methods and data

    RefCountPtr<Epetra_Map> ElementMap_;
    RefCountPtr<Epetra_Map> VertexMap_;
    RefCountPtr<Epetra_Map> linearVertexMap_;
    RefCountPtr<Epetra_Export> linearExporter_;
    RefCountPtr<grid::Element> GridElement_;

    RefCountPtr<Epetra_MultiVector> COO_;
    RefCountPtr<Epetra_MultiVector> linearCOO_;
    RefCountPtr<EpetraExt::DistArray<int> > ADJ_;

    string elementType_;

    int numElementData_;
    int numVertexData_;
    RefCountPtr<Epetra_MultiVector> elementData_;
    RefCountPtr<Epetra_MultiVector> vertexData_;
    // @}

}; // class Loadable

} // namespace grid

}; // namespace Galeri
#endif
