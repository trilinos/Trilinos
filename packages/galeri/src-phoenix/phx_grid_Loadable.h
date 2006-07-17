// @HEADER
// ************************************************************************
//
//                  Galeri Matrix Generation Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef PHX_LOADABLE_GRID_H
#define PHX_LOADABLE_GRID_H

#include "../src/Galeri_ConfigDefs.h"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Hashtable.hpp"

#include "phx_core_Object.h"
#include "phx_core_Utils.h"
#include "phx_grid_Element.h"
#include "phx_grid_Point.h"
#include "phx_grid_Segment.h"
#include "phx_grid_Triangle.h"
#include "phx_grid_Quad.h"
#include "phx_grid_Tet.h"
#include "phx_grid_Hex.h"

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Export.h"

#include "EpetraExt_DistArray.h"

using namespace Teuchos;

namespace phx {

namespace grid {

class Loadable : public core::Object
{
  public:
    // @{ \name Constructors and destructors.
    //! Constructor.
    Loadable()
    {}

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
        GridElement_ = rcp(new phx::grid::Point);
      else if (elementType == "Segment")
        GridElement_ = rcp(new phx::grid::Segment);
      else if (elementType == "Triangle")
        GridElement_ = rcp(new phx::grid::Triangle);
      else if (elementType == "Quad")
        GridElement_ = rcp(new phx::grid::Quad);
      else if (elementType == "Tet")
        GridElement_ = rcp(new phx::grid::Tet);
      else if (elementType == "Hex")
        GridElement_ = rcp(new phx::grid::Hex);
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
    
    inline int getNumGlobalElements() const 
    {
      return(ElementMap_->NumGlobalElements());
    }

    inline int getNumMyElements() const 
    {
      return(ElementMap_->NumMyElements());
    }

    inline int getNumGlobalVertices() const 
    {
      // NOTE: this requires IndexBase == 0
      return(VertexMap_->MaxAllGID() + 1);
    }

    inline int getNumMyVertices() const 
    {
      return(VertexMap_->NumMyElements());
    }

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

    inline int getGEID(const int LEID) const
    {
      return(ElementMap_->GID(LEID));
    }

    inline int getGVID(const int LVID) const
    {
      return(VertexMap_->GID(LVID));
    }

    inline int getLEID(const int GEID) const
    {
      return(ElementMap_->LID(GEID));
    }

    inline int getLVID(const int GVID) const
    {
      return(VertexMap_->LID(GVID));
    }

    const grid::Element getElement() const
    {
      return(*GridElement_);
    }

    inline int getNumVertexData() const
    {
      return(numVertexData_);
    }

    inline int getNumElementData() const
    {
      return(numElementData_);
    }

    inline double getVertexData(const int GVID, const int which) const
    {
      int LVID = getLVID(GVID);
      return((*vertexData_)[which][LVID]);
    }

    inline void setVertexData(const int GVID, const int which, const double val)
    {
      int LVID = getLVID(GVID);
      (*vertexData_)[which][LVID] = val;
    }

    inline double getElementData(const int GEID, const int which) const
    {
      int LEID = getLEID(GEID);
      return((*elementData_)[which][LEID]);
    }

    inline void setElementData(const int GEID, const int which, const double val)
    {
      int LEID = getLEID(GEID);
      (*elementData_)[which][LEID] = val;
    }

    // @}
    // @{ \name Data access methods
    
    inline void setGlobalCoordinates(const int GID, const int index, const double value)
    {
      COO_->ReplaceGlobalValue(GID, index, value);
    }

    inline double& getGlobalCoordinates(const int GID, const int index)
    {
      int LID = VertexMap_->LID(GID);
      return((*COO_)[index][LID]);
    }

    inline double& getMyCoordinates(const int LID, const int index)
    {
      return((*COO_)[index][LID]);
    }

    inline void setGlobalConnectivity(const int GID, const int index, const int what)
    {
      int LID = ElementMap_->LID(GID);
      assert (LID != -1);
      (*ADJ_)(LID, index) = what;
    }

    inline int& getGlobalConnectivity(const int GID, const int index)
    {
      int LID = ElementMap_->LID(GID);
      return ((*ADJ_)(LID, index));
    }

    // FIXME???
    inline int& getMyConnectivity(const int LID, const int index)
    {
      return ((*ADJ_)(LID, index));
    }

    inline int& ADJ(const int LVID, const int index)
    {
      return((*ADJ_)(LVID, index));
    }

    // @} 
    // @{ \name 
    
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
                                        phx::core::Utils::getNumDimensions()));

      if (numVertexData_ > 0)
        vertexData_ = rcp(new Epetra_MultiVector(*VertexMap_, numVertexData_));
    }

    void freezeCoordinates()
    {
      // do-nothing at this point
    }

    virtual void print(ostream & os) const
    {
      cout << *ElementMap_;

      cout << *VertexMap_;

      cout << *ADJ_;

      cout << *COO_;
    }

    Epetra_Map getLinearVertexMap()
    {
      if (linearVertexMap_ == Teuchos::null)
      {
        linearVertexMap_ = rcp(new Epetra_Map(getNumGlobalVertices(), 0, ElementMap_->Comm()));
      }
      return(*linearVertexMap_);
    }

    const Epetra_MultiVector& getLinearCoordinates()
    {
      if (linearCOO_ == Teuchos::null)
      {
        Epetra_Map linearVertexMap = getLinearVertexMap();
        linearCOO_ = rcp(new Epetra_MultiVector(linearVertexMap, 
                                                phx::core::Utils::getNumDimensions()));
        linearExporter_ = rcp(new Epetra_Export(getVertexMap(), linearVertexMap));
        linearCOO_->Export(*COO_, *linearExporter_, Insert);
      }
      return(*linearCOO_);
    }

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

}; // namespace phx
#endif
