#ifndef PHX_LOADABLE_GRID_H
#define PHX_LOADABLE_GRID_H

#include "Galeri_ConfigDefs.h"
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
#include "Epetra_FEVector.h"
#include "EpetraExt_DistArray.h"

using namespace Teuchos;

namespace phx {

namespace grid {

class Loadable : public core::Object
{
  public:
    // @{ \name Constructors and destructors.
    //! Constructor.
    Loadable(const Epetra_Comm& comm,
             const int numGlobalElements,
             const int numMyElements,
             const string& elementType)
    { 
      ElementMap_ = rcp(new Epetra_Map(numGlobalElements, numMyElements, 0, comm));

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
    }

    Loadable(RefCountPtr<Epetra_Map> ElementMap, 
             RefCountPtr<grid::Element> Element) :
      ElementMap_(ElementMap),
      GridElement_(Element)
    { 
      ADJ_ = rcp(new EpetraExt::DistArray<int>(*ElementMap_, GridElement_->getNumVertices()));
    }

    //! Destructor.
    ~Loadable() {}

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

    // FIXME: THIS IS THE WRONG NUMBER *OVERLAPPING*
    inline int getNumGlobalVertices() const 
    {
      return(VertexMap_->NumGlobalElements());
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
    inline const RefCountPtr<Epetra_Map> getElementMap() const
    {
      return(ElementMap_);
    }

    //! Returns the Epetra_Map associated to the vertex distribution.
    inline const RefCountPtr<Epetra_Map> getVertexMap() const
    {
      return(VertexMap_);
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

    const RefCountPtr<grid::Element> getElement() const
    {
      return(GridElement_);
    }

    //! Returns the number of integer data associated with each element.
    int getIntDataSize() const
    {
      return(IntData_.size());
    }

    void setIntDataSize(const int size)
    {
      IntData_.resize(size);
    }

    //! Returns the \c what integer data vector associated with elements.
    RefCountPtr<Epetra_IntVector> getIntData(const int what)
    {
      return(IntData_[what]);
    }

    void setIntData(const int what, RefCountPtr<Epetra_IntVector> IntVector)
    {
      IntData_[what] = IntVector;
    }

    int getDoubleDataSize() const
    {
      return(DoubleData_.size());
    }

    void setDoubleDataSize(const int size)
    {
      DoubleData_.resize(size);
    }

    //! Returns the \c what double data vector associated with vertices.
    RefCountPtr<Epetra_Vector> getDoubleData(const int what)
    {
      return(DoubleData_[what]);
    }

    void setDoubleData(const int what, RefCountPtr<Epetra_Vector> Vector)
    {
      DoubleData_[what] = Vector;
    }

    // @}
    // @{ \name Data access methods
    
    inline void setGlobalCoordinates(const int GID, const int index, const double value)
    {
      COO_->SumIntoGlobalValue(GID, index, value);
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

      int MaxSize = getElementMap()->NumMyElements() * getElement()->getNumVertices();

      vector<int> MyGlobalElements(MaxSize);
      int count = 0;

      // insert all elements in a hash table
      Teuchos::Hashtable<int, short int> hash(MaxSize * 2);

      for (int i = 0; i < getElementMap()->NumMyElements(); ++i)
      {
        for (int j = 0; j < getElement()->getNumVertices(); ++j)
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

  private:

    // @}
    // @{ \name private methods and data

    RefCountPtr<Epetra_Map> ElementMap_;
    RefCountPtr<Epetra_Map> VertexMap_;
    RefCountPtr<grid::Element> GridElement_;

    vector<RefCountPtr<Epetra_IntVector> > IntData_;
    vector<RefCountPtr<Epetra_Vector> > DoubleData_;

    RefCountPtr<Epetra_MultiVector> COO_;
    RefCountPtr<EpetraExt::DistArray<int> > ADJ_;

    // @}

}; // class Loadable

} // namespace grid

}; // namespace phx
#endif
