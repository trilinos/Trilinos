#ifndef PHX_LOADABLE_GRID_H
#define PHX_LOADABLE_GRID_H

#include "Galeri_ConfigDefs.h"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Hashtable.hpp"

#include "phx_core_Object.h"
#include "phx_grid_Element.h"

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
    Loadable(RefCountPtr<Epetra_Map> ElementMap, 
             RefCountPtr<grid::Element> Element) :
      ElementMap_(ElementMap),
      GridElement_(Element)
    { 
      ADJ_ = rcp(new EpetraExt::DistArray<int>(*ElementMap_, GridElement_->getNumVertices()));
      COO_.resize(GridElement_->getNumDimensions());
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

    inline int getNumGlobalVertices() const 
    {
      return(VertexMap_->NumGlobalElements());
    }

    inline int getNumMyVertices() const 
    {
      return(VertexMap_->NumMyElements());
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
    
    inline void setGlobalCoordinates(const int index, const int NumEntries,
                                     const int* GVID, const double* value)
    {
      Epetra_FEVector& V = *COO_[index];
      V.SumIntoGlobalValues(NumEntries, GVID, value);
    }

    inline void setGlobalCoordinates(const int GID, const int index, const double value)
    {
      Epetra_FEVector& V = *COO_[index];
      V.SumIntoGlobalValues(1, &GID, &value);
    }

    inline double& getGlobalCoordinates(const int GID, const int index)
    {
      Epetra_FEVector& V = *COO_[index];
      int LID = VertexMap_->LID(GID);
      assert (LID != -1);
      return(V[0][LID]);
    }

    inline double& getMyCoordinates(const int LID, const int index)
    {
      Epetra_FEVector& V = *COO_[index];
      return(V[0][LID]);
    }

    inline double& COO(const int GVID, const int index)
    {
      Epetra_FEVector& V = *COO_[index];
      return(V[0][VertexMap_->LID(GVID)]);
    }

    inline void setGlobalConnectivity(const int GID, const int index, const int what)
    {
      int LID = ElementMap_->LID(GID);
      (*ADJ_)(LID, index) = what;
    }

    inline int& getGlobalConnectivity(const int GID, const int index)
    {
      int LID = ElementMap_->LID(GID);
      return ((*ADJ_)(LID, index));
    }

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

      for (int i = 0; i < COO_.size(); ++i)
        COO_[i] = rcp(new Epetra_FEVector(*VertexMap_));
    }

    void freezeCoordinates()
    {
      const Epetra_Comm& Comm = ADJ_->Comm();

      for (int i = 0; i < COO_.size(); ++i)
      {
        COO_[i]->GlobalAssemble(Insert);
      }
    }

    const Epetra_FEVector& COOPtr(const int i) const
    {
      return(*(COO_[i].get()));
    }

    virtual void print(ostream & os) const
    {
      cout << *ADJ_;

      for (int i = 0; i < COO_.size(); ++i)
        cout << *(COO_[i]);
    }

  private:

    // @}
    // @{ \name private methods and data

    RefCountPtr<Epetra_Map> ElementMap_;
    RefCountPtr<Epetra_Map> VertexMap_;
    RefCountPtr<grid::Element> GridElement_;

    vector<RefCountPtr<Epetra_IntVector> > IntData_;
    vector<RefCountPtr<Epetra_Vector> > DoubleData_;

    vector<RefCountPtr<Epetra_FEVector> > COO_;
    RefCountPtr<EpetraExt::DistArray<int> > ADJ_;

    // @}

}; // class Loadable

} // namespace grid

}; // namespace phx
#endif
