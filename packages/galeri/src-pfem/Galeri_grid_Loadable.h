// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
#include "Galeri_core_Object.h"
#include "Galeri_core_Workspace.h"
#include "Galeri_grid_Element.h"

#include "Teuchos_Assert.hpp"
#include "Teuchos_RefCountPtr.hpp"

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Export.h"
#include "Epetra_MultiVector.h"

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
    Loadable() :
      status_(core::Workspace::UNINITIALIZED)
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
    Loadable(const Epetra_Comm& comm,
             const int numGlobalElements,
             const int numMyElements,
             const std::string& elementType,
             const int* myGlobalElements = 0,
             const int numElementData = 0,
             const int numVertexData = 0);

    //! Destructor.
    ~Loadable() 
    {
      status_ = core::Workspace::UNINITIALIZED;
    }

    //! Initialization method.
    /*! @param comm [In] communicator object
     *
     *  @param numGlobalElements [In] number of global elements in \c this
     *                              grid object.
     *
     *  @param numMyElements [In] number of local elements in \c this
     *                          grid object.
     *
     *  @param element  [In] element to be used.
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
     *
     *  \warning both vertexMap and elementMap are built with a base index of
     *  0. The code may not work otherwise.
     */
    void initialize(const Epetra_Comm& comm,
                    const int numGlobalElements,
                    const int numMyElements,
                    const Galeri::grid::Element& element,
                    const int* myGlobalElements = 0,
                    const int numElementData = 0,
                    const int numVertexData = 0);

    // @}
    // @{ \name Get/Set Methods.
    
    //! Returns the communicator of \c this object.
    inline const Epetra_Comm& getComm() const
    {
#ifdef GALERI_CHECK
      TEUCHOS_TEST_FOR_EXCEPTION(status_ == core::Workspace::INITIALIZED, std::logic_exception,
                         "method getComm() called, but the object is " <<
                         "uninitialized");
#endif
      return(elementMap_->Comm());
    }

    //! Returns the global number of grid elements in \c this object.
    inline int getNumGlobalElements() const 
    {
#ifdef GALERI_CHECK
      TEUCHOS_TEST_FOR_EXCEPTION(status_ == core::Workspace::INITIALIZED, std::logic_exception,
                         "method getNumGlobalElements() called, but the object is " <<
                         "uninitialized");
#endif
      return(elementMap_->NumGlobalElements());
    }

    //! Returns the local number of grid elements in \c this object.
    inline int getNumMyElements() const 
    {
#ifdef GALERI_CHECK
      TEUCHOS_TEST_FOR_EXCEPTION(status_ == core::Workspace::INITIALIZED, std::logic_exception,
                         "method getNumMyElements() called, but the object is " <<
                         "uninitialized");
#endif
      return(elementMap_->NumMyElements());
    }

    //! Returns the global number of grid vertices in \c this object.
    inline int getNumGlobalVertices() const 
    {
#ifdef GALERI_CHECK
      TEUCHOS_TEST_FOR_EXCEPTION(status_ == core::Workspace::INITIALIZED, std::logic_exception,
                         "method getNumGlobalVertices() called, but the object is " <<
                         "uninitialized");
#endif
      return(vertexMap_->MaxAllGID() + 1);
    }

    //! Returns the local number of grid vertices in \c this object.
    inline int getNumMyVertices() const 
    {
#ifdef GALERI_CHECK
      TEUCHOS_TEST_FOR_EXCEPTION(status_ == core::Workspace::INITIALIZED, std::logic_exception,
                         "method getNumMyVertices() called, but the object is " <<
                         "uninitialized");
#endif
      return(vertexMap_->NumMyElements());
    }

    /*! \brief Returns the number of vertices per element, which is
     *  constant value across all grid elements.
     */
    inline int getNumVerticesPerElement() const
    {
#ifdef GALERI_CHECK
      TEUCHOS_TEST_FOR_EXCEPTION(status_ == core::Workspace::INITIALIZED, std::logic_exception,
                         "method getNumVerticesPerElement() called, but the object is " <<
                         "uninitialized");
#endif
      return(element_.getNumVertices());
    }

    //! Returns the Epetra_Map associated to the element distribution.
    inline const Epetra_Map getElementMap() const
    {
#ifdef GALERI_CHECK
      TEUCHOS_TEST_FOR_EXCEPTION(status_ == core::Workspace::INITIALIZED, std::logic_exception,
                         "method getElementMap() called, but the object is " <<
                         "uninitialized");
#endif
      return(*elementMap_);
    }

    //! Returns the Epetra_Map associated to the vertex distribution.
    inline const Epetra_Map getVertexMap() const
    {
#ifdef GALERI_CHECK
      TEUCHOS_TEST_FOR_EXCEPTION(status_ == core::Workspace::INITIALIZED, std::logic_exception,
                         "method getVertexMap() called, but the object is " <<
                         "uninitialized");
#endif
      return(*vertexMap_);
    }

    //! Returns the global grid element ID for the specified local (and locally owned) grid element ID.
    inline int getGEID(const int LEID) const
    {
#ifdef GALERI_CHECK
      TEUCHOS_TEST_FOR_EXCEPTION(status_ == core::Workspace::INITIALIZED, std::logic_exception,
                         "method getGEID() called, but the object is " <<
                         "uninitialized");
#endif
      return(elementMap_->GID(LEID));
    }

    //! Returns the global grid vertex ID for the specified (and locally owned) local grid vertex ID.
    inline int getGVID(const int LVID) const
    {
#ifdef GALERI_CHECK
      TEUCHOS_TEST_FOR_EXCEPTION(status_ == core::Workspace::INITIALIZED, std::logic_exception,
                         "method getGVID() called, but the object is " <<
                         "uninitialized");
#endif
      return(vertexMap_->GID(LVID));
    }

    //! Returns the local grid element ID for the specified (and locally owned) global grid element ID.
    inline int getLEID(const int GEID) const
    {
#ifdef GALERI_CHECK
      TEUCHOS_TEST_FOR_EXCEPTION(status_ == core::Workspace::INITIALIZED, std::logic_exception,
                         "method getLEID() called, but the object is " <<
                         "uninitialized");
#endif
      return(elementMap_->LID(GEID));
    }

    //! Returns the local grid vertex ID for the specified (and locally owned) global grid vertex ID.
    inline int getLVID(const int GVID) const
    {
#ifdef GALERI_CHECK
      TEUCHOS_TEST_FOR_EXCEPTION(status_ == core::Workspace::INITIALIZED, std::logic_exception,
                         "method getLVID() called, but the object is " <<
                         "uninitialized");
#endif
      return(vertexMap_->LID(GVID));
    }

    //! Returns the Galeri::grid::Element object of \c this object.
    const grid::Element getElement() const
    {
#ifdef GALERI_CHECK
      TEUCHOS_TEST_FOR_EXCEPTION(status_ == core::Workspace::INITIALIZED, std::logic_exception,
                         "method getElement() called, but the object is " <<
                         "uninitialized");
#endif
      return(element_);
    }

    //! Returns the number of optional double-typed data associated to each grid vertex.
    inline int getNumVertexData() const
    {
#ifdef GALERI_CHECK
      TEUCHOS_TEST_FOR_EXCEPTION(status_ == core::Workspace::INITIALIZED, std::logic_exception,
                         "method getNumVertexData() called, but the object is " <<
                         "uninitialized");
#endif
      return(numVertexData_);
    }

    //! Returns the number of optional double-typed data associated to each grid element.
    inline int getNumElementData() const
    {
#ifdef GALERI_CHECK
      TEUCHOS_TEST_FOR_EXCEPTION(status_ == core::Workspace::INITIALIZED, std::logic_exception,
                         "method getNumElementData() called, but the object is " <<
                         "uninitialized");
#endif
      return(numElementData_);
    }

    // @}
    // @{ Optional element and vertex data

    //! Returns the optional data associated to the specified (and locally owned) global grid element ID, stored in position \c which in the data array.
    inline double getElementData(const int GEID, const int which) const
    {
#ifdef GALERI_CHECK
      TEUCHOS_TEST_FOR_EXCEPTION(status_ == core::Workspace::INITIALIZED, std::logic_exception,
                         "method getElementData() called, but the object is " <<
                         "uninitialized");
#endif
      int LEID = getLEID(GEID);
      // FIXME
      return((*elementData_)[which][LEID]);
    }

    //! Sets the optional data associated to the specified (and locally owned) global grid element ID, are stores it in position \c which in the data array.
    inline void setElementData(const int GEID, const int which, const double val)
    {
#ifdef GALERI_CHECK
      TEUCHOS_TEST_FOR_EXCEPTION(status_ == core::Workspace::INITIALIZED, std::logic_exception,
                         "method setElementData() called, but the object is " <<
                         "uninitialized");
#endif
      int LEID = getLEID(GEID);
      // FIXME
      (*elementData_)[which][LEID] = val;
    }

    //! Returns the optional data associated to the specified (and locally owned) global grid vertex ID, stored in position \c which in the data array.
    inline double getVertexData(const int GVID, const int which) const
    {
#ifdef GALERI_CHECK
      TEUCHOS_TEST_FOR_EXCEPTION(status_ != core::Workspace::CONNECTIVITY_FREEZED, std::logic_exception,
                         "method setVertexData() called, but freezeConnectivity() has not been called");
#endif
      int LVID = getLVID(GVID);
      // FIXME
      return((*vertexData_)[which][LVID]);
    }

    //! Sets the optional data associated to the specified (and locally owned) global grid vertex ID, are stores it in position \c which in the data array.
    inline void setVertexData(const int GVID, const int which, const double val)
    {
#ifdef GALERI_CHECK
      TEUCHOS_TEST_FOR_EXCEPTION(status_ == core::Workspace::INITIALIZED, std::logic_exception,
                         "method setVertexData() called, but the object is " <<
                         "uninitialized");
#endif
      int LVID = getLVID(GVID);
      // FIXME
      (*vertexData_)[which][LVID] = val;
    }

    // @}
    // @{ \name Data access methods
    
    //! Sets the \c index coordinate of the specified (and locally owned) global grid vertex ID to \c value.
    inline void setGlobalCoordinates(const int GID, const int index, const double value)
    {
#ifdef GALERI_CHECK
      TEUCHOS_TEST_FOR_EXCEPTION(status_ == core::Workspace::INITIALIZED, std::logic_exception,
                         "method setGlobalCoordinates() called, but the object is " <<
                         "uninitialized");

      TEUCHOS_TEST_FOR_EXCEPTION(status_ == core::Workspace::COORDINATES_FREEZED, std::logic_exception,
                         "method setGlobalCoordinates() called, but the coordinates are " <<
                         "already freezed");
#endif
      COO_->ReplaceGlobalValue(GID, index, value);
    }

    //! Sets the coordinates of the specified (and locally owned) global grid vertex ID to \c value.
    inline double& getGlobalCoordinates(const int GID, const int index)
    {
#ifdef GALERI_CHECK
      TEUCHOS_TEST_FOR_EXCEPTION(status_ == core::Workspace::INITIALIZED, std::logic_exception,
                         "method getGlobalCoordinates() called, but the object is " <<
                         "uninitialized");
#endif
      int LID = vertexMap_->LID(GID);
      return((*COO_)[index][LID]);
    }

    //! Sets the \c index coordinate of the specified (and locally owned) local grid vertex ID to \c value.
    inline double& getMyCoordinates(const int LID, const int index)
    {
#ifdef GALERI_CHECK
      TEUCHOS_TEST_FOR_EXCEPTION(status_ == core::Workspace::INITIALIZED, std::logic_exception,
                         "method getMyCoordinates() called, but the object is " <<
                         "uninitialized");
#endif
      return((*COO_)[index][LID]);
    }

    //! Sets the \c index coordinate of the specified (and locally owned) local grid vertex ID to \c value. (const version)
    inline const double& getMyCoordinates(const int LID, const int index) const
    {
#ifdef GALERI_CHECK
      TEUCHOS_TEST_FOR_EXCEPTION(status_ == core::Workspace::INITIALIZED, std::logic_exception,
                         "method getMyCoordinates() called, but the object is " <<
                         "uninitialized");
#endif
      return((*COO_)[index][LID]);
    }

    //! Sets the \c index coordinate of the specified (and locally owned) local grid vertex ID to \c value.
    inline void setGlobalConnectivity(const int GID, const int index, const int what)
    {
#ifdef GALERI_CHECK
      TEUCHOS_TEST_FOR_EXCEPTION(status_ == core::Workspace::INITIALIZED, std::logic_exception,
                         "method setGlobalConnectivity() called, but the object is " <<
                         "uninitialized");
#endif
      int LID = elementMap_->LID(GID);
#ifdef GALERI_CHECK
      TEUCHOS_TEST_FOR_EXCEPTION(LID == -1, std::logic_exception,
                         "Requested GID (" << GID << ") is not locally owned");
#endif
      (*CON_)(LID, index) = what;
    }

    //! Sets the \c index-th component of the specified (and locally owned) global grid element ID to \c value.
    inline int& getGlobalConnectivity(const int GID, const int index)
    {
#ifdef GALERI_CHECK
      TEUCHOS_TEST_FOR_EXCEPTION(status_ == core::Workspace::INITIALIZED, std::logic_exception,
                         "method getGlobalConnectivity() called, but the object is " <<
                         "uninitialized");
#endif
      int LID = elementMap_->LID(GID);
#ifdef GALERI_CHECK
      TEUCHOS_TEST_FOR_EXCEPTION(LID == -1, std::logic_exception,
                         "Requested GID (" << GID << ") is not locally owned");
#endif
      return ((*CON_)(LID, index));
    }

    //! Sets the \c index-th component of the specified (and locally owned) local grid element ID to \c value.
    inline int& getMyConnectivity(const int LID, const int index)
    {
#ifdef GALERI_CHECK
      TEUCHOS_TEST_FOR_EXCEPTION(status_ == core::Workspace::INITIALIZED, std::logic_exception,
                         "method getMyConnectivity() called, but the object is " <<
                         "uninitialized");
#endif
      return ((*CON_)(LID, index));
    }

    //! Sets the \c index-th component of the specified (and locally owned) local grid element ID to \c value. (const version)
    inline const int& getMyConnectivity(const int LID, const int index) const
    {
#ifdef GALERI_CHECK
      TEUCHOS_TEST_FOR_EXCEPTION(status_ == core::Workspace::INITIALIZED, std::logic_exception,
                         "method getMyConnectivity() called, but the object is " <<
                         "uninitialized");
#endif
      return ((*CON_)(LID, index));
    }

    // @} 
    // @{ \name 
    
    //! Freezes the grid connectivity, which cannot be modified any longer.
    /*! 
     * This method constructs the set of locally owned vertices, by looking
     * over all local elements, and building the vertexMap of the grid.
     * After having called this method, it is no longer possible to modify the
     * grid connectivity, because this would change the vertexMap as well.
     *
     * \warning Performances to be improved for large problems.
     */
    virtual void freezeConnectivity();

    //! Freezes the grid coordinates, which cannot be modified any longer.
    void freezeCoordinates();

    //! Prints the grid on \c os.
    virtual void print(ostream & os) const;

    //! Returns the Epetra_Map associated with grid vertices; each vertex is owned by exactly one processor.
    const Epetra_Map& getNonOverlappingVertexMap();

    //! Returns the Epetra_MultiVector containing the coordinates of vertices; each vertex is owned by exactly one processor.
    const Epetra_MultiVector& getNonOverlappingCoordinates();

    //! Returns a linear Epetra_Map for grid vertices.
    const Epetra_Map& getLinearVertexMap();

    //! Returns the coordinates as a vector based on linearVertexMap.
    const Epetra_MultiVector& getLinearCoordinates();

  private:

    // @}
    // @{ \name private methods and data

    //! Status of the object (UNINITILAZED, INITIALIZED, CONNECTIVITY_FREEZED, COORDINATES_FREEZED).
    int status_;
    //! Map for element distribution. A given element is assigned to exactly one processor.
    RefCountPtr<Epetra_Map> elementMap_;
    //! Map for vertex distribution. A given vertex may be assigned to more than one processor.
    RefCountPtr<Epetra_Map> vertexMap_;

    //! A map for vertices. A given vertex is assigned to exactly one processor.
    RefCountPtr<Epetra_Map> nonOverlappingVertexMap_;
    //! Exporter to nonOverlappingVertexMap.
    RefCountPtr<Epetra_Export> nonOverlappingExporter_;

    //! A linear map for vertices. A given vertex is assigned to exactly one processor.
    RefCountPtr<Epetra_Map> linearVertexMap_;
    //! Exporter to linearVertexMap.
    RefCountPtr<Epetra_Export> linearExporter_;

    //! Element used in \c this object.
    grid::Element element_;
    //! Container for vertex coordinates, based on vertexMap_.
    RefCountPtr<Epetra_MultiVector> COO_;
    //! Container for vertex coordinates, based on nonOverlappingVertexMap_.
    RefCountPtr<Epetra_MultiVector> nonOverlappingCOO_;
    //! Container for vertex coordinates, based on linearVertexMap_.
    RefCountPtr<Epetra_MultiVector> linearCOO_;

    //! Container for element connectivity, based on elementMap_.
    RefCountPtr<EpetraExt::DistArray<int> > CON_;
    //! Amount of additional data assigned to each grid element.
    int numElementData_;
    //! Amount of additional data assigned to each grid vertex.
    int numVertexData_;
    //! Additional data assigned to grid elements.
    RefCountPtr<Epetra_MultiVector> elementData_;
    //! Additional data assigned to grid vertices.
    RefCountPtr<Epetra_MultiVector> vertexData_;

    // @}

}; // class Loadable

} // namespace grid

}; // namespace Galeri
#endif
