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

#include "Galeri_grid_Loadable.h"
#include "Galeri_core_Object.h"
#include "Galeri_core_Workspace.h"
#include "Galeri_grid_Element.h"
#include "Galeri_grid_Point.h"
#include "Galeri_grid_Segment.h"
#include "Galeri_grid_Triangle.h"
#include "Galeri_grid_Quad.h"
#include "Galeri_grid_Tet.h"
#include "Galeri_grid_Hex.h"

#include "Teuchos_HashTable.hpp"

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Export.h"

#include "EpetraExt_DistArray.h"

using namespace Teuchos;

// ============================================================================ 
Galeri::grid::Loadable::Loadable(const Epetra_Comm& comm,
                                 const int numGlobalElements,
                                 const int numMyElements,
                                 const string& elementType,
                                 const int numElementData,
                                 const int numVertexData) :
  status_(core::Workspace::UNINITIALIZED)
{ 
  Element element;

  if (elementType == "Point")
  {
    Point point;
    element = point;
  }
  else if (elementType == "Segment")
  {
    Segment segment;
    element = segment;
  }
  else if (elementType == "Triangle")
  {
    Triangle triangle;
    element = triangle;
  }
  else if (elementType == "Quad")
  {
    Quad quad;
    element = quad;
  }
  else if (elementType == "Tet")
  {
    Tet tet;
    element = tet;
  }
  else if (elementType == "Hex")
  {
    Hex hex;
    element = hex;
  }
  else
    TEST_FOR_EXCEPTION(true, std::logic_error,
                       "input elementType not recognized, " << elementType);

  initialize(comm, numGlobalElements, numMyElements, element, 0,
             numElementData, numVertexData);
}

// ============================================================================
void Galeri::grid::Loadable::initialize(const Epetra_Comm& comm,
                                        const int numGlobalElements,
                                        const int numMyElements,
                                        const Galeri::grid::Element& element,
                                        const int* myGlobalElements,
                                        const int numElementData,
                                        const int numVertexData)
{ 
  status_ = core::Workspace::UNINITIALIZED;

  element_ = element;

  if (myGlobalElements != 0)
    elementMap_ = rcp(new Epetra_Map(numGlobalElements, numMyElements, 
                                     myGlobalElements, 0, comm));
  else if (numMyElements != -1)
    elementMap_ = rcp(new Epetra_Map(numGlobalElements, numMyElements, 0, comm));
  else
    elementMap_ = rcp(new Epetra_Map(numGlobalElements, 0, comm));

  CON_ = rcp(new EpetraExt::DistArray<int>(*elementMap_, element_.getNumVertices()));

  numElementData_ = numElementData;
  numVertexData_ = numVertexData;

  if (numElementData_ > 0)
    elementData_ = rcp(new Epetra_MultiVector(*elementMap_, numElementData_));

  status_ = core::Workspace::INITIALIZED;
}

// ============================================================================ 
void Galeri::grid::Loadable::freezeConnectivity()
{
#ifdef GALERI_CHECK
  TEST_FOR_EXCEPTION(status_ == core::Workspace::INITIALIZED, std::logic_exception,
                     "method freezeConnectivity() called, but the object is " <<
                     "uninitialized");
#endif

  const Epetra_Comm& Comm = CON_->Comm();

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
      const int& GVID = getMyConnectivity(i, j);
      if (hash.containsKey(GVID) == false)
      {
        MyGlobalElements[count++] = GVID;
        hash.put(GVID, 1);
      }
    }
  }

  vertexMap_ = rcp(new Epetra_Map(-1, count, &MyGlobalElements[0], 0,  Comm));

  COO_ = rcp(new Epetra_MultiVector(*vertexMap_, 
                                    Galeri::core::Workspace::getNumDimensions()));

  if (numVertexData_ > 0)
    vertexData_ = rcp(new Epetra_MultiVector(*vertexMap_, numVertexData_));

  status_ = core::Workspace::CONNECTIVITY_FREEZED;
}

// ============================================================================ 
void Galeri::grid::Loadable::freezeCoordinates()
{
#ifdef GALERI_CHECK
  TEST_FOR_EXCEPTION(status_ != core::Workspace::CONNECTIVITY_FREEZED, std::logic_exception,
                     "method freezeCoordinates() called, but freezeCoordinates() " <<
                     "has not been called");
#endif
  // do-nothing at this point
  status_ = core::Workspace::COORDINATES_FREEZED;
}

// ============================================================================
void Galeri::grid::Loadable::print(ostream & os) const
{
  if (status_ == core::Workspace::UNINITIALIZED)
  {
    cout << endl;
    cout << "** Galeri::grid::Loadable object, not initialized" << endl;
    cout << endl;
    return;
  }
  else
  {
    int k = getElement().getNumVertices();

    const Epetra_Comm& comm = elementMap_->Comm();
    if (comm.MyPID() == 0)
    {
      cout << endl;
      cout << "** Galeri::grid::Loadable object" << endl;
      cout << "** Label: " << getLabel() << endl;
      cout << "** Number of processors: " << comm.NumProc() << endl;
      cout << "** Global number of elements: " << getNumMyElements() << endl;
      cout << "** Global number of vertices: " << getNumMyVertices() << endl;
      cout << "** Label of grid element: " << getElement().getLabel() << endl;
      cout << "** Number of data assigned to elements: " << getNumElementData() << endl;
      cout << "** Number of data assigned to vertices: " << getNumVertexData() << endl;
      cout << endl;
    }

    for (int iproc = 0; iproc < comm.NumProc(); ++iproc)
    { 
      if (iproc == comm.MyPID())
      {
        if (comm.MyPID() == 0)
        {
          cout << "** Element connectivity ** " << endl;
          cout << endl;
          cout << setw(8) << "ProcID";
          cout << setw(10) << "LEID";
          cout << setw(10) << "GEID" << "    ";
          for (int i = 0; i < k; ++i)
            cout << setw(8) << "GVID";
          cout << endl;
        }

        for (int i = 0; i < getNumMyElements(); ++i)
        {
          cout << setw(8) << comm.MyPID();
          cout << setw(10) << i;
          cout << setw(10) << getGEID(i) << "    ";
          for (int j = 0; j < k; ++j)
          {
            cout << setw(8) << getMyConnectivity(i, j);
          }
          cout << endl;
        }
      }
      comm.Barrier();
    }

    int d = core::Workspace::getNumDimensions();

    for (int iproc = 0; iproc < comm.NumProc(); ++iproc)
    { 
      if (iproc == comm.MyPID())
      {
        if (comm.MyPID() == 0)
        {
          cout << endl;
          cout << "** Vertex coordinates ** " << endl;
          cout << endl;
          cout << setw(8) << "ProcID" << setw(10) << "LVID" << setw(10) << "GVID";
          if (d == 1)
            cout << setw(15) << "x-coord" << endl;
          else if (d == 2)
            cout << setw(15) << "x-coord" << setw(15) << "y-coord" << endl;
          else if (d == 3)
            cout << setw(15) << "x-coord" << setw(15) << "y-coord" << setw(15) << "z-coord" << endl;
          cout << endl;
        }

        for (int i = 0; i < getNumMyVertices(); ++i)
        {
          cout << setw(8) << comm.MyPID();
          cout << setw(10) << i;
          cout << setw(10) << getGVID(i);
          if (d == 1)
            cout << setw(15) << getMyCoordinates(i, 0) << endl;
          else if (d == 2)
            cout << setw(15) << getMyCoordinates(i, 0) 
              << setw(15) << getMyCoordinates(i, 1) << endl;
          else if (d == 3)
            cout << setw(15) << getMyCoordinates(i, 0) 
              << setw(15) << getMyCoordinates(i, 1)
              << setw(15) << getMyCoordinates(i, 2) << endl;
        }
      }
      comm.Barrier();
    }
  }
}

// ============================================================================ 
const Epetra_Map& Galeri::grid::Loadable::getLinearVertexMap()
{
#ifdef GALERI_CHECK
  TEST_FOR_EXCEPTION(status_ == core::Workspace::INITIALIZED, std::logic_exception,
                     "method getLinearVertexMap() called, but the object is " <<
                     "uninitialized");
#endif

  if (linearVertexMap_ == Teuchos::null)
  {
    linearVertexMap_ = rcp(new Epetra_Map(getNumGlobalVertices(), 0, elementMap_->Comm()));
  }
  return(*linearVertexMap_);
}

// ============================================================================ 
const Epetra_MultiVector& Galeri::grid::Loadable::getLinearCoordinates()
{
#ifdef GALERI_CHECK
  TEST_FOR_EXCEPTION(status_ == core::Workspace::INITIALIZED, std::logic_exception,
                     "method getLinearVertexMap() called, but the object is " <<
                     "uninitialized");
#endif

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
