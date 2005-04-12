#ifndef ML_BASICGRID_H
#define ML_BASICGRID_H

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_IntSerialDenseMatrix.h"
#include "Teuchos_RefCountPtr.hpp"
#include "ml_AbstractGrid.h"
#include "MLAPI_Workspace.h"
#include <vector>
#include <algorithm>

using namespace std;
using namespace Teuchos;
using namespace ML_Epetra;
using namespace MLAPI;

/*
 * \class TriangleRectangleGrid
 *
 * \brief Creates a grid with triangles in a rectangle.
 *
 * \author Marzio Sala, SNL 9214.
 *
 * \date Last updated on 03-Apr-05.
 */
class BasicGrid : public AbstractGrid
{

public:

  //! Constructor.
  BasicGrid(Epetra_Comm& Comm, const string ElementType,
            const int EstimatedNumMyElements = 0,
            const int EstimatedNumMyVertices = 0,
            const int EstimatedNumMyBoundaryFaces = 0) :
    ElementType_(ElementType),
    Comm_(Comm),
    EstimatedNumMyElements_(EstimatedNumMyElements),
    EstimatedNumMyVertices_(EstimatedNumMyVertices),
    EstimatedNumMyBoundaryFaces_(EstimatedNumMyBoundaryFaces)
  {
    if (ElementType_ == "ML_TRIANGLE")
    {
      NumDimensions_ = 2;
      NumVerticesPerElement_ = 3;
      NumFacesPerElement_ = 3;
      NumVerticesPerFace_ = 2;
    }
    else if (ElementType_ == "ML_QUAD")
    {
      NumDimensions_ = 2;
      NumVerticesPerElement_ = 4;
      NumFacesPerElement_ = 4;
      NumVerticesPerFace_ = 2;
    }
    else if (ElementType_ == "ML_TET")
    {
      NumDimensions_ = 3;
      NumVerticesPerElement_ = 4;
      NumFacesPerElement_ = 4;
      NumVerticesPerFace_ = 3;
    }
    else if (ElementType_ == "ML_HEX")
    {
      NumDimensions_ = 3;
      NumVerticesPerElement_ = 8;
      NumFacesPerElement_ = 6;
      NumVerticesPerFace_ = 4;
    }
  }

  //! Destructor
  virtual ~BasicGrid() {}

  virtual int NumDimensions() const
  {
    return(NumDimensions_);
  }

  virtual int NumVerticesPerElement() const
  {
    return(NumVerticesPerElement_);
  }

  virtual int NumFacesPerElement() const
  {
    return(NumFacesPerElement_);
  }

  virtual int NumVerticesPerFace() const
  {
    return(NumVerticesPerFace_);
  }

  virtual string ElementType() const
  {
    return(ElementType_);
  }

  virtual const Epetra_Comm& Comm() const
  {
    return(Comm_);
  }

  virtual int NumMyElements() const
  {
    return(NumMyElements_);
  }

  virtual int NumGlobalElements() const
  {
    return(NumGlobalElements_);
  }

  virtual int NumMyVertices() const
  {
    return(NumMyVertices_);
  }

  virtual int NumGlobalVertices() const
  {
    return(NumGlobalVertices_);
  }

  virtual int NumMyBoundaryFaces() const
  {
    return(NumMyBoundaryFaces_);
  }

  virtual int NumGlobalBoundaryFaces() const
  {
    return(NumGlobalBoundaryFaces_);
  }

  virtual void VertexCoord(const int LocalID, double* coord) const
  {
  }

  virtual void VertexCoord(const int Length, const int* IDs, 
                           double* x, double* y, double* z) const
  {
  }

  virtual void ElementVertices(const int LocalID, int* elements) const
  {
  }

  virtual double ElementMinLength(const int LocalElement) const
  {
    ML_THROW("NOT IMPL", -1);
  }

  virtual double ElementMaxLength(const int LocalElement) const
  {
    ML_THROW("NOT IMPL", -1);
  }

  virtual const RefCountPtr<Epetra_Map> RCPVertexMap() const
  {
    return(VertexMap_);
  }

  virtual const RefCountPtr<Epetra_Map> RCPElementMap() const
  {
    return(ElementMap_);
  }

  virtual const Epetra_Map& VertexMap() const
  {
    return(*(VertexMap_.get()));
  }

  virtual const Epetra_Map& ElementMap() const
  {
    return(*(ElementMap_.get()));
  }

  virtual const Epetra_Map& RowMap() const
  {
    return(*(RowMap_.get()));
  }

  virtual const Epetra_Import& Importer() const
  {
    return(*(Importer_.get()));
  }

  virtual int ElementTag(const int ID) const
  {
    return(1);
  }

  virtual int VertexTag(const int ID) const
  {
    return(1);
  }

  virtual double ElementVolume() const
  {
    return(DeltaX() * DeltaY());
  }

  virtual void FaceVertices(const int LocalFace, int& tag, int* IDs) const
  {
  }

  inline int FacePatch(const int LocalFace) const
  {
  }

  virtual double ElementVolume(const int LocalElement) const
  {
    ML_THROW("NOT IMPLT", -1);
  }

  virtual double FaceArea(const int LocalFace) const
  {
  }

  virtual double MyVolume() const
  {
    ML_THROW("NOT IMPLT", -1);
  }

  virtual double GlobalVolume() const
  {
    ML_THROW("NOT IMPLT", -1);
  }

  void ExportToVertexMap(const Epetra_DistObject& RowObject,
                         Epetra_DistObject& VertexObject) const
  {
    VertexObject.Import(RowObject, Importer(), Insert);
  }

  void ExportToRowMap(const Epetra_DistObject& VertexObject,
                            Epetra_DistObject& RowObject) const
  {
    RowObject.Export(VertexObject, Importer(), Insert);
  }

private:

  Epetra_Comm& Comm_;

  int NumMyVertices_;
  int NumGlobalVertices_;
  int NumMyElements_;
  int NumGlobalElements_;
  int NumMyBoundaryFaces_;
  int NumGlobalBoundaryFaces_;

  RefCountPtr<Epetra_Map> VertexMap_;
  RefCountPtr<Epetra_Map> ElementMap_;
  RefCountPtr<Epetra_Map> RowMap_;
  RefCountPtr<Epetra_Import> Importer_;
};

#endif
