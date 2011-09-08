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

#ifndef GALERI_FILEGRID_H
#define GALERI_FILEGRID_H

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_DistObject.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_IntSerialDenseMatrix.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Galeri_Utils.h"
#include "Galeri_Workspace.h"
#include "Galeri_AbstractGrid.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>

using namespace Teuchos;

namespace Galeri {
namespace FiniteElements {

/*!
 * \class Galeri_FileGrid
 *
 * \brief Reads a grid from file, in Galeri format.
 *
 * \author Marzio Sala, ETHZ/COLAB.
 *
 * \date Last updated on 15-Sep-05.
 */
class FileGrid : public AbstractGrid
{

public:

  //! Constructor.
  /*! 
   * \param Comm - (In) Communicator object.
   *
   * \param FileName - (In) Name of grid file.
   */
  FileGrid(const Epetra_Comm& Comm, const string FileName) :
    Comm_(Comm),
    FileName_(FileName),
    NumDimensions_(-1),
    NumVerticesPerElement_(-1),
    NumVerticesPerFace_(-1),
    NumFacesPerElement_(-1),
    NumMyVertices_(-1),
    NumGlobalVertices_(-1),
    NumMyElements_(-1),
    NumGlobalElements_(-1),
    NumMyBoundaryFaces_(-1),
    NumGlobalBoundaryFaces_(-1),
    ElementType_("NOT SET")
  {
    // works only in serial at this point
    if (Comm.NumProc() != 1)
      throw(Exception(__FILE__, __LINE__,
                      "FileGrid can be used only in serial"));

    // read the file. A description of the file format is reported above.

    std::ifstream GridFile;
    GridFile.open(FileName.c_str());
    if (!GridFile.good())
      throw(Exception(__FILE__, __LINE__,
                      "Error opening file `" + FileName + ")"));

    bool IsVertexCoordRead = false;
    bool IsElementVerticesRead = false;
    bool IsFaceVerticesRead = false;

    // STILL TO DO: ELEMENT TAG, VERTEX TAG, FACE TAG
    string what;
    while (!GridFile.eof())
    {
      GridFile >> what;
      if (what == "NumDimensions:") 
        GridFile >> NumDimensions_;
      else if (what == "NumVerticesPerElement:") 
        GridFile >> NumVerticesPerElement_;
      else if (what == "NumVerticesPerFace:") 
        GridFile >> NumVerticesPerFace_;
      else if (what == "NumFacesPerElement:") 
        GridFile >> NumFacesPerElement_;
      else if (what == "ElementType:") 
        GridFile >> ElementType_;
      else if (what == "NumMyElements:") 
        GridFile >> NumMyElements_;
      else if (what == "NumMyVertices:") 
        GridFile >> NumMyVertices_;
      else if (what == "NumMyBoundaryFaces:") 
        GridFile >> NumMyBoundaryFaces_;

      else if (what == "VertexCoord:")
      {
        if (NumMyVertices_ == -1)
          throw(Exception(__FILE__, __LINE__,
                          "NumMyVertices not set when reading VertexCoord"));
        if (NumDimensions_ == -1)
          throw(Exception(__FILE__, __LINE__,
                          "NumDimensions not set when reading VertexCoord"));

        // last entry is the tag
        VertexCoord_.Shape(NumMyVertices_, NumDimensions_ + 1);

        for (int i = 0 ; i < NumMyVertices_ ; ++i)
          for (int j = 0 ; j < NumDimensions_ + 1 ; ++j)
            GridFile >> VertexCoord_(i, j);

        IsVertexCoordRead = true;
      }
      else if (what == "ElementVertices:")
      {
        if (NumMyElements_ == -1)
          throw(Exception(__FILE__, __LINE__,
                          "NumMyElements not set when reading ElementVertices"));
        if (NumVerticesPerElement_ == -1)
          throw(Exception(__FILE__, __LINE__,
                          "NumVerticesPerElement not set when reading ElementVertices"));

        // last entry is the tag
        ElementVertices_.Shape(NumMyElements_, NumVerticesPerElement_ + 1);

        for (int i = 0 ; i < NumMyElements_ ; ++i)
          for (int j = 0 ; j < NumVerticesPerElement_ + 1 ; ++j)
            GridFile >> ElementVertices_(i, j);

        IsElementVerticesRead = true;
      }

      else if (what == "FaceVertices:")
      {
        if (NumMyBoundaryFaces_ == -1)
          throw(Exception(__FILE__, __LINE__,
                          "NumMyBoundaryFaces not set when reading FaceVertices"));
        if (NumVerticesPerFace_ == -1)
          throw(Exception(__FILE__, __LINE__,
                          "NumVerticesPerFace not set when reading FaceVertices"));

        // last entry is the tag
        FaceVertices_.Shape(NumMyBoundaryFaces_, NumVerticesPerFace_ + 1);

        for (int i = 0 ; i < NumMyBoundaryFaces_ ; ++i)
        {
          for (int j = 0 ; j < NumVerticesPerFace_ + 1 ; ++j)
            GridFile >> FaceVertices_(i, j);
        }

        IsFaceVerticesRead = true;
      }
    } // till file end

    if (!IsVertexCoordRead || !IsElementVerticesRead ||
        !IsFaceVerticesRead)
      throw(Exception(__FILE__, __LINE__,
                      "One among VertexCoord, ElementVertices and",
                      "FaceVertices has not been read"));

    if (NumDimensions_ == -1)
      throw(Exception(__FILE__, __LINE__, 
                      "NumDimensions not read from file"));
    if (NumVerticesPerElement_ == -1)
      throw(Exception(__FILE__, __LINE__, 
                      "NumVerticesPerElement not read from file"));
    if (NumVerticesPerFace_ == -1)
      throw(Exception(__FILE__, __LINE__, 
                      "NumVerticesPerFace not read from file"));
    if (NumFacesPerElement_ == -1)
      throw(Exception(__FILE__, __LINE__, 
                      "NumFacesPerElement not read from file"));
    if (NumMyVertices_ == -1)
      throw(Exception(__FILE__, __LINE__, 
                      "NumMyVertices not read from file"));
    if (NumMyElements_ == -1)
      throw(Exception(__FILE__, __LINE__, 
                      "NumMyElements not read from file"));
    if (NumMyBoundaryFaces_ == -1)
      throw(Exception(__FILE__, __LINE__, 
                      "NumMyBoundaryFaces not read from file"));
    if (ElementType_ == "NOT SET")
      throw(Exception(__FILE__, __LINE__, 
                      "ElementType not read from file"));
    
    NumGlobalVertices_ = NumMyVertices_;
    NumGlobalElements_ = NumMyElements_;
    NumGlobalBoundaryFaces_ = NumMyBoundaryFaces_;

    // This maps are all trivial here
    ElementMap_ = rcp(new Epetra_Map(-1, NumMyElements(), 0, Comm));
    VertexMap_ = rcp(new Epetra_Map(-1, NumMyVertices(), 0, Comm));
    RowMap_ = rcp(new Epetra_Map(-1, NumMyVertices(), 0, Comm));

    Importer_ = rcp(new Epetra_Import(VertexMap(), RowMap()));
    
    // computes the length and other element quantities

    ElementMaxLength_.Shape(NumMyElements_, 1);
    ElementMinLength_.Shape(NumMyElements_, 1);
    ElementVolume_.Shape(NumMyElements_, 1);

    if (ElementType_ == "GALERI_TRIANGLE")
    {
      int IDs[3];
      double x[3], y[3], z[3];
      z[0] = 0.0, z[1] = 0.0, z[2] = 0.0;

      for (int ie = 0 ; ie < NumMyElements_ ; ++ie)
      {
        ElementVertices(ie, IDs);

        for (int i = 0 ; i < 3 ; ++i)
        {
          x[i] = VertexCoord_(IDs[i], 0);
          y[i] = VertexCoord_(IDs[i], 1);
        }

        double l0 = Length(x[0], y[0], 0.0, x[1], y[1], 0.0);
        double l1 = Length(x[1], y[1], 0.0, x[2], y[2], 0.0);
        double l2 = Length(x[2], y[2], 0.0, x[0], y[0], 0.0);

        double max = l0;
        if (l1 > max) max = l1;
        if (l2 > max) max = l2;

        double min = l0;
        if (l1 < min) min = l1;
        if (l2 < min) min = l2;

        ElementMaxLength_[ie] = max;
        ElementMinLength_[ie] = min;
        ElementVolume_[ie] = AreaOfTriangle(x, y, z);
      }
    }
    else if (ElementType_ == "GALERI_QUAD")
    {
      int IDs[4];
      double x[4], y[4], z[4];
      z[0] = 0.0, z[1] = 0.0, z[2] = 0.0, z[3] = 0.0;

      for (int ie = 0 ; ie < NumMyElements_ ; ++ie)
      {
        ElementVertices(ie, IDs);

        for (int i = 0 ; i < 4 ; ++i)
        {
          x[i] = VertexCoord_(IDs[i], 0);
          y[i] = VertexCoord_(IDs[i], 1);
        }

        double l0 = Length(x[0], y[0], 0.0, x[1], y[1], 0.0);
        double l1 = Length(x[1], y[1], 0.0, x[2], y[2], 0.0);
        double l2 = Length(x[2], y[2], 0.0, x[3], y[3], 0.0);
        double l3 = Length(x[3], y[3], 0.0, x[0], y[0], 0.0);

        double max = l0;
        if (l1 > max) max = l1;
        if (l2 > max) max = l2;
        if (l3 > max) max = l3;

        double min = l0;
        if (l1 < min) min = l1;
        if (l2 < min) min = l2;
        if (l3 < max) max = l3;

        ElementMaxLength_[ie] = max;
        ElementMinLength_[ie] = min;
        ElementVolume_[ie] = AreaOfQuad(x, y, z);
      }
    }
    else if (ElementType_ == "GALERI_TET")
    {
      for (int ie = 0 ; ie < NumMyElements_ ; ++ie)
      {
        ElementMaxLength_[ie] = -1.;
        ElementMinLength_[ie] = -1.;
        ElementVolume_[ie] = -1;
      }
    }
    else if (ElementType_ == "GALERI_HEX")
    {
      for (int ie = 0 ; ie < NumMyElements_ ; ++ie)
      {
        ElementMaxLength_[ie] = -1.;
        ElementMinLength_[ie] = -1.;
        ElementVolume_[ie] = -1.;
      }
    }
  }

  virtual ~FileGrid() {}

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
    if (LocalID >= NumMyVertices_)
      throw(Exception(__FILE__, __LINE__, 
                      "LocalID = " + toString(LocalID) + ", " +
                      "NumMyVertices = " + toString(NumMyVertices_)));

    for (int i = 0 ; i < NumDimensions_ ; ++i)
      coord[i] = VertexCoord_(LocalID, i);
    if (NumDimensions_ == 2)
      coord[2] = 0.0;
  }

  virtual void VertexCoord(const int Length, const int* IDs, 
                           double* x, double* y, double* z) const
  {
    for (int i = 0 ; i < Length ; ++i)
    {
      if (IDs[i] >= NumMyVertices_)
        throw(Exception(__FILE__, __LINE__, 
                        "LocalID = " + toString(IDs[i]) + ", " +
                        "NumMyVertices = " + toString(NumMyVertices_)));

      x[i] = VertexCoord_(IDs[i], 0);
      y[i] = VertexCoord_(IDs[i], 1);
      if (NumDimensions_ == 3)
        z[i] = VertexCoord_(IDs[i], 2);
      else
        z[i] = 0.0;
    }
  }

  virtual void ElementVertices(const int LocalID, int* elements) const
  {
    for (int i = 0 ; i < NumVerticesPerElement_ ; ++i)
      elements[i] = ElementVertices_(LocalID, i);
  }

  virtual double ElementMinLength(const int LocalElement) const
  {
    return(ElementMinLength_[LocalElement]);
  }

  virtual double ElementMaxLength(const int LocalElement) const
  {
    return(ElementMaxLength_[LocalElement]);
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

  virtual int ElementTag(const int LocalID) const
  {
    return(ElementVertices_(LocalID, NumVerticesPerElement_));
  }

  virtual int VertexTag(const int LocalID) const
  {
    return((int)(VertexCoord_(LocalID, NumDimensions_)));
  }

  virtual double ElementVolume() const
  {
    throw(Exception(__FILE__, __LINE__, "Not implemented"));
  }

  virtual void FaceVertices(const int LocalFace, int& tag, int* IDs) const
  {
    for (int i = 0 ; i < NumVerticesPerFace_ ; ++i)
      IDs[i] = FaceVertices_(LocalFace, i);
    tag = FaceVertices_(LocalFace, NumVerticesPerFace_);
  }

  inline int FacePatch(const int LocalFace) const
  {
    return(FaceVertices_(LocalFace, NumVerticesPerFace_));
  }

  virtual double ElementVolume(const int LocalElement) const
  {
    return(ElementVolume_(LocalElement));
  }

  virtual double FaceArea(const int LocalFace) const
  {
    throw(Exception(__FILE__, __LINE__, "Not implemented"));
  }

  virtual double MyVolume() const
  {
    throw(Exception(__FILE__, __LINE__, "Not implemented"));
  }

  virtual double GlobalVolume() const
  {
    throw(Exception(__FILE__, __LINE__, "Not implemented"));
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

   int NumNeighborsPerElement() const
   {
     throw(Exception(__FILE__, __LINE__, "Not implemented"));
   }

   void ElementNeighbors(int, int*) const
   {
     throw(Exception(__FILE__, __LINE__, "Not implemented"));
   }

private:
  const Epetra_Comm& Comm_;
  string FileName_;

  int NumDimensions_;
  int NumVerticesPerElement_;
  int NumVerticesPerFace_;
  int NumFacesPerElement_;

  int NumMyVertices_;
  int NumGlobalVertices_;

  int NumMyElements_;
  int NumGlobalElements_;

  int NumMyBoundaryFaces_;
  int NumGlobalBoundaryFaces_;

  string ElementType_;

  Epetra_SerialDenseMatrix VertexCoord_;
  Epetra_IntSerialDenseMatrix ElementVertices_;
  Epetra_IntSerialDenseMatrix FaceVertices_;
  Epetra_SerialDenseVector ElementMaxLength_;
  Epetra_SerialDenseVector ElementMinLength_;
  Epetra_SerialDenseVector ElementVolume_;

  RefCountPtr<Epetra_Map> VertexMap_;
  RefCountPtr<Epetra_Map> ElementMap_;
  RefCountPtr<Epetra_Map> RowMap_;
  RefCountPtr<Epetra_Import> Importer_;
}; // class FileGrid

} // namespace FiniteElements
} // namespace Galeri
#endif
