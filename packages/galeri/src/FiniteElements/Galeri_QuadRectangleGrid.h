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

#ifndef GALERI_QUADRECTANGLEGRID_H
#define GALERI_QUADRECTANGLEGRID_H

/*! 
 * \file Galeri_QuadRectangleGrid.h
 */

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Galeri_AbstractGrid.h"
#include "Galeri_Workspace.h"
#include <vector>
#include <algorithm>

using namespace std;
using namespace Teuchos;

namespace Galeri {
namespace FiniteElements {

/*!
 * \class QuadRectangleGrid
 *
 * \brief Creates a grid with quadrilaterals on a rectangle.
 *
 * \author Marzio Sala, SNL 9214.
 *
 * \date Last updated on 31-Mar-05.
 */
class QuadRectangleGrid : public AbstractGrid
{

public:

  //! Constructor.
  /*! 
   * \param Comm - (In) Communicator object.
   *
   * \param nx - (In) number of elements along the X-axis.
   *
   * \param ny - (In) number of elements along the Y-axis.
   *
   * \param mx - (In) Number of subdomains along the X-axis.
   *
   * \param my - (In) Number of subdomains along the Y-axis.
   *
   * \param lx - (In) Length of the rectangle along the X-axis.
   *
   * \param ly - (In) Length of the rectangle along the Y-axis.
   *
   * \note The total number of processors must equal mx * my.
   */
  QuadRectangleGrid(const Epetra_Comm& Comm, const int nx, const int ny, 
                    const int mx, const int my, 
                    const double lx = 1.0, const double ly = 1.0) :
    Comm_(Comm),
    nx_(nx),
    ny_(ny),
    lx_(lx),
    ly_(ly),
    mx_(mx),
    my_(my)
  {
    // check input
    if (lx <= 0.0 || ly <= 0.0)
    {
      cerr << "Invalid length, lx = " << lx << ", ly = " << ly << endl;
      cerr << "File " << __FILE__ << ", line " << __LINE__ << endl;
      throw(-1);
    }

    if (mx * my != Comm.NumProc())
    {
      cerr << "Incorrect processor subdivision, mx = " << mx
           << ", my = " << my << endl;
      cerr << "File " << __FILE__ << ", line " << __LINE__ << endl;
      throw(-1);
    }

    NumGlobalElements_ = nx_ * ny_;
    NumGlobalVertices_ = (nx_ + 1) * (ny_ + 1);

    deltax_ = lx_ / nx_;
    deltay_ = ly_ / ny_;

    CreateElementMap();
    CreateVertexMap();
    CreateBoundaryFaceMap();
    CreateRowMap();

    Importer_ = rcp(new Epetra_Import(VertexMap(), RowMap()));
  }

  //! Destructor
  virtual ~QuadRectangleGrid() {}

  virtual int NumDimensions() const
  {
    return(2);
  }

  virtual int NumVerticesPerElement() const
  {
    return(4);
  }

  virtual int NumFacesPerElement() const
  {
    return(4);
  }

  virtual int NumVerticesPerFace() const
  {
    return(2);
  }

  //! Returns \c GALERI_QUAD
  virtual string ElementType() const
  {
    return("GALERI_QUAD");
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
    int GlobalID;
    GlobalID = VertexMap_->GID(LocalID);

    if (GlobalID == -1)
    {
      cerr << "Vertex ID out of bound or non local" << endl;
      cerr << "File " << __FILE__ << ", line " << __LINE__ << endl;
      throw(-1);
    }

    int ix, iy;
    ix = GlobalID % NumGlobalVerticesX();
    iy = GlobalID / NumGlobalVerticesX();

    coord[0] = DeltaX() * ix;
    coord[1] = DeltaY() * iy;
  }

  virtual void VertexCoord(const int Length, const int* IDs, 
                           double* x, double* y, double* z) const
  {
    for (int i = 0 ; i < Length ; ++i)
    {
      int ID;
      ID = VertexMap_->GID(IDs[i]);

      int ix, iy;
      ix = ID % NumGlobalVerticesX();
      iy = ID / NumGlobalVerticesX();

      x[i] = DeltaX() * ix;
      y[i] = DeltaY() * iy;
      z[i] = 0.0;
    }
  }

  virtual void ElementVertices(const int LocalID, int* elements) const
  {
    IL_ElementVertices(LocalID, elements, false);
  }

  virtual double ElementMinLength(const int LocalElement) const
  {
    if (DeltaX() < DeltaY())
      return(DeltaX());
    else 
      return(DeltaY());
  }

  virtual double ElementMaxLength(const int LocalElement) const
  {
    return(sqrt(DeltaX() * DeltaX() + DeltaY()*DeltaY()));
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

  virtual const Epetra_Map& BoundaryFaceMap() const
  {
    return(*(BoundaryFaceMap_.get()));
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
    // FIXME: the tag is not correct for parallel runs
    int face = BoundaryFaceMap_->GID(LocalFace);

    if (face < NumGlobalElementsX())
    {
      IDs[0] = face;
      IDs[1] = face + 1;
      tag = GALERI_BOTTOM;
    }
    else if (face < NumGlobalElementsX() + NumGlobalElementsY())
    {
      int mod = face - NumGlobalElementsX();
      IDs[0] = (mod + 1) * NumGlobalVerticesX() - 1;
      IDs[1] = IDs[0] + NumGlobalVerticesX();
      tag = GALERI_RIGHT;
    }
    else if (face < 2 * NumGlobalElementsX() + NumGlobalElementsY())
    {
      int mod = face - NumGlobalElementsX() - NumGlobalElementsY();
      IDs[0] = NumGlobalVerticesX() * (NumGlobalVerticesY() - 1) + mod;
      IDs[1] = IDs[0] + 1;
      tag = GALERI_TOP;
    }
    else 
    {
      int mod = face - 2 * NumGlobalElementsX() - NumGlobalElementsY();
      IDs[0] = NumGlobalVerticesX() * mod;
      IDs[1] = IDs[0] + NumGlobalVerticesX();
      tag = GALERI_LEFT;
    }

    IDs[0] = VertexMap_->LID(IDs[0]);
    IDs[1] = VertexMap_->LID(IDs[1]);

    if (IDs[0] == -1 || IDs[1] == 0)
    {
      cerr << "Internal error in FaceVertices() for face " << LocalFace << endl;
      cerr << "IDs[0] = " << IDs[0] << ", IDs[1] = " << IDs[1] << endl;
      cerr << "(file " << __FILE__ << ", line " << __LINE__ << ")" << endl;
      throw(-1);
    }
  }

  inline int FacePatch(const int LocalFace) const
  {
    int GlobalFace = BoundaryFaceMap_->GID(LocalFace);

    int Patch;

    if (GlobalFace < NumGlobalElementsX())
      Patch = GALERI_BOTTOM;
    else if (GlobalFace < NumGlobalElementsX() + NumGlobalElementsY())
      Patch = GALERI_RIGHT;
    else if (GlobalFace < 2 * NumGlobalElementsX() + NumGlobalElementsY())
      Patch = GALERI_TOP;
    else 
      Patch = GALERI_LEFT;

    return(Patch);
  }

  inline int NumMyElementsX() const
  {
    return(NumMyElementsX_);
  }

  inline int NumMyElementsY() const
  {
    return(NumMyElementsY_);
  }

  inline int NumMyVerticesX() const
  {
    return(NumMyVerticesX_);
  }

  inline int NumMyVerticesY() const
  {
    return(NumMyVerticesY_);
  }

  inline int NumGlobalElementsX() const
  {
    return(nx_);
  }

  inline int NumGlobalElementsY() const
  {
    return(ny_);
  }

  inline int NumGlobalVerticesX() const
  {
    return(nx_ + 1);
  }

  inline int NumGlobalVerticesY() const
  {
    return(ny_ + 1);
  }

  inline double LengthX() const
  {
    return(lx_);
  }

  inline double LengthY() const
  {
    return(ly_);
  }

  inline double DeltaX() const
  {
    return(deltax_);
  }

  inline double DeltaY() const
  {
    return(deltay_);
  }
  
  virtual double ElementVolume(const int LocalElement) const
  {
    return(DeltaX() * DeltaY());
  }

  virtual double FaceArea(const int LocalFace) const
  {
    int tag = FacePatch(LocalFace);
    if (tag == 0 || tag == 2)
      return (DeltaX());
    else
      return (DeltaY());
  }

  virtual double MyVolume() const
  {
    return(LengthX() * LengthY());
  }

  virtual double GlobalVolume() const
  {
    return(LengthX() * LengthY());
  }

  int NumDomainsX() const
  {
    return(mx_);
  }

  int NumDomainsY() const
  {
    return(my_);
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
     return(4);
   }

   void ElementNeighbors(int, int*) const
   {
     cerr << "ElementNeighbors() is not yet implemented" << endl;
     cerr << "File " << __FILE__ << ", line " << __LINE__ << endl;
     throw(-1);
   }

private:

  inline void IL_ElementVertices(const int LocalID, int* elements,
                                 const bool ReturnGlobal = false) const
  {
    int GlobalID;
    GlobalID = ElementMap_->GID(LocalID);

    int ix, iy;
    ix = GlobalID % NumGlobalElementsX();
    iy = GlobalID / NumGlobalElementsX();

    elements[0] = ix + iy * NumGlobalVerticesX();
    elements[1] = elements[0] + 1;  
    elements[2] = elements[1] + NumGlobalVerticesX();
    elements[3] = elements[2] - 1;

    if (ReturnGlobal == false)
      for (int i = 0 ; i < 4 ; ++i)
        elements[i] = VertexMap_->LID(elements[i]);
  }

  void CreateRowMap()
  {
    int modx = NumGlobalVerticesX() / NumDomainsX(); 
    int resx = NumGlobalVerticesX() % NumDomainsX();
    int mody = NumGlobalVerticesY() / NumDomainsY(); 
    int resy = NumGlobalVerticesY() % NumDomainsY();

    int startx, starty, endx, endy;
    int xpid = Comm().MyPID() % NumDomainsX();
    int ypid = Comm().MyPID() / NumDomainsX();

    startx = xpid * modx;
    endx   = (xpid + 1) * modx;
    if (xpid == NumDomainsX() - 1) endx += resx;

    starty = ypid * mody;
    endy   = (ypid + 1) * mody;
    if (ypid == NumDomainsY() - 1) endy += resy;

    int size = (endx - startx) * (endy - starty);

    int count = 0;
    vector<int> itmp(size);
    for (int j = starty ; j < endy ; ++j) 
    {
      for (int i = startx ; i < endx ; ++i) 
      {
        itmp[count++] = i + j * NumGlobalVerticesX();
      }
    }

    RowMap_ = rcp(new Epetra_Map(-1, count, &itmp[0], 0, Comm()));

    return;
  }

  void CreateElementMap()
  {
    int modx = NumGlobalElementsX() / NumDomainsX(); 
    int resx = NumGlobalElementsX() % NumDomainsX();
    int mody = NumGlobalElementsY() / NumDomainsY(); 
    int resy = NumGlobalElementsY() % NumDomainsY();

    int startx, starty, endx, endy;
    int xpid = Comm().MyPID() % NumDomainsX();
    int ypid = Comm().MyPID() / NumDomainsX();

    startx = xpid * modx;
    endx   = (xpid + 1) * modx;
    if (xpid == NumDomainsX() - 1) endx += resx;

    starty = ypid * mody;
    endy   = (ypid + 1) * mody;
    if (ypid == NumDomainsY() - 1) endy += resy;

    NumMyElementsX_ = endx - startx;
    NumMyElementsY_ = endy - starty;
    NumMyElements_ = (endx - startx) * (endy - starty);

    vector<int> itmp(NumMyElements());
    int count = 0;

    for (int j = starty ; j < endy ; ++j) 
    {
      for (int i = startx ; i < endx ; ++i) 
      {
        itmp[count++] = i + j * NumGlobalElementsX();
      }
    }

    ElementMap_ = rcp(new Epetra_Map(-1, count, &itmp[0], 0, Comm()));

    return;
  }

  void CreateBoundaryFaceMap()
  {
    int xpid = Comm().MyPID() % NumDomainsX();
    int ypid = Comm().MyPID() / NumDomainsX();

    NumMyBoundaryFaces_ = 0;

    if (ypid == 0) 
      NumMyBoundaryFaces_ += NumMyElementsX();
    if (xpid == NumDomainsX() - 1) 
      NumMyBoundaryFaces_ += NumMyElementsY();
    if (ypid == NumDomainsY() - 1) 
      NumMyBoundaryFaces_ += NumMyElementsX();
    if (xpid == 0) 
      NumMyBoundaryFaces_ += NumMyElementsY();

    vector<int> itmp(NumMyBoundaryFaces());
    int count = 0;

    if (ypid == 0)
    {
      int offset = xpid * (NumGlobalElementsX() / NumDomainsX());
      for (int i = 0 ; i < NumMyElementsX() ; ++i)
        itmp[count++] = offset + i;
    }

    if (xpid == NumDomainsX() - 1)
    {
      int offset = ypid * (NumGlobalElementsY() / NumDomainsY())
        + NumGlobalElementsX();
      for (int i = 0 ; i < NumMyElementsY() ; ++i)
        itmp[count++] = offset + i;
    }

    if (ypid == NumDomainsY() - 1)
    {
      int offset = NumGlobalElementsX() + NumGlobalElementsY() + 
        xpid * (NumGlobalElementsX() / NumDomainsX());

      for (int i = 0 ; i < NumMyElementsX() ; ++i)
        itmp[count++] = offset + i;
    }

    if (xpid == 0)
    {
      int offset = ypid * (NumGlobalElementsY() / NumDomainsY())
        + 2 * NumGlobalElementsX() + NumGlobalElementsY();
      for (int i = 0 ; i < NumMyElementsY() ; ++i)
        itmp[count++] = offset + i;
    }

    BoundaryFaceMap_ = rcp(new Epetra_Map(-1, count, &itmp[0], 0, Comm()));

    if (NumMyBoundaryFaces_ != count)
    {
      cerr << "Internal error, NumMyBoundaryFaces_ != count,";
      cerr << NumMyBoundaryFaces_ << " vs. " << count << endl;
      cerr << "(file " << __FILE__ << ", line " << __LINE__ << ")" << endl;
      throw(-1);
    }

    NumGlobalBoundaryFaces_ = BoundaryFaceMap_->NumGlobalElements();
    return;
  }

  void CreateVertexMap()
  {
    vector<int> tmp;
    vector<int>::iterator where;
    int Vertices[4];

    for (int i = 0 ; i < NumMyElements() ; ++i)
    {
      IL_ElementVertices(i, Vertices, true);
      for (int j = 0 ; j < NumVerticesPerElement() ; ++j)
      {
        where = find(tmp.begin(), tmp.end(), Vertices[j]);
        if (where == tmp.end())
          tmp.push_back(Vertices[j]);
      }
    }

    VertexMap_ = rcp(new Epetra_Map(-1, (int)tmp.size(), &tmp[0], 0, Comm()));
    NumMyVertices_ = VertexMap_->NumMyElements();

    return;
  }

  const Epetra_Comm& Comm_;

  int NumMyVertices_;
  int NumMyVerticesX_;
  int NumMyVerticesY_;
  int NumGlobalVertices_;
  int NumMyElements_;
  int NumMyElementsX_;
  int NumMyElementsY_;
  int NumGlobalElements_;
  int NumMyBoundaryFaces_;
  int NumGlobalBoundaryFaces_;

  int nx_;
  int ny_;
  double lx_;
  double ly_;
  int mx_;
  int my_;
  double deltax_;
  double deltay_;

  RefCountPtr<Epetra_Map> VertexMap_;
  RefCountPtr<Epetra_Map> ElementMap_;
  RefCountPtr<Epetra_Map> BoundaryFaceMap_;
  RefCountPtr<Epetra_Map> RowMap_;
  RefCountPtr<Epetra_Import> Importer_;

};

} // namespace FiniteElements
} // namespace Galeri
#endif
