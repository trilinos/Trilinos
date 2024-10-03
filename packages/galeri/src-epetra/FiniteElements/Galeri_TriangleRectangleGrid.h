// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_TRIANGLERECTANGLEGRID_H
#define GALERI_TRIANGLERECTANGLEGRID_H

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_DistObject.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_IntSerialDenseMatrix.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Galeri_Workspace.h"
#include "Galeri_AbstractGrid.h"
#include <vector>
#include <algorithm>
#include <limits>

using namespace Teuchos;

namespace Galeri {
namespace FiniteElements {

/*!
 * \class TriangleRectangleGrid
 *
 * \brief Creates a grid composed by triangles, the domain is a rectangle.
 *
 * This class defined, on-the-fly, the triangulation of a 2D rectangular
 * domain. The elements are all triangles. For parallel run, the rectangle
 * is subdivided along the X- and Y-axis, as specified by the user.
 *
 * \author Marzio Sala, SNL 9214.
 *
 * \date Last updated on 03-Apr-05.
 */
class TriangleRectangleGrid : public AbstractGrid
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
  TriangleRectangleGrid(const Epetra_Comm& Comm, const int nx, const int ny, 
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
    if (lx <= 0.0 || ly <= 0.0 )
    {
      cerr << "Invalid length, lx = " << lx << ", ly = " << ly << endl;
      cerr << "File " << __FILE__ << ", line " << __LINE__ << endl;
      throw(-1);
    }

    if (mx * my != Comm.NumProc())
    {
      cerr << "Incorrect processor subdivision, mx = " << mx
           << ", my = " << my << endl;
      cerr << "(file " << __FILE__ << ", line " << __LINE__ << endl;
      throw(-1);
    }

    int px, py;
    GetProcessorXY(px, py);

    NumGlobalElements_ = 2 * nx_ * ny_;
    NumGlobalVertices_ = (nx_ + 1) * (ny_ + 1);

    NumMyElementsX_ = nx_ / mx_;
    NumMyElementsY_ = ny_ / my_;
    if (px == mx_ - 1) NumMyElementsX_ += nx_ % mx_;
    if (py == my_ - 1) NumMyElementsY_ += ny_ % my_;
    NumMyElements_ = 2 * NumMyElementsX_ * NumMyElementsY_;

    NumMyVerticesX_ = NumMyElementsX_ + 1;
    NumMyVerticesY_ = NumMyElementsY_ + 1;
    NumMyVertices_ = NumMyVerticesX_ * NumMyVerticesY_;

    deltax_ = lx_ / nx_;
    deltay_ = ly_ / ny_;

    CreateElementMap();
    CreateVertexMap();
    CreateBoundaryFaces();
    CreateRowMap();

    Importer_ = rcp(new Epetra_Import(VertexMap(), RowMap()));
  }

  virtual ~TriangleRectangleGrid() {}

  virtual int NumDimensions() const
  {
    return(2);
  }

  virtual int NumVerticesPerElement() const
  {
    return(3);
  }

  virtual int NumFacesPerElement() const
  {
    return(3);
  }

  virtual int NumVerticesPerFace() const
  {
    return(2);
  }

  //! Returns \c GALERI_TRIANGLE
  virtual std::string ElementType() const
  {
    return("GALERI_TRIANGLE");
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
    long long LLGlobalID = VertexMap_->GID64(LocalID);

    if(LLGlobalID > std::numeric_limits<int>::max())
    {
      cerr << "Vertex ID out of int bound" << endl;
      cerr << "File " << __FILE__ << ", line " << __LINE__ << endl;
      throw(-1);
    }

    int GlobalID = (int) LLGlobalID;

    int ix, iy;
    GetVertexXY(GlobalID, ix, iy);

    coord[0] = DeltaX() * ix;
    coord[1] = DeltaY() * iy;
    coord[2] = 0.0;
  }

  virtual void VertexCoord(const int Length, const int* IDs, 
                           double* x, double* y, double* z) const
  {
    for (int i = 0 ; i < Length ; ++i)
    {
      int ID;
      long long LLID = VertexMap_->GID64(IDs[i]);

      if(LLID > std::numeric_limits<int>::max())
      {
        cerr << "LLID ID out of int bound" << endl;
        cerr << "File " << __FILE__ << ", line " << __LINE__ << endl;
        throw(-1);
      }

      ID = (int) LLID;

      int ix, iy;
      GetVertexXY(ID, ix, iy);

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
      return DeltaX();
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
    for (int i = 0 ; i < 2 ; ++i)
      IDs[i] = BF_(LocalFace, i);
    tag = BF_(LocalFace, 2);
  }

  inline int FacePatch(const int LocalFace) const
  {
    return(BF_(LocalFace, 2));
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
    int patch = BF_(LocalFace, 2);

    if (patch == GALERI_LEFT || patch == GALERI_RIGHT)
      return(DeltaY());
    else
      return(DeltaX());
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
     return(3);
   }

   void ElementNeighbors(int, int*) const
   {
     cerr << "ElementNeighbors() not yet implemented" << endl;
     throw(-1);
   }

private:

  inline void GetVertexXY(const int& GlobalID, int& ix, int& iy) const
  {
    // FIXME: add NumGlobalVerticesXY();
    iy = GlobalID / NumGlobalVerticesX();
    ix = GlobalID % NumGlobalVerticesX();
  }

  inline void GetElementXY(const int& GlobalID, int& ix, int& iy) const
  {
    iy = GlobalID / NumGlobalElementsX();
    ix = GlobalID % NumGlobalElementsX();
  }

  inline void GetLocalElementXY(const int& LocalID, int& ix, int& iy) const
  {
    iy = LocalID / NumMyElementsX();
    ix = LocalID % NumMyElementsX();
  }

  inline void GetProcessorXY(int& ix, int& iy) const
  {
    iy = Comm().MyPID() / mx_;
    ix = Comm().MyPID() % mx_;
  }

  inline void IL_ElementVertices(const int LocalID, int* elements,
                                 const bool ReturnGlobal = false) const
  {
    int ix, iy;
    GetLocalElementXY(LocalID / 2, ix, iy);

    if (LocalID % 2 == 0)
    {
      elements[0] = ix + iy * NumMyVerticesX();
      elements[1] = elements[0] + 1;
      elements[2] = elements[0] + NumMyVerticesX();
    }
    else
    {
      elements[0] = ix + iy * NumMyVerticesX() + 1;
      elements[1] = elements[0] + NumMyVerticesX();
      elements[2] = elements[1] - 1;
    }
  }

  void CreateElementMap()
  {
    ElementMap_ = rcp(new Epetra_Map(-1, NumMyElements(), 0, Comm()));
    return;
  }

  void CreateBoundaryFaces()
  {

    /* I decompose the square in the following way:

                  ML_TOP  
             +--------------+
	     |              |
   ML_LEFT   |              |  ML_RIGHT 
	     |              |
	     +--------------+
	         ML_BOTTOM
    */

    int px, py;
    GetProcessorXY(px, py);

    NumMyBoundaryFaces_ = 0;
    if (px == 0)
      NumMyBoundaryFaces_ += NumMyElementsY();
    if (px == mx_ - 1)
      NumMyBoundaryFaces_ += NumMyElementsY();
    if (py == 0)
      NumMyBoundaryFaces_ += NumMyElementsX();
    if (py == my_ - 1)
      NumMyBoundaryFaces_ += NumMyElementsX();

    BF_.Shape(NumMyBoundaryFaces(), 3);

    int count = 0;

    int nx = NumMyVerticesX();
    int ny = NumMyVerticesY();

    // GALERI_BOTTOM
    if (py == 0)
    {
      for (int ix = 0 ; ix < NumMyElementsX() ; ix++) 
      {
        BF_(count, 0) = ix;
        BF_(count, 1) = ix + 1;
        ++count;
      }
    }
  
    // GALERI_RIGHT
    if (px == mx_ - 1)
    {
      for (int iy = 0 ; iy < NumMyElementsY() ; iy++) 
      {
        BF_(count, 0) = nx * (iy + 1) -1;
        BF_(count, 1) = nx * (iy + 1) + nx - 1;
        ++count;
      }
    }
    
    // GALERI_TOP
    if (py == my_ - 1)
    {
      for (int ix = 0 ; ix < NumMyElementsX() ; ix++)
      {
        BF_(count, 0) = nx * (ny - 1) + ix;
        BF_(count, 1) = nx * (ny - 1) + ix + 1;
        ++count;
      }
    }

    // GALERI_LEFT
    if (px == 0)
    {
      for (int iy = 0 ; iy < NumMyElementsY() ; iy++) 
      {
        BF_(count, 0) = iy * nx;
        BF_(count, 1) = iy + 1 * nx;
        ++count;
      }
    }

    if (count != NumMyBoundaryFaces())
    {
      cerr << "Internal error, count != NumMyBoundaryFaces(), "
           << count << " vs. " << NumMyBoundaryFaces() << endl;
      cerr << "File " << __FILE__<< ", line " << __LINE__ << endl;
      throw(-1);
    }

    return;
  }

  void CreateVertexMap()
  {
#ifdef EPETRA_NO_32BIT_GLOBAL_INDICES
    std::vector<long long> itmp(NumMyVertices());
#else
    std::vector<int> itmp(NumMyVertices());
#endif

    int count = 0;
    int px, py;
    GetProcessorXY(px, py);
    int startx = px * (NumGlobalElementsX() / NumDomainsX());
    int starty = py * (NumGlobalElementsY() / NumDomainsY());
    int endx = startx + NumMyVerticesX();
    int endy = starty + NumMyVerticesY();

    for (int iy = starty ; iy < endy ; ++iy)
    {
      for (int ix = startx ; ix < endx ; ++ix)
      {
        itmp[count++] = ix + iy * NumGlobalVerticesX();
      }
    }
    assert (count == NumMyVertices());

    VertexMap_ = rcp(new Epetra_Map(-1, NumMyVertices(), &itmp[0], 0, Comm()));

    return;
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
#ifdef EPETRA_NO_32BIT_GLOBAL_INDICES
    std::vector<long long> itmp(size);
#else
    std::vector<int> itmp(size);
#endif
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
  RefCountPtr<Epetra_Map> RowMap_;
  RefCountPtr<Epetra_Import> Importer_;

  Epetra_IntSerialDenseMatrix BF_;
};

} // namespace FiniteElements
} // namespace Galeri
#endif
