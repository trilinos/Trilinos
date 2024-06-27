// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_TETCUBEGRID_H
#define GALERI_TETCUBEGRID_H

/*!
 * \file Galeri_TetCubeGrid.h
 */
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_IntSerialDenseMatrix.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Galeri_Workspace.h"
#include "Galeri_AbstractGrid.h"
#include <vector>

using namespace Teuchos;

namespace Galeri {
namespace FiniteElements {

/*!
 * \class TetCubeGrid
 *
 * \brief Creates a grid with tetrahedral elements in a cube.
 *
 * \author Marzio Sala, SNL 9214.
 *
 * \date Last updated on 31-Mar-05.
 */

class TetCubeGrid : public AbstractGrid
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
   * \param mz - (In) Number of subdomains along the Z-axis.
   *
   * \param lx - (In) Length of the cube along the X-axis.
   *
   * \param ly - (In) Length of the cube along the Y-axis.
   *
   * \param lz - (In) Length of the cube along the Z-axis.
   *
   * \note The total number of processors must equal mx * my * mz.
   */
  TetCubeGrid(const Epetra_Comm& Comm, const int nx, const int ny, const int nz,
              const int mx, const int my, const int mz,
              const double lx = 1.0, const double ly = 1.0, const double lz = 1.0) :
    Comm_(Comm),
    nx_(nx),
    ny_(ny),
    nz_(nz),
    lx_(lx),
    ly_(ly),
    lz_(lz),
    mx_(mx),
    my_(my),
    mz_(mz)
  {
    // check input
    if (lx <= 0.0 || ly <= 0.0 || lz <= 0.0)
    {
      cerr << "Invalid length, lx = " << lx << ", ly = " << ly
           << ", lz = " << lz << endl;
      throw(-1);
    }

    if (mx * my * mz != Comm.NumProc())
    {
      cerr << "Incorrect processor subdivision, mx = " << mx 
           << "my = " << my << ", mz = " << mz << endl;
      throw(-1);
    }

    int px, py, pz;
    GetProcessorXYZ(px, py, pz);

    NumGlobalElements_ = 6* nx_ * ny_ * nz_;
    NumGlobalVertices_ = (nx_ + 1) * (ny_ + 1) * (nz + 1);

    NumMyElementsX_ = nx_ / mx_;
    NumMyElementsY_ = ny_ / my_;
    NumMyElementsZ_ = nz_ / mz_;
    if (px == mx_ - 1) NumMyElementsX_ += nx_ % mx_;
    if (py == my_ - 1) NumMyElementsY_ += ny_ % my_;
    if (pz == mz_ - 1) NumMyElementsZ_ += nz_ % mz_;
    NumMyElements_ = 6 * NumMyElementsX_ * NumMyElementsY_ * NumMyElementsZ_;

    NumMyVerticesX_ = NumMyElementsX_ + 1;
    NumMyVerticesY_ = NumMyElementsY_ + 1;
    NumMyVerticesZ_ = NumMyElementsZ_ + 1;
    NumMyVertices_ = NumMyVerticesX_ * NumMyVerticesY_ * NumMyVerticesZ_;

    deltax_ = lx_ / nx_;
    deltay_ = ly_ / ny_;
    deltaz_ = lz_ / nz_;

    CreateT();
    CreateBF();

    CreateVertexMap();
    CreateRowMap();

    Importer_ = rcp(new Epetra_Import(VertexMap(), RowMap()));
  }

  //! Destructor
  virtual ~TetCubeGrid() {}

  virtual int NumDimensions() const
  {
    return(3);
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
    return(3);
  }

  //! Returns \c GALERI_TET
  virtual std::string ElementType() const
  {
    return("GALERI_TET");
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
    int GlobalID = VertexMap_->GID(LocalID);

    int ix, iy, iz;
    GetVertexXYZ(GlobalID, ix, iy, iz);

    coord[0] = DeltaX() * ix;
    coord[1] = DeltaY() * iy;
    coord[2] = DeltaZ() * iz;
  }

  virtual void VertexCoord(const int Length, const int* IDs, 
                           double* x, double* y, double* z) const
  {
    for (int i = 0 ; i < Length ; ++i)
    {
      int ID = VertexMap_->GID(IDs[i]);

      int ix, iy, iz;
      GetVertexXYZ(ID, ix, iy, iz);

      x[i] = DeltaX() * ix;
      y[i] = DeltaY() * iy;
      z[i] = DeltaZ() * iz;
    }
  }

  virtual void ElementVertices(const int LocalID, int* elements) const
  {
    for (int i = 0 ; i < 4 ; ++i)
      elements[i] = T_(LocalID, i);
  }

  virtual double ElementMinLength(const int LocalElement) const
  {
    double min = DeltaX();
    if (DeltaY() < min) min = DeltaY();
    if (DeltaZ() < min) min = DeltaZ();
    return(min);
  }

  virtual double ElementMaxLength(const int LocalElement) const
  {
    return(sqrt(DeltaX() * DeltaX() + DeltaY()*DeltaY()) + DeltaZ() * DeltaZ());
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
    return((1.0 / 6) * DeltaX() * DeltaY() * DeltaZ());
  }

  virtual void FaceVertices(const int LocalFace, int& tag, int* IDs) const
  {
    for (int i = 0 ; i < 4 ; ++i)
      IDs[i] = BF_(LocalFace, i);
    tag = BF_(LocalFace, 3);
  }

  inline int FacePatch(const int LocalFace) const
  {
    return(BF_(LocalFace, 3));
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

  inline int NumMyElementsX() const
  {
    return(NumMyElementsX_);
  }

  inline int NumMyElementsY() const
  {
    return(NumMyElementsY_);
  }

  inline int NumMyElementsZ() const
  {
    return(NumMyElementsZ_);
  }

  inline int NumMyVerticesX() const
  {
    return(NumMyVerticesX_);
  }

  inline int NumMyVerticesY() const
  {
    return(NumMyVerticesY_);
  }

  inline int NumMyVerticesZ() const
  {
    return(NumMyVerticesZ_);
  }

  inline int NumGlobalElementsX() const
  {
    return(nx_);
  }

  inline int NumGlobalElementsY() const
  {
    return(ny_);
  }

  inline int NumGlobalElementsZ() const
  {
    return(nz_);
  }

  inline int NumGlobalVerticesX() const
  {
    return(nx_ + 1);
  }

  inline int NumGlobalVerticesY() const
  {
    return(ny_ + 1);
  }

  inline int NumGlobalVerticesZ() const
  {
    return(nz_ + 1);
  }

  inline double LengthX() const
  {
    return(lx_);
  }

  inline double LengthY() const
  {
    return(ly_);
  }

  inline double LengthZ() const
  {
    return(lz_);
  }

  inline double DeltaX() const
  {
    return(deltax_);
  }

  inline double DeltaY() const
  {
    return(deltay_);
  }
  
  inline double DeltaZ() const
  {
    return(deltaz_);
  }
  
  virtual double ElementVolume(const int LocalElement) const
  {
    return(DeltaX() * DeltaY() * DeltaZ());
  }

  virtual double FaceArea(const int LocalFace) const
  {
    int patch = BF_(LocalFace, 4);

    if (patch == GALERI_BOTTOM || patch == GALERI_TOP)
      return(0.5 * DeltaX() * DeltaY());
    else if (patch == GALERI_LEFT || patch == GALERI_RIGHT)
      return(0.5 * DeltaY() * DeltaZ());
    else
      return(0.5 * DeltaX() * DeltaZ());
  }

  virtual double MyVolume() const
  {
    return((1.0 / 6) * LengthX() * LengthY() * LengthZ());
  }

  virtual double GlobalVolume() const
  {
    return(LengthX() * LengthY() * LengthZ());
  }

  int NumDomainsX() const
  {
    return(mx_);
  }

  int NumDomainsY() const
  {
    return(my_);
  }

  int NumDomainsZ() const
  {
    return(mz_);
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

  void CreateT()
  {
    T_.Shape(NumMyElements(), 4);

    int nx = NumMyVerticesX();
    int ny = NumMyVerticesY();
    int nz = NumMyVerticesZ();

    int count = 0;
    for (int iy = 1 ; iy < ny ; iy++) 
    {
      for (int ix = 1 ; ix < nx ; ix++) 
      {
        for (int iz = 1 ; iz < nz ; iz++) 
        {
          int i1 = ( ix + (iy-1)*nx  - 1 ) + nx*ny *(iz - 1);
          int i2 = ( ix + (iy-1)*nx      ) + nx*ny *(iz - 1);
          int i3 = ( ix + iy * nx        ) + nx*ny *(iz - 1);
          int i4 = ( ix + iy * nx    - 1 ) + nx*ny *(iz - 1);
          int i5 = i1 + nx*ny;
          int i6 = i2 + nx*ny;
          int i7 = i3 + nx*ny;
          int i8 = i4 + nx*ny;

          T_(count,0) = i1;
          T_(count,1) = i2;
          T_(count,2) = i4;
          T_(count,3) = i5; count++;

          T_(count,0) = i2;
          T_(count,1) = i4;
          T_(count,2) = i5;
          T_(count,3) = i6; count++;

          T_(count,0) = i4;
          T_(count,1) = i6;
          T_(count,2) = i8;
          T_(count,3) = i5; count++;

          T_(count,0) = i4;
          T_(count,1) = i3;
          T_(count,2) = i8;
          T_(count,3) = i6; count++;

          T_(count,0) = i3;
          T_(count,1) = i7;
          T_(count,2) = i8;
          T_(count,3) = i6; count++;

          T_(count,0) = i2;
          T_(count,1) = i3;
          T_(count,2) = i4;
          T_(count,3) = i6; count++;
        }
      }
    }

    assert (NumMyElements() == count);
  }

  void CreateBF()
  {

    /* I decompose the cube in the following way:
     
                 +--------------+     
                /              /|
               / upper_side   / |      
              /              /  |
 face_5 ---> +--------------+ <-+--- face_3
	     |              |   +
             |              |  /
             |  face_2      | / 
	     |              |/
	     +--------------+
                lower_side
   z
   ^  + y 		
   | /
   |/
   +----> x

    */

    int px, py, pz;
    GetProcessorXYZ(px, py, pz);

    NumMyBoundaryFaces_ = 0;
    if (px == 0)
      NumMyBoundaryFaces_ += 2 * NumMyElementsY() * NumMyElementsZ();
    if (px == mx_ - 1)
      NumMyBoundaryFaces_ += 2 * NumMyElementsY() * NumMyElementsZ();
    if (py == 0)
      NumMyBoundaryFaces_ += 2 * NumMyElementsX() * NumMyElementsZ();
    if (py == my_ - 1)
      NumMyBoundaryFaces_ += 2 * NumMyElementsX() * NumMyElementsZ();
    if (pz == 0)
      NumMyBoundaryFaces_ += 2 * NumMyElementsX() * NumMyElementsY();
    if (pz == mz_ - 1)
      NumMyBoundaryFaces_ += 2 * NumMyElementsX() * NumMyElementsY();

    int nx = NumMyVerticesX();
    int ny = NumMyVerticesY();
    int nz = NumMyVerticesZ();
    
    BF_.Shape(NumMyBoundaryFaces(), 4);

    int count = 0;

    int offset;
    if (pz == 0)
    {
      for (int iy = 0 ; iy < NumMyElementsY() ; iy++) 
      {
        offset = iy * NumMyVerticesX();
        for (int ix = 0 ; ix < NumMyElementsX() ; ix++) 
        {
          BF_(count,0) = offset + ix;
          BF_(count,1) = offset + ix + 1;
          BF_(count,2) = offset + ix + NumMyVerticesX();
          BF_(count,3) = GALERI_BOTTOM;
          ++count;

          BF_(count,0) = offset + ix + 1;
          BF_(count,1) = offset + ix + NumMyVerticesX() + 1;
          BF_(count,2) = offset + ix + NumMyVerticesX();
          BF_(count,3) = GALERI_BOTTOM;
          ++count;
        }
      }
    }

    if (pz == mz_ - 1)
    {
      for (int iy = 0 ; iy < NumMyElementsY() ; iy++) 
      {
        for (int ix = 0 ; ix < NumMyElementsX() ; ix++) 
        {
          offset = NumMyVerticesX() * NumMyVerticesY() *
            (NumMyVerticesZ() - 1) + iy * NumMyVerticesX();

          BF_(count,0) = offset + ix;
          BF_(count,1) = offset + ix + 1;
          BF_(count,2) = offset + ix + NumMyVerticesX();
          BF_(count,3) = GALERI_TOP;
          ++count;

          BF_(count,0) = offset + ix + 1;
          BF_(count,1) = offset + ix + 1 + NumMyVerticesX();
          BF_(count,2) = offset + ix + NumMyVerticesX();
          BF_(count,3) = GALERI_TOP;
          ++count;
        }
      }
    }

    if (py == 0)
    {
      for (int iz = 0 ; iz < NumMyElementsZ() ; iz++) 
      {
        offset = nx*ny*iz;
        for (int ix = 0 ; ix < NumMyElementsX() ; ix++)
        {
          BF_(count,0) = ix + offset;
          BF_(count,1) = ix + offset + 1;
          BF_(count,2) = ix + offset + nx*ny;
          BF_(count,3) = GALERI_FRONT;
          ++count;

          BF_(count,0) = ix + offset + 1;
          BF_(count,1) = ix + offset + nx*ny + 1;
          BF_(count,2) = ix + offset + nx*ny;
          BF_(count,3) = GALERI_FRONT;
          ++count;
        }
      }
    }
  
    if (py == my_ - 1)
    {
      for (int iz = 0 ; iz < nz-1 ; iz++)
      {
        offset = nx*ny*iz + nx*(ny-1) ;
        for (int ix = 0 ; ix < nx-1 ; ix++)
        {
          BF_(count,0) = ix + offset;
          BF_(count,1) = ix + offset + 1;
          BF_(count,2) = ix + offset + nx * ny;
          BF_(count,3) = GALERI_REAR;
          ++count;

          BF_(count,0) = ix + offset + 1;
          BF_(count,1) = ix + offset + nx * ny + 1;
          BF_(count,2) = ix + offset + nx * ny;
          BF_(count,3) = GALERI_REAR;
          ++count;

        }
      }
    }

    if (px == 0)
    {
      for (int iz = 0 ; iz < nz-1 ; iz++) 
      {
        for (int iy = 0 ; iy < ny-1 ; iy++) 
        {
          offset = nx*ny*iz + (iy+1)*nx - 1;

          BF_(count,0) = offset;
          BF_(count,1) = offset + nx;
          BF_(count,2) = offset + nx*ny;
          BF_(count,3) = GALERI_RIGHT;
          ++count;

          BF_(count,0) = offset + nx;
          BF_(count,1) = offset + nx + nx*ny;
          BF_(count,2) = offset + nx*ny;
          BF_(count,3) = GALERI_RIGHT;
          ++count;
        }
      }
    }

    if (px == mx_ - 1)
    {
      for (int iz=0 ; iz<nz-1 ; iz++ ) 
      {
        for (int iy=0 ; iy<ny-1 ; iy++ ) 
        {
          offset = nx*ny*iz + iy*nx;

          BF_(count,0) = offset;
          BF_(count,1) = offset + nx;
          BF_(count,2) = offset + nx*ny;
          BF_(count,3) = GALERI_LEFT;
          ++count;

          BF_(count,0) = offset + nx;
          BF_(count,1) = offset + nx + nx*ny;
          BF_(count,2) = offset + nx*ny;
          BF_(count,3) = GALERI_LEFT;
          ++count;
        }
      }
    }

    assert (NumMyBoundaryFaces() == count);
  }

  inline void GetVertexXYZ(const int& GlobalID, int& ix, int& iy, int& iz) const
  {
    // FIXME: add NumGlobalVerticesXY();
    int ixy = GlobalID % (NumGlobalVerticesX() * NumGlobalVerticesY());
    iz = GlobalID / (NumGlobalVerticesX() * NumGlobalVerticesY());
    iy = ixy / NumGlobalVerticesX();
    ix = ixy % NumGlobalVerticesX();
  }

  inline void GetElementXYZ(const int& GlobalID, int& ix, int& iy, int& iz) const
  {
    int ixy = GlobalID % (NumGlobalElementsX() * NumGlobalElementsY());
    iz = GlobalID / (NumGlobalElementsX() * NumGlobalElementsY());
    iy = ixy / NumGlobalElementsX();
    ix = ixy % NumGlobalElementsX();
  }

  inline void GetLocalElementXYZ(const int& MyID, int& ix, int& iy, int& iz) const
  {
    int ixy = MyID % (NumMyElementsX() * NumMyElementsY());
    iz = MyID / (NumMyElementsX() * NumMyElementsY());
    iy = ixy / NumMyElementsX();
    ix = ixy % NumMyElementsX();
  }

  inline void GetProcessorXYZ(int& ix, int& iy, int& iz) const
  {
    int ixy = Comm().MyPID() % (mx_ * my_);
    iz = Comm().MyPID() / (mx_ * my_);
    iy = ixy / mx_;
    ix = ixy % mx_;
  }

  inline void IL_ElementVertices(const int LocalID, int* elements,
                                 const bool ReturnGlobal = false) const
  {
    int ix, iy, iz;
    GetLocalElementXYZ(LocalID, ix, iy, iz);
    int nx = NumMyVerticesX();
    int ny = NumMyVerticesY();
    //int nz = NumMyVerticesZ();
    elements[0] = (ix + iy * nx          ) + nx * ny * iz;
    elements[1] = (ix + iy * nx + 1      ) + nx * ny * iz;
    elements[2] = (ix + (iy + 1) * nx + 1) + nx * ny * iz;
    elements[3] = (ix + (iy + 1) * nx    ) + nx * ny * iz;
    elements[4] = elements[0] + nx * ny;
    elements[5] = elements[1] + nx * ny;
    elements[6] = elements[2] + nx * ny;
    elements[7] = elements[3] + nx * ny;
  }

  void CreateElementMap()
  {
    ElementMap_ = rcp(new Epetra_Map(-1, NumGlobalElements(), 0, Comm()));
    return;
  }

  void CreateBoundaryFaces()
  {
    /* I decompose the cube in the following way:
     
                 +--------------+     
                /              /|
               /   ML_TOP     / |      
              /              /  |
 ML_LEFT --> +--------------+ <-+--- ML_RIGHT
	     |              |   +
             |              |  /
             |   ML_FRONT   | / 
	     |              |/
	     +--------------+
                ML_BOTTOM 
   z
   ^  + y 		
   | /
   |/
   +----> x

    */

    int px, py, pz;
    GetProcessorXYZ(px, py, pz);

    NumMyBoundaryFaces_ = 0;
    if (px == 0)
      NumMyBoundaryFaces_ += NumMyElementsY() * NumMyElementsZ();
    if (px == mx_ - 1)
      NumMyBoundaryFaces_ += NumMyElementsY() * NumMyElementsZ();
    if (py == 0)
      NumMyBoundaryFaces_ += NumMyElementsX() * NumMyElementsZ();
    if (py == my_ - 1)
      NumMyBoundaryFaces_ += NumMyElementsX() * NumMyElementsZ();
    if (pz == 0)
      NumMyBoundaryFaces_ += NumMyElementsX() * NumMyElementsY();
    if (pz == mz_ - 1)
      NumMyBoundaryFaces_ += NumMyElementsX() * NumMyElementsY();

    BF_.Shape(NumMyBoundaryFaces(), 5);

    int count = 0, offset = 0;
    int nx = NumMyVerticesX();
    int ny = NumMyVerticesY();
    int nz = NumMyVerticesZ();

    // lower side
    
    if (pz == 0)
    {
      for (int iy = 0 ; iy < NumMyElementsY() ; iy++) 
      {
        for (int ix = 0 ; ix < NumMyElementsX() ; ix++) 
        {
          BF_(count,0) = ix + iy*nx;
          BF_(count,1) = ix + iy*nx+1;
          BF_(count,2) = ix + (iy+1)*nx + 1;
          BF_(count,3) = ix + (iy+1)*nx;
          BF_(count,4) = GALERI_BOTTOM;
          ++count;
        }
      }
    }

    // upper_side 

    if (pz == mz_ - 1)
    {
      offset = nx*ny*(nz-1);
      for (int iy = 0 ; iy < NumMyElementsY() ; iy++) 
      {
        for (int ix = 0 ; ix < NumMyElementsX() ; ix++) 
        {
          BF_(count,0) = offset + ix + iy*nx;
          BF_(count,1) = offset + ix + iy*nx+1;
          BF_(count,2) = offset + ix + (iy+1)*nx + 1;
          BF_(count,3) = offset + ix + (iy+1)*nx;
          BF_(count,4) = GALERI_TOP;
          ++count;
        }
      }
    }

    // face_2 -- front

    if (py == 0)
    {
      for (int iz = 0 ; iz < NumMyElementsZ() ; iz++) 
      {
        offset = nx*ny*iz;
        for (int ix = 0 ; ix < NumMyElementsX() ; ix++) 
        {
          BF_(count,0) = ix + offset;
          BF_(count,1) = ix + offset + 1;
          BF_(count,2) = ix + offset + nx*ny + 1;
          BF_(count,3) = ix + offset + nx*ny;
          BF_(count,4) = GALERI_FRONT;
          ++count;
        }
      }
    }
  
    // face_3 -- right

    if (px == mx_ - 1)
    {
      for (int iz = 0 ; iz < NumMyElementsZ() ; iz++) 
      {
        for (int iy = 0 ; iy < NumMyElementsY() ; iy++) 
        {
          offset = nx*ny*iz + (iy+1)*nx - 1;

          BF_(count,0) = offset;
          BF_(count,1) = offset + nx;
          BF_(count,2) = offset + nx + nx*ny;
          BF_(count,3) = offset + nx*ny;
          BF_(count,4) = GALERI_RIGHT;
          ++count;
        }
      }
    }

    // face_4 -- rear
    
    if (py == my_ - 1)
    {
      for (int iz = 0 ; iz < NumMyElementsZ() ; iz++) 
      {
        offset = nx*ny*(iz+1) - nx ;
        for (int ix = 0 ; ix < NumMyElementsX() ; ix++) 
        {

          BF_(count,0) = ix + offset;
          BF_(count,1) = ix + offset + 1;
          BF_(count,2) = ix + offset + nx*ny + 1;
          BF_(count,3) = ix + offset + nx*ny;
          BF_(count,4) = GALERI_REAR;
          ++count;
        }
      }
    }
    
    // face_5 -- left
    
    if (px == 0)
    {
      for (int iz = 0 ; iz < NumMyElementsZ() ; iz++) 
      {
        for (int iy = 0 ; iy < NumMyElementsY() ; iy++) 
        {
          offset = nx*ny*iz;
          BF_(count,0) = offset + iy*nx;
          BF_(count,1) = offset + (iy+1)*nx;
          BF_(count,2) = offset + iy*nx + nx*ny;
          BF_(count,3) = offset + (iy+1)*ny + nx*ny;
          BF_(count,4) = GALERI_LEFT;
          ++count;
        }
      }
    }

    if (count != NumMyBoundaryFaces())
    {
      cerr << "internal error, count != NumMyBoundaryFaces(), " 
           << count << " vs. " << NumMyBoundaryFaces() << endl;
      throw(-1);
    }
  }

  void CreateVertexMap()
  {
    vector<int> itmp(NumMyVertices());

    int count = 0;
    int px, py, pz;
    GetProcessorXYZ(px, py, pz);
    int startx = px * (NumGlobalElementsX() / mx_);
    int starty = py * (NumGlobalElementsY() / my_);
    int startz = pz * (NumGlobalElementsZ() / mz_);
    int endx = startx + NumMyVerticesX();
    int endy = starty + NumMyVerticesY();
    int endz = startz + NumMyVerticesZ();

    for (int iz = startz ; iz < endz ; ++iz)
    {
      for (int iy = starty ; iy < endy ; ++iy)
      {
        for (int ix = startx ; ix < endx ; ++ix)
        {
          itmp[count++] = ix + iy * NumGlobalVerticesY() + 
            iz * NumGlobalVerticesX() * NumGlobalVerticesY();
        }
      }
    }
    assert (count == NumMyVertices());

    VertexMap_ = rcp(new Epetra_Map(-1, NumMyVertices(), &itmp[0], 0, Comm()));

    return;
  }

  void CreateRowMap()
  {
    vector<int> itmp(NumMyVertices());

    int count = 0;
    int px, py, pz;
    GetProcessorXYZ(px, py, pz);
    int startx = px * (NumGlobalVerticesX() / mx_);
    int starty = py * (NumGlobalVerticesY() / my_);
    int startz = pz * (NumGlobalVerticesZ() / mz_);
    int endx = startx + NumMyVerticesX() - (px == mx_ - 1 ? 0 : 1);
    int endy = starty + NumMyVerticesY() - (py == my_ - 1 ? 0 : 1);
    int endz = startz + NumMyVerticesZ() - (pz == mz_ - 1 ? 0 : 1);

    for (int iz = startz ; iz < endz ; ++iz)
    {
      for (int iy = starty ; iy < endy ; ++iy)
      {
        for (int ix = startx ; ix < endx ; ++ix)
        {
          itmp[count++] = ix + iy * NumGlobalVerticesY() + 
            iz * NumGlobalVerticesX() * NumGlobalVerticesY();
        }
      }
    }
    assert (count <= NumMyVertices());

    RowMap_ = rcp(new Epetra_Map(-1, count, &itmp[0], 0, Comm()));

    return;
  }

  const Epetra_Comm& Comm_;

  int NumMyVertices_;
  int NumMyVerticesX_;
  int NumMyVerticesY_;
  int NumMyVerticesZ_;
  int NumGlobalVertices_;
  int NumMyElements_;
  int NumMyElementsX_;
  int NumMyElementsY_;
  int NumMyElementsZ_;
  int NumGlobalElements_;
  int NumMyBoundaryFaces_;
  int NumGlobalBoundaryFaces_;

  int nx_;
  int ny_;
  int nz_;
  double lx_;
  double ly_;
  double lz_;
  int mx_;
  int my_;
  int mz_;
  double deltax_;
  double deltay_;
  double deltaz_;

  RefCountPtr<Epetra_Map> VertexMap_;
  RefCountPtr<Epetra_Map> ElementMap_;
  RefCountPtr<Epetra_Map> RowMap_;
  RefCountPtr<Epetra_Import> Importer_;

  Epetra_IntSerialDenseMatrix T_;
  Epetra_IntSerialDenseMatrix BF_;
};

} // namespace FiniteElements
} // namespace Galeri
#endif
