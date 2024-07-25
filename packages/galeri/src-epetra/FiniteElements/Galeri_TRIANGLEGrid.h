// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_TRIANGLEGRID_H
#define GALERI_TRIANGLEGRID_H

#include "Galeri_AbstractGrid.h"
extern "C" {
#define ANSI_DECLARATORS
#define REAL double
#include "triangle.h"
}

namespace Galeri {
namespace FiniteElements {

class TRIANGLEGrid : public AbstractGrid 
{
public:

  // @{ Constructor and Destructor
  TRIANGLEGrid(const Epetra_Comm& Comm, 
               const int NumPoints, const double* x, const double* y,
               const double MaxArea) :
    Comm_(Comm)
  {

    if (Comm_.NumProc() != 1)
      throw(Exception(__FILE__, __LINE__,
                      "TRIANGLEGrid can be used w/ 1 processor only"));

    /* Define input points. */

    in_.numberofpoints = NumPoints;
    in_.numberofpointattributes = 1;
    in_.pointlist = (double *) malloc(in_.numberofpoints * 2 * sizeof(double));
    for (int i = 0 ; i < NumPoints ; ++i)
    {
      in_.pointlist[2 * i    ] = x[i];
      in_.pointlist[2 * i + 1] = y[i];
    }
    in_.pointattributelist = (double *) malloc(in_.numberofpoints *
                                               in_.numberofpointattributes *
                                               sizeof(double));
    // not so sure of the commands below
    in_.pointattributelist[0] = 0.0;
    in_.pointattributelist[1] = 1.0;
    in_.pointattributelist[2] = 11.0;
    in_.pointattributelist[3] = 10.0;
    in_.pointmarkerlist = (int *) malloc(in_.numberofpoints * sizeof(int));
    in_.pointmarkerlist[0] = 0;
    in_.pointmarkerlist[1] = 2;
    in_.pointmarkerlist[2] = 0;
    in_.pointmarkerlist[3] = 0;

    in_.numberofsegments = 0;
    in_.numberofholes = 0;
    in_.numberofregions = 1;
    in_.regionlist = (double *) malloc(in_.numberofregions * 4 * sizeof(double));
    in_.regionlist[0] = 1.0;
    in_.regionlist[1] = 1.0;
    in_.regionlist[2] = 1.0;            /* Regional attribute (for whole mesh). */
    in_.regionlist[3] = MaxArea;

    out_.pointlist = (double *) NULL;            /* Not needed if -N switch used. */
    out_.pointattributelist = (double *) NULL;
    out_.pointmarkerlist = (int *) NULL; /* Not needed if -N or -B switch used. */
    out_.trianglelist = (int *) NULL;          /* Not needed if -E switch used. */
    /* Not needed if -E switch used or number of triangle attributes is zero: */
    out_.triangleattributelist = (double *) NULL;
    out_.neighborlist = (int *) NULL;         /* Needed only if -n switch used. */
    /* Needed only if segments are output (-p or -c) and -P not used: */
    out_.segmentlist = (int *) NULL;
    /* Needed only if segments are output (-p or -c) and -P and -B not used: */
    out_.segmentmarkerlist = (int *) NULL;
    out_.edgelist = (int *) NULL;             /* Needed only if -e switch used. */
    out_.edgemarkerlist = (int *) NULL;   /* Needed if -e used and -B not used. */

    /* Triangulate the points.  Switches are chosen to read and write a  */
    /*   PSLG (p), preserve the convex hull (c), number everything from  */
    /*   zero (z), assign a regional attribute to each element (A), and  */
    /*   produce an edge list (e), a Voronoi diagram (v), and a triangle */
    /*   neighbor list (n).                                              */

    triangulate("apczAen", &in_, &out_, NULL);

    VertexMap_ = new Epetra_Map(NumMyVertices(), 0, Comm_);
    RowMap_    = new Epetra_Map(NumMyVertices(), 0, Comm_);

    min_h_.resize(NumMyElements());
    max_h_.resize(NumMyElements());

    for (int i = 0 ; i < NumMyVertices() ; ++i)
    {
      min_h_[i] = 1000000000.0;
      max_h_[i] = 0.0;
    }

    // computes the length of each element
    for (int ie = 0 ; ie < NumMyElements() ; ++ie)
    {
      int vertices[3];
      double x[3], y[3], z[3];
      double h_0, h_1, h_2;

      ElementVertices(ie, vertices);
      VertexCoord(3, vertices, x, y, z);

      h_0 = sqrt((x[0] - x[1]) * (x[0] - x[1]) + 
                 (y[0] - y[1]) * (y[0] - y[1]));
      h_1 = sqrt((x[1] - x[2]) * (x[1] - x[2]) + 
                 (y[1] - y[2]) * (y[1] - y[2]));
      h_2 = sqrt((x[0] - x[2]) * (x[0] - x[2]) + 
                 (y[0] - y[2]) * (y[0] - y[2]));

      if (h_0 < min_h_[ie]) min_h_[ie] = h_0;
      if (h_1 < min_h_[ie]) min_h_[ie] = h_1;
      if (h_2 < min_h_[ie]) min_h_[ie] = h_2;

      if (h_0 > max_h_[ie]) max_h_[ie] = h_0;
      if (h_1 > max_h_[ie]) max_h_[ie] = h_1;
      if (h_2 > max_h_[ie]) max_h_[ie] = h_2;
    }
  }

  ~TRIANGLEGrid()
  {
    if (in_.pointlist != NULL) free(in_.pointlist);
    if (in_.pointattributelist != NULL) free(in_.pointattributelist);
    if (in_.pointmarkerlist != NULL) free(in_.pointmarkerlist);
    if (in_.regionlist != NULL) free(in_.regionlist);
    if (out_.pointlist != NULL) free(out_.pointlist);
    if (out_.pointattributelist != NULL) free(out_.pointattributelist);
    if (VertexMap_) delete VertexMap_;
    if (RowMap_) delete RowMap_;
  }

  virtual int NumDimensions() const
  {
    return(2);
  }

  //! Returns the number of vertices contained in each element.
  virtual int NumVerticesPerElement() const
  {
    return(3);
  }

  //! Returns the number of faces contained in each element.
  virtual int NumFacesPerElement() const 
  {
    return(3);
  }

  //! Returns the number of vertices contained in each face.
  virtual int NumVerticesPerFace() const
  {
    return(2);
  }

  //! Returns a string containing the element type.
  virtual std::string ElementType() const
  {
    return("ML_TRIANGLE");
  }

  virtual int NumNeighborsPerElement() const
  {
    return(3);
  }

  //! Returns the number of finite elements on the calling process.
  virtual int NumMyElements() const
  {
    return(out_.numberoftriangles);
  }

  //! Returns the global number of finite elements.
  virtual int NumGlobalElements() const
  {
    return(out_.numberoftriangles);
  }

  //! Returns the number of vertices on the calling process.
  virtual int NumMyVertices() const
  {
    return(out_.numberofpoints);
  }

  //! Returns the global number of vertices.
  virtual int NumGlobalVertices() const
  {
    return(out_.numberofpoints);
  }

  //! Returns the number of boundary faces on the calling process.
  virtual int NumMyBoundaryFaces() const
  {
    return(out_.numberofsegments);
  }

  //! Returns the global number of boundary faces.
  virtual int NumGlobalBoundaryFaces() const
  {
    return(out_.numberofsegments);
  }

  //! Returns the volume of all local elements.
  virtual double MyVolume() const
  {
    throw(Exception(__FILE__, __LINE__,
                    "Feature not implemented"));
  }

  //! Returns the global volume of the grid.
  virtual double GlobalVolume() const
  {
    throw(Exception(__FILE__, __LINE__,
                    "Feature not implemented"));
  }

  //! Returns the coordinates of local vertex \c LocalVertex in vector \c coord.
  virtual void VertexCoord(const int LocalVertex, double* coord) const
  {
    coord[0] = out_.pointlist[2 * LocalVertex];
    coord[1] = out_.pointlist[2 * LocalVertex + 1];
  }

  //! Returns the coordinates of specified local vertices.
  virtual void VertexCoord(const int Length, const int* IDs, double* x, 
                           double* y, double* z) const
  {
    for (int i = 0 ; i < Length ; ++i)
    {
      int ID = IDs[i];
      x[i] = out_.pointlist[2 * ID];
      y[i] = out_.pointlist[2 * ID + 1];
      z[i] = 0.0;
    }
  }

  //! Returns the local vertex IDs of the specified local finite element.
  virtual void ElementVertices(const int LocalElement, int* elements) const
  {
    for (int j = 0; j < out_.numberofcorners; j++)
      elements[j] = out_.trianglelist[LocalElement * out_.numberofcorners + j];
  }

  //! Returns the local vertex IDs of vertices contained in the specified boundary face.
  virtual void FaceVertices(const int LocalFace, int& tag, int* IDs) const
  {
    for (int j = 0; j < 2; j++) {
      IDs[j] = out_.segmentlist[LocalFace * 2 + j];
      tag = 0;
    }
  }

  virtual void ElementNeighbors(const int LocalElement, int* elements) const
  {
    for (int j = 0 ; j < 3; ++j)
      elements[j] = out_.neighborlist[LocalElement * 3 + j];
  }

  //! Returns the patch ID of the specified face.
  virtual int FacePatch(const int LocalFace) const
  {
    return(0);
  }

  //! Returns the volume of the specified local finite element.
  virtual double ElementMinLength(const int LocalElement) const
  {
    return(min_h_[LocalElement]);
  }

  //! Returns the volume of the specified local finite element.
  virtual double ElementMaxLength(const int LocalElement) const
  {
    return(max_h_[LocalElement]);
  }

  //! Returns the volume of the specified local finite element.
  virtual double ElementVolume(const int LocalElement) const
  {
    throw(Exception(__FILE__, __LINE__,
                    "Feature not implemented"));
  }

  //! Returns the area of the specified local face.
  virtual double FaceArea(const int LocalFace) const
  {
    throw(Exception(__FILE__, __LINE__,
                    "Feature not implemented"));
  }

  // @}
  // @{ Maps and import/export

  //! Returns a reference to the map representing the vertex distribution.
  virtual const Epetra_Map& VertexMap() const
  {
    return(*VertexMap_);
  }

  //! Returns a reference to the map representing the distribution of rows.
  virtual const Epetra_Map& RowMap() const
  {
    return(*RowMap_);
  }

  //! Exports distributed object from RowMap() to VertexMap().
  virtual void ExportToVertexMap(const Epetra_DistObject& RowObject,
                                 Epetra_DistObject& VertexObject) const
  {
    throw(Exception(__FILE__, __LINE__,
                    "Feature not implemented"));
  }

  //! Exports distributed object from VertexMap() to RowMap().
  virtual void ExportToRowMap(const Epetra_DistObject& RowObject,
                              Epetra_DistObject& VertexObject) const
  {
    throw(Exception(__FILE__, __LINE__,
                    "Feature not implemented"));
  }

  //! Returns a reference to the communicator object.
  virtual const Epetra_Comm& Comm() const
  {
    return(Comm_);
  }

  // @}
private:

  const Epetra_Comm& Comm_;
  struct triangulateio in_, out_;
  Epetra_Map* VertexMap_;
  Epetra_Map* RowMap_;
  vector<double> min_h_;
  vector<double> max_h_;

}; // class TRIANGLEGrid

} // namespace FiniteElements
} // namespace Galeri
#endif
