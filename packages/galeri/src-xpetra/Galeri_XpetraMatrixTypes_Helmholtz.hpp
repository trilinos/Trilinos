// @HEADER
//
// ***********************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
/*
  Direct translation of parts of Galeri to use Tpetra or Xpetra rather than Epetra. Epetra also supported.
*/

// TODO: rename variables (camelCase)

#ifndef GALERI_XPETRAMATRIXTYPES_HELMHOLTZ_HPP
#define GALERI_XPETRAMATRIXTYPES_HELMHOLTZ_HPP
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayView.hpp>
#include "Galeri_config.h"
#include "Galeri_MatrixTraits.hpp"
#include "Galeri_Problem.hpp"

namespace Galeri {

  namespace Xpetra {

    template <typename Scalar, typename GlobalOrdinal>
    void GetPMLvalues(const GlobalOrdinal i,
		      const GlobalOrdinal nx,          const GlobalOrdinal ny,
		      const double h,                  const double delta,
		      const double Dx,                 const double Dy,
		      const double LBx,                const double RBx,
		      const double LBy,                const double RBy,
		      Scalar& sx_left,        Scalar& sx_center,        Scalar& sx_right,
		      Scalar& sy_left,        Scalar& sy_center,        Scalar& sy_right);

    template <typename Scalar, typename GlobalOrdinal>
    void GetPMLvalues(const GlobalOrdinal i,
                      const GlobalOrdinal nx,          const GlobalOrdinal ny,            const GlobalOrdinal nz,
		      const double h,                  const double delta,
		      const double Dx,                 const double Dy,                   const double Dz,
		      const double LBx,                const double RBx,
		      const double LBy,                const double RBy,
		      const double LBz,                const double RBz,
		      Scalar& sx_left,        Scalar& sx_center,        Scalar& sx_right,
		      Scalar& sy_left,        Scalar& sy_center,        Scalar& sy_right,
		      Scalar& sz_left,        Scalar& sz_center,        Scalar& sz_right);

    /* ****************************************************************************************************** *
     *    Helmholtz 1D
     * ****************************************************************************************************** */
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
    Teuchos::RCP<Matrix>
    TriDiag_Helmholtz(const Teuchos::RCP<const Map> & map, const GlobalOrdinal nx,
		      const double h, const double omega, const Scalar shift) {

      Teuchos::RCP<Matrix> mtx = MatrixTraits<Map,Matrix>::Build(map, 3);
      LocalOrdinal NumMyElements = map->getNodeNumElements();
      Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getNodeElementList();
      Teuchos::RCP<const Teuchos::Comm<int> > comm = map->getComm();
      GlobalOrdinal NumGlobalElements = map->getGlobalNumElements();
      GlobalOrdinal NumEntries;
      LocalOrdinal nnz=2;
      std::vector<Scalar> Values(nnz);
      std::vector<GlobalOrdinal> Indices(nnz);
      Scalar one = (Scalar) 1.0;
      Scalar two = (Scalar) 2.0;
      comm->barrier();
      if (comm->getRank() == 0) {
        std::cout << "starting global insert" << std::endl;
      }
      Teuchos::RCP<Teuchos::Time> timer = rcp(new Teuchos::Time("TriDiag global insert"));
      timer->start(true);
      for (LocalOrdinal i = 0; i < NumMyElements; ++i) {
	if (MyGlobalElements[i] == 0) {
	  // off-diagonal for first row
	  Indices[0] = 1;
	  NumEntries = 1;
	  Values[0] = -one;
	}
	else if (MyGlobalElements[i] == NumGlobalElements - 1) {
	  // off-diagonal for last row
	  Indices[0] = NumGlobalElements - 2;
	  NumEntries = 1;
	  Values[0] = -one;
	}
	else {
	  // off-diagonal for internal row
	  Indices[0] = MyGlobalElements[i] - 1;
	  Indices[1] = MyGlobalElements[i] + 1;
	  Values[0] = -one;
	  Values[1] = -one;
	  NumEntries = 2;
	}
	// put the off-diagonal entries
	// Xpetra wants ArrayViews (sigh)
	Teuchos::ArrayView<Scalar> av(&Values[0],NumEntries);
	Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0],NumEntries);
	mtx->insertGlobalValues(MyGlobalElements[i], iv, av);
	// Put in the diagonal entry
	mtx->insertGlobalValues(MyGlobalElements[i],
				Teuchos::tuple<GlobalOrdinal>(MyGlobalElements[i]),
				Teuchos::tuple<Scalar>(two-shift*omega*omega*h*h) );
      }
      timer->stop();
      timer = rcp(new Teuchos::Time("TriDiag fillComplete"));
      timer->start(true);
      mtx->fillComplete();
      timer->stop();
      return mtx;

    } //TriDiag_Helmholtz

    /* ****************************************************************************************************** *
     *    Helmholtz 2D
     * ****************************************************************************************************** */
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
    Teuchos::RCP<Matrix>
    Cross2D_Helmholtz(const Teuchos::RCP<const Map> & map,
		      const GlobalOrdinal nx,      const GlobalOrdinal ny,
		      const double h,              const double delta,
		      const int PMLgridptsx_left,  const int PMLgridptsx_right,
		      const int PMLgridptsy_left,  const int PMLgridptsy_right,
                      const double omega,          const Scalar shift) {

      Teuchos::RCP<Matrix> mtx = MatrixTraits<Map,Matrix>::Build(map, 5);
      LocalOrdinal NumMyElements = map->getNodeNumElements();
      Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getNodeElementList();
      GlobalOrdinal left, right, lower, upper, center;
      LocalOrdinal nnz=5;
      std::vector<Scalar> Values(nnz);
      std::vector<GlobalOrdinal> Indices(nnz);

      double LBx, RBx, LBy, RBy, Dx, Dy;
      Scalar sx_left, sx_center, sx_right;
      Scalar sy_left, sy_center, sy_right;
      // Calculate some parameters
      Dx = ((double) nx-1)*h;
      Dy = ((double) ny-1)*h;
      LBx = ((double) PMLgridptsx_left)*h;
      RBx = Dx-((double) PMLgridptsx_right)*h;
      LBy = ((double) PMLgridptsy_left)*h;
      RBy = Dy-((double) PMLgridptsy_right)*h;

      for (LocalOrdinal i = 0; i < NumMyElements; ++i)  {

        size_t numEntries = 0;
        center = MyGlobalElements[i];
        GetNeighboursCartesian2d(center, nx, ny, left, right, lower, upper);
	GetPMLvalues(center, nx, ny, h, delta,
		     Dx, Dy, LBx, RBx, LBy, RBy,
		     sx_left, sx_center, sx_right,
		     sy_left, sy_center, sy_right);

	if (left != -1) {
	  Indices[numEntries] = left;
	  Values [numEntries] = -(sy_center/sx_center + sy_center/sx_left)/2.0;
	  numEntries++;
	}
	if (right != -1) {
	  Indices[numEntries] = right;
	  Values [numEntries] = -(sy_center/sx_center + sy_center/sx_right)/2.0;
	  numEntries++;
	}
	if (lower != -1) {
	  Indices[numEntries] = lower;
	  Values [numEntries] = -(sx_center/sy_center + sx_center/sy_left)/2.0;
	  numEntries++;
	}
	if (upper != -1) {
	  Indices[numEntries] = upper;
	  Values [numEntries] = -(sx_center/sy_center + sx_center/sy_right)/2.0;
	  numEntries++;
	}

	// diagonal
	Scalar z = (Scalar) 0.0;
	for (size_t j = 0; j < numEntries; j++)
	  z -= Values[j];

	// mass matrix term
	z -= shift*omega*omega*h*h*sx_center*sy_center;

	Indices[numEntries] = center;
	Values [numEntries] = z;
	numEntries++;

	Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0], numEntries);
	Teuchos::ArrayView<Scalar>        av(&Values[0],  numEntries);
	mtx->insertGlobalValues(center, iv, av);

      }

      mtx->fillComplete();

      return mtx;

    } //Cross2D_Helmholtz

    /* ****************************************************************************************************** *
     *    Helmholtz 3D
     * ****************************************************************************************************** */
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
    Teuchos::RCP<Matrix>
    Cross3D_Helmholtz(const Teuchos::RCP<const Map> & map,
		      const GlobalOrdinal nx,      const GlobalOrdinal ny,        const GlobalOrdinal nz,
                      const double h,              const double delta,
		      const int PMLgridptsx_left,  const int PMLgridptsx_right,
		      const int PMLgridptsy_left,  const int PMLgridptsy_right,
		      const int PMLgridptsz_left,  const int PMLgridptsz_right,
                      const double omega,          const Scalar shift) {

      Teuchos::RCP<Matrix> mtx = MatrixTraits<Map,Matrix>::Build(map, 7);

      LocalOrdinal                               NumMyElements = map->getNodeNumElements();
      Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getNodeElementList();

      GlobalOrdinal left, right, bottom, top, front, back, center;
      std::vector<Scalar> Values(7);
      std::vector<GlobalOrdinal> Indices(7);

      double LBx, RBx, LBy, RBy, LBz, RBz, Dx, Dy, Dz;
      Scalar sx_left, sx_center, sx_right;
      Scalar sy_left, sy_center, sy_right;
      Scalar sz_left, sz_center, sz_right;
      // Calculate some parameters
      Dx = ((double) nx-1)*h;
      Dy = ((double) ny-1)*h;
      Dz = ((double) nz-1)*h;
      LBx = ((double) PMLgridptsx_left)*h;
      RBx = Dx-((double) PMLgridptsx_right)*h;
      LBy = ((double) PMLgridptsy_left)*h;
      RBy = Dy-((double) PMLgridptsy_right)*h;
      LBz = ((double) PMLgridptsz_left)*h;
      RBz = Dz-((double) PMLgridptsz_right)*h;

      for (GlobalOrdinal i = 0; i < NumMyElements; ++i) {

        size_t numEntries = 0;
        center = MyGlobalElements[i];
        GetNeighboursCartesian3d(center, nx, ny, nz, left, right, front, back, bottom, top);
	GetPMLvalues(center, nx, ny, nz, h, delta, Dx, Dy, Dz,
		     LBx, RBx, LBy, RBy, LBz, RBz,
		     sx_left, sx_center, sx_right,
		     sy_left, sy_center, sy_right,
		     sz_left, sz_center, sz_right);

	if (left != -1) {
	  Indices[numEntries] = left;
	  Values [numEntries] = -(sy_center*sz_center/sx_center + sy_center*sz_center/sx_left)/2.0;
	  numEntries++;
	}
	if (right != -1) {
	  Indices[numEntries] = right;
	  Values [numEntries] = -(sy_center*sz_center/sx_center + sy_center*sz_center/sx_right)/2.0;
	  numEntries++;
	}
	if (front != -1) {
	  Indices[numEntries] = front;
	  Values [numEntries] = -(sx_center*sz_center/sy_center + sx_center*sz_center/sy_right)/2.0;
	  numEntries++;
	}
	if (back != -1) {
	  Indices[numEntries] = back;
	  Values [numEntries] = -(sx_center*sz_center/sy_center + sx_center*sz_center/sy_left)/2.0;
	  numEntries++;
	}
	if (bottom != -1) {
	  Indices[numEntries] = bottom;
	  Values [numEntries] = -(sx_center*sy_center/sz_center + sx_center*sy_center/sz_left)/2.0;
	  numEntries++;
	}
	if (top != -1) {
	  Indices[numEntries] = top;
	  Values [numEntries] = -(sx_center*sy_center/sz_center + sx_center*sy_center/sz_right)/2.0;
	  numEntries++;
	}

	// diagonal
	Scalar z = (Scalar) 0.0;
	for (size_t j = 0; j < numEntries; j++)
	  z -= Values[j];

	// mass matrix term
	z -= shift*omega*omega*h*h*sx_center*sy_center*sz_center;

	Indices[numEntries] = center;
	Values [numEntries] = z;
	numEntries++;

	Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0], numEntries);
	Teuchos::ArrayView<Scalar>        av(&Values[0],  numEntries);
	mtx->insertGlobalValues(center, iv, av);

      }

      mtx->fillComplete();

      return mtx;

    } //Cross3D_Helmholtz

    /* ****************************************************************************************************** *
     *    Utilities
     * ****************************************************************************************************** */
    template <typename Scalar, typename GlobalOrdinal>
    void GetPMLvalues(const GlobalOrdinal i,
		      const GlobalOrdinal nx,          const GlobalOrdinal ny,
		      const double h,                  const double delta,
		      const double Dx,                 const double Dy,
		      const double LBx,                const double RBx,
		      const double LBy,                const double RBy,
		      Scalar& sx_left,        Scalar& sx_center,        Scalar& sx_right,
		      Scalar& sy_left,        Scalar& sy_center,        Scalar& sy_right)
    {
      GlobalOrdinal ix, iy;
      ix = i % nx;
      iy = (i - ix) / nx;
      double xcoord, ycoord;
      xcoord=((double) ix)*h;
      ycoord=((double) iy)*h;
      double Lx=max(LBx,Dx-RBx); if(Lx==0.0) { Lx=1.0; }
      double Ly=max(LBy,Dy-RBy); if(Ly==0.0) { Ly=1.0; }
      GetStretch(xcoord-h,delta,LBx,RBx,Lx,sx_left);
      GetStretch(xcoord+0,delta,LBx,RBx,Lx,sx_center);
      GetStretch(xcoord+h,delta,LBx,RBx,Lx,sx_right);
      GetStretch(ycoord-h,delta,LBy,RBy,Ly,sy_left);
      GetStretch(ycoord+0,delta,LBy,RBy,Ly,sy_center);
      GetStretch(ycoord+h,delta,LBy,RBy,Ly,sy_right);
    }

    template <typename Scalar, typename GlobalOrdinal>
    void GetPMLvalues(const GlobalOrdinal i,
                      const GlobalOrdinal nx,          const GlobalOrdinal ny,            const GlobalOrdinal nz,
		      const double h,                  const double delta,
		      const double Dx,                 const double Dy,                   const double Dz,
		      const double LBx,                const double RBx,
		      const double LBy,                const double RBy,
		      const double LBz,                const double RBz,
		      Scalar& sx_left,        Scalar& sx_center,        Scalar& sx_right,
		      Scalar& sy_left,        Scalar& sy_center,        Scalar& sy_right,
		      Scalar& sz_left,        Scalar& sz_center,        Scalar& sz_right)
    {
      GlobalOrdinal ixy, ix, iy, iz;
      ixy = i % (nx * ny);
      iz = (i - ixy) / (nx * ny);
      ix = ixy % nx;
      iy = (ixy - ix) / nx;
      double xcoord, ycoord, zcoord;
      xcoord=((double) ix)*h;
      ycoord=((double) iy)*h;
      zcoord=((double) iz)*h;
      double Lx=max(LBx,Dx-RBx); if(Lx==0.0) { Lx=1.0; }
      double Ly=max(LBy,Dy-RBy); if(Ly==0.0) { Ly=1.0; }
      double Lz=max(LBz,Dz-RBz); if(Lz==0.0) { Lz=1.0; }
      GetStretch(xcoord-h,delta,LBx,RBx,Lx,sx_left);
      GetStretch(xcoord+0,delta,LBx,RBx,Lx,sx_center);
      GetStretch(xcoord+h,delta,LBx,RBx,Lx,sx_right);
      GetStretch(ycoord-h,delta,LBy,RBy,Ly,sy_left);
      GetStretch(ycoord+0,delta,LBy,RBy,Ly,sy_center);
      GetStretch(ycoord+h,delta,LBy,RBy,Ly,sy_right);
      GetStretch(zcoord-h,delta,LBz,RBz,Lz,sz_left);
      GetStretch(zcoord+0,delta,LBz,RBz,Lz,sz_center);
      GetStretch(zcoord+h,delta,LBz,RBz,Lz,sz_right);
    }

    template <typename Scalar>
    void GetStretch(double x, const double& delta, const double& LB, const double& RB, double& L, Scalar& stretch) {
      double cpxpart;
      if(x<LB)        { cpxpart = delta*pow((x-LB)/L,2);  }
      else if(x>RB)   { cpxpart = delta*pow((x-RB)/L,2);  }
      else            { cpxpart = 0.0;                    }
      Scalar sx(1.0,cpxpart);
      stretch = sx;
    }

  } // namespace Xpetra
} // namespace Galeri


#endif //ifndef GALERI_XPETRAMATRIXTYPES_HELMHOLTZ_HPP
