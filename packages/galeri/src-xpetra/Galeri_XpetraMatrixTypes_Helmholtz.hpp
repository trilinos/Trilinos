// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
#include "Galeri_XpetraMatrixTypes.hpp"
#include "Galeri_VelocityModel.hpp"

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

    template <typename Scalar>
    void GetStretch(double x, const double& delta, const double& LB, const double& RB, double& L, Scalar& stretch);

    /* ****************************************************************************************************** *
     *    Helmholtz 1D
     * ****************************************************************************************************** */
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
    Teuchos::RCP<Matrix>
    TriDiag_Helmholtz(const Teuchos::RCP<const Map> & map, const GlobalOrdinal nx,
                      const double h, const double omega, const Scalar shift) {

      Teuchos::RCP<Matrix> mtx = MatrixTraits<Map,Matrix>::Build(map, 3);
      LocalOrdinal NumMyElements = map->getLocalNumElements();
      Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getLocalElementList();
      Teuchos::RCP<const Teuchos::Comm<int> > comm = map->getComm();
      GlobalOrdinal NumGlobalElements = map->getGlobalNumElements();
      GlobalOrdinal NumEntries;
      LocalOrdinal nnz=2;
      std::vector<Scalar> Values(nnz);
      std::vector<GlobalOrdinal> Indices(nnz);
      Scalar one = (Scalar) 1.0;
      Scalar two = (Scalar) 2.0;
      comm->barrier();
      Teuchos::RCP<Teuchos::Time> timer = Teuchos::rcp(new Teuchos::Time("TriDiag global insert"));
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
      timer = Teuchos::rcp(new Teuchos::Time("TriDiag fillComplete"));
      timer->start(true);
      mtx->fillComplete();
      timer->stop();
      return mtx;

    } //TriDiag_Helmholtz

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
    std::pair< Teuchos::RCP<Matrix>, Teuchos::RCP<Matrix> >
    TriDiag_Helmholtz_Pair(const Teuchos::RCP<const Map> & map, const GlobalOrdinal nx,
                           const double h, const double omega, const Scalar shift) {

      Teuchos::RCP<Matrix> ktx = MatrixTraits<Map,Matrix>::Build(map, 3);
      Teuchos::RCP<Matrix> mtx = MatrixTraits<Map,Matrix>::Build(map, 1);
      LocalOrdinal NumMyElements = map->getLocalNumElements();
      Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getLocalElementList();
      Teuchos::RCP<const Teuchos::Comm<int> > comm = map->getComm();
      GlobalOrdinal NumGlobalElements = map->getGlobalNumElements();
      GlobalOrdinal NumEntries;
      LocalOrdinal nnz=2;
      std::vector<Scalar> Values(nnz);
      std::vector<GlobalOrdinal> Indices(nnz);
      Scalar one = (Scalar) 1.0;
      comm->barrier();
      if (comm->getRank() == 0) {
        std::cout << "starting global insert" << std::endl;
      }
      Teuchos::RCP<Teuchos::Time> timer = Teuchos::rcp(new Teuchos::Time("TriDiag global insert"));
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
        ktx->insertGlobalValues(MyGlobalElements[i], iv, av);
        // Put in the diagonal entry
        mtx->insertGlobalValues(MyGlobalElements[i],
                                Teuchos::tuple<GlobalOrdinal>(MyGlobalElements[i]),
                                Teuchos::tuple<Scalar>(h*h) );
      }
      timer->stop();
      timer = Teuchos::rcp(new Teuchos::Time("TriDiag fillComplete"));
      timer->start(true);
      ktx->fillComplete();
      mtx->fillComplete();
      timer->stop();
      std::pair< Teuchos::RCP<Matrix>, Teuchos::RCP<Matrix> > system;
      system=std::make_pair(ktx,mtx);
      return system;

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
                      const double omega,          const Scalar shift,
		      const int model) {

      Teuchos::RCP<Matrix> mtx = MatrixTraits<Map,Matrix>::Build(map, 5);
      LocalOrdinal NumMyElements = map->getLocalNumElements();
      Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getLocalElementList();
      GlobalOrdinal left, right, lower, upper, center;
      LocalOrdinal nnz=5;
      std::vector<Scalar> Values(nnz);
      std::vector<GlobalOrdinal> Indices(nnz);
      VelocityModel<Scalar,LocalOrdinal,GlobalOrdinal> velocitymodel(2,model);

      double LBx, RBx, LBy, RBy, Dx, Dy;
      Scalar sx_left, sx_center, sx_right;
      Scalar sy_left, sy_center, sy_right;
      Scalar c;
      // Calculate some parameters
      Dx = ((double) nx-1)*h;
      Dy = ((double) ny-1)*h;
      LBx = ((double) PMLgridptsx_left)*h;
      RBx = Dx-((double) PMLgridptsx_right)*h;
      LBy = ((double) PMLgridptsy_left)*h;
      RBy = Dy-((double) PMLgridptsy_right)*h;

      for (LocalOrdinal i = 0; i < NumMyElements; ++i)  {

	// calculate PML functions
        size_t numEntries = 0;
        center = MyGlobalElements[i];
        GetNeighboursCartesian2d(center, nx, ny, left, right, lower, upper);
        GetPMLvalues(center, nx, ny, h, delta,
                     Dx, Dy, LBx, RBx, LBy, RBy,
                     sx_left, sx_center, sx_right,
                     sy_left, sy_center, sy_right);
	// get velocity
	GlobalOrdinal ix, iy;
	ix = i % nx;
	iy = (i - ix) / nx;
	double xcoord, ycoord;
	xcoord=((double) ix)*h;
	ycoord=((double) iy)*h;
	c = velocitymodel.getVelocity(xcoord,ycoord,0.0);

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
        z -= shift*omega*omega*h*h*sx_center*sy_center/(c*c);

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

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
    std::pair< Teuchos::RCP<Matrix>, Teuchos::RCP<Matrix> >
    Cross2D_Helmholtz_Pair(const Teuchos::RCP<const Map> & map,
                           const GlobalOrdinal nx,      const GlobalOrdinal ny,
                           const double h,              const double delta,
                           const int PMLgridptsx_left,  const int PMLgridptsx_right,
                           const int PMLgridptsy_left,  const int PMLgridptsy_right,
                           const double omega,          const Scalar shift,
			   const int model) {

      Teuchos::RCP<Matrix> ktx = MatrixTraits<Map,Matrix>::Build(map, 5);
      Teuchos::RCP<Matrix> mtx = MatrixTraits<Map,Matrix>::Build(map, 1);
      LocalOrdinal NumMyElements = map->getLocalNumElements();
      Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getLocalElementList();
      GlobalOrdinal left, right, lower, upper, center;
      LocalOrdinal nnz=5;
      std::vector<Scalar> Values(nnz);
      std::vector<GlobalOrdinal> Indices(nnz);
      VelocityModel<Scalar,LocalOrdinal,GlobalOrdinal> velocitymodel(2,model);

      double LBx, RBx, LBy, RBy, Dx, Dy;
      Scalar sx_left, sx_center, sx_right;
      Scalar sy_left, sy_center, sy_right;
      Scalar c;
      // Calculate some parameters
      Dx = ((double) nx-1)*h;
      Dy = ((double) ny-1)*h;
      LBx = ((double) PMLgridptsx_left)*h;
      RBx = Dx-((double) PMLgridptsx_right)*h;
      LBy = ((double) PMLgridptsy_left)*h;
      RBy = Dy-((double) PMLgridptsy_right)*h;

      for (LocalOrdinal i = 0; i < NumMyElements; ++i)  {

	// calculate PML functions
        size_t numEntries = 0;
        center = MyGlobalElements[i];
        GetNeighboursCartesian2d(center, nx, ny, left, right, lower, upper);
        GetPMLvalues(center, nx, ny, h, delta,
                     Dx, Dy, LBx, RBx, LBy, RBy,
                     sx_left, sx_center, sx_right,
                     sy_left, sy_center, sy_right);
	// get velocity
	GlobalOrdinal ix, iy;
	ix = i % nx;
	iy = (i - ix) / nx;
	double xcoord, ycoord;
	xcoord=((double) ix)*h;
	ycoord=((double) iy)*h;
	c = velocitymodel.getVelocity(xcoord,ycoord,0.0);

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

        Indices[numEntries] = center;
        Values [numEntries] = z;
        numEntries++;

        Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0], numEntries);
        Teuchos::ArrayView<Scalar>        av(&Values[0],  numEntries);
        ktx->insertGlobalValues(center, iv, av);
        mtx->insertGlobalValues(center,
                                Teuchos::tuple<GlobalOrdinal>(center),
                                Teuchos::tuple<Scalar>(h*h*sx_center*sy_center/(c*c)) );

      }

      ktx->fillComplete();
      mtx->fillComplete();
      std::pair< Teuchos::RCP<Matrix>, Teuchos::RCP<Matrix> > system;
      system=std::make_pair(ktx,mtx);
      return system;

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
                      const double omega,          const Scalar shift,
		      const int model) {

      Teuchos::RCP<Matrix> mtx = MatrixTraits<Map,Matrix>::Build(map, 7);

      LocalOrdinal                               NumMyElements = map->getLocalNumElements();
      Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getLocalElementList();

      GlobalOrdinal left, right, bottom, top, front, back, center;
      std::vector<Scalar> Values(7);
      std::vector<GlobalOrdinal> Indices(7);
      VelocityModel<Scalar,LocalOrdinal,GlobalOrdinal> velocitymodel(3,model);

      double LBx, RBx, LBy, RBy, LBz, RBz, Dx, Dy, Dz;
      Scalar sx_left, sx_center, sx_right;
      Scalar sy_left, sy_center, sy_right;
      Scalar sz_left, sz_center, sz_right;
      Scalar c;
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

	// calculate PML functions
        size_t numEntries = 0;
        center = MyGlobalElements[i];
        GetNeighboursCartesian3d(center, nx, ny, nz, left, right, front, back, bottom, top);
        GetPMLvalues(center, nx, ny, nz, h, delta, Dx, Dy, Dz,
                     LBx, RBx, LBy, RBy, LBz, RBz,
                     sx_left, sx_center, sx_right,
                     sy_left, sy_center, sy_right,
                     sz_left, sz_center, sz_right);
	// get velocity
	GlobalOrdinal ixy, ix, iy, iz;
	ixy = i % (nx * ny);
	iz = (i - ixy) / (nx * ny);
	ix = ixy % nx;
	iy = (ixy - ix) / nx;
	double xcoord, ycoord, zcoord;
	xcoord=((double) ix)*h;
	ycoord=((double) iy)*h;
	zcoord=((double) iz)*h;
	c = velocitymodel.getVelocity(xcoord,ycoord,zcoord);

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
        z -= shift*omega*omega*h*h*sx_center*sy_center*sz_center/(c*c);

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

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
        std::pair< Teuchos::RCP<Matrix>, Teuchos::RCP<Matrix> >
        Cross3D_Helmholtz_Pair(const Teuchos::RCP<const Map> & map,
                               const GlobalOrdinal nx,      const GlobalOrdinal ny,        const GlobalOrdinal nz,
                               const double h,              const double delta,
                               const int PMLgridptsx_left,  const int PMLgridptsx_right,
                               const int PMLgridptsy_left,  const int PMLgridptsy_right,
                               const int PMLgridptsz_left,  const int PMLgridptsz_right,
                               const double omega,          const Scalar shift,
			       const int model) {

          Teuchos::RCP<Matrix> ktx = MatrixTraits<Map,Matrix>::Build(map, 7);
          Teuchos::RCP<Matrix> mtx = MatrixTraits<Map,Matrix>::Build(map, 1);

          LocalOrdinal                               NumMyElements = map->getLocalNumElements();
          Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getLocalElementList();

          GlobalOrdinal left, right, bottom, top, front, back, center;
          std::vector<Scalar> Values(7);
          std::vector<GlobalOrdinal> Indices(7);
	  VelocityModel<Scalar,LocalOrdinal,GlobalOrdinal> velocitymodel(3,model);

          double LBx, RBx, LBy, RBy, LBz, RBz, Dx, Dy, Dz;
          Scalar sx_left, sx_center, sx_right;
          Scalar sy_left, sy_center, sy_right;
          Scalar sz_left, sz_center, sz_right;
	  Scalar c;
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

	    // calculate PML functions
            size_t numEntries = 0;
            center = MyGlobalElements[i];
            GetNeighboursCartesian3d(center, nx, ny, nz, left, right, front, back, bottom, top);
            GetPMLvalues(center, nx, ny, nz, h, delta, Dx, Dy, Dz,
                         LBx, RBx, LBy, RBy, LBz, RBz,
                         sx_left, sx_center, sx_right,
                         sy_left, sy_center, sy_right,
                         sz_left, sz_center, sz_right);
	    // get velocity
	    GlobalOrdinal ixy, ix, iy, iz;
	    ixy = i % (nx * ny);
	    iz = (i - ixy) / (nx * ny);
	    ix = ixy % nx;
	    iy = (ixy - ix) / nx;
	    double xcoord, ycoord, zcoord;
	    xcoord=((double) ix)*h;
	    ycoord=((double) iy)*h;
	    zcoord=((double) iz)*h;
	    c = velocitymodel.getVelocity(xcoord,ycoord,zcoord);
	    
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

            Indices[numEntries] = center;
            Values [numEntries] = z;
            numEntries++;

            Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0], numEntries);
            Teuchos::ArrayView<Scalar>        av(&Values[0],  numEntries);
            ktx->insertGlobalValues(center, iv, av);
            mtx->insertGlobalValues(center,
                                    Teuchos::tuple<GlobalOrdinal>(center),
                                    Teuchos::tuple<Scalar>(h*h*sx_center*sy_center*sz_center/(c*c)) );

          }

          ktx->fillComplete();
          mtx->fillComplete();
          std::pair< Teuchos::RCP<Matrix>, Teuchos::RCP<Matrix> > system;
          system=std::make_pair(ktx,mtx);
          return system;

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
      double Lx=std::max(LBx,Dx-RBx); if(Lx==0.0) { Lx=1.0; }
      double Ly=std::max(LBy,Dy-RBy); if(Ly==0.0) { Ly=1.0; }
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
      double Lx=std::max(LBx,Dx-RBx); if(Lx==0.0) { Lx=1.0; }
      double Ly=std::max(LBy,Dy-RBy); if(Ly==0.0) { Ly=1.0; }
      double Lz=std::max(LBz,Dz-RBz); if(Lz==0.0) { Lz=1.0; }
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
