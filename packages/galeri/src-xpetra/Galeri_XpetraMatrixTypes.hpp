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

#ifndef GALERI_XPETRAMATRIXTYPES_HPP
#define GALERI_XPETRAMATRIXTYPES_HPP
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "Galeri_config.h"
#include "Galeri_MatrixTraits.hpp"
#include "Galeri_Problem.hpp"

namespace Galeri {

  namespace Xpetra {

    /* prototypes */
    template  <typename GlobalOrdinal>
    bool IsBoundary2d(const GlobalOrdinal i, const GlobalOrdinal nx, const GlobalOrdinal ny);
    template  <typename GlobalOrdinal>
    bool IsBoundary3d(const GlobalOrdinal i, const GlobalOrdinal nx, const GlobalOrdinal ny, const GlobalOrdinal nz);

    template <typename GlobalOrdinal>
    void GetNeighboursCartesian2d(const GlobalOrdinal i,
                                  const GlobalOrdinal nx, const GlobalOrdinal ny,
                                  GlobalOrdinal& left,  GlobalOrdinal& right,
                                  GlobalOrdinal& lower, GlobalOrdinal& upper);

    template <typename GlobalOrdinal>
    void GetNeighboursCartesian2d(GlobalOrdinal i,      GlobalOrdinal nx,      const GlobalOrdinal ny,
                                  GlobalOrdinal& left,  GlobalOrdinal& right,  GlobalOrdinal& lower,  GlobalOrdinal& upper,
                                  GlobalOrdinal& left2, GlobalOrdinal& right2, GlobalOrdinal& lower2, GlobalOrdinal& upper2);

    template <typename GlobalOrdinal>
    void GetNeighboursCartesian3d(const GlobalOrdinal i,
                                  const GlobalOrdinal nx, const GlobalOrdinal ny, const GlobalOrdinal nz,
                                  GlobalOrdinal& left,   GlobalOrdinal& right,
                                  GlobalOrdinal& front,  GlobalOrdinal& back,
                                  GlobalOrdinal& bottom, GlobalOrdinal& top);

    template <typename GlobalOrdinal, typename Scalar>
    void Fill9PointStencil(const GlobalOrdinal center,
                           std::vector<Scalar>& Values, std::vector<GlobalOrdinal>& Indices, size_t& numEntries,
                           const GlobalOrdinal nx, const GlobalOrdinal ny, const GlobalOrdinal nz,
                           const Scalar b,  const Scalar c,  const Scalar d,  const Scalar e,
                           const Scalar z1, const Scalar z2, const Scalar z3, const Scalar z4,
                           GlobalOrdinal left  = -2, GlobalOrdinal right = -2,
                           GlobalOrdinal lower = -2, GlobalOrdinal upper = -2);
    /* end of prototypes */

    /* ****************************************************************************************************** *
     *    (Scaled) Identity
     * ****************************************************************************************************** */
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
    Teuchos::RCP<Matrix>
    Identity(const Teuchos::RCP<const Map>& map, const Scalar a)
    {
      using Teuchos::TimeMonitor;
      using Teuchos::RCP;
      using Teuchos::rcp;

      Teuchos::RCP<Matrix> mtx = MatrixTraits<Map,Matrix>::Build(map, 1);

      LocalOrdinal NumMyElements = map->getLocalNumElements();
      Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getLocalElementList();

      {
	Teuchos::RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Galeri: Scaled Identity Generation")));

	for (LocalOrdinal i = 0; i < NumMyElements; i++)
	  mtx->insertGlobalValues(MyGlobalElements[i],
				  Teuchos::tuple<GlobalOrdinal>(MyGlobalElements[i]),
				  Teuchos::tuple<Scalar>(a));
      }
      {
	Teuchos::RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Galeri: Scaled Identity fillComplete")));
	mtx->fillComplete();
      }
      return mtx;
    }

    /* ****************************************************************************************************** *
     *    Laplace 1D
     * ****************************************************************************************************** */
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
    Teuchos::RCP<Matrix>
    TriDiag(const Teuchos::RCP<const Map> & map,
            const GlobalOrdinal nx,                            // note: nx unused
            const Scalar a, const Scalar b, const Scalar c)
    {
      using Teuchos::TimeMonitor;
      using Teuchos::RCP;
      using Teuchos::rcp;

      Teuchos::RCP<Matrix> mtx = MatrixTraits<Map,Matrix>::Build(map, 3);

      LocalOrdinal NumMyElements = map->getLocalNumElements();
      Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getLocalElementList();
      GlobalOrdinal indexBase = map->getIndexBase();

      Teuchos::RCP<const Teuchos::Comm<int> > comm = map->getComm();

      GlobalOrdinal NumGlobalElements = map->getGlobalNumElements();

      GlobalOrdinal NumEntries;
      LocalOrdinal  nnz = 2;
      std::vector<Scalar>        Values(nnz);
      std::vector<GlobalOrdinal> Indices(nnz);

      comm->barrier();

      // c a b
      {
	Teuchos::RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Galeri: Laplace 1D Generation")));

	for (LocalOrdinal i = 0; i < NumMyElements; i++) {
	  if (MyGlobalElements[i] == indexBase) {
	    // off-diagonal for first row
	    Indices[0] = 1 + indexBase;
	    Values [0] = c;
	    NumEntries = 1;
	    
	  } else if (MyGlobalElements[i] == NumGlobalElements + indexBase - 1) {
	    // off-diagonal for last row
	    Indices[0] = NumGlobalElements - 2 + indexBase;
	    Values [0] = b;
	    NumEntries = 1;

	  } else {
	    // off-diagonal for internal row
	    Indices[0] = MyGlobalElements[i] - 1;
	    Values [0] = b;
	    Indices[1] = MyGlobalElements[i] + 1;
	    Values [1] = c;
	    NumEntries = 2;
	  }

	  // put the off-diagonal entries
	  // Xpetra wants ArrayViews (sigh)
	  Teuchos::ArrayView<Scalar>        av(&Values [0], NumEntries);
	  Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0], NumEntries);
	  mtx->insertGlobalValues(MyGlobalElements[i], iv, av);

	  // Put in the diagonal entry
	  mtx->insertGlobalValues(MyGlobalElements[i],
				  Teuchos::tuple<GlobalOrdinal>(MyGlobalElements[i]),
				  Teuchos::tuple<Scalar>(a));
	}
      }
      
      {
	Teuchos::RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Galeri: Laplace 1D fillComplete")));
	mtx->fillComplete();
      }

      return mtx;
    }

    /* ****************************************************************************************************** *
     *    Laplace 2D
     * ****************************************************************************************************** */
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
    Teuchos::RCP<Matrix>
    Cross2D(const Teuchos::RCP<const Map>& map,
            const GlobalOrdinal nx, const GlobalOrdinal ny,
            const Scalar a, const Scalar b, const Scalar c,
            const Scalar d, const Scalar e,
            const DirBC DirichletBC = 0, const bool keepBCs = false)
    {
      using Teuchos::TimeMonitor;
      using Teuchos::RCP;
      using Teuchos::rcp;
      
      LocalOrdinal nnz = 5;

      RCP<Matrix> mtx = MatrixTraits<Map,Matrix>::Build(map, nnz);

      LocalOrdinal  numMyElements = map->getLocalNumElements();
      GlobalOrdinal indexBase     = map->getIndexBase();

      Teuchos::ArrayView<const GlobalOrdinal> myGlobalElements = map->getLocalElementList();

      GlobalOrdinal center, left, right, lower, upper;
      std::vector<Scalar>        vals(nnz);
      std::vector<GlobalOrdinal> inds(nnz);

      //    e
      //  b a c
      //    d
      {
	Teuchos::RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Galeri: Laplace 2D Generation")));
	for (LocalOrdinal i = 0; i < numMyElements; ++i)  {
	  size_t n = 0;

	  center = myGlobalElements[i] - indexBase;
	  GetNeighboursCartesian2d(center, nx, ny, left, right, lower, upper);

	  bool isDirichlet = (left  == -1 && (DirichletBC & DIR_LEFT))   ||
                             (right == -1 && (DirichletBC & DIR_RIGHT))  ||
                             (lower == -1 && (DirichletBC & DIR_BOTTOM)) ||
                             (upper == -1 && (DirichletBC & DIR_TOP));

	  if (isDirichlet && keepBCs) {
	    // Dirichlet unknown we want to keep
	    inds[n]   = center;
	    vals[n++] = Teuchos::ScalarTraits<Scalar>::one();

	  } else {
	    // The Neumann b.c. are treated in a sane way. The Dirichlet b.c., however, are treated
	    // insane when the option keepBCs=false. Speicifically, in this case we don't want to keep
	    // Dirichlet b.c., but that would result in inconsistency between the map and the number of
	    // degrees of freedom, plus the problem with GIDs. Therefore, we virtually expand domain by
	    // one node in the direction of the Dirichlet b.c., and then assume that that node was
	    // not kept. But we use an old GIDs. So yes, that's weird.

	    if (left  != -1) { inds[n] = left;  vals[n++] = b; }
	    if (right != -1) { inds[n] = right; vals[n++] = c; }
	    if (lower != -1) { inds[n] = lower; vals[n++] = d; }
	    if (upper != -1) { inds[n] = upper; vals[n++] = e; }

	    // diagonal
	    Scalar z = a;
	    if (IsBoundary2d(center, nx, ny) && !isDirichlet) {
	      // Neumann boundary unknown (diagonal = sum of all offdiagonal)
	      z = Teuchos::ScalarTraits<Scalar>::zero();
	      for (size_t j = 0; j < n; j++)
		z -= vals[j];
	    }
	    inds[n]   = center;
	    vals[n++] = z;
	  }

	  for (size_t j = 0; j < n; j++)
	    inds[j] += indexBase;

	  Teuchos::ArrayView<GlobalOrdinal> iv(&inds[0], n);
	  Teuchos::ArrayView<Scalar>        av(&vals[0], n);
	  mtx->insertGlobalValues(myGlobalElements[i], iv, av);
	}
      }

      {
	Teuchos::RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Galeri: Laplace 2D FillComplete")));
	mtx->fillComplete();
      }

      return mtx;
    }

    /* ****************************************************************************************************** *
     *    Star2D
     * ****************************************************************************************************** */
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
    Teuchos::RCP<Matrix>
    Star2D(const Teuchos::RCP<const Map> & map,
           const GlobalOrdinal nx, const GlobalOrdinal ny,
           const Scalar a,  const Scalar b, const Scalar c,
           const Scalar d,  const Scalar e,
           const Scalar z1, const Scalar z2,
           const Scalar z3, const Scalar z4,
           const DirBC DirichletBC = 0, const bool keepBCs = false)
    {
      using Teuchos::TimeMonitor;
      using Teuchos::RCP;
      using Teuchos::rcp;

      LocalOrdinal nnz = 9;

      Teuchos::RCP<Matrix> mtx = MatrixTraits<Map,Matrix>::Build(map, nnz);

      LocalOrdinal  numMyElements = map->getLocalNumElements();
      GlobalOrdinal indexBase     = map->getIndexBase();

      Teuchos::ArrayView<const GlobalOrdinal> myGlobalElements = map->getLocalElementList();

      GlobalOrdinal center, left, right, lower, upper;
      std::vector<Scalar>        vals(nnz);
      std::vector<GlobalOrdinal> inds(nnz);

      //  z3  e  z4
      //   b  a  c
      //  z1  d  z2
      {
	Teuchos::RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Galeri: Star2D generation")));

	for (LocalOrdinal i = 0; i < numMyElements; ++i) {
	  size_t n = 0;
	  
	  center = myGlobalElements[i] - indexBase;
	  GetNeighboursCartesian2d(center, nx, ny, left, right, lower, upper);

	  bool isDirichlet = (left  == -1 && (DirichletBC & DIR_LEFT))   ||
	                     (right == -1 && (DirichletBC & DIR_RIGHT))  ||
                             (lower == -1 && (DirichletBC & DIR_BOTTOM)) ||
                             (upper == -1 && (DirichletBC & DIR_TOP));

	  if (isDirichlet && keepBCs) {
	    // Dirichlet unknown we want to keep
	    inds[n]   = center;
	    vals[n++] = Teuchos::ScalarTraits<Scalar>::one();

	  } else {
	    // See comments about weird in Cross2D
	    if (left  != -1)                { inds[n] = left;    vals[n++] = b;  }
	    if (right != -1)                { inds[n] = right;   vals[n++] = c;  }
	    if (lower != -1)                { inds[n] = lower;   vals[n++] = d;  }
	    if (upper != -1)                { inds[n] = upper;   vals[n++] = e;  }
	    if (left  != -1 && lower != -1) { inds[n] = lower-1; vals[n++] = z1; }
	    if (right != -1 && lower != -1) { inds[n] = lower+1; vals[n++] = z2; }
	    if (left  != -1 && upper != -1) { inds[n] = upper-1; vals[n++] = z3; }
	    if (right != -1 && upper != -1) { inds[n] = upper+1; vals[n++] = z4; }

	    // diagonal
	    Scalar z = a;
	    if (IsBoundary2d(center, nx, ny) && !isDirichlet) {
	      // Neumann boundary unknown (diagonal = sum of all offdiagonal)
	      z = Teuchos::ScalarTraits<Scalar>::zero();
	      for (size_t j = 0; j < n; j++)
		z -= vals[j];
	    }
	    inds[n]   = center + indexBase;
	    vals[n++] = z;
	  }

	  for (size_t j = 0; j < n; j++)
	    inds[j] += indexBase;

	  Teuchos::ArrayView<GlobalOrdinal> iv(&inds[0], n);
	  Teuchos::ArrayView<Scalar>        av(&vals[0], n);
	  mtx->insertGlobalValues(myGlobalElements[i], iv, av);
	}
      }

      {
	Teuchos::RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Galeri: star2D fillComplete")));
	mtx->fillComplete();
      }

      return mtx;
    }

    /* ****************************************************************************************************** *
     *    BigStar2D (2D Biharmonic operator, for example)
     * ****************************************************************************************************** */
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
    Teuchos::RCP<Matrix>
    BigStar2D(const Teuchos::RCP<const Map> & map,
              const GlobalOrdinal nx, const GlobalOrdinal ny,
              const Scalar a, const Scalar b, const Scalar c,
              const Scalar d, const Scalar e,
              const Scalar z1, const Scalar z2,
              const Scalar z3, const Scalar z4,
              const Scalar bb, const Scalar cc, const Scalar dd, const Scalar ee,
              const DirBC DirichletBC = 0, const bool keepBCs = false)
    {
      using Teuchos::TimeMonitor;
      using Teuchos::RCP;
      using Teuchos::rcp;

      LocalOrdinal nnz = 13;

      Teuchos::RCP<Matrix> mtx = MatrixTraits<Map,Matrix>::Build(map, nnz);

      LocalOrdinal  numMyElements = map->getLocalNumElements();
      GlobalOrdinal indexBase     = map->getIndexBase();

      Teuchos::ArrayView<const GlobalOrdinal> myGlobalElements = map->getLocalElementList();

      GlobalOrdinal center, left, right, lower, upper, left2, right2, lower2, upper2;
      std::vector<Scalar>        vals(nnz);
      std::vector<GlobalOrdinal> inds(nnz);

      //        ee
      //    z3  e  z4
      // bb  b  a  c  cc
      //    z1  d  z2
      //        dd
      {
	Teuchos::RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Galeri: bigStar2D Generation")));

	for (LocalOrdinal i = 0; i < numMyElements; ++i) {
	  size_t n = 0;

	  center = myGlobalElements[i] - indexBase;
	  GetNeighboursCartesian2d(center, nx, ny, left, right, lower, upper, left2, right2, lower2, upper2);

	  bool isDirichlet = (left  == -1 && (DirichletBC & DIR_LEFT))   ||
                             (right == -1 && (DirichletBC & DIR_RIGHT))  ||
                             (lower == -1 && (DirichletBC & DIR_BOTTOM)) ||
                             (upper == -1 && (DirichletBC & DIR_TOP));

	  if (isDirichlet && keepBCs) {
	    // Dirichlet unknown we want to keep
	    inds[n]   = center;
	    vals[n++] = Teuchos::ScalarTraits<Scalar>::one();

	  } else {
	    // See comments about weird in Cross2D
	    if (left   != -1)                { inds[n] = left;      vals[n++] = b ; }
	    if (right  != -1)                { inds[n] = right;     vals[n++] = c ; }
	    if (lower  != -1)                { inds[n] = lower;     vals[n++] = d ; }
	    if (upper  != -1)                { inds[n] = upper;     vals[n++] = e ; }
	    if (left   != -1 && lower != -1) { inds[n] = lower - 1; vals[n++] = z1; }
	    if (right  != -1 && lower != -1) { inds[n] = lower + 1; vals[n++] = z2; }
	    if (left   != -1 && upper != -1) { inds[n] = upper - 1; vals[n++] = z3; }
	    if (right  != -1 && upper != -1) { inds[n] = upper + 1; vals[n++] = z4; }
	    if (left2  != -1)                { inds[n] = left2;     vals[n++] = bb; }
	    if (right2 != -1)                { inds[n] = right2;    vals[n++] = cc; }
	    if (lower2 != -1)                { inds[n] = lower2;    vals[n++] = dd; }
	    if (upper2 != -1)                { inds[n] = upper2;    vals[n++] = ee; }

	    // diagonal
	    Scalar z = a;
	    if (IsBoundary2d(center, nx, ny) && !isDirichlet) {
	      // Neumann boundary unknown (diagonal = sum of all offdiagonal)
	      z = Teuchos::ScalarTraits<Scalar>::zero();
	      for (size_t j = 0; j < n; j++)
		z -= vals[j];
	    }
	    inds[n]   = center + indexBase;
	    vals[n++] = z;
	  }

	  for (size_t j = 0; j < n; j++)
	    inds[j] += indexBase;

	  Teuchos::ArrayView<GlobalOrdinal> iv(&inds[0], n);
	  Teuchos::ArrayView<Scalar>        av(&vals[0], n);
	  mtx->insertGlobalValues(myGlobalElements[i], iv, av);
	}
      }

      {
	Teuchos::RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Galeri: bigStar2D fillComplete")));

	mtx->fillComplete();
      }

      return mtx;
    }

    /* ****************************************************************************************************** *
     *    Laplace 3D
     * ****************************************************************************************************** */
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
    Teuchos::RCP<Matrix>
    Cross3D(const Teuchos::RCP<const Map>& map,
            const GlobalOrdinal nx, const GlobalOrdinal ny, const GlobalOrdinal nz,
            const Scalar a, const Scalar b, const Scalar c,
            const Scalar d, const Scalar e,
            const Scalar f, const Scalar g,
            const DirBC DirichletBC = 0, const bool keepBCs = false)
    {
      using Teuchos::TimeMonitor;
      using Teuchos::RCP;
      using Teuchos::rcp;

      LocalOrdinal nnz = 7;

      Teuchos::RCP<Matrix> mtx = MatrixTraits<Map,Matrix>::Build(map, nnz);

      LocalOrdinal  numMyElements = map->getLocalNumElements();
      GlobalOrdinal indexBase     = map->getIndexBase();

      Teuchos::ArrayView<const GlobalOrdinal> myGlobalElements = map->getLocalElementList();

      GlobalOrdinal center, left, right, bottom, top, front, back;
      std::vector<GlobalOrdinal> inds(nnz);
      std::vector<Scalar>        vals(nnz);

      //    e
      //  b a c
      //    d
      // + f bottom and g top
      {
	Teuchos::RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Galeri: Laplace 3D generation")));
	for (LocalOrdinal i = 0; i < numMyElements; ++i) {
	  size_t n = 0;

	  center = myGlobalElements[i] - indexBase;
	  GetNeighboursCartesian3d(center, nx, ny, nz,
				   left, right, front, back, bottom, top);

	  bool isDirichlet = (left   == -1 && (DirichletBC & DIR_LEFT))   ||
                             (right  == -1 && (DirichletBC & DIR_RIGHT))  ||
                             (bottom == -1 && (DirichletBC & DIR_BOTTOM)) ||
                             (top    == -1 && (DirichletBC & DIR_TOP))    ||
                             (front  == -1 && (DirichletBC & DIR_FRONT))  ||
                             (back   == -1 && (DirichletBC & DIR_BACK));

	  if (isDirichlet && keepBCs) {
	    // Dirichlet unknown we want to keep
	    inds[n]   = center;
	    vals[n++] = Teuchos::ScalarTraits<Scalar>::one();

	  } else {
	    // See comments about weird in Cross2D
	    if (left   != -1) { inds[n] = left;   vals[n++] = b; }
	    if (right  != -1) { inds[n] = right;  vals[n++] = c; }
	    if (front  != -1) { inds[n] = front;  vals[n++] = d; }
	    if (back   != -1) { inds[n] = back;   vals[n++] = e; }
	    if (bottom != -1) { inds[n] = bottom; vals[n++] = f; }
	    if (top    != -1) { inds[n] = top;    vals[n++] = g; }

	    // diagonal
	    Scalar z = a;
	    if (IsBoundary3d(center, nx, ny, nz) && !isDirichlet) {
	      // Neumann boundary unknown (diagonal = sum of all offdiagonal)
	      z = Teuchos::ScalarTraits<Scalar>::zero();
	      for (size_t j = 0; j < n; j++)
		z -= vals[j];
	    }
	    inds[n]   = center;
	    vals[n++] = z;
	  }

	  for (size_t j = 0; j < n; j++)
	    inds[j] += indexBase;

	  Teuchos::ArrayView<GlobalOrdinal> iv(&inds[0], n);
	  Teuchos::ArrayView<Scalar>        av(&vals[0], n);
	  mtx->insertGlobalValues(myGlobalElements[i], iv, av);
	}
      }

      {
	Teuchos::RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Galeri: Laplace 3D fillComplete")));
	mtx->fillComplete();
      }

      return mtx;
    }

    /* ****************************************************************************************************** *
     *    3D 27-point stencil
     * ****************************************************************************************************** */
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
    Teuchos::RCP<Matrix>
    Brick3D(const Teuchos::RCP<const Map> & map,
            const GlobalOrdinal nx, const GlobalOrdinal ny, const GlobalOrdinal nz,
            const Scalar a, const Scalar b, const Scalar c,
            const Scalar d, const Scalar e,
            const Scalar f, const Scalar g,
            const DirBC DirichletBC = 0, const bool keepBCs = false)
    {
      using Teuchos::TimeMonitor;
      using Teuchos::RCP;
      using Teuchos::rcp;

      LocalOrdinal nnz = 27;

      Teuchos::RCP<Matrix> mtx = MatrixTraits<Map,Matrix>::Build(map, nnz);

      LocalOrdinal  numMyElements = map->getLocalNumElements();
      GlobalOrdinal indexBase     = map->getIndexBase();

      Teuchos::ArrayView<const GlobalOrdinal> myGlobalElements = map->getLocalElementList();

      GlobalOrdinal center, left, right, front, back, below, above;
      std::vector<Scalar>        vals(nnz);
      std::vector<GlobalOrdinal> inds(nnz);

      // upper plane
      //   e  d  e
      //   d  b  d
      //   e  d  e

      // middle plane
      //   c  b  c
      //   b  a  b
      //   c  b  c

      // lower plane
      //   e  d  e
      //   d  b  d
      //   e  d  e

      {
	Teuchos::RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Galeri: 3D 27 point stencil generation")));
	for (LocalOrdinal i = 0; i < numMyElements; ++i) {
	  size_t n = 0;

	  center = myGlobalElements[i] - indexBase;
	  GetNeighboursCartesian3d(center, nx, ny, nz,
				   left, right, front, back, below, above);

	  bool isDirichlet = (left  == -1 && (DirichletBC & DIR_LEFT))   ||
                             (right == -1 && (DirichletBC & DIR_RIGHT))  ||
                             (below == -1 && (DirichletBC & DIR_BOTTOM)) ||
                             (above == -1 && (DirichletBC & DIR_TOP))    ||
                             (front == -1 && (DirichletBC & DIR_FRONT))  ||
                             (back  == -1 && (DirichletBC & DIR_BACK));

	  if (isDirichlet && keepBCs) {
	    // Dirichlet unknown we want to keep
	    inds[n]   = center;
	    vals[n++] = Teuchos::ScalarTraits<Scalar>::one();

	  } else {
	    // See comments about weird in Cross2D

	    // center plance (centered on center)
	    Fill9PointStencil(center, vals, inds,
			      n, nx, ny, nz,
			      b, b, b, b, c, c, c, c,
			      left, right, front, back);
	    // lower plane (centered on "below")
	    if (below != -1) {
	      inds[n]   = below;
	      vals[n++] = b;
	      Fill9PointStencil(below, vals, inds, n, nx, ny, nz,
				d, d, d, d, e, e, e, e);
	    }
	    // upper plane (centered on "upper")
	    if (above != -1) {
	      inds[n]   = above;
	      vals[n++] = b;
	      Fill9PointStencil(above, vals, inds, n, nx, ny, nz,
				d, d, d, d, e, e, e, e);
	    }

	    // diagonal
	    Scalar z = a;
	    if (IsBoundary3d(center, nx, ny, nz) && !isDirichlet) {
	      // Neumann boundary unknown (diagonal = sum of all offdiagonal)
	      z = Teuchos::ScalarTraits<Scalar>::zero();
	      for (size_t j = 0; j < n; j++)
		z -= vals[j];
	    }
	    inds[n]   = center;
	    vals[n++] = z;
	  }

	  for (size_t j = 0; j < n; j++)
	    inds[j] += indexBase;

	  Teuchos::ArrayView<GlobalOrdinal> iv(&inds[0], n);
	  Teuchos::ArrayView<Scalar>        av(&vals[0], n);
	  mtx->insertGlobalValues(myGlobalElements[i], iv, av);
	}
      }

      {
	Teuchos::RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Galeri: 3D 27 point stencil fillComplete")));
	mtx->fillComplete();
      }
      return mtx;
    }

    /* ****************************************************************************************************** *
     *    Utilities
     * ****************************************************************************************************** */

    /* IsBoundary2d */
    template  <typename GlobalOrdinal>
    bool IsBoundary2d(const GlobalOrdinal i, const GlobalOrdinal nx, const GlobalOrdinal ny) {
      GlobalOrdinal ix, iy;
      ix = i % nx;
      iy = (i - ix) / nx;

      return (ix == 0 || ix == nx-1 || iy == 0 || iy == ny-1);
    }

    /* IsBoundary3d */
    template  <typename GlobalOrdinal>
    bool IsBoundary3d(const GlobalOrdinal i, const GlobalOrdinal nx, const GlobalOrdinal ny, const GlobalOrdinal nz) {
      GlobalOrdinal ix, iy, ixy, iz;
      ix  = i % nx;
      ixy = i % (nx * ny);
      iy  = (ixy - ix) / nx;
      iz  = (i - ixy) / (nx * ny);

      return (ix == 0 || ix == nx-1 || iy == 0 || iy == ny-1 || iz == 0 || iz == nz-1);
    }

    /* GetNeighboursCartesian2d */
    template <typename GlobalOrdinal>
    void GetNeighboursCartesian2d(const GlobalOrdinal i, const GlobalOrdinal nx, const GlobalOrdinal ny,
                                  GlobalOrdinal& left,  GlobalOrdinal& right,
                                  GlobalOrdinal& lower, GlobalOrdinal& upper)
    {
      GlobalOrdinal ix, iy;
      ix = i % nx;
      iy = (i - ix) / nx;

      if (ix == 0)      left  = -1;
      else              left  = i - 1;
      if (ix == nx - 1) right = -1;
      else              right = i + 1;
      if (iy == 0)      lower = -1;
      else              lower = i - nx;
      if (iy == ny - 1) upper = -1;
      else              upper = i + nx;
    }

    /* GetNeighboursCartesian2d */
    template <typename GlobalOrdinal>
    void GetNeighboursCartesian2d(const GlobalOrdinal i, const GlobalOrdinal nx, const GlobalOrdinal ny,
                                  GlobalOrdinal& left,  GlobalOrdinal& right,  GlobalOrdinal& lower,  GlobalOrdinal& upper,
                                  GlobalOrdinal& left2, GlobalOrdinal& right2, GlobalOrdinal& lower2, GlobalOrdinal& upper2)
    {
      GlobalOrdinal ix, iy;
      ix = i % nx;
      iy = (i - ix) / nx;

      if (ix == 0)      left  = -1;
      else              left  = i - 1;
      if (ix == nx - 1) right = -1;
      else              right = i + 1;
      if (iy == 0)      lower = -1;
      else              lower = i - nx;
      if (iy == ny - 1) upper = -1;
      else              upper = i + nx;

      if (ix <= 1)      left2  = -1;
      else              left2  = i - 2;
      if (ix >= nx - 2) right2 = -1;
      else              right2 = i + 2;
      if (iy <= 1)      lower2 = -1;
      else              lower2 = i - 2 * nx;
      if (iy >= ny - 2) upper2 = -1;
      else              upper2 = i + 2 * nx;
    }

    /* GetNeighboursCartesian3d */
    template <typename GlobalOrdinal>
    void GetNeighboursCartesian3d(const GlobalOrdinal i,
                                  const GlobalOrdinal nx, const GlobalOrdinal ny, const GlobalOrdinal nz,
                                  GlobalOrdinal& left,   GlobalOrdinal& right,
                                  GlobalOrdinal& front,  GlobalOrdinal& back,
                                  GlobalOrdinal& bottom, GlobalOrdinal& top)
    {
      GlobalOrdinal ixy, iz;
      ixy = i % (nx * ny);

      iz = (i - ixy) / (nx * ny);

      if (iz == 0)      bottom = -1;
      else              bottom = i - nx * ny;
      if (iz == nz - 1) top    = -1;
      else              top    = i + nx * ny;

      GetNeighboursCartesian2d(ixy, nx, ny, left, right, front, back);

      if (left  != -1) left  += iz * (nx * ny);
      if (right != -1) right += iz * (nx * ny);
      if (front != -1) front += iz * (nx * ny);
      if (back  != -1) back  += iz * (nx * ny);
    }

    /* Fill9PointStencil */
    template <typename GlobalOrdinal, typename Scalar>
    void Fill9PointStencil(const GlobalOrdinal center,
                           std::vector<Scalar>& vals, std::vector<GlobalOrdinal>& inds, size_t& n,
                           const GlobalOrdinal nx, const GlobalOrdinal ny, const GlobalOrdinal nz,
                           const Scalar b,  const Scalar c,  const Scalar d,  const Scalar e,
                           const Scalar z1, const Scalar z2, const Scalar z3, const Scalar z4,
                           GlobalOrdinal left,  GlobalOrdinal right,
                           GlobalOrdinal lower, GlobalOrdinal upper)
    {
      //  z3  e  z4
      //   b  .  c
      //  z1  d  z2
      GlobalOrdinal below, above;
      if (left == -2)
        GetNeighboursCartesian3d(center, nx, ny, nz, left, right, lower, upper, below, above);

      if (left  != -1)                { inds[n] = left;      vals[n++] = b;  }
      if (right != -1)                { inds[n] = right;     vals[n++] = c;  }
      if (lower != -1)                { inds[n] = lower;     vals[n++] = d;  }
      if (upper != -1)                { inds[n] = upper;     vals[n++] = e;  }
      if (left  != -1 && lower != -1) { inds[n] = lower - 1; vals[n++] = z1; }
      if (right != -1 && lower != -1) { inds[n] = lower + 1; vals[n++] = z2; }
      if (left  != -1 && upper != -1) { inds[n] = upper - 1; vals[n++] = z3; }
      if (right != -1 && upper != -1) { inds[n] = upper + 1; vals[n++] = z4; }
    }

  } // namespace Xpetra

} // namespace Galeri


#endif //ifndef GALERI_XPETRAMATRIXTYPES_HPP
