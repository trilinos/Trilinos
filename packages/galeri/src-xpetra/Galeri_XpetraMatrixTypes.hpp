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

#ifndef GALERI_XPETRAMATRIXTYPES_HPP
#define GALERI_XPETRAMATRIXTYPES_HPP
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayView.hpp>

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
                                  GlobalOrdinal & left, GlobalOrdinal & right,
                                  GlobalOrdinal & lower, GlobalOrdinal & upper);

    template <typename GlobalOrdinal>
    void GetNeighboursCartesian2d(const GlobalOrdinal i, const GlobalOrdinal nx, const GlobalOrdinal ny,
                                  GlobalOrdinal& left, GlobalOrdinal& right, GlobalOrdinal& lower, GlobalOrdinal& upper,
                                  GlobalOrdinal& left2, GlobalOrdinal& right2, GlobalOrdinal& lower2, GlobalOrdinal& upper2);

    template <typename GlobalOrdinal>
    void GetNeighboursCartesian3d(const GlobalOrdinal i,
                                  const GlobalOrdinal nx, const GlobalOrdinal ny, const GlobalOrdinal nz,
                                  GlobalOrdinal& left, GlobalOrdinal& right,
                                  GlobalOrdinal& front, GlobalOrdinal& back,
                                  GlobalOrdinal& bottom, GlobalOrdinal& top);

    template <typename GlobalOrdinal, typename Scalar>
    void Fill9PointStencil(const GlobalOrdinal center,
                           std::vector<Scalar>  &Values,
                           std::vector<GlobalOrdinal>  &Indices,
                           GlobalOrdinal &NumEntries,
                           const GlobalOrdinal nx, const GlobalOrdinal ny,
                           const GlobalOrdinal nz,
                           const Scalar b,  const Scalar c,  const Scalar d,  const Scalar e,
                           const Scalar z1, const Scalar z2, const Scalar z3, const Scalar z4,
                           GlobalOrdinal left=-2, GlobalOrdinal right=-2,
                           GlobalOrdinal lower=-2, GlobalOrdinal upper=-2);
    /* end of prototypes */

    /* ****************************************************************************************************** *
     *    (Scaled) Identity
     * ****************************************************************************************************** */
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
    Teuchos::RCP<Matrix>
    Identity(const Teuchos::RCP<const Map> & map,
            const Scalar a)
    {

      Teuchos::RCP<Matrix> mtx = MatrixTraits<Map,Matrix>::Build(map, 1);

      LocalOrdinal NumMyElements = map->getNodeNumElements();
      Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getNodeElementList();

      for (LocalOrdinal i = 0; i < NumMyElements; ++i) {
          mtx->insertGlobalValues(MyGlobalElements[i],
                                  Teuchos::tuple<GlobalOrdinal>(MyGlobalElements[i]),
                                  Teuchos::tuple<Scalar>(a) );
      }

      mtx->fillComplete();

      return mtx;
    } //Identity

    /* ****************************************************************************************************** *
     *    Laplace 1D
     * ****************************************************************************************************** */
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
    Teuchos::RCP<Matrix>
    TriDiag(const Teuchos::RCP<const Map> & map,
            const GlobalOrdinal nx, // note: nx unused
            const Scalar a, const Scalar b, const Scalar c)
    {

      Teuchos::RCP<Matrix> mtx = MatrixTraits<Map,Matrix>::Build(map, 3);

      LocalOrdinal NumMyElements = map->getNodeNumElements();
      Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getNodeElementList();
      GlobalOrdinal indexBase = map->getIndexBase();

      Teuchos::RCP<const Teuchos::Comm<int> > comm = map->getComm();

      GlobalOrdinal NumGlobalElements = map->getGlobalNumElements();

      GlobalOrdinal NumEntries;
      LocalOrdinal nnz=2;
      std::vector<Scalar> Values(nnz);
      std::vector<GlobalOrdinal> Indices(nnz);

      comm->barrier();
      if (comm->getRank() == 0) {
        std::cout << "starting global insert" << std::endl;
      }

      Teuchos::RCP<Teuchos::Time> timer = rcp(new Teuchos::Time("TriDiag global insert"));
      timer->start(true);

      for (LocalOrdinal i = 0; i < NumMyElements; ++i)
        {
          if (MyGlobalElements[i] == indexBase)
            {
              // off-diagonal for first row
              Indices[0] = 1 + indexBase;
              NumEntries = 1;
              Values[0] = c;
            }
          else if (MyGlobalElements[i] == NumGlobalElements + indexBase - 1)
            {
              // off-diagonal for last row
              Indices[0] = NumGlobalElements - 2 + indexBase;
              NumEntries = 1;
              Values[0] = b;
            }
          else
            {
              // off-diagonal for internal row
              Indices[0] = MyGlobalElements[i] - 1;
              Values[1] = b;
              Indices[1] = MyGlobalElements[i] + 1;
              Values[0] = c;
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
                                  Teuchos::tuple<Scalar>(a) );

        } //for (LocalOrdinal i = 0; i < NumMyElements; ++i)

        timer->stop();

      timer = rcp(new Teuchos::Time("TriDiag fillComplete"));
      timer->start(true);

      mtx->fillComplete();

      timer->stop();

      return mtx;
    } //TriDiag

    /* ****************************************************************************************************** *
     *    Laplace 2D
     * ****************************************************************************************************** */
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
    Teuchos::RCP<Matrix>
    Cross2D(const Teuchos::RCP<const Map> & map,
            const GlobalOrdinal nx, const GlobalOrdinal ny,
            const Scalar a, const Scalar b, const Scalar c,
            const Scalar d, const Scalar e, const DirBC DirichletBC = 0, const bool keepBCs=false)
    {

      Teuchos::RCP<Matrix> mtx = MatrixTraits<Map,Matrix>::Build(map, 5);

      LocalOrdinal NumMyElements = map->getNodeNumElements();
      Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getNodeElementList();
      GlobalOrdinal indexBase = map->getIndexBase();

      GlobalOrdinal left, right, lower, upper, center;
      LocalOrdinal nnz=5;
      std::vector<Scalar> Values(nnz);
      std::vector<GlobalOrdinal> Indices(nnz);

      //    e
      //  b a c
      //    d

      for (LocalOrdinal i = 0; i < NumMyElements; ++i)  {
        size_t numEntries = 0;

        center = MyGlobalElements[i];
        //GetNeighboursCartesian2d is zero-based, so shift the center point to get the correct neighbors
        GetNeighboursCartesian2d(center-indexBase, nx, ny,
                                 left, right, lower, upper);

        bool isDirichlet = (left  == -1 && (DirichletBC & DIR_LEFT))   ||
                           (right == -1 && (DirichletBC & DIR_RIGHT))  ||
                           (lower == -1 && (DirichletBC & DIR_BOTTOM)) ||
                           (upper == -1 && (DirichletBC & DIR_TOP));

        if (isDirichlet && keepBCs) {
          // Dirichlet unknown we want to keep
          mtx->insertGlobalValues(center,
                                  Teuchos::tuple<GlobalOrdinal>(center),
                                  Teuchos::tuple<Scalar>(Teuchos::ScalarTraits<Scalar>::one()) );
        } else {
          // The Neumann b.c. are treated in a sane way. The Dirichlet b.c., however, are treated
          // insane when the option keepBCs=false. Speicifically, in this case we don't want to keep
          // Dirichlet b.c., but that would result in inconsistency between the map and the number of
          // degrees of freedom, plus the problem with GIDs. Therefore, we virtually expand domain by
          // one node in the direction of the Dirichlet b.c., and then assume that that node was
          // not kept. But we use an old GIDs. So yes, that's weird.

          if (left != -1) {
            Indices[numEntries] = left+indexBase;
            Values [numEntries] = b;
            numEntries++;
          }
          if (right != -1) {
            Indices[numEntries] = right+indexBase;
            Values [numEntries] = c;
            numEntries++;
          }
          if (lower != -1) {
            Indices[numEntries] = lower+indexBase;
            Values [numEntries] = d;
            numEntries++;
          }
          if (upper != -1) {
            Indices[numEntries] = upper+indexBase;
            Values [numEntries] = e;
            numEntries++;
          }
          // diagonal
          Scalar z = a;
          if (IsBoundary2d(center, nx, ny) && !isDirichlet) {
            // Neumann boundary unknown
            // Diagonal = sum of all offdiagonal
            z = Teuchos::ScalarTraits<Scalar>::zero();
            for (size_t j = 0; j < numEntries; j++)
              z -= Values[j];
          }
          Indices[numEntries] = center;
          Values [numEntries] = z;
          numEntries++;

          Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0], numEntries);
          Teuchos::ArrayView<Scalar>        av(&Values[0],  numEntries);
          mtx->insertGlobalValues(center, iv, av);
        }
      }

      mtx->fillComplete();

      return mtx;
    } //Cross2D

    /* ****************************************************************************************************** *
     *    Star2D
     * ****************************************************************************************************** */
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
    Teuchos::RCP<Matrix>
    Star2D(const Teuchos::RCP<const Map> & map,
           const GlobalOrdinal nx, const GlobalOrdinal ny,
           const Scalar a, const Scalar b, const Scalar c,
           const Scalar d, const Scalar e,
           const Scalar z1, const Scalar z2,
           const Scalar z3, const Scalar z4)
    {

      Teuchos::RCP<Matrix> mtx = MatrixTraits<Map,Matrix>::Build(map, 9);

      LocalOrdinal NumMyElements = map->getNodeNumElements();
      Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getNodeElementList();

      GlobalOrdinal left, right, lower, upper;
      //Scalar Values[9];
      //GlobalOrdinal    Indices[9];
      std::vector<Scalar> Values(9);
      std::vector<GlobalOrdinal>    Indices(9);

      //  z3  e  z4
      //   b  a  c
      //  z1  d  z2
      for (LocalOrdinal i = 0; i < NumMyElements; ++i)
        {
          GlobalOrdinal NumEntries = 0;
          left = right = lower = upper = -1;  //FIXME JJH: I don't remember why the next lines are commented out....
                                              //FIXME but this is just to get rid of an annoying compiler warning
          /*
            GetNeighboursCartesian2d(MyGlobalElements[i], nx, ny,
            left, right, lower, upper);
          */
#define OLD_CODE
#ifndef OLD_CODE
          Fill9PointStencil<GlobalOrdinal,Scalar>(MyGlobalElements[i], Values, Indices, NumEntries, nx, ny,
                                                  b, c, d, e, z1, z2, z3, z4);
#else
          if (left != -1)
            {
              Values[NumEntries] = b;
              Indices[NumEntries] = left;
              ++NumEntries;
            }
          if (right != -1)
            {
              Values[NumEntries] = c;
              Indices[NumEntries] = right;
              ++NumEntries;
            }
          if (lower != -1)
            {
              Values[NumEntries] = d;
              Indices[NumEntries] = lower;
              ++NumEntries;
            }
          if (upper != -1)
            {
              Values[NumEntries] = e;
              Indices[NumEntries] = upper;
              ++NumEntries;
            }
          if (left != -1 && lower != -1)
            {
              Values[NumEntries] = z1;
              Indices[NumEntries] = lower - 1;
              ++NumEntries;
            }
          if (right != -1 && lower != -1)
            {
              Values[NumEntries] = z2;
              Indices[NumEntries] = lower + 1;
              ++NumEntries;
            }
          if (left != -1 && upper != -1)
            {
              Values[NumEntries] = z3;
              Indices[NumEntries] = upper - 1;
              ++NumEntries;
            }
          if (right != -1 && upper != -1)
            {
              Values[NumEntries] = z4;
              Indices[NumEntries] = upper + 1;
              ++NumEntries;
            }
#endif

          Values[NumEntries] = a;
          Indices[NumEntries] = MyGlobalElements[i];
          ++NumEntries;

          Teuchos::ArrayView<Scalar> av(&Values[0],NumEntries);
          Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0],NumEntries);
          mtx->insertGlobalValues(MyGlobalElements[i], iv, av);
        }

      mtx->fillComplete();

      return mtx;
    } //Star2D

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
              const Scalar bb, const Scalar cc, const Scalar dd, const Scalar ee)
    {

      Teuchos::RCP<Matrix> mtx= MatrixTraits<Map,Matrix>::Build(map, 13);

      LocalOrdinal NumMyElements = map->getNodeNumElements();
      Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getNodeElementList();

      GlobalOrdinal left, right, lower, upper;
      GlobalOrdinal left2, right2, lower2, upper2;
      Scalar Values[13];
      GlobalOrdinal Indices[13];

      //        ee
      //    z3  e  z4
      // bb  b  a  c  cc
      //    z1  d  z2
      //        dd

      for (GlobalOrdinal i = 0; i < NumMyElements; ++i)
        {
          GlobalOrdinal NumEntries = 0;
          GetNeighboursCartesian2d(MyGlobalElements[i], nx, ny, left, right, lower, upper,
                                   left2, right2, lower2, upper2);

          if (left != -1)
            {
              Values[NumEntries] = b;
              Indices[NumEntries] = left;
              ++NumEntries;
            }
          if (right != -1)
            {
              Values[NumEntries] = c;
              Indices[NumEntries] = right;
              ++NumEntries;
            }
          if (lower != -1)
            {
              Values[NumEntries] = d;
              Indices[NumEntries] = lower;
              ++NumEntries;
            }
          if (upper != -1)
            {
              Values[NumEntries] = e;
              Indices[NumEntries] = upper;
              ++NumEntries;
            }
          if (left != -1 && lower != -1)
            {
              Values[NumEntries] = z1;
              Indices[NumEntries] = lower - 1;
              ++NumEntries;
            }
          if (right != -1 && lower != -1)
            {
              Values[NumEntries] = z2;
              Indices[NumEntries] = lower + 1;
              ++NumEntries;
            }
          if (left != -1 && upper != -1)
            {
              Values[NumEntries] = z3;
              Indices[NumEntries] = upper - 1;
              ++NumEntries;
            }
          if (right != -1 && upper != -1)
            {
              Values[NumEntries] = z4;
              Indices[NumEntries] = upper + 1;
              ++NumEntries;
            }
          if (left2 != -1)
            {
              Values[NumEntries] = bb;
              Indices[NumEntries] = left2;
              ++NumEntries;
            }
          if (right2 != -1)
            {
              Values[NumEntries] = cc;
              Indices[NumEntries] = right2;
              ++NumEntries;
            }
          if (lower2 != -1)
            {
              Values[NumEntries] = dd;
              Indices[NumEntries] = lower2;
              ++NumEntries;
            }
          if (upper2 != -1)
            {
              Values[NumEntries] = ee;
              Indices[NumEntries] = upper2;
              ++NumEntries;
            }

          Values[NumEntries] = a;
          Indices[NumEntries] = MyGlobalElements[i];
          ++NumEntries;

          Teuchos::ArrayView<Scalar> av(&Values[0],NumEntries);
          Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0],NumEntries);
          mtx->insertGlobalValues(MyGlobalElements[i], iv, av);
        }

      mtx->fillComplete();

      return mtx;
    } //BigStar2D

    /* ****************************************************************************************************** *
     *    Laplace 3D
     * ****************************************************************************************************** */
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
    Teuchos::RCP<Matrix>
    Cross3D(const Teuchos::RCP<const Map> & map,
            const GlobalOrdinal nx, const GlobalOrdinal ny, const GlobalOrdinal nz,
            const Scalar a, const Scalar b, const Scalar c,
            const Scalar d, const Scalar e,
            const Scalar f, const Scalar g,
            const DirBC DirichletBC = 0, const bool keepBCs = false)
    {

      Teuchos::RCP<Matrix> mtx = MatrixTraits<Map,Matrix>::Build(map, 7);

      LocalOrdinal                               NumMyElements = map->getNodeNumElements();
      Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getNodeElementList();

      GlobalOrdinal left, right, bottom, top, front, back, center;
      std::vector<Scalar> Values(7);
      std::vector<GlobalOrdinal> Indices(7);

      //    e
      //  b a c
      //    d
      // + f bottom and g top

      for (GlobalOrdinal i = 0; i < NumMyElements; ++i) {
        size_t numEntries = 0;

        center = MyGlobalElements[i];
        GetNeighboursCartesian3d(center, nx, ny, nz,
                                 left, right, front, back, bottom, top);

        bool isDirichlet = (left  == -1 && (DirichletBC & DIR_LEFT))   ||
                           (right == -1 && (DirichletBC & DIR_RIGHT))  ||
                           (front == -1 && (DirichletBC & DIR_BOTTOM)) ||
                           (back  == -1 && (DirichletBC & DIR_TOP))    ||
                           (front == -1 && (DirichletBC & DIR_FRONT))  ||
                           (back  == -1 && (DirichletBC & DIR_BACK));

        if (isDirichlet && keepBCs) {
          // Dirichlet unknown we want to keep
          mtx->insertGlobalValues(center,
                                  Teuchos::tuple<GlobalOrdinal>(center),
                                  Teuchos::tuple<Scalar>(Teuchos::ScalarTraits<Scalar>::one()) );
        } else {
          // See comments about weird in Cross2D

          if (left != -1) {
            Indices[numEntries] = left;
            Values [numEntries] = b;
            numEntries++;
          }
          if (right != -1) {
            Indices[numEntries] = right;
            Values [numEntries] = c;
            numEntries++;
          }
          if (front != -1) {
            Indices[numEntries] = front;
            Values [numEntries] = d;
            numEntries++;
          }
          if (back != -1) {
            Indices[numEntries] = back;
            Values [numEntries] = e;
            numEntries++;
          }
          if (bottom != -1) {
            Indices[numEntries] = bottom;
            Values [numEntries] = f;
            numEntries++;
          }
          if (top != -1) {
            Indices[numEntries] = top;
            Values [numEntries] = g;
            numEntries++;
          }
          // diagonal
          Scalar z = a;
          if (IsBoundary3d(center, nx, ny, nz) && !isDirichlet) {
            // Neumann boundary unknown
            // Diagonal = sum of all offdiagonal
            z = Teuchos::ScalarTraits<Scalar>::zero();
            for (size_t j = 0; j < numEntries; j++)
              z -= Values[j];
          }
          Indices[numEntries] = center;
          Values [numEntries] = z;
          numEntries++;

          Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0], numEntries);
          Teuchos::ArrayView<Scalar>        av(&Values[0],  numEntries);
          mtx->insertGlobalValues(center, iv, av);
        }
      }

      mtx->fillComplete();

      return mtx;
    } //Cross3D

    /* ****************************************************************************************************** *
     *    3D 27-point stencil
     * ****************************************************************************************************** */
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
    Teuchos::RCP<Matrix>
    Brick3D(const Teuchos::RCP<const Map> & map,
            const GlobalOrdinal nx, const GlobalOrdinal ny, const GlobalOrdinal nz,
            const Scalar a, const Scalar b, const Scalar c,
            const Scalar d, const Scalar e,
            const Scalar f, const Scalar g)
    {

      Teuchos::RCP<Matrix> mtx = MatrixTraits<Map,Matrix>::Build(map, 27);

      LocalOrdinal NumMyElements = map->getNodeNumElements();
      Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getNodeElementList();

      GlobalOrdinal left, right, lower, upper, below, above;
      std::vector<Scalar> Values(26);
      std::vector<GlobalOrdinal> Indices(26);

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

      for (GlobalOrdinal i = 0; i < NumMyElements; ++i)
        {
          GlobalOrdinal NumEntries = 0;
          // First process the middle plane and figure out the
          // centers of the upper and lower planes
          GetNeighboursCartesian3d(MyGlobalElements[i], nx, ny, nz,
                                   left, right, lower, upper, below, above);

#ifdef DEBUG_STENCIL
          std::cout << "middle plane, center " << MyGlobalElements[i]  << std::endl;
#endif
          Fill9PointStencil(MyGlobalElements[i],Values, Indices,
                            NumEntries, nx, ny, nz,
                            b, b, b, b, c, c, c, c,
                            left, right, lower, upper);

          //process lower plane, centered on "below"
          if (below != -1)
            {
              Indices[NumEntries] = below;
              Values[NumEntries] = b;
              ++NumEntries;
#ifdef DEBUG_STENCIL
              std::cout << "lower plane, center " << below  << std::endl;
#endif
              Fill9PointStencil(below, Values, Indices, NumEntries, nx, ny, nz,
                                d, d, d, d, e, e, e, e);
            }
          //process upper plane, centered on "upper"
          if (above != -1)
            {
              Indices[NumEntries] = above;
              Values[NumEntries] = b;
              ++NumEntries;
#ifdef DEBUG_STENCIL
              std::cout << "upper plane, center " << above  << std::endl;
#endif
              Fill9PointStencil(above, Values, Indices, NumEntries, nx, ny, nz,
                                d, d, d, d, e, e, e, e);
            }


          // put the off-diagonal entries
          Teuchos::ArrayView<Scalar> av(&Values[0],NumEntries);
          Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0],NumEntries);
          mtx->insertGlobalValues(MyGlobalElements[i], iv, av);

          // Put in the diagonal entry
          mtx->insertGlobalValues(MyGlobalElements[i],
                                  Teuchos::tuple<GlobalOrdinal>(MyGlobalElements[i]),
                                  Teuchos::tuple<Scalar>(a) );
        } //for (GlobalOrdinal i = 0; i < NumMyElements; ++i)
      mtx->fillComplete();

      return mtx;
    } //Brick3D

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
                                  GlobalOrdinal & left, GlobalOrdinal & right,
                                  GlobalOrdinal & lower, GlobalOrdinal & upper)
    {
      GlobalOrdinal ix, iy;
      ix = i % nx;
      iy = (i - ix) / nx;

      if (ix == 0)      left = -1;
      else              left = i - 1;
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
                                  GlobalOrdinal& left, GlobalOrdinal& right, GlobalOrdinal& lower, GlobalOrdinal& upper,
                                  GlobalOrdinal& left2, GlobalOrdinal& right2, GlobalOrdinal& lower2, GlobalOrdinal& upper2)
    {
      GlobalOrdinal ix, iy;
      ix = i % nx;
      iy = (i - ix) / nx;

      if (ix == 0)      left = -1;
      else              left = i - 1;
      if (ix == nx - 1) right = -1;
      else              right = i + 1;
      if (iy == 0)      lower = -1;
      else              lower = i - nx;
      if (iy == ny - 1) upper = -1;
      else              upper = i + nx;

      if (ix <= 1)      left2 = -1;
      else              left2 = i - 2;
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
                                  GlobalOrdinal& left, GlobalOrdinal& right,
                                  GlobalOrdinal& front, GlobalOrdinal& back,
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
                           std::vector<Scalar>  &Values,
                           std::vector<GlobalOrdinal>  &Indices,
                           GlobalOrdinal &NumEntries,
                           const GlobalOrdinal nx, const GlobalOrdinal ny, const GlobalOrdinal nz,
                           const Scalar b,  const Scalar c,  const Scalar d,  const Scalar e,
                           const Scalar z1, const Scalar z2, const Scalar z3, const Scalar z4,
                           GlobalOrdinal left, GlobalOrdinal right,
                           GlobalOrdinal lower, GlobalOrdinal upper)
    {
      //  z3  e  z4
      //   b  .  c
      //  z1  d  z2
      GlobalOrdinal below,above;
#ifdef DEBUG_STENCIL
      std::cout << "   ";
#endif
      if (left == -2) {
        GetNeighboursCartesian3d(center, nx, ny, nz, left, right, lower, upper, below, above);
      }

      if (left != -1)
        {
#ifdef DEBUG_STENCIL
          std::cout << left << " ";
#endif
          Values[NumEntries] = b;
          Indices[NumEntries] = left;
          ++NumEntries;
        }
      if (right != -1)
        {
#ifdef DEBUG_STENCIL
          std::cout << right << " ";
#endif
          Values[NumEntries] = c;
          Indices[NumEntries] = right;
          ++NumEntries;
        }
      if (lower != -1)
        {
#ifdef DEBUG_STENCIL
          std::cout << lower << " ";
#endif
          Values[NumEntries] = d;
          Indices[NumEntries] = lower;
          ++NumEntries;
        }
      if (upper != -1)
        {
#ifdef DEBUG_STENCIL
          std::cout << upper << " ";
#endif
          Values[NumEntries] = e;
          Indices[NumEntries] = upper;
          ++NumEntries;
        }
      if (left != -1 && lower != -1)
        {
#ifdef DEBUG_STENCIL
          std::cout << lower-1 << " ";
#endif
          Values[NumEntries] = z1;
          Indices[NumEntries] = lower - 1;
          ++NumEntries;
        }
      if (right != -1 && lower != -1)
        {
#ifdef DEBUG_STENCIL
          std::cout << lower+1 << " ";
#endif
          Values[NumEntries] = z2;
          Indices[NumEntries] = lower + 1;
          ++NumEntries;
        }
      if (left != -1 && upper != -1)
        {
#ifdef DEBUG_STENCIL
          std::cout << upper-1 << " ";
#endif
          Values[NumEntries] = z2;
          Values[NumEntries] = z3;
          Indices[NumEntries] = upper - 1;
          ++NumEntries;
        }
      if (right != -1 && upper != -1)
        {
#ifdef DEBUG_STENCIL
          std::cout << upper+1 << " ";
#endif
          Values[NumEntries] = z4;
          Indices[NumEntries] = upper + 1;
          ++NumEntries;
        }
#ifdef DEBUG_STENCIL
      std::cout << std::endl;
#endif
    } //Fill9PointStencil()


  } // namespace Xpetra
} // namespace Galeri


#endif //ifndef GALERI_XPETRAMATRIXTYPES_HPP
