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

namespace Galeri {
  
  namespace Xpetra {
    
    /* prototypes */
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
                                  GlobalOrdinal& lower, GlobalOrdinal& upper,
                                  GlobalOrdinal& below, GlobalOrdinal& above);

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

      Teuchos::RCP<const Teuchos::Comm<int> > comm = map->getComm();

      GlobalOrdinal NumGlobalElements = map->getGlobalNumElements();

      GlobalOrdinal NumEntries;
      LocalOrdinal nnz=2;
      std::vector<Scalar> Values(nnz);
      std::vector<GlobalOrdinal> Indices(nnz);

      comm->barrier();
      if (comm->getRank() == 0) {
        std::cout << "starting global insert" << std::endl;
	//        std::cout << MemUtils::PrintMemoryUsage() << std::endl;
      }

/*
      double t0 = MPI_Wtime();
      double t1,t2;
*/

      Teuchos::RCP<Teuchos::Time> timer = rcp(new Teuchos::Time("TriDiag global insert"));
      timer->start(true);

      for (LocalOrdinal i = 0; i < NumMyElements; ++i)
        {
          if (MyGlobalElements[i] == 0)
            {
              // off-diagonal for first row
              Indices[0] = 1;
              NumEntries = 1;
              Values[0] = c;
            }
          else if (MyGlobalElements[i] == NumGlobalElements - 1)
            {
              // off-diagonal for last row
              Indices[0] = NumGlobalElements - 2;
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

/*
        if ( (comm->getRank() == 0) && (NumMyElements >= 10) &&  (i % (NumMyElements / 10) == 0) ) {
            int percDone = (int) floor((((double)i)/NumMyElements)*100);
            t1 = MPI_Wtime() - t0;
            std::cout << percDone << "% done (" << i << " rows) in " << t1 << " seconds, [pid 0]" << std::endl;
            std::cout << MemUtils::PrintMemoryUsage() << std::endl;
          }
*/
        } //for (LocalOrdinal i = 0; i < NumMyElements; ++i)

        timer->stop();
        //MemUtils::ReportTimeAndMemory(*timer, *comm);

/*
      t1 = MPI_Wtime() - t0;
      if (comm->getRank() == 0) {
        std::cout << "100% done in " << t1 << " seconds [pid 0]" << std::endl;
        std::cout << MemUtils::PrintMemoryUsage() << std::endl;
      }
      if (comm->getRank() == 0) {
        std::cout << "starting fill complete" << std::endl;
        std::cout << MemUtils::PrintMemoryUsage() << std::endl;
      }
*/

      //t2 = MPI_Wtime();
      timer = rcp(new Teuchos::Time("TriDiag fillComplete"));
      timer->start(true);

      mtx->fillComplete();

      timer->stop();
      //MemUtils::ReportTimeAndMemory(*timer, *comm);
      /*
      t2 = MPI_Wtime() - t2;
      t1 = MPI_Wtime() - t0;
      if (comm->getRank() == 0) {
        std::cout << "time to FillComplete = " << t2 << " seconds [pid 0]" << std::endl;
        std::cout << "total time = " << t1 << " seconds [pid 0]" << std::endl;
        std::cout << MemUtils::PrintMemoryUsage() << std::endl;
      }
      */

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
            const Scalar d, const Scalar e, const bool keepBCs=false)
    {

      Teuchos::RCP<Matrix> mtx = MatrixTraits<Map,Matrix>::Build(map, 5);

      LocalOrdinal NumMyElements = map->getNodeNumElements();
      Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getNodeElementList();

      GlobalOrdinal left, right, lower, upper;
      LocalOrdinal nnz=5;
      std::vector<Scalar> Values(nnz);
      std::vector<GlobalOrdinal> Indices(nnz);

      //    e
      //  b a c
      //    d
      
      if (keepBCs) {
        for (LocalOrdinal i = 0; i < NumMyElements; ++i)  {
          GlobalOrdinal NumEntries = 0;
          GetNeighboursCartesian2d(MyGlobalElements[i], nx, ny, 
                                   left, right, lower, upper);
          if (left == -1 || right == -1 || lower == -1 || upper == -1) {
            //Dirichlet unknown
            mtx->insertGlobalValues(MyGlobalElements[i],
                                    Teuchos::tuple<GlobalOrdinal>(MyGlobalElements[i]),
                                    Teuchos::tuple<Scalar>(Teuchos::ScalarTraits<Scalar>::one()) );
          } else {
            Indices[0] = left;
            Values[0] = b;
            Indices[1] = right;
            Values[1] = c;
            Indices[2] = lower;
            Values[2] = d;
            Indices[3] = upper;
            Values[3] = e;
            Indices[4] = MyGlobalElements[i]; //diagonal
            Values[4] = a;
            NumEntries=5;
            Teuchos::ArrayView<Scalar> av(&Values[0],NumEntries);
            Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0],NumEntries);
            mtx->insertGlobalValues(MyGlobalElements[i], iv, av);
          }
            
        }
 
      } else {
        for (LocalOrdinal i = 0; i < NumMyElements; ++i) 
        {
          GlobalOrdinal NumEntries = 0;
          GetNeighboursCartesian2d(MyGlobalElements[i], nx, ny, 
                                   left, right, lower, upper);

          if (left != -1) {
             Indices[NumEntries] = left;
             Values[NumEntries] = b;
             ++NumEntries;
          }
          if (right != -1) {
             Indices[NumEntries] = right;
             Values[NumEntries] = c;
             ++NumEntries;
          }
          if (lower != -1) {
             Indices[NumEntries] = lower;
             Values[NumEntries] = d;
             ++NumEntries;
          }
          if (upper != -1) {
             Indices[NumEntries] = upper;
             Values[NumEntries] = e;
             ++NumEntries;
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
        } //for (LocalOrdinal i=0...
      } //if (keepBCs)

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
     *    Laplace3D
     * ****************************************************************************************************** */
    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
    Teuchos::RCP<Matrix>
    Cross3D(const Teuchos::RCP<const Map> & map,
            const GlobalOrdinal nx, const GlobalOrdinal ny, const GlobalOrdinal nz,
            const Scalar a, const Scalar b, const Scalar c, 
            const Scalar d, const Scalar e,
            const Scalar f, const Scalar g)
    {

      Teuchos::RCP<Matrix> mtx = MatrixTraits<Map,Matrix>::Build(map, 7);

      LocalOrdinal NumMyElements = map->getNodeNumElements();
      Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getNodeElementList();

      GlobalOrdinal left, right, lower, upper, below, above;
      std::vector<Scalar> Values(6);
      std::vector<GlobalOrdinal> Indices(6);

      //    e
      //  b a c
      //    d
      // + f below and g above
  
      for (GlobalOrdinal i = 0; i < NumMyElements; ++i) 
        {
          GlobalOrdinal NumEntries = 0;
          GetNeighboursCartesian3d(MyGlobalElements[i], nx, ny, nz,
                                   left, right, lower, upper, below, above);

          if (left != -1) 
            {
              Indices[NumEntries] = left;
              Values[NumEntries] = b;
              ++NumEntries;
            }
          if (right != -1) 
            {
              Indices[NumEntries] = right;
              Values[NumEntries] = c;
              ++NumEntries;
            }
          if (lower != -1) 
            {
              Indices[NumEntries] = lower;
              Values[NumEntries] = d;
              ++NumEntries;
            }
          if (upper != -1) 
            {
              Indices[NumEntries] = upper;
              Values[NumEntries] = e;
              ++NumEntries;
            }
          if (below != -1) 
            {
              Indices[NumEntries] = below;
              Values[NumEntries] = f;
              ++NumEntries;
            }
          if (above != -1) 
            {
              Indices[NumEntries] = above;
              Values[NumEntries] = g;
              ++NumEntries;
            }
          // put the off-diagonal entries
          Teuchos::ArrayView<Scalar> av(&Values[0],NumEntries);
          Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0],NumEntries);
          mtx->insertGlobalValues(MyGlobalElements[i], iv, av);

          // Put in the diagonal entry
          mtx->insertGlobalValues(MyGlobalElements[i],
                                  Teuchos::tuple<GlobalOrdinal>(MyGlobalElements[i]),
                                  Teuchos::tuple<Scalar>(a) );
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
                                  GlobalOrdinal& lower, GlobalOrdinal& upper,
                                  GlobalOrdinal& below, GlobalOrdinal& above)
    {
      GlobalOrdinal ixy, iz;
      ixy = i % (nx * ny);

      iz = (i - ixy) / (nx * ny); 

      if (iz == 0)      below = -1;
      else              below = i - nx * ny;
      if (iz == nz - 1) above = -1;
      else              above = i + nx * ny;

      GetNeighboursCartesian2d(ixy, nx, ny, left, right, lower, upper);

      if (left  != -1) left  += iz * (nx * ny);
      if (right != -1) right += iz * (nx * ny);
      if (lower != -1) lower += iz * (nx * ny);
      if (upper != -1) upper += iz * (nx * ny);
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
