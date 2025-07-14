// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// This software is a result of the research described in the report
//
//     "A comparison of algorithms for modal analysis in the absence
//     of a sparse direct method", P. Arbenz, R. Lehoucq, and U. Hetmaniuk,
//     Sandia National Laboratories, Technical report SAND2003-1028J.
//

#ifndef ANASAZI_TPETRA_MODE_LAPLACE_1D_Q1_H
#define ANASAZI_TPETRA_MODE_LAPLACE_1D_Q1_H

#include "Anasazitpetra_ModeLaplace_DLLExportMacro.h"

#include <Teuchos_RCP.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_MultiVector.hpp>

#include "TpetraCheckingTools.hpp"
#include "TpetraSortingTools.hpp"


template<class Scalar, class LO, class GO, class Node>
class ANASAZITPETRA_MODELAPLACE_LIB_DLL_EXPORT ModeLaplace1DQ1 {

  private:

    const TpetraCheckingTools<Scalar,LO,GO,Node> myVerify;
    Teuchos::RCP<const Teuchos::Comm<int> > myComm;
    const TpetraSortingTools<Scalar> mySort;

    Teuchos::RCP<Tpetra::Map<LO,GO,Node> > Map;
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LO,GO,Node> > K;
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LO,GO,Node> > M;

    Scalar Lx;
    GO nX;

    std::vector<Scalar> x;

    static const int dofEle;
    static const int maxConnect;
#ifndef M_PI
    static const double M_PI;
#endif

    // Private member functions
    void preProcess();
    LO countElements(std::vector<bool>& isTouched);
    void makeMyElementsTopology(std::vector<GO>& elemTopo, const std::vector<bool>& isTouched);
    void makeMyConnectivity(std::vector<GO>& elemTopo, LO numEle, std::vector<GO>& connectivity, std::vector<size_t>& numNz);
    void makeStiffness(std::vector<GO>& elemTopo, LO numEle, std::vector<GO>& connectivity, const std::vector<size_t>& numNz);
    void makeMass(std::vector<GO>& elemTopo, LO numEle, std::vector<GO>& connectivity, const std::vector<size_t>& numNz);

    // Don't define these functions
    ModeLaplace1DQ1(const ModeLaplace1DQ1<Scalar,LO,GO,Node> &ref);
    ModeLaplace1DQ1& operator=(const ModeLaplace1DQ1<Scalar,LO,GO,Node> &ref);

  public:

    ModeLaplace1DQ1(Teuchos::RCP<const Teuchos::Comm<int> >& _Comm, Scalar _Lx, GO _nX)
    : myVerify(_Comm),
      myComm(_Comm),
      mySort(),
      Map(Teuchos::null),
      K(Teuchos::null),
      M(Teuchos::null),
      Lx(_Lx),
      nX(_nX)
    {
      preProcess();
    }

    ~ModeLaplace1DQ1() {}

    Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LO,GO,Node> > getStiffness() const { return K; }
    Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LO,GO,Node> > getMass()      const { return M; }
};

template<class Scalar, class LO, class GO, class Node>
void ModeLaplace1DQ1<Scalar, LO, GO, Node>::preProcess()
{
  GO globalSize = nX - 1;
  TEUCHOS_TEST_FOR_EXCEPTION(globalSize <= myComm->getSize(),std::logic_error,"Parameter error in ModeLaplace1DQ1.");

  // Create a uniform distribution of the unknowns across processors
  Map = Teuchos::rcp( new Tpetra::Map<LO,GO,Node>(globalSize, 0, myComm) );

  // Count the number of elements touched by this processor
  std::vector<bool> isTouched(nX, false);
  LO numEle = countElements(isTouched);

  // Create the mesh
  std::vector<GO> elemTopo(dofEle*numEle);
  makeMyElementsTopology(elemTopo, isTouched);

  // Get the number of nonzeros per row
  LO localSize = Map->getLocalNumElements();
  std::vector<size_t> numNz(localSize);
  std::vector<GO> connectivity(localSize*maxConnect);
  makeMyConnectivity(elemTopo, numEle, connectivity, numNz);

  // Make the stiffness matrix
  makeStiffness(elemTopo, numEle, connectivity, numNz);

  // Assemble the mass matrix
  makeMass(elemTopo, numEle, connectivity, numNz);

  // Get the geometrical coordinates of the managed nodes
  Scalar hx = Lx/nX;
  x.resize(localSize);

  globalSize = Map->getGlobalNumElements();
  for (GO i=0; i<globalSize; ++i) {
    if (Map->getLocalElement(i) > -1) {
      x[Map->getLocalElement(i)] = (i+1)*hx;
    }
  }
}

template<class Scalar, class LO, class GO, class Node>
LO ModeLaplace1DQ1<Scalar, LO, GO, Node>::countElements(std::vector<bool>& isTouched) 
{

  // This routine counts and flags the elements that contain the nodes
  // on this processor.

  GO i;
  LO numEle = 0;

  for (i=0; i<nX; ++i) {
    GO node;
    node = (i==0)  ? -1 : i-1;
    if ((node > -1) && (Map->getLocalElement(node) > -1)) {
      isTouched[i] = true;
      numEle++;
      continue;
    }
    node = (i==nX-1) ? -1 : i;
    if ((node > -1) && (Map->getLocalElement(node) > -1)) {
      isTouched[i] = true;
      numEle++;
      continue;
    }
  }

  return numEle;
}

template<class Scalar, class LO, class GO, class Node>
void ModeLaplace1DQ1<Scalar,LO,GO,Node>::makeMyElementsTopology(std::vector<GO>& elemTopo, const std::vector<bool>& isTouched) 
{

  // Create the element topology for the elements containing nodes for this processor
  // Note: Put the flag -1 when the node has a Dirichlet boundary condition

  LO i;
  LO numEle = 0;

  for (i=0; i<nX; ++i) {
    if (isTouched[i] == false)
      continue;
    elemTopo[dofEle*numEle]   = (i==0)  ? -1 : i-1;
    elemTopo[dofEle*numEle+1] = (i==nX-1) ? -1 : i;
    numEle++;
  }
}

template<class Scalar, class LO, class GO, class Node>
void ModeLaplace1DQ1<Scalar,LO,GO,Node>::makeMyConnectivity(std::vector<GO>& elemTopo, LO numEle, std::vector<GO>& connectivity, std::vector<size_t>& numNz) 
{
  // This routine creates the connectivity of each managed node
  // from the element topology.

  LO i, j;
  LO localSize = Map->getLocalNumElements();

  for (i=0; i<localSize; ++i) {
    numNz[i] = 0;
    for (j=0; j<maxConnect; ++j) {
      connectivity[i*maxConnect + j] = -1;
    }
  }

  for (i=0; i<numEle; ++i) {
    for (j=0; j<dofEle; ++j) {
      if (elemTopo[dofEle*i+j] == -1)
        continue;
      int node = Map->getLocalElement(elemTopo[dofEle*i+j]);
      if (node > -1) {
        int k;
        for (k=0; k<dofEle; ++k) {
          GO neighbor = elemTopo[dofEle*i+k];
          if (neighbor == -1)
            continue;
          // Check if this neighbor is already stored
          int l;
          for (l=0; l<maxConnect; ++l) {
            if (neighbor == connectivity[node*maxConnect + l])
              break;
            if (connectivity[node*maxConnect + l] == -1) {
              connectivity[node*maxConnect + l] = neighbor;
              numNz[node]++;
              break;
            }
          } // for (l = 0; l < maxConnect; ++l)
        } // for (k = 0; k < dofEle; ++k)
      } // if (node > -1)
    } // for (j = 0; j < dofEle; ++j)
  } // for (i = 0; i < numEle; ++i)
}

template<class Scalar, class LO, class GO, class Node>
void ModeLaplace1DQ1<Scalar,LO,GO,Node>::makeStiffness(std::vector<GO>& elemTopo, LO numEle, std::vector<GO>& connectivity, const std::vector<size_t>& numNz) 
{
  // Create Tpetra::CrsMatrix for stiffness
  K = Teuchos::rcp( new Tpetra::CrsMatrix<Scalar,LO,GO,Node>(Map, Teuchos::ArrayView<const size_t>(numNz)) );

  LO i;
  LO localSize = Map->getLocalNumElements();

  std::vector<Scalar> values(maxConnect, Teuchos::ScalarTraits<Scalar>::zero());

  for (i=0; i<localSize; ++i) 
    K->insertGlobalValues(Map->getGlobalElement(i), numNz[i], values.data(), connectivity.data()+(maxConnect*i));

  // Define the elementary matrix
  double hx = Lx/nX;
  std::vector<double> kel(dofEle*dofEle);
  kel[0] = 1.0/hx; kel[1] = -1.0/hx;
  kel[2] = -1.0/hx; kel[3] = 1.0/hx;

  // Assemble the matrix
  std::vector<GO> indices(dofEle);
  LO numEntries;
  LO j;
  for (i=0; i<numEle; ++i) {
    for (j=0; j<dofEle; ++j) {
      if (elemTopo[dofEle*i + j] == -1)
        continue;
      if (Map->getLocalElement(elemTopo[dofEle*i+j]) == -1)
        continue;
      numEntries = 0;
      LO k;
      for (k=0; k<dofEle; ++k) {
        if (elemTopo[dofEle*i+k] == -1)
          continue;
        indices[numEntries] = elemTopo[dofEle*i+k];
        values[numEntries] = kel[dofEle*j + k];
        numEntries++;
      }
      K->sumIntoGlobalValues(elemTopo[dofEle*i+j], numEntries, values.data(), indices.data());
    }
  }

  K->fillComplete();
}

template<class Scalar, class LO, class GO, class Node>
void ModeLaplace1DQ1<Scalar,LO,GO,Node>::makeMass(std::vector<GO>& elemTopo, LO numEle, std::vector<GO>& connectivity, const std::vector<size_t>& numNz) 
{
  // Create Tpetra::CrsMatrix for mass
  M = Teuchos::rcp( new Tpetra::CrsMatrix<Scalar,LO,GO,Node>(Map, Teuchos::ArrayView<const size_t>(numNz)) );

  LO i;
  LO localSize = Map->getLocalNumElements();

  std::vector<Scalar> values(maxConnect, Teuchos::ScalarTraits<Scalar>::zero());
  for (i=0; i<maxConnect; ++i)
    values[i] = 0.0;
  for (i=0; i<localSize; ++i)
    M->insertGlobalValues(Map->getGlobalElement(i), numNz[i], values.data(), connectivity.data()+(maxConnect*i));

  // Define the elementary matrix
  Scalar hx = Lx/nX;

  std::vector<Scalar> mel(dofEle*dofEle);
  mel[0] = hx/3.0; mel[1] = hx/6.0;
  mel[2] = hx/6.0; mel[3] = hx/3.0;

  // Assemble the matrix
  std::vector<GO> indices(dofEle);
  LO numEntries;
  LO j;
  for (i=0; i<numEle; ++i) {
    for (j=0; j<dofEle; ++j) {
      if (elemTopo[dofEle*i + j] == -1)
        continue;
      if (Map->getLocalElement(elemTopo[dofEle*i+j]) == -1)
        continue;
      numEntries = 0;
      LO k;
      for (k=0; k<dofEle; ++k) {
        if (elemTopo[dofEle*i+k] == -1)
          continue;
        indices[numEntries] = elemTopo[dofEle*i+k];
        values[numEntries] = mel[dofEle*j + k];
        numEntries++;
      }
      M->sumIntoGlobalValues(elemTopo[dofEle*i+j], numEntries, values.data(), indices.data());
    }
  }

  M->fillComplete();
}

#endif
