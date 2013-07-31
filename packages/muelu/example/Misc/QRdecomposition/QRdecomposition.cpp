/*
 * QRdecomposition.cpp
 *
 *  Created on: Jan 31, 2013
 *      Author: tobias
 */




#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_FancyOStream.hpp>

// MueLu
#include "MueLu_QR_InterfaceEx.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

int main(int argc, char *argv[]) {
  using Teuchos::RCP;
  using std::cout;
  using std::endl;

  RCP<Teuchos::FancyOStream> fos = getFancyOStream(Teuchos::rcpFromRef(std::cout));


  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  *fos << "Hello world" << std::endl;

  size_t nullSpaceDimension = 3;
  size_t myAggSize = 6;
  Teuchos::ArrayRCP<double> localQR(myAggSize*nullSpaceDimension); // reserve memory for local null space

  // fill local fine level nullspace
  for (size_t j=0; j<nullSpaceDimension; ++j) { // loop over columns
    for (int k=0; k<myAggSize; ++k) {            // loop over rows
      if(k % nullSpaceDimension == j)
        localQR[j* myAggSize + k] = 1.0;
      else
        localQR[j* myAggSize + k] = 0.0;
    }
  }
  localQR[1] = 0.5;
  localQR[2] = 0.2;
  localQR[6] = 0.0;
  localQR[7] = 0.0;
  localQR[10] = 0.0;

  // print input matrix
  for (int k=0; k<myAggSize; ++k) {            // loop over rows
    for (size_t j=0; j<nullSpaceDimension; ++j) { // loop over columns
      *fos << localQR[j*myAggSize + k] << "\t";
    }
    *fos << std::endl;
  }

  MueLu::QR_InterfaceEx<double,int> qrtest(nullSpaceDimension);

  qrtest.Compute(myAggSize, localQR);

  // define output array for coarse null space
  Teuchos::ArrayRCP< Teuchos::ArrayRCP<double> > coarseNS(nullSpaceDimension); // array of arrays!
  for (size_t j=0; j<nullSpaceDimension; ++j) {
    coarseNS[j] = Teuchos::ArrayRCP<double>(nullSpaceDimension); // max size is nullSpaceDimension, since we're only interested in upper triangle
  }

  int offset = 0;
  for (size_t j=0; j<nullSpaceDimension; ++j) {
    for (size_t k=0; k<=j; ++k) {
      // fill j-th null space vector
      coarseNS[j][offset+k] = localQR[ myAggSize*j + k ]; //TODO is offset+k the correct local ID?!
    }
  }

  *fos << std::endl;
  *fos << "coarse level nullspace part R=" << std::endl;

  // print coarse level null space
  for (size_t k=0; k<nullSpaceDimension; ++k) {
    for (size_t j=0; j<nullSpaceDimension; ++j) {
      *fos << coarseNS[j][k] << "\t";
    }
    *fos << std::endl;
  }

  // Calculate Q, the tentative prolongator.
  // The Lapack GEQRF call only works for myAggsize >= NSDim
  qrtest.ExtractQ(myAggSize, localQR);

  // define array for Ptent output
  Teuchos::ArrayRCP< Teuchos::ArrayRCP<double> > Ptent(nullSpaceDimension); // array of arrays!
  for (size_t j=0; j<nullSpaceDimension; ++j) {
    Ptent[j] = Teuchos::ArrayRCP<double>(myAggSize);
  }

  for (size_t r=0; r<myAggSize; ++r) {

      for (size_t j=0; j<nullSpaceDimension; ++j) {
        // fill j-th column of Ptent output array
        Ptent[j][r] = localQR[j*myAggSize+r];
      }
  }

  // print prolongator part
  *fos << std::endl;
  *fos << " Q= " << std::endl;
  for (size_t k=0; k<myAggSize; ++k) {
    for (size_t j=0; j<nullSpaceDimension; ++j) {
      *fos << Ptent[j][k] << "\t";
    }
    *fos << std::endl;
  }

  return EXIT_SUCCESS;

}
