// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_01.cpp
    \brief Test SROMvector interface.
*/


#include "ROL_SROMVector.hpp"
#include "ROL_TeuchosBatchManager.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_DefaultComm.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv,0);
  ROL::Ptr<const Teuchos::Comm<int> > comm
    = ROL::toPtr(Teuchos::DefaultComm<int>::getComm());

  int iprint = argc - 1;
  ROL::nullstream bhs; // outputs nothing
  std::ostream& outStream = (iprint > 0 && !Teuchos::rank<int>(*comm)) ? std::cout : bhs;

  int errorFlag = 0;

  RealT errtol = ROL::ROL_THRESHOLD<RealT>();

  try {
    // Batch manager
    ROL::Ptr<ROL::BatchManager<RealT> > bman =
      ROL::makePtr<ROL::TeuchosBatchManager<RealT,int>>(comm);

    // Dimension of the optimization vector
    int dimension = 5, numMyAtoms = 10;
    int size = dimension*numMyAtoms;

    // Create batch std vectors 
    ROL::Ptr<std::vector<RealT> > b1_ptr = ROL::makePtr<std::vector<RealT>>(size);
    ROL::Ptr<std::vector<RealT> > b2_ptr = ROL::makePtr<std::vector<RealT>>(size);
    ROL::Ptr<std::vector<RealT> > b3_ptr = ROL::makePtr<std::vector<RealT>>(size);
    for (int i = 0; i < size; ++i) {
      (*b1_ptr)[i] = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
      (*b2_ptr)[i] = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
      (*b3_ptr)[i] = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
    }
    ROL::Ptr<ROL::BatchStdVector<RealT> > b1
      = ROL::makePtr<ROL::BatchStdVector<RealT>>(b1_ptr,bman);
    ROL::Ptr<ROL::BatchStdVector<RealT> > b2
      = ROL::makePtr<ROL::BatchStdVector<RealT>>(b2_ptr,bman);
    ROL::Ptr<ROL::BatchStdVector<RealT> > b3
      = ROL::makePtr<ROL::BatchStdVector<RealT>>(b3_ptr,bman);

    // Create atom vectors 
    ROL::Ptr<std::vector<RealT> > a1_ptr = ROL::makePtr<std::vector<RealT>>(size);
    ROL::Ptr<std::vector<RealT> > a2_ptr = ROL::makePtr<std::vector<RealT>>(size);
    ROL::Ptr<std::vector<RealT> > a3_ptr = ROL::makePtr<std::vector<RealT>>(size);
    ROL::Ptr<std::vector<RealT> > aW_ptr = ROL::makePtr<std::vector<RealT>>(size);
    for (int i = 0; i < size; ++i) {
      (*a1_ptr)[i] = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
      (*a2_ptr)[i] = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
      (*a3_ptr)[i] = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
      (*aW_ptr)[i] = static_cast<RealT>(2);
    }
    ROL::Ptr<ROL::PrimalAtomVector<RealT> > a1
      = ROL::makePtr<ROL::PrimalAtomVector<RealT>>(a1_ptr,bman,numMyAtoms,dimension,aW_ptr);
    ROL::Ptr<ROL::PrimalAtomVector<RealT> > a2
      = ROL::makePtr<ROL::PrimalAtomVector<RealT>>(a2_ptr,bman,numMyAtoms,dimension,aW_ptr);
    ROL::Ptr<ROL::PrimalAtomVector<RealT> > a3
      = ROL::makePtr<ROL::PrimalAtomVector<RealT>>(a3_ptr,bman,numMyAtoms,dimension,aW_ptr);

    // Create probability vectors
    ROL::Ptr<std::vector<RealT> > p1_ptr = ROL::makePtr<std::vector<RealT>>(numMyAtoms);
    ROL::Ptr<std::vector<RealT> > p2_ptr = ROL::makePtr<std::vector<RealT>>(numMyAtoms);
    ROL::Ptr<std::vector<RealT> > p3_ptr = ROL::makePtr<std::vector<RealT>>(numMyAtoms);
    ROL::Ptr<std::vector<RealT> > pW_ptr = ROL::makePtr<std::vector<RealT>>(numMyAtoms);
    for (int i = 0; i < numMyAtoms; ++i) {
      (*p1_ptr)[i] = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
      (*p2_ptr)[i] = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
      (*p3_ptr)[i] = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
      (*pW_ptr)[i] = static_cast<RealT>(2);
    }
    ROL::Ptr<ROL::PrimalProbabilityVector<RealT> > p1
      = ROL::makePtr<ROL::PrimalProbabilityVector<RealT>>(p1_ptr,bman,pW_ptr);
    ROL::Ptr<ROL::PrimalProbabilityVector<RealT> > p2
      = ROL::makePtr<ROL::PrimalProbabilityVector<RealT>>(p2_ptr,bman,pW_ptr);
    ROL::Ptr<ROL::PrimalProbabilityVector<RealT> > p3
      = ROL::makePtr<ROL::PrimalProbabilityVector<RealT>>(p3_ptr,bman,pW_ptr);

    // Create SROM vectors
    ROL::SROMVector<RealT> x1(p1,a1);
    ROL::SROMVector<RealT> x2(p2,a2);
    ROL::SROMVector<RealT> x3(p3,a3);

    // Standard tests.
    std::vector<RealT> consistencyBMAN = b1->checkVector(*b2, *b3, true, outStream);
    ROL::StdVector<RealT> checkvecBMAN( ROL::makePtrFromRef(consistencyBMAN) );
    if (checkvecBMAN.norm() > std::sqrt(errtol)) {
      errorFlag++;
    }
    std::vector<RealT> consistencyAtom = a1->checkVector(*a2, *a3, true, outStream);
    ROL::StdVector<RealT> checkvecAtom( ROL::makePtrFromRef(consistencyAtom) );
    if (checkvecAtom.norm() > std::sqrt(errtol)) {
      errorFlag++;
    }
    std::vector<RealT> consistencyProb = p1->checkVector(*p2, *p3, true, outStream);
    ROL::StdVector<RealT> checkvecProb( ROL::makePtrFromRef(consistencyProb) );
    if (checkvecProb.norm() > std::sqrt(errtol)) {
      errorFlag++;
    }
    std::vector<RealT> consistencySROM = x1.checkVector(x2, x3, true, outStream);
    ROL::StdVector<RealT> checkvecSROM( ROL::makePtrFromRef(consistencySROM) );
    if (checkvecSROM.norm() > std::sqrt(errtol)) {
      errorFlag++;
    }

    RealT numProcs = static_cast<RealT>(Teuchos::size<int>(*comm));
    RealT anorm = std::sqrt(numProcs*size), pnorm = std::sqrt(numProcs*numMyAtoms);
    RealT norm = std::sqrt(anorm*anorm + pnorm*pnorm);
    RealT sqrt2 = static_cast<RealT>(std::sqrt(2.));

    // Create batch std vectors 
    ROL::Ptr<std::vector<RealT> > b_ptr = ROL::makePtr<std::vector<RealT>>(size,1);
    ROL::Ptr<ROL::BatchStdVector<RealT> > b
      = ROL::makePtr<ROL::BatchStdVector<RealT>>(b_ptr,bman);
    RealT bnorm = b->norm();
    outStream << "BatchStdVector Norm Error:          "
              << std::abs(bnorm - anorm) << std::endl;
    if ( std::abs(bnorm - anorm) > std::sqrt(errtol) ) {
      errorFlag++;
    }

    // Create atom vectors 
    ROL::Ptr<std::vector<RealT> > ap_ptr = ROL::makePtr<std::vector<RealT>>(size,1);
    ROL::Ptr<ROL::PrimalAtomVector<RealT> > ap
      = ROL::makePtr<ROL::PrimalAtomVector<RealT>>(ap_ptr,bman,numMyAtoms,dimension,aW_ptr);
    RealT apnorm = ap->norm();
    outStream << "PrimalAtomVector Norm Error:        "
              << std::abs(apnorm - sqrt2*anorm) << std::endl;
    if ( std::abs(apnorm - sqrt2*anorm) > std::sqrt(errtol) ) {
      errorFlag++;
    }
    ROL::Ptr<std::vector<RealT> > ad_ptr = ROL::makePtr<std::vector<RealT>>(size,1);
    ROL::Ptr<ROL::DualAtomVector<RealT> > ad
      = ROL::makePtr<ROL::DualAtomVector<RealT>>(ad_ptr,bman,numMyAtoms,dimension,aW_ptr);
    RealT adnorm = ad->norm();
    outStream << "DualAtomVector Norm Error:          "
              << std::abs(adnorm - anorm/sqrt2) << std::endl;
    if ( std::abs(adnorm - anorm/sqrt2) > std::sqrt(errtol) ) {
      errorFlag++;
    }

    // Create probability vectors
    ROL::Ptr<std::vector<RealT> > pp_ptr = ROL::makePtr<std::vector<RealT>>(numMyAtoms,1);
    ROL::Ptr<ROL::PrimalProbabilityVector<RealT> > pp
      = ROL::makePtr<ROL::PrimalProbabilityVector<RealT>>(pp_ptr,bman,pW_ptr);
    RealT ppnorm = pp->norm();
    outStream << "PrimalProbabilityVector Norm Error: "
              << std::abs(ppnorm - sqrt2*pnorm) << std::endl;
    if ( std::abs(ppnorm - sqrt2*pnorm) > std::sqrt(errtol) ) {
      errorFlag++;
    }
    ROL::Ptr<std::vector<RealT> > pd_ptr = ROL::makePtr<std::vector<RealT>>(numMyAtoms,1);
    ROL::Ptr<ROL::DualProbabilityVector<RealT> > pd
      = ROL::makePtr<ROL::DualProbabilityVector<RealT>>(pd_ptr,bman,pW_ptr);
    RealT pdnorm = pd->norm();
    outStream << "DualProbabilityVector Norm Error:   "
              << std::abs(pdnorm - pnorm/sqrt2) << std::endl;
    if ( std::abs(pdnorm - pnorm/sqrt2) > std::sqrt(errtol) ) {
      errorFlag++;
    }

    // Create SROM vectors
    ROL::SROMVector<RealT> xp(pp,ap);
    RealT xpnorm = xp.norm();
    outStream << "PrimalSROMVector Norm Error:        "
              << std::abs(xpnorm - sqrt2*norm) << std::endl;
    if ( std::abs(xpnorm - sqrt2*norm) > std::sqrt(errtol) ) {
      errorFlag++;
    }
    ROL::SROMVector<RealT> xd(pd,ad);
    RealT xdnorm = xd.norm();
    outStream << "DualSROMVector Norm Error:          "
              << std::abs(xdnorm - norm/sqrt2) << std::endl;
    if ( std::abs(xdnorm - norm/sqrt2) > std::sqrt(errtol) ) {
      errorFlag++;
    }
    outStream << std::endl;
  }

  catch (std::logic_error& err) {
    outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}
