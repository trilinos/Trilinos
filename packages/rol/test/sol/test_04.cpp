// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#include "ROL_ParameterList.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"

#include "ROL_TeuchosBatchManager.hpp"
#include "ROL_SROMGenerator.hpp"
#include "ROL_DistributionFactory.hpp"

#include <ctime>

typedef double RealT;

int main(int argc, char* argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  ROL::Ptr<const Teuchos::Comm<int> > commptr =
    ROL::toPtr(Teuchos::DefaultComm<int>::getComm());

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0 && commptr->getRank() == 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  try {
    /**********************************************************************************************/
    /************************* CONSTRUCT SOL COMPONENTS *******************************************/
    /**********************************************************************************************/
    // Set random seed
    srand(123456789);
    // Build samplers
    size_t dimension = 1;

    // Initialize distribution
    ROL::Ptr<ROL::Distribution<RealT> > dist;
    std::vector<ROL::Ptr<ROL::Distribution<RealT> > > distVec(dimension);
    ROL::ParameterList Dlist;
    Dlist.sublist("SOL").sublist("Distribution").set("Name","Beta");
    RealT alpha = 1., beta = 4.;
    // Fill moment vector and initial guess
    for (size_t d = 0; d < dimension; d++) {
      // Build distribution for dimension d
      alpha++; beta++;
      Dlist.sublist("SOL").sublist("Distribution").sublist("Beta").set("Shape 1",alpha);
      Dlist.sublist("SOL").sublist("Distribution").sublist("Beta").set("Shape 2",beta);
      dist = ROL::DistributionFactory<RealT>(Dlist);
      distVec[d] = ROL::DistributionFactory<RealT>(Dlist);
    }

    // Get ROL parameterlist
    std::string filename = "input_04.xml";
    
    auto parlist = ROL::getParametersFromXmlFile( filename );

    ROL::ParameterList &list = parlist->sublist("SOL").sublist("Sample Generator").sublist("SROM");
    std::vector<int> moments = ROL::getArrayFromStringParameter<int>(list,"Moments");
    size_t numMoments = static_cast<size_t>(moments.size());

    std::clock_t timer = std::clock();
    ROL::Ptr<ROL::BatchManager<RealT> > bman =
      ROL::makePtr<ROL::TeuchosBatchManager<RealT,int>>(commptr);
    ROL::Ptr<ROL::SampleGenerator<RealT> > sampler =
      ROL::makePtr<ROL::SROMGenerator<RealT>>(*parlist,bman,distVec,*outStream);
    *outStream << std::endl << "Sample Time: "
               << (std::clock()-timer)/(RealT)CLOCKS_PER_SEC << " seconds"
               << std::endl;

    RealT val = 0., error = 0., data = 0., sum = 0.;
    *outStream << std::endl;
    *outStream << std::scientific << std::setprecision(11);
    *outStream << std::right << std::setw(20) << "Computed Moment"
                             << std::setw(20) << "True Moment"
                             << std::setw(20) << "Relative Error"
                             << std::endl;
    for (size_t m = 0; m < numMoments; m++) {
      for (size_t d = 0; d < dimension; d++) {
        val = 0.; data = distVec[d]->moment(moments[m]);
        for (size_t k = 0; k < (size_t)sampler->numMySamples(); k++) {
          val += sampler->getMyWeight(k)*std::pow((sampler->getMyPoint(k))[d],moments[m]);
        }
        bman->sumAll(&val,&sum,1);
        error = std::abs(sum-data)/std::abs(data);
        if ( error > 1.e-1 ) {
          errorFlag++;
        }
        *outStream << std::right << std::setw(20) << sum
                                 << std::setw(20) << data
                                 << std::setw(20) << error
                                 << std::endl;
      }
    }
    *outStream << std::endl;

//    std::ofstream file;
//    std::stringstream name;
//    name << "samples." << commptr->getRank() << ".txt";
//    file.open(name.str().c_str());
//    for (size_t k = 0; k < (size_t)sampler->numMySamples(); k++) {
//      for (size_t d = 0; d < dimension; d++) {
//        file << std::setprecision(std::numeric_limits<RealT>::digits10)
//             << std::scientific
//             << (sampler->getMyPoint(k))[d];
//        file << "  ";
//      }
//      file << std::setprecision(std::numeric_limits<RealT>::digits10)
//           << std::scientific
//           << sampler->getMyWeight(k) << std::endl;
//    }
//    file.close();
//    commptr->barrier();

  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}
