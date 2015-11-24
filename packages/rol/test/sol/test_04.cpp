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

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_BatchManager.hpp"
#include "ROL_SROMGenerator.hpp"
#include "ROL_MomentObjective.hpp"
#include "ROL_CDFObjective.hpp"
#include "ROL_LinearCombinationObjective.hpp"
//#include "ROL_SROMBoundConstraint.hpp"
#include "ROL_SROMVector.hpp"
#include "ROL_DistributionFactory.hpp"

int main(int argc, char* argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

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
    Teuchos::RCP<ROL::Distribution<double> > dist;
    std::vector<Teuchos::RCP<ROL::Distribution<double> > > distVec(dimension);
    Teuchos::ParameterList Dlist;
    Dlist.sublist("Distribution").set("Name","Beta");
    double alpha = 1., beta = 4.;
    // Fill moment vector and initial guess
    for (size_t d = 0; d < dimension; d++) {
      // Build distribution for dimension d
      alpha++; beta++;
      Dlist.sublist("Distribution").sublist("Beta").set("Shape 1",alpha);
      Dlist.sublist("Distribution").sublist("Beta").set("Shape 2",beta);
      dist = ROL::DistributionFactory<double>(Dlist);
      distVec[d] = ROL::DistributionFactory<double>(Dlist);
    }

    // Get ROL parameterlist
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    Teuchos::ParameterList &list = parlist->sublist("SOL").sublist("Sample Generator").sublist("SROM");
    Teuchos::Array<int> moments = Teuchos::getArrayFromStringParameter<int>(list,"Moments");
    size_t numMoments = static_cast<size_t>(moments.size());

    Teuchos::RCP<ROL::BatchManager<double> > bman =
      Teuchos::rcp(new ROL::BatchManager<double>());
    Teuchos::RCP<ROL::SampleGenerator<double> > sampler =
      Teuchos::rcp(new ROL::SROMGenerator<double>(*parlist,bman,distVec));

    double val = 0., error = 0., data = 0.;
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
        error = std::abs(val-data)/std::abs(data);
        if ( error > 1.e-2 ) {
          errorFlag++;
        }
        *outStream << std::right << std::setw(20) << val
                                 << std::setw(20) << data
                                 << std::setw(20) << error
                                 << std::endl;
      }
    }
    *outStream << std::endl;

    for (size_t k = 0; k < (size_t)sampler->numMySamples(); k++) {
      for (size_t d = 0; d < dimension; d++) {
        *outStream << std::setw(20) << (sampler->getMyPoint(k))[d] << "  ";
      }
      *outStream << std::setw(20) << sampler->getMyWeight(k) << std::endl;
    }
    *outStream << std::endl;
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
