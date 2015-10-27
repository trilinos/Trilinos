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
#include "ROL_SROMBoundConstraint.hpp"
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
    size_t nSamp = 5;
    size_t numMoments = 5;

    // Initialize sample and weight vector
    Teuchos::RCP<std::vector<double> > xpt
      = Teuchos::rcp(new std::vector<double>(dimension*nSamp,0.));
    Teuchos::RCP<std::vector<double> > xwt
      = Teuchos::rcp(new std::vector<double>(nSamp,1./(double)nSamp));
    ROL::SROMVector<double> x(xpt,xwt);
    Teuchos::RCP<ROL::Vector<double> > xptr = Teuchos::rcp(&x,false);
    // Initialize distribution
    Teuchos::RCP<ROL::Distribution<double> > dist;
    std::vector<Teuchos::RCP<ROL::Distribution<double> > > distVec(dimension);
    Teuchos::ParameterList Dlist;
    Dlist.sublist("SOL").sublist("Distribution").set("Name","Beta");
    double alpha = 1., beta = 4.;
    std::vector<std::vector<std::pair<size_t,double> > > moments(dimension);
    std::vector<std::pair<size_t,double> > data(numMoments);
    // Fill moment vector and initial guess
    for (size_t d = 0; d < dimension; d++) {
      // Build distribution for dimension d
      alpha++; beta++;
      Dlist.sublist("SOL").sublist("Distribution").sublist("Beta").set("Shape 1",alpha);
      Dlist.sublist("SOL").sublist("Distribution").sublist("Beta").set("Shape 2",beta);
      dist = ROL::DistributionFactory<double>(Dlist);
      distVec[d] = ROL::DistributionFactory<double>(Dlist);
      // Compute moments
      for (size_t m = 0; m < numMoments; m++) {
        data[m] = std::make_pair(m+1,dist->moment(m+1));
      }
      moments[d].assign(data.begin(),data.end());
      // Set initial sample guess to random samples
      for (size_t k = 0; k < nSamp; k++) {
        (*xpt)[k*dimension + d] = dist->invertCDF((double)rand()/(double)RAND_MAX);
      }
    }
    // Initialize bound constraints
    std::vector<double> xlo(dimension,0.), xup(dimension,1.);
    Teuchos::RCP<ROL::BoundConstraint<double> > bnd
      = Teuchos::rcp(new ROL::SROMBoundConstraint<double>(xlo,xup));
    // Initialize objective function
    Teuchos::RCP<ROL::Objective<double> > obj_moment
      = Teuchos::rcp(new ROL::MomentObjective<double>(moments));
    double scale = 1.e-1;
    Teuchos::RCP<ROL::Objective<double> > obj_CDF
      = Teuchos::rcp(new ROL::CDFObjective<double>(distVec,xlo,xup,scale));
    std::vector<double> weights(2,0.);
    weights[0] = 1.; weights[1] = 1.;
    std::vector<Teuchos::RCP<ROL::Objective<double> > > objVec(2);
    objVec[0] = obj_moment; objVec[1] = obj_CDF;
    Teuchos::RCP<ROL::Objective<double> > obj
      = Teuchos::rcp(new ROL::LinearCombinationObjective<double>(weights,objVec));

    bool derivCheck = true;
    if ( derivCheck ) {
      Teuchos::RCP<std::vector<double> > ypt
        = Teuchos::rcp(new std::vector<double>(dimension*nSamp,0.));
      Teuchos::RCP<std::vector<double> > ywt
        = Teuchos::rcp(new std::vector<double>(nSamp,0.));
      ROL::SROMVector<double> y(ypt,ywt);
      for (size_t k = 0; k < nSamp; k++) {
        for (size_t d = 0; d < dimension; d++) {
          (*ypt)[k*dimension + d] = 2.*(double)rand()/(double)RAND_MAX - 1.;
        }
        (*ywt)[k] = 2.*(double)rand()/(double)RAND_MAX - 1.;
      }
      *outStream << "\n  CHECK MOMENT OBJECTIVE DERIVATIVES" << std::endl;
      obj_moment->checkGradient(x,y,true);
      obj_moment->checkHessVec(x,y,true);
      *outStream << "\n  CHECK CDF OBJECTIVE DERIVATIVES" << std::endl;
      obj_CDF->checkGradient(x,y,true);
      obj_CDF->checkHessVec(x,y,true);
      *outStream << "\n  CHECK COMBINED OBJECTIVE DERIVATIVES" << std::endl;
      obj->checkGradient(x,y,true);
      obj->checkHessVec(x,y,true);
    }

    // Get ROL parameterlist
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    Teuchos::RCP<ROL::BatchManager<double> > bman =
      Teuchos::rcp(new ROL::BatchManager<double>());
    Teuchos::RCP<ROL::SampleGenerator<double> > sampler =
      Teuchos::rcp(new ROL::SROMGenerator<double>(*parlist,bman,obj,bnd,xptr,dimension,nSamp));

    double val = 0., error = 0.;
    *outStream << std::endl;
    *outStream << std::scientific << std::setprecision(11);
    *outStream << std::right << std::setw(20) << "Computed Moment"
                             << std::setw(20) << "True Moment"
                             << std::setw(20) << "Relative Error"
                             << std::endl;
    for (size_t m = 0; m < numMoments; m++) {
      for (size_t d = 0; d < dimension; d++) {
        val = 0.; data = moments[d];
        for (size_t k = 0; k < (size_t)sampler->numMySamples(); k++) {
          val += sampler->getMyWeight(k)*std::pow((sampler->getMyPoint(k))[d],data[m].first);
        }
        error = std::abs(val-data[m].second)/std::abs(data[m].second);
        if ( error > 1.e-2 ) {
          errorFlag++;
        }
        *outStream << std::right << std::setw(20) << val
                                 << std::setw(20) << data[m].second
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
