// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file MultiJagged.cpp
    \brief An example of partitioning coordinates with MultiJagged.
    \todo add more cases to this test.
 */

#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>
#include <Zoltan2_BasicVectorAdapter.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <GeometricGenerator.hpp>

#include <Zoltan2_EvaluatePartition.hpp>

#include "Teuchos_XMLParameterListHelpers.hpp"

#include <Teuchos_LAPACK.hpp>
#include <fstream>
#include <string>
using Teuchos::RCP;
using Teuchos::rcp;


//#define hopper_separate_test
#ifdef hopper_separate_test
#include "stdio.h"
#endif
#define CATCH_EXCEPTIONS_AND_RETURN(pp) \
        catch (std::runtime_error &e) { \
            std::cout << "Runtime exception returned from " << pp << ": " \
            << e.what() << " FAIL" << std::endl; \
            return -1; \
        } \
        catch (std::logic_error &e) { \
            std::cout << "Logic exception returned from " << pp << ": " \
            << e.what() << " FAIL" << std::endl; \
            return -1; \
        } \
        catch (std::bad_alloc &e) { \
            std::cout << "Bad_alloc exception returned from " << pp << ": " \
            << e.what() << " FAIL" << std::endl; \
            return -1; \
        } \
        catch (std::exception &e) { \
            std::cout << "Unknown exception returned from " << pp << ": " \
            << e.what() << " FAIL" << std::endl; \
            return -1; \
        }

#define CATCH_EXCEPTIONS_WITH_COUNT(ierr, pp) \
        catch (std::runtime_error &e) { \
            std::cout << "Runtime exception returned from " << pp << ": " \
            << e.what() << " FAIL" << std::endl; \
            (ierr)++; \
        } \
        catch (std::logic_error &e) { \
            std::cout << "Logic exception returned from " << pp << ": " \
            << e.what() << " FAIL" << std::endl; \
            (ierr)++; \
        } \
        catch (std::bad_alloc &e) { \
            std::cout << "Bad_alloc exception returned from " << pp << ": " \
            << e.what() << " FAIL" << std::endl; \
            (ierr)++; \
        } \
        catch (std::exception &e) { \
            std::cout << "Unknown exception returned from " << pp << ": " \
            << e.what() << " FAIL" << std::endl; \
            (ierr)++; \
        }

/*! \test MultiJaggedTest.cpp
    An example of the use of the MultiJagged algorithm to partition coordinate data.
 */


const char param_comment = '#';

string trim_right_copy(
        const string& s,
        const string& delimiters = " \f\n\r\t\v" )
{
    return s.substr( 0, s.find_last_not_of( delimiters ) + 1 );
}

string trim_left_copy(
        const string& s,
        const string& delimiters = " \f\n\r\t\v" )
{
    return s.substr( s.find_first_not_of( delimiters ) );
}

string trim_copy(
        const string& s,
        const string& delimiters = " \f\n\r\t\v" )
{
    return trim_left_copy( trim_right_copy( s, delimiters ), delimiters );
}

template <typename Adapter>
void print_boxAssign_result(
  const char *str,
  int dim,
  typename Adapter::scalar_t *lower,
  typename Adapter::scalar_t *upper,
  size_t nparts,
  typename Adapter::part_t *parts
)
{
  std::cout << "boxAssign test " << str << ":  Box (";
  for (int j = 0; j < dim; j++) std::cout << lower[j] << " ";
  std::cout << ") x (";
  for (int j = 0; j < dim; j++) std::cout << upper[j] << " ";

  if (nparts == 0)
    std::cout << ") does not overlap any parts" << std::endl;
  else {
    std::cout << ") overlaps parts ";
    for (size_t k = 0; k < nparts; k++) std::cout << parts[k] << " ";
    std::cout << std::endl;
  }
}

template <typename Adapter>
int run_pointAssign_tests(
  Zoltan2::PartitioningProblem<Adapter> *problem,
  RCP<tMVector_t> &coords,
  bool print_details)
{
    int ierr = 0;

    // pointAssign tests
    int coordDim = coords->getNumVectors();
    zscalar_t *pointDrop = new zscalar_t[coordDim];
    typename Adapter::part_t part = -1;

    char mechar[10];
    sprintf(mechar, "%d", problem->getComm()->getRank());
    string me(mechar);

    // test correctness of pointAssign for owned points
    {
      const typename Adapter::part_t *solnPartView =
                                      problem->getSolution().getPartListView();

      size_t numPoints = coords->getLocalLength();
      for (size_t localID = 0; localID < numPoints; localID++) {

        typename Adapter::part_t solnPart = solnPartView[localID];

        for (int i = 0; i < coordDim; i++)
          pointDrop[i] = coords->getData(i)[localID];

        try {
          part = problem->getSolution().pointAssign(coordDim, pointDrop);
        }
        CATCH_EXCEPTIONS_WITH_COUNT(ierr, me + ": pointAssign -- OwnedPoints");

        if(print_details) {
          std::cout << me << " Point " << localID
                    << " gid " << coords->getMap()->getGlobalElement(localID)
                    << " (" << pointDrop[0];
          if (coordDim > 1) std::cout << " " << pointDrop[1];
          if (coordDim > 2) std::cout << " " << pointDrop[2];
          std::cout << ") in boxPart " << part
                    << "  in solnPart " << solnPart
                    << std::endl;
        }

// this error test does not work for points that fall on the cuts.
// like Zoltan's RCB, pointAssign arbitrarily picks a part along the cut.
// the arbitrarily chosen part will not necessarily be the one to which
// the coordinate was assigned in partitioning.
//
//        if (part != solnPart) {
//          std::cout << me << " pointAssign:  incorrect part " << part
//                    << " found; should be " << solnPart
//                    << " for point " << j << std::endl;
//          ierr++;
//        }
      }
    }

    {
      const std::vector<Zoltan2::coordinateModelPartBox>
            pBoxes = problem->getSolution().getPartBoxesView();
      if(print_details) {
        for (size_t i = 0; i < pBoxes.size(); i++) {
          typename Zoltan2::coordinateModelPartBox::coord_t *lmin = pBoxes[i].getlmins();
          typename Zoltan2::coordinateModelPartBox::coord_t *lmax = pBoxes[i].getlmaxs();;
          std::cout << me << " pBox " << i << " pid " << pBoxes[i].getpId()
                    << " (" << lmin[0] << "," << lmin[1] << ","
                    << (coordDim > 2 ? lmin[2] : 0) << ") x "
                    << " (" << lmax[0] << "," << lmax[1] << ","
                    << (coordDim > 2 ? lmax[2] : 0) << ")" << std::endl;
        }
      }
    }

    // test the origin
    {
      for (int i = 0; i < coordDim; i++) pointDrop[i] = 0.;
      try {
        part = problem->getSolution().pointAssign(coordDim, pointDrop);
      }
      CATCH_EXCEPTIONS_WITH_COUNT(ierr, me + " pointAssign -- Origin");

      std::cout << me << " OriginPoint (" << pointDrop[0];
      if (coordDim > 1) std::cout << " " << pointDrop[1];
      if (coordDim > 2) std::cout << " " << pointDrop[2];
      std::cout << ")  part " << part << std::endl;
    }

    // test point with negative coordinates
    {
      for (int i = 0; i < coordDim; i++) pointDrop[i] = -100.+i;
      try {
        part = problem->getSolution().pointAssign(coordDim, pointDrop);
      }
      CATCH_EXCEPTIONS_WITH_COUNT(ierr, me + " pointAssign -- Negative Point");
      std::cout << me << " NegativePoint (" << pointDrop[0];
      if (coordDim > 1) std::cout << " " << pointDrop[1];
      if (coordDim > 2) std::cout << " " << pointDrop[2];
      std::cout << ")  part " << part << std::endl;
    }

    // test a point that's way out there
    {
      for (int i = 0; i < coordDim; i++) pointDrop[i] = i*5;
      try {
        part = problem->getSolution().pointAssign(coordDim, pointDrop);
      }
      CATCH_EXCEPTIONS_WITH_COUNT(ierr, me + " pointAssign -- i*5 Point");
      std::cout << me << " i*5-Point (" << pointDrop[0];
      if (coordDim > 1) std::cout << " " << pointDrop[1];
      if (coordDim > 2) std::cout << " " << pointDrop[2];
      std::cout << ")  part " << part << std::endl;
    }

    // test a point that's way out there
    {
      for (int i = 0; i < coordDim; i++) pointDrop[i] = 10+i*5;
      try {
        part = problem->getSolution().pointAssign(coordDim, pointDrop);
      }
      CATCH_EXCEPTIONS_WITH_COUNT(ierr, me + " pointAssign -- WoopWoop");
      std::cout << me << " WoopWoop-Point (" << pointDrop[0];
      if (coordDim > 1) std::cout << " " << pointDrop[1];
      if (coordDim > 2) std::cout << " " << pointDrop[2];
      std::cout << ")  part " << part << std::endl;
    }

    delete [] pointDrop;
    return ierr;
}

template <typename Adapter>
int run_boxAssign_tests(
  Zoltan2::PartitioningProblem<Adapter> *problem,
  RCP<tMVector_t> &coords)
{
    int ierr = 0;

    // boxAssign tests
    int coordDim = coords->getNumVectors();
    zscalar_t *lower = new zscalar_t[coordDim];
    zscalar_t *upper = new zscalar_t[coordDim];

    char mechar[10];
    sprintf(mechar, "%d", problem->getComm()->getRank());
    string me(mechar);

    const std::vector<Zoltan2::coordinateModelPartBox>
          pBoxes = problem->getSolution().getPartBoxesView();
    size_t nBoxes = pBoxes.size();

    // test a box that is smaller than a part
    {
      size_t nparts;
      typename Adapter::part_t *parts;
      size_t pickabox = nBoxes / 2;
      for (int i = 0; i < coordDim; i++) {
        zscalar_t dd = 0.2 * (pBoxes[pickabox].getlmaxs()[i] -
                              pBoxes[pickabox].getlmins()[i]);
        lower[i] = pBoxes[pickabox].getlmins()[i] + dd;
        upper[i] = pBoxes[pickabox].getlmaxs()[i] - dd;
      }
      try {
        problem->getSolution().boxAssign(coordDim, lower, upper,
                                         nparts, &parts);
      }
      CATCH_EXCEPTIONS_WITH_COUNT(ierr, me + " boxAssign -- smaller");
      if (nparts > 1) {
        std::cout << me << " FAIL boxAssign error: smaller test, nparts > 1"
                  << std::endl;
        ierr++;
      }
      print_boxAssign_result<Adapter>("smallerbox", coordDim,
                                      lower, upper, nparts, parts);
      delete [] parts;
    }

    // test a box that is larger than a part
    {
      size_t nparts;
      typename Adapter::part_t *parts;
      size_t pickabox = nBoxes / 2;
      for (int i = 0; i < coordDim; i++) {
        zscalar_t dd = 0.2 * (pBoxes[pickabox].getlmaxs()[i] -
                              pBoxes[pickabox].getlmins()[i]);
        lower[i] = pBoxes[pickabox].getlmins()[i] - dd;
        upper[i] = pBoxes[pickabox].getlmaxs()[i] + dd;
      }
      try {
        problem->getSolution().boxAssign(coordDim, lower, upper,
                                         nparts, &parts);
      }
      CATCH_EXCEPTIONS_WITH_COUNT(ierr, me + " boxAssign -- larger");

      // larger box should have at least two parts in it for k > 1.
      if ((nBoxes > 1) && (nparts < 2)) {
        std::cout << me << " FAIL boxAssign error: "
                  << "larger test, nparts < 2"
                  << std::endl;
        ierr++;
      }

      // parts should include pickabox's part
      bool found_pickabox = 0;
      for (size_t i = 0; i < nparts; i++)
        if (parts[i] == pBoxes[pickabox].getpId()) {
           found_pickabox = 1;
           break;
        }
      if (!found_pickabox) {
        std::cout << me << " FAIL boxAssign error: "
                  << "larger test, pickabox not found"
                  << std::endl;
        ierr++;
      }

      print_boxAssign_result<Adapter>("largerbox", coordDim,
                                      lower, upper, nparts, parts);
      delete [] parts;
    }

    // test a box that includes all parts
    {
      size_t nparts;
      typename Adapter::part_t *parts;
      for (int i = 0; i < coordDim; i++) {
        lower[i] = std::numeric_limits<zscalar_t>::max();
        upper[i] = std::numeric_limits<zscalar_t>::min();
      }
      for (size_t j = 0; j < nBoxes; j++) {
        for (int i = 0; i < coordDim; i++) {
          if (pBoxes[j].getlmins()[i] <= lower[i])
            lower[i] = pBoxes[j].getlmins()[i];
          if (pBoxes[j].getlmaxs()[i] >= upper[i])
            upper[i] = pBoxes[j].getlmaxs()[i];
        }
      }
      try {
        problem->getSolution().boxAssign(coordDim, lower, upper,
                                         nparts, &parts);
      }
      CATCH_EXCEPTIONS_WITH_COUNT(ierr, me + " boxAssign -- global");

      // global box should have all parts
      if (nparts != nBoxes) {
        std::cout << me << " FAIL boxAssign error: "
                  << "global test, nparts found " << nparts
                  << " != num global parts " << nBoxes
                  << std::endl;
        ierr++;
      }
      print_boxAssign_result<Adapter>("globalbox", coordDim,
                                      lower, upper, nparts, parts);
      delete [] parts;
    }

    // test a box that is bigger than the entire domain
    // Assuming lower and upper are still set to the global box boundary
    // from the previous test
    {
      size_t nparts;
      typename Adapter::part_t *parts;
      for (int i = 0; i < coordDim; i++) {
        lower[i] -= 2.;
        upper[i] += 2.;
      }

      try {
        problem->getSolution().boxAssign(coordDim, lower, upper,
                                         nparts, &parts);
      }
      CATCH_EXCEPTIONS_WITH_COUNT(ierr, me + " boxAssign -- bigdomain");

      // bigdomain box should have all parts
      if (nparts != nBoxes) {
        std::cout << me << " FAIL boxAssign error: "
                  << "bigdomain test, nparts found " << nparts
                  << " != num global parts " << nBoxes
                  << std::endl;
        ierr++;
      }
      print_boxAssign_result<Adapter>("bigdomainbox", coordDim,
                                      lower, upper, nparts, parts);
      delete [] parts;
    }

    // test a box that is way out there
    // Assuming lower and upper are still set to at least the global box
    // boundary from the previous test
    {
      size_t nparts;
      typename Adapter::part_t *parts;
      for (int i = 0; i < coordDim; i++) {
        lower[i] = upper[i] + 10;
        upper[i] += 20;
      }

      try {
        problem->getSolution().boxAssign(coordDim, lower, upper,
                                         nparts, &parts);
      }
      CATCH_EXCEPTIONS_WITH_COUNT(ierr, me + " boxAssign -- out there");

      // For now, boxAssign returns zero if there is no box overlap.
      // TODO:  this result should be changed in boxAssign definition
      if (nparts != 0) {
        std::cout << me << " FAIL boxAssign error: "
                  << "outthere test, nparts found " << nparts
                  << " != zero"
                  << std::endl;
        ierr++;
      }
      print_boxAssign_result<Adapter>("outthere box", coordDim,
                                      lower, upper, nparts, parts);
      delete [] parts;
    }

    delete [] lower;
    delete [] upper;
    return ierr;
}

void readGeoGenParams(string paramFileName, Teuchos::ParameterList &geoparams, const RCP<const Teuchos::Comm<int> > & comm){
    std::string input = "";
    char inp[25000];
    for(int i = 0; i < 25000; ++i){
        inp[i] = 0;
    }

    bool fail = false;
    if(comm->getRank() == 0){

        std::fstream inParam(paramFileName.c_str());
        if (inParam.fail())
        {
            fail = true;
        }
        if(!fail)
        {
            std::string tmp = "";
            getline (inParam,tmp);
            while (!inParam.eof()){
                if(tmp != ""){
                    tmp = trim_copy(tmp);
                    if(tmp != ""){
                        input += tmp + "\n";
                    }
                }
                getline (inParam,tmp);
            }
            inParam.close();
            for (size_t i = 0; i < input.size(); ++i){
                inp[i] = input[i];
            }
        }
    }



    int size = input.size();
    if(fail){
        size = -1;
    }
    comm->broadcast(0, sizeof(int), (char*) &size);
    if(size == -1){
        throw "File " + paramFileName + " cannot be opened.";
    }
    comm->broadcast(0, size, inp);
    std::istringstream inParam(inp);
    string str;
    getline (inParam,str);
    while (!inParam.eof()){
        if(str[0] != param_comment){
            size_t pos = str.find('=');
            if(pos == string::npos){
                throw  "Invalid Line:" + str  + " in parameter file";
            }
            string paramname = trim_copy(str.substr(0,pos));
            string paramvalue = trim_copy(str.substr(pos + 1));
            geoparams.set(paramname, paramvalue);
        }
        getline (inParam,str);
    }
}

template<class bv_use_node_t>
int compareWithBasicVectorAdapterTest(RCP<const Teuchos::Comm<int> > &comm,
  Teuchos::RCP<Teuchos::ParameterList> params,
  Zoltan2::PartitioningProblem<Zoltan2::XpetraMultiVectorAdapter<tMVector_t>> *problem,
  RCP<tMVector_t> coords,
  Zoltan2::XpetraMultiVectorAdapter<tMVector_t>::scalar_t ** weights = NULL, int numWeightsPerCoord = 0) {

  typedef Zoltan2::XpetraMultiVectorAdapter<tMVector_t> inputAdapter_t;

  // Run a test with BasicVectorAdapter and xyzxyz format coordinates
  const int bvme = comm->getRank();
  const inputAdapter_t::lno_t bvlen =
                        inputAdapter_t::lno_t(coords->getLocalLength());
  const size_t bvnvecs = coords->getNumVectors();
  const size_t bvsize = coords->getNumVectors() * coords->getLocalLength();

  ArrayRCP<inputAdapter_t::scalar_t> *bvtpetravectors =
          new ArrayRCP<inputAdapter_t::scalar_t>[bvnvecs];
  for (size_t i = 0; i < bvnvecs; i++)
    bvtpetravectors[i] = coords->getDataNonConst(i);
  int idx = 0;
  inputAdapter_t::gno_t *bvgids = new
                         inputAdapter_t::gno_t[coords->getLocalLength()];
  inputAdapter_t::scalar_t *bvcoordarr = new inputAdapter_t::scalar_t[bvsize];
  for (inputAdapter_t::lno_t j = 0; j < bvlen; j++) {
    bvgids[j] = coords->getMap()->getGlobalElement(j);
    for (size_t i = 0; i < bvnvecs; i++) {
      bvcoordarr[idx++] = bvtpetravectors[i][j];
    }
  }

  // my test node type
  typedef Zoltan2::BasicUserTypes<inputAdapter_t::scalar_t,
                                  inputAdapter_t::lno_t,
                                  inputAdapter_t::gno_t,
                                  bv_use_node_t> bvtypes_t;
  typedef Zoltan2::BasicVectorAdapter<bvtypes_t> bvadapter_t;
  std::vector<const inputAdapter_t::scalar_t *> bvcoords(bvnvecs);
  std::vector<int> bvstrides(bvnvecs);
  for (size_t i = 0; i < bvnvecs; i++) {
    bvcoords[i] = &bvcoordarr[i];
    bvstrides[i] = bvnvecs;
  }
  std::vector<const inputAdapter_t::scalar_t *> bvwgts;
  std::vector<int> bvwgtstrides;

  if(numWeightsPerCoord > 0) {
    bvwgts = std::vector<const inputAdapter_t::scalar_t *>(numWeightsPerCoord);
    bvwgtstrides = std::vector<int>(coords->getLocalLength());
    for (size_t i = 0; i < coords->getLocalLength(); i++) {
      bvwgtstrides[i] = numWeightsPerCoord;
    }
    for (int i = 0; i < numWeightsPerCoord; i++) {
      bvwgts[i] = weights[i];
    }
  }

  bvadapter_t bvia(bvlen, bvgids, bvcoords, bvstrides,
                     bvwgts, bvwgtstrides);

  Zoltan2::PartitioningProblem<bvadapter_t> *bvproblem;
  try {
    bvproblem = new Zoltan2::PartitioningProblem<bvadapter_t>(&bvia,
                                               params.getRawPtr(),
                                               comm);
  }
  CATCH_EXCEPTIONS_AND_RETURN("PartitioningProblem()")

  try {
      bvproblem->solve();
  }
  CATCH_EXCEPTIONS_AND_RETURN("solve()")

  int ierr = 0;

  // Compare with MultiVectorAdapter result
  for (inputAdapter_t::lno_t i = 0; i < bvlen; i++) {
    if (problem->getSolution().getPartListView()[i] !=
        bvproblem->getSolution().getPartListView()[i]) {
      std::cout << bvme << " " << i << " "
           << coords->getMap()->getGlobalElement(i) << " " << bvgids[i]
           << ": XMV " << problem->getSolution().getPartListView()[i]
           << "; BMV " << bvproblem->getSolution().getPartListView()[i]
           << "  :  FAIL" << std::endl;
        ++ierr;
    }
    /*  For debugging - plot all success as well
    else {
      std::cout << bvme << " " << i << " "
           << coords->getMap()->getGlobalElement(i) << " " << bvgids[i]
           << ": XMV " << problem->getSolution().getPartListView()[i]
           << "; BMV " << bvproblem->getSolution().getPartListView()[i]
           << "  :  PASS" << std::endl;
    }
    */
  }

  delete [] bvgids;
  delete [] bvcoordarr;
  delete [] bvtpetravectors;
  delete bvproblem;

  if (coords->getGlobalLength() < 40) {
      int len = coords->getLocalLength();
      const inputAdapter_t::part_t *zparts =
            problem->getSolution().getPartListView();
      for (int i = 0; i < len; i++)
          std::cout << comm->getRank()
          << " lid " << i
          << " gid " << coords->getMap()->getGlobalElement(i)
          << " part " << zparts[i] << std::endl;
  }

  return ierr;
}

template<class bv_use_node_t>
int GeometricGenInterface(RCP<const Teuchos::Comm<int> > &comm,
        int numTeams, int numParts, float imbalance,
        std::string paramFile, std::string pqParts,
        std::string pfname,
        int k,
        int migration_check_option,
        int migration_all_to_all_type,
        zscalar_t migration_imbalance_cut_off,
        int migration_processor_assignment_type,
        int migration_doMigration_type,
        bool uvm,
        bool print_details,
        bool test_boxes,
        bool rectilinear,
        int  mj_premigration_option,
        int  mj_premigration_coordinate_cutoff   
)
{
    int ierr = 0;
    Teuchos::ParameterList geoparams("geo params");
    readGeoGenParams(paramFile, geoparams, comm);

    GeometricGen::GeometricGenerator<zscalar_t, zlno_t, zgno_t, znode_t> *gg =
      new GeometricGen::GeometricGenerator<zscalar_t,zlno_t,zgno_t, znode_t>(geoparams,
                                                                      comm);

    int coord_dim = gg->getCoordinateDimension();
    int numWeightsPerCoord = gg->getNumWeights();
    zlno_t numLocalPoints = gg->getNumLocalCoords();
    zgno_t numGlobalPoints = gg->getNumGlobalCoords();
    zscalar_t **coords = new zscalar_t * [coord_dim];
    for(int i = 0; i < coord_dim; ++i){
        coords[i] = new zscalar_t[numLocalPoints];
    }
    gg->getLocalCoordinatesCopy(coords);
    zscalar_t **weight = NULL;
    if (numWeightsPerCoord) {
        weight= new zscalar_t * [numWeightsPerCoord];
        for(int i = 0; i < numWeightsPerCoord; ++i){
            weight[i] = new zscalar_t[numLocalPoints];
        }
        gg->getLocalWeightsCopy(weight);
    }

    delete gg;

    // Run 1st test with MV which always runs UVM on
    RCP<Tpetra::Map<zlno_t, zgno_t, znode_t> > mp = rcp(
                new Tpetra::Map<zlno_t, zgno_t, znode_t>(numGlobalPoints,
                                                      numLocalPoints, 0, comm));
    Teuchos::Array<Teuchos::ArrayView<const zscalar_t> > coordView(coord_dim);
    for (int i=0; i < coord_dim; i++){
        if(numLocalPoints > 0){
            Teuchos::ArrayView<const zscalar_t> a(coords[i], numLocalPoints);
            coordView[i] = a;
        }
        else {
            Teuchos::ArrayView<const zscalar_t> a;
            coordView[i] = a;
        }
    }
    RCP<tMVector_t> tmVector = RCP<tMVector_t>(new
                                   tMVector_t(mp, coordView.view(0, coord_dim),
                                              coord_dim));
    std::vector<const zscalar_t *> weights;
    if(numWeightsPerCoord){
        for (int i = 0; i < numWeightsPerCoord;++i){
          weights.push_back(weight[i]);
        }
    }
    std::vector<int> stride;
    typedef Zoltan2::XpetraMultiVectorAdapter<tMVector_t> inputAdapter_t;
    typedef Zoltan2::EvaluatePartition<inputAdapter_t> quality_t;
    //inputAdapter_t ia(coordsConst);
    inputAdapter_t *ia = new inputAdapter_t(tmVector, weights, stride);

    Teuchos::RCP<Teuchos::ParameterList> params;

    //Teuchos::ParameterList params("test params");
    if(pfname != ""){
        params = Teuchos::getParametersFromXmlFile(pfname);
    }
    else {
        params =RCP<Teuchos::ParameterList>(new Teuchos::ParameterList, true);
    }
/*
    params->set("memory_output_stream" , "std::cout");
    params->set("memory_procs" , 0);
    */

    params->set("timer_output_stream" , "std::cout");

    params->set("algorithm", "multijagged");
    if (test_boxes)
        params->set("mj_keep_part_boxes", true); // bool parameter
    if (rectilinear)
        params->set("rectilinear", true ); // bool parameter

    if(imbalance > 1)
        params->set("imbalance_tolerance", double(imbalance));
    params->set("mj_premigration_option", mj_premigration_option);
    if (mj_premigration_coordinate_cutoff > 0){
        params->set("mj_premigration_coordinate_count", 
                     mj_premigration_coordinate_cutoff);
    }

    if(pqParts != "")
        params->set("mj_parts", pqParts);
    if(numTeams > 0) {
        params->set("mj_num_teams", numTeams);
    }
    if(numParts > 0)
        params->set("num_global_parts", numParts);
    if (k > 0)
        params->set("mj_concurrent_part_count", k);
    if(migration_check_option >= 0)
        params->set("mj_migration_option", migration_check_option);
    if(migration_imbalance_cut_off >= 0)
        params->set("mj_minimum_migration_imbalance",
                    double(migration_imbalance_cut_off));

    Zoltan2::PartitioningProblem<inputAdapter_t> *problem;
    try {
        problem = new Zoltan2::PartitioningProblem<inputAdapter_t>(ia,
                                                   params.getRawPtr(),
                                                   comm);
    }
    CATCH_EXCEPTIONS_AND_RETURN("PartitioningProblem()")

    try {
        problem->solve();
    }
    CATCH_EXCEPTIONS_AND_RETURN("solve()")
    {
      ierr += compareWithBasicVectorAdapterTest<bv_use_node_t>(
        comm, params, problem, tmVector,
        weight, numWeightsPerCoord);
    }

    // create metric object

    RCP<quality_t> metricObject =
      rcp(new quality_t(ia,params.getRawPtr(),comm,&problem->getSolution()));

    if (comm->getRank() == 0){
      metricObject->printMetrics(std::cout);
    }

    problem->printTimers();

    // run pointAssign tests
    if (test_boxes) {
      ierr += run_pointAssign_tests<inputAdapter_t>(problem, tmVector,
        print_details);
      ierr += run_boxAssign_tests<inputAdapter_t>(problem, tmVector);
    }

    if(numWeightsPerCoord){
        for(int i = 0; i < numWeightsPerCoord; ++i) {
          delete [] weight[i];
        }
        delete [] weight;
    }
    if(coord_dim){
        for(int i = 0; i < coord_dim; ++i) {
          delete [] coords[i];
        }
        delete [] coords;
    }

    delete problem;
    delete ia;
    return ierr;
}

template<class bv_use_node_t>
int testFromDataFile(
        RCP<const Teuchos::Comm<int> > &comm,
        int numTeams,
        int numParts,
        float imbalance,
        std::string fname,
        std::string pqParts,
        std::string pfname,
        int k,
        int migration_check_option,
        int migration_all_to_all_type,
        zscalar_t migration_imbalance_cut_off,
        int migration_processor_assignment_type,
        int migration_doMigration_type,
        bool uvm,
        bool print_details,
        bool test_boxes,
        bool rectilinear,
        int  mj_premigration_option, 
        int mj_premigration_coordinate_cutoff
)
{
    int ierr = 0;
    //std::string fname("simple");
    //std::cout << "running " << fname << std::endl;

    UserInputForTests uinput(testDataFilePath, fname, comm, true);

    RCP<tMVector_t> coords = uinput.getUICoordinates();

    typedef Zoltan2::XpetraMultiVectorAdapter<tMVector_t> inputAdapter_t;
    typedef Zoltan2::EvaluatePartition<inputAdapter_t> quality_t;
    inputAdapter_t *ia = new inputAdapter_t(coords);

    Teuchos::RCP <Teuchos::ParameterList> params ;

    //Teuchos::ParameterList params("test params");
    if(pfname != ""){
        params = Teuchos::getParametersFromXmlFile(pfname);
    }
    else {
        params =RCP <Teuchos::ParameterList> (new Teuchos::ParameterList, true);
    }

    //params->set("timer_output_stream" , "std::cout");
    if (test_boxes)
        params->set("mj_keep_part_boxes", true); // bool parameter
    if (rectilinear)
        params->set("rectilinear", true); // bool parameter
    params->set("algorithm", "multijagged");
    if(imbalance > 1){
        params->set("imbalance_tolerance", double(imbalance));
    }

    if(pqParts != ""){
        params->set("mj_parts", pqParts);
    }
    params->set("mj_premigration_option", mj_premigration_option);

    if(numParts > 0){
        params->set("num_global_parts", numParts);
    }
    if (k > 0){
        params->set("mj_concurrent_part_count", k);
    }
    if(migration_check_option >= 0){
        params->set("mj_migration_option", migration_check_option);
    }
    if(migration_imbalance_cut_off >= 0){
        params->set("mj_minimum_migration_imbalance",
                    double (migration_imbalance_cut_off));
    }
    if (mj_premigration_coordinate_cutoff > 0){
        params->set("mj_premigration_coordinate_count",
                     mj_premigration_coordinate_cutoff);
    }

    Zoltan2::PartitioningProblem<inputAdapter_t> *problem;
    try {
        problem = new Zoltan2::PartitioningProblem<inputAdapter_t>(ia,
                                                   params.getRawPtr(),
                                                   comm);
    }
    CATCH_EXCEPTIONS_AND_RETURN("PartitioningProblem()")
    try {
        problem->solve();
    }
    CATCH_EXCEPTIONS_AND_RETURN("solve()")
    {
      ierr += compareWithBasicVectorAdapterTest<bv_use_node_t>(
        comm, params, problem, coords);
    }

    // create metric object

    RCP<quality_t> metricObject =
      rcp(new quality_t(ia,params.getRawPtr(),comm,&problem->getSolution()));

    if (comm->getRank() == 0){
      metricObject->printMetrics(std::cout);
        std::cout << "testFromDataFile is done " << std::endl;
    }

    problem->printTimers();

    // run pointAssign tests
    if (test_boxes) {
      ierr += run_pointAssign_tests<inputAdapter_t>(problem, coords,
        print_details);
      ierr += run_boxAssign_tests<inputAdapter_t>(problem, coords);
    }

    delete problem;
    delete ia;
    return ierr;
}

#ifdef hopper_separate_test

template <typename zscalar_t, typename zlno_t>
void getCoords(zscalar_t **&coords, zlno_t &numLocal, int &dim, string fileName){
    FILE *f = fopen(fileName.c_str(), "r");
    if (f == NULL){
        std::cout << fileName << " cannot be opened" << std::endl;
        std::terminate();
    }
    fscanf(f, "%d", &numLocal);
    fscanf(f, "%d", &dim);
    coords = new zscalar_t *[ dim];
    for (int i = 0; i < dim; ++i){
        coords[i] = new zscalar_t[numLocal];
    }
    for (int i = 0; i < dim; ++i){
        for (zlno_t j = 0; j < numLocal; ++j){
            fscanf(f, "%lf", &(coords[i][j]));
        }
    }
    fclose(f);
    std::cout << "done reading:" << fileName << std::endl;
}

int testFromSeparateDataFiles(
        RCP<const Teuchos::Comm<int> > &comm,
        int numTeams,
        int numParts,
        float imbalance,
        std::string fname,
        std::string pqParts,
        std::string pfname,
        int k,
        int migration_check_option,
        int migration_all_to_all_type,
        zscalar_t migration_imbalance_cut_off,
        int migration_processor_assignment_type,
        int migration_doMigration_type,
        bool uvm,
        bool test_boxes,
        bool rectilinear,
        int  mj_premigration_option
)
{
    //std::string fname("simple");
    //std::cout << "running " << fname << std::endl;

    int ierr = 0;
    int mR = comm->getRank();
    if (mR == 0) std::cout << "size of zscalar_t:" << sizeof(zscalar_t) << std::endl;
    string tFile = fname +"_" + std::to_string(mR) + ".mtx";
    zscalar_t **double_coords;
    zlno_t numLocal = 0;
    int dim = 0;
    getCoords<zscalar_t, zlno_t>(double_coords, numLocal, dim, tFile);
    //UserInputForTests uinput(testDataFilePath, fname, comm, true);
    Teuchos::Array<Teuchos::ArrayView<const zscalar_t> > coordView(dim);
    for (int i=0; i < dim; i++){
        if(numLocal > 0){
            Teuchos::ArrayView<const zscalar_t> a(double_coords[i], numLocal);
            coordView[i] = a;
        } else{
            Teuchos::ArrayView<const zscalar_t> a;
            coordView[i] = a;
        }
    }

    zgno_t numGlobal;
    zgno_t nL = numLocal;
    Teuchos::Comm<int> *tcomm =  (Teuchos::Comm<int> *)comm.getRawPtr();

    reduceAll<int, zgno_t>(
            *tcomm,
            Teuchos::REDUCE_SUM,
            1,
            &nL,
            &numGlobal
    );


    RCP<Tpetra::Map<zlno_t, zgno_t, znode_t> > mp = rcp(
            new Tpetra::Map<zlno_t, zgno_t, znode_t> (numGlobal, numLocal, 0, comm));
    RCP< Tpetra::MultiVector<zscalar_t, zlno_t, zgno_t, znode_t> >coords = RCP< Tpetra::MultiVector<zscalar_t, zlno_t, zgno_t, znode_t> >(
            new Tpetra::MultiVector<zscalar_t, zlno_t, zgno_t, znode_t>( mp, coordView.view(0, dim), dim));


    RCP<const tMVector_t> coordsConst = rcp_const_cast<const tMVector_t>(coords);

    typedef Zoltan2::XpetraMultiVectorAdapter<tMVector_t> inputAdapter_t;
    typedef Zoltan2::EvaluatePartition<inputAdapter_t> quality_t;

    inputAdapter_t *ia = new inputAdapter_t(coordsConst);

    
   
    Teuchos::RCP <Teuchos::ParameterList> params ;

    //Teuchos::ParameterList params("test params");
    if(pfname != ""){
        params = Teuchos::getParametersFromXmlFile(pfname);
    }
    else {
        params =RCP <Teuchos::ParameterList> (new Teuchos::ParameterList, true);
    }

    //params->set("timer_output_stream" , "std::cout");
    params->set("algorithm", "multijagged");
    if(imbalance > 1){
        params->set("imbalance_tolerance", double(imbalance));
    }
 
    params->set("mj_premigration_option", mj_premigration_option);
    if(pqParts != ""){
        params->set("mj_parts", pqParts);
    }
    if(numTeams > 0){
        params->set("mj_num_teams", numTeams);
    }
    if(numParts > 0){
        params->set("num_global_parts", numParts);
    }
    if (k > 0){
        params->set("parallel_part_calculation_count", k);
    }
    if(migration_processor_assignment_type >= 0){
        params->set("migration_processor_assignment_type", migration_processor_assignment_type);
    }
    if(migration_check_option >= 0){
        params->set("migration_check_option", migration_check_option);
    }
    if(migration_all_to_all_type >= 0){
        params->set("migration_all_to_all_type", migration_all_to_all_type);
    }
    if(migration_imbalance_cut_off >= 0){
        params->set("migration_imbalance_cut_off",
                    double (migration_imbalance_cut_off));
    }
    if (migration_doMigration_type >= 0){
        params->set("migration_doMigration_type", int (migration_doMigration_type));
    }
    if (test_boxes)
        params->set("mj_keep_part_boxes", true); // bool parameter
    if (rectilinear)
        params->set("rectilinear", true); // bool parameter

    Zoltan2::PartitioningProblem<inputAdapter_t> *problem;
    try {
        problem =
          new Zoltan2::PartitioningProblem<inputAdapter_t>(ia,
                                                           params.getRawPtr(),
                                                           comm);
    }
    CATCH_EXCEPTIONS_AND_RETURN("PartitioningProblem()")

    try {
        problem->solve();
    }
    CATCH_EXCEPTIONS_AND_RETURN("solve()")
    {
      ierr += compareWithBasicVectorAdapterTest<bv_use_node_t>(
        comm, params, problem, coords);
    }

    if (coordsConst->getGlobalLength() < 40) {
        int len = coordsConst->getLocalLength();
        const inputAdapter_t::part_t *zparts =
                                      problem->getSolution().getPartListView();
        for (int i = 0; i < len; i++)
            std::cout << comm->getRank()
            << " gid " << coords->getMap()->getGlobalElement(i)
            << " part " << zparts[i] << std::endl;
    }

    //create metric object

    RCP<quality_t> metricObject =
      rcp(new quality_t(ia,params.getRawPtr(),comm,&problem->getSolution()));

    if (comm->getRank() == 0){
      metricObject->printMetrics(std::cout);
        std::cout << "testFromDataFile is done " << std::endl;
    }

    problem->printTimers();

    // run pointAssign tests
    if (test_boxes) {
      ierr += run_pointAssign_tests<inputAdapter_t>(problem, coords,
        print_details);
      ierr += run_boxAssign_tests<inputAdapter_t>(problem, coords);
    }

    delete problem;
    delete ia;
    return ierr;
}
#endif



string convert_to_string(char *args){
    string tmp = "";
    for(int i = 0; args[i] != 0; i++)
        tmp += args[i];
    return tmp;
}
bool getArgumentValue(string &argumentid, double &argumentValue, string argumentline){
    std::stringstream stream(std::stringstream::in | std::stringstream::out);
    stream << argumentline;
    getline(stream, argumentid, '=');
    if (stream.eof()){
        return false;
    }
    stream >> argumentValue;
    return true;
}

void getArgVals(
        int narg,
        char **arg,
        int &numTeams,
        int &numParts,
        float &imbalance ,
        string &pqParts,
        int &opt,
        std::string &fname,
        std::string &pfname,
        int &k,
        int &migration_check_option,
        int &migration_all_to_all_type,
        zscalar_t &migration_imbalance_cut_off,
        int &migration_processor_assignment_type,
        int &migration_doMigration_type,
        bool &uvm,
        bool &print_details,
        bool &test_boxes,
        bool &rectilinear,
        int  &mj_premigration_option,
	int &mj_coordinate_cutoff 
)
{
    bool isCset = false;
    bool isPset = false;
    bool isFset = false;
    bool isPFset = false;

    for(int i = 0; i < narg; ++i){
        string tmp = convert_to_string(arg[i]);
        string identifier = "";
        long long int value = -1; double fval = -1;
        if(!getArgumentValue(identifier, fval, tmp)) continue;
        value = (long long int) (fval);

        if(identifier == "W"){
            if(value == 0 || value == 1){
                print_details = (value == 0 ? false : true);
            } else {
                throw "Invalid argument at " + tmp;
            }
        }
        else if(identifier == "UVM"){
            if(value == 0 || value == 1){
                uvm = (value == 0 ? false : true);
            } else {
                throw "Invalid argument at " + tmp;
            }
        }
        else if(identifier == "T"){
            if(value > 0){
                numTeams=value;
            } else {
                throw  "Invalid argument at " + tmp;
            }
        } else if(identifier == "C"){
            if(value > 0){
                numParts=value;
                isCset = true;
            } else {
                throw  "Invalid argument at " + tmp;
            }
        } else if(identifier == "P"){
            std::stringstream stream(std::stringstream::in | std::stringstream::out);
            stream << tmp;
            string ttmp;
            getline(stream, ttmp, '=');
            stream >> pqParts;
            isPset = true;
        }else if(identifier == "I"){
            if(fval > 0){
                imbalance=fval;
            } else {
                throw "Invalid argument at " + tmp;
            }
        } else if(identifier == "MI"){
            if(fval > 0){
                migration_imbalance_cut_off=fval;
            } else {
                throw "Invalid argument at " + tmp;
            }
        } else if(identifier == "MO"){
            if(value >=0 ){
                migration_check_option = value;
            } else {
                throw "Invalid argument at " + tmp;
            }
        } else if(identifier == "AT"){
            if(value >=0 ){
                migration_processor_assignment_type = value;
            } else {
                throw "Invalid argument at " + tmp;
            }
        }

        else if(identifier == "MT"){
            if(value >=0 ){
                migration_all_to_all_type = value;
            } else {
                throw "Invalid argument at " + tmp;
            }
        }
        else if(identifier == "PCC"){
            if(value >=0 ){
                mj_coordinate_cutoff = value;
            } else {
                throw "Invalid argument at " + tmp;
            }
        }

        else if(identifier == "PM"){
            if(value >=0 ){
		mj_premigration_option = value;
            } else {
                throw "Invalid argument at " + tmp;
            }
        }

        else if(identifier == "DM"){
            if(value >=0 ){
                migration_doMigration_type = value;
            } else {
                throw "Invalid argument at " + tmp;
            }
        }
        else if(identifier == "F"){
            std::stringstream stream(std::stringstream::in | std::stringstream::out);
            stream << tmp;
            getline(stream, fname, '=');

            stream >> fname;
            isFset = true;
        }
        else if(identifier == "PF"){
            std::stringstream stream(std::stringstream::in | std::stringstream::out);
            stream << tmp;
            getline(stream, pfname, '=');

            stream >> pfname;
            isPFset = true;
        }

        else if(identifier == "O"){
            if(value >= 0 && value <= 3){
                opt = value;
            } else {
                throw "Invalid argument at " + tmp;
            }
        }
        else if(identifier == "K"){
            if(value >=0 ){
                k = value;
            } else {
                throw "Invalid argument at " + tmp;
            }
        }
        else if(identifier == "TB"){
            if(value >=0 ){
                test_boxes = (value == 0 ? false : true);
            } else {
                throw "Invalid argument at " + tmp;
            }
        }
        else if(identifier == "R"){
            if(value >=0 ){
                rectilinear = (value == 0 ? false : true);
            } else {
                throw "Invalid argument at " + tmp;
            }
        }
        else {
            throw "Invalid argument at " + tmp;
        }

    }
    if(!( (isCset || isPset || isPFset) && isFset)){
        throw "(C || P || PF) && F are mandatory arguments.";
    }

}

void print_usage(char *executable){
    std::cout << "\nUsage:" << std::endl;
    std::cout << executable << " arglist" << std::endl;
    std::cout << "arglist:" << std::endl;
    std::cout << "\tT=numTeams: numTeams > 0" << std::endl;
    std::cout << "\tC=numParts: numParts > 0" << std::endl;
    std::cout << "\tP=MultiJaggedPart: Example: P=512,512" << std::endl;
    std::cout << "\tI=imbalance: Example I=1.03 (ignored for now.)" << std::endl;
    std::cout << "\tF=filePath: When O=0 the path of the coordinate input file, for O>1 the path to the geometric generator parameter file." << std::endl;
    std::cout << "\tO=input option: O=0 for reading coordinate from file, O>0 for generating coordinate from coordinate generator file. Default will run geometric generator." << std::endl;
    std::cout << "\tK=concurrent part calculation input: K>0." << std::endl;
    std::cout << "\tMI=migration_imbalance_cut_off: MI=1.35. " << std::endl;
    std::cout << "\tMT=migration_all_to_all_type: 0 for alltoallv, 1 for Zoltan_Comm, 2 for Zoltan2 Distributor object(Default 1)." << std::endl;
    std::cout << "\tMO=migration_check_option: 0 for decision on imbalance, 1 for forcing migration, >1 for avoiding migration. (Default-0)" << std::endl;
    std::cout << "\tAT=migration_processor_assignment_type. 0-for assigning procs with respect to proc ownment, otherwise, assignment with respect to proc closeness." << std::endl;
    std::cout << "Example:\n" << executable << " P=2,2,2 C=8 F=simple O=0" << std::endl;
}

int main(int narg, char *arg[])
{
    Tpetra::ScopeGuard tscope(&narg, &arg);
    Teuchos::RCP<const Teuchos::Comm<int> > tcomm = Tpetra::getDefaultComm();

    int rank = tcomm->getRank();


    int numTeams = 0; // will use default if not set
    int numParts = -10;
    float imbalance = -1.03;
    int k = -1;

    string pqParts = "";
    int opt = 1;
    std::string fname = "";
    std::string paramFile = "";


    int migration_check_option = -2;
    int migration_all_to_all_type = -1;
    zscalar_t migration_imbalance_cut_off = -1.15;
    int migration_processor_assignment_type = -1;
    int migration_doMigration_type = -1;
    int  mj_premigration_option = 0;
    int mj_premigration_coordinate_cutoff = 0;

    bool uvm = true;
    bool print_details = true;
    bool test_boxes = false;
    bool rectilinear = false;

#ifdef KOKKOS_ENABLE_CUDA
    // make a new node type so we can run BasicVectorAdapter with UVM off
    // The Tpetra MV will still run with UVM on and we'll compare the results.
    // For Serial/OpenMP the 2nd test will be turned off at the CMake level.
    // For CUDA we control uvm on/off with parameter uvm set to 0 or 1.
    // TODO: Probably this should all change eventually so we don't have a node
    // declared like this.
    typedef Kokkos::Compat::KokkosDeviceWrapperNode<
      Kokkos::Cuda, Kokkos::CudaSpace>  uvm_off_node_t;
#endif

    try{
        try {
            getArgVals(
                    narg,
                    arg,
                    numTeams,
                    numParts,
                    imbalance ,
                    pqParts,
                    opt,
                    fname,
                    paramFile,
                    k,
                    migration_check_option,
                    migration_all_to_all_type,
                    migration_imbalance_cut_off,
                    migration_processor_assignment_type,
                    migration_doMigration_type,
                    uvm,
                    print_details,
                    test_boxes,
                    rectilinear,
                    mj_premigration_option, mj_premigration_coordinate_cutoff);
        }
        catch(std::string s){
            if(tcomm->getRank() == 0){
                print_usage(arg[0]);
            }
            throw s;
        }

        catch(char * s){
            if(tcomm->getRank() == 0){
                print_usage(arg[0]);
            }
            throw s;
        }

        int ierr = 0;

        switch (opt){

        case 0:
          if(uvm == true) { // true by default, if not CUDA it should be unset
            ierr = testFromDataFile<znode_t>(tcomm, numTeams, numParts,
                    imbalance, fname, pqParts, paramFile, k,
                    migration_check_option,
                    migration_all_to_all_type,
                    migration_imbalance_cut_off,
                    migration_processor_assignment_type,
                    migration_doMigration_type, uvm, print_details, test_boxes,
                    rectilinear,
                    mj_premigration_option, mj_premigration_coordinate_cutoff);
          }
          else {
#ifdef KOKKOS_ENABLE_CUDA
            ierr = testFromDataFile<uvm_off_node_t>(tcomm, numTeams, numParts,
                    imbalance, fname, pqParts, paramFile, k,
                    migration_check_option,
                    migration_all_to_all_type,
                    migration_imbalance_cut_off,
                    migration_processor_assignment_type,
                    migration_doMigration_type, uvm, print_details, test_boxes,
                    rectilinear,
                    mj_premigration_option, mj_premigration_coordinate_cutoff);
#else
            throw std::logic_error("uvm set off but this is not a cuda test.");
#endif
          }
          break;

#ifdef hopper_separate_test
        case 1:
            // TODO: Note made changes here but did not actually run this
            // method.
            ierr = testFromSeparateDataFiles(tcomm, numTeams, numParts,
                    imbalance, fname, pqParts, paramFile, k,
                    migration_check_option,
                    migration_all_to_all_type,
                    migration_imbalance_cut_off,
                    migration_processor_assignment_type,
                    migration_doMigration_type,uvm, print_details, test_boxes,
                    rectilinear,
                    mj_premigration_option, mj_premigration_coordinate_cutoff);
            break;
#endif
        default:
          if(uvm == true) { // true by default, if not CUDA it should be unset
            ierr = GeometricGenInterface<znode_t>(tcomm, numTeams, numParts,
                    imbalance, fname, pqParts, paramFile, k,
                    migration_check_option,
                    migration_all_to_all_type,
                    migration_imbalance_cut_off,
                    migration_processor_assignment_type,
                    migration_doMigration_type, uvm, print_details, test_boxes,
                    rectilinear,
                    mj_premigration_option, mj_premigration_coordinate_cutoff);
          }
          else {
#ifdef KOKKOS_ENABLE_CUDA
            ierr = GeometricGenInterface<uvm_off_node_t>(tcomm, numTeams,
                    numParts, imbalance, fname, pqParts, paramFile, k,
                    migration_check_option,
                    migration_all_to_all_type,
                    migration_imbalance_cut_off,
                    migration_processor_assignment_type,
                    migration_doMigration_type, uvm, print_details, test_boxes,
                    rectilinear,
                    mj_premigration_option, mj_premigration_coordinate_cutoff);
#else
            throw std::logic_error("uvm set off but this is not a cuda test.");
#endif
          }
          break;
        }

        if (rank == 0) {
            if (ierr == 0) std::cout << "PASS" << std::endl;
            else std::cout << "FAIL" << std::endl;
        }
    }


    catch(std::string &s){
        if (rank == 0)
            std::cerr << s << std::endl;
    }

    catch(char * s){
        if (rank == 0)
            std::cerr << s << std::endl;
    }

    return 0;
}
