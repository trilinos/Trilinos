// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Teuchos_XMLParameterListCoreHelpers.hpp>

#include "MueLu_Exceptions.hpp"
#include "MueLu_MasterList.hpp"

namespace MueLu {

  Teuchos::RCP<const Teuchos::ParameterList> MasterList::List() {
    if (masterList_.is_null()) {
      masterList_ = Teuchos::getParametersFromXmlString(stringList_);
    }

    return masterList_;
  }

  Teuchos::RCP<Teuchos::ParameterList> MasterList::GetProblemSpecificList(std::string const & problemType) {

    if ( (problemType != problemType_) || problemSpecificList_.is_null() ) {
      if (DefaultProblemTypeLists_.find(problemType) != DefaultProblemTypeLists_.end()) {
        problemType_ = problemType;
        problemSpecificList_ = Teuchos::getParametersFromXmlString(DefaultProblemTypeLists_[problemType]);
      } else {
        //TODO provide valid problem types
        TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "Invalid problem type " << problemType << ".");
      }
    }
    return problemSpecificList_;
  }

  Teuchos::RCP<Teuchos::ParameterList> MasterList::masterList_ = Teuchos::null;
  Teuchos::RCP<Teuchos::ParameterList> MasterList::problemSpecificList_ = Teuchos::null;
  std::string                          MasterList::problemType_ = "unknown";
  const std::string                    MasterList::stringList_ =
"<ParameterList name=\"MueLu\">"
  "<Parameter name=\"problem: type\" type=\"string\" value=\"unknown\"/>"
  "<Parameter name=\"verbosity\" type=\"string\" value=\"high\"/>"
  "<Parameter name=\"number of equations\" type=\"int\" value=\"1\"/>"
  "<Parameter name=\"max levels\" type=\"int\" value=\"10\"/>"
  "<Parameter name=\"cycle type\" type=\"string\" value=\"V\"/>"
  "<Parameter name=\"problem: symmetric\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"smoother: pre or post\" type=\"string\" value=\"both\"/>"
  "<Parameter name=\"smoother: type\" type=\"string\" value=\"gs\"/>"
  "<Parameter name=\"smoother: pre type\" type=\"string\" value=\"gs\"/>"
  "<Parameter name=\"smoother: post type\" type=\"string\" value=\"gs\"/>"
  "<ParameterList name=\"smoother: params\"/>"
  "<ParameterList name=\"smoother: pre params\"/>"
  "<ParameterList name=\"smoother: post params\"/>"
  "<Parameter name=\"smoother: overlap\" type=\"int\" value=\"0\"/>"
  "<Parameter name=\"smoother: pre overlap\" type=\"int\" value=\"0\"/>"
  "<Parameter name=\"smoother: post overlap\" type=\"int\" value=\"0\"/>"
  "<Parameter name=\"coarse: max size\" type=\"int\" value=\"2000\"/>"
  "<Parameter name=\"coarse: type\" type=\"string\" value=\"SuperLU\"/>"
  "<ParameterList name=\"coarse: params\"/>"
  "<Parameter name=\"coarse: overlap\" type=\"int\" value=\"0\"/>"
  "<Parameter name=\"aggregation: type\" type=\"string\" value=\"uncoupled\"/>"
  "<Parameter name=\"aggregation: ordering\" type=\"string\" value=\"natural\"/>"
  "<Parameter name=\"aggregation: drop scheme\" type=\"string\" value=\"classical\"/>"
  "<Parameter name=\"aggregation: drop tol\" type=\"double\" value=\"0.0\"/>"
  "<Parameter name=\"aggregation: min agg size\" type=\"int\" value=\"2\"/>"
  "<Parameter name=\"aggregation: max agg size\" type=\"int\" value=\"2147483647\"/>"
  "<Parameter name=\"aggregation: max selected neighbors\" type=\"int\" value=\"0\"/>"
  "<Parameter name=\"aggregation: Dirichlet threshold\" type=\"double\" value=\"0.0\"/>"
  "<Parameter name=\"aggregation: enable phase 1\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"aggregation: enable phase 2a\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"aggregation: enable phase 2b\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"aggregation: enable phase 3\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"aggregation: preserve Dirichlet points\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"aggregation: export visualization data\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"aggregation: mode\" type=\"string\" value=\"old\"/>"
  "<ParameterList name=\"export data\"/>"
  "<Parameter name=\"print initial parameters\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"print unused parameters\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"transpose: use implicit\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"multigrid algorithm\" type=\"string\" value=\"sa\"/>"
  "<Parameter name=\"semicoarsen: coarsen rate\" type=\"int\" value=\"3\"/>"
  "<Parameter name=\"sa: damping factor\" type=\"double\" value=\"1.33333333\"/>"
  "<Parameter name=\"sa: use filtered matrix\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"sa: calculate eigenvalue estimate\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"sa: eigenvalue estimate num iterations\" type=\"int\" value=\"10\"/>"
  "<Parameter name=\"filtered matrix: use lumping\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"filtered matrix: reuse eigenvalue\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"filtered matrix: reuse graph\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"emin: iterative method\" type=\"string\" value=\"cg\"/>"
  "<Parameter name=\"emin: num iterations\" type=\"int\" value=\"2\"/>"
  "<Parameter name=\"emin: num reuse iterations\" type=\"int\" value=\"1\"/>"
  "<Parameter name=\"emin: pattern\" type=\"string\" value=\"AkPtent\"/>"
  "<Parameter name=\"emin: pattern order\" type=\"int\" value=\"1\"/>"
  "<Parameter name=\"repartition: enable\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"repartition: partitioner\" type=\"string\" value=\"zoltan2\"/>"
  "<ParameterList name=\"repartition: params\"/>"
  "<Parameter name=\"repartition: start level\" type=\"int\" value=\"2\"/>"
  "<Parameter name=\"repartition: min rows per proc\" type=\"int\" value=\"800\"/>"
  "<Parameter name=\"repartition: max imbalance\" type=\"double\" value=\"1.2\"/>"
  "<Parameter name=\"repartition: remap parts\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"repartition: remap num values\" type=\"int\" value=\"4\"/>"
  "<Parameter name=\"repartition: print partition distribution\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"repartition: rebalance P and R\" type=\"bool\" value=\"true\"/>"
"</ParameterList>"
;
  std::map<std::string,std::string> MasterList::DefaultProblemTypeLists_ = DefaultProblemStrings<std::string,std::string>
("Poisson-2D",

    "<ParameterList name=\"MueLu\">"
    
            "<Parameter name=\"number of equations\" type=\"int\" value=\"1\"/>"
          
            "<Parameter name=\"smoother: type\" type=\"string\" value=\"CHEBYSHEV\"/>"
          
            "<Parameter name=\"multigrid algorithm\" type=\"string\" value=\"sa\"/>"
          
    "</ParameterList>"
  )
("Poisson-3D",

    "<ParameterList name=\"MueLu\">"
    
            "<Parameter name=\"number of equations\" type=\"int\" value=\"1\"/>"
          
            "<Parameter name=\"smoother: type\" type=\"string\" value=\"CHEBYSHEV\"/>"
          
            "<Parameter name=\"multigrid algorithm\" type=\"string\" value=\"sa\"/>"
          
    "</ParameterList>"
  )
("Elasticity-2D",

    "<ParameterList name=\"MueLu\">"
    
            "<Parameter name=\"number of equations\" type=\"int\" value=\"3\"/>"
          
            "<Parameter name=\"smoother: type\" type=\"string\" value=\"CHEBYSHEV\"/>"
          
            "<Parameter name=\"multigrid algorithm\" type=\"string\" value=\"sa\"/>"
          
    "</ParameterList>"
  )
("Elasticity-3D",

    "<ParameterList name=\"MueLu\">"
    
            "<Parameter name=\"number of equations\" type=\"int\" value=\"6\"/>"
          
            "<Parameter name=\"smoother: type\" type=\"string\" value=\"CHEBYSHEV\"/>"
          
            "<Parameter name=\"multigrid algorithm\" type=\"string\" value=\"sa\"/>"
          
    "</ParameterList>"
  )
("MHD",

    "<ParameterList name=\"MueLu\">"
    
            "<Parameter name=\"smoother: type\" type=\"string\" value=\"SCHWARZ\"/>"
          
    "<ParameterList name=\"smoother: params\">"
    
        "<Parameter name=\"schwarz: overlap level\" type=\"int\" value=\"1\"/>"
        
        "<Parameter name=\"schwarz: combine mode\" type=\"string\" value=\"Zero\"/>"
        
        "<Parameter name=\"schwarz: use reordering\" type=\"bool\" value=\"false\"/>"
        
        "<Parameter name=\"subdomain solver name\" type=\"string\" value=\"RILUK\"/>"
        
    "<ParameterList name=\"subdomain solver parameters\">"
    
        "<Parameter name=\"fact: iluk level-of-fill\" type=\"int\" value=\"0\"/>"
        
        "<Parameter name=\"fact: absolute threshold\" type=\"double\" value=\"0.\"/>"
        
        "<Parameter name=\"fact: relative threshold\" type=\"double\" value=\"1.\"/>"
        
        "<Parameter name=\"fact: relax value\" type=\"double\" value=\"0.\"/>"
        
    "</ParameterList>"
  
    "</ParameterList>"
  
            "<Parameter name=\"aggregation: mode\" type=\"string\" value=\"new\"/>"
          
            "<Parameter name=\"transpose: use implicit\" type=\"bool\" value=\"true\"/>"
          
            "<Parameter name=\"multigrid algorithm\" type=\"string\" value=\"unsmoothed\"/>"
          
    "</ParameterList>"
  )
("ConvectionDiffusion",

    "<ParameterList name=\"MueLu\">"
    
            "<Parameter name=\"problem: symmetric\" type=\"bool\" value=\"false\"/>"
          
            "<Parameter name=\"smoother: type\" type=\"string\" value=\"RELAXATION\"/>"
          
    "<ParameterList name=\"smoother: params\">"
    
        "<Parameter name=\"relaxation: type\" type=\"string\" value=\"Gauss-Seidel\"/>"
        
    "</ParameterList>"
  
            "<Parameter name=\"multigrid algorithm\" type=\"string\" value=\"pg\"/>"
          
            "<Parameter name=\"sa: use filtered matrix\" type=\"bool\" value=\"true\"/>"
          
    "</ParameterList>"
  )
;

}

