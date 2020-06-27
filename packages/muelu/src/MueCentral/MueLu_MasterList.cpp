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

   std::string MasterList::interpretParameterName(const std::string& name, const std::string& value) {

    // used to concatenate the return string
    std::stringstream ss;

    // put in short cuts here!

    if (name == "verbosity") {
      std::string verb = "none";
      if (value == "\"0\"") verb = "none";
      if (value == "\"1\"" || value == "\"2\"" || value == "\"3\"") verb = "low";
      if (value == "\"4\"" || value == "\"5\"" || value == "\"6\"") verb = "medium";
      if (value == "\"7\"" || value == "\"8\"") verb = "high";
      if (value == "\"9\"") verb = "extreme";
      if (value == "\"10\"") verb = "test";
      verb = "\"" + verb + "\"";
      ss << "<Parameter name=\"verbosity\" type=\"string\" value=" << verb << "/>";
      return ss.str();
    }

    if (name == "cycle type") {
      std::stringstream temp1; temp1 << "\"" << "MGV" << "\"";
      std::stringstream temp2; temp2 << "\"" << "MGV" << "\"";
      if (value == temp1.str() ) { ss << "<Parameter name=\"cycle type\" type=\"string\" value=\"V\"/>"; }
      else if (value == temp2.str()) { ss << "<Parameter name=\"cycle type\" type=\"string\" value=\"W\"/>"; }
      else TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "MasterList::interpretParameterName, Line " << __LINE__ << ". "
                                           << "The parameter " << value << " is not supported by MueLu.");
      return ss.str();
    }

    // energy minimization is enabled
    if (name == "multigrid algorithm") {
      std::stringstream temp; temp << "\"" << "1" << "\"";
      if (value == temp.str() ) { ss << "<Parameter name=\"multigrid algorithm\" type=\"string\" value=\"pg\"/>"; return ss.str(); }
    }

    if (name == "repartition: enable") {
      std::stringstream temp1; temp1 << "\"" << "1" << "\"";
      if (value == temp1.str()) {
        RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
        *out << "WARNING: repartitioning in MueLu is different to MLs. Please refer to the MueLu users Manual for more information." << std::endl;
      }
    }

    // put in auto-generated code here


    if (name == "output filename") { ss << "<Parameter name=\"output filename\" type=\"string\" value=" << value << "/>"; return ss.str(); }      
    if (name == "number of equations") { ss << "<Parameter name=\"number of equations\" type=\"int\" value=" << value << "/>"; return ss.str(); }      
    if (name == "max levels") { ss << "<Parameter name=\"max levels\" type=\"int\" value=" << value << "/>"; return ss.str(); }      
    if (name == "W cycle start level") { ss << "<Parameter name=\"W cycle start level\" type=\"int\" value=" << value << "/>"; return ss.str(); }      
    if (name == "coarse grid correction scaling factor") { ss << "<Parameter name=\"coarse grid correction scaling factor\" type=\"double\" value=" << value << "/>"; return ss.str(); }      
    if (name == "fuse prolongation and update") { ss << "<Parameter name=\"fuse prolongation and update\" type=\"bool\" value=" << value << "/>"; return ss.str(); }      
    if (name == "number of vectors") { ss << "<Parameter name=\"number of vectors\" type=\"int\" value=" << value << "/>"; return ss.str(); }      
    if (name == "problem: symmetric") { ss << "<Parameter name=\"problem: symmetric\" type=\"bool\" value=" << value << "/>"; return ss.str(); }      
    if (name == "hierarchy label") { ss << "<Parameter name=\"hierarchy label\" type=\"string\" value=" << value << "/>"; return ss.str(); }      
    if (name == "aggregation: drop tol") { ss << "<Parameter name=\"aggregation: drop tol\" type=\"double\" value=" << value << "/>"; return ss.str(); }      
    if (name == "print initial parameters") { ss << "<Parameter name=\"print initial parameters\" type=\"bool\" value=" << value << "/>"; return ss.str(); }      
    if (name == "print unused parameters") { ss << "<Parameter name=\"print unused parameters\" type=\"bool\" value=" << value << "/>"; return ss.str(); }      
    if (name == "sa: damping factor") { ss << "<Parameter name=\"sa: damping factor\" type=\"double\" value=" << value << "/>"; return ss.str(); }      
    if (name == "sa: eigenvalue estimate num iterations") { ss << "<Parameter name=\"sa: eigenvalue estimate num iterations\" type=\"int\" value=" << value << "/>"; return ss.str(); }      
    if (name == "pcoarsen: element") { ss << "<Parameter name=\"pcoarsen: element\" type=\"string\" value=" << value << "/>"; return ss.str(); }      
    if (name == "pcoarsen: schedule") { ss << "<Parameter name=\"pcoarsen: schedule\" type=\"string\" value=" << value << "/>"; return ss.str(); }      
    if (name == "pcoarsen: hi basis") { ss << "<Parameter name=\"pcoarsen: hi basis\" type=\"string\" value=" << value << "/>"; return ss.str(); }      
    if (name == "pcoarsen: lo basis") { ss << "<Parameter name=\"pcoarsen: lo basis\" type=\"string\" value=" << value << "/>"; return ss.str(); }      
    if (name == "smoother: neighborhood type") { ss << "<Parameter name=\"smoother: neighborhood type\" type=\"string\" value=" << value << "/>"; return ss.str(); }      
    if (name == "tentative: calculate qr") { ss << "<Parameter name=\"tentative: calculate qr\" type=\"bool\" value=" << value << "/>"; return ss.str(); }      
    if (name == "tentative: constant column sums") { ss << "<Parameter name=\"tentative: constant column sums\" type=\"bool\" value=" << value << "/>"; return ss.str(); }      
    if (name == "repartition: enable") { ss << "<Parameter name=\"repartition: enable\" type=\"bool\" value=" << value << "/>"; return ss.str(); }      
    if (name == "repartition: start level") { ss << "<Parameter name=\"repartition: start level\" type=\"int\" value=" << value << "/>"; return ss.str(); }      
    if (name == "repartition: use map") { ss << "<Parameter name=\"repartition: use map\" type=\"bool\" value=" << value << "/>"; return ss.str(); }      
    if (name == "repartition: node repartition level") { ss << "<Parameter name=\"repartition: node repartition level\" type=\"int\" value=" << value << "/>"; return ss.str(); }      
    if (name == "repartition: node id") { ss << "<Parameter name=\"repartition: node id\" type=\"int\" value=" << value << "/>"; return ss.str(); }      
    if (name == "repartition: min rows per proc") { ss << "<Parameter name=\"repartition: min rows per proc\" type=\"int\" value=" << value << "/>"; return ss.str(); }      
    if (name == "repartition: max imbalance") { ss << "<Parameter name=\"repartition: max imbalance\" type=\"double\" value=" << value << "/>"; return ss.str(); }      
    if (name == "use external multigrid package") { ss << "<Parameter name=\"use external multigrid package\" type=\"string\" value=" << value << "/>"; return ss.str(); }      
    if (name == "refmaxwell: mode") { ss << "<Parameter name=\"refmaxwell: mode\" type=\"string\" value=" << value << "/>"; return ss.str(); }      
    if (name == "refmaxwell: disable addon") { ss << "<Parameter name=\"refmaxwell: disable addon\" type=\"bool\" value=" << value << "/>"; return ss.str(); }      
    if (name == "refmaxwell: use as preconditioner") { ss << "<Parameter name=\"refmaxwell: use as preconditioner\" type=\"bool\" value=" << value << "/>"; return ss.str(); }      
    if (name == "refmaxwell: dump matrices") { ss << "<Parameter name=\"refmaxwell: dump matrices\" type=\"bool\" value=" << value << "/>"; return ss.str(); }      
    if (name == "refmaxwell: subsolves on subcommunicators") { ss << "<Parameter name=\"refmaxwell: subsolves on subcommunicators\" type=\"bool\" value=" << value << "/>"; return ss.str(); }      
    if (name == "refmaxwell: ratio AH / A22 subcommunicators") { ss << "<Parameter name=\"refmaxwell: ratio AH / A22 subcommunicators\" type=\"double\" value=" << value << "/>"; return ss.str(); }      
    return "";
  }

  Teuchos::RCP<Teuchos::ParameterList> MasterList::masterList_ = Teuchos::null;
  Teuchos::RCP<Teuchos::ParameterList> MasterList::problemSpecificList_ = Teuchos::null;
  std::string                          MasterList::problemType_ = "unknown";
  const std::string                    MasterList::stringList_ =
"<ParameterList name=\"MueLu\">"
  "<Parameter name=\"problem: type\" type=\"string\" value=\"unknown\"/>"
  "<Parameter name=\"verbosity\" type=\"string\" value=\"high\"/>"
  "<Parameter name=\"output filename\" type=\"string\" value=\"\"/>"
  "<Parameter name=\"number of equations\" type=\"int\" value=\"1\"/>"
  "<Parameter name=\"max levels\" type=\"int\" value=\"10\"/>"
  "<Parameter name=\"cycle type\" type=\"string\" value=\"V\"/>"
  "<Parameter name=\"W cycle start level\" type=\"int\" value=\"0\"/>"
  "<Parameter name=\"coarse grid correction scaling factor\" type=\"double\" value=\"1.0\"/>"
  "<Parameter name=\"fuse prolongation and update\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"number of vectors\" type=\"int\" value=\"1\"/>"
  "<Parameter name=\"problem: symmetric\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"xml parameter file\" type=\"string\" value=\"\"/>"
  "<Parameter name=\"parameterlist: syntax\" type=\"string\" value=\"muelu\"/>"
  "<Parameter name=\"hierarchy label\" type=\"string\" value=\"\"/>"
  "<ParameterList name=\"matvec params\"/>"
  "<Parameter name=\"smoother: pre or post\" type=\"string\" value=\"both\"/>"
  "<Parameter name=\"smoother: type\" type=\"string\" value=\"RELAXATION\"/>"
  "<Parameter name=\"smoother: pre type\" type=\"string\" value=\"RELAXATION\"/>"
  "<Parameter name=\"smoother: post type\" type=\"string\" value=\"RELAXATION\"/>"
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
  "<Parameter name=\"aggregation: mode\" type=\"string\" value=\"uncoupled\"/>"
  "<Parameter name=\"aggregation: ordering\" type=\"string\" value=\"natural\"/>"
  "<Parameter name=\"aggregation: drop scheme\" type=\"string\" value=\"classical\"/>"
  "<Parameter name=\"aggregation: number of random vectors\" type=\"int\" value=\"10\"/>"
  "<Parameter name=\"aggregation: number of times to pre or post smooth\" type=\"int\" value=\"10\"/>"
  "<Parameter name=\"aggregation: penalty parameters\" type=\"Array(double)\" value=\"{12.,-.2,0,0,0} \"/>"
  "<Parameter name=\"aggregation: distance laplacian algo\" type=\"string\" value=\"default\"/>"
  "<Parameter name=\"aggregation: classical algo\" type=\"string\" value=\"default\"/>"
  "<Parameter name=\"aggregation: drop tol\" type=\"double\" value=\"0.0\"/>"
  "<Parameter name=\"aggregation: min agg size\" type=\"int\" value=\"2\"/>"
  "<Parameter name=\"aggregation: max agg size\" type=\"int\" value=\"-1\"/>"
  "<Parameter name=\"aggregation: compute aggregate qualities\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"aggregation: brick x size\" type=\"int\" value=\"2\"/>"
  "<Parameter name=\"aggregation: brick y size\" type=\"int\" value=\"2\"/>"
  "<Parameter name=\"aggregation: brick z size\" type=\"int\" value=\"2\"/>"
  "<Parameter name=\"aggregation: brick x Dirichlet\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"aggregation: brick y Dirichlet\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"aggregation: brick z Dirichlet\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"aggregation: max selected neighbors\" type=\"int\" value=\"0\"/>"
  "<Parameter name=\"aggregation: Dirichlet threshold\" type=\"double\" value=\"0.0\"/>"
  "<Parameter name=\"aggregation: deterministic\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"aggregation: coloring algorithm\" type=\"string\" value=\"serial\"/>"
  "<Parameter name=\"aggregation: enable phase 1\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"aggregation: enable phase 2a\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"aggregation: enable phase 2b\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"aggregation: enable phase 3\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"aggregation: phase2a include root\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"aggregation: error on nodes with no on-rank neighbors\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"aggregation: phase3 avoid singletons\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"aggregation: allow empty prolongator columns\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"aggregation: preserve Dirichlet points\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"aggregation: allow user-specified singletons\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"aggregation: use interface aggregation\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"aggregation: export visualization data\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"aggregation: output filename\" type=\"string\" value=\"\"/>"
  "<Parameter name=\"aggregation: output file: time step\" type=\"int\" value=\"0\"/>"
  "<Parameter name=\"aggregation: output file: iter\" type=\"int\" value=\"0\"/>"
  "<Parameter name=\"aggregation: output file: agg style\" type=\"string\" value=\"Point Cloud\"/>"
  "<Parameter name=\"aggregation: output file: fine graph edges\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"aggregation: output file: coarse graph edges\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"aggregation: output file: build colormap\" type=\"bool\" value=\"false\"/>"
  "<ParameterList name=\"aggregation: params\"/>"
  "<ParameterList name=\"strength-of-connection: params\"/>"
  "<Parameter name=\"aggregation: mesh layout\" type=\"string\" value=\"Global Lexicographic\"/>"
  "<Parameter name=\"aggregation: output type\" type=\"string\" value=\"Aggregates\"/>"
  "<Parameter name=\"aggregation: coarsening rate\" type=\"string\" value=\"{3}\"/>"
  "<Parameter name=\"aggregation: number of spatial dimensions\" type=\"int\" value=\"3\"/>"
  "<Parameter name=\"aggregation: coarsening order\" type=\"int\" value=\"0\"/>"
  "<Parameter name=\"aggregation: pairwise: size\" type=\"int\" value=\"8\"/>"
  "<Parameter name=\"aggregation: pairwise: tie threshold\" type=\"double\" value=\"1e-6\"/>"
  "<Parameter name=\"aggregate qualities: check symmetry\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"aggregate qualities: good aggregate threshold\" type=\"double\" value=\"100.0\"/>"
  "<Parameter name=\"aggregate qualities: file output\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"aggregate qualities: file base\" type=\"string\" value=\"agg_qualities\"/>"
  "<Parameter name=\"aggregate qualities: algorithm\" type=\"string\" value=\"forward\"/>"
  "<Parameter name=\"aggregate qualities: zero threshold\" type=\"double\" value=\"1e-12\"/>"
  "<Parameter name=\"aggregate qualities: percentiles\" type=\"Array(double)\" value=\"{}\"/>"
  "<ParameterList name=\"export data\"/>"
  "<Parameter name=\"print initial parameters\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"print unused parameters\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"transpose: use implicit\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"restriction: scale nullspace\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"use kokkos refactor\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"synchronize factory timers\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"rap: triple product\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"multigrid algorithm\" type=\"string\" value=\"sa\"/>"
  "<Parameter name=\"toggle: mode\" type=\"string\" value=\"semicoarsen\"/>"
  "<Parameter name=\"semicoarsen: coarsen rate\" type=\"int\" value=\"3\"/>"
  "<Parameter name=\"semicoarsen: number of levels\" type=\"int\" value=\"3\"/>"
  "<Parameter name=\"linedetection: orientation\" type=\"string\" value=\"vertical\"/>"
  "<Parameter name=\"linedetection: num layers\" type=\"int\" value=\"-1\"/>"
  "<Parameter name=\"sa: damping factor\" type=\"double\" value=\"1.33\"/>"
  "<Parameter name=\"sa: use filtered matrix\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"sa: calculate eigenvalue estimate\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"sa: eigenvalue estimate num iterations\" type=\"int\" value=\"10\"/>"
  "<Parameter name=\"interp: interpolation order\" type=\"int\" value=\"1\"/>"
  "<Parameter name=\"interp: build coarse coordinates\" type=\"bool\" value=\"true\"/>"
  "<ParameterList name=\"transfer: params\"/>"
  "<Parameter name=\"pcoarsen: element\" type=\"string\" value=\"\"/>"
  "<Parameter name=\"pcoarsen: schedule\" type=\"string\" value=\"\"/>"
  "<Parameter name=\"pcoarsen: hi basis\" type=\"string\" value=\"\"/>"
  "<Parameter name=\"pcoarsen: lo basis\" type=\"string\" value=\"\"/>"
  "<Parameter name=\"smoother: neighborhood type\" type=\"string\" value=\"\"/>"
  "<Parameter name=\"filtered matrix: use lumping\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"filtered matrix: reuse eigenvalue\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"filtered matrix: reuse graph\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"emin: iterative method\" type=\"string\" value=\"cg\"/>"
  "<Parameter name=\"emin: num iterations\" type=\"int\" value=\"2\"/>"
  "<Parameter name=\"emin: num reuse iterations\" type=\"int\" value=\"1\"/>"
  "<Parameter name=\"emin: pattern\" type=\"string\" value=\"AkPtent\"/>"
  "<Parameter name=\"emin: pattern order\" type=\"int\" value=\"1\"/>"
  "<Parameter name=\"emin: use filtered matrix\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"tentative: calculate qr\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"tentative: constant column sums\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"tentative: build coarse coordinates\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"repartition: enable\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"repartition: partitioner\" type=\"string\" value=\"zoltan2\"/>"
  "<ParameterList name=\"repartition: params\"/>"
  "<Parameter name=\"repartition: start level\" type=\"int\" value=\"2\"/>"
  "<Parameter name=\"repartition: use map\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"repartition: node repartition level\" type=\"int\" value=\"-1\"/>"
  "<Parameter name=\"repartition: node id\" type=\"int\" value=\"-1\"/>"
  "<Parameter name=\"repartition: min rows per proc\" type=\"int\" value=\"800\"/>"
  "<Parameter name=\"repartition: target rows per proc\" type=\"int\" value=\"0\"/>"
  "<Parameter name=\"repartition: min rows per thread\" type=\"int\" value=\"0\"/>"
  "<Parameter name=\"repartition: target rows per thread\" type=\"int\" value=\"0\"/>"
  "<Parameter name=\"repartition: max imbalance\" type=\"double\" value=\"1.2\"/>"
  "<Parameter name=\"repartition: remap parts\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"repartition: remap num values\" type=\"int\" value=\"4\"/>"
  "<Parameter name=\"repartition: remap accept partition\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"repartition: print partition distribution\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"repartition: rebalance P and R\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"repartition: rebalance Nullspace\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"repartition: use subcommunicators\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"rap: relative diagonal floor\" type=\"Array(double)\" value=\"{}\"/>"
  "<Parameter name=\"rap: fix zero diagonals\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"rap: fix zero diagonals threshold\" type=\"double\" value=\"0.\"/>"
  "<Parameter name=\"rap: shift\" type=\"double\" value=\"0.0\"/>"
  "<Parameter name=\"rap: shift diagonal M\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"rap: shift low storage\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"rap: shift array\" type=\"Array(double)\" value=\"{}\"/>"
  "<Parameter name=\"rap: cfl array\" type=\"Array(double)\" value=\"{}\"/>"
  "<Parameter name=\"rap: algorithm\" type=\"string\" value=\"galerkin\"/>"
  "<ParameterList name=\"matrixmatrix: kernel params\"/>"
  "<Parameter name=\"matrixmatrix: kernel params:MM_TAFC_OptimizationCoreCount \" type=\"int\" value=\"3000 \"/>"
  "<Parameter name=\"isMatrixMatrix_TransferAndFillComplete \" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"reuse: type\" type=\"string\" value=\"none\"/>"
  "<Parameter name=\"use external multigrid package\" type=\"string\" value=\"none\"/>"
  "<ParameterList name=\"amgx:params\"/>"
  "<Parameter name=\"debug: graph level\" type=\"int\" value=\"-1\"/>"
  "<Parameter name=\"refmaxwell: mode\" type=\"string\" value=\"additive\"/>"
  "<Parameter name=\"refmaxwell: disable addon\" type=\"bool\" value=\"true\"/>"
  "<ParameterList name=\"refmaxwell: 11list\"/>"
  "<ParameterList name=\"refmaxwell: 22list\"/>"
  "<Parameter name=\"refmaxwell: use as preconditioner\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"refmaxwell: dump matrices\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"refmaxwell: subsolves on subcommunicators\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"refmaxwell: ratio AH / A22 subcommunicators\" type=\"double\" value=\"1.0\"/>"
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
("Poisson-2D-complex",

    "<ParameterList name=\"MueLu\">"
    
            "<Parameter name=\"number of equations\" type=\"int\" value=\"1\"/>"
          
            "<Parameter name=\"smoother: type\" type=\"string\" value=\"RELAXATION\"/>"
          
    "<ParameterList name=\"smoother: params\">"
    
        "<Parameter name=\"relaxation: type\" type=\"string\" value=\"Symmetric Gauss-Seidel\"/>"
        
    "</ParameterList>"
  
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
("Poisson-3D-complex",

    "<ParameterList name=\"MueLu\">"
    
            "<Parameter name=\"number of equations\" type=\"int\" value=\"1\"/>"
          
            "<Parameter name=\"smoother: type\" type=\"string\" value=\"RELAXATION\"/>"
          
    "<ParameterList name=\"smoother: params\">"
    
        "<Parameter name=\"relaxation: type\" type=\"string\" value=\"Symmetric Gauss-Seidel\"/>"
        
    "</ParameterList>"
  
            "<Parameter name=\"multigrid algorithm\" type=\"string\" value=\"sa\"/>"
          
    "</ParameterList>"
  )
("Elasticity-2D",

    "<ParameterList name=\"MueLu\">"
    
            "<Parameter name=\"number of equations\" type=\"int\" value=\"2\"/>"
          
            "<Parameter name=\"smoother: type\" type=\"string\" value=\"CHEBYSHEV\"/>"
          
            "<Parameter name=\"multigrid algorithm\" type=\"string\" value=\"sa\"/>"
          
    "</ParameterList>"
  )
("Elasticity-2D-complex",

    "<ParameterList name=\"MueLu\">"
    
            "<Parameter name=\"number of equations\" type=\"int\" value=\"2\"/>"
          
            "<Parameter name=\"smoother: type\" type=\"string\" value=\"RELAXATION\"/>"
          
    "<ParameterList name=\"smoother: params\">"
    
        "<Parameter name=\"relaxation: type\" type=\"string\" value=\"Symmetric Gauss-Seidel\"/>"
        
    "</ParameterList>"
  
            "<Parameter name=\"multigrid algorithm\" type=\"string\" value=\"sa\"/>"
          
    "</ParameterList>"
  )
("Elasticity-3D",

    "<ParameterList name=\"MueLu\">"
    
            "<Parameter name=\"number of equations\" type=\"int\" value=\"3\"/>"
          
            "<Parameter name=\"smoother: type\" type=\"string\" value=\"CHEBYSHEV\"/>"
          
            "<Parameter name=\"multigrid algorithm\" type=\"string\" value=\"sa\"/>"
          
    "</ParameterList>"
  )
("Elasticity-3D-complex",

    "<ParameterList name=\"MueLu\">"
    
            "<Parameter name=\"number of equations\" type=\"int\" value=\"3\"/>"
          
            "<Parameter name=\"smoother: type\" type=\"string\" value=\"RELAXATION\"/>"
          
    "<ParameterList name=\"smoother: params\">"
    
        "<Parameter name=\"relaxation: type\" type=\"string\" value=\"Symmetric Gauss-Seidel\"/>"
        
    "</ParameterList>"
  
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
          
            "<Parameter name=\"emin: use filtered matrix\" type=\"bool\" value=\"true\"/>"
          
    "</ParameterList>"
  )
;
  std::map<std::string,std::string> MasterList::ML2MueLuLists_ = DefaultProblemStrings<std::string,std::string>

         ("default values","problem: type")
      
         ("ML output","verbosity")
      
         ("output filename","output filename")
      
         ("PDE equations","number of equations")
      
         ("max levels","max levels")
      
         ("prec type","cycle type")
      
         ("W cycle start level","W cycle start level")
      
         ("coarse grid correction scaling factor","coarse grid correction scaling factor")
      
         ("fuse prolongation and update","fuse prolongation and update")
      
         ("number of vectors","number of vectors")
      
         ("problem: symmetric","problem: symmetric")
      
         ("xml parameter file","xml parameter file")
      
         ("parameterlist: syntax","parameterlist: syntax")
      
         ("hierarchy label","hierarchy label")
      
         ("matvec params","matvec params")
      
         ("smoother: pre or post","smoother: pre or post")
      
         ("smoother: type","smoother: type")
      
         ("smoother: pre type","smoother: pre type")
      
         ("smoother: post type","smoother: post type")
      
         ("smoother: params","smoother: params")
      
         ("smoother: pre params","smoother: pre params")
      
         ("smoother: post params","smoother: post params")
      
         ("smoother: overlap","smoother: overlap")
      
         ("smoother: pre overlap","smoother: pre overlap")
      
         ("smoother: post overlap","smoother: post overlap")
      
         ("max size","coarse: max size")
      
         ("coarse: type","coarse: type")
      
         ("coarse: params","coarse: params")
      
         ("coarse: overlap","coarse: overlap")
      
         ("aggregation: type","aggregation: type")
      
         ("aggregation: mode","aggregation: mode")
      
         ("aggregation: ordering","aggregation: ordering")
      
         ("aggregation: drop scheme","aggregation: drop scheme")
      
         ("aggregation: number of random vectors","aggregation: number of random vectors")
      
         ("aggregation: number of times to pre or post smooth","aggregation: number of times to pre or post smooth")
      
         ("aggregation: penalty parameters","aggregation: penalty parameters")
      
         ("aggregation: distance laplacian algo","aggregation: distance laplacian algo")
      
         ("aggregation: classical algo","aggregation: classical algo")
      
         ("aggregation: threshold","aggregation: drop tol")
      
         ("aggregation: min agg size","aggregation: min agg size")
      
         ("aggregation: max agg size","aggregation: max agg size")
      
         ("aggregation: compute aggregate qualities","aggregation: compute aggregate qualities")
      
         ("aggregation: brick x size","aggregation: brick x size")
      
         ("aggregation: brick y size","aggregation: brick y size")
      
         ("aggregation: brick z size","aggregation: brick z size")
      
         ("aggregation: brick x Dirichlet","aggregation: brick x Dirichlet")
      
         ("aggregation: brick y Dirichlet","aggregation: brick y Dirichlet")
      
         ("aggregation: brick z Dirichlet","aggregation: brick z Dirichlet")
      
         ("aggregation: max selected neighbors","aggregation: max selected neighbors")
      
         ("aggregation: Dirichlet threshold","aggregation: Dirichlet threshold")
      
         ("aggregation: deterministic","aggregation: deterministic")
      
         ("aggregation: coloring algorithm","aggregation: coloring algorithm")
      
         ("aggregation: enable phase 1","aggregation: enable phase 1")
      
         ("aggregation: enable phase 2a","aggregation: enable phase 2a")
      
         ("aggregation: enable phase 2b","aggregation: enable phase 2b")
      
         ("aggregation: enable phase 3","aggregation: enable phase 3")
      
         ("aggregation: phase2a include root","aggregation: phase2a include root")
      
         ("aggregation: error on nodes with no on-rank neighbors","aggregation: error on nodes with no on-rank neighbors")
      
         ("aggregation: phase3 avoid singletons","aggregation: phase3 avoid singletons")
      
         ("aggregation: allow empty prolongator columns","aggregation: allow empty prolongator columns")
      
         ("aggregation: preserve Dirichlet points","aggregation: preserve Dirichlet points")
      
         ("aggregation: allow user-specified singletons","aggregation: allow user-specified singletons")
      
         ("aggregation: use interface aggregation","aggregation: use interface aggregation")
      
         ("aggregation: export visualization data","aggregation: export visualization data")
      
         ("aggregation: output filename","aggregation: output filename")
      
         ("aggregation: output file: time step","aggregation: output file: time step")
      
         ("aggregation: output file: iter","aggregation: output file: iter")
      
         ("aggregation: output file: agg style","aggregation: output file: agg style")
      
         ("aggregation: output file: fine graph edges","aggregation: output file: fine graph edges")
      
         ("aggregation: output file: coarse graph edges","aggregation: output file: coarse graph edges")
      
         ("aggregation: output file: build colormap","aggregation: output file: build colormap")
      
         ("aggregation: params","aggregation: params")
      
         ("strength-of-connection: params","strength-of-connection: params")
      
         ("aggregation: mesh layout","aggregation: mesh layout")
      
         ("aggregation: output type","aggregation: output type")
      
         ("aggregation: coarsening rate","aggregation: coarsening rate")
      
         ("aggregation: number of spatial dimensions","aggregation: number of spatial dimensions")
      
         ("aggregation: coarsening order","aggregation: coarsening order")
      
         ("aggregation: pairwise: size","aggregation: pairwise: size")
      
         ("aggregation: pairwise: tie threshold","aggregation: pairwise: tie threshold")
      
         ("aggregate qualities: check symmetry","aggregate qualities: check symmetry")
      
         ("aggregate qualities: good aggregate threshold","aggregate qualities: good aggregate threshold")
      
         ("aggregate qualities: file output","aggregate qualities: file output")
      
         ("aggregate qualities: file base","aggregate qualities: file base")
      
         ("aggregate qualities: algorithm","aggregate qualities: algorithm")
      
         ("aggregate qualities: zero threshold","aggregate qualities: zero threshold")
      
         ("aggregate qualities: percentiles","aggregate qualities: percentiles")
      
         ("export data","export data")
      
         ("ML print initial list","print initial parameters")
      
         ("print unused","print unused parameters")
      
         ("transpose: use implicit","transpose: use implicit")
      
         ("restriction: scale nullspace","restriction: scale nullspace")
      
         ("use kokkos refactor","use kokkos refactor")
      
         ("synchronize factory timers","synchronize factory timers")
      
         ("rap: triple product","rap: triple product")
      
         ("energy minimization: enable","multigrid algorithm")
      
         ("toggle: mode","toggle: mode")
      
         ("semicoarsen: coarsen rate","semicoarsen: coarsen rate")
      
         ("semicoarsen: number of levels","semicoarsen: number of levels")
      
         ("linedetection: orientation","linedetection: orientation")
      
         ("linedetection: num layers","linedetection: num layers")
      
         ("aggregation: damping factor","sa: damping factor")
      
         ("sa: use filtered matrix","sa: use filtered matrix")
      
         ("sa: calculate eigenvalue estimate","sa: calculate eigenvalue estimate")
      
         ("eigen-analysis: iterations","sa: eigenvalue estimate num iterations")
      
         ("interp: interpolation order","interp: interpolation order")
      
         ("interp: build coarse coordinates","interp: build coarse coordinates")
      
         ("transfer: params","transfer: params")
      
         ("pcoarsen: element","pcoarsen: element")
      
         ("pcoarsen: schedule","pcoarsen: schedule")
      
         ("pcoarsen: hi basis","pcoarsen: hi basis")
      
         ("pcoarsen: lo basis","pcoarsen: lo basis")
      
         ("smoother: neighborhood type","smoother: neighborhood type")
      
         ("filtered matrix: use lumping","filtered matrix: use lumping")
      
         ("filtered matrix: reuse eigenvalue","filtered matrix: reuse eigenvalue")
      
         ("filtered matrix: reuse graph","filtered matrix: reuse graph")
      
         ("emin: iterative method","emin: iterative method")
      
         ("emin: num iterations","emin: num iterations")
      
         ("emin: num reuse iterations","emin: num reuse iterations")
      
         ("emin: pattern","emin: pattern")
      
         ("emin: pattern order","emin: pattern order")
      
         ("emin: use filtered matrix","emin: use filtered matrix")
      
         ("tentative: calculate qr","tentative: calculate qr")
      
         ("tentative: constant column sums","tentative: constant column sums")
      
         ("tentative: build coarse coordinates","tentative: build coarse coordinates")
      
         ("repartition: enable","repartition: enable")
      
         ("repartition: partitioner","repartition: partitioner")
      
         ("repartition: params","repartition: params")
      
         ("repartition: start level","repartition: start level")
      
         ("repartition: use map","repartition: use map")
      
         ("repartition: node repartition level","repartition: node repartition level")
      
         ("repartition: node id","repartition: node id")
      
         ("repartition: min per proc","repartition: min rows per proc")
      
         ("repartition: target rows per proc","repartition: target rows per proc")
      
         ("repartition: min rows per thread","repartition: min rows per thread")
      
         ("repartition: target rows per thread","repartition: target rows per thread")
      
         ("repartition: max min ratio","repartition: max imbalance")
      
         ("repartition: remap parts","repartition: remap parts")
      
         ("repartition: remap num values","repartition: remap num values")
      
         ("repartition: remap accept partition","repartition: remap accept partition")
      
         ("repartition: print partition distribution","repartition: print partition distribution")
      
         ("repartition: rebalance P and R","repartition: rebalance P and R")
      
         ("repartition: rebalance Nullspace","repartition: rebalance Nullspace")
      
         ("repartition: use subcommunicators","repartition: use subcommunicators")
      
         ("rap: relative diagonal floor","rap: relative diagonal floor")
      
         ("rap: fix zero diagonals","rap: fix zero diagonals")
      
         ("rap: fix zero diagonals threshold","rap: fix zero diagonals threshold")
      
         ("rap: shift","rap: shift")
      
         ("rap: shift diagonal M","rap: shift diagonal M")
      
         ("rap: shift low storage","rap: shift low storage")
      
         ("rap: shift array","rap: shift array")
      
         ("rap: cfl array","rap: cfl array")
      
         ("rap: algorithm","rap: algorithm")
      
         ("matrixmatrix: kernel params","matrixmatrix: kernel params")
      
         ("matrixmatrix: kernel params:MM_TAFC_OptimizationCoreCount ","matrixmatrix: kernel params:MM_TAFC_OptimizationCoreCount ")
      
         ("isMatrixMatrix_TransferAndFillComplete ","isMatrixMatrix_TransferAndFillComplete ")
      
         ("reuse: type","reuse: type")
      
         ("use external multigrid package","use external multigrid package")
      
         ("amgx:params","amgx:params")
      
         ("debug: graph level","debug: graph level")
      
         ("refmaxwell: mode","refmaxwell: mode")
      
         ("refmaxwell: disable addon","refmaxwell: disable addon")
      
         ("refmaxwell: 11list","refmaxwell: 11list")
      
         ("refmaxwell: 22list","refmaxwell: 22list")
      
         ("zero starting solution","refmaxwell: use as preconditioner")
      
         ("refmaxwell: dump matrices","refmaxwell: dump matrices")
      
         ("refmaxwell: subsolves on subcommunicators","refmaxwell: subsolves on subcommunicators")
      
         ("refmaxwell: ratio AH / A22 subcommunicators","refmaxwell: ratio AH / A22 subcommunicators")
      ;

}

