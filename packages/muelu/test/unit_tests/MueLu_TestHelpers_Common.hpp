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
#ifndef MUELU_TEST_HELPERS_COMMON_HPP
#define MUELU_TEST_HELPERS_COMMON_HPP

#define MUELU_DEBUGGER_MACRO \
  GO tmpIndex; out << "debug: scalar = " << typeid(tmpIndex).name() << std::endl; \
  NO tmpNode;  out << "debug: node   = " << typeid(tmpNode).name() << std::endl; \
  if (TestHelpers::Parameters::getLib() == Xpetra::UseEpetra) out << "debug: linAlgebra = Epetra" << std::endl; \
  else                                                        out << "debug: linAlgebra = Tpetra" << std::endl;

//#define MUELU_DEBUGGER_MACRO

// If linAlgebra==Epetra, this macro will return early from the test
// if SC!=double, GO!={int}, or NO!=Serial.
#define MUELU_TESTING_LIMIT_EPETRA_SCOPE(SC, GO, NO) \
  MUELU_DEBUGGER_MACRO \
  if (TestHelpers::Parameters::getLib() == Xpetra::UseEpetra) { \
    NO nodeCheck; \
    std::string nodeName = typeid(nodeCheck).name(); \
    if (nodeName.find("Serial") == std::string::npos) { \
      out << "Skipping Epetra for non-Serial nodes" << std::endl; \
      return; \
    } \
    GO index; \
    int epetraIntegerCheck; \
    std::string goName = typeid(index).name(); \
    if (goName != typeid(epetraIntegerCheck).name()) { \
       out << "Skipping Epetra for GO other than \"int\" and \"long long\"" << std::endl; \
       return; \
    } \
    SC epetraScalarCheck; \
    double doubleScalar; \
    if (typeid(doubleScalar).name() != typeid(epetraScalarCheck).name()) { \
       out << "Skipping Epetra for SC other than \"double\"" << std::endl; \
       return; \
    } \
  }

// If linAlgebra==Tpetra, but the test also requires Epetra, this macro will cause the test
// to return early if SC!=double, GO!={int}, or NO!=Serial.
#define MUELU_TESTING_LIMIT_EPETRA_SCOPE_TPETRA_IS_DEFAULT(SC, GO, NO) \
    Node nodeCheck; \
    std::string nodeName = typeid(nodeCheck).name();  \
    if (nodeName.find("Serial") == std::string::npos) { \
      out << "Skipping Epetra for non-Serial nodes" << std::endl; \
      return; \
    } \
    if (Teuchos::OrdinalTraits<GlobalOrdinal>::name() != std::string("int")) { \
       out << "Skipping Epetra for GO other than \"int\" and \"long long\"" << std::endl; \
       return; \
    } \
    if (Teuchos::ScalarTraits<Scalar>::name() != std::string("double")) { \
       out << "Skipping Epetra for SC other than \"double\"" << std::endl; \
       return; \
    }

//Macro to set MueLu's internal oh-so FancyOStream to be the same as the one used by Teuchos' unit testing framework.
//This prevents MueLu's output from intermingling with with the unit test pass/fail summary lines.
#define MUELU_TESTING_SET_OSTREAM \
   MueLu::VerboseObject::SetDefaultOStream(Teuchos::fancyOStream(out.getOStream()));


#endif // ifndef MUELU_TEST_HELPERS_COMMON_HPP
