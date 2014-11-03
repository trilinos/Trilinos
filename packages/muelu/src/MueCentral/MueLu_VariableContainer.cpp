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
#include "Kokkos_DefaultNode.hpp"
#include "Tpetra_ETIHelperMacros.h"

#include "MueLu_VariableContainer.hpp"

namespace MueLu {

#define MUELU_INST1_SC_LO_GO_NO(SC,LO,GO,NO) \
  template<> \
  const Teuchos::RCP<Xpetra::Operator<SC,LO,GO,NO> >& VariableContainer::GetData<Teuchos::RCP<Xpetra::Operator<SC,LO,GO,NO> > >() const { \
    typedef Xpetra::Operator<SC,LO,GO,NO> Operator; \
    typedef Xpetra::Matrix  <SC,LO,GO,NO> Matrix; \
    \
    typedef const Teuchos::RCP<Operator> TO; \
    typedef const Teuchos::RCP<Matrix>   TM; \
    \
    const std::string typeTOName = Teuchos::TypeNameTraits<TO>::name(); \
    const std::string typeTMName = Teuchos::TypeNameTraits<TM>::name(); \
    TEUCHOS_TEST_FOR_EXCEPTION(data_ == NULL, Teuchos::bad_any_cast, \
                               "Error, cast to type Data<" << typeTOName << "> failed since the content is NULL"); \
    if (data_->type() == typeid(TO)) { \
      Data<TO>* data = dynamic_cast<Data<TO>*>(data_); \
      TEUCHOS_TEST_FOR_EXCEPTION(!data, std::logic_error, \
                                 "Error, cast to type Data<" << typeTOName << "> failed but should not have and the actual underlying type is " \
                                 "\'" << data_->typeName() << "! The problem might be related to incompatible RTTI systems in static and shared libraries!"); \
      return data->data_; \
    } \
    \
    TEUCHOS_TEST_FOR_EXCEPTION(data_->type() != typeid(TM), Teuchos::bad_any_cast, \
                               "Error, cast to type Data<" << typeTMName << "> failed since the actual underlying type is " \
                               "\'" << data_->typeName() << "!"); \
    Data<TM>* data = dynamic_cast<Data<TM>*>(data_); \
    TEUCHOS_TEST_FOR_EXCEPTION(!data, std::logic_error, \
                               "Error, cast to type Data<" << typeTMName << "> failed but should not have and the actual underlying type is " \
                               "\'" << data_->typeName() << "! The problem might be related to incompatible RTTI systems in static and shared libraries!"); \
    if (datah_ == NULL) \
      datah_ = new Data<TO>(rcp_dynamic_cast<Operator>(data->data_)); \
    Data<TO>* datah = dynamic_cast<Data<TO>*>(datah_); \
    TEUCHOS_TEST_FOR_EXCEPTION(!datah, std::logic_error, \
                               "Error, cast to type Data<" << typeTOName << "> failed but should not have and the actual underlying type is " \
                               "\'" << datah_->typeName() << "! The problem might be related to incompatible RTTI systems in static and shared libraries!"); \
    return datah->data_; \
  }

#define MUELU_INST2_SC_LO_GO_NO(SC,LO,GO,NO) \
  template<> \
  Teuchos::RCP<Xpetra::Operator<SC,LO,GO,NO> >& VariableContainer::GetData<Teuchos::RCP<Xpetra::Operator<SC,LO,GO,NO> > >() { \
    typedef Xpetra::Operator<SC,LO,GO,NO> Operator; \
    typedef Xpetra::Matrix  <SC,LO,GO,NO> Matrix; \
    \
    typedef Teuchos::RCP<Operator> TO; \
    typedef Teuchos::RCP<Matrix>   TM; \
    \
    const std::string typeTOName = Teuchos::TypeNameTraits<TO>::name(); \
    const std::string typeTMName = Teuchos::TypeNameTraits<TM>::name(); \
    TEUCHOS_TEST_FOR_EXCEPTION(data_ == NULL, Teuchos::bad_any_cast, \
                               "Error, cast to type Data<" << typeTOName << "> failed since the content is NULL"); \
    if (data_->type() == typeid(TO)) { \
      Data<TO>* data = dynamic_cast<Data<TO>*>(data_); \
      TEUCHOS_TEST_FOR_EXCEPTION(!data, std::logic_error, \
                                 "Error, cast to type Data<" << typeTOName << "> failed but should not have and the actual underlying type is " \
                                 "\'" << data_->typeName() << "! The problem might be related to incompatible RTTI systems in static and shared libraries!"); \
      return data->data_; \
    } \
    \
    TEUCHOS_TEST_FOR_EXCEPTION(data_->type() != typeid(TM), Teuchos::bad_any_cast, \
                               "Error, cast to type Data<" << typeTMName << "> failed since the actual underlying type is " \
                               "\'" << data_->typeName() << "!"); \
    Data<TM>* data = dynamic_cast<Data<TM>*>(data_); \
    TEUCHOS_TEST_FOR_EXCEPTION(!data, std::logic_error, \
                               "Error, cast to type Data<" << typeTMName << "> failed but should not have and the actual underlying type is " \
                               "\'" << data_->typeName() << "! The problem might be related to incompatible RTTI systems in static and shared libraries!"); \
    if (datah_ == NULL) \
      datah_ = new Data<TO>(rcp_dynamic_cast<Operator>(data->data_)); \
    Data<TO>* datah = dynamic_cast<Data<TO>*>(datah_); \
    TEUCHOS_TEST_FOR_EXCEPTION(!datah, std::logic_error, \
                               "Error, cast to type Data<" << typeTOName << "> failed but should not have and the actual underlying type is " \
                               "\'" << datah_->typeName() << "! The problem might be related to incompatible RTTI systems in static and shared libraries!"); \
    return datah->data_; \
  }

#ifdef HAVE_MUELU_INST_DOUBLE_INT_INT
  MUELU_INST1_SC_LO_GO_NO(double,int,int,KokkosClassic::DefaultNode::DefaultNodeType);
  MUELU_INST2_SC_LO_GO_NO(double,int,int,KokkosClassic::DefaultNode::DefaultNodeType);
#endif

#ifdef HAVE_MUELU_INST_DOUBLE_INT_LONGLONGINT
# ifdef HAVE_TEUCHOS_LONG_LONG_INT
  MUELU_INST1_SC_LO_GO_NO(double,int,long long int,KokkosClassic::DefaultNode::DefaultNodeType);
  MUELU_INST2_SC_LO_GO_NO(double,int,long long int,KokkosClassic::DefaultNode::DefaultNodeType);
# else
# warning To compile MueLu with 'long long int' support, please turn on Teuchos_ENABLE_LONG_LONG_INT
# endif
#endif

#ifdef HAVE_MUELU_INST_COMPLEX_INT_INT
# ifdef HAVE_TEUCHOS_COMPLEX
#include <complex>
  MUELU_INST1_SC_LO_GO_NO(std::complex<double>,int,long long int,KokkosClassic::DefaultNode::DefaultNodeType);
  MUELU_INST2_SC_LO_GO_NO(std::complex<double>,int,long long int,KokkosClassic::DefaultNode::DefaultNodeType);
# else
# warning To compile MueLu with 'complex' support, please turn on Teuchos_ENABLE_COMPLEX
# endif
#endif

#if defined(HAVE_KOKKOSCLASSIC_KOKKOSCOMPAT) && defined(KOKKOS_HAVE_PTHREAD) && defined(HAVE_MUELU_INST_DOUBLE_INT_INT) && !defined(HAVE_KOKKOSCLASSIC_DEFAULTNODE_THREADSWRAPPERNODE)
  MUELU_INST1_SC_LO_GO_NO(double,int,int,Kokkos_Compat_KokkosThreadsWrapperNode);
  MUELU_INST2_SC_LO_GO_NO(double,int,int,Kokkos_Compat_KokkosThreadsWrapperNode);
#endif

#if defined(HAVE_KOKKOSCLASSIC_THREADPOOL) && defined(HAVE_MUELU_INST_DOUBLE_INT_INT) && !defined(HAVE_KOKKOSCLASSIC_DEFAULTNODE_TPINODE)
  MUELU_INST1_SC_LO_GO_NO(double,int,int,KokkosClassic_TPINode);
  MUELU_INST2_SC_LO_GO_NO(double,int,int,KokkosClassic_TPINode);
#endif

#if defined(HAVE_KOKKOSCLASSIC_KOKKOSCOMPAT) && defined(KOKKOS_HAVE_OPENMP) && defined(HAVE_MUELU_INST_DOUBLE_INT_INT) && !defined(HAVE_KOKKOSCLASSIC_DEFAULTNODE_OPENMPWRAPPERNODE)
  MUELU_INST1_SC_LO_GO_NO(double,int,int,Kokkos_Compat_KokkosOpenMPWrapperNode);
  MUELU_INST2_SC_LO_GO_NO(double,int,int,Kokkos_Compat_KokkosOpenMPWrapperNode);
#endif

#undef MUELU_INST2_SC_LO_GO_NO
#undef MUELU_INST1_SC_LO_GO_NO

}
