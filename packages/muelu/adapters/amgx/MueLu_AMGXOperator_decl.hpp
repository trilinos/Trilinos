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
#ifndef MUELU_AMGXOPERATOR_DECL_HPP
#define MUELU_AMGXOPERATOR_DECL_HPP

#if defined (HAVE_MUELU_EXPERIMENTAL) and defined (HAVE_MUELU_AMGX)
#include <Teuchos_ParameterList.hpp>

#include <Tpetra_Operator.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_MultiVector.hpp>

#include "MueLu_Exceptions.hpp"
#include "MueLu_TpetraOperator.hpp"
#include "MueLu_VerboseObject.hpp"

#include <cuda_runtime.h>
#include <amgx_c.h>

namespace MueLu {


  /*! @class TemplatedAMGXOperator
      This templated version of the class throws errors in all methods as AmgX is not implemented for datatypes where scalar!=double/float and ordinal !=int
 */
  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node>
  class AMGXOperator : public TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
  private:
    typedef Scalar          SC;
    typedef LocalOrdinal    LO;
    typedef GlobalOrdinal   GO;
    typedef Node            NO;

  public:

    //! @name Constructor/Destructor
    //@{

    //! Constructor
    AMGXOperator(const Teuchos::RCP<Tpetra::CrsMatrix<SC,LO,GO,NO> > &InA, Teuchos::ParameterList &paramListIn) { }

    //! Destructor.
    virtual ~AMGXOperator() { }

    //@}

    //! Returns the Tpetra::Map object associated with the domain of this operator.
    Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > getDomainMap() const{
      throw Exceptions::RuntimeError("Cannot use AMGXOperator with scalar != double and/or global ordinal != int \n");
    }

    //! Returns the Tpetra::Map object associated with the range of this operator.
    Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > getRangeMap() const{
      throw Exceptions::RuntimeError("Cannot use AMGXOperator with scalar != double and/or global ordinal != int \n");
    }

    //! Returns a solution for the linear system AX=Y in the  Tpetra::MultiVector X.
    /*!
      \param[in]  X - Tpetra::MultiVector of dimension NumVectors that contains the solution to the linear system.
      \param[out] Y -Tpetra::MultiVector of dimension NumVectors containing the RHS of the linear system.
    */
    void apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
                                         Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
                                         Teuchos::ETransp mode = Teuchos::NO_TRANS,
                                         Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
                                         Scalar beta  = Teuchos::ScalarTraits<Scalar>::one()) const {
      throw Exceptions::RuntimeError("Cannot use AMGXOperator with scalar != double and/or global ordinal != int \n");
    }

    //! Indicates whether this operator supports applying the adjoint operator
    bool hasTransposeApply() const{
      throw Exceptions::RuntimeError("Cannot use AMGXOperator with scalar != double and/or global ordinal != int \n");
    }

    RCP<MueLu::Hierarchy<SC,LO,GO,NO> > GetHierarchy() const {
        throw Exceptions::RuntimeError("AMGXOperator does not hold a MueLu::Hierarchy object \n");
    }

  private:
  };

  /*! @class AMGXOperator
      Creates and AmgX Solver object with a Tpetra Matrix. Partial specialization of the template for data types supported by AmgX.
  */
  template<class Node>
  class AMGXOperator<double, int, int, Node> : public TpetraOperator<double, int, int, Node> {
  private:
    typedef double  SC;
    typedef int     LO;
    typedef int     GO;
    typedef Node    NO;

  public:

    //! @name Constructor/Destructor
    //@{
    AMGXOperator(const Teuchos::RCP<Tpetra::CrsMatrix<SC,LO,GO,NO> > &inA, Teuchos::ParameterList &paramListIn){
      AMGX_Mode mode;
      /* init */
      AMGX_SAFE_CALL(AMGX_initialize());
      AMGX_SAFE_CALL(AMGX_initialize_plugins());
      /*system*/
      //AMGX_SAFE_CALL(AMGX_register_print_callback(&print_callback));
      AMGX_SAFE_CALL(AMGX_install_signal_handler());

      mode = AMGX_mode_dDDI;
      Teuchos::ArrayRCP<const size_t> row_ptr_t;
      Teuchos::ArrayRCP<const int> col_ind;
      Teuchos::ArrayRCP<const double> data;

      domainMap_ = inA->getDomainMap();
      rangeMap_ = inA->getRangeMap();

      N = inA->getGlobalNumRows();
      int nnz = inA->getGlobalNumEntries();

      inA->getAllValues(row_ptr_t, col_ind, data);


      Teuchos::ArrayRCP<int> row_ptr(row_ptr_t.size(), 0);
      for(int i = 0; i < row_ptr_t.size(); i++){
        row_ptr[i] = Teuchos::as<int, size_t> (row_ptr_t[i]);
      }
      Teuchos::ParameterList configs = paramListIn.sublist("amgx:params", true);
      if (configs.isParameter("json file")) {
        AMGX_SAFE_CALL(AMGX_config_create_from_file(&Config_, (const char *) &configs.get<std::string>("json file")[0]));
      } else {
        std::ostringstream oss;
        oss << "";
        ParameterList::ConstIterator itr;
        for( itr = configs.begin(); itr != configs.end(); ++itr){
          const std::string & paramName = configs.name(itr);
          const ParameterEntry &value = configs.entry(itr);
          oss << paramName << "=" << filterValueToString(value) << ", ";
        }
        oss << "\0";
        std::string configString = oss.str();
        if (configString == "") {
          //print msg that using defaults
    //      GetOStream(Warnings0) << "Warning: No configuration parameters specified, using default AMGX configuration parameters. \n";
        }
        AMGX_SAFE_CALL(AMGX_config_create(&Config_, (const char *) &configString[0]));
      }
      AMGX_resources_create_simple(&Resources_, Config_);
      AMGX_matrix_create(&A_, Resources_, mode);
      AMGX_solver_create(&Solver_, Resources_, mode, Config_);

      AMGX_matrix_upload_all(A_, N, nnz, 1, 1, &row_ptr[0], &col_ind[0], &data[0], NULL);
      AMGX_solver_setup(Solver_, A_);
      AMGX_vector_create(&X_, Resources_, mode);
      AMGX_vector_create(&Y_, Resources_, mode);
    }

    //! Destructor.
    virtual ~AMGXOperator() { }

    //! Returns the Tpetra::Map object associated with the domain of this operator.
    Teuchos::RCP<const Tpetra::Map<LO,GO,NO> > getDomainMap() const;

    //! Returns the Tpetra::Map object associated with the range of this operator.
    Teuchos::RCP<const Tpetra::Map<LO,GO,NO> > getRangeMap() const;

    //! Returns in X the solution to the linear system AX=Y.
    /*!
       \param[in]  X - Tpetra::MultiVector of dimension NumVectors containing the solution to the linear system
       \param[out] Y -Tpetra::MultiVector of dimension NumVectors containing the RHS of the linear system.                 */
    void apply(const Tpetra::MultiVector<SC,LO,GO,NO>& X, Tpetra::MultiVector<SC,LO,GO,NO>& Y,
               Teuchos::ETransp mode = Teuchos::NO_TRANS,
               SC alpha = Teuchos::ScalarTraits<SC>::one(),                                                                           SC beta  = Teuchos::ScalarTraits<SC>::one()) const;

    //! Indicates whether this operator supports applying the adjoint operator.
    bool hasTransposeApply() const;

    RCP<MueLu::Hierarchy<SC,LO,GO,NO> > GetHierarchy() const {
        throw Exceptions::RuntimeError("AMGXOperator does not hold a MueLu::Hierarchy object \n");
    }

    std::string filterValueToString(const Teuchos::ParameterEntry& entry )
    {
      return ( entry.isList() ? std::string("...") : toString(entry.getAny()) );
    }

    int sizeA() {
      int sizeX, sizeY, n;
      AMGX_matrix_get_size(A_, &n, &sizeX, &sizeY);
      return n;
    }

    int iters() {
      int it;
      AMGX_solver_get_iterations_number(Solver_, &it);
      return it;
    }
 
    AMGX_SOLVE_STATUS getStatus() {
      AMGX_SOLVE_STATUS status;
      AMGX_solver_get_status(Solver_, &status);
      return status;
    }
      

  private:
    AMGX_solver_handle     Solver_;
    AMGX_resources_handle  Resources_;
    AMGX_config_handle     Config_;
    AMGX_matrix_handle     A_;
    AMGX_vector_handle     X_;
    AMGX_vector_handle     Y_;
    int N;
    Teuchos::RCP<const Tpetra::Map<LO,GO,NO> > domainMap_;
    Teuchos::RCP<const Tpetra::Map<LO,GO,NO> > rangeMap_;
  };

} // namespace

#endif //HAVE_MUELU_EXPERIMENTAL && HAVE_MUELU_EXPERIMENTAL
#endif // MUELU_AMGXOPERATOR_DECL_HPP
