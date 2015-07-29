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
#include <Tpetra_Distributor.hpp>
#include <Tpetra_HashTable.hpp>
#include <Tpetra_Import.hpp>
#include <Tpetra_Import_Util.hpp>

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
    AMGXOperator(const Teuchos::RCP<Tpetra::CrsMatrix<SC,LO,GO,NO> > &inA, Teuchos::ParameterList &paramListIn) {
      RCP<const Teuchos::Comm<LO> > comm = inA->getCrsGraph()->getComm();

      // Initialize
      AMGX_SAFE_CALL(AMGX_initialize());
      AMGX_SAFE_CALL(AMGX_initialize_plugins());

      /*system*/
      //AMGX_SAFE_CALL(AMGX_register_print_callback(&print_callback));
      AMGX_SAFE_CALL(AMGX_install_signal_handler());
      Teuchos::ParameterList configs = paramListIn.sublist("amgx:params", true);
      if (configs.isParameter("json file")) {
        AMGX_SAFE_CALL(AMGX_config_create_from_file(&Config_, (const char *) &configs.get<std::string>("json file")[0]));
      } else {
        std::ostringstream oss;
        oss << "";
        ParameterList::ConstIterator itr;
        for (itr = configs.begin(); itr != configs.end(); ++itr) {
          const std::string&    name  = configs.name(itr);
          const ParameterEntry& entry = configs.entry(itr);
          oss << name << "=" << filterValueToString(entry) << ", ";
        }
        oss << "\0";
        std::string configString = oss.str();
        if (configString == "") {
          //print msg that using defaults
          //GetOStream(Warnings0) << "Warning: No configuration parameters specified, using default AMGX configuration parameters. \n";
        }
        AMGX_SAFE_CALL(AMGX_config_create(&Config_, (const char *) &configString[0]));
      }

      AMGX_Mode mode = AMGX_mode_dDDI;

      int nnz = -1;
      if (comm->getSize() > 1) {
        // FIXME: the below code handling MPI communication is wrong
        // Specifically, we completely missed the fact the matrix itself has to be reordered (packed) before passing it to AMGX.
        // AMGX_matrix_upload_all documentation states this:
        //     1. Select the rows (with global column indices) belonging to rank i and read them into memory
        //     2. Pack the columns:
        //        a) Reorder the column indices such that the columns corresponding to the global diagonal elements present
        //           on rank i come ﬁrst.
        //        b) Then, reorder the remaining column indices, in the order to which they belong (have connections) to
        //           rank’s i neighbors. Columns belonging to neighbor 0 first, neighbor 1 second, etc., leaving the natural
        //           ordering of columns within the same neighbor.
        //        c) Notice that we now have a rectangular matrix with reordered column indices that are numbered locally.
        //           To keep track of this reordering create a local to global map which indicates how the locally renumbered
        //           columns map into their global counterparts.
        //     3. Reorder the rows and columns (such that the rows with no connections to neighbors come ﬁrst):
        //        a) Find the reordering that moves rows with no connections to neighbors ﬁrst
        //        b) Apply it to rows and columns (belonging to the global diagonal elements present on rank i only, in other
        //           words, columns up to n).
        //     4. Pass the resulting packed local matrix to this routine.
        // This is non-trivial and expensive.
        RCP<const Tpetra::Import<LO,GO> > importer = inA->getCrsGraph()->getImporter();

        TEUCHOS_TEST_FOR_EXCEPTION(importer.is_null(), MueLu::Exceptions::RuntimeError, "The matrix A has no Import object.");

        Tpetra::Distributor distributor = importer->getDistributor();

        Array<int> sendRanks = distributor.getImagesTo();
        Array<int> recvRanks = distributor.getImagesFrom();

        std::sort(sendRanks.begin(), sendRanks.end());
        std::sort(recvRanks.begin(), recvRanks.end());

        bool match = true;
        if (sendRanks.size() != recvRanks.size()) {
          match = false;
        } else {
          for (int i = 0; i < sendRanks.size(); i++) {
            if (recvRanks[i] != sendRanks[i])
              match = false;
              break;
          }
        }
        TEUCHOS_TEST_FOR_EXCEPTION(!match, MueLu::Exceptions::RuntimeError, "AMGX requires that the processors that we send to and receive from are the same. "
                                   "This is not the case: we send to {" << sendRanks << "} and receive from {" << recvRanks << "}");

        int        num_neighbors = sendRanks.size();  // does not include the calling process
        const int* neighbors     = &sendRanks[0];

        // Later on, we'll have to organize the send and recv data by PIDs,
        // i.e, a vector V of vectors, where V[i] is PID i's vector of data.
        // Hence we need to be able to quickly look up  an array index
        // associated with each PID.
        Tpetra::Details::HashTable<int,int> hashTable(3*num_neighbors);
        for (int i = 0; i < num_neighbors; i++)
          hashTable.add(neighbors[i], i);

        // Construct send arrays
        ArrayView<const int> exportLIDs = importer->getExportLIDs();
        ArrayView<const int> exportPIDs = importer->getExportPIDs();

        std::vector<std::vector<int> > sendDatas (num_neighbors);
        std::vector<int>               send_sizes(num_neighbors, 0);
        for (int i = 0; i < exportPIDs.size(); i++) {
          int index = hashTable.get(exportPIDs[i]);
          sendDatas [index].push_back(exportLIDs[i]);
          send_sizes[index]++;
        }

        std::vector<const int*> send_maps(num_neighbors);
        for (int i = 0; i < num_neighbors; i++)
          send_maps[i] = &(sendDatas[i][0]);

        // Construct recv arrays
        Array<int> importPIDs;
        Tpetra::Import_Util::getPids(*importer, importPIDs, true/* make local -1 */);

        std::vector<std::vector<int> > recvDatas (num_neighbors);
        std::vector<int>               recv_sizes(num_neighbors, 0);
        for (int i = 0; i < importPIDs.size(); i++)
          if (importPIDs[i] != -1) {
            int index = hashTable.get(importPIDs[i]);
            recvDatas [index].push_back(i);
            recv_sizes[index]++;
        }
        std::vector<const int*> recv_maps(num_neighbors);
        for (int i = 0; i < num_neighbors; i++)
          recv_maps[i] = &(recvDatas[i][0]);

        RCP<const Teuchos::MpiComm<int> > tmpic = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm->duplicate());
        TEUCHOS_TEST_FOR_EXCEPTION(tmpic.is_null(), Exceptions::RuntimeError, "Communicator is not MpiComm");
        RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > rawMpiComm = tmpic->getRawMpiComm();
        MPI_Comm mpiComm = *rawMpiComm;

        int numGPUDevices;
        cudaGetDeviceCount(&numGPUDevices);
        int device[] = {(comm->getRank() % numGPUDevices)};

        // TODO: we probably need to add "exception_handling=1" to the parameter list
        // to switch on internal error handling (with no need for AMGX_SAFE_CALL)
        AMGX_config_add_parameters(&Config_, "communicator=MPI");
        AMGX_resources_create(&Resources_, Config_, &mpiComm, 1/* number of GPU devices utilized by this rank */, device);

        AMGX_solver_create(&Solver_, Resources_, mode,  Config_);
        AMGX_matrix_create(&A_,      Resources_, mode);
        AMGX_vector_create(&X_,      Resources_, mode);
        AMGX_vector_create(&Y_,      Resources_, mode);

        AMGX_matrix_comm_from_maps_one_ring(A_, 1, num_neighbors, neighbors, &send_sizes[0], &send_maps[0], &recv_sizes[0], &recv_maps[0]);

        AMGX_vector_bind(X_, A_);
        AMGX_vector_bind(Y_, A_);

        N_  = inA->getNodeNumRows();
        nnz = inA->getNodeNumEntries();

      } else {
        AMGX_resources_create_simple(&Resources_, Config_);

        AMGX_solver_create(&Solver_, Resources_, mode,  Config_);
        AMGX_matrix_create(&A_,      Resources_, mode);
        AMGX_vector_create(&X_,      Resources_, mode);
        AMGX_vector_create(&Y_,      Resources_, mode);

        N_  = inA->getGlobalNumRows();
        nnz = inA->getGlobalNumEntries();
      }

      domainMap_ = inA->getDomainMap();
      rangeMap_  = inA->getRangeMap();

      ArrayRCP<const size_t> row_ptr_t;
      ArrayRCP<const int>    col_ind;
      ArrayRCP<const double> data;
      inA->getAllValues(row_ptr_t, col_ind, data);

      Teuchos::ArrayRCP<int> row_ptr(row_ptr_t.size(), 0);
      for (int i = 0; i < row_ptr_t.size(); i++)
        row_ptr[i] = Teuchos::as<int, size_t>(row_ptr_t[i]);

      // TODO We also need to pin memory here through AMGX_pin_memory
      AMGX_matrix_upload_all(A_, N_, nnz, 1, 1, &row_ptr[0], &col_ind[0], &data[0], NULL);
      AMGX_solver_setup(Solver_, A_);
    }

    //! Destructor.
    virtual ~AMGXOperator() { }

    //! Returns the Tpetra::Map object associated with the domain of this operator.
    Teuchos::RCP<const Tpetra::Map<LO,GO,NO> > getDomainMap() const;

    //! Returns the Tpetra::Map object associated with the range of this operator.
    Teuchos::RCP<const Tpetra::Map<LO,GO,NO> > getRangeMap() const;

    //! Returns in X the solution to the linear system AX=Y.
    /*!
       \param[out] X - Tpetra::MultiVector of dimension NumVectors containing the RHS of the linear system
       \param[in]  Y - Tpetra::MultiVector of dimension NumVectors containing the solution to the linear system
       */
    void apply(const Tpetra::MultiVector<SC,LO,GO,NO>& X, Tpetra::MultiVector<SC,LO,GO,NO>& Y,
               Teuchos::ETransp mode = Teuchos::NO_TRANS,
               SC alpha = Teuchos::ScalarTraits<SC>::one(), SC beta = Teuchos::ScalarTraits<SC>::one()) const;

    //! Indicates whether this operator supports applying the adjoint operator.
    bool hasTransposeApply() const;

    RCP<MueLu::Hierarchy<SC,LO,GO,NO> > GetHierarchy() const {
      throw Exceptions::RuntimeError("AMGXOperator does not hold a MueLu::Hierarchy object \n");
    }

    std::string filterValueToString(const Teuchos::ParameterEntry& entry ) {
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
    int                    N_;

    Teuchos::RCP<const Tpetra::Map<LO,GO,NO> > domainMap_;
    Teuchos::RCP<const Tpetra::Map<LO,GO,NO> > rangeMap_;
  };

} // namespace

#endif //HAVE_MUELU_EXPERIMENTAL && HAVE_MUELU_EXPERIMENTAL
#endif // MUELU_AMGXOPERATOR_DECL_HPP
