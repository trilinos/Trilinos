#ifndef MUELU_CREATE_AMGX_PRECONDITIONER_HPP
#define MUELU_CREATE_AMGX_PRECONDITIONER_HPP

#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Tpetra_CrsMatrix_decl.hpp>
#include "cuda_runtime.h"
#include <amgx_c.h>

//! @file MueLu_CreateAMGXPreconditioner.hpp

namespace MueLu {

  /*! \fn CreateTpetraPreconditioner
    @brief Helper function to create a MueLu preconditioner that can be used by Tpetra.

    Given a Tpetra matrix, this function returns a constructed MueLu preconditioner.

    @param[in] inA Matrix
    @param[in] paramList Parameter list
    @param[in] inCoords (optional) Coordinates.  The first vector is x, the second (if necessary) y, the third (if necessary) z.
    @param[in] inNullspace (optional) Near nullspace of the matrix.
    */
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
 // Teuchos::RCP<MueLu::TpetraOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  Teuchos::RCP<MueLu::AMGXOperator>
  CreateAMGXPreconditioner(const Teuchos::RCP<Tpetra::CrsMatrix  <Scalar, LocalOrdinal, GlobalOrdinal, Node> >& inA,
                             Teuchos::ParameterList& paramListIn,
                             const Teuchos::RCP<Tpetra::MultiVector<double, LocalOrdinal, GlobalOrdinal, Node> >& inCoords    = Teuchos::null,
                             const Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& inNullspace = Teuchos::null)
  {

    using   Teuchos::ParameterList;

    bool hasParamList=paramListIn.numParams();

	    
    AMGX_Mode mode;
    AMGX_config_handle cfg;
    AMGX_resources_handle rsrc;
    AMGX_matrix_handle A;
    AMGX_solver_handle solver;
    //status handling
    AMGX_SOLVE_STATUS status;

    /* init */
    AMGX_SAFE_CALL(AMGX_initialize());
    AMGX_SAFE_CALL(AMGX_initialize_plugins());

    mode = AMGX_mode_dDDI;

    /*transform matrix into crs arrays*/
    /*need to transform size_t  into int*/ 
    Teuchos::ArrayRCP<size_t> row_ptr;
    Teuchos::ArrayRCP<int> col_ind;
    Teuchos::ArrayRCP<double> data;
 
    int N = inA->getGlobalNumRows();
 
    int nnz = inA->getGlobalNumEntries();

    
    inA->getAllValues(row_ptr, col_ind, data);


     

    /*need to use paramlist to set up config*/
    if(hasParamList && paramListIn.isParameter("config string")){	
       AMGX_SAFE_CALL(AMGX_config_create(&cfg, (const string) paramListIn.get("config string"));

    }
    else if(hasParamList && paramListIn.isParameter("config file"){
	
       AMGX_SAFE_CALL(AMGX_config_create_from_file(&cfg, (const string) paramListIn.get("config file"));

    }
    else{
       AMGX_SAFE_CALL(AMGX_config_create(&cfg, "exception_handling=1");
    }
    


    /*create resources -- for one device on one node*/
    AMGX_resources_create_simple(&rsrc, cfg);

    AMGX_matrix_create(&A, rsrc, mode);
    AMGX_solver_create(&solver, rsrc, mode, cfg);

    AMGX_matrix_upload_all(A, N, nnz, 1, 1, &row_ptr[0], &col_ind[0], &data[0], NULL);
    
    AMGX_solver_setup(solver, A);


    
    
    return rcp(new MueLu::AMGXOperator(solver, rsrc, cfg, A, N);

  }

  /*! \fn CreateTpetraPreconditioner
    @brief Helper function to create a MueLu preconditioner that can be used by Tpetra.

    Given a Tpetra matrix, this function returns a constructed MueLu preconditioner.

    @param[in] inA Matrix
    @param[in] inCoords (optional) Coordinates.  The first vector is x, the second (if necessary) y, the third (if necessary) z.
    @param[in] inNullspace (optional) Near nullspace of the matrix.
    */
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
 // Teuchos::RCP<MueLu::TpetraOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > 
  Teuchos::RCP<MueLu::AMGXOperator
  CreateAMGXPreconditioner(const Teuchos::RCP<Tpetra::CrsMatrix  <Scalar, LocalOrdinal, GlobalOrdinal, Node> >& inA,
                             const Teuchos::RCP<Tpetra::MultiVector<double, LocalOrdinal, GlobalOrdinal, Node> >& inCoords    = Teuchos::null,
                             const Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& inNullspace = Teuchos::null) {
    /*create void paramList resulting in default config settings*/
    Teuchos::ParameterList paramList;
    return CreateAMGXPreconditioner(inA, paramList, inCoords, inNullspace);
  }

  /*! \fn CreateTpetraPreconditioner
    @brief Helper function to create a MueLu preconditioner that can be used by Tpetra.

    Given a Tpetra matrix, this function returns a constructed MueLu preconditioner.

    @param[in] inA Matrix
    @param[in] xmlFileName XML file containing MueLu options
    @param[in] inCoords (optional) Coordinates.  The first vector is x, the second (if necessary) y, the third (if necessary) z.
    @param[in] inNullspace (optional) Near nullspace of the matrix.
    */




  /*Is this method necessary?*/
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ReuseAMGXPreconditioner(const Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& inA,
                                 MueLu::AMGXOperator& Op) {
    typedef Scalar          SC;
    typedef LocalOrdinal    LO;
    typedef GlobalOrdinal   GO;
    typedef Node            NO;
    
    /*need to get arrays of inA matrix and reset solver*/
    
    Teuchos::ArrayRCP<size_t> row_ptr;
    Teuchos::ArrayRCP<int> col_ind;
    Teuchos::ArrayRCP<double> data;
 
    int N = inA->getGlobalNumRows();
 
    int nnz = inA->getGlobalNumEntries();

    
    inA->getAllValues(row_ptr, col_ind, data);
  
    AMGX_Mode mode = dDDI;
  
    AMGX_matrix_handle A;		
    AMGX_matrix_create(&A, Op.get_rsrc(), mode);

    AMGX_matrix_upload_all(A, N, nnz, 1, 1, &row_ptr[0], &col_ind[0], &data[0], NULL);
   
    /*need to get solver object (and other objects from AMGXOperator) - need to write methods to accomplish this*/ 
    AMGX_solver_setup(solver, A);

    Op = rcp(new MueLu::AMGXPreconditioner(solver, Op.get_rsrc(), Op.get_cfg(), A, N));
  }

} //namespace

#endif //ifndef MUELU_CREATE_AMGX_PRECONDITIONER_HPP
