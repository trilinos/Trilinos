#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_CUSPARSEOps.hpp"

Teuchos::RCP<cusparseHandle_t> Kokkos::CUSPARSEdetails::Session::session_handle_ = Teuchos::null;

void Kokkos::CUSPARSEdetails::Session::init() 
{
  if (session_handle_ == null) {
    cusparseHandle_t *hndl = new cusparseHandle_t;
    cusparseStatus_t status = cusparseCreate(hndl);
    TEUCHOS_TEST_FOR_EXCEPTION(status == CUSPARSE_STATUS_NOT_INITIALIZED, std::runtime_error, 
        "Kokkos::CUSPARSEdetails::init(): the CUDA Runtime initialization failed.")
    TEUCHOS_TEST_FOR_EXCEPTION(status == CUSPARSE_STATUS_ALLOC_FAILED, std::runtime_error, 
        "Kokkos::CUSPARSEdetails::init(): the resources could not be allocated.")
    TEUCHOS_TEST_FOR_EXCEPTION(status == CUSPARSE_STATUS_ARCH_MISMATCH, std::runtime_error, 
        "Kokkos::CUSPARSEdetails::init(): the device compute capability (CC) is less than 1.1.")
    TEUCHOS_TEST_FOR_EXCEPTION(status != CUSPARSE_STATUS_SUCCESS, std::runtime_error, 
        "Kokkos::CUSPARSEdetails::init(): unspecified error.")
    session_handle_ = Teuchos::rcpWithDealloc(hndl,CUSPARSESessionDestroyer(),true);
  }
}

Teuchos::RCP<const cusparseHandle_t> Kokkos::CUSPARSEdetails::Session::getHandle() {
  return session_handle_;
}

Teuchos::RCP<cusparseMatDescr_t> Kokkos::CUSPARSEdetails::createMatDescr() 
{
  cusparseMatDescr_t *desc = new cusparseMatDescr_t;
  cusparseStatus_t status = cusparseCreateMatDescr(desc);
  TEUCHOS_TEST_FOR_EXCEPTION(status == CUSPARSE_STATUS_ALLOC_FAILED, std::runtime_error, 
      "Kokkos::CUSPARSEdetails::createMatDescr(): the resources could not be allocated.")
  TEUCHOS_TEST_FOR_EXCEPTION(status != CUSPARSE_STATUS_SUCCESS, std::runtime_error, 
      "Kokkos::CUSPARSEdetails::createMatDescr(): unspecified error.")
  return Teuchos::rcpWithDealloc(desc,CUSPARSEMatDescDestroyer(),true);
}

Teuchos::RCP<cusparseSolveAnalysisInfo_t> Kokkos::CUSPARSEdetails::createSolveAnalysisInfo()
{
  cusparseSolveAnalysisInfo_t *info = new cusparseSolveAnalysisInfo_t;
  cusparseStatus_t status = cusparseCreateSolveAnalysisInfo(info);
  TEUCHOS_TEST_FOR_EXCEPTION(status == CUSPARSE_STATUS_ALLOC_FAILED, std::runtime_error, 
      "Kokkos::CUSPARSEdetails::createSolveAnalysisInfo(): the resources could not be allocated.")
  TEUCHOS_TEST_FOR_EXCEPTION(status != CUSPARSE_STATUS_SUCCESS, std::runtime_error, 
      "Kokkos::CUSPARSEdetails::createSolveAnalysisInfo(): unspecified error.")
  return Teuchos::rcpWithDealloc(info,CUSPARSESolveAnalysisDestroyer(),true);
}

//////// deallocator for cusparse sessions
Kokkos::CUSPARSEdetails::CUSPARSESessionDestroyer::CUSPARSESessionDestroyer() {}
void Kokkos::CUSPARSEdetails::CUSPARSESessionDestroyer::free(cusparseHandle_t *ptr) 
{
  cusparseStatus_t status = cusparseDestroy( *ptr );
  delete ptr;
  TEUCHOS_TEST_FOR_EXCEPTION(status != CUSPARSE_STATUS_SUCCESS, std::runtime_error,
      "Kokkos::CUSPARSEdetails::CUSPARSESessionDestroyer::free(): library was never initialized; we should not have been called.")
}

//////// deallocator for cusparse matrix description structs
Kokkos::CUSPARSEdetails::CUSPARSEMatDescDestroyer::CUSPARSEMatDescDestroyer() {}
void Kokkos::CUSPARSEdetails::CUSPARSEMatDescDestroyer::free(cusparseMatDescr_t *ptr) 
{
  cusparseStatus_t status = cusparseDestroyMatDescr( *ptr );
  delete ptr;
  TEUCHOS_TEST_FOR_EXCEPTION(status != CUSPARSE_STATUS_SUCCESS, std::runtime_error,
      "Kokkos::CUSPARSEdetails::CUSPARSEMatDescDestroyer::free(): unknown error.")
}

//////// deallocator for cusparse solve analysis structs
Kokkos::CUSPARSEdetails::CUSPARSESolveAnalysisDestroyer::CUSPARSESolveAnalysisDestroyer() {}
void Kokkos::CUSPARSEdetails::CUSPARSESolveAnalysisDestroyer::free(cusparseSolveAnalysisInfo_t *ptr) 
{
  cusparseStatus_t status = cusparseDestroySolveAnalysisInfo( *ptr );
  delete ptr;
  TEUCHOS_TEST_FOR_EXCEPTION(status != CUSPARSE_STATUS_SUCCESS, std::runtime_error,
      "Kokkos::CUSPARSEdetails::CUSPARSESolveAnalysisDestroyer::free(): unknown error.")
}

