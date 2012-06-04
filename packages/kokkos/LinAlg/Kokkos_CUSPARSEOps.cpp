#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_CUSPARSEOps.hpp"

Kokkos::CUSPARSEdetails::CUSPARSEDestroyer::CUSPARSEDestroyer() {}

void Kokkos::CUSPARSEdetails::CUSPARSEDestroyer::free(void *ptr) 
{
  if (ptr) {
    cusparseStatus_t status = cusparseDestroy( (cusparseHandle_t)ptr );
    TEUCHOS_TEST_FOR_EXCEPTION(status == CUSPARSE_STATUS_SUCCESS, std::runtime_error,
        "Kokkos::CUSPARSEdetails::CUSPARSEDestroyer::free(): library was never initialized; we should not have been called.")
  }
}

void Kokkos::CUSPARSEdetails::initCUSPARSEsession() 
{
  if (session_handle == null) {
    cusparseHandle_t hndl;
    cusparseStatus_t status = cusparseCreate(hndl);
    TEUCHOS_TEST_FOR_EXCEPTION(status == CUSPARSE_STATUS_NOT_INITIALIZED, std::runtime_error, 
        "Kokkos::CUSPARSEdetails::initCUSPARSEsession(): the CUDA Runtime initialization failed.")
    TEUCHOS_TEST_FOR_EXCEPTION(status == CUSPARSE_STATUS_ALLOC_FAILED, std::runtime_error, 
        "Kokkos::CUSPARSEdetails::initCUSPARSEsession(): the resources could not be allocated.")
    TEUCHOS_TEST_FOR_EXCEPTION(status == CUSPARSE_STATUS_ARCH_MISMATCH, std::runtime_error, 
        "Kokkos::CUSPARSEdetails::initCUSPARSEsession(): the device compute capability (CC) is less than 1.1.")
    TEUCHOS_TEST_FOR_EXCEPTION(status != CUSPARSE_STATUS_SUCCESS, std::runtime_error, 
        "Kokkos::CUSPARSEdetails::initCUSPARSEsession(): unspecified error.")
    session_handle = Teuchos::rcpWithDealloc(hndl,CUSPARSEDestroyer(),true);
  }
}
