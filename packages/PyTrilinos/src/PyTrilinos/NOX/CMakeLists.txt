# The NOX module employs namespaces, so include the Epetra directory
IF(NOX_ENABLE_Epetra)
  ADD_SUBDIRECTORY(Epetra)
ENDIF(NOX_ENABLE_Epetra)
