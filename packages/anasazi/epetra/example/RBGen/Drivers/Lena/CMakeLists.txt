

TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyLenaExampleFiles
  SOURCE_DIR ${Anasazi_SOURCE_DIR}/epetra/example/RBGen/Drivers/Lena
  SOURCE_FILES lenaInputAnasazi.xml lenaInputIncSVD.xml lenaInputStSVD.xml lenaInput.xml lena.mtx README
  )
