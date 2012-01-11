
INCLUDE(TribitsLibraryMacros)
INCLUDE(GenerateXDR)

#  This CMakeLists.txt file is only used to generate XDR header and 
#  source code using rpcgen 
#

GenerateXDR(nnti_xdr)
GenerateXDR(nssi_types_xdr)
GenerateXDR(nssi_service_args)
