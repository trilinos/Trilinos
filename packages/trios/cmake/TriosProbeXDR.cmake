
########## CHECK FOR HEADER FILES ############

INCLUDE(CheckIncludeFiles)

########## CHECK FOR FUNCTIONS ############

INCLUDE(CheckLibraryExists)
INCLUDE(CheckFunctionExists)
INCLUDE(CheckCSourceCompiles)
INCLUDE(CheckCXXSourceCompiles)

CHECK_FUNCTION_EXISTS(xdr_u_int16_t HAVE_TRIOS_XDR_U_INT16_T)
CHECK_FUNCTION_EXISTS(xdr_u_int32_t HAVE_TRIOS_XDR_U_INT32_T)
CHECK_FUNCTION_EXISTS(xdr_u_int64_t HAVE_TRIOS_XDR_U_INT64_T)

# XDR_SIZEOF
check_c_source_compiles(
    "#include <rpc/xdr.h>\nint main(){xdr_sizeof(NULL,NULL);return 0;}"
    HAVE_TRIOS_XDR_SIZEOF
)
