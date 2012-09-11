SET(SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
        commsplitter    libraries/commsplitter         EX  OPTIONAL
        support         libraries/support              EX  OPTIONAL
        nnti            libraries/nessie/nnti          EX  OPTIONAL
        nssi            libraries/nessie/nssi          EX  OPTIONAL
        programs        programs                       EX  OPTIONAL
        examples        examples                       EX  OPTIONAL
        tests           tests                          EX  OPTIONAL
        netcdf-service  services/netcdf                EX  OPTIONAL
)

SET(LIB_REQUIRED_DEP_PACKAGES )
SET(LIB_OPTIONAL_DEP_PACKAGES Teuchos)
SET(TEST_REQUIRED_DEP_PACKAGES)
SET(TEST_OPTIONAL_DEP_PACKAGES)
SET(LIB_REQUIRED_DEP_TPLS)
SET(LIB_OPTIONAL_DEP_TPLS MPI Netcdf Pnetcdf Pthread CrayPortals Portals Gemini InfiniBand Pablo HPCToolkit)
SET(TEST_REQUIRED_DEP_TPLS)
SET(TEST_OPTIONAL_DEP_TPLS)
