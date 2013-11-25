INCLUDE(TribitsSubPackageMacros)
INCLUDE(TribitsLibraryMacros)
INCLUDE(FortranSettings)

TRIBITS_SUBPACKAGE(SVDI)

SET(PLT_VERSION_MAJOR "1")
SET(PLT_VERSION_MINOR "0")
SET(PLT_VERSION_PATCH "0")
SET(PLT_VERSION "${PLT_VERSION_MAJOR}.${PLT_VERSION_MINOR}")
SET(PLT_VERSION_FULL "${PLT_VERSION}.${PLT_VERSION_PATCH}")

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})

SET(CGI_SOURCES "")
APPEND_GLOB(CGI_SOURCES
	${DIR}/cgi/cgisdum.c
	${DIR}/cgi/vdicgi.c
        ${DIR}/cgi/sdcgif.f
        ${DIR}/cgi/vdicgi_char.f
        ${DIR}/cgi/mdcgi.c
        ${DIR}/cgi/sdcgi.c
        ${DIR}/cgi/x11_vdix11.c
        ${DIR}/cgi/x11_x11xlate.c
        ${DIR}/cgi/pst_vdipst.f
        ${DIR}/cgi/pst_pstxlate.c
        ${DIR}/cgi/met_vdimet.f
        ${DIR}/cgi/met_metxlate.c
        ${DIR}/cgi/svdi_addrwrap.F)

TRIBITS_ADD_LIBRARY(
	svdi_cgi
	SOURCES ${CGI_SOURCES}
	DEPLIBS svdi_cdr
)

SET(CDR_SOURCES "")
APPEND_GLOB(CDR_SOURCES ${DIR}/cdr/*.f ${DIR}/cdr/*.F ${DIR}/cdr/*.c)

TRIBITS_ADD_LIBRARY(
	svdi_cdr
	SOURCES ${CDR_SOURCES}
)

TRIBITS_ADD_LIBRARY(svdi_cgi_x11
	SOURCES ${DIR}/cgi/x11_cgisx11.c
	DEPLIBS svdi_cgi
)

TRIBITS_ADD_LIBRARY(svdi_cgi_pst
	SOURCES ${DIR}/cgi/pst_cgispst.c
	DEPLIBS svdi_cgi
)

TRIBITS_ADD_LIBRARY(svdi_cgi_met
	SOURCES ${DIR}/cgi/met_cgismet.c
	DEPLIBS svdi_cgi
)

TRIBITS_ADD_LIBRARY(vdx11cps
	SOURCES ${DIR}/vdi_drivers/vdx11cps.F
	DEPLIBS vdicps_dual svdi_cgi
)

TRIBITS_ADD_LIBRARY(vdicps_dual
	SOURCES ${DIR}/vdi_drivers/vdicps_dual.f
	DEPLIBS svdi_cgi
)

TRIBITS_ADD_LIBRARY(svdi_cgi_cps
	SOURCES ${DIR}/vdi_drivers/vdicps.f
	DEPLIBS svdi_cgi
)

TRIBITS_SUBPACKAGE_POSTPROCESS()

