INCLUDE(SubPackageMacros)
INCLUDE(PackageLibraryMacros)

SUBPACKAGE(Supes)

INCLUDE(FortranSettings)

SET(SUPES_VERSION_MAJOR "4")
SET(SUPES_VERSION_MINOR "4")
SET(SUPES_VERSION_PATCH "0")
SET(SUPES_VERSION "${SUPES_VERSION_MAJOR}.${SUPES_VERSION_MINOR}")
SET(SUPES_VERSION_FULL "${SUPES_VERSION}.${SUPES_VERSION_PATCH}")

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(SOURCES
	ext_lib/exmemy.c      ext_lib/extime.c      ext_lib/ixlchr.c
	ext_lib/README        ext_lib/exname.c      ext_lib/exupcs.c  ext_lib/ixlnum.c
	ext_lib/excpus.c      ext_lib/exparm.c      ext_lib/getline_int.c
	ext_lib/exdate.c      ext_lib/exread.c      
	fre_fld/ffistr.f      fre_fld/frefld.f      fre_fld/getinp.f  fre_fld/gsupev.f
	fre_fld/quoted.f      fre_fld/stripb.f	    mem_mgr/mxexec.f  mem_mgr/mxlong.f	 mem_mgr/mycomp.f   mem_mgr/mymemy.f   mem_mgr/shftc.f
	mem_mgr/mdinit.F      mem_mgr/mxfind.f	    mem_mgr/mxlook.f  mem_mgr/mydel.f    mem_mgr/mynsrt.f   mem_mgr/shfti.f
	mem_mgr/mxcomp.f      mem_mgr/mxget.f	    mem_mgr/mxnsrt.f  mem_mgr/myfind.f   mem_mgr/myprnt.f   mem_mgr/srchc.f
	mem_mgr/mxdel.f	      mem_mgr/mxgive.f	    mem_mgr/mxprnt.f  mem_mgr/myget.f    mem_mgr/myrsrv.f   mem_mgr/srchi.f
	mem_mgr/mxeror.f      mem_mgr/mxlist.f	    mem_mgr/mxrsrv.f  mem_mgr/mylong.f   mem_mgr/vtable.f)

PACKAGE_ADD_LIBRARY(
	supes
	SOURCES ${SOURCES}
)

ADD_EXECUTABLE(exttest exttest.f)
ADD_EXECUTABLE(ffrtest ffrtest.f)
ADD_EXECUTABLE(memtest memtest.f)

TARGET_LINK_LIBRARIES(exttest supes)
TARGET_LINK_LIBRARIES(ffrtest supes)
TARGET_LINK_LIBRARIES(memtest supes)

SUBPACKAGE_POSTPROCESS()
