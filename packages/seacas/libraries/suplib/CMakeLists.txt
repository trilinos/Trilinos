INCLUDE(SubPackageMacros)
INCLUDE(PackageLibraryMacros)

SUBPACKAGE(Suplib)
INCLUDE(FortranSettings)

SET(SUPLIB_VERSION_MAJOR "3")
SET(SUPLIB_VERSION_MINOR "08")
SET(SUPLIB_VERSION_PATCH "0")
SET(SUPLIB_VERSION "${SUPLIB_VERSION_MAJOR}.${SUPLIB_VERSION_MINOR}")
SET(SUPLIB_VERSION_FULL "${SUPLIB_VERSION}.${SUPLIB_VERSION_PATCH}")

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(SOURCES
	cpyrgt.f          dboelb.f          flush.f           locint.f          pltond.f
	cpystr.f          dboess.f          f2kcli.F          gslibv.f          locrea.f          pltone.f
	dberr.f           dboini.f          f2kcli.inc        help.F            locstr.f          prterr.f
	dbieb1.f          dbomap.f          ffaddc.f          help.h            logerr.f          scanner.c
	abrstr.f          dbielb.f          dbonam.f          ffaddi.f          help_c.c          lowstr.f          scneof.f
	addlog.c          dbiess.f          dbonm1.f          ffaddo.f          idblnk.f          matstr.f          screrr.f
	addv.f            dbiini.f          dbonps.f          ffaddr.f          iengrx.f          memerr.f          shocmd.f
	banner.f          dbimap.f          dboqa.f           ffaddv.f          indexr.f          minmax.f          sqzstr.f
	bannr2.f          dbinam.f          dboste.f          ffchar.f          inigrd.f          minmxl.f          strtup.f
	batch.f           dbinm1.f          dboxyz.f          ffexst.f          iniint.f          minmxs.f          subv.f
	c_xhelp.c         dbinps.f          dbpelb.f          ffintg.f          inilog.f          num2ix.f          util.c
	ckexod.f          dbiqa.f           dbpini.f          ffmatc.f          inirea.f          numeqi.f          util.h
	cknone.f          dbist1.f          dbpnam.f          ffneed.f          inistr.f          numeql.f          wheneq.f
	ckwhol.f          dbist2.f          dbpqa.f           ffnumb.f          intadd.f          numstr.f
	copint.f          dbiste.f          dbptim.f          ffonof.f          intstr.f          okabrt.f          wrapup.f
	cpuifc.f          dbitim.f          dbvini.f          ffreal.f          iowdsz.f          opnfil.f          zeroin.f
	cpyint.f          dbiv0.f           do_help.c         ffrrng.f          ismax.f           packt.f
	cpylog.f          dbixyz.f          engnot.f          ismin.f           pcklab.f 
	cpyrea.f          dboeb1.f          engstr.f          ffvrng.f          lenstr.f          pckstr.f)

PACKAGE_ADD_LIBRARY(
	suplib
	SOURCES ${SOURCES}
	DEPLIBS supes
)

SUBPACKAGE_POSTPROCESS()
