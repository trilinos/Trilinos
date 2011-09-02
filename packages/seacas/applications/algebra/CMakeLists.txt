INCLUDE(SubPackageMacros)
INCLUDE(PackageAddExecutable)
INCLUDE(FortranSettings)

SUBPACKAGE(Algebra)

SET(CMAKE_EXECUTABLE_SUFFIX "")

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
SET(SOURCES 
   addvar.f   chkpar.f   dbinam.f   dbpqa.f    dooper.f   lnksto.f   myamod.f   postfx.f   rweval.f   sordel.f   unode.f
   chksyn.f   dbiqa.f    dbptim.f   eval.f     lnkvar.f   mymax.f    princ3.f   rwmap.f    sortid.f   wconat.f
   alicmd.f   cmdtim.f   dbivtt.f   dbvini.f   filtel.f   loceqv.f   mymin.f    prtdeb.f   rwnmap.f   store.f    welb.f
   aparse.f   cnvnum.f   dblist.f   delcmd.f   fixone.f   locnam.f   mysign.f   pxn.f      rwxyz.f    svcord.f   wnam.f
   caltim.f   comand.f   dboqa.f    dofnc1.f   idblnk.f   locrl.f    mysum.f    pxn2.f     savcmd.f   tmag.f     wrstep.f
   chkdis.f   cpyvar.f   dbpelb.f   dofnc2.f   initim.f   locstr.f   nwhsel.f   rdeqns.f   shotsp.f   uelem.f    zmess.f
   chkerr.f   dbicon.f   dbpini.f   dofncg.f   isize.f    makevo.f   nxttim.f   rdstep.f   show.f     uglob.f    zmfixd.f
   chkfnc.f   dbielb.f   dbpnam.f   doif.f     lnkfnc.f   mlist.f    opnlog.f   remel.f    shptim.f   uhist.f    zmnps.f)

PACKAGE_ADD_LIBRARY(algebra_app_lib SOURCES ${SOURCES} DEPLIBS supes suplib exodus_for)
PACKAGE_ADD_EXECUTABLE(algebra NOEXEPREFIX NOEXESUFFIX HEADERS ${HEADERS} SOURCES algebra.f DEPLIBS algebra_app_lib COMM serial)
install_executable(algebra)

SUBPACKAGE_POSTPROCESS()

