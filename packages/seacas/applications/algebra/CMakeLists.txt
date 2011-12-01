INCLUDE(TribitsSubPackageMacros)
INCLUDE(TribitsAddExecutable)
INCLUDE(FortranSettings)

TRIBITS_SUBPACKAGE(Algebra)

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
SET(SOURCES 
   addvar.f   chkpar.f   ag_dbinam.f   dooper.f   lnksto.f   myamod.f   postfx.f   rweval.f   sordel.f   unode.f
   chksyn.f   ag_dbiqa.f    ag_dbptim.f   eval.f     lnkvar.f   mymax.f    princ3.f   rwmap.f    sortid.f   wconat.f
   alicmd.f   cmdtim.f   ag_dbivtt.f   ag_dbvini.f   filtel.f   loceqv.f   mymin.f    prtdeb.f   rwnmap.f   store.f    welb.f
   aparse.f   cnvnum.f   ag_dblist.f   delcmd.f   fixone.f   locnam.f   mysign.f   pxn.f      rwxyz.f    svcord.f   wnam.f
   caltim.f   comand.f   ag_dboqa.f    dofnc1.f   locrl.f    mysum.f    pxn2.f     savcmd.f   tmag.f     wrstep.f
   chkdis.f   cpyvar.f   ag_dbpelb.f   dofnc2.f   initim.f   nwhsel.f   rdeqns.f   shotsp.f   uelem.f    zmess.f
   chkerr.f   ag_dbicon.f   ag_dbpini.f   dofncg.f   isize.f    makevo.f   nxttim.f   rdstep.f   show.f     uglob.f    zmfixd.f
   chkfnc.f   ag_dbielb.f   ag_dbpnam.f   doif.f     lnkfnc.f   mlist.f    opnlog.f   remel.f    shptim.f   uhist.f    zmnps.f)

TRIBITS_ADD_LIBRARY(algebra_app_lib SOURCES ${SOURCES} DEPLIBS supes suplib exodus_for)
TRIBITS_ADD_EXECUTABLE(algebra NOEXEPREFIX NOEXESUFFIX HEADERS ${HEADERS} LINKER_LANGUAGE Fortran SOURCES algebra.f DEPLIBS algebra_app_lib COMM serial mpi)
install_executable(algebra)

TRIBITS_SUBPACKAGE_POSTPROCESS()

