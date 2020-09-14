/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#include "fmt/ostream.h"
#include "nem_spread.h"     // for NemSpread
#include "ps_pario_const.h" // for PIO_Info, Parallel_IO
#include "rf_io_const.h"    // for ExoFile, Exo_LB_File, etc
#include <copy_string_cpp.h>
#include <cstdio>     // for stderr
#include <cstring>    // for strlen
#include <exodusII.h> // for ex_close, ex_open, EX_READ

template int NemSpread<double, int>::check_inp(void);
template int NemSpread<float, int>::check_inp(void);
template int NemSpread<double, int64_t>::check_inp(void);
template int NemSpread<float, int64_t>::check_inp(void);

template <typename T, typename INT> int NemSpread<T, INT>::check_inp()
{
  int   exid;
  int   icpu_ws = 0;
  int   iio_ws  = 0;
  float vers    = 0.0;

  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  /*                 Check the input and output files                          */
  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  /* check if the Mesh file was specified */
  if (ExoFile.empty()) {
    fmt::print(stderr, "{}: fatal - must specify input FEM file.\n", __func__);
    return 0;
  }

  /* check for the existence of a readable FEM file */
  int mode = EX_READ | int64api;
  if ((exid = ex_open(ExoFile.c_str(), mode, &icpu_ws, &iio_ws, &vers)) < 0) {
    fmt::print(stderr, "{}: fatal - unable to open input FEM file, {}.\n", __func__,
               ExoFile.c_str());
    return 0;
  }
  ex_close(exid);

  /* check that there is a load balance file specified */
  if (Exo_LB_File.empty()) {
    fmt::print(stderr, "{}: fatal - must specify input FEM file.\n", __func__);
    return 0;
  }

  /* check for the existence of a readable load balance file */
  icpu_ws = 0;
  iio_ws  = 0;
  if ((exid = ex_open(Exo_LB_File.c_str(), mode, &icpu_ws, &iio_ws, &vers)) < 0) {
    fmt::print(stderr, "{}: fatal - unable to open load balance file, {}.\n", __func__,
               Exo_LB_File.c_str());
    return 0;
  }
  ex_close(exid);

  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  /*               Check the result spreading specifications                   */
  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

  /* check if anything was specified for restart information */
  if (Restart_Info.Flag < 0) {
    /* default is now to spread results if they exist */
    Restart_Info.Flag      = 1;
    Restart_Info.Num_Times = -1; /* -1 means spread all results */
  }

  /* check to see if there is a separate restart file */
  if (Restart_Info.Flag > 0) {
    if (Exo_Res_File.empty()) {
      Exo_Res_File = ExoFile; /* if not use the input FEM file */
    }
  }

  /* check if space is to be reserved for variables in the parallel files */
  if (Num_Glob_Var < 0) {
    Num_Glob_Var = 0;
  }
  if (Num_Nod_Var < 0) {
    Num_Nod_Var = 0;
  }
  if (Num_Elem_Var < 0) {
    Num_Elem_Var = 0;
  }
  if (Num_Nset_Var < 0) {
    Num_Nset_Var = 0;
  }
  if (Num_Sset_Var < 0) {
    Num_Sset_Var = 0;
  }

  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  /*                 Check the parallel IO specifications                      */
  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

  /* default is not to have preceding 0's in the disk names */
  if (PIO_Info.Zeros < 0) {
    PIO_Info.Zeros = 0;
  }
  /* most systems that we deal with start their files systems with 1 not 0 */
  if (PIO_Info.PDsk_Add_Fact < 0) {
    PIO_Info.PDsk_Add_Fact = 1;
  }

  /* check that there is a list of disks, or a number of raids */
  if ((PIO_Info.Dsk_List_Cnt <= 0) && (PIO_Info.Num_Dsk_Ctrlrs <= 0)) {
    fmt::print(stderr,
               "{}: fatal - must specify a number of raids, or a disk"
               " list.\n",
               __func__);
    return 0;
  }

  if (PIO_Info.Par_Dsk_Root.empty()) {
    fmt::print(stderr,
               "{}: Error - Root directory for parallel files must"
               " be specified.\n",
               __func__);
    return 0;
  }

  if (PIO_Info.Par_Dsk_SubDirec.empty()) {
    fmt::print(stderr,
               "{}: Error - Subdirectory for parallel files must"
               " be specified.\n",
               __func__);
    return 0;
  }

  return 1;
}
