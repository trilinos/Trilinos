/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*
 *
 *
 *
 *      Include file for I/O global variables used in FEM
 *      problem specification
 *
 */

std::string ExoFile;               /* Exodus II File containing problem definition. */
                                   /* This name is the root name.                   */
std::string Output_File_Base_Name; /* Base name of output file. If it has a suffix, it will be
                                        stripped */

std::string Exo_LB_File;  /* Exodus II file containing the mesh
                           * load-balanceinformation                       */
std::string Exo_Res_File; /* Exodus II file containing the mesh results  */
                          /* information                                   */

int Debug_Flag = 1; /* Flag to specify debug info is to be printed out.
                       The value of this flag determines the level of
                       diagnostic info which is printed to stdout
                       Debug_Flag == 0  No debug output
                                     .
                                     .
                                     9  maximum output               */
int Gen_Flag = 1;   /* Flag used by nem_join to determine if the user
                       wants to use an existing genesis file rather
                       than generating one from the parallel files */

int Num_Nod_Var = 0;  /* The number of nodal variables to reserve */
                      /* space in the output file for. */
int Num_Elem_Var = 0; /* The number of elemental variables to */
                      /* reserve space in the output file for. */
int Num_Glob_Var = 0; /* The number of global variables to reserve */
                      /* space in the output file for. */
int Num_Nset_Var = 0; /* The number of nodeset variables to reserve */
                      /* space in the output file for. */
int Num_Sset_Var = 0; /* The number of sideset variables to reserve */
                      /* space in the output file for. */
