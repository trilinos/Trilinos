/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#ifndef _ELB_INP_CONST_H_
#define _ELB_INP_CONST_H_

#include <string> // for string
struct Machine_Description;
struct Problem_Description;
struct Solver_Description;
template <typename INT> struct LB_Description;
template <typename INT> struct Weight_Description;

/* Prototype for command-line parsing function */
template <typename INT>
int cmd_line_arg_parse(
    int                  argc,           /* The command line argument count */
    char *               argv[],         /* The command line arguments array */
    std::string &        exoII_inp_file, /* The ExodusII input FEM file name */
    std::string &        ascii_inp_file, /* The ASCII input file name */
    std::string &        nemI_out_file,  /* The output NemesisI file name */
    Machine_Description *machine,        /* Pointer to structure in which to place machine
                                          * information */
    LB_Description<INT> *lb,             /* Pointer to structure in which to place load
                                          * balance parameters */
    Problem_Description *prob,           /* Pointer to structure in which to place general
                                          * information about the run */
    Solver_Description *solver,          /* Pointer to structure in which to place parameters
                                          * for the eigensolver */
    Weight_Description<INT> *weight      /* Pointer to structure in which to place parameters
                                          * for the graph weighting scheme */
);

/* Prototype for function which reads in the ASCII input file */
template <typename INT>
int read_cmd_file(std::string &        ascii_inp_file, /* The ASCII input file name */
                  std::string &        exoII_inp_file, /* The ExodusII input FEM file name */
                  std::string &        nemI_out_file,  /* The output NemesisI file name */
                  Machine_Description *machine, /* Pointer to structure in which to place machine
                                                 * information */
                  LB_Description<INT> *lb,      /* Pointer to structure in which to place load
                                                 * balance parameters */
                  Problem_Description *problem, /* Pointer to structure in which to place general
                                                 * information about the run */
                  Solver_Description *solver,   /* Pointer to structure in which to place parameters
                                                 * for the eigensolver */
                  Weight_Description<INT> *weight /* Pointer to structure in which to place
                                                   * parameters for the eigensolver */
);

/* Prototype for function which checks the user specified input */
template <typename INT>
int check_inp_specs(std::string &        exoII_inp_file, /* The ExodusII input FEM file name */
                    std::string &        nemI_out_file,  /* The output NemesisI file name */
                    Machine_Description *machine, /* Pointer to structure in which to place machine
                                                   * information */
                    LB_Description<INT> *lb,      /* Pointer to structure in which to place load
                                                   * balance parameters */
                    Problem_Description *prob,    /* Pointer to structure in which to place general
                                                   * information about the run */
                    Solver_Description *solver, /* Pointer to structure in which to place parameters
                                                 * for the eigensolver */
                    Weight_Description<INT> *weight /* Pointer to structure in which to place
                                                     * parameters for the weighting scheme */
);

/* Various defines used by the input routines */
#define NS_NONE -1

#endif /* _ELB_INP_CONST_H_ */
