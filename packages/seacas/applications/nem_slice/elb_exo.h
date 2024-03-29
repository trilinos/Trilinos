/*
 * Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#pragma once

#include "elb.h"

/* Function prototypes */
template <typename INT>
int read_exo_weights(Problem_Description *prob,   /* Pointer to problem info structure */
                     Weight_Description  *weight, /* Pointer to weight info structure */
                     INT                  dummy);

template <typename INT>
int read_mesh_params(const std::string     &exo_file, /* Name of ExodusII geometry file */
                     Problem_Description   *problem,  /* Pointer to problem info structure */
                     Mesh_Description<INT> *mesh,     /* Mesh information structure */
                     Sphere_Info           *sphere);            /* Sphere element info structure */

template <typename INT>
int read_mesh(const std::string     &exo_file, /* Name of ExodusII geometry file */
              Problem_Description   *problem,  /* Problem information */
              Mesh_Description<INT> *mesh,     /* Mesh information structure */
              Weight_Description    *weight);     /* Weight specification structure */

template <typename INT>
int init_weight_struct(Problem_Description   *problem, /* Problem information */
                       Mesh_Description<INT> *mesh,    /* Mesh information structure */
                       Weight_Description    *weight);    /* Weight specification structure */
