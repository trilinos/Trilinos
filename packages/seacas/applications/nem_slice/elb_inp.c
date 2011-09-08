/*
 * Copyright (C) 2009 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 *----------------------------------------------------------------------------
 * Functions contained in this file:
 *	cmd_line_arg_parse()
 *	read_cmd_file()
 *	check_inp_specs()
 *+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <stdlib.h>
#include <exodusII.h>

#include "getopt.h"
#include "md_getsubopt.h"
#include "elb_inp_const.h"
#include "elb_err_const.h"
#include "elb_util_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* This function parses the command line options and stores the information
 * in the appropriate data locations.
 *---------------------------------------------------------------------------*/
int cmd_line_arg_parse(
  int argc, char *argv[],	/* Args as passed by main() */
  char exoII_inp_file[],	/* The input ExodusII file name */
  char ascii_inp_file[],	/* The ASCII input file name */
  char nemI_out_file[],		/* Output NemesisI file name */
  MACHINE_PTR machine,		/* Structure for machine description */
  LB_INFO_PTR lb,		/* Structure for load balance description */
  PROB_INFO_PTR prob,		/* Structure for various problem params */
  SOLVE_INFO_PTR solver,	/* Structure for eigen solver params */
  WEIGHT_INFO_PTR weight	/* Structure for weighting graph */
  )
{
  int   opt_let, iret, el_blk, wgt, max_dim=0, i;
  char *sub_opt=NULL, *value=NULL, *cptr=NULL, *cptr2=NULL, ctemp[MAX_ERR_MSG+1];

  extern char *optarg;
  extern int   optind;

  /* see NOTE in elb_const.h about the order of the following array */
  char *weight_subopts[] = {
     "none",
     "read",
     "eb",
     "var_index",
     "edges",
     "time_index",
     "var_name",
     NULL
  };

  char *mach_subopts[] = {
    "mesh",
    "hcube",
    "hypercube",
    "cluster",
    NULL
  };

  char *lb_subopts[] = {
    "multikl",
    "spectral",
    "inertial",
    "linear",
    "random",
    "scattered",
    "infile",
    "kl",
    "none",
    "num_sects",
    "cnctd_dom",
    "outfile",
    "zpinch",
    "brick",
    "rcb",
    "rib",
    "hsfc",
    NULL
  };

  char *solve_subopts[] = {
    "tolerance",
    "use_rqi",
    "vmax",
    NULL
  };

/*---------------------------Execution Begins--------------------------------*/

  /*
   * Make sure there were command line options given. If not assign the
   * name of the default ascii input file.
   */
  if(argc <= 1)
  {
    print_usage();
    exit(0);
  }

  /* Loop over each command line option */
  while((opt_let=getopt(argc, argv, "a:hm:l:nes:x:w:vo:cg:fp")) != EOF) {

    /* case over the option letter */
    switch(opt_let)
    {
    case 'v':
      /* Should an ouput visualization file be output */
      if(prob->vis_out == ELB_TRUE)
        prob->vis_out = ELB_FALSE;
      else
        prob->vis_out =  ELB_TRUE;

      break;

    case 'x':
      /* Undocumented flag for setting the error level */
      if(optarg == NULL)
        error_lev = 1;
      else
      {
        iret = sscanf(optarg, "%d", &error_lev);
        if(iret != 1)
          error_lev = 1;

        if(error_lev > 3)
          error_lev = 3;
      }
      break;

    case 'c':
      /* flag to allow the user to skip some error checks */
      prob->skip_checks = 1;
      break;

    case 'f':
      /*
       * use the face method to calculate adjacencies for
       * elemental decompositions
       */
      prob->face_adj = 1;
      break;

    case 'p':
      /*
       * use the partial method to determine adjacencies:
       * only 3/4 of face nodes must match instead of all
       */
       prob->partial_adj = 1;
       break;

    case 'w':
      /* Weighting options */
      sub_opt = optarg;
      while(*sub_opt != '\0')
      {
        switch(md_getsubopt(&sub_opt, weight_subopts, &value))
        {
        case READ_EXO:
          if(value == NULL)
          {
            sprintf(ctemp, "fatal: must specify a file name with %s",
                    weight_subopts[READ_EXO]);
            Gen_Error(0, ctemp);
            return 0;
          }
          if(strlen(value) == 0)
          {
            sprintf(ctemp, "fatal: must specify a file name with %s",
                    weight_subopts[READ_EXO]);
            Gen_Error(0, ctemp);
            return 0;
          }

          if (weight->type < 0)                weight->type = READ_EXO;
          else if (!(weight->type & READ_EXO)) weight->type += READ_EXO;

          /* check if the read is after an element block weight */
          if (weight->type & EL_BLK) weight->ow_read = 0;

          strcpy(weight->exo_filename, value);

          break; /* End "case READ_EXO" */

        case VAR_INDX:
          if(value == NULL)
          {
            sprintf(ctemp, "fatal: must specify a value with %s",
                    weight_subopts[VAR_INDX]);
            Gen_Error(0, ctemp);
            return 0;
          }
          if(strlen(value) == 0)
          {
            sprintf(ctemp, "fatal: must specify a value with %s",
                    weight_subopts[VAR_INDX]);
            Gen_Error(0, ctemp);
            return 0;
          }

          iret = sscanf(value, "%d", &(weight->exo_vindx));
          if(iret != 1)
          {
            sprintf(ctemp, "fatal: invalid value specified for %s",
                    weight_subopts[VAR_INDX]);
            Gen_Error(0, ctemp);
            return 0;
          }
          break;

        case TIME_INDX:
          if(value == NULL)
          {
            sprintf(ctemp, "fatal: must specify a value with %s",
                    weight_subopts[TIME_INDX]);
            Gen_Error(0, ctemp);
            return 0;
          }
          if(strlen(value) == 0)
          {
            sprintf(ctemp, "fatal: must specify a value with %s",
                    weight_subopts[TIME_INDX]);
            Gen_Error(0, ctemp);
            return 0;
          }

          iret = sscanf(value, "%d", &(weight->exo_tindx));
          if(iret != 1)
          {
            sprintf(ctemp, "fatal: invalid value specified for %s",
                    weight_subopts[TIME_INDX]);
            Gen_Error(0, ctemp);
            return 0;
          }
          break;

        case VAR_NAME:
          if(value == NULL)
          {
            sprintf(ctemp, "fatal: must specify a value with %s",
                    weight_subopts[VAR_NAME]);
            Gen_Error(0, ctemp);
            return 0;
          }
          if(strlen(value) == 0)
          {
            sprintf(ctemp, "fatal: must specify a value with %s",
                    weight_subopts[VAR_NAME]);
            Gen_Error(0, ctemp);
            return 0;
          }

          strcpy(weight->exo_varname, value);

          break;

        case EL_BLK:
          if(value == NULL)
          {
            sprintf(ctemp, "fatal: must specify an element block and weight with %s",
                    weight_subopts[EL_BLK]);
            Gen_Error(0, ctemp);
            return 0;
          }
          el_blk = -1;
          wgt = -1;

          iret = sscanf(value, "%d:%d", &el_blk, &wgt);
          if(iret != 2)
          {
            Gen_Error(0, "invalid element block weight");
            return 0;
          }
          if(el_blk <= 0)
          {
            printf(ctemp, "invalid element block, %d", el_blk);
            Gen_Error(0, ctemp);
            return 0;
          }
          if(wgt < 0)
          {
            printf(ctemp, "invalid weight, %d", wgt);
            Gen_Error(0, ctemp);
            return 0;
          }

          /* check if this is the first element block weight */
          weight->num_ebw++;
          if (weight->num_ebw == 1) {
            if (weight->type < 0) weight->type = EL_BLK;
            else                  weight->type += EL_BLK;

            /* since this is the first time through allocate array */
            weight->elemblk = (int *) malloc(sizeof(int));
            weight->elemblk_wgt = (int *) malloc(sizeof(int));
            if(!weight->elemblk || !weight->elemblk_wgt)
            {
              Gen_Error(0, "fatal: insufficient memory");
              return 0;
            }
          }
          else { /* just realloc arrays */
            weight->elemblk = (int *) realloc(weight->elemblk,
                                              weight->num_ebw*sizeof(int));
            weight->elemblk_wgt = (int *) realloc(weight->elemblk_wgt,
                                                  weight->num_ebw*sizeof(int));
            if(!weight->elemblk || !weight->elemblk_wgt)
            {
              Gen_Error(0, "fatal: insufficient memory");
              return 0;
            }
          }
          /* now set the values */
          weight->elemblk[weight->num_ebw-1] = el_blk;
          weight->elemblk_wgt[weight->num_ebw-1] = wgt;

          if (weight->type < 0)              weight->type = EL_BLK;
          else if (!(weight->type & EL_BLK)) weight->type += EL_BLK;

          /* check if the element block weight needs to over write the read */
          if (weight->type & READ_EXO) weight->ow_read = 1;

          break;

        case EDGE_WGT:
          if (weight->type < 0)                weight->type = EDGE_WGT;
          else if (!(weight->type & EDGE_WGT)) weight->type += EDGE_WGT;
          break;

        case NO_WEIGHT:
          weight->type = NO_WEIGHT;
          break;

        default:
          sprintf(ctemp, "fatal: unknown suboption %s specified", value);
          Gen_Error(0, ctemp);
          return 0;

        } /* End "switch(md_getsubopt(&sub_opt, weight_subopts, &value))" */

      } /* End "while(*sub_opt != '\0')" */

      break; /* End "case 'w'" */

    case 'a':
      /* Only an ASCII input file name */
      strcpy(ascii_inp_file, optarg);
      break;

    case 'o':
      /* Output NemesisI file name */
      strcpy(nemI_out_file, optarg);
      break;

    case 'n':
      /* Nodal decomposition */
      if(prob->type == ELEMENTAL)
      {
        Gen_Error(0, "fatal: -e and -n are mutually exclusive");
        return 0;
      }
      prob->type = NODAL;
      break;

    case 'e':
      /* Elemental decomposition */
      if(prob->type == NODAL)
      {
        Gen_Error(0, "fatal: -e and -n are mutually exclusive\n");
        return 0;
      }
      prob->type = ELEMENTAL;
      break;

    case 'm':
      /* Machine type */
      sub_opt = optarg;
      string_to_lower(sub_opt, '\0');
      while(*sub_opt != '\0')
      {

        /* Switch over the machine description */
        switch(md_getsubopt(&sub_opt, mach_subopts, &value))
        {
        case HCUBE:
        case HYPERCUBE:
          if (machine->type < 0)
          {
            machine->type = HCUBE;
            max_dim = 1;
          }
          /* fall thru */

        case MESH:
          if (machine->type < 0)
          {
            machine->type = MESH;
            max_dim = 3;
          }

          cptr = value;  /* want to set this for both mesh and hcube */

          /* fall thru */

        case CLUSTER:
          if (machine->type < 0)  /* so, get the number of boxes */
          {
            if(value == NULL || strlen(value) == 0)
            {
              Gen_Error(0, "fatal: need to specify number of boxes");
              return 0;
            }

            /* now need to find what each box consists of */
            cptr = strpbrk(value, "mMhH");
            if (*cptr == 'm' || *cptr == 'M') {
              machine->type = MESH;
              max_dim = 3;
            }
            else if (*cptr == 'h' || *cptr == 'H') {
              machine->type = HCUBE;
              max_dim = 1;
            }
            else
            {
              Gen_Error(0, "fatal: unknown type specified with cluster");
              return 0;
            }
            /* blank out character and move cptr to next char */
            *cptr = '\0';
            cptr++;

            /* get the number of boxes from value */
            iret = sscanf(value, "%d", &(machine->num_boxes));
            if(iret <= 0 || machine->num_boxes <= 0)
            {
              Gen_Error(0, "fatal: invalid number of boxes");
              return 0;
            }
          }

          if(cptr == NULL || strlen(cptr) == 0)
          {
            Gen_Error(0, "fatal: need to specify dimension");
            return 0;
          }
          cptr2 = strtok(cptr, "xX");
          if(cptr2 == NULL)
          {
            Gen_Error(0, "fatal: bad size for dimension specification");
            return 0;
          }
          machine->num_dims = 0;
          for (i = 0; i < max_dim; i++) machine->dim[i] = 1;
          while(cptr2)
          {
            iret = sscanf(cptr2, "%d", &(machine->dim[machine->num_dims]));
            if(iret <= 0 || machine->dim[machine->num_dims] <= 0)
            {
              Gen_Error(0, "fatal: invalid dimension specification");
              return 0;
            }

            machine->num_dims++;
            cptr2 = strtok(NULL, "xX");

            /* Only up to three-dimensional allowed */
            if(machine->num_dims == max_dim && cptr2 != NULL)
            {
              Gen_Error(0, "fatal: maximum number of dimensions exceeded");
              return 0;
            }
          }

          break; /* End "case MESH or HCUBE or CLUSTER" */

        default:
          Gen_Error(0, "fatal: unknown machine type");
          return 0;

        } /* End "switch(md_getsubopt(&sub_opt, mach_subopts, &value))" */

      } /* End "while(*sub_opt != '\0')" */

      break; /* End "case 'm'" */

    case 'l':
      /* Load balance information */
      sub_opt = optarg;
      string_to_lower(sub_opt, '\0');
      while(*sub_opt != '\0')
      {
        switch(md_getsubopt(&sub_opt, lb_subopts, &value))
        {
        case MULTIKL:
          lb->type = MULTIKL;
          break;

        case SPECTRAL:
          lb->type = SPECTRAL;
          break;

        case INERTIAL:
          lb->type = INERTIAL;
          break;

        case ZPINCH:
          lb->type = ZPINCH;
          break;

        case BRICK:
          lb->type = BRICK;
          break;

        case ZOLTAN_RCB:
          lb->type = ZOLTAN_RCB;
          break;

        case ZOLTAN_RIB:
          lb->type = ZOLTAN_RIB;
          break;

        case ZOLTAN_HSFC:
          lb->type = ZOLTAN_HSFC;
          break;

        case LINEAR:
          lb->type = LINEAR;
          break;

        case RANDOM:
          lb->type = RANDOM;
          break;

        case SCATTERED:
          lb->type = SCATTERED;
          break;

        case INFILE:
          if(value == NULL)
          {
            Gen_Error(0, "fatal: need to specify a value with file");
            return 0;
          }
          iret = sscanf(value, "%s", lb->file);
          if(iret != 1)
          {
            Gen_Error(0, "fatal: invalid value associated with file");
            return 0;
          }
          lb->type = INFILE;
          break;

        case KL_REFINE:
          lb->refine = KL_REFINE;
          break;

        case NO_REFINE:
          lb->refine = NO_REFINE;
          break;

        case NUM_SECTS:
          if(value == NULL)
          {
            Gen_Error(0, "fatal: need to specify a value with num_sects");
            return 0;
          }
          iret = sscanf(value, "%d", &(lb->num_sects));
          if(iret != 1)
          {
            Gen_Error(0, "fatal: invalid value associated with num_sects");
            return 0;
          }
          break;

        case CNCT_DOM:
          lb->cnctd_dom = 1;
          break;

        case OUTFILE:
          if(value == NULL)
          {
            Gen_Error(0, "fatal: need to specify a value with outfile");
            return 0;
          }
          iret = sscanf(value, "%s", lb->file);
          if(iret != 1)
          {
            Gen_Error(0, "fatal: invalid value associated with outfile");
            return 0;
          }
          lb->outfile = ELB_TRUE;
          break;

        default:
          sprintf(ctemp, "fatal: unknown lb param \"%s\"", value);
          Gen_Error(0, ctemp);
          return 0;

        } /* End "switch(md_getsubopt(&sup_opt, mach_subopts, &value))" */

      } /* End "while(*sup_opt != '\0')" */

      break; /* End "case 'l'" */

    case 's':
      /* Eigen solver options */
      sub_opt = optarg;
      string_to_lower(sub_opt, '\0');
      while(*sub_opt != '\0')
      {
        switch(md_getsubopt(&sub_opt, solve_subopts, &value))
        {
        case TOLER:
          if(value == NULL)
          {
            fprintf(stderr, "fatal: tolerance specification requires \
value\n");
            return 0;
          }
          iret = sscanf(value, "%le", &(solver->tolerance));
          if(iret != 1)
          {
            fprintf(stderr, "fatal: incorrect value for tolerance\n");
            return 0;
          }

          break;

        case USE_RQI:
          if(solver->rqi_flag == USE_RQI)
            solver->rqi_flag = -1;
          else
            solver->rqi_flag = USE_RQI;

          break;

        case VMAX:
          if(value == NULL)
          {
            fprintf(stderr, "fatal: must specify a value with %s\n",
                    solve_subopts[VMAX]);
            return 0;
          }
          iret = sscanf(value, "%d", &(solver->vmax));
          if(iret != 1)
          {
            fprintf(stderr, "fatal: invalid value read for %s\n",
                    solve_subopts[VMAX]);
            return 0;
          }

          break;

        default:
          fprintf(stderr, "fatal: unknown solver option\n");
          return 0;

        } /* End "switch(md_getsubopt(&sub_opt, solve_subopts, &value))" */

      } /* End "while(sub_opt != '\0')" */

      break; /* End "case 's'" */

    case 'g':
      /* group designations */
      /* allocate string to hold designation */
      prob->groups = malloc(strlen(optarg) + 1);
      strcpy(prob->groups, optarg);
      break;

    case 'h':
      /* Usage info was requested */
      print_usage();
      exit(0);

    default:
      /* Default case. Error on unknown argument. */
      return 0;

    } /* End "switch(opt_let)" */

  } /* End "while((opt_let=getopt(argc, argv, "i")) != EOF)" */

  /* Get the input file name, if specified on the command line */
  if((argc-optind) >= 1)
    strcpy(exoII_inp_file, argv[optind]);
  else
    exoII_inp_file[0] = '\0';

  return 1;

} /*-------End cmd_line_arg_parse()--------*/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* This function reads in the ASCII command file.
 *---------------------------------------------------------------------------*/
int read_cmd_file(
  char ascii_inp_file[],
  char exoII_inp_file[],
  char nemI_out_file[],
  MACHINE_PTR machine,
  LB_INFO_PTR lb,
  PROB_INFO_PTR problem,
  SOLVE_INFO_PTR solver,
  WEIGHT_INFO_PTR weight
  )
{
  FILE *inp_fd;
  char  ctemp[MAX_ERR_MSG+1], inp_line[MAX_INP_LINE];
  char  inp_copy[MAX_INP_LINE];
  char *cptr, *cptr2;

  int iret, el_blk, wgt, i, ilen, max_dim;
/*-----------------------------Execution Begins------------------------------*/
  if(!(inp_fd=fopen(ascii_inp_file, "r")))
  {
    sprintf(ctemp, "fatal: unable to open ASCII input file %s",
            ascii_inp_file);
    Gen_Error(0, ctemp);
    return 0;
  }

  /* Begin parsing the input file */
  while(fgets(inp_line, MAX_INP_LINE, inp_fd))
  {
    if(inp_line[0] != '#')
    {
      strcpy(inp_copy, inp_line);
      clean_string(inp_line, " \t");
      cptr = strtok(inp_line, "\t=");
      if(token_compare(cptr, "input exodusii file"))
      {
        /* The input ExodusII file name */
        if(strlen(exoII_inp_file) == 0)
        {
          cptr = strtok(NULL, "\t=");
          strip_string(cptr, " \t\n");
          strcpy(exoII_inp_file, cptr);
        }
      }
      else if(token_compare(cptr, "output visualization file"))
      {
        if(problem->vis_out < 0)
        {
          /* Output a visualization file */
          cptr = strtok(NULL, "\t=");
          strip_string(cptr, " \t\n");
          if(strcasecmp(cptr, "yes") == 0 || strcasecmp(cptr, "true") == 0)
            problem->vis_out = 1;
          else if(strcasecmp(cptr, "no") == 0 ||
                  strcasecmp(cptr, "false") == 0)
          {
            problem->vis_out = 0;
          }
          else
          {
            iret = sscanf(cptr, "%d", &problem->vis_out);
            if(iret != 1)
            {
              Gen_Error(1, "warning: unknown visualization output flag");
              problem->vis_out = 0;
            }
            else
            {
              if(problem->vis_out != 1)
                problem->vis_out = 0;
            }
          }
        }
      }
      else if(token_compare(cptr, "output nemesisi file"))
      {
        /* The NemesisI output file name */
        if(strlen(nemI_out_file) == 0)
        {
          cptr = strtok(NULL, "\t=");
          strip_string(cptr, " \t\n");
          strcpy(nemI_out_file, cptr);
        }
      }
      else if(token_compare(cptr, "decomposition method"))
      {
        /* The method to use for decomposing the graph */
        if(lb->type < 0 || lb->refine < 0 || lb->num_sects < 0)
        {

          /* Search to the first null character */
          cptr = strchr(cptr, '\0');
          cptr++;
          strip_string(cptr, " \t\n=");
          cptr = strtok(cptr, ",");
          while(cptr != NULL)
          {
            strip_string(cptr, " \t\n");
            string_to_lower(cptr, '\0');
            if(strcmp(cptr, "multikl") == 0)
            {
              if(lb->type < 0)
                lb->type = MULTIKL;
            }
            else if(strcmp(cptr, "spectral") == 0)
            {
              if(lb->type < 0)
                lb->type = SPECTRAL;
            }
            else if(strcmp(cptr, "scattered") == 0)
            {
              if(lb->type < 0)
                lb->type = SCATTERED;
            }
            else if(strcmp(cptr, "linear") == 0)
            {
              if(lb->type < 0)
                lb->type = LINEAR;
            }
            else if(strcmp(cptr, "inertial") == 0)
            {
              if(lb->type < 0)
                lb->type = INERTIAL;
            }
            else if(strcmp(cptr, "zpinch") == 0)
            {
              if(lb->type < 0)
                lb->type = ZPINCH;
            }
            else if(strcmp(cptr, "brick") == 0)
            {
              if(lb->type < 0)
                lb->type = BRICK;
            }
            else if(strcmp(cptr, "rcb") == 0)
            {
              if(lb->type < 0)
                lb->type = ZOLTAN_RCB;
            }
            else if(strcmp(cptr, "rib") == 0)
            {
              if(lb->type < 0)
                lb->type = ZOLTAN_RIB;
            }
            else if(strcmp(cptr, "hsfc") == 0)
            {
              if(lb->type < 0)
                lb->type = ZOLTAN_HSFC;
            }
            else if(strcmp(cptr, "random") == 0)
            {
              if(lb->type < 0)
                lb->type = RANDOM;
            }
            else if(strstr(cptr, "infile"))
            {
              if(lb->type < 0)
              {
                lb->type = INFILE;
                cptr2 = strchr(cptr, '=');
                if(cptr2 == NULL)
                {
                  Gen_Error(0, "fatal: need to specify a value with infile");
                  return 0;
                }

                cptr2++;
                iret = sscanf(cptr2, "%s", lb->file);
                if(iret != 1)
                {
                  Gen_Error(0, "fatal: invalid value for infile");
                  return 0;
                }
              }
            }
            else if(strcmp(cptr, "kl") == 0)
            {
              if(lb->refine < 0)
                lb->refine = KL_REFINE;
            }
            else if(strcmp(cptr, "none") == 0)
            {
              if(lb->refine < 0)
                lb->refine = NONE;
            }
            else if(strstr(cptr, "num_sects"))
            {
              if(lb->num_sects < 0)
              {
                cptr2 = strchr(cptr, '=');
                if(cptr2 == NULL)
                {
                  Gen_Error(0,
                            "fatal: need to specify a value with num_sects");
                  return 0;
                }

                cptr2++;
                iret = sscanf(cptr2, "%d", &(lb->num_sects));
                if(iret != 1)
                {
                  Gen_Error(0, "fatal: invalid value for num_sects");
                  return 0;
                }
              }
            }
            else if(strcmp(cptr, "cnctd_dom") == 0)
            {
              if(lb->cnctd_dom < 0)
                lb->cnctd_dom = ELB_TRUE;
            }
            else if(strstr(cptr, "outfile"))
            {
              if(lb->outfile < 0)
              {
                lb->outfile = ELB_TRUE;
                cptr2 = strchr(cptr, '=');
                if(cptr2 == NULL)
                {
                  Gen_Error(0, "fatal: need to specify a value with outfile");
                  return 0;
                }

                cptr2++;
                iret = sscanf(cptr2, "%s", lb->file);
                if(iret != 1)
                {
                  Gen_Error(0, "fatal: invalid value for outfile");
                  return 0;
                }
              }
            }
            else
            {
              sprintf(ctemp,
                      "fatal: unknown LB method \"%s\" specified in command"
                      " file",
                      cptr);
              Gen_Error(0, ctemp);
            }
            cptr = strtok(NULL, ",");
          }
        }
      }
      else if(token_compare(cptr, "solver specifications"))
      {
        /* Solver specs */

        /* Search to the first null character */
        cptr = strchr(cptr, '\0');
        cptr++;
        strip_string(cptr, " \t\n=");
        cptr = strtok(cptr, ",");

        /* Loop until all the suboptions have been specified */
        while(cptr)
        {
          strip_string(cptr, " \t\n");
          string_to_lower(cptr, '\0');

          /* Check to see if this is the "tolerance" suboption */
          if(strstr(cptr, "tolerance"))
          {
            if(solver->tolerance < 0.0)
            {
              cptr2 = strchr(cptr, '=');
              if(cptr2 == NULL)
              {
                Gen_Error(0,
                          "fatal: tolerance specification requires a value");
                return 0;
              }

              cptr2++;
              iret = sscanf(cptr2, "%le", &(solver->tolerance));
              if(iret != 1)
              {
                Gen_Error(0, "fatal: invalid value for tolerance");
                return 0;
              }
            }
          }
          else if(strcmp(cptr, "use_rqi") == 0)
          {
            if(solver->rqi_flag == USE_RQI)
              solver->rqi_flag = -1;
            else
              solver->rqi_flag = USE_RQI;
          }
          else if(strstr(cptr, "vmax"))
          {
            if(solver->vmax < 0)
            {
              cptr2 = strchr(cptr, '=');
              if(cptr2 == NULL)
              {
                Gen_Error(0, "fatal: vmax must have a value");
                return 0;
              }

              cptr2++;
              iret = sscanf(cptr2, "%d", &(solver->vmax));
              if(iret != 1)
              {
                Gen_Error(0, "fatal: invalid value for vmax");
                return 0;
              }
            }
          }
          else
          {
            sprintf(ctemp, "warning: unknown solver suboption %s",
                    cptr);
            Gen_Error(1, ctemp);
          }

          cptr = strtok(NULL, ",");
        }
      }
      else if(token_compare(cptr, "graph type"))
      {
        if(problem->type < 0)
        {
          cptr = strtok(NULL, "\t=");
          strip_string(cptr, " \t\n");
          string_to_lower(cptr, '\0');
          if(strcmp(cptr, "nodal") == 0)
            problem->type = NODAL;
          else if(strcmp(cptr, "node") == 0)
            problem->type = NODAL;
          else if(strcmp(cptr, "elemental") == 0)
            problem->type = ELEMENTAL;
          else if(strcmp(cptr, "element") == 0)
            problem->type = ELEMENTAL;
          else
          {
            Gen_Error(0, "fatal: unknown graph type specified");
            return 0;
          }

        }
      }
      else if(token_compare(cptr, "machine description"))
      {
        /* Machine specs */
        if(machine->num_dims < 0)
        {
          /* Search to first null character */
          cptr = strchr(cptr, '\0');
          cptr++;
          strip_string(cptr, " \t\n=");

          /* Search to equal sign */
          cptr2 = strchr(cptr, '=');
          if(cptr2 == NULL)
          {
            Gen_Error(0, "fatal: machine must have a dimension specified");
            return 0;
          }

          *cptr2 = '\0';

          /* Find out the machine type */
          strip_string(cptr, " \t\n");
          if(strcasecmp(cptr, "mesh") == 0)
          {
            machine->type = MESH;
            max_dim = 3;
          }
          else if(strcasecmp(cptr, "hcube") == 0 ||
                  strcasecmp(cptr, "hypercube") == 0)
          {
            machine->type = HCUBE;
            max_dim = 1;
          }
          else if(strcasecmp(cptr, "cluster") == 0)
          {
            /* now need to find what each box consists of */
            cptr = cptr2 + 1;
            cptr2 = strpbrk(cptr, "mMhH");
            if (*cptr2 == 'm' || *cptr2 == 'M') {
              machine->type = MESH;
              max_dim = 3;
            }
            else if (*cptr2 == 'h' || *cptr2 == 'H') {
              machine->type = HCUBE;
              max_dim = 1;
            }
            else
            {
              Gen_Error(0, "fatal: unknown type specified with cluster");
              return 0;
            }
            /* blank out character and move cptr to next char */
            *cptr2 = '\0';

            /* get the number of boxes from value */
            iret = sscanf(cptr, "%d", &(machine->num_boxes));
            if(iret <= 0 || machine->num_boxes <= 0)
            {
              Gen_Error(0, "fatal: invalid number of boxes");
              return 0;
            }
          }
          else
          {
            Gen_Error(0, "fatal: unknown machine type specified");
            return 0;
          }

          machine->num_dims = 0;
          cptr = cptr2 + 1;
          cptr = strtok(cptr, "xX");
          while(cptr)
          {
            iret = sscanf(cptr, "%d", &(machine->dim[machine->num_dims]));
            if(iret != 1)
            {
              Gen_Error(0, "fatal: invalid dimension specified for machine");
              return 0;
            }

            machine->num_dims++;
            cptr = strtok(NULL, "xX");

            /* Check how many dimensions there are */
            if(machine->num_dims == max_dim && cptr != NULL)
            {
              Gen_Error(0, "fatal: maximum number of dimensions exceeded");
              return 0;
            }
          }
        }
      }
      else if(token_compare(cptr, "weighting specifications"))
      {
        /* Parameters for weighting the graph */
        if(weight->type < 0)
        {
          cptr = strchr(cptr, '\0');
          cptr++;
          strip_string(cptr," \t\n=");
          cptr = strtok(cptr, ",");

          while(cptr != NULL)
          {
            strip_string(cptr, " \t\n");
            string_to_lower(cptr, '\0');
            if(strstr(cptr, "read"))
            {
              cptr2 = strchr(cptr, '=');
              if(cptr2 == NULL)
              {
                Gen_Error(0, "fatal: must specify file name with \"read\"");
                return 0;
              }
              cptr2++;
              if(strlen(cptr2) == 0)
              {
                Gen_Error(0, "fatal: invalid file name with \"read\"");
                return 0;
              }
              strcpy(weight->exo_filename, cptr2);
              if (weight->type < 0)                weight->type = READ_EXO;
              else if (!(weight->type & READ_EXO)) weight->type += READ_EXO;

              /* check if the read is after an element block weight */
              if (weight->type & EL_BLK) weight->ow_read = 0;
            }
            else if(strstr(cptr, "var_name"))
            {
              cptr2 = strchr(cptr, '=');
              if(cptr2 == NULL)
              {
                Gen_Error(0, "fatal: must specify a name with \"var_name\"");
                return 0;
              }
              cptr2++;
              if(strlen(cptr2) == 0)
              {
                Gen_Error(0, "fatal: invalid variable name specified with"
                          " \"var_name\"");
                return 0;
              }
              strcpy(weight->exo_varname, cptr2);
            }
            else if(strstr(cptr, "var_index"))
            {
              cptr2 = strchr(cptr, '=');
              if(cptr2 == NULL)
              {
                Gen_Error(0, "fatal: must specify a value with \"var_index\"");
                return 0;
              }
              cptr2++;

              iret = sscanf(cptr2, "%d", &(weight->exo_vindx));
              if(iret != 1)
              {
                Gen_Error(0, "fatal: invalid value with \"var_index\"");
                return 0;
              }

            }
            else if(strstr(cptr, "time_index"))
            {
              cptr2 = strchr(cptr, '=');
              if(cptr2 == NULL)
              {
                Gen_Error(0,
                          "fatal: must specify a value with \"time_index\"");
                return 0;
              }
              cptr2++;

              iret = sscanf(cptr2, "%d", &(weight->exo_tindx));
              if(iret != 1)
              {
                Gen_Error(0, "fatal: invalid value with \"time_index\"");
                return 0;
              }
            }
            else if(strstr(cptr, "eb"))
            {
              cptr2 = strchr(cptr, '=');
              if(cptr2 == NULL)
              {
                Gen_Error(0,
                          "fatal: must specify a value with \"eb\"");
                return 0;
              }
              cptr2++;

              el_blk = -1;
              wgt = -1;

              iret = sscanf(cptr2, "%d:%d", &el_blk, &wgt);
              if(iret != 2)
              {
                Gen_Error(0, "fatal: invalid value with \"eb\"");
                return 0;
              }
              if(el_blk <= 0)
              {
                printf(ctemp, "invalid element block, %d", el_blk);
                Gen_Error(0, ctemp);
                return 0;
              }
              if(wgt < 1)
              {
                printf(ctemp, "invalid weight, %d", wgt);
                Gen_Error(0, ctemp);
                return 0;
              }

              if (weight->type < 0)              weight->type = EL_BLK;
              else if (!(weight->type & EL_BLK)) weight->type += EL_BLK;

              /* check if this is the first element block weight */
              weight->num_ebw++;
              if (weight->num_ebw == 1) {

                /* since this is the first time through allocate array */
                weight->elemblk = (int *) malloc(sizeof(int));
                weight->elemblk_wgt = (int *) malloc(sizeof(int));
                if(!weight->elemblk || !weight->elemblk_wgt)
                {
                  Gen_Error(0, "fatal: insufficient memory");
                  return 0;
                }
              }
              else { /* just realloc arrays */
                weight->elemblk = (int *) realloc(weight->elemblk,
                                              weight->num_ebw*sizeof(int));
                weight->elemblk_wgt = (int *) realloc(weight->elemblk_wgt,
                                                  weight->num_ebw*sizeof(int));
                if(!weight->elemblk || !weight->elemblk_wgt)
                {
                  Gen_Error(0, "fatal: insufficient memory");
                  return 0;
                }
              }
              /* now set the values */
              weight->elemblk[weight->num_ebw-1] = el_blk;
              weight->elemblk_wgt[weight->num_ebw-1] = wgt;

              /* check if the elem block weight needs to overwrite the read */
              if (weight->type & READ_EXO) weight->ow_read = 1;

            }
            else if(strstr(cptr, "edges"))
            {
              if (weight->type < 0)                weight->type = EDGE_WGT;
              else if (!(weight->type & EDGE_WGT)) weight->type += EDGE_WGT;
            }
            else
            {
              sprintf(ctemp, "fatal: unknown suboption \"%s\" specified",
                      cptr);
              Gen_Error(0, ctemp);
              return 0;
            }

            cptr = strtok(NULL, ",");
          }

        } /* End "if(weight->type < 0)" */

      }
      else if(token_compare(cptr, "misc options"))
      {
        /* Misc Options */

        /* Search to the first null character */
        cptr = strchr(cptr, '\0');
        cptr++;
        /*
         * Need to handle case where users have put comma's in
         * the group descriptor. This will mess up getting the
         * tokens using strtok(). So, search for commas between
         * the beginning delimiter, "{", and the end delimiter,
         * "}", and change them to blank spaces.
         */
        cptr2 = strchr(cptr, '{');
        if (cptr2 != NULL) {
          ilen = strlen(cptr2);
          for (i = 0; i < ilen; i++) {
            if (*cptr2 == '}') break;
            if (*cptr2 == ',') *cptr2 = ' ';
            cptr2++;
          }
        }

        strip_string(cptr, " \t\n=");
        cptr = strtok(cptr, ",");

        /* Loop until all the suboptions have been specified */
        while(cptr != NULL)
        {
          strip_string(cptr, " \t\n");
          string_to_lower(cptr, '\0');

          /* Check to see if the side id error checks need to be skipped */
          if(strstr(cptr, "checks_off"))
          {
            if(problem->skip_checks < 0)
              problem->skip_checks = 1;
          }
          /* Check to see if using face definition of adjacency */
          else if(strstr(cptr, "face_adj"))
          {
            if(problem->face_adj < 0)
              problem->face_adj = 1;
          }
          /* Check to see if looking for global mechanisms */
          else if(strstr(cptr, "global_mech"))
          {
            if(problem->global_mech < 0)
              problem->global_mech = 1;
          }
          /* Check to see if looking for introduced mechanisms */
          else if(strstr(cptr, "local_mech"))
          {
            if(problem->local_mech < 0)
              problem->local_mech = 1;
          }
          /* Check to see if looking for connected domains */
          else if(strstr(cptr, "find_cnt_domains"))
          {
            if(problem->find_cnt_domains < 0)
              problem->find_cnt_domains = 1;
          }
          /* Check to see if user wants to add processors to take care of
           * introduced mechanisms 
           */
          else if(strstr(cptr, "mech_add_procs"))
          {
            if(problem->mech_add_procs < 0)
              problem->mech_add_procs = 1;
          }
          /* Check to see if user wants to add processors to take care of
           * introduced mechanisms 
           */
          else if(strstr(cptr, "dsd_add_procs"))
          {
            if(problem->dsd_add_procs < 0)
              problem->dsd_add_procs = 1;
          }
          /* Check for treating spheres as concentrated masses*/
          else if(strstr(cptr, "no_sph"))
          {
            if(problem->no_sph < 0)
              problem->no_sph = 1;
          }
          /* Check for group designation sub-option */
          else if (strstr(cptr, "groups"))
          {
            /* "{" defines the beginning of the group designator */
            cptr2 = strchr(cptr, '{');
            if (cptr2== NULL) {
              Gen_Error(0, "fatal: group start designator \"}\" not found");
              return 0;
            }
            cptr2++;
            /* allocate space to hold the group designator */
            problem->groups = malloc(strlen(cptr2) + 1);
            strcpy(problem->groups, cptr2);
            /* get rid of ending bracket */
            cptr2 = strchr(problem->groups, '}');
            *cptr2 = '\0';
          }
          else
          {
            sprintf(ctemp, "warning: unknown miscellaneous suboption %s",
                    cptr);
            Gen_Error(1, ctemp);
          }
          cptr = strtok(NULL, ",");
        }
      }
      else
      {
        /* Generate an error, but continue reading for an unknown key */
        strip_string(inp_copy, " #\t");
        if(strlen(inp_copy) > 5)
        {
          sprintf(ctemp,
                  "warning: don't know how to process line: \n%s\nin command"
                  " file, ignored",
                  inp_copy);
          Gen_Error(1, ctemp);
        }
      }

    } /* End "if(inp_line[0] != '#')" */

  } /* End "while(fgets(inp_line, MAX_INP_LINE, inp_fd))" */

  /* Close the input ascii command file */
  fclose(inp_fd);

  return 1;

} /*------------End read_cmd_file()---------------*/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* This function performs error checks on the user input.
 *---------------------------------------------------------------------------*/
int check_inp_specs(
  char exoII_inp_file[],
  char nemI_out_file[],
  MACHINE_PTR machine,
  LB_INFO_PTR lb,
  PROB_INFO_PTR prob,
  SOLVE_INFO_PTR solver,
  WEIGHT_INFO_PTR weight
  )
{
  char  ctemp[MAX_ERR_MSG+1];
  char  ctemp2[MAX_FNL+1];
  char *cptr;
  char **var_names;
  int   cnt, exoid, cpu_ws=0, io_ws=0, nvars, tmp_vindx=0;
  int   ntimes;
  float version;
  FILE *inp_fd;
  int   exid_inp, icpu_ws=0, iio_ws=0;
  float vers;

  /* Check that an input ExodusII file name was specified */
  if(strlen(exoII_inp_file) <= 0)
  {
    Gen_Error(0, "fatal: no input ExodusII file specified");
    return 0;
  }

  /* Check for the existence and readability of the input file */
  if((exid_inp=ex_open(exoII_inp_file, EX_READ, &icpu_ws, &iio_ws, &vers)) < 0)
  {
    sprintf(ctemp, "fatal: unable to open input ExodusII file %s",
            exoII_inp_file);
    Gen_Error(0, ctemp);
    return 0;
  }
  ex_close(exid_inp);

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                       Check the machine specification                     */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  if(machine->type != MESH && machine->type != HCUBE)
  {
    Gen_Error(0, "fatal: machine type not properly set");
    return 0;
  }
  if(machine->type == HCUBE && machine->num_dims != 1)
  {
    Gen_Error(0, "fatal: improper number of dimension for a hypercube, only"
              " 1 allowed");
    return 0;
  }
  if(machine->type == MESH && machine->num_dims > 3)
  {
    Gen_Error(0, "fatal: maximum of 3 dimensions for a mesh exceeded");
    return 0;
  }

  /* non-cluster machines have only one box */
  if (machine->num_boxes < 0) machine->num_boxes = 1;

  /* Find out the number of processors */
  if(machine->type == HCUBE)
    machine->procs_per_box = 1 << machine->dim[0];
  else
  {
    machine->procs_per_box = machine->dim[0];
    for(cnt=1; cnt < machine->num_dims; cnt++)
      machine->procs_per_box *= machine->dim[cnt];
  }

  /* now calculate the total number of processors */
  machine->num_procs = machine->num_boxes * machine->procs_per_box;

  /*
   * currently, do not allow groups and clusters since the
   * loops to chaco get a bit too confusing
   */
  if (machine->num_boxes > 1 && prob->groups != NULL)
  {
    Gen_Error(0, "fatal: groups cannot be designated for a cluster machine");
    return 0;
  }

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                      Check the problem specifications                     */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  if(prob->type != ELEMENTAL && prob->type != NODAL)
  {
    Gen_Error(0, "fatal: unknown problem type specified");
    return 0;
  }

  if (prob->skip_checks < 0) prob->skip_checks = 0;

  if (prob->face_adj < 0) prob->face_adj = 0;

  /*
   * using face definition of adjacencies only makes sense
   * with an elemental decomposition
   */
  if (prob->type != ELEMENTAL && prob->face_adj) {
    Gen_Error(1, "warning: can only use face definition of");
    Gen_Error(1, "warning: adjacency with elemental decomposition");
    Gen_Error(1, "warning: face definition turned off");
    prob->face_adj = 0;
  }

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                      Check the load balance parameters                    */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  if(lb->type != MULTIKL && lb->type != SPECTRAL && lb->type != INERTIAL &&
     lb->type != LINEAR  && lb->type != RANDOM   && lb->type != SCATTERED &&
     lb->type != INFILE  && lb->type != ZPINCH   && lb->type != BRICK &&
     lb->type != ZOLTAN_RCB && lb->type != ZOLTAN_RIB && 
     lb->type != ZOLTAN_HSFC)
  {
    Gen_Error(0, "fatal: unknown load balance type requested");
    return 0;
  }

  if(lb->type == MULTIKL)
    lb->refine = KL_REFINE;

  if(lb->refine != KL_REFINE && lb->refine != NO_REFINE)
    lb->refine = NO_REFINE;	/* Default if not specified */

  if(lb->num_sects <= 0)
      lb->num_sects = 1;	/* Default if not specified */

  if(lb->cnctd_dom < 0)
    lb->cnctd_dom = 0;
  else if(!prob->face_adj) {
    Gen_Error(1, "warning: can only set connected domain");
    Gen_Error(1, "warning: when using face definition of adjacency");
    Gen_Error(1, "warning: connected domain turned off");
    lb->cnctd_dom = 0;
  }

  if(lb->outfile < 0)
    lb->outfile = ELB_FALSE;

  if(lb->type == INFILE)
  {
    if(lb->outfile)
    {
      Gen_Error(0, "fatal: both infile and outfile cannot be specified");
      return 0;
    }

    if (lb->refine != NO_REFINE) {
      Gen_Error(1, "warning: no refinement can be specified with infile");
      return 0;
    }
  }


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                   Check the eigensolver parameters                        */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  if(lb->type == SPECTRAL || lb->type == MULTIKL)
  {
    if(solver->tolerance < 0.0)
    {
      Gen_Error(1, "warning: using default value for eigensolver"
                " tolerance");
      solver->tolerance = 1.0e-3;
    }

    if(solver->rqi_flag < 0)
      solver->rqi_flag = 0;

    if(solver->vmax < 0)
    {
      Gen_Error(1, "warning: no value for vmax specified,"
                " using default");
      solver->vmax = 200;
    }

  }

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                     Check the output file name                            */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  if(strlen(nemI_out_file) <= 0)
  {
    /*
     * Generate the file name from the input file name and the requested
     * load balance method.
     */
    switch(machine->type)
    {
    case MESH:
      sprintf(ctemp2, "-m%d-", machine->num_procs);
      break;

    case HCUBE:
      sprintf(ctemp2, "-h%d-", machine->num_procs);
      break;
    }

    switch(lb->type)
    {
    case MULTIKL:
    case SPECTRAL:
      if(lb->num_sects == 1)
        strcat(ctemp2, "b");
      else if(lb->num_sects == 2)
        strcat(ctemp2, "q");
      else if(lb->num_sects == 3)
        strcat(ctemp2, "o");

      break;

    case INERTIAL:
      strcat(ctemp2, "i");
      break;

    case ZPINCH:
      strcat(ctemp2, "z");
      break;

    case BRICK:
      strcat(ctemp2, "x");
      break;

    case ZOLTAN_RCB:
    case ZOLTAN_RIB:
    case ZOLTAN_HSFC:
      strcat(ctemp2, "Z");
      break;

    case SCATTERED:
      strcat(ctemp2, "s");
      break;

    case RANDOM:
      strcat(ctemp2, "r");
      break;

    case LINEAR:
      strcat(ctemp2, "l");
      break;
    }

    if(lb->refine == KL_REFINE)
      strcat(ctemp2, "KL");

    /* Generate the complete file name */
    strcpy(nemI_out_file, exoII_inp_file);
    cptr = strrchr(nemI_out_file, '.');
    *cptr = '\0';
    strcat(nemI_out_file, ctemp2);
    strcat(nemI_out_file, ".nemI");

  } /* End "if(strlen(nemI_out_file) <= 0)" */

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                 Check the weighting specifications                        */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  if(strlen(weight->exo_filename) > 0)
  {

    /* Check for the existence of the specified file */
    if(!(inp_fd=fopen(weight->exo_filename, "r")))
    {
      sprintf(ctemp, "fatal: unable to open ExodusII weighting file %s",
              weight->exo_filename);
      Gen_Error(0, ctemp);
      return 0;
    }
    fclose(inp_fd);

    /* Check that a variable name and/or index was specified. */
    if(strlen(weight->exo_varname) == 0 && weight->exo_vindx <= 0)
    {
      Gen_Error(0, "fatal: must specify an index and/or a name for weighting"
                " variable");
      return 0;
    }

    /*
     * If a variable name and or index was specified then open the ExodusII
     * file and compare the specified name against what exists in the file.
     */
    if((exoid=ex_open(weight->exo_filename, EX_READ, &cpu_ws,
                      &io_ws, &version)) < 0)
    {
      sprintf(ctemp, "fatal: failed to open ExodusII weighting file %s",
              weight->exo_filename);
      Gen_Error(0, ctemp);
      return 0;
    }

    if(ex_inquire(exoid, EX_INQ_TIME, &ntimes, NULL, NULL) < 0)
    {
      Gen_Error(0, "fatal: can't get number of time steps in weighting file");
      ex_close(exoid);
      return 0;
    }

    /* Check the time index */
    if(weight->exo_tindx <= 0)
      weight->exo_tindx = 1;	/* Defaults to 1 */

    if(weight->exo_tindx > ntimes)
    {
      sprintf(ctemp,
              "fatal: requested time index %d not available in weighting file",
              weight->exo_tindx);
      Gen_Error(0, ctemp);
      return 0;
    }

    if(prob->type == NODAL)
      strcpy(ctemp, "n");
    else
      strcpy(ctemp, "e");

    /*
     * First check that there are variables of the requested type in the
     * specified ExodusII file.
     */
    if(ex_get_var_param(exoid, ctemp, &nvars) < 0)
    {
      Gen_Error(0, "fatal: unable to get variable params from ExodusII"
                " weighting file");
      return 0;
    }
    if(nvars <= 0)
    {
      Gen_Error(0, "fatal: no variables found in the ExodusII weighting file");
      return 0;
    }

    /* Read the variable names from the requested file */
    var_names = malloc(nvars*sizeof(char *));
    if(!var_names)
    {
      Gen_Error(0, "fatal: insufficient memory");
      return 0;
    }

    {
      int max_name_length = ex_inquire_int(exoid, EX_INQ_DB_MAX_USED_NAME_LENGTH);
      ex_set_max_name_length(exoid, max_name_length);
      
      for(cnt=0; cnt < nvars; cnt++) {
	var_names[cnt] = malloc((max_name_length+1)*sizeof(char));
	if(!var_names[cnt]) {
	  Gen_Error(0, "fatal: insufficient memory");
	  return 0;
	}
      }
    }

    if(ex_get_var_names(exoid, ctemp, nvars, var_names) < 0) {
      Gen_Error(0, "fatal: unable to obtain variable names for weighting");
      return 0;
    }

    /*
     * If there was a variable name specified but no index then get the
     * index. If a variable name AND an index were specified then make
     * sure they match.
     */
    if(strlen(weight->exo_varname) > 0) {
      for(cnt=0; cnt < nvars; cnt++) {
        if(strcmp(var_names[cnt], weight->exo_varname) == 0) {
          tmp_vindx = cnt+1;

          break; /* out of "for" loop */
        }
      }

      if(weight->exo_vindx <= 0)
        weight->exo_vindx = tmp_vindx;
      else if(weight->exo_vindx != tmp_vindx) {
        Gen_Error(0, "fatal: requested weight variable index doesn't match"
                  " the variable name in the ExodusII file");
        return 0;
      }
    }

    /* Free up memory */
    if(nvars > 0) {
      for(cnt=0; cnt < nvars; cnt++)
        free(var_names[cnt]);

      free(var_names);
    }

    /*
     * If there is still no valid index then the variable name does
     * not exist in the specified file.
     */
    if(weight->exo_vindx <= 0) {
      sprintf(ctemp,
              "fatal: requested weighting variable %s not found in ExodusII"
              " file", weight->exo_varname);
      Gen_Error(0, ctemp);
      return 0;
    }

    if(nvars < weight->exo_vindx) {
      Gen_Error(0, "fatal: requested variable index is larger than number in"
                " ExodusII weighting file");
      return 0;
    }

    if(weight->exo_vindx <= 0) {
      sprintf(ctemp, "fatal: variable index must be in the range [1,%d]",
              nvars);
      Gen_Error(0, ctemp);
      return 0;
    }

    ex_close(exoid);

  } /* End "if(strlen(weight->exo_filename) > 0)" */

  if ((weight->type & EL_BLK) && (weight->ow_read)) {
    if (weight->num_ebw > 1) {
      /* start by sorting the two arrays by the element block number */
      sort2_int_int(weight->num_ebw, weight->elemblk, weight->elemblk_wgt);

      /* now loop through, and make sure that we don't have multiple values */
      for (cnt=1; cnt < weight->num_ebw; cnt++)
        if (weight->elemblk[cnt] == weight->elemblk[cnt-1]) {
          sprintf(ctemp, "warning: multiple weight specified for block %d",
                  weight->elemblk[cnt]);
          Gen_Error(1, ctemp);
        }
    }

  } /* if ((weight->type & EL_BLK) && (weight->ow_read)) */

  return 1;

} /*-------------------End check_inp_specs()-----------------*/
