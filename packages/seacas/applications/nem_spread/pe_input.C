/*
 * Copyright(C) 1999-2020, 2023 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*--------------------------------------------------------------------------*/
/* Purpose: Determine file types for command files and read in the parallel */
/*          ExodusII command file.                                          */
/*--------------------------------------------------------------------------*/

#include "copy_string_cpp.h"
#include "fmt/ostream.h"
#include "nem_spread.h"
#include "ps_pario_const.h" // for PIO_Info, Parallel_IO, etc
#include "rf_allo.h"        // for array_alloc
#include "rf_io_const.h"    // for ExoFile, etc
#include "scopeguard.h"
#include <cstdio>  // for nullptr, stderr, etc
#include <cstdlib> // for exit, realloc
#include <cstring> // for strtok, strchr, strstr, etc

#define TLIST_CNT 5

constexpr int MAX_INPUT_STR_LN = 4096; /* maximum string length for read_string()  */

/*****************************************************************************/
int read_mesh_file_name(const char *filename)
{
  /* local declarations */
  FILE *file_cmd = nullptr;
  char  inp_line[MAX_INPUT_STR_LN + 1];
  char  inp_copy[MAX_INPUT_STR_LN + 1];

  /* Open the file */
  if ((file_cmd = fopen(filename, "r")) == nullptr) {
    return -1;
  }
  ON_BLOCK_EXIT(fclose, file_cmd);

  /* Begin parsing the input file */
  while (fgets(inp_line, MAX_INPUT_STR_LN, file_cmd) != nullptr) {
    /* skip any line that is a comment */
    if ((inp_line[0] != '#') && (inp_line[0] != '\n')) {

      copy_string(inp_copy, inp_line);
      clean_string(inp_line, " \t");
      char *cptr = strtok(inp_line, "\t=");
      /****** The input ExodusII file name ******/
      if (token_compare(cptr, "input fem file") != 0) {
        if (ExoFile.empty()) {
          cptr = strtok(nullptr, "\t=");
          strip_string(cptr, " \t\n");
          ExoFile = std::string(cptr);
          break;
        }
      }
    }
  }
  return 0;
}

/*****************************************************************************/

template <typename T, typename INT>
int read_pexoII_info(NemSpread<T, INT> &spreader, const char *filename)

/*
 *          This function reads the ASCII parallel-exodus command file.
 *
 *   Input
 *   -----
 *   filename - The name of the command file.
 */
{
  /* local declarations */
  FILE *file_cmd = nullptr;
  char  inp_line[MAX_INPUT_STR_LN + 1];
  char  inp_copy[MAX_INPUT_STR_LN + 1];
  char *cptr;
  char *cptr2;
  char *cptr3;
  int   i;
  int   icnt;
  int   tlist_alloc;

  /***************************** BEGIN EXECUTION ******************************/

  /* Open the file */
  if ((file_cmd = fopen(filename, "r")) == nullptr) {
    return -1;
  }
  ON_BLOCK_EXIT(fclose, file_cmd);

  /* Begin parsing the input file */
  while (fgets(inp_line, MAX_INPUT_STR_LN, file_cmd)) {
    /* skip any line that is a comment */
    if ((inp_line[0] != '#') && (inp_line[0] != '\n')) {

      copy_string(inp_copy, inp_line);
      clean_string(inp_line, " \t");
      cptr = strtok(inp_line, "\t=");
      /****** The input ExodusII file name ******/
      if (token_compare(cptr, "input fem file")) {
        if (ExoFile.empty()) {
          cptr = strtok(nullptr, "\t=");
          strip_string(cptr, " \t\n");
          ExoFile = std::string(cptr);
        }
      }
      /****** The input NemesisI load balance file name ******/
      else if (token_compare(cptr, "lb file")) {
        if (Exo_LB_File.empty()) {
          cptr = strtok(nullptr, "\t=");
          strip_string(cptr, " \t\n");
          Exo_LB_File = std::string(cptr);
        }
      }
      /****** The scalar results ExodusII file name ******/
      else if (token_compare(cptr, "scalar results fem file")) {
        if (Exo_Res_File.empty()) {
          cptr = strtok(nullptr, "\t=");
          strip_string(cptr, " \t\n");
          Exo_Res_File = std::string(cptr);
        }
      }
      /****** The parallel results ExodusII file name ******/
      else if (token_compare(cptr, "parallel results file base name")) {
        if (Output_File_Base_Name.empty()) {
          cptr = strtok(nullptr, "\t=");
          strip_string(cptr, " \t\n");
          Output_File_Base_Name = std::string(cptr);
        }
      }
      /****** The Number of Processors ******/
      else if (token_compare(cptr, "number of processors")) {
        if (spreader.Proc_Info[0] < 0) {
          cptr = strtok(nullptr, "\t=");
          strip_string(cptr, " \t\n");
          if (sscanf(cptr, "%d", &(spreader.Proc_Info[0])) != 1) {
            fmt::print(stderr,
                       "{}: ERROR, can\'t interpret int for number of"
                       " Processors.\n",
                       __func__);
            exit(1);
          }
        }
      }
      /****** The File extension to use for spread files ******/
      else if (token_compare(cptr, "file extension for spread files")) {
        cptr = strtok(nullptr, "\t=");
        strip_string(cptr, " \t\n");
        PIO_Info.Exo_Extension = std::string(cptr);
      }

      /****** Is There a Scalar Mesh File to Use ******/
      else if (token_compare(cptr, "use scalar mesh file")) {
        if (Gen_Flag < 0) {
          cptr = strtok(nullptr, "\t=");
          strip_string(cptr, " \t\n");
          if (Gen_Flag < 0) {
            if (token_compare(cptr, "yes")) {
              Gen_Flag = 1;
            }
            else {
              Gen_Flag = 0;
            }
          }
        }
      }
      /****** The Debug reporting level ******/
      else if (token_compare(cptr, "debug")) {
        if (Debug_Flag < 0) {
          cptr = strtok(nullptr, "\t=");
          strip_string(cptr, " \t\n");
          if (sscanf(cptr, "%d", &Debug_Flag) != 1) {
            fmt::print(stderr, "{}: ERROR, can\'t interpret int for Debug_Flag\n", __func__);
            exit(1);
          }
        }
      }
      /****** Restart Time List ******/
      else if (token_compare(cptr, "restart info")) {

        cptr = strchr(cptr, '\0');
        cptr++;

        /*
         * If there is a list, then need to change the comma's in
         * the list to blank spaces so that the strtok below works
         * correctly. So, search for commas between the beginning
         * delimiter, "{", and the end delimiter, "}", and change
         * them to blank spaces.
         */
        cptr2 = strchr(cptr, '{');
        if (cptr2 != nullptr) {
          icnt = strlen(cptr2);
          for (i = 0; i < icnt; i++) {
            if (*cptr2 == '}') {
              break;
            }
            if (*cptr2 == ',') {
              *cptr2 = ' ';
            }
            cptr2++;
          }
        }

        strip_string(cptr, " \t\n=");
        cptr = strtok(cptr, ",");

        /* Loop until all the suboptions have been specified */
        while (cptr != nullptr) {
          strip_string(cptr, " \t\n");
          string_to_lower(cptr, '\0');

          /* check if the user wants to specifically turn this off */
          if (strcmp(cptr, "off") == 0) {
            if (spreader.Restart_Info.Flag < 0) {
              spreader.Restart_Info.Flag      = 0;
              spreader.Restart_Info.Num_Times = 0;
            }
          }
          /* check if the user wants all of the time steps */
          else if (strcmp(cptr, "all") == 0) {
            if (spreader.Restart_Info.Flag < 0) {
              spreader.Restart_Info.Flag      = 1;
              spreader.Restart_Info.Num_Times = -1; /* -1 designates read all times */
            }
          }
          /* IGNORED check if the user wants move variable in blocks IGNORED */
          else if (strstr(cptr, "block")) {
            cptr2 = strchr(cptr, '=');
            if (cptr2 == nullptr) {
              fmt::print(stderr, "fatal: must specify a value with \"block\"");
              exit(1);
            }
            cptr2++;
          }
          /* check if the user has a list of time indices to get */
          else if (strstr(cptr, "list")) {
            /* "{" defines the beginning of the group designator */
            cptr2 = strchr(cptr, '{');
            if (cptr2 == nullptr) {
              fmt::print(stderr, fmt::runtime("fatal: list start designator \"{\" not found"));
              exit(1);
            }
            cptr2++;
            cptr3 = strchr(cptr, '}');
            if (cptr3 == nullptr) {
              fmt::print(stderr, fmt::runtime("fatal: list end designator \"}\" not found"));
              exit(1);
            }
            *cptr3 = '\0';

            /* Allocate the time list */
            spreader.Restart_Info.Time_Idx.resize(TLIST_CNT);
            tlist_alloc                     = TLIST_CNT;
            spreader.Restart_Info.Num_Times = 0;

            while (cptr2) {
              /* first check to see if they want to get the last time index */
              if (strncmp(cptr2, "last", 4) == 0) {
                icnt                                                            = 0;
                spreader.Restart_Info.Time_Idx[spreader.Restart_Info.Num_Times] = 0;
              }
              else {
                icnt = sscanf(cptr2, "%d",
                              &(spreader.Restart_Info.Time_Idx[spreader.Restart_Info.Num_Times]));
              }

              if (icnt >= 0) {
                (spreader.Restart_Info.Num_Times)++;
              }

              if (spreader.Restart_Info.Num_Times >= tlist_alloc) {
                tlist_alloc += TLIST_CNT;
                spreader.Restart_Info.Time_Idx.resize(tlist_alloc);
              }
              /* move to the next blank space */
              cptr3 = strchr(cptr2, ' ');
              if (cptr3) {
                /* find the next non-blank space */
                while (*cptr3 == ' ') {
                  cptr3++;
                }
              }
              cptr2 = cptr3;
            }
          }
          else {
            fmt::print(stderr, "warning: unknown restart info suboption {}", cptr);
            exit(1);
          }
          cptr = strtok(nullptr, ",");
        }
      } /* End "if (token_compare(cptr, "restart time list"))" */
      /****** Reserved Space for Variables ******/
      else if (token_compare(cptr, "reserve space")) {

        cptr = strchr(cptr, '\0');
        cptr++;
        strip_string(cptr, " \t\n=");
        cptr = strtok(cptr, ",");

        while (cptr != nullptr) {
          strip_string(cptr, " \t\n");
          string_to_lower(cptr, '=');
          if (strstr(cptr, "nodal")) {
            cptr2 = strchr(cptr, '=');
            if (cptr2 == nullptr) {
              fmt::print(stderr, "Error: integer value must be specified for"
                                 " reserve space.\n");
              return 0;
            }
            cptr2++;
            icnt = sscanf(cptr2, "%d", &Num_Nod_Var);
            if ((icnt <= 0) || (Num_Nod_Var < 0)) {
              fmt::print(stderr, "Error: Invalid value for nodal variable\n");
              return 0;
            }
          }
          else if (strstr(cptr, "elemental")) {
            cptr2 = strchr(cptr, '=');
            if (cptr2 == nullptr) {
              fmt::print(stderr, "Error: integer value must be specified for"
                                 " reserve space.\n");
              return 0;
            }
            cptr2++;
            icnt = sscanf(cptr2, "%d", &Num_Elem_Var);
            if ((icnt <= 0) || (Num_Elem_Var < 0)) {
              fmt::print(stderr, "Error: Invalid value for elemental variable\n");
              return 0;
            }
          }
          else if (strstr(cptr, "global")) {
            cptr2 = strchr(cptr, '=');
            if (cptr2 == nullptr) {
              fmt::print(stderr, "Error: integer value must be specified for"
                                 " reserve space.\n");
              return 0;
            }
            cptr2++;
            icnt = sscanf(cptr2, "%d", &Num_Glob_Var);
            if ((icnt <= 0) || (Num_Glob_Var < 0)) {
              fmt::print(stderr, "Error: Invalid value for global variable\n");
              return 0;
            }
          }
          else if (strstr(cptr, "nodeset")) {
            cptr2 = strchr(cptr, '=');
            if (cptr2 == nullptr) {
              fmt::print(stderr, "Error: integer value must be specified for"
                                 " reserve space.\n");
              return 0;
            }
            cptr2++;
            icnt = sscanf(cptr2, "%d", &Num_Nset_Var);
            if ((icnt <= 0) || (Num_Nset_Var < 0)) {
              fmt::print(stderr, "Error: Invalid value for nodeset variable\n");
              return 0;
            }
          }
          else if (strstr(cptr, "sideset")) {
            cptr2 = strchr(cptr, '=');
            if (cptr2 == nullptr) {
              fmt::print(stderr, "Error: integer value must be specified for"
                                 " reserve space.\n");
              return 0;
            }
            cptr2++;
            icnt = sscanf(cptr2, "%d", &Num_Sset_Var);
            if ((icnt <= 0) || (Num_Sset_Var < 0)) {
              fmt::print(stderr, "Error: Invalid value for sideset variable\n");
              return 0;
            }
          }

          cptr = strtok(nullptr, ",");

        } /* End "while (cptr != nullptr)" */
      } /* End "else if (token_compare(cptr, "reserve space"))" */
      /****** Parallel Disk Information ******/
      else if (token_compare(cptr, "parallel disk info")) {

        cptr = strchr(cptr, '\0');
        cptr++;
        strip_string(cptr, " \t\n=");
        cptr = strtok(cptr, ",");
        strip_string(cptr, " \t\n");
        string_to_lower(cptr, '=');

        /* the first sub-option must be "number" */
        if (!strstr(cptr, "number")) {
          fmt::print(stderr, "Error: First sup-option for disk info must be "
                             "\"number\"\n");
          return 0;
        }

        cptr2 = strchr(cptr, '=');
        if (cptr2 == nullptr) {
          fmt::print(stderr, "Error: integer value must be specified for"
                             " reserve space.\n");
          return 0;
        }
        cptr2++;
        icnt = sscanf(cptr2, "%d", &(PIO_Info.Num_Dsk_Ctrlrs));
        if ((icnt <= 0) || (PIO_Info.Num_Dsk_Ctrlrs <= 0)) {
          fmt::print(stderr, "Error: Invalid value for # of raid controllers\n");
          return 0;
        }

        cptr = strtok(nullptr, ",");
        while (cptr != nullptr) {
          strip_string(cptr, " \t\n");
          string_to_lower(cptr, '=');
          if (strstr(cptr, "list")) {
            /*
             * So, "number" references the length of the list, and
             * I need to do some shuffling to make the new form
             * work with the old code.
             */
            PIO_Info.Dsk_List_Cnt   = PIO_Info.Num_Dsk_Ctrlrs;
            PIO_Info.Num_Dsk_Ctrlrs = 0;

            /* "{" defines the beginning of the list */
            cptr = strchr(cptr, '{');
            if (cptr == nullptr) {
              fmt::print(stderr, "Error: disk list must be specified\n");
              return 0;
            }
            cptr++;

            /* allocate memory for to hold the values */
            PIO_Info.Dsk_List = reinterpret_cast<int *>(
                array_alloc(__FILE__, __LINE__, 1, PIO_Info.Dsk_List_Cnt, sizeof(int)));
            for (i = 0; i < (PIO_Info.Dsk_List_Cnt - 1); i++) {
              sscanf(cptr, "%d", &(PIO_Info.Dsk_List[i]));
              cptr = strtok(nullptr, ", \t;");
            }
            /* last one is a special case */
            sscanf(cptr, "%d}", &(PIO_Info.Dsk_List[i]));
          }
          else if (strstr(cptr, "offset")) {
            cptr2 = strchr(cptr, '=');
            if (cptr2 == nullptr) {
              fmt::print(stderr, "Error: value must be specified with the "
                                 "\"offset\" option.\n");
              return 0;
            }
            cptr2++;
            icnt = sscanf(cptr2, "%d", &(PIO_Info.PDsk_Add_Fact));
            if ((icnt <= 0) || (PIO_Info.PDsk_Add_Fact < 0)) {
              fmt::print(stderr, "Error: Invalid value for offset\n");
              return 0;
            }
          }
          else if (strstr(cptr, "zeros")) {
            PIO_Info.Zeros = 1;
          }
          else if (strstr(cptr, "nosubdirectory")) {
            PIO_Info.NoSubdirectory = 1;
          }
          else if (strstr(cptr, "stage_off")) {
            PIO_Info.Staged_Writes = false;
          }
          else if (strstr(cptr, "stage_on")) {
            PIO_Info.Staged_Writes = true;
          }

          cptr = strtok(nullptr, ",");
        }
      } /* End "else if (token_compare(cptr, "parallel disk info"))" */
      else if (token_compare(cptr, "parallel file location")) {
        cptr = strchr(cptr, '\0');
        cptr++;
        strip_string(cptr, " \t\n=");
        cptr = strtok(cptr, ",");

        while (cptr != nullptr) {
          strip_string(cptr, " \t\n");
          string_to_lower(cptr, '=');
          if (strstr(cptr, "root")) {
            cptr2 = strchr(cptr, '=');
            if (cptr2 == nullptr) {
              fmt::print(stderr, "fatal: must specify a path with \"root\"");
              return 0;
            }
            cptr2++;
            if (strlen(cptr2) == 0) {
              fmt::print(stderr, "fatal: invalid path name specified with"
                                 " \"root\"");
              return 0;
            }
            PIO_Info.Par_Dsk_Root = std::string(cptr2);
          }
          if (strstr(cptr, "subdir")) {
            cptr2 = strchr(cptr, '=');
            if (cptr2 == nullptr) {
              fmt::print(stderr, "fatal: must specify a path with \"subdir\"");
              return 0;
            }
            cptr2++;
            if (strlen(cptr2) == 0) {
              fmt::print(stderr, "fatal: invalid path name specified with"
                                 " \"subdir\"");
              return 0;
            }
            PIO_Info.Par_Dsk_SubDirec = std::string(cptr2);
            if (PIO_Info.Par_Dsk_SubDirec.back() != '/') {
              PIO_Info.Par_Dsk_SubDirec += '/';
            }
          }

          cptr = strtok(nullptr, ",");
        }
      }

    } /* End "if(inp_line[0] != '#')" */
  } /* End "while(fgets(inp_line, MAX_INPUT_STR_LN, file_cmd))" */

  if (Output_File_Base_Name.empty() && !Exo_LB_File.empty()) {
    // User did not specify a base name.  Use the basenmae of the
    // Exo_LB_File instead.
    Output_File_Base_Name = Exo_LB_File;

    // If there is an extension, strip it off...
    size_t found = Output_File_Base_Name.find_last_of('.');
    if (found != std::string::npos) {
      Output_File_Base_Name = Output_File_Base_Name.substr(0, found);
    }
  }
  return 0;
}

template int read_pexoII_info(NemSpread<float, int> &spreader, const char *filename);
template int read_pexoII_info(NemSpread<double, int> &spreader, const char *filename);

template int read_pexoII_info(NemSpread<double, int64_t> &spreader, const char *filename);
template int read_pexoII_info(NemSpread<float, int64_t> &spreader, const char *filename);
