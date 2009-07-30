/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#include "Parser_dh.h"
#include "Mem_dh.h"
 
typedef struct _optionsNode OptionsNode;

struct _parser_dh {
  OptionsNode *head;  
  OptionsNode *tail;  
};

struct _optionsNode {
  char *name;
  char *value;
  OptionsNode *next;  
};

static bool find(Parser_dh p, char *option, OptionsNode** ptr);
static void init_from_default_settings_private(Parser_dh p);


#undef __FUNC__
#define __FUNC__ "Parser_dhCreate"
void Parser_dhCreate(Parser_dh *p)
{
  START_FUNC_DH
  OptionsNode *ptr;
 
  /* allocate storage for object */
  struct _parser_dh* tmp = (struct _parser_dh*)MALLOC_DH(sizeof(struct _parser_dh)); CHECK_V_ERROR;
  *p = tmp;

  /* consruct header node */
  tmp->head = tmp->tail = (OptionsNode*)MALLOC_DH(sizeof(OptionsNode)); CHECK_V_ERROR;
  ptr = tmp->head;
  ptr->next = NULL;
  ptr->name  = (char*)MALLOC_DH(6*sizeof(char)); CHECK_V_ERROR;
  ptr->value = (char*)MALLOC_DH(6*sizeof(char)); CHECK_V_ERROR;
  strcpy(ptr->name, "JUNK");
  strcpy(ptr->value, "JUNK");
  END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "Parser_dhDestroy"
void Parser_dhDestroy(Parser_dh p)
{
  START_FUNC_DH
  OptionsNode *ptr2 = p->head, *ptr1 = ptr2;
  if (ptr1 != NULL) {
    do {
      ptr2 = ptr2->next;
      FREE_DH(ptr1->name);
      FREE_DH(ptr1->value);
      FREE_DH(ptr1);
      ptr1 = ptr2;
    } while (ptr1 != NULL);
  }
  FREE_DH(p);
  END_FUNC_DH
}


#undef __FUNC__
#define __FUNC__ "Parser_dhUpdateFromFile"
void Parser_dhUpdateFromFile(Parser_dh p, char *filename)
{
  START_FUNC_DH_2
  char line[80], name[80], value[80];
  FILE *fp;

  if ((fp = fopen(filename, "r")) == NULL) {
    sprintf(msgBuf_dh, "can't open >>%s<< for reading", filename);
    SET_INFO(msgBuf_dh);
  } else {
    sprintf(msgBuf_dh, "updating parser from file: >>%s<<", filename);
    SET_INFO(msgBuf_dh);
    while (!feof(fp)) {
      if (fgets(line, 80, fp) == NULL) break;
      if (line[0] != '#') { 
        if (sscanf(line, "%s %s", name, value) != 2) break;
        Parser_dhInsert(p, name, value);   
      }
    }
    fclose(fp);
  }
  END_FUNC_DH_2
}

#undef __FUNC__
#define __FUNC__ "Parser_dhInit"
void Parser_dhInit(Parser_dh p, int argc, char *argv[])
{
  START_FUNC_DH_2
  int j;

  /* read option names and values from default database */
/*  Parser_dhUpdateFromFile(p, MASTER_OPTIONS_LIST); CHECK_V_ERROR;
*/
  init_from_default_settings_private(p); CHECK_V_ERROR;

  /* attempt to update from "./database" in local directory */
  Parser_dhUpdateFromFile(p, "./database"); CHECK_V_ERROR;

  /* attempt to update from specified file */
  for (j=1; j<argc; ++j) {
    if (strcmp(argv[j],"-db_filename") == 0) {  
       ++j;
      if (j < argc) {
        Parser_dhUpdateFromFile(p, argv[j]); CHECK_V_ERROR;
      }
    }
  }

  /* update from command-line options and values */
  {
  int i = 0;
  while (i < argc) {
    if (argv[i][0] == '-') {
      char value[] = { "1" };  /* option's default value */
      bool flag = false;       /* yuck! flag for negative numbers */
      if (i+1 < argc && argv[i+1][0] == '-' && argv[i+1][1] == '-') {
        flag = true;
      }

      if ( (i+1 == argc || argv[i+1][0] == '-') && !flag ) {
        Parser_dhInsert(p, argv[i], value);
      } else if (flag) {
        Parser_dhInsert(p, argv[i], argv[i+1]+1); /* insert a negative number */
      } else {
        Parser_dhInsert(p, argv[i], argv[i+1]);
      }
    }
    ++i;
  }}
  END_FUNC_DH_2
}


#undef __FUNC__
#define __FUNC__ "Parser_dhHasSwitch"
bool Parser_dhHasSwitch(Parser_dh p, char* s)
{
  START_FUNC_DH_2
  bool has_switch = false;
  OptionsNode *node;

  if (p != NULL && find(p,s,&node)) {
    if (! strcmp(node->value, "0")) {
      has_switch = false;
    } else if  (! strcmp(node->value, "false")) {
      has_switch = false;
    } else if  (! strcmp(node->value, "False")) {
      has_switch = false;
    } else if  (! strcmp(node->value, "FALSE")) {
      has_switch = false;
    } else {
      has_switch = true;
    }
  }
  END_FUNC_VAL_2(has_switch)
}

/* returns false if option isn't found, or if
 * its value is zero.
 */
#undef __FUNC__
#define __FUNC__ "Parser_dhReadInt"
bool Parser_dhReadInt(Parser_dh p, char* in, int* out)
{
  START_FUNC_DH_2
  bool has_switch = false;
  OptionsNode *node;

  if (p != NULL && find(p,in,&node)) {
    *out = atoi(node->value);
    if (! strcmp(node->value, "0")) {
      has_switch = false;
    } else {
      has_switch = true;
    }
  }
  END_FUNC_VAL_2(has_switch)
}


#undef __FUNC__
#define __FUNC__ "Parser_dhReadDouble"
bool Parser_dhReadDouble(Parser_dh p, char* in, double *out)
{
  START_FUNC_DH_2
  bool optionExists = false;
  OptionsNode *node;

  if (p != NULL && find(p,in,&node)) {
    *out = atof(node->value);
    optionExists = true;
  }
  END_FUNC_VAL_2(optionExists)
}

#undef __FUNC__
#define __FUNC__ "Parser_dhReadString"
bool Parser_dhReadString(Parser_dh p, char* in, char **out)
{
  START_FUNC_DH_2
  bool optionExists = false;
  OptionsNode *node;

  if (p != NULL && find(p,in,&node)) {
    *out = node->value;
    optionExists = true;
  } 
  END_FUNC_VAL_2(optionExists)
}


#undef __FUNC__
#define __FUNC__ "Parser_dhPrint"
void Parser_dhPrint(Parser_dh p, FILE *fp, bool allPrint)
{
  START_FUNC_DH_2
  OptionsNode *ptr = p->head;

  if (fp == NULL) SET_V_ERROR("fp == NULL");

  if (myid_dh == 0 || allPrint) {
    fprintf(fp, "------------------------ registered options:\n");
    if (ptr == NULL) {
      fprintf(fp, "Parser object is invalid; nothing to print!\n");
    } else {
      ptr = ptr->next;
      while (ptr != NULL) {
        fprintf(fp, "   %s  %s\n", ptr->name, ptr->value);
        fflush(fp);
        ptr = ptr->next;
      } 
    } 
    fprintf(fp, "\n");
    fflush(fp);
  }
  END_FUNC_DH_2
}

#undef __FUNC__
#define __FUNC__ "Parser_dhInsert"
void Parser_dhInsert(Parser_dh p, char *option, char *value)
{
  START_FUNC_DH_2
  OptionsNode *node;
  int length;

  if (p == NULL) goto PARSER_NOT_INITED;

  /* if option is already in the list, update its value */
  if (find(p, option,&node)) {
    int length2 = strlen(node->value)+1;
    length = strlen(value)+1;
    if (length2 < length) {
      FREE_DH(node->value);
      node->value  =  (char*)MALLOC_DH(length*sizeof(char)); CHECK_V_ERROR;
    }
    strcpy(node->value, value);
  }
  /* otherwise, add a new node to the list */
  else {
    node = p->tail;
    p->tail = node->next = (OptionsNode*)MALLOC_DH(sizeof(OptionsNode)); CHECK_V_ERROR;
    node = node->next;
    length = strlen(option)+1;
    node->name = (char*)MALLOC_DH(length*sizeof(char)); CHECK_V_ERROR;
    strcpy(node->name, option);
    length = strlen(value)+1;
    node->value = (char*)MALLOC_DH(length*sizeof(char)); CHECK_V_ERROR;
    strcpy(node->value, value);
    node->next = NULL;
  } 

PARSER_NOT_INITED:
      ;

  END_FUNC_DH_2
}

#undef __FUNC__
#define __FUNC__ "find"
bool find(Parser_dh p, char *option, OptionsNode** ptr)
{
  START_FUNC_DH_2
  OptionsNode *tmpPtr = p->head;
  bool foundit = false;
  while (tmpPtr != NULL) {
    if (strcmp(tmpPtr->name,option) == 0) {
      foundit = true;
      *ptr = tmpPtr;
      break;
    }
    tmpPtr = tmpPtr->next;
  }
  END_FUNC_VAL_2(foundit)
}


#undef __FUNC__
#define __FUNC__ "init_from_default_settings_private"
void init_from_default_settings_private(Parser_dh p)
{
  START_FUNC_DH_2
  /* default is to intercept certain signals
     (floating point error, segmentation violation, etc.)
  */
  Parser_dhInsert(p, "-sig_dh", "1"); CHECK_V_ERROR;

  /* used by MetGenFD */
  Parser_dhInsert(p, "-px", "1"); CHECK_V_ERROR;
  Parser_dhInsert(p, "-py", "1"); CHECK_V_ERROR;
  Parser_dhInsert(p, "-pz", "0"); CHECK_V_ERROR;
  Parser_dhInsert(p, "-m", "4"); CHECK_V_ERROR;

  Parser_dhInsert(p, "-xx_coeff", "-1.0"); CHECK_V_ERROR;
  Parser_dhInsert(p, "-yy_coeff", "-1.0"); CHECK_V_ERROR;
  Parser_dhInsert(p, "-zz_coeff", "-1.0"); CHECK_V_ERROR;

  Parser_dhInsert(p, "-level", "1"); CHECK_V_ERROR;

  Parser_dhInsert(p, "-printStats", "0"); CHECK_V_ERROR;
  END_FUNC_DH_2
}
