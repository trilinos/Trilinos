/*
 * Copyright(C) 2011 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 * * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *           
 * * Redistributions in binary form must reproduce the above
 *   copyright notice, this list of conditions and the following
 *   disclaimer in the documentation and/or other materials provided
 *   with the distribution.
 *                         
 * * Neither the name of Sandia Corporation nor the names of its
 *   contributors may be used to endorse or promote products derived
 *   from this software without specific prior written permission.
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
/* exodus II to matlab m file, copy of
   exo2mat. 
   exo2mat was written by mrtabba
   exo2m modifications by gmreese to permit usage on machines without
     the matlab libraries.

   modified by D. Todd Griffith, 01/12/2006
      to include changes made in versions 1.4 through 1.6 on the 
      SEACAS tools repository as the previous modifications dated
      12/08 and 12/15/2005 were made starting with Version 1.3.
      In particular, the only changes are those made here are those
      made in Version 1.6 which include a special provision 
      for exodus files which contain no distribution factors 
   modified by D. Todd Griffith, 12/15/2005
      to write distribution factors as double precision type
   modified by D. Todd Griffith  12/08/2005
      to include complete writing of side set and node set information 
      so it will be available for the mat2exo conversion stage

*/

#include <vector>
#include <algorithm>
#include <iostream>

#include <assert.h>                     // for assert
#include <stddef.h>                     // for size_t
#include <stdio.h>                      // for fprintf, printf, sprintf, etc
#include <stdlib.h>                     // for free, calloc, exit, malloc
#include <string.h>                     // for strcat, strlen, strcpy, etc
#include "add_to_log.h"                 // for add_to_log
#include "exodusII.h"                   // for ex_get_variable_param, etc
#include "matio.h"                      // for Mat_VarCreate, Mat_VarFree, etc

#if __cplusplus > 199711L
#define TOPTR(x) x.data()
#else
#define TOPTR(x) (x.empty() ? NULL : &x[0])
#endif

#define EXT ".mat"
int textfile=0;

FILE* m_file=NULL;     /* file for m file output */
mat_t *mat_file=NULL;  /* file for binary .mat output */

static const char *qainfo[] =
{
  "exo2mat",
  "2015/10/20",
  "3.00",
};


void usage()
{
  std::cout << "exo2mat [options] exodus_file_name.\n"
	    << "   the exodus_file_name is required (exodus only).\n"
	    << "   Options:\n"
	    << "   -t   write a text (.m) file rather than a binary .mat\n"
	    << "   -o   output file name (rather than auto generate)\n"
	    << "   -v5  output version 5 mat file (default)\n"
	    << "   -v73 output version 7.3 mat file (hdf5-based)\n"
	    << "   -v7.3 output version 7.3 mat file (hdf5-based)\n"
	    << " ** note **\n"
	    << "Binary files are written by default on all platforms.\n";
}

/* put a string into an m file. If the string has
   line feeds, we put it as ints, and use 'char()' to convert it */
void mPutStr (const char *name, char *str)
{
  assert(m_file!=0);
  if (strchr(str,'\n')==0)
    fprintf(m_file,"%s='%s';\n",name,str);
  else {
    fprintf(m_file,"%s=[",name);
    size_t i;
    size_t j;
    for (j=i=0;i<strlen(str);i++,j++){
      if (j>=20){
	j=0;
	fprintf(m_file,"...\n");
      }
      fprintf(m_file,"%d ",str[i]);
    }
    fprintf(m_file,"];\n");
    fprintf(m_file,"%s=char(%s);\n",name,name);
  }
}

/* put double array in m file */
void mPutDbl (const char *name,int n1,int n2,double *pd)
{
  assert(m_file != NULL);
  if ( n1==1 && n2 ==1 ){
    fprintf(m_file,"%s=%15.8e;\n",name,*pd);
    return;
  }
  fprintf(m_file,"%s=zeros(%d,%d);\n",name,n1,n2);
  for (int i=0;i<n1;i++)
    for (int j=0;j<n2;j++)
      fprintf(m_file,"%s(%d,%d)=%15.8e;\n",name,i+1,j+1,pd[i*n2+j]);
}


/* put integer array in m file */
void mPutInt (const char *name, int pd)
{
  assert(m_file != NULL);
  fprintf(m_file,"%s=%d;\n",name,pd);
  return;
}

/* put integer array in m file */
void mPutInt (const char *name,int n1,int n2, int *pd)
{
  assert(m_file != NULL);
  if ( n1==1 && n2 ==1 ){
    fprintf(m_file,"%s=%d;\n",name,*pd);
    return;
  }
  fprintf(m_file,"%s=zeros(%d,%d);\n",name,n1,n2);
  for (int i=0;i<n1;i++)
    for (int j=0;j<n2;j++)
      fprintf(m_file,"%s(%d,%d)=%d;\n",name,i+1,j+1,pd[i*n2+j]);
}


/* put string in mat file*/
void matPutStr (const char *name, char *str)
{
  matvar_t *matvar = NULL;
  size_t dims[2];

  dims[0] = 1;
  dims[1] = strlen(str);

  matvar = Mat_VarCreate(name, MAT_C_CHAR, MAT_T_UINT8, 2, dims, str, MAT_F_DONT_COPY_DATA);
  Mat_VarWrite(mat_file, matvar, MAT_COMPRESSION_NONE);
  Mat_VarFree(matvar);
}

/* put double in mat file*/
void matPutDbl (const char *name,int n1,int n2,double *pd)
{
  matvar_t *matvar = NULL;
  
  size_t dims[2];
  dims[0] = n1;
  dims[1] = n2;
  
  matvar = Mat_VarCreate(name, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, pd, MAT_F_DONT_COPY_DATA);
  Mat_VarWrite(mat_file, matvar, MAT_COMPRESSION_ZLIB);
  Mat_VarFree(matvar);
}

/* put integer in mat file*/
void matPutInt (const char *name,int n1,int n2, int *pd)
{
  matvar_t *matvar = NULL;
  
  size_t dims[2];
  dims[0] = n1;
  dims[1] = n2;
  
  matvar = Mat_VarCreate(name, MAT_C_INT32, MAT_T_INT32, 2, dims, pd, MAT_F_DONT_COPY_DATA);
  Mat_VarWrite(mat_file, matvar, MAT_COMPRESSION_ZLIB);
  Mat_VarFree(matvar);
}

/*----------------------*/

/* wrappers for the output routine types */
void PutStr (const char *name, char *str)
{
  if ( textfile )
    mPutStr(name,str);
  else
    matPutStr(name,str);
}

void PutInt (const char *name,int pd)
{
  if ( textfile )
    mPutInt(name,pd);
  else
    matPutInt(name,1,1,&pd);
}
    
void PutInt (const char *name,int n1,int n2,int *pd)
{
  if ( textfile )
    mPutInt(name,n1,n2,pd);
  else
    matPutInt(name,n1,n2,pd);
}
    
void PutDbl (const char *name,int n1,int n2,double *pd)
{
  if ( textfile )
    mPutDbl(name,n1,n2,pd);
  else
    matPutDbl(name,n1,n2,pd);
}
    

/**********************************************************************/
/* remove an argument from the list */
void del_arg(int *argc, char* argv[], int j)
{
  for (int jj=j+1;jj<*argc;jj++)
    argv[jj-1]=argv[jj];
  (*argc)--;
  argv[*argc]=0;
}
/**********************************************************************/
int main (int argc, char *argv[])
{

  char  
    *str,**str2,*line, *oname, *dot, *filename;

  const char* ext=EXT;

  int   
    i,j,k,n,n1,n2,cpu_word_size,io_word_size,exo_file,err,
    num_axes,num_nodes,num_elements,num_blocks,
    num_side_sets,num_node_sets,num_time_steps,
    num_info_lines,num_global_vars,
    num_nodal_vars,num_element_vars,num_nodeset_vars, num_sideset_vars,
    nstr2, has_ss_dfac;
    
  float
    exo_version;

  int mat_version = 50;
  
  oname=0;

  /* process arguments */
  for (j=1; j< argc; j++){
    if ( strcmp(argv[j],"-t")==0){    /* write text file (*.m) */
      del_arg(&argc,argv,j);
      textfile=1;
      j--;
      continue;
    }
    if ( strcmp(argv[j],"-h")==0){    /* write help info */
      del_arg(&argc,argv,j);
      usage();
      exit(1);
    }
    if ( strcmp(argv[j],"-v73")==0){    /* Version 7.3 */
      del_arg(&argc,argv,j);
      mat_version = 73;
      j--;
      continue;
    }
    // This matches the option used in matlab
    if ( (strcmp(argv[j],"-v7.3")==0) || (strcmp(argv[j],"-V7.3")==0)){    /* Version 7.3 */
      del_arg(&argc,argv,j);
      mat_version = 73;
      j--;
      continue;
    }
    if ( strcmp(argv[j],"-v5")==0){    /* Version 5 (default) */
      del_arg(&argc,argv,j);
      mat_version = 50;
      j--;
      continue;
    }
    if ( strcmp(argv[j],"-o")==0){    /* specify output file name */
      del_arg(&argc,argv,j);
      if ( argv[j] ){
	oname=(char*)calloc(strlen(argv[j])+10,sizeof(char));
	strcpy(oname,argv[j]);
	del_arg(&argc,argv,j);
	std::cout << "output file: " << oname << "\n";
      }
      else {
	std::cerr << "ERROR: Invalid output file specification.\n";
	return 2;
      }
      j--;

      continue;
    }
  }

  /* QA Info */
  printf("%s: %s, %s\n", qainfo[0], qainfo[2], qainfo[1]);
  
  /* usage message*/
  if(argc != 2){
    usage();
    exit(1);
  }

  /* open output file */
  if ( textfile )
    ext=".m";

  if ( !oname ){
    filename = (char*)malloc( strlen(argv[1])+10);
    strcpy(filename,argv[1]);
    dot=strrchr(filename,'.');
    if ( dot ) *dot='\0';
    strcat(filename,ext);
  }
  else {
    filename=oname;
  }

  if ( textfile ){
    m_file = fopen(filename,"w");
    if (!m_file ){
      std::cerr << "ERROR: Unable to open " << filename << "\n";
      exit(1);
    }
  }
  else {
    if (mat_version == 50) {
      mat_file = Mat_CreateVer(filename, NULL, MAT_FT_MAT5);
    } else if (mat_version == 73) {
      mat_file = Mat_CreateVer(filename, NULL, MAT_FT_MAT73);
    }
    if (mat_file == NULL) {
      std::cerr << "ERROR: Unable to create matlab file " << filename << "\n";
      exit(1);
    }
  }

  /* word sizes */
  cpu_word_size=sizeof(double);
  io_word_size=0;

  /* open exodus file */
  exo_file=ex_open(argv[1],EX_READ,&cpu_word_size,&io_word_size,&exo_version);
  if (exo_file < 0){
    std::cerr << "ERROR: Cannot open " << argv[1] << "\n";
    exit(1);
  }

  /* print */
  std::cout << "\ttranslating " << argv[1] << " to " << filename << "...\n";

  /* read database paramters */
  line=(char *) calloc ((MAX_LINE_LENGTH+1),sizeof(char));
  ex_get_init(exo_file,line,
	      &num_axes,&num_nodes,&num_elements,&num_blocks,
	      &num_node_sets,&num_side_sets);
  num_info_lines = ex_inquire_int(exo_file,EX_INQ_INFO);
  num_time_steps = ex_inquire_int(exo_file,EX_INQ_TIME);
  ex_get_variable_param(exo_file,EX_GLOBAL,&num_global_vars);
  ex_get_variable_param(exo_file,EX_NODAL,&num_nodal_vars);
  ex_get_variable_param(exo_file,EX_ELEM_BLOCK,&num_element_vars);
  ex_get_variable_param(exo_file,EX_NODE_SET,&num_nodeset_vars);
  ex_get_variable_param(exo_file,EX_SIDE_SET,&num_sideset_vars);


  /* export paramters */
  PutInt("naxes",  num_axes);
  PutInt("nnodes", num_nodes);
  PutInt("nelems", num_elements);
  PutInt("nblks",  num_blocks);
  PutInt("nnsets", num_node_sets);
  PutInt("nssets", num_side_sets);
  PutInt("nsteps", num_time_steps);
  PutInt("ngvars", num_global_vars);
  PutInt("nnvars", num_nodal_vars);
  PutInt("nevars", num_element_vars);
  PutInt("nnsvars",num_nodeset_vars);
  PutInt("nssvars",num_sideset_vars);

  /* allocate -char- scratch space*/
  n = num_info_lines;
  n = std::max(n, num_global_vars);
  n = std::max(n, num_nodal_vars);
  n = std::max(n, num_element_vars);
  n = std::max(n, num_blocks);
  nstr2 = n;
  str2= (char **) calloc (n,sizeof(char *));
  for (i=0;i<nstr2;i++)
    str2[i]=(char *) calloc ((MAX_LINE_LENGTH+1),sizeof(char));
  str= (char *) calloc ((MAX_LINE_LENGTH+1)*n,sizeof(char));

  /* title */
  PutStr("Title",line);

  /* information records */
  if (num_info_lines > 0 ){
    ex_get_info(exo_file,str2);
    str[0]='\0';
    for (i=0;i<num_info_lines;i++) {
      if (strlen(str2[i]) > 0) {
	strcat(str, str2[i]);
	strcat(str, "\n");
      }
    }
    PutStr("info",str);
    str[0]='\0';
    for (i=0;i<num_info_lines;i++) {
      if (strlen(str2[i]) > 0 && strncmp(str2[i],"cavi",4)==0) {
	strcat(str, str2[i]);
	strcat(str, "\n");
      }
    }
    PutStr("cvxp",str);
  }

  /* nodal coordinates */
  {
    std::vector<double> x, y, z;
    x.resize(num_nodes);
    if (num_axes >= 2) y.resize(num_nodes);
    if (num_axes == 3) z.resize(num_nodes);
    ex_get_coord(exo_file,TOPTR(x), TOPTR(y), TOPTR(z));
    PutDbl("x0", num_nodes, 1, TOPTR(x));
    if (num_axes >= 2) {
      PutDbl("y0", num_nodes, 1, TOPTR(y));
    }
    if (num_axes == 3){ 
      PutDbl("z0",num_nodes,1, TOPTR(z));
    }
  }  

  /* side sets */
  std::vector<int> ids;
  if(num_side_sets > 0) {
    ids.resize(num_side_sets);
    ex_get_ids(exo_file,EX_SIDE_SET,TOPTR(ids));
    PutInt( "ssids",num_side_sets, 1,TOPTR(ids));
    std::vector<int> nsssides(num_side_sets);
    std::vector<int> nssdfac(num_side_sets);
    std::vector<int> iscr;
    std::vector<int> jscr;
    std::vector<double> scr;
    std::vector<int> elem_list;
    std::vector<int> side_list;
    std::vector<int> junk;
    for (i=0;i<num_side_sets;i++) {
      ex_get_set_param(exo_file,EX_SIDE_SET, ids[i],&n1,&n2);
      nsssides[i]=n1; /* dgriffi */
      nssdfac[i]=n2;  /* dgriffi */
      /*
       * the following provision is from Version 1.6 when there are no
       * distribution factors in exodus file
       */
      has_ss_dfac = (n2 != 0);
      if(n2==0 || n1==n2){
	
	std::cerr << "WARNING: Exodus II file does not contain distribution factors.\n";
	
	/* n1=number of faces, n2=number of df */
	/* using distribution factors to determine number of nodes in the sideset
	   causes a lot grief since some codes do not output distribution factors
	   if they are all equal to 1. mkbhard: I am using the function call below
	   to figure out the total number of nodes in this sideset. Some redundancy
	   exists, but it works for now */

	junk.resize(n1);
	ex_get_side_set_node_count(exo_file,ids[i],TOPTR(junk));
	n2=0; /* n2 will be equal to the total number of nodes in the sideset */
	for (j=0;j<n1;j++)
	  n2+=junk[j];
      }
	
      iscr.resize(n1);
      jscr.resize(n2);
      ex_get_side_set_node_list(exo_file,ids[i],TOPTR(iscr),TOPTR(jscr));
      /* number-of-nodes-per-side list */
      sprintf(str,"ssnum%02d",i+1);
      PutInt(str,n1,1,TOPTR(iscr)); 
      /* nodes list */
      sprintf(str,"ssnod%02d",i+1);
      PutInt(str,n2,1,TOPTR(jscr));

      /* distribution-factors list */
      scr.resize(n2);
      if (has_ss_dfac) {
	ex_get_side_set_dist_fact(exo_file,ids[i], TOPTR(scr));
      } else {
	for (j=0; j<n2; j++) {
	  scr[j] = 1.0;
	}
      }
      sprintf(str,"ssfac%02d",i+1);
      PutDbl(str,n2,1,TOPTR(scr));

      /* element and side list for side sets (dgriffi) */
      elem_list.resize(n1);
      side_list.resize(n1);
      ex_get_set(exo_file,EX_SIDE_SET,ids[i],TOPTR(elem_list),TOPTR(side_list));
      sprintf(str,"ssside%02d",i+1);
      PutInt(str,n1,1,TOPTR(side_list));
      sprintf(str,"sselem%02d",i+1);
      PutInt(str,n1,1,TOPTR(elem_list));
    }
    /* Store # sides and # dis. factors per side set (dgriffi) */
    PutInt("nsssides",num_side_sets,1,TOPTR(nsssides));
    PutInt("nssdfac",num_side_sets,1,TOPTR(nssdfac));
  }

  /* node sets (section by dgriffi) */
  if(num_node_sets > 0){
    std::vector<int> iscr;
    std::vector<double> scr;
    std::vector<int> ids(num_node_sets);
    ex_get_ids(exo_file,EX_NODE_SET, TOPTR(ids));
    PutInt( "nsids",num_node_sets, 1,TOPTR(ids));
    std::vector<int> nnsnodes(num_node_sets);
    std::vector<int> nnsdfac(num_node_sets);
    for (i=0;i<num_node_sets;i++){
      ex_get_set_param(exo_file,EX_NODE_SET,ids[i],&n1,&n2);
      iscr.resize(n1);
      ex_get_node_set(exo_file,ids[i],TOPTR(iscr));
      /* nodes list */
      sprintf(str,"nsnod%02d",i+1);
      PutInt(str,n1,1,TOPTR(iscr));
      {
	/* distribution-factors list */
	scr.resize(n2);
	ex_get_node_set_dist_fact(exo_file,ids[i],TOPTR(scr));  
	sprintf(str,"nsfac%02d",i+1);
	PutDbl(str,n2,1,TOPTR(scr));
      }
      nnsnodes[i]=n1;
      nnsdfac[i]=n2;
    }

    /* Store # nodes and # dis. factors per node set */
    PutInt("nnsnodes",num_node_sets,1,TOPTR(nnsnodes));
    PutInt("nnsdfac",num_node_sets,1,TOPTR(nnsdfac));
  }

  /* element blocks */
  std::vector<int> num_elem_in_block(num_blocks);
  {
    ids.resize(num_blocks);
    std::vector<int> iscr;
    ex_get_ids(exo_file,EX_ELEM_BLOCK,TOPTR(ids));
    PutInt( "blkids",num_blocks, 1,TOPTR(ids));
    for (i=0;i<num_blocks;i++) {
      ex_get_elem_block(exo_file,ids[i],str2[i],&n,&n1,&n2);
      num_elem_in_block[i]=n;
      iscr.resize(n*n1);
      ex_get_conn(exo_file,EX_ELEM_BLOCK,ids[i],TOPTR(iscr), NULL, NULL);
      sprintf(str,"blk%02d",i+1);
      PutInt(str,n1,n,TOPTR(iscr));
    }
    str[0]='\0';
    for (i=0;i<num_blocks;i++) {
      strcat(str, str2[i]);
      strcat(str, "\n");
    }
    PutStr("blknames",str);
  }

  /* time values */
  if (num_time_steps > 0 ) {
    std::vector<double> scr(num_time_steps);
    ex_get_all_times (exo_file, TOPTR(scr));
    PutDbl( "time", num_time_steps, 1, TOPTR(scr));
  }

  /* global variables */
  if (num_global_vars > 0 ) {
    ex_get_variable_names(exo_file,EX_GLOBAL,num_global_vars,str2);
    str[0]='\0';
    for (i=0;i<num_global_vars;i++) {
      strcat(str, str2[i]);
      strcat(str, "\n");
    }
    PutStr("gnames",str);
    std::vector<double> scr(num_time_steps);
    for (i=0;i<num_global_vars;i++){
      sprintf(str,"gvar%02d",i+1);
      ex_get_glob_var_time(exo_file,i+1,1,num_time_steps,TOPTR(scr));
      PutDbl(str,num_time_steps,1,TOPTR(scr));
    }
  }

  /* nodal variables */
  if (num_nodal_vars > 0 ) {
    ex_get_variable_names(exo_file,EX_NODAL,num_nodal_vars,str2);
    str[0]='\0';
    for (i=0;i<num_nodal_vars;i++) {
      strcat(str, str2[i]);
      strcat(str, "\n");
    }
    PutStr("nnames",str);
    std::vector<double> scr(num_nodes*num_time_steps);
    for (int i=0; i<num_nodal_vars; i++){
      sprintf(str,"nvar%02d",i+1);
      for (int j=0; j<num_time_steps; j++)
	ex_get_nodal_var(exo_file,j+1,i+1,num_nodes,
			 &scr[num_nodes*j]);
      PutDbl(str,num_nodes,num_time_steps,TOPTR(scr));
    }
  }

  /* element variables */
  if (num_element_vars > 0 ) {
    ex_get_variable_names(exo_file,EX_ELEM_BLOCK,num_element_vars,str2);
    str[0]='\0';
    for (i=0;i<num_element_vars;i++) {
      strcat(str, str2[i]);
      strcat(str, "\n");
    }
    PutStr("enames",str);
    /* truth table */
    std::vector<int> iscr(num_element_vars*num_blocks);
    ex_get_elem_var_tab(exo_file,num_blocks,num_element_vars,TOPTR(iscr));
    std::vector<double> scr(num_elements * num_time_steps);
    for (i=0;i<num_element_vars;i++){
      std::fill(scr.begin(), scr.end(), 0.0);
      n=0;
      sprintf(str,"evar%02d",i+1);
      for (j=0;j<num_time_steps;j++){
	for (k=0;k<num_blocks;k++){ 
          if(iscr[num_element_vars*k+i]==1)
	    ex_get_elem_var(exo_file,j+1,i+1,ids[k],num_elem_in_block[k],&scr[n]);
	  n=n+num_elem_in_block[k];
	      
	}
      }
      PutDbl(str,num_elements,num_time_steps,TOPTR(scr));
    }
  }
  // Clear out num_elem_in_block...
 
  /* node and element number maps */
  ex_opts(0);  /* turn off error reporting. It is not an error to have no map*/
  ids.resize(num_nodes);
  err = ex_get_node_num_map(exo_file,TOPTR(ids));
  if ( err==0 ){
    PutInt("node_num_map",num_nodes,1,TOPTR(ids));
  }

  ids.resize(num_elements);
  err = ex_get_elem_num_map(exo_file,TOPTR(ids));
  if ( err==0 ){
    PutInt("elem_num_map",num_elements,1,TOPTR(ids));
  }


  /* close exo file */
  ex_close(exo_file);
  
  /* close mat file */
  if ( textfile )
    fclose(m_file);
  else
    Mat_Close(mat_file);

  /* */
  std::cout << "done...\n";

  free(filename);
  free(line);
  
  free(str);
  for (i=0;i<nstr2;i++)
    free(str2[i]);
  free(str2);
  

  /* exit status */
  add_to_log("exo2mat", 0);
  return(0);
}
