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

   Note:
     If exo2mat is available, it should be used. However,
     the mathsoft libraries are a moving target. This crutch
     delivers almost the same functionality, but does so by
     writing an .m file rather than the binary .mat file.
     gmreese. April 1, 2003.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <assert.h>
#include "exodusII.h"

#include "add_to_log.h"
#include "matio.h"

#define EXT ".mat"
int textfile=0;

FILE* m_file=0;     /* file for m file output */
mat_t *mat_file=0;  /* file for binary .mat output */

static char *qainfo[] =
{
  "exo2mat",
  "2012/07/09",
  "2.01",
};


/* put a string into an m file. If the string has
   line feeds, we put it as ints, and use 'char()' to convert it */
void mPutStr (char *name,char *str)
{
  unsigned int i,j;
  assert(m_file!=0);
  if (strchr(str,'\n')==0)
    fprintf(m_file,"%s='%s';\n",name,str);
  else {
    fprintf(m_file,"%s=[",name);
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
void mPutDbl (char *name,int n1,int n2,double *pd)
{
  int i,j;
  assert(m_file!=0);
  if ( n1==1 && n2 ==1 ){
    fprintf(m_file,"%s=%15.8e;\n",name,*pd);
    return;
  }
  fprintf(m_file,"%s=zeros(%d,%d);\n",name,n1,n2);
  for (i=0;i<n1;i++)
    for (j=0;j<n2;j++)
      fprintf(m_file,"%s(%d,%d)=%15.8e;\n",name,i+1,j+1,pd[i*n2+j]);
}


/* put integer array in m file */
void mPutInt (char *name,int n1,int n2, int *pd)
{
  int i,j;
  assert(m_file!=0);
  if ( n1==1 && n2 ==1 ){
    fprintf(m_file,"%s=%d;\n",name,*pd);
    return;
  }
  fprintf(m_file,"%s=zeros(%d,%d);\n",name,n1,n2);
  for (i=0;i<n1;i++)
    for (j=0;j<n2;j++)
      fprintf(m_file,"%s(%d,%d)=%d;\n",name,i+1,j+1,pd[i*n2+j]);
}


/* put string in mat file*/
void matPutStr (char *name,char *str)
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
void matPutDbl (char *name,int n1,int n2,double *pd)
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
void matPutInt (char *name,int n1,int n2, int *pd)
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
void PutStr (char *name,char *str)
{
  if ( textfile )
    mPutStr(name,str);
  else
    matPutStr(name,str);
}

void PutDbl (char *name,int n1,int n2,double *pd)
{
  if ( textfile )
    mPutDbl(name,n1,n2,pd);
  else
    matPutDbl(name,n1,n2,pd);
}
    
void PutInt (char *name,int n1,int n2,int *pd)
{
  if ( textfile )
    mPutInt(name,n1,n2,pd);
  else
    matPutInt(name,n1,n2,pd);
}
    

/**********************************************************************/
/* remove an argument from the list */
void del_arg(int *argc, char* argv[], int j)
{
  int jj;
  for (jj=j+1;jj<*argc;jj++)
    argv[jj-1]=argv[jj];
  (*argc)--;
  argv[*argc]=0;
}
/**********************************************************************/
int main (int argc, char *argv[])
{

  char  
    *str,**str2,*(*qa_records)[4],*line, *oname, *dot, *filename;

  const char* ext=EXT;

  int   
    i,j,k,n,n1,n2,cpu_word_size,io_word_size,exo_file,err,
    num_axes,num_nodes,num_elements,num_blocks,
    num_side_sets,num_node_sets,num_time_steps,
    num_qa_lines,num_info_lines,num_global_vars,
    num_nodal_vars,num_element_vars,num_nodeset_vars, num_sideset_vars,
    *ids,*iscr,*num_elem_in_block,*junk,
    *elem_list,*side_list,
    *nsssides,*nssdfac,
    *nnsnodes,*nnsdfac,
    nstr2, has_ss_dfac;
    
  float
    exo_version;

  double
    *scr,*x,*y,*z;

  oname=0;

  /* process arguments */
  for (j=1; j< argc; j++){
    if ( strcmp(argv[j],"-t")==0){    /* write text file (*.m) */
      del_arg(&argc,argv,j);
      textfile=1;
      j--;
      continue;
    }
    if ( strcmp(argv[j],"-o")==0){    /* specify output file name */
      del_arg(&argc,argv,j);
      if ( argv[j] ){
         oname=(char*)calloc(strlen(argv[j])+10,sizeof(char));
	 strcpy(oname,argv[j]);
	 del_arg(&argc,argv,j);
	 printf("output file: %s\n",oname);
      }
      else {
         fprintf(stderr,"Invalid output file specification.\n");
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
    printf("%s [options] exodus_file_name.\n",argv[0]);
    printf("   the exodus_file_name is required (exodusII only).\n");
    printf("   Options:\n");
    printf("     -t write a text (.m) file rather than a binary .mat\n");
    printf("     -o output file name (rather than auto generate)\n");
    printf(" ** note **\n");
    printf("Binary files are written by default on all platforms with");
    printf(" available libraries.\n");
    exit(1);
  }

  /* open output file */
  if ( textfile )
    ext=".m";

  if ( !oname ){
      filename = (char*)malloc( strlen(argv[1])+10);
      strcpy(filename,argv[1]);
      dot=strrchr(filename,'.');
      if ( dot ) *dot=0;
      strcat(filename,ext);
  }
  else {
      filename=oname;
  }

  if ( textfile ){
    m_file = fopen(filename,"w");
    if (!m_file ){
      fprintf(stderr,"Unable to open %s\n",filename);
      exit(1);
    }
  }
  else {
    mat_file = Mat_CreateVer(filename, NULL, MAT_FT_MAT5);
    if (mat_file == NULL) {
      fprintf(stderr,"Unable to create matlab file %s\n",filename);
      exit(1);
    }
  }

  /* word sizes */
  cpu_word_size=sizeof(double);
  io_word_size=0;

  /* open exodus file */
  exo_file=ex_open(argv[1],EX_READ,&cpu_word_size,&io_word_size,&exo_version);
  if (exo_file < 0){
    printf("error opening %s\n",argv[1]);
    exit(1);
  }

  /* print */
  fprintf(stderr,"translating %s to %s ...\n",argv[1],filename);

  /* read database paramters */
  line=(char *) calloc ((MAX_LINE_LENGTH+1),sizeof(char *));
  err = ex_get_init(exo_file,line,
	&num_axes,&num_nodes,&num_elements,&num_blocks,
        &num_node_sets,&num_side_sets);
  num_qa_lines   = ex_inquire_int(exo_file,EX_INQ_QA);
  num_info_lines = ex_inquire_int(exo_file,EX_INQ_INFO);
  num_time_steps = ex_inquire_int(exo_file,EX_INQ_TIME);
  err=ex_get_variable_param(exo_file,EX_GLOBAL,&num_global_vars);
  err=ex_get_variable_param(exo_file,EX_NODAL,&num_nodal_vars);
  err=ex_get_variable_param(exo_file,EX_ELEM_BLOCK,&num_element_vars);
  err=ex_get_variable_param(exo_file,EX_NODE_SET,&num_nodeset_vars);
  err=ex_get_variable_param(exo_file,EX_SIDE_SET,&num_sideset_vars);


  /* export paramters */
  PutInt("naxes",  1, 1,&num_axes);
  PutInt("nnodes", 1, 1,&num_nodes);
  PutInt("nelems", 1, 1,&num_elements);
  PutInt("nblks",  1, 1,&num_blocks);
  PutInt("nnsets", 1, 1,&num_node_sets);
  PutInt("nssets", 1, 1,&num_side_sets);
  PutInt("nsteps", 1, 1,&num_time_steps);
  PutInt("ngvars", 1, 1,&num_global_vars);
  PutInt("nnvars", 1, 1,&num_nodal_vars);
  PutInt("nevars", 1, 1,&num_element_vars);
  PutInt("nnsvars", 1, 1,&num_nodeset_vars);
  PutInt("nssvars", 1, 1,&num_sideset_vars);

  /* allocate -char- scratch space*/
  n =                              num_info_lines;
  n = (n > num_global_vars) ?  n : num_global_vars;
  n = (n > num_nodal_vars) ?   n : num_nodal_vars;
  n = (n > num_element_vars) ? n : num_element_vars;
  n = (n > num_blocks) ?       n : num_blocks;
  nstr2 = n;
  str2= (char **) calloc (n,sizeof(char *));
  for (i=0;i<nstr2;i++)
    str2[i]=(char *) calloc ((MAX_LINE_LENGTH+1),sizeof(char));
  str= (char *) calloc ((MAX_LINE_LENGTH+1)*n,sizeof(char *));

  /* title */
  PutStr("Title",line);

  /* QA records */
  if (num_qa_lines > 0 ){
    qa_records  =(char *(*)[4]) calloc (num_qa_lines*4,sizeof(char **));
    for (i=0;i<num_qa_lines;i++) 
      for (j=0;j<4;j++)
	qa_records[i][j]=(char *) calloc ((MAX_STR_LENGTH+1),sizeof(char));
    err=ex_get_qa(exo_file,qa_records);
    str[0]='\0';
    for (i=0;i<num_qa_lines;i++){
      for (j=0;j<4;j++)
	sprintf(str+strlen(str),"%s ",qa_records[i][j]);
      strcat(str,"\n");
    }
    for (i=0;i<num_qa_lines;i++){
        for (j=0;j<4;j++)
	  free(qa_records[i][j]);
    }
    free(qa_records);
  }

  /* information records */
  if (num_info_lines > 0 ){
    err = ex_get_info(exo_file,str2);
    str[0]='\0';
    for (i=0;i<num_info_lines;i++)
      sprintf(str+strlen(str),"%s\n",str2[i]);
    PutStr("info",str);
    str[0]='\0';
    for (i=0;i<num_info_lines;i++)
      if (strncmp(str2[i],"cavi",4)==0)
	sprintf(str+strlen(str),"%s\n",str2[i]);
    PutStr("cvxp",str);
  }

  /* nodal coordinates */
  x = (double *) calloc(num_nodes,sizeof(double));
  y = (double *) calloc(num_nodes,sizeof(double));
  if (num_axes == 3) 
    z = (double *) calloc(num_nodes,sizeof(double));
  else 
    z = NULL;
  err = ex_get_coord(exo_file,x,y,z);
  PutDbl("x0", num_nodes, 1, x);
  PutDbl("y0", num_nodes, 1, y);
  free(x);
  free(y);
  if (num_axes == 3){ 
    PutDbl("z0",num_nodes,1, z);
    free(z);
  }
  
   /* side sets */
  if(num_side_sets > 0){
    ids=(int *) calloc(num_side_sets,sizeof(int));
    err = ex_get_ids(exo_file,EX_SIDE_SET,ids);
    PutInt( "ssids",num_side_sets, 1,ids);
    nsssides = (int *) calloc(num_side_sets,sizeof(int)); /*dgriffi */
    nssdfac  = (int *) calloc(num_side_sets,sizeof(int)); /*dgriffi */
    for (i=0;i<num_side_sets;i++){
      err = ex_get_set_param(exo_file,EX_SIDE_SET, ids[i],&n1,&n2);
      nsssides[i]=n1; /* dgriffi */
      nssdfac[i]=n2;  /* dgriffi */
      /*
       * the following provision is from Version 1.6 when there are no
       * distribution factors in exodus file
       */
      has_ss_dfac = (n2 != 0);
      if(n2==0 || n1==n2){
	
	printf(" WARNING: Exodus II file does not contain distribution factors.\n");
	
	/* n1=number of faces, n2=number of df */
	/* using distribution factors to determine number of nodes in the sideset
         causes a lot grief since some codes do not output distribution factors
         if they are all equal to 1. mkbhard: I am using the function call below
         to figure out the total number of nodes in this sideset. Some redundancy
         exists, but it works for now */

	junk = (int*) calloc(n1,sizeof(int)); 
	err = ex_get_side_set_node_count(exo_file,ids[i],junk);
	n2=0; /* n2 will be equal to the total number of nodes in the sideset */
	for (j=0;j<n1;j++) n2+=junk[j];
	free(junk);

      }
	
      iscr = (int *) calloc(n1+n2,sizeof(int));
      err = ex_get_side_set_node_list(exo_file,ids[i],iscr,iscr+n1);
      /* number-of-nodes-per-side list */
      sprintf(str,"ssnum%02d",i+1);
      PutInt(str,n1,1,iscr); 
      /* nodes list */
      sprintf(str,"ssnod%02d",i+1);
      PutInt(str,n2,1,iscr+n1);
      free(iscr);
      /* distribution-factors list */
      scr = (double *) calloc (n2,sizeof(double));
      if (has_ss_dfac) {
	ex_get_side_set_dist_fact(exo_file,ids[i],scr);
      } else {
	for (j=0; j<n2; j++) {
	  scr[j] = 1.0;
	}
      }
      sprintf(str,"ssfac%02d",i+1);
      PutDbl(str,n2,1,scr);
      free(scr);
      /* element and side list for side sets (dgriffi) */
      elem_list = (int *) calloc(n1, sizeof(int));
      side_list = (int *) calloc(n1, sizeof(int));
      err = ex_get_set(exo_file,EX_SIDE_SET,ids[i],elem_list,side_list);
      sprintf(str,"ssside%02d",i+1);
      PutInt(str,n1,1,side_list);
      sprintf(str,"sselem%02d",i+1);
      PutInt(str,n1,1,elem_list);
      free(elem_list);
      free(side_list);

    }
    /* Store # sides and # dis. factors per side set (dgriffi) */
    PutInt("nsssides",num_side_sets,1,nsssides);
    PutInt("nssdfac",num_side_sets,1,nssdfac);
    free(ids);
    free(nsssides);
    free(nssdfac);
  }

  /* node sets (section by dgriffi) */
  if(num_node_sets > 0){
    ids=(int *) calloc(num_node_sets,sizeof(int));
    err = ex_get_ids(exo_file,EX_NODE_SET, ids);
    PutInt( "nsids",num_node_sets, 1,ids);
    nnsnodes = (int *) calloc(num_node_sets,sizeof(int)); 
    nnsdfac  = (int *) calloc(num_node_sets,sizeof(int));
    for (i=0;i<num_node_sets;i++){
      err = ex_get_set_param(exo_file,EX_NODE_SET,ids[i],&n1,&n2);
      iscr = (int *) calloc(n1,sizeof(int));
      err = ex_get_node_set(exo_file,ids[i],iscr);
      /* nodes list */
      sprintf(str,"nsnod%02d",i+1);
      PutInt(str,n1,1,iscr);
      free(iscr);
      /* distribution-factors list */
      scr = (double *) calloc (n2,sizeof(double));
      ex_get_node_set_dist_fact(exo_file,ids[i],scr);  
      sprintf(str,"nsfac%02d",i+1);
      PutDbl(str,n2,1,scr);
      free(scr);

      nnsnodes[i]=n1;
      nnsdfac[i]=n2;

    }

      /* Store # nodes and # dis. factors per node set */
      PutInt("nnsnodes",num_node_sets,1,nnsnodes);
      PutInt("nnsdfac",num_node_sets,1,nnsdfac);
      free(ids);
   
    free(nnsdfac);
    free(nnsnodes);
    
  }

  /* element blocks */
  ids=(int *) calloc(num_blocks,sizeof(int));
  num_elem_in_block=(int *) calloc(num_blocks,sizeof(int));
  err = ex_get_ids(exo_file,EX_ELEM_BLOCK,ids);
  PutInt( "blkids",num_blocks, 1,ids);
  for (i=0;i<num_blocks;i++) {
    err = ex_get_elem_block(exo_file,ids[i],str2[i],&n,&n1,&n2);
    num_elem_in_block[i]=n;
    iscr = (int *) calloc(n*n1,sizeof(int));
    err = ex_get_conn(exo_file,EX_ELEM_BLOCK,ids[i],iscr, NULL, NULL);
    sprintf(str,"blk%02d",i+1);
    PutInt(str,n1,n,iscr);
    free(iscr);
  }
  str[0]='\0';
  for (i=0;i<num_blocks;i++)
    sprintf(str+strlen(str),"%s\n",str2[i]);
  PutStr("blknames",str);  

  /* time values */
  if (num_time_steps > 0 ) {
    scr = (double *) calloc (num_time_steps,sizeof(double));
    err= ex_get_all_times (exo_file,scr);
    PutDbl( "time", num_time_steps, 1,scr);
    free(scr); 
  }

  /* global variables */
  if (num_global_vars > 0 ) {
    err = ex_get_variable_names(exo_file,EX_GLOBAL,num_global_vars,str2);
    str[0]='\0';
    for (i=0;i<num_global_vars;i++)
      sprintf(str+strlen(str),"%s\n",str2[i]);
    PutStr("gnames",str);
    scr = (double *) calloc (num_time_steps,sizeof(double));
    for (i=0;i<num_global_vars;i++){
      sprintf(str,"gvar%02d",i+1);
      err=ex_get_glob_var_time(exo_file,i+1,1,num_time_steps,scr);
      PutDbl(str,num_time_steps,1,scr);
    }
    free(scr);
  }

  /* nodal variables */
  if (num_nodal_vars > 0 ) {
    err = ex_get_variable_names(exo_file,EX_NODAL,num_nodal_vars,str2);
    str[0]='\0';
    for (i=0;i<num_nodal_vars;i++)
      sprintf(str+strlen(str),"%s\n",str2[i]);
    PutStr("nnames",str);
    scr = (double *) calloc (num_nodes*num_time_steps,sizeof(double));
    for (i=0;i<num_nodal_vars;i++){
      sprintf(str,"nvar%02d",i+1);
      for (j=0;j<num_time_steps;j++)
	err=ex_get_nodal_var(exo_file,j+1,i+1,num_nodes,
                                  scr+num_nodes*j);
      PutDbl(str,num_nodes,num_time_steps,scr);
    }
    free(scr);
  }

  /* element variables */
  if (num_element_vars > 0 ) {
    err = ex_get_variable_names(exo_file,EX_ELEM_BLOCK,num_element_vars,str2);
    str[0]='\0';
    for (i=0;i<num_element_vars;i++)
      sprintf(str+strlen(str),"%s\n",str2[i]);
    PutStr("enames",str);
    /* truth table */
    iscr = (int *) calloc(num_element_vars*num_blocks, sizeof(int));
    ex_get_elem_var_tab(exo_file,num_blocks,num_element_vars,iscr);
    for (i=0;i<num_element_vars;i++){
      scr = (double *) calloc (num_elements*num_time_steps,sizeof(double));
      n=0;
      sprintf(str,"evar%02d",i+1);
      for (j=0;j<num_time_steps;j++){
	for (k=0;k<num_blocks;k++){ 
          if(iscr[num_element_vars*k+i]==1)
	      ex_get_elem_var(exo_file,j+1,i+1,ids[k],num_elem_in_block[k],scr+n);
	      n=n+num_elem_in_block[k];
	      
	}
      }
      PutDbl(str,num_elements,num_time_steps,scr);
      free(scr);
    }
    free(iscr);
  }
  free(num_elem_in_block);
  free(ids);
 
  /* node and element number maps */
  ex_opts(0);  /* turn off error reporting. It is not an error to have no map*/
  ids = (int *)malloc(num_nodes*sizeof(int));
  err = ex_get_node_num_map(exo_file,ids);
  if ( err==0 ){
    PutInt("node_num_map",num_nodes,1,ids);
  }
  free(ids);

  ids = (int *)malloc(num_elements*sizeof(int));
  err = ex_get_elem_num_map(exo_file,ids);
  if ( err==0 ){
    PutInt("elem_num_map",num_elements,1,ids);
  }
  free(ids);


  /* close exo file */
  ex_close(exo_file);
  
  /* close mat file */
  if ( textfile )
    fclose(m_file);
  else
    Mat_Close(mat_file);

  /* */
  fprintf(stderr,"done.\n");

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
