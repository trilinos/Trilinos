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

/**************** conditional compilation ************************/
#ifdef HAS_MATLAB
#include "mat.h"
#define EXT ".mat"
int textfile=0;
#else


/* no matlab libraries */
#define MATFile FILE
#define matOpen fopen
#define matClose fclose
#define matPutStr mPutStr
#define matPutDbl mPutDbl
#define EXT ".m"
int textfile=1;
#endif
/************* end conditional compilation ***********************/

FILE* m_file=0;       /* file for m file output */
MATFile *mat_file=0;  /* file for binary .mat output */

extern void add_to_log(const char *my_name);

static char *qainfo[] =
{
  "exo2mat",
  "2010/09/28",
  "1.11",
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




/*----------------------*/
#ifdef HAS_MATLAB

/* put string in mat file*/
void matPutStr (char *name,char *str)
{
  mxArray *pa;
  pa=mxCreateString(str);

  /* MATLAB 5: 
   * mxSetName(pa,name);
   * matPutArray(mat_file,pa);
   */

  /* MATLAB 6: */
  matPutVariable(mat_file,name,pa);

  mxDestroyArray(pa);
}

/* put double in mat file*/
void matPutDbl (char *name,int n1,int n2,double *pd)
{
  mxArray *pa;
  pa=mxCreateDoubleMatrix(n1,n2,mxREAL);
  memcpy(mxGetPr(pa),pd,n1*n2*sizeof(double));

  /* MATLAB 5 */
  /* mxSetName(pa,name); */
  /* matPutArray(mat_file,pa); */
  
  /* MATLAB 6: */
  matPutVariable(mat_file,name,pa);

  mxDestroyArray(pa);
}
#endif
/*----------------------*/

/* wrappers for the two output routine types */
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
    


/* int-to-double copy (since MATLAB needs double)*/
double *i2d (int *iscr, int n)
{
  int i;
  double *scr;
  scr = (double *) calloc(n,sizeof(double));
  for (i=0;i<n;i++)
    scr[i]=iscr[i];
  return scr;
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
    *str,**str2,c[1],*(*qa_records)[4],*line, *oname, *dot, *filename;

  const char* ext=EXT;

  int   
    i,j,k,n,n1,n2,cpu_word_size,io_word_size,exo_file,err,
    num_axes,num_nodes,num_elements,num_blocks,
    num_side_sets,num_node_sets,num_time_steps,
    num_qa_lines,num_info_lines,num_global_vars,
    num_nodal_vars,num_element_vars,*ids,*iscr,*num_elem_in_block,*junk,
    *elem_list,*side_list,
    *nsssides,*nssdfac,
    *nnsnodes,*nnsdfac,
    nstr2, has_ss_dfac;
    
  float
    exo_version;

  double
    f,*scr,*x,*y,*z;

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
  else 
    if ( (mat_file = matOpen(filename, "w"))==0 ){
      fprintf(stderr,"Unable to open %s\n",filename);
      exit(1);
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
  err=ex_inquire(exo_file,EX_INQ_QA,&num_qa_lines,&f,c);
  err=ex_inquire(exo_file,EX_INQ_INFO,&num_info_lines,&f,c);
  err=ex_inquire(exo_file,EX_INQ_TIME,&num_time_steps,&f,c);
  err=ex_get_var_param(exo_file,"g",&num_global_vars);
  err=ex_get_var_param(exo_file,"n",&num_nodal_vars);
  err=ex_get_var_param(exo_file,"e",&num_element_vars);


  /* export paramters */
  f=num_axes         ;PutDbl("naxes",  1, 1,&f);
  f=num_nodes        ;PutDbl("nnodes", 1, 1,&f);
  f=num_elements     ;PutDbl("nelems", 1, 1,&f);
  f=num_blocks       ;PutDbl("nblks",  1, 1,&f);
  f=num_node_sets    ;PutDbl("nnsets", 1, 1,&f);
  f=num_side_sets    ;PutDbl("nssets", 1, 1,&f);
  f=num_time_steps   ;PutDbl("nsteps", 1, 1,&f);
  f=num_global_vars  ;PutDbl("ngvars", 1, 1,&f);
  f=num_nodal_vars   ;PutDbl("nnvars", 1, 1,&f);
  f=num_element_vars ;PutDbl("nevars", 1, 1,&f);

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
    err = ex_get_side_set_ids(exo_file,ids);
    scr = i2d(ids,num_side_sets);
    PutDbl( "ssids",num_side_sets, 1,scr);
    free(scr);
    nsssides = (int *) calloc(num_side_sets,sizeof(int)); /*dgriffi */
    nssdfac  = (int *) calloc(num_side_sets,sizeof(int)); /*dgriffi */
    for (i=0;i<num_side_sets;i++){
      err = ex_get_side_set_param(exo_file,ids[i],&n1,&n2);
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
	err = ex_get_side_set_param(exo_file,ids[i],&n1,&n2);
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
      scr = i2d(iscr,n1);
      PutDbl(str,n1,1,scr); 
      free(scr);
      /* nodes list */
      sprintf(str,"ssnod%02d",i+1);
      scr = i2d(iscr+n1,n2);
      PutDbl(str,n2,1,scr);
      free(scr);
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
      err = ex_get_side_set(exo_file,ids[i],elem_list,side_list);
      sprintf(str,"ssside%02d",i+1);
      scr =i2d(side_list,n1);
      PutDbl(str,n1,1,scr);
      free(scr);
      sprintf(str,"sselem%02d",i+1);
      scr =i2d(elem_list,n1);
      PutDbl(str,n1,1,scr);
      free(scr);
      free(elem_list);
      free(side_list);

    }
    /* Store # sides and # dis. factors per side set (dgriffi) */
    scr = i2d(nsssides,num_side_sets);
    PutDbl("nsssides",num_side_sets,1,scr);
    free(scr);
    scr = i2d(nssdfac,num_side_sets);
    PutDbl("nssdfac",num_side_sets,1,scr);
    free(scr);
    free(ids);
    free(nsssides);
    free(nssdfac);
  }

  /* node sets (section by dgriffi) */
  if(num_node_sets > 0){
    ids=(int *) calloc(num_node_sets,sizeof(int));
    err = ex_get_node_set_ids(exo_file,ids);
    scr = i2d(ids,num_node_sets);
    PutDbl( "nsids",num_node_sets, 1,scr);
    free(scr);
    nnsnodes = (int *) calloc(num_node_sets,sizeof(int)); 
    nnsdfac  = (int *) calloc(num_node_sets,sizeof(int));
    for (i=0;i<num_node_sets;i++){
      err = ex_get_node_set_param(exo_file,ids[i],&n1,&n2);
      iscr = (int *) calloc(n1,sizeof(int));
      err = ex_get_node_set(exo_file,ids[i],iscr);
      /* nodes list */
      sprintf(str,"nsnod%02d",i+1);
      scr = i2d(iscr,n1);
      PutDbl(str,n1,1,scr);
      free(scr);
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
      scr = i2d(nnsnodes,num_node_sets);
      PutDbl("nnsnodes",num_node_sets,1,scr);
      free(scr);
      scr = i2d(nnsdfac,num_node_sets);
      PutDbl("nnsdfac",num_node_sets,1,scr);
      free(scr);
      free(ids);
   
    free(nnsdfac);
    free(nnsnodes);
    
  }

  /* element blocks */
  ids=(int *) calloc(num_blocks,sizeof(int));
  num_elem_in_block=(int *) calloc(num_blocks,sizeof(int));
  err = ex_get_elem_blk_ids(exo_file,ids);
  scr = i2d(ids,num_blocks);
  PutDbl( "blkids",num_blocks, 1,scr);
  free(scr);
  for (i=0;i<num_blocks;i++) {
    err = ex_get_elem_block(exo_file,ids[i],str2[i],&n,&n1,&n2);
    num_elem_in_block[i]=n;
    iscr = (int *) calloc(n*n1,sizeof(int));
    err = ex_get_elem_conn(exo_file,ids[i],iscr);
    sprintf(str,"blk%02d",i+1);
    scr = i2d(iscr,n*n1);
    PutDbl(str,n1,n,scr);
    free(iscr);
    free(scr);
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
    err = ex_get_var_names(exo_file,"g",num_global_vars,str2);
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
    err = ex_get_var_names(exo_file,"n",num_nodal_vars,str2);
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
    err = ex_get_var_names(exo_file,"e",num_element_vars,str2);
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
    scr = i2d( ids, num_nodes );
    PutDbl("node_num_map",num_nodes,1,scr);
    free(scr);
  }
  free(ids);

  ids = (int *)malloc(num_elements*sizeof(int));
  err = ex_get_elem_num_map(exo_file,ids);
  if ( err==0 ){
    scr = i2d( ids, num_elements );
    PutDbl("elem_num_map",num_elements,1,scr);
    free(scr);
  }
  free(ids);


  /* close exo file */
  ex_close(exo_file);
  
  /* close mat file */
  if ( textfile )
    fclose(m_file);
  else
    matClose(mat_file);

  /* */
  fprintf(stderr,"done.\n");

  free(filename);
  free(line);
  
  free(str);
  for (i=0;i<nstr2;i++)
    free(str2[i]);
  free(str2);
  

  /* exit status */
  add_to_log("exo2mat");
  return(0);
}
