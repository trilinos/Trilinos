#include<stdio.h>
/* converts a Sundance-style matrix file to an aztec-style matrix file */


int convertSund2AZdatafile( char *infile, char *outfile)
{
	FILE *infp, *outfp;
	int row, col, nrows, prevrow;
	double val;

	nrows=0;  

	infp=fopen(infile, "r");
	outfp=fopen(outfile, "w");

	if ((infp==NULL) || (outfp==NULL))
		printf("convertSund2AZdatafile couldn't open one of the files\n");

	while (fscanf(infp,"%d %d %lf", &row, &col, &val) != EOF)
		nrows=row+1; /* zero-based array, so need to add 1 */
 	
	fprintf(outfp, "%d\n", nrows);

	rewind(infp);

	prevrow=0;

	while (fscanf(infp,"%d %d %lf", &row, &col, &val) != EOF) {
		if (row !=prevrow)
			fprintf(outfp,"  -1");
		if (val != 0.0) 
			fprintf(outfp,"\n%d   %g", col, val);
		prevrow=row;
	}

	fprintf(outfp,"  -1\n");

	fclose(infp);
	fclose(outfp);
	return (nrows);
}




main()
{
	char infile[80], outfile[80];
	int ret;
	 
	sprintf(infile, "/home/people/dawn/Software/Sundance/examples/Amatrix.dat\0");
	sprintf(outfile, "AmatrixAZ.dat\0");
	
	ret=convertSund2AZdatafile(infile, outfile);

	sprintf(infile, "/home/people/dawn/Software/Sundance/examples/Mpmatrix.dat\0");
	sprintf(outfile, "MpmatrixAZ.dat\0");
	
	ret=convertSund2AZdatafile(infile, outfile);
	
}





