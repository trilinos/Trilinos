/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <limits>
#include <cmath>
#include <stdexcept>

#include <snl_fei_Utils.hpp>
#include "fei_Record.hpp"
#include <fei_MatrixGraph.hpp>
#include <fei_SparseRowGraph.hpp>
#include <fei_Matrix_Impl.hpp>
#include <fei_ParameterSet.hpp>

#include <fei_CommUtils.hpp>
#include <fei_TemplateUtils.hpp>
#include <fei_CSVec.hpp>
#include <fei_chk_mpi.hpp>

#undef fei_file
#define fei_file "snl_fei_Utils.cpp"
#include <fei_ErrMacros.hpp>

//----------------------------------------------------------------------------
const char* snl_fei::getParam(const char* key,
			      int numParams,
			      const char* const* paramStrings)
{
  int foundOffset = -1;
  return( getParam(key, numParams, paramStrings, foundOffset) );
}

//----------------------------------------------------------------------------
const char* snl_fei::getParamValue(const char* key,
				   int numParams,
				   const char* const* paramStrings,
				   char separator)
{
  const char* param = getParam(key, numParams, paramStrings);

  return( param != NULL ? skipSeparator(param, separator) : NULL );
}

//----------------------------------------------------------------------------
const char* snl_fei::getParamValue(const char* key,
				   int numParams,
				   const char* const* paramStrings,
				   int& foundOffset,
				   char separator)
{
  const char* param = getParam(key, numParams, paramStrings, foundOffset);

  return( param != NULL ? skipSeparator(param, separator) : NULL );
}

//----------------------------------------------------------------------------
int snl_fei::getIntParamValue(const char* key,
			      int numParams,
			      const char*const* params,
			      int& paramValue)
{
  const char* tempstr = getParamValue(key, numParams, params);
  if (tempstr==NULL) return(-1);

  std::string strg(tempstr);
  FEI_ISTRINGSTREAM isstr(strg);
  isstr >> paramValue;
  return( isstr.fail() ? -1 : 0 );
}

//----------------------------------------------------------------------------
int snl_fei::getDoubleParamValue(const char* key,
				 int numParams,
				 const char*const* params,
				 double& paramValue)
{
  const char* tempstr = getParamValue(key, numParams, params);
  if (tempstr==NULL) return(-1);

  std::string strg(tempstr);
  FEI_ISTRINGSTREAM isstr(strg);
  isstr >> paramValue;
  return( isstr.fail() ? -1 : 0 );
}

//----------------------------------------------------------------------------
int snl_fei::getDoubleParamValue(const char* key,
				 std::vector<std::string>& params,
				 double& paramValue)
{
  const char* tempstr = getParamValue(key, params);
  if (tempstr==NULL) return(-1);

  std::string strg(tempstr);
  FEI_ISTRINGSTREAM isstr(strg);
  isstr >> paramValue;
  return( isstr.fail() ? -1 : 0 );
}

//----------------------------------------------------------------------------
const char* snl_fei::getParam(const char* key,
			      int numParams,
			      const char* const* paramStrings,
			      int& foundOffset)
{
  const char* returnPtr = NULL;
  foundOffset = -1;

  //if the search-key is null, or if the list of strings to be searched is null,
  //or if the number of strings to be searched is null, then return now.
  //
  if (key == NULL || paramStrings == NULL || numParams == 0) {
    return(returnPtr);
  }

  unsigned keyLen = strlen(key);

  for(int i=numParams-1; i>=0; --i) {
    const char* paramStr = paramStrings[i];

    //if the i-th paramString is null, skip it.
    if (paramStr == NULL) continue;

    unsigned paramStr_len = leading_substring_length(paramStr);
    if (paramStr_len != keyLen) continue;

    if (strncmp(key, paramStr, keyLen) == 0) {
      returnPtr = paramStr;
      foundOffset = i;
      return(returnPtr);
    }
  }

  return(returnPtr);
}

//----------------------------------------------------------------------------
const char* snl_fei::getParam(const char* key,
			      std::vector<std::string>& params,
			      int& foundOffset)
{
  std::vector<std::string>::iterator
    p_iter = params.begin(),
    p_end = params.end();

  int offset = 0;
  for(; p_iter != p_end; ++p_iter, ++offset) {
    std::string& pstr = *p_iter;
    int ssize = pstr.size();
    if (ssize < 1) continue;

    std::string::size_type i = pstr.find(key);
    if (i == 0) {
      foundOffset = offset;
      return( pstr.c_str() );
    }
  }

  return( NULL );
}

//----------------------------------------------------------------------------
const char* snl_fei::getParamValue(const char* key,
				   std::vector<std::string>& params,
				   char separator)
{
  int offset = 0;
  return( skipSeparator( getParam(key, params, offset), separator) );
}

//----------------------------------------------------------------------------
void snl_fei::separate_string(const char* input_string,
			      const char* substring,
			      const char*& before_substring,
			      int& len_before_substring,
			      const char*& after_substring,
			      int& len_after_substring)
{
  if (input_string == NULL) {
    before_substring = NULL;
    len_before_substring = 0;
    after_substring = NULL;
    len_after_substring = 0;
    return;
  }

  int len_input_string = strlen(input_string);
  before_substring = input_string;

  if (substring == NULL) {
    len_before_substring = len_input_string;
    after_substring = NULL;
    len_after_substring = 0;
    return;
  }

  int len_substring = strlen(substring);

  const char* s1 = strstr(input_string, substring);
  if (s1 == NULL) {
    len_before_substring = len_input_string;
    after_substring = NULL;
    len_after_substring = 0;
    return;
  }

  after_substring = skipSeparator(s1, substring[len_substring-1]);
  len_before_substring = s1 - input_string;
  len_after_substring = len_input_string - len_before_substring - len_substring;
}

//----------------------------------------------------------------------------
int snl_fei::storeNamedAttribute(const char* name,
				 void* attribute,
				 std::vector<char*>& attributeNames,
				 std::vector<void*>& attributes)
{
  int offset = -1;
  const char* foundname = attributeNames.empty() ?
    0 : snl_fei::getParam(name, attributeNames.size(),
			  &(attributeNames[0]), offset);

  if (foundname != 0) {
    attributes[offset] = attribute;
  }
  else {
    char* newname = new char[strlen(name)+1];
    strcpy(newname, name);
    attributeNames.push_back(newname);
    attributes.push_back(attribute);
  }

  return(0);
}

//----------------------------------------------------------------------------
void* snl_fei::retrieveNamedAttribute(const char* name,
				      std::vector<char*>& attributeNames,
				      std::vector<void*>& attributes)
{
  int offset = -1;
  const char* foundname = attributeNames.empty() ?
    0 : snl_fei::getParam(name, attributeNames.size(),
			  &(attributeNames[0]), offset);

  void* returnVal = NULL;

  if (foundname != 0) {
    returnVal = attributes[offset];
  }
  else {
    return( returnVal );
  }

  return( returnVal );
}

//----------------------------------------------------------------------------
unsigned snl_fei::leading_substring_length(const char* string)
{
  if (string == NULL) {
    return(0);
  }

  const char* lastchar = string;
  unsigned len = 0;

  while(*lastchar != ' ' && *lastchar != '\t' && *lastchar != '\0') {
    ++len;
    ++lastchar;
  }

  return(len);
}

//----------------------------------------------------------------------------
const char* snl_fei::skipSeparator(const char* paramString,
				   char separator)
{
  if (paramString == NULL) return(NULL);

  const char* result = strchr(paramString, separator);

  if (result == NULL) return(result);

  //allow for the possibility that separator is repeated
  while(result[0] == separator) result++;

  return( result );
}

//----------------------------------------------------------------------------
int snl_fei::mergeStringLists(char**& strings,
			      int& numStrings,
			      const char*const* stringsToMerge,
			      int numStringsToMerge)
{
  int i;
  if (numStrings == 0) {
    strings = new char*[numStringsToMerge];

    for(i=0; i<numStringsToMerge; ++i){
      strings[i] = new char[strlen(stringsToMerge[i])+1];

      strcpy(strings[i], stringsToMerge[i]);
    }

    numStrings = numStringsToMerge;
  }
  else {
    int numMerged = 0;
    bool* merged = new bool[numStringsToMerge];
    for(i=0; i<numStringsToMerge; ++i) {
      const char* temp;
      const char* temp2;
      int len_temp, len_temp2;
      separate_string(stringsToMerge[i], " ",
                      temp, len_temp, temp2, len_temp2);
      std::string strtemp(temp, len_temp);
      int foundOffset = -1;
      const char* matchingString = 
	snl_fei::getParam(strtemp.c_str(), numStrings, strings, foundOffset);
      if (matchingString != NULL) {
        int len = strlen(stringsToMerge[i])+1;
	delete [] strings[foundOffset];
	strings[foundOffset] = new char[len];
	strcpy(strings[foundOffset], stringsToMerge[i]);
	merged[i] = true;
	++numMerged;
      }
      else merged[i] = false;
    }

    if (numMerged == numStringsToMerge) {
      delete [] merged;
      return(0);
    }

    char** newStrings = new char*[numStrings+numStringsToMerge-numMerged];
    for(i=0; i<numStrings; ++i) {
      newStrings[i] = strings[i];
    }
    int offset = numStrings;
    for(i=0; i<numStringsToMerge; ++i) {
      if (!merged[i]) {
	int len = strlen(stringsToMerge[i])+1;
	newStrings[offset] = new char[len];
	strcpy(newStrings[offset++], stringsToMerge[i]);
      }
    }
    delete [] merged;

    delete [] strings;
    strings = newStrings;
    numStrings+= numStringsToMerge-numMerged;
  }

  return(0);
}

//----------------------------------------------------------------------------
int snl_fei::resolveConflictingCRs(fei::MatrixGraph& matrixGraph,
				   fei::Matrix& bcEqns,
                                   const std::vector<int>& bcEqnNumbers)
{
  int numLagrangeConstraints = matrixGraph.getLocalNumLagrangeConstraints();
  if (numLagrangeConstraints < 1) {
    return(0);
  }

  double coefs[3];
  int indices[3];
  indices[0] = 0;
  indices[1] = 1;
  indices[2] = 2;
  coefs[0] = 1.0;
  coefs[1] = 0.0;
  coefs[2] = 1.0;

  double* coefPtr = &(coefs[0]);
  int* indicesPtr = &(indices[0]);

  double fei_eps = 1.e-49;

  std::vector<int> cr_indices;
  std::map<int,Constraint<fei::Record<int>*>*>& lagrangeConstraints =
    matrixGraph.getLagrangeConstraints();

  std::map<int,Constraint<fei::Record<int>*>*>::const_iterator
    cr_iter = lagrangeConstraints.begin(),
    cr_end  = lagrangeConstraints.end();

  while(cr_iter != cr_end) {
    Constraint<fei::Record<int>*>* cr = (*cr_iter).second;

    CHK_ERR( matrixGraph.getConstraintConnectivityIndices(cr, cr_indices) );

    std::vector<double>& weights = cr->getMasterWeights();
    double* weightsPtr = &weights[0];

    int len = cr_indices.size();
    int* cr_indPtr = &(cr_indices[0]);
    for(int j=0; j<len; ++j) {
      if (std::abs(weightsPtr[j] + 1.0) > fei_eps) continue;

      int eqn = cr_indPtr[j];
      if (fei::binarySearch(eqn, bcEqnNumbers) > -1) {
	int cr_eqn = cr->getEqnNumber();

	CHK_ERR( bcEqns.copyIn(1, &cr_eqn, 3, indicesPtr,
			       (double**)(&coefPtr)) );
      }
    }
    ++cr_iter;
  }

  return(0);
}

//----------------------------------------------------------------------------
int snl_fei::gatherRemoteEssBCs(fei::CSVec& essBCs,
				fei::SparseRowGraph* remoteGraph,
				fei::Matrix& matrix)
{
  std::vector<int>& rrowOffs = remoteGraph->rowOffsets;
  std::vector<int>& rcols = remoteGraph->packedColumnIndices;

  int numEssEqns = essBCs.size();
  if (numEssEqns > 0) {
    int* essEqns = &(essBCs.indices()[0]);
    double* coefs = &(essBCs.coefs()[0]);

    if (rrowOffs.size() > 0 && rcols.size() > 0) {

      int* rowOffsPtr = &(rrowOffs[0]);
      int* rcolsPtr = &(rcols[0]);

      for(int j=0; j<numEssEqns; ++j) {

        int eqn = essEqns[j];

        for(unsigned i=0; i<rrowOffs.size()-1; ++i) {
          int len = rowOffsPtr[i+1]-rowOffsPtr[i];
          int* colsPtr = &(rcolsPtr[rowOffsPtr[i]]);

          if (fei::binarySearch(eqn, colsPtr, len) > -1) {
            double coef = coefs[j];
            double* coefPtr = &coef;

            for(int k=0; k<len; ++k) {
              CHK_ERR( matrix.copyIn(1, &(colsPtr[k]),
                       1, &(eqn), &coefPtr) );
            }
          }
        }
      }
    }
  }

  CHK_ERR( matrix.gatherFromOverlap(false) );

  return(0);
}

fei::SharedPtr<fei::SparseRowGraph>
snl_fei::mergeSparseRowGraphs(const fei::SparseRowGraph* srg1,
                              const fei::SparseRowGraph* srg2)
{
  fei::SharedPtr<fei::SparseRowGraph> newgraph(new fei::SparseRowGraph);

  int numrows = srg1->rowNumbers.size() + srg2->rowNumbers.size();
  int nnz = srg1->packedColumnIndices.size() + srg2->packedColumnIndices.size();

  newgraph->rowNumbers.resize(numrows);
  newgraph->rowOffsets.resize(numrows+1);
  newgraph->packedColumnIndices.resize(nnz);

  std::map<int,int> rowmap;

  for(unsigned i=0; i<srg1->rowNumbers.size(); ++i) {
    int rowlen = srg1->rowOffsets[i+1]-srg1->rowOffsets[i];
    rowmap.insert(std::make_pair(srg1->rowNumbers[i], rowlen));
  }

  for(unsigned i=0; i<srg2->rowNumbers.size(); ++i) {
    int rowlen = srg2->rowOffsets[i+1]-srg2->rowOffsets[i];
    rowmap.insert(std::make_pair(srg2->rowNumbers[i], rowlen));
  }

  std::map<int,int>::iterator
    r_iter = rowmap.begin(),
    r_end  = rowmap.end();

  int offset = 0;
  for(unsigned i=0; r_iter != r_end; ++r_iter, ++i) {
    newgraph->rowNumbers[i] = r_iter->first;
    int rowlen = r_iter->second;
    newgraph->rowOffsets[i] = offset;
    offset += rowlen;
    r_iter->second = i;
  }
  newgraph->rowOffsets[numrows] = offset;

  for(unsigned i=0; i<srg1->rowNumbers.size(); ++i) {
    r_iter = rowmap.find(srg1->rowNumbers[i]);

    int newcoloffset = newgraph->rowOffsets[r_iter->second];

    int jbegin = srg1->rowOffsets[i];
    int jend   = srg1->rowOffsets[i+1];
    for(int j=jbegin; j<jend; ++j) {
      newgraph->packedColumnIndices[newcoloffset++] =
        srg1->packedColumnIndices[j];
    }
  }

  for(unsigned i=0; i<srg2->rowNumbers.size(); ++i) {
    r_iter = rowmap.find(srg2->rowNumbers[i]);

    int newcoloffset = newgraph->rowOffsets[r_iter->second];

    int jbegin = srg2->rowOffsets[i];
    int jend   = srg2->rowOffsets[i+1];
    for(int j=jbegin; j<jend; ++j) {
      newgraph->packedColumnIndices[newcoloffset++] =
        srg2->packedColumnIndices[j];
    }
  }

  return(newgraph);
}

void snl_fei::copy2DBlockDiagToColumnContig(int numBlocks,
					    const int* blockSizes,
					    const double*const* values2d,
					    int format,
					    double* colcontigvalues)
{
  int i, j, k, offset, coffset;

  switch(format) {
  case FEI_BLOCK_DIAGONAL_ROW:
    offset = 0;
    coffset = 0;
    for(i=0; i<numBlocks; ++i) {
      int bsize = blockSizes[i];
      for(j=0; j<bsize; ++j) {
	for(k=0; k<bsize; ++k) {
	  colcontigvalues[coffset+j*bsize+k] = values2d[offset+j][k];
	}
      }
      offset += bsize;
      coffset += bsize*bsize;
    }
    break;
  default:
    FEI_OSTRINGSTREAM osstr;
    osstr << "copy2DBlockDiagToColumnContig ERROR, format="<<format
	  << " not recognized.";
    throw std::runtime_error(osstr.str());
  }
}

void snl_fei::copy2DToColumnContig(int numrows,
				   int numcols,
				   const double*const* values2d,
				   int format,
				   double* colcontigvalues)
{
  int i, j;

  switch(format) {

  case FEI_DENSE_ROW:
    for(i=0; i<numrows; ++i) {
      for(j=0; j<numcols; ++j) {
	colcontigvalues[j*numrows+i] = values2d[i][j];
      }
    }
    break;

  case FEI_DENSE_COL:
    for(j=0; j<numcols; ++j) {
      for(i=0; i<numrows; ++i) {
	colcontigvalues[j*numrows+i] = values2d[j][i];
      }
    }
    break;

  default:
    FEI_OSTRINGSTREAM osstr;
    osstr << "copy2DToColumnContig ERROR, format="<<format<<" not recognized.";
    throw std::runtime_error(osstr.str());
  }
}
