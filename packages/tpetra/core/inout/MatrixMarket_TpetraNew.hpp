// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __MATRIXMARKET_TPETRA_NEW_HPP
#define __MATRIXMARKET_TPETRA_NEW_HPP

///////////////////////////////////////////////////////////////////////////
// Functions to read a MatrixMarket file and load it into a matrix.
// Adapted from 
//    Trilinos/packages/epetraext/src/inout/EpetraExt_CrsMatrixIn.cpp
// Modified by Jon Berry and Karen Devine to make matrix reallocations 
// more efficient; rather than insert each non-zero separately, we 
// collect rows of non-zeros for insertion.
// Modified by Karen Devine and Steve Plimpton to prevent all processors 
// from reading the same file at the same time (ouch!); chunks of the file 
// are read and broadcast by processor zero; each processor then extracts 
// its nonzeros from the chunk, sorts them by row, and inserts them into 
// the matrix.
// The variable "chunkSize" can be changed to modify the size of the chunks
// read from the file.  Larger values of chunkSize lead to faster execution
// and greater memory use.
///////////////////////////////////////////////////////////////////////////

private:

template <typename global_ordinal_type, typename scalar_type>
static 
Teuchos::RCP<Distribution<global_ordinal_type,scalar_type> > 
buildDistribution(
  std::string &distribution,  // string indicating whether to use 1D, 2D or
                              // file-based matrix distribution
  size_t nRow,                // Number of global matrix rows
  size_t nCol,                // Number of global matrix rows
  const Teuchos::ParameterList &params, // Parameters to the file reading
  const Teuchos::RCP<const Teuchos::Comm<int> > &comm
                            // communicator to be used in maps
)
{
  // Function to build the sets of GIDs for 1D or 2D matrix distribution.
  // Output is a Distribution object.

  int me = comm->getRank();

  using basedist_t = Distribution<global_ordinal_type,scalar_type>;
  basedist_t *retval;

  bool randomize = false;   // Randomly permute rows and columns before storing
  {
  const Teuchos::ParameterEntry *pe = params.getEntryPtr("randomize");
  if (pe != NULL) 
    randomize = pe->getValue<bool>(&randomize);
  }

  std::string partitionFile = "";  // File to reassign rows to parts
  {
  const Teuchos::ParameterEntry *pe = params.getEntryPtr("partitionFile");
  if (pe != NULL) 
    partitionFile = pe->getValue<std::string>(&partitionFile);
  }

  if (distribution == "2D") {   // Generate 2D distribution
    if (partitionFile != "") {
      // Generate 2D distribution with vector partition file
      TEUCHOS_TEST_FOR_EXCEPTION(randomize, std::logic_error,
              "Randomization not available for 2DVec distributions.");
      Distribution2DVec<global_ordinal_type,scalar_type> *dist =
              new Distribution2DVec<global_ordinal_type, scalar_type>(
                                    nRow, comm, params, partitionFile);
      retval = dynamic_cast<basedist_t *>(dist);
    }
    else if (randomize) {
      // Randomly assign rows/columns to processors
      Distribution2DRandom<global_ordinal_type,scalar_type> *dist =
              new Distribution2DRandom<global_ordinal_type,scalar_type>(
                                       nRow, comm, params);
      retval = dynamic_cast<basedist_t *>(dist);
    }
    else {
      // The vector will be distributed linearly, with extra vector entries
      // spread among processors to maintain balance in rows and columns.
      Distribution2DLinear<global_ordinal_type,scalar_type> *dist =
              new Distribution2DLinear<global_ordinal_type,scalar_type>(
                                       nRow, comm, params);
      retval = dynamic_cast<basedist_t *>(dist);
    }
  }
  else if (distribution == "1D") {
    if (partitionFile != "") {
      // Generate 1D row-based distribution with vector partition file
      Distribution1DVec<global_ordinal_type,scalar_type> *dist =
              new Distribution1DVec<global_ordinal_type,scalar_type>(
                                    nRow, comm, params, partitionFile);
      retval = dynamic_cast<basedist_t*>(dist);
    }
    else if (randomize) {
      // Randomly assign rows to processors.
      Distribution1DRandom<global_ordinal_type,scalar_type> *dist =
              new Distribution1DRandom<global_ordinal_type,scalar_type>(
                                       nRow, comm, params);
      retval = dynamic_cast<basedist_t *>(dist);
    }
    else {
      // Linear map similar to Trilinos default.
      Distribution1DLinear<global_ordinal_type,scalar_type> *dist =
              new Distribution1DLinear<global_ordinal_type,scalar_type>(
                                       nRow, comm, params);
      retval = dynamic_cast<basedist_t *>(dist);
    }
  }
  else if (distribution == "LowerTriangularBlock") {
    if (randomize && me == 0) {
      std::cout << "WARNING: Randomization not available for "
                << "LowerTriangularBlock distributions." << std::endl;
    }

    DistributionLowerTriangularBlock<global_ordinal_type,scalar_type> *dist = 
          new DistributionLowerTriangularBlock<global_ordinal_type,scalar_type>(
                                       nRow, comm, params);
    retval = dynamic_cast<basedist_t*>(dist);
  }
  else if (distribution == "MMFile") {
    // Generate user-defined distribution from Matrix-Market file
    if (randomize && me == 0) {
      std::cout << "WARNING: Randomization not available for MMFile "
                << "distributions." << std::endl;
    }
    DistributionMMFile<global_ordinal_type,scalar_type> *dist =
            new DistributionMMFile<global_ordinal_type,scalar_type>(
                                   nRow, comm, params);
    retval = dynamic_cast<basedist_t*>(dist);
  }

  else {
    std::cout << "ERROR:  Invalid distribution method " << distribution 
              << std::endl;
    exit(-1);
  }
  return Teuchos::rcp<basedist_t>(retval);
}

static
void 
readMatrixMarket(  
  const std::string &filename,    // MatrixMarket file to read
  const Teuchos::RCP<const Teuchos::Comm<int> > &comm,  
  const Teuchos::ParameterList &params,
  size_t &nRow,
  size_t &nCol,
  typename Distribution<global_ordinal_type,scalar_type>::LocalNZmap_t &localNZ,
  Teuchos::RCP<Distribution<global_ordinal_type,scalar_type> > &dist
)
{
  bool useTimers = false;   // should we time various parts of the reader?
  {
  const Teuchos::ParameterEntry *pe = params.getEntryPtr("useTimers");
  if (pe != NULL) 
    useTimers = pe->getValue<bool>(&useTimers);
  }

  Teuchos::RCP<Teuchos::TimeMonitor> timer = Teuchos::null;
  if (useTimers) 
    timer = rcp(new Teuchos::TimeMonitor(
                    *Teuchos::TimeMonitor::getNewTimer("RMM parameterSetup")));

  int me = comm->getRank();

  bool verbose = false;   // Print status as reading
  {
  const Teuchos::ParameterEntry *pe = params.getEntryPtr("verbose");
  if (pe != NULL) 
    verbose = pe->getValue<bool>(&verbose);
  }

  size_t chunkSize = 500000;   // Number of lines to read / broadcast at once
  {
  const Teuchos::ParameterEntry *pe = params.getEntryPtr("chunkSize");
  if (pe != NULL) 
    chunkSize = pe->getValue<size_t>(&chunkSize);
  }

  bool symmetrize = false;  // Symmetrize the matrix
  {
  const Teuchos::ParameterEntry *pe = params.getEntryPtr("symmetrize");
  if (pe != NULL) 
    symmetrize = pe->getValue<bool>(&symmetrize);
  }

  bool transpose = false;   // Store the transpose
  {
  const Teuchos::ParameterEntry *pe = params.getEntryPtr("transpose");
  if (pe != NULL) 
    transpose = pe->getValue<bool>(&transpose);
  }

  std::string diagonal = "";  // Are diagonal entries required (so add them)
                              // or ignored (so delete them) or neither?
                              // Default is neither.
  {
  const Teuchos::ParameterEntry *pe = params.getEntryPtr("diagonal");
  if (pe != NULL) 
    diagonal = pe->getValue<std::string>(&diagonal);
  }
  bool ignoreDiagonal = (diagonal == "exclude");
  bool requireDiagonal = (diagonal == "require");

  std::string distribution = "1D";  // Default distribution is 1D row-based
  {
  const Teuchos::ParameterEntry *pe = params.getEntryPtr("distribution");
  if (pe != NULL) 
    distribution = pe->getValue<std::string>(&distribution);
  }

  if (useTimers) {
    timer = Teuchos::null;
    timer = rcp(new Teuchos::TimeMonitor(
                    *Teuchos::TimeMonitor::getNewTimer("RMM header")));
  }

  if (verbose && me == 0)
    std::cout << "Reading matrix market file... " << filename << std::endl;

  FILE *fp = NULL;
  size_t dim[3] = {0,0,0};  // nRow, nCol, nNz as read from MatrixMarket
  MM_typecode mmcode;

  if (me == 0) {

    fp = fopen(filename.c_str(), "r");

    if (fp == NULL) {
      std::cout << "Error:  cannot open file " << filename << std::endl;
    }
    else {
      // Read MatrixMarket banner to get type of matrix
      if (mm_read_banner(fp, &mmcode) != 0) {
        std::cout << "Error:  bad MatrixMarket banner " << std::endl;
      }
      else {
        // Error check for file types that this function can read
        if (!mm_is_valid(mmcode) ||
             mm_is_dense(mmcode) || mm_is_array(mmcode) ||
             mm_is_complex(mmcode) || 
             mm_is_skew(mmcode) || mm_is_hermitian(mmcode)) {
          std::cout << "Error:  matrix type is not supported by this reader\n "
                    << "This reader supports sparse, coordinate, "
                    << "real, pattern, integer, symmetric, general"
                    << std::endl;
        }
        else {
          // Read nRow, nCol, nNz from MatrixMarket file
          if (mm_read_mtx_crd_size(fp, &dim[0], &dim[1], &dim[2]) != 0) {
            std::cout << "Error:  bad matrix dimensions " << std::endl;
            dim[0] = dim[1] = dim[2] = 0;
          }
        }
      }
    }
  }

  // Broadcast relevant info
  // Bad input if dim[0] or dim[1] still is zero after broadcast; 
  // all procs throw together
  Teuchos::broadcast<int, size_t>(*comm, 0, 3, dim);
  if (dim[0] == 0 || dim[1] == 0) {
    throw std::runtime_error("Error:  bad matrix header information");
  }
  Teuchos::broadcast<int, char>(*comm, 0, sizeof(MM_typecode), mmcode);

  nRow = dim[0];
  nCol = dim[1];

  TEUCHOS_TEST_FOR_EXCEPTION(nRow != nCol, std::logic_error,
        "This overload of readSparseFile requires nRow == nCol " 
         << "(nRow = " << nRow << ", nCol = " << nCol << "); "
         << "For now, use a different overload of readSparseFile until this "
         << "overload is fixed to read rectangular matrices. "
         << "See Trilinos github issues #7045 and #8472.");
  
  size_t nNz = dim[2];
  bool patternInput = mm_is_pattern(mmcode);
  bool symmetricInput = mm_is_symmetric(mmcode);
  if (symmetricInput) symmetrize = true;

  if (verbose && me == 0)
    std::cout << "Matrix market file... " 
              << "\n  symmetrize          = " << symmetrize
              << "\n  transpose           = " << transpose
              << "\n  change diagonal     = " << diagonal
              << "\n  distribution        = " << distribution
              << std::endl;

  if (useTimers) {
    timer = Teuchos::null;
    timer = rcp(new Teuchos::TimeMonitor(
                    *Teuchos::TimeMonitor::getNewTimer("RMM distribution")));
  }
  
  // Create distribution based on nRow, nCol, npRow, npCol
  dist = buildDistribution<global_ordinal_type,scalar_type>(distribution,
							    nRow, nCol, params,
							    comm);
  if (useTimers) {
    timer = Teuchos::null;
    timer = rcp(new Teuchos::TimeMonitor(
                    *Teuchos::TimeMonitor::getNewTimer("RMM readChunks")));
  }

  std::set<global_ordinal_type> diagset;  
                            // If diagonal == require, this set keeps track of 
                            // which diagonal entries already existing so we can
                            // add those that don't

  using nzindex_t = 
        typename Distribution<global_ordinal_type,scalar_type>::NZindex_t;

  // Chunk information and buffers
  const int maxLineLength = 81;
  char *buffer = new char[chunkSize*maxLineLength+1];
  size_t nChunk; 
  size_t nMillion = 0;
  size_t nRead = 0;
  size_t rlen;

  auto timerRead = Teuchos::TimeMonitor::getNewTimer("RMM    readNonzeros");
  auto timerSelect = Teuchos::TimeMonitor::getNewTimer("RMM    selectNonzeros");
  // Read chunks until the entire file is read
  Teuchos::RCP<Teuchos::TimeMonitor> innerTimer = Teuchos::null;
  while (nRead < nNz) {

    innerTimer = rcp(new Teuchos::TimeMonitor(*timerRead));

    if (nNz-nRead > chunkSize) nChunk = chunkSize;
    else nChunk = (nNz - nRead);

    // Processor 0 reads a chunk
    if (me == 0) {
      char *eof;
      rlen = 0;
      for (size_t i = 0; i < nChunk; i++) {
        eof = fgets(&buffer[rlen],maxLineLength,fp);
        if (eof == NULL) {
          std::cout << "Unexpected end of matrix file." << std::endl;
          std::cout.flush();
          exit(-1);
        }
        rlen += strlen(&buffer[rlen]);
      }
      buffer[rlen++]='\n';
    }

    // Processor 0 broadcasts the chunk
    Teuchos::broadcast<int, size_t>(*comm, 0, 1, &rlen);
    Teuchos::broadcast<int, char>(*comm, 0, rlen, buffer);

    buffer[rlen++] = '\0';
    nRead += nChunk;

    innerTimer = Teuchos::null;
    innerTimer = rcp(new Teuchos::TimeMonitor(*timerSelect));

    // All processors check the received data, saving nonzeros belonging to them
    char *lineptr = buffer;
    for (rlen = 0; rlen < nChunk; rlen++) {

      char *next = strchr(lineptr,'\n');
      global_ordinal_type I = atol(strtok(lineptr," \t\n")) 
              - 1; // since MatrixMarket files are one-based
      global_ordinal_type J = atol(strtok(NULL," \t\n")) 
              - 1; // since MatrixMarket files are one-based

      scalar_type V = (patternInput ? -1. : atof(strtok(NULL," \t\n")));
      lineptr = next + 1;

      // Special processing of nonzero
      if ((I == J) && ignoreDiagonal) continue;

      if (transpose) std::swap<global_ordinal_type>(I,J);

      // Add nonzero (I,J) to the map if it should be on this processor
      // Some file-based distributions have processor assignment stored as 
      // the non-zero's value, so pass the value to Mine.
      if (dist->Mine(I,J,int(V))) {
        nzindex_t idx = std::make_pair(I,J);
        localNZ[idx] = V;   
        if (requireDiagonal && (I == J)) diagset.insert(I);
      }

      // If symmetrizing, add (J,I) to the map if it should be on this processor
      // Some file-based distributions have processor assignment stored as 
      // the non-zero's value, so pass the value to Mine.
      if (symmetrize && (I != J) && dist->Mine(J,I,int(V))) {
        //  Add entry (J, I) if need to symmetrize
        //  This processor keeps this non-zero.
        nzindex_t idx = std::make_pair(J,I);
        localNZ[idx] = V;   
      }
    }

    // Status check
    if (verbose) {
      if (nRead / 1000000 > nMillion) {
        nMillion++;
        if (me == 0) std::cout << nMillion << "M ";
      }
    }

    innerTimer = Teuchos::null;
  }

  if (verbose && me == 0) 
    std::cout << std::endl << nRead << " nonzeros read " << std::endl;

  if (fp != NULL) fclose(fp);
  delete [] buffer;

  if (useTimers) {
    timer = Teuchos::null;
    timer = rcp(new Teuchos::TimeMonitor(
                    *Teuchos::TimeMonitor::getNewTimer("RMM diagonal")));
  }

  //  Add diagonal entries if they are required.
  //  Check to make sure they are all here; add them if they are not.
  if (requireDiagonal) {
    if (dist->DistType() == MMFile) {
      // All diagonal entries should be present in the file; we cannot create
      // them for file-based data distributions.  
      // Do an error check to ensure they are all there.
      size_t localss = diagset.size();
      size_t globalss;
      Teuchos::reduceAll<int, size_t>(*comm, Teuchos::REDUCE_SUM, 1,
                                      &localss, &globalss);
      TEUCHOS_TEST_FOR_EXCEPTION(globalss != nRow, std::logic_error,
        "File-based nonzero distribution requires all diagonal "
         << "entries to be present in the file.  \n" 
         << "Number of diagonal entries in file = " << globalss << "\n"
         << "Number of matrix rows = " << nRow << "\n");
    }
    else {
      for (global_ordinal_type i = 0;
                               i < static_cast<global_ordinal_type>(nRow); i++) 
      {
        if (dist->Mine(i,i)) {
          if (diagset.find(i) == diagset.end()) {
            nzindex_t idx = std::make_pair(i,i);
            localNZ[idx] = 0;   
          }
        }
      }
    }
  }
  // Done with diagset; free its memory
  std::set<global_ordinal_type>().swap(diagset);

  if (useTimers) 
    timer = Teuchos::null;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////  BINARY FILE READER & RESPECTIVE FORMAT ///////////////////
//        (This format is used by the upcoming readBinary function)         //
//////////////////////////////////////////////////////////////////////////////
// 
// File format:
//   #vertices #edges src_1 dst_1 src_2 dst_2 ... src_#edges dst_#edges
// 
// Types:
//   #edges:            unsigned long long
//   everything else:   unsigned int
//
// Base of indexing: 
//   1
// 
//////////////////////////////////////////////////////////////////////////////
//
// A code example that converts MatrixMarket format into this format:
//
//   typedef unsigned int ord_type;
//   typedef unsigned long long size_type;
//
//   std::string line;
//   std::ifstream in(infilename);
//   std::ofstream out(outfilename,  std::ios::out | std::ios::binary);
//
//  ord_type nv, dummy, edge[2];
//  size_type ne;
//
//   do
//     std::getline (in, line);
//   while(line[0] == '%');
//
//   std::stringstream ss(line);
//   ss >> nv >> dummy >> ne;
//   out.write((char *)&nv, sizeof(ord_type));
//   out.write((char *)&ne, sizeof(size_type));
//
//   while(std::getline(in, line)) {
//     std::stringstream ss2(line);
//     ss2 >> edge[0] >> edge[1];
//     out.write((char *)edge, sizeof(ord_type)*2);
//
//    }
//   in.close();
//   out.close();
//
//////////////////////////////////////////////////////////////////////////////
  
static
void 
readBinary(  
  const std::string &filename,    // MatrixMarket file to read
  const Teuchos::RCP<const Teuchos::Comm<int> > &comm,  
  const Teuchos::ParameterList &params,
  size_t &nRow,
  size_t &nCol,
  typename Distribution<global_ordinal_type,scalar_type>::LocalNZmap_t &localNZ,
  Teuchos::RCP<Distribution<global_ordinal_type,scalar_type> > &dist
)
{

  int me = comm->getRank();

  bool verbose = false;   // Print status as reading
  {
  const Teuchos::ParameterEntry *pe = params.getEntryPtr("verbose");
  if (pe != NULL) 
    verbose = pe->getValue<bool>(&verbose);
  }

  size_t chunkSize = 500000;   // Number of lines to read / broadcast at once
  {
  const Teuchos::ParameterEntry *pe = params.getEntryPtr("chunkSize");
  if (pe != NULL) 
    chunkSize = pe->getValue<size_t>(&chunkSize);
  }

  bool symmetrize = false;  // Symmetrize the matrix
  {
  const Teuchos::ParameterEntry *pe = params.getEntryPtr("symmetrize");
  if (pe != NULL) 
    symmetrize = pe->getValue<bool>(&symmetrize);
  }

  bool transpose = false;   // Store the transpose
  {
  const Teuchos::ParameterEntry *pe = params.getEntryPtr("transpose");
  if (pe != NULL) 
    transpose = pe->getValue<bool>(&transpose);
  }

  std::string diagonal = "";  // Are diagonal entries required (so add them)
                              // or ignored (so delete them) or neither?
                              // Default is neither.
  {
  const Teuchos::ParameterEntry *pe = params.getEntryPtr("diagonal");
  if (pe != NULL) 
    diagonal = pe->getValue<std::string>(&diagonal);
  }
  bool ignoreDiagonal = (diagonal == "exclude");
  bool requireDiagonal = (diagonal == "require");

  std::string distribution = "1D";  // Default distribution is 1D row-based
  {
  const Teuchos::ParameterEntry *pe = params.getEntryPtr("distribution");
  if (pe != NULL) 
    distribution = pe->getValue<std::string>(&distribution);
  }

  if (verbose && me == 0)
    std::cout << "Reading binary file... " << filename << std::endl;

  FILE *fp = NULL;
  size_t dim[3] = {0,0,0};  // Expected to read nRow and nNz, nCol = nRow

  if (me == 0) {

    fp = fopen(filename.c_str(), "rb");

    if (fp == NULL) {
      std::cout << "Error:  cannot open file " << filename << std::endl;
    }
    else {
      // The header in the binary file contains only numVertices and numEdges
      unsigned int nv = 0;
      unsigned long long ne = 0;
      if (fread(&nv, sizeof(unsigned int), 1, fp) != 1)
           throw std::runtime_error("Error: bad number of read values.");
      if (fread(&ne, sizeof(unsigned long long), 1, fp) != 1)
           throw std::runtime_error("Error: bad number of read values.");
      dim[0] = nv;
      dim[1] = dim[0];
      dim[2] = ne;
    }
  }

  // Broadcast relevant info
  // Bad input if dim[0] or dim[1] still is zero after broadcast; 
  // all procs throw together
  Teuchos::broadcast<int, size_t>(*comm, 0, 3, dim);
  if (dim[0] == 0 || dim[1] == 0) {
    throw std::runtime_error("Error:  bad matrix header information");
  }

  nRow = dim[0];
  nCol = dim[1];
  size_t nNz = dim[2];

  if (verbose && me == 0)
    std::cout << "Binary file... " 
              << "\n  symmetrize          = " << symmetrize
              << "\n  transpose           = " << transpose
              << "\n  change diagonal     = " << diagonal
              << "\n  distribution        = " << distribution
              << std::endl;
  
  // Create distribution based on nRow, nCol, npRow, npCol
  dist = buildDistribution<global_ordinal_type,scalar_type>(distribution,
							    nRow, nCol, params,
							    comm);

  std::set<global_ordinal_type> diagset;  
                            // If diagonal == require, this set keeps track of 
                            // which diagonal entries already existing so we can
                            // add those that don't

  using nzindex_t = 
        typename Distribution<global_ordinal_type,scalar_type>::NZindex_t;

  // Chunk information and buffers
  unsigned int *buffer = new unsigned int[chunkSize*2];
  size_t nChunk; 
  size_t nMillion = 0;
  size_t nRead = 0;
  size_t rlen;
  const scalar_type ONE = Teuchos::ScalarTraits<scalar_type>::one();


  // Read chunks until the entire file is read
  while (nRead < nNz) {
    if (nNz-nRead > chunkSize) nChunk = chunkSize;
    else nChunk = (nNz - nRead);

    // Processor 0 reads a chunk
    if (me == 0) {
      size_t ret = 0;
      rlen = 0;
      for (size_t i = 0; i < nChunk; i++) {
	ret = fread(&buffer[rlen], sizeof(unsigned int), 2, fp);
        if (ret == 0) {
          std::cout << "Unexpected end of matrix file." << std::endl;
          std::cout.flush();
          exit(-1);
        }
        rlen += 2;
      }
    }

    // Processor 0 broadcasts the chunk
    Teuchos::broadcast<int, size_t>(*comm, 0, 1, &rlen);
    Teuchos::broadcast<int, unsigned int>(*comm, 0, rlen, buffer);

    nRead += nChunk;

    // All processors check the received data, saving nonzeros belonging to them
    for (rlen = 0; rlen < nChunk; rlen++) {

      global_ordinal_type I = buffer[2*rlen]-1;
      global_ordinal_type J = buffer[2*rlen+1]-1;
      
      // Special processing of nonzero
      if ((I == J) && ignoreDiagonal) continue;

      if (transpose) std::swap<global_ordinal_type>(I,J);

      // Add nonzero (I,J) to the map if it should be on this processor
      // Some file-based distributions have processor assignment stored as 
      // the non-zero's value, so pass the value to Mine.
      if (dist->Mine(I,J,ONE)) {
        nzindex_t idx = std::make_pair(I,J);
        localNZ[idx] = ONE;  // For now, the input binary format does not 
	                     // support numeric values, so we insert one.  
        if (requireDiagonal && (I == J)) diagset.insert(I);
      }

      // If symmetrizing, add (J,I) to the map if it should be on this processor
      // Some file-based distributions have processor assignment stored as 
      // the non-zero's value, so pass the value to Mine.
      if (symmetrize && (I != J) && dist->Mine(J,I,ONE)) {
        //  Add entry (J, I) if need to symmetrize
        //  This processor keeps this non-zero.
        nzindex_t idx = std::make_pair(J,I);
        localNZ[idx] = ONE;   
      }
    }

    // Status check
    if (verbose) {
      if (nRead / 1000000 > nMillion) {
        nMillion++;
        if (me == 0) std::cout << nMillion << "M ";
      }
    }
  }

  if (verbose && me == 0) 
    std::cout << std::endl << nRead << " nonzeros read " << std::endl;

  if (fp != NULL) fclose(fp);
  delete [] buffer;

  //  Add diagonal entries if they are required.
  //  Check to make sure they are all here; add them if they are not.
  if (requireDiagonal) {
    if (dist->DistType() == MMFile) {
      // All diagonal entries should be present in the file; we cannot create
      // them for file-based data distributions.  
      // Do an error check to ensure they are all there.
      size_t localss = diagset.size();
      size_t globalss;
      Teuchos::reduceAll<int, size_t>(*comm, Teuchos::REDUCE_SUM, 1,
                                      &localss, &globalss);
      TEUCHOS_TEST_FOR_EXCEPTION(globalss != nRow, std::logic_error,
        "File-based nonzero distribution requires all diagonal "
         << "entries to be present in the file.  \n" 
         << "Number of diagonal entries in file = " << globalss << "\n"
         << "Number of matrix rows = " << nRow << "\n");
    }
    else {
      for (global_ordinal_type i = 0;
                               i < static_cast<global_ordinal_type>(nRow); i++)
      {
        if (dist->Mine(i,i)) {
          if (diagset.find(i) == diagset.end()) {
            nzindex_t idx = std::make_pair(i,i);
            localNZ[idx] = 0;   
          }
        }
      }
    }
  }
  // Done with diagset; free its memory
  std::set<global_ordinal_type>().swap(diagset);

}

//////////////////////////////////////////////////////////////////////////////
//////////// PER-PROCESS BINARY FILE READER & RESPECTIVE FORMAT //////////////
//   (This format is used by the upcoming readPerProcessBinary function)    //
//////////////////////////////////////////////////////////////////////////////
// 
// 
// FILE FORMAT:
//     globalNumRows globalNumCols localNumNnzs row_1 col_1 ... row_n col_n, 
//  where n = localNumNnzs
// 
//
// ASSUMPTIONS:
//  - The nonzeros should be sorted by their row indices within each file.
//  - Nonzeros have global row and column indices.
//  - The user specifies the base <filename>, where rank i reads <filename>.i.cooBin.  
// 
// 
// TYPES:
//  localNumNnzs:      unsigned long long
//  everything else:   unsigned int
//
//
// BASE OF INDEXING: 1
//
// 
//////////////////////////////////////////////////////////////////////////////
static
void 
readPerProcessBinary(
  const std::string &filename,
  const Teuchos::RCP<const Teuchos::Comm<int> > &comm,  
  const Teuchos::ParameterList &params,
  size_t &nRow,
  size_t &nCol,
  typename Distribution<global_ordinal_type,scalar_type>::LocalNZmap_t &localNZ,
  Teuchos::RCP<Distribution<global_ordinal_type,scalar_type> > &dist,
  unsigned int* &buffer,
  size_t &nNz
)
{

  int me = comm->getRank();

  bool verbose = false;   // Print status as reading
  {
  const Teuchos::ParameterEntry *pe = params.getEntryPtr("verbose");
  if (pe != NULL) 
    verbose = pe->getValue<bool>(&verbose);
  }

  std::string distribution = "1D";  // Default distribution is 1D row-based
  {
  const Teuchos::ParameterEntry *pe = params.getEntryPtr("distribution");
  if (pe != NULL) 
    distribution = pe->getValue<std::string>(&distribution);
  }

  if (verbose && me == 0)
    std::cout << "Reading per-process binary files... " << filename << std::endl;


  std::string rankFileName = filename + "." + std::to_string(me) + ".cooBin";
  FILE *fp = NULL;

  fp = fopen(rankFileName.c_str(), "rb");
  if (fp == NULL) {
    std::cout << "Error:  cannot open file " << filename << std::endl;
    throw std::runtime_error("Error:  non-existing input file: " + rankFileName);
  }
 
  // The header in each per-process file:  globalNumRows globalNumCols localNumNonzeros
  unsigned int globalNumRows = 0, globalNumCols = 0;
  unsigned long long localNumNonzeros = 0;
  if (fread(&globalNumRows, sizeof(unsigned int), 1, fp) != 1)
    throw std::runtime_error("Error: bad number of read values.");
  if (fread(&globalNumCols, sizeof(unsigned int), 1, fp) != 1)
    throw std::runtime_error("Error: bad number of read values.");
  if (fread(&localNumNonzeros, sizeof(unsigned long long), 1, fp) != 1)
    throw std::runtime_error("Error: bad number of read values.");
 
  nRow = static_cast<size_t>(globalNumRows);
  nCol = static_cast<size_t>(globalNumCols);
  nNz = static_cast<size_t>(localNumNonzeros);

  // Fill the simple buffer array instead of a std::map
  // S. Acer: With large graphs, we can't afford std::map
  buffer = new unsigned int[nNz*2];

  if(nNz > 0) {
    size_t ret = fread(buffer, sizeof(unsigned int), 2*nNz, fp);
    if (ret == 0) {
      std::cout << "Unexpected end of matrix file: " << rankFileName << std::endl;
      std::cout.flush();
      delete [] buffer;
      exit(-1);
    }
  }
  if (fp != NULL) fclose(fp);

  // This barrier is not necessary but useful for debugging
  comm->barrier();
  if(verbose && me == 0)
    std::cout << "All ranks finished reading their nonzeros from their individual files\n";

  // Create distribution based on nRow, nCol, npRow, npCol
  dist = buildDistribution<global_ordinal_type,scalar_type>(distribution,
							    nRow, nCol, params,
							    comm);
}

public:

// This is the default interface.
static Teuchos::RCP<sparse_matrix_type>
readSparseFile(
  const std::string &filename,    // MatrixMarket file to read
  const Teuchos::RCP<const Teuchos::Comm<int> > &comm,  
  const Teuchos::ParameterList &params
)
{
  Teuchos::RCP<Distribution<global_ordinal_type,scalar_type> > dist;
  return readSparseFile(filename, comm, params, dist);
}

// This version has the Distribution object as an output parameter.
// S. Acer needs the distribution object to get the chunk cuts from
// LowerTriangularBlock distribution.  
static Teuchos::RCP<sparse_matrix_type>
readSparseFile(
  const std::string &filename,    // MatrixMarket file to read
  const Teuchos::RCP<const Teuchos::Comm<int> > &comm,  
  const Teuchos::ParameterList &params,
  Teuchos::RCP<Distribution<global_ordinal_type,scalar_type> > &dist 
)
{
  bool useTimers = false;   // should we time various parts of the reader?
  {
  const Teuchos::ParameterEntry *pe = params.getEntryPtr("useTimers");
  if (pe != NULL) 
    useTimers = pe->getValue<bool>(&useTimers);
  }

  Teuchos::RCP<Teuchos::TimeMonitor> timer = Teuchos::null;
  if (useTimers) 
    timer = rcp(new Teuchos::TimeMonitor(
                    *Teuchos::TimeMonitor::getNewTimer("RSF parameterSetup")));

  int me = comm->getRank();

  // Check parameters to determine how to process the matrix while reading
  // TODO:  Add validators for the parameters

  bool verbose = false;   // Print status as reading
  {
  const Teuchos::ParameterEntry *pe = params.getEntryPtr("verbose");
  if (pe != NULL) 
    verbose = pe->getValue<bool>(&verbose);
  }
  
  bool callFillComplete = true;   // should we fillComplete the new CrsMatrix?
  {
  const Teuchos::ParameterEntry *pe = params.getEntryPtr("callFillComplete");
  if (pe != NULL) 
    callFillComplete = pe->getValue<bool>(&callFillComplete);
  }

  // Don't want to use MatrixMarket's coordinate reader, because don't want
  // entire matrix on one processor.
  // Instead, Proc 0 reads nonzeros in chunks and broadcasts chunks to all 
  // processors.
  // All processors insert nonzeros they own into a std::map

  // Storage for this processor's nonzeros.
  using localNZmap_t = 
        typename Distribution<global_ordinal_type,scalar_type>::LocalNZmap_t;
  localNZmap_t localNZ;

  bool binary = false;   // should we read a binary file?
  {
  const Teuchos::ParameterEntry *pe = params.getEntryPtr("binary");
  if (pe != NULL) 
    binary = pe->getValue<bool>(&binary);
  }

  bool readPerProcess = false;   // should we read a separate file per process?
  {
  const Teuchos::ParameterEntry *pe = params.getEntryPtr("readPerProcess");
  if (pe != NULL) 
    readPerProcess = pe->getValue<bool>(&readPerProcess);
  }

  if (useTimers) {
    const char *timername = (binary?"RSF readBinary":"RSF readMatrixMarket");
    timer = Teuchos::null;
    timer = rcp(new Teuchos::TimeMonitor(
                   *Teuchos::TimeMonitor::getNewTimer(timername)));
  }
  
  // Read nonzeros from the given file(s)
  size_t nRow = 0, nCol = 0;
  unsigned int *buffer=0; size_t nNz = 0;
  if(binary){
    if(readPerProcess)
      readPerProcessBinary(filename, comm, params, nRow, nCol, localNZ, dist, buffer, nNz);
    else 
      readBinary(filename, comm, params, nRow, nCol, localNZ, dist);
  }
  else
    readMatrixMarket(filename, comm, params, nRow, nCol, localNZ, dist);

  if(readPerProcess == false){

    // Redistribute nonzeros as needed to satisfy the Distribution
    // For most Distributions, this is a no-op    
    if (useTimers) {
      timer = Teuchos::null;
      timer = rcp(new Teuchos::TimeMonitor(
                   *Teuchos::TimeMonitor::getNewTimer("RSF redistribute")));
    }
  
    dist->Redistribute(localNZ);
  }

  if (useTimers) {
    timer = Teuchos::null;
    timer = rcp(new Teuchos::TimeMonitor(
                   *Teuchos::TimeMonitor::getNewTimer("RSF nonzerosConstruction")));
  }

  //  Now construct the matrix.
  //  Count number of entries in each row for best use of StaticProfile
  if (verbose && me == 0)
    std::cout << std::endl << "Constructing the matrix" << std::endl;

  Teuchos::Array<global_ordinal_type> rowIdx;
  Teuchos::Array<size_t> nnzPerRow;
  Teuchos::Array<global_ordinal_type> colIdx;
  Teuchos::Array<scalar_type> val;
  Teuchos::Array<size_t> offsets;

  if(readPerProcess) {
    global_ordinal_type prevI = std::numeric_limits<global_ordinal_type>::max();
    for (size_t it = 0; it < nNz; it++){
      global_ordinal_type I = buffer[2*it]-1;
      global_ordinal_type J = buffer[2*it+1]-1;
      if (prevI != I) {
	prevI = I;
	rowIdx.push_back(I);
	nnzPerRow.push_back(0);
      }
      nnzPerRow.back()++;
      colIdx.push_back(J);
    }
    delete [] buffer;      
  }
  else {
    // Exploit fact that map has entries sorted by I, then J
    global_ordinal_type prevI = std::numeric_limits<global_ordinal_type>::max();
    for (auto it = localNZ.begin(); it != localNZ.end(); it++) {
      global_ordinal_type I = it->first.first;
      global_ordinal_type J = it->first.second;
      scalar_type V = it->second;
      if (prevI != I) {
	prevI = I;
	rowIdx.push_back(I);
	nnzPerRow.push_back(0);
      }
      nnzPerRow.back()++;
      colIdx.push_back(J);
      val.push_back(V);
    }

    // Done with localNZ map; free it.
    localNZmap_t().swap(localNZ);
  }

  // Compute prefix sum in offsets array
  offsets.resize(rowIdx.size() + 1);
  offsets[0] = 0;
  for (size_t row = 0; row < static_cast<size_t>(rowIdx.size()); row++)
    offsets[row+1] = offsets[row] + nnzPerRow[row];

  if (useTimers) {
    timer = Teuchos::null;
    timer = rcp(new Teuchos::TimeMonitor(
                   *Teuchos::TimeMonitor::getNewTimer("RSF insertNonzeros")));
  }

  // Create a new RowMap with only rows having non-zero entries.
  size_t dummy = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  Teuchos::RCP<const Tpetra::Map<> > rowMap = 
       Teuchos::rcp(new Tpetra::Map<>(dummy, rowIdx(), 0, comm));

  Teuchos::RCP<sparse_matrix_type> A = 
           rcp(new sparse_matrix_type(rowMap, nnzPerRow()));

  //  Insert the global values into the matrix row-by-row.
  if (verbose && me == 0) 
    std::cout << "Inserting global values" << std::endl;

  if(readPerProcess){
    const scalar_type ONE = Teuchos::ScalarTraits<scalar_type>::one();
    for (int i = 0; i < rowIdx.size(); i++) {
      size_t nnz = nnzPerRow[i];
      size_t off = offsets[i];
      val.assign(nnz, ONE);
      // ReadPerProcess routine does not read any numeric values from the file,
      // So we insert dummy values here. 
      A->insertGlobalValues(rowIdx[i], colIdx(off, nnz), val());
    }
  }
  else {
    for (int i = 0; i < rowIdx.size(); i++) {
      size_t nnz = nnzPerRow[i];
      size_t off = offsets[i];
      A->insertGlobalValues(rowIdx[i], colIdx(off, nnz), val(off, nnz));
    }
  }

  // free memory that is no longer needed
  Teuchos::Array<size_t>().swap(nnzPerRow);
  Teuchos::Array<size_t>().swap(offsets);
  Teuchos::Array<global_ordinal_type>().swap(rowIdx);
  Teuchos::Array<global_ordinal_type>().swap(colIdx);
  Teuchos::Array<scalar_type>().swap(val);

  if (useTimers)
    timer = Teuchos::null;

  if (callFillComplete) {

    if (verbose && me == 0) 
      std::cout << "Building vector maps" << std::endl;

    if (useTimers) {
      timer = Teuchos::null;
      timer = rcp(new Teuchos::TimeMonitor(
                     *Teuchos::TimeMonitor::getNewTimer("RSF buildVectorMaps")));
    }

    // Build domain map that corresponds to distribution 
    Teuchos::Array<global_ordinal_type> vectorSet;
    for (global_ordinal_type i = 0;
                             i < static_cast<global_ordinal_type>(nCol); i++) 
      if (dist->VecMine(i)) vectorSet.push_back(i);

    dummy = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    Teuchos::RCP<const Tpetra::Map<> > domainMap = 
           Teuchos::rcp(new Tpetra::Map<>(dummy, vectorSet(), 0, comm));
  
    Teuchos::Array<global_ordinal_type>().swap(vectorSet);

    // Build range map that corresponds to distribution
    for (global_ordinal_type i = 0; 
                             i < static_cast<global_ordinal_type>(nRow); i++) 
      if (dist->VecMine(i)) vectorSet.push_back(i);
  
    dummy = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    Teuchos::RCP<const Tpetra::Map<> > rangeMap = 
           Teuchos::rcp(new Tpetra::Map<>(dummy, vectorSet(), 0, comm));

    Teuchos::Array<global_ordinal_type>().swap(vectorSet);

    // FillComplete the matrix
    if (useTimers) {
      timer = Teuchos::null;
      timer = rcp(new Teuchos::TimeMonitor(
                     *Teuchos::TimeMonitor::getNewTimer("RSF fillComplete")));
    }

    if (verbose && me == 0) 
      std::cout << "Calling fillComplete" << std::endl;

    A->fillComplete(domainMap, rangeMap);

    if (useTimers) 
      timer = Teuchos::null;

    if (verbose) {
      std::cout << "\nRank " << A->getComm()->getRank() 
                << "  nRows " << A->getLocalNumRows()
                << "  minRow " << A->getRowMap()->getMinGlobalIndex()
                << "  maxRow " << A->getRowMap()->getMaxGlobalIndex()
                << "\nRank " << A->getComm()->getRank() 
                << "  nCols " << A->getLocalNumCols()
                << "  minCol " << A->getColMap()->getMinGlobalIndex()
                << "  maxCol " << A->getColMap()->getMaxGlobalIndex()
                << "\nRank " << A->getComm()->getRank() 
                << "  nDomain " << A->getDomainMap()->getLocalNumElements()
                << "  minDomain " << A->getDomainMap()->getMinGlobalIndex()
                << "  maxDomain " << A->getDomainMap()->getMaxGlobalIndex()
                << "\nRank " << A->getComm()->getRank() 
                << "  nRange " << A->getRangeMap()->getLocalNumElements()
                << "  minRange " << A->getRangeMap()->getMinGlobalIndex()
                << "  maxRange " << A->getRangeMap()->getMaxGlobalIndex()
                << "\nRank " << A->getComm()->getRank() 
                << "  nEntries " << A->getLocalNumEntries()
                << std::endl;
    }
  }

  return A;
}



#endif

