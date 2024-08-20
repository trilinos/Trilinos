// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_Array.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"

#include <iostream>
#include <unordered_map>

// Test to determine appropriate setting of Threshold in insert_crs_indices 
// Threshold is number of nonzeros in a row; it determines whether to do
// linear search for duplicates or use an unordered_map .  #8794
// Experiments on 2/22/21 on ghost show cross-over in runtime from serial 
// search to unordered_map occurs at roughly 400 entries being inserted.

// Insert using linear search
size_t insert_crs_indices_linearsearch(
    Teuchos::Array<Tpetra::Map<>::global_ordinal_type> &cur_indices,
    size_t &num_assigned,
    const Teuchos::Array<Tpetra::Map<>::global_ordinal_type> &new_indices)
{
  const size_t start = 0;
  size_t end = start + num_assigned;
  const size_t num_avail = cur_indices.size() - end;

  const size_t num_new_indices = new_indices.size ();
  size_t num_inserted = 0;

  {
    for (size_t k = 0; k < num_new_indices; ++k) {
      const Tpetra::Map<>::global_ordinal_type idx = new_indices[k];
      size_t row_offset = start;
      for (; row_offset < end; ++row_offset) {
        if (idx == cur_indices[row_offset]) {
          break;
        }
      }
      if (row_offset == end) {
        if (num_inserted >= num_avail) { // not enough room
          return Teuchos::OrdinalTraits<size_t>::invalid();
        }
        // This index is not yet in indices
        cur_indices[end++] = idx;
        num_inserted++;
      }
    }
  }
  return num_inserted++;
}

// Insert using unordered_map
size_t insert_crs_indices_table(
    Teuchos::Array<Tpetra::Map<>::global_ordinal_type> &cur_indices,
    size_t &num_assigned,
    const Teuchos::Array<Tpetra::Map<>::global_ordinal_type> &new_indices)
{
  const size_t start = 0;
  size_t end = start + num_assigned;
  const size_t num_avail = cur_indices.size() - end;

  const size_t num_new_indices = new_indices.size ();
  size_t num_inserted = 0;

  {
    std::unordered_map<Tpetra::Map<>::global_ordinal_type, size_t> 
         idxLookup(num_assigned+num_new_indices);

    // Put existing indices into the lookup table
    for (size_t k = 0; k < num_assigned; k++) {
      idxLookup[cur_indices[k]] = start+k;
    }

    // Check for new indices in table; insert if not there yet
    for (size_t k = 0; k < num_new_indices; k++) {
      const Tpetra::Map<>::global_ordinal_type idx = new_indices[k];
      size_t row_offset;

      auto it = idxLookup.find(idx);
      if (it == idxLookup.end()) {
        if (num_inserted >= num_avail) { // not enough room
          return Teuchos::OrdinalTraits<size_t>::invalid();
        }
        // index not found; insert it
        row_offset = end;
        cur_indices[end++] = idx;
        idxLookup[idx] = row_offset;
        num_inserted++;
      }
      else {
        // index found; note its position
        row_offset = it->second;
      }
    }
  }
  num_assigned += num_inserted;
  return num_inserted;
}

////////////////////////////////////////////////////////////////////////
int main(int narg, char **arg)
{
  Tpetra::ScopeGuard scope(&narg, &arg);;

  using gno_t = Tpetra::Map<>::global_ordinal_type;
  int nIter = 10;

  // Test various number of nonzeros per row to insert
  for (int nnz = 1; nnz < 4000; nnz+=100) {
    std::cout << "NNZ = " << nnz << std::endl;

    // Run each test for several iterations
    for (int iter = 0; iter < nIter; iter++) {

      Teuchos::Array<gno_t> curInd_LS(nnz);
      Teuchos::Array<gno_t> newInd_LS(nnz);

      Teuchos::Array<gno_t> curInd_T(nnz);
      Teuchos::Array<gno_t> newInd_T(nnz);

      for (int i = 0; i < nnz; i++) 
        newInd_LS[i] = std::rand() % (nnz * 10000);
      for (int i = 0; i < nnz; i++) 
        newInd_T[i] = newInd_LS[i];

      // Insert entries using linear search
      size_t nInserted_LS; 
      {
        char name[25];
        sprintf(name, "nz=%04d search", nnz);
        Teuchos::RCP<Teuchos::TimeMonitor> timeMonitor(
          new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(name)));
        size_t num_assigned = 0;
        nInserted_LS =
          insert_crs_indices_linearsearch(curInd_LS, num_assigned, newInd_LS);
      }
      
      // Insert entries using unordered_map
      size_t nInserted_T;
      {
        char name[25];
        sprintf(name, "nz=%04d table", nnz);
        Teuchos::RCP<Teuchos::TimeMonitor> timeMonitor(
          new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(name)));
        size_t num_assigned = 0;
        nInserted_T =
          insert_crs_indices_table(curInd_T, num_assigned, newInd_T);
      }

      // Make sure the number of insertions matched for both methods
      if (nInserted_LS != nInserted_T) 
        std::cout << "FAIL nInserted_LS " << nInserted_LS << " != " 
                 << nInserted_T << " nInserted_T" << std::endl;
    }
  }

  // Report the times
  Teuchos::TimeMonitor::summarize();
  std::cout << "PASS" << std::endl;
  return 0;
}
