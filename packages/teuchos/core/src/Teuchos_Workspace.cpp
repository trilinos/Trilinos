// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_Workspace.hpp"

namespace {
Teuchos::RCP<Teuchos::WorkspaceStore>  default_workspace_store(Teuchos::null);
}

// Global functions

void Teuchos::set_default_workspace_store( const Teuchos::RCP<WorkspaceStore> &default_workspace_store_in )
{
  default_workspace_store = default_workspace_store_in;
}

Teuchos::RCP<Teuchos::WorkspaceStore> Teuchos::get_default_workspace_store()
{
  return default_workspace_store;
}

void Teuchos::print_memory_usage_stats( const WorkspaceStore* workspace_store, std::ostream& out )
{
  if( workspace_store ) {
    out
      << "\n*** Statistics for autmatic array workspace:"
      << "\n  Number of megabytes of preallocated workspace                = "
      << (workspace_store->num_bytes_total()*1e-6)
      << "\n  Number of megabytes needed                                   = "
      << (workspace_store->num_max_bytes_needed()*1e-6)
      << "\n  Number of allocations using preallocated workspace           = "
      << workspace_store->num_static_allocations()
      << "\n  Number of dynamic allocations beyond preallocated workspace  = "
      << workspace_store->num_dyn_allocations()
      << "\n";
  }
  else {
    out
      << "\n*** Statistics for autmatic array workspace:"
      << "\n  No workspace storage was allocated!\n";
  }
}

namespace Teuchos {

// WorkspaceStore

WorkspaceStore::WorkspaceStore(size_t num_bytes)
  : workspace_begin_(NULL)
  , workspace_end_(NULL)
  , curr_ws_ptr_(NULL)
  , num_static_allocations_(0)
  , num_dyn_allocations_(0)
  , num_current_bytes_total_(0)
  , num_max_bytes_needed_(0)
{
  if(num_bytes)
    protected_initialize(num_bytes);
}

WorkspaceStore::~WorkspaceStore() {
  if(workspace_begin_) delete [] workspace_begin_;
}

void WorkspaceStore::protected_initialize(size_t num_bytes)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    curr_ws_ptr_ != workspace_begin_, std::logic_error
    ,"WorkspaceStore::set_workspace_size(...) : Error, "
    "You can not reset the workspace size when any RawWorkspace objects "
    "are using workspace!" );
  if(workspace_begin_) delete [] workspace_begin_;
  workspace_begin_        = ::new char[num_bytes];
  workspace_end_          = workspace_begin_ + num_bytes;
  curr_ws_ptr_            = workspace_begin_;
  num_static_allocations_ = 0;
  num_dyn_allocations_    = 0;
  num_current_bytes_total_= 0;
  num_max_bytes_needed_   = 0;
}

// RawWorkspace

RawWorkspace::RawWorkspace(WorkspaceStore* workspace_store, size_t num_bytes_in)
{
  if(num_bytes_in) {
    workspace_store_ = workspace_store;
    if( !workspace_store_ || workspace_store_->num_bytes_remaining() < num_bytes_in ) {
      workspace_begin_ = ::new char[num_bytes_in];
      workspace_end_   = workspace_begin_ + num_bytes_in;
      owns_memory_     = true;
      if(workspace_store_)
        workspace_store_->num_dyn_allocations_++;
    }
    else {
      workspace_begin_ = workspace_store_->curr_ws_ptr_;
      workspace_end_   = workspace_begin_ + num_bytes_in;
      owns_memory_     = false;
      workspace_store_->curr_ws_ptr_ += num_bytes_in;
      workspace_store_->num_static_allocations_++;
    }
  }
  else {
    workspace_store_ = NULL;
    workspace_begin_ = NULL;
    workspace_end_   = NULL;
    owns_memory_     = false;
  }
  if(workspace_store_) {
    workspace_store_->num_current_bytes_total_ += num_bytes_in;
    if( workspace_store_->num_current_bytes_total_ > workspace_store_->num_max_bytes_needed_ )
      workspace_store_->num_max_bytes_needed_ = workspace_store_->num_current_bytes_total_;
  }
}

RawWorkspace::~RawWorkspace()
{
  if(workspace_store_)
    workspace_store_->num_current_bytes_total_ -= this->num_bytes();
  if(owns_memory_) {
    if(workspace_begin_) delete [] workspace_begin_;
  }
  else {
    if(workspace_store_) {
      TEUCHOS_TEST_FOR_TERMINATION(
        workspace_store_->curr_ws_ptr_ != workspace_end_
        ,"RawWorkspace::~RawWorkspace(...): Error, "
        "Invalid usage of RawWorkspace class, corrupted WorspaceStore object!" );
      workspace_store_->curr_ws_ptr_ = workspace_begin_;
    }
  }
}

#ifdef __PGI // Should not have to define this since it should not be called!
void* RawWorkspace::operator new(size_t)
{
  assert(0);
  return NULL;
}
#endif

} // end namespace Teuchos
