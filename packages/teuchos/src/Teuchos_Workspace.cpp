// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
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
  TEST_FOR_EXCEPTION(
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
      TEST_FOR_EXCEPTION(
        workspace_store_->curr_ws_ptr_ != workspace_end_, std::logic_error
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
