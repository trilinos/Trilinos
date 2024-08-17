// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// ************************************************************************
// 
// C++ wrappers for Zoltan's Distributed Directory utility
//
// Two styles of initialization: 
//
//   C++ style:  Zoltan_DD dd(comm, num_gid, num_lid, len1, len2, debug);
//
//   C style: Zoltan_DD dd;
//            dd.Create(comm, num_gid, num_lid, len1, len2, debug);
//
// ************************************************************************

#ifndef ZOLTAN_DD_CPP_H_
#define ZOLTAN_DD_CPP_H_

#include "zoltan_dd.h"

class Zoltan_DD {

public:

  Zoltan_DD(const MPI_Comm &comm, const int &num_gid, const int &num_lid, 
    const int &user_length,  const int &table_length, const int &debug_level) 
    {
    Zoltan_DD_Create (&this->DD, comm, num_gid, 
                  num_lid, user_length,  table_length, debug_level);
    }

  Zoltan_DD()
    {
    this->DD = NULL;

    // Creator of this object must call Zoltan_DD::Create to finish
    // initialization.
    }

  int Create(const MPI_Comm &comm, const int &num_gid, const int &num_lid, 
    const int &user_length_in_chars,  const int &table_length, const int &debug_level) 
    {
    if (this->DD)
      {
      Zoltan_DD_Destroy(&this->DD);
      this->DD = NULL;
      }

    int rc =  Zoltan_DD_Create (&this->DD, comm, num_gid, 
                  num_lid, user_length_in_chars,  table_length, debug_level);

    return rc;
    }

  ~Zoltan_DD()
    {
    Zoltan_DD_Destroy (&this->DD) ;
    }

  Zoltan_DD (const Zoltan_DD &dd) // Copy constructor
   {
   this->DD = Zoltan_DD_Copy(dd.DD);
   }

  Zoltan_DD & operator= (const Zoltan_DD &dd) // Copy operator
  {
    Zoltan_DD_Copy_To(&this->DD, dd.DD);

    return *this;
  }
 
  int Update (ZOLTAN_ID_PTR gid, ZOLTAN_ID_PTR lid, 
    char *user, int *partition, const int &count) 
    {
    return Zoltan_DD_Update (this->DD, gid, lid, user, partition, count) ;
    }
  
  int Find (ZOLTAN_ID_PTR gid, ZOLTAN_ID_PTR lid, char *data, 
                   int *partition, const int &count, int *owner) const
    {
    return Zoltan_DD_Find (this->DD, gid, lid, data, partition, count, owner);
    }
  
  int Remove (ZOLTAN_ID_PTR gid, const int &count)
    {
    return Zoltan_DD_Remove (this->DD, gid, count);
    }
  
  int Set_Hash_Fn (unsigned int (*hash) (ZOLTAN_ID_PTR, int, unsigned int))
    {
    return Zoltan_DD_Set_Hash_Fn (this->DD, hash);
    } 

  void Stats () const
    {
    Zoltan_DD_Stats (this->DD) ;
    }
  
  int Print () const
    {
    return Zoltan_DD_Print (this->DD) ;
    }

private:

  Zoltan_DD_Directory *DD;
};

#endif
