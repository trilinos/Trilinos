/*
  This is an implementation of the EntityList interface

  Jaideep Ray, SNL, Livermore, 08/23/02
*/

#ifndef EntityListImplHSeen
#define  EntityListImplHSeen

#include "EntityList.h"

namespace ZoltanSpace
{
  class EntityListImpl : public virtual ::LoadPartitionerSpace::EntityList
  {
    public :

    EntityListImpl(int a) : ::LoadPartitionerSpace::EntityList()
    {
      
      proc_size = a ;  gid_size = 0 ;  lid_size = 0 ;  user_def_size = 0;
      gid_list = lid_list = user_def_list = 0 ; proc_list = 0 ;
      ref = 1 ;
    }

    ~EntityListImpl() 
    { 
      if (gid_list != 0) gid_list = 0 ;  
      if (lid_list != 0) lid_list = 0 ;
      if (proc_list != 0) proc_list = 0 ;
      if (user_def_list != 0) delete [] user_def_list;
    }

    virtual unsigned int GetListLength() { return(proc_size) ; } 

    virtual int GetGIDListLength() { return(gid_size) ; }

    virtual unsigned int *GetAllGIDs() { return(gid_list) ; }

    virtual int GetLIDListLength() { return(lid_size) ; }

    virtual unsigned int *GetAllLIDs() { return(lid_list) ; } 

    virtual int *GetResidentProcsList() { return(proc_list) ; }

    virtual int GetExtraDataSize() { return(user_def_size) ; } 

    virtual unsigned int *GetExtraData() { return(user_def_list) ; } 

    virtual void addRef() { ref++ ; }

    virtual void deleteRef() { ref--; if (ref == 0) delete (this); }

    int create_gid_list(int n, unsigned int *glist) 
    { 
      gid_size = n ; gid_list = glist ;
    }

    int create_lid_list(int n, unsigned int *llist) 
    { 
      lid_size = n ; lid_list = llist ; 
    }

    int create_user_def_list(int n, unsigned int *udlist) 
    { 
      user_def_size = n ; user_def_list = udlist ;
    }

    int create_proc_list(int n, int *plist) 
    { 
      proc_size = n ; proc_list = plist ;
    }

    int zero_out_lists() 
    {
      gid_list = lid_list = 0 ;  proc_list = 0 ;
      if (user_def_list != 0) delete [] user_def_list ; user_def_list = 0 ;
      gid_size = lid_size = proc_size = user_def_size = 0 ;
    }

    private :

    int proc_size, gid_size, lid_size, user_def_size, ref ;
    unsigned int *gid_list, *lid_list, *user_def_list ;
    int *proc_list ;
  } ;
};
#endif
