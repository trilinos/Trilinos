/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#include "lb_const.h"
#include "octree_const.h"
#include "migoct_const.h"
#include "comm_const.h"
#include "all_allo_const.h"

/* function prototypes */

static int Zoltan_Oct_Update_Connections(ZZ *zz, pOctant *octs, int *newpids, pOctant *newocts, int nocts);

static int Zoltan_Oct_Final_Migration(ZZ *zz, pOctant *octs, int *newpids, pOctant *newocts, int nocts, int nrecocts);
static int Zoltan_Oct_Update_Map(ZZ *zz);

/* Use high tag numbers. The MPI standard says all tags should be
   <= 32767. Note that tags 32766-32767 are used elsewhere. */
#define MigOctCommCreate 32760
#define MigOctCommDo 32761
#define MigOctCommReverse 32762
#define MigUpdCommCreate 32763
#define MigUpdCommDo 32764
#define MigFinCommCreate 32765
#define MigFinCommDo 32759
#define RootListCommCreate 32758
#define RootListCommDo  32757

typedef struct
{
  int num;           /* requestor's number for this request */
  pOctant ptr;       /* pointer that requestor asked for */
} OCTNEW_msg;

typedef struct
{
  pOctant oct;                /* Update link of oct */
  int     childnum;           /* Update child (0-7) or parent (-1) */
  pOctant newptr;             /* New pointer for child or parent */
  int     newpid;             /* New pid for child or parent */
} Update_msg;

#define FILLUPDATEMSG(msg, octant, chnum, nptr, npid) { msg.oct = octant; \
                                                        msg.childnum = chnum; \
                                                        msg.newptr = nptr; \
                                                        msg.newpid = npid; }

typedef struct
{
  pOctant ptr;            /* New location */

  int ppid;
  pOctant parent;
  int childnum;
  int cpids[8];
  pOctant children[8];

  double min[3];
  double max[3];

  int id;              
  int dir;
  int mapidx;
  int from;
} Migrate_msg;

#define FILLMIGRATEMSG(oct, NEW_OCTANTPTR, msg, proc)  { msg.ptr=NEW_OCTANTPTR; \
                                            msg.parent = Zoltan_Oct_parent(oct); \
				            msg.ppid   = Zoltan_Oct_Ppid(oct); \
				            msg.childnum = Zoltan_Oct_childnum(oct);  \
				            msg.id = Zoltan_Oct_id(oct); \
				            msg.dir = Zoltan_Oct_dir(oct); \
				            Zoltan_Oct_children(oct,msg.children);  \
				            Zoltan_Oct_cpids(oct,msg.cpids);  \
				            Zoltan_Oct_bounds(oct,msg.min,msg.max); \
                                            msg.mapidx = Zoltan_Oct_mapidx(oct); \
                                            msg.from = proc; } 

#define SETOCTFROMMIGRATEMSG(OCT_info, msg)  { Zoltan_Oct_setID(msg.ptr,msg.id); \
                                               Zoltan_Oct_POct_setparent(OCT_info, msg.ptr, msg.parent, msg.ppid); \
                                               Zoltan_Oct_setchildnum(msg.ptr,msg.childnum); \
                                               Zoltan_Oct_setchildren(msg.ptr, msg.children, msg.cpids); \
                                               Zoltan_Oct_setbounds(msg.ptr, msg.min, msg.max); \
                                               Zoltan_Oct_setID(msg.ptr,msg.id); \
                                               Zoltan_Oct_setDir(msg.ptr,msg.dir); \
                                               Zoltan_Oct_setMapIdx(msg.ptr,msg.mapidx); }
      
int Zoltan_Oct_migrate_octants(ZZ *zz, int *newpids, pOctant *octs, int nocts, int *nrecocts) {
  int i,j = 0;
  int nsends = 0;
  int nreceives = 0;
  int *despid = NULL;
  OCTNEW_msg *snd_reply = NULL;
  OCTNEW_msg *rcv_reply = NULL;
  int ierr = ZOLTAN_OK;
  ZOLTAN_COMM_OBJ *comm_plan;        /* Object returned by communication routines */
  char *yo = "Zoltan_Oct_migrate_octants";
  pOctant *newocts = NULL;                          /* New foreign octant pointers */

  if((newocts = (pOctant *) ZOLTAN_MALLOC(sizeof(pOctant)*(nocts+10))) == NULL) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ZOLTAN_MEMERR;
  }

  for (i=0; i<nocts; i++)
    newocts[i]=NULL;


  /* count number of sends */

  nsends=0;
  for (i=0; i<nocts; i++)            
    if (newpids[i]!=zz->Proc) 
      nsends++;

  /* build message array */

/*   if(nsends > 0) { */
    if((snd_reply = (OCTNEW_msg *) ZOLTAN_MALLOC((nsends + 1) * sizeof(OCTNEW_msg))) == NULL) {
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }
    if((despid = (int *) ZOLTAN_MALLOC((nsends+10) * sizeof(int))) == NULL) {
      ZOLTAN_TRACE_EXIT(zz, yo);
      ZOLTAN_FREE(&snd_reply);
      return ZOLTAN_MEMERR;
    }
/*   } */
/*   else { */
/*     snd_reply = NULL; */
/*     despid = NULL; */
/*   } */

  j = 0;
  for (i=0; i<nocts; i++)                    
    if (newpids[i]!=zz->Proc) {  
      snd_reply[j].num = i;
      despid[j++] = newpids[i];
    }
   
  /* send messages */

  ierr = Zoltan_Comm_Create(&comm_plan, nsends, despid, zz->Communicator, MigOctCommCreate, &nreceives);
  if(ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    ZOLTAN_FREE(&snd_reply);
    ZOLTAN_FREE(&despid);
    return (ierr);
  }

  if((rcv_reply = (OCTNEW_msg *) ZOLTAN_MALLOC((nreceives + 1) * sizeof(OCTNEW_msg))) == NULL) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    ZOLTAN_FREE(&snd_reply);
    ZOLTAN_FREE(&despid);
    ZOLTAN_FREE(&rcv_reply);
    return ZOLTAN_MEMERR;
  }

  ierr = Zoltan_Comm_Do(comm_plan, MigOctCommDo, (char *) snd_reply,
		    sizeof(OCTNEW_msg), (char *) rcv_reply);
  if(ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    ZOLTAN_FREE(&snd_reply);
    ZOLTAN_FREE(&despid);
    ZOLTAN_FREE(&rcv_reply);
    return (ierr);
  }
    /* Reply to malloc requests and Receive malloc replies */
    
  for (i=0; i< nreceives; i++) {  
    rcv_reply[i].ptr = Zoltan_Oct_POct_new((OCT_Global_Info *) (zz->Data_Structure)); 
  }
;
  ierr = Zoltan_Comm_Do_Reverse(comm_plan, MigOctCommReverse, (char *) rcv_reply,
			    sizeof(OCTNEW_msg), NULL, (char *) snd_reply);
  if(ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    ZOLTAN_FREE(&snd_reply);
    ZOLTAN_FREE(&despid);
    ZOLTAN_FREE(&rcv_reply);
    return (ierr);
  }
  


  /* store remote pointers locally for future migration */
  for (i=0; i<nsends; i++) {                  
    newocts[snd_reply[i].num] = snd_reply[i].ptr;
  }
    
  ierr = Zoltan_Comm_Destroy(&comm_plan);
  if(ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    ZOLTAN_FREE(&snd_reply);
    ZOLTAN_FREE(&despid);
    ZOLTAN_FREE(&rcv_reply);
    return (ierr);
  }


  ZOLTAN_FREE(&snd_reply);
  ZOLTAN_FREE(&despid);
  ZOLTAN_FREE(&rcv_reply);

  /* set return value */

  *nrecocts = nreceives;
  ierr = Zoltan_Oct_Update_Connections(zz, octs, newpids, newocts, nocts);
  if(ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    abort();
    return (ierr);
  }
  ierr = Zoltan_Oct_Final_Migration(zz, octs,newpids,newocts,nocts, *nrecocts);
  if(ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    abort();
    return (ierr);
  }
  
  Zoltan_Oct_Update_Map(zz);

  ZOLTAN_FREE(&newocts);
  return ierr;
}

/*
 * update_connections(octs,newpids,newocts,nocts)
 *
 * Now that the new pid and pointer for each octant is known,
 * send that data to the parent and children so that they can
 * update their links.  Note: must delay handling of local pointers
 * until all sends are complete, otherwise necessary info wil
 * be overwritten.  After this function returns, all octants
 * will have correct pointers for the post-migrated locations,
 * and all that remains is to actually migrate the octants.
 * 
 */


static int Zoltan_Oct_Update_Connections(zz, octs,newpids,newocts,nocts)
ZZ *zz;
pOctant *octs;      /* octs[nocts]    */
int *newpids;       /* newpids[nocts] */
pOctant *newocts;   /* newocts[nocts] */
int nocts;          /* number of octants leaving this processor */
{
  int i, j;
  int nsends;
  int nreceives;
  pOctant parent;
  pOctant child;
  int ppid;
  int cpid;
  int childnum;
  int *despid = NULL;
  Update_msg umsg;
  Update_msg *localumsg = NULL;
  Update_msg *remoteumsg = NULL;
  Update_msg *rcv_umsg = NULL;
  int localcount;
  int remotecount;
  int ierr = ZOLTAN_OK;


  ZOLTAN_COMM_OBJ *comm_plan;           /* Object returned by communication routines */
  char *yo = "Zoltan_Oct_Update_Connections";
  OCT_Global_Info *OCT_info = (OCT_Global_Info *) zz->Data_Structure;
  localcount=0;
  remotecount=0;

  /* count number of sends */
  nsends = 0;
  for (i=0; i<nocts; i++)              
    if (newpids[i]!=zz->Proc)
      nsends++;

  if(nocts > 0) {
    if((remoteumsg = (Update_msg *) ZOLTAN_MALLOC((nocts+1) * sizeof(Update_msg)*9)) == NULL) {
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }
    
    if((localumsg  = (Update_msg *) ZOLTAN_MALLOC((nocts+1) * sizeof(Update_msg)*9)) == NULL) {
      ZOLTAN_TRACE_EXIT(zz, yo);
      ZOLTAN_FREE(&remoteumsg);
      return ZOLTAN_MEMERR;
    }
    
    if((despid = (int *) ZOLTAN_MALLOC((nocts+1) * sizeof(int)*9)) == NULL) {
      ZOLTAN_TRACE_EXIT(zz, yo);
      ZOLTAN_FREE(&remoteumsg);
      ZOLTAN_FREE(&localumsg);
      return ZOLTAN_MEMERR;
    }
  }
  else {
    remoteumsg = NULL;
    localumsg = NULL;
    despid = NULL;
  }
  localcount = 0;
  remotecount = 0;

  for (i=0; i<nocts; i++)                       /* Send connection updates */
    if (newpids[i]!=zz->Proc) {
	parent = Zoltan_Oct_parent(octs[i]); 
        ppid   = Zoltan_Oct_Ppid(octs[i]); 
        childnum = Zoltan_Oct_childnum(octs[i]);
	if (parent) {      /* Let parent of oct[i] know that it's moving   */
	  if (ppid==zz->Proc) {
	    FILLUPDATEMSG(localumsg[localcount], parent, childnum, newocts[i], newpids[i]);
	    localcount++;
	  }
	  else {
	    FILLUPDATEMSG(remoteumsg[remotecount], parent, childnum, newocts[i], newpids[i]);
	    despid[remotecount++] = ppid;
	  }
	}
	for (j=0; j<8; j++) {
	  child = Zoltan_Oct_child(octs[i],j);
	  cpid = octs[i]->cpid[j];
	  /* Tell child of oct[i] that it is moving */
	  if (child) {
	    if (cpid==zz->Proc) {
	      /* NOTE: -1 signals PARENT   */
	      FILLUPDATEMSG(localumsg[localcount], child, -1, newocts[i], newpids[i]);
	      localcount++;
	    }
	    else {
	      /* NOTE: -1 signals PARENT   */
	      FILLUPDATEMSG(remoteumsg[remotecount], child, -1, newocts[i], newpids[i]); 
	      despid[remotecount++] = cpid;
	    }
	  }
	}
    }

  ierr = Zoltan_Comm_Create(&comm_plan, remotecount, despid, zz->Communicator,
			MigUpdCommCreate, &nreceives);
  if(ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    ZOLTAN_FREE(&remoteumsg);
    ZOLTAN_FREE(&localumsg);
    ZOLTAN_FREE(&despid);
    return (ierr);
  }


/*   if(nreceives > 0) { */
    if((rcv_umsg = (Update_msg *) ZOLTAN_MALLOC((nreceives +1) * sizeof(Update_msg)*9)) == NULL) {
      ZOLTAN_TRACE_EXIT(zz, yo);
      ZOLTAN_FREE(&remoteumsg);
      ZOLTAN_FREE(&localumsg);
      ZOLTAN_FREE(&despid);
      return ZOLTAN_MEMERR;
    }

    
    ierr = Zoltan_Comm_Do(comm_plan, MigUpdCommDo, (char *) remoteumsg,
		      sizeof(Update_msg), (char *) rcv_umsg);
    if(ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_TRACE_EXIT(zz, yo);
      ZOLTAN_FREE(&remoteumsg);
      ZOLTAN_FREE(&localumsg);
      ZOLTAN_FREE(&despid);
      ZOLTAN_FREE(&rcv_umsg);
      return (ierr);
    }

/*   } */
/*   else { */
/*     rcv_umsg = NULL; */
/*   } */
  /* update new octants */
  for (i=0; i< (localcount+nreceives); i++)  {   
    if (i<localcount) 
      umsg=localumsg[i];
    else 
      umsg=rcv_umsg[i-localcount];
    if (umsg.childnum>=0) {
      Zoltan_Oct_setchild(umsg.oct,umsg.childnum,umsg.newptr);
      Zoltan_Oct_setCpid(umsg.oct,umsg.childnum,umsg.newpid);
    }
    else {
      if((Zoltan_Oct_data_newpid(umsg.oct) ==  OCT_info->OCT_localpid) ||
	 ((Zoltan_Oct_data_newpid(umsg.oct) !=  OCT_info->OCT_localpid) && (umsg.newpid == OCT_info->OCT_localpid)))
	Zoltan_Oct_POct_setparent(OCT_info, umsg.oct,umsg.newptr,umsg.newpid);
      else {
	umsg.oct->ppid = umsg.newpid;
	umsg.oct->parent = umsg.newptr;
      }
    }
  }

  ierr = Zoltan_Comm_Destroy(&comm_plan);
  if(ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    ZOLTAN_FREE(&remoteumsg);
    ZOLTAN_FREE(&localumsg);
    ZOLTAN_FREE(&despid);
    ZOLTAN_FREE(&rcv_umsg);
    return (ierr);
  }

  ZOLTAN_FREE(&remoteumsg);
  ZOLTAN_FREE(&localumsg);
  ZOLTAN_FREE(&rcv_umsg);
  ZOLTAN_FREE(&despid);
  return ierr;
}

/*
 * final_migration(octs,newpids,newocts,nocts,nrecocts)
 *
 * send updated octants to their new processors and delete them
 * locally
 *
 */

static int Zoltan_Oct_Final_Migration(zz, octs,newpids,newocts,nocts,nrecocts)
ZZ *zz;
pOctant *octs;      /* octs[nocts]    */
int *newpids;       /* newpids[nocts] */
pOctant *newocts;   /* newocts[nocts] */
int nocts;          /* number of octants leaving this processor */
int nrecocts;       /* number of octants received in this processor */
{
  int i;
  Migrate_msg *msnd = NULL, *mrcv = NULL;
  int remotecount;
  int nsends;
  int nreceives;
  int *despid = NULL;
  int ierr = ZOLTAN_OK;
  ZOLTAN_COMM_OBJ *comm_plan;           /* Object returned by communication routines */
  char *yo = "Zoltan_Oct_Final_Migration";
  OCT_Global_Info *OCT_info = (OCT_Global_Info *) zz->Data_Structure;

  /* count number of sends */
  nsends=0;
  for (i=0; i<nocts; i++)              
    if (newpids[i]!=zz->Proc)
      nsends++;

  if((despid = (int *) ZOLTAN_MALLOC((nocts+10) * sizeof(int))) == NULL) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ZOLTAN_MEMERR;
  }

  if((msnd   = (Migrate_msg *) ZOLTAN_MALLOC((nocts+10) * sizeof(Migrate_msg))) == NULL) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    ZOLTAN_FREE(&despid);
    return ZOLTAN_MEMERR;
  }

  remotecount = 0;
  for (i=0; i<nocts; i++)                          /* Send and free */
    if (newpids[i]!=zz->Proc)
      { 
	FILLMIGRATEMSG(octs[i], newocts[i], msnd[remotecount], zz->Proc); /* bug */
	despid[remotecount++] = newpids[i];
	Zoltan_Oct_clearRegions(octs[i]);
        /* KDDKDDFREE Change oct to &oct to allow NULL from ZOLTAN_FREE 
         * KDDKDDFREE to propagate back. */
	Zoltan_Oct_POct_free(OCT_info, &(octs[i]));
      }

  ierr = Zoltan_Comm_Create(&comm_plan, remotecount, despid, zz->Communicator,
			MigFinCommCreate, &nreceives);
  if(ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    ZOLTAN_FREE(&despid);
    ZOLTAN_FREE(&msnd);
    return (ierr);
  }

  if((mrcv = (Migrate_msg *) ZOLTAN_MALLOC((nreceives+10) * sizeof(Migrate_msg))) == NULL) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    ZOLTAN_FREE(&despid);
    ZOLTAN_FREE(&msnd);
    return ZOLTAN_MEMERR;
  }

  ierr = Zoltan_Comm_Do(comm_plan, MigFinCommDo, (char *) msnd,
		    sizeof(Migrate_msg), (char *) mrcv);
  if(ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    ZOLTAN_FREE(&despid);
    ZOLTAN_FREE(&msnd);
    ZOLTAN_FREE(&mrcv);
    return (ierr);
  }

  for (i=0; i<nreceives; i++) {                     /* Receive new parocts */
/*     Zoltan_Oct_setID(mrcv[i].ptr,mrcv[i].id); */
/*     Zoltan_Oct_POct_setparent(OCT_info, mrcv[i].ptr, mrcv[i].parent, mrcv[i].ppid);  */
    SETOCTFROMMIGRATEMSG(OCT_info, mrcv[i]); 
/*     Zoltan_Oct_setchildnum(mrcv[i].ptr,mrcv[i].childnum);  */
/*     Zoltan_Oct_setchildren(mrcv[i].ptr, mrcv[i].children, mrcv[i].cpids);  */
/*     Zoltan_Oct_setbounds(mrcv[i].ptr, mrcv[i].min, mrcv[i].max);   */
/*     Zoltan_Oct_setDir(mrcv[i].ptr,mrcv[i].dir); */
/*     Zoltan_Oct_setMapIdx(mrcv[i].ptr,mrcv[i].mapidx); */
  }

  ierr = Zoltan_Comm_Destroy(&comm_plan);
  if(ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    ZOLTAN_FREE(&despid);
    ZOLTAN_FREE(&msnd);
    ZOLTAN_FREE(&mrcv);
    return (ierr);
  }

  ZOLTAN_FREE(&despid);
  ZOLTAN_FREE(&msnd);
  ZOLTAN_FREE(&mrcv);

  return ierr;
}


static int Zoltan_Oct_build_global_rootlist(ZZ *zz,Migrate_msg  **ret_rmsg, int *size) {
  int j, k;
  int *despid = NULL;
  int nroots, nreceives;
  pRList  RootList;                  /* list of the local roots */
  pOctant RootOct;
  Migrate_msg *snd_rmsg = NULL;
  Migrate_msg *rcv_rmsg = NULL;
  OCT_Global_Info *OCT_info = (OCT_Global_Info *)(zz->Data_Structure);
/*Map *array = OCT_info->map;*/
  ZOLTAN_COMM_OBJ *comm_plan;                /* Object returned by communication routines */

  int ierr = ZOLTAN_OK;
  char *yo = "Zoltan_Oct_build_global_rootlist";
 
  nroots = RL_numRootOctants(Zoltan_Oct_POct_localroots(OCT_info));

  if((despid = (int *) ZOLTAN_MALLOC((zz->Num_Proc)*nroots * sizeof(int))) == NULL) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ZOLTAN_MEMERR;
  }

  if((snd_rmsg = (Migrate_msg *) ZOLTAN_MALLOC((zz->Num_Proc)*nroots * sizeof(Migrate_msg))) == NULL) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    ZOLTAN_FREE(&despid);
    return ZOLTAN_MEMERR;
  }
  
  k = 0;
  for (j=0; j<zz->Num_Proc; j++) {
    RootList = Zoltan_Oct_POct_localroots(OCT_info);
    while((RootOct = RL_nextRootOctant(&RootList))) {	
 /*      if(array[Zoltan_Oct_mapidx(RootOct)].npid > 0) { */
	FILLMIGRATEMSG(RootOct, RootOct, snd_rmsg[k], zz->Proc);
	despid[k] = j;
	k++;
/*       } */
    }
  }
  
  ierr = Zoltan_Comm_Create(&comm_plan, k, despid, zz->Communicator,
			RootListCommCreate, &nreceives);

  if(ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    ZOLTAN_FREE(&despid);
    ZOLTAN_FREE(&snd_rmsg);
    return (ierr);
  }

  if((rcv_rmsg = (Migrate_msg *) ZOLTAN_MALLOC(nreceives * sizeof(Migrate_msg))) == NULL) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    ZOLTAN_FREE(&despid);
    ZOLTAN_FREE(&snd_rmsg);
    return ZOLTAN_MEMERR;
  }

  

  ierr = Zoltan_Comm_Do(comm_plan, RootListCommDo, (char *) snd_rmsg, sizeof(Migrate_msg), (char *) rcv_rmsg);
  if(ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    ZOLTAN_FREE(&despid);
    ZOLTAN_FREE(&snd_rmsg);
    ZOLTAN_FREE(&rcv_rmsg);
    return (ierr);
  }

  ZOLTAN_FREE(&despid);
  ZOLTAN_FREE(&snd_rmsg);


  ierr = Zoltan_Comm_Destroy(&comm_plan);
  if(ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ierr);
  }

  *ret_rmsg = rcv_rmsg;
  *size = nreceives;

  return ierr;
}

static int Zoltan_Oct_Update_Map(ZZ *zz) {
  int i;
  double x,y;
  pRList  RootList;                
  pOctant RootOct;
  pOctant remoteoctant;
  Migrate_msg *rootlists = NULL;
  OCT_Global_Info *OCT_info = (OCT_Global_Info *)(zz->Data_Structure);
  Map *array = OCT_info->map;
  int mapsize = OCT_info->mapsize;
  int rlsize = 0;
  int ierr = ZOLTAN_OK;
  char *yo = "Zoltan_Oct_Update_Map";

  if((ierr = Zoltan_Oct_build_global_rootlist(zz, &rootlists, &rlsize)) != ZOLTAN_OK) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ierr;
  }


  for(i = 0; i < mapsize; i++) {
    RootList = array[i].list;
    while((RootOct = RL_nextRootOctant(&RootList)))  {
      Zoltan_Oct_free(OCT_info, &RootOct);
      /* KDDKDDFREE set oct pointer of RootList to NULL. */
      RootList->oct = NULL;
    }
    RL_clearRootOctants(&(array[i].list));
  }
  
  for(i = 0; i < rlsize; i++) {
    remoteoctant = Zoltan_Oct_newremote();
    remoteoctant->remoteptr = rootlists[i].ptr;
    x = rootlists[i].max[0] - rootlists[i].min[0];
    y = rootlists[i].max[1] - rootlists[i].min[1];
    remoteoctant->area = x*y;
    remoteoctant->ppid = rootlists[i].ppid;
    remoteoctant->npid = rootlists[i].from;
    Zoltan_Oct_setID(remoteoctant,rootlists[i].id); 
    Zoltan_Oct_setchildnum(remoteoctant,rootlists[i].childnum); 
    Zoltan_Oct_setchildren(remoteoctant, rootlists[i].children, rootlists[i].cpids); 
    Zoltan_Oct_setbounds(remoteoctant, rootlists[i].min, rootlists[i].max);  
    Zoltan_Oct_setDir(remoteoctant,rootlists[i].dir);
    Zoltan_Oct_setMapIdx(remoteoctant,rootlists[i].mapidx);
    RL_addRootOctant(array[Zoltan_Oct_mapidx(remoteoctant)].list, remoteoctant);	
  }

  if(rlsize > 0)
    ZOLTAN_FREE(&rootlists);
  return ierr;
}
