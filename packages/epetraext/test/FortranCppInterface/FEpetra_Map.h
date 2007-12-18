/* True C header file! */


#ifndef FEPETRA_MAP_H
#define FEPETRA_MAP_H


#ifdef __cplusplus
extern "C" {
#endif


typedef int MapID;

MapID FEpetra_Map_Create( int numGlobalElements );

void FEpetra_Map_Destroy( MapID mapID );

int FEpetra_Map_NumGlobalElements( MapID mapID );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* FEPETRA_MAP_H */
