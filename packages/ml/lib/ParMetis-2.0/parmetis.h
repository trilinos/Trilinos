typedef int idxtype;

void ParMETIS_PartKway(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *, MPI_Comm *);
void ParMETIS_RefineKway(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, idxtype *, MPI_Comm *);
void ParMETIS_RepartLDiffusion(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, idxtype *, MPI_Comm *);
void ParMETIS_RepartGDiffusion(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, idxtype *, MPI_Comm *);
void ParMETIS_RepartRemap(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, idxtype *, MPI_Comm *);
void ParMETIS_RepartMLRemap(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, idxtype *, MPI_Comm *);
void ParMETIS_PartGraphGeomKway(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, int *, idxtype *, MPI_Comm *);
void ParMETIS_PartGraphGeomRefine(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *, MPI_Comm *);
void ParMETIS_PartGraphGeom(idxtype *, int *, float *, idxtype *, MPI_Comm *);
void ParMETIS_NodeND(idxtype *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *, MPI_Comm *);

void PARUAMETIS(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, MPI_Comm);
void PARDAMETIS(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, MPI_Comm);
void PARKMETIS(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, MPI_Comm);
void PARRMETIS(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, MPI_Comm);
void PAROMETIS(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, MPI_Comm);
void PARGKMETIS(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int, float *, idxtype *, int *, MPI_Comm);
void PARGRMETIS(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int, float *, idxtype *, int *, MPI_Comm);
void PARGMETIS(idxtype *, idxtype *, idxtype *, int, float *, idxtype *, int *, MPI_Comm);


