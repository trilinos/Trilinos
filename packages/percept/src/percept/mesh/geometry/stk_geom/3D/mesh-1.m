indexQuad[i_,j_,n_] := i + n j;
indexTri[i_,j_,n_] := ((n-j)(n - j + 1))/2 + 1 + i

edgeIndicesTri[iedge_ , n_] :=
    If[iedge == 0,
       Table[indexTri[i,0,n], {i,0,n}],
       If[iedge == 1,
          Table[indexTri[n-j,j,n], {j, 0, n}],
          Table[indexTri[0,j,n], {j, n, 0}]
         ]
      ]

edgeIndicesQuad[iedge_ , n_] :=
    If[iedge == 0,
       Table[indexQuad[i,0,n], {i,0,n}],
       If[iedge == 1,
          Table[indexQuad[n,j,n], {j, 0, n}],
          If[iedge == 2,
             Table[indexQuad[i,n,n], {i, n, 0}],
             If[iedge == 3,
                Table[indexQuad[0,j,n], {j, n, 0}]
               ]
            ]
         ]
      ]


FitCubics[mesh_] :=
    Module[{nodes = mesh[[1]], elements = mesh[[2]], normals = mesh[[3]]},
           ret = {};
           (* we double up the points to account for extra Gregory points *)
           Do[ ret = Join[ret, Table[{{0,0,0},{0,0,0}},Ncp[Length[elements[[ielem]] ] ] ] ], {ielem,Length[elements]}];

           Do[
               elemI = elements[[Ielem]];
               Do[
                   elemJ = elements[[Jelem]];
                   If[Ielem != Jelem,
                      (* find if the two elements share an edge *)
                      edge = Intersection[nodes[[elemI]], nodes[[elemJ]]  ];
                      If[Length[edge] == 2,
                         indexI = Position[elemI, edge[[1]] ][[1, 1]];
                         indexJ = Position[elemJ, edge[[2]] ][[1, 1]];
                         If[Length[elemI] == 3,
                            idsI = edgeIndicesTri[indexI - 1, 3];
                            idsJ = edgeIndicesTri[indexJ - 1, 3],

                            idsI = edgeIndicesQuad[indexI - 1, 3]
                            idsJ = edgeIndicesQuad[indexJ - 1, 3]
                           ]
                         rule = Flatten[ {  Thread[V[pi]->nodes[[edge[[1]] ]] ], Thread[V[pj]->nodes[[edge[[2]] ]] ],
                                            Thread[V[ni]->normals[[edge[[1]] ]] ], Thread[V[nj]->normals[[edge[[2]] ]] ] } ];
                         sol = solCubic /. rule;
                         sol = Partition[sol, 3];
                         Do[
                             Do[  ret[[ Ielem, idsI[[ KK ]] , 1, LL ]] = sol[[ KK , LL]] , {LL,3} ]
                             ,{KK,Length[idsI]}]

                         rule = Flatten[ {  Thread[V[pi]->nodes[[edge[[2]] ]] ], Thread[V[pj]->nodes[[edge[[1]] ]] ],
                                            Thread[V[ni]->normals[[edge[[2]] ]] ], Thread[V[nj]->normals[[edge[[1]] ]] ] } ];
                         sol = solCubic /. rule;
                         sol = Partition[sol, 3];
                         Do[
                             Do[  ret[[ Jelem, idsJ[[ KK ]] , 1, LL ]] = sol[[ KK , LL]] , {LL,3} ]
                             ,{KK,Length[idsJ]}]

                        ]
                     ]
                 ]
             ]
           ret
          ]
