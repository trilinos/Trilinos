(* Notes: use iimage-mode to see embedded images; use "C-c n" to turn on removal of ^M from Mathematica output *)

(* ::Package:: *)

(* ::Section::Closed:: *)
(*Background	*)


(* Triangular Bezier patch by Garry Helzer *)
baryRules=Solve[{u{a1,b1}+v{a2,b2}+w{a3,b3}=={x,y},u+v+w==1},{u,v,w}] ;
BarycentricCoordinates[Polygon[{{a1_,b1_},{a2_,b2_},{a3_,b3_}}]][{x_,y_}] = {u,v,w} /. baryRules //Flatten ;
Subdivide[l_]:= l /. Polygon[{p_,q_,r_}] :> Polygon /@ ({{p+p,p+q,p+r},{p+q,q+q,q+r},{p+r,q+r,r+r},{p+q,q+r,r+p}}/2) ;
Transform[F_][L_] := L /. Polygon[l_] :> Polygon[F /@ l] ;
TriEval[ControlPts_][{u_,v_,w_}]:= Module[{x,y,z,n=(Sqrt[8 Length[ControlPts]+1]-3)/2},
  ((List @@ Expand[(x+y+z)^n]) /. {x->u,y->v,z->w}).ControlPts] ;
Param[Tri_,CP_][{x_,y_}]:=With[{p=BarycentricCoordinates[Tri][{x, y}]},TriEval[CP][p]]
(* Triangular bezier patch for n=3 *)
Tri=Polygon[{{1, 0}, {0, 1}, {0, 0}}];
CP={P300,P210,P120,P030, P201,P111,P021, P102,P012, P003} ={{3,0,0},{2.5,1,.5},{2,2,0},{1.5,3,0}, {2,0,1},{1.5,1,2},{1,2,.5}, {1,0,1},{.5,1,.5}, {0,0,0}};
SubT=Nest[Subdivide, Tri, 3];
Patch=Transform[Param[Tri, CP]][SubT];
cpts={PointSize[0.02], Point/@CP};
coord={AbsoluteThickness[1],Line/@{{{0,0,0},{3.2,0,0}},{{0,0,0},{0,3.4,0}},{{0,0,0},{0,0,1.3}}}};
cpolygon={AbsoluteThickness[2],Line[{P300,P210,P120,P030,P021,P012,P003,P102,P201,P300}],Line[{P012,P102,P111,P120,P021,P111,P201,P210,P111,P012}]};
Show[Graphics3D[{cpolygon,cpts,coord,Patch}], Boxed->False, PlotRange->All,ViewPoint->{2.620, -3.176, 2.236}]



pnts={{3,3,0}, {2,2,0},{4,2,1}, {1,1,0},{3,1,1},{5,1,2},
         {0,0,0},{2,0,1},{4,0,2},{6,0,3}};


BezierP[i_,j_,k_,n_,u_,v_,w_]:=(n!/(i! j! k!))u^i v^j w^k;
BezierP[i_,j_,n_,u_,v_,w_] := BezierP[i,j,n-i-j,n,u,v,w]
n=3; u=1/6; v=2/6; w=3/6; Tsrpt={0,0,0};
TriIndex[i_,j_,n_] := (n-j)(n-j+1)/2+1+i;
EvalTri[pnts_,n_,u_,v_]:=Module[{Tsrp={0,0,0},w=1-u-v,k},
  Do[{k=n-i-j, Tsrpt=Tsrpt+BezierP[i,j,k,n,u,v,w] pnts[[ TriIndex[i,j,n] ]] }, {j,0,n}, {i,0,n-j}];
  Tsrpt
]

EvalTri[pnts,n,u,v]

TriEval[pnts][{u,v,w}]

TriIndex[0,0,3]

Clear[pnts,u,v,w]


(* ::Section::Closed:: *)
(*More background*)


(* biquadratic bezier surface patch *)
Clear[u,v,pwr,bern,spnts,n,bzSurf,g1,g2]; n=2;

spnts={{{0,0,0},{1,0,1},{0,0,2}},{{1,1,0},{4,1,1},{1,1,2}}, {{0,2,0},{1,2,1},{0,2,2}}};

(* Handle Indeterminate condition *)
pwr[x_,y_]:=If[x==0 && y==0, 1, x^y];
bern[n_,i_,u_]:=Binomial[n,i]pwr[u,i]pwr[1-u,n-i] ;
bzSurf[u_,w_]:=Sum[bern[n,i,u] spnts[[i+1,j+1]] bern[n,j,w],{i,0,n}, {j,0,n}] ;
g1=ParametricPlot3D[bzSurf[u,w],{u,0,1}, {w,0,1}, Ticks->{{0,1,4},{0,1,2},{0,1,2}}, Compiled->False, DisplayFunction->Identity];
g2=Graphics3D[{AbsolutePointSize[3],Table[Point[spnts[[i,j]]],{i,1,n+1},{j,1,n+1}]}];
Show[g1,g2, ViewPoint->{2.783, -3.090, 1.243}, PlotRange->All, DisplayFunction->$DisplayFunction]


(* ::Section:: *)
(*Implementation / Literature	*)


(* ::Text:: *)
(*
file://images/Screenshot-1.png
*)

(* ::Subsection:: *)
(*Code Gen Helpers*)

powerRules:={Power[x_,2]->MyPow2[x],Power[x_,3]->MyPow3[x],Power[x_,-1]->MyInverse[x]};
Clear[MyPow]
MyPow[x_,2]:=MyPow2[x];
MyPow[x_,0]:=1;
MyPow[x_,1]:=x;
MyPow[x_,(n_/;IntegerQ[n]&&n<0)]:=MyPow[MyInverse[x],-n];
powerRules:={Power[x_,n_]:>MyPow[x,n]};
Formatter:=CForm;
IndentString:="      ";
ConvertRule:=(lhs_->rhs_):>"\n" <> IndentString<>ToString[Formatter[lhs]]<>" = "<>ToString[Formatter[rhs]]<>" ;";
CodeGen[listIn_,file_]:=Module[{list=listIn//.powerRules/.ConvertRule},Do[PutAppend[OutputForm[list[[i]]],file],{i,Length[list]}]];
CodeGen[listIn_]:=Module[{list=listIn//.powerRules/.ConvertRule},ColumnForm[list]];
SetOptions[Splice,FormatType->OutputForm,PageWidth->120];

(* ::Subsection:: *)
(*Degree Elevation*)

(* ::Text:: *)
(*
file://images/degree-elev.jpg
*)


DegreeElevateLine[b_,i_,n_,u_]:=i/(n+1) b[i-1][u]+(1-i/(n+1))b[i][u];
DegreeElevateLine[b_,i_,n_]:=i/(n+1) b[i-1]+(1-i/(n+1))b[i];
DegreeElevateTri[b_,i_,j_,n_,u_,v_]:=1/(n+1) (i b[i-1,j][u,v]+j b[i,j-1][u,v]+(n+1-i-j)b[i,j][u,v]);


DegreeElevateTri[P,#[[1]],#[[2]],2,u,v]&/@{{0,0},{1,0},{2,0},{3,0},{0,1},{1,1},{2,1}}//ColumnForm
eqqh=Table[qh[i]==DegreeElevateLine[q,i,3],{i,0,4}];
eqqhRule=Table[qh[i]->DegreeElevateLine[q,i,3],{i,0,4}];
Solve[Take[eqqh,4],Table[q[i],{i,0,3}]][[1]];
qDL:=Solve[Take[eqqh,4],Table[q[i],{i,0,3}]][[1]]/.q[i_]->V[q,i]/.qh[i_]->V[qh,i];
qDegreeLower=Flatten[Thread/@qDL];

(* ::Subsection:: *)
(*Quad Patch*)


(* ::Text:: *)
(*
file://images/Screenshot-10.png
*)

nq:=3


bt[u_,v_]:=Table[b[i,j][u,v],{j,0,nq},{i,0,nq}]


bt[u,v][[1,2]]


bt[u,v]//MatrixForm


btv[u_,v_]:=Table[Flatten[Table[b[i,j,k,l][u,v],{l,0,1},{k,0,1-l}]],{j,0,nq},{i,0,nq}]


btv[u,v]


Eb[u_,v_]:=Sum[b[i,j][u,v]36/(i!j!(3-i)!(3-j)!) (1-u)^(3-i) u^i (1-v)^(3-j) v^j,{i,0,nq},{j,0,nq}]


Ebv[u_,v_]:=Eb[uu,vv]/.{uu->u,vv->v}


ModB[u_,v_]:={
b[1,1][u,v]->u/(u+v) b[1,1,1,0][u,v]+v/(u+v) b[1,1,0,1][u,v],
b[2,1][u,v]->(1-u)/(1-u+v) b[2,1,1,0][u,v]+v/(1-u+v) b[2,1,0,1][u,v],
b[1,2][u,v]->u/(1-v+u) b[1,2,1,0][u,v]+(1-v)/(1-v+u) b[1,2,0,1][u,v],
b[2,2][u,v]->(1-u)/(2-u-v) b[2,2,1,0][u,v]+(1-v)/(2-u-v) b[2,2,0,1][u,v]
}


ModBtRHS[u_,v_]:={
{b[1,1,1,0][u,v],b[1,1,0,1][u,v]},
{b[2,1,1,0][u,v],b[2,1,0,1][u,v]},
{b[1,2,1,0][u,v],b[1,2,0,1][u,v]},
{b[2,2,1,0][u,v],b[2,2,0,1][u,v]}
}


ModBtLHS[u_,v_]:={
b[1,1][u,v],
b[2,1][u,v],
b[1,2][u,v],
b[2,2][u,v]
}


ModBt[u_,v_]:=bt[u,v]/.Thread[ModBtLHS[u,v]->ModBtRHS[u,v]]


ModB[u,v]


Clear[u,v,w]


Gpq[u_,v_]:=Eb[u,v]/.ModB[u,v]//Simplify


Gpq[u,v]


bcl[u_,v_]:=Flatten[ModBt[u,v]]


bcl[u,v]//TableForm


bcl[u,v]//CForm


spts:={{0,0,0},{1,0,0},{2,0,0},{3,0,0},
{0,1,0},{1,1,1},{1,1,1},{2,1,1},{2,1,1},{3,1,0},
{0,2,0},{1,2,1},{1,2,1},{2,2,1},{2,2,1},{3,2,0},
{0,3,0},{1,3,0},{2,3,0},{3,3,0}}


spts//Length


cpts:={PointSize[0.02], Point/@spts};
cpolygon:={AbsoluteThickness[2],Line[Take[spts,4]],Line[Take[spts,{5,10}]],Line[Take[spts,{11,16}]],Line[Take[spts,{17,20}]]};


plt=Gpq[u,v]/.Thread[bcl[u,v]->spts];


Show[ParametricPlot3D[plt,{u,0,1},{v,0,1}],Graphics3D[{cpts,cpolygon}],PlotRange->All]


spts:={{0,0,0},{1,0,0},{2,0,0},{3,0,0},
{0,1,0},{.5,0.2,1},{1.3,1.2,1},{2,1,1},{2,1,1},{3,1,0},
{0,2,0},{1,2,1},{1,2,1},{2,2,1},{2,2,1},{3,2,0},
{0,3,0},{1,3,0},{2,3,0},{3,3,0}}


spts:={{0,0,0},{1,0,0},{2,0,0},{3,0,0},
{0,1,0},{8.7,.7,1.4},{1.2,1.2,1},{2,1,1},{2,1,1},{3,1,0},
{0,2,0},{1,2,1},{1,2,1},{2,2,1},{2,2,1},{3,2,0},
{0,3,0},{1,3,0},{2,3,0},{3,3,0}}


plt=Gpq[u,v]/.Thread[bcl[u,v]->spts];


ParametricPlot3D[plt,{u,0,1},{v,0,1}]


bcl[u,v]


pts:=Table[Cp[i],{i,0,Length[bcl[u,v]]-1}]


Table[ToString["double Cpl[]"]<>"="<>ToString[InputForm[pts[[jc]]]]<>";"<>"\n",{jc,3}]//OutputForm>>>"quad_codegen.cpp"


EvalPatch[u_,v_]:=Simplify[Gpq[uu,vv]/.uu->u]/.{vv->v}/.Thread[bcl[u,v]->pts]


quadEvalList:={genUV->EvalPatch[u,v],
u0v->EvalPatch[0,v],
u1v->EvalPatch[1,v],
uv0->EvalPatch[u,0],
uv1->EvalPatch[u,1]
}


EvalPatch[u,v]


EvalPatchQuad[u_,v_]=EvalPatch[u,v]


bclQuad[u_,v_]=bcl[u,v]


Splice["eval_quad.mcpp","eval_quad.cpp"]


CodeGen[quadEvalList,"quad.cpp"]


(* ::Text:: *)
(*Image[CompressedData["*)


(* ::Subsection:: *)
(*Tri Patch*)


(* ::Text:: *)
(*Image[CompressedData["*)


nq:=4


bt[u_,v_]:=Table[b[i,j][u,v],{j,0,nq},{i,0,nq-j}]


bt[u,v][[1,2]]


bt[u,v]//MatrixForm


Eb[u_,v_]:=Sum[b[i,j][u,v]24/(i!j!(nq-i-j)!) u^i v^j (1-u-v)^(nq-i-j),{j,0,nq},{i,0,nq-j}]


Ebv[u_,v_]:=Eb[uu,vv]/.{uu->u,vv->v}


Ebv[uu,vv]


(* ::Text:: *)
(*Image[CompressedData["*)


ModB[u_,v_]:={
b[1,1][u,v]->u/(u+v) b[1,1,1,0][u,v]+v/(u+v) b[1,1,0,1][u,v],
b[2,1][u,v]->(1-u-v)/(1-u) b[2,1,1,0][u,v]+v/(1-u) b[2,1,0,1][u,v],
b[1,2][u,v]->u/(1-v) b[1,2,1,0][u,v]+(1-u-v)/(1-v) b[1,2,0,1][u,v]
}


ModBtRHS[u_,v_]:={
{b[1,1,1,0][u,v],b[1,1,0,1][u,v]},
{b[2,1,1,0][u,v],b[2,1,0,1][u,v]},
{b[1,2,1,0][u,v],b[1,2,0,1][u,v]}
}


ModBtLHS[u_,v_]:={
b[1,1][u,v],
b[2,1][u,v],
b[1,2][u,v]
}


ModBt[u_,v_]:=bt[u,v]/.Thread[ModBtLHS[u,v]->ModBtRHS[u,v]]


ModB[u,v]


ModBt[u,v]


Gpq[u_,v_]:=Eb[u,v]/.ModB[u,v]//Simplify


Gpq[u,v]


bcl[u_,v_]:=Flatten[ModBt[u,v]]


bcl[u,v]//TableForm


spts:={{0,0,0},{1,0,0},{2,0,0},{3,0,0},{4,0,0},
{0,1,0},{1,1,1},{1,1,1},{2,1,1},{2,1,1},{3,1,0},
{0,2,0},{1,2,1},{1,2,1},{2,2,0},
{0,3,0},{1,3,0},
{0,4,0}}


spts//Length


cpts:={PointSize[0.02], Point/@spts};
cpolygon:={AbsoluteThickness[2],Line[Take[spts,5]],Line[Take[spts,{6,11}]],Line[Take[spts,{12,15}]],Line[Take[spts,{16,17}]]};


xx=Gpq[u,v]/.Thread[bcl[u,v]->spts]//Simplify


g1=ParametricPlot3D[xx,{v,0,1},{u,0,1-v}];


Show[g1,Graphics3D[{cpts,cpolygon}],PlotRange->All]


spts:={{0,0,0},{1,0,0},{2,0,0},{3,0,0},{4,0,0},
{0,1,0},{0.5,0.5,1},{1.5,1.5,3},{2,1,3},{2,1,1},{3,1,0},
{0,2,0},{1,2,4},{1,2,1},{2,2,0},
{0,3,0},{1,3,0},
{0,4,0}}


xx=Gpq[u,v]/.Thread[bcl[u,v]->spts]//Simplify;


xx/.u->0.999/.v->0





g1=ParametricPlot3D[xx,{v,0,1},{u,0,1-v}];


Show[g1,Graphics3D[{cpts,cpolygon}],PlotRange->All]


EvalPatch[u_,v_]:=Simplify[Gpq[uu,vv]/.uu->u]/.{vv->v}/.Thread[bcl[u,v]->pts]


EvalPatch[u,v]


EvalPatch[u,v]/.{v->1-u}//Simplify


list:={genUV->EvalPatch[u,v],
u0v->EvalPatch[0,v],
vEqual1Mu->EvalPatch[u,1-u],
uv0->EvalPatch[u,0]
}


bclTri[u_,v_]=bcl[u,v]


EvalPatchTri[u_,v_]=EvalPatch[u,v]


Splice["eval_tri.mcpp","eval_tri.cpp"]


CodeGen[list,"tri.cpp"]


(* ::Text:: *)
(*Image[CompressedData["*)


(* ::Subsection:: *)
(*Enforcing G^1 Conditions*)


(* ::Text:: *)
(*Image[CompressedData["*)


(* ::Text:: *)
(*Image[CompressedData["*)


A:={{-3(1-\[Lambda]0),-3\[Lambda]0,0,0},
{(1-\[Lambda]1),\[Lambda]1,(1-\[Lambda]0),\[Lambda]0},
{0,0,-3(1-\[Lambda]1),-3\[Lambda]1}}


A//MatrixForm


greekRules:={\[Lambda]0->Lam0,\[Lambda]1->Lam1,\[Mu]0->Mu0,\[Mu]1->Mu1}


AA:=A/.greekRules


AA


MakeVec[p_,i_]:=Table[p[i,j],{j,0,2}]


MakeVec[p_]:=Table[p[j],{j,0,2}]


V[p_,i_]:=MakeVec[p,i]


V[p_]:=MakeVec[p]


V[p,0]


X:={V[p,1],V[r,1],V[p,2],V[r,2]}


A.X//Dimensions


B:={((1-\[Lambda]1)V[p,0]+\[Lambda]1 V[r,0])-((1-\[Mu]1)V[qh,0]+\[Mu]1 V[qh,1])-3((1-\[Mu]0)V[qh,1]+\[Mu]0 V[qh,2]),
((1-\[Mu]1)V[qh,1]+\[Mu]1 V[qh,2])+((1-\[Mu]0)V[qh,2]+\[Mu]0 V[qh,3]),
((1-\[Lambda]0)V[p,3]-\[Lambda]0 V[r,3]) - 3((1-\[Mu]1)V[qh,2]+\[Mu]1 V[qh,3]) - ((1-\[Mu]0) V[qh,3]+\[Mu]0 V[qh,4])
}


BB:=B/.greekRules


eq0:=(1-\[Lambda]0)V[p,0]+\[Lambda]0 V[r,0]-(1-\[Mu]0)V[q,0]+\[Mu]0 V[q,1]


eq1:=(1-\[Lambda]1)V[p,3]+\[Lambda]1 V[r,3]-(1-\[Mu]1)V[q,2]+\[Mu]1 V[q,3]


sol0:=Solve[Take[eq0,2]=={0,0},{\[Lambda]0,\[Mu]0}]//Simplify


sol1:=Solve[Take[eq1,2]=={0,0},{\[Lambda]1,\[Mu]1}]//Simplify


mat0:=Transpose[{Coefficient[eq0,\[Lambda]0],Coefficient[eq0,\[Mu]0]}]


mat1:=Transpose[{Coefficient[eq1,\[Lambda]1],Coefficient[eq1,\[Mu]1]}]


rhs0:=-(eq0-mat0.{\[Lambda]0,\[Mu]0})//Simplify


rhs1:=-(eq1-mat1.{\[Lambda]1,\[Mu]1})//Simplify


m0t:=Table[m0[i,j],{i,0,2},{j,0,1}]


r0t:=Table[R0[i],{i,0,2}]


m1t:=Table[m1[i,j],{i,0,2},{j,0,1}]


r1t:=Table[R1[i],{i,0,2}]


solRhs0:=Inverse[Transpose[m0t].m0t].Transpose[m0t].r0t//Simplify


solRhs1:=Inverse[Transpose[m1t].m1t].Transpose[m1t].r1t//Simplify


Ro10:=Thread[Flatten[Table[m0[i,j],{i,0,1},{j,0,2}]]->Flatten[mat0]]


Ro11:=Thread[Flatten[Table[m1[i,j],{i,0,1},{j,0,2}]]->Flatten[mat1]]


Ro12:=Thread[r0t->rhs0]


Ro12


Ro13:=Thread[r1t->rhs1]


Ro13


Ro20:=Thread[{\[Lambda]0,\[Mu]0}->solRhs0]/.greekRules


Ro21:=Thread[{\[Lambda]1,\[Mu]1}->solRhs0]/.greekRules


Ro20


(* ::Text:: *)
(*Image[CompressedData["*)


(* given pi, pj and ni, nj, find cubic *)


cubicFitRules:=Flatten[Thread/@{
V[ci]->(2V[pi]+V[pj])/3,
V[cj]->(V[pi]+2V[pj])/3,
V[c,0]->V[pi],
V[c,1]->(V[ci]-V[pi]).V[ni],
V[c,2]->(V[cj]-V[pj]).V[nj],
V[c,3]->V[pj]}]


cubicFitRules


ptriRules:=Flatten[Thread/@{
V[pt,0]->(V[q,0]+3V[p,0])/4,
V[pt,3]->(V[q,3]+3V[p,3])/4
}]


rtriRules:=Flatten[Thread/@{
V[rt,0]->(V[q,0]+3V[r,0])/4,
V[rt,3]->(V[q,3]+3V[r,3])/4
}]


pquadRules:=Flatten[Thread/@{
V[pt,0]->V[p,0],V[pt,3]->V[p,3]}]


rquadRules:=Flatten[Thread/@{
V[rt,0]->V[r,0],V[rt,3]->V[r,3]}]


estRules:=Flatten[Thread/@{
V[pe,1]->V[q,1]+2(V[pt,0]-V[q,0])/3+(V[pt,3]-V[q,3])/3,
V[pe,2]->V[q,2]+(V[pt,0]-V[q,0])/3+2(V[pt,3]-V[q,3])/3,
V[re,1]->V[q,1]+2(V[rt,0]-V[q,0])/3+(V[rt,3]-V[q,3])/3,
V[re,2]->V[q,2]+(V[rt,0]-V[q,0])/3+2(V[rt,3]-V[q,3])/3}]


(* ::Text:: *)
(*Image[CompressedData["*)


Xe:={V[pe,1],V[pe,2],V[re,1],V[re,2]}


A.Xe


Dimensions[A]


ALhs:=Table[a[i,j],{i,0,2},{j,0,3}]


ALhs//Dimensions


AAtLhs:=Table[aat[i,j],{i,0,2},{j,0,2}]


AAtInvLhs:=Table[aatI[i,j],{i,0,2},{j,0,2}]


AAtRhs:=ALhs.Transpose[ALhs]


DetAAtRhs:=Det[AAtLhs]


AAtInvRhs:=Simplify[DetAAtRhs Inverse[AAtLhs]]/DetAAtLhs


AAtInvRhs


Ro1:=Join[Thread[Flatten[AAtLhs]->Flatten[AAtRhs]],{DetAAtLhs->DetAAtRhs},Thread[Flatten[AAtInvLhs]->Flatten[AAtInvRhs]]]


powerRules


Thread[Flatten[AAtInvLhs]->Flatten[AAtInvRhs]][[1]]/.powerRules


Ro1


Ro2:=Thread[Flatten[ALhs]->Flatten[AA]]


BLhs:=Table[b[i,j],{i,0,2},{j,0,2}]


Ro3:=Thread[Flatten[BLhs]->Flatten[BB]]


DLhs:=Table[d[i,j],{i,0,2},{j,0,2}]


D1Lhs:=Table[d1[i,j],{i,0,2},{j,0,2}]


D1Rhs:=(BB-ALhs.Xe)//Simplify


DRhs:=AAtInvLhs.D1Lhs//Simplify


DRhs


Ro40:=Join[Thread[Flatten[D1Lhs]->Flatten[D1Rhs]],Thread[Flatten[DLhs]->Flatten[DRhs]]]


Ro40


?DegreeElevateLine


degQ:=Table[DegreeElevateLine[q,i,3],{i,0,4}]/.q[i_]->V[q,i]


Ro30:=Thread[Flatten[Table[V[qh,i],{i,0,4}]]->Flatten[degQ]]


Ro4


Dimensions[DLhs]


xRhs:=Transpose[ALhs].DLhs+Xe


xRhs


X


Ro50:=Thread[Flatten[X]->Flatten[xRhs]]


Ro50


cgList:={triRules,quadRules,estRules,Ro30,Ro10,Ro11,Ro12,Ro13,Ro20,Ro21,Ro3,Ro2,Ro1,Ro40,Ro50}


{triRules,quadRules,estRules};


cgList//ColumnForm;


CodeGen[#,"codegen.cpp"]&/@cgList


Splice["codegen.mcpp","codegen.cpp"]

(* ::Section:: *)
(*Mesh	*)
