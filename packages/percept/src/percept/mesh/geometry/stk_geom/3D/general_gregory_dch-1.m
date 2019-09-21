debug = False
(* ::Package:: *)

(* ::Subtitle:: *)
(*G1 Gregory Patches (agnostic)*)


(* ::Section:: *)
(**)
(*Graphics Directives*)


(* ::Input:: *)

(* Viewpoint *)
(* gf had 0, -1, 2 *)
(* myviewpoint = {0.2, -1.2, 2.0}; *) (* for sphere usedata = 5 *)
(* myviewpoint = {1.3, -2.4, 2.2}; *) (* sinxsiny usedata = 8 *)
 myviewpoint = {-0.1, -2, 3.1}; (* tri&rec patch usedata = 6 *)
(* myviewpoint = {4, -1,0}; *) (* for cone usedata = 9 *)

myviewvertical = {0,0,1}  (* default -- unless specified, us this one for data sets*)
(* myviewvertical = {0,-1,0} *)  (* cone data set for upright view *)

(* Turn on/off the control points on the shaded images *)
withcontrolstructure =1;

(* rendermethod = 0 smooth shading, 1 = colormap *)
rendermethod = 0;

(* This should be the number of evals -- but not sure *)
myplotpoints = 5;

colorsurface = RGBColor[1,.8,.2];
myopacity = 1.0; (* was 0.9 *)
myspecularity = 10;

(* Create two colors for the rectangle patch *)
(* One for main control net and one for extra Gregory points *)
colorrec1 = RGBColor[1, 0.5, 0.3];
colorrec2 = RGBColor[1, 0.8, 0.5];
colortri1 = RGBColor[0.7, 0.8, 0.3];
colortri2 = RGBColor[0.8, 0.9, 0.6];
colorpolygon = Gray;
colordatapoints = Black;
colornormals = Black;
colorbndry = Black

pointsizemedium = 0.025;
pointsizelarge = 0.03;

linewidththin = 0.003;
linewidthmedium = 0.004;
linewidththick = 0.007;

normallengthfactor = 1.0;



(* ::Section:: *)
(**)
(*Input Data*)


(* ::Input:: *)

(SetDirectory[ToFileName[Extract["FileName"/.NotebookInformation[EvaluationNotebook[]],{1},FrontEnd`FileName]]];)

(* Flag to check G1 continuity. *)
(* Added this with thought that it would speed up calculations, but it doesn't *)
checkg1condition = 1;
(* check occurs at t=0.5 -- this is same point for opposing patches. *)


(* Initialize dimension of bb1 -- the working Bezier patch *)
bb1=Table[{0,0,0},{i,1,4},{j,1,4}];
(* And the Bezier patch for the triangle data *)
tt1 = Table[{0,0,0},{i,1,5},{j,1,5}];

(* ================== Data Sets ===================== *)(* Convention:faces are ordered counterclockwise.
Neighbors are given across an edge.
faces lists the number of vertices,index into vertex list,and neighbors
*)

(* usedata set
0 = simple two rectangles
1 = two rectangles
2 = 2x2 rectangular grid
3 = 3x3 rectangular grid
4 = one triangle
5 = spherical cap
6 = one triangle and one rectangle
7 = two triangles
8 = 3x3 grid with triangle in middle
9 = cone
*)
usedataset = 9;

(* -----------Data Set: Simple Two rectangle------------- *)

If[usedataset == 0, (
data=Table[{0,0,0},{i,1,6}];
ndata=Length[data];

data[[1]]={0,0.0,0};
data[[2]]={1.0,0,0};
data[[3]]={1,1,0};
data[[4]]={0,1,0};
data[[5]]={0,-1,0};
data[[6]]={1,-1,0};

nn=Table[{0,0,1},{i,1,ndata}];


faces={
{4,{1,2,3,4},{2,0,0,0}},
{4,{6,2,1, 5},{0,1,0,0}}
};
nfaces=Length[faces];

)];

(* -----------Data Set: Two rectangle as spaceship ------------- *)

If[usedataset == 1, (
data=Table[{0,0,0},{i,1,6}];
ndata=Length[data];

data[[1]]={-1.5,0.0,0};
data[[2]]={1.5,0,0};
data[[3]]={1,1,0};
data[[4]]={0,1,0};
data[[5]]={0,-1,0};
data[[6]]={1,-1,0};

(* nn=Table[{0,0,1},{i,1,ndata}]; *)

nn=Table[(data[[i]]-{0, 0.5, -2})/Norm[(data[[i]]-{0, 0.5, -2})],{i,1,ndata}];



faces={
{4,{1,2,3,4},{2,0,0,0}},
{4,{2,1,5,6},{1,0,0,0}}
};
nfaces=Length[faces];

)];

(*------------Data Set:2x2 rectangles--------------*)

If[usedataset == 2,(
data=Table[{0,0,0},{i,1,9}];
	ndata=Length[data];

Do[data[[1+j]]={j,2,0},{j,0,2}];
	Do[data[[4+j]]={j,1,0},{j,0,2}];
Do[data[[7+j]]={j,0,0},{j,0,2}];


data[[5]] = data[[5]] + {0.5, 0.0, 0.25};


(* nn=Table[{0,0,1},{i,1,ndata}]; *)
nn=Table[(data[[i]]-{1.5, 1.5, -2})/Norm[(data[[i]]-{1.5, 1.5, -2})],{i,1,ndata}];


faces={{4,{1,4,5,2},{0,3,2,0}},
            {4,{2,5,6,3},{1,4,0,0}},
	  {4,{4,7,8,5},{0,0,4,1}},
	  {4,{5,8,9,6},{3,0,0,2}}};
nfaces = Length[faces];

)];


(*------------Data Set:3x3 rectangles--------------*)

If[usedataset == 3,(
data=Table[{0,0,0},{i,1,16}];
	ndata=Length[data];

Do[data[[1+j]]={j,3,0},{j,0,3}];
	Do[data[[5+j]]={j,2,0},{j,0,3}];
Do[data[[9+j]]={j,1,0},{j,0,3}];
	Do[data[[13+j]]={j,0,0},{j,0,3}];



movez = 0.5;

data[[6]] = data[[6]] + {0.3,0.3, movez};
data[[7]] = data[[7]] + {0, 0,- movez};
data[[10]] = data[[10]] + {0.0, 0.0, -movez};
data[[11]] = data[[11]] + {0, 0, movez};

nn=Table[{0,0,1},{i,1,ndata}];
(*nn=Table[(data[[i]]-{1.5, 1.5, -2})/Norm[(data[[i]]-{1.5, 1.5, -2})],{i,1,ndata}]; *)


faces={{4,{5,6,2,1},{4,2,0,0}},
            {4,{2,6,7,3},{1,5,3,0}},
	  {4,{3,7,8,4},{2,6,0,0}},
	  {4,{5,9,10,6},{0,7,5,1}},
	  {4,{6,10,11,7},{4,8,6,2}},
	  {4,{7,11,12,8},{5,9,0,3}},
	  {4,{9,13,14,10},{0,0,8,4}},
	  {4,{10,14,15,11},{7,0,9,5}},
	  {4,{11,15,16,12},{8,0,0,6}}};
nfaces = Length[faces];

)];

(* -----------Data Set: One triangle  ------------- *)

If[usedataset == 4, (
data=Table[{0,0,0},{i,1,3}];
ndata=Length[data];

data[[1]]={0,0.0,0.5};
data[[2]]={1.0,0,0.0};
data[[3]]={0,1,0};

nn=Table[{0,0,1},{i,1,ndata}];

faces={
{3,{1,2,3},{0,0,0}}
};
nfaces=Length[faces];

)];

(* -----------Data Set: Spherical Cap ------------- *)

If[usedataset == 5, (
data=Table[{0,0,0},{i,1,7}]; (* initialize *)
ndata=Length[data];

data[[1]]={0,2,0}//N;

(* -0.5 *)
cc={0,0,-0.25}; (* center of sphere that data should be on *)
radius=Norm[data[[1]]-cc]; (* radius of that sphere *)
data[[6]]= cc+{0,0,radius};
data[[7]]=cc+radius*({0,-0.5,0}-cc)/Norm[{0,-0.5,0}-cc];

deg=72*Pi/180//N;
rotmat={{Cos[deg],-Sin[deg],0},
	{Sin[deg], Cos[deg],0},
	{0,0,1}};

Do[data[[i+1]]=rotmat.data[[i]],{i,1,4}]; (* create pentagon *)


(* 3D unit normals come from sphere with center cc *)
nn=Table[(data[[i]]-cc)/Norm[(data[[i]]-cc)],{i,1,7}];

arrows=Table[{data[[i]],data[[i]]+nn[[i]]},{i,1,7}]; (* for Graphics *)


(* gf's original *)
(*
faces={
{4,{1,2,6,5},{0,2,5,0}},
{4,{2,3,7,6},{0,3,5,1}},
{3,{3,4,7},{0,4,2}},
{3,{7,4,5},{3,0,5}},
{3,{5,6,7},{1,2,4}}
};
*)

(* slightly different connectivity for triangles *)

faces={
{4,{1,2,6,5},{0,2,5,0}},
{4,{2,3,7,6},{0,3,4,1}},
{3,{3,4,7},{0,4,2}},
{3,{7,4,6},{3,5,2}},
{3,{4,5,6},{0,1,4}}
};

(* same as previous but with e2 to outside if possible *)
(*
faces={
{4,{1,2,6,5},{0,2,5,0}},
{4,{2,3,7,6},{0,3,4,1}},
{3,{7,3,4},{2,0,4}},
{3,{4,6,7},{5,2,3}},
{3,{6,4,5},{4,0,1}}
};
*)

(* all triangles *)
(*
faces={
{3,{1,2,6},{0,2,7}},
{3,{2,7,6},{3,6,1}},
{3,{2,3,7},{0,4,2}},
{3,{3,4,7},{0,5,3}},
{3,{7,4,5},{4,0,6}},
{3,{7,5,6},{5,7,2}},
{3,{5,1,6},{0,1,6}}
};
*)

(* just two triangles *)

(*
faces={
{3,{3,4,7},{0,2,0}},
{3,{7,4,5},{1,0,0}}
};
*)
(*just rectangles *)
(*
faces={
{4,{2,3,7,6},{0,0,2,3}},
{4,{4,5,6,7},{0,3,1,0}},
{4,{1,2,6,5},{0,1,2,0}}
};
*)

nfaces=Length[faces];

)];

(* -----  End of Spherical cap data set -----  *)
(* -----------Data Set: One triangle and one rectangle ------------- *)

If[usedataset == 6, (
data=Table[{0,0,0},{i,1,5}];
ndata=Length[data];

data[[1]]={0,0.0,0.0};
data[[2]]={1.0,0,0.25};
data[[3]]={0,1,0.25};
data[[4]]={1.5,0.5,0};
data[[5]]={1,1.5,0.0};

cc={0.75,0.75,-0.5};
(* nn=Table[{0,0,1},{i,1,ndata}]; *)
nn=Table[(data[[i]]-cc)/Norm[(data[[i]]-cc)],{i,1,ndata}];

(*
faces={
{3,{1,2,3},{0,2,0}},
{4,{2,4,5,3}, {0, 0, 0, 1}}
};
*)
faces={
{3,{2,3,1},{2,0,0}},
{4,{2,4,5,3}, {0, 0, 0, 1}}
};

nfaces=Length[faces];

)];

(* -----------Data Set: Two triangles ------------- *)

If[usedataset == 7, (
data=Table[{0,0,0},{i,1,4}];
ndata=Length[data];

data[[1]]={0.0,0.0,0.0};
data[[2]]={1.0,0,0.0};
data[[3]]={1,1,0.0};
data[[4]]={0,1,0};

(* nn=Table[{0,0,1},{i,1,ndata}]; *)

nn=Table[(data[[i]]-{0.5, 0.5, -2.0})/Norm[(data[[i]]-{0.5, 0.5, -2.0})],{i,1,ndata}];


faces={
{3,{2,4,1},{2,0,0}},
{3,{4,2,3},{1,0,0}}
};
nfaces=Length[faces];

)];

(*------------Data Set:3x3 rectangles with middle rec as 2 triangles --------------*)

If[usedataset == 8,(
data=Table[{0,0,0},{i,1,16}];
	ndata=Length[data];

(* Data over (-3,-3) to (3,3) grid *)

Do[data[[1+j]]={-3+(j/3)*6,3,0},{j,0,3}];
	Do[data[[5+j]]={-3+(j/3)*6,1,0},{j,0,3}];
Do[data[[9+j]]={-3+(j/3)*6,-1,0},{j,0,3}];
	Do[data[[13+j]]={-3+(j/3)*6,-3,0},{j,0,3}];

(* function x,y, sinx*siny *)

 Do[data[[j,3]] = Sin[data[[j,1]]]*Sin[data[[j,2]]], {j,1,ndata}];

nn=Table[({-Sin[data[[i,2]]]*Cos[data[[i,1]]], -Sin[data[[i,1]]]*Cos[data[[i,2]]], 1.0})/Norm[{-Sin[data[[i,2]]]*Cos[data[[i,1]]], -Sin[data[[i,1]]]*Cos[data[[i,2]]], 1.0}],{i,1,ndata}];


(*
faces={{4,{5,6,2,1},{4,2,0,0}},
            {4,{2,6,7,3},{1,6,3,0}},
	  {4,{3,7,8,4},{2,7,0,0}},
	  {4,{5,9,10,6},{0,8,5,1}},
	  {3,{6,10,11},{4,9,6}},
	  {3,{11,7,6},{7,2,5}},
            {4,{7,11,12,8},{6,10, 0, 3}},
	  {4,{9,13,14,10},{0,0,9,4}},
	  {4,{10,14,15,11},{8,0,10,5}},
	  {4,{11,15,16,12},{9,0,0,7}}};
nfaces = Length[faces];
*)

faces={{4,{5,6,2,1},{4,2,0,0}},
            {4,{2,6,7,3},{1,5,3,0}},
	  {4,{3,7,8,4},{2,7,0,0}},
	  {4,{5,9,10,6},{0,8,5,1}},
	  {3,{6,10,7},{4,6,2}},
	  {3,{11,7,10},{7,5,9}},
            {4,{7,11,12,8},{6,10, 0, 3}},
	  {4,{9,13,14,10},{0,0,9,4}},
	  {4,{10,14,15,11},{8,0,10,6}},
	  {4,{11,15,16,12},{9,0,0,7}}};
nfaces = Length[faces];

)];

(*------------Data  around  --------------*)

If[usedataset == 9,(
data=Table[{0,0,0},{i,1,18}];
	ndata=Length[data];

(* basic ring is at 2. Needs to be greater than 1 *)
r2 = 1.7;

data[[1]] = {-1,1,r2};
data[[2]] = {-1,-1,2};
data[[3]] = {1,-1,2};
data[[4]] = {1,1,r2};

data[[5]] = {-r2,1,1};
data[[6]] = {-2,-1,1};
data[[7]] = {-2,-1,-1};
data[[8]] = {-r2,1,-1};

data[[9]] = {-1,-1,-2};
data[[10]] = {1,-1,-2};
data[[11]] = {1,1,-r2};
data[[12]] = {-1,1,-r2};

data[[13]] = {2,-1,-1};
data[[14]] = {r2,1,-1};
data[[15]] = {r2,1,1};
data[[16]] = {2,-1,1};

data[[17]] = {0, -4, 0};
data[[18]] = {0, 4, 0};

cc={0,0,0};
nn=Table[(data[[i]]-cc)/Norm[(data[[i]]-cc)],{i,1,ndata}];

(* just around *)
(*
faces={{4,{1,2,3,4},{8,0,2,0}},
            {4,{3,16,15,4},{0,3,0,1}},
{4,{13,14,15,16},{4,0,2,0}},
{4,{11,14,13,10},{0,3,0,5}},
{4,{11,10,9,12},{4,0,6,0}},
{4,{12,9,7,8},{5,0,7,0}},
{4,{7,6,5,8},{0,8,0,6}},
	  {4,{6,2,1,5},{0,1,0,7}}};
*)

(* with one end of trinagles *)

faces={{4,{1,2,3,4},{8,9,2,0}},
            {4,{3,16,15,4},{10,3,0,1}},
{4,{13,14,15,16},{4,0,2,11}},
{4,{11,14,13,10},{0,3,12,5}},
{4,{11,10,9,12},{4,13,6,0}},
{4,{12,9,7,8},{5,14,7,0}},
{4,{7,6,5,8},{15,8,0,6}},
	  {4,{6,2,1,5},{16,1,0,7}},
{3,{17,3,2},{10,1,16}},
{3,{17,16,3},{11,2,9}},
{3,{17,13,16},{12,3,10}},
{3,{17,10,13},{13,4,11}},
{3,{17,9,10},{14,5,12}},
{3,{17,7,9},{15,6,13}},
{3,{17,6,7},{16,7,14}},
{3,{17,2,6},{9,8,15}}
};


(* faces with cone at each end *)
(*
faces={{4,{1,2,3,4},{8,9,2,17}},
            {4,{3,16,15,4},{10,3,18,1}},
{4,{13,14,15,16},{4,19,2,11}},
{4,{11,14,13,10},{20,3,12,5}},
{4,{11,10,9,12},{4,13,6,21}},
{4,{12,9,7,8},{5,14,7,22}},
{4,{7,6,5,8},{15,8,23,6}},
	  {4,{6,2,1,5},{16,1,24,7}},
{3,{17,3,2},{10,1,16}},
{3,{17,16,3},{11,2,9}},
{3,{17,13,16},{12,3,10}},
{3,{17,10,13},{13,4,11}},
{3,{17,9,10},{14,5,12}},
{3,{17,7,9},{15,6,13}},
{3,{17,6,7},{16,7,14}},
{3,{17,2,6},{9,8,15}},
{3, {1,4,18},{1,18,24}},
{3, {4,15,18},{2,19,17}},
{3, {15,14,18},{3,20,18}},
{3, {14,11,18},{4,21,19}},
{3, {11,12,18},{5,22,20}},
{3, {12,8,18},{6,23,21}},
{3, {8,5,18},{7,24,22}},
{3, {5,1,18},{8,17,23}}
};
*)


nfaces = Length[faces];

)];





(* ----------- Print Input Data-------------- *)Print["Input data points  ndata =",ndata,"  data = ",MatrixForm[data]];
Print["Input normals nn = ",MatrixForm[nn]];
Print["Input faces  nfaces =",nfaces];
Do[Print["face ",i," #verts =",faces[[i,1]],
                                          "   vertices =",faces[[i,2]],
                                          "   neighbors=",faces[[i,3]]],
{i,1,nfaces}];
Print["Created working Bezier patch bb1 = ",MatrixForm[bb1]];

(* ------------------------------------------------ *)






(* ::Section:: *)
(**)
(*Core G1 Functions*)


(* ::Input:: *)
project[aa_,nn_,xx_]:=xx+((aa-xx).nn)*nn;



(* project xx into plane through aa with unit normal nn *)
(* This ratio (1-t):t should depend on the normal angle ?? *)
(* Sphere can be good with t=0.5, but this causes problems in planar *)
(* t=0.4 also good for curved *)
(* Piper = t=1/3 *)

aux[i_,j_]:= (2/3)*data[[i]]+(1/3)*data[[j]];



(*boundary[i_,j_]:={data[[i]],project[data[[i]],nn[[i]],aux[i,j]],project[data[[j]],nn[[j]],aux[j,i]],data[[j]]};
*)
boundary[i_,j_]:={data[[i]], project[data[[i]],nn[[i]],aux[i,j]], project[data[[j]],nn[[j]],aux[j,i]], data[[j]]};




columnfix[a_,j_,vec_] :=(auxx=a;
Do[auxx[[i,j]]=vec[[i]],{i,1,Length[vec]}];
auxx
);

guesstri[pp_,qq_]:=(
(* "intelligent" guesses for the unknowns pp[[2]], pp[[3]]*)
(* shrink to deg 4 patch *)
ribbon1 = 3/4(pp[[1]] - qq[[1]]);
ribbon2 = 3/4(pp[[4]] - qq[[4]]);
{pp[[1]],
 qq[[2]] + (2ribbon1+ ribbon2 )/3,
 qq[[3]] + (ribbon1   + 2ribbon2)/3,
pp[[4]]}
)//N;

guessleft[pp_,qq_]:={
(* "intelligent" guesses for the unknowns pp[[2]], pp[[3]]*)
pp[[1]],
 qq[[2]] + (2(pp[[1]]-qq[[1]])+(pp[[4]]-qq[[4]]) )/3,
 qq[[3]] + ((pp[[1]]-qq[[1]])   + 2(pp[[4]]-qq[[4]]))/3,
pp[[4]]
}//N;

guessright[qq_,rr_]:=
{
(* "intelligent" guesses for the unknowns rr[[2]], rr[[3]]*)
rr[[1]],
 qq[[2]] + (2(rr[[1]]-qq[[1]])+(rr[[4]]-qq[[4]]) )/3,
 qq[[3]] + ((rr[[1]]-qq[[1]])   + 2(rr[[4]]-qq[[4]]))/3,
rr[[4]]
}//N;

elevate[bb_]:={bb[[1]],(bb[[1]]+3bb[[2]])/4,(bb[[2]]+bb[[3]])/2, (3bb[[3]]+bb[[4]])/4,
bb[[4]]}//N;



bern[n_,i_,t_]:=(
If[i==0, ti = 1.0, ti = t^i];
If[(n-i)== 0, nmit = 1.0, nmit = (1-t)^(n-i)];
Binomial[n,i]*nmit*ti
);

tribern[n_,i_,j_,u_,v_]:=(

facn = Factorial[n];
faci = Factorial[i];
facj = Factorial[j];
facnij = Factorial[n-i-j];

If[i == 0, ui = 1.0, ui = u^i];
If[j==0, vj = 1.0, vj = v^j];
If[(n-i-j) == 0, uvij = 1.0, uvij = (1.0-u-v)^(n-i-j)];

val = (facn/(faci*facj*facnij))*(ui*vj*uvij)
);

(*Gregory quad: *)
(* old gquad[bb_,c22_,c23_,c32_,c33_,u_,v_] *)
gquad[u_,v_]:=
(*the second and third rows of bb are from ribbons across v=0 and v=1 *)
(
aa=bb1;

aa[[2,2]]=(u*aa[[2,2]]+v*cc22)/(u+v);
aa[[2,3]]=(u*aa[[2,3]]+(1-v)*cc23)/(1+u-v);
aa[[3,2]]=((1-u)*aa[[3,2]]+v*cc32)/(1-u+v);
aa[[3,3]]=((1-u)*aa[[3,3]]+(1-v)*cc33)/(2-u-v);


Sum[aa[[i+1,j+1]]*bern[3,i,u]*bern[3,j,v],{i,0,3},{j,0,3}]
);

(*Gregory triangle: *)

gtri[u_,v_]:=(
ss=tt1;

eps = 0.00000001;
upv = u+v;
u1 = 1.0 -u;
v1 = 1.0 - v;
(*
If[Abs[upv] \[LessEqual] eps, If[upv< 0, upv = upv - eps, upv = upv + eps]];
If[Abs[u1]<= eps, If[u1 < 0, u1 =  u1 - eps, u1 = u1 + eps]];
If[Abs[v1] <= eps, If[v1 < 0, v1 =  v1 - eps, v1 = v1 + eps]];
*)

(*  This is the old one *)
(*
ss[[2,2]]=(u*ss[[2,2]]+v*dd22)/(u+v);
ss[[2,3]]=(u*ss[[2,3]]+(1-v)*dd23)/(1-v+u);
ss[[3,2]]=((1-u)*ss[[3,2]]+v*dd32)/(1-u+v);
*)

ss[[2,2]]=(u*ss[[2,2]]+v*dd22)/upv;
ss[[2,3]]=(u*ss[[2,3]]+(1-u-v)*dd23)/v1;
ss[[3,2]]=((1-u-v)*ss[[3,2]]+ v*dd32)/u1;

Sum[ss[[i+1,j+1]]*tribern[4,i,j,u,v],{j,0,4},{i,0,4-j}]
);

(* ************************************************* *)
(* Create G1 tangent ribbons. Numbering of topleft etc is *away* from qq. *)
interiorG1[qq_,topleft_,topright_,bottomleft_,bottomright_, nume_, numenbr_]:=
(
If[debug,
Print["in interiorG1"];
Print["qq = ",MatrixForm[qq],"  topl=",MatrixForm[topleft],"  topr = ",MatrixForm[topright],"botl=",MatrixForm[bottomleft],"  botr=",MatrixForm[bottomright]];
Print["qq = ",Reverse[qq],"\n  topl=",topleft,"\n  topr = ",topright,"\n  botl=",bottomleft,"\n  botr=",bottomright];
Print["numedges=",nume,"  numedges nbr = ",numenbr];
  ];

(* First, determine geometry parameters l0, m0  and l1 m1 (lambda's and mu's in text). *)

(* Create degree 4 boundary data that is needed *)
qqhat = elevate[qq];
If[debug, Print["After elevate: qqhat=",qqhat]; ];

If[nume == 4,
(bndrytl = topleft[[2]];
bndrybl = bottomleft[[2]]),
(bndrytl = (3*topleft[[2]] + topleft[[1]])/4.0;
bndrybl = (3*bottomleft[[2]] + bottomleft[[1]])/4.0)
];
If[numenbr == 4,
(bndrytr = topright[[2]];
bndrybr = bottomright[[2]]),
(bndrytr = (3*topright[[2]] + topright[[1]])/4.0;
bndrybr = (3*bottomright[[2]] + bottomright[[1]])/4.0)
];
If[debug,
Print["bndry points to use below depend on num edges"];
Print["bndrytl=",bndrytl,"bndrybl=",bndrybl," bndrytr=",bndrytr," bndrybr=",bndrybr];
  ];
(* Is it better to use the degree 4 boundary for the triangle in guess ?*)

pp=guessleft[{bndrytl,{0,0,0},{0,0,0},bndrybl},qq];
rr=guessright[qq,{bndrytr,{0,0,0},{0,0,0},bndrybr}];

If[debug,Print["Guess pp = ",MatrixForm[pp]," Guess rr ",MatrixForm[rr]]; ];


lmat = {-bndrytl + bndrytr, qqhat[[1]] - qqhat[[2]]};
ltrans = Transpose[lmat];
temp0= LeastSquares[ltrans, -bndrytl + qqhat[[1]] ];
l0=temp0[[1]];
m0 = temp0[[2]];

mmat = {-bndrybl + bndrybr, qqhat[[4]] - qqhat[[5]]};
mtrans = Transpose[mmat];
temp1 = LeastSquares[mtrans,-bndrybl + qqhat[[4]]];
l1=temp1[[1]];
m1 = temp1[[2]];
If[debug,
Print["Factors from degree 4 bndry: l0=",l0," m0=", m0,"  l1=", l1,"  m1=",m1, " m0= " , Transpose[lmat], " R0= " ,
-bndrytl + qqhat[[1]], " bndrytl = ", bndrytl, " qqhat[[1]]= ", qqhat[[1]], " m1= " , Transpose[mmat], " R1= " , -bndrybl + qqhat[[4]],
" bndrybl= " , bndrybl, " qqhat[4]= " , qqhat[[4]]
     ];
  ];

(* Having found the geometry params, compute the unknowns: *)
   mat = {{-3*(1 - l0),  -3*l0, 0, 0}, {(1 - l1), l1, (1 - l0), l0}, {0, 0, 3*(1 - l1), 3*l1}};

 lc[k_, i_] := (1 - k)*qqhat[[i]]  +  k*qqhat[[i + 1]] ;

rhs = {(1 - l1)*bndrytl + l1*bndrytr - lc[m1, 1] - 3*lc[m0, 2],
        lc[m0, 3] + lc[m1, 2],
        3*lc[m1, 3] - ((1 - l0)*bndrybl + l0*bndrybr - lc[m0, 4])};
If[debug,
Print["Matrix = ",MatrixForm[mat]];
Print["rhs = ",MatrixForm[rhs]];
Print["Matrix = ",Identity[mat]];
Print["rhs = ",Identity[rhs]];
  ];
approx = {pp[[2]], rr[[2]], pp[[3]], rr[[3]] };

yy = LeastSquares[mat.Transpose[mat], rhs - mat.approx];
solution = approx + Transpose[mat].yy;

(* For the tri patch, the first and last entries should be deg elevated *)
(* But these aren't used *)
outpp = {pp[[1]], solution[[1]], solution[[3]], pp[[4]]};outrr ={rr[[1]], solution[[2]], solution[[4]], rr[[4]]};

(*
Print["Leaving inGbot1 pp = ",MatrixForm[outpp]];
Print["Leaving inGbot1 rr = ",MatrixForm[outrr]];
*)
If[nume == 3,
If[debug, Print["Leaving inGbot1 \npp = ",outpp, "\nqq= ", qq, "\nrr= ", outrr]; ];
If[debug, Print["Leaving inGbot1 rev \npp = ",Reverse[outpp], "\nqq= ", Reverse[qq], "\nrr= ", Reverse[outrr]]; ];
  ];

outpp
);

(* ****************************************************************** *)
(* ****************************************************************** *)
(* ****************************************************************** *)

g1check[pp_,qq_,rr_]:=Max[
Table[
Abs[
Det[{
(1-t)^3*rr[[1]]
+ 3(1-t)^2*t *rr[[2]] + 3 (1-t)*t^2*rr[[3]] + t^3*rr[[4]]
  - ((1-t)^2*qq[[1]] + 2 (1-t)*t*qq[[2]] + t^2*qq[[3]]),
  (1-t)^2*qq[[2]] + 2 (1-t)*t*qq[[3]] + t^2*qq[[4]] - ((1-t)^2*qq[[1]] + 2 (1-t)*t*qq[[2]] + t^2*qq[[3]]),
(1-t)^3*pp[[1]]
+ 3(1-t)^2*t *pp[[2]] + 3 (1-t)*t^2*pp[[3]] + t^3*pp[[4]]-
  ((1-t)^2*qq[[1]] + 2 (1-t)*t*qq[[2]] + t^2*qq[[3]])
}]
],{t,0,1,.1}]
];



(* ::Section:: *)
(**)
(*Data Structure Related Functions*)
(**)


(* ::Input:: *)

(* ===============Rectangle patch: Load Default bb1 and cc data code===========================*)

(*
Input:global data,global normals,and faces info, global current iface
Output:global bb1,cc
Module
*)
loadbb1patch[iface_]:=(

(* Print["loadbb1patch: iface = ",iface]; *)

(* Load the boundaries *)

temp = boundary[faces[[iface]][[2]][[1]], faces[[iface]][[2]][[2]]];
Do[bb1[[i,1]] = temp[[i]],{i,1,4}];

temp = boundary[faces[[iface]][[2]][[4]], faces[[iface]][[2]][[3]]];
Do[bb1[[i,4]] = temp[[i]],{i,1,4}];

temp = boundary[faces[[iface]][[2]][[1]], faces[[iface]][[2]][[4]]];
Do[bb1[[1,i]] = temp[[i]],{i,1,4}];

temp = boundary[faces[[iface]][[2]][[2]], faces[[iface]][[2]][[3]]];
Do[bb1[[4,i]] = temp[[i]],{i,1,4}];


(* Load interior data *)

interiorrow = {bb1[[1,2]],{0,0,0},{0,0,0},bb1[[4,2]]};
bndry = {bb1[[1,1]], bb1[[2,1]],bb1[[3,1]],bb1[[4,1]]};
newrow = guessleft[interiorrow, bndry];
bb1[[2,2]] = newrow[[2]];
bb1[[3,2]] = newrow[[3]];

interiorrow = {bb1[[1,3]],{0,0,0},{0,0,0},bb1[[4,3]]};
bndry = {bb1[[1,4]], bb1[[2,4]],bb1[[3,4]],bb1[[4,4]]};
newrow = guessleft[interiorrow, bndry];
bb1[[2,3]] = newrow[[2]];
bb1[[3,3]] = newrow[[3]];

interiorrow = {bb1[[2,1]],{0,0,0},{0,0,0},bb1[[2,4]]};
bndry = {bb1[[1,1]], bb1[[1,2]],bb1[[1,3]],bb1[[1,4]]};
newrow = guessleft[interiorrow, bndry];
cc22 = newrow[[2]];
cc23 = newrow[[3]];

interiorrow = {bb1[[3,1]],{0,0,0},{0,0,0},bb1[[3,4]]};
bndry = {bb1[[4,1]], bb1[[4,2]],bb1[[4,3]],bb1[[4,4]]};
newrow = guessleft[interiorrow, bndry];
cc32 = newrow[[2]];
cc33 = newrow[[3]];

If[debug,
Print["loadbb1patch: loaded boundaries and interior guesses = ",MatrixForm[bb1]];
Print["  cc22=",cc22," cc23=",cc23,"  cc32=",cc32,"  cc33=",cc33];
  ];
)



(* ****************************************************************** *)
(* ****************************************************************** *)
(* ****************************************************************** *)


(* Triangle patch: Load Default bb1, tt1 and dd data code *)

(*
Input:global data,global normals,and faces info, global current iface
Output:global bb1,tt1,dd
Module
*)
loadpatchtricase[iface_]:=(

(* Print["loadpatchtricase: iface = ",iface]; *)

(* Load the boundaries *)

tedge1 = boundary[faces[[iface]][[2]][[1]], faces[[iface]][[2]][[2]]];

tedge2 = boundary[faces[[iface]][[2]][[2]], faces[[iface]][[2]][[3]]];

tedge3 = boundary[faces[[iface]][[2]][[3]], faces[[iface]][[2]][[1]]];



(* Load interior data *)

interiorrow = {tedge3[[3]],{0,0,0},{0,0,0},tedge2[[2]]};
bndry = tedge1;
newrow = guesstri[interiorrow, bndry];
tt1[[2,2]] = newrow[[2]];
tt1[[3,2]] = newrow[[3]];


interiorrow = {tedge2[[3]],{0,0,0},{0,0,0},tedge1[[2]]};
bndry = tedge3;
newrow = guesstri[interiorrow, bndry];
dd23 = newrow[[2]];
dd22 = newrow[[3]];


interiorrow = {tedge1[[3]],{0,0,0},{0,0,0},tedge3[[2]]};
bndry = tedge2;
newrow = guesstri[interiorrow, bndry];
dd32 = newrow[[2]];
tt1[[2,3]] = newrow[[3]];

deg4curve = elevate[tedge1];
Do[tt1[[i,1]] = deg4curve[[i]],{i,1,5}];


deg4curve = elevate[tedge2];
tt1[[4,2]] = deg4curve[[2]];
tt1[[3,3]] = deg4curve[[3]];
tt1[[2,4]] = deg4curve[[4]];


deg4curve = elevate[tedge3];
Do[tt1[[1, 6-i]] = deg4curve[[i]],{i,1,4}];

)


(* ::Input:: *)

(* ======================================= *)
(* ****************************************************************** *)
(* ****************************************************************** *)
(* ****************************************************************** *)

(*-----------------Rectangular Patch Rules--------------------*)(*
Rectangular patch should be in ccw order. Vertices v1,v2,v3,v4

v1 to v2-- called E1-- has Bezier points 11,12,13,14
v1 to v4-- called E4-- has Bezier points 11,21,31,41
Gregory points definition/convention:
Working across E1 and E3 load Gregory points into cc points.
Working across E2 and E4 load Gregory points into bb1 positions
*)


(* ===============Load G1 Input code===========================*)
(*
Input:global faces,inbr,iedge,bb1,nbrverts
Output:q,topl,topr,botl, botr
inbr is the patch neighboring the current patch's iedge bb1 should have correct boundary information loaded nbrverts take values 1 to 4:

nbrvert[[2]] has the face index for the neighboring patch that is at the start of q
nbrvert[[3]] is at the end of q
nbrvert[[1]] is nbrvert[[2]]-1 mod 4
nbrvert[[4]] is nbrvert[[3]]+1 mod 4
*)



(* ::Input:: *)
loadG1input[iface_, iedge_, nbr_]:= (

(* Print["loadG1input: iface=",iface,"  iedge=",iedge,"  nbr=",nbr]; *)

If[debug, Print["loadG1input 1"];];

If[iedge==1,(
	q={bb1[[4,1]],bb1[[3,1]],bb1[[2,1]],bb1[[1,1]]};
	topl={bb1[[4,1]],bb1[[4,2]],bb1[[4,3]],bb1[[4,4]]};
	botl={bb1[[1,1]],bb1[[1,2]],bb1[[1,3]],bb1[[1,4]]};)
];
If[debug, Print["loadG1input 2"];];

If[iedge==2, (
	q={bb1[[4,4]],bb1[[4,3]],bb1[[4,2]],bb1[[4,1]]};
	topl={bb1[[4,4]],bb1[[3,4]],bb1[[2,4]],bb1[[1,4]]};
	botl={bb1[[4,1]],bb1[[3,1]],bb1[[2,1]],bb1[[1,1]]};)
];
If[debug, Print["loadG1input 3"];];

If[iedge==3, (
	q={bb1[[1,4]],bb1[[2,4]],bb1[[3,4]],bb1[[4,4]]};
	topl={bb1[[1,4]],bb1[[1,3]],bb1[[1,2]],bb1[[1,1]]};
	botl={bb1[[4,4]],bb1[[4,3]],bb1[[4,2]],bb1[[4,1]]}; )
];
If[debug, Print["loadG1input 4"];];

If[iedge==4, (
	q={bb1[[1,1]],bb1[[1,2]],bb1[[1,3]],bb1[[1,4]]};
	topl={bb1[[1,1]],bb1[[2,1]],bb1[[3,1]],bb1[[4,1]]};
	botl={bb1[[1,4]],bb1[[2,4]],bb1[[3,4]],bb1[[4,4]]};)
];


If[faces[[nbr,1]] == 4,
nbrverts = loadnghbrinfo[iface, nbr],
nbrverts = loadnghbrinfotri[iface, nbr]
];
If[debug, Print["loadG1input 5"];];

(* Print["loadG1input: nbrverts = ",nbrverts]; *)

(*load data starting at common boundary*)

topr=boundary[
faces[[nbr,2,nbrverts[[2]]]],
faces[[nbr,2,nbrverts[[1]]]]];

botr=boundary[
faces[[nbr,2,nbrverts[[3]]]],
faces[[nbr,2,nbrverts[[4]]]]];

(*
Print["loadG1input: data loaded"];
Print["    q = ",q];
Print["    topl = ",topl];
Print["    botl = ",botl];
Print["    topr = ",topr];
Print["    botr = ",botr];
*)
)



loadG1inputtri[iface_, iedge_, nbr_]:= (

If[debug, Print["loadG1input tri case: iface=",iface,"  iedge=",iedge,"  nbr=",nbr]; ];

If[iedge==1,(
	q={tedge1[[4]], tedge1[[3]], tedge1[[2]],tedge1[[1]]};
	topl= tedge2;
	botl={tedge3[[4]], tedge3[[3]], tedge3[[2]],tedge3[[1]]};)
];

If[iedge==2, (
	q={tedge2[[4]], tedge2[[3]], tedge2[[2]],tedge2[[1]]};
	topl=tedge3;
	botl={tedge1[[4]], tedge1[[3]], tedge1[[2]],tedge1[[1]]};)
];

If[iedge==3, (
	q={tedge3[[4]], tedge3[[3]], tedge3[[2]],tedge3[[1]]};
	topl=tedge1;
	botl={tedge2[[4]], tedge2[[3]], tedge2[[2]],tedge2[[1]]}; )
];


If[faces[[nbr,1]] == 4,
nbrverts = loadnghbrinfo[iface, nbr],
nbrverts = loadnghbrinfotri[iface, nbr]
];



(*load data starting at common boundary*)

topr=boundary[
faces[[nbr,2,nbrverts[[2]]]],
faces[[nbr,2,nbrverts[[1]]]]];

botr=boundary[
faces[[nbr,2,nbrverts[[3]]]],
faces[[nbr,2,nbrverts[[4]]]]];

 If[debug,
Print["loadG1inputtri: data loaded"];
Print["    q = ",q];
Print["    topl = ",topl];
Print["    botl = ",botl];
Print["    topr = ",topr];
Print["    botr = ",botr];
   ];
)


(* ::Input:: *)

(* ===============Load Gregory Points code===========================*)
(*
Input:iedge,p
Output:bb1,cc

Load Gregory points into bb1 or cc arrays depending on edge working on.
Input array p is the row/column of Bezier points that is adjacent to boundary.
Don't need to reload the control points on the boundaries-- they didn't get modified.
Modify interior control points in the following way:
Edges 1 and 3 load into cc and Edges 2 and 4 load into bb1.

FIX FOR TRIANGLES!!!!
*)


(* ::Input:: *)
loadG1output[iedge_,p_] :=(

If[iedge==1,(bb1[[3,2]]=p[[2]]; bb1[[2,2]]=p[[3]])];

If[iedge==2,(cc33=p[[2]];cc32=p[[3]])];

If[iedge==3,(bb1[[2,3]]=p[[2]];bb1[[3,3]]=p[[3]])];

     If[iedge==4,(cc22=p[[2]];cc23=p[[3]])];

(*
Print["loadG1output: loaded p using case iedge=",iedge];
Print["loadG1output: bb1 = ",MatrixForm[bb1]];
Print["  cc22=",cc22," cc23=",cc23,"  cc32=",cc32,"  cc33=",cc33];
*)
)

loadG1outputtri[iedge_,p_] :=(

If[iedge==1,(tt1[[3,2]]=p[[2]]; tt1[[2,2]]=p[[3]])];

If[iedge==2,(tt1[[2,3]]=p[[2]];dd32=p[[3]])];

If[iedge==3,(dd22=p[[2]];dd23=p[[3]])];


(*
Print["loadG1outputtri: loaded p using case iedge=",iedge];
Print["loadG1outputtri: tt1 = ",MatrixForm[tt1]];
Print["  dd22=",dd22," dd23=",dd23,"  dd32=",dd32,"  dd33=",dd33];
*)
)


(* ::Input:: *)

(* ===============Construct neighbor vertices===========================*)
(*
Input:gobal faces,iface,inbr
Output:nbrverts which holds 4 vertex pointers This routine assumes that the patches are oriented ccw.
nbrverts = {l-1, l, l+1, l+2} where l is the vertex of the neighboring patch at the start of q
l is [1,4] -- a pointer into the faces list
*)

loadnghbrinfo[iface_, nbr_]:=(

nbrverts={0,0,0,0};

(* Print["loadnghbrinfo: neighbors of neighbor -- faces[nbr,3] = ",faces[[nbr,3]]]; *)

Do[If[faces[[nbr,3,ii]] == iface, nbrverts[[2]]=ii],{ii,1,4}];

If[nbrverts[[2]]==0, Print["ERROR: loadnghbrinfo nbrverts[[2]] = 0"]];

nbrverts[[1]]=nbrverts[[2]]-1;
If[nbrverts[[1]]== 0,nbrverts[[1]]=4];

nbrverts[[3]]=nbrverts[[2]]+1;
    If[nbrverts[[3]]==5,nbrverts[[3]]=1];

nbrverts[[4]]=nbrverts[[3]]+1;
    If[nbrverts[[4]]==5,nbrverts[[4]]=1];

(* Print["loadnbrinfo: {lm1, l, l+1, l+2} == nbrverts = ",nbrverts]; *)

Return[nbrverts]
)




(* Construct neighbor vertices for a neighbor that is a tri patch *)
(*
Input:gobal faces,iface,inbr
Output:nbrverts which holds 4 vertex pointers This routine assumes that the patches are oriented ccw.
nbrverts = {l-1, l, l+1, l+2} where l is the vertex of the neighboring patch at the start of q
l is [1,4] -- a pointer into the faces list

---- fix comment ... l+2 is same as l-1
*)

loadnghbrinfotri[iface_, nbr_]:=(

nbrverts={0,0,0,0};

Do[If[faces[[nbr,3,ii]] == iface, nbrverts[[2]]=ii],{ii,1,3}];

If[nbrverts[[2]]==0, Print["ERROR: loadnghbrinfo nbrverts[[2]] = 0"]];

nbrverts[[1]]=nbrverts[[2]]-1;
If[nbrverts[[1]]== 0,nbrverts[[1]]=3];

nbrverts[[3]]=nbrverts[[2]]+1;
    If[nbrverts[[3]]==4,nbrverts[[3]]=1];

nbrverts[[4]]=nbrverts[[3]]+1;
    If[nbrverts[[4]]==4,nbrverts[[4]]=1];

Return[nbrverts];
)







maketrinet[iface_] := (

trinet  = {{tt1[[1,1]], tt1[[2,1]], tt1[[3,1]], tt1[[4,1]], tt1[[5,1]]}, {tt1[[1,2]], tt1[[2,2]], tt1[[3,2]], tt1[[4,2]]}, {tt1[[1,3]], tt1[[2,3]], tt1[[3,3]]},
{tt1[[1,4]], tt1[[2,4]]},
{tt1[[1,5]]}};

trinetu  = {{tt1[[1,1]], tt1[[2,1]], tt1[[3,1]], tt1[[4,1]], tt1[[5,1]]}, {tt1[[1,2]], tt1[[2,2]], tt1[[3,2]], tt1[[4,2]]}, {tt1[[1,3]], tt1[[2,3]], tt1[[3,3]]}, {tt1[[1,4]], tt1[[2,4]]}};

trinetv  = {{tt1[[1,1]], tt1[[1,2]], tt1[[1,3]], tt1[[1,4]], tt1[[1,5]]}, {tt1[[2,1]], tt1[[2,2]], tt1[[2,3]], tt1[[2,4]]}, {tt1[[3,1]], tt1[[3,2]], tt1[[3,3]]}, {tt1[[4,1]], tt1[[4,2]]}};

trinetw  = {{tt1[[5,1]], tt1[[4,2]], tt1[[3,3]], tt1[[2,4]], tt1[[1,5]]}, {tt1[[4,1]], tt1[[3,2]], tt1[[2,3]], tt1[[1,4]]}, {tt1[[3,1]], tt1[[2,2]], tt1[[1,3]]}, {tt1[[2,1]], tt1[[1,2]]}};

bndrytrisv0 = {tt1[[2,1]], tt1[[2,2]], tt1[[3,1]], tt1[[3,2]],tt1[[4,1]]};
bndrytrisu0 = {tt1[[1,2]], dd22, tt1[[1,3]], dd23, tt1[[1,4]]};
bndrytrisw0 = {tt1[[4,2]], dd32, tt1[[3,3]], tt1[[2,3]], tt1[[2,4]]};


)


(* ::Input:: *)

(* Compute the partials du and dv of the rectangular patch at the edge specified *)
(* iedge = 1 (u,v=0)   edge=2 (u=1, v) *)

secantpartialsrec[iedge_, param_] := (

eps = 0.000000005;

If[iedge == 1,(
partialdu = (gquad[param+eps, 0.0] - gquad[param, 0.0])/eps;
partialdv = (gquad[param, eps] - gquad[param, 0.0])/eps)
];

If[iedge == 2, (
partialdu = (gquad[1.0, param] - gquad[1.0 - eps, param])/eps;
partialdv = (gquad[1.0, param + eps] - gquad[1.0, param])/eps)
];

If[iedge == 3, (
partialdu = (gquad[param + eps, 1.0] - gquad[param, 1.0])/eps;
partialdv = (gquad[param, 1.0-eps] - gquad[param, 1.0])/eps)
];

If[iedge == 4, (
partialdu = (gquad[eps, param] - gquad[0.0, param])/eps;
partialdv = (gquad[0.0, param + eps] - gquad[0.0, param])/eps)
];

Print["secant partial rectangle: iedge=",iedge,"  du = ",partialdu,"  dv = ",partialdv];

);



(* ::Input:: *)

(* Compute the partials du and dv of the triangular patch at the edge specified *)
(* (u,v) origin at v1. e1 corresponds to u-dir. for e1 and e3 du and dv are consistent with param space. For e2 boundary curve deriv is stored in du and cross bndry is in dv  *)

secantpartialstri[iedge_, param_] := (

eps = 1.0*10^(-10);
eps2 = eps/2;

If[iedge == 1,(
partialdu = (gtri[param+eps, 0.0] - gtri[param, 0.0])/eps;
partialdv = (gtri[param-eps2, eps] - gtri[param, 0.0])/eps)
];

If[iedge == 2, (
partialdu = (gtri[param-eps, param+eps] - gtri[param, param])/eps;
partialdv = (gtri[param-eps2, param -eps2] - gtri[param, param])/eps)
];

If[iedge == 3, (
partialdu = (gtri[eps, param-eps2] - gtri[0.0, param])/eps;
partialdv = (gtri[0.0, param + eps] - gtri[0.0, param])/eps)
];
(*
Print["secant partial triangle: iedge=",iedge,"  du = ",partialdu,"  dv = ",partialdv];
*)
);



(* ::Section:: *)
(**)
(*Color Functions*)


(* ::Input:: *)


zebracolorfct[val_] := (
onstripe = 0;
vizeps = 0.05;
isoval = 0.4;
If[val >= (isoval-vizeps) && val <= (isoval + vizeps),onstripe = 1];
isoval = 0.8;
If[val >= (isoval-vizeps) && val <= (isoval + vizeps),onstripe = 1];
isoval = 1.2;
If[val >= (isoval-vizeps) && val <= (isoval + vizeps),onstripe = 1];
If[onstripe == 1, RGBColor[0,0,0], RGBColor[1,1,1]]
);





(* ::Section:: *)
(*Gregory coefficients computed + graphics*)


(* ::Input:: *)

Print["Run through the input to create a input data plot"];



(* u- and v-partials for four edges of each face. If triangle face, don't use fourth set *)
If[checkg1condition == 1,(
partials = Table[{{{0,0,0},{0,0,0}},
                                   {{0,0,0},{0,0,0}},
                                  {{0,0,0},{0,0,0}},
                                  {{0,0,0},{0,0,0}}}, {i,1,nfaces}];
)];



(* Generate a plot of the input data *)

plotcmatinput = Table[0,{i,1,nfaces}];
plotbb1input= Table[0,{i,1,nfaces}];
quadpicinput= Table[0,{i,1,nfaces}];
plotdata = Table[0,{i,1,nfaces}];
plotfaces = Table[0,{i,1,nfaces}];
plotnormals = Table[0,{i,1,nfaces}];
plotbndrypoly = Table[0,{i,1,nfaces}];
plotbndrycurve = Table[0,{i,1,nfaces}];


(* Process each face to create a display list *)
Do[(
(*Extract the number of edges*)
numedges = faces[[iface,1]];

(*Print["face ",iface,"  number of edges = ",faces[[iface,1]]]; *)
If[numedges < 3 || numedges > 4, Print["ERROR -- numfaces out of range -- numfaces=",numedges]];

(* RECTANGLE *)
If[numedges == 4,(

loadbb1patch[iface];
(*
Print["Initial iface=",iface,"   bb1 = ",MatrixForm[bb1]];
Print["  cc22=",cc22," cc23=",cc23,"  cc32=",cc32,"  cc33=",cc33];
*)
facedata ={data[[faces[[iface,2,1]]]], data[[faces[[iface,2,2]]]], data[[faces[[iface,2,3]]]], data[[faces[[iface,2,4]]]]};
plotdata[[iface]] =Graphics3D[{PointSize[pointsizelarge],colordatapoints,Point[facedata]}];
plotfaces[[iface]] = Graphics3D[{colorpolygon, Thickness[linewidththick], Line[facedata], Line[{facedata[[1]],facedata[[4]]}]}];

facenormals = {nn[[faces[[iface,2,1]]]], nn[[faces[[iface,2,2]]]], nn[[faces[[iface,2,3]]]], nn[[faces[[iface,2,4]]]]};
plotnormals[[iface]] = Graphics3D[{colornormals, Thickness[linewidththick], Arrow[{facedata[[1]], facedata[[1]] + normallengthfactor*facenormals[[1]]}],
Arrow[{facedata[[2]], facedata[[2]] + normallengthfactor*facenormals[[2]]}],
Arrow[{facedata[[3]], facedata[[3]] + normallengthfactor*facenormals[[3]]}],
Arrow[{facedata[[4]], facedata[[4]] + normallengthfactor*facenormals[[4]]}]}];

plotbb1input[[iface]]=Graphics3D[{colorrec1,PointSize[pointsizemedium],Map[Point,bb1],colorpolygon,Thickness[linewidththin],
Line[{bb1[[2,1]],bb1[[2,2]]}], Line[{bb1[[3,1]],bb1[[3,2]]}],
Line[{bb1[[2,4]],bb1[[2,3]]}], Line[{bb1[[3,4]],bb1[[3,3]]}],
Line[{bb1[[1,2]],cc22}], Line[{bb1[[1,3]],cc23}],
Line[{bb1[[4,2]],cc32}], Line[{bb1[[4,3]],cc33}]}];

cmat1={{cc22,cc23},{cc32,cc33}};
plotcmatinput[[iface]]=Graphics3D[{colorrec1,PointSize[pointsizemedium],Map[Point,cmat1]}];

(* boundary curves *)
bndry1 = {bb1[[1,1]],bb1[[2,1]], bb1[[3,1]],bb1[[4,1]]};
bndry2 = {bb1[[4,1]],bb1[[4,2]], bb1[[4,3]],bb1[[4,4]]};
bndry3 = {bb1[[1,4]],bb1[[2,4]], bb1[[3,4]],bb1[[4,4]]};
bndry4 = {bb1[[1,1]],bb1[[1,2]], bb1[[1,3]],bb1[[1,4]]};

plotbndrypoly[[iface]] = Graphics3D[{colorrec1,Thickness[linewidththick],Line[bndry1],
Line[bndry2],Line[bndry3],Line[bndry4],
colorrec1,PointSize[pointsizemedium], Point[{bb1[[1,1]],bb1[[2,1]], bb1[[3,1]],bb1[[4,1]],
bb1[[1,4]],bb1[[2,4]], bb1[[3,4]],bb1[[4,4]],bb1[[1,2]], bb1[[1,3]],bb1[[4,2]], bb1[[4,3]]}]}];

plotbndrycurve[[iface]] = Graphics3D[{Black, Thickness[linewidththick],BezierCurve[bndry1], BezierCurve[bndry2], BezierCurve[bndry3],BezierCurve[bndry4]}];

If[rendermethod == 0,
quadpicinput [[iface]]=
ParametricPlot3D[gquad[u,v],{v,0,1},{u,0,1},PlotStyle->{colorsurface,Specularity[White,myspecularity], Opacity[myopacity]},Mesh->None]
];
If[rendermethod == 1,
quadpicinput [[iface]]=
ParametricPlot3D[gquad[u,v],{v,0,1},{u,0,1},Mesh->False,PlotPoints->myplotpoints,ColorFunction->Function[{x,y,z,v,u},zebracolorfct[z]],ColorFunctionScaling->False]
];


(* plot arrows on first patch *)
If[iface == 1,(
axisc1input=Graphics3D[{Green,Arrow[{cc22,cc23}]}];
axisbb1input=Graphics3D[{Red,Arrow[{bb1[[1,1]],bb1[[1,2]]}]}];
)];

(* Get partials to check G1 condition *)
If[checkg1condition == 1,(
Do[(
nbr=faces[[iface,3,iedge]];
If[nbr != 0,(
(* store partials at this edge *)

If[iedge == 1, uu=0.5; vv=0.0];
If[iedge == 2, uu=1.0; vv=0.5];
If[iedge == 3, uu=0.5; vv=1.0];
If[iedge == 4, uu=0.0; vv=0.5];
 gotpartial= Derivative[1,0][gquad][uu,vv];
partials[[iface,iedge,1]] = (gotpartial/Norm[gotpartial]);
gotpartial = Derivative[0,1][gquad][uu,vv];
partials[[iface,iedge,2]] = (gotpartial/Norm[gotpartial]);
(*  Derivative command works well for rectangle patches
secantpartialsrec[iedge, 0.5];
partials[[iface,iedge,1]] = partialdu;
partials[[iface,iedge,2]] = partialdv;
*)

)]
),{iedge,1,4}]
)] (* if checkg1 *)
)]; (* If numedges = 4 *)


(* TRIANGLE *)
If[numedges == 3,(
loadpatchtricase[iface];

facedata ={data[[faces[[iface,2,1]]]], data[[faces[[iface,2,2]]]], data[[faces[[iface,2,3]]]]};
plotdata[[iface]] =Graphics3D[{colordatapoints,PointSize[pointsizelarge],Point[facedata]}];
plotfaces[[iface]] = Graphics3D[{colorpolygon, Thickness[linewidththick], Line[facedata], Line[{facedata[[1]],facedata[[3]]}]}];

facenormals = {nn[[faces[[iface,2,1]]]], nn[[faces[[iface,2,2]]]], nn[[faces[[iface,2,3]]]]};
plotnormals[[iface]] = Graphics3D[{colornormals, Thickness[linewidththick], Arrow[{facedata[[1]], facedata[[1]] + normallengthfactor*facenormals[[1]]}],
Arrow[{facedata[[2]], facedata[[2]] + normallengthfactor*facenormals[[2]]}],
Arrow[{facedata[[3]], facedata[[3]] + normallengthfactor*facenormals[[3]]}]}];

maketrinet[iface];
plotbb1input[[iface]]=Graphics3D[{colortri1,PointSize[pointsizemedium], Map[Point,trinet],colorpolygon, Thickness[linewidththin],Line[bndrytrisu0],Line[bndrytrisv0],Line[bndrytrisw0]}];


dmat1 = {dd22, dd32, dd23, dd22};
plotcmatinput[[iface]]=Graphics3D[{colortri1,PointSize[pointsizemedium],Map[Point,dmat1],colorpolygon,Thickness[linewidththin]}];


(* boundary curves *)
bndry1 = {tt1[[1,1]],tt1[[2,1]],tt1[[3,1]], tt1[[4,1]], tt1[[5,1]]};
bndry2 = {tt1[[5,1]], tt1[[4,2]], tt1[[3,3]], tt1[[2,4]], tt1[[1,5]]};
bndry3 = {tt1[[1,1]], tt1[[1,2]], tt1[[1,3]], tt1[[1,4]], tt1[[1,5]]};

plotbndrypoly[[iface]] = Graphics3D[{colortri1,Thickness[linewidththick],Line[{tt1[[1,1]],tt1[[2,1]],tt1[[3,1]], tt1[[4,1]], tt1[[5,1]]}],
Line[{tt1[[1,1]], tt1[[1,2]], tt1[[1,3]], tt1[[1,4]], tt1[[1,5]]}],
Line[{tt1[[5,1]], tt1[[4,2]], tt1[[3,3]], tt1[[2,4]], tt1[[1,5]]}],
colortri1,PointSize[pointsizemedium], Point[{tt1[[1,1]],tt1[[2,1]],tt1[[3,1]], tt1[[4,1]], tt1[[5,1]],tt1[[1,2]], tt1[[1,3]], tt1[[1,4]], tt1[[1,5]],tt1[[4,2]], tt1[[3,3]], tt1[[2,4]]}]}];

plotbndrycurve[[iface]] = Graphics3D[{Black, Thickness[linewidththick],BezierCurve[bndry1, SplineDegree->4], BezierCurve[bndry2, SplineDegree->4], BezierCurve[bndry3, SplineDegree->4]}];


(* Problems with this on w edge
RegionFunction\[Rule]Function[{x,y,z,u,v}, 0 <= u+v <= 1.0]
*)
If[rendermethod == 0,
quadpicinput [[iface]]=
ParametricPlot3D[gtri[u,v],{v,0.0, 1.0},{u,0.0,1.0-v},
Mesh->None,PlotStyle->{colorsurface,Specularity[White,myspecularity],Opacity[myopacity]}]
];


If[rendermethod == 1,
quadpicinput [[iface]]=
ParametricPlot3D[gtri[u,v],{v,0.0, 1.0},{u,0.0,1.0-v},Mesh->False,PlotPoints->myplotpoints,ColorFunction->Function[{x,y,z,v,u},zebracolorfct[z]],ColorFunctionScaling->False]
];

(* ColorFunction\[Rule]Function[{x,y,z,v,u},ColorData["DarkBands"][z/1.6]] *)


(* plot arrows on first patch *)
If[iface == 1,(
axisc1input=Graphics3D[{Green,Arrow[{dd22,dd23}]}];
axisbb1input=Graphics3D[{Red,Arrow[{tt1[[1,1]],tt1[[1,2]]}]}];
)]


(* Get partials to check G1 condition *)
If[checkg1condition == 1, (
Do[(
nbr=faces[[iface,3,iedge]];
If[nbr != 0,(
(* store partials at this edge *)
(*
If[iedge \[Equal] 1, uu=0.5; vv=0.0];
If[iedge \[Equal] 2, uu=0.5; vv=0.5];
If[iedge \[Equal] 3, uu=0.0; vv=0.5];
 gotpartial= Derivative[1,0][gtri][uu,vv];
partials[[iface,iedge,1]] = gotpartial;
gotpartial = Derivative[0,1][gtri][uu,vv];
partials[[iface,iedge,2]] = gotpartial;
Print["Derivative tri: iedge=",iedge,"  du=",partials[[iface,iedge,1]],"  dv=",partials[[iface,iedge,2]]];
*)

secantpartialstri[iedge, 0.5];
partials[[iface,iedge,1]] = (partialdu/Norm[partialdu]);
partials[[iface,iedge,2]] = partialdv/Norm[partialdv];
)];
),{iedge,1,3}]
)] (* if checkg1 *)
)] (* if triangle case *)

),{iface,1,nfaces}];



(* Check G1 conditions across faces *)
If[checkg1condition == 1, (
Print["Check G1 conditions **before ** G1 modifications "];
(*Print["Partials = ",partials]; *)

Do[(
numedges = faces[[iface,1]];
Do[(
nbr=faces[[iface,3,iedge]];
If[nbr > iface,(
numedgesnbr = faces[[nbr,1]];
If[numedgesnbr == 4,
loadnghbrinfo[iface, nbr],
loadnghbrinfotri[iface,nbr]
];
nbredge = nbrverts[[2]];

If[numedges == 4,(
If[iedge == 1 || iedge == 3,iacross=2, iacross=1];
volume = Det[{partials[[iface,iedge,iacross]],partials[[nbr,nbredge,1]], partials[[nbr,nbredge,2]]}]),

(If[iedge == 1 || iedge == 2,iacross=2, iacross=1];
volume = Det[{partials[[iface,iedge,iacross]],partials[[nbr,nbredge,1]], partials[[nbr,nbredge,2]]}];
)];
(*
Print["Main: Face #",iface,"  iedge=",iedge," u-partial=",partials[[iface]][[iedge]][[1]],"  v-partial=",partials[[iface]][[iedge]][[2]]];
Print["Main: Face #",nbr,"  nbredge=",nbredge," u-partial=",partials[[nbr]][[nbredge]][[1]],"  v-partial=",partials[[nbr]][[nbredge]][[2]]];
*)
Print["                   Across face ",iface," and face",nbr," is it G1 at mid-param of edge?  volume of partials = ",volume];

)] (* if nbrface index > current face index *)
),{iedge,1,numedges}]
),{iface,1,nfaces}]
)]; (* if checkg1 *)

(* ============================================================\[Equal] *)
(* ============================================================\[Equal] *)
(* ============================================================\[Equal] *)
(* ============================================================\[Equal] *)

Print["Run the the data to do G1 correction "];

(* Now do the correction and generate a plot *)

plotcmat1 = Table[0,{i,1,nfaces}];
plotbb1= Table[0,{i,1,nfaces}];
quadpic1= Table[0,{i,1,nfaces}];

Do[(


(*Extract the number of edges*)
numedges = faces[[iface,1]];

If[debug,Print["face ",iface,"  number of edges = ",faces[[iface,1]]]; ];
If[numedges < 3 || numedges > 4, Print["ERROR -- numfaces out of range -- numfaces=",numedges]];

(* RECTANGLE *)
If[numedges == 4,(

If[debug, Print["iface= " ,iface , " numedges= " , numedges];];

loadbb1patch[iface];
Do[(
If[debug, Print["Main: Face #",iface,"  iedge=",iedge]; ];
(*Extract the neighbor across edge*)
nbr=faces[[iface,3,iedge]];
If[debug, Print["Main: nbr = ",nbr]; ];

If[nbr != 0,(
If[debug,
Print["       Face #",iface,"  iedge=",iedge," nbr=",nbr,"  Computing G1 condition"];
  ];

numedgesnbr = faces[[nbr,1]];

If[debug,
Print["    loadG1input..."];
  ];

loadG1input[iface, iedge, nbr];

p=interiorG1[q,topl,topr,botl,botr, numedges, numedgesnbr] ;
loadG1output[iedge,p];

(* store partials at this edge *)
If[checkg1condition == 1, (
If[iedge == 1, uu=0.5; vv=0.0];
If[iedge == 2, uu=1.0; vv=0.5];
If[iedge == 3, uu=0.5; vv=1.0];
If[iedge == 4, uu=0.0; vv=0.5];
 gotpartial= Derivative[1,0][gquad][uu,vv];
partials[[iface,iedge,1]] = gotpartial/Norm[gotpartial];
gotpartial = Derivative[0,1][gquad][uu,vv];
partials[[iface,iedge,2]] = gotpartial/Norm[gotpartial];
(* for rectangles, Derivative seems to work fine *)
(*
secantpartialsrec[iedge, 0.5];
partials[[iface,iedge,1]] = partialdu;
partials[[iface,iedge,2]] = partialdv;
*)
)]
 If[debug,
Print["Main: Updated bb1 = ",MatrixForm[bb1]];
   Print["Main: check updates to cc  "];
  Print["  cc22=",cc22," cc23=",cc23,"  cc32=",cc32,"  cc33=",cc33];
   ];
)] (* nbr <> 0 *)
),{iedge,1,4}];

If[debug,
Print["Main: All done with patch   bb1 = ",MatrixForm[bb1]];
Print["  cc22=",cc22," cc23=",cc23,"  cc32=",cc32,"  cc33=",cc33];
  ];

(* plot arrows on first patch *)
If[iface == 1,(
axisc1=Graphics3D[{Green,Arrow[{cc22,cc23}]}];
axisbb1=Graphics3D[{Red,Arrow[{bb1[[1,1]],bb1[[1,2]]}]}];
)];

plotbb1[[iface]]=Graphics3D[{colorrec1,PointSize[pointsizemedium],Map[Point,bb1],colorpolygon, Thickness[linewidththin],Line[{bb1[[2,1]],bb1[[2,2]]}], Line[{bb1[[3,1]],bb1[[3,2]]}],
Line[{bb1[[2,4]],bb1[[2,3]]}], Line[{bb1[[3,4]],bb1[[3,3]]}],
Line[{bb1[[1,2]],cc22}], Line[{bb1[[1,3]],cc23}],
Line[{bb1[[4,2]],cc32}], Line[{bb1[[4,3]],cc33}]}];

cmat1={{cc22,cc23},{cc32,cc33}};

plotcmat1[[iface]]=Graphics3D[{colorrec1,PointSize[pointsizemedium],Map[Point,cmat1]}];



If[rendermethod == 0,
quadpic1 [[iface]]=
ParametricPlot3D[gquad[u,v],{v,0,1},{u,0,1},PlotStyle->{colorsurface,Specularity[White,myspecularity],Opacity[myopacity]},Mesh->None]
];
If[rendermethod == 1,
quadpic1 [[iface]]=
ParametricPlot3D[gquad[u,v],{v,0,1},{u,0,1},Mesh->False,PlotPoints->myplotpoints,ColorFunction->Function[{x,y,z,v,u},zebracolorfct[z]],ColorFunctionScaling->False]
]
(*ColorFunction\[Rule]Function[{x,y,z,v,u},ColorData["DarkBands"][z/1.6]],ColorFunctionScaling\[Rule]False] *)

(*quadrefcurve=ParametricPlot[quadrefline[l1_,l2_,t_,bb_,c22_,c23_,c32_,c33_],{t,0,1}];
*)
(quadData[iface] = Join[Flatten[bb1,1], {cc22, cc23, cc32, cc33}];)

)]; (* rectangle case *)



(* TRIANGLE *)
If[numedges == 3,(
loadpatchtricase[iface];
If[debug,
Print["Main: Initialized tri patch   tt1 = ",MatrixForm[tt1]];
Print["  dd22=",dd22," dd23=",dd23,"  dd32=",dd32];
  ];
Do[(
If[debug, Print["Main: Face #",iface,"  iedge=",iedge]; ];

(*Extract the neighbor across edge*)
nbr=faces[[iface,3,iedge]];
If[debug, Print["Main: nbr = ",nbr]; ];


If[nbr != 0,(
If[debug,
Print["       Face #",iface,"  iedge=",iedge," nbr=",nbr,"  Computing G1 condition"];
  ];

numedgesnbr = faces[[nbr,1]];

loadG1inputtri[iface, iedge, nbr];

p=interiorG1[q,topl,topr,botl,botr, numedges, numedgesnbr] ;

loadG1outputtri[iedge,p];

(* store partials at this edge *)
If[checkg1condition == 1, (
(*
If[iedge \[Equal] 1, uu=0.5; vv=0.0];
If[iedge \[Equal] 2, uu=0.5; vv=0.5];
If[iedge \[Equal] 3, uu=0.0; vv=0.5];

 gotpartial= Derivative[1,0][gtri][uu,vv];
partials[[iface,iedge,1]] = gotpartial;
gotpartial = Derivative[0,1][gtri][uu,vv];
partials[[iface,iedge,2]] = gotpartial;

Print["Derivative tri: iedge=",iedge,"  du=",partials[[iface,iedge,1]],"  dv=",partials[[iface,iedge,2]]];
*)
secantpartialstri[iedge, 0.5];
partials[[iface,iedge,1]] = partialdu/Norm[partialdu];
partials[[iface,iedge,2]] = partialdv/Norm[partialdv];

)]
)] (* nbr <> 0 *)
),{iedge,1,3}];

If[debug,
Print["Main: All done with patch   tt1 = ",MatrixForm[tt1]];
Print["  dd22=",dd22," dd23=",dd23,"  dd32=",dd32];
  ];


maketrinet[iface];
plotbb1[[iface]]=Graphics3D[{colortri1,PointSize[pointsizemedium],Map[Point,trinet],colorpolygon,Thickness[linewidththin], Line[bndrytrisu0], Line[bndrytrisv0], Line[bndrytrisw0]}];

dmat1 = {dd22, dd32, dd23, dd22};
plotcmat1[[iface]]=Graphics3D[{colortri1,PointSize[pointsizemedium],Map[Point,dmat1]}];


If[rendermethod ==0,
quadpic1 [[iface]]=
ParametricPlot3D[gtri[u,v],{v,0,1},{u,0,1-v},PlotStyle->{colorsurface,Specularity[White,myspecularity],Opacity[myopacity]},Mesh->None]
];
If[rendermethod == 1,
quadpic1 [[iface]]=
ParametricPlot3D[gtri[u,v],{v,0,1},{u,0,1-v},Mesh->False,PlotPoints->myplotpoints,ColorFunction->Function[{x,y,z,v,u},zebracolorfct[z]],ColorFunctionScaling->False]
];



(* plot arrows on first patch *)
If[iface == 1,(
axisc1=Graphics3D[{Green,Arrow[{dd22,dd23}]}];
axisbb1=Graphics3D[{Red,Arrow[{tt1[[1,1]],tt1[[1,2]]}]}];
)]

(
 triData[iface] = Join[Flatten[Table[Table[tt1[[i,j]],{j,6-i}],{i,5}],1], {dd22, dd23, dd32}];
)


)]; (* Triangle case *)

),{iface,1,nfaces}];


(* Check G1 conditions across faces *)
If[checkg1condition == 1, (
Print["Check partials *after* G1 conditions: "];
Do[(
numedges = faces[[iface,1]];
Do[(
nbr=faces[[iface,3,iedge]];
If[nbr > iface,(
numedgesnbr = faces[[nbr,1]];
If[numedgesnbr == 4,
loadnghbrinfo[iface, nbr],
loadnghbrinfotri[iface,nbr]
];
nbredge = nbrverts[[2]];

If[numedges == 4,(
If[iedge == 1 || iedge == 3,iacross=2, iacross=1];
volume = Det[{partials[[iface]][[iedge]][[iacross]],partials[[nbr]][[nbredge]][[1]], partials[[nbr]][[nbredge]][[2]]}]),

(If[iedge == 1 || iedge == 2,iacross=2, iacross=1];
volume = Det[{partials[[iface]][[iedge]][[iacross]],partials[[nbr]][[nbredge]][[1]], partials[[nbr]][[nbredge]][[2]]}];
)];
(*
Print["Main: Face #",iface,"  iedge=",iedge," u-partial=",partials[[iface]][[iedge]][[1]],"  v-partial=",partials[[iface]][[iedge]][[2]]];
Print["Main: Face #",nbr,"  nbredge=",nbredge," u-partial=",partials[[nbr]][[nbredge]][[1]],"  v-partial=",partials[[nbr]][[nbredge]][[2]]];
*)

Print["                   Across face ",iface," and face",nbr," is it G1 at mid-param of edge?  volume of partials = ",volume];
)];
),{iedge,1,numedges}];
),{iface,1,nfaces}];
)];




If[withcontrolstructure == 1,picc0 = Show[quadpicinput,  plotbb1input,plotcmatinput, plotbndrypoly, plotdata, Axes->False, Boxed -> False,ViewPoint-> myviewpoint,ViewVertical->myviewvertical,PlotRange->Automatic,PreserveImageOptions ->False]]
If[withcontrolstructure == 1,picg1 = Show[quadpic1, plotbb1,plotcmat1,plotbndrypoly,plotdata, Axes->False, Boxed -> False,ViewPoint-> myviewpoint,ViewVertical->myviewvertical,PlotRange->Automatic,PreserveImageOptions ->False]]

If[withcontrolstructure == 0,picc0 = Show[quadpicinput,Axes->False, Boxed -> False, ViewPoint-> myviewpoint,ViewVertical->myviewvertical,PlotRange->Automatic,PreserveImageOptions ->False]]

If[withcontrolstructure == 0,picg1 = Show[quadpic1,Axes->False, Boxed -> False, ViewPoint-> myviewpoint, ViewVertical->myviewvertical,PlotRange->Automatic,PreserveImageOptions ->False]]

picinput = Show[plotdata,plotfaces,plotnormals, Axes->False, Boxed -> False, ViewPoint-> myviewpoint,ViewVertical->myviewvertical,PlotRange->Automatic,PreserveImageOptions ->False]

picbndry = Show[plotbndrycurve,plotbndrypoly, plotdata,Axes->False, Boxed -> False, ViewPoint->myviewpoint,ViewVertical->myviewvertical,PlotRange->Automatic,PreserveImageOptions ->False]

(*
Export["picc0.eps", picc0];
Export["picg1.eps", picg1];
Export["picbndry.eps", picbndry];
Export["picinput.eps", picinput];
*)



