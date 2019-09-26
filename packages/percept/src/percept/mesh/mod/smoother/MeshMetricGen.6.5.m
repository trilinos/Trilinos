(* ::Package:: *)

(* ::Subsection:: ==================================================================================================================================================== *)
(*Code Gen Helpers*)

powerRules:={Power[x_,2]->MyPow2[x],Power[x_,3]->MyPow3[x],Power[x_,-1]->MyInverse[x]};

MyPow[x_,2]:=MyPow2[x];
MyPow[x_,3]:=MyPow3[x];
MyPow[x_,4]:=MyPow4[x];
MyPow[x_,0]:=1;
MyPow[x_,1]:=x;
MyPow[x_,(n_/;IntegerQ[n]&&n<0)]:=MyPow[MyInverse[x],-n];
powerRules:={Power[x_,n_]:>MyPow[x,n]};

MakeVec[p_,i_]:=Table[p[i,j],{j,0,2}];
MakeVec[p_]:=Table[p[j],{j,0,2}];

V[p_,i_]:=MakeVec[p,i];
V[p_]:=MakeVec[p];

Formatter:=CForm;
IndentString:="      ";
Clear[ConvertRule];
ConvertRule[assign_:" = ", ltype_:" double "]:=(lhs_->rhs_):>"\n" <> IndentString<>(If[MatchQ[lhs,x_[i__]]," ",ltype,""])<>ToString[Formatter[lhs]]<> assign <>ToString[Formatter[rhs]]<>" ;";

Clear[CodeGen];
CodeGenFile[listIn_,file_,assign_:" = ", ltype_ : " double "]:=Module[{list=listIn//.powerRules/.ConvertRule[assign, ltype]},Do[PutAppend[OutputForm[list[[i]]],file],{i,Length[list]}]];
CodeGen[listIn_,assign_:" = ", ltype_ : " double "]:=Module[{list=listIn//.powerRules/.ConvertRule[assign, ltype] },list=StringReplace[list,"$"->"_"];ColumnForm[list]];
SetOptions[Splice,FormatType->OutputForm,PageWidth->120];

(* ::Subsection:: ==================================================================================================================================================== *)
(*More Helpers*)

CF=ColumnForm;
Si=Simplify;
Fl=Flatten;
Jo=Join;
MF=MatrixForm;

DerivRule=Derivative[x__][zz__][xx__]:>ToExpression[StringJoin@@Join[{"D$",ToString[zz],"$"},ToString/@List[x]]][xx];
RemoveDeps := { f_@@IndVar:> f, f_@@fx:> f};

OComplement[all_List, i__List] :=
  DeleteDuplicates[Join @ ##] ~Drop~ Length[#] &[Union @ i, DeleteDuplicates @ all];

ZeroRules[rule__]:=OComplement[rule,DeleteCases[rule,r_->0]];

DeleteUnused[r_]:=Module[{rl={}},
                         Do[
                             Do[
                                 If[(Last[r[[i]]]/.r[[j]])!=Last[r[[i]]],True,True,rl=Jo[rl,{r[[j]]}]]
                                 ,{i,Length[r]}]
                             ,{j,Length[r]}];
                         Complement[r,rl]
                        ];


ProcessRules[r__]:=Module[{r0,r1=r,zrl,zr={}},
                          Print["ProcessRules start... "];

                          Do[
                              r0=r1;
                              r1=DeleteCases[r1,0->0];
                              zrl=ZeroRules[r1];
                              zr=Union[zr,zrl];
                              r1=OComplement[r1,zrl];
                              If[False,Print["ProcessRules iter= ",i];];
                              r1=r1//.zrl;
                              r1=DeleteCases[r1,0->0];

                              If[r0==r1,Break[]];
                              ,{i,10}];
                          Print["ProcessRules ... done "];

                          {r1,zr}
                         ];

rtest={xx->1,aa->0,cc->aa,dd->cc,ee->dd,mm->2ee+xx};
ProcessRules[rtest];

FirstDeriv[xx_,aa_]:=D[xx,#1]&/@aa;
SecondDeriv[xx_,aa_,bb_]:=Outer[D[xx,#1,#2]&,aa,bb];
break1stRule = "/* aa */" -> CHECKGRAD[getGrad];
break2ndRule = "/* aa */" -> CHECKHESS[getHess];
makeFirst[rl_]:=Fl[Jo[{rl,FirstDeriv[#,IndVar]&/@rl}]];
makeSecond[rl_]:=Fl[Jo[{makeFirst[rl],SecondDeriv[#,IndVar,IndVar]&/@rl}]];

frob[A_]:=Plus@@Flatten[A*A];

Print["Page 1"];
(* ::Subsection:: ==================================================================================================================================================== *)
(*Jacobians for each element type*)

<<"jac.m";

Print["Page 2"];
epr=True;

(* ::Subsection:: ==================================================================================================================================================== *)
(*Metrics Setup and Definition*)

(* set useLambda in the includer of this file Main.nb *)
lambdaList = {};
extraDof = 0;
If[useLambda, lambdaList = {Lambda0, Lambda1, Lambda2, Lambda3}; extraDof = 1; ];
If[useLambda,
   x=Transpose[Join[Table[{X0[j],X1[j],X2[j],X3[j]},{j,0,2}], {lambdaList}]];
   ,
   x=Transpose[Table[{X0[j],X1[j],X2[j],X3[j]},{j,0,2}]];
  ];
fx=Flatten[x];
IndVar:=fx;
alhs=Table[A[i,j]@@IndVar,{i,0,2},{j,0,2}];
alhsNoDep=Table[A[i,j],{i,0,2},{j,0,2}];
a=Table[ToExpression["A"<>ToString[i]<>ToString[j]]@@IndVar,{i,0,2},{j,0,2}];
fa=Fl[a];
w=Table[W[i,j],{i,0,2},{j,0,2}];
w=Table[ToExpression["W"<>ToString[i]<>ToString[j]],{i,0,2},{j,0,2}];
fw=Fl[w];
wi=Table[WI[i,j],{i,0,2},{j,0,2}];
wi=Table[ToExpression["WI"<>ToString[i]<>ToString[j]],{i,0,2},{j,0,2}];

If[epr, Print["Page 2a"]; ];

lhsTofA=Table[T[i,j]@@IndVar,{i,0,2},{j,0,2}];
lhsTofA=Table[ToExpression["T"<>ToString[i]<>ToString[j]]@@IndVar,{i,0,2},{j,0,2}];

If[epr, Print["Page 2b"]; ];

RD1m=Fl[Table[GRAD[II[inode], jc]->D[met@@IndVar,IndVar[[inode*(3 + extraDof)+jc+1]] ],{inode,0,3},{jc,0,2+extraDof}]];
If[epr, Print["Page 2b1"]; ];
RD2m=Fl[Table[HESS[II[inode],jc, II[jnode],kc]->D[met@@IndVar,IndVar[[inode*(3 + extraDof)+jc+1]],IndVar[[jnode*(3 + extraDof)+kc+1]] ],{inode,0,3},{jc,0,2 + extraDof},{jnode,0,3},{kc,0,2 + extraDof}]];

If[epr, Print["Page 2c"]; ];

Scmrhs=Table[Sc[i,j],{i,0,2},{j,0,2}];
Scm=Table[ToExpression["Sc"<>ToString[i]<>ToString[j]],{i,0,2},{j,0,2}];
RScm=Thread[Fl[Scm]->Fl[Scmrhs]];
Rbase=Fl[Jo[{
(*    Thread[fa->Fl[JM[Hex].Scm]],
    Thread[fw->Fl[WJM[Hex].Scm]], *)
    Thread[fa->Fl[JM[Hex]]],
    Thread[fw->Fl[WJM[Hex]]],
    detA@@IndVar->Si[Det[a]],
    detW->Si[Det[w]],
    detWI->1/detW,
    Thread[Fl[wi]->Fl[detWI Simplify[Det[w]Inverse[w]]]],
    Thread[Fl[lhsTofA]->Fl[a.wi]]
          }]];

Ra=Thread[Fl[alhsNoDep]->fa];

If[epr, Print["Page 2d"]; ];

Rdm=Fl[{RD1m,RD2m}];
ff=FirstDeriv[met@@IndVar,IndVar];

extraTerm = 0;
If[useLambda, extraTerm = onBoundary*(Lambda0 * Sum[(X0[j] - Y0[j])*normal[j], {j, 0, 2}]) + (1 - onBoundary)*Lambda0^2 ; ];

RUntangle=Fl[Jo[{
    Rbase,
    vU@@IndVar->beta detW - detA@@IndVar,
    hv->MyHeaviside[vU,0],
    met@@IndVar-> hv (vU@@IndVar)^2 + extraTerm
 }]];

If[epr, Print["Page 2e"]; ];

DimFac[2]=2;
DimFac[3]=3 Sqrt[3];

RB1[dim_]:=Fl[Jo[{
    Rbase,
    fT@@IndVar->frob[lhsTofA],

    fTDimO2@@IndVar -> Power[fT@@IndVar,dim/2],

    dimFacI -> 1/DimFac[dim],

    detAI@@IndVar -> 1/detA@@IndVar,

    met@@IndVar-> detW ( fTDimO2@@IndVar dimFacI detAI@@IndVar detW - 1 + extraTerm)

  (*    met@@IndVar->Power[fT@@IndVar,dim/2]/(DimFac[dim]detA@@IndVar/detW)-1 + extraTerm *)

               }]];

(* ::Text:: *)
(* dm/dx_ij = dm/dA_pq dA_pq/dx_ij *)

(* ::Text:: *)
(* d2m/dx_ij dx_lm = d2m/dA_pq dA_rs dA_pq/dx_ij dA_rs/dx_lm *)

If[False,

   dr=makeSecond[RUntangle];
   dru=RUntangle;
   {dr1,zr}=ProcessRules[dr];
   dr2=dr1/.DerivRule/.RemoveDeps;
  ];

Print["Page 3"];
(* ::Subsection:: ==================================================================================================================================================== *)
(*Metrics*)
If[True,

   {RUVal,RUValZR}=ProcessRules[RUntangle];
   {RU1stDeriv,RU1stDerivZR}=ProcessRules[makeFirst[RUntangle]];
   {RU2ndDeriv,RU2ndDerivZR}=ProcessRules[makeSecond[RUntangle]];

   {RB1Val,RB1ValZR}=ProcessRules[RB1[2]];
   {RB11stDeriv,RB11stDerivZR}=ProcessRules[makeFirst[RB1[2]]];
   {RB12ndDeriv,RB12ndDerivZR}=ProcessRules[makeSecond[RB1[2]]];

   {R3B1Val,R3B1ValZR}=ProcessRules[RB1[3]];
   {R3B11stDeriv,R3B11stDerivZR}=ProcessRules[makeFirst[RB1[3]]];
   {R3B12ndDeriv,R3B12ndDerivZR}=ProcessRules[makeSecond[RB1[3]]];

,
   {RUVal,RUValZR}=ProcessRules[RUntangle];
   {RU1stDeriv,RU1stDerivZR}=ProcessRules[makeFirst[RUntangle]];
   {RU2ndDeriv,RU2ndDerivZR}=ProcessRules[makeSecond[RUntangle]];

  ];


Print["Page 4"];

(* ::Subsection:: ==================================================================================================================================================== *)
(*ElimZeros*)

ElimZeros[R1_,ZR_]:=Module[{r1=R1,zr=ZR},
r1=r1//.zr;
zr=ZeroRules[r1];
r1=OComplement[r1,zr];
r1];

ElimZeros[Rdm,RU2ndDerivZR];
Print["Page 5"];

(* ::Subsection:: ==================================================================================================================================================== *)
(*Splice*)

Print["Splicing..."];
Splice["ScalingMatricesGen.mhpp", "ScalingMatricesGen.hpp"];
Print["Splicing...done"];

Print["Page 6"];
