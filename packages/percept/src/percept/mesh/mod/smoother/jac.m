

rr={X1[0]->1,X2[0]->0,X3[0]->0,X1[1]->0,X2[1]->1,X3[1]->0,X1[2]->0,X2[2]->0,X3[2]->1,X0[0]->0,X0[1]->0,X0[2]->0}


(* find S for each element type s.t. Atype = Abase . Stype *)
(* hex *)
Abase[0,0] = (X1[0] - X0[0]);
Abase[0,1] = (X2[0] - X0[0]);
Abase[0,2] = (X3[0] - X0[0]);

Abase[1,0] = (X1[1] - X0[1]);
Abase[1,1] = (X2[1] - X0[1]);
Abase[1,2] = (X3[1] - X0[1]);

Abase[2,0] = (X1[2] - X0[2]);
Abase[2,1] = (X2[2] - X0[2]);
Abase[2,2] = (X3[2] - X0[2]);
JM[Hex]      = Array[Abase,{3,3},{0,0}];
JMrr[Hex]      = Array[Abase,{3,3},{0,0}] /. rr;
SM[Hex]      = Inverse[JMrr[Hex]].JM[Hex] /. rr;

Wbase[0,0] = (WX1[0] - WX0[0]);
Wbase[0,1] = (WX2[0] - WX0[0]);
Wbase[0,2] = (WX3[0] - WX0[0]);

Wbase[1,0] = (WX1[1] - WX0[1]);
Wbase[1,1] = (WX2[1] - WX0[1]);
Wbase[1,2] = (WX3[1] - WX0[1]);

Wbase[2,0] = (WX1[2] - WX0[2]);
Wbase[2,1] = (WX2[2] - WX0[2]);
Wbase[2,2] = (WX3[2] - WX0[2]);
WJM[Hex]      = Array[Wbase,{3,3},{0,0}];
WJMrr[Hex]      = Array[Wbase,{3,3},{0,0}] /. rr;

(* pyramid *)
Apyr[0,0] = X1[0] - X0[0];
Apyr[0,1] = X2[0] - X0[0];
Apyr[0,2] = (2*X3[0] - X1[0] - X2[0])*0.5;

Apyr[1,0] = X1[1] - X0[1];
Apyr[1,1] = X2[1] - X0[1];
Apyr[1,2] = (2*X3[1] - X1[1] - X2[1])*0.5;

Apyr[2,0] = X1[2] - X0[2];
Apyr[2,1] = X2[2] - X0[2];
Apyr[2,2] = (2*X3[2] - X1[2] - X2[2])*0.5;
JM[Pyr]      = Array[Apyr,{3,3},{0,0}];
SM[Pyr]      = Inverse[JMrr[Hex]].JM[Pyr] /. rr;

(* wedge *)
Awedge[0,0] = X1[0] - X0[0];
Awedge[0,1] = isqrt3 * (2 * X2[0] - X1[0] - X0[0]);
Awedge[0,2] = X3[0] - X0[0];

Awedge[1,0] = X1[1] - X0[1];
Awedge[1,1] = isqrt3 * (2 * X2[1] - X1[1] - X0[1]);
Awedge[1,2] = X3[1] - X0[1];

Awedge[2,0] = X1[2] - X0[2];
Awedge[2,1] = isqrt3 * (2 * X2[2] - X1[2] - X0[2]);
Awedge[2,2] = X3[2] - X0[2];
JM[Wedge]      = Array[Awedge,{3,3},{0,0}];
SM[Wedge]      = Inverse[JMrr[Hex]].JM[Wedge] /. rr;

(* tet *)
Atet[0,0] = X1[0] - X0[0];
Atet[0,1] = (2*X2[0] - (X1[0] + X0[0]))*isqrt3;
Atet[0,2] = (3*X3[0] - X2[0] - (X1[0] + X0[0]))*isqrt6;

Atet[1,0] = X1[1] - X0[1];
Atet[1,1] = (2*X2[1] - (X1[1] + X0[1]))*isqrt3;
Atet[1,2] = (3*X3[1] - X2[1] - (X1[1] + X0[1]))*isqrt6;

Atet[2,0] = X1[2] - X0[2];
Atet[2,1] = (2*X2[2] - (X1[2] + X0[2]))*isqrt3;
Atet[2,2] = (3*X3[2] - X2[2] - (X1[2] + X0[2]))*isqrt6;
JM[Tet]      = Array[Atet,{3,3},{0,0}];
SM[Tet]      = Inverse[JMrr[Hex]].JM[Tet] /. rr;

(* tri *)
Atri[0,0] = X1[0] - X0[0];
Atri[0,1] = (2*X2[0] - X1[0] - X0[0])*isqrt3;
Atri[0,2] = 0;

Atri[1,0] = X1[1] - X0[1];
Atri[1,1] = (2*X2[1] - X1[1] - X0[1])*isqrt3;
Atri[1,2] = 0;

Atri[2,0] = 0;
Atri[2,1] = 0;
Atri[2,2] = 1;
JM[Tri]      = Array[Atri,{3,3},{0,0}];
SM[Tri]      = Inverse[JMrr[Hex]].JM[Tri] /. rr;

(* quad *)
Aquad[0,0] = (X1[0] - X0[0]);
Aquad[0,1] = (X2[0] - X0[0]);
Aquad[0,2] = 0  (* (X3[0] - X0[0]); *)

Aquad[1,0] = (X1[1] - X0[1]);
Aquad[1,1] = (X2[1] - X0[1]);
Aquad[1,2] = 0;  (* (X3[1] - X0[1]); *)

Aquad[2,0] = 0;  (* (X1[2] - X0[2]); *)
Aquad[2,1] = 0;  (* (X2[2] - X0[2]); *)
Aquad[2,2] = 1;  (* (X3[2] - X0[2]); *)
JM[Quad]      = Array[Aquad,{3,3},{0,0}];
SM[Quad]      = Inverse[JMrr[Hex]].JM[Quad] /. rr;

(* test1 = JM[Tri] -  JM[Hex].SM[Tri]  // Si *)
