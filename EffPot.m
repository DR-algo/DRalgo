(* ::Package:: *)

(* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ *)

(* :Title: EffPot                                                          	*)

(*
       This software is covered by the GNU General Public License 3.
       Copyright (C) 2021-2022 Andreas Ekstedt
       Copyright (C) 2021-2022 Philipp Schicho
       Copyright (C) 2021-2022 Tuomas V.I. Tenkanen

*)

(* :Summary:	Compute the effecitve potential.        					*)	

(* ------------------------------------------------------------------------ *)


(*
	User-defined model.
*)
DefineNewTensorsUS[\[Mu]ijI_,\[Lambda]4I_,\[Lambda]3I_,gvssI_,gvvvI_]:=Module[{\[Mu]ijP=\[Mu]ijI,\[Lambda]4P=\[Lambda]4I,\[Lambda]3P=\[Lambda]3I,gvssP=gvssI,gvvvP=gvvvI},
\[Mu]ijEP=\[Mu]ijP//SparseArray;
gvvvEP=gvvvP//SparseArray;
gvssEP=gvssP//SparseArray;
\[Lambda]4EP=\[Lambda]4P//SparseArray;
\[Lambda]3EP=\[Lambda]3P//SparseArray;
nsEP=Length[gvssEP[[1]]];
nvEP=Length[gvvvEP];

];


(*
	Defines parameters directly that are used by the effective potential
*)
DefineTensorsUS[]:=Module[{},
(*We only need to modify the gauge couplings*)
SubGauge=Table[c->Symbol[ToString[c]<>ToString["3dUS"]],{c,GaugeCouplingNames}];
gvvvEP=gvvv//Normal//ReplaceAll[#,SubGauge]&//SparseArray;

(*We only need to modify the gauge couplings*)
SubGauge=Table[c->Symbol[ToString[c]<>ToString["3dUS"]],{c,Variables[Normal[\[Lambda]4Light]]}];
\[Lambda]4EP=\[Lambda]4Light//Normal//ReplaceAll[#,SubGauge]&//SparseArray;
(*We only need to modify the gauge couplings*)
SubGauge=Table[c->Symbol[ToString[c]<>ToString["3dUS"]],{c,Variables[Normal[\[Lambda]3CLight]]}];
\[Lambda]3EP=\[Lambda]3CLight//Normal//ReplaceAll[#,SubGauge]&//SparseArray;
(*We only need to modify the gauge couplings*)
SubGauge=Table[c->Symbol[ToString[c]<>ToString["3dUS"]],{c,Variables[Normal[\[Mu]ijLight]]}];
\[Mu]ijEP=\[Mu]ijLight//Normal//ReplaceAll[#,SubGauge]&//SparseArray;

(*We only need to modify the gauge couplings*)
VarGauge=Table[Symbol[ToString[c]<>ToString["3d"]],{c,GaugeCouplingNames}];
SubGauge=Table[c->Symbol[ToString[c]<>ToString["US"]],{c,VarGauge}];
gvssEP=gvssL//Normal//ReplaceAll[#,SubGauge]&//SparseArray;

nsEP=Length[gvssEP[[1]]];
nvEP=Length[gvvvEP];

];


(*
	DiagonalMatrix was only introduced in Mathematica 12. So this is a hack for older versions.
	Did I mention that I am a fan of Mathematica? 
*)
DiagonalMatrixQAE[matI_]:=Module[{matP=matI},
VarMat=#->RandomInteger[10000]&/@(Normal[matP]//Variables);
matDia=Normal[DiagonalMatrix[Diagonal[matP]]]/.VarMat;
matT=Normal[matP]/.VarMat;
Return[matDia==matT]
];


{\[Mu]ijEP\[Phi],\[Mu]ijVec\[Phi],Gvvs\[Phi],\[Lambda]3\[Phi],\[Lambda]4\[Phi],gvvv\[Phi],gvss\[Phi],\[Phi]Vev};


(*
	Creates field-dependent couplings. These are created from the user's \[Phi]Vecp input.
*)
DefineVEVS[\[Phi]Vecp_]:=Module[{\[Phi]Vec=\[Phi]Vecp},

HabijVEP=Transpose[Activate@TensorContract[Inactive@TensorProduct[gvssEP,gvssEP],{{3,5}}],{1,3,2,4}]+Transpose[Transpose[Activate@TensorContract[Inactive@TensorProduct[gvssEP,gvssEP],{{3,5}}],{1,3,2,4}],{2,1,3,4}];

\[Phi]Vev=\[Phi]Vecp//SparseArray;
\[Mu]ijEP\[Phi]=\[Mu]ijEP+1/2 Activate@TensorContract[Inactive@TensorProduct[\[Lambda]4EP,\[Phi]Vec,\[Phi]Vec],{{3,5},{4,6}}]+ Activate@TensorContract[Inactive@TensorProduct[\[Lambda]3EP,\[Phi]Vec],{{3,4}}]//SparseArray;
\[Mu]ijVec\[Phi]=-1/2 Activate@TensorContract[Inactive@TensorProduct[HabijVEP,\[Phi]Vec,\[Phi]Vec],{{3,5},{4,6}}]//SparseArray;
Gvvs\[Phi]=-Activate@TensorContract[Inactive@TensorProduct[HabijVEP,\[Phi]Vec],{{4,5}}];
\[Lambda]3\[Phi]=\[Lambda]3EP+Activate@TensorContract[Inactive@TensorProduct[\[Lambda]4EP,\[Phi]Vec],{{4,5}}];
\[Lambda]4\[Phi]=\[Lambda]4EP//SparseArray;
gvvv\[Phi]=gvvvEP//SparseArray;
gvss\[Phi]=gvssEP//SparseArray;

If[DiagonalMatrixQAE[\[Mu]ijVec\[Phi]]==False,
Print["The Vector Mass-Matrix is not Diagonal"];
];


If[DiagonalMatrixQAE[\[Mu]ijEP\[Phi]]==False,
Print["The Scalar Mass-Matrix is not Diagonal"];
];

];


Options[CalculatePotentialUS] = {CustomMasses -> False}


(*
	Calculates the effective potential with custom masses
*)
CalculatePotentialUS[ScalMassI_,VecMassI_,OptionsPattern[]]:=Module[{ScalMassP=ScalMassI,VecMassP=VecMassI},

CustomMass= OptionValue[CustomMasses];
If[CustomMass==True,
\[Mu]ijEP\[Phi]=ScalMassP//SparseArray;
\[Mu]ijVec\[Phi]=VecMassP//SparseArray;

If[DiagonalMatrixQAE[\[Mu]ijEP\[Phi]]==True && DiagonalMatrixQAE[\[Mu]ijVec\[Phi]]==True,
CalculateLOPotentialSS[];
CalculateNLOPotentialSS[];

VTot={VLO,VNLO};
,
Print["The Mass matrices are not diagonal. Please rotate to the mass-basis using RotateTensorsUSPostVeV[]"];
];
,
Print["Please set CustomMasses->True"]
];


];


(*
	Calculates the effective potential.
*)
CalculatePotentialUS[OptionsPattern[]]:=Module[{},

CustomMass= OptionValue[CustomMasses];
If[CustomMass==False,
If[DiagonalMatrixQAE[\[Mu]ijEP\[Phi]]==True && DiagonalMatrixQAE[\[Mu]ijVec\[Phi]]==True,
CalculateLOPotentialSS[];
CalculateNLOPotentialSS[];
CalculateNNLOPotentialSS[];

VTot={VLO,VNLO,VNNLO};
,
Print["The Mass matrices are not diagonal. Please rotate to the mass-basis using RotateTensorsUSPostVeV[]"];
];
,
Print["Please supply scalar and vector mass matrices"]
];


];


(*
	Prints the effective potential.
*)
PrintEffectivePotential[optP_]:=Module[{opt=optP},
EffPotPrint=Switch[opt,"LO",VTot[[1]],"NLO",VTot[[2]],"NNLO",VTot[[3]]];

(*Printing Result*)
ToExpression[StringReplace[ToString[StandardForm[EffPotPrint]],"DRalgo`Private`"->""]]
];


{VLO,VNLO,VNNLO};


(*
	Calculates the two-loop effective potential.
*)
CalculateNNLOPotentialSS[]:=Module[{},
If[verbose==True,Print["Calculating the 2-Loop Effective Potential"]];

(*The notation follows Martin's convention*)
(*Please see arXiv:1808.07615*)

Q=\[Mu]3US; (*RG scale*)

(*Definitions*)
\[Sigma]=1/(16 \[Pi]^2);
f[x_]:=-(1/(12 \[Pi])) x^(3/2);
A[0]=0;
A[x_]:=-(1/(4 \[Pi])) Sqrt[x];
I2[0,0,0]=0;
I2[x_,y_,z_]:=1/(4 \[Pi])^2 (Log[Q/(Sqrt[x]+Sqrt[y]+Sqrt[z])]+1/2);
I2Div[x_,y_,z_]:=1/(4 \[Pi])^2 1/4;
\[CapitalLambda][x_,y_,z_]:=x^2+y^2+z^2-2x y-2x z-2 y z;
fsss[x_,y_,z_]:=-I2[x,y,z];
fss[x_,y_]:=A[x]A[y];
fvs[x_,y_]:=2 A[x]A[y];
fssv[x_,y_,z_]:=1/z (-\[CapitalLambda][x,y,z]I2[x,y,z]+(x-y)^2 I2[x,y,0]+z A[x] A[y]+(y-x-z)A[x] A[z]+(x-y-z)A[y] A[z]+(x-y)A[x] A[0]+(y-x)A[y] A[0]);
fssv[x_,y_,0]:=(x+4 Sqrt[x] Sqrt[y]+y+4 (x+y) Log[Q/(Sqrt[x]+Sqrt[y])])/(32 \[Pi]^2);
fssv[0,0,0]:=0;

fvbvbs[x_,y_,z_]:=1/(4 x y) (-(x+y-z)^2I2[x,y,z]+(x-z)^2 I2[x,0,z]+(y-z)^2 I2[0,y,z]-z^2 I2[0,0,z]+(z-x-y)A[x] A[y]+y A[x] A[z]+x A[y] A[z]+z A[0] A[0]+(x-z)A[x]A[0]+(y-z)A[y]A[0]-(x+y)A[z]A[0]);
fvbvbs[0,y_,z_]:=1/(128 \[Pi]^2 y) (-3 y+4 Sqrt[y] Sqrt[z]-4 y Log[Q/(Sqrt[y]+Sqrt[z])]+4 z Log[Q/(Sqrt[y]+Sqrt[z])]-4 z Log[Q/Sqrt[z]]);
fvbvbs[x_,0,z_]:=fvbvbs[0,x,z];
fvbvbs[0,0,z_]:=(-1-4 Log[Q/Sqrt[z]])/(128 \[Pi]^2);
fvbvbs[0,0,0]:=0;
fvvs[x_,y_,z_]:=fvbvbs[x,y,z]-I2[x,y,z]+2I2Div[x,y,z];
fvv[x_,y_]:=8/3 A[x]A[y];
fvvv[x_,y_,z_]:=1/(4x y z) (-\[CapitalLambda][x,y,z](\[CapitalLambda][x,y,z] I2[x,y,z]+4(x y+ x z+ y z)(2 I2[x,y,z]-2 I2Div[x,y,z]))+(x-y)^2 ((x^2+y^2+6x y)I2[x,y,0]-8 x y I2Div[x,y,0])-z^4 I2[z,0,0]+(x-z)^2 ((x^2+z^2+6x z)I2[x,z,0]-8 x z I2Div[x,0,z])-y^4 I2[y,0,0]+(y-z)^2 ((y^2+z^2+6y z)I2[y,z,0]-8y z I2Div[0,y,z])-x^4 I2[x,0,0]+ z A[x]A[y](z^2-5(x^2+y^2-x z -y z)-26/3x y)+ y A[x]A[z](y^2-5(x^2+z^2-x y -y z)-26/3x z)+ x A[z]A[y](x^2-5(z^2+y^2-x z -y x)-26/3z y)-A[x]A[0](-9 x^2 (y+z)+9x (y^2+z^2)+y^3+z^3)-A[y]A[0](-9 y^2 (x+z)+9y (x^2+z^2)+x^3+z^3)-A[z]A[0](-9 z^2 (x+y)+9z (x^2+y^2)+x^3+y^3)+(x^3+y^3+z^3)A[0]A[0]);
fvvv[x_,y_,0]:=1/(384 \[Pi]^2 x y) (-48 x^(5/2) Sqrt[y]+27 x^2 y-18 x^(3/2) y^(3/2)+27 x y^2-48 Sqrt[x] y^(5/2)-136 \[Pi]^2 x^2 y \[Sigma]-136 \[Pi]^2 x y^2 \[Sigma]+26 x y (x+y) Log[(E^EulerGamma Glaisher^12 Q^2)/(16 \[Pi]^2 T^2)]+48 x^3 Log[Q/Sqrt[x]]-48 x^3 Log[Q/(Sqrt[x]+Sqrt[y])]+192 x^2 y Log[Q/(Sqrt[x]+Sqrt[y])]+192 x y^2 Log[Q/(Sqrt[x]+Sqrt[y])]-48 y^3 Log[Q/(Sqrt[x]+Sqrt[y])]+48 y^3 Log[Q/Sqrt[y]]);
fvvv[x_,0,z_]:=fvvv[x,z,0];
fvvv[0,y_,z_]:=fvvv[y,z,0];
fvvv[0,0,0]:=0;
fvvv[x_,0,0]:=-(1/(384 \[Pi]^2)) x (-3+136 \[Pi]^2 \[Sigma]-26 Log[(E^EulerGamma Glaisher^12 Q^2)/(16 \[Pi]^2 T^2)]-192 Log[Q/Sqrt[x]]);
fvvv[0,y_,0]:=fvvv[y,0,0];
fvvv[0,0,z_]:=fvvv[z,0,0];
fvvv[0,0,0]:=0;
f\[Eta]\[Eta]v[x_,y_,z_]:=2 1/(4z) (\[CapitalLambda][x,y,z]I2[x,y,z]-(x-y)^2 I2[x,y,0]-z A[x]A[y]+(x-y+z) A[x]A[z]+(y-x+z)A[y]A[z]-(x-y)A[0](A[x]-A[y]));
f\[Eta]\[Eta]v[x_,y_,0]:=-(((x+4 Sqrt[x] Sqrt[y]+y+4 (x+y) Log[Q/(Sqrt[x]+Sqrt[y])]))/(64 \[Pi]^2));
f\[Eta]\[Eta]v[0,0,0]:=0;

f\[Eta]\[Eta]v[0,0,z]:=(z (Log[Q/Sqrt[z]]+1/2))/(32 \[Pi]^2);

aS=Table[\[Mu]ijEP\[Phi][[i,i]],{i,1,nsEP}]//SparseArray;
av=Table[\[Mu]ijVec\[Phi][[i,i]],{i,1,nvEP}]//SparseArray;
(*Potential*)
ss=1/8 TensorProduct[\[Lambda]4\[Phi]];
Vss=Sum[ss[[j,j,k,k]]fss[aS[[j]],aS[[k]]],{j,nsEP},{k,nsEP}];
sss=1/12 TensorProduct[\[Lambda]3\[Phi],\[Lambda]3\[Phi]];
Vsss=Sum[sss[[i,j,k,i,j,k]]fsss[aS[[i]],aS[[j]],aS[[k]]],{j,nsEP},{k,nsEP},{i,nsEP}];
vvs=1/4 TensorProduct[Gvvs\[Phi],Gvvs\[Phi]];
Vvvs=Sum[vvs[[a,b,i,a,b,i]]fvvs[av[[a]],av[[b]],aS[[i]]],{a,nvEP},{b,nvEP},{i,nsEP}];
ssv=1/4 TensorProduct[gvss\[Phi],gvss\[Phi]];
Vssv=Sum[ssv[[a,i,j,a,i,j]]fssv[aS[[i]],aS[[j]],av[[a]]],{a,nvEP},{j,nsEP},{i,nsEP}];
vs=1/2 TensorProduct[gvss\[Phi],gvss\[Phi]];
Vvs=Sum[vs[[a,i,j,a,i,j]]fvs[aS[[i]],av[[a]]],{a,nvEP},{j,nsEP},{i,nsEP}];
vv=1/4 TensorProduct[gvvv\[Phi],gvvv\[Phi]];
Vvv=Sum[vv[[a,b,c,a,b,c]]fvv[av[[a]],av[[b]]],{a,nvEP},{b,nvEP},{c,nvEP}];
vvv=1/12 TensorProduct[gvvv\[Phi],gvvv\[Phi]];
Vvvv=Sum[vvv[[a,b,c,a,b,c]]fvvv[av[[a]],av[[b]],av[[b]]],{a,nvEP},{b,nvEP},{c,nvEP}];
ggv=1/4 TensorProduct[gvvv\[Phi],gvvv\[Phi]];
V\[Eta]\[Eta]v=Sum[ggv[[a,b,c,a,b,c]]f\[Eta]\[Eta]v[0,0,av[[b]]],{a,nvEP},{b,nvEP},{c,nvEP}];

VNNLO= Vss+ Vsss+ Vvvs+ Vssv+ Vvs+ Vvvv+ Vvv+ V\[Eta]\[Eta]v;

];


(*
	Calculates the one-loop effective potential.
*)
CalculateNLOPotentialSS[]:=Module[{},
If[verbose==True,Print["Calculating the 1-Loop Effective Potential"]];

ALog[x_]:=-(x^(3/2)/(12 \[Pi]));

V1=Sum[ALog[\[Mu]ijEP\[Phi][[i,i]]],{i,1,nsEP}];
V2=2*Sum[ALog[\[Mu]ijVec\[Phi][[i,i]]],{i,1,nvEP}];

VNLO=V1+V2;
];


(*
	Calculates the tree-level effective potential.
*)
CalculateLOPotentialSS[]:=Module[{},
If[verbose==True,Print["Calculating the Tree-Level Effective Potential"]];

V1=\[Lambda]4EP . \[Phi]Vev . \[Phi]Vev . \[Phi]Vev . \[Phi]Vev;
V2=\[Mu]ijEP . \[Phi]Vev . \[Phi]Vev;
V3=\[Lambda]3EP . \[Phi]Vev . \[Phi]Vev . \[Phi]Vev;

VLO=1/4! V1+V2/2!+V3/3!;
];


(*
	Prints field-dependent tensors.
*)
PrintTensorsVEV[]:=Module[{},
Print["(1) Scalar Mass, (2) Vector Mass, (3) Scalar Quartic, (4) Scalar Cubic, (5) VSS, (6) VVS, (7) VVV"];

ToExpression[StringReplace[ToString[StandardForm[Join[\[Mu]ijEP\[Phi],\[Mu]ijVec\[Phi],\[Lambda]4\[Phi],\[Lambda]3\[Phi],gvss\[Phi],Gvvs\[Phi],gvvv\[Phi]]]],"DRalgo`Private`"->""]]

];


(*
	Rotates to a diagonal-mass basis.
*)
RotateTensorsUSPostVEV[DScalarsp_,DVectorsp_]:=Module[{DS=DScalarsp,DV=DVectorsp},
DS=DS//SparseArray;
DV=DV//SparseArray;

\[Lambda]4\[Phi]=Transpose[DS] . \[Lambda]4\[Phi] . DS//Activate@TensorContract[Inactive@TensorProduct[DS,#],{{1,4}}]&//Transpose[#,{2,1,3,4}]&//Activate@TensorContract[Inactive@TensorProduct[DS,#],{{1,5}}]&//Transpose[#,{3,2,1,4}]&//SimplifySparse;
\[Lambda]3\[Phi]=Transpose[DS] . \[Lambda]3\[Phi] . DS//Activate@TensorContract[Inactive@TensorProduct[DS,#],{{1,4}}]&//Transpose[#,{2,1,3}]&//SimplifySparse;
gvss\[Phi]=Transpose[DV] . gvss\[Phi] . DS//Activate@TensorContract[Inactive@TensorProduct[DS,#],{{1,4}}]&//Transpose[#,{2,1,3}]&//SimplifySparse;
Gvvs\[Phi]=Transpose[DV] . Gvvs\[Phi] . DS//Activate@TensorContract[Inactive@TensorProduct[DV,#],{{1,4}}]&//Transpose[#,{2,1,3}]&//SimplifySparse;
gvvv\[Phi]=Transpose[DV] . gvvv\[Phi] . DV//Activate@TensorContract[Inactive@TensorProduct[DV,#],{{1,4}}]&//Transpose[#,{2,1,3}]&//SimplifySparse;
\[Mu]ijVec\[Phi]=Activate@TensorContract[Inactive@TensorProduct[DV,DV,\[Mu]ijVec\[Phi]],{{1,5},{3,6}}]//SimplifySparse;
\[Mu]ijEP\[Phi]=Activate@TensorContract[Inactive@TensorProduct[DS,DS,\[Mu]ijEP\[Phi]],{{1,5},{3,6}}]//SimplifySparse;


If[DiagonalMatrixQAE[\[Mu]ijVec\[Phi]]==False,
Print["The Vector mass-Matrix is not diagonal"];
];


If[DiagonalMatrixQAE[\[Mu]ijEP\[Phi]]==False,
Print["The Scalar mass-Matrix is not diagonal"];
];

];


(*
	Rotates to a diagonal-mass basis.
*)
RotateTensorsCustomMass[DScalarsp_,DVectorsp_,ScalarMass_,vectorMass_]:=Module[{DS=DScalarsp,DV=DVectorsp},
DS=DS//SparseArray;
DV=DV//SparseArray;

\[Lambda]4\[Phi]=Transpose[DS] . \[Lambda]4\[Phi] . DS//Activate@TensorContract[Inactive@TensorProduct[DS,#],{{1,4}}]&//Transpose[#,{2,1,3,4}]&//Activate@TensorContract[Inactive@TensorProduct[DS,#],{{1,5}}]&//Transpose[#,{3,2,1,4}]&//SimplifySparse;
\[Lambda]3\[Phi]=Transpose[DS] . \[Lambda]3\[Phi] . DS//Activate@TensorContract[Inactive@TensorProduct[DS,#],{{1,4}}]&//Transpose[#,{2,1,3}]&//SimplifySparse;
gvss\[Phi]=Transpose[DV] . gvss\[Phi] . DS//Activate@TensorContract[Inactive@TensorProduct[DS,#],{{1,4}}]&//Transpose[#,{2,1,3}]&//SimplifySparse;
Gvvs\[Phi]=Transpose[DV] . Gvvs\[Phi] . DS//Activate@TensorContract[Inactive@TensorProduct[DV,#],{{1,4}}]&//Transpose[#,{2,1,3}]&//SimplifySparse;
gvvv\[Phi]=Transpose[DV] . gvvv\[Phi] . DV//Activate@TensorContract[Inactive@TensorProduct[DV,#],{{1,4}}]&//Transpose[#,{2,1,3}]&//SimplifySparse;

\[Mu]ijVec\[Phi]=vectorMass//SparseArray;
\[Mu]ijEP\[Phi]=ScalarMass//SparseArray;


If[DiagonalMatrixQAE[\[Mu]ijVec\[Phi]]==False,
Print["The Vector mass-Matrix is not diagonal"];
];


If[DiagonalMatrixQAE[\[Mu]ijEP\[Phi]]==False,
Print["The Scalar mass-Matrix is not diagonal"];
];

];
