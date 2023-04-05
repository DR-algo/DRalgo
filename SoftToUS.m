(* ::Package:: *)

(* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ *)

(* :Title: SoftToUS                                                         *)

(*
       This software is covered by the GNU General Public License 3.
       Copyright (C) 2021-2022 Andreas Ekstedt
       Copyright (C) 2021-2022 Philipp Schicho
       Copyright (C) 2021-2022 Tuomas V.I. Tenkanen

*)

(* :Summary:	Dimensonal reduction from soft to supersoft scale.            *)	

(* ------------------------------------------------------------------------ *)


(*
	Prints the result from SymmEnergy.
*)
PrintPressureUS[optP_]:=Module[{opt=optP},
SymmPrint=Switch[opt,"LO",SymmEnergyUS[[1]],"NLO",SymmEnergyUS[[2]]];


(*Printing Result*)
ToExpression[StringReplace[ToString[StandardForm[SymmPrint]],"DRalgo`Private`"->""]]
];


(*
	Calculates the preassure in the ultrasoft theory. Only the preassure in the symmetric
	phae is calculated.
*)
SymmetricPhaseEnergyUS[]:=Module[{},
(*
	Counterterms are needed to calculate
	SymmetricPhaseNLO and SymmetricPhaseNNLO
*)
(*The minus signs is a convention to get the pressure*)
Tot={-SymmetricPhaseUSLO[],-SymmetricPhaseUSNLO[]};
SymmEnergyUS=Tot;
];


(*
	Calculates the 1-loop pressure in the ultrasoft theory.
*)
SymmetricPhaseUSLO[]:=Module[{},
If[verbose,Print["Calculating Leading-Order \!\(\*SuperscriptBox[\(T\), \(4\)]\) Terms"]];

ContriScalars=Sum[-1/(12 \[Pi]) \[Mu]ijL[[i,i]]^3,{i,1,Length[\[Mu]ijL]}];

ToExpression[StringReplace[ToString[StandardForm[ContriScalars]],"DRalgo`Private`"->""]]
];



(*
	Calculates the 2-loop pressure in the ultrasoft theory.
*)
SymmetricPhaseUSNLO[]:=Module[{},
If[verbose,Print["Calculating NLO \!\(\*SuperscriptBox[\(T\), \(4\)]\) Terms"]];
(*Definitions*)
fSSV[x_,y_]:=(4 (x^2+y^2) Log[\[Mu]3/(x+y)]+4 Sqrt[x^2] Sqrt[y^2]+x^2+y^2)/(32 \[Pi]^2);

I1Temp=1/(16^2 \[Pi]);
Vss=1/8*I1Temp*TensorContract[\[Mu]ijL . \[Lambda]KTotal . \[Mu]ijL,{{1,2},{3,4}}];

ssv=1/4 TensorProduct[gvssVTot,gvssVTot];
Vssv=Sum[ssv[[a,i,j,a,i,j]]fSSV[\[Mu]ijL[[i,i]],\[Mu]ijL[[j,j]]],{a,nv},{i,nSH},{j,nSH}];


VNLO= Vss+Vssv;

ToExpression[StringReplace[ToString[StandardForm[VNLO]],"DRalgo`Private`"->""]]
];



(*
	Prints all the couplings in the supersoft theory.	
*)




(*
	Prints all the couplings in the supersoft theory.	
*)

PrintCouplingsUS[]:=Module[{},

(*Non-Abelian Couplings*)

NonAbelianCouplingSS[]; (*Calculates non-abelian couplings*)

GabcdTemp=GgvvvSS . gvvvSS+gvvvSS . GgvvvSS;
 Temp=Activate @ TensorContract[
        Inactive[TensorProduct][gvvvSS,gvvvSS], {{3, 6}}]    ; 
GabVTree=TensorContract[Temp,{{2,3}}]//Normal;
GabVLoop=TensorContract[GabcdTemp,{{2,3}}]//Normal;


HelpList=DeleteDuplicates@Flatten@FullSimplify[GabVTree+GabVLoop ]//Sort; (*Adds the tree-level result*)
HelpVar=Table[ \[Lambda]VNASS[a],{a,1,Delete[HelpList,1]//Length}];
HelpVarMod=RelationsBVariables[HelpList,HelpVar]; (*Finds a minimal basis of couplings*)
HelpSolveNASS=Table[{Delete[HelpList,1][[a]]->HelpVarMod[[a]]},{a,1,Delete[HelpList,1]//Length}]//Flatten;
\[Lambda]VecNASS=GabVTree+GabVLoop//Normal//FullSimplify//ReplaceAll[#,HelpSolveNASS]&//SparseArray;
IdentMatNASS=List/@HelpSolveNASS/.{b_->a_}:>a->b//Flatten[#,1]&;



VarGauge=Table[Symbol[ToString[c]<>ToString["3d"]],{c,GaugeCouplingNames}];
SubGauge=Table[c->Symbol[ToString[c]<>ToString["US"]],{c,VarGauge}];


A1=\[Lambda]VecNASS//Normal;
A2=GabVTree//Normal//ReplaceAll[#,SubGauge]&;
Var3D=A2//Variables;
RepVar3D=#->Sqrt[#]&/@Var3D;
A2Mod=A2/.RepVar3D;
Sol1=Solve[A2Mod==A1,Var3D]/.IdentMatNASS//Flatten[#,1]&//FullSimplify; (*Solves ultrasoft couplings in terms of soft ones*)
ResGaugeNASS=Table[List[Sol1[[c]]]/.{b_->a_}:>b^2->a,{c,1,Length[Sol1]}]//Simplify;




(* Gauge couplings*)
A1=TensorContract[\[Lambda]KVecTSS,{{3,4}}]//Normal;
A2=TensorContract[HabijVL,{{3,4}}]//Normal//ReplaceAll[#,SubGauge]&;
(*Trick to avoid problems when kinetic mixing*)
A1=DiagonalMatrix[Diagonal[A1]];
A2=DiagonalMatrix[Diagonal[A2]];


Var3D=VarGauge//ReplaceAll[#,SubGauge]&//Variables;
RepVar3D=#->Sqrt[#]&/@Var3D;
A2Mod=A2/.RepVar3D;
Sol1=Solve[A2Mod==A1,Var3D]/.IdentMatSS//Flatten[#,1]&;

ResGauge=Table[List[Sol1[[c]]]/.{b_->a_}:>b^2->a,{c,1,Length[Sol1]}];



(* Scalar Quartics*)
QuarticVar=\[Lambda]4S//Normal//Variables;
SubGauge=Table[c->Symbol[ToString[c]<>ToString["US"]],{c,QuarticVar}];
NonZeroPos=SparseArray[\[Lambda]4S]["NonzeroPositions"];
SolVar=Extract[\[Lambda]4S-\[Lambda]3DSS,{#}]&/@NonZeroPos//DeleteDuplicates;
ResScalp=Reduce[SolVar==0,QuarticVar]//ToRules[#]&;
SolveTemp=QuarticVar/.ResScalp;
ResScal=Table[{ReplaceAll[QuarticVar[[i]],SubGauge]->SolveTemp[[i]]},{i,1,Length@QuarticVar}]//Flatten[#,1]&//ReplaceAll[#,IdentMatSS]&;


(* Scalar Cubics*)
VarGauge=Join[\[Lambda]3CLight//Normal//Variables]//DeleteDuplicates;
SubGauge=Table[c->Symbol[ToString[c]<>ToString["US"]],{c,VarGauge}];
\[Lambda]3p=\[Lambda]3CLight//Normal//ReplaceAll[#,SubGauge]&;
SolVar=\[Lambda]3CSRedSS-\[Lambda]3p//Normal;
CubicVar=\[Lambda]3p//Normal//Variables;
ResCubicSS=Solve[SolVar==0,CubicVar]/.IdentMatSS//Flatten[#,1]&;




(*Printing Result*)
PrintPre=Join[ResScal,ResGauge,ResGaugeNASS,ResCubicSS]//Normal//FullSimplify//DeleteDuplicates;
PrintPre=DeleteCases[PrintPre,0->0];

ToExpression[StringReplace[ToString[StandardForm[PrintPre]],"DRalgo`Private`"->""]]

];


(*
	Creates auxiliary tensors that appear in the matching.
*)
CreateHelpTensorsSS[]:=Module[{},
If[verbose,Print["Creating Help Tensors"]];


HabijL=Transpose[Activate @ TensorContract[
        Inactive[TensorProduct][gvssL,gvssL], { {3, 5}}],{1,3,2,4}]//SimplifySparse;
HabijVL=HabijL+Transpose[HabijL,{2,1,3,4}]//SimplifySparse;
SelfEnergySS=  Inactivate[TensorProduct[gAvss,gAvss]]//SimplifySparse;
HabijA=Transpose[Activate@TensorContract[SelfEnergySS,{{3,5}}],{1,3,2,4}]//SimplifySparse;
HabijVA=HabijA+Transpose[HabijA,{1,2,4,3}]//SimplifySparse;
];


(*
	Scalar self-energy in the effective theory.
*)
ScalarSelfEnergySS[]:=Module[{},
If[verbose,Print["Calculating Scalar Self-Energy"]];

SelfEnergySS=-1/(12\[Pi]);
ContriSS=SelfEnergySS/2*Simplify[Table[Sum[\[Lambda]3Cx[[i,ii,jj]]\[Lambda]3Cx[[j,ii,jj]]/(\[Mu]ijL[[ii,ii]]+\[Mu]ijL[[jj,jj]])^3,{ii,1,nSH},{jj,1,nSH}],{i,1,nSL},{j,1,nSL}]];
ContriSS2=SelfEnergySS*Simplify[Table[Sum[\[Lambda]3Cy[[i,ii,jj]]\[Lambda]3Cy[[j,ii,jj]]/(\[Mu]ijL[[jj,jj]])^3,{ii,1,nSL},{jj,1,nSH}],{i,1,nSL},{j,1,nSL}]];


ZSij=-(ContriSS+ContriSS2)/2;
];


(*
	Matching of scalar-cubic couplings. \[Lambda]3Cy corresponds to light*light*heavy scalar cubic, and \[Lambda]3Cx corresponds to light*heavy*heavy scalar coupling.
*)

ScalarCubicsSS[]:=Module[{},
If[verbose,Print["Calculating Scalar Cubic Couplings"]];

If[Length[\[Lambda]3//Normal//Variables]==0,
\[Lambda]3CSSS=0;
,
MassVec=Table[\[Mu]ijL[[n,n]],{n,1,nSH}];

TensHelp=Table[1/n^2,{n,MassVec}]//SparseArray;
TensHelp2=TensorProduct[\[Lambda]3Cy,\[Lambda]3CHeavy]//DiagonalTensor2[#,3,6]&;
Temp=-TensHelp . TensHelp2//Flatten[#,{{3},{4},{1},{2}}]&;
\[Lambda]KTemp=\[Lambda]K+Temp;

TensHelp=Table[1/n^2,{n,MassVec}]//SparseArray;
TensHelp2=TensorProduct[\[Lambda]3Cx,\[Lambda]3Cy]//DiagonalTensor2[#,3,6]&;
Temp=-TensHelp . TensHelp2//Transpose[#,{1,4,3,2}]&;
Temp=3Symmetrize[Temp,Symmetric[{1,2,3}]]//SparseArray//SimplifySparse;
\[Lambda]yTemp=\[Lambda]y+Temp;

TensHelp=Table[1/(n+m),{n,MassVec},{m,MassVec}]//SparseArray;
TensHelp2=TensorProduct[\[Lambda]KTemp,\[Lambda]3Cx]//DiagonalTensor2[#,1,6]&//DiagonalTensor2[#,2,6]&;
Contri1Pre=Total[TensHelp TensHelp2,-4];
Contri1=1/(4 \[Pi]) /2*3*Symmetrize[Contri1Pre,Symmetric[{1,2,3}]]//SparseArray//SimplifySparse;

TensHelp=Table[1/(n),{n,MassVec}]//SparseArray;
TensHelp2=TensorProduct[\[Lambda]yTemp,\[Lambda]3Cy]//TensorProduct[#,{3,6}]&//DiagonalTensor2[#,4,7]&;
Contri2Pre=TensHelp . TensHelp2;
Contri2=1/(4 \[Pi])*3*Symmetrize[Contri2Pre,Symmetric[{1,2,3}]]//SparseArray//SimplifySparse;

MassVec2=Table[\[Mu]ijLS[[n,n]],{n,1,ns}];
TensHelp=Table[MassHelpTriangle[n,m,l],{n,MassVec2},{m,MassVec2},{l,MassVec2}]//SparseArray;
TensHelp2=TensorProduct[\[Lambda]3CTot,\[Lambda]3CTot]//DiagonalTensor2[#,2,6]&//TensorProduct[#,\[Lambda]3CTot]&//DiagonalTensor2[#,5,7]&//DiagonalTensor2[#,3,7]&;
Contri3Pre=Total[TensHelp TensHelp2,-4][[LightScalar[[;;,1]],LightScalar[[;;,1]],LightScalar[[;;,1]]]];
Contri3=-3*Symmetrize[Contri3Pre,Symmetric[{1,2,3}]]//SparseArray//SimplifySparse;

TensHelp=Table[n/m^2,{n,MassVec},{m,MassVec}]//SparseArray;
TensHelp2=TensorProduct[\[Lambda]3Cy,\[Lambda]x]//DiagonalTensor2[#,5,6]&//DiagonalTensor2[#,4,6]&;
Temp=TensHelp . TensHelp2//TensorContract[#,{{1,2}}]&;
ContriMixed=-1/(4 \[Pi])*1/2*3 Symmetrize[Temp,Symmetric[{1,2,3}]]//SparseArray//SimplifySparse;

ContriSE=3 Symmetrize[ZSij . \[Lambda]3CLight,Symmetric[{1,2,3}]]//SparseArray//SimplifySparse;

\[Lambda]3CSSS=-Contri3-Contri1-Contri2-ContriSE-ContriMixed;
];

];


(*
	Calculates the 2-loop scalar mass in the ultrasoft theory.
*)
ScalarMass2LoopSS[]:=Module[{},
If[verbose,Print["Calculating 2-Loop Scalar Mass"]];

MassVector=Table[\[Mu]ijL[[n,n]],{n,1,nSH}];

TensHelp=Table[(1/2+Log[\[Mu]3/( a+b)]),{a,MassVector},{b,MassVector}]//SparseArray;
AE1=\[Lambda]K . Transpose[\[Lambda]K,{4,2,3,1}]//DiagonalTensor[#,1,4]&//DiagonalTensor[#,1,5]&;
Contri1=-1/(16 \[Pi]^2)1/2AE1 . TensHelp//TensorContract[#,{2,4}]&;

MassHelp=Table[(1/2+2Log[\[Mu]3/(2 a)]),{a,MassVector}]//SparseArray;
TensHelp=Transpose[Transpose[gAvss,{2,1,3}],{1,3,2}] . gAvss//TensorContract[#,{1,3}]&;
TensHelp2=MassHelp TensHelp . \[Lambda]K//TensorContract[#,{1,2}]&;
Contri2=1/(16 \[Pi]^2)*(1/2)TensHelp2;

MassHelp=Table[(-Log[\[Mu]3/(2 a)]),{a,MassVector}]//SparseArray;
TensHelp=gAvss . Transpose[gAvss,{2,1,3}]//SparseArray //DiagonalTensor[#,2,4]&//Transpose[#,{3,2,1}]&;
TensHelp2=TensHelp . MassHelp . HabijVL//TensorContract[#,{1,2}]&;
Contri3=1/(16 \[Pi]^2)*(1/4)TensHelp2;

MassHelp=Table[1/(a+b),{a,MassVector},{b,MassVector}]//SparseArray;
TensHelp=TensorProduct[DiagonalTensor[\[Lambda]4K,3,4] . MassVector,MassHelp]//DiagonalTensor[#,1,3]&//DiagonalTensor[#,1,3]&;
Contri4=1/4*(1/(16 \[Pi]^2))*TensHelp . \[Lambda]K//TensorContract[#,{1,2}]&;

(*
Contri1=-1/(16 \[Pi]^2)1/2Simplify[Table[Sum[ \[Lambda]K[[a,b,i,n]]\[Lambda]K[[a,b,n,j]](1/2+Log[\[Mu]3US/( \[Mu]ijL[[a,a]]+ \[Mu]ijL[[b,b]])]),{a,1,nSH},{b,1,nSH},{n,1,nSL}],{i,1,nSL},{j,1,nSL}]];
Contri2=1/(16 \[Pi]^2)*(-1/2)Simplify[Table[Sum[\[Lambda]K[[n,m,i,j]]gAvss[[a,n,l]]gAvss[[a,l,m]](1/2+2Log[\[Mu]3US/(2 \[Mu]ijL[[n,n]])]),{a,1,nv},{l,1,nSH},{n,1,nSH},{m,1,nSH}],{i,1,nSL},{j,1,nSL}]];
Contri3=1/(16 \[Pi]^2)*(1/4)Simplify[Table[Sum[HabijVL[[a,b,i,j]]gAvss[[a,n,l]]gAvss[[b,l,n]](-Log[\[Mu]3US/(2 \[Mu]ijL[[n,n]])]),{a,1,nv},{b,1,nv},{l,1,nSH},{n,1,nSH}],{i,1,nSL},{j,1,nSL}]];
Contri4=1/4*(1/(16 \[Pi]^2))Simplify[Table[Sum[\[Lambda]K[[a,b,i,j]]\[Lambda]4K[[a,b,c,c]]\[Mu]ijL[[c,c]]/(\[Mu]ijL[[a,a]]+\[Mu]ijL[[b,b]]),{a,1,nSH},{b,1,nSH},{c,1,nSH}],{i,1,nSL},{j,1,nSL}]];
*)

If[Length[\[Lambda]x//Normal//Variables]==0&&Length[\[Lambda]y//Normal//Variables]==0,
Contri5=0;
Contri6=0;
Contri7=0;
,
TensHelp=Table[(1/2+Log[\[Mu]3/(a+b+n)]) ,{a,MassVector},{b,MassVector},{n,MassVector}]//SparseArray;
TensHelp2=TensorProduct[\[Lambda]x,\[Lambda]x]//DiagonalTensor2[#,2,6]&//DiagonalTensor2[#,3,6]&//DiagonalTensor2[#,4,6]&;
Contri5=-1/(16 \[Pi]^2)/3!*Total[TensHelp TensHelp2,-3];

TensHelp=Table[(1/2+Log[\[Mu]3/(a)]) ,{a,MassVector}]//SparseArray;
TensHelp2=TensorProduct[\[Lambda]y,\[Lambda]y]//TensorContract[#,{2,6}]&//TensorContract[#,{2,5}]&//DiagonalTensor2[#,2,4]&;
Contri6=-1/(16 \[Pi]^2)/2!*Total[TensHelp TensHelp2,-3];

TensHelp=Table[c/b ,{c,MassVector},{b,MassVector}]//SparseArray;
TensHelp2=Transpose[\[Lambda]y,{1,2,4,3}] . \[Lambda]x//DiagonalTensor2[#,3,4]&//DiagonalTensor2[#,4,5]&;
Contri7=(1/(16 \[Pi]^2))/2*Total[TensHelp TensHelp2,-3];

];
TensHelp=Table[n/m^2*l ,{n,MassVector},{m,MassVector},{l,MassVector}]//SparseArray;
TensHelp2=TensorProduct[DiagonalTensor2[\[Lambda]x,2,3],DiagonalTensor2[\[Lambda]x,2,3]]//DiagonalTensor2[#,3,6]&//Flatten[#,{{2},{1},{4}}]&;
ContriMix1=-1/4*1/(4 \[Pi])^2*Total[TensHelp TensHelp2,-3];

ContriMix2=-Simplify[Table[Sum[(\[Mu]ijMix[[i,m]])\[Mu]ijL[[m,m]]^-2 \[Mu]ijMix[[j,m]],{m,1,nSH}],{i,1,nSL},{j,1,nSL}]];

If[Length[\[Lambda]3//Normal//Variables]||nv>=nSH,
ContriC1=0;
ContriC2=0;
ContriC3=0;
ContriC4=0;
\[Mu]ijSSNLO2=0;
,
(*Check signs*)
MassVector=Table[\[Mu]ijLS[[n,n]],{n,1,ns}];
TensHelp=Table[MassHelp1[j,i,k,l,m]  ,{j,MassVector},{i,MassVector},{k,MassVector},{l,MassVector},{m,MassVector}]//SparseArray//SimplifySparse&;
TensHelp2=TensorProduct[\[Lambda]3CTot,\[Lambda]3CTot]//DiagonalTensor2[#,1,4]&//TensorProduct[#,\[Lambda]3CTot]&//DiagonalTensor2[#,2,6]&//DiagonalTensor2[#,4,6]&;
TensHelp3=TensorProduct[TensHelp2,\[Lambda]3CTot]//DiagonalTensor2[#,4,7]&//DiagonalTensor2[#,5,7]&;
ContriC1=-1/2*Total[TensHelp TensHelp3,-3][[LightScalar[[;;,1]],LightScalar[[;;,1]]]]//SimplifySparse;

TensHelp=Table[MassHelp2[i,j,m,n]  ,{i,MassVector},{j,MassVector},{m,MassVector},{n,MassVector}]//SparseArray//SimplifySparse;
TensHelp2=TensorProduct[\[Lambda]4Tot,\[Lambda]3CTot]//DiagonalTensor2[#,4,6]&//DiagonalTensor2[#,4,5]&;
TensHelp3=TensorProduct[TensHelp2,\[Lambda]3CTot]//DiagonalTensor2[#,4,7]&//DiagonalTensor2[#,5,7]&;
ContriC2=1/2*Total[TensHelp TensHelp3,-3][[LightScalar[[;;,1]],LightScalar[[;;,1]]]]//SimplifySparse;

TensHelp=Table[MassHelp2[i,j,m,n]  ,{i,MassVector},{j,MassVector},{m,MassVector},{n,MassVector}]//SparseArray//SimplifySparse;
TensHelp2=TensorProduct[\[Lambda]3CTot,\[Lambda]3CTot]//DiagonalTensor2[#,3,4]&//DiagonalTensor2[#,3,4]&;
TensHelp3=TensorProduct[TensHelp2,\[Lambda]4Tot]//DiagonalTensor2[#,3,5]&//DiagonalTensor2[#,4,7]&;
ContriC4=-1/4*Total[TensHelp TensHelp3,-3][[LightScalar[[;;,1]],LightScalar[[;;,1]]]]//SimplifySparse;

\[Mu]ijSSNLO2=(ContriC1+2*ContriC2+ContriC4)//Simplify//SparseArray;
];

\[Mu]ijTemp=\[Mu]ijLight+\[Mu]ijSSLO//SparseArray;
ContriF=-ZSij . \[Mu]ijTemp-\[Mu]ijTemp . ZSij;

\[Mu]ijSSNLO=(Contri1+Contri2+Contri3+ Contri4+Contri5+Contri6+Contri7+ContriMix1+ContriMix2+ContriF)//Simplify//SparseArray;

];



(*
	These functions appear in the matching.
*)
MassHelp1[0,0,0,0,0]:=0
MassHelp1[x_,0,0,0,0]:=1/(16 \[Pi]^2 x^4)
MassHelp1[0,x_,0,0,0]:=1/(16 \[Pi]^2 x^4)
MassHelp1[0,0,x_,0,0]:=1/(16 \[Pi]^2 x^4)
MassHelp1[0,0,0,x_,0]:=1/(16 \[Pi]^2 x^4)
MassHelp1[y_,y_,w_,z_,t_]:=1/(16 \[Pi]^2 (t^2+w^2+y^2) (t^2+y^2+z^2))
MassHelp1[x_,y_,z_,z_,t_]:=1/(16 \[Pi]^2 (t^2+x^2+z^2) (t^2+y^2+z^2))
MassHelp1[0,y_,0,0,t_]:=1/(16 \[Pi]^2 (t^2) (t^2+y^2))
MassHelp1[x_,y_,w_,z_,t_]:=(Log[\[Mu]3/(t+w+x)]-Log[\[Mu]3/(t+w+y)]-Log[\[Mu]3/(t+x+z)]+Log[\[Mu]3/(t+y+z)])/(16 \[Pi]^2 (w^2-z^2) (x^2-y^2))
MassHelp1[0,y_,0,z_,0]:=-((2 Log[\[Mu]3/y]-2 Log[\[Mu]3/(y+z)]+2 Log[\[Mu]3/z]+1)/(32 \[Pi]^2 y^2 z^2))
MassHelp1[x_,0,w_,0,0]:=-((2 Log[\[Mu]3/x]-2 Log[\[Mu]3/(x+w)]+2 Log[\[Mu]3/w]+1)/(32 \[Pi]^2 x^2 w^2))
MassHelp1[0,y_,w_,0,0]:=-((2 Log[\[Mu]3/y]-2 Log[\[Mu]3/(y+w)]+2 Log[\[Mu]3/w]+1)/(32 \[Pi]^2 y^2 w^2))
MassHelp1[x_,0,0,z_,0]:=-((2 Log[\[Mu]3/x]-2 Log[\[Mu]3/(x+z)]+2 Log[\[Mu]3/z]+1)/(32 \[Pi]^2 x^2 z^2))


(*
	These functions appear in the matching.
*)
MassHelpTriangle[0,0,0]:=0
MassHelpTriangle[x_,0,0]:=-(1/(4 \[Pi] x^3))
MassHelpTriangle[0,y_,0]:=-(1/(4 \[Pi] y^3))
MassHelpTriangle[0,0,z_]:=-(1/(4 \[Pi] z^3))
MassHelpTriangle[x_,y_,z_]:=1/(4 \[Pi] (x+y) (x+z) (y+z))


(*
	These functions appear in the matching.
*)
MassHelp2[0,0,0,0]:=0
MassHelp2[0,0,z_,t_]:=1/(16 \[Pi]^2 (t^2+z^2))
MassHelp2[x_,0,0,0]:=(-(2 Log[\[Mu]3/x])-1)/(32 \[Pi]^2 x^2)
MassHelp2[0,y_,0,0]:=(-(2 Log[\[Mu]3/y])-1)/(32 \[Pi]^2 y^2)
MassHelp2[y_,y_,z_,t_]:=1/(16 \[Pi]^2 (t^2+y^2+z^2))
MassHelp2[x_,y_,z_,t_]:=(Log[\[Mu]3/(t+y+z)]-Log[\[Mu]3/(t+x+z)])/(16 \[Pi]^2 (x^2-y^2))


(*
	Adds tree-level contributions before the matching.
*)
IdentifyTensorsPreSSDRalgo[]:=Module[{},


(*Quartic Tensor*)
HelpList=DeleteDuplicates@SparseArray[Flatten@Simplify[ (\[Lambda]4S)]]//Sort//FullSimplify;
HelpVar=Table[\[Lambda]PSS[a],{a,1,Delete[HelpList,1]//Length}];
HelpVarMod=RelationsBVariables[HelpList,HelpVar];
 HelpSolveQuarticPreS=Table[{Delete[HelpList,1][[a]]->HelpVarMod[[a]]},{a,1,Delete[HelpList,1]//Length}]//Flatten;
\[Lambda]4S=(\[Lambda]4S)//Normal//Simplify//FullSimplify//ReplaceAll[#,HelpSolveQuarticPreS]&//SparseArray;


HelpList=DeleteDuplicates@SparseArray[Flatten@Simplify[ (\[Lambda]K)]]//Sort//FullSimplify;
HelpVar=Table[\[Lambda]PSK[a],{a,1,Delete[HelpList,1]//Length}];
HelpVarMod=RelationsBVariables[HelpList,HelpVar];
 HelpSolveQuarticPreK=Table[{Delete[HelpList,1][[a]]->HelpVarMod[[a]]},{a,1,Delete[HelpList,1]//Length}]//Flatten;
\[Lambda]K=(\[Lambda]K)//Normal//Simplify//FullSimplify//ReplaceAll[#,HelpSolveQuarticPreK]&//SparseArray;


HelpList=DeleteDuplicates@SparseArray[Flatten@Simplify[ (\[Lambda]y)]]//Sort//FullSimplify;
HelpVar=Table[\[Lambda]PSY[a],{a,1,Delete[HelpList,1]//Length}];
HelpVarMod=RelationsBVariables[HelpList,HelpVar];
 HelpSolveQuarticPreY=Table[{Delete[HelpList,1][[a]]->HelpVarMod[[a]]},{a,1,Delete[HelpList,1]//Length}]//Flatten;
\[Lambda]y=(\[Lambda]y)//Normal//Simplify//FullSimplify//ReplaceAll[#,HelpSolveQuarticPreY]&//SparseArray;


HelpList=DeleteDuplicates@SparseArray[Flatten@Simplify[ \[Lambda]x]]//Sort//FullSimplify;
HelpVar=Table[\[Lambda]PSX[a],{a,1,Delete[HelpList,1]//Length}];
HelpVarMod=RelationsBVariables[HelpList,HelpVar];
 HelpSolveQuarticPreX=Table[{Delete[HelpList,1][[a]]->HelpVarMod[[a]]},{a,1,Delete[HelpList,1]//Length}]//Flatten;
\[Lambda]x=(\[Lambda]x)//Normal//Simplify//FullSimplify//ReplaceAll[#,HelpSolveQuarticPreX]&//SparseArray;


HelpList=DeleteDuplicates@SparseArray[Flatten@Simplify[ \[Lambda]4Tot]]//Sort//FullSimplify;
HelpVar=Table[\[Lambda]PST[a],{a,1,Delete[HelpList,1]//Length}];
HelpVarMod=RelationsBVariables[HelpList,HelpVar];
 HelpSolveQuarticPreTot=Table[{Delete[HelpList,1][[a]]->HelpVarMod[[a]]},{a,1,Delete[HelpList,1]//Length}]//Flatten;
\[Lambda]4Tot=(\[Lambda]4Tot)//Normal//Simplify//FullSimplify//ReplaceAll[#,HelpSolveQuarticPreTot]&//SparseArray;

IdentMatPre=List/@Join[HelpSolveQuarticPreS,HelpSolveQuarticPreK,HelpSolveQuarticPreY,HelpSolveQuarticPreX,HelpSolveQuarticPreTot]/.{b_->a_}:>a->b//Flatten[#,1]&;

];


(*
	Writes ultrasoft parameters in terms of soft ones.
*)
IdentifyTensorsSSDRalgo[]:=Module[{},

If[verbose,Print["Identifying Components"]];
(*Quartic Tensor*)


HelpList=DeleteDuplicates@SparseArray[Simplify@Normal@Flatten[ (\[Lambda]4S+\[Lambda]3DSS)]]//Sort;
If[HelpList[[1]]==0&&Length[HelpList]>1,
	HelpList=Delete[HelpList,1];
];
HelpVar=Table[\[Lambda]SS[a],{a,1,HelpList//Length}];
HelpVarMod=RelationsBVariables[HelpList,HelpVar];
HelpSolveQuarticS=Table[{HelpList[[a]]->HelpVarMod[[a]]},{a,1,HelpList//Length}]//Flatten;
\[Lambda]3DSS=(\[Lambda]4S+\[Lambda]3DSS)//Normal//Simplify//ReplaceAll[#,HelpSolveQuarticS]&//SparseArray;



(*Cubic Tensor*)
HelpList=DeleteDuplicates@SparseArray[Flatten@Simplify[(\[Lambda]3CSSS+\[Lambda]3CLight)]]//Sort//FullSimplify;
If[HelpList[[1]]==0&&Length[HelpList]>1,
	HelpList=Delete[HelpList,1];
];
HelpVar=Table[cSSSS[a],{a,1,HelpList//Length}];
HelpVarMod=RelationsBVariables[HelpList,HelpVar];
HelpSolveCubicSSS=Table[{HelpList[[a]]->HelpVarMod[[a]]},{a,1,HelpList//Length}]//Flatten//Simplify;
\[Lambda]3CSRedSS=(\[Lambda]3CSSS+\[Lambda]3CLight)//Normal//Simplify//FullSimplify//ReplaceAll[#,HelpSolveCubicSSS]&//SparseArray;
(*Scalar-Vector Tensor*)

HelpList=DeleteDuplicates@Simplify@Flatten[(HabijVL+GvvssTSS)]//Sort;
If[HelpList[[1]]==0&&Length[HelpList]>1,
	HelpList=Delete[HelpList,1];
];
HelpVar=Table[ \[Lambda]VTSS[a],{a,1,HelpList//Length}];
HelpVarMod=RelationsBVariables[HelpList,HelpVar];
HelpSolveVecTS=Table[{HelpList[[a]]->HelpVarMod[[a]]},{a,1,HelpList//Length}]//Flatten;
\[Lambda]KVecTSS= (HabijVL+GvvssTSS)//Normal//Simplify//ReplaceAll[#,HelpSolveVecTS]&//SparseArray;


(*Scalar Mass*)
If[mode>=2,
HelpList=DeleteDuplicates@Flatten@Simplify[ xLO \[Mu]ijLight+xLO \[Mu]ijSSLO+xNLO(\[Mu]ijSSNLO+\[Mu]ijSSNLO2)]//Sort;
If[HelpList[[1]]==0&&Length[HelpList]>1,
	HelpList=Delete[HelpList,1];
];
HelpVar=Table[ \[Mu]ijSS[a],{a,1,HelpList//Length}];
HelpVarMod=RelationsBVariables[HelpList,HelpVar];
HelpSolveMassS=Table[{HelpList[[a]]->HelpVarMod[[a]]},{a,1,HelpList//Length}]//Simplify//Flatten;
\[Mu]ijSNLOSS=xLO \[Mu]ijLight+xLO \[Mu]ijSSLO+xNLO(\[Mu]ijSSNLO+\[Mu]ijSSNLO2)//Normal//Simplify//ReplaceAll[#,HelpSolveMassS]&//SparseArray;
,
HelpList=DeleteDuplicates@Flatten@Simplify[ xLO \[Mu]ijLight+xLO \[Mu]ijSSLO]//Sort;
If[HelpList[[1]]==0&&Length[HelpList]>1,
	HelpList=Delete[HelpList,1];
];
HelpVar=Table[ \[Mu]ijSS[a],{a,1,HelpList//Length}];
HelpVarMod=RelationsBVariables[HelpList,HelpVar];
HelpSolveMassS=Table[{HelpList[[a]]->HelpVarMod[[a]]},{a,1,HelpList//Length}]//Simplify//Flatten;
\[Mu]ijSNLOSS=xLO \[Mu]ijLight+xLO \[Mu]ijSSLO//Normal//Simplify//ReplaceAll[#,HelpSolveMassS]&//SparseArray;
];

(*Scalar Tadpoles*)
HelpList=DeleteDuplicates@Flatten@Simplify[ TadPoleLightSS]//Sort;
If[HelpList[[1]]==0&&Length[HelpList]>1,
	HelpList=Delete[HelpList,1];
];
HelpVar=Table[ dSS[a],{a,1,HelpList//Length}];
HelpVarMod=RelationsBVariables[HelpList,HelpVar];
HelpSolveTadpoleSS=Table[{HelpList[[a]]->HelpVarMod[[a]]},{a,1,HelpList//Length}]//Flatten;
TadPoleSSSLO=TadPoleLightSS//Normal//Simplify//ReplaceAll[#,HelpSolveTadpoleSS]&//SparseArray;



IdentMatSS=List/@Join[HelpSolveQuarticS,HelpSolveVecTS,HelpSolveMassS,HelpSolveCubicSSS,HelpSolveTadpoleSS]/.{b_->a_}:>a->b//Flatten[#,1]&;

];


{TadPoleSSSLO};


(*
	Prints tadpoles.
*)
PrintTadpolesUS[optP_]:=Module[{opt=optP},
If[verbose,Print["Printing Scalar Masses"]];

VarGauge=Join[TadPoleLight//Normal//Variables]//DeleteDuplicates;
SubGauge=Table[c->Symbol[ToString[c]<>ToString["US"]],{c,VarGauge}];


\[Lambda]1p=TadPoleLight//Normal//ReplaceAll[#,SubGauge]&;
var=Normal[\[Lambda]1p]//Variables;
helpMass=Normal[\[Lambda]1p-TadPoleSSSLO];
ResScalp=Reduce[helpMass==0,var]//ToRules[#]&;
SolveTemp=var/.ResScalp;
SolTadpole=Table[{var[[i]]->SolveTemp[[i]]},{i,1,Length@var}]//Flatten[#,1]&//ReplaceAll[#,IdentMatSS]&;

ToExpression[StringReplace[ToString[StandardForm[Join[SolTadpole]]],"DRalgo`Private`"->""]]

];


(*
	Prints scalar masses.
*)
PrintScalarMassUS[optP_]:=Module[{opt=optP},
If[verbose,Print["Printing Scalar Masses"]];

VarGauge=Join[\[Mu]ijLight//Normal//Variables]//DeleteDuplicates;
SubGauge=Table[c->Symbol[ToString[c]<>ToString["US"]],{c,VarGauge}];

\[Mu]ijp=\[Mu]ijLight//Normal//ReplaceAll[#,SubGauge]&;
var=Normal[\[Mu]ijp]//Variables;
helpMass=Normal[\[Mu]ijp-\[Mu]ijSNLOSS];
ResScalp=Reduce[helpMass==0,var]//ToRules[#]&;
SolveTemp=var/.ResScalp;
SolMassPre=Table[{var[[i]]->SolveTemp[[i]]},{i,1,Length@var}]//Flatten[#,1]&//ReplaceAll[#,IdentMatSS]&;


If[opt=="LO",
SolMass=SolMassPre/.xLO->1/.xNLO->0;
,
SolMass=SolMassPre/.xLO->0/.xNLO->1;
];
(*Printing Result*)
ToExpression[StringReplace[ToString[StandardForm[Join[SolMass]]],"DRalgo`Private`"->""]]

];


PrintIdentificationUS[]:=Module[{},
ToExpression[StringReplace[ToString[StandardForm[Join[IdentMatSS]]],"DRalgo`Private`"->""]]
];


PrintTensorUSDRalgo[]:=Module[{},
Print["Order of Tensors: (1) Scalar Quartic, (2) Vector-Scalar Couplings, (3) Scalar Mass"];
ToExpression[StringReplace[ToString[StandardForm[Join[\[Lambda]3DSS,\[Lambda]KVecTSS,\[Mu]ijSNLOSS]]],"DRalgo`Private`"->""]]
];


(*
	Calculates the 1-loop scalar mass.
*)
ScalarMassSS[]:=Module[{},
If[verbose,Print["Calculating 1-Loop Scalar Mass"]];

SelfEnergySS=1/(4\[Pi]);
ContriSS=SelfEnergySS/2*Simplify[Table[Sum[\[Lambda]3Cx[[i,ii,jj]]\[Lambda]3Cx[[j,ii,jj]]/(\[Mu]ijL[[ii,ii]]+\[Mu]ijL[[jj,jj]]),{ii,1,nSH},{jj,1,nSH}],{i,1,nSL},{j,1,nSL}]];
ContriSS2=SelfEnergySS*Simplify[Table[Sum[\[Lambda]3Cy[[i,ii,jj]]\[Lambda]3Cy[[j,ii,jj]]/(\[Mu]ijL[[jj,jj]]),{ii,1,nSL},{jj,1,nSH}],{i,1,nSL},{j,1,nSL}]];

ContriTadpole=Table[Sum[\[Lambda]3Cy[[i,j,ll]]TadPoleHeavy[[ll]]/(\[Mu]ijL[[ll,ll]]^2),{ll,1,nSH}],{i,1,nSL},{j,1,nSL}];

ContriSS3=1/(4 \[Pi])/2 Simplify[Table[Sum[ \[Mu]ijL[[a,a]]\[Lambda]K[[a,a,i,j]],{a,1,nSH}],{i,1,nSL},{j,1,nSL}]];


\[Mu]ijSSLO=-ContriSS3-ContriSS-ContriSS2-ContriTadpole//SparseArray;

];


(*
	Adds tree-level Corrections to Scalar-quartic couplings.
*)
ScalarModifiedSS[]:=Module[{},
If[verbose,Print["Adding tree-level Corrections to Scalar Quartics"]];

ContriTLTemp=-Simplify[Table[Sum[1/(\[Mu]ijL[[n,n]]^2)\[Lambda]3Cy[[i,j,n]]\[Lambda]3Cy[[k,l,n]],{n,1,nSH}],{i,1,nSL},{j,1,nSL},{k,1,nSL},{l,1,nSL}]];
ContriTL=Table[ContriTLTemp[[i,j,k,l]]+ContriTLTemp[[i,k,j,l]]+ContriTLTemp[[i,l,j,k]],{i,1,nSL},{j,1,nSL},{k,1,nSL},{l,1,nSL}];
\[Lambda]Mod=\[Lambda]4S+ContriTL;
\[Lambda]4S=\[Lambda]Mod;

ContriTLTemp=-Simplify[Table[Sum[1/(\[Mu]ijL[[n,n]]^2)\[Lambda]3CHeavy[[ii,jj,n]]\[Lambda]3Cy[[k,l,n]],{n,1,nSH}],{ii,1,nSH},{jj,1,nSH},{k,1,nSL},{l,1,nSL}]];
ContriTLTemp2=-Simplify[Table[Sum[1/(\[Mu]ijL[[n,n]]^2)\[Lambda]3Cx[[k,ii,n]]\[Lambda]3Cx[[l,jj,n]],{n,1,nSH}],{ii,1,nSH},{jj,1,nSH},{k,1,nSL},{l,1,nSL}]];
ContriTL=Table[ContriTLTemp[[ii,jj,k,l]]+ContriTLTemp2[[ii,jj,k,l]]+ContriTLTemp2[[jj,ii,k,l]],{ii,1,nSH},{jj,1,nSH},{k,1,nSL},{l,1,nSL}];
\[Lambda]Mod=\[Lambda]K+ContriTL;
\[Lambda]K=\[Lambda]Mod;

ContriTLTemp=-Simplify[Table[Sum[1/(\[Mu]ijL[[n,n]]^2)\[Lambda]3Cx[[i,jj,n]]\[Lambda]3Cy[[k,l,n]],{n,1,nSH}],{i,1,nSL},{k,1,nSL},{l,1,nSL},{jj,1,nSH}]];
ContriTL=Table[ContriTLTemp[[i,k,l,jj]]+ContriTLTemp[[k,i,l,jj]]+ContriTLTemp[[l,k,i,jj]],{i,1,nSL},{k,1,nSL},{l,1,nSL},{jj,1,nSH}];
\[Lambda]Mod=\[Lambda]y+ContriTL;
\[Lambda]y=\[Lambda]Mod;


ContriTLTemp=-Simplify[Table[Sum[1/(\[Mu]ijL[[n,n]]^2)\[Lambda]3CHeavy[[ii,jj,n]]\[Lambda]3Cx[[k,ll,n]],{n,1,nSH}],{k,1,nSL},{ii,1,nSH},{jj,1,nSH},{ll,1,nSH}]];
ContriTL=Table[ContriTLTemp[[k,ii,jj,ll]]+ContriTLTemp[[k,ii,ll,jj]]+ContriTLTemp[[k,ll,jj,ii]],{k,1,nSL},{ii,1,nSH},{jj,1,nSH},{ll,1,nSH}];
\[Lambda]Mod=\[Lambda]x+ContriTL;
\[Lambda]x=\[Lambda]Mod;

MassHelp[0]:=0;
MassHelp[x_]:=1/(x);

ContriTLTemp=-Simplify[Table[Sum[MassHelp[\[Mu]ijLS[[n,n]]]^2 \[Lambda]3CTot[[i,j,n]]\[Lambda]3CTot[[k,l,n]],{n,1,ns}],{i,1,ns},{j,1,ns},{k,1,ns},{l,1,ns}]];
ContriTL=Table[ContriTLTemp[[i,j,k,l]]+ContriTLTemp[[i,k,j,l]]+ContriTLTemp[[i,l,k,j]],{i,1,ns},{j,1,ns},{k,1,ns},{l,1,ns}];
\[Lambda]4Tot=\[Lambda]3DSp+ContriTL;


\[Lambda]4Tot[[HeavyScalars[[;;,1]],HeavyScalars[[;;,1]],HeavyScalars[[;;,1]],HeavyScalars[[;;,1]]]]=\[Lambda]3DSp[[HeavyScalars[[;;,1]],HeavyScalars[[;;,1]],HeavyScalars[[;;,1]],HeavyScalars[[;;,1]]]];


];


(*
	Calculates the scalar quartic in the ultrasoft theory.
*)
ScalarQuarticSS[]:=Module[{},
If[verbose,Print["Calculating Scalar Quartics"]];



(*Maybe include Cubics later*)
(*Temp=0*Simplify[Table[Sum[(1/(\[Mu]ijL[[n,n]]^2))\[Lambda]3Cy[[i,j,n]]\[Lambda]3CHeavy[[k,l,n]],{n,1,nSH}],{k,1,nSH},{l,1,nSH},{i,1,nSL},{j,1,nSL}]];*)
\[Lambda]KEff=\[Lambda]K//SparseArray;
(*Temp=0*Simplify[Table[Sum[(1/(\[Mu]ijL[[n,n]]^2))\[Lambda]3Cy[[i,j,n]]\[Lambda]3Cx[[k,l,n]],{n,1,nSH}],{i,1,nSL},{j,1,nSL},{k,1,nSL},{l,1,nSH}]];*)
\[Lambda]yEff=\[Lambda]y//SparseArray;
(*\[Lambda]yEff2=Table[\[Lambda]y[[i,j,k,l]]-0*Temp[[i,j,k,l]]-0*Temp[[i,k,j,l]]-0*Temp[[k,j,i,l]],{i,1,nSL},{j,1,nSL},{k,1,nSL},{l,1,nSH}]//SparseArray;*)
\[Lambda]yEff2=\[Lambda]y;

(*Loop level*)
MassVec=Table[\[Mu]ijL[[n,n]],{n,1,nSH}];
TensHelp=Table[1/(n+m),{n,MassVec},{m,MassVec}]//SparseArray//Flatten[#,{1,2}]&;
AE2=Flatten[\[Lambda]KEff,{{1,2}}]//SparseArray;
Temp=1/2*1/(4 \[Pi])*Transpose[AE2,{3,2,1}] . (TensHelp\[NonBreakingSpace]AE2)//SparseArray;
ContriSS=3Symmetrize[Temp,Symmetric[{1,2,3,4}]]//SparseArray//SimplifySparse;


If[Length[\[Lambda]y//Normal//Variables]==0,
ContriSS2=EmptyArray[{nSL,nSL,nSL,nSL}];
ContriSS3=EmptyArray[{nSL,nSL,nSL,nSL}];
,
MassVec=Table[\[Mu]ijL[[n,n]],{n,1,nSH}];
TensHelp=Table[1/(n),{n,MassVec}]//SparseArray;
AE2=TensHelp Transpose[\[Lambda]yEff,{4,2,3,1 }]//Flatten[#,{2,1}]& ;
Temp=1/(4 \[Pi]) Transpose[Flatten[\[Lambda]yEff,{3,4}],{3,2,1}] . AE2;
ContriSS2=3Symmetrize[Temp,Symmetric[{1,2,3,4}]]//SparseArray//SimplifySparse;

TensHelp=Table[n/m^2,{m,MassVec},{n,MassVec}]//SparseArray;
TensHelp2=Transpose[\[Lambda]x,{4,2,3,1}]//Table[i,{i,#}]&//Table[#[[i,i]],{i,1,Length[#]}]&//SparseArray;
TensHelp3=TensHelp . TensHelp2//Table[i,{i,#}]&//Table[#[[i,i]],{i,1,Length[#]}]&//SparseArray;
ContriSS3=-1/(4 \[Pi])*1/2*4*Symmetrize[\[Lambda]yEff2 . TensHelp3,Symmetric[{1,2,3,4}]]//SparseArray//SimplifySparse;

];

	
If[Length[\[Lambda]3//Normal//Variables]==0,
\[Lambda]3DSS=-ContriSS-ContriSS2-ContriSS3;
,
ContriSETemp=ZSij . \[Lambda]4S//SparseArray;
ContriSE=ContriSETemp+Transpose[ContriSETemp,{2,1,3,4}]+Transpose[ContriSETemp,{3,1,2,4}]+Transpose[ContriSETemp,{4,1,2,3}]//SparseArray//SimplifySparse;

TensHelp=Table[1/(n^2),{n,MassVec}]//SparseArray;
TensHelp2=TensorProduct[\[Lambda]3Cy,\[Lambda]3Cy]//DiagonalTensor2[#,3,6]&;
ContriTLTemp1=TensHelp . TensHelp2;
TensHelp=Table[1/(n^2*m^2),{n,MassVec},{m,MassVec}]//SparseArray;
TensHelp=\[Mu]HEff TensHelp;
ContriTLTemp2=\[Lambda]3Cy . TensHelp . Transpose[\[Lambda]3Cy,{3,2,1}];
ContriTL=3 Symmetrize[ContriTLTemp1+ContriTLTemp2,Symmetric[{1,2,3,4}]]//SparseArray//SimplifySparse;

TensHelp=Table[1/(n^2*m^2),{n,MassVec},{m,MassVec}]//SparseArray;
TensHelp2=TensorProduct[\[Lambda]yEff2,\[Lambda]3Cx,TadPoleHeavySS]//DiagonalTensor2[#,4,6]&//DiagonalTensor2[#,6,7]&;
Temp=-TensHelp . TensHelp2//TensorContract[#,{1,2}]&;
CouplingSSTadpole=4 Symmetrize[Temp,Symmetric[{1,2,3,4}]]//SparseArray//SimplifySparse;


(*TadPole*)
TensHelp=Table[1/(n^2*m^2),{n,MassVec},{m,MassVec}]//SparseArray;
TensHelp2=TensorProduct[\[Lambda]K,\[Lambda]3Cy,TadPoleHeavySS]//DiagonalTensor2[#,1,7]&//DiagonalTensor2[#,2,7]&;
Temp=-TensHelp . TensHelp2//TensorContract[#,{1,2}]&;
ContriTadpole=2*3 Symmetrize[Temp,Symmetric[{1,2,3,4}]]//SparseArray//SimplifySparse;

Clear[MassHelp,MassHelp5,MassHelp2,MassHelp3];
MassHelp5[x_]:=-(1/(4 \[Pi] x^3));
MassHelp4[x_,y_,z_]:=1/(4 \[Pi] (x+y) (x+z) (y+z));
MassHelp[x_,y_,w_,z_]:=(w+x+y+z)/(4 \[Pi] (w+x) (w+y) (w+z) (x+y) (x+z) (y+z));
MassHelp2[w_,z_]:=-((w^2+w z+z^2)/(4 \[Pi] w^3 z^3 (w+z)));
MassHelp3[w_]:=1/(4 \[Pi] w^5);

If[nSH>nv,
MassVec2=Table[\[Mu]ijL[[n,n]],{n,nv+1,nSH}];
TensHelp=Table[MassHelp[n,m,l,o],{l,MassVec2},{m,MassVec2},{n,MassVec2},{o,MassVec2}]//SparseArray;
TensHelp=Delete[TensHelp//ArrayRules,-1]//ReplaceAll[#,{x_,y_,z_,w_}->{x+nv,y+nv,z+nv,w+nv}]&//SparseArray;
TensHelp2=TensorProduct[\[Lambda]3Cx,\[Lambda]3Cx]//DiagonalTensor2[#,3,5]&//TensorProduct[#,\[Lambda]3Cx]&//DiagonalTensor2[#,5,7]&//TensorProduct[#,\[Lambda]3Cx]&//DiagonalTensor2[#,7,10]&//DiagonalTensor2[#,5,9]&;
ContriSSSSTemp=Total[TensHelp TensHelp2,-5];
ContriCubic1=2*3 Symmetrize[ContriSSSSTemp,Symmetric[{1,2,3,4}]]//SparseArray//SimplifySparse;

TensHelp=Table[MassHelp[n,0,l,o],{l,MassVec2},{n,MassVec2},{o,MassVec2}]//SparseArray;
TensHelp=Delete[TensHelp//ArrayRules,-1]//ReplaceAll[#,{x_,y_,z_}->{x+nv,y+nv,z+nv}]&//SparseArray;
TensHelp2=TensorProduct[\[Lambda]3Cy,\[Lambda]3Cy]//TensorContract[#,{2,5}]&//TensorProduct[#,\[Lambda]3Cx]&//DiagonalTensor2[#,4,6]&//TensorProduct[#,\[Lambda]3Cx]&//DiagonalTensor2[#,6,9]&//DiagonalTensor2[#,4,8]&;
ContriSSSSTemp=Total[TensHelp TensHelp2,-5];
ContriCubic2=12 Symmetrize[ContriSSSSTemp,Symmetric[{1,2,3,4}]]//SparseArray//SimplifySparse;
 
TensHelp=Table[MassHelp2[n,o],{n,MassVec2},{o,MassVec2}]//SparseArray;
TensHelp=Delete[TensHelp//ArrayRules,-1]//ReplaceAll[#,{x_,y_}->{x+nv,y+nv}]&//SparseArray;
TensHelp2=Transpose[\[Lambda]3Cy,{3,2,1}] . \[Lambda]3CLight . \[Lambda]3Cy//TensorProduct[#,\[Lambda]3Cx]&//DiagonalTensor2[#,1,8]&//DiagonalTensor2[#,5,7]&;
ContriSSSSTemp=Total[TensHelp TensHelp2,-5];
ContriCubic3=4! Symmetrize[ContriSSSSTemp,Symmetric[{1,2,3,4}]]//SparseArray//SimplifySparse;

TensHelp=Table[MassHelp2[n,o],{n,MassVec2},{o,MassVec2}]//SparseArray;
TensHelp=Delete[TensHelp//ArrayRules,-1]//ReplaceAll[#,{x_,y_}->{x+nv,y+nv}]&//SparseArray;
TensHelp2=TensorProduct[Transpose[\[Lambda]3Cy,{3,2,1}] . \[Lambda]3Cy, Transpose[\[Lambda]3Cy,{3,2,1}] . \[Lambda]3Cy]//DiagonalTensor2[#,1,8]&//DiagonalTensor2[#,4,5]&;
ContriSSSSTemp=Total[TensHelp TensHelp2,-5];
ContriCubic4=3! Symmetrize[ContriSSSSTemp,Symmetric[{1,2,3,4}]]//SparseArray//SimplifySparse;

TensHelp=Table[MassHelp3[n],{n,MassVec2}]//SparseArray;
TensHelp=ReplaceAll[#,{y_}->{y+nv}]&/@Delete[TensHelp//ArrayRules,-1]//SparseArray;
TensHelp2=Transpose[\[Lambda]3Cy,{3,2,1}] . \[Lambda]3CLight . \[Lambda]3CLight . \[Lambda]3Cy//DiagonalTensor2[#,1,6]&;
ContriSSSSTemp=TensHelp . TensHelp2;
ContriCubic5=3! Symmetrize[ContriSSSSTemp,Symmetric[{1,2,3,4}]]//SparseArray//SimplifySparse;

TensHelp=Table[MassHelp4[n,m,l],{n,MassVec2},{m,MassVec2},{l,MassVec2}]//SparseArray;
TensHelp=Delete[TensHelp//ArrayRules,-1]//ReplaceAll[#,{x_,y_,z_}->{x+nv,y+nv,z+nv}]&//SparseArray;
TensHelp2=TensorProduct[\[Lambda]K,\[Lambda]3Cx]//DiagonalTensor2[#,1,6]&//TensorProduct[#,\[Lambda]3Cx]&//DiagonalTensor2[#,2,8]&//DiagonalTensor2[#,6,8]&;
ContriSSSSTemp=Total[TensHelp TensHelp2,-5];
ContriCubic6=-3! Symmetrize[ContriSSSSTemp,Symmetric[{1,2,3,4}]]//SparseArray//SimplifySparse;

TensHelp=Table[MassHelp4[n,m,0],{n,MassVec2},{m,MassVec2}]//SparseArray;
TensHelp=Delete[TensHelp//ArrayRules,-1]//ReplaceAll[#,{x_,y_}->{x+nv,y+nv}]&//SparseArray;
TensHelp2=Transpose[\[Lambda]3Cy,{3,2,1}] . \[Lambda]3Cy//TensorProduct[\[Lambda]K,#]&//DiagonalTensor2[#,1,5]&//DiagonalTensor2[#,2,7]&;
ContriSSSSTemp=Total[TensHelp TensHelp2,-5];
ContriCubic7=-3! Symmetrize[ContriSSSSTemp,Symmetric[{1,2,3,4}]]//SparseArray//SimplifySparse;

TensHelp=Table[MassHelp5[n],{n,MassVec2}]//SparseArray;
TensHelp=ReplaceAll[#,{y_}->{y+nv}]&/@Delete[TensHelp//ArrayRules,-1]//SparseArray;
TensHelp2=\[Lambda]3CLight . \[Lambda]3Cy//TensorProduct[\[Lambda]y,#]&//TensorContract[#,{3,5}]&//DiagonalTensor2[#,3,6]&;
ContriSSSSTemp=TensHelp . TensHelp2;
ContriCubic8=-3! Symmetrize[ContriSSSSTemp,Symmetric[{1,2,3,4}]]//SparseArray//SimplifySparse;

TensHelp=Table[MassHelp4[n,m,0],{n,MassVec2},{m,MassVec2}]//SparseArray;
TensHelp=Delete[TensHelp//ArrayRules,-1]//ReplaceAll[#,{x_,y_}->{x+nv,y+nv}]&//SparseArray;
TensHelp2=Transpose[\[Lambda]y,{1,2,4,3}] . \[Lambda]3Cy//TensorProduct[#,\[Lambda]3Cx]&//DiagonalTensor2[#,3,8]&//DiagonalTensor2[#,5,7]&;
ContriSSSSTemp=Total[TensHelp TensHelp2,-5];
ContriCubic9=-3! Symmetrize[ContriSSSSTemp,Symmetric[{1,2,3,4}]]//SparseArray//SimplifySparse;

(*
ContriSSSSTemp=-Table[Sum[\[Lambda]4S[[i,j,n,m]]\[Lambda]3Cy[[k,n,l]]\[Lambda]3Cx[[t,m,l]] MassHelp5[\[Mu]ijL[[l,l]]],{n,1,nSL},{m,1,nSL},{l,nv+1,nSH}],{i,1,nSL},{j,1,nSL},{k,1,nSL},{t,1,nSL}]//Simplify//SparseArray;
ContriCubic10=Table[(ContriSSSSTemp[[i,j,k,t]])+(ContriSSSSTemp[[i,k,j,t]])+(ContriSSSSTemp[[i,t,k,j]])+(ContriSSSSTemp[[k,t,i,j]])+(ContriSSSSTemp[[k,j,i,t]])+(ContriSSSSTemp[[j,t,i,k]]),{i,1,nSL},{j,1,nSL},{k,1,nSL},{t,1,nSL}];
*)
\[Lambda]3DSS=ContriSE-ContriTL-CouplingSSTadpole- ContriSS-ContriSS2-ContriSS3-ContriCubic1-ContriCubic2-ContriCubic3-ContriCubic4-ContriCubic5-ContriCubic6-ContriCubic7-ContriCubic8-ContriCubic9;

,

\[Lambda]3DSS=ContriSE-ContriTL-CouplingSSTadpole- ContriSS-ContriSS2-ContriSS3;
];

];

];


(*
	Calculates non-abelian couplings in the ultrasoft theory.
*)
NonAbelianCouplingSS[]:=Module[{},
ContriAnomVV= gvvvSS . Transpose[ZLij]//SparseArray//SimplifySparse;
(*Simplify[Table[Sum[ ZLij[[c,d]]gvvvSS[[a,b,d]],{d,1,nv}],{a,1,nv},{b,1,nv},{c,1,nv}]];*)
GgvvvSS=-ContriAnomVV;
];


(*
	Calculates vector couplings in the ultrasoft theory.
*)
ScalarVectorCouplingSS[]:=Module[{},
If[verbose,Print["Calculating Vector-Scalar Couplings"]];


(* Self-Energy contribution*)
ContriSEVector=-ZLij . HabijVL-Transpose[ZLij . HabijVL,{2,1,3,4}];
(*ContriSEVector=-Simplify[Table[Sum[ZLij[[a,c]](HabijVL[[c,b,i,j]])+ZLij[[b,c]](HabijVL[[a,c,i,j]]),{c,1,nv}],{a,1,nv},{b,1,nv},{i,1,nSL},{j,1,nSL}]];*)


GvvssTSS= ContriSEVector//Simplify;

];


(*
	Calculates the transverse-vector self energy.
*)
VectorSelfEnergySS[]:=Module[{},
If[verbose,Print["Calculating Vector Self-Energy"]];

SelfEnergySS=-1/2*(-1/(24 \[Pi]));
ContriSS=SelfEnergySS*Simplify[Table[Sum[HabijA[[a,b,i,i]]/\[Mu]ijL[[i,i]],{i,1,nSH}],{a,1,nv},{b,1,nv}]];
ZLij=-ContriSS/2;
];


{DiaRules,RotRules};(*Diagonalization matrix and masses*)


(*
	Prepares the reduction to the ultrasoft scale. Here scalars are divided into heavy and soft depending on the user's input.
	Heavy scalars are then grouped with temporal scalars and integrated out.
*)
PrepareSoftToSuperSoft[ListHardI_]:=Module[{ListP=ListHardI},

(*Fix for Non-diagonal Debye masses*)
(*************************)

(*
	Here we take the "Damn the torpedo, full speed ahead" approach.
*)

HeavyScalars=Transpose[List@ListP]; (*The list of all heavy scalars*)

(*
	Renames 4d couplings to 3d ones
*)

VarQuartic=Join[\[Lambda]4//Normal//Variables]//DeleteDuplicates;
SubQuartic=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarQuartic}];
\[Lambda]3DSp= \[Lambda]4//Normal//ReplaceAll[#,SubQuartic]&//SparseArray;

VarCubic=Join[\[Lambda]3//Normal//Variables]//DeleteDuplicates;
SubCubic=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarCubic}];
\[Lambda]3CSRedUS=\[Lambda]3//Normal//ReplaceAll[#,SubCubic]&//SparseArray;

VarMass=Join[\[Mu]ij//Normal//Variables]//DeleteDuplicates;
SubMass=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarMass}];
\[Mu]ijSNLOP=\[Mu]ij//Normal//ReplaceAll[#,SubMass]&//SparseArray;



If[Length[HeavyScalars]<1,HScal=False,HScal=True];

TotScalar=Table[{i},{i,1,ns}];
LightScalar=Delete[TotScalar,HeavyScalars];
nSH=nv+(HeavyScalars//Length); (*Adds heavy scalars to the temporal scalars*)
nSL=ns-(HeavyScalars//Length); (*Removes heavy scalars from the scalars*)

\[Mu]ijLS=Sqrt[\[Mu]ijSNLOP]//SparseArray; (*Rewrites squared masses as masses. Just a trick to make things easier.*)
\[Mu]ijLS[[LightScalar[[;;,1]],LightScalar[[;;,1]]]]=0;

VarGauge=GaugeCouplingNames;
SubGauge=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarGauge}];
gvvvSS=gvvv//Normal//ReplaceAll[#,SubGauge]&//SparseArray;


If[HScal,
(*Quartic Tensor*)
(*
	Here we create couplings between light and heavy scalarws.
*)
\[Lambda]4Light=Table[\[Lambda]3DSp[[a,b,c,d]],{a,LightScalar[[;;,1]]},{b,LightScalar[[;;,1]]},{c,LightScalar[[;;,1]]},{d,LightScalar[[;;,1]]}]//SparseArray; (*Light-scalar part of scalar quartic*)
\[Lambda]4Heavy=Table[\[Lambda]3DSp[[a,b,c,d]],{a,HeavyScalars[[;;,1]]},{b,HeavyScalars[[;;,1]]},{c,HeavyScalars[[;;,1]]},{d,HeavyScalars[[;;,1]]}]//SparseArray; (*Heavy-scalar part of scalar quartic*)
\[Lambda]KHeavy=Table[\[Lambda]KVec[[a,b,c,d]],{a,1,nv},{b,1,nv},{c,HeavyScalars[[;;,1]]},{d,HeavyScalars[[;;,1]]}]//SparseArray;  (*Heavy-scalar part of the scalar-temporalScal couplings*)


A1=Delete[\[Lambda]4Heavy//ArrayRules,-1]//ReplaceAll[#,{x_,y_,z_,w_}->{x+nv,y+nv,z+nv,w+nv}]&;
A2=ReplaceAll[Delete[\[Lambda]KHeavy//ArrayRules,-1],({x_,y_,z_,w_}->a_)->{x,y,z+nv,w+nv}->a];
A3=Delete[\[Lambda]AAS//ArrayRules,-1];
\[Lambda]KTotal=SymmetrizedArray[Join[A1,A2,A3],{nSH,nSH,nSH,nSH},Symmetric[{1,2,3,4}]]//SparseArray; (*Total heavy-scalar tensors: Heavy scalars+Temporal scalars*)

(*TadPoles*)
VarTadpole=Join[\[Lambda]1//Normal//Variables]//DeleteDuplicates;
SubTadpole=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarTadpole}];

TadPoleLight=Table[\[Lambda]1[[a]],{a,LightScalar[[;;,1]]}]//ReplaceAll[#,SubTadpole]&//SparseArray;
TadPoleTemp=Table[\[Lambda]1[[a]],{a,HeavyScalars[[;;,1]]}]//ReplaceAll[#,SubTadpole]&//SparseArray;


If[Length[TadPoleTemp//ArrayRules]==1,
TadPoleHeavy=SymmetrizedArray[{{1}->0},{nSH},Symmetric[{1}]]//SparseArray;
,
TadPoleHeavy=Delete[TadPoleTemp//ArrayRules//ReplaceAll[#,{x_}->{x+nv}]&,-1]//SparseArray;
];


\[Lambda]x1=Table[\[Lambda]3DSp[[a,b,c,d]],{a,LightScalar[[;;,1]]},{b,HeavyScalars[[;;,1]]},{c,HeavyScalars[[;;,1]]},{d,HeavyScalars[[;;,1]]}]//SparseArray;
\[Lambda]y1=Table[\[Lambda]3DSp[[a,b,c,d]],{a,LightScalar[[;;,1]]},{b,LightScalar[[;;,1]]},{c,LightScalar[[;;,1]]},{d,HeavyScalars[[;;,1]]}]//SparseArray;


If[Length[\[Lambda]x1//ArrayRules]==1,
\[Lambda]x=SymmetrizedArray[{{1,1,1,1}->0},{nSL,nSH,nSH,nSH},Symmetric[{2,3,4}]]//SparseArray;
,
\[Lambda]x=Delete[\[Lambda]x1//ArrayRules//ReplaceAll[#,{x_,y_,z_,w_}->{x,y+nv,z+nv,w+nv}]&,-1]//SparseArray;
];

If[Length[\[Lambda]y1//ArrayRules]==1,
\[Lambda]y=SymmetrizedArray[{{1,1,1,1}->0},{nSL,nSL,nSL,nSH},Symmetric[{1,2,3}]]//SparseArray;
,
\[Lambda]y=Delete[\[Lambda]y1//ArrayRules//ReplaceAll[#,{x_,y_,z_,w_}->{x,y,z,w+nv}]&,-1]//SparseArray;
];


(*Mixed Cubic tensors*)
\[Lambda]3CLight=Table[\[Lambda]3CSRedUS[[a,b,c]],{a,LightScalar[[;;,1]]},{b,LightScalar[[;;,1]]},{c,LightScalar[[;;,1]]}]//SparseArray;
\[Lambda]3CHeavy1=Table[\[Lambda]3CSRedUS[[a,b,c]],{a,HeavyScalars[[;;,1]]},{b,HeavyScalars[[;;,1]]},{c,HeavyScalars[[;;,1]]}]//SparseArray;
\[Lambda]3CTot=\[Lambda]3CSRedUS//SparseArray;

If[Length[\[Lambda]3CHeavy1//ArrayRules]==1,
\[Lambda]3CHeavy=SymmetrizedArray[{{1,1,1}->0},{nSH,nSH,nSH},Symmetric[{2,3}]]//SparseArray;
,
\[Lambda]3CHeavy=Delete[\[Lambda]3CHeavy1//ArrayRules//ReplaceAll[#,{x_,y_,z_}->{x+nv,y+nv,z+nv}]&,-1]//SparseArray;
];


\[Lambda]3Cx1=Table[\[Lambda]3CSRedUS[[a,b,c]],{a,LightScalar[[;;,1]]},{b,HeavyScalars[[;;,1]]},{c,HeavyScalars[[;;,1]]}]//SparseArray;
\[Lambda]3Cy1=Table[\[Lambda]3CSRedUS[[a,b,c]],{a,LightScalar[[;;,1]]},{b,LightScalar[[;;,1]]},{c,HeavyScalars[[;;,1]]}]//SparseArray;


If[Length[\[Lambda]3Cx1//ArrayRules]==1,
\[Lambda]3Cx=SymmetrizedArray[{{1,1,1}->0},{nSL,nSH,nSH},Symmetric[{2,3}]]//SparseArray;
,
\[Lambda]3Cx=Delete[\[Lambda]3Cx1//ArrayRules//ReplaceAll[#,{x_,y_,z_}->{x,y+nv,z+nv}]&,-1]//SparseArray;
];



If[Length[\[Lambda]3Cy1//ArrayRules]==1,
\[Lambda]3Cy=SymmetrizedArray[{{1,1,1}->0},{nSL,nSL,nSH},Symmetric[{1,2}]]//SparseArray;
,
\[Lambda]3Cy=Delete[\[Lambda]3Cy1//ArrayRules//ReplaceAll[#,{x_,y_,z_}->{x,y,z+nv}]&,-1]//SparseArray;
];

(*Mixed Hard-Soft mass matrix*)
\[Mu]ijMix1=Table[\[Mu]ij[[a,b]],{a,LightScalar[[;;,1]]},{b,HeavyScalars[[;;,1]]}]//SparseArray;

If[Length[\[Mu]ijMix1//ArrayRules]==1,
\[Mu]ijMix=SymmetrizedArray[{{1,1}->0},{nSL,nSH},Symmetric[{1}]]//SparseArray;
,
\[Mu]ijMix=Delete[\[Mu]ijMix1//ArrayRules//ReplaceAll[#,{x_,y_}->{x,y+nv}]&,-1]//SparseArray;
];


(*Scalar-Temporal-Vector Tensor*)
\[Lambda]KVLight=Table[\[Lambda]KVec[[a,b,c,d]],{a,1,nv},{b,1,nv},{c,LightScalar[[;;,1]]},{d,LightScalar[[;;,1]]}]//SparseArray;
\[Lambda]KVHeavy=Table[\[Lambda]3DSp[[a,b,c,d]],{a,HeavyScalars[[;;,1]]},{b,HeavyScalars[[;;,1]]},{c,LightScalar[[;;,1]]},{d,LightScalar[[;;,1]]}]//SparseArray;
A1=Delete[\[Lambda]KVLight//ArrayRules,-1];
A2=Delete[\[Lambda]KVHeavy//ArrayRules//ReplaceAll[#,{x_,y_,z_,w_}->{x+nv,y+nv,z,w}]&,-1];
\[Lambda]kVTot=SparseArray[Join[A1,A2]]//SparseArray;


VarGauge=GaugeCouplingNames;
SubGauge=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarGauge}];


(*Heavy-Scalar-Vector Tensor*)
gvssHeavy=Table[gvss[[a,b,c]],{a,1,nv},{b,HeavyScalars[[;;,1]]},{c,HeavyScalars[[;;,1]]}]//ReplaceAll[#,SubGauge]&//SparseArray;
gvssVec=gvvv//Normal//ReplaceAll[#,SubGauge]&//SparseArray;
gvssVTot=Table[ArrayFlatten[{{gvssVec[[a]],0},{0,gvssHeavy[[a]]}}],{a,1,Length[gvssVec]}]//SparseArray;

(*Light-Scalar-Vector Tensor*)
gvssL=Table[gvss[[a,c,d]],{a,1,nv},{c,LightScalar[[;;,1]]},{d,LightScalar[[;;,1]]}]//Normal//ReplaceAll[#,SubGauge]&//SparseArray;

(*Heavy-Mass matrix*)
MassHeavy=Table[\[Mu]ijSNLOP[[a,b]],{a,HeavyScalars[[;;,1]]},{b,HeavyScalars[[;;,1]]}]//SparseArray;
\[Mu]ijL=ArrayFlatten[{{Sqrt[\[Mu]abDef],0},{0,Sqrt[MassHeavy]}}]//SparseArray;

(*Light-Mass matrix*)
MassLight=Table[\[Mu]ijSNLOP[[a,b]],{a,LightScalar[[;;,1]]},{b,LightScalar[[;;,1]]}]//SparseArray;
\[Mu]ijLight=ArrayFlatten[MassLight]//SparseArray;

\[Lambda]K=\[Lambda]kVTot//SparseArray;
\[Lambda]4S=\[Lambda]4Light//SparseArray;
\[Lambda]4K=\[Lambda]KTotal//SparseArray;
gAvss=gvssVTot//SparseArray;
ToExpression[StringReplace[ToString[StandardForm[Join[\[Lambda]K,\[Lambda]4S,\[Lambda]4K,\[Lambda]x,\[Lambda]y,gAvss,gvssL,\[Mu]ijL]]],"DRalgo`Private`"->""]];

,
\[Lambda]3CTot=\[Lambda]3CSRedUS//SparseArray;
A3=Delete[\[Lambda]AAS//ArrayRules,-1];
\[Lambda]KTotal=SymmetrizedArray[A3,{nSH,nSH,nSH,nSH},Symmetric[{1,2,3,4}]]//SparseArray;
\[Lambda]3Cx=SymmetrizedArray[{{1,1,1}->0},{nSL,nSH,nSH},Symmetric[{2,3}]]//SparseArray;
\[Lambda]3Cx=Table[\[Lambda]vvsLS[[b,c,a]],{a,LightScalar[[;;,1]]},{b,1,nv},{c,nv}]//SparseArray;
\[Lambda]3CHeavy=SymmetrizedArray[{{1,1,1}->0},{nSH,nSH,nSH},Symmetric[{1,2,3}]]//SparseArray;
\[Lambda]3Cy=SymmetrizedArray[{{1,1,1}->0},{nSL,nSL,nSH},Symmetric[{1,2}]]//SparseArray;
\[Lambda]x=SymmetrizedArray[{{1,1,1,1}->0},{nSL,nSH,nSH,nSH},Symmetric[{2,3,4}]]//SparseArray;
\[Lambda]y=SymmetrizedArray[{{1,1,1,1}->0},{nSL,nSL,nSL,nSH},Symmetric[{1,2,3}]]//SparseArray;
\[Mu]ijMix=SymmetrizedArray[{{1,1}->0},{nSL,nSH},Symmetric[{1}]]//SparseArray;
\[Lambda]3CLight=Table[\[Lambda]3CSRedUS[[a,b,c]],{a,LightScalar[[;;,1]]},{b,LightScalar[[;;,1]]},{c,LightScalar[[;;,1]]}]//SparseArray;
\[Lambda]4Light=Table[\[Lambda]3DSp[[a,b,c,d]],{a,LightScalar[[;;,1]]},{b,LightScalar[[;;,1]]},{c,LightScalar[[;;,1]]},{d,LightScalar[[;;,1]]}]//SparseArray;

(*TadPoles*)
VarTadpole=Join[\[Lambda]1//Normal//Variables]//DeleteDuplicates;
SubTadpole=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarTadpole}];

TadPoleLight=Table[\[Lambda]1[[a]],{a,LightScalar[[;;,1]]}]//ReplaceAll[#,SubTadpole]&//SparseArray;
TadPoleHeavy=SymmetrizedArray[{{1}->0},{nSH},Symmetric[{1}]]//SparseArray;

\[Lambda]KVLight=Table[\[Lambda]KVec[[a,b,c,d]],{a,1,nv},{b,1,nv},{c,LightScalar[[;;,1]]},{d,LightScalar[[;;,1]]}]//SparseArray;
\[Lambda]kVTot=\[Lambda]KVLight//SparseArray;

VarGauge=GaugeCouplingNames;
SubGauge=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarGauge}];

gvssVec=gvvv//Normal//ReplaceAll[#,SubGauge]&//SparseArray;
gvssVTot=gvssVec//SparseArray; 
gvssL=Table[gvss[[a,c,d]],{a,1,nv},{c,LightScalar[[;;,1]]},{d,LightScalar[[;;,1]]}]//Normal//ReplaceAll[#,SubGauge]&//SparseArray;

\[Mu]ijL=Sqrt[\[Mu]abDef]//Normal//SparseArray;
\[Mu]ijLight=\[Mu]ijSNLOP//Normal//SparseArray;

\[Lambda]K=\[Lambda]kVTot//SparseArray;
\[Lambda]4S=\[Lambda]4Light//SparseArray;
\[Lambda]4K=\[Lambda]KTotal//SparseArray;
gAvss=gvssVTot//SparseArray;

];


];


(*
	Prints the beta functions in the ultrasoft theory.
*)
BetaFunctions3DUS[]:=Module[{},
If[verbose,Print["Finding SuperSoft \[Beta]-functions"]];


VarGauge=Join[\[Mu]ijLight//Normal//Variables]//DeleteDuplicates;
SubGauge=Table[c->Symbol[ToString[c]<>ToString["US"]],{c,VarGauge}];

(*This is the epsilon pole of the 3d sunset diagram. This is the only divergence in 3d.*)
I111=1/(4)1/(16 \[Pi]^2);

Coeff=1/2*2;
Contri1=Coeff*Flatten[TensorContract[HabijL,{1,2}]] . Flatten[\[Lambda]4Light,{1,2}];

Coeff=-1/4*(-1);
Contri2=Coeff*Flatten[TensorContract[HabijL,{3,4}]] . Flatten[HabijVL,{1,2}];

Coeff=1/3!;
Contri3=Coeff *Flatten[\[Lambda]4Light,{{1},{2,3,4}}] . Flatten[\[Lambda]4Light,{{1,2,3}}];

Coeff=1/2*((1+1/2));
Contri4=Coeff *Flatten[HabijVL,{{3},{1,2,4}}] . Flatten[HabijVL,{{1,2,4},{3}}];

Coeff=(-1)(-1/4) /2;
GabcdV2=TensorContract[GabcdV,{2,4}];
Contri5=Coeff*Flatten[GabcdV2] . Flatten[HabijVL,{1,2}];
 
Coeff=(-1)(-1)1/4*(20/4);
Contri6=Coeff*Flatten[GabcdV2] . Flatten[HabijVL,{1,2}];


SubGauge=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,GaugeCouplingNames}];

\[Delta]\[Mu]3dS=Contri1+   Contri2+ Contri3+  Contri4+ Contri5+  Contri6;(*Sum of all the diagrams*)
\[Beta]\[Mu]3ijS=4*I111* \[Delta]\[Mu]3dS//Normal //ReplaceAll[#,SubGauge]&//Simplify; (*The scalar-mass counter-term*)

(*Printing Result*)
(* Scalar Mass*)
VarGauge=Join[\[Mu]ijLight//Normal//Variables]//DeleteDuplicates;
SubGauge=Table[c->Symbol[ToString[c]<>ToString["temp"]],{c,VarGauge}];
SubGauge2=Table[Symbol[ToString[c]<>ToString["temp"]]->Symbol[ToString[c]<>ToString["US"]],{c,VarGauge}];


GaugeHelp=Table[Symbol[ToString[c]<>ToString["3d"]],{c,GaugeCouplingNames}];
VarUS=Join[{\[Lambda]3CLight//Normal//Variables,\[Lambda]4Light//Normal//Variables,GaugeHelp//Normal//Variables}]//DeleteDuplicates//Flatten[#,1]&;
SubUS=Table[c->Symbol[ToString[c]<>ToString["US"]],{c,VarUS}];


\[Lambda]4p=\[Mu]ijLight//Normal//ReplaceAll[#,SubGauge]&//SparseArray;
SolVar=\[Beta]\[Mu]3ijS-\[Lambda]4p//Normal;
QuarticVar=\[Lambda]4p//Normal//Variables;
ResMass=Solve[SolVar==0,QuarticVar]/.SubGauge2/.SubUS//Flatten[#,1]&; (*Finds the beta-function for each scalar mass*)


(*Printing Result*)
PrintPre=ResMass//Normal//FullSimplify//DeleteDuplicates;




(*Tadpoles*)
\[CapitalLambda]gHelp=TensorContract[HabijL,{{1,2}}]//SimplifySparse; 
Contri10=1/3!*I111*Activate@TensorContract[Inactive@TensorProduct[\[Lambda]4Light//SparseArray,\[Lambda]3CLight//SparseArray],{{2,5},{3,6},{4,7}}];
Contri11=1/2*I111*2Activate@TensorContract[Inactive@TensorProduct[\[Lambda]3CLight//SparseArray,\[CapitalLambda]gHelp//SparseArray],{{2,4},{3,5}}];

VarGauge=GaugeCouplingNames;
SubGauge=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarGauge}];

\[Delta]\[Mu]Tadpole3dUS=Contri10+Contri11//Normal//ReplaceAll[#,SubGauge]&//SparseArray;(*Same deal for tadpoles*)
\[Beta]\[Mu]TadpoleUS=4*\[Delta]\[Mu]Tadpole3dUS //Simplify;

VarGauge=Join[TadPoleLight//Normal//Variables]//DeleteDuplicates;
SubGauge=Table[c->Symbol[ToString[c]<>ToString["temp"]],{c,VarGauge}];
SubGauge2=Table[Symbol[ToString[c]<>ToString["temp"]]->Symbol[ToString[c]<>ToString["US"]],{c,VarGauge}];

\[Lambda]4p=TadPoleLight//Normal//ReplaceAll[#,SubGauge]&//SparseArray;
SolVar=\[Beta]\[Mu]TadpoleUS-\[Lambda]4p//Normal;
QuarticVar=\[Lambda]4p//Normal//Variables;
ResTadpole=Solve[SolVar==0,QuarticVar]/.SubGauge2/.SubUS//Flatten[#,1]&; (*Finds the beta-function for each scalar mass*)

PrintPre=Join[PrintPre,ResTadpole]//Normal//FullSimplify//DeleteDuplicates;

ToExpression[StringReplace[ToString[StandardForm[PrintPre]],"DRalgo`Private`"->""]]


];


{Contri\[Beta]UltraSoftToSoft,ContriTadPoleUltraSoftToSoft};


(*
	Calculates the 1-loop tadpole in the ultrasoft theory.
*)
TadPoleSS[]:=Module[{},
If[verbose,Print["Calculating 1-Loop Tadpoles"]];


fac=-1/(4 \[Pi]);
ContriL1=-fac/2*Simplify[Table[Sum[\[Mu]ijL[[ii,ii]] \[Lambda]3Cx[[i,ii,ii]],{ii,1,nSH}],{i,1,nSL}]];
ContriLSE=-ZSij . TadPoleLight;

TadPoleLightSS=TadPoleLight+ContriLSE-ContriL1;

ContriH1=-fac/2*Simplify[Table[Sum[\[Mu]ijL[[ii,ii]] \[Lambda]3CHeavy[[i,ii,ii]],{ii,1,nSH}],{i,1,nSH}]];

TadPoleHeavySS=TadPoleHeavy-ContriH1;
];


(*
	This is the mass that pops up if the heavy-scalar line carry a soft momenta
*)
HeavyScalarMassSS[]:=Module[{},
fac=1/(4 \[Pi]);
Contri1=1/2  *fac*Simplify[Table[Sum[\[Lambda]3CHeavy[[i,ii,jj]]\[Lambda]3CHeavy[[j,ii,jj]]/(\[Mu]ijL[[ii,ii]]+\[Mu]ijL[[jj,jj]]),{ii,1,nSH},{jj,1,nSH}],{i,1,nSH},{j,1,nSH}]];
fac=-1/(4 \[Pi]);
Contri2=-1/2*fac*Simplify[Table[Sum[\[Lambda]KTotal[[i,j,ii,ii]] \[Mu]ijL[[ii,ii]],{ii,1,nSH}],{i,1,nSH},{j,1,nSH}]]//FullSimplify;
fac=1/(4 \[Pi]);
Contri3=fac*Simplify[Table[Sum[\[Lambda]3Cx[[jj,ii,i]]\[Lambda]3Cx[[jj,ii,j]]/(\[Mu]ijL[[ii,ii]]),{ii,1,nSH},{jj,1,nSL}],{i,1,nSH},{j,1,nSH}]];
Contri4=-Simplify[Table[Sum[TadPoleHeavySS[[ii]]\[Lambda]3CHeavy[[i,j,ii]]/(\[Mu]ijL[[ii,ii]]^2),{ii,1,nSH}],{i,1,nSH},{j,1,nSH}]];

HeavyScalarMass=Contri1+Contri2+Contri3+Contri4;

HelpList=DeleteDuplicates@Flatten@Simplify[ HeavyScalarMass]//Sort;
HelpVar=Table[ \[Mu]H[a],{a,1,Delete[HelpList,1]//Length}];
HelpVarMod=RelationsBVariables[HelpList,HelpVar];
HelpSolveEffectiveHardM=Table[{Delete[HelpList,1][[a]]->HelpVarMod[[a]]},{a,1,Delete[HelpList,1]//Length}]//Flatten;
\[Mu]HEff=HeavyScalarMass//Normal//ReplaceAll[#,HelpSolveEffectiveHardM]&//SparseArray;

];


DiagonalTensor2[s_SparseArray,a_Integer,b_Integer] := With[
    {
    s1=Flatten[s,{{a},{b}}]
    },
     Table[i,{i,s1}]//Table[#[[i,i]],{i,1,Length[#]}]&//SparseArray
    ]


DiagonalTensor[s_SparseArray,a_Integer,b_Integer]:= Module[{},
If[a==1&&b==2,
    PermTens=Table[i,{i,s}]//Table[#[[i,i]],{i,1,Length[#]}]&//SparseArray
,
If[a==1,
   PermTens=Transpose[s,b<->2];
    PermTens=Table[i,{i,PermTens}]//Table[#[[i,i]],{i,1,Length[#]}]&//SparseArray;
   Transpose[PermTens,(b-1)<->1]
,
	   PermTens=Transpose[s,a<->1]//Transpose[#,b<->2]&;
    PermTens=Table[i,{i,PermTens}]//Table[#[[i,i]],{i,1,Length[#]}]&//SparseArray;
   Transpose[PermTens,(a-1)<->1]//Transpose[#,(b-1)<->2]&
]
]
    ]
