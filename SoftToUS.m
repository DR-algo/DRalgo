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
SubGauge=Table[c->Symbol[ToString[c]<>ToString["SS"]],{c,VarGauge}];


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
SubGauge=Table[c->Symbol[ToString[c]<>ToString["SS"]],{c,QuarticVar}];
NonZeroPos=SparseArray[\[Lambda]4S]["NonzeroPositions"];
SolVar=Extract[\[Lambda]4S-\[Lambda]3DSS,{#}]&/@NonZeroPos//DeleteDuplicates;

ResScalp=Reduce[SolVar==0,QuarticVar]//ToRules[#]&;
SolveTemp=QuarticVar/.ResScalp;
ResScal=Table[{ReplaceAll[QuarticVar[[i]],SubGauge]->SolveTemp[[i]]},{i,1,Length@QuarticVar}]//Flatten[#,1]&//ReplaceAll[#,IdentMatSS]&;



(* Scalar Cubics*)
VarGauge=Join[\[Lambda]3CLight//Normal//Variables]//DeleteDuplicates;
SubGauge=Table[c->Symbol[ToString[c]<>ToString["SS"]],{c,VarGauge}];
\[Lambda]3p=\[Lambda]3CLight//Normal//ReplaceAll[#,SubGauge]&;
SolVar=\[Lambda]3CSRedSS-\[Lambda]3p//Normal;
CubicVar=\[Lambda]3p//Normal//Variables;
ResCubicSS=Solve[SolVar==0,CubicVar]/.IdentMatSS//Flatten[#,1]&;




(*Printing Result*)
PrintPre=Join[ResScal,ResGauge,ResGaugeNASS,ResCubicSS]//Normal//FullSimplify//DeleteDuplicates;

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

Temp=-Simplify[Table[Sum[1/(\[Mu]ijL[[n,n]]^2)\[Lambda]3Cy[[i,j,n]]\[Lambda]3CHeavy[[k,l,n]],{n,1,nSH},{m,1,nSH}],{k,1,nSH},{l,1,nSH},{i,1,nSL},{j,1,nSL}]]//SparseArray;
\[Lambda]KTemp=\[Lambda]K+Temp;


Temp=-Simplify[Table[Sum[1/(\[Mu]ijL[[n,n]]^2)\[Lambda]3Cx[[i,jj,n]]\[Lambda]3Cy[[k,l,n]],{n,1,nSH}],{i,1,nSL},{k,1,nSL},{l,1,nSL},{jj,1,nSH}]];
\[Lambda]yTemp=\[Lambda]y+Temp;

Contri1Pre=1/(4 \[Pi]) /2Simplify[Table[Sum[\[Lambda]KTemp[[ii,jj,i,j]]\[Lambda]3Cx[[k,ii,jj]]1/(\[Mu]ijL[[ii,ii]]+\[Mu]ijL[[jj,jj]]),{ii,1,nSH},{jj,1,nSH}],{i,1,nSL},{j,1,nSL},{k,1,nSL}]];
Contri1=Contri1Pre+Transpose[Contri1Pre,{3,2,1}]+Transpose[Contri1Pre,{1,3,2}];

Contri2Pre=1/(4 \[Pi]) Simplify[Table[Sum[\[Lambda]yTemp[[i,j,ii,jj]]\[Lambda]3Cy[[k,ii,jj]]1/(\[Mu]ijL[[jj,jj]]),{ii,1,nSL},{jj,1,nSH}],{i,1,nSL},{j,1,nSL},{k,1,nSL}]];
Contri2=Contri2Pre+Transpose[Contri2Pre,{3,2,1}]+Transpose[Contri2Pre,{1,3,2}];

MassTemp=Table[MassHelpTriangle[\[Mu]ijLS[[n,n]],\[Mu]ijLS[[l,l]],\[Mu]ijLS[[m,m]]],{n,1,ns},{m,1,ns},{l,1,ns}]//SparseArray;
Contri3Pre=-Table[Sum[\[Lambda]3CTot[[i,n,m]]\[Lambda]3CTot[[j,l,n]]\[Lambda]3CTot[[k,l,m]]MassTemp[[l,m,n]],{l,1,ns},{m,1,ns},{n,1,ns}],{i,LightScalar[[;;,1]]},{j,LightScalar[[;;,1]]},{k,LightScalar[[;;,1]]}];
Contri3=Contri3Pre+Transpose[Contri3Pre,{3,2,1}]+Transpose[Contri3Pre,{1,3,2}];

ContriMixed=-1/(4 \[Pi])*1/2* Simplify[Table[Sum[(\[Mu]ijL[[m,m]])/\[Mu]ijL[[l,l]]^2(\[Lambda]3Cy[[i,j,l]]*\[Lambda]x[[k,m,m,l]]+\[Lambda]3Cy[[i,k,l]]*\[Lambda]x[[j,m,m,l]]+\[Lambda]3Cy[[k,j,l]]*\[Lambda]x[[i,m,m,l]]),{l,1,nSH},{m,1,nSH}],{i,1,nSL},{j,1,nSL},{k,1,nSL}]];
ContriSE=Simplify[Table[Sum[ZSij[[i,l]]\[Lambda]3CLight[[l,j,k]]+ZSij[[j,l]]\[Lambda]3CLight[[l,i,k]]+ZSij[[k,l]]\[Lambda]3CLight[[l,j,i]],{l,1,nSL}],{i,1,nSL},{j,1,nSL},{k,1,nSL}]];


\[Lambda]3CSSS=-Contri3-Contri1-Contri2-ContriSE-ContriMixed;

];


(*
	Calculates the 2-loop scalar mass in the ultrasoft theory.
*)
ScalarMass2LoopSS[]:=Module[{},
If[verbose,Print["Calculating 2-Loop Scalar Mass"]];

Contri1=-1/(16 \[Pi]^2)1/2Simplify[Table[Sum[ \[Lambda]K[[a,b,i,n]]\[Lambda]K[[a,b,n,j]](1/2+Log[\[Mu]3/( \[Mu]ijL[[a,a]]+ \[Mu]ijL[[b,b]])]),{a,1,nSH},{b,1,nSH},{n,1,nSL}],{i,1,nSL},{j,1,nSL}]];
Contri2=1/(16 \[Pi]^2)*(-1/2)Simplify[Table[Sum[\[Lambda]K[[n,m,i,j]]gAvss[[a,n,l]]gAvss[[a,l,m]](1/2+2Log[\[Mu]3/(2 \[Mu]ijL[[n,n]])]),{a,1,nv},{l,1,nSH},{n,1,nSH},{m,1,nSH}],{i,1,nSL},{j,1,nSL}]];
Contri3=1/(16 \[Pi]^2)*(1/4)Simplify[Table[Sum[HabijVL[[a,b,i,j]]gAvss[[a,n,l]]gAvss[[b,l,n]](-Log[\[Mu]3/(2 \[Mu]ijL[[n,n]])]),{a,1,nv},{b,1,nv},{l,1,nSH},{n,1,nSH}],{i,1,nSL},{j,1,nSL}]];
Contri4=1/4*(1/(16 \[Pi]^2))Simplify[Table[Sum[\[Lambda]K[[a,b,i,j]]\[Lambda]4K[[a,b,c,c]]\[Mu]ijL[[c,c]]/(\[Mu]ijL[[a,a]]+\[Mu]ijL[[b,b]]),{a,1,nSH},{b,1,nSH},{c,1,nSH}],{i,1,nSL},{j,1,nSL}]];
Contri5=-1/(16 \[Pi]^2)/3!Simplify[Table[Sum[ \[Lambda]x[[i,a,b,n]]\[Lambda]x[[j,a,b,n]](1/2+Log[\[Mu]3/( \[Mu]ijL[[a,a]]+ \[Mu]ijL[[b,b]]+ \[Mu]ijL[[n,n]])]),{a,1,nSH},{b,1,nSH},{n,1,nSH}],{i,1,nSL},{j,1,nSL}]];
Contri6=-1/(16 \[Pi]^2)/2Simplify[Table[Sum[ \[Lambda]y[[i,a,b,n]]\[Lambda]y[[j,a,b,n]](1/2+Log[\[Mu]3/( \[Mu]ijL[[n,n]])]),{a,1,nSL},{b,1,nSL},{n,1,nSH}],{i,1,nSL},{j,1,nSL}]];
Contri7=(1/(16 \[Pi]^2))/2Simplify[Table[Sum[\[Lambda]y[[i,j,a,b]]\[Lambda]x[[a,b,c,c]]\[Mu]ijL[[c,c]]/(\[Mu]ijL[[b,b]]),{a,1,nSL},{b,1,nSH},{c,1,nSH}],{i,1,nSL},{j,1,nSL}]];
Coupling=-1/4*1/(4 \[Pi])^2;
ContriMix1=Coupling*Simplify[Table[Sum[(\[Mu]ijL[[n,n]])\[Mu]ijL[[m,m]]^-2 \[Mu]ijL[[l,l]]\[Lambda]x[[i,n,n,m]]\[Lambda]x[[j,l,l,m]],{l,1,nSH},{m,1,nSH},{n,1,nSH}],{i,1,nSL},{j,1,nSL}]];
ContriMix2=-Simplify[Table[Sum[(\[Mu]ijMix[[i,m]])\[Mu]ijL[[m,m]]^-2 \[Mu]ijMix[[j,m]],{m,1,nSH}],{i,1,nSL},{j,1,nSL}]];

ContriC1=-1/2Simplify[Table[Sum[MassHelp1[\[Mu]ijLS[[jj,jj]],\[Mu]ijLS[[ii,ii]],\[Mu]ijLS[[kk,kk]],\[Mu]ijLS[[ll,ll]],\[Mu]ijLS[[mm,mm]]] \[Lambda]3CTot[[i,jj,ii]]\[Lambda]3CTot[[j,kk,ll]]\[Lambda]3CTot[[jj,mm,kk]]\[Lambda]3CTot[[ii,mm,ll]],{jj,1,ns},{ii,1,ns},{mm,1,ns},{kk,1,ns},{ll,1,ns}],{i,LightScalar[[;;,1]]},{j,LightScalar[[;;,1]]}]];

ContriC2=1/2*Simplify[Table[Sum[MassHelp2[\[Mu]ijLS[[ii,ii]],\[Mu]ijLS[[jj,jj]],\[Mu]ijLS[[mm,mm]],\[Mu]ijLS[[nn,nn]]] \[Lambda]4Tot[[i,ii,nn,mm]]\[Lambda]3CTot[[mm,nn,jj]]\[Lambda]3CTot[[j,ii,jj]],{jj,1,ns},{ii,1,ns},{mm,1,ns},{nn,1,ns}],{i,LightScalar[[;;,1]]},{j,LightScalar[[;;,1]]}]];
ContriC3=1/2*Simplify[Table[Sum[MassHelp2[\[Mu]ijLS[[ii,ii]],\[Mu]ijLS[[jj,jj]],\[Mu]ijLS[[mm,mm]],\[Mu]ijLS[[nn,nn]]] \[Lambda]4Tot[[j,ii,nn,mm]]\[Lambda]3CTot[[mm,nn,jj]]\[Lambda]3CTot[[i,ii,jj]],{jj,1,ns},{ii,1,ns},{mm,1,ns},{nn,1,ns}],{i,LightScalar[[;;,1]]},{j,LightScalar[[;;,1]]}]];
ContriC4=-1/4*Simplify[Table[Sum[MassHelp2[\[Mu]ijLS[[ii,ii]],\[Mu]ijLS[[jj,jj]],\[Mu]ijLS[[mm,mm]],\[Mu]ijLS[[nn,nn]]] \[Lambda]4Tot[[i,j,ii,jj]]\[Lambda]3CTot[[mm,jj,nn]]\[Lambda]3CTot[[mm,ii,nn]],{jj,1,ns},{ii,1,ns},{mm,1,ns},{nn,1,ns}],{i,LightScalar[[;;,1]]},{j,LightScalar[[;;,1]]}]];


\[Mu]ijTemp=\[Mu]ijLight+\[Mu]ijSSLO;
ContriF=-Table[Sum[ZSij[[i,n]]\[Mu]ijTemp[[n,j]]+ZSij[[j,n]]\[Mu]ijTemp[[n,i]],{n,1,nSL}],{i,1,nSL},{j,1,nSL}];

          


\[Mu]ijSSNLO2=(ContriC1+ContriC2+ContriC3+ContriC4)//Simplify//SparseArray;
\[Mu]ijSSNLO=(Contri1+Contri2+Contri3+ Contri4+Contri5+Contri6+Contri7+ContriMix1+ContriMix2+ContriF)//Simplify//SparseArray;

];



(*
	These functions appear in the matching.
*)
MassHelp1[0,0,0,0,0]:=0
MassHelp1[x_,0,0,0,0]:=1/(16 \[Pi]^2 x^2)
MassHelp1[0,x_,0,0,0]:=1/(16 \[Pi]^2 x^2)
MassHelp1[0,0,x_,0,0]:=1/(16 \[Pi]^2 x^2)
MassHelp1[0,0,0,x_,0]:=1/(16 \[Pi]^2 x^2)
MassHelp1[y_,y_,w_,z_,t_]:=1/(16 \[Pi]^2 (t+w+y) (t+y+z))
MassHelp1[x_,y_,z_,z_,t_]:=1/(16 \[Pi]^2 (t+x+z) (t+y+z))
MassHelp1[0,y_,0,0,t_]:=1/(16 \[Pi]^2 (t) (t+y))
MassHelp1[x_,y_,w_,z_,t_]:=(Log[\[Mu]/(t+w+x)]-Log[\[Mu]/(t+w+y)]-Log[\[Mu]/(t+x+z)]+Log[\[Mu]/(t+y+z)])/(16 \[Pi]^2 (w-z) (x-y))
MassHelp1[0,y_,0,z_,0]:=-((2 Log[\[Mu]/y]-2 Log[\[Mu]/(y+z)]+2 Log[\[Mu]/z]+1)/(32 \[Pi]^2 y z))
MassHelp1[x_,0,w_,0,0]:=-((2 Log[\[Mu]/x]-2 Log[\[Mu]/(x+w)]+2 Log[\[Mu]/w]+1)/(32 \[Pi]^2 x w))
MassHelp1[0,y_,w_,0,0]:=-((2 Log[\[Mu]/y]-2 Log[\[Mu]/(y+w)]+2 Log[\[Mu]/w]+1)/(32 \[Pi]^2 y w))
MassHelp1[x_,0,0,z_,0]:=-((2 Log[\[Mu]/x]-2 Log[\[Mu]/(x+z)]+2 Log[\[Mu]/z]+1)/(32 \[Pi]^2 x z))


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
MassHelp2[0,0,z_,t_]:=1/(16 \[Pi]^2 (t+z))
MassHelp2[x_,0,0,0]:=(-(2 Log[\[Mu]/x])-1)/(32 \[Pi]^2 x)
MassHelp2[0,y_,0,0]:=(-(2 Log[\[Mu]/y])-1)/(32 \[Pi]^2 y)
MassHelp2[y_,y_,z_,t_]:=1/(16 \[Pi]^2 (t+y+z))
MassHelp2[x_,y_,z_,t_]:=(Log[\[Mu]/(t+y+z)]-Log[\[Mu]/(t+x+z)])/(16 \[Pi]^2 (x-y))


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



HelpList=DeleteDuplicates@SparseArray[Simplify@Flatten[ (\[Lambda]4S+\[Lambda]3DSS)]]//Sort;
HelpVar=Table[\[Lambda]SS[a],{a,1,Delete[HelpList,1]//Length}];
HelpVarMod=RelationsBVariables[HelpList,HelpVar];
 HelpSolveQuarticS=Table[{Delete[HelpList,1][[a]]->HelpVarMod[[a]]},{a,1,Delete[HelpList,1]//Length}]//Flatten;
\[Lambda]3DSS=(\[Lambda]4S+\[Lambda]3DSS)//Normal//Simplify//ReplaceAll[#,HelpSolveQuarticS]&//SparseArray;



(*Cubic Tensor*)
HelpList=DeleteDuplicates@SparseArray[Flatten@Simplify[(\[Lambda]3CSSS+\[Lambda]3CLight)]]//Sort//FullSimplify;
HelpVar=Table[cSSSS[a],{a,1,Delete[HelpList,1]//Length}];
HelpVarMod=RelationsBVariables[HelpList,HelpVar];
HelpSolveCubicSSS=Table[{Delete[HelpList,1][[a]]->HelpVarMod[[a]]},{a,1,Delete[HelpList,1]//Length}]//Flatten//Simplify;
\[Lambda]3CSRedSS=(\[Lambda]3CSSS+\[Lambda]3CLight)//Normal//Simplify//FullSimplify//ReplaceAll[#,HelpSolveCubicSSS]&//SparseArray;

(*Scalar-Vector Tensor*)

HelpList=DeleteDuplicates@Simplify@Flatten[(HabijVL+GvvssTSS)]//Sort;
HelpVar=Table[ \[Lambda]VTSS[a],{a,1,Delete[HelpList,1]//Length}];
HelpVarMod=RelationsBVariables[HelpList,HelpVar];
HelpSolveVecTS=Table[{Delete[HelpList,1][[a]]->HelpVarMod[[a]]},{a,1,Delete[HelpList,1]//Length}]//Flatten;
\[Lambda]KVecTSS= (HabijVL+GvvssTSS)//Normal//Simplify//ReplaceAll[#,HelpSolveVecTS]&//SparseArray;


(*Scalar Mass*)
If[mode>=2,

HelpList=DeleteDuplicates@Flatten@Simplify[ xLO \[Mu]ijLight+xLO \[Mu]ijSSLO+xNLO(\[Mu]ijSSNLO+\[Mu]ijSSNLO2)]//Sort;
HelpVar=Table[ \[Mu]ijSS[a],{a,1,Delete[HelpList,1]//Length}];
HelpVarMod=RelationsBVariables[HelpList,HelpVar];
HelpSolveMassS=Table[{Delete[HelpList,1][[a]]->HelpVarMod[[a]]},{a,1,Delete[HelpList,1]//Length}]//Simplify//Flatten;
\[Mu]ijSNLOSS=xLO \[Mu]ijLight+xLO \[Mu]ijSSLO+xNLO(\[Mu]ijSSNLO+\[Mu]ijSSNLO2)//Normal//Simplify//ReplaceAll[#,HelpSolveMassS]&//SparseArray;
,
HelpList=DeleteDuplicates@Flatten@Simplify[ xLO \[Mu]ijLight+xLO \[Mu]ijSSLO]//Sort;
HelpVar=Table[ \[Mu]ijSS[a],{a,1,Delete[HelpList,1]//Length}];
HelpVarMod=RelationsBVariables[HelpList,HelpVar];
HelpSolveMassS=Table[{Delete[HelpList,1][[a]]->HelpVarMod[[a]]},{a,1,Delete[HelpList,1]//Length}]//Simplify//Flatten;
\[Mu]ijSNLOSS=xLO \[Mu]ijLight+xLO \[Mu]ijSSLO//Normal//Simplify//ReplaceAll[#,HelpSolveMassS]&//SparseArray;
];

(*Scalar Tadpoles*)
HelpList=DeleteDuplicates@Flatten@Simplify[ TadPoleLightSS]//Sort;
HelpVar=Table[ dSS[a],{a,1,Delete[HelpList,1]//Length}];
HelpVarMod=RelationsBVariables[HelpList,HelpVar];
HelpSolveTadpoleSS=Table[{Delete[HelpList,1][[a]]->HelpVarMod[[a]]},{a,1,Delete[HelpList,1]//Length}]//Flatten;
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
SubGauge=Table[c->Symbol[ToString[c]<>ToString["SS"]],{c,VarGauge}];


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
SubGauge=Table[c->Symbol[ToString[c]<>ToString["SS"]],{c,VarGauge}];

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

ContriTLTemp1=Simplify[Table[Sum[(1/(\[Mu]ijL[[n,n]]^2))\[Lambda]3Cy[[i,j,n]]\[Lambda]3Cy[[k,l,n]],{n,1,nSH}],{i,1,nSL},{j,1,nSL},{k,1,nSL},{l,1,nSL}]];

ContriTLTemp2=Simplify[Table[Sum[1/(\[Mu]ijL[[n,n]]^2)\[Mu]HEff[[n,m]]1/(\[Mu]ijL[[m,m]]^2)\[Lambda]3Cy[[i,j,n]]\[Lambda]3Cy[[k,l,m]],{n,1,nSH},{m,1,nSH}],{i,1,nSL},{j,1,nSL},{k,1,nSL},{l,1,nSL}]];
ContriTLTemp=ContriTLTemp1+ContriTLTemp1//Simplify;
ContriTL=Table[ContriTLTemp[[i,j,k,l]]+ContriTLTemp[[i,k,j,l]]+ContriTLTemp[[i,l,j,k]],{i,1,nSL},{j,1,nSL},{k,1,nSL},{l,1,nSL}];
(*TadPole*)
Temp=-Simplify[Table[Sum[1/(\[Mu]ijL[[n,n]]^2)1/(\[Mu]ijL[[m,m]]^2)TadPoleHeavySS[[m]](\[Lambda]K[[n,m,i,j]]\[Lambda]3Cy[[k,l,n]]+\[Lambda]K[[n,m,k,l]]\[Lambda]3Cy[[j,i,n]]),{n,1,nSH},{m,1,nSH}],{i,1,nSL},{j,1,nSL},{k,1,nSL},{l,1,nSL}]]//Simplify;
ContriTadpole=Table[Temp[[i,j,k,l]]+Temp[[i,k,j,l]]+Temp[[i,l,j,k]],{i,1,nSL},{j,1,nSL},{k,1,nSL},{l,1,nSL}];

(*Maybe include Cubics later*)
Temp=0*Simplify[Table[Sum[(1/(\[Mu]ijL[[n,n]]^2))\[Lambda]3Cy[[i,j,n]]\[Lambda]3CHeavy[[k,l,n]],{n,1,nSH}],{k,1,nSH},{l,1,nSH},{i,1,nSL},{j,1,nSL}]];
\[Lambda]KEff=\[Lambda]K-Temp//SparseArray;
Temp=0*Simplify[Table[Sum[(1/(\[Mu]ijL[[n,n]]^2))\[Lambda]3Cy[[i,j,n]]\[Lambda]3Cx[[k,l,n]],{n,1,nSH}],{i,1,nSL},{j,1,nSL},{k,1,nSL},{l,1,nSH}]];
\[Lambda]yEff=\[Lambda]y-Temp//SparseArray;
\[Lambda]yEff2=Table[\[Lambda]y[[i,j,k,l]]-Temp[[i,j,k,l]]-Temp[[i,k,j,l]]-Temp[[k,j,i,l]],{i,1,nSL},{j,1,nSL},{k,1,nSL},{l,1,nSH}]//SparseArray;

(*Loop level*)
CouplingSSSS=1/2*1/(4 \[Pi])TensorProduct[\[Lambda]KEff,\[Lambda]KEff];
Temp=Table[Sum[1/(\[Mu]ijL[[n,n]]+\[Mu]ijL[[m,m]])(CouplingSSSS[[n,m,i,j,n,m,k,l]]),{n,1,nSH},{m,1,nSH}],{i,1,nSL},{j,1,nSL},{k,1,nSL},{l,1,nSL}]//SparseArray;
ContriSS=Table[Temp[[i,j,k,l]]+Temp[[i,k,j,l]]+Temp[[i,l,j,k]],{i,1,nSL},{j,1,nSL},{k,1,nSL},{l,1,nSL}];
CouplingSSSS=1/(4 \[Pi])TensorProduct[\[Lambda]yEff,\[Lambda]yEff];
Temp=Simplify[Table[Sum[1/(\[Mu]ijL[[m,m]])(CouplingSSSS[[i,j,n,m,k,l,n,m]]),{n,1,nSL},{m,1,nSH}],{i,1,nSL},{j,1,nSL},{k,1,nSL},{l,1,nSL}]]//SparseArray;
ContriSS2=Table[Temp[[i,j,k,l]]+Temp[[i,k,j,l]]+Temp[[i,l,j,k]],{i,1,nSL},{j,1,nSL},{k,1,nSL},{l,1,nSL}];
CouplingSSSS=-1/(4 \[Pi])*1/2TensorProduct[\[Lambda]yEff2,\[Lambda]x]//SparseArray;
Temp=Simplify[Table[Sum[1/(\[Mu]ijL[[m,m]]^2)\[Mu]ijL[[n,n]](CouplingSSSS[[i,j,k,m,l,n,n,m]]),{n,1,nSH},{m,1,nSH}],{i,1,nSL},{j,1,nSL},{k,1,nSL},{l,1,nSL}]];
ContriSS3=Table[Temp[[i,j,k,l]]+Temp[[i,j,l,k]]+Temp[[i,l,k,k]]+Temp[[l,j,k,i]],{i,1,nSL},{j,1,nSL},{k,1,nSL},{l,1,nSL}]//SparseArray;
Temp=-Simplify[Table[Sum[1/(\[Mu]ijL[[m,m]]^2)1/(\[Mu]ijL[[n,n]]^2)\[Lambda]yEff2[[i,j,k,m]]\[Lambda]3Cx[[l,m,n]]TadPoleHeavySS[[n]],{n,1,nSH},{m,1,nSH}],{i,1,nSL},{j,1,nSL},{k,1,nSL},{l,1,nSL}]]//SparseArray;
CouplingSSTadpole=Table[Temp[[i,j,k,l]]+Temp[[i,j,l,k]]+Temp[[i,l,k,k]]+Temp[[l,j,k,i]],{i,1,nSL},{j,1,nSL},{k,1,nSL},{l,1,nSL}]//SparseArray;

MassHelp5[x_]:=-(1/(4 \[Pi] x^3));
MassHelp4[x_,y_,z_]:=1/(4 \[Pi] (x+y) (x+z) (y+z));
MassHelp[x_,y_,w_,z_]:=(w+x+y+z)/(4 \[Pi] (w+x) (w+y) (w+z) (x+y) (x+z) (y+z));
MassHelp2[w_,z_]:=-((w^2+w z+z^2)/(4 \[Pi] w^3 z^3 (w+z)));
MassHelp3[w_]:=1/(4 \[Pi] w^5);

ContriSSSSTemp=Table[Sum[\[Lambda]3Cx[[i,n,m]]\[Lambda]3Cx[[j,m,l]]\[Lambda]3Cx[[k,l,o]]\[Lambda]3Cx[[t,n,o]] MassHelp[\[Mu]ijL[[n,n]],\[Mu]ijL[[m,m]],\[Mu]ijL[[l,l]],\[Mu]ijL[[o,o]]],{n,nv+1,nSH},{m,nv+1,nSH},{l,nv+1,nSH},{o,nv+1,nSH}],{i,1,nSL},{j,1,nSL},{k,1,nSL},{t,1,nSL}]//Simplify;
ContriCubic1=2Table[ContriSSSSTemp[[i,j,k,t]]+ContriSSSSTemp[[i,k,j,t]]+ContriSSSSTemp[[i,t,k,j]],{i,1,nSL},{j,1,nSL},{k,1,nSL},{t,1,nSL}];
ContriSSSSTemp=Table[Sum[\[Lambda]3Cy[[i,m,n]]\[Lambda]3Cy[[j,m,l]]\[Lambda]3Cx[[k,l,o]]\[Lambda]3Cx[[t,n,o]] MassHelp[\[Mu]ijL[[n,n]],0,\[Mu]ijL[[l,l]],\[Mu]ijL[[o,o]]],{n,nv+1,nSH},{m,1,nSL},{l,nv+1,nSH},{o,nv+1,nSH}],{i,1,nSL},{j,1,nSL},{k,1,nSL},{t,1,nSL}]//Simplify;
ContriCubic2=Table[(ContriSSSSTemp[[i,j,k,t]]+ContriSSSSTemp[[j,i,k,t]])+(ContriSSSSTemp[[i,k,j,t]]+ContriSSSSTemp[[k,i,j,t]])+(ContriSSSSTemp[[i,t,k,j]]+ContriSSSSTemp[[t,i,k,j]])+(ContriSSSSTemp[[k,t,i,j]]+ContriSSSSTemp[[t,k,i,j]])+(ContriSSSSTemp[[k,j,i,t]]+ContriSSSSTemp[[j,k,i,t]])+(ContriSSSSTemp[[j,t,i,k]]+ContriSSSSTemp[[t,j,i,k]]),{i,1,nSL},{j,1,nSL},{k,1,nSL},{t,1,nSL}];
 
ContriSSSSTemp=Table[Sum[\[Lambda]3Cy[[i,m,n]]\[Lambda]3CLight[[j,m,l]]\[Lambda]3Cy[[k,l,o]]\[Lambda]3Cx[[t,n,o]] MassHelp2[\[Mu]ijL[[n,n]],\[Mu]ijL[[o,o]]],{n,nv+1,nSH},{m,1,nSL},{l,1,nSL},{o,nv+1,nSH}],{i,1,nSL},{j,1,nSL},{k,1,nSL},{t,1,nSL}]//Simplify//SparseArray;
ContriCubic3H=2Table[ContriSSSSTemp[[i,j,k,t]]+ContriSSSSTemp[[k,j,t,i]]+ContriSSSSTemp[[t,j,t,k]],{i,1,nSL},{j,1,nSL},{k,1,nSL},{t,1,nSL}];
ContriCubic3=Table[ContriCubic3H[[i,j,k,t]]+ContriCubic3H[[j,i,k,t]]+ContriCubic3H[[i,k,j,t]]+ContriCubic3H[[i,t,k,j]],{i,1,nSL},{j,1,nSL},{k,1,nSL},{t,1,nSL}];
ContriSSSSTemp=Table[Sum[\[Lambda]3Cy[[i,m,n]]\[Lambda]3Cy[[j,m,l]]\[Lambda]3Cy[[k,o,l]]\[Lambda]3Cy[[t,o,n]] MassHelp2[\[Mu]ijL[[n,n]],\[Mu]ijL[[l,l]]],{n,nv+1,nSH},{m,1,nSL},{l,nv+1,nSH},{o,1,nSL}],{i,1,nSL},{j,1,nSL},{k,1,nSL},{t,1,nSL}]//Simplify//SparseArray;
ContriCubic3H=Table[ContriSSSSTemp[[i,j,k,t]]+ContriSSSSTemp[[j,i,k,t]],{i,1,nSL},{j,1,nSL},{k,1,nSL},{t,1,nSL}];
ContriCubic4=Table[ContriCubic3H[[i,j,k,t]]+ContriCubic3H[[i,k,j,t]]+ContriCubic3H[[i,t,k,j]],{i,1,nSL},{j,1,nSL},{k,1,nSL},{t,1,nSL}];
ContriSSSSTemp=Table[Sum[\[Lambda]3Cy[[i,m,n]]\[Lambda]3CLight[[j,m,l]]\[Lambda]3CLight[[k,l,o]]\[Lambda]3Cy[[t,o,n]] MassHelp3[\[Mu]ijL[[n,n]]],{n,nv+1,nSH},{m,1,nSL},{l,1,nSL},{o,1,nSL}],{i,1,nSL},{j,1,nSL},{k,1,nSL},{t,1,nSL}]//Simplify//SparseArray;
ContriCubic5=2Table[ContriSSSSTemp[[i,j,k,t]]+ContriSSSSTemp[[i,k,j,t]]+ContriSSSSTemp[[i,t,k,j]],{i,1,nSL},{j,1,nSL},{k,1,nSL},{t,1,nSL}];


ContriSSSSTemp=-Table[Sum[\[Lambda]K[[n,m,i,j]]\[Lambda]3Cx[[k,n,l]]\[Lambda]3Cx[[t,m,l]] MassHelp4[\[Mu]ijL[[n,n]],\[Mu]ijL[[m,m]],\[Mu]ijL[[l,l]]],{n,nv+1,nSH},{m,nv+1,nSH},{l,nv+1,nSH}],{i,1,nSL},{j,1,nSL},{k,1,nSL},{t,1,nSL}]//Simplify//SparseArray;
ContriCubic6=Table[(ContriSSSSTemp[[i,j,k,t]])+(ContriSSSSTemp[[i,k,j,t]])+(ContriSSSSTemp[[i,t,k,j]])+(ContriSSSSTemp[[k,t,i,j]])+(ContriSSSSTemp[[k,j,i,t]])+(ContriSSSSTemp[[j,t,i,k]]),{i,1,nSL},{j,1,nSL},{k,1,nSL},{t,1,nSL}];
ContriSSSSTemp=-Table[Sum[\[Lambda]K[[n,m,i,j]]\[Lambda]3Cy[[k,l,n]]\[Lambda]3Cy[[t,l,m]] MassHelp4[\[Mu]ijL[[n,n]],\[Mu]ijL[[m,m]],0],{n,nv+1,nSH},{m,nv+1,nSH},{l,1,nSL}],{i,1,nSL},{j,1,nSL},{k,1,nSL},{t,1,nSL}]//Simplify//SparseArray;
ContriCubic7=Table[(ContriSSSSTemp[[i,j,k,t]])+(ContriSSSSTemp[[i,k,j,t]])+(ContriSSSSTemp[[i,t,k,j]])+(ContriSSSSTemp[[k,t,i,j]])+(ContriSSSSTemp[[k,j,i,t]])+(ContriSSSSTemp[[j,t,i,k]]),{i,1,nSL},{j,1,nSL},{k,1,nSL},{t,1,nSL}];
ContriSSSSTemp=-Table[Sum[\[Lambda]y[[i,j,n,m]]\[Lambda]3CLight[[k,l,n]]\[Lambda]3Cy[[t,l,m]] MassHelp5[\[Mu]ijL[[m,m]]],{n,1,nSL},{m,nv+1,nSH},{l,1,nSL}],{i,1,nSL},{j,1,nSL},{k,1,nSL},{t,1,nSL}]//Simplify//SparseArray;
ContriCubic8=Table[(ContriSSSSTemp[[i,j,k,t]])+(ContriSSSSTemp[[i,k,j,t]])+(ContriSSSSTemp[[i,t,k,j]])+(ContriSSSSTemp[[k,t,i,j]])+(ContriSSSSTemp[[k,j,i,t]])+(ContriSSSSTemp[[j,t,i,k]]),{i,1,nSL},{j,1,nSL},{k,1,nSL},{t,1,nSL}];
ContriSSSSTemp=-Table[Sum[\[Lambda]y[[i,j,n,m]]\[Lambda]3Cy[[k,n,l]]\[Lambda]3Cx[[t,l,m]] MassHelp4[\[Mu]ijL[[l,l]],\[Mu]ijL[[m,m]],0],{n,1,nSL},{m,nv+1,nSH},{l,nv+1,nSH}],{i,1,nSL},{j,1,nSL},{k,1,nSL},{t,1,nSL}]//Simplify//SparseArray;
ContriCubic9=Table[(ContriSSSSTemp[[i,j,k,t]])+(ContriSSSSTemp[[i,k,j,t]])+(ContriSSSSTemp[[i,t,k,j]])+(ContriSSSSTemp[[k,t,i,j]])+(ContriSSSSTemp[[k,j,i,t]])+(ContriSSSSTemp[[j,t,i,k]]),{i,1,nSL},{j,1,nSL},{k,1,nSL},{t,1,nSL}];
ContriSSSSTemp=-Table[Sum[\[Lambda]4S[[i,j,n,m]]\[Lambda]3Cy[[k,n,l]]\[Lambda]3Cx[[t,m,l]] MassHelp5[\[Mu]ijL[[l,l]]],{n,1,nSL},{m,1,nSL},{l,nv+1,nSH}],{i,1,nSL},{j,1,nSL},{k,1,nSL},{t,1,nSL}]//Simplify//SparseArray;
ContriCubic10=Table[(ContriSSSSTemp[[i,j,k,t]])+(ContriSSSSTemp[[i,k,j,t]])+(ContriSSSSTemp[[i,t,k,j]])+(ContriSSSSTemp[[k,t,i,j]])+(ContriSSSSTemp[[k,j,i,t]])+(ContriSSSSTemp[[j,t,i,k]]),{i,1,nSL},{j,1,nSL},{k,1,nSL},{t,1,nSL}];


CouplingSE=-Inactivate[TensorProduct[ZSij,\[Lambda]4S]];
ContriSETemp=Simplify[Activate @ TensorContract[CouplingSE, {{2,3}}]]//SparseArray;
ContriSE=ContriSETemp+Transpose[ContriSETemp,{2,1,3,4}]+Transpose[ContriSETemp,{3,1,2,4}]+Transpose[ContriSETemp,{4,1,2,3}]//Simplify;

\[Lambda]3DSS=ContriSE-ContriTL-CouplingSSTadpole- ContriSS-ContriSS2-ContriSS3-ContriCubic1-ContriCubic2-ContriCubic3-ContriCubic4-ContriCubic5-ContriCubic6-ContriCubic7-ContriCubic8-ContriCubic9-ContriCubic10;
];


(*
	Calculates non-abelian couplings in the ultrasoft theory.
*)
NonAbelianCouplingSS[]:=Module[{},
ContriAnomVV= Simplify[Table[Sum[ ZLij[[c,d]]gvvvSS[[a,b,d]],{d,1,nv}],{a,1,nv},{b,1,nv},{c,1,nv}]];
GgvvvSS=-ContriAnomVV;
];


(*
	Calculates vector couplings in the ultrasoft theory.
*)
ScalarVectorCouplingSS[]:=Module[{},
If[verbose,Print["Calculating Vector-Scalar Couplings"]];

(*Cubic contribution*)
VarGauge=GaugeCouplingNames;
SubGauge=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarGauge}];
gvssM=gvss//Normal//ReplaceAll[#,SubGauge]&//SparseArray;
HabijM=Transpose[Activate @ TensorContract[
        Inactive[TensorProduct][gvssM,gvssM], { {3, 5}}],{1,3,2,4}]//SparseArray;
HabijMV=HabijM+Transpose[HabijM,{2,1,3,4}]//SparseArray;
MassTemp=Table[MassHelpTriangle[\[Mu]ijLS[[n,n]],\[Mu]ijLS[[l,l]],\[Mu]ijLS[[m,m]]],{n,1,ns},{m,1,ns},{l,1,ns}]//SparseArray;

ContriC1=Table[Sum[\[Lambda]3CTot[[i,l,n]]\[Lambda]3CTot[[j,l,m]]MassTemp[[l,m,n]]HabijMV[[a,b,n,m]],{n,1,ns},{m,1,ns},{l,1,ns}],{a,1,nv},{b,1,nv},{i,LightScalar[[;;,1]]},{j,LightScalar[[;;,1]]}]//SparseArray;
(* Self-Energy contribution*)



ContriSEVector=-Simplify[Table[Sum[ZLij[[a,c]](HabijVL[[c,b,i,j]])+ZLij[[b,c]](HabijVL[[a,c,i,j]]),{c,1,nv}],{a,1,nv},{b,1,nv},{i,1,nSL},{j,1,nSL}]];
ContriSEScalar=-Simplify[Table[Sum[ZSij[[i,k]](HabijVL[[a,b,k,j]])+ZSij[[j,k]](HabijVL[[a,b,i,k]]),{k,1,nSL}],{a,1,nv},{b,1,nv},{i,1,nSL},{j,1,nSL}]];


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

\[Lambda]3CSRed=\[Lambda]3//Normal//ReplaceAll[#,SubCubic]&//SparseArray;


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
A1=Delete[\[Lambda]4Heavy//ArrayRules//ReplaceAll[#,{x_,y_,z_,w_}->{x+nv,y+nv,z+nv,w+nv}]&,-1];
A2=Delete[\[Lambda]KHeavy//ArrayRules//ReplaceAll[#,{x_,y_,z_,w_}->{x,y,z+nv,w+nv}]&,-1];
\[Lambda]KTotal=SymmetrizedArray[Join[A1,A2],{nSH,nSH,nSH,nSH},Symmetric[{1,2,3,4}]]//SparseArray; (*Total heavy-scalar tensors: Heavy scalars+Temporal scalars*)

(*TadPoles*)
VarTadpole=Join[\[Lambda]1//Normal//Variables]//DeleteDuplicates;
SubTadpole=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarTadpole}];

TadPoleLight=Table[\[Lambda]1[[a]],{a,LightScalar[[;;,1]]}]//ReplaceAll[#,SubMass]&//SparseArray;
TadPoleTemp=Table[\[Lambda]1[[a]],{a,HeavyScalars[[;;,1]]}]//ReplaceAll[#,SubMass]&//SparseArray;


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
\[Lambda]3CLight=Table[\[Lambda]3CSRed[[a,b,c]],{a,LightScalar[[;;,1]]},{b,LightScalar[[;;,1]]},{c,LightScalar[[;;,1]]}]//SparseArray;
\[Lambda]3CHeavy1=Table[\[Lambda]3CSRed[[a,b,c]],{a,HeavyScalars[[;;,1]]},{b,HeavyScalars[[;;,1]]},{c,HeavyScalars[[;;,1]]}]//SparseArray;
\[Lambda]3CTot=\[Lambda]3CSRed//SparseArray;

If[Length[\[Lambda]3CHeavy1//ArrayRules]==1,
\[Lambda]3CHeavy=SymmetrizedArray[{{1,1,1}->0},{nSH,nSH,nSH},Symmetric[{2,3}]]//SparseArray;
,
\[Lambda]3CHeavy=Delete[\[Lambda]3CHeavy1//ArrayRules//ReplaceAll[#,{x_,y_,z_}->{x+nv,y+nv,z+nv}]&,-1]//SparseArray;
];


\[Lambda]3Cx1=Table[\[Lambda]3CSRed[[a,b,c]],{a,LightScalar[[;;,1]]},{b,HeavyScalars[[;;,1]]},{c,HeavyScalars[[;;,1]]}]//SparseArray;
\[Lambda]3Cy1=Table[\[Lambda]3CSRed[[a,b,c]],{a,LightScalar[[;;,1]]},{b,LightScalar[[;;,1]]},{c,HeavyScalars[[;;,1]]}]//SparseArray;


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
\[Lambda]3CTot=\[Lambda]3CSRed//SparseArray;
\[Lambda]KTotal=SymmetrizedArray[{{1,1,1,1}->0},{nSH,nSH,nSH,nSH},Symmetric[{1,2,3,4}]]//SparseArray;
\[Lambda]3Cx=SymmetrizedArray[{{1,1,1}->0},{nSL,nSH,nSH},Symmetric[{2,3}]]//SparseArray;
\[Lambda]3Cx=Table[\[Lambda]vvsLS[[b,c,a]],{a,LightScalar[[;;,1]]},{b,1,nv},{c,nv}]//SparseArray;
\[Lambda]3CHeavy=SymmetrizedArray[{{1,1,1}->0},{nSH,nSH,nSH},Symmetric[{1,2,3}]]//SparseArray;
\[Lambda]3Cy=SymmetrizedArray[{{1,1,1}->0},{nSL,nSL,nSH},Symmetric[{1,2}]]//SparseArray;
\[Lambda]x=SymmetrizedArray[{{1,1,1,1}->0},{nSL,nSH,nSH,nSH},Symmetric[{2,3,4}]]//SparseArray;
\[Lambda]y=SymmetrizedArray[{{1,1,1,1}->0},{nSL,nSL,nSL,nSH},Symmetric[{1,2,3}]]//SparseArray;
\[Mu]ijMix=SymmetrizedArray[{{1,1}->0},{nSL,nSH},Symmetric[{1}]]//SparseArray;
\[Lambda]3CLight=Table[\[Lambda]3CSRed[[a,b,c]],{a,LightScalar[[;;,1]]},{b,LightScalar[[;;,1]]},{c,LightScalar[[;;,1]]}]//SparseArray;
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


ToExpression[StringReplace[ToString[StandardForm[Join[\[Lambda]K,\[Lambda]4S,\[Lambda]4K,gAvss,gvssL,\[Mu]ijL]]],"DRalgo`Private`"->""]];
];


];


(*
	Prints the beta functions in the ultrasoft theory.
*)
BetaFunctions3DUS[]:=Module[{},
If[verbose,Print["Finding SuperSoft \[Beta]-functions"]];


VarGauge=Join[\[Mu]ijLight//Normal//Variables]//DeleteDuplicates;
SubGauge=Table[c->Symbol[ToString[c]<>ToString["SS"]],{c,VarGauge}];

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

Coeff=(-1/4) /2;
GabcdV2=TensorContract[GabcdV,{2,4}];
Contri5=Coeff*Flatten[GabcdV2] . Flatten[HabijVL,{1,2}];
 
Coeff=(-1)1/4*(20/4);
Contri6=Coeff*Flatten[GabcdV2] . Flatten[HabijVL,{1,2}];


VarGauge=GaugeCouplingNames;
SubGauge=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarGauge}];

\[Delta]\[Mu]3dS=Contri1+  Contri2+ Contri3+ Contri4+Contri5+ Contri6;(*Sum of all the diagrams*)

\[Beta]\[Mu]3ijS=4*I111* \[Delta]\[Mu]3dS //ReplaceAll[#,SubGauge]&//Simplify; (*The scalar-mass counter-term*)

(* Scalar Mass*)
VarGauge=Join[\[Mu]ij//Normal//Variables]//DeleteDuplicates;
SubGauge=Table[c->Symbol[ToString[c]<>ToString["temp"]],{c,VarGauge}];
SubGauge2=Table[Symbol[ToString[c]<>ToString["temp"]]->c,{c,VarGauge}];


\[Lambda]4p=\[Mu]ijLight//Normal//ReplaceAll[#,SubGauge]&//SparseArray;
SolVar=\[Beta]\[Mu]3ijS-\[Lambda]4p//Normal;
QuarticVar=\[Lambda]4p//Normal//Variables;
ResMass=Solve[SolVar==0,QuarticVar]/.SubGauge2//Flatten[#,1]&; (*Finds the beta-function for each scalar mass*)


(*Printing Result*)
PrintPre=ResMass//Normal//FullSimplify//DeleteDuplicates;

ToExpression[StringReplace[ToString[StandardForm[PrintPre]],"DRalgo`Private`"->""]]


];


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
\[Mu]HEff=HeavyScalarMass//Normal//Simplify//ReplaceAll[#,HelpSolveEffectiveHardM]&//SparseArray;

];
