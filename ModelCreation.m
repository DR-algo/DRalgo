(* ::Package:: *)

(* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ *)

(* :Title: ModelCreation													*)

(*
	This software is covered by the GNU General Public License 3.
	Copyright (C) 2021-2022
					Andreas Ekstedt,
					Philipp Schicho,
					Tuomas V.I. Tenkanen
*)

(* :Summary:	Creating models from model file.							*)

(*This file uses functions from GroupMath. A citation to arXiv:2011.01764 should*)
(*be used if any of these functions are used.*)
(* ------------------------------------------------------------------------ *)


(*
	This module creates representation matrices for the gauge sector. Other coupling tensors are allocated as zero-tensors.
*)
AllocateTensors[GroupI_,RepNameI_,CouplingNameI_,FermionRepI_,ScalarRepI_]:=Module[{GroupP=GroupI,RepNameP=RepNameI,CouplingNameP=CouplingNameI,FermionRepP=FermionRepI,ScalarRepP=ScalarRepI},
Clear[a,b,c,d];
Needs["GroupMath`"];
GroupMathCleared=False;

GroupP=ToGM[GroupP]; (*Converts the Group from DR input to GroupMath input*)
gen=RepMatrices[GroupP,RepNameP]; (*The representation matrices*)

SizeGroups= Quiet[Table[Length[RepMatrices[GroupP[[i]],RepNameP[[i]]]],{i,1,Length[GroupP]}]]; (*The number of components for each group*)
SizeGroupsMod=SizeGroups//ReplaceAll[#,0->1]&;  (*Fix for U1 groups*)
TotalComponents=Total[SizeGroupsMod]; (**)
GaugeIndices=ConstantArray[0,{Length[GroupP]}]; (*Position of every gauge rep*)


GaugeIndices[[1]]=1;;SizeGroupsMod[[1]];
Do[
helpInd=GaugeIndices[[i-1]]/.(a_;;b_)-> b;
GaugeIndices[[i]]=(helpInd+1);;(helpInd+SizeGroupsMod[[i]]);
,{i,2,Length[GroupP]}
];
  

(*Create coupling matrix. This matrix contains the name of coupling constants, and incorporates correct normalization of generators*)
(*This normalization is required when finding structure constants through commutation relations*)
Components=ConstantArray[0,Length[GroupP]+1];
Components[[1]]=1;
Do[
Components[[i]]=Components[[i-1]]+SizeGroupsMod[[i-1]]
,
{i,2,Length[GroupP]+1}
];

Do[
Components[[i]]=Components[[i-1]]+SizeGroupsMod[[i-1]]
,
{i,2,Length[GroupP]+1}
];


CouplingArray=ConstantArray[0,{TotalComponents,TotalComponents}];
Do[
(*Multiplies the relevant copupling-constant for each group*)
CouplingArray[[Components[[i]];;Components[[i+1]]-1,Components[[i]];;Components[[i+1]]-1]]=CouplingNameP[[i]]*IdentityMatrix[SizeGroupsMod[[i]]];
,
{i,1,Length[GroupP]}
];
(*Structure constants*)
Temp=I GaugeRep[GroupP]//SparseArray;
fabc=Table[ Temp[[a]] CouplingArray[[a,a]],{a,1,TotalComponents}]//SparseArray;

(*Fix for non-numerical U1s*)
PosU1=Position[GroupP,{}]//Flatten[#]&;
PosU1Components=ConstantArray[0,Length[PosU1]];
Do[
PosU1Components[[i-PosU1[[1]]+1]]=Total[SizeGroupsMod[[1;;i]]];
,{i,PosU1}];

ChargeU1F=ConstantArray[0,{Length[FermionRepP],Length[GroupP]}];
Do[
Do[
If[NumericQ[FermionRepP[[j]][[1]][[i]]]==False,
ChargeU1F[[j,i]]=FermionRepP[[j]][[1]][[i]];
FermionRepP[[j]][[1]][[i]]=1;
]
,{i,PosU1}],{j,1,Length[FermionRepP]}];

(********************)


(*Fermion Reps*)
If[Length[FermionRepP]<1,
gvffP=SparseArray[{{1,1,1}->0},{TotalComponents,1,1}]//SparseArray;
,
gvffP=FermionRep[GroupP,FermionRepP[[1]]];
(*Arbitrary U1 charges*)
Do[
If[NumericQ[ChargeU1F[[1]][[i]]]==False,
h=PosU1Components[[i+1-PosU1[[1]]]];
gvffP[[h]]=gvffP[[h]]*ChargeU1F[[1,i]];
]
,{i,PosU1}];
(*********************)

If[Length[FermionRepP]>1,
Do[
Temp=FermionRep[GroupP,FermionRepP[[i]]];

(*Arbitrary U1 charges*)
Do[
If[NumericQ[ChargeU1F[[i,j]]]==False,
h=PosU1Components[[j+1-PosU1[[1]]]];
Temp[[h]]=Temp[[h]]*ChargeU1F[[i,j]];
]
,{j,PosU1}];
(*********************)


gvffTemp=Table[ArrayFlatten[{{gvffP[[a]],0},{0,Temp[[a]]}}],{a,1,Length[gvffP]}];(*Stacks all fermions on the diagonal*)
gvffP=gvffTemp;
,
{i,2,Length[FermionRepP]}]

];

];

CouplingArray=ConstantArray[0,{TotalComponents,TotalComponents}]; (*Only contains coupling constants. No Dynkin index as generators are properly normalized*)
Do[
CouplingArray[[Components[[i]];;Components[[i+1]]-1,Components[[i]];;Components[[i+1]]-1]]=CouplingNameP[[i]]IdentityMatrix[SizeGroupsMod[[i]]];
,
{i,1,Length[GroupP]}
];
gvffTemp=Table[CouplingArray[[a,a]]gvffP[[a]],{a,1,Length[gvffP]}];
gvffP=gvffTemp;


(*Scalar Reps*)
ChargeU1=ConstantArray[0,{Length[ScalarRepP],Length[GroupP]}];
Do[
Do[
If[NumericQ[ScalarRepP[[j]][[1]][[i]]]==False,
ChargeU1[[j,i]]=ScalarRepP[[j]][[1]][[i]];
ScalarRepP[[j]][[1]][[i]]=1;
]
,{i,PosU1}],{j,1,Length[ScalarRepP]}];
(*The above takes care of cases with arbitrary U1 charges*)

If[Length[ScalarRepP]<1,
gvssP=SparseArray[{{1,1,1}->0},{TotalComponents,1,1}]//SparseArray;
,

gvssP=ScalarRepCreation[GroupP,ScalarRepP[[1]]];
(*Arbitrary U1 charges*)
Do[
If[NumericQ[ChargeU1[[1]][[i]]]==False,
h=PosU1Components[[i+1-PosU1[[1]]]];
gvssP[[h]]=gvssP[[h]]*ChargeU1[[1,i]];
]
,{i,PosU1}];
(*********************)
(*********************)

If[Length[ScalarRepP]>1,
Do[
Temp=ScalarRepCreation[GroupP,ScalarRepP[[i]]];

(*********************)
(*Arbitrary U1 charges*)
Do[
If[NumericQ[ChargeU1[[i]][[j]]]==False,
h=PosU1Components[[j+1-PosU1[[1]]]];
Temp[[h]]=Temp[[h]]*ChargeU1[[i,j]];
]
,{j,PosU1}];
(*********************)


gvssTemp=Table[ArrayFlatten[{{gvssP[[a]],0},{0,Temp[[a]]}}],{a,1,Length[gvssP]}]; (*Stack scalar reps on the diagonal*)
gvssP=gvssTemp;
,
{i,2,Length[ScalarRepP]}]

];
];

gvssTemp=Table[CouplingArray[[a,a]]gvssP[[a]],{a,1,Length[gvssP]}];
gvssP=gvssTemp//SparseArray;

(* Creates All tensors*)
(*By default all tensors barring gauge ones are empty*)
nsP=Length[gvssP[[1]]];
nvP=Length[fabc];
nfP=Length[gvffP[[1]]];
\[Lambda]1P=EmptyArray[{nsP}];
\[Lambda]3P=EmptyArray[{nsP,nsP,nsP}];
\[Lambda]4P=EmptyArray[{nsP,nsP,nsP,nsP}];
\[Mu]ijP=EmptyArray[{nsP,nsP}];
\[Mu]IJCP=EmptyArray[{nfP,nfP}];
\[Mu]IJP=EmptyArray[{nfP,nfP}];
YsffP=EmptyArray[{nsP,nfP,nfP}];
YsffCP=EmptyArray[{nsP,nfP,nfP}];
(* Returns*)

(*Creates variables for constructing invariants. Basically naming all components*)
CreateScalarComponents[GroupP,ScalarRepP];
CreateFermionComponents[GroupP,FermionRepP];



GaugeCouplingNames=CouplingNameI;

Return[Join[{fabc//SparseArray//SimplifySparse,gvffP//SparseArray//SimplifySparse,gvssP//SparseArray//SimplifySparse,\[Lambda]1P//SparseArray,\[Lambda]3P//SparseArray,\[Lambda]4P//SparseArray,\[Mu]ijP//SparseArray,\[Mu]IJP//SparseArray,\[Mu]IJCP//SparseArray,YsffP//SparseArray,YsffCP//SparseArray}]];
];


{GaugeCouplingNames};


(*
	Creates representation matrices for scalars. It is possible for the representation to be either real or complex.
*)
ScalarRepCreation[GroupI_,ScalarRepI_]:=Module[{GroupP=GroupI,ScalarRepP=ScalarRepI},
Temp=RepMatrices[GroupP,ScalarRepP[[1]]];
If[ScalarRepP[[2]]=="C",
(*For complex representations we need to rewrite everything in a real basis*)
help=- Table[{{Im[Temp[[a]]],Re[Temp[[a]]]},{-Re[Temp[[a]]],Im[Temp[[a]]]}}//ArrayFlatten,{a,1,Length[Temp]}]//ArrayFlatten;
,
 (*By default real representation matrices are not anti-symmetric. So we need to change basis*)
U=RotationToRealBasis[GroupP,ScalarRepP[[1]]];
help=   I  U . # . ConjugateTranspose[U]&/@Temp;

];
Return[help]
];


GradH[f_,x_List?VectorQ]:=D[f,{x}];


(*
	Creates representation matrices for fermions. It is possible for the representation to be either left- or right-handed.
*)
FermionRep[GroupI_,FermionRepI_]:=Module[{GroupP=GroupI,FermionRepP=FermionRepI},
If[FermionRepP[[2]]=="L",
	help=RepMatrices[GroupP,FermionRepP[[1]]]; (*No change for left-handed fermions*)
,
	help=-Conjugate[RepMatrices[GroupP,FermionRepP[[1]]]]; (*Right-handed matrices are charge-conjugated*)
];
Return[help]
];


(*
	Redunant module. This should be removed.
*)
RotateToAdjBasis[GroupI_,ScalarRepI_]:=Module[{GroupP=GroupI,ScalarRepP=ScalarRepI},
(*This Module maps components from a real basis to the adjoint basis defined by gvss=gvvv*)
(*This mapping is a similarity transformation*)
(*This Module is no longer necessary. Should probably delete it eventually.*)
posGroup=FirstPosition[ScalarRepP, _?(# != 0 &)][[1]]; (*Only simple adjoint reps are allowed*)
OrigRep=RepMatrices[GroupP[[posGroup]],ScalarRepP[[posGroup]]]//ArrayFlatten[#]&//Normal; (*This is the standard rep from GroupMath.*)

PreFacHelp=DynkinIndex[GroupP[[posGroup]],ScalarRepP[[posGroup]]]; (*Normalization for structure constants*)
fabc=-I PreFacHelp^-1 Table[Tr[(b . c-c . b) . a],{a,OrigRep},{b,OrigRep},{c, OrigRep}]//SparseArray;

U=RotationToRealBasis[GroupP,ScalarRepP]; (*Rotation matrix from OrigRep to a real basis*)

(*Choose to find the transformation, for simplicity, via a diagonal matrix.*)
Do[
If[DiagonalMatrixQAE[OrigRep[[i]]]==True,
RefMat=OrigRep[[i]]//Normal;
AdjMat=fabc[[i]]//Normal;
Break[];
];
,{i,1,Length[OrigRep]}];

ListHelp=Array[a[##]&,Length[RefMat]];
U2=DiagonalMatrix[ListHelp]; (**)
E1=U2 . Normalize/@Eigenvectors[AdjMat]; (*The diagonalization matrix is only defined up to a diagonal U(1) matrix*)
E2=Normalize/@Eigenvectors[RefMat];
ETot=ConjugateTranspose[E1] . E2; (*This matrix maps from fabc to OrigRep up to a U(1) tf*)

SolveFac=-I ConjugateTranspose[ETot] . #  . ETot&/@fabc;(*We now choose U2 to ensure that this is a similarity tf*)
solU=FindInstance[Normal[SolveFac]== Normal[OrigRep],ListHelp][[1]];
EFinal=ETot/.solU; (*Complete tf from fabc to OrigRep*)

P=U . ConjugateTranspose[EFinal];(*This is the tf that maps from the real basis to the adjoint basis defined by fabc*)
Return[P]
];


{ScalarVariablesIndices,GaugeIndices};


PrintGaugeRepPositions[]:=Module[{},
Return[GaugeIndices]
];


PrintFermionRepPositions[]:=Module[{},
Return[FermionVariablesIndices]
];


PrintScalarRepPositions[]:=Module[{},
Return[ScalarVariablesIndices]
];


(*
	This Module names all scalar components.
*)
CreateScalarComponents[GroupI_,ScalarRepI_]:=Module[{GroupP=GroupI,ScalarRepP=ScalarRepI},

GroupHelp=DeleteCases[GroupP,{}]; (*Removes all U1 factors since they don't have individual components*)
PosU1=Position[GroupP,{}];(*Position of all U1:s*)
ScalarComponents=ConstantArray[0,{Length[ScalarRepP]}];(*Contains all subsitution rules*)
ScalarComponentsC=ConstantArray[0,{Length[ScalarRepP]}];(*Conjugated components*)
ScalarVariables=ConstantArray[0,{Length[ScalarRepP]}]; (*Names of all components*)
ScalarVariablesIndices=ConstantArray[0,{Length[ScalarRepP]}]; (*Position of every scalar rep*);

Do[
SizeRep=DimR[GroupHelp,ScalarRepP[[i]][[1]]//Delete[#,PosU1]&];
(*By default GroupMath creates invariants with names a[x], b[x], ...*)
If[ScalarRepP[[i]][[2]]=="C",
(*\[Phi] are real components and \[Psi] ar the imagionary ones*)
ScalarComponents[[i]]=Array[a[##]->1/Sqrt[2]Symbol[ToString[\[Phi]]<>ToString[i]][##]+I 1/Sqrt[2] Symbol[ToString[\[Psi]]<>ToString[i]][##]&,SizeRep]//Flatten[#]&;
ScalarComponentsC[[i]]=Array[a[##]->1/Sqrt[2] Symbol[ ToString[\[Phi]]<>ToString[i]][##]-I 1/Sqrt[2] Symbol[ToString[\[Psi]]<>ToString[i]][##]&,SizeRep]//Flatten[#]&;
ScalarVariables[[i]]={#}&/@Flatten[ScalarComponents[[i]]]/.{a_->b_}->b//Variables;
If[i==1,
ScalarVariablesIndices[[i]]=1;;(Length[ScalarVariables[[i]]]);
,
helpInd=ScalarVariablesIndices[[i-1]]/.(a_;;b_)-> b;
ScalarVariablesIndices[[i]]=(helpInd+1);;(helpInd+Length[ScalarVariables[[i]]]);
];
,
(*For a real rep we need to rotate it to a real basis*)
If[ScalarRepP[[i]][[2]]=="R",
RotMat=RotationToRealBasis[GroupP,ScalarRepP[[i]][[1]]]; (*Maps to a real basis*)
Components=Array[Symbol[ToString[\[Phi]]<>ToString[i]][##]&,SizeRep]//Flatten[#]&;
BasisChange=MapThread[Rule,{Components,ConjugateTranspose[RotMat] . Components}];
,
(*This option should be removed*)
RotMat1=RotationToRealBasis[GroupP,ScalarRepP[[i]][[1]]]; (*Maps to a real basis*)
RotMat=RotateToAdjBasis[GroupP,ScalarRepP[[i]][[1]]]; (*Maps to an adjoint basis. Only for simple reps*)
Components=Array[Symbol[ToString[\[Phi]]<>ToString[i]][##]&,SizeRep]//Flatten[#]&;
BasisChange=MapThread[Rule,{Components,ConjugateTranspose[RotMat1] . RotMat . Components}];
];

ScalarComponents[[i]]=Array[a[##]->Symbol[ToString[\[Phi]]<>ToString[i]][##]&,SizeRep]//Flatten[#]&//ReplaceAll[#,BasisChange]&;
ScalarComponentsC[[i]]=Array[a[##]->Symbol[ToString[\[Phi]]<>ToString[i]][##]&,SizeRep]//Flatten[#]&//ReplaceAll[#,BasisChange]&;
ScalarVariables[[i]]={#}&/@Flatten[ScalarComponents[[i]]]/.{a_->b_}->b//Variables;
If[i==1,
ScalarVariablesIndices[[i]]=1;;(Length[ScalarVariables[[i]]]);
,
helpInd=ScalarVariablesIndices[[i-1]]/.(a_;;b_)-> b;
ScalarVariablesIndices[[i]]=(helpInd+1);;(helpInd+Length[ScalarVariables[[i]]]);
];
];
,
{i,1,Length[ScalarRepI]}];


(*Check if the group only contains U1s*)
If[Total[Length[#]&/@GroupP]<1,
(*If the group only contains U1s the components don't have an index. so a[x]-> a. This is fixed below*)
ScalarComponents=ScalarComponents//ReplaceAll[#,a[1]->a]&; 
ScalarComponentsC=ScalarComponentsC//ReplaceAll[#,a[1]->a]&;
];
ScalarVariables=ScalarVariables//Flatten[#]&;
];


{ScalarComponents,ScalarComponentsC,ScalarVariables};


(*
	Names fermion components.
*)
CreateFermionComponents[GroupI_,FermionRepI_]:=Module[{GroupP=GroupI,FermionRepP=FermionRepI},

GroupHelp=DeleteCases[GroupP,{}]; (*Removes all U1 factors since they don't have individual components*)
PosU1=Position[GroupP,{}];(*Position of all U1:s*)

FermionComponents=ConstantArray[0,{Length[FermionRepP]}]; (*Contains all subsitution rules*)
FermionVariables=ConstantArray[0,{Length[FermionRepP]}];(*Contains the name of all fermion variables*)
FermionVariablesIndices=ConstantArray[0,{Length[FermionRepP]}]; (*Position of every scalar rep*);


Do[
SizeRep=DimR[GroupHelp,FermionRepP[[i]][[1]]//Delete[#,PosU1]&];
FermionComponents[[i]]=Array[a[##]->Symbol[ToString[\[CapitalPsi]]<>ToString[i]][##]&,SizeRep]//Flatten[#]&;
FermionVariables[[i]]={#}&/@Flatten[FermionComponents[[i]]]/.{a_->b_}->b//Variables;
If[i==1,
FermionVariablesIndices[[i]]=1;;(Length[FermionVariables[[i]]]);
,
helpInd=FermionVariablesIndices[[i-1]]/.(a_;;b_)-> b;
FermionVariablesIndices[[i]]=(helpInd+1);;(helpInd+Length[FermionVariables[[i]]]);
];
,
{i,1,Length[FermionRepP]}];

(*Check if the group only contains U1s*)
If[Total[Length[#]&/@GroupP]<1,
FermionComponents=FermionComponents//ReplaceAll[#,a[1]->a]&;
];

FermionVariables=FermionVariables//Flatten[#]&;

];


{FermionComponents,FermionVariables,FermionVariablesIndices};


(*
	Creates invariant scalar operators with ordered components.
*)
CreateInvariant[GroupI_,ScalarRepI_,InvariantI_]:=Module[{GroupP=GroupI,ScalarRepP=ScalarRepI,InvariantP=InvariantI},

SubArray=ConstantArray[0,Length[InvariantP[[1]]]];

GroupP=ToGM[GroupP]; (*Converts to GroupMath input*)

(*Fix for non-numeric U1 charges*)
PosU1=Position[GroupP,{}]//Flatten[#]&;

If[Length[Delete[GroupP,{#}&/@PosU1]]>0,
GroupP=Delete[GroupP,{#}&/@PosU1];
Do[
ScalarRepP[[i]][[1]]=Delete[ScalarRepP[[i]][[1]],{#}&/@PosU1];
,{i,1,Length[ScalarRepP]}];
,
Do[
ScalarRepP[[j]][[1]][[i]]=0; (*For arbitrary gauge-charges we set all charges to 0 when creating an invariant. Up to user to ensure that the term is gauge invariant*)
,{j,1,Length[ScalarRepP]},{i,1,Length[GroupP]}]
];
(********************)

ComponentsP={ScalarComponents,ScalarComponentsC,ScalarVariables};
(*This just creates a Sub list*)
j=1;
Do[
If[ScalarRepP[[i]][[2]]=="C",
If[InvariantP[[2]][[j]]==True,
SubArray[[j]]=ComponentsP[[1]][[i]]/.a->ToExpression[Alphabet[][[j]]] 
,
SubArray[[j]]=ComponentsP[[2]][[i]]/.a->ToExpression[Alphabet[][[j]]](*Complex conjugated components*)
];
,
SubArray[[j]]=ComponentsP[[1]][[i]]/.a->ToExpression[Alphabet[][[j]]](*No conjugated components for real reps*)
];
j=j+1;
,
{i,InvariantP[[1]]}];


SubArray=SubArray//Flatten[#]&;
InvRepsP=ScalarRepP[[InvariantP[[1]]]];(*These are the reps specified by the user*)
InvReps=Table[InvRepsP[[i]][[1]],{i,1,Length[InvariantP[[1]]]}];


Temp=Invariants[GroupP,InvReps,Conjugations->InvariantP[[2]]];
InvOperator=ToExpression[StringReplace[ToString[StandardForm[Temp]],"GroupMath`"->""]]//ReplaceAll[#,SubArray]&;
Return[InvOperator]
];


(*
	Creates invariant pure-fermion operators.
*)
CreateInvariantFermion[GroupI_,FermionRepI_,InvariantI_]:=Module[{GroupP=GroupI,FermionRepP=FermionRepI,InvariantP=InvariantI},
(*Here I create an invariant, with ordered components, using the Invariant function from Groupmath*)
SubArray=ConstantArray[0,Length[InvariantP[[1]]]];

GroupP=ToGM[GroupP]; (*Converts to groupmath input*)

(*Fix for non-numeric U1 charges*)
PosU1=Position[GroupP,{}]//Flatten[#]&;
If[Length[Delete[GroupP,{#}&/@PosU1]]>0,
GroupP=Delete[GroupP,{#}&/@PosU1];
Do[
FermionRepP[[i]][[1]]=Delete[FermionRepP[[i]][[1]],{#}&/@PosU1];
,{i,1,Length[FermionRepP]}];
,
Do[
FermionRepP[[j]][[1]][[i]]=0;
,{j,1,Length[FermionRepP]},{i,1,Length[GroupP]}]
];
(********************)

ComponentsP={FermionComponents,FermionVariables};
(*This just creates a Sub list*)
j=1;
Do[
SubArray[[j]]=ComponentsP[[1]][[i]]/.a->ToExpression[Alphabet[][[j]]];
j=j+1;
,
{i,InvariantP[[1]]}];


SubArray=SubArray//Flatten[#]&;
InvRepsP=FermionRepP[[InvariantP[[1]]]];(*These are the reps specified by the user*)
InvReps=Table[InvRepsP[[i]][[1]],{i,1,Length[InvariantP[[1]]]}];



Temp=Invariants[GroupP,InvReps,Conjugations->InvariantP[[2]]];


InvOperator=ToExpression[StringReplace[ToString[StandardForm[Temp]],"GroupMath`"->""]]//ReplaceAll[#,SubArray]&;
Return[InvOperator]
];


(*
	Create Yukawa couplings.
*)
CreateInvariantYukawa[GroupI_,ScalarRepI_,FermionRepI_,InvariantI_]:=Module[{GroupP=GroupI,ScalarRepP=ScalarRepI,FermionRepP=FermionRepI,InvariantP=InvariantI},
(*Here I create an invariant, with ordered components, using the Invariant function from Groupmath*)
(*This just creates a Sub list*)

GroupP=ToGM[GroupP]; (*Converts to GroupMath input*)


ComponentsP={ScalarComponents,ScalarComponentsC,ScalarVariables};
ComponentsFermionP={FermionComponents,FermionVariables};



(*Fix for non-numeric U1 charges*)
PosU1=Position[GroupP,{}]//Flatten[#]&;
If[Length[Delete[GroupP,{#}&/@PosU1]]>0,
GroupP=Delete[GroupP,{#}&/@PosU1];
Do[
ScalarRepP[[i]][[1]]=Delete[ScalarRepP[[i]][[1]],{#}&/@PosU1];
,{i,1,Length[ScalarRepP]}];

Do[
FermionRepP[[i]][[1]]=Delete[FermionRepP[[i]][[1]],{#}&/@PosU1];
,{i,1,Length[FermionRepP]}];


,
Do[
ScalarRepP[[j]][[1]][[i]]=0;
,{j,1,Length[ScalarRepP]},{i,1,Length[GroupP]}];

Do[
FermionRepP[[j]][[1]][[i]]=0;
,{j,1,Length[FermionRepP]},{i,1,Length[GroupP]}];

];



(********************)

j=1;
i=InvariantP[[1]][[1]];
If[ScalarRepP[[i]][[2]]=="C",
If[InvariantP[[2]][[j]]==True,
SubArray=ComponentsP[[1]][[i]]/.a->ToExpression[Alphabet[][[j]]]
,
SubArray=ComponentsP[[2]][[i]]/.a->ToExpression[Alphabet[][[j]]]
];
,
SubArray=ComponentsP[[1]][[i]]/.a->ToExpression[Alphabet[][[j]]]
];


SubArrayFermion=ConstantArray[0,2];

i=InvariantP[[1]][[2]];
SubArrayFermion[[1]]=ComponentsFermionP[[1]][[i]]/.a->ToExpression[Alphabet[][[2]]];
i=InvariantP[[1]][[3]];
SubArrayFermion[[2]]=ComponentsFermionP[[1]][[i]]/.a->ToExpression[Alphabet[][[3]]];


SubArray=SubArray//Flatten[#]&;
SubArrayFermion=SubArrayFermion//Flatten[#]&;


InvRepsScalar=ScalarRepP[[InvariantP[[1]][[1]]]];
(*These are the reps specified by the user*)
InvRepsFemion=FermionRepP[[InvariantP[[1]][[2;;3]]]];


InvReps={InvRepsScalar[[1]],InvRepsFemion[[1]][[1]],InvRepsFemion[[2]][[1]]};
(*//ReplaceAll[#,SubArray]&//ReplaceAll[#,SubArrayFermion]&*)


Temp=Invariants[GroupP,InvReps,Conjugations->InvariantP[[2]]]//ReplaceAll[#,SubArray]&//ReplaceAll[#,SubArrayFermion]&;


InvOperator=ToExpression[StringReplace[ToString[StandardForm[Temp]],"GroupMath`"->""]]//ReplaceAll[#,SubArray]&//ReplaceAll[#,SubArrayFermion]&;
Return[InvOperator];
];


GradMass[PotentialI_]:=Module[{PotentialP=PotentialI},
SparseArray[GradH[PotentialP,ScalarVariables]]//GradH[#,ScalarVariables]&
];


GradMassFermion[PotentialI_]:=Module[{PotentialP=PotentialI},
SparseArray[GradH[PotentialP,FermionVariables]]//GradH[#,FermionVariables]&//SparseArray
];


GradQuartic[PotentialI_]:=Module[{PotentialP=PotentialI},
SparseArray[GradH[PotentialP,ScalarVariables]]//GradH[#,ScalarVariables]&//GradH[#,ScalarVariables]&//GradH[#,ScalarVariables]&//SparseArray
];


GradSextic[PotentialI_]:=Module[{PotentialP=PotentialI},
GradQuartic[PotentialP]//GradH[#,ScalarVariables]&//GradH[#,ScalarVariables]&//SparseArray
];


GradCubic[PotentialI_]:=Module[{PotentialP=PotentialI},
SparseArray[GradH[PotentialP,ScalarVariables]]//GradH[#,ScalarVariables]&//GradH[#,ScalarVariables]&//SparseArray
];


GradTadpole[PotentialI_]:=Module[{PotentialP=PotentialI},
SparseArray[GradH[PotentialP,ScalarVariables]]//SparseArray
];


GradYukawa[PotentialI_]:=Module[{PotentialP=PotentialI},
SparseArray[GradH[PotentialP,ScalarVariables]]//GradH[#,FermionVariables]&//GradH[#,FermionVariables]&//SparseArray
];


(*
	Creates an empty array.
*)
EmptyArray[ComponentsI_]:=Module[{ComponentsP=ComponentsI},
ComponentsMod=ComponentsP;
Do[
If[ComponentsP[[i]]<1,ComponentsMod[[i]]=1,ComponentsMod[[i]]=ComponentsP[[i]]
];
,{i,Length[ComponentsP]}];
SparseArray[{ComponentsMod->0},ComponentsMod]
];


(*
	Finds the dimension of all groups.
*)
ToGM[GroupI_]:=Module[{GroupP=GroupI},
Temp=StringReplace[GroupP,{"SU"->""}];
SUTensor=DeleteCases[GroupP,Alternatives@@Temp];(*Identifies SU groups*)
numTensorSU=CartanMatrix["SU",#]&/@(Read[StringToStream[#],Number]&/@StringReplace[SUTensor,{"SU"->""}]);(*Replace the group with the dimension of the adjoint rep*)

Temp=StringReplace[GroupP,{"SO"->""}];
SUTensor=DeleteCases[GroupP,Alternatives@@Temp];(*Identifies SU groups*)
numTensorSO=CartanMatrix["SO",#]&/@(Read[StringToStream[#],Number]&/@StringReplace[SUTensor,{"SO"->""}]);(*Replace the group with the dimension of the adjoint rep*)

Temp=StringReplace[GroupP,{"SP"->""}];
SUTensor=DeleteCases[GroupP,Alternatives@@Temp];(*Identifies SU groups*)
numTensorSP=CartanMatrix["SP",#]&/@(Read[StringToStream[#],Number]&/@StringReplace[SUTensor,{"SP"->""}]);(*Replace the group with the dimension of the adjoint rep*)

Temp=StringReplace[GroupP,{"E"->""}];
SUTensor=DeleteCases[GroupP,Alternatives@@Temp];(*Identifies SU groups*)
numTensorE=CartanMatrix["E",#]&/@(Read[StringToStream[#],Number]&/@StringReplace[SUTensor,{"E"->""}]);(*Replace the group with the dimension of the adjoint rep*)

Temp=StringReplace[GroupP,{"G"->""}];
SUTensor=DeleteCases[GroupP,Alternatives@@Temp];(*Identifies SU groups*)
numTensorG=CartanMatrix["G",#]&/@(Read[StringToStream[#],Number]&/@StringReplace[SUTensor,{"G"->""}]);(*Replace the group with the dimension of the adjoint rep*)

Temp=StringReplace[GroupP,{"F"->""}];
SUTensor=DeleteCases[GroupP,Alternatives@@Temp];(*Identifies SU groups*)
numTensorF=CartanMatrix["F",#]&/@(Read[StringToStream[#],Number]&/@StringReplace[SUTensor,{"F"->""}]);(*Replace the group with the dimension of the adjoint rep*)

Temp=StringReplace[GroupP,{"U1"->""}];
SUTensor=DeleteCases[GroupP,Alternatives@@Temp];(*Identifies SU groups*)
numTensorU1=CartanMatrix["U",#]&/@(Read[StringToStream[#],Number]&/@StringReplace[SUTensor,{"U"->""}]);(*Replace the group with the dimension of the adjoint rep*)


{SoRes,SuRes,SpRes,G2Res,F4Res,E6Res,E7Res,E8Res,U1Res}//Flatten[#,1]&;
Return[Join[numTensorSO,numTensorSU,numTensorSP,numTensorG,numTensorF,numTensorE,numTensorU1]];
];
