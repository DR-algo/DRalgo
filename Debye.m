(* ::Package:: *)

(* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ *)

(* :Title: Debye															*)

(*
	This software is covered by the GNU General Public License 3.
	Copyright (C) 2021-2022 Andreas Ekstedt
	Copyright (C) 2021-2022 Philipp Schicho
	Copyright (C) 2021-2022 Tuomas V.I. Tenkanen

*)

(* :Summary:	Names Debye masses for the dimensionally reduced theory.	*)

(* ------------------------------------------------------------------------ *)


(*
	This module takes a list of groups and names Debye masses for each group.
	If there are multiple groups with the same names, the masses are indexed.
	In addition, for such cases there are additional masses created to 
	represent mixing between the groups.

*)

CreateDebyeMasses[GroupI_]:=Module[{GroupP=GroupI},

Components=ComponentsOfGroup[GroupP];(*The number of components for each rep*)
Indexing=IndexDuplicateGroups[GroupP]; (*Position of groups with multiple factors*)
If[Length[GroupP]<1,
Return[{}] (*No mass is created if there are no group*)
,
Debye=StringJoin["\[Mu]sq",Indexing[[1]]]IdentityMatrix[Components[[1]]]; (*names each mass \[Mu]sq+"group name"*)
If[Length[GroupP]>1,
Do[
Temp=StringJoin["\[Mu]sq",Indexing[[i]]]IdentityMatrix[Components[[i]]];
DebyeTemp=ArrayFlatten[{{Debye,0},{0,Temp}}]; (*Stack the masses on the diagonal*)
Debye=DebyeTemp;
,
{i,2,Length[GroupP]}]
];
];

(*

	Create off-diagonal masses. The convention is that for the group {"U1","U1","U1"}, \[Mu]sqU1Mix1 is the mixing between the first and second group.
	\[Mu]sqU1Mix2 is the mixing between the first and third group. And \[Mu]sqU1Mix3 is the mixing between the second and third group etc.
*)
PosDup=Select[positionDuplicates[GroupP],Length[#]>1&];
Do[
k=1;
nameHelp=Reps[[l]][[1]];
PosDupHelp=PosDup[[l]];
TempList=Drop[PosDupHelp,-1];
TempList2=Drop[PosDupHelp,1];
MixVec={};
Do[
Do[
NameTemp=StringJoin[StringJoin[nameHelp,"Mix"],ToString[k]];
Size=Components[[i]];
(*Determines where in the mass matrix each group should go*)
HorInd=(Components[[1;;j]]//Total[#]&)-Size+1;
VertInd=(Components[[1;;i]]//Total[#]&)-Size+1;
VecTemp={NameTemp,HorInd,VertInd,Size};
AppendTo[MixVec,VecTemp];
k=k+1;
,{j,TempList2}];
TempList2=Drop[TempList2,1];
,{i,TempList}];
Do[
MixVecHelp=MixVec[[n]];
MatTemp=StringJoin["\[Mu]",MixVecHelp[[1]]]IdentityMatrix[MixVecHelp[[4]]];

Debye[[MixVecHelp[[2]];;MixVecHelp[[2]]+MixVecHelp[[4]]-1,MixVecHelp[[3]];;MixVecHelp[[3]]+MixVecHelp[[4]]-1]]=MatTemp;
Debye[[MixVecHelp[[3]];;MixVecHelp[[3]]+MixVecHelp[[4]]-1,MixVecHelp[[2]];;MixVecHelp[[2]]+MixVecHelp[[4]]-1]]=MatTemp;
,{n,1,Length[MixVec]}];

,{l,1,Length[PosDup]}];
Return[SparseArray[ToExpression[Debye]]]
];


(*
	Old routine.This routine should be removed.
*)
CreateDebyeMassesOld[GroupI_]:=Module[{GroupP=GroupI},
Components=ComponentsOfGroup[GroupP];(*The number of components for each rep*)
Indexing=IndexDuplicateGroups[GroupP];
If[Length[GroupP]<1,
Return[{}]
,
Debye=StringJoin["\[Mu]",Indexing[[1]]]IdentityMatrix[Components[[1]]];
If[Length[GroupP]>1,
Do[
Temp=StringJoin["\[Mu]",Indexing[[i]]]IdentityMatrix[Components[[i]]];
DebyeTemp=ArrayFlatten[{{Debye,0},{0,Temp}}];
Debye=DebyeTemp;
,
{i,2,Length[GroupP]}]

];
];
Return[SparseArray[ToExpression[Debye]]]
];


(*
	Counts the number of components for the adjoint reps. That is SUReps["SUN"] would correspond to N^2-1.
*)
ComponentsOfGroup[GroupI_]:=Module[{GroupP=GroupI},
SpRes=SPReps[GroupP];
SuRes=SUReps[GroupP];
SoRes=SOReps[GroupP];
G2Res=G2Reps[GroupP];
F4Res=F4Reps[GroupP];
E6Res=E6Reps[GroupP];
E7Res=E7Reps[GroupP];
E8Res=E8Reps[GroupP];
U1Res=U1Reps[GroupP];
TotRes={SoRes,SuRes,SpRes,G2Res,F4Res,E6Res,E7Res,E8Res,U1Res}//Flatten[#,1]&;
Return[TotRes]
];


U1Reps[GroupI_]:=Module[{GroupP=GroupI},
Temp=StringReplace[GroupP,{"U1"->""}];
U1Tensor=DeleteCases[GroupP,Alternatives@@Temp];(*Identifies U1 groups*)
numTensor=1&/@(Read[StringToStream[#],Number]&/@StringReplace[U1Tensor,{"U1"->""}]);(*Replace the group with the dimension of the adjoint rep*)
Return[numTensor]
];


E8Reps[GroupI_]:=Module[{GroupP=GroupI},
Temp=StringReplace[GroupP,{"E8"->""}];
E8Tensor=DeleteCases[GroupP,Alternatives@@Temp];(*Identifies E8 groups*)
numTensor=248&/@(Read[StringToStream[#],Number]&/@StringReplace[E8Tensor,{"E8"->""}]);(*Replace the group with the dimension of the adjoint rep*)
Return[numTensor]
];


E7Reps[GroupI_]:=Module[{GroupP=GroupI},
Temp=StringReplace[GroupP,{"E7"->""}];
E7Tensor=DeleteCases[GroupP,Alternatives@@Temp];(*Identifies E7 groups*)
numTensor=133&/@(Read[StringToStream[#],Number]&/@StringReplace[E7Tensor,{"E7"->""}]);(*Replace the group with the dimension of the adjoint rep*)
Return[numTensor]
];


E6Reps[GroupI_]:=Module[{GroupP=GroupI},
Temp=StringReplace[GroupP,{"E6"->""}];
E6Tensor=DeleteCases[GroupP,Alternatives@@Temp];(*Identifies E6 groups*)
numTensor=78&/@(Read[StringToStream[#],Number]&/@StringReplace[E6Tensor,{"E6"->""}]);(*Replace the group with the dimension of the adjoint rep*)
Return[numTensor]
];


F4Reps[GroupI_]:=Module[{GroupP=GroupI},
Temp=StringReplace[GroupP,{"F4"->""}];
F4Tensor=DeleteCases[GroupP,Alternatives@@Temp];(*Identifies F4 groups*)
numTensor=52&/@(Read[StringToStream[#],Number]&/@StringReplace[F4Tensor,{"F4"->""}]);(*Replace the group with the dimension of the adjoint rep*)
Return[numTensor]
];


G2Reps[GroupI_]:=Module[{GroupP=GroupI},
Temp=StringReplace[GroupP,{"G2"->""}];
G2Tensor=DeleteCases[GroupP,Alternatives@@Temp];(*Identifies G2 groups*)
numTensor=14&/@(Read[StringToStream[#],Number]&/@StringReplace[G2Tensor,{"G2"->""}]);(*Replace the group with the dimension of the adjoint rep*)
Return[numTensor]
];


SPReps[GroupI_]:=Module[{GroupP=GroupI},
Temp=StringReplace[GroupP,{"SP"->""}];
SPTensor=DeleteCases[GroupP,Alternatives@@Temp];(*Identifies SP groups*)
numTensor=SPAdjoint[#]&/@(Read[StringToStream[#],Number]&/@StringReplace[SPTensor,{"SP"->""}]);(*Replace the group with the dimension of the adjoint rep*)
numTensorS=ToString[#]&/@numTensor;(*Converts to string*)
Return[numTensor]
];


SPAdjoint[n_]:=1/2n(n+1);


SUReps[GroupI_]:=Module[{GroupP=GroupI},
Temp=StringReplace[GroupP,{"SU"->""}];
SUTensor=DeleteCases[GroupP,Alternatives@@Temp];(*Identifies SU groups*)
numTensor=SUAdjoint[#]&/@(Read[StringToStream[#],Number]&/@StringReplace[SUTensor,{"SU"->""}]);(*Replace the group with the dimension of the adjoint rep*)
numTensorS=ToString[#]&/@numTensor;(*Converts to string*)
Return[numTensor]
];


SUAdjoint[n_]:=n^2-1;


SOReps[GroupI_]:=Module[{GroupP=GroupI},
Temp=StringReplace[GroupP,{"SO"->""}];
SOTensor=DeleteCases[GroupP,Alternatives@@Temp];(*Identifies SO groups*)
numTensor=SOAdjoint[#]&/@(Read[StringToStream[#],Number]&/@StringReplace[SOTensor,{"SO"->""}]);(*Replace the group with the dimension of the adjoint rep*)
numTensorS=ToString[#]&/@numTensor;(*Converts to string*)
Return[numTensor]
];


SOAdjoint[n_]:=n (n-1)/2;


positionDuplicates[list_] := GatherBy[Range@Length[list], list[[#]] &]


IndexDuplicateGroups[GroupI_]:=Module[{GroupP=GroupI},
Duplicates=Select[positionDuplicates[GroupP],Length[#]>1&];(*Selects all reccuring groups*)
Reps=Table[GroupP[[Duplicates[[i]]]],{i,Length[Duplicates]}];
RepsRenamed=Table[IndexDuplicates[Reps[[i]]],{i,Length[Duplicates]}];
ResTemp=GroupP;
Do[
Temp=RepsRenamed[[i,;;]];
ResTemp[[Duplicates[[i]]]]=Temp;
,{i,Length[Duplicates]}];
Return[ResTemp]
];


IndexDuplicates[GroupI_]:=Module[{GroupP=GroupI},
Dup=ToString[#]&/@(Table[i,{i,Length[GroupP]}]);
DupFinal=StringJoin["P",#]&/@Dup;
Return[Table[StringJoin[GroupP[[i]],DupFinal[[i]]],{i,Length[GroupP]}]]
];
