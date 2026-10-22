(* ::Package:: *)

If[$InputFileName=="",
	SetDirectory[NotebookDirectory[]],
	SetDirectory[DirectoryName[$InputFileName]]
];
(*Several models are created in the same kernel*)
DRalgo`DRalgo`$GroupMathMultipleModels=True;

DRalgo`DRalgo`$LoadGroupMath=True;
(*DRalgo`DRalgo`$InstallGroupMath=True;*)

Check[
    Get["../Kernel/DRalgo.m"],
    Message[Get::noopen, "DRalgo` at "<>ToString[$UserBaseDirectory]<>"/Applications"];
    Abort[];
]


(* ::Chapter:: *)
(*Gauge invariance of the Yukawa couplings of the examples*)


(*
	DRalgo converts right-handed fermions to left-handed ones with the conjugate generators, and Ysff
	multiplies the conjugates of these fields. In CreateInvariantYukawa, left-handed fermions enter
	conjugated and right-handed fermions not, as in sm.m with {False,False,True} for (\[Phi],QL,tR).
	Then gvss[c].Ysff+I(gvff[c].Ysff+(gvff[c].Ysff)^T)=0 for every generator c, where gvff[c] acts on the
	first fermion index and the transpose is in the fermion indices.

	For the opposite flags, Ysff contains the conjugate invariants, which are not invariant under gvff.
	The fermion loops in the one-loop couplings of the temporal scalars to the scalars then break the
	gauge symmetry, and PrintDebyeMass fails, since the Debye masses of one group run differently.

	Checked here:
	- the Yukawa couplings of sm.m, SU5.m, and LRSymmetric.m are invariant, the conjugate ones of LRSymmetric.m not;
	- the one-loop couplings of the temporal scalars to the scalars of LRSymmetric.m are gauge invariant,
	  and PrintDebyeMass gives one Debye mass per group without messages.
*)


testList={};

private[name_String]:=ToExpression["DRalgo`DRalgo`Private`"<>name];

(*The model part of an example, up to and including the call of ImportModelDRalgo*)
exampleModel[file_]:=Module[{text,start,stop},
	text=Import[FileNameJoin[{"..","examples",file}],"Text"];
	text=StringReplace[text,{"SetDirectory[NotebookDirectory[]];"->"","<<../Kernel/DRalgo.m"->"","<<DRalgo`DRalgo`"->""}];
	start=StringPosition[text,"ImportModelDRalgo["][[1,1]];
	stop=SelectFirst[StringPosition[text,"];"],#[[1]]>start&][[2]];
	StringTake[text,stop]
];

(*Variation of the Yukawa couplings y_iIJ under the generator c*)
yukawaVariation[y_,c_]:=Module[{yN=Normal[y],fermionY},
	fermionY=Transpose[Normal[gvff][[c]] . Transpose[yN,{2,1,3}],{2,1,3}];
	Normal[gvss][[c]] . yN+I (fermionY+Transpose[fermionY,{1,3,2}])
];
yukawaResidual[y_]:=Union@Flatten@Simplify@Table[yukawaVariation[y,c],{c,Length[gvss]}];

(*Variation of the couplings G_abij of two temporal scalars to two scalars under the generator c*)
temporalVariation[g_,c_]:=Module[{gN=Normal[g],f=Normal[gvvv][[c]],s=Normal[gvss][[c]]},
	f . gN+Transpose[f . Transpose[gN,{2,1,3,4}],{2,1,3,4}]
	+Transpose[s . Transpose[gN,{2,3,1,4}],{3,1,2,4}]+Transpose[s . Transpose[gN,{2,3,4,1}],{4,1,2,3}]
];


(* ::Section:: *)
(*Standard Model*)


ToExpression[exampleModel["sm.m"]];

AppendTo[testList,
VerificationTest[
	yukawaResidual[Ysff],
	{0},
	TestID->"sm-yukawa-couplings-are-invariant"
]];


(* ::Section:: *)
(*SU(5)*)


ToExpression[exampleModel["SU5.m"]];

AppendTo[testList,
VerificationTest[
	yukawaResidual[Ysff],
	{0},
	TestID->"SU5-yukawa-couplings-are-invariant"
]];


(* ::Section:: *)
(*Left-right symmetric model*)


ToExpression[exampleModel["LRSymmetric.m"]];

AppendTo[testList,
VerificationTest[
	yukawaResidual[Ysff],
	{0},
	TestID->"LRSymmetric-yukawa-couplings-are-invariant"
]];

AppendTo[testList,
VerificationTest[
	yukawaResidual[YsffC]==={0},
	False,
	TestID->"LRSymmetric-conjugate-yukawa-couplings-are-not-invariant"
]];

PerformDRhard[];

AppendTo[testList,
VerificationTest[
	Union@Flatten@Simplify@Table[temporalVariation[private["GvvssL"],c],{c,Length[gvss]}],
	{0},
	TestID->"LRSymmetric-temporal-scalar-couplings-are-gauge-invariant"
]];

AppendTo[testList,
VerificationTest[
	AllTrue[PrintDebyeMass["NLO"][[All,2]],FreeQ[#,_Rule|_List]&],
	True,
	TestID->"LRSymmetric-debye-masses"
]];


(* ::Section:: *)
(*Report*)


report=TestReport[testList];
Print[report];
Print[Column[{#["TestID"],#["Outcome"]}&/@Values[report["TestResults"]]]];
Exit[If[report["AllTestsSucceeded"],0,1]]
