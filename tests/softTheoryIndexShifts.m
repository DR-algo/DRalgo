(* ::Package:: *)

If[$InputFileName=="",
	SetDirectory[NotebookDirectory[]],
	SetDirectory[DirectoryName[$InputFileName]]
];
(*Several models are created in the same kernel*)
DRalgo`DRalgo`$GroupMathMultipleModels=True;

DRalgo`DRalgo`$LoadGroupMath=True;

Check[
    Get["../Kernel/DRalgo.m"],
    Message[Get::noopen, "DRalgo` at "<>ToString[$UserBaseDirectory]<>"/Applications"];
    Abort[];
]


(* ::Chapter:: *)
(*Temporal scalars in the effective potential of the soft theory*)


(*
	UseSoftTheory (EffPot.m) appends the temporal scalars to the scalars, shifting their indices by the number of
	scalars ns. The couplings of the temporal scalars must be shifted entry by entry, also if a tensor has as many
	nonzero entries as indices, and keep the dimensions of all fields.

	Models: scalar QED with two scalars of charges 1 and 2, where the couplings of two temporal scalars to two
	scalars have four entries, and SU(2) with a doublet and a singlet with the cubic muHS, where the couplings of
	two temporal scalars to the singlet have three entries.
*)


testList={};
private[s_]:=ToExpression["DRalgo`DRalgo`Private`"<>s];


(* ::Section:: *)
(*Scalar QED with charges 1 and 2*)


Group={"U1"};RepAdjoint={0};CouplingName={g1};RepFermion={};
RepScalar={{{1},"C"},{{2},"C"}};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
phi1=CreateInvariant[Group,RepScalar,{{1,1},{True,False}}][[1]];
phi2=CreateInvariant[Group,RepScalar,{{2,2},{True,False}}][[1]];
\[Mu]ij=GradMass[m1 phi1+m2 phi2]//SparseArray;
\[Lambda]4=GradQuartic[l1 phi1^2+l2 phi2^2+l12 phi1 phi2];
ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PerformDRhard[];

AppendTo[testList,
VerificationTest[
	Check[UseSoftTheory[];{Length[Most[ArrayRules[private["lamKVec"]]]],Dimensions[private["lam4EP"]]},$Failed],
	{4,{5,5,5,5}},
	TestID->"QED-no-messages"
]];

(*The temporal scalar is the fifth field*)
AppendTo[testList,
VerificationTest[
	Normal[private["lam4EP"]][[5,5,1;;4,1;;4]]-Normal[private["lamKVec"]][[1,1]]//Simplify//Flatten//Union,
	{0},
	TestID->"QED-temporal-scalar-couplings"
]];


(* ::Section:: *)
(*SU(2) with a doublet and a singlet*)


Group={"SU2"};RepAdjoint={{2}};CouplingName={gw};RepFermion={};
RepScalar={{{{1}},"C"},{{{0}},"R"}};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
phi2=CreateInvariant[Group,RepScalar,{{1,1},{True,False}}][[1]]//Simplify;
s2=CreateInvariant[Group,RepScalar,{{2,2},{True,True}}][[1]]//Simplify;
s1=CreateInvariant[Group,RepScalar,{{2},{True}}][[1]];
\[Mu]ij=GradMass[m2 phi2+m2S/2 s2]//Simplify//SparseArray;
\[Lambda]3=GradCubic[muHS phi2 s1];
\[Lambda]4=GradQuartic[lH phi2^2+lS s2^2+lHS phi2 s2];
ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PerformDRhard[];

AppendTo[testList,
VerificationTest[
	Check[UseSoftTheory[];{Length[Most[ArrayRules[private["lamvvsLS"]]]],Dimensions[private["lam3EP"]]},$Failed],
	{3,{8,8,8}},
	TestID->"SU2-no-messages"
]];

(*The temporal scalars are the fields 6 to 8, and the singlet is the field 5*)
AppendTo[testList,
VerificationTest[
	Normal[private["lam3EP"]][[6;;8,6;;8,1;;5]]-Normal[private["lamvvsLS"]]//Simplify//Flatten//Union,
	{0},
	TestID->"SU2-temporal-scalar-cubics"
]];


(* ::Section:: *)
(*Report*)


report=TestReport[testList];
Print[report];
Print[Column[{#["TestID"],#["Outcome"]}&/@Values[report["TestResults"]]]];
Exit[If[report["AllTestsSucceeded"],0,1]]
