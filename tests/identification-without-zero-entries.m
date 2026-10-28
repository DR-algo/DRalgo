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
(*Identification of the soft couplings of tensors without zero entries*)


(*
	ProcessCoupling (HardToSoft.m) names the distinct entries of a tensor, e.g. of the scalar masses or tadpoles,
	by variables that keep the numerical relations between them. The entries are sorted, so that a vanishing entry
	comes first and is dropped. For a tensor without vanishing entries, each entry got the variable of the previous
	entry, which gave wrong masses and tadpoles if entries were related by a number.

	Model: two real singlets f1 and f2 with all mass terms and tadpoles; in Mode->0 the soft masses and tadpoles
	are those of the tree level.
*)


testList={};

Group={"U1"};RepAdjoint={0};CouplingName={g1};RepFermion={};
RepScalar={{{0},"R"},{{0},"R"}};


(* ::Section:: *)
(*Entries related by a number*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
{f1,f2}=Table[CreateInvariant[Group,RepScalar,{{i},{True}}][[1]],{i,2}];
\[Mu]ij=GradMass[m/2 f1^2+(2m)/2 f2^2+m12 f1 f2]//SparseArray;
\[Lambda]1=GradTadpole[t f1+2t f2];
ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False,Mode->0];
PerformDRhard[];

AppendTo[testList,
VerificationTest[
	Sort[PrintScalarMass["LO"]],
	Sort[{m3d->m,m123d->m12}],
	TestID->"masses-related-by-a-number"
]];

AppendTo[testList,
VerificationTest[
	PrintTadpoles["LO"],
	{t3d->t/Sqrt[T]},
	TestID->"tadpoles-related-by-a-number"
]];


(* ::Section:: *)
(*Independent entries*)


(*The one-loop tadpoles from the cubics, 1/2 lam3_ijj times the hard tadpole integral T^2/12 of a massless scalar line, are T^2/24 times the cubics*)
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
\[Mu]ij=GradMass[m1/2 f1^2+m2/2 f2^2+m12 f1 f2]//SparseArray;
\[Lambda]1=GradTadpole[t1 f1+t2 f2];
\[Lambda]3=GradCubic[c1/6 f1^3+c2/2 f1^2 f2+c3/2 f1 f2^2+c4/6 f2^3];
ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
(*Without quartics, the identification of the quartics emits messages, which do not affect the other couplings*)
Quiet[PerformDRhard[]];

AppendTo[testList,
VerificationTest[
	Sort[PrintCouplings[]],
	Sort[{c13d->Sqrt[T]c1,c23d->Sqrt[T]c2,c33d->Sqrt[T]c3,c43d->Sqrt[T]c4}],
	TestID->"independent-cubics"
]];

AppendTo[testList,
VerificationTest[
	Sort[PrintScalarMass["LO"]],
	Sort[{m13d->m1,m23d->m2,m123d->m12}],
	TestID->"independent-masses"
]];

AppendTo[testList,
VerificationTest[
	Simplify[Sort[PrintTadpoles["LO"]][[;;,2]]-{(t1+(c1+c3)T^2/24)/Sqrt[T],(t2+(c2+c4)T^2/24)/Sqrt[T]}],
	{0,0},
	TestID->"independent-tadpoles"
]];


(* ::Section:: *)
(*Report*)


report=TestReport[testList];
Print[report];
Print[Column[{#["TestID"],#["Outcome"]}&/@Values[report["TestResults"]]]];
Exit[If[report["AllTestsSucceeded"],0,1]]
