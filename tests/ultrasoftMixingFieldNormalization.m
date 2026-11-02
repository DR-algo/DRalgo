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
(*Light-heavy mixing in the ultrasoft quartics*)


(*
	A light real singlet f and a heavy real singlet s with the mixing m12 f s; s is integrated out at the soft scale,
	and the results are of first order in m12.

	At tree level, s=-(m12 f+l31 f^3/6)/M with M=ms3d, so that V ⊃ -(m12 f+l31 f^3/6)^2/(2M) and the quartic of f gets
	-4 m12 l31/M, also without cubic couplings.

	The light two-point function contains the exchange of s between the mixing m12 and the one-loop light-heavy
	self-energy Pi(p), -2 m12 Pi(p)/(p^2+M), whose p^2 term 2 m12 (Pi(0)/M^2-Pi'/M) enters the field normalization Z, and
	the quartic l4 of f gets -2(Z-1) l4. With the tadpole integral -Sqrt[M]/(4Pi) and the bubble 1/(8Pi Sqrt[M]) of two
	lines of mass Sqrt[M], whose p^2 term is -1/(96Pi M^(3/2)):
	- with the quartic l13/6 f s^3, Pi=-l13 Sqrt[M]/(8Pi), so that the quartic of f gets m12 l13 l4/(2Pi M^(3/2));
	- with the cubics cm/2 f s^2 and cs/6 s^3, Pi(0)=-cm cs/(16Pi Sqrt[M])+cm cs/(8Pi Sqrt[M]) from the bubble and from
	  the shift of s by its one-loop tadpole -cs Sqrt[M]/(8Pi), and Pi'=cm cs/(192Pi M^(3/2)), so that the quartic of f
	  gets -11 m12 cm cs l4/(48Pi M^(5/2)).
*)


testList={};

Group={"U1"};RepAdjoint={0};CouplingName={g1};RepFermion={};
RepScalar={{{0},"R"},{{0},"R"}};


(* ::Section:: *)
(*Quartic couplings only*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
{f1,s1}=Table[CreateInvariant[Group,RepScalar,{{i},{True}}][[1]],{i,2}];
\[Mu]ij=GradMass[mf/2 f1^2+ms/2 s1^2+m12 f1 s1]//SparseArray;
\[Lambda]4=GradQuartic[l4/24 f1^4+l31/6 f1^3 s1+l13/6 f1 s1^3+l22/4 f1^2 s1^2+ls/24 s1^4];
ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False,Mode->1];
PerformDRhard[];
PerformDRsoft[{2}];
quartic=Expand[l43dUS/.PrintCouplingsUS[]];

AppendTo[testList,
VerificationTest[
	Coefficient[quartic/.{l223d->0,l133d->0,ls3d->0,l43d->0},m123d l313d],
	-4/ms3d,
	TestID->"tree-level-mixing-without-cubics"
]];

AppendTo[testList,
VerificationTest[
	Simplify[Coefficient[quartic,m123d l133d l43d]],
	1/(2\[Pi] ms3d^(3/2)),
	TestID->"field-normalization-from-a-quartic"
]];


(* ::Section:: *)
(*Cubic couplings*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
\[Mu]ij=GradMass[mf/2 f1^2+ms/2 s1^2+m12 f1 s1]//SparseArray;
\[Lambda]3=GradCubic[cm/2 f1 s1^2+cs/6 s1^3];
\[Lambda]4=GradQuartic[l4/24 f1^4];
ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False,Mode->1];
PerformDRhard[];
PerformDRsoft[{2}];
quartic=Expand[l43dUS/.PrintCouplingsUS[]];

AppendTo[testList,
VerificationTest[
	Simplify[Coefficient[quartic,m123d cm3d cs3d l43d]],
	-11/(48\[Pi] ms3d^(5/2)),
	TestID->"field-normalization-from-cubics"
]];


(* ::Section:: *)
(*Report*)


report=TestReport[testList];
Print[report];
Print[Column[{#["TestID"],#["Outcome"]}&/@Values[report["TestResults"]]]];
Exit[If[report["AllTestsSucceeded"],0,1]]
