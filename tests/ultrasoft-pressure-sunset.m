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
(*Ultrasoft pressure: scalar sunset with cubic couplings*)


(*
	Issue #66: the scalar sunset of the soft theory,
	-1/12 \[Lambda]3^2 (1/(4\[Epsilon])+Log[\[Mu]3/(m1+m2+m3)]+1/2)/(16\[Pi]^2) in the free energy,
	contributes +T \[Lambda]3^2 Log[\[Mu]3]/(192\[Pi]^2) to the pressure, which the running from \[Mu] to \[Mu]3
	in PrintPressure["NNLO"] cancels. When heavy scalars are integrated out, the sunsets with at least
	one heavy line belong to the ultrasoft pressure, with the light masses set to zero, while the sunset
	with only light lines belongs to the ultrasoft theory.

	Checked here:
	- one light and one heavy real scalar with all cubic couplings: the sunsets in PrintPressureUS["NLO"],
	  and the \[Mu]3-dependence of the hard and ultrasoft pressure, which is that of the light sunset only;
	- the scalar of examples/ScalYukawa.m integrated out, the example of issue #66: the \[Mu]3-dependence
	  of the hard pressure is cancelled.
*)


testList={};

Group={"U1"};RepAdjoint={0};CouplingName={g1};RepFermion={};


(* ::Section:: *)
(*One light and one heavy real scalar with all cubic couplings*)


(*The field f is light and s is heavy*)
RepScalar={{{0},"R"},{{0},"R"}};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
fField=CreateInvariant[Group,RepScalar,{{1},{True}}][[1]];
sField=CreateInvariant[Group,RepScalar,{{2},{True}}][[1]];
\[Mu]ij=GradMass[mf/2 fField^2+ms/2 sField^2]//SparseArray;
\[Lambda]3=GradCubic[cf/6 fField^3+cs/6 sField^3+ya/2 fField^2 sField+yb/2 fField sField^2];
\[Lambda]4=GradQuartic[lf/24 fField^4+ls/24 sField^4+l22/4 fField^2 sField^2];
ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
Quiet[PerformDRhard[]];
PerformDRsoft[{2}];

(*Sunsets with the lines sss, ffs, and fss, where the light line f is massless*)
sunsets=1/(192\[Pi]^2)(
	cs3d^2 (1/2+Log[\[Mu]3/(3Sqrt[ms3d])])
	+3ya3d^2 (1/2+Log[\[Mu]3/Sqrt[ms3d]])
	+3yb3d^2 (1/2+Log[\[Mu]3/(2Sqrt[ms3d])]));
pressureUS=PrintPressureUS["NLO"];

AppendTo[testList,
VerificationTest[
	Simplify[pressureUS-(pressureUS/.{cs3d->0,ya3d->0,yb3d->0})-sunsets],
	0,
	TestID->"light-and-heavy-scalar-sunsets"
]];

(*The 3d couplings at LO; the sunset of the light cubic cf remains for the ultrasoft theory*)
lo3d={cf3d->Sqrt[T]cf,cs3d->Sqrt[T]cs,ya3d->Sqrt[T]ya,yb3d->Sqrt[T]yb};
AppendTo[testList,
VerificationTest[
	Simplify[\[Mu]3 D[PrintPressure["NNLO"]+T (PrintPressureUS["LO"]+pressureUS),\[Mu]3]/.lo3d],
	-T^2 cf^2/(192\[Pi]^2),
	TestID->"light-and-heavy-scalar-mu3-dependence"
]];


(* ::Section:: *)
(*Scalar Yukawa model with the scalar integrated out*)


(*The model of examples/ScalYukawa.m, see tests/ScalYukawa-rg-invariance.m*)
RepScalar={{{0},"R"}};
RepFermion={{{0},"L"},{{0},"R"}};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJC,\[Mu]IJ,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];

\[Phi]=CreateInvariant[Group,RepScalar,{{1},{True}}][[1]];
\[Lambda]1=GradTadpole[\[Sigma] \[Phi]];
\[Mu]ij=GradMass[msq/2 \[Phi]^2]//SparseArray;
\[Lambda]3=GradCubic[\[Gamma]/6 \[Phi]^3];
\[Lambda]4=GradQuartic[\[Lambda]/24 \[Phi]^4];
\[Mu]IJ=m\[Psi] GradMassFermion[CreateInvariantFermion[Group,RepFermion,{{2,1},{False,True}}][[1]]];
\[Mu]IJC=m\[Psi] GradMassFermion[CreateInvariantFermion[Group,RepFermion,{{1,2},{False,True}}][[1]]];
Ysff=y GradYukawa[CreateInvariantYukawa[Group,RepScalar,RepFermion,{{1,2,1},{True,False,True}}][[1]]];
YsffC=y GradYukawa[CreateInvariantYukawa[Group,RepScalar,RepFermion,{{1,1,2},{True,False,True}}][[1]]];

ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PerformDRhard[];
(*Without light scalars, the matching of the light sector emits messages, which do not affect the pressure*)
Quiet[PerformDRsoft[{1}]];

AppendTo[testList,
VerificationTest[
	Simplify[\[Mu]3 D[PrintPressure["NNLO"]+T (PrintPressureUS["LO"]+PrintPressureUS["NLO"]),\[Mu]3]/.\[Gamma]3d->Sqrt[T]\[Gamma]],
	0,
	TestID->"ScalYukawa-mu3-dependence-cancels"
]];


(* ::Section:: *)
(*Report*)


report=TestReport[testList];
Print[report];
Print[Column[{#["TestID"],#["Outcome"]}&/@Values[report["TestResults"]]]];
Exit[If[report["AllTestsSucceeded"],0,1]]
