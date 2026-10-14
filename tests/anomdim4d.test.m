(* ::Package:: *)

(*
	AnomDim4D returns the one-loop anomalous dimensions of the 4d fields.

	Before the fix every Switch branch in AnomDim4D ended in ";", so the function
	returned Null for scalars, fermions and vectors.

	The model is the Standard Model with one generation and the top Yukawa coupling.
	The tests check
	- that a value (not Null) is returned,
	- the Landau-gauge values: in DRalgo's convention the anomalous dimension is the
	  1/\[Epsilon] coefficient of the one-loop field renormalisation, Z = 1 + \[Gamma]/\[Epsilon],
	  which is minus the usual \[Gamma] = 1/2 d log Z/d log \[Mu]; the values agree with
	  anom_dim_4d of DRalgo 2.0 for the same model,
	- two relations to the beta functions of BetaFunctions4D[]:
	  the U(1) Ward identity \[Beta](g1^2) = -2 g1^2 \[Gamma]_B, and the running of the Higgs mass
	  parameter at \[Lambda]1H = 0, \[Beta](m1) = -2 \[Gamma]_H m1.
*)


If[$InputFileName=="",
	SetDirectory[NotebookDirectory[]],
	SetDirectory[DirectoryName[$InputFileName]]
];
(*Put this if you want to create multiple model-files with the same kernel*)
(*DRalgo`$GroupMathMultipleModels=True;*)

DRalgo`DRalgo`$LoadGroupMath=True;
(*DRalgo`DRalgo`$InstallGroupMath=True;*)

Check[
    Get["../Kernel/DRalgo.m"],
    Message[Get::noopen, "DRalgo` at "<>ToString[$UserBaseDirectory]<>"/Applications"];
    Abort[];
]


(* ::Section:: *)
(*Model: Standard Model with one generation and the top Yukawa coupling*)


Group={"SU3","SU2","U1"};
RepAdjoint={{1,1},{2},0};
CouplingName={g3,g2,g1};
RepScalar={{{{0,0},{1},1/2},"C"}};
RepFermion={
	{{{1,0},{1},1/6},"L"},
	{{{1,0},{0},2/3},"R"},
	{{{1,0},{0},-1/3},"R"},
	{{{0,0},{1},-1/2},"L"},
	{{{0,0},{0},-1},"R"}
	};


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];


InputInv={{1,1},{True,False}};
MassTerm=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
\[Mu]ij=GradMass[m1*MassTerm]//SparseArray;
\[Lambda]4=GradQuartic[\[Lambda]1H*MassTerm^2];


InputInv={{1,1,2},{False,False,True}};
YukawaTop=CreateInvariantYukawa[Group,RepScalar,RepFermion,InputInv][[1]]//Simplify;
Ysff=-yt1*GradYukawa[YukawaTop];
YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt1>0}]];


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False,Mode->2];


PosScalar=PrintScalarRepPositions[];
PosVector=PrintGaugeRepPositions[];
PosFermion=PrintFermionRepPositions[];
{posH}=PosScalar;
{posGluon,posW,posB}=PosVector;
{posQ,posU,posD,posL,posE}=PosFermion;


tests={};


(* ::Section:: *)
(*A value is returned, not Null*)


AppendTo[tests,VerificationTest[
	FreeQ[{AnomDim4D["S",{posH,posH}],AnomDim4D["V",{posB,posB}],AnomDim4D["F",{posQ,posQ}]},Null],
	True,
	TestID->"anomdim4d-returns-a-value"
]];


(* ::Section:: *)
(*Scalar: Higgs doublet*)


(*Landau gauge: \[Gamma]_H = (3 g2^2 C2(SU2) + 3 g1^2 Y^2 - 2 Y2(S))/(16 \[Pi]^2), with C2 = 3/4, Y = 1/2, 2 Y2(S) = 3 yt1^2*)
AppendTo[tests,VerificationTest[
	AnomDim4D["S",{posH,posH}]//Simplify,
	(9 g2^2+3 g1^2-12 yt1^2)/(64 \[Pi]^2)//Simplify,
	TestID->"anomdim4d-higgs"
]];


(* ::Section:: *)
(*Vectors: gluon, W, and B*)


(*Landau gauge: \[Gamma]_A = g^2 (13/6 C_A - 1/3 T(complex scalars) - 2/3 T(Weyl fermions))/(16 \[Pi]^2)*)
AppendTo[tests,VerificationTest[
	Simplify[{
		AnomDim4D["V",{posGluon,posGluon}],
		AnomDim4D["V",{posW,posW}],
		AnomDim4D["V",{posB,posB}]
		}],
	Simplify[{31 g3^2/(96 \[Pi]^2), 17 g2^2/(96 \[Pi]^2), -43 g1^2/(288 \[Pi]^2)}],
	TestID->"anomdim4d-gauge-fields"
]];


AppendTo[tests,VerificationTest[
	{AnomDim4D["V",{posGluon,posW}],AnomDim4D["V",{posW,posB}],AnomDim4D["V",{posB,posGluon}]},
	{0,0,0},
	TestID->"anomdim4d-gauge-fields-do-not-mix"
]];


(* ::Section:: *)
(*Fermions*)


(*Landau gauge: only the Yukawa coupling contributes, -Y2(F)/(16 \[Pi]^2): yt1^2/2 for the quark doublet, yt1^2 for the up quark*)
AppendTo[tests,VerificationTest[
	Simplify[AnomDim4D["F",{#,#}]&/@{posQ,posU,posD,posL,posE}],
	{-yt1^2/(32 \[Pi]^2), -yt1^2/(16 \[Pi]^2), 0, 0, 0},
	TestID->"anomdim4d-fermions"
]];


(* ::Section:: *)
(*U(1) Ward identity: \[Beta](g1^2) = -2 g1^2 \[Gamma]_B*)


beta4D=BetaFunctions4D[];


AppendTo[tests,VerificationTest[
	Simplify[(g1^2/.beta4D)+2 g1^2 AnomDim4D["V",{posB,posB}]],
	0,
	TestID->"anomdim4d-u1-ward-identity"
]];


(* ::Section:: *)
(*Higgs mass running at \[Lambda]1H = 0: \[Beta](m1) = -2 \[Gamma]_H m1*)


AppendTo[tests,VerificationTest[
	Simplify[(m1/.beta4D/.\[Lambda]1H->0)+2 m1 AnomDim4D["S",{posH,posH}]],
	0,
	TestID->"anomdim4d-higgs-mass-running"
]];


report=TestReport[tests];
Print[report];
Exit[If[report["AllTestsSucceeded"],0,1]]
