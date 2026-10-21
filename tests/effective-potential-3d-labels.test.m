(* ::Package:: *)

If[$InputFileName=="",
	SetDirectory[NotebookDirectory[]],
	SetDirectory[DirectoryName[$InputFileName]]
];
(*Put this if you want to create multiple model-files with the same kernel*)
(*DRalgo`DRalgo`$GroupMathMultipleModels=True;*)

DRalgo`DRalgo`$LoadGroupMath=True;
(*DRalgo`DRalgo`$InstallGroupMath=True;*)

Check[
    Get["../Kernel/DRalgo.m"],
    Message[Get::noopen, "DRalgo` at "<>ToString[$UserBaseDirectory]<>"/Applications"];
    Abort[];
]


(* ::Chapter:: *)
(*Effective potential of the soft theory: 3d parameters*)


(*
	Issue #57: with UseSoftTheory[], the NNLO effective potential of the SU(2)+Higgs model
	contained the 4d gauge coupling gw next to gw3d. UseSoftTheory builds the couplings of
	the temporal scalars to the spatial vectors from the structure constants, which now
	carry the 3d gauge couplings, as gvvvEP.

	Checked here: at LO, NLO, and NNLO, with and without PerturbativeDiagonalization, the
	effective potential depends only on the parameters of the soft theory.
*)


(* ::Section:: *)
(*Model*)


Group={"SU2"};
RepAdjoint={{2}};
RepScalar={{{{1}},"C"}};
CouplingName={gw};
RepFermion={};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
\[Phi]sq=CreateInvariant[Group,RepScalar,{{1,1},{True,False}}][[1]]//Simplify;
\[Mu]ij=GradMass[msq \[Phi]sq]//Simplify//SparseArray;
\[Lambda]4=GradQuartic[\[Lambda] \[Phi]sq^2];
ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PerformDRhard[];

testList={};


(* ::Section:: *)
(*Effective potential*)


(*The four components of the doublet and the three temporal scalars*)
UseSoftTheory[];
DefineVEVS[{0,0,0,\[Phi]0,0,0,0}];

(*The 4d parameters of the model*)
fourDParametersQ[expr_]:=!FreeQ[expr,gw|msq|\[Lambda]];

Do[
	CalculatePotential[PerturbativeDiagonalization->diag];
	AppendTo[testList,
	VerificationTest[
		fourDParametersQ[PrintEffectivePotential[order]],
		False,
		TestID->"only-3d-parameters-"<>order<>If[diag,"-perturbative-diagonalization",""]
	]];
,{diag,{False,True}},{order,{"LO","NLO","NNLO"}}];

(*The NNLO potential contains the 3d gauge coupling*)
AppendTo[testList,
VerificationTest[
	FreeQ[PrintEffectivePotential["NNLO"],gw3d],
	False,
	TestID->"NNLO-potential-contains-gw3d"
]];


(* ::Section:: *)
(*Report*)


report=TestReport[testList];
Print[report];
Print[Column[{#["TestID"],#["Outcome"]}&/@Values[report["TestResults"]]]];
Exit[If[report["AllTestsSucceeded"],0,1]]
