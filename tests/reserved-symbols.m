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
(*Couplings with the names of symbols in the results*)


(*
	Issue #32: a cubic coupling \[Mu]3, as in examples/dark-su2-higgs-singlet.m and examples/xsm.m, is the
	same symbol as the 3d renormalization scale \[Mu]3 in the results. ImportModelDRalgo warns if a coupling
	has the name of the temperature T, of a renormalization scale \[Mu], \[Mu]3, or \[Mu]3US, of Lb, Lf, \[Epsilon],
	of the gauge parameter xi, or of a Debye mass.

	Checked here:
	- a real scalar with a cubic coupling \[Mu]3 or Lb gives the warning, with c3 not;
	- SU(2) with a Higgs doublet whose mass is called \[Mu]sqSU2, the name of the Debye mass, gives the warning;
	- the examples xsm.m and dark-su2-higgs-singlet.m, whose cubic coupling is now c3, give no warning,
	  nor LRSymmetric.m, whose mass of the triplets is now \[Mu]\[CapitalDelta].
*)


testList={};

realScalar[cubic_]:=Module[{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,\[Phi]},
	{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=
		AllocateTensors[{"U1"},{0},{g1},{},{{{0},"R"}}];
	\[Phi]=CreateInvariant[{"U1"},{{{0},"R"}},{{1},{True}}][[1]];
	\[Mu]ij=GradMass[msq/2 \[Phi]^2]//SparseArray;
	\[Lambda]3=GradCubic[cubic/6 \[Phi]^3];
	\[Lambda]4=GradQuartic[\[Lambda]/24 \[Phi]^4];
	ImportModelDRalgo[{"U1"},gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False]
];

AppendTo[testList,
VerificationTest[
	realScalar[\[Mu]3],
	Null,
	{ImportModelDRalgo::reserved},
	TestID->"cubic-coupling-mu3-warns"
]];

AppendTo[testList,
VerificationTest[
	realScalar[Lb],
	Null,
	{ImportModelDRalgo::reserved},
	TestID->"cubic-coupling-Lb-warns"
]];

AppendTo[testList,
VerificationTest[
	realScalar[c3],
	Null,
	TestID->"cubic-coupling-c3-does-not-warn"
]];

higgsDoublet[mass_]:=Module[{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,\[Phi]sq},
	{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=
		AllocateTensors[{"SU2"},{{2}},{gw},{},{{{1},"C"}}];
	\[Phi]sq=CreateInvariant[{"SU2"},{{{1},"C"}},{{1,1},{True,False}}][[1]];
	\[Mu]ij=GradMass[mass \[Phi]sq]//SparseArray;
	\[Lambda]4=GradQuartic[\[Lambda] \[Phi]sq^2];
	ImportModelDRalgo[{"SU2"},gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False]
];

AppendTo[testList,
VerificationTest[
	higgsDoublet[\[Mu]sqSU2],
	Null,
	{ImportModelDRalgo::reserved},
	TestID->"mass-with-the-name-of-the-Debye-mass-warns"
]];

AppendTo[testList,
VerificationTest[
	higgsDoublet[msq],
	Null,
	TestID->"higgs-doublet-does-not-warn"
]];


(* ::Section:: *)
(*Examples*)


(*The model part of an example, up to and including the call of ImportModelDRalgo*)
exampleModel[file_]:=Module[{text,start,stop},
	text=Import[FileNameJoin[{"..","examples",file}],"Text"];
	text=StringReplace[text,{"SetDirectory[NotebookDirectory[]];"->"","<<../Kernel/DRalgo.m"->"","<<DRalgo`DRalgo`"->""}];
	start=StringPosition[text,"ImportModelDRalgo["][[1,1]];
	stop=SelectFirst[StringPosition[text,"];"],#[[1]]>start&][[2]];
	StringTake[text,stop]
];

AppendTo[testList,
VerificationTest[
	ToExpression[exampleModel["xsm.m"]],
	Null,
	TestID->"example-xsm-does-not-warn"
]];

AppendTo[testList,
VerificationTest[
	ToExpression[exampleModel["dark-su2-higgs-singlet.m"]],
	Null,
	TestID->"example-dark-su2-higgs-singlet-does-not-warn"
]];

AppendTo[testList,
VerificationTest[
	ToExpression[exampleModel["LRSymmetric.m"]],
	Null,
	TestID->"example-LRSymmetric-does-not-warn"
]];


(* ::Section:: *)
(*Report*)


report=TestReport[testList];
Print[report];
Print[Column[{#["TestID"],#["Outcome"]}&/@Values[report["TestResults"]]]];
Exit[If[report["AllTestsSucceeded"],0,1]]
