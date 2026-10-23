(* ::Package:: *)

If[$InputFileName=="",
	SetDirectory[NotebookDirectory[]],
	SetDirectory[DirectoryName[$InputFileName]]
];

DRalgo`DRalgo`$LoadGroupMath=True;
(*DRalgo`DRalgo`$InstallGroupMath=True;*)

Check[
    Get["../Kernel/DRalgo.m"],
    Message[Get::noopen, "DRalgo` at "<>ToString[$UserBaseDirectory]<>"/Applications"];
    Abort[];
]


(* ::Chapter:: *)
(*Options and terminal output*)


(*
	Issue #35: the options CustomMasses and PerturbativeDiagonalization of CalculatePotential and
	CalculatePotentialUS, and FastRotation of RotateTensorsCustomMass, are public symbols with usage
	messages, like the options of ImportModelDRalgo. OptionValue compares the names of options, so
	that they can also be given as strings.

	Issue #50: with wolframscript, the header of DRalgo, the banner of GroupMath, and the model information
	of LoadModelDRalgo were printed as boxes. Without a front end, they are now printed as plain text, and
	DRalgo`DRalgo`$DRalgoBanner=False suppresses the header and the banner of GroupMath.

	Checked here:
	- the options are in the context DRalgo`DRalgo` and have usage messages;
	- PerturbativeDiagonalization->True changes the effective potential of a model with off-diagonal masses,
	  and the option given as a string gives the same potential;
	- CustomMasses->True given as a string gives the same potential as the symbol;
	- loading DRalgo with wolframscript prints a few lines of plain text, and nothing with $DRalgoBanner=False.
*)


testList={};


(* ::Section:: *)
(*Public options*)


AppendTo[testList,
VerificationTest[
	{Context[#],StringQ[MessageName[#,"usage"]]}&/@{CustomMasses,PerturbativeDiagonalization,FastRotation},
	Table[{"DRalgo`DRalgo`",True},{3}],
	TestID->"options-are-public"
]];


(*Two real scalars with a mass mixing m12, so that the mass matrix in the background is not diagonal*)
Group={"U1"};RepAdjoint={0};CouplingName={g1};
RepScalar={{{0},"R"},{{0},"R"}};RepFermion={};
{gvvv,gvff,gvss,lam1,lam3,lam4,muij,muIJ,muIJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
{f1,f2}=Table[CreateInvariant[Group,RepScalar,{{i},{True}}][[1]],{i,2}];
muij=GradMass[m1/2 f1^2+m2/2 f2^2+m12 f1 f2]//SparseArray;
lam4=GradQuartic[l1/24 f1^4+l2/24 f2^4];
ImportModelDRalgo[Group,gvvv,gvff,gvss,lam1,lam3,lam4,muij,muIJ,muIJC,Ysff,YsffC,Verbose->False];
PerformDRhard[];
UseSoftTheory[];
Quiet[DefineVEVS[{v,0}]];

(*With the mass mixing as an insertion, the background is defined again for each potential*)
potential[opts___]:=(Quiet[DefineVEVS[{v,0}]];CalculatePotential[opts];PrintEffectivePotential["NNLO"]);
perturbative=potential[PerturbativeDiagonalization->True];

AppendTo[testList,
VerificationTest[
	{FreeQ[perturbative,m123d],Simplify[potential["PerturbativeDiagonalization"->True]-perturbative]},
	{False,0},
	TestID->"perturbative-diagonalization-as-symbol-and-string"
]];


(*Diagonal mass matrices in the background, given as custom masses; the temporal scalar is massless*)
massS=DiagonalMatrix[{m13d+l13d v^2/2,m23d,0}]//SparseArray;
massV={{0}}//SparseArray;
CalculatePotential[massS,massV,CustomMasses->True];
custom=PrintEffectivePotential["NLO"];
CalculatePotential[massS,massV,"CustomMasses"->True];

AppendTo[testList,
VerificationTest[
	{FreeQ[custom,m123d],Simplify[PrintEffectivePotential["NLO"]-custom],Simplify[custom+(m13d+l13d v^2/2)^(3/2)/(12\[Pi])+m23d^(3/2)/(12\[Pi])]},
	{True,0,0},
	TestID->"custom-masses-as-symbol-and-string"
]];


(* ::Section:: *)
(*Terminal output*)


(*Loads DRalgo with wolframscript in a separate process and returns what it prints*)
loadOutput[banner_]:=RunProcess[{"wolframscript","-code",
	"DRalgo`DRalgo`$DRalgoBanner="<>ToString[banner]<>";Get[\""<>FileNameJoin[{Directory[],"..","Kernel","DRalgo.m"}]<>"\"];Exit[0]"},
	"StandardOutput"];
withBanner=loadOutput[True];
withoutBanner=loadOutput[False];

AppendTo[testList,
VerificationTest[
	{StringContainsQ[withBanner,"DRalgo "<>$DRalgoVersion],StringContainsQ[withBanner,"Please cite GroupMath"],
		StringFreeQ[withBanner,"Graphics"|"TemplateBox"|"Grid["|"Panel["],Length[StringSplit[withBanner,"\n"]]<=6,
		StringTrim[withoutBanner]},
	{True,True,True,True,""},
	TestID->"plain-text-header-in-a-terminal"
]];


(* ::Section:: *)
(*Report*)


report=TestReport[testList];
Print[report];
Print[Column[{#["TestID"],#["Outcome"]}&/@Values[report["TestResults"]]]];
Exit[If[report["AllTestsSucceeded"],0,1]]
