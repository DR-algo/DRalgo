(* ::Package:: *)

Quit[];


If[$InputFileName=="",
	SetDirectory[NotebookDirectory[]],
	SetDirectory[DirectoryName[$InputFileName]]
];

DRalgo`DRalgo`$LoadGroupMath=True;
DRalgo`DRalgo`$GroupMathMultipleModels=True;
<<../Kernel/DRalgo.m


(* ::Chapter:: *)
(*SU2+Higgs+singlet higher-dimensional operators*)


(*
	Regression test for the physical dimension-5 operator \[Phi]^3 W0^2,
	ODIM5["PHYS",3,d]. The cubic couplings of the real singlet together with
	the gauge coupling make the cubic-coupling contractions of this operator
	nonzero.
*)


(* ::Section::Closed:: *)
(*Model*)


Group={"SU2"};
RepAdjoint={{2}};
HiggsDoublet={{{1}},"C"};
RealSinglet={{{0}},"R"};
RepScalar={HiggsDoublet,RealSinglet};
CouplingName={gw};
RepFermion={};


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];

InputInv={{1,1},{True,False}};
MassTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{2,2},{True,True}};
MassTerm2=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
VMass=m2*MassTerm1+m2S/2*MassTerm2;
\[Mu]ij=GradMass[VMass]//Simplify//SparseArray;

InputInv={{1,1,2},{True,False,True}};
CubicTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{2,2,2},{True,True,True}};
CubicTerm2=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
VCubic=\[Mu]HS*CubicTerm1+\[Mu]S*CubicTerm2;
\[Lambda]3=GradCubic[VCubic];

VQuartic=\[Lambda]1H*MassTerm1^2+\[Lambda]S*MassTerm2^2+\[Lambda]HS*MassTerm1*MassTerm2;
\[Lambda]4=GradQuartic[VQuartic];


(* ::Section::Closed:: *)
(*Dimensional Reduction*)


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PerformDRhard[];


(* ::Section:: *)
(*Physical Dimension 5 Matching*)


(*Nonzero components of the physical-basis tensor of \[Phi]^3 W0^2*)
OP3=DeleteCases[Union[Flatten[Simplify[Normal[ODIM5["PHYS",3,3]]]]],0];

(*Powers of Sqrt[T] in each term of these components*)
TPowers=Union[Flatten[
	Exponent[#,sT]&/@Flatten[
		With[{e=Expand[#/.T->sT^2]},If[Head[e]===Plus,List@@e,{e}]]&/@OP3]]];


(* ::Section:: *)
(*TESTS*)


testList={};

(*The cubic-coupling contractions are nonzero in this model*)
AppendTo[testList,
TestCreate[
	FreeQ[OP3,\[Mu]HS]||FreeQ[OP3,gw],
	False
]];

(*Each term scales as T^(3/2)*)
AppendTo[testList,
TestCreate[
	TPowers,
	{3}
]];

AppendTo[testList,
TestCreate[
	Union[Simplify[D[OP3,xi]]],
	{0}
]];

report=TestReport[testList];
Print[report["ResultsDataset"]]
