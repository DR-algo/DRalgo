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
(*Ultrasoft theory without light scalars*)


(*
	If all scalars are integrated out at the soft scale, the ultrasoft theory contains only the gauge fields:
	PerformDRsoft computes the pressure and the gauge couplings, and the ultrasoft scalar masses, tadpoles, and
	couplings are empty, also after an earlier reduction with light scalars.

	Model: SU(2) with a Higgs doublet and a real singlet with the cubic muHS; the US gauge coupling gets the
	one-loop terms of the heavy fields, -g^4 T(R)/(24 Pi m) for a complex representation R of mass m and half of it
	for a real one, from the doublet (T=1/2) and the temporal scalars (adjoint, T=2).
*)


testList={};


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

(*An earlier reduction with the light doublet*)
PerformDRsoft[{5}];
massLight=PrintScalarMassUS["LO"];

AppendTo[testList,
VerificationTest[
	Check[PerformDRsoft[{1,2,3,4,5}];True,False],
	True,
	TestID->"no-messages"
]];

AppendTo[testList,
VerificationTest[
	PrintCouplingsUS[]/.(gw3dUS^2->x_):>(gw3dUS^2->Simplify[x-(gw3d^2-gw3d^4/(24\[Pi]) (1/2/Sqrt[m23d]+1/Sqrt[musqSU2]))]),
	{gw3dUS^2->0},
	TestID->"gauge-coupling"
]];

AppendTo[testList,
VerificationTest[
	{massLight=!={},PrintScalarMassUS["LO"],PrintScalarMassUS["NLO"],PrintTadpolesUS["LO"],PrintIdentificationUS[]},
	{True,{},{},{},{}},
	TestID->"no-scalar-masses-and-tadpoles"
]];

(*Four components of the doublet, the singlet, and three temporal scalars*)
AppendTo[testList,
VerificationTest[
	Simplify[PrintPressureUS["LO"]-(4m23d^(3/2)+m2S3d^(3/2)+3musqSU2^(3/2))/(12\[Pi])],
	0,
	TestID->"LO-pressure"
]];


(* ::Section:: *)
(*Scalar QED with the scalar integrated out*)


(*The abelian gauge couplings are identified from the couplings of the light scalars, so that none remain*)
Group={"U1"};RepAdjoint={0};CouplingName={g1};RepFermion={};
RepScalar={{{1},"C"}};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
phi2=CreateInvariant[Group,RepScalar,{{1,1},{True,False}}][[1]];
\[Mu]ij=GradMass[msq phi2]//SparseArray;
\[Lambda]4=GradQuartic[lam phi2^2];
ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PerformDRhard[];

AppendTo[testList,
VerificationTest[
	Check[PerformDRsoft[{1,2}];{PrintCouplingsUS[],PrintScalarMassUS["LO"]},$Failed],
	{{},{}},
	TestID->"scalar-QED-no-messages"
]];


(* ::Section:: *)
(*Report*)


report=TestReport[testList];
Print[report];
Print[Column[{#["TestID"],#["Outcome"]}&/@Values[report["TestResults"]]]];
Exit[If[report["AllTestsSucceeded"],0,1]]
