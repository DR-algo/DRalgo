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
(*LO pressure: sectors that couple to nothing*)


(*
	The interface needs a gauge group and allocates at least one scalar, so a model
	without gauge interactions gets a U(1) under which every field is neutral, and a
	model without scalars gets a placeholder scalar. SymmetricPhaseLO drops such a
	sector from the ideal-gas pressure when all its couplings vanish.

	The check of the vector sector, and of the scalar sector, compared gvss (shape
	{nv,ns,ns}) with an empty array of shape {nv,nf,nf}. For ns!=nf the comparison
	always failed, so the decoupled photon of ScalYukawa.m and WessZumino.m added
	-Pi^2 T^4/45 to the free energy, and the placeholder scalar of a model without
	scalars added -Pi^2 T^4/90, while the photon of RealScalar.m (ns=nf=1 after
	padding) was dropped.

	Expected: the LO pressure is the ideal gas of the fields that take part in the
	model, Pi^2 T^4/90 (Nb + 7/8 Nf), with Nb bosonic and Nf fermionic degrees of
	freedom (two per vector, one per real scalar, two per Weyl fermion).
*)


idealGas[nb_, nf_] := Pi^2 T^4/90 (nb + 7/8 nf);


(*The potential of RealScalar.m and ScalYukawa.m for the first scalar*)
realScalarPotential[] := Module[{MassTerm},
	\[Lambda]1=GradTadpole[\[Sigma] CreateInvariant[Group,RepScalar,{{1},{True}}][[1]]];
	MassTerm=CreateInvariant[Group,RepScalar,{{1,1},{True,True}}][[1]];
	\[Mu]ij=GradMass[msq/2 MassTerm]//SparseArray;
	\[Lambda]3=GradCubic[\[Gamma]/6 CreateInvariant[Group,RepScalar,{{1,1,1},{True,True,True}}][[1]]];
	\[Lambda]4=GradQuartic[\[Lambda]/24 MassTerm^2];
];


(*
	The LO pressure of the model stored in the global tensors. Mode->0 computes only
	what it needs: PrintPressure then evaluates SymmetricPhaseLO directly, whereas
	the default mode stores the same value in SymmEnergy (checked for ScalYukawa).
*)
pressureLO[mode_:0] := (
	ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False,Mode->mode];
	PerformDRhard[];
	PrintPressure["LO"]
);


tests={};


(* ::Section:: *)
(*A U(1) that couples to nothing does not contribute (ns!=nf)*)


(*ScalYukawa.m: a real scalar and a Dirac fermion (ns=1, nf=2)*)
Group={"U1"};
CouplingName={g};
RepAdjoint={0};
RepScalar={{{0},"R"}};
RepFermion={{{0},"L"},{{0},"R"}};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJC,\[Mu]IJ,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
realScalarPotential[];
\[Mu]IJ=m\[Psi]*GradMassFermion[CreateInvariantFermion[Group,RepFermion,{{2,1},{False,True}}][[1]]];
\[Mu]IJC=m\[Psi]*GradMassFermion[CreateInvariantFermion[Group,RepFermion,{{1,2},{False,True}}][[1]]];
Ysff=y*GradYukawa[CreateInvariantYukawa[Group,RepScalar,RepFermion,{{1,2,1},{True,False,True}}][[1]]];
YsffC=y*GradYukawa[CreateInvariantYukawa[Group,RepScalar,RepFermion,{{1,1,2},{True,False,True}}][[1]]];
pScalYukawa=pressureLO[];

AppendTo[tests, VerificationTest[
	pressureLO[2],
	pScalYukawa,
	TestID->"same-LO-pressure-in-default-mode"
]];

AppendTo[tests, VerificationTest[
	pScalYukawa,
	idealGas[1, 4],
	SameTest->(Simplify[#1-#2]==0&),
	TestID->"decoupled-U1-not-counted-ScalYukawa"
]];
(*Reference value: one real scalar and one Dirac fermion, (1+7/8*4)/90*)
AppendTo[tests, VerificationTest[
	pScalYukawa,
	\[Pi]^2 T^4/20,
	TestID->"ScalYukawa-LO-pressure-value"
]];


(*WessZumino.m: two real scalars and one Weyl fermion (ns=2, nf=1)*)
Group={"U1"};
CouplingName={gY};
RepAdjoint={0};
RepScalar={{{0},"R"},{{0},"R"}};
RepFermion={{{0},"L"}};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
MassTerm1=CreateInvariant[Group,RepScalar,{{1,1},{True,True}}][[1]];
MassTerm2=CreateInvariant[Group,RepScalar,{{2,2},{True,True}}][[1]];
\[Mu]ij=GradMass[1/2 (mA MassTerm1+mB MassTerm2)]//SparseArray;
\[Lambda]4=GradQuartic[\[Lambda]/4 (MassTerm1+MassTerm2)^2];
\[Lambda]3=GradCubic[\[Lambda]c (CreateInvariant[Group,RepScalar,{{1,1,1},{True,True,True}}][[1]]
	+CreateInvariant[Group,RepScalar,{{1,2,2},{True,True,True}}][[1]])];
\[Lambda]1=GradTadpole[\[Lambda]Tad1 CreateInvariant[Group,RepScalar,{{1},{True}}][[1]]
	+\[Lambda]Tad2 CreateInvariant[Group,RepScalar,{{2},{True}}][[1]]];
\[Mu]IJ=mF/2*GradMassFermion[CreateInvariantFermion[Group,RepFermion,{{1,1},{True,False}}][[1]]];
\[Mu]IJC=\[Mu]IJ;
Ysff=yt/Sqrt[2]*GradYukawa[
	CreateInvariantYukawa[Group,RepScalar,RepFermion,{{1,1,1},{True,True,True}}][[1]]
	+I CreateInvariantYukawa[Group,RepScalar,RepFermion,{{2,1,1},{True,True,True}}][[1]]];
YsffC=Ysff/.Complex[0,a_]->Complex[0,-a];
pWessZumino=pressureLO[];

AppendTo[tests, VerificationTest[
	pWessZumino,
	idealGas[2, 2],
	SameTest->(Simplify[#1-#2]==0&),
	TestID->"decoupled-U1-not-counted-WessZumino"
]];


(* ::Section:: *)
(*The result does not depend on whether ns equals nf*)


(*RealScalar.m: a real scalar and no fermions (ns=nf=1 after padding)*)
Group={"U1"};
CouplingName={g};
RepAdjoint={0};
RepScalar={{{0},"R"}};
RepFermion={};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJC,\[Mu]IJ,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
realScalarPotential[];
pRealScalar=pressureLO[];

AppendTo[tests, VerificationTest[
	pRealScalar,
	idealGas[1, 0],
	SameTest->(Simplify[#1-#2]==0&),
	TestID->"decoupled-U1-not-counted-RealScalar"
]];
(*Adding a neutral Dirac fermion adds only its own ideal-gas term*)
AppendTo[tests, VerificationTest[
	pScalYukawa-pRealScalar,
	idealGas[0, 4],
	SameTest->(Simplify[#1-#2]==0&),
	TestID->"neutral-fermion-adds-only-its-own-term"
]];


(* ::Section:: *)
(*A model without scalars has no scalar contribution*)


(*QCD with one quark flavour: 8 gluons and 2x3 Weyl fermions, no scalars*)
Group={"SU3"};
CouplingName={gs};
RepAdjoint={{1,1}};
RepScalar={};
RepFermion={{{{1,0}},"L"},{{{1,0}},"R"}};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
pQCD=pressureLO[];

AppendTo[tests, VerificationTest[
	pQCD,
	idealGas[16, 12],
	SameTest->(Simplify[#1-#2]==0&),
	TestID->"no-placeholder-scalar-QCD"
]];


(* ::Section:: *)
(*Gauge fields that couple are counted*)


(*Abelian Higgs model: a charged scalar (2 real components) and a photon*)
Group={"U1"};
CouplingName={g1};
RepAdjoint={0};
RepScalar={{{Y\[Phi]},"C"}};
RepFermion={};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
MassTerm=CreateInvariant[Group,RepScalar,{{1,1},{True,False}}][[1]];
\[Mu]ij=GradMass[msq MassTerm]//SparseArray;
\[Lambda]4=GradQuartic[\[Lambda] MassTerm^2];
pAH=pressureLO[];

AppendTo[tests, VerificationTest[
	pAH,
	idealGas[4, 0],
	SameTest->(Simplify[#1-#2]==0&),
	TestID->"coupled-photon-counted-AbelianHiggs"
]];


(* ::Section:: *)
(*Report*)


report=TestReport[tests];
Print[report];
Exit[If[report["AllTestsSucceeded"], 0, 1]]
