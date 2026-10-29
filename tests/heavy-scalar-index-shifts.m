(* ::Package:: *)

If[$InputFileName=="",
	SetDirectory[NotebookDirectory[]],
	SetDirectory[DirectoryName[$InputFileName]]
];

DRalgo`DRalgo`$LoadGroupMath=True;

Check[
    Get["../Kernel/DRalgo.m"],
    Message[Get::noopen, "DRalgo` at "<>ToString[$UserBaseDirectory]<>"/Applications"];
    Abort[];
]


(* ::Chapter:: *)
(*Index shifts of the heavy-scalar tensors*)


(*
	PrepareSoftToSuperSoft (SoftToUS.m) shifts the indices of the heavy scalars by the number of temporal
	scalars, so that the heavy scalars follow the temporal scalars, and builds the tensors of the light and heavy
	fields with fixed dimensions. A tensor with few nonzero components, e.g. the light-light-heavy cubic with as many
	components as indices, must be shifted as well, and must keep its dimensions.

	Model: scalar QED with a light complex scalar \[Phi] and two heavy real singlets S and C, with the cubics
	mum/2 S \[Phi]^+\[Phi], c3/3 S^3, and kap/2 S C^2; S and C are integrated out at the soft scale.
*)


(* ::Section:: *)
(*Model*)


Group={"U1"};
CouplingName={g1};
RepAdjoint={0};
RepScalar={{{1},"C"},{{0},"R"},{{0},"R"}}; (*\[Phi] (light), C and S (heavy)*)
RepFermion={};

{gvvv,gvff,gvss,lam1,lam3,lam4,muij,muIJ,muIJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];

P1=CreateInvariant[Group,RepScalar,{{1,1},{True,False}}][[1]]//Simplify; (*\[Phi]^+\[Phi]*)
C1=CreateInvariant[Group,RepScalar,{{2},{True}}][[1]];
S1=CreateInvariant[Group,RepScalar,{{3},{True}}][[1]];

muij=GradMass[m1 P1+mc/2 C1^2+ms/2 S1^2]//SparseArray;
lam4=GradQuartic[lam P1^2+lc/2 C1^2 P1+ls/2 S1^2 P1+lss/4 S1^4+lcc/4 C1^4+lsc/4 S1^2 C1^2];
lam3=GradCubic[mum/2 S1 P1+c3/3 S1^3+kap/2 S1 C1^2];

ImportModelDRalgo[Group,gvvv,gvff,gvss,lam1,lam3,lam4,muij,muIJ,muIJC,Ysff,YsffC,Verbose->False];
PerformDRhard[];


tests={};


(* ::Section:: *)
(*Ultrasoft mass*)


messages=Reap[Check[PerformDRsoft[{3,4}],Sow["failed"]]][[2]];

AppendTo[tests,
VerificationTest[
	messages,
	{},
	TestID->"PerformDRsoft-without-messages"
]];


(*
	At LO, the light mass gets the tadpoles of S and C, the A0 tadpole, the S exchange between the light line and
	the one-loop tadpole of S, and the bubble of the light scalar with S, with <S^2>=-Sqrt[ms3d]/(4\[Pi]) and <C^2>=-Sqrt[mc3d]/(4\[Pi]).
	The tadpole of S has the loop of the A0, with the coupling kapA/2 A0^2 S of the soft theory, which arises at one loop.
*)
kapA=Normal[DRalgo`DRalgo`Private`lamvvsLS][[1,1,4]]/.DRalgo`DRalgo`Private`lamVVSL->lamVVSL;
expected=(m13d
	-ls3d Sqrt[ms3d]/(8\[Pi])-lc3d Sqrt[mc3d]/(8\[Pi])
	-mum3d/(2ms3d) (-c33d Sqrt[ms3d]/(4\[Pi])-kap3d Sqrt[mc3d]/(8\[Pi])-kapA Sqrt[musqU1]/(8\[Pi]))
	-mum3d^2/(16\[Pi] Sqrt[ms3d])
	-lamVL[1] Sqrt[musqU1]/(8\[Pi]));

AppendTo[tests,
VerificationTest[
	Simplify[(m13dUS/.PrintScalarMassUS["LO"])-expected,{ms3d>0,mc3d>0,musqU1>0}],
	0,
	TestID->"LO-ultrasoft-mass-with-two-heavy-singlets"
]];


AppendTo[tests,
VerificationTest[
	FreeQ[PrintScalarMassUS["NLO"],_Part|_Dot|_Thread|_SparseArray],
	True,
	TestID->"NLO-ultrasoft-mass-is-evaluated"
]];


(* ::Section:: *)
(*Report*)


report=TestReport[tests];
Print[report];
Print[Column[{#["TestID"],#["Outcome"]}&/@Values[report["TestResults"]]]];
Exit[If[report["AllTestsSucceeded"],0,1]]
