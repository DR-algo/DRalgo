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
(*Temporal-scalar x heavy-scalar quartics when heavy scalars are integrated out*)


(*
	When heavy scalars are integrated out, PrepareSoftToSuperSoft (SoftToUS.m) collects
	all quartics of the heavy fields (temporal scalars A0 and heavy scalars H) in
	lamKTotal, including the A0 A0 H H quartics, which enter the heavy-field masses,
	the ultrasoft masses, and the NLO ultrasoft pressure.

	Model: scalar QED with a light scalar \[Phi] (charge 1) and a heavy scalar H (charge 2),
	H integrated out at the soft scale.
*)


(* ::Section:: *)
(*Model*)


Group={"U1"};
CouplingName={g1};
RepAdjoint={0};
RepScalar={{{1},"C"},{{2},"C"}}; (*\[Phi] (light), H (heavy)*)
RepFermion={};


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];


MassTerm1=CreateInvariant[Group,RepScalar,{{1,1},{True,False}}][[1]]//Simplify; (*\[Phi]^+\[Phi]*)
MassTerm2=CreateInvariant[Group,RepScalar,{{2,2},{True,False}}][[1]]//Simplify; (*H^+H*)
\[Mu]ij=GradMass[m1 MassTerm1+m2 MassTerm2]//Simplify//SparseArray;
\[Lambda]4=GradQuartic[l1 MassTerm1^2+l2 MassTerm2^2+l12 MassTerm1 MassTerm2];


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PerformDRhard[];
temporalCouplings=PrintTemporalScalarCouplings[];
PerformDRsoft[{3,4}];


tests={};


(* ::Section:: *)
(*The A0 A0 H H block of the heavy-field quartic equals the soft temporal-scalar couplings*)


(*
	lamKTotal has the indices {A0 (nv components), heavy scalars}; its A0 A0 H H block
	must be the A0 A0 S S coupling lamKVec of the soft theory restricted to H.
*)
nv=Length[gvvv];
heavy={3,4};
nH=Length[heavy];
lambdaKTotal=DRalgo`DRalgo`Private`lamKTotal;
lambdaKVec=DRalgo`DRalgo`Private`lamKVec;


AppendTo[tests,VerificationTest[
	Normal[lambdaKTotal[[1;;nv,1;;nv,nv+1;;nv+nH,nv+1;;nv+nH]]],
	Normal[lambdaKVec[[1;;nv,1;;nv,heavy,heavy]]],
	TestID->"A0A0HH-block-equals-temporal-scalar-couplings"
]];


AppendTo[tests,VerificationTest[
	FreeQ[Normal[lambdaKVec[[1;;nv,1;;nv,heavy,heavy]]],DRalgo`DRalgo`Private`lamVL],
	False,
	TestID->"A0A0HH-coupling-is-nonzero-in-this-model"
]];


(* ::Section:: *)
(*NLO ultrasoft pressure: figure-eight diagram with one A0 and one H loop*)


(*
	The term e3^2 q^2 A0^2 H^+H of |D0 H|^2 in the 3d Lagrangian makes the figure-eight
	diagram of the A0 and H tadpoles contribute the free energy e3^2 q^2 mD mH/(16\[Pi]^2),
	i.e. the pressure -e3^2 q^2 mD mH/(16\[Pi]^2). Here q=2 and e3^2=g1^2 T at leading order,
	so the mixed derivative of the pressure with respect to mD and mH is -g1^2 T/(4\[Pi]^2) at
	leading order. No other NLO term depends on both mD and mH.
*)
pressureNLO=PrintPressureUS["NLO"]/.{musqU1->mD^2,m23d->mH^2};
mixed=Simplify[D[pressureNLO,mD,mH],mD>0&&mH>0];
(*Leading order in the power counting g1^2 ~ l1 ~ l2 ~ l12*)
leadingOrder[expr_]:=Normal@Series[
	expr/.{g1->\[Epsilon] g1,l1->\[Epsilon]^2 l1,l2->\[Epsilon]^2 l2,l12->\[Epsilon]^2 l12},
	{\[Epsilon],0,2}]/.\[Epsilon]->1;


AppendTo[tests,VerificationTest[
	FreeQ[mixed,mD|mH],
	True,
	TestID->"pressure-cross-term-is-bilinear-in-mD-mH"
]];


AppendTo[tests,VerificationTest[
	leadingOrder[mixed/.temporalCouplings]//Simplify,
	-g1^2 T/(4\[Pi]^2),
	TestID->"pressure-cross-term-A0-H-figure-eight"
]];


(* ::Section:: *)
(*Report*)


report=TestReport[tests];
Print[report];
Print[Column[{#["TestID"],#["Outcome"]}&/@Values[report["TestResults"]]]];
Exit[If[report["AllTestsSucceeded"],0,1]]
