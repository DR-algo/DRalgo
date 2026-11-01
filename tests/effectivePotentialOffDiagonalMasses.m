(* ::Package:: *)

If[$InputFileName=="",
	SetDirectory[NotebookDirectory[]],
	SetDirectory[DirectoryName[$InputFileName]]
];
(*Two models are created in this file*)
DRalgo`DRalgo`$GroupMathMultipleModels=True;

DRalgo`DRalgo`$LoadGroupMath=True;
(*DRalgo`DRalgo`$InstallGroupMath=True;*)

Check[
    Get["../Kernel/DRalgo.m"],
    Message[Get::noopen, "DRalgo` at "<>ToString[$UserBaseDirectory]<>"/Applications"];
    Abort[];
]


(*
	Off-diagonal masses as insertions in the effective potential:
	CalculatePotentialUS[PerturbativeDiagonalization->True] and, for masses
	that involve heavy fields, CalculatePotentialHET[].

	The off-diagonal part P of the field-dependent mass matrix is inserted on
	the one-loop diagram. To second order in P this gives
		-1/4 Sum[BubbleIntegral3d[m_i,m_l] P_il P_li], with the bubble 1/(4\[Pi](Sqrt[x]+Sqrt[y])),
	per real scalar, the O(P^2) term of 1/2 Tr ln(p^2+m^2+P), and twice that per
	vector (Landau gauge). It must agree with the exact one-loop potential
	-Sum[\[Lambda]^(3/2)]/(12\[Pi]) over the eigenvalues \[Lambda] of the full
	mass matrix, expanded to second order in P.
*)


(* ::Section:: *)
(*Helpers*)


(*The one-loop potential of real scalars with these squared masses*)
oneLoop[eigenvalues_List]:=Total[-eigenvalues^(3/2)/(12 \[Pi])];

(*The O(P^2) term of the exact one-loop potential of real scalars with the mass matrix D+P*)
exactSecondOrder[mass_]:=Module[{diag,eps},
	diag=DiagonalMatrix[Diagonal[mass]];
	SeriesCoefficient[oneLoop[Eigenvalues[diag+eps (mass-diag)]],{eps,0,2}]
];

tests={};


(* ::Section:: *)
(*Model: two real scalars with a mass mixing and no interactions*)


Group={"U1"};
RepAdjoint={0};
CouplingName={g1};
RepScalar={{{0},"R"},{{0},"R"}};
RepFermion={};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
\[Mu]ij=SparseArray[{{m1,m12},{m12,m2}}];

ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False,Mode->0];
DefineNewTensorsUS[\[Mu]ij,\[Lambda]4,\[Lambda]3,gvss,gvvv];
DefineVEVS[{0,0}];
CalculatePotentialUS[PerturbativeDiagonalization->True];
(*Without interactions the two-loop order is the insertion alone*)
insertionScalars=PrintEffectivePotential["NNLO"];


(* ::Section:: *)
(*Two scalars: the insertion is the O(P^2) term of the exact diagonalization*)


(*m1 = m2 + c with c > 0, so that Sqrt[(m1-m2)^2] = c*)
AppendTo[tests,VerificationTest[
	Simplify[(insertionScalars/.m1->m2+c)-exactSecondOrder[{{m2+c,m12},{m12,m2}}],
		Assumptions->c>0&&m2>0&&m12>0],
	0,
	TestID->"two-scalars-match-exact-diagonalization"
]];

AppendTo[tests,VerificationTest[
	Simplify[insertionScalars+m12^2/(8 \[Pi] (Sqrt[m1]+Sqrt[m2]))],
	0,
	TestID->"two-scalars-reference-value"
]];


(* ::Section:: *)
(*HEFT: a heavy scalar mixing with a light one*)


(*CalculatePotentialUS has kept only the diagonal masses, so the background is set again*)
DefineVEVS[{0,0}];
PrepareHET[{2},{}];
CalculatePotentialHET[];
insertionHET=PrintActionHET["NLO"];

(*The light line is massless (the hard region): the heavy eigenvalue of {{0,P},{P,m2}}*)
heavyEigenvalue=Last[SortBy[Eigenvalues[{{0,eps m12},{eps m12,m2}}],(#/.{eps->0,m2->1})&]];

AppendTo[tests,VerificationTest[
	Simplify[insertionHET-SeriesCoefficient[oneLoop[{heavyEigenvalue}],{eps,0,2}],
		Assumptions->m2>0&&m12>0],
	0,
	TestID->"het-heavy-light-mixing-is-the-hard-region"
]];

AppendTo[tests,VerificationTest[
	Simplify[insertionHET+m12^2/(8 \[Pi] Sqrt[m2])],
	0,
	TestID->"het-reference-value"
]];


(* ::Section:: *)
(*Model: a complex scalar with charges 1 and 2 under U(1)xU(1), which mixes the vectors*)


Group={"U1","U1"};
RepAdjoint={0,0};
CouplingName={g1,g2};
RepScalar={{{1,2},"C"}};
RepFermion={};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
MassTerm=CreateInvariant[Group,RepScalar,{{1,1},{True,False}}][[1]];
\[Mu]ij=GradMass[msq MassTerm]//SparseArray;

ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False,Mode->0];
DefineNewTensorsUS[\[Mu]ij,\[Lambda]4,\[Lambda]3,gvss,gvvv];
DefineVEVS[{\[CurlyPhi],0}];
vectorMass=PrintTensorsVEV[2]//Normal;

CalculatePotentialUS[PerturbativeDiagonalization->True];
withInsertions=PrintEffectivePotential["NNLO"];
(*The masses are diagonal now: the same two-loop potential without the insertions*)
CalculatePotentialUS[];
insertionVectors=withInsertions-PrintEffectivePotential["NNLO"];


(* ::Section:: *)
(*Vectors: twice the insertion of a real scalar (two polarizations in Landau gauge)*)


AppendTo[tests,VerificationTest[
	vectorMass[[1,2]]=!=0,
	True,
	TestID->"vectors-mix"
]];

(*At rational values with rational vector masses; the comparison is exact*)
values={g1->2/3,g2->1/5,\[CurlyPhi]->7/4};
AppendTo[tests,VerificationTest[
	Simplify[(insertionVectors/.values)-2 exactSecondOrder[vectorMass/.values]],
	0,
	TestID->"vectors-match-exact-diagonalization"
]];

AppendTo[tests,VerificationTest[
	Simplify[insertionVectors+g1^2 g2^2 \[CurlyPhi]^3/(\[Pi] (g1+2 g2)),
		Assumptions->g1>0&&g2>0&&\[CurlyPhi]>0],
	0,
	TestID->"vectors-reference-value"
]];


report=TestReport[tests];
Print[report];
Print[Column[{#["TestID"],#["Outcome"]}&/@Values[report["TestResults"]]]];
Exit[If[report["AllTestsSucceeded"],0,1]]
