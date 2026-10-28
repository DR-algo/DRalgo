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
(*Light-mass insertion in the two-loop ultrasoft mass*)


(*
	When a heavy scalar S with a cubic coupling S \[Phi]^+\[Phi] is integrated out, the bubble with one light and one heavy line
	contributes to the two-loop light mass with one insertion of the light mass on the light line, from the hard region of the
	light line. The light mass inserted there is the tree-level mass and the momentum-independent one-loop terms: the A0 tadpole
	and the exchange of S between the light line and the tadpole of S. Together with the field normalization, all of them come
	with the same factor as the tree-level mass, -mum3d^2/(12\[Pi] ms3d^(3/2)), of which the field normalization is one quarter.

	Model: scalar QED with a light complex scalar \[Phi] and a heavy real singlet S, with the cubic mum/2 S \[Phi]^+\[Phi] and the
	tadpole t1 S; S is integrated out at the soft scale.
*)


(* ::Section:: *)
(*Model*)


Group={"U1"};
CouplingName={g1};
RepAdjoint={0};
RepScalar={{{1},"C"},{{0},"R"}}; (*\[Phi] (light), S (heavy)*)
RepFermion={};

{gvvv,gvff,gvss,lam1,lam3,lam4,muij,muIJ,muIJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];

P1=CreateInvariant[Group,RepScalar,{{1,1},{True,False}}][[1]]//Simplify; (*\[Phi]^+\[Phi]*)
S1=CreateInvariant[Group,RepScalar,{{2},{True}}][[1]];

muij=GradMass[m1 P1+ms/2 S1^2]//SparseArray;
lam4=GradQuartic[lam P1^2+ls/2 S1^2 P1+lss/4 S1^4];
lam3=GradCubic[mum/2 S1 P1];
lam1=GradTadpole[t1 S1];

ImportModelDRalgo[Group,gvvv,gvff,gvss,lam1,lam3,lam4,muij,muIJ,muIJC,Ysff,YsffC,Verbose->False];
PerformDRhard[];
PerformDRsoft[{3}];

lo=Expand[m13dUS/.PrintScalarMassUS["LO"]];
nlo=Expand[m13dUS/.PrintScalarMassUS["NLO"]];
temporal=Cases[lo/.t13d->0,lamVL[_],Infinity]//Union; (*the coupling of the A0 to the light scalar*)

(*The tadpole of S also enters the LO mass at one loop, with the couplings of S; only the terms without them are compared*)
noS={ls3d->0,lss3d->0,lamVL[_]->0};

(*The factor of the tree-level mass in the two-loop mass*)
factor=Coefficient[nlo,m13d mum3d^2];


tests={};


(* ::Section:: *)
(*Two-loop mass*)


AppendTo[tests,
VerificationTest[
	{Simplify[factor+1/(12\[Pi] ms3d^(3/2)),ms3d>0],Length[temporal]},
	{0,1},
	TestID->"tree-level-mass-insertion"
]];


(*The A0 tadpole on the light line*)
AppendTo[tests,
VerificationTest[
	Simplify[Coefficient[nlo,mum3d^2 Sqrt[musqU1] temporal[[1]]]-factor Coefficient[lo,Sqrt[musqU1] temporal[[1]]],ms3d>0],
	0,
	TestID->"A0-tadpole-insertion"
]];


(*The exchange of S between the light line and the tree-level tadpole of S*)
AppendTo[tests,
VerificationTest[
	{Simplify[(Coefficient[nlo,mum3d^3 t13d]-factor Coefficient[lo,mum3d t13d])/.noS,ms3d>0],(Coefficient[lo,mum3d t13d]/.noS)===0},
	{0,False},
	TestID->"tadpole-exchange-insertion"
]];


(* ::Section:: *)
(*Report*)


report=TestReport[tests];
Print[report];
Print[Column[{#["TestID"],#["Outcome"]}&/@Values[report["TestResults"]]]];
Exit[If[report["AllTestsSucceeded"],0,1]]
