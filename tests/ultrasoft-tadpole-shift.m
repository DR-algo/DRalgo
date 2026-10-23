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
(*Tree-level tadpoles of heavy scalars in the ultrasoft mass*)


(*
	A tree-level tadpole t1 of a heavy scalar S shifts S by v0=-t1/M^2. The LO ultrasoft mass of the light scalar contains,
	besides the tree-level exchange of S, the one-loop terms at first order in v0: the one-loop light-light-heavy vertex times v0,
	the shift of the one-loop tadpole of S by v0, and the H^2 S^2 quartic with v0 and the one-loop shift of S.

	The shift S -> S+s of the 3d theory changes the parameters and creates the tadpole t1=ms s+c3 s^2+lss s^3 from t1=0.
	The ultrasoft theory does not depend on s, so that the LO ultrasoft mass has no term linear in s. The couplings of the
	temporal scalars are set to zero, since the one-loop couplings of two temporal scalars to a heavy scalar are not shifted.

	Model: scalar QED with a light complex scalar \[Phi] and a heavy real singlet S with the cubics mum/2 S \[Phi]^+\[Phi]
	and c3/3 S^3, the quartics ls/2 S^2 \[Phi]^+\[Phi] and lss/4 S^4, and the tadpole t1 S; S is integrated out at the soft scale.
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
lam3=GradCubic[mum/2 S1 P1+c3/3 S1^3];
lam1=GradTadpole[t1 S1];

ImportModelDRalgo[Group,gvvv,gvff,gvss,lam1,lam3,lam4,muij,muIJ,muIJC,Ysff,YsffC,Verbose->False];
PerformDRhard[];
PerformDRsoft[{3}];

noTemporal={lamVL[_]->0,lamVLL[_]->0,lamVVSL[_]->0};
lo=(m13dUS/.PrintScalarMassUS["LO"])/.noTemporal;

(*The 3d parameters after S -> S+s, starting from t1=0*)
shift={
	m13d->m13d+mum3d s/2+ls3d s^2/2,
	ms3d->ms3d+2c33d s+3lss3d s^2,
	t13d->ms3d s+c33d s^2+lss3d s^3,
	c33d->c33d+3lss3d s,
	mum3d->mum3d+2ls3d s};


tests={};


(* ::Section:: *)
(*Shift of the heavy scalar*)


AppendTo[tests,
VerificationTest[
	{FreeQ[lo,t13d],Simplify[(lo/.t13d->0)+(D[lo,t13d]/.t13d->0) t13d-lo]},
	{False,0},
	TestID->"LO-mass-is-linear-in-the-tadpole"
]];


AppendTo[tests,
VerificationTest[
	Simplify[D[lo/.shift,s]/.s->0,ms3d>0],
	0,
	TestID->"LO-mass-is-invariant-under-a-shift-of-the-heavy-scalar"
]];


(* ::Section:: *)
(*Report*)


report=TestReport[tests];
Print[report];
Print[Column[{#["TestID"],#["Outcome"]}&/@Values[report["TestResults"]]]];
Exit[If[report["AllTestsSucceeded"],0,1]]
