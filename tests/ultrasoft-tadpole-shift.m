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
(*Tree-level tadpoles of heavy scalars in the ultrasoft mass and quartic*)


(*
	A tree-level tadpole t1 of a heavy scalar S shifts S by v0=-t1/M^2. The LO ultrasoft mass and quartic of the light
	scalar contain, besides the tree-level diagrams, the one-loop one-particle-reducible diagrams at first order in v0
	(OneLoopOnePRSS), and the quartic contains the field normalization at first order in v0 (TadpoleZSS).

	The shift S -> S+s of the 3d theory changes the parameters and creates the tadpole t1=ms s+c3 s^2+lss s^3 from t1=0.
	The ultrasoft theory does not depend on s, so that the LO ultrasoft mass and quartic have no terms linear in s. The
	couplings of the temporal scalars are set to zero, since the one-loop couplings of two temporal scalars to a heavy scalar
	are not shifted.

	The 1PR diagrams are also compared with the one-loop effective potential of the hard modes of \[Phi]1 and S, with S
	eliminated at tree level, computed here from 1/2 Tr ln of the mass matrix with the light line expanded in the mixing.

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
quartic=(lam3dUS/.PrintCouplingsUS[])/.noTemporal;

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


AppendTo[tests,
VerificationTest[
	Simplify[D[quartic/.shift,s]/.s->0,ms3d>0],
	0,
	TestID->"quartic-is-invariant-under-a-shift-of-the-heavy-scalar"
]];


(* ::Section:: *)
(*1PR diagrams and the effective potential*)


(*
	The hard one-loop effective potential of the real component f of \[Phi] and of S: 1/2 Tr ln(p^2+M^2) of the mass matrix M^2,
	with the light line expanded in the light-light entry L and in the mixing X/(p^2+A), at the tree-level solution S(f).
	As in OneLoopOnePRSS, L has no field-independent part, which is an insertion of the light mass.
*)
effectivePotential=Module[{int,V0,ans,eqs,cs,star,Mss,Msf,Mff,V1,L,trunc,tot,d,t,f,sv,c},
	int[a_,b_,A_]:=A^(3/2-a-b) Gamma[3/2-a]Gamma[a+b-3/2]/((4Pi)^(3/2)Gamma[3/2]Gamma[b]);
	V0[sv_,f_]:=m13d/2 f^2+lam3d/4 f^4+ms3d/2 sv^2+t t13d sv+c33d/3 sv^3+lss3d/4 sv^4+mum3d/4 sv f^2+ls3d/4 sv^2 f^2;
	trunc[x_]:=Normal[Series[x,{f,0,4},{t,0,1}]];
	ans=Sum[c[i,j] f^i t^j,{i,0,4,2},{j,0,1}]/.c[0,0]->0;
	eqs=Flatten[CoefficientList[Normal[Series[D[V0[sv,f],sv]/.sv->ans,{f,0,4},{t,0,1}]],{f,t}]];
	cs=Solve[Thread[eqs==0],Cases[ans,_c,Infinity]][[1]];
	star=ans/.cs;
	Mss=D[V0[sv,f],{sv,2}]/.sv->star; Msf=D[V0[sv,f],sv,f]/.sv->star; Mff=D[V0[sv,f],{f,2}]/.sv->star;
	L=Mff-(Mff/.f->0);
	V1=1/2(-Mss^(3/2)/(6Pi)+Sum[(-1)^(n+1)/n Sum[Binomial[n,k]L^(n-k)(-Msf^2)^k int[n,k,Mss],{k,1,n}],{n,1,5}]);
	tot=Simplify[trunc[V1]/.t->1,ms3d>0];
	(*The one-loop mass and quartic of the tensors, d^2/df^2 and d^4/df^4*)
	{2 Coefficient[tot,f,2],24 Coefficient[tot,f,4]}
];

{onePRMass,onePRQuartic}=DRalgo`DRalgo`Private`OneLoopOnePRSS[];
(*The tensors of DRalgo contain the internal names of the couplings of the temporal scalars*)
noTemporalPrivate={DRalgo`DRalgo`Private`lamVL[_]->0,DRalgo`DRalgo`Private`lamVLL[_]->0,DRalgo`DRalgo`Private`lamVVSL[_]->0};

AppendTo[tests,
VerificationTest[
	Table[Simplify[{
		(Normal[onePRMass[n]][[1,1]]/.noTemporalPrivate)-Coefficient[effectivePotential[[1]],t13d,n]t13d^n,
		(Normal[onePRQuartic[n]][[1,1,1,1]]/.noTemporalPrivate)-Coefficient[effectivePotential[[2]],t13d,n]t13d^n},ms3d>0],{n,0,1}],
	{{0,0},{0,0}},
	TestID->"1PR-diagrams-equal-the-effective-potential"
]];


(* ::Section:: *)
(*Report*)


report=TestReport[tests];
Print[report];
Print[Column[{#["TestID"],#["Outcome"]}&/@Values[report["TestResults"]]]];
Exit[If[report["AllTestsSucceeded"],0,1]]
