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
(*Tree-level tadpole and mixing of a heavy scalar with a light scalar*)


(*
	A light real singlet s and a heavy real singlet S with all cubic and quartic couplings, the mixing mx s S, and a tadpole
	t1 S; S is integrated out at the soft scale. The shift v0=-t1/M^2 of S gives the mixing c5 v0 s S through the cubic
	c5/2 s S^2, which connects the light leg to the one-loop mixing and to the one-loop tadpole of S. The mixing and the
	tadpole are counted as small, and the results are of first order in both.

	The LO ultrasoft mass of s is compared with the tree-level and the hard one-loop effective potential of s and S,
	with S eliminated at tree level, computed here from 1/2 Tr ln of the mass matrix with the light line expanded in the
	mixing, at zeroth and first order in the tadpole and the mixing. The one-loop 1PR diagrams of the mass and the quartic
	(OneLoopOnePRSS) are compared with the same potential.

	The shift S -> S+x of the 3d theory changes the parameters and creates the tadpole t1=ms x+c4 x^2/2+lss x^3 and the
	mixing mx=c5 x+l5 x^2/2 from t1=mx=0. The ultrasoft theory does not depend on x, so that the LO ultrasoft mass and
	quartic have no terms linear in x.
*)


(* ::Section:: *)
(*Model*)


Group={"U1"};
CouplingName={g1};
RepAdjoint={0};
RepScalar={{{0},"R"},{{0},"R"}}; (*s (light), S (heavy)*)
RepFermion={};

{gvvv,gvff,gvss,lam1,lam3,lam4,muij,muIJ,muIJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];

s1=CreateInvariant[Group,RepScalar,{{1},{True}}][[1]];
S1=CreateInvariant[Group,RepScalar,{{2},{True}}][[1]];

muij=GradMass[m1/2 s1^2+ms/2 S1^2+mx s1 S1]//SparseArray;
lam4=GradQuartic[lam/4 s1^4+lss/4 S1^4+lm/4 s1^2 S1^2+l5/6 s1 S1^3+l4/6 s1^3 S1];
lam3=GradCubic[c1/6 s1^3+c4/6 S1^3+c5/2 s1 S1^2+c6/2 s1^2 S1];
lam1=GradTadpole[t1 S1];

ImportModelDRalgo[Group,gvvv,gvff,gvss,lam1,lam3,lam4,muij,muIJ,muIJC,Ysff,YsffC,Verbose->False];
PerformDRhard[];
PerformDRsoft[{2}];

lo=m13dUS/.PrintScalarMassUS["LO"];
quartic=lam3dUS/.PrintCouplingsUS[];

(*The terms of order n in the tadpole and the mixing*)
order[expr_,n_]:=Coefficient[expr/.{t13d->tau t13d,mx3d->tau mx3d},tau,n];

(*The 3d parameters after S -> S+x, starting from t1=mx=0*)
shift={
	m13d->m13d+c63d x+lm3d x^2/2,
	ms3d->ms3d+c43d x+3lss3d x^2,
	mx3d->c53d x+l53d x^2/2,
	t13d->ms3d x+c43d x^2/2+lss3d x^3,
	c13d->c13d+l43d x,
	c43d->c43d+6lss3d x,
	c53d->c53d+l53d x,
	c63d->c63d+lm3d x};


(* ::Section:: *)
(*Effective potential*)


(*
	The tree-level mass, and the one-loop mass and quartic of the tensors, d^2/ds^2 and d^4/ds^4, as in
	ultrasoft-tadpole-shift.m; t counts the order of the tadpole and the mixing
*)
{treeMass,loopMass,loopQuartic}=Module[{int,V0,ans,eqs,cs,star,Mss,Msf,Mff,V1,L,trunc,tree,tot,t,f,sv,c},
	int[a_,b_,A_]:=A^(3/2-a-b) Gamma[3/2-a]Gamma[a+b-3/2]/((4Pi)^(3/2)Gamma[3/2]Gamma[b]);
	V0[sv_,f_]:=m13d/2 f^2+lam3d/4 f^4+ms3d/2 sv^2+lss3d/4 sv^4+lm3d/4 f^2 sv^2+l53d/6 f sv^3+l43d/6 f^3 sv
		+c13d/6 f^3+c43d/6 sv^3+c53d/2 f sv^2+c63d/2 f^2 sv+t t13d sv+t mx3d f sv;
	trunc[x_]:=Normal[Series[x,{f,0,4},{t,0,1}]];
	ans=Sum[c[i,j] f^i t^j,{i,0,4},{j,0,1}]/.c[0,0]->0;
	eqs=Flatten[CoefficientList[Normal[Series[D[V0[sv,f],sv]/.sv->ans,{f,0,4},{t,0,1}]],{f,t}]];
	cs=Solve[Thread[eqs==0],Cases[ans,_c,Infinity]][[1]];
	star=ans/.cs;
	tree=Simplify[trunc[V0[star,f]]/.t->1];
	Mss=D[V0[sv,f],{sv,2}]/.sv->star; Msf=D[V0[sv,f],sv,f]/.sv->star; Mff=D[V0[sv,f],{f,2}]/.sv->star;
	L=Mff-(Mff/.f->0);
	V1=1/2(-Mss^(3/2)/(6Pi)+Sum[(-1)^(n+1)/n Sum[Binomial[n,k]L^(n-k)(-Msf^2)^k int[n,k,Mss],{k,1,n}],{n,1,5}]);
	tot=Simplify[trunc[V1]/.t->1,ms3d>0];
	{2 Coefficient[tree,f,2],2 Coefficient[tot,f,2],24 Coefficient[tot,f,4]}
];

{onePRMass,onePRQuartic}=DRalgo`DRalgo`Private`OneLoopOnePRSS[];


tests={};


(* ::Section:: *)
(*LO mass and 1PR diagrams*)


AppendTo[tests,
VerificationTest[
	Table[Simplify[order[lo-treeMass-loopMass,n],ms3d>0],{n,0,1}],
	{0,0},
	TestID->"LO-mass-equals-the-effective-potential"
]];


AppendTo[tests,
VerificationTest[
	Table[Simplify[{
		Normal[onePRMass[n]][[1,1]]-order[loopMass,n],
		Normal[onePRQuartic[n]][[1,1,1,1]]-order[loopQuartic,n]},ms3d>0],{n,0,1}],
	{{0,0},{0,0}},
	TestID->"1PR-diagrams-equal-the-effective-potential"
]];


(* ::Section:: *)
(*Shift of the heavy scalar*)


AppendTo[tests,
VerificationTest[
	{Simplify[D[lo/.shift,x]/.x->0,ms3d>0],Simplify[D[quartic/.shift,x]/.x->0,ms3d>0]},
	{0,0},
	TestID->"LO-mass-and-quartic-are-invariant-under-a-shift-of-the-heavy-scalar"
]];


(* ::Section:: *)
(*Report*)


report=TestReport[tests];
Print[report];
Print[Column[{#["TestID"],#["Outcome"]}&/@Values[report["TestResults"]]]];
Exit[If[report["AllTestsSucceeded"],0,1]]
