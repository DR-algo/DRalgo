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
(*Two-loop ultrasoft masses: double bubbles and insertions on heavy lines*)


(*
	The two-loop mass of a light real singlet f with heavy real singlets, from the second derivative of the two-loop
	potential with the masses shifted by the cubics: the figure-8 diagram (1/8)lam4[a,b,c,d]G[a,b]G[c,d] gives
	- the double bubble 1/4 lam3[f,a,b]B[a,b]lam4[a,b,c,d]B[c,d]lam3[c,d,f], with the bubble B[x,y]=1/(4Pi(Sqrt[x]+Sqrt[y]));
	- the tadpole 1/2 lam4[a,a,c,d]A[a] on a line of the one-loop bubbles lam3 lam3 B, with the tadpole A[x]=-Sqrt[x]/(4Pi)
	  and the triangle T[x,y,z] at zero momenta.
	Likewise, the shift of a heavy field by its one-loop tadpole changes its mass through its cubic, on a line of the bubbles.

	Checked here, with the 3d masses M of the heavy singlets:
	- cubics c1/2 f s1^2 and c2/2 f s2^2, quartic l12/4 s1^2 s2^2: the double bubble gives c1 c2 l12 B[M1,M1]B[M2,M2]/2;
	- the same with l11/24 s1^4: the double bubble and the tadpole on a line cancel, since the second derivative of
	  A[M+c1 f]^2 vanishes for A ~ Sqrt[M];
	- cubic ca/2 f^2 s and quartic dd/24 s^4: the tadpole dd A[M]/2 on the heavy line of the bubble with f and s gives
	  ca^2 dd A[M]T[0,M,M]/2=-ca^2 dd/(64Pi^2 M);
	- cubics ca/2 f^2 s and cs/6 s^3: the one-loop tadpole cs A[M]/2 shifts the mass of s by -cs^2 A[M]/(2M), which gives
	  ca^2 cs^2/(64Pi^2 M^2) on the heavy line of the bubble, and the bubble cs^2 of two s on the line of s between the cubics of
	  the light legs gives -ca^2 cs^2 I3/2 with the sunset I3=(6 Log[3/2]-1)/(96Pi^2 M^2) of the masses {0,m,m} on one line
	  and m, m, by partial fractions, so that ca^2 cs^2 has (2-3 Log[3/2])/(96Pi^2 M^2);
	- with these couplings, the dependence of the two-loop mass on mu3 cancels the running of the ultrasoft and soft masses.
	The sunset with two cubics on one line, i.e. a bubble on a line of the bubble of the light legs, has three propagators on
	this line.
*)


testList={};

Group={"U1"};RepAdjoint={0};CouplingName={g1};RepFermion={};


(* ::Section:: *)
(*Two heavy singlets*)


RepScalar={{{0},"R"},{{0},"R"},{{0},"R"}};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
{f1,s1,s2}=Table[CreateInvariant[Group,RepScalar,{{i},{True}}][[1]],{i,3}];
\[Mu]ij=GradMass[mf/2 f1^2+m1/2 s1^2+m2/2 s2^2]//SparseArray;
\[Lambda]3=GradCubic[c1/2 f1 s1^2+c2/2 f1 s2^2];
\[Lambda]4=GradQuartic[l12/4 s1^2 s2^2+l11/24 s1^4];
ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
(*Without the quartics of f, the identification of the quartics emits messages, which do not affect the masses*)
Quiet[PerformDRhard[]];
PerformDRsoft[{2,3}];
nlo=Expand[mf3dUS/.PrintScalarMassUS["NLO"]];

AppendTo[testList,
VerificationTest[
	Simplify[Coefficient[nlo,c13d c23d l123d]-1/(128\[Pi]^2 Sqrt[m13d]Sqrt[m23d])],
	0,
	TestID->"double-bubble"
]];

AppendTo[testList,
VerificationTest[
	Simplify[Coefficient[nlo,c13d^2 l113d]],
	0,
	TestID->"double-bubble-and-tadpole-on-a-heavy-line-cancel"
]];


(* ::Section:: *)
(*One heavy singlet*)


RepScalar={{{0},"R"},{{0},"R"}};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
{f1,s1}=Table[CreateInvariant[Group,RepScalar,{{i},{True}}][[1]],{i,2}];
\[Mu]ij=GradMass[mf/2 f1^2+ms/2 s1^2]//SparseArray;
\[Lambda]3=GradCubic[ca/2 f1^2 s1+cs/6 s1^3];
\[Lambda]4=GradQuartic[dd/24 s1^4+lf/24 f1^4];
ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
Quiet[PerformDRhard[]];
betaSoft=mf/.BetaFunctions3DS[];
PerformDRsoft[{2}];
nlo=Expand[mf3dUS/.PrintScalarMassUS["NLO"]];

AppendTo[testList,
VerificationTest[
	Simplify[Coefficient[nlo,ca3d^2 dd3d]],
	-1/(64\[Pi]^2 ms3d),
	TestID->"tadpole-on-a-heavy-line"
]];

AppendTo[testList,
VerificationTest[
	Simplify[Coefficient[nlo,ca3d^2 cs3d^2]-(2-3Log[3/2])/(96\[Pi]^2 ms3d^2)],
	0,
	TestID->"one-loop-tadpole-shift-and-bubble-on-a-heavy-line"
]];

(*
	The dependence on mu3 of the two-loop mass equals the running in the ultrasoft theory, with the ultrasoft couplings at tree
	level, minus that in the soft theory; the sunset with the bubble of two f on the line of s, between the cubics of the
	light legs, gives the term ca^4/(32Pi^2 M^2) of the running
*)
tree=Table[r[[1]]->(Normal[Series[r[[2]]/.Pi->Pi/eps,{eps,0,0}]]/.eps->0/.muH[_]->0),{r,PrintCouplingsUS[]}];
AppendTo[testList,
VerificationTest[
	Simplify[mu3 D[nlo,mu3]-((mf3dUS/.BetaFunctions3DUS[]/.tree)-betaSoft)],
	0,
	TestID->"mu3-dependence-of-the-two-loop-mass"
]];


(* ::Section:: *)
(*A charged light scalar*)


(*
	Scalar QED with a light complex scalar phi and a heavy real singlet S with the cubic mum/2 S phi^+phi and a tadpole t1 S: the
	terms with g1^2 of the dependence on mu3 of the two-loop mass, from the vector loop on the line of phi between the cubics, equal
	the running in the ultrasoft theory, with the quartic -mum^2/(8M) from the exchange of S, minus that of the soft mass and of
	the tadpole, which enters the tree-level mass as -mum t1/(2M) and runs with the loop of phi with the vector
*)
RepScalar={{{1},"C"},{{0},"R"}};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
p1=CreateInvariant[Group,RepScalar,{{1,1},{True,False}}][[1]]//Simplify;
s1=CreateInvariant[Group,RepScalar,{{2},{True}}][[1]];
\[Mu]ij=GradMass[m1 p1+ms/2 s1^2]//SparseArray;
\[Lambda]1=GradTadpole[t1 s1];
\[Lambda]3=GradCubic[mum/2 s1 p1];
\[Lambda]4=GradQuartic[lam p1^2+ls/2 s1^2 p1+lss/4 s1^4];
ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PerformDRhard[];
betaSoftQED=BetaFunctions3DS[];
PerformDRsoft[{3}];
nloQED=m13dUS/.PrintScalarMassUS["NLO"];
treeQED=Table[r[[1]]->(Normal[Series[r[[2]]/.Pi->Pi/eps,{eps,0,0}]]/.eps->0/.muH[_]->0),{r,PrintCouplingsUS[]}];
treeMass=m13d-mum3d t13d/(2ms3d);
mismatchQED=mu3 D[nloQED,mu3]-((m13dUS/.BetaFunctions3DUS[]/.treeQED/.g13dUS->g13d)-(m1/.betaSoftQED)-D[treeMass,t13d](t1/.betaSoftQED));

AppendTo[testList,
VerificationTest[
	{Simplify[treeMass-(Normal[Series[(m13dUS/.PrintScalarMassUS["LO"])/.Pi->Pi/eps,{eps,0,0}]]/.eps->0)],Simplify[Coefficient[Expand[mismatchQED/.t13d->0],g13d,2]]},
	{0,0},
	TestID->"mu3-dependence-with-a-vector-loop"
]];


(* ::Section:: *)
(*Report*)


report=TestReport[testList];
Print[report];
Print[Column[{#["TestID"],#["Outcome"]}&/@Values[report["TestResults"]]]];
Exit[If[report["AllTestsSucceeded"],0,1]]
