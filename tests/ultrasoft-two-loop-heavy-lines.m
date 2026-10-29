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
	  ca^2 cs^2/(64Pi^2 M^2) on the heavy line of the bubble.
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
	Simplify[Coefficient[nlo,ca3d^2 cs3d^2]],
	1/(64\[Pi]^2 ms3d^2),
	TestID->"one-loop-tadpole-shift-on-a-heavy-line"
]];


(* ::Section:: *)
(*Report*)


report=TestReport[testList];
Print[report];
Print[Column[{#["TestID"],#["Outcome"]}&/@Values[report["TestResults"]]]];
Exit[If[report["AllTestsSucceeded"],0,1]]
