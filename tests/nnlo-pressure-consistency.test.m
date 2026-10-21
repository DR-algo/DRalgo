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
(*NNLO pressure: RG invariance and comparison with the literature*)


(*
	PrintPressure["NNLO"] is the 3-loop pressure of the hard modes. Together with the
	pressure of the soft theory it reproduces the full pressure, and at fixed \[Mu]3 it does
	not depend on \[Mu] at O(g^4): its explicit \[Mu]-dependence cancels the running of the
	NLO pressure (and of the vacuum energy, which DRalgo does not track).

	Up to DRalgo 1.5.0 (issue #37, and beyond the scalar Yukawa model):
	- the m^4 term of the one-loop diagram had the wrong sign;
	- the fermion masses at O(m^2 g^2) were missing;
	- the fermion loop with crossed scalar lines was missing (it vanishes for complex
	  scalars with holomorphic Yukawa couplings, as in the SM);
	- the running from \[Mu] to \[Mu]3 and the pole of the soft two-loop diagrams times the O(\[Epsilon])
	  parts of the LO masses were missing. For QCD this is the term -d_A C_A \[Alpha]_E5/4
	  of hep-ph/0211321, so that PrintPressure["NNLO"] plus PrintPressureUS["NLO"] did not
	  agree with the g^4 term of Arnold and Zhai, hep-ph/9410360.

	Checked here:
	- QCD with one quark flavour: hard plus soft pressure agrees with hep-ph/9410360;
	- the SM with the top Yukawa coupling: RG invariance, and the y_t^4, y_t^2 \[Lambda] and
	  y_t^2 \[Nu]^2 terms of Gynther and Vepsalainen, hep-ph/0510375, appendix A;
	- a U(1) with a charged Dirac fermion of mass m\[Psi] and a neutral real scalar: RG
	  invariance, which tests the fermion masses with gauge and Yukawa couplings.
*)


testList={};

L=Log[\[Mu]/(4\[Pi] T)];
z1=Derivative[1][Zeta][-1]/Zeta[-1];
z3=Derivative[1][Zeta][-3]/Zeta[-3];
toL[expr_]:=Expand[expr/.{Lb->2L+2EulerGamma,Lf->2L+2EulerGamma+4Log[2]}];

(*Explicit \[Mu]-dependence at fixed \[Mu]3, and running with rules coupling->beta*)
explicitMu[expr_]:=Simplify[\[Mu] D[expr/.{Lb->2Log[\[Mu]/(4\[Pi] T)]+2EulerGamma,Lf->2Log[\[Mu]/(\[Pi] T)]+2EulerGamma},\[Mu]]];
running[expr_,rules_]:=Sum[r[[2]] D[expr,r[[1]]],{r,rules}];
(*BetaFunctions4D returns the gauge couplings as g^2->beta(g^2)*)
betaRules[]:=BetaFunctions4D[]/.(Power[c_,2]->b_):>(c->b/(2c));


(* ::Section:: *)
(*QCD with one quark flavour: hard plus soft pressure*)


Group={"SU3"};
CouplingName={gs};
RepAdjoint={{1,1}};
RepScalar={};
RepFermion={{{{1,0}},"L"},{{{1,0}},"R"}};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PerformDRhard[];
pHard=PrintPressure["NNLO"];
pNLO=PrintPressure["NLO"];
debyeLO=PrintDebyeMass["LO"];
(*Without scalars, the soft step issues messages from the empty scalar sector*)
Quiet[PerformDRsoft[{}]];
(*The A0^4 couplings contribute at O(g^6)*)
pSoft=PrintPressureUS["NLO"]/.debyeLO/.{gs3d->gs Sqrt[T],\[Lambda]VLL[_]->0};

(*The free energy F=-p of hep-ph/9410360 at O(g^4), for SU(3) with one Dirac quark*)
dA=8; CA=3; SF=1/2; S2F=2/3;
pAZ=-dA T^4 \[Pi]^2/9 (gs/(4\[Pi]))^4 (
	-48 CA (CA+SF) Log[gs/(2\[Pi]) Sqrt[(CA+SF)/3]]
	+CA^2 (22/3 L+38/3 z3-148/3 z1-4EulerGamma+64/5)
	+CA SF (47/3 L+1/3 z3-74/3 z1-8EulerGamma+1759/60+37/5 Log[2])
	+SF^2 (-20/3 L+8/3 z3-16/3 z1-4EulerGamma-1/3+88/5 Log[2])
	+S2F (-105/4+24 Log[2]));

AppendTo[testList,
VerificationTest[
	FullSimplify[PowerExpand[toL[pHard+T pSoft-pAZ]/.\[Mu]3->\[Mu]],{gs>0,T>0,\[Mu]>0}],
	0,
	TestID->"QCD-hard-plus-soft-pressure-agrees-with-Arnold-Zhai"
]];

(*The soft pressure contains -24 g3d^2 \[Mu]sqSU3 Log[\[Mu]3]/(16\[Pi]^2); the hard pressure cancels it*)
AppendTo[testList,
VerificationTest[
	Simplify[\[Mu]3 D[pHard+T pSoft,\[Mu]3]],
	0,
	TestID->"QCD-pressure-independent-of-mu3"
]];

(*Two-loop running of g^2 with one Dirac quark, beta(g^2) = -2(11-2/3) g^4/(16\[Pi]^2)*)
AppendTo[testList,
VerificationTest[
	Simplify[running[pNLO,{gs->-(11-2/3) gs^3/(16\[Pi]^2)}]+explicitMu[pHard]],
	0,
	TestID->"QCD-hard-pressure-is-rg-invariant"
]];


(* ::Section:: *)
(*Standard Model with the top Yukawa coupling*)


Group={"SU3","SU2","U1"};
RepAdjoint={{1,1},{2},0};
RepScalar={{{{0,0},{1},1/2},"C"}};
CouplingName={gs,gw,gY};
RepFermion={{{{1,0},{1},1/6},"L"},{{{1,0},{0},2/3},"R"},{{{1,0},{0},-1/3},"R"},{{{0,0},{1},-1/2},"L"},{{{0,0},{0},-1},"R"}};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
\[Phi]sq=CreateInvariant[Group,RepScalar,{{1,1},{True,False}}][[1]]//Simplify;
\[Mu]ij=GradMass[m2 \[Phi]sq]//Simplify//SparseArray;
\[Lambda]4=GradQuartic[\[Lambda]1H \[Phi]sq^2];
Ysff=-GradYukawa[yt CreateInvariantYukawa[Group,RepScalar,RepFermion,{{1,1,2},{False,False,True}}][[1]]];
YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt>0}]];
ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PerformDRhard[];
pNNLO=PrintPressure["NNLO"];

(*Four real scalars of mass m2: d\[CapitalLambda]/dln\[Mu] = 4 m2^2/2/(16\[Pi]^2)*)
AppendTo[testList,
VerificationTest[
	Simplify[running[PrintPressure["NLO"],betaRules[]]+explicitMu[pNNLO]-2m2^2/(16\[Pi]^2)],
	0,
	TestID->"SM-hard-pressure-is-rg-invariant"
]];

(*The coefficient of a monomial in the couplings, all other couplings set to zero*)
smCoefficient[mono_]:=Coefficient[toL[pNNLO]/.{gs->0,gw->0,gY->0},mono]/.{yt->0,\[Lambda]1H->0,m2->0};

(*\[Alpha]_EYY, \[Alpha]_EY\[Lambda] and \[Alpha]_EY\[Nu] of hep-ph/0510375, appendix A, with g_Y=yt, \[Lambda]=\[Lambda]1H, \[Nu]^2=-m2 and Nc=3*)
Nc=3;
aEYY=-1/32 Nc (L-239/120-11/5 Log[2]+2z1-z3-Nc (10/9 L+53/90-106/45 Log[2]+4/9 EulerGamma+4/3 z1-2/3 z3));
aEYl=1/6 Nc (L-Log[2]+EulerGamma);
aEYnu=-1/3 Nc (L-Log[2]+EulerGamma);

AppendTo[testList,
VerificationTest[
	FullSimplify[PowerExpand[{
		smCoefficient[yt^4]-T^4 aEYY/(16\[Pi]^2),
		smCoefficient[yt^2 \[Lambda]1H]-T^4 aEYl/(16\[Pi]^2),
		smCoefficient[yt^2 m2]+T^2 aEYnu/(16\[Pi]^2)}],{T>0,\[Mu]>0}],
	{0,0,0},
	TestID->"SM-top-Yukawa-terms-agree-with-Gynther-Vepsalainen"
]];

(*-m^4/4 Sum[1/P^4] for each of the four real scalars; same sign as \[Alpha]_E\[Nu]\[Nu] of hep-ph/0510375*)
AppendTo[testList,
VerificationTest[
	Simplify[smCoefficient[m2^2]],
	(L+EulerGamma)/(8\[Pi]^2),
	TestID->"SM-scalar-mass-squared-term-of-one-loop-diagram"
]];


(* ::Section:: *)
(*U(1) with a charged Dirac fermion of mass m\[Psi] and a neutral real scalar*)


Group={"U1"};
RepAdjoint={0};
CouplingName={g};
RepScalar={{{0},"R"}};
RepFermion={{{q},"L"},{{q},"R"}};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
\[Phi]=CreateInvariant[Group,RepScalar,{{1},{True}}][[1]];
\[Mu]ij=GradMass[msq/2 \[Phi]^2]//SparseArray;
\[Lambda]3=GradCubic[\[Gamma]/6 \[Phi]^3];
\[Lambda]4=GradQuartic[\[Lambda]/24 \[Phi]^4];
\[Mu]IJ=m\[Psi] GradMassFermion[CreateInvariantFermion[Group,RepFermion,{{2,1},{False,True}}][[1]]];
\[Mu]IJC=m\[Psi] GradMassFermion[CreateInvariantFermion[Group,RepFermion,{{1,2},{False,True}}][[1]]];
Ysff=y GradYukawa[CreateInvariantYukawa[Group,RepScalar,RepFermion,{{1,2,1},{True,False,True}}][[1]]];
YsffC=y GradYukawa[CreateInvariantYukawa[Group,RepScalar,RepFermion,{{1,1,2},{True,False,True}}][[1]]];
ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PerformDRhard[];

(*
	BetaFunctions4D does not return beta(g^2) without charged scalars; for one Dirac fermion
	of charge q it is 8/3 q^2 g^4/(16\[Pi]^2). The vacuum energy runs with (msq^2/2-2m\[Psi]^4)/(16\[Pi]^2).
*)
AppendTo[testList,
VerificationTest[
	Simplify[running[PrintPressure["NLO"],Append[betaRules[],g->4/3 q^2 g^3/(16\[Pi]^2)]]
		+explicitMu[PrintPressure["NNLO"]]-(msq^2/2-2m\[Psi]^4)/(16\[Pi]^2)],
	0,
	TestID->"charged-massive-fermion-hard-pressure-is-rg-invariant"
]];


(* ::Section:: *)
(*Report*)


report=TestReport[testList];
Print[report];
Print[Column[{#["TestID"],#["Outcome"]}&/@Values[report["TestResults"]]]];
Exit[If[report["AllTestsSucceeded"],0,1]]
