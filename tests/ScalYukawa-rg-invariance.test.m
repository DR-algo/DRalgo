(* ::Package:: *)

If[$InputFileName=="",
	SetDirectory[NotebookDirectory[]],
	SetDirectory[DirectoryName[$InputFileName]]
];
(*Put this if you want to create multiple model-files with the same kernel*)
(*DRalgo`DRalgo`$GroupMathMultipleModels=True;*)

DRalgo`DRalgo`$LoadGroupMath=True;
(*DRalgo`DRalgo`$InstallGroupMath=True;*)

Check[
    Get["../Kernel/DRalgo.m"],
    Message[Get::noopen, "DRalgo` at "<>ToString[$UserBaseDirectory]<>"/Applications"];
    Abort[];
]


(* ::Chapter:: *)
(*Real scalar + Dirac fermion: RG invariance of the matching relations and of the pressure*)


(*
	The model of examples/ScalYukawa.m (see 2108.04377 and 2310.02308).

	The parameters of the 3d theory do not depend on the 4d scale \[Mu] at the order
	computed: the explicit \[Mu]-dependence of the NLO terms, through Lb and Lf, cancels
	the running of the LO terms with BetaFunctions4D[]. The 3d scale \[Mu]3 is held
	fixed. The same holds for the hard pressure, together with the running of the
	vacuum energy, which DRalgo does not track.

	Up to DRalgo 1.5.0 this failed in four places (issues #36, #37, #38):
	- the fermion-mass term of \[Beta]\[Gamma] was -24 m\[Psi] y^3 instead of -48 m\[Psi] y^3;
	- the NLO tadpole lacked the field normalisation of \[Sigma] and had a wrong m\[Psi] y^3 term;
	- the LO pressure counted the photon of the auxiliary U(1);
	- the NNLO pressure had the wrong sign of msq^2, lacked the O(m\[Psi]^2 y^2) terms,
	  the fermion loop with crossed scalar lines (y^4), and the running from \[Mu] to \[Mu]3 (\[Gamma]^2).

	The tests check the RG invariance, and reference values for the \[Beta]-functions
	(2310.02308, appendix D), the NLO tadpole (the one-point function \[CapitalGamma]^(1) of 2310.02308,
	appendix B, with the field normalisation), and the LO pressure (ideal gas).
*)


(* ::Section:: *)
(*Model*)


Group={"U1"};
CouplingName={g};
RepAdjoint={0};
RepScalar={{{0},"R"}};
RepFermion={{{0},"L"},{{0},"R"}};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJC,\[Mu]IJ,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];

\[Phi]=CreateInvariant[Group,RepScalar,{{1},{True}}][[1]];
\[Lambda]1=GradTadpole[\[Sigma] \[Phi]];
\[Mu]ij=GradMass[msq/2 \[Phi]^2]//SparseArray;
\[Lambda]3=GradCubic[\[Gamma]/6 \[Phi]^3];
\[Lambda]4=GradQuartic[\[Lambda]/24 \[Phi]^4];
\[Mu]IJ=m\[Psi] GradMassFermion[CreateInvariantFermion[Group,RepFermion,{{2,1},{False,True}}][[1]]];
\[Mu]IJC=m\[Psi] GradMassFermion[CreateInvariantFermion[Group,RepFermion,{{1,2},{False,True}}][[1]]];
Ysff=y GradYukawa[CreateInvariantYukawa[Group,RepScalar,RepFermion,{{1,2,1},{True,False,True}}][[1]]];
YsffC=y GradYukawa[CreateInvariantYukawa[Group,RepScalar,RepFermion,{{1,1,2},{True,False,True}}][[1]]];

ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PerformDRhard[];

testList={};


(* ::Section:: *)
(*Helpers*)


betas=BetaFunctions4D[];

(*Explicit \[Mu]-dependence at fixed \[Mu]3*)
explicitMu[expr_]:=Simplify[\[Mu] D[expr/.{Lb->2Log[\[Mu]/(4\[Pi] T)]+2EulerGamma,Lf->2Log[\[Mu]/(\[Pi] T)]+2EulerGamma},\[Mu]]];
running[expr_]:=Sum[r[[2]] D[expr,r[[1]]],{r,betas}];

(*The 3d couplings in the logarithms of \[Mu]3 at LO*)
lo3d={\[Lambda]3d->T \[Lambda],\[Gamma]3d->Sqrt[T] \[Gamma]};
rgResidual[lo_,nlo_]:=Simplify[running[lo]+explicitMu[nlo]/.lo3d];


(* ::Section:: *)
(*\[Beta]-functions*)


(*2310.02308, appendix D, in units of 1/(16\[Pi]^2)*)
expandRules[rules_]:=Sort[rules/.(a_->b_):>(a->Expand[b])];

AppendTo[testList,
VerificationTest[
	expandRules[betas],
	expandRules[{
		\[Lambda]->(3\[Lambda]^2-48y^4+8\[Lambda] y^2)/(16\[Pi]^2),
		\[Gamma]->(3\[Gamma] \[Lambda]-48y^3 m\[Psi]+6\[Gamma] y^2)/(16\[Pi]^2),
		y->5y^3/(16\[Pi]^2),
		msq->(\[Gamma]^2+\[Lambda] msq+4y^2 msq-24y^2 m\[Psi]^2)/(16\[Pi]^2),
		m\[Psi]->3y^2 m\[Psi]/(16\[Pi]^2),
		\[Sigma]->(\[Gamma] msq-8y m\[Psi]^3+2\[Sigma] y^2)/(16\[Pi]^2)}],
	TestID->"beta-functions-agree-with-2310.02308"
]];


(* ::Section:: *)
(*RG invariance of the 3d parameters*)


couplings=PrintCouplings[];
coupling[name_]:=name/.couplings;

AppendTo[testList,
VerificationTest[
	rgResidual[coupling[\[Gamma]3d]/.{Lb->0,Lf->0},coupling[\[Gamma]3d]],
	0,
	TestID->"cubic-coupling-is-rg-invariant"
]];

AppendTo[testList,
VerificationTest[
	rgResidual[coupling[\[Lambda]3d]/.{Lb->0,Lf->0},coupling[\[Lambda]3d]],
	0,
	TestID->"quartic-coupling-is-rg-invariant"
]];

AppendTo[testList,
VerificationTest[
	rgResidual[PrintScalarMass["LO"][[1,2]],PrintScalarMass["NLO"][[1,2]]],
	0,
	TestID->"scalar-mass-is-rg-invariant"
]];

AppendTo[testList,
VerificationTest[
	rgResidual[PrintTadpoles["LO"][[1,2]],PrintTadpoles["NLO"][[1,2]]],
	0,
	TestID->"tadpole-is-rg-invariant"
]];


(* ::Section:: *)
(*NLO tadpole*)


(*
	\[CapitalGamma]^(1) of 2310.02308, appendix B, contains m\[Psi] T^2 y^3/(6(4\[Pi])^2) (Lf+2Lb), and \[Sigma]
	enters only at tree level. Dividing by Sqrt[T] and by the field normalisation
	Sqrt[1+2y^2 Lf/(4\[Pi])^2] gives the two terms below. (The field normalisation of
	2310.02308 has an additional constant 4y^2/(3(4\[Pi])^2); averaging the fermion loop over
	d-1 spatial directions, its p^2 term is proportional to Lf alone.)
*)
tadpoleNLO=Expand[PrintTadpoles["NLO"][[1,2]]];

AppendTo[testList,
VerificationTest[
	Coefficient[tadpoleNLO,\[Sigma]],
	-Lf y^2/(16\[Pi]^2 Sqrt[T]),
	TestID->"tadpole-field-normalisation-of-sigma"
]];

AppendTo[testList,
VerificationTest[
	Simplify[Coefficient[tadpoleNLO,m\[Psi] y^3]],
	Lb T^(3/2)/(48\[Pi]^2),
	TestID->"tadpole-fermion-mass-yukawa-term"
]];


(* ::Section:: *)
(*Pressure*)


(*One real scalar and one Dirac fermion; the photon of the auxiliary U(1) does not count*)
AppendTo[testList,
VerificationTest[
	PrintPressure["LO"],
	\[Pi]^2 T^4/90 (1+7/8 4),
	TestID->"LO-pressure-is-ideal-gas-of-scalar-and-Dirac-fermion"
]];

(*
	The vacuum energy \[CapitalLambda] runs with d\[CapitalLambda]/dln\[Mu] = (msq^2/2-2 m\[Psi]^4)/(16\[Pi]^2), from the
	Coleman-Weinberg potential of one real scalar and one Dirac fermion; the pressure
	contains -\[CapitalLambda].
*)
betaVacuum=(msq^2/2-2m\[Psi]^4)/(16\[Pi]^2);

AppendTo[testList,
VerificationTest[
	Simplify[running[PrintPressure["NLO"]]+explicitMu[PrintPressure["NNLO"]]-betaVacuum],
	0,
	TestID->"NNLO-pressure-is-rg-invariant"
]];

(*
	The scalar sunset of the 3d theory, -1/12 \[Gamma]3d^2 (1/(4\[Epsilon])+Log[\[Mu]3/(3 m3)]+1/2)/(16\[Pi]^2),
	contributes +T \[Gamma]3d^2 Log[\[Mu]3]/(192\[Pi]^2) to the pressure; the hard pressure cancels it.
*)
AppendTo[testList,
VerificationTest[
	Simplify[\[Mu]3 D[PrintPressure["NNLO"],\[Mu]3]],
	-T^2 \[Gamma]^2/(192\[Pi]^2),
	TestID->"NNLO-pressure-cancels-mu3-of-soft-cubic-sunset"
]];


(* ::Section:: *)
(*Report*)


report=TestReport[testList];
Print[report];
Print[Column[{#["TestID"],#["Outcome"]}&/@Values[report["TestResults"]]]];
Exit[If[report["AllTestsSucceeded"],0,1]]
