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
(*Two Dirac flavors with complex Yukawa and mass matrices*)


(*
	A U(1) with a neutral real scalar \[Phi] and two Dirac fermions of charge q, whose Yukawa and
	mass matrices are complex. In contrast to one flavor, the matrices do not commute, so
	that the order of the factors in the \[Beta]-functions matters.

	A constant shift \[Phi] -> \[Phi]+v maps the Lagrangian onto itself, with the shifted parameters
	\[Sigma]+msq v+\[Gamma] v^2/2+\[Lambda] v^3/6, msq+\[Gamma] v+\[Lambda] v^2/2, \[Gamma]+\[Lambda] v, and m+y v. Since v runs with
	the anomalous dimension of \[Phi], the \[Beta]-functions of the shifted parameters follow from
	the \[Beta]-functions of the unshifted ones. This relates the \[Beta]-functions of the cubic
	coupling and of the fermion masses to those of the quartic and Yukawa couplings.

	Checked here:
	- the shift relations for the tadpole, the scalar mass, the cubic coupling, and the fermion masses;
	- the \[Beta]-function of the Yukawa couplings against the general one-loop formula
	  (Machacek and Vaughn; Luo, Wang and Xiao, hep-ph/0211440), with the vertex
	  correction Ysff[[b]].YsffC[[a]].Ysff[[b]];
	- the RG invariance of the 3d couplings, the scalar mass, the tadpole, the Debye mass,
	  and the NNLO pressure, which involve the Yukawa and gauge corrections to the fermion loops;
	- the reality of the NNLO pressure.
*)


testList={};
private[name_String]:=ToExpression["DRalgo`DRalgo`Private`"<>name];


(* ::Section:: *)
(*Model*)


Group={"U1"};
RepAdjoint={0};
CouplingName={g1};
RepScalar={{{0},"R"}};
RepFermion={{{q},"L"},{{q},"R"},{{q},"L"},{{q},"R"}};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
{gvvv0,gvff0,gvss0}={gvvv,gvff,gvss};
\[Phi]=CreateInvariant[Group,RepScalar,{{1},{True}}][[1]];

(*y[[i,j]] and m[[i,j]] couple the Weyl fermions L_i (2i-1) and R_j (2j); yC and mC are their conjugates*)
y={{y11,y12},{y21,y22}};
yC={{yC11,yC12},{yC21,yC22}};
m={{m11,m12},{m21,m22}};
mC={{mC11,mC12},{mC21,mC22}};
yukawa[i_,j_]:=GradYukawa[CreateInvariantYukawa[Group,RepScalar,RepFermion,{{1,2j,2i-1},{True,False,True}}][[1]]];
yukawaC[i_,j_]:=GradYukawa[CreateInvariantYukawa[Group,RepScalar,RepFermion,{{1,2i-1,2j},{True,False,True}}][[1]]];
fermionMass[i_,j_]:=GradMassFermion[CreateInvariantFermion[Group,RepFermion,{{2j,2i-1},{False,True}}][[1]]];
fermionMassC[i_,j_]:=GradMassFermion[CreateInvariantFermion[Group,RepFermion,{{2i-1,2j},{False,True}}][[1]]];
yukawaTensor=Sum[y[[i,j]] yukawa[i,j],{i,2},{j,2}]//SparseArray;
yukawaTensorC=Sum[yC[[i,j]] yukawaC[i,j],{i,2},{j,2}]//SparseArray;
massTensor[mm_]:=Sum[mm[[i,j]] fermionMass[i,j],{i,2},{j,2}]//SparseArray;
massTensorC[mmC_]:=Sum[mmC[[i,j]] fermionMassC[i,j],{i,2},{j,2}]//SparseArray;

(*Imports the model with the scalar parameters {\[Sigma],msq,\[Gamma],\[Lambda]} and the fermion masses {m,mC}*)
import[{sig_,msq_,gam_,lam_},{mm_,mmC_}]:=ImportModelDRalgo[Group,gvvv0,gvff0,gvss0,
	GradTadpole[sig \[Phi]],GradCubic[gam/6 \[Phi]^3],GradQuartic[lam/24 \[Phi]^4],GradMass[msq/2 \[Phi]^2]//SparseArray,
	massTensor[mm],massTensorC[mmC],yukawaTensor,yukawaTensorC,Verbose->False];

(*The \[Beta]-function tensors of the tadpole, scalar mass, cubic, quartic, fermion masses, and Yukawa couplings*)
betaTensors[]:=(
	private["CounterTerm"][];
	Normal/@{private["\[Beta]lam1"],private["\[Beta]mij"],private["\[Beta]lamijk"],private["\[Beta]lamijkl"],
		private["\[Beta]muIJF"],private["\[Beta]muIJFC"],private["\[Beta]Ysij"],private["\[Beta]YsijC"],private["\[Gamma]ij"]}
);


(* ::Section:: *)
(*Shift relations of the \[Beta]-functions*)


import[{\[Sigma],msq,\[Gamma],\[Lambda]},{m,mC}];
{bT,bM,bC,bQ,bF,bFC,bY,bYC,anomalous}=betaTensors[];
import[{\[Sigma]+msq v+\[Gamma] v^2/2+\[Lambda] v^3/6,msq+\[Gamma] v+\[Lambda] v^2/2,\[Gamma]+\[Lambda] v,\[Lambda]},{m+y v,mC+yC v}];
{bTs,bMs,bCs,bQs,bFs,bFCs,bYs,bYCs,anomalousS}=betaTensors[];

(*The running of v*)
betaV=anomalous[[1,1]] v;
zeroQ[expr_]:=Union[Flatten[{Expand[expr]}]]==={0};

AppendTo[testList,
VerificationTest[
	zeroQ/@{
		bTs[[1]]-(bT[[1]]+v bM[[1,1]]+v^2/2 bC[[1,1,1]]+v^3/6 bQ[[1,1,1,1]]+(msq+\[Gamma] v+\[Lambda] v^2/2) betaV),
		bMs[[1,1]]-(bM[[1,1]]+v bC[[1,1,1]]+v^2/2 bQ[[1,1,1,1]]+(\[Gamma]+\[Lambda] v) betaV),
		bCs[[1,1,1]]-(bC[[1,1,1]]+v bQ[[1,1,1,1]]+\[Lambda] betaV),
		bFs-(bF+v bY[[1]]+Normal[yukawaTensor[[1]]] betaV),
		bFCs-(bFC+v bYC[[1]]+Normal[yukawaTensorC[[1]]] betaV)},
	{True,True,True,True,True},
	TestID->"beta-functions-satisfy-shift-relations"
]];

(*General one-loop formula, in units of 1/(16\[Pi]^2), with Y2 = Sum[YsffC[[b]].Ysff[[b]]]*)
{yP,yCP,gP}=Normal/@{yukawaTensor,yukawaTensorC,gvff0};
Y2=Sum[yCP[[l]] . yP[[l]],{l,Length[yP]}];
betaYukawa[k_]:=(
	(1/2 (Transpose[Y2] . yP[[k]]+yP[[k]] . Y2)
	+2 Sum[yP[[l]] . yCP[[k]] . yP[[l]],{l,Length[yP]}]
	+6 Sum[gP[[n]] . yP[[k]] . Transpose[gP[[n]]],{n,Length[gP]}])/(16\[Pi]^2)
	-Sum[anomalous[[k,l]] yP[[l]],{l,Length[yP]}]);

AppendTo[testList,
VerificationTest[
	zeroQ[bY-Table[betaYukawa[k],{k,Length[yP]}]],
	True,
	TestID->"yukawa-beta-agrees-with-general-one-loop-formula"
]];


(* ::Section:: *)
(*RG invariance*)


import[{\[Sigma],msq,\[Gamma],\[Lambda]},{m,mC}];
PerformDRhard[];

(*
	The \[Beta]-functions from the tensors, since BetaFunctions4D does not solve for the independent
	parameters m and mC, and y and yC. For two Dirac fermions of charge q, \[Beta](g1) = 8/3 q^2 g1^3/(16\[Pi]^2).
*)
{bT,bM,bC,bQ,bF,bFC,bY,bYC,anomalous}=betaTensors[];
{yT,yCT,mT,mCT}=Normal/@{yukawaTensor,yukawaTensorC,massTensor[m],massTensorC[mC]};
component[tensor_,par_,i_,j_]:=tensor[[2i-1,2j]]/D[{yT[[1]],yCT[[1]],mT,mCT}[[par]][[2i-1,2j]],{y,yC,m,mC}[[par,i,j]]];
betas=Join[
	Flatten@Table[{
		y[[i,j]]->component[bY[[1]],1,i,j],yC[[i,j]]->component[bYC[[1]],2,i,j],
		m[[i,j]]->component[bF,3,i,j],mC[[i,j]]->component[bFC,4,i,j]},{i,2},{j,2}],
	{\[Sigma]->bT[[1]],msq->bM[[1,1]],\[Gamma]->bC[[1,1,1]],\[Lambda]->bQ[[1,1,1,1]],g1->8/3 q^2 g1^3/(16\[Pi]^2)}];

explicitMu[expr_]:=Simplify[mu4 D[expr/.{Lb->2Log[mu4/(4\[Pi] T)]+2EulerGamma,Lf->2Log[mu4/(\[Pi] T)]+2EulerGamma},mu4]];
running[expr_]:=Sum[r[[2]] D[expr,r[[1]]],{r,betas}];
lo3d={\[Lambda]3d->T \[Lambda],\[Gamma]3d->Sqrt[T] \[Gamma]};
rgResidual[lo_,nlo_]:=Expand[running[lo]+explicitMu[nlo]/.lo3d];
couplings=PrintCouplings[];

AppendTo[testList,
VerificationTest[
	{
		rgResidual[\[Lambda]3d/.couplings/.{Lb->0,Lf->0},\[Lambda]3d/.couplings],
		rgResidual[\[Gamma]3d/.couplings/.{Lb->0,Lf->0},\[Gamma]3d/.couplings],
		rgResidual[PrintTadpoles["LO"][[1,2]],PrintTadpoles["NLO"][[1,2]]],
		rgResidual[PrintDebyeMass["LO"][[1,2]],PrintDebyeMass["NLO"][[1,2]]]},
	{0,0,0,0},
	TestID->"couplings-tadpole-and-debye-mass-are-rg-invariant"
]];

(*
	The running from mu4 to mu3 contains the square of the one-loop coupling lamVL. This term of
	O(g1^4 y^4/(16\[Pi]^2)^3) lies beyond the order computed and remains in the residual.
*)
AppendTo[testList,
VerificationTest[
	FreeQ[Expand[rgResidual[PrintScalarMass["LO"][[1,2]],PrintScalarMass["NLO"][[1,2]]]\[Pi]^6],\[Pi]],
	True,
	TestID->"scalar-mass-is-rg-invariant-at-the-order-computed"
]];

(*The vacuum energy runs with (msq^2/2-Tr[(\[Mu]IJF.\[Mu]IJFC)^2])/(16\[Pi]^2)*)
mmC=mT . mCT;
pressureNNLO=PrintPressure["NNLO"];
AppendTo[testList,
VerificationTest[
	Expand[running[PrintPressure["NLO"]]+explicitMu[pressureNNLO]-(msq^2/2-Tr[mmC . mmC])/(16\[Pi]^2)],
	0,
	TestID->"NNLO-pressure-is-rg-invariant"
]];

conjugate[expr_]:=expr/.Join[Thread[Flatten[{y,yC,m,mC}]->Flatten[{yC,y,mC,m}]],{Complex[re_,im_]:>Complex[re,-im]}];
AppendTo[testList,
VerificationTest[
	Expand[conjugate[pressureNNLO]-pressureNNLO],
	0,
	TestID->"NNLO-pressure-is-real"
]];


(* ::Section:: *)
(*Report*)


report=TestReport[testList];
Print[report];
Print[Column[{#["TestID"],#["Outcome"]}&/@Values[report["TestResults"]]]];
Exit[If[report["AllTestsSucceeded"],0,1]]
