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
(*Soft to ultrasoft: integrating out a heavy scalar*)


(*
	Issue #28: integrating out the singlet of xSM.m, the ultrasoft Higgs self-coupling contains
	\[Lambda]1H3d \[Mu]m3d^2/(\[Pi] \[Mu]\[Sigma]3d^(3/2)) with the coefficient -19/48 (+1/24 in earlier versions).
	The triangle diagram gives -3/8, see Eq. (3.30) of 1609.06230, and the field normalization of the Higgs,
	1+\[Mu]m3d^2/(48 \[Pi] \[Mu]\[Sigma]3d^(3/2)), adds -1/24, so that the coefficient is -5/12.

	The matching is computed here independently of DRalgo. The one-loop effective potential is the hard region
	of 1/2 Tr ln(p^2+M^2) in a background of the light field h, with the heavy field at its tree-level minimum.
	This contains all one-point reducible diagrams with heavy lines, and dividing by powers of the field
	normalization gives the parameters of the ultrasoft theory. Light masses are kept to first order.

	Checked here:
	- one light and one heavy real scalar with all couplings: tadpole, LO and NLO mass, cubic, and quartic;
	- the scalar sector of xSM.m, the Higgs doublet and a singlet: the coefficient of issue #28, the quartic,
	  and the mass.
*)


testList={};

(*3d integral of 1/((p^2)^k (p^2+m2)^j) in dimensional regularization*)
integral[k_,j_,m2_]:=Gamma[k+j-3/2]Gamma[3/2-k]/((4Pi)^(3/2)Gamma[j]Gamma[3/2]) m2^(3/2-k-j);

(*
	Effective potential of the light background h, for V[h,s] with the light field h and the heavy field s:
	tree level at the minimum s0 = c1 h + ... + c4 h^4 and the hard region of the one-loop determinant of the
	(h,s) block. The parameter w counts insertions of the light mass.
*)
effectivePotential[V_]:=Module[{d2h,d2hs,d2s,oneLoop,sMin,coeffs,sol},
	d2h=D[V[h,s],h,h];
	d2hs=D[V[h,s],h,s];
	d2s=D[V[h,s],s,s];
	oneLoop=-d2s^(3/2)/(12Pi)
		+1/2 Sum[(-1)^(k+1)/k Binomial[k,j] d2h^(k-j) (-d2hs^2)^j integral[k,j,d2s],{k,1,6},{j,1,k}];
	sMin=c1 h+c2 h^2+c3 h^3+c4 h^4;
	coeffs=CoefficientList[Normal[Series[D[V[h,s],s]/.s->sMin,{h,0,4}]],h];
	sol=First@Solve[Thread[PadRight[coeffs,5][[2;;5]]==0],{c1,c2,c3,c4}];
	Normal[Series[Normal[Series[V[h,sMin/.sol]+eps oneLoop/.s->(sMin/.sol),{w,0,1}]],{h,0,4}]]
];

(*Expansion to one loop*)
oneLoopExpand[expr_]:=Expand[Normal[Series[expr,{eps,0,1}]]/.w->1];


(* ::Section:: *)
(*One light and one heavy real scalar with all couplings*)


Vgen[f_,s_]:=tf f+w mf/2 f^2+cf/6 f^3+lf/24 f^4+ms/2 s^2+cs/6 s^3+ls/24 s^4+ya/2 f^2 s+yb/2 f s^2+l22/4 f^2 s^2+l31/6 f^3 s+l13/6 f s^3;
Veff=effectivePotential[Vgen];
coef[n_]:=Coefficient[Veff,h,n] n!;
(*Field normalization from the bubbles with one light and one heavy line (ya) and with two heavy lines (yb)*)
dZ=ya^2/(12Pi ms^(3/2))+1/2 yb^2/(12Pi (2Sqrt[ms])^3);
expected=oneLoopExpand/@{coef[1](1-eps dZ/2),coef[2](1-eps dZ),coef[3](1-3eps dZ/2),coef[4](1-2eps dZ)};

Group={"U1"};RepAdjoint={0};CouplingName={g1};
RepScalar={{{0},"R"},{{0},"R"}};RepFermion={};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
fField=CreateInvariant[Group,RepScalar,{{1},{True}}][[1]];
sField=CreateInvariant[Group,RepScalar,{{2},{True}}][[1]];
\[Lambda]1=GradTadpole[tf fField];
\[Mu]ij=GradMass[mf/2 fField^2+ms/2 sField^2]//SparseArray;
\[Lambda]3=GradCubic[Vgen[fField,sField]-(Vgen[fField,sField]/.{cf->0,cs->0,ya->0,yb->0})/.{tf->0,w->0,ms->0}];
\[Lambda]4=GradQuartic[Vgen[fField,sField]/.{tf->0,w->0,ms->0,cf->0,cs->0,ya->0,yb->0}];
ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
Quiet[PerformDRhard[]];
PerformDRsoft[{2}];

to3d={tf->tf3d,mf->mf3d,cf->cf3d,lf->lf3d,ms->ms3d,cs->cs3d,ls->ls3d,ya->ya3d,yb->yb3d,l22->l223d,l31->l313d,l13->l133d};
couplingsUS=PrintCouplingsUS[];
massLO=mf3dUS/.PrintScalarMassUS["LO"];
massNLO=mf3dUS/.PrintScalarMassUS["NLO"];
(*The terms of O(eps^0) and O(eps) without the light mass are the LO, those with the light mass the NLO*)
{tadpoleExp,massExp,cubicExp,quarticExp}=expected/.to3d;

AppendTo[testList,
VerificationTest[
	Simplify/@{
		(tf3dUS/.PrintTadpolesUS["All"])-(tadpoleExp/.eps->1),
		massLO-(massExp/.eps->1/.mf3d->0)-mf3d,
		Coefficient[Expand[massNLO],mf3d]-Coefficient[Coefficient[massExp,eps],mf3d],
		(cf3dUS/.couplingsUS)-(cubicExp/.eps->1/.mf3d->0),
		(lf3dUS/.couplingsUS)-(quarticExp/.eps->1/.mf3d->0)},
	{0,0,0,0,0},
	TestID->"light-and-heavy-scalar-with-all-couplings"
]];


(* ::Section:: *)
(*Scalar sector of xSM.m: Higgs doublet and singlet*)


(*The four components of the doublet with phi^+phi = h^2/2 along the background, and the singlet s*)
VxSM[hh_,s_]:=w m1 hh^2/2+lh hh^4/4+mus/2 s^2+cs/3 s^3+ls/4 s^4+mum/4 s hh^2+lm/4 s^2 hh^2;
(*Only the (h,s) block has a hard region; the Goldstone modes contain light lines only*)
Veff=effectivePotential[VxSM];
dZ=(mum/2)^2/(12Pi mus^(3/2));
{massExp,quarticExp}=oneLoopExpand/@{2Coefficient[Veff,h,2](1-eps dZ),4Coefficient[Veff,h,4](1-2eps dZ)};

RepScalar=Table[{{0},"R"},{5}];
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
components=Table[CreateInvariant[Group,RepScalar,{{i},{True}}][[1]],{i,5}];
phiSq=Total[components[[1;;4]]^2]/2;
sField=components[[5]];
\[Mu]ij=GradMass[m1 phiSq+mus/2 sField^2]//SparseArray;
\[Lambda]3=GradCubic[mum/2 sField phiSq+cs/3 sField^3];
\[Lambda]4=GradQuartic[lh phiSq^2+lm/2 sField^2 phiSq+ls/4 sField^4];
ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
Quiet[PerformDRhard[]];
PerformDRsoft[{5}];

to3d={m1->m13d,lh->lh3d,mus->mus3d,cs->cs3d,ls->ls3d,mum->mum3d,lm->lm3d};
quarticUS=lh3dUS/.PrintCouplingsUS[];
{massExp,quarticExp}={massExp,quarticExp}/.to3d;

AppendTo[testList,
VerificationTest[
	Simplify[Coefficient[Expand[quarticUS],lh3d mum3d^2]Pi mus3d^(3/2)],
	-5/12,
	TestID->"xSM-coefficient-of-issue-28"
]];

AppendTo[testList,
VerificationTest[
	Simplify/@{
		quarticUS-(quarticExp/.eps->1/.m13d->0),
		(m13dUS/.PrintScalarMassUS["LO"])-(massExp/.eps->1/.m13d->0)-m13d,
		Coefficient[Expand[m13dUS/.PrintScalarMassUS["NLO"]],m13d]-Coefficient[Coefficient[massExp,eps],m13d]},
	{0,0,0},
	TestID->"xSM-quartic-and-mass"
]];


(* ::Section:: *)
(*Report*)


report=TestReport[testList];
Print[report];
Print[Column[{#["TestID"],#["Outcome"]}&/@Values[report["TestResults"]]]];
Exit[If[report["AllTestsSucceeded"],0,1]]
