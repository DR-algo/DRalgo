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
(*Counterterms of the 4d theory with renormalized fields*)


(*
	Issue #47: CounterTerms4D[] gives the counterterms \[Delta]x=\[Beta]x/(2\[Epsilon]) of the bare couplings, in which
	the field renormalization is absorbed, and no field counterterms. Usually, see 1802.05056,
	appendix C.2, the fields are renormalized, \[Phi]_b=(1+\[Delta]Z)^(1/2)\[Phi], and the couplings as
	\[Lambda]_b=Z^(-2)\[Mu]^(2\[Epsilon])(\[Lambda]+\[Delta]\[Lambda]), so that \[Delta]x does not contain the field counterterms \[Delta]Z/2 x of the legs.

	CounterTerms4D["S"|"F"|"V",{pos1,pos2}] gives the field counterterms \[Delta]Z, and
	CounterTerms4D[FieldRenormalization->True] the counterterms with renormalized fields, both in
	Landau gauge. The gauge couplings carry no Z factor, g_b=\[Mu]^\[Epsilon](g+\[Delta]g).

	Checked here:
	- the 2HDM of 1802.05056 with real couplings, Nd=2 doublets, the top quark coupled to \[Phi]2, and nF
	  generations: \[Delta]Z_A, \[Delta]Z_B, \[Delta]Z_\[Phi]1, \[Delta]Z_\[Phi]2, \[Delta]Z_q, \[Delta]Z_t, \[Delta]g, \[Delta]g', \[Delta]g_Y, \[Delta]\[Mu]^2, and \[Delta]\[Lambda] of appendix C.2;
	  the default output against the \[Beta]-functions of appendix C.2 and BetaFunctions4D[];
	  the relation \[Delta]x(FieldRenormalization->True) = \[Delta]x + \[Sum]_legs \[Delta]Z/2 x;
	- a 2HDM with the top quark coupled to both doublets, which mixes the doublets, \[Delta]Z_12!=0:
	  the same relation for the masses and the Yukawa couplings;
	- the real scalar and Dirac fermion of examples/ScalYukawa.m: \[Delta]Z of Yukawa theory, and the
	  same relation for the tadpole, the scalar mass, the cubic and quartic couplings, the Yukawa
	  coupling, and the fermion mass.

	The notation of 1802.05056 in terms of the couplings here:
	g=gw, g'=gY, g_s=gs, g_Y=yt, \[Mu]^2_11=m11, \[Mu]^2_22=m22, \[Mu]^2_12=m12, \[Lambda]_i=\[Lambda]iH, and N_f=nF, see DefineNF.
	The potential and the Yukawa coupling are those of eq. (2.3) and the Lagrangian of section 2.
	The hypercharges of 1802.05056 are twice those here, Y_\[Phi]=1 with Q=I_3+Y/2, and the gauge coupling
	is g' Y/2 instead of gY Y here; the gauge coupling g'=gY is the same.
*)


testList={};
\[Kappa]=1/(16\[Pi]^2);

(*\[Delta]x(FieldRenormalization->True) - \[Delta]x - (\[Sum]_legs \[Delta]Z/2) x for the rules x->(\[Sum]_legs \[Delta]Z/2) x*)
legResidual[ct_,ctFR_,legs_]:=Simplify[(#[[1]]/.ctFR)-(#[[1]]/.ct)-#[[2]]]&/@legs;


(* ::Section:: *)
(*2HDM of 1802.05056*)


Group={"SU3","SU2","U1"};
RepAdjoint={{1,1},{2},0};
RepScalar={{{{0,0},{1},1/2},"C"},{{{0,0},{1},1/2},"C"}};
CouplingName={gs,gw,gY};
RepFermion={
	{{{1,0},{1},1/6},"L"},
	{{{1,0},{0},2/3},"R"},
	{{{1,0},{0},-1/3},"R"},
	{{{0,0},{1},-1/2},"L"},
	{{{0,0},{0},-1},"R"}
	};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];

(*p[i,j] = \[Phi]i^+ \[Phi]j*)
p[i_,j_]:=CreateInvariant[Group,RepScalar,{{j,i},{True,False}}][[1]]//Simplify;
\[Mu]ij=GradMass[m11 p[1,1]+m22 p[2,2]+m12 (p[1,2]+p[2,1])]//Simplify//SparseArray;
\[Lambda]4=GradQuartic[
	+\[Lambda]1H p[1,1]^2+\[Lambda]2H p[2,2]^2+\[Lambda]3H p[1,1]p[2,2]+\[Lambda]4H p[1,2]p[2,1]
	+\[Lambda]5H/2 (p[1,2]^2+p[2,1]^2)+\[Lambda]6H p[1,1](p[1,2]+p[2,1])+\[Lambda]7H p[2,2](p[1,2]+p[2,1])];
Ysff=-GradYukawa[yt CreateInvariantYukawa[Group,RepScalar,RepFermion,{{2,1,2},{False,False,True}}][[1]]];
YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt>0}]];

ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
(*nF generations in the gauge-field self-energies; only one top quark couples to \[Phi]2, as in 1802.05056*)
DefineNF[Table[{nF,pos},{pos,PrintFermionRepPositions[]}]];

{pos\[Phi]1,pos\[Phi]2}=PrintScalarRepPositions[];
{posC,posA,posB}=PrintGaugeRepPositions[];
{posQ,posU,posD,posL,posE}=PrintFermionRepPositions[];

ct=CounterTerms4D[];
ctFR=CounterTerms4D[FieldRenormalization->True];
\[Delta]Z\[Phi]1=CounterTerms4D["S",{pos\[Phi]1,pos\[Phi]1}];
\[Delta]Z\[Phi]2=CounterTerms4D["S",{pos\[Phi]2,pos\[Phi]2}];
\[Delta]Zq=CounterTerms4D["F",{posQ,posQ}];
\[Delta]Zt=CounterTerms4D["F",{posU,posU}];


(* ::Subsection:: *)
(*Field counterterms*)


(*1802.05056, appendix C.2, with Nd=2*)
AppendTo[testList,
VerificationTest[
	Simplify[{
		CounterTerms4D["V",{posA,posA}],
		CounterTerms4D["V",{posB,posB}],
		\[Delta]Z\[Phi]1,
		\[Delta]Z\[Phi]2,
		CounterTerms4D["S",{pos\[Phi]1,pos\[Phi]2}],
		\[Delta]Zq,
		\[Delta]Zt
		}],
	Simplify[{
		gw^2/(16\[Pi]^2 \[Epsilon]) ((26-2)/6-4/3 nF),
		-gY^2/(96\[Pi]^2 \[Epsilon]) (2+40/3 nF),
		1/(16\[Pi]^2 \[Epsilon]) (9/4 gw^2+3/4 gY^2),
		1/(16\[Pi]^2 \[Epsilon]) (9/4 gw^2+3/4 gY^2-3 yt^2),
		0,
		-yt^2/(32\[Pi]^2 \[Epsilon]),
		-yt^2/(16\[Pi]^2 \[Epsilon])
		}],
	TestID->"2hdm-field-counterterms-agree-with-1802.05056"
]];


(*In DRalgo's convention, Z=1+\[Gamma]/\[Epsilon] with the anomalous dimension \[Gamma] of AnomDim4D*)
AppendTo[testList,
VerificationTest[
	Simplify[{
		CounterTerms4D["V",{posC,posC}]-AnomDim4D["V",{posC,posC}]/\[Epsilon],
		CounterTerms4D["F",{posQ,posQ}]-AnomDim4D["F",{posQ,posQ}]/\[Epsilon],
		CounterTerms4D["S",{pos\[Phi]2,pos\[Phi]2}]-AnomDim4D["S",{pos\[Phi]2,pos\[Phi]2}]/\[Epsilon]
		}],
	{0,0,0},
	TestID->"field-counterterms-are-anomalous-dimensions-over-epsilon"
]];


(* ::Subsection:: *)
(*Counterterms of the bare couplings*)


(*The default output is unchanged by the option, and \[Delta]x=\[Beta]x/(2\[Epsilon]) of BetaFunctions4D[]*)
AppendTo[testList,
VerificationTest[
	CounterTerms4D[FieldRenormalization->False],
	ct,
	TestID->"2hdm-default-counterterms-do-not-depend-on-the-option"
]];


AppendTo[testList,
VerificationTest[
	Simplify[(#[[1]]/.ct)-#[[2]]/(2\[Epsilon])&/@BetaFunctions4D[]],
	ConstantArray[0,Length[ct]],
	TestID->"2hdm-default-counterterms-are-beta-over-2-epsilon"
]];


(*The \[Beta]-functions of 1802.05056, appendix C.2, divided by 2\[Epsilon]; \[Beta](yt)=\[Beta](yt^2)/(2yt)*)
paperBeta={
	gw^2->-gw^4/(8\[Pi]^2) (22/3-2/6-4/3 nF),
	gY^2->gY^4/(8\[Pi]^2) (2/6+20/9 nF),
	yt->1/(2yt) yt^2/(8\[Pi]^2) (9/2 yt^2-9/4 gw^2-17/12 gY^2-8 gs^2),
	m11->\[Kappa] (-3 m11 (3/2 gw^2+1/2 gY^2-4\[Lambda]1H)+2 m22 (2\[Lambda]3H+\[Lambda]4H)+12 \[Lambda]6H m12),
	m22->\[Kappa] (-3 m22 (3/2 gw^2+1/2 gY^2-2 yt^2-4\[Lambda]2H)+2 m11 (2\[Lambda]3H+\[Lambda]4H)+12 \[Lambda]7H m12),
	m12->\[Kappa] (-3 m12 (3/2 gw^2+1/2 gY^2-yt^2)+6 (m11 \[Lambda]6H+m22 \[Lambda]7H)+6 \[Lambda]5H m12+2 (2\[Lambda]4H+\[Lambda]3H) m12),
	\[Lambda]1H->\[Kappa]/2 (48 \[Lambda]1H^2+4 \[Lambda]3H^2+4 \[Lambda]3H \[Lambda]4H+2 \[Lambda]4H^2+2 \[Lambda]5H^2+24 \[Lambda]6H^2+3/4 (3 gw^4+gY^4+2 gw^2 gY^2)
		-6 \[Lambda]1H (3 gw^2+gY^2)),
	\[Lambda]2H->\[Kappa]/2 (48 \[Lambda]2H^2+4 \[Lambda]3H^2+4 \[Lambda]3H \[Lambda]4H+2 \[Lambda]4H^2+2 \[Lambda]5H^2+24 \[Lambda]7H^2+3/4 (3 gw^4+gY^4+2 gw^2 gY^2)
		-12 yt^4-6 \[Lambda]2H (3 gw^2+gY^2-4 yt^2)),
	\[Lambda]3H->2\[Kappa] (2 (\[Lambda]1H+\[Lambda]2H) (3 \[Lambda]3H+\[Lambda]4H)+2 \[Lambda]3H^2+\[Lambda]4H^2+\[Lambda]5H^2+3/8 (3 gw^4+gY^4-2 gw^2 gY^2)
		+2 (\[Lambda]6H^2+\[Lambda]7H^2)+8 \[Lambda]6H \[Lambda]7H-3/2 \[Lambda]3H (3 gw^2+gY^2-2 yt^2)),
	\[Lambda]4H->2\[Kappa] (2 (\[Lambda]1H+\[Lambda]2H) \[Lambda]4H+2 \[Lambda]4H^2+4 \[Lambda]3H \[Lambda]4H+4 \[Lambda]5H^2+3/2 gw^2 gY^2
		+5 (\[Lambda]6H^2+\[Lambda]7H^2)+2 \[Lambda]6H \[Lambda]7H-3/2 \[Lambda]4H (3 gw^2+gY^2-2 yt^2)),
	\[Lambda]5H->2\[Kappa] (2 (\[Lambda]1H+\[Lambda]2H+2 \[Lambda]3H+3 \[Lambda]4H) \[Lambda]5H+5 (\[Lambda]6H^2+\[Lambda]7H^2)+2 \[Lambda]6H \[Lambda]7H
		-3/2 \[Lambda]5H (3 gw^2+gY^2-2 yt^2)),
	\[Lambda]6H->2\[Kappa] (12 \[Lambda]1H \[Lambda]6H+(3 \[Lambda]3H+2 \[Lambda]4H) \[Lambda]7H+\[Lambda]5H \[Lambda]7H+(3 \[Lambda]3H+4 \[Lambda]4H) \[Lambda]6H+5 \[Lambda]5H \[Lambda]6H
		-3/2 \[Lambda]6H (3 gw^2+gY^2-yt^2)),
	\[Lambda]7H->2\[Kappa] (12 \[Lambda]2H \[Lambda]7H+(3 \[Lambda]3H+2 \[Lambda]4H) \[Lambda]6H+\[Lambda]5H \[Lambda]6H+(3 \[Lambda]3H+4 \[Lambda]4H) \[Lambda]7H+5 \[Lambda]5H \[Lambda]7H
		-3/2 \[Lambda]7H (3 gw^2+gY^2-3 yt^2))
	};

AppendTo[testList,
VerificationTest[
	Simplify[(#[[1]]/.ct)-#[[2]]/(2\[Epsilon])&/@paperBeta],
	ConstantArray[0,Length[paperBeta]],
	TestID->"2hdm-default-counterterms-agree-with-1802.05056-beta-functions"
]];


(* ::Subsection:: *)
(*Counterterms with renormalized fields*)


(*1802.05056, appendix C.2, with Nd=2 and \[Delta]g^2=2g \[Delta]g*)
paperCT={
	gw^2->2gw (-gw^3/(16\[Pi]^2 \[Epsilon]) ((44-2)/12-2/3 nF)),
	gY^2->2gY (gY^3/(192\[Pi]^2 \[Epsilon]) (2+40/3 nF)),
	yt->-yt/(16\[Pi]^2 \[Epsilon]) (1/3 gY^2+4 gs^2),
	m11->\[Kappa]/\[Epsilon] (6 \[Lambda]1H m11+(2\[Lambda]3H+\[Lambda]4H) m22+6 \[Lambda]6H m12),
	m22->\[Kappa]/\[Epsilon] (6 \[Lambda]2H m22+(2\[Lambda]3H+\[Lambda]4H) m11+6 \[Lambda]7H m12),
	m12->\[Kappa]/\[Epsilon] (3 \[Lambda]6H m11+3 \[Lambda]7H m22+3 \[Lambda]5H m12+(2\[Lambda]4H+\[Lambda]3H) m12),
	\[Lambda]1H->\[Kappa]/\[Epsilon]/2 (24 \[Lambda]1H^2+2 \[Lambda]3H^2+2 \[Lambda]3H \[Lambda]4H+\[Lambda]4H^2+\[Lambda]5H^2+12 \[Lambda]6H^2+3/8 (3 gw^4+gY^4+2 gw^2 gY^2)),
	\[Lambda]2H->\[Kappa]/\[Epsilon]/2 (24 \[Lambda]2H^2+2 \[Lambda]3H^2+2 \[Lambda]3H \[Lambda]4H+\[Lambda]4H^2+\[Lambda]5H^2+12 \[Lambda]7H^2+3/8 (3 gw^4+gY^4+2 gw^2 gY^2)-6 yt^4),
	\[Lambda]3H->\[Kappa]/\[Epsilon] (2 (\[Lambda]1H+\[Lambda]2H) (3 \[Lambda]3H+\[Lambda]4H)+2 \[Lambda]3H^2+\[Lambda]4H^2+\[Lambda]5H^2+3/8 (3 gw^4+gY^4-2 gw^2 gY^2)
		+2 (\[Lambda]6H^2+\[Lambda]7H^2)+8 \[Lambda]6H \[Lambda]7H),
	\[Lambda]4H->\[Kappa]/\[Epsilon] (2 (\[Lambda]1H+\[Lambda]2H) \[Lambda]4H+2 \[Lambda]4H^2+4 \[Lambda]3H \[Lambda]4H+4 \[Lambda]5H^2+3/2 gw^2 gY^2+5 (\[Lambda]6H^2+\[Lambda]7H^2)+2 \[Lambda]6H \[Lambda]7H),
	\[Lambda]5H->\[Kappa]/\[Epsilon] (2 (\[Lambda]1H+\[Lambda]2H+2 \[Lambda]3H+3 \[Lambda]4H) \[Lambda]5H+5 (\[Lambda]6H^2+\[Lambda]7H^2)+2 \[Lambda]6H \[Lambda]7H),
	\[Lambda]6H->\[Kappa]/\[Epsilon] (12 \[Lambda]1H \[Lambda]6H+(3 \[Lambda]3H+2 \[Lambda]4H) \[Lambda]7H+\[Lambda]5H \[Lambda]7H+(3 \[Lambda]3H+4 \[Lambda]4H) \[Lambda]6H+5 \[Lambda]5H \[Lambda]6H),
	\[Lambda]7H->\[Kappa]/\[Epsilon] (12 \[Lambda]2H \[Lambda]7H+(3 \[Lambda]3H+2 \[Lambda]4H) \[Lambda]6H+\[Lambda]5H \[Lambda]6H+(3 \[Lambda]3H+4 \[Lambda]4H) \[Lambda]7H+5 \[Lambda]5H \[Lambda]7H)
	};

AppendTo[testList,
VerificationTest[
	Simplify[(#[[1]]/.ctFR)-#[[2]]&/@paperCT],
	ConstantArray[0,Length[paperCT]],
	TestID->"2hdm-counterterms-with-renormalized-fields-agree-with-1802.05056"
]];


(*The gauge couplings carry no Z factor*)
AppendTo[testList,
VerificationTest[
	Simplify[{gs^2,gw^2,gY^2}/.ctFR],
	Simplify[{gs^2,gw^2,gY^2}/.ct],
	TestID->"2hdm-gauge-counterterms-do-not-depend-on-the-option"
]];


AppendTo[testList,
VerificationTest[
	Sort[Keys[ctFR]],
	Sort[Keys[ct]],
	TestID->"2hdm-same-couplings-with-and-without-field-renormalization"
]];


(*\[Lambda]6H multiplies (\[Phi]1^+\[Phi]1)(\[Phi]1^+\[Phi]2+h.c.), with three legs \[Phi]1 and one leg \[Phi]2*)
AppendTo[testList,
VerificationTest[
	legResidual[ct,ctFR,{
		\[Lambda]1H->2\[Delta]Z\[Phi]1 \[Lambda]1H,
		\[Lambda]2H->2\[Delta]Z\[Phi]2 \[Lambda]2H,
		\[Lambda]3H->(\[Delta]Z\[Phi]1+\[Delta]Z\[Phi]2) \[Lambda]3H,
		\[Lambda]4H->(\[Delta]Z\[Phi]1+\[Delta]Z\[Phi]2) \[Lambda]4H,
		\[Lambda]5H->(\[Delta]Z\[Phi]1+\[Delta]Z\[Phi]2) \[Lambda]5H,
		\[Lambda]6H->(3\[Delta]Z\[Phi]1+\[Delta]Z\[Phi]2)/2 \[Lambda]6H,
		\[Lambda]7H->(\[Delta]Z\[Phi]1+3\[Delta]Z\[Phi]2)/2 \[Lambda]7H,
		m11->\[Delta]Z\[Phi]1 m11,
		m22->\[Delta]Z\[Phi]2 m22,
		m12->(\[Delta]Z\[Phi]1+\[Delta]Z\[Phi]2)/2 m12,
		yt->(\[Delta]Z\[Phi]2+\[Delta]Zq+\[Delta]Zt)/2 yt
		}],
	ConstantArray[0,11],
	TestID->"2hdm-counterterms-differ-by-the-field-counterterms-of-the-legs"
]];


(* ::Section:: *)
(*2HDM with the top quark coupled to both doublets*)


(*The top loop mixes the doublets, and the field counterterms are a matrix in (\[Phi]1,\[Phi]2)*)
Ysff=-GradYukawa[
	+yt1 CreateInvariantYukawa[Group,RepScalar,RepFermion,{{1,1,2},{False,False,True}}][[1]]
	+yt2 CreateInvariantYukawa[Group,RepScalar,RepFermion,{{2,1,2},{False,False,True}}][[1]]];
YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt1>0,yt2>0}]];

ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];

ct=CounterTerms4D[];
ctFR=CounterTerms4D[FieldRenormalization->True];
\[Delta]Z11=CounterTerms4D["S",{pos\[Phi]1,pos\[Phi]1}];
\[Delta]Z22=CounterTerms4D["S",{pos\[Phi]2,pos\[Phi]2}];
\[Delta]Z12=CounterTerms4D["S",{pos\[Phi]1,pos\[Phi]2}];
\[Delta]Zq=CounterTerms4D["F",{posQ,posQ}];
\[Delta]Zt=CounterTerms4D["F",{posU,posU}];


AppendTo[testList,
VerificationTest[
	Simplify[{\[Delta]Z12,CounterTerms4D["S",{pos\[Phi]2,pos\[Phi]1}],\[Delta]Z11-(\[Delta]Z22/.{yt2->yt1})}],
	{-3yt1 yt2/(16\[Pi]^2 \[Epsilon]),-3yt1 yt2/(16\[Pi]^2 \[Epsilon]),0},
	TestID->"2hdm-mixing-field-counterterms"
]];


(*For the mass matrix \[Mu]^2: (\[Delta]Z.\[Mu]^2+\[Mu]^2.\[Delta]Z)/2. For the Yukawa coupling of \[Phi]1: \[Delta]Z11/2 yt1+\[Delta]Z12/2 yt2*)
AppendTo[testList,
VerificationTest[
	legResidual[ct,ctFR,{
		m11->\[Delta]Z11 m11+\[Delta]Z12 m12,
		m22->\[Delta]Z22 m22+\[Delta]Z12 m12,
		m12->(\[Delta]Z11+\[Delta]Z22)/2 m12+\[Delta]Z12 (m11+m22)/2,
		yt1->(\[Delta]Z11+\[Delta]Zq+\[Delta]Zt)/2 yt1+\[Delta]Z12/2 yt2,
		yt2->(\[Delta]Z22+\[Delta]Zq+\[Delta]Zt)/2 yt2+\[Delta]Z12/2 yt1
		}],
	ConstantArray[0,5],
	TestID->"2hdm-mixing-counterterms-differ-by-the-field-counterterms-of-the-legs"
]];


(* ::Section:: *)
(*Real scalar and Dirac fermion of examples/ScalYukawa.m*)


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

{pos\[Phi]}=PrintScalarRepPositions[];
{pos\[Psi]L,pos\[Psi]R}=PrintFermionRepPositions[];
ct=CounterTerms4D[];
ctFR=CounterTerms4D[FieldRenormalization->True];
\[Delta]Z\[Phi]=CounterTerms4D["S",{pos\[Phi],pos\[Phi]}];
\[Delta]Z\[Psi]L=CounterTerms4D["F",{pos\[Psi]L,pos\[Psi]L}];
\[Delta]Z\[Psi]R=CounterTerms4D["F",{pos\[Psi]R,pos\[Psi]R}];


(*Yukawa theory: the anomalous dimensions are 2y^2/(16\[Pi]^2) for \[Phi] and y^2/(32\[Pi]^2) for \[Psi], in the usual sign convention*)
AppendTo[testList,
VerificationTest[
	Simplify[{\[Delta]Z\[Phi],\[Delta]Z\[Psi]L,\[Delta]Z\[Psi]R}],
	{-2y^2/(16\[Pi]^2 \[Epsilon]),-y^2/(32\[Pi]^2 \[Epsilon]),-y^2/(32\[Pi]^2 \[Epsilon])},
	TestID->"scalar-yukawa-field-counterterms"
]];


(*The Yukawa vertex correction y^3/(16\[Pi]^2 \[Epsilon]) and the quartic (3\[Lambda]^2-48y^4)/(32\[Pi]^2 \[Epsilon]) without the leg corrections*)
AppendTo[testList,
VerificationTest[
	Simplify[{y,\[Lambda]}/.ctFR],
	{y^3/(16\[Pi]^2 \[Epsilon]),(3\[Lambda]^2-48y^4)/(32\[Pi]^2 \[Epsilon])},
	TestID->"scalar-yukawa-vertex-counterterms"
]];


AppendTo[testList,
VerificationTest[
	legResidual[ct,ctFR,{
		\[Sigma]->\[Delta]Z\[Phi]/2 \[Sigma],
		msq->\[Delta]Z\[Phi] msq,
		\[Gamma]->3\[Delta]Z\[Phi]/2 \[Gamma],
		\[Lambda]->2\[Delta]Z\[Phi] \[Lambda],
		y->(\[Delta]Z\[Phi]+\[Delta]Z\[Psi]L+\[Delta]Z\[Psi]R)/2 y,
		m\[Psi]->(\[Delta]Z\[Psi]L+\[Delta]Z\[Psi]R)/2 m\[Psi]
		}],
	ConstantArray[0,6],
	TestID->"scalar-yukawa-counterterms-differ-by-the-field-counterterms-of-the-legs"
]];


(* ::Section:: *)
(*Sextic coupling*)


(*
	The same model with c6/6! \[Phi]^6 in Mode 3. The scalar loop with one quartic and one sextic coupling gives
	15 \[Lambda] c6/(16\[Pi]^2), and the field normalization of the six legs 6\[Times]2y^2 c6/(16\[Pi]^2), as 4\[Times]2y^2 \[Lambda] in \[Beta](\[Lambda]).
*)
ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False,Mode->3];
DefineDim6[GradSextic[c6/720 \[Phi]^6]];
ct=CounterTerms4D[];
ctFR=CounterTerms4D[FieldRenormalization->True];

AppendTo[testList,
VerificationTest[
	Simplify[{(c6/.BetaFunctions4D[])-(15\[Lambda]+12y^2)c6/(16\[Pi]^2),(c6/.ct)-(15\[Lambda]+12y^2)c6/(32\[Pi]^2 \[Epsilon]),(c6/.ctFR)-(c6/.ct)-3\[Delta]Z\[Phi] c6}],
	{0,0,0},
	TestID->"sextic-beta-function-and-counterterms"
]];


(* ::Section:: *)
(*Report*)


report=TestReport[testList];
Print[report];
Print[Column[{#["TestID"],#["Outcome"]}&/@Values[report["TestResults"]]]];
Exit[If[report["AllTestsSucceeded"],0,1]]
