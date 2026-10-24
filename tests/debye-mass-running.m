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
(*Running of the Debye masses in the soft theory*)


(*
	Issue #59: the Debye masses of the soft theory run at two loops, as the scalar masses do.
	BetaFunctions3DS gives their \[Beta]-functions from the temporal block of the two-loop mass counterterm
	of the soft theory, which includes the temporal-scalar quartics.

	The lattice mass counterterms of hep-lat/9504001, appendix, contain -C/(16\[Pi]^2) Log[6/(a\[Mu])], where C is
	the coefficient of the MSbar counterterm, so that mu3 dm^2/dmu3 = -C/(16\[Pi]^2). For the Debye mass of SU(2)
	with an adjoint and a fundamental Higgs, (A3) and (A11) give C = 20\[Lambda]A g3^2-10\[Lambda]A^2-g3^4+6h3 g3^2-8h3^2,
	and the U(1) results give C = -6\[Lambda]A^2+4h3 e3^2-4h3^2. There, the Lagrangian contains (1/4)\[Lambda]A (A0^a A0^a)^2
	and h3 \[Phi]^\[Dagger]\[Phi] A0^a A0^a with the doublet \[Phi], or h3 \[Phi]^*\[Phi] A0^2 for U(1). The soft tensors enter as
	(1/4!)\[Lambda]AAS_abcd A0^a A0^b A0^c A0^d and (1/4)\[Lambda]KVec_abij A0^a A0^b \[Phi]_i \[Phi]_j with real components \[Phi]_i,
	so that \[Lambda]A=\[Lambda]AAS_1111/6 and h3=\[Lambda]KVec_1111/2.

	At the tree-level values h3=g3^2/4 and \[Lambda]A=0, the O(g^4) part of the running vanishes, as noted
	below (A11). The running is therefore of O(g^6), beyond the order of the NLO Debye masses, which are
	RG invariant at O(g^4) and do not depend on mu3.

	Checked here:
	- SU(2)+Higgs: the normalization of h3, the \[Beta]-function against (A3) and (A11), that it starts at O(g^6),
	  the RG invariance of the NLO Debye mass at O(g^4), its independence of mu3 and AutoRG, and Mode->1;
	- scalar QED: the \[Beta]-function against the U(1) results of hep-lat/9504001, and that it starts at O(g^6);
	- two U(1) with a scalar of charges (1,1): the off-diagonal Debye masses and temporal-scalar quartics
	  are those of scalar QED with ge^2=g1^2+g2^2, rotated by (g1,g2)/ge;
	- QCD, without scalars: the \[Beta]-function against the SU(N) form of (A3), 2(dA+2)\[Lambda]A(\[Lambda]A-CA g3^2),
	  and the RG invariance of the NLO Debye mass at O(g^4).
*)


testList={};


(* ::Section:: *)
(*Helpers*)


(*mu4 d/dmu4 at fixed mu3 and the 4d running with BetaFunctions4D, which gives \[Beta] of g^2 for gauge couplings*)
explicitMu[expr_]:=Simplify[mu4 D[expr/.{Lb->2Log[mu4/(4\[Pi] T)]+2EulerGamma,Lf->2Log[mu4/(\[Pi] T)]+2EulerGamma},mu4]];
running[expr_,betas_]:=Sum[r[[2]] D[expr,r[[1]]],{r,betas/.(Power[cc_,2]->bb_):>(cc->bb/(2cc))}];

(*The coefficient of x^n, with gauge couplings of O(x) and quartics and masses of O(x^2)*)
order[expr_,n_,scaling_]:=Simplify[SeriesCoefficient[expr/.scaling,{x,0,n}]];

(*Placeholders in terms of 4d parameters*)
resolve[expr_]:=expr/.PrintIdentification[];

(*The temporal-scalar couplings in the normalization of hep-lat/9504001*)
temporalCouplings[]:={PrintTensorDRalgo[3][[1,1,1,1]]/2,PrintTensorDRalgo[4][[1,1,1,1]]/6};

(*Running of the Debye mass m in the soft theory: mu3 dm/dmu3*)
debyeBeta[m_]:=resolve[m/.BetaFunctions3DS[]];


(* ::Section:: *)
(*SU(2)+Higgs*)


importSU2Higgs[opts___]:=Module[{\[Phi]sq},
	Group={"SU2"};
	RepAdjoint={{2}};
	CouplingName={gw};
	RepScalar={{{{1}},"C"}};
	RepFermion={};
	{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
	\[Phi]sq=CreateInvariant[Group,RepScalar,{{1,1},{True,False}}][[1]]//Simplify;
	\[Mu]ij=GradMass[m2 \[Phi]sq]//Simplify//SparseArray;
	\[Lambda]4=GradQuartic[\[Lambda] \[Phi]sq^2];
	ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False,opts];
	PerformDRhard[];
];
scalingSU2={gw->x gw,\[Lambda]->x^2 \[Lambda],m2->x^2 m2};

importSU2Higgs[];
{h3,\[Lambda]A}=temporalCouplings[];
(*(A3) and (A11) of hep-lat/9504001*)
betaSU2=-1/(16\[Pi]^2)((20\[Lambda]A g3^2-10\[Lambda]A^2)+(-g3^4+6h3 g3^2-8h3^2))/.g3->gw3d;
runningSU2=debyeBeta[musqSU2];
(*With the NLO gauge coupling, the coefficient of g^6 is complete*)
runningSU2NLO=runningSU2/.gw3d->Sqrt[gw3d^2/.PrintCouplings[]];
nloSU2Output=PrintDebyeMass["NLO"];
nloSU2=resolve[musqSU2/.nloSU2Output];

(*h3=g3^2/4 at tree level, see the remark below (A11)*)
AppendTo[testList,
VerificationTest[
	order[resolve[h3],2,scalingSU2],
	gw^2 T/4,
	TestID->"SU2-temporal-scalar-coupling-in-normalization-of-hep-lat-9504001"
]];

AppendTo[testList,
VerificationTest[
	Simplify[runningSU2-resolve[betaSU2]],
	0,
	TestID->"SU2-debye-mass-runs-as-A3-and-A11-of-hep-lat-9504001"
]];

AppendTo[testList,
VerificationTest[
	{order[runningSU2NLO,4,scalingSU2],order[runningSU2NLO,6,scalingSU2]===0},
	{0,False},
	TestID->"SU2-debye-mass-running-is-of-O(g^6)"
]];

(*With the 3d couplings at LO, the NLO Debye mass cancels the running of the LO Debye mass at O(g^4)*)
AppendTo[testList,
VerificationTest[
	order[running[musqSU2/.PrintDebyeMass["LO"],BetaFunctions4D[]]+explicitMu[nloSU2/.{gw3d->gw Sqrt[T],\[Lambda]3d->\[Lambda] T}],4,scalingSU2],
	0,
	TestID->"SU2-debye-mass-is-rg-invariant-at-O(g^4)"
]];

(*The running is beyond the order of the NLO Debye mass, which therefore does not depend on mu3 or AutoRG*)
importSU2Higgs[AutoRG->False];
nloSU2NoRG=PrintDebyeMass["NLO"];

AppendTo[testList,
VerificationTest[
	{FreeQ[nloSU2Output,mu3],nloSU2NoRG===nloSU2Output},
	{True,True},
	TestID->"SU2-nlo-debye-mass-does-not-depend-on-mu3-or-AutoRG"
]];

importSU2Higgs[Mode->1];

AppendTo[testList,
VerificationTest[
	PrintDebyeMass["All"],
	{musqSU2->5/6 gw^2 T^2},
	TestID->"SU2-debye-mass-in-mode-1-is-LO"
]];


(* ::Section:: *)
(*Scalar QED*)


importU1[charges_,couplings_]:=Module[{\[Phi]sq},
	Group=Table["U1",{Length[charges]}];
	RepAdjoint=Table[0,{Length[charges]}];
	CouplingName=couplings;
	RepScalar={{charges,"C"}};
	RepFermion={};
	{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
	\[Phi]sq=CreateInvariant[Group,RepScalar,{{1,1},{True,False}}][[1]]//Simplify;
	\[Mu]ij=GradMass[m2 \[Phi]sq]//Simplify//SparseArray;
	\[Lambda]4=GradQuartic[\[Lambda] \[Phi]sq^2];
	ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
	PerformDRhard[];
];

importU1[{1},{ge}];
{h3,\[Lambda]A}=temporalCouplings[];
(*U(1)+fundamental Higgs+adjoint Higgs of hep-lat/9504001*)
betaU1=-1/(16\[Pi]^2)(-6\[Lambda]A^2+4h3 e3^2-4h3^2)/.e3->ge3d;

runningU1=debyeBeta[musqU1];
runningU1NLO=runningU1/.ge3d->Sqrt[ge3d^2/.PrintCouplings[]];

AppendTo[testList,
VerificationTest[
	Simplify[runningU1-resolve[betaU1]],
	0,
	TestID->"U1-debye-mass-runs-as-hep-lat-9504001"
]];

AppendTo[testList,
VerificationTest[
	With[{scaling={ge->x ge,\[Lambda]->x^2 \[Lambda],m2->x^2 m2}},
		{order[runningU1NLO,4,scaling],order[runningU1NLO,6,scaling]===0}],
	{0,False},
	TestID->"U1-debye-mass-running-is-of-O(g^6)"
]];

(*Results of scalar QED in terms of 4d parameters, for the comparison below*)
debyeQED=resolve[musqU1/.PrintDebyeMass["All"]];
quarticQED=resolve[PrintTensorDRalgo[4][[1,1,1,1]]];


(* ::Section:: *)
(*Two U(1) with off-diagonal Debye masses*)


(*
	For the charges (1,1), the vector (g1 B1+g2 B2)/ge with ge^2=g1^2+g2^2 couples like the photon of
	scalar QED, and the orthogonal one decouples. The same holds for the 3d gauge couplings.
*)
importU1[{1,1},{g1,g2}];
rotation={ge->Sqrt[g1^2+g2^2],ge3d->Sqrt[g13d^2+g23d^2]};
gList={g1,g2};

AppendTo[testList,
VerificationTest[
	Simplify[{
		resolve[musqU1P1/.PrintDebyeMass["All"]]-g1^2/ge^2 debyeQED,
		resolve[musqU1P2/.PrintDebyeMass["All"]]-g2^2/ge^2 debyeQED,
		resolve[muU1Mix1/.PrintDebyeMass["All"]]-g1 g2/ge^2 debyeQED}/.rotation,
		Assumptions->{g1>0,g2>0,g13d>0,g23d>0}],
	{0,0,0},
	TestID->"two-U1-debye-masses-are-rotated-scalar-QED"
]];

AppendTo[testList,
VerificationTest[
	Simplify[resolve[Normal[PrintTensorDRalgo[4]]]-Array[gList[[#1]]gList[[#2]]gList[[#3]]gList[[#4]]&,{2,2,2,2}]/ge^4 quarticQED/.rotation],
	Array[0&,{2,2,2,2}],
	TestID->"two-U1-temporal-scalar-quartics-are-rotated-scalar-QED"
]];


(* ::Section:: *)
(*QCD*)


(*Without scalars, the model gets a placeholder scalar that does not couple*)
Group={"SU3"};
RepAdjoint={{1,1}};
CouplingName={gs};
RepScalar={};
RepFermion={{{{1,0}},"L"},{{{1,0}},"R"}};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PerformDRhard[];
\[Lambda]A=PrintTensorDRalgo[4][[1,1,1,1]]/6;
(*For SU(3), (A0^a A0^a)^2 is the only quartic invariant, with dA=8 and CA=3*)
betaQCD=2(8+2)\[Lambda]A(\[Lambda]A-3gs3d^2)/(16\[Pi]^2);
(*Outside of the tests, since the placeholder scalar leads to messages*)
runningQCD=debyeBeta[musqSU3]-resolve[betaQCD];
loQCD=musqSU3/.PrintDebyeMass["LO"];
nloQCD=resolve[musqSU3/.PrintDebyeMass["NLO"]];
betasQCD=BetaFunctions4D[];

AppendTo[testList,
VerificationTest[
	Simplify[runningQCD],
	0,
	TestID->"QCD-debye-mass-runs-as-SU(N)-form-of-A3"
]];

AppendTo[testList,
VerificationTest[
	order[running[loQCD,betasQCD]+explicitMu[nloQCD/.gs3d->gs Sqrt[T]],4,{gs->x gs}],
	0,
	TestID->"QCD-debye-mass-is-rg-invariant-at-O(g^4)"
]];


(* ::Section:: *)
(*Report*)


report=TestReport[testList];
Print[report];
Print[Column[{#["TestID"],#["Outcome"]}&/@Values[report["TestResults"]]]];
Exit[If[report["AllTestsSucceeded"],0,1]]
