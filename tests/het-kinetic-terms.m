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
(*Scalar kinetic terms in the HET*)


(*
	Issue #42: PrintScalarKineticHET[] for examples/tripletHET.m, with the field normalization from the
	scalar-scalar, scalar-vector, and vector-vector bubbles of ScalarSelfEnergyHET[]. The kinetic term is
	(1+2 \[Delta]Z) (\[PartialD]R)^2/2 with \[Delta]Z=PrintScalarKineticHET[], in Landau gauge.

	Checked here:
	- tripletHET.m with the SU(2) vectors integrated out: for the Higgs, 2\[Delta]Z=-11 g/(16\[Pi] \[CurlyPhi]) is the
	  term Z_NLO of 2512.05565 at g'=0; for the Goldstone modes and the triplet, the scalar-vector bubble
	  with a massless scalar gives 2\[Delta]Z=-g/(\[Pi] \[CurlyPhi]) and -8g/(3\[Pi] \[CurlyPhi]);
	- tripletHET.m with the triplet integrated out: the scalar-scalar bubble with two heavy lines;
	- tripletHET.m in the soft theory with the SU(2) vectors and the temporal scalars integrated out,
	  also with the vev {0,0,0,\[CurlyPhi]} of the 4d doublet only, as in issue #42, which DefineVEVS pads with zeros;
	- the scalar-scalar bubble with one light and one heavy line, -1/(12\[Pi] M^3) in the p^2 term of
	  1/(4\[Pi] p) ArcTan[p/M], for three real scalars with the cubic coupling c h \[Chi] s.
*)


testList={};


(* ::Section:: *)
(*tripletHET.m*)


Group={"SU3","SU2","U1"};
RepAdjoint={{1,1},{2},0};
HiggsDoublet1={{{0,0},{1},1/2},"C"};
HiggsTriplet1={{{0,0},{2},0},"R"};
RepScalar={HiggsDoublet1,HiggsTriplet1};
CouplingName={g3,gw,g1};

Rep1={{{1,0},{1},1/6},"L"};
Rep2={{{1,0},{0},2/3},"R"};
Rep3={{{1,0},{0},-1/3},"R"};
Rep4={{{0,0},{1},-1/2},"L"};
Rep5={{{0,0},{0},-1},"R"};
RepFermion1Gen={Rep1,Rep2,Rep3,Rep4,Rep5};
RepFermion3Gen={RepFermion1Gen,RepFermion1Gen,RepFermion1Gen}//Flatten[#,1]&;

{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion3Gen,RepScalar];

MassTerm1=CreateInvariant[Group,RepScalar,{{1,1},{True,False}}][[1]]//Simplify//FullSimplify;
MassTerm2=CreateInvariant[Group,RepScalar,{{2,2},{True,True}}][[1]]//Simplify//FullSimplify;
\[Mu]ij=GradMass[mH*MassTerm1+m\[CapitalSigma]/2*MassTerm2]//Simplify//SparseArray;
\[Lambda]4=GradQuartic[\[Lambda]*MassTerm1^2+b4/4*MassTerm2^2+a2/2*MassTerm1*MassTerm2];
YukawaDoublet1=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,{{1,1,2},{False,False,True}}]//Simplify;
Ysff=-GradYukawa[yt*YukawaDoublet1[[1]]];
YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt>0}]];

ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PerformDRhard[];
PerformDRsoft[{}];

positive={gw3dUS>0,gw3d>0,\[CurlyPhi]>0};


(* ::Subsection:: *)
(*Ultrasoft theory*)


UseUltraSoftTheory[];
DefineVEVS[{0,0,0,\[CurlyPhi],0,0,0}//SparseArray];

(*The SU(2) vectors live on 9, 10, 11*)
PrepareHET[{},{9,10,11}];
CalculatePotentialHET[];
kinetic=Normal[PrintScalarKineticHET[]];

AppendTo[testList,
VerificationTest[
	Simplify[2Diagonal[kinetic],positive],
	gw3dUS/(\[Pi] \[CurlyPhi]){-1,-1,-1,-11/16,-8/3,-8/3,-8/3},
	TestID->"heavy-SU2-vectors-ultrasoft"
]];

(*The triplet lives on 5, 6, 7*)
PrepareHET[{5,6,7},{}];
CalculatePotentialHET[];
kinetic=Normal[PrintScalarKineticHET[]];
mTriplet=m\[CapitalSigma]3dUS+a23dUS \[CurlyPhi]^2/2;

(*Bubble of two heavy lines with the same mass M^2, whose p^2 term is -1/(96\[Pi] M^3): 2\[Delta]Z=1/2 Sum[\[Lambda]3^2]/(96\[Pi] M^3)*)
AppendTo[testList,
VerificationTest[
	Simplify[2kinetic-DiagonalMatrix[{0,0,0,1/2 3(a23dUS \[CurlyPhi])^2/(96\[Pi] mTriplet^(3/2))}]],
	ConstantArray[0,{4,4}],
	TestID->"heavy-triplet-ultrasoft"
]];


(* ::Subsection:: *)
(*Soft theory*)


UseSoftTheory[];
(*The temporal scalars follow the 7 scalars: SU(3) on 8,...,15, SU(2) on 16, 17, 18, and U(1) on 19*)
heavySoft={Range[8,19],{9,10,11}};

DefineVEVS[{0,0,0,\[CurlyPhi],0*Range[3],0*Range[8],0*Range[3],0}//Flatten//SparseArray];
PrepareHET@@heavySoft;
CalculatePotentialHET[];
kinetic=Normal[PrintScalarKineticHET[]];

(*Without the couplings of the doublet to the SU(2) and U(1) temporal scalars, only the SU(3) bubble remains*)
mGluon=musqSU3+lamVL[4]\[CurlyPhi]^2/2;
AppendTo[testList,
VerificationTest[
	Simplify[(2Diagonal[kinetic]/.{lamVL[5]->0,lamVL[6]->0,lamVL[7]->0})
		-gw3d/(\[Pi] \[CurlyPhi]){-1,-1,-1,-11/16,-8/3,-8/3,-8/3}
		-{0,0,0,1/2 8(lamVL[4]\[CurlyPhi])^2/(96\[Pi] mGluon^(3/2)),0,0,0},Join[positive,{mGluon>0}]],
	ConstantArray[0,7],
	TestID->"heavy-SU2-vectors-and-temporal-scalars-soft"
]];

(*The vev of the 4d doublet only, which DefineVEVS pads with zeros*)
Quiet[DefineVEVS[{0,0,0,\[CurlyPhi]}//SparseArray]];
PrepareHET@@heavySoft;
CalculatePotentialHET[];

AppendTo[testList,
VerificationTest[
	Simplify[Normal[PrintScalarKineticHET[]]-kinetic],
	ConstantArray[0,{7,7}],
	TestID->"vev-of-the-4d-doublet-only"
]];


(* ::Section:: *)
(*Scalar-scalar bubble with one light and one heavy line*)


(*
	The light fields h and \[Chi] and the heavy field s; with the background h=\[CurlyPhi], the cubic coupling c h \[Chi] s
	gives the bubbles (h,s) for \[Chi] and (\[Chi],s) for h, each with one massless line. The p^2 term of the bubble
	is -c^2/(12\[Pi] M^3), and \[Delta]Z=c^2/(24\[Pi] M^3) with M^2=ms.
*)
\[Lambda]3Light=SparseArray[Thread[Permutations[{1,2,3}]->c],{3,3,3}];
DefineNewTensorsUS[DiagonalMatrix[{mh,m\[Chi],ms}]//SparseArray,SparseArray[{},{3,3,3,3}],\[Lambda]3Light,SparseArray[{},{1,3,3}],SparseArray[{},{1,1,1}]];
Quiet[DefineVEVS[{\[CurlyPhi],0,0}//SparseArray]];
PrepareHET[{3},{}];
CalculatePotentialHET[];

AppendTo[testList,
VerificationTest[
	Simplify[Normal[PrintScalarKineticHET[]],ms>0],
	c^2/(24\[Pi] ms^(3/2)) IdentityMatrix[2],
	TestID->"bubble-with-one-light-and-one-heavy-line"
]];


(* ::Section:: *)
(*Report*)


report=TestReport[testList];
Print[report];
Print[Column[{#["TestID"],#["Outcome"]}&/@Values[report["TestResults"]]]];
Exit[If[report["AllTestsSucceeded"],0,1]]
