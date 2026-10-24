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
(*Gauge invariance of the couplings*)


(*
	DRalgo converts right-handed fermions to left-handed ones with the conjugate generators, and Ysff
	multiplies the conjugates of these fields. In CreateInvariantYukawa, left-handed fermions enter
	conjugated and right-handed fermions not, as in sm.m with {False,False,True} for (\[Phi],QL,tR).
	Then gvss[c].Ysff+I(gvff[c].Ysff+(gvff[c].Ysff)^T)=0 for every generator c, where gvff[c] acts on the
	first fermion index and the transpose is in the fermion indices.

	For the opposite flags, Ysff contains the conjugate invariants, which are not invariant under gvff.
	The fermion loops in the one-loop couplings of the temporal scalars to the scalars then break the
	gauge symmetry, and PrintDebyeMass fails, since the Debye masses of one group run differently.
	CreateInvariantYukawa sets the U(1) charges to zero, so the conjugation of the scalar has to fit them:
	in the Standard Model, the up-type quarks couple to the conjugate Higgs and the down-type quarks to the Higgs.

	ImportModelDRalgo checks the invariance of all couplings and warns with ImportModelDRalgo::noninvariant.
	With symbolic U(1) charges, it gives the condition on the charges with ImportModelDRalgo::invariantif.

	Checked here:
	- the Yukawa couplings of sm.m, SU5.m, and LRSymmetric.m are invariant, the conjugate ones of LRSymmetric.m not;
	- in the Standard Model with top and bottom quarks, BetaFunctions4D gives the one-loop \[Beta]-functions
	  of the Yukawa couplings;
	- smCKM.m imports without a warning, but with the flags of the up-type quarks for the down-type quarks not;
	- LRSymmetric.m with the conjugate flags gives the warning;
	- the one-loop couplings of the temporal scalars to the scalars of LRSymmetric.m are gauge invariant,
	  and PrintDebyeMass gives one Debye mass per group without messages;
	- xsm.m, with symbolic hypercharges, is invariant only for Yu=Yq+Y\[Phi];
	- a Dirac fermion with a mass and a Yukawa coupling to a real scalar is invariant if both chiralities
	  have the same U(1) charge, and, with symbolic charges qL and qR, only for qL=qR;
	- the cubic coupling d_abc of an SU(3) octet is invariant for the type "R", but not for the type "A",
	  whose components and gauge couplings are in different bases.
*)


testList={};

private[name_String]:=ToExpression["DRalgo`DRalgo`Private`"<>name];

(*The model part of an example, up to and including the call of ImportModelDRalgo*)
exampleModel[file_]:=Module[{text,start,stop},
	text=Import[FileNameJoin[{"..","examples",file}],"Text"];
	text=StringReplace[text,{"SetDirectory[NotebookDirectory[]];"->"","<<../Kernel/DRalgo.m"->"","<<DRalgo`DRalgo`"->""}];
	start=StringPosition[text,"ImportModelDRalgo["][[1,1]];
	stop=SelectFirst[StringPosition[text,"];"],#[[1]]>start&][[2]];
	StringTake[text,stop]
];

(*Variation of the Yukawa couplings y_iIJ under the generator c*)
yukawaVariation[y_,c_]:=Module[{yN=Normal[y],fermionY},
	fermionY=Transpose[Normal[gvff][[c]] . Transpose[yN,{2,1,3}],{2,1,3}];
	Normal[gvss][[c]] . yN+I (fermionY+Transpose[fermionY,{1,3,2}])
];
yukawaResidual[y_]:=Union@Flatten@Simplify@Table[yukawaVariation[y,c],{c,Length[gvss]}];

(*Variation of the couplings G_abij of two temporal scalars to two scalars under the generator c*)
temporalVariation[g_,c_]:=Module[{gN=Normal[g],f=Normal[gvvv][[c]],s=Normal[gvss][[c]]},
	f . gN+Transpose[f . Transpose[gN,{2,1,3,4}],{2,1,3,4}]
	+Transpose[s . Transpose[gN,{2,3,1,4}],{3,1,2,4}]+Transpose[s . Transpose[gN,{2,3,4,1}],{4,1,2,3}]
];


(* ::Section:: *)
(*Standard Model*)


ToExpression[exampleModel["sm.m"]];

AppendTo[testList,
VerificationTest[
	yukawaResidual[Ysff],
	{0},
	TestID->"sm-yukawa-couplings-are-invariant"
]];

(*One generation of quarks with top and bottom Yukawa couplings*)
Module[{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,group,scalars,fermions,mass,top,bottom},
	group={"SU3","SU2","U1"};
	scalars={{{{0,0},{1},1/2},"C"}};
	fermions={{{{1,0},{1},1/6},"L"},{{{1,0},{0},2/3},"R"},{{{1,0},{0},-1/3},"R"}};
	{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=
		AllocateTensors[group,{{1,1},{2},0},{gs,gw,gY},fermions,scalars];
	mass=CreateInvariant[group,scalars,{{1,1},{True,False}}][[1]]//Simplify;
	\[Mu]ij=GradMass[msq*mass]//SparseArray;
	\[Lambda]4=GradQuartic[\[Lambda]*mass^2];
	top=CreateInvariantYukawa[group,scalars,fermions,{{1,1,2},{False,False,True}}][[1]]//Simplify; (*\[Phi]^+ QL^+ tR*)
	bottom=CreateInvariantYukawa[group,scalars,fermions,{{1,1,3},{True,False,True}}][[1]]//Simplify; (*\[Phi] QL^+ bR*)
	Ysff=-GradYukawa[yt*top+yb*bottom];
	YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt>0,yb>0}]];
	ImportModelDRalgo[group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
];
\[Beta]4D=BetaFunctions4D[];

AppendTo[testList,
VerificationTest[
	Simplify[{yt,yb}/.\[Beta]4D]-{
		yt/(16\[Pi]^2) (9/2 yt^2+3/2 yb^2-8gs^2-9/4 gw^2-17/12 gY^2),
		yb/(16\[Pi]^2) (9/2 yb^2+3/2 yt^2-8gs^2-9/4 gw^2-5/12 gY^2)}//Simplify,
	{0,0},
	TestID->"sm-top-and-bottom-yukawa-beta-functions"
]];


(* ::Section:: *)
(*SU(5)*)


ToExpression[exampleModel["SU5.m"]];

AppendTo[testList,
VerificationTest[
	yukawaResidual[Ysff],
	{0},
	TestID->"SU5-yukawa-couplings-are-invariant"
]];


(* ::Section:: *)
(*Standard Model with the CKM matrix*)


AppendTo[testList,
VerificationTest[
	ToExpression[exampleModel["smCKM.m"]],
	Null,
	TestID->"smCKM-imports-without-warning"
]];

AppendTo[testList,
VerificationTest[
	ToExpression[StringReplace[exampleModel["smCKM.m"],
		"{{1,1+(i-1)*5,3+(j-1)*5},{True,False,True}}"->"{{1,1+(i-1)*5,3+(j-1)*5},{False,False,True}}"]],
	Null,
	{ImportModelDRalgo::noninvariant},
	TestID->"smCKM-with-the-flags-of-the-up-type-quarks-warns"
]];


(* ::Section:: *)
(*Left-right symmetric model*)


AppendTo[testList,
VerificationTest[
	ToExpression[StringReplace[exampleModel["LRSymmetric.m"],{
		"InputInv={{1,1,2},{False,False,True}};"->"InputInv={{1,1,2},{True,True,False}};",
		"InputInv={{1,1,2},{True,False,True}};"->"InputInv={{1,1,2},{False,True,False}};"}]],
	Null,
	{ImportModelDRalgo::noninvariant},
	TestID->"LRSymmetric-with-the-conjugate-flags-warns"
]];

ToExpression[exampleModel["LRSymmetric.m"]];

AppendTo[testList,
VerificationTest[
	yukawaResidual[Ysff],
	{0},
	TestID->"LRSymmetric-yukawa-couplings-are-invariant"
]];

AppendTo[testList,
VerificationTest[
	yukawaResidual[YsffC]==={0},
	False,
	TestID->"LRSymmetric-conjugate-yukawa-couplings-are-not-invariant"
]];

PerformDRhard[];

AppendTo[testList,
VerificationTest[
	Union@Flatten@Simplify@Table[temporalVariation[private["GvvssL"],c],{c,Length[gvss]}],
	{0},
	TestID->"LRSymmetric-temporal-scalar-couplings-are-gauge-invariant"
]];

AppendTo[testList,
VerificationTest[
	AllTrue[PrintDebyeMass["NLO"][[All,2]],FreeQ[#,_Rule|_List]&],
	True,
	TestID->"LRSymmetric-debye-masses"
]];


(* ::Section:: *)
(*Symbolic U(1) charges*)


AppendTo[testList,
VerificationTest[
	ToExpression[exampleModel["xsm.m"]],
	Null,
	{ImportModelDRalgo::invariantif},
	TestID->"xsm-warns-with-the-condition-on-the-hypercharges"
]];

AppendTo[testList,
VerificationTest[
	Union[Cases[private["InvarianceConditions"][],(_->c_)/;!TrueQ[c]:>Simplify[c/.Yu->Yq+Y\[Phi]]]],
	{True},
	TestID->"xsm-is-invariant-for-the-standard-model-hypercharges"
]];

(*A Dirac fermion with the charges qL and qR, a mass, and a Yukawa coupling to a real scalar*)
diracFermion[qL_,qR_]:=Module[{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,group,scalars,fermions,\[Phi]},
	group={"U1"};
	scalars={{{0},"R"}};
	fermions={{{qL},"L"},{{qR},"R"}};
	{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[group,{0},{g1},fermions,scalars];
	\[Phi]=CreateInvariant[group,scalars,{{1},{True}}][[1]];
	\[Mu]ij=GradMass[msq/2 \[Phi]^2]//SparseArray;
	\[Lambda]4=GradQuartic[\[Lambda]/24 \[Phi]^4];
	\[Mu]IJ=m\[Psi]*GradMassFermion[CreateInvariantFermion[group,fermions,{{1,2},{False,True}}][[1]]]; (*\[Psi]L^+ \[Psi]R*)
	\[Mu]IJC=m\[Psi]*GradMassFermion[CreateInvariantFermion[group,fermions,{{2,1},{False,True}}][[1]]]; (*\[Psi]R^+ \[Psi]L*)
	Ysff=y*GradYukawa[CreateInvariantYukawa[group,scalars,fermions,{{1,1,2},{True,False,True}}][[1]]]; (*\[Phi] \[Psi]L^+ \[Psi]R*)
	YsffC=y*GradYukawa[CreateInvariantYukawa[group,scalars,fermions,{{1,2,1},{True,False,True}}][[1]]]; (*\[Phi] \[Psi]R^+ \[Psi]L*)
	ImportModelDRalgo[group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False]
];

AppendTo[testList,
VerificationTest[
	diracFermion[1,1],
	Null,
	TestID->"dirac-fermion-is-invariant"
]];

AppendTo[testList,
VerificationTest[
	diracFermion[1,2],
	Null,
	{ImportModelDRalgo::noninvariant},
	TestID->"dirac-fermion-with-different-charges-warns"
]];

AppendTo[testList,
VerificationTest[
	diracFermion[qL,qR],
	Null,
	{ImportModelDRalgo::invariantif},
	TestID->"dirac-fermion-with-symbolic-charges-warns"
]];

AppendTo[testList,
VerificationTest[
	With[{conditions=private["InvarianceConditions"][]},
		{Cases[conditions,(name_->c_)/;!TrueQ[c]:>name],Union[Cases[conditions,(_->c_)/;!TrueQ[c]:>Simplify[c/.qR->qL]]]}],
	{{"muIJ","muIJC","Ysff","YsffC"},{True}},
	TestID->"dirac-fermion-is-invariant-only-for-equal-charges"
]];


(* ::Section:: *)
(*Adjoint scalars*)


(*An SU(3) octet of the type "R" or "A" with a mass, a cubic, and a quartic coupling*)
octet[type_]:=Module[{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,group,scalars,square},
	group={"SU3"};
	scalars={{{{1,1}},type}};
	{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[group,{{1,1}},{gs},{},scalars];
	square=CreateInvariant[group,scalars,{{1,1},{True,True}}][[1]];
	\[Mu]ij=GradMass[msq/2 square]//SparseArray;
	\[Lambda]3=GradCubic[c3/6 CreateInvariant[group,scalars,{{1,1,1},{True,True,True}}][[1]]];
	\[Lambda]4=GradQuartic[\[Lambda] square^2];
	ImportModelDRalgo[group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False]
];

AppendTo[testList,
VerificationTest[
	octet["R"],
	Null,
	TestID->"octet-of-type-R-is-invariant"
]];

AppendTo[testList,
VerificationTest[
	octet["A"],
	Null,
	{ImportModelDRalgo::noninvariant},
	TestID->"octet-of-type-A-warns"
]];

AppendTo[testList,
VerificationTest[
	private["InvarianceConditions"][],
	{"lam3"->False,"lam4"->True,"muij"->True},
	TestID->"octet-of-type-A-has-a-noninvariant-cubic"
]];


(* ::Section:: *)
(*Report*)


report=TestReport[testList];
Print[report];
Print[Column[{#["TestID"],#["Outcome"]}&/@Values[report["TestResults"]]]];
Exit[If[report["AllTestsSucceeded"],0,1]]
