(* ::Package:: *)

If[$InputFileName=="",
	SetDirectory[NotebookDirectory[]],
	SetDirectory[DirectoryName[$InputFileName]]
];
DRalgo`DRalgo`$LoadGroupMath=True;

Check[
    Get["../Kernel/DRalgo.m"],
    Message[Get::noopen, "DRalgo` at "<>ToString[$UserBaseDirectory]<>"/Applications"];
    Abort[];
]


(* ::Chapter:: *)
(*Operators with six fields*)


(*
	The tensors of the operators with six fields, phi^6, phi^4 A0^2, phi^2 A0^4, and A0^6, are sums of contractions
	with up to nv^6 components. SymmetricTensorSum (HigherDimensionalOperators.m) computes each contraction once up to
	the relabeling of identical tensors and sums the projections onto the tensors symmetric under the permutations of
	the scalars and of the temporal scalars, which the operators project onto.
*)


testList={};
private[s_]:=ToExpression["DRalgo`DRalgo`Private`"<>s];


(* ::Section:: *)
(*Sums of contractions*)


(*
	Contractions of a tensor g with itself and with a tensor h, of which some are equal up to relabeling the copies of g,
	with transposes and symbolic coefficients: the sum equals the projection of the sum of the evaluated terms.
*)
SeedRandom[6];
g=SparseArray[RandomInteger[{-3,3},{3,3,3}]];
h=SparseArray[RandomInteger[{-3,3},{3,3}]];
heldSum=Hold[
	+a Contract[g,g,g,{{2,4},{5,7},{3,8}}]
	+b Contract[g,g,g,{{2,7},{3,5},{6,8}}]
	-2 Transpose[Contract[g,g,g,{{3,4},{6,7},{2,9}}],{3,1,2}]
	+(a+b)/3 Transpose[Contract[g,g,g,{{2,5},{4,8},{3,7}}],{2,3,1}]
	+c Contract[g,h,{{3,4}}]
	+Transpose[Contract[h,g,{{1,4}}],{2,1,3}]
];
symmetrize[t_,groups_]:=Fold[Symmetrize[#1,Symmetric[#2]]&,t,Select[groups,Length[#]>1&]];

Do[
	AppendTo[testList,
	VerificationTest[
		Normal[heldSum/.Hold[s_]:>private["SymmetricTensorSum"][groups,s]]-Normal[symmetrize[ReleaseHold[heldSum],groups]]//Simplify//Flatten//Union,
		{0},
		TestID->"symmetric-sum-"<>StringRiffle[ToString/@Length/@groups,"-"]
	]],
{groups,{{{1,2,3}},{{1,2},{3}},{{1},{2},{3}}}}];


(* ::Section:: *)
(*A0^6 of SU(3) in d dimensions*)


(*
	The coefficients of A0^6 of the pure SU(3) gauge theory, as listed in 2605.15176 for d spatial dimensions.
	The 1194 terms of V6B contain 804 contractions of six gvvv, which are 4 distinct ones up to relabeling.
*)
Group={"SU3"};CouplingName={gs};RepAdjoint={{1,1}};RepScalar={};RepFermion={};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PerformDRhard[];

Tadj=-I gvvv/gs;
X2=Contract[Tadj,Tadj,{{2,6},{3,5}}];
X4=Contract[Tadj,Tadj,Tadj,Tadj,{{2,12},{3,5},{6,8},{9,11}}];
X6=Contract[Tadj,Tadj,Tadj,Tadj,Tadj,Tadj,{{2,18},{3,5},{6,8},{9,11},{12,14},{15,17}}];
OT612=SymmetrizeTensor[\[Alpha][A^6,1]X6+\[Alpha][A^6,2]TensorProduct[X4,X2],{{1,2,3,4,5,6}}];

AppendTo[testList,
VerificationTest[
	Simplify[{\[Alpha][A^6,1],\[Alpha][A^6,2]}/.Dimension6Matching[{OT612},{12},{\[Alpha][A^6,1],\[Alpha][A^6,2]},d][[1]]],
	{8(-5+d)(-3+d)(-1+d)^2 gs^6 T^2 Zb[3,0],0},
	SameTest->(Simplify[#1-#2]==={0,0}&),
	TestID->"A0^6-of-SU3-in-d-dimensions"
]];


(* ::Section:: *)
(*Report*)


report=TestReport[testList];
Print[report];
Print[Column[{#["TestID"],#["Outcome"]}&/@Values[report["TestResults"]]]];
Exit[If[report["AllTestsSucceeded"],0,1]]
