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
(*BetaFunctions4D for Lagrangians that are not closed under renormalisation*)


(*
	Issue #54: a 3HDM with the \[CapitalSigma](36)-symmetric potential of examples/3hdm.m in which
	the top quark couples only to \[Phi]3. BetaFunctions4D gave \[Beta](\[Lambda]11) without a term yt3^4.

	The Yukawa coupling breaks \[CapitalSigma](36), and the fermion loop generates (\[Phi]3^+ \[Phi]3)^2 alone,
	which the two couplings \[Lambda]11, \[Lambda]12 cannot express. The \[Beta]-function tensor contains the
	yt3^4 term; BetaFunctions4D projects it onto the tree-level couplings, and its
	components are inconsistent: some contain yt3^4, others do not. Up to DRalgo 1.5.0 the
	printed result was one of them, without any warning.

	BetaFunctions4D now checks that the \[Beta]-function tensors are reproduced by the printed
	\[Beta]-functions, and warns if not. examples/3hdm.m had the same problem: its top Yukawa
	couplings broke \[CapitalSigma](36), and its \[Lambda]H3 term was not invariant even under the cyclic
	permutations of the doublets.

	Checked here:
	- the issue model: BetaFunctions4D warns about the scalar quartics, and the fermion
	  loop is present in the (\[Phi]3^+ \[Phi]3)^2 component of the \[Beta]-function tensor;
	- the SM: no warning;
	- examples/3hdm.m: the potential is \[CapitalSigma](36)-invariant, and no warning.
*)


testList={};
private[name_String]:=ToExpression["DRalgo`DRalgo`Private`"<>name];

(*The lines that an expression prints*)
SetAttributes[printedBy,HoldAll];
printedBy[expr_]:=Reap[Block[{Print=Sow[StringJoin[ToString/@{##}]]&},expr]][[2]]/.{{s___}}:>{s};
warningQ[lines_,sector_String]:=AnyTrue[lines,StringContainsQ[#,"\[Beta]-functions not defined at tree-level"]&&StringContainsQ[#,sector]&];

(*BetaFunctions4D needs only the counterterms, so the models are not reduced*)
import[]:=ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];


(* ::Section:: *)
(*Issue model: \[CapitalSigma](36) potential, top quark coupled to \[Phi]3 only*)


Group={"SU3","SU2","U1"};
RepAdjoint={{1,1},{2},0};
RepScalar=Table[{{{0,0},{1},1/2},"C"},{3}];
CouplingName={g3,g2,g1};
RepFermion={{{{1,0},{1},1/6},"L"},{{{1,0},{0},2/3},"R"},{{{1,0},{0},-1/3},"R"},{{{0,0},{1},-1/2},"L"},{{{0,0},{0},-1},"R"}};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];

(*p[i,j] = \[Phi]i^+ \[Phi]j*)
p[i_,j_]:=CreateInvariant[Group,RepScalar,{{j,i},{True,False}}][[1]]//Simplify;
\[Mu]ij=GradMass[-\[Mu]sq (p[1,1]+p[2,2]+p[3,3])]//Simplify//SparseArray;
\[Lambda]4=GradQuartic[
	+\[Lambda]11 (p[1,1]+p[2,2]+p[3,3])^2
	+\[Lambda]12 (2(p[1,2]p[2,1]+p[2,3]p[3,2]+p[3,1]p[1,3])
		-p[1,2]p[2,3]-p[2,3]p[3,1]-p[3,1]p[1,2]-p[2,1]p[3,2]-p[3,2]p[1,3]-p[1,3]p[2,1])];
Ysff=-GradYukawa[yt3 CreateInvariantYukawa[Group,RepScalar,RepFermion,{{3,1,2},{False,False,True}}][[1]]];
YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt3>0}]];
import[];

lines=printedBy[BetaFunctions4D[]];

AppendTo[testList,
VerificationTest[
	warningQ[lines,"scalar quartics"],
	True,
	TestID->"issue-model-warns-about-quartics"
]];

(*
	The real components of \[Phi]i have the scalar indices 4i-3,...,4i. As in the SM, the
	fermion loop gives \[Beta](\[Lambda]) = -6 yt^4/(16\[Pi]^2) for \[Lambda](\[Phi]^+\[Phi])^2, whose component (a,a,a,a)
	is 6\[Lambda]; here for \[Phi]3 only.
*)
\[Beta]Quartic=Normal[private["\[Beta]\[Lambda]ijkl"]];
AppendTo[testList,
VerificationTest[
	Simplify[Coefficient[Expand[{\[Beta]Quartic[[9,9,9,9]],\[Beta]Quartic[[1,1,1,1]],\[Beta]Quartic[[5,5,5,5]]}],yt3^4]],
	{-36/(16\[Pi]^2),0,0},
	TestID->"issue-model-beta-tensor-contains-top-loop-for-phi3"
]];


(* ::Section:: *)
(*SM: closed under renormalisation*)


RepScalar={{{{0,0},{1},1/2},"C"}};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
\[Phi]sq=CreateInvariant[Group,RepScalar,{{1,1},{True,False}}][[1]]//Simplify;
\[Mu]ij=GradMass[m2 \[Phi]sq]//Simplify//SparseArray;
\[Lambda]4=GradQuartic[\[Lambda]1H \[Phi]sq^2];
Ysff=-GradYukawa[yt CreateInvariantYukawa[Group,RepScalar,RepFermion,{{1,1,2},{False,False,True}}][[1]]];
YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt>0}]];
import[];

AppendTo[testList,
VerificationTest[
	printedBy[BetaFunctions4D[]],
	{},
	TestID->"SM-beta-functions-without-warning"
]];


(* ::Section:: *)
(*examples/3hdm.m*)


(*Evaluates the model part of an example file, without the lines that load DRalgo*)
loadExampleModel[file_String]:=Module[{text,lines,tmp},
	text=Import[file,"Text"];
	text=StringTake[text,First[First[StringPosition[text,"ImportModelDRalgo["]]]-1];
	lines=Select[StringSplit[text,"\n"],
		!StringContainsQ[#,"SetDirectory"|"NotebookDirectory"|"Kernel/DRalgo.m"|"$LoadGroupMath"]&];
	tmp=FileNameJoin[{$TemporaryDirectory,FileBaseName[file]<>"-model-"<>ToString[$ProcessID]<>".m"}];
	Export[tmp,StringRiffle[lines,"\n"],"Text"];
	Get[tmp];
	DeleteFile[tmp];
];
loadExampleModel["../examples/3hdm.m"];

(*
	Generators of \[CapitalSigma](36) on the doublets (\[Phi]1,\[Phi]2,\[Phi]3), see 1410.6139: the cyclic
	permutation, the phase rotation and the matrix \[Omega]^((i-1)(j-1))/(I Sqrt[3]).
*)
\[Omega]3=Exp[2\[Pi] I/3];
sigma36={
	{{0,0,1},{1,0,0},{0,1,0}},
	DiagonalMatrix[{1,\[Omega]3,\[Omega]3^2}],
	Table[\[Omega]3^((i-1)(j-1)),{i,3},{j,3}]/(I Sqrt[3])
};
potential=VMass+VQuartic;

(*The real components \[Phi]a[1,k] and \[Psi]a[1,k] of the doublets (in DRalgo's private context)*)
components=Cases[Variables[potential],_Symbol[__]];
component[name_String,a_,k_]:=First[Select[components,SymbolName[Head[#]]==name<>ToString[a]&&Last[#]==k&]];

(*A unitary U acting on the complex components (\[Phi]a + I \[Psi]a)/Sqrt[2], written for the real components*)
doubletRotation[U_]:=Module[{re=ComplexExpand[Re[U]],im=ComplexExpand[Im[U]],phi,psi},
	phi[a_,k_]:=component["\[Phi]",a,k];
	psi[a_,k_]:=component["\[Psi]",a,k];
	Flatten@Table[{
		phi[a,k]->Sum[re[[a,b]]phi[b,k]-im[[a,b]]psi[b,k],{b,3}],
		psi[a,k]->Sum[im[[a,b]]phi[b,k]+re[[a,b]]psi[b,k],{b,3}]},
	{a,3},{k,2}]
];
invariantUnderQ[V_,U_]:=Simplify[Expand[(V/.doubletRotation[U])-V]]===0;

AppendTo[testList,
VerificationTest[
	Length[components],
	12,
	TestID->"3hdm-potential-depends-on-12-real-components"
]];

(*Control: the sign flip of one doublet is no symmetry of the \[Lambda]H3 term*)
AppendTo[testList,
VerificationTest[
	invariantUnderQ[potential,DiagonalMatrix[{1,-1,1}]],
	False,
	TestID->"3hdm-potential-not-invariant-under-sign-flip-of-phi2"
]];

AppendTo[testList,
VerificationTest[
	AllTrue[Join[sigma36,Conjugate[sigma36]],invariantUnderQ[potential,#]&],
	True,
	TestID->"3hdm-potential-invariant-under-sigma36"
]];

import[];
AppendTo[testList,
VerificationTest[
	printedBy[BetaFunctions4D[]],
	{},
	TestID->"3hdm-beta-functions-without-warning"
]];


(* ::Section:: *)
(*Report*)


report=TestReport[testList];
Print[report];
Print[Column[{#["TestID"],#["Outcome"]}&/@Values[report["TestResults"]]]];
Exit[If[report["AllTestsSucceeded"],0,1]]
