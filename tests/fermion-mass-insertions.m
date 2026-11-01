(* ::Package:: *)

If[$InputFileName=="",
	SetDirectory[NotebookDirectory[]],
	SetDirectory[DirectoryName[$InputFileName]]
];
DRalgo`DRalgo`$LoadGroupMath=True;
DRalgo`DRalgo`$GroupMathMultipleModels=True;

Check[
    Get["../Kernel/DRalgo.m"],
    Message[Get::noopen, "DRalgo` at "<>ToString[$UserBaseDirectory]<>"/Applications"];
    Abort[];
]


(* ::Chapter:: *)
(*Dimension-5 operators with fermion-mass insertions*)


(*
	A fermion loop with an odd number of scalar legs needs an odd number of insertions of the fermion mass.
	The fermion loops of the dimension-5 operators are therefore proportional to the fermion masses, and
	ODIM5Fermion (HigherDimensionalOperators.m) obtains them from the fermion loops of the dimension-6
	operators, with the masses as the Yukawa couplings of a constant scalar.
*)


testList={};


(* ::Section:: *)
(*U(1) with a charged Dirac fermion*)


(*
	A Dirac fermion of charge 1 and a neutral real scalar, with V = lam phi^4, the fermion mass m\[Psi], and the Yukawa
	couplings y + I y5 of phi psi_R^+ psi_L and y - I y5 of its conjugate, i.e. a scalar coupling y and a pseudoscalar
	coupling y5. At d = 3, the coefficients agree with the matching of Matchotter (2604.21972), both off-shell and
	on-shell.
*)
Group={"U1"};CouplingName={g1};RepAdjoint={0};RepScalar={{{0},"R"}};RepFermion={{{1},"L"},{{1},"R"}};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJC,\[Mu]IJ,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
\[Lambda]4=GradQuartic[\[Lambda] CreateInvariant[Group,RepScalar,{{1,1},{True,True}}][[1]]^2];
\[Mu]IJ=m\[Psi] GradMassFermion[CreateInvariantFermion[Group,RepFermion,{{2,1},{False,True}}][[1]]];
\[Mu]IJC=m\[Psi] GradMassFermion[CreateInvariantFermion[Group,RepFermion,{{1,2},{False,True}}][[1]]];
Ysff=(y+I y5) GradYukawa[CreateInvariantYukawa[Group,RepScalar,RepFermion,{{1,2,1},{True,False,True}}][[1]]];
YsffC=(y-I y5) GradYukawa[CreateInvariantYukawa[Group,RepScalar,RepFermion,{{1,1,2},{True,False,True}}][[1]]];
ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];

(*one scalar and one temporal scalar: each tensor has a single component*)
one[n_]:=SparseArray[ConstantArray[1,ConstantArray[1,n]]];
sol5=Dimension5Matching[Table[\[Alpha][k] one[{3,5,5,5,3,3,3,3}[[k]]],{k,8}],Range[8],Array[\[Alpha],8],3][[1]];
sol5P=Dimension5Matching["PHYS",Table[\[Beta][k] one[{3,5,5,5,3}[[k]]],{k,5}],Range[5],Array[\[Beta],5],3][[1]];

AppendTo[testList,
VerificationTest[
	Array[\[Alpha],8]/.sol5,
	m\[Psi] Zf[3,0] {
		-8/3 g1^2 y Sqrt[T],
		-480 y (y^2+y5^2)^2 T^(3/2),
		192 g1^2 y (y^2+y5^2) T^(3/2),
		0,
		20/3 y (y^2+y5^2) Sqrt[T],
		-16/3 g1^2 y Sqrt[T],
		0,
		-8 I g1^2 y5 Sqrt[T]},
	SameTest->(Union[Simplify[#1-#2]]==={0}&),
	TestID->"U1-off-shell"
]];

AppendTo[testList,
VerificationTest[
	Array[\[Beta],5]/.sol5P,
	m\[Psi] Zf[3,0] {
		-8/3 g1^2 y Sqrt[T],
		-160 y (y^2+y5^2) (3 y^2+3 y5^2-10 \[Lambda]) T^(3/2),
		64 g1^2 y (3 y^2+3 y5^2-2 \[Lambda]) T^(3/2),
		0,
		-8 I g1^2 y5 Sqrt[T]},
	SameTest->(Union[Simplify[#1-#2]]==={0}&),
	TestID->"U1-on-shell"
]];


(* ::Section:: *)
(*U(1) without fermion masses*)


(*
	Without fermion masses, ODIM5Fermion does not evaluate the fermion loops.
*)
\[Mu]IJ=0 \[Mu]IJ;\[Mu]IJC=0 \[Mu]IJC;
ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];

AppendTo[testList,
VerificationTest[
	Union[Flatten[Normal[ODIM5[#,3]]&/@Range[8]]],
	{0},
	TestID->"U1-without-masses"
]];


(* ::Section:: *)
(*SU(3) with a Dirac triplet and a real octet*)


(*
	A Dirac triplet with the mass m\[Psi] and the Yukawa coupling yo Sigma^a T^a of a real octet Sigma, whose components
	are those of the gauge fields. At d = 3, the coefficients of d^abc Sigma^a F^b F^c and d^abc (D^2 Sigma)^a A0^b A0^c
	agree with the matching of Matchotter, and those of Sigma^a (D^2 A0)^b A0^c vanish.
*)
Group={"SU3"};CouplingName={gs};RepAdjoint={{1,1}};RepScalar={{{{1,1}},"R"}};RepFermion={{{{1,0}},"L"},{{{1,0}},"R"}};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJC,\[Mu]IJ,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
gvss=gvvv;
tGen=Normal[Coefficient[gvff,gs]][[All,1;;3,1;;3]];
fG=Normal[Coefficient[gvvv,gs]];
dG=Table[2 Tr[(tGen[[i]] . tGen[[j]]+tGen[[j]] . tGen[[i]]) . tGen[[k]]],{i,8},{j,8},{k,8}]//Simplify;
(*the fermions 1 to 3 are psi_L, 4 to 6 the conjugates of psi_R*)
\[Mu]IJ=SparseArray[m\[Psi] ArrayFlatten[{{0,IdentityMatrix[3]},{IdentityMatrix[3],0}}]];
\[Mu]IJC=\[Mu]IJ;
Ysff=SparseArray[Table[yo ArrayFlatten[{{0,tGen[[i]]},{Transpose[tGen[[i]]],0}}],{i,8}]];
YsffC=SparseArray[Conjugate[Normal[Ysff]]/.Conjugate[yo]->yo];
ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
sol5=Dimension5Matching[{\[Alpha][1] dG,\[Alpha][6] dG,\[Alpha][7,1] fG+\[Alpha][7,2] dG},{1,6,7},{\[Alpha][1],\[Alpha][6],\[Alpha][7,1],\[Alpha][7,2]},3][[1]];

AppendTo[testList,
VerificationTest[
	{\[Alpha][1],\[Alpha][6],\[Alpha][7,1],\[Alpha][7,2]}/.sol5,
	gs^2 m\[Psi] yo Sqrt[T] Zf[3,0] {-2/3,-4/3,0,0},
	SameTest->(Union[Simplify[#1-#2]]==={0}&),
	TestID->"SU3-octet"
]];


(* ::Section:: *)
(*Report*)


report=TestReport[testList];
Print[report];
Print[Column[{#["TestID"],#["Outcome"]}&/@Values[report["TestResults"]]]];
Exit[If[report["AllTestsSucceeded"],0,1]]
