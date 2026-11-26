(* ::Package:: *)

(*Quit[];*)


If[$InputFileName=="",
	SetDirectory[NotebookDirectory[]],
	SetDirectory[DirectoryName[$InputFileName]]
];
DRalgo`DRalgo`$LoadGroupMath=True;
DRalgo`DRalgo`$GroupMathMultipleModels=True; (*Put this if you want to create multiple model-files with the same kernel*)
<<../Kernel/DRalgo.m


(* ::Chapter:: *)
(*Three-Higgs-doublet model: the operator H^6*)


(*
	See 2605.15176 [hep-ph] for the dimension-6 operator matching and its physical basis.
	Three Higgs doublets with a general quartic,
		V \[Superset] (1/4!) Sum L4[i,j,k,l] (\[Phi]i^\[Dagger]\[Phi]j)(\[Phi]k^\[Dagger]\[Phi]l),
	with random integer entries, and one quark generation whose top Yukawa coupling yt3 couples to \[Phi]3.
	The one-loop Wilson coefficient of the six-doublet operator in the physical basis, ODIM6["PHYS", 8, d],
	is matched onto the flavor structures C6 of L4, the identity and the Yukawa vector,
		O = Sum C6[i1,j1,i2,j2,i3,j3] (\[Phi]i1^\[Dagger]\[Phi]j1)(\[Phi]i2^\[Dagger]\[Phi]j2)(\[Phi]i3^\[Dagger]\[Phi]j3).
	A random L4 makes the structures of one order in the couplings independent unless they are identical for every
	L4, so that their coefficients are those of a general quartic.
	The 12 real scalars give a tensor of the operator with 12^6 components, of which 12376 are independent:
	Dimension6Matching compares only those. Its run time is that of the operator, 8 to 12 minutes and up to 12 GB of
	memory on a laptop, mostly the contractions of the loops with six scalars; the matching takes 15 seconds.
*)


(* ::Section::Closed:: *)
(*Model*)


Group={"SU3","SU2","U1"};
RepAdjoint={{1,1},{2},0};
HiggsDoublet1={{{0,0},{1},1/2},"C"};
HiggsDoublet2={{{0,0},{1},1/2},"C"};
HiggsDoublet3={{{0,0},{1},1/2},"C"};
RepScalar={HiggsDoublet1,HiggsDoublet2,HiggsDoublet3};

GaugeCouplings={g3,g2,g1};


Rep1={{{1,0},{1},1/6},"L"};
Rep2={{{1,0},{0},2/3},"R"};
RepFermion1Gen={Rep1,Rep2};


{gvvv,gvff,gvss,\[CapitalLambda]1,\[CapitalLambda]3,\[CapitalLambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,GaugeCouplings,RepFermion1Gen,RepScalar];


(*
	The notation for the doublet products is \[Phi]ij = \[Phi]i^\[Dagger]\[Phi]j: the first index is the conjugated doublet.
*)
InputInv={{1,1},{True,False}};
\[Phi]11=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;
InputInv={{2,2},{True,False}};
\[Phi]22=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;
InputInv={{3,3},{True,False}};
\[Phi]33=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;
InputInv={{1,2},{True,False}};
\[Phi]21=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;
InputInv={{2,1},{True,False}};
\[Phi]12=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;


(* complex parameters with their real and imaginary parts *)
\[Mu]12sq = \[Mu]12sqRe + I*\[Mu]12sqIm;
\[Mu]12sqConj = Conjugate[\[Mu]12sq]//ComplexExpand;

VMass=-\[Mu]1sq*\[Phi]11 - \[Mu]2sq*\[Phi]22 - \[Mu]3sq*\[Phi]33 - \[Mu]12sq*\[Phi]12 - \[Mu]12sqConj*\[Phi]21 // Simplify // Expand;
\[Mu]ij=GradMass[VMass]//Simplify//SparseArray;


InputInv={{3,1,2},{False,False,True}};
YukawaDoublet3=CreateInvariantYukawa[Group,RepScalar,RepFermion1Gen,InputInv][[1]]//Simplify;
Ysff=-GradYukawa[yt3*YukawaDoublet3];
YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt3>0}]];


(*
	The quartic: L4[i,j,k,l] = L4[k,l,i,j] with random integers. The real components of the doublets are
	\[Phi]i = Viab[[i,\[Alpha],a]] x_a, so that \[Phi]i^\[Dagger]\[Phi]j = \[Phi]ijab[[i,j,a,b]] x_a x_b.
*)
canon[i1_, i2_, i3_, i4_] := If[OrderedQ[{{i1, i2}, {i3, i4}}], l[i1, i2, i3, i4], l[i3, i4, i1, i2]];
L4 = Array[canon, {3,3,3,3}];
allLSymbols = DeleteDuplicates@Cases[L4, _l, Infinity];
SeedRandom[1234];
L4 = L4/.Thread[allLSymbols -> RandomInteger[{-9, 9}, Length[allLSymbols]]];

Viab={{{1,0,I,0,0,0,0,0,0,0,0,0},{0,1,0,I,0,0,0,0,0,0,0,0}},{{0,0,0,0,1,0,I,0,0,0,0,0},{0,0,0,0,0,1,0,I,0,0,0,0}},{{0,0,0,0,0,0,0,0,1,0,I,0},{0,0,0,0,0,0,0,0,0,1,0,I}}}/Sqrt[2];
ViabConj=Conjugate[Viab];
\[Phi]ijab=Transpose[Contract[ViabConj,Viab,{{2,5}}],{1,3,2,4}];
\[CapitalLambda]4=SparseArray[SymmetrizeTensor[Contract[L4,\[Phi]ijab,\[Phi]ijab,{{1,5},{2,6},{3,9},{4,10}}],{{1,2,3,4}}]];


(* ::Section:: *)
(*Dimensional Reduction*)


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[CapitalLambda]1,\[CapitalLambda]3,\[CapitalLambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False,Mode->2,Normalization4D->False,AutoRG->False];
PerformDRhard[];


(* ::Section:: *)
(*Operator basis*)


(*
	The flavor structures C6[i1,j1,i2,j2,i3,j3] of the one-loop operator by order in the couplings: products of L4
	(indices c u c u, c for the conjugated doublet of a bilinear), of \[Delta] (c u) and of the Yukawa vector e = (0,0,1),
	whose indices are joined by lines from u to c. A line is the propagator of a doublet, \[Delta], summed over m, n, p,
	or for a line through the top loop the projector onto \[Phi]3, which sets its doublet to 3. The six free indices,
	three c and three u, form the three bilinears.
		family   structures            lines                 coefficients
		g6       \[Delta] \[Delta] \[Delta]                 -                     g^6
		y6       e (x6)                -                     yt^6
		Lg4      L4 \[Delta]                  -                     \[Lambda] g^4
		Ly4      L4 e (x3)             L4 - top loop         \[Lambda] yt^4
		LLg2     L4 L4                 one \[Delta] line            \[Lambda]^2 g^2
		LLy2     L4 L4                 one line through \[Phi]3  \[Lambda]^2 yt^2
		LLL      L4 L4 L4              three \[Delta] lines,         \[Lambda]^3 (triangles, and the chains L4 - L4 = L4
		                               connected             from the equations of motion)
	Each family lists its structures that are linearly independent for this L4 (94 in all). The basis of the first
	version of this example had 45 structures: it lacked the chains of \[Lambda]^3, the families LLy2 and Ly4 (O(6,15) and
	O(6,4) of the top loop), and could not reproduce \[Lambda]^3, \[Lambda]^2 yt^2 and \[Lambda] yt^4.
*)

structures = {
	"g6" -> {
		\[Delta][i1,j1] \[Delta][i2,j2] \[Delta][i3,j3],
		\[Delta][i1,j1] \[Delta][i2,j3] \[Delta][i3,j2]
	},
	"y6" -> {
		e[i1] e[i2] e[i3] e[j1] e[j2] e[j3]
	},
	"Lg4" -> {
		L[i1,j1,i2,j2] \[Delta][i3,j3],
		L[i1,j1,i2,j3] \[Delta][i3,j2],
		L[i1,j2,i2,j1] \[Delta][i3,j3]
	},
	"Ly4" -> {
		L[3,j1,i1,j2] e[i2] e[i3] e[j3],
		L[3,j2,i1,j1] e[i2] e[i3] e[j3],
		L[3,j2,i1,j3] e[i2] e[i3] e[j1],
		L[i1,3,i2,j1] e[i3] e[j2] e[j3],
		L[i1,3,i2,j2] e[i3] e[j1] e[j3],
		L[i1,3,i2,j3] e[i3] e[j1] e[j2]
	},
	"LLg2" -> {
		L[i1,m,i2,j1] L[m,j2,i3,j3],
		L[i1,m,i2,j1] L[m,j3,i3,j2],
		L[i1,m,i2,j2] L[m,j1,i3,j3],
		L[i1,m,i2,j3] L[m,j1,i3,j2],
		L[i1,m,i2,j2] L[m,j3,i3,j1]
	},
	"LLy2" -> {
		L[i1,3,i2,j1] L[3,j2,i3,j3],
		L[i1,3,i2,j1] L[3,j3,i3,j2],
		L[i1,3,i2,j2] L[3,j1,i3,j3],
		L[i1,3,i2,j3] L[3,j1,i3,j2],
		L[i1,3,i2,j2] L[3,j3,i3,j1]
	},
	"LLL" -> {
		L[i1,m,i2,n] L[m,p,n,j1] L[p,j2,i3,j3],
		L[i1,m,i2,n] L[m,p,n,j1] L[p,j3,i3,j2],
		L[i1,m,i2,n] L[m,p,n,j2] L[p,j1,i3,j3],
		L[i1,m,i2,n] L[m,p,n,j3] L[p,j1,i3,j2],
		L[i1,m,i2,n] L[m,p,n,j2] L[p,j3,i3,j1],
		L[p,m,i1,n] L[m,j1,n,j2] L[i2,p,i3,j3],
		L[p,m,i1,n] L[m,j1,n,j3] L[i2,p,i3,j2],
		L[p,m,i1,n] L[m,j2,n,j1] L[i2,p,i3,j3],
		L[p,m,i1,n] L[m,j3,n,j1] L[i2,p,i3,j2],
		L[p,m,i1,n] L[m,j2,n,j3] L[i2,p,i3,j1],
		L[p,m,i1,n] L[m,p,i2,j1] L[n,j2,i3,j3],
		L[p,m,i1,n] L[m,p,i2,j1] L[n,j3,i3,j2],
		L[p,m,i1,n] L[m,p,i2,j2] L[n,j1,i3,j3],
		L[p,m,i1,n] L[m,p,i2,j3] L[n,j1,i3,j2],
		L[p,m,i1,n] L[m,p,i2,j2] L[n,j3,i3,j1],
		L[i1,m,p,n] L[m,p,i2,j1] L[n,j2,i3,j3],
		L[i1,m,p,n] L[m,p,i2,j1] L[n,j3,i3,j2],
		L[i1,m,p,n] L[m,p,i2,j2] L[n,j1,i3,j3],
		L[i1,m,p,n] L[m,p,i2,j3] L[n,j1,i3,j2],
		L[i1,m,p,n] L[m,p,i2,j2] L[n,j3,i3,j1],
		L[i1,m,i2,n] L[m,p,i3,j1] L[n,j2,p,j3],
		L[i1,m,i2,n] L[m,p,i3,j1] L[n,j3,p,j2],
		L[i1,m,i2,n] L[m,p,i3,j2] L[n,j1,p,j3],
		L[i1,m,i2,n] L[m,p,i3,j3] L[n,j1,p,j2],
		L[i1,m,i2,n] L[m,p,i3,j2] L[n,j3,p,j1],
		L[p,m,i1,n] L[m,j1,i2,p] L[n,j2,i3,j3],
		L[p,m,i1,n] L[m,j1,i2,p] L[n,j3,i3,j2],
		L[p,m,i1,n] L[m,j2,i2,p] L[n,j1,i3,j3],
		L[p,m,i1,n] L[m,j3,i2,p] L[n,j1,i3,j2],
		L[p,m,i1,n] L[m,j2,i2,p] L[n,j3,i3,j1],
		L[i1,m,p,n] L[m,j1,i2,p] L[n,j2,i3,j3],
		L[i1,m,p,n] L[m,j1,i2,p] L[n,j3,i3,j2],
		L[i1,m,p,n] L[m,j2,i2,p] L[n,j1,i3,j3],
		L[i1,m,p,n] L[m,j3,i2,p] L[n,j1,i3,j2],
		L[i1,m,p,n] L[m,j2,i2,p] L[n,j3,i3,j1],
		L[i1,m,i2,n] L[m,j1,i3,p] L[n,j2,p,j3],
		L[i1,m,i2,n] L[m,j1,i3,p] L[n,j3,p,j2],
		L[i1,m,i2,n] L[m,j2,i3,p] L[n,j1,p,j3],
		L[i1,m,i2,n] L[m,j3,i3,p] L[n,j1,p,j2],
		L[i1,m,i2,n] L[m,j2,i3,p] L[n,j3,p,j1],
		L[n,m,p,j1] L[m,n,i1,j2] L[i2,p,i3,j3],
		L[n,m,p,j1] L[m,n,i1,j3] L[i2,p,i3,j2],
		L[n,m,p,j2] L[m,n,i1,j1] L[i2,p,i3,j3],
		L[n,m,p,j3] L[m,n,i1,j1] L[i2,p,i3,j2],
		L[n,m,p,j2] L[m,n,i1,j3] L[i2,p,i3,j1],
		L[p,m,n,j1] L[m,n,i1,j2] L[i2,p,i3,j3],
		L[p,m,n,j1] L[m,n,i1,j3] L[i2,p,i3,j2],
		L[p,m,n,j2] L[m,n,i1,j1] L[i2,p,i3,j3],
		L[p,m,n,j3] L[m,n,i1,j1] L[i2,p,i3,j2],
		L[p,m,n,j2] L[m,n,i1,j3] L[i2,p,i3,j1],
		L[i1,m,n,j1] L[m,n,p,j2] L[i2,p,i3,j3],
		L[i1,m,n,j1] L[m,n,p,j3] L[i2,p,i3,j2],
		L[i1,m,n,j2] L[m,n,p,j1] L[i2,p,i3,j3],
		L[i1,m,n,j3] L[m,n,p,j1] L[i2,p,i3,j2],
		L[i1,m,n,j2] L[m,n,p,j3] L[i2,p,i3,j1],
		L[p,m,i1,j1] L[m,n,i2,j2] L[n,p,i3,j3],
		L[p,m,i1,j1] L[m,n,i2,j3] L[n,p,i3,j2],
		L[p,m,i1,j3] L[m,n,i2,j1] L[n,p,i3,j2],
		L[i1,m,p,j1] L[m,n,i2,j2] L[n,p,i3,j3],
		L[i1,m,p,j1] L[m,n,i2,j3] L[n,p,i3,j2],
		L[i1,m,p,j2] L[m,n,i2,j1] L[n,p,i3,j3],
		L[i1,m,p,j3] L[m,n,i2,j1] L[n,p,i3,j2],
		L[i1,m,p,j2] L[m,n,i2,j3] L[n,p,i3,j1],
		L[i1,m,p,j1] L[m,n,i2,j2] L[n,j3,i3,p],
		L[i1,m,p,j1] L[m,n,i2,j3] L[n,j2,i3,p],
		L[i1,m,p,j2] L[m,n,i2,j1] L[n,j3,i3,p],
		L[i1,m,p,j3] L[m,n,i2,j1] L[n,j2,i3,p],
		L[i1,m,p,j2] L[m,n,i2,j3] L[n,j1,i3,p],
		L[i1,m,i2,j1] L[m,n,p,j2] L[n,j3,i3,p],
		L[i1,m,i2,j3] L[m,n,p,j1] L[n,j2,i3,p],
		L[i1,m,p,j1] L[m,j2,i2,n] L[n,j3,i3,p],
		L[i1,m,p,j1] L[m,j3,i2,n] L[n,j2,i3,p]
	}
};

structures=Association[structures];


\[Delta][a_Integer,b_Integer]:=KroneckerDelta[a,b];
e[a_Integer]:=KroneckerDelta[a,3];
L[a_Integer,b_Integer,c_Integer,k_Integer]:=L4[[a,b,c,k]];
FlavorTensor[expr_]:=With[{internal=Union[Cases[expr,m|n|p,{0,Infinity}]]},
	Table[If[internal==={},expr,Sum@@Prepend[{#,3}&/@internal,expr]],{i1,3},{j1,3},{i2,3},{j2,3},{i3,3},{j3,3}]
];

factor=<|"y6"->yt3^6,"Ly4"->yt3^4,"LLy2"->yt3^2|>;
coefficients=Flatten[Table[\[Beta][family,k],{family,Keys[structures]},{k,Length[structures[[Key[family]]]]}]];
C6tensor=Sum[\[Beta][family,k]Lookup[factor,family,1]FlavorTensor[structures[[Key[family],k]]],{family,Keys[structures]},{k,Length[structures[[Key[family]]]]}];


(*
	The group tensor of the operator: C6 averaged over the permutations of the bilinears, contracted with the
	bilinears of the real components, and symmetrized.
*)
perms=Permutations[{{1,2},{3,4},{5,6}}];
C6tensorSymm=Total[TensorTranspose[C6tensor,Flatten[#]]&/@perms]/Length[perms];
Tens8NS=Contract[C6tensorSymm,\[Phi]ijab,\[Phi]ijab,\[Phi]ijab,{{1,7},{2,8},{3,11},{4,12},{5,15},{6,16}}]//SparseArray;
Tens8=Factorial[6]SymmetrizeTensor[Tens8NS,{{1,2,3,4,5,6}}];


(* ::Section:: *)
(*Matching*)


(*
	For a numeric quartic, structures of different orders in the couplings can be linearly dependent: the
	coefficients of 7 L4^3 structures are left free, and the other coefficients depend on them. They are fixed by
	requiring that each coefficient \[Beta][family, k], times the couplings of its structure (factor), contains only
	the monomials of its family.
*)
sol=First@Quiet[Dimension6Matching["PHYS",{Tens8},{8},coefficients,3],Solve::svars];
free=Complement[coefficients,sol[[All,1]]]


(*
	The family of a monomial by its powers of the gauge couplings and of yt3. The gauge couplings are taken from
	GaugeCouplings: after AllocateTensors, g2 refers to a symbol of GroupMath.
*)
familyOf[monomial_]:=Lookup[
	<|{0,0}->"LLL",{2,0}->"LLg2",{4,0}->"Lg4",{6,0}->"g6",{0,2}->"LLy2",{0,4}->"Ly4",{0,6}->"y6"|>,
	Key[{Total[Exponent[monomial,#]&/@GaugeCouplings],Exponent[monomial,yt3]}]
];
terms[expr_]:=Replace[Expand[expr],{s_Plus:>List@@s,0->{},s_:>{s}}];
monomialOf[term_]:=Times@@Cases[Replace[term,{t_Times:>List@@t,t_:>{t}}],Except[_?NumericQ|_\[Tau]]];
structureFactor[c_]:=Lookup[factor,c[[1]],1];
monomials=GroupBy[DeleteDuplicates[monomialOf/@Flatten[terms[(#[[2]]/.Thread[free->0])structureFactor[#[[1]]]]&/@sol]],familyOf];
ansatz=Table[With[{own=Lookup[monomials,f[[1]],{}]},f->Sum[\[Tau][f,k]own[[k]],{k,Length[own]}]/structureFactor[f]],{f,free}];
conditions=Flatten[Table[
	Values[GroupBy[
		Select[{monomialOf[#],#}&/@terms[(rule[[2]]/.ansatz)structureFactor[rule[[1]]]],familyOf[#[[1]]]=!=rule[[1,1]]&],
		First->(#[[2]]/#[[1]]&),Total]],
	{rule,sol}]];
fixed=ansatz/.First[Solve[Thread[conditions==0],Cases[ansatz,_\[Tau],Infinity]//DeleteDuplicates]];
result=Factor/@Association[Join[sol/.fixed,fixed]];


Normal[result]//TableForm
