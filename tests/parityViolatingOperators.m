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
(*Parity-violating operators*)


(*
	The operators O(6,18), O(6,19), O(6,20) and O(5,8) contain the dual field strength Ftilde_i = (1/2) eps_ijk F_jk
	of eq. (dim_6_operators) of 2605.15176. Their coefficients agree with a direct evaluation of the fermion loops
	with explicit Dirac matrices, for static plane waves and for constant non-abelian backgrounds, and those of
	O(6,19) and O(5,8) agree with the matching of Matchotter (2604.21972).
*)


testList={};


(* ::Section:: *)
(*Chiral U(1)*)


(*
	A left-handed fermion psi of charge 1, a right-handed fermion chi of charge 0, and a complex scalar H of
	charge 1, with V = lam |H|^4 and the Yukawa coupling y H psi^+ chi + h.c. The real components of H enter
	O(6,19) through eps2. Matchotter gives no O(6,18), whose coefficient is that of the direct evaluation of
	the A0 A_j A_l function of psi at third order in the momenta.
*)
Group={"U1"};CouplingName={e1};RepAdjoint={0};RepScalar={{{1},"C"}};RepFermion={{{1},"L"},{{0},"R"}};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJC,\[Mu]IJ,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
\[Lambda]4=GradQuartic[\[Lambda] CreateInvariant[Group,RepScalar,{{1,1},{True,False}}][[1]]^2];
Ysff=-GradYukawa[y CreateInvariantYukawa[Group,RepScalar,RepFermion,{{1,1,2},{True,False,True}}][[1]]];
YsffC=SparseArray[Simplify[Conjugate[Normal[Ysff]],Assumptions->{y>0}]];
ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];

one[n_]:=SparseArray[ConstantArray[1,ConstantArray[1,n]]];
eps2={{0,1},{-1,0}};
O19=SparseArray[\[Alpha][19,1] TensorProduct[IdentityMatrix[2],{{1}}]+\[Alpha][19,2] TensorProduct[eps2,{{1}}]];
sol6=Dimension6Matching[{\[Alpha][18] one[3],O19,\[Alpha][20] one[4]},{18,19,20},{\[Alpha][18],\[Alpha][19,1],\[Alpha][19,2],\[Alpha][20]},3][[1]];

AppendTo[testList,
VerificationTest[
	{\[Alpha][18],\[Alpha][19,1],\[Alpha][19,2],\[Alpha][20]}/.sol6,
	Zf[3,0] {2/3 I e1^3 Sqrt[T],0,4/3 I e1^2 y^2 T,0},
	SameTest->(Union[Simplify[#1-#2]]==={0}&),
	TestID->"U1-chiral"
]];


(* ::Section:: *)
(*U(1) with a scalar and a pseudoscalar Yukawa coupling*)


(*
	A Dirac fermion of charge 1 and two neutral real scalars, with the Yukawa couplings y1 phi1 psibar psi and
	y5 phi2 psibar i gamma5 psi, i.e. y1 + I y5 for phi_a psi_R^+ psi_L. The coefficient of O(6,19) agrees with
	Matchotter. Away from d = 3, the fermion loops of O(6,2), O(6,5) and O(6,7) are not symmetric in the two
	scalars, and only their symmetric parts, which enter the operators, are matched onto the symmetric tensors.
	Their values at d = 3 are those of the fermion loops without this projection.
*)
RepScalar={{{0},"R"},{{0},"R"}};RepFermion={{{1},"L"},{{1},"R"}};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJC,\[Mu]IJ,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
yuk[s_,f_]:=GradYukawa[CreateInvariantYukawa[Group,RepScalar,RepFermion,{{s,f[[1]],f[[2]]},{True,False,True}}][[1]]];
Ysff=y1 yuk[1,{2,1}]+I y5 yuk[2,{2,1}];
YsffC=y1 yuk[1,{1,2}]-I y5 yuk[2,{1,2}];
ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];

O19=SparseArray[TensorProduct[Array[a19,{2,2}],{{1}}]];
OS=SparseArray[TensorProduct[{{a11,a12},{a12,a22}},{{1}}]];

AppendTo[testList,
VerificationTest[
	Flatten[Array[a19,{2,2}]]/.Dimension6Matching[{O19},{19},Flatten[Array[a19,{2,2}]],3][[1]],
	8 I e1^2 T y1 y5 Zf[3,0] {0,1,-1,0},
	SameTest->(Union[Simplify[#1-#2]]==={0}&),
	TestID->"U1-pseudoscalar-O(6,19)"
]];

symmetricPair[dim_]:=Table[{a11,a12,a22}/.Dimension6Matching[{OS},{x},{a11,a12,a22},dim][[1]],{x,{2,5,7}}];

AppendTo[testList,
VerificationTest[
	symmetricPair[d],
	e1^2 T Zf[3,0] {-4/3 (3d-7),-2/3 (d^2+8d-57),-2/3 (d^2-9)} ConstantArray[{y1^2,0,y5^2},3],
	SameTest->(Union[Flatten[Simplify[#1-#2]]]==={0}&),
	TestID->"U1-pseudoscalar-symbolic-d"
]];

AppendTo[testList,
VerificationTest[
	symmetricPair[3],
	e1^2 T Zf[3,0] {-8/3,16,0} ConstantArray[{y1^2,0,y5^2},3],
	SameTest->(Union[Flatten[Simplify[#1-#2]]]==={0}&),
	TestID->"U1-pseudoscalar-d=3"
]];


(* ::Section:: *)
(*Chiral SU(2)xU(1) doublet on constant backgrounds*)


(*
	For constant gauge fields A^a_mu, the fermion loops of the operators without scalars give the sixth order of
	-Tr log(Dslash) in a_mu = A^a_mu gvff^a, i.e. (1/6) Sum-Integral_P tr (S V)^6 with S = -I Pslash/P^2 and
	V = gamma_mu PL (-I a_mu) for the left-handed Weyl fermions, at T = 1 and with the couplings set to 1.
	For a left-handed doublet of hypercharge 1/2, both O(6,18) and O(6,20) contribute, and the comparison fixes
	their normalization independently of the tensors of DRalgo. The covariant derivatives and field strengths
	follow eq. (dim_6_operators) of 2605.15176: (D_k X)^a = g^abc A^b_k X^c and F^a_ij = g^abc A^b_i A^c_j.
*)
Group={"SU2","U1"};CouplingName={gw,gY};RepAdjoint={{2},0};RepScalar={};RepFermion={{{{1},1/2},"L"}};
pureOps={1,3,8,12,13,16,17,18,20};
numeric[x_]:=x/.{s_Symbol/;MemberQ[{"gw","gY","T","xi"},SymbolName[s]]:>1}/.
	{z_[3,0]/;SymbolName[z]==="Zb":>Zeta[3]/(128 Pi^4),z_[3,0]/;SymbolName[z]==="Zf":>7 Zeta[3]/(128 Pi^4)};
pureGaugeTensors[fermions_]:=(
	{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJC,\[Mu]IJ,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,fermions,RepScalar];
	ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
	Association[Table[x->numeric[Normal[ODIM6[x,3]]],{x,pureOps}]]
);
withFermions=pureGaugeTensors[RepFermion];
gens=numeric[Normal[gvff]];
fabc=numeric[Normal[gvvv]];
nV=Length[fabc];
withoutFermions=pureGaugeTensors[{}];
Oc=Association[Table[x->withFermions[x]-withoutFermions[x],{x,pureOps}]];

eps3=Normal[LeviCivitaTensor[3]];
dralgoLagrangian[A_]:=Module[{A0,Ai,cov,F,Ft,DA0,DDA0,DF},
	A0=A[[All,4]];Ai=Table[A[[All,k1]],{k1,3}];
	cov[k_,X_]:=Table[Sum[fabc[[i1,i2,i3]]Ai[[k,i2]]X[[i3]],{i2,nV},{i3,nV}],{i1,nV}];
	F=Table[Table[Sum[fabc[[i1,i2,i3]]Ai[[k1,i2]]Ai[[k2,i3]],{i2,nV},{i3,nV}],{i1,nV}],{k1,3},{k2,3}];
	Ft=Table[1/2 Sum[eps3[[k1,k2,k3]]F[[k2,k3]],{k2,3},{k3,3}],{k1,3}];
	DA0=Table[cov[k1,A0],{k1,3}];
	DDA0=Sum[cov[k1,DA0[[k1]]],{k1,3}];
	DF=Table[Sum[cov[k1,F[[k1,k2]]],{k1,3}],{k2,3}];
	<|
	1->1/6 Sum[Oc[1][[i1,i2,i3]]F[[k1,k2,i1]]F[[k2,k3,i2]]F[[k3,k1,i3]],{i1,nV},{i2,nV},{i3,nV},{k1,3},{k2,3},{k3,3}],
	3->1/4 Sum[Oc[3][[i1,i2,i3,i4]]F[[k1,k2,i1]]F[[k1,k2,i2]]A0[[i3]]A0[[i4]],{i1,nV},{i2,nV},{i3,nV},{i4,nV},{k1,3},{k2,3}],
	8->1/4 Sum[Oc[8][[i1,i2,i3,i4]]DA0[[k1,i1]]DA0[[k1,i2]]A0[[i3]]A0[[i4]],{i1,nV},{i2,nV},{i3,nV},{i4,nV},{k1,3}],
	12->1/720 Total[Flatten[Oc[12] TensorProduct@@ConstantArray[A0,6]]],
	13->1/2 Sum[Oc[13][[i1,i2]]DF[[k2,i1]]DF[[k2,i2]],{i1,nV},{i2,nV},{k2,3}],
	16->-1/2 Sum[Oc[16][[i1,i2,i3]]DF[[k1,i1]](DA0[[k1,i2]]A0[[i3]]-DA0[[k1,i3]]A0[[i2]]),{i1,nV},{i2,nV},{i3,nV},{k1,3}],
	17->1/2 Sum[Oc[17][[i1,i2]]DDA0[[i1]]DDA0[[i2]],{i1,nV},{i2,nV}],
	18->Sum[Oc[18][[i1,i2,i3]]A0[[i1]]Ft[[k1,i2]]DF[[k1,i3]],{i1,nV},{i2,nV},{i3,nV},{k1,3}],
	20->1/2 Sum[Oc[20][[i1,i2,i3,i4]]A0[[i1]]A0[[i2]]DA0[[k1,i3]]Ft[[k1,i4]],{i1,nV},{i2,nV},{i3,nV},{i4,nV},{k1,3}]
	|>
];

(*Sum-Integral_{P} p0^(2 al) |p|^(2 m)/(P^2)^n over fermionic modes in d = 3, with the angular average of the p_i*)
zfI[n_,al_]:=(2^(2n-2al-3)-1) 2 (2 Pi)^(3-2n+2al)/(4 Pi)^(3/2) Gamma[n-3/2]/Gamma[n] Zeta[2n-2al-3];
sumIntegral[{k1_,k2_,k3_,k0_},n_]:=If[OddQ[k0]||OddQ[k1]||OddQ[k2]||OddQ[k3],0,
	Module[{al=k0/2,m=(k1+k2+k3)/2},
		Factorial2[k1-1]Factorial2[k2-1]Factorial2[k3-1]/Factorial2[2m+1] Sum[Binomial[m,j](-1)^j zfI[n-m+j,al+j],{j,0,m}]]];
pauli={{{0,1},{1,0}},{{0,-I},{I,0}},{{1,0},{0,-1}}};
gamma[k_]:=If[k<=3,ArrayFlatten[{{0,-I pauli[[k]]},{I pauli[[k]],0}}],ArrayFlatten[{{0,IdentityMatrix[2]},{IdentityMatrix[2],0}}]];
gamma5=gamma[1] . gamma[2] . gamma[3] . gamma[4];
fermionLoop[a_]:=Module[{p={p1,p2,p3,p0},SV},
	SV=Sum[KroneckerProduct[(-I Sum[p[[nu]]gamma[nu],{nu,4}]) . gamma[mu] . (IdentityMatrix[4]+gamma5)/2,-I a[[mu]]],{mu,4}];
	1/6 Total[(#2 sumIntegral[#1,6])&@@@CoefficientRules[Expand[Tr[MatrixPower[SV,6]]],p]]
];

backgrounds={{{1,0,-1,2},{-2,1,0,1},{0,2,1,-1},{1,1,-2,1}},{{0,-1,2,1},{1,2,-1,0},{2,0,1,-2},{-1,1,1,2}}};
backgroundTest=Table[
	Module[{L=dralgoLagrangian[A]},
		{Simplify[Total[L]-fermionLoop[Table[Sum[A[[i1,mu]]gens[[i1]],{i1,nV}],{mu,4}]]],L[18]=!=0,L[20]=!=0}],
	{A,backgrounds}];

AppendTo[testList,
VerificationTest[
	backgroundTest,
	{{0,True,True},{0,True,True}},
	TestID->"SU2xU1-chiral-doublet-backgrounds"
]];


(* ::Section:: *)
(*Report*)


report=TestReport[testList];
Print[report];
Print[Column[{#["TestID"],#["Outcome"]}&/@Values[report["TestResults"]]]];
Exit[If[report["AllTestsSucceeded"],0,1]]
