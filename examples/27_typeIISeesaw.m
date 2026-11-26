(* ::Package:: *)

(*Quit[];*)


If[$InputFileName=="",
	SetDirectory[NotebookDirectory[]],
	SetDirectory[DirectoryName[$InputFileName]]
];
DRalgo`DRalgo`$LoadGroupMath=True;
<<../Kernel/DRalgo.m


(* ::Title:: *)
(*Type-II seesaw model*)


(*
	See 2611.01634 [hep-ph]: the Standard Model with a complex SU(2) triplet \[CapitalDelta] of hypercharge 1
	(Y=2 in the convention Q=T3+Y/2 of the paper). With the doublet \[CapitalPhi] and the 2x2 matrix
	\[CapitalDelta]={{\[CapitalDelta]^+/Sqrt[2],H^++},{\[CapitalDelta]^0,-\[CapitalDelta]^+/Sqrt[2]}}, the potential of eq. (2.3) reads
		V = -msq \[CapitalPhi]^\[Dagger]\[CapitalPhi] + Msq Tr(\[CapitalDelta]^\[Dagger]\[CapitalDelta]) + [\[Mu] \[CapitalPhi]^T i\[Tau]2 \[CapitalDelta]^\[Dagger] \[CapitalPhi] + h.c.]
		  + \[Lambda]s1 (\[CapitalPhi]^\[Dagger]\[CapitalPhi])^2 + \[Lambda]s2 [Tr(\[CapitalDelta]^\[Dagger]\[CapitalDelta])]^2 + \[Lambda]s3 Tr(\[CapitalDelta]^\[Dagger]\[CapitalDelta]\[CapitalDelta]^\[Dagger]\[CapitalDelta])
		  + \[Lambda]s4 (\[CapitalPhi]^\[Dagger]\[CapitalPhi]) Tr(\[CapitalDelta]^\[Dagger]\[CapitalDelta]) + \[Lambda]s5 \[CapitalPhi]^\[Dagger]\[CapitalDelta]\[CapitalDelta]^\[Dagger]\[CapitalPhi],
	where msq, Msq, \[Lambda]s1,...,\[Lambda]s5 are m^2, M^2, \[Lambda]1,...,\[Lambda]5 of the paper.
	The triplet-lepton Yukawa couplings are neglected, as in the paper.
*)


(* ::Section::Closed:: *)
(*Model*)


Group={"SU3","SU2","U1"};
RepAdjoint={{1,1},{2},0};
HiggsDoublet={{{0,0},{1},1/2},"C"};
HiggsTriplet={{{0,0},{2},1},"C"};
RepScalar={HiggsDoublet,HiggsTriplet};
CouplingName={g3,g2,g1};


Rep1={{{1,0},{1},1/6},"L"};
Rep2={{{1,0},{0},2/3},"R"};
Rep3={{{1,0},{0},-1/3},"R"};
Rep4={{{0,0},{1},-1/2},"L"};
Rep5={{{0,0},{0},-1},"R"};
RepFermion1Gen={Rep1,Rep2,Rep3,Rep4,Rep5};


RepFermion3Gen={RepFermion1Gen,RepFermion1Gen,RepFermion1Gen}//Flatten[#,1]&;


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion3Gen,RepScalar];


InputInv={{1,1},{True,False}}; (*\[CapitalPhi]^\[Dagger]\[CapitalPhi]*)
\[CapitalPhi]\[CapitalPhi]=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{2,2},{True,False}}; (*Tr(\[CapitalDelta]^\[Dagger]\[CapitalDelta])*)
\[CapitalDelta]\[CapitalDelta]=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;


VMass=-msq*\[CapitalPhi]\[CapitalPhi]+Msq*\[CapitalDelta]\[CapitalDelta];
\[Mu]ij=GradMass[VMass]//Simplify//SparseArray;


(*
	GroupMath writes the triplet in the weight basis, here (-H^++,\[CapitalDelta]^+,\[CapitalDelta]^0), and its invariants
	are normalized differently from the traces of the paper. The factors below are fixed by
	writing the operators of the paper out in components:
		\[CapitalPhi]^T i\[Tau]2 \[CapitalDelta]^\[Dagger] \[CapitalPhi] = -3^(1/4)/Sqrt[2] (\[CapitalPhi]\[CapitalPhi]\[CapitalDelta]^\[Dagger])_GroupMath,
		Tr(\[CapitalDelta]^\[Dagger]\[CapitalDelta]\[CapitalDelta]^\[Dagger]\[CapitalDelta]) = [Tr(\[CapitalDelta]^\[Dagger]\[CapitalDelta])]^2 - 1/2 |Tr(\[CapitalDelta]\[CapitalDelta])|^2,
		\[CapitalPhi]^\[Dagger]\[CapitalDelta]\[CapitalDelta]^\[Dagger]\[CapitalPhi] = 1/2 (\[CapitalPhi]^\[Dagger]\[CapitalPhi]) Tr(\[CapitalDelta]^\[Dagger]\[CapitalDelta]) + 1/2 \[CapitalPhi]^\[Dagger][\[CapitalDelta],\[CapitalDelta]^\[Dagger]]\[CapitalPhi].
	With this sign of the cubic term, the neutral components are the real parts
	\[CapitalPhi]^0=\[Phi]/Sqrt[2] and \[CapitalDelta]^0=\[Delta]/Sqrt[2] of the scalars 2 and 7, and the potential on these
	backgrounds is eq. (2.4) of the paper.
*)
InputInv={{1,1,2},{True,True,False}}; (*\[CapitalPhi]^T i\[Tau]2 \[CapitalDelta]^\[Dagger] \[CapitalPhi]*)
\[CapitalPhi]\[CapitalPhi]\[CapitalDelta]c=-3^(1/4)/Sqrt[2]*CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{1,1,2},{False,False,True}}; (*h.c.*)
\[CapitalPhi]c\[CapitalPhi]c\[CapitalDelta]=-3^(1/4)/Sqrt[2]*CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;


VCubic=\[Mu]*(\[CapitalPhi]\[CapitalPhi]\[CapitalDelta]c+\[CapitalPhi]c\[CapitalPhi]c\[CapitalDelta]);
\[Lambda]3=GradCubic[VCubic];


InputInv={{2,2},{True,True}}; (*Tr(\[CapitalDelta]\[CapitalDelta]) up to a sign*)
\[CapitalDelta]Sq=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{2,2},{False,False}}; (*Tr(\[CapitalDelta]^\[Dagger]\[CapitalDelta]^\[Dagger]) up to a sign*)
\[CapitalDelta]Sqc=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{1,1,2,2},{True,False,True,False}}; (*(\[CapitalPhi]^\[Dagger]\[CapitalPhi])Tr(\[CapitalDelta]^\[Dagger]\[CapitalDelta]) and \[CapitalPhi]^\[Dagger][\[CapitalDelta],\[CapitalDelta]^\[Dagger]]\[CapitalPhi]/Sqrt[2]*)
\[CapitalPhi]\[CapitalPhi]\[CapitalDelta]\[CapitalDelta]=CreateInvariant[Group,RepScalar,InputInv]//Simplify;


TrDDDD=\[CapitalDelta]\[CapitalDelta]^2-1/2*\[CapitalDelta]Sq*\[CapitalDelta]Sqc; (*Tr(\[CapitalDelta]^\[Dagger]\[CapitalDelta]\[CapitalDelta]^\[Dagger]\[CapitalDelta])*)
\[CapitalPhi]DD\[CapitalPhi]=1/2*\[CapitalPhi]\[CapitalPhi]\[CapitalDelta]\[CapitalDelta][[1]]+1/Sqrt[2]*\[CapitalPhi]\[CapitalPhi]\[CapitalDelta]\[CapitalDelta][[2]]; (*\[CapitalPhi]^\[Dagger]\[CapitalDelta]\[CapitalDelta]^\[Dagger]\[CapitalPhi]*)


VQuartic=(
	+\[Lambda]s1*\[CapitalPhi]\[CapitalPhi]^2
	+\[Lambda]s2*\[CapitalDelta]\[CapitalDelta]^2
	+\[Lambda]s3*TrDDDD
	+\[Lambda]s4*\[CapitalPhi]\[CapitalPhi]*\[CapitalDelta]\[CapitalDelta]
	+\[Lambda]s5*\[CapitalPhi]DD\[CapitalPhi]
	);
\[Lambda]4=GradQuartic[VQuartic];


InputInv={{1,1,2},{False,False,True}};
YukawaDoublet=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv][[1]]//Simplify;


Ysff=-yt*GradYukawa[YukawaDoublet];
YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt>0}]];


(* ::Section:: *)
(*Dimensional Reduction*)


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PerformDRhard[]


PrintCouplings[]


(*The beta functions and the anomalous dimensions (Landau gauge) of eqs. (3.7)-(3.11) of the paper*)
BetaFunctions4D[]


PosScalar=PrintScalarRepPositions[];
Table[AnomDim4D["S",{a,a}],{a,PosScalar}]


(*The thermal masses of eq. (4.7) are the T^2 terms of the LO masses*)
PrintScalarMass["LO"]
PrintScalarMass["NLO"]


(*The Debye masses of eq. (4.8)*)
PrintDebyeMass["LO"]
PrintDebyeMass["NLO"]


PrintTemporalScalarCouplings[]


(* ::Subsection:: *)
(*Comparison with 2611.01634*)


(*
	BetaFunctions4D[] prints dX/dlog\[Mu] for the squared gauge couplings and for all other parameters,
	which are eqs. (3.7) and (3.9)-(3.11) of the paper with g'=g1.
	AnomDim4D[] prints dlogR/dlog\[Mu] of the renormalized fields, which is minus the \[Gamma] of eq. (3.8).
	The thermal masses of eq. (4.7) are msq-msq3d and Msq3d-Msq at LO,
	and the Debye masses \[CapitalPi]_WL and \[CapitalPi]_BL of eq. (4.8) are musqSU2 and musqU1.
	All differences vanish.
*)
betaPaper={
	g1^2->2*47/6*g1^4,
	g2^2->2*(-5/2)*g2^4,
	g3^2->2*(-7)*g3^4,
	yt->yt*(9/2*yt^2-17/12*g1^2-9/4*g2^2-8*g3^2),
	\[Lambda]s1->3/8*g1^4+3/4*g1^2*g2^2+9/8*g2^4-(3*g1^2+9*g2^2)*\[Lambda]s1+24*\[Lambda]s1^2+3*\[Lambda]s4^2+3*\[Lambda]s4*\[Lambda]s5
		+5/4*\[Lambda]s5^2+12*\[Lambda]s1*yt^2-6*yt^4,
	\[Lambda]s2->6*g1^4-12*g1^2*g2^2+15*g2^4-(12*g1^2+24*g2^2)*\[Lambda]s2+28*\[Lambda]s2^2+24*\[Lambda]s2*\[Lambda]s3
		+6*\[Lambda]s3^2+2*\[Lambda]s4^2+2*\[Lambda]s4*\[Lambda]s5,
	\[Lambda]s3->24*g1^2*g2^2-6*g2^4-(12*g1^2+24*g2^2)*\[Lambda]s3+24*\[Lambda]s2*\[Lambda]s3+18*\[Lambda]s3^2+\[Lambda]s5^2,
	\[Lambda]s4->3*g1^4-6*g1^2*g2^2+6*g2^4-(15/2*g1^2+33/2*g2^2)*\[Lambda]s4+12*\[Lambda]s1*\[Lambda]s4+4*\[Lambda]s1*\[Lambda]s5
		+16*\[Lambda]s2*\[Lambda]s4+12*\[Lambda]s3*\[Lambda]s4+4*\[Lambda]s4^2+6*\[Lambda]s2*\[Lambda]s5+2*\[Lambda]s3*\[Lambda]s5+\[Lambda]s5^2+6*\[Lambda]s4*yt^2,
	\[Lambda]s5->12*g1^2*g2^2-(15/2*g1^2+33/2*g2^2)*\[Lambda]s5+4*\[Lambda]s5^2+4*\[Lambda]s1*\[Lambda]s5+4*\[Lambda]s2*\[Lambda]s5
		+8*\[Lambda]s3*\[Lambda]s5+8*\[Lambda]s4*\[Lambda]s5+6*\[Lambda]s5*yt^2,
	msq->msq*(12*\[Lambda]s1+6*yt^2-9/2*g2^2-3/2*g1^2)-3*Msq*(2*\[Lambda]s4+\[Lambda]s5)-12*\[Mu]^2,
	Msq->Msq*(16*\[Lambda]s2+12*\[Lambda]s3-12*g2^2-6*g1^2)-2*msq*(2*\[Lambda]s4+\[Lambda]s5)+4*\[Mu]^2,
	\[Mu]->\[Mu]*(4*\[Lambda]s1+4*\[Lambda]s4+6*\[Lambda]s5+6*yt^2-21/2*g2^2-9/2*g1^2)
	};
Simplify[Lookup[BetaFunctions4D[],Keys[betaPaper]]-Values[betaPaper]/(16\[Pi]^2)]


\[Gamma]Paper={3*yt^2-9/4*g2^2-3/4*g1^2,-6*g2^2-3*g1^2}/(16\[Pi]^2);
Simplify[Table[AnomDim4D["S",{a,a}],{a,PosScalar}]+\[Gamma]Paper]


\[CapitalPi]Paper={
	T^2/16*(3*g2^2+g1^2+8*\[Lambda]s1+4*\[Lambda]s4+2*\[Lambda]s5+4*yt^2),
	T^2/12*(6*g2^2+3*g1^2+8*\[Lambda]s2+6*\[Lambda]s3+2*\[Lambda]s4+\[Lambda]s5)
	};
Simplify[({msq-msq3d,Msq3d-Msq}/.PrintScalarMass["LO"])-\[CapitalPi]Paper]


Simplify[({musqSU2,musqU1}/.PrintDebyeMass["LO"])-{5/2*g2^2*T^2,17/6*g1^2*T^2}]


(* ::Section:: *)
(*Integrate out the temporal scalars*)


PerformDRsoft[{}];


PrintCouplingsUS[]


PrintScalarMassUS["LO"]
PrintScalarMassUS["NLO"]


BetaFunctions3DUS[]


(* ::Section:: *)
(*Integrate out the temporal scalars and the triplet*)


(*
	For a soft triplet mass, Msq3d ~ (gT)^2, only the doublet remains in the ultrasoft theory.
	The tree-level exchange of the triplet through the cubic coupling \[Mu] shifts the doublet quartic.
	The triplet lives on 5,...,10.
*)
PerformDRsoft[{5,6,7,8,9,10}];


PrintCouplingsUS[]


PrintScalarMassUS["LO"]
PrintScalarMassUS["NLO"]


(* ::Section:: *)
(*Effective potential on the neutral backgrounds*)


(*
	Computations are done in the soft theory.
	The doublet lives on 1,...,4 and the triplet on 5,...,10, followed by the temporal scalars:
	SU(3) on 11,...,18, SU(2) on 19, 20, 21, and U(1) on 22
*)
UseSoftTheory[];
PrintScalarRepPositions[]


(*The background fields \[Phi] of \[CapitalPhi]^0 and \[Delta] of \[CapitalDelta]^0, as in eq. (2.4) of the paper*)
vev=SparseArray[{2->\[Phi],7->\[Delta]},{22}];
DefineVEVS[vev];


ScalarMass=PrintTensorsVEV[1]//Normal;
VectorMass=PrintTensorsVEV[2]//Normal;


(*
	The mass matrices are block diagonal: the background fields mix the components {1,6}, {2,7}, {3,9},
	and {4,10} of the doublet and the triplet, the temporal scalars {21,22}, and the neutral vectors {11,12}.
	The rotation to the mass eigenstates diagonalizes each block, whose eigenvalues are the masses.
*)
RotationToMassBasis[mass_,blocks_]:=Module[
	{rotation=IdentityMatrix[Length[mass]],masses=Diagonal[mass],values,vectors},
	Do[
		{values,vectors}=Eigensystem[mass[[b,b]]];
		rotation[[b,b]]=Transpose[Simplify[#/Sqrt[# . #]]&/@vectors];
		masses[[b]]=values;
	,{b,blocks}];
	{SparseArray[rotation],SparseArray[DiagonalMatrix[masses]]}
];


{DSRot,ScalarMassDiag}=RotationToMassBasis[ScalarMass,{{1,6},{2,7},{3,9},{4,10},{21,22}}];
{DVRot,VectorMassDiag}=RotationToMassBasis[VectorMass,{{11,12}}];


(*The two-loop potential needs only the components \[Lambda]_iikk of the rotated quartic couplings*)
RotateTensorsCustomMass[DSRot,DVRot,ScalarMassDiag,VectorMassDiag,FastRotation->True];


CalculatePotential[];


PrintEffectivePotential["LO"]
PrintEffectivePotential["NLO"]


(*The two-loop potential is a large expression*)
VNNLO=PrintEffectivePotential["NNLO"];
