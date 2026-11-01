(* ::Package:: *)

(*Quit[];*)


If[$InputFileName=="",
	SetDirectory[NotebookDirectory[]],
	SetDirectory[DirectoryName[$InputFileName]]
];
DRalgo`DRalgo`$LoadGroupMath=True;
<<../Kernel/DRalgo.m


(* ::Title:: *)
(*Scotogenic model*)


(*
	See hep-ph/0601225: the Standard Model with a second Higgs doublet \[Eta] and
	three right-handed neutrinos N with Majorana masses.
	The background fields of both doublets are kept in the soft theory,
	and the heavy fields are integrated out with their field-dependent masses.
*)


(* ::Section::Closed:: *)
(*Model*)


Group={"SU3","SU2","U1"};
RepAdjoint={{1,1},{2},0};
Higgs={{{0,0},{1},1/2},"C"};
Eta={{{0,0},{1},1/2},"C"};
RepScalar={Higgs,Eta};
CouplingName={g3,g2,g1};


Rep1={{{1,0},{1},1/6},"L"};
Rep2={{{1,0},{0},2/3},"R"};
Rep3={{{1,0},{0},-1/3},"R"};
Rep4={{{0,0},{1},-1/2},"L"};
Rep5={{{0,0},{0},-1},"R"};
Rep6={{{0,0},{0},0},"R"};
RepFermion1Gen={Rep1,Rep2,Rep3,Rep4,Rep5,Rep6};


RepFermion3Gen={RepFermion1Gen,RepFermion1Gen,RepFermion1Gen}//Flatten[#,1]&;


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion3Gen,RepScalar];


InputInv={{1,1},{True,False}}; (*\[Phi]^+\[Phi]*)
\[Phi]\[Phi]=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{2,2},{True,False}}; (*\[Eta]^+\[Eta]*)
\[Eta]\[Eta]=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{1,2},{True,False}}; (*\[Phi]^+\[Eta]*)
\[Phi]\[Eta]=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{2,1},{True,False}}; (*\[Eta]^+\[Phi]*)
\[Eta]\[Phi]=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;


VMass=msq1*\[Phi]\[Phi]+msq2*\[Eta]\[Eta];
\[Mu]ij=GradMass[VMass]//Simplify;


VQuartic=(
	+\[Lambda]s1*\[Phi]\[Phi]^2
	+\[Lambda]s2*\[Eta]\[Eta]^2
	+\[Lambda]s3*\[Phi]\[Phi]*\[Eta]\[Eta]
	+\[Lambda]s4*\[Phi]\[Eta]*\[Eta]\[Phi]
	+\[Lambda]s5/2*(\[Phi]\[Eta]^2+\[Eta]\[Phi]^2)
	);
\[Lambda]4=GradQuartic[VQuartic];


(*The right-handed neutrinos N1, N2, N3 are the fermions 6, 12, 18*)
InputInv={{6,6},{True,True}};
NN1=CreateInvariantFermion[Group,RepFermion3Gen,InputInv][[1]]//Simplify;
InputInv={{12,12},{True,True}};
NN2=CreateInvariantFermion[Group,RepFermion3Gen,InputInv][[1]]//Simplify;
InputInv={{18,18},{True,True}};
NN3=CreateInvariantFermion[Group,RepFermion3Gen,InputInv][[1]]//Simplify;


FermionMasses=1/2*MJ1*NN1+1/2*MJ2*NN2+1/2*MJ3*NN3;
\[Mu]IJ=GradMassFermion[FermionMasses];
\[Mu]IJC=SparseArray[Simplify[Conjugate[\[Mu]IJ]//Normal,Assumptions->{MJ1>0,MJ2>0,MJ3>0}]];


(*The top-quark Yukawa coupling, and the coupling of \[Eta] to the third-generation lepton doublet L3 (16) and N3 (18)*)
InputInv={{1,1,2},{False,False,True}};
YukawaDoublet1=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv][[1]]//Simplify;
InputInv={{2,16,18},{False,False,True}};
YukawaScoto=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv][[1]]//Simplify;


Ysff=-yt*GradYukawa[YukawaDoublet1]-yN3*GradYukawa[YukawaScoto];
YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt>0,yN3>0}]];


(* ::Section:: *)
(*Dimensional Reduction*)


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PerformDRhard[]


PrintCouplings[]
BetaFunctions4D[]


PrintScalarMass["LO"]
PrintScalarMass["NLO"]


PrintDebyeMass["LO"]
PrintDebyeMass["NLO"]


PrintTemporalScalarCouplings[]


(* ::Section:: *)
(*Effective potential with both background fields*)


(*
	Computations are done in the soft theory.
	The doublets live on 1,...,4 and 5,...,8, followed by the temporal scalars:
	SU(3) on 9,...,16, SU(2) on 17, 18, 19, and U(1) on 20
*)
UseSoftTheory[];
PrintScalarRepPositions[]


(*The background fields \[Phi]1 of \[Phi] and \[Phi]2 of \[Eta]*)
vev=SparseArray[{4->\[Phi]1,8->\[Phi]2},{20}];
DefineVEVS[vev];


ScalarMass=PrintTensorsVEV[1]//Normal;
VectorMass=PrintTensorsVEV[2]//Normal;


(*
	The mass matrices are block diagonal: the background fields mix the components
	{1,5}, {2,6}, {3,7}, and {4,8} of the two doublets, the temporal scalars {19,20},
	and the neutral vectors {11,12}. The rotation to the mass eigenstates diagonalizes each block,
	whose eigenvalues are the masses.
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


{DSRot,ScalarMassDiag}=RotationToMassBasis[ScalarMass,{{1,5},{2,6},{3,7},{4,8},{19,20}}];
{DVRot,VectorMassDiag}=RotationToMassBasis[VectorMass,{{11,12}}];


(*
	The two-loop potential needs only the components \[Lambda]_iikk of the rotated quartic couplings.
	With FastRotation->True, only these are rotated, which is much faster for this rotation.
*)
RotateTensorsCustomMass[DSRot,DVRot,ScalarMassDiag,VectorMassDiag,FastRotation->True];


CalculatePotential[];


PrintEffectivePotential["LO"]
PrintEffectivePotential["NLO"]


(*The two-loop potential is a large expression*)
VNNLO=PrintEffectivePotential["NNLO"];


(* ::Section:: *)
(*Integrate out the heavy fields in a Higgs Effective-type Theory (HET)*)


(*
	The temporal scalars, the W bosons, and the Z boson are integrated out
	with their field-dependent masses. The photon is massless.
*)
HeavyVectors=Flatten@Position[Diagonal[VectorMassDiag],m_/;m=!=0,{1},Heads->False]


PrepareHET[Range[9,20],HeavyVectors];
CalculatePotentialHET[];


PrintActionHET["LO"]
PrintActionHET["NLO"]
