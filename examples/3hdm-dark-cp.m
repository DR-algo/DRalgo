(* ::Package:: *)

(*Quit[];*)


SetDirectory[NotebookDirectory[]];
DRalgo`DRalgo`$LoadGroupMath=True;
<<../Kernel/DRalgo.m


(* ::Chapter:: *)
(*Z2-symmetric 3HDM with a CP-violating dark sector*)


(*See 2511.04636 [hep-ph]*)


(* ::Section::Closed:: *)
(*Model*)


Group={"SU3","SU2","U1"};
RepAdjoint={{1,1},{2},0};
HiggsDoublet1={{{0,0},{1},1/2},"C"};
HiggsDoublet2={{{0,0},{1},1/2},"C"};
HiggsDoublet3={{{0,0},{1},1/2},"C"};
RepScalar={HiggsDoublet1,HiggsDoublet2,HiggsDoublet3};
CouplingName={g3,g2,g1};


Rep1={{{1,0},{1},1/6},"L"};
Rep2={{{1,0},{0},2/3},"R"};
Rep3={{{1,0},{0},-1/3},"R"};
Rep4={{{0,0},{1},-1/2},"L"};
Rep5={{{0,0},{0},-1},"R"};
RepFermion1Gen={Rep1,Rep2,Rep3,Rep4,Rep5};


RepFermion3Gen={RepFermion1Gen,RepFermion1Gen,RepFermion1Gen}//Flatten[#,1]&;


(* ::Text:: *)
(*The input for the gauge interactions to DRalgo are then given by*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion3Gen,RepScalar];


(*
	The potential of Eq. (2.1) of 2511.04636. Under Z2 = diag(-1,-1,1), the inert doublets \[Phi]1, \[Phi]2 are odd and \[Phi]3 is even.
	The complex parameters \[Mu]12^2, \[Lambda]1, \[Lambda]2, \[Lambda]3 of the paper are split into real and imaginary parts, \[Mu]12sqR+I \[Mu]12sqI etc.
	The names \[Lambda]1R, \[Lambda]3R, ... avoid the tensors \[Lambda]1 and \[Lambda]3, and \[Mu]3sq avoids the 3d scale mu3.
	The dark democracy limit of the paper is \[Mu]2sq=\[Mu]1sq, \[Lambda]3=\[Lambda]2, \[Lambda]31=\[Lambda]23, and \[Lambda]31p=\[Lambda]23p.
*)
(*bilinear[i,j] = \[Phi]i^+ \[Phi]j*)
bilinear[i_,j_]:=CreateInvariant[Group,RepScalar,{{j,i},{True,False}}][[1]]//Simplify//FullSimplify;


VMass=(
	-\[Mu]1sq bilinear[1,1]-\[Mu]2sq bilinear[2,2]-\[Mu]3sq bilinear[3,3]
	-\[Mu]12sqR (bilinear[1,2]+bilinear[2,1])-I \[Mu]12sqI (bilinear[1,2]-bilinear[2,1])
	);


\[Mu]ij=GradMass[VMass]//Simplify//SparseArray;


(*
	As in 2511.04636, the Z2-invariant operators (\[Phi]1^+\[Phi]2)(\[Phi]1^+\[Phi]1), (\[Phi]1^+\[Phi]2)(\[Phi]2^+\[Phi]2), (\[Phi]3^+\[Phi]1)(\[Phi]2^+\[Phi]3), and (\[Phi]1^+\[Phi]2)(\[Phi]3^+\[Phi]3)
	are not included. Under the phase rotations of the three doublets, they carry the odd charges of \[Mu]12^2, while all
	quartic and Yukawa couplings carry even charges. Therefore, the \[Beta]-functions do not generate these operators,
	and the potential is closed under renormalization.
*)
VQuartic=(
	+\[Lambda]11 bilinear[1,1]^2+\[Lambda]22 bilinear[2,2]^2+\[Lambda]33 bilinear[3,3]^2
	+\[Lambda]12 bilinear[1,1]bilinear[2,2]+\[Lambda]23 bilinear[2,2]bilinear[3,3]+\[Lambda]31 bilinear[3,3]bilinear[1,1]
	+\[Lambda]12p bilinear[1,2]bilinear[2,1]+\[Lambda]23p bilinear[2,3]bilinear[3,2]+\[Lambda]31p bilinear[3,1]bilinear[1,3]
	+\[Lambda]1R (bilinear[1,2]^2+bilinear[2,1]^2)+I \[Lambda]1I (bilinear[1,2]^2-bilinear[2,1]^2)
	+\[Lambda]2R (bilinear[2,3]^2+bilinear[3,2]^2)+I \[Lambda]2I (bilinear[2,3]^2-bilinear[3,2]^2)
	+\[Lambda]3R (bilinear[3,1]^2+bilinear[1,3]^2)+I \[Lambda]3I (bilinear[3,1]^2-bilinear[1,3]^2)
	);


\[Lambda]4=GradQuartic[VQuartic];


(*Type-I Yukawa sector of Eq. (2.4): only the active doublet \[Phi]3 couples to fermions, here the top quark*)
InputInv={{3,1,2},{False,False,True}};
YukawaDoublet3=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv][[1]]//Simplify;


Ysff=-GradYukawa[yt YukawaDoublet3];


YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt>0}]];


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
(*Integrating out the temporal scalars*)


(*As in 2511.04636, the temporal scalars are integrated out and all three doublets remain light*)
PerformDRsoft[{}];


PrintCouplingsUS[]


PrintScalarMassUS["LO"]
PrintScalarMassUS["NLO"]
