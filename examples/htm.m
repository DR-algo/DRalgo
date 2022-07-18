(* ::Package:: *)

(*Quit[];*)


SetDirectory[NotebookDirectory[]];
$LoadGroupMath=True;
$GroupMathMultipleModels=True; (*Put this if you want to create multiple model-files with the same kernel*)
<<../DRalgo.m


(* ::Chapter:: *)
(*SM+Real SU(2)Triplet*)


(*See 1802.10500 [hep-ph], *)


(* ::Section:: *)
(*Model*)


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


(* ::Text:: *)
(*The input for the gauge interactions toDRalgo are then given by*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion3Gen,RepScalar];


(* ::Text:: *)
(*The first element is the vector self - interaction matrix :*)


InputInv={{1,1},{True,False}}; (*\[Phi]\[Phi]^+*)
MassTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;
InputInv={{2,2},{True,True}}; (*\[CapitalSigma]^a\[CapitalSigma]^a*)
MassTerm2=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;


VMass=(
	+mH*MassTerm1
	+m\[CapitalSigma]/2*MassTerm2
	);


\[Mu]ij=GradMass[VMass]//Simplify//SparseArray;


QuarticTerm1=MassTerm1^2; (*[(\[Phi]\[Phi]^+)]^2*)
QuarticTerm2=MassTerm2^2; (*[\[CapitalSigma]^a\[CapitalSigma]^a]^2*)
QuarticTerm3=MassTerm1*MassTerm2; (*[\[Phi]\[Phi]^+][\[CapitalSigma]^a\[CapitalSigma]^a]*)


VQuartic=(
	+\[Lambda]*QuarticTerm1
	+b4/4*QuarticTerm2
	+a2/2*QuarticTerm3
	);


\[Lambda]4=GradQuartic[VQuartic];


InputInv={{1,1,2},{False,False,True}}; 
YukawaDoublet1=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv]//Simplify;


Ysff=-GradYukawa[yt*YukawaDoublet1[[1]]];


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


(* ::Section:: *)
(*Identify Temporal Couplings:*)


PosScalarsDRalgo=PrintScalarRepPositions[]; (*The position of all scalar reps*)
PosGaugeDRalgo=PrintGaugeRepPositions[]; (*The position of all scalar reps*)


(* ::Subsection:: *)
(*3 d Model*)


(*The convention follows eq 11 in 1802.10500*)


Group={"SU3","SU2","U1"};
RepAdjoint={{1,1},{2},0};


CouplingName={g3,gw,g1};


Higgs\[CapitalSigma]={{{0,0},{2},0},"R"}; (*This is the triplet rep*)
HiggsDoublet={{{0,0},{1},1/2},"C"}; (*This is the Higgs doublet*) 
HiggsA0={{{0,0},{2},0},"R"}; (*This is the SU(2) Temporal Scalar*)
HiggsG0={{{1,1},{0},0},"R"}; (*This is the SU(3) Temporal Scalar*)
HiggsB0={{{0,0},{0},0},"R"}; (*This is the U1(1) Temporal Scalar*)
RepScalar={HiggsG0,HiggsA0,HiggsB0,Higgs\[CapitalSigma],HiggsDoublet};


{gvvv2,gvff,gvss2,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJC,\[Mu]IJ,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,{},RepScalar];


PosScalarsNew=PrintScalarRepPositions[]; (*The position of the new scalar reps*)


InputInv={{2,4},{True,True}}; (*\[CapitalSigma]aAa0*)
MassTerm1=CreateInvariant[Group, RepScalar, InputInv][[1]] //Simplify;
InputInv={{4,4},{True,True}}; (*\[CapitalSigma]a\[CapitalSigma]a*)
MassTerm2=CreateInvariant[Group, RepScalar, InputInv][[1]] //Simplify;
InputInv={{2,2},{True,True}}; (*Aa0Aa0*)
MassTerm3=CreateInvariant[Group, RepScalar, InputInv][[1]] //Simplify;
VQuartic1=(
	+\[Delta]3*MassTerm2*MassTerm3
	+\[Delta]3p*MassTerm1^2
	);


InputInv={{1,1},{True,True}}; (*G0^aG0^a*)
MassTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{5,5},{True,False}}; (*\[Phi]\[Phi]^+*)
MassTerm2=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
VQuartic2=\[Omega]3*MassTerm1*MassTerm2;


InputInv={{2,3,5,5},{True,True,True,False}}; (*A0^a B0 \[Phi]\[Tau]^a\[Phi]^+*)
MassTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
VQuartic3=h3pp*MassTerm1;


InputInv={{3,3,5,5}, {True,True,True,False}}; (*\[Phi]^+\[Phi]^+B0^2*)
MassTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
VQuartic4=h3p*MassTerm1;


InputInv={{2,2,5,5},{True,True,True,False}}; (*\[Phi]^+\[Phi]^+A0^aA0^a*)
MassTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
VQuartic5=h3*MassTerm1;


\[Lambda]KCompare=GradQuartic[VQuartic1+VQuartic2+VQuartic3+VQuartic4+VQuartic5];


(* ::Subsection:: *)
(*Comparison*)


(* ::Subsubsection:: *)
(*\[Delta]3 and \[Delta]3p Couplings*)


(*First we need to pick out the components*)


(*The convention, in DRalgo, is that the first two indices are temporal bosons. And the last two scalar particles*)


(*Because the ordering is different we need to use the Postions printed before*)


(*We are after the rep A0^2\[CapitalSigma]^2*)
(*In DRalgo output this corresponds to gauge rep 2 and scalar rep 2*)
(*While with our rep {HiggsG0,HiggsA0,HiggsB0,Higgs\[CapitalSigma],HiggsDoublet}, this corresponds to reps 2 and 4*)


\[Lambda]KComparep=\[Lambda]KCompare[[PosScalarsNew[[2]],PosScalarsNew[[2]],PosScalarsNew[[4]],PosScalarsNew[[4]]]];


\[Lambda]KDRalgo=PrintTensorDRalgo[][[3]][[PosGaugeDRalgo[[2]],PosGaugeDRalgo[[2]],PosScalarsDRalgo[[2]],PosScalarsDRalgo[[2]]]];


CompareInvariants[\[Lambda]KComparep,\[Lambda]KDRalgo]


(* ::Subsection:: *)
(*\[Omega]3*)


(*First we need to pick out the components*)


(*The convention, in DRalgo, is that the first two indices are temporal bosons. And the last two scalar particles*)


(*Because the ordering is different we need to use the Postions printed before*)


(*We are after the rep G0^2(\[Phi]\[Phi]^+)*)
(*In DRalgo output this corresponds to gauge rep 1 and scalar rep 1*)
(*While with our rep {HiggsG0,HiggsA0,HiggsB0,Higgs\[CapitalSigma],HiggsDoublet}, this corresponds to reps 1 and 5*)


\[Lambda]KComparep=\[Lambda]KCompare[[PosScalarsNew[[1]],PosScalarsNew[[1]],PosScalarsNew[[5]],PosScalarsNew[[5]]]];


\[Lambda]KDRalgo=PrintTensorDRalgo[][[3]][[PosGaugeDRalgo[[1]],PosGaugeDRalgo[[1]],PosScalarsDRalgo[[1]],PosScalarsDRalgo[[1]]]];


CompareInvariants[\[Lambda]KComparep,\[Lambda]KDRalgo]//Expand


(* ::Subsection:: *)
(*h3pp*)


(*First we need to pick out the components*)


(*The convention, in DRalgo, is that the first two indices are temporal bosons. And the last two scalar particles*)


(*Because the ordering is different we need to use the Postions printed before*)


(*We are after the rep A0^a B0 \[Phi]\[Tau]^a\[Phi]^+*)
(*In DRalgo output this corresponds to gauge reps 2,3 and scalar rep 1*)
(*While with our rep {HiggsG0,HiggsA0,HiggsB0,Higgs\[CapitalSigma],HiggsDoublet}, this corresponds to reps 2,3 and 5*)


\[Lambda]KComparep=\[Lambda]KCompare[[PosScalarsNew[[2]],PosScalarsNew[[3]],PosScalarsNew[[5]],PosScalarsNew[[5]]]];


\[Lambda]KDRalgo=PrintTensorDRalgo[][[3]][[PosGaugeDRalgo[[2]],PosGaugeDRalgo[[3]],PosScalarsDRalgo[[1]],PosScalarsDRalgo[[1]]]];


I*Sqrt[3]^(-1/2)*h3pp/.CompareInvariants[\[Lambda]KComparep,\[Lambda]KDRalgo]//Expand


(* ::Subsection:: *)
(*h3p*)


(*First we need to pick out the components*)


(*The convention, in DRalgo, is that the first two indices are temporal bosons. And the last two scalar particles*)


(*Because the ordering is different we need to use the Postions printed before*)


(*We are after the rep B0^2 \[Phi]\[Phi]^+*)
(*In DRalgo output this corresponds to gauge rep 3 and scalar rep 1*)
(*While with our rep {HiggsG0,HiggsA0,HiggsB0,Higgs\[CapitalSigma],HiggsDoublet}, this corresponds to reps 3 and 5*)


\[Lambda]KComparep=\[Lambda]KCompare[[PosScalarsNew[[3]],PosScalarsNew[[3]],PosScalarsNew[[5]],PosScalarsNew[[5]]]];


\[Lambda]KDRalgo=PrintTensorDRalgo[][[3]][[PosGaugeDRalgo[[3]],PosGaugeDRalgo[[3]],PosScalarsDRalgo[[1]],PosScalarsDRalgo[[1]]]];


CompareInvariants[\[Lambda]KComparep,\[Lambda]KDRalgo]


(* ::Subsection:: *)
(*h3*)


(*First we need to pick out the components*)


(*The convention, in DRalgo, is that the first two indices are temporal bosons. And the last two scalar particles*)


(*Because the ordering is different we need to use the Postions printed before*)


(*We are after the rep A0^2 \[Phi]\[Phi]^+*)
(*In DRalgo output this corresponds to gauge rep 2 and scalar rep 1*)
(*While with our rep {HiggsG0,HiggsA0,HiggsB0,Higgs\[CapitalSigma],HiggsDoublet}, this corresponds to reps 2 and 5*)


\[Lambda]KComparep=\[Lambda]KCompare[[PosScalarsNew[[2]],PosScalarsNew[[2]],PosScalarsNew[[5]],PosScalarsNew[[5]]]];


\[Lambda]KDRalgo=PrintTensorDRalgo[][[3]][[PosGaugeDRalgo[[2]],PosGaugeDRalgo[[2]],PosScalarsDRalgo[[1]],PosScalarsDRalgo[[1]]]];


CompareInvariants[\[Lambda]KComparep,\[Lambda]KDRalgo]
