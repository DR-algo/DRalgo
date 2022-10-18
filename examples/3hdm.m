(* ::Package:: *)

(*Quit[];*)


SetDirectory[NotebookDirectory[]];
$LoadGroupMath=True;
<<../DRalgo.m


(* ::Chapter:: *)
(*3HDM with \[CapitalSigma](36) symmetry*)


(*See 1410.6139 [hep-ph]*)


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
(*The input for the gauge interactions toDRalgo are then given by*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion3Gen,RepScalar];


(* ::Text:: *)
(*The first element is the vector self - interaction matrix :*)


InputInv={{1,1},{True,False}}; (*\[Phi]1 \[Phi]1^+*)
MassTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;
InputInv={{2,2},{True,False}}; (*\[Phi]2 \[Phi]2^+*)
MassTerm2=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;
InputInv={{3,3},{True,False}}; (*\[Phi]3 \[Phi]3^+*)
MassTerm3=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;
InputInv={{1,2},{True,False}}; (*\[Phi]1\[Phi]2^+*)
MassTerm4Re=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;
InputInv={{2,1},{True,False}};(*\[Phi]2\[Phi]1^+*)
MassTerm4Im=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;
InputInv={{1,3},{True,False}}; (*\[Phi]1\[Phi]3^+*)
MassTerm5Re=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;
InputInv={{3,1},{True,False}};(*\[Phi]3\[Phi]1^+*)
MassTerm5Im=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;
InputInv={{2,3},{True,False}}; (*\[Phi]2\[Phi]3^+*)
MassTerm6Re=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;
InputInv={{3,2},{True,False}};(*\[Phi]3\[Phi]2^+*)
MassTerm6Im=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;


VMass=m1(MassTerm1+MassTerm2+MassTerm3);


\[Mu]ij=GradMass[VMass]//Simplify//SparseArray;


QuarticTerm1=\[Lambda]H1 (MassTerm1+MassTerm2+MassTerm3)^2;
QuarticTerm2=\[Lambda]H2(
	+MassTerm4Re*MassTerm4Im
	+MassTerm5Re*MassTerm5Im
	+MassTerm6Re*MassTerm6Im
	);
QuarticTerm3=-\[Lambda]H2(
	+MassTerm1*MassTerm2
	+MassTerm1*MassTerm3
	+MassTerm2*MassTerm3
	);
QuarticTerm4=\[Lambda]H3(
	+(MassTerm4Re-MassTerm5Re)*(MassTerm4Im-MassTerm5Im)
	+(MassTerm4Re-MassTerm6Re)*(MassTerm4Im-MassTerm6Im)
	+(MassTerm5Re-MassTerm6Re)*(MassTerm5Im-MassTerm6Im)
	);


VQuartic=QuarticTerm1+QuarticTerm2+QuarticTerm3+QuarticTerm4;


\[Lambda]4=GradQuartic[VQuartic];


InputInv={{1,1,2},{False,False,True}}; 
YukawaDoublet1=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv][[1]]//Simplify;


InputInv={{2,1,2},{False,False,True}}; 
YukawaDoublet2=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv][[1]]//Simplify;


InputInv={{3,1,2},{False,False,True}}; 
YukawaDoublet3=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv][[1]]//Simplify;


Ysff=-GradYukawa[yt1*YukawaDoublet1+yt2*YukawaDoublet2+yt3*YukawaDoublet3];


YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt1>0,yt2>0,yt3>0}]];


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
