(* ::Package:: *)

(*Quit[];*)


SetDirectory[NotebookDirectory[]];
$LoadGroupMath=True;
<<../DRalgo.m


(* ::Chapter:: *)
(*SM+2 real singlets*)


(* ::Section::Closed:: *)
(*Model*)


Group={"SU3","SU2","U1"};
RepAdjoint={{1,1},{2},0};
HiggsDoublet1={{{0,0},{1},1/2},"C"};
HiggsSinglet1={{{0,0},{0},0},"R"};
HiggsSinglet2={{{0,0},{0},0},"R"};
RepScalar={HiggsDoublet1,HiggsSinglet1,HiggsSinglet2};
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


InputInv={{1,1},{True,False}}; (*\[Phi]\[Phi]^+*)
MassTerm1=CreateInvariant[Group,RepScalar,InputInv]//Simplify;
InputInv={{2,2},{True,True}}; (*S1^2*)
MassTerm2=CreateInvariant[Group,RepScalar,InputInv]//Simplify;
InputInv={{3,3},{True,True}}; (*S2^2*)
MassTerm3=CreateInvariant[Group,RepScalar,InputInv]//Simplify;


VMass=(
	+m1*MassTerm1
	+1/2*mS1*MassTerm2
	+1/2*mS2*MassTerm3
	);


\[Mu]ij=GradMass[VMass[[1]]]//Simplify; 


QuarticTerm1=MassTerm1[[1]]^2; (*[(\[Phi]\[Phi]^+)]^2*)
QuarticTerm2=MassTerm2[[1]]^2; (* S1^4*)
QuarticTerm3=MassTerm3[[1]]^2; (* S2^4*)
QuarticTerm4=MassTerm1[[1]]*MassTerm2[[1]]; (* [\[Phi]\[Phi]^+]S1^2*)
QuarticTerm5=MassTerm1[[1]]*MassTerm3[[1]]; (* [\[Phi]\[Phi]^+]S2^2*)
QuarticTerm6=MassTerm2[[1]]*MassTerm3[[1]]; (* S1^2S2^2*)


InputInv={{2,2,2,3},{True,True,True,True}}; (*S1^3S2*)
QuarticTerm7=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{3,3,3,2},{True,True,True,True}}; (*S2^3S1*)
QuarticTerm8=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;


VQuartic=(
	+\[Lambda]H*QuarticTerm1
	+\[Lambda]S1/4!*QuarticTerm2
	+\[Lambda]S2/4!*QuarticTerm3
	+\[Lambda]M1/2*QuarticTerm4
	+\[Lambda]M2/2*QuarticTerm5
	+\[Lambda]M3/2*QuarticTerm6
	+\[Lambda]M4/4!*QuarticTerm7
	+\[Lambda]M5/4!*QuarticTerm8
	);


\[Lambda]4=GradQuartic[VQuartic];


InputInv={{1,1,2},{True,False,True}}; (*\[Phi]\[Phi]^+S1*)
CubicTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{1,1,3},{True,False,True}}; (*\[Phi]\[Phi]^+S2*)
CubicTerm2=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{2,2,2},{True,True,True}}; (*S1^3*)
CubicTerm3=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{3,3,3},{True,True,True}}; (*S2^3*)
CubicTerm4=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{2,3,3},{True,True,True}}; (*S1*S2^2*)
CubicTerm5=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{2,2,3},{True,True,True}}; (*S1^2*S2*)
CubicTerm6=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;


VCubic=(
	+\[Lambda]C1*CubicTerm1
	+\[Lambda]C2*CubicTerm2
	+\[Lambda]C3/3!*CubicTerm3
	+\[Lambda]C4/3!*CubicTerm4
	+\[Lambda]C5/2!*CubicTerm5
	+\[Lambda]C6/2!*CubicTerm6
	);


\[Lambda]3=GradCubic[VCubic];


InputInv={{2},{True}}; (*S1*)
TadpoleTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{3},{True}}; (*S2*)
TadpoleTerm2=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;


VTadpole=(
	+\[Lambda]Tad1*TadpoleTerm1
	+\[Lambda]Tad2*TadpoleTerm2
	);


\[Lambda]1=GradTadpole[VTadpole];


InputInv={{1,1,2},{False,False,True}}; 
YukawaDoublet1=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv]//Simplify;


Ysff=-yt*GradYukawa[YukawaDoublet1[[1]]];


YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt>0}]];


(* ::Section:: *)
(*Dimensional Reduction*)


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PerformDRhard[]


BetaFunctions4D[]


PrintCouplings[]


PrintTadpoles["LO"]


PrintTemporalScalarCouplings[]


PrintDebyeMass["LO"]
PrintDebyeMass["NLO"]
PrintScalarMass["LO"]
PrintScalarMass["LO"]


PrintPressure["LO"]
PrintPressure["NLO"]
PrintPressure["NNLO"]


(* ::Text:: *)
(*All Scalars active:*)


PerformDRsoft[{}];


PrintCouplingsUS[]


PrintTadpolesUS["LO"]


PrintScalarMassUS["LO"]//Simplify


(* ::Text:: *)
(*2 Scalars active:*)


PerformDRsoft[{6}]


PrintCouplingsUS[]


PrintScalarMassUS["LO"]//Simplify


(* ::Text:: *)
(*1 Scalar active:*)


PerformDRsoft[{5,6}];


PrintCouplingsUS[]


PrintScalarMassUS["LO"]//Simplify
