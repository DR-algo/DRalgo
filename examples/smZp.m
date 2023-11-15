(* ::Package:: *)

(*Quit[];*)


SetDirectory[NotebookDirectory[]];
$LoadGroupMath=True;
<<../DRalgo.m


(* ::Chapter:: *)
(*SM+Gauge B-L+Right-handed Neutrinos	*)


(*See 0212073 [hep-ph]*)


(* ::Section::Closed:: *)
(*Model*)


Group={"SU3","SU2","U1","U1"};
RepAdjoint={{1,1},{2},0,0};
HiggsDoublet={{{0,0},{1},1/2,0},"C"};
HiggsSinglet={{{0,0},{0},0,1/2},"C"};
RepScalar={HiggsDoublet,HiggsSinglet};
CouplingName={g3,g2,g1,gZ};


Rep1={{{1,0},{1},1/6,1/6},"L"};
Rep2={{{1,0},{0},2/3,1/6},"R"};
Rep3={{{1,0},{0},-1/3,1/6},"R"};
Rep4={{{0,0},{1},-1/2,-1/2},"L"};
Rep5={{{0,0},{0},-1,-1/2},"R"};
Rep6={{{0,0},{0},0,1/4},"R"}; (*Right-handed Neutrino*)
RepFermion1Gen={Rep1,Rep2,Rep3,Rep4,Rep5,Rep6};


RepFermion3Gen={RepFermion1Gen,RepFermion1Gen,RepFermion1Gen}//Flatten[#,1]&;


(* ::Text:: *)
(*The input for the gauge interactions toDRalgo are then given by*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion3Gen,RepScalar];


InputInv={{1,1},{True,False}}; (*\[Phi]\[Phi]^+*)
MassTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{2,2},{True,False}}; (*S*S^+*)
MassTerm2=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;


VMass=(
	+m1*MassTerm1
	+mS1*MassTerm2
	);


\[Mu]ij=GradMass[VMass]//Simplify; 


QuarticTerm1=MassTerm1^2; (*[(\[Phi]\[Phi]^+)]^2*)
QuarticTerm2=MassTerm2^2; (* (S*S^+)^2*)
QuarticTerm3=MassTerm1*MassTerm2; (* (\[Phi]\[Phi]^+)(S*S^+)*)


VQuartic=(
	+\[Lambda]H*QuarticTerm1
	+\[Lambda]S*QuarticTerm2
	+\[Lambda]SH*QuarticTerm3
	);


\[Lambda]4=GradQuartic[VQuartic];


InputInv={{1,1,2},{False,False,True}}; 
YukawaDoublet1=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv][[1]]//Simplify;


Ysff=-yt*GradYukawa[YukawaDoublet1];


YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt>0}]];


(* ::Section:: *)
(*Dimensional Reduction*)


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PerformDRhard[]


BetaFunctions4D[]


PrintCouplings[]


PrintScalarMass["LO"]
PrintScalarMass["NLO"]


PrintDebyeMass["LO"]
PrintDebyeMass["NLO"]


PosScalar=PrintScalarRepPositions[];
PosVector=PrintGaugeRepPositions[];
PosFermion=PrintFermionRepPositions[];


Table[AnomDim4D["S",{i,i}],{i,PosScalar}]


Table[AnomDim4D["V",{i,i}],{i,PosVector}]


Table[AnomDim4D["F",{i,i}],{i,PosFermion}]
