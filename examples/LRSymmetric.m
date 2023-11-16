(* ::Package:: *)

(*Quit[];*)


SetDirectory[NotebookDirectory[]];
$LoadGroupMath=True;
<<../DRalgo.m


(* ::Chapter:: *)
(*Left-Right Symmetric Model*)


(*See 1811.06869 [hep-ph]*)


(* ::Section:: *)
(*Model*)


Group={"SU3","SU2","SU2","U1"};
RepAdjoint={{1,1},{2},{2},0};
HiggsBiDoublet={{{0,0},{1},{1},0},"C"};
HiggsLTriplet={{{0,0},{2},{0},1},"C"};
HiggsRTriplet={{{0,0},{0},{2},1},"C"};
RepScalar={HiggsBiDoublet,HiggsLTriplet,HiggsRTriplet};
CouplingName={gs,gL,gR,gBL};
(*For a LR symmetric model we should instead use*)
CouplingName={gs,gL,gL,gBL}; (*So gR=gL*)


Rep1={{{1,0},{1},{0},1/6},"L"}; (*QL*)
Rep2={{{1,0},{0},{1},1/6},"R"}; (*QR*)
Rep3={{{0,0},{1},{0},-1/2},"L"}; (*lL*)
Rep4={{{0,0},{0},{1},-1/2},"R"}; (*lR*)
RepFermion1Gen={Rep1,Rep2,Rep3,Rep4};


RepFermion3Gen={RepFermion1Gen,RepFermion1Gen,RepFermion1Gen}//Flatten[#,1]&;


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion3Gen,RepScalar];


Table[Tr[a . (b . c+c . b)],{a,gvff},{b,gvff},{c,gvff}]//SparseArray (*Anomaly check*)


InputInv={{1,1},{True,False}}; (*Tr \[Phi]\[Phi]^+*)
MassTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{1,1},{True,True}}; (*Tr \[Phi]\[Phi]^~*)
MassTerm2=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{1,1},{False,False}}; (*Tr \[Phi]^+\[Phi]^+~*)
MassTerm3=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{2,2},{True,False}}; (*Tr\[CapitalDelta]L \[CapitalDelta]L^+*)
MassTerm4=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{3,3},{True,False}}; (*Tr\[CapitalDelta]R \[CapitalDelta]R^+*)
MassTerm5=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;


VMass=(
	+\[Mu]1*MassTerm1
	+\[Mu]2(MassTerm2+MassTerm3)
	+\[Mu]3(MassTerm4+MassTerm5)
	);


\[Mu]ij=GradMass[VMass]//SparseArray;


(*First Quartics without ~*)


Quartic1=\[Lambda]H1 MassTerm1^2+\[Rho]1(MassTerm4^2+MassTerm5^2)+\[Rho]3 MassTerm4*MassTerm5;
InputInv={{1,1,2,2},{True,False,True,False}}; (*Terms with \[Phi]\[Phi]^+\[CapitalDelta]L\[CapitalDelta]L^+*)
QuarticTemp1=CreateInvariant[Group,RepScalar,InputInv]//Simplify;
InputInv={{1,1,3,3},{True,False,True,False}}; (*Terms with \[Phi]\[Phi]^+\[CapitalDelta]R\[CapitalDelta]R^+*)
QuarticTemp2=CreateInvariant[Group,RepScalar,InputInv]//Simplify;
InputInv={{1,1,2,3},{True,False,True,False}}; (*Terms with \[Phi]\[Phi]^+\[CapitalDelta]L\[CapitalDelta]R^+*)
QuarticTemp3=CreateInvariant[Group,RepScalar,InputInv]//Simplify;
InputInv={{1,1,3,2},{True,False,True,False}}; (*Terms with \[Phi]\[Phi]^+\[CapitalDelta]R\[CapitalDelta]L^+*)
QuarticTemp4=CreateInvariant[Group,RepScalar,InputInv]//Simplify;


Quartic2=(
	+\[Alpha]1(QuarticTemp1[[1]]+QuarticTemp2[[1]])
	+\[Alpha]3(QuarticTemp1[[2]]+QuarticTemp2[[2]])
	+\[Beta]1(QuarticTemp3[[1]]+QuarticTemp4[[1]])
	);


(*Now terms with Tr\[CapitalDelta]\[CapitalDelta]Tr \[CapitalDelta]^+(\[CapitalDelta]^+)~*)


InputInv={{2,2},{True,True}}; (*Tr \[CapitalDelta]L\[CapitalDelta]L*)
\[CapitalDelta]LSq=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{2,2},{False,False}}; (*Tr \[CapitalDelta]L^+\[CapitalDelta]L^+*)
\[CapitalDelta]LSqhc=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{3,3},{True,True}}; (*Tr \[CapitalDelta]R\[CapitalDelta]R*)
\[CapitalDelta]RSq=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{3,3},{False,False}}; (*Tr \[CapitalDelta]R^+\[CapitalDelta]R^+*)
\[CapitalDelta]RSqhc=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;


Quartic3=(
	+\[Rho]2*(\[CapitalDelta]LSq*\[CapitalDelta]LSqhc+\[CapitalDelta]RSq*\[CapitalDelta]RSqhc)
	+\[Rho]4*(\[CapitalDelta]LSq*\[CapitalDelta]RSqhc+\[CapitalDelta]LSqhc*\[CapitalDelta]RSq)
	);


(*Now terms with \[Phi]^~*)


Quartic4=(
	+\[Lambda]H2*(MassTerm2^2+MassTerm3^2)
	+\[Lambda]H3*MassTerm2*MassTerm3
	+\[Lambda]H4*MassTerm1*(MassTerm2+MassTerm3)
	);


Quartic5=\[Alpha]2(
	+(MassTerm4+MassTerm5)*MassTerm2
	+(MassTerm4+MassTerm5)*MassTerm3
	);


InputInv={{1,1,2,3},{True,True,True,False}}; (*Terms with \[Phi]\[Phi]^~\[CapitalDelta]L\[CapitalDelta]R^+*)
QuarticTemp3=CreateInvariant[Group,RepScalar,InputInv]//Simplify;
InputInv={{1,1,3,2},{True,True,True,False}}; (*Terms with \[Phi]\[Phi]^~\[CapitalDelta]R\[CapitalDelta]L^+*)
QuarticTemp4=CreateInvariant[Group,RepScalar,InputInv]//Simplify;
InputInv={{1,1,2,3},{False,False,True,False}}; (*Terms with \[Phi]\[Phi]^~\[CapitalDelta]L\[CapitalDelta]R^+*)
QuarticTemp5=CreateInvariant[Group,RepScalar,InputInv]//Simplify;
InputInv={{1,1,3,2},{False,False,True,False}}; (*Terms with \[Phi]\[Phi]^~\[CapitalDelta]R\[CapitalDelta]L^+*)
QuarticTemp6=CreateInvariant[Group,RepScalar,InputInv]//Simplify;


Quartic6=\[Beta]2(
	+QuarticTemp3[[1]]+QuarticTemp4[[1]]
	+QuarticTemp5[[1]]+QuarticTemp6[[1]]
	);


VQuartic=Quartic1+Quartic2+Quartic3+Quartic4+Quartic5+Quartic6;


\[Lambda]4=GradQuartic[VQuartic]//SparseArray;


InputInv={{1,1,2},{True,True,False}};  
YukawaTerm1=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv]//Simplify;


InputInv={{1,1,2},{False,True,False}};  
YukawaTerm2=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv]//Simplify;


Ysff=-GradYukawa[yt1*YukawaTerm1[[1]]+yt2*YukawaTerm2[[1]]];


YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt1>0,yt2>0}]];


(* ::Section:: *)
(*Dimensional Reduction*)


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False,Mode->2];
PerformDRhard[]


PrintCouplings[]


BetaFunctions4D[]


PrintTemporalScalarCouplings[]


PrintScalarMass["LO"]
PrintScalarMass["NLO"]


PrintDebyeMass["LO"]
PrintDebyeMass["NLO"]
