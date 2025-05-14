(* ::Package:: *)

(*Quit[];*)


SetDirectory[NotebookDirectory[]];
DRalgo`DRalgo`$LoadGroupMath=True;
<<../Kernel/DRalgo.m


(* ::Title:: *)
(*Wess-Zumino Model*)


(*See 0303260 [hep-th]*)


(* ::Section::Closed:: *)
(*Model*)


Group={"U1"};
RepAdjoint={0};
HiggsA={{0},"R"};
HiggsB={{0},"R"};
RepScalar={HiggsA,HiggsB};
CouplingName={gY};


Rep1={{0},"L"};
RepFermion={Rep1};


(* ::Text:: *)
(*The input for the gauge interactions toDRalgo are then given by*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];


(* ::Text:: *)
(*The first element is the vector self - interaction matrix :*)


InputInv={{1,1},{True,True}}; (*A^2*)
MassTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;
InputInv={{2,2},{True,True}}; (*B^2*)
MassTerm2=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;


VMass=1/2(mA*MassTerm1+mB*MassTerm2);


\[Mu]ij=GradMass[VMass]//Simplify//SparseArray;


QuarticTerm1=MassTerm1^2; (*A^4*)
QuarticTerm2=MassTerm2^2; (*B^4*)
QuarticTerm3=MassTerm1*MassTerm2; (*A^2B^2*)


VQuartic=1/4*\[Lambda](QuarticTerm1+QuarticTerm2+2*QuarticTerm3)//Simplify;


\[Lambda]4=GradQuartic[VQuartic];


InputInv={{1,1,1},{True,True,True}}; (*A^3*)
CubicTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;
InputInv={{1,2,2},{True,True,True}}; (*A^3*)
CubicTerm2=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;
VCubic=\[Lambda]c(CubicTerm1+CubicTerm2);


\[Lambda]3=GradCubic[VCubic];


InputInv={{1},{True}}; (*A*)
TadPoleTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;
InputInv={{2},{True}}; (*B*)
TadPoleTerm2=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;
VTadpole=\[Lambda]Tad1*TadPoleTerm1+\[Lambda]Tad2*TadPoleTerm2;


\[Lambda]1=GradTadpole[VTadpole];


InputInv={{1,1},{True,False}}; (*\[CapitalPsi]^+ \[CapitalPsi]*)
MassTerm1=CreateInvariantFermion[Group,RepFermion,InputInv][[1]]//Simplify//FullSimplify;


\[Mu]IJ=mF/2*GradMassFermion[MassTerm1];
\[Mu]IJC=SparseArray[Simplify[Conjugate[\[Mu]IJ]//Normal,Assumptions->{mF>0}]];


InputInv={{1,1,1},{True,True,True}}; 
YukawaDoublet1=CreateInvariantYukawa[Group,RepScalar,RepFermion,InputInv][[1]]//Simplify;
InputInv={{2,1,1},{True,True,True}}; 
YukawaDoublet2=CreateInvariantYukawa[Group,RepScalar,RepFermion,InputInv][[1]]//Simplify;


Ysff= yt/Sqrt[2]*GradYukawa[YukawaDoublet1+I*YukawaDoublet2 ];


YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt>0}]];


(* ::Section:: *)
(*Dimensional Reduction*)


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PerformDRhard[]


PrintCouplings[]


PrintScalarMass["LO"]
PrintScalarMass["NLO"]


BetaFunctions4D[]


PrintTadpoles["LO"]
PrintTadpoles["NLO"]


PrintPressure["LO"]
PrintPressure["NLO"]
PrintPressure["NNLO"]

