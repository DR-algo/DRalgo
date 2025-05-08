(* ::Package:: *)

(*Quit[];*)


SetDirectory[NotebookDirectory[]];
DRalgo`$LoadGroupMath=True;
<<../DRalgo.m


(* ::Chapter:: *)
(*2HDM-Two Higgs doublet model*)


(*See 1106.0034 [hep-ph] for a review*)


(* ::Section::Closed:: *)
(*Model*)


Group={"SU3","SU2","U1"};
RepAdjoint={{1,1},{2},0};
HiggsDoublet1={{{0,0},{1},1/2},"C"};
HiggsDoublet2={{{0,0},{1},1/2},"C"};
RepScalar={HiggsDoublet1,HiggsDoublet2};
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


InputInv={{1,1},{True,False}};
MassTerm1=CreateInvariant[Group,RepScalar,InputInv]//Simplify//FullSimplify;
InputInv={{2,2},{True,False}};
MassTerm2=CreateInvariant[Group,RepScalar,InputInv]//Simplify//FullSimplify;
InputInv={{1,2},{True,False}};
MassTerm3=CreateInvariant[Group,RepScalar,InputInv]//Simplify//FullSimplify;
InputInv={{2,1},{True,False}};
MassTerm4=CreateInvariant[Group,RepScalar,InputInv]//Simplify//FullSimplify;


VMass=(
	+m1*MassTerm1
	+m2*MassTerm2
	-(m12R+I*m12I)(MassTerm3)
	-(m12R-I*m12I)(MassTerm4)
	);


\[Mu]ij=GradMass[VMass[[1]]]//Simplify//SparseArray;


QuarticTerm1=MassTerm1[[1]]^2;
QuarticTerm2=MassTerm2[[1]]^2;
QuarticTerm3=MassTerm1[[1]]*MassTerm2[[1]];
QuarticTerm4=MassTerm3[[1]]*MassTerm4[[1]];
QuarticTerm5=(MassTerm3[[1]]^2+MassTerm4[[1]]^2)//Simplify;
QuarticTerm6=MassTerm1[[1]]MassTerm3[[1]]+MassTerm1[[1]]MassTerm4[[1]];
QuarticTerm7=MassTerm2[[1]]MassTerm3[[1]]+MassTerm2[[1]]MassTerm4[[1]];


VQuartic=(
	+\[Lambda]1H*QuarticTerm1
	+\[Lambda]2H*QuarticTerm2
	+\[Lambda]3H*QuarticTerm3
	+\[Lambda]4H*QuarticTerm4
	+\[Lambda]5H/2*QuarticTerm5
	+\[Lambda]6H*QuarticTerm6
	+\[Lambda]7H*QuarticTerm7
	);


\[Lambda]4=GradQuartic[VQuartic];


InputInv={{1,1,2},{False,False,True}}; 
YukawaDoublet1=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv]//Simplify;


InputInv={{2,1,2},{False,False,True}}; 
YukawaDoublet2=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv]//Simplify;


Ysff=-GradYukawa[yt1*YukawaDoublet1[[1]]+yt2*YukawaDoublet2[[1]]];


YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt1>0,yt2>0}]];


(* ::Section:: *)
(*Dimensional Reduction*)


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PerformDRhard[]


BetaFunctions4D[]


PrintScalarMass["LO"]
PrintScalarMass["NLO"]


PrintDebyeMass["LO"]
PrintDebyeMass["NLO"]


PrintCouplings[]
PrintTemporalScalarCouplings[]


(* ::Text:: *)
(*Two active doublets:*)


PerformDRsoft[{}];


PrintCouplingsUS[]//Simplify


PrintScalarMassUS["LO"]
PrintScalarMassUS["NLO"]


BetaFunctions3DUS[]


(* ::Text:: *)
(*One active doublet :*)


PerformDRsoft[{5,6,7,8}];


PrintCouplingsUS[]


PrintScalarMassUS["LO"]
PrintScalarMassUS["NLO"]


BetaFunctions3DUS[]


PrintPressureUS["LO"]
PrintPressureUS["NLO"]


(* ::Chapter:: *)
(*Saving the model*)


result={};
AppendTo[result,Row[{
	TexFor["DRDRDRDRDRDRDRDRDRDRDRDRDRDR "],
	TexFor["DRalgo"],
	TexFor[" DRDRDRDRDRDRDRDRDRDRDRDRDRDRD"]}]];
AppendTo[result,Row[{"Model: "//TexFor,"Two-Higgs doublet model. See hep-ph:1106.0034 for further details"//TexFor}]];
AppendTo[result,Row[{"Version: "//TexFor,"1.0 beta (16-05-2022)"//TexFor}]];
AppendTo[result,Row[{"Authors: "//TexFor,"Andreas Ekstedt, Philipp Schicho, Tuomas V.I. Tenkanen"//TexFor}]];
AppendTo[result,Row[{"Reference: "//TexFor,"2205.08815 [hep-ph]"//TexFor}]];
AppendTo[result,Row[{"Repository link: "//TexFor,
	Hyperlink[Mouseover[TexFor["github.com/DR-algo/DRalgo"],Style["github.com/DR-algo/DRalgo",Bold]],
	"https://github.com/DR-algo/DRalgo"]}]];
AppendTo[result,Style["DRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRD",{GrayLevel[0.3]}]];


ModelInfo=result;
SaveModelDRalgo[ModelInfo,"2hdm.txt"]


(* ::Chapter::Closed:: *)
(*Loading the model*)


{Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=LoadModelDRalgo["2hdm.txt"];
ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PerformDRhard[]
