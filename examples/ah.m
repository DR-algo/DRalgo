(* ::Package:: *)

(*Quit[];*)


SetDirectory[NotebookDirectory[]];
$LoadGroupMath=True;
<<../DRalgo.m


(* ::Chapter:: *)
(*Abelian Higgs*)


(* ::Section::Closed:: *)
(*Model*)


Group={"U1"};
CouplingName={g1};
RepAdjoint={0};
Higgs={{Y\[Phi]},"C"}; (* Y\[Phi] = 1 *)
RepScalar={Higgs};


RepFermion={};


(* ::Text:: *)
(*The first element is the vector self-interaction matrix:*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];


InputInv={{1,1},{True,False}}; (*This specifies that we want a \[Phi]^+\[Phi] term*)
MassTerm1=CreateInvariant[Group,RepScalar,InputInv]//Simplify//FullSimplify;


VMass=msq*MassTerm1[[1]];(*This is the \[Phi]^+\[Phi] term written in component form*)


\[Mu]ij=GradMass[VMass]//Simplify//SparseArray;


QuarticTerm1=MassTerm1[[1]]^2; (*Because MassTerm1=\[Phi]^+\[Phi], we can write (\[Phi]^+\[Phi])^2=MassTerm1^2*)


VQuartic=\[Lambda]*QuarticTerm1;


\[Lambda]4=GradQuartic[VQuartic];


(* ::Section:: *)
(*Dimensional Reduction*)


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PerformDRhard[]
PrintCouplings[]


PrintConstants[]


PrintScalarMass["LO"]


PrintScalarMass["NLO"]



PrintDebyeMass["LO"]
PrintDebyeMass["NLO"]


PrintTensorDRalgo[];


PrintTemporalScalarCouplings[]//Simplify


PrintPressure["LO"]
PrintPressure["NLO"]
PrintPressure["NNLO"]


BetaFunctions4D[]


posScalar=PrintScalarRepPositions[];
posVector=PrintGaugeRepPositions[];


AnomDim4D["S",{posScalar[[1]],posScalar[[1]]}]
AnomDim4D["V",{posVector[[1]],posVector[[1]]}]


(* ::Text:: *)
(*Integrating out temporal scalars:*)


PerformDRsoft[{}];


PrintCouplingsUS[]


PrintScalarMassUS["LO"]
PrintScalarMassUS["NLO"]


BetaFunctions3DUS[]


PrintPressureUS["LO"]
PrintPressureUS["NLO"]


(* ::Text:: *)
(*Effective potential:*)


DefineNewTensorsUS[\[Mu]ij,\[Lambda]4,\[Lambda]3,gvss,gvvv];
\[Phi]VeV={0,\[Phi]}//SparseArray;
DefineVEVS[\[Phi]VeV];
FieldMasses = PrintTensorsVEV[];


CalculatePotentialUS[];


PrintEffectivePotential["LO"]
PrintEffectivePotential["NLO"]
PrintEffectivePotential["NNLO"]


(* ::Section:: *)
(*Saving the model*)


(*result={};
AppendTo[result,Row[{
	TexFor["DRDRDRDRDRDRDRDRDRDRDRDRDRDR "],
	TexFor["DRalgo"],
	TexFor[" DRDRDRDRDRDRDRDRDRDRDRDRDRDRD"]}]];
	AppendTo[result,Row[{"Model: "//TexFor,"Abelian Higgs. See hep-ph:9709418 for further details"//TexFor}]];
AppendTo[result,Row[{"Version: "//TexFor,"1.0 beta (16-05-2022)"//TexFor}]];
AppendTo[result,Row[{"Authors: "//TexFor,"Andreas Ekstedt, Philipp Schicho, Tuomas V.I. Tenkanen"//TexFor}]];
AppendTo[result,Row[{"Reference: "//TexFor,"2205.xxxxx [hep-ph]"//TexFor}]];
AppendTo[result,Row[{"Repository link: "//TexFor,
	Hyperlink[Mouseover[TexFor["github.com/DR-algo/DRalgo"],Style["github.com/DR-algo/DRalgo",Bold]],
	"https://github.com/DR-algo/DRalgo"]}]];
AppendTo[result,Style["DRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRD",{GrayLevel[0.3]}]];*)


(*ModelInfo=result;*)


(*SaveModelDRalgo[ModelInfo,"ah.txt"]*)


(* ::Section:: *)
(*Loading the model*)


(*{Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=LoadModelDRalgo["ah.txt"];*)
