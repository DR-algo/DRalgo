(* ::Package:: *)

Quit[];


SetDirectory[NotebookDirectory[]];
$LoadGroupMath=True;
<<../DRalgo.m


(* ::Chapter:: *)
(*2 phi 4*)


(* ::Section::Closed:: *)
(*Model*)


Group={"SU3","SU2","U1"};
(*Group={};*)
RepAdjoint={{1,1},{2},0};
(*RepAdjoint={};*)
HiggsSinglet1={{{0,0},{0},0},"R"};
(*HiggsSinglet1={"R"};*)
HiggsSinglet2={{{0,0},{0},0},"R"};
(*HiggsSinglet2={"R"};*)
RepScalar={HiggsSinglet1,HiggsSinglet2};
CouplingName={g3,g2,g1};
(*CouplingName={};*)


(*Rep1={{{1,0},{1},1/6},"L"};
Rep2={{{1,0},{0},2/3},"R"};
Rep3={{{1,0},{0},-1/3},"R"};
Rep4={{{0,0},{1},-1/2},"L"};
Rep5={{{0,0},{0},-1},"R"};*)
(*RepFermion1Gen={Rep1,Rep2,Rep3,Rep4,Rep5};*)
RepFermion1Gen={};


RepFermion3Gen={RepFermion1Gen,RepFermion1Gen,RepFermion1Gen}//Flatten[#,1]&;


(* ::Text:: *)
(*The input for the gauge interactions toDRalgo are then given by*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion3Gen,RepScalar];


(* ::Text:: *)
(*The first element is the vector self - interaction matrix :*)


InputInv={{1,1},{True,True}}; (*S1^2*)
MassTerm1=CreateInvariant[Group,RepScalar,InputInv]//Simplify;
InputInv={{2,2},{True,True}}; (*S2^2*)
MassTerm2=CreateInvariant[Group,RepScalar,InputInv]//Simplify;


VMass=(
	+1/2*mS1*MassTerm1
	+1/2*mS2*MassTerm2
	);


\[Mu]ij=GradMass[VMass[[1]]]//Simplify; 


QuarticTerm1=MassTerm1[[1]]^2; (* S1^4*)
QuarticTerm2=MassTerm2[[1]]^2; (* S2^4*)
QuarticTerm3=MassTerm1[[1]]*MassTerm2[[1]] (* S1^2S2^2*)


VQuartic=(
	+\[Lambda]S1/4!*QuarticTerm1
	+\[Lambda]S2/4!*QuarticTerm2
	+\[Lambda]S12/4*QuarticTerm3
		);


\[Lambda]4=GradQuartic[VQuartic];


InputInv={{1},{True}}; (*S1*)
TadpoleTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{2},{True}}; (*S2*)
TadpoleTerm2=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;


VTadpole=(
	+\[Lambda]Tad1*TadpoleTerm1
	+\[Lambda]Tad2*TadpoleTerm2
	);


\[Lambda]1=GradTadpole[VTadpole];


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
PrintScalarMass["NLO"]


PrintPressure["LO"]
PrintPressure["NLO"]
PrintPressure["NNLO"]


(* ::Subsubsection:: *)
(*2 Scalars active:*)


PerformDRsoft[{}]


PrintCouplingsUS[]


PrintScalarMassUS["LO"]//Simplify


DefineTensorsUS[]


DefineTensorsUS[]
\[CurlyPhi]vev={\[CurlyPhi],0}//SparseArray;
DefineVEVS[\[CurlyPhi]vev];


FieldMasses=PrintTensorsVEV[]
CalculatePotentialUS[];


V3dLO=PrintEffectivePotential["LO"]
V3dNLO=PrintEffectivePotential["NLO"]
V3dNNLO=PrintEffectivePotential["NNLO"]


(* ::Subsubsection:: *)
(*1 Scalar active:*)


PerformDRsoft[{2}]
PrintCouplingsUS[]


PrintScalarMassUS["LO"]//Simplify
PrintScalarMassUS["NLO"]//Simplify


DefineTensorsUS[]


\[CurlyPhi]vev={\[CurlyPhi]}//SparseArray;
DefineVEVS[\[CurlyPhi]vev];
FieldMasses=PrintTensorsVEV[];


\[Lambda]4US=PrintTensorUSDRalgo[][[1]]
gvssUS=PrintTensorUSDRalgo[][[2]]
\[Mu]ijUS=PrintTensorUSDRalgo[][[3]]


CalculatePotentialUS[];


V3dLO=PrintEffectivePotential["LO"]
V3dNLO=PrintEffectivePotential["NLO"]
V3dNNLO=PrintEffectivePotential["NNLO"]


(* ::Subsection:: *)
(*Plotting the potential*)


(*For now only use US potential with two scalars active*)


Rlrest={\[CurlyPhi]3d-> \[Phi] /Sqrt[T], \[CurlyPhi]0-> 300/Sqrt[T],Lb-> Log[\[Mu]^2/T^2]+2EulerGamma-2 Log[4\[Pi]]};
V3d=V3dLO+V3dNLO+V3dNNLO;
V3dprep=V3d/.{mS1-> mS13d,\[CurlyPhi]-> \[CurlyPhi]3d,mS2-> mS23d,\[Lambda]S1-> \[Lambda]S13d,\[Lambda]S2-> \[Lambda]S23d,\[Lambda]S12-> \[Lambda]S123d};
Vfull=V3dprep/.PrintCouplings[];
Vfull=Vfull/.Thread[PrintScalarMass["LO"][[;;,1]]->PrintScalarMass["LO"][[;;,2]]+PrintScalarMass["NLO"][[;;,2]]];
Vfull=Vfull/.Rlrest;


(*Potential*)
(*looks strange since two minimia where origin is a local maximum*)
Vpar=T Vfull/.{\[Lambda]S1-> 6*0.048828125`,\[Lambda]S12-> 4.4,\[Lambda]S2-> 0.1*6,mS2-> 8000,mS1-> -7812.5};
Plot[Vpar/.{T-> 295,\[Mu]-> 500,\[Mu]3-> 500},{\[Phi],0,350}]
