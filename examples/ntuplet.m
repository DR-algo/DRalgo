(* ::Package:: *)

Quit[];


SetDirectory[NotebookDirectory[]];
DRalgo`DRalgo`$GroupMathMultipleModels=True; (*Put this if you want to create multiple model-files with the same kernel*)
<<DRalgo`DRalgo`


(* ::Chapter:: *)
(*SM+ 5tuplet*)


(*See 1812.07829, 1802.10500 [hep-ph], *)


(* ::Section:: *)
(*Model*)


Group={"SU3","SU2","U1"};
RepAdjoint={{1,1},{2},0};
HiggsDoublet1={{{0,0},{1},Y\[Phi]/2},"C"};
(*ndim, triplet (ndim=3), ...*)
ndim=7;
Scalarntuplet={{{0,0},{(ndim-1)},0},"R"};
RepScalar={HiggsDoublet1,Scalarntuplet};
CouplingName={g3,gw,g1};


(* Show different dimensional representations of SU2 *)
su2Reps = RepsUpToDimN[SU2,9];
Grid[Prepend[{#,RepName[SU2,#]}&/@ su2Reps,{"Dynkin coefficients","Name"}],
Frame->All,FrameStyle->LightGray]


Rep1={{{1,0},{1},1/6},"L"};
Rep2={{{1,0},{0},2/3},"R"};
Rep3={{{1,0},{0},-1/3},"R"};
Rep4={{{0,0},{1},-1/2},"L"};
Rep5={{{0,0},{0},-1},"R"};
RepFermion1Gen={Rep1,Rep2,Rep3,Rep4,Rep5};


RepFermion3Gen={RepFermion1Gen}//Flatten[#,1]&;


(* ::Text:: *)
(*The input for the gauge interactions toDRalgo are then given by*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion3Gen,RepScalar];


(* ::Text:: *)
(*The first element is the vector self - interaction matrix :*)


InputInv={{1,1},{True,False}}; (*\[Phi]\[Phi]^+*)
MassTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;
(* two real ntuplets *)
InputInv={{2,2},{True,True}}; (*\[CapitalSigma]^a\[CapitalSigma]^a*)
MassTerm2=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;


VMass=(
	+mHsq*MassTerm1
	+mNTsq/2*MassTerm2
	);


\[Mu]ij=GradMass[VMass]//Simplify//SparseArray;


(*Higgs self interaction*)
InputInv={{1,1,1,1},{True,False,True,False}};
QuarticTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]];
(*Portal*)
InputInv={{1,1,2,2},{True,False,True,True}};
QuarticTerm2=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;

(*nTuplet self interaction*)
InputInv={{2,2,2,2},{True,True,True,True}};
QuarticTerm3=CreateInvariant[Group,RepScalar,InputInv]//DeleteCases[#,0]&;


VQuartic=(
	+\[Lambda]*QuarticTerm1
	+a2/2*QuarticTerm2
	+QuarticTerm3 . Table[ToExpression["lambda"<>ToString[i]],{i,0,2*(Length[QuarticTerm3]-1),2}]
	);


\[Lambda]4=GradQuartic[VQuartic];


InputInv={{1,1,2},{False,False,True}}; 
YukawaDoublet1=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv]//Simplify;


Ysff=-GradYukawa[yt*YukawaDoublet1[[1]]];


YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt>0}]];


(* ::Section:: *)
(*Dimensional Reduction*)


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PosFermion=PrintFermionRepPositions[];
FermionMat=Table[{nF,i},{i,PosFermion}];
DefineNF[FermionMat]
PerformDRhard[]


BetaFunctions4D[]


PrintCouplings[]


(*Higgs thermal mass*)
mHsq3d/.PrintScalarMass["LO"]//Expand
mHsq3d/.PrintScalarMass["NLO"]//Expand;


(*ntuplet thermal mass*)
mNTsq3d/.PrintScalarMass["LO"]//Expand
mNTsq3d/.PrintScalarMass["NLO"]//Expand


(* ::Subsubsection:: *)
(*Check Debye Mass*)


dY=1/3*j*(j+1)*(2j+1);


PrintDebyeMass["LO"]
(* SU2 Debye Mass can be written in terms of j *)
mDsq = gw^2*T^2*((4+Nd+dY)/6+nF/3)/.{Nd->1,j->(ndim-1)/2}//Simplify

PrintDebyeMass["NLO"]


(* ::Subsection:: *)
(*Reduce to US theory*)


PerformDRsoft[Range[5,5+ndim-1]]


PrintCouplingsUS[]


PrintScalarMassUS["LO"]
PrintScalarMassUS["NLO"]


BetaFunctions3DUS[]


PrintTadpolesUS["LO"]


(* ::Section:: *)
(*Effective potential*)


(*
Computations are done in the soft theory
*)
UseUltraSoftTheory[];


\[CurlyPhi]vev={\[CurlyPhi],0,0,0}//SparseArray; 
DefineVEVS[\[CurlyPhi]vev];


ScalarMass=PrintTensorsVEV[1]//Normal
VectorMass=PrintTensorsVEV[2]//Normal


VectorMass=PrintTensorsVEV[2]//Normal;
VectorEigenvectors=FullSimplify[
    Transpose[Normalize/@Eigenvectors[VectorMass[[11;;12,11;;12]]]],
Assumptions->{gw>0,gY>0,\[CurlyPhi]>0}]; DVRot={{IdentityMatrix[10],0},{0,VectorEigenvectors}}//ArrayFlatten; DSRot=IdentityMatrix[4];
RotateTensorsUSPostVEV[DSRot,DVRot];


PrintTensorsVEV[2]//Normal//Simplify//MatrixForm


CalculatePotential[];


PrintEffectivePotential["LO"]
PrintEffectivePotential["NLO"]
PrintEffectivePotential["NNLO"]



