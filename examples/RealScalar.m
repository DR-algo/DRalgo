(* ::Package:: *)

(*Quit[];*)


SetDirectory[NotebookDirectory[]];
$LoadGroupMath=True;
<<../DRalgo.m


(* ::Title:: *)
(*A single real-scalar field; model created by Oliver Gould (28-06-2022)*)


(*See 2101.05528 [hep-th] for further details and for independent calculations*)


(* ::Section::Closed:: *)
(*Model*)


Group={"U1"};
CouplingName={g};
RepAdjoint={0};
RepScalar={{{0},"R"}};
RepFermion={};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJC,\[Mu]IJ,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];


(* ::Text:: *)
(*Tadpoles*)


(* \[Sigma] \[Phi] *)
InputInv={{1},{True}};
LinearTerm=CreateInvariant[Group,RepScalar,InputInv][[1]];
VLinear=\[Sigma] LinearTerm;
\[Lambda]1=GradTadpole[VLinear];


(* ::Text:: *)
(*Scalar-Mass terms*)


(* 1/2m^2\[Phi]^2 *)
InputInv={{1,1},{True,True}}; 
MassTerm=CreateInvariant[Group,RepScalar,InputInv][[1]];
VMass=msq/2 MassTerm;
\[Mu]ij=GradMass[VMass]//SparseArray;


(* ::Text:: *)
(*Scalar - Cubic terms*)


(* 1/6\[Gamma]\[Phi]^3 *)
InputInv={{1,1,1},{True,True,True}};
CubicTerm=CreateInvariant[Group,RepScalar,InputInv][[1]];
VCubic=\[Gamma]/6 CubicTerm;
\[Lambda]3=GradCubic[VCubic];


(* ::Text:: *)
(*Scalar - Quartic terms*)


(* 1/24\[Lambda]\[Phi]^4 *)
InputInv={{1,1,1,1},{True,True,True,True}};
QuarticTerm=CreateInvariant[Group,RepScalar,InputInv][[1]];
VQuartic=\[Lambda]/24 MassTerm^2;
\[Lambda]4=GradQuartic[VQuartic];


(* ::Section::Closed:: *)
(*Dimensional Reduction*)


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PerformDRhard[]


PrintCouplings[]


PrintScalarMass["LO"]
PrintScalarMass["NLO"]


BetaFunctions4D[]


PrintTadpoles["LO"]//Expand
PrintTadpoles["NLO"]


PrintPressure["LO"]
PrintPressure["NLO"]
PrintPressure["NNLO"]



