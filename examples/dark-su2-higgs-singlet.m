(* ::Package:: *)

(*Quit[];*)


SetDirectory[NotebookDirectory[]];
$LoadGroupMath=True;
<<../DRalgo.m


(* ::Chapter:: *)
(*SU(2) + Doublet + Singlet Dark Sector*)


(* 
This is a DRalgo model file for the
SU(2) + doublet + singlet dark sector 'toy model' used in appendix C of
"Cosmological phase transitions at three loops: the final verdict on perturbation theory"
by Andreas Ekstedt, Philipp Schicho and Tuomas V. I. Tenkanen, [2405.18349]. 

The model is a simplified toy setup, where only
the SU(2) gauge fields and doublet and singlet scalars are added.
The Standard Model field content (U(1) and SU(3) gauge sectors, and fermion sector)
is not included. For effects of these sectors, 
see e.g. DRalgo model file sm.m. 
*)


(* ::Section::Closed:: *)
(*Model*)


(*
	Running this cell defines the model.
	For more documentation, see the DRalgo manual [2205.08815].
*)
(* SU(2) gauge group *)
Group={"SU2"};
(* gauge field in adjoint representation (rep) *)
RepAdjoint={{2}};
(* scalar in fundamental rep *)
Doublet={{{1}},"C"};
(* real scalar, in trivial rep *)
Singlet={{{0}},"R"};
RepScalar={Doublet,Singlet};
(* gauge coupling *)
CouplingName={g};


(* no fermions *)
RepFermion = {}; 


(* ::Text:: *)
(*The input for the gauge interactions to DRalgo are then given by*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];


(* ::Text:: *)
(*Define scalar group invariants for scalar potential:*)


InputInv={{1,1},{True,False}};
MassTerm1=CreateInvariant[Group,RepScalar,InputInv]//Simplify//FullSimplify;
InputInv={{2,2},{True,True}};
MassTerm2=CreateInvariant[Group,RepScalar,InputInv]//Simplify//FullSimplify;


(* mass terms *)
(* here m\[Phi]sq = m\[Phi]^2 and  mssq = ms^2 *)
VMass=(m\[Phi]sq*MassTerm1+1/2 mssq*MassTerm2); 
\[Mu]ij=GradMass[VMass[[1]]]//Simplify//SparseArray;


(* ::Text:: *)
(*Add quartic operators:*)


QuarticTerm1=MassTerm1[[1]]^2;
QuarticTerm2=MassTerm2[[1]]^2;
QuarticTerm3=MassTerm1[[1]]*MassTerm2[[1]];


VQuartic=(
	+\[Lambda]\[Phi]*QuarticTerm1 (* doublet self-interaction *)
	+\[Lambda]s/4*QuarticTerm2 (* singlet self-interaction *)
	+\[Lambda]m/2*QuarticTerm3 (* quartic portal interaction *)
	);


\[Lambda]4=GradQuartic[VQuartic];


(* ::Text:: *)
(*Add cubic operators:*)


InputInv={{1,1,2},{True,False,True}};
CubicTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{2,2,2},{True,True,True}};
CubicTerm2=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;


VCubic=(
	+\[Mu]m/2*CubicTerm1
	+\[Mu]3/3*CubicTerm2
	);


\[Lambda]3=GradCubic[VCubic];


(* ::Text:: *)
(*Add tadpole operator:*)


InputInv={{2},{True}};
TadpoleTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;


VTadpole=\[Mu]1*TadpoleTerm1;


\[Lambda]1=GradTadpole[VTadpole];


(* ::Text:: *)
(*Add some dimension six operators:*)


(*
The leading marginal operators could be included by defining following Dim6Term,
and by using below the ImportModelDRalgo with option Mode -> 3.
Here, we do not include them for simplicity.
*)	
(*
Dim6Term=(
	+ c1 MassTerm1[[1]]^3
	+ c3 MassTerm1[[1]]^2 MassTerm2[[1]]
	+ c3 MassTerm1[[1]] MassTerm2[[1]]^2
	+ c4 MassTerm2[[1]]^3);
\[Lambda]6=GradSextic[Dim6Term];
*)


(* ::Section:: *)
(*Dimensional Reduction*)


(* ::Subsection::Closed:: *)
(*Model Import*)


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False]; 


(* ::Text:: *)
(* Mode -> 3 option is used for marginal operators if desired, in this case one also includes:*)


(*ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False,Mode->3];*)
(*DefineDim6[\[Lambda]6]*)


(* ::Subsection::Closed:: *)
(*Hard-to-soft EFT construction*)


PerformDRhard[] 


(*
the Lb logarithm is defnined as
Lb = Log[\[Mu]^2/T^2]+2*EulerGamma-2*Log[4\[Pi]];
*)
PrintConstants[]


(* ::Subsubsection::Closed:: *)
(*Couplings  in  the  3D EFT*)


(*
The 3D couplings are functions of
temperature T and parameters of the parent theory
*)
PrintCouplings[]


(* ::Subsubsection::Closed:: *)
(*Couplings to temporal gauge field components:*)


PrintTemporalScalarCouplings[]


(* ::Subsubsection::Closed:: *)
(*Debye  masses  of  temporal  gauge  field  component  to 2-loop  order*)


(*
(here \[Mu]sqSU2 == mD^2)
*)
PrintDebyeMass["LO"]
PrintDebyeMass["NLO"]


(* ::Subsubsection::Closed:: *)
(*Print 3D masses in the EFT*)


(*
These are the thermal masses,
in terms of temperature T and parameters of the parent theory
*)
PrintScalarMass["LO"] (* one-loop thermal masses *)
PrintScalarMass["NLO"] (* two-loop correction to thermal masses *)


(* ::Subsubsection::Closed:: *)
(*Singlet tadpole to 2-loops:*)


PrintTadpoles["LO"] 
PrintTadpoles["NLO"]


(* ::Subsubsection::Closed:: *)
(*1-loop parent  theory  beta  functions:*)


(*
Notation defined as
beta[g^2] = g^2 /. g^2\[Rule] -((43 g^4)/(48 \[Pi]^2))
etc.
*)
BetaFunctions4D[]


(* ::Subsubsection::Closed:: *)
(*Unit operator pE:*)


(*
this is the hard mode contribution to symmetric phase pressure
*)
PrintPressure["LO"]
PrintPressure["NLO"] 
PrintPressure["NNLO"]


(*
indices 1 to 4 are doublet components, and index 5 a singlet
*)
PosScalar=PrintScalarRepPositions[]
(*PosVector=PrintGaugeRepPositions[]
PosFermion=PrintFermionRepPositions[]*)


(* ::Subsection::Closed:: *)
(*Soft-to-softer EFT construction*)


(*
by integrating out temporal gauge field components, and
the singlet (scalar list index 5)
*)
PerformDRsoft[{5}];


(* ::Subsubsection::Closed:: *)
(*Scalar and gauge couplings in softer 3D EFT:*)


PrintCouplingsUS[]


(* ::Subsubsection::Closed:: *)
(*Scalar 3D mass in softer 3D EFT:*)


PrintScalarMassUS["LO"]
PrintScalarMassUS["NLO"]


(* ::Subsubsection::Closed:: *)
(*Softer 3D EFT beta functions:*)


BetaFunctions3DUS[]


(* ::Subsubsection::Closed:: *)
(*The  unit  operator  pM:*)


(*
soft mode constributions to unit operator (pM)
*)
PrintPressureUS["LO"] 
PrintPressureUS["NLO"]

(*
Currently, the "NNLO" part cannot be computed by DRalgo
*) 
PrintPressureUS["NNLO"]



