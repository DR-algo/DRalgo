(* ::Package:: *)

Quit[]


SetDirectory[NotebookDirectory[]];
$LoadGroupMath=True;
<<../DRalgo.m


(* ::Chapter:: *)
(*SU(2) + Doublet + Singlet Dark Sector*)


(* 
This is a DRalgo model file for the
SU(2) + doublet + singlet dark sector 'toy model' used in appendix C of
"Cosmological phase transitions at three loops: the final verdict on perturbation theory"
by Andreas Ekstedt, Philipp Schicho and Tuomas V. I. Tenkanen, [2405.xxxxx]. 

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
Group={"SU2"}; (* SU(2) gauge group *)
RepAdjoint={{2}}; (* gauge field in adjoint representation (rep) *)
Doublet={{{1}},"C"}; (* scalar in fundamental rep *)
Singlet={{{0}},"R"}; (* real scalar, in trivial rep *)
RepScalar={Doublet,Singlet}; 
CouplingName={g}; (* gauge coupling *)

(* no fermions *)
RepFermion = {}; 

(* allocating tensor structures: *)
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];

(* construct scalar group invariants for scalar potential: *)
InputInv={{1,1},{True,False}};
MassTerm1=CreateInvariant[Group,RepScalar,InputInv]//Simplify//FullSimplify;
InputInv={{2,2},{True,True}};
MassTerm2=CreateInvariant[Group,RepScalar,InputInv]//Simplify//FullSimplify;

(* mass terms *)
VMass=(m\[Phi]sq*MassTerm1+1/2 mssq*MassTerm2); (* here m\[Phi]sq = m\[Phi]^2 and  mssq = ms^2 *)
\[Mu]ij=GradMass[VMass[[1]]]//Simplify//SparseArray;

(* quartic terms *)
QuarticTerm1=MassTerm1[[1]]^2;
QuarticTerm2=MassTerm2[[1]]^2;
QuarticTerm3=MassTerm1[[1]]*MassTerm2[[1]];

VQuartic=(
	+\[Lambda]\[Phi]*QuarticTerm1 (* doublet self-interaction *)
	+\[Lambda]s/4*QuarticTerm2 (* singlet self-interaction *)
	+\[Lambda]m/2*QuarticTerm3 (* quartic portal interaction *)
	);
	
\[Lambda]4=GradQuartic[VQuartic];
	
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

(*
define group invariants for cubic operators:
*)
InputInv={{1,1,2},{True,False,True}};
CubicTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{2,2,2},{True,True,True}};
CubicTerm2=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;

(* and add cubic operators: *)
VCubic=(
	+\[Mu]m/2*CubicTerm1
	+\[Mu]3/3*CubicTerm2
	);

\[Lambda]3=GradCubic[VCubic];

(* add tadpole: *)
InputInv={{2},{True}};
TadpoleTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;

VTadpole=\[Mu]1*TadpoleTerm1;

\[Lambda]1=GradTadpole[VTadpole];


(* ::Section:: *)
(*Dimensional Reduction*)


(* import model defined above *)
ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False(*,Mode->3*)]; 
(* Mode -> 3 option is used for marginal operators if desired, in this case one also includes: *)
(*DefineDim6[\[Lambda]6]*)



(* ::Subsection:: *)
(*Hard to soft EFT construction*)


(*
'hard to soft' EFT construction
*)
PerformDRhard[] 


(*
the Lb logarithm is defnined as
Lb = Log[\[Mu]^2/T^2]+2*EulerGamma-2*Log[4\[Pi]];
*)
PrintConstants[]


(*
print the 3D couplings in the EFT in terms of
temperature T and parameters of the parent theory
*)
PrintCouplings[]


(*
couplings to temporal gauge field components
*)
PrintTemporalScalarCouplings[]


(*
Debye masses of temporal gauge field component to
2-loop order
(here \[Mu]sqSU2 == mD^2)
*)
PrintDebyeMass["LO"]
PrintDebyeMass["NLO"]


(*
print 3D masses in the EFT,i.e. thermal masses,
in terms of temperature T and parameters of the parent theory
*)
PrintScalarMass["LO"] (* one-loop thermal masses *)
PrintScalarMass["NLO"] (* two-loop correction to thermal masses *)


(*
singlet tadpole to two-loops:
*)
PrintTadpoles["LO"] 
PrintTadpoles["NLO"]


(* =
parent theory beta functions at one-loop:
*)
(*
Notation defined as
beta[g^2] = g^2 /. g^2\[Rule] -((43 g^4)/(48 \[Pi]^2))
etc.
*)
BetaFunctions4D[]


(*
unit operator pE, or hard mode contribution to symmetric phase pressure
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


(* ::Subsection:: *)
(*Soft to softer EFT construction*)


(*
'soft to softer' EFT construction,
by integrating out temporal gauge field components, and
the singlet (scalar list index 5)
*)
PerformDRsoft[{5}];


(*
scalar and gauge couplings in 'softer' 3D EFT:
*)
PrintCouplingsUS[]


(*
scalar 3d mass in 'softer' 3D EFT:
*)
PrintScalarMassUS["LO"]
PrintScalarMassUS["NLO"]


(*
the running of mass in 3d is governed by
*)
BetaFunctions3DUS[]


(*
soft mode constributions to unit operator (pM)
*)
PrintPressureUS["LO"] 
PrintPressureUS["NLO"]
(*
Currently, the "NNLO" part cannot be computed by DRalgo
*) 



