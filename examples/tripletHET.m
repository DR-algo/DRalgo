(* ::Package:: *)

(*Quit[];*)


SetDirectory[NotebookDirectory[]];
$LoadGroupMath=True;
<<../DRalgo.m


(* ::Chapter:: *)
(*SM+Real SU(2)Triplet*)


(*See 1802.10500 [hep-ph], *)


(* ::Section:: *)
(*Model*)


Group={"SU3","SU2","U1"};
RepAdjoint={{1,1},{2},0};
HiggsDoublet1={{{0,0},{1},1/2},"C"};
HiggsTriplet1={{{0,0},{2},0},"R"};
RepScalar={HiggsDoublet1,HiggsTriplet1};
CouplingName={g3,gw,g1};


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


InputInv={{1,1},{True,False}}; (*\[Phi]\[Phi]^+*)
MassTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;
InputInv={{2,2},{True,True}}; (*\[CapitalSigma]^a\[CapitalSigma]^a*)
MassTerm2=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;


VMass=(
	+mH*MassTerm1
	+m\[CapitalSigma]/2*MassTerm2
	);


\[Mu]ij=GradMass[VMass]//Simplify//SparseArray;


QuarticTerm1=MassTerm1^2; (*[(\[Phi]\[Phi]^+)]^2*)
QuarticTerm2=MassTerm2^2; (*[\[CapitalSigma]^a\[CapitalSigma]^a]^2*)
QuarticTerm3=MassTerm1*MassTerm2; (*[\[Phi]\[Phi]^+][\[CapitalSigma]^a\[CapitalSigma]^a]*)


VQuartic=(
	+\[Lambda]*QuarticTerm1
	+b4/4*QuarticTerm2
	+a2/2*QuarticTerm3
	);


\[Lambda]4=GradQuartic[VQuartic];


InputInv={{1,1,2},{False,False,True}}; 
YukawaDoublet1=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv]//Simplify;


Ysff=-GradYukawa[yt*YukawaDoublet1[[1]]];


YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt>0}]];


(* ::Section:: *)
(*Dimensional Reduction*)


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PerformDRhard[]


PrintCouplings[]
BetaFunctions4D[]


PrintScalarMass["LO"]
PrintScalarMass["NLO"]


PrintDebyeMass["LO"]
PrintDebyeMass["NLO"]


PerformDRsoft[{}];


PrintCouplingsUS[]


(* ::Section:: *)
(*Integrate out triplet in a Higgs Effective-type Theory (HET)*)


(*
Computations are done in the ultrasoft theory
*)
UseUltraSoftTheory[];


(*
In HEFT-like theories the vev of the Higgs or any other scalar
contributes substationally to the heavy particles mass.
Thus, the Higgs-vev contribution must be kept when integrating out the heavy particle.
*)


\[CurlyPhi]vev={0,0,0,\[CurlyPhi],0,0,0}//SparseArray;
DefineVEVS[\[CurlyPhi]vev];


PrintTensorsVEV[1]//Normal


(* ::Chapter:: *)
(*Integrating out SU (2) bosons*)


PrepareHET[{},{9,10,11}]


CalculatePotentialHET[]


(*
When integrating out SU(2) bosons the LO potential is the same,
but we get a \[Phi]^3 term at NLO; this term is responsible for
a barrier in models with radiative symmetry breaking
*)


PrintActionHET["LO"]


PrintActionHET["NLO"]


(*
The contribution at NNLO is
*)
PrintActionHET["NNLO"]


(*
The SU(2) gauge bosons also change the scalar fields' kinetic term,
printed below is the Z factor
*)
(*
that is PrintScalarKineticHET[][[i,j]]=\[Delta]Z^ij; (\[Delta]^ij+2 \[Delta]Z^ij)\!\(
\*SubscriptBox[\(\[Del]\), \(\[Mu]\)]\(R[i]\)\) \[Del]^\[Mu]R[j]
*)
(*
Notice that the kinetic term is different for the Higgs components;
this reflects that \[CurlyPhi]!=0 breaks the symmetry;
so the Higgs and Goldstone terms become different
*)

PrintScalarKineticHET[]//Normal


(* ::Chapter:: *)
(*Integrating out Triplet bosons*)


(*
Integrate out the triplet scalars
*)


PrepareHET[{5,6,7},{}]


CalculatePotentialHET[]


PrintActionHET["LO"]


PrintActionHET["NLO"]


PrintActionHET["NNLO"]


PrintScalarKineticHET[]//Normal
