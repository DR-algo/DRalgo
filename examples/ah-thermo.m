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


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJC,\[Mu]IJ,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];


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


PrintScalarMassUS["NLO"]


(* ::Section:: *)
(*Thermodynamics*)


(* 

A possible application of DRalgo output, by Andreas Ekstedt, Philipp Schicho and Tuomas V. I. Tenkanen, May 2022.
DRalgo tutorial Appendix A.1 covers this example.  

In this pedagogic example, we use the output of DRalgo to determine selected
thermodynamic properties  of the Abelian Higgs Model at high temperature,
by using perturbation theory within the 3d EFT. 

Required output are 

1) 4d theory beta functions, 
2) 3d EFT mathing relations,
3) 3d EFT effective potential.  

We implement these in functions

1) solveBetas[]
2) DRstep[]
3) Veff3d[]

Selected thermodynamic properties, the critical temperature Tc, critical field value at Tc (\[Phi]c/Tc) 
and released latent heat at Tc (L/Tc^4) are determined by a function 

4) findThermo[]

We emphasize that following perturbative determination is done at Landau gauge, and is gauge dependent.
This serves as a familiar and simple example;
a recipe for more sophisticated, gauge invariant determination can be found e.g. in [2203.04284 [hep-ph]]
{also cf. references therein), but is not implemented here, in order to maintain minimality.
We encourage DRalgo users to develop and optimise their own implementations for algorithms to determine 
thermodynamic properties. For future versions of DRalgo, we envision a possibility to include efficient
algorithms for determination of thermodynamics, but in the current \[Beta]-version,
these features are not implemented.

Below we read the output of DRalgo, and use it to construct the aformentioned functions.
For detailed documentation, see commentary within the implementation below,
as well as Appendix A.1 of the DRalgo manual.  
*)


(*
	Beta functions in 4d theory can be read from
*)
BetaFunctions4D[];
(*
	Results can be ported as
*)
sub1={g1->Sqrt[gsq[\[CapitalLambda]]],\[Lambda]->\[Lambda]h[\[CapitalLambda]],msq->\[Mu]hsq[\[CapitalLambda]]}; 
(* 
	Above is simply a substition rule to apply preferred notation. 
	Parameters are here functions of 4d renormalisation scale, denoted here by \[CapitalLambda]. 
*)
\[Beta]gsq=g1^2/.BetaFunctions4D[]/.sub1;
\[Beta]\[Lambda]=\[Lambda]/.BetaFunctions4D[]/.sub1//Expand;
\[Beta]msq=msq/.BetaFunctions4D[]/.sub1//Expand;

(*
	We implement these in a following function, that takes for
	an argument a list of three 4d MSbar parameters,
	defined at a given initial scale.
*)
solveBetas[{gsqInit_,\[Mu]hsqInit_,\[Lambda]hInit_}]:=Module[{initialScale,\[CapitalLambda],gsq,\[Lambda]h,\[Mu]hsq,eqgsq,eqg1sq,eq\[Mu]hsq,eq\[Lambda]h,Y\[Phi]},

Y\[Phi]=1;

eqgsq=(Y\[Phi]^2 gsq[\[CapitalLambda]]^2)/(24 \[Pi]^2); 
eq\[Mu]hsq=-((3 Y\[Phi]^2 gsq[\[CapitalLambda]] \[Mu]hsq[\[CapitalLambda]])/(8 \[Pi]^2))+(\[Lambda]h[\[CapitalLambda]] \[Mu]hsq[\[CapitalLambda]])/(2 \[Pi]^2);
eq\[Lambda]h=(3 Y\[Phi]^4 gsq[\[CapitalLambda]]^2)/(8 \[Pi]^2)-(3 Y\[Phi]^2 gsq[\[CapitalLambda]] \[Lambda]h[\[CapitalLambda]])/(4 \[Pi]^2)+(5 \[Lambda]h[\[CapitalLambda]]^2)/(4 \[Pi]^2);

(*
	Solve the system of beta functions using the initial condition.
*)
(* an arbitrary choice, hard-coded *)
initialScale=100; 

{gsq,\[Mu]hsq,\[Lambda]h}=NDSolveValue[{
(* beta functions *)
\[CapitalLambda]*D[gsq[\[CapitalLambda]],\[CapitalLambda]]==eqgsq,
\[CapitalLambda]*D[\[Mu]hsq[\[CapitalLambda]],\[CapitalLambda]]==eq\[Mu]hsq,
\[CapitalLambda]*D[\[Lambda]h[\[CapitalLambda]],\[CapitalLambda]]==eq\[Lambda]h,

(* run from initialScale *)
gsq[initialScale]==gsqInit,
\[Mu]hsq[initialScale]==\[Mu]hsqInit,
\[Lambda]h[initialScale]==\[Lambda]hInit},

{gsq,\[Mu]hsq,\[Lambda]h},
(*
	the lower and upper limits for renormalisation scale are tuned for
	the problem at hand, but otherwise arbitrary.
*)
{\[CapitalLambda],5,5000} 
];

(* Function returns MSbar parameters as Interpolating function of \[CapitalLambda]. *)
Return[{gsq,\[Mu]hsq,\[Lambda]h}]; 
];

(* 
Next, we can read 3d matching relations from Dralgo, using 

PrintCouplings[],
PrintScalarMass[],
PrintDebyeMass[],
PrintTemporalScalarCouplings[],

with following hacks for preferred notation: 
*)

(* Soft scale theory parameters: *)

sub2={g1->Sqrt[gsq],\[Lambda]->\[Lambda]h,msq->\[Mu]hsq}; (* preferred notation *)
gsq3d=g13d^2/.PrintCouplings[]/.sub2;
\[Lambda]h3d=\[Lambda]3d/.PrintCouplings[]/.sub2;
msq3dLO=msq3d/.PrintScalarMass["LO"]/.sub2;
msq3dNLO=msq3d/.PrintScalarMass["NLO"]/.sub2;
\[Mu]hsq3d=msq3dLO+msq3dNLO;
\[Mu]sqU1LO=\[Mu]sqU1/.PrintDebyeMass["LO"]/.sub2;
\[Mu]sqU1NLO=\[Mu]sqU1/.PrintDebyeMass["NLO"]/.sub2;
mDsq3d=\[Mu]sqU1LO+\[Mu]sqU1NLO;
\[Kappa]33d=\[Lambda]VLL/.PrintTemporalScalarCouplings[]/.sub2;
h33d=\[Lambda]VL[1]/.PrintTemporalScalarCouplings[]/.sub2;

(*
In a similar manner, we can use 

PrintCouplingsUS[],
PrintScalarMassUS[]

*)

(* to get ultra soft theory parameters: *)

bar\[Lambda]h3d=\[Lambda]3dUS/.PrintCouplingsUS[]/.\[Lambda]3d->\[Lambda]h3d/.\[Lambda]VL[1]->h33d/.\[Mu]sqU1->mDsq3d/.sub2;
bargsq3d=g13dUS^2/.PrintCouplingsUS[]/.g13d^2->gsq3d;

PrintScalarMassUS["LO"];

msq3dUSLO=msq3dUS/.PrintScalarMassUS["LO"]/.\[Lambda]3d->\[Lambda]h3d/.\[Lambda]VL[1]->h33d/.\[Mu]sqU1->mDsq3d/.sub2;
msq3dUSNLO=msq3dUS/.PrintScalarMassUS["NLO"]/.\[Lambda]3d->\[Lambda]h3d/.\[Lambda]VL[1]->h33d/.\[Mu]sqU1->mDsq3d/.sub2;

bar\[Mu]hsq3d=msq3dSSLO+msq3dSSNLO/.msq3d->\[Mu]hsq3d/.g13d->Sqrt[gsq3d]/.\[Lambda]3d->\[Lambda]h3d/.\[Lambda]VL[1]->h33d/.\[Mu]sqU1->mDsq3d/.sub2;

(* Expressions for ultra soft theory parameters bargsq3d, bar\[Lambda]h3d and bar\[Mu]hsq3d are implemented below in DRstep[].  *)

(* Two-loop effective potential in 3d EFT can be directly read from: *)

PrintEffectivePotential["LO"];
PrintEffectivePotential["NLO"];
PrintEffectivePotential["NNLO"];
(*
	we implement these expressions below in Veff3d[].
	Note that here we simply use predefined notation of DRalgo.
*)


(* Soft scale theory parameters: *)

sub2={g1->Sqrt[gsq],\[Lambda]->\[Lambda]h,msq->\[Mu]hsq}; (* preferred notation *)
gsq3d=g13d^2/.PrintCouplings[]/.sub2;
\[Lambda]h3d=\[Lambda]3d/.PrintCouplings[]/.sub2;
msq3dLO=msq3d/.PrintScalarMass["LO"]/.sub2;
msq3dNLO=msq3d/.PrintScalarMass["NLO"]/.sub2;
\[Mu]hsq3d=msq3dLO+msq3dNLO;
\[Mu]sqU1LO=\[Mu]sqU1/.PrintDebyeMass["LO"]/.sub2;
\[Mu]sqU1NLO=\[Mu]sqU1/.PrintDebyeMass["NLO"]/.sub2;
mDsq3d=\[Mu]sqU1LO+\[Mu]sqU1NLO;
\[Kappa]33d=\[Lambda]VLL/.PrintTemporalScalarCouplings[]/.sub2;
h33d=\[Lambda]VL[1]/.PrintTemporalScalarCouplings[]/.sub2;

(*
In a similar manner, we can use 

PrintCouplingsUS[],
PrintScalarMassUS[]

*)

(* to get ultra soft theory parameters: *)

bar\[Lambda]h3d=\[Lambda]3dUS/.PrintCouplingsUS[]/.\[Lambda]3d->\[Lambda]h3d/.\[Lambda]VL[1]->h33d/.\[Mu]sqU1->mDsq3d/.sub2;
bargsq3d=g13dUS^2/.PrintCouplingsUS[]/.g13d^2->gsq3d;

PrintScalarMassUS["LO"];

msq3dSSLO=msq3dUS/.PrintScalarMassUS["LO"]/.\[Lambda]3d->\[Lambda]h3d/.\[Lambda]VL[1]->h33d/.\[Mu]sqU1->mDsq3d/.sub2;
msq3dSSNLO=msq3dUS/.PrintScalarMassUS["NLO"]/.\[Lambda]3d->\[Lambda]h3d/.\[Lambda]VL[1]->h33d/.\[Mu]sqU1->mDsq3d/.sub2;

(*We choose \[Mu]3=mD to minimize logarithms*)
bar\[Mu]hsq3dpre=msq3dSSLO+msq3dSSNLO/.msq3d->\[Mu]hsq3d/.g13d->Sqrt[gsq3d]/.\[Lambda]3d->\[Lambda]h3d/.\[Lambda]VL[1]->h33d/.\[Mu]sqU1->mDsq3d/.sub2/.\[Mu]3->Sqrt[mDsq3d];
bar\[Mu]hsq3dRG=(msq3dUS/.BetaFunctions3DUS[])Log[\[Mu]3US/Sqrt[mDsq3d]]/.g13dUS^4->(g13dUS^2)^a/.PrintCouplingsUS[]/.a->1/.\[Lambda]3d->\[Lambda]h3d/.\[Lambda]VL[1]->h33d/.\[Mu]sqU1->mDsq3d/.sub2/.\[Mu]3->Sqrt[mDsq3d];
bar\[Mu]hsq3d=bar\[Mu]hsq3dpre+bar\[Mu]hsq3dRG/.PrintTemporalScalarCouplings[]/.g13d->Sqrt[gsq3d]/.\[Lambda]3d->\[Lambda]h3d/.\[Lambda]VL[1]->h33d/.\[Mu]sqU1->mDsq3d/.sub2;
(* Expressions for ultra soft theory parameters bargsq3d, bar\[Lambda]h3d and bar\[Mu]hsq3d are implemented below in DRstep[].  *)

(* Two-loop effective potential in 3d EFT can be directly read from: *)

PrintEffectivePotential["LO"];
PrintEffectivePotential["NLO"];
PrintEffectivePotential["NNLO"];
(*
	we implement these expressions below in Veff3d[].
	Note that here we simply use predefined notation of DRalgo.
*)


(*
	NLO dimensional reduction for the Abelian Higgs Model
*)
(* 
	Here function arguments are temperature T and 4d renormalisation scale (denoted herein by \[Mu]), 
	given as numbers, and MSbar parameters as outputted by solveBetas[]. 
*)
DRstep[T_,\[Mu]_,gsq_,\[Mu]hsq_,\[Lambda]h_]:=Module[{Lb,\[Mu]3US,bargsq3,bar\[Mu]hsq3,bar\[Lambda]h3,Y\[Phi],Y\[Phi]3d},

Y\[Phi]=1;
Y\[Phi]3d=1;

Lb=2 Log[\[Mu]/T]-2(Log[4\[Pi]]-EulerGamma); 

bargsq3=
	+gsq*T
	-(gsq^2 Lb T Y\[Phi]^2)/(48 \[Pi]^2);
bar\[Lambda]h3=
	+(T(gsq^2 (2-3 Lb) Y\[Phi]^4+6 gsq Lb Y\[Phi]^2 \[Lambda]h+2 \[Lambda]h (8 \[Pi]^2-5 Lb \[Lambda]h)))/(16 \[Pi]^2)
	-(gsq^2 T^2 Y\[Phi]^4 (48 \[Pi]^2-gsq (-4+Lb) Y\[Phi]^2+24 \[Lambda]h)^2)/(18432 \[Pi]^5 Sqrt[1/3 gsq T^2 Y\[Phi]^2+(gsq Y\[Phi]^2 (T^2 (-2 gsq (-7+Lb) Y\[Phi]^2+24 \[Lambda]h)+72 \[Mu]hsq))/(288 \[Pi]^2)]);
	
\[Mu]3US=bargsq3;


bar\[Mu]hsq3=1/12 T^2 (3 gsq Y\[Phi]^2+4 \[Lambda]h)+\[Mu]hsq-(gsq T Y\[Phi]^2 (48 \[Pi]^2-gsq (-4+Lb) Y\[Phi]^2+24 \[Lambda]h) Sqrt[1/3 gsq T^2 Y\[Phi]^2+(gsq Y\[Phi]^2 (T^2 (-2 gsq (-7+Lb) Y\[Phi]^2+24 \[Lambda]h)+72 \[Mu]hsq))/(288 \[Pi]^2)])/(192 \[Pi]^3)+(gsq T Y\[Phi]^2 (48 \[Pi]^2-gsq (-4+Lb) Y\[Phi]^2+24 \[Lambda]h) ((gsq^2 T Y\[Phi]^4)/\[Pi]^2-(gsq T Y\[Phi]^2 (48 \[Pi]^2-gsq (-4+Lb) Y\[Phi]^2+24 \[Lambda]h))/(12 \[Pi]^2)+(gsq T Y\[Phi]^2 (48 \[Pi]^2-gsq (-4+Lb) Y\[Phi]^2+24 \[Lambda]h) Log[2])/(6 \[Pi]^2)))/(3072 \[Pi]^4)+1/(4 \[Pi]^2) (Y\[Phi]^4 (gsq T-(gsq^2 Lb T Y\[Phi]^2)/(48 \[Pi]^2))-2 Y\[Phi]^2 (gsq T-(gsq^2 Lb T Y\[Phi]^2)/(48 \[Pi]^2)) ((T (gsq^2 (2-3 Lb) Y\[Phi]^4+6 gsq Lb Y\[Phi]^2 \[Lambda]h+2 \[Lambda]h (8 \[Pi]^2-5 Lb \[Lambda]h)))/(16 \[Pi]^2)-(gsq^2 T^2 Y\[Phi]^4 (48 \[Pi]^2-gsq (-4+Lb) Y\[Phi]^2+24 \[Lambda]h)^2)/(18432 \[Pi]^5 Sqrt[1/3 gsq T^2 Y\[Phi]^2+(gsq Y\[Phi]^2 (T^2 (-2 gsq (-7+Lb) Y\[Phi]^2+24 \[Lambda]h)+72 \[Mu]hsq))/(288 \[Pi]^2)]))+2 ((T (gsq^2 (2-3 Lb) Y\[Phi]^4+6 gsq Lb Y\[Phi]^2 \[Lambda]h+2 \[Lambda]h (8 \[Pi]^2-5 Lb \[Lambda]h)))/(16 \[Pi]^2)-(gsq^2 T^2 Y\[Phi]^4 (48 \[Pi]^2-gsq (-4+Lb) Y\[Phi]^2+24 \[Lambda]h)^2)/(18432 \[Pi]^5 Sqrt[1/3 gsq T^2 Y\[Phi]^2+(gsq Y\[Phi]^2 (T^2 (-2 gsq (-7+Lb) Y\[Phi]^2+24 \[Lambda]h)+72 \[Mu]hsq))/(288 \[Pi]^2)]))^2) Log[\[Mu]3US/Sqrt[1/3 gsq T^2 Y\[Phi]^2+(gsq Y\[Phi]^2 (T^2 (-2 gsq (-7+Lb) Y\[Phi]^2+24 \[Lambda]h)+72 \[Mu]hsq))/(288 \[Pi]^2)]]+1/(576 \[Pi]^2) (12 gsq Y\[Phi]^2 (Lb (-6 T^2 \[Lambda]h+9 \[Mu]hsq)+2 T^2 \[Lambda]h (1+6 EulerGamma-72 Log[Glaisher]))+24 \[Lambda]h (Lb (T^2 \[Lambda]h-6 \[Mu]hsq)-6 T^2 \[Lambda]h (EulerGamma-12 Log[Glaisher]))+gsq^2 T^2 Y\[Phi]^4 (-8-108 EulerGamma+69 Lb+1296 Log[Glaisher])+18 (8 Y\[Phi]^4 (gsq T-(gsq^2 Lb T Y\[Phi]^2)/(48 \[Pi]^2))^2+(gsq^2 T^2 Y\[Phi]^4 (48 \[Pi]^2-gsq (-4+Lb) Y\[Phi]^2+24 \[Lambda]h)^2)/(576 \[Pi]^4)-(T Y\[Phi]^2 (gsq T-(gsq^2 Lb T Y\[Phi]^2)/(48 \[Pi]^2)) (gsq^2 (2-3 Lb) Y\[Phi]^4+6 gsq Lb Y\[Phi]^2 \[Lambda]h+2 \[Lambda]h (8 \[Pi]^2-5 Lb \[Lambda]h)))/\[Pi]^2+(T^2 (gsq^2 (2-3 Lb) Y\[Phi]^4+6 gsq Lb Y\[Phi]^2 \[Lambda]h+2 \[Lambda]h (8 \[Pi]^2-5 Lb \[Lambda]h))^2)/(16 \[Pi]^4)) Log[Sqrt[1/3 gsq T^2 Y\[Phi]^2+(gsq Y\[Phi]^2 (T^2 (-2 gsq (-7+Lb) Y\[Phi]^2+24 \[Lambda]h)+72 \[Mu]hsq))/(288 \[Pi]^2)]/\[Mu]]);(* return *)
{bargsq3,bar\[Mu]hsq3,bar\[Lambda]h3} 
(* 
Function returns numeric values for ultra soft scale theory parameters, at fixed input temperature T 
and renormalisation scale \[Mu].
 *)

];

(* 
Landau gauge 3d Veff at 2-loops, as function of a generic background field \[Phi] and
ultrasoft theory parameters (at fixed temperature). 
Arguments LO, NLO and NNLO are assumed to be flags that control if corresponding pieces are included, 
for example for Veff at 1-loop level one sets LO=1, NLO=1 and NNLO=0.
Background field \[Phi] is assumed to be given as undetermined variable, 
and resulting function is a single-variable function of \[Phi].
*)
Veff3d[\[Phi]_,bargsq3_,bar\[Mu]hsq3_,bar\[Lambda]h3_,LO_,NLO_,NNLO_]:=Module[{veff,veffLO,veffNLO,veffNNLO,msq,\[Lambda],g1,\[Mu]3US, Y\[Phi]},

Y\[Phi]=1;
(*
	for simplicity, we hardcode here 3d RG scale to be fixed as \[Mu]3US = bargsq3.
	Naturally, other options are possible.
*)
\[Mu]3US=bargsq3; 

g1=Sqrt[bargsq3];
\[Lambda]=bar\[Lambda]h3;
msq=bar\[Mu]hsq3;

veffLO=
	+(msq \[Phi]^2)/2
	+(\[Lambda] \[Phi]^4)/4;


veffNLO=
	-((g1^2 Y\[Phi]^2 \[Phi]^2)^(3/2)/(6 \[Pi]))
	-(msq+\[Lambda] \[Phi]^2)^(3/2)/(12 \[Pi])
	-(msq+3\[Lambda] \[Phi]^2)^(3/2)/(12 \[Pi]);  

veffNNLO=
	+(g1^2 Y\[Phi]^2 Sqrt[g1^2 Y\[Phi]^2 \[Phi]^2] Sqrt[msq+\[Lambda] \[Phi]^2])/(16 \[Pi]^2)
	+(3\[Lambda](msq+\[Lambda] \[Phi]^2))/(64 \[Pi]^2)
	+(g1^2 Y\[Phi]^2 Sqrt[g1^2 Y\[Phi]^2 \[Phi]^2] Sqrt[msq+3 \[Lambda] \[Phi]^2])/(16 \[Pi]^2)
	+(\[Lambda] Sqrt[msq+\[Lambda] \[Phi]^2] Sqrt[msq+3 \[Lambda] \[Phi]^2])/(32 \[Pi]^2)
	+(3 \[Lambda] (msq+3 \[Lambda] \[Phi]^2))/(64 \[Pi]^2)
	-(3 \[Lambda]^2 \[Phi]^2 (1/2+Log[\[Mu]3US/(3 Sqrt[msq+3 \[Lambda] \[Phi]^2])]))/(16 \[Pi]^2)
	+g1^4 Y\[Phi]^4 \[Phi]^2 (1/(32 \[Pi]^2)
	-(1/2+Log[\[Mu]3US/(2 Sqrt[g1^2 Y\[Phi]^2 \[Phi]^2]+Sqrt[msq+3 \[Lambda] \[Phi]^2])])/(16 \[Pi]^2)
	+1/(4 g1^4 Y\[Phi]^4 \[Phi]^4) ((g1^2 Y\[Phi]^2 \[Phi]^2 Sqrt[g1^2 Y\[Phi]^2 \[Phi]^2] Sqrt[msq+3 \[Lambda] \[Phi]^2])/(8 \[Pi]^2)+(g1^2 Y\[Phi]^2 \[Phi]^2 (msq-2 g1^2 Y\[Phi]^2 \[Phi]^2+3 \[Lambda] \[Phi]^2))/(16 \[Pi]^2)-((msq+3 \[Lambda] \[Phi]^2)^2 (1/2+Log[\[Mu]3US/Sqrt[msq+3 \[Lambda] \[Phi]^2]]))/(16 \[Pi]^2)+((-msq+g1^2 Y\[Phi]^2 \[Phi]^2-3 \[Lambda] \[Phi]^2)^2 (1/2+Log[\[Mu]3US/(Sqrt[g1^2 Y\[Phi]^2 \[Phi]^2]+Sqrt[msq+3 \[Lambda] \[Phi]^2])]))/(8 \[Pi]^2)-((-msq+2 g1^2 Y\[Phi]^2 \[Phi]^2-3 \[Lambda] \[Phi]^2)^2 (1/2+Log[\[Mu]3US/(2 Sqrt[g1^2 Y\[Phi]^2 \[Phi]^2]+Sqrt[msq+3 \[Lambda] \[Phi]^2])]))/(16 \[Pi]^2)))+1/(2 \[Phi]^2) ((Sqrt[g1^2 Y\[Phi]^2 \[Phi]^2] Sqrt[msq+\[Lambda] \[Phi]^2] (-g1^2 Y\[Phi]^2 \[Phi]^2+2 \[Lambda] \[Phi]^2))/(16 \[Pi]^2)+(Sqrt[g1^2 Y\[Phi]^2 \[Phi]^2] (-g1^2 Y\[Phi]^2 \[Phi]^2-2 \[Lambda] \[Phi]^2) Sqrt[msq+3 \[Lambda] \[Phi]^2])/(16 \[Pi]^2)+(g1^2 Y\[Phi]^2 \[Phi]^2 Sqrt[msq+\[Lambda] \[Phi]^2] Sqrt[msq+3 \[Lambda] \[Phi]^2])/(16 \[Pi]^2)+(\[Lambda]^2 \[Phi]^4 (1/2+Log[\[Mu]3US/(Sqrt[msq+\[Lambda] \[Phi]^2]+Sqrt[msq+3 \[Lambda] \[Phi]^2])]))/(4 \[Pi]^2)-((g1^4 Y\[Phi]^4 \[Phi]^4-2 g1^2 Y\[Phi]^2 \[Phi]^2 (msq+\[Lambda] \[Phi]^2)+(msq+\[Lambda] \[Phi]^2)^2-2 g1^2 Y\[Phi]^2 \[Phi]^2 (msq+3 \[Lambda] \[Phi]^2)-2 (msq+\[Lambda] \[Phi]^2) (msq+3 \[Lambda] \[Phi]^2)+(msq+3 \[Lambda] \[Phi]^2)^2) (1/2+Log[\[Mu]3US/(Sqrt[g1^2 Y\[Phi]^2 \[Phi]^2]+Sqrt[msq+\[Lambda] \[Phi]^2]+Sqrt[msq+3 \[Lambda] \[Phi]^2])]))/(16 \[Pi]^2))-(\[Lambda]^2 \[Phi]^2 (1/2+Log[\[Mu]3US/(2 Sqrt[msq+\[Lambda] \[Phi]^2]+Sqrt[msq+3 \[Lambda] \[Phi]^2])]))/(16 \[Pi]^2) ;

veff=
	+LO*veffLO
	+NLO*veffNLO
	+NNLO*veffNNLO;

(* return *)
{veff}
(*
	function returns a list of single element,
	which is a single variable function of \[Phi] at fixed T.
*)
];

(*
Finally, we present a simple algorithm to find Tc, \[Phi]c/Tc and L/Tc^4. 
Critical temperature is determined from the condition that potential has degenerate minima. 
Since we normalise the effective potential to zero at the origin, i.e. at the symmetric phase \[Phi]=0, 
we need to find temperature where the effective potential vanishes at the broken minimum \[Phi]c\[NotEqual]0.
We implement simple binary search to find this temperature.       
*)

(* Find Tc from degenerate minima, using binary search:  *)
(* 
Here function arguments are "physical", tree-level Higgs mass, and \[Lambda] and g^2 at chosen initial scale, 
and scaleFactor N that controls 4d RG scale as \[CapitalLambda] = N*\[Pi]*T. This implements natural choice for hard scale ~\[Pi]*T.  
Flag 'plotting' can be set to plotting=1, to inspect that binary
search is actually converging to degenerate minima. Otherwise, this flag can be set to zero.  

Note that since we are working on a toy model instead of an actual BSM theory, here mass M can be in arbitrary units, 
and other dimensionfull quantities are assumed to be in same units.
*)
findThermo[M_,\[Lambda]_,gsq_,scaleFactor_,plotting_]:=Module[
	{Tc,phicTc,Tmax,Tmin,Tm,Tp,TT,Tloopmax,v1Limit,veff,veff1,veff2,veff0,\[Phi],minSoln,Vc,plot,veffNorm,\[Mu]sq,scale,\[Phi]c,Tnext,
	dT,val1,val2,dVbro,dVsym,L,lat,gsqInit,\[Mu]sqInit,\[Lambda]Init,gsqInt,\[Mu]sqInt,\[Lambda]Int,ngsq,n\[Mu]sq,n\[Lambda],T,bargsq3,bar\[Mu]hsq3,bar\[Lambda]h3,LO,NLO,NNLO,
	\[Phi]3,\[Phi]0,method,constraints,acc,scalingFactor},

(* binary search for critical temperature: *)

(* range for temperature interval. These numbers depend on a chosen parameter points (in this case fixed M=100 and varying \[Lambda] from 0.003 to 0.03), 
and are not optimised beyond our simple illustration below. *)
Tmax=700; 
Tmin=200;

(* Auxilliary varibles Tm, Tp, and TT: *)
Tm=Tmin;
Tp=Tmax;
TT=(Tm+Tp)/2; (* we start our binary search from average temperature of given T range. *)
Tloopmax=20; (* determines how many iterations we do in binary search. This must be large enough for decent accuracy.  *)

(* limits for background field range: *)

(*v1Limit = 10;*)
(* this number has to be adjusted for binary search to work. *)
v1Limit=200;

(* tree-level relation of MSbar mass parameter and "physical mass" M *)
\[Mu]sq=-((M)^2/2); 

{gsqInit,\[Mu]sqInit,\[Lambda]Init}={gsq,\[Mu]sq,\[Lambda]};
{gsqInt,\[Mu]sqInt,\[Lambda]Int}=solveBetas[{gsqInit,\[Mu]sqInit,\[Lambda]Init}];
(* we start by solving the beta functions *)

(* Here we set 3d EFT Veff accuracy flags to maximum available accuracy including the 2-loop contributions. *)
LO=1;
NLO=1;
NNLO=1;

(* Next, we do for-loop with lenght Tloopmax, for the binary search: *)
Do[
scale=scaleFactor*\[Pi]*TT; (* 4d scale is fixed to a number * (\[Pi]*T) *)
{ngsq,n\[Mu]sq,n\[Lambda]}={gsqInt[scale],\[Mu]sqInt[scale],\[Lambda]Int[scale]}; (* this fixes numerical values for MSbar parameters at this thermal scale. *)

{bargsq3,bar\[Mu]hsq3,bar\[Lambda]h3}=DRstep[TT,scale,ngsq,n\[Mu]sq,n\[Lambda]]; (* solves 3d EFT parameters at ultra soft scale. *)

\[Phi]3=\[Phi]/Sqrt[TT];  (* We express 3d field \[Phi]3, that has mass dimension [\[Phi]3]=1/2 , in terms of more familiar 4d field \[Phi], with dimension [\[Phi]]=1. *)
veff=Veff3d[\[Phi]3,bargsq3,bar\[Mu]hsq3,bar\[Lambda]h3,LO,NLO,NNLO]; (* construct and store 3d Veff. *)
veff1=veff[[1]]; (* here we pick the expression from the list, and store it to variable veff1. *)
\[Phi]0=0.0001; 
(* 
we normalise Veff to zero at close to the origin. 
This is a carefree implementation where we have not worked out actual expression for the 2-loop part of the Veff at limit \[Phi]\[Rule]0. 
Resorting to such shortcuts or hacks, one has to make sure the final numerical results are not sensitive to any similar shennanigans.        
*)
veff0=veff1/.{\[Phi]->\[Phi]0};
veffNorm=veff1-veff0;

veff2=veffNorm;
acc=10; (* set accuracy for DifferentialEvolution method *)
scalingFactor=0.6; (* Mathematica default is 0.6 *)
minSoln=NMinimize[{Re[veff2],\[Phi]>v1Limit},{\[Phi],1,2000},
	Method->{"DifferentialEvolution","SearchPoints"->acc,"ScalingFactor"->scalingFactor}
	];
(* 
For the minimisation of the potential, we use Mathematica's build in function NMinimize, with "DifferentialEvolution" method. 
Here we limit the search to (4d) field values between 1 and 2000.  
Alternative, and potentially faster, methods could be available, for example "RandomSearch" but in this simple example we not aim to optimise 
the performance of this algorithm, but instead provide something simple. 
However, such optimisations become crucially important for BSM theory applications, that contain multiple free
parameters and background fields. 
*)

Vc=minSoln[[1]]; (* value of the Veff at the broken phase minimum. *)

If[plotting==1,
(* note relation between 3d EFT and 4d Veffs: Veff4d = T Veff3d. *)
plot=Plot[{TT*Re[veff2]},{\[Phi],0,1500},
	AxesLabel->{"\[Phi]","Re[\!\(\*SubscriptBox[\(V\),\(eff\)]\)[\[Phi]]]"}];
Print[plot];
(* if plotting flag is set to unity, Veff is printed at each value of T in binary search. *)
];

(* 
Then the actual binary search:  
We aim to find Veff[\[Phi]broken] \[Equal] 0, so if at given T (denoted by auxilliary variable TT here) Vc is positive, 
we choose next guess between TT and Tm, and if Vc is negative, we choose next guess between TT and Tp, 
and adjust auxilliary variables TT, Tp and Tm accordingly, as follows:      
*)
If[Vc>0,
	Tnext=(TT+Tm)/2//N;
	Tp=TT;
	TT=Tnext;
];
If[Vc<0,
	Tnext=(TT+Tp)/2//N;
	Tm=TT;
	TT=Tnext;
];

(* 
At the end, we asumme to have landed close enough to Tc (naturally, Tloopmax has to be large enough here for any decent accuracy) 
*)
If[Tloop==Tloopmax,
\[Phi]c=Abs[\[Phi]]/.minSoln[[2]]; (* save the critical field value *)

Tc=TT//N; (* save the critical temperature *)

phicTc=\[Phi]c/Tc; (* Ratio \[Phi]c/Tc correlates with phase transition strenght, albeit being a gauge dependent quantity. *)

(* Keep in mind that also our estimate for Tc here is gauge dependent, and thereof not well-defined physical quantity.  *)


];
,{Tloop,1,Tloopmax,1}];

(* Finally, we also determine the latent heat, defined as  *)

(* 
L = \[CapitalDelta] T dVeff/dT , where we denote \[CapitalDelta] \[Equal] [low T phase] - [high T phase]
as difference of a quantity between two phases.
Note that above Veff is in 4d units.
*)

(* 
For simplicity, we compute the T-derivative as a simple finite difference.
For this to be numerically accurate,  we must be close enough to Tc and use small enough
finite difference T interval dT. In this kind of a hack implementation, 
one needs to make sure that final numerical results are not sensitive to such choices of
the implementation.
Again, for actual, reliable, BSM application, more careful implementation is a must. 
*)

(* we think to be happy with this T-interval for finite difference. *)
dT = 0.1; 

T=Tc-dT; 
(* repeat minimisation at T= Tc - dT: *)
scale=scaleFactor(\[Pi]*T);
{ngsq,n\[Mu]sq,n\[Lambda]}={gsqInt[scale],\[Mu]sqInt[scale],\[Lambda]Int[scale]};
{bargsq3,bar\[Mu]hsq3,bar\[Lambda]h3}=DRstep[T,scale,ngsq,n\[Mu]sq,n\[Lambda]];
\[Phi]3=\[Phi]/Sqrt[TT];  
veff=Veff3d[\[Phi]3,bargsq3,bar\[Mu]hsq3,bar\[Lambda]h3,LO,NLO,NNLO];
(*veff4d = T veff[[1]];*)
veff1=veff[[1]];
\[Phi]0=0.0001;
veff0=veff1/.{\[Phi]-> \[Phi]0};
veffNorm=veff1-veff0;
veff2=veffNorm;
acc=10;
scalingFactor=0.6;
(* Mathematica default is 0.6 *)
minSoln=NMinimize[{Re[veff2],\[Phi]>v1Limit},{\[Phi],1,2000},
	Method->{"DifferentialEvolution","SearchPoints"->acc,"ScalingFactor"->scalingFactor}
	];
(* alternative method: "RandomSearch" *)
Vc=minSoln[[1]]; (* value of the Veff at the minimum. *)
val2=T*Vc; (* convert to 4d quantity *)

T=Tc-2*dT; 
(* we then repeat minimisation at T= Tc - 2dT:  *)
scale=scaleFactor(\[Pi]*T);
{ngsq,n\[Mu]sq,n\[Lambda]}={gsqInt[scale],\[Mu]sqInt[scale],\[Lambda]Int[scale]};
{bargsq3,bar\[Mu]hsq3,bar\[Lambda]h3}=DRstep[T,scale,ngsq,n\[Mu]sq,n\[Lambda]];
\[Phi]3=\[Phi]/Sqrt[TT];  
veff=Veff3d[\[Phi]3,bargsq3,bar\[Mu]hsq3,bar\[Lambda]h3,LO,NLO,NNLO];
(*veff4d=T veff[[1]];*)
veff1=veff[[1]];
\[Phi]0=0.0001;
veff0=veff1/.{\[Phi]->\[Phi]0};
veffNorm=veff1-veff0;

veff2=veffNorm;
acc=10;
scalingFactor=0.6; (* Mathematica default is 0.6 *)
minSoln=NMinimize[{Re[veff2],\[Phi]>v1Limit},{\[Phi],1,2000},
	Method->{"DifferentialEvolution","SearchPoints"->acc,"ScalingFactor"->scalingFactor}
	];
(* alternative method: "RandomSearch" *)
Vc=minSoln[[1]]; (* value of the Veff at the minimum. *)
val1=T*Vc; (* convert to 4d quantity *)

(* Now, we estimate the T-derivative by finite difference:  *)
dVbro=(val2-val1)/dT;
dVsym=0; (* Note that this is zero since Vsym is normalised to zero for every T. *)

L=Tc(dVbro-dVsym)//Abs; (* result for the latent heat near Tc *)

lat=L/Tc^4; (* dimensionless ratio *)

(* return *)
{Tc,phicTc,lat}
(* 
Function returns estimates for Tc, \[Phi]c/Tc and L/Tc^4 at Tc. 
DISCLAIMER: note that here we wrote our algorithm in relatively carefree attitude, 
since we work in a simple toy-model setup, whereas for actual BSM application similar 
algoritm has to be optimised for speed, and reliability and numerical accuracy ensured carefully!
*)
];


(* Numerical example, for a fixed benchmark point: *)

M=100; (* We fix Higgs field mass to some arbitrary number (understood to be in some units of mass). *)
gsqInit=0.42; (* We fix the gauge coupling, again, this choice is arbitrary. *)
\[Lambda]hInit=0.005; (* Finally, we fine tune scalar self coupling to find strong first order phase transition, in perturbation theory. *)

(*
	Note that 3d EFT dynamics of this model are fully described by two dimensionless ratios
	x = \[Lambda]3/gsq3 and y = msq3/gsq3^2,
	cf. eg. [2112.08912 [hep-ph]] and references therein.  
*)

(*
	We fix renormalisation scale to be 1.0*\[Pi]*T
*)
scaleFactor=1; 
plotting=0;

Print["{Tc,phicTc,lat}"]
findThermo[M,\[Lambda]hInit,gsqInit,scaleFactor,plotting]


(*
	Finally, loop over \[Lambda], to mimick a scan over BSM input parameters:
*)

(*
In order to take a theoretical uncertainty related to varying RG scale seriously, 
(c.f. [2009.10080 [hep-ph]] and [2104.04399 [hep-ph]]) we perform the scan with two 
choices of RG scale, 0.5*\[Pi]*T and 2.0*\[Pi]*T and interpet the difference to be the minimal
intrinsic uncertainty in our perturbative computation.   
*)

M=100;
gsqInit=0.42;

\[Lambda]max=0.03;
\[Lambda]min=0.003;
d\[Lambda]=0.001; (* step size in \[Lambda] *)

index1=1; (* running index for arrays/lists *)
len=Ceiling[((\[Lambda]max-\[Lambda]min)/d\[Lambda]+1)]; (* length for arrays/lists to store results *)

(* Allocate lists to store results: *)

(* To store results with RG scale 0.5*\[Pi]*T *)
TcList1=ConstantArray[0,len]; 
phicTcList1=ConstantArray[0,len]; 
latentList1=ConstantArray[0,len];

(* To store results with RG scale 2.0*\[Pi]*T *)
TcList2=ConstantArray[0,len]; 
phicTcList2=ConstantArray[0,len]; 
latentList2=ConstantArray[0,len];

Do[
(* auxilliary printing to showcase the progress *)
Print[\[Lambda]];

\[Lambda]hInit=\[Lambda];
plotting=0;

scaleFactor=0.5;
{Tc,phicTc,lat}=findThermo[M,\[Lambda]hInit,gsqInit,scaleFactor,plotting];
TcList1[[index1]]={\[Lambda],Tc};
phicTcList1[[index1]]={\[Lambda],phicTc};
latentList1[[index1]]={\[Lambda],lat};

scaleFactor=2.0;
{Tc,phicTc,lat}=findThermo[M,\[Lambda]hInit,gsqInit,scaleFactor,plotting];
TcList2[[index1]]={\[Lambda],Tc};
phicTcList2[[index1]]={\[Lambda],phicTc};
latentList2[[index1]]={\[Lambda],lat};

index1=index1+1; 
,{\[Lambda],\[Lambda]min,\[Lambda]max,d\[Lambda]}];





(* Plot the results: *)
ListLinePlot[{TcList1,TcList2},AxesLabel->{"\[Lambda]","Tc"},PlotStyle->Black]
ListLinePlot[{phicTcList1,phicTcList2},AxesLabel->{"\[Lambda]","\!\(\*FractionBox[\(\[Phi]c\),\(Tc\)]\)"},PlotStyle->Black]
ListLinePlot[{latentList1,latentList2},AxesLabel->{"\[Lambda]","\!\(\*FractionBox[\(L\),SuperscriptBox[\(Tc\),\(4\)]]\)"},PlotStyle->Black]


(*
	Export the data to produce publication quality plots with Python
	by .examples/ah-thermo-python-plots/plot.py
*)
Export["./ah-thermo-python-plots/TcList1.dat",TcList1]
Export["./ah-thermo-python-plots/phicTcList1.dat",phicTcList1]
Export["./ah-thermo-python-plots/latentList1.dat",latentList1]

Export["./ah-thermo-python-plots/TcList2.dat",TcList2]
Export["./ah-thermo-python-plots/phicTcList2.dat",phicTcList2]
Export["./ah-thermo-python-plots/latentList2.dat",latentList2]



