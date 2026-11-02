(* ::Package:: *)

(* ::Subsection:: *)
(*Basic Example of the Abelian Higgs Model*)


(* ::Subsubsection:: *)
(*Model*)


(* ::Text:: *)
(*We consider a simple model with a U(1) gauge group. The adjoint representation is trivial for abelian groups and has Dynkin index {0}. We also include a single complex scalar in the fundamental representation (R = 1):*)


Group={"U1"};
RepAdjoint={0}; (*Dynkin index of the adjoint*)
CouplingName={g1}; (*Gauge coupling*)


(* ::Text:: *)
(*One complex scalar (interpreted as a Higgs-like field) is present:*)


Higgs={{Y\[Phi]},"C"}; (*Complex scalar with Hypercharge Y\[Phi]*)
RepScalar={Higgs}; (*Representation of the scalar field*)


(* ::Text:: *)
(*Fermions are treated as Weyl spinors in DRalgo. Therefore, to define a Dirac fermion, one must include both a left-handed and a right-handed Weyl spinor. In this minimal model, no fermions are included:*)


RepFermion={}; (*No fermions in the model*)


(* ::Text:: *)
(*We allocate all necessary coupling tensors for the model:*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];


(* ::Text:: *)
(*To build the Lagrangian, we construct the mass term m\.b2 \[Phi]\[Dagger]\[Phi]:*)


InputInv={{1,1},{True,False}}; (*Indicates a \[Phi]\[Dagger]\[Phi]-type invariant*)
MassTerm1=CreateInvariant[Group,RepScalar,InputInv]//Simplify//FullSimplify;
VMass=msq*MassTerm1[[1]]; (*Component-level expression for \[Phi]\[Dagger]\[Phi]*)
\[Mu]ij=GradMass[VMass]//Simplify//SparseArray;


(* ::Text:: *)
(*We add the quartic interaction of the complex scalar:*)


QuarticTerm1=MassTerm1[[1]]^2; (*Because MassTerm1=\[Phi]^+\[Phi], we can write (\[Phi]^+\[Phi])^2=MassTerm1^2*)
VQuartic=\[Lambda]*QuarticTerm1;
\[Lambda]4=GradQuartic[VQuartic];


(* ::Text:: *)
(*To invoke the model in DRalgo, it needs to be imported:*)


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->True];


(* ::Subsubsection:: *)
(*Dimensional reduction from hard to soft scale*)


(* ::Text:: *)
(*Integrating out the hard-scale modes in verbose mode:*)


PerformDRhard[]


(* ::Text:: *)
(*After the reduction completes, the effective couplings are given by:*)


PrintCouplings[]


(* ::Text:: *)
(*Internally, a list of constants is used and they can be printed via:*)


PrintConstants[]


(* ::Text:: *)
(*Both the leading order (LO) and next-to-leading order (NLO) scalar masses are given by:*)


PrintScalarMass["LO"]


PrintScalarMass["NLO"]


(* ::Text:: *)
(*Similarly, the Debye masses are given by:*)


PrintDebyeMass["LO"]


PrintDebyeMass["NLO"]


(* ::Text:: *)
(*Often one needs the explicit coupling tensors of DRalgo whose indices can be accessed via:*)


PrintTensorDRalgo[];


(* ::Text:: *)
(*For a specific tensor, e.g. for the scalars quartics with index "1", we can print it directly:*)


PrintTensorDRalgo[1]//MatrixForm


(* ::Text:: *)
(*The hard-scale pressure contribution in the high-temperature expansion is obtained by:*)


PrintPressure["LO"]


PrintPressure["NLO"]


PrintPressure["NNLO"]


(* ::Text:: *)
(*Finally, the \[Beta]-functions of the model are:*)


BetaFunctions4D[]


(* ::Subsubsection:: *)
(*Dimensional reduction from hard to ''ultrasoft'' scale*)


(* ::Text:: *)
(*Integrate out the soft-scale modes in verbose mode:*)


PerformDRsoft[{}];


(* ::Text:: *)
(*After the reduction completes, the effective couplings are given by:*)


PrintCouplingsUS[]


(* ::Text:: *)
(*Both the leading order (LO) and next-to-leading order (NLO) scalar masses are given by:*)


PrintScalarMassUS["LO"]


PrintScalarMassUS["NLO"]


(* ::Text:: *)
(*The soft-scale pressure contribution in the high-temperature expansion is obtained by:*)


PrintPressureUS["LO"]


PrintPressureUS["NLO"]


(* ::Text:: *)
(*The NNLO soft-scale pressure is not yet included in DRalgo:*)


PrintPressureUS["NNLO"]


(* ::Text:: *)
(*Finally, the \[Beta]-functions of the model are:*)


BetaFunctions3DUS[]


(* ::Subsubsection:: *)
(*Effective potential directly*)


(* ::Text:: *)
(*First, the theory scale is chosen:*)


UseUltraSoftTheory[];


(* ::Text:: *)
(*Then, a constant background is introduced to the scalars and the potential is calculated:*)


\[CurlyPhi]vev={0,\[CurlyPhi]}//SparseArray;
DefineVEVS[\[CurlyPhi]vev];
CalculatePotential[]


(* ::Text:: *)
(*Finally, the different orders of the potential can be computed:*)


PrintEffectivePotential["LO"]


PrintEffectivePotential["NLO"]


PrintEffectivePotential["NNLO"]


(* ::Subsubsection:: *)
(*Effective potential directly using HEFT*)


(* ::Text:: *)
(*In this example, we assume that the temporal and vector degrees of freedom are heavy compared to the scalar sector. Thus, we integrate them out starting from the soft theory, ultimately constructing a High-Energy Effective Theory (HEFT).*)


(* ::Text:: *)
(*Start from the soft theory:*)


UseSoftTheory[];


(* ::Text:: *)
(*We introduce a background field, \[CurlyPhi], for the scalars and calculate the effective potential:*)


\[CurlyPhi]vev={0,\[CurlyPhi],0}//SparseArray;
DefineVEVS[\[CurlyPhi]vev];
CalculatePotential[]


(* ::Text:: *)
(*Before computing the final potential, we must identify which temporal scalars and gauge bosons are heavy and are to be integrated out. Use `PrepareHET[{scalarIndices}, {vectorIndices}]` to specify them:*)


PrepareHET[{3},{1}]


(* ::Text:: *)
(*This works similarly to `PerformDRsoft[{}]`, and the indices of the scalar and gauge fields are found with:*)


PrintScalarRepPositions[]


PrintGaugeRepPositions[]


(* ::Text:: *)
(*Afterwards, the effective potential is calculated in the HEFT setup:*)


CalculatePotentialHET[]


(* ::Text:: *)
(*Finally, the different orders of the potential can be computed:*)


PrintActionHET["LO"]


PrintActionHET["NLO"]


(* ::Text:: *)
(*The U(1) gauge bosons modify the scalar kinetic terms. The renormalization (wavefunction normalization) is encoded in the Z-factor:*)


PrintScalarKineticHET[]//Normal


(* ::Subsection:: *)
(*Higher-Dimensional Operators in the Abelian Higgs Model*)


(* ::Text:: *)
(*DRalgo matches the dimension-5 and dimension-6 operators of the three-dimensional theory at one loop. The matching starts from the four-dimensional model above, imported once more without verbose output:*)


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];


(* ::Text:: *)
(*Each operator enters through a group tensor that acts on the real components of the fields. The complex scalar decomposes as \[Phi] = V_i \[Phi]_i with V = (1, i)/\[Sqrt]2, and the single U(1) gauge field is projected by `id = {1}`:*)


V=(1/Sqrt[2]){1,I};
Vs=(1/Sqrt[2]){1,-I};
id={1};


(* ::Subsubsection:: *)
(*Matching onto the physical basis*)


(* ::Text:: *)
(*The physical (on-shell) basis contains no redundant operators, so its Wilson coefficients are gauge independent. Each coefficient `\[Beta][\[Ellipsis]]` is labeled by the field content of its operator, with the magnetic field strength F, the temporal photon B0, and the covariant derivative D. The group tensors of the operators read:*)


OTP[2]=4 \[Beta][\[Phi]^2 F^2] SymmetrizeTensor[TensorProduct[V,Vs,id,id],{{1,2}}];
OTP[3]=4 \[Beta][B0^2 F^2] TensorProduct[id,id,id,id];
OTP[4]=SymmetrizeTensor[\[Beta][D^2 \[Phi]^4] SymmetrizeTensor[TensorProduct[V,Vs,V,Vs],{{1,2},{3,4}}],"P"];
OTP[5]=SymmetrizeTensor[4 \[Beta][D^2 \[Phi]^2 B0^2] SymmetrizeTensor[TensorProduct[V,Vs,id,id],{{1,2}}],"SS"];
OTP[8]=6! \[Beta][\[Phi]^6] Symmetrize[TensorProduct[Vs,Vs,Vs,V,V,V]];
OTP[9]=2*4! \[Beta][B0^2 \[Phi]^4] TensorProduct[Symmetrize[TensorProduct[Vs,Vs,V,V]],id,id];
OTP[10]=2*4! \[Beta][B0^4 \[Phi]^2] TensorProduct[Symmetrize[TensorProduct[Vs,V]],id,id,id,id];
OTP[11]=6! \[Beta][B0^6] TensorProduct[id,id,id,id,id,id];


(* ::Text:: *)
(*The list of tensors, the positions of the operators in the basis, and the Wilson coefficients define the matching problem:*)


TensorListP={OTP[2],OTP[3],OTP[4],OTP[5],OTP[8],OTP[9],OTP[10],OTP[11]};
NListP={2,3,4,5,8,9,10,11};
WCsP=DeleteDuplicates[Cases[Normal[TensorListP],\[Beta][__],\[Infinity]]]


(* ::Text:: *)
(*`Dimension6Matching` solves for the coefficients in `d` spatial dimensions. They are proportional to the bosonic one-loop sum-integral `Zb[3,0]`:*)


sol6P=Dimension6Matching["PHYS",TensorListP,NListP,WCsP,d][[1]];
TableForm[Collect[sol6P,{_Zb},Factor]]


(* ::Text:: *)
(*In three spatial dimensions, the sum-integral evaluates to:*)


HardThermal1LoopInt["B",3,0,3]


(* ::Text:: *)
(*which fixes the Wilson coefficients of the three-dimensional theory:*)


TableForm[Factor[sol6P/.d->3/.Zb[3,0]->HardThermal1LoopInt["B",3,0,3]]]


(* ::Subsubsection:: *)
(*Matching onto the off-shell basis*)


(* ::Text:: *)
(*Individual operators can also be matched onto the off-shell basis of [arXiv:2503.18904](https://arxiv.org/abs/2503.18904), which retains the redundant operators. There, the coefficient `\[Alpha][\[Phi]^6]` of the sextic operator depends on the gauge parameter `xi`, in contrast to the physical coefficient `\[Beta][\[Phi]^6]` above:*)


Tens9=6! \[Alpha][\[Phi]^6] Symmetrize[TensorProduct[Vs,Vs,Vs,V,V,V]];
Collect[Dimension6Matching[{Tens9},{9},{\[Alpha][\[Phi]^6]},d][[1]],{_Zb},Factor]
