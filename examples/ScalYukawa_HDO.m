(* ::Package:: *)

Quit[];


SetDirectory[NotebookDirectory[]];
DRalgo`DRalgo`$LoadGroupMath=True;
DRalgo`DRalgo`$GroupMathMultipleModels=True; (*Put this if you want to create multiple model-files with the same kernel*)
<<../Kernel/DRalgo.m


(* ::Chapter:: *)
(*A real-scalar + Dirac-fermion model with higher-order operator matching*)


(*
	See 2108.04377 [hep-th] for further details and for independent calculations
*)


(* ::Section::Closed:: *)
(*Model*)


Group={"U1"};
CouplingName={g};
RepAdjoint={0};
RepScalar={{{0},"R"}};
RepFermion={{{0},"L"},{{0},"R"}};
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


(*\[Kappa] \[Phi]^3 *)
InputInv={{1,1,1},{True,True,True}};
CubicTerm=CreateInvariant[Group,RepScalar,InputInv][[1]];
VCubic=\[Kappa] CubicTerm;
\[Lambda]3=GradCubic[VCubic];


(* ::Text:: *)
(*Scalar - Quartic terms*)


(* \[Lambda]\[Phi]^4 *)
InputInv={{1,1,1,1},{True,True,True,True}};
QuarticTerm=CreateInvariant[Group,RepScalar,InputInv][[1]];
VQuartic=\[Lambda] MassTerm^2;
\[Lambda]4=GradQuartic[VQuartic];


(* ::Text:: *)
(*Fermion-mass terms*)


(* m(Subscript[\[Psi], R]^+Subscript[\[Psi], L]+Subscript[\[Psi]^+, L]Subscript[\[Psi], R])*)
InputInv={{2,1},{False,True}}; (*Subscript[\[Psi], R]^+Subscript[\[Psi], L]*)
MassTerm1=CreateInvariantFermion[Group,RepFermion,InputInv][[1]];
InputInv={{1,2},{False,True}};  (*Subscript[\[Psi]^+, L]Subscript[\[Psi], R]*)
MassTerm2=CreateInvariantFermion[Group,RepFermion,InputInv][[1]];
\[Mu]IJ=m\[Psi]*GradMassFermion[MassTerm1];
\[Mu]IJC=m\[Psi]*GradMassFermion[MassTerm2];


(* ::Text:: *)
(*Yukawa terms*)


InputInv={{1,2,1},{True,False,True}};
YukawaDoublet1=CreateInvariantYukawa[Group,RepScalar,RepFermion,InputInv][[1]]//Simplify;
InputInv={{1,1,2},{True,False,True}}; 
YukawaDoublet2=CreateInvariantYukawa[Group,RepScalar,RepFermion,InputInv][[1]]//Simplify;
Ysff= y*GradYukawa[YukawaDoublet1];
YsffC=y*GradYukawa[YukawaDoublet2];


Ysff//MatrixForm
YsffC//MatrixForm
\[Mu]IJ//MatrixForm
\[Mu]IJC//MatrixForm


(* ::Section:: *)
(*Dimensional Reduction*)


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PerformDRhard[]


(* ::Section:: *)
(*DIMENSION 5 MATCHING*)


(*Here we construct the group tensors of the higher dimensional operators*)

OT[5,2]=\[Alpha][5,1]*TensorProduct[{1},{1},{1},{1},{1}];
OT[5,5]=\[Beta][5,1]*TensorProduct[{1},{1},{1}];

(*
	TensorList is the array of the various group tensors refered to
	higher dimensional operators
*)
TensorList5={OT[5,2],OT[5,5]};
(*
	NList is the array that indicate to which operator
	the group tensors in TensorList refer to
*)
NList5={2,5};
(*
	WC is the array of the various Wilson Coefficients \[Alpha][...]
*)
WC5 = {\[Alpha][5,1],\[Beta][5,1]};

sol5=Dimension5Matching[TensorList5,NList5,WC5,d];//AbsoluteTiming
Collect[sol5,{_Zb,_Zf},Factor]//TableForm


(* ::Section:: *)
(*DIMENSION 6 MATCHING*)


OT[6,4]=\[Beta][6,2]*TensorProduct[{1},{1},{1},{1}];
OT[6,9]=\[Alpha][6,1]*TensorProduct[{1},{1},{1},{1},{1},{1}];
OT[6,15]=\[Beta][6,1]*TensorProduct[{1},{1}];

(*
	TensorList is the array of the various group tensors refered to
	higher dimensional operators
*)
TensorList6={OT[6,4],OT[6,9],OT[6,15]};
(*
	NList is the array that indicate to which operator
	the group tensors in TensorList refer to
*)
NList6={4,9,15};
(*
	WC is the array of the various Wilson Coefficients \[Alpha][...]
*)
WC6 = {\[Alpha][6,1],\[Beta][6,1],\[Beta][6,2]};
(*
	TensorList is the array of the various group tensors refered to
	higher dimensional operators
*)
(*
	Dimension6Matching and Dimension5Matching find
	the values of the Wilson coefficients listed in WC,
	d is the number of spatial dimensions,
	Zb and Zf are 1 loop master integral
*)
sol6=Dimension6Matching[TensorList6,NList6,WC6,d][[1]];
Collect[sol6,{_Zb,_Zf},Factor]//TableForm


(*
	The matching can be done individually for each group tensor
*)
sol=Dimension5Matching[{OT[5,2]},{2},{\[Alpha][5,1]},d][[1]]; 
sol//Factor//TableForm


(*
	The matching can be done individually for each group tensor
*)
sol=Dimension6Matching[{OT[6,4],OT[6,15]},{4,15},{\[Beta][6,1],\[Beta][6,2]},d][[1]];
sol//Factor//TableForm


(*
	The functions ODIM6 and ODIM5 return the group tensors of the various operators
*)
ODIM5[5,d]//MatrixForm


(*
	The functions ODIM6 and ODIM5 return the group tensors of the various operators
*)
ODIM6[4,d]//MatrixForm


(*
	The solutions returned by the functions Dimension5Matching and Dimension6Matching are
	consistent with the ones reported in 2406.02667
*)
