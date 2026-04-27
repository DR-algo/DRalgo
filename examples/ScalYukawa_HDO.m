(* ::Package:: *)

Quit[];


SetDirectory[NotebookDirectory[]];
DRalgo`DRalgo`$LoadGroupMath=True;
<<../Kernel/DRalgo.m


(* ::Title:: *)
(*A real-scalar + Dirac-fermion model; model created by Oliver Gould  and Joonas Hirvonen (07-07-2022)*)


(*See 2108.04377 [hep-th] for further details and for independent calculations*)


(*<<../Higher_Dimensional_Operators_v1.wl (*HIGHER DIMENSIONAL OPERATORS ROUTINES*)*)


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


id={1};


(*Here we construct the group tensors of the higher dimensional operators*)

Tens2=Factorial[5]\[Alpha][\[Phi]^5]TensorProduct[id,id,id,id,id];
Tens5=Factorial[2]\[Alpha][\[Phi]^3 D^2]TensorProduct[id,id,id];


TensorList5={Tens2,Tens5}; (*TensorList is the array of the various group tensors refered to higher dimensional operators*)
NList5={2,5}; (*NList is the array that indicate to which operator the group tensors in TensorList refer to*)
WC5 = DeleteDuplicates[Cases[TensorList5//Normal, \[Alpha][__], \[Infinity]]]; (*WC is the array of the various Wilson Coefficients \[Alpha][...]*)


sol=DIMENSION5MATCHING[TensorList5,NList5,WC5,d](*[[1]]*) (*DIMENSION6MATCHING and DIMENSION5MATCHING find the values of the Wilson coefficients listed in WC, d is the number of spatial dimensions, Zb and Zf are 1 loop master integral *)
sol//Factor//TableForm


sol=DIMENSION5MATCHING[{Tens2},{2},{\[Alpha][\[Phi]^5]},d][[1]];  (*The matching can be done singularly for each group tensor*)
sol//Factor//TableForm


ODIM5[2,d]//MatrixForm (*The functions ODIM6 and ODIM5 return the group tensors of the various operators*)


(* ::Section:: *)
(*DIMENSION 6 MATCHING*)


id={1};


(*Here we construct the group tensors of the higher dimensional operators*)

Tens4=-3*Factorial[2]Factorial[2]\[Alpha][\[Phi]^4D^2]TensorProduct[id,id,id,id];
Tens9=Factorial[6]\[Alpha][\[Phi]^6]TensorProduct[id,id,id,id,id,id];
Tens15=Factorial[2]\[Alpha][\[Phi]^2D^4]TensorProduct[id,id];


TensorList6={Tens4,Tens9,Tens15}; (*TensorList is the array of the various group tensors refered to higher dimensional operators*)
NList6={4,9,15}; (*NList is the array that indicate to which operator the group tensors in TensorList refer to*)
WC6 = DeleteDuplicates[Cases[TensorList6//Normal, \[Alpha][__], \[Infinity]]]; (*WC is the array of the various Wilson Coefficients \[Alpha][...]*)


sol=DIMENSION6MATCHING[TensorList6,NList6,WC6,d][[1]]; (*DIMENSION6MATCHING and DIMENSION5MATCHING find the values of the Wilson coefficients listed in WC, d is the number of spatial dimensions, Zb and Zf are 1 loop master integral *)
sol//Factor//TableForm


sol=DIMENSION6MATCHING[{Tens4},{4},{\[Alpha][\[Phi]^4 D^2]},d][[1]]; (*The matching can be done singularly for each group tensor*)
sol//Factor//TableForm


ODIM6[4,d]//MatrixForm (*The functions ODIM6 and ODIM5 return the group tensors of the various operators*)


(* ::Section:: *)
(*CROSSCHECK WITH LITERATURE*)


(* ::Text:: *)
(*The effects of Higher Dimensional Operators for Scalar-Yukawa model has already been investigated in  arxiv 2406.02667*)
