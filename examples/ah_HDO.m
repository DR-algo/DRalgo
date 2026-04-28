(* ::Package:: *)

Quit[];


SetDirectory[NotebookDirectory[]];
DRalgo`DRalgo`$LoadGroupMath=True;
<<../Kernel/DRalgo.m


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


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];


InputInv={{1,1},{True,False}}; (*This specifies that we want a \[Phi]^+\[Phi] term*)
MassTerm1=CreateInvariant[Group,RepScalar,InputInv]//Simplify//FullSimplify;


VMass=msq*MassTerm1[[1]];(*This is the \[Phi]^+\[Phi] term written in component form*)


\[Mu]ij=GradMass[VMass]//Simplify//SparseArray;


QuarticTerm1=MassTerm1[[1]]^2; (*Because MassTerm1=\[Phi]^+\[Phi], we can write (\[Phi]^+\[Phi])^2=MassTerm1^2*)


VQuartic=\[Lambda]*QuarticTerm1;


\[Lambda]4=GradQuartic[VQuartic];


Y\[Phi]=1;


(* ::Section:: *)
(*Dimensional Reduction*)


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PerformDRhard[]


(* ::Section:: *)
(*Dimension 6 Matching*)

(*The complex scalar field \phi is decomposed in terms of two real scalar field \\phi_1 and \phi_2 as follows \phi=1/Sqrt[2](\phi_1+i\phi_2) which can be written as \phi=V_i \phi_i,
therefore we can construct all the tensors using the vector V={1,I}/Sqrt[2], its complex conjugate Vs={1,-I}/Sqrt[2] and the projector on the gauge field, which, since it's the only gauge field of the model, it's just id={1}*)

V=(1/Sqrt[2]){1,I}; 
Vs=(1/Sqrt[2]){1,-I};
id={1};

(*We define \[Alpha][__] the Wilson coefficients associated to the dimension-6 operator basis defined in arxiv 2503.18904, c[__] are just combination of them*)

Subscript[c, 1]=\[Alpha][ D^2 \[Phi]^4,1]-4\[Alpha][ D^2 \[Phi]^4,2];
Subscript[c, 2]=-\[Alpha][ D^2 \[Phi]^4,2];
Subscript[c, 3]=-2 \[Alpha][ D^2 \[Phi]^2 B0^2,2]+\[Alpha][ D^2 \[Phi]^2 B0^2,3];
Subscript[c, 4]=-\[Alpha][ D^2 \[Phi]^2 B0^2,1];
Subscript[c, 5]=-2 \[Alpha][ D^2 \[Phi]^2 B0^2,2]-\[Alpha][ D^2 \[Phi]^2 B0^2,1];
Subscript[c, 6]=-3 \[Alpha][ D^2 B0^4];


(*Here we construct the group tensors of the higher dimensional operators*)

Tens2=4\[Alpha][ \[Phi]^2 F^2]SymmetrizeTensor[TensorProduct[V,Vs,id,id],{{1,2}}];
Tens3=4\[Alpha][ B0^2 F^2]TensorProduct[id,id,id,id];
Tens4=4(Subscript[c, 1]SymmetrizeTensor[TensorProduct[V,Vs,V,Vs],{{1,2},{3,4}}]+Subscript[c, 2](SymmetrizeTensor[TensorProduct[V,V,Vs,Vs],{{1,2},{3,4}}]+SymmetrizeTensor[TensorProduct[Vs,Vs,V,V],{{1,2},{3,4}}]));
Tens5=4Subscript[c, 3]SymmetrizeTensor[TensorProduct[V,Vs,id,id],{{1,2}}];
Tens6=2Subscript[c, 5]SymmetrizeTensor[TensorProduct[V,Vs,id,id],{{1,2}}];
Tens7=4Subscript[c, 4]SymmetrizeTensor[TensorProduct[V,Vs,id,id],{{1,2}}];
Tens8=4Subscript[c, 6]TensorProduct[id,id,id,id];
Tens9=\[Alpha][ \[Phi]^6]Factorial[6]Symmetrize[TensorProduct[Vs,Vs,Vs,V,V,V]];
Tens10=2*Factorial[4]\[Alpha][ B0^2 \[Phi]^4]TensorProduct[Symmetrize[TensorProduct[Vs,Vs,V,V]],id,id];
Tens11=2*Factorial[4]\[Alpha][ B0^4 \[Phi]^2]TensorProduct[Symmetrize[TensorProduct[Vs,V]],id,id,id,id];
Tens12=Factorial[6] \[Alpha][ B0^6]TensorProduct[id,id,id,id,id,id];
Tens13=2\[Alpha][ D^2 F^2]TensorProduct[id,id];
Tens14=-I*\[Alpha][ D^2 \[Phi]^2 F](TensorProduct[V,Vs,id]-TensorProduct[Vs,V,id])/2;
Tens15=2\[Alpha][ D^4 \[Phi]^2]Symmetrize[TensorProduct[V,Vs]];
Tens17=2\[Alpha][ D^4 B0^2]TensorProduct[id,id];


TensorList={Tens2,Tens3,Tens4,Tens5,Tens6,Tens7,Tens8,Tens9,Tens10,Tens11,Tens12,Tens13,Tens14,Tens15,Tens17}; (*TensorList is the array of the various group tensors refered to higher dimensional operators*)
NList={2,3,4,5,6,7,8,9,10,11,12,13,14,15,17}; (*NList is the array that indicate to which operator the group tensors in TensorList refer to*)
WC = DeleteDuplicates[Cases[TensorList//Normal, \[Alpha][__], \[Infinity]]]; (*WC is the array of the various Wilson Coefficients \[Alpha][...]*)


sol=DIMENSION6MATCHING[TensorList,NList,WC,d][[1]]; (*DIMENSION6MATCHING and DIMENSION5MATCHING find the values of the Wilson coefficients listed in WC, d is the number of spatial dimensions, Zb and Zf are 1 loop master integral *)
sol//Factor//TableForm

(*The solutions of the Wilson Coefficient is left in generic R-\xi gauge, the user can easily check that this gauge dependence is removed once redundancy of the operator basis is removed as well *)
(*Zb[3,0] is an Hard Thermal 1 loop Integral, its value in d=3 can be reproduce with the function HardThermal1LoopInt["B",3,0,3]*)

sol=DIMENSION6MATCHING[{Tens2},{2},{\[Alpha][\[Phi]^2 F^2]},d][[1]];  (*The matching can be done singularly for each group tensor*)
sol//Factor//TableForm


ODIM6[2,d]//MatrixForm (*The functions ODIM6 and ODIM5 return the group tensors of the various operators*)


(* ::Section:: *)
(*CROSSCHECK WITH LITERATURE*)

(*The solutions reported by the function DIMENSION6MATCHING are the same as in arxiv: 2503.18904, where the Abelian-Higss model has already been investigated *)
