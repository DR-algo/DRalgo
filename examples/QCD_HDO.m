(* ::Package:: *)

Quit[];


SetDirectory[NotebookDirectory[]];
DRalgo`DRalgo`$LoadGroupMath=True;
<<../Kernel/DRalgo.m


(* ::Chapter:: *)
(*QCD*)


(* ::Section::Closed:: *)
(*Model*)


Group={"SU3"};
CouplingName={gs};
RepAdjoint={{1,1}};
RepScalar={};
RepFermionL={{{1,0}},"L"};
RepFermionR={{{1,0}},"R"};
RepFermion={RepFermionL,RepFermionR};


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];


(* ::Section::Closed:: *)
(*Dimensional Reduction*)


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PerformDRhard[]


(* ::Section:: *)
(*Dimension 6 Matching*)

(*Here we construct tensors associated to the group structure constants and its contraction*)

Tadj=-I gvvv/gs;
X2=Contract[Tadj,Tadj,{{2,6},{3,5}}];
X3=Contract[Tadj,Tadj,Tadj,{{2,9},{3,5},{6,8}}];
X4=Contract[Tadj,Tadj,Tadj,Tadj,{{2,12},{3,5},{6,8},{9,11}}];
X6=Contract[Tadj,Tadj,Tadj,Tadj,Tadj,Tadj,{{2,18},{3,5},{6,8},{9,11},{12,14},{15,17}}];

(*Below we construct the group tensors of the higher dimensional operators, having defined \[Alpha][__] our Wilson Coefficients*)

OT[6,1]=\[Alpha][F^3] X3;
OT[6,3]=SymmetrizeTensor[\[Alpha][A^2F^2,1]X4+\[Alpha][A^2F^2,2]Transpose[X4,{1,3,2,4}]+\[Alpha][A^2F^2,3]TensorProduct[X2,X2],{{1,2},{3,4}}];
OT[6,8]=SymmetrizeTensor[\[Alpha][D^2A^4,1]X4+\[Alpha][D^2A^4,2]Transpose[X4,{1,3,2,4}]+\[Alpha][D^2A^4,3]TensorProduct[X2,X2],{{1,2},{3,4}}];
OT[6,12]=SymmetrizeTensor[\[Alpha][A^6,1]X6+\[Alpha][A^6,2]TensorProduct[X4,X2],{{1,2,3,4,5,6}}];
OT[6,13]=\[Alpha][D^2F^2]X2;
OT[6,16]=\[Alpha][D^2 A^2 F]X3;
OT[6,17]=\[Alpha][D^4A^2]X2;


TensorList={OT[6,1],OT[6,3],OT[6,8],OT[6,12],OT[6,13],OT[6,16],OT[6,17]};
NList={1,3,8,12,13,16,17};
WCs=DeleteDuplicates[Cases[TensorList//Normal, \[Alpha][__], \[Infinity]]]; (*WC is the array of the various Wilson Coefficients \[Alpha][...]*)

sol6=Dimension6Matching[TensorList,NList,WCs,3][[1]]; (*Dimension6Matching and Dimension5Matching find the values of the Wilson coefficients listed in WC, d is the number of spatial dimensions, Zb and Zf are 1 loop master integral *)
Collect[sol6,{_Zb,_Zf},Factor]//TableForm

(*The matching of the tensor Tens12 takes approximately 10 minutes to complete, as it entails handling an 8\[Times]8\[Times]8\[Times]8\[Times]8\[Times]8 tensor. It is thus advantageous to first perform 
the matching for the remaining tensors and subsequently address the A0^6 operator separately. This procedure can be implemented as follows*)

TensorListWithoutA6={OT[6,1],OT[6,3],OT[6,8],OT[6,13],OT[6,16],OT[6,17]};
NListWithoutA6={1,3,8,13,16,17};
WCsWithoutA6=DeleteDuplicates[Cases[TensorListWithoutA6//Normal, \[Alpha][__], \[Infinity]]]; (*WC is the array of the various Wilson Coefficients \[Alpha][...]*)

solWithoutA6=Dimension6Matching[TensorListWithoutA6,NListWithoutA6,WCsWithoutA6,3][[1]]; (*Dimension6Matching and Dimension5Matching find the values of the Wilson coefficients listed in WC, d is the number of spatial dimensions, Zb and Zf are 1 loop master integral *)
Collect[solWithoutA6,{_Zb,_Zf},Factor]//TableForm


solA6=Dimension6Matching[{OT[12]},{12},{\[Alpha][A^6,1],\[Alpha][A^6,2]},3][[1]];  (*The matching can be done singularly for each group tensor*)
solA6//Factor//TableForm



(*The matching can be performed on individual operators*)
sol=Dimension6Matching[{OT[6,1]},{1},{\[Alpha][F^3]},3][[1]];
Collect[sol,{_Zb,_Zf},Factor]//TableForm

ODIM6[1,3]//MatrixForm (*The functions ODIM6 and ODIM5 return the group tensors of the various operators*)
ODIM6[17,3]//MatrixForm (*The functions ODIM6 and ODIM5 return the group tensors of the various operators*)


(* ::Section:: *)
(*FIELD REDEFINITIONS*)

(*The current version of Higher Dimensional Operators routines desn't include field redefinitions*)
(*Through field redefinition is it possible to check gauge dependence cancellation*)

gE=gs Sqrt[T];
solDRalgo=Join[solA6,solWithoutA6];

\[Alpha][F^3]/.solDRalgo//Factor
\[Alpha][A^2F^2,1]/.solDRalgo//Factor
\[Alpha][A^2F^2,2]/.solDRalgo//Factor
(\[Alpha][D^2A^4,2]-1/2 \[Alpha][D^2A^4,1]-6 I gE \[Alpha][D^2A^2F]-6gE^2 \[Alpha][D^2F^2])/.solDRalgo//Factor
\[Alpha][A^6,1]/.solDRalgo//Factor
\[Alpha][A^6,2]/.solDRalgo//Factor
