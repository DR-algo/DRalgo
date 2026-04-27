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
RepFermion={};


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];


(* ::Section::Closed:: *)
(*Dimensional Reduction*)


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PerformDRhard[]


(* ::Section:: *)
(*Dimension 6 Matching*)


fvvv=gvvv/gs;
Tvvv=-I*fvvv;
id=IdentityMatrix[8];
X2=Contract[Tvvv,Tvvv,{{2,6},{3,5}}];
X3=Contract[Tvvv,Tvvv,Tvvv,{{2,9},{3,5},{6,8}}];
X4=Contract[Tvvv,Tvvv,Tvvv,Tvvv,{{2,12},{3,5},{6,8},{9,11}}];
X6=Contract[Tvvv,Tvvv,Tvvv,Tvvv,Tvvv,Tvvv,{{2,18},{3,5},{6,8},{9,11},{12,14},{15,17}}];

gE=gs;


(*Here we construct the group tensors of the higher dimensional operators*)

Tens1=Factorial[3]I \[Alpha][F^3] gE^3 X3;
Tens3=SymmetrizeTensor[4(\[Alpha][A^2F^2,1]X4+\[Alpha][A^2F^2,2]Transpose[X4,{1,3,2,4}]),{{1,2},{3,4}}] gE^4;
Tens8=SymmetrizeTensor[4(\[Alpha][D^2A^4,1]X4+\[Alpha][D^2A^4,2]Transpose[X4,{1,3,2,4}]),{{1,2},{3,4}}] gE^4;
Tens12=(\[Alpha][A^6,1]SymmetrizeTensor[X6,{{1,2,3,4,5,6}}]+\[Alpha][A^6,2]SymmetrizeTensor[TensorProduct[X4,X2],{{1,2,3,4,5,6}}])Factorial[6] gE^6 ;
Tens13=2 \[Alpha][D^2F^2]X2 gE^2 ;
Tens16=2 I \[Alpha][D^2 A^2 F]X3 gE^3;
Tens17=2 \[Alpha][D^4A^2]X2 gE^2;


TensorList={Tens1,Tens3,Tens8,Tens12,Tens13,Tens16,Tens17}; (*TensorList is the array of the various group tensors refered to higher dimensional operators*)
NList={1,3,8,12,13,16,17}; (*NList is the array that indicate to which operator the group tensors in TensorList refer to*)
WCs=DeleteDuplicates[Cases[TensorList//Normal, \[Alpha][__], \[Infinity]]]; (*WC is the array of the various Wilson Coefficients \[Alpha][...]*)


soltot=DIMENSION6MATCHING[TensorList,NList,WCs,3][[1]]; (*DIMENSION6MATCHING and DIMENSION5MATCHING find the values of the Wilson coefficients listed in WC, d is the number of spatial dimensions, Zb and Zf are 1 loop master integral *)
soltot//Factor//TableForm


sol=DIMENSION6MATCHING[{Tens1},{1},{\[Alpha][F^3]},d][[1]];  (*The matching can be done singularly for each group tensor*)
sol//Factor//TableForm


ODIM6[1,d]//MatrixForm (*The functions ODIM6 and ODIM5 return the group tensors of the various operators*)


(* ::Section:: *)
(*FIELD REDEFINITIONS*)


SOST1={\[Alpha][D^2F^2]->c[1],\[Alpha][D^4A^2]->c[1]+c[2],\[Alpha][F^3]->c[3],\[Alpha][D^2A^2F]->+2c[1]-3/2 c[3]-1/2 c[4]+1/2 c[5],\[Alpha][A^2F^2,1]->3/2 c[3]+1/2 c[4]+c[6],\[Alpha][A^2F^2,2]->-3/2 c[3]-1/2 c[4]+c[7],\[Alpha][D^2A^4,1]->2c[1]+c[5]+2c[6]+c[8],\[Alpha][D^2A^4,2]->-2c[1]-c[5]+2c[7]+c[9],\[Alpha][A^6,1]->c[10],\[Alpha][A^6,2]->c[11]};
solRB={c[1]->2(41-d)/120+2(8-\[Alpha])\[Alpha]/48,c[2]->-2(d-1)(5-d)/(120)+2(d-5)(4+\[Alpha])\[Alpha]/48,c[3]->2(1-d)/180,c[4]->c[5]-2(d-1)(d-5)/60-(d-5)\[Alpha]/6,c[5]->2c[7]+2(d-21)(d-5)/30+(d-5)\[Alpha]/6,c[6]->-c[7]-2(d-5)(d-25)/24,c[8]->2(d-1)(d-3)(d-5)/20+(d-5)(d-3)\[Alpha]/3,c[9]->2(d-1)(d-3)(d-5)/30+(d-5)(d-3)\[Alpha]/6,c[10]->2(d-1)^2 (d-3)(d-5)/180,c[11]->0,c[7]->2(21-d)(5-d)/60};


(*The current version of Higher Dimensional Operators routines desn't include field redefinitions*)


(*Through field redefinition is it possible to check gauge dependence cancellation*)


(*Wih field redefinition, is it possible to compare with results already existing in literature, for this model with arXiv: 1803.08689*)


gE=gs Sqrt[T];
T=1;

solDRalgo=soltot; (*{\[Alpha][F^3]\[Rule]-(1/90) (-1+d) Sqrt[T] Zb[3,0],\[Alpha][A^2 F^2,1]\[Rule]1/60 (-419 Zb[3,0]+103 d Zb[3,0]-4 d^2 Zb[3,0]),\[Alpha][A^2 F^2,2]\[Rule]1/60 (-206 Zb[3,0]+47 d Zb[3,0]-d^2 Zb[3,0]),\[Alpha][A^4 D^2,1]\[Rule]1/60 T (-303+32 d-20 d^2+6 d^3-530 \[Xi]+300 d \[Xi]-40 d^2 \[Xi]-5 \[Xi]^2) Zb[3,0],\[Alpha][A^4 D^2,2]\[Rule]1/60 (-197 T Zb[3,0]+18 d T Zb[3,0]-20 d^2 T Zb[3,0]+4 d^3 T Zb[3,0]-370 T \[Xi] Zb[3,0]+180 d T \[Xi] Zb[3,0]-20 d^2 T \[Xi] Zb[3,0]+5 T \[Xi]^2 Zb[3,0]),\[Alpha][A^6,1]\[Rule]1/90 (-5+d) (-3+d) (-1+d)^2 T^2 Zb[3,0],\[Alpha][A^6,2]\[Rule]0,\[Alpha][D^2 F^2]\[Rule]-(1/120) (-117+2 d+30 \[Xi]+5 \[Xi]^2) Zb[3,0],\[Alpha][A^2 D^2 F]\[Rule]1/60 Sqrt[T] (71+3 d+d^2+20 \[Xi]-10 d \[Xi]-5 \[Xi]^2) Zb[3,0],\[Alpha][A^2 D^4]\[Rule]1/120 (2+11 d+2 d^2+120 \[Xi]-30 d \[Xi]-30 \[Xi]^2+5 d \[Xi]^2) Zb[3,0]}*) 
solLit=Join[SOST1,solRB];


\[Alpha][F^3]/.solDRalgo//Factor
\[Alpha][F^3]//.solLit//Factor


\[Alpha][A^2F^2,1]/.solDRalgo//Factor
\[Alpha][A^2F^2,1]//.solLit//Factor


\[Alpha][A^2F^2,2]/.solDRalgo//Factor
\[Alpha][A^2F^2,2]//.solLit//Factor


(\[Alpha][D^2A^4,2]-1/2 \[Alpha][D^2A^4,1]+3\[Alpha][D^2A^2F]-3\[Alpha][D^2F^2])/.solDRalgo//Factor
(\[Alpha][D^2A^4,2]-1/2 \[Alpha][D^2A^4,1]+3\[Alpha][D^2A^2F]-3\[Alpha][D^2F^2])//.solLit//Factor


\[Alpha][A^6,1]/.solDRalgo//Factor
\[Alpha][A^6,1]//.solLit


\[Alpha][A^6,2]/.solDRalgo//Factor
\[Alpha][A^6,2]//.solLit

