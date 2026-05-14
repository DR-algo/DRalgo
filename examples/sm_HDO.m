(* ::Package:: *)

Quit[];


SetDirectory[NotebookDirectory[]];
DRalgo`DRalgo`$LoadGroupMath=True;
DRalgo`DRalgo`$GroupMathMultipleModels=True; (*Put this if you want to create multiple model-files with the same kernel*)
<<../Kernel/DRalgo.m


(* ::Chapter:: *)
(*Standard Model*)


(*See 1106.0034 [hep-ph] for a review*)


(* ::Section:: *)
(*Model*)


Group={"SU3","SU2","U1"};
RepAdjoint={{1,1},{2},0};
HiggsDoublet={{{0,0},{1},1/2},"C"};
RepScalar={HiggsDoublet};
CouplingName={gs,gw,gY};


(*Rep1={{{1,0},{1},1/6},"L"};
Rep2={{{1,0},{0},2/3},"R"};
Rep3={{{1,0},{0},-1/3},"R"};
Rep4={{{0,0},{1},-1/2},"L"};
Rep5={{{0,0},{0},-1},"R"};
RepFermion1Gen={Rep1,Rep2,Rep3,Rep4,Rep5};*)


Rep1={{{1,0},{1},Yq/2},"L"};
Rep2={{{1,0},{0},Yu/2},"R"};
Rep3={{{1,0},{0},Yd/2},"R"};
Rep4={{{0,0},{1},Yl/2},"L"};
Rep5={{{0,0},{0},Ye/2},"R"};
RepFermion1Gen={Rep1,Rep2,Rep3,Rep4,Rep5};


RepFermion3Gen={RepFermion1Gen,RepFermion1Gen,RepFermion1Gen}//Flatten[#,1]&;


(* ::Text:: *)
(*The input for the gauge interactions toDRalgo are then given by*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion3Gen,RepScalar];


(* ::Text:: *)
(*The first element is the vector self - interaction matrix :*)


InputInv={{1,1},{True,False}};
MassTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;


VMass=m2*MassTerm1;


\[Mu]ij=GradMass[VMass]//Simplify//SparseArray;


QuarticTerm1=MassTerm1^2;


VQuartic=\[Lambda]1H*QuarticTerm1;


\[Lambda]4=GradQuartic[VQuartic];


InputInv={{1,1,2},{False,False,True}}; 
YukawaDoublet=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv]//Simplify;


Ysff=-GradYukawa[yt*YukawaDoublet[[1]]];


YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt>0}]];


(*Llmat={{0,0,0,0,1,0,0},{0,0,0,0,0,1,0}};
Ermat={0,0,0,0,0,0,1};
Qlmat={{1,0,0,0,0,0,0},{0,1,0,0,0,0,0}};
Urmat={0,0,1,0,0,0,0};
Drmat={0,0,0,1,0,0,0};
V={{1,0,I,0},{0,1,0,I}}/Sqrt[2];
Vs={{1,0,-I,0},{0,1,0,-I}}/Sqrt[2];

Ysff=2SymmetrizeTensor[(yl Contract[V,Llmat,Ermat,{{1,3}}]+yd Contract[V,Qlmat,Drmat,{{1,3}}]+ yu Contract[Vtilde,Qlmat,Urmat,{{1,3}}]),{{2,3}}];
YsffC=Conjugate[Ysff];*)


HYPERCHARGE={Yq->1/3,Yu->4/3,Yd->-2/3,Yl->-1,Ye->-2};


(* ::Section:: *)
(*Dimensional Reduction*)


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False,Mode->3];
PerformDRhard[]


(* ::Section:: *)
(*Dimension-6 Matching*)


(* ::Subsection:: *)
(*Definitions*)


(*
	The Higgs field H is internally decomposed in single real scalar fields
	as follows H=(\phi_1+I\phi_3,\phi_2+I\phi_4)^T,
	which can be compactly be written has H_I=V_{I i}\phi_i with i=1,..,4. 
	All the operator tensors associated to higher dimensional operators
	can be constructed in terms of the matrix V, and its complex conjugate Vs,
	the identity projected on the strong gauge fields G idG,
	the identity projected on the Weak gauge fields W idW,
	the identity projected on the hypercharge gauge field B idB,
	the vector of the field B ProjB,
	the pauli matrices (namely the generator in the fundamental representation of SU(2)),
	the group structure constant of SU(2)_L \[Epsilon]W,
	the group structure constants of SU(3)_c fG, and
	the structure constant d^{abc} of SU(3)
*)

idW = SparseArray[{
   {i_, i_} /; 9 <= i <= 11 -> 1
   }, {12, 12}];

idG = SparseArray[{
   {i_, i_} /; 1 <= i <= 8 -> 1
   }, {12, 12}];
   
idB = SparseArray[{
   {i_, i_} /; 12 <= i <= 12 -> 1
   }, {12, 12}];

ProjB=SparseArray[{
   {i_} /; 12 <= i <= 12 -> 1
   }, {12}];
   
V={{1,0,I,0},{0,1,0,I}}/Sqrt[2];
Vs={{1,0,-I,0},{0,1,0,-I}}/Sqrt[2];

pauli = {
   {{0, 1}, {1, 0}},
   {{0, -I}, {I, 0}},
   {{1, 0}, {0, -1}}
};


Sigma = Table[
   If[9 <= i <= 11,
      pauli[[i - 8]],
      ConstantArray[0, {2, 2}]
   ],
   {i, 12}
];

\[Epsilon]W=Coefficient[gvvv,gw];
fG=Coefficient[gvvv,gs];

HdH=Contract[Vs,V,{{1,3}}];
Hd\[Sigma]H=Contract[Vs,V,Sigma,{{1,6},{3,7}}];

TGfund=Coefficient[gvff,gs];
TGadj=-I fG;
dG=SymmetrizeTensor[Contract[TGfund,TGfund,TGfund,{{2,9},{3,5},{6,8}}],{{1,2,3}}];

(* Gell-Mann Matrices *)
lambda[1] = {{0, 1, 0}, {1, 0, 0}, {0, 0, 0}};
lambda[2] = {{0, -I, 0}, {I, 0, 0}, {0, 0, 0}};
lambda[3] = {{1, 0, 0}, {0, -1, 0}, {0, 0, 0}};
lambda[4] = {{0, 0, 1}, {0, 0, 0}, {1, 0, 0}};
lambda[5] = {{0, 0, -I}, {0, 0, 0}, {I, 0, 0}};
lambda[6] = {{0, 0, 0}, {0, 0, 1}, {0, 1, 0}};
lambda[7] = {{0, 0, 0}, {0, 0, -I}, {0, I, 0}};
lambda[8] = (1/Sqrt[3]) {{1, 0, 0}, {0, 1, 0}, {0, 0, -2}};

(* d^{abc} Tensor*)
dG = ConstantArray[0, {12, 12, 12}];
dG[[1 ;; 8, 1 ;; 8, 1 ;; 8]] = Table[
  (1/4) Tr[(lambda[a] . lambda[b] + lambda[b] . lambda[a]) . lambda[c]],
  {a, 8}, {b, 8}, {c, 8}
];


(* ::Subsection:: *)
(*Individual operator matching*)


(*
	the matching can be performed on individual operators
*)
OT[6,18]=(
	+I \[Alpha][B0 B^2 D]TensorProduct[idB,ProjB]
	+I \[Alpha][B0 W^2 D]TensorProduct[ProjB,idW]
	+I \[Alpha][W0 B W D,1]Transpose[TensorProduct[ProjB,idW],{2,1,3}]
	+I \[Alpha][W0 B W D,2]TensorProduct[idW,ProjB]
	+I \[Alpha][W0 W^2 D]\[Epsilon]W
	+\[Alpha][B0 G^2 D]TensorProduct[ProjB,idG]
	+I \[Alpha][G0 B G D,1]Transpose[TensorProduct[ProjB,idG],{2,1,3}]
	+I \[Alpha][G0 B G D,2]TensorProduct[idG,ProjB]
	+I \[Alpha][G0 G^2 D,1]fG+\[Alpha][G0 G^2 D,2]dG
	);	

WCs = Cases[OT[6,18]//Normal, \[Alpha][__], \[Infinity]];

sol=Dimension6Matching[{OT[6,18]},{18},3][[1]];
Collect[sol,{_Zb,_Zf},Factor]//TableForm


(* ::Subsection:: *)
(*Full matching*)


(*
	below we construct the full list of tensors
*)
X2=Contract[TGadj,TGadj,{{2,6},{3,5}}];
X3=Contract[TGadj,TGadj,{{2,6},{3,5}}];
X4=Contract[TGadj,TGadj,TGadj,TGadj,{{2,12},{3,5},{6,8},{9,11}}];
X6=Contract[TGadj,TGadj,TGadj,TGadj,TGadj,TGadj,{{2,18},{3,5},{6,8},{9,11},{12,14},{15,17}}];

(*
	We define \[Alpha][__] our Wilson coefficients
*)
Tens1=Factorial[3](\[Alpha][W^3]\[Epsilon]W+\[Alpha][G^3]fG);

Tens2=Factorial[2]^2 SymmetrizeTensor[
	+\[Alpha][H^2B^2]TensorProduct[HdH,idB]
	+\[Alpha][H^2W^2]TensorProduct[HdH,idW]
	+\[Alpha][H^2 B W]TensorProduct[Hd\[Sigma]H,ProjB]
	+\[Alpha][H^2G^2]TensorProduct[HdH,idG],{{1,2},{3,4}}
	];

Tens3=Factorial[2]^2 SymmetrizeTensor[
	+\[Alpha][B0^2 B^2]TensorProduct[idB,idB]
	+\[Alpha][W0^2 W^2,1]TensorProduct[idW,idW]
	+\[Alpha][W0^2 W^2,2]Transpose[TensorProduct[idW,idW],{1,3,2,4}]
	+\[Alpha][B0^2 W^2]TensorProduct[idB,idW]
	+\[Alpha][B0 W0 B W]Transpose[TensorProduct[idB,idW],{1,3,2,4}]
	+\[Alpha][W0^2 B^2]TensorProduct[idW,idB]
	+\[Alpha][G0^2 G^2,1]TensorProduct[idG,idG]
	+\[Alpha][G0^2 G^2,2]Transpose[TensorProduct[idG,idG],{1,3,2,4}]
	+\[Alpha][G0^2 G^2,3]Transpose[X4,{1,3,2,4}]
	+\[Alpha][G0^2 G B]TensorProduct[dG,ProjB]
	+\[Alpha][G0 B0 G^2]TensorProduct[ProjB,dG]
	+\[Alpha][G0^2 B^2]TensorProduct[idG,idB]
	+\[Alpha][G0 B0 G B]Transpose[TensorProduct[idG,idB],{1,3,2,4}]
	+\[Alpha][B0^2 G^2]TensorProduct[idB,idG]
	+\[Alpha][G0^2 W^2]TensorProduct[idG,idW]
	+\[Alpha][G0 W0 G W]Transpose[TensorProduct[idG,idW],{1,3,2,4}]
	+\[Alpha][W0^2 G^2]TensorProduct[idW,idG],{{1,2},{3,4}}
	];
						
Tens4=Factorial[2]^2 SymmetrizeTensor[
	+\[Alpha][H^4D^2,1]Re[Transpose[TensorProduct[HdH,HdH],{1,3,2,4}]]
	+\[Alpha][H^4D^2,2]Im[Transpose[TensorProduct[HdH,HdH],{1,3,2,4}]]
	+\[Alpha][H^4D^2,3]Transpose[TensorProduct[HdH,HdH],{3,2,1,4}]
	+\[Alpha][H^4D^2,4]TensorProduct[HdH,HdH],{{1,2},{3,4}}
	];
						
Tens5=Factorial[2]^2 SymmetrizeTensor[
	+\[Alpha][H^2 B0^2 D^2,1]TensorProduct[HdH,idB]
	+\[Alpha][H^2 W0^2 D^2,1]TensorProduct[HdH,idW]
	+\[Alpha][H^2 W0 B0 D^2,1]TensorProduct[Hd\[Sigma]H,ProjB]
	+\[Alpha][H^2 G0^2 D^2,1]TensorProduct[HdH,idG],{{1,2},{3,4}}
	];

Tens6=(
	+\[Alpha][H^2 B0^2 D^2,2]TensorProduct[Re[HdH],idB]
	+\[Alpha][H^2 B0^2 D^2,3]TensorProduct[Im[HdH],idB]
	+\[Alpha][H^2 W0^2 D^2,2]TensorProduct[Re[HdH],idW]
	+\[Alpha][H^2 W0^2 D^2,3]TensorProduct[Im[HdH],idW]
	+\[Alpha][H^2 W0^2 D^2,4]Contract[Re[Hd\[Sigma]H],\[Epsilon]W,{{3,4}}]
	+\[Alpha][H^2 W0^2 D^2,5]Contract[Im[Hd\[Sigma]H],\[Epsilon]W,{{3,4}}]
	+\[Alpha][H^2 B0 W0 D^2,2]TensorProduct[Re[Hd\[Sigma]H],ProjB]
	+\[Alpha][H^2 B0 W0 D^2,3]TensorProduct[Im[Hd\[Sigma]H],ProjB]
	+\[Alpha][H^2 B0 W0 D^2,4]Transpose[TensorProduct[Re[Hd\[Sigma]H],ProjB],{1,2,4,3}]
	+\[Alpha][H^2 B0 W0 D^2,5]Transpose[TensorProduct[Im[Hd\[Sigma]H],ProjB],{1,2,4,3}]
	+\[Alpha][H^2 G0^2 D^2,2]TensorProduct[Re[HdH],idG]
	+\[Alpha][H^2 G0^2 D^2,3]TensorProduct[Im[HdH],idG]
	); 
	   
Tens7=Factorial[2]^2 SymmetrizeTensor[
	+\[Alpha][H^2 B0^2 D^2,4]TensorProduct[HdH,idB]
	+\[Alpha][H^2 W0^2 D^2,6]TensorProduct[HdH,idW]
	+\[Alpha][H^2 W0 B0 D^2,6]TensorProduct[Hd\[Sigma]H,ProjB]
	+\[Alpha][H^2 G0^2 D^2,6]TensorProduct[HdH,idG],{{1,2},{3,4}}
	];

Tens8=Factorial[2]^2 SymmetrizeTensor[
	+\[Alpha][B0^4 D^2]TensorProduct[idB,idB]
	+\[Alpha][B0^2 W0^2 D^2,1]TensorProduct[idB,idW]
	+\[Alpha][B0^2 W0^2 D^2,2]Transpose[TensorProduct[idB,idW],{1,3,2,4}]
	+\[Alpha][B0^2 W0^2 D^2,3]TensorProduct[idW,idB]
	+\[Alpha][W0^4 D^2,1]TensorProduct[idW,idW]
	+\[Alpha][W0^4 D^2,2]Transpose[TensorProduct[idW,idW],{1,3,2,4}]
	+\[Alpha][G0^4 D^2,1]TensorProduct[idG,idG]
	+\[Alpha][G0^4 D^2,2]Transpose[TensorProduct[idG,idG],{1,3,2,4}]
	+\[Alpha][G0^4 D^2,3]Transpose[X4,{1,3,2,4}]
	+\[Alpha][G0^3 B0 D^2,1]TensorProduct[ProjB,dG]
	+\[Alpha][G0^3 B0 D^2,2]TensorProduct[dG,ProjB]
	+\[Alpha][G0^2 B0^2 D^2,1]TensorProduct[idB,idG]
	+\[Alpha][G0^2 B0^2 D^2,2]Transpose[TensorProduct[idG,idB],{1,3,2,4}]
	+\[Alpha][G0^2 B0^2 D^2,3]TensorProduct[idG,idB]
	+\[Alpha][G0^2 W0^2 D^2,1]TensorProduct[idW,idG]
	+\[Alpha][G0^2 W0^2 D^2,2]Transpose[TensorProduct[idG,idW],{1,3,2,4}]
	+\[Alpha][G0^2 W0^2 D^2,3]TensorProduct[idG,idW],{{1,2},{3,4}}
	];

Tens9=Factorial[6] SymmetrizeTensor[
	+\[Alpha][H^6]TensorProduct[HdH,HdH,HdH],{{1,2,3,4,5,6}}
	];

Tens10=Factorial[4]Factorial[2] SymmetrizeTensor[
	+\[Alpha][H^4 B0^2]TensorProduct[HdH,HdH,idB]
	+\[Alpha][H^4 W0^2,1]TensorProduct[HdH,HdH,idW]
	+\[Alpha][H^4 W0^2,2]Transpose[TensorProduct[Hd\[Sigma]H,Hd\[Sigma]H],{1,2,5,4,3,6}]
	+\[Alpha][H^4 B0 W0]TensorProduct[HdH,Hd\[Sigma]H,ProjB]
	+\[Alpha][H^4 G0^2]TensorProduct[HdH,HdH,idG],{{1,2,3,4},{5,6}}
	];
												
Tens11=Factorial[4]Factorial[2] SymmetrizeTensor[
	+\[Alpha][H^2 B0^4]TensorProduct[HdH,idB,idB]
	+\[Alpha][H^2 B0^3 W0]TensorProduct[Hd\[Sigma]H,idB,ProjB]
	+\[Alpha][H^2 B0^2 W0^2]TensorProduct[HdH,idW,idB]
	+\[Alpha][H^2 B0 W0^3]TensorProduct[Hd\[Sigma]H,idW,ProjB]
	+\[Alpha][H^2 W0^4]TensorProduct[HdH,idW,idW]
	+\[Alpha][H^2 G0^2 B0^2]TensorProduct[HdH,idG,idB]
	+\[Alpha][H^2 G0^2 W0^2]TensorProduct[HdH,idG,idW]
	+\[Alpha][H^2 G0^2 W0 B0]TensorProduct[Hd\[Sigma]H,idG,ProjB]
	+\[Alpha][H^2 G0^3 B0]TensorProduct[HdH,dG,ProjB]
	+\[Alpha][H^2 G0^3 W0]TensorProduct[Hd\[Sigma]H,dG]
	+\[Alpha][H^2 G0^4]TensorProduct[HdH,idG,idG],{{1,2},{3,4,5,6}}
	];	
											  
Tens12=Factorial[6]SymmetrizeTensor[
	+\[Alpha][B0^6]TensorProduct[idB,idB,idB]
	+\[Alpha][B0^4 W0^2]TensorProduct[idB,idB,idW]
	+\[Alpha][B0^2 W0^4]TensorProduct[idB,idW,idW]
	+\[Alpha][W0^6]TensorProduct[idW,idW,idW]
	+\[Alpha][G0^2 B0^4]TensorProduct[idG,idB,idB]
	+\[Alpha][G0^2 B0^2 W0^2]TensorProduct[idG,idW,idB]
	+\[Alpha][G0^2 W0^4]TensorProduct[idG,idW,idW]
	+\[Alpha][G0^3 B0^3]TensorProduct[dG,idB,ProjB]
	+\[Alpha][G0^3 B0 W0^2]TensorProduct[dG,ProjB,idW]
	+\[Alpha][G0^4 B0^2]TensorProduct[idG,idG,idB]
	+\[Alpha][G0^4 W0^2]TensorProduct[idG,idG,idW]
	+\[Alpha][G0^5 B0]TensorProduct[dG,idG,ProjB]
	+\[Alpha][G0^6,1]TensorProduct[dG,dG]
	+\[Alpha][G0^6,2]TensorProduct[idG,idG,idG],{{1,2,3,4,5,6}}
	];	
								   
Tens13=Factorial[2]SymmetrizeTensor[
	+\[Alpha][W^2D^2]TensorProduct[idW]
	+\[Alpha][B^2D^2]TensorProduct[idB]
	+\[Alpha][G^2D^2]TensorProduct[idG],{{1,2}}
	];

Tens14=(
	+\[Alpha][H^2 W D^2]Im[Hd\[Sigma]H]
	+\[Alpha][H^2 B D^2]TensorProduct[Im[HdH],ProjB]
	);	

Tens15=Factorial[2]\[Alpha][H^2 D^4]SymmetrizeTensor[HdH,{{1,2}}];

Tens16=(
	+\[Alpha][W0^2 W D^2]\[Epsilon]W
	+Factorial[2]\[Alpha][W0 B0 W D^2](TensorProduct[idW,ProjB]-Transpose[TensorProduct[idW,ProjB],{1,3,2}])
	+\[Alpha][G0^2 G D^2]fG+Factorial[2]\[Alpha][G0 B0 G D^2](TensorProduct[idG,ProjB]-Transpose[TensorProduct[idG,ProjB],{1,3,2}])
	);

Tens17=Factorial[2](
	+\[Alpha][W0^2 D^4]idW
	+\[Alpha][B0^2 D^4]idB
	+\[Alpha][G0^2 D^4]idG
	);												   					  																											  																						   					  																											  																								   					  																											  																																											  																								   					  																											  																						   					  																											  																								   					  																											  																						

Tens18=(
	+I \[Alpha][B0 B^2 D]TensorProduct[idB,ProjB]
	+I \[Alpha][B0 W^2 D]TensorProduct[ProjB,idW]
	+I \[Alpha][W0 B W D,1]Transpose[TensorProduct[ProjB,idW],{2,1,3}]
	+I \[Alpha][W0 B W D,2]TensorProduct[idW,ProjB]
	+I \[Alpha][W0 W^2 D]\[Epsilon]W
	+\[Alpha][B0 G^2 D]TensorProduct[ProjB,idG]
	+I \[Alpha][G0 B G D,1]Transpose[TensorProduct[ProjB,idG],{2,1,3}]
	+I \[Alpha][G0 B G D,2]TensorProduct[idG,ProjB]
	+I \[Alpha][G0 G^2 D,1]fG+\[Alpha][G0 G^2 D,2]dG
	);	

Tens19=(
	+I \[Alpha][H^2 B0 B D,1]TensorProduct[Re[HdH],idB]
	+I \[Alpha][H^2 B0 B D,2]TensorProduct[Im[HdH],idB]
	+I \[Alpha][H^2 B0 W D,1]Transpose[TensorProduct[Re[Hd\[Sigma]H],ProjB],{1,2,4,3}]
	+I \[Alpha][H^2 B0 W D,2]Transpose[TensorProduct[Im[Hd\[Sigma]H],ProjB],{1,2,4,3}]
	+I \[Alpha][H^2 W0 B D,1]TensorProduct[Re[Hd\[Sigma]H],ProjB]
	+I \[Alpha][H^2 W0 B D,2]TensorProduct[Im[Hd\[Sigma]H],ProjB]
	+I \[Alpha][H^2 W0 W D,1]TensorProduct[Re[HdH],idW]
	+I \[Alpha][H^2 W0 W D,2]TensorProduct[Im[HdH],idW]
	+I \[Alpha][H^2 W0 W D,3]Contract[Re[Hd\[Sigma]H],\[Epsilon]W,{{3,4}}]
	+I \[Alpha][H^2 W0 W D,4]Contract[Im[Hd\[Sigma]H],\[Epsilon]W,{{3,4}}]
	+I \[Alpha][H^2 G0 G D,1]TensorProduct[Re[HdH],idG]
	+I \[Alpha][H^2 G0 G D,2]TensorProduct[Im[HdH],idG]
	);

Tens20NS=(
	+I \[Alpha][B0^2 W0 W D]TensorProduct[idB,idW]
	+I \[Alpha][B0 W0^2 B D]Transpose[TensorProduct[idB,idW],{1,4,3,2}]
	+I \[Alpha][B0 W0^2 W D]TensorProduct[ProjB,\[Epsilon]W]
	+I \[Alpha][W0^3 W D]TensorProduct[idW,idW]
	+I \[Alpha][B0^2 G0 G D]TensorProduct[idB,idG]
	+I \[Alpha][B0 G0^2 G D,1]Transpose[TensorProduct[idB,idG],{1,4,3,2}]
	+I \[Alpha][B0 G0^2 G D,2]TensorProduct[ProjB,fG]
	+I \[Alpha][G0^3 G D]TensorProduct[idG,idG]
	+I \[Alpha][B0 G0^2 G D,3]Transpose[TensorProduct[dG,ProjB],{1,2,4,3}]
	+I \[Alpha][W0^2 G0 G D]TensorProduct[idW,idG]
	+I \[Alpha][G0^2 W0 W D]TensorProduct[idG,idW]
	);

Tens20=(
	+(1/2)(
		+Tens20NS+
		Transpose[Tens20NS,{2,1,3,4}]
		)
	-(1/6)(
		+Tens20NS
		+Transpose[Tens20NS,{2,1,3,4}]
		+Transpose[Tens20NS,{2,3,1,4}]
		+Transpose[Tens20NS,{3,2,1,4}]
		+Transpose[Tens20NS,{1,3,2,4}]
		+Transpose[Tens20NS,{3,1,2,4}]
		)
	);


TensorList={
	Tens1,Tens2,Tens3,Tens4,Tens5,Tens6,Tens7,Tens8,Tens9,Tens10,
	Tens11,Tens12,Tens13,Tens14,Tens15,Tens16,Tens17,Tens18,Tens19,Tens20
	};
NList={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
WC = DeleteDuplicates[Cases[TensorList//Normal, \[Alpha][__], \[Infinity]]];

(*
	The matching of operators with 6 fields, namely
	phi^6, phi^4 B0^2, phi^2 B0^4, and B0^6,
	requires a significant amount of running time.
	It is therefore convenient to perform the matching of these operators separately.
*)
n={1,2,3,4,5,6,7,8,13,14,15,16,17,18,19,20};
sol=Dimension6Matching[TensorList[[n]],NList[[n]],WC,3][[1]];
Collect[sol,{_Zb,_Zf},Factor]/.HYPERCHARGE//TableForm


n={9,10,11,12};

sol6=Dimension6Matching[TensorList[[n]],NList[[n]],WC,3][[1]];//AbsoluteTiming
Collect[sol6,{_Zb,_Zf},Factor]/.HYPERCHARGE//TableForm


solTex=Collect[Join[sol,sol6]/.HYPERCHARGE/.{\[Lambda]1H->l4,\[Xi]->xi},{_Zb,_Zf},Factor]//TableForm


(* ::Section:: *)
(*COMPARISON WITH LITERATURE*)


(*
	The Wilson coefficients related to the electro-weak sector only
	agree with what is reported in arxiv 2503.20016
*)
