(* ::Package:: *)

Quit[];


SetDirectory[NotebookDirectory[]];
DRalgo`DRalgo`$LoadGroupMath=True;
<<../Kernel/DRalgo.m


(* ::Chapter:: *)
(*Standard Model*)


(*See 1106.0034 [hep-ph] for a review*)


(* ::Section::Closed:: *)
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


IPERCHARGE={Yq->1/3,Yu->4/3,Yd->-2/3,Yl->-1,Ye->-2};


(* ::Section:: *)
(*Dimensional Reduction*)


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False,Mode->3];
PerformDRhard[]


(* ::Section:: *)
(*Dimension-6 Matching*)


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

X2=Contract[TGadj,TGadj,{{2,6},{3,5}}];
X3=Contract[TGadj,TGadj,{{2,6},{3,5}}];
X4=Contract[TGadj,TGadj,TGadj,TGadj,{{2,12},{3,5},{6,8},{9,11}}];
X6=Contract[TGadj,TGadj,TGadj,TGadj,TGadj,TGadj,{{2,18},{3,5},{6,8},{9,11},{12,14},{15,17}}];


Tens1=Factorial[3](\[Alpha][W^3]\[Epsilon]W+\[Alpha][G^3]fG);

Tens2=Factorial[2]^2 SymmetrizeTensor[\[Alpha][H^2B^2]TensorProduct[HdH,idB]+\[Alpha][H^2W^2]TensorProduct[HdH,idW]+\[Alpha][H^2 B W]TensorProduct[Hd\[Sigma]H,ProjB]+\[Alpha][H^2G^2]TensorProduct[HdH,idG],{{1,2},{3,4}}];

Tens3=Factorial[2]^2 SymmetrizeTensor[\[Alpha][B0^2 B^2]TensorProduct[idB,idB]+\[Alpha][W0^2 W^2,1]TensorProduct[idW,idW]+\[Alpha][W0^2 W^2,2]Transpose[TensorProduct[idW,idW],{1,3,2,4}]+\
						\[Alpha][B0^2 W^2]TensorProduct[idB,idW]+\[Alpha][B0 W0 B W]Transpose[TensorProduct[idB,idW],{1,3,2,4}]+\[Alpha][W0^2 B^2]TensorProduct[idW,idB]+\
						\[Alpha][G0^2 G^2,1]TensorProduct[idG,idG]+\[Alpha][G0^2 G^2,2]Transpose[TensorProduct[idG,idG],{1,3,2,4}]+\[Alpha][G0^2 G^2,3]Transpose[X4,{1,3,2,4}]+\[Alpha][G0^2 G B]TensorProduct[dG,ProjB]+\
						\[Alpha][G0 B0 G^2]TensorProduct[ProjB,dG]+\[Alpha][G0^2 B^2]TensorProduct[idG,idB]+\[Alpha][G0 B0 G B]Transpose[TensorProduct[idG,idB],{1,3,2,4}]+\
						\[Alpha][B0^2 G^2]TensorProduct[idB,idG]+\[Alpha][G0^2 W^2]TensorProduct[idG,idW]+\[Alpha][G0 W0 G W]Transpose[TensorProduct[idG,idW],{1,3,2,4}]+\[Alpha][W0^2 G^2]TensorProduct[idW,idG],{{1,2},{3,4}}];
						
Tens4=Factorial[2]^2 SymmetrizeTensor[\[Alpha][H^4D^2,1]Re[Transpose[TensorProduct[HdH,HdH],{1,3,2,4}]]+\[Alpha][H^4D^2,2]Im[Transpose[TensorProduct[HdH,HdH],{1,3,2,4}]]+\[Alpha][H^4D^2,3]Transpose[TensorProduct[HdH,HdH],{3,2,1,4}]+\
					 \[Alpha][H^4D^2,4]TensorProduct[HdH,HdH],{{1,2},{3,4}}];
						
Tens5=Factorial[2]^2 SymmetrizeTensor[\[Alpha][H^2 B0^2 D^2,1]TensorProduct[HdH,idB]+\[Alpha][H^2 W0^2 D^2,1]TensorProduct[HdH,idW]+\[Alpha][H^2 W0 B0 D^2,1]TensorProduct[Hd\[Sigma]H,ProjB]+\[Alpha][H^2 G0^2 D^2,1]TensorProduct[HdH,idG],{{1,2},{3,4}}];

Tens6=\[Alpha][H^2 B0^2 D^2,2]TensorProduct[Re[HdH],idB]+\[Alpha][H^2 B0^2 D^2,3]TensorProduct[Im[HdH],idB]+\[Alpha][H^2 W0^2 D^2,2]TensorProduct[Re[HdH],idW]+\
	  \[Alpha][H^2 W0^2 D^2,3]TensorProduct[Im[HdH],idW]+\[Alpha][H^2 W0^2 D^2,4]Contract[Re[Hd\[Sigma]H],\[Epsilon]W,{{3,4}}]+\[Alpha][H^2 W0^2 D^2,5]Contract[Im[Hd\[Sigma]H],\[Epsilon]W,{{3,4}}]+\
	  \[Alpha][H^2 B0 W0 D^2,2]TensorProduct[Re[Hd\[Sigma]H],ProjB]+\[Alpha][H^2 B0 W0 D^2,3]TensorProduct[Im[Hd\[Sigma]H],ProjB]+\[Alpha][H^2 B0 W0 D^2,4]Transpose[TensorProduct[Re[Hd\[Sigma]H],ProjB],{1,2,4,3}]+\
	  \[Alpha][H^2 B0 W0 D^2,5]Transpose[TensorProduct[Im[Hd\[Sigma]H],ProjB],{1,2,4,3}]+\[Alpha][H^2 G0^2 D^2,2]TensorProduct[Re[HdH],idG]+\[Alpha][H^2 G0^2 D^2,3]TensorProduct[Im[HdH],idG]; 
	   
Tens7=Factorial[2]^2 SymmetrizeTensor[\[Alpha][H^2 B0^2 D^2,4]TensorProduct[HdH,idB]+\[Alpha][H^2 W0^2 D^2,6]TensorProduct[HdH,idW]+\[Alpha][H^2 W0 B0 D^2,6]TensorProduct[Hd\[Sigma]H,ProjB]+\[Alpha][H^2 G0^2 D^2,6]TensorProduct[HdH,idG],{{1,2},{3,4}}];

Tens8=Factorial[2]^2 SymmetrizeTensor[\[Alpha][B0^4 D^2]TensorProduct[idB,idB]+\[Alpha][B0^2 W0^2 D^2,1]TensorProduct[idB,idW]+\[Alpha][B0^2 W0^2 D^2,2]Transpose[TensorProduct[idB,idW],{1,3,2,4}]+\
									\[Alpha][B0^2 W0^2 D^2,3]TensorProduct[idW,idB]+\[Alpha][W0^4 D^2,1]TensorProduct[idW,idW]+\[Alpha][W0^4 D^2,2]Transpose[TensorProduct[idW,idW],{1,3,2,4}]+\
									\[Alpha][G0^4 D^2,1]TensorProduct[idG,idG]+\[Alpha][G0^4 D^2,2]Transpose[TensorProduct[idG,idG],{1,3,2,4}]+\[Alpha][G0^4 D^2,3]Transpose[X4,{1,3,2,4}]+\[Alpha][G0^3 B0 D^2,1]TensorProduct[ProjB,dG]+\
									\[Alpha][G0^3 B0 D^2,2]TensorProduct[dG,ProjB]+\[Alpha][G0^2 B0^2 D^2,1]TensorProduct[idB,idG]+\[Alpha][G0^2 B0^2 D^2,2]Transpose[TensorProduct[idG,idB],{1,3,2,4}]+\
									\[Alpha][G0^2 B0^2 D^2,3]TensorProduct[idG,idB]+\[Alpha][G0^2 W0^2 D^2,1]TensorProduct[idW,idG]+\[Alpha][G0^2 W0^2 D^2,2]Transpose[TensorProduct[idG,idW],{1,3,2,4}]+\
									\[Alpha][G0^2 W0^2 D^2,3]TensorProduct[idG,idW],{{1,2},{3,4}}];

Tens9=Factorial[6] SymmetrizeTensor[\[Alpha][H^6]TensorProduct[HdH,HdH,HdH],{{1,2,3,4,5,6}}];

Tens10=Factorial[4]Factorial[2] SymmetrizeTensor[\[Alpha][H^4 B0^2]TensorProduct[HdH,HdH,idB]+\[Alpha][H^4 W0^2,1]TensorProduct[HdH,HdH,idW]+\[Alpha][H^4 W0^2,2]Transpose[TensorProduct[Hd\[Sigma]H,Hd\[Sigma]H],{1,2,5,4,3,6}]+\
												\[Alpha][H^4 B0 W0]TensorProduct[HdH,Hd\[Sigma]H,ProjB]+\[Alpha][H^4 G0^2]TensorProduct[HdH,HdH,idG],{{1,2,3,4},{5,6}}];
												
Tens11=Factorial[4]Factorial[2] SymmetrizeTensor[\[Alpha][H^2 B0^4]TensorProduct[HdH,idB,idB]+\[Alpha][H^2 B0^3 W0]TensorProduct[Hd\[Sigma]H,idB,ProjB]+\[Alpha][H^2 B0^2 W0^2]TensorProduct[HdH,idW,idB]+\
											    \[Alpha][H^2 B0 W0^3]TensorProduct[Hd\[Sigma]H,idW,ProjB]+\[Alpha][H^2 W0^4]TensorProduct[HdH,idW,idW]+\[Alpha][H^2 G0^2 B0^2]TensorProduct[HdH,idG,idB]+\
											    \[Alpha][H^2 G0^2 W0^2]TensorProduct[HdH,idG,idW]+\[Alpha][H^2 G0^2 W0 B0]TensorProduct[Hd\[Sigma]H,idG,ProjB]+\[Alpha][H^2 G0^3 B0]TensorProduct[HdH,dG,ProjB]+\
											    \[Alpha][H^2 G0^3 W0]TensorProduct[Hd\[Sigma]H,dG]+\[Alpha][H^2 G0^4]TensorProduct[HdH,idG,idG],{{1,2},{3,4,5,6}}];	
											  
Tens12=Factorial[6]SymmetrizeTensor[\[Alpha][B0^6]TensorProduct[idB,idB,idB]+\[Alpha][B0^4 W0^2]TensorProduct[idB,idB,idW]+\[Alpha][B0^2 W0^4]TensorProduct[idB,idW,idW]+\[Alpha][W0^6]TensorProduct[idW,idW,idW]+\
								   \[Alpha][G0^2 B0^4]TensorProduct[idG,idB,idB]+\[Alpha][G0^2 B0^2 W0^2]TensorProduct[idG,idW,idB]+\[Alpha][G0^2 W0^4]TensorProduct[idG,idW,idW]+\
								   \[Alpha][G0^3 B0^3]TensorProduct[dG,idB,ProjB]+\[Alpha][G0^3 B0 W0^2]TensorProduct[dG,ProjB,idW]+\[Alpha][G0^4 B0^2]TensorProduct[idG,idG,idB]+\
								   \[Alpha][G0^4 W0^2]TensorProduct[idG,idG,idW]+\[Alpha][G0^5 B0]TensorProduct[dG,idG,ProjB]+\[Alpha][G0^6,1]TensorProduct[dG,dG]+\[Alpha][G0^6,2]TensorProduct[idG,idG,idG],{{1,2,3,4,5,6}}];	
								   
Tens13=Factorial[2]SymmetrizeTensor[\[Alpha][W^2D^2]TensorProduct[idW]+\[Alpha][B^2D^2]TensorProduct[idB]+\[Alpha][G^2D^2]TensorProduct[idG],{{1,2}}];

Tens14=\[Alpha][H^2 W D^2]Im[Hd\[Sigma]H]+\[Alpha][H^2 B D^2]TensorProduct[Im[HdH],ProjB];	

Tens15=Factorial[2]\[Alpha][H^2 D^4]SymmetrizeTensor[HdH,{{1,2}}];

Tens16=\[Alpha][W0^2 W D^2]\[Epsilon]W+Factorial[2]\[Alpha][W0 B0 W D^2](TensorProduct[idW,ProjB]-Transpose[TensorProduct[idW,ProjB],{1,3,2}])+\
	   \[Alpha][G0^2 G D^2]fG+Factorial[2]\[Alpha][G0 B0 G D^2](TensorProduct[idG,ProjB]-Transpose[TensorProduct[idG,ProjB],{1,3,2}]);

Tens17=Factorial[2](\[Alpha][W0^2 D^4]idW+\[Alpha][B0^2 D^4]idB+\[Alpha][G0^2 D^4]idG);												   					  																											  																						   					  																											  																								   					  																											  																																											  																								   					  																											  																						   					  																											  																								   					  																											  																						

Tens18=I \[Alpha][B0 B^2 D]TensorProduct[idB,ProjB]+I \[Alpha][B0 W^2 D]TensorProduct[ProjB,idW]+I \[Alpha][W0 B W D,1]Transpose[TensorProduct[ProjB,idW],{2,1,3}]+\
		I \[Alpha][W0 B W D,2]TensorProduct[idW,ProjB]+I \[Alpha][W0 W^2 D]\[Epsilon]W+\[Alpha][B0 G^2 D]TensorProduct[ProjB,idG]+I \[Alpha][G0 B G D,1]Transpose[TensorProduct[ProjB,idG],{2,1,3}]+\
		I \[Alpha][G0 B G D,2]TensorProduct[idG,ProjB]+I \[Alpha][G0 G^2 D,1]fG+\[Alpha][G0 G^2 D,2]dG;	

Tens19=I \[Alpha][H^2 B0 B D,1]TensorProduct[Re[HdH],idB]+I \[Alpha][H^2 B0 B D,2]TensorProduct[Im[HdH],idB]+I \[Alpha][H^2 B0 W D,1]Transpose[TensorProduct[Re[Hd\[Sigma]H],ProjB],{1,2,4,3}]+\
	   I \[Alpha][H^2 B0 W D,2]Transpose[TensorProduct[Im[Hd\[Sigma]H],ProjB],{1,2,4,3}]+I \[Alpha][H^2 W0 B D,1]TensorProduct[Re[Hd\[Sigma]H],ProjB]+I \[Alpha][H^2 W0 B D,2]TensorProduct[Im[Hd\[Sigma]H],ProjB]+\
	   I \[Alpha][H^2 W0 W D,1]TensorProduct[Re[HdH],idW]+I \[Alpha][H^2 W0 W D,2]TensorProduct[Im[HdH],idW]+I \[Alpha][H^2 W0 W D,3]Contract[Re[Hd\[Sigma]H],\[Epsilon]W,{{3,4}}]+\
	   I \[Alpha][H^2 W0 W D,4]Contract[Im[Hd\[Sigma]H],\[Epsilon]W,{{3,4}}]+I \[Alpha][H^2 G0 G D,1]TensorProduct[Re[HdH],idG]+I \[Alpha][H^2 G0 G D,2]TensorProduct[Im[HdH],idG];

Tens20NS=I \[Alpha][B0^2 W0 W D]TensorProduct[idB,idW]+I \[Alpha][B0 W0^2 B D]Transpose[TensorProduct[idB,idW],{1,4,3,2}]+I \[Alpha][B0 W0^2 W D]TensorProduct[ProjB,\[Epsilon]W]+\
	   I \[Alpha][W0^3 W D]TensorProduct[idW,idW]+I \[Alpha][B0^2 G0 G D]TensorProduct[idB,idG]+I \[Alpha][B0 G0^2 G D,1]Transpose[TensorProduct[idB,idG],{1,4,3,2}]+\
	   I \[Alpha][B0 G0^2 G D,2]TensorProduct[ProjB,fG]+I \[Alpha][G0^3 G D]TensorProduct[idG,idG]+I \[Alpha][B0 G0^2 G D,3]Transpose[TensorProduct[dG,ProjB],{1,2,4,3}]+\
	   I \[Alpha][W0^2 G0 G D]TensorProduct[idW,idG]+I \[Alpha][G0^2 W0 W D]TensorProduct[idG,idW];

Tens20=(1/2)(Tens20NS+Transpose[Tens20NS,{2,1,3,4}])-(1/6)(Tens20NS+Transpose[Tens20NS,{2,1,3,4}]+Transpose[Tens20NS,{2,3,1,4}]+Transpose[Tens20NS,{3,2,1,4}]+Transpose[Tens20NS,{1,3,2,4}]+Transpose[Tens20NS,{3,1,2,4}]);


TensorList={Tens1,Tens2,Tens3,Tens4,Tens5,Tens6,Tens7,Tens8,Tens9,Tens10,Tens11,Tens12,Tens13,Tens14,Tens15,Tens16,Tens17,Tens18,Tens19,Tens20};
NList={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
WC = DeleteDuplicates[Cases[TensorList//Normal, \[Alpha][__], \[Infinity]]];


n={1,2,3,4,5,6,7,8,13,14,15,16,17,18,19,20};
sol=DIMENSION6MATCHING[TensorList[[n]],NList[[n]],WC,3][[1]];
Collect[sol,{_Zb,_Zf},Factor]/.IPERCHARGE//TableForm


n={9,10,11,12};

sol6=DIMENSION6MATCHING[TensorList[[n]],NList[[n]],WC,3][[1]];
Collect[sol6,{_Zb,_Zf},Factor]/.IPERCHARGE//TableForm


solTex=Collect[Join[sol,sol6]/.IPERCHARGE/.{\[Lambda]1H->l4,\[Xi]->xi},{_Zb,_Zf},Factor]//TableForm


(* ::Section::Closed:: *)
(*SOLUTIONS MIKI*)


solMFer={mH2 -> -2 (Abs[yl]^2 + Nc yd Conjugate[yd] + Nc yu Conjugate[yu]) I4[
    1], KH -> (Abs[yl]^2 + Nc yd Conjugate[yd] + 
     Nc yu Conjugate[yu]) I4[2], 
 mB0 -> -(1/3) (-1 + d) g1^2 (27 + 11 NC) I4[1], 
 KB0 -> 1/18 (-1 + d) g1^2 (27 + 11 NC) I4[2], 
 KB -> 1/9 g1^2 (27 + 11 NC) I4[2], 
 lmbdH2B02 -> 
  1/36 g1^2 ((-6 + 5 d) Nc yd Conjugate[yd] I4[2] + 
     9 (-14 + 5 d) yl Conjugate[yl] I4[
       2] + (-42 + 17 d) Nc yu Conjugate[yu] I4[2] - 
     24 (-1 + d) I4[
       1] (3 (Tr[aHe] + Tr[aHl1]) + 
        NC (Tr[aHd] - Tr[aHq1] - 2 Tr[aHu]))), 
 lmbdB04 -> ((-3 + d) (-1 + d) g1^4 (729 + 137 NC) I4[2])/1296, 
 lmbdH4 -> 
  2 (aeH Conjugate[yl] + yl Conjugate[aeH] + 
      NC (adH Conjugate[yd] + auH Conjugate[yu] + 
         yd Conjugate[adH] + yu Conjugate[auH])) I4[
     1] + (Abs[yl]^4 + Nc (Abs[yd]^4 + Abs[yu]^4)) I4[2], 
 aH4D21 -> 
  1/2 (yl Conjugate[
        aeH] + (aeH + 2 (-aHe + aHl1 + aHl3) yl) Conjugate[yl] + 
      NC (yd Conjugate[adH] + 
         yu Conjugate[
           auH] + (adH + 2 (-aHd + aHq1 + aHq3) yd) Conjugate[
           yd] + (auH + 2 (-aHq1 + aHq3 + aHu) yu) Conjugate[yu])) I4[
     2] + 1/3 (Abs[yl]^4 + Nc (Abs[yd]^4 + Abs[yu]^4)) I4[3], 
 aH4D22 -> -2 (NC (2 aHd yd - 2 aHq1 yd - aHud yu) Conjugate[yd] + 
      2 (aHe - aHl1) yl Conjugate[yl] + 
      NC (2 (aHq1 - aHu) yu - yd Conjugate[aHud]) Conjugate[yu]) I4[
     2] - 2/3 (Abs[yl]^4 + 
      Nc (Abs[yd]^4 + Abs[yu]^4 - 
         2 yd yu Conjugate[yd] Conjugate[yu])) I4[3], 
 rH4D23 -> (-yl Conjugate[aeH] - (aeH - 4 aHl3 yl) Conjugate[yl] - 
      NC (yd Conjugate[adH] + 
         yu Conjugate[auH] + (adH - 4 aHq3 yd + 2 aHud yu) Conjugate[
           yd] + (auH - 4 aHq3 yu + 2 yd Conjugate[aHud]) Conjugate[
           yu])) I4[2] - 
   4/3 (Abs[yl]^4 + 
      Nc (Abs[yd]^4 + Abs[yu]^4 + 
         yd yu Conjugate[yd] Conjugate[yu])) I4[3], 
 rH4D24 -> 
  1/2 I (-yl Conjugate[aeH] + aeH Conjugate[yl] + 
     NC (-yd Conjugate[adH] + yu Conjugate[auH] + adH Conjugate[yd] - 
        auH Conjugate[yu])) I4[2], 
 rH2BD2 -> 
  -1/18 g1 (3 NC yd Conjugate[yd] I4[3] + 3 yl Conjugate[yl] I4[3] + 
     3 NC yu Conjugate[yu] I4[3] - 
     4 I4[2] (3 (Tr[aHe] + Tr[aHl1]) + 
        NC (Tr[aHd] - Tr[aHq1] - 2 Tr[aHu]))), 
 rB2D2 -> -(1/45) g1^2 (27 + 11 NC) I4[3], 
 rB04D2 -> -(1/648) (-5 + d) (-1 + d) g1^4 (729 + 137 NC) I4[3], 
 aB02B2 -> -(1/648) (-5 + d) g1^4 (729 + 137 NC), 
 aH2B2 -> -(1/216)
     g1^2 (NC yd Conjugate[yd] + 81 yl Conjugate[yl] + 
     25 NC yu Conjugate[yu]) I4[3], 
 aH2B02D21 -> -(1/18)
     g1^2 ((-2 + d) NC yd Conjugate[yd] + 
     3 (-34 + 7 d) yl Conjugate[yl] + (-32 + 7 d) NC yu Conjugate[
       yu]) I4[3], 
 rH2B02D22 -> 
  1/108 g1^2 ((6 - 5 d) NC yd Conjugate[yd] I4[3] + 
     9 (14 - 5 d) yl Conjugate[yl] I4[
       3] + (42 - 17 d) NC yu Conjugate[yu] I4[3] + 
     12 (-1 + d) I4[
       2] (3 (Tr[aHe] + Tr[aHl1]) + 
        NC (Tr[aHd] - Tr[aHq1] - 2 Tr[aHu]))), 
 rH2B02D23 -> 
  1/54 g1^2 ((7 - 5 d) NC yd Conjugate[yd] I4[3] + 
     9 (23 - 5 d) yl Conjugate[yl] I4[
       3] + (67 - 17 d) NC yu Conjugate[yu] I4[3] + 
     6 (-1 + d) I4[
       2] (3 (Tr[aHe] + Tr[aHl1]) + 
        NC (Tr[aHd] - Tr[aHq1] - 2 Tr[aHu]))), 
 aH6 -> -2 (NC yd Conjugate[
        yd] (yd Conjugate[adH] + adH Conjugate[yd]) + 
      yl Conjugate[yl] (yl Conjugate[aeH] + aeH Conjugate[yl]) + 
      NC yu Conjugate[yu] (yu Conjugate[auH] + auH Conjugate[yu])) I4[
     2] - 2/3 (Abs[yl]^6 + NC (Abs[yd]^6 + Abs[yu]^6)) I4[3], 
 aB06 -> -(((-5 + d) (-3 + d) (-1 + d) g1^6 (24057 + 2081 NC) I4[3])/
   349920), 
 aH2B04 -> -(1/3888)
   g1^4 ((66 + d (-82 + 17 d)) Nc yd Conjugate[yd] I4[3] + 
      81 (186 + d (-122 + 17 d)) yl Conjugate[yl] I4[
        3] + (2046 + d (-1678 + 257 d)) Nc yu Conjugate[yu] I4[3] - 
      24 (-3 + d) (-1 + d) I4[
        2] (27 (4 Tr[aHe] + Tr[aHl1]) + 
         Nc (4 Tr[aHd] - Tr[aHq1] - 32 Tr[aHu]))), 
 aH4B02 -> 
  1/36 g1^2 ((6 - 5 d) NC yd Conjugate[adH] I4[2] + 
     9 (14 - 5 d) yl Conjugate[aeH] I4[
       2] + (42 - 17 d) NC yu Conjugate[auH] I4[2] + 
     NC (adH (6 - 5 d) + 36 aHd yd + 
        12 (-2 aHd + aHq1 + aHq3) d yd) Conjugate[yd] I4[2] + 
     126 aeH Conjugate[yl] I4[2] - 45 aeH d Conjugate[yl] I4[2] + 
     180 aHe yl Conjugate[yl] I4[2] + 
     144 aHl1 yl Conjugate[yl] I4[2] + 
     144 aHl3 yl Conjugate[yl] I4[2] - 
     72 aHe d yl Conjugate[yl] I4[2] - 
     36 aHl1 d yl Conjugate[yl] I4[2] - 
     36 aHl3 d yl Conjugate[yl] I4[2] + 
     42 auH NC Conjugate[yu] I4[2] - 
     17 auH d NC Conjugate[yu] I4[2] - 
     72 aHq1 NC yu Conjugate[yu] I4[2] + 
     72 aHq3 NC yu Conjugate[yu] I4[2] - 
     108 aHu NC yu Conjugate[yu] I4[2] + 
     12 aHq1 d NC yu Conjugate[yu] I4[2] - 
     12 aHq3 d NC yu Conjugate[yu] I4[2] + 
     48 aHu d NC yu Conjugate[yu] I4[2] + 7 NC Abs[yd]^4 I4[3] - 
     5 d NC Abs[yd]^4 I4[3] + 207 Abs[yl]^4 I4[3] - 
     45 d Abs[yl]^4 I4[3] + 67 NC Abs[yu]^4 I4[3] - 
     17 d NC Abs[yu]^4 I4[3]), KW -> g2^2 (1 + NC) I4[2], 
 KW0 -> 1/2 (-1 + d) g2^2 (1 + NC) I4[2], 
 mW02 -> -3 (-1 + d) g2^2 (1 + NC) I4[1], 
 lmbdH2W02 -> 
  1/4 (-2 + d) g2^2 (Abs[yl]^2 + NC yd Conjugate[yd] + 
     NC yu Conjugate[yu]) I4[2], 
 lmbdW04 -> 1/16 (-3 + d) (-1 + d) g2^4 (1 + NC) I4[2], 
 rW2D2 -> -(1/5) g2^2 (1 + NC) I4[3], 
 aW3 -> -(1/60) g2^3 (1 + NC) I4[3], 
 rW02D4 -> -(1/20) (-1 + d) g2^2 (1 + NC) I4[3], 
 rH2D4 -> -(1/
    3) (Abs[yl]^2 + NC yd Conjugate[yd] + NC yu Conjugate[yu]) I4[3], 
 rH2WD2 -> -(1/6)
     g2 (NC yd Conjugate[yd] I4[3] + yl Conjugate[yl] I4[3] + 
     NC yu Conjugate[yu] I4[3] + 4 I4[2] (Tr[aHl3] + NC Tr[aHq3])), 
 aH2W2 -> -(1/24)
     g2^2 (Abs[yl]^2 + NC yd Conjugate[yd] + NC yu Conjugate[yu]) I4[
    3], aH2W02D21 -> -(1/
    12) (-4 + d) g2^2 (Abs[yl]^2 + NC yd Conjugate[yd] + 
     NC yu Conjugate[yu]) I4[3], 
 aH2W02D22 -> -(1/6)
     g2^2 (NC yd Conjugate[yd] I4[3] + yl Conjugate[yl] I4[3] + 
     NC yu Conjugate[yu] I4[3] + 
     2 (-1 + d) I4[2] (Tr[aHl3] + NC Tr[aHq3])), 
 rH2W02D24 -> 
  1/12 (-4 + d) g2^2 (Abs[yl]^2 + NC yd Conjugate[yd] + 
     NC yu Conjugate[yu]) I4[3], 
 rH2W02D25 -> 
  1/6 g2^2 (NC yd Conjugate[yd] I4[3] + yl Conjugate[yl] I4[3] + 
     NC yu Conjugate[yu] I4[3] + 
     2 (-1 + d) I4[2] (Tr[aHl3] + NC Tr[aHq3])), 
 rW02D4 -> -(1/20) (-1 + d) g2^2 (1 + NC) I4[3], 
 rW02WD2 -> -(1/20) (4 + d) g2^3 (1 + NC) I4[3], 
 aW02W21 -> -(1/20) (-1 + d) g2^4 (1 + NC) I4[3], 
 aW02W22 -> -(1/40) (-23 + 3 d) g2^4 (1 + NC) I4[3], 
 aW04D21 -> -(1/40) (-3 + d) (-1 + d) g2^4 (1 + NC) I4[3], 
 rW04D22 -> -(1/20) (-1 + d) (-11 + 2 d) g2^4 (1 + NC) I4[3], 
 aH4W021 -> -(1/12)
     g2^2 (3 ((-2 + d) NC yd Conjugate[adH] + (-2 + d) yl Conjugate[
          aeH] + (-2 + d) NC yu Conjugate[auH] + 
        NC ((-2 + d) (adH - 4 aHq3 yd) + 2 aHud yu) Conjugate[yd] - 
        2 aeH Conjugate[yl] + aeH d Conjugate[yl] + 
        8 aHl3 yl Conjugate[yl] - 4 aHl3 d yl Conjugate[yl] - 
        2 auH NC Conjugate[yu] + auH d NC Conjugate[yu] + 
        8 aHq3 NC yu Conjugate[yu] - 4 aHq3 d NC yu Conjugate[yu] + 
        2 NC yd Conjugate[aHud] Conjugate[yu]) I4[2] + 
     2 (-3 + d) (Abs[yl]^4 + 
        NC (Abs[yd]^4 + Abs[yu]^4 + 
           yd yu Conjugate[yd] Conjugate[yu])) I4[3]), 
 aH4W022 -> (aHu - aHq1 (-2 + d)) g2^2 NC yu Conjugate[yu] I4[2] - 
   1/12 (-3 + d) g2^2 (Abs[yl]^4 + 
      NC (Abs[yd]^4 + Abs[yu]^4 - 
         2 yd yu Conjugate[yd] Conjugate[yu])) I4[3], 
 aH2W04 -> -(1/48)
     g2^4 ((6 + (-6 + d) d) NC yd Conjugate[yd] I4[
       3] + (6 + (-6 + d) d) yl Conjugate[yl] I4[
       3] + (6 + (-6 + d) d) NC yu Conjugate[yu] I4[3] + 
     8 (-3 + d) (-1 + d) I4[2] (Tr[aHl3] + NC Tr[aHq3])), 
 aW06 -> -(1/480) (-5 + d) (-3 + d) (-1 + d) g2^6 (1 + NC) I4[3], 
 lmbdB02W02 -> 1/24 (-3 + d) (-1 + d) g1^2 g2^2 (9 + NC) I4[2], 
 lmbdH2B0W0 -> 
  1/6 g1 g2 (d NC yd Conjugate[yd] I4[2] - 
     3 (-4 + d) yl Conjugate[yl] I4[2] - (-6 + d) NC yu Conjugate[
       yu] I4[2] - 
     4 (-1 + d) I4[
       1] (3 (Tr[aHe] + Tr[aHl1] - Tr[aHl3]) + 
        NC (Tr[aHd] - Tr[aHq1] - 3 Tr[aHq3] - 2 Tr[aHu]))), 
 rB02W02D22 -> -(1/72) (-5 + d) (-1 + d) g1^2 g2^2 (9 + NC) I4[3], 
 rB02W02D23 -> -(1/18) (-5 + d) (-1 + d) g1^2 g2^2 (9 + NC) I4[3], 
 aB02W02D21 -> -(1/72) (-5 + d) (-1 + d) g1^2 g2^2 (9 + NC) I4[3], 
 aB02W2 -> -(1/72) (-5 + d) g1^2 g2^2 (9 + NC) I4[3], 
 aB0BW0W -> -(1/18) (-5 + d) g1^2 g2^2 (9 + NC) I4[3], 
 aB2W02 -> -(1/72) (-5 + d) g1^2 g2^2 (9 + NC) I4[3], 
 aH2BW -> 1/
   36 g1 g2 (NC yd Conjugate[yd] + 9 yl Conjugate[yl] + 
     5 NC yu Conjugate[yu]) I4[3], 
 aH2B0W0D21 -> 
  1/9 g1 g2 ((-2 + d) NC yd Conjugate[yd] I4[3] + 
     3 (-4 + d) yl Conjugate[yl] I4[3] + (-7 + 2 d) NC yu Conjugate[
       yu] I4[3] - 
     2 (-1 + d) I4[
       2] (3 (Tr[aHe] + Tr[aHl1]) + 
        NC (Tr[aHd] - Tr[aHq1] - 2 Tr[aHu]))), 
 aH2B0W0D22 -> 
  1/18 g1 g2 (d NC yd Conjugate[yd] I4[3] - 
     3 (-4 + d) yl Conjugate[yl] I4[3] - (-6 + d) NC yu Conjugate[
       yu] I4[3] - 
     4 (-1 + d) I4[
       2] (3 (Tr[aHe] + Tr[aHl1]) + 
        NC (Tr[aHd] - Tr[aHq1] - 2 Tr[aHu]))), 
 rH2B0W0D24 -> 
  1/18 g1 g2 (d NC yd Conjugate[yd] I4[3] - 
     3 (-4 + d) yl Conjugate[yl] I4[3] - (-6 + d) NC yu Conjugate[
       yu] I4[3] - 
     2 (-1 + d) I4[
       2] (3 (Tr[aHe] + Tr[aHl1] - Tr[aHl3]) + 
        NC (Tr[aHd] - Tr[aHq1] - 3 Tr[aHq3] - 2 Tr[aHu]))), 
 rH2B0W0D25 -> 
  1/18 g1 g2 ((1 + d) NC yd Conjugate[yd] I4[3] - 
     3 (-7 + d) yl Conjugate[yl] I4[3] - (-11 + d) NC yu Conjugate[
       yu] I4[3] - 
     2 (-1 + d) I4[
       2] (3 (Tr[aHe] + Tr[aHl1]) + 
        NC (Tr[aHd] - Tr[aHq1] - 2 Tr[aHu]))), 
 aH2B02W02 -> -(1/72)
     g1^2 g2^2 ((-2 + d)^2 NC yd Conjugate[yd] I4[3] + 
     9 (-4 + d)^2 yl Conjugate[yl] I4[
       3] + (34 + (-10 + d) d) NC yu Conjugate[yu] I4[3] - 
     4 (-3 + d) (-1 + d) I4[
       2] (9 Tr[aHl1] - 9 Tr[aHl3] - NC (3 Tr[aHq1] + Tr[aHq3]))), 
 aH2B03W0 -> 
  1/324 g1^3 g2 (-((-3 + d + d^2) NC yd Conjugate[yd] I4[3]) + 
     27 (39 + (-13 + d) d) yl Conjugate[yl] I4[
       3] + (168 + (-44 + d) d) NC yu Conjugate[yu] I4[3] + 
     2 (-3 + d) (-1 + d) I4[
       2] (27 (4 Tr[aHe] + Tr[aHl1] - 3 Tr[aHl3]) + 
        NC (4 Tr[aHd] - Tr[aHq1] - 9 Tr[aHq3] - 32 Tr[aHu]))), 
 aH2B0W03 -> 
  1/36 g1 g2^3 (-((3 + (-5 + d) d) NC yd Conjugate[yd] I4[3]) + 
     3 (9 + (-7 + d) d) yl Conjugate[yl] I4[
       3] + (-6 + d) (-2 + d) NC yu Conjugate[yu] I4[3] + 
     6 (-3 + d) (-1 + d) I4[
       2] (3 Tr[aHl1] - Tr[aHl3] - NC (Tr[aHq1] + Tr[aHq3]))), 
 aH4B0W0 -> 
  1/6 g1 g2 ((-6 + d) NC yu Conjugate[auH] I4[2] + 
     NC (auH (-6 + d) + 
        4 (-3 aHu + (-aHq1 + aHq3 + 2 aHu) d) yu) Conjugate[yu] I4[
       2] - ((1 + d) NC Abs[yd]^4 - 
        3 (-7 + d) Abs[yl]^4 - (-11 + d) NC Abs[yu]^4) I4[3]), 
 aB04W02 -> -(((-5 + d) (-3 + d) (-1 + d) g1^4 g2^2 (81 + NC) I4[3])/
   2592), aB02W04 -> -(1/
    288) (-5 + d) (-3 + d) (-1 + d) g1^2 g2^4 (9 + NC) I4[3], 
 rB02D4 -> -(1/180) (-1 + d) g1^2 (27 + 11 NC) I4[3],
 aH2W0WD-> 1/6 g2^2 I4[3] NC yd^2-1/6 g2^2 I4[3] NC yu^2+1/6 g2^2 I4[3] yl^2,
 aH2W0BD->1/6 g1 g2 I4[3] NC yu^2-1/3 g1 g2 I4[3] yl^2,
 aH2B0WD->1/6 g1 g2 I4[3] NC yu^2-1/3 g1 g2 I4[3] yl^2,
 aH2B0BD->1/18 g1^2 I4[3] NC yd^2-7/18 g1^2 I4[3] NC yu^2+7/6 g1^2 I4[3] yl^2
};

(* Note: WCs are given in 4-dimensional units. Multiply by the corresponding powers of T to bring them to 3-dimensional ones *)

solMBos={mH2 -> ((d*(g1^2 + 3*g2^2) + 24*lam)*I4[1] - 
    4*aHD*(-(mu^2*I4[1]) + mu^4*I4[2]) + 
    8*aHDD*(-(mu^2*I4[1]) + mu^4*I4[2]) + (g1^2 + 3*(g2^2 - 8*lam))*mu^2*
     (I4[2] - mu^2*I4[3]))/4, KH -> aHD*(I4[1] - mu^2*I4[2]) - 
   2*aHDD*(I4[1] - mu^2*I4[2]) - ((g1^2 + 3*g2^2)*(3*I4[2] - 2*mu^2*I4[3]))/
    6, lmbdH4 -> (-(d*(g1^4 + 2*g1^2*g2^2 + 3*g2^4)*I4[2]) - 
    192*aH*(I4[1] - mu^2*I4[2]) + 4*aHD*(d*(g1^2 + g2^2)*I4[1] - 
      16*lam*(I4[1] - 2*mu^2*I4[2])) - 
    8*lam*(-24*aHDD*(I4[1] - 2*mu^2*I4[2]) - (g1^2 + 3*g2^2 - 24*lam)*
       (I4[2] - 2*mu^2*I4[3])))/16, 
 lmbdH2B02 -> -1/48*(g1^2*(-3*g1^2*I4[2] + 3*d*g1^2*I4[2] - 9*g2^2*I4[2] + 
     9*d*g2^2*I4[2] - 216*lam*I4[2] + 72*d*lam*I4[2] - 
     24*aHDD*((-3 + 2*d)*I4[1] - (-10 + 3*d)*mu^2*I4[2]) - 
     12*aHD*((-3 + 4*d)*I4[1] - (-8 + 3*d)*mu^2*I4[2]) + 6*g1^2*mu^2*I4[3] - 
     2*d*g1^2*mu^2*I4[3] + 18*g2^2*mu^2*I4[3] - 6*d*g2^2*mu^2*I4[3] + 
     720*lam*mu^2*I4[3] - 144*d*lam*mu^2*I4[3])), 
 mB02 -> g1^2*((-1 + d)*I4[1] - (-3 + d)*mu^2*I4[2] + mu^4*I4[3]), 
 KB0 -> -1/6*(g1^2*((-4 + d)*I4[2] - 2*(-6 + d)*mu^2*I4[3] )), 
 KB -> (g1^2*(I4[2] - 2*mu^2*I4[3]))/6, 
 lmbdB04 -> -1/48*((-3 + d)*g1^4*((-1 + d)*I4[2] - 2*(-5 + d)*mu^2*I4[3])), 
 KW -> (g2^2*((-49 + 2*d)*I4[2] - 2*mu^2*I4[3] ))/6, 
 KW0 -> -1/6*(g2^2*((52 + d*(-9 + 2*d))*I4[2] - 2*(-6 + d)*mu^2*I4[3])), 
 mW02 -> g2^2*((-1 + d)*(-1 + 2*d)*I4[1] - (-3 + d)*mu^2*I4[2] + 
    (-5 + d)*mu^4*I4[3]), rH2D4 -> ((g1^2 + 3*g2^2)*I4[3])/6, 
 rH2BD2 -> (2*(aHD + aHDD)*g1*I4[2] - g1*(g1^2 + 3*g2^2)*I4[3])/12, 
 aH2B2 -> (g1^2*(7*g1^2 + 21*g2^2 - 24*lam)*I4[3])/96, 
 aH2B02D21 -> (g1^2*(-3*(aHD - 2*aHDD)*(-3 + d)*I4[2] + 
     (-4 + d)*(g1^2 + 3*g2^2)*I4[3]))/12, 
 rH2B02D22 -> (g1^2*(-8*(aHD + aHDD)*(-4 + d)*I4[2] + 
     ((-2 + d)*(g1^2 + 3*g2^2) + 24*(-6 + d)*lam)*I4[3]))/48, 
 rH2B02D23 -> -1/24*(g1^2*(4*aHDD*(32 - 11*d)*I4[2] + 
     2*aHD*(-17 + 5*d)*I4[2] + (-5 + d)*(g1^2 + 3*g2^2 - 24*lam)*I4[3])), 
 aH4D21 -> -1/4*(aHDD*(g1^2 - 3*g2^2 + 56*lam)*I4[2]) + 
   aHD*((-3*g2^2*I4[2])/4 + lam*I4[2]) - 
   (((5 + 2*d)*g1^4 + (37 + 6*d)*g2^4 + 2*g1^2*((5 + 2*d)*g2^2 - 64*lam) - 
      224*g2^2*lam + 320*lam^2)*I4[3])/96, 
 aH4D22 -> (-3*aHD*(g1^2 + 21*g2^2 + 24*lam)*I4[2] + 
    g1^2*(-36*aHDD*I4[2] + ((15 - 2*d)*g2^2 + 40*lam)*I4[3]))/12, 
 rH4D23 -> aHDD*(-3*g2^2 + 4*lam)*I4[2] - 
   (aHD*(3*g1^2 - 3*g2^2 + 4*lam)*I4[2])/2 + 
   ((5*g1^4 + g1^2*((-5 + 2*d)*g2^2 + 8*lam) + 
      4*(g2^4 + 16*g2^2*lam + 8*lam^2))*I4[3])/12, 
 aH6 -> (-6*(aHD*d*(g1^2 + g2^2)^2 + 6*aH*(g1^2 + 3*g2^2 - 72*lam) + 
      32*(-5*aHD + 18*aHDD)*lam^2)*I4[2] + 
    (d*(g1^6 + 3*g1^4*g2^2 + 3*g1^2*g2^4 + 3*g2^6) - 
      48*(g1^2 + 3*g2^2 - 40*lam)*lam^2)*I4[3])/48, 
 rB2D2 -> -1/60*(g1^2*I4[3]), rB04D2 -> ((-5 + d)*(-4 + d)*g1^4*I4[3])/24, 
 aB02B2 -> -1/48*((-5 + d)*g1^4*I4[3]), 
 aB06 -> ((-5 + d)*(-3 + d)*(-1 + d)*g1^6*I4[3])/1440, 
 aH2B04 -> (g1^4*(-4*(-3 + d)*(-aHD - 10*aHDD + 4*(aHD + aHDD)*d)*I4[2] + 
     ((-3 + (-4 + d)*d)*(g1^2 + 3*g2^2) + 24*(-5 + d)*(-3 + d)*lam)*I4[3]))/
   192, aH4B02 -> 
  (g1^2*(-3*(-48*aH*(-3 + d) + aHD*((-2 + 3*d)*g1^2 + (-6 + 7*d)*g2^2 + 
         32*(-2 + d)*lam) + 4*aHDD*(-88*lam + d*(g1^2 + g2^2 + 28*lam)))*
      I4[2] + (d*g1^4 + 2*d*g1^2*g2^2 - 6*g2^4 + 2*d*g2^4 + 
       4*(3*(-1 + d)*g1^2 + (-9 + 5*d)*g2^2)*lam + 32*(-39 + 8*d)*lam^2)*
      I4[3]))/48, rW2D2 -> ((81 - 2*d)*g2^2*I4[3])/60, 
 aW3 -> ((-1 + 2*d)*g2^3*I4[3])/180, 
 lmbdW04 -> -1/48*((-3 + d)*g2^4*((7 - 15*d + 8*d^2)*I4[2] - 
     2*(-5 + d)*mu^2*I4[3])), lmbdH2W02 -> 
  -1/48*(g2^2*(-3*g1^2*I4[2] + 3*d*g1^2*I4[2] + 63*g2^2*I4[2] - 
     87*d*g2^2*I4[2] + 24*d^2*g2^2*I4[2] - 216*lam*I4[2] + 72*d*lam*I4[2] + 
     6*g1^2*mu^2*I4[3] - 2*d*g1^2*mu^2*I4[3] - 222*g2^2*mu^2*I4[3] + 
     42*d*g2^2*mu^2*I4[3] + 720*lam*mu^2*I4[3] - 144*d*lam*mu^2*I4[3] - 
     12*aHD*(I4[1] + (-4 + d)*mu^2*I4[2] + (11 - 2*d)*mu^4*I4[3]) - 
     24*aHDD*((-3 + 2*d)*I4[1] - (-10 + 3*d)*mu^2*I4[2] + 
       (-21 + 4*d)*mu^4*I4[3]) + 3*(g1^2 + (115 - 16*d)*g2^2 + 
       72*(-7 + d)*lam)*mu^4*I4[4])), 
 aW06 -> ((-5 + d)*(-3 + d)*(-1 + d)*(-31 + 32*d)*g2^6*I4[3])/1440, 
 rH2WD2 -> (g2*(2*aHDD*I4[2] - (g1^2 - 10*g2^2)*I4[3]))/12, 
 aH2W2 -> (g2^2*(7*g1^2 + (125 - 8*d)*g2^2 - 24*lam)*I4[3])/96, 
 aH4W021 -> 
  (g2^2*(-6*(-24*aH*(-3 + d) + aHD*d*g1^2 + d*(2*aHDD + aHD*(-7 + 2*d))*
        g2^2 + 4*(aHD*(10 - 3*d) + 2*aHDD*(-22 + 7*d))*lam)*I4[2] + 
     (4*d^2*g2^2*(g1^2 + 3*g2^2) + d*(g1^4 - 58*g2^4 - 76*g2^2*lam + 
         256*lam^2 + g1^2*(-19*g2^2 + 4*lam)) - 
       6*(g2^4 - 74*g2^2*lam + 2*lam*(g1^2 + 104*lam)))*I4[3]))/48, 
 aH4W022 -> (g2^2*(-3*aHD*((-2 + d)*g1^2 + (18 + (13 - 4*d)*d)*g2^2 + 
       24*(-2 + d)*lam)*I4[2] + d*g1^2*(-12*aHDD*I4[2] + 
       ((21 - 4*d)*g2^2 + 8*lam)*I4[3])))/48, 
 aH2W04 -> (g2^4*(-4*(-3 + d)*(3*aHD + 2*aHDD*(-5 + 2*d))*I4[2] + 
     ((-3 + (-4 + d)*d)*g1^2 + (-729 + d*(772 + d*(-285 + 32*d)))*g2^2 + 
       24*(-5 + d)*(-3 + d)*lam)*I4[3]))/192, 
 rW02D4 -> ((86 + d*(-13 + 2*d))*g2^2*I4[3])/60, 
 rW02WD2 -> ((166 + d*(-13 + 2*d))*g2^3*I4[3])/60, 
 aW02W21 -> -1/240*((1331 + d*(-383 + 12*d))*g2^4*I4[3]), 
 aW02W22 -> -1/60*((911 + 7*(-29 + d)*d)*g2^4*I4[3]), 
 aW04D21 -> ((-1658 + d*(655 + d*(-119 + 12*d)))*g2^4*I4[3])/120, 
 rW04D22 -> ((-1921 + 2*d*(355 + d*(-69 + 7*d)))*g2^4*I4[3])/60, 
 aH2W02D21 -> ((-4 + d)*g2^2*(7*g1^2 + (125 - 8*d)*g2^2 - 24*lam)*I4[3])/48, 
 rH2W02D24 -> -1/24*(g2^2*((-4 + d)*g1^2 + 3*(-24 + 5*d)*g2^2)*I4[3]), 
 rH2W02D25 -> (g2^2*(4*aHDD*(-4 + d)*I4[2] + 
     ((-13 + 3*d)*g1^2 + (-295 + (95 - 8*d)*d)*g2^2 - 24*(-5 + d)*lam)*
      I4[3]))/24, aH2W02D22 -> 
  (g2^2*(-2*aHDD*(-4 + d)*I4[2] - (g1^2 - 10*g2^2)*I4[3]))/12, 
 lmbdB02W02 -> -1/8*((-3 + d)*g1^2*g2^2*((-1 + d)*I4[2] - 
     2*(-5 + d)*mu^2*I4[3])), rH2W02D24 -> 
  -1/24*(g2^2*((-4 + d)*g1^2 + 3*(-24 + 5*d)*g2^2)*I4[3]), 
 aB02W2 -> -1/48*((-5 + d)*g1^2*g2^2*I4[3]), 
 rB02W02D22 -> ((-5 + d)*(-4 + d)*g1^2*g2^2*I4[3])/24, 
 rB02W02D23 -> ((-5 + d)*(-4 + d)*g1^2*g2^2*I4[3])/6, 
 aB02W02D21 -> ((-5 + d)*(-4 + d)*g1^2*g2^2*I4[3])/24, 
 aB0BW0W -> -1/12*((-5 + d)*g1^2*g2^2*I4[3]), 
 aB2W02 -> -1/48*((-5 + d)*g1^2*g2^2*I4[3]), 
 aH2BW -> (g1*g2*(7*g1^2 - 23*g2^2 - 8*lam)*I4[3])/48, 
 aH2B0W0D21 -> (g1*g2*(4*(aHD*(5 - 2*d) + aHDD*(-13 + 4*d))*I4[2] + 
     ((-13 + 3*d)*g1^2 + (61 - 11*d)*g2^2 - 8*(-5 + d)*lam)*I4[3]))/12, 
 aH2B0W0D22 -> (g1*g2*(8*(aHD + aHDD)*(-4 + d)*I4[2] + 
     (-((-2 + d)*g1^2) + (30 + d)*g2^2 - 8*(-6 + d)*lam)*I4[3]))/24, 
 rH2B0W0D24 -> (g1*g2*(4*(aHD + 2*aHDD)*(-4 + d)*I4[2] + 
     (-((-2 + d)*(g1^2 - 3*g2^2)) - 8*(-6 + d)*lam)*I4[3]))/24, 
 rH2B0W0D25 -> (g1*g2*(2*(aHD + 2*aHDD*(-4 + d) - aHD*d)*I4[2] + 
     (-5 + d)*(g1^2 - 5*g2^2 - 8*lam)*I4[3]))/24, 
 aH2B02W02 -> (g1^2*g2^2*(-4*(-3 + d)*(aHD + 2*aHD*d + 2*aHDD*(-5 + 2*d))*
      I4[2] + ((-3 + (-4 + d)*d)*g1^2 + (31 + d*(-20 + 3*d))*g2^2 + 
       24*(-5 + d)*(-3 + d)*lam)*I4[3]))/32, 
 aH2B03W0 -> (g1^3*g2*(-4*(-3 + d)*(aHD*(2 + d) + 2*aHDD*(-5 + 2*d))*I4[2] + 
     ((-3 + (-4 + d)*d)*g1^2 - (9 + (-8 + d)*d)*g2^2 + 
       8*(-5 + d)*(-3 + d)*lam)*I4[3]))/48, 
 aH2B0W03 -> (g1*g2^3*(-4*(-3 + d)*(3*aHD*d + 2*aHDD*(-5 + 2*d))*I4[2] + 
     ((-3 + (-4 + d)*d)*g1^2 - (129 + (-32 + d)*d)*g2^2 + 
       8*(-5 + d)*(-3 + d)*lam)*I4[3]))/48, 
 aH4B0W0 -> 
  (g1*g2*(-3*(-24*aH*(-3 + d) + aHD*((-2 + 3*d)*g1^2 + 3*(2 + d)*g2^2 + 
         16*(-1 + d)*lam) + 4*aHDD*(-88*lam + d*(g1^2 + g2^2 + 28*lam)))*
      I4[2] + (-6*(g2^4 + 6*g2^2*lam + 2*lam*(g1^2 + 64*lam)) + 
       d*(g1^4 + 2*g1^2*(g2^2 + 6*lam) + 2*(g2^4 + 10*g2^2*lam + 80*lam^2)))*
      I4[3]))/24, aB04W02 -> ((-5 + d)*(-3 + d)*(-1 + d)*g1^4*g2^2*I4[3])/96, 
 aB02W04 -> ((-5 + d)*(-3 + d)*(-1 + d)*g1^2*g2^4*I4[3])/96, 
 rW02D4 -> ((86 + d*(-13 + 2*d))*g2^2*I4[3])/60, 
 lmbdH2B0W0 -> -1/24*(g1*g2*(-3*g1^2*I4[2] + 3*d*g1^2*I4[2] + 27*g2^2*I4[2] - 
     3*d*g2^2*I4[2] - 72*lam*I4[2] + 24*d*lam*I4[2] + 6*g1^2*mu^2*I4[3] - 
     2*d*g1^2*mu^2*I4[3] + 18*g2^2*mu^2*I4[3] - 6*d*g2^2*mu^2*I4[3] + 
     240*lam*mu^2*I4[3] - 48*d*lam*mu^2*I4[3] - 
     12*aHD*((-1 + 2*d)*I4[1] - (-2 + d)*mu^2*I4[2] + mu^4*I4[3]) - 
     24*aHDD*((-3 + 2*d)*I4[1] - (-10 + 3*d)*mu^2*I4[2] + 
       (-21 + 4*d)*mu^4*I4[3]) + 3*(g1^2 + (-25 + 4*d)*g2^2 + 
       24*(-7 + d)*lam)*mu^4*I4[4])), rB02D4 -> ((-6 + d)*g1^2*I4[3])/60};


(* ::Section:: *)
(*CROSSCHECK*)


Expand[(-2(rH2B0W0D25)B/.solMBos)]//.{d->3,Abs[yt]->yt,g1->gY,g2->gw,lam->\[Lambda]1H,NC->3,Nc->3,yu->yt,yl->0,yd->0,Conjugate[yt]->yt,I4[2]->0,I4[1]->0}/.{F I4[3]->Zf[3,0],B I4[3]->Zb[3,0]}/.{B->0,F->0}//Factor


WilCoef=aW06;
Collect[(Expand[(WilCoef)B /.solMBos])+(Expand[(WilCoef)F/.solMFer])//.{d->3,Abs[yt]->yt,g1->gY,g2->gw,lam->\[Lambda]1H,NC->3,Nc->3,yu->yt,yl->0,yd->0,Conjugate[yt]->yt,I4[2]->0,I4[1]->0}/.{F I4[3]->Zf[3,0],B I4[3]->Zb[3,0]}/.{B->0,F->0},{_Zb,_Zf},Factor]


Collect[\[Alpha][BW0^6]/.sol6/.{\[Xi]->1}/.IPERCHARGE,{_Zb,_Zf},Factor]
