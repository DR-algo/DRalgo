(* ::Package:: *)

Quit[];


If[$InputFileName=="",
	SetDirectory[NotebookDirectory[]],
	SetDirectory[DirectoryName[$InputFileName]]
];

DRalgo`DRalgo`$LoadGroupMath=True;
DRalgo`DRalgo`$GroupMathMultipleModels=True;
<<../Kernel/DRalgo.m


(* ::Chapter:: *)
(*SU2+Higgs higher-dimensional operators*)


(* ::Section::Closed:: *)
(*Model*)


Group={"SU2"};
RepAdjoint={{2}};
HiggsDoublet={{{1}},"C"};
RepScalar={HiggsDoublet};
CouplingName={gw};
RepFermion={};


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];

InputInv={{1,1},{True,False}};
MassTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;
VMass=m2*MassTerm1;
\[Mu]ij=GradMass[VMass]//Simplify//SparseArray;

QuarticTerm1=MassTerm1^2;
VQuartic=\[Lambda]1H*QuarticTerm1;
\[Lambda]4=GradQuartic[VQuartic];


(* ::Section::Closed:: *)
(*Dimensional Reduction*)


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PerformDRhard[];


(* ::Section:: *)
(*Dimension 6 Matching*)


id=IdentityMatrix[3];
V={{1,0,I,0},{0,1,0,I}}/Sqrt[2];
Vs={{1,0,-I,0},{0,1,0,-I}}/Sqrt[2];
Sigma={{{0,1},{1,0}},{{0,-I},{I,0}},{{1,0},{0,-1}}};
LC=LeviCivitaTensor[3];

ge=gw*Sqrt[T];
Subscript[C, 1]=\[Alpha][W0^2 W^2,1]-(1/2) ge*\[Alpha][W0^2 W D^2];
Subscript[C, 2]=\[Alpha][W0^2 W^2,2]+(1/2) ge*\[Alpha][W0^2 W D^2];
Subscript[C, 3]=\[Alpha][\[Phi]^4 D^2,2]-2 \[Alpha][\[Phi]^4 D^2,1];
Subscript[C, 4]=-(\[Alpha][\[Phi]^4 D^2,1]+I \[Alpha][\[Phi]^4 D^2,4]);
Subscript[Cs, 4]=-(\[Alpha][\[Phi]^4 D^2,1]-I \[Alpha][\[Phi]^4 D^2,4]);
Subscript[C, 5]=-2 \[Alpha][\[Phi]^2 W0^2 D^2,4]-\[Alpha][\[Phi]^2 W0^2 D^2,5]+I \[Alpha][\[Phi]^2 W0^2 D^2,3];
Subscript[Cs, 5]=-2 \[Alpha][\[Phi]^2 W0^2 D^2,4]-\[Alpha][\[Phi]^2 W0^2 D^2,5]-I \[Alpha][\[Phi]^2 W0^2 D^2,3];
Subscript[C, 6]=-\[Alpha][\[Phi]^2 W0^2 D^2,6]-I \[Alpha][\[Phi]^2 W0^2 D^2,2];
Subscript[Cs, 6]=-\[Alpha][\[Phi]^2 W0^2 D^2,6]+I \[Alpha][\[Phi]^2 W0^2 D^2,2];
Subscript[C, 7]=\[Alpha][\[Phi]^2 W0^2 D^2,1]-\[Alpha][\[Phi]^2 W0^2 D^2,5];

Tens1=Factorial[3] \[Alpha][W^3] LC;
Tens2=4 \[Alpha][\[Phi]^2 W^2] SymmetrizeTensor[Contract[V,Vs,id,{{1,3}}],{{1,2}}];
Tens3=4 Subscript[C, 1] TensorProduct[id,id]+4 Subscript[C, 2] SymmetrizeTensor[Transpose[TensorProduct[id,id],{1,3,4,2}],{{1,2},{3,4}}];
Tens4=4 (
	Subscript[C, 3] SymmetrizeTensor[Contract[Vs,V,V,Vs,{{1,5},{3,7}}],{{1,2},{3,4}}]
	+Subscript[C, 4] SymmetrizeTensor[Contract[Vs,Vs,V,V,{{1,5},{3,7}}],{{1,2},{3,4}}]
	+Subscript[Cs, 4] SymmetrizeTensor[Contract[V,V,Vs,Vs,{{1,5},{3,7}}],{{1,2},{3,4}}]
	+\[Alpha][\[Phi]^4 D^2,3] SymmetrizeTensor[Contract[Vs,V,Vs,V,{{1,3},{5,7}}],{{1,2},{3,4}}]
);
Tens5=-8 \[Alpha][\[Phi]^2 W0^2 D^2,4] SymmetrizeTensor[Contract[Vs,V,id,{{1,3}}],{{1,2}}];
Tens6=(Subscript[C, 5]) Contract[Vs,V,id,{{1,3}}]+(Subscript[Cs, 5]) Contract[V,Vs,id,{{1,3}}]+(Subscript[C, 6]) Contract[Sigma,V,Vs,LC,{{1,8},{2,6},{3,4}}]+(Subscript[Cs, 6]) Conjugate[Contract[Sigma,V,Vs,LC,{{1,8},{2,6},{3,4}}]];
Tens7=4 Subscript[C, 7] SymmetrizeTensor[Contract[Vs,V,id,{{1,3}}],{{1,2},{3,4}}];
Tens8=4 \[Alpha][W0^4 D^2,1] SymmetrizeTensor[TensorProduct[id,id],{{1,2},{3,4}}]+4 \[Alpha][W0^4 D^2,2] SymmetrizeTensor[Transpose[TensorProduct[id,id],{1,3,4,2}],{{1,2},{3,4}}];
Tens9=Factorial[6] \[Alpha][\[Phi]^6] SymmetrizeTensor[Contract[V,Vs,V,Vs,V,Vs,{{1,3},{5,7},{9,11}}],{{1,2,3,4,5,6}}];
Tens10=Factorial[2] Factorial[4] (
	\[Alpha][\[Phi]^4 W0^2,1] SymmetrizeTensor[Contract[V,Vs,V,Vs,id,{{1,3},{5,7}}],{{1,2,3,4},{5,6}}]
	+\[Alpha][\[Phi]^4 W0^2,2] SymmetrizeTensor[Contract[Vs,V,Vs,V,Sigma,Sigma,{{1,10},{3,11},{5,13},{7,14}}],{{1,2,3,4},{5,6}}]
);
Tens11=Factorial[2] Factorial[4] \[Alpha][\[Phi]^2 W0^4] SymmetrizeTensor[Contract[V,Vs,id,id,{{1,3}}],{{1,2},{3,4,5,6}}];
Tens12=Factorial[6] \[Alpha][W0^6] SymmetrizeTensor[TensorProduct[id,id,id],{{1,2,3,4,5,6}}];
Tens13=2 \[Alpha][W^2 D^2] id;
Tens14=I \[Alpha][\[Phi]^2 W D^2] (Contract[V,Vs,Sigma,{{1,7},{3,6}}]-Contract[Vs,V,Sigma,{{1,6},{3,7}}]);
Tens15=2 \[Alpha][\[Phi]^2 D^4] SymmetrizeTensor[Contract[V,Vs,{{1,3}}],{{1,2}}];
Tens16=-\[Alpha][W0^2 W D^2] LC;
Tens17=2 \[Alpha][W0^2 D^4] id;

TensorList={
	Tens1,Tens2,Tens3,Tens4,Tens5,Tens6,Tens7,Tens8,Tens9,Tens10,
	Tens11,Tens12,Tens13,Tens14,Tens15,Tens16,Tens17
};
NList={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17};
WC=DeleteDuplicates[Cases[TensorList//Normal, \[Alpha][__], \[Infinity]]];

soltot=Dimension6Matching[TensorList,NList,WC,3][[1]];

(* ::Section:: *)
(*Physical Dimension 6 Matching*)

Xb=Contract[V,Vs,{{1,3}}];
Qb=I (Contract[V,Vs,{{1,3}}]-Contract[Vs,V,{{1,3}}]);
Zeps=Contract[Sigma,V,Vs,LC,{{1,8},{2,6},{3,4}}];

OTP[1]=Factorial[3] \[Beta][W^3] LC;
OTP[2]=4 \[Beta][\[Phi]^2 W^2] SymmetrizeTensor[Contract[V,Vs,id,{{1,3}}],{{1,2}}];
OTP[3]=4 \[Beta][W0^2 W^2,1] TensorProduct[id,id]+4 \[Beta][W0^2 W^2,2] SymmetrizeTensor[Transpose[TensorProduct[id,id],{1,3,4,2}],{{1,2},{3,4}}];
OTP[4]=SymmetrizeTensor[
    \[Beta][D^2 \[Phi]^4,1] SymmetrizeTensor[Transpose[TensorProduct[Qb,Qb],{1,3,2,4}],{{1,2},{3,4}}]
  + \[Beta][D^2 \[Phi]^4,2] (SymmetrizeTensor[TensorProduct[Xb,Xb],{{1,2},{3,4}}]
      - SymmetrizeTensor[Transpose[TensorProduct[Xb,Xb],{1,3,4,2}],{{1,2},{3,4}}]),"P"];
OTP[5]=SymmetrizeTensor[4 \[Beta][D^2 \[Phi]^2 W0^2,1] SymmetrizeTensor[Contract[V,Vs,id,{{1,3}}],{{1,2}}],"SS"];
OTP[6]=SymmetrizeTensor[-I \[Beta][D^2 \[Phi]^2 W0^2,2] (Zeps-Conjugate[Zeps]),"AA"];
OTP[7]=SymmetrizeTensor[\[Beta][D^2 W0^4] (SymmetrizeTensor[Transpose[TensorProduct[id,id],{1,3,4,2}],{{1,2},{3,4}}]-TensorProduct[id,id]),"P"];
OTP[8]=Factorial[6] \[Beta][\[Phi]^6] SymmetrizeTensor[Contract[V,Vs,V,Vs,V,Vs,{{1,3},{5,7},{9,11}}],{{1,2,3,4,5,6}}];
OTP[9]=Factorial[2] Factorial[4] (
    \[Beta][\[Phi]^4 W0^2,1] SymmetrizeTensor[Contract[V,Vs,V,Vs,id,{{1,3},{5,7}}],{{1,2,3,4},{5,6}}]
  + \[Beta][\[Phi]^4 W0^2,2] SymmetrizeTensor[Contract[Vs,V,Vs,V,Sigma,Sigma,{{1,10},{3,11},{5,13},{7,14}}],{{1,2,3,4},{5,6}}]);
OTP[10]=Factorial[2] Factorial[4] \[Beta][\[Phi]^2 W0^4] SymmetrizeTensor[Contract[V,Vs,id,id,{{1,3}}],{{1,2},{3,4,5,6}}];
OTP[11]=Factorial[6] \[Beta][W0^6] SymmetrizeTensor[TensorProduct[id,id,id],{{1,2,3,4,5,6}}];

TensorListP={OTP[1],OTP[2],OTP[3],OTP[4],OTP[5],OTP[6],OTP[7],OTP[8],OTP[9],OTP[10],OTP[11]};
NListP={1,2,3,4,5,6,7,8,9,10,11};
WCsP=DeleteDuplicates[Cases[TensorListP//Normal,\[Beta][__],\[Infinity]]];

solP=Dimension6Matching["PHYS",TensorListP,NListP,WCsP,3][[1]];

solP

(* ::Section:: *)
(*TESTS*)


testList={};

AppendTo[testList,
TestCreate[
	\[Alpha][W^3]/.soltot,
	gw^3 Sqrt[T] Zb[3,0]/36
]];

AppendTo[testList,
TestCreate[
	\[Alpha][\[Phi]^2 W^2]/.soltot,
	(101 gw^4 T Zb[3,0]-24 gw^2 T \[Lambda]1H Zb[3,0])/96
]];

AppendTo[testList,
TestCreate[
	\[Alpha][W^2 D^2]/.soltot,
	(22 gw^2 Zb[3,0]-6 gw^2 xi Zb[3,0]-gw^2 xi^2 Zb[3,0])/12
]];

AppendTo[testList,
TestCreate[
	\[Alpha][\[Phi]^6]/.soltot,
	(3 gw^6 T^2 Zb[3,0]-48 gw^2 T^2 xi \[Lambda]1H^2 Zb[3,0]+640 T^2 \[Lambda]1H^3 Zb[3,0])/16
]];

AppendTo[testList,
TestCreate[
	\[Alpha][W0^6]/.soltot,
	0
]];

AppendTo[testList,
TestCreate[
	ODIM6[1,3][[1,2,3]],
	gw^3 Sqrt[T] Zb[3,0]/6
]];

AppendTo[testList,
TestCreate[
	ODIM6[13,3][[1,1]],
	-1/6 (gw^2 (-22+6 xi+xi^2) Zb[3,0])
]];

expectedP={
 \[Beta][W^3]->gw^3 Sqrt[T] Zb[3,0]/36,
 \[Beta][\[Phi]^2 W^2]->(101 gw^4 T Zb[3,0]-24 gw^2 T \[Lambda]1H Zb[3,0])/96,
 \[Beta][W0^2 W^2,1]->-(29/12) gw^4 T Zb[3,0],
 \[Beta][W0^2 W^2,2]->-(39/8) gw^4 T Zb[3,0],
 \[Beta][D^2 \[Phi]^4,1]->-(1/96) T (23 gw^4-160 gw^2 \[Lambda]1H+64 \[Lambda]1H^2) Zb[3,0],
 \[Beta][D^2 \[Phi]^4,2]->-(1/24) T (23 gw^4-160 gw^2 \[Lambda]1H+64 \[Lambda]1H^2) Zb[3,0],
 \[Beta][D^2 \[Phi]^2 W0^2,1]->(1/96) gw^2 T (101 gw^2-24 \[Lambda]1H) Zb[3,0],
 \[Beta][D^2 \[Phi]^2 W0^2,2]->(1/24) gw^4 T Zb[3,0],
 \[Beta][D^2 W0^4]->-(1/36) gw^4 T Zb[3,0],
 \[Beta][\[Phi]^6]->(1/48) T^2 (9 gw^6-84 gw^4 \[Lambda]1H+80 gw^2 \[Lambda]1H^2+1408 \[Lambda]1H^3) Zb[3,0],
 \[Beta][W0^2 \[Phi]^4,1]->(1/48) gw^2 T^2 (-101 gw^4+508 gw^2 \[Lambda]1H-544 \[Lambda]1H^2) Zb[3,0],
 \[Beta][W0^2 \[Phi]^4,2]->0,
 \[Beta][W0^4 \[Phi]^2]->(175/72) gw^6 T^2 Zb[3,0],
 \[Beta][W0^6]->0};

AppendTo[testList,
TestCreate[Union[Simplify[(First/@expectedP/.solP)-(Last/@expectedP)]],
	{0}
]];

AppendTo[testList,
TestCreate[
	ODIM6["PHYS",1,3][[1,2,3]],
	gw^3 Sqrt[T] Zb[3,0]/6
]];

AppendTo[testList,
TestCreate[
	Table[
		With[{T4=Normal[ODIM6["PHYS",n,3]]},
			Union[Flatten[Simplify[{
				T4-Transpose[T4,{2,1,3,4}],
				T4-Transpose[T4,{1,2,4,3}],
				T4-Transpose[T4,{3,4,1,2}],
				T4+Transpose[T4,{1,3,2,4}]+Transpose[T4,{1,4,3,2}]}]]]],
		{n,{4,7}}],
	{{0},{0}}
]];

AppendTo[testList,
TestCreate[
	With[{T5=Normal[ODIM6["PHYS",5,3]],T6=Normal[ODIM6["PHYS",6,3]]},
		Union[Flatten[Simplify[{
			T5-Transpose[T5,{2,1,3,4}],T5-Transpose[T5,{1,2,4,3}],
			T6+Transpose[T6,{2,1,3,4}],T6+Transpose[T6,{1,2,4,3}]}]]]],
	{0}
]];

AppendTo[testList,
TestCreate[
	Union[Flatten[MapThread[
		Simplify[(Normal[#1]/.solP)-Normal[ODIM6["PHYS",#2,3]]]&,
		{TensorListP,NListP}]]],
	{0}
]];

AppendTo[testList,
TestCreate[
	Union[Simplify[D[solP[[All,2]],xi]]],
	{0}
]];

report=TestReport[testList];
Print[report["ResultsDataset"]]
