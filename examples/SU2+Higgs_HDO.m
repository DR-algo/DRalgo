(* ::Package:: *)

Quit[];


SetDirectory[NotebookDirectory[]];
DRalgo`DRalgo`$LoadGroupMath=True;
DRalgo`DRalgo`$GroupMathMultipleModels=True; (*Put this if you want to create multiple model-files with the same kernel*)
<<../Kernel/DRalgo.m


(* ::Chapter:: *)
(*SU2+Higgs*)


(*
	See 2503.20016 [hep-ph] for the dimension-6 operator matching and
	field-redefinition comparisons
*)


(* ::Section::Closed:: *)
(*Model*)


Group={"SU2"};
RepAdjoint={{2}};
HiggsDoublet={{{1}},"C"};
RepScalar={HiggsDoublet};
CouplingName={g2};


RepFermion={};


(* ::Text:: *)
(*The input for the gauge interactions toDRalgo are then given by*)


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
PerformDRhard[]


(* ::Section:: *)
(*Dimension 6 Matching*)


(*
	The scalar field \phi in the fundamental representation is
	internally decomposed into single real scalar fields as follows
	H=(\phi_1+I\phi_3,\phi_2+I\phi_4)^T,
	which can be compactly written as H_I=V_{I i}\phi_i with i=1,..,4. 
	All the operator tensors associated to higher-dimensional operators
	can be constructed in terms of
	the matrix V, and its complex conjugate Vs,
	the identity on the gauge fields id, which,
	since they are 3, it is a 3x3 identity matrix, 
	the generators in the fundamental representation, namely
	the Pauli matrices Sigma, and the group structure constants of SU(2) LC
*)
id=IdentityMatrix[3];
V={{1,0,I,0},{0,1,0,I}}/Sqrt[2];
Vs={{1,0,-I,0},{0,1,0,-I}}/Sqrt[2];
Sigma={{{0,1},{1,0}},{{0,-I},{I,0}},{{1,0},{0,-1}}};
LC=LeviCivitaTensor[3];


ge=g2*Sqrt[T];
Subscript[C, 1]=\[Alpha][W0^2 W^2,1]-(1/2)ge*\[Alpha][W0^2 W D^2];
Subscript[C, 2]=\[Alpha][W0^2 W^2,2]+(1/2)ge*\[Alpha][W0^2 W D^2];
Subscript[C, 3]=(\[Alpha][\[Phi]^4 D^2,2]-2\[Alpha][\[Phi]^4 D^2,1]);
Subscript[C, 4]=-(\[Alpha][\[Phi]^4 D^2,1]+I \[Alpha][\[Phi]^4 D^2,4]);
Subscript[Cs, 4]=-(\[Alpha][\[Phi]^4 D^2,1]-I \[Alpha][\[Phi]^4 D^2,4]);
Subscript[C, 5]=-2 \[Alpha][\[Phi]^2 W0^2 D^2,4]-\[Alpha][\[Phi]^2 W0^2 D^2,5]+I \[Alpha][\[Phi]^2 W0^2 D^2,3];
Subscript[Cs, 5]=-2 \[Alpha][\[Phi]^2 W0^2 D^2,4]-\[Alpha][\[Phi]^2 W0^2 D^2,5]-I \[Alpha][\[Phi]^2 W0^2 D^2,3];
Subscript[C, 6]=-\[Alpha][\[Phi]^2 W0^2 D^2,6]-I \[Alpha][\[Phi]^2 W0^2 D^2,2];
Subscript[Cs, 6]=-\[Alpha][\[Phi]^2 W0^2 D^2,6]+I \[Alpha][\[Phi]^2 W0^2 D^2,2];
Subscript[C, 7]=\[Alpha][\[Phi]^2 W0^2 D^2,1]-\[Alpha][\[Phi]^2 W0^2 D^2,5];


(*
	Here we construct the group tensors of the higher dimensional operators
*)
Tens1=Factorial[3]\[Alpha][W^3]LC;
Tens2=4 \[Alpha][\[Phi]^2 W^2]SymmetrizeTensor[Contract[V,Vs,id,{{1,3}}],{{1,2}}];
Tens3=4 Subscript[C, 1]TensorProduct[id,id]+4 Subscript[C, 2]SymmetrizeTensor[Transpose[TensorProduct[id,id],{1,3,4,2}],{{1,2},{3,4}}];
Tens4=4(Subscript[C, 3]SymmetrizeTensor[Contract[Vs,V,V,Vs,{{1,5},{3,7}}],{{1,2},{3,4}}]+Subscript[C, 4]SymmetrizeTensor[Contract[Vs,Vs,V,V,{{1,5},{3,7}}],{{1,2},{3,4}}]+Subscript[Cs, 4]SymmetrizeTensor[Contract[V,V,Vs,Vs,{{1,5},{3,7}}],{{1,2},{3,4}}]+\[Alpha][ \[Phi]^4 D^2,3]SymmetrizeTensor[Contract[Vs,V,Vs,V,{{1,3},{5,7}}],{{1,2},{3,4}}]);
Tens5=-8 \[Alpha][\[Phi]^2 W0^2 D^2,4]SymmetrizeTensor[Contract[Vs,V,id,{{1,3}}],{{1,2}}];
Tens6=(Subscript[C, 5])Contract[Vs,V,id,{{1,3}}]+(Subscript[Cs, 5])Contract[V,Vs,id,{{1,3}}]+(Subscript[C, 6])Contract[Sigma,V,Vs,LC,{{1,8},{2,6},{3,4}}]+(Subscript[Cs, 6])Conjugate[Contract[Sigma,V,Vs,LC,{{1,8},{2,6},{3,4}}]];
Tens7=4Subscript[C, 7]SymmetrizeTensor[Contract[Vs,V,id,{{1,3}}],{{1,2},{3,4}}];
Tens8=4 \[Alpha][W0^4 D^2,1]SymmetrizeTensor[TensorProduct[id,id],{{1,2},{3,4}}]+4\[Alpha][ W0^4 D^2,2]SymmetrizeTensor[Transpose[TensorProduct[id,id],{1,3,4,2}],{{1,2},{3,4}}];
Tens9=Factorial[6]\[Alpha][\[Phi]^6]SymmetrizeTensor[Contract[V,Vs,V,Vs,V,Vs,{{1,3},{5,7},{9,11}}],{{1,2,3,4,5,6}}];
Tens10=Factorial[2]Factorial[4](\[Alpha][\[Phi]^4 W0^2,1]SymmetrizeTensor[Contract[V,Vs,V,Vs,id,{{1,3},{5,7}}],{{1,2,3,4},{5,6}}]+\[Alpha][ \[Phi]^4 W0^2,2]SymmetrizeTensor[Contract[Vs,V,Vs,V,Sigma,Sigma,{{1,10},{3,11},{5,13},{7,14}}],{{1,2,3,4},{5,6}}]);
Tens11=Factorial[2]Factorial[4]\[Alpha][\[Phi]^2 W0^4]SymmetrizeTensor[Contract[V,Vs,id,id,{{1,3}}],{{1,2},{3,4,5,6}}];
Tens12=Factorial[6]\[Alpha][W0^6]SymmetrizeTensor[TensorProduct[id,id,id],{{1,2,3,4,5,6}}];
Tens13=2\[Alpha][W^2 D^2]id;
Tens14=I \[Alpha][\[Phi]^2 W D^2](Contract[V,Vs,Sigma,{{1,7},{3,6}}]-Contract[Vs,V,Sigma,{{1,6},{3,7}}]);
Tens15=2\[Alpha][\[Phi]^2 D^4]SymmetrizeTensor[Contract[V,Vs,{{1,3}}],{{1,2}}];
Tens16=-\[Alpha][W0^2 W D^2]LC;
Tens17=2\[Alpha][W0^2 D^4]id;


(*
	TensorList is the array of the various group tensors refered to
	higher dimensional operators
*)
TensorList={
	Tens1,Tens2,Tens3,Tens4,Tens5,Tens6,Tens7,Tens8,Tens9,Tens10,
	Tens11,Tens12,Tens13,Tens14,Tens15,Tens16,Tens17
	};
(*
	NList is the array that indicate to which operator
	the group tensors in TensorList refer to
*)
NList={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17};
(*
	WC is the array of the various Wilson Coefficients \[Alpha][...]
*)
WC = DeleteDuplicates[Cases[TensorList//Normal, \[Alpha][__], \[Infinity]]];


(*
	Dimension6Matching and Dimension5Matching find the values of
	the Wilson coefficients listed in WC,
	d is the number of spatial dimensions,
	Zb and Zf are 1-loop master integrals
*)
soltot=Dimension6Matching[TensorList,NList,WC,d][[1]];
soltot//Factor//TableForm


(*
	The matching can be done individually for each group tensor
*)
sol=Dimension6Matching[{Tens1},{1},{\[Alpha][W^3]},d][[1]]; 
sol//Factor//TableForm


(*
	The functions ODIM6 and ODIM5
	return the group tensors of the various operators
*)
ODIM6[1,d]//MatrixForm


(* ::Section:: *)
(*FIELD REDEFINITIONS (BOSONIC CONTR)*)


(*
	The current version of Higher Dimensional Operators routines
	does not include field redefinitions
*)


(*
	Through field redefinition is it possible
	to check gauge dependence cancellation
*)


(*
	With field redefinition,
	is it possible to compare with results already existing in literature,
	for this model with arXiv:2503.20016
*)


ge=g2*Sqrt[T];
(*
	Higher Dimensional operators found with DRalgo
*)
solDralgo=Join[{\[Lambda][1]->\[Lambda]1H*T,\[Lambda][2]->1/4 ge^2},soltot];
(*
	Higher Dimensional operators from arXiv:2503.20016
*)
solLit={
	aH->0,
	aHD->0,
	aHDD->0,
	g1->0,
	lam->\[Lambda]1H,
	\[Lambda][1]->\[Lambda]1H,
	\[Lambda][2]->1/4 g2^2,
	mH2->((d*(g1^2+3*g2^2)+24*lam)*I4[1]-4*aHD*(-(mu^2*I4[1])+mu^4*I4[2])+8*aHDD*(-(mu^2*I4[1])+mu^4*I4[2])+(g1^2+3*(g2^2-8*lam))*mu^2*(I4[2]-mu^2*I4[3]))/4,
	KH->aHD*(I4[1]-mu^2*I4[2])-2*aHDD*(I4[1]-mu^2*I4[2])-((g1^2+3*g2^2)*(3*I4[2]-2*mu^2*I4[3]))/6,
	lmbdH4->(-(d*(g1^4+2*g1^2*g2^2+3*g2^4)*I4[2])-192*aH*(I4[1]-mu^2*I4[2])+4*aHD*(d*(g1^2+g2^2)*I4[1]-16*lam*(I4[1]-2*mu^2*I4[2]))-8*lam*(-24*aHDD*(I4[1]-2*mu^2*I4[2])-(g1^2+3*g2^2-24*lam)*(I4[2]-2*mu^2*I4[3])))/16,
	KW->(g2^2*((-49+2*d)*I4[2]-2*mu^2*I4[3]))/6,
	KW0->-1/6*(g2^2*((52+d*(-9+2*d))*I4[2]-2*(-6+d)*mu^2*I4[3])),
	mW02->g2^2*((-1+d)*(-1+2*d)*I4[1]-(-3+d)*mu^2*I4[2]+(-5+d)*mu^4*I4[3]),
	\[Alpha][\[Phi]^2D^4]->((g1^2+3*g2^2)*I4[3])/6,
	\[Alpha][\[Phi]^4D^2,1]->-1/4*(aHDD*(g1^2-3*g2^2+56*lam)*I4[2])+aHD*((-3*g2^2*I4[2])/4+lam*I4[2])-(((5+2*d)*g1^4+(37+6*d)*g2^4+2*g1^2*((5+2*d)*g2^2-64*lam)-224*g2^2*lam+320*lam^2)*I4[3])/96,
	\[Alpha][\[Phi]^4D^2,2]->(-3*aHD*(g1^2+21*g2^2+24*lam)*I4[2]+g1^2*(-36*aHDD*I4[2]+((15-2*d)*g2^2+40*lam)*I4[3]))/12,
	\[Alpha][\[Phi]^4D^2,3]->aHDD*(-3*g2^2+4*lam)*I4[2]-(aHD*(3*g1^2-3*g2^2+4*lam)*I4[2])/2+((5*g1^4+g1^2*((-5+2*d)*g2^2+8*lam)+4*(g2^4+16*g2^2*lam+8*lam^2))*I4[3])/12,
	\[Alpha][\[Phi]^6]->(-6*(aHD*d*(g1^2+g2^2)^2+6*aH*(g1^2+3*g2^2-72*lam)+32*(-5*aHD+18*aHDD)*lam^2)*I4[2]+(d*(g1^6+3*g1^4*g2^2+3*g1^2*g2^4+3*g2^6)-48*(g1^2+3*g2^2-40*lam)*lam^2)*I4[3])/48,
	\[Alpha][W^2D^2]->((81-2*d)*g2^2*I4[3])/60,
	\[Alpha][W^3]->((-1+2*d)*g2^3*I4[3])/180,
	lmbdW04->-1/48*((-3+d)*g2^4*((7-15*d+8*d^2)*I4[2]-2*(-5+d)*mu^2*I4[3])),
	lmbdH2W02->-1/48*(g2^2*(-3*g1^2*I4[2]+3*d*g1^2*I4[2]+63*g2^2*I4[2]-87*d*g2^2*I4[2]+24*d^2*g2^2*I4[2]-216*lam*I4[2]+72*d*lam*I4[2]+6*g1^2*mu^2*I4[3]-2*d*g1^2*mu^2*I4[3]-222*g2^2*mu^2*I4[3]+42*d*g2^2*mu^2*I4[3]+720*lam*mu^2*I4[3]-144*d*lam*mu^2*I4[3]-12*aHD*(I4[1]+(-4+d)*mu^2*I4[2]+(11-2*d)*mu^4*I4[3])-24*aHDD*((-3+2*d)*I4[1]-(-10+3*d)*mu^2*I4[2]+(-21+4*d)*mu^4*I4[3])+3*(g1^2+(115-16*d)*g2^2+72*(-7+d)*lam)*mu^4*I4[4])),
	\[Alpha][W0^6]->((-5+d)*(-3+d)*(-1+d)*(-31+32*d)*g2^6*I4[3])/1440,
	\[Alpha][\[Phi]^2W D^2]->(g2*(2*aHDD*I4[2]-(g1^2-10*g2^2)*I4[3]))/12,
	\[Alpha][\[Phi]^2W^2]->(g2^2*(7*g1^2+(125-8*d)*g2^2-24*lam)*I4[3])/96,
	\[Alpha][\[Phi]^4W0^2,1]->(g2^2*(-6*(-24*aH*(-3+d)+aHD*d*g1^2+d*(2*aHDD+aHD*(-7+2*d))*g2^2+4*(aHD*(10-3*d)+2*aHDD*(-22+7*d))*lam)*I4[2]+(4*d^2*g2^2*(g1^2+3*g2^2)+d*(g1^4-58*g2^4-76*g2^2*lam+256*lam^2+g1^2*(-19*g2^2+4*lam))-6*(g2^4-74*g2^2*lam+2*lam*(g1^2+104*lam)))*I4[3]))/48,
	\[Alpha][\[Phi]^4W0^2,2]->(g2^2*(-3*aHD*((-2+d)*g1^2+(18+(13-4*d)*d)*g2^2+24*(-2+d)*lam)*I4[2]+d*g1^2*(-12*aHDD*I4[2]+((21-4*d)*g2^2+8*lam)*I4[3])))/48,
	\[Alpha][\[Phi]^2W0^4]->(g2^4*(-4*(-3+d)*(3*aHD+2*aHDD*(-5+2*d))*I4[2]+((-3+(-4+d)*d)*g1^2+(-729+d*(772+d*(-285+32*d)))*g2^2+24*(-5+d)*(-3+d)*lam)*I4[3]))/192,
	\[Alpha][W0^2D^4]->((86+d*(-13+2*d))*g2^2*I4[3])/60,
	\[Alpha][W0^2W D^2]->((166+d*(-13+2*d))*g2^3*I4[3])/60,
	\[Alpha][W0^2W^2,1]->-1/240*((1331+d*(-383+12*d))*g2^4*I4[3]),
	\[Alpha][W0^2W^2,2]->-1/60*((911+7*(-29+d)*d)*g2^4*I4[3]),
	\[Alpha][W0^4D^2,1]->((-1658+d*(655+d*(-119+12*d)))*g2^4*I4[3])/120,
	\[Alpha][W0^4D^2,2]->((-1921+2*d*(355+d*(-69+7*d)))*g2^4*I4[3])/60,
	\[Alpha][\[Phi]^2W0^2D^2,1]->((-4+d)*g2^2*(7*g1^2+(125-8*d)*g2^2-24*lam)*I4[3])/48,
	\[Alpha][\[Phi]^2W0^2D^2,4]->-1/24*(g2^2*((-4+d)*g1^2+3*(-24+5*d)*g2^2)*I4[3]),
	\[Alpha][\[Phi]^2W0^2D^2,5]->(g2^2*(4*aHDD*(-4+d)*I4[2]+((-13+3*d)*g1^2+(-295+(95-8*d)*d)*g2^2-24*(-5+d)*lam)*I4[3]))/24,
	\[Alpha][\[Phi]^2W0^2D^2,2]->(g2^2*(-2*aHDD*(-4+d)*I4[2]-(g1^2-10*g2^2)*I4[3]))/12,
	\[Alpha][\[Phi]^2W0^2D^2,4]->-1/24*(g2^2*((-4+d)*g1^2+3*(-24+5*d)*g2^2)*I4[3]),
	\[Alpha][W0^2D^4]->((86+d*(-13+2*d))*g2^2*I4[3])/60
};


\[Alpha][W^3]/.solDralgo/.{Global`g2->g2}//Factor
\[Alpha][W^3]/.solLit//Factor


\[Alpha][\[Phi]^6]+\[Lambda][1](4*\[Alpha][\[Phi]^4 D^2,1]+ge^2 \[Alpha][W^2D^2]-2*ge*\[Alpha][\[Phi]^2W D^2]+4*\[Lambda][1]\[Alpha][\[Phi]^2D^4])/.solDralgo/.{Global`g2->g2}//Factor
\[Alpha][\[Phi]^6]+\[Lambda][1](4\[Alpha][\[Phi]^4 D^2,1]+g2^2 \[Alpha][W^2D^2]-2*g2*\[Alpha][\[Phi]^2W D^2]+4\[Lambda][1]\[Alpha][\[Phi]^2D^4])//.solLit//Factor


\[Alpha][\[Phi]^4 W0^2,1]+4\[Lambda][1]\[Alpha][\[Phi]^2W0^2D^2,4]+2\[Alpha][\[Phi]^2 W0^2 D^2,1](\[Lambda][1]-\[Lambda][2])+2\[Lambda][2]\[Alpha][\[Phi]^4 D^2,1]+2\[Lambda][2]\[Alpha][\[Phi]^2W0^2D^2,5]-ge*\[Lambda][2]\[Alpha][\[Phi]^2 W D^2]+(1/2)\[Lambda][2]ge^2 \[Alpha][W^2 D^2]+4\[Lambda][1]\[Lambda][2]\[Alpha][\[Phi]^2D^4]+4\[Lambda][2]^2 \[Alpha][W0^2D^4]/.solDralgo/.{Global`g2->g2}//Factor
\[Alpha][\[Phi]^4 W0^2,1]+4\[Lambda][1]\[Alpha][\[Phi]^2W0^2D^2,4]+2\[Alpha][\[Phi]^2 W0^2 D^2,1](\[Lambda][1]-\[Lambda][2])+2\[Lambda][2]\[Alpha][\[Phi]^4 D^2,1]+2\[Lambda][2]\[Alpha][\[Phi]^2W0^2D^2,5]-g2*\[Lambda][2]\[Alpha][\[Phi]^2 W D^2]+(1/2)\[Lambda][2]g2^2 \[Alpha][W^2 D^2]+4\[Lambda][1]\[Lambda][2]\[Alpha][\[Phi]^2D^4]+4\[Lambda][2]^2 \[Alpha][W0^2D^4]//.solLit//Factor


\[Alpha][\[Phi]^4 W0^2,2]/.solDralgo//Factor
\[Alpha][\[Phi]^4 W0^2,2]//.solLit//Factor


\[Alpha][\[Phi]^2 W0^4]+\[Lambda][2](-2\[Alpha][W0^4 D^2,1]+\[Alpha][\[Phi]^2W0^2 D^2,1]+2\[Alpha][\[Phi]^2W0^2D^2,4]+2ge*\[Alpha][W0^2W D^2]-2ge^2 \[Alpha][W^2D^2]+\[Lambda][2]\[Alpha][\[Phi]^2D^4])/.solDralgo/.{Global`g2->g2}//Factor
\[Alpha][\[Phi]^2 W0^4]+\[Lambda][2](-2\[Alpha][W0^4 D^2,1]+\[Alpha][\[Phi]^2W0^2 D^2,1]+2\[Alpha][\[Phi]^2W0^2D^2,4]+2g2*\[Alpha][W0^2W D^2]-2g2^2 \[Alpha][W^2D^2]+\[Lambda][2]\[Alpha][\[Phi]^2D^4])//.solLit//Factor


\[Alpha][W0^6]/.solDralgo//Factor
\[Alpha][W0^6]//.solLit//Factor


\[Alpha][\[Phi]^4 D^2,2]/.solDralgo//Factor
\[Alpha][\[Phi]^4 D^2,2]//.solLit//Factor


2\[Alpha][\[Phi]^4 D^2,1]+\[Alpha][\[Phi]^4 D^2,3]-(3/2)ge (2\[Alpha][\[Phi]^2 W D^2]-ge \[Alpha][W^2 D^2])/.solDralgo/.{Global`g2->g2}//Factor
2\[Alpha][\[Phi]^4 D^2,1]+\[Alpha][\[Phi]^4 D^2,3]-(3/2)g2 (2\[Alpha][\[Phi]^2 W D^2]-g2 \[Alpha][W^2 D^2])//.solLit//Factor


\[Alpha][\[Phi]^2 W0^2 D^2,1]/.solDralgo/.{Global`g2->g2} //Factor
\[Alpha][\[Phi]^2 W0^2 D^2,1]//.solLit//Factor


2\[Alpha][\[Phi]^2 W0^2 D^2,2]-ge*(\[Alpha][W0^2 W D^2]+2(\[Alpha][\[Phi]^2 W D^2]-ge*\[Alpha][W^2D^2]))/.solDralgo/.{Global`g2->g2}//Factor
2\[Alpha][\[Phi]^2 W0^2 D^2,2]-g2*(\[Alpha][W0^2 W D^2]+2(\[Alpha][\[Phi]^2 W D^2]-g2*\[Alpha][W^2D^2]))//.solLit//Factor


\[Alpha][\[Phi]^2W^2]/.solDralgo/.{Global`g2->g2}//Factor
\[Alpha][\[Phi]^2W^2]//.solLit//Factor


\[Alpha][W0^2W^2,2]+(1/2)ge*\[Alpha][W0^2W D^2]/.solDralgo/.{Global`g2->g2}//Factor
\[Alpha][W0^2W^2,2]+(1/2)g2*\[Alpha][W0^2W D^2]//.solLit//Factor


\[Alpha][W0^2W^2,1]-(1/2)ge*\[Alpha][W0^2W D^2]/.solDralgo/.{Global`g2->g2}//Factor
\[Alpha][W0^2W^2,1]-(1/2)g2*\[Alpha][W0^2W D^2]//.solLit//Factor



