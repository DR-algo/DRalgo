(* ::Package:: *)

(*Quit[];*)


SetDirectory[NotebookDirectory[]];
$LoadGroupMath=True;
$GroupMathMultipleModels=True; (*Put this if you want to create multiple model-files with the same kernel*)
<<../DRalgo.m


(* ::Chapter:: *)
(*SU5*)


(*see 9702255 [hep-ph]*)


(* ::Section:: *)
(*Model*)


Group={"SU5"};
RepAdjoint={{1,0,0,1}};
HiggsAdjoint={{{1,0,0,1}},"R"};
HiggsFundamental={{{1,0,0,0}},"C"};
RepScalar={HiggsAdjoint,HiggsFundamental};
CouplingName={g1};


Rep1={{{0,1,0,0}},"L"};
Rep2={{{0,0,0,1}},"L"};


RepFermion1Gen={Rep1,Rep2};
RepFermion3Gen={RepFermion1Gen,RepFermion1Gen,RepFermion1Gen}//Flatten[#,1]&;


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion3Gen,RepScalar];


InputInv={{1,1},{True,True}}; (*Tr \[CapitalPhi]^2*)
MassTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{2,2},{True,False}}; (*H H^+*)
MassTerm2=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;


VMass=(
	+m\[CapitalPhi]/2*MassTerm1
	+mH*MassTerm2
	);


\[Mu]ij=GradMass[VMass]//SparseArray;


InputInv={{1,1,1,1},{True,True,True,True}}; (*For this model there are two possible quartics for the adjoint*)
QuarticTemp=CreateInvariant[Group,RepScalar,InputInv];


InputInv={{2,1,1,2},{True,True,True,False}}; (*There are also two mixed terms*)
QuarticMix=CreateInvariant[Group,RepScalar,InputInv];


VQuartic=(
	+\[Lambda]H*QuarticTemp[[1]]
	+\[Lambda]S*QuarticTemp[[2]]
	+\[Lambda]*MassTerm2^2
	+\[Lambda]M1*QuarticMix[[1]]
	+\[Lambda]M2*QuarticMix[[2]]
	);


\[Lambda]4=GradQuartic[VQuartic]//SparseArray;


InputInv={{2,1,2},{False,True,True}};  (*H^+ Subscript[5, F] Subscript[10, F]*)
YukawaTerm1=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv][[1]]//Simplify;


InputInv={{2,1,1},{True,True,True}};  (*H^+ Subscript[10, F] Subscript[10, F]*)
YukawaTerm2=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv][[1]]//Simplify;


Ysff=-GradYukawa[yt1*YukawaTerm1+yt2*YukawaTerm2];


YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt1>0,yt2>0}]];


(* ::Section:: *)
(*Dimensional Reduction*)


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->True];
PerformDRhard[]


PrintCouplings[]


BetaFunctions4D[]


PrintTemporalScalarCouplings[]


PrintScalarMass["LO"]
PrintScalarMass["NLO"]


PrintDebyeMass["LO"]
PrintDebyeMass["NLO"]


(* ::Section:: *)
(*Comparing with the literature *)


GradH[f_,x_List?VectorQ]:=D[f,{x}];


(*In the litterature it's common to represent the adjoint representation as \[CapitalPhi]=T^a\[Phi]^2 where T^a are traceless hermitian matrices
satisfying Tr T^aT^b=1/2 \[Delta]^2*)


(*Yet groupmath uses a different basis. So we need to find the relation between them.*)


<<GroupMath` (*Here we use GroupMath to find the hermitian matrices. Yet the user could very-well supply these themselves*)


Mat=RepMatrices[SU5,{1,0,0,0}]; (*Defines 5x5 matrices of the fundamental rep. These are equivalent to the T^a above*)


\[Phi]Mat=Sum[Mat[[i]]\[Phi][i],{i,1,Length[Mat]}]; (*\[CapitalPhi]=T^a\[Phi]^2*)


var=Normal[\[Phi]Mat]//Variables;
QuarticMat=\[Delta]1*Tr[\[Phi]Mat . \[Phi]Mat]^2+\[Delta]2 Tr[\[Phi]Mat . \[Phi]Mat . \[Phi]Mat . \[Phi]Mat]//Simplify; (*\[Delta]1 (Tr\[CapitalPhi]^2)^2+\[Delta]2 Tr \[CapitalPhi]^4*)
QuarticA=GradH[QuarticMat,var]//GradH[#,var]&//GradH[#,var]&//GradH[#,var]&//SparseArray;


(*Now we create various invariants in the basis above, and the DRalgo basis.*)


v1=Tr[TensorContract[QuarticA,{1,2}]]-Tr[TensorContract[\[Lambda]4,{1,2}]]//Simplify;
v2=Total[QuarticA*QuarticA,-1]-Total[\[Lambda]4 \[Lambda]4,-1]//Simplify;
v3=Tr[Flatten[QuarticA,{{1,2},{3,4}}] . Flatten[QuarticA,{{1,2},{3,4}}] . Flatten[QuarticA,{{1,2},{3,4}}]]-Tr[Flatten[\[Lambda]4,{{1,2},{3,4}}] . Flatten[\[Lambda]4,{{1,2},{3,4}}] . Flatten[\[Lambda]4,{{1,2},{3,4}}]]//Simplify;


(*Demanding that the invariants are equal gives the relation between couplings*)


(*These replacement rules singles out pure adjoint quartics, and removes mixing to the 10plt*)


v1=v1/.\[Lambda]M1->0/.\[Lambda]M2->0/.\[Lambda]->0;
v2=v2/.\[Lambda]M1->0/.\[Lambda]M2->0/.\[Lambda]->0;
v3=v3/.\[Lambda]M1->0/.\[Lambda]M2->0/.\[Lambda]->0;


Solve[v1==0&&v2==0&&v3==0,{\[Lambda]H,\[Lambda]S}]//Simplify


(* ::Section:: *)
(*Other Types of Representations*)


(*Let's now see how to transform between DRalgo's, and conventions in the litterature, for other reps*)


(* ::Subtitle:: *)
(*Antisymmetric *)


(*First Let's define the quartics in DRalgo*)


Group={"SU5"};
RepAdjoint={{1,0,0,1}};
HiggsSymmetric={{{0,1,0,0}},"C"}; (*This is the symmetric rep*)
RepScalar={HiggsSymmetric};
CouplingName={g};


RepFermion3Gen={}//Flatten[#,1]&;


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJC,\[Mu]IJ,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion3Gen,RepScalar];


InputInv={{1,1,1,1},{True,True,False,False}}; (*(Tr\[CapitalPhi]\[CapitalPhi]^+)^2+(Tr\[CapitalPhi]\[CapitalPhi]^+\[CapitalPhi]\[CapitalPhi]^+)*)
QuarticTemp=CreateInvariant[Group,RepScalar,InputInv];


VQuartic=\[Lambda]H*QuarticTemp[[1]]+\[Lambda]S*QuarticTemp[[2]];


\[Lambda]4=GradQuartic[VQuartic];


(*Now let's define the same term with another convention. The 10 of SU(5) is a two-index symmetric tensor*)


n=5;


Sym=SymmetrizedArray[{i_,j_}->A[i,j],{n,n},Symmetric[{1,2}]];


Sym=SymmetrizedArray[{i_,j_}->A[i,j],{n,n},Antisymmetric[{1,2}]];


SymMod=Sym-0*IdentityMatrix[5]/5*Tr[Sym]; (*Makes sure that the matrix is symmetric*)


SymRe=SparseArray[Normal[SymMod]//ReplaceAll[#,A[x___]:>(AR[x]+I AI[x])/Sqrt[2]]&];
SymIm=SparseArray[Normal[SymMod]//ReplaceAll[#,A[x___]:>(AR[x]-I AI[x])/Sqrt[2]]&];
Vars=SymRe//Normal//Variables;
varAssump=#>0&/@Vars;


V=(
	+\[Delta]1*Tr[SymRe . SymIm]^2
	+\[Delta]2*Tr[SymRe . SymIm . SymRe . SymIm]
	);


QuarticA=GradH[V,Vars]//GradH[#,Vars]&//GradH[#,Vars]&//GradH[#,Vars]&//Normal//Simplify//SparseArray;


(*Now let's compare invariants like before*)


v1=Tr[TensorContract[QuarticA,{1,2}]]-Tr[TensorContract[\[Lambda]4,{1,2}]]//Simplify;


v2=Total[QuarticA*QuarticA,-1]-Total[\[Lambda]4*\[Lambda]4,-1]//Simplify;


v3=Tr[Flatten[QuarticA,{{1,2},{3,4}}] . Flatten[QuarticA,{{1,2},{3,4}}] . Flatten[QuarticA,{{1,2},{3,4}}]]-Tr[Flatten[\[Lambda]4,{{1,2},{3,4}}] . Flatten[\[Lambda]4,{{1,2},{3,4}}] . Flatten[\[Lambda]4,{{1,2},{3,4}}]]//Simplify;


Solve[v1==0&&v2==0&&v3==0,{\[Lambda]H,\[Lambda]S}]//Simplify
