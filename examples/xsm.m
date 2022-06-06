(* ::Package:: *)

(*Quit[];*)


SetDirectory[NotebookDirectory[]];
$LoadGroupMath=True;
<<../DRalgo.m


(* ::Chapter:: *)
(*SM+sr1*)


(*see 2102.11145 [hep-ph]*)


(* ::Section::Closed:: *)
(*Model*)


Group={"SU3","SU2","U1"};
RepAdjoint={{1,1},{2},0};
scalar1={{{0,0},{1},Y\[Phi]/2},"C"};
scalar2={{{0,0},{0},0},"R"};
RepScalar={scalar1,scalar2};
CouplingName={g3,g2,g1};


Rep1={{{1,0},{1},Yq/2},"L"};
Rep2={{{1,0},{0},Yu/2},"R"};
Rep3={{{1,0},{0},Yd/2},"R"};
Rep4={{{0,0},{1},Yl/2},"L"};
Rep5={{{0,0},{0},Ye/2},"R"};
RepFermion1Gen={Rep1,Rep2,Rep3,Rep4,Rep5};


(* ::Text:: *)
(*The input for the gauge interactions to DRalgo are then given by*)


RepFermion1Gen={Rep1,Rep2,Rep3,Rep4,Rep5};
RepFermion3Gen={RepFermion1Gen}//Flatten[#,1]&;


(* ::Text:: *)
(*The first element is the vector self-interaction matrix:*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion3Gen,RepScalar];


InputInv={{1,1},{True,False}};
MassTerm1=CreateInvariant[Group,RepScalar,InputInv]//Simplify//FullSimplify;
InputInv={{2,2},{True,True}};
MassTerm2=CreateInvariant[Group,RepScalar,InputInv]//Simplify//FullSimplify;


VMass=(
	+m1*MassTerm1
	+\[Mu]\[Sigma]/2*MassTerm2
	);


\[Mu]ij=GradMass[VMass[[1]]]//Simplify//SparseArray;


QuarticTerm1=MassTerm1[[1]]^2;
QuarticTerm2=MassTerm2[[1]]^2;
QuarticTerm3=MassTerm1[[1]]*MassTerm2[[1]];


VQuartic=(
	+\[Lambda]1H*QuarticTerm1
	+\[Lambda]\[Sigma]/4*QuarticTerm2
	+\[Lambda]m/2*QuarticTerm3
	);


\[Lambda]4=GradQuartic[VQuartic];


InputInv={{1,1,2},{True,False,True}};
CubicTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{2,2,2},{True,True,True}};
CubicTerm2=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;


VCubic=(
	+\[Mu]m/2*CubicTerm1
	+\[Mu]3/3*CubicTerm2
	);


\[Lambda]3=GradCubic[VCubic];


InputInv={{2},{True}};
TadpoleTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;


VTadpole=\[Mu]1*TadpoleTerm1;


\[Lambda]1=GradTadpole[VTadpole];


InputInv={{1,1,2},{False,False,True}}; 
YukawaDoublet1=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv]//Simplify;


Ysff=-yt1*GradYukawa[YukawaDoublet1[[1]]];


YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt1>0}]];


(* ::Section:: *)
(*Dimensional Reduction*)


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PosFermion=PrintFermionRepPositions[];
FermionMat=Table[{nF,i},{i,PosFermion}];
DefineNF[FermionMat]
PerformDRhard[]


BetaFunctions4D[]
PrintCouplings[]//Simplify


(4\[Pi])^2 m1/.BetaFunctions4D[]


(4\[Pi])^2  \[Mu]m3d/.PrintCouplings[]
1/T^(3/2) (4\[Pi])^2  \[Mu]33d/.PrintCouplings[]//Expand
1/T^(2/2) (4\[Pi])^2  \[Lambda]\[Sigma]3d/.PrintCouplings[]//Expand
1/T^(2/2) (4\[Pi])^2  \[Lambda]m3d/.PrintCouplings[]//Expand


PrintTadpoles["LO"]
\[Mu]13d/.PrintTadpoles["NLO"]


PrintTemporalScalarCouplings[]


PrintDebyeMass["LO"]
PrintDebyeMass["NLO"]


PrintScalarMass["LO"]//Simplify
PrintScalarMass["NLO"]//Simplify


PerformDRsoft[{}];


PrintCouplingsUS[]


PrintScalarMassUS["LO"]
PrintScalarMassUS["NLO"]


BetaFunctions3DUS[]


PerformDRsoft[{5}];


PrintCouplingsUS[]


PrintScalarMassUS["LO"]
PrintScalarMassUS["LO"]


BetaFunctions3DUS[]


(* ::Section:: *)
(*Effective potential*)


DefineNewTensorsUS[\[Mu]ij,\[Lambda]4,\[Lambda]3,gvss,gvvv]; \[CurlyPhi]vev={0,0,0,\[CurlyPhi],s}//SparseArray; DefineVEVS[\[CurlyPhi]vev];
PrintTensorsVEV[];


(*The vector-mass matrix is not diagonal. Let's first diagonalize this matrix*)


MassMatrix=PrintTensorsVEV[];
VectorMass=MassMatrix[[2]]//Normal;
VectorEigenvectors=Simplify[
    Transpose[Normalize/@Eigenvectors[VectorMass[[11;;12,11;;12]]]],
Assumptions->{g1>0,g2>0,\[CurlyPhi]>0,Y\[Phi]>0}]; DVRot={{IdentityMatrix[10],0},{0,VectorEigenvectors}}//ArrayFlatten//Simplify;


(*We now need to diagonalize the scalar-mass matrix*)


PrintTensorsVEV[][[1]]//Normal


ScalarMass=MassMatrix[[1]]//Normal;


(*We see that the scalar-masses are quite complicated. So it is easier to do the rotation to the diagonal basis by hand.*)


(*If we only want the 1-loop effective potential it suffices to use the following method:*)


DSRot={{IdentityMatrix[5]}}//ArrayFlatten//Simplify;


\[Mu]vecDia=Simplify[Transpose[DVRot] . VectorMass . DVRot,
Assumptions->{g1>0,g2>0,\[CurlyPhi]>0,Y\[Phi]>0}];


\[Mu]ijDiag={{ScalarMass[[1;;3,1;;3]],0},{0,DiagonalMatrix[Eigenvalues[ScalarMass[[4;;5,4;;5]]]]}}//ArrayFlatten//Simplify//FullSimplify;


DiagonalMatrixQAE[\[Mu]ijDiag]


RotateTensorsCustomMass[DSRot,DVRot,\[Mu]vecDia,\[Mu]ijDiag];


CalculatePotentialUS[];


PrintEffectivePotential["LO"]


PrintEffectivePotential["NLO"]


(*If we want the 2-loop effective potential we also need the scalar-mass diagonalization matrix:*)


(*Since the masses are quite complicated it is useful to use the replacement*)


HelpReplace={ScalarMass[[4;;5,4;;5]][[1,1]]->M11,ScalarMass[[4;;5,4;;5]][[1,2]]->M12,ScalarMass[[4;;5,4;;5]][[2,1]]->M12,ScalarMass[[4;;5,4;;5]][[2,2]]->M22};


ScalarMassHelp=ScalarMass[[4;;5,4;;5]]//ReplaceAll[#,HelpReplace]&;


HelpRevert=HelpReplace/.(a_->b_)->(b->a);


ScalarEigenvectors=FullSimplify[
    Transpose[Normalize/@Eigenvectors[ScalarMassHelp]]];


(*Since we know that the matrix should be SO(2), let's write it as*)


ScalDiaMatrix={{c\[Theta],-s\[Theta]},{s\[Theta],c\[Theta]}};


(*We can then write the total rotation matrix as (first 3 components are already diagonal)*)


 DSRot={{IdentityMatrix[3],0},{0,ScalDiaMatrix}}//ArrayFlatten


(*We might as well simplify the scalar masses as well*)


DiagonalMatrix[Eigenvalues[ScalarMass[[4;;5,4;;5]]]]


\[Mu]ijDiag={{ScalarMass[[1;;3,1;;3]],0},{0,DiagonalMatrix[{\[Mu]S11,\[Mu]S22}]}}//ArrayFlatten//Simplify//FullSimplify


(*And finally we get the effective potential from*)


RotateTensorsCustomMass[DSRot,DVRot,\[Mu]vecDia,\[Mu]ijDiag];


CalculatePotentialUS[];


PrintEffectivePotential["LO"]/.{c\[Theta]->1,s\[Theta]->0}//Simplify (*The replacement rule is just because the rotation should not act on the tree-level potential*)


PrintEffectivePotential["NLO"]


PrintEffectivePotential["NNLO"]
