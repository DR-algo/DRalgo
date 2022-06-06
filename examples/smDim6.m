(* ::Package:: *)

(*Quit[];*)


SetDirectory[NotebookDirectory[]];
$LoadGroupMath=True;
<<../DRalgo.m


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


Rep1={{{1,0},{1},1/6},"L"};
Rep2={{{1,0},{0},2/3},"R"};
Rep3={{{1,0},{0},-1/3},"R"};
Rep4={{{0,0},{1},-1/2},"L"};
Rep5={{{0,0},{0},-1},"R"};
RepFermion1Gen={Rep1,Rep2,Rep3,Rep4,Rep5};


(*Rep1={{{1,0},{1},Yu/2},"L"};
Rep2={{{1,0},{0},Yd/2},"R"};
Rep3={{{1,0},{0},Yd/2},"R"};
Rep4={{{0,0},{1},Yl/2},"L"};
Rep5={{{0,0},{0},Ye/2},"R"};
RepFermion1Gen={Rep1,Rep2,Rep3,Rep4,Rep5};*)


RepFermion3Gen={RepFermion1Gen}//Flatten[#,1]&;


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


Dim6Term=MassTerm1^3;


\[Lambda]6=c6 GradSextic[Dim6Term];


(* ::Section:: *)
(*Dimensional Reduction*)


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False,Mode->3];
PosFermion=PrintFermionRepPositions[];
FermionMat=Table[{nF,i},{i,PosFermion}];
DefineNF[FermionMat]
PerformDRhard[]


PrintCouplings[]


PrintScalarMass["LO"]
PrintScalarMass["NLO"]


PrintDebyeMass["LO"]
PrintDebyeMass["NLO"]


PrintTemporalScalarCouplings[]


BetaFunctions4D[]


PrintPressure["LO"]
PrintPressure["NLO"]
PrintPressure["NNLO"]


BetaFunctions4D[]


PosScalar=PrintScalarRepPositions[];
PosVector=PrintGaugeRepPositions[];
PosFermion=PrintFermionRepPositions[];


Table[AnomDim4D["S",{a,b}],{a,PosScalar},{b,PosScalar}]
Table[AnomDim4D["V",{a,b}],{a,PosVector},{b,PosVector}]
Table[AnomDim4D["F",{a,b}],{a,PosFermion},{b,PosFermion}]


(* ::Text:: *)
(*One active doublets:*)


PerformDRsoft[{}];


PrintCouplingsUS[]


PrintScalarMassUS["LO"]
PrintScalarMassUS["NLO"]


BetaFunctions3DUS[]


(* ::Text:: *)
(*2 - Loop Effective potential*)


DefineNewTensorsUS[\[Mu]ij,\[Lambda]4,\[Lambda]3,gvss,gvvv]; \[CurlyPhi]vev={0,0,0,\[CurlyPhi]}//SparseArray; DefineVEVS[\[CurlyPhi]vev];
PrintTensorsVEV[];


(*The vector-mass matrix is not diagonal. This can be seen by looking at*)


PrintTensorsVEV[][[2]]//Normal


(*So the problem is the A^3-B mixing as usual.*)


(*We can diagonalize the mass-matrix via*)


MassMatrix=PrintTensorsVEV[];
VectorMass=MassMatrix[[2]]//Normal;
VectorEigenvectors=FullSimplify[
    Transpose[Normalize/@Eigenvectors[VectorMass[[11;;12,11;;12]]]],
Assumptions->{gw>0,gY>0,\[CurlyPhi]>0}]; DVRot={{IdentityMatrix[10],0},{0,VectorEigenvectors}}//ArrayFlatten; DSRot=IdentityMatrix[4];
RotateTensorsUSPostVEV[DSRot,DVRot];


(*We now see the new diagonal masses by writing*)


PrintTensorsVEV[][[2]]//Normal//Simplify;


(*The effective potential is now given by*)


CalculatePotentialUS[];


PrintEffectivePotential["LO"]
PrintEffectivePotential["NLO"]
PrintEffectivePotential["NNLO"]





(* ::Title:: *)
(*Testing remove !*)


Habij=Transpose[Activate @ TensorContract[
        Inactive[TensorProduct][gvss,gvss], { {3, 5}}],{1,3,2,4}]//SparseArray;   
        
HabijV=Habij+Transpose[Habij,{2,1,3,4}]//SparseArray;


YTemp=Transpose[Transpose[Simplify[ Activate @ TensorContract[
        Inactive[TensorProduct][Ysff,YsffC], {{3, 5}}]],{1,3,2,4}],{1,2,4,3}]//SparseArray;


\[CapitalLambda]\[Lambda]6tem=Flatten[\[Lambda]4 . \[Lambda]4,{{1},{2},{4},{5},{3,6}}];
\[CapitalLambda]\[Lambda]6tot=\[CapitalLambda]\[Lambda]6tem . Flatten[\[Lambda]4,{1,2}];
Prefac=-Zeta[3]/(128 \[Pi]^4 T^2);
SymHelp=Table[\[CapitalLambda]\[Lambda]6tot["NonzeroPositions"][[i]]->\[CapitalLambda]\[Lambda]6tot["NonzeroValues"][[i]],{i,1,Length[\[CapitalLambda]\[Lambda]6tot["NonzeroPositions"]]}];
ContriSS=Prefac*SymmetrizedArray[SymHelp,Dimensions[\[CapitalLambda]\[Lambda]6tot],Symmetric[{1,2,3,4,5,6}]]//SparseArray//SimplifySparse;


\[CapitalLambda]\[Lambda]6tem=Flatten[Transpose[HabijV,{1,4,3,2}] . HabijV,{{2},{3},{5},{6},{1,4}}];
\[CapitalLambda]\[Lambda]6tot=\[CapitalLambda]\[Lambda]6tem . Flatten[HabijV,{{1,2},{3},{4}}];
Prefac=-3*Zeta[3]/(128 \[Pi]^4 T^2);
SymHelp=Table[\[CapitalLambda]\[Lambda]6tot["NonzeroPositions"][[i]]->\[CapitalLambda]\[Lambda]6tot["NonzeroValues"][[i]],{i,1,Length[\[CapitalLambda]\[Lambda]6tot["NonzeroPositions"]]}];
ContriVV=Prefac*SymmetrizedArray[SymHelp,Dimensions[\[CapitalLambda]\[Lambda]6tot],Symmetric[{1,2,3,4,5,6}]]//SparseArray//SimplifySparse;


\[CapitalLambda]\[Lambda]6tot1=Table[Tr[a . b . c . d . e . f],{a,Ysff},{b,YsffC},{c,Ysff},{d,YsffC},{e,Ysff},{f,YsffC}]//SparseArray;
\[CapitalLambda]\[Lambda]6totC=Table[Tr[a . b . c . d . e . f],{a,YsffC},{b,Ysff},{c,YsffC},{d,Ysff},{e,YsffC},{f,Ysff}]//SparseArray;
\[CapitalLambda]\[Lambda]6tot=\[CapitalLambda]\[Lambda]6tot1+ \[CapitalLambda]\[Lambda]6totC;
Prefac=((7 Zeta[3])/(64 \[Pi]^4 T^2));
SymHelp=Table[\[CapitalLambda]\[Lambda]6tot["NonzeroPositions"][[i]]->\[CapitalLambda]\[Lambda]6tot["NonzeroValues"][[i]],{i,1,Length[\[CapitalLambda]\[Lambda]6tot["NonzeroPositions"]]}];
ContriFF=Prefac*SymmetrizedArray[SymHelp,Dimensions[\[CapitalLambda]\[Lambda]6tot],Symmetric[{1,2,3,4,5,6}]]//SparseArray//SimplifySparse;

ContriSETemp=-0*ZijS . \[Lambda]6;
ContriSE=ContriSETemp+Transpose[ContriSETemp,{2,1,3,4,5,6}]+Transpose[ContriSETemp,{3,1,2,4,5,6}]+Transpose[ContriSETemp,{4,1,2,3,5,6}]+Transpose[ContriSETemp,{5,1,2,3,4,6}]+Transpose[ContriSETemp,{6,1,2,3,5,4}]//SimplifySparse;

(*Minus sign from matching*)
\[Lambda]6D=- ContriSS- ContriVV- ContriFF//SimplifySparse; 


\[CapitalLambda]\[Lambda]6tem=Flatten[Transpose[HabijV,{1,4,3,2}] . HabijV,{{2},{3},{5},{6},{1,4}}];
\[CapitalLambda]\[Lambda]6tot=\[CapitalLambda]\[Lambda]6tem . Flatten[HabijV,{{1,2},{3},{4}}];
Prefac=-3*Zeta[3]/(128 \[Pi]^4 T^2);
SymHelp=Table[\[CapitalLambda]\[Lambda]6tot["NonzeroPositions"][[i]]->\[CapitalLambda]\[Lambda]6tot["NonzeroValues"][[i]],{i,1,Length[\[CapitalLambda]\[Lambda]6tot["NonzeroPositions"]]}];
ContriVV=Prefac*SymmetrizedArray[SymHelp,Dimensions[\[CapitalLambda]\[Lambda]6tot],Symmetric[{1,2,3,4,5,6}]]//SparseArray//SimplifySparse;



\[Lambda]6D//Flatten[#]&//DeleteDuplicates[#]&//Sort//Simplify


GradH[f_,x_List?VectorQ]:=D[f,{x}];


AEHelp=Dim6Term//Variables;


\[Lambda]6=c6*SparseArray[GradH[Dim6Term,AEHelp]]//GradH[#,AEHelp]&//GradH[#,AEHelp]&//GradH[#,AEHelp]&//GradH[#,AEHelp]&//GradH[#,AEHelp]&;


Solve[Normal[ContriFF]==Normal[\[Lambda]6],c6]


\[Lambda]6//Flatten[#]&//DeleteDuplicates[#]&//Sort//Simplify


(*
	Trick to simplify tensors because Mathematica 13 sucks.
*)
SimplifySparse[s_SparseArray] := With[
    {
    elem =Simplify[s["NonzeroValues"]],
    pos=s["NonzeroPositions"],
    dim = Dimensions[s]
    },
SparseArray[pos->elem,dim,0]
    ]

