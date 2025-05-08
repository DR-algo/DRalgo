(* ::Package:: *)

Quit[];


If[$InputFileName=="",
	SetDirectory[NotebookDirectory[]],
	SetDirectory[DirectoryName[$InputFileName]]
];
(*Put this if you want to create multiple model-files with the same kernel*)
(*DRalgo`$GroupMathMultipleModels=True;*)

DRalgo`$LoadGroupMath=True;
DRalgo`$InstallGroupMath=True;

Check[
    Get["../DRalgo.m"],
    Message[Get::noopen, "DRalgo` at "<>ToString[$UserBaseDirectory]<>"/Applications"];
    Abort[];
]


(* ::Chapter:: *)
(*SM+sr1*)


(*see 2102.11145 [hep-ph]*)


(* ::Section::Closed:: *)
(*Model*)


Group={"SU3","SU2","U1"};
RepAdjoint={{1,1},{2},0};
HiggsDoublet1={{{0,0},{1},1/2},"C"};
HiggsSinglet1={{{0,0},{0},0},"R"};
HiggsSinglet2={{{0,0},{0},0},"R"};
RepScalar={HiggsDoublet1,HiggsSinglet1,HiggsSinglet2};
CouplingName={g3,g2,g1};


Rep1={{{1,0},{1},1/6},"L"};
Rep2={{{1,0},{0},2/3},"R"};
Rep3={{{1,0},{0},-1/3},"R"};
Rep4={{{0,0},{1},-1/2},"L"};
Rep5={{{0,0},{0},-1},"R"};
RepFermion1Gen={Rep1,Rep2,Rep3,Rep4,Rep5};


RepFermion3Gen={RepFermion1Gen,RepFermion1Gen,RepFermion1Gen}//Flatten[#,1]&;


(* ::Text:: *)
(*The input for the gauge interactions toDRalgo are then given by*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion3Gen,RepScalar];


(* ::Text:: *)
(*The first element is the vector self - interaction matrix :*)


InputInv={{1,1},{True,False}}; (*\[Phi]\[Phi]^+*)
MassTerm1=CreateInvariant[Group,RepScalar,InputInv]//Simplify;
InputInv={{2,2},{True,True}}; (*S1^2*)
MassTerm2=CreateInvariant[Group,RepScalar,InputInv]//Simplify;
InputInv={{3,3},{True,True}}; (*S2^2*)
MassTerm3=CreateInvariant[Group,RepScalar,InputInv]//Simplify;


VMass=(
	+m1*MassTerm1
	+1/2*mS1*MassTerm2
	+1/2*mS2*MassTerm3
	);


\[Mu]ij=GradMass[VMass[[1]]]//Simplify; 


QuarticTerm1=MassTerm1[[1]]^2; (*[(\[Phi]\[Phi]^+)]^2*)
QuarticTerm2=MassTerm2[[1]]^2; (* S1^4*)
QuarticTerm3=MassTerm3[[1]]^2; (* S2^4*)
QuarticTerm4=MassTerm1[[1]]*MassTerm2[[1]]; (* [\[Phi]\[Phi]^+]S1^2*)
QuarticTerm5=MassTerm1[[1]]*MassTerm3[[1]]; (* [\[Phi]\[Phi]^+]S2^2*)
QuarticTerm6=MassTerm2[[1]]*MassTerm3[[1]]; (* S1^2S2^2*)


InputInv={{2,2,2,3},{True,True,True,True}}; (*S1^3S2*)
QuarticTerm7=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{3,3,3,2},{True,True,True,True}}; (*S2^3S1*)
QuarticTerm8=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;


VQuartic=(
	+\[Lambda]H*QuarticTerm1
	+\[Lambda]S1/4!*QuarticTerm2
	+\[Lambda]S2/4!*QuarticTerm3
	+\[Lambda]M1/2*QuarticTerm4
	+\[Lambda]M2/2*QuarticTerm5
	+\[Lambda]M3/2*QuarticTerm6
	+\[Lambda]M4/4!*QuarticTerm7
	+\[Lambda]M5/4!*QuarticTerm8
	);


\[Lambda]4=GradQuartic[VQuartic];


InputInv={{1,1,2},{True,False,True}}; (*\[Phi]\[Phi]^+S1*)
CubicTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{1,1,3},{True,False,True}}; (*\[Phi]\[Phi]^+S2*)
CubicTerm2=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{2,2,2},{True,True,True}}; (*S1^3*)
CubicTerm3=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{3,3,3},{True,True,True}}; (*S2^3*)
CubicTerm4=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{2,3,3},{True,True,True}}; (*S1*S2^2*)
CubicTerm5=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{2,2,3},{True,True,True}}; (*S1^2*S2*)
CubicTerm6=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;


VCubic=(
	+\[Lambda]C1*CubicTerm1
	+\[Lambda]C2*CubicTerm2
	+\[Lambda]C3/3!*CubicTerm3
	+\[Lambda]C4/3!*CubicTerm4
	+\[Lambda]C5/2!*CubicTerm5
	+\[Lambda]C6/2!*CubicTerm6
	);


\[Lambda]3=GradCubic[VCubic];


InputInv={{2},{True}}; (*S1*)
TadpoleTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;
InputInv={{3},{True}}; (*S2*)
TadpoleTerm2=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;


VTadpole=(
	+\[Lambda]Tad1*TadpoleTerm1
	+\[Lambda]Tad2*TadpoleTerm2
	);


\[Lambda]1=GradTadpole[VTadpole];


InputInv={{1,1,2},{False,False,True}}; 
YukawaDoublet1=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv]//Simplify;


Ysff=-yt*GradYukawa[YukawaDoublet1[[1]]];


YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt>0}]];


(* ::Section:: *)
(*Dimensional Reduction*)


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->True];
PerformDRhard[]


(* ::Section:: *)
(*Tests*)


(* ::Subsection:: *)
(*Test hard*)


testList={};


(* ::Subsubsection:: *)
(*Beta functions*)


AppendTo[testList,
TestCreate[BetaFunctions4D[],
	{g1^2->(41 g1^4)/(48 \[Pi]^2),g2^2->-((19 g2^4)/(48 \[Pi]^2)),g3^2->-((7 g3^4)/(8 \[Pi]^2)),\[Lambda]H->(3 g1^4+9 g2^4+6 g1^2 (g2^2-4 \[Lambda]H)-72 g2^2 \[Lambda]H+96 \[Lambda]H (yt^2+2 \[Lambda]H)+4 (-12 yt^4+\[Lambda]M1^2+\[Lambda]M2^2))/(128 \[Pi]^2),\[Lambda]M1->(-3 g1^2 \[Lambda]M1-9 g2^2 \[Lambda]M1+4 \[Lambda]M2 \[Lambda]M3+2 \[Lambda]M1 (6 yt^2+12 \[Lambda]H+4 \[Lambda]M1+\[Lambda]S1))/(32 \[Pi]^2),\[Lambda]M2->(-3 g1^2 \[Lambda]M2-9 g2^2 \[Lambda]M2+4 \[Lambda]M1 \[Lambda]M3+2 \[Lambda]M2 (6 yt^2+12 \[Lambda]H+4 \[Lambda]M2+\[Lambda]S2))/(32 \[Pi]^2),\[Lambda]M3->(32 \[Lambda]M1 \[Lambda]M2+\[Lambda]M4^2+\[Lambda]M4 \[Lambda]M5+\[Lambda]M5^2+16 \[Lambda]M3 (8 \[Lambda]M3+\[Lambda]S1+\[Lambda]S2))/(256 \[Pi]^2),\[Lambda]M4->(3 (2 \[Lambda]M3 (2 \[Lambda]M4+\[Lambda]M5)+\[Lambda]M4 \[Lambda]S1))/(16 \[Pi]^2),\[Lambda]M5->(3 (2 \[Lambda]M3 (\[Lambda]M4+2 \[Lambda]M5)+\[Lambda]M5 \[Lambda]S2))/(16 \[Pi]^2),\[Lambda]S1->(3 (32 \[Lambda]M1^2+32 \[Lambda]M3^2+\[Lambda]M4^2+8 \[Lambda]S1^2))/(128 \[Pi]^2),\[Lambda]S2->(3 (32 \[Lambda]M2^2+32 \[Lambda]M3^2+\[Lambda]M5^2+8 \[Lambda]S2^2))/(128 \[Pi]^2),\[Lambda]C1->(-3 \[Lambda]C1 (g1^2+3 g2^2-4 (yt^2+2 \[Lambda]H))+2 (4 \[Lambda]C1+\[Lambda]C3) \[Lambda]M1+2 \[Lambda]C5 \[Lambda]M2)/(32 \[Pi]^2),\[Lambda]C2->(-3 \[Lambda]C2 (g1^2+3 g2^2-4 (yt^2+2 \[Lambda]H))+2 \[Lambda]C6 \[Lambda]M1+2 (4 \[Lambda]C2+\[Lambda]C4) \[Lambda]M2)/(32 \[Pi]^2),\[Lambda]C3->(3 (8 \[Lambda]C1 \[Lambda]M1+4 \[Lambda]C5 \[Lambda]M3+\[Lambda]C6 \[Lambda]M4+2 \[Lambda]C3 \[Lambda]S1))/(32 \[Pi]^2),\[Lambda]C4->(3 (8 \[Lambda]C2 \[Lambda]M2+4 \[Lambda]C6 \[Lambda]M3+\[Lambda]C5 \[Lambda]M5+2 \[Lambda]C4 \[Lambda]S2))/(32 \[Pi]^2),\[Lambda]C5->(8 \[Lambda]C1 \[Lambda]M2+4 \[Lambda]C3 \[Lambda]M3+\[Lambda]C6 \[Lambda]M4+(\[Lambda]C4+\[Lambda]C6) \[Lambda]M5+2 \[Lambda]C5 (8 \[Lambda]M3+\[Lambda]S2))/(32 \[Pi]^2),\[Lambda]C6->(8 \[Lambda]C2 \[Lambda]M1+4 \[Lambda]C4 \[Lambda]M3+(\[Lambda]C3+\[Lambda]C5) \[Lambda]M4+\[Lambda]C5 \[Lambda]M5+2 \[Lambda]C6 (8 \[Lambda]M3+\[Lambda]S1))/(32 \[Pi]^2),yt->(yt (-17 g1^2-27 g2^2-96 g3^2+54 yt^2))/(192 \[Pi]^2),\[Lambda]Tad1->(4 m1 \[Lambda]C1+mS1 \[Lambda]C3+mS2 \[Lambda]C5)/(16 \[Pi]^2),\[Lambda]Tad2->(4 m1 \[Lambda]C2+mS2 \[Lambda]C4+mS1 \[Lambda]C6)/(16 \[Pi]^2)}
]];


(* ::Subsubsection:: *)
(*Dimension 4 matching relations*)


AppendTo[testList,
TestCreate[PrintCouplings[],
	{g13d^2->g1^2 T-(g1^4 (Lb+40 Lf) T)/(96 \[Pi]^2),g23d^2->g2^2 T+(g2^4 (4+43 Lb-24 Lf) T)/(96 \[Pi]^2),g33d^2->g3^2 T+(g3^4 (1+11 Lb-4 Lf) T)/(16 \[Pi]^2),\[Lambda]H3d->1/(256 \[Pi]^2) T (g2^4 (6-9 Lb)+g1^4 (2-3 Lb)+72 g2^2 Lb \[Lambda]H+256 \[Pi]^2 \[Lambda]H+g1^2 (g2^2 (4-6 Lb)+24 Lb \[Lambda]H)+48 Lf (yt^4-2 yt^2 \[Lambda]H)-4 Lb (48 \[Lambda]H^2+\[Lambda]M1^2+\[Lambda]M2^2)),\[Lambda]M13d->1/(64 \[Pi]^2) T (3 g1^2 Lb \[Lambda]M1+9 g2^2 Lb \[Lambda]M1-2 (-32 \[Pi]^2 \[Lambda]M1+6 Lf yt^2 \[Lambda]M1+2 Lb \[Lambda]M2 \[Lambda]M3+Lb \[Lambda]M1 (12 \[Lambda]H+4 \[Lambda]M1+\[Lambda]S1))),\[Lambda]M23d->1/(64 \[Pi]^2) T (3 g1^2 Lb \[Lambda]M2+9 g2^2 Lb \[Lambda]M2-2 (-32 \[Pi]^2 \[Lambda]M2+6 Lf yt^2 \[Lambda]M2+2 Lb \[Lambda]M1 \[Lambda]M3+Lb \[Lambda]M2 (12 \[Lambda]H+4 \[Lambda]M2+\[Lambda]S2))),\[Lambda]M33d->T \[Lambda]M3-(Lb T (32 \[Lambda]M1 \[Lambda]M2+\[Lambda]M4^2+\[Lambda]M4 \[Lambda]M5+\[Lambda]M5^2+16 \[Lambda]M3 (8 \[Lambda]M3+\[Lambda]S1+\[Lambda]S2)))/(512 \[Pi]^2),\[Lambda]M43d->T \[Lambda]M4-(3 Lb T (2 \[Lambda]M3 (2 \[Lambda]M4+\[Lambda]M5)+\[Lambda]M4 \[Lambda]S1))/(32 \[Pi]^2),\[Lambda]M53d->T \[Lambda]M5-(3 Lb T (2 \[Lambda]M3 (\[Lambda]M4+2 \[Lambda]M5)+\[Lambda]M5 \[Lambda]S2))/(32 \[Pi]^2),\[Lambda]S13d->T (\[Lambda]S1-(3 Lb (32 \[Lambda]M1^2+32 \[Lambda]M3^2+\[Lambda]M4^2+8 \[Lambda]S1^2))/(256 \[Pi]^2)),\[Lambda]S23d->T (\[Lambda]S2-(3 Lb (32 \[Lambda]M2^2+32 \[Lambda]M3^2+\[Lambda]M5^2+8 \[Lambda]S2^2))/(256 \[Pi]^2)),\[Lambda]C13d->Sqrt[T] (\[Lambda]C1+(3 g1^2 Lb \[Lambda]C1+9 g2^2 Lb \[Lambda]C1-2 (6 Lf yt^2 \[Lambda]C1+Lb (\[Lambda]C3 \[Lambda]M1+4 \[Lambda]C1 (3 \[Lambda]H+\[Lambda]M1)+\[Lambda]C5 \[Lambda]M2)))/(64 \[Pi]^2)),\[Lambda]C23d->Sqrt[T] (\[Lambda]C2+(3 g1^2 Lb \[Lambda]C2+9 g2^2 Lb \[Lambda]C2-2 (6 Lf yt^2 \[Lambda]C2+Lb (\[Lambda]C6 \[Lambda]M1+\[Lambda]C4 \[Lambda]M2+4 \[Lambda]C2 (3 \[Lambda]H+\[Lambda]M2))))/(64 \[Pi]^2)),\[Lambda]C33d->Sqrt[T] (\[Lambda]C3-(3 Lb (8 \[Lambda]C1 \[Lambda]M1+4 \[Lambda]C5 \[Lambda]M3+\[Lambda]C6 \[Lambda]M4+2 \[Lambda]C3 \[Lambda]S1))/(64 \[Pi]^2)),\[Lambda]C43d->Sqrt[T] (\[Lambda]C4-(3 Lb (8 \[Lambda]C2 \[Lambda]M2+4 \[Lambda]C6 \[Lambda]M3+\[Lambda]C5 \[Lambda]M5+2 \[Lambda]C4 \[Lambda]S2))/(64 \[Pi]^2)),\[Lambda]C53d->Sqrt[T] (\[Lambda]C5-(Lb (8 \[Lambda]C1 \[Lambda]M2+4 \[Lambda]C3 \[Lambda]M3+\[Lambda]C6 \[Lambda]M4+(\[Lambda]C4+\[Lambda]C6) \[Lambda]M5+2 \[Lambda]C5 (8 \[Lambda]M3+\[Lambda]S2)))/(64 \[Pi]^2)),\[Lambda]C63d->Sqrt[T] (\[Lambda]C6-(Lb (8 \[Lambda]C2 \[Lambda]M1+4 \[Lambda]C4 \[Lambda]M3+(\[Lambda]C3+\[Lambda]C5) \[Lambda]M4+\[Lambda]C5 \[Lambda]M5+2 \[Lambda]C6 (8 \[Lambda]M3+\[Lambda]S1)))/(64 \[Pi]^2))}
]];


(* ::Subsubsection:: *)
(*Report*)


report=TestReport[testList]
report["ResultsDataset"]



