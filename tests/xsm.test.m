(* ::Package:: *)

Quit[];


If[$InputFileName=="",
	SetDirectory[NotebookDirectory[]],
	SetDirectory[DirectoryName[$InputFileName]]
];
(*Put this if you want to create multiple model-files with the same kernel*)
(*DRalgo`$GroupMathMultipleModels=True;*)

DRalgo`$LoadGroupMath=True;
(*DRalgo`$InstallGroupMath=True;*)

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


Dim6Term=(
	+ c61*MassTerm1[[1]]^3
	+ c62*MassTerm2[[1]]^3
	+ c63*MassTerm1[[1]]^2*MassTerm2[[1]]
	+ c64*MassTerm2[[1]]^2*MassTerm1[[1]]
	);
\[Lambda]6=GradSextic[Dim6Term];


(* ::Section:: *)
(*Dimensional Reduction*)


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->True,Mode->3];
DefineDim6[\[Lambda]6];
PosFermion=PrintFermionRepPositions[];
FermionMat=Table[{nF,i},{i,PosFermion}];
DefineNF[FermionMat]
PerformDRhard[]


(* ::Section:: *)
(*Tests*)


(* ::Subsection:: *)
(*Test hard*)


Replacements={
	Thread[{c61,c62,c63,c64}->0]
}//Flatten;


testList={};


(* ::Subsubsection::Closed:: *)
(*Beta functions*)


AppendTo[testList,
TestCreate[BetaFunctions4D[]/.Replacements//DeleteCases[#,0->0]&,
	{g1^2->(g1^2 (-18 yt1^2 (Yq-Yu)^2+(18 yt1^2+g1^2 nF (3 Yd^2+Ye^2+2 Yl^2+6 Yq^2+3 Yu^2)) Y\[Phi]^2+g1^2 Y\[Phi]^4))/(48 \[Pi]^2 Y\[Phi]^2),g2^2->(g2^4 (-43+8 nF))/(48 \[Pi]^2),g3^2->(g3^4 (-33+4 nF))/(24 \[Pi]^2),\[Lambda]1H->(9 g2^4-48 yt1^4+3 g1^4 Y\[Phi]^4+6 g2^2 (g1^2 Y\[Phi]^2-12 \[Lambda]1H)-24 g1^2 Y\[Phi]^2 \[Lambda]1H+96 \[Lambda]1H (yt1^2+2 \[Lambda]1H)+4 \[Lambda]m^2)/(128 \[Pi]^2),\[Lambda]m->(\[Lambda]m (-9 g2^2+12 yt1^2-3 g1^2 Y\[Phi]^2+24 \[Lambda]1H+8 \[Lambda]m+12 \[Lambda]\[Sigma]))/(32 \[Pi]^2),\[Lambda]\[Sigma]->(\[Lambda]m^2+9 \[Lambda]\[Sigma]^2)/(8 \[Pi]^2),\[Mu]3->(3 (6 \[Lambda]\[Sigma] \[Mu]3+\[Lambda]m \[Mu]m))/(16 \[Pi]^2),\[Mu]m->(-3 (3 g2^2-4 yt1^2+g1^2 Y\[Phi]^2-8 \[Lambda]1H) \[Mu]m+8 \[Lambda]m (\[Mu]3+\[Mu]m))/(32 \[Pi]^2),yt1->(yt1 (-9 g2^2-32 g3^2+18 yt1^2-3 g1^2 (2 Yq Yu+Y\[Phi]^2)))/(64 \[Pi]^2),m1->(-3 m1 (3 g2^2-4 yt1^2+g1^2 Y\[Phi]^2-8 \[Lambda]1H)+\[Mu]m^2+2 \[Lambda]m \[Mu]\[Sigma])/(32 \[Pi]^2),\[Mu]\[Sigma]->(4 m1 \[Lambda]m+4 \[Mu]3^2+\[Mu]m^2+6 \[Lambda]\[Sigma] \[Mu]\[Sigma])/(16 \[Pi]^2),\[Mu]1->(m1 \[Mu]m+\[Mu]3 \[Mu]\[Sigma])/(8 \[Pi]^2)}
]];


(* ::Subsubsection:: *)
(*Dimension 4 matching relations*)


AppendTo[testList,
TestCreate[PrintCouplings[]/.Replacements//Simplify,
	{g13d^2->1/(96 \[Pi]^2 Y\[Phi]^2) g1^2 T (18 Lf yt1^2 (Yq-Yu)^2+(96 \[Pi]^2-Lf (18 yt1^2+g1^2 nF (3 Yd^2+Ye^2+2 Yl^2+6 Yq^2+3 Yu^2))) Y\[Phi]^2-g1^2 Lb Y\[Phi]^4),g23d^2->g2^2 T+(g2^4 (4+43 Lb-8 Lf nF) T)/(96 \[Pi]^2),g33d^2->g3^2 T+(g3^4 (3+33 Lb-4 Lf nF) T)/(48 \[Pi]^2),\[Lambda]1H3d->1/(256 \[Pi]^2) T (48 Lf yt1^4+(2-3 Lb) (3 g2^4+2 g1^2 g2^2 Y\[Phi]^2+g1^4 Y\[Phi]^4)+256 \[Pi]^2 \[Lambda]1H+24 (3 g2^2 Lb-4 Lf yt1^2+g1^2 Lb Y\[Phi]^2) \[Lambda]1H-4 Lb (48 \[Lambda]1H^2+\[Lambda]m^2)),\[Lambda]m3d->(T \[Lambda]m (9 g2^2 Lb+64 \[Pi]^2-12 Lf yt1^2+3 g1^2 Lb Y\[Phi]^2-4 Lb (6 \[Lambda]1H+2 \[Lambda]m+3 \[Lambda]\[Sigma])))/(64 \[Pi]^2),\[Lambda]\[Sigma]3d->T \[Lambda]\[Sigma]-(Lb T (\[Lambda]m^2+9 \[Lambda]\[Sigma]^2))/(16 \[Pi]^2),\[Mu]33d->1/2 Sqrt[T] (2 \[Mu]3-(3 Lb (6 \[Lambda]\[Sigma] \[Mu]3+\[Lambda]m \[Mu]m))/(16 \[Pi]^2)),\[Mu]m3d->(Sqrt[T] (4 (16 \[Pi]^2-3 Lf yt1^2) \[Mu]m+3 Lb (3 g2^2+g1^2 Y\[Phi]^2-8 \[Lambda]1H) \[Mu]m-8 Lb \[Lambda]m (\[Mu]3+\[Mu]m)))/(64 \[Pi]^2)}
]];


AppendTo[testList,
TestCreate[PrintTadpoles["LO"]/.Replacements,
	{\[Mu]13d->(12 \[Mu]1+T^2 (\[Mu]3+\[Mu]m))/(12 Sqrt[T])}
]];
AppendTo[testList,
TestCreate[PrintTadpoles["NLO"]/.Replacements//Expand,
	{\[Mu]13d->1/(768 \[Pi]^2 Sqrt[T]) (-Lb (T^2 (-12 \[Lambda]\[Sigma] \[Mu]3+2 \[Lambda]m (4 \[Mu]3-5 \[Mu]m)+3 (9 g2^2+6 yt1^2+3 g1^2 Y\[Phi]^2+8 \[Lambda]1H) \[Mu]m)+48 (m1 \[Mu]m+\[Mu]3 \[Mu]\[Sigma]))+2 T^2 (3 Lf yt1^2 \[Mu]m+g1^2 Y\[Phi]^2 \[Mu]m+6 EulerGamma (-4 \[Lambda]\[Sigma] \[Mu]3+(3 g2^2+g1^2 Y\[Phi]^2-2 \[Lambda]m) \[Mu]m)+3 g2^2 \[Mu]m (1-72 Log[Glaisher])+288 \[Lambda]\[Sigma] \[Mu]3 Log[Glaisher]-72 g1^2 Y\[Phi]^2 \[Mu]m Log[Glaisher]+144 \[Lambda]m \[Mu]m Log[Glaisher])+24 Sqrt[T] (4 \[Lambda]\[Sigma]3d \[Mu]33d-(3 g23d^2+g13d^2 Y\[Phi]^2-2 \[Lambda]m3d) \[Mu]m3d) Log[\[Mu]3/\[Mu]])}//Expand
]];


AppendTo[testList,
TestCreate[PrintTemporalScalarCouplings[]/.Replacements,
	{\[Lambda]VLL[1]->-((g2^2 g3^2 nF T)/(4 \[Pi]^2)),\[Lambda]VLL[2]->(g3^4 (9-2 nF) T)/(4 \[Pi]^2),\[Lambda]VLL[3]->(g2^4 (17-4 nF) T)/(8 \[Pi]^2),\[Lambda]VLL[4]->-((I g1 g2^3 nF T (Yl+3 Yq))/(24 \[Pi]^2)),\[Lambda]VLL[5]->-((I g1 g3^3 nF T (Yd-2 Yq+Yu))/(48 \[Pi]^2)),\[Lambda]VLL[6]->-((g1 g3^3 nF T (Yd+2 Yq+Yu))/(16 \[Pi]^2)),\[Lambda]VLL[7]->-((g1^2 g3^2 nF T (Yd^2+2 Yq^2+Yu^2))/(8 \[Pi]^2)),\[Lambda]VLL[8]->(g1^2 g2^2 T (-nF (Yl^2+3 Yq^2)+Y\[Phi]^2))/(8 \[Pi]^2),\[Lambda]VLL[9]->-((g1^4 T (nF (3 Yd^4+Ye^4+2 Yl^4+6 Yq^4+3 Yu^4)-2 Y\[Phi]^4))/(16 \[Pi]^2)),\[Lambda]VVSL[1]->(g2^2 Sqrt[T] \[Mu]m)/(16 \[Pi]^2),\[Lambda]VVSL[2]->(g1^2 Sqrt[T] Y\[Phi]^2 \[Mu]m)/(16 \[Pi]^2),\[Lambda]VL[1]->-((g3^2 T yt1^2)/(4 \[Pi]^2)),\[Lambda]VL[2]->1/192 g1^2 T (96 Y\[Phi]^2+1/\[Pi]^2 (-36 Lf Yq yt1^2 Yu+18 yt1^2 ((-2+Lf) Yu^2-Lf Y\[Phi]^2)+6 Yq^2 (3 (-2+Lf) yt1^2-g1^2 (-1+Lf) nF Y\[Phi]^2)+Y\[Phi]^2 (9 g2^2+g1^2 (-((-1+Lf) nF (3 Yd^2+Ye^2+2 Yl^2+3 Yu^2))-(-1+Lb) Y\[Phi]^2)+72 \[Lambda]1H))),\[Lambda]VL[3]->(g2^2 T (g2^2 (51+43 Lb+8 nF-8 Lf nF)+96 \[Pi]^2-36 yt1^2+3 g1^2 Y\[Phi]^2+72 \[Lambda]1H))/(192 \[Pi]^2),\[Lambda]VL[4]->1/384 g1 g2 T (192 Y\[Phi]+1/\[Pi]^2 (-36 (-2+Lf) Yq yt1^2+36 Lf yt1^2 (Yu-Y\[Phi])+g2^2 (-12+43 Lb+8 nF) Y\[Phi]-6 g1^2 (-1+Lf) nF Yq^2 Y\[Phi]-Lf nF (8 g2^2+g1^2 (3 Yd^2+Ye^2+2 Yl^2+3 Yu^2)) Y\[Phi]+g1^2 Y\[Phi] (nF (3 Yd^2+Ye^2+2 Yl^2+3 Yu^2)-(-4+Lb) Y\[Phi]^2)+48 Y\[Phi] \[Lambda]1H)),\[Lambda]VL[5]->(g2^2 T \[Lambda]m)/(8 \[Pi]^2),\[Lambda]VL[6]->(g1^2 T Y\[Phi]^2 \[Lambda]m)/(8 \[Pi]^2)}
]];


AppendTo[testList,
TestCreate[PrintDebyeMass["LO"]/.Replacements,
	{\[Mu]sqSU2->1/6 g2^2 (5+2 nF) T^2,\[Mu]sqSU3->1/3 g3^2 (3+nF) T^2,\[Mu]sqU1->1/24 g1^2 T^2 (nF (3 Yd^2+Ye^2+2 Yl^2+6 Yq^2+3 Yu^2)+4 Y\[Phi]^2)}
]];
AppendTo[testList,
TestCreate[PrintDebyeMass["NLO"]/.Replacements,
	{\[Mu]sqSU2->1/(1152 \[Pi]^2) g2^2 (144 m1+T^2 (-144 g3^2 nF+9 g1^2 (-nF (Yl^2+3 Yq^2)+Y\[Phi]^2)+6 (-3 yt1^2+12 \[Lambda]1H+\[Lambda]m)+g2^2 (207-4 EulerGamma (5+2 nF) (-43+8 nF)+4 nF (11+8 nF-192 Log[2])))+4 g2^2 T^2 (8 nF (-7+2 nF) Log[\[Pi] T]-39 Log[4 \[Pi] T]+(39+8 (7-2 nF) nF) Log[\[Mu]]+2 (88-5 nF) Log[\[Mu]/(4 \[Pi] T)])),\[Mu]sqSU3->1/(1152 \[Pi]^2) g3^2 T^2 (-9 (6 g2^2 nF+8 yt1^2+g1^2 nF (Yd^2+2 Yq^2+Yu^2))-8 g3^2 (-45+2 EulerGamma (3+nF) (-33+4 nF)+nF (-3-4 nF+132 Log[2]))+16 g3^2 (nF (-21+4 nF) (Log[\[Pi]]+Log[T]-Log[\[Mu]])+99 Log[\[Mu]/(4 \[Pi] T)])),\[Mu]sqU1->1/(576 \[Pi]^2) (-27 g1^2 T^2 yt1^2 (Yq^2+Yu^2)-9/4 g1^2 nF T^2 (6 g2^2 (Yl^2+3 Yq^2)+16 g3^2 (Yd^2+2 Yq^2+Yu^2)+g1^2 (3 Yd^4+Ye^4+2 Yl^4+6 Yq^4+3 Yu^4))+72 g1^2 m1 Y\[Phi]^2+18 g1^2 T^2 yt1^2 Y\[Phi]^2+9/2 T^2 (3 g1^2 g2^2 Y\[Phi]^2+g1^4 Y\[Phi]^4)+3 g1^2 T^2 Y\[Phi]^2 (12 \[Lambda]1H+\[Lambda]m)-2 g1^4 T^2 Y\[Phi]^4 (1+EulerGamma-Log[4 \[Pi] T]+Log[\[Mu]])-1/4 g1^4 nF^2 T^2 (3 Yd^2+Ye^2+2 Yl^2+6 Yq^2+3 Yu^2)^2 (-1+2 EulerGamma-2 Log[\[Pi] T]+2 Log[\[Mu]])-1/2 g1^4 nF T^2 (3 Yd^2+Ye^2+2 Yl^2+6 Yq^2+3 Yu^2) Y\[Phi]^2 (-1+5 EulerGamma+Log[1/(4 \[Pi]^5)]+5 Log[\[Mu]/T]))}
]];


AppendTo[testList,
TestCreate[PrintScalarMass["LO"]/.Replacements//Simplify,
	{m13d->m1+1/48 T^2 (9 g2^2+12 yt1^2+3 g1^2 Y\[Phi]^2+24 \[Lambda]1H+2 \[Lambda]m),\[Mu]\[Sigma]3d->1/12 T^2 (2 \[Lambda]m+3 \[Lambda]\[Sigma])+\[Mu]\[Sigma]}
]];
AppendTo[testList,
TestCreate[PrintScalarMass["NLO"]/.Replacements//Simplify,
	{m13d->1/(9216 \[Pi]^2) (-1728 Lf m1 yt1^2-1152 g3^2 T^2 yt1^2-384 g3^2 Lb T^2 yt1^2+1536 g3^2 Lf T^2 yt1^2-36 g1^2 T^2 Yq^2 yt1^2+18 g1^2 Lb T^2 Yq^2 yt1^2+90 g1^2 Lf T^2 Yq^2 yt1^2+216 Lb T^2 yt1^4-144 g1^2 T^2 Yq yt1^2 Yu-108 g1^2 Lb T^2 Yq yt1^2 Yu+108 g1^2 Lf T^2 Yq yt1^2 Yu-36 g1^2 T^2 yt1^2 Yu^2+18 g1^2 Lb T^2 yt1^2 Yu^2+90 g1^2 Lf T^2 yt1^2 Yu^2+432 g1^2 Lb m1 Y\[Phi]^2+6 g1^4 nF T^2 Yd^2 Y\[Phi]^2-27 g1^4 Lb nF T^2 Yd^2 Y\[Phi]^2+9 g1^4 Lf nF T^2 Yd^2 Y\[Phi]^2+2 g1^4 nF T^2 Ye^2 Y\[Phi]^2-9 g1^4 Lb nF T^2 Ye^2 Y\[Phi]^2+3 g1^4 Lf nF T^2 Ye^2 Y\[Phi]^2+4 g1^4 nF T^2 Yl^2 Y\[Phi]^2-18 g1^4 Lb nF T^2 Yl^2 Y\[Phi]^2+6 g1^4 Lf nF T^2 Yl^2 Y\[Phi]^2+12 g1^4 nF T^2 Yq^2 Y\[Phi]^2-54 g1^4 Lb nF T^2 Yq^2 Y\[Phi]^2+18 g1^4 Lf nF T^2 Yq^2 Y\[Phi]^2+108 g1^2 Lb T^2 yt1^2 Y\[Phi]^2-108 g1^2 Lf T^2 yt1^2 Y\[Phi]^2+6 g1^4 nF T^2 Yu^2 Y\[Phi]^2-27 g1^4 Lb nF T^2 Yu^2 Y\[Phi]^2+9 g1^4 Lf nF T^2 Yu^2 Y\[Phi]^2+2 g1^4 T^2 Y\[Phi]^4-126 EulerGamma g1^4 T^2 Y\[Phi]^4+66 g1^4 Lb T^2 Y\[Phi]^4-3456 Lb m1 \[Lambda]1H-1296 Lb T^2 yt1^2 \[Lambda]1H-432 Lf T^2 yt1^2 \[Lambda]1H+144 g1^2 T^2 Y\[Phi]^2 \[Lambda]1H+864 EulerGamma g1^2 T^2 Y\[Phi]^2 \[Lambda]1H-432 g1^2 Lb T^2 Y\[Phi]^2 \[Lambda]1H-3456 EulerGamma T^2 \[Lambda]1H^2-72 Lf T^2 yt1^2 \[Lambda]m+18 g1^2 Lb T^2 Y\[Phi]^2 \[Lambda]m-144 Lb T^2 \[Lambda]1H \[Lambda]m-144 EulerGamma T^2 \[Lambda]m^2+24 Lb T^2 \[Lambda]m^2-72 Lb T^2 \[Lambda]m \[Lambda]\[Sigma]-144 Lb \[Mu]m^2-288 Lb \[Lambda]m \[Mu]\[Sigma]+54 g2^2 (Lb (24 m1+T^2 (7 yt1^2+8 g1^2 Y\[Phi]^2-24 \[Lambda]1H+\[Lambda]m))-T^2 ((2+Lf) yt1^2-8 \[Lambda]1H (1+6 EulerGamma-72 Log[Glaisher])+2 g1^2 Y\[Phi]^2 (1+5 EulerGamma-60 Log[Glaisher])))+6 g2^4 T^2 (167+243 EulerGamma+8 nF+12 Lf nF-3 Lb (47+12 nF)-2916 Log[Glaisher])+1512 g1^4 T^2 Y\[Phi]^4 Log[Glaisher]-10368 g1^2 T^2 Y\[Phi]^2 \[Lambda]1H Log[Glaisher]+41472 T^2 \[Lambda]1H^2 Log[Glaisher]+1728 T^2 \[Lambda]m^2 Log[Glaisher]-36 Log[\[Mu]3/\[Mu]] (39 g23d^4-5 g13d^4 Y\[Phi]^4+48 g13d^2 Y\[Phi]^2 \[Lambda]1H3d+g23d^2 (-18 g13d^2 Y\[Phi]^2+144 \[Lambda]1H3d+96 \[Lambda]VL[3])-8 (24 \[Lambda]1H3d^2+\[Lambda]m3d^2-48 g33d^2 \[Lambda]VL[1]+8 \[Lambda]VL[1]^2+\[Lambda]VL[2]^2+3 \[Lambda]VL[3]^2+6 \[Lambda]VL[4]^2))),\[Mu]\[Sigma]3d->-(1/(384 \[Pi]^2))(Lb (48 m1 \[Lambda]m+27 g2^2 T^2 \[Lambda]m+18 T^2 yt1^2 \[Lambda]m+9 g1^2 T^2 Y\[Phi]^2 \[Lambda]m+24 T^2 \[Lambda]1H \[Lambda]m-10 T^2 \[Lambda]m^2+12 T^2 \[Lambda]m \[Lambda]\[Sigma]-18 T^2 \[Lambda]\[Sigma]^2+48 \[Mu]3^2+12 \[Mu]m^2+72 \[Lambda]\[Sigma] \[Mu]\[Sigma])-2 T^2 (3 Lf yt1^2 \[Lambda]m+g1^2 Y\[Phi]^2 \[Lambda]m+6 EulerGamma g1^2 Y\[Phi]^2 \[Lambda]m-12 EulerGamma \[Lambda]m^2-36 EulerGamma \[Lambda]\[Sigma]^2+3 g2^2 \[Lambda]m (1+6 EulerGamma-72 Log[Glaisher])-72 g1^2 Y\[Phi]^2 \[Lambda]m Log[Glaisher]+144 \[Lambda]m^2 Log[Glaisher]+432 \[Lambda]\[Sigma]^2 Log[Glaisher])+12 Log[\[Mu]3/\[Mu]] (2 g13d^2 Y\[Phi]^2 \[Lambda]m3d-4 \[Lambda]m3d^2-12 \[Lambda]\[Sigma]3d^2-3 \[Lambda]VL[5]^2+6 g23d^2 (\[Lambda]m3d+2 \[Lambda]VL[5])-\[Lambda]VL[6]^2))}
]];


(* ::Subsubsection::Closed:: *)
(*Dimension 6 matching relations*)


AppendTo[testList,
TestCreate[PrintCouplingsEffective[],
	{c613d->1/(6144 \[Pi]^4) (96 \[Pi]^2 T^2 (c61 (64 \[Pi]^2-36 Lf yt1^2+9 Lb (3 g2^2+g1^2 Y\[Phi]^2-24 \[Lambda]1H))-4 c63 Lb \[Lambda]m)+(9 g2^6+9 g1^2 g2^4 Y\[Phi]^2+9 g1^4 g2^2 Y\[Phi]^4+3 g1^6 Y\[Phi]^6+8 (-84 yt1^6+240 \[Lambda]1H^3+\[Lambda]m^3)) Zeta[3]),c623d->c62 T^2-(Lb T^2 (c64 \[Lambda]m+45 c62 \[Lambda]\[Sigma]))/(16 \[Pi]^2)+((\[Lambda]m^3+54 \[Lambda]\[Sigma]^3) Zeta[3])/(1536 \[Pi]^4),c633d->1/(256 \[Pi]^4) (-8 \[Pi]^2 T^2 (12 (c61+c64) Lb \[Lambda]m+c63 (-9 g2^2 Lb-32 \[Pi]^2+12 Lf yt1^2+Lb (-3 g1^2 Y\[Phi]^2+48 \[Lambda]1H+16 \[Lambda]m+6 \[Lambda]\[Sigma])))+\[Lambda]m (24 \[Lambda]1H^2+12 \[Lambda]1H \[Lambda]m+\[Lambda]m (2 \[Lambda]m+3 \[Lambda]\[Sigma])) Zeta[3]),c643d->1/(256 \[Pi]^4) (-4 \[Pi]^2 T^2 (12 (5 c62+c63) Lb \[Lambda]m+c64 (-64 \[Pi]^2+12 Lf yt1^2+Lb (-9 g2^2-3 g1^2 Y\[Phi]^2+24 \[Lambda]1H+32 \[Lambda]m+72 \[Lambda]\[Sigma])))+\[Lambda]m (3 \[Lambda]1H \[Lambda]m+(\[Lambda]m+3 \[Lambda]\[Sigma])^2) Zeta[3])}
]];


(* ::Subsubsection::Closed:: *)
(*Symmetric pressure*)


PrintPressure["NLO"]


AppendTo[testList,
TestCreate[PrintPressure["LO"],
	(29*Pi^2*T^4)/90 + (7*nF*Pi^2*T^4)/24
]];
AppendTo[testList,
TestCreate[PrintPressure["NLO"]/.Replacements//Simplify,
	-1/2304*(T^2*(384*m1 + T^2*(64*g3^2*(6 + 5*nF) + 12*g2^2*(13 + 10*nF) + 
     5*g1^2*(nF*(3*Yd^2 + Ye^2 + 2*Yl^2 + 6*Yq^2 + 3*Yu^2) + 4*Y\[Phi]^2) + 
     4*(30*yt1^2 + 24*\[Lambda]1H + 4*\[Lambda]m + 3*\[Lambda]\[Sigma])) + 96*\[Mu]\[Sigma]))
]];
AppendTo[testList,
TestCreate[PrintPressure["NNLO"],
	-1/4423680*(69120*Lb*(4*m1^2 + \[Mu]\[Sigma]^2) + 480*T^4*(3*g2^2 + g1^2*Y\[Phi]^2)*(12*\[Lambda]1H + \[Lambda]m)*
    (2 + 12*EulerGamma - 9*Lb - 144*Log[Glaisher]) + 2880*T^4*yt1^2*(12*\[Lambda]1H + \[Lambda]m)*
    (8*EulerGamma - 7*Lb + Lf - 96*Log[Glaisher]) + 
   69120*m1*T^2*yt1^2*(4*EulerGamma - 5*Lb + Lf - 48*Log[Glaisher]) + 
   23040*m1*T^2*(3*g2^2 + g1^2*Y\[Phi]^2)*(1 + 3*EulerGamma - 3*Lb - 36*Log[Glaisher]) + 
   5760*T^2*(18*g2^2*m1 - 2*m1*(12*yt1^2 - 3*g1^2*Y\[Phi]^2 + 24*\[Lambda]1H + 2*\[Lambda]m) - 4*\[Mu]3^2 - 
     3*\[Mu]m^2 - 4*\[Lambda]m*\[Mu]\[Sigma] - 6*\[Lambda]\[Sigma]*\[Mu]\[Sigma])*(2*EulerGamma - Lb - 24*Log[Glaisher]) + 
   23040*T^2*(24*m1*\[Lambda]1H + 2*m1*\[Lambda]m + 2*\[Lambda]m*\[Mu]\[Sigma] + 3*\[Lambda]\[Sigma]*\[Mu]\[Sigma])*
    (EulerGamma - Lb - 12*Log[Glaisher]) + 20*T^4*(10320*g2^4 + 26832*EulerGamma*g2^4 + 
     50688*g3^4 + 101376*EulerGamma*g3^4 - 13416*g2^4*Lb - 50688*g3^4*Lb + 
     3240*g2^4*nF + 15648*EulerGamma*g2^4*nF + 14976*g3^4*nF + 
     72192*EulerGamma*g3^4*nF - 10920*g2^4*Lb*nF - 48768*g3^4*Lb*nF + 3096*g2^4*Lf*nF + 
     12672*g3^4*Lf*nF - 960*g2^4*nF^2 - 3840*EulerGamma*g2^4*nF^2 - 2560*g3^4*nF^2 - 
     10240*EulerGamma*g3^4*nF^2 + 2496*g2^4*Lb*nF^2 + 6656*g3^4*Lb*nF^2 - 
     576*g2^4*Lf*nF^2 - 1536*g3^4*Lf*nF^2 - 45*g1^4*nF^2*Yd^4 - 
     180*EulerGamma*g1^4*nF^2*Yd^4 + 117*g1^4*Lb*nF^2*Yd^4 - 27*g1^4*Lf*nF^2*Yd^4 - 
     30*g1^4*nF^2*Yd^2*Ye^2 - 120*EulerGamma*g1^4*nF^2*Yd^2*Ye^2 + 
     78*g1^4*Lb*nF^2*Yd^2*Ye^2 - 18*g1^4*Lf*nF^2*Yd^2*Ye^2 - 5*g1^4*nF^2*Ye^4 - 
     20*EulerGamma*g1^4*nF^2*Ye^4 + 13*g1^4*Lb*nF^2*Ye^4 - 3*g1^4*Lf*nF^2*Ye^4 - 
     60*g1^4*nF^2*Yd^2*Yl^2 - 240*EulerGamma*g1^4*nF^2*Yd^2*Yl^2 + 
     156*g1^4*Lb*nF^2*Yd^2*Yl^2 - 36*g1^4*Lf*nF^2*Yd^2*Yl^2 - 20*g1^4*nF^2*Ye^2*Yl^2 - 
     80*EulerGamma*g1^4*nF^2*Ye^2*Yl^2 + 52*g1^4*Lb*nF^2*Ye^2*Yl^2 - 
     12*g1^4*Lf*nF^2*Ye^2*Yl^2 - 20*g1^4*nF^2*Yl^4 - 80*EulerGamma*g1^4*nF^2*Yl^4 + 
     52*g1^4*Lb*nF^2*Yl^4 - 12*g1^4*Lf*nF^2*Yl^4 - 180*g1^4*nF^2*Yd^2*Yq^2 - 
     720*EulerGamma*g1^4*nF^2*Yd^2*Yq^2 + 468*g1^4*Lb*nF^2*Yd^2*Yq^2 - 
     108*g1^4*Lf*nF^2*Yd^2*Yq^2 - 60*g1^4*nF^2*Ye^2*Yq^2 - 
     240*EulerGamma*g1^4*nF^2*Ye^2*Yq^2 + 156*g1^4*Lb*nF^2*Ye^2*Yq^2 - 
     36*g1^4*Lf*nF^2*Ye^2*Yq^2 - 120*g1^4*nF^2*Yl^2*Yq^2 - 
     480*EulerGamma*g1^4*nF^2*Yl^2*Yq^2 + 312*g1^4*Lb*nF^2*Yl^2*Yq^2 - 
     72*g1^4*Lf*nF^2*Yl^2*Yq^2 - 180*g1^4*nF^2*Yq^4 - 720*EulerGamma*g1^4*nF^2*Yq^4 + 
     468*g1^4*Lb*nF^2*Yq^4 - 108*g1^4*Lf*nF^2*Yq^4 + 6480*EulerGamma*g2^2*yt1^2 + 
     23040*EulerGamma*g3^2*yt1^2 - 4212*g2^2*Lb*yt1^2 - 14976*g3^2*Lb*yt1^2 + 
     972*g2^2*Lf*yt1^2 + 3456*g3^2*Lf*yt1^2 + 288*g1^2*Yq^2*yt1^2 + 
     1440*EulerGamma*g1^2*Yq^2*yt1^2 - 720*g1^2*Lb*Yq^2*yt1^2 - 
     12960*EulerGamma*yt1^4 + 8424*Lb*yt1^4 - 1944*Lf*yt1^4 - 576*g1^2*Yq*yt1^2*Yu + 
     1440*EulerGamma*g1^2*Yq*yt1^2*Yu - 1368*g1^2*Lb*Yq*yt1^2*Yu + 
     648*g1^2*Lf*Yq*yt1^2*Yu - 90*g1^4*nF^2*Yd^2*Yu^2 - 360*EulerGamma*g1^4*nF^2*Yd^2*
      Yu^2 + 234*g1^4*Lb*nF^2*Yd^2*Yu^2 - 54*g1^4*Lf*nF^2*Yd^2*Yu^2 - 
     30*g1^4*nF^2*Ye^2*Yu^2 - 120*EulerGamma*g1^4*nF^2*Ye^2*Yu^2 + 
     78*g1^4*Lb*nF^2*Ye^2*Yu^2 - 18*g1^4*Lf*nF^2*Ye^2*Yu^2 - 60*g1^4*nF^2*Yl^2*Yu^2 - 
     240*EulerGamma*g1^4*nF^2*Yl^2*Yu^2 + 156*g1^4*Lb*nF^2*Yl^2*Yu^2 - 
     36*g1^4*Lf*nF^2*Yl^2*Yu^2 - 180*g1^4*nF^2*Yq^2*Yu^2 - 
     720*EulerGamma*g1^4*nF^2*Yq^2*Yu^2 + 468*g1^4*Lb*nF^2*Yq^2*Yu^2 - 
     108*g1^4*Lf*nF^2*Yq^2*Yu^2 + 288*g1^2*yt1^2*Yu^2 + 1440*EulerGamma*g1^2*yt1^2*
      Yu^2 - 720*g1^2*Lb*yt1^2*Yu^2 - 45*g1^4*nF^2*Yu^4 - 180*EulerGamma*g1^4*nF^2*
      Yu^4 + 117*g1^4*Lb*nF^2*Yu^4 - 27*g1^4*Lf*nF^2*Yu^4 - 63*g1^4*nF*Yd^2*Y\[Phi]^2 - 
     300*EulerGamma*g1^4*nF*Yd^2*Y\[Phi]^2 + 159*g1^4*Lb*nF*Yd^2*Y\[Phi]^2 - 
     9*g1^4*Lf*nF*Yd^2*Y\[Phi]^2 - 21*g1^4*nF*Ye^2*Y\[Phi]^2 - 100*EulerGamma*g1^4*nF*Ye^2*
      Y\[Phi]^2 + 53*g1^4*Lb*nF*Ye^2*Y\[Phi]^2 - 3*g1^4*Lf*nF*Ye^2*Y\[Phi]^2 - 
     42*g1^4*nF*Yl^2*Y\[Phi]^2 - 200*EulerGamma*g1^4*nF*Yl^2*Y\[Phi]^2 + 
     106*g1^4*Lb*nF*Yl^2*Y\[Phi]^2 - 6*g1^4*Lf*nF*Yl^2*Y\[Phi]^2 - 126*g1^4*nF*Yq^2*Y\[Phi]^2 - 
     600*EulerGamma*g1^4*nF*Yq^2*Y\[Phi]^2 + 318*g1^4*Lb*nF*Yq^2*Y\[Phi]^2 - 
     18*g1^4*Lf*nF*Yq^2*Y\[Phi]^2 - 288*g1^2*yt1^2*Y\[Phi]^2 + 720*EulerGamma*g1^2*yt1^2*
      Y\[Phi]^2 - 684*g1^2*Lb*yt1^2*Y\[Phi]^2 + 324*g1^2*Lf*yt1^2*Y\[Phi]^2 - 
     63*g1^4*nF*Yu^2*Y\[Phi]^2 - 300*EulerGamma*g1^4*nF*Yu^2*Y\[Phi]^2 + 
     159*g1^4*Lb*nF*Yu^2*Y\[Phi]^2 - 9*g1^4*Lf*nF*Yu^2*Y\[Phi]^2 - 16*g1^4*Y\[Phi]^4 - 
     80*EulerGamma*g1^4*Y\[Phi]^4 + 40*g1^4*Lb*Y\[Phi]^4 - 306432*g2^4*Log[Glaisher] - 
     1216512*g3^4*Log[Glaisher] - 187776*g2^4*nF*Log[Glaisher] - 
     866304*g3^4*nF*Log[Glaisher] + 46080*g2^4*nF^2*Log[Glaisher] + 
     122880*g3^4*nF^2*Log[Glaisher] + 2160*g1^4*nF^2*Yd^4*Log[Glaisher] + 
     1440*g1^4*nF^2*Yd^2*Ye^2*Log[Glaisher] + 240*g1^4*nF^2*Ye^4*Log[Glaisher] + 
     2880*g1^4*nF^2*Yd^2*Yl^2*Log[Glaisher] + 960*g1^4*nF^2*Ye^2*Yl^2*Log[Glaisher] + 
     960*g1^4*nF^2*Yl^4*Log[Glaisher] + 8640*g1^4*nF^2*Yd^2*Yq^2*Log[Glaisher] + 
     2880*g1^4*nF^2*Ye^2*Yq^2*Log[Glaisher] + 5760*g1^4*nF^2*Yl^2*Yq^2*Log[Glaisher] + 
     8640*g1^4*nF^2*Yq^4*Log[Glaisher] - 77760*g2^2*yt1^2*Log[Glaisher] - 
     276480*g3^2*yt1^2*Log[Glaisher] - 17280*g1^2*Yq^2*yt1^2*Log[Glaisher] + 
     72576*yt1^4*Log[Glaisher] - 17280*g1^2*Yq*yt1^2*Yu*Log[Glaisher] + 
     4320*g1^4*nF^2*Yd^2*Yu^2*Log[Glaisher] + 1440*g1^4*nF^2*Ye^2*Yu^2*Log[Glaisher] + 
     2880*g1^4*nF^2*Yl^2*Yu^2*Log[Glaisher] + 8640*g1^4*nF^2*Yq^2*Yu^2*Log[Glaisher] - 
     17280*g1^2*yt1^2*Yu^2*Log[Glaisher] + 2160*g1^4*nF^2*Yu^4*Log[Glaisher] + 
     10368*g1^2*g2^2*Y\[Phi]^2*Log[Glaisher] + 3600*g1^4*nF*Yd^2*Y\[Phi]^2*Log[Glaisher] + 
     1200*g1^4*nF*Ye^2*Y\[Phi]^2*Log[Glaisher] + 2400*g1^4*nF*Yl^2*Y\[Phi]^2*Log[Glaisher] + 
     7200*g1^4*nF*Yq^2*Y\[Phi]^2*Log[Glaisher] - 8640*g1^2*yt1^2*Y\[Phi]^2*Log[Glaisher] + 
     3600*g1^4*nF*Yu^2*Y\[Phi]^2*Log[Glaisher] + 6144*g1^4*Y\[Phi]^4*Log[Glaisher] - 
     124416*g2^2*\[Lambda]1H*Log[Glaisher] + 165888*yt1^2*\[Lambda]1H*Log[Glaisher] - 
     41472*g1^2*Y\[Phi]^2*\[Lambda]1H*Log[Glaisher] + 331776*\[Lambda]1H^2*Log[Glaisher] - 
     10368*g2^2*\[Lambda]m*Log[Glaisher] + 13824*yt1^2*\[Lambda]m*Log[Glaisher] - 
     3456*g1^2*Y\[Phi]^2*\[Lambda]m*Log[Glaisher] + 27648*\[Lambda]1H*\[Lambda]m*Log[Glaisher] + 
     19584*\[Lambda]m^2*Log[Glaisher] + 13824*\[Lambda]m*\[Lambda]\[Sigma]*Log[Glaisher] + 
     31104*\[Lambda]\[Sigma]^2*Log[Glaisher] - 48*(27*g2^4 - 144*yt1^4 + 9*g1^4*Y\[Phi]^4 - 
       72*g1^2*Y\[Phi]^2*\[Lambda]1H + 576*\[Lambda]1H^2 + 18*g2^2*(g1^2*Y\[Phi]^2 - 12*\[Lambda]1H - \[Lambda]m) - 
       6*g1^2*Y\[Phi]^2*\[Lambda]m + 48*\[Lambda]1H*\[Lambda]m + 34*\[Lambda]m^2 + 24*yt1^2*(12*\[Lambda]1H + \[Lambda]m) + 
       24*\[Lambda]m*\[Lambda]\[Sigma] + 54*\[Lambda]\[Sigma]^2)*Log[4*Pi*T] + 24*(27*g2^4 - 144*yt1^4 + 9*g1^4*Y\[Phi]^4 - 
       72*g1^2*Y\[Phi]^2*\[Lambda]1H + 576*\[Lambda]1H^2 + 18*g2^2*(g1^2*Y\[Phi]^2 - 12*\[Lambda]1H - \[Lambda]m) - 
       6*g1^2*Y\[Phi]^2*\[Lambda]m + 48*\[Lambda]1H*\[Lambda]m + 34*\[Lambda]m^2 + 24*yt1^2*(12*\[Lambda]1H + \[Lambda]m) + 
       24*\[Lambda]m*\[Lambda]\[Sigma] + 54*\[Lambda]\[Sigma]^2)*Log[\[Mu]^2]) - 
   216*T^4*yt1^4*(211 + 760*EulerGamma - 100*Lb - 200*Lf - 1944*Log[2] - 720*Log[4] + 
     9600*Log[Glaisher] - 1440*Log[Pi] + 1440*Log[\[Mu]/T] - 
     14400*Derivative[1][Zeta][-3]) + 
   2*nF*T^4*(24*g2^4 + g1^4*(3*Yd^2 + Ye^2 + 2*Yl^2 + 6*Yq^2 + 3*Yu^2)*Y\[Phi]^2)*
    (997 - 1920*EulerGamma + 840*Lb - 120*Lf + 1188*Log[2] + 2940*Log[4] - 
     12480*Log[Glaisher] + 2940*Log[Pi] - 2940*Log[\[Mu]/T] + 
     2400*Derivative[1][Zeta][-3]) + 9216*(g2^4 + 6*g3^4)*T^4*
    (107 - 40*EulerGamma + 15*Lb + 310*Log[2] - 1760*Log[Glaisher] + 155*Log[Pi] - 
     155*Log[\[Mu]/T] + 3800*Derivative[1][Zeta][-3]) + 
   96*T^4*(1392*\[Lambda]1H^2 - 1440*EulerGamma*\[Lambda]1H^2 - 240*EulerGamma*\[Lambda]1H*\[Lambda]m + 116*\[Lambda]m^2 - 
     50*EulerGamma*\[Lambda]m^2 - 120*EulerGamma*\[Lambda]m*\[Lambda]\[Sigma] + 174*\[Lambda]\[Sigma]^2 - 
     90*EulerGamma*\[Lambda]\[Sigma]^2 + 17280*\[Lambda]1H^2*Log[2] + 1440*\[Lambda]1H*\[Lambda]m*Log[2] + 
     1020*\[Lambda]m^2*Log[2] + 720*\[Lambda]m*\[Lambda]\[Sigma]*Log[2] + 1620*\[Lambda]\[Sigma]^2*Log[2] - 
     103680*\[Lambda]1H^2*Log[Glaisher] - 5760*\[Lambda]1H*\[Lambda]m*Log[Glaisher] - 
     6960*\[Lambda]m^2*Log[Glaisher] - 2880*\[Lambda]m*\[Lambda]\[Sigma]*Log[Glaisher] - 
     10800*\[Lambda]\[Sigma]^2*Log[Glaisher] + 4320*\[Lambda]1H^2*Log[Pi] + 360*\[Lambda]m^2*Log[Pi] + 
     540*\[Lambda]\[Sigma]^2*Log[Pi] + 30*(144*\[Lambda]1H^2 + 24*\[Lambda]1H*\[Lambda]m + 5*\[Lambda]m^2 + 12*\[Lambda]m*\[Lambda]\[Sigma] + 
       9*\[Lambda]\[Sigma]^2)*Log[Pi*T] - 180*(24*\[Lambda]1H^2 + 2*\[Lambda]m^2 + 3*\[Lambda]\[Sigma]^2)*Log[\[Mu]/T] - 
     2160*\[Lambda]1H^2*Log[\[Mu]^2] - 360*\[Lambda]1H*\[Lambda]m*Log[\[Mu]^2] - 75*\[Lambda]m^2*Log[\[Mu]^2] - 
     180*\[Lambda]m*\[Lambda]\[Sigma]*Log[\[Mu]^2] - 135*\[Lambda]\[Sigma]^2*Log[\[Mu]^2] + 
     172800*\[Lambda]1H^2*Derivative[1][Zeta][-3] + 14400*\[Lambda]m^2*Derivative[1][Zeta][-3] + 
     21600*\[Lambda]\[Sigma]^2*Derivative[1][Zeta][-3]) - 
   24*T^4*yt1^2*(27136*g3^2 + 34560*EulerGamma*g3^2 - 13440*g3^2*Lb - 3840*g3^2*Lf + 
     1617*g1^2*Yq^2 + 3240*EulerGamma*g1^2*Yq^2 - 1260*g1^2*Lb*Yq^2 - 
     360*g1^2*Lf*Yq^2 + 1854*g1^2*Yq*Yu + 1617*g1^2*Yu^2 + 3240*EulerGamma*g1^2*Yu^2 - 
     1260*g1^2*Lb*Yu^2 - 360*g1^2*Lf*Yu^2 - 183*g1^2*Yq*Y\[Phi] + 183*g1^2*Yu*Y\[Phi] - 
     606*g1^2*Y\[Phi]^2 - 1440*EulerGamma*g1^2*Y\[Phi]^2 + 1260*g1^2*Lb*Y\[Phi]^2 - 
     180*g1^2*Lf*Y\[Phi]^2 + 48768*g3^2*Log[2] + 1512*g1^2*Yq^2*Log[2] + 
     6120*g1^2*Yq*Yu*Log[2] + 1512*g1^2*Yu^2*Log[2] - 612*g1^2*Yq*Y\[Phi]*Log[2] + 
     612*g1^2*Yu*Y\[Phi]*Log[2] - 1224*g1^2*Y\[Phi]^2*Log[2] - 11520*g3^2*Log[4] - 
     2160*g1^2*Yq^2*Log[4] + 2160*g1^2*Yq*Yu*Log[4] - 2160*g1^2*Yu^2*Log[4] - 
     540*g1^2*Yq*Y\[Phi]*Log[4] + 540*g1^2*Yu*Y\[Phi]*Log[4] - 1080*g1^2*Y\[Phi]^2*Log[4] - 
     460800*g3^2*Log[Glaisher] - 17280*g1^2*Yq^2*Log[Glaisher] - 
     51840*g1^2*Yq*Yu*Log[Glaisher] - 17280*g1^2*Yu^2*Log[Glaisher] + 
     8640*g1^2*Yq*Y\[Phi]*Log[Glaisher] - 8640*g1^2*Yu*Y\[Phi]*Log[Glaisher] + 
     34560*g1^2*Y\[Phi]^2*Log[Glaisher] - 5760*g3^2*Log[Pi] - 2160*g1^2*Yq^2*Log[Pi] + 
     3240*g1^2*Yq*Yu*Log[Pi] - 2160*g1^2*Yu^2*Log[Pi] - 540*g1^2*Yq*Y\[Phi]*Log[Pi] + 
     540*g1^2*Yu*Y\[Phi]*Log[Pi] - 1080*g1^2*Y\[Phi]^2*Log[Pi] + 
     180*(45*g2^2 + 32*g3^2 + 3*g1^2*(4*Yq^2 - 6*Yq*Yu + 4*Yu^2 + Yq*Y\[Phi] - Yu*Y\[Phi] + 
         2*Y\[Phi]^2))*Log[\[Mu]/T] + 1152000*g3^2*Derivative[1][Zeta][-3] + 
     43200*g1^2*Yq^2*Derivative[1][Zeta][-3] + 129600*g1^2*Yq*Yu*
      Derivative[1][Zeta][-3] + 43200*g1^2*Yu^2*Derivative[1][Zeta][-3] - 
     21600*g1^2*Yq*Y\[Phi]*Derivative[1][Zeta][-3] + 21600*g1^2*Yu*Y\[Phi]*
      Derivative[1][Zeta][-3] - 43200*g1^2*Y\[Phi]^2*Derivative[1][Zeta][-3] + 
     18*g2^2*(199 + 300*EulerGamma - 90*Lf + 150*Log[2] - 450*Log[4] + 
       1440*Log[Glaisher] - 450*Log[Pi] + 3600*Derivative[1][Zeta][-3])) - 
   8*T^4*(30*(2025*g2^4 + 594*g1^2*g2^2*Y\[Phi]^2 + 169*g1^4*Y\[Phi]^4)*Log[\[Mu]/T] + 
     3*g2^4*(-12589 + 7020*EulerGamma - 2025*Lb - 40500*Log[2] + 187920*Log[Glaisher] - 
       20250*Log[Pi] - 291600*Derivative[1][Zeta][-3]) + 
     g1^4*Y\[Phi]^4*(-2543 + 420*EulerGamma + 45*Lb - 10140*Log[2] + 70320*Log[Glaisher] - 
       5070*Log[Pi] - 145200*Derivative[1][Zeta][-3]) - 
     18*g1^2*g2^2*Y\[Phi]^2*(499 + 180*EulerGamma - 135*Lb + 1980*Log[2] - 
       18000*Log[Glaisher] + 990*Log[Pi] + 39600*Derivative[1][Zeta][-3])) + 
   nF*T^4*(544512*g3^4 - 2580480*EulerGamma*g3^4 + 1075200*g3^4*Lb + 30720*g3^4*Lf + 
     128000*g3^4*nF - 245760*EulerGamma*g3^4*nF + 122880*g3^4*Lb*nF - 
     30720*g3^4*Lf*nF - 50400*g1^2*g3^2*Yd^2 - 103680*EulerGamma*g1^2*g3^2*Yd^2 + 
     40320*g1^2*g3^2*Lb*Yd^2 + 11520*g1^2*g3^2*Lf*Yd^2 - 4725*g1^4*Yd^4 - 
     9720*EulerGamma*g1^4*Yd^4 + 3780*g1^4*Lb*Yd^4 + 1080*g1^4*Lf*Yd^4 + 
     2250*g1^4*nF*Yd^4 - 4320*EulerGamma*g1^4*nF*Yd^4 + 2160*g1^4*Lb*nF*Yd^4 - 
     540*g1^4*Lf*nF*Yd^4 + 1500*g1^4*nF*Yd^2*Ye^2 - 2880*EulerGamma*g1^4*nF*Yd^2*Ye^2 + 
     1440*g1^4*Lb*nF*Yd^2*Ye^2 - 360*g1^4*Lf*nF*Yd^2*Ye^2 - 1575*g1^4*Ye^4 - 
     3240*EulerGamma*g1^4*Ye^4 + 1260*g1^4*Lb*Ye^4 + 360*g1^4*Lf*Ye^4 + 
     250*g1^4*nF*Ye^4 - 480*EulerGamma*g1^4*nF*Ye^4 + 240*g1^4*Lb*nF*Ye^4 - 
     60*g1^4*Lf*nF*Ye^4 + 3000*g1^4*nF*Yd^2*Yl^2 - 5760*EulerGamma*g1^4*nF*Yd^2*Yl^2 + 
     2880*g1^4*Lb*nF*Yd^2*Yl^2 - 720*g1^4*Lf*nF*Yd^2*Yl^2 + 1000*g1^4*nF*Ye^2*Yl^2 - 
     1920*EulerGamma*g1^4*nF*Ye^2*Yl^2 + 960*g1^4*Lb*nF*Ye^2*Yl^2 - 
     240*g1^4*Lf*nF*Ye^2*Yl^2 - 3150*g1^4*Yl^4 - 6480*EulerGamma*g1^4*Yl^4 + 
     2520*g1^4*Lb*Yl^4 + 720*g1^4*Lf*Yl^4 + 1000*g1^4*nF*Yl^4 - 
     1920*EulerGamma*g1^4*nF*Yl^4 + 960*g1^4*Lb*nF*Yl^4 - 240*g1^4*Lf*nF*Yl^4 - 
     100800*g1^2*g3^2*Yq^2 - 207360*EulerGamma*g1^2*g3^2*Yq^2 + 
     80640*g1^2*g3^2*Lb*Yq^2 + 23040*g1^2*g3^2*Lf*Yq^2 + 9000*g1^4*nF*Yd^2*Yq^2 - 
     17280*EulerGamma*g1^4*nF*Yd^2*Yq^2 + 8640*g1^4*Lb*nF*Yd^2*Yq^2 - 
     2160*g1^4*Lf*nF*Yd^2*Yq^2 + 3000*g1^4*nF*Ye^2*Yq^2 - 5760*EulerGamma*g1^4*nF*Ye^2*
      Yq^2 + 2880*g1^4*Lb*nF*Ye^2*Yq^2 - 720*g1^4*Lf*nF*Ye^2*Yq^2 + 
     6000*g1^4*nF*Yl^2*Yq^2 - 11520*EulerGamma*g1^4*nF*Yl^2*Yq^2 + 
     5760*g1^4*Lb*nF*Yl^2*Yq^2 - 1440*g1^4*Lf*nF*Yl^2*Yq^2 - 9450*g1^4*Yq^4 - 
     19440*EulerGamma*g1^4*Yq^4 + 7560*g1^4*Lb*Yq^4 + 2160*g1^4*Lf*Yq^4 + 
     9000*g1^4*nF*Yq^4 - 17280*EulerGamma*g1^4*nF*Yq^4 + 8640*g1^4*Lb*nF*Yq^4 - 
     2160*g1^4*Lf*nF*Yq^4 - 50400*g1^2*g3^2*Yu^2 - 103680*EulerGamma*g1^2*g3^2*Yu^2 + 
     40320*g1^2*g3^2*Lb*Yu^2 + 11520*g1^2*g3^2*Lf*Yu^2 + 4500*g1^4*nF*Yd^2*Yu^2 - 
     8640*EulerGamma*g1^4*nF*Yd^2*Yu^2 + 4320*g1^4*Lb*nF*Yd^2*Yu^2 - 
     1080*g1^4*Lf*nF*Yd^2*Yu^2 + 1500*g1^4*nF*Ye^2*Yu^2 - 2880*EulerGamma*g1^4*nF*Ye^2*
      Yu^2 + 1440*g1^4*Lb*nF*Ye^2*Yu^2 - 360*g1^4*Lf*nF*Ye^2*Yu^2 + 
     3000*g1^4*nF*Yl^2*Yu^2 - 5760*EulerGamma*g1^4*nF*Yl^2*Yu^2 + 
     2880*g1^4*Lb*nF*Yl^2*Yu^2 - 720*g1^4*Lf*nF*Yl^2*Yu^2 + 9000*g1^4*nF*Yq^2*Yu^2 - 
     17280*EulerGamma*g1^4*nF*Yq^2*Yu^2 + 8640*g1^4*Lb*nF*Yq^2*Yu^2 - 
     2160*g1^4*Lf*nF*Yq^2*Yu^2 - 4725*g1^4*Yu^4 - 9720*EulerGamma*g1^4*Yu^4 + 
     3780*g1^4*Lb*Yu^4 + 1080*g1^4*Lf*Yu^4 + 2250*g1^4*nF*Yu^4 - 
     4320*EulerGamma*g1^4*nF*Yu^4 + 2160*g1^4*Lb*nF*Yu^4 - 540*g1^4*Lf*nF*Yu^4 + 
     709632*g3^4*Log[2] + 516096*g3^4*nF*Log[2] + 69120*g1^2*g3^2*Yd^2*Log[2] + 
     6480*g1^4*Yd^4*Log[2] + 9072*g1^4*nF*Yd^4*Log[2] + 6048*g1^4*nF*Yd^2*Ye^2*Log[2] + 
     2160*g1^4*Ye^4*Log[2] + 1008*g1^4*nF*Ye^4*Log[2] + 12096*g1^4*nF*Yd^2*Yl^2*
      Log[2] + 4032*g1^4*nF*Ye^2*Yl^2*Log[2] + 4320*g1^4*Yl^4*Log[2] + 
     4032*g1^4*nF*Yl^4*Log[2] + 138240*g1^2*g3^2*Yq^2*Log[2] + 
     36288*g1^4*nF*Yd^2*Yq^2*Log[2] + 12096*g1^4*nF*Ye^2*Yq^2*Log[2] + 
     24192*g1^4*nF*Yl^2*Yq^2*Log[2] + 12960*g1^4*Yq^4*Log[2] + 
     36288*g1^4*nF*Yq^4*Log[2] + 69120*g1^2*g3^2*Yu^2*Log[2] + 
     18144*g1^4*nF*Yd^2*Yu^2*Log[2] + 6048*g1^4*nF*Ye^2*Yu^2*Log[2] + 
     12096*g1^4*nF*Yl^2*Yu^2*Log[2] + 36288*g1^4*nF*Yq^2*Yu^2*Log[2] + 
     6480*g1^4*Yu^4*Log[2] + 9072*g1^4*nF*Yu^4*Log[2] + 3409920*g3^4*Log[4] + 
     491520*g3^4*nF*Log[4] + 69120*g1^2*g3^2*Yd^2*Log[4] + 6480*g1^4*Yd^4*Log[4] + 
     8640*g1^4*nF*Yd^4*Log[4] + 5760*g1^4*nF*Yd^2*Ye^2*Log[4] + 2160*g1^4*Ye^4*Log[4] + 
     960*g1^4*nF*Ye^4*Log[4] + 11520*g1^4*nF*Yd^2*Yl^2*Log[4] + 
     3840*g1^4*nF*Ye^2*Yl^2*Log[4] + 4320*g1^4*Yl^4*Log[4] + 3840*g1^4*nF*Yl^4*Log[4] + 
     138240*g1^2*g3^2*Yq^2*Log[4] + 34560*g1^4*nF*Yd^2*Yq^2*Log[4] + 
     11520*g1^4*nF*Ye^2*Yq^2*Log[4] + 23040*g1^4*nF*Yl^2*Yq^2*Log[4] + 
     12960*g1^4*Yq^4*Log[4] + 34560*g1^4*nF*Yq^4*Log[4] + 69120*g1^2*g3^2*Yu^2*Log[4] + 
     17280*g1^4*nF*Yd^2*Yu^2*Log[4] + 5760*g1^4*nF*Ye^2*Yu^2*Log[4] + 
     11520*g1^4*nF*Yl^2*Yu^2*Log[4] + 34560*g1^4*nF*Yq^2*Yu^2*Log[4] + 
     6480*g1^4*Yu^4*Log[4] + 8640*g1^4*nF*Yu^4*Log[4] - 9584640*g3^4*Log[Glaisher] - 
     3440640*g3^4*nF*Log[Glaisher] - 60480*g1^4*nF*Yd^4*Log[Glaisher] - 
     40320*g1^4*nF*Yd^2*Ye^2*Log[Glaisher] - 6720*g1^4*nF*Ye^4*Log[Glaisher] - 
     80640*g1^4*nF*Yd^2*Yl^2*Log[Glaisher] - 26880*g1^4*nF*Ye^2*Yl^2*Log[Glaisher] - 
     26880*g1^4*nF*Yl^4*Log[Glaisher] - 241920*g1^4*nF*Yd^2*Yq^2*Log[Glaisher] - 
     80640*g1^4*nF*Ye^2*Yq^2*Log[Glaisher] - 161280*g1^4*nF*Yl^2*Yq^2*Log[Glaisher] - 
     241920*g1^4*nF*Yq^4*Log[Glaisher] - 120960*g1^4*nF*Yd^2*Yu^2*Log[Glaisher] - 
     40320*g1^4*nF*Ye^2*Yu^2*Log[Glaisher] - 80640*g1^4*nF*Yl^2*Yu^2*Log[Glaisher] - 
     241920*g1^4*nF*Yq^2*Yu^2*Log[Glaisher] - 60480*g1^4*nF*Yu^4*Log[Glaisher] - 
     540*g2^2*(16*g3^2 + g1^2*(Yl^2 + 3*Yq^2))*(35 + 72*EulerGamma - 28*Lb - 8*Lf - 
       48*Log[2] - 48*Log[4] - 72*Log[Pi]) + 3363840*g3^4*Log[Pi] + 
     491520*g3^4*nF*Log[Pi] + 103680*g1^2*g3^2*Yd^2*Log[Pi] + 9720*g1^4*Yd^4*Log[Pi] + 
     8640*g1^4*nF*Yd^4*Log[Pi] + 5760*g1^4*nF*Yd^2*Ye^2*Log[Pi] + 
     3240*g1^4*Ye^4*Log[Pi] + 960*g1^4*nF*Ye^4*Log[Pi] + 11520*g1^4*nF*Yd^2*Yl^2*
      Log[Pi] + 3840*g1^4*nF*Ye^2*Yl^2*Log[Pi] + 6480*g1^4*Yl^4*Log[Pi] + 
     3840*g1^4*nF*Yl^4*Log[Pi] + 207360*g1^2*g3^2*Yq^2*Log[Pi] + 
     34560*g1^4*nF*Yd^2*Yq^2*Log[Pi] + 11520*g1^4*nF*Ye^2*Yq^2*Log[Pi] + 
     23040*g1^4*nF*Yl^2*Yq^2*Log[Pi] + 19440*g1^4*Yq^4*Log[Pi] + 
     34560*g1^4*nF*Yq^4*Log[Pi] + 103680*g1^2*g3^2*Yu^2*Log[Pi] + 
     17280*g1^4*nF*Yd^2*Yu^2*Log[Pi] + 5760*g1^4*nF*Ye^2*Yu^2*Log[Pi] + 
     11520*g1^4*nF*Yl^2*Yu^2*Log[Pi] + 34560*g1^4*nF*Yq^2*Yu^2*Log[Pi] + 
     9720*g1^4*Yu^4*Log[Pi] + 8640*g1^4*nF*Yu^4*Log[Pi] - 
     120*(128*g3^4*(219 + 32*nF) + 24*g2^4*(277 + 64*nF) + 
       324*g2^2*(16*g3^2 + g1^2*(Yl^2 + 3*Yq^2)) + 864*g1^2*g3^2*
        (Yd^2 + 2*Yq^2 + Yu^2) + g1^4*(9*(9 + 8*nF)*Yd^4 + (27 + 8*nF)*Ye^4 + 54*Yl^4 + 
         32*nF*Yl^4 + 192*nF*Yl^2*Yq^2 + 162*Yq^4 + 288*nF*Yq^4 + 96*nF*Yl^2*Yu^2 + 
         288*nF*Yq^2*Yu^2 + 81*Yu^4 + 72*nF*Yu^4 + 16*nF*Ye^2*(2*Yl^2 + 6*Yq^2 + 
           3*Yu^2) + 48*nF*Yd^2*(Ye^2 + 2*Yl^2 + 6*Yq^2 + 3*Yu^2)))*Log[\[Mu]/T] + 
     1843200*g3^4*Derivative[1][Zeta][-3] + 4915200*g3^4*nF*Derivative[1][Zeta][-3] + 
     86400*g1^4*nF*Yd^4*Derivative[1][Zeta][-3] + 57600*g1^4*nF*Yd^2*Ye^2*
      Derivative[1][Zeta][-3] + 9600*g1^4*nF*Ye^4*Derivative[1][Zeta][-3] + 
     115200*g1^4*nF*Yd^2*Yl^2*Derivative[1][Zeta][-3] + 38400*g1^4*nF*Ye^2*Yl^2*
      Derivative[1][Zeta][-3] + 38400*g1^4*nF*Yl^4*Derivative[1][Zeta][-3] + 
     345600*g1^4*nF*Yd^2*Yq^2*Derivative[1][Zeta][-3] + 115200*g1^4*nF*Ye^2*Yq^2*
      Derivative[1][Zeta][-3] + 230400*g1^4*nF*Yl^2*Yq^2*Derivative[1][Zeta][-3] + 
     345600*g1^4*nF*Yq^4*Derivative[1][Zeta][-3] + 172800*g1^4*nF*Yd^2*Yu^2*
      Derivative[1][Zeta][-3] + 57600*g1^4*nF*Ye^2*Yu^2*Derivative[1][Zeta][-3] + 
     115200*g1^4*nF*Yl^2*Yu^2*Derivative[1][Zeta][-3] + 345600*g1^4*nF*Yq^2*Yu^2*
      Derivative[1][Zeta][-3] + 86400*g1^4*nF*Yu^4*Derivative[1][Zeta][-3] + 
     24*g2^4*(6547 + 120*Lf + 2000*nF - 480*Lf*nF + 60*Lb*(175 + 32*nF) - 
       120*EulerGamma*(209 + 32*nF) + 6192*Log[2] + 8064*nF*Log[2] + 34320*Log[4] + 
       7680*nF*Log[4] - 99840*Log[Glaisher] - 53760*nF*Log[Glaisher] + 33240*Log[Pi] + 
       7680*nF*Log[Pi] + 19200*Derivative[1][Zeta][-3] + 
       76800*nF*Derivative[1][Zeta][-3])))/Pi^2
]];


(* ::Subsubsection:: *)
(*Report*)


report=TestReport[testList]
report["ResultsDataset"]


(* ::Subsection:: *)
(*Test scalar 1 and 2 softer*)


PerformDRsoft[{}]


testList={};


AppendTo[testList,
TestCreate[PrintCouplingsUS[],
	{\[Lambda]1H3dUS->\[Lambda]1H3d-((8 \[Lambda]VL[1]^2)/Sqrt[\[Mu]sqSU3]+\[Lambda]VL[2]^2/Sqrt[\[Mu]sqU1]+(3 \[Lambda]VL[3]^2)/Sqrt[\[Mu]sqSU2]+(4 \[Lambda]VL[4]^2)/(Sqrt[\[Mu]sqSU2]+Sqrt[\[Mu]sqU1]))/(32 \[Pi]),\[Lambda]m3dUS->-((-192 \[Pi] \[Lambda]m3d+(36 \[Lambda]VL[3] \[Lambda]VL[5])/Sqrt[\[Mu]sqSU2]+(3 \[Lambda]m3d \[Lambda]VVSL[1]^2)/\[Mu]sqSU2^(3/2)+(12 \[Mu]sqU1 \[Lambda]VL[2] \[Lambda]VL[6]+\[Lambda]m3d \[Lambda]VVSL[2]^2)/\[Mu]sqU1^(3/2))/(192 \[Pi])),\[Lambda]\[Sigma]3dUS->-((-96 \[Pi] \[Lambda]\[Sigma]3d+(9 \[Lambda]VL[5]^2)/Sqrt[\[Mu]sqSU2]+(3 \[Lambda]\[Sigma]3d \[Lambda]VVSL[1]^2)/\[Mu]sqSU2^(3/2)+(3 \[Mu]sqU1 \[Lambda]VL[6]^2+\[Lambda]\[Sigma]3d \[Lambda]VVSL[2]^2)/\[Mu]sqU1^(3/2))/(96 \[Pi])),g13dUS^2->g13d^2,g23dUS^2->g23d^2-g23d^4/(24 \[Pi] Sqrt[\[Mu]sqSU2]),g33dUS^2->g33d^2-g33d^4/(16 \[Pi] Sqrt[\[Mu]sqSU3]),\[Mu]33dUS->(128 \[Pi] \[Mu]33d-(36 \[Lambda]VL[5] \[Lambda]VVSL[1])/Sqrt[\[Mu]sqSU2]+(3 \[Mu]33d \[Lambda]VVSL[1]^2)/\[Mu]sqSU2^(3/2)+(\[Lambda]VVSL[2] (-12 \[Mu]sqU1 \[Lambda]VL[6]+\[Mu]33d \[Lambda]VVSL[2]))/\[Mu]sqU1^(3/2))/(128 \[Pi]),\[Mu]m3dUS->(384 \[Pi] \[Mu]m3d+(3 \[Lambda]VVSL[1] (-48 \[Mu]sqSU2 \[Lambda]VL[3]+\[Mu]m3d \[Lambda]VVSL[1]))/\[Mu]sqSU2^(3/2)-(48 \[Lambda]VL[2] \[Lambda]VVSL[2])/Sqrt[\[Mu]sqU1]+(\[Mu]m3d \[Lambda]VVSL[2]^2)/\[Mu]sqU1^(3/2))/(384 \[Pi])}
]];


AppendTo[testList,
TestCreate[PrintScalarMassUS["LO"],
	{m13dUS->m13d-(8 Sqrt[\[Mu]sqSU3] \[Lambda]VL[1]+Sqrt[\[Mu]sqU1] \[Lambda]VL[2]+3 Sqrt[\[Mu]sqSU2] \[Lambda]VL[3])/(8 \[Pi]),\[Mu]\[Sigma]3dUS->\[Mu]\[Sigma]3d-(6 Sqrt[\[Mu]sqSU2] \[Lambda]VL[5]+(3 \[Lambda]VVSL[1]^2)/Sqrt[\[Mu]sqSU2]+(2 \[Mu]sqU1 \[Lambda]VL[6]+\[Lambda]VVSL[2]^2)/Sqrt[\[Mu]sqU1])/(16 \[Pi])}
]];
AppendTo[testList,
TestCreate[PrintScalarMassUS["NLO"],
	{m13dUS->1/(128 \[Pi]^2) (48 g33d^2 \[Lambda]VL[1]+32 Log[\[Mu]3/(2 Sqrt[\[Mu]sqSU3])] (6 g33d^2-\[Lambda]VL[1]) \[Lambda]VL[1]-16 \[Lambda]VL[1]^2-2 \[Lambda]VL[2]^2-4 Log[\[Mu]3/(2 Sqrt[\[Mu]sqU1])] \[Lambda]VL[2]^2+12 g23d^2 \[Lambda]VL[3]-6 \[Lambda]VL[3]^2-6 Log[\[Mu]3/(2 Sqrt[\[Mu]sqSU2])] (g23d^4-8 g23d^2 \[Lambda]VL[3]+2 \[Lambda]VL[3]^2)-12 \[Lambda]VL[4]^2-24 Log[\[Mu]3/(Sqrt[\[Mu]sqSU2]+Sqrt[\[Mu]sqU1])] \[Lambda]VL[4]^2+(24 Sqrt[\[Mu]sqSU2] \[Lambda]VL[1] \[Lambda]VLL[1])/Sqrt[\[Mu]sqSU3]+(24 Sqrt[\[Mu]sqSU3] \[Lambda]VL[3] \[Lambda]VLL[1])/Sqrt[\[Mu]sqSU2]+80/3 \[Lambda]VL[1] \[Lambda]VLL[2]+5 \[Lambda]VL[3] \[Lambda]VLL[3]+(8 Sqrt[\[Mu]sqU1] \[Lambda]VL[1] \[Lambda]VLL[7])/Sqrt[\[Mu]sqSU3]+(8 Sqrt[\[Mu]sqSU3] \[Lambda]VL[2] \[Lambda]VLL[7])/Sqrt[\[Mu]sqU1]+(3 Sqrt[\[Mu]sqSU2] \[Lambda]VL[2] \[Lambda]VLL[8])/Sqrt[\[Mu]sqU1]+(3 Sqrt[\[Mu]sqU1] \[Lambda]VL[3] \[Lambda]VLL[8])/Sqrt[\[Mu]sqSU2]+\[Lambda]VL[2] \[Lambda]VLL[9]),\[Mu]\[Sigma]3dUS->1/(384 \[Pi]^2) (36 g23d^2 (1+4 Log[\[Mu]3/(2 Sqrt[\[Mu]sqSU2])]) \[Lambda]VL[5]-18 (1+2 Log[\[Mu]3/(2 Sqrt[\[Mu]sqSU2])]) \[Lambda]VL[5]^2-6 (1+2 Log[\[Mu]3/(2 Sqrt[\[Mu]sqU1])]) \[Lambda]VL[6]^2+(3 \[Lambda]VL[5] (24 Sqrt[\[Mu]sqSU3] \[Lambda]VLL[1]+5 Sqrt[\[Mu]sqSU2] \[Lambda]VLL[3]+3 Sqrt[\[Mu]sqU1] \[Lambda]VLL[8]))/Sqrt[\[Mu]sqSU2]+(3 \[Lambda]VL[6] (8 Sqrt[\[Mu]sqSU3] \[Lambda]VLL[7]+3 Sqrt[\[Mu]sqSU2] \[Lambda]VLL[8]+Sqrt[\[Mu]sqU1] \[Lambda]VLL[9]))/Sqrt[\[Mu]sqU1]+1/8 ((3 \[Lambda]VVSL[1]^2)/\[Mu]sqSU2^(3/2)+\[Lambda]VVSL[2]^2/\[Mu]sqU1^(3/2)) (16 \[Pi] \[Mu]\[Sigma]3d-6 Sqrt[\[Mu]sqSU2] \[Lambda]VL[5]-(3 \[Lambda]VVSL[1]^2)/Sqrt[\[Mu]sqSU2]+(-2 \[Mu]sqU1 \[Lambda]VL[6]-\[Lambda]VVSL[2]^2)/Sqrt[\[Mu]sqU1]))}
]];


AppendTo[testList,
TestCreate[BetaFunctions3DUS[],
	{m13dUS->1/(256 \[Pi]^2) (-51 g23dUS^4+5 g13dUS^4 Y\[Phi]^4+18 g23dUS^2 (g13dUS^2 Y\[Phi]^2-8 \[Lambda]1H3dUS)-48 g13dUS^2 Y\[Phi]^2 \[Lambda]1H3dUS+8 (24 \[Lambda]1H3dUS^2+\[Lambda]m3dUS^2)),\[Mu]\[Sigma]3dUS->(-((3 g23dUS^2+g13dUS^2 Y\[Phi]^2-2 \[Lambda]m3dUS) \[Lambda]m3dUS)+6 \[Lambda]\[Sigma]3dUS^2)/(16 \[Pi]^2),\[Mu]13dUS->(4 \[Lambda]\[Sigma]3dUS \[Mu]33dUS-3 g23dUS^2 \[Mu]m3dUS-g13dUS^2 Y\[Phi]^2 \[Mu]m3dUS+2 \[Lambda]m3dUS \[Mu]m3dUS)/(32 \[Pi]^2)}
]];


AppendTo[testList,
TestCreate[PrintTadpolesUS["LO"],
	{\[Mu]13dUS->\[Mu]13d-(3 Sqrt[\[Mu]sqSU2] \[Lambda]VVSL[1]+Sqrt[\[Mu]sqU1] \[Lambda]VVSL[2])/(8 \[Pi])-(\[Mu]13d ((3 \[Lambda]VVSL[1]^2)/\[Mu]sqSU2^(3/2)+\[Lambda]VVSL[2]^2/\[Mu]sqU1^(3/2)))/(384 \[Pi])}
]];


report=TestReport[testList]
report["ResultsDataset"]


(* ::Subsection:: *)
(*Test scalar 2 softer*)


PerformDRsoft[{5}];


testList={};


AppendTo[testList,
TestCreate[PrintCouplingsUS[],
	{\[Lambda]1H3dUS->1/768 (8 \[Lambda]1H3d (96-(38 \[Mu]m3d^2)/(\[Pi] \[Mu]\[Sigma]3d^(3/2)))-(96 \[Mu]m3d^2)/\[Mu]\[Sigma]3d+(192 \[Mu]13d \[Mu]m3d (\[Mu]33d \[Mu]m3d+2 \[Lambda]m3d \[Mu]\[Sigma]3d))/\[Mu]\[Sigma]3d^3-(3 \[Mu]m3d (24 \[Mu]33d^2 \[Mu]m3d+9 \[Mu]m3d^3+32 \[Lambda]m3d \[Mu]33d \[Mu]\[Sigma]3d))/(\[Pi] \[Mu]\[Sigma]3d^(5/2))-(24 ((-((\[Lambda]m3d-3 \[Lambda]\[Sigma]3d) \[Mu]m3d^2)+\[Lambda]m3d^2 \[Mu]\[Sigma]3d)/\[Mu]\[Sigma]3d^(3/2)+(8 \[Lambda]VL[1]^2)/Sqrt[\[Mu]sqSU3]+\[Lambda]VL[2]^2/Sqrt[\[Mu]sqU1]+(3 \[Lambda]VL[3]^2)/Sqrt[\[Mu]sqSU2]+(4 \[Lambda]VL[4]^2)/(Sqrt[\[Mu]sqSU2]+Sqrt[\[Mu]sqU1])))/\[Pi]),g13dUS^2->g13d^2,g23dUS^2->g23d^2-g23d^4/(24 \[Pi] Sqrt[\[Mu]sqSU2]),g33dUS^2->g33d^2-g33d^4/(16 \[Pi] Sqrt[\[Mu]sqSU3])}
]];


AppendTo[testList,
TestCreate[PrintScalarMassUS["LO"],
	{m13dUS->m13d-1/(16 \[Pi] \[Mu]\[Sigma]3d) (8 \[Pi] \[Mu]13d \[Mu]m3d+\[Mu]m3d^2 Sqrt[\[Mu]\[Sigma]3d]+2 \[Mu]\[Sigma]3d (\[Lambda]m3d Sqrt[\[Mu]\[Sigma]3d]+8 Sqrt[\[Mu]sqSU3] \[Lambda]VL[1]+Sqrt[\[Mu]sqU1] \[Lambda]VL[2]+3 Sqrt[\[Mu]sqSU2] \[Lambda]VL[3]))}
]];
AppendTo[testList,
TestCreate[PrintScalarMassUS["NLO"],
	{m13dUS->1/(1536 \[Pi]^2) (72 \[Lambda]m3d \[Lambda]\[Sigma]3d-(8 \[Pi] \[Mu]13d \[Mu]m3d^3)/\[Mu]\[Sigma]3d^(5/2)-\[Mu]m3d^4/\[Mu]\[Sigma]3d^2+(16 m13d \[Pi] \[Mu]m3d^2)/\[Mu]\[Sigma]3d^(3/2)-(144 \[Lambda]1H3d \[Mu]m3d^2)/\[Mu]\[Sigma]3d+(6 \[Lambda]m3d \[Mu]m3d (\[Mu]33d+\[Mu]m3d))/\[Mu]\[Sigma]3d+(8 (4 \[Lambda]m3d \[Mu]33d^2+18 \[Lambda]1H3d \[Mu]m3d^2+3 \[Lambda]m3d \[Mu]m3d^2))/\[Mu]\[Sigma]3d+(48 \[Lambda]m3d \[Mu]m3d (\[Mu]m3d Log[2]+\[Mu]33d Log[9/4]))/\[Mu]\[Sigma]3d-72 g23d^4 Log[\[Mu]3/(2 Sqrt[\[Mu]sqSU2])]-24 \[Lambda]m3d^2 (1+2 Log[\[Mu]3/(2 Sqrt[\[Mu]\[Sigma]3d])])-(288 \[Lambda]1H3d \[Mu]m3d^2 Log[\[Mu]3/Sqrt[\[Mu]\[Sigma]3d]])/\[Mu]\[Sigma]3d+(3 \[Mu]m3d^3 (\[Mu]m3d-8 \[Mu]33d Log[4/3]+\[Mu]m3d Log[4]+2 \[Mu]m3d Log[\[Mu]3/Sqrt[\[Mu]\[Sigma]3d]]))/\[Mu]\[Sigma]3d^2+576 g33d^2 (1+4 Log[\[Mu]3/(2 Sqrt[\[Mu]sqSU3])]) \[Lambda]VL[1]-192 (1+2 Log[\[Mu]3/(2 Sqrt[\[Mu]sqSU3])]) \[Lambda]VL[1]^2-24 (1+2 Log[\[Mu]3/(2 Sqrt[\[Mu]sqU1])]) \[Lambda]VL[2]^2+144 g23d^2 (1+4 Log[\[Mu]3/(2 Sqrt[\[Mu]sqSU2])]) \[Lambda]VL[3]-72 (1+2 Log[\[Mu]3/(2 Sqrt[\[Mu]sqSU2])]) \[Lambda]VL[3]^2-(2 \[Mu]m3d^2 (\[Lambda]m3d Sqrt[\[Mu]\[Sigma]3d]+8 Sqrt[\[Mu]sqSU3] \[Lambda]VL[1]+Sqrt[\[Mu]sqU1] \[Lambda]VL[2]+3 Sqrt[\[Mu]sqSU2] \[Lambda]VL[3]))/\[Mu]\[Sigma]3d^(3/2)-144 (1+2 Log[\[Mu]3/(Sqrt[\[Mu]sqSU2]+Sqrt[\[Mu]sqU1])]) \[Lambda]VL[4]^2+(32 \[Lambda]VL[1] (9 Sqrt[\[Mu]sqSU2] \[Lambda]VLL[1]+10 Sqrt[\[Mu]sqSU3] \[Lambda]VLL[2]+3 Sqrt[\[Mu]sqU1] \[Lambda]VLL[7]))/Sqrt[\[Mu]sqSU3]+(12 \[Lambda]VL[3] (24 Sqrt[\[Mu]sqSU3] \[Lambda]VLL[1]+5 Sqrt[\[Mu]sqSU2] \[Lambda]VLL[3]+3 Sqrt[\[Mu]sqU1] \[Lambda]VLL[8]))/Sqrt[\[Mu]sqSU2]+(12 \[Lambda]VL[2] (8 Sqrt[\[Mu]sqSU3] \[Lambda]VLL[7]+3 Sqrt[\[Mu]sqSU2] \[Lambda]VLL[8]+Sqrt[\[Mu]sqU1] \[Lambda]VLL[9]))/Sqrt[\[Mu]sqU1])}
]];


AppendTo[testList,
TestCreate[BetaFunctions3DUS[],
	{m13dUS->1/(256 \[Pi]^2) (-51 g23dUS^4+5 g13dUS^4 Y\[Phi]^4+18 g23dUS^2 (g13dUS^2 Y\[Phi]^2-8 \[Lambda]1H3dUS)-48 g13dUS^2 Y\[Phi]^2 \[Lambda]1H3dUS+192 \[Lambda]1H3dUS^2)}
]];


AppendTo[testList,
TestCreate[PrintTadpolesUS["LO"],
	{}
]];


report=TestReport[testList]
report["ResultsDataset"]


(* ::Subsection:: *)
(*Test scalar 1 softer*)


PerformDRsoft[{1,2,3,4}]


testList={};


AppendTo[testList,
TestCreate[PrintCouplingsUS[],
	{\[Lambda]\[Sigma]3dUS->-(((192 m13d^2 \[Lambda]m3d^2+16 m13d (-3 \[Lambda]m3d+\[Lambda]\[Sigma]3d) \[Mu]m3d^2-3 \[Mu]m3d^4)/m13d^(5/2)+48 (-32 \[Pi] \[Lambda]\[Sigma]3d+(3 \[Lambda]VL[5]^2)/Sqrt[\[Mu]sqSU2]+\[Lambda]VL[6]^2/Sqrt[\[Mu]sqU1]))/(1536 \[Pi])),g23dUS^2->g23d^2-(g23d^4 (1/Sqrt[m13d]+2/Sqrt[\[Mu]sqSU2]))/(48 \[Pi]),g33dUS^2->g33d^2-g33d^4/(16 \[Pi] Sqrt[\[Mu]sqSU3]),\[Mu]33dUS->(128 m13d^(3/2) \[Pi] \[Mu]33d-24 m13d \[Lambda]m3d \[Mu]m3d+\[Mu]m3d^2 (\[Mu]33d+3 \[Mu]m3d))/(128 m13d^(3/2) \[Pi])}
]];


AppendTo[testList,
TestCreate[PrintScalarMassUS["LO"],
	{\[Mu]\[Sigma]3dUS->\[Mu]\[Sigma]3d-(8 m13d \[Lambda]m3d+\[Mu]m3d^2+2 Sqrt[m13d] (3 Sqrt[\[Mu]sqSU2] \[Lambda]VL[5]+Sqrt[\[Mu]sqU1] \[Lambda]VL[6]))/(16 Sqrt[m13d] \[Pi])}
]];
AppendTo[testList,
TestCreate[PrintScalarMassUS["NLO"],
	{\[Mu]\[Sigma]3dUS->1/(1536 \[Pi]^2) (576 \[Lambda]1H3d \[Lambda]m3d+(12 \[Lambda]m3d \[Mu]m3d^2)/m13d+(24 (\[Lambda]m3d+3 \[Lambda]\[Sigma]3d) \[Mu]m3d^2)/m13d-\[Mu]m3d^4/(2 m13d^2)+(96 \[Lambda]m3d \[Mu]m3d (2 \[Mu]33d+\[Mu]m3d))/m13d-(3 \[Mu]m3d^3 (8 \[Mu]33d+\[Mu]m3d))/m13d^2+(8 \[Pi] \[Mu]m3d^2 \[Mu]\[Sigma]3d)/m13d^(3/2)-96 \[Lambda]m3d^2 (1+2 Log[\[Mu]3/(2 Sqrt[m13d])])+24 (3 g23d^2+g13d^2 Y\[Phi]^2) \[Lambda]m3d (1+4 Log[\[Mu]3/(2 Sqrt[m13d])])+144 g23d^2 (1+4 Log[\[Mu]3/(2 Sqrt[\[Mu]sqSU2])]) \[Lambda]VL[5]-72 (1+2 Log[\[Mu]3/(2 Sqrt[\[Mu]sqSU2])]) \[Lambda]VL[5]^2-24 (1+2 Log[\[Mu]3/(2 Sqrt[\[Mu]sqU1])]) \[Lambda]VL[6]^2-(\[Mu]m3d^2 (4 Sqrt[m13d] \[Lambda]m3d+3 Sqrt[\[Mu]sqSU2] \[Lambda]VL[5]+Sqrt[\[Mu]sqU1] \[Lambda]VL[6]))/m13d^(3/2)+(12 \[Lambda]VL[5] (24 Sqrt[\[Mu]sqSU3] \[Lambda]VLL[1]+5 Sqrt[\[Mu]sqSU2] \[Lambda]VLL[3]+3 Sqrt[\[Mu]sqU1] \[Lambda]VLL[8]))/Sqrt[\[Mu]sqSU2]+(12 \[Lambda]VL[6] (8 Sqrt[\[Mu]sqSU3] \[Lambda]VLL[7]+3 Sqrt[\[Mu]sqSU2] \[Lambda]VLL[8]+Sqrt[\[Mu]sqU1] \[Lambda]VLL[9]))/Sqrt[\[Mu]sqU1])}
]];


AppendTo[testList,
TestCreate[BetaFunctions3DUS[],
	{\[Mu]\[Sigma]3dUS->(3 \[Lambda]\[Sigma]3dUS^2)/(8 \[Pi]^2),\[Mu]13dUS->(\[Lambda]\[Sigma]3dUS \[Mu]33dUS)/(8 \[Pi]^2)}
]];


AppendTo[testList,
TestCreate[PrintTadpolesUS["LO"],
	{\[Mu]13dUS->\[Mu]13d-(Sqrt[m13d] \[Mu]m3d)/(4 \[Pi])-(\[Mu]13d \[Mu]m3d^2)/(384 m13d^(3/2) \[Pi])}
]];


AppendTo[testList,
TestCreate[PrintPressureUS["LO"],
	m13d^(3/2)/(3 \[Pi])+\[Mu]sqSU2^(3/2)/(4 \[Pi])+(2 \[Mu]sqSU3^(3/2))/(3 \[Pi])+\[Mu]sqU1^(3/2)/(12 \[Pi])
]];
AppendTo[testList,
TestCreate[PrintPressureUS["NLO"],
	-((3 g23d^2 (6 m13d+8 m13d Log[\[Mu]3/(2 Sqrt[m13d])]))/(128 \[Pi]^2))-(g13d^2 Y\[Phi]^2 (6 m13d+8 m13d Log[\[Mu]3/(2 Sqrt[m13d])]))/(128 \[Pi]^2)-(3 g23d^2 (6 \[Mu]sqSU2+8 \[Mu]sqSU2 Log[\[Mu]3/(2 Sqrt[\[Mu]sqSU2])]))/(64 \[Pi]^2)-(3 g33d^2 (6 \[Mu]sqSU3+8 \[Mu]sqSU3 Log[\[Mu]3/(2 Sqrt[\[Mu]sqSU3])]))/(16 \[Pi]^2)-1/(128 \[Pi]^2) (48 m13d \[Lambda]1H3d+48 Sqrt[\[Mu]sqSU2] Sqrt[\[Mu]sqSU3] \[Lambda]VLL[1]+80/3 \[Mu]sqSU3 \[Lambda]VLL[2]+5 \[Mu]sqSU2 \[Lambda]VLL[3]+16 Sqrt[\[Mu]sqSU3] Sqrt[\[Mu]sqU1] \[Lambda]VLL[7]+6 Sqrt[\[Mu]sqSU2] Sqrt[\[Mu]sqU1] \[Lambda]VLL[8]+\[Mu]sqU1 \[Lambda]VLL[9])
]];


report=TestReport[testList]
report["ResultsDataset"]




