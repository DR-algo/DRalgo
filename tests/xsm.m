(* ::Package:: *)

Quit[];


If[$InputFileName=="",
	SetDirectory[NotebookDirectory[]],
	SetDirectory[DirectoryName[$InputFileName]]
];
(*Put this if you want to create multiple model-files with the same kernel*)
(*DRalgo`$GroupMathMultipleModels=True;*)

DRalgo`DRalgo`$LoadGroupMath=True;
(*DRalgo`DRalgo`$InstallGroupMath=True;*)

Check[
    Get["../Kernel/DRalgo.m"],
    Message[Get::noopen, "DRalgo` at "<>ToString[$UserBaseDirectory]<>"/Applications"];
    Abort[];
]


(* ::Chapter:: *)
(*SM+sr1*)


(*see 2102.11145 [hep-ph]*)


(* ::Section:: *)
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
	+c3/3*CubicTerm2
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
(*Dimensional Reduction Mode 0*)


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->True,Mode->0];
PosFermion=PrintFermionRepPositions[];
FermionMat=Table[{nF,i},{i,PosFermion}];
DefineNF[FermionMat]
PerformDRhard[]


(* ::Section:: *)
(*Test Mode 0*)


(* ::Subsection:: *)
(*Test hard*)


testList={};


(* ::Subsubsection:: *)
(*Beta functions*)


AppendTo[testList,
TestCreate[BetaFunctions4D[],
	{g3^2->0,g2^2->0,g1^2->0,\[Lambda]1H->0,\[Lambda]m->0,\[Lambda]\[Sigma]->0,c3->0,\[Mu]m->0,yt1->0,m1->0,\[Mu]\[Sigma]->0,\[Mu]1->0}
]];


(* ::Subsubsection:: *)
(*Dimension 4 matching relations*)


AppendTo[testList,
TestCreate[PrintCouplings[],
	{g13d^2->g1^2 T,g23d^2->g2^2 T,g33d^2->g3^2 T,\[Lambda]1H3d->T \[Lambda]1H,\[Lambda]m3d->T \[Lambda]m,\[Lambda]\[Sigma]3d->T \[Lambda]\[Sigma],c33d->Sqrt[T] c3,\[Mu]m3d->Sqrt[T] \[Mu]m}
]];


AppendTo[testList,
TestCreate[PrintTadpoles["LO"],
	{\[Mu]13d -> (\[Mu]1 + (T^2*(c3 + \[Mu]m))/12)/Sqrt[T]}
]];
AppendTo[testList,
TestCreate[PrintTadpoles["NLO"],
	{\[Mu]13d->0}
]];


AppendTo[testList,
TestCreate[PrintTemporalScalarCouplings[],
	{lamVL[1]->(g2^2 T)/2,lamVL[2]->-(1/2) g1 g2 T Y\[Phi],lamVL[3]->1/2 g1^2 T Y\[Phi]^2}
]];


AppendTo[testList,
TestCreate[PrintDebyeMass["LO"],
	{musqSU2->1/6 g2^2 (5+2 nF) T^2,musqSU3->1/3 g3^2 (3+nF) T^2,musqU1->1/24 g1^2 T^2 (nF (3 Yd^2+Ye^2+2 Yl^2+6 Yq^2+3 Yu^2)+4 Y\[Phi]^2)}
]];
AppendTo[testList,
TestCreate[PrintDebyeMass["NLO"],
	{musqSU2->0,musqSU3->0,musqU1->0}
]];


AppendTo[testList,
TestCreate[PrintScalarMass["LO"],
	{m13d->m1+1/48 T^2 (9 g2^2+12 yt1^2+3 g1^2 Y\[Phi]^2+24 \[Lambda]1H+2 \[Lambda]m),\[Mu]\[Sigma]3d->1/12 T^2 (2 \[Lambda]m+3 \[Lambda]\[Sigma])+\[Mu]\[Sigma]}
]];
AppendTo[testList,
TestCreate[PrintScalarMass["NLO"],
	{m13d->0,\[Mu]\[Sigma]3d->0}
]];


(* ::Subsubsection:: *)
(*Symmetric pressure*)


AppendTo[testList,
TestCreate[PrintPressure["LO"],
	(29*Pi^2*T^4)/90 + (7*nF*Pi^2*T^4)/24
]];


(* ::Subsubsection:: *)
(*Report*)


report=TestReport[testList]
report["ResultsDataset"]


(* ::Section:: *)
(*Dimensional Reduction Mode 3*)


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->True,Mode->3];
DefineDim6[\[Lambda]6];
PosFermion=PrintFermionRepPositions[];
FermionMat=Table[{nF,i},{i,PosFermion}];
DefineNF[FermionMat]
PerformDRhard[]


(* ::Section:: *)
(*Tests Mode 3*)


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
	{g1^2->(g1^2 (-18 yt1^2 (Yq-Yu)^2+(18 yt1^2+g1^2 nF (3 Yd^2+Ye^2+2 Yl^2+6 Yq^2+3 Yu^2)) Y\[Phi]^2+g1^2 Y\[Phi]^4))/(48 \[Pi]^2 Y\[Phi]^2),g2^2->(g2^4 (-43+8 nF))/(48 \[Pi]^2),g3^2->(g3^4 (-33+4 nF))/(24 \[Pi]^2),\[Lambda]1H->(9 g2^4-48 yt1^4+3 g1^4 Y\[Phi]^4+6 g2^2 (g1^2 Y\[Phi]^2-12 \[Lambda]1H)-24 g1^2 Y\[Phi]^2 \[Lambda]1H+96 \[Lambda]1H (yt1^2+2 \[Lambda]1H)+4 \[Lambda]m^2)/(128 \[Pi]^2),\[Lambda]m->(\[Lambda]m (-9 g2^2+12 yt1^2-3 g1^2 Y\[Phi]^2+24 \[Lambda]1H+8 \[Lambda]m+12 \[Lambda]\[Sigma]))/(32 \[Pi]^2),\[Lambda]\[Sigma]->(\[Lambda]m^2+9 \[Lambda]\[Sigma]^2)/(8 \[Pi]^2),c3->(3 (6 \[Lambda]\[Sigma] c3+\[Lambda]m \[Mu]m))/(16 \[Pi]^2),\[Mu]m->(-3 (3 g2^2-4 yt1^2+g1^2 Y\[Phi]^2-8 \[Lambda]1H) \[Mu]m+8 \[Lambda]m (c3+\[Mu]m))/(32 \[Pi]^2),yt1->(yt1 (-9 g2^2-32 g3^2+18 yt1^2-3 g1^2 (2 Yq Yu+Y\[Phi]^2)))/(64 \[Pi]^2),m1->(-3 m1 (3 g2^2-4 yt1^2+g1^2 Y\[Phi]^2-8 \[Lambda]1H)+\[Mu]m^2+2 \[Lambda]m \[Mu]\[Sigma])/(32 \[Pi]^2),\[Mu]\[Sigma]->(4 m1 \[Lambda]m+4 c3^2+\[Mu]m^2+6 \[Lambda]\[Sigma] \[Mu]\[Sigma])/(16 \[Pi]^2),\[Mu]1->(m1 \[Mu]m+c3 \[Mu]\[Sigma])/(8 \[Pi]^2)}
]];


(* ::Subsubsection:: *)
(*Dimension 4 matching relations*)


AppendTo[testList,
TestCreate[PrintCouplings[]/.Replacements//Simplify,
	{g13d^2 -> (g1^2*T*(18*Lf*yt1^2*(Yq - Yu)^2 + (96*Pi^2 - Lf*(18*yt1^2 + g1^2*nF*(3*Yd^2 + Ye^2 + 2*Yl^2 + 6*Yq^2 + 3*Yu^2)))*Y\[Phi]^2 - g1^2*Lb*Y\[Phi]^4))/(96*Pi^2*Y\[Phi]^2), g23d^2 -> g2^2*T + (g2^4*(4 + 43*Lb - 8*Lf*nF)*T)/(96*Pi^2), g33d^2 -> g3^2*T + (g3^4*(3 + 33*Lb - 4*Lf*nF)*T)/(48*Pi^2), \[Lambda]1H3d -> (T*(48*Lf*yt1^4 + (2 - 3*Lb)*(3*g2^4 + 2*g1^2*g2^2*Y\[Phi]^2 + g1^4*Y\[Phi]^4) + 256*Pi^2*\[Lambda]1H + 24*(3*g2^2*Lb - 4*Lf*yt1^2 + g1^2*Lb*Y\[Phi]^2)*\[Lambda]1H - 4*Lb*(48*\[Lambda]1H^2 + \[Lambda]m^2)))/(256*Pi^2), \[Lambda]m3d -> (T*\[Lambda]m*(9*g2^2*Lb + 64*Pi^2 - 12*Lf*yt1^2 + 3*g1^2*Lb*Y\[Phi]^2 - 4*Lb*(6*\[Lambda]1H + 2*\[Lambda]m + 3*\[Lambda]\[Sigma])))/(64*Pi^2), \[Lambda]\[Sigma]3d -> T*\[Lambda]\[Sigma] - (Lb*T*(\[Lambda]m^2 + 9*\[Lambda]\[Sigma]^2))/(16*Pi^2), c33d -> (Sqrt[T]*(2*c3 - (3*Lb*(6*c3*\[Lambda]\[Sigma] + \[Lambda]m*\[Mu]m))/(16*Pi^2)))/2, \[Mu]m3d -> (Sqrt[T]*(-8*c3*Lb*\[Lambda]m + (9*g2^2*Lb + 64*Pi^2 - 12*Lf*yt1^2 + 3*g1^2*Lb*Y\[Phi]^2 - 8*Lb*(3*\[Lambda]1H + \[Lambda]m))*\[Mu]m))/(64*Pi^2)}
]];


AppendTo[testList,
TestCreate[PrintTadpoles["LO"]/.Replacements,
	{\[Mu]13d -> (c3*T^2 + 12*\[Mu]1 + T^2*\[Mu]m)/(12*Sqrt[T])}
]];
AppendTo[testList,
TestCreate[PrintTadpoles["NLO"]/.Replacements//Expand,
	{\[Mu]13d->1/(768 \[Pi]^2 Sqrt[T]) (-Lb (T^2 (-12 \[Lambda]\[Sigma] c3+2 \[Lambda]m (4 c3-5 \[Mu]m)+3 (9 g2^2+6 yt1^2+3 g1^2 Y\[Phi]^2+8 \[Lambda]1H) \[Mu]m)+48 (m1 \[Mu]m+c3 \[Mu]\[Sigma]))+2 T^2 (3 Lf yt1^2 \[Mu]m+g1^2 Y\[Phi]^2 \[Mu]m+6 EulerGamma (-4 \[Lambda]\[Sigma] c3+(3 g2^2+g1^2 Y\[Phi]^2-2 \[Lambda]m) \[Mu]m)+3 g2^2 \[Mu]m (1-72 Log[Glaisher])+288 \[Lambda]\[Sigma] c3 Log[Glaisher]-72 g1^2 Y\[Phi]^2 \[Mu]m Log[Glaisher]+144 \[Lambda]m \[Mu]m Log[Glaisher])+24 Sqrt[T] (4 \[Lambda]\[Sigma]3d c33d-(3 g23d^2+g13d^2 Y\[Phi]^2-2 \[Lambda]m3d) \[Mu]m3d) Log[mu3/mu4])}//Expand
]];


AppendTo[testList,
TestCreate[PrintTemporalScalarCouplings[]/.Replacements,
	{lamVLL[1]->-((g2^2 g3^2 nF T)/(4 \[Pi]^2)),lamVLL[2]->(g3^4 (9-2 nF) T)/(4 \[Pi]^2),lamVLL[3]->(g2^4 (17-4 nF) T)/(8 \[Pi]^2),lamVLL[4]->-((I g1 g2^3 nF T (Yl+3 Yq))/(24 \[Pi]^2)),lamVLL[5]->-((I g1 g3^3 nF T (Yd-2 Yq+Yu))/(48 \[Pi]^2)),lamVLL[6]->-((g1 g3^3 nF T (Yd+2 Yq+Yu))/(16 \[Pi]^2)),lamVLL[7]->-((g1^2 g3^2 nF T (Yd^2+2 Yq^2+Yu^2))/(8 \[Pi]^2)),lamVLL[8]->(g1^2 g2^2 T (-nF (Yl^2+3 Yq^2)+Y\[Phi]^2))/(8 \[Pi]^2),lamVLL[9]->-((g1^4 T (nF (3 Yd^4+Ye^4+2 Yl^4+6 Yq^4+3 Yu^4)-2 Y\[Phi]^4))/(16 \[Pi]^2)),lamVVSL[1]->(g2^2 Sqrt[T] \[Mu]m)/(16 \[Pi]^2),lamVVSL[2]->(g1^2 Sqrt[T] Y\[Phi]^2 \[Mu]m)/(16 \[Pi]^2),lamVL[1]->-((g3^2 T yt1^2)/(4 \[Pi]^2)),lamVL[2]->1/192 g1^2 T (96 Y\[Phi]^2+1/\[Pi]^2 (-36 Lf Yq yt1^2 Yu+18 yt1^2 ((-2+Lf) Yu^2-Lf Y\[Phi]^2)+6 Yq^2 (3 (-2+Lf) yt1^2-g1^2 (-1+Lf) nF Y\[Phi]^2)+Y\[Phi]^2 (9 g2^2+g1^2 (-((-1+Lf) nF (3 Yd^2+Ye^2+2 Yl^2+3 Yu^2))-(-1+Lb) Y\[Phi]^2)+72 \[Lambda]1H))),lamVL[3]->(g2^2 T (g2^2 (51+43 Lb+8 nF-8 Lf nF)+96 \[Pi]^2-36 yt1^2+3 g1^2 Y\[Phi]^2+72 \[Lambda]1H))/(192 \[Pi]^2),lamVL[4]->1/384 g1 g2 T (192 Y\[Phi]+1/\[Pi]^2 (-36 (-2+Lf) Yq yt1^2+36 Lf yt1^2 (Yu-Y\[Phi])+g2^2 (-12+43 Lb+8 nF) Y\[Phi]-6 g1^2 (-1+Lf) nF Yq^2 Y\[Phi]-Lf nF (8 g2^2+g1^2 (3 Yd^2+Ye^2+2 Yl^2+3 Yu^2)) Y\[Phi]+g1^2 Y\[Phi] (nF (3 Yd^2+Ye^2+2 Yl^2+3 Yu^2)-(-4+Lb) Y\[Phi]^2)+48 Y\[Phi] \[Lambda]1H)),lamVL[5]->(g2^2 T \[Lambda]m)/(8 \[Pi]^2),lamVL[6]->(g1^2 T Y\[Phi]^2 \[Lambda]m)/(8 \[Pi]^2)}
]];


AppendTo[testList,
TestCreate[PrintDebyeMass["LO"]/.Replacements,
	{musqSU2->1/6 g2^2 (5+2 nF) T^2,musqSU3->1/3 g3^2 (3+nF) T^2,musqU1->1/24 g1^2 T^2 (nF (3 Yd^2+Ye^2+2 Yl^2+6 Yq^2+3 Yu^2)+4 Y\[Phi]^2)}
]];
AppendTo[testList,
TestCreate[PrintDebyeMass["NLO"]/.Replacements,
	{musqSU2 -> (g2^2*(144*m1 + T^2*(-144*g3^2*nF + 9*g1^2*(-(nF*(Yl^2 + 3*Yq^2)) + Y\[Phi]^2) + 6*(-3*yt1^2 + 12*\[Lambda]1H + \[Lambda]m) + g2^2*(207 - 4*EulerGamma*(5 + 2*nF)*(-43 + 8*nF) + 4*nF*(11 + 8*nF - 192*Log[2]))) + 4*g2^2*T^2*((39 + 8*(7 - 2*nF)*nF)*Log[mu4] + 2*(88 - 5*nF)*Log[mu4/(4*Pi*T)] + 8*nF*(-7 + 2*nF)*Log[Pi*T] - 39*Log[4*Pi*T])))/(1152*Pi^2), musqSU3 -> (g3^2*T^2*(-9*(6*g2^2*nF + 8*yt1^2 + g1^2*nF*(Yd^2 + 2*Yq^2 + Yu^2)) - 8*g3^2*(-45 + 2*EulerGamma*(3 + nF)*(-33 + 4*nF) + nF*(-3 - 4*nF + 132*Log[2])) + 16*g3^2*(99*Log[mu4/(4*Pi*T)] + nF*(-21 + 4*nF)*(-Log[mu4] + Log[Pi*T]))))/(1152*Pi^2), musqU1 -> (-27*g1^2*T^2*yt1^2*(Yq^2 + Yu^2) - (9*g1^2*nF*T^2*(6*g2^2*(Yl^2 + 3*Yq^2) + 16*g3^2*(Yd^2 + 2*Yq^2 + Yu^2) + g1^2*(3*Yd^4 + Ye^4 + 2*Yl^4 + 6*Yq^4 + 3*Yu^4)))/4 + 72*g1^2*m1*Y\[Phi]^2 + 18*g1^2*T^2*yt1^2*Y\[Phi]^2 + (9*T^2*(3*g1^2*g2^2*Y\[Phi]^2 + g1^4*Y\[Phi]^4))/2 + 3*g1^2*T^2*Y\[Phi]^2*(12*\[Lambda]1H + \[Lambda]m) - (g1^4*nF*T^2*(3*Yd^2 + Ye^2 + 2*Yl^2 + 6*Yq^2 + 3*Yu^2)*Y\[Phi]^2*(-1 + 5*EulerGamma + Log[1/(4*Pi^5)] + 5*Log[mu4/T]))/2 - (g1^4*nF^2*T^2*(3*Yd^2 + Ye^2 + 2*Yl^2 + 6*Yq^2 + 3*Yu^2)^2*(-1 + 2*EulerGamma + 2*Log[mu4] - 2*Log[Pi*T]))/4 - 2*g1^4*T^2*Y\[Phi]^4*(1 + EulerGamma + Log[mu4] - Log[4*Pi*T]))/(576*Pi^2)}
]];


AppendTo[testList,
TestCreate[PrintScalarMass["LO"]/.Replacements//Simplify,
	{m13d->m1+1/48 T^2 (9 g2^2+12 yt1^2+3 g1^2 Y\[Phi]^2+24 \[Lambda]1H+2 \[Lambda]m),\[Mu]\[Sigma]3d->1/12 T^2 (2 \[Lambda]m+3 \[Lambda]\[Sigma])+\[Mu]\[Sigma]}
]];
AppendTo[testList,
TestCreate[PrintScalarMass["NLO"]/.Replacements//Simplify,
	{m13d -> (-1728*Lf*m1*yt1^2 - 1152*g3^2*T^2*yt1^2 - 384*g3^2*Lb*T^2*yt1^2 + 1536*g3^2*Lf*T^2*yt1^2 - 36*g1^2*T^2*Yq^2*yt1^2 + 18*g1^2*Lb*T^2*Yq^2*yt1^2 + 90*g1^2*Lf*T^2*Yq^2*yt1^2 + 216*Lb*T^2*yt1^4 - 144*g1^2*T^2*Yq*yt1^2*Yu - 108*g1^2*Lb*T^2*Yq*yt1^2*Yu + 108*g1^2*Lf*T^2*Yq*yt1^2*Yu - 36*g1^2*T^2*yt1^2*Yu^2 + 18*g1^2*Lb*T^2*yt1^2*Yu^2 + 90*g1^2*Lf*T^2*yt1^2*Yu^2 + 432*g1^2*Lb*m1*Y\[Phi]^2 + 6*g1^4*nF*T^2*Yd^2*Y\[Phi]^2 - 27*g1^4*Lb*nF*T^2*Yd^2*Y\[Phi]^2 + 9*g1^4*Lf*nF*T^2*Yd^2*Y\[Phi]^2 + 2*g1^4*nF*T^2*Ye^2*Y\[Phi]^2 - 9*g1^4*Lb*nF*T^2*Ye^2*Y\[Phi]^2 + 3*g1^4*Lf*nF*T^2*Ye^2*Y\[Phi]^2 + 4*g1^4*nF*T^2*Yl^2*Y\[Phi]^2 - 18*g1^4*Lb*nF*T^2*Yl^2*Y\[Phi]^2 + 6*g1^4*Lf*nF*T^2*Yl^2*Y\[Phi]^2 + 12*g1^4*nF*T^2*Yq^2*Y\[Phi]^2 - 54*g1^4*Lb*nF*T^2*Yq^2*Y\[Phi]^2 + 18*g1^4*Lf*nF*T^2*Yq^2*Y\[Phi]^2 + 108*g1^2*Lb*T^2*yt1^2*Y\[Phi]^2 - 108*g1^2*Lf*T^2*yt1^2*Y\[Phi]^2 + 6*g1^4*nF*T^2*Yu^2*Y\[Phi]^2 - 27*g1^4*Lb*nF*T^2*Yu^2*Y\[Phi]^2 + 9*g1^4*Lf*nF*T^2*Yu^2*Y\[Phi]^2 + 2*g1^4*T^2*Y\[Phi]^4 - 126*EulerGamma*g1^4*T^2*Y\[Phi]^4 + 66*g1^4*Lb*T^2*Y\[Phi]^4 - 3456*Lb*m1*\[Lambda]1H - 1296*Lb*T^2*yt1^2*\[Lambda]1H - 432*Lf*T^2*yt1^2*\[Lambda]1H + 144*g1^2*T^2*Y\[Phi]^2*\[Lambda]1H + 864*EulerGamma*g1^2*T^2*Y\[Phi]^2*\[Lambda]1H - 432*g1^2*Lb*T^2*Y\[Phi]^2*\[Lambda]1H - 3456*EulerGamma*T^2*\[Lambda]1H^2 - 72*Lf*T^2*yt1^2*\[Lambda]m + 18*g1^2*Lb*T^2*Y\[Phi]^2*\[Lambda]m - 144*Lb*T^2*\[Lambda]1H*\[Lambda]m - 144*EulerGamma*T^2*\[Lambda]m^2 + 24*Lb*T^2*\[Lambda]m^2 - 72*Lb*T^2*\[Lambda]m*\[Lambda]\[Sigma] - 144*Lb*\[Mu]m^2 - 288*Lb*\[Lambda]m*\[Mu]\[Sigma] + 54*g2^2*(Lb*(24*m1 + T^2*(7*yt1^2 + 8*g1^2*Y\[Phi]^2 - 24*\[Lambda]1H + \[Lambda]m)) - T^2*((2 + Lf)*yt1^2 - 8*\[Lambda]1H*(1 + 6*EulerGamma - 72*Log[Glaisher]) + 2*g1^2*Y\[Phi]^2*(1 + 5*EulerGamma - 60*Log[Glaisher]))) + 6*g2^4*T^2*(167 + 243*EulerGamma + 8*nF + 12*Lf*nF - 3*Lb*(47 + 12*nF) - 2916*Log[Glaisher]) + 1512*g1^4*T^2*Y\[Phi]^4*Log[Glaisher] - 10368*g1^2*T^2*Y\[Phi]^2*\[Lambda]1H*Log[Glaisher] + 41472*T^2*\[Lambda]1H^2*Log[Glaisher] + 1728*T^2*\[Lambda]m^2*Log[Glaisher] - 36*Log[mu3/mu4]*(39*g23d^4 - 5*g13d^4*Y\[Phi]^4 + 48*g13d^2*Y\[Phi]^2*\[Lambda]1H3d + g23d^2*(-18*g13d^2*Y\[Phi]^2 + 144*\[Lambda]1H3d + 96*lamVL[3]) - 8*(24*\[Lambda]1H3d^2 + \[Lambda]m3d^2 - 48*g33d^2*lamVL[1] + 8*lamVL[1]^2 + lamVL[2]^2 + 3*lamVL[3]^2 + 6*lamVL[4]^2)))/(9216*Pi^2), \[Mu]\[Sigma]3d -> -1/384*(48*c3^2*Lb + 48*Lb*m1*\[Lambda]m - 6*g2^2*T^2*\[Lambda]m - 36*EulerGamma*g2^2*T^2*\[Lambda]m + 27*g2^2*Lb*T^2*\[Lambda]m + 18*Lb*T^2*yt1^2*\[Lambda]m - 6*Lf*T^2*yt1^2*\[Lambda]m - 2*g1^2*T^2*Y\[Phi]^2*\[Lambda]m - 12*EulerGamma*g1^2*T^2*Y\[Phi]^2*\[Lambda]m + 9*g1^2*Lb*T^2*Y\[Phi]^2*\[Lambda]m + 24*Lb*T^2*\[Lambda]1H*\[Lambda]m + 24*EulerGamma*T^2*\[Lambda]m^2 - 10*Lb*T^2*\[Lambda]m^2 + 12*Lb*T^2*\[Lambda]m*\[Lambda]\[Sigma] + 72*EulerGamma*T^2*\[Lambda]\[Sigma]^2 - 18*Lb*T^2*\[Lambda]\[Sigma]^2 + 12*Lb*\[Mu]m^2 + 72*Lb*\[Lambda]\[Sigma]*\[Mu]\[Sigma] + 432*g2^2*T^2*\[Lambda]m*Log[Glaisher] + 144*g1^2*T^2*Y\[Phi]^2*\[Lambda]m*Log[Glaisher] - 288*T^2*\[Lambda]m^2*Log[Glaisher] - 864*T^2*\[Lambda]\[Sigma]^2*Log[Glaisher] + 12*Log[mu3/mu4]*(2*g13d^2*Y\[Phi]^2*\[Lambda]m3d - 4*\[Lambda]m3d^2 - 12*\[Lambda]\[Sigma]3d^2 - 3*lamVL[5]^2 + 6*g23d^2*(\[Lambda]m3d + 2*lamVL[5]) - lamVL[6]^2))/Pi^2}
]];


(* ::Subsubsection:: *)
(*Dimension 6 matching relations*)


AppendTo[testList,
TestCreate[PrintCouplingsEffective[],
	{c613d->1/(6144 \[Pi]^4) (96 \[Pi]^2 T^2 (c61 (64 \[Pi]^2-36 Lf yt1^2+9 Lb (3 g2^2+g1^2 Y\[Phi]^2-24 \[Lambda]1H))-4 c63 Lb \[Lambda]m)+(9 g2^6+9 g1^2 g2^4 Y\[Phi]^2+9 g1^4 g2^2 Y\[Phi]^4+3 g1^6 Y\[Phi]^6+8 (-84 yt1^6+240 \[Lambda]1H^3+\[Lambda]m^3)) Zeta[3]),c623d->c62 T^2-(Lb T^2 (c64 \[Lambda]m+45 c62 \[Lambda]\[Sigma]))/(16 \[Pi]^2)+((\[Lambda]m^3+54 \[Lambda]\[Sigma]^3) Zeta[3])/(1536 \[Pi]^4),c633d->1/(256 \[Pi]^4) (-8 \[Pi]^2 T^2 (12 (c61+c64) Lb \[Lambda]m+c63 (-9 g2^2 Lb-32 \[Pi]^2+12 Lf yt1^2+Lb (-3 g1^2 Y\[Phi]^2+48 \[Lambda]1H+16 \[Lambda]m+6 \[Lambda]\[Sigma])))+\[Lambda]m (24 \[Lambda]1H^2+12 \[Lambda]1H \[Lambda]m+\[Lambda]m (2 \[Lambda]m+3 \[Lambda]\[Sigma])) Zeta[3]),c643d->1/(256 \[Pi]^4) (-4 \[Pi]^2 T^2 (12 (5 c62+c63) Lb \[Lambda]m+c64 (-64 \[Pi]^2+12 Lf yt1^2+Lb (-9 g2^2-3 g1^2 Y\[Phi]^2+24 \[Lambda]1H+32 \[Lambda]m+72 \[Lambda]\[Sigma])))+\[Lambda]m (3 \[Lambda]1H \[Lambda]m+(\[Lambda]m+3 \[Lambda]\[Sigma])^2) Zeta[3])}
]];


(* ::Subsubsection::Closed:: *)
(*Symmetric pressure*)


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
	-1/4423680*(-276480*Lb*m1^2 + 69120*g2^2*m1*T^2 + 414720*EulerGamma*g2^2*m1*T^2 - 
	311040*g2^2*Lb*m1*T^2 + 846648*g2^4*T^4 - 769440*EulerGamma*g2^4*T^4 + 3612672*g3^4*T^4 -
	3502080*EulerGamma*g3^4*T^4 + 303000*g2^4*Lb*T^4 + 1474560*g3^4*Lb*T^4 + 131544*g2^4*nF*T^4 -
	657600*EulerGamma*g2^4*nF*T^4 - 302400*g2^2*g3^2*nF*T^4 - 622080*EulerGamma*g2^2*g3^2*nF*T^4 +
	291072*g3^4*nF*T^4 - 2242560*EulerGamma*g3^4*nF*T^4 + 281280*g2^4*Lb*nF*T^4 +
	241920*g2^2*g3^2*Lb*nF*T^4 + 929280*g3^4*Lb*nF*T^4 - 10080*g2^4*Lf*nF*T^4 +
	69120*g2^2*g3^2*Lf*nF*T^4 + 7680*g3^4*Lf*nF*T^4 + 28800*g2^4*nF^2*T^4 -
	168960*EulerGamma*g2^4*nF^2*T^4 + 76800*g3^4*nF^2*T^4 - 450560*EulerGamma*g3^4*nF^2*T^4 +
	96000*g2^4*Lb*nF^2*T^4 + 256000*g3^4*Lb*nF^2*T^4 - 23040*g2^4*Lf*nF^2*T^4 -
	61440*g3^4*Lf*nF^2*T^4 - 50400*g1^2*g3^2*nF*T^4*Yd^2 - 103680*EulerGamma*g1^2*g3^2*nF*T^4*Yd^2 +
	40320*g1^2*g3^2*Lb*nF*T^4*Yd^2 + 11520*g1^2*g3^2*Lf*nF*T^4*Yd^2 - 4725*g1^4*nF*T^4*Yd^4 -
	9720*EulerGamma*g1^4*nF*T^4*Yd^4 + 3780*g1^4*Lb*nF*T^4*Yd^4 + 1080*g1^4*Lf*nF*T^4*Yd^4 +
	1350*g1^4*nF^2*T^4*Yd^4 - 7920*EulerGamma*g1^4*nF^2*T^4*Yd^4 + 4500*g1^4*Lb*nF^2*T^4*Yd^4 -
	1080*g1^4*Lf*nF^2*T^4*Yd^4 + 900*g1^4*nF^2*T^4*Yd^2*Ye^2 - 5280*EulerGamma*g1^4*nF^2*T^4*Yd^2*
	Ye^2 + 3000*g1^4*Lb*nF^2*T^4*Yd^2*Ye^2 - 720*g1^4*Lf*nF^2*T^4*Yd^2*Ye^2 -
	1575*g1^4*nF*T^4*Ye^4 - 3240*EulerGamma*g1^4*nF*T^4*Ye^4 + 1260*g1^4*Lb*nF*T^4*Ye^4 +
	360*g1^4*Lf*nF*T^4*Ye^4 + 150*g1^4*nF^2*T^4*Ye^4 - 880*EulerGamma*g1^4*nF^2*T^4*Ye^4 +
	500*g1^4*Lb*nF^2*T^4*Ye^4 - 120*g1^4*Lf*nF^2*T^4*Ye^4 - 18900*g1^2*g2^2*nF*T^4*Yl^2 -
	38880*EulerGamma*g1^2*g2^2*nF*T^4*Yl^2 + 15120*g1^2*g2^2*Lb*nF*T^4*Yl^2 +
	4320*g1^2*g2^2*Lf*nF*T^4*Yl^2 + 1800*g1^4*nF^2*T^4*Yd^2*Yl^2 -
	10560*EulerGamma*g1^4*nF^2*T^4*Yd^2*Yl^2 + 6000*g1^4*Lb*nF^2*T^4*Yd^2*Yl^2 -
	1440*g1^4*Lf*nF^2*T^4*Yd^2*Yl^2 + 600*g1^4*nF^2*T^4*Ye^2*Yl^2 -
	3520*EulerGamma*g1^4*nF^2*T^4*Ye^2*Yl^2 + 2000*g1^4*Lb*nF^2*T^4*Ye^2*Yl^2 -
	480*g1^4*Lf*nF^2*T^4*Ye^2*Yl^2 - 3150*g1^4*nF*T^4*Yl^4 - 6480*EulerGamma*g1^4*nF*T^4*Yl^4 +
	2520*g1^4*Lb*nF*T^4*Yl^4 + 720*g1^4*Lf*nF*T^4*Yl^4 + 600*g1^4*nF^2*T^4*Yl^4 -
	3520*EulerGamma*g1^4*nF^2*T^4*Yl^4 + 2000*g1^4*Lb*nF^2*T^4*Yl^4 - 480*g1^4*Lf*nF^2*T^4*Yl^4 -
	56700*g1^2*g2^2*nF*T^4*Yq^2 - 116640*EulerGamma*g1^2*g2^2*nF*T^4*Yq^2 -
	100800*g1^2*g3^2*nF*T^4*Yq^2 - 207360*EulerGamma*g1^2*g3^2*nF*T^4*Yq^2 +
	45360*g1^2*g2^2*Lb*nF*T^4*Yq^2 + 80640*g1^2*g3^2*Lb*nF*T^4*Yq^2 +
	12960*g1^2*g2^2*Lf*nF*T^4*Yq^2 + 23040*g1^2*g3^2*Lf*nF*T^4*Yq^2 + 5400*g1^4*nF^2*T^4*Yd^2*Yq^2 -
	31680*EulerGamma*g1^4*nF^2*T^4*Yd^2*Yq^2 + 18000*g1^4*Lb*nF^2*T^4*Yd^2*Yq^2 -
	4320*g1^4*Lf*nF^2*T^4*Yd^2*Yq^2 + 1800*g1^4*nF^2*T^4*Ye^2*Yq^2 -
	10560*EulerGamma*g1^4*nF^2*T^4*Ye^2*Yq^2 + 6000*g1^4*Lb*nF^2*T^4*Ye^2*Yq^2 -
	1440*g1^4*Lf*nF^2*T^4*Ye^2*Yq^2 + 3600*g1^4*nF^2*T^4*Yl^2*Yq^2 -
	21120*EulerGamma*g1^4*nF^2*T^4*Yl^2*Yq^2 + 12000*g1^4*Lb*nF^2*T^4*Yl^2*Yq^2 -
	2880*g1^4*Lf*nF^2*T^4*Yl^2*Yq^2 - 9450*g1^4*nF*T^4*Yq^4 - 19440*EulerGamma*g1^4*nF*T^4*Yq^4 +
	7560*g1^4*Lb*nF*T^4*Yq^4 + 2160*g1^4*Lf*nF*T^4*Yq^4 + 5400*g1^4*nF^2*T^4*Yq^4 -
	31680*EulerGamma*g1^4*nF^2*T^4*Yq^4 + 18000*g1^4*Lb*nF^2*T^4*Yq^4 - 4320*g1^4*Lf*nF^2*T^4*Yq^4 -
	207360*Lb*m1*T^2*yt1^2 + 69120*Lf*m1*T^2*yt1^2 - 85968*g2^2*T^4*yt1^2 -
	103680*EulerGamma*g2^2*T^4*yt1^2 - 651264*g3^2*T^4*yt1^2 - 368640*EulerGamma*g3^2*T^4*yt1^2 -
	6480*g2^2*Lb*T^4*yt1^2 + 23040*g3^2*Lb*T^4*yt1^2 + 32400*g2^2*Lf*T^4*yt1^2 +
	161280*g3^2*Lf*T^4*yt1^2 - 33048*g1^2*T^4*Yq^2*yt1^2 - 48960*EulerGamma*g1^2*T^4*Yq^2*yt1^2 +
	15840*g1^2*Lb*T^4*Yq^2*yt1^2 + 8640*g1^2*Lf*T^4*Yq^2*yt1^2 - 45576*T^4*yt1^4 -
	423360*EulerGamma*T^4*yt1^4 + 190080*Lb*T^4*yt1^4 + 4320*Lf*T^4*yt1^4 -
	56016*g1^2*T^4*Yq*yt1^2*Yu + 28800*EulerGamma*g1^2*T^4*Yq*yt1^2*Yu -
	27360*g1^2*Lb*T^4*Yq*yt1^2*Yu + 12960*g1^2*Lf*T^4*Yq*yt1^2*Yu - 50400*g1^2*g3^2*nF*T^4*Yu^2 -
	103680*EulerGamma*g1^2*g3^2*nF*T^4*Yu^2 + 40320*g1^2*g3^2*Lb*nF*T^4*Yu^2 +
	11520*g1^2*g3^2*Lf*nF*T^4*Yu^2 + 2700*g1^4*nF^2*T^4*Yd^2*Yu^2 -
	15840*EulerGamma*g1^4*nF^2*T^4*Yd^2*Yu^2 + 9000*g1^4*Lb*nF^2*T^4*Yd^2*Yu^2 -
	2160*g1^4*Lf*nF^2*T^4*Yd^2*Yu^2 + 900*g1^4*nF^2*T^4*Ye^2*Yu^2 -
	5280*EulerGamma*g1^4*nF^2*T^4*Ye^2*Yu^2 + 3000*g1^4*Lb*nF^2*T^4*Ye^2*Yu^2 -
	720*g1^4*Lf*nF^2*T^4*Ye^2*Yu^2 + 1800*g1^4*nF^2*T^4*Yl^2*Yu^2 -
	10560*EulerGamma*g1^4*nF^2*T^4*Yl^2*Yu^2 + 6000*g1^4*Lb*nF^2*T^4*Yl^2*Yu^2 -
	1440*g1^4*Lf*nF^2*T^4*Yl^2*Yu^2 + 5400*g1^4*nF^2*T^4*Yq^2*Yu^2 -
	31680*EulerGamma*g1^4*nF^2*T^4*Yq^2*Yu^2 + 18000*g1^4*Lb*nF^2*T^4*Yq^2*Yu^2 -
	4320*g1^4*Lf*nF^2*T^4*Yq^2*Yu^2 - 33048*g1^2*T^4*yt1^2*Yu^2 -
	48960*EulerGamma*g1^2*T^4*yt1^2*Yu^2 + 15840*g1^2*Lb*T^4*yt1^2*Yu^2 +
	8640*g1^2*Lf*T^4*yt1^2*Yu^2 - 4725*g1^4*nF*T^4*Yu^4 - 9720*EulerGamma*g1^4*nF*T^4*Yu^4 +
	3780*g1^4*Lb*nF*T^4*Yu^4 + 1080*g1^4*Lf*nF*T^4*Yu^4 + 1350*g1^4*nF^2*T^4*Yu^4 -
	7920*EulerGamma*g1^4*nF^2*T^4*Yu^4 + 4500*g1^4*Lb*nF^2*T^4*Yu^4 - 1080*g1^4*Lf*nF^2*T^4*Yu^4 +
	4392*g1^2*T^4*Yq*yt1^2*Y\[Phi] - 4392*g1^2*T^4*yt1^2*Yu*Y\[Phi] + 23040*g1^2*m1*T^2*Y\[Phi]^2 +
	138240*EulerGamma*g1^2*m1*T^2*Y\[Phi]^2 - 103680*g1^2*Lb*m1*T^2*Y\[Phi]^2 +
	54576*g1^2*g2^2*T^4*Y\[Phi]^2 - 25920*EulerGamma*g1^2*g2^2*T^4*Y\[Phi]^2 +
	6480*g1^2*g2^2*Lb*T^4*Y\[Phi]^2 + 4722*g1^4*nF*T^4*Yd^2*Y\[Phi]^2 -
	17520*EulerGamma*g1^4*nF*T^4*Yd^2*Y\[Phi]^2 + 8220*g1^4*Lb*nF*T^4*Yd^2*Y\[Phi]^2 -
	900*g1^4*Lf*nF*T^4*Yd^2*Y\[Phi]^2 + 1574*g1^4*nF*T^4*Ye^2*Y\[Phi]^2 -
	5840*EulerGamma*g1^4*nF*T^4*Ye^2*Y\[Phi]^2 + 2740*g1^4*Lb*nF*T^4*Ye^2*Y\[Phi]^2 -
	300*g1^4*Lf*nF*T^4*Ye^2*Y\[Phi]^2 + 3148*g1^4*nF*T^4*Yl^2*Y\[Phi]^2 -
	11680*EulerGamma*g1^4*nF*T^4*Yl^2*Y\[Phi]^2 + 5480*g1^4*Lb*nF*T^4*Yl^2*Y\[Phi]^2 -
	600*g1^4*Lf*nF*T^4*Yl^2*Y\[Phi]^2 + 9444*g1^4*nF*T^4*Yq^2*Y\[Phi]^2 -
	35040*EulerGamma*g1^4*nF*T^4*Yq^2*Y\[Phi]^2 + 16440*g1^4*Lb*nF*T^4*Yq^2*Y\[Phi]^2 -
	1800*g1^4*Lf*nF*T^4*Yq^2*Y\[Phi]^2 + 8784*g1^2*T^4*yt1^2*Y\[Phi]^2 +
	14400*EulerGamma*g1^2*T^4*yt1^2*Y\[Phi]^2 - 18000*g1^2*Lb*T^4*yt1^2*Y\[Phi]^2 +
	2160*g1^2*Lf*T^4*yt1^2*Y\[Phi]^2 + 4722*g1^4*nF*T^4*Yu^2*Y\[Phi]^2 -
	17520*EulerGamma*g1^4*nF*T^4*Yu^2*Y\[Phi]^2 + 8220*g1^4*Lb*nF*T^4*Yu^2*Y\[Phi]^2 -
	900*g1^4*Lf*nF*T^4*Yu^2*Y\[Phi]^2 + 17144*g1^4*T^4*Y\[Phi]^4 -
	13600*EulerGamma*g1^4*T^4*Y\[Phi]^4 + 4760*g1^4*Lb*T^4*Y\[Phi]^4 -
	276480*Lb*m1*T^2*\[Lambda]1H + 34560*g2^2*T^4*\[Lambda]1H - 51840*g2^2*Lb*T^4*\[Lambda]1H +
	276480*EulerGamma*T^4*yt1^2*\[Lambda]1H - 241920*Lb*T^4*yt1^2*\[Lambda]1H +
	34560*Lf*T^4*yt1^2*\[Lambda]1H + 11520*g1^2*T^4*Y\[Phi]^2*\[Lambda]1H -
	17280*g1^2*Lb*T^4*Y\[Phi]^2*\[Lambda]1H + 133632*T^4*\[Lambda]1H^2 -
	138240*EulerGamma*T^4*\[Lambda]1H^2 - 23040*Lb*m1*T^2*\[Lambda]m + 2880*g2^2*T^4*\[Lambda]m -
	4320*g2^2*Lb*T^4*\[Lambda]m + 23040*EulerGamma*T^4*yt1^2*\[Lambda]m -
	20160*Lb*T^4*yt1^2*\[Lambda]m + 2880*Lf*T^4*yt1^2*\[Lambda]m +
	960*g1^2*T^4*Y\[Phi]^2*\[Lambda]m - 1440*g1^2*Lb*T^4*Y\[Phi]^2*\[Lambda]m -
	23040*EulerGamma*T^4*\[Lambda]1H*\[Lambda]m + 11136*T^4*\[Lambda]m^2 -
	4800*EulerGamma*T^4*\[Lambda]m^2 - 11520*EulerGamma*T^4*\[Lambda]m*\[Lambda]\[Sigma] +
	16704*T^4*\[Lambda]\[Sigma]^2 - 8640*EulerGamma*T^4*\[Lambda]\[Sigma]^2 -
	46080*EulerGamma*T^2*c3^2 + 23040*Lb*T^2*c3^2 - 34560*EulerGamma*T^2*\[Mu]m^2 +
	17280*Lb*T^2*\[Mu]m^2 - 23040*Lb*T^2*\[Lambda]m*\[Mu]\[Sigma] -
	34560*Lb*T^2*\[Lambda]\[Sigma]*\[Mu]\[Sigma] - 69120*Lb*\[Mu]\[Sigma]^2 +
	3828960*g2^4*T^4*Log[2] + 17141760*g3^4*T^4*Log[2] + 2135232*g2^4*nF*T^4*Log[2] +
	1244160*g2^2*g3^2*nF*T^4*Log[2] + 7529472*g3^4*nF*T^4*Log[2] + 562176*g2^4*nF^2*T^4*Log[2] +
	1499136*g3^4*nF^2*T^4*Log[2] + 207360*g1^2*g3^2*nF*T^4*Yd^2*Log[2] +
	19440*g1^4*nF*T^4*Yd^4*Log[2] + 26352*g1^4*nF^2*T^4*Yd^4*Log[2] +
	17568*g1^4*nF^2*T^4*Yd^2*Ye^2*Log[2] + 6480*g1^4*nF*T^4*Ye^4*Log[2] +
	2928*g1^4*nF^2*T^4*Ye^4*Log[2] + 77760*g1^2*g2^2*nF*T^4*Yl^2*Log[2] +
	35136*g1^4*nF^2*T^4*Yd^2*Yl^2*Log[2] + 11712*g1^4*nF^2*T^4*Ye^2*Yl^2*Log[2] +
	12960*g1^4*nF*T^4*Yl^4*Log[2] + 11712*g1^4*nF^2*T^4*Yl^4*Log[2] +
	233280*g1^2*g2^2*nF*T^4*Yq^2*Log[2] + 414720*g1^2*g3^2*nF*T^4*Yq^2*Log[2] +
	105408*g1^4*nF^2*T^4*Yd^2*Yq^2*Log[2] + 35136*g1^4*nF^2*T^4*Ye^2*Yq^2*Log[2] +
	70272*g1^4*nF^2*T^4*Yl^2*Yq^2*Log[2] + 38880*g1^4*nF*T^4*Yq^4*Log[2] +
	105408*g1^4*nF^2*T^4*Yq^4*Log[2] + 324000*g2^2*T^4*yt1^2*Log[2] - 617472*g3^2*T^4*yt1^2*Log[2] +
	67392*g1^2*T^4*Yq^2*yt1^2*Log[2] + 730944*T^4*yt1^4*Log[2] -
	250560*g1^2*T^4*Yq*yt1^2*Yu*Log[2] + 207360*g1^2*g3^2*nF*T^4*Yu^2*Log[2] +
	52704*g1^4*nF^2*T^4*Yd^2*Yu^2*Log[2] + 17568*g1^4*nF^2*T^4*Ye^2*Yu^2*Log[2] +
	35136*g1^4*nF^2*T^4*Yl^2*Yu^2*Log[2] + 105408*g1^4*nF^2*T^4*Yq^2*Yu^2*Log[2] +
	67392*g1^2*T^4*yt1^2*Yu^2*Log[2] + 19440*g1^4*nF*T^4*Yu^4*Log[2] +
	26352*g1^4*nF^2*T^4*Yu^4*Log[2] + 40608*g1^2*T^4*Yq*yt1^2*Y\[Phi]*Log[2] -
	40608*g1^2*T^4*yt1^2*Yu*Y\[Phi]*Log[2] + 285120*g1^2*g2^2*T^4*Y\[Phi]^2*Log[2] +
	42408*g1^4*nF*T^4*Yd^2*Y\[Phi]^2*Log[2] + 14136*g1^4*nF*T^4*Ye^2*Y\[Phi]^2*Log[2] +
	28272*g1^4*nF*T^4*Yl^2*Y\[Phi]^2*Log[2] + 84816*g1^4*nF*T^4*Yq^2*Y\[Phi]^2*Log[2] +
	81216*g1^2*T^4*yt1^2*Y\[Phi]^2*Log[2] + 42408*g1^4*nF*T^4*Yu^2*Y\[Phi]^2*Log[2] +
	81120*g1^4*T^4*Y\[Phi]^4*Log[2] + 1658880*T^4*\[Lambda]1H^2*Log[2] +
	138240*T^4*\[Lambda]1H*\[Lambda]m*Log[2] + 97920*T^4*\[Lambda]m^2*Log[2] +
	69120*T^4*\[Lambda]m*\[Lambda]\[Sigma]*Log[2] + 155520*T^4*\[Lambda]\[Sigma]^2*Log[2] -
	4976640*g2^2*m1*T^2*Log[Glaisher] - 17631360*g2^4*T^4*Log[Glaisher] -
	81838080*g3^4*T^4*Log[Glaisher] - 3432960*g2^4*nF*T^4*Log[Glaisher] -
	13639680*g3^4*nF*T^4*Log[Glaisher] - 368640*g2^4*nF^2*T^4*Log[Glaisher] -
	983040*g3^4*nF^2*T^4*Log[Glaisher] - 17280*g1^4*nF^2*T^4*Yd^4*Log[Glaisher] -
	11520*g1^4*nF^2*T^4*Yd^2*Ye^2*Log[Glaisher] - 1920*g1^4*nF^2*T^4*Ye^4*Log[Glaisher] -
	23040*g1^4*nF^2*T^4*Yd^2*Yl^2*Log[Glaisher] - 7680*g1^4*nF^2*T^4*Ye^2*Yl^2*Log[Glaisher] -
	7680*g1^4*nF^2*T^4*Yl^4*Log[Glaisher] - 69120*g1^4*nF^2*T^4*Yd^2*Yq^2*Log[Glaisher] -
	23040*g1^4*nF^2*T^4*Ye^2*Yq^2*Log[Glaisher] - 46080*g1^4*nF^2*T^4*Yl^2*Yq^2*Log[Glaisher] -
	69120*g1^4*nF^2*T^4*Yq^4*Log[Glaisher] - 933120*g2^2*T^4*yt1^2*Log[Glaisher] +
	5529600*g3^2*T^4*yt1^2*Log[Glaisher] + 69120*g1^2*T^4*Yq^2*yt1^2*Log[Glaisher] -
	622080*T^4*yt1^4*Log[Glaisher] + 898560*g1^2*T^4*Yq*yt1^2*Yu*Log[Glaisher] -
	34560*g1^4*nF^2*T^4*Yd^2*Yu^2*Log[Glaisher] - 11520*g1^4*nF^2*T^4*Ye^2*Yu^2*Log[Glaisher] -
	23040*g1^4*nF^2*T^4*Yl^2*Yu^2*Log[Glaisher] - 69120*g1^4*nF^2*T^4*Yq^2*Yu^2*Log[Glaisher] +
	69120*g1^2*T^4*yt1^2*Yu^2*Log[Glaisher] - 17280*g1^4*nF^2*T^4*Yu^4*Log[Glaisher] -
	207360*g1^2*T^4*Yq*yt1^2*Y\[Phi]*Log[Glaisher] + 207360*g1^2*T^4*yt1^2*Yu*Y\[Phi]*
	Log[Glaisher] - 1658880*g1^2*m1*T^2*Y\[Phi]^2*Log[Glaisher] -
	1762560*g1^2*g2^2*T^4*Y\[Phi]^2*Log[Glaisher] - 2880*g1^4*nF*T^4*Yd^2*Y\[Phi]^2*Log[Glaisher] -
	960*g1^4*nF*T^4*Ye^2*Y\[Phi]^2*Log[Glaisher] - 1920*g1^4*nF*T^4*Yl^2*Y\[Phi]^2*Log[Glaisher] -
	5760*g1^4*nF*T^4*Yq^2*Y\[Phi]^2*Log[Glaisher] - 587520*g1^2*T^4*yt1^2*Y\[Phi]^2*Log[Glaisher] -
	2880*g1^4*nF*T^4*Yu^2*Y\[Phi]^2*Log[Glaisher] - 336000*g1^4*T^4*Y\[Phi]^4*Log[Glaisher] -
	2488320*g2^2*T^4*\[Lambda]1H*Log[Glaisher] - 829440*g1^2*T^4*Y\[Phi]^2*\[Lambda]1H*
	Log[Glaisher] - 3317760*T^4*\[Lambda]1H^2*Log[Glaisher] -
	207360*g2^2*T^4*\[Lambda]m*Log[Glaisher] - 69120*g1^2*T^4*Y\[Phi]^2*\[Lambda]m*Log[Glaisher] -
	276480*T^4*\[Lambda]m^2*Log[Glaisher] - 414720*T^4*\[Lambda]\[Sigma]^2*Log[Glaisher] +
	552960*T^2*c3^2*Log[Glaisher] + 414720*T^2*\[Mu]m^2*Log[Glaisher] +
	1914480*g2^4*T^4*Log[Pi] + 8570880*g3^4*T^4*Log[Pi] + 938880*g2^4*nF*T^4*Log[Pi] +
	622080*g2^2*g3^2*nF*T^4*Log[Pi] + 3363840*g3^4*nF*T^4*Log[Pi] + 184320*g2^4*nF^2*T^4*Log[Pi] +
	491520*g3^4*nF^2*T^4*Log[Pi] + 103680*g1^2*g3^2*nF*T^4*Yd^2*Log[Pi] +
	9720*g1^4*nF*T^4*Yd^4*Log[Pi] + 8640*g1^4*nF^2*T^4*Yd^4*Log[Pi] +
	5760*g1^4*nF^2*T^4*Yd^2*Ye^2*Log[Pi] + 3240*g1^4*nF*T^4*Ye^4*Log[Pi] +
	960*g1^4*nF^2*T^4*Ye^4*Log[Pi] + 38880*g1^2*g2^2*nF*T^4*Yl^2*Log[Pi] +
	11520*g1^4*nF^2*T^4*Yd^2*Yl^2*Log[Pi] + 3840*g1^4*nF^2*T^4*Ye^2*Yl^2*Log[Pi] +
	6480*g1^4*nF*T^4*Yl^4*Log[Pi] + 3840*g1^4*nF^2*T^4*Yl^4*Log[Pi] +
	116640*g1^2*g2^2*nF*T^4*Yq^2*Log[Pi] + 207360*g1^2*g3^2*nF*T^4*Yq^2*Log[Pi] +
	34560*g1^4*nF^2*T^4*Yd^2*Yq^2*Log[Pi] + 11520*g1^4*nF^2*T^4*Ye^2*Yq^2*Log[Pi] +
	23040*g1^4*nF^2*T^4*Yl^2*Yq^2*Log[Pi] + 19440*g1^4*nF*T^4*Yq^4*Log[Pi] +
	34560*g1^4*nF^2*T^4*Yq^4*Log[Pi] + 194400*g2^2*T^4*yt1^2*Log[Pi] +
	138240*g3^2*T^4*yt1^2*Log[Pi] + 51840*g1^2*T^4*Yq^2*yt1^2*Log[Pi] + 311040*T^4*yt1^4*Log[Pi] -
	77760*g1^2*T^4*Yq*yt1^2*Yu*Log[Pi] + 103680*g1^2*g3^2*nF*T^4*Yu^2*Log[Pi] +
	17280*g1^4*nF^2*T^4*Yd^2*Yu^2*Log[Pi] + 5760*g1^4*nF^2*T^4*Ye^2*Yu^2*Log[Pi] +
	11520*g1^4*nF^2*T^4*Yl^2*Yu^2*Log[Pi] + 34560*g1^4*nF^2*T^4*Yq^2*Yu^2*Log[Pi] +
	51840*g1^2*T^4*yt1^2*Yu^2*Log[Pi] + 9720*g1^4*nF*T^4*Yu^4*Log[Pi] +
	8640*g1^4*nF^2*T^4*Yu^4*Log[Pi] + 12960*g1^2*T^4*Yq*yt1^2*Y\[Phi]*Log[Pi] -
	12960*g1^2*T^4*yt1^2*Yu*Y\[Phi]*Log[Pi] + 142560*g1^2*g2^2*T^4*Y\[Phi]^2*Log[Pi] +
	17640*g1^4*nF*T^4*Yd^2*Y\[Phi]^2*Log[Pi] + 5880*g1^4*nF*T^4*Ye^2*Y\[Phi]^2*Log[Pi] +
	11760*g1^4*nF*T^4*Yl^2*Y\[Phi]^2*Log[Pi] + 35280*g1^4*nF*T^4*Yq^2*Y\[Phi]^2*Log[Pi] +
	25920*g1^2*T^4*yt1^2*Y\[Phi]^2*Log[Pi] + 17640*g1^4*nF*T^4*Yu^2*Y\[Phi]^2*Log[Pi] +
	40560*g1^4*T^4*Y\[Phi]^4*Log[Pi] + 414720*T^4*\[Lambda]1H^2*Log[Pi] +
	34560*T^4*\[Lambda]m^2*Log[Pi] + 51840*T^4*\[Lambda]\[Sigma]^2*Log[Pi] +
	414720*T^4*\[Lambda]1H^2*Log[Pi*T] + 69120*T^4*\[Lambda]1H*\[Lambda]m*Log[Pi*T] +
	14400*T^4*\[Lambda]m^2*Log[Pi*T] + 34560*T^4*\[Lambda]m*\[Lambda]\[Sigma]*Log[Pi*T] +
	25920*T^4*\[Lambda]\[Sigma]^2*Log[Pi*T] - 25920*g2^4*T^4*Log[4*Pi*T] +
	138240*T^4*yt1^4*Log[4*Pi*T] - 17280*g1^2*g2^2*T^4*Y\[Phi]^2*Log[4*Pi*T] -
	8640*g1^4*T^4*Y\[Phi]^4*Log[4*Pi*T] + 207360*g2^2*T^4*\[Lambda]1H*Log[4*Pi*T] -
	276480*T^4*yt1^2*\[Lambda]1H*Log[4*Pi*T] + 69120*g1^2*T^4*Y\[Phi]^2*\[Lambda]1H*Log[4*Pi*T] -
	552960*T^4*\[Lambda]1H^2*Log[4*Pi*T] + 17280*g2^2*T^4*\[Lambda]m*Log[4*Pi*T] -
	23040*T^4*yt1^2*\[Lambda]m*Log[4*Pi*T] + 5760*g1^2*T^4*Y\[Phi]^2*\[Lambda]m*Log[4*Pi*T] -
	46080*T^4*\[Lambda]1H*\[Lambda]m*Log[4*Pi*T] - 32640*T^4*\[Lambda]m^2*Log[4*Pi*T] -
	23040*T^4*\[Lambda]m*\[Lambda]\[Sigma]*Log[4*Pi*T] - 51840*T^4*\[Lambda]\[Sigma]^2*Log[4*Pi*T] -
	1914480*g2^4*T^4*Log[mu4/T] - 8570880*g3^4*T^4*Log[mu4/T] -
	938880*g2^4*nF*T^4*Log[mu4/T] - 622080*g2^2*g3^2*nF*T^4*Log[mu4/T] -
	3363840*g3^4*nF*T^4*Log[mu4/T] - 184320*g2^4*nF^2*T^4*Log[mu4/T] -
	491520*g3^4*nF^2*T^4*Log[mu4/T] - 103680*g1^2*g3^2*nF*T^4*Yd^2*Log[mu4/T] -
	9720*g1^4*nF*T^4*Yd^4*Log[mu4/T] - 8640*g1^4*nF^2*T^4*Yd^4*Log[mu4/T] -
	5760*g1^4*nF^2*T^4*Yd^2*Ye^2*Log[mu4/T] - 3240*g1^4*nF*T^4*Ye^4*Log[mu4/T] -
	960*g1^4*nF^2*T^4*Ye^4*Log[mu4/T] - 38880*g1^2*g2^2*nF*T^4*Yl^2*Log[mu4/T] -
	11520*g1^4*nF^2*T^4*Yd^2*Yl^2*Log[mu4/T] - 3840*g1^4*nF^2*T^4*Ye^2*Yl^2*Log[mu4/T] -
	6480*g1^4*nF*T^4*Yl^4*Log[mu4/T] - 3840*g1^4*nF^2*T^4*Yl^4*Log[mu4/T] -
	116640*g1^2*g2^2*nF*T^4*Yq^2*Log[mu4/T] - 207360*g1^2*g3^2*nF*T^4*Yq^2*Log[mu4/T] -
	34560*g1^4*nF^2*T^4*Yd^2*Yq^2*Log[mu4/T] - 11520*g1^4*nF^2*T^4*Ye^2*Yq^2*Log[mu4/T] -
	23040*g1^4*nF^2*T^4*Yl^2*Yq^2*Log[mu4/T] - 19440*g1^4*nF*T^4*Yq^4*Log[mu4/T] -
	34560*g1^4*nF^2*T^4*Yq^4*Log[mu4/T] - 194400*g2^2*T^4*yt1^2*Log[mu4/T] -
	138240*g3^2*T^4*yt1^2*Log[mu4/T] - 51840*g1^2*T^4*Yq^2*yt1^2*Log[mu4/T] -
	311040*T^4*yt1^4*Log[mu4/T] + 77760*g1^2*T^4*Yq*yt1^2*Yu*Log[mu4/T] -
	103680*g1^2*g3^2*nF*T^4*Yu^2*Log[mu4/T] - 17280*g1^4*nF^2*T^4*Yd^2*Yu^2*Log[mu4/T] -
	5760*g1^4*nF^2*T^4*Ye^2*Yu^2*Log[mu4/T] - 11520*g1^4*nF^2*T^4*Yl^2*Yu^2*Log[mu4/T] -
	34560*g1^4*nF^2*T^4*Yq^2*Yu^2*Log[mu4/T] - 51840*g1^2*T^4*yt1^2*Yu^2*Log[mu4/T] -
	9720*g1^4*nF*T^4*Yu^4*Log[mu4/T] - 8640*g1^4*nF^2*T^4*Yu^4*Log[mu4/T] -
	12960*g1^2*T^4*Yq*yt1^2*Y\[Phi]*Log[mu4/T] + 12960*g1^2*T^4*yt1^2*Yu*Y\[Phi]*Log[mu4/T] -
	142560*g1^2*g2^2*T^4*Y\[Phi]^2*Log[mu4/T] - 17640*g1^4*nF*T^4*Yd^2*Y\[Phi]^2*Log[mu4/T] -
	5880*g1^4*nF*T^4*Ye^2*Y\[Phi]^2*Log[mu4/T] - 11760*g1^4*nF*T^4*Yl^2*Y\[Phi]^2*Log[mu4/T] -
	35280*g1^4*nF*T^4*Yq^2*Y\[Phi]^2*Log[mu4/T] - 25920*g1^2*T^4*yt1^2*Y\[Phi]^2*Log[mu4/T] -
	17640*g1^4*nF*T^4*Yu^2*Y\[Phi]^2*Log[mu4/T] - 40560*g1^4*T^4*Y\[Phi]^4*Log[mu4/T] -
	414720*T^4*\[Lambda]1H^2*Log[mu4/T] - 34560*T^4*\[Lambda]m^2*Log[mu4/T] -
	51840*T^4*\[Lambda]\[Sigma]^2*Log[mu4/T] + 12960*g2^4*T^4*Log[mu4^2] -
	69120*T^4*yt1^4*Log[mu4^2] + 8640*g1^2*g2^2*T^4*Y\[Phi]^2*Log[mu4^2] +
	4320*g1^4*T^4*Y\[Phi]^4*Log[mu4^2] - 103680*g2^2*T^4*\[Lambda]1H*Log[mu4^2] +
	138240*T^4*yt1^2*\[Lambda]1H*Log[mu4^2] - 34560*g1^2*T^4*Y\[Phi]^2*\[Lambda]1H*Log[mu4^2] +
	69120*T^4*\[Lambda]1H^2*Log[mu4^2] - 8640*g2^2*T^4*\[Lambda]m*Log[mu4^2] +
	11520*T^4*yt1^2*\[Lambda]m*Log[mu4^2] - 2880*g1^2*T^4*Y\[Phi]^2*\[Lambda]m*Log[mu4^2] -
	11520*T^4*\[Lambda]1H*\[Lambda]m*Log[mu4^2] + 9120*T^4*\[Lambda]m^2*Log[mu4^2] -
	5760*T^4*\[Lambda]m*\[Lambda]\[Sigma]*Log[mu4^2] + 12960*T^4*\[Lambda]\[Sigma]^2*
	Log[mu4^2] - 829440*g2^2*m1*T^2*Log[mu3/mu4] - 1537920*g2^4*T^4*Log[mu3/mu4] -
	6635520*g3^4*T^4*Log[mu3/mu4] - 552960*g2^4*nF*T^4*Log[mu3/mu4] -
	2211840*g3^4*nF*T^4*Log[mu3/mu4] - 207360*c61*g2^2*T^6*Log[mu3/mu4] -
	34560*c63*g2^2*T^6*Log[mu3/mu4] - 17280*c64*g2^2*T^6*Log[mu3/mu4] -
	207360*g2^2*T^4*yt1^2*Log[mu3/mu4] - 276480*g1^2*m1*T^2*Y\[Phi]^2*Log[mu3/mu4] -
	103680*g1^2*g2^2*T^4*Y\[Phi]^2*Log[mu3/mu4] - 69120*c61*g1^2*T^6*Y\[Phi]^2*
	Log[mu3/mu4] - 11520*c63*g1^2*T^6*Y\[Phi]^2*Log[mu3/mu4] -
	5760*c64*g1^2*T^6*Y\[Phi]^2*Log[mu3/mu4] - 69120*g1^2*T^4*yt1^2*Y\[Phi]^2*
	Log[mu3/mu4] - 17280*g1^4*T^4*Y\[Phi]^4*Log[mu3/mu4] -
	414720*g2^2*T^4*\[Lambda]1H*Log[mu3/mu4] - 138240*g1^2*T^4*Y\[Phi]^2*\[Lambda]1H*
	Log[mu3/mu4] - 34560*g2^2*T^4*\[Lambda]m*Log[mu3/mu4] -
	11520*g1^2*T^4*Y\[Phi]^2*\[Lambda]m*Log[mu3/mu4] + 92160*T^2*c3^2*Log[mu3/mu4] +
	69120*T^2*\[Mu]m^2*Log[mu3/mu4] + 42019200*g2^4*T^4*Derivative[1][Zeta][-3] +
	210124800*g3^4*T^4*Derivative[1][Zeta][-3] + 576000*g2^4*nF*T^4*Derivative[1][Zeta][-3] +
	1843200*g3^4*nF*T^4*Derivative[1][Zeta][-3] + 1843200*g2^4*nF^2*T^4*Derivative[1][Zeta][-3] +
	4915200*g3^4*nF^2*T^4*Derivative[1][Zeta][-3] + 86400*g1^4*nF^2*T^4*Yd^4*
	Derivative[1][Zeta][-3] + 57600*g1^4*nF^2*T^4*Yd^2*Ye^2*Derivative[1][Zeta][-3] +
	9600*g1^4*nF^2*T^4*Ye^4*Derivative[1][Zeta][-3] + 115200*g1^4*nF^2*T^4*Yd^2*Yl^2*
	Derivative[1][Zeta][-3] + 38400*g1^4*nF^2*T^4*Ye^2*Yl^2*Derivative[1][Zeta][-3] +
	38400*g1^4*nF^2*T^4*Yl^4*Derivative[1][Zeta][-3] + 345600*g1^4*nF^2*T^4*Yd^2*Yq^2*
	Derivative[1][Zeta][-3] + 115200*g1^4*nF^2*T^4*Ye^2*Yq^2*Derivative[1][Zeta][-3] +
	230400*g1^4*nF^2*T^4*Yl^2*Yq^2*Derivative[1][Zeta][-3] + 345600*g1^4*nF^2*T^4*Yq^4*
	Derivative[1][Zeta][-3] - 1555200*g2^2*T^4*yt1^2*Derivative[1][Zeta][-3] -
	27648000*g3^2*T^4*yt1^2*Derivative[1][Zeta][-3] - 1036800*g1^2*T^4*Yq^2*yt1^2*
	Derivative[1][Zeta][-3] + 3110400*T^4*yt1^4*Derivative[1][Zeta][-3] -
	3110400*g1^2*T^4*Yq*yt1^2*Yu*Derivative[1][Zeta][-3] + 172800*g1^4*nF^2*T^4*Yd^2*Yu^2*
	Derivative[1][Zeta][-3] + 57600*g1^4*nF^2*T^4*Ye^2*Yu^2*Derivative[1][Zeta][-3] +
	115200*g1^4*nF^2*T^4*Yl^2*Yu^2*Derivative[1][Zeta][-3] + 345600*g1^4*nF^2*T^4*Yq^2*Yu^2*
	Derivative[1][Zeta][-3] - 1036800*g1^2*T^4*yt1^2*Yu^2*Derivative[1][Zeta][-3] +
	86400*g1^4*nF^2*T^4*Yu^4*Derivative[1][Zeta][-3] + 518400*g1^2*T^4*Yq*yt1^2*Y\[Phi]*
	Derivative[1][Zeta][-3] - 518400*g1^2*T^4*yt1^2*Yu*Y\[Phi]*Derivative[1][Zeta][-3] +
	5702400*g1^2*g2^2*T^4*Y\[Phi]^2*Derivative[1][Zeta][-3] + 14400*g1^4*nF*T^4*Yd^2*Y\[Phi]^2*
	Derivative[1][Zeta][-3] + 4800*g1^4*nF*T^4*Ye^2*Y\[Phi]^2*Derivative[1][Zeta][-3] +
	9600*g1^4*nF*T^4*Yl^2*Y\[Phi]^2*Derivative[1][Zeta][-3] + 28800*g1^4*nF*T^4*Yq^2*Y\[Phi]^2*
	Derivative[1][Zeta][-3] + 1036800*g1^2*T^4*yt1^2*Y\[Phi]^2*Derivative[1][Zeta][-3] +
	14400*g1^4*nF*T^4*Yu^2*Y\[Phi]^2*Derivative[1][Zeta][-3] +
	1161600*g1^4*T^4*Y\[Phi]^4*Derivative[1][Zeta][-3] + 16588800*T^4*\[Lambda]1H^2*
	Derivative[1][Zeta][-3] + 1382400*T^4*\[Lambda]m^2*Derivative[1][Zeta][-3] +
	2073600*T^4*\[Lambda]\[Sigma]^2*Derivative[1][Zeta][-3])/Pi^2
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
	{\[Lambda]1H3dUS -> \[Lambda]1H3d - ((8*lamVL[1]^2)/Sqrt[musqSU3] + lamVL[2]^2/Sqrt[musqU1] + (3*lamVL[3]^2)/Sqrt[musqSU2] + (4*lamVL[4]^2)/(Sqrt[musqSU2] + Sqrt[musqU1]))/(32*Pi), \[Lambda]m3dUS -> -1/192*(-192*Pi*\[Lambda]m3d + (36*lamVL[3]*lamVL[5])/Sqrt[musqSU2] + (3*\[Lambda]m3d*lamVVSL[1]^2)/musqSU2^(3/2) + (12*musqU1*lamVL[2]*lamVL[6] + \[Lambda]m3d*lamVVSL[2]^2)/musqU1^(3/2))/Pi, \[Lambda]\[Sigma]3dUS -> -1/96*(-96*Pi*\[Lambda]\[Sigma]3d + (9*lamVL[5]^2)/Sqrt[musqSU2] + (3*\[Lambda]\[Sigma]3d*lamVVSL[1]^2)/musqSU2^(3/2) + (3*musqU1*lamVL[6]^2 + \[Lambda]\[Sigma]3d*lamVVSL[2]^2)/musqU1^(3/2))/Pi, g13dUS^2 -> g13d^2, g23dUS^2 -> g23d^2 - g23d^4/(24*Pi*Sqrt[musqSU2]), g33dUS^2 -> g33d^2 - g33d^4/(16*Pi*Sqrt[musqSU3]), c33dUS -> c33d + (-36*musqSU2*musqU1^(3/2)*lamVL[5]*lamVVSL[1] + 3*musqU1^(3/2)*lamVVSL[1]^2*(-c33d + 2*lamVVSL[1]) - musqSU2^(3/2)*lamVVSL[2]*(12*musqU1*lamVL[6] + (c33d - 2*lamVVSL[2])*lamVVSL[2]))/(128*Pi*musqSU2^(3/2)*musqU1^(3/2)), \[Mu]m3dUS -> \[Mu]m3d - ((3*lamVVSL[1]*(48*musqSU2*lamVL[3] + \[Mu]m3d*lamVVSL[1]))/musqSU2^(3/2) + (48*lamVL[2]*lamVVSL[2])/Sqrt[musqU1] + (\[Mu]m3d*lamVVSL[2]^2)/musqU1^(3/2))/(384*Pi)}
]];


AppendTo[testList,
TestCreate[PrintScalarMassUS["LO"],
	{m13dUS->m13d-(8 Sqrt[musqSU3] lamVL[1]+Sqrt[musqU1] lamVL[2]+3 Sqrt[musqSU2] lamVL[3])/(8 \[Pi]),\[Mu]\[Sigma]3dUS->\[Mu]\[Sigma]3d-(6 Sqrt[musqSU2] lamVL[5]+(3 lamVVSL[1]^2)/Sqrt[musqSU2]+(2 musqU1 lamVL[6]+lamVVSL[2]^2)/Sqrt[musqU1])/(16 \[Pi])}
]];
AppendTo[testList,
TestCreate[PrintScalarMassUS["NLO"],
	{m13dUS -> (48*g33d^2*lamVL[1] - 16*lamVL[1]^2 - 2*lamVL[2]^2 + 12*g23d^2*lamVL[3] - 6*lamVL[3]^2 - 12*lamVL[4]^2 + (24*Sqrt[musqSU2]*lamVL[1]*lamVLL[1])/Sqrt[musqSU3] + (24*Sqrt[musqSU3]*lamVL[3]*lamVLL[1])/Sqrt[musqSU2] + (80*lamVL[1]*lamVLL[2])/3 + 5*lamVL[3]*lamVLL[3] + (8*Sqrt[musqU1]*lamVL[1]*lamVLL[7])/Sqrt[musqSU3] + (8*Sqrt[musqSU3]*lamVL[2]*lamVLL[7])/Sqrt[musqU1] + (3*Sqrt[musqSU2]*lamVL[2]*lamVLL[8])/Sqrt[musqU1] + (3*Sqrt[musqU1]*lamVL[3]*lamVLL[8])/Sqrt[musqSU2] + lamVL[2]*lamVLL[9] - 6*(g23d^4 - 8*g23d^2*lamVL[3] + 2*lamVL[3]^2)*Log[mu3/(2*Sqrt[musqSU2])] + 32*(6*g33d^2 - lamVL[1])*lamVL[1]*Log[mu3/(2*Sqrt[musqSU3])] - 24*lamVL[4]^2*Log[mu3/(Sqrt[musqSU2] + Sqrt[musqU1])] - 4*lamVL[2]^2*Log[mu3/(2*Sqrt[musqU1])])/(128*Pi^2), \[Mu]\[Sigma]3dUS -> ((3*lamVL[5]*(24*Sqrt[musqSU3]*lamVLL[1] + 5*Sqrt[musqSU2]*lamVLL[3] + 3*Sqrt[musqU1]*lamVLL[8]))/Sqrt[musqSU2] + (3*lamVL[6]*(8*Sqrt[musqSU3]*lamVLL[7] + 3*Sqrt[musqSU2]*lamVLL[8] + Sqrt[musqU1]*lamVLL[9]))/Sqrt[musqU1] + ((3*musqU1^(3/2)*lamVVSL[1]^2 + musqSU2^(3/2)*lamVVSL[2]^2)*(6*musqSU2*Sqrt[musqU1]*lamVL[5] + 3*Sqrt[musqU1]*lamVVSL[1]^2 + Sqrt[musqSU2]*(-16*Sqrt[musqU1]*Pi*\[Mu]\[Sigma]3d + 2*musqU1*lamVL[6] + lamVVSL[2]^2)))/(8*musqSU2^2*musqU1^2) - 18*lamVL[5]^2*(1 + 2*Log[mu3/(2*Sqrt[musqSU2])]) + 36*g23d^2*lamVL[5]*(1 + 4*Log[mu3/(2*Sqrt[musqSU2])]) - 6*lamVL[6]^2*(1 + 2*Log[mu3/(2*Sqrt[musqU1])]))/(384*Pi^2)}
]];


AppendTo[testList,
TestCreate[BetaFunctions3DUS[],
	{m13dUS->1/(256 \[Pi]^2) (-51 g23dUS^4+5 g13dUS^4 Y\[Phi]^4+18 g23dUS^2 (g13dUS^2 Y\[Phi]^2-8 \[Lambda]1H3dUS)-48 g13dUS^2 Y\[Phi]^2 \[Lambda]1H3dUS+8 (24 \[Lambda]1H3dUS^2+\[Lambda]m3dUS^2)),\[Mu]\[Sigma]3dUS->(-((3 g23dUS^2+g13dUS^2 Y\[Phi]^2-2 \[Lambda]m3dUS) \[Lambda]m3dUS)+6 \[Lambda]\[Sigma]3dUS^2)/(16 \[Pi]^2),\[Mu]13dUS->(4 \[Lambda]\[Sigma]3dUS c33dUS-3 g23dUS^2 \[Mu]m3dUS-g13dUS^2 Y\[Phi]^2 \[Mu]m3dUS+2 \[Lambda]m3dUS \[Mu]m3dUS)/(32 \[Pi]^2)}
]];


AppendTo[testList,
TestCreate[PrintTadpolesUS["LO"],
	{\[Mu]13dUS->\[Mu]13d-(3 Sqrt[musqSU2] lamVVSL[1]+Sqrt[musqU1] lamVVSL[2])/(8 \[Pi])-(\[Mu]13d ((3 lamVVSL[1]^2)/musqSU2^(3/2)+lamVVSL[2]^2/musqU1^(3/2)))/(384 \[Pi])}
]];


report=TestReport[testList]
report["ResultsDataset"]


(* ::Subsection:: *)
(*Test scalar 2 softer*)


PerformDRsoft[{5}];


testList={};


AppendTo[testList,
TestCreate[PrintCouplingsUS[],
	{\[Lambda]1H3dUS -> (384*\[Lambda]1H3d + (9*\[Mu]m3d^4)/(Pi*\[Mu]\[Sigma]3d^(5/2)) + (48*\[Lambda]m3d*\[Mu]m3d*(4*Pi*\[Mu]13d - c33d*Sqrt[\[Mu]\[Sigma]3d]))/(Pi*\[Mu]\[Sigma]3d^2) + (12*(-12*\[Lambda]1H3d + \[Lambda]m3d)*\[Mu]m3d^2)/(Pi*\[Mu]\[Sigma]3d^(3/2)) + (2*\[Mu]m3d*(5*c33d*\[Mu]13d*\[Mu]m3d^3 - 8*\[Mu]13d*\[Mu]m3d*(3*c33d*\[Lambda]1H3d + \[Lambda]m3d*\[Mu]m3d)*\[Mu]\[Sigma]3d + (32*\[Lambda]1H3d*\[Lambda]m3d*\[Mu]13d + \[Mu]m3d^3)*\[Mu]\[Sigma]3d^2 - 8*\[Lambda]1H3d*\[Mu]m3d*\[Mu]\[Sigma]3d^3))/(Pi*\[Mu]\[Sigma]3d^(9/2)) - (12*(\[Lambda]m3d^2/Sqrt[\[Mu]\[Sigma]3d] + (8*lamVL[1]^2)/Sqrt[musqSU3] + lamVL[2]^2/Sqrt[musqU1] + (3*lamVL[3]^2)/Sqrt[musqSU2] + (4*lamVL[4]^2)/(Sqrt[musqSU2] + Sqrt[musqU1])))/Pi - (6*\[Mu]m3d*(-2*c33d^2*\[Mu]m3d*Sqrt[\[Mu]\[Sigma]3d] - \[Mu]m3d^3*Sqrt[\[Mu]\[Sigma]3d] + 2*c33d*(8*Pi*\[Mu]13d*\[Mu]m3d + Sqrt[\[Mu]\[Sigma]3d]*(\[Mu]m3d^2 - 2*\[Lambda]m3d*\[Mu]\[Sigma]3d)) + \[Mu]m3d*\[Mu]\[Sigma]3d*(2*(-4*\[Lambda]m3d + 3*\[Lambda]\[Sigma]3d)*Sqrt[\[Mu]\[Sigma]3d] + 8*Pi*\[Mu]\[Sigma]3d + 3*Sqrt[musqSU2]*lamVL[5] + Sqrt[musqU1]*lamVL[6])))/(Pi*\[Mu]\[Sigma]3d^3) + (3*\[Mu]13d*(20*c33d^3*Sqrt[musqSU2]*Sqrt[musqU1]*\[Mu]m3d^2 - 20*c33d^2*Sqrt[musqSU2]*Sqrt[musqU1]*\[Mu]m3d*(\[Mu]m3d^2 + 2*\[Lambda]m3d*\[Mu]\[Sigma]3d) + 4*\[Mu]m3d*\[Mu]\[Sigma]3d*(Sqrt[musqU1]*(Sqrt[musqSU2]*((-10*\[Lambda]m3d + 3*\[Lambda]\[Sigma]3d)*\[Mu]m3d^2 + 2*\[Lambda]m3d*(24*\[Lambda]1H3d - 10*\[Lambda]m3d + 9*\[Lambda]\[Sigma]3d)*\[Mu]\[Sigma]3d) - 3*Sqrt[\[Mu]\[Sigma]3d]*(-4*musqSU2*\[Lambda]m3d + \[Mu]\[Sigma]3d*lamVL[3])*lamVL[5]) + Sqrt[musqSU2]*Sqrt[\[Mu]\[Sigma]3d]*(4*musqU1*\[Lambda]m3d - \[Mu]\[Sigma]3d*lamVL[2])*lamVL[6]) + c33d*Sqrt[musqSU2]*Sqrt[musqU1]*(25*\[Mu]m3d^4 + 12*\[Lambda]m3d^2*\[Mu]\[Sigma]3d^2 - 12*\[Mu]m3d^2*(12*\[Lambda]1H3d*\[Mu]\[Sigma]3d - 7*\[Lambda]m3d*\[Mu]\[Sigma]3d + 5*\[Lambda]\[Sigma]3d*\[Mu]\[Sigma]3d + 3*Sqrt[musqSU2]*Sqrt[\[Mu]\[Sigma]3d]*lamVL[5] + Sqrt[musqU1]*Sqrt[\[Mu]\[Sigma]3d]*lamVL[6]))))/(Sqrt[musqSU2]*Sqrt[musqU1]*Pi*\[Mu]\[Sigma]3d^(9/2)))/384, g13dUS^2 -> g13d^2, g23dUS^2 -> g23d^2 - g23d^4/(24*Sqrt[musqSU2]*Pi), g33dUS^2 -> g33d^2 - g33d^4/(16*Sqrt[musqSU3]*Pi)}
]];


AppendTo[testList,
TestCreate[PrintScalarMassUS["LO"],
	{m13dUS -> m13d - (-2*c33d^2*\[Mu]13d*\[Mu]m3d + c33d*(-2*\[Mu]m3d*\[Mu]\[Sigma]3d^2 + \[Mu]13d*(\[Mu]m3d^2 + 2*\[Lambda]m3d*\[Mu]\[Sigma]3d)) + Sqrt[\[Mu]\[Sigma]3d]*(6*\[Lambda]\[Sigma]3d*\[Mu]13d*\[Mu]m3d*Sqrt[\[Mu]\[Sigma]3d] + 8*Pi*\[Mu]13d*\[Mu]m3d*\[Mu]\[Sigma]3d + \[Mu]m3d^2*\[Mu]\[Sigma]3d^(3/2) + 2*\[Lambda]m3d*Sqrt[\[Mu]\[Sigma]3d]*(-2*\[Mu]13d*\[Mu]m3d + \[Mu]\[Sigma]3d^2) + 16*Sqrt[musqSU3]*\[Mu]\[Sigma]3d^2*lamVL[1] + 2*Sqrt[musqU1]*\[Mu]\[Sigma]3d^2*lamVL[2] + 6*Sqrt[musqSU2]*\[Mu]\[Sigma]3d^2*lamVL[3] + 3*Sqrt[musqSU2]*\[Mu]13d*\[Mu]m3d*lamVL[5] + Sqrt[musqU1]*\[Mu]13d*\[Mu]m3d*lamVL[6]))/(16*Pi*\[Mu]\[Sigma]3d^(5/2))}
]];
AppendTo[testList,
TestCreate[PrintScalarMassUS["NLO"],
	{m13dUS -> ((-4*\[Mu]m3d^3*(-4*Pi*\[Mu]13d + c33d*Sqrt[\[Mu]\[Sigma]3d]))/\[Mu]\[Sigma]3d^(5/2) + (2*\[Mu]m3d^4)/\[Mu]\[Sigma]3d^2 - (32*m13d*Pi*\[Mu]m3d^2)/\[Mu]\[Sigma]3d^(3/2) + (16*c33d^2*\[Lambda]m3d)/\[Mu]\[Sigma]3d + (24*\[Lambda]m3d*\[Mu]m3d^2)/\[Mu]\[Sigma]3d + (12*(-6*\[Lambda]1H3d + \[Lambda]m3d)*\[Mu]m3d^2)/\[Mu]\[Sigma]3d + (24*c33d^2*(\[Mu]m3d^2 - 2*\[Lambda]m3d*\[Mu]\[Sigma]3d))/\[Mu]\[Sigma]3d^2 + (4*\[Mu]m3d^2*(\[Lambda]m3d*Sqrt[\[Mu]\[Sigma]3d] + 8*Sqrt[musqSU3]*lamVL[1] + Sqrt[musqU1]*lamVL[2] + 3*Sqrt[musqSU2]*lamVL[3]))/\[Mu]\[Sigma]3d^(3/2) + (12*\[Mu]m3d^2*(4*Pi*\[Mu]13d*\[Mu]m3d - c33d*\[Mu]m3d*Sqrt[\[Mu]\[Sigma]3d] - 8*m13d*Pi*\[Mu]\[Sigma]3d + \[Lambda]m3d*\[Mu]\[Sigma]3d^(3/2) + 8*Sqrt[musqSU3]*\[Mu]\[Sigma]3d*lamVL[1] + Sqrt[musqU1]*\[Mu]\[Sigma]3d*lamVL[2] + 3*Sqrt[musqSU2]*\[Mu]\[Sigma]3d*lamVL[3]))/\[Mu]\[Sigma]3d^(5/2) - (6*\[Mu]m3d^2*(6*\[Lambda]\[Sigma]3d*Sqrt[\[Mu]\[Sigma]3d] + 3*Sqrt[musqSU2]*lamVL[5] + Sqrt[musqU1]*lamVL[6]))/\[Mu]\[Sigma]3d^(3/2) + (12*\[Lambda]m3d*(6*\[Lambda]\[Sigma]3d*Sqrt[\[Mu]\[Sigma]3d] + 3*Sqrt[musqSU2]*lamVL[5] + Sqrt[musqU1]*lamVL[6]))/Sqrt[\[Mu]\[Sigma]3d] + (2*\[Mu]13d*\[Mu]m3d^2*(-2*c33d^2*\[Mu]m3d - 4*\[Lambda]m3d*\[Mu]m3d*\[Mu]\[Sigma]3d + c33d*(\[Mu]m3d^2 + 2*\[Lambda]m3d*\[Mu]\[Sigma]3d) + \[Mu]m3d*Sqrt[\[Mu]\[Sigma]3d]*(6*\[Lambda]\[Sigma]3d*Sqrt[\[Mu]\[Sigma]3d] + 3*Sqrt[musqSU2]*lamVL[5] + Sqrt[musqU1]*lamVL[6])))/\[Mu]\[Sigma]3d^4 + (32*lamVL[1]*(9*Sqrt[musqSU2]*lamVLL[1] + 10*Sqrt[musqSU3]*lamVLL[2] + 3*Sqrt[musqU1]*lamVLL[7]))/Sqrt[musqSU3] + (36*lamVL[3]*(Sqrt[\[Mu]\[Sigma]3d]*lamVL[5] + 8*Sqrt[musqSU3]*lamVLL[1] + (5*Sqrt[musqSU2]*lamVLL[3])/3 + Sqrt[musqU1]*lamVLL[8]))/Sqrt[musqSU2] + (12*lamVL[2]*(Sqrt[\[Mu]\[Sigma]3d]*lamVL[6] + 8*Sqrt[musqSU3]*lamVLL[7] + 3*Sqrt[musqSU2]*lamVLL[8] + Sqrt[musqU1]*lamVLL[9]))/Sqrt[musqU1] - 72*g23d^4*Log[mu3/(2*Sqrt[musqSU2])] - 72*lamVL[3]^2*(1 + 2*Log[mu3/(2*Sqrt[musqSU2])]) + 144*g23d^2*lamVL[3]*(1 + 4*Log[mu3/(2*Sqrt[musqSU2])]) - 192*lamVL[1]^2*(1 + 2*Log[mu3/(2*Sqrt[musqSU3])]) + 576*g33d^2*lamVL[1]*(1 + 4*Log[mu3/(2*Sqrt[musqSU3])]) - 144*lamVL[4]^2*(1 + 2*Log[mu3/(Sqrt[musqSU2] + Sqrt[musqU1])]) - 24*lamVL[2]^2*(1 + 2*Log[mu3/(2*Sqrt[musqU1])]) - 24*\[Lambda]m3d^2*(1 + 2*Log[mu3/(2*Sqrt[\[Mu]\[Sigma]3d])]) + (3*\[Mu]m3d^3*(\[Mu]m3d - 8*c33d*Log[4/3] + \[Mu]m3d*Log[4] + 2*\[Mu]m3d*Log[mu3/Sqrt[\[Mu]\[Sigma]3d]]))/\[Mu]\[Sigma]3d^2 + (48*\[Mu]m3d*(-3*\[Lambda]1H3d*\[Mu]m3d + \[Lambda]m3d*\[Mu]m3d*Log[2] + c33d*\[Lambda]m3d*Log[9/4] - 6*\[Lambda]1H3d*\[Mu]m3d*Log[mu3/Sqrt[\[Mu]\[Sigma]3d]]))/\[Mu]\[Sigma]3d)/(1536*Pi^2)}
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
	{\[Lambda]\[Sigma]3dUS -> -1/1536*((192*m13d^2*\[Lambda]m3d^2 + 16*m13d*(-3*\[Lambda]m3d + \[Lambda]\[Sigma]3d)*\[Mu]m3d^2 + 3*\[Mu]m3d^4)/m13d^(5/2) + 48*(-32*Pi*\[Lambda]\[Sigma]3d + (3*lamVL[5]^2)/Sqrt[musqSU2] + lamVL[6]^2/Sqrt[musqU1]))/Pi, g23dUS^2 -> g23d^2 - (g23d^4*(1/Sqrt[m13d] + 2/Sqrt[musqSU2]))/(48*Pi), g33dUS^2 -> g33d^2 - g33d^4/(16*Pi*Sqrt[musqSU3]), c33dUS -> c33d + (\[Mu]m3d*(-24*m13d*\[Lambda]m3d + \[Mu]m3d*(-c33d + \[Mu]m3d)))/(128*m13d^(3/2)*Pi)}
]];


AppendTo[testList,
TestCreate[PrintScalarMassUS["LO"],
	{\[Mu]\[Sigma]3dUS->\[Mu]\[Sigma]3d-(8 m13d \[Lambda]m3d+\[Mu]m3d^2+2 Sqrt[m13d] (3 Sqrt[musqSU2] lamVL[5]+Sqrt[musqU1] lamVL[6]))/(16 Sqrt[m13d] \[Pi])}
]];
AppendTo[testList,
TestCreate[PrintScalarMassUS["NLO"],
	{\[Mu]\[Sigma]3dUS -> ((288*\[Lambda]1H3d*\[Mu]m3d^2)/m13d + (24*(2*\[Lambda]m3d - 3*\[Lambda]\[Sigma]3d)*\[Mu]m3d^2)/m13d + (3*(8*c33d - \[Mu]m3d)*\[Mu]m3d^3)/m13d^2 + (2*\[Mu]m3d^4)/m13d^2 + (192*\[Lambda]m3d*\[Mu]m3d*(-c33d + \[Mu]m3d))/m13d - (32*Pi*\[Mu]m3d^2*\[Mu]\[Sigma]3d)/m13d^(3/2) + (192*\[Lambda]m3d*(12*Sqrt[m13d]*\[Lambda]1H3d + 8*Sqrt[musqSU3]*lamVL[1] + Sqrt[musqU1]*lamVL[2] + 3*Sqrt[musqSU2]*lamVL[3]))/Sqrt[m13d] - (24*\[Mu]m3d^2*(12*Sqrt[m13d]*\[Lambda]1H3d + 8*Sqrt[musqSU3]*lamVL[1] + Sqrt[musqU1]*lamVL[2] + 3*Sqrt[musqSU2]*lamVL[3]))/m13d^(3/2) + (4*\[Mu]m3d^2*(4*Sqrt[m13d]*\[Lambda]m3d + 3*Sqrt[musqSU2]*lamVL[5] + Sqrt[musqU1]*lamVL[6]))/m13d^(3/2) + (48*lamVL[5]*(12*Sqrt[m13d]*lamVL[3] + 24*Sqrt[musqSU3]*lamVLL[1] + 5*Sqrt[musqSU2]*lamVLL[3] + 3*Sqrt[musqU1]*lamVLL[8]))/Sqrt[musqSU2] + (48*lamVL[6]*(4*Sqrt[m13d]*lamVL[2] + 8*Sqrt[musqSU3]*lamVLL[7] + 3*Sqrt[musqSU2]*lamVLL[8] + Sqrt[musqU1]*lamVLL[9]))/Sqrt[musqU1] - 384*\[Lambda]m3d^2*(1 + 2*Log[mu3/(2*Sqrt[m13d])]) + 96*(3*g23d^2 + g13d^2*Y\[Phi]^2)*\[Lambda]m3d*(1 + 4*Log[mu3/(2*Sqrt[m13d])]) - 288*lamVL[5]^2*(1 + 2*Log[mu3/(2*Sqrt[musqSU2])]) + 576*g23d^2*lamVL[5]*(1 + 4*Log[mu3/(2*Sqrt[musqSU2])]) - 96*lamVL[6]^2*(1 + 2*Log[mu3/(2*Sqrt[musqU1])]))/(6144*Pi^2)}
]];


AppendTo[testList,
TestCreate[BetaFunctions3DUS[],
	{\[Mu]\[Sigma]3dUS->(3 \[Lambda]\[Sigma]3dUS^2)/(8 \[Pi]^2),\[Mu]13dUS->(\[Lambda]\[Sigma]3dUS c33dUS)/(8 \[Pi]^2)}
]];


AppendTo[testList,
TestCreate[PrintTadpolesUS["LO"],
	{\[Mu]13dUS->\[Mu]13d-(Sqrt[m13d] \[Mu]m3d)/(4 \[Pi])-(\[Mu]13d \[Mu]m3d^2)/(384 m13d^(3/2) \[Pi])}
]];


AppendTo[testList,
TestCreate[PrintPressureUS["LO"],
	m13d^(3/2)/(3 \[Pi])+musqSU2^(3/2)/(4 \[Pi])+(2 musqSU3^(3/2))/(3 \[Pi])+musqU1^(3/2)/(12 \[Pi])
]];
AppendTo[testList,
TestCreate[PrintPressureUS["NLO"],
	-1/128*(48*m13d*\[Lambda]1H3d + 64*Sqrt[m13d]*Sqrt[musqSU3]*lamVL[1] + 8*Sqrt[m13d]*Sqrt[musqU1]*lamVL[2] + 24*Sqrt[m13d]*Sqrt[musqSU2]*lamVL[3] + 48*Sqrt[musqSU2]*Sqrt[musqSU3]*lamVLL[1] + (80*musqSU3*lamVLL[2])/3 + 5*musqSU2*lamVLL[3] + 16*Sqrt[musqSU3]*Sqrt[musqU1]*lamVLL[7] + 6*Sqrt[musqSU2]*Sqrt[musqU1]*lamVLL[8] + musqU1*lamVLL[9])/Pi^2 + (\[Mu]m3d^2*(1/2 + Log[mu3/(2*Sqrt[m13d])]))/(64*Pi^2) - (3*g23d^2*(-1/16*m13d/Pi^2 - 2*(-1/16*m13d/Pi^2 - (m13d*(1/2 + Log[mu3/(2*Sqrt[m13d])]))/(8*Pi^2))))/4 - (g13d^2*Y\[Phi]^2*(-1/16*m13d/Pi^2 - 2*(-1/16*m13d/Pi^2 - (m13d*(1/2 + Log[mu3/(2*Sqrt[m13d])]))/(8*Pi^2))))/4 - (3*g23d^2*(-1/16*musqSU2/Pi^2 - 2*(-1/16*musqSU2/Pi^2 - (musqSU2*(1/2 + Log[mu3/(2*Sqrt[musqSU2])]))/(8*Pi^2))))/2 - 6*g33d^2*(-1/16*musqSU3/Pi^2 - 2*(-1/16*musqSU3/Pi^2 - (musqSU3*(1/2 + Log[mu3/(2*Sqrt[musqSU3])]))/(8*Pi^2)))
]];


report=TestReport[testList]
report["ResultsDataset"]
