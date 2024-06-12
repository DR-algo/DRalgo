(* ::Package:: *)

(*Quit[];*)


SetDirectory[NotebookDirectory[]];
$LoadGroupMath=True;
<<..//DRalgo.m


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


(* ::Section:: *)
(*Tests*)


(* ::Subsection:: *)
(*Test hard*)


testList={};


AppendTo[testList,
TestCreate[BetaFunctions4D[],
	{g1^2->(g1^2 (-18 yt1^2 (Yq-Yu)^2+(18 yt1^2+g1^2 nF (3 Yd^2+Ye^2+2 Yl^2+6 Yq^2+3 Yu^2)) Y\[Phi]^2+g1^2 Y\[Phi]^4))/(48 \[Pi]^2 Y\[Phi]^2),g2^2->(g2^4 (-43+8 nF))/(48 \[Pi]^2),g3^2->(g3^4 (-33+4 nF))/(24 \[Pi]^2),\[Lambda]1H->(9 g2^4-48 yt1^4+3 g1^4 Y\[Phi]^4+6 g2^2 (g1^2 Y\[Phi]^2-12 \[Lambda]1H)-24 g1^2 Y\[Phi]^2 \[Lambda]1H+96 \[Lambda]1H (yt1^2+2 \[Lambda]1H)+4 \[Lambda]m^2)/(128 \[Pi]^2),\[Lambda]m->(\[Lambda]m (-9 g2^2+12 yt1^2-3 g1^2 Y\[Phi]^2+24 \[Lambda]1H+8 \[Lambda]m+12 \[Lambda]\[Sigma]))/(32 \[Pi]^2),\[Lambda]\[Sigma]->(\[Lambda]m^2+9 \[Lambda]\[Sigma]^2)/(8 \[Pi]^2),\[Mu]3->(3 (6 \[Lambda]\[Sigma] \[Mu]3+\[Lambda]m \[Mu]m))/(16 \[Pi]^2),\[Mu]m->(-3 (3 g2^2-4 yt1^2+g1^2 Y\[Phi]^2-8 \[Lambda]1H) \[Mu]m+8 \[Lambda]m (\[Mu]3+\[Mu]m))/(32 \[Pi]^2),yt1->(yt1 (-9 g2^2-32 g3^2+18 yt1^2-3 g1^2 (2 Yq Yu+Y\[Phi]^2)))/(64 \[Pi]^2),m1->(-3 m1 (3 g2^2-4 yt1^2+g1^2 Y\[Phi]^2-8 \[Lambda]1H)+\[Mu]m^2+2 \[Lambda]m \[Mu]\[Sigma])/(32 \[Pi]^2),\[Mu]\[Sigma]->(4 m1 \[Lambda]m+4 \[Mu]3^2+\[Mu]m^2+6 \[Lambda]\[Sigma] \[Mu]\[Sigma])/(16 \[Pi]^2),\[Mu]1->(m1 \[Mu]m+\[Mu]3 \[Mu]\[Sigma])/(8 \[Pi]^2)}
]];


AppendTo[testList,
TestCreate[PrintCouplings[],
	{g13d^2->1/(96 \[Pi]^2 Y\[Phi]^2) g1^2 T (18 Lf yt1^2 (Yq-Yu)^2+(96 \[Pi]^2-Lf (18 yt1^2+g1^2 nF (3 Yd^2+Ye^2+2 Yl^2+6 Yq^2+3 Yu^2))) Y\[Phi]^2-g1^2 Lb Y\[Phi]^4),g23d^2->g2^2 T+(g2^4 (4+43 Lb-8 Lf nF) T)/(96 \[Pi]^2),g33d^2->g3^2 T+(g3^4 (3+33 Lb-4 Lf nF) T)/(48 \[Pi]^2),\[Lambda]1H3d->1/(256 \[Pi]^2) T (48 Lf yt1^4+(2-3 Lb) (3 g2^4+2 g1^2 g2^2 Y\[Phi]^2+g1^4 Y\[Phi]^4)+256 \[Pi]^2 \[Lambda]1H+24 (3 g2^2 Lb-4 Lf yt1^2+g1^2 Lb Y\[Phi]^2) \[Lambda]1H-4 Lb (48 \[Lambda]1H^2+\[Lambda]m^2)),\[Lambda]m3d->(T \[Lambda]m (9 g2^2 Lb+64 \[Pi]^2-12 Lf yt1^2+3 g1^2 Lb Y\[Phi]^2-4 Lb (6 \[Lambda]1H+2 \[Lambda]m+3 \[Lambda]\[Sigma])))/(64 \[Pi]^2),\[Lambda]\[Sigma]3d->T \[Lambda]\[Sigma]-(Lb T (\[Lambda]m^2+9 \[Lambda]\[Sigma]^2))/(16 \[Pi]^2),\[Mu]33d->1/2 Sqrt[T] (2 \[Mu]3-(3 Lb (6 \[Lambda]\[Sigma] \[Mu]3+\[Lambda]m \[Mu]m))/(16 \[Pi]^2)),\[Mu]m3d->(Sqrt[T] (4 (16 \[Pi]^2-3 Lf yt1^2) \[Mu]m+3 Lb (3 g2^2+g1^2 Y\[Phi]^2-8 \[Lambda]1H) \[Mu]m-8 Lb \[Lambda]m (\[Mu]3+\[Mu]m)))/(64 \[Pi]^2)}
]];


AppendTo[testList,
TestCreate[PrintTadpoles["LO"],
	{\[Mu]13d->(12 \[Mu]1+T^2 (\[Mu]3+\[Mu]m))/(12 Sqrt[T])}
]];
AppendTo[testList,
TestCreate[PrintTadpoles["NLO"],
	{\[Mu]13d->1/(768 \[Pi]^2 Sqrt[T]) (-Lb (T^2 (-12 \[Lambda]\[Sigma] \[Mu]3+2 \[Lambda]m (4 \[Mu]3-5 \[Mu]m)+3 (9 g2^2+6 yt1^2+3 g1^2 Y\[Phi]^2+8 \[Lambda]1H) \[Mu]m)+48 (m1 \[Mu]m+\[Mu]3 \[Mu]\[Sigma]))+2 T^2 (3 Lf yt1^2 \[Mu]m+g1^2 Y\[Phi]^2 \[Mu]m+6 EulerGamma (-4 \[Lambda]\[Sigma] \[Mu]3+(3 g2^2+g1^2 Y\[Phi]^2-2 \[Lambda]m) \[Mu]m)+3 g2^2 \[Mu]m (1-72 Log[Glaisher])+288 \[Lambda]\[Sigma] \[Mu]3 Log[Glaisher]-72 g1^2 Y\[Phi]^2 \[Mu]m Log[Glaisher]+144 \[Lambda]m \[Mu]m Log[Glaisher])+24 Sqrt[T] (4 \[Lambda]\[Sigma]3d \[Mu]33d-(3 g23d^2+g13d^2 Y\[Phi]^2-2 \[Lambda]m3d) \[Mu]m3d) Log[\[Mu]3/\[Mu]])}
]];


AppendTo[testList,
TestCreate[PrintTemporalScalarCouplings[],
	{\[Lambda]VLL[1]->-((g2^2 g3^2 nF T)/(4 \[Pi]^2)),\[Lambda]VLL[2]->(g3^4 (9-2 nF) T)/(4 \[Pi]^2),\[Lambda]VLL[3]->(g2^4 (17-4 nF) T)/(8 \[Pi]^2),\[Lambda]VLL[4]->-((I g1 g2^3 nF T (Yl+3 Yq))/(24 \[Pi]^2)),\[Lambda]VLL[5]->-((I g1 g3^3 nF T (Yd-2 Yq+Yu))/(48 \[Pi]^2)),\[Lambda]VLL[6]->-((g1 g3^3 nF T (Yd+2 Yq+Yu))/(16 \[Pi]^2)),\[Lambda]VLL[7]->-((g1^2 g3^2 nF T (Yd^2+2 Yq^2+Yu^2))/(8 \[Pi]^2)),\[Lambda]VLL[8]->(g1^2 g2^2 T (-nF (Yl^2+3 Yq^2)+Y\[Phi]^2))/(8 \[Pi]^2),\[Lambda]VLL[9]->-((g1^4 T (nF (3 Yd^4+Ye^4+2 Yl^4+6 Yq^4+3 Yu^4)-2 Y\[Phi]^4))/(16 \[Pi]^2)),\[Lambda]VVSL[1]->(g2^2 Sqrt[T] \[Mu]m)/(16 \[Pi]^2),\[Lambda]VVSL[2]->(g1^2 Sqrt[T] Y\[Phi]^2 \[Mu]m)/(16 \[Pi]^2),\[Lambda]VL[1]->-((g3^2 T yt1^2)/(4 \[Pi]^2)),\[Lambda]VL[2]->1/192 g1^2 T (96 Y\[Phi]^2+1/\[Pi]^2 (-36 Lf Yq yt1^2 Yu+18 yt1^2 ((-2+Lf) Yu^2-Lf Y\[Phi]^2)+6 Yq^2 (3 (-2+Lf) yt1^2-g1^2 (-1+Lf) nF Y\[Phi]^2)+Y\[Phi]^2 (9 g2^2+g1^2 (-((-1+Lf) nF (3 Yd^2+Ye^2+2 Yl^2+3 Yu^2))-(-1+Lb) Y\[Phi]^2)+72 \[Lambda]1H))),\[Lambda]VL[3]->(g2^2 T (g2^2 (51+43 Lb+8 nF-8 Lf nF)+96 \[Pi]^2-36 yt1^2+3 g1^2 Y\[Phi]^2+72 \[Lambda]1H))/(192 \[Pi]^2),\[Lambda]VL[4]->1/384 g1 g2 T (192 Y\[Phi]+1/\[Pi]^2 (-36 (-2+Lf) Yq yt1^2+36 Lf yt1^2 (Yu-Y\[Phi])+g2^2 (-12+43 Lb+8 nF) Y\[Phi]-6 g1^2 (-1+Lf) nF Yq^2 Y\[Phi]-Lf nF (8 g2^2+g1^2 (3 Yd^2+Ye^2+2 Yl^2+3 Yu^2)) Y\[Phi]+g1^2 Y\[Phi] (nF (3 Yd^2+Ye^2+2 Yl^2+3 Yu^2)-(-4+Lb) Y\[Phi]^2)+48 Y\[Phi] \[Lambda]1H)),\[Lambda]VL[5]->(g2^2 T \[Lambda]m)/(8 \[Pi]^2),\[Lambda]VL[6]->(g1^2 T Y\[Phi]^2 \[Lambda]m)/(8 \[Pi]^2)}
]];


AppendTo[testList,
TestCreate[PrintDebyeMass["LO"],
	{\[Mu]sqSU2->1/6 g2^2 (5+2 nF) T^2,\[Mu]sqSU3->1/3 g3^2 (3+nF) T^2,\[Mu]sqU1->1/24 g1^2 T^2 (nF (3 Yd^2+Ye^2+2 Yl^2+6 Yq^2+3 Yu^2)+4 Y\[Phi]^2)}
]];
AppendTo[testList,
TestCreate[PrintDebyeMass["NLO"],
{\[Mu]sqSU2->1/(1152 \[Pi]^2) g2^2 (144 m1+T^2 (-144 g3^2 nF+9 g1^2 (-nF (Yl^2+3 Yq^2)+Y\[Phi]^2)+6 (-3 yt1^2+12 \[Lambda]1H+\[Lambda]m)+g2^2 (207-4 EulerGamma (5+2 nF) (-43+8 nF)+4 nF (11+8 nF-192 Log[2])))+4 g2^2 T^2 (8 nF (-7+2 nF) Log[\[Pi] T]-39 Log[4 \[Pi] T]+(39+8 (7-2 nF) nF) Log[\[Mu]]+2 (88-5 nF) Log[\[Mu]/(4 \[Pi] T)])),\[Mu]sqSU3->1/(1152 \[Pi]^2) g3^2 T^2 (-9 (6 g2^2 nF+8 yt1^2+g1^2 nF (Yd^2+2 Yq^2+Yu^2))-8 g3^2 (-45+2 EulerGamma (3+nF) (-33+4 nF)+nF (-3-4 nF+132 Log[2]))+16 g3^2 (nF (-21+4 nF) (Log[\[Pi]]+Log[T]-Log[\[Mu]])+99 Log[\[Mu]/(4 \[Pi] T)])),\[Mu]sqU1->1/(576 \[Pi]^2) (-27 g1^2 T^2 yt1^2 (Yq^2+Yu^2)-9/4 g1^2 nF T^2 (6 g2^2 (Yl^2+3 Yq^2)+16 g3^2 (Yd^2+2 Yq^2+Yu^2)+g1^2 (3 Yd^4+Ye^4+2 Yl^4+6 Yq^4+3 Yu^4))+72 g1^2 m1 Y\[Phi]^2+18 g1^2 T^2 yt1^2 Y\[Phi]^2+9/2 T^2 (3 g1^2 g2^2 Y\[Phi]^2+g1^4 Y\[Phi]^4)+3 g1^2 T^2 Y\[Phi]^2 (12 \[Lambda]1H+\[Lambda]m)-2 g1^4 T^2 Y\[Phi]^4 (1+EulerGamma-Log[4 \[Pi] T]+Log[\[Mu]])-1/4 g1^4 nF^2 T^2 (3 Yd^2+Ye^2+2 Yl^2+6 Yq^2+3 Yu^2)^2 (-1+2 EulerGamma-2 Log[\[Pi] T]+2 Log[\[Mu]])-1/2 g1^4 nF T^2 (3 Yd^2+Ye^2+2 Yl^2+6 Yq^2+3 Yu^2) Y\[Phi]^2 (-1+5 EulerGamma+Log[1/(4 \[Pi]^5)]+5 Log[\[Mu]/T]))}
]];


AppendTo[testList,
TestCreate[PrintScalarMass["LO"],
	{m13d->m1+1/48 T^2 (9 g2^2+12 yt1^2+3 g1^2 Y\[Phi]^2+24 \[Lambda]1H+2 \[Lambda]m),\[Mu]\[Sigma]3d->1/12 T^2 (2 \[Lambda]m+3 \[Lambda]\[Sigma])+\[Mu]\[Sigma]}
]];
AppendTo[testList,
TestCreate[PrintScalarMass["NLO"],
{m13d->1/(9216 \[Pi]^2) (-1728 Lf m1 yt1^2-1152 g3^2 T^2 yt1^2-384 g3^2 Lb T^2 yt1^2+1536 g3^2 Lf T^2 yt1^2-36 g1^2 T^2 Yq^2 yt1^2+18 g1^2 Lb T^2 Yq^2 yt1^2+90 g1^2 Lf T^2 Yq^2 yt1^2+216 Lb T^2 yt1^4-144 g1^2 T^2 Yq yt1^2 Yu-108 g1^2 Lb T^2 Yq yt1^2 Yu+108 g1^2 Lf T^2 Yq yt1^2 Yu-36 g1^2 T^2 yt1^2 Yu^2+18 g1^2 Lb T^2 yt1^2 Yu^2+90 g1^2 Lf T^2 yt1^2 Yu^2+432 g1^2 Lb m1 Y\[Phi]^2+6 g1^4 nF T^2 Yd^2 Y\[Phi]^2-27 g1^4 Lb nF T^2 Yd^2 Y\[Phi]^2+9 g1^4 Lf nF T^2 Yd^2 Y\[Phi]^2+2 g1^4 nF T^2 Ye^2 Y\[Phi]^2-9 g1^4 Lb nF T^2 Ye^2 Y\[Phi]^2+3 g1^4 Lf nF T^2 Ye^2 Y\[Phi]^2+4 g1^4 nF T^2 Yl^2 Y\[Phi]^2-18 g1^4 Lb nF T^2 Yl^2 Y\[Phi]^2+6 g1^4 Lf nF T^2 Yl^2 Y\[Phi]^2+12 g1^4 nF T^2 Yq^2 Y\[Phi]^2-54 g1^4 Lb nF T^2 Yq^2 Y\[Phi]^2+18 g1^4 Lf nF T^2 Yq^2 Y\[Phi]^2+108 g1^2 Lb T^2 yt1^2 Y\[Phi]^2-108 g1^2 Lf T^2 yt1^2 Y\[Phi]^2+6 g1^4 nF T^2 Yu^2 Y\[Phi]^2-27 g1^4 Lb nF T^2 Yu^2 Y\[Phi]^2+9 g1^4 Lf nF T^2 Yu^2 Y\[Phi]^2+2 g1^4 T^2 Y\[Phi]^4-126 EulerGamma g1^4 T^2 Y\[Phi]^4+66 g1^4 Lb T^2 Y\[Phi]^4-3456 Lb m1 \[Lambda]1H-1296 Lb T^2 yt1^2 \[Lambda]1H-432 Lf T^2 yt1^2 \[Lambda]1H+144 g1^2 T^2 Y\[Phi]^2 \[Lambda]1H+864 EulerGamma g1^2 T^2 Y\[Phi]^2 \[Lambda]1H-432 g1^2 Lb T^2 Y\[Phi]^2 \[Lambda]1H-3456 EulerGamma T^2 \[Lambda]1H^2-72 Lf T^2 yt1^2 \[Lambda]m+18 g1^2 Lb T^2 Y\[Phi]^2 \[Lambda]m-144 Lb T^2 \[Lambda]1H \[Lambda]m-144 EulerGamma T^2 \[Lambda]m^2+24 Lb T^2 \[Lambda]m^2-72 Lb T^2 \[Lambda]m \[Lambda]\[Sigma]-144 Lb \[Mu]m^2-288 Lb \[Lambda]m \[Mu]\[Sigma]+54 g2^2 (Lb (24 m1+T^2 (7 yt1^2+8 g1^2 Y\[Phi]^2-24 \[Lambda]1H+\[Lambda]m))-T^2 ((2+Lf) yt1^2-8 \[Lambda]1H (1+6 EulerGamma-72 Log[Glaisher])+2 g1^2 Y\[Phi]^2 (1+5 EulerGamma-60 Log[Glaisher])))+6 g2^4 T^2 (167+243 EulerGamma+8 nF+12 Lf nF-3 Lb (47+12 nF)-2916 Log[Glaisher])+1512 g1^4 T^2 Y\[Phi]^4 Log[Glaisher]-10368 g1^2 T^2 Y\[Phi]^2 \[Lambda]1H Log[Glaisher]+41472 T^2 \[Lambda]1H^2 Log[Glaisher]+1728 T^2 \[Lambda]m^2 Log[Glaisher]-36 Log[\[Mu]3/\[Mu]] (39 g23d^4-5 g13d^4 Y\[Phi]^4+48 g13d^2 Y\[Phi]^2 \[Lambda]1H3d+g23d^2 (-18 g13d^2 Y\[Phi]^2+144 \[Lambda]1H3d+96 \[Lambda]VL[3])-8 (24 \[Lambda]1H3d^2+\[Lambda]m3d^2-48 g33d^2 \[Lambda]VL[1]+8 \[Lambda]VL[1]^2+\[Lambda]VL[2]^2+3 \[Lambda]VL[3]^2+6 \[Lambda]VL[4]^2))),\[Mu]\[Sigma]3d->-(1/(384 \[Pi]^2))(Lb (48 m1 \[Lambda]m+27 g2^2 T^2 \[Lambda]m+18 T^2 yt1^2 \[Lambda]m+9 g1^2 T^2 Y\[Phi]^2 \[Lambda]m+24 T^2 \[Lambda]1H \[Lambda]m-10 T^2 \[Lambda]m^2+12 T^2 \[Lambda]m \[Lambda]\[Sigma]-18 T^2 \[Lambda]\[Sigma]^2+48 \[Mu]3^2+12 \[Mu]m^2+72 \[Lambda]\[Sigma] \[Mu]\[Sigma])-2 T^2 (3 Lf yt1^2 \[Lambda]m+g1^2 Y\[Phi]^2 \[Lambda]m+6 EulerGamma g1^2 Y\[Phi]^2 \[Lambda]m-12 EulerGamma \[Lambda]m^2-36 EulerGamma \[Lambda]\[Sigma]^2+3 g2^2 \[Lambda]m (1+6 EulerGamma-72 Log[Glaisher])-72 g1^2 Y\[Phi]^2 \[Lambda]m Log[Glaisher]+144 \[Lambda]m^2 Log[Glaisher]+432 \[Lambda]\[Sigma]^2 Log[Glaisher])+12 Log[\[Mu]3/\[Mu]] (2 g13d^2 Y\[Phi]^2 \[Lambda]m3d-4 \[Lambda]m3d^2-12 \[Lambda]\[Sigma]3d^2-3 \[Lambda]VL[5]^2+6 g23d^2 (\[Lambda]m3d+2 \[Lambda]VL[5])-\[Lambda]VL[6]^2))}
]];


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



