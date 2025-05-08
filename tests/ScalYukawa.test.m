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


Group={"U1"};
CouplingName={g};
RepAdjoint={0};
RepScalar={{{0},"R"}};
RepFermion={{{0},"L"},{{0},"R"}};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJC,\[Mu]IJ,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];


(* ::Text:: *)
(*Tadpoles*)


(* \[Sigma] \[Phi] *)
InputInv={{1},{True}};
LinearTerm=CreateInvariant[Group,RepScalar,InputInv][[1]];
VLinear=\[Sigma] LinearTerm;
\[Lambda]1=GradTadpole[VLinear];


(* ::Text:: *)
(*Scalar-Mass terms*)


(* 1/2m^2\[Phi]^2 *)
InputInv={{1,1},{True,True}}; 
MassTerm=CreateInvariant[Group,RepScalar,InputInv][[1]];
VMass=msq/2 MassTerm;
\[Mu]ij=GradMass[VMass]//SparseArray;


(* ::Text:: *)
(*Scalar - Cubic terms*)


(* 1/6\[Gamma]\[Phi]^3 *)
InputInv={{1,1,1},{True,True,True}};
CubicTerm=CreateInvariant[Group,RepScalar,InputInv][[1]];
VCubic=\[Gamma]/6 CubicTerm;
\[Lambda]3=GradCubic[VCubic];


(* ::Text:: *)
(*Scalar - Quartic terms*)


(* 1/24\[Lambda]\[Phi]^4 *)
InputInv={{1,1,1,1},{True,True,True,True}};
QuarticTerm=CreateInvariant[Group,RepScalar,InputInv][[1]];
VQuartic=\[Lambda]/24 MassTerm^2;
\[Lambda]4=GradQuartic[VQuartic];


(* ::Text:: *)
(*Fermion-mass terms*)


(* m(Subscript[\[Psi], R]^+Subscript[\[Psi], L]+Subscript[\[Psi]^+, L]Subscript[\[Psi], R])*)
InputInv={{2,1},{False,True}}; (*Subscript[\[Psi], R]^+Subscript[\[Psi], L]*)
MassTerm1=CreateInvariantFermion[Group,RepFermion,InputInv][[1]];
InputInv={{1,2},{False,True}};  (*Subscript[\[Psi]^+, L]Subscript[\[Psi], R]*)
MassTerm2=CreateInvariantFermion[Group,RepFermion,InputInv][[1]];
\[Mu]IJ=m\[Psi]*GradMassFermion[MassTerm1];
\[Mu]IJC=m\[Psi]*GradMassFermion[MassTerm2];


(* ::Text:: *)
(*Yukawa terms*)


InputInv={{1,2,1},{True,False,True}};
YukawaDoublet1=CreateInvariantYukawa[Group,RepScalar,RepFermion,InputInv][[1]]//Simplify;
InputInv={{1,1,2},{True,False,True}}; 
YukawaDoublet2=CreateInvariantYukawa[Group,RepScalar,RepFermion,InputInv][[1]]//Simplify;
Ysff= y*GradYukawa[YukawaDoublet1];
YsffC=y*GradYukawa[YukawaDoublet2];


(* ::Section:: *)
(*Dimensional Reduction*)


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->True];
PerformDRhard[]


(* ::Section:: *)
(*Tests*)


(* ::Subsection:: *)
(*Test hard*)


(*Replacements={
	Thread[{c61,c62,c63,c64}->0]
}//Flatten;*)


testList={};


(* ::Subsubsection:: *)
(*Beta functions*)


AppendTo[testList,
TestCreate[BetaFunctions4D[],
	{\[Lambda]->(-48 y^4+8 y^2 \[Lambda]+3 \[Lambda]^2)/(16 \[Pi]^2),\[Gamma]->(-24 m\[Psi] y^3+3 \[Gamma] (2 y^2+\[Lambda]))/(16 \[Pi]^2),y->(5 y^3)/(16 \[Pi]^2),msq->(4 (msq-6 m\[Psi]^2) y^2+\[Gamma]^2+msq \[Lambda])/(16 \[Pi]^2),m\[Psi]->(3 m\[Psi] y^2)/(16 \[Pi]^2),\[Sigma]->(-8 m\[Psi]^3 y+msq \[Gamma]+2 y^2 \[Sigma])/(16 \[Pi]^2)}
]];


(* ::Subsubsection:: *)
(*Dimension 4 matching relations*)


AppendTo[testList,
TestCreate[PrintCouplings[],
	{\[Lambda]3d->(T (48 Lf y^4+32 \[Pi]^2 \[Lambda]-8 Lf y^2 \[Lambda]-3 Lb \[Lambda]^2))/(32 \[Pi]^2),\[Gamma]3d->(Sqrt[T] (6 Lf y^2 (8 m\[Psi] y-\[Gamma])+32 \[Pi]^2 \[Gamma]-3 Lb \[Gamma] \[Lambda]))/(32 \[Pi]^2)}
]];


AppendTo[testList,
TestCreate[PrintTadpoles["LO"]//Expand,
	{\[Sigma]3d->1/6 m\[Psi] T^(3/2) y+1/24 T^(3/2) \[Gamma]+\[Sigma]/Sqrt[T]}
]];
AppendTo[testList,
TestCreate[PrintTadpoles["NLO"],
	{\[Sigma]3d->1/(768 \[Pi]^2 Sqrt[T]) (192 Lf m\[Psi]^3 y-32 EulerGamma m\[Psi] T^2 y^3+28 Lb m\[Psi] T^2 y^3+4 Lf m\[Psi] T^2 y^3-24 Lb msq \[Gamma]-6 Lb T^2 y^2 \[Gamma]-4 EulerGamma T^2 \[Gamma] \[Lambda]+Lb T^2 \[Gamma] \[Lambda]+384 m\[Psi] T^2 y^3 Log[Glaisher]+48 T^2 \[Gamma] \[Lambda] Log[Glaisher]+8 Sqrt[T] \[Gamma]3d \[Lambda]3d Log[\[Mu]3/\[Mu]])}
]];


AppendTo[testList,
TestCreate[PrintScalarMass["LO"],
	{msq3d->msq+1/24 T^2 (4 y^2+\[Lambda])}
]];
AppendTo[testList,
TestCreate[PrintScalarMass["NLO"]//Simplify,
	{msq3d->1/(768 \[Pi]^2) (16 Lb T^2 y^4-24 Lb \[Gamma]^2-24 Lb msq \[Lambda]-6 Lb T^2 y^2 \[Lambda]-4 EulerGamma T^2 \[Lambda]^2+Lb T^2 \[Lambda]^2-2 Lf y^2 (48 msq-288 m\[Psi]^2+T^2 (4 y^2+\[Lambda]))+48 T^2 \[Lambda]^2 Log[Glaisher]+8 \[Lambda]3d^2 Log[\[Mu]3/\[Mu]])}
]];


(* ::Subsubsection:: *)
(*Dimension 6 matching relations*)


(*AppendTo[testList,
TestCreate[PrintCouplingsEffective[],
	{c613d->1/(6144 \[Pi]^4) (96 \[Pi]^2 T^2 (c61 (64 \[Pi]^2-36 Lf yt1^2+9 Lb (3 g2^2+g1^2 Y\[Phi]^2-24 \[Lambda]1H))-4 c63 Lb \[Lambda]m)+(9 g2^6+9 g1^2 g2^4 Y\[Phi]^2+9 g1^4 g2^2 Y\[Phi]^4+3 g1^6 Y\[Phi]^6+8 (-84 yt1^6+240 \[Lambda]1H^3+\[Lambda]m^3)) Zeta[3]),c623d->c62 T^2-(Lb T^2 (c64 \[Lambda]m+45 c62 \[Lambda]\[Sigma]))/(16 \[Pi]^2)+((\[Lambda]m^3+54 \[Lambda]\[Sigma]^3) Zeta[3])/(1536 \[Pi]^4),c633d->1/(256 \[Pi]^4) (-8 \[Pi]^2 T^2 (12 (c61+c64) Lb \[Lambda]m+c63 (-9 g2^2 Lb-32 \[Pi]^2+12 Lf yt1^2+Lb (-3 g1^2 Y\[Phi]^2+48 \[Lambda]1H+16 \[Lambda]m+6 \[Lambda]\[Sigma])))+\[Lambda]m (24 \[Lambda]1H^2+12 \[Lambda]1H \[Lambda]m+\[Lambda]m (2 \[Lambda]m+3 \[Lambda]\[Sigma])) Zeta[3]),c643d->1/(256 \[Pi]^4) (-4 \[Pi]^2 T^2 (12 (5 c62+c63) Lb \[Lambda]m+c64 (-64 \[Pi]^2+12 Lf yt1^2+Lb (-9 g2^2-3 g1^2 Y\[Phi]^2+24 \[Lambda]1H+32 \[Lambda]m+72 \[Lambda]\[Sigma])))+\[Lambda]m (3 \[Lambda]1H \[Lambda]m+(\[Lambda]m+3 \[Lambda]\[Sigma])^2) Zeta[3])}
]];*)


(* ::Subsubsection:: *)
(*Symmetric pressure*)


AppendTo[testList,
TestCreate[PrintPressure["LO"],
	(13 \[Pi]^2 T^4)/180
]];
AppendTo[testList,
TestCreate[PrintPressure["NLO"],
	-((T^2 (48 msq+96 m\[Psi]^2+T^2 (10 y^2+\[Lambda])))/1152)
]];
AppendTo[testList,
TestCreate[PrintPressure["NNLO"],
	-(1/(552960 \[Pi]^2))(8640 Lb msq^2+34560 Lf m\[Psi]^4-4320 Lb msq T^2 y^2+1440 Lf msq T^2 y^2+34560 EulerGamma m\[Psi]^2 T^2 y^2-17280 Lb m\[Psi]^2 T^2 y^2-633 T^4 y^4-8280 EulerGamma T^4 y^4+4200 Lb T^4 y^4-300 Lf T^4 y^4-1440 EulerGamma T^2 \[Gamma]^2+720 Lb T^2 \[Gamma]^2-720 Lb msq T^2 \[Lambda]+480 EulerGamma T^4 y^2 \[Lambda]-420 Lb T^4 y^2 \[Lambda]+60 Lf T^4 y^2 \[Lambda]+58 T^4 \[Lambda]^2-30 EulerGamma T^4 \[Lambda]^2+5832 T^4 y^4 Log[2]+540 T^4 \[Lambda]^2 Log[2]+2160 T^4 y^4 Log[4]-414720 m\[Psi]^2 T^2 y^2 Log[Glaisher]+8640 T^4 y^4 Log[Glaisher]+17280 T^2 \[Gamma]^2 Log[Glaisher]-1440 T^4 \[Lambda]^2 Log[Glaisher]+4320 T^4 y^4 Log[\[Pi]]+180 T^4 \[Lambda]^2 Log[\[Pi]]+90 T^4 \[Lambda]^2 Log[\[Pi] T]+2880 T^4 y^4 Log[4 \[Pi] T]-480 T^4 y^2 \[Lambda] Log[4 \[Pi] T]-180 T^4 \[Lambda]^2 Log[4 \[Pi] T]-4320 T^4 y^4 Log[\[Mu]/T]-180 T^4 \[Lambda]^2 Log[\[Mu]/T]-1440 T^4 y^4 Log[\[Mu]^2]+240 T^4 y^2 \[Lambda] Log[\[Mu]^2]+45 T^4 \[Lambda]^2 Log[\[Mu]^2]+43200 T^4 y^4 Derivative[1][Zeta][-3]+7200 T^4 \[Lambda]^2 Derivative[1][Zeta][-3])
]];


(* ::Subsubsection:: *)
(*Report*)


report=TestReport[testList]
report["ResultsDataset"]



