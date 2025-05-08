(* ::Package:: *)

(* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ *)

(* :Title: HardToSoft							    *)

(*
      	This software is covered by the GNU General Public License 3.
       	Copyright (C) 2021-2025 Andreas Ekstedt
       	Copyright (C) 2021-2025 Philipp Schicho
       	Copyright (C) 2021-2025 Tuomas V.I. Tenkanen

*)

(* :Summary:	Dimensonal reduction from hard to soft scale.		    *)

(* ------------------------------------------------------------------------ *)


(* ::Section::Closed:: *)
(*Help routines*)


(* Helper functions *)
VerbosePrint[msg_] := If[verbose, Print[msg]]
LogPrint[msg_] := If[verbose, Print[DateString[{"Hour", ":", "Minute"}], " \[LongDash] ", msg]];


(*
	Rewrites short-hand constants with Lb and Lf.
*)
ReplaceLb={
	Lbb->2Lb+2(12Log[Glaisher]-EulerGamma),
	Lbbb->Lb/12+1/6(12Log[Glaisher]-EulerGamma),
	LBF->5/2Lb-1/2Lf+2(12Log[Glaisher]-EulerGamma),
	LFF->1/2Lf+3/2Lb+2(12Log[Glaisher]-EulerGamma),
	LFB->Lb+Lf+2(12Log[Glaisher]-EulerGamma),
	Lfff->1/2Lf-3/2Lb-2(12Log[Glaisher]-EulerGamma),
	LbbM->2Lb+2(12Log[Glaisher]-EulerGamma)};


(*
  Contract:
  Performs a TensorContract after forming the TensorProduct of any number of tensors.
  Uses Inactive TensorProduct to delay evaluation until after contraction.
  Automatically activates the resulting expression.
*)
Contract[tensors__SparseArray, indices_] := Module[
  {
    product,
    tensorList = {tensors}
  },
  If[Length[tensorList] < 2,
    Message[Contract::argcount, Length[tensorList]];
    Return[$Failed];
  ];
  
  (* Build the TensorProduct without evaluation *)
  product = Inactive[TensorProduct][tensors];
  
  (* Perform contraction and activate result *)
  Activate @ TensorContract[product, indices]
];


(*
  SimplifySparse:
  Simplifies the non-zero elements of a SparseArray.
  Useful workaround for limitations in Simplify behavior with SparseArray in Mathematica 13.
*)
SimplifySparse[s_SparseArray] := Module[
  {
    simplifiedValues,
    positions,
    shape
  },
  
  (* Extract and simplify the non-zero values *)
  simplifiedValues = Simplify[s["NonzeroValues"]];
  
  (* Extract positions and array dimensions *)
  positions = s["NonzeroPositions"];
  shape = Dimensions[s];
  
  (* Return a new SparseArray with simplified entries *)
  SparseArray[positions -> simplifiedValues, shape, 0]
];


(*
	Finds relations between parameters in two different basis.
*)
CompareInvariants[Tens1I_,Tens2I_]:=Module[{Tens1P=Tens1I,Tens2P=Tens2I},
	Sol1=Flatten[Tens1P,{{3},{4},{1,2}}] . Flatten[Tens1P,{{1,2}}]//Normal;
	Sol2=Flatten[Tens2P,{{3},{4},{1,2}}] . Flatten[Tens2P,{{1,2}}]//Normal//ReplaceAll[#,PrintIdentification[]]&;

	dimHelp=Dimensions[Tens1P];

	If[dimHelp[[1]]==dimHelp[[2]],
		Sol3=TensorContract[Tens1P,{1,2}]//Normal;
		Sol4=TensorContract[Tens2P,{1,2}]//Normal//ReplaceAll[#,PrintIdentification[]]&;
		Sol5=Activate@TensorContract[Inactive@TensorProduct[Tens1P,Tens1P,Tens1P],{{2,5},{6,9},{10,1}}]//Normal;
		Sol6=Activate@TensorContract[Inactive@TensorProduct[Tens2P,Tens2P,Tens2P],{{2,5},{6,9},{10,1}}]//Normal//ReplaceAll[#,PrintIdentification[]]&;

		SolInvariants=Solve[Sol1==Sol2&&Sol3==Sol4&&Sol5==Sol6,Sol5//Variables][[1]];
	,
		SolInvariants=Solve[Sol1==Sol2,Sol1//Variables][[1]];
	];
	
	Return[SolInvariants]
];


(*
	Comparing symbolic with numeric tensors
*)
CompareTensors[mat1I_,mat2I_]:=Module[{mat1P=mat1I,mat2P=mat2I},
	VarMat1=#->RandomInteger[10000]&/@(Normal[mat1P]//Variables);
	VarMat2=#->RandomInteger[10000]&/@(Normal[mat2P]//Variables);
	mat1P=Normal[mat1P]/.VarMat1;
	mat2P=Normal[mat2P]/.VarMat2;
Return[mat1P==mat2P]
];


(*
  CreateHelpTensors:
  Constructs intermediate tensors from scalar-vector, fermion-vector, Yukawa, and structure constant couplings.
  These are used in higher-loop pressure calculations and symmetry analyses.
*)
ClearAll[CreateHelpTensors];
CreateHelpTensors[] := Module[{},
  VerbosePrint["Creating Help Tensors"];
  
  (* --- Tensors from two scalar-vector Trilinear Couplings --- *)
  Habij = Contract[gvss, gvss, {{3, 5}}] // Transpose[#, {1, 3, 2, 4}] & // SimplifySparse;
  Hg = TensorContract[Habij, {{3, 4}}] // SimplifySparse;
  \[CapitalLambda]g = TensorContract[Habij, {{1, 2}}] // SimplifySparse;
  HabijV = (Habij + Transpose[Habij, {2, 1, 3, 4}]) // SparseArray // SimplifySparse;

  (* --- Tensors from two structure constants --- *)
  GabcdV = gvvv . gvvv // SparseArray // SimplifySparse;

  (* --- Tensors from two fermion-vector trilinear couplings --- *)
  HabIJF = Contract[gvff, gvff, {{3, 5}}] // Transpose[#, {1, 3, 2, 4}] & // SimplifySparse;

  (* --- Tensors from two Yukawa coupling contractions --- *)
  Ysij = Contract[Ysff, YsffC, {{2, 5}, {3, 6}}] // SimplifySparse;
  YsijC = Contract[YsffC, Ysff, {{2, 5}, {3, 6}}] // SimplifySparse;

  (* --- Mode-dependent contributions --- *)
  If[mode >= 1,
    
    (* Tensors from two scalar quartic contractions *)
    \[CapitalLambda]\[Lambda] = Flatten[\[Lambda]4, {{1}, {2}, {3, 4}}] . Flatten[\[Lambda]4, {1, 2}] // SparseArray // SimplifySparse;

    (* Tensors from Yukawa \[Times] Yukawa\:1d9c *)
    YTemp = Ysff . Transpose[YsffC] 
      // Transpose[#, {1, 3, 2, 4}] & 
      // Transpose[#, {1, 2, 4, 3}] & 
      // SimplifySparse;

    YTempC = YsffC . Transpose[Ysff] 
      // Transpose[#, {1, 3, 2, 4}] & 
      // Transpose[#, {1, 2, 4, 3}] & 
      // SimplifySparse;

    (* Flatten to contract internal indices and obtain invariants *)
    Yhelp = Flatten[YTemp, {{1}, {2}, {3, 4}}] . Flatten[YTemp, {4, 3}] // SparseArray // SimplifySparse;
    YhelpC = Flatten[YTempC, {{1}, {2}, {3, 4}}] . Flatten[YTempC, {4, 3}] // SparseArray // SimplifySparse;
  ];
];


(* ::Section::Closed:: *)
(*Pressure calculations*)


(*
  SymmetricPhaseEnergy:
  Computes the total pressure in the symmetric phase of the soft theory.
  Includes contributions up to NNLO.
*)
SymmetricPhaseEnergy[] := Module[
  {
    totalPressure
  },

  (* Counterterms are required for higher-order SymmetricPhaseNLO and SymmetricPhaseNNLO contributions *)
  CounterTerm[];

  (* Convention: pressure is minus the free energy *)
  totalPressure = {
    -SymmetricPhaseLO[],
    -SymmetricPhaseNLO[],
    -SymmetricPhaseNNLO[]
  };

  (* Store result *)
  SymmEnergy = totalPressure;
];


(* ::Subsection:: *)
(*LO symmetric pressure*)


(*
  SymmetricPhaseLO:
  Calculates the 1-loop pressure in the symmetric phase of the soft theory.
*)
SymmetricPhaseLO[] := Module[
  {
    isEmptyTensorQ,                  (* helper function *)
    scalarEmptyQ, vectorEmptyQ, fermionEmptyQ,
    scalarContribution, vectorContribution, fermionContribution
  },
  
  VerbosePrint["Calculating Leading-Order \!\(\*SuperscriptBox[\(T\), \(4\)]\) Terms"];

  (* --- Helper to check whether a tensor is numerically empty --- *)
  isEmptyTensorQ[tensor_, shape_] := (
    AllTrue[Flatten[tensor], NumericQ] &&
    CompareTensors[tensor, EmptyArray[shape]]
  );

  (* --- Scalar contribution --- *)
  scalarEmptyQ = And @@ {
    isEmptyTensorQ[gvss, {nv, nf, nf}],
    isEmptyTensorQ[Ysff, {ns, nf, nf}],
    isEmptyTensorQ[\[Mu]ij, {ns, ns}],
    isEmptyTensorQ[\[Lambda]4, {ns, ns, ns, ns}],
    isEmptyTensorQ[\[Lambda]3, {ns, ns, ns}],
    isEmptyTensorQ[\[Lambda]1, {ns}]
  };
  scalarContribution = If[scalarEmptyQ, 0, Sum[-(\[Pi]^2/90) T^4, {a, 1, ns}]];

  (* --- Vector contribution --- *)
  vectorEmptyQ = And @@ {
    isEmptyTensorQ[gvss, {nv, nf, nf}],
    isEmptyTensorQ[gvvv, {nv, nv, nv}],
    isEmptyTensorQ[gvff, {nv, nf, nf}]
  };
  vectorContribution = If[vectorEmptyQ, 0, Sum[-(\[Pi]^2/90) T^4 * 2, {a, 1, nv}]];

  (* --- Fermion contribution --- *)
  fermionEmptyQ = And @@ {
    isEmptyTensorQ[gvff, {nv, nf, nf}],
    isEmptyTensorQ[Ysff, {ns, nf, nf}],
    isEmptyTensorQ[\[Mu]IJF, {nf, nf}]
  };
  fermionContribution = If[
    fermionEmptyQ, 
    0, 
    Sum[-(7 \[Pi]^2)/720 T^4 * 2 * NFMat[[a, a]], {a, 1, nf}]
  ];

  (* --- Return total contribution, cleaning internal formatting --- *)
  ToExpression @ StringReplace[
    ToString[StandardForm[scalarContribution + vectorContribution + fermionContribution]],
    "DRalgo`Private`" -> ""
  ]
];


(* ::Subsection:: *)
(*NLO symmetric pressure*)


(*
  SymmetricPhaseNLO:
  Calculates the 2-loop pressure in the symmetric phase of the soft theory.
  Uses notation following Martin (arXiv:1808.07615).
*)
SymmetricPhaseNLO[] := Module[
	{
	tempCoeff, Vss, Vssv, Vvs, Vvv, Vvvv, Vggv, VFFs, VFFv,
	massInsertionScalar, massInsertionFermion,
	VssLambda6, HabIJFnF, totalNLOPressure
	},

  VerbosePrint["Calculating NLO \!\(\*SuperscriptBox[\(T\), \(4\)]\) Terms"];

  (* --- Two-loop diagram contributions --- *)
  (* Scalar quartic self-interaction term *)
  tempCoeff = (1/144) * T^4;
  Vss = (1/8) * tempCoeff * TensorContract[\[Lambda]4, {{1, 2}, {3, 4}}];

  (* Scalar-vector interaction (type I) *)
  tempCoeff = -(1/144) * T^4;
  Vssv = (1/4) * tempCoeff * Total[Flatten[gvss] Flatten[gvss]];

  (* Scalar-vector interaction (type II) *)
  tempCoeff = (1/48) * T^4;
  Vvs = (1/2) * tempCoeff * Total[Flatten[gvss] Flatten[gvss]];

  (* Vector cubic self-interaction *)
  tempCoeff = (3/64) * T^4;
  Vvv = (1/4) * tempCoeff * Total[Flatten[gvvv] Flatten[gvvv]];

  (* Vector quartic self-interaction *)
  tempCoeff = -(13/192) * T^4;
  Vvvv = (1/12) * tempCoeff * Total[Flatten[gvvv] Flatten[gvvv]];

  (* Ghost-vector interaction *)
  tempCoeff = (1/288) * T^4;
  Vggv = (1/4) * tempCoeff * Total[Flatten[gvvv] Flatten[gvvv]];

  (* Fermion-scalar Yukawa interaction *)
  tempCoeff = (5/576) * T^4;
  VFFs = (1/2) * tempCoeff * TensorContract[Ysij, {{1, 2}}];

  (* --- Fermion-vector coupling correction with nF factor --- *)
  tempCoeff = 5/288 T^4;
  HabIJFnF = HabIJF . NFMat;
  VFFv = (1/2) * tempCoeff * TensorContract[HabIJFnF, {{1, 2}, {3, 4}}];

  (* --- Mass insertions in 1-loop diagrams --- *)
  massInsertionScalar = (1/2) * Tr[\[Mu]ij] * T^2 / 12;
  massInsertionFermion = (1/2) * T^2 / 12 * Tr[\[Mu]IJF . \[Mu]IJFC];

  (* --- Higher-dimensional operator (\[Lambda]6) contribution --- *)
  VssLambda6 = If[
    mode >= 3,
    (1/82944) * T^6 * TensorContract[\[Lambda]6, {{1, 2}, {3, 4}, {5, 6}}],
    0
  ];

  (* --- Total expression and cleanup --- *)
  totalNLOPressure = (
    + massInsertionScalar
    + massInsertionFermion
    + Vss + Vssv + Vvs + Vvv + Vvvv + Vggv + VFFs + VFFv
    + VssLambda6
    );

  ToExpression @ StringReplace[
    ToString[StandardForm[FullSimplify[totalNLOPressure]]],
    "DRalgo`Private`" -> ""
  ]
];


(* ::Subsection:: *)
(*NNLO symmetric pressure*)


(*
  Computes the 3-loop pressure in the soft theory.
  Diagrammatic contributions follow the naming and structure of arXiv:1709.02397.
*)
SymmetricPhaseNNLO[]:=Module[
	{
	  \[Kappa], I1, I2M2, I1I2, dI1I2, I211M020, I2M2I2, I3M2I1, I1p, 
	  I3M2p, I2p, I4M4p, I2M2p, IF1p, IF2p, IF3M2p, IF2M2p,
	  EPre, M00, LInt, EFPre, EFInt, N00, MF1M1, MF00, MFM2P2, N2M2, MP2M2,
	  Temp, Temp1, Temp2, Temp3, Temp4, Temp5, 
	  Help1, Help2, Help3, GHelp, GHelp2, I1Temp, HFab, GabV, gvffnF, gvffNf, \[Lambda]4Eff,
	  HabIJFnF,
	  DiagramsSimplified, DiagramsTotal
	},
	
	VerbosePrint["Calculating NNLO \!\(\*SuperscriptBox[\(T\), \(4\)]\) Terms"];

	(* --- Precompute loop integrals --- *)
	(* One- and two-loop integrals appearing in scalar, gauge and fermionic diagrams *)
	\[Kappa]=1/(16 \[Pi]^2);
	I1=T^2/(12 \[Epsilon])+T^2 Lbbb;
	I2M2=T^2/(24 )(-1/\[Epsilon]-  12Lbbb+2);
	I1I2=\[Kappa] T^2/(12 \[Epsilon]b); 
	dI1I2=3\[Kappa] T^2/(12 \[Epsilon]b)-\[Kappa] T^2/6;
	I211M020=T^2/(192 \[Pi]^2);
	I2M2I2=T^2/(192 \[Pi]^2)-T^2/(384 \[Pi]^2 \[Epsilon]bbM);
	I3M2I1=T^2/(768 \[Pi]^2 \[Epsilon]bbM)+T^2/(384 \[Pi]^2);
	Clear[LbbM];
	I1p=T^2/(12 )+\[Epsilon] T^2 Lbbb;
	I3M2p=1/(64 \[Pi]^2 \[Epsilon]bp)+1/(32 \[Pi]^2);
	I2p=1/(16 \[Pi]^2 \[Epsilon]bp);
	I4M4p=1/(128 \[Pi]^2 \[Epsilon]bp)+1/(48 \[Pi]^2);
	I2M2p=-T^2/24+T^2 (((-1/24)LbbM - (-1)*1/24 Lb )+1/12)\[Epsilon];
	IF1p=-T^2/(24)+1/24 T^2 Lfff \[Epsilon];
	IF2p=1/(16 (\[Pi]^2) ) (1/\[Epsilon]+Lf);
	IF3M2p=1/(64 \[Pi]^2)(1/\[Epsilon]+Lf+2);
	IF2M2p=(-T^2/24-T^2 1/48* Lfff)\[Epsilon]+T^2/48;


	(* These integrals are given in 9408276 and 9410360 *)
	EPre=1/(24*16 \[Pi]^2) T^4*(\[Mu]/(4 \[Pi] T))^(6 \[Epsilon]) (1/\[Epsilon]+91/15+8( Derivative[1][Zeta][-1])/Zeta[-1]-2 Derivative[1][Zeta][-3]/Zeta[-3]);
	Series[EPre,{\[Epsilon],0,0}]//Normal//FullSimplify;
	M00=Series[EPre,{\[Epsilon],0,0}]//Normal;
	LInt=(T^4 Log[Glaisher])/(48 \[Pi]^2)+(T^4 Log[\[Mu]^2])/(768 \[Pi]^2)+T^4/(2304 \[Pi]^2 \[Epsilon])+(EulerGamma T^4)/(1152 \[Pi]^2)-(T^4 Log[\[Pi] T])/(384 \[Pi]^2)-(T^4 Log[2])/(192 \[Pi]^2);
	EFPre=1/(96*16 \[Pi]^2) T^4*(\[Mu]/(4 \[Pi] T))^(6 \[Epsilon]) (1/\[Epsilon]+173/30-42/5 Log[2]+8( Derivative[1][Zeta][-1])/Zeta[-1]-2 Derivative[1][Zeta][-3]/Zeta[-3]);
	EFInt=Series[EFPre,{\[Epsilon],0,0}]//Normal//FullSimplify//Expand;
	N00=EFInt;
	MF1M1=Series[-1/(192*16 \[Pi]^2) T^4*(\[Mu]/(4 \[Pi] T))^(6 \[Epsilon]) (1/\[Epsilon]+361/60+76/5 Log[2]+6EulerGamma-4( Derivative[1][Zeta][-1])/Zeta[-1]+4Derivative[1][Zeta][-3]/Zeta[-3]),{\[Epsilon],0,0}]//Normal;
	MF00=Series[-1/(192*16 \[Pi]^2) T^4*(\[Mu]/(4 \[Pi] T))^(6 \[Epsilon]) (1/\[Epsilon]+179/30-34/5 Log[2]+8( Derivative[1][Zeta][-1])/Zeta[-1]-2Derivative[1][Zeta][-3]/Zeta[-3]),{\[Epsilon],0,0}]//Normal;
	MFM2P2=Series[-29/(1728*16 \[Pi]^2) T^4*(\[Mu]/(4 \[Pi] T))^(6 \[Epsilon]) (1/\[Epsilon]+89/29-90/29 Log[2]+48/29EulerGamma+136/29( Derivative[1][Zeta][-1])/Zeta[-1]-10/29Derivative[1][Zeta][-3]/Zeta[-3]),{\[Epsilon],0,0}]//Normal;
	N2M2=Series[1/(108*16 \[Pi]^2) T^4*(\[Mu]/(4 \[Pi] T))^(6 \[Epsilon]) (1/\[Epsilon]+35/8+3/2EulerGamma-63/10 Log[2]+5( Derivative[1][Zeta][-1])/Zeta[-1]-1/2Derivative[1][Zeta][-3]/Zeta[-3]),{\[Epsilon],0,0}]//Normal;
	MP2M2=Series[11/(216*16 \[Pi]^2) T^4*(\[Mu]/(4 \[Pi] T))^(6 \[Epsilon]) (1/\[Epsilon]+73/22+12/11EulerGamma+64/11( Derivative[1][Zeta][-1])/Zeta[-1]-10/11Derivative[1][Zeta][-3]/Zeta[-3]),{\[Epsilon],0,0}]//Normal;


	(*Note that some contractions are inefficent. Will fix later but they are anyway not bottlenecks*)
	(*The names of all diagrams follow Martin's notation arXiv:1709.02397*)

	(* --- Pure scalar sector --- *)
	Temp=TensorContract[\[Lambda]4,{3,4}];
	I1Temp=-LInt;
	LSSSS=1/16  *I1Temp*Activate@TensorContract[Inactive@TensorProduct[Temp,Temp],{{1,3},{2,4}}];
	
	I1Temp=-M00;
	ESSSS=1/48I1Temp*TensorContract[\[CapitalLambda]\[Lambda],{{1,3},{2,4}}];

	(* --- Scalar vector sector --- *)
	I1Temp=1/2Activate@TensorContract[Inactive@TensorProduct[\[CapitalLambda]g,\[Lambda]4],{{1,3},{2,4},{5,6}}];
	LSSVS=I1Temp/2*D I1p^2 I2p;
	
	I1Temp= 1/2Activate@TensorContract[Inactive@TensorProduct[\[CapitalLambda]g,\[Lambda]4],{{1,3},{2,4},{5,6}}];
	JSSVSS=-1/2I1Temp I1p^2 I2p;

	(* --- Yukawa scalar sector --- *)
	Temp1=Activate@TensorContract[Inactive@TensorProduct[Ysij,Ysij],{{1,3},{2,4}}]+Activate@TensorContract[Inactive@TensorProduct[Ysij,YsijC],{{1,3},{2,4}}];
	KSSFFFF=-Temp1/8*(4*IF1p^2 I2p  + N00);
	
	Help1=Activate@TensorContract[Inactive@TensorProduct[Ysff,YsffC],{{1,4},{3,5}}];
	Temp1=Activate@TensorContract[Inactive@TensorProduct[Help1,Help1],{{1,4},{2,3}}];
	KFFSFSF=-1/2*(-1)*Temp1*(     MF1M1+    MF00+(I1p-IF1p)^2 IF2p);

	(* --- Yukawa scalar sector with Lambda4 --- *)
	help1=TensorContract[\[Lambda]4,{{3,4}}];
	help2=Activate@TensorContract[Inactive@TensorProduct[Ysff,YsffC],{{2,5},{3,6}}]//SparseArray;
	Temp1=Activate@TensorContract[Inactive@TensorProduct[help1,help2],{{1,3},{2,4}}];
	JSSFFS=1/4*2Temp1*I1p IF1p I2p;


	(* --- Yukawa vector sector --- *)
	Temp1=Activate@TensorContract[Inactive@TensorProduct[Ysij,\[CapitalLambda]g],{{1,3},{2,4}}];
	KSSSVFF=-Temp1/2*2*(MF00-I1p IF1p I2p);
	
	JSSFFV=-Temp1*D I1p IF1p I2p/.D->4-2\[Epsilon];
	
	Help1=Flatten[Ysff,{{3},{1,2}}] . Flatten[YsffC,{1,2}];
	Temp1=Total[Flatten[TensorContract[HabIJF,{1,2}]]Flatten[Help1]];
	KFFFSVF=Temp1/2(2D-4)(MF1M1+MF00+(I1p-IF1p)^2 IF2p)/.D->4-2\[Epsilon];
	
	Help1=Transpose[Transpose[gvff,{2,1,3}],{1,3,2}] . gvff;
	Help2=Transpose[Transpose[Ysff,{2,1,3}],{1,3,2}] . YsffC;
	Temp1=Total[Flatten[Help1]Flatten[Help2,{4,1,3,2}]];
	HFFSVFF=Temp1/2*((2D-4)MF00+(3-D)N00)/.D->4-2\[Epsilon];
	
	Help1=Flatten[TensorProduct[YsffC,Ysff],{{2},{3},{5},{6},{1,4}}] . Flatten[gvss,{2,3}];
	Temp1=-I Flatten[Help1,{{1,4},{5,2,3}}] . Flatten[gvff]//Tr[#]&//Simplify;
	HSSFVFF=Temp1*(MF00);


	(* --- Fermion vector sector --- *)
	(*General nF modification*)
	gvffnF=gvff . NFMat;
	(***********)
	Help1=Flatten[HabIJF,{{1},{2},{4,3}}] . Flatten[gvffnF,{2,3}];
	Temp1=6*I/3*Total[Flatten[Help1]Flatten[gvvv]];
	HFFFVVV=1/2(D-2)Temp1 MF00/.D->4-2\[Epsilon];

	(*General nF modification*)
	gvffnF=gvff . NFMat;
	(***********)
	Temp1=8*1/4Sum[Tr[gvffnF[[a]] . gvff[[b]] . gvff[[a]] . gvff[[b]]],{a,1,nv},{b,1,nv}]//Simplify//Expand;
	HFFVVFF=1/8*(D-2)*(2*(4-D)*MF00+(D-6)N00)*Temp1/.D->4-2\[Epsilon];
	GabV=TensorContract[GabcdV,{{2,4}}]//SparseArray;

	(*General nF modification*)
	HabIJFnF=HabIJF . NFMat;
	(***********)
	Temp1=-4*1/4*Tr[Flatten[GabV] . Flatten[HabIJFnF,{1,2}]];
	KGaugeFF=-1/2*Temp1*(D-2)*(2 MFM2P2+MF00+2(D-6)IF1p I2p I1p)/.D->4-2\[Epsilon];

	(*General nF modification*)
	HabIJFnF=HabIJF . NFMat;
	(***********)
	Help1=TensorContract[HabIJFnF,{3,4}]//SparseArray;
	Temp1=-1*Flatten[Help1] . Flatten[Help1];
	KVVFFFF=Temp1/4*(4N2M2+(D-4)N00-4(6-D)IF1p^2 I2p)/.D->4-2\[Epsilon];

	(*General nF modification*)
	gvffNf=gvff . NFMat;
	(***********)
	Temp1=-1/2Sum[Tr[gvffNf[[a]] . gvff[[a]] . gvff[[b]] . gvff[[b]]],{a,1,nv},{b,1,nv}]//Simplify//Expand;
	KFFFVVF=-(2-D)^2Temp1*(MF1M1+MF00+(IF1p-I1p)^2 IF2p)/.D->4-2\[Epsilon];


	(* --- Scalar vector sector --- *)
	Help1=Table[Tr[gvss[[a]] . gvss[[b]] . gvss[[c]]],{a,1,nv},{b,1,nv},{c,1,nv}]//SparseArray;
	Temp1=-1*Total[Flatten[Help1]Flatten[gvvv]];
	HSSSVVV=5/8Temp1 M00;
	
	Temp1=-Sum[Tr[gvss[[a]] . gvss[[b]] . gvss[[a]] . gvss[[b]]],{a,1,nv},{b,1,nv}];
	HSSVVSS=Temp1*5/8*M00;
	
	Temp1=-1/2*Total[Flatten[HabijV]Flatten[HabijV]];
	EVVSS=1/4 D Temp1*M00;
	
	Temp1=-1/2Total[Flatten[HabijV]Flatten[HabijV]];
	GSSVVS=-9/8*Temp1*M00;
	
	GHelp=TensorContract[GabcdV,{2,4}];
	Temp1=1/2Total[Flatten[GHelp]Flatten[Hg]];
	KGaugeSS=Temp1/2*(2*(D-2)MP2M2-(D+2)/2M00-8(D-2)I1p^2 I2p);
	
	KGaugeS=Temp1*(D-2)^2 I1p^2 I2p;
	
	Temp1=-1/4Total[Flatten[Hg]Flatten[Hg]];
	KVVSSSS=Temp1*1/4*(4MP2M2-M00+4(D-6)I1p^2 I2p);
	
	Temp1=-1/2Total[Flatten[\[CapitalLambda]g]Flatten[\[CapitalLambda]g]];
	JSSSVV=Temp1*(-D I1p^2 I2p+2 M00+1/2 I1p^2 I2p+D^2/2 I1p^2 I2p);
	
	(* --- Pure vector sector --- *)
	Help1=Flatten[GabcdV,{{1},{3},{2,4}}] . Flatten[gvvv,{2,3}];
	Temp1=-2*Total[Flatten[Help1]Flatten[gvvv]];
	HGauge=SerEnergyHelp[(1/8*(5D-5-3/4)-1/16-1/32+3/16(D-1)*D-27/16*(D-1))*M00]*Temp1;
	
	GHelp=TensorContract[GabcdV,{2,4}];
	Temp1=-Total[Flatten[GHelp]Flatten[GHelp]];
	KGauge=Temp1*SerEnergyHelp[1/4*((D-2)^2 MP2M2-((D+2)^2/4-4D)M00+(D-6)(D-2)^2 I1p^2 I2p)];
	
	KGhost=-1/8 M00*Temp1;

	(* --- Fermion vector scalar sector --- *)
	(* General nF modification *)
	HabIJFnF=HabIJF . NFMat;
	(***********)
	HFab=TensorContract[HabIJFnF,{{3,4}}];
	Temp1=2*Total[Flatten[Hg]Flatten[HFab]];
	KVVFFSS=-1/2*Temp1*(MFM2P2-1/2 MF00+(D-6)I1p IF1p I2p);


	(* --- Counterterms --- *)
	Temp=1/2*\[Gamma]ij . \[Lambda]4;
	\[Lambda]4Eff=\[Epsilon]*Z\[Lambda]ijkl;
	I1Temp=1/3 T^4 \[Epsilon] Log[Glaisher]+1/72 T^4 \[Epsilon] Log[\[Mu]^2]-1/36 T^4 \[Epsilon] Log[4 \[Pi] T]+T^4/144;
	VssZ=1/8* \[Epsilon]^-1 I1Temp*TensorContract[\[Lambda]4Eff,{{1,2},{3,4}}];

	V2SV=I1p^2 (D+\[Xi]-1)/.\[Xi]->1;
	
	V2SSV=-I1p^2 (2 \[Xi]+1)/.\[Xi]->1;
	
	V2VV=1/2 (2 (I1p (D+\[Xi]-2)+I2M2p (-\[Xi])+I2M2p) (I1p ((D+\[Xi]-2)D -2 \[Xi]+2)+D I2M2p (\[Xi]-1)))/(D-1)/.\[Xi]->1;
	
	V2VVV=-((3 (I1p^2 (2 (D-2) \[Xi]^2+ (4 D-11) \[Xi] D+ (2 D-3)D+11 \[Xi]-1)-2 D I2M2p^2 (\[Xi]-1)^2+4 I1p I2M2p (\[Xi]-1)^2))/(2 (D-1)))/.\[Xi]->1;
	
	V2ggV=2 1/4 I1p^2 (\[Xi]+1)/.\[Xi]->1;
	
	V2FFV=-(D-2) IF1p (2 I1p-IF1p)/.\[Xi]->1;
	
	V2FFS=IF1p (IF1p-2 I1p);
	
	VssvZ=1/4(-1)*2/4V2SSV*TensorContract[Zgvvss,{{1,2},{3,4}}];
	
	VvsZ=1/2*(-1)*2/4*V2SV*TensorContract[Zgvvss,{{1,2},{3,4}}];
	
	Temp=2Total[Flatten[Zgvvv]Flatten[gvvv]];
	VvvZ=1/4*\[Epsilon]^-1 V2VV*Temp;
	
	VvvvZ=1/12*\[Epsilon]^-1 V2VVV*Temp;
	
	VggvZ=1/4*\[Epsilon]^-1 V2ggV*Temp;
	
	Temp=Flatten[ZYsij,{{1},{2,3}}] . Flatten[YsffC,{2,3}]+Flatten[Ysff,{{1},{2,3}}] . Flatten[ZYsijC,{2,3}];
	VFFsZ=1/2*\[Epsilon]^-1 V2FFS*Tr[Temp];

	(*General nF modification*)
	gvffnF=gvff . NFMat;
	(*************************)
	Temp=Transpose[Zgvff . Transpose[gvffnF,{2,1,3}]+gvffnF . Transpose[Zgvff,{2,1,3}],{1,3,2,4}];
	VFFvZ=1/2*\[Epsilon]^-1*V2FFV*Tr[Flatten[Temp,{{1,3},{2,4}}]];

	(* --- Scalar mass contributions --- *)
	\[Lambda]Help=TensorContract[\[Lambda]4,{3,4}];
	V\[Mu]SS=Tr[\[Lambda]Help . \[Mu]ij]/4* I1p I2p*(-1);
	
	V\[Mu]FFS=-2*1/2 IF1p I2p*Tr[Ysij . \[Mu]ij]*(-1);
	
	V\[Mu]SV=-2/4*(D-1)*Tr[\[Mu]ij . \[CapitalLambda]g]I2p I1p*(-1);
	
	VLOtoNLOCT=-1/2Tr[\[Beta]mij]/2*I1p*\[Epsilon]^-1*(-1);
	
	VLOtoNNLO=-1/4Tr[\[Mu]ij . \[Mu]ij]*1/(16 (\[Pi]^2) ) Lb*(-1);

	ContriMass=VLOtoNNLO+ VLOtoNLOCT+V\[Mu]SS+V\[Mu]FFS+ V\[Mu]SV;

	(* --- Fermion mass contributions --- *)
	VLOtoNNLOF=1/2*Tr[\[Mu]IJ . \[Mu]IJC . \[Mu]IJ . \[Mu]IJC]*1/(16 (\[Pi]^2) ) Lf;

	(* Total NNLO symmetric pressure contributions, organized by sector *)
	(* Counterterm diagrams *)
	DiagramsCT = SerEnergyHelp[
	    VssvZ + VvsZ + VssZ + VvvZ + VvvvZ + VggvZ + VFFsZ + VFFvZ
	] // Simplify;
	
	(* Pure scalar sector (\[Lambda]-dependent) *)
	DiagramsScalar\[Lambda] = SerEnergyHelp[LSSSS + ESSSS] // Simplify;
	
	(* Scalar-vector mixing (\[Lambda]-dependent) *)
	DiagramsScalarVector\[Lambda] = SerEnergyHelp[LSSVS + JSSVSS] // Simplify;
	
	(* Yukawa: scalar couplings *)
	DiagramsYukawaScalar = SerEnergyHelp[KSSFFFF + KFFSFSF] // Simplify;
	
	(* Yukawa: \[Lambda]-dependent scalar-fermion *)
	DiagramsYukawaScalar\[Lambda] = SerEnergyHelp[JSSFFS] // Simplify;
	
	(* Yukawa: vector couplings *)
	DiagramsYukawaVector = SerEnergyHelp[
	    KSSSVFF + JSSFFV + KFFFSVF + HFFSVFF + HSSFVFF
	] // Simplify;
	
	(* Fermion-vector interactions *)
	DiagramsFermionVector = SerEnergyHelp[
	    HFFFVVV + HFFVVFF + KGaugeFF + KVVFFFF + KFFFVVF
	] // Simplify;
	
	(* Scalar-vector diagrams not tied to \[Lambda] *)
	DiagramsScalarVector = SerEnergyHelp[
	    HSSSVVV + EVVSS + GSSVVS + KGaugeSS + KGaugeS + 
	    KVVSSSS + JSSSVV + HSSVVSS
	] // Simplify;
	
	(* Pure gauge and ghost contributions *)
	DiagramsPureVector = SerEnergyHelp[
	    HGauge + KGauge + KGhost
	] // Simplify;
	
	(* Mixed fermion-scalar-vector contribution *)
	DiagramsFermionScalarVector = SerEnergyHelp[KVVFFSS] // Simplify;
	
	(* One-loop and mass-related term *)
	DiagramsFermionMass = VLOtoNNLOF;
	
	(* Total NNLO contribution *)
	DiagramsTotal = (
	      DiagramsFermionMass
	    + ContriMass
	    + DiagramsCT
	    + DiagramsFermionScalarVector
	    + DiagramsPureVector
	    + DiagramsScalarVector
	    + DiagramsFermionVector
	    + DiagramsYukawaVector
	    + DiagramsYukawaScalar\[Lambda]
	    + DiagramsScalar\[Lambda]
	    + DiagramsScalarVector\[Lambda]
	    + DiagramsYukawaScalar
	);

	(* \[Epsilon]-expansion and simplification *)
	DiagramsSimplified=Series[
		DiagramsTotal/.D->4-2\[Epsilon] 
		/.\[Epsilon]bp->(1/\[Epsilon]+Lb)^-1
		/.\[Epsilon]b->(1/\[Epsilon]+Lbb)^-1
		/.\[Epsilon]BF->(1/\[Epsilon]+LBF)^-1 
		/.\[Epsilon]F->(1/\[Epsilon]+LFF)^-1
		/.\[Epsilon]FB->(1/\[Epsilon]+LFB)^-1
		/.\[Epsilon]bbM->(1/\[Epsilon]+LbbM)^-1
		/.ReplaceLb,
		{\[Epsilon],0,0}]//Normal;

	(* Extract \[Epsilon]^0 coefficient and clean namespace *)
	ToExpression[
		StringReplace[
			ToString[StandardForm[Coefficient[Simplify[DiagramsSimplified],\[Epsilon],0]]],
			"DRalgo`Private`"->""]]
];


(* ::Section::Closed:: *)
(*Scalar masses*)


(* ::Subsection:: *)
(*1loop scalar mass*)


(* 
    Computes the 1-loop scalar mass correction (aS3D) in the dimensionally reduced 
    soft theory, incorporating contributions from scalar self-interactions, 
    gauge boson loops, fermion self-energies, and optionally 6-point scalar couplings.
*)
ScalarMass[] := Module[
  {
    selfEnergyFF, contriSS, contriVV, contriFF,
    contriSSLambda6, result
  },

  VerbosePrint["Calculating 1-Loop Scalar Mass"];

  (* Scalar 4-point interaction contribution *)
  contriSS = -T^2/24 TensorContract[\[Lambda]4, {{3, 4}}] // SimplifySparse;

  (* Gauge boson loop contribution *)
  contriVV = T^2/4 \[CapitalLambda]g;

  (* Fermion self-energy contribution *)
  selfEnergyFF = -T^2/12;
  contriFF = 1/2 selfEnergyFF (Ysij + YsijC) // SimplifySparse;

  (* If mode \[GreaterEqual] 3, include 6-point scalar coupling contribution *)
  If[mode >= 3,
    contriSSLambda6 = T^4/1152 TensorContract[\[Lambda]6, {{1, 2}, {3, 4}}];
    result = \[Mu]ij - contriSS - contriVV - contriFF + contriSSLambda6,
    
    (* Otherwise, only include up to 4-point interactions *)
    (* Minus signs from the matching *)
    result = \[Mu]ij - contriSS - contriVV - contriFF
  ];

  (* Store the final corrected scalar mass, fully simplified *)
  aS3D = result // Normal // FullSimplify // Expand;
];


(* ::Subsection:: *)
(*1loop scalar self-energy*)


(*
  Computes the 1-loop scalar self-energy (wavefunction renormalization)
  in the dimensionally reduced soft theory.
*)
ScalarSelfEnergy[] := Module[
	{
		ContriVV, SelfEnergyFF, ContriFF
	},
	
  VerbosePrint["Calculating Scalar-Field Renormalization"];

  (* Gauge boson loop contribution *)
  ContriVV = -3/(16 \[Pi]^2) Lb \[CapitalLambda]g;

  (* Fermion loop contribution via Yukawa couplings *)
  SelfEnergyFF = -Lf/(16 \[Pi]^2);
  ContriFF = 1/2 SelfEnergyFF (Ysij + YsijC);

  (* Total wavefunction renormalization from gauge + fermion loops *)
  ZijS = -ContriVV/2 - ContriFF/2 // SimplifySparse;
];


(* ::Subsection:: *)
(*2loop scalar mass*)


(* Calculates the scalar mass to 2 loops in the soft theory *)

ScalarMass2Loop[] := Module[
	{
		I1I2, dI1I2, I1, I2I3m1, I1FI2F, I1BI2F, I1FI2B, 
		I1FI2FD, I1BI2FD, I1FI2BT1, I1F, Contri1, Contri2, Contri21, Contri3, 
		Contri31, Contri32, Contri4, Contri5, Contri6, Contri61, Contri62, Contri7, 
		Contri8, Contri9, help, Contri10, Contri11, Contri12, Contri13, Contri14, 
		Contri141, Contri15, Contri16, ContriF, ContriCubic, Fac, Temp1, Temp1C,
		simplifyContract
	},
  
  VerbosePrint["Calculating 2-Loop Scalar Mass"];
  
  (* Constants *)
  \[Kappa] = 1/(16 \[Pi]^2);
  I1I2 = \[Kappa] T^2/(12 \[Epsilon]b);
  dI1I2 = 3 \[Kappa] T^2/(12 \[Epsilon]b) - \[Kappa] T^2/6;
  I1 = T^2/(12 \[Epsilon]) + T^2 Lbbb;
  I2I3m1 = T^2/(512 \[Pi]^2 \[Epsilon]b) - T^2/(768 \[Pi]^2);
  I1FI2F = -((T^2)/(384 \[Pi]^2))/\[Epsilon]F;
  I1BI2F = ((T^2)/(192 \[Pi]^2))/\[Epsilon]FB;
  I1FI2B = -((T^2)/(384 \[Pi]^2))/\[Epsilon]BF;
  I1FI2FD = (3 - 2 \[Epsilon]F) I1FI2F;
  I1BI2FD = (3 - 2 \[Epsilon]FB) I1BI2F;
  I1FI2BT1 = I1FI2B(-4) - T^2/(96 \[Pi]^2);
  I1F = -T^2/(24 \[Epsilon]) + 1/24 T^2 Lfff;

  (* Helper Functions *)
  simplifyContract[expr_] := Simplify[expr] // SimplifySparse;

  (* Contributions *)
  Contri1 = (1/4) I1I2 Simplify[TensorContract[\[CapitalLambda]\[Lambda], {3, 4}]];
  Contri2 = (1/2) I1 Flatten[\[Gamma]ij] . Flatten[\[Lambda]4, {1, 2}];
  Contri21 = (1/(16 \[Pi]^2)) Lb * (1/2) Flatten[\[Mu]ij] . Flatten[\[Lambda]4, {1, 2}];
  Contri3 = -(1/2) I1 \[Epsilon] TensorContract[Z\[Lambda]ijkl, {3, 4}];
  Contri31 = -(1/2) I1 Flatten[\[Gamma]ij] . Flatten[\[Lambda]4, {1, 2}];
  Contri32 = -(1/4) I1*(
    + \[Gamma]ij . TensorContract[\[Lambda]4, {1, 2}]
    + TensorContract[\[Lambda]4, {1, 2}] . \[Gamma]ij
    );
  Contri4 = -(1/2) dI1I2 Flatten[\[CapitalLambda]g] . Flatten[\[Lambda]4, {1, 2}];
  Contri5 = (1/2) dI1I2 Contract[Hg, HabijV, {{1, 3}, {2, 4}}] // simplifyContract;
  Contri6 = (3 - 2 \[Epsilon]) I1 \[Epsilon] (1/2) TensorContract[Zgvvss, {1, 2}];
  Contri61 = (3 - 2 \[Epsilon]) I1 (1/2) Contract[\[Gamma]ab, HabijV, {{1, 3}, {2, 4}}] // simplifyContract;
  Contri62 = (3 - 2 \[Epsilon]) I1 (1/2) (
     + Contract[\[Gamma]ij, Habij, {{3, 4}, {2, 5}}]
     + Contract[\[Gamma]ij, Habij, {{3, 4}, {2, 6}}]) // simplifyContract;
  Contri7 = -(1/2) I1 (3 - 2 \[Epsilon]) Contract[\[Gamma]ab, HabijV, {{1, 3}, {2, 4}}] // simplifyContract;
  Contri8 = -(1/2) I1I2 Contract[Hg, HabijV, {{1, 3}, {2, 4}}] // simplifyContract;
  Contri9 = T^2 (1/4) (-1) (9/(8 \[Epsilon]b) - 13/8) (1/(16 \[Pi]^2)) Contract[GabcdV, HabijV, {{2, 4}, {1, 5}, {3, 6}}] // simplifyContract;
  
  help = Contract[GabcdV, HabijV, {{1, 5}, {2, 3}, {4, 6}}] // simplifyContract;
  Contri10 = (I1I2/4) help;
  
  Contri11 = T^2 (-1) (1/4) (13/(24 \[Epsilon]b) - 7/24)/(16 \[Pi]^2) help;
  
  Contri12 = (-1/2) (-1) 2 (I1FI2F - I1BI2F)*(
    + Flatten[TensorContract[YTemp, {1, 2}]] . Flatten[YTemp, {4, 3}]
    + Flatten[TensorContract[YTempC, {1, 2}]] . Flatten[YTempC, {4, 3}]
    );
  
  Contri13 = (-2 I1FI2B) (1/4)*(
    + Flatten[Ysij] . Flatten[\[Lambda]4, {1, 2}]
    + Flatten[YsijC] . Flatten[\[Lambda]4, {1, 2}]
    );

  (* Kosm or "Kos" contributions *)
  help1 = (
    + Flatten[Ysff, {{1}, {2, 3}}] . Flatten[ZYsijC, {2, 3}]
    + Flatten[YsffC, {{1}, {2, 3}}] . Flatten[ZYsij, {2, 3}]
    );
  Contri14 = I1F (help1 + Transpose[help1, {2, 1}]) // Simplify;
  
  help1 = Ysij . \[Gamma]ij;
  Contri141 = I1F (help1 + Transpose[help1]) // Simplify;
  
  (* General nF modification *)
  help = TensorContract[HabIJF . NFMat, {3, 4}];
  
  (* Contributions from self-energy, cubics, and fermion masses *)
  Contri15 = 2 I1FI2BT1 (-1) (1/4) (Flatten[help] . Flatten[HabijV, {1, 2}]);
  
  help = TensorContract[HabIJF, {1, 2}];
  help2 = (
    + Flatten[help] . Flatten[YTemp, {3, 4}]
    + Flatten[help] . Flatten[YTempC, {3, 4}]
    );
  Contri16 = (-1) (I1FI2FD - I1BI2FD) (-1) (2/3) (1 - 1 \[Epsilon]/3) help2;

  (* Self-energy contribution *)
  ContriF = -(ZijS . aS3D + Transpose[ZijS . aS3D]);

  (* Contribution from Cubics *)
  ContriCubic = (1/2) (1/(16 \[Pi]^2)) Lb Contract[\[Lambda]3, \[Lambda]3, {{2, 5}, {3, 6}}] // simplifyContract;

  (* Contribution from Fermion Masses *)
  Fac = 1/(8 \[Pi]^2) Lf;
  Temp1 = Activate @ TensorContract[Inactive @ TensorProduct[Ysff, \[Mu]IJFC, \[Mu]IJFC, Ysff], {{2, 4}, {3, 6}, {5, 9}, {7, 10}}] // Normal;
  Temp1C = Activate @ TensorContract[Inactive @ TensorProduct[YsffC, \[Mu]IJF, \[Mu]IJF, YsffC], {{2, 4}, {3, 6}, {5, 9}, {7, 10}}] // Normal;
  ContriFFMass1 = -(1/2) Fac (Temp1 + Temp1C) // Simplify // PowerExpand;
  
  Temp1 = Activate @ TensorContract[Inactive @ TensorProduct[Ysff, YsffC, \[Mu]IJFC, \[Mu]IJF], {{3, 5}, {2, 7}, {6, 10}, {8, 9}}] // Normal;
  Temp1C = Activate @ TensorContract[Inactive @ TensorProduct[YsffC, Ysff, \[Mu]IJF, \[Mu]IJFC], {{3, 5}, {2, 7}, {6, 10}, {8, 9}}] // Normal;
  ContriFFMass2 = -Fac (Temp1 + Temp1C) // Simplify // PowerExpand;

  (* Final scalar mass *)
  \[Mu]SijNLO = (
    + SerEnergyHelp[ContriF // Normal]
    - SerEnergyHelp[ContriCubic // Normal]
    - SerEnergyHelp[(ContriFFMass1 + ContriFFMass2) // Normal]
    - SerEnergyHelp[(Contri1 + Contri2 + Contri3 + Contri31 + Contri32 + Contri4 + Contri5) // Normal]
    - SerEnergyHelp[(Contri6 + Contri61 + Contri62 + Contri7 + Contri8 + Contri9 + Contri10 + Contri11) // Normal]
    - SerEnergyHelp[(Contri21 + Contri12 + Contri13 + Contri14 + Contri141 + Contri15 + Contri16) // Normal]
    ) // Simplify;
];


(* ::Section::Closed:: *)
(*Temporal masses*)


(* ::Subsection:: *)
(*1loop temporal mass*)


(* 
    1-loop Debye mass in the soft theory.
*)
VectorMass[] := Module[
  {
    ContriSS, ContriVV, fac, ContriVVV, SelfEnergyFF, HabIJFnF, ContriFF
  },

  VerbosePrint["Calculating 1-Loop Vector Mass"];

  (* Scalar and vector boson contributions *)
  ContriSS = (T^2/12) Hg;
  ContriVV = (T^2/12) Hg;

  (* Pure gauge 3-vertex contribution *)
  fac = Simplify[T^2/24 - T^2/4 - (1/2) T^2/4];
  ContriVVV = fac TensorContract[GabcdV, {{2, 4}}];

  (* Fermion contribution *)
  SelfEnergyFF = -(T^2/6);
  HabIJFnF = HabIJF . NFMat;  (* Incorporate general nF modification *)
  ContriFF = SelfEnergyFF Simplify[TensorContract[HabIJFnF, {{3, 4}}] // Normal];

  (* Result: Minus sign due to matching convention *)
  aV3D = -ContriSS - ContriVV - ContriVVV - ContriFF // Normal // FullSimplify // Expand;
];


(* ::Subsection:: *)
(*1loop vector self-energy*)


(* 
    Vector self-energy in the soft theory.
*)
VectorSelfEnergy[] := Module[
  {
    SelfEnergySS, ContriSS, fac, ContriVVV, SelfEnergyFF,
    HabIJFnF, ContriFF
  },

  VerbosePrint["Calculating Vector-Field Renormalization"];

  (* Contribution from scalar self-energy (transverse part) *)
  SelfEnergySS = -1/2 * 1/(16 \[Pi]^2) * (-1/3 Lb);
  ContriSS = SelfEnergySS * Hg;

  (* Vector self-energy contribution (transverse part) *)
  fac = 1/(16 \[Pi]^2) * (1/12 Lb + 1/2 (25/6 Lb + 2/3));
  ContriVVV = fac * Simplify[TensorContract[GabcdV, {{2, 4}}]];

  (* Fermion self-energy contribution *)
  SelfEnergyFF = - (2/3 Lf)/(16 \[Pi]^2);
  HabIJFnF = HabIJF . NFMat;  (* General nF modification *)
  ContriFF = SelfEnergyFF * Simplify[TensorContract[HabIJFnF, {{3, 4}}] // Normal];

  (* Transverse vector contribution *)
  ZabT = -(ContriSS + ContriVVV + ContriFF)/2 // Normal // Simplify;

  (* Contribution from scalar self-energy (longitudinal part) *)
  SelfEnergySS = -1/2 * 1/(16 \[Pi]^2) * (-1) * (1/3 Lb + 2/3);
  ContriSS = SelfEnergySS * Hg;

  (* Vector self-energy contribution (longitudinal part) *)
  fac = 1/(16 \[Pi]^2) * (1/12 Lb + 1/6 + 1/2 (25/6 Lb - 3));
  ContriVVV = fac * Simplify[TensorContract[GabcdV, {{2, 4}}]];

  (* Fermion self-energy contribution (longitudinal part) *)
  SelfEnergyFF = - (2/3 Lf - 2/3)/(16 \[Pi]^2);
  HabIJFnF = HabIJF . NFMat;  (* General nF modification *)
  ContriFF = SelfEnergyFF * Simplify[TensorContract[HabIJFnF, {{3, 4}}] // Normal];

  (* Longitudinal vector contribution *)
  ZabL = -(ContriSS + ContriVVV + ContriFF)/2 // Normal // Simplify;

];


(* ::Subsection:: *)
(*2loop temporal mass*)


(* 
    Calculations of the 2-loop Debye mass in the soft theory.
*)
VectorMass2Loop[] := Module[
  {
    gvffnF, TrSS, TrFF, TrVV, Gabij, Gabcd, GabIJ,
    GroupFacVVVV, GroupTheoryFacm2, GroupTheoryFac\[Lambda], 
    GroupTheoryFacTrSSTrSS, GroupTheoryFacSSSS, GroupTheoryFacVVTrSS, 
    GroupTheoryFacTrVVTrSS, GroupTheoryFacVVTrFF, GroupTheoryFacTrVVTrFF, 
    GroupTheoryFacFFFF, GroupFacTrSSTrFF, GroupTheoryYuk1, GroupTheoryYuk2, 
    \[CapitalPi]V, \[CapitalPi]S, \[CapitalPi]F, \[CapitalPi]SF
  },

  VerbosePrint["Calculating 2-Loop Debye Mass"];

  (* General nF modification *)
  gvffnF = gvff . NFMat;

  (* Compute trace terms *)
  TrSS = Table[Tr[a . b], {a, gvss}, {b, gvss}] // SparseArray;
  TrFF = Table[Tr[a . b], {a, gvff}, {b, gvffnF}] // SparseArray;
  TrVV = Table[Tr[a . b], {a, gvvv}, {b, gvvv}] // SparseArray;

  (* Compute Gab terms *)
  Gabij = Table[a . b, {a, gvss}, {b, gvss}] // SparseArray;
  Gabcd = gvvv . gvvv // SparseArray;
  GabIJ = Table[a . b, {a, gvff}, {b, gvff}] // SparseArray;

  (* Group theory factors *)
  GroupFacVVVV = -Flatten[Gabcd, {{1}, {2, 3, 4}}] . Flatten[Gabcd, {1, 3, 2}];
  GroupTheoryFacm2 = Flatten[\[Mu]ij] . Flatten[Gabij, {3, 4}];
  GroupTheoryFac\[Lambda] = Flatten[TensorContract[\[Lambda]4, {{1, 2}}]] . Flatten[Gabij, {3, 4}];

  (* Various contributions based on traces and group theory *)
  GroupTheoryFacTrSSTrSS = TrSS . TrSS;
  GroupTheoryFacSSSS = Table[Sum[Tr[c . d . a . a], {a, gvss}], {c, gvss}, {d, gvss}];
  GroupTheoryFacVVTrSS = Flatten[TrSS] . Flatten[Gabcd, {2, 4}];
  GroupTheoryFacTrVVTrSS = TrSS . TrVV + TrVV . TrSS;
  GroupTheoryFacVVTrFF = Flatten[TrFF] . Flatten[Gabcd, {2, 4}];
  GroupTheoryFacTrVVTrFF = TrFF . TrVV + TrVV . TrFF;
  GroupTheoryFacFFFF = Table[Sum[Tr[a . a . c . d], {a, gvff}], {c, gvff}, {d, gvffnF}] // SparseArray;
  GroupFacTrSSTrFF = (TrSS . TrFF + TrFF . TrSS);
  GroupTheoryFacTrFFTrFF = TrFF . TrFF;

  (* Yukawa terms *)
  GroupTheoryYuk1 = Flatten[TensorContract[TensorProduct[Ysff, YsffC], {{1, 4}, {3, 5}}]] . Flatten[GabIJ, {3, 4}];
  GroupTheoryYuk2 = Flatten[Table[Tr[a . b + b . a], {a, Ysff}, {b, YsffC}] // SparseArray] . Flatten[Gabij, {3, 4}];

  (* Contributions to the Debye mass *)
  \[CapitalPi]V = GroupFacVVVV * (
    T^2 * (11 Log[\[Mu]/(4 \[Pi] T)] + 11 EulerGamma - 1/2) / (36 \[Pi]^2) + T^2 / (12 \[Pi]^2)
  );

  \[CapitalPi]S = -GroupTheoryFac\[Lambda] * T^2 / (192 \[Pi]^2) 
    - GroupTheoryFacm2 / (8 \[Pi]^2) 
    - GroupTheoryFacTrSSTrSS * (T^2 * (Log[\[Mu]/(4 \[Pi])] - Log[T] + EulerGamma + 1)) / (288 \[Pi]^2)
    + GroupTheoryFacSSSS * T^2 / (32 \[Pi]^2)
    - (T^2 * (Log[\[Mu]] - Log[T] + EulerGamma - Log[\[Pi]] - 2 Log[2])) / (24 \[Pi]^2) * GroupTheoryFacVVTrSS
    - GroupTheoryFacVVTrSS * T^2 / (48 \[Pi]^2)
    + (T^2 * (8 Log[\[Mu]] - 8 Log[T] + 8 EulerGamma - 3 - 8 Log[\[Pi]] - 16 Log[2])) / (576 \[Pi]^2) * GroupTheoryFacTrVVTrSS;

  \[CapitalPi]F = (T^2 * (Log[\[Mu]] - Log[T] + EulerGamma - Log[\[Pi]] - 2 Log[2])) / (24 \[Pi]^2) * GroupTheoryFacVVTrFF
    - (T^2 * (2 Log[\[Mu]] - 2 Log[T] + 2 EulerGamma - 1 - 2 Log[\[Pi]])) / (144 \[Pi]^2) * GroupTheoryFacTrFFTrFF
    - (T^2 * (2 Log[\[Mu]] - 2 Log[T] + 2 EulerGamma + 3 - 2 Log[\[Pi]] - 20 Log[2])) / (576 \[Pi]^2) * GroupTheoryFacTrVVTrFF
    - GroupTheoryFacFFFF * T^2 / (16 \[Pi]^2)
    + GroupTheoryFacVVTrFF * T^2 / (48 \[Pi]^2);

  \[CapitalPi]SF = (T^2 * (5 Log[\[Mu]/(4 \[Pi] T)] + 5 EulerGamma - 1 + 8 Log[2])) / (576 \[Pi]^2) * GroupFacTrSSTrFF
    - 1/(32 \[Pi]^2) * T^2 * GroupTheoryYuk1
    - T^2 / (192 \[Pi]^2) * GroupTheoryYuk2;

  (* Final expression for the 2-loop Debye mass *)
  \[Mu]VabNLO = \[CapitalPi]V + \[CapitalPi]S + \[CapitalPi]F + \[CapitalPi]SF // Normal // Simplify // FullSimplify;

];


(* ::Section::Closed:: *)
(*Effective couplings*)


(* ::Subsection::Closed:: *)
(*Gauge couplings*)


(* 
    Calculates non-abelian couplings from ghost-renormalization.
*)
NonAbelianCoupling[] := Module[
  {
    fac, ContriVVV, Zab\[Eta], ContriAnomVV
  },

  (* Define factor for the contribution *)
  fac = 3/4 * Lb / (16 \[Pi]^2);

  (* Contribution from the vector interactions *)
  ContriVVV = fac * Table[-Tr[a . b], {a, gvvv}, {b, gvvv}] // SparseArray // SimplifySparse;

  (* Calculate Zab\[Eta] term *)
  Zab\[Eta] = -1/2 * ContriVVV;

  (* Anomalous dimension contribution to the vector interaction *)
  ContriAnomVV = gvvv . ZabT + Zab\[Eta] . gvvv + Table[Zab\[Eta] . a, {a, gvvv}] // SparseArray;

  (* Final result *)
  Ggvvv = -ContriAnomVV;
];


(* ::Subsection:: *)
(*SSS*)


(* 
    Calculates 1-loop scalar cubics in the soft theory.
*)
ScalarCubic[] := Module[
  {
    SelfEnergySSC, ContriSSCTemp, ContriSSC,
    Contri1, Contri2, ContriFF,
    ContriSETemp, ContriSE
  },

  VerbosePrint["Calculating Scalar-Cubic Couplings"];

  (* Scalar self-energy contribution *)
  SelfEnergySSC = 1 / (16 \[Pi]^2) * Lb * 1 / 2;
  ContriSSCTemp = Contract[\[Lambda]4, \[Lambda]3, {{3, 5}, {4, 6}}];
  ContriSSC = SelfEnergySSC * (
    ContriSSCTemp + Transpose[ContriSSCTemp, {1, 3, 2}] + Transpose[ContriSSCTemp, {3, 2, 1}]
  );

  (* Fermion mass contributions *)
  Contri1 = Table[Tr[a . \[Mu]IJFC . b . c], {a, Ysff}, {b, Ysff}, {c, YsffC}];
  Contri2 = Table[Tr[a . \[Mu]IJF . b . c], {a, YsffC}, {b, YsffC}, {c, Ysff}];
  ContriFF = -1 / (8 \[Pi]^2) * Lf * 3 * Symmetrize[(Contri1 + Contri2), Symmetric[{1, 2, 3}]];

  (* Self-energy contribution *)
  ContriSETemp = -ZijS . \[Lambda]3;
  ContriSE = ContriSETemp + Transpose[ContriSETemp, {2, 1, 3}] + Transpose[ContriSETemp, {3, 1, 2}] // Simplify;

  (* Final result for \[Lambda]3CS *)
  \[Lambda]3CS = -ContriSSC - ContriFF + ContriSE // Simplify;
];


(* ::Subsection:: *)
(*VVVV*)


(*
    Calculates Temporal-Vector Quartics ~ (V^4). Note that all terms are finite.
*)
LongitudinalVVVV[] := Module[
  {
    CouplingSSSS, ContriSSSSTemp, Help, ContriSSSS, CouplingVV, ContriVVTemp, ContriVVTemp2, 
    ContriVV, HabIJFnF, helpF, ContriFFFF
  },

  (* Verbose output for debugging *)
  VerbosePrint["Calculating Temporal-Vector Quartics"];

  (* Temporal-Scalar Scalar contributions (Sum of Bubbles, Triangles, and Boxes with internal scalars) *)
  CouplingSSSS = -(1 / (24 \[Pi]^2));
  ContriSSSSTemp = Flatten[HabijV, {{1}, {2}, {3, 4}}] . Flatten[HabijV, {3, 4}];
  Help = Simplify[
    ContriSSSSTemp + 
    Transpose[ContriSSSSTemp, {1, 3, 2, 4}] + 
    Transpose[ContriSSSSTemp, {1, 4, 3, 2}]
  ] // SparseArray // SimplifySparse;
  ContriSSSS = CouplingSSSS * Help // SparseArray;

  (* Vector-Vector contributions (General nF modification) *)
  CouplingVV = -(1 / (6 \[Pi]^2));
  ContriVVTemp = Flatten[GabcdV, {{1}, {3}, {2, 4}}] . Flatten[GabcdV, {2, 4}];
  ContriVVTemp2 = ContriVVTemp + Transpose[ContriVVTemp, {1, 2, 4, 3}] // Simplify;
  Help = Simplify[
    (ContriVVTemp2 
    + Transpose[ContriVVTemp2, {1, 3, 2, 4}] 
    + Transpose[ContriVVTemp2, {1, 4, 3, 2}])
  ] // SparseArray // SimplifySparse;
  ContriVV = CouplingVV * Help // SparseArray;

  (* General nF modification (adjusting fermionic contributions) *)
  HabIJFnF = HabIJF . NFMat;

  (* Fermion-Fermion contributions (Sum of Boxes with Internal Fermions) *)
  helpF = Flatten[HabIJFnF, {{1}, {2}, {3, 4}}] . Flatten[Transpose[HabIJF, {1, 2, 4, 3}], {3, 4}] // SimplifySparse;
  ContriFFFF = 1 / (3 \[Pi]^2) * (
    + helpF 
    + Transpose[helpF, {1, 3, 2, 4}] 
    + Transpose[helpF, {1, 2, 4, 3}]) // SparseArray; (* Check this one after modification *)

  (* Final result with minus sign for matching *)
  \[Lambda]AA = -(ContriSSSS + ContriVV + ContriFFFF) // SimplifySparse // SparseArray;
];


(* ::Subsection:: *)
(*VVS*)


(* 
    Calculates cubics between two temporal-vectors and one scalar.
    No tree-level contribution.
*)
LongitudinalVVS[] := Module[
  {
    CouplingSSSS, ContriSSSSTemp
  },

  VerbosePrint["Calculating Temporal-Vector-Scalar Cubics"];

  (* Cubic contribution with internal scalars *)
  CouplingSSSS = 1 / (8 \[Pi]^2);
  ContriSSSSTemp = CouplingSSSS * Contract[Habij, \[Lambda]3, {{3, 5}, {4, 6}}];

  (* Minus sign from matching *)
  GvvsL = -ContriSSSSTemp;
];


(* ::Subsection:: *)
(*SSSS*)


(* 
    Calculates Scalar Quartics in the soft theory.
*)
ScalarQuartic[] := Module[
  {
    ContriSS, ContriVVTemp, ContriVVTemp2, ContriVV, 
    CouplingFF, ContriFF, ContriSETemp, ContriSE, ContriSS\[Lambda]6
  },

  (* Verbose output for debugging *)
  VerbosePrint["Calculating Scalar Quartic"];

  (* Scalar contribution to the quartic terms *)
  ContriSS = (1 / (16 \[Pi]^2)) * Lb * 1 / 2 * Simplify[
    (\[CapitalLambda]\[Lambda] + Transpose[\[CapitalLambda]\[Lambda], {1, 4, 3, 2}] + 
     Transpose[\[CapitalLambda]\[Lambda], {1, 3, 2, 4}])
  ];

  (* Temporal-vector contribution *)
  ContriVVTemp = Transpose[Flatten[HabijV, {1, 2}], {3, 2, 1}] . Flatten[HabijV, {1, 2}] // SimplifySparse;
  ContriVVTemp2 = ContriVVTemp + Transpose[ContriVVTemp, {1, 3, 2, 4}] + Transpose[ContriVVTemp, {1, 4, 2, 3}];
  ContriVV = (1 / (16 \[Pi]^2)) * (3 Lb - 2) / 2 * ContriVVTemp2;

  (* Fermion contribution *)
  CouplingFF = 2 * 1 / (16 \[Pi]^2) * (2 Lf) * (-1) * 1 / 4;
  ContriFF = CouplingFF * 12 * Symmetrize[Yhelp, Symmetric[{1, 2, 3, 4}]] // SparseArray // SimplifySparse;

  (* Self-energy contribution *)
  ContriSETemp = -ZijS . \[Lambda]4;
  ContriSE = ContriSETemp + Transpose[ContriSETemp, {2, 1, 3, 4}] + Transpose[ContriSETemp, {3, 1, 2, 4}] + 
    Transpose[ContriSETemp, {4, 1, 2, 3}] // Simplify;

  (* Conditional handling based on mode *)
  If[mode >= 3,
    (* Additional contribution based on the mode *)
    ContriSS\[Lambda]6 = -T^2 / 24 * TensorContract[\[Lambda]6, {1, 2}];
    \[Lambda]3D = -ContriSS - ContriVV + ContriSE - ContriFF - ContriSS\[Lambda]6,
    
    (* Default contribution *)
    \[Lambda]3D = -ContriSS - ContriVV + ContriSE - ContriFF;
  ];
];


(* ::Subsection:: *)
(*SSSSSS*)


(*
    Calculates Scalar Sextic contributions in the soft theory
*)
ScalarSextic[] := Module[
  {
    LambdaLambda6Temp, LambdaLambda6Tot, Prefac, ContriSS, SymHelp, 
    ContriSSS, ContriVV, ContriFF, ContriSETemp, ContriSE, LambdaLambda6Tot1, 
    LambdaLambda6TotC, LambdaLambda6TotFermion
  },

  (* Verbose output for debugging *)
  VerbosePrint["Calculating Scalar Sextic"];

  (* Scalar loop contributions *)
  LambdaLambda6Temp = Flatten[\[Lambda]4 . \[Lambda]4, {{1}, {2}, {4}, {5}, {3, 6}}];
  LambdaLambda6Tot = LambdaLambda6Temp . Flatten[\[Lambda]4, {1, 2}];
  Prefac = -Zeta[3] * 15 / (128 \[Pi]^4 T^2);
  ContriSS = Prefac * Symmetrize[LambdaLambda6Tot, Symmetric] // SparseArray // SimplifySparse;

  (* Scalar loop with mixed \[Lambda]6 and \[Lambda]4 vertices *)
  LambdaLambda6Tot = Flatten[\[Lambda]4, {{1}, {2}, {3, 4}}] . Flatten[\[Lambda]6, {1, 2}] // SimplifySparse;
  Prefac = 15 / 2 * 1 / (16 \[Pi]^2) * Lb;
  SymHelp = Symmetrize[LambdaLambda6Tot, Symmetric] // SparseArray;
  ContriSSS = Prefac * SymHelp;

  (* Vector loop contributions *)
  LambdaLambda6Temp = Flatten[Transpose[HabijV, {1, 4, 3, 2}] . HabijV, {{2}, {3}, {5}, {6}, {1, 4}}];
  LambdaLambda6Tot = LambdaLambda6Temp . Flatten[HabijV, {{1, 2}, {3}, {4}}];
  Prefac = 3 * 15 * Zeta[3] / (128 \[Pi]^4 T^2);
  ContriVV = Prefac * Symmetrize[LambdaLambda6Tot, Symmetric] // SparseArray // SimplifySparse;

  (* Fermion loop contributions *)
  LambdaLambda6Tot1 = Table[Tr[a . b . c . d . e . f], {a, Ysff}, {b, YsffC}, {c, Ysff}, {d, YsffC}, {e, Ysff}, {f, YsffC}] // SparseArray;
  LambdaLambda6TotC = Table[Tr[a . b . c . d . e . f], {a, YsffC}, {b, Ysff}, {c, YsffC}, {d, Ysff}, {e, YsffC}, {f, Ysff}] // SparseArray;
  LambdaLambda6TotFermion = LambdaLambda6Tot1 + LambdaLambda6TotC;
  Prefac = ((7 * Zeta[3]) * 15 * 4) / (64 \[Pi]^4 T^2);
  ContriFF = Prefac * Symmetrize[LambdaLambda6TotFermion, Symmetric] // SparseArray // SimplifySparse;

  (* Field-strength renormalization *)
  ContriSETemp = -ZijS . \[Lambda]6;
  ContriSE = 6 * Symmetrize[ContriSETemp, Symmetric] // SparseArray // SimplifySparse;

  (* Minus sign from matching *)
  \[Lambda]6D = ContriSE - ContriSS - ContriVV - ContriFF - ContriSSS // SparseArray;
];


(* ::Subsection:: *)
(*SSVV transverse*)


(*
    Calculates Scalar-Vector gauge couplings in the soft theory.
*)
TransverseSSVV[] := Module[
  {
    CouplingSV, ContriSVTemp, ContriSV, CouplingVVVV, ContriVVVV, CouplingFFFF,
    HabIJFAE, YTemp2, ContriFFFFTemp, ContriFFFFTemp2, ContriFFFF, HgY, HgY4, AE1, AE2,
    CouplingFFFF2, ContriFFFF2, ContriSEScalar, ContriSEVector
  },

  (* Verbose output for debugging *)
  VerbosePrint["Calculating Transverse-Vector Couplings"];

  (* Scalar-Vector coupling *)
  CouplingSV = 1 / (16 \[Pi]^2) * 3 / 4 * Lb;
  ContriSVTemp = Transpose[Simplify[Transpose[Transpose[Flatten[HabijV, {2, 3}], {2, 1, 3}], {1, 3, 2}] . Flatten[HabijV, {2, 4}]], {1, 3, 2, 4}];
  ContriSV = CouplingSV * Simplify[(ContriSVTemp + Transpose[ContriSVTemp, {2, 1, 3, 4}])] // SimplifySparse;

  (* Vector-Vector coupling *)
  CouplingVVVV = -((3 * Lb) / (64 \[Pi]^2));
  ContriVVVV = CouplingVVVV * Transpose[Transpose[Flatten[GabcdV, {2, 3}], {2, 1, 3}], {1, 3, 2}] . Flatten[HabijV, {1, 2}] // SimplifySparse;

  (* Fermion-Fermion coupling *)
  CouplingFFFF = -1 * 1 / (16 \[Pi]^2) * Lf;
  HabIJFAE = gvff . Transpose[gvff, {2, 1, 3}];
  YTemp2 = Ysff . Transpose[YsffC, {2, 1, 3}];
  
  ContriFFFFTemp = Transpose[Transpose[Flatten[HabIJFAE, {2, 4}], {2, 1, 3}], {1, 3, 2}] . Flatten[Transpose[YTemp2, {1, 4, 3, 2}], {2, 4}];
  ContriFFFFTemp2 = ContriFFFFTemp + Transpose[Transpose[ContriFFFFTemp, {1, 2, 4, 3}], {2, 1, 3, 4}];
  ContriFFFF = CouplingFFFF * Simplify[ContriFFFFTemp2 + Transpose[ContriFFFFTemp2, {2, 1, 3, 4}]];

  (* Additional Fermion-Fermion contribution *)
  HgY = gvff . Transpose[Ysff, {2, 1, 3}];
  HgY4 = Transpose[gvff, {1, 3, 2}] . Transpose[YsffC, {2, 1, 3}];
  
  AE1 = Transpose[Transpose[Flatten[HgY, {2, 4}], {2, 1, 3}], {1, 3, 2}] . Flatten[Transpose[HgY4, {1, 4, 3, 2}], {2, 4}] // Transpose[#, {1, 3, 2, 4}] &;
  AE2 = Transpose[Transpose[Flatten[HgY4, {2, 4}], {2, 1, 3}], {1, 3, 2}] . Flatten[Transpose[HgY, {1, 4, 3, 2}], {2, 4}] // Transpose[#, {1, 3, 2, 4}] &;
  
  CouplingFFFF2 = 2 * (-1) * 2 * 1 / 2 * 1 / (16 \[Pi]^2) * Lf;
  ContriFFFF2 = CouplingFFFF2 * Simplify[AE1 + AE2];

  (* Self-energy contributions *)
  ContriSEScalar = -HabijV . ZijS - Transpose[Transpose[HabijV, {1, 2, 4, 3}] . ZijS, {1, 2, 4, 3}] // SparseArray;
  ContriSEVector = -ZabT . HabijV - Transpose[ZabT . Transpose[HabijV, {2, 1, 3, 4}], {2, 1, 3, 4}] // SparseArray // SimplifySparse;

  (* Final result: total transverse-Scalar-Vector couplings *)
  GvvssT = ContriSV + ContriVVVV + ContriFFFF + ContriFFFF2 + ContriSEScalar + ContriSEVector // SparseArray;
];


(* ::Subsection:: *)
(*SSVV longitudinal*)


(*
    Calculates Scalar and Temporal-Scalar Quartics in the Soft Theory.
*)
LongitudinalSSVV[] := Module[
  {
    CouplingSV, ContriSVTemp, ContriSV, CouplingSSS, ContriSSS, CouplingVVVV, ContriVVVV,
    CouplingFFFF, HabIJFAE, YTemp2, ContriFFFFTemp, ContriFFFFTemp2, ContriFFFF,
    HgY, HgY4, AE1, AE2, CouplingFFFF2, ContriFFFF2, ContriSEScalar, ContriSEVector
  },

  (* Verbose output for debugging *)
  VerbosePrint["Calculating Scalar-Temporal-Vector Couplings"];

  (* Scalar-Vector coupling contribution *)
  CouplingSV = 1 / (16 \[Pi]^2) * (3 / 4 * Lb - 1 / 2);
  ContriSVTemp = Transpose[Simplify[Transpose[Transpose[Flatten[HabijV, {2, 3}], {2, 1, 3}], {1, 3, 2}] . Flatten[HabijV, {2, 4}]], {1, 3, 2, 4}];
  ContriSV = CouplingSV * Simplify[(ContriSVTemp + Transpose[ContriSVTemp, {2, 1, 3, 4}])] // SimplifySparse;

  (* Scalar-Scalar coupling contribution *)
  CouplingSSS = 1 / (8 \[Pi]^2);
  ContriSSS = CouplingSSS * Transpose[Transpose[Flatten[Habij, {3, 4}], {2, 1, 3}], {1, 3, 2}] . Flatten[\[Lambda]4, {1, 2}] // SimplifySparse;

  (* Vector-Vector coupling contribution *)
  CouplingVVVV = -((3 * Lb + 22) / (64 \[Pi]^2));
  ContriVVVV = CouplingVVVV * Transpose[Transpose[Flatten[GabcdV, {2, 3}], {2, 1, 3}], {1, 3, 2}] . Flatten[HabijV, {1, 2}] // SimplifySparse;

  (* Fermion-Fermion coupling contribution *)
  CouplingFFFF = -1 * 1 / (16 \[Pi]^2) * (Lf - 2);
  HabIJFAE = gvff . Transpose[gvff, {2, 1, 3}];
  YTemp2 = Ysff . Transpose[YsffC, {2, 1, 3}];
  ContriFFFFTemp = Transpose[Transpose[Flatten[HabIJFAE, {2, 4}], {2, 1, 3}], {1, 3, 2}] . Flatten[Transpose[YTemp2, {1, 4, 3, 2}], {2, 4}];
  ContriFFFFTemp2 = ContriFFFFTemp + Transpose[Transpose[ContriFFFFTemp, {1, 2, 4, 3}], {2, 1, 3, 4}] // SimplifySparse;
  ContriFFFF = CouplingFFFF * Simplify[ContriFFFFTemp2 + Transpose[ContriFFFFTemp2, {2, 1, 3, 4}]];

  (* Additional Fermion-Fermion contribution *)
  HgY = gvff . Transpose[Ysff, {2, 1, 3}];
  HgY4 = Transpose[gvff, {1, 3, 2}] . Transpose[YsffC, {2, 1, 3}];
  
  AE1 = Transpose[Transpose[Flatten[HgY, {2, 4}], {2, 1, 3}], {1, 3, 2}] . Flatten[Transpose[HgY4, {1, 4, 3, 2}], {2, 4}] // Transpose[#, {1, 3, 2, 4}] &;
  AE2 = Transpose[Transpose[Flatten[HgY4, {2, 4}], {2, 1, 3}], {1, 3, 2}] . Flatten[Transpose[HgY, {1, 4, 3, 2}], {2, 4}] // Transpose[#, {1, 3, 2, 4}] &;
  
  CouplingFFFF2 = 2 * (-1) * 2 * 1 / 2 * 1 / (16 \[Pi]^2) * Lf;
  ContriFFFF2 = CouplingFFFF2 * Simplify[AE1 + AE2];

  (* Self-energy contributions *)
  ContriSEScalar = -HabijV . ZijS - Transpose[Transpose[HabijV, {1, 2, 4, 3}] . ZijS, {1, 2, 4, 3}] // SimplifySparse // SparseArray;
  ContriSEVector = -ZabL . HabijV - Transpose[ZabL . Transpose[HabijV, {2, 1, 3, 4}], {2, 1, 3, 4}] // SparseArray // SimplifySparse;

  (* Final result: Total scalar-temporal-vector couplings *)
  GvvssL = ContriSV + ContriSSS + ContriVVVV + ContriFFFF + ContriFFFF2 + ContriSEScalar + ContriSEVector // SimplifySparse // SparseArray;
];


(* ::Section:: *)
(*Tadpoles*)


(* ::Subsection:: *)
(*1loop scalar tadpole*)


(*
    Calculates the 1-loop tadpole in the soft theory.
*)
TadPole[] := Module[
  {
    ContriS, fac, Temp1, Temp1C, ContriFF
  },

  VerbosePrint["Calculating 1-Loop Tadpoles"];

  (* Scalar bubble contribution (from one-loop with cubic coupling).
     Minus sign comes from the Feynman rule. *)
  ContriS = -1/2 * T^2 / 12 * TensorContract[\[Lambda]3, {{2, 3}}];

  (* Fermion mass-insertion contribution.
     Includes both Y and Y\[Dagger] contractions with fermion masses. *)
  fac = -(T^2 / 12);
  Temp1 = Contract[Ysff, \[Mu]IJFC, {{2, 4}, {3, 5}}] // SimplifySparse;
  Temp1C = Contract[YsffC, \[Mu]IJF, {{2, 4}, {3, 5}}] // SimplifySparse;

  (* Symmetrized fermionic contribution; factor (-1)^2 is from Wick contraction and symmetry. *)
  ContriFF = (1/2) * fac * (Temp1 + Temp1C);

  (* Final result; minus sign from matching to effective potential. *)
  TadPoleLO = -ContriS - ContriFF;
];


(* ::Subsection:: *)
(*2loop scalar tadpole*)


(*
    Calculates the 2-loop tadpole in the soft theory.
*)
TadPole2Loop[] := Module[
	{
	    (* General variables *)
	    \[Kappa], I1, I2M2, I1I2, dI1I2, I211M020, I2M2I2, I3M2I1,
	    I1p, I3M2p, I2p, I4M4p, I2M2p, IF1p, IF2p, IF3M2p, IF2M2p,
	
	    (* Contractions & results *)
	    \[Lambda]Help, Yhelp1, Yhelp2, YHelp2, YHelpC2, FHelp,
	    temp1, temp1C, temp2, temp3, temp4, temp5, temp6,
	    Contri1, Contri2, Contri3, Contri4, Contri5, Contri6,
	    Contri7, Contri71, Contri72,
	    ContriCTS, ContriCTFF, ContriF, TotalResult
	},

	VerbosePrint["Calculating 2-Loop Tadpoles"];

	(* --- Master integrals --- *)
	\[Kappa]=1/(16 \[Pi]^2);
	I1=T^2/(12 \[Epsilon])+T^2 Lbbb;
	I2M2=T^2/(24 )(-1/\[Epsilon]-  12Lbbb+2);
	I1I2=\[Kappa] T^2/(12 \[Epsilon]b); 
	dI1I2=3\[Kappa] T^2/(12 \[Epsilon]b)-\[Kappa] T^2/6;
	I211M020=T^2/(192 \[Pi]^2);
	I2M2I2=T^2/(192 \[Pi]^2)-T^2/(384 \[Pi]^2 \[Epsilon]bbM);
	I3M2I1=T^2/(768 \[Pi]^2 \[Epsilon]bbM)+T^2/(384 \[Pi]^2);
	
	Clear[LbbM];
	I1p=T^2/(12 )+\[Epsilon] T^2 Lbbb;
	I3M2p=1/(64 \[Pi]^2 \[Epsilon]bp)+1/(32 \[Pi]^2);
	I2p=1/(16 \[Pi]^2 \[Epsilon]bp);
	I4M4p=1/(128 \[Pi]^2 \[Epsilon]bp)+1/(48 \[Pi]^2);
	I2M2p=-T^2/24+T^2 (((-1/24)LbbM - (-1)*1/24 Lb )+1/12)\[Epsilon];
	IF1p=-T^2/(24)+1/24 T^2 Lfff \[Epsilon];
	IF2p=1/(16 (\[Pi]^2) ) (1/\[Epsilon]+Lf);
	IF3M2p=1/(64 \[Pi]^2)(1/\[Epsilon]+Lf+2);
	IF2M2p=(-T^2/24-T^2 1/48* Lfff)\[Epsilon]+T^2/48;

	\[Lambda]Help=TensorContract[\[Lambda]4,{3,4}];

	(* --- Pure diagram contributions --- *)
	Contri1 = (1/4) * I2p * I1p * Contract[\[Lambda]3, \[Lambda]Help, {{2, 4}, {3, 5}}] // SimplifySparse;
	Contri2 = -1/4 * (D - 1) * I2p * I1p * Contract[\[Lambda]3, \[CapitalLambda]g, {{2, 4}, {3, 5}}] // SimplifySparse;
	Contri3 = (-2)/2 * I2p * IF1p * Contract[\[Lambda]3, Ysij, {{2, 4}, {3, 5}}] // SimplifySparse;

	(* --- Scalar Mass insertion --- *)
	Contri4 = 1/2 * I2p * Contract[\[Lambda]3, \[Mu]ij, {{2, 4}, {3, 5}}] // SimplifySparse;


	(* --- Fermion Mass insertion --- *)
	Yhelp1 = \[Mu]IJF . \[Mu]IJFC . \[Mu]IJF;
	Yhelp2 = \[Mu]IJFC . \[Mu]IJF . \[Mu]IJFC;
	Contri71 = Contract[YsffC, Yhelp1, {{2, 4}, {3, 5}}] // SimplifySparse;
	Contri72 = Contract[Ysff, Yhelp2, {{2, 4}, {3, 5}}] // SimplifySparse;
	Contri7 = -IF2p * (Contri71 + Contri72);

	(* --- Fermion loop contributions --- *)
	YHelp2 = Contract[YsffC, Ysff, {{1, 4}, {2, 5}}] // SimplifySparse;
	YHelpC2 = Contract[Ysff, YsffC, {{1, 4}, {2, 5}}] // SimplifySparse;
	temp1 = Contract[Ysff, \[Mu]IJFC, YHelp2, {{2, 4}, {5, 7}, {3, 6}}] // SimplifySparse;
	temp1C = Contract[YsffC, \[Mu]IJF, YHelpC2, {{2, 4}, {5, 7}, {3, 6}}] // SimplifySparse;
	Contri5 = -(-1)*(IF1p * IF2p - IF2p * I1p) * (temp1 + temp1C) // SimplifySparse;

	FHelp = TensorContract[HabIJF, {1, 2}] // SimplifySparse;
	temp1 = Contract[Ysff, \[Mu]IJFC, FHelp, {{2, 4}, {5, 6}, {3, 7}}] // SimplifySparse;
	temp1C = Contract[YsffC, \[Mu]IJF, FHelp, {{2, 4}, {5, 6}, {3, 7}}] // SimplifySparse;
	Contri6 = -(D - 2) * (IF1p * IF2p - IF2p * I1p) * (temp1 + temp1C) // SimplifySparse;

	(* --- Counterterm contribution --- *)
	
	ContriCTS=-1/2 I1p \[Epsilon]^-1*(
		+ TensorContract[Z\[Lambda]ijk,{{2,3}}]
		+ 1/2*TensorContract[\[Gamma]ij . \[Lambda]3,{2,3}]
		)//SparseArray//SimplifySparse;
		
	temp1=Contract[ZYsij,\[Mu]IJFC,{{2,4},{3,5}}]//SimplifySparse;
	temp2=Contract[ZYsijC,\[Mu]IJF,{{2,4},{3,5}}]//SimplifySparse;
	temp3=Contract[Ysff,Z\[Mu]IJFC,{{2,4},{3,5}}]//SimplifySparse;
	temp4=Contract[YsffC,Z\[Mu]IJF,{{2,4},{3,5}}]//SimplifySparse;
	temp5=1/2Contract[\[Gamma]ij,Ysff,\[Mu]IJFC,{{2,3},{4,6},{5,7}}]//SimplifySparse;
	temp6=1/2Contract[\[Gamma]ij,YsffC,\[Mu]IJF,{{2,3},{4,6},{5,7}}]//SimplifySparse;
	ContriCTFF=1/2*2*IF1p*\[Epsilon]^-1 (-1)^2*(temp1+temp2+temp3+temp4+temp5+temp6);
	
	(* --- Self-Energy correction --- *)
	ContriF=-TadPoleLO . ZijS//SimplifySparse;
	
	(* --- Final result --- *)
	TotalResult = ContriF-Contri1-Contri2-Contri3-ContriCTS-Contri4-ContriCTFF-Contri5-Contri6-Contri7//Normal; (*minus signs from matching*)
	TadPoleNLO=Series[
		TotalResult/.{
			D->4-2\[Epsilon],
			\[Epsilon]bp->(1/\[Epsilon]+Lb)^-1,
			\[Epsilon]b->(1/\[Epsilon]+Lbb)^-1,
			\[Epsilon]BF->(1/\[Epsilon]+LBF)^-1,
			\[Epsilon]F->(1/\[Epsilon]+LFF)^-1,
			\[Epsilon]FB->(1/\[Epsilon]+LFB)^-1,
			\[Epsilon]bbM->(1/\[Epsilon]+LbbM)^-1
		},{\[Epsilon],0,0}
	]/.ReplaceLb//Normal//Coefficient[#,\[Epsilon],0]&//Simplify//FullSimplify;
];


(* ::Section:: *)
(*Counter-terms and beta functions*)


(* ::Subsection:: *)
(*Running hard to soft*)


(*
	Calculates counter-terms, and beta functions, in the soft theory.
*)
RGRunningHardToSoft[]:=Module[{},
	VerbosePrint["RG-evolving from \[Mu] to \[Mu]3"];

(*
	All calculations are done within the 3d theory.
	So effective couplings are used
	The module finds all \[Epsilon] poles to determine the scalar-mass counterterm
*)

	I111=1/(4)1/(16 \[Pi]^2);(*All possible divergences comes from I111, which is the scalar sunset in 3-D*)

(*Help Tensors*)
	VarGauge=GaugeCouplingNames//DeleteDuplicates;
	SubGauge=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarGauge}]; 
(*VarGauge and SubGauge are used because the running happens in the soft theory*)


	gvssHeavy=gvvv//Normal//ReplaceAll[#,SubGauge]&//SparseArray; (*Heavy meaning temporal vectors*)
	gvssVTot=Table[ArrayFlatten[{{gvssHeavy[[a]],0},{0,gvss[[a]]//Normal//ReplaceAll[#,SubGauge]&}}],{a,1,Length[gvssHeavy]}]//SparseArray;
(* gvssVTot combines temporal vectors with scalars into a single vector-scalar^2 coupling.*)

	VarGauge=Join[\[Lambda]4//Normal//Variables]//DeleteDuplicates;
	SubGauge=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarGauge}];
(*VarGauge and SubGauge are used because the running happens in the soft theory*)

(*
	Because temporal vectors also talk with scalars, via quartics,
	they need to be treated on the same footing.
	There are ns+nv scalars in the soft theory
	First nv indices are allocated to the temporal scalars
*)

	\[Lambda]4p=\[Lambda]4//Normal//ReplaceAll[#,SubGauge]&//SparseArray; 
(*Changes made to \[Lambda]KVLight*)
	\[Lambda]KVLight=\[Lambda]KVec//SparseArray;(*Temporal vector^2 times scalar^2 couplings*)
	\[Lambda]KVHeavy=\[Lambda]4p//SparseArray;(*Original scalars*)
	\[Lambda]4Tot=SymmetrizedArray[{{1,1,1,1}->0},{nv+ns,nv+ns,nv+ns,nv+ns},Symmetric[{1,2,3,4}]]//SparseArray;
	\[Lambda]4Tot[[1;;nv,1;;nv,nv+1;;ns+nv,nv+1;;ns+nv]]=\[Lambda]KVLight//SparseArray; (*Populates the scalar quartic and ensures the symmetry*)
	\[Lambda]4Tot[[nv+1;;ns+nv,1;;nv,1;;nv,nv+1;;ns+nv]]=Transpose[Transpose[\[Lambda]KVLight[[1;;nv,1;;nv,1;;ns,1;;ns]],{1,3,2,4}],{2,1,3,4}]//SparseArray;
	\[Lambda]4Tot[[1;;nv,nv+1;;ns+nv,1;;nv,nv+1;;ns+nv]]=Transpose[\[Lambda]KVLight[[1;;nv,1;;nv,1;;ns,1;;ns]],{1,3,2,4}]//SparseArray;
	\[Lambda]4Tot[[1;;nv,nv+1;;ns+nv,nv+1;;ns+nv,1;;nv]]=Transpose[Transpose[\[Lambda]KVLight[[1;;nv,1;;nv,1;;ns,1;;ns]],{1,3,2,4}],{1,2,4,3}]//SparseArray;
	\[Lambda]4Tot[[nv+1;;ns+nv,1;;nv,nv+1;;ns+nv,1;;nv]]=Transpose[Transpose[Transpose[\[Lambda]KVLight[[1;;nv,1;;nv,1;;ns,1;;ns]],{1,3,2,4}],{1,2,4,3}],{2,1,3,4}]//SparseArray;
	\[Lambda]4Tot[[nv+1;;ns+nv,nv+1;;ns+nv,1;;nv,1;;nv]]=Transpose[Transpose[Transpose[Transpose[\[Lambda]KVLight[[1;;nv,1;;nv,1;;ns,1;;ns]],{1,3,2,4}],{1,2,4,3}],{2,1,3,4}],{1,3,2,4}]//SparseArray;
	\[Lambda]4Tot[[nv+1;;ns+nv,nv+1;;ns+nv,nv+1;;ns+nv,nv+1;;ns+nv]]=\[Lambda]KVHeavy//SparseArray;

(*\[Lambda]4ijkl: S^4,\[Lambda]KVec_abij; V^2(S^2)., ab=1,nv ijkl=1, ns*)
(*\[Lambda]4TotIJKL, I=1,...,nv,....,nv+ns*)

(*Help tensors*)
	\[CapitalLambda]\[Lambda]RG =Transpose[Flatten[\[Lambda]4Tot,{3,4}],{3,2,1}] . Flatten[\[Lambda]4Tot,{1,2}]//SparseArray;
	HabijRG=Transpose[Activate @ TensorContract[
        Inactive[TensorProduct][gvssVTot,gvssVTot], { {3, 5}}],{1,3,2,4}]//SparseArray;
	HabijVRG=HabijRG+Transpose[HabijRG,{2,1,3,4}]//SparseArray;
	GabcdVp=Simplify[ Activate @ TensorContract[
        Inactive[TensorProduct][gvssHeavy,gvssHeavy], {{3, 6}}]]//SparseArray;
	HijRG=TensorContract[TensorProduct[HabijRG],{{1,2}}]//SparseArray;


(**)
	Coeff=1/2*2*I111;
	Contri1RG=Coeff*Flatten[HijRG] . Flatten[\[Lambda]4Tot,{3,4}];

	Coeff=-1/4*(-1)*I111;
	help1=TensorContract[HabijRG,{3,4}]//SparseArray;
	Contri2RG=Coeff*Flatten[help1] . Flatten[HabijVRG,{1,2}];

	Coeff=1/3!*I111;
	Contri3RG=Coeff *TensorContract[\[CapitalLambda]\[Lambda]RG,{{2,3}}];


	HijV=Transpose[Transpose[Flatten[HabijVRG,{1,2}],{2,1,3}],{1,3,2}] . Flatten[HabijVRG,{1,2}];
	Coeff=1/2*(I111(1+1/2));
	Contri4RG=Coeff *TensorContract[HijV,{{2,4}}]//Simplify;

	GabcdV2=TensorContract[GabcdVp,{2,4}];
	Contri5RG=(-1/4*I111) /2*(-1)*Flatten[GabcdV2] . Flatten[HabijVRG,{1,2}];
        
    Contri6RG=(-1)1/4*(20/4 I111) *(-1)*Flatten[GabcdV2] . Flatten[HabijVRG,{1,2}];

	\[Delta]\[Mu]3d=Contri1RG+  Contri2RG+ Contri3RG+ Contri4RG+Contri5RG+ Contri6RG; (*All mass-counterterms in the soft theory*)       
        
	\[Beta]\[Mu]3ij=4 \[Delta]\[Mu]3d //Simplify; (*Finds the beta function in the soft theory*)
	Contri\[Beta]SoftToHard=\[Beta]\[Mu]3ij[[nv+1;;ns+nv,nv+1;;ns+nv]] Log[\[Mu]3/\[Mu]]//SimplifySparse;(*Running from the hard matching scale to the soft scale \[Mu]3*)
(*In the soft theory temporal scalars don't run. So only the original-scalar running is needed*)


(*Tadpoles*)

	Contri10RG=1/3!*I111*Activate@TensorContract[Inactive@TensorProduct[\[Lambda]4//SparseArray,\[Lambda]3//SparseArray],{{2,5},{3,6},{4,7}}];
	Contri11RG=1/2*I111*2Activate@TensorContract[Inactive@TensorProduct[\[Lambda]3//SparseArray,\[CapitalLambda]g//SparseArray],{{2,4},{3,5}}];

	VarTadpole=(Contri10RG+Contri11RG)//Normal//Variables;

(*Help Tensors*)
	VarGauge=GaugeCouplingNames//DeleteDuplicates;
	VarQuartic=\[Lambda]4//Normal//Variables//DeleteDuplicates;
	VarCubic=\[Lambda]3//Normal//Variables//DeleteDuplicates;
	VarTot=Join[VarGauge,VarQuartic,VarCubic];
	SubGauge=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarTot}]; 


	\[Delta]\[Mu]Tadpole3d=Contri10RG+Contri11RG//Normal//ReplaceAll[#,SubGauge]&//SparseArray;(*Same deal for tadpoles*)
	\[Beta]\[Mu]Tadpole=4  \[Delta]\[Mu]Tadpole3d //Simplify;
	ContriTadPoleSoftToHard=\[Beta]\[Mu]Tadpole Log[\[Mu]3/\[Mu]]//SimplifySparse;

];


(* ::Subsection:: *)
(*Counter terms*)


(*
	Calculates all counterterms required for two-loop calculations.
*)
CounterTerm[]:=Module[{},
(*Only performs the calculation once*)
If[CT==False,
	CT=True;
	CreateBasisVanDeVis[];
		
	VerbosePrint["Calculating CounterTerms"];
		\[Kappa]=1/(16 \[Pi]^2);

(*
	Scalar-field renormalization
*)
	SelfEnergySS=-3 \[CapitalLambda]g;
	ContriSS=SelfEnergySS;

	SelfEnergyFF=-1;
	ContriFF=1/2SelfEnergyFF (Ysij+YsijC);
	\[Gamma]ij=(ContriSS+ContriFF) \[Kappa]//SparseArray; (*Scalar anomalous dimension*)


(*
	Scalar mass renormalization
*)

(*Scalar-quartic contribution*)
	ContriSS=Flatten[\[Mu]ij] . Flatten[\[Lambda]4,{3,4}];

(*Fermion-mass contribution*)
	Fac=1/(8 \[Pi]^2)*2;
	Temp1=Contract[Ysff,\[Mu]IJFC,\[Mu]IJFC,Ysff,{{2,4},{3,6},{5,9},{7,10}}]//Normal;
	Temp1C=Contract[YsffC,\[Mu]IJF,\[Mu]IJF,YsffC,{{2,4},{3,6},{5,9},{7,10}}]//Normal;
	ContriFFMass1=-1/2Fac*(Temp1+Temp1C)//Simplify//PowerExpand;
	
	Temp1=Contract[Ysff,YsffC,\[Mu]IJFC,\[Mu]IJF,{{3,5},{2,7},{6,10},{8,9}}]//Normal;
	Temp1C=Contract[YsffC,Ysff,\[Mu]IJF,\[Mu]IJFC,{{3,5},{2,7},{6,10},{8,9}}]//Normal;
	ContriFFMass2=-Fac*(Temp1+Temp1C)//Simplify//PowerExpand;
  
(*Cubic-coupling contribution*)
	ContriAnom=-\[Mu]ij . \[Gamma]ij-Transpose[\[Mu]ij . \[Gamma]ij];
	ContriCubic=1/2*1/(16 \[Pi]^2)*2 Activate@TensorContract[Inactive@TensorProduct[\[Lambda]3,\[Lambda]3],{{2,5},{3,6}}];
	\[Beta]mij=(\[Kappa] ContriSS+ContriAnom+ContriFFMass1+ContriFFMass2+ContriCubic)//SparseArray;

(*
	Scalar quartic renormalization
*)

(*Scalar-quartic contribution*)
	ContriSS=\[Kappa](\[CapitalLambda]\[Lambda]+Transpose[\[CapitalLambda]\[Lambda],{1,3,2,4}]+Transpose[\[CapitalLambda]\[Lambda],{1,4,3,2}]);

(*Gauge contribution*)
	CouplingVV2=2*1/(16 \[Pi]^2)(3 1)/2 ;
	ContriVVTemp=CouplingVV2*Transpose[Transpose[Flatten[HabijV,{1,2}],{2,1,3}],{1,3,2}] . Flatten[HabijV,{1,2}];
	ContriVV=ContriVVTemp+Transpose[ContriVVTemp,{1,3,2,4}]+Transpose[ContriVVTemp,{1,4,2,3}]//Simplify;


(*Yukawa contribution*)
	CouplingFF=2*2*1/(16 \[Pi]^2)(2 )(-1)*1/4;
	ContriFF=CouplingFF*12Symmetrize[Yhelp,Symmetric[{1,2,3,4}]]//SparseArray;

(*Anomalous dimension contribution*)
	ContriSETemp=\[Gamma]ij . \[Lambda]4;
	ContriAnom=ContriSETemp+Transpose[ContriSETemp,{2,1,3,4}]+Transpose[ContriSETemp,{3,1,2,4}]+Transpose[ContriSETemp,{4,1,2,3}]//Simplify;

	\[Beta]\[Lambda]ijkl= ContriSS+ ContriVV+ ContriFF-  ContriAnom//Simplify//SparseArray;
	Z\[Lambda]ijkl=\[Beta]\[Lambda]ijkl/(2 \[Epsilon])//SparseArray;(*renormalization constants*)


(*
	Scalar cubic renormalization
*)
	SelfEnergySSC=1/(16 \[Pi]^2) ;
	ContriSSCTemp=Simplify[Activate @ TensorContract[Inactive[TensorProduct][\[Lambda]4,\[Lambda]3], {{3, 5},{4,6}}]];
	ContriSSC=SelfEnergySSC(ContriSSCTemp+Transpose[ContriSSCTemp,{1,3,2}]+Transpose[ContriSSCTemp,{3,2,1}]);


(*Fermion mass contributions*)
	Contri1=Table[Tr[a . \[Mu]IJFC . b . c],{a,Ysff},{b,Ysff},{c,YsffC}];
	Contri2=Table[Tr[a . \[Mu]IJF . b . c],{a,YsffC},{b,YsffC},{c,Ysff}];
	ContriFF=-1/(8 \[Pi]^2) 3 Symmetrize[(Contri1+Contri2),Symmetric[{1,2,3}]]//SparseArray;
	
(*Anomalous dimension contribution*)
	ContriAnom=Simplify[Activate @ TensorContract[Inactive[TensorProduct][\[Gamma]ij,\[Lambda]3], {{2, 3}}]];
	\[Beta]\[Lambda]ijk=(-ContriAnom-Transpose[ContriAnom,{2,1,3}]-Transpose[ContriAnom,{3,2,1}]+ContriSSC+ContriFF)//SparseArray//SimplifySparse; (*beta function*)
	Z\[Lambda]ijk=\[Beta]\[Lambda]ijk/2;(*renormalization constants*)

(*
	Vector-field renormalization
*)
	SelfEnergySS=1/6 Hg;
	fac=13/6;
	ContriVVV=fac*TensorContract[GabcdV,{2,4}];
	SelfEnergyFF=(-1)(2/3 );
(*General nF modification*)
	HabIJFnF=HabIJF . NFMat;
(************************)
	ContriFF=SelfEnergyFF*TensorContract[HabIJFnF,{3,4}];

	\[Gamma]ab=( SelfEnergySS+1ContriVVV+  ContriFF) \[Kappa]//SparseArray; (*vector anomalous dimension*)


(*
	Gauge-coupling renormalization
*)
	ContriSS=-1/(16 \[Pi]^2)1*1/2*Simplify[Transpose[Transpose[Flatten[HabijV,{3,4}],{2,1,3}],{1,3,2}] . Flatten[\[Lambda]4,{1,2}]];
	CouplingSV= 1/(16 \[Pi]^2)*3/4*1 ;
	ContriSVTemp=Transpose[Transpose[Flatten[HabijV,{2,4}],{2,1,3}],{1,3,2}] . Flatten[HabijV,{2,3}]//Transpose[#,{1,3,2,4}]&;
	ContriSV=CouplingSV Simplify[(ContriSVTemp+Transpose[ContriSVTemp,{2,1,3,4}])];

	CouplingSSS=1/(16 \[Pi]^2)1/2 1;
	ContriSSS=2*CouplingSSS*Simplify[Transpose[Transpose[Flatten[Habij,{3,4}],{2,1,3}],{1,3,2}] . Flatten[\[Lambda]4,{1,2}]];

	CouplingVVV= 1/(16 \[Pi]^2)*(-1)/2(9/2 1) ;
	ContriVVV=CouplingVVV*Simplify[Transpose[Transpose[Flatten[GabcdV,{2,4}],{2,1,3}],{1,3,2}] . Flatten[HabijV,{1,2}]];
 
	CouplingVVVV= 1/(16 \[Pi]^2)*(-1)(3 1) ;
	ContriVVVV=CouplingVVVV*Simplify[Transpose[Transpose[Flatten[GabcdV,{2,3}],{2,1,3}],{1,3,2}] . Flatten[HabijV,{1,2}]];
 
 (*Fermion contributions*)
	CouplingFFFF=(-1) 1/(16 \[Pi]^2)( 1) ;
	HabIJFAE=Simplify[ Activate @ TensorContract[
        Inactive[TensorProduct][gvff,gvff], {{3, 5}}]];
	YTemp2=Simplify[ Activate @ TensorContract[
        Inactive[TensorProduct][Ysff,YsffC], {{3, 5}}]];
	ContriFFFFTemp=Transpose[Transpose[Flatten[HabIJFAE,{2,4}],{2,1,3}],{1,3,2}] . Flatten[Transpose[YTemp2,{1,4,3,2}],{2,4}];
	ContriFFFFTemp2=ContriFFFFTemp+Transpose[Transpose[ContriFFFFTemp,{1,2,4,3}],{2,1,3,4}];
	ContriFFFF=CouplingFFFF*Simplify[ContriFFFFTemp2+Transpose[ContriFFFFTemp2,{2,1,3,4}]];


	HgY=Simplify[ Activate @ TensorContract[
        Inactive[TensorProduct][gvff,Ysff], {{3, 5}}]]//SparseArray;
	HgY4=Simplify[ Activate @ TensorContract[
        Inactive[TensorProduct][gvff,YsffC], {{2, 5}}]]//SparseArray;
	AE1=Transpose[Transpose[Flatten[Transpose[HgY,{1,4,3,2}],{2,4}],{2,1,3}],{1,3,2}] . Flatten[HgY4,{2,4}]//Transpose[#,{1,3,2,4}]&;
	AE2=Transpose[Transpose[Flatten[Transpose[HgY4,{1,4,3,2}],{2,4}],{2,1,3}],{1,3,2}] . Flatten[HgY,{2,4}]//Transpose[#,{1,3,2,4}]&;
	CouplingFFFF2=2* (-1)*2*1/2*1/(16 \[Pi]^2)( 1) ;
	ContriFFFF2=CouplingFFFF2*Simplify[ AE1 +  AE2 ];

(*Anomalous-dimension contributions*)
	ContriSEScalar=-HabijV . \[Gamma]ij-Transpose[Transpose[HabijV,{1,2,4,3}] . \[Gamma]ij,{1,2,4,3}]//SparseArray//Simplify;
	ContriSEVector=-\[Gamma]ab . HabijV-Transpose[\[Gamma]ab . Transpose[HabijV,{2,1,3,4}],{2,1,3,4}]//SparseArray//Simplify;

	\[Beta]vvss=-2(ContriSS+   ContriSV+ContriSSS+  ContriVVV+ ContriVVVV+ ContriFFFF+ContriFFFF2+   -1/2( ContriSEScalar   + ContriSEVector))//Simplify//SparseArray; (*beta function*)
	Zgvvss=\[Beta]vvss/(2 \[Epsilon])//SparseArray;(*renormalization constant*)

(*
	Ghost-field renormalization
*)
	fac=3/4;
	ContriVVV=fac TensorContract[GabcdV,{{2,4}}];
	\[Gamma]ab\[Eta]=(1ContriVVV) \[Kappa]//SparseArray; (*Who you gonna call?*)


(*
	Non-abelian coupling renormalization
*)
(*There are only anomalous-dimension contributions*)
	ContriAnomVV=gvvv . \[Gamma]ab+\[Gamma]ab\[Eta] . gvvv+Transpose[\[Gamma]ab\[Eta] . Transpose[gvvv,{2,1,3}],{2,1,3}]//SparseArray//Simplify;

	\[Beta]gvvv=-ContriAnomVV//SparseArray; (*beta functions*)
	Zgvvv=\[Beta]gvvv/2//SparseArray; (*Renormalization constants*)
(*For some renormalization constants I decided in my infinite wisdom to not include \[Epsilon] poles. Should probably make the same convention everywhere.*)
(*I blame this inconsistency on flouridation of water.*)

(*
	Fermion-field renormalization
*)
	ContriSS=(-1)/2*TensorContract[YTemp,{1,2}];
	\[Gamma]IJF=\[Kappa] (ContriSS)//SparseArray; (*Anomalous dimension*)



(*
	Yukawa-coupling renormalization
*)
	ContriSS=2/(16 \[Pi]^2)*Transpose[Transpose[Flatten[YTemp,{1,4}],{2,1,3}],{1,3,2}] . Flatten[Ysff,{1,2}];
	HelpAE=Transpose[Transpose[gvff,{2,1,3}],{1,3,2}] . gvff; 
	ContriVV=6/(16 \[Pi]^2)*Transpose[Flatten[Ysff,{2,3}]] . Flatten[HelpAE,{2,4}];

(*Anomalous dimension contribution*)
	ContriAnomSS=-\[Gamma]ij . Ysff;
	ContriAnomFF=-(Ysff . \[Gamma]IJF+Transpose[Transpose[Ysff,{1,3,2}] . \[Gamma]IJF,{1,3,2}]);

	\[Beta]Ysij=   ContriSS+  ContriVV+  ContriAnomSS+ ContriAnomFF//Simplify//SparseArray; (*beta function*)
	ZYsij=\[Beta]Ysij/2//SparseArray;(*renormalization constants*)

(*
	Yukawa-coupling renormalization for conjugated tensor.
	A bit redundant, but it's not even close to being a bottleneck.
*)
	ContriSS=2/(16 \[Pi]^2)Simplify[ Activate @ TensorContract[
        Inactive[TensorProduct][YTempC,YsffC], {{1, 5},{4,6}}]];
	ContriVV=6/(16 \[Pi]^2)Simplify[ Activate @ TensorContract[
        Inactive[TensorProduct][YsffC,gvff,gvff], {{4, 7},{2,6},{3,9}}]];
	ContriAnomSS=-\[Gamma]ij . YsffC;
	ContriAnomFF=-(Transpose[YsffC,{1,3,2}] . Transpose[\[Gamma]IJF]+Transpose[YsffC . Transpose[\[Gamma]IJF],{1,3,2}]);        

	\[Beta]YsijC=  ContriSS+ContriVV+ ContriAnomSS+ContriAnomFF//Simplify//SparseArray;
	ZYsijC=\[Beta]YsijC/2//SparseArray;


	Zij=2*\[Gamma]ij/\[Epsilon] /2//SparseArray  ;
	Zmij=\[Beta]mij/2//SparseArray;

(*
	Gauge-coupling renormalization. Just in case there are U1's not interacting with the scalars.
*)
	Contri2=1/(16 \[Pi]^2) Transpose[Flatten[gvvv,{2,3}]] . Flatten[HabIJF,{1,2}];
(*Anomalous-dimension contribution*)
	ContriSE= \[Gamma]ab . gvff;

	\[Beta]gvff=(I 3 Contri2- ContriSE)//Simplify//SparseArray; (*beta function*)
	Zgvff=\[Beta]gvff/2//SparseArray;(*Renormalization-constant*)


(*
	Fermion-mass renormalization.
*)
	YtempMass= Simplify[ Activate @ TensorContract[
        Inactive[TensorProduct][Ysff,\[Mu]IJFC], {{3, 5}}]];
	ContriSS=2/(16 \[Pi]^2)Simplify[ Activate @ TensorContract[
        Inactive[TensorProduct][YtempMass,Ysff], {{1,4},{3,5}}]];
	ContriVV=6/(16 \[Pi]^2)Simplify[ Activate @ TensorContract[
        Inactive[TensorProduct][\[Mu]IJF,gvff,gvff], {{3,6},{1,5},{2,8}}]];
        
        (*Anomalous-dimension contribution*)
	ContriAnomFF=-(\[Gamma]IJF . \[Mu]IJF+Transpose[\[Gamma]IJF . \[Mu]IJF]);
	
	\[Beta]\[Mu]IJF=   ContriSS+  ContriVV+  ContriAnomFF//Simplify//SparseArray; (*beta function*)
	Z\[Mu]IJF=\[Beta]\[Mu]IJF/2//SparseArray;(*renormalization constant*)


(*
	Conjugated Fermion-mass renormalization.
*)
	YtempMass= Simplify[ Activate @ TensorContract[
        Inactive[TensorProduct][YsffC,\[Mu]IJF], {{3, 5}}]];
	ContriSS=(-1)2/(16 \[Pi]^2)Simplify[ Activate @ TensorContract[
        Inactive[TensorProduct][YtempMass,YsffC], {{1,4},{3,5}}]];
	ContriVV=6/(16 \[Pi]^2)Simplify[ Activate @ TensorContract[
        Inactive[TensorProduct][\[Mu]IJFC,gvff,gvff], {{3,6},{1,5},{2,8}}]];
	ContriAnomFF=-(\[Gamma]IJF . \[Mu]IJFC+Transpose[\[Gamma]IJF . \[Mu]IJFC]);

	\[Beta]\[Mu]IJFC=   ContriSS+  ContriVV+  ContriAnomFF//Simplify//SparseArray;
	Z\[Mu]IJFC=\[Beta]\[Mu]IJFC/2//SparseArray;


(*
	Tadpole renormalization
*)
	Contri4=1/2*I2p*Activate@TensorContract[Inactive@TensorProduct[\[Lambda]3,\[Mu]ij],{{2,4},{3,5}}];


(*Fermion Mass insertion*)
	Yhelp1=\[Mu]IJF . \[Mu]IJFC . \[Mu]IJF;
	Yhelp2=\[Mu]IJFC . \[Mu]IJF . \[Mu]IJFC;
	Contri71=Activate@TensorContract[Inactive@TensorProduct[YsffC,Yhelp1],{{2,4},{3,5}}];
	Contri72=Activate@TensorContract[Inactive@TensorProduct[Ysff,Yhelp2],{{2,4},{3,5}}];
	ContriF=-1/(16 \[Pi]^2)*2(Contri71+Contri72);

	ContriS=1/(16 \[Pi]^2)*Activate@TensorContract[Inactive@TensorProduct[\[Lambda]3,\[Mu]ij],{{2,4},{3,5}}];
	
	ContriAnom=-\[Lambda]1 . \[Gamma]ij;

	\[Beta]\[Lambda]1=ContriF+ContriS+ContriAnom;
	Z\[Lambda]1=\[Beta]\[Lambda]1/2;

	If[mode>=3,
(*Effective higher-order couplings*)

(*Scalar sextic operator*)
(*Scalar loops*)

(*Scalar loop with mixed \[Lambda]6 and \[Lambda]4 vertices*)
		\[CapitalLambda]\[Lambda]6tot=Flatten[\[Lambda]4,{{1},{2},{3,4}}] . Flatten[\[Lambda]6,{1,2}]//SimplifySparse;
		Prefac=15/2*1/(16 \[Pi]^2);
		SymHelp=Symmetrize[\[CapitalLambda]\[Lambda]6tot,Symmetric]//SparseArray;
		ContriSSS=Prefac*SymHelp;

(*Field-strength renormalization*)

		ContriSETemp=-1/2*\[Gamma]ij . \[Lambda]6;
		ContriSE=6*Symmetrize[ContriSETemp,Symmetric]//SparseArray//SimplifySparse;

(*Minus sign from matching*)
		\[Beta]\[Lambda]6ijklnm=-2*(ContriSE+ContriSSS)//SparseArray; 
	];



	]

];


(* ::Section:: *)
(*Identify tensors*)


(* Helper functions *)
CleanedList[list_] := If[Length[list] > 1 && list[[1]] === 0, Rest[list], list];

MakeReplacementRules[list_, vars_] := Module[{modVars},
  modVars = RelationsBVariables[list, vars];
  Table[Delete[list, 1][[a]] -> modVars[[a]], {a, Length[Delete[list, 1]]}] // Flatten // Simplify
];

ProcessCoupling[expr_, sym_] := Module[{flat, cleaned, vars, mods, rules},
  flat = DeleteDuplicates@Flatten@Simplify[expr] // Sort;
  cleaned = CleanedList[flat];
  vars = Table[sym[a], {a, Length[cleaned]}];
  mods = RelationsBVariables[flat, vars];
  rules = MakeReplacementRules[flat, vars];
  {cleaned, vars, mods, rules}
];


(*
	Identifies couplings and masses in the soft theory.
*)

IdentifyTensorsDRalgo[]:=Module[
	{
		couplingTensor,
		HelpList, HelpListCleaned, HelpVar, HelpVarMod
	},

(*Option to keep 4d-normalization in the matching relations*)
	Tfac=T;
	If[normalization4D, Tfac=1];

	If[mode>=1,
		(* --- Scalar quartic couplings --- *)
		VerbosePrint["Calculating Quartic Tensor"];
		
		couplingTensor = Tfac(\[Lambda]4+\[Lambda]3D);
		HelpList=DeleteDuplicates@Flatten[couplingTensor]//Sort//Simplify;
		HelpVarMod=RelationsBVariables3[HelpList]//ReplaceAll[#,\[Lambda]VL[v1_]->\[Lambda][v1]]&;
		HelpListCleaned=CleanedList[HelpList];
		HelpSolveQuartic=Table[{HelpListCleaned[[a]]->HelpVarMod[[a]]},{a,1,HelpListCleaned//Length}]//Flatten//Simplify;
		\[Lambda]3DS=couplingTensor//SimplifySparse//Normal//ReplaceAll[#,HelpSolveQuartic]&//SparseArray;


		If[Length[SparseArray[\[Lambda]3DS]["NonzeroValues"]]!=Length[SparseArray[\[Lambda]4]["NonzeroValues"]],
			Print["Detected 1-loop Scalar Quartics not defined at tree-level"];
			Print["Please Check if you defined all couplings allowed by symmetry"];
		];

		(* --- Scalar cubic couplings --- *)
		VerbosePrint["Calculating Cubic Tensor "];
		
		couplingTensor = Tfac^(1/2)*(\[Lambda]3CS+\[Lambda]3);
		{HelpListCleaned, HelpVar, HelpVarMod, HelpSolveCubicS} = ProcessCoupling[couplingTensor, cSS];	
		\[Lambda]3CSRed=couplingTensor//Normal//Simplify//FullSimplify//ReplaceAll[#,HelpSolveCubicS]&//SparseArray;

		If[Length[SparseArray[\[Lambda]3CSRed]["NonzeroValues"]]!=Length[SparseArray[\[Lambda]3]["NonzeroValues"]],
			Print["Detected 1-loop Scalar Cubic not defined at tree-level"];
			Print["Please Check if you defined all couplings allowed by symmetry"];
		];
	
		(* --- Vector-scalar couplings --- *)
		VerbosePrint["Calculating Vector-Scalar Interaction Tensor "];
		
		couplingTensor = Tfac*(HabijV+GvvssT);
		HelpList=DeleteDuplicates@Flatten@SimplifySparse[couplingTensor]//Sort;
		HelpListCleaned=CleanedList[HelpList];
		If[Length[HelpListCleaned]<1,(*Fix for when the tensor is empty*)
			\[Lambda]KVecT=couplingTensor//SparseArray;
			HelpSolveVecT={};
		,
			HelpVarMod=RelationsBVariables3[HelpList]//ReplaceAll[#,\[Lambda]VL[v1_]->\[Lambda]VT[v1]]&;
			HelpSolveVecT=Table[{HelpListCleaned[[a]]->HelpVarMod[[a]]},{a,1,HelpListCleaned//Length}]//Flatten//Simplify;
			\[Lambda]KVecT=couplingTensor//SimplifySparse//Normal//ReplaceAll[#,HelpSolveVecT]&//SparseArray;
		];

		(* --- Temporal-scalar/scalar cross couplings --- *)
		
		couplingTensor = -Tfac(HabijV+GvvssL);
		HelpList=DeleteDuplicates@Simplify@Flatten[couplingTensor]//Sort;
		HelpListCleaned=CleanedList[HelpList];
		If[Length[HelpListCleaned]<1,
			\[Lambda]KVec=couplingTensor//Normal//Simplify//SparseArray;
			HelpSolveVecL={};
		,
			HelpVarMod=RelationsBVariables3[HelpList];
			HelpSolveVecL=Table[{HelpListCleaned[[a]]->HelpVarMod[[a]]},{a,1,HelpListCleaned//Length}]//Flatten//Simplify;
			\[Lambda]KVec=couplingTensor//Normal//Simplify//ReplaceAll[#,HelpSolveVecL]&//SparseArray;
		];

		(* --- Temporal-Scalar quartics --- *)
		VerbosePrint["Calculating Temporal-Vector Quartics "];
		
		couplingTensor=Tfac*\[Lambda]AA;
		HelpList=DeleteDuplicates@SparseArray[Flatten@Simplify[couplingTensor]]//Sort//FullSimplify;
		HelpListCleaned=CleanedList[HelpList];
		HelpVar=Table[\[Lambda]VLL[a],{a,1,Rest[HelpList]//Length}];
		If[Length[HelpVar]<1,
			HelpVar=\[Lambda]VLL[1];
			HelpVarMod=RelationsBVariables[HelpList,HelpVar];
			HelpSolveQuarticL={HelpList[[1]]->HelpVarMod}//Flatten;

			\[Lambda]AAS=couplingTensor//Normal//Simplify//ReplaceAll[#,HelpSolveQuarticL]&//SparseArray;
		,
			HelpVarMod=RelationsBVariables3[HelpList]//ReplaceAll[#,\[Lambda]VL[v1_]->\[Lambda]VLL[v1]]&;
			HelpSolveQuarticL=Table[{HelpListCleaned[[a]]->HelpVarMod[[a]]},{a,1,HelpListCleaned//Length}]//Flatten//Simplify;
			\[Lambda]AAS=couplingTensor//Normal//FullSimplify//ReplaceAll[#,HelpSolveQuarticL]&//SparseArray;
		];

		(* --- Temporal-Scalar/scalar cross cubic couplings --- *)
		couplingTensor = Sqrt[Tfac]*GvvsL;
		{HelpListCleaned, HelpVar, HelpVarMod, HelpSolveCubicL} = ProcessCoupling[couplingTensor, \[Lambda]VVSL];	
		\[Lambda]vvsLS=couplingTensor//Normal//Simplify//FullSimplify//ReplaceAll[#,HelpSolveCubicL]&//SparseArray;
	];
		
	(* --- TemporalVector-scalar couplings --- *)	
	If[mode==0,
		GvvssL=EmptyArray[{nv,nv,ns,ns}]
	];
	
	couplingTensor = -Tfac(HabijV+GvvssL);
	HelpList=DeleteDuplicates@Simplify@Flatten[couplingTensor]//Sort;
	HelpListCleaned=CleanedList[HelpList];
	If[Length[HelpListCleaned]<1,
		\[Lambda]KVec=couplingTensor//Normal//Simplify//SparseArray;
		HelpSolveVecL={};
	,
		HelpVarMod=RelationsBVariables3[HelpList];
		HelpSolveVecL=Table[{HelpListCleaned[[a]]->HelpVarMod[[a]]},{a,1,HelpListCleaned//Length}]//Flatten//Simplify;
		\[Lambda]KVec=couplingTensor//Normal//Simplify//ReplaceAll[#,HelpSolveVecL]&//SparseArray;
	];			
			
	If[mode>=3,
		(* --- Scalar sextic couplings --- *)
		VerbosePrint["Calculating Sextic Tensor"];
		
		couplingTensor = Tfac^2*(\[Lambda]6+\[Lambda]6D);
		{HelpListCleaned, HelpVar, HelpVarMod, HelpSolveSextic} = ProcessCoupling[couplingTensor, \[Lambda]6d];	
		\[Lambda]6DS=couplingTensor//Normal//Simplify//ReplaceAll[#,HelpSolveSextic]&//SparseArray;

		If[Length[SparseArray[\[Lambda]6DS]["NonzeroValues"]]!=Length[SparseArray[\[Lambda]6]["NonzeroValues"]],
			Print["Detected 1-loop Scalar Sextic not defined at tree-level"];
			Print["Please Check if you defined all couplings allowed by symmetry"];
		];
	];
	
	(* --- Debye masses --- *);
	If[mode>=2,
		couplingTensor = xLO*aV3D+ xNLO*\[Mu]VabNLO;
	,
		couplingTensor = xLO*aV3D;	
	];
	{HelpListCleaned, HelpVar, HelpVarMod, HelpSolveVectorMass} = ProcessCoupling[couplingTensor, \[Mu]ijV];	
	\[Mu]ijVNLO=couplingTensor//Normal//Simplify//FullSimplify//ReplaceAll[#,HelpSolveVectorMass]&//SparseArray;

	If[DiagonalMatrixQAE[Normal[\[Mu]ijVNLO]]==False,
		Print["Off-Diagonal Debye Matrices Detected"]
	];

	(* --- Scalar masses --- *)
	If[mode>=2,
		RGRunningHardToSoft[]; (*Runs from the hard to the soft scale in the effective theory*)
		couplingTensor = xLO*aS3D+xNLO*(\[Mu]SijNLO+rgFac*Contri\[Beta]SoftToHard);
	,
		couplingTensor = xLO*aS3D;
	];
	{HelpListCleaned, HelpVar, HelpVarMod, HelpSolveMass} = ProcessCoupling[couplingTensor, \[Mu]ijS];	
	\[Mu]ijSNLO=couplingTensor//Normal//Simplify//ReplaceAll[#,HelpSolveMass]&//SparseArray;

	(* --- Scalar tadpoles --- *)
	If[mode>=2,
		couplingTensor = xLO*Tfac^(-1/2)(\[Lambda]1+TadPoleLO)+xNLO(TadPoleNLO*Tfac^(-1/2) +rgFac*ContriTadPoleSoftToHard);
	,
		couplingTensor = xLO*Tfac^(-1/2)(\[Lambda]1+TadPoleLO);
	];
	{HelpListCleaned, HelpVar, HelpVarMod, HelpSolveTadpole} = ProcessCoupling[couplingTensor, dS];	
	TadPoleS=couplingTensor//Normal//Simplify//ReplaceAll[#,HelpSolveTadpole]&//SparseArray;

	If[Length[SparseArray[TadPoleS]["NonzeroValues"]]!=Length[SparseArray[\[Lambda]1]["NonzeroValues"]],
		Print["Detected 1-loop Scalar Tadpoles not defined at tree-level"];
		Print["Please Check if you defined all couplings allowed by symmetry"];
	];

	If[mode>=1,
		IdentMat=List/@Join[
			HelpSolveCubicS,HelpSolveVecT,HelpSolveVecL,
			HelpSolveQuarticL,HelpSolveVectorMass,HelpSolveMass,
			HelpSolveCubicL,HelpSolveTadpole,HelpSolveQuartic
		]/.{b_->a_}:>a->b//Flatten[#,1]&;
	,
		IdentMat=List/@Join[HelpSolveVectorMass,HelpSolveMass,HelpSolveVecL]/.{b_->a_}:>a->b//Flatten[#,1]&;
	];
	
	If[mode>=3,
		IdentMatEff=List/@Join[HelpSolveSextic]/.{b_->a_}:>a->b//Flatten[#,1]&;
	];
];


(* ::Section:: *)
(*Printing the results*)


(* ::Subsection:: *)
(*Print Debye mass*)


(*
    Prints the Debye mass in the soft theory.
    Can be customized with options for LO, NLO, or All.
*)
PrintDebyeMass[optP_:"All"] := Module[
	{
		opt = optP,
		VarGauge, \[Mu]ijp, var, helpMass, SolMassPre, SolMass
	},
    
    (* Verbose output *)
    VerbosePrint["Printing Debye Masses"];
    
    (* Prepare gauge couplings and variables *)
    VarGauge = Join[\[Mu]abDef // Normal // Variables] // DeleteDuplicates;
    
    \[Mu]ijp = \[Mu]abDef // Normal;
    var = Variables[Normal[\[Mu]ijp]];  (* Get variables from the expression *)
    
    (* Calculate the mass term difference *)
    helpMass = Normal[\[Mu]ijp - \[Mu]ijVNLO];
    
    (* Solve for mass and flatten the result *)
    SolMassPre = Solve[helpMass == 0, var] /. IdentMat // Flatten[#, 1] &;
    SolMass = SolMassPre;
    
    (* Apply options to SolMass based on optP *)
    Switch[opt,
        "All", SolMass = SolMassPre /. xLO -> 1 /. xNLO -> 1 /. ReplaceLb // Simplify,
        "LO", SolMass = SolMassPre /. xLO -> 1 /. xNLO -> 0 /. ReplaceLb // Simplify,
        "NLO", SolMass = SolMassPre /. xLO -> 0 /. xNLO -> 1 /. ReplaceLb // Simplify,
        _, Message[PrintResults::badopt, opt]; Return[$Failed];
    ];
    
    (* Return the result in readable form *)
    Return[ToExpression[StringReplace[ToString[StandardForm[Join[SolMass]]], "DRalgo`Private`" -> ""]]]
];


(* ::Subsection:: *)
(*Print higher-dimensional effective couplings*)


(*
	Prints higher-order couplings.
*)
PrintCouplingsEffective[]:=Module[{},
	VerbosePrint["Printing higher-dimension couplings"];

(*Scalar Sextic*);
	VarGauge=Join[\[Lambda]6//Normal//Variables]//DeleteDuplicates;
	SubGauge=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarGauge}];
	\[Lambda]6p=\[Lambda]6//Normal//ReplaceAll[#,SubGauge]&;
	SolVar=\[Lambda]6DS-\[Lambda]6p//Normal;
	SexticVar=\[Lambda]6p//Normal//Variables;
	ResScalp=Reduce[SolVar==0,SexticVar]//ToRules[#]&;
	SolveTemp=SexticVar/.ResScalp;
	ResScal=Table[{SexticVar[[i]]->SolveTemp[[i]]},{i,1,Length@SexticVar}]//Flatten[#,1]&//ReplaceAll[#,IdentMatEff]&//Simplify;


(*Printing Result*)
	PrintPre=Join[ResScal]//Normal//FullSimplify//DeleteDuplicates;

	ToExpression[StringReplace[ToString[StandardForm[PrintPre]],"DRalgo`Private`"->""]]

];


(* ::Subsection:: *)
(*Print effective couplings*)


(*
	Prints effective couplings in the soft theory. The module calculates all tensors in the soft theory and matches them with corresponding tensors in the original 4d theory.
*)
PrintCouplings[]:=Module[
	{
		VarGauge, SubGauge,
		A1, A2
	},
	
	VerbosePrint["Printing 3D vector and quartic couplings in terms of 4D couplings"];

	(*VarGauge=Join[gvvv//Normal//Variables,gvss//Normal//Variables,gvff//Normal//Variables]//DeleteDuplicates;*)
	VarGauge=GaugeCouplingNames//Variables;
	SubGauge=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarGauge}];
	(*VarGauge denotes all possible vector couplings*)

	(* Gauge couplings*)
	If[mode==0,\[Lambda]KVecT=Tfac HabijV]; (*If mode=0 no couplings are calculated*)
	
	A1=TensorContract[\[Lambda]KVecT,{{3,4}}]//Normal; 
	A2=TensorContract[HabijV,{{3,4}}]//Normal//ReplaceAll[#,SubGauge]&;
	
	(* --- Trick to avoid problems when kinetic mixing --- *)
	A1=DiagonalMatrix[Diagonal[A1]];
	A2=DiagonalMatrix[Diagonal[A2]];
	(* --- end of trick --- *)
	
	Var3D=VarGauge//ReplaceAll[#,SubGauge]&//Variables;
	RepVar3D=#->Sqrt[#]&/@Var3D;
	A2Mod=A2/.RepVar3D;
	Sol1=Solve[A2Mod==A1,Var3D]/.IdentMat//Flatten[#,1]&//FullSimplify;
	ResGauge=Table[List[Sol1[[c]]]/.{b_->a_}:>b^2->a,{c,1,Length[Sol1]}];
	

	(* --- Non-Abelian Couplings --- *)
	NonAbelianCoupling[];
	
	If[mode==0,Ggvvv=EmptyArray[{nv,nv,nv}]]; (*If mode=0 no couplings are calculated*)
	GabcdTemp=Ggvvv . gvvv+gvvv . Ggvvv;
	GabVTree=TensorContract[GabcdV,{{2,3}}]//Normal;
	GabVLoop=TensorContract[GabcdTemp,{{2,3}}]//Normal;

	HelpList=DeleteDuplicates@Flatten@FullSimplify[T (GabVTree+GabVLoop) ]//Sort;
	HelpVar=Table[ \[Lambda]VNA[a],{a,1,Rest[HelpList]//Length}];
	HelpVarMod=RelationsBVariables[HelpList,HelpVar];
	HelpSolveNA=Table[{Rest[HelpList][[a]]->ReplaceAll[HelpVarMod[[a]],SubGauge]},{a,1,Rest[HelpList]//Length}]//Flatten;
	\[Lambda]VecNA=Tfac (GabVTree+GabVLoop)//Normal//FullSimplify//ReplaceAll[#,HelpSolveNA]&//SparseArray;
	IdentMatNA=List/@Join[HelpSolveNA]/.{b_->a_}:>a->b//Flatten[#,1]&;

	A1=\[Lambda]VecNA//Normal;
	A2=GabVTree//Normal//ReplaceAll[#,SubGauge]&;
	Var3D=A2//Variables;
	RepVar3D=#->Sqrt[#]&/@Var3D;
	A2Mod=A2/.RepVar3D;
	Sol1=Solve[A2Mod==A1,Var3D]/.IdentMatNA//Flatten[#,1]&//FullSimplify;
	ResGaugeNA=Table[List[Sol1[[c]]]/.{b_->a_}:>b^2->a,{c,1,Length[Sol1]}]//Simplify;


	(* --- Scalar quartics --- *)
	(*If mode=0 no couplings are calculated*)
	VarGauge=Join[\[Lambda]4//Normal//Variables]//DeleteDuplicates;
	SubGauge=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarGauge}];
	
	If[mode==0,
			\[Lambda]4p= \[Lambda]4//Normal//ReplaceAll[#,SubGauge]&;
			SolVar=Tfac*\[Lambda]4-\[Lambda]4p//Normal;
			QuarticVar=\[Lambda]4p//Normal//Variables;
			ResQuartic=Solve[SolVar==0,QuarticVar]/.IdentMat//Flatten[#,1]&;
	,
			SolVar=SparseArray[\[Lambda]3DS-\[Lambda]4]["NonzeroValues"]//DeleteDuplicates//ReplaceAll[#,SubGauge]&;
			QuarticVar=\[Lambda]4//Normal//ReplaceAll[#,SubGauge]&//Variables;
			ResScalp=Reduce[SolVar==0,QuarticVar]//ToRules[#]&;
			SolveTemp=QuarticVar/.ResScalp;
			ResQuartic=Table[{QuarticVar[[i]]->SolveTemp[[i]]},{i,1,Length@QuarticVar}]//Flatten[#,1]&//ReplaceAll[#,IdentMat]&//Simplify;
	]; 
	

	(* --- Scalar Cubics --- *)

	VarGauge=Join[\[Lambda]3//Normal//Variables]//DeleteDuplicates;
	SubGauge=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarGauge}];
	
	If[mode==0,\[Lambda]3CSRed=Tfac^(1/2) \[Lambda]3]; (*If mode=0 no couplings are calculated*)
	\[Lambda]3p=\[Lambda]3//Normal//ReplaceAll[#,SubGauge]&;
	SolVar=\[Lambda]3CSRed-\[Lambda]3p//Normal;
	CubicVar=\[Lambda]3p//Normal//Variables;
	ResCubic=Solve[SolVar==0,CubicVar]/.IdentMat//Flatten[#,1]&;

	(* --- Printing Result --- *)
	PrintPre=Join[ResGauge,ResGaugeNA,ResQuartic,ResCubic]//Normal//FullSimplify//DeleteDuplicates;

	ToExpression[StringReplace[ToString[StandardForm[PrintPre]],"DRalgo`Private`"->""]]

];


(* ::Subsection::Closed:: *)
(*Print symmetric pressure*)


(*
    Prints the result from SymmEnergy.
    Optional argument selects the loop order: "LO", "NLO", or "NNLO".
*)
PrintPressure[optP_: "All"] := Module[{opt = optP, SymmPrint},

    If[mode > 0,
        SymmPrint = Switch[opt,
            "LO", SymmEnergy[[1]],
            "NLO", SymmEnergy[[2]],
            "NNLO", SymmEnergy[[3]],
            "All", Total[SymmEnergy[[1 ;; 3]]],
            _, Message[PrintResults::badopt, opt]; Return[$Failed];
        ],
        (* Symmetric phase when mode \[LessEqual] 0 *)
        SymmPrint = -SymmetricPhaseLO[];
    ];

    (* Format and return output *)
    ToExpression[
        StringReplace[
            ToString[StandardForm[SymmPrint]],
            "DRalgo`Private`" -> ""
        ]
    ]
];


(* ::Subsection:: *)
(*Print scalar mass*)


(*
    Prints 1-loop and 2-loop effective scalar masses in the soft theory.
    Optional argument: "LO", "NLO", or "All" (default: "All").
*)
PrintScalarMass[optP_: "All"] := Module[
    {opt = optP, VarGauge, SubGauge, \[Mu]ijp, var, helpMass, ResScalp, SolveTemp, SolMassPre, SolMass},

    VerbosePrint["Printing Scalar Masses"];

    (* Collect gauge variables and apply 3d replacements *)
    VarGauge = DeleteDuplicates[Variables[Normal[\[Mu]ij]]];
    SubGauge = Rule @@@ (List @@@ {#, Symbol[ToString[#] <> "3d"]} & /@ VarGauge);

    \[Mu]ijp = Normal[\[Mu]ij] /. SubGauge;
    var = Variables[Normal[\[Mu]ijp]];
    helpMass = Normal[\[Mu]ijp - \[Mu]ijSNLO];

    (* Solve for scalar masses *)
    ResScalp = ToRules[Reduce[helpMass == 0, var]];
    SolveTemp = var /. ResScalp;
    SolMassPre = ReplaceAll[Flatten[Table[{var[[i]] -> SolveTemp[[i]]}, {i, Length[var]}]], IdentMat];

    (* Apply loop order option *)
    SolMass = Switch[opt,
        "LO",    SolMassPre /. xLO -> 1 /. xNLO -> 0 /. ReplaceLb // Simplify,
        "NLO",   SolMassPre /. xLO -> 0 /. xNLO -> 1 /. ReplaceLb // Simplify,
        "All",   SolMassPre /. xLO -> 1 /. xNLO -> 1 /. ReplaceLb // Simplify,
        _,      Message[PrintResults::badopt, opt]; Return[$Failed];
    ];

    (* Final print with cleanup *)
    ToExpression[
        StringReplace[
            ToString[StandardForm[Join[SolMass]]],
            "DRalgo`Private`" -> ""
        ]
    ]
];


(* ::Subsection:: *)
(*Print tadpole*)


(*
    Prints 1-loop and 2-loop effective tadpoles in the soft theory.
    Valid options: "LO", "NLO", or "All". Default: "All".
*)
PrintTadpoles[optP_: "All"] := Module[
    {
    opt = optP,
    VarGauge, SubGauge, \[Lambda]1p, var, helpMass, ResScalp, SolveTemp, SolMassPre, SolTadpole
    },

    VerbosePrint["Printing Tadpoles"];

    (* Gauge substitutions for 3D theory *)
    VarGauge = DeleteDuplicates[Variables[Normal[\[Lambda]1]]];
    SubGauge = Rule @@@ (List @@@ {#, Symbol[ToString[#] <> "3d"]} & /@ VarGauge);

    (* Setup and solve tadpole equations *)
    \[Lambda]1p = Normal[\[Lambda]1] /. SubGauge;
    var = Variables[Normal[\[Lambda]1p]];
    helpMass = Normal[\[Lambda]1p - TadPoleS];
    ResScalp = ToRules[Reduce[helpMass == 0, var]];
    SolveTemp = var /. ResScalp;
    SolMassPre = ReplaceAll[Flatten[Table[{var[[i]] -> SolveTemp[[i]]}, {i, Length[var]}]], IdentMat];

    (* Loop-order selection *)
    SolTadpole = Switch[opt,
        "LO",   SolMassPre /. xLO -> 1 /. xNLO -> 0,
        "NLO",  SolMassPre /. xLO -> 0 /. xNLO -> 1,
        "All",  SolMassPre /. xLO -> 1 /. xNLO -> 1,
        _,      Message[PrintResults::badopt, opt]; Return[$Failed];
    ];

    (* Output final expression *)
    ToExpression[
        StringReplace[
            ToString[StandardForm[Join[SolTadpole]]],
            "DRalgo`Private`" -> ""
        ]
    ]
];


(*
	Rewrites internal names of effective couplings in terms of 4d parameters.
*)
PrintIdentification[]:=Module[{},
	ToExpression[StringReplace[ToString[StandardForm[IdentMat]],"DRalgo`Private`"->""]]
];


(* ::Subsection::Closed:: *)
(*Print anomalous dimensions*)


(*
	Prints anomalous dimensions of 4d fields.
*)
AnomDim4D[ParticleI_,ComponentsI_]:=Module[{ParticleP=ParticleI,ComponentsP=ComponentsI},
	If[mode>=2,
		CounterTerm[];(*Calculates counterterms if not already done*)
	,
		\[Gamma]ij=EmptyArray[{ns,ns}];
		\[Gamma]IJF=EmptyArray[{nf,nf}];
		\[Gamma]ab=EmptyArray[{nv,nv}];
		
	]; 
	Switch[ParticleP,"S",
		Ret=\[Gamma]ij[[ComponentsP[[1]][[1]],ComponentsP[[2]][[1]]]];
	,"F",
		Ret=\[Gamma]IJF[[ComponentsP[[1]][[1]],ComponentsP[[2]][[1]]]];
	,"V",
		Ret=\[Gamma]ab[[ComponentsP[[1]][[1]],ComponentsP[[2]][[1]]]];
	];

	ToExpression[StringReplace[ToString[StandardForm[Ret]],"DRalgo`Private`"->""]]
];


(*
	Finds the \[Epsilon]^0 coefficient of an expression
*)
SerEnergyHelp[optP_]:=Module[{opt=optP},
	HelpSE=Series[opt/.D->4-2\[Epsilon]/.\[Epsilon]bp->(1/\[Epsilon]+Lb)^-1/.\[Epsilon]b->(1/\[Epsilon]+Lbb)^-1/.\[Epsilon]BF->(1/\[Epsilon]+LBF)^-1/.\[Epsilon]F->(1/\[Epsilon]+LFF)^-1/.\[Epsilon]FB->(1/\[Epsilon]+LFB)^-1/.\[Epsilon]bbM->(1/\[Epsilon]+LbbM)^-1,{\[Epsilon],0,0}]//Normal;
	Return[Coefficient[HelpSE//Normal//Expand,\[Epsilon],0]]
];


(* ::Subsection:: *)
(*Print temporal scalar couplings*)


(*
	Prints temporal couplings
*)
PrintTemporalScalarCouplings[]:=Module[
	{
	quarticVars, cubicVars, quarticTempVars,
	quarticList, cubicList, quarticTempList,
	temporalCouplings, AE
	},
	
	VerbosePrint["Printing Temporal Vector Couplings"];
	If[verbose,
		Print[
			"\[Lambda]VLL[a] are \!\(\*SuperscriptBox[\(V\), \(4\)]\) Couplings, ",
			"\[Lambda]vvsLS[a] are \!\(\*SuperscriptBox[\(V\), \(2\)]\)S Couplings, ",
			"\[Lambda]VL[a] are \!\(\*SuperscriptBox[\(V\), \(2\)]\)\!\(\*SuperscriptBox[\(S\), \(2\)]\) Couplings."
		];
	];

	(* Quartic: V^4 *)
	If[mode==0,\[Lambda]AAS=EmptyArray[{nv,nv,nv,nv}]]; (*If mode=0 no couplings are calculated*)
	quarticVars=\[Lambda]AAS//Normal//Variables;
	AE=ToExpression[StringReplace[ToString[StandardForm[quarticVars]],"DRalgo`Private`"->""]];
	quarticList=#->(ReplaceAll[#,PrintIdentification[]])&/@AE;

	(* Cubic: V^2S *)
	If[mode==0,\[Lambda]vvsLS=EmptyArray[{nv,nv,ns}]]; (*If mode=0 no couplings are calculated*)
	cubicVars=\[Lambda]vvsLS//Normal//Variables;
	AE=ToExpression[StringReplace[ToString[StandardForm[cubicVars]],"DRalgo`Private`"->""]];
	cubicList=#->(ReplaceAll[#,PrintIdentification[]])&/@AE;

	(* Mixed: V^2S^2 *)
	quarticTempVars=\[Lambda]KVec//Normal//Variables;
	AE=ToExpression[StringReplace[ToString[StandardForm[quarticTempVars]],"DRalgo`Private`"->""]];
	quarticTempList=#->(ReplaceAll[#,PrintIdentification[]])&/@AE;

	(* Combine and print all temporal couplings *)
	temporalCouplings=Join[quarticList, cubicList, quarticTempList]//Flatten;
	
	ToExpression[StringReplace[
		ToString[StandardForm[temporalCouplings]],
		"DRalgo`Private`"->""]]
];


errPrintTens="Order of Tensors:\n\
	(1) Scalar Quartic,\n\
	(2) Vector-Scalar Couplings,\n\
	(3) Temporal-Vector Scalar Couplings,\n\
	(4) Temporal-Vector Quartics,\n\
	(5) Temporal-Vector Mass,\n\
	(6) Scalar Mass,\n\
	(7) Scalar Cubics,\n\
	(8) Temporal Vector-Scalar Cubics,\n\
	(9) Tadpoles";


PrintTensorDRalgo[]:="nothing"/;Print[errPrintTens];


PrintTensorDRalgo[ind_]:=Module[{},	
	
	Switch[ind, 
    1, Return[OutputFormatDR[\[Lambda]3DS]], 
    2, Return[OutputFormatDR[\[Lambda]KVecT]], 
    3, Return[OutputFormatDR[\[Lambda]KVec]], 
    4, Return[OutputFormatDR[\[Lambda]AAS]], 
    5, Return[OutputFormatDR[\[Mu]ijVNLO]], 
    6, Return[OutputFormatDR[\[Mu]ijSNLO]], 
    7, Return[OutputFormatDR[\[Lambda]3CSRed]], 
    8, Return[OutputFormatDR[\[Lambda]vvsLS]],   
    9, Return[OutputFormatDR[TadPoleS]],    
    _, 
    Print[errPrintTens];
  ];
      Return[""]

];


(* ::Subsection:: *)
(*Print 4D counter terms*)


(*
	Prints Counterterms of 4d couplings and masses.
*)
CounterTerms4D[]:=Module[{},
	VerbosePrint["Finding \[Beta]-functions"];
	CounterTerm[];

(*To make the comparisons easier*)
		VarGauge=GaugeCouplingNames//Variables;
		SubGauge=Table[c->Symbol[ToString[c]<>ToString["temp"]],{c,VarGauge}];
		SubGauge2=Table[Symbol[ToString[c]<>ToString["temp"]]->c,{c,VarGauge}];

(* 
	Gauge couplings
*)
(*Anomalous-dimension contributions*)

		A1=TensorContract[\[Beta]vvss//Normal,{{3,4}}];
		A2=TensorContract[HabijV,{{3,4}}]//Normal//ReplaceAll[#,SubGauge]&;
(*Trick to avoid problems when kinetic mixing*)
		A1=DiagonalMatrix[Diagonal[A1]];
		A2=DiagonalMatrix[Diagonal[A2]];
(*end of trick*)
		Var3D=VarGauge//ReplaceAll[#,SubGauge]&//Variables;
		RepVar3D=#->Sqrt[#]&/@Var3D;
		A2Mod=A2/.RepVar3D;

		Sol1=Solve[A2Mod==A1,Var3D]//Flatten[#,1]&//FullSimplify;
		ResGauge=Table[List[Sol1[[c]]]/.{b_->a_}:>b^2->a,{c,1,Length[Sol1]}]/.SubGauge2;

(* 
	Non-abelian couplings
*)

		GabcdTemp=(\[Beta]gvvv . gvvv+gvvv . \[Beta]gvvv)//SparseArray;
		GabVTree=TensorContract[GabcdV,{{2,3}}]//Normal;
		GabVLoop=TensorContract[GabcdTemp,{{2,3}}]//Normal;

		A1=GabVLoop//Normal;
		A2=GabVTree//Normal//ReplaceAll[#,SubGauge]&;
		Var3D=A2//Variables;
		RepVar3D=#->Sqrt[#]&/@Var3D;
		A2Mod=A2/.RepVar3D;
		Sol1=Solve[A2Mod==A1,Var3D]//Flatten[#,1]&//FullSimplify;
		ResGaugeNA=Table[List[Sol1[[c]]]/.{b_->a_}:>b^2->a,{c,1,Length[Sol1]}]/.SubGauge2//Simplify;


(* 
	Scalar-quartic couplings
*)

		VarGauge=Join[\[Lambda]4//Normal//Variables]//DeleteDuplicates;
		SubGauge=Table[c->Symbol[ToString[c]<>ToString["temp"]],{c,VarGauge}];
		SubGauge2=Table[Symbol[ToString[c]<>ToString["temp"]]->c,{c,VarGauge}];


		HelpList=DeleteDuplicates@Flatten@SimplifySparse[\[Beta]\[Lambda]ijkl]//Sort;
		HelpVarMod=RelationsBVariables3[HelpList]//ReplaceAll[#,\[Lambda]VL[v1_]->\[Lambda]Beta[v1]]&;
		If[HelpList[[1]]==0&&Length[HelpList]>1,
				HelpList=Rest[HelpList];
		];
		HelpSolve\[Beta]=Table[{HelpList[[a]]->HelpVarMod[[a]]},{a,1,HelpList//Length}]//Flatten//Simplify;
		HelpSolve\[Beta]2=Table[{HelpVarMod[[a]]->HelpList[[a]]},{a,1,HelpList//Length}]//Flatten//Simplify;
		\[Lambda]4\[Beta]=\[Beta]\[Lambda]ijkl//SimplifySparse//Normal//ReplaceAll[#,HelpSolve\[Beta]]&//SparseArray;

		\[Lambda]4p=\[Lambda]4//Normal//ReplaceAll[#,SubGauge]&;
		SolVar=\[Lambda]4\[Beta]-\[Lambda]4p//Normal;
		QuarticVar=\[Lambda]4p//Normal//Variables;
		ResScalp=Reduce[SolVar==0,QuarticVar]//ToRules[#]&;
		SolveTemp=QuarticVar/.ResScalp;
		ResScal=Table[{QuarticVar[[i]]->SolveTemp[[i]]},{i,1,Length@QuarticVar}]/.SubGauge2//Flatten[#,1]&//ReplaceAll[#,HelpSolve\[Beta]2]&//Simplify;

(* 
	Scalar-cubic couplings
*)
		VarGauge=Join[\[Lambda]3//Normal//Variables]//DeleteDuplicates;
		SubGauge=Table[c->Symbol[ToString[c]<>ToString["temp"]],{c,VarGauge}];
		SubGauge2=Table[Symbol[ToString[c]<>ToString["temp"]]->c,{c,VarGauge}];

		\[Lambda]3p=\[Lambda]3//Normal//ReplaceAll[#,SubGauge]&;
		SolVar=\[Beta]\[Lambda]ijk-\[Lambda]3p//Normal;
		CubicVar=\[Lambda]3p//Normal//Variables;
		ResCubic=Solve[SolVar==0,CubicVar]/.SubGauge2//Flatten[#,1]&;

(* 
	Yukawa couplings
*)
		VarGauge=Join[Ysff//Normal//Variables]//DeleteDuplicates;
		SubGauge=Table[c->Symbol[ToString[c]<>ToString["temp"]],{c,VarGauge}];
		SubGauge2=Table[Symbol[ToString[c]<>ToString["temp"]]->c,{c,VarGauge}];

		\[Lambda]4p=Ysff//Normal//ReplaceAll[#,SubGauge]&;
		SolVar=\[Beta]Ysij-\[Lambda]4p//Normal;
		QuarticVar=\[Lambda]4p//Normal//Variables;
		ResYuk=Solve[SolVar==0,QuarticVar]/.SubGauge2//Flatten[#,1]&;



(* 
	Scalar masses
*)

		VarGauge=Join[\[Mu]ij//Normal//Variables]//DeleteDuplicates;
		SubGauge=Table[c->Symbol[ToString[c]<>ToString["temp"]],{c,VarGauge}];
		SubGauge2=Table[Symbol[ToString[c]<>ToString["temp"]]->c,{c,VarGauge}];

		\[Lambda]4p=\[Mu]ij//Normal//ReplaceAll[#,SubGauge]&;
		SolVar=\[Beta]mij-\[Lambda]4p//Normal;
		QuarticVar=\[Lambda]4p//Normal//Variables;
		ResMass=Solve[SolVar==0,QuarticVar]/.SubGauge2//Flatten[#,1]&;

(* 
	Fermion masses
*)

		VarGauge=Join[\[Mu]IJF//Normal//Variables]//DeleteDuplicates;
		SubGauge=Table[c->Symbol[ToString[c]<>ToString["temp"]],{c,VarGauge}];
		SubGauge2=Table[Symbol[ToString[c]<>ToString["temp"]]->c,{c,VarGauge}];

		\[Lambda]4p=\[Mu]IJF//Normal//ReplaceAll[#,SubGauge]&;
		SolVar=\[Beta]\[Mu]IJF-\[Lambda]4p//Normal;
		QuarticVar=\[Lambda]4p//Normal//Variables;
		ResMassF=Solve[SolVar==0,QuarticVar]/.SubGauge2//Flatten[#,1]&;



(* 
	Scalar tadpoles
*)
(*Scalar Mass insertion*)
		Contri4=1/2*I2p*Activate@TensorContract[Inactive@TensorProduct[\[Lambda]3,\[Mu]ij],{{2,4},{3,5}}];


(*Fermion Mass insertion*)
		VarGauge=Join[\[Lambda]1//Normal//Variables]//DeleteDuplicates;
		SubGauge=Table[c->Symbol[ToString[c]<>ToString["temp"]],{c,VarGauge}];
		SubGauge2=Table[Symbol[ToString[c]<>ToString["temp"]]->c,{c,VarGauge}];

		\[Lambda]4p=\[Lambda]1//Normal//ReplaceAll[#,SubGauge]&;
		SolVar=\[Beta]\[Lambda]1-\[Lambda]4p//Normal;
		QuarticVar=\[Lambda]4p//Normal//Variables;
		ResTadpole=Solve[SolVar==0,QuarticVar]/.SubGauge2//Flatten[#,1]&;

(*Effective couplings*)
		If[mode>=3,

(* 
	Scalar-cubic couplings
*)

			VarGauge=Join[\[Lambda]6//Normal//Variables]//DeleteDuplicates;
			SubGauge=Table[c->Symbol[ToString[c]<>ToString["temp"]],{c,VarGauge}];
			SubGauge2=Table[Symbol[ToString[c]<>ToString["temp"]]->c,{c,VarGauge}];

			\[Lambda]6p=\[Lambda]6//Normal//ReplaceAll[#,SubGauge]&;
			SolVar=\[Beta]\[Lambda]6ijklnm-\[Lambda]6p//Normal;
			SexticVar=\[Lambda]6p//Normal//Variables;
			ResSextic=Solve[SolVar==0,SexticVar]/.SubGauge2//Flatten[#,1]&;



(*Printing Result*)
			PrintPre=Join[ResGauge,ResGaugeNA,ResScal,ResCubic,ResYuk,ResMass,ResMassF,ResTadpole,ResSextic]//Normal//FullSimplify//DeleteDuplicates;
		,
(*Printing Result*)
			PrintPre=Join[ResGauge,ResGaugeNA,ResScal,ResCubic,ResYuk,ResMass,ResMassF,ResTadpole]//Normal//FullSimplify//DeleteDuplicates;
		];

		PrintPre=ReplaceAll[PrintPre,(a_->b_):>a->b/(2 \[Epsilon])];
		ToExpression[StringReplace[ToString[StandardForm[PrintPre]],"DRalgo`Private`"->""]]

];


(* ::Subsection:: *)
(*Print 4D beta functions*)


(*
	Prints beta-functions of 4d couplings and masses.
*)
BetaFunctions4D[]:=Module[{},
	VerbosePrint["Finding \[Beta]-functions"];
	If[mode>=2,
		CounterTerm[];


(*To make the comparisons easier*)
		VarGauge=GaugeCouplingNames//Variables;
		SubGauge=Table[c->Symbol[ToString[c]<>ToString["temp"]],{c,VarGauge}];
		SubGauge2=Table[Symbol[ToString[c]<>ToString["temp"]]->c,{c,VarGauge}];

(* 
	Gauge couplings
*)
		A1=TensorContract[\[Beta]vvss//Normal,{{3,4}}];
		A2=TensorContract[HabijV,{{3,4}}]//Normal//ReplaceAll[#,SubGauge]&;
(*Trick to avoid problems when kinetic mixing*)
		A1=DiagonalMatrix[Diagonal[A1]];
		A2=DiagonalMatrix[Diagonal[A2]];
(*end of trick*)
		Var3D=VarGauge//ReplaceAll[#,SubGauge]&//Variables;
		RepVar3D=#->Sqrt[#]&/@Var3D;
		A2Mod=A2/.RepVar3D;

		Sol1=Solve[A2Mod==A1,Var3D]//Flatten[#,1]&//FullSimplify;
		ResGauge=Table[List[Sol1[[c]]]/.{b_->a_}:>b^2->a,{c,1,Length[Sol1]}]/.SubGauge2;

(* 
	Non-abelian couplings
*)
		GabcdTemp=\[Beta]gvvv . gvvv+gvvv . \[Beta]gvvv//SparseArray;
		GabVTree=TensorContract[GabcdV,{{2,3}}]//Normal;
		GabVLoop=TensorContract[GabcdTemp,{{2,3}}]//Normal;

		A1=GabVLoop//Normal;
		A2=GabVTree//Normal//ReplaceAll[#,SubGauge]&;
		Var3D=A2//Variables;
		RepVar3D=#->Sqrt[#]&/@Var3D;
		A2Mod=A2/.RepVar3D;
		Sol1=Solve[A2Mod==A1,Var3D]//Flatten[#,1]&//FullSimplify;
		ResGaugeNA=Table[List[Sol1[[c]]]/.{b_->a_}:>b^2->a,{c,1,Length[Sol1]}]/.SubGauge2//Simplify;

(* 
	Scalar-quartic couplings
*)
		VarGauge=Join[\[Lambda]4//Normal//Variables]//DeleteDuplicates;
		SubGauge=Table[c->Symbol[ToString[c]<>ToString["temp"]],{c,VarGauge}];
		SubGauge2=Table[Symbol[ToString[c]<>ToString["temp"]]->c,{c,VarGauge}];

	
		HelpList=DeleteDuplicates@Flatten@SimplifySparse[\[Beta]\[Lambda]ijkl]//Sort;
		HelpVarMod=RelationsBVariables3[HelpList]//ReplaceAll[#,\[Lambda]VL[v1_]->\[Lambda]Beta[v1]]&;

		If[HelpList[[1]]==0&&Length[HelpList]>1,
			HelpList=Rest[HelpList];
		];
		HelpSolve\[Beta]=Table[{HelpList[[a]]->HelpVarMod[[a]]},{a,1,HelpList//Length}]//Flatten//Simplify;
		HelpSolve\[Beta]2=Table[{HelpVarMod[[a]]->HelpList[[a]]},{a,1,HelpList//Length}]//Flatten//Simplify;
		\[Lambda]4\[Beta]=\[Beta]\[Lambda]ijkl//SimplifySparse//Normal//ReplaceAll[#,HelpSolve\[Beta]]&//SparseArray;

		\[Lambda]4p=\[Lambda]4//Normal//ReplaceAll[#,SubGauge]&;
		SolVar=\[Lambda]4\[Beta]-\[Lambda]4p//Normal;
		QuarticVar=\[Lambda]4p//Normal//Variables;
		ResScalp=Reduce[SolVar==0,QuarticVar]//ToRules[#]&;
		SolveTemp=QuarticVar/.ResScalp;
		ResScal=Table[{QuarticVar[[i]]->SolveTemp[[i]]},{i,1,Length@QuarticVar}]/.SubGauge2//Flatten[#,1]&//ReplaceAll[#,HelpSolve\[Beta]2]&//Simplify;

(* 
	Scalar-cubic couplings
*)
		VarGauge=Join[\[Lambda]3//Normal//Variables]//DeleteDuplicates;
		SubGauge=Table[c->Symbol[ToString[c]<>ToString["temp"]],{c,VarGauge}];
		SubGauge2=Table[Symbol[ToString[c]<>ToString["temp"]]->c,{c,VarGauge}];

	
		\[Lambda]3p=\[Lambda]3//Normal//ReplaceAll[#,SubGauge]&;
		SolVar=\[Beta]\[Lambda]ijk-\[Lambda]3p//Normal;
		CubicVar=\[Lambda]3p//Normal//Variables;
		ResCubic=Solve[SolVar==0,CubicVar]/.SubGauge2//Flatten[#,1]&;

(* 
	Yukawa couplings
*)
		VarGauge=Join[Ysff//Normal//Variables]//DeleteDuplicates;
		SubGauge=Table[c->Symbol[ToString[c]<>ToString["temp"]],{c,VarGauge}];
		SubGauge2=Table[Symbol[ToString[c]<>ToString["temp"]]->c,{c,VarGauge}];

	
		\[Lambda]4p=Ysff//Normal//ReplaceAll[#,SubGauge]&;
		SolVar=\[Beta]Ysij-\[Lambda]4p//Normal;
		QuarticVar=\[Lambda]4p//Normal//Variables;
		ResYuk=Solve[SolVar==0,QuarticVar]/.SubGauge2//Flatten[#,1]&;

(* 
	Scalar masses
*)
		VarGauge=Join[\[Mu]ij//Normal//Variables]//DeleteDuplicates;
		SubGauge=Table[c->Symbol[ToString[c]<>ToString["temp"]],{c,VarGauge}];
		SubGauge2=Table[Symbol[ToString[c]<>ToString["temp"]]->c,{c,VarGauge}];

		\[Lambda]4p=\[Mu]ij//Normal//ReplaceAll[#,SubGauge]&;
		SolVar=\[Beta]mij-\[Lambda]4p//Normal;
		QuarticVar=\[Lambda]4p//Normal//Variables;
		ResMass=Solve[SolVar==0,QuarticVar]/.SubGauge2//Flatten[#,1]&;

(* 
	Fermion masses
*)
		VarGauge=Join[\[Mu]IJF//Normal//Variables]//DeleteDuplicates;
		SubGauge=Table[c->Symbol[ToString[c]<>ToString["temp"]],{c,VarGauge}];
		SubGauge2=Table[Symbol[ToString[c]<>ToString["temp"]]->c,{c,VarGauge}];

	
		\[Lambda]4p=\[Mu]IJF//Normal//ReplaceAll[#,SubGauge]&;
		SolVar=\[Beta]\[Mu]IJF-\[Lambda]4p//Normal;
		QuarticVar=\[Lambda]4p//Normal//Variables;
		ResMassF=Solve[SolVar==0,QuarticVar]/.SubGauge2//Flatten[#,1]&;

(* 
	Scalar tadpoles
*)
		If[Length[\[Lambda]1//Normal//Variables]==0,
			ResTadpole={};
		,
			VarGauge=Join[\[Lambda]1//Normal//Variables]//DeleteDuplicates;
			SubGauge=Table[c->Symbol[ToString[c]<>ToString["temp"]],{c,VarGauge}];
			SubGauge2=Table[Symbol[ToString[c]<>ToString["temp"]]->c,{c,VarGauge}];

			\[Lambda]4p=\[Lambda]1//Normal//ReplaceAll[#,SubGauge]&;
			SolVar=\[Beta]\[Lambda]1-\[Lambda]4p//Normal;
			QuarticVar=\[Lambda]4p//Normal//Variables;
			ResTadpole=Solve[SolVar==0,QuarticVar]/.SubGauge2//Flatten[#,1]&;
		];
(*Dim-6 couplings*)
	If[mode>=3,
		VarGauge=Join[\[Lambda]6//Normal//Variables]//DeleteDuplicates;
		SubGauge=Table[c->Symbol[ToString[c]<>ToString["temp"]],{c,VarGauge}];
		SubGauge2=Table[Symbol[ToString[c]<>ToString["temp"]]->c,{c,VarGauge}];
		
		\[Lambda]6p=\[Lambda]6//Normal//ReplaceAll[#,SubGauge]&;
		SolVar=\[Beta]\[Lambda]6ijklnm-\[Lambda]6p//Normal;
		SexticVar=\[Lambda]6p//Normal//Variables;
		ResSextic=Solve[SolVar==0,SexticVar]/.SubGauge2//Flatten[#,1]&;
(*Printing Result*)
		PrintPre=Join[ResGauge,ResGaugeNA,ResScal,ResCubic,ResYuk,ResMass,ResMassF,ResTadpole,ResSextic]//Normal//FullSimplify//DeleteDuplicates;
	,
(*Printing Result*)
		PrintPre=Join[ResGauge,ResGaugeNA,ResScal,ResCubic,ResYuk,ResMass,ResMassF,ResTadpole]//Normal//FullSimplify//DeleteDuplicates;
	];
,
	ResGauge=#^2->0&/@GaugeCouplingNames;
	ResScal=#->0&/@Variables[\[Lambda]4//Normal];
	ResCubic=#->0&/@Variables[\[Lambda]3//Normal];
	ResYuk=#->0&/@Variables[Ysff//Normal];
	ResMass=#->0&/@Variables[\[Mu]ij//Normal];
	ResMassF=#->0&/@Variables[\[Mu]IJF//Normal];
	ResTadpole=#->0&/@Variables[\[Lambda]1//Normal];
	PrintPre=Join[ResGauge,ResScal,ResCubic,ResYuk,ResMass,ResMassF,ResTadpole]//Normal//FullSimplify//DeleteDuplicates;
];
ToExpression[StringReplace[ToString[StandardForm[PrintPre]],"DRalgo`Private`"->""]]


];


(* ::Subsection:: *)
(*Print 3D beta functions*)


(*
	Prints beta-functions of the soft masses and tadpoles.
*)
BetaFunctions3DS[]:=Module[{},
	VerbosePrint["Finding \[Beta]-functions"];
	If[mode>0,

(* 
	Scalar masses
*)
		VarGauge=Join[\[Mu]ij//Normal//Variables]//DeleteDuplicates;
		SubGauge=Table[c->Symbol[ToString[c]<>ToString["temp"]],{c,VarGauge}];
		SubGauge2=Table[Symbol[ToString[c]<>ToString["temp"]]->c,{c,VarGauge}];

		\[Lambda]4p=\[Mu]ij//Normal//ReplaceAll[#,SubGauge]&;
		SolVar=Contri\[Beta]SoftToHard-\[Lambda]4p//Normal;
		QuarticVar=\[Lambda]4p//Normal//Variables;
		ResMass=Solve[SolVar==0,QuarticVar]/.SubGauge2//Flatten[#,1]&;
(* 
	Scalar tadpoles
*)
		If[Length[\[Lambda]1//Normal//Variables]==0,
			ResTadpole={};
		,
			VarGauge=Join[\[Lambda]1//Normal//Variables]//DeleteDuplicates;
			SubGauge=Table[c->Symbol[ToString[c]<>ToString["temp"]],{c,VarGauge}];
			SubGauge2=Table[Symbol[ToString[c]<>ToString["temp"]]->c,{c,VarGauge}];

			\[Lambda]4p=\[Lambda]1//Normal//ReplaceAll[#,SubGauge]&;
			SolVar=ContriTadPoleSoftToHard-\[Lambda]4p//Normal;
			QuarticVar=\[Lambda]4p//Normal//Variables;
			ResTadpole=Solve[SolVar==0,QuarticVar]/.SubGauge2//Flatten[#,1]&;
		];

		PrintPre=Join[ResMass,ResTadpole]/.\[Mu]3->\[Mu] Exp[1]//Normal//FullSimplify//DeleteDuplicates;
,
	ResMass=#->0&/@Variables[\[Mu]ij//Normal];
	ResTadpole=#->0&/@Variables[\[Lambda]1//Normal];
	PrintPre=Join[ResMass,ResTadpole]/.\[Mu]3->\[Mu] Exp[1]//Normal//FullSimplify//DeleteDuplicates;
];
ToExpression[StringReplace[ToString[StandardForm[PrintPre]],"DRalgo`Private`"->""]]


];
