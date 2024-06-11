(* ::Package:: *)

(* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ *)

(* :Title: EffPot                                                          	*)

(*
      	This software is covered by the GNU General Public License 3.
       	Copyright (C) 2021-2023 Andreas Ekstedt
       	Copyright (C) 2021-2023 Philipp Schicho
       	Copyright (C) 2021-2023 Tuomas V.I. Tenkanen

*)

(* :Summary:	Compute the effective potential.        					*)	

(* ------------------------------------------------------------------------ *)


(* ::Section::Closed:: *)
(*Definition of model*)


(*
	User-defined model.
*)
DefineNewTensorsUS[\[Mu]ijI_,\[Lambda]4I_,\[Lambda]3I_,gvssI_,gvvvI_]:=Module[{\[Mu]ijP=\[Mu]ijI,\[Lambda]4P=\[Lambda]4I,\[Lambda]3P=\[Lambda]3I,gvssP=gvssI,gvvvP=gvvvI},
	\[Mu]ijEP=\[Mu]ijP//SparseArray;
	gvvvEP=gvvvP//SparseArray;
	gvssEP=gvssP//SparseArray;
	\[Lambda]4EP=\[Lambda]4P//SparseArray;
	\[Lambda]3EP=\[Lambda]3P//SparseArray;
	nsEP=Length[gvssEP[[1]]];
	nvEP=Length[gvvvEP];

];


(*
	Defines parameters directly that are used by the effective potential
*)
UseSoftTheory[]:=Module[{},
(*Here the user wants to use the soft-theory to calculate the effective potential*)
(*We then need to make sure that all couplings are defined in the soft theory,
	and that temporal-vector bosons are treated as scalars*)

(*We only need to change the gauge coupling name*)
	SubGauge=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,GaugeCouplingNames}];
	gvvvEP=gvvv//Normal//ReplaceAll[#,SubGauge]&//SparseArray;

(*We first change name of the tree-level quartic couplings*)
	VarQuartic=Join[\[Lambda]4//Normal//Variables]//DeleteDuplicates;
	SubQuartic=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarQuartic}];
	\[Lambda]4Temp= \[Lambda]4//Normal//ReplaceAll[#,SubQuartic]&//SparseArray;
	
	(*We then combine normal and temporal scalars*)
	ns=(\[Lambda]4//Length);
	nSH=nv+ns; (*Adds normal scalars to the temporal scalars*)
	A1=Delete[\[Lambda]AAS//ArrayRules,-1]//ReplaceAll[#,{x_,y_,z_,w_}->{x+ns,y+ns,z+ns,w+ns}]&;
	A2=ReplaceAll[Delete[\[Lambda]KVec//ArrayRules,-1],{x_,y_,z_,w_}->{x+ns,y+ns,z,w}];
	A3=Delete[\[Lambda]4Temp//ArrayRules,-1];
	\[Lambda]4EP=SymmetrizedArray[Join[A1,A2,A3],{nSH,nSH,nSH,nSH},Symmetric[{1,2,3,4}]]//SparseArray;
	
(*We first change name of the tree-level cubic couplings*)
	VarCubic=Join[\[Lambda]3//Normal//Variables]//DeleteDuplicates;
	SubCubic=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarCubic}];
	\[Lambda]3Temp= \[Lambda]3//Normal//ReplaceAll[#,SubCubic]&//SparseArray;
	
	(*We then combine normal and temporal scalars*)
	A1=ReplaceAll[Delete[\[Lambda]vvsLS//ArrayRules,-1],{x_,y_,z_}->{x+ns,y+ns,z}];
	A2=Delete[\[Lambda]3Temp//ArrayRules,-1];
	\[Lambda]3EP=SymmetrizedArray[Join[A1,A2],{nSH,nSH,nSH},Symmetric[{1,2,3}]]//SparseArray;

(*We first change name of the tree-level masses*)
	VarMass=Join[\[Mu]ij//Normal//Variables]//DeleteDuplicates;
	SubMass=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarMass}];
	\[Mu]ijTemp=\[Mu]ij//Normal//ReplaceAll[#,SubMass]&//SparseArray;

	(*We then combine normal and temporal scalar masses*)
	A1=ReplaceAll[Delete[\[Mu]abDef//SparseArray//ArrayRules,-1],{x_,y_}->{x+ns,y+ns}];
	A2=Delete[\[Mu]ijTemp//ArrayRules,-1];

	\[Mu]ijEP=SymmetrizedArray[Join[A1,A2],{nSH,nSH},Symmetric[{1,2}]]//SparseArray;


(*We first change the name of trillinear scalar-vector couplings*)
	VarGauge=Join[gvss//Normal//Variables]//DeleteDuplicates;
	SubGauge=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarGauge}];
	gvssTEMP=gvss//Normal//ReplaceAll[#,SubGauge]&//SparseArray;

	(*we then add vector-temporal trilinear couplings*)
	If[Length[gvvv//ArrayRules]==1,
		gvssEP=SymmetrizedArray[{{1,1,1}->0},{nv,nSH,nSH},Antisymmetric[{2,3}]]//SparseArray;
		gvssEP[[;;,1;;ns,1;;ns]]=gvssTEMP;
	,
		A1=ReplaceAll[Delete[gvvv//ArrayRules,-1],{x_,y_,z_}->{x,y+ns,z+ns}];
		A2=Delete[gvssTEMP//ArrayRules,-1];
		gvssEP=SparseArray[Join[A1,A2],{nv,nSH,nSH}];
	];
	
	nsEP=Length[gvssEP[[1]]];
	nvEP=Length[gvvvEP];
	
	Print["The theory now has ", nsEP, " scalar degrees of freedom"];
	


];


(*
	Defines parameters directly that are used by the effective potential
*)
UseUltraSoftTheory[]:=Module[{},
(*We only need to modify the gauge couplings*)
	SubGauge=Table[c->Symbol[ToString[c]<>ToString["3dUS"]],{c,GaugeCouplingNames}];
	gvvvEP=gvvv//Normal//ReplaceAll[#,SubGauge]&//SparseArray;

(*We only need to modify the gauge couplings*)
	SubGauge=Table[c->Symbol[ToString[c]<>ToString["US"]],{c,Variables[Normal[\[Lambda]4Light]]}];
	\[Lambda]4EP=\[Lambda]4Light//Normal//ReplaceAll[#,SubGauge]&//SparseArray;
(*We only need to modify the gauge couplings*)
	SubGauge=Table[c->Symbol[ToString[c]<>ToString["US"]],{c,Variables[Normal[\[Lambda]3CLight]]}];
	\[Lambda]3EP=\[Lambda]3CLight//Normal//ReplaceAll[#,SubGauge]&//SparseArray;
(*We only need to modify the gauge couplings*)
	SubGauge=Table[c->Symbol[ToString[c]<>ToString["US"]],{c,Variables[Normal[\[Mu]ijLight]]}];
	\[Mu]ijEP=\[Mu]ijLight//Normal//ReplaceAll[#,SubGauge]&//SparseArray;

(*We only need to modify the gauge couplings*)
	VarGauge=Table[Symbol[ToString[c]<>ToString["3d"]],{c,GaugeCouplingNames}];
	SubGauge=Table[c->Symbol[ToString[c]<>ToString["US"]],{c,VarGauge}];
	gvssEP=gvssL//Normal//ReplaceAll[#,SubGauge]&//SparseArray;

	nsEP=Length[gvssEP[[1]]];
	nvEP=Length[gvvvEP];
	
	Print["The theory now has ", nsEP, " scalar degrees of freedom"];

];


(*
	DiagonalMatrix was only introduced in Mathematica 12. So this is a hack for older versions.
	Did I mention that I am a fan of Mathematica? 
*)
DiagonalMatrixQAE[matI_]:=Module[{matP=matI},
	VarMat=#->RandomInteger[10000]&/@(Normal[matP]//Variables);
	matDia=Normal[DiagonalMatrix[Diagonal[matP]]]/.VarMat;
	matT=Normal[matP]/.VarMat;
	Return[matDia==matT]
];


{\[Mu]ij\[Phi],\[Mu]ab\[Phi],Gvvs\[Phi],\[Lambda]3\[Phi],\[Lambda]4\[Phi],gvvv\[Phi],gvss\[Phi],\[Phi]Vev};


(*
	Creates field-dependent couplings. These are created from the user's \[Phi]Vecp input.
*)
DefineVEVS[\[Phi]Vecp_]:=Module[{\[Phi]Vec=\[Phi]Vecp},

	If[Length[\[Phi]Vec]!=nsEP,
		Print["Wrong dimension on the vacuum expectation value"];
		Print["Setting all unspecified VeVs to zero"];
		\[Phi]Vec=SparseArray[\[Phi]Vec//ArrayRules,{nsEP}];
	];
	
	HabijVEP=Transpose[Activate@TensorContract[Inactive@TensorProduct[gvssEP,gvssEP],{{3,5}}],{1,3,2,4}]+Transpose[Transpose[Activate@TensorContract[Inactive@TensorProduct[gvssEP,gvssEP],{{3,5}}],{1,3,2,4}],{2,1,3,4}];
	\[Phi]Vev=\[Phi]Vec//SparseArray;
	\[Mu]ij\[Phi]=\[Mu]ijEP+1/2 Activate@TensorContract[Inactive@TensorProduct[\[Lambda]4EP,\[Phi]Vec,\[Phi]Vec],{{3,5},{4,6}}]+ Activate@TensorContract[Inactive@TensorProduct[\[Lambda]3EP,\[Phi]Vec],{{3,4}}]//SparseArray;
	\[Mu]ab\[Phi]=-1/2 Activate@TensorContract[Inactive@TensorProduct[HabijVEP,\[Phi]Vec,\[Phi]Vec],{{3,5},{4,6}}]//SparseArray;
	Gvvs\[Phi]=-Activate@TensorContract[Inactive@TensorProduct[HabijVEP,\[Phi]Vec],{{4,5}}];
	\[Lambda]3\[Phi]=\[Lambda]3EP+Activate@TensorContract[Inactive@TensorProduct[\[Lambda]4EP,\[Phi]Vec],{{4,5}}];
	\[Lambda]4\[Phi]=\[Lambda]4EP//SparseArray;
	gvvv\[Phi]=gvvvEP//SparseArray;
	gvss\[Phi]=gvssEP//SparseArray;

	If[DiagonalMatrixQAE[\[Mu]ab\[Phi]]==False,
		Print["The Vector Mass-Matrix is not Diagonal"];
	];


	If[DiagonalMatrixQAE[\[Mu]ij\[Phi]]==False,
		Print["The Scalar Mass-Matrix is not Diagonal"];
	];

];


(* ::Section:: *)
(*Calculation of the effective potential*)


	Options[CalculatePotentialUS] = {CustomMasses -> False,PerturbativeDiagonalization->False}
	Options[CalculatePotential] = {CustomMasses -> False,PerturbativeDiagonalization->False}


(*
	Calculates the effective potential with custom masses
*)
CalculatePotentialUS[ScalMassI_,VecMassI_,Opt_]:=Module[{},

	CalculatePotential[ScalMassI,VecMassI,Opt]
];

(*
	Calculates the effective potential with custom masses (just a copy of CalculatePotentialUS)
*)
CalculatePotential[ScalMassI_,VecMassI_,OptionsPattern[]]:=Module[{ScalMassP=ScalMassI,VecMassP=VecMassI},

	CustomMass= OptionValue[CustomMasses];
	If[CustomMass==True,
		\[Mu]ij\[Phi]=ScalMassP//SparseArray;
		\[Mu]ab\[Phi]=VecMassP//SparseArray;

			If[DiagonalMatrixQAE[\[Mu]ij\[Phi]]==True && DiagonalMatrixQAE[\[Mu]ab\[Phi]]==True,
				CalculateLOPotentialSS[];
				CalculateNLOPotentialSS[];
				VTot={VLO,VNLO};
			,
			Print["The Mass matrices are not diagonal. Please rotate to the mass-basis using RotateTensorsUSPostVeV[]"];
			];
	,
	Print["Please set CustomMasses->True"]
	];
];


(*
	Calculates the effective potential.
*)
CalculatePotentialUS[OptionsPattern[]]:=Module[{},

	CustomMass= OptionValue[CustomMasses]; (*Whether the user wants to use their own masses*)
	PertDia=OptionValue[PerturbativeDiagonalization]; (*If a perturbative diagonalization in terms of off-diagonal masses should be carried out*)
	
	If[PertDia==False,
		If[CustomMass==False,
			If[DiagonalMatrixQAE[\[Mu]ij\[Phi]]==True && DiagonalMatrixQAE[\[Mu]ab\[Phi]]==True,
				CalculateLOPotentialSS[];
				CalculateNLOPotentialSS[];
				CalculateNNLOPotentialSS[];

				VTot={VLO,VNLO,VNNLO};
			,
			Print["The Mass matrices are not diagonal. Please rotate to the mass-basis using RotateTensorsUSPostVeV[]"];
			];
		,
			Print["Please supply scalar and vector mass matrices"]
		];
	,
		(*We now treat off-diagonal masses as small*)
		\[Mu]ijPert=ArrayRules[\[Mu]ij\[Phi]]/.({x_Integer,y_Integer}->a_)/;Equal[x,y]->{x,y}->0//SparseArray;
		\[Mu]ij\[Phi]=ArrayRules[\[Mu]ij\[Phi]]/.({x_Integer,y_Integer}->a_)/;Unequal[x,y]->{x,y}->0//SparseArray;
		
		\[Mu]ab\[Phi]Pert=ArrayRules[\[Mu]ab\[Phi]]/.({x_Integer,y_Integer}->a_)/;Equal[x,y]->{x,y}->0//SparseArray;
		\[Mu]ab\[Phi]=ArrayRules[\[Mu]ab\[Phi]]/.({x_Integer,y_Integer}->a_)/;Unequal[x,y]->{x,y}->0//SparseArray;
		
		CalculateLOPotentialSS[];
		CalculateNLOPotentialSS[];
		CalculateNNLOPotentialSS[];
		CalculateOffDiagPotentialSS[];

		VTot={VLO,VNLO,VNNLO};
	
	];


];

(*
	Calculates the effective potential. (just a copy of CalculatePotentialUS)
*)
CalculatePotential[OptionsPattern[]]:=Module[{},

	CustomMass= OptionValue[CustomMasses]; (*Whether the user wants to use their own masses*)
	PertDia=OptionValue[PerturbativeDiagonalization]; (*If a perturbative diagonalization in terms of off-diagonal masses should be carried out*)
	
	If[PertDia==False,
		If[CustomMass==False,
			If[DiagonalMatrixQAE[\[Mu]ij\[Phi]]==True && DiagonalMatrixQAE[\[Mu]ab\[Phi]]==True,
				CalculateLOPotentialSS[];
				CalculateNLOPotentialSS[];
				CalculateNNLOPotentialSS[];

				VTot={VLO,VNLO,VNNLO};
			,
			Print["The Mass matrices are not diagonal. Please rotate to the mass-basis using RotateTensorsUSPostVeV[]"];
			];
		,
			Print["Please supply scalar and vector mass matrices"]
		];
	,
		(*We now treat off-diagonal masses as small*)
		\[Mu]ijPert=ArrayRules[\[Mu]ij\[Phi]]/.({x_Integer,y_Integer}->a_)/;Equal[x,y]->{x,y}->0//SparseArray;
		\[Mu]ij\[Phi]=ArrayRules[\[Mu]ij\[Phi]]/.({x_Integer,y_Integer}->a_)/;Unequal[x,y]->{x,y}->0//SparseArray;
		
		\[Mu]ab\[Phi]Pert=ArrayRules[\[Mu]ab\[Phi]]/.({x_Integer,y_Integer}->a_)/;Equal[x,y]->{x,y}->0//SparseArray;
		\[Mu]ab\[Phi]=ArrayRules[\[Mu]ab\[Phi]]/.({x_Integer,y_Integer}->a_)/;Unequal[x,y]->{x,y}->0//SparseArray;
		
		CalculateLOPotentialSS[];
		CalculateNLOPotentialSS[];
		CalculateNNLOPotentialSS[];
		CalculateOffDiagPotentialSS[];

		VTot={VLO,VNLO,VNNLO};
	
	];


];


(*
	Prints the effective potential.
*)
PrintEffectivePotential[]:=Module[{},
	EffPotPrint=VTot[[1]]+VTot[[2]]+VTot[[3]];

(*Printing Result*)
	OutputFormatDR[EffPotPrint]
];

PrintEffectivePotential[optP_]:=Module[{opt=optP},
	EffPotPrint=Switch[opt,"LO",VTot[[1]],"NLO",VTot[[2]],"NNLO",VTot[[3]],__,VTot[[1]]+VTot[[2]]+VTot[[3]]];

(*Printing Result*)
	OutputFormatDR[EffPotPrint]
];


{VLO,VNLO,VNNLO};


(*
	Calculates the tree-level effective potential.
*)
CalculateOffDiagPotentialSS[]:=Module[{V1,V2,helpTens,AD},
		
(*Scalar integrals*)	
	AD[x_,y_]:=-(1/(4 \[Pi]))( Sqrt[x]-Sqrt[y])/(y-x);
	AD[0,0]:=0;
	AD[x_,x_]:=(1/(8 \[Pi]))/Sqrt[x];
	
(*Scalar contribution*)
		helpTens=Table[AD[a,b],{a,Diagonal[\[Mu]ij\[Phi]]},{b,Diagonal[\[Mu]ij\[Phi]]}]//SparseArray;
		helpTens=TensorProduct[helpTens,\[Mu]ijPert]//SparseArray//DiagonalTensor2[#,1,3]&;
		V1=Tr[helpTens . \[Mu]ijPert];

(*Vector contribution*)
		helpTens=Table[AD[a,b],{a,Diagonal[\[Mu]ab\[Phi]]},{b,Diagonal[\[Mu]ab\[Phi]]}]//SparseArray;
		helpTens=TensorProduct[helpTens,\[Mu]ab\[Phi]Pert]//SparseArray//DiagonalTensor2[#,1,3]&;
		V2=2*Tr[helpTens . \[Mu]ab\[Phi]Pert];
		
		
		VNNLO+=V1+V2;
		
	];


(*
	Calculates the two-loop effective potential.
*)
CalculateNNLOPotentialSS[]:=Module[{Vss, Vsss, Vvvs, Vssv, Vvs, Vvvv, Vvv, V\[Eta]\[Eta]v
									,fvvv,f\[Eta]\[Eta]v,fssv,fvvs,fsss,fss,fvs,fvv,aS,av,
									ss,sss,vvs,ssv,vs,vv,vvv,ggv},
If[verbose==True,Print["Calculating the 2-Loop Effective Potential"]];

	
(*Contraction with coupling-tensors and masses*)

	(*Assuming that the mass matrices are diagonal*)
	aS=Table[\[Mu]ij\[Phi][[i,i]],{i,1,nsEP}]//SparseArray;
	av=Table[\[Mu]ab\[Phi][[i,i]],{i,1,nvEP}]//SparseArray;

(*Loading two-loop master integrals*)
	{fvvv,f\[Eta]\[Eta]v,fssv,fvvs,fsss,fss,fvs,fvv}=TwoLoopFunctions[];
(*Potential*)

	If[Length[TensorDimensions[\[Lambda]4\[Phi]]]<4,
	(*This occurs when the FastRotation->True option is used for mass diagonalization*)
		ss=1/8 TensorProduct[\[Lambda]4\[Phi]];
		Vss=Sum[ss[[j,k]]fss[aS[[j]],aS[[k]]],{j,nsEP},{k,nsEP}];
	,
		ss=1/8 TensorProduct[\[Lambda]4\[Phi]];
		Vss=Sum[ss[[j,j,k,k]]fss[aS[[j]],aS[[k]]],{j,nsEP},{k,nsEP}];
	];
	sss=1/12 TensorProduct[\[Lambda]3\[Phi],\[Lambda]3\[Phi]];
	Vsss=Sum[sss[[i,j,k,i,j,k]]fsss[aS[[i]],aS[[j]],aS[[k]]],{j,nsEP},{k,nsEP},{i,nsEP}];
	vvs=1/4 TensorProduct[Gvvs\[Phi],Gvvs\[Phi]];
	Vvvs=Sum[vvs[[a,b,i,a,b,i]]fvvs[av[[a]],av[[b]],aS[[i]]],{a,nvEP},{b,nvEP},{i,nsEP}];
	ssv=1/4 TensorProduct[gvss\[Phi],gvss\[Phi]];
	Vssv=Sum[ssv[[a,i,j,a,i,j]]fssv[aS[[i]],aS[[j]],av[[a]]],{a,nvEP},{j,nsEP},{i,nsEP}];
	vs=1/2 TensorProduct[gvss\[Phi],gvss\[Phi]];
	Vvs=Sum[vs[[a,i,j,a,i,j]]fvs[aS[[i]],av[[a]]],{a,nvEP},{j,nsEP},{i,nsEP}];
	vv=1/4 TensorProduct[gvvv\[Phi],gvvv\[Phi]];
	Vvv=Sum[vv[[a,b,c,a,b,c]]fvv[av[[a]],av[[b]]],{a,nvEP},{b,nvEP},{c,nvEP}];
	vvv=1/12 TensorProduct[gvvv\[Phi],gvvv\[Phi]];
	Vvvv=Sum[vvv[[a,b,c,a,b,c]]fvvv[av[[a]],av[[b]],av[[c]]],{a,nvEP},{b,nvEP},{c,nvEP}];
	ggv=1/4 TensorProduct[gvvv\[Phi],gvvv\[Phi]];
	V\[Eta]\[Eta]v=Sum[ggv[[a,b,c,a,b,c]]f\[Eta]\[Eta]v[0,0,av[[b]]],{a,nvEP},{b,nvEP},{c,nvEP}];

	VNNLO= Vss+ Vsss+Vvvs+ Vssv+ Vvs+ Vvvv+ Vvv+ V\[Eta]\[Eta]v;

];


(*
	Calculates the one-loop effective potential.
*)
CalculateNLOPotentialSS[]:=Module[{ALog},
If[verbose==True,Print["Calculating the 1-Loop Effective Potential"]];
	
	ALog[x_]:=-(x^(3/2)/(12 \[Pi]));
		
	V1=Sum[ALog[\[Mu]ij\[Phi][[i,i]]],{i,1,nsEP}];
	V2=2*Sum[ALog[\[Mu]ab\[Phi][[i,i]]],{i,1,nvEP}];

	VNLO=V1+V2;
];


(*
	Calculates the tree-level effective potential.
*)
CalculateLOPotentialSS[]:=Module[{},
	If[verbose==True,Print["Calculating the Tree-Level Effective Potential"]];

	V1=\[Lambda]4EP . \[Phi]Vev . \[Phi]Vev . \[Phi]Vev . \[Phi]Vev;
	V2=\[Mu]ijEP . \[Phi]Vev . \[Phi]Vev;
	V3=\[Lambda]3EP . \[Phi]Vev . \[Phi]Vev . \[Phi]Vev;
	VLO=1/4! V1+V2/2!+V3/3!;
	];


(* ::Section::Closed:: *)
(*Rotations and diagonalization*)


DiagonalTensor2[s_SparseArray,a_Integer,b_Integer] := With[
    {
    s1=Flatten[s,{{a},{b}}]
    },
     Table[i,{i,s1}]//Table[#[[i,i]],{i,1,Length[#]}]&//SparseArray
    ]


errPrintTensVEV="(1) Scalar Mass, (2) Vector Mass, (3) Scalar Quartic, (4) Scalar Cubic, (5) VSS, (6) VVS, (7) VVV";


PrintTensorsVEV[]:="nothing"/;Print[errPrintTensVEV];


PrintTensorsVEV[ind_]:=Module[{},	
	
	Switch[ind, 
    1, Return[OutputFormatDR[\[Mu]ij\[Phi]]], 
    2, Return[OutputFormatDR[\[Mu]ab\[Phi]]], 
    3, Return[OutputFormatDR[\[Lambda]4\[Phi]]], 
    4, Return[OutputFormatDR[\[Lambda]3\[Phi]]], 
    5, Return[OutputFormatDR[gvss\[Phi]]], 
    6, Return[OutputFormatDR[Gvvs\[Phi]]], 
    7, Return[OutputFormatDR[gvvv\[Phi]]],   
    _, 
    Print[errPrintTensVEV];
  ];
      Return[""]

];


(*
	Rotates to a diagonal-mass basis.
*)
RotateTensorsUSPostVEV[DScalarsp_,DVectorsp_]:=Module[{DS=DScalarsp,DV=DVectorsp},
	DS=DS//SparseArray;
	DV=DV//SparseArray;

	\[Lambda]4\[Phi]=Transpose[DS] . \[Lambda]4\[Phi] . DS//Activate@TensorContract[Inactive@TensorProduct[DS,#],{{1,4}}]&//Transpose[#,{2,1,3,4}]&//Activate@TensorContract[Inactive@TensorProduct[DS,#],{{1,5}}]&//Transpose[#,{3,2,1,4}]&//SimplifySparse;
	\[Lambda]3\[Phi]=Transpose[DS] . \[Lambda]3\[Phi] . DS//Activate@TensorContract[Inactive@TensorProduct[DS,#],{{1,4}}]&//Transpose[#,{2,1,3}]&//SimplifySparse;
	gvss\[Phi]=Transpose[DV] . gvss\[Phi] . DS//Activate@TensorContract[Inactive@TensorProduct[DS,#],{{1,4}}]&//Transpose[#,{2,1,3}]&//SimplifySparse;
	Gvvs\[Phi]=Transpose[DV] . Gvvs\[Phi] . DS//Activate@TensorContract[Inactive@TensorProduct[DV,#],{{1,4}}]&//Transpose[#,{2,1,3}]&//SimplifySparse;
	gvvv\[Phi]=Transpose[DV] . gvvv\[Phi] . DV//Activate@TensorContract[Inactive@TensorProduct[DV,#],{{1,4}}]&//Transpose[#,{2,1,3}]&//SimplifySparse;
	\[Mu]ab\[Phi]=Activate@TensorContract[Inactive@TensorProduct[DV,DV,\[Mu]ab\[Phi]],{{1,5},{3,6}}]//SimplifySparse;
	\[Mu]ij\[Phi]=Activate@TensorContract[Inactive@TensorProduct[DS,DS,\[Mu]ij\[Phi]],{{1,5},{3,6}}]//SimplifySparse;


	If[DiagonalMatrixQAE[\[Mu]ab\[Phi]]==False,
		Print["The Vector mass-Matrix is not diagonal"];
	];


	If[DiagonalMatrixQAE[\[Mu]ij\[Phi]]==False,
		Print["The Scalar mass-Matrix is not diagonal"];
	];

];


	Options[RotateTensorsCustomMass] = {FastRotation->False}


(*
	Rotates to a diagonal-mass basis.
*)
RotateTensorsCustomMass[DScalarsp_,DVectorsp_,ScalarMass_,vectorMass_,OptionsPattern[]]:=Module[{DS=DScalarsp,DV=DVectorsp},
	
	DS=DS//SparseArray;
	DV=DV//SparseArray;
	
	If[OptionValue[FastRotation],
		(*For large rotation matrices we can simplify the two-loop calculation by only doing the rotation for Subscript[\[Lambda], iikk] components*)
		RotHelp=Flatten[TensorProduct[DS,DS],{{2},{4},{1,3}}]//DiagonalTensor2[#,1,2]&;
		\[Lambda]4\[Phi]=RotHelp . Flatten[\[Lambda]4\[Phi] ,{{1,2},{3,4}}] . Transpose[RotHelp];
		
		RotHelp=Flatten[TensorProduct[DS,DS,DS],{{2},{4},{6},{1,3,5}}];
		\[Lambda]3\[Phi]=RotHelp . Flatten[\[Lambda]3\[Phi]];

		RotHelp=Flatten[TensorProduct[DV,DS,DS],{{2},{4},{6},{1,3,5}}];
		gvss\[Phi]=RotHelp . Flatten[gvss\[Phi]];
		
		RotHelp=Flatten[TensorProduct[DV,DV,DS],{{2},{4},{6},{1,3,5}}];
		Gvvs\[Phi]=RotHelp . Flatten[Gvvs\[Phi]];
		
		RotHelp=Flatten[TensorProduct[DV,DV,DV],{{2},{4},{6},{1,3,5}}];
		gvvv\[Phi]=RotHelp . Flatten[gvvv\[Phi]];
	,
		\[Lambda]4\[Phi]=Transpose[DS] . \[Lambda]4\[Phi] . DS//Activate@TensorContract[Inactive@TensorProduct[DS,#],{{1,4}}]&//Transpose[#,{2,1,3,4}]&//Activate@TensorContract[Inactive@TensorProduct[DS,#],{{1,5}}]&//Transpose[#,{3,2,1,4}]&//SimplifySparse;
		\[Lambda]3\[Phi]=Transpose[DS] . \[Lambda]3\[Phi] . DS//Activate@TensorContract[Inactive@TensorProduct[DS,#],{{1,4}}]&//Transpose[#,{2,1,3}]&//SimplifySparse;
		gvss\[Phi]=Transpose[DV] . gvss\[Phi] . DS//Activate@TensorContract[Inactive@TensorProduct[DS,#],{{1,4}}]&//Transpose[#,{2,1,3}]&//SimplifySparse;
		Gvvs\[Phi]=Transpose[DV] . Gvvs\[Phi] . DS//Activate@TensorContract[Inactive@TensorProduct[DV,#],{{1,4}}]&//Transpose[#,{2,1,3}]&//SimplifySparse;
		gvvv\[Phi]=Transpose[DV] . gvvv\[Phi] . DV//Activate@TensorContract[Inactive@TensorProduct[DV,#],{{1,4}}]&//Transpose[#,{2,1,3}]&//SimplifySparse;
	];
	


	\[Mu]ab\[Phi]=vectorMass//SparseArray;
	\[Mu]ij\[Phi]=ScalarMass//SparseArray;


	If[DiagonalMatrixQAE[\[Mu]ab\[Phi]]==False,
		Print["The Vector mass-Matrix is not diagonal"];
	];


	If[DiagonalMatrixQAE[\[Mu]ij\[Phi]]==False,
		Print["The Scalar mass-Matrix is not diagonal"];
	];

];


(* ::Section:: *)
(*Scalar master integrals*)


TwoLoopFunctions[]:=Module[{f,A,I2,I2Div,AD,ALog,fvvv,f\[Eta]\[Eta]v,fssv,fvvs,fsss,fss,fvs,fvv},
(*The notation follows Martin's convention*)
(*Please see arXiv:1808.07615*)
	f[x_]:=-(1/(12 \[Pi])) x^(3/2);
	A[0]=0;
	A[x_]:=-(1/(4 \[Pi])) Sqrt[x];
	I2[0,0,0]=0;
	I2[x_,y_,z_]:=1/(4 \[Pi])^2 (Log[\[Mu]3US/(Sqrt[x]+Sqrt[y]+Sqrt[z])]+1/2);
	I2Div=1/(4 \[Pi])^2 1/4;

(*Scalar integrals*)	
	AD[x_,y_]:=-(1/(4 \[Pi]))( Sqrt[x]-Sqrt[y])/(y-x);
	AD[0,0]:=0;
	AD[x_,x_]:=(1/(8 \[Pi]))/Sqrt[x];
	ALog[x_]:=-(x^(3/2)/(12 \[Pi]));
		
		
(*Sunset with three vector bosons*)
	fvvv[x_,y_,z_]:=1/(4 x y z) (1/3 (x^4 (2 I2Div-3 I2[0,0,x])+(x^2-2 x (y+z)+(y-z)^2) (-(3 (x^2+6 x (y+z)+y^2+6 y z+z^2) I2[x,y,z])
							-I2Div (3 (-(8 x (y+z))-8 y z)-2 (x^2+6 x (y+z)+y^2+6 y z+z^2)))+(x-y)^2 (3 (x^2+6 x y+y^2) I2[0,x,y]
							+I2Div (-(2 (x^2+6 x y+y^2))-24 x y))+(x-z)^2 (3 (x^2+6 x z+z^2) I2[0,x,z]+I2Div (-(2 (x^2+6 x z+z^2))-24 x z))
							+y^4 (2 I2Div-3 I2[0,0,y])+(y-z)^2 (3 (y^2+6 y z+z^2) I2[0,y,z]+I2Div (-(2 (y^2+6 y z+z^2))-24 y z))+z^4 (2 I2Div-3 I2[0,0,z])
							+x A[y] A[z] (3 x^2+15 x y+15 x z-15 y^2-26 y z-15 z^2)+z A[x] A[y] (-(15 x^2)-26 x y+15 x z-15 y^2+15 y z+3 z^2)
							+y A[x] A[z] (-(15 x^2)+15 x y-26 x z+3 y^2+15 y z-15 z^2))+40/3 (I2Div x^2 y z+I2Div x y^2 z+I2Div x y z^2));
	(*Special cases where some massesa are 0*)	
		fvvv[0,0,0]:=0;
	
		fvvv[0,y_,z_]:=1/(6 y z) (-(6 y^3 I2[0,y,z])+6 y^3 I2[0,0,y]+30 y^2 z I2[0,y,z]-6 z^3 I2[0,y,z]+30 y z^2 I2[0,y,z]
					+6 z^3 I2[0,0,z]-6 y^2 A[y] A[z]-6 z^2 A[y] A[z]-4 y z A[y] A[z]-57 I2Div y^2 z-57 I2Div y z^2);
		fvvv[x_,0,z_]:=fvvv[0,x,z];
		fvvv[x_,y_,0]:=fvvv[0,x,y];
		fvvv[0,0,z_]:=1/4 z (20 I2[0,0,z]-46 I2Div);
		fvvv[0,y_,0]:=fvvv[0,0,y];
		fvvv[x_,0,0]:=fvvv[0,0,x];

(*Ghost-vector sunset*)

		f\[Eta]\[Eta]v[x_,y_,z_]:=((-x+y+z) A[y] A[z]+A[x] (-z A[y]+(x-y+z) A[z])-(x-y)^2 I2[0,x,y]+(y^2+(x-z)^2-2 y (x+z)) I2[x,z,y])/(2 z);
	(*Special cases where vector-mass is zero*)
		f\[Eta]\[Eta]v[x_,y_,0]:=1/2 (2 I2Div (x+y)-2 (A[x] A[y]+(x+y) I2[0,x,y]));
		f\[Eta]\[Eta]v[0,0,0]:=0;
		f\[Eta]\[Eta]v[0,0,z]:=1/2 z I2[0,0,z];

(*Scalar-scalar-vector sunset*)
	fssv[x_,y_,z_]:=1/z (-((-x+y+z) A[y] A[z])-A[x] (-z A[y]+(x-y+z) A[z])+(x-y)^2 I2[0,x,y]-((x-y)^2-2 (x+y) z+z^2) I2[x,y,z]);
	(*Special case where vector-mass is zero*)
		fssv[x_,y_,0]:=-2 I2Div (x+y)-2 (-A[x] A[y]-(x+y) I2[0,x,y]);
		fssv[0,0,0]:=0;

(*Vector-vector-scalar sunset*)
	fvvs[x_,y_,z_]:=1/(4 x y) (8 I2Div x y-(x+y-z) A[x] A[y]+y A[x] A[z]+x A[y] A[z]
						-z^2 I2[0,0,z]+(x-z)^2 I2[0,x,z]+(y-z)^2 I2[0,y,z]-(x^2+6 x y+(y-z)^2-2 x z) I2[x,y,z]);
	(*Special case where vector-masses are zero*)
		fvvs[0,y_,z_]:=(6 I2Div y+2 (A[y] A[z]-z I2[0,0,z]+(-3 y+z) I2[0,y,z]))/(4 y);
		fvvs[x_,0,z_]:=fvvs[0,x,z];
		fvvs[0,0,z_]:=1/4 (10 I2Div-6 I2[0,0,z]);
		fvvs[0,0,0]:=0;
		
(*scalar-scalar-scalar sunset*)		
	fsss[x_,y_,z_]:=-I2[x,y,z];
	
(*scalar-scalar bubble*)
	fss[x_,y_]:=A[x]A[y];
	
(*scalar-vector bubble*)
	fvs[x_,y_]:=2 A[x]A[y];
	
(*vector-vector bubble*)	
	fvv[x_,y_]:=8/3 A[x]A[y];
	
	Return[{fvvv,f\[Eta]\[Eta]v,fssv,fvvs,fsss,fss,fvs,fvv}];
]

