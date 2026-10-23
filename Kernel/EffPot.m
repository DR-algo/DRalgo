(* ::Package:: *)

(* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ *)

(* :Title: EffPot                                                          	*)

(*
      	This software is covered by the GNU General Public License 3.
       	Copyright (C) 2021-2026 Andreas Ekstedt
       	Copyright (C) 2021-2026 Philipp Schicho
       	Copyright (C) 2021-2026 Tuomas V.I. Tenkanen

*)

(* :Summary:	Compute the effective potential.        					*)	

(* ------------------------------------------------------------------------ *)


(* ::Section::Closed:: *)
(*Definition of model*)


(*
	User-defined model.
*)
DefineNewTensorsUS[muijI_,lam4I_,lam3I_,gvssI_,gvvvI_]:=Module[{muijP=muijI,lam4P=lam4I,lam3P=lam3I,gvssP=gvssI,gvvvP=gvvvI},
	muijEP=muijP//SparseArray;
	gvvvEP=gvvvP//SparseArray;
	gvssEP=gvssP//SparseArray;
	lam4EP=lam4P//SparseArray;
	lam3EP=lam3P//SparseArray;
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
	VarQuartic=Join[lam4//Normal//Variables]//DeleteDuplicates;
	SubQuartic=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarQuartic}];
	lam4Temp= lam4//Normal//ReplaceAll[#,SubQuartic]&//SparseArray;
	
	(*We then combine normal and temporal scalars*)
	ns=(lam4//Length);
	nSH=nv+ns; (*Adds normal scalars to the temporal scalars*)
	A1=Delete[lamAAS//ArrayRules,-1]//ReplaceAll[#,{x_,y_,z_,w_}->{x+ns,y+ns,z+ns,w+ns}]&;
	A2=ReplaceAll[Delete[lamKVec//ArrayRules,-1],{x_,y_,z_,w_}->{x+ns,y+ns,z,w}];
	A3=Delete[lam4Temp//ArrayRules,-1];
	lam4EP=SymmetrizedArray[Join[A1,A2,A3],{nSH,nSH,nSH,nSH},Symmetric[{1,2,3,4}]]//SparseArray;
	
(*We first change name of the tree-level cubic couplings*)
	VarCubic=Join[lam3//Normal//Variables]//DeleteDuplicates;
	SubCubic=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarCubic}];
	lam3Temp= lam3//Normal//ReplaceAll[#,SubCubic]&//SparseArray;
	
	(*We then combine normal and temporal scalars*)
	A1=ReplaceAll[Delete[lamvvsLS//ArrayRules,-1],{x_,y_,z_}->{x+ns,y+ns,z}];
	A2=Delete[lam3Temp//ArrayRules,-1];
	lam3EP=SymmetrizedArray[Join[A1,A2],{nSH,nSH,nSH},Symmetric[{1,2,3}]]//SparseArray;

(*We first change name of the tree-level masses*)
	VarMass=Join[muij//Normal//Variables]//DeleteDuplicates;
	SubMass=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarMass}];
	muijTemp=muij//Normal//ReplaceAll[#,SubMass]&//SparseArray;

	(*We then combine normal and temporal scalar masses*)
	A1=ReplaceAll[Delete[Normal[muabDef]/.SubDebyeDecoupled//SparseArray//ArrayRules,-1],{x_,y_}->{x+ns,y+ns}];
	A2=Delete[muijTemp//ArrayRules,-1];

	muijEP=SymmetrizedArray[Join[A1,A2],{nSH,nSH},Symmetric[{1,2}]]//SparseArray;


(*We first change the name of trillinear scalar-vector couplings*)
	VarGauge=Join[gvss//Normal//Variables]//DeleteDuplicates;
	SubGauge=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarGauge}];
	gvssTEMP=gvss//Normal//ReplaceAll[#,SubGauge]&//SparseArray;

	(*we then add vector-temporal trilinear couplings, with the 3d gauge couplings of gvvvEP*)
	If[Length[gvvv//ArrayRules]==1,
		gvssEP=SymmetrizedArray[{{1,1,1}->0},{nv,nSH,nSH},Antisymmetric[{2,3}]]//SparseArray;
		gvssEP[[;;,1;;ns,1;;ns]]=gvssTEMP;
	,
		A1=ReplaceAll[Delete[gvvvEP//ArrayRules,-1],{x_,y_,z_}->{x,y+ns,z+ns}];
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
	SubGauge=Table[c->Symbol[ToString[c]<>ToString["US"]],{c,Variables[Normal[lam4Light]]}];
	lam4EP=lam4Light//Normal//ReplaceAll[#,SubGauge]&//SparseArray;
(*We only need to modify the gauge couplings*)
	SubGauge=Table[c->Symbol[ToString[c]<>ToString["US"]],{c,Variables[Normal[lam3CLight]]}];
	lam3EP=lam3CLight//Normal//ReplaceAll[#,SubGauge]&//SparseArray;
(*We only need to modify the gauge couplings*)
	SubGauge=Table[c->Symbol[ToString[c]<>ToString["US"]],{c,Variables[Normal[muijLight]]}];
	muijEP=muijLight//Normal//ReplaceAll[#,SubGauge]&//SparseArray;

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


{muij\[Phi],muab\[Phi],Gvvs\[Phi],lam3\[Phi],lam4\[Phi],gvvv\[Phi],gvss\[Phi],\[Phi]Vev};


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
	muij\[Phi]=muijEP+1/2 Activate@TensorContract[Inactive@TensorProduct[lam4EP,\[Phi]Vec,\[Phi]Vec],{{3,5},{4,6}}]+ Activate@TensorContract[Inactive@TensorProduct[lam3EP,\[Phi]Vec],{{3,4}}]//SparseArray;
	muab\[Phi]=-1/2 Activate@TensorContract[Inactive@TensorProduct[HabijVEP,\[Phi]Vec,\[Phi]Vec],{{3,5},{4,6}}]//SparseArray;
	Gvvs\[Phi]=-Activate@TensorContract[Inactive@TensorProduct[HabijVEP,\[Phi]Vec],{{4,5}}];
	lam3\[Phi]=lam3EP+Activate@TensorContract[Inactive@TensorProduct[lam4EP,\[Phi]Vec],{{4,5}}];
	lam4\[Phi]=lam4EP//SparseArray;
	gvvv\[Phi]=gvvvEP//SparseArray;
	gvss\[Phi]=gvssEP//SparseArray;

	If[DiagonalMatrixQAE[muab\[Phi]]==False,
		Print["The Vector Mass-Matrix is not Diagonal"];
	];


	If[DiagonalMatrixQAE[muij\[Phi]]==False,
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
		muij\[Phi]=ScalMassP//SparseArray;
		muab\[Phi]=VecMassP//SparseArray;

			If[DiagonalMatrixQAE[muij\[Phi]]==True && DiagonalMatrixQAE[muab\[Phi]]==True,
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
			If[DiagonalMatrixQAE[muij\[Phi]]==True && DiagonalMatrixQAE[muab\[Phi]]==True,
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
		muijPert=ArrayRules[muij\[Phi]]/.({x_Integer,y_Integer}->a_)/;Equal[x,y]->{x,y}->0//SparseArray;
		muij\[Phi]=ArrayRules[muij\[Phi]]/.({x_Integer,y_Integer}->a_)/;Unequal[x,y]->{x,y}->0//SparseArray;
		
		muab\[Phi]Pert=ArrayRules[muab\[Phi]]/.({x_Integer,y_Integer}->a_)/;Equal[x,y]->{x,y}->0//SparseArray;
		muab\[Phi]=ArrayRules[muab\[Phi]]/.({x_Integer,y_Integer}->a_)/;Unequal[x,y]->{x,y}->0//SparseArray;
		
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
			If[DiagonalMatrixQAE[muij\[Phi]]==True && DiagonalMatrixQAE[muab\[Phi]]==True,
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
		muijPert=ArrayRules[muij\[Phi]]/.({x_Integer,y_Integer}->a_)/;Equal[x,y]->{x,y}->0//SparseArray;
		muij\[Phi]=ArrayRules[muij\[Phi]]/.({x_Integer,y_Integer}->a_)/;Unequal[x,y]->{x,y}->0//SparseArray;
		
		muab\[Phi]Pert=ArrayRules[muab\[Phi]]/.({x_Integer,y_Integer}->a_)/;Equal[x,y]->{x,y}->0//SparseArray;
		muab\[Phi]=ArrayRules[muab\[Phi]]/.({x_Integer,y_Integer}->a_)/;Unequal[x,y]->{x,y}->0//SparseArray;
		
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
CalculateOffDiagPotentialSS[]:=Module[{V1,V2,helpTens},
		
(*Scalar contribution*)
		helpTens=Table[BubbleIntegral3d[a,b],{a,Diagonal[muij\[Phi]]},{b,Diagonal[muij\[Phi]]}]//SparseArray;
		helpTens=TensorProduct[helpTens,muijPert]//SparseArray//DiagonalTensor2[#,1,3]&;
		V1=Tr[helpTens . muijPert];

(*Vector contribution*)
		helpTens=Table[BubbleIntegral3d[a,b],{a,Diagonal[muab\[Phi]]},{b,Diagonal[muab\[Phi]]}]//SparseArray;
		helpTens=TensorProduct[helpTens,muab\[Phi]Pert]//SparseArray//DiagonalTensor2[#,1,3]&;
		V2=2*Tr[helpTens . muab\[Phi]Pert];
		
		
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
	aS=Table[muij\[Phi][[i,i]],{i,1,nsEP}]//SparseArray;
	av=Table[muab\[Phi][[i,i]],{i,1,nvEP}]//SparseArray;

(*Loading two-loop master integrals*)
	{fvvv,f\[Eta]\[Eta]v,fssv,fvvs,fsss,fss,fvs,fvv}=TwoLoopFunctions[mu3US];
(*Potential*)

	If[Length[TensorDimensions[lam4\[Phi]]]<4,
	(*This occurs when the FastRotation->True option is used for mass diagonalization*)
		ss=1/8 TensorProduct[lam4\[Phi]];
		Vss=Sum[ss[[j,k]]fss[aS[[j]],aS[[k]]],{j,nsEP},{k,nsEP}];
	,
		ss=1/8 TensorProduct[lam4\[Phi]];
		Vss=Sum[ss[[j,j,k,k]]fss[aS[[j]],aS[[k]]],{j,nsEP},{k,nsEP}];
	];
	sss=1/12 TensorProduct[lam3\[Phi],lam3\[Phi]];
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
CalculateNLOPotentialSS[]:=Module[{},
If[verbose==True,Print["Calculating the 1-Loop Effective Potential"]];
		
	V1=Sum[OneLoopPotential3d[muij\[Phi][[i,i]]],{i,1,nsEP}];
	V2=2*Sum[OneLoopPotential3d[muab\[Phi][[i,i]]],{i,1,nvEP}];

	VNLO=V1+V2;
];


(*
	Calculates the tree-level effective potential.
*)
CalculateLOPotentialSS[]:=Module[{},
	If[verbose==True,Print["Calculating the Tree-Level Effective Potential"]];

	V1=lam4EP . \[Phi]Vev . \[Phi]Vev . \[Phi]Vev . \[Phi]Vev;
	V2=muijEP . \[Phi]Vev . \[Phi]Vev;
	V3=lam3EP . \[Phi]Vev . \[Phi]Vev . \[Phi]Vev;
	VLO=1/4! V1+V2/2!+V3/3!;
	];


(* ::Section:: *)
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
    1, Return[OutputFormatDR[muij\[Phi]]], 
    2, Return[OutputFormatDR[muab\[Phi]]], 
    3, Return[OutputFormatDR[lam4\[Phi]]], 
    4, Return[OutputFormatDR[lam3\[Phi]]], 
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

	lam4\[Phi]=Transpose[DS] . lam4\[Phi] . DS//Activate@TensorContract[Inactive@TensorProduct[DS,#],{{1,4}}]&//Transpose[#,{2,1,3,4}]&//Activate@TensorContract[Inactive@TensorProduct[DS,#],{{1,5}}]&//Transpose[#,{3,2,1,4}]&//SimplifySparse;
	lam3\[Phi]=Transpose[DS] . lam3\[Phi] . DS//Activate@TensorContract[Inactive@TensorProduct[DS,#],{{1,4}}]&//Transpose[#,{2,1,3}]&//SimplifySparse;
	gvss\[Phi]=Transpose[DV] . gvss\[Phi] . DS//Activate@TensorContract[Inactive@TensorProduct[DS,#],{{1,4}}]&//Transpose[#,{2,1,3}]&//SimplifySparse;
	Gvvs\[Phi]=Transpose[DV] . Gvvs\[Phi] . DS//Activate@TensorContract[Inactive@TensorProduct[DV,#],{{1,4}}]&//Transpose[#,{2,1,3}]&//SimplifySparse;
	gvvv\[Phi]=Transpose[DV] . gvvv\[Phi] . DV//Activate@TensorContract[Inactive@TensorProduct[DV,#],{{1,4}}]&//Transpose[#,{2,1,3}]&//SimplifySparse;
	muab\[Phi]=Activate@TensorContract[Inactive@TensorProduct[DV,DV,muab\[Phi]],{{1,5},{3,6}}]//SimplifySparse;
	muij\[Phi]=Activate@TensorContract[Inactive@TensorProduct[DS,DS,muij\[Phi]],{{1,5},{3,6}}]//SimplifySparse;


	If[DiagonalMatrixQAE[muab\[Phi]]==False,
		Print["The Vector mass-Matrix is not diagonal"];
	];


	If[DiagonalMatrixQAE[muij\[Phi]]==False,
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
		(*For large rotation matrices we can simplify the two-loop calculation by only doing the rotation for Subscript[lam, iikk] components*)
		RotHelp=Flatten[TensorProduct[DS,DS],{{2},{4},{1,3}}]//DiagonalTensor2[#,1,2]&;
		lam4\[Phi]=RotHelp . Flatten[lam4\[Phi] ,{{1,2},{3,4}}] . Transpose[RotHelp];
		
		RotHelp=Flatten[TensorProduct[DS,DS,DS],{{2},{4},{6},{1,3,5}}];
		lam3\[Phi]=RotHelp . Flatten[lam3\[Phi]];

		RotHelp=Flatten[TensorProduct[DV,DS,DS],{{2},{4},{6},{1,3,5}}];
		gvss\[Phi]=RotHelp . Flatten[gvss\[Phi]];
		
		RotHelp=Flatten[TensorProduct[DV,DV,DS],{{2},{4},{6},{1,3,5}}];
		Gvvs\[Phi]=RotHelp . Flatten[Gvvs\[Phi]];
		
		RotHelp=Flatten[TensorProduct[DV,DV,DV],{{2},{4},{6},{1,3,5}}];
		gvvv\[Phi]=RotHelp . Flatten[gvvv\[Phi]];
	,
		lam4\[Phi]=Transpose[DS] . lam4\[Phi] . DS//Activate@TensorContract[Inactive@TensorProduct[DS,#],{{1,4}}]&//Transpose[#,{2,1,3,4}]&//Activate@TensorContract[Inactive@TensorProduct[DS,#],{{1,5}}]&//Transpose[#,{3,2,1,4}]&//SimplifySparse;
		lam3\[Phi]=Transpose[DS] . lam3\[Phi] . DS//Activate@TensorContract[Inactive@TensorProduct[DS,#],{{1,4}}]&//Transpose[#,{2,1,3}]&//SimplifySparse;
		gvss\[Phi]=Transpose[DV] . gvss\[Phi] . DS//Activate@TensorContract[Inactive@TensorProduct[DS,#],{{1,4}}]&//Transpose[#,{2,1,3}]&//SimplifySparse;
		Gvvs\[Phi]=Transpose[DV] . Gvvs\[Phi] . DS//Activate@TensorContract[Inactive@TensorProduct[DV,#],{{1,4}}]&//Transpose[#,{2,1,3}]&//SimplifySparse;
		gvvv\[Phi]=Transpose[DV] . gvvv\[Phi] . DV//Activate@TensorContract[Inactive@TensorProduct[DV,#],{{1,4}}]&//Transpose[#,{2,1,3}]&//SimplifySparse;
	];

	muab\[Phi]=vectorMass//SparseArray;
	muij\[Phi]=ScalarMass//SparseArray;

	If[DiagonalMatrixQAE[muab\[Phi]]==False,
		Print["The Vector mass-Matrix is not diagonal"];
	];


	If[DiagonalMatrixQAE[muij\[Phi]]==False,
		Print["The Scalar mass-Matrix is not diagonal"];
	];

];


