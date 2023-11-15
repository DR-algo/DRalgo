(* ::Package:: *)

(* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ *)

(* :Title: HEFT                                                          	*)

(*
      	This software is covered by the GNU General Public License 3.
       	Copyright (C) 2021-2023 Andreas Ekstedt
       	Copyright (C) 2021-2023 Philipp Schicho
       	Copyright (C) 2021-2023 Tuomas V.I. Tenkanen

*)

(* :Summary:	Integrates out particles with that get large field-dependent masses					*)	

(* ------------------------------------------------------------------------ *)


(* ::Section:: *)
(*Help functions*)


(*
	Trick to simplify tensors because Mathematica 13 sucks.
*)
SparseZero[s_,dims_] := Module[{},
	If[Length[s]==0,
(*If there are no array elements, as if for example when no hard scalars have been defined,
we want to return an array with 0 elements; this makes it easier to avoid unnecessary if statements down the line*)
		Return[SparseArray[{Table[1,{i,1,Length[dims]}]->0},dims]]
	,
		Return[SparseArray[s]];
	];
];


(*
	Trick to simplify tensors because Mathematica 13 sucks.
*)
HeavyTensor[s_,dims_] := Module[{},
	If[SparseArrayQ[s]==False,
(*If there are no array elements, as if for example when no hard scalars have been defined,
we want to return an array with 0 elements; this makes it easier to avoid unnecessary if statements down the line*)
		Return[SparseArray[{Table[1,{i,1,Length[dims]}]->0},dims]]
	,
		Return[SparseArray[s]];
	];
];


(* ::Section:: *)
(*Definition of model*)


(*
	Prepares the effective semi-soft/supersoft theory by creating hard and soft coupling tensors
*)
PrepareHET[HardScalarI_,HardSVectorI_]:=Module[{ListScalar=HardScalarI,ListVector=HardSVectorI},
	
	If[ValueQ[gvvvEP]==False,
		Print["You have to define the model, see UseUltraSoftTheory[], UseSoftTheory[], or DefineNewTensorsUS[]"];
		Print["Please also define the relevant VEVs with DefineVEVS[] before calculating the HET theory"];
		Return[];
	];
	
(*These particles will be integrated out*)
	HeavyScalarsHET=Transpose[List@ListScalar]; (*The list of all heavy scalars*)
	HeavyVectorsHET=Transpose[List@ListVector]; (*The list of all heavy scalars*)

	TotScalar=Table[{i},{i,1,nsEP}];
	TotVector=Table[{i},{i,1,nvEP}];
(*These are the particles that should not be integrated out*)	
	If[Length[HeavyScalarsHET]<1,
		LightScalarHET=TotScalar[[;;,1]];
		nHETS=1; (*Trick to avoid an exessive amount of if statements*)
		nHETSlight=nsEP;
	,
		LightScalarHET=Delete[TotScalar,HeavyScalarsHET][[;;,1]];
		nHETS=Length[HeavyScalarsHET];
		nHETSlight=nsEP-nHETS;
		HeavyScalarsHET=HeavyScalarsHET[[;;,1]];
	];
	If[Length[HeavyVectorsHET]<1,
		LightVectorHET=TotVector[[;;,1]];
		nHETV=1;(*Trick to avoid an exessive amount of if statements*)
		nHETVlight=nvEP;
	,
		LightVectorHET=Delete[TotVector,HeavyVectorsHET][[;;,1]];
		nHETV=Length[HeavyVectorsHET];
		nHETVlight=nvEP-nHETV;
		HeavyVectorsHET=HeavyVectorsHET[[;;,1]];
	];
	If[nHETVlight<1,
		nHETVlight=1;
		LightVectorHET={};
	];(*Trick to avoid errors in situations where all vector bosons are integrated out*)
(*We now extract masses for the hard particles*)
	\[Mu]ijHET=Table[\[Mu]ij\[Phi][[a,b]],{a,HeavyScalarsHET},{b,HeavyScalarsHET}]//SparseZero[#,{nHETS,nHETS}]&;
	\[Mu]abHET=Table[\[Mu]ab\[Phi][[a,b]],{a,HeavyVectorsHET},{b,HeavyVectorsHET}]//SparseZero[#,{nHETV,nHETV}]&;
	
	
(*Defining couplings between hard and soft particles*)
(*Convention: Capital letters are hard, otherwise a soft particle*)

(*Contraction with coupling-tensors and masses*)
(*We now have to integrate out hard particles in each topology*)
(*We treat off-diagonal masses perturbatively, so only the diagonal elements are used here*)	
	aS=Table[\[Mu]ijHET[[i,i]],{i,1,nHETS}]//SparseArray;
	av=Table[\[Mu]abHET[[i,i]],{i,1,nHETV}]//SparseArray;

(*For the effective potential we can loop over all particles (both light and heavy) if we set the light masses to zero*)
	avF=Table[0,{i,nvEP}]//SparseArray;
	asF=Table[0,{i,nsEP}]//SparseArray;
	
	If[Length[HeavyVectorsHET]>0,
		avF[[HeavyVectorsHET]]=av;
	];
	
	If[Length[HeavyScalarsHET]>0,
		asF[[HeavyScalarsHET]]=aS;
	];
	
	Return[]
];


(* ::Section::Closed:: *)
(*Effective potential*)


(*
	Calculates the tree-level effective potential.
*)
CalculateLOPotentialHET[]:=Module[{},
	If[verbose==True,Print["Calculating the Tree-Level Effective Potential"]];
	
	V1=\[Lambda]4EP . \[Phi]Vev . \[Phi]Vev . \[Phi]Vev . \[Phi]Vev;
	V2=\[Mu]ijEP . \[Phi]Vev . \[Phi]Vev;
	V3=\[Lambda]3EP . \[Phi]Vev . \[Phi]Vev . \[Phi]Vev;
	VHETLO=1/4! V1+V2/2!+V3/3!;
	];


(*
	Calculates the one-loop effective potential.
*)
CalculateNLOPotentialHET[]:=Module[{},
If[verbose==True,Print["Calculating the 1-Loop Effective Potential"]];

	V1=Sum[ALog[\[Mu]ijHET[[i,i]]],{i,1,nHETS}];
	V2=2*Sum[ALog[\[Mu]abHET[[i,i]]],{i,1,nHETV}];

	VHETNLO=V1+V2;
];


(*
	Calculates the two-loop effective potential.
*)
CalculateNNLOPotentialHET[]:=Module[{},
If[verbose==True,Print["Calculating the 2-Loop Effective Potential"]];


	
(*scalar-scalar bubble*)	
	ss=1/8 TensorProduct[\[Lambda]4\[Phi]];
	Vss=Sum[ss[[j,j,k,k]]fss[asF[[j]],asF[[k]]],{j,nsEP},{k,nsEP}];

(*scalar-scalar-scalar sunset*)	
	sss=1/12 TensorProduct[\[Lambda]3\[Phi],\[Lambda]3\[Phi]];
	Vsss=Sum[sss[[i,j,k,i,j,k]]fsss[asF[[i]],asF[[j]],asF[[k]]],{j,nsEP},{k,nsEP},{i,nsEP}];
	
(*scalar-vector bubble*)
	vs=1/2 TensorProduct[gvss\[Phi],gvss\[Phi]];
	Vvs=Sum[vs[[a,i,j,a,i,j]]fvs[asF[[i]],avF[[a]]],{a,nvEP},{j,nsEP},{i,nsEP}];


(*Vector-vector-scalar sunset diagrams*)
	vvs=1/4 TensorProduct[Gvvs\[Phi],Gvvs\[Phi]];
	Vvvs=Sum[vvs[[a,b,i,a,b,i]]fvvs[avF[[a]],avF[[b]],asF[[i]]],{a,nvEP},{b,nvEP},{i,nsEP}];
	
(*Scalar-Scalar-vector sunset diagrams*)									
	ssv=1/4 TensorProduct[gvss\[Phi],gvss\[Phi]];
	Vssv=Sum[ssv[[a,i,j,a,i,j]]fssv[asF[[i]],asF[[j]],avF[[a]]],{a,nvEP},{j,nsEP},{i,nsEP}];

(*vector-vector-vector sunset diagrams*)	
	vvv=1/12 TensorProduct[gvvv\[Phi],gvvv\[Phi]];
	Vvvv=Sum[vvv[[a,b,c,a,b,c]]fvvv[avF[[a]],avF[[b]],avF[[c]]],{a,nvEP},{b,nvEP},{c,nvEP}];
	
(*vector-vector bubble diagrams*)	
	vv=1/4 TensorProduct[gvvv\[Phi],gvvv\[Phi]];
	Vvv=Sum[vv[[a,b,c,a,b,c]]fvv[avF[[a]],avF[[b]]],{a,nvEP},{b,nvEP},{c,nvEP}];
	
(*ghost diagrams*)
	ggv=1/4 TensorProduct[gvvv\[Phi],gvvv\[Phi]];
	V\[Eta]\[Eta]v=Sum[ggv[[a,b,c,a,b,c]]f\[Eta]\[Eta]v[0,0,avF[[b]]],{a,nvEP},{b,nvEP},{c,nvEP}];
	

(*We now treat off-diagonal masses as small*)
	\[Mu]ijPert=\[Mu]ij\[Phi];
	\[Mu]ijPert[[LightScalarHET,LightScalarHET]]=0;
	\[Mu]ijPert=ArrayRules[\[Mu]ijPert]/.({x_Integer,y_Integer}->a_)/;Equal[x,y]->{x,y}->0//SparseArray[#,{nsEP,nsEP}]&;
	
	\[Mu]ab\[Phi]Pert=\[Mu]ab\[Phi];
	\[Mu]ab\[Phi]Pert[[LightVectorHET,LightVectorHET]]=0;
	\[Mu]ab\[Phi]Pert=ArrayRules[\[Mu]ab\[Phi]Pert]/.({x_Integer,y_Integer}->a_)/;Equal[x,y]->{x,y}->0//SparseArray[#,{nvEP,nvEP}]&;


(*Scalar contribution*)
	helpTens=Table[AD[a,b],{a,asF},{b,asF}]//SparseArray;
	helpTens=TensorProduct[helpTens,\[Mu]ijPert]//SparseArray//DiagonalTensor2[#,1,3]&;
	V1=Tr[helpTens . \[Mu]ijPert];

(*Vector contribution*)
	helpTens=Table[AD[a,b],{a,avF},{b,avF}]//SparseArray;
	helpTens=TensorProduct[helpTens,\[Mu]ab\[Phi]Pert]//SparseArray//DiagonalTensor2[#,1,3]&;
	V2=2*Tr[helpTens . \[Mu]ab\[Phi]Pert];
		
	VMix=V1+V2;
	
	
(*Potential*)
(*


	
	*)
	VHETNNLO= Vss+Vsss+Vvs+ Vvvs+  Vssv+ Vvvv+ Vvv+ V\[Eta]\[Eta]v+VMix;

];


(*
	Prints the effective potential.
*)
PrintActionHET[optP_]:=Module[{opt=optP},
	EffActionPrint=Switch[opt,"LO",VTotHET[[1]],"NLO",VTotHET[[2]],"NNLO",VTotHET[[3]]];

(*Printing Result*)
	OutputFormatDR[EffActionPrint]
];


(*
	Calculates the effective potential.
*)
CalculatePotentialHET[]:=Module[{},

	CalculateLOPotentialHET[];
	CalculateNLOPotentialHET[];
	CalculateNNLOPotentialHET[];

	VTotHET={VHETLO,VHETNLO,VHETNNLO};
	

];


(* ::Section:: *)
(*Scalar mass*)


(*
	Prints the effective potential.
*)
PrintScalarKineticHET[]:=Module[{},
	ScalarSelfEnergyHET[];

(*Printing Result*)
	Return[OutputFormatDR[ZSijHET]]
];


(*
	Scalar self-energy in the effective theory.
*)
ScalarSelfEnergyHET[]:=Module[{},
If[verbose,Print["Calculating Scalar Self-Energy"]];

(*Scalar-scalar bubble*)
	\[Lambda]3IJi=\[Lambda]3[[;;,;;,LightScalarHET]]//HeavyTensor[#,{nsEP,nsEP,nHETSlight}]&;
	TensHelp=Table[LSS[i,a],{i,asF},{a,asF}]//SparseArray;
	Temp=TensorProduct[TensHelp,\[Lambda]3IJi]//DiagonalTensor[#,2,3]&//DiagonalTensor[#,2,3]&;
	ContriSS=1/2*Contract[\[Lambda]3IJi,Temp,{{1,4},{2,5}}]//SimplifySparse;

(*Scalar-vector bubble*)
	gAIj=gvss\[Phi][[;;,LightScalarHET,;;]]//HeavyTensor[#,{nvEP,nHETSlight,nsEP}]&;
	TensHelp=Table[LSV[i,a],{i,asF},{a,avF}]//SparseArray;
	Temp=TensorProduct[TensHelp,gAIj]//DiagonalTensor[#,2,3]&//DiagonalTensor[#,2,4]&;
	ContriSV=TensorContract[gAIj . Temp,{1,3}];

(*Vector-vector bubble*)	
	GABi=Gvvs\[Phi][[;;,;;,LightScalarHET]]//HeavyTensor[#,{nvEP,nvEP,nHETSlight}]&;
	TensHelp=Table[LVV[i,a],{i,avF},{a,avF}]//SparseArray;
	Temp=TensorProduct[TensHelp,GABi]//DiagonalTensor[#,2,3]&//DiagonalTensor[#,2,3]&;
	ContriVV=1/2*Contract[GABi,Temp,{{1,4},{2,5}}]//SimplifySparse;
	
	ZSijHET=(- ContriSV- ContriVV-ContriSS)/2;
	

];


(* ::Section:: *)
(*Master Integrals*)


(*scalar-vector 2-point*)


LSV[x_,y_]:=(8 (-A[x]+A[y]))/(3 (x-y));


LSV[x_,x_]:=-((4 A[x] )/(3 x));
LSV[x_,0]:=-((8 A[x] )/(3 x))
LSV[0,0]=0;


(*vector-vector 2-point*)


LVV[x_,y_]:=(2 ((5 x-y) y^2 A[x]+x^2 (x-5 y) A[y]) )/(3 x (x-y)^3 y) ;


LVV[x_,x_]:=(5 A[x] )/(12 x^2);


LVV[x_,0]:=0;


LVV[0,x_]:=LVV[x,0];
LVV[0,0]:=0;


(*Scalar-scalar 2-point*)


LSS[x_,y_]:=((x+3 y) A[x]-(3 x+y) A[y]) /(3 (x-y)^3);


LSS[x_,x_]:=A[x]/(24 x^2);


LSS[x_,0]:=A[x] /(24 x^2);
LSS[0,x_]:=LSS[x,0];


LSS[0,0]:=0;
