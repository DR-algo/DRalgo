(* ::Package:: *)

(* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ *)

(* :Title: HEFT                                                          	*)

(*
       This software is covered by the GNU General Public License 3.
       Copyright (C) 2021-2026 Andreas Ekstedt
       Copyright (C) 2021-2026 Philipp Schicho
       Copyright (C) 2021-2026 Tuomas V.I. Tenkanen

*)

(* :Summary:	Integrates out particles with that get large field-dependent masses					*)	

(* ------------------------------------------------------------------------ *)


(* ::Section::Closed:: *)
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
PrepareHET[HardScalarIndices_,HardVectorIndices_]:=Module[
	{
		ListScalar=HardScalarIndices,
		ListVector=HardVectorIndices
	},
	
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
	muijHET=Table[muij\[Phi][[a,b]],{a,HeavyScalarsHET},{b,HeavyScalarsHET}]//SparseZero[#,{nHETS,nHETS}]&;
	muabHET=Table[muab\[Phi][[a,b]],{a,HeavyVectorsHET},{b,HeavyVectorsHET}]//SparseZero[#,{nHETV,nHETV}]&;
	
	
(*Defining couplings between hard and soft particles*)
(*Convention: Capital letters are hard, otherwise a soft particle*)

(*Contraction with coupling-tensors and masses*)
(*We now have to integrate out hard particles in each topology*)
(*We treat off-diagonal masses perturbatively, so only the diagonal elements are used here*)	
	aS=Table[muijHET[[i,i]],{i,1,nHETS}]//SparseArray;
	av=Table[muabHET[[i,i]],{i,1,nHETV}]//SparseArray;

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


(* ::Section:: *)
(*Effective potential*)


(*
	Calculates the tree-level effective potential.
*)
CalculateLOPotentialHET[]:=Module[{V1,V2,V3},
	If[verbose==True,Print["Calculating the Tree-Level Effective Potential"]];
	
	V1=lam4EP . \[Phi]Vev . \[Phi]Vev . \[Phi]Vev . \[Phi]Vev;
	V2=muijEP . \[Phi]Vev . \[Phi]Vev;
	V3=lam3EP . \[Phi]Vev . \[Phi]Vev . \[Phi]Vev;
	VHETLO=1/4! V1+V2/2!+V3/3!;
	];


(*
	Calculates the one-loop effective potential.
*)
CalculateNLOPotentialHET[]:=Module[{V1,V2},
If[verbose==True,Print["Calculating the 1-Loop Effective Potential"]];
	
	V1=Sum[OneLoopPotential3d[muijHET[[i,i]]],{i,1,nHETS}];
	V2=2*Sum[OneLoopPotential3d[muabHET[[i,i]]],{i,1,nHETV}];

	VHETNLO=V1+V2;
];


(*
	Calculates the two-loop effective potential.
*)
CalculateNNLOPotentialHET[]:=Module[{fvvv,f\[Eta]\[Eta]v,fssv,fvvs,fsss,fss,fvs,fvv,ss,sss,vs,vvs,ssv,vv,vvv,ggv,
							Vss,Vsss,Vvs,Vvvs,Vssv,Vvvv,Vvv,V\[Eta]\[Eta]v,
							muab\[Phi]Pert,muijPert,helpTens,V1,V2},
If[verbose==True,Print["Calculating the 2-Loop Effective Potential"]];


	(*Loading two-loop master integrals*)
	{fvvv,f\[Eta]\[Eta]v,fssv,fvvs,fsss,fss,fvs,fvv}=TwoLoopFunctions[mu3US];
(*scalar-scalar bubble*)	
	ss=1/8 TensorProduct[lam4\[Phi]];
	Vss=Sum[ss[[j,j,k,k]]fss[asF[[j]],asF[[k]]],{j,nsEP},{k,nsEP}];

(*scalar-scalar-scalar sunset*)	
	sss=1/12 TensorProduct[lam3\[Phi],lam3\[Phi]];
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
	muijPert=muij\[Phi];
	muijPert[[LightScalarHET,LightScalarHET]]=0;
	muijPert=ArrayRules[muijPert]/.({x_Integer,y_Integer}->a_)/;Equal[x,y]:>{x,y}->0//SparseArray[#,{nsEP,nsEP}]&;
	
	muab\[Phi]Pert=muab\[Phi];
	muab\[Phi]Pert[[LightVectorHET,LightVectorHET]]=0;
	muab\[Phi]Pert=ArrayRules[muab\[Phi]Pert]/.({x_Integer,y_Integer}->a_)/;Equal[x,y]:>{x,y}->0//SparseArray[#,{nvEP,nvEP}]&;

(*The one-loop diagram with two insertions of the off-diagonal masses P, -1/4 Sum[BubbleIntegral3d[m_i,m_l] P_il P_li],
	the O(P^2) term of the one-loop potential, with massless light lines (the hard region). Vectors give twice that.*)
(*Scalar contribution*)
	helpTens=Table[BubbleIntegral3d[a,b],{a,asF},{b,asF}]//SparseArray;
	V1=-1/4 Tr[(helpTens muijPert) . muijPert];

(*Vector contribution*)
	helpTens=Table[BubbleIntegral3d[a,b],{a,avF},{b,avF}]//SparseArray;
	V2=2*(-1/4) Tr[(helpTens muab\[Phi]Pert) . muab\[Phi]Pert];
		
	
(*Potential*)
(*	*)
	VHETNNLO= Vss+Vsss+Vvs+ Vvvs+  Vssv+ Vvvv+ Vvv+ V\[Eta]\[Eta]v+V1+V2;

];


(*
   Prints the effective potential at a given order.

   Options:
     "LO"   \[Dash] leading order (tree + 1-loop)
     "NLO"  \[Dash] next-to-leading order
     "All"  \[Dash] (default) sum of all available orders
*)
(* no argument \[RightArrow] same as "All" *)
PrintActionHET[] := PrintActionHET["All"];

PrintActionHET[opt_String] := Module[{effAction},

  effAction = Switch[opt,
    "LO",   Total[Take[VTotHET, UpTo[2]]],
    "NLO",  If[Length[VTotHET] >= 3, VTotHET[[3]], 0],
    "All",  If[VTotHET === {}, 0, Total[VTotHET]],
    _,      Message[PrintActionHET::badopt, opt]; Return[$Failed]
  ];

  OutputFormatDR[effAction]
];

(* custom message for invalid option *)
PrintActionHET::badopt = 
  "Invalid option `1`. Expected \"LO\", \"NLO\", or \"All\".";


(*
   Calculates the effective potential up to NNLO 
   and stores the results in VTotHET.
*)
CalculatePotentialHET[] := Module[{},
  
  (* compute contributions at each order *)
  CalculateLOPotentialHET[];
  CalculateNLOPotentialHET[];
  CalculateNNLOPotentialHET[];
  
  (* store results in global container *)
  VTotHET = {VHETLO, VHETNLO, VHETNNLO};
  
  Null
];


(* ::Section:: *)
(*Scalar mass*)


(*
   Prints the scalar kinetic terms in the effective theory.
*)
PrintScalarKineticHET[] := Module[{},
  
  (* compute scalar self-energy corrections *)
  ScalarSelfEnergyHET[];
  
  (* format and return kinetic term matrix *)
  OutputFormatDR[ZSijHET]
];


(*
	Scalar self-energy in the effective theory.
*)
ScalarSelfEnergyHET[]:=Module[{lam3IJi,gAIj,GABi,TensHelp,ContriSV,ContriVV,ContriSS},
If[verbose,Print["Calculating Scalar Self-Energy"]];

(*Scalar-scalar bubble*)
	lam3IJi=lam3\[Phi][[;;,;;,LightScalarHET]]//HeavyTensor[#,{nsEP,nsEP,nHETSlight}]&;
	TensHelp=Table[BubbleIntegral3dP2[i,a],{i,asF},{a,asF}]//SparseArray;
	TensHelp=TensorProduct[TensHelp,lam3IJi]//DiagonalTensor[#,2,3]&//DiagonalTensor[#,2,3]&;
	ContriSS=1/2*Contract[lam3IJi,TensHelp,{{1,4},{2,5}}]//SimplifySparse;

(*Scalar-vector bubble*)
	gAIj=gvss\[Phi][[;;,LightScalarHET,;;]]//HeavyTensor[#,{nvEP,nHETSlight,nsEP}]&;
	TensHelp=Table[BubbleIntegral3dP2SV[i,a],{i,asF},{a,avF}]//SparseArray;
	TensHelp=TensorProduct[TensHelp,gAIj]//DiagonalTensor[#,2,3]&//DiagonalTensor[#,2,4]&;
	ContriSV=TensorContract[gAIj . TensHelp,{1,3}];

(*Vector-vector bubble*)	
	GABi=Gvvs\[Phi][[;;,;;,LightScalarHET]]//HeavyTensor[#,{nvEP,nvEP,nHETSlight}]&;
	TensHelp=Table[BubbleIntegral3dP2VV[i,a],{i,avF},{a,avF}]//SparseArray;
	TensHelp=TensorProduct[TensHelp,GABi]//DiagonalTensor[#,2,3]&//DiagonalTensor[#,2,3]&;
	ContriVV=1/2*Contract[GABi,TensHelp,{{1,4},{2,5}}]//SimplifySparse;
	
	ZSijHET=(- ContriSV- ContriVV-ContriSS)/2;
	

];


