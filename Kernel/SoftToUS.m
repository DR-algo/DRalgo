(* ::Package:: *)

(* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ *)

(* :Title: SoftToUS                                                         *)

(*
      	This software is covered by the GNU General Public License 3.
       	Copyright (C) 2021-2026 Andreas Ekstedt
       	Copyright (C) 2021-2026 Philipp Schicho
       	Copyright (C) 2021-2026 Tuomas V.I. Tenkanen

*)

(* :Summary:	Dimensonal reduction from soft to supersoft scale.            *)	

(* ------------------------------------------------------------------------ *)


(* ::Section:: *)
(*Pressure calculation*)


(*
	Prints the result from SymmEnergy.
*)
PrintPressureUS[optP_]:=Module[{opt=optP},
	If[mode<2,
		SymmPrint=-SymmetricPhaseUSLO[]
	,
		SymmPrint=Switch[opt,
			"LO",SymmEnergyUS[[1]],
			"NLO",SymmEnergyUS[[2]],
			"NNLO",Message[DRalgo::failmsg, "No NNLO ultrasoft pressure implemented in DRalgo."]
		];
	];
(*Printing Result*)
	OutputFormatDR[SymmPrint]
];

PrintPressureUS[]:=Module[{},
	If[mode<2,
		SymmPrint=-SymmetricPhaseUSLO[]
	,
		SymmPrint=SymmEnergyUS[[1]]+SymmEnergyUS[[2]];
	];
(*Printing Result*)
	OutputFormatDR[SymmPrint]
];


(*
	Calculates the preassure in the ultrasoft theory. Only the preassure in the symmetric
	phae is calculated.
*)
SymmetricPhaseEnergyUS[]:=Module[{},
(*
	Counterterms are needed to calculate
	SymmetricPhaseNLO and SymmetricPhaseNNLO
*)
(*The minus signs is a convention to get the pressure*)
	Tot={-SymmetricPhaseUSLO[],-SymmetricPhaseUSNLO[]};
	SymmEnergyUS=Tot;
];


(*
	Calculates the 1-loop pressure in the ultrasoft theory.
*)
SymmetricPhaseUSLO[]:=Module[{ContriScalars,VLO},
	If[verbose,Print["Calculating Leading-Order \!\(\*SuperscriptBox[\(T\), \(4\)]\) Terms"]];
	
(*This just adds the m^3 term for all heavy scalars*)
	VLO=Sum[-1/(12 \[Pi]) muijL[[i,i]]^3,{i,1,Length[muijL]}];

	OutputFormatDR[VLO/.SubDebyeDecoupled]
];



(*
	Calculates the 2-loop pressure in the ultrasoft theory.
*)
SymmetricPhaseUSNLO[]:=Module[{fssv,Vss,TensHelp,Vssvs,Vssv,Vsss,heavyIndices,masses,VNLO},
	If[verbose,Print["Calculating NLO \!\(\*SuperscriptBox[\(T\), \(4\)]\) Terms"]];
(*Scalar-Scalar-Vector sunset diagram*)
	fssv=TwoLoopFunctions[mu3][[3]];

	Vss=1/8/(16 \[Pi]^2)*TensorContract[muijL . lamKTotal . muijL,{{1,2},{3,4}}];

	TensHelp=1/4 TensorProduct[gvssVTot,gvssVTot];
	Vssv=Sum[TensHelp[[a,i,j,a,i,j]]fssv[muijL[[i,i]]^2,muijL[[j,j]]^2,0],{a,nv},{i,nSH},{j,nSH}];

(*
	Scalar sunset diagram with two cubic couplings, -1/12 lam3^2 times the sunset integral,
	for the lines with at least one heavy scalar, where the light masses are zero.
	The sunset with only light lines belongs to the ultrasoft theory.
	The cubic couplings of two temporal scalars to a scalar are of one-loop order, and their sunset is of higher order.
*)
	heavyIndices=Flatten[HeavyScalars];
	masses=Table[If[MemberQ[heavyIndices,i],muijSNLOP[[i,i]],0],{i,1,ns}];
	Vsss=-1/12*Total[
		If[ContainsAny[#[[1]],heavyIndices],#[[2]]^2 SunsetIntegral3d[Sequence@@masses[[#[[1]]]],mu3],0]&/@Most[ArrayRules[lam3CSRedUS]]];

	VNLO= Vss+Vssv+Vsss;

	OutputFormatDR[VNLO/.SubDebyeDecoupled]
];



(* ::Section::Closed:: *)
(*Scalar masses*)


(*
	Scalar self-energy in the effective theory.
*)
ScalarSelfEnergySS[]:=Module[{ContriSS,ContriSS2},
If[verbose,Print["Calculating Scalar Self-Energy"]];

(*The p^2 terms of the bubbles; two heavy lines, with symmetry factor 1/2*)
	ContriSS=1/2*Simplify[Table[Sum[lam3Cx[[i,ii,jj]]lam3Cx[[j,ii,jj]]BubbleIntegral3dP2[muijL[[ii,ii]]^2,muijL[[jj,jj]]^2],{ii,1,nSH},{jj,1,nSH}],{i,1,nSL},{j,1,nSL}]];
(*One light and one heavy line, without symmetry factor as in ScalarMassSS*)
	ContriSS2=Simplify[Table[Sum[lam3Cy[[i,ii,jj]]lam3Cy[[j,ii,jj]]BubbleIntegral3dP2[0,muijL[[jj,jj]]^2],{ii,1,nSL},{jj,1,nSH}],{i,1,nSL},{j,1,nSL}]];

(*The field normalization is 1+2 ZSij; each light leg of a coupling contributes -ZSij*)
	ZSij=-(ContriSS+ContriSS2)/2;
];


(*
	Calculates the 2-loop scalar mass in the ultrasoft theory.
*)


ScalarMass2LoopSS[]:=Module[{MassHelp,TensHelp,TensHelp2,HeavyMasses,AllMasses
						,ContriInsertion,Contri1,Contri2,Contri3,Contri4,Contri5,Contri6,Contri7,ContriMix1,ContriMix2,ContriSE,muijSSCubic,
						ContriC1,ContriC2,ContriC3,ContriC4,muijTemp},
If[verbose,Print["Calculating 2-Loop Scalar Mass"]];

(*List of scalar masses*)
	HeavyMasses=Table[muijL[[n,n]],{n,1,nSH}]; (*A list of all the heavy-scalar masses*)
	AllMasses=Table[muijLS[[n,n]],{n,1,ns}]; (*A list of all the scalar masses*)

(*Diagrams without cubic couplings*)
(*Sunset diagram with two S^2H^2 quartic couplings*)
	MassHelp=Table[SunsetIntegral3d[a^2,b^2,0,mu3],{a,HeavyMasses},{b,HeavyMasses}]//SparseArray;
	TensHelp=lamK . Transpose[lamK,{4,2,3,1}]//DiagonalTensor[#,1,4]&//DiagonalTensor[#,1,5]&;
	Contri1=-1/2TensHelp . MassHelp//TensorContract[#,{2,4}]&;

(*Sunset diagram with a S^2H^2 quartic coupling and two H^2V couplings*)
	MassHelp=Table[(1/2+2Log[mu3/(2 a)]),{a,HeavyMasses}]//SparseArray;
	TensHelp=Transpose[Transpose[gAvss,{2,1,3}],{1,3,2}] . gAvss//TensorContract[#,{1,3}]&;
	TensHelp=MassHelp TensHelp . lamK//TensorContract[#,{1,2}]&;
	Contri2=1/(16 \[Pi]^2)*(1/2)TensHelp;

(*Sunset diagram with one S^2V^2 and two H^2V couplings*)
	MassHelp=Table[(-Log[mu3/(2 a)]),{a,HeavyMasses}]//SparseArray;
	TensHelp=gAvss . Transpose[gAvss,{2,1,3}]//SparseArray //DiagonalTensor[#,2,4]&//Transpose[#,{3,2,1}]&;
	TensHelp=TensHelp . MassHelp . HabijVL//TensorContract[#,{1,2}]&;
	Contri3=1/(16 \[Pi]^2)*(1/4)TensHelp;

(*Figure 8 diagram with two quartic couplings*)
	MassHelp=Table[1/(a+b),{a,HeavyMasses},{b,HeavyMasses}]//SparseArray;
	TensHelp=TensorProduct[DiagonalTensor[lam4K,3,4] . HeavyMasses,MassHelp]//DiagonalTensor[#,1,3]&//DiagonalTensor[#,1,3]&;
	Contri4=1/4*(1/(16 \[Pi]^2))*TensHelp . lamK//TensorContract[#,{1,2}]&;


	If[Length[lamx//Normal//Variables]==0&&Length[lamy//Normal//Variables]==0,
		Contri5=0;
		Contri6=0;
		Contri7=0;
	,
(*Sunset with two H^3S quartics*)
		TensHelp=Table[SunsetIntegral3d[a^2,b^2,n^2,mu3],{a,HeavyMasses},{b,HeavyMasses},{n,HeavyMasses}]//SparseArray;
		TensHelp2=TensorProduct[lamx,lamx]//DiagonalTensor2[#,2,6]&//DiagonalTensor2[#,3,6]&//DiagonalTensor2[#,4,6]&;
		Contri5=-1/3!*Total[TensHelp TensHelp2,-3];

(*Sunset with two HS^3 quartics*)
		TensHelp=Table[SunsetIntegral3d[a^2,0,0,mu3],{a,HeavyMasses}]//SparseArray;
		TensHelp2=TensorProduct[lamy,lamy]//TensorContract[#,{2,6}]&//TensorContract[#,{2,5}]&//DiagonalTensor2[#,2,4]&;
		Contri6=-1/2!*Total[TensHelp TensHelp2,-3];

(*Figure 8 bubbles with one HS^3quartic, and one H^3S quartic*)
		TensHelp=Table[c/b ,{c,HeavyMasses},{b,HeavyMasses}]//SparseArray;
		TensHelp2=Transpose[lamy,{1,2,4,3}] . lamx//DiagonalTensor2[#,3,4]&//DiagonalTensor2[#,4,5]&;
		Contri7=(1/(16 \[Pi]^2))/2*Total[TensHelp TensHelp2,-3];
	];

(*These two mixed contributions come from one-point reducible mixing between hard and soft particles*)

(*												 ------
Mixing diagram with two H^3S quartics: ----O---O------
													---*)
	TensHelp=Table[n/m^2*l ,{n,HeavyMasses},{m,HeavyMasses},{l,HeavyMasses}]//SparseArray;
	TensHelp2=TensorProduct[DiagonalTensor2[lamx,2,3],DiagonalTensor2[lamx,2,3]]//DiagonalTensor2[#,3,6]&//Flatten[#,{{2},{1},{4}}]&;
	ContriMix1=-1/4*1/(4 \[Pi])^2*Total[TensHelp TensHelp2,-3];

(*The contribution from off-diagonal HS masses*)
	ContriMix2=-Simplify[Table[Sum[(muijMix[[i,m]])muijL[[m,m]]^-2 muijMix[[j,m]],{m,1,nSH}],{i,1,nSL},{j,1,nSL}]];

	If[Length[lam3//Normal//Variables]==0||nv>=nSH,
		muijSSCubic=0;
	,
(*The contribution from cubic couplings. There are so many diagrams here that it is not worth dividing the cubics into
	H^2S HS^2 etc. Insdead lam3CTot is a master tensor that has all scalar (heavy and light) interactions. The
	sunset integrals, which vanish if all lines are massless, then remove the pure-light contributions.
	*)

(*Sunset diagram with four cubic couplings*)		
		MassHelp=Table[SunsetIntegral3d[{j,i},{k,l},m,mu3],{j,AllMasses^2},{i,AllMasses^2},{k,AllMasses^2},{l,AllMasses^2},{m,AllMasses^2}]//SparseArray//SimplifySparse;
		TensHelp=TensorProduct[lam3CTot,lam3CTot]//DiagonalTensor2[#,1,4]&//TensorProduct[#,lam3CTot]&//DiagonalTensor2[#,2,6]&//DiagonalTensor2[#,4,6]&;
		TensHelp=TensorProduct[TensHelp,lam3CTot]//DiagonalTensor2[#,4,7]&//DiagonalTensor2[#,5,7]&;
		ContriC1=-1/2*Total[MassHelp TensHelp,-3][[LightScalar[[;;,1]],LightScalar[[;;,1]]]]//SimplifySparse;

(*Sunset diagram with one quartic and two cubic couplings*)		
		MassHelp=Table[SunsetIntegral3d[{i,j},m,n,mu3],{i,AllMasses^2},{j,AllMasses^2},{m,AllMasses^2},{n,AllMasses^2}]//SparseArray//SimplifySparse;
		TensHelp=TensorProduct[lam4Tot,lam3CTot]//DiagonalTensor2[#,4,6]&//DiagonalTensor2[#,4,5]&;
		TensHelp=TensorProduct[TensHelp,lam3CTot]//DiagonalTensor2[#,4,7]&//DiagonalTensor2[#,5,7]&;
		ContriC2=1/2*Total[MassHelp TensHelp,-3][[LightScalar[[;;,1]],LightScalar[[;;,1]]]]//SimplifySparse;

(*Sunset diagram with one quartic (light scalars both connect to the same quartic) and two cubic couplings*)
		MassHelp=Table[SunsetIntegral3d[{i,j},m,n,mu3],{i,AllMasses^2},{j,AllMasses^2},{m,AllMasses^2},{n,AllMasses^2}]//SparseArray//SimplifySparse;
		TensHelp=TensorProduct[lam3CTot,lam3CTot]//DiagonalTensor2[#,3,4]&//DiagonalTensor2[#,3,4]&;
		TensHelp=TensorProduct[TensHelp,lam4Tot]//DiagonalTensor2[#,3,5]&//DiagonalTensor2[#,4,7]&;
		ContriC3=1/4*Total[MassHelp TensHelp,-3][[LightScalar[[;;,1]],LightScalar[[;;,1]]]]//SimplifySparse;
		
(*Two bubbles connected in the middle by a quartic, and connecting to the external lines with a cubic coupling*)
		MassHelp=Table[BubbleIntegral3d[i,j] BubbleIntegral3d[k,l],{i,AllMasses^2},{j,AllMasses^2},{k,AllMasses^2},{l,AllMasses^2}]//SparseArray//SimplifySparse;
		TensHelp=TensorProduct[lam3CTot,lam4Tot]//DiagonalTensor2[#,1,4]&//DiagonalTensor2[#,2,4]&;
		TensHelp=TensorProduct[TensHelp, lam3CTot]//DiagonalTensor2[#,4,6]&//DiagonalTensor2[#,4,6]&;
		ContriC4=1/8*Total[MassHelp TensHelp,-3][[LightScalar[[;;,1]],LightScalar[[;;,1]]]]//SimplifySparse;
		
	
		muijSSCubic=(ContriC1+2*ContriC2+ContriC3+ContriC4)//Simplify//SparseArray;
	];

(*This is the contribution from field-renormalization, where the Z factor multiplies the entire LO (tree-level+1-loop) mass*)
	muijTemp=muijLight+muijSSLO//SparseArray;
	ContriSE=-(ZSij . muijTemp+muijTemp . ZSij);

(*Bubble with one light and one heavy line, with one insertion of the light mass: the tree-level mass and the momentum-independent
	one-loop terms, from the tadpoles of the heavy fields and the temporal scalars*)
	ContriInsertion=Simplify[Table[Sum[lam3Cy[[i,a,hh]](muijLight+muijSSLOTadpoles)[[a,b]]lam3Cy[[j,b,hh]]TriangleIntegral3d[0,0,muijL[[hh,hh]]^2],{a,1,nSL},{b,1,nSL},{hh,1,nSH}],{i,1,nSL},{j,1,nSL}]];

	muijSSNLO=(Contri1+Contri2+Contri3+ Contri4+Contri5+Contri6+Contri7+ContriMix1+ContriMix2+ContriSE+ContriInsertion+muijSSCubic)//Simplify//SparseArray;

];



(*
	Calculates the 1-loop scalar mass.
*)
ScalarMassSS[]:=Module[{ContriSS,ContriSS2,ContriTadpole,ContriSS3},
If[verbose,Print["Calculating 1-Loop Scalar Mass"]];

(*Two H^2S cubics*)
	ContriSS=1/2*Simplify[Table[Sum[lam3Cx[[i,ii,jj]]lam3Cx[[j,ii,jj]]BubbleIntegral3d[muijL[[ii,ii]]^2,muijL[[jj,jj]]^2],{ii,1,nSH},{jj,1,nSH}],{i,1,nSL},{j,1,nSL}]];

(*Two HS^2 cubics*)
	ContriSS2=Simplify[Table[Sum[lam3Cy[[i,ii,jj]]lam3Cy[[j,ii,jj]]BubbleIntegral3d[0,muijL[[jj,jj]]^2],{ii,1,nSL},{jj,1,nSH}],{i,1,nSL},{j,1,nSL}]];

(*One HS^2 cubic with the tree-level and one-loop tadpoles of the heavy fields*)
	ContriTadpole=Table[Sum[lam3Cy[[i,j,ll]]TadPoleHeavySS[[ll]]/(muijL[[ll,ll]]^2),{ll,1,nSH}],{i,1,nSL},{j,1,nSL}];

(*A H^2S^2 quartic*)
	ContriSS3=-1/2 Simplify[Table[Sum[TadpoleIntegral3d[muijL[[a,a]]^2]lamK[[a,a,i,j]],{a,1,nSH}],{i,1,nSL},{j,1,nSL}]];

	muijSSLO=-ContriSS3-ContriSS-ContriSS2-ContriTadpole+TadpoleShiftMassSS[]//SparseArray;
(*The terms without momentum dependence, which are inserted on the light lines of the two-loop masses*)
	muijSSLOTadpoles=-ContriSS3-ContriTadpole//SparseArray;

];


(*
	This is the mass that pops up if the heavy-scalar line carry a soft momenta
*)
HeavyScalarMassSS[]:=Module[{},
	fac=1/(4 \[Pi]);
	Contri1=1/2  *fac*Simplify[Table[Sum[lam3CHeavy[[i,ii,jj]]lam3CHeavy[[j,ii,jj]]/(muijL[[ii,ii]]+muijL[[jj,jj]]),{ii,1,nSH},{jj,1,nSH}],{i,1,nSH},{j,1,nSH}]];
	fac=-1/(4 \[Pi]);
	Contri2=-1/2*fac*Simplify[Table[Sum[lamKTotal[[i,j,ii,ii]] muijL[[ii,ii]],{ii,1,nSH}],{i,1,nSH},{j,1,nSH}]]//FullSimplify;
	fac=1/(4 \[Pi]);
	Contri3=fac*Simplify[Table[Sum[lam3Cx[[jj,ii,i]]lam3Cx[[jj,ii,j]]/(muijL[[ii,ii]]),{ii,1,nSH},{jj,1,nSL}],{i,1,nSH},{j,1,nSH}]];
(*The tadpole shifts the heavy field by -TadPoleHeavySS/M^2; like the other terms, Contri4 is minus the mass correction*)
	Contri4=Simplify[Table[Sum[TadPoleHeavySS[[ii]]lam3CHeavy[[i,j,ii]]/(muijL[[ii,ii]]^2),{ii,1,nSH}],{i,1,nSH},{j,1,nSH}]];

	HeavyScalarMass=Contri1+Contri2+Contri3+Contri4;

	HelpList=DeleteDuplicates@Flatten@Simplify[ HeavyScalarMass]//Sort;
	HelpVar=Table[ muH[a],{a,1,Delete[HelpList,1]//Length}];
	HelpVarMod=RelationsBVariables[HelpList,HelpVar];
	HelpSolveEffectiveHardM=Table[{Delete[HelpList,1][[a]]->HelpVarMod[[a]]},{a,1,Delete[HelpList,1]//Length}]//Flatten;
	muHEff=HeavyScalarMass//Normal//ReplaceAll[#,HelpSolveEffectiveHardM]&//SparseArray;

];


(* ::Section:: *)
(*Effective couplings*)


(*
	Matching of scalar-cubic couplings. lam3Cy corresponds to light*light*heavy scalar cubic, and lam3Cx corresponds to light*heavy*heavy scalar coupling.
*)

ScalarCubicsSS[]:=Module[{HeavyMasses,ContriLoop,ContriMixed,ContriSE,cubicTree,mixingCubic,tadpoleCubic},
If[verbose,Print["Calculating Scalar Cubic Couplings"]];

	If[Length[lam3//Normal//Variables]==0,
		lam3CSSS=0;
	,
		HeavyMasses=Table[muijL[[n,n]],{n,1,nSH}]; (*A list of all the heavy-particle masses*)

(*Bubble and triangle diagrams with at least one heavy line*)
		ContriLoop=OneLoopCubicSS[Range[nSL],Range[nSL],Range[nSL]];

(*
	One-point reducible diagrams: a light-heavy mixing whose heavy field ends at a H^2S cubic, and the
	tadpoles of the heavy fields at a HS^3 quartic.
*)
		mixingCubic[mix_]:=-3 Symmetrize[Table[Sum[mix[[i,h]]lam3Cy[[j,k,h]]/HeavyMasses[[h]]^2,{h,1,nSH}],{i,1,nSL},{j,1,nSL},{k,1,nSL}],Symmetric[{1,2,3}]];
		tadpoleCubic[tadpole_]:=-Table[Sum[lamy[[i,j,k,h]]tadpole[[h]]/HeavyMasses[[h]]^2,{h,1,nSH}],{i,1,nSL},{j,1,nSL},{k,1,nSL}];
		ContriMixed=mixingCubic[MixingLH[TadPoleHeavySS]]+tadpoleCubic[TadPoleHeavySS];

(*The field normalization acts on the tree-level cubic of the light fields, including the tree-level tadpoles*)
		cubicTree=lam3CLight+mixingCubic[TadpoleShiftLH[TadPoleHeavy]]+tadpoleCubic[TadPoleHeavy];
		ContriSE=-3 Symmetrize[ZSij . cubicTree,Symmetric[{1,2,3}]];

		lam3CSSS=ContriLoop+ContriMixed+ContriSE//SparseArray//SimplifySparse;
	];

];


(*
	Adds tree-level Corrections to Scalar-quartic couplings.
*)
ScalarModifiedSS[]:=Module[{},
If[verbose,Print["Adding tree-level Corrections to Scalar Quartics"]];

	ContriTLTemp=-Simplify[Table[Sum[1/(muijL[[n,n]]^2)lam3Cy[[i,j,n]]lam3Cy[[k,l,n]],{n,1,nSH}],{i,1,nSL},{j,1,nSL},{k,1,nSL},{l,1,nSL}]];
	ContriTL=Table[ContriTLTemp[[i,j,k,l]]+ContriTLTemp[[i,k,j,l]]+ContriTLTemp[[i,l,j,k]],{i,1,nSL},{j,1,nSL},{k,1,nSL},{l,1,nSL}];
	lamMod=lam4S+ContriTL;
	lam4S=lamMod;

	ContriTLTemp=-Simplify[Table[Sum[1/(muijL[[n,n]]^2)lam3CHeavy[[ii,jj,n]]lam3Cy[[k,l,n]],{n,1,nSH}],{ii,1,nSH},{jj,1,nSH},{k,1,nSL},{l,1,nSL}]];
	ContriTLTemp2=-Simplify[Table[Sum[1/(muijL[[n,n]]^2)lam3Cx[[k,ii,n]]lam3Cx[[l,jj,n]],{n,1,nSH}],{ii,1,nSH},{jj,1,nSH},{k,1,nSL},{l,1,nSL}]];
	ContriTL=Table[ContriTLTemp[[ii,jj,k,l]]+ContriTLTemp2[[ii,jj,k,l]]+ContriTLTemp2[[jj,ii,k,l]],{ii,1,nSH},{jj,1,nSH},{k,1,nSL},{l,1,nSL}];
	lamMod=lamK+ContriTL;
	lamK=lamMod;

	ContriTLTemp=-Simplify[Table[Sum[1/(muijL[[n,n]]^2)lam3Cx[[i,jj,n]]lam3Cy[[k,l,n]],{n,1,nSH}],{i,1,nSL},{k,1,nSL},{l,1,nSL},{jj,1,nSH}]];
	ContriTL=Table[ContriTLTemp[[i,k,l,jj]]+ContriTLTemp[[k,i,l,jj]]+ContriTLTemp[[l,k,i,jj]],{i,1,nSL},{k,1,nSL},{l,1,nSL},{jj,1,nSH}];
	lamMod=lamy+ContriTL;
	lamy=lamMod;


	ContriTLTemp=-Simplify[Table[Sum[1/(muijL[[n,n]]^2)lam3CHeavy[[ii,jj,n]]lam3Cx[[k,ll,n]],{n,1,nSH}],{k,1,nSL},{ii,1,nSH},{jj,1,nSH},{ll,1,nSH}]];
	ContriTL=Table[ContriTLTemp[[k,ii,jj,ll]]+ContriTLTemp[[k,ii,ll,jj]]+ContriTLTemp[[k,ll,jj,ii]],{k,1,nSL},{ii,1,nSH},{jj,1,nSH},{ll,1,nSH}];
	lamMod=lamx+ContriTL;
	lamx=lamMod;

	MassHelp[0]:=0;
	MassHelp[x_]:=1/(x);

	ContriTLTemp=-Simplify[Table[Sum[MassHelp[muijLS[[n,n]]]^2 lam3CTot[[i,j,n]]lam3CTot[[k,l,n]],{n,1,ns}],{i,1,ns},{j,1,ns},{k,1,ns},{l,1,ns}]];
	ContriTL=Table[ContriTLTemp[[i,j,k,l]]+ContriTLTemp[[i,k,j,l]]+ContriTLTemp[[i,l,k,j]],{i,1,ns},{j,1,ns},{k,1,ns},{l,1,ns}];
	lam3DSp+ContriTL;

	lam4Tot[[HeavyScalars[[;;,1]],HeavyScalars[[;;,1]],HeavyScalars[[;;,1]],HeavyScalars[[;;,1]]]]=lam3DSp[[HeavyScalars[[;;,1]],HeavyScalars[[;;,1]],HeavyScalars[[;;,1]],HeavyScalars[[;;,1]]]];


];


(*
	Calculates the scalar quartic in the ultrasoft theory.
*)
ScalarQuarticSS[]:=Module[{lam3CyLoop,ContriTree,ContriTadpole,ContriMixing,lamKEff,lamyEff,lamyEff2,HeavyMasses,
					AllMasses,TensHelp,MassHelp,HeavyScal,ContriSE,
					ContriCubic,ContriQuarCubic,ContriSS,ContriSS2,ContriSS3,ContriTL,ContriTadpole1,ContriTadpole2},
If[verbose,Print["Calculating Scalar Quartics"]];

(*Quartics with a various number of heavy-scalar legs*)
	lamKEff=lamK//SparseArray;
	lamyEff=lamy//SparseArray;
	lamyEff2=lamy;

(*List of scalar masses*)
	HeavyMasses=Table[muijL[[n,n]],{n,1,nSH}]; (*A list of all the heavy-particle masses*)
	HeavyScal=Table[muijL[[n,n]],{n,nv+1,nSH}]; (*A list of all the heavy-scalar masses*)
	AllMasses=Table[muijLS[[n,n]],{n,1,ns}]; (*A list of all the scalar masses*)
	
(*A bubble with two H^2S^2 quartics*)
	MassHelp=Table[1/(n+m),{n,HeavyMasses},{m,HeavyMasses}]//SparseArray//Flatten[#,{1,2}]&;
	TensHelp=Flatten[lamKEff,{{1,2}}]//SparseArray;
	TensHelp=1/2*1/(4 \[Pi])*Transpose[TensHelp,{3,2,1}] . (MassHelp\[NonBreakingSpace]TensHelp)//SparseArray;
	ContriSS=3Symmetrize[TensHelp,Symmetric[{1,2,3,4}]]//SparseArray//SimplifySparse;
	
	
	If[Length[lamy//Normal//Variables]==0,
		ContriSS2=EmptyArray[{nSL,nSL,nSL,nSL}];
		ContriSS3=EmptyArray[{nSL,nSL,nSL,nSL}];
	,

(*A bubble with two two HS^3 quartics*)
		MassHelp=Table[1/(n),{n,HeavyMasses}]//SparseArray;
		TensHelp=MassHelp Transpose[lamyEff,{4,2,3,1 }]//Flatten[#,{2,1}]& ;
		TensHelp=1/(4 \[Pi]) Transpose[Flatten[lamyEff,{3,4}],{3,2,1}] . TensHelp;
		ContriSS2=3Symmetrize[TensHelp,Symmetric[{1,2,3,4}]]//SparseArray//SimplifySparse;

(*A one-point reducible diagram where (for one external line) a H^3S quartic forms a bubble and where the
	remaining H line connects to a HS^3 quartic*)
	
		MassHelp=Table[n/m^2,{m,HeavyMasses},{n,HeavyMasses}]//SparseArray;
		TensHelp=Transpose[lamx,{4,2,3,1}]//Table[i,{i,#}]&//Table[#[[i,i]],{i,1,Length[#]}]&//SparseArray;
		TensHelp=MassHelp . TensHelp//Table[i,{i,#}]&//Table[#[[i,i]],{i,1,Length[#]}]&//SparseArray;
		ContriSS3=-1/(4 \[Pi])*1/2*4*Symmetrize[lamyEff2 . TensHelp,Symmetric[{1,2,3,4}]]//SparseArray//SimplifySparse;
	];


	
	If[Length[lam3//Normal//Variables]==0,
		lam3DSS=-ContriSS-ContriSS2-ContriSS3;
	,
(*A tree-level two-cubic diagram*)
		MassHelp=Table[1/(n^2),{n,HeavyMasses}]//SparseArray;
		TensHelp=TensorProduct[lam3Cy,lam3Cy]//DiagonalTensor2[#,3,6]&;
		TensHelp=MassHelp . TensHelp;
		ContriTree=3 Symmetrize[TensHelp,Symmetric[{1,2,3,4}]]//SparseArray//SimplifySparse;
		
(*A tree-level two-cubic diagram with a loop correction*)		
		MassHelp=Table[1/(n^2*m^2),{n,HeavyMasses},{m,HeavyMasses}]//SparseArray;
		MassHelp=muHEff MassHelp;
		TensHelp+=lam3Cy . MassHelp . Transpose[lam3Cy,{3,2,1}];

(*A tree-level two-cubic diagram with a loop correction to one of the light-light-heavy vertices*)
		lam3CyLoop=VertexCorrectionLLH[];
		MassHelp=DiagonalMatrix[Table[1/(n^2),{n,HeavyMasses}]]//SparseArray;
		TensHelp+=lam3CyLoop . MassHelp . Transpose[lam3Cy,{3,2,1}]+lam3Cy . MassHelp . Transpose[lam3CyLoop,{3,2,1}];
		ContriTL=3 Symmetrize[TensHelp,Symmetric[{1,2,3,4}]]//SparseArray//SimplifySparse;

(*A one-point reducible diagram where one leg has a tadpole insertion, with HS^2 and H^2S^2 couplings*)
		ContriTadpole[tadpole_]:=Module[{help},
			MassHelp=Table[1/(n^2*m^2),{n,HeavyMasses},{m,HeavyMasses}]//SparseArray;
			help=TensorProduct[lamK,lam3Cy,tadpole]//DiagonalTensor2[#,1,7]&//DiagonalTensor2[#,2,7]&;
			help=-MassHelp . help//TensorContract[#,{1,2}]&;
			2*3 Symmetrize[help,Symmetric[{1,2,3,4}]]//SparseArray//SimplifySparse
		];
		ContriTadpole2=ContriTadpole[TadPoleHeavySS];

(*One-point reducible diagrams with a light-heavy mixing from loops and from the tadpoles of the heavy fields*)
		ContriMixing=MixingQuarticSS[MixingLH[TadPoleHeavySS]];

(*
	The self-energy contribution. The field normalization acts on the tree-level quartic of the light fields,
	including the tree-level exchange of heavy fields and the tree-level tadpoles.
*)
		ContriSE=-ZSij . (lam4S-ContriTree-ContriTadpole[TadPoleHeavy]+MixingQuarticSS[TadpoleShiftLH[TadPoleHeavy]])//SparseArray;
		ContriSE=ContriSE+Transpose[ContriSE,{2,1,3,4}]+Transpose[ContriSE,{3,1,2,4}]+Transpose[ContriSE,{4,1,2,3}]//SparseArray//SimplifySparse;

	If[nSH>nv,
(*Here I I don't separate the quartics/cubic couplings into light and heavy. Instead I do the most general calculation
	for all the scalars and then project out the relevant contributions in the end.*)	

(*A box diagram with 4 cubic couplings*)
		MassHelp=Table[BoxIntegral3d[n,m,l,w],{n,AllMasses^2},{m,AllMasses^2},{l,AllMasses^2},{w,AllMasses^2}]//SparseArray;	
		TensHelp=TensorProduct[lam3CTot,lam3CTot,lam3CTot]//DiagonalTensor2[#,2,5]&//DiagonalTensor2[#,5,7]&;
		TensHelp=TensorProduct[TensHelp,lam3CTot]//DiagonalTensor2[#,4,9]&//DiagonalTensor2[#,7,8]&;
		TensHelp=Total[MassHelp TensHelp,-5][[LightScalar[[;;,1]],LightScalar[[;;,1]],LightScalar[[;;,1]],LightScalar[[;;,1]]]];
(*The three orderings of the external legs; the box enters with the opposite sign of the triangle*)
		ContriCubic=3 Symmetrize[TensHelp,Symmetric[{1,2,3,4}]]//SparseArray//SimplifySparse;

(*A triangle diagram with two cubic and one quartic couplings*)				
		MassHelp=Table[TriangleIntegral3d[n,m,l],{n,AllMasses^2},{m,AllMasses^2},{l,AllMasses^2}]//SparseArray;	
		TensHelp=TensorProduct[lam4Tot,lam3CTot,lam3CTot]//DiagonalTensor2[#,3,5]&//DiagonalTensor2[#,4,7]&//DiagonalTensor2[#,6,7]&;
		TensHelp=Total[MassHelp TensHelp,-5][[LightScalar[[;;,1]],LightScalar[[;;,1]],LightScalar[[;;,1]],LightScalar[[;;,1]]]];
		ContriQuarCubic=-3!Symmetrize[TensHelp,Symmetric[{1,2,3,4}]]//SparseArray//SimplifySparse;
		
		lam3DSS=  ContriSE-ContriTL-ContriTadpole2+ContriMixing- ContriSS-ContriSS2- ContriQuarCubic- ContriCubic;
	,
		lam3DSS=ContriSE-ContriTL-ContriTadpole2+ContriMixing- ContriSS-ContriSS2;
	];

	];

];


(*
	Calculates non-abelian couplings in the ultrasoft theory.
*)
NonAbelianCouplingSS[]:=Module[{},
	ContriAnomVV= gvvvSS . Transpose[ZLij]//SparseArray//SimplifySparse;
(*Simplify[Table[Sum[ ZLij[[c,d]]gvvvSS[[a,b,d]],{d,1,nv}],{a,1,nv},{b,1,nv},{c,1,nv}]];*)
	GgvvvSS=-ContriAnomVV;
];


(*
	Calculates vector couplings in the ultrasoft theory.
*)
ScalarVectorCouplingSS[]:=Module[{},
If[verbose,Print["Calculating Vector-Scalar Couplings"]];


(* Self-Energy contribution*)
	ContriSEVector=-ZLij . HabijVL-Transpose[ZLij . HabijVL,{2,1,3,4}];
(*ContriSEVector=-Simplify[Table[Sum[ZLij[[a,c]](HabijVL[[c,b,i,j]])+ZLij[[b,c]](HabijVL[[a,c,i,j]]),{c,1,nv}],{a,1,nv},{b,1,nv},{i,1,nSL},{j,1,nSL}]];*)


	GvvssTSS= ContriSEVector//Simplify;

];


(*
	Calculates the transverse-vector self energy.
*)
VectorSelfEnergySS[]:=Module[{},
If[verbose,Print["Calculating Vector Self-Energy"]];

	SelfEnergySS=-1/2*(-1/(24 \[Pi]));
	ContriSS=SelfEnergySS*Simplify[Table[Sum[HabijA[[a,b,i,i]]/muijL[[i,i]],{i,1,nSH}],{a,1,nv},{b,1,nv}]];
	ZLij=-ContriSS/2;
];


(*
	Prints the beta functions in the ultrasoft theory.
*)
BetaFunctions3DUS[]:=Module[{},
If[verbose,Print["Finding SuperSoft \[Beta]-functions"]];


	VarGauge=Join[muijLight//Normal//Variables]//DeleteDuplicates;
	SubGauge=Table[c->Symbol[ToString[c]<>ToString["US"]],{c,VarGauge}];

(*This is the epsilon pole of the 3d sunset diagram. This is the only divergence in 3d.*)
	I111=1/(4)1/(16 \[Pi]^2);

	Coeff=1/2*2;
	Contri1=Coeff*Flatten[TensorContract[HabijL,{1,2}]] . Flatten[lam4Light,{1,2}];

	Coeff=-1/4*(-1);
	Contri2=Coeff*Flatten[TensorContract[HabijL,{3,4}]] . Flatten[HabijVL,{1,2}];

	Coeff=1/3!;
	Contri3=Coeff *Flatten[lam4Light,{{1},{2,3,4}}] . Flatten[lam4Light,{{1,2,3}}];

	Coeff=1/2*((1+1/2));
	Contri4=Coeff *Flatten[HabijVL,{{3},{1,2,4}}] . Flatten[HabijVL,{{1,2,4},{3}}];

	Coeff=(-1)(-1/4) /2;
	GabcdV2=TensorContract[GabcdV,{2,4}];
	Contri5=Coeff*Flatten[GabcdV2] . Flatten[HabijVL,{1,2}];
 
	Coeff=(-1)(-1)1/4*(20/4);
	Contri6=Coeff*Flatten[GabcdV2] . Flatten[HabijVL,{1,2}];

	SubGauge=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,GaugeCouplingNames}];

	\[Delta]mu3dS=Contri1+   Contri2+ Contri3+  Contri4+ Contri5+  Contri6;(*Sum of all the diagrams*)
	\[Beta]mu3ijS=4*I111* \[Delta]mu3dS//Normal //ReplaceAll[#,SubGauge]&//Simplify; (*The scalar-mass counter-term*)

(*Printing Result*)
(* Scalar Mass*)
	VarGauge=Join[muijLight//Normal//Variables]//DeleteDuplicates;
	SubGauge=Table[c->Symbol[ToString[c]<>ToString["temp"]],{c,VarGauge}];
	SubGauge2=Table[Symbol[ToString[c]<>ToString["temp"]]->Symbol[ToString[c]<>ToString["US"]],{c,VarGauge}];

	GaugeHelp=Table[Symbol[ToString[c]<>ToString["3d"]],{c,GaugeCouplingNames}];
	VarUS=Join[{lam3CLight//Normal//Variables,lam4Light//Normal//Variables,GaugeHelp//Normal//Variables}]//DeleteDuplicates//Flatten[#,1]&;
	SubUS=Table[c->Symbol[ToString[c]<>ToString["US"]],{c,VarUS}];

	lam4p=muijLight//Normal//ReplaceAll[#,SubGauge]&//SparseArray;
	SolVar=\[Beta]mu3ijS-lam4p//Normal;
	QuarticVar=lam4p//Normal//Variables;
	ResMass=Solve[SolVar==0,QuarticVar]/.SubGauge2/.SubUS//Flatten[#,1]&; (*Finds the beta-function for each scalar mass*)

(*Printing Result*)
	PrintPre=ResMass//Normal//FullSimplify//DeleteDuplicates;
	
(*Tadpoles*)
	\[CapitalLambda]gHelp=TensorContract[HabijL,{{1,2}}]//SimplifySparse; 
	Contri10=1/3!*I111*Activate@TensorContract[Inactive@TensorProduct[lam4Light//SparseArray,lam3CLight//SparseArray],{{2,5},{3,6},{4,7}}];
	Contri11=1/2*I111*2Activate@TensorContract[Inactive@TensorProduct[lam3CLight//SparseArray,\[CapitalLambda]gHelp//SparseArray],{{2,4},{3,5}}];

	VarGauge=GaugeCouplingNames;
	SubGauge=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarGauge}];

	\[Delta]muTadpole3dUS=Contri10+Contri11//Normal//ReplaceAll[#,SubGauge]&//SparseArray;(*Same deal for tadpoles*)
	\[Beta]muTadpoleUS=4*\[Delta]muTadpole3dUS //Simplify;

	VarGauge=Join[TadPoleLight//Normal//Variables]//DeleteDuplicates;
	SubGauge=Table[c->Symbol[ToString[c]<>ToString["temp"]],{c,VarGauge}];
	SubGauge2=Table[Symbol[ToString[c]<>ToString["temp"]]->Symbol[ToString[c]<>ToString["US"]],{c,VarGauge}];

	lam4p=TadPoleLight//Normal//ReplaceAll[#,SubGauge]&//SparseArray;
	SolVar=\[Beta]muTadpoleUS-lam4p//Normal;
	QuarticVar=lam4p//Normal//Variables;
	ResTadpole=Solve[SolVar==0,QuarticVar]/.SubGauge2/.SubUS//Flatten[#,1]&; (*Finds the beta-function for each scalar mass*)

	PrintPre=Join[PrintPre,ResTadpole]//Normal//FullSimplify//DeleteDuplicates;

	OutputFormatDR[PrintPre]
];


(*
	Calculates the 1-loop tadpole in the ultrasoft theory.
*)
TadPoleSS[]:=Module[{},
	If[verbose,Print["Calculating 1-Loop Tadpoles"]];
	fac=-1/(4 \[Pi]);
	ContriL1=-fac/2*Simplify[Table[Sum[muijL[[ii,ii]] lam3Cx[[i,ii,ii]],{ii,1,nSH}],{i,1,nSL}]];
	ContriLSE=-ZSij . TadPoleLight;

	TadPoleLightSS=TadPoleLight+ContriLSE-ContriL1;

	ContriH1=-fac/2*Simplify[Table[Sum[muijL[[ii,ii]] lam3CHeavy[[i,ii,ii]],{ii,1,nSH}],{i,1,nSH}]];

	TadPoleHeavySS=TadPoleHeavy-ContriH1;
];


(* ::Section:: *)
(*Help functions*)


(*
	The squared masses and the cubic and quartic couplings of all scalars, light and heavy, as one tensor each.
	The light fields have the indices 1,...,nSL, the heavy fields nSL+1,...,nSL+nSH, and the light masses are zero.
*)
AllScalarCouplingsSS[]:=Module[{nAll,masses,blockRules,allPermutations,lam3All,lam4All},
	nAll=nSL+nSH;
	masses=Join[Table[0,{nSL}],Table[muijL[[n,n]]^2,{n,1,nSH}]];
	blockRules[tensor_,offset_]:=Delete[SparseArray[tensor]//ArrayRules,-1]/.(x_List->v_):>(x+offset->v);
	allPermutations[rules_]:=DeleteDuplicates[Flatten[Table[Thread[Permutations[r[[1]]]->r[[2]]],{r,rules}]],#1[[1]]==#2[[1]]&];
	lam3All=SparseArray[allPermutations[Join[
		blockRules[lam3CLight,{0,0,0}],
		blockRules[lam3Cy,{0,0,nSL}],
		blockRules[lam3Cx,{0,nSL,nSL}],
		blockRules[lam3CHeavy,{nSL,nSL,nSL}]]],{nAll,nAll,nAll}];
	lam4All=SparseArray[allPermutations[Join[
		blockRules[lam4Light,{0,0,0,0}],
		blockRules[lamy,{0,0,0,nSL}],
		blockRules[lamK,{nSL,nSL,0,0}],
		blockRules[lamx,{0,nSL,nSL,nSL}],
		blockRules[lamKTotal,{nSL,nSL,nSL,nSL}]]],{nAll,nAll,nAll,nAll}];
	{masses,lam3All,lam4All}
];


(*
	One-loop cubic couplings of the fields in list1, list2, list3, from loops with at least one heavy line.
	It is the third derivative of 1/2 Tr ln(p^2+M^2) with respect to the fields: bubbles with one quartic and one
	cubic coupling, -1/2 lam4.lam3 B for each choice of the two legs at the quartic, and triangles with three cubic couplings, lam3.lam3.lam3 T.
*)
OneLoopCubicSS[list1_,list2_,list3_]:=Module[{masses,lam3All,lam4All,bubbleMat,triangle,nonzero,dLam},
	{masses,lam3All,lam4All}=AllScalarCouplingsSS[];
	bubbleMat=SparseArray[Table[BubbleIntegral3d[x,y],{x,masses},{y,masses}]];

	nonzero[tensor_]:=Delete[ArrayRules[tensor],-1];
	triangle[i_,j_,k_]:=Sum[
		If[c[[1,1]]==0,0,
			r[[2]]c[[2]]lam3All[[k,c[[1,1]],r[[1,1]]]]TriangleIntegral3d[masses[[r[[1,1]]]],masses[[r[[1,2]]]],masses[[c[[1,1]]]]]],
		{r,nonzero[lam3All[[i]]]},{c,Append[nonzero[lam3All[[j,r[[1,2]]]]],{0}->0]}];

	dLam[i_,j_,k_]:=-1/2*(
			Total[lam4All[[i,j]] lam3All[[k]] bubbleMat,2]
			+Total[lam4All[[i,k]] lam3All[[j]] bubbleMat,2]
			+Total[lam4All[[j,k]] lam3All[[i]] bubbleMat,2])
		+triangle[i,j,k];

	SparseArray[Table[dLam[i,j,k],{i,list1},{j,list2},{k,list3}]]//SimplifySparse
];


(*One-loop correction to the light-light-heavy cubic couplings lam3Cy*)
VertexCorrectionLLH[]:=OneLoopCubicSS[Range[nSL],Range[nSL],nSL+Range[nSH]];


(*
	One-loop light masses at first order in the tree-level tadpoles of the heavy fields, which shift the heavy fields by
	v0=-tadpole/M^2: the one-loop light-light-heavy vertex times v0, without the triangle with two light lines, which is the
	insertion of the light mass lam3Cy.v0 in ScalarMass2LoopSS; the shift of the one-loop tadpoles of the heavy fields by
	the one-loop heavy masses and the H^3 cubics, v01, at the H S^2 cubic; and the H^2 S^2 quartic with v0 and the shift v1
	by the one-loop tadpoles.
*)
TadpoleShiftMassSS[]:=Module[{M2,v0,v1,v01,masses,lam3All,lam4All,tadpoleVec,bubbleMat,massLoop,triangle,vertex},
	M2=Table[muijL[[h,h]]^2,{h,1,nSH}];
	v0=-Normal[TadPoleHeavy]/M2;
	If[nSL==0||AllTrue[v0,#===0&],Return[0]];
	v1=-Normal[TadPoleHeavySS-TadPoleHeavy]/M2;

	{masses,lam3All,lam4All}=AllScalarCouplingsSS[];
	tadpoleVec=Table[TadpoleIntegral3d[x],{x,masses}];
	bubbleMat=Table[BubbleIntegral3d[x,y],{x,masses},{y,masses}];
	massLoop=Table[1/2 Total[Diagonal[lam4All[[nSL+h,nSL+k]]] tadpoleVec]-1/2 Total[lam3All[[nSL+h]] lam3All[[nSL+k]] bubbleMat,2],{h,1,nSH},{k,1,nSH}];
	v01=-(massLoop . v0+Table[v0 . Normal[lam3CHeavy[[h]]] . v1,{h,1,nSH}])/M2;

	triangle=Table[Sum[lam3Cy[[i,a,b]]lam3Cy[[j,c,b]]lam3Cy[[a,c,h]]TriangleIntegral3d[0,0,M2[[b]]],{a,1,nSL},{c,1,nSL},{b,1,nSH}],{i,1,nSL},{j,1,nSL},{h,1,nSH}];
	vertex=Normal[VertexCorrectionLLH[]]-triangle;

	SparseArray[vertex . v0+Normal[lam3Cy] . v01+Table[v0 . Normal[lamK[[;;,;;,i,j]]] . v1,{i,1,nSL},{j,1,nSL}]]//SimplifySparse
];


(*Mixing between light and heavy fields from the shift of the heavy fields by their tadpoles, lam3Cx.(-tadpole/M^2)*)
TadpoleShiftLH[tadpole_]:=SparseArray[Table[-Sum[lam3Cx[[i,k,l]]tadpole[[l]]/muijL[[l,l]]^2,{l,1,nSH}],{i,1,nSL},{k,1,nSH}]];


(*
	Mixing between light and heavy fields, as a correction to the mass matrix at zero momentum: the one-loop
	diagrams with at least one heavy line, 1/2 lam4.A-1/2 lam3.lam3 B, and the shift of the heavy fields by their tadpoles.
*)
MixingLH[tadpole_]:=Module[{masses,lam3All,lam4All,tadpoleVec,bubbleMat,loop},
	{masses,lam3All,lam4All}=AllScalarCouplingsSS[];
	tadpoleVec=SparseArray[Table[TadpoleIntegral3d[x],{x,masses}]];
	bubbleMat=SparseArray[Table[BubbleIntegral3d[x,y],{x,masses},{y,masses}]];
	loop=Table[1/2 Total[Diagonal[lam4All[[i,nSL+k]]] tadpoleVec]-1/2 Total[lam3All[[i]] lam3All[[nSL+k]] bubbleMat,2],{i,1,nSL},{k,1,nSH}];
	SparseArray[loop]+TadpoleShiftLH[tadpole]//SimplifySparse
];


(*
	Contributions of a light-heavy mixing to the light quartic couplings, at first order in the mixing:
	the mixing converts one light leg into a heavy field, which ends at a H S^3 quartic or at a H^2 S cubic
	that is connected to a H S^2 cubic by the tree-level exchange of a heavy field.
*)
MixingQuarticSS[mix_]:=Module[{inv,help1,help2,result},
	inv=DiagonalMatrix[Table[1/muijL[[n,n]]^2,{n,1,nSH}]]//SparseArray;
	help1=Table[Sum[mix[[i,k]]inv[[k,k]]lamy[[j,l,m,k]],{k,1,nSH}],{i,1,nSL},{j,1,nSL},{l,1,nSL},{m,1,nSL}];
	help2=Table[Sum[mix[[i,k]]inv[[k,k]]lam3Cx[[j,k,kk]]inv[[kk,kk]]lam3Cy[[l,m,kk]],{k,1,nSH},{kk,1,nSH}],{i,1,nSL},{j,1,nSL},{l,1,nSL},{m,1,nSL}];
	result=-4 Symmetrize[help1,Symmetric[{1,2,3,4}]]+12 Symmetrize[help2,Symmetric[{1,2,3,4}]];
	result//SparseArray//SimplifySparse
];




(*
	Adds tree-level contributions before the matching.
*)
IdentifyTensorsPreSSDRalgo[]:=Module[{},

(*Quartic Tensor*)
	HelpList=DeleteDuplicates@SparseArray[Flatten@Simplify[ (lam4S)]]//Sort//FullSimplify;
	HelpVar=Table[lamPSS[a],{a,1,Delete[HelpList,1]//Length}];
	HelpVarMod=RelationsBVariables[HelpList,HelpVar];
	HelpSolveQuarticPreS=Table[{Delete[HelpList,1][[a]]->HelpVarMod[[a]]},{a,1,Delete[HelpList,1]//Length}]//Flatten;
	lam4S=(lam4S)//Normal//Simplify//FullSimplify//ReplaceAll[#,HelpSolveQuarticPreS]&//SparseArray;

	HelpList=DeleteDuplicates@SparseArray[Flatten@Simplify[ (lamK)]]//Sort//FullSimplify;
	HelpVar=Table[lamPSK[a],{a,1,Delete[HelpList,1]//Length}];
	HelpVarMod=RelationsBVariables[HelpList,HelpVar];
	HelpSolveQuarticPreK=Table[{Delete[HelpList,1][[a]]->HelpVarMod[[a]]},{a,1,Delete[HelpList,1]//Length}]//Flatten;
	lamK=(lamK)//Normal//Simplify//FullSimplify//ReplaceAll[#,HelpSolveQuarticPreK]&//SparseArray;

	HelpList=DeleteDuplicates@SparseArray[Flatten@Simplify[ (lamy)]]//Sort//FullSimplify;
	HelpVar=Table[lamPSY[a],{a,1,Delete[HelpList,1]//Length}];
	HelpVarMod=RelationsBVariables[HelpList,HelpVar];
	HelpSolveQuarticPreY=Table[{Delete[HelpList,1][[a]]->HelpVarMod[[a]]},{a,1,Delete[HelpList,1]//Length}]//Flatten;
	lamy=(lamy)//Normal//Simplify//FullSimplify//ReplaceAll[#,HelpSolveQuarticPreY]&//SparseArray;

	HelpList=DeleteDuplicates@SparseArray[Flatten@Simplify[ lamx]]//Sort//FullSimplify;
	HelpVar=Table[lamPSX[a],{a,1,Delete[HelpList,1]//Length}];
	HelpVarMod=RelationsBVariables[HelpList,HelpVar];
	HelpSolveQuarticPreX=Table[{Delete[HelpList,1][[a]]->HelpVarMod[[a]]},{a,1,Delete[HelpList,1]//Length}]//Flatten;
	lamx=(lamx)//Normal//Simplify//FullSimplify//ReplaceAll[#,HelpSolveQuarticPreX]&//SparseArray;
	
	HelpList=DeleteDuplicates@SparseArray[Flatten@Simplify[ lam4Tot]]//Sort//FullSimplify;
	HelpVar=Table[lamPST[a],{a,1,Delete[HelpList,1]//Length}];
	HelpVarMod=RelationsBVariables[HelpList,HelpVar];
	HelpSolveQuarticPreTot=Table[{Delete[HelpList,1][[a]]->HelpVarMod[[a]]},{a,1,Delete[HelpList,1]//Length}]//Flatten;
	(lam4Tot)//Normal//Simplify//FullSimplify//ReplaceAll[#,HelpSolveQuarticPreTot]&//SparseArray;

	IdentMatPre=List/@Join[HelpSolveQuarticPreS,HelpSolveQuarticPreK,HelpSolveQuarticPreY,HelpSolveQuarticPreX,HelpSolveQuarticPreTot]/.{b_->a_}:>a->b//Flatten[#,1]&;

];


(*
	Writes ultrasoft parameters in terms of soft ones.
*)
IdentifyTensorsSSDRalgo[]:=Module[{},

If[verbose,Print["Identifying Components"]];

	If[mode>=1,

(*Quartic Tensor*)
		HelpList=SparseArray[lam4S+lam3DSS]["NonzeroValues"]//Simplify//DeleteDuplicates//Sort;
		If[HelpList[[1]]==0&&Length[HelpList]>1,
			HelpList=Delete[HelpList,1];
		];
		HelpVar=Table[lamSS[a],{a,1,HelpList//Length}];
		HelpVarMod=RelationsBVariables[HelpList,HelpVar];
		HelpSolveQuarticS=Table[{HelpList[[a]]->HelpVarMod[[a]]},{a,1,HelpList//Length}]//Flatten;
		lam4US=(lam4S+lam3DSS)//Normal//Simplify//ReplaceAll[#,HelpSolveQuarticS]&//SparseArray;
	
		If[Length[Variables[lam4US["NonzeroValues"]/.lamSS[x_]->0]]>0,
			HelpSolveQuarticS2=HelpSolveQuarticS/.lamSS[x_]->lamSSS[x];
			HelpList=DeleteDuplicates[lam4US["NonzeroValues"]//Simplify]/.lamSS[x_]->lamSSS[x]//Sort;
			If[HelpList[[1]]==0&&Length[HelpList]>1,
				HelpList=Delete[HelpList,1];
			];
			HelpVar=Table[lamSS[a],{a,1,HelpList//Length}];
			HelpVarMod=RelationsBVariables[HelpList,HelpVar];
			HelpSolveQuarticS=Table[{HelpList[[a]]->HelpVarMod[[a]]},{a,1,HelpList//Length}]//Flatten;
			lam4US=Normal[lam4US]/.lamSS[x_]->lamSSS[x]//Simplify//ReplaceAll[#,HelpSolveQuarticS]&//SparseArray;
			
			HelpSolveQuarticS=HelpSolveQuarticS//ReplaceAll[#,HelpSolveQuarticS2/.(b_->a_):>a->b]&
		];
	
	
(*Cubic Tensor*)
		HelpList=DeleteDuplicates@SparseArray[Flatten@Simplify[(lam3CSSS+lam3CLight)]]//Sort//FullSimplify;
		If[HelpList[[1]]==0&&Length[HelpList]>1,
			HelpList=Delete[HelpList,1];
		];
		HelpVar=Table[cSSSS[a],{a,1,HelpList//Length}];
		HelpVarMod=RelationsBVariables[HelpList,HelpVar];
		HelpSolveCubicSSS=Table[{HelpList[[a]]->HelpVarMod[[a]]},{a,1,HelpList//Length}]//Flatten//Simplify;
		lam3CSRedSS=(lam3CSSS+lam3CLight)//Normal//Simplify//FullSimplify//ReplaceAll[#,HelpSolveCubicSSS]&//SparseArray;
(*Scalar-Vector Tensor*)

		HelpList=DeleteDuplicates@Simplify@Flatten[(HabijVL+GvvssTSS)]//Sort;
		If[HelpList[[1]]==0&&Length[HelpList]>1,
			HelpList=Delete[HelpList,1];
		];
		HelpVar=Table[ lamVTSS[a],{a,1,HelpList//Length}];
		HelpVarMod=RelationsBVariables[HelpList,HelpVar];
		HelpSolveVecTS=Table[{HelpList[[a]]->HelpVarMod[[a]]},{a,1,HelpList//Length}]//Flatten;
		lamKVecTSS= (HabijVL+GvvssTSS)//Normal//Simplify//ReplaceAll[#,HelpSolveVecTS]&//SparseArray;		
	];
(*Scalar Mass*)
	If[mode>=2,
		HelpList=DeleteDuplicates@Flatten@Simplify[ xLO muijLight+xLO muijSSLO+xNLO muijSSNLO]//Sort;
		If[HelpList[[1]]==0&&Length[HelpList]>1,
				HelpList=Delete[HelpList,1];
			];
		HelpVar=Table[ muijSS[a],{a,1,HelpList//Length}];
		HelpVarMod=RelationsBVariables[HelpList,HelpVar];
		HelpSolveMassS=Table[{HelpList[[a]]->HelpVarMod[[a]]},{a,1,HelpList//Length}]//Simplify//Flatten;
		muijSNLOSS=xLO muijLight+xLO muijSSLO+xNLO muijSSNLO//Normal//Simplify//ReplaceAll[#,HelpSolveMassS]&//SparseArray;
	,
		HelpList=DeleteDuplicates@Flatten@Simplify[ xLO muijLight+xLO muijSSLO]//Sort;
		If[HelpList[[1]]==0&&Length[HelpList]>1,
			HelpList=Delete[HelpList,1];
		];
		HelpVar=Table[ muijSS[a],{a,1,HelpList//Length}];
		HelpVarMod=RelationsBVariables[HelpList,HelpVar];
		HelpSolveMassS=Table[{HelpList[[a]]->HelpVarMod[[a]]},{a,1,HelpList//Length}]//Simplify//Flatten;
		muijSNLOSS=xLO muijLight+xLO muijSSLO//Normal//Simplify//ReplaceAll[#,HelpSolveMassS]&//SparseArray;
	];

(*Scalar Tadpoles*)
	HelpList=DeleteDuplicates@Flatten@Simplify[ TadPoleLightSS]//Sort;
	If[HelpList[[1]]==0&&Length[HelpList]>1,
		HelpList=Delete[HelpList,1];
	];
	HelpVar=Table[ dSS[a],{a,1,HelpList//Length}];
	HelpVarMod=RelationsBVariables[HelpList,HelpVar];
	HelpSolveTadpoleSS=Table[{HelpList[[a]]->HelpVarMod[[a]]},{a,1,HelpList//Length}]//Flatten;
	TadPoleSSSLO=TadPoleLightSS//Normal//Simplify//ReplaceAll[#,HelpSolveTadpoleSS]&//SparseArray;
	
	If[mode>=1,
		IdentMatSS=List/@Join[HelpSolveQuarticS,HelpSolveVecTS,HelpSolveMassS,HelpSolveCubicSSS,HelpSolveTadpoleSS]/.{b_->a_}:>a->b//Flatten[#,1]&;
	,
		IdentMatSS=List/@Join[HelpSolveMassS,HelpSolveTadpoleSS]/.{b_->a_}:>a->b//Flatten[#,1]&;
	];

];


(*
	Creates auxiliary tensors that appear in the matching.
*)
CreateHelpTensorsSS[]:=Module[{},
If[verbose,Print["Creating Help Tensors"]];

	HabijL=Transpose[Activate @ TensorContract[Inactive[TensorProduct][gvssL,gvssL], { {3, 5}}],{1,3,2,4}]//SimplifySparse;
	HabijVL=HabijL+Transpose[HabijL,{2,1,3,4}]//SimplifySparse;
	SelfEnergySS=  Inactivate[TensorProduct[gAvss,gAvss]]//SimplifySparse;
	HabijA=Transpose[Activate@TensorContract[SelfEnergySS,{{3,5}}],{1,3,2,4}]//SimplifySparse;
	HabijVA=HabijA+Transpose[HabijA,{1,2,4,3}]//SimplifySparse;
];


{TadPoleSSSLO};


{DiaRules,RotRules};(*Diagonalization matrix and masses*)


(*
	The tensor with all indices shifted by offset, in an array of dimensions dims.
	The rules are replaced one by one, so that the list of rules itself is never matched, and dims is
	given explicitly, since the largest index of the nonzero entries can be smaller.
*)
ShiftedTensorSS[tensor_,offset_,dims_]:=SparseArray[Delete[ArrayRules[SparseArray[tensor]],-1]/.(x_List->v_):>(x+offset->v),dims];


(*
	Prepares the reduction to the ultrasoft scale. Here scalars are divided into heavy and soft depending on the user's input.
	Heavy scalars are then grouped with temporal scalars and integrated out.
*)
PrepareSoftToSuperSoft[ListHardI_]:=Module[{ListP=ListHardI},

(*Fix for Non-diagonal Debye masses*)
(*************************)

(*
	Here we take the "Damn the torpedo, full speed ahead" approach.
*)

	HeavyScalars=Transpose[List@ListP]; (*The list of all heavy scalars*)

(*
	Renames 4d couplings to 3d ones
*)

	VarQuartic=Join[lam4//Normal//Variables]//DeleteDuplicates;
	SubQuartic=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarQuartic}];
	lam3DSp= lam4//Normal//ReplaceAll[#,SubQuartic]&//SparseArray;

	VarCubic=Join[lam3//Normal//Variables]//DeleteDuplicates;
	SubCubic=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarCubic}];
	lam3CSRedUS=lam3//Normal//ReplaceAll[#,SubCubic]&//SparseArray;

	VarMass=Join[muij//Normal//Variables]//DeleteDuplicates;
	SubMass=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarMass}];
	muijSNLOP=muij//Normal//ReplaceAll[#,SubMass]&//SparseArray;



	If[Length[HeavyScalars]<1,HScal=False,HScal=True];

	TotScalar=Table[{i},{i,1,ns}];
	LightScalar=Delete[TotScalar,HeavyScalars];
	nSH=nv+(HeavyScalars//Length); (*Adds heavy scalars to the temporal scalars*)
	nSL=ns-(HeavyScalars//Length); (*Removes heavy scalars from the scalars*)

	muijLS=Sqrt[muijSNLOP]//SparseArray; (*Rewrites squared masses as masses. Just a trick to make things easier.*)
	muijLS[[LightScalar[[;;,1]],LightScalar[[;;,1]]]]=0;

	VarGauge=GaugeCouplingNames;
	SubGauge=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarGauge}];
	gvvvSS=gvvv//Normal//ReplaceAll[#,SubGauge]&//SparseArray;


	If[HScal,
(*Quartic Tensor*)
(*
	Here we create couplings between light and heavy scalarws.
*)
			lam4Light=Table[lam3DSp[[a,b,c,d]],{a,LightScalar[[;;,1]]},{b,LightScalar[[;;,1]]},{c,LightScalar[[;;,1]]},{d,LightScalar[[;;,1]]}]//SparseArray; (*Light-scalar part of scalar quartic*)
			lam4Heavy=Table[lam3DSp[[a,b,c,d]],{a,HeavyScalars[[;;,1]]},{b,HeavyScalars[[;;,1]]},{c,HeavyScalars[[;;,1]]},{d,HeavyScalars[[;;,1]]}]//SparseArray; (*Heavy-scalar part of scalar quartic*)
			lamKHeavy=Table[lamKVec[[a,b,c,d]],{a,1,nv},{b,1,nv},{c,HeavyScalars[[;;,1]]},{d,HeavyScalars[[;;,1]]}]//SparseArray;  (*Heavy-scalar part of the scalar-temporalScal couplings*)
			lam4Tot=lam3DSp//SparseArray;


			A1=Delete[ShiftedTensorSS[lam4Heavy,{nv,nv,nv,nv},{nSH,nSH,nSH,nSH}]//ArrayRules,-1];
			A2=Delete[ShiftedTensorSS[lamKHeavy,{0,0,nv,nv},{nSH,nSH,nSH,nSH}]//ArrayRules,-1];
			A3=Delete[lamAAS//ArrayRules,-1];
			lamKTotal=SymmetrizedArray[Join[A1,A2,A3],{nSH,nSH,nSH,nSH},Symmetric[{1,2,3,4}]]//SparseArray; (*Total heavy-scalar tensors: Heavy scalars+Temporal scalars*)

(*TadPoles*)
			VarTadpole=Join[lam1//Normal//Variables]//DeleteDuplicates;
			SubTadpole=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarTadpole}];

			TadPoleLight=Table[lam1[[a]],{a,LightScalar[[;;,1]]}]//ReplaceAll[#,SubTadpole]&//SparseArray;
			TadPoleTemp=Table[lam1[[a]],{a,HeavyScalars[[;;,1]]}]//ReplaceAll[#,SubTadpole]&//SparseArray;


			TadPoleHeavy=ShiftedTensorSS[TadPoleTemp,{nv},{nSH}];


			lamx1=Table[lam3DSp[[a,b,c,d]],{a,LightScalar[[;;,1]]},{b,HeavyScalars[[;;,1]]},{c,HeavyScalars[[;;,1]]},{d,HeavyScalars[[;;,1]]}]//SparseArray;
			lamy1=Table[lam3DSp[[a,b,c,d]],{a,LightScalar[[;;,1]]},{b,LightScalar[[;;,1]]},{c,LightScalar[[;;,1]]},{d,HeavyScalars[[;;,1]]}]//SparseArray;


			lamx=ShiftedTensorSS[lamx1,{0,nv,nv,nv},{nSL,nSH,nSH,nSH}];
			lamy=ShiftedTensorSS[lamy1,{0,0,0,nv},{nSL,nSL,nSL,nSH}];


(*Mixed Cubic tensors*)
			lam3CLight=Table[lam3CSRedUS[[a,b,c]],{a,LightScalar[[;;,1]]},{b,LightScalar[[;;,1]]},{c,LightScalar[[;;,1]]}]//SparseArray;
			lam3CHeavy1=Table[lam3CSRedUS[[a,b,c]],{a,HeavyScalars[[;;,1]]},{b,HeavyScalars[[;;,1]]},{c,HeavyScalars[[;;,1]]}]//SparseArray;
			lam3CTot=lam3CSRedUS//SparseArray;
			
			lam3CHeavy=ShiftedTensorSS[lam3CHeavy1,{nv,nv,nv},{nSH,nSH,nSH}];


			lam3Cx1=Table[lam3CSRedUS[[a,b,c]],{a,LightScalar[[;;,1]]},{b,HeavyScalars[[;;,1]]},{c,HeavyScalars[[;;,1]]}]//SparseArray;
			lam3Cy1=Table[lam3CSRedUS[[a,b,c]],{a,LightScalar[[;;,1]]},{b,LightScalar[[;;,1]]},{c,HeavyScalars[[;;,1]]}]//SparseArray;


			lam3Cx=ShiftedTensorSS[lam3Cx1,{0,nv,nv},{nSL,nSH,nSH}];
			lam3Cy=ShiftedTensorSS[lam3Cy1,{0,0,nv},{nSL,nSL,nSH}];

(*Mixed Hard-Soft mass matrix*)
			muijMix1=Table[muij[[a,b]],{a,LightScalar[[;;,1]]},{b,HeavyScalars[[;;,1]]}]//SparseArray;

			muijMix=ShiftedTensorSS[muijMix1,{0,nv},{nSL,nSH}];


(*Scalar-Temporal-Vector Tensor*)
		lamKVLight=Table[lamKVec[[a,b,c,d]],{a,1,nv},{b,1,nv},{c,LightScalar[[;;,1]]},{d,LightScalar[[;;,1]]}]//SparseArray;
		lamKVHeavy=Table[lam3DSp[[a,b,c,d]],{a,HeavyScalars[[;;,1]]},{b,HeavyScalars[[;;,1]]},{c,LightScalar[[;;,1]]},{d,LightScalar[[;;,1]]}]//SparseArray;
		lamkVTot=ShiftedTensorSS[lamKVLight,{0,0,0,0},{nSH,nSH,nSL,nSL}]+ShiftedTensorSS[lamKVHeavy,{nv,nv,0,0},{nSH,nSH,nSL,nSL}];


		VarGauge=GaugeCouplingNames;
		SubGauge=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarGauge}];

(*Heavy-Scalar-Vector Tensor*)
		gvssHeavy=Table[gvss[[a,b,c]],{a,1,nv},{b,HeavyScalars[[;;,1]]},{c,HeavyScalars[[;;,1]]}]//ReplaceAll[#,SubGauge]&//SparseArray;
		gvssVec=gvvv//Normal//ReplaceAll[#,SubGauge]&//SparseArray;
		gvssVTot=Table[ArrayFlatten[{{gvssVec[[a]],0},{0,gvssHeavy[[a]]}}],{a,1,Length[gvssVec]}]//SparseArray;

(*Light-Scalar-Vector Tensor*)
		gvssL=Table[gvss[[a,c,d]],{a,1,nv},{c,LightScalar[[;;,1]]},{d,LightScalar[[;;,1]]}]//Normal//ReplaceAll[#,SubGauge]&//SparseArray;

(*Heavy-Mass matrix*)
		MassHeavy=Table[muijSNLOP[[a,b]],{a,HeavyScalars[[;;,1]]},{b,HeavyScalars[[;;,1]]}]//SparseArray;
		muijL=ArrayFlatten[{{Sqrt[muabDef],0},{0,Sqrt[MassHeavy]}}]//SparseArray;

(*Light-Mass matrix*)
		MassLight=Table[muijSNLOP[[a,b]],{a,LightScalar[[;;,1]]},{b,LightScalar[[;;,1]]}]//SparseArray;
		muijLight=ArrayFlatten[MassLight]//SparseArray;

		lamK=lamkVTot//SparseArray;
		lam4S=lam4Light//SparseArray;
		lam4K=lamKTotal//SparseArray;
		gAvss=gvssVTot//SparseArray;
		OutputFormatDR[Join[lamK,lam4S,lam4K,lamx,lamy,gAvss,gvssL,muijL]];
	,
		lam3CTot=lam3CSRedUS//SparseArray;
		A3=Delete[lamAAS//ArrayRules,-1];
		lamKTotal=SymmetrizedArray[A3,{nSH,nSH,nSH,nSH},Symmetric[{1,2,3,4}]]//SparseArray;
		lam3Cx=SymmetrizedArray[{{1,1,1}->0},{nSL,nSH,nSH},Symmetric[{2,3}]]//SparseArray;
		lam3Cx=Table[lamvvsLS[[b,c,a]],{a,LightScalar[[;;,1]]},{b,1,nv},{c,nv}]//SparseArray;
		lam3CHeavy=SymmetrizedArray[{{1,1,1}->0},{nSH,nSH,nSH},Symmetric[{1,2,3}]]//SparseArray;
		lam3Cy=SymmetrizedArray[{{1,1,1}->0},{nSL,nSL,nSH},Symmetric[{1,2}]]//SparseArray;
		lamx=SymmetrizedArray[{{1,1,1,1}->0},{nSL,nSH,nSH,nSH},Symmetric[{2,3,4}]]//SparseArray;
		lamy=SymmetrizedArray[{{1,1,1,1}->0},{nSL,nSL,nSL,nSH},Symmetric[{1,2,3}]]//SparseArray;
		muijMix=SymmetrizedArray[{{1,1}->0},{nSL,nSH},Symmetric[{1}]]//SparseArray;
		lam3CLight=Table[lam3CSRedUS[[a,b,c]],{a,LightScalar[[;;,1]]},{b,LightScalar[[;;,1]]},{c,LightScalar[[;;,1]]}]//SparseArray;
		lam4Light=Table[lam3DSp[[a,b,c,d]],{a,LightScalar[[;;,1]]},{b,LightScalar[[;;,1]]},{c,LightScalar[[;;,1]]},{d,LightScalar[[;;,1]]}]//SparseArray;

(*TadPoles*)
		VarTadpole=Join[lam1//Normal//Variables]//DeleteDuplicates;
		SubTadpole=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarTadpole}];

		TadPoleLight=Table[lam1[[a]],{a,LightScalar[[;;,1]]}]//ReplaceAll[#,SubTadpole]&//SparseArray;
		TadPoleHeavy=SymmetrizedArray[{{1}->0},{nSH},Symmetric[{1}]]//SparseArray;

		lamKVLight=Table[lamKVec[[a,b,c,d]],{a,1,nv},{b,1,nv},{c,LightScalar[[;;,1]]},{d,LightScalar[[;;,1]]}]//SparseArray;
		lamkVTot=lamKVLight//SparseArray;

		VarGauge=GaugeCouplingNames;
		SubGauge=Table[c->Symbol[ToString[c]<>ToString["3d"]],{c,VarGauge}];

		gvssVec=gvvv//Normal//ReplaceAll[#,SubGauge]&//SparseArray;
		gvssVTot=gvssVec//SparseArray; 
		gvssL=Table[gvss[[a,c,d]],{a,1,nv},{c,LightScalar[[;;,1]]},{d,LightScalar[[;;,1]]}]//Normal//ReplaceAll[#,SubGauge]&//SparseArray;

		muijL=Sqrt[muabDef]//Normal//SparseArray;
		muijLight=muijSNLOP//Normal//SparseArray;

		lamK=lamkVTot//SparseArray;
		lam4S=lam4Light//SparseArray;
		lam4K=lamKTotal//SparseArray;
		gAvss=gvssVTot//SparseArray;

	];


];


{Contri\[Beta]UltraSoftToSoft,ContriTadPoleUltraSoftToSoft};


DiagonalTensor2[s_SparseArray,a_Integer,b_Integer] := With[
(*
	Takes a A_i1,...ia,...ib,... and creates_ia,i1,...,...; so the ia and ib indices are collapsed to one diagonal index.
*)    
    {
    s1=Flatten[s,{{a},{b}}]
    },
     Table[i,{i,s1}]//Table[#[[i,i]],{i,1,Length[#]}]&//SparseArray
    ]


DiagonalTensor[s_SparseArray,a_Integer,b_Integer]:= Module[{},
(*
	Takes a A_i1,...ia,...ib,... and creates_ia,i1,...,...; so the ia and ib indices are collapsed to one diagonal index.
*)   
If[a==1&&b==2,
    PermTens=Table[i,{i,s}]//Table[#[[i,i]],{i,1,Length[#]}]&//SparseArray
,
If[a==1,
   PermTens=Transpose[s,b<->2];
    PermTens=Table[i,{i,PermTens}]//Table[#[[i,i]],{i,1,Length[#]}]&//SparseArray;
   Transpose[PermTens,(b-1)<->1]
,
	   PermTens=Transpose[s,a<->1]//Transpose[#,b<->2]&;
    PermTens=Table[i,{i,PermTens}]//Table[#[[i,i]],{i,1,Length[#]}]&//SparseArray;
   Transpose[PermTens,(a-1)<->1]//Transpose[#,(b-1)<->2]&
]
]
    ]


(* ::Section:: *)
(*Printing functions*)


(*
	Prints scalar masses.
*)
PrintScalarMassUS[optP_]:=Module[{opt=optP},
If[verbose,Print["Printing Scalar Masses"]];

	VarGauge=Join[muijLight//Normal//Variables]//DeleteDuplicates;
	SubGauge=Table[c->Symbol[ToString[c]<>ToString["US"]],{c,VarGauge}];

	muijp=muijLight//Normal//ReplaceAll[#,SubGauge]&;
	var=Normal[muijp]//Variables;
	helpMass=Normal[muijp-muijSNLOSS];
	ResScalp=Reduce[helpMass==0,var]//ToRules[#]&;
	SolveTemp=var/.ResScalp;
	SolMassPre=Table[{var[[i]]->SolveTemp[[i]]},{i,1,Length@var}]//Flatten[#,1]&//ReplaceAll[#,IdentMatSS]&;

	If[opt=="LO",
		SolMass=SolMassPre/.xLO->1/.xNLO->0;
	,
		SolMass=SolMassPre/.xLO->0/.xNLO->1;
	];
(*Printing Result*)
	OutputFormatDR[Join[SolMass]]
];

PrintScalarMassUS[]:=Module[{},
If[verbose,Print["Printing Scalar Masses"]];

	VarGauge=Join[muijLight//Normal//Variables]//DeleteDuplicates;
	SubGauge=Table[c->Symbol[ToString[c]<>ToString["US"]],{c,VarGauge}];

	muijp=muijLight//Normal//ReplaceAll[#,SubGauge]&;
	var=Normal[muijp]//Variables;
	helpMass=Normal[muijp-muijSNLOSS];
	ResScalp=Reduce[helpMass==0,var]//ToRules[#]&;
	SolveTemp=var/.ResScalp;
	SolMassPre=Table[{var[[i]]->SolveTemp[[i]]},{i,1,Length@var}]//Flatten[#,1]&//ReplaceAll[#,IdentMatSS]&;
	
	SolMass=SolMassPre/.xLO->1/.xNLO->1;
	
(*Printing Result*)
	OutputFormatDR[Join[SolMass]]
];


(*
	Prints all the couplings in the supersoft theory.	
*)

PrintCouplingsUS[]:=Module[{},

(*Non-Abelian Couplings*)

	NonAbelianCouplingSS[]; (*Calculates non-abelian couplings*)
	
	If[mode==0,GgvvvSS=EmptyArray[{nv,nv,nv}]];
	GabcdTemp=GgvvvSS . gvvvSS+gvvvSS . GgvvvSS;
	Temp=Activate @ TensorContract[Inactive[TensorProduct][gvvvSS,gvvvSS], {{3, 6}}]; 
	GabVTree=TensorContract[Temp,{{2,3}}]//Normal;
	GabVLoop=TensorContract[GabcdTemp,{{2,3}}]//Normal;

	HelpList=DeleteDuplicates@Flatten@FullSimplify[GabVTree+GabVLoop ]//Sort; (*Adds the tree-level result*)
	HelpVar=Table[ lamVNASS[a],{a,1,Delete[HelpList,1]//Length}];
	HelpVarMod=RelationsBVariables[HelpList,HelpVar]; (*Finds a minimal basis of couplings*)
	HelpSolveNASS=Table[{Delete[HelpList,1][[a]]->HelpVarMod[[a]]},{a,1,Delete[HelpList,1]//Length}]//Flatten;
	lamVecNASS=GabVTree+GabVLoop//Normal//FullSimplify//ReplaceAll[#,HelpSolveNASS]&//SparseArray;
	IdentMatNASS=List/@HelpSolveNASS/.{b_->a_}:>a->b//Flatten[#,1]&;

	VarGauge=Table[Symbol[ToString[c]<>ToString["3d"]],{c,GaugeCouplingNames}];
	SubGauge=Table[c->Symbol[ToString[c]<>ToString["US"]],{c,VarGauge}];


	A1=lamVecNASS//Normal;
	A2=GabVTree//Normal//ReplaceAll[#,SubGauge]&;
	Var3D=A2//Variables;
	RepVar3D=#->Sqrt[#]&/@Var3D;
	A2Mod=A2/.RepVar3D;
	Sol1=Solve[A2Mod==A1,Var3D]/.IdentMatNASS//Flatten[#,1]&//FullSimplify; (*Solves ultrasoft couplings in terms of soft ones*)
	ResGaugeNASS=Table[List[Sol1[[c]]]/.{b_->a_}:>b^2->a,{c,1,Length[Sol1]}]//Simplify;

(* Gauge couplings*)
	If[mode==0,lamKVecTSS=HabijVL];
	A1=TensorContract[lamKVecTSS,{{3,4}}]//Normal;
	A2=TensorContract[HabijVL,{{3,4}}]//Normal//ReplaceAll[#,SubGauge]&;
(*Trick to avoid problems when kinetic mixing*)
	A1=DiagonalMatrix[Diagonal[A1]];
	A2=DiagonalMatrix[Diagonal[A2]];


	Var3D=VarGauge//ReplaceAll[#,SubGauge]&//Variables;
	RepVar3D=#->Sqrt[#]&/@Var3D;
	A2Mod=A2/.RepVar3D;
	Sol1=Solve[A2Mod==A1,Var3D]/.IdentMatSS//Flatten[#,1]&;
	ResGauge=Table[List[Sol1[[c]]]/.{b_->a_}:>b^2->a,{c,1,Length[Sol1]}];

(* Scalar Quartics*)
	If[mode==0,lam4US=lam4S];
	QuarticVar=lam4S//Normal//Variables;
	SubGauge=Table[c->Symbol[ToString[c]<>ToString["US"]],{c,QuarticVar}];
	NonZeroPos=SparseArray[lam4S]["NonzeroPositions"];
	SolVar=Extract[lam4S-lam4US,{#}]&/@NonZeroPos//DeleteDuplicates;
	ResScalp=Reduce[SolVar==0,QuarticVar]//ToRules[#]&;
	SolveTemp=QuarticVar/.ResScalp;
	ResScal=Table[{ReplaceAll[QuarticVar[[i]],SubGauge]->SolveTemp[[i]]},{i,1,Length@QuarticVar}]//Flatten[#,1]&//ReplaceAll[#,IdentMatSS]&;


(* Scalar Cubics*)
	If[mode==0,lam3CSRedSS=lam3CLight];
	VarGauge=Join[lam3CLight//Normal//Variables]//DeleteDuplicates;
	SubGauge=Table[c->Symbol[ToString[c]<>ToString["US"]],{c,VarGauge}];
	lam3p=lam3CLight//Normal//ReplaceAll[#,SubGauge]&;
	SolVar=lam3CSRedSS-lam3p//Normal;
	CubicVar=lam3p//Normal//Variables;
	ResCubicSS=Solve[SolVar==0,CubicVar]/.IdentMatSS//Flatten[#,1]&;

(*Printing Result*)
	PrintPre=Join[ResScal,ResGauge,ResGaugeNASS,ResCubicSS]//Normal//FullSimplify//DeleteDuplicates;
	PrintPre=DeleteCases[PrintPre,0->0];

	OutputFormatDR[PrintPre]
];


(*
	Prints tadpoles.
*)
PrintTadpolesUS[optP_]:=Module[{opt=optP},
If[verbose,Print["Printing Scalar Masses"]];

	VarGauge=Join[TadPoleLight//Normal//Variables]//DeleteDuplicates;
	SubGauge=Table[c->Symbol[ToString[c]<>ToString["US"]],{c,VarGauge}];

	lam1p=TadPoleLight//Normal//ReplaceAll[#,SubGauge]&;
	var=Normal[lam1p]//Variables;
	helpMass=Normal[lam1p-TadPoleSSSLO];
	ResScalp=Reduce[helpMass==0,var]//ToRules[#]&;
	SolveTemp=var/.ResScalp;
	SolTadpole=Table[{var[[i]]->SolveTemp[[i]]},{i,1,Length@var}]//Flatten[#,1]&//ReplaceAll[#,IdentMatSS]&;

	OutputFormatDR[Join[SolTadpole]]
];


PrintIdentificationUS[]:=Module[{},
	OutputFormatDR[Join[IdentMatSS]]
];


errPrintTensUS="Order of Tensors: (1) Scalar Quartic, (2) Vector-Scalar Couplings, (3) Scalar Mass";


PrintTensorUSDRalgo[]:="nothing"/;Print[errPrintTensUS];


PrintTensorUSDRalgo[ind_]:=Module[{},	
	
	Switch[ind, 
    1, Return[OutputFormatDR[lam4US]], 
    2, Return[OutputFormatDR[lamKVecTSS]], 
    3, Return[OutputFormatDR[muijSNLOSS]], 
    _, 
    Print[errPrintTensUS];
  ];
      Return[""]

];
