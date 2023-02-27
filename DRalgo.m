(* ::Package:: *)

(* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ *)

(* :Title: DRalgo                                                          	*)

(*
       This software is covered by the GNU General Public License 3.
       Copyright (C) 2021-2022 Andreas Ekstedt
       Copyright (C) 2021-2022 Philipp Schicho
       Copyright (C) 2021-2022 Tuomas V.I. Tenkanen

*)

(* :Summary:	DRalgo is an algorithm that constructs
				the effective high-temperature field theory for generic models.	*)	

(* ------------------------------------------------------------------------ *)


BeginPackage["DRalgo`"] (* digital history: at its early days in 2021, DRalgo was called FireStorm *)


Unprotect@Definition;
Definition[x_Symbol] /; StringMatchQ[Context[x], "Package`" ~~ ___] := 
    StringReplace[ToString@FullDefinition[x], 
        (WordCharacter .. ~~ DigitCharacter ... ~~ "`") .. ~~ s_ :> s
    ];
Protect@Definition;


(*
	Welcome banner: All credit for this part to GroupMath
*)
TexFor[text_]:=Style[text,{GrayLevel[0.3]}]
result={};
AppendTo[result,Row[{
	TexFor["DRDRDRDRDRDRDRDRDRDRDRDRDRDR "],
	TexFor["DRalgo"],
	TexFor[" DRDRDRDRDRDRDRDRDRDRDRDRDRDRD"]}]];
AppendTo[result,Row[{"Version: "//TexFor,"1.02 beta (16-05-2022)"//TexFor}]];
AppendTo[result,Row[{"Authors: "//TexFor,"Andreas Ekstedt, Philipp Schicho, Tuomas V.I. Tenkanen"//TexFor}]];
AppendTo[result,Row[{"Reference: "//TexFor,"2205.08815 [hep-ph]"//TexFor}]];
AppendTo[result,Row[{"Repository link: "//TexFor,
	Hyperlink[Mouseover[TexFor["github.com/DR-algo/DRalgo"],Style["github.com/DR-algo/DRalgo",Bold]],
	"https://github.com/DR-algo/DRalgo"]}]];
(* DRalgoLoad=Import[FileNameJoin[{DirectoryName[$InputFileName],"logo.eps"}],ImageSize->{200.,300.}]; *)
DRalgoLoad=Image[CompressedData["1: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"],"Byte",ColorSpace->"RGB",Interleaving->True];

AppendTo[result,Style["DRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRD",{GrayLevel[0.3]}]];
Print[Grid[{{Image[DRalgoLoad,ImageSize->140],Row[result,"\n",BaseStyle->(FontFamily->"Consolas")]}},Alignment->{Left,Center}]]
(**)

(*
	List of public functions.
*)
ImportModelDRalgo::usage="Loads the model and creates help tensors"
DefineVectorMass::usuage="Defines Custom Vector-Boson masses"
PrintDebyeMass::usuage="Prints Debye Masses in the soft theory"
PrintTemporalScalarCouplings::usage="Prints Temporal-Vector Couplings in the soft theory"
PerformDRhard::usage="Performs the dimensional reduction from hard to soft"
PrintIdentification::usage="Shows dim-red parameters in terms of 4-d ones"
PrintTensorDRalgo::usage="Shows dim-red parameters in terms of 4-d ones"
PrintCouplings::usage="Prints all couplings in the soft theory"
PrintScalarMass::usage="Prints scalar masses in the soft theory"
BetaFunctions3DUS::usage="Prints beta functions in the supersoft theory"
PrintPressure::usage="Prints the preassure in the soft theory"
PrintEffectivePotential::usage="Prints the effective potential in the supersoft theory"
PrintScalarMassUS::usage="Prints scalar masses in the supersoft theory"
AnomDim4D::usage="Prints anomalous dimensions in the 4d theory"
BetaFunctions4D::usage="Prints beta-functions in the 4d theory"
DefineNewTensorsUS::usage="Defines Custom Tensors for Effective-Potential Computations"
SymmetricPhaseEnergy::usage="Prints the energy of the Symmetric Phase"
PrintCouplingsUS::usage="Prints couplings in the supersoft theory"
PrintIdentificationUS::usage="Shows supersoft parameters in terms of soft parameters"
PrintTensorUSDRalgo::usage="Prints tensors in the supersoft theory"
PrintConstants::usage="Prints constants used in the matching"
PrintGenericBasis::usage="Rewrites the results in a different basis"
PrintTensorsVEV::usage="Prints background-field dependent couplings and masses"
DefineVEVS::usuage="Defines background fields for scalars"
RotateTensorsUSPostVEV::usuage="Rotate the field basis. Both for scalar and vector."
RotateTensorsCustomMass::usuage="Creates custom field-dependent masses"
DefineGroup::usuage="Loads the group and names Debye masses"
PerformDRsoft::usage="Perform the reduction from soft to supersoft"
CalculatePotentialUS::usage="Calculates the effective potential"
CompareInvariants::usuage="Finds relations between couplings by calculating basis-invariant tensors"
AllocateTensors::usage="Creates gauge generators"
GradQuartic::usage="Creates Quartic tensors"
GradCubic::usage="Creates Cubic tensors"
GradTadpole::usage="Creates Tadpole tensors"
GradSextic::usage="Creates dim 6 tensors"
GradMass::usage="Creates Mass tensors"
CreateInvariant::usage="Creates an invariant"
CreateInvariantYukawa::usage="Creates Yukawa Tensor"
GradYukawa::usuage="Creates Yukawa tensor"
PrintTadpoles::usuage="Prints Tadpoles"
PrintGaugeRepPositions::usuage="Prints the indices of Gauge reps"
PrintScalarRepPositions::usuage="Prints the indices of Scalar reps"
PrintFermionRepPositions::usuage="Prints the indices of Fermion reps"
DefineNF::usuage="Allows the user to add an arbitrary number of fermion families"
PrintTadpolesUS::usuage="Prints tadpoles in the supersoft theory"
GradMassFermion::usuage="Creates Fermion Invariants"
CreateInvariantFermion::usuage="Creates Fermion Invariants"
SaveModelDRalgo::usuage="Saves a model to a file"
LoadModelDRalgo::usuage="Loads a model from file"
DefineDim6::usuage="Defines a dimension 6 operator"
PrintPressureUS::usuage="Calculates the preassure in the ultrasoft theory"
PrintCouplingsEffective::usuage="Prints higher-order couplings"
CounterTerms4D::usuage="Prints 4d CounterTerms"
DefineTensorsUS::usuage="Uses ultrasoft couplings to construct the potential"

(* end of public functions*)


$DRalgoDic=DirectoryName[$InputFileName];


(*
	Functions from groupmath are used to create the model.
*)
If[ Global`$LoadGroupMath,
Get["GroupMath`"];
Print["GroupMath is an independent package, and is not part of DRalgo"];
Print["Please Cite GroupMath: Comput.Phys.Commun. 267 (2021) 108085 \[Bullet] e-Print: 2011.01764 [hep-th]
"];
];


(*
	Verbose=True removes progress messages.
	Mode=2 calculates everything, Mode=1 only calculates LO masses and couplings
	 Mode=0 only calculates LO masses
*)
Options[ImportModelDRalgo] = {Verbose -> False,Mode->2,Dim6->False}


Begin["`Private`"]


(*
	Loads all functions.
*)
Get[FileNameJoin[{$DRalgoDic,"Debye.m"}]]; (*Loads additional functions*)
Get[FileNameJoin[{$DRalgoDic,"HardToSoft.m"}]];(*Loads Hard -> Soft functions*)
Get[FileNameJoin[{$DRalgoDic,"SoftToUS.m"}]];(*Loads Soft -> SS functions*)
Get[FileNameJoin[{$DRalgoDic,"EffPot.m"}]];(*Loads Effective Potential Functions*)
Get[FileNameJoin[{$DRalgoDic,"ModelCreation.m"}]];(*Loads Effective Potential Functions*)


(*
	Defines internal tensors from the loaded model. Also creates help-tensors used for
	intermediate calculations.
*)
ImportModelDRalgo[GroupI_,gvvvI_,gvffI_,gvssI_,\[Lambda]1I_,\[Lambda]3I_,\[Lambda]4I_,\[Mu]ijI_,\[Mu]IJFI_,\[Mu]IJFCI_,YsffI_,YsffCI_, OptionsPattern[]]:=Module[{GroupP=GroupI,gvvvP=gvvvI,gvffP=gvffI,gvssP=gvssI,\[Lambda]1IP=\[Lambda]1I,\[Lambda]3P=\[Lambda]3I,\[Lambda]4P=\[Lambda]4I,\[Mu]ijP=\[Mu]ijI,\[Mu]IJFP=\[Mu]IJFI,\[Mu]IJFCP=\[Mu]IJFCI,YsffP=YsffI,YsffCP=YsffCI},


If[ Global`$LoadGroupMath,
If[!GroupMathCleared && !ValueQ[Global`$GroupMathMultipleModels],
Remove["GroupMath`*"];
GroupMathCleared=True;
];
];

\[Mu]ij=\[Mu]ijP//SparseArray//SimplifySparse;
gvvv=gvvvP//SparseArray//SimplifySparse;
gvss=gvssP//SparseArray//SimplifySparse;
\[Lambda]4=\[Lambda]4P//SparseArray//SimplifySparse;
\[Lambda]3=\[Lambda]3P//SparseArray//SimplifySparse;
Ysff=YsffP//SparseArray//SimplifySparse;
YsffC=YsffCP//SparseArray//SimplifySparse;
gvff=gvffP//SparseArray//SimplifySparse;
\[Mu]IJF=\[Mu]IJFP//SparseArray//SimplifySparse;
\[Mu]IJFC=\[Mu]IJFCP//SparseArray//SimplifySparse;
\[Lambda]1=\[Lambda]1IP//SparseArray//SimplifySparse;
ns=Length[gvss[[1]]];
nv=Length[gvvv];
nf=Length[gvff[[1]]];
\[Lambda]6=EmptyArray[{ns,ns,ns,ns,ns,ns}];

(*Options*)
verbose = OptionValue[Verbose];
mode = OptionValue[Mode]; (*If 2 everthing is calculated. And if 2 only 1-loop contributions are calculated*)
NFMat=IdentityMatrix[nf]//SparseArray; (*This matrix is only relevant if the user wants an arbitrary number of fermion families*)
(*End of Options*)

CT=False; (*Checks if counter-terms have already been calculated*)
DefineGroup[GroupP]; (*Names Debye masses*)
GroupDR=GroupP; (*For saving purposes*)
CreateHelpTensors[] (*Creates recurring tensors*)
];


(*
	Defines a \[Phi]^6 operator
*)
DefineDim6[\[Lambda]6I_]:=Module[{\[Lambda]6P=\[Lambda]6I},
If[mode>=3,
\[Lambda]6=\[Lambda]6P//SparseArray;
,
Print["Please set mode=3 to use this feature"];
];
];


(*
	Takes the user-defined group and names debye masses.
*)
DefineGroup[GroupI_]:=Module[{GroupP=GroupI},
(*The Debye mass matrix is \[Mu]abDef*)
(*At tree-level this matrix is 0. But after DR, thermal masses are named accoording to \[Mu]abDef.*)
\[Mu]abDef=CreateDebyeMasses[GroupP];
VecMassDefined=True;
];


(*
	Enables the user to add an arbitrary number of fermions.
	This works by creating a diagonal matrix. Each element NFMatP[[i,i]] corresponds to how many times the
	fermion labled by ii appears.
*)
DefineNF[NFMatP_]:=Module[{NFMatI=NFMatP},
Do[NFMat[[i[[2]],i[[2]]]]*=i[[1]],{i,NFMatI}];(*Each element is multiplied with nF_i*)
];


(*
	Performs the dimensional reduction for all couplings and masses.
*)

PerformDRhard[]:=Module[{},

If[mode>=0,
ScalarMass[];
VectorMass[];
];

If[mode>=1,
CreateBasisVanDeVis[];

ScalarSelfEnergy[];
VectorSelfEnergy[];
ScalarCubic[];
ScalarQuartic[];
TransverseSSVV[];
LongitudionalSSVV[];
LongitudionalVVVV[];
LongitudionalVVS[];
TadPole[];
];


If[mode>=2,
CounterTerm[];
VectorMass2Loop[];
ScalarMass2Loop[];
TadPole2Loop[];
SymmetricPhaseEnergy[];
];

If[mode>=3,
(*Calculates effective dim 6 operators*)
ScalarSextic[];
];


(*
	This step takes all calculations and removes redundancies. For example, if one element is (g^4 Lb)
	and another 3(g^4 Lb), the function replaces element 2 by three times the first element. This also
	works for linear combinations of elements.
*)

IdentifyTensorsDRalgo[];

];



(*
	Performs the reduction from soft to supersoft.
	ListHardP is a list which tells the code which scalars should be integrated out.
	By default all Debye-vectors (temporal-vectors) are integrated out.
*)

PerformDRsoft[ListHardP_]:=Module[{ListHardI=ListHardP},
(*Options*)

PrepareSoftToSuperSoft[ListHardI];
CreateHelpTensorsSS[];

ScalarSelfEnergySS[];
TadPoleSS[];(*Calculates tadpoles*)
HeavyScalarMassSS[];(*Includes self-energies for heavy-lines with soft momenta*)
VectorSelfEnergySS[];
ScalarVectorCouplingSS[];
ScalarQuarticSS[];
ScalarCubicsSS[];
ScalarMassSS[];
If[mode>=2,
ScalarMass2LoopSS[];
SymmetricPhaseEnergyUS[];
];

(*
	This step takes all calculations and removes redundancies. For example, if one element is (g^4 Lb)
	and another 3(g^4 Lb), the function replaces element 2 by three times the first element. This also
	works for linear combinations of elements.
*)
IdentifyTensorsSSDRalgo[]
];


(*
	Prints constants that appear in the effective theory.
*)
PrintConstants[]:=Module[{},
ToExpression[StringReplace[ToString[StandardForm[{Lb->(Log[\[Mu]^2/T^2]-2 Log[4 \[Pi]]+2EulerGamma),Lf->(Log[\[Mu]^2/T^2]-2 Log[4 \[Pi]]+2EulerGamma+4 Log[2])}]],"DRalgo`Private`"->""]]
];


(*
	Replaces the Glaisher constant by c, which is sometimes used in the litterature. See for example hep-ph: 9508379
*)
PrintGenericBasis[]:=Module[{},
ToExpression[StringReplace[ToString[StandardForm[{Log[Glaisher]->-1/12 (Lb+2cplus-EulerGamma)}]],"DRalgo`Private`"->""]]
];


(*
	Checks if two variables are identical up to a numerical factor
	For example if a=x^2 y and b=y, then there's no relation, and the function returns nothing.
	While if a=x^2 and b=5 x^2 the function returns 5
*)
CompExp2[a0_,b0_]:=Module[{a=a0,b=b0},
Comps=Solve[a v[1]- b ==0,v[1]]//Simplify//DeleteDuplicates//Select[#,UnsameQ[#,{}]&]&;
Comps/. {v[x_]->a_}:>a
]


(*
	Old function: Should maybe be removed.
	This function finds linear dependencies between the variables in list.
	Only numerical factors are included.
	The result is given in Mat. Where Mat[[i,j]]=0 if list[[i]]!=N list[[j]] where N is a number.
	Otherwise Mat[[i,j]]=N if list[[i]]=Nlist[[j]]
*)
OverallFac[list_,varMat_]:=Module[{L=list,Mat=varMat},

For[i=1,i<Length[list],i++,
For[j=i+1,j<Length[list]+1,j++,
help=CompExp2[list[[i]],list[[j]]];
TempIf=If[Length[help]>0,{NumericQ/@help}[[1,1]]&&Equal@@NumericQ/@help,False];
If[TempIf,a=help[[1]],a=0];
Mat[[i,j]]=a;
];
];
Mat
]


(*
	Old function:Should maybe be removed.
	This function finds linear dependences between the variables in listvar
	Only relations between two variables are considered.
	This function is only used for the soft->supersoft step where the faster, and more general,
	RelationsBVariables3 has problems due to inverse powers of masses.
*)
RelationsBVariables[list_,listVar_]:=Module[{L=list,LV=listVar},
LVTemp=LV;
If[L[[1]]==0,
	L=Delete[L,1];
];

Mat=ConstantArray[0,{Length[L],Length[L]}];(*Linear-dependency matrix*)
Mat1=OverallFac[L,Mat];
TempVar=ConstantArray[0,{Length[LV]}];
For[i=1,i<Length[LV],i++,
For[j=i+1,j<Length[LV]+1,j++,
IfTemp=Mat1[[1;;i,j]]//DeleteDuplicates//DeleteCases[#,0,Infinity]&; (*Removes cases when the numerical factor is 0 or infinity*)
If[Length[IfTemp]==1&&Mat1[[i,j]]!=0,
LVTemp[[j]]=Mat1[[i,j]]LV[[i]];(*Creates a list with all linear relations*)
];
];
];
LVTemp
]


(*
	Creates a list of all possible couplings that can appear in
	1-loop matching relations
	Note that this does not include couplings in debye/scalar masses
	at 1-loop level
*)
CreateBasisVanDeVis[]:=Module[{},
(* This module received founding from the fish *)

FermionFamVar=Normal[NFMat]//Variables;
ScalVar=\[Lambda]4//Normal//Variables;


GaugeVarPre=Join[Normal[gvvv]//Variables,Normal[gvff]//Variables,Normal[gvss]//Variables]//DeleteDuplicates; (*Includes possible non-numeric charges*)
GaugeCharge=Complement[GaugeVarPre,GaugeCouplingNames];(*Possible non-numeric gauge charges*)

GaugeVarHelp=Table[i*j,{i,GaugeCouplingNames},{j,GaugeCharge}];
GaugeVar=Join[GaugeCouplingNames,GaugeVarHelp];

YukVar=Normal[Ysff]//Variables;
AuxVar={1,Lb,Lf};
AuxVar=Join[AuxVar,FermionFamVar]//DeleteDuplicates;
varHelp=Join[ScalVar,GaugeVar,YukVar];
t2=Table[i*j,{i,varHelp},{j,varHelp}]//Flatten[#]&//DeleteDuplicates;
t3G=Table[i*j*k,{i,ScalVar},{j,GaugeVar},{k,GaugeVar}]//Flatten[#]&//DeleteDuplicates;
t3F=Table[i*j*k,{i,ScalVar},{j,YukVar},{k,YukVar}]//Flatten[#]&//DeleteDuplicates;
t4=Table[i*j*k*l,{i,GaugeVar},{j,GaugeVar},{k,GaugeVar},{l,GaugeVar}]//Flatten[#]&//DeleteDuplicates;
t4F=Table[i*j*k*l,{i,YukVar},{j,YukVar},{k,YukVar},{l,YukVar}]//Flatten[#]&//DeleteDuplicates;
t4FG=Table[i*j*k*l,{i,YukVar},{j,YukVar},{k,GaugeVar},{l,GaugeVar}]//Flatten[#]&//DeleteDuplicates;
basPre=Join[varHelp,t2,t3G,t3F,t4,t4F,t4FG];

basDR=Table[i*j*k,{i,basPre},{j,AuxVar},{k,{1,T,T^2}}]//Flatten[#]&//DeleteDuplicates;
]


{basDR};


(*
	This function finds linear dependences between the variables in list.
	This works by treating each element in list as a vector in the space spanned by basDR.
	All vectors are then rowreduced, and a minimal set of basis vectors (in list) are found.
*)
RelationsBVariables3[list_]:=Module[{L=list},
(*Creates a vector-basis*)

(*One could say that v3 and v2 are identical. With v3 being almost twice as identical as v2*)
If[L[[1]]==0&&Length[L]>1,
	Lp=Delete[L,1];
,
	Lp=L;
];


varHelp=Lp//Variables;
varFix=#->0&/@varHelp; (*Trick to ensure that vectors are expanded properly*)

(*
Expands all elements in Lp in terms of basDR.
*)
setVecs=Table[Coefficient[i,basDR],{i,Lp}]/.varFix;  
(*Delete columns with only 0s*)
setVecs=Transpose[DeleteCases[Transpose[setVecs], {0 ..}, Infinity]];

(*Finds independent basis*)
rr = setVecs // Transpose // RowReduce;
rr=DeleteCases[rr, {0 ..}, Infinity];
(*Basis elements*)
basisElements = Flatten[FirstPosition[#, 1, Nothing] & /@ rr];
varBasis=Table[ \[Lambda]VL[a],{a,1,Length[basisElements]}];

(*Puts everything together*)
LVTemp=LVTemp=ConstantArray[0,Length[Lp]];
Do[LVTemp[[i]]=rr[[;;,i]] . varBasis,{i,1,Length[LVTemp]}];

Return[LVTemp];
]


myPrint[args__,{style__}]:=Print[Row[{args},BaseStyle->{style}]]


(*
	Old function: Should maybe be removed.
	This function finds linear dependencies between the variables in list.
	Only numerical factors are included.
	The result is given in Mat. Where Mat[[i,j]]=0 if list[[i]]!=N list[[j]] where N is a number.
	Otherwise Mat[[i,j]]=N if list[[i]]=Nlist[[j]]
*)
OverallFac2[list_,varMat_]:=Module[{L=list,Mat=varMat},

TotVar=L//Variables;
DoneList=ConstantArray[0,1];

For[i=1,i<Length[list],i++,
Var=list[[i]]//Variables;
varCompliment=Complement[TotVar,Var];
For[j=i+1,j<Length[list]+1,j++,
If[!MemberQ[DoneList, j],
If[!CheckVariables[list[[j]],varCompliment],
h=CompExp3[list[[i]],list[[j]]];
Mat[[i,j]]=h;
If[h!=0,AppendTo[DoneList,j]];
];

];
];
];
Mat
]


(*
	This function finds linear dependencies between the variables in list.
	Only numerical factors are included.
	The result is given in Mat. Where Mat[[i,j]]=0 if list[[i]]!=N list[[j]] where N is a number.
	Otherwise Mat[[i,j]]=N if list[[i]]=Nlist[[j]]
*)
RelationsBVariables2[list_,listVar_]:=Module[{L=list,LV=listVar},
LVTemp=LV//FullSimplify;
Mat=ConstantArray[0,{Length[Delete[L,1]],Length[Delete[L,1]]}];(*Linear-dependency matrix*)
Mat1=OverallFac2[Delete[L,1],Mat];
TempVar=ConstantArray[0,{Length[LV]}];
For[i=1,i<Length[LV],i++,
For[j=i+1,j<Length[LV]+1,j++,
IfTemp=Mat1[[1;;i,j]]//DeleteDuplicates//DeleteCases[#,0,Infinity]&; (*Removes cases when the numerical factor is 0 or infinity*)
If[Length[IfTemp]==1&&Mat1[[i,j]]!=0,
LVTemp[[j]]=Mat1[[i,j]]LV[[i]](*Creates a list with all linear relations*)
];
];
];
LVTemp
]


CheckVariables[a_,Vars_]:=MemberQ[Boole@(MemberQ[a, #, {0, -1}, Heads -> True]&/@Vars),1, {0, -1}, Heads -> True];


(*
	Checks if two variables are identical up to a numerical factor
	For example if a=x^2 y and b=y, then there's no relation, and the function returns nothing.
	 While if a=x^2 and b=5 x^2 the function returns 5.
*)
CompExp3[a0_,b0_]:=Module[{a=a0,b=b0},
Temp=Simplify[b/a];
If[NumericQ[Temp],Return[Temp],Return[0]];
]


(*
	Converts an array to a saveable form
*)
ConvertToSaveFormat[tens_SparseArray]:=Module[{},
Return[{tens["NonzeroPositions"],tens["NonzeroValues"],Dimensions[tens]}]
];


(*
	Converts imported data, as defined by ConvertToSaveFormat, to a sparse array
*)
ConvertToSparse[arr_]:=Module[{},
Return[SparseArray[arr[[1]]->arr[[2]],arr[[3]]]]
];


(*
	Loads the position of all scalar,gauge, and fermion representations.
*)
LoadRepPositions[repPos_]:=Module[{repPosP=repPos},
ScalarVariablesIndices=repPosP[[1]];
GaugeIndices=repPosP[[2]];
FermionVariablesIndices=repPosP[[3]];
];


(*
	Loads the names of gauge couplings.
*)
LoadCouplingNames[couplingNamesI_]:=Module[{couplingNamesP=couplingNamesI},
GaugeCouplingNames=couplingNamesP;
];


(*
	Saves a model by converting all coupling-tensors to a list.
*)
SaveModelDRalgo[modelInfo_,fileName_]:=Module[{modelInfoP=modelInfo},

PosScalar=PrintScalarRepPositions[];
PosVector=PrintGaugeRepPositions[];
PosFermion=PrintFermionRepPositions[];
PosReps={PosScalar,PosVector,PosFermion};
tensP={modelInfoP,PosReps,GaugeCouplingNames,GroupDR,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJFC,\[Mu]IJF,Ysff,YsffC};
SaveFile={tensP[[1]],tensP[[2]],tensP[[3]],tensP[[4]]}; (*The fourth element is the group*)
tensP=Delete[Delete[Delete[Delete[tensP,1],1],1],1];
Do[
AppendTo[SaveFile,ConvertToSaveFormat[i]];
,{i,tensP}];
Export[fileName,SaveFile];
];


(*
	Loads tensors that are saved by SaveModelDRalgo
*)
LoadModelDRalgo[fileName_]:=Module[{},
arrImp=ReadList[fileName];
InfoText=arrImp[[1]];(*The first element is the info*)
arrImp=Delete[arrImp,1];

LoadRepPositions[arrImp[[1]]];(*The Second element is the repPositions*)
arrImp=Delete[arrImp,1];

LoadCouplingNames[arrImp[[1]]];(*The Third element is the gauge-coupling names*)
arrImp=Delete[arrImp,1];

ImportFile={arrImp[[1]]};(*The fourth element is the group*)
arrImp=Delete[arrImp,1];

Do[
AppendTo[ImportFile,ConvertToSparse[i]];
,{i,arrImp}];
(*Prints the info text*)
Print[Grid[{{Row[InfoText,"\n",BaseStyle->(FontFamily->"Consolas")]}},Alignment->{Left,Center}]];
(**)
Return[ImportFile]
];


(*
	All private constants.
*)


{ZLij,GvvssTSS,\[Lambda]3DSS,\[Mu]ijSSLO,\[Mu]ijSSNLO,\[Mu]ijSNLOSS,\[Lambda]3DSS,\[Lambda]KVecTSS,\[Lambda]3CTot,\[Lambda]3CSSS,ZSij,\[Mu]ijSSNLO2,Ggvvv};


{TadPoleS,ContriTadPoleSoftToHard,GgvvvSS,\[Lambda]4SMod,\[Lambda]4Tot,IdentMatPre};


{\[Beta]gvff,Zgvff,\[Mu]ijEP,gvvvEP,gvssEP,\[Lambda]4EP,\[Lambda]3EP,nsEP,nvEP,CT,HelpSolveEffectiveHardM};


{aS3D,ZijS,aV3D,ZabT,ZabL,\[Lambda]3D,GvvssL,GvvssT,\[Lambda]AA,\[Lambda]3CS,\[Mu]SijNLO,GvvsL};(*DimRed Results*)


{\[CapitalLambda]\[Lambda],\[CapitalLambda]g,Hg,Habij,HabIJF,HabIJFC,Ysij,YsijC,YTemp,YTempC,Yhelp,YhelpC};(*Private Variables*)


{\[Gamma]ij,\[Beta]mij,\[Beta]\[Lambda]ijkl,Z\[Lambda]ijkl,\[Gamma]ab,\[Beta]vvss,Zgvvss,\[Beta]gvvv,Zgvvv,\[Gamma]IJF,\[Beta]Ysij,ZYsij,\[Beta]YsijC,ZYsijC};(*CounterTerms*)


{\[Lambda]3DS,\[Lambda]KVecT,\[Lambda]KVec,\[Lambda]AAS,IdentMat,\[Mu]ijVNLO,\[Mu]ijSNLO,\[Lambda]3CSRed,\[Lambda]1IP,NFMat,NSMat};(*Private Variables*)


{nsH,nSl,\[Lambda]K,\[Lambda]4S,\[Lambda]4K,\[Lambda]x,\[Lambda]y,gAvss,gvssL,\[Mu]ijL,\[Mu]ijLight};


{HabijL,HabijVL,HabijA,HabijVA,\[Lambda]3Cx,\[Lambda]3Cy,\[Lambda]3CLight,\[Lambda]3CHeavy,\[Mu]IJF,\[Mu]IJFC,\[Mu]VabNLO,\[Mu]abDef,GroupMathCleared}


End[]
EndPackage[]
