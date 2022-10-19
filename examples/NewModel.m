(* ::Package:: *)

(*Quit[];*)


SetDirectory[NotebookDirectory[]];
$LoadGroupMath=True;
<<../DRalgo.m


(* ::Chapter:: *)
(*Fun with DRalgo*)


(* ::Section::Closed:: *)
(*Model for Standard-Model like theory*)


Group={"SU3","SU2","U1"}; (*This specifies the gauge group. So for the Standard Model SU3xSU2xU1 *)


(*To add particles we need to specifiy the representation*)
(*This is done by giving Dynking indices. For a particular rep: Rep={DynkinSU3,DynkinSU2,DynkinU1};*)


(*These Dynkin indices can be found in various group theory books, and also the paper by Slansky*)
(*Yet for practical purposes they are directly given by GroupMath. See the accompanying "RepNames.m" file*)


RepVector={{1,1},{2},0}; (*For the vector bosons we have a colour octed {1,1}, a Weak Triplet {2}, and no charge under U1*)


HiggsDoublet={{{0,0},{1},1/2},"C"}; (*The Higgs doublet is uncharged under colour {0,0}, transforms as a doublet under SU2 {2}*)
									(* and has Hypercharge 1/2*)
									(*The Higgs is complex "C". For a real representation a "R" should be given*)


Rdoublet={{{1,0},{1},1/6},"C"}; 


SSinglet={{{0,1},{0},1/3},"C"}; 


RepScalar={HiggsDoublet,Rdoublet,SSinglet};
CouplingName={gs,gw,gY};


(*We now need to specify the fermions*)


Rep1={{{1,0},{1},1/6},"L"};  (*QL- Left-handed quark doublet with hypercharge 1/6*)
Rep2={{{1,0},{0},2/3},"R"};  (*tR- Right-handed top quark with hypercharge 2/3*)
Rep3={{{1,0},{0},-1/3},"R"}; (*bR- Right-handed bottom quark with hypercharge -1/3*)
Rep4={{{0,0},{1},-1/2},"L"}; (*L- Left-handed lepton doublet with hypercharge -1/2*)
Rep5={{{0,0},{0},-1},"R"};   (*L- Right-handed lepton with hypercharge -1*)
RepFermion1Gen={Rep1,Rep2,Rep3,Rep4,Rep5}; (*This specifies one generation of fermions*)


(*Rep1={{{1,0},{1},Yu/2},"L"};
Rep2={{{1,0},{0},Yd/2},"R"};
Rep3={{{1,0},{0},Yd/2},"R"};
Rep4={{{0,0},{1},Yl/2},"L"};
Rep5={{{0,0},{0},Ye/2},"R"};
RepFermion1Gen={Rep1,Rep2,Rep3,Rep4,Rep5};*)


(*W can then stack the generations together if we want multiple generations. In this case 3*)
(*See the manual for how to create an arbitrary number of NF generations*)


RepFermion3Gen={RepFermion1Gen,RepFermion1Gen,RepFermion1Gen}//Flatten[#,1]&; 


(* ::Text:: *)
(*The input for the gauge interactions toDRalgo are then given by*)


(*The tensors on the left-hand side are used in many applications. For example general 2-loop beta functions*)
(*and two/three loop effective potentials for any models*)


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepVector,CouplingName,RepFermion3Gen,RepScalar];


(* ::Text:: *)
(*Specifying the scalar potential :*)


(*You don't have to use DRalgo to define your model. You can also write down everything yourself, or import*)
(*the result from another code.*)


(*Let us start with the mass. To get a mass term with two Higgs fields we specify that we want two RepScalar[[1]] fields*)


RepScalar={HiggsDoublet,Rdoublet,SSinglet};


InputInv={{1,1},{True,False}};  (*\[CapitalPhi]\[CapitalPhi]^+. False-> Hermitian conjugate*)
MassTerm1=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;


InputInv={{2,2},{True,False}};
MassTerm2=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;


InputInv={{3,3},{True,False}};
MassTerm3=CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify;


(*To create the mass matrix we then add all our mass terms together*)


VMass=m2*MassTerm1+\[Mu]R*MassTerm2+\[Mu]S*MassTerm3;


\[Mu]ij=GradMass[VMass]//Simplify//SparseArray; (*GradMass[] takes the Lagrangian terms and extracts the mass matrix*)


(*Let us now move on to quartic couplings*)
(*We can either define these by writing InputInv={{1,1,1,1},{True,True,False,False}}; as above*)
(*Or, since we already have the mass terms we can just use those.*)


QuarticTerm1=\[Lambda]H*MassTerm1^2;
QuarticTerm2=gHR MassTerm1*MassTerm2;
QuarticTerm3=gHS MassTerm1*MassTerm3;


InputInv={{1,1,2,2},{True,False,True,False}}; 
QuarticTerm4=gHp CreateInvariant[Group,RepScalar,InputInv][[2]]//Simplify;


InputInv={{2,2,2,2},{True,True,False,False}}; 
QuarticTerm4=\[Lambda]S1 CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;


InputInv={{2,2,2,2},{True,True,False,False}}; 
QuarticTerm5=\[Lambda]S2 CreateInvariant[Group,RepScalar,InputInv][[2]]//Simplify//FullSimplify;


InputInv={{3,3,3,3},{True,True,False,False}}; 
QuarticTerm6=\[Lambda]R1 CreateInvariant[Group,RepScalar,InputInv][[1]]//Simplify//FullSimplify;


InputInv={{2,2,3,3},{True,False,True,False}}; 
QuarticTerm7=gRS CreateInvariant[Group,RepScalar,InputInv][[2]]//Simplify;


(*If we have a a Higgs doublet: \[CapitalPsi]=1/(Sqrt[2])(\[Phi]1+\[ImaginaryI]\[Psi]1,\[Phi]2+\[ImaginaryI] \[Psi]2) the code stores this as (\[Phi]1,\[Phi]2,\[Psi]1,\[Psi]2)*)


VQuartic=QuarticTerm1+QuarticTerm2+QuarticTerm3+0*(QuarticTerm4+QuarticTerm5+QuarticTerm6+QuarticTerm7); (*We then add all quartic terms together*)


\[Lambda]4=GradQuartic[VQuartic]; (*GradQuartic[] then extracts the quartic couplings*)


(*To define a Yukawa-coupling we need to specify the Scalar+First fermion+Second fermion*)
(*So for a \[CapitalPsi]^+ QL^+ tR term we would write:*)


(*We now need to specify the fermions*)


InputInv={{1,1,2},{False,False,True}};  
YukawaDoublet=yt*CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv][[1]]//Simplify;


InputInv={{3,1,4},{True,True,True}};  
YukawaTerm2=y\[CapitalTheta]*CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv][[1]]//Simplify;


InputInv={{2,4,3},{False,False,True}};  
YukawaTerm3=y\[CapitalOmega]*CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv][[1]]//Simplify;


InputInv={{3,2,5},{True,True,True}};  
YukawaTerm4=yY*CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,InputInv][[1]]//Simplify;


Ysff=-GradYukawa[YukawaDoublet+YukawaTerm2+YukawaTerm3+YukawaTerm4]; (*GradYukawa[] then extracts the Yukawa coupling*)


YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt>0,y\[CapitalTheta]>0,y\[CapitalOmega]>0,yY>0}]]; (*YsffC is always the complex conjugate of Ysff*)


(* ::Section:: *)
(*Dimensional Reduction*)


(*We can now find the effective couplings. Verbose->False turns of messages*)
(*1-loop thermal masses: Mode->0; 1-loop masses and couplings: Mode->1; 2-loop masses and couplings: Mode->2*)
(*Dimension 6 operators: Mode->3*)


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->True,Mode->2];


(*PerformDRhard[]  performs the matching from 4d to 3d*)


PerformDRhard[] 


(*PrintCouplings[]shows the effective couplings*)


PrintCouplings[] 


(*Note that Lb and Lf contain RG-dependent terms. They can be shown by writing PrintConstants[]*)


PrintConstants[]


(*PrintScalarMass["LO"] shows the tree-level mass+ the 1-loop thermal mass*)


PrintScalarMass["LO"]


PrintTemporalScalarCouplings[]


PrintPressure["LO"]


(*PrintScalarMass["NLO"] shows the two-loop thermal mass*)


PrintScalarMass["NLO"]


(*So for the complete NLO mass you would need PrintScalarMass["LO"]+PrintScalarMass["NLO"]*)


(*Next debye masses*)


(*PrintDebyeMass["LO"] shows the 1-loop thermal mass*)


PrintDebyeMass["LO"]


(*PrintScalarMass["NLO"] shows the two-loop thermal mass*)


PrintDebyeMass["NLO"]


(*So for the complete NLO mass you would need PrintDebyeMass["LO"]+PrintDebyeMass["NLO"]*)


(*We can also find various quartic couplings for the temporal gauge field*)


PrintTemporalScalarCouplings[]


(*DRalgo also calculates normal beta functions *)


BetaFunctions4D[]


(*Anomalous dimensions can also be found for all fields, fermions included*)
(*The notation is d /dlog\[Mu]Subscript[R, a]=Subscript[\[Gamma], ab]Subscript[R, b]*)


(*First we find where all fields live.*)


PosScalar=PrintScalarRepPositions[];
PosVector=PrintGaugeRepPositions[];
PosFermion=PrintFermionRepPositions[];


(*We can then find the anomalous dimensions by writing*)


(*Here AnomDim4D["S",{a,b}] corresponds to Subscript[(\[Gamma]^s), ab]. So Subscript[(\[Gamma]^s), 11] is the anomalous dimension for the Higgs field*)


Table[AnomDim4D["S",{a,b}],{a,PosScalar},{b,PosScalar}]
Table[AnomDim4D["V",{a,b}],{a,PosVector},{b,PosVector}]
Table[AnomDim4D["F",{a,b}],{a,PosFermion},{b,PosFermion}]


(* ::Section:: *)
(*Integrating out the temporal component*)


(*We can also integrate out the temporal vector component by writing PerformDRsoft[{}]*)


PerformDRsoft[{},IncludeCubics->"False"]


(*The resulting theory is known as the "ultrasoft" theory and all the commands are as above*)


PrintCouplingsUS[]


PrintScalarMassUS["LO"]


PrintScalarMassUS["NLO"]


BetaFunctions3DUS[]


(*In addition, if there are heavy scalars, these can also be integrated out simultaneously with the temporal vector bosons*)


CommunicatingRowAndColumns[matrix_]:=Module[{allIndices,idxs,continue,newIdxs,newGroup,result},
allIndices=Range[Length[matrix]];

result=Reap[While[Length[allIndices]>0,
idxs={allIndices[[1]]};
continue=True;
While[continue,
newIdxs=Sort[DeleteDuplicates[Join[idxs,ArrayRules[matrix[[idxs]]][[1;;-2,1,2]]]]];
If[Sort[idxs]==Sort[newIdxs],
continue=False;,
continue=True;
];
idxs=newIdxs;
];
newGroup=idxs;
Sow[newGroup];
allIndices=Complement[allIndices,newGroup];
]][[2,1]];

Return[result];
]


(* ::Section:: *)
(*2 - Loop Effective potential*)


vev=Table[0,{i,1,22}];


vev[[1]]=\[Phi];


vev[[5]]=0;


(*There are two options to define your model*)


(*First, if you have used for example "PerformDRsoft[{}]" the couplings can be defined with the command:*)


DefineTensorsUS[]


(*Second, you can define your custom model by writing "DefineNewTensorsUS[\[Mu]ij,\[Lambda]4,\[Lambda]3,gvss,gvvv];"(here just taking the original model as an example)*)


DefineNewTensorsUS[\[Mu]ij,\[Lambda]4,\[Lambda]3,gvss,gvvv];
 \[CurlyPhi]vev=vev//SparseArray;
  DefineVEVS[\[CurlyPhi]vev];
PrintTensorsVEV[];


(*The vector-mass matrix is not diagonal. This can be seen by looking at*)


PrintTensorsVEV[][[2]]//Normal


(*So the problem is the A^3-B mixing as usual.*)


(*We can diagonalize the mass-matrix via*)


MassMatrix=PrintTensorsVEV[];
VectorMass=MassMatrix[[2]]//Normal;
VectorEigenvectors=FullSimplify[
    Transpose[Normalize/@Eigenvectors[VectorMass[[11;;12,11;;12]]]],
Assumptions->{gw>0,gY>0,\[CurlyPhi]>0}]; DVRot={{IdentityMatrix[10],0},{0,VectorEigenvectors}}//ArrayFlatten; 
DSRot=IdentityMatrix[22];
RotateTensorsUSPostVEV[DSRot,DVRot];


PrintTensorsVEV[][[2]]//Normal


(* ::Subsection:: *)
(*Calculating the potential*)


(*After the massmatrices are diagonal, the potential is given by CalculatePotentialUS[];*)


CalculatePotentialUS[]


(*We can first give the tree-level potential*)


PrintEffectivePotential["LO"]


(*Next the 1-loop potential*)


PrintEffectivePotential["NLO"]


(*We can set the hypercharge coupling to zero to make the result neater*)


PrintEffectivePotential["NLO"]/.gY3dUS->0


(*And finally, the two-loop effective potential is given by*)


PrintEffectivePotential["NNLO"]/.gY3dUS->0
