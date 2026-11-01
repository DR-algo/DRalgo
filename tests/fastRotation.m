(* ::Package:: *)

If[$InputFileName=="",
	SetDirectory[NotebookDirectory[]],
	SetDirectory[DirectoryName[$InputFileName]]
];
DRalgo`DRalgo`$LoadGroupMath=True;
(*DRalgo`DRalgo`$InstallGroupMath=True;*)

Check[
    Get["../Kernel/DRalgo.m"],
    Message[Get::noopen, "DRalgo` at "<>ToString[$UserBaseDirectory]<>"/Applications"];
    Abort[];
]


(* ::Chapter:: *)
(*RotateTensorsCustomMass with FastRotation*)


(*
	Issue #52: with FastRotation->True, RotateTensorsCustomMass rotates only the components lam4_iikk
	of the quartic couplings, which the two-loop potentials need. The rotated quartic is a rank-4 tensor
	whose only non-zero components are lam4_iikk, which CalculatePotentialUS[] and CalculatePotentialHET[]
	take as they take the fully rotated one.

	Checked here, for the model of examples/25_tripletHET.m in the ultrasoft theory and an exact orthogonal
	rotation that mixes the doublet, the triplet, and the neutral vectors:
	- the rotated quartic is a rank-4 tensor, and its components lam4_iikk agree with the full rotation;
	- the two-loop potential of CalculatePotentialUS[] agrees with the full rotation;
	- CalculatePotentialHET[] gives no messages, and its potential and PrintScalarKineticHET[] agree with
	  the full rotation, with the triplet and the SU(2) vectors integrated out.
*)


testList={};


(* ::Section:: *)
(*Model*)


Group={"SU3","SU2","U1"};
RepAdjoint={{1,1},{2},0};
HiggsDoublet1={{{0,0},{1},1/2},"C"};
HiggsTriplet1={{{0,0},{2},0},"R"};
RepScalar={HiggsDoublet1,HiggsTriplet1};
CouplingName={g3,gw,g1};

Rep1={{{1,0},{1},1/6},"L"};
Rep2={{{1,0},{0},2/3},"R"};
Rep3={{{1,0},{0},-1/3},"R"};
Rep4={{{0,0},{1},-1/2},"L"};
Rep5={{{0,0},{0},-1},"R"};
RepFermion1Gen={Rep1,Rep2,Rep3,Rep4,Rep5};
RepFermion3Gen={RepFermion1Gen,RepFermion1Gen,RepFermion1Gen}//Flatten[#,1]&;

{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion3Gen,RepScalar];

MassTerm1=CreateInvariant[Group,RepScalar,{{1,1},{True,False}}][[1]]//Simplify//FullSimplify;
MassTerm2=CreateInvariant[Group,RepScalar,{{2,2},{True,True}}][[1]]//Simplify//FullSimplify;
\[Mu]ij=GradMass[mH*MassTerm1+m\[CapitalSigma]/2*MassTerm2]//Simplify//SparseArray;
\[Lambda]4=GradQuartic[\[Lambda]*MassTerm1^2+b4/4*MassTerm2^2+a2/2*MassTerm1*MassTerm2];
YukawaDoublet1=CreateInvariantYukawa[Group,RepScalar,RepFermion3Gen,{{1,1,2},{False,False,True}}]//Simplify;
Ysff=-GradYukawa[yt*YukawaDoublet1[[1]]];
YsffC=SparseArray[Simplify[Conjugate[Ysff]//Normal,Assumptions->{yt>0}]];

ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PerformDRhard[];
PerformDRsoft[{}];


(* ::Section:: *)
(*Rotation*)


UseUltraSoftTheory[];
DefineVEVS[{0,0,0,\[CurlyPhi],0,0,0}//SparseArray];
scalarMass=PrintTensorsVEV[1]//Normal;
vectorMass=PrintTensorsVEV[2]//Normal;

(*Exact rotations of the pairs {i,j} of fields by the angles with cosine c and sine s*)
rotation[n_,pairs_]:=Module[{d=IdentityMatrix[n]},
	Do[d[[p[[1]],p[[1]]]]={{p[[2]],-p[[3]]},{p[[3]],p[[2]]}},{p,pairs}];
	SparseArray[d]
];
(*The Goldstone mode 1 with 2, the Higgs 4 with the triplet 5, and the neutral vectors 11 and 12*)
DS=rotation[7,{{{1,2},3/5,4/5},{{4,5},5/13,12/13}}];
DV=rotation[12,{{{11,12},8/17,15/17}}];

(*The diagonal of the rotated mass matrices serves as the custom masses*)
scalarMassRotated=DiagonalMatrix[Diagonal[Transpose[DS] . scalarMass . DS]]//SparseArray;
vectorMassRotated=DiagonalMatrix[Diagonal[Transpose[DV] . vectorMass . DV]]//SparseArray;

(*Random positive values of the parameters with 40 digits, for numerical comparisons*)
SeedRandom[52];
numerical[expr_]:=Module[{vars=Variables[Level[expr,{-1}]~Join~Cases[expr,_lamVL|_lamVLL,Infinity]]},
	expr/.Thread[vars->N[RandomInteger[{1,100},Length[vars]]/100,40]]
];
results[fast_]:=Module[{lam4,potential,het,kinetic},
	DefineVEVS[{0,0,0,\[CurlyPhi],0,0,0}//SparseArray];
	RotateTensorsCustomMass[DS,DV,scalarMassRotated,vectorMassRotated,FastRotation->fast];
	lam4=PrintTensorsVEV[3];
	CalculatePotentialUS[];
	potential=PrintEffectivePotential["NNLO"];
	PrepareHET[{5,6,7},{9,10,11}];
	het=Check[CalculatePotentialHET[];PrintActionHET[],$Failed];
	kinetic=Normal[PrintScalarKineticHET[]];
	{lam4,potential,het,kinetic}
];


(* ::Section:: *)
(*Comparison*)


{lam4Full,potentialFull,hetFull,kineticFull}=results[False];
{lam4Fast,potentialFast,hetFast,kineticFast}=results[True];
values=numerical[{potentialFull-potentialFast,hetFull-hetFast,kineticFull-kineticFast}/.\[CurlyPhi]->N[7/10,40]];

AppendTo[testList,
VerificationTest[
	{Dimensions[lam4Fast],Simplify[Table[lam4Full[[i,i,k,k]]-lam4Fast[[i,i,k,k]],{i,7},{k,7}]]},
	{{7,7,7,7},ConstantArray[0,{7,7}]},
	TestID->"rank-4-quartic-with-lam4-iikk"
]];

AppendTo[testList,
VerificationTest[
	Chop[values[[1]],10^-20],
	0,
	TestID->"two-loop-potential-ultrasoft"
]];

AppendTo[testList,
VerificationTest[
	{hetFast=!=$Failed,Chop[values[[2]],10^-20],Chop[values[[3]],10^-20]},
	{True,0,ConstantArray[0,{4,4}]},
	TestID->"het-potential-and-kinetic-terms"
]];


(* ::Section:: *)
(*Report*)


report=TestReport[testList];
Print[report];
Print[Column[{#["TestID"],#["Outcome"]}&/@Values[report["TestResults"]]]];
Exit[If[report["AllTestsSucceeded"],0,1]]
