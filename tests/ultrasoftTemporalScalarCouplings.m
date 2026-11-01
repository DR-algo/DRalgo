(* ::Package:: *)

If[$InputFileName=="",
	SetDirectory[NotebookDirectory[]],
	SetDirectory[DirectoryName[$InputFileName]]
];
(*Several models are created in the same kernel*)
DRalgo`DRalgo`$GroupMathMultipleModels=True;

DRalgo`DRalgo`$LoadGroupMath=True;

Check[
    Get["../Kernel/DRalgo.m"],
    Message[Get::noopen, "DRalgo` at "<>ToString[$UserBaseDirectory]<>"/Applications"];
    Abort[];
]


(* ::Chapter:: *)
(*Couplings of the temporal scalar to light and heavy scalars in the ultrasoft theory*)


(*
	The couplings of the soft theory that arise at one loop, of two temporal scalars to a scalar, kap3/2 A0^2 R_i, and to two
	scalars of different flavors, kap4/4 A0^2 R_i R_j, enter the matching to the ultrasoft theory like the other couplings, also
	when heavy scalars are integrated out. With the tadpole integral A[x]=-Sqrt[x]/(4Pi) and the bubble B[x,x]=1/(8Pi Sqrt[x]):
	- scalar QED with a light complex scalar phi1 and a heavy one phi2 of the same charge, with l6 (phi1^+phi1)(phi1^+phi2+h.c.):
	  the A0 loop mixes the real components of phi1 and phi2 by Pi=kap4 A[mD^2]/2, so that phi2=-(Pi+l6 |phi1|^2)phi1/M and
	  V ⊃ -2 Pi l6 |phi1|^4/M, and the quartic l1 (phi1^+phi1)^2 gets -2 Pi l6/M=kap4 l6 mD/(4Pi M);
	- scalar QED with a light real singlet f and a heavy complex scalar chi with cf f chi^+chi: the A0 loop of the coupling
	  to f gives the mass of f the bubble -kap3^2 B[mD^2,mD^2]/2=-kap3^2/(16Pi mD), as without heavy scalars.
*)


testList={};
private[s_]:=ToExpression["DRalgo`DRalgo`Private`"<>s];
toGlobal[x_]:=x/.{private["lamVL"]->lamVL,private["lamVVSL"]->lamVVSL};

Group={"U1"};RepAdjoint={0};CouplingName={g1};RepFermion={};


(* ::Section:: *)
(*Temporal scalars coupled to a light and a heavy scalar*)


RepScalar={{{1},"C"},{{1},"C"}};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
p11=CreateInvariant[Group,RepScalar,{{1,1},{True,False}}][[1]];
p22=CreateInvariant[Group,RepScalar,{{2,2},{True,False}}][[1]];
p12=CreateInvariant[Group,RepScalar,{{1,2},{True,False}}][[1]];
p21=CreateInvariant[Group,RepScalar,{{2,1},{True,False}}][[1]];
\[Mu]ij=GradMass[m1 p11+m2 p22]//Simplify//SparseArray;
\[Lambda]4=GradQuartic[l1 p11^2+l2 p22^2+l6 p11(p12+p21)]//Simplify;
ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PerformDRhard[];
PerformDRsoft[{3,4}];

(*The coupling of the temporal scalar to the real parts of phi1 and phi2*)
kap4=toGlobal[Normal[private["lamKVec"]][[1,1,1,3]]];
quartic=Expand[l13dUS/.PrintCouplingsUS[]];

AppendTo[testList,
VerificationTest[
	{kap4=!=0,Simplify[Coefficient[quartic,l63d]/.{l13d->0,l23d->0,g13d->0}]},
	{True,kap4 Sqrt[musqU1]/(4\[Pi] m23d)},
	TestID->"mixing-from-the-temporal-scalar-loop"
]];


(* ::Section:: *)
(*Temporal scalars coupled to a light singlet, with a heavy scalar*)


RepScalar={{{0},"R"},{{1},"C"}};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
f1=CreateInvariant[Group,RepScalar,{{1},{True}}][[1]];
chi2=CreateInvariant[Group,RepScalar,{{2,2},{True,False}}][[1]];
\[Mu]ij=GradMass[mf/2 f1^2+mc chi2]//Simplify//SparseArray;
\[Lambda]3=GradCubic[cf f1 chi2];
\[Lambda]4=GradQuartic[lf/24 f1^4+lc chi2^2]//Simplify;
ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
(*Without the mixed quartics, the identification of the quartics emits messages, which do not affect the masses*)
Quiet[PerformDRhard[]];
PerformDRsoft[{2,3}];

kap3=toGlobal[Normal[private["lamvvsLS"]][[1,1,1]]];
kap3Symbol=First[Cases[kap3,_lamVVSL,{0,Infinity}]];
massLO=Expand[mf3dUS/.PrintScalarMassUS["LO"]];

AppendTo[testList,
VerificationTest[
	{kap3=!=0,Simplify[Coefficient[massLO,kap3Symbol,2]/(kap3/kap3Symbol)^2]},
	{True,-1/(16\[Pi] Sqrt[musqU1])},
	TestID->"bubble-of-temporal-scalars-with-a-heavy-scalar"
]];


(* ::Section:: *)
(*Report*)


report=TestReport[testList];
Print[report];
Print[Column[{#["TestID"],#["Outcome"]}&/@Values[report["TestResults"]]]];
Exit[If[report["AllTestsSucceeded"],0,1]]
