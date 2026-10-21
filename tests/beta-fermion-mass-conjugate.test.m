(* ::Package:: *)

If[$InputFileName=="",
	SetDirectory[NotebookDirectory[]],
	SetDirectory[DirectoryName[$InputFileName]]
];
(*Two models are created in this file*)
DRalgo`DRalgo`$GroupMathMultipleModels=True;

DRalgo`DRalgo`$LoadGroupMath=True;
(*DRalgo`DRalgo`$InstallGroupMath=True;*)

Check[
    Get["../Kernel/DRalgo.m"],
    Message[Get::noopen, "DRalgo` at "<>ToString[$UserBaseDirectory]<>"/Applications"];
    Abort[];
]


(* ::Chapter:: *)
(*One-loop beta functions of the fermion masses and of their conjugates*)


(*
	The Lagrangian contains the Weyl-fermion masses m (\[Mu]IJ) and their
	complex conjugates mC = Conjugate[m] (\[Mu]IJC), and likewise Ysff and
	YsffC = Conjugate[Ysff]. Hermiticity requires beta(mC) = Conjugate[beta(m)],
	with Ysff and m replaced by their conjugates. The one-loop result (Machacek
	and Vaughn; Luo, Wang and Xiao, hep-ph/0211440), in units of 1/(16 Pi^2), is

		beta(m) = 1/2 (Transpose[Y2].m + m.Y2) + 2 Y_b.mC.Y_b - 3 g^2 {C2(F), m},
		Y2 = Sum[YsffC_b.Ysff_b],

	where the order of Y2 and Transpose[Y2] follows from covariance under
	unitary rotations of the fermions, as in the beta function of Ysff.

	Up to DRalgo 1.5.0, CounterTerm had an extra overall (-1) in the Yukawa
	term of beta(mC), and the field-renormalisation term of beta(m) in the
	transposed order, 1/2 (Y2.m + m.Transpose[Y2]), which differs when Y2 is
	not real. The counterterms Z\[Mu]IJF and Z\[Mu]IJFC enter the two-loop
	tadpoles.

	Checked here, for two Dirac flavours with complex Yukawa and mass matrices
	(the conjugates are independent symbols yC and mC):
	- beta(mC) is the complex conjugate of beta(m);
	- beta(m) is the formula above;
	- the two-loop tadpole of the real scalar is real;
	and the textbook value beta(m) = 3 m y^2 for one Dirac fermion with a real
	Yukawa coupling, for the mass and its conjugate.
*)


testList={};
private[name_String]:=ToExpression["DRalgo`DRalgo`Private`"<>name];


(* ::Section:: *)
(*Model: two Dirac flavours with complex Yukawa and mass matrices*)


(*U(1) with a neutral real scalar and two Dirac fermions of charge q: Weyl fermions L1, R1, L2, R2*)
Group={"U1"};
RepAdjoint={0};
CouplingName={g};
RepScalar={{{0},"R"}};
RepFermion={{{q},"L"},{{q},"R"},{{q},"L"},{{q},"R"}};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];

\[Phi]=CreateInvariant[Group,RepScalar,{{1},{True}}][[1]];
\[Lambda]1=GradTadpole[\[Sigma] \[Phi]];
\[Mu]ij=GradMass[msq/2 \[Phi]^2]//SparseArray;
\[Lambda]3=GradCubic[\[Gamma]/6 \[Phi]^3];
\[Lambda]4=GradQuartic[\[Lambda]/24 \[Phi]^4];

(*m[f,h] and y[f,h] couple R_h (fermion 2h) and L_f (fermion 2f-1); mC and yC are their conjugates*)
fermionMass[f_,h_]:=GradMassFermion[CreateInvariantFermion[Group,RepFermion,{{2h,2f-1},{False,True}}][[1]]];
fermionMassC[f_,h_]:=GradMassFermion[CreateInvariantFermion[Group,RepFermion,{{2f-1,2h},{False,True}}][[1]]];
yukawa[f_,h_]:=GradYukawa[CreateInvariantYukawa[Group,RepScalar,RepFermion,{{1,2h,2f-1},{True,False,True}}][[1]]];
yukawaC[f_,h_]:=GradYukawa[CreateInvariantYukawa[Group,RepScalar,RepFermion,{{1,2f-1,2h},{True,False,True}}][[1]]];
\[Mu]IJ=Sum[m[f,h] fermionMass[f,h],{f,2},{h,2}]//SparseArray;
\[Mu]IJC=Sum[mC[f,h] fermionMassC[f,h],{f,2},{h,2}]//SparseArray;
Ysff=Sum[y[f,h] yukawa[f,h],{f,2},{h,2}]//SparseArray;
YsffC=Sum[yC[f,h] yukawaC[f,h],{f,2},{h,2}]//SparseArray;

conjugate[expr_]:=expr/.{
	y[f_,h_]:>yC[f,h],yC[f_,h_]:>y[f,h],
	m[f_,h_]:>mC[f,h],mC[f_,h_]:>m[f,h],
	Complex[a_,b_]:>Complex[a,-b]};

ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
(*Only the steps of PerformDRhard that the checks need, to keep the test fast*)
private["TadPole"][];
private["ScalarSelfEnergy"][];
private["CounterTerm"][];
private["TadPole2Loop"][];

AppendTo[testList,
VerificationTest[
	{conjugate[Normal[\[Mu]IJ]]-Normal[\[Mu]IJC],conjugate[Normal[Ysff]]-Normal[YsffC]}//Expand//Flatten//Union,
	{0},
	TestID->"model-masses-and-yukawas-come-with-their-conjugates"
]];


(* ::Section:: *)
(*beta(mC) is the complex conjugate of beta(m)*)


\[Beta]Mass=Normal[private["\[Beta]\[Mu]IJF"]];
\[Beta]MassC=Normal[private["\[Beta]\[Mu]IJFC"]];

AppendTo[testList,
VerificationTest[
	Expand[\[Beta]MassC-conjugate[\[Beta]Mass]]//Flatten//Union,
	{0},
	TestID->"conjugate-mass-beta-is-complex-conjugate-of-mass-beta"
]];


(* ::Section:: *)
(*beta(m) agrees with the general one-loop formula*)


{yP,yCP,mP,mCP,gP}=Normal/@{Ysff,YsffC,\[Mu]IJ,\[Mu]IJC,gvff};
Y2=Sum[yCP[[b]] . yP[[b]],{b,Length[yP]}];
fieldTerm=1/2 (Transpose[Y2] . mP+mP . Y2);
vertexTerm=2 Sum[yP[[b]] . Transpose[mCP] . yP[[b]],{b,Length[yP]}];
(*-3 g^2 {C2(F), m} = 6 g_a m g_a^T for a gauge-invariant mass*)
gaugeTerm=6 Sum[gP[[a]] . mP . Transpose[gP[[a]]],{a,Length[gP]}];

AppendTo[testList,
VerificationTest[
	Expand[\[Beta]Mass-(fieldTerm+vertexTerm+gaugeTerm)/(16\[Pi]^2)]//Flatten//Union,
	{0},
	TestID->"mass-beta-agrees-with-general-one-loop-formula"
]];


(* ::Section:: *)
(*The two-loop tadpole of the real scalar is real*)


(*The counterterms Z\[Mu]IJF and Z\[Mu]IJFC enter the two-loop tadpole TadPoleNLO*)
tadpoleNLO=Normal[private["TadPoleNLO"]];

AppendTo[testList,
VerificationTest[
	{FreeQ[tadpoleNLO,m],Expand[conjugate[tadpoleNLO]-tadpoleNLO]},
	{False,{0}},
	TestID->"two-loop-tadpole-of-real-scalar-is-real"
]];


(* ::Section:: *)
(*Textbook value for one Dirac fermion with a real Yukawa coupling*)


(*The model of ScalYukawa.test.m: beta(m\[Psi]) = 3 m\[Psi] y^2/(16 Pi^2) (e.g. Srednicki, ch. 52), for the mass and its conjugate*)
Group={"U1"};
RepAdjoint={0};
CouplingName={g};
RepScalar={{{0},"R"}};
RepFermion={{{0},"L"},{{0},"R"}};
{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJC,\[Mu]IJ,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];
\[Phi]=CreateInvariant[Group,RepScalar,{{1},{True}}][[1]];
\[Mu]ij=GradMass[msq/2 \[Phi]^2]//SparseArray;
\[Lambda]4=GradQuartic[\[Lambda]/24 \[Phi]^4];
\[Mu]IJ=m\[Psi] GradMassFermion[CreateInvariantFermion[Group,RepFermion,{{2,1},{False,True}}][[1]]];
\[Mu]IJC=m\[Psi] GradMassFermion[CreateInvariantFermion[Group,RepFermion,{{1,2},{False,True}}][[1]]];
Ysff=y GradYukawa[CreateInvariantYukawa[Group,RepScalar,RepFermion,{{1,2,1},{True,False,True}}][[1]]];
YsffC=y GradYukawa[CreateInvariantYukawa[Group,RepScalar,RepFermion,{{1,1,2},{True,False,True}}][[1]]];

ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
private["CounterTerm"][];

AppendTo[testList,
VerificationTest[
	{Normal[private["\[Beta]\[Mu]IJF"]],Normal[private["\[Beta]\[Mu]IJFC"]]}//Expand,
	Table[3 m\[Psi] y^2/(16\[Pi]^2) {{0,1},{1,0}},{2}]//Expand,
	TestID->"scalar-yukawa-mass-beta-is-textbook-value"
]];


(* ::Section:: *)
(*Report*)


report=TestReport[testList];
Print[report];
Print[Column[{#["TestID"],#["Outcome"]}&/@Values[report["TestResults"]]]];
Exit[If[report["AllTestsSucceeded"],0,1]]
