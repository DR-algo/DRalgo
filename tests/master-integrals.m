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
(*Master integrals*)


(*
	The master integrals of MasterIntegrals.m, which HardToSoft.m, SoftToUS.m, EffPot.m, and HEFT.m share.
	The arguments of the 3d integrals are squared masses.

	Checked here:
	- the one-loop integrals against their closed forms and partial fractions, and the bubble and its
	  p^2 term against the integral 1/(4\[Pi] p) ArcTan[p/(m1+m2)];
	- the p^2 terms of the HEFT bubbles with vector lines against their form in squared masses;
	- the sunsets with several propagators on a line against numerical integrals, including equal masses on a line
	  and two massless propagators on a line, whose hard part is the integral with the IR power divergences removed;
	- the two-loop functions of the effective potential, fssv with a massless vector against the sunset of SoftToUS.m;
	- the one-loop sum-integral I1p against HardThermal1LoopInt.
*)


testList={};
tadpole=DRalgo`DRalgo`Private`TadpoleIntegral3d;
bubble=DRalgo`DRalgo`Private`BubbleIntegral3d;
bubbleP2=DRalgo`DRalgo`Private`BubbleIntegral3dP2;
bubbleP2SV=DRalgo`DRalgo`Private`BubbleIntegral3dP2SV;
bubbleP2VV=DRalgo`DRalgo`Private`BubbleIntegral3dP2VV;
triangle=DRalgo`DRalgo`Private`TriangleIntegral3d;
box=DRalgo`DRalgo`Private`BoxIntegral3d;
sunset=DRalgo`DRalgo`Private`SunsetIntegral3d;
positive=x>0&&y>0&&z>0&&w>0;


(* ::Section:: *)
(*One loop*)


(*The bubble with external momentum k, 1/(4\[Pi] k) ArcTan[k/(m1+m2)], at k=0 and its k^2 term*)
bubbleK[k_,x_,y_]:=ArcTan[k/(Sqrt[x]+Sqrt[y])]/(4\[Pi] k);

AppendTo[testList,
VerificationTest[
	Simplify[{
		bubble[x,y]-Limit[bubbleK[k,x,y],k->0],
		bubble[x,x]-Limit[bubbleK[k,x,x],k->0],
		bubble[x,0]-Limit[bubbleK[k,x,0],k->0],
		bubbleP2[x,y]-SeriesCoefficient[bubbleK[k,x,y],{k,0,2}],
		bubbleP2[x,x]-SeriesCoefficient[bubbleK[k,x,x],{k,0,2}],
		bubbleP2[x,0]-SeriesCoefficient[bubbleK[k,x,0],{k,0,2}],
		bubble[0,0],bubbleP2[0,0]},positive],
	{0,0,0,0,0,0,0,0},
	TestID->"bubble-and-its-p2-term"
]];


AppendTo[testList,
VerificationTest[
	N[{bubble[1.3,0.2],
		NIntegrate[4\[Pi] k^2/(2\[Pi])^3 1/((k^2+1.3)(k^2+0.2)),{k,0,Infinity}]}]//Apply[Subtract]//Abs//#<10^-8&,
	True,
	TestID->"bubble-numerically"
]];


(*Partial fractions, 1/((p^2+x)(p^2+y))=(1/(p^2+y)-1/(p^2+x))/(x-y)*)
AppendTo[testList,
VerificationTest[
	Simplify[{
		bubble[x,y]-(tadpole[y]-tadpole[x])/(x-y),
		triangle[x,y,z]-(bubble[y,z]-bubble[x,z])/(x-y),
		triangle[x,y,0]-(bubble[y,0]-bubble[x,0])/(x-y),
		triangle[x,0,0]-(bubble[0,0]-bubble[x,0])/(x-0),
		triangle[x,x,y]-Limit[(bubble[z,y]-bubble[x,y])/(x-z),z->x],
		box[x,y,z,w]-(triangle[y,z,w]-triangle[x,z,w])/(x-y),
		box[x,y,z,0]-(triangle[y,z,0]-triangle[x,z,0])/(x-y),
		box[x,y,0,0]-(triangle[y,0,0]-triangle[x,0,0])/(x-y),
		box[0,x,0,y]-(triangle[x,0,y]-triangle[0,0,y])/(0-x),
		box[x,0,0,0]-(triangle[0,0,0]-triangle[x,0,0])/(x-0),
		triangle[0,0,0],box[0,0,0,0]},positive],
	{0,0,0,0,0,0,0,0,0,0,0,0},
	TestID->"triangle-and-box-from-partial-fractions"
]];


(*The p^2 terms of the HEFT bubbles with vector lines, in the form of squared masses*)
AppendTo[testList,
VerificationTest[
	Simplify[{
		bubbleP2SV[x,y]-(8 (-tadpole[x]+tadpole[y]))/(3 (x-y)),
		bubbleP2SV[x,x]+(4 tadpole[x])/(3 x),
		bubbleP2SV[x,0]+(8 tadpole[x])/(3 x),
		bubbleP2VV[x,y]-(2 ((5 x-y) y^2 tadpole[x]+x^2 (x-5 y) tadpole[y]))/(3 x (x-y)^3 y),
		bubbleP2VV[x,x]-(5 tadpole[x])/(12 x^2),
		bubbleP2SV[0,0],bubbleP2VV[x,0],bubbleP2VV[0,x],bubbleP2VV[0,0]},positive],
	{0,0,0,0,0,0,0,0,0},
	TestID->"heft-bubbles-with-vector-lines"
]];


(* ::Section:: *)
(*Two loops*)


(*
	Numerical sunsets at the scale mu: the q-line and the (p+q)-line give the bubble of the external momentum p,
	and two propagators on the q-line its divided difference. Two massless propagators on a line give the hard part,
	where the IR power divergence \[Integral]1/p^4 vanishes, and for the q-line minus the coefficient of w^2 of the bubble,
	which is even in w. A line with two massless propagators is taken as the p-line, since the q-integral alone
	vanishes for a massless (p+q)-line, while the p-integral then has an IR divergence.
*)
mu=1.7;
bubbleNum[p_,a_,b_]:=If[a+b==0,1/(8 p),ArcTan[p/(a+b)]/(4\[Pi] p)];
qLine[p_,{c_},t_]:=bubbleNum[p,c,t];
qLine[p_,{c1_,c2_},t_]:=Which[
	c1!=c2,(bubbleNum[p,c2,t]-bubbleNum[p,c1,t])/(c1^2-c2^2),
	c1!=0,-1/(2c1) D[bubbleNum[p,cc,t],cc]/.cc->c1,
	True,-t/(4\[Pi] (t^2+p^2)^2)];
sunsetNum[{a__},{c__},t_]:=Module[{f,x,lines={a}},
	f[p_]:=qLine[p,{c},t];
	Which[
		Length[lines]==1,
			Infinity,
		lines[[1]]!=lines[[2]],
			NIntegrate[p^2/(2\[Pi]^2) f[p]/((p^2+lines[[1]]^2)(p^2+lines[[2]]^2)),{p,0,Infinity},WorkingPrecision->20,PrecisionGoal->10],
		lines[[1]]!=0,
			NIntegrate[p^2/(2\[Pi]^2) f[p]/(p^2+lines[[1]]^2)^2,{p,0,Infinity},WorkingPrecision->20,PrecisionGoal->10],
		True,
			With[{f0=Limit[f[x],x->0]},
				NIntegrate[p^2/(2\[Pi]^2) (f[p]-f0)/p^4,{p,0,Infinity},WorkingPrecision->20,PrecisionGoal->10]]]
];

(*Linear masses {line 1}, {line 2}, line 3, as rationals for the numerical integration*)
cases={
	{{1,3/5},{1/2},3/10},
	{{1,1},{1/2},3/10},
	{{1,1},{0},0},
	{{4/5,4/5},{1},0},
	{{0,0},{1/2},3/10},
	{{0,0},{1},0},
	{{1,3/5},{7/10,1/2},3/10},
	{{1,1},{7/10,1/2},3/10},
	{{1,3/5},{1/2,1/2},3/10},
	{{1,3/5},{1/2,1/2},0},
	{{1,1},{1,1},1},
	{{0,0},{1/2,1/2},3/10},
	{{0,0},{1,1},0},
	{{0,0},{0,0},1},
	{{1,0},{0,0},1/2}};

sunsetError[{a_,c_,t_}]:=Abs[N[sunset[a^2,c^2,t^2,mu],20]-sunsetNum[a,c,t]]/Abs[sunsetNum[a,c,t]];

AppendTo[testList,
VerificationTest[
	Max[Quiet[sunsetError/@cases]]<10^-7,
	True,
	TestID->"sunsets-with-several-propagators-on-a-line-numerically"
]];


(*
	Three propagators on the p-line and one on the q-line, as in the two-loop masses with a bubble on a line: with n0 massless
	propagators on the p-line, the Taylor terms of order below p^(2 n0-2) of the rest of the integrand are scaleless and make
	the integral IR divergent, and are dropped.
*)
sunsetNum3[{a__},{c_},t_]:=Module[{x,lines={a},n0,massive,rest,taylor},
	n0=Count[lines,0];
	massive=DeleteCases[lines,0];
	rest[p_]:=qLine[p,{c},t] Times@@(1/(p^2+#^2)&/@massive);
	taylor=If[n0<2,0,Normal[Series[rest[x],{x,0,2n0-4}]]];
	NIntegrate[p^2/(2\[Pi]^2) (rest[p]-(taylor/.x->p))/p^(2n0),{p,0,Infinity},WorkingPrecision->20,PrecisionGoal->10]
];

casesTripled={
	{{1,3/5,1/2},{7/10},3/10},
	{{1,1,3/5},{1/2},3/10},
	{{1,1,1},{1/2},3/10},
	{{0,1,3/5},{1/2},3/10},
	{{0,1,1},{1/2},0},
	{{0,0,1},{1/2},3/10},
	{{0,0,1},{1},0},
	{{0,0,0},{1/2},3/10}};

sunsetError3[{a_,c_,t_}]:=Abs[N[sunset[a^2,c^2,t^2,mu],20]-sunsetNum3[a,c,t]]/Abs[sunsetNum3[a,c,t]];

AppendTo[testList,
VerificationTest[
	Max[Quiet[sunsetError3/@casesTripled]]<10^-7,
	True,
	TestID->"sunsets-with-three-propagators-on-a-line-numerically"
]];


AppendTo[testList,
VerificationTest[
	{sunset[0,0,0,mu],sunset[{0,0},0,0,mu],sunset[{0,0},{0,0},0,mu],sunset[{0,0},{1,1},0,mu]-sunset[{1,1},{0,0},0,mu],
		Simplify[sunset[x,y,z,mu3]-1/(16\[Pi]^2)(Log[mu3/(Sqrt[x]+Sqrt[y]+Sqrt[z])]+1/2)],
		Simplify[sunset[{0,0},{x,0},0,mu3]-1/(32\[Pi]^2 x^2),positive]},
	{0,0,0,0,0,0},
	TestID->"sunset-at-zero-masses-and-symmetric"
]];


(*The scalar-scalar-vector sunset of the effective potential with a massless vector, as in the NLO ultrasoft pressure*)
fssv=DRalgo`DRalgo`Private`TwoLoopFunctions[mu3][[3]];
AppendTo[testList,
VerificationTest[
	Simplify[fssv[x^2,y^2,0]-(4 (x^2+y^2) Log[mu3/(x+y)]+4 x y+x^2+y^2)/(32 \[Pi]^2),positive],
	0,
	TestID->"effective-potential-fssv-with-a-massless-vector"
]];


(* ::Section:: *)
(*Sum-integrals*)


AppendTo[testList,
VerificationTest[
	Module[{I1p=DRalgo`DRalgo`Private`OneLoopSumIntegrals[]["I1p"],eps=DRalgo`DRalgo`Private`eps,
			Lb=DRalgo`DRalgo`Private`Lb,mu=DRalgo`DRalgo`Private`\[CapitalLambda]},
		FullSimplify[
			(Series[HardThermal1LoopInt["B",1,0,3-2eps],{eps,0,1}]//Normal)
			-(I1p/.DRalgo`DRalgo`Private`ReplaceLb/.Lb->Log[mu^2/T^2]-2 Log[4 \[Pi]]+2EulerGamma),
		T>0&&mu>0]],
	0,
	TestID->"one-loop-sum-integral-I1p"
]];


(* ::Section:: *)
(*Report*)


report=TestReport[testList];
Print[report];
Print[Column[{#["TestID"],#["Outcome"]}&/@Values[report["TestResults"]]]];
Exit[If[report["AllTestsSucceeded"],0,1]]
