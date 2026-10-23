(* ::Package:: *)

(* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ *)

(* :Title: MasterIntegrals                                                  *)

(*
       This software is covered by the GNU General Public License 3.
       Copyright (C) 2021-2026 Andreas Ekstedt
       Copyright (C) 2021-2026 Philipp Schicho
       Copyright (C) 2021-2026 Tuomas V.I. Tenkanen

*)

(* :Summary:	One-, two-, and three-loop master integrals of the hard modes and of the 3d theories	*)

(* ------------------------------------------------------------------------ *)


(* ::Section:: *)
(*Sum-integrals of the hard modes*)


(*
	Rewrites short-hand constants with Lb and Lf.
*)
ReplaceLb={
	Lbb->2Lb+2(12Log[Glaisher]-EulerGamma),
	Lbbb->Lb/12+1/6(12Log[Glaisher]-EulerGamma),
	LBF->5/2Lb-1/2Lf+2(12Log[Glaisher]-EulerGamma),
	LFF->1/2Lf+3/2Lb+2(12Log[Glaisher]-EulerGamma),
	LFB->Lb+Lf+2(12Log[Glaisher]-EulerGamma),
	Lfff->1/2Lf-3/2Lb-2(12Log[Glaisher]-EulerGamma),
	LbbM->2Lb+2(12Log[Glaisher]-EulerGamma)};


(*
	HardThermal1LoopInt["B", a, b, d]
   Bosonic hard-thermal 1-loop sum-integral I_B(a,b) in d dimensions.
   a : power of the loop momentum in the denominator (propagator power),
   b : power of the spatial loop momentum in the numerator,
   d : spacetime dimension (d->3 for the 3d EFT).
   Evaluates to the closed form in terms of the Riemann Zeta function
   and Gamma functions after performing the bosonic Matsubara frequency sum.
*)
HardThermal1LoopInt["B",a_,b_,d_]:=((\[CapitalLambda]^2 Exp[EulerGamma])/(4\[Pi]))^((3-d)/2) 2 T (2\[Pi] T)^(d-2a+2b)/(4\[Pi])^(d/2) Gamma[a-d/2]/Gamma[a] Zeta[2a-2b-d];

(*
	HardThermal1LoopInt["F", a, b, d]
   Fermionic hard-thermal 1-loop sum-integral I_F(a,b) in d dimensions.
   Related to the bosonic case by the overall factor (2^(2a-2b-d)-1) that
   accounts for the half-integer fermionic Matsubara frequencies.
*)
HardThermal1LoopInt["F",a_,b_,d_]:=(2^(2a-2b-d)-1)HardThermal1LoopInt["B",a,b,d];


(*
	One-loop sum-integrals of the massless hard modes in d=3-2eps, with I for bosonic and IF for fermionic ones.
	The ones with a trailing p are those of the 3-loop pressure and the 2-loop tadpoles.
*)
OneLoopSumIntegrals[]:=<|
	"I1"->T^2/(12 eps)+T^2 Lbbb,
	"I2M2"->T^2/(24 )(-1/eps-  12Lbbb+2),
	"I1p"->T^2/(12 )+eps T^2 Lbbb,
	"I3M2p"->1/(64 \[Pi]^2 epsbp)+1/(32 \[Pi]^2),
	"I2p"->1/(16 \[Pi]^2 epsbp),
	"I4M4p"->1/(128 \[Pi]^2 epsbp)+1/(48 \[Pi]^2),
	"I2M2p"->-T^2/24+T^2 (((-1/24)LbbM - (-1)*1/24 Lb )+1/12)eps,
	"IF1p"->-T^2/(24)+1/24 T^2 Lfff eps,
	"IF2p"->1/(16 (\[Pi]^2) ) (1/eps+Lf),
	"IF3M2p"->1/(64 \[Pi]^2)(1/eps+Lf+2),
	"IF2M2p"->(-T^2/24-T^2 1/48* Lfff)eps+T^2/48,
	"I1F"->-T^2/(24 eps)+1/24 T^2 Lfff|>;


(*
	Two-loop sum-integrals of the massless hard modes, with the shorthand epsb, epsF, ... for their divergences,
	1/epsb=1/eps+Lbb, see ReplaceLb.
*)
TwoLoopSumIntegrals[]:=Module[{\[Kappa]=1/(16 \[Pi]^2)},<|
	"I1I2"->\[Kappa] T^2/(12 epsb),
	"dI1I2"->3\[Kappa] T^2/(12 epsb)-\[Kappa] T^2/6,
	"I211M020"->T^2/(192 \[Pi]^2),
	"I2M2I2"->T^2/(192 \[Pi]^2)-T^2/(384 \[Pi]^2 epsbbM),
	"I3M2I1"->T^2/(768 \[Pi]^2 epsbbM)+T^2/(384 \[Pi]^2),
	"I2I3m1"->T^2/(512 \[Pi]^2 epsb)-T^2/(768 \[Pi]^2),
	"I1FI2F"->-((T^2)/(384 \[Pi]^2))/epsF,
	"I1BI2F"->((T^2)/(192 \[Pi]^2))/epsFB,
	"I1FI2B"->-((T^2)/(384 \[Pi]^2))/epsBF|>
];


(*
	Three-loop sum-integrals of the massless hard modes, given in 9408276 and 9410360.
	The basketball M00 and the fermionic N00, and the sum-integrals of the remaining 3-loop topologies of the pressure.
*)
ThreeLoopSumIntegrals[]:=Module[{prefactor,expand},
	prefactor=T^4*(mu/(4 \[Pi] T))^(6 eps);
	expand[x_]:=Series[x,{eps,0,0}]//Normal;
	<|
	"M00"->expand[1/(24*16 \[Pi]^2) prefactor (1/eps+91/15+8( Derivative[1][Zeta][-1])/Zeta[-1]-2 Derivative[1][Zeta][-3]/Zeta[-3])],
	"LInt"->(T^4 Log[Glaisher])/(48 \[Pi]^2)+(T^4 Log[mu^2])/(768 \[Pi]^2)+T^4/(2304 \[Pi]^2 eps)+(EulerGamma T^4)/(1152 \[Pi]^2)-(T^4 Log[\[Pi] T])/(384 \[Pi]^2)-(T^4 Log[2])/(192 \[Pi]^2),
	"N00"->expand[1/(96*16 \[Pi]^2) prefactor (1/eps+173/30-42/5 Log[2]+8( Derivative[1][Zeta][-1])/Zeta[-1]-2 Derivative[1][Zeta][-3]/Zeta[-3])]//FullSimplify//Expand,
	"MF1M1"->expand[-1/(192*16 \[Pi]^2) prefactor (1/eps+361/60+76/5 Log[2]+6EulerGamma-4( Derivative[1][Zeta][-1])/Zeta[-1]+4Derivative[1][Zeta][-3]/Zeta[-3])],
	"MF00"->expand[-1/(192*16 \[Pi]^2) prefactor (1/eps+179/30-34/5 Log[2]+8( Derivative[1][Zeta][-1])/Zeta[-1]-2Derivative[1][Zeta][-3]/Zeta[-3])],
	"MFM2P2"->expand[-29/(1728*16 \[Pi]^2) prefactor (1/eps+89/29-90/29 Log[2]+48/29EulerGamma+136/29( Derivative[1][Zeta][-1])/Zeta[-1]-10/29Derivative[1][Zeta][-3]/Zeta[-3])],
	"N2M2"->expand[1/(108*16 \[Pi]^2) prefactor (1/eps+35/8+3/2EulerGamma-63/10 Log[2]+5( Derivative[1][Zeta][-1])/Zeta[-1]-1/2Derivative[1][Zeta][-3]/Zeta[-3])],
	"MP2M2"->expand[11/(216*16 \[Pi]^2) prefactor (1/eps+73/22+12/11EulerGamma+64/11( Derivative[1][Zeta][-1])/Zeta[-1]-10/11Derivative[1][Zeta][-3]/Zeta[-3])]
	|>
];


(* ::Section:: *)
(*Integrals of the 3d theories*)


(*
	The arguments are the squared masses x of the propagators 1/(p^2+x), and the integrals are \[Integral]d^dp/(2\[Pi])^d
	in d=3-2eps at their finite part. Light fields are massless, and scaleless integrals vanish.
*)


(* ::Subsection:: *)
(*One loop*)


(*Tadpole, \[Integral]1/(p^2+x)*)
TadpoleIntegral3d[x_]:=-Sqrt[x]/(4 \[Pi]);


(*One-loop effective potential of a real scalar, 1/2\[Integral]Log[p^2+x]*)
OneLoopPotential3d[x_]:=-x^(3/2)/(12 \[Pi]);


(*Bubble at zero external momentum, \[Integral]1/((p^2+x)(p^2+y))*)
BubbleIntegral3d[0,0]=0;
BubbleIntegral3d[x_,y_]:=1/(4 \[Pi] (Sqrt[x]+Sqrt[y]));


(*The k^2 term of the bubble \[Integral]1/((p^2+x)((p+k)^2+y))=ArcTan[k/(Sqrt[x]+Sqrt[y])]/(4\[Pi] k)*)
BubbleIntegral3dP2[0,0]=0;
BubbleIntegral3dP2[x_,y_]:=-1/(12 \[Pi] (Sqrt[x]+Sqrt[y])^3);


(*
	The k^2 terms of the scalar self-energy in HEFT from a bubble with one scalar and one vector line
	and with two vector lines, with the coupling-independent factors of the diagrams. A massless vector line gives no contribution.
*)
BubbleIntegral3dP2SV[0,0]=0;
BubbleIntegral3dP2SV[x_,y_]:=2/(3 \[Pi] (Sqrt[x]+Sqrt[y]));

BubbleIntegral3dP2VV[x_,y_]:=If[x===0||y===0,0,-(x+3 Sqrt[x] Sqrt[y]+y)/(6 \[Pi] Sqrt[x] Sqrt[y] (Sqrt[x]+Sqrt[y])^3)];


(*Triangle at zero external momenta, \[Integral]1/((p^2+x)(p^2+y)(p^2+z))*)
TriangleIntegral3d[x_,y_,z_]:=Module[{m=Sqrt[{x,y,z}]},
	Switch[Count[m,0],
		3,0,
		2,-1/(4 \[Pi] Total[m]^3),
		_,1/(4 \[Pi] Times@@(Plus@@@Subsets[m,{2}]))]
];


(*Box at zero external momenta, \[Integral]1/((p^2+x)(p^2+y)(p^2+z)(p^2+w))*)
BoxIntegral3d[x_,y_,z_,w_]:=Module[{m=Sqrt[{x,y,z,w}],a,b},
	Switch[Count[m,0],
		4,0,
		3,1/(4 \[Pi] Total[m]^5),
		2,{a,b}=DeleteCases[m,0];-(a^2+a b+b^2)/(4 \[Pi] a^3 b^3 (a+b)),
		_,Total[m]/(4 \[Pi] Times@@(Plus@@@Subsets[m,{2}]))]
];


(* ::Subsection:: *)
(*Two loops*)


(*
	Sunset, \[Integral]\[Integral]1/((p^2+x)(q^2+y)((p+q)^2+z)), at the renormalization scale scale. Its divergence is SunsetDivergence3d/eps.
	A list of squared masses as an argument gives several propagators on the same line, e.g. SunsetIntegral3d[{x1,x2},y,z,scale]
	is \[Integral]\[Integral]1/((p^2+x1)(p^2+x2)(q^2+y)((p+q)^2+z)); lists are allowed on the first two lines, with at most two masses each.
*)
SunsetDivergence3d=1/(4 \[Pi])^2 1/4;

SunsetIntegral3d[x_,y_,z_,scale_]/;!MemberQ[{x,y,z},_List]:=sunsetDoubled[{Sqrt[x]},{Sqrt[y]},Sqrt[z],scale];
SunsetIntegral3d[x_List,y_,z_,scale_]:=sunsetDoubled[Sqrt[x],Sqrt[Flatten[{y}]],Sqrt[z],scale];
SunsetIntegral3d[x_,y_List,z_,scale_]:=SunsetIntegral3d[y,x,z,scale];


(*
	The sunset with the masses a and c, and the lists of masses {a1,a2} and {c1,c2} of the propagators on one line.
	Distinct masses on a line are separated by partial fractions, 1/((p^2+a1^2)(p^2+a2^2))=(1/(p^2+a2^2)-1/(p^2+a1^2))/(a1^2-a2^2),
	and equal ones are a derivative with respect to the squared mass. Two massless propagators on a line give the hard
	part of the integral, the coefficient of a^2 in its expansion in the mass a, which is even in a.
*)
sunsetDoubled[{a_},{c_},t_,scale_]:=If[a===0&&c===0&&t===0,0,1/(4 \[Pi])^2 (Log[scale/(a+c+t)]+1/2)];

sunsetDoubled[{c__},{a1_,a2_},t_,scale_]/;Length[{c}]==1:=sunsetDoubled[{a1,a2},{c},t,scale];

sunsetDoubled[{a1_,a2_},{c_},t_,scale_]:=Module[{s=c+t},
	Which[
		a1=!=a2&&s===0&&a2===0,-(Log[scale/a1]+1/2)/(16 \[Pi]^2 a1^2),
		a1=!=a2&&s===0&&a1===0,-(Log[scale/a2]+1/2)/(16 \[Pi]^2 a2^2),
		a1=!=a2,Log[(a1+s)/(a2+s)]/(16 \[Pi]^2 (a1^2-a2^2)),
		a1=!=0,1/(32 \[Pi]^2 a1 (a1+s)),
		s=!=0,-1/(32 \[Pi]^2 s^2),
		True,0]
];

sunsetDoubled[{a1_,a2_},{c1_,c2_},t_,scale_]:=Module[{u},
	Which[
		a1=!=a2,(sunsetDoubled[{c1,c2},{a2},t,scale]-sunsetDoubled[{c1,c2},{a1},t,scale])/(a1^2-a2^2),
		a1=!=0,-D[sunsetDoubled[{c1,c2},{u},t,scale],u]/(2u)/.u->a1,
		c1=!=c2||c1=!=0,sunsetDoubled[{c1,c2},{a1,a2},t,scale],
		t=!=0,3/(32 \[Pi]^2 t^4),
		True,0]
];


(*
	Two-loop functions of the effective potential in Landau gauge, at the renormalization scale scale.
	The notation follows Martin's convention, see arXiv:1808.07615.
*)
TwoLoopFunctions[scale_]:=Module[{A,I2,I2Div,fvvv,f\[Eta]\[Eta]v,fssv,fvvs,fsss,fss,fvs,fvv},
	A=TadpoleIntegral3d;
	I2[x_,y_,z_]:=SunsetIntegral3d[x,y,z,scale];
	I2Div=SunsetDivergence3d;

(*Sunset with three vector bosons*)
	fvvv[x_,y_,z_]:=1/(4 x y z) (1/3 (x^4 (2 I2Div-3 I2[0,0,x])+(x^2-2 x (y+z)+(y-z)^2) (-(3 (x^2+6 x (y+z)+y^2+6 y z+z^2) I2[x,y,z])
							-I2Div (3 (-(8 x (y+z))-8 y z)-2 (x^2+6 x (y+z)+y^2+6 y z+z^2)))+(x-y)^2 (3 (x^2+6 x y+y^2) I2[0,x,y]
							+I2Div (-(2 (x^2+6 x y+y^2))-24 x y))+(x-z)^2 (3 (x^2+6 x z+z^2) I2[0,x,z]+I2Div (-(2 (x^2+6 x z+z^2))-24 x z))
							+y^4 (2 I2Div-3 I2[0,0,y])+(y-z)^2 (3 (y^2+6 y z+z^2) I2[0,y,z]+I2Div (-(2 (y^2+6 y z+z^2))-24 y z))+z^4 (2 I2Div-3 I2[0,0,z])
							+x A[y] A[z] (3 x^2+15 x y+15 x z-15 y^2-26 y z-15 z^2)+z A[x] A[y] (-(15 x^2)-26 x y+15 x z-15 y^2+15 y z+3 z^2)
							+y A[x] A[z] (-(15 x^2)+15 x y-26 x z+3 y^2+15 y z-15 z^2))+40/3 (I2Div x^2 y z+I2Div x y^2 z+I2Div x y z^2));
	(*Special cases where some masses are 0*)
		fvvv[0,0,0]:=0;

		fvvv[0,y_,z_]:=1/(6 y z) (-(6 y^3 I2[0,y,z])+6 y^3 I2[0,0,y]+30 y^2 z I2[0,y,z]-6 z^3 I2[0,y,z]+30 y z^2 I2[0,y,z]
					+6 z^3 I2[0,0,z]-6 y^2 A[y] A[z]-6 z^2 A[y] A[z]-4 y z A[y] A[z]-57 I2Div y^2 z-57 I2Div y z^2);
		fvvv[x_,0,z_]:=fvvv[0,x,z];
		fvvv[x_,y_,0]:=fvvv[0,x,y];
		fvvv[0,0,z_]:=1/4 z (20 I2[0,0,z]-46 I2Div);
		fvvv[0,y_,0]:=fvvv[0,0,y];
		fvvv[x_,0,0]:=fvvv[0,0,x];

(*Ghost-vector sunset*)

		f\[Eta]\[Eta]v[x_,y_,z_]:=((-x+y+z) A[y] A[z]+A[x] (-z A[y]+(x-y+z) A[z])-(x-y)^2 I2[0,x,y]+(y^2+(x-z)^2-2 y (x+z)) I2[x,z,y])/(2 z);
	(*Special cases where vector-mass is zero*)
		f\[Eta]\[Eta]v[x_,y_,0]:=1/2 (2 I2Div (x+y)-2 (A[x] A[y]+(x+y) I2[0,x,y]));
		f\[Eta]\[Eta]v[0,0,0]:=0;
		f\[Eta]\[Eta]v[0,0,z_]:=1/2 z I2[0,0,z];

(*Scalar-scalar-vector sunset*)
	fssv[x_,y_,z_]:=1/z (-((-x+y+z) A[y] A[z])-A[x] (-z A[y]+(x-y+z) A[z])+(x-y)^2 I2[0,x,y]-((x-y)^2-2 (x+y) z+z^2) I2[x,y,z]);
	(*Special case where vector-mass is zero*)
		fssv[x_,y_,0]:=-2 I2Div (x+y)-2 (-A[x] A[y]-(x+y) I2[0,x,y]);
		fssv[0,0,0]:=0;

(*Vector-vector-scalar sunset*)
	fvvs[x_,y_,z_]:=1/(4 x y) (8 I2Div x y-(x+y-z) A[x] A[y]+y A[x] A[z]+x A[y] A[z]
						-z^2 I2[0,0,z]+(x-z)^2 I2[0,x,z]+(y-z)^2 I2[0,y,z]-(x^2+6 x y+(y-z)^2-2 x z) I2[x,y,z]);
	(*Special case where vector-masses are zero*)
		fvvs[0,y_,z_]:=(6 I2Div y+2 (A[y] A[z]-z I2[0,0,z]+(-3 y+z) I2[0,y,z]))/(4 y);
		fvvs[x_,0,z_]:=fvvs[0,x,z];
		fvvs[0,0,z_]:=1/4 (10 I2Div-6 I2[0,0,z]);
		fvvs[0,0,0]:=0;

(*scalar-scalar-scalar sunset*)
	fsss[x_,y_,z_]:=-I2[x,y,z];

(*scalar-scalar bubble*)
	fss[x_,y_]:=A[x]A[y];

(*scalar-vector bubble*)
	fvs[x_,y_]:=2 A[x]A[y];

(*vector-vector bubble*)
	fvv[x_,y_]:=8/3 A[x]A[y];

	Return[{fvvv,f\[Eta]\[Eta]v,fssv,fvvs,fsss,fss,fvs,fvv}];
]
