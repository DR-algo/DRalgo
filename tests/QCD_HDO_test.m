(* ::Package:: *)

Quit[];


SetDirectory[NotebookDirectory[]];
DRalgo`DRalgo`$LoadGroupMath=True;
<<../Kernel/DRalgo.m


(* ::Chapter:: *)
(*QCD*)


(* ::Section::Closed:: *)
(*Model*)


Group={"SU3"};
CouplingName={gs};
RepAdjoint={{1,1}};
RepScalar={};
RepFermionL={{{1,0}},"L"};
RepFermionR={{{1,0}},"R"};
RepFermion={RepFermionL,RepFermionR};


{gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC}=AllocateTensors[Group,RepAdjoint,CouplingName,RepFermion,RepScalar];


(* ::Section::Closed:: *)
(*Dimensional Reduction*)


ImportModelDRalgo[Group,gvvv,gvff,gvss,\[Lambda]1,\[Lambda]3,\[Lambda]4,\[Mu]ij,\[Mu]IJ,\[Mu]IJC,Ysff,YsffC,Verbose->False];
PerformDRhard[]


(* ::Section:: *)
(*Dimension 6 Matching*)


(*
	Here we construct tensors associated to the group structure constants and
	their contractions
*)
Tadj=-I gvvv/gs;
X2=Contract[Tadj,Tadj,{{2,6},{3,5}}];
X3=Contract[Tadj,Tadj,Tadj,{{2,9},{3,5},{6,8}}];
X4=Contract[Tadj,Tadj,Tadj,Tadj,{{2,12},{3,5},{6,8},{9,11}}];
X6=Contract[Tadj,Tadj,Tadj,Tadj,Tadj,Tadj,{{2,18},{3,5},{6,8},{9,11},{12,14},{15,17}}];

(*
	Below we construct the group tensors of the higher dimensional operators,
	having defined \[Alpha][__] our Wilson Coefficients
*)
OT[6,1]=\[Alpha][F^3] X3;
OT[6,3]=SymmetrizeTensor[
	+\[Alpha][A^2F^2,1]X4
	+\[Alpha][A^2F^2,2]Transpose[X4,{1,3,2,4}]
	+\[Alpha][A^2F^2,3]TensorProduct[X2,X2],{{1,2},{3,4}}
	];
OT[6,8]=SymmetrizeTensor[
	+\[Alpha][D^2A^4,1]X4
	+\[Alpha][D^2A^4,2]Transpose[X4,{1,3,2,4}]
	+\[Alpha][D^2A^4,3]TensorProduct[X2,X2],{{1,2},{3,4}}
	];
OT[6,12]=SymmetrizeTensor[
	+\[Alpha][A^6,1]X6
	+\[Alpha][A^6,2]TensorProduct[X4,X2],{{1,2,3,4,5,6}}
	];
OT[6,13]=\[Alpha][D^2F^2]X2;
OT[6,16]=\[Alpha][D^2 A^2 F]X3;
OT[6,17]=\[Alpha][D^4A^2]X2;


TensorListWithoutA6={OT[6,1],OT[6,3],OT[6,8],OT[6,13],OT[6,16],OT[6,17]};
NListWithoutA6={1,3,8,13,16,17};
(*
	WC is the array of the various Wilson Coefficients \[Alpha][...]
*)
WCsWithoutA6=DeleteDuplicates[Cases[TensorListWithoutA6//Normal, \[Alpha][__], \[Infinity]]];

(*
	Dimension6Matching and Dimension5Matching find the values of
	the Wilson coefficients listed in WC,
	d is the number of spatial dimensions,
	Zb and Zf are 1-loop master integrals
*)
solWithoutA6=Dimension6Matching[TensorListWithoutA6,NListWithoutA6,WCsWithoutA6,3][[1]];


(* ::Section:: *)
(*TESTS*)


testList={};


AppendTo[testList,
TestCreate[
	solWithoutA6,
	{
		\[Alpha][F^3]->-(2/45) I (3 gs^3 Sqrt[T] Zb[3,0]-gs^3 Sqrt[T] Zf[3,0]),
		\[Alpha][A^2 F^2,1]->-(2/405) (1971 gs^4 Zb[3,0]-37 gs^4 Zf[3,0]),
		\[Alpha][A^2 F^2,2]->-(2/405) (999 gs^4 Zb[3,0]-23 gs^4 Zf[3,0]),
		\[Alpha][A^2 F^2,3]->0,
		\[Alpha][A^4 D^2,1]->1/81 (-1215 gs^4 T Zb[3,0]+54 gs^4 T xi Zb[3,0]-27 gs^4 T xi^2 Zb[3,0]+8 gs^4 T Zf[3,0]),
		\[Alpha][A^4 D^2,2]->1/81 (-1161 gs^4 T Zb[3,0]-54 gs^4 T xi Zb[3,0]+27 gs^4 T xi^2 Zb[3,0]+40 gs^4 T Zf[3,0]),
		\[Alpha][A^4 D^2,3]->0,\[Alpha][D^2 F^2]->1/180 (333 gs^2 Zb[3,0]-90 gs^2 xi Zb[3,0]-15 gs^2 xi^2 Zb[3,0]-16 gs^2 Zf[3,0]),
		\[Alpha][A^2 D^2 F]->-(1/90) I (-267 gs^3 Sqrt[T] Zb[3,0]+30 gs^3 Sqrt[T] xi Zb[3,0]+15 gs^3 Sqrt[T] xi^2 Zb[3,0]+14 gs^3 Sqrt[T] Zf[3,0]),
		\[Alpha][A^2 D^4]->1/180 (159 gs^2 Zb[3,0]+90 gs^2 xi Zb[3,0]-45 gs^2 xi^2 Zb[3,0]-8 gs^2 Zf[3,0])
		}
]];


AppendTo[testList,
TestCreate[
	ODIM6[1,3][[1,2,3]],
	1/15 gs^3 Sqrt[T] (3 Zb[3,0]-Zf[3,0])
]];

AppendTo[testList,
TestCreate[
	ODIM6[8,3][[1,1,1,1]],
	2/3 gs^4 T (-99 Zb[3,0]+2 Zf[3,0])
]];

AppendTo[testList,
TestCreate[
	HardThermal1LoopInt["B",3,0,3],
	Zeta[3]/(128 \[Pi]^4 T^2)
]];

AppendTo[testList,
TestCreate[
	HardThermal1LoopInt["F",3,0,3],
	(7 Zeta[3])/(128 \[Pi]^4 T^2)
]];


report=TestReport[testList]
report["ResultsDataset"]



