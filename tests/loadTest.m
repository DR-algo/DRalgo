(* ::Package:: *)

Quit[];


(* ::Section:: *)
(*Testing package loading*)


If[$InputFileName=="",
	SetDirectory[NotebookDirectory[]],
	SetDirectory[DirectoryName[$InputFileName]]
];


(* ::Subsection:: *)
(*Test direct paclet installation*)


PacletUninstall["DRalgo/DRalgo"]
(*PacletInstall["../build/DRalgo__DRalgo-1.3.0.paclet", ForceVersionInstall -> True]*)
PacletInstall["https://github.com/DR-algo/DRalgo/releases/latest/download/DRalgo.paclet", ForceVersionInstall -> True]


(* Delte GroupMath first *)
GroupMathPath = FileNameJoin[{$UserBaseDirectory, "Applications", "GroupMath"}];
If[FileExistsQ[GroupMathPath],
  DeleteDirectory[FileNameJoin[{$UserBaseDirectory, "Applications", "GroupMath"}], DeleteContents -> True];
  Print["Deleted GroupMath"]
];

DRalgo`DRalgo`$LoadGroupMath=True;
DRalgo`DRalgo`$InstallGroupMath=True;

Check[
    Get["DRalgo`DRalgo`"],
    Message[Get::noopen, "DRalgo` at "<>ToString[$UserBaseDirectory]<>"/Applications"];
    Abort[];
]


(* ::Subsection:: *)
(*Test manual installation*)


PacletUninstall["DRalgo/DRalgo"]


(* Delte GroupMath first *)
GroupMathPath = FileNameJoin[{$UserBaseDirectory, "Applications", "GroupMath"}];
If[FileExistsQ[GroupMathPath],
  DeleteDirectory[FileNameJoin[{$UserBaseDirectory, "Applications", "GroupMath"}], DeleteContents -> True];
  Print["Deleted GroupMath"]
];

(*Put this if you want to create multiple model-files with the same kernel*)
(*DRalgo`$GroupMathMultipleModels=True;*)

DRalgo`DRalgo`$LoadGroupMath=True;
DRalgo`DRalgo`$InstallGroupMath=True;

Check[
    Get["../Kernel/DRalgo.m"],
    Message[Get::noopen, "DRalgo` at "<>ToString[$UserBaseDirectory]<>"/Applications"];
    Abort[];
]

