#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Title:: *)
(*Release tooling for the DRalgo paclet*)


(* ::Text:: *)
(*This script builds the paclet, evaluates the examples of its landing page, and submits the version in PacletInfo.m to the Wolfram Language Paclet Repository. The landing page is defined by two text files, paclet/Resource.m and paclet/Examples.m; the definition notebook of the Paclet Repository is generated from them in build/ and never edited by hand. The paclet contains PacletInfo.m and the files in Kernel/ and FrontEnd/ that git tracks, so there is no manifest to maintain.*)
(*Run "wolframscript -file paclet.wls help" for the list of commands.*)


Begin["DRalgoPaclet`"];


(* ::Section:: *)
(*Paths and messages*)


$root = DirectoryName[ExpandFileName[$InputFileName]];
$pacletInfoFile = FileNameJoin[{$root, "PacletInfo.m"}];
$kernelFile = FileNameJoin[{$root, "Kernel", "DRalgo.m"}];
$resourceFile = FileNameJoin[{$root, "paclet", "Resource.m"}];
$examplesFile = FileNameJoin[{$root, "paclet", "Examples.m"}];
$buildDirectory = FileNameJoin[{$root, "build"}];
$stagingDirectory = FileNameJoin[{$buildDirectory, "DRalgo"}];
$definitionNotebook = FileNameJoin[{$stagingDirectory, "ResourceDefinition.nb"}];
$examplesNotebook = FileNameJoin[{$buildDirectory, "Examples.nb"}];
$pacletName = "DRalgo/DRalgo";


$usage = "\
Usage: wolframscript -file paclet.wls <command> [options]

Commands:
  version           Compare the version in PacletInfo.m with Kernel/DRalgo.m,
                    the Paclet Repository and the git tags.
  build             Build build/DRalgo__DRalgo-<version>.paclet from PacletInfo.m
                    and the tracked files in Kernel/ and FrontEnd/, and load it
                    in a fresh kernel.
  examples          Evaluate paclet/Examples.m with the staged paclet and write
                    the result to build/Examples.nb.
  notebook          Generate the definition notebook build/DRalgo/ResourceDefinition.nb
                    from paclet/Resource.m and the evaluated examples, e.g. to
                    preview the landing page from Mathematica.
  check [--all]     Run the checks of the Paclet Repository on the generated
                    definition notebook; --all lists every code-inspection hint.
  submit [--dry-run] [--yes]
                    Submit the version in PacletInfo.m to the Paclet Repository.
                    --dry-run stops before the submission, --yes skips the
                    confirmation. Without a cloud login, the publisher token in
                    the environment variable RESOURCE_PUBLISHER_TOKEN is used.
  release [--yes]   Run build, check and submit in this order.
";


say[args___] := WriteString["stdout", args, "\n"];
warn[args___] := WriteString["stderr", "Warning: ", args, "\n"];
die[args___] := (WriteString["stderr", "Error: ", args, "\n"]; Exit[1]);
relative[file_] := StringReplace[ToString[First[Flatten[{file}]]], StartOfString ~~ $root -> ""];


(* ::Section:: *)
(*Versions and git*)


pacletInfo[] := Get[$pacletInfoFile];
localVersion[] := pacletInfo[]["Version"];


kernelVersion[] := First[
	StringCases[ReadString[$kernelFile], "$DRalgoVersion = \"" ~~ v : Except["\""] .. ~~ "\"" :> v],
	Missing["NotFound"]
];


checkVersions[] := If[kernelVersion[] =!= localVersion[],
	die["PacletInfo.m has version ", localVersion[], " but Kernel/DRalgo.m sets $DRalgoVersion = ", kernelVersion[], "."]
];


(* The highest version of DRalgo/DRalgo on the paclet sites, or None. *)
publishedVersion[] := Module[{remote},
	remote = Quiet[PacletFindRemote[$pacletName, UpdatePacletSites -> True]];
	If[MatchQ[remote, {__PacletObject}],
		Last[SortBy[#["Version"] & /@ remote, versionNumbers]],
		None
	]
];


versionNumbers[v_String] := PadRight[ToExpression /@ StringSplit[v, "."], 3];
newerQ[v1_String, v2_String] := Order[versionNumbers[v1], versionNumbers[v2]] === -1;


git[args__String] := Module[{r = RunProcess[{"git", "-C", $root, args}]},
	If[r["ExitCode"] === 0, StringTrim[r["StandardOutput"]], $Failed]
];


warnUncommitted[] := If[git["status", "--porcelain", "--", "Kernel", "FrontEnd", "PacletInfo.m", "paclet"] =!= "",
	warn["Kernel/, FrontEnd/, PacletInfo.m or paclet/ contain uncommitted changes."]
];


warnUntagged[version_] := If[!MemberQ[StringSplit[git["tag", "--points-at", "HEAD"]], "v" <> version],
	warn["HEAD is not tagged v", version, "; the GitHub release is created from this tag."]
];


versionCommand[] := Module[{version = localVersion[], published = publishedVersion[]},
	say["PacletInfo.m:       ", version];
	say["Kernel/DRalgo.m:    ", kernelVersion[]];
	say["Paclet Repository:  ", Replace[published, None -> "not published"]];
	say["Latest git tag:     ", git["describe", "--tags", "--abbrev=0"]];
	If[kernelVersion[] =!= version, warn["PacletInfo.m and Kernel/DRalgo.m disagree on the version."]];
	warnUntagged[version];
	If[StringQ[published] && !newerQ[version, published],
		warn["the version in PacletInfo.m is not newer than the published version ", published, "."]
	];
];


(* ::Section:: *)
(*Build*)


(*
	The staging directory build/DRalgo is the paclet directory: PacletInfo.m and
	the files in Kernel/ and FrontEnd/ that git tracks. Local files such as
	.DS_Store or untracked drafts never enter the paclet.
*)
stage[] := Module[{files},
	files = Append[StringSplit[git["ls-files", "--", "Kernel", "FrontEnd"], "\n"], "PacletInfo.m"];
	files = Select[files, FileExistsQ[FileNameJoin[{$root, #}]] &];
	If[Length[files] < 2, die["git lists no files in Kernel/; run the script inside the git repository."]];
	Quiet[DeleteDirectory[$stagingDirectory, DeleteContents -> True]];
	Scan[copyInto[$stagingDirectory], files];
	files
];


copyInto[stage_][file_] := With[{target = FileNameJoin[{stage, file}]},
	Quiet[CreateDirectory[DirectoryName[target], CreateIntermediateDirectories -> True]];
	CopyFile[FileNameJoin[{$root, file}], target, OverwriteTarget -> True]
];


buildCommand[] := Module[{version, files, result, archive},
	Needs["PacletTools`"];
	checkVersions[];
	warnUncommitted[];
	version = localVersion[];
	files = stage[];
	Quiet[DeleteDirectory[FileNameJoin[{$buildDirectory, "DRalgo__DRalgo"}], DeleteContents -> True]];
	result = PacletTools`PacletBuild[$stagingDirectory, $buildDirectory];
	If[!MatchQ[result, _Success], die["PacletBuild failed: ", ToString[result, InputForm]]];
	archive = result["PacletArchive"];
	say["Built ", relative[archive], " with ", Length[files], " files."];
	loadTest[archive, version];
	archive
];


(* Loads the built archive in a fresh kernel and checks that DRalgo comes from it. *)
loadTest[archive_, version_] := Module[{dir, code, r},
	dir = FileNameJoin[{$buildDirectory, "loadTest"}];
	Quiet[DeleteDirectory[dir, DeleteContents -> True]];
	ExtractArchive[archive, dir];
	dir = First[Select[FileNames["*", dir], DirectoryQ], die["cannot extract ", archive]];
	code = StringJoin[
		"PacletDirectoryLoad[", ToString[dir, InputForm], "];",
		"DRalgo`DRalgo`$InstallGroupMath = True;",
		"Block[{Print}, Needs[\"DRalgo`DRalgo`\"]];",
		"Print[FindFile[\"DRalgo`DRalgo`\"]];"
	];
	r = Quiet[RunProcess[{"wolframscript", "-code", code}]];
	DeleteDirectory[DirectoryName[dir], DeleteContents -> True];
	If[!AssociationQ[r], die["cannot run wolframscript for the load test."]];
	If[r["ExitCode"] =!= 0 || !StringContainsQ[r["StandardOutput"], dir],
		die["the built paclet does not load:\n", r["StandardOutput"], r["StandardError"]]
	];
	say["Loaded DRalgo ", version, " from the archive in a fresh kernel."];
];


(* ::Section:: *)
(*Examples*)


(* ::Subsection:: *)
(*Reading paclet/Examples.m*)


(*
	The front end parses the package file into cells; Code cells become Input
	cells, and text cells take the markup of markdownText, e.g. `code`.
*)
exampleCells[] := Module[{cells},
	cells = UsingFrontEnd[Module[{nbo, expr},
		nbo = NotebookOpen[$examplesFile, Visible -> False];
		expr = NotebookGet[nbo];
		NotebookClose[nbo];
		expr
	]];
	If[Head[cells] =!= Notebook, die["the front end cannot read ", relative[$examplesFile], "."]];
	First[cells] /. {
		Cell[c_, "Code", opts___] :> Cell[c, "Input", Sequence @@ DeleteCases[{opts}, InitializationCell -> _]],
		Cell[s_String, "Text", opts___] :> Cell[TextData[markdownText[s]], "Text", opts]
	}
];


(* DRalgo is loaded from the staging directory, so the examples test the paclet that is submitted. *)
loadStagedDRalgo[] := Module[{loaded},
	PacletDirectoryLoad[$stagingDirectory];
	DRalgo`DRalgo`$InstallGroupMath = True;
	Block[{Print}, Needs["DRalgo`DRalgo`"]];
	loaded = FindFile["DRalgo`DRalgo`"];
	If[!StringQ[loaded] || !StringStartsQ[loaded, $stagingDirectory],
		die["the examples need DRalgo from ", relative[$stagingDirectory], " but load ", loaded, "."]
	];
];


(* ::Subsection:: *)
(*Evaluation*)


(*
	The examples are evaluated cell by cell in this kernel, which has loaded
	nothing else. Each top-level expression of an input cell increments the line
	number, as in a notebook session. Print and Message output is collected in the
	order in which it occurs and stored in Print and Message cells. A headless
	front end cannot insert results with NotebookEvaluate, hence this loop.
*)
$outputStyles = "Output" | "Print" | "Message" | "Echo" | "PrintTemporary";
inputQ[c_] := MatchQ[c, Cell[_BoxData, "Input" | "Code", ___]] && FreeQ[c, Evaluatable -> False, {1}];
outputQ[c_] := MatchQ[c, Cell[_, $outputStyles, ___]];


processCell[Cell[CellGroupData[{first_?inputQ, rest___}, _], ___]] :=
	Sequence @@ Flatten[{runCell[first], processCell /@ Select[{rest}, !outputQ[#] &]}];
processCell[Cell[CellGroupData[cells_List, state_], opts___]] :=
	Cell[CellGroupData[Flatten[{processCell /@ cells}], state], opts];
processCell[c_?inputQ] := runCell[c];
processCell[_?outputQ] := Sequence[];
processCell[c_] := c;


runCell[Cell[BoxData[boxes_], style_, rest___]] := Module[{extra, opts, lines, first, out},
	extra = Cases[{rest}, _String];
	opts = DeleteCases[Cases[{rest}, _Rule | _RuleDelayed], CellLabel -> _];
	lines = If[ListQ[boxes], DeleteCases[boxes, s_String /; StringMatchQ[s, WhitespaceCharacter ...]], {boxes}];
	first = $line + 1;
	out = Flatten[runLine[#, extra] & /@ lines];
	With[{input = Cell[BoxData[boxes], style, Sequence @@ extra, Sequence @@ opts, CellLabel -> "In[" <> ToString[first] <> "]:="]},
		If[out === {}, input, Cell[CellGroupData[Prepend[out, input], Open]]]
	]
];


runLine[box_, extra_] := Module[{held, result, n = ++$line},
	held = Quiet[ToExpression[box, StandardForm, HoldComplete]];
	If[!MatchQ[held, _HoldComplete], die["cannot parse the input ", ToString[box, InputForm]]];
	If[held === HoldComplete[], Return[{}]];
	$log = {};
	result = Internal`HandlerBlock[{"Message", messageHandler[n, extra]},
		Block[{Print = printHandler[n, extra], $Messages = {}},
			CheckAbort[ReleaseHold[held], $Aborted]
		]
	];
	If[result === Null, $log, Append[$log, outputCell[n, extra, result]]]
];


printHandler[n_, extra_][args___] := ($log = Append[$log,
	Cell[BoxData[ToBoxes[If[Length[{args}] === 1, args, SequenceForm[args]], StandardForm]], "Print",
		Sequence @@ extra, CellLabel -> duringLabel[n]]
]; Null);


messageHandler[n_, extra_][Hold[Message[MessageName[s_, t_String], args___], True]] := ($log = Append[$log,
	messageCell[n, extra, HoldComplete[s], t, HoldComplete[args]]
]);
messageHandler[___][_] := Null;


messageCell[n_, extra_, HoldComplete[s_], t_, HoldComplete[args___]] := Module[{template, text},
	template = Which[
		StringQ[MessageName[s, t]], MessageName[s, t],
		StringQ[MessageName[General, t]], MessageName[General, t],
		True, "-- Message text not found --"
	];
	text = ToBoxes[StringForm[template, Sequence @@ (Short[HoldForm[#], 3] & /@ {args})], StandardForm];
	Cell[BoxData[RowBox[{SymbolName[Unevaluated[s]] <> "::" <> t, ":", " ", text}]], "Message", "MSG",
		Sequence @@ extra, CellLabel -> duringLabel[n]]
];


outputCell[n_, extra_, result_] :=
	Cell[BoxData[ToBoxes[result, StandardForm]], "Output", Sequence @@ extra, CellLabel -> "Out[" <> ToString[n] <> "]="];


duringLabel[n_] := "During evaluation of In[" <> ToString[n] <> "]:=";


(* Evaluates paclet/Examples.m once per run; requires the staging directory. *)
evaluatedExamples[] := evaluatedExamples[] = Module[{cells, messages},
	cells = exampleCells[];
	loadStagedDRalgo[];
	$line = 0;
	cells = processCell /@ cells;
	messages = Cases[cells, Cell[BoxData[RowBox[{name_String, ___}]], "Message", ___, CellLabel -> label_, ___] :> label <> " " <> name, Infinity];
	say["Evaluated ", $line, " example inputs with ", Length[messages], " messages."];
	Scan[say["  ", #] &, messages];
	Notebook[cells]
];


examplesCommand[] := (
	stage[];
	writeNotebook[evaluatedExamples[], $examplesNotebook];
	say["Wrote ", relative[$examplesNotebook], "."]
);


(* ::Section:: *)
(*Definition notebook*)


(*
	The definition notebook is the blank template of the Paclet Repository that
	ships with Mathematica, filled with the content of PacletInfo.m,
	paclet/Resource.m and the evaluated examples in the cell forms that the
	definition-notebook tools read. It lies in the staging directory, which is
	therefore its paclet directory: the manifest lists exactly the files of the
	paclet.
*)
ensureNotebook[] := ensureNotebook[] = Module[{info, resource, nb},
	stage[];
	info = pacletInfo[];
	resource = Get[$resourceFile];
	If[!AssociationQ[resource], die["cannot read ", relative[$resourceFile], "."]];
	nb = UsingFrontEnd[
		Needs["DefinitionNotebookClient`"];
		DefinitionNotebookClient`GenerateDefinitionNotebook["Paclet", <||>]
	];
	If[Head[nb] =!= Notebook, die["the template of the definition notebook is not available."]];
	nb = fillTemplate[nb, info, resource, First[evaluatedExamples[]]];
	writeNotebook[nb, $definitionNotebook];
	say["Generated ", relative[$definitionNotebook], "."];
	$definitionNotebook
];


fillTemplate[nb0_, info_, resource_, examples_] := Module[{nb = nb0},
	nb = nb /. Cell[_, "Title", opts___] :> Cell[info["Name"], "Title", opts];
	nb = nb /. c : Cell[_, "Text", opts___] /; MemberQ[tagsOf[c], "Description"] :>
		Cell[info["Description"], "Text", Sequence @@ withoutDefaultTag[{opts}]];
	nb = fillGroup[nb, "HeroImage", {Cell[BoxData[ToBoxes[heroImage[FileNameJoin[{$root, resource["HeroImage"]}]]]], "Input"]}];
	nb = fillGroup[nb, "LongDescription", {Cell[resource["LongDescription"], "Text"]}];
	nb = fillGroup[nb, "Details", {Cell[TextData[Flatten[Riffle[markdownText /@ resource["Details"], "\n\n"]]], "Notes"]}];
	nb = fillGroup[nb, "PrimaryContext", {Cell[info["PrimaryContext"], "Text"]}];
	nb = fillGroup[nb, "Contributed By", {Cell[StringRiffle[resource["ContributedBy"], "\n"], "Text"]}];
	nb = fillGroup[nb, "SourceControlURL", {Cell[resource["Source"], "Text"]}];
	nb = fillGroup[nb, "Keywords", Cell[#, "Item"] & /@ resource["Keywords"]];
	nb = fillGroup[nb, "Links", Cell[TextData[hyperlinkBox[#1, #2]], "Item"] & @@@ resource["Links"]];
	nb = fillGroup[nb, "CompatibilityWolframLanguageVersionRequired", {Cell[info["WolframVersion"], "Text"]}];
	nb = checkBoxes[nb, "Categories-Checkboxes", resource["Categories"]];
	nb = setLicense[nb, info["License"]];
	nb = setExamples[nb, info["PrimaryContext"], examples];
	addCellIDs[nb]
];


(* The definition-notebook tools address cells by their CellID, which new cells lack. *)
addCellIDs[Notebook[cells_, opts___]] := Module[{next},
	next = Max[Append[Cases[cells, (CellID -> id_Integer) :> id, Infinity], 0]] + 1;
	Notebook[cells /. Cell[c_, style_String, o___] /; FreeQ[{o}, CellID] :> Cell[c, style, o, CellID -> next++], opts]
];


(*
	The Paclet Repository expects a headline image of 400 to 1500 pixels with an
	aspect ratio between 0.75 and 1.25; the logo is rendered 600 pixels high and
	padded to a white square.
*)
heroImage[file_String] := Module[{logo = Import[file], image},
	image = Switch[logo,
		_Graphics, Rasterize[Show[logo, ImageSize -> {Automatic, 600}], Background -> White],
		_Image, ImageResize[logo, {Automatic, 600}],
		_, die["cannot import the headline image ", relative[file], "."]
	];
	ImageCrop[RemoveAlphaChannel[image, White], {600, 600}, Padding -> White]
];


tagsOf[Cell[___, CellTags -> tags_, ___]] := Flatten[{tags}];
tagsOf[_] := {};
withoutDefaultTag[opts_List] := Replace[opts, (CellTags -> tags_) :> (CellTags -> DeleteCases[Flatten[{tags}], "DefaultContent"]), {1}];
headerQ[name_][c_] := MatchQ[c, Cell[___, TaggingRules -> rules_, ___] /; MemberQ[Flatten[{rules}], "TemplateGroupName" -> name]];


(* The placeholder cells of a template section are replaced by the given cells. *)
fillGroup[nb_, name_, cells_List] := Module[{done = False, result},
	result = nb /. Cell[CellGroupData[{h_?(headerQ[name]), rest___}, state_], opts___] :> (done = True;
		Cell[CellGroupData[Join[{h}, replaceDefaults[{rest}, cells]], state], opts]);
	If[!done, die["the template of the definition notebook has no section ", name, "."]];
	result
];


replaceDefaults[rest_, cells_] := Module[{pos = FirstPosition[rest, c_ /; MemberQ[tagsOf[c], "DefaultContent"], None, {1}, Heads -> False]},
	If[pos === None,
		Join[rest, cells],
		Join[Take[rest, pos[[1]] - 1], cells, Select[Drop[rest, pos[[1]] - 1], !MemberQ[tagsOf[#], "DefaultContent"] &]]
	]
];


hyperlinkBox[label_, url_] := ButtonBox[label, BaseStyle -> "Hyperlink", ButtonData -> {URL[url], None}, ButtonNote -> url];


(* [text](url) becomes a hyperlink, `code` is set as code and *text* in italics. *)
markdownText[s_String] := StringSplit[s, {
	"[" ~~ t : Except["]"] .. ~~ "](" ~~ u : Except[")"] .. ~~ ")" :> hyperlinkBox[t, u],
	"`" ~~ t : Except["`"] .. ~~ "`" :> StyleBox[t, "Input"],
	"*" ~~ t : Except["*"] .. ~~ "*" :> StyleBox[t, FontSlant -> "Italic"]
}];


(* A checked box of a checkbox section carries its label as value. *)
checkBoxes[nb_, tag_, labels_List] := Module[{known},
	known = Cases[nb, c : Cell[__] /; MemberQ[tagsOf[c], tag] :> Cases[c, CheckboxBox[_, {_, label_String}] :> label, Infinity], Infinity];
	If[!SubsetQ[Flatten[known], labels], die["unknown entries for ", tag, ": ", StringRiffle[Complement[labels, Flatten[known]], ", "]]];
	nb /. c : Cell[__] /; MemberQ[tagsOf[c], tag] :>
		(c /. CheckboxBox[_, {f_, label_String}] :> CheckboxBox[If[MemberQ[labels, label], label, f], {f, label}])
];


(*
	The License section is a set of radio buttons; the license of PacletInfo.m is
	selected as "Other" by its SPDX name, as the notebook stores a license chosen
	by hand.
*)
setLicense[nb_Notebook, id_String] := Module[{data, name, result},
	Needs["ResourceSystemClient`"];
	data = Quiet[ResourceSystemClient`Private`validateLicensingInformation[id]];
	name = If[MatchQ[data, KeyValuePattern["Name" -> _String]], data["Name"], id];
	result = nb /. Cell[c_, style_, a___, TaggingRules -> rules_, b___, CellTags -> tags_List, d___] /; MemberQ[tags, "LicensingInformation-RadioButtons"] :>
		Cell[c, style, a, TaggingRules -> Join[
			DeleteCases[Flatten[{rules}], ("RadioButtonValue" | "SpecifiedLicenseID") -> _],
			{"RadioButtonValue" -> "Other", "SpecifiedLicenseID" -> name}
		], b, CellTags -> tags, d];
	If[FreeQ[result, "SpecifiedLicenseID" -> name], die["the template of the definition notebook has no license buttons."]];
	result
];


(*
	The Examples section keeps the excluded initialization cells, which load the
	staged paclet when the notebook is evaluated in Mathematica, and replaces the
	sample examples of the template.
*)
setExamples[nb_, context_, examples_List] := Module[{done = False, init, result},
	init = {
		Cell[BoxData[RowBox[{RowBox[{"PacletDirectoryLoad", "[", RowBox[{"NotebookDirectory", "[", "]"}], "]"}], ";"}]],
			"Input", "Excluded", InitializationCell -> True],
		Cell[BoxData[{
			RowBox[{RowBox[{"DRalgo`DRalgo`$InstallGroupMath", "=", "True"}], ";"}], "\n",
			RowBox[{RowBox[{"Needs", "[", "\"" <> context <> "\"", "]"}], ";"}]}],
			"Input", "Excluded", InitializationCell -> True]
	};
	result = nb /. Cell[CellGroupData[{h_?(headerQ["ExampleNotebook"]), rest___}, state_], opts___] :> (done = True;
		Cell[CellGroupData[Join[
			{h},
			Cases[{rest}, Cell[CellGroupData[{ih_?(headerQ["ExampleInitialization"]), ___}, is_], io___] :> Cell[CellGroupData[Prepend[init, ih], is], io]],
			examples
		], state], opts]);
	If[!done, die["the template of the definition notebook has no Examples section."]];
	result
];


(*
	The front end serializes the notebook; its header and outline cache are
	dropped.
*)
writeNotebook[nb_Notebook, file_String] := Module[{string, a, b, stream},
	string = UsingFrontEnd[First[FrontEndExecute[FrontEnd`ExportPacket[nb, "NotebookString"]]]];
	If[!StringQ[string], die["the front end could not serialize ", relative[file], "."]];
	a = StringPosition[string, "(* Beginning of Notebook Content *)\n", 1];
	b = StringPosition[string, "(* End of Notebook Content *)", 1];
	If[a === {} || b === {}, die["unexpected notebook format from the front end."]];
	string = StringTake[string, {a[[1, 2]] + 1, b[[1, 1]] - 1}] <> "\n";
	Quiet[CreateDirectory[DirectoryName[file], CreateIntermediateDirectories -> True]];
	WithCleanup[
		stream = OpenWrite[file, CharacterEncoding -> "ASCII"],
		WriteString[stream, string],
		Close[stream]
	];
];


notebookCommand[] := say["Open ", ensureNotebook[], " in Mathematica to preview or deploy the landing page."];


(* ::Section:: *)
(*Checks of the Paclet Repository*)


checkCommand[all_] := Module[{file, hints, nb, files, inspection, shown, counts},
	file = ensureNotebook[];
	{hints, nb} = UsingFrontEnd[
		Needs["DefinitionNotebookClient`"];
		Block[{$PublisherID = pacletInfo[]["PublisherID"]},
			Module[{nbo, h, n},
				nbo = NotebookOpen[file, Visible -> False];
				DefinitionNotebookClient`CheckDefinitionNotebook[nbo];
				h = DefinitionNotebookClient`HintData["Paclet", None, {"Tag", "Level", "Message", "CellID"}];
				n = NotebookGet[nbo];
				NotebookClose[nbo];
				{h, n}
			]
		]
	];
	If[!ListQ[hints], die["the check of the definition notebook failed."]];
	files = manifestFiles[nb];
	inspection = StringStartsQ[#Tag, "CodeInspection"] &;
	shown = Select[hints, all || !inspection[#] || #Level =!= "Suggestion" &];
	Scan[say[
		StringPadRight[#Level, 11],
		StringPadRight[StringDelete[#Tag, "CodeInspectionFileIssue/"], 26], " ", hintText[#Message],
		If[inspection[#], " [" <> Lookup[files, #CellID, "?"] <> "]", ""]
	] &, SortBy[shown, {Position[{"Error", "Warning", "Suggestion"}, #Level] &}]];
	counts = Counts[StringDelete[#Tag, "CodeInspectionFileIssue/"] & /@ Complement[hints, shown]];
	If[counts =!= <||>, say["Code-inspection suggestions: ", ToString[Normal[counts]], "; list them with --all."]];
	Select[hints, #Level === "Error" &]
];


(* Paths of the files in the manifest of the definition notebook, keyed by CellID. *)
manifestFiles[nb_] := Association[Cases[nb,
	Cell[content_, "FileManagerFile" | "FileManagerDirectory", ___, CellID -> id_, ___] :>
		(id -> FirstCase[content, HoldPattern[FileNameJoin[{_, path___String}]] :> StringRiffle[{path}, "/"], "", Infinity]),
	Infinity]];


hintText[s_String] := StringReplace[s, "\n" -> " "];
hintText[other_] := ToString[other, OutputForm, PageWidth -> Infinity];


(* ::Section:: *)
(*Submission*)


(*
	The generated definition notebook builds the paclet, collects the resource
	data and submits them with ResourceSubmit, as its Submit button does. The
	publisher ID of PacletInfo.m is used, since a cloud account may own several.
	Without a cloud login, ResourceSystemClient takes the publisher token from the
	environment variable RESOURCE_PUBLISHER_TOKEN, as in a GitHub workflow.
*)
submitCommand[dryRun_, yes_] := Module[{version, published, publisher, token, file, result},
	checkVersions[];
	version = localVersion[];
	publisher = pacletInfo[]["PublisherID"];
	published = publishedVersion[];
	If[StringQ[published] && !newerQ[version, published],
		die["PacletInfo.m has version ", version, ", which is not newer than the published version ", published, "."]
	];
	warnUncommitted[];
	warnUntagged[version];
	Needs["ResourceSystemClient`"];
	token = ResourceSystemClient`$PublisherToken;
	If[!TrueQ[$CloudConnected] && !MatchQ[token, _String | _ByteArray],
		die["neither a cloud login nor a publisher token in RESOURCE_PUBLISHER_TOKEN is available."]
	];
	file = ensureNotebook[];
	say["Submitting ", $pacletName, " ", version, " with publisher ID ", publisher, " as ",
		If[TrueQ[$CloudConnected], $CloudUserID, "the holder of the publisher token"],
		" (published: ", Replace[published, None -> "none"], ")."];
	result = UsingFrontEnd[
		Needs["DefinitionNotebookClient`"];
		Needs["PacletResource`DefinitionNotebook`"];
		Block[{$PublisherID = publisher},
			Module[{nbo, built, ro, examples, submission},
				nbo = NotebookOpen[file, Visible -> False];
				(* The resource refers to notebooks of this session; they stay open until the submission. *)
				WithCleanup[
					built = PacletResource`DefinitionNotebook`BuildPaclet[nbo, Interactive -> False];
					If[!MatchQ[built, KeyValuePattern["Result" -> _Success]],
						die["the definition notebook could not build the paclet: ", ToString[built, InputForm]]
					];
					say["Built ", relative[built["Result"]["PacletArchive"]], " through the definition notebook."];
					ro = DefinitionNotebookClient`ScrapeResource["Paclet", nbo,
						Interactive -> False, "ClickedButton" -> "Submit", "FailureCondition" -> None
					];
					If[Head[ro] =!= ResourceObject, die["the definition notebook did not produce a resource: ", ToString[ro, InputForm]]];
					examples = Replace[Quiet[ro["ExampleNotebook"]], nb_NotebookObject :> NotebookGet[nb]];
					say["Resource ", ro["Name"], " with ", Length[Cases[examples, Cell[_, "Input", ___], Infinity]], " example inputs."];
					submission = Which[
						dryRun, say["Dry run: the paclet was not submitted."]; Null,
						!yes && !confirmedQ["Submit " <> $pacletName <> " " <> version <> " to the Paclet Repository?"], say["Not submitted."]; Null,
						True, ResourceSubmit[ro, "Version" -> version]
					];
					If[submission =!= Null,
						say[ToString[DefinitionNotebookClient`SubmitResultMessage["Paclet", submission], OutputForm]]
					];
					submission,
					NotebookClose[nbo]
				]
			]
		]
	];
	If[FailureQ[result], die["the submission failed: ", ToString[result, InputForm]]];
	result
];


confirmedQ[question_] := Module[{answer},
	WriteString["stdout", question, " [y/N] "];
	answer = Quiet[ReadLine["stdin"]];
	StringQ[answer] && StringMatchQ[StringTrim[answer], "y" | "yes", IgnoreCase -> True]
];


(* ::Section:: *)
(*Command line*)


checkOrDie[all_] := If[checkCommand[all] =!= {}, die["the definition notebook has errors."]];


main[{"version"}] := versionCommand[];
main[{"build"}] := buildCommand[];
main[{"examples"}] := examplesCommand[];
main[{"notebook"}] := notebookCommand[];
main[{"check", opts___}] := checkOrDie[MemberQ[{opts}, "--all"]];
main[{"submit", opts___}] := submitCommand[MemberQ[{opts}, "--dry-run"], MemberQ[{opts}, "--yes"]];
main[{"release", opts___}] := (
	buildCommand[];
	checkOrDie[False];
	submitCommand[False, MemberQ[{opts}, "--yes"]]
);
main[{"help" | "--help" | "-h"} | {}] := say[$usage];
main[args_] := (WriteString["stderr", "Unknown arguments: ", StringRiffle[args, " "], "\n\n", $usage]; Exit[2]);


End[];


DRalgoPaclet`main[Rest[$ScriptCommandLine]];
Exit[0];
