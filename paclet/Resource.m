(* ::Package:: *)

(*
	Landing page of DRalgo in the Wolfram Language Paclet Repository, read by
	paclet.wls. The name, version, description, publisher, primary context,
	license and Wolfram Language version are taken from PacletInfo.m. The
	creators of the paclet are listed in "ContributedBy", all contributors in
	"Contributors", which the last paragraph of the details names. In
	"Details", [text](url) is a hyperlink, `code` is set as code and *text* in
	italics.
*)
<|
	"LongDescription" -> "DRalgo is an algorithmic implementation that constructs an effective, dimensionally reduced, high-temperature field theory for generic models.",
	"Details" -> {
		"DRalgo is an algorithmic implementation that constructs an effective, dimensionally reduced, high-temperature field theory for generic models. The corresponding Mathematica package automatically performs the matching up to next-to-leading order, as described in the public release [arXiv:2205.08815](https://arxiv.org/abs/2205.08815). DRalgo also matches the dimension-5 and dimension-6 operators of the three-dimensional effective theory at one loop, as described in [arXiv:2605.15176](https://arxiv.org/abs/2605.15176).",
		"To create model files, DRalgo uses functions from [GroupMath](https://renatofonseca.net/groupmath). Since GroupMath is an external package, any use of DRalgo's model-creation features should be accompanied by a corresponding citation:",
		"R. M. Fonseca, GroupMath: *A Mathematica package for group theory calculations*, Comput. Phys. Commun. 267 (2021) 108085, [arXiv:2011.01764](https://arxiv.org/abs/2011.01764)."
	},
	"ContributedBy" -> {"Andreas Ekstedt", "Philipp Schicho", "Tuomas V.I. Tenkanen", "Fabio Bernardo"},
	"Contributors" -> {"Andreas Ekstedt", "Philipp Schicho", "Tuomas V.I. Tenkanen", "Fabio Bernardo", "Romain Guillermo Reinle"},
	"Categories" -> {"Symbolic & Numeric Computation"},
	"Keywords" -> {"Dimensional reduction", "Thermal field theory"},
	"Links" -> {
		"Publication of DRalgo 1.0" -> "https://arxiv.org/abs/2205.08815",
		"Publication of the higher-dimensional operators in DRalgo" -> "https://arxiv.org/abs/2605.15176",
		"DRalgo example files" -> "https://github.com/DR-algo/DRalgo/tree/main/examples"
	},
	"Source" -> "https://github.com/DR-algo/DRalgo",
	"HeroImage" -> "FrontEnd/logo.eps"
|>
