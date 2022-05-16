# DRalgo
**DRalgo** is an algorithmic implementation that constructs an effective,
dimensionally reduced, high-temperature field theory for generic models.
The corresponding Mathematica package automatically performs the matching to next-to-leading order. 
Public release of
[https://arxiv.org/abs/22xx.xxxx](https://arxiv.org/abs/22xx.xxxx).

## Installation

To load DRalgo the following commands need to be executed from the root directory:

	SetDirectory[NotebookDirectory[]]; 
	$LoadGroupMath=True;
	<<DRalgo.m

### Dependencies

To create model files, DRalgo uses functions from **GroupMath**
[https://renatofonseca.net/groupmath](https://renatofonseca.net/groupmath). Since GroupMath is an external package, any use of the model-creation features in DRalgo should be accompanied by a corresponding citation of:
*R. M. Fonseca, GroupMath: A Mathematica package for group theory calculations, Comput. Phys. Commun. 267 (2021) 108085 [[2011.01764]](https://arxiv.org/abs/2011.01764)*


## New models

Users are encouraged to make their own models available for the community. This is possible by submitting the model file via the
[Issue Tracker](https://github.com/dr-algo/DRalgo/issues) of Github.
The model will then be verified and added to the model repository. When submitting a model please refer to a paper or explicitly write out the Lagrangian in an accompanying notebook.

## Issues and licensing provisions
Consult "[COPYING](COPYING)" for information about copying and licencing of this software.
Bugs and issues can be reported via the
[Issue Tracker](https://github.com/dr-algo/DRalgo/issues) of Github.# DRalgo
