<img src="https://raw.githubusercontent.com/DR-algo/DRalgo/refs/heads/main/FrontEnd/logo.svg" alt="DRalgoLogo" width="100"/>

# DRalgo

**DRalgo** is an algorithmic implementation that constructs an effective,
dimensionally reduced, high-temperature field theory for generic models.
The corresponding Mathematica package automatically performs the matching to next-to-leading order. 
Public release of
[https://arxiv.org/abs/2205.08815](https://arxiv.org/abs/2205.08815).

## Status

[![license: GPLv3](https://img.shields.io/badge/license-GPLv3-brightgreen.svg)](https://github.com/DR-algo/DRalgo/blob/master/LICENSE)
[![Version](https://img.shields.io/github/v/tag/DR-algo/DRalgo?label=Version)](https://github.com/DR-algo/DRalgo/releases/latest/)
![compatibility](https://img.shields.io/badge/Mathematica-12.x_13.x_14.x-brightgreen.svg)

## Installation

### Paclet Installation

**DRalgo** can be installed as a Wolfram Paclet by running one of the following commands in Mathematica:

#### From the Wolfram Repository
```mathematica
PacletInstall["DRalgo/DRalgo"]
```
[Visit the Wolfram Repository](https://resources.wolframcloud.com/PacletRepository/resources/DRalgo/DRalgo/)

#### From the GitHub Repository
```mathematica
PacletInstall["https://github.com/DR-algo/DRalgo/releases/latest/download/DRalgo.paclet"]
```

> **Note:** Ensure that all dependencies of **DRalgo** are installed. Refer to the [Requirements section](#requirements) for details.

### Manual Installation

**DRalgo** can also be installed manually by cloning the repository into the **Applications** folder within either the base or user-specific **Mathematica Applications** directory. These directories store Mathematica packages and can be located by evaluating the variables **`$BaseDirectory`** and **`$UserBaseDirectory`** in a Mathematica session. To identify these directories, run the following commands in Mathematica:

```mathematica
Print["Base Directory: ", FileNameJoin[{$BaseDirectory, "Applications"}]]
Print["User Base Directory: ", FileNameJoin[{$UserBaseDirectory, "Applications"}]]
```

#### Common Installation Paths

For versions of Mathematica prior to 14.1, the typical paths for these directories are as follows:

**Linux**
- `/usr/share/Mathematica/Applications`
- `~/.Mathematica/Applications`

**macOS**
- `~/Library/Mathematica/Applications`

**Windows**
- `C:\ProgramData\Mathematica\Applications`
- `C:\Users\<username>\AppData\Roaming\Mathematica\Applications`

Starting with Mathematica 14.1, these directories have been renamed to include "Wolfram" instead of "Mathematica." For example, on macOS, the path becomes:

- `~/Library/Wolfram/Applications`

Ensure that the **DRalgo** package is placed in the appropriate directory for your system and Mathematica version.

### Requirements

**DRalgo** is written in the *Wolfram Mathematica* language and relies on the external package **GroupMath**, available here:  

- [**Mathematica**](https://www.wolfram.com/mathematica/): versions 12.x, 13.x, and 14.x  
- [**GroupMath**](https://renatofonseca.net/groupmath): version 1.1.3

GroupMath can be installed manually from the link above or automatically by setting the following flag **before** loading DRalgo in Mathematica:
```mathematica
DRalgo`DRalgo`$InstallGroupMath = True;
```
By default, DRalgo will automatically load GroupMath when needed.
If you wish to disable automatic loading, set the following flag instead:
```mathematica
DRalgo`DRalgo`$LoadGroupMath = False;
```

Since GroupMath is an external package, any use of the model-creation features
in DRalgo should be accompanied by a corresponding citation of:
*R. M. Fonseca, GroupMath: A Mathematica package for group theory calculations,
Comput. Phys. Commun. 267 (2021) 108085 [[2011.01764]](https://arxiv.org/abs/2011.01764)*

## Running

Once the **DRalgo** package is installed, it can be loaded in Mathematica by executing the following command:

```mathematica
<<DRalgo`DRalgo`
```

This command initializes the package, making its functions and features available for use within your Mathematica session.

### Running the examples

To explore how **DRalgo** works in practice, we recommend reviewing the provided [examples](https://github.com/DR-algo/DRalgo/tree/main/examples).

These examples can be executed directly within a Mathematica notebook or via [Wolframscript](https://www.wolfram.com/wolframscript/). Wolframscript enables the execution of Wolfram Language scripts without requiring a full Mathematica installation, offering the core computational capabilities of Wolfram Mathematica.

To run an example file using Wolframscript, use the following command:

```bash
$ wolframscript -file examples/xsm.m
```

Ensure that an active `WolframKernel` is available for the above command to work.

## New models
Users are encouraged to contribute their own models to the community by submitting the model file through the [Issue Tracker](https://github.com/DR-algo/DRalgo/issues) on GitHub. Submitted models will undergo verification before being added to the official model repository. 

When submitting a model, please include a reference to a relevant paper or provide the explicit Lagrangian in an accompanying Mathematica notebook to ensure clarity and reproducibility.

## Issues and licensing provisions
For details regarding the licensing and terms of use for this software, please refer to the [LICENSE](LICENSE) file.  

If you encounter any bugs or issues, please report them through the [GitHub Issue Tracker](https://github.com/DR-algo/DRalgo/issues). Your feedback is invaluable in improving the software.
