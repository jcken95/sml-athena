# sml-code
E: j.c.kennedy1@ncl.ac.uk
jcken95.github.io

#### summary of files and data etc ####
* emulator-data - contains the cheap, expensive and HetGP data for the Athena example

* .stan files - stan files used to fit HetGP and SML emulators. hetGP.stan fits a `standard` hetGP emulator, het-SML fits a stochastic multilevel emulator

* GPfunctionOptim.R - used to construct simple emulators
hetGPfunctions.R - used for `full MAP` estimation

* emulator-integrated-means.R - used to fit and validate the Athena emulators
* diagnostics.R - contains some function for hetgp type emulator diagnostics


* psa - conducts a probabalistic sensitivity analysis

* run-first-order.R will run the analysis for both emulators

* <emulator>-<component>.R runs the analysis for the mean OR variance for HetGP OR SML

* toy-example contains MAP estimate code to run the simple sinusoidal example




