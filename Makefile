## Makefile for reproducing automatically the results of this repository

# Type 'make help' for a summary of available targets

##
## TARGETS:
##

.PHONY: all 

## all    : (default) run the main analysis, executing HBGWL_analysis.
all: HBGWL_analysis.nbconvert.ipynb


.PHONY : env clean help
## env    : Create a conda environment with all the necessary dependencies.
##          Don't forget to activate it with `source activate ewind`!
env:
	conda env create -f environment.yml

## clean  : Remove auto-generated notebook from execution.
clean:
	@rm *.nbconvert.ipynb

## help   : Provide information about available targets.
help : Makefile
	@sed -n 's/^##//p' $<


# Pattern rule to execute notebooks. Our main notebook has cells that can
# take a long time to run, so we extend the default timeout to 10 minutes.

%.nbconvert.ipynb : %.ipynb
	time jupyter nbconvert --execute --to notebook --ExecutePreprocessor.timeout=600 $<
