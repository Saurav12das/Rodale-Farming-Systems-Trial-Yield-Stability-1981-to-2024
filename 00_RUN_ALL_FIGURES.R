# FST manuscript supplementary: rebuild all figure tables and six main figures and one supplemental figure.
# Put raw yield Excel and PRISM CSV in the parent directory of this code folder.
# Command line: Rscript 00_RUN_ALL_FIGURES.R
# For a figures-only replot using bundled derived CSVs:
#   Rscript 00_RUN_ALL_FIGURES.R --plots-only
script_arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
CODE_DIR <- if(length(script_arg)) dirname(normalizePath(sub("^--file=","",script_arg[1]))) else getwd()
OUTPUT_DIR <- CODE_DIR
ROOT_DIR <- dirname(CODE_DIR)
source(file.path(CODE_DIR,"01_PREPARE_FIGURE_DATA.R"),local=FALSE)
source(file.path(CODE_DIR,"02_MAKE_FIGURES.R"),local=FALSE)
if(!("--plots-only" %in% commandArgs(TRUE))) fst_prepare()
fst_draw_figures()
