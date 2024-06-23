# Created on 2019/10/11

########################### DEFAULT PARAMETER VALUES ###########################

ID_SLOT <- "orig.ident"


################################# OPTION MENU ##################################

library("optparse")

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")
cat("\n")
option_list <- list( 
    make_option("--object", type = "character",
        help="[REQUIRED] .Robj file containing the Seurat object"),
    make_option("--out_dir", type = "character",
        help="[REQUIRED] directory where output will be stored"),
    make_option("--id_slot", type = "character",
        help="[REQUIRED] mode ID for clustering (see object@meta.data), to plot the tSNE coloured by cluster")
)


################################ PARSE OPTIONS #################################

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$object)) {
	write("Option --object is required\nTry --help for help", stderr()) 
	q()
} else {
	OBJECT <- opt$object
}

if (is.null(opt$out_dir)) {
	write("Option --out_dir is required\nTry --help for help", stderr()) 
	q()
} else {
	OUT_DIR <- opt$out_dir
}

if (!is.null(opt$id_slot)) 
	ID_SLOT <- opt$id_slot


############################### EXPORT FUNCTIONS ###############################

library("funr")

WORKING_DIR <- getwd()
SCRIPT_PATH <- dirname(sys.script())
SCRIPT_NAME <- basename(sys.script())
setwd(SCRIPT_PATH)
FUNCTIONS <- "functions.R"
source(FUNCTIONS)
setwd(WORKING_DIR)

######################### PRINT THE PARAMETERS TO FILE #########################

sink(file = paste(OUT_DIR, "parameters.txt", sep = "/"), append = TRUE)
cat("======== clusters composition parameters ========\n")
for (i in 1:length(opt))
	cat(paste(names(opt)[i], "=", opt[[i]], "\n", sep=""))
sink()


################################## EXECUTION ###################################

object <- eval(parse(text=load(OBJECT)))

SamplesComposition(object = object, out.dir = OUT_DIR, id_slot = ID_SLOT)	
BarplotSampleComposition(out.dir = OUT_DIR)
BarplotClusterComposition(out.dir = OUT_DIR)


sessionInfo()

q()

