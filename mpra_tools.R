# mpra tools
# Stephen P. Plassmeyer
# updated November 18, 2020
# Tony's Notes 25Nov2020 - added script to auto install packages and load
# if they aren't found, otherwise load - modified from Stephen's install script
# On volcano plot, changed to geom_text_repel to read easier, added subtitle,
# included an option to color points by density of surrounding points
# rather than sigificance arguments - see function for details,
# centered both title and subtitle
# Added enrichment plot based on logFC 07Dec2020 - Stephen's code
##########################################################################################################################################################
# To Do's
# * get all the functions on board with the proper data nomenclature
# * get_count_threshold: give flexibility to use alternative hypothesis testing functions. Also, consider possible issue of libraries of varying depth
#     - i.e. might have to use different thresholds for libraries with different depth if DNA and RNA are sequenced at significant different depths
##########################################################################################################################################################

packs <- c("R.utils","data.table","plyr","tidyverse","dplyr","lme4",
           "lmerTest","lattice","VennDiagram",
           "patchwork","progress","BiocManager","ggrepel","MASS","ggpubr",
           "nortest","ggrepel","rstatix","psych","gglorenz","ggVennDiagram",
           "remotes","latex2exp","ggpointdensity","umap")
Bioc.packs <- c("edgeR","MPRAnalyze","mpra","M3C")

package.check <- lapply(
  packs,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

Bioc.package.check <- lapply(
  Bioc.packs,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      BiocManager::install(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)
rm(Bioc.packs,Bioc.package.check,packs,package.check)

rm()
remotes::install_github("r-link/corrmorant")
library(corrmorant)
theme_set(theme_bw(base_size = 18))

##########################################################################################################################################################
# Expected Format for data.table-type Data Matrices
# 
# The functions provided within this r script use variable nomenclature in a certain way to refer to the columns and values within a data matrix.
# There are a number expected columns, or fields, that are used by a number of functions, although many should provide users with the flexibility 
# to use alternate fields (e.g. taking log-ratios across experiments rather than conditions within an experiment).
# 
# The nomenclature for the data table structure is as follows:
#   * Columns are called 'fields'
#   * Rows are called 'records' -- each record is a measurement for a single barcode that came from one condition in one replicate in one experiment.
#   * Values may be referred to as 'values' if they are a measurement, such as counts or a transformed metric. Otherwise, values for categorical 
#     data fields will be referred to as 'levels' of that field.
#
# The standard format and expected count matrix is as follows:
# 
#        Field        |                               Description                              |               Example Levels     
# --------------------+------------------------------------------------------------------------+------------------------------------------
#                     |      String indicating what experiment the counts originated from.     |
#      Experiment     |   If libraries from a single experiment were sequenced multiple time,  |            "Vglut_410_NextSeq"
#                     |    for example once as a spike-in and again as a full lane, unique     |             "SY5Y_455_NextSeq"
#                     |       experiment names are used to distinguish the measurements.       |
# --------------------+------------------------------------------------------------------------+-------------------------------------------
#       Element       | String value, interpreted as factor, providing a unique identifier for |
#                     |   each library element. Note that for allelic comparison designs,      |                'Rbfox2_wgs'
#                     |  all variants for a particular element share the same 'Element' value  | Variant;chr3:123954485;CA|C;Family=14179
#                     |  and are distinguished as unique library elements by the Allele field  |                                                        
# --------------------+------------------------------------------------------------------------+-------------------------------------------
#       Allele        | String value, interpreted as factor, denoting the allele type for a    |                ref, alt, shuf      
#                     |       a particular Element for allelic comparison studies.             |  
# --------------------+------------------------------------------------------------------------+-------------------------------------------
#         BC          | Integer value, interpreted as a factor, denoting which unique Barcode  |                1, 2, 3, 4, 5, 6
#                     |        a measurement came from for a particular Element.               |
# --------------------+------------------------------------------------------------------------+-------------------------------------------
#                     | Integer value, interpreted as a factor, signifying a unique biological |
#                     |  replicate from which the measurement was taken. Typically, if a DNA   |               For RNA Libraries:
#       BiolRep       |  library for a viral or maxi pool is used to normalize counts across   |                    1, 2, 3, ...
#                     |   multiple RNA libraries, the DNA replicate is coded as '0' while the  |        For maxi or viral DNA libraries:
#                     |       RNA levels are coded with a positive integer starting at 1.      |                         0        
# --------------------+------------------------------------------------------------------------+-------------------------------------------
#                     |  Integer value, intepreted as a factor, signifying a unique technical  |
#       TechRep       |  replicate for the given measurement. All samples are coded with 1 by  |           1 by default, then 2, 3, ...
#                     |       default, even if technical replicates were not prepared.         |
# --------------------+------------------------------------------------------------------------+-------------------------------------------
#                     |   Short string value denoting the nature of the sequencing library,    |              Strings formated as 
#       Fraction      | whether it's an RNA-seq or DNA-seq library and the source material for |               <rna/dna>.<source>
#                     |  that library, whether input vs TRAP IP RNA or plasmid vs maxi DNA     |        Example: 'rna.trap','dna.maxi'
# --------------------+------------------------------------------------------------------------+-------------------------------------------
#        Counts       |        Integer value for the number of counts for each record.         |
# --------------------+------------------------------------------------------------------------+-------------------------------------------
#         CPM         | Counts per million normalized values within each sequencing library    |
#

##########################################################################################################################################################
# Format Data
# This function takes in a matrix of counts as output by a barcode counting script, where each row is a particular barcode with associated element
# name information for that barcode. Each column is a separate sequencing index with count values for each barcode down each row.
# A sample info table is provided with a specified column of values that match the column labels in the count matrix
# The counts matrix is expanded and applies the specified sample info to genarate a fully annotated count data.table which is returned to the user
# All added sample information is appended as factors to the Count/CPM matrix

# To-Do's: 
# * Make the number of barcode labels generalizeable to cases where not all elements have same number of barcodes in the library design
# * Allow this to be run without splitting a column into element\allele if allele is submitted in the original count matrix as a separate field
# * Build in flexibility to select subset of fields to return and in what order.

format.counts<-function(counts, # count matrix of the format that each field has a sample with per barcode counts, each row corresponds to a specific barcode
                        sample_info, # matrix with a SampleID field matching the field labels in the count matrix containing additional metadata for each sample
                        sample.id = 'SampleID', # designates which field in the sample_info table corresponds to the field labels in the counts matrix
                        element.name = 'Element', # designates which field in the counts matrix contains information about the barcode/element
                        split.element.into = c('Element','Allele'), # specifies what fields to split the element into, useful for pulling out the allele info if not explicitly stated in count matrix 
                        split.char = '_', # character to split element name on to separate different information fields
                        n.bc = 6, # number of barcodes asigned to each element
                        experiment.label = NA, # provides an experiment specific label
                        return.fields = c('Experiment', 'Element', 'Allele', 'BC', 'BiolRep', 'TechRep', 'Fraction', 'Counts', 'CPM')) # what fields to return and in what order
{ 
  
  # Get rid of unneeded fields in the count matrix
  counts<-subset(counts, select = c(element.name, unlist(sample_info[,..sample.id])))%>%as.data.table() # remove columns that aren't the element name or the counts
  
  # sort counts numerically, then alphabetically by element name -- groups all elements together.
  # all samples will be kept in-line at this point so there shouldn't be a concern of swapping of barcodes for the same element between different samples
  counts<-counts[order(counts[,..element.name]),]
  
  counts[,BC := as.factor(rep(1:n.bc, nrow(counts)/n.bc))] # label barcodes as a number
  
  # if labels for splitting the element name in the count matrix are provided, then split data fields accordingly
  if (!is.na(element.name) & length(split.element.into) > 1){
    counts<-separate(counts, c(element.name), split.element.into, sep = split.char)%>%as.data.table()}
  
  data<-c() # dummy variable for holding matrix as we go
  for (s in unlist(sample_info[,..sample.id])){ # for each sample in the sample info, filter down to just those counts in count matrix, add meta data and append to 'data' table
    # if a sample in the count matrix is not listed the sample_info table, it does not get added to the final output data table
    tmp<-subset(counts, select = c(split.element.into, 'BC', s))%>%as.data.table() # filter down to just the element info and the counts for that one sample
    tmp<-setnames(tmp, old = s, new = 'Counts') # rename the sample-specific counts column as 'Counts'
    tmp[,CPM := Counts / sum(Counts) * 1e6] # compute CPM for that sample
    for (col in colnames(sample_info)[colnames(sample_info) != sample.id]){ # for each additional field provided for that sample in the 
      tmp[, c(col) := factor(rep(filter(sample_info, sample_info[,get(sample.id)] == s)[[col]], nrow(tmp)))] # all added sample info will be made into factors
    }
    data<-rbind(data, tmp)
  }
  
  if (!is.na(experiment.label)){data[,Experiment:=rep(experiment.label, nrow(data))]} # add on an experiment label for the entire data set.
  # note that if multiple experiments are represented within a single count matrix, the Experiment value can be specified in the sample_info and left NA in the arguments
  
  data<-subset(data, select = return.fields)%>%as.data.table() # subset and reorder the columsn to be returned
  return(data)
}

##########################################################################################################################################################
# Aggregate Counts
# This is a post-formatting, pre-expression calculation step that allows counts from various levels of replication (Technical, Biological, etc.) to 
# be aggregated prior to any form of normalization, averaging, or log transformation.

aggregate.counts<-function(data, # formatted Count/CPM table output from the format.counts function
                           by.field, # column in the count table on which to filter for the samples to aggregate counts from
                           agg.level, # the level within the specified column that will be filtered (ex. dna.maxi designated in the Fraction field)
                           new.level = NA, # what the new aggregated counts will be labeled in the returned counts table
                           rep.fields, # other fields for which unique information is removed to facilitate aggregation. i.e. label for technical replicates
                           new.rep.levels, # ordered list corresponding with the values specified in rep.field for the values that will be overwritten for all samples
                           sample.fields, # vector of fields which after aggregation will uniquely identify samples for CPM normalization (ex. Fraction and BiolRep)
                           remove.old = F){ # Optionally, the source counts can be removed from the 
  
  # subset to the counts to be aggregated.
  tmp<-filter(data, UQ(as.symbol(by.field)) == agg.level)%>%as.data.table()
  
  # Re-label any unique identifying information -- only values with identical identifying info will be aggregated at the summation step
  # if sets of counts for two technical replicates are to be summed, the technical replicate label must be made the same across all counts
  # this can be done for example as rep.fields = c('TechRep'), new.rep.levels = c(1)
  for (i in 1:length(rep.fields)){tmp[, c(rep.fields[i]) := rep(new.rep.levels[i], nrow(tmp))]}
  
  # drop CPM, and sum the Counts
  tmp<-ddply(tmp, colnames(tmp)[!grepl('Counts|CPM', colnames(tmp))], summarise, Counts = sum(Counts))%>%as.data.table()
  # Calculate CPM
  tmp<-ddply(tmp, sample.fields, transform, CPM = Counts / sum(Counts) * 1e6)%>%as.data.table()
  
  # if alternate name for the aggregated count entries is given, apply it to the newly aggregated counts
  if (!is.na(new.level)){tmp[,c(by.field) := rep(new.level, nrow(tmp))]}
  
  # Remove the entries for the unaggregated counts, keep only the new aggregated values
  if (remove.old){tmp<-rbind(as.data.table(filter(data, UQ(as.symbol(by.field)) != agg.level)), tmp)}
  else{tmp<-rbind(data, tmp)}
  
  return(tmp)
}

##########################################################################################################################################################
# Compute Expression
# Everything below is either an accessory function or a function used to compute expression
# To Do's:
# * Is error propagation done correctly? We assume uncertainties are uncorrelated, which may be violated for uncertainties of the same element across reps
# 

normalization.step<-function(data, # count matrix
                             normalization.field, # field by which numerator and denominator levels are defined
                             numerator.level, # level for values serving as the numerator
                             denominator.level, # level for values serving as the denominator
                             replicate.field, # field by which duplicate denominator records should be sorted to match numerator records
                             propagate.error, # if True, then propagated error will be returned
                             remove.na){
  # Updated 07-28-2020
  # Rather than checking if number of numerator and denominator records are equal or a multiple, checks whether set of replicates are the same.
  # Now, the user must always specify a replicate field!
  # Additionally, it removes any records which are not represented in both the numerator and denominator set for each replicate. This handles 
  # the case where elements may drop out in the RNA or DNA measurement prior to normalization, which occurs if log is taken before normalization.
  
  # Takes a data.table, and, for a specified field, divides the values for the 
  # numerator level by the corresponding values for second level within that field.
  # If the two sets of values are unequal in size, then the size of the numerator set 
  # must be a multiple of the denominator set. Otherwise, an error will be thrown.
  # * To handle the use of the same denominator records for multiple sets of numerator
  #   records, a field must be designated on which to sort the values so all the measurements
  #   will be properly paired. Typically values are repeated by number of biological replicates
  #   so the matrix will be sorted first based on BiolRep then all other informative fields.
  # If, propagate.error == T, error will be propagated for the resultant divided values. 
  # Error for f(x,y) = x/y --> delta_f = sqrt( (delta_x/y)^2 + (delta_y * -x / y^2)^2)
  # Any infinite or NA values will be removed if remove.na == T 
  
  numerator<-filter(data, UQ(as.symbol(normalization.field)) == numerator.level)%>%as.data.table() # data table with just numerator records
  denominator<-filter(data, UQ(as.symbol(normalization.field)) == denominator.level)%>%as.data.table() # data table with just denominator records
  
  numerator.replicates <- unique(numerator[,get(replicate.field)]) # get set of numerator replicates
  denominator.replicates <- unique(denominator[,get(replicate.field)]) # get set of denominator replicates
  
  if (!setequal(numerator.replicates, denominator.replicates)){ # if there isn't exact one-to-one correspondence of numerator and denominator replicates...
    if (length(denominator.replicates) == 1) { # then there should be only one denominator replicate
      denominator.tmp <- c()
      # if there's only one set of DNA counts, simply use those for each RNA
      for (i in 1:as.integer(nrow(numerator) / nrow(denominator))){
        denominator[,c(replicate.field) := rep(unique(numerator[,get(replicate.field)])[i], nrow(denominator))]
        denominator.tmp<-rbind(denominator.tmp, denominator)
      }
      denominator <- denominator.tmp
      
      # set the order for sorting to data to sort by the field you want to duplicate the denominator for (e.g. duplicate a single biological replicate)
      sort.by = colnames(numerator)[!colnames(numerator)%in%c('ERROR', 'VALUE', normalization.field)]
      setorderv(numerator, sort.by)
      setorderv(denominator, sort.by)
    } else {stop('Numerator and denominator records must have either 1) identical replicate labels or 2) there must be a single denomintor replicate.')}
    # else, throw an error and stop
  }
  
  # for each replicate, the records in both the numerator and denominator sets must be the same (i.e. the same elements are present in both)
  tmp.num<-as.data.table(subset(numerator, 
                                select = colnames(numerator)[!grepl(paste(c('VALUE', 'ERROR', normalization.field), sep = '', collapse = '|'), colnames(numerator))]))
  
  tmp.den<-as.data.table(subset(denominator, 
                                select = colnames(denominator)[!grepl(paste(c('VALUE', 'ERROR', normalization.field), sep = '', collapse = '|'), colnames(denominator))]))
  
  common.records <- dplyr::intersect(tmp.num, tmp.den) # identifying information for common records in both numerator and denominator tables
  
  # create a dummy vector of strings which have all of the identifying information for all records in numerator and denominator, matching info in common.records
  numerator.records   <- tidyr::unite(tmp.num,     'RECORD', colnames(common.records), sep = "::")$RECORD; rm(tmp.num)
  denominator.records <- tidyr::unite(tmp.den, 'RECORD', colnames(common.records), sep = "::")$RECORD; rm(tmp.den)
  # do same for common.records
  common.records <- tidyr::unite(common.records, 'RECORD', colnames(common.records), sep = "::")$RECORD
  # filter records
  numerator   <- filter(numerator,   numerator.records%in%common.records); numerator <- as.data.table(numerator)
  denominator <- filter(denominator, denominator.records%in%common.records); denominator <- as.data.table(denominator)
  #re-sort to ensure everything is in the same order in both tables
  sort.by = colnames(numerator)[!colnames(numerator)%in%c('ERROR', 'VALUE', normalization.field)]
  setorderv(numerator, sort.by)
  setorderv(denominator, sort.by)
  
  # create a temporary table with the same record identifiers as the numerator table, to which we will write the normalized values and error
  tmp<-as.data.table(subset(numerator, 
                            select = colnames(numerator)[!grepl(paste(c('VALUE', 'ERROR', normalization.field), sep = '', collapse = '|'), colnames(numerator))]))
  
  # if the values are not in log space, normalize by division
  if (!unique(numerator$LOG_SCALE) & !unique(denominator$LOG_SCALE)){
    tmp[,VALUE := numerator$VALUE / denominator$VALUE]
    if (propagate.error){
      tmp[,ERROR := sqrt(((numerator$ERROR / denominator$VALUE) ** 2) + ((denominator$ERROR * numerator$VALUE / (denominator$VALUE ** 2)) ** 2))]
    }
    else{tmp[,ERROR := NA]}
  }else if (unique(numerator$LOG_SCALE) & unique(denominator$LOG_SCALE)){ # else, if values are log, subtract 
    tmp[,VALUE := numerator$VALUE - denominator$VALUE]
    if (propagate.error){tmp[,ERROR := sqrt((numerator$ERROR ** 2) + (denominator$ERROR ** 2))]}
    else{tmp[,ERROR := NA]}
  }else{stop('Numerator and demoninator values must both be in either Log or Linear scales before normalization.')}
  
  if(remove.na){tmp<-filter(tmp, !is.na(tmp$VALUE), is.finite(tmp$VALUE))%>%as.data.table()}
  
  return(tmp)
}

log.step<-function(data, # count matrix
                   subset.field = NA, # field for which a subset of measurements at some specified level will be log transformed
                   subset.level = NA, # log transformation will be done for only the measurements corresponding to this level within the specified field
                   propagate.error, # if True, propagated error for the measurement will be returned
                   remove.na){
  
  # For the specified 'VALUE' value in the provided data.table, a log2 transform value will be computed
  # If, propagate.error == T, error will be propagated for the resultant log values 
  # Error for f(x) = log2(x) --> sd_f = sd_x / x / ln(2) 
  # Any infinite or NA values will be removed if remove.na == T 
  
  if (!is.na(subset.field) & !is.na(subset.level)){
    data.remainder <- filter(data, UQ(as.symbol(subset.field) != subset.level))
    data <- filter(data, UQ(as.symbol(subset.field) == subset.level))
  }
  
  tmp<-as.data.table(subset(data, select = colnames(data)[!grepl('VALUE|ERROR', colnames(data))]))
  
  if (propagate.error){tmp[,ERROR := data$ERROR / data$VALUE / log(2)]} # must calculate erorr first, using untransformed measures
  else{tmp[,ERROR := NA]}
  
  tmp[,VALUE := log2(data$VALUE)] # 
  
  if(remove.na){tmp<-filter(tmp, !is.na(tmp$VALUE), is.finite(tmp$VALUE))%>%as.data.table()}
  
  tmp<-tmp%>%as.data.table()
  
  tmp[,LOG_SCALE := rep(TRUE, nrow(tmp))]
  
  if (!is.na(subset.field) & !is.na(subset.level)){tmp <- rbind(data.remainder, tmp)}
  
  return(tmp)
}

averaging.step<-function(data, 
                         average.field, 
                         subset.field = NA, # if you want to average a subset of levels, say average all TRAP counts together, specify 'Fraction' for instance
                         subset.level = NA, # give a level within the subset.field to average across.
                         new.rep.level = NA, # alternate level to assign to the newly averaged subet
                         propagate.error, # if True, will return a propagated error
                         remove.na){
  
  # For a provided data.table, this function averages the values contained in the 'VALUE' column acrossed a specified 'average.field'
  #   i.e. specifying average.field = BC will average all of the VALUE values from Barcodes, while preserving other existing structure (biological replicates, experiment, etc.)
  # If, propagate.error == T, error will be propagated for the resultant averaged values 
  # Error for f(X) =  1/N sum(X) (X is a vector of values, of size N), delta_f = sqrt(sum((delta_xi / N)^2)) -- equivalent to averaging the sd for each element xi of X 
  # Any infinite or NA values will be removed if remove.na == T 
  
  if (!is.na(subset.field) & !is.na(subset.level)){ # if you want to only average over a subset, you must specificy which levels to average and in which field 
    data.remainder <- filter(data, UQ(as.symbol(subset.field)) != subset.level)
    data <- filter(data, UQ(as.symbol(subset.field)) == subset.level)
    data.remainder <- as.data.table(data.remainder)
    data <- as.data.table(data)
  }
  
  
  
  if (propagate.error & !all(is.na(data$ERROR))){
    tmp<-as.data.table(ddply(data, colnames(data)[!grepl(paste(c('VALUE', 'ERROR', average.field), collapse = '|'), colnames(data))], 
                             summarise, 
                             ERROR = sqrt(sum((ERROR / length(ERROR)) ** 2)), 
                             VALUE = mean(VALUE)))
  }
  else{
    tmp<-as.data.table(ddply(data, colnames(data)[!grepl(paste(c('VALUE', 'ERROR', average.field), collapse = '|'), colnames(data))], 
                             summarise, 
                             ERROR = sd(VALUE), # initializes error as the sample standard deviation 
                             VALUE = mean(VALUE)))
  }
  
  tmp<-as.data.table(tmp)
  
  if (!is.na(subset.field) & !is.na(subset.level) & !is.na(new.rep.level)){tmp[,c(average.field) := rep(new.rep.level, nrow(tmp))]}
  
  if(remove.na){tmp<-filter(tmp, !is.na(tmp$VALUE), is.finite(tmp$VALUE))%>%as.data.table()}
  
  if (!is.na(subset.field) & !is.na(subset.level)){tmp <- rbind(data.remainder, tmp)}
  
  return(tmp)
}

sum.step<-function(data, 
                   sum.field, 
                   subset.field = NA, # if you want to average a subset of levels, say average all TRAP counts together, specify 'Fraction' for instance
                   subset.level = NA, # give a level within the subset.field to average across.
                   new.rep.level = NA, # alternate level to assign to the newly averaged subet
                   propagate.error, # if True, will return a propagated error
                   remove.na){
  
  # For a provided data.table, this function returns the sum of the values contained in the 'VALUE' column acrossed the specified 'sum.field'
  #   i.e. specifying sum.field = BC will sum all of the VALUE values from Barcodes, while preserving other existing structure (biological replicates, experiment, etc.)
  # If, propagate.error == T, error will be propagated for the resultant averaged values 
  # Error for f(X) =  sum(X) (X is a vector of values, of size N), delta_f = sqrt(sum((delta_xi)^2)) -- equivalent to adding the sd for each element xi of X 
  # Any infinite or NA values will be removed if remove.na == T 
  
  if (!is.na(subset.field) & !is.na(subset.level)){ # if you want to only average over a subset, you must specificy which levels to average and in which field 
    data.remainder <- filter(data, UQ(as.symbol(subset.field)) != subset.level)
    data <- filter(data, UQ(as.symbol(subset.field)) == subset.level)
    data.remainder <- as.data.table(data.remainder)
    data <- as.data.table(data)
  }
  
  if (propagate.error){
    tmp<-as.data.table(ddply(data, colnames(data)[!grepl(paste(c('VALUE', 'ERROR', sum.field), collapse = '|'), colnames(data))], 
                             summarise, 
                             ERROR = sqrt(sum(ERROR ** 2)), 
                             VALUE = sum(VALUE)))
  }
  else{
    tmp<-as.data.table(ddply(data, colnames(data)[!grepl(paste(c('VALUE', 'ERROR', sum.field), collapse = '|'), colnames(data))], 
                             summarise, 
                             ERROR = NA,
                             VALUE = sum(VALUE)))
  }
  
  tmp<-tmp%>%as.data.table()
  
  if (!is.na(subset.field) & !is.na(subset.level) & !is.na(new.rep.level)){tmp[,c(sum.field) := rep(new.rep.level, nrow(tmp))]}
  
  if(remove.na){tmp<-filter(tmp, !is.na(tmp$VALUE), is.finite(tmp$VALUE))%>%as.data.table()}
  
  if (!is.na(subset.field) & !is.na(subset.level)){tmp <- rbind(data.remainder, tmp)}
  
  return(tmp)
}

compute.expression<-function(data, # count matrix
                             steps = list(c('average','BC'), c('normalize','Fraction', 'rna.input', 'dna.maxi', 'BiolRep'), 
                                          c('log'), c('average','BiolRep')), # list of steps and arguments for transforming the provided input data
                             propagate.error = T, # if True, standard deviation will be propagated through all steps
                             remove.na = T, # removes Nan values at each step
                             output.value.label = 'Expression',  # final column label for the returned transformed values
                             output.error.label = 'SD') # the associated error of the returned transformed values
{
  
  # Takes a data.table of Counts/CPM values and performs a series of transformation steps (averaging, normalizing, and log transforming).
  # Starting with the first averaging step, standard deviation can be propagated through each of the remaining transformations 
  # to provide a final estimate for the error.
  
  tmp<-as.data.table(subset(data, select = colnames(data)[!grepl('CPM|Counts',colnames(data))]))
  tmp[,VALUE := data$CPM]
  tmp[,ERROR := rep(NA, nrow(tmp))]
  tmp[,LOG_SCALE := rep(FALSE, nrow(tmp))]
  
  for (s in steps){
    if (s[1] == 'sum'){
      if (length(s) == 5) {tmp<-sum.step(tmp, s[2], s[3], s[4], s[5], propagate.error = propagate.error, remove.na = remove.na)}
      else {tmp<-sum.step(tmp, s[2], propagate.error = propagate.error, remove.na = remove.na)}}
    if (s[1] == 'average'){
      if (length(s) == 5) {tmp<-averaging.step(tmp, s[2], s[3], s[4], s[5], propagate.error = propagate.error, remove.na = remove.na)}
      else {tmp<-averaging.step(tmp, s[2], propagate.error = propagate.error, remove.na = remove.na)}}
    if (s[1] == 'normalize'){tmp<-normalization.step(tmp, s[2], s[3], s[4], s[5], propagate.error = propagate.error, remove.na = remove.na)}
    if (s[1] == 'log'){
      if (length(s) == 3) {tmp<-log.step(tmp, s[2], s[3], propagate.error = propagate.error, remove.na = remove.na)}
      else {tmp<-log.step(tmp, propagate.error = propagate.error, remove.na = remove.na)}}
  }
  
  tmp<-tmp%>%subset(select = colnames(.)[!grepl('LOG_SCALE',colnames(.))])%>%as.data.table()
  tmp<-setnames(tmp, old = 'VALUE', new = output.value.label)
  tmp<-setnames(tmp, old = 'ERROR', new = output.error.label)
  tmp<-tmp%>%as.data.table()
  
  return(tmp)
}

##########################################################################################################################################################
# Apply Additional Label Values Expression
# Takes two data frames, one consisting of MPRA records, the second consisting of value matched labels to be appended to the table
# Values in both table by which the value labels will be applied to the MPRA table are give as a vector, matching.fields
# The fields in the value table are provided as value.fields

apply.values <- function(d, # data.table, table of MPRA records
                         values,# data.table, table with value labels which will be applied to the table 'd', creating a new column with the labels
                         matching.fields, # vector, listing fields in both tables by which values will be assigned from table 'values' to table 'd'
                         value.fields){ # vector, listing fields in table 'values' which will be appended to table 'd' based on the matching criterion
  
  # Takes a data.table of MPRA records or similar measurements, where each record is associated with a hierarchical organization of replicate levels
  # (i.e. Biological sample, technical replicate, sample fraction, element, barcode, etc.). A second table contains a subset of the distinguishing
  # fields in the data.table d. A specified set of matching.fields is used to assign values from the value table to the table d by creating a new
  # column for each specified value.fields and writing that value to the susbet of d matching the values in each row of the value table.
  # Returns a data table with the additional value columns.
  
  d <- as.data.table(d)
  values <- as.data.table(values)
  
  if (length(matching.fields) >  1){
    d.names <- tidyr::unite(d, NAMES, all_of(matching.fields), sep = '::')$NAMES # join the specified matching fields into a single list for indexing
    values.names <- tidyr::unite(values, NAMES,  all_of(matching.fields), sep = '::')$NAMES # join the specified matching fields into a list for indexing
  } else {
    d.names <- as.character(d[,get(matching.fields[1])])
    values.names <- as.character(values[,get(matching.fields[1])])
  }
  d.new <- c() # initialize collector variable
  for (i in 1:nrow(values)){  # iterate through each row in the value table
    tmp <- filter(d, d.names == values.names[i]); tmp <- as.data.table(tmp) # filter down to the records in d which matched the fiels in values
    for (v in value.fields){ tmp[, c(v) := values[,get(v)][i]] } # assign the values to subset of d, creating new column
    d.new <- rbind(d.new, tmp) # add subset to collector variable
  }
  
  d.new <- as.data.table(d.new)
  return(d.new)
}

##########################################################################################################################################################
# edgeR Count Normalization
# Takes a data.table containing counts data and applies edgeR's normalization factor calculation procedure to generate normalized CPM values
# which reflect differences in library composition in addition to sequencing depth.
# Takes a data.table with one barcode measurement per experiment, per replicate as each row
# Takes a specified set of identifying sample.fields and groups measurements accordingly. 
# Applies calcNormFactors, gets new CPM values
# returns data.table with new CPM values
# Note, any fields which are not specified in the sample.fields, element.fields, or count.field will not be returned in the final output.... 

normalize.counts.edgeR <- function(d, # data.table, MPRA counts table, with each row as a single barcode measure within a distinct experiment, replicate
                                   sample.fields = c('Experiment','Fraction','BiolRep','TechRep'), # vector, set of fields which distinguish replicates
                                   element.fields = c('Element','Allele','BC'), # vector, set of fields which distinguish each element
                                   count.field = 'Counts', # field containing count values
                                   sep_pattern = '::'){ # used to keep sample and element information separate through intermediate steps. 
  # should not be present in sample or element names
  
  d.counts <- tidyr::unite(d, 'Sample', all_of(sample.fields), sep = sep_pattern)           # join columns for sample information into single column
  d.counts <- tidyr::unite(d.counts, 'Element', all_of(element.fields), sep = sep_pattern)  # join columns for element information into single column
  
  d.counts <- subset(d.counts, select = c('Element','Sample','Counts')) # subset to just sample, element, and count information
  
  d.counts <- spread(d.counts, Sample, Counts) # spread into an N element x M Sample matrix
  
  d.counts.matrix <- as.matrix(d.counts[,c(-1)]) # convert to matrix
  
  row.names(d.counts.matrix) <- d.counts[,1] # set row names to element names
  
  y <- DGEList(counts = d.counts.matrix) # convert to a edgeR DGEList object
  
  y <- calcNormFactors(y) # calculate the normalization factors
  
  d.cpm <- cpm(y) # get the cpm values
  
  element.names <- row.names(d.cpm) # get element.names back out
  d.cpm <- as.data.table(d.cpm) # convert to a data.table
  
  d.cpm[, Element := element.names] # set element names as column
  
  d.cpm <- gather(d.cpm, Sample, CPM, colnames(d.cpm)[colnames(d.cpm) != 'Element']) # gather records back to one bc measurement per row
  
  d.cpm <- separate(d.cpm, 'Element', element.fields, sep = sep_pattern) # extract element information back out
  d.cpm <- separate(d.cpm, 'Sample', sample.fields, sep = sep_pattern) # extract sample information back out
  
  return(d.cpm)
}
##########################################################################################################################################################
# Realistic MPRA Simulation
# To Do's:
# * Think of better way to handle elements with nan values -- how to handle replicates with different numbers of measuable (non-zero count) elements.
# * For repetitive sampling/simulating to low p-values, integrate a way to repeat the shuffling procedures N times to return more simulation data 
#   without having to re-do all of the expression calculations and standardization over and over again.

filter.multiple <- function(d, filter.dt, exclude.match = T, exclude.fields = c()){
  filter.match <- c()
  for (i in 1:nrow(filter.dt)){
    check_match <- rep(T, nrow(d))
    for (f in colnames(filter.dt)){if(f%in%colnames(d) & !f%in%exclude.fields){check_match <- (d[,get(f)] == filter.dt[,get(f)][i]) & check_match}}
    filter.match <- cbind(filter.match, check_match)
  }
  filter.match <- rowSums(filter.match) > 0L
  if(exclude.match){filter.match <- !filter.match}
  return(d[filter.match,]%>%as.data.table())
}


# NOTE: The simulation method restricts the values returned to those derived from elements which have a full set of non-nan values across all replicates
# This should still be fine for calibrating Type I error rates, but the comparison with real data may be less relevant if the real data contains many 
# elements with Nan values. I would be especially wary of uncertainty and false positives from elements which are low abundance in counts 
# and only have 2-3 measurements.

simulate.mpra <- function(
  d, # counts data table
  barcode.field = 'BC', # field denoting which barcode each measurement is assigned to for a particular element
  comparison.field = 'Allele', # which field should be used for distinguishing the comparison groups for the count simulation
  comparison.levels = c('ref','alt'), # which levels in the comparison.field will be used. The first level provided will serve as the reference
  replicate.field = 'BiolRep', # which field denotes the replicate label
  fraction.field = 'Fraction', # field denoting the library type ('rna.total', 'dna.maxi', etc)
  rna.level = 'rna.input', # which level in the fraction.field denotes the RNA counts
  dna.level = 'dna.aggmaxi', # which level in the fraction.field denotes the DNA counts
  generate.counts = F){ # if TRUE, modelled CPM values will be converted to Counts by sampling barcodes using their CPM as their sampling probability
  # to a number of sampling events equal to the size of each library
  
  # if we're asking the function to output integer count values, we need to store the actual number of counts from each library in the original dataset.
  if (generate.counts){total.counts <- ddply(filter(d, UQ(as.symbol(fraction.field)) == rna.level), replicate.field, summarise, TOTAL_COUNTS = sum(Counts))%>%as.data.table()}
  
  # filter down to just the dna counts, remove the actual counts column, just preserve the CPM
  dna.counts <- filter(d, UQ(as.symbol(fraction.field)) == dna.level)%>%subset(select = colnames(d)[colnames(d) != 'Counts'])%>%as.data.table()
  
  # check how many replicates of DNA and RNA libraries were suplied
  dna.reps <- as.data.table(filter(d, UQ(as.symbol(fraction.field)) == dna.level))[,get(replicate.field)]%>%unique()
  rna.reps <- as.data.table(filter(d, UQ(as.symbol(fraction.field)) == rna.level))[,get(replicate.field)]%>%unique()
  
  # Check to see if the RNA libraries match the provided DNA libraries. If there isn't perfect 1to1 correspondence, was a common DNA library provided?
  if (all(dna.reps != rna.reps) & length(dna.reps) == 1){
    tmp <- c()
    for (i in rna.reps){ # go through each RNA replicate
      dna.counts[,(replicate.field) := rep(i, nrow(dna.counts))] # copy the single DNA counts and re-label them for each RNA library
      tmp <- rbind(tmp, dna.counts) # add the re-labeled counts to a collection variable
    }
    dna.counts <- tmp; rm(tmp) # Re-assign collection variable as the DNA counts
  }else if (all(dna.reps != rna.reps) | (length(dna.reps) != length(rna.reps))){
    stop('The DNA replicate levels do not match the RNA replicate levels. Please provide paired DNA and RNA replicates or a single DNA replicate for all RNA replicates.')
  }; rm(dna.reps, rna.reps)
  
  # Steps for aggregating CPM values into per-replicate estimates of each element's log2 activity
  agg.steps = list(c('average', barcode.field), c('normalize', fraction.field, rna.level, dna.level, replicate.field), c('log')) # steps for calculating aggregated log-ratio activity estimates
  # run the the calculation to get log2 activity estimates
  d.expression <- compute.expression(d, steps = agg.steps, output.value.label = 'VALUE', output.error.label = 'ERROR') # calculate the log-ratio activity estimate per biological replicate
  # what the informative columns in our activity data table that distinguish element/allele/replicate/etc
  summarise.fields <- colnames(d.expression)[!colnames(d.expression)%in%c(replicate.field, 'VALUE', 'ERROR')] 
  # record just the values for the designated reference allele as the true activity of each element
  group.1.activity <- ddply(filter(d.expression, UQ(as.symbol(comparison.field)) == comparison.levels[1]), summarise.fields, summarise, ACTIVITY = mean(VALUE))%>%as.data.table() 
  # For each expression value across replicates, mean center and record the standard deviation across all of the replicates
  std.residuals <- ddply(d.expression, summarise.fields, summarise, SD = sd(VALUE - mean(VALUE)))%>%as.data.table() # the standard deviation across replicates for each element-allele
  # standardize the activity for each element across all replicates. Scales the distribution for all elements to have mean 0 and unit variance; these values are shuffled amongst elements later
  std.activity <- ddply(d.expression, summarise.fields, transform, ACTIVITY = (VALUE - mean(VALUE)) / sd(VALUE))%>%subset(select = c(summarise.fields, replicate.field, 'ACTIVITY'))%>%as.data.table()
  # re-format the standardized activity to have paired activities for the two comparison groups in each row
  paired.std.activity <- spread(std.activity, comparison.field, ACTIVITY)%>%as.data.table()
  # re-order the table based on the informative columns
  setorderv(paired.std.activity, c(summarise.fields[summarise.fields != comparison.field], replicate.field))
  
  # determine which ref-alt pairs have at least one value in one replicate that had a nan value.
  # Should I just exclude the replicate field from this list? 
  nan.pairs <- subset(paired.std.activity, select = c(summarise.fields[summarise.fields != comparison.field], replicate.field))%>%filter(!complete.cases(paired.std.activity))%>%as.data.table()
  
  # filter dna counts to the just the allele-element values for which there was a non-nan value 
  dna.counts<-filter.multiple(dna.counts, nan.pairs, exclude.fields = c('BiolRep')) # only filter to elements which had a complete set of non-nan values across all replicates.
  tmp <- c() # collector variable
  for (i in unique(dna.counts[,get(replicate.field)])){ # iterate through each replicate and renormalize the counts to CPM after removing elements with nan values
    rep.tmp<-filter(dna.counts, UQ(as.symbol(replicate.field)) == i)%>%as.data.table()
    rep.tmp[, CPM := CPM / sum(CPM) * 1e6]
    tmp<-rbind(tmp, rep.tmp); rm(rep.tmp)
  }
  dna.counts<-tmp; rm(tmp, i)
  
  # filter out nan-value elements from the activitiy and residuals 
  group.1.activity<-filter.multiple(group.1.activity, nan.pairs) # contains one value for each allele for each element
  std.residuals<-filter.multiple(std.residuals, nan.pairs) # contains one value for each allele for each element
  paired.std.activity <-filter.multiple(paired.std.activity, nan.pairs) # contains one value for each allele of each element in each replicate.
  # Filter only by Allele and Element to exclude any elements that don't have a full set of replicates for comparison
  # store the replicate levels for the activity
  replicate.levels <- unique(paired.std.activity[,get(replicate.field)])
  # initialize a storage variable for resampled data
  resampled.std.activity<-c()
  # for each replicate
  for (rl in replicate.levels){
    tmp <- filter(paired.std.activity, UQ(as.symbol(replicate.field)) == rl)%>%as.data.table() # temporarily grab just the standardized activities for that replicate
    resample.index <- sample(1:nrow(tmp), nrow(tmp)) # shuffle by sampling without replacement, based on the number index
    
    tmp[, (comparison.levels[1]) := get(comparison.levels[1])[resample.index]] # reorder the values by index
    tmp[, (comparison.levels[2]) := get(comparison.levels[2])[resample.index]]
    
    resampled.std.activity <- rbind(resampled.std.activity, tmp) # append the resampled values
  }; rm(tmp, rl, replicate.levels)
  # reformat the newly resampled activity values
  synthetic.activity<-gather(resampled.std.activity, UQ(as.symbol(comparison.field)), RESIDUAL_ACTIVITY, all_of(comparison.levels))%>%as.data.table()
  setorderv(synthetic.activity, c(summarise.fields, replicate.field)) # reorder them by the informative columns
  #return(list(synthetic.activity, std.residuals))
  tmp<-c() # expand the per-element standardized residuals (standard deviation across replicates) to fit the synthetic activity
  for (i in unique(synthetic.activity[,get(replicate.field)])){
    std.residuals[, (replicate.field) := rep(i, nrow(std.residuals))]
    tmp <- rbind(tmp, std.residuals)
  }
  std.residuals<-as.data.table(tmp); rm(tmp, i)
  setorderv(std.residuals, c(summarise.fields, replicate.field)) # reorder everything for consistency
  
  #synthetic.activity[, SD_ACTIVITY := std.residuals$SD] # assign the std.residuals to the synthetic activity table
  synthetic.activity <- apply.values(synthetic.activity, std.residuals, c('Element','Allele','BiolRep'),'SD')
  synthetic.activity <- setnames(synthetic.activity, old = "SD", new = "SD_ACTIVITY")
  
  tmp<-group.1.activity # assign the reference activities to a temporary variable
  tmp[,(comparison.field) := rep(comparison.levels[2], nrow(tmp))] # re-label these as the secondary comparison group
  mean.activity<-rbind(group.1.activity, tmp) # combine the repeated activities
  tmp<-c()
  for (i in unique(synthetic.activity[,get(replicate.field)])){ # copy mean activity to fit number of replicates
    mean.activity[, (replicate.field) := rep(i, nrow(mean.activity))]
    tmp <- rbind(tmp, mean.activity)
  }
  mean.activity<-tmp; rm(tmp, i)
  setorderv(mean.activity, c(summarise.fields, replicate.field)) # reorder everything acording to the informative columns
  
  synthetic.activity[, MEAN_ACTIVITY := mean.activity$ACTIVITY] # assign mean activities to synthetic activity 
  
  synthetic.activity[, ACTIVITY := RESIDUAL_ACTIVITY * SD_ACTIVITY + MEAN_ACTIVITY] # construct actual synthetic activity
  
  tmp<-c() # re-assign barcodes to the element-allele-replicate synthetic activity values, expand synth activity to fit no. barcodes
  for (i in unique(dna.counts[,get(barcode.field)])){
    synthetic.activity[,(barcode.field) := rep(i, nrow(synthetic.activity))]
    tmp<-rbind(tmp, synthetic.activity)
  }
  synthetic.activity<-tmp; rm(tmp, i)
  
  setorderv(synthetic.activity, c(summarise.fields, barcode.field, replicate.field)) # order synthetic activity and dna counts with same columns
  setorderv(dna.counts, c(summarise.fields, barcode.field, replicate.field))
  
  synthetic.counts <- c() # collector for synthetic counts
  for (i in unique(synthetic.activity[,get(replicate.field)])){
    tmp.activity <- filter(synthetic.activity, UQ(as.symbol(replicate.field)) == i)%>%as.data.table()
    # only pull the counts for measurements in this particular replicate, for elements which are included in the activity table
    tmp.dna <- filter(dna.counts, UQ(as.symbol(replicate.field)) == i, Element%in%tmp.activity$Element)%>%as.data.table() 
    tmp.activity[, CPM := (2 ** ACTIVITY) * tmp.dna$CPM] # use synthetic activity and barcode count abundance (as CPM) to generate synthetic RNA CPM
    if (generate.counts){ # if generate.counts == T, turn CPM into simulated integer count variables
      library.size <- total.counts[get(replicate.field) == i,]$TOTAL_COUNTS # library size matches the provided size (will be slightly different if nan values were removed)
      sample.counts <- sample(1:nrow(tmp.activity), library.size, replace = T, prob = tmp.activity$CPM / sum(tmp.activity$CPM)) # sum over range of bc indexes, for as many times as counts in library
      tmp.activity[, Counts := as.numeric(table(factor(sample.counts, levels = c(1:nrow(tmp.activity)))))] # turn list of sampled indexes into table with same levels as barcode indexes
      tmp.activity[, CPM := Counts / sum(Counts) * 1e6] # write CPM as well
    } else { tmp.activity[, CPM := CPM / sum(CPM) * 1e6] } # else if not returning counts, re-normalize the simulated values to CPM
    synthetic.counts <- rbind(synthetic.counts, tmp.activity) # add counts to collector
  }; rm(tmp.activity, tmp.dna, i)
  # label these counts with the same rna.label as provided
  synthetic.counts[, (fraction.field) := rep(rna.level, nrow(synthetic.counts))]
  
  # ensure new synthetic rna counts and original dna counts are formated in tables with the same field information
  if (!generate.counts){ # if not generating pseduo counts, subset both the synthetic rna counts and original dna counts to all fields except for count
    synthetic.counts<-synthetic.counts%>%subset(select = colnames(d)[colnames(d) != 'Counts'])%>%as.data.table()
    dna.counts <- filter(d, UQ(as.symbol(fraction.field)) == dna.level)%>%subset(select = colnames(d)[colnames(d) != 'Counts'])%>%as.data.table()}
  else{ # else, do include counts in both tables
    synthetic.counts<-synthetic.counts%>%subset(select = colnames(d))%>%as.data.table()
    dna.counts <- filter(d, UQ(as.symbol(fraction.field)) == dna.level)%>%subset(select = colnames(d))%>%as.data.table()}
  # filter down to just the dna elements which were modeled in the simulation, after removing elements with nan values
  
  # tmp <- c() # re-normalize the dna counts after removing elements with nan-values
  # for (i in unique(dna.counts[,get(replicate.field)])){
  #   # get the elements in the corresponding synthetic counts table for this replicate 'i'
  #   non.nan.elements <- filter(synthetic.counts, UQ(as.symbol(replicate.field)) == i)$Element; print(length(non.nan.elements))
  #   rep.tmp<-filter(dna.counts, UQ(as.symbol(replicate.field)) == i, Element%in%non.nan.elements)%>%as.data.table()
  #   if (generate.counts){rep.tmp[, CPM := Counts / sum(Counts) * 1e6]} # normalize the remaining counts after filtering
  #   else {rep.tmp[, CPM := CPM / sum(CPM) * 1e6]} # CPM must also be normalized if elements are removed
  #   tmp<-rbind(tmp, rep.tmp); rm(rep.tmp)
  # }
  # dna.counts<-tmp; rm(tmp, i)
  # recombine the synthetic and dna counts into a single table
  synthetic.counts<-rbind(synthetic.counts, dna.counts)
  
  return(synthetic.counts)
}

##########################################################################################################################################################
# Outlier Detection

# To Do's
# * Fix the entire function apparently -- it doesn't work for 410 data and it doesn't work for passing column/field names to ddply as strings through get()
#   - Those strings must be specified as global variables for it to work
#   - Essentially, this comes down to figuring out how to use ddply properly

identify.outlier <- function(d, value.field = 'Expression', by.field = 'BC', within.field = 'BiolRep', exclude.fields = c('TechRep','SEM'), return.direction = F){
  d.expression <- copy(d)
  d.expression <- setnames(d.expression, value.field, 'VALUE'); value.field <- 'VALUE'
  std.bc.activity <- ddply(d.expression, colnames(d.expression)[!colnames(d.expression)%in%c(by.field, value.field)], transform, STD_ACTIVITY = (VALUE - mean(VALUE)) / sd(VALUE))
  std.bc.activity <- subset(std.bc.activity, select = colnames(std.bc.activity)[!colnames(std.bc.activity)%in%c(value.field, exclude.fields, within.field)])
  std.bc.average.activity <- ddply(std.bc.activity, colnames(std.bc.activity)[!colnames(std.bc.activity)%in%c('STD_ACTIVITY')], summarise, AVG_STD_ACTIVITY = mean(STD_ACTIVITY))
  std.bc.average.activity <- std.bc.average.activity%>%filter(complete.cases(.))
  std.bc.average.activity <- ddply(std.bc.average.activity, colnames(std.bc.activity)[!colnames(std.bc.activity)%in%c(by.field, 'STD_ACTIVITY')], 
                                   transform, Q1 = as.numeric(quantile(AVG_STD_ACTIVITY, 0.25)), Q3 = as.numeric(quantile(AVG_STD_ACTIVITY, 0.75)))
  std.bc.average.activity <- as.data.table(std.bc.average.activity)
  std.bc.average.activity[, OUTLIER := (AVG_STD_ACTIVITY < (Q1 - 1.5 * (Q3 - Q1))) | (AVG_STD_ACTIVITY > (Q3 + 1.5 * (Q3 - Q1)))]
  
  if (return.direction){
    std.bc.average.activity[, HI := (AVG_STD_ACTIVITY > (Q3 + 1.5 * (Q3 - Q1)))]
    std.bc.average.activity[, LO := (AVG_STD_ACTIVITY < (Q1 - 1.5 * (Q3 - Q1)))]
  }
  
  return.fields <- colnames(d.expression)[!colnames(d.expression)%in%c(exclude.fields, within.field, value.field)]
  if (return.direction){return.fields <- c(return.fields, 'HI', 'LO')}
  outlier.bcs <- filter(std.bc.average.activity, OUTLIER)%>%subset(select = return.fields)%>%as.data.table()
  
  return(outlier.bcs)
}

# Filter Outliers
remove.outliers <- function(counts, outliers){
  # "data poisioning approach"
  # strip original data and outliers to common features, label one as outliers, the other as not
  # add the two tables together
  # use a tranform() function to change any similar entries that match any of the outliers to T 
  #   Transform will group similar entries by Element, Allele, Barcode number etc. and if any 
  #   in the group were True (ie from the table of outliers appended), then all will be labeld true.
  
  common.fields <- colnames(outliers)[colnames(outliers)%in%colnames(counts)]
  
  tmp<-subset(counts, select = common.fields)
  
  tmp.outliers <-subset(outliers, select = common.fields)
  
  tmp[, OUTLIER := rep(F, nrow(tmp))]
  tmp[, id := 1:nrow(tmp)]
  tmp.outliers[, OUTLIER := rep(T, nrow(tmp.outliers))]
  tmp.outliers[, id := rep(NA, nrow(tmp.outliers))]
  
  tmp <- rbind(tmp, tmp.outliers)
  
  tmp <- ddply(tmp, common.fields, transform, OUTLIER = any(OUTLIER))
  
  return(counts[filter(tmp, !OUTLIER)$id,])
}

##########################################################################################################################################################
# Well Measured
# This function filters elements based on representation in a specified sequencing library pool

# To-Do's:
# * 

well.measured<-function(counts, fraction.field, fraction.level, return.levels, min.counts = 20, min.bc = 3){
  
  d <- filter(counts, UQ(as.symbol(fraction.field)) == fraction.level)
  
  d <- ddply(d, c('Element', 'Allele', 'BC'), Counts = mean(Counts))
  
  d <- filter(d, Counts >= min.counts)%>%as.data.table()
  
  d <- ddply(d, return.levels, summarise, N = length(BC))
  
  d <- filter(d, N >= min.bc)
  
  d <- subset(d, select = return.levels)
  
  d <- as.data.table(d)
  
  return(d)
}

##########################################################################################################################################################
# CPM Threshold Determination

# Filter set of elements based on available count data to determine which to include in hypothesis testing
#   Takes a matrix of counts used to calculate log activity ratios for a set of elements and parameters to be use in downstream tests, returns which
#     elements are measured sufficiently, thereby minimizing the total number of tests performed.
# Specify parameters of how the elements will be tested (i.e by what Fractions measurements will be normalized to calculate the log activity score)

# Arguments
#   normalization.field -- field in counts table used to specify which records to normalize by to obtain the log activity ration 
#     (specifies a numerator and denominator field)
#   num.level -- numerator level in the normalization.field column
#   den.level -- denominator level in normalization.field column
#   comparison.field -- field by which elements are compared for the hypothesis testing
#   comparison.levels -- levels in the comparison.field column which are actually tested by hypothesis testing
#   element.fields -- fields which uniquely identify each element, independent of the allele
#   replicate.field -- field identifying which replicate each record came from
#   min.replicate, min.bc, min.count -- integer values specifying the minimum values each element can have 
#     for inclusion in the final set of testable elements. (i.e. an element must have the minimum number of counts in at least 
#     the minimum number of barcodes in at least the minimum number of replicates, in that order)
#   paired -- specifies if both alleles/elements for a comparison must reach barcode and count thresholds for that replicate to be counted for the pair,
#     which is necessary to specify if a paired hypothesis testing model is used (i.e. for paired t-test). Paired is true by default

filter_elements <- function(counts, 
                            normalization.field = 'Fraction', 
                            num.level = 'rna.input', den.level = 'dna.maxi', 
                            comparison.field = 'Allele', comparison.levels = c('ref','alt'), 
                            element.fields = c('Element'),
                            replicate.field = 'BiolRep',
                            min.replicate = 4, min.bc = 3, min.count = 10, paired = T){
  
  d <- filter(counts, UQ(as.symbol(normalization.field)) == num.level)
  d <- filter(d, Counts >= min.count)
  d <- ddply(d, c(element.fields, comparison.field, replicate.field), summarise, N_BC = length(Counts))
  d <- filter(d, N_BC >= min.bc)
  if (paired){
    d <- setnames(d, old = comparison.field, new = 'COMPARISON.FIELD')
    d <- ddply(d, c(element.fields, replicate.field), here(summarise), Paired = all(comparison.levels%in%COMPARISON.FIELD))
    d <- filter(d, Paired)
  }
  d <- setnames(d, old = replicate.field, new = 'REPLICATE.FIELD')
  d <- ddply(d, c(element.fields), here(summarise), N_BiolRep = length(REPLICATE.FIELD))
  d <- filter(d, N_BiolRep >= min.replicate)
  
  num.elements <- subset(d, select = element.fields)
  
  d <- filter(counts, UQ(as.symbol(normalization.field)) == num.level)
  d <- filter(d, Counts >= min.count)
  d <- ddply(d, c(element.fields, comparison.field, replicate.field), summarise, N_BC = length(Counts))
  d <- filter(d, N_BC >= min.bc)
  if (paired){
    d <- setnames(d, old = comparison.field, new = 'COMPARISON.FIELD')
    d <- ddply(d, c(element.fields, replicate.field), here(summarise), Paired = all(comparison.levels%in%COMPARISON.FIELD))
    d <- filter(d, Paired)
  }
  d  <- setnames(d, old = replicate.field, new = 'REPLICATE.FIELD')
  if (length(unique(d$REPLICATE.FIELD)) > 1){
    d <- ddply(d, c(element.fields), here(summarise), N_BiolRep = length(REPLICATE.FIELD))
    d <- filter(d, N_BiolRep >= min.replicate)
  }
  
  den.elements <- subset(d, select = element.fields)
  
  return(dplyr::intersect(num.elements, den.elements))
}

# Determine Count Threshold for Minimization of Multiple Tests
# default method uses t.test for hypothesis test
# prioritizes finding 1) highest number of hits with 2) highest number of barcodes per element and 3) highest number of replicates. 
# then optimizes count levels
get_count_threshold <- function(d, replicate.field = 'BiolRep', 
                                fraction.field = 'Fraction', rna.level = 'rna.input', dna.level = 'dna.maxi', 
                                aggregation.steps = list(c('average', 'BC'), 
                                                         c('normalize', 'Fraction', 'rna.input', 'dna.maxi', 'BiolRep'), 
                                                         c('log')),
                                sample.fields = c('BiolRep','Fraction'), 
                                element.fields = c('Element','Allele','BC'),
                                rep.levels = c(4),
                                bc.levels = c(3),
                                count.resolution = 20,
                                hit.threshold = 0.05,
                                count.minimum = 10,
                                count.max = NA,
                                show.plot = T){
  
  # calculate p-values for entire set
  d.norm <- normalize.counts.edgeR(d, sample.fields = sample.fields, element.fields = c('Element','Allele','BC'), count.field = 'Counts')
  d.expr <- compute.expression(d.norm, 
                               steps = aggregation.steps)
  res <- t.test.by.group(d.expr, show.warnings = F)
  
  if (is.na(count.max)) {count.max  <- quantile(d$Counts, 0.75) + 1}
  count.levels <- seq(count.minimum, count.max, length.out = count.resolution) # set the discrete values at which to test count thresholds
  
  if (length(rep.levels) != 1 | length(bc.levels) != 1){ # if there are other parameters besides count threshold to optimize
    
    # create table for recording results
    results.by.threshold <- data.table(CountThreshold = NA, RepThreshold = NA, BCThreshold = NA, hits = NA)
    
    # populate table by iterating through all combinations of thresholds at count, replicate, and barcode level
    tmp <- c()
    for (i in count.levels){
      sub.tmp <-copy(results.by.threshold)
      sub.tmp[, CountThreshold := i]
      tmp <- rbind(tmp, sub.tmp)
    }
    results.by.threshold <- tmp
    
    # populate replication thresholds
    tmp <- c()
    for (i in rep.levels){
      sub.tmp <-copy(results.by.threshold)
      sub.tmp[, RepThreshold := rep(i, nrow(sub.tmp))]
      tmp <- rbind(tmp, sub.tmp)
    }
    results.by.threshold <- tmp
    
    # populate barcode thresholds
    tmp <- c()
    for (i in bc.levels){
      sub.tmp <-copy(results.by.threshold)
      sub.tmp[, BCThreshold := rep(i, nrow(sub.tmp))]
      tmp <- rbind(tmp, sub.tmp)
    }
    results.by.threshold <- tmp; rm(tmp, sub.tmp, i)
    
    # iterate through each set of thresholds, get set of elements which meet criterion, re-adjust p-values and see which met criterion
    pb <- pb <- progress_bar$new(format = " performing grid search [:bar] :percent eta: :eta",total = nrow(results.by.threshold), clear = FALSE, width= 60)
    for (i in 1:nrow(results.by.threshold)){
      
      elements.to.test <- filter_elements(d, 
                                          normalization.field = fraction.field, 
                                          num.level = rna.level, 
                                          den.level = dna.level, 
                                          min.replicate = results.by.threshold$RepThreshold[i], 
                                          min.bc = results.by.threshold$BCThreshold[i], 
                                          min.count = results.by.threshold$CountThreshold[i])
      
      tmp <- filter(res, Element%in%elements.to.test$Element)
      tmp <- as.data.table(tmp)
      tmp[, fdr := p.adjust(pval, method = 'BH')]
      
      results.by.threshold$hits[i] <- sum(tmp$fdr < hit.threshold)
      pb$tick()
    }
    
    if (show.plot) { 
      bc.labels <- paste0("Min BC = ", results.by.threshold$BCThreshold)
      names(bc.labels) <- results.by.threshold$BCThreshold
      
      p <- ggplot(results.by.threshold, aes(x = CountThreshold, y = RepThreshold, fill = hits)) +
        geom_tile(color = 'black', size=0.5) + 
        facet_grid(.~BCThreshold, labeller = labeller(BCThreshold = bc.labels)) + 
        theme(text = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        scale_y_continuous(breaks = rep.levels) +
        scale_x_continuous(breaks = round(count.levels)) + ggtitle('Parameter Grid Search') + xlab('Count Threshold') + ylab('Replicate Threshold')
      
      plot(p)}
    
    max.hits <- max(results.by.threshold$hits)
    
    thresholds.by.max.hits <- filter(results.by.threshold, hits == max.hits)
    
    bc.threshold <- max(thresholds.by.max.hits$BCThreshold) # prioritize maximizing the number of barcodes over the number of replicates
    
    thresholds.by.max.hits <- filter(thresholds.by.max.hits, BCThreshold == bc.threshold)
    
    rep.threshold <- max(thresholds.by.max.hits$RepThreshold)
    
    thresholds.by.max.hits <- filter(thresholds.by.max.hits, RepThreshold == rep.threshold)
    
  } 
  else { bc.threshold <- bc.levels[1]; rep.threshold <- rep.levels[1] } # default to provide values if only one
  
  results.by.threshold <- data.table(CountThreshold = NA, RepThreshold = rep.threshold, BCThreshold = bc.threshold, hits = NA)
  
  tmp <- c()
  for (i in count.levels){
    sub.tmp <-copy(results.by.threshold)
    sub.tmp[, CountThreshold := i]
    tmp <- rbind(tmp, sub.tmp)
  }
  results.by.threshold <- tmp; rm(tmp, sub.tmp, i)
  
  pb <- pb <- progress_bar$new(format = " evaluating count threshold [:bar] :percent eta: :eta",total = nrow(results.by.threshold), clear = FALSE, width= 60)
  for (i in 1:nrow(results.by.threshold)){
    elements.to.test <- filter_elements(d, normalization.field = fraction.field, num.level = rna.level, den.level = dna.level, min.replicate = results.by.threshold$RepThreshold[i], min.bc = results.by.threshold$BCThreshold[i], min.count = results.by.threshold$CountThreshold[i])
    tmp <- filter(res, Element%in%elements.to.test$Element); tmp <- as.data.table(tmp)
    tmp[, fdr := p.adjust(pval, method = 'BH')]
    results.by.threshold$hits[i] <- sum(tmp$fdr < hit.threshold)
    pb$tick()
  }
  
  fit <- loess(hits ~ CountThreshold, results.by.threshold)
  if (show.plot){
    p <- ggplot(results.by.threshold, aes(x = CountThreshold, y = hits)) + 
      stat_smooth(method = 'loess') + 
      geom_point() + 
      coord_equal(ratio = 1.5 * max(results.by.threshold$CountThreshold) / max(results.by.threshold$hits)) + 
      geom_vline(xintercept = round(fit$x[(fit$fitted == max(fit$fitted))] - fit$s), color = 'red') + 
      xlab('Count Threshold') + 
      ylab(paste('Number Significnat at FDR < ', hit.threshold, sep = '')) + theme(text = element_text(size=12))
    plot(p)
  }
  
  count.threshold <- max(c(round(fit$x[(fit$fitted == max(fit$fitted))] - fit$s), 1))
  
  final_thresholds <- c(count.threshold, bc.threshold, rep.threshold)
  names(final_thresholds) <- c('Count Threshold', 'BC Threshold', 'Replicate Threshold')
  
  return(final_thresholds)
}

##########################################################################################################################################################
# Assign Library Annotations and Test for Enrichment/Burden

# Test for Enrichment
# Test for the independence of two categorical variables for a set of variants
# Specifically, a set of variants annotated with two boolean variables is submitted, denoting binary class membership for each variant in two classes
# The enrichment of members of the v2 class within the v1 class is evaluated using one of three tests (hypergeometric, Fisher exact, or binomal)
# returns a p-value, and if possible, and odds ratio
test.for.enrichment <- function(d, v1, v2, test.type = 'hyper.test'){
  
  n <- nrow(filter(d, UQ(as.symbol(v1)),  UQ(as.symbol(v2))))
  k <- nrow(filter(d, !UQ(as.symbol(v1)),  UQ(as.symbol(v2))))
  N <- nrow(filter(d, UQ(as.symbol(v1))))
  K <- nrow(filter(d, !UQ(as.symbol(v1))))
  
  if (test.type == 'binom.test'){p <-  binom.test(n, n + k, N / (N + K))$p.value}
  else if (test.type == 'hyper.test'){p <- phyper(n - 1, N, K, n + k, lower.tail = F)}
  else if (toupper(test.type) == 'FET'){p <- fisher.test(matrix(c(n, N - n, k, K - k), 2, 2), alternative = 'greater')$p.value}
  else {stop('test.type must be binom.test, hyper.test, or FET')}
  
  if (k > 0 & N > 0 & K > 0){OR <- (n / k) / (N / K)} else { OR <- NA }
  
  res <- c(p, OR)
  names(res) <- c('P-Value', 'Odds Ratio')
  return(res)
}

# the pipeline should be to 1) assign annotations and 2) send separate lists of annotated into the enrichment function
# test.for.enrichment <- function(significant, background, test.type = 'hyper.test'){
#   
#   total.case <- sum(background$Pheno == 'case')
#   total.control <- sum(background$Pheno == 'control')
#   
#   significant.case <- sum(significant$Pheno == 'case')
#   significant.control <- sum(significant$Pheno == 'control')
#   
#   if (test.type == 'binom.test'){
#     p <- binom.test(significant.case, nrow(significant), p = total.case / nrow(background), alternative = 'greater')$p.value # test that proportion of cases in significant results is greater than
#   } 
#   else if (test.type == 'hyper.test') {
#     p <- phyper(significant.case - 1, total.case, total.control, nrow(significant), lower.tail = F)
#   } else if (toupper(test.type) == 'FET') {
#     p <- fisher.test(matrix(c(significant.case, 
#                               total.case - significant.case, 
#                               significant.control, 
#                               total.control - significant.control), 2, 2), alternative = 'greater')$p.value
#   } else {
#     p <- phyper(significant.case - 1, total.case, total.control, nrow(significant), lower.tail = F)
#   }
#   
#   
#   if (significant.control > 0 & total.case > 0 & total.control > 0){
#     OR <- (significant.case / significant.control) / (total.case / total.control) # odds ratio of being proband with functional 3'UTR mutation
#   } else { OR <- NA }
#   
#   res <- c(p, OR)
#   names(res) <- c('P-Value', 'Odds Ratio')
#   
#   return(res)
# }

# slightly deprecated, not relevant for 2.0
assign.annotations <- function(results, annotation){
  
  tmp <- copy (results)
  
  tmp[,Chr := NA]
  tmp[,hg19Start := NA]
  tmp[,Strand := NA]
  tmp[,Ref := NA]
  tmp[,Alt := NA]
  tmp[,Pheno := NA]
  tmp[,pLI := NA]
  tmp[,LoF_Intolerant := NA]
  
  
  for (i in 1:nrow(tmp)){
    element.annotation <- filter(annotation, seqName == tmp$Element[i])
    tmp$Chr[i] <- element.annotation$Chr[1]
    tmp$hg19Start[i] <- element.annotation$hg19Start[1]
    tmp$Strand[i] <- element.annotation$Strand[1]
    tmp$Ref[i] <- element.annotation$Ref[1]
    tmp$Alt[i] <- element.annotation$Alt[1]
    tmp$Pheno[i] <- element.annotation$Pheno[1]
    tmp$pLI[i] <- element.annotation$pLI[1]
    tmp$LoF_Intolerant[i] <- element.annotation$LoF_Intolerant[1]
  }
  
  return(tmp)
}

##########################################################################################################################################################
# Linear Mixed Model for Hypothesis testing

# To Do's
# * provide option to use more exhaustive ways for estimating p-value than LRT and Satterthwaite
# * devise a way to return more than just the single FC and pvalue for one coefficient of the model
# Tony's notes - this is crazy.  anova(model1,model2) gives the same result.  Consider using this instead.

# Likelihood Ratio Test using the difference in deviance between two linear models and the difference in degrees of freedom to 
# obtain p-value from a Chi-squared distribution.
likelihood_ratio_test<-function(model1, model2){
  # Model1 -- complex model (fewer degrees of freedom, aka more terms)
  # Model2 -- reduced model (more degrees of freedom, aka fewer terms)
  # When comparing a LMM as the complex model, run the model with ML, i.e. set lmer(, REML = FALSE)
  # When comparing random effects structures between LMM, use REML, i.e. use lmer(, REML = TRUE) (which is the default)
  # Performs a likelihood ratio test by comparing the deviance for two models
  # Difference in deviance and degress of freedom for each model 
  
  dev1 <- -2*logLik(model1)
  dev0 <- -2*logLik(model2)
  
  devdiff <- as.numeric(dev0-dev1)
  dfdiff <- attr(dev1,"df")-attr(dev0,"df")
  if (dfdiff < 1){stop("The first argument must be a model with at least 1 additional DF as the second argument's model.")}
  #cat('Chi-square =', devdiff, '(df=', dfdiff,'), p =', pchisq(devdiff,dfdiff,lower.tail=FALSE))
  return(pchisq(devdiff,dfdiff,lower.tail=FALSE)) # returns p-value for the LRT
}

# Performs statistical inference by a linear mixed effects model for given group comparison amongst library elements based on a log-ratio activity value
# A mixed effects model formula is specified by the 'model' argument
# The library units to be compared across groups is specified by the element.field
# The field across which library members are compared is specified by the comparison.field argument
# Levels for comparison within the comparison.field are specified as ref.level and alt.level 
# The model estimate and p-value from the resulting lmm summary table for the coefficient of interest is specified by result.index
# pvalue.method indicates the method for estimating degrees of freedom or performing an LRT to obtain a p-value for the model parameter of interest
# If pvalue.method is set to 'LRT', the likelihood.ratio.test method will be used to compute the p-value by comparing the specified model to a 
# provided null model. The null model will be identical to the full model with the exclusion of the single fixed effect for you inference is performed.
#
# Returns a results table with a columns:
# * FC -- gives the log2 fold change for the effect size of the comparison. (we assume values provided are already log transformed)
# * pval -- nominal pvalue for the comparison
# * fdr -- Benjamini-Hochberg adjusted pvalues, obtained using p.adjust function.
# * bonferroni -- bonferroni adjusted pvalues, obtained using p.adjust function.

mpratools.lmm<-function(d, # data.table with log-ratio activitiy values for fitting to a specified linear mixed model
                        model = Expression ~ Allele + (1|BiolRep), # linear mixed model formula
                        reduced.model = Expression ~ (1|BiolRep), # if pvalue.method == 'LRT', must provide an appropriate null model for comparsion
                        element.field = 'Element', # field designating within what variable the group comparison should be made
                        comparison.field = 'Allele',  # field designating the grouping variable
                        ref.level = 'ref',  # level within the comparison field denoting the reference group for the comparison, level set to '0'
                        alt.level = 'alt', # level within the comparison field denoting the alternative group for the comparison, level set to '1'
                        result.index = 2, # which row in the summary table should the effect size and pvalue be taken
                        pvalue.method = 'LRT',  # method for estimating the pvalue, 
                        # defaults to 'Satterthwaite' method for computing denominator degrees of freedom for anova table
                        run.diagnostic = F,
                        diagnostic.model = Expression ~ Allele,
                        remove.na = T){ # removes nan results from final table, prior to p-value correction.
  
  # filter provided data to just the two factors for comparison, typically by Allele
  d<-filter(d, UQ(as.symbol(comparison.field)) == ref.level | UQ(as.symbol(comparison.field)) == alt.level)%>%as.data.table()
  d[, (comparison.field) := factor(get(comparison.field), levels = c(ref.level, alt.level)) ] # re-order the factor
  
  results<-data.table(logFC = rep(NA, length(unique(d[,get(element.field)]))), # initialize a table for the results
                      pval = rep(NA, length(unique(d[,get(element.field)]))))
  
  results[,(element.field) := d[,unique(d[,get(element.field)])]] # label each row for each element comparison
  
  if(run.diagnostic){
    results[,singular.fit := rep(NA, nrow(results))]
    results[,diagnostic.pval := rep(NA, nrow(results))]}
  
  for (i in 1:nrow(results)){
    
    # filter down to data for a particular element
    d.element<-filter(d, UQ(as.symbol(element.field)) == results[,get(element.field)][i])
    
    # make sure there are at least two observations of the ref and alt labels
    if(nrow(filter(d.element, UQ(as.symbol(comparison.field)) == ref.level)) <  2 | nrow(filter(d.element, UQ(as.symbol(comparison.field)) == alt.level)) < 2){
      warning(paste('Not enough levels to compare for element:', results[,get(element.field)][i], sep = ' ')); next}
    
    tryCatch({
      res<-lmer(model, d.element)
      
      full<-lmer(model, REML = FALSE, d.element)
      if (pvalue.method == 'LRT'){null<-lmer(reduced.model, REML = FALSE, d.element)}
      
      results$logFC[i]<- summary(res)$coefficients[result.index, 1] # pull the relative effect size 
      
      if (pvalue.method == 'Satterthwaite') {results$pval[i]<- summary(res)$coefficients[result.index, 5]}
      else if (pvalue.method == 'LRT'){results$pval[i] <- likelihood_ratio_test(full, null)}
      else {results$pval[i]<- summary(res)$coefficients[result.index, 5]}
      
      if (run.diagnostic){
        results$singular.fit[i] <- isSingular(res)
        full <- lmer(model, REML = FALSE, d.element)
        null <- lm(diagnostic.model, d.element)
        results$diagnostic.pval[i] <- likelihood_ratio_test(full, null)
      }
      
    },
    # case that only one replicate for a given term is measured for an Element
    error = function(err){print(paste("Error:  ", err))})
    
  }
  
  if(remove.na){results<-results[complete.cases(results),]}
  
  results[,fdr:=p.adjust(pval, method = 'BH')]
  results[,bonferroni:=p.adjust(pval, method = 'bonferroni')]
  
  if (run.diagnostic){return.fields <- c(element.field, 'logFC', 'pval', 'fdr', 'bonferroni', 'singular.fit', 'diagnostic.pval')}
  else{return.fields <- c(element.field, 'logFC', 'pval', 'fdr', 'bonferroni')}
  
  results<-results%>%subset(select = return.fields)%>%as.data.table()
  
  return(results)
}

#Tony's inclusion, interaction term, right now uses restricted maximum likelihood
#and Satterthwaite F tests for model significance - this is safest, but need
#better explanation/ML testing of reduced models using LRT (ChiSquare) -
#there is a reduced model in there for this, I will incorporate later
#Pulls logFC from the coefficients using summary() and p-values using anova()
mpratools.lmm.int<-function(d, # data.table with log-ratio activitiy values for fitting to a specified linear mixed model
                            model = Expression ~ Allele*SpliceStatus + (1|BC), # linear mixed model formula
                            reduced.model= Expression ~ Allele + SpliceStatus + (1|BC),
                            element.field = 'Element', # field designating within what variable the group comparison should be made
                            comparison.field = 'Allele',  # field designating the grouping variable
                            ref.level = 'ref',  # level within the comparison field denoting the reference group for the comparison, level set to '0'
                            alt.level = 'alt', # level within the comparison field denoting the alternative group for the comparison, level set to '1'
                            result.index.allele = 2, # which row in the summary table should the effect size and pvalue be taken
                            result.index.secondterm = 3, # which row in the summary table should the effect size and pvalue for the second term be taken
                            result.index.interaction = 4, # which row in the summary table should the effect size and pvalue for the interaction be taken
                            method = "Satterthwaite") #alternate hypothesis testing here is "LRT", which uses anova(model,reduced.model)
{
  
  
  
  
  # filter provided data to just the two factors for comparison, typically by Allele
  d<-filter(d, UQ(as.symbol(comparison.field)) == ref.level | UQ(as.symbol(comparison.field)) == alt.level)%>%as.data.table()
  d[, (comparison.field) := factor(get(comparison.field), levels = c(ref.level, alt.level)) ] # re-order the factor
  #For the second interaction, it will be alphabetical unless you make them factors and reorder them ahead of time
  #(may not be a bad idea to do anyways) - this will matter when you grab logFC to make it go the right way, although the logFC is
  #not quite as interpretable here
  results<-data.table(logFC_allele = rep(NA, length(unique(d[,get(element.field)]))), # initialize a table for the results
                      pval_allele = rep(NA, length(unique(d[,get(element.field)]))),
                      logFC_second = rep(NA, length(unique(d[,get(element.field)]))),
                      pval_second = rep(NA, length(unique(d[,get(element.field)]))),
                      logFC_interaction = rep(NA, length(unique(d[,get(element.field)]))),
                      pval_interaction = rep(NA, length(unique(d[,get(element.field)])))
  )
  
  results[,(element.field) := d[,unique(d[,get(element.field)])]] # label each row for each element comparison
  
  for (i in 1:nrow(results)){
    
    # filter down to data for a particular element
    d.element<-filter(d, UQ(as.symbol(element.field)) == results[,get(element.field)][i])
    
    # make sure there are at least two observations of the ref and alt labels
    if(nrow(filter(d.element, UQ(as.symbol(comparison.field)) == ref.level)) <  2 | nrow(filter(d.element, UQ(as.symbol(comparison.field)) == alt.level)) < 2){
      warning(paste('Not enough levels to compare for element:', results[,get(element.field)][i], sep = ' ')); next}
    
    tryCatch({
      res<-lmer(model, d.element)
      results$logFC_allele[i]<- summary(res)$coefficients[result.index.allele, 1] # pull the relative allele effect size 
      results$pval_allele[i]<- anova(res)[1,6] # allele p-val
      results$logFC_secondterm[i]<- summary(res)$coefficients[result.index.secondterm, 1] # second term logFC
      results$pval_secondterm[i]<- anova(res)[2,6] # second term pval
      results$logFC_interaction[i]<- summary(res)$coefficients[result.index.interaction, 1] # interaction logFC
      if(method=="LRT"){
      full <- res
      null <- lmer(reduced.model,d.element)
      results$pval_interaction[i]<- anova(full,null)[2,8] # interaction pval
      } else {results$pval_interaction[i]<- anova(res)[3,6]} # interaction pval
    },
    # case that only one replicate for a given term is measured for an Element
    error = function(err){print(paste("Error:  ", err))})
    
  }
  
  results[,fdr_allele:=p.adjust(pval_allele, method = 'BH')]
  results[,bonferroni_allele:=p.adjust(pval_allele, method = 'bonferroni')]
  results[,fdr_secondterm:=p.adjust(pval_secondterm, method = 'BH')]
  results[,bonferroni_secondterm:=p.adjust(pval_secondterm, method = 'bonferroni')]
  results[,fdr_interaction:=p.adjust(pval_interaction, method = 'BH')]
  results[,bonferroni_interaction:=p.adjust(pval_interaction, method = 'bonferroni')]
  
  
  return.fields <- c(element.field, 'logFC_allele', 'pval_allele', 'fdr_allele', 'bonferroni_allele', 'logFC_secondterm', 'pval_secondterm', 'fdr_secondterm', 'bonferroni_secondterm','logFC_interaction', 'pval_interaction', 'fdr_interaction', 'bonferroni_interaction')
  
  results<-results%>%subset(select = return.fields)%>%as.data.table()
  
  return(results)
}
##########################################################################################################################################################
# Linear Model for Hypothesis testing

# To-Do's
# * Apply a better way to handle cases where less than two observations are retained for one or more terms in the specified model

# Performs statistical inference by a linear model for given group comparison amongst library elements based on a log-ratio activity value
# A linear model formula is specified by the 'model' argument
# The library units to be compared across groups is specified by the element.field
# The field across which library members are compared is specified by the comparison.field argument
# Levels for comparison within the comparison.field are specified as ref.level and alt.level 
# The model estimate and p-value from the resulting lmm summary table for the coefficient of interet is specified by result.index
#
# Returns a results table with a columns:
# * FC -- gives the log2 fold change for the effect size of the comparison. (we assume values provided are already log transformed)
# * pval -- nominal pvalue for the comparison
# * fdr -- Benjamini-Hochberg adjusted pvalues, obtained using p.adjust function.
# * bonferroni -- bonferroni adjusted pvalues, obtained using p.adjust function.

mpratools.lm<-function(d, # data.table with log-ratio activitiy values for fitting to a specified linear model 
                       model = Expression ~ Allele,  # linear model formula
                       element.field = 'Element', # field designating within what variable the group comparison should be made
                       comparison.field = 'Allele',  # field designating the grouping variable
                       ref.level = 'ref',  # level within the comparison field denoting the reference group for the comparison, level set to '0'
                       alt.level = 'alt', # level within the comparison field denoting the alternative group for the comparison, level set to '1'
                       result.index = 2, # which row in the summary table should the effect size and pvalue be taken
                       remove.na = T){
  
  # The hypothesis test results (effect size and pval) will be returned for the first term in the provided model formula
  
  # filter provided data to just the two factors for comparison, typically by Allele
  d<-filter(d, UQ(as.symbol(comparison.field)) == ref.level | UQ(as.symbol(comparison.field)) == alt.level)%>%as.data.table()
  d[, (comparison.field) := factor(get(comparison.field), levels = c(ref.level, alt.level)) ] # re-order the factor
  
  results<-data.table(logFC = rep(NA, length(unique(d[,get(element.field)]))), # initialize a table for the results
                      pval = rep(NA, length(unique(d[,get(element.field)]))))
  
  results[,(element.field) := d[,unique(d[,get(element.field)])]] # label each row for each element comparison
  
  for (i in 1:nrow(results)){
    
    # filter down to data for a particular element
    d.element<-filter(d, UQ(as.symbol(element.field)) == results[,get(element.field)][i])
    
    # make sure there are at least two observations of the ref and alt labels
    if(nrow(filter(d.element, UQ(as.symbol(comparison.field)) == ref.level)) <  2| nrow(filter(d.element, UQ(as.symbol(comparison.field)) == alt.level)) < 2){next}
    
    tryCatch({
      res<-lm(model, d.element)
      results$logFC[i]<- summary(res)$coefficients[result.index, 1] # pull's the relative logFC/Effect size and the pvalue
      results$pval[i]<- summary(res)$coefficients[result.index, 4] # could provide option to select the index (2 by default) to select which term to record effect for
    },
    # case that only one replicate for a given term is measured for an Element
    error = function(err){print(paste("Error:  ", err))})
  }
  
  if(remove.na){results<-results[complete.cases(results),]}
  
  results[,fdr:=p.adjust(pval, method = 'BH')]
  results[,bonferroni:=p.adjust(pval, method = 'bonferroni')]
  
  # re-order the table columns
  results<-results%>%subset(select = c(element.field, 'logFC', 'pval', 'fdr', 'bonferroni'))%>%as.data.table()
  return(results)
}

#Here also is a function of just linear modeling with interaction, no random effects, very similar to the above
mpratools.lm.int<-function(d, # data.table with log-ratio activitiy values for fitting to a specified linear mixed model
                           model = Expression ~ Allele*SpliceStatus, # linear mixed model formula
                           element.field = 'Element', # field designating within what variable the group comparison should be made
                           comparison.field = 'Allele',  # field designating the grouping variable
                           ref.level = 'ref',  # level within the comparison field denoting the reference group for the comparison, level set to '0'
                           alt.level = 'alt', # level within the comparison field denoting the alternative group for the comparison, level set to '1'
                           result.index.allele = 2, # which row in the summary table should the effect size and pvalue be taken
                           result.index.secondterm = 3, # which row in the summary table should the effect size and pvalue for the second term be taken
                           result.index.interaction = 4) # which row in the summary table should the effect size and pvalue for the interaction be taken
{
  
  
  
  
  # filter provided data to just the two factors for comparison, typically by Allele
  d<-filter(d, UQ(as.symbol(comparison.field)) == ref.level | UQ(as.symbol(comparison.field)) == alt.level)%>%as.data.table()
  d[, (comparison.field) := factor(get(comparison.field), levels = c(ref.level, alt.level)) ] # re-order the factor
  #For the second interaction, it will be alphabetical unless you make them factors and reorder them ahead of time
  #(may not be a bad idea to do anyways) - this will matter when you grab logFC to make it go the right way
  results<-data.table(logFC_allele = rep(NA, length(unique(d[,get(element.field)]))), # initialize a table for the results
                      pval_allele = rep(NA, length(unique(d[,get(element.field)]))),
                      logFC_secondterm = rep(NA, length(unique(d[,get(element.field)]))),
                      pval_secondterm = rep(NA, length(unique(d[,get(element.field)]))),
                      logFC_interaction = rep(NA, length(unique(d[,get(element.field)]))),
                      pval_interaction = rep(NA, length(unique(d[,get(element.field)])))
  )
  
  results[,(element.field) := d[,unique(d[,get(element.field)])]] # label each row for each element comparison
  
  for (i in 1:nrow(results)){
    
    # filter down to data for a particular element
    d.element<-filter(d, UQ(as.symbol(element.field)) == results[,get(element.field)][i])
    
    # make sure there are at least two observations of the ref and alt labels
    if(nrow(filter(d.element, UQ(as.symbol(comparison.field)) == ref.level)) <  2 | nrow(filter(d.element, UQ(as.symbol(comparison.field)) == alt.level)) < 2){
      warning(paste('Not enough levels to compare for element:', results[,get(element.field)][i], sep = ' ')); next}
    
    tryCatch({
      res<-lm(model, d.element)
      results$logFC_allele[i]<- summary(res)$coefficients[result.index.allele, 1] # pull the relative allele effect size 
      results$pval_allele[i]<- summary(res)$coefficients[result.index.allele, 4] # allele p-val
      results$logFC_secondterm[i]<- summary(res)$coefficients[result.index.secondterm, 1] # second term logFC
      results$pval_secondterm[i]<- summary(res)$coefficients[result.index.secondterm, 4] # second term pval
      results$logFC_interaction[i]<- summary(res)$coefficients[result.index.interaction, 1] # interaction logFC
      results$pval_interaction[i]<- summary(res)$coefficients[result.index.interaction, 4] # interaction pval
      
    },
    # case that only one replicate for a given term is measured for an Element
    error = function(err){print(paste("Error:  ", err))})
    
  }
  
  results[,fdr_allele:=p.adjust(pval_allele, method = 'BH')]
  results[,bonferroni_allele:=p.adjust(pval_allele, method = 'bonferroni')]
  results[,fdr_secondterm:=p.adjust(pval_secondterm, method = 'BH')]
  results[,bonferroni_secondterm:=p.adjust(pval_secondterm, method = 'bonferroni')]
  results[,fdr_interaction:=p.adjust(pval_interaction, method = 'BH')]
  results[,bonferroni_interaction:=p.adjust(pval_interaction, method = 'bonferroni')]
  
  
  return.fields <- c(element.field, 'logFC_allele', 'pval_allele', 'fdr_allele', 'bonferroni_allele', 'logFC_secondterm', 'pval_secondterm', 'fdr_secondterm', 'bonferroni_secondterm','logFC_interaction', 'pval_interaction', 'fdr_interaction', 'bonferroni_interaction')
  
  results<-results%>%subset(select = return.fields)%>%as.data.table()
  
  return(results)
}
##########################################################################################################################################################
# T-test for Alt-Ref Hypothesis Testing
# 
# Performs a T-test between two comparison levels for some field for the provided data.table
# * Set which field denotes the elements for which a comparison will be made across groups
# * Set the field denoting the comparison groups and which levels will serve as the ref and alt group
# * Set the field for the numerical value used in the comparison. Assumes the values provided are log2 transformed
# * Select T/F if the t-test should be paired
# * Selected the field denoting replicates for the paired t-test
# If the provided data has less than 2 observations for either the ref or the alt, the comparison
# for that element will be excluded from the final results. If the t-test is to be paired, only
# replicates common to both the ref and alt will be used for the comparison, in which case, if 
# there are less than 2 replicates common to both ref and alt, the element will be excluded.
# 
# Returns a results table with a columns:
# * FC -- gives the log2 fold change for the effect size of the comparison. (we assume values provided are already log transformed)
# * pval -- nominal pvalue for the comparison
# * fdr -- Benjamini-Hochberg adjusted pvalues, obtained using p.adjust function.
# * bonferroni -- bonferroni adjusted pvalues, obtained using p.adjust function.

t.test.by.group<-function(d,  # data table containing calculated log-ratio activity measurements
                          element.field = 'Element', # field designating within what variable the group comparison should be made
                          comparison.field = 'Allele',  # field designating the grouping variable
                          ref.level = 'ref',  # level within the comparison field denoting the reference group for the comparison, level set to '0'
                          alt.level = 'alt', # level within the comparison field denoting the alternative group for the comparison, level set to '1'
                          value.field = 'Expression', # field denoting the value to be compared in the t-test
                          paired = T, # run the t-test as paired or un-paired, True by default
                          replicate.field = 'BiolRep', # what replicates will be compared across 
                          remove.na = T, # removes nan results from final table, prior to p-value correction.
                          show.warnings = T){  # if too few observations for performing test, warning is printed, othewise silent
  
  # filter provided data to just the two levels for comparison, typically by Allele
  d<-filter(d, UQ(as.symbol(comparison.field)) == ref.level | UQ(as.symbol(comparison.field)) == alt.level)%>%as.data.table()
  d[, (comparison.field) := factor(get(comparison.field), levels = c(ref.level, alt.level)) ] # re-order the factor, set ref.level to be reference
  d[, (replicate.field) := factor(get(replicate.field)) ] # make the replicate level a factor
  
  results<-data.table(logFC   = rep(NA, length(unique(d[,get(element.field)]))), # initialize a table for the results
                      pval = rep(NA, length(unique(d[,get(element.field)]))))
  
  results[,(element.field) := d[,unique(d[,get(element.field)])]] # label each row for each element comparison
  
  for (i in 1:nrow(results)){
    # filter the ref and alt data for a particular element
    d.ref <- filter(d, UQ(as.symbol(element.field)) == results[,get(element.field)][i] & UQ(as.symbol(comparison.field)) == ref.level)%>%as.data.table()
    d.alt <- filter(d, UQ(as.symbol(element.field)) == results[,get(element.field)][i] & UQ(as.symbol(comparison.field)) == alt.level)%>%as.data.table()
    
    # if using a paired t-test
    if(paired){
      common.reps <- base::intersect(d.ref[, get(replicate.field)], d.alt[, get(replicate.field)]) # see which measurable replicates intersect ref and alt
      d.ref <- filter(d.ref, UQ(as.symbol(replicate.field))%in%common.reps)%>%as.data.table() # filter ref and alt data to the shared replicates
      d.alt <- filter(d.alt, UQ(as.symbol(replicate.field))%in%common.reps)%>%as.data.table()
      
      setorderv(d.ref, replicate.field) # order the measurements by replicate for each ref and alt data table
      setorderv(d.alt, replicate.field)
    }
    
    # if there aren't enough observations of either ref or alt, then throw a warning and skip the t-test
    if(nrow(d.ref) <  2 | nrow(d.alt) < 2){if(show.warnings){warning(paste('Too few observations for at least one level to compare for element:', results[,get(element.field)][i], sep = ' '))}; next}
    
    # perform the t-test
    tryCatch({res<-t.test(d.ref[, get(value.field)], d.alt[, get(value.field)], paired = paired)},
             error = function(err, show.warnings){if(show.warnings){print(err, paste(as.character(results[,get(element.field)][i]), sep = '; On Element: '))}})
    
    # write the fold change and pvalue to results table
    if (paired){results$logFC[i] <- as.numeric(-1 * res$estimate[1])} # t.test only gives a single estimate for difference of means
    else{results$logFC[i] <- as.numeric(res$estimate[2] - res$estimate[1])}
    results$pval[i] <- res$p.value
  }
  
  # drop the nan results if remove.na == T
  if(remove.na){results<-results[complete.cases(results),]}
  
  # generate fdr and bonferroni corrected pvalues
  results[,fdr:=p.adjust(pval, method = 'BH')]
  results[,bonferroni:=p.adjust(pval, method = 'bonferroni')]
  
  # re-organize the table
  results<-results%>%subset(select = c(element.field, 'logFC', 'pval', 'fdr', 'bonferroni'))%>%as.data.table()
  # return results
  return(results)
}

##########################################################################################################################################################
# Wilcoxon Rank Sum and Rank Signed Test for Allelic/Categorical Effects
# 
# DESCRIPTION
# 
# Returns a results table with a columns:
# * FC -- gives the log2 fold change for the effect size of the comparison. (we assume values provided are already log transformed)
# * pval -- nominal pvalue for the comparison
# * fdr -- Benjamini-Hochberg adjusted pvalues, obtained using p.adjust function.
# * bonferroni -- bonferroni adjusted pvalues, obtained using p.adjust function.

wilcox.test.by.group<-function(d,  # data table containing calculated log-ratio activity measurements
                               element.field = 'Element', # field designating within what variable the group comparison should be made
                               comparison.field = 'Allele',  # field designating the grouping variable
                               ref.level = 'ref',  # level within the comparison field denoting the reference group for the comparison, level set to '0'
                               alt.level = 'alt', # level within the comparison field denoting the alternative group for the comparison, level set to '1'
                               value.field = 'Expression', # field denoting the value to be compared in the t-test
                               paired = T, # run the t-test as paired or un-paired, True by default
                               replicate.field = 'BiolRep', # what replicates will be compared across 
                               remove.na = T, # removes nan results from final table, prior to p-value correction.
                               show.warnings = T){  # if too few observations for performing test, warning is printed, othewise silent
  
  # filter provided data to just the two levels for comparison, typically by Allele
  d<-filter(d, UQ(as.symbol(comparison.field)) == ref.level | UQ(as.symbol(comparison.field)) == alt.level)%>%as.data.table()
  d[, (comparison.field) := factor(get(comparison.field), levels = c(ref.level, alt.level)) ] # re-order the factor, set ref.level to be reference
  d[, (replicate.field) := factor(get(replicate.field)) ] # make the replicate level a factor
  
  results<-data.table(logFC   = rep(NA, length(unique(d[,get(element.field)]))), # initialize a table for the results
                      pval = rep(NA, length(unique(d[,get(element.field)]))))
  
  results[,(element.field) := d[,unique(d[,get(element.field)])]] # label each row for each element comparison
  
  for (i in 1:nrow(results)){
    # filter the ref and alt data for a particular element
    d.ref <- filter(d, UQ(as.symbol(element.field)) == results[,get(element.field)][i] & UQ(as.symbol(comparison.field)) == ref.level)%>%as.data.table()
    d.alt <- filter(d, UQ(as.symbol(element.field)) == results[,get(element.field)][i] & UQ(as.symbol(comparison.field)) == alt.level)%>%as.data.table()
    
    # if using a paired t-test
    if(paired){
      common.reps <- base::intersect(d.ref[, get(replicate.field)], d.alt[, get(replicate.field)]) # see which measurable replicates intersect ref and alt
      d.ref <- filter(d.ref, UQ(as.symbol(replicate.field))%in%common.reps)%>%as.data.table() # filter ref and alt data to the shared replicates
      d.alt <- filter(d.alt, UQ(as.symbol(replicate.field))%in%common.reps)%>%as.data.table()
      
      setorderv(d.ref, replicate.field) # order the measurements by replicate for each ref and alt data table
      setorderv(d.alt, replicate.field)
    }
    
    # if there aren't enough observations of either ref or alt, then throw a warning and skip the t-test
    if(nrow(d.ref) <  2 | nrow(d.alt) < 2){if(show.warnings){warning(paste('Too few observations for at least one level to compare for element:', results[,get(element.field)][i], sep = ' '))}; next}
    
    # perform the t-test
    tryCatch({res<-wilcox.test(d.ref[, get(value.field)], d.alt[, get(value.field)], paired = paired)}, # does rank signed if paired, else rank sum
             error = function(err, show.warnings){if(show.warnings){print(err, paste(as.character(results[,get(element.field)][i]), sep = '; On Element: '))}})
    
    # write the fold change and pvalue to results table
    results$logFC[i] <- mean(d.alt[, get(value.field)]) - mean(d.ref[, get(value.field)]) # t.test only gives a single estimate for difference of means
    results$pval[i] <- res$p.value
  }
  
  # drop the nan results if remove.na == T
  if(remove.na){results<-results[complete.cases(results),]}
  
  # generate fdr and bonferroni corrected pvalues
  results[,fdr:=p.adjust(pval, method = 'BH')]
  results[,bonferroni:=p.adjust(pval, method = 'bonferroni')]
  
  # re-organize the table
  results<-results%>%subset(select = c(element.field, 'logFC', 'pval', 'fdr', 'bonferroni'))%>%as.data.table()
  # return results
  return(results)
}

##########################################################################################################################################################
# Run MPRAnalyze

# Explicitly sets the 'fraction' to select 'rna' from 'dna' by a given label
# Explicity removes shuffled sequences
# Explicity uses Counts

run_mpranalyze<-function(data, 
                         rna.label,
                         dna.label,
                         dna.levels = c('Allele', 'BC'), 
                         rna.levels = c('Allele', 'BC', 'BiolRep'), 
                         dna.model = ~ Allele + BC, # cases where only a single DNA replicate is provided, use BiolRep as DNA predictor
                         rna.model = ~ Allele,
                         comparison.field = 'Allele', 
                         ref.label = 'ref', 
                         alt.label = 'alt'){
  
  # filter provided data to just the two factors for comparison, typically by Allele
  data<-filter(data, UQ(as.symbol(comparison.field)) == ref.label | UQ(as.symbol(comparison.field)) == alt.label)%>%as.data.table()
  data[, (comparison.field) := factor(get(comparison.field), levels = c(ref.label, alt.label)) ] # re-order the factor
  
  dna.tmp<-filter(data, Fraction == dna.label)%>%subset(select = c('Element', dna.levels, 'Counts'))%>%as.data.table()
  rna.tmp<-filter(data, Fraction == rna.label)%>%subset(select = c('Element', rna.levels, 'Counts'))%>%as.data.table()
  
  element_comparison_list<-intersect(filter(dna.tmp, UQ(as.symbol(comparison.field)) == ref.label)$Element, 
                                     filter(dna.tmp, UQ(as.symbol(comparison.field)) == alt.label)$Element)
  
  dna.tmp<-filter(dna.tmp, Element%in%element_comparison_list)
  rna.tmp<-filter(rna.tmp, Element%in%element_comparison_list)
  
  dna.annot<-dna.tmp%>%subset(select = dna.levels)%>%unique()%>%as.data.table()
  rna.annot<-rna.tmp%>%subset(select = rna.levels)%>%unique()%>%as.data.table()
  
  dna.counts<-unique(dna.tmp$Element)
  rna.counts<-unique(rna.tmp$Element)
  
  for (i in 1:nrow(dna.annot)){
    tmp<-dna.tmp
    for (l in dna.levels){tmp<-filter(tmp, UQ(as.symbol(l)) == dna.annot[,..l][i][[1]])}
    dna.counts<-cbind(dna.counts, as.numeric(tmp$Counts))
  }
  rownames(dna.counts) <- dna.counts[,1]
  dna.counts<-dna.counts[,-1]%>%as.data.table()
  dna.counts[, names(dna.counts) := lapply(.SD, as.numeric)]
  dna.counts<-as.matrix(dna.counts)
  rownames(dna.counts)<-element_comparison_list
  
  for (i in 1:nrow(rna.annot)){
    tmp<-rna.tmp
    for (l in rna.levels){tmp<-filter(tmp, UQ(as.symbol(l)) == rna.annot[,..l][i][[1]])}
    rna.counts<-cbind(rna.counts, as.numeric(tmp$Counts))
  }
  rownames(rna.counts) <- rna.counts[,1]
  rna.counts<-rna.counts[,-1]%>%as.data.table()
  rna.counts[, names(rna.counts) := lapply(.SD, as.numeric)]
  rna.counts<-as.matrix(rna.counts)
  rownames(rna.counts)<-element_comparison_list
  
  obj <- MpraObject(dnaCounts = dna.counts, rnaCounts = rna.counts, 
                    dnaAnnot = dna.annot, rnaAnnot = rna.annot)
  
  obj <- estimateDepthFactors(obj, lib.factor = dna.levels, which.lib = "dna", depth.estimator = "uq")
  obj <- estimateDepthFactors(obj, lib.factor = rna.levels, which.lib = "rna", depth.estimator = "uq")
  
  obj <- analyzeQuantification(obj = obj, dnaDesign = dna.model, rnaDesign = rna.model)
  
  # reduced design is compared to rnaDesign, essentially testing whether the condition factor matters
  obj <- analyzeComparative(obj = obj, dnaDesign = dna.model, rnaDesign = rna.model, reducedDesign = ~ 1) 
  res <- testLrt(obj)
  
  element.ids <- rownames(res)
  
  res <- as.data.table(res)
  
  res[, Element := element.ids]
  res[,bonferroni:=p.adjust(pval, method = 'bonferroni')]
  
  return(res)
}

##########################################################################################################################################################
# Run mpralm

run_mpralm <- function(d,
                       element.field = 'Element',
                       bc.field = 'BC',
                       fraction.field = 'Fraction',
                       rna.level = 'rna.input',
                       dna.level = 'dna.aggmaxi',
                       comparison.field = 'Allele',
                       comparison.levels = c('ref','alt'),
                       replicate.field = 'BiolRep',
                       value.field = 'Counts',
                       aggregate.method = 'sum',
                       show.plot = F){
  
  rna.counts <- filter(d, UQ(as.symbol(fraction.field)) == rna.level)%>%subset(select = c(element.field, comparison.field, bc.field, replicate.field, value.field))%>%as.data.table()
  dna.counts <- filter(d, UQ(as.symbol(fraction.field)) == dna.level)%>%subset(select = c(element.field, comparison.field, bc.field, replicate.field, value.field))%>%as.data.table()
  
  rna.reps <- as.data.table(filter(d, UQ(as.symbol(fraction.field)) == rna.level))[,get(replicate.field)]%>%unique()
  dna.reps <- as.data.table(filter(d, UQ(as.symbol(fraction.field)) == dna.level))[,get(replicate.field)]%>%unique()
  
  if (all(dna.reps != rna.reps) & length(dna.reps) == 1){
    tmp <- c()
    for (i in rna.reps){
      dna.counts[,(replicate.field) := rep(i, nrow(dna.counts))]
      tmp <- rbind(tmp, dna.counts)
    }
    dna.counts <- tmp; rm(tmp)
  }else if (all(dna.reps != rna.reps) | (length(dna.reps) != length(rna.reps))){
    stop('The DNA replicate levels do not match the RNA replicate levels. Please provide paired DNA and RNA replicates or a single DNA replicate for all RNA replicates.')
  }
  
  rna.ref.counts <- filter(rna.counts, UQ(as.symbol(comparison.field)) == comparison.levels[1])
  rna.alt.counts <- filter(rna.counts, UQ(as.symbol(comparison.field)) == comparison.levels[2])
  dna.ref.counts <- filter(dna.counts, UQ(as.symbol(comparison.field)) == comparison.levels[1])
  dna.alt.counts <- filter(dna.counts, UQ(as.symbol(comparison.field)) == comparison.levels[2])
  
  rna.ref.counts<-spread(rna.ref.counts, (replicate.field), (value.field))%>%as.data.table()
  rna.alt.counts<-spread(rna.alt.counts, (replicate.field), (value.field))%>%as.data.table()
  dna.ref.counts<-spread(dna.ref.counts, (replicate.field), (value.field))%>%as.data.table()
  dna.alt.counts<-spread(dna.alt.counts, (replicate.field), (value.field))%>%as.data.table()
  
  eids <- base::intersect(dna.ref.counts[, get(element.field)], dna.alt.counts[, get(element.field)])
  eids <- as.data.table(filter(dna.ref.counts, UQ(as.symbol(element.field))%in%eids))[,get(element.field)]%>%as.character()
  
  rnames <- filter(dna.ref.counts, UQ(as.symbol(element.field))%in%eids)%>%tidyr::unite(ROWNAME, c(element.field, bc.field))
  rnames <- rnames$ROWNAME
  
  rna.ref.counts <- rna.ref.counts%>%filter(UQ(as.symbol(element.field))%in%eids)%>%subset(select = colnames(.)[!colnames(.)%in%c(element.field, bc.field, comparison.field)])
  rna.alt.counts <- rna.alt.counts%>%filter(UQ(as.symbol(element.field))%in%eids)%>%subset(select = colnames(.)[!colnames(.)%in%c(element.field, bc.field, comparison.field)])
  dna.ref.counts <- dna.ref.counts%>%filter(UQ(as.symbol(element.field))%in%eids)%>%subset(select = colnames(.)[!colnames(.)%in%c(element.field, bc.field, comparison.field)])
  dna.alt.counts <- dna.alt.counts%>%filter(UQ(as.symbol(element.field))%in%eids)%>%subset(select = colnames(.)[!colnames(.)%in%c(element.field, bc.field, comparison.field)])
  
  colnames(rna.ref.counts) <- paste0(paste0('replicate_', colnames(rna.ref.counts)), paste('_', comparison.levels[1], sep = ''))
  colnames(rna.alt.counts) <- paste0(paste0('replicate_', colnames(rna.alt.counts)), paste('_', comparison.levels[1], sep = ''))
  colnames(dna.ref.counts) <- paste0(paste0('replicate_', colnames(dna.ref.counts)), paste('_', comparison.levels[1], sep = ''))
  colnames(dna.alt.counts) <- paste0(paste0('replicate_', colnames(dna.alt.counts)), paste('_', comparison.levels[1], sep = ''))
  
  rownames(rna.ref.counts) <- rnames
  rownames(rna.alt.counts) <- rnames
  rownames(dna.ref.counts) <- rnames
  rownames(dna.alt.counts) <- rnames
  
  rna <- cbind(rna.ref.counts, rna.alt.counts)
  dna <- cbind(dna.ref.counts, dna.alt.counts)
  
  design <- data.frame(intcpt = 1, TMP = c(rep(F, length(rna.reps)), rep(T, length(rna.reps))))%>%setnames('TMP', comparison.levels[2])
  block_vector <- rep(c(1:length(rna.reps)), length(comparison.levels))
  
  mpra_object <- MPRASet(DNA = dna, RNA = rna, eid = eids, eseq = NULL, barcode = NULL)
  
  mpralm_allele_fit <- mpralm(object = mpra_object, design = design, aggregate = aggregate.method, normalize = TRUE, block = block_vector, model_type = "corr_groups", plot = show.plot)
  
  results <- topTable(mpralm_allele_fit, coef = 2, number = Inf)
  element.names <- rownames(results)
  results<-subset(results, select = c('logFC','P.Value'))%>%setnames('P.Value','pval')%>%as.data.table()
  results[,(element.field) := element.names]
  results[,fdr:=p.adjust(pval, method = 'BH')]
  results[,bonferroni:=p.adjust(pval, method = 'bonferroni')]
  results<-subset(results, select = c(element.field, 'logFC','pval','fdr','bonferroni'))
  
  return(results)
}

##########################################################################################################################################################
# Volcano plotting

# updated 10-15-2020: color and added effect size thresholds
# updated 12Dec2020 by Tony with an option to color points by
# local plotting density - can help to display local density of points
# and mitigate/visualize overplotting, nice visual
# default is to color by significance, to change it to density change
# argument to color.by="density"
volcano.plot <- function(results, sig.field = 'fdr', sig.threshold = 0.05, effect.field = 'logFC', effect.threshold = 0,title = '', subtitle = '', xlabel = 'log2 FC', ylabel = '-log10 FDR', show.sig.names = T, names.list = c(), color.by = 'significance'){
  
  results <- copy(results)
  results<-separate(results, Element, c('Name','varset'), sep='_')
  results <- as.data.table(results)
  results[,significance := get(sig.field) < sig.threshold & abs(get(effect.field)) >= effect.threshold]
  results[,name_listed:=sapply(Name, function(x){if(x%in%names.list){T}else{F}})]
  results$plotcolor <- "A"
  results$plotcolor[results[[sig.field]]<sig.threshold & results[[effect.field]] < -(effect.threshold)] <- "B"
  results$plotcolor[results[[sig.field]]<sig.threshold & results[[effect.field]] > effect.threshold] <- "C"
  
 #max.x <- max(c(max(results$logFC,-min(results$logFC))))
  max.x <- max(c(max(results[,get(effect.field)],-min(results[,get(effect.field)]))))
  
  
  if (length(names.list)>0 & color.by=="sigandeffect"){
    p<-ggplot(results, aes(x = get(effect.field), y = -log10(get(sig.field)),color=plotcolor)) +
      geom_point() + scale_color_manual(values=c("#B9B9B9","#AA4499","#117733"),labels=c("Not Significant","Significantly Downregulated","Significantly Upregulated"))+
      geom_hline(yintercept = -log10(sig.threshold), linetype="dashed", color = "black") +
      ggtitle(title,subtitle) + xlim(-(max.x),max.x) +
      xlab(xlabel) +
      ylab(ylabel) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5),legend.title.align=0.5) + 
      geom_text_repel(aes(label=ifelse(show.sig.names & significance, as.character(Name),'')), hjust=0.5, vjust=0,show.legend = F)+
      labs(color="Significance")+
      geom_text_repel(aes(label=ifelse(name_listed,as.character(Name),'')),hjust=0.5, vjust=0)
  } else{}
  
   if (length(names.list) > 0 & color.by=="significance"){
      p<-ggplot(results, aes(x = get(effect.field), y = -log10(get(sig.field)), color = name_listed)) +
        geom_point() +
        geom_hline(yintercept = -log10(sig.threshold), linetype="dashed", color = "black") +
        ggtitle(title,subtitle) +
        xlab(xlabel) +
        ylab(ylabel) +
        theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
        geom_text_repel(aes(label=ifelse(name_listed,as.character(Name),'')),hjust=0.5, vjust=0) + scale_color_manual(values=c("#9C9C9C", "#C637FF"))
    } else{}
  
  if (length(names.list)==0 & color.by=="sigandeffect"){
    p<-ggplot(results, aes(x = get(effect.field), y = -log10(get(sig.field)),color=plotcolor)) +
      geom_point() + scale_color_manual(values=c("#B9B9B9","#AA4499","#117733"),labels=c("Not Significant","Significantly Downregulated","Significantly Upregulated"))+
      geom_hline(yintercept = -log10(sig.threshold), linetype="dashed", color = "black") +
      ggtitle(title,subtitle) + xlim(-(max.x),max.x) +
      xlab(xlabel) +
      ylab(ylabel) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5),legend.title.align=0.5) + 
      geom_text_repel(aes(label=ifelse(show.sig.names & significance, as.character(Name),'')), hjust=0.5, vjust=0,show.legend = F)+
      labs(color="Significance")
  } else{}
  
  if (length(names.list)==0 & color.by=="density"){
    p<-ggplot(results, aes(x = get(effect.field), y = -log10(get(sig.field)))) +
      geom_pointdensity() + scale_color_viridis_c("Local\nPlotting\nDensity\n(Nearest Neighbors)",option="plasma") +
      geom_hline(yintercept = -log10(sig.threshold), linetype="dashed", color = "black") +
      ggtitle(title,subtitle) + xlim(-(max.x),max.x) +
      xlab(xlabel) +
      ylab(ylabel) +
      #scale_color_continuous("Local\nPlotting\nDensity",breaks=c(min(results$density),max(results$density)), labels = c("Low", "High"),type="viridis") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5),legend.title.align=0.5) + 
      geom_text_repel(aes(label=ifelse(show.sig.names & significance, as.character(Name),'')), hjust=0.5, vjust=0)
  } else{}
  
  if (length(names.list)>0 & color.by=="density"){
    p<-ggplot(results, aes(x = get(effect.field), y = -log10(get(sig.field)))) +
      geom_pointdensity() + scale_color_viridis_c("Local\nPlotting\nDensity\n(Nearest Neighbors)",option="plasma") +
      geom_hline(yintercept = -log10(sig.threshold), linetype="dashed", color = "black") +
      ggtitle(title,subtitle) + xlim(-(max.x),max.x) +
      xlab(xlabel) +
      ylab(ylabel) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5),legend.title.align=0.5) + 
      geom_text_repel(aes(label=ifelse(name_listed,as.character(Name),'')),hjust=0.5, vjust=0) + geom_text_repel(aes(label=ifelse(show.sig.names & significance, as.character(Name),'')), hjust=0, vjust=0)
  } else{}
  
  if (length(names.list)==0 & color.by=="significance"){
  p<-ggplot(results, aes(x = get(effect.field), y = -log10(get(sig.field)),color=significance)) +
      geom_point() +
      geom_hline(yintercept = -log10(sig.threshold), linetype="dashed", color = "black") +
      ggtitle(title,subtitle) + xlim(-(max.x),max.x) +
      xlab(xlabel) +
      ylab(ylabel) +
      theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5)) +
      geom_text_repel(aes(label=ifelse(show.sig.names & significance, as.character(Name),'')), hjust=0.5, vjust=0) + scale_color_manual(values=c("#9C9C9C", "#C637FF"))
  } else{}
  
  
  if (effect.threshold != 0){
    p <- p + geom_vline(xintercept = effect.threshold, linetype="dashed", color = "black") +
      geom_vline(xintercept = -(effect.threshold), linetype="dashed", color = "black")}

  return(p)
}

##########################################################################################################################################################
plot.venn.diagram <- function(res.list, names.list, filename){
  venn.diagram(x = res.list, category.names = names.list, filename = filename, 
               output = TRUE ,
               imagetype="png" ,
               height = 480 , 
               width = 480 , 
               resolution = 300, compression = "lzw",
               lwd = 1,
               cex = 0.5,
               fontfamily = "sans",
               cat.cex = 0.3,
               cat.default.pos = "outer",
               cat.fontfamily = "sans")
}

##########################################################################################################################################################
# edgeR MDS Plot
# Plots a given set of counts using the edgeR MDS function
# By default, counts are normalized using calcNormFactors
# Samples are grouped by given sample.fields, into measurements for elements specified by element.fields

mpratools.plot.mds <- function(d, # data.table, MPRA counts table, with each row as a single barcode measure within a distinct experiment, replicate
                               sample.fields = c('Experiment','Fraction','BiolRep','TechRep'), # vector, set of fields which distinguish replicates
                               element.fields = c('Element','Allele','BC'), # vector, set of fields which distinguish each element
                               count.field = 'Counts', # field containing count values
                               labels = c(), # list of sample labels
                               plot.title = '',
                               normalize.counts = T, 
                               sep_pattern = '::'){ # used to keep sample and element information separate through intermediate steps. 
  # should not be present in sample or element names
  
  d.counts <- unite(d, 'Sample', all_of(sample.fields), sep = sep_pattern)           # join columns for sample information into single column
  d.counts <- unite(d.counts, 'Element', all_of(element.fields), sep = sep_pattern)  # join columns for element information into single column
  
  d.counts <- subset(d.counts, select = c('Element','Sample',count.field)) # subset to just sample, element, and count information
  
  d.counts <- spread(d.counts, Sample, sym(count.field)) # spread into an N element x M Sample matrix
  
  d.counts.matrix <- as.matrix(d.counts[,c(-1)]) # convert to matrix
  
  row.names(d.counts.matrix) <- d.counts[,1][[1]] # set row names to element names
  
  y <- DGEList(counts = d.counts.matrix) # convert to a edgeR DGEList object
  
  if (normalize.counts){ y <- calcNormFactors(y); y <- cpm(y)} # calculate the normalization factors, obtain normalized counts matrix
  
  if (length(labels) > 0){ plotMDS(y, labels = labels, main = plot.title)} else {plotMDS(y, main = plot.title)}
}

##########################################################################################################################################################
# Plot PCA
# Produces a PCA graph for a provided set of counts or other metric using the native R prcomp function
# Adapted from StatQuest (Josh Starmer)
# Samples are grouped by given sample.fields, into measurements for elements specified by element.fields

mpratools.plot.pca <- function(d, # data.table, MPRA counts table, with each row as a single barcode measure within a distinct experiment, replicate
                               sample.fields = c('Experiment','Fraction','BiolRep','TechRep'), # vector, set of fields which distinguish replicates
                               element.fields = c('Element','Allele','BC'), # vector, set of fields which distinguish each element
                               count.field = 'Counts', # field containing count values
                               labels = c(), # list of sample labels
                               plot.title = '',
                               normalize.counts = T, 
                               sep_pattern = '::', # used to keep sample and element information separate through intermediate steps
                               pc_x = 1,
                               pc_y = 2,
                               return_pca = F){ # boolean, will return the pca object from prcomp if true
  # should not be present in sample or element names
  
  d.counts <- unite(d, 'Sample', all_of(sample.fields), sep = sep_pattern)           # join columns for sample information into single column
  d.counts <- unite(d.counts, 'Element', all_of(element.fields), sep = sep_pattern)  # join columns for element information into single column
  
  d.counts <- subset(d.counts, select = c('Element','Sample', count.field)) # subset to just sample, element, and count information
  
  d.counts <- spread(d.counts, Sample, sym(count.field)) # spread into an N element x M Sample matrix
  
  d.counts.matrix <- as.matrix(d.counts[apply(d.counts[,-1], 1, function(x) !all(x==0)),-1]) # convert to matrix
  
  row.names(d.counts.matrix) <- d.counts[apply(d.counts[,-1], 1, function(x) !all(x==0)),][[1]] # set row names to element names
  
  pca <- prcomp(t(d.counts.matrix), scale = TRUE) # run PCA on transpose, 'scale' standardizes all barcodes to have unit variance
  
  if (length(labels) > 0) {pca.data <- data.frame(Sample=labels, x = pca$x[,pc_x], y = pca$x[,pc_y])} 
  else { pca.data <- data.frame(Sample=rownames(pca$x), x = pca$x[,pc_x], y = pca$x[,pc_y]) }
  
  pca.var <- pca$sdev^2
  pca.var.per <- round(pca.var / sum(pca.var) * 100, 1)
  
  barplot(pca.var.per, xlab = "Principal Component", ylab = "Percent Variation (%)")
  
  p<-ggplot(pca.data, aes(x=x, y=y, label=Sample)) + 
    geom_text() +
    xlab(paste("PC", pc_x, " - ", pca.var.per[pc_x], "%", sep = "")) +
    ylab(paste("PC", pc_y, " - ", pca.var.per[pc_y], "%", sep = "")) +
    ggtitle(plot.title) + 
    theme_bw()
  plot(p)
  
  if(return_pca){return(pca)}
}

##########################################################################################################################################################
# Plot Hierarchical Clustering


mpratools.plot.hclust <- function(d, # data.table, MPRA counts table, with each row as a single barcode measure within a distinct experiment, replicate
                                  sample.fields = c('Experiment','Fraction','BiolRep','TechRep'), # vector, set of fields which distinguish replicates
                                  element.fields = c('Element','Allele','BC'), # vector, set of fields which distinguish each element
                                  count.field = 'Counts', # field containing count values
                                  labels = c(), # list of sample labels
                                  plot.title = '',
                                  normalize.counts = T, 
                                  sep_pattern = '::' # used to keep sample and element information separate through intermediate steps
){ # boolean, will return the pca object from prcomp if true
  # should not be present in sample or element names
  
  d.counts <- unite(d, 'Sample', all_of(sample.fields), sep = sep_pattern)           # join columns for sample information into single column
  d.counts <- unite(d.counts, 'Element', all_of(element.fields), sep = sep_pattern)  # join columns for element information into single column
  
  d.counts <- subset(d.counts, select = c('Element','Sample', count.field)) # subset to just sample, element, and count information
  
  d.counts <- spread(d.counts, Sample, sym(count.field)) # spread into an N element x M Sample matrix
  
  d.counts.matrix <- as.matrix(d.counts[apply(d.counts[,-1], 1, function(x) !all(x==0)),-1]) # convert to matrix
  
  row.names(d.counts.matrix) <- d.counts[apply(d.counts[,-1], 1, function(x) !all(x==0)),][[1]] # set row names to element names
  
  hc <- hclust(dist(t(d.counts.matrix)))
  
  plot(hc)
}
##########################################################################################################################################################
# Plot Barcode Recovery

plot.bc.recovery <- function(counts, sample.fields = c('BiolRep','Fraction'), element.fields = c('Element','Allele'),plot.title='Barcode Recovery by Sample'){
  
  d <- copy(counts)
  
  if (length(sample.fields) >  1){d <- unite(d, Sample, all_of(sample.fields), sep = '::')
  } else {d[, Sample := get(sample.fields[1])]}
  if (length(element.fields) >  1){d <- unite(d, Element, all_of(element.fields), sep = '::')
  } else {d[, Element := get(element.fields[1])]}
  
  barcode_recovery <- ddply(d, c('Sample', 'Element'), summarise, Recovered_BC = sum(Counts > 0))
  barcode_recovery <- ddply(barcode_recovery, c('Sample','Recovered_BC'), summarise, Count = length(Recovered_BC))
  barcode_recovery <- ddply(barcode_recovery, c('Sample'), transform, Percent = Count / sum(Count) * 100)
  barcode_recovery <- setnames(barcode_recovery, old = 'Recovered_BC', new = 'Number BCs Recovered')
  
  barcode_recovery <- as.data.table(barcode_recovery)
  barcode_recovery[, `Number BCs Recovered` := factor(`Number BCs Recovered`)]
  
  p <- ggplot(barcode_recovery, aes(x = Percent, y = Sample, fill = `Number BCs Recovered`)) + 
    geom_col() + 
    ggtitle(plot.title) + 
    scale_fill_brewer(palette="Blues") + 
    xlab('% of Elements') + 
    theme(text = element_text(size=12), plot.title = element_text(hjust = 0.5))
  
  plot(p)
  
  print(spread(subset(barcode_recovery, select = c('Sample','Percent','Number BCs Recovered')), Sample, Percent))
}

##########################################################################################################################################################
# Borrowed Functions
panel.cor <- function(x, y, digits=2, prefix="", cex.cor){
  
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(0, 1, 0, 1)) 
  r <- cor(x, y) 
  txt <- format(c(r, 0.123456789), digits=digits)[1] 
  txt <- paste(prefix, txt, sep="") 
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
  
  test <- cor.test(x,y) 
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")) 
  
  text(0.5, 0.5, txt, cex = cex * r) 
  text(.8, .8, Signif, cex=cex, col=2) 
}

qqunif.plot<-function(pvalues, 
                      should.thin=T, thin.obs.places=2, thin.exp.places=2, 
                      xlab=expression(paste("Expected (",-log[10], " p-value)")),
                      ylab=expression(paste("Observed (",-log[10], " p-value)")), 
                      draw.conf=TRUE, conf.points=1000, conf.col="lightgray", conf.alpha=.05,
                      already.transformed=FALSE, pch=20, aspect="iso", prepanel=prepanel.qqunif,
                      par.settings=list(superpose.symbol=list(pch=pch)), ...) {
  #error checking
  if (length(pvalues)==0) stop("pvalue vector is empty, can't draw plot")
  if(!(class(pvalues)=="numeric" || 
       (class(pvalues)=="list" && all(sapply(pvalues, class)=="numeric"))))
    stop("pvalue vector is not numeric, can't draw plot")
  if (any(is.na(unlist(pvalues)))) stop("pvalue vector contains NA values, can't draw plot")
  if (already.transformed==FALSE) {
    if (any(unlist(pvalues)==0)) stop("pvalue vector contains zeros, can't draw plot")
  } else {
    if (any(unlist(pvalues)<0)) stop("-log10 pvalue vector contains negative values, can't draw plot")
  }
  
  grp<-NULL
  n<-1
  exp.x<-c()
  if(is.list(pvalues)) {
    nn<-sapply(pvalues, length)
    rs<-cumsum(nn)
    re<-rs-nn+1
    n<-min(nn)
    if (!is.null(names(pvalues))) {
      grp=factor(rep(names(pvalues), nn), levels=names(pvalues))
      names(pvalues)<-NULL
    } else {
      grp=factor(rep(1:length(pvalues), nn))
    }
    pvo<-pvalues
    pvalues<-numeric(sum(nn))
    exp.x<-numeric(sum(nn))
    for(i in 1:length(pvo)) {
      if (!already.transformed) {
        pvalues[rs[i]:re[i]] <- -log10(pvo[[i]])
        exp.x[rs[i]:re[i]] <- -log10((rank(pvo[[i]], ties.method="first")-.5)/nn[i])
      } else {
        pvalues[rs[i]:re[i]] <- pvo[[i]]
        exp.x[rs[i]:re[i]] <- -log10((nn[i]+1-rank(pvo[[i]], ties.method="first")-.5)/(nn[i]+1))
      }
    }
  } else {
    n <- length(pvalues)+1
    if (!already.transformed) {
      exp.x <- -log10((rank(pvalues, ties.method="first")-.5)/n)
      pvalues <- -log10(pvalues)
    } else {
      exp.x <- -log10((n-rank(pvalues, ties.method="first")-.5)/n)
    }
  }
  
  #this is a helper function to draw the confidence interval
  panel.qqconf<-function(n, conf.points=1000, conf.col="gray", conf.alpha=.05, ...) {
    require(grid)
    conf.points = min(conf.points, n-1);
    mpts<-matrix(nrow=conf.points*2, ncol=2)
    for(i in seq(from=1, to=conf.points)) {
      mpts[i,1]<- -log10((i-.5)/n)
      mpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i))
      mpts[conf.points*2+1-i,1]<- -log10((i-.5)/n)
      mpts[conf.points*2+1-i,2]<- -log10(qbeta(conf.alpha/2, i, n-i))
    }
    grid.polygon(x=mpts[,1],y=mpts[,2], gp=gpar(fill=conf.col, lty=0), default.units="native")
  }
  
  #reduce number of points to plot
  if (should.thin==T) {
    if (!is.null(grp)) {
      thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                exp.x = round(exp.x, thin.exp.places),
                                grp=grp))
      grp = thin$grp
    } else {
      thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                exp.x = round(exp.x, thin.exp.places)))
    }
    pvalues <- thin$pvalues
    exp.x <- thin$exp.x
  }
  gc()
  
  prepanel.qqunif= function(x,y,...) {
    A = list()
    A$xlim = range(x, y)*1.02
    A$xlim[1]=0
    A$ylim = A$xlim
    return(A)
  }
  
  #draw the plot12-5-18_biological reps.xls
  xyplot(pvalues~exp.x, groups=grp, xlab=xlab, ylab=ylab, aspect=aspect,
         prepanel=prepanel, scales=list(axs="i"), pch=pch,
         panel = function(x, y, ...) {
           if (draw.conf) {
             panel.qqconf(n, conf.points=conf.points, 
                          conf.col=conf.col, conf.alpha=conf.alpha)
           };
           panel.xyplot(x,y, ...);
           panel.abline(0,1);
         }, par.settings=par.settings, ...
  )
}

##########################################################################################################################################################
#Enrichment Plot by logFC
#Stephen's script for plotting significant enrichment of case vs. control
#in the ref/alt expression by incremental steps in logFC
#The function takes two arguments, the logFC table generated by a statistical
#test (i.e., t-test or linear/mixed modeling), and a title for the graph
#You can also prefilter the table to only analyze elements with positive logFC
#(stabilizing) or negative logFC (destabilizing)

plot_enrichment <- function(mpra_res, sig.field, effect.field, title = ""){
  
  mpra_res <- as.data.table(mpra_res)
  
  points <- abs(filter(mpra_res, get(sig.field) < 0.05)[[effect.field]])
  
  enrichment <- c()
  pval <- c()
  n_variants <- c()
  
  mpra_res[, CASE := Pheno == 'case']
  
  for (i in points){
    
    mpra_res[, SIGNIFICANT := (get(sig.field) < 0.05)  & (abs(get(effect.field)) > i)]
    enrich_test <- test.for.enrichment(mpra_res, "SIGNIFICANT", "CASE", test.type = 'hyper.test')
    enrichment <- c(enrichment, enrich_test[2])
    n_variants <- c(n_variants, nrow(filter(mpra_res, SIGNIFICANT)))
    pval <- c(pval, enrich_test[1])
  }
  
  res <- data.table(logFC_Threshold = points, enrichment = enrichment, pval = pval, n = n_variants)
  res[, enrichment := as.double(gsub("NaN", NA, enrichment))]
  res[, significant := pval < 0.05]
  
  p1 <- ggplot(res, aes(x = logFC_Threshold, y = enrichment)) +
    geom_point() +
    geom_hline(yintercept = 1, linetype = "dashed") +
    ylab('Odds Ratio') + 
    ggtitle(title) + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  p2 <- ggplot(res, aes(x = logFC_Threshold, y = pval, color = significant)) + 
    geom_point() + 
    ylab('Hypergeometric\nP-Value') + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "none") + 
    scale_y_continuous(trans = "log10") + 
    scale_color_manual(values=c("#9C9C9C", "#379eff"))
  
  p3 <- ggplot(res, aes(x = logFC_Threshold, y = n)) + 
    geom_point() + 
    ylab('N Variants') + 
    xlab('Absolute log2FC Threshold')
  
  p <- (p1 / p2 / p3)
  
  return(p)
}

################################################################################
#Here are Tony's scripts for pulling out random effects conditional modes (means)
#from a linear mixed model
#The arguments are the same table you gave the lmm, with the same model description
#The function returns the same table (minus any shuffled alleles) plus a new column
#of conditional means for that random effect (usually barcode)
random.effects <- function(d,model=Expression ~ Allele*SpliceStatus + (1|BC)){
  filter(d,Element!="shuf")
  big_ol_data <- list()
  uniques <- unique(d$Element)
  for(i in 1:length(uniques)){
    rand <- as.data.table(ranef(lmer(formula=model,filter(d,Element==uniques[i]))))
    rand[,Element:=uniques[i]]
    big_ol_data[[i]] <- rand
  }
  big_call <- do.call(rbind,big_ol_data)
  d$ConditionalMode <- 0
  for(i in 1:nrow(big_call)){
    d$ConditionalMode[which(d$Element==big_call$Element[i] & d$BC==big_call$grp[i])] <- big_call$condval[i]
  }
  return(d)
}

################################################################################
#A function for displaying Pearson's correlation on correlation matrices using
#the ggcorrm package
corrs <- function(d, title="You forgot to title this graph.", method, rescale,
                  boxtext="I'm not sure what I'm looking at here."){
  
  corrs.text <- function(x,y,method=method){
    res <- cor.test(x,y,method=method)
    if (res$p.value<=0.0001){
      stars <- "****"
    } else
      if (res$p.value<=0.001){
        stars <- "***"
      } else
        if (res$p.value<=0.01){
          stars <- "**"
        } else
          if (res$p.value<=0.05){
            stars <- "*"
          } else
            stars <- "ns"
          return(paste(boxtext,"\n",round(res$estimate,4),stars))
  }
  
  p <- ggcorrm(d, mapping = aes(col = .corr, fill = .corr), rescale = rescale) +
    lotri(geom_smooth(method = "lm", size = .3,color="black")) +
    lotri(geom_point(alpha = 0.5)) +
    utri_funtext(fun = ~ corrs.text(.x,.y, method = method)) +
    dia_names(y_pos = 0.15, size = 2) +
    dia_histogram(bins=50,color="black") +
    dia_density(color="black",fill="transparent") +
    scale_color_corr(option="A",aesthetics = c("fill", "color"))
  
  p <- p + labs(title=title) + theme(plot.title = element_text(hjust=0.5))
  return(p)
}
################################################################################
#Here is a function to apply library 1 and 2 annotations
get1and2annotations <- function(d){
  #Split the libraries
  lib2res <- filter(d,grepl("^chr",Element))
  lib1res <- filter(d,!grepl("^chr",Element))
  
  
  libann <- fread("./lib1annotations.csv") %>% setnames(.,"seqName","Element")
  libann$Element <- str_replace(libann$Element,"_ref","")
  #Apply Library label and put it first
  libann$Library <- 1
  setcolorder(libann, c("Library", names(libann)[-ncol(libann)]))
  libann$Strand <- ifelse(libann$Strand=="+",1,-1)
  
  lib1res <- full_join(lib1res,libann[,c(1,2,3,4,5,6,9,10,11,12,14)],by="Element")
  
  lib1res$Type <- ifelse((str_length(lib1res$Ref) > 1 | str_length(lib1res$Alt) > 1), "Indel", "SNV")
  lib1res <- tidyr::separate(lib1res,hg38positionTag,"Chr",sep=";")
  lib1res <- tidyr::separate(lib1res,Chr,c("Chr","Pos"),sep=":")
  #lib1res$logFC[is.na(lib1res$logFC)] <- "Not Tested"
  #lib1res$pval[is.na(lib1res$pval)] <- "Not Tested"
  #lib1res$fdr[is.na(lib1res$fdr)] <- "Not Tested"
  #lib1res$bonferroni[is.na(lib1res$bonferroni)] <- "Not Tested"
  lib1res$Source <- str_replace_all(lib1res$Source,
                                    "Exome:SSC_Iossifov_2014",
                                    "WES_Iossifov_2014")
  lib1res$Source <- str_replace_all(lib1res$Source,
                                    "Exome:ASC_DeRubeis",
                                    "WES_DeRubeis_2014")
  lib1res$Source <- str_replace_all(lib1res$Source,
                                    "WGS",
                                    "WGS_Werling_2018")
  
  
  libann2 <- fread("./lib2annotations.csv") %>% setnames(.,"PositionTag","Element")
  libann2$Library <- 2
  setcolorder(libann2, c("Library", names(libann2)[-ncol(libann2)]))
  setnames(libann2,"Fam","familyID")
  setnames(libann2,"SYMBOL","Gene")
  setnames(libann2,"gene_strand","Strand")
  libann2$Source <- "WGS_An_2018"
  libann2$sampleID <- ifelse(libann2$Pheno=="case",
                             paste0(libann2$familyID,".p1"),
                             paste0(libann2$familyID,".s1"))
  
  lib2res <- full_join(lib2res,libann2[,c(1,2,3,4,37,8,45,13,5,6,10,44,7)],by="Element")
  
  #lib2res$logFC[is.na(lib2res$logFC)] <- "Not Tested"
  #lib2res$pval[is.na(lib2res$pval)] <- "Not Tested"
  #lib2res$fdr[is.na(lib2res$fdr)] <- "Not Tested"
  #lib2res$bonferroni[is.na(lib2res$bonferroni)] <- "Not Tested"
  
  
  results.all <- rbind(lib1res,lib2res)
  
  return(results.all)
}