################################################################################
#' Bin quality scores from a non-binned FASTQ file.
#' 
#' This Function takes a raw FASTQ file with unbinned quality scores and outputs
#' a new FASTQ file with binned quality scores based on the selected binning
#' scheme
#' 
#' 
#' @param infastq (Required). Character vector.
#'  Path to the raw FASTQ file with unbinned quality scores
#' 
#' @param outfastq (Required). Character vector.
#'  Path for the output binned FASTQ file
#'   
#' @param scheme (Optional). Character vector. 
#'  Default is "novaseqc1.3".
#'  Binning scheme to be used. Options include ("novaseqc1.3","novaseqc1.2")
#'   
#' @param bins (Optional). Numeric vector. 
#'  If bins and binlabs given this overrides
#'  any scheme given. If 0 is not included as the bottom of the first bin, it will be
#'  added.
#' 
#' @param binlabs (Optional) Character vector. If bins and binlabs given this
#'  overrides any scheme given. Must be n-1 of bins (including the 0 minimum)
#' 
#' 
#' @return new_fastq ShortReadQ object .
#'  ShortRead object. Contains original sreads, new binned quality scores, and
#'  original ids.
#'
#' 
#' @examples
#' input="pathtoinput"
#' output="pathtooutput"
#' binqualscores(input,output)
#' 
#' 


binqualscores<- function(infastq,outfastq,scheme="novaseqc1.3",bins=NULL,binlabs=NULL){
  #Check input file exists
  if (!file.exists(infastq)){
    stop("Cannot find input FASTQ. Ensure correct path and file name")
  }
  #Warn if user does not specify binning scheme
  if (missing(scheme) & (missing(bins)&missing(binlabs))){
    warning("No binning scheme specified. Assuming NovaSeq Control v 1.3")
  }
  #Warn user about override
  if (!missing(scheme) & (!missing(bins)&!missing(binlabs)) ) {
    warning("Custom Binning overrides preset binning schemes.")
  }
  
  # Set bins and bin labels
  
  #Check for manual binnning (overrides scheme)
  if (!missing(bins)&!missing(binlabs)){
    #Set binning scheme to custom
    scheme="CUSTOM"
    
    #check for 0 as proper starting point of bins
    if(!(0 %in% bins)){
      bins<-append(bins,0)
    }
    
    #Ensure bins are sorted smallest to largest
    bins<-sort(bins)
    
    #check number of bins against bin labels
    if(length(binlabs)!=length(bins)-1){
      stop("Incorrect amount of bin labels for given bins.")
    }
  } else {
    if(!exists(scheme,binschemes)){
      stop("Premade Binning scheme does not exist. Please choose from: novaseqc1.3,
             or enter manual bin and binning label information")
    }
    else {
      bins<-binschemes[[scheme]]$bins
      binlabs<-binschemes[[scheme]]$binlabs
    }
  }
  
  #Read in input FASTQ file
  og_fastq<- readFastq(infastq)

  #Pull quality scores into a data frame
  Quality = (quality(og_fastq))
  quality_scores<-as(Quality,"matrix")
  quality_scores<-as.data.frame(quality_scores)
  
  #Bins to match binning scheme for quality scores
  df_binned <- quality_scores %>%
    mutate(across(everything(), ~ cut(.,breaks = bins,labels = binlabs)))
  
  #Turn to character and then integer
  df_binned<- df_binned %>% mutate_all(as.character) %>% mutate_all(as.integer)
  
  #Back to matrix
  temp<-as.matrix(df_binned)
  dimnames(temp) <- NULL
  
  #Convert to string for recreating quality score object
  phred_to_string <- function(x, offset = 33L) {
    paste0(intToUtf8(x + offset, multiple = TRUE), collapse = "")
  }
  qual_strings <- apply(temp, 1, phred_to_string)
  
  #Create Quality Score Object
  temp_qual <- PhredQuality(qual_strings)
  
  #ShortRead object
  new_fastq <- ShortReadQ(
    sread   = sread(og_fastq),
    quality = temp_qual,
    id      = ShortRead:::id(og_fastq)
  )
  
  #Write out for later use
  writeFastq(new_fastq, outfastq)
  
  return(new_fastq)
}

################################################################################

#Binning Schemes

binschemes<-list("novaseqc1.3"=
                      list("bins"=c(0,2,17,29,100),
                           "binlabs"=c("2","9","24","40")),
                 "novaseqc1.2"=
                      list("bins"=c(0,2,17,29,100),
                           "binlabs"=c("2","12","24","40")
                    )
)







