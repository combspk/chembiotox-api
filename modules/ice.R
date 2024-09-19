# Load required packages for querying the ICE API
library(jsonlite)

### ---------------------------------------- ###
###            Module functions              ###
### ---------------------------------------- ###
# Query ICE API for a list of chemicals
# Returns a dataframe
query_ice <- function(chemicals=c()){
    res <- tryCatch({
        body_content <- paste0('{"chemids": ["', paste0(chemicals, collapse='", "'), '"]}')
        res <- httr::POST(url="https://ice.ntp.niehs.nih.gov/api/v1/search", httr::content_type_json(), body=body_content)
        json <- fromJSON(httr::content(res, "text"))$endPoints
        
        return(json)
    }, error=function(cond){
        print(cond)
        return(data.frame(matrix(ncol=23, nrow=0)))
    })
}

aggregate_ice_assays <- function(x=data.frame()){
    if(!("assay" %in% colnames(x))){
        return(NULL)
    }
    assay_counts <- split(x, x$assay)
    assay_counts <- unlist(lapply(assay_counts, function(assay_count) length(unique(assay_count$dtxsid))))
    assay_counts <- sort(assay_counts, decreasing=TRUE)
    assay_counts <- data.frame(assay_name=names(assay_counts), assay_num=assay_counts)
    
    return(assay_counts)
}