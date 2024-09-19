get_metabolite_leadscope_models <- function(smi_ids){
    model_tables <- lapply(as.double(smi_ids), function(x){
        for(i in seq_len(13)){
            metabolite_models <- run_query_metabolites(paste0("
                SELECT
                    *
                FROM
                    metabolite_", i-1, "
                WHERE
                    smi_id = $1
            "), args=list(x))
            if(nrow(metabolite_models) > 0){
                return(metabolite_models)
            }
        }
        return(NULL)
    })
    names(model_tables) <- as.double(smi_ids)
    return(model_tables)
}