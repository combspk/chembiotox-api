process_chemical_file <- function(path, delim=",", identifier="DTXSID"){
    tryCatch({
        chemical_input <- fread(path, header=FALSE, fill=TRUE, sep=delim)
        
        print("==chemical_input==")
        print(chemical_input)
        
        if(identifier == "DTXSID"){
            # Check if DTXSID, find DTXSID column
            dtxsid_col <- ""
            first_row <- chemical_input[1, ]
            for (col in colnames(first_row)){
                if(grepl("^DTXSID", first_row[1, ..col])){
                    dtxsid_col <- col
                }
            }
            
            setkeyv(chemical_input, dtxsid_col)
            
            dtxsids <- unlist(chemical_input[, ..dtxsid_col])
            matched_chemicals <- run_query(paste0("
                SELECT
                    bcc.preferred_name,
                    bcc.casrn,
                    bc.dsstox_substance_id,
                    bc.epa_id--,
                    --bcs.smi_id,
                    --ss.smiles
                FROM
                    base_chemicals bc,
                    base_chemical_compounds bcc--,
                    --base_chemical_to_smiles bcs,
                    --smiles_strings ss
                WHERE
                    bc.dsstox_substance_id IN (", paste0(lapply(seq_len(length(dtxsids)), function(x) paste0("$", x)), collapse=", "), ")
                AND bc.epa_id = bcc.epa_id
                --AND bc.epa_id = bcs.epa_id
                --AND bcs.smi_id = ss.smi_id
            "), args=as.list(unname(dtxsids)))
            
            print("==matched_chemicals==")
            print(matched_chemicals)
    
            missing_dtxsids <- setdiff(dtxsids, matched_chemicals$dsstox_substance_id)
            missing_dtxsids <- unique(chemical_input[missing_dtxsids, ])
            
            return(list(
                matched=matched_chemicals,
                missing=missing_dtxsids
            ))
        } else {
            
            # Check if CASRN, find CASRN column
            casrn_col <- ""
            first_row <- chemical_input[1, ]
            for (col in colnames(first_row)){
                if(grepl("^[0-9]+-[0-9]+-[0-9]+", first_row[1, ..col])){
                    casrn_col <- col
                }
            }
            
            setkeyv(chemical_input, casrn_col)
            
            casrns <- unlist(chemical_input[, ..casrn_col])
            matched_chemicals <- run_query(paste0("
                SELECT
                    bcc.preferred_name,
                    bcc.casrn,
                    bc.dsstox_substance_id,
                    bc.epa_id--,
                    --bcs.smi_id,
                    --ss.smiles
                FROM
                    base_chemicals bc,
                    base_chemical_compounds bcc--,
                    --base_chemical_to_smiles bcs,
                    --smiles_strings ss
                WHERE
                    bcc.casrn IN (", paste0(lapply(seq_len(length(casrns)), function(x) paste0("$", x)), collapse=", "), ")
                AND bc.epa_id = bcc.epa_id
                --AND bc.epa_id = bcs.epa_id
                --AND bcs.smi_id = ss.smi_id
            "), args=as.list(unname(casrns)))
            
            print("==matched_chemicals==")
            print(matched_chemicals)
            
            missing_casrns <- setdiff(casrns, matched_chemicals$casrn)
            missing_casrns <- unique(chemical_input[missing_casrns, ])
            
            return(list(
                matched=matched_chemicals,
                missing=missing_casrns
            ))
            
        }
        
    }, error=function(err){
        print(err)
        return(NULL)
    })
}