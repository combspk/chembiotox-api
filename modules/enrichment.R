library(data.table)
source("modules/postgres.R")

base_chembl_alerts_query <- paste0("
    SELECT DISTINCT
        chemblsa_id,
        CONCAT(alert_name, ' | ', smarts) AS alert_name
    FROM
        chembl_structural_alerts_annotated
")
base_chembl_alerts <- run_query(query=base_chembl_alerts_query)
base_chembl_alerts <- as.data.table(base_chembl_alerts)


# Convert annotation table to heatmap
annotation_heatmap <- function(data){
    if(nrow(data) < 2){
        return(NULL)
    }
    heatmaply(
        x=data,
        grid_gap=0.5,
        width=1000,
        height=1000,
        
    ) %>% layout(
        width=1000,
        height=1000
    )
    
}


# Load NTP Studies
get_ntp_studies <- function(epa_ids){
    ntp_studies <- run_query(paste0("
        SELECT DISTINCT
            bc.dsstox_substance_id,
            bcc.preferred_name,
            bcc.casrn,
            bcta.testarticle,
            bcta.ntpid,
            bcta.study_area,
            bcta.study_number,
            CONCAT('<a target=_blank rel=noopener noreferrer href=', ntau.web_url, '>View study</a>') AS url
        FROM
            base_chemical_test_articles bcta,
            ntp_test_article_urls ntau,
            base_chemicals bc,
            base_chemical_compounds bcc
        WHERE
            bcta.epa_id IN (", paste0(lapply(seq_len(length(epa_ids)), function(x) paste0("$", x, "::bigint")), collapse=","), ")
        AND bcta.ntpid = ntau.ntpid
        AND bcta.epa_id = bc.epa_id
        AND bcta.epa_id = bcc.epa_id
    "), args=epa_ids)
    if(ncol(ntp_studies) == 8){
        colnames(ntp_studies) <- c("DTXSID", "Chemical Name", "CASRN", "Article", "NTP ID", "Study Area", "Study Number", "URL")
    }
    return(ntp_studies)
}

# Load chemical availability information
get_chemical_availability <- function(epa_ids){
    chemical_availability <- run_query(paste0("
        SELECT DISTINCT
            bc.dsstox_substance_id,
            bcc.preferred_name,
            bcc.casrn,
            pca.source_name,
            --CONCAT('<a target=_blank rel=noopener noreferrer href=', pca.source_record_url, '>View at source</a>') AS url
            pca.source_record_url AS url
        FROM
            base_chemical_to_pubchem_cid bpc,
            pubchem_chemical_availability pca,
            base_chemicals bc,
            base_chemical_compounds bcc
        WHERE
            bc.epa_id IN (", paste0(lapply(seq_len(length(epa_ids)), function(x) paste0("$", x, "::bigint")), collapse=","), ")
        AND bc.epa_id = bcc.epa_id
        AND bc.epa_id = bpc.epa_id
        AND bpc.pubchem_cid = pca.pubchem_cid
    "), args=epa_ids)
    if(ncol(chemical_availability) == 5){
        
        chemical_availability$url <- lapply(chemical_availability$url, function(x){
            if(is.na(x)){
                return("No link available")
            }
            return(paste0("<a target=_blank rel=noopener noreferrer href='", x, "'>View at source</a>"))
        })
        
        colnames(chemical_availability) <- c("DTXSID", "Chemical Name", "CASRN", "Source Name", "URL")
    }
    return(chemical_availability)
}

# Load chemical annotations for a given cluster
chemical_annotations <- function(selected, input, selected_node, dtxsids, switch_mode=""){
    input_all_annotations <- list(
        "switch__admet_binr",
        "switch__admet_catg",
        "switch__admet_cont",
        "switch__chembl",
        "switch__cpd",
        "switch__ctd_bioprocess",
        "switch__ctd_cellcomp",
        "switch__ctd_diseases",
        "switch__ctd_genes",
        "switch__ctd_molfunct",
        "switch__ctd_phenotypes",
        "switch__drugbank_atccodes",
        "switch__drugbank_carriers",
        "switch__drugbank_enzymes",
        "switch__drugbank_targets",
        "switch__drugbank_transporters",
        "switch__hmdb_biospecimenlocations",
        "switch__hmdb_cellularlocations",
        "switch__hmdb_diseases",
        "switch__hmdb_genes",
        "switch__hmdb_tissuelocations",
        "switch__ice",
        "switch__invitrodb",
        "switch__leadscope",
        "switch__ochem",
        "switch__opera_adme",
        "switch__opera_environmental_fate",
        "switch__opera_physicochem",
        "switch__opera_structural_properties",
        "switch__opera_toxicology_endpoints",
        "switch__saagar",
        "switch__t3db",
        "switch__toxrefdb_neoplastic",
        "switch__toxrefdb_nonneoplastic",
        "switch__pubchem",
        "switch__pubchem_props",
        "switch__gras",
        "switch__foodb_enzymes",
        "switch__foodb_flavors",
        "switch__foodb_foodcontent",
        "switch__foodb_healtheffects",
        "switch__foodb_ontology",
        "switch__epa_props"
    )
    
    input <- lapply(input_all_annotations, function(x){
        if(x %in% input){
            return(TRUE)
        }
        return(FALSE)
    })
    names(input) <- input_all_annotations

    
    anno_admet_binr <- data.frame(matrix(ncol=3, nrow=0))
    colnames(anno_admet_binr) <- c("model_name", "interpretation", "description")
    anno_admet_catg <- data.frame(matrix(ncol=3, nrow=0))
    colnames(anno_admet_catg) <- c("model_name", "interpretation", "description")
    anno_admet_cont <- data.frame(matrix(ncol=3, nrow=0))
    colnames(anno_admet_cont) <- c("model_name", "interpretation", "description")
    
    anno_admet_binr_hm <- plotly_empty()
    anno_admet_catg_hm <- plotly_empty()
    anno_admet_cont_hm <- plotly_empty()
    
    anno_chembl <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_chembl) <- c("alert_name", "alert_num")
    anno_cpd <- data.frame(matrix(ncol=3, nrow=0))
    colnames(anno_cpd) <- c("gen_cat", "prod_fam", "prod_type")
    
    anno_ctd_bioprocess <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_ctd_bioprocess) <- c("term_name", "term_num")
    anno_ctd_cellcomp <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_ctd_cellcomp) <- c("term_name", "term_num")
    anno_ctd_diseases <- data.frame(matrix(ncol=7, nrow=0))
    colnames(anno_ctd_diseases) <- c("disease_name", "disease_id", "direct_evidence", "inference_gene_symbol", "inference_score", "omimids", "pubmed_ids")
    anno_ctd_genes <- data.frame(matrix(ncol=7, nrow=0))
    colnames(anno_ctd_genes) <- c("gene_symbol", "gene_id", "gene_forms", "organism", "interaction", "interaction_actions", "pubmed_ids")
    anno_ctd_molfunct <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_ctd_molfunct) <- c("term_name", "term_num")
    anno_ctd_phenotypes <- data.frame(matrix(ncol=9, nrow=0))
    colnames(anno_ctd_phenotypes) <- c("phenotype_name", "phenotype_id", "comentioned_items", "organism", "interaction", "interaction_actions", "anatomy_terms", "inference_gene_symbols", "pubmed_ids")

    anno_drugbank_atccodes <- data.frame(matrix(ncol=3, nrow=0))
    colnames(anno_drugbank_atccodes) <- c("atc_code", "atc_annotation", "atc_level")
    anno_drugbank_carriers <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_drugbank_carriers) <- c("annotation_name", "annotation_num")
    anno_drugbank_enzymes <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_drugbank_enzymes) <- c("annotation_name", "annotation_num")
    anno_drugbank_targets <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_drugbank_targets) <- c("annotation_name", "annotation_num")
    anno_drugbank_transporters <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_drugbank_transporters) <- c("annotation_name", "annotation_num")

    anno_hmdb_biospecimenlocations <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_hmdb_biospecimenlocations) <- c("annotation_name", "annotation_num")
    anno_hmdb_cellularlocations <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_hmdb_cellularlocations) <- c("annotation_name", "annotation_num")
    anno_hmdb_diseases <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_hmdb_diseases) <- c("annotation_name", "annotation_num")
    anno_hmdb_genes <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_hmdb_genes) <- c("annotation_name", "annotation_num")
    anno_hmdb_tissuelocations <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_hmdb_tissuelocations) <- c("annotation_name", "annotation_num")

    anno_ice <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_ice) <- c("assay_name", "assay_num")

    anno_invitrodb <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_invitrodb) <- c("assay_endpoint_name", "assay_num")
    anno_leadscope <- data.frame(matrix(ncol=3, nrow=0))
    colnames(anno_leadscope) <- c("model_name", "interpretation", "description")

    anno_ochem <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_ochem) <- c("alert_name", "alert_num")

    anno_opera_adme <- data.frame(matrix(ncol=6, nrow=0))
    colnames(anno_opera_adme) <- c("preferred_name", "dsstox_substance_id", "model_name", "model_value", "category", "model_description")

    anno_opera_environmental_fate <- data.frame(matrix(ncol=6, nrow=0))
    colnames(anno_opera_environmental_fate) <- c("preferred_name", "dsstox_substance_id", "model_name", "model_value", "category", "model_description")

    anno_opera_physicochem <- data.frame(matrix(ncol=6, nrow=0))
    colnames(anno_opera_physicochem) <- c("preferred_name", "dsstox_substance_id", "model_name", "model_value", "category", "model_description")

    anno_opera_structural_properties <- data.frame(matrix(ncol=6, nrow=0))
    colnames(anno_opera_structural_properties) <- c("preferred_name", "dsstox_substance_id", "model_name", "model_value", "category", "model_description")

    anno_opera_toxicology_endpoints <- data.frame(matrix(ncol=6, nrow=0))
    colnames(anno_opera_toxicology_endpoints) <- c("preferred_name", "dsstox_substance_id", "model_name", "model_value", "category", "model_description")

    anno_saagar <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_saagar) <- c("descript_name", "descript_num")
    anno_t3db <- data.frame(matrix(ncol=3, nrow=0))
    colnames(anno_t3db) <- c("target_id", "target_name", "uniprot_id")
    
    anno_toxrefdb_neoplastic <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_toxrefdb_neoplastic) <- c("smi_id", "annotation")
    
    anno_toxrefdb_nonneoplastic <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_toxrefdb_nonneoplastic) <- c("smi_id", "annotation")
    
    anno_pubchem <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_pubchem) <- c("smi_id", "annotation")
    
    anno_pubchem_props <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_pubchem_props) <- c("smi_id", "annotation")
    
    anno_gras <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_gras) <- c("smi_id", "annotation")
    
    anno_foodb_enzymes <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_foodb_enzymes) <- c("smi_id", "annotation")
    anno_foodb_flavors <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_foodb_flavors) <- c("smi_id", "annotation")
    anno_foodb_foodcontent <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_foodb_foodcontent) <- c("smi_id", "annotation")
    anno_foodb_healtheffects <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_foodb_healtheffects) <- c("smi_id", "annotation")
    anno_foodb_ontology <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_foodb_ontology) <- c("smi_id", "annotation")
    
    anno_epa_props <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_epa_props) <- c("smi_id", "annotation")

    # ADMET
    if(input[[paste0("switch__admet_binr", switch_mode)]] == TRUE){
        anno_admet_binr <- run_query(paste0("
                        SELECT DISTINCT
                            bc.dsstox_substance_id,
                            abcp.model_results
                        FROM
                            base_chemical_to_smiles bcs,
                            base_chemicals bc,
                            admet_base_chemical_predictions_interpretation_binary_v11_json abcp
                        WHERE
                            bcs.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                        AND bcs.epa_id = bc.epa_id
                        AND bcs.smi_id = abcp.smi_id
                    "), args=as.list(selected))
        anno_admet_binr <- fromJSON(anno_admet_binr$model_results[1])
        
        anno_admet_binr <- lapply(anno_admet_binr, function(x){
            if(as.numeric(x) == 0) return("in domain, inactive")
            else if(as.numeric(x) == 1) return("in domain, active")
            return("out of domain")
        })
        
        anno_admet_binr <- data.frame(model_name=names(anno_admet_binr), interpretation=unname(unlist(anno_admet_binr)))
    }

    if(input[[paste0("switch__admet_catg", switch_mode)]] == TRUE){
        anno_admet_catg_plot <- run_query(paste0("
                        SELECT DISTINCT
                            bc.dsstox_substance_id,
                            abcp.model_results
                        FROM
                            base_chemical_to_smiles bcs,
                            base_chemicals bc,
                            admet_base_chemical_predictions_interpretation_catego_v11_json abcp
                        WHERE
                            bcs.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                        AND bcs.epa_id = bc.epa_id
                        AND bcs.smi_id = abcp.smi_id
                    "), args=as.list(selected))
        anno_admet_catg_plot_rows <- lapply(anno_admet_catg_plot$model_results, function(x){
            fromJSON(x)
        })
        anno_admet_catg_plot_rows <- rbindlist(anno_admet_catg_plot_rows) %>% mutate_if(is.character, as.numeric)
        anno_admet_catg_plot_rows <- cbind(anno_admet_catg_plot$dsstox_substance_id, anno_admet_catg_plot_rows)
        colnames(anno_admet_catg_plot_rows)[1] <- "dsstox_substance_id"

        #n_electr messes up scale
        #anno_admet_catg_plot_rows$N_Electr <- NULL

        #anno_admet_catg_hm <- annotation_heatmap(anno_admet_catg_plot_rows)
    }

    if(input[[paste0("switch__admet_cont", switch_mode)]] == TRUE){
        anno_admet_cont_plot <- run_query(paste0("
                        SELECT DISTINCT
                            bc.dsstox_substance_id,
                            abcp.model_results
                        FROM
                            base_chemical_to_smiles bcs,
                            base_chemicals bc,
                            admet_base_chemical_predictions_interpretation_contin_v11_json abcp
                        WHERE
                            bcs.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                        AND bcs.epa_id = bc.epa_id
                        AND bcs.smi_id = abcp.smi_id
                    "), args=as.list(selected))
        anno_admet_cont_plot_rows <- lapply(anno_admet_cont_plot$model_results, function(x){
            fromJSON(x)
        })
        anno_admet_cont_plot_rows <- rbindlist(anno_admet_cont_plot_rows) %>% mutate_if(is.character, as.numeric)
        anno_admet_cont_plot_rows <- cbind(anno_admet_cont_plot$dsstox_substance_id, anno_admet_cont_plot_rows)
        colnames(anno_admet_cont_plot_rows)[1] <- "dsstox_substance_id"
        anno_admet_cont_plot_rows <- anno_admet_cont_plot_rows %>% select_if(~ !any(is.na(.)))
        #anno_admet_cont_hm <- annotation_heatmap(anno_admet_cont_plot_rows)
    }
    
    # ChEMBL
    if(input[[paste0("switch__chembl", switch_mode)]] == TRUE){
        anno_chembl <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            csa.structural_alerts_results
                        FROM
                            base_chemicals bc,
                            base_chemical_compounds bcc,
                            base_chemical_to_smiles bcs,
                            base_chemical_chembl_structural_alerts_hstore csa
                        WHERE
                            bc.epa_id = bcc.epa_id
                        AND bc.epa_id = bcs.epa_id
                        AND bcs.smi_id = csa.smi_id
                        AND bc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                    "), args=as.list(selected))
        chembl_alert_names <- apply(anno_chembl, 1, function(i){
            x <- i["structural_alerts_results"]
            tmp <- gsub("\"", "", unlist(str_split(x, ", ")))
            tmp2 <- lapply(tmp, function(y) unlist(str_split(y, "=>")))
            tmp3 <- lapply(tmp2, function(y) y[2])
            names(tmp3) <- lapply(tmp2, function(y) y[1])
            return(
                as.list(c(unlist(i[c("preferred_name", "casrn", "dsstox_substance_id")]), unlist(tmp3)))
            )
        })
        anno_chembl <- rbindlist(chembl_alert_names)
        
        colnames(anno_chembl) <- unlist(lapply(colnames(anno_chembl), function(x){
            if(x != "preferred_name" & x != "casrn" & x != "dsstox_substance_id"){
                return(base_chembl_alerts[as.character(base_chembl_alerts$chemblsa_id) == x, "alert_name"])
            } 
                return(x)
        }))
        lapply(colnames(anno_chembl), function(xx){
            tmpp <- anno_chembl[, ..xx]
        })
    }
    
    # CPD
    if(input[[paste0("switch__cpd", switch_mode)]] == TRUE){
        anno_cpd <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            cic.gen_cat,
                            cic.prod_fam,
                            cic.prod_type
                        FROM
                            cpd_ice_clean cic,
                            base_chemicals bc,
                            base_chemical_compounds bcc
                        WHERE
                            bc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                        AND bc.dsstox_substance_id = cic.dtxsid
                        AND bc.epa_id = bcc.epa_id
                    "), args=as.list(selected))
    }
    
    # CTD
    if(input[[paste0("switch__ctd_bioprocess", switch_mode)]] == TRUE){
        anno_ctd_bioprocess <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            mv.go_term_name,
                            mv.corrected_pvalue
                        FROM
                            base_chemicals bc,
                            base_chemical_compounds bcc,
                            ctd_biological_process_to_base_chemical_saagar_clusters mv
                        WHERE
                            bc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                        AND bc.epa_id = bcc.epa_id
                        AND bc.epa_id = mv.epa_id
                    "), args=as.list(selected))
    }
    
    if(input[[paste0("switch__ctd_cellcomp", switch_mode)]] == TRUE){
        anno_ctd_cellcomp <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            mv.go_term_name,
                            mv.corrected_pvalue
                        FROM
                            base_chemicals bc,
                            base_chemical_compounds bcc,
                            ctd_cellular_component_to_base_chemical_saagar_clusters mv
                        WHERE
                            bc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                        AND bc.epa_id = bcc.epa_id
                        AND bc.epa_id = mv.epa_id
                    "), args=as.list(selected))
    }
    
    if(input[[paste0("switch__ctd_diseases", switch_mode)]] == TRUE){
        anno_ctd_diseases <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            ccd.disease_name,
                            ccd.disease_id,
                            ccd.direct_evidence,
                            ccd.inference_gene_symbol,
                            ccd.inference_score,
                            ccd.omimids,
                            ccd.pubmed_ids
                        FROM
                            ctd_chemicals_to_diseases ccd,
                            ctd_to_base_chemicals cbc,
                            base_chemical_compounds bcc,
                            base_chemicals bc
                        WHERE
                            cbc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                        AND cbc.epa_id = bc.epa_id
                        AND cbc.ctd_id = ccd.ctd_id
                        AND cbc.epa_id = bcc.epa_id
                    "), args=as.list(selected))
    }
    if(input[[paste0("switch__ctd_genes", switch_mode)]] == TRUE){
        anno_ctd_genes <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            ccg.gene_symbol,
                            ccg.gene_id,
                            ccg.gene_forms,
                            ccg.organism,
                            ccg.interaction,
                            ccg.interaction_actions,
                            ccg.pubmed_ids
                        FROM
                            ctd_chemicals_to_genes ccg,
                            ctd_to_base_chemicals cbc,
                            base_chemical_compounds bcc,
                            base_chemicals bc
                        WHERE
                            cbc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                        AND cbc.epa_id = bc.epa_id
                        AND cbc.ctd_id = ccg.ctd_id
                        AND cbc.epa_id = bcc.epa_id
                    "), args=as.list(selected))
    }
    
    if(input[[paste0("switch__ctd_molfunct", switch_mode)]] == TRUE){
        anno_ctd_molfunct <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            mv.go_term_name,
                            mv.corrected_pvalue
                        FROM
                            base_chemicals bc,
                            base_chemical_compounds bcc,
                            ctd_molecular_function_to_base_chemical_saagar_clusters mv
                        WHERE
                            bc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                        AND bc.epa_id = bcc.epa_id
                        AND bc.epa_id = mv.epa_id
                    "), args=as.list(selected))
    }
    
    if(input[[paste0("switch__ctd_phenotypes", switch_mode)]] == TRUE){
        anno_ctd_phenotypes <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            ccp.phenotype_name,
                            ccp.phenotype_id,
                            ccp.comentioned_terms,
                            ccp.organism,
                            ccp.interaction,
                            ccp.interaction_actions,
                            ccp.anatomy_terms,
                            ccp.inference_gene_symbols,
                            ccp.pubmed_ids
                        FROM
                            ctd_chemicals_to_phenotypes ccp,
                            ctd_to_base_chemicals cbc,
                            base_chemical_compounds bcc,
                            base_chemicals bc
                        WHERE
                            cbc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                        AND cbc.epa_id = bc.epa_id
                        AND cbc.ctd_id = ccp.ctd_id
                        AND cbc.epa_id = bcc.epa_id
                    "), args=as.list(selected))
    }

    # DrugBank
    if(input[[paste0("switch__drugbank_atccodes", switch_mode)]] == TRUE){
        anno_drugbank_atccodes <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            dac.atc_code,
                            dac.atc_annotation,
                            dac.atc_level
                        FROM
                            drugbank_chemicals_to_atccodes dac,
                            drugbank_to_base_chemicals dbc,
                            drugbank_curated_chemicals dcc,
                            base_chemical_compounds bcc,
                            base_chemicals bc
                        WHERE
                            dbc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                        AND dbc.epa_id = bc.epa_id
                        AND dbc.drugbank_id = dcc.drugbank_id
                        AND dcc.db_id = dac.db_id
                        AND dbc.epa_id = bcc.epa_id
                    "), args=as.list(selected))
    }
    
    if(input[[paste0("switch__drugbank_carriers", switch_mode)]] == TRUE){
        anno_drugbank_carriers <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            mv.annotation
                        FROM
                            base_chemicals bc,
                            base_chemical_compounds bcc,
                            base_chemical_drugbank_carriers_saagar_clusters mv
                        WHERE
                            bc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                        AND bc.epa_id = bcc.epa_id
                        AND bc.epa_id = mv.epa_id
                    "), args=as.list(selected))
    }
    
    if(input[[paste0("switch__drugbank_enzymes", switch_mode)]] == TRUE){
        anno_drugbank_enzymes <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            mv.annotation
                        FROM
                            base_chemicals bc,
                            base_chemical_compounds bcc,
                            base_chemical_drugbank_enzymes_saagar_clusters mv
                        WHERE
                            bc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                        AND bc.epa_id = bcc.epa_id
                        AND bc.epa_id = mv.epa_id
                    "), args=as.list(selected))
    }
    
    
    if(input[[paste0("switch__drugbank_targets", switch_mode)]] == TRUE){
        anno_drugbank_targets <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            mv.annotation
                        FROM
                            base_chemicals bc,
                            base_chemical_compounds bcc,
                            base_chemical_drugbank_targets_saagar_clusters mv
                        WHERE
                            bc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                        AND bc.epa_id = bcc.epa_id
                        AND bc.epa_id = mv.epa_id
                    "), args=as.list(selected))
    }
    
    
    if(input[[paste0("switch__drugbank_transporters", switch_mode)]] == TRUE){
        anno_drugbank_transporters <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            mv.annotation
                        FROM
                            base_chemicals bc,
                            base_chemical_compounds bcc,
                            base_chemical_drugbank_transporters_saagar_clusters mv
                        WHERE
                            bc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                        AND bc.epa_id = bcc.epa_id
                        AND bc.epa_id = mv.epa_id
                    "), args=as.list(selected))
    }
    
    # HMDB
    if(input[[paste0("switch__hmdb_biospecimenlocations", switch_mode)]] == TRUE){
         anno_hmdb_biospecimenlocations <- run_query(paste0("
                         SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            mv.annotation
                         FROM
                             base_chemical_hmdb_biospecimen_locations_saagar_clusters mv,
                             base_chemicals bc,
                             base_chemical_compounds bcc
                         WHERE
                            bc.epa_id = bcc.epa_id
                        AND bc.epa_id = mv.epa_id
                        AND bc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                     "), args=as.list(selected))
     }

    if(input[[paste0("switch__hmdb_cellularlocations", switch_mode)]] == TRUE){
         anno_hmdb_cellularlocations <- run_query(paste0("
                         SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            mv.annotation
                         FROM
                             base_chemical_hmdb_cellular_locations_saagar_clusters mv,
                             base_chemicals bc,
                             base_chemical_compounds bcc
                         WHERE
                            bc.epa_id = bcc.epa_id
                        AND bc.epa_id = mv.epa_id
                        AND bc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                     "), args=as.list(selected))
     }
    
    if(input[[paste0("switch__hmdb_diseases", switch_mode)]] == TRUE){
         anno_hmdb_diseases <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            mv.annotation
                         FROM
                             base_chemical_hmdb_diseases_saagar_clusters mv,
                             base_chemicals bc,
                             base_chemical_compounds bcc
                         WHERE
                            bc.epa_id = bcc.epa_id
                        AND bc.epa_id = mv.epa_id
                        AND bc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                     "), args=as.list(selected))
     }
    
    if(input[[paste0("switch__hmdb_genes", switch_mode)]] == TRUE){
         anno_hmdb_genes <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            mv.annotation
                         FROM
                             base_chemical_hmdb_genes_saagar_clusters mv,
                             base_chemicals bc,
                             base_chemical_compounds bcc
                         WHERE
                            bc.epa_id = bcc.epa_id
                        AND bc.epa_id = mv.epa_id
                        AND bc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                     "), args=as.list(selected))
     }
    
    if(input[[paste0("switch__hmdb_tissuelocations", switch_mode)]] == TRUE){
         anno_hmdb_tissuelocations <- run_query(paste0("
                         SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            mv.annotation
                         FROM
                             base_chemical_hmdb_tissue_locations_saagar_clusters mv,
                             base_chemicals bc,
                             base_chemical_compounds bcc
                         WHERE
                            bc.epa_id = bcc.epa_id
                        AND bc.epa_id = mv.epa_id
                        AND bc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                     "), args=as.list(selected))
     }
    
    # ICE
    if(input[[paste0("switch__ice", switch_mode)]] == TRUE){
        anno_ice <- query_ice(chemicals=dtxsids)
    }
    
    # InVitroDB
    if(input[[paste0("switch__invitrodb", switch_mode)]] == TRUE){
        anno_invitrodb <- run_query(paste0("
            SELECT DISTINCT
                bcc.preferred_name,
                bcc.casrn,
                bc.dsstox_substance_id,
                isr.invitrodb_values AS assay_name,
                iaea.assay_endpoint_attribute,
                iaea.assay_endpoint_value,
                x.invitrodb_values AS hit_call
            FROM
                invitrodb_summarized_results isr,
                invitrodb_annotation ia,
                invitrodb_assay_endpoint_annotation iaea,
                base_chemicals bc,
                base_chemical_to_smiles bcs,
                base_chemical_compounds bcc,
                (SELECT i.row_index, i.invitrodb_attributes, i.invitrodb_values FROM invitrodb_summarized_results i WHERE i.invitrodb_attributes='hitc') AS x
            WHERE
                bc.epa_id = bcs.epa_id
            AND bc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
            AND bc.epa_id = bcc.epa_id
            AND bcs.smi_id = isr.smi_id
            AND isr.invitrodb_attributes = ia.invitrodb_attributes
            AND (isr.invitrodb_attributes = 'aenm' AND isr.invitrodb_values = iaea.aenm)
            AND isr.row_index = x.row_index
        "), args=as.list(selected))
    }
    
    # Leadscope
    if(input[[paste0("switch__leadscope", switch_mode)]] == TRUE){
        anno_leadscope <- run_query(paste0("
                        SELECT DISTINCT
                            --bcc.preferred_name,
                            --bcc.casrn,
                            --bc.dsstox_substance_id,
                            lpm.model_name,
                            lcpi.interpretation,
                            --lpm.short_description
                            lqd.endpoint_desc AS short_description
                        FROM
                            leadscope_chemical_predictions_interpretation lcpi,
                            base_chemical_to_smiles bcs,
                            leadscope_predictive_models lpm,
                            base_chemical_compounds bcc,
                            base_chemicals bc,
                            leadscope_qmrf_descriptions lqd
                        WHERE
                            bcs.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                        AND bcs.epa_id = bc.epa_id
                        AND bcs.epa_id = bcc.epa_id
                        AND bcs.smi_id = lcpi.smi_id
                        AND lcpi.lmodel_id = lpm.lmodel_id
                        AND lcpi.lmodel_id = lqd.lmodel_id
                    "), args=as.list(selected))
        # anno_admet_binr <- fromJSON(anno_admet_binr$model_results[1])
        anno_leadscope$interpretation <- lapply(anno_leadscope$interpretation, function(x){
            if(as.numeric(x) == 0) return("in domain, inactive")
            else if(as.numeric(x) == 1) return("in domain, active")
            return("out of domain")
        })
    }
    
    # OChem
    if(input[[paste0("switch__ochem", switch_mode)]] == TRUE){
        anno_ochem <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            mv.alert_name,
                            mv.saresult
                        FROM
                            mv_ochem_sa_positives_saagar mv,
                            base_chemical_to_smiles bcs,
                            base_chemical_compounds bcc,
                            base_chemicals bc
                        WHERE
                            bcs.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                        AND bcs.epa_id = bc.epa_id
                        AND bcs.epa_id = bcc.epa_id
                        AND bcs.epa_id = mv.epa_id
                    "), args=as.list(selected))
    }
    
    # Opera
    if(input[[paste0("switch__opera_adme", switch_mode)]] == TRUE){
        anno_opera_adme <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            opm.model_name,
                            ocp.model_value,
                            opm.category,
                            opm.model_description
                        FROM
                            base_chemicals bc,
                            base_chemical_to_smiles bcs,
                            base_chemical_compounds bcc,
                            opera_predictive_models opm,
                            opera_chemical_predictions ocp
                        WHERE
                            bc.epa_id = bcs.epa_id
                        AND bcs.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                        AND bc.epa_id = bcc.epa_id
                        AND bcs.smi_id = ocp.smi_id
                        AND ocp.omodel_id = opm.omodel_id
                        AND opm.category = 'ADME endpoints'
                    "), args=as.list(selected))
    }

    if(input[[paste0("switch__opera_environmental_fate", switch_mode)]] == TRUE){
        anno_opera_environmental_fate <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            opm.model_name,
                            ocp.model_value,
                            opm.category,
                            opm.model_description
                        FROM
                            base_chemicals bc,
                            base_chemical_to_smiles bcs,
                            base_chemical_compounds bcc,
                            opera_predictive_models opm,
                            opera_chemical_predictions ocp
                        WHERE
                            bc.epa_id = bcs.epa_id
                        AND bcs.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                        AND bc.epa_id = bcc.epa_id
                        AND bcs.smi_id = ocp.smi_id
                        AND ocp.omodel_id = opm.omodel_id
                        AND opm.category = 'Environmental fate'
                    "), args=as.list(selected))
    }

    if(input[[paste0("switch__opera_physicochem", switch_mode)]] == TRUE){
        anno_opera_physicochem <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            opm.model_name,
                            ocp.model_value,
                            opm.category,
                            opm.model_description
                        FROM
                            base_chemicals bc,
                            base_chemical_to_smiles bcs,
                            base_chemical_compounds bcc,
                            opera_predictive_models opm,
                            opera_chemical_predictions ocp
                        WHERE
                            bc.epa_id = bcs.epa_id
                        AND bcs.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                        AND bc.epa_id = bcc.epa_id
                        AND bcs.smi_id = ocp.smi_id
                        AND ocp.omodel_id = opm.omodel_id
                        AND opm.category = 'Physicochemical properties'
                    "), args=as.list(selected))
    }

    if(input[[paste0("switch__opera_structural_properties", switch_mode)]] == TRUE){
        anno_opera_structural_properties <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            opm.model_name,
                            ocp.model_value,
                            opm.category,
                            opm.model_description
                        FROM
                            base_chemicals bc,
                            base_chemical_to_smiles bcs,
                            base_chemical_compounds bcc,
                            opera_predictive_models opm,
                            opera_chemical_predictions ocp
                        WHERE
                            bc.epa_id = bcs.epa_id
                        AND bcs.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                        AND bc.epa_id = bcc.epa_id
                        AND bcs.smi_id = ocp.smi_id
                        AND ocp.omodel_id = opm.omodel_id
                        AND opm.category = 'Structural properties'
                    "), args=as.list(selected))
    }

    if(input[[paste0("switch__opera_toxicology_endpoints", switch_mode)]] == TRUE){
        anno_opera_toxicology_endpoints <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            opm.model_name,
                            ocp.model_value,
                            opm.category,
                            opm.model_description
                        FROM
                            base_chemicals bc,
                            base_chemical_to_smiles bcs,
                            base_chemical_compounds bcc,
                            opera_predictive_models opm,
                            opera_chemical_predictions ocp
                        WHERE
                            bc.epa_id = bcs.epa_id
                        AND bcs.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                        AND bc.epa_id = bcc.epa_id
                        AND bcs.smi_id = ocp.smi_id
                        AND ocp.omodel_id = opm.omodel_id
                        AND opm.category = 'Toxicology endpoints'
                    "), args=as.list(selected))
    }
    
    # Saagar
    if(input[[paste0("switch__saagar", switch_mode)]] == TRUE){
        anno_saagar <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            mv.alert_name,
                            mv.saresult
                        FROM
                            mv_saagar_de_positives_saagar mv,
                            base_chemical_to_smiles bcs,
                            base_chemical_compounds bcc,
                            base_chemicals bc
                        WHERE
                            bcs.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                        AND bcs.epa_id = bc.epa_id
                        AND bcs.epa_id = bcc.epa_id
                        AND bcs.epa_id = mv.epa_id
                    "), args=as.list(selected))
    }
    
    # T3DB
    if(input[[paste0("switch__t3db", switch_mode)]] == TRUE){
        anno_t3db <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            tct.target_id,
                            tct.target_name,
                            tct.uniprot_id
                        FROM
                            t3db_chemicals_to_targets tct,
                            t3db_to_base_chemicals tbc,
                            t3db_curated_chemicals tcc,
                            base_chemical_compounds bcc,
                            base_chemicals bc
                        WHERE
                            tbc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                        AND bc.epa_id = bcc.epa_id
                        AND tbc.epa_id = bc.epa_id
                        AND tbc.t3db_id = tcc.t3db_id
                        AND tcc.t3_id = tct.t3_id
                    "), args=as.list(selected))
    }
    
    # ToxRefDB
    if(input[[paste0("switch__toxrefdb_neoplastic", switch_mode)]] == TRUE){
        anno_toxrefdb_neoplastic <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            CONCAT(mv.disease_type, '_',
                            mv.study_type, '_',
                            mv.species, '_',
                            mv.sex, '_',
                            mv.life_stage, '_',
                            mv.endpoint_type, '_',
                            mv.endpoint_target) AS annotation
                        FROM
                            mv_toxrefdb_cancer_saagar mv,
                            base_chemical_compounds bcc,
                            base_chemicals bc
                        WHERE
                            mv.epa_id = bcc.epa_id
                        AND mv.epa_id = bc.epa_id
                        AND mv.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                    "), args=as.list(selected))
    }
    
    if(input[[paste0("switch__toxrefdb_nonneoplastic", switch_mode)]] == TRUE){
        anno_toxrefdb_nonneoplastic <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            CONCAT(mv.disease_type, '_',
                            mv.study_type, '_',
                            mv.species, '_',
                            mv.sex, '_',
                            mv.life_stage, '_',
                            mv.endpoint_type, '_',
                            mv.endpoint_target) AS annotation
                        FROM
                            mv_toxrefdb_nonneoplastic_saagar mv,
                            base_chemical_compounds bcc,
                            base_chemicals bc
                        WHERE
                            mv.epa_id = bcc.epa_id
                        AND mv.epa_id = bc.epa_id
                        AND mv.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                    "), args=as.list(selected))
    }

    if(input[[paste0("switch__pubchem", switch_mode)]] == TRUE){
        anno_pubchem <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            pb.bioassay_name,
                            pb.source_name,
                            pb.source_id,
                            pb.substance_type,
                            pb.outcome_type,
                            pb.project_category,
                            pb.bioassay_group,
                            pb.bioassay_types,
                            pb.protein_accessions,
                            pbm.aid,
                            pbm.sid,
                            pbm.sid_group,
                            pbm.activity_outcome,
                            pbm.activity_name,
                            pbm.activity_qualifier,
                            pbm.activity_value,
                            pbm.protein_accession,
                            pbm.gene_id,
                            pbm.target_tax_id,
                            pbm.pmid
                        FROM
                            base_chemical_to_pubchem_cid bcpc,
                            pubchem_bioactivity_mapping pbm,
                            pubchem_bioassays pb,
                            base_chemicals bc,
                            base_chemical_compounds bcc
                        WHERE
                            bcpc.pubchem_cid = pbm.cid
                        AND bcpc.epa_id = bc.epa_id
                        AND bcpc.epa_id = bcc.epa_id
                        AND pbm.aid = pb.aid
                        AND bcpc.epa_id IN (
                            ", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                    "), args=as.list(selected))
    }
    
    if(input[[paste0("switch__pubchem_props", switch_mode)]] == TRUE){
        anno_pubchem_props <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            pba.pubchem_attribute,
                            pba.pubchem_value
                        FROM
                            base_chemical_to_pubchem_cid bcpc,
                            pubchem_attributes pba,
                            base_chemicals bc,
                            base_chemical_compounds bcc
                        WHERE
                            bcpc.pubchem_cid = pba.cid
                        AND bcpc.epa_id = bc.epa_id
                        AND bcpc.epa_id = bcc.epa_id
                        AND bcpc.epa_id IN (
                            ", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                    "), args=as.list(selected))
    }
    
    if(input[[paste0("switch__gras", switch_mode)]] == TRUE){
        anno_gras <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            bcg.gras_substance,
                            bcg.scogs_report_number,
                            bcg.year,
                            bcg.scogs_conclusion,
                            CONCAT(gsc.scogs_conclusion, ': ', gsc.scogs_interpretation) AS scogs_interpretation,
                            bcg.\"21_cfr_regulation\",
                            bcg.ntis_accession
                        FROM
                            base_chemical_to_gras bcg,
                            base_chemicals bc,
                            base_chemical_compounds bcc,
                            gras_scogs_conclusions gsc
                        WHERE
                            bcg.epa_id = bc.epa_id
                        AND bcg.scogs_conclusion LIKE '%' || gsc.scogs_conclusion || '%'
                        AND bcg.epa_id = bcc.epa_id
                        AND bcg.epa_id IN (
                            ", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                    "), args=as.list(selected))
    }

    if(input[[paste0("switch__foodb_enzymes", switch_mode)]] == TRUE){
        anno_foodb_enzymes <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            fce.enzyme_name,
                            fce.gene_name,
                            fce.kingdom,
                            fce.superclass,
                            fce.class,
                            fce.subclass,
                            fce.citations
                        FROM
                            base_chemical_to_foodb_compound bcfc,
                            foodb_compound_enzymes fce,
                            base_chemicals bc,
                            base_chemical_compounds bcc
                        WHERE
                            bcfc.epa_id = bc.epa_id
                        AND bcfc.foodb_compound_id = fce.compound_id
                        AND bcfc.epa_id = bcc.epa_id
                        AND bcfc.epa_id IN (
                            ", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                    "), args=as.list(selected))
    }
    
    if(input[[paste0("switch__foodb_flavors", switch_mode)]] == TRUE){
        anno_foodb_flavors <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            fcf.flavor_name,
                            fcf.flavor_group,
                            fcf.category,
                            fcf.kingdom,
                            fcf.superclass,
                            fcf.class,
                            fcf.subclass,
                            fcf.citations
                        
                        FROM
                            base_chemical_to_foodb_compound bcfc,
                            foodb_compound_flavors fcf,
                            base_chemicals bc,
                            base_chemical_compounds bcc
                        WHERE
                            bcfc.epa_id = bc.epa_id
                        AND bcfc.foodb_compound_id = fcf.compound_id
                        AND bcfc.epa_id = bcc.epa_id
                        AND bcfc.epa_id IN (
                            ", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                    "), args=as.list(selected))
    }
    
    if(input[[paste0("switch__foodb_foodcontent", switch_mode)]] == TRUE){
        anno_foodb_foodcontent <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            fcfc.name,
                            fcfc.name_scientific,
                            fcfc.food_group,
                            fcfc.food_subgroup,
                            fcfc.food_type,
                            fcfc.description,
                            fcfc.kingdom,
                            fcfc.superclass,
                            fcfc.class,
                            fcfc.subclass,
                            fcfc.orig_content,
                            fcfc.orig_unit
                        FROM
                            base_chemical_to_foodb_compound bcfc,
                            foodb_compound_food_content fcfc,
                            base_chemicals bc,
                            base_chemical_compounds bcc
                        WHERE
                            bcfc.epa_id = bc.epa_id
                        AND bcfc.foodb_compound_id = fcfc.compound_id
                        AND bcfc.epa_id = bcc.epa_id
                        AND fcfc.orig_content > 0
                        AND bcfc.epa_id IN (
                            ", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                    "), args=as.list(selected))
    }
    
    if(input[[paste0("switch__foodb_healtheffects", switch_mode)]] == TRUE){
        anno_foodb_healtheffects <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            fche.health_effect_name,
                            fche.kingdom,
                            fche.superclass,
                            fche.class,
                            fche.subclass,
                            fche.description
                        FROM
                            base_chemical_to_foodb_compound bcfc,
                            foodb_compound_health_effects fche,
                            base_chemicals bc,
                            base_chemical_compounds bcc
                        WHERE
                            bcfc.epa_id = bc.epa_id
                        AND bcfc.foodb_compound_id = fche.compound_id
                        AND bcfc.epa_id = bcc.epa_id
                        AND bcfc.epa_id IN (
                            ", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                    "), args=as.list(selected))
    }
    
    
    if(input[[paste0("switch__foodb_ontology", switch_mode)]] == TRUE){
        anno_foodb_ontology <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            fcot.term,
                            fcot.definition,
                            fcot.level,
                            fcot.kingdom,
                            fcot.superclass,
                            fcot.class,
                            fcot.subclass
                        FROM
                            base_chemical_to_foodb_compound bcfc,
                            foodb_compound_ontology_terms fcot,
                            base_chemicals bc,
                            base_chemical_compounds bcc
                        WHERE
                            bcfc.epa_id = bc.epa_id
                        AND bcfc.foodb_compound_id = fcot.compound_id
                        AND bcfc.epa_id = bcc.epa_id
                        AND bcfc.epa_id IN (
                            ", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                    "), args=as.list(selected))
    }
    
    # EPA
    if(input[[paste0("switch__epa_props", switch_mode)]] == TRUE){
        anno_epa_props <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            ep.epa_attribute,
                            ep.epa_attribute_value
                        FROM
                            epa_properties ep,
                            base_chemicals bc,
                            base_chemical_compounds bcc
                        WHERE
                            ep.epa_id = bc.epa_id
                        AND ep.epa_id = bcc.epa_id
                        AND ep.epa_id IN (
                            ", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                    "), args=as.list(selected))
    }
    return(list(
        anno_admet_binr=anno_admet_binr,
        anno_admet_catg=anno_admet_catg,
        anno_admet_cont=anno_admet_cont,
        
        anno_chembl=anno_chembl,
        anno_cpd=anno_cpd,
        
        anno_ctd_bioprocess=anno_ctd_bioprocess,
        anno_ctd_cellcomp=anno_ctd_cellcomp,
        anno_ctd_diseases=anno_ctd_diseases,
        anno_ctd_genes=anno_ctd_genes,
        anno_ctd_molfunct=anno_ctd_molfunct,
        anno_ctd_phenotypes=anno_ctd_phenotypes,

        anno_drugbank_atccodes=anno_drugbank_atccodes,
        anno_drugbank_carriers=anno_drugbank_carriers,
        anno_drugbank_enzymes=anno_drugbank_enzymes,
        anno_drugbank_targets=anno_drugbank_targets,
        anno_drugbank_transporters=anno_drugbank_transporters,

        anno_hmdb_biospecimenlocations=anno_hmdb_biospecimenlocations,
        anno_hmdb_cellularlocations=anno_hmdb_cellularlocations,
        anno_hmdb_diseases=anno_hmdb_diseases,
        anno_hmdb_genes=anno_hmdb_genes,
        anno_hmdb_tissuelocations=anno_hmdb_tissuelocations,

        anno_ice=anno_ice,

        anno_invitrodb=anno_invitrodb,
        anno_leadscope=anno_leadscope,
        anno_ochem=anno_ochem,

        anno_opera_adme=anno_opera_adme,
        anno_opera_environmental_fate=anno_opera_environmental_fate,
        anno_opera_physicochem=anno_opera_physicochem,
        anno_opera_structural_properties=anno_opera_structural_properties,
        anno_opera_toxicology_endpoints=anno_opera_toxicology_endpoints,

        anno_saagar=anno_saagar,
        anno_t3db=anno_t3db,
        anno_toxrefdb_neoplastic=anno_toxrefdb_neoplastic,
        anno_toxrefdb_nonneoplastic=anno_toxrefdb_nonneoplastic,
        anno_pubchem=anno_pubchem,
        anno_pubchem_props=anno_pubchem_props,
        anno_gras=anno_gras,

        anno_foodb_enzymes=anno_foodb_enzymes,
        anno_foodb_flavors=anno_foodb_flavors,
        anno_foodb_foodcontent=anno_foodb_foodcontent,
        anno_foodb_healtheffects=anno_foodb_healtheffects,
        anno_foodb_ontology=anno_foodb_ontology,
        anno_epa_props=anno_epa_props
    ))
}








# Load chemical annotations for a chemical for use with LLMs
llm_annotations <- function(selected, dtxsids){
    
    admet_desc <- run_query(paste0("
                        SELECT
                            model_name,
                            description
                        FROM
                            admet_predictive_models
                    "))
    
    
    anno_admet_binr <- data.frame(matrix(ncol=3, nrow=0))
    anno_admet_catg <- data.frame(matrix(ncol=3, nrow=0))
    anno_admet_cont <- data.frame(matrix(ncol=3, nrow=0))
    anno_chembl <- data.frame(matrix(ncol=2, nrow=0))
    anno_cpd <- data.frame(matrix(ncol=3, nrow=0))
    anno_ctd_bioprocess <- data.frame(matrix(ncol=2, nrow=0))
    anno_ctd_cellcomp <- data.frame(matrix(ncol=2, nrow=0))
    anno_ctd_diseases <- data.frame(matrix(ncol=7, nrow=0))
    anno_ctd_genes <- data.frame(matrix(ncol=7, nrow=0))
    anno_ctd_molfunct <- data.frame(matrix(ncol=2, nrow=0))
    anno_ctd_phenotypes <- data.frame(matrix(ncol=9, nrow=0))
    anno_drugbank_atccodes <- data.frame(matrix(ncol=3, nrow=0))
    anno_drugbank_carriers <- data.frame(matrix(ncol=2, nrow=0))
    anno_drugbank_enzymes <- data.frame(matrix(ncol=2, nrow=0))
    anno_drugbank_targets <- data.frame(matrix(ncol=2, nrow=0))
    anno_drugbank_transporters <- data.frame(matrix(ncol=2, nrow=0))
    anno_hmdb_biospecimenlocations <- data.frame(matrix(ncol=2, nrow=0))
    anno_hmdb_cellularlocations <- data.frame(matrix(ncol=2, nrow=0))
    anno_hmdb_diseases <- data.frame(matrix(ncol=2, nrow=0))
    anno_hmdb_genes <- data.frame(matrix(ncol=2, nrow=0))
    anno_hmdb_tissuelocations <- data.frame(matrix(ncol=2, nrow=0))
    anno_ice <- data.frame(matrix(ncol=2, nrow=0))
    anno_invitrodb <- data.frame(matrix(ncol=2, nrow=0))
    anno_leadscope <- data.frame(matrix(ncol=3, nrow=0))
    anno_ochem <- data.frame(matrix(ncol=2, nrow=0))
    anno_opera_adme <- data.frame(matrix(ncol=6, nrow=0))
    anno_opera_environmental_fate <- data.frame(matrix(ncol=6, nrow=0))
    anno_opera_physicochem <- data.frame(matrix(ncol=6, nrow=0))
    anno_opera_structural_properties <- data.frame(matrix(ncol=6, nrow=0))
    anno_opera_toxicology_endpoints <- data.frame(matrix(ncol=6, nrow=0))
    anno_saagar <- data.frame(matrix(ncol=2, nrow=0))
    anno_t3db <- data.frame(matrix(ncol=3, nrow=0))
    anno_toxrefdb_neoplastic <- data.frame(matrix(ncol=2, nrow=0))
    anno_toxrefdb_nonneoplastic <- data.frame(matrix(ncol=2, nrow=0))
    anno_pubchem <- data.frame(matrix(ncol=2, nrow=0))
    anno_pubchem_props <- data.frame(matrix(ncol=2, nrow=0))
    anno_gras <- data.frame(matrix(ncol=2, nrow=0))
    anno_foodb_enzymes <- data.frame(matrix(ncol=2, nrow=0))
    anno_foodb_flavors <- data.frame(matrix(ncol=2, nrow=0))
    anno_foodb_foodcontent <- data.frame(matrix(ncol=2, nrow=0))
    anno_foodb_healtheffects <- data.frame(matrix(ncol=2, nrow=0))
    anno_foodb_ontology <- data.frame(matrix(ncol=2, nrow=0))
    anno_epa_props <- data.frame(matrix(ncol=2, nrow=0))

        # anno_admet_binr_plot <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bc.dsstox_substance_id,
        #                     abcp.model_results
        #                 FROM
        #                     base_chemical_to_smiles bcs,
        #                     base_chemicals bc,
        #                     admet_base_chemical_predictions_interpretation_binary_v11_json abcp
        #                 WHERE
        #                     bcs.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #                 AND bcs.epa_id = bc.epa_id
        #                 AND bcs.smi_id = abcp.smi_id
        #             "), args=as.list(selected))
        # anno_admet_binr_plot_rows <- lapply(anno_admet_binr_plot$model_results, function(x){
        #     fromJSON(x)
        # })
        # anno_admet_binr_plot_rows <- rbindlist(anno_admet_binr_plot_rows) %>% mutate_if(is.character, as.numeric)
        # anno_admet_binr_plot_cols <- colnames(anno_admet_binr_plot_rows)
        # anno_admet_binr_plot_rows <- transpose(anno_admet_binr_plot_rows)
        # anno_admet_binr_plot <- cbind(anno_admet_binr_plot_cols, anno_admet_binr_plot_rows)
        # colnames(anno_admet_binr_plot) <- c("model", "value")
        # anno_admet_binr_plot$description <- apply(anno_admet_binr_plot, 1, function(x){
        #     admet_desc[admet_desc$model_name == x, "description"]
        # })
        # anno_admet_binr <- anno_admet_binr_plot
        # 
        # 
        # 
        # 
        # 
        # anno_admet_catg_plot <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bc.dsstox_substance_id,
        #                     abcp.model_results
        #                 FROM
        #                     base_chemical_to_smiles bcs,
        #                     base_chemicals bc,
        #                     admet_base_chemical_predictions_interpretation_catego_v11_json abcp
        #                 WHERE
        #                     bcs.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #                 AND bcs.epa_id = bc.epa_id
        #                 AND bcs.smi_id = abcp.smi_id
        #             "), args=as.list(selected))
        # anno_admet_catg_plot_rows <- lapply(anno_admet_catg_plot$model_results, function(x){
        #     fromJSON(x)
        # })
        # anno_admet_catg_plot_rows <- rbindlist(anno_admet_catg_plot_rows) %>% mutate_if(is.character, as.numeric)
        # anno_admet_catg_plot_cols <- colnames(anno_admet_catg_plot_rows)
        # anno_admet_catg_plot_rows <- transpose(anno_admet_catg_plot_rows)
        # anno_admet_catg_plot <- cbind(anno_admet_catg_plot_cols, anno_admet_catg_plot_rows)
        # colnames(anno_admet_catg_plot) <- c("model", "value")
        # anno_admet_catg_plot$description <- apply(anno_admet_catg_plot, 1, function(x){
        #     admet_desc[admet_desc$model_name == x, "description"]
        # })
        # anno_admet_catg <- anno_admet_catg_plot
        # 
        # anno_admet_cont_plot <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bc.dsstox_substance_id,
        #                     abcp.model_results
        #                 FROM
        #                     base_chemical_to_smiles bcs,
        #                     base_chemicals bc,
        #                     admet_base_chemical_predictions_interpretation_contin_v11_json abcp
        #                 WHERE
        #                     bcs.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #                 AND bcs.epa_id = bc.epa_id
        #                 AND bcs.smi_id = abcp.smi_id
        #             "), args=as.list(selected))
        # anno_admet_cont_plot_rows <- lapply(anno_admet_cont_plot$model_results, function(x){
        #     fromJSON(x)
        # })
        # anno_admet_cont_plot_rows <- rbindlist(anno_admet_cont_plot_rows) %>% mutate_if(is.character, as.numeric)
        # anno_admet_cont_plot_cols <- colnames(anno_admet_cont_plot_rows)
        # anno_admet_cont_plot_rows <- transpose(anno_admet_cont_plot_rows)
        # anno_admet_cont_plot <- cbind(anno_admet_cont_plot_cols, anno_admet_cont_plot_rows)
        # colnames(anno_admet_cont_plot) <- c("model", "value")
        # anno_admet_cont_plot$description <- apply(anno_admet_cont_plot, 1, function(x){
        #     admet_desc[admet_desc$model_name == x, "description"]
        # })
        # anno_admet_cont <- anno_admet_cont_plot
        # 
        # anno_chembl <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     csa.structural_alerts_results
        #                 FROM
        #                     base_chemicals bc,
        #                     base_chemical_compounds bcc,
        #                     base_chemical_to_smiles bcs,
        #                     base_chemical_chembl_structural_alerts_hstore csa
        #                 WHERE
        #                     bc.epa_id = bcc.epa_id
        #                 AND bc.epa_id = bcs.epa_id
        #                 AND bcs.smi_id = csa.smi_id
        #                 AND bc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #             "), args=as.list(selected))
        # chembl_alert_names <- apply(anno_chembl, 1, function(i){
        #     x <- i["structural_alerts_results"]
        #     tmp <- gsub("\"", "", unlist(str_split(x, ", ")))
        #     tmp2 <- lapply(tmp, function(y) unlist(str_split(y, "=>")))
        #     tmp3 <- lapply(tmp2, function(y) y[2])
        #     names(tmp3) <- lapply(tmp2, function(y) y[1])
        #     return(
        #         #as.list(c(unlist(i[c("preferred_name", "casrn", "dsstox_substance_id")]), unlist(tmp3)))
        #         as.list(unlist(tmp3))
        #     )
        # })
        # anno_chembl <- rbindlist(chembl_alert_names)
        # 
        # anno_chembl_cols <- unlist(lapply(colnames(anno_chembl), function(x){
        #     if(x != "preferred_name" & x != "casrn" & x != "dsstox_substance_id"){
        #         return(base_chembl_alerts[as.character(base_chembl_alerts$chemblsa_id) == x, "alert_name"])
        #     } 
        #     return(x)
        # }))
        # anno_chembl <- transpose(anno_chembl)
        # anno_chembl <- cbind(anno_chembl_cols, anno_chembl)
        # colnames(anno_chembl) <- c("feature", "present")
        # 
        # anno_cpd <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     cic.gen_cat,
        #                     cic.prod_fam,
        #                     cic.prod_type
        #                 FROM
        #                     cpd_ice_clean cic,
        #                     base_chemicals bc,
        #                     base_chemical_compounds bcc
        #                 WHERE
        #                     bc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #                 AND bc.dsstox_substance_id = cic.dtxsid
        #                 AND bc.epa_id = bcc.epa_id
        #             "), args=as.list(selected))
        # 
        # anno_ctd_bioprocess <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     mv.go_term_name,
        #                     mv.corrected_pvalue
        #                 FROM
        #                     base_chemicals bc,
        #                     base_chemical_compounds bcc,
        #                     ctd_biological_process_to_base_chemical_saagar_clusters mv
        #                 WHERE
        #                     bc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #                 AND bc.epa_id = bcc.epa_id
        #                 AND bc.epa_id = mv.epa_id
        #             "), args=as.list(selected))
        # 
        # anno_ctd_cellcomp <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     mv.go_term_name,
        #                     mv.corrected_pvalue
        #                 FROM
        #                     base_chemicals bc,
        #                     base_chemical_compounds bcc,
        #                     ctd_cellular_component_to_base_chemical_saagar_clusters mv
        #                 WHERE
        #                     bc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #                 AND bc.epa_id = bcc.epa_id
        #                 AND bc.epa_id = mv.epa_id
        #             "), args=as.list(selected))
        # 
        # anno_ctd_diseases <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     ccd.disease_name,
        #                     ccd.disease_id,
        #                     ccd.direct_evidence,
        #                     ccd.inference_gene_symbol,
        #                     ccd.inference_score,
        #                     ccd.omimids,
        #                     ccd.pubmed_ids
        #                 FROM
        #                     ctd_chemicals_to_diseases ccd,
        #                     ctd_to_base_chemicals cbc,
        #                     base_chemical_compounds bcc,
        #                     base_chemicals bc
        #                 WHERE
        #                     cbc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #                 AND cbc.epa_id = bc.epa_id
        #                 AND cbc.ctd_id = ccd.ctd_id
        #                 AND cbc.epa_id = bcc.epa_id
        #             "), args=as.list(selected))
        # 
        # anno_ctd_genes <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     ccg.gene_symbol,
        #                     ccg.gene_id,
        #                     ccg.gene_forms,
        #                     ccg.organism,
        #                     ccg.interaction,
        #                     ccg.interaction_actions,
        #                     ccg.pubmed_ids
        #                 FROM
        #                     ctd_chemicals_to_genes ccg,
        #                     ctd_to_base_chemicals cbc,
        #                     base_chemical_compounds bcc,
        #                     base_chemicals bc
        #                 WHERE
        #                     cbc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #                 AND cbc.epa_id = bc.epa_id
        #                 AND cbc.ctd_id = ccg.ctd_id
        #                 AND cbc.epa_id = bcc.epa_id
        #             "), args=as.list(selected))
        # 
        # anno_ctd_molfunct <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     mv.go_term_name,
        #                     mv.corrected_pvalue
        #                 FROM
        #                     base_chemicals bc,
        #                     base_chemical_compounds bcc,
        #                     ctd_molecular_function_to_base_chemical_saagar_clusters mv
        #                 WHERE
        #                     bc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #                 AND bc.epa_id = bcc.epa_id
        #                 AND bc.epa_id = mv.epa_id
        #             "), args=as.list(selected))
        # 
        # anno_ctd_phenotypes <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     ccp.phenotype_name,
        #                     ccp.phenotype_id,
        #                     ccp.comentioned_terms,
        #                     ccp.organism,
        #                     ccp.interaction,
        #                     ccp.interaction_actions,
        #                     ccp.anatomy_terms,
        #                     ccp.inference_gene_symbols,
        #                     ccp.pubmed_ids
        #                 FROM
        #                     ctd_chemicals_to_phenotypes ccp,
        #                     ctd_to_base_chemicals cbc,
        #                     base_chemical_compounds bcc,
        #                     base_chemicals bc
        #                 WHERE
        #                     cbc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #                 AND cbc.epa_id = bc.epa_id
        #                 AND cbc.ctd_id = ccp.ctd_id
        #                 AND cbc.epa_id = bcc.epa_id
        #             "), args=as.list(selected))
        # 
        # anno_drugbank_atccodes <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     dac.atc_code,
        #                     dac.atc_annotation,
        #                     dac.atc_level
        #                 FROM
        #                     drugbank_chemicals_to_atccodes dac,
        #                     drugbank_to_base_chemicals dbc,
        #                     drugbank_curated_chemicals dcc,
        #                     base_chemical_compounds bcc,
        #                     base_chemicals bc
        #                 WHERE
        #                     dbc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #                 AND dbc.epa_id = bc.epa_id
        #                 AND dbc.drugbank_id = dcc.drugbank_id
        #                 AND dcc.db_id = dac.db_id
        #                 AND dbc.epa_id = bcc.epa_id
        #             "), args=as.list(selected))
        # 
        # anno_drugbank_carriers <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     mv.annotation
        #                 FROM
        #                     base_chemicals bc,
        #                     base_chemical_compounds bcc,
        #                     base_chemical_drugbank_carriers_saagar_clusters mv
        #                 WHERE
        #                     bc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #                 AND bc.epa_id = bcc.epa_id
        #                 AND bc.epa_id = mv.epa_id
        #             "), args=as.list(selected))
        # 
        # anno_drugbank_enzymes <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     mv.annotation
        #                 FROM
        #                     base_chemicals bc,
        #                     base_chemical_compounds bcc,
        #                     base_chemical_drugbank_enzymes_saagar_clusters mv
        #                 WHERE
        #                     bc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #                 AND bc.epa_id = bcc.epa_id
        #                 AND bc.epa_id = mv.epa_id
        #             "), args=as.list(selected))
        # 
        # anno_drugbank_targets <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     mv.annotation
        #                 FROM
        #                     base_chemicals bc,
        #                     base_chemical_compounds bcc,
        #                     base_chemical_drugbank_targets_saagar_clusters mv
        #                 WHERE
        #                     bc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #                 AND bc.epa_id = bcc.epa_id
        #                 AND bc.epa_id = mv.epa_id
        #             "), args=as.list(selected))
        # 
        # anno_drugbank_transporters <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     mv.annotation
        #                 FROM
        #                     base_chemicals bc,
        #                     base_chemical_compounds bcc,
        #                     base_chemical_drugbank_transporters_saagar_clusters mv
        #                 WHERE
        #                     bc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #                 AND bc.epa_id = bcc.epa_id
        #                 AND bc.epa_id = mv.epa_id
        #             "), args=as.list(selected))
        # 
        # anno_hmdb_biospecimenlocations <- run_query(TRUE, paste0("
        #                  SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     mv.annotation
        #                  FROM
        #                      base_chemical_hmdb_biospecimen_locations_saagar_clusters mv,
        #                      base_chemicals bc,
        #                      base_chemical_compounds bcc
        #                  WHERE
        #                     bc.epa_id = bcc.epa_id
        #                 AND bc.epa_id = mv.epa_id
        #                 AND bc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #              "), args=as.list(selected))
        # 
        # anno_hmdb_cellularlocations <- run_query(TRUE, paste0("
        #                  SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     mv.annotation
        #                  FROM
        #                      base_chemical_hmdb_cellular_locations_saagar_clusters mv,
        #                      base_chemicals bc,
        #                      base_chemical_compounds bcc
        #                  WHERE
        #                     bc.epa_id = bcc.epa_id
        #                 AND bc.epa_id = mv.epa_id
        #                 AND bc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #              "), args=as.list(selected))
        # 
        # anno_hmdb_diseases <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     mv.annotation
        #                  FROM
        #                      base_chemical_hmdb_diseases_saagar_clusters mv,
        #                      base_chemicals bc,
        #                      base_chemical_compounds bcc
        #                  WHERE
        #                     bc.epa_id = bcc.epa_id
        #                 AND bc.epa_id = mv.epa_id
        #                 AND bc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #              "), args=as.list(selected))
        # 
        # anno_hmdb_genes <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     mv.annotation
        #                  FROM
        #                      base_chemical_hmdb_genes_saagar_clusters mv,
        #                      base_chemicals bc,
        #                      base_chemical_compounds bcc
        #                  WHERE
        #                     bc.epa_id = bcc.epa_id
        #                 AND bc.epa_id = mv.epa_id
        #                 AND bc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #              "), args=as.list(selected))
        # 
        # anno_hmdb_tissuelocations <- run_query(TRUE, paste0("
        #                  SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     mv.annotation
        #                  FROM
        #                      base_chemical_hmdb_tissue_locations_saagar_clusters mv,
        #                      base_chemicals bc,
        #                      base_chemical_compounds bcc
        #                  WHERE
        #                     bc.epa_id = bcc.epa_id
        #                 AND bc.epa_id = mv.epa_id
        #                 AND bc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #              "), args=as.list(selected))
        # 
        # anno_ice <- query_ice(chemicals=dtxsids)
        # 
        # anno_invitrodb <- run_query(TRUE, paste0("
        #     SELECT DISTINCT
        #         bcc.preferred_name,
        #         bcc.casrn,
        #         bc.dsstox_substance_id,
        #         isr.invitrodb_values AS assay_name,
        #         iaea.assay_endpoint_attribute,
        #         iaea.assay_endpoint_value
        #     FROM
        #         invitrodb_summarized_results isr,
        #         invitrodb_annotation ia,
        #         invitrodb_assay_endpoint_annotation iaea,
        #         base_chemicals bc,
        #         base_chemical_to_smiles bcs,
        #         base_chemical_compounds bcc
        #     WHERE
        #         bc.epa_id = bcs.epa_id
        #     AND bc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #     AND bc.epa_id = bcc.epa_id
        #     AND bcs.smi_id = isr.smi_id
        #     AND isr.invitrodb_attributes = ia.invitrodb_attributes
        #     AND isr.invitrodb_attributes = 'aenm'
        #     AND isr.invitrodb_values = iaea.aenm
        # "), args=as.list(selected))
        # 
        # anno_leadscope <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     lpm.model_name,
        #                     lcpi.interpretation,
        #                     --lpm.short_description
        #                     lqd.endpoint_desc AS short_description
        #                 FROM
        #                     leadscope_chemical_predictions_interpretation lcpi,
        #                     base_chemical_to_smiles bcs,
        #                     leadscope_predictive_models lpm,
        #                     base_chemical_compounds bcc,
        #                     base_chemicals bc,
        #                     leadscope_qmrf_descriptions lqd
        #                 WHERE
        #                     bcs.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #                 AND bcs.epa_id = bc.epa_id
        #                 AND bcs.epa_id = bcc.epa_id
        #                 AND bcs.smi_id = lcpi.smi_id
        #                 AND lcpi.lmodel_id = lpm.lmodel_id
        #                 AND lcpi.lmodel_id = lqd.lmodel_id
        #             "), args=as.list(selected))
        # 
        # anno_ochem <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     mv.alert_name,
        #                     mv.saresult
        #                 FROM
        #                     mv_ochem_sa_positives_saagar mv,
        #                     base_chemical_to_smiles bcs,
        #                     base_chemical_compounds bcc,
        #                     base_chemicals bc
        #                 WHERE
        #                     bcs.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #                 AND bcs.epa_id = bc.epa_id
        #                 AND bcs.epa_id = bcc.epa_id
        #                 AND bcs.epa_id = mv.epa_id
        #             "), args=as.list(selected))
        # 
        # anno_opera_adme <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     opm.model_name,
        #                     ocp.model_value,
        #                     opm.category,
        #                     opm.model_description
        #                 FROM
        #                     base_chemicals bc,
        #                     base_chemical_to_smiles bcs,
        #                     base_chemical_compounds bcc,
        #                     opera_predictive_models opm,
        #                     opera_chemical_predictions ocp
        #                 WHERE
        #                     bc.epa_id = bcs.epa_id
        #                 AND bcs.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #                 AND bc.epa_id = bcc.epa_id
        #                 AND bcs.smi_id = ocp.smi_id
        #                 AND ocp.omodel_id = opm.omodel_id
        #                 AND opm.category = 'ADME endpoints'
        #             "), args=as.list(selected))
        # 
        # anno_opera_environmental_fate <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     opm.model_name,
        #                     ocp.model_value,
        #                     opm.category,
        #                     opm.model_description
        #                 FROM
        #                     base_chemicals bc,
        #                     base_chemical_to_smiles bcs,
        #                     base_chemical_compounds bcc,
        #                     opera_predictive_models opm,
        #                     opera_chemical_predictions ocp
        #                 WHERE
        #                     bc.epa_id = bcs.epa_id
        #                 AND bcs.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #                 AND bc.epa_id = bcc.epa_id
        #                 AND bcs.smi_id = ocp.smi_id
        #                 AND ocp.omodel_id = opm.omodel_id
        #                 AND opm.category = 'Environmental fate'
        #             "), args=as.list(selected))
        # 
        # anno_opera_physicochem <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     opm.model_name,
        #                     ocp.model_value,
        #                     opm.category,
        #                     opm.model_description
        #                 FROM
        #                     base_chemicals bc,
        #                     base_chemical_to_smiles bcs,
        #                     base_chemical_compounds bcc,
        #                     opera_predictive_models opm,
        #                     opera_chemical_predictions ocp
        #                 WHERE
        #                     bc.epa_id = bcs.epa_id
        #                 AND bcs.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #                 AND bc.epa_id = bcc.epa_id
        #                 AND bcs.smi_id = ocp.smi_id
        #                 AND ocp.omodel_id = opm.omodel_id
        #                 AND opm.category = 'Physicochemical properties'
        #             "), args=as.list(selected))
        # 
        # anno_opera_structural_properties <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     opm.model_name,
        #                     ocp.model_value,
        #                     opm.category,
        #                     opm.model_description
        #                 FROM
        #                     base_chemicals bc,
        #                     base_chemical_to_smiles bcs,
        #                     base_chemical_compounds bcc,
        #                     opera_predictive_models opm,
        #                     opera_chemical_predictions ocp
        #                 WHERE
        #                     bc.epa_id = bcs.epa_id
        #                 AND bcs.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #                 AND bc.epa_id = bcc.epa_id
        #                 AND bcs.smi_id = ocp.smi_id
        #                 AND ocp.omodel_id = opm.omodel_id
        #                 AND opm.category = 'Structural properties'
        #             "), args=as.list(selected))
        # 
        # anno_opera_toxicology_endpoints <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     opm.model_name,
        #                     ocp.model_value,
        #                     opm.category,
        #                     opm.model_description
        #                 FROM
        #                     base_chemicals bc,
        #                     base_chemical_to_smiles bcs,
        #                     base_chemical_compounds bcc,
        #                     opera_predictive_models opm,
        #                     opera_chemical_predictions ocp
        #                 WHERE
        #                     bc.epa_id = bcs.epa_id
        #                 AND bcs.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #                 AND bc.epa_id = bcc.epa_id
        #                 AND bcs.smi_id = ocp.smi_id
        #                 AND ocp.omodel_id = opm.omodel_id
        #                 AND opm.category = 'Toxicology endpoints'
        #             "), args=as.list(selected))
        # 
        # anno_saagar <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     mv.alert_name,
        #                     mv.saresult
        #                 FROM
        #                     mv_saagar_de_positives_saagar mv,
        #                     base_chemical_to_smiles bcs,
        #                     base_chemical_compounds bcc,
        #                     base_chemicals bc
        #                 WHERE
        #                     bcs.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #                 AND bcs.epa_id = bc.epa_id
        #                 AND bcs.epa_id = bcc.epa_id
        #                 AND bcs.epa_id = mv.epa_id
        #             "), args=as.list(selected))
        # 
        # anno_t3db <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     tct.target_id,
        #                     tct.target_name,
        #                     tct.uniprot_id
        #                 FROM
        #                     t3db_chemicals_to_targets tct,
        #                     t3db_to_base_chemicals tbc,
        #                     t3db_curated_chemicals tcc,
        #                     base_chemical_compounds bcc,
        #                     base_chemicals bc
        #                 WHERE
        #                     tbc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #                 AND bc.epa_id = bcc.epa_id
        #                 AND tbc.epa_id = bc.epa_id
        #                 AND tbc.t3db_id = tcc.t3db_id
        #                 AND tcc.t3_id = tct.t3_id
        #             "), args=as.list(selected))
        # 
        # anno_toxrefdb_neoplastic <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     CONCAT(mv.disease_type, '_',
        #                     mv.study_type, '_',
        #                     mv.species, '_',
        #                     mv.sex, '_',
        #                     mv.life_stage, '_',
        #                     mv.endpoint_type, '_',
        #                     mv.endpoint_target) AS annotation
        #                 FROM
        #                     mv_toxrefdb_cancer_saagar mv,
        #                     base_chemical_compounds bcc,
        #                     base_chemicals bc
        #                 WHERE
        #                     mv.epa_id = bcc.epa_id
        #                 AND mv.epa_id = bc.epa_id
        #                 AND mv.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #             "), args=as.list(selected))
        # 
        # anno_toxrefdb_nonneoplastic <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     CONCAT(mv.disease_type, '_',
        #                     mv.study_type, '_',
        #                     mv.species, '_',
        #                     mv.sex, '_',
        #                     mv.life_stage, '_',
        #                     mv.endpoint_type, '_',
        #                     mv.endpoint_target) AS annotation
        #                 FROM
        #                     mv_toxrefdb_nonneoplastic_saagar mv,
        #                     base_chemical_compounds bcc,
        #                     base_chemicals bc
        #                 WHERE
        #                     mv.epa_id = bcc.epa_id
        #                 AND mv.epa_id = bc.epa_id
        #                 AND mv.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #             "), args=as.list(selected))
        # 
        # anno_pubchem <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     pb.bioassay_name,
        #                     pb.source_name,
        #                     pb.source_id,
        #                     pb.substance_type,
        #                     pb.outcome_type,
        #                     pb.project_category,
        #                     pb.bioassay_group,
        #                     pb.bioassay_types,
        #                     pb.protein_accessions,
        #                     pbm.aid,
        #                     pbm.sid,
        #                     pbm.sid_group,
        #                     pbm.activity_outcome,
        #                     pbm.activity_name,
        #                     pbm.activity_qualifier,
        #                     pbm.activity_value,
        #                     pbm.protein_accession,
        #                     pbm.gene_id,
        #                     pbm.target_tax_id,
        #                     pbm.pmid
        #                 FROM
        #                     base_chemical_to_pubchem_cid bcpc,
        #                     pubchem_bioactivity_mapping pbm,
        #                     pubchem_bioassays pb,
        #                     base_chemicals bc,
        #                     base_chemical_compounds bcc
        #                 WHERE
        #                     bcpc.pubchem_cid = pbm.cid
        #                 AND bcpc.epa_id = bc.epa_id
        #                 AND bcpc.epa_id = bcc.epa_id
        #                 AND pbm.aid = pb.aid
        #                 AND bcpc.epa_id IN (
        #                     ", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #             "), args=as.list(selected))
   
        anno_pubchem_props <- run_query(TRUE, paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            pba.pubchem_attribute,
                            pba.pubchem_value
                        FROM
                            base_chemical_to_pubchem_cid bcpc,
                            pubchem_attributes pba,
                            base_chemicals bc,
                            base_chemical_compounds bcc
                        WHERE
                            bcpc.pubchem_cid = pba.cid
                        AND bcpc.epa_id = bc.epa_id
                        AND bcpc.epa_id = bcc.epa_id
                        AND bcpc.epa_id IN (
                            ", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                    "), args=as.list(selected))
    
        # anno_gras <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     bcg.gras_substance,
        #                     bcg.scogs_report_number,
        #                     bcg.year,
        #                     bcg.scogs_conclusion,
        #                     CONCAT(gsc.scogs_conclusion, ': ', gsc.scogs_interpretation) AS scogs_interpretation,
        #                     bcg.\"21_cfr_regulation\",
        #                     bcg.ntis_accession
        #                 FROM
        #                     base_chemical_to_gras bcg,
        #                     base_chemicals bc,
        #                     base_chemical_compounds bcc,
        #                     gras_scogs_conclusions gsc
        #                 WHERE
        #                     bcg.epa_id = bc.epa_id
        #                 AND bcg.scogs_conclusion LIKE '%' || gsc.scogs_conclusion || '%'
        #                 AND bcg.epa_id = bcc.epa_id
        #                 AND bcg.epa_id IN (
        #                     ", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #             "), args=as.list(selected))
        # 
        # anno_foodb_enzymes <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     fce.enzyme_name,
        #                     fce.gene_name,
        #                     fce.kingdom,
        #                     fce.superclass,
        #                     fce.class,
        #                     fce.subclass,
        #                     fce.citations
        #                 FROM
        #                     base_chemical_to_foodb_compound bcfc,
        #                     foodb_compound_enzymes fce,
        #                     base_chemicals bc,
        #                     base_chemical_compounds bcc
        #                 WHERE
        #                     bcfc.epa_id = bc.epa_id
        #                 AND bcfc.foodb_compound_id = fce.compound_id
        #                 AND bcfc.epa_id = bcc.epa_id
        #                 AND bcfc.epa_id IN (
        #                     ", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #             "), args=as.list(selected))
        # 
        # anno_foodb_flavors <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     fcf.flavor_name,
        #                     fcf.flavor_group,
        #                     fcf.category,
        #                     fcf.kingdom,
        #                     fcf.superclass,
        #                     fcf.class,
        #                     fcf.subclass,
        #                     fcf.citations
        #                 
        #                 FROM
        #                     base_chemical_to_foodb_compound bcfc,
        #                     foodb_compound_flavors fcf,
        #                     base_chemicals bc,
        #                     base_chemical_compounds bcc
        #                 WHERE
        #                     bcfc.epa_id = bc.epa_id
        #                 AND bcfc.foodb_compound_id = fcf.compound_id
        #                 AND bcfc.epa_id = bcc.epa_id
        #                 AND bcfc.epa_id IN (
        #                     ", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #             "), args=as.list(selected))
        # 
        # anno_foodb_foodcontent <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     fcfc.name,
        #                     fcfc.name_scientific,
        #                     fcfc.food_group,
        #                     fcfc.food_subgroup,
        #                     fcfc.food_type,
        #                     fcfc.description,
        #                     fcfc.kingdom,
        #                     fcfc.superclass,
        #                     fcfc.class,
        #                     fcfc.subclass
        #                 FROM
        #                     base_chemical_to_foodb_compound bcfc,
        #                     foodb_compound_food_content fcfc,
        #                     base_chemicals bc,
        #                     base_chemical_compounds bcc
        #                 WHERE
        #                     bcfc.epa_id = bc.epa_id
        #                 AND bcfc.foodb_compound_id = fcfc.compound_id
        #                 AND bcfc.epa_id = bcc.epa_id
        #                 AND bcfc.epa_id IN (
        #                     ", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #             "), args=as.list(selected))
        # 
        # anno_foodb_healtheffects <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     fche.health_effect_name,
        #                     fche.kingdom,
        #                     fche.superclass,
        #                     fche.class,
        #                     fche.subclass,
        #                     fche.description
        #                 FROM
        #                     base_chemical_to_foodb_compound bcfc,
        #                     foodb_compound_health_effects fche,
        #                     base_chemicals bc,
        #                     base_chemical_compounds bcc
        #                 WHERE
        #                     bcfc.epa_id = bc.epa_id
        #                 AND bcfc.foodb_compound_id = fche.compound_id
        #                 AND bcfc.epa_id = bcc.epa_id
        #                 AND bcfc.epa_id IN (
        #                     ", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #             "), args=as.list(selected))
        # 
        # anno_foodb_ontology <- run_query(TRUE, paste0("
        #                 SELECT DISTINCT
        #                     bcc.preferred_name,
        #                     bcc.casrn,
        #                     bc.dsstox_substance_id,
        #                     fcot.term,
        #                     fcot.definition,
        #                     fcot.level,
        #                     fcot.kingdom,
        #                     fcot.superclass,
        #                     fcot.class,
        #                     fcot.subclass
        #                 FROM
        #                     base_chemical_to_foodb_compound bcfc,
        #                     foodb_compound_ontology_terms fcot,
        #                     base_chemicals bc,
        #                     base_chemical_compounds bcc
        #                 WHERE
        #                     bcfc.epa_id = bc.epa_id
        #                 AND bcfc.foodb_compound_id = fcot.compound_id
        #                 AND bcfc.epa_id = bcc.epa_id
        #                 AND bcfc.epa_id IN (
        #                     ", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
        #             "), args=as.list(selected))
    
        anno_epa_props <- run_query(TRUE, paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
                            ep.epa_attribute,
                            ep.epa_attribute_value
                        FROM
                            epa_properties ep,
                            base_chemicals bc,
                            base_chemical_compounds bcc
                        WHERE
                            ep.epa_id = bc.epa_id
                        AND ep.epa_id = bcc.epa_id
                        AND ep.epa_id IN (
                            ", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                    "), args=as.list(selected))
    
    
    return(list(
        anno_admet_binr=anno_admet_binr,
        anno_admet_catg=anno_admet_catg,
        anno_admet_cont=anno_admet_cont,
        anno_chembl=anno_chembl,
        anno_cpd=anno_cpd,
        anno_ctd_bioprocess=anno_ctd_bioprocess,
        anno_ctd_cellcomp=anno_ctd_cellcomp,
        anno_ctd_diseases=anno_ctd_diseases,
        anno_ctd_genes=anno_ctd_genes,
        anno_ctd_molfunct=anno_ctd_molfunct,
        anno_ctd_phenotypes=anno_ctd_phenotypes,
        anno_drugbank_atccodes=anno_drugbank_atccodes,
        anno_drugbank_carriers=anno_drugbank_carriers,
        anno_drugbank_enzymes=anno_drugbank_enzymes,
        anno_drugbank_targets=anno_drugbank_targets,
        anno_drugbank_transporters=anno_drugbank_transporters,
        anno_hmdb_biospecimenlocations=anno_hmdb_biospecimenlocations,
        anno_hmdb_cellularlocations=anno_hmdb_cellularlocations,
        anno_hmdb_diseases=anno_hmdb_diseases,
        anno_hmdb_genes=anno_hmdb_genes,
        anno_hmdb_tissuelocations=anno_hmdb_tissuelocations,
        anno_ice=anno_ice,
        anno_invitrodb=anno_invitrodb,
        anno_leadscope=anno_leadscope,
        anno_ochem=anno_ochem,
        anno_opera_adme=anno_opera_adme,
        anno_opera_environmental_fate=anno_opera_environmental_fate,
        anno_opera_physicochem=anno_opera_physicochem,
        anno_opera_structural_properties=anno_opera_structural_properties,
        anno_opera_toxicology_endpoints=anno_opera_toxicology_endpoints,
        anno_saagar=anno_saagar,
        anno_t3db=anno_t3db,
        anno_toxrefdb_neoplastic=anno_toxrefdb_neoplastic,
        anno_toxrefdb_nonneoplastic=anno_toxrefdb_nonneoplastic,
        anno_pubchem=anno_pubchem,
        anno_pubchem_props=anno_pubchem_props,
        anno_gras=anno_gras,
        anno_foodb_enzymes=anno_foodb_enzymes,
        anno_foodb_flavors=anno_foodb_flavors,
        anno_foodb_foodcontent=anno_foodb_foodcontent,
        anno_foodb_healtheffects=anno_foodb_healtheffects,
        anno_foodb_ontology=anno_foodb_ontology,
        anno_epa_props=anno_epa_props
    ))
}