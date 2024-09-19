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



# load_chembl <- function(epa_ids){
#     epa_df  <- run_query_chembl(paste0("
#         SELECT
# 	XCR.molregno,
# 	chemical_name,
# 	mw_freebase,
# 	alogp,
# 	hba,
# 	hbd,
# 	psa,
# 	rtb,
# 	ro3_pass,
# 	num_ro5_violations,
# 	cx_most_apka,
# 	cx_most_bpka,
# 	cx_logp,
# 	cx_logd,
# 	molecular_species,
# 	full_mwt,
# 	aromatic_rings,
# 	heavy_atoms,
# 	qed_weighted,
# 	mw_monoisotopic,
# 	full_molformula,
# 	hba_lipinski,
# 	hbd_lipinski,
# 	num_lipinski_ro5_violations,
# 	np_likeness_score,
# 	standard_inchi,
# 	standard_inchi_key,
# 	canonical_smiles,
# 	structural_alerts,
# 	study_title,
# 	study_doi,
# 	study_doctype,
# 	study_authors,
# 	study_abstract,
# 	mechanism_of_action,
# 	mechanism_comment,
# 	efo_term,
# 	warning_type,
# 	warning_class,
# 	src_description,
# 	chembl_id,
# 	activities
# FROM
# 	(
# 		SELECT cr.molregno, cr.compound_name AS chemical_name, cr.doc_id, cr.record_id, cr.src_id
# 		FROM compound_records cr
# 		WHERE cr.record_id = 1343079
# 	) XCR
# 	
# 	LEFT JOIN
# 	(
# 		SELECT 
# 			cp.molregno,
# 			cp.mw_freebase,
# 			cp.alogp,
# 			cp.hba,
# 			cp.hbd,
# 			cp.psa,
# 			cp.rtb,
# 			cp.ro3_pass,
# 			cp.num_ro5_violations,
# 			cp.cx_most_apka,
# 			cp.cx_most_bpka,
# 			cp.cx_logp,
# 			cp.cx_logd,
# 			cp.molecular_species,
# 			cp.full_mwt,
# 			cp.aromatic_rings,
# 			cp.heavy_atoms,
# 			cp.qed_weighted,
# 			cp.mw_monoisotopic,
# 			cp.full_molformula,
# 			cp.hba_lipinski,
# 			cp.hbd_lipinski,
# 			cp.num_lipinski_ro5_violations,
# 			cp.np_likeness_score
# 		FROM compound_properties cp
# 	) XCP ON XCR.molregno = XCP.molregno
# 	
# 	LEFT JOIN
# 	(
# 		SELECT
# 			cs.molregno,
# 			cs.standard_inchi,
# 			cs.standard_inchi_key,
# 			cs.canonical_smiles
# 		FROM compound_structures cs
# 	) XCS ON XCR.molregno = XCS.molregno
# 	
# 	LEFT JOIN
# 	(
# 		SELECT
# 			csa.molregno,
# 			string_agg(sa.alert_name, '; ') AS structural_alerts
# 		FROM
# 			compound_structural_alerts csa,
# 			structural_alerts sa
# 		WHERE csa.alert_id = sa.alert_id
# 		GROUP BY 1
# 	) XSA ON XCR.molregno = XSA.molregno
# 	
# 	LEFT JOIN
# 	(
# 		SELECT
# 			d.doc_id,
# 			d.title AS study_title,
# 			d.doi AS study_doi,
# 			d.doc_type AS study_doctype,
# 			d.authors AS study_authors,
# 			d.abstract AS study_abstract
# 		FROM docs d
# 	) XD ON XCR.doc_id = XD.doc_id
# 	
# 	LEFT JOIN
# 	(
# 		SELECT
# 			dm.record_id,
# 			dm.mechanism_of_action,
# 			dm.mechanism_comment
# 		FROM drug_mechanism dm
# 	) XDM ON XCR.record_id = XDM.record_id
# 	
# 	LEFT JOIN
# 	(
# 		SELECT
# 			di.record_id,
# 			di.efo_term
# 		FROM drug_indication di
# 	) XDI ON XCR.record_id = XDI.record_id
# 	
# 	LEFT JOIN
# 	(
# 		SELECT
# 			dw.record_id,
# 			dw.warning_type,
# 			dw.warning_class
# 		FROM drug_warning dw
# 	) XDW ON XCR.record_id = XDW.record_id
# 	
# 	LEFT JOIN
# 	(
# 		SELECT
# 			s.src_id,
# 			s.src_description
# 		FROM "source" s
# 	) XS ON XCR.src_id = XS.src_id
# 	
# 	LEFT JOIN
# 	(
# 		SELECT
# 			md.molregno,
# 			md.chembl_id
# 		FROM molecule_dictionary md
# 	) XMD ON XCR.molregno = XMD.molregno
# 	
# 	LEFT JOIN
# 	(
# 		SELECT
# 			a.record_id,
# 			string_agg(CONCAT(a.standard_type, ' ', a.standard_relation, ' ', a.standard_value, a.standard_units), '; ') AS activities
# 		FROM activities a
# 		GROUP BY 1
# 	) XA ON XCR.record_id = XA.record_id
# 
# 
# 
#         
#     "), args=epa_ids)
# }


# Convenient handler for loading annotations
load_annotations <- function(selected_saagar_clusters, selected_node, input, reactive_clusters, reactive_selected_cluster, mode, type=""){
    tmp_usage_df <- unique(selected_saagar_clusters[[selected_node]][, c("preferred_name")])
    epa_ids <- selected_saagar_clusters[[selected_node]]
    num_chemicals <- length(unique(epa_ids[, "epa_id"]))
    dtxsids <- unique(epa_ids[, "dsstox_substance_id"])
    tmp <- NULL
    final_plots <- NULL
    
    if(mode == "cluster"){
        tmp <- cluster_annotations(selected=as.integer(epa_ids$epa_id), clust_l3=epa_ids[1, "clust_l3"], input=input, reactive_clusters=reactive_clusters, reactive_selected_cluster=reactive_selected_cluster, dtxsids=dtxsids, switch_mode=type)
        topn_admet_binr <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_admet_binr) <- c("annotation", "percentage", "num_of_chemicals")
        topn_chembl <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_chembl) <- c("annotation", "percentage", "num_of_chemicals")
        topn_cpd <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_cpd) <- c("annotation", "percentage", "num_of_chemicals")
        topn_ctd_bioprocess <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_ctd_bioprocess) <- c("annotation", "percentage", "num_of_chemicals")
        topn_ctd_cellcomp <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_ctd_cellcomp) <- c("annotation", "percentage", "num_of_chemicals")
        topn_ctd_diseases <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_ctd_diseases) <- c("annotation", "percentage", "num_of_chemicals")
        topn_ctd_genes <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_ctd_genes) <- c("annotation", "percentage", "num_of_chemicals")
        topn_ctd_molfunct <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_ctd_molfunct) <- c("annotation", "percentage", "num_of_chemicals")
        topn_ctd_phenotypes <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_ctd_phenotypes) <- c("annotation", "percentage", "num_of_chemicals")
        topn_drugbank_atccodes <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_drugbank_atccodes) <- c("annotation", "percentage", "num_of_chemicals")
        topn_drugbank_carriers <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_drugbank_carriers) <- c("annotation", "percentage", "num_of_chemicals")
        topn_drugbank_enzymes <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_drugbank_enzymes) <- c("annotation", "percentage", "num_of_chemicals")
        topn_drugbank_targets <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_drugbank_targets) <- c("annotation", "percentage", "num_of_chemicals")
        topn_drugbank_transporters <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_drugbank_transporters) <- c("annotation", "percentage", "num_of_chemicals")
        topn_hmdb_biospecimenlocations <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_hmdb_biospecimenlocations) <- c("annotation", "percentage", "num_of_chemicals")
        topn_hmdb_cellularlocations <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_hmdb_cellularlocations) <- c("annotation", "percentage", "num_of_chemicals")
        topn_hmdb_diseases <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_hmdb_diseases) <- c("annotation", "percentage", "num_of_chemicals")
        topn_hmdb_genes <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_hmdb_genes) <- c("annotation", "percentage", "num_of_chemicals")
        topn_hmdb_tissuelocations <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_hmdb_tissuelocations) <- c("annotation", "percentage", "num_of_chemicals")
        topn_ice <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_ice) <- c("annotation", "percentage", "num_of_chemicals")
        topn_invitrodb <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_invitrodb) <- c("annotation", "percentage", "num_of_chemicals")
        topn_leadscope <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_leadscope) <- c("annotation", "percentage", "num_of_chemicals")
        topn_ochem <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_ochem) <- c("annotation", "percentage", "num_of_chemicals")
        topn_opera_adme <- data.frame(matrix(ncol=6, nrow=0))
        colnames(topn_opera_adme) <- c("annotation", "percentage", "num_of_chemicals")
        topn_opera_environmental_fate <- data.frame(matrix(ncol=6, nrow=0))
        colnames(topn_opera_environmental_fate) <- c("annotation", "percentage", "num_of_chemicals")
        topn_opera_physicochem <- data.frame(matrix(ncol=6, nrow=0))
        colnames(topn_opera_physicochem) <- c("annotation", "percentage", "num_of_chemicals")
        topn_opera_structural_properties <- data.frame(matrix(ncol=6, nrow=0))
        colnames(topn_opera_structural_properties) <- c("annotation", "percentage", "num_of_chemicals")
        topn_opera_toxicology_endpoints <- data.frame(matrix(ncol=6, nrow=0))
        colnames(topn_opera_toxicology_endpoints) <- c("annotation", "percentage", "num_of_chemicals")
        topn_saagar <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_saagar) <- c("annotation", "percentage", "num_of_chemicals")
        topn_t3db <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_t3db) <- c("annotation", "percentage", "num_of_chemicals")
        topn_toxrefdb_neoplastic <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_toxrefdb_neoplastic) <- c("annotation", "percentage", "num_of_chemicals")
        topn_toxrefdb_nonneoplastic <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_toxrefdb_nonneoplastic) <- c("annotation", "percentage", "num_of_chemicals")
        
        topn_pubchem <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_pubchem) <- c("annotation", "percentage", "num_of_chemicals")
        
        topn_pubchem_props <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_pubchem_props) <- c("annotation", "percentage", "num_of_chemicals")
        
        topn_gras <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_gras) <- c("annotation", "percentage", "num_of_chemicals")
        
        topn_foodb_enzymes <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_foodb_enzymes) <- c("annotation", "percentage", "num_of_chemicals")
        topn_foodb_flavors <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_foodb_flavors) <- c("annotation", "percentage", "num_of_chemicals")
        topn_foodb_foodcontent <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_foodb_foodcontent) <- c("annotation", "percentage", "num_of_chemicals")
        topn_foodb_healtheffects <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_foodb_healtheffects) <- c("annotation", "percentage", "num_of_chemicals")
        topn_foodb_ontology <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_foodb_ontology) <- c("annotation", "percentage", "num_of_chemicals")
        
        topn_epa_props <- data.frame(matrix(ncol=3, nrow=0))
        colnames(topn_epa_props) <- c("annotation", "percentage", "num_of_chemicals")
        
        if(nrow(tmp[["anno_admet_binr"]]) > 0){
            topn_admet_binr <- tmp[["anno_admet_binr"]][, c("model_name", "positives")]
            topn_admet_binr$num_of_chemicals <- lapply(topn_admet_binr$positives, function(x) paste0(x, "/", num_chemicals))
            topn_admet_binr$positives <- lapply(topn_admet_binr$positives, function(x) x/num_chemicals)
            colnames(topn_admet_binr) <- c("annotation", "percentage", "num_of_chemicals")
        }
        
        if(nrow(tmp[["anno_chembl"]]) > 0){
            topn_chembl <- tmp[["anno_chembl"]][, c("alert_name", "alert_num")]
            topn_chembl_split <- split(topn_chembl, topn_chembl$alert_name)
            topn_chembl_split_num_of_chemicals <- unlist(lapply(topn_chembl_split, function(x) paste0("Avg. ", sum(x$alert_num)/nrow(x), "/", num_chemicals)))
            topn_chembl_split <- unlist(lapply(topn_chembl_split, function(x) sum(unlist(lapply(x$alert_num, function(y) y/num_chemicals)))/nrow(x)))
            topn_chembl <- data.frame(annotation=names(topn_chembl_split), percentage=topn_chembl_split, num_of_chemicals=topn_chembl_split_num_of_chemicals, stringsAsFactors=FALSE)
        }
        
        if(nrow(tmp[["anno_cpd"]]) > 0){
            topn_cpd <- tmp[["anno_cpd"]][, c("product_name", "product_num")]
            topn_cpd$num_of_chemicals <- unlist(lapply(topn_cpd$product_num, function(x) paste0(x, "/", num_chemicals)))
            topn_cpd$product_num <- lapply(topn_cpd$product_num, function(x) x/num_chemicals)
            colnames(topn_cpd) <- c("annotation", "percentage", "num_of_chemicals")
        }
        
        if(nrow(tmp[["anno_ctd_bioprocess"]]) > 0){
            topn_ctd_bioprocess <- tmp[["anno_ctd_bioprocess"]][, c("term_name", "term_num")]
            topn_ctd_bioprocess$num_of_chemicals <- unlist(lapply(topn_ctd_bioprocess$term_num, function(x) paste0(x, "/", num_chemicals)))
            topn_ctd_bioprocess$term_num <- lapply(topn_ctd_bioprocess$term_num, function(x) x/num_chemicals)
            colnames(topn_ctd_bioprocess) <- c("annotation", "percentage", "num_of_chemicals")
        }
        if(nrow(tmp[["anno_ctd_cellcomp"]]) > 0){
            topn_ctd_cellcomp <- tmp[["anno_ctd_cellcomp"]][, c("term_name", "term_num")]
            topn_ctd_cellcomp$num_of_chemicals <- unlist(lapply(topn_ctd_cellcomp$term_num, function(x) paste0(x, "/", num_chemicals)))
            topn_ctd_cellcomp$term_num <- lapply(topn_ctd_cellcomp$term_num, function(x) x/num_chemicals)
            colnames(topn_ctd_cellcomp) <- c("annotation", "percentage", "num_of_chemicals")
        }
        if(nrow(tmp[["anno_ctd_diseases"]]) > 0){
            topn_ctd_diseases <- tmp[["anno_ctd_diseases"]][, c("term_name", "term_num")]
            topn_ctd_diseases$num_of_chemicals <- unlist(lapply(topn_ctd_diseases$term_num, function(x) paste0(x, "/", num_chemicals)))
            topn_ctd_diseases$term_num <- lapply(topn_ctd_diseases$term_num, function(x) x/num_chemicals)
            colnames(topn_ctd_diseases) <- c("annotation", "percentage", "num_of_chemicals")
        }
        if(nrow(tmp[["anno_ctd_genes"]]) > 0){
            topn_ctd_genes <- tmp[["anno_ctd_genes"]][, c("term_name", "term_num")]
            topn_ctd_genes$num_of_chemicals <- unlist(lapply(topn_ctd_genes$term_num, function(x) paste0(x, "/", num_chemicals)))
            topn_ctd_genes$term_num <- lapply(topn_ctd_genes$term_num, function(x) x/num_chemicals)
            colnames(topn_ctd_genes) <- c("annotation", "percentage", "num_of_chemicals")
        }
        if(nrow(tmp[["anno_ctd_molfunct"]]) > 0){
            topn_ctd_molfunct <- tmp[["anno_ctd_molfunct"]][, c("term_name", "term_num")]
            topn_ctd_molfunct$num_of_chemicals <- unlist(lapply(topn_ctd_molfunct$term_num, function(x) paste0(x, "/", num_chemicals)))
            topn_ctd_molfunct$term_num <- lapply(topn_ctd_molfunct$term_num, function(x) x/num_chemicals)
            colnames(topn_ctd_molfunct) <- c("annotation", "percentage", "num_of_chemicals")
        }
        if(nrow(tmp[["anno_ctd_phenotypes"]]) > 0){
            topn_ctd_phenotypes <- tmp[["anno_ctd_phenotypes"]][, c("term_name", "term_num")]
            topn_ctd_phenotypes$num_of_chemicals <- unlist(lapply(topn_ctd_phenotypes$term_num, function(x) paste0(x, "/", num_chemicals)))
            topn_ctd_phenotypes$term_num <- lapply(topn_ctd_phenotypes$term_num, function(x) x/num_chemicals)
            colnames(topn_ctd_phenotypes) <- c("annotation", "percentage", "num_of_chemicals")
        }
        
        if(nrow(tmp[["anno_drugbank_atccodes"]]) > 0){
            topn_drugbank_atccodes <- tmp[["anno_drugbank_atccodes"]][, c("annotation_name", "annotation_num")]
            topn_drugbank_atccodes_split <- split(topn_drugbank_atccodes, topn_drugbank_atccodes$annotation_name)
            topn_drugbank_atccodes_num_of_chemicals <- unlist(lapply(topn_drugbank_atccodes_split, function(x) paste0("Avg. ", sum(x$annotation_num)/nrow(x), "/", num_chemicals)))
            topn_drugbank_atccodes_split <- unlist(lapply(topn_drugbank_atccodes_split, function(x) sum(unlist(lapply(x$annotation_num, function(y) y/num_chemicals)))/nrow(x)))
            topn_drugbank_atccodes <- data.frame(annotation=names(topn_drugbank_atccodes_split), percentage=topn_drugbank_atccodes_split, num_of_chemicals=topn_drugbank_atccodes_num_of_chemicals, stringsAsFactors=FALSE)
        }
        
        if(nrow(tmp[["anno_drugbank_carriers"]]) > 0){
            topn_drugbank_carriers <- tmp[["anno_drugbank_carriers"]][, c("annotation_name", "annotation_num")]
            topn_drugbank_carriers_split <- split(topn_drugbank_carriers, topn_drugbank_carriers$annotation_name)
            topn_drugbank_carriers_num_of_chemicals <- unlist(lapply(topn_drugbank_carriers_split, function(x) paste0("Avg. ", sum(x$annotation_num)/nrow(x), "/", num_chemicals)))
            topn_drugbank_carriers_split <- unlist(lapply(topn_drugbank_carriers_split, function(x) sum(unlist(lapply(x$annotation_num, function(y) y/num_chemicals)))/nrow(x)))
            topn_drugbank_carriers <- data.frame(annotation=names(topn_drugbank_carriers_split), percentage=topn_drugbank_carriers_split, num_of_chemicals=topn_drugbank_carriers_split, stringsAsFactors=FALSE)
        }
        
        if(nrow(tmp[["anno_drugbank_enzymes"]]) > 0){
            topn_drugbank_enzymes <- tmp[["anno_drugbank_enzymes"]][, c("annotation_name", "annotation_num")]
            topn_drugbank_enzymes_split <- split(topn_drugbank_enzymes, topn_drugbank_enzymes$annotation_name)
            topn_drugbank_enzymes_num_of_chemicals <- unlist(lapply(topn_drugbank_enzymes_split, function(x) paste0("Avg. ", sum(x$annotation_num)/nrow(x), "/", num_chemicals)))
            topn_drugbank_enzymes_split <- unlist(lapply(topn_drugbank_enzymes_split, function(x) sum(unlist(lapply(x$annotation_num, function(y) y/num_chemicals)))/nrow(x)))
            topn_drugbank_enzymes <- data.frame(annotation=names(topn_drugbank_enzymes_split), percentage=topn_drugbank_enzymes_split, num_of_chemicals=topn_drugbank_enzymes_num_of_chemicals, stringsAsFactors=FALSE)
        }
        
        if(nrow(tmp[["anno_drugbank_targets"]]) > 0){
            topn_drugbank_targets <- tmp[["anno_drugbank_targets"]][, c("annotation_name", "annotation_num")]
            topn_drugbank_targets_split <- split(topn_drugbank_targets, topn_drugbank_targets$annotation_name)
            topn_drugbank_targets_num_of_chemicals <- unlist(lapply(topn_drugbank_targets_split, function(x) paste0("Avg. ", sum(x$annotation_num)/nrow(x), "/", num_chemicals)))
            topn_drugbank_targets_split <- unlist(lapply(topn_drugbank_targets_split, function(x) sum(unlist(lapply(x$annotation_num, function(y) y/num_chemicals)))/nrow(x)))
            topn_drugbank_targets <- data.frame(annotation=names(topn_drugbank_targets_split), percentage=topn_drugbank_targets_split, num_of_chemicals=topn_drugbank_targets_num_of_chemicals, stringsAsFactors=FALSE)
        }
        if(nrow(tmp[["anno_drugbank_transporters"]]) > 0){
            topn_drugbank_transporters <- tmp[["anno_drugbank_transporters"]][, c("annotation_name", "annotation_num")]
            topn_drugbank_transporters_split <- split(topn_drugbank_transporters, topn_drugbank_transporters$annotation_name)
            topn_drugbank_transporters_num_of_chemicals <- unlist(lapply(topn_drugbank_transporters_split, function(x) paste0("Avg. ", sum(x$annotation_num)/nrow(x), "/", num_chemicals)))
            topn_drugbank_transporters_split <- unlist(lapply(topn_drugbank_transporters_split, function(x) sum(unlist(lapply(x$annotation_num, function(y) y/num_chemicals)))/nrow(x)))
            topn_drugbank_transporters <- data.frame(annotation=names(topn_drugbank_transporters_split), percentage=topn_drugbank_transporters_split, num_of_chemicals=topn_drugbank_transporters_num_of_chemicals, stringsAsFactors=FALSE)
        }
        
        if(nrow(tmp[["anno_hmdb_biospecimenlocations"]]) > 0){
            topn_hmdb_biospecimenlocations <- tmp[["anno_hmdb_biospecimenlocations"]][, c("annotation_name", "annotation_num")]
            topn_hmdb_biospecimenlocations_split <- split(topn_hmdb_biospecimenlocations, topn_hmdb_biospecimenlocations$annotation_name)
            topn_hmdb_biospecimenlocations_num_of_chemicals <- unlist(lapply(topn_hmdb_biospecimenlocations_split, function(x) paste0("Avg. ", sum(x$annotation_num)/nrow(x), "/", num_chemicals)))
            topn_hmdb_biospecimenlocations_split <- unlist(lapply(topn_hmdb_biospecimenlocations_split, function(x) sum(unlist(lapply(x$annotation_num, function(y) y/num_chemicals)))/nrow(x)))
            topn_hmdb_biospecimenlocations <- data.frame(annotation=names(topn_hmdb_biospecimenlocations_split), percentage=topn_hmdb_biospecimenlocations_split, num_of_chemicals=topn_hmdb_biospecimenlocations_num_of_chemicals, stringsAsFactors=FALSE)
        }
        if(nrow(tmp[["anno_hmdb_cellularlocations"]]) > 0){
            topn_hmdb_cellularlocations <- tmp[["anno_hmdb_cellularlocations"]][, c("annotation_name", "annotation_num")]
            topn_hmdb_cellularlocations_split <- split(topn_hmdb_cellularlocations, topn_hmdb_cellularlocations$annotation_name)
            topn_hmdb_cellularlocations_num_of_chemicals <- unlist(lapply(topn_hmdb_cellularlocations_split, function(x) paste0("Avg. ", sum(x$annotation_num)/nrow(x), "/", num_chemicals)))
            topn_hmdb_cellularlocations_split <- unlist(lapply(topn_hmdb_cellularlocations_split, function(x) sum(unlist(lapply(x$annotation_num, function(y) y/num_chemicals)))/nrow(x)))
            topn_hmdb_cellularlocations <- data.frame(annotation=names(topn_hmdb_cellularlocations_split), percentage=topn_hmdb_cellularlocations_split, num_of_chemicals=topn_hmdb_cellularlocations_num_of_chemicals, stringsAsFactors=FALSE)
        }
        if(nrow(tmp[["anno_hmdb_diseases"]]) > 0){
            topn_hmdb_diseases <- tmp[["anno_hmdb_diseases"]][, c("annotation_name", "annotation_num")]
            topn_hmdb_diseases_split <- split(topn_hmdb_diseases, topn_hmdb_diseases$annotation_name)
            topn_hmdb_diseases_num_of_chemicals <- unlist(lapply(topn_hmdb_diseases_split, function(x) paste0("Avg. ", sum(x$annotation_num)/nrow(x), "/", num_chemicals)))
            topn_hmdb_diseases_split <- unlist(lapply(topn_hmdb_diseases_split, function(x) sum(unlist(lapply(x$annotation_num, function(y) y/num_chemicals)))/nrow(x)))
            topn_hmdb_diseases <- data.frame(annotation=names(topn_hmdb_diseases_split), percentage=topn_hmdb_diseases_split, num_of_chemicals=topn_hmdb_diseases_num_of_chemicals, stringsAsFactors=FALSE)
        }
        if(nrow(tmp[["anno_hmdb_genes"]]) > 0){
            topn_hmdb_genes <- tmp[["anno_hmdb_genes"]][, c("annotation_name", "annotation_num")]
            topn_hmdb_genes_split <- split(topn_hmdb_genes, topn_hmdb_genes$annotation_name)
            topn_hmdb_genes_num_of_chemicals <- unlist(lapply(topn_hmdb_genes_split, function(x) paste0("Avg. ", sum(x$annotation_num)/nrow(x), "/", num_chemicals)))
            topn_hmdb_genes_split <- unlist(lapply(topn_hmdb_genes_split, function(x) sum(unlist(lapply(x$annotation_num, function(y) y/num_chemicals)))/nrow(x)))
            topn_hmdb_genes <- data.frame(annotation=names(topn_hmdb_genes_split), percentage=topn_hmdb_genes_split, num_of_chemicals=topn_hmdb_genes_num_of_chemicals, stringsAsFactors=FALSE)
        }
        if(nrow(tmp[["anno_hmdb_tissuelocations"]]) > 0){
            topn_hmdb_tissuelocations <- tmp[["anno_hmdb_tissuelocations"]][, c("annotation_name", "annotation_num")]
            topn_hmdb_tissuelocations_split <- split(topn_hmdb_tissuelocations, topn_hmdb_tissuelocations$annotation_name)
            topn_hmdb_tissuelocations_num_of_chemicals <- unlist(lapply(topn_hmdb_tissuelocations_split, function(x) paste0("Avg. ", sum(x$annotation_num)/nrow(x), "/", num_chemicals)))
            topn_hmdb_tissuelocations_split <- unlist(lapply(topn_hmdb_tissuelocations_split, function(x) sum(unlist(lapply(x$annotation_num, function(y) y/num_chemicals)))/nrow(x)))
            topn_hmdb_tissuelocations <- data.frame(annotation=names(topn_hmdb_tissuelocations_split), percentage=topn_hmdb_tissuelocations_split, num_of_chemicals=topn_hmdb_tissuelocations_num_of_chemicals, stringsAsFactors=FALSE)
        }
        
        if(nrow(tmp[["anno_ice"]]) > 0){
            topn_ice <- tmp[["anno_ice"]]
            topn_ice$num_of_chemicals <- lapply(topn_ice$assay_num, function(x) paste0(x, "/", num_chemicals))
            topn_ice$assay_num <- lapply(topn_ice$assay_num, function(x) x/num_chemicals)
            colnames(topn_ice) <- c("annotation", "percentage", "num_of_chemicals")
        }
        
        if(nrow(tmp[["anno_invitrodb"]]) > 0){
            topn_invitrodb <- tmp[["anno_invitrodb"]][, c("assay_endpoint_name", "assay_num")]
            topn_invitrodb$num_of_chemicals <- lapply(topn_invitrodb$assay_num, function(x) paste0(x, "/", num_chemicals))
            topn_invitrodb$assay_num <- lapply(topn_invitrodb$assay_num, function(x) x/num_chemicals)
            colnames(topn_invitrodb) <- c("annotation", "percentage", "num_of_chemicals")
        }
        
        if(nrow(tmp[["anno_leadscope"]]) > 0){
            topn_leadscope <- tmp[["anno_leadscope"]][, c("model_name", "id_count", "id_mean")]
            topn_leadscope$num_of_chemicals <- lapply(topn_leadscope$id_count, function(x) paste0(x, "/", num_chemicals))
            topn_leadscope$id_count <- lapply(topn_leadscope$id_count, function(x) x/num_chemicals)
            colnames(topn_leadscope) <- c("annotation", "percentage", "num_of_chemicals")
        }
        
        if(nrow(tmp[["anno_ochem"]]) > 0){
            topn_ochem <- tmp[["anno_ochem"]][, c("alert_name", "alert_num")]
            topn_ochem_split <- split(topn_ochem, topn_ochem$alert_name)
            topn_ochem_split_num_of_chemicals <- unlist(lapply(topn_ochem_split, function(x) paste0("Avg. ", sum(x$alert_num)/nrow(x), "/", num_chemicals)))
            topn_ochem_split <- unlist(lapply(topn_ochem_split, function(x) sum(unlist(lapply(x$alert_num, function(y) y/num_chemicals)))/nrow(x)))
            topn_ochem <- data.frame(annotation=names(topn_ochem_split), percentage=topn_ochem_split, num_of_chemicals=topn_ochem_split_num_of_chemicals, stringsAsFactors=FALSE)
        }
        
        # TODO: handle opera differently for now
        if(nrow(tmp[["anno_opera_adme"]]) > 0){
            topn_opera_adme <- tmp[["anno_opera_adme"]]
        }
        if(nrow(tmp[["anno_opera_environmental_fate"]]) > 0){
            topn_opera_environmental_fate <- tmp[["anno_opera_environmental_fate"]]
        }
        if(nrow(tmp[["anno_opera_physicochem"]]) > 0){
            topn_opera_physicochem <- tmp[["anno_opera_physicochem"]]
        }
        if(nrow(tmp[["anno_opera_structural_properties"]]) > 0){
            topn_opera_structural_properties <- tmp[["anno_opera_structural_properties"]]
        }
        if(nrow(tmp[["anno_opera_toxicology_endpoints"]]) > 0){
            topn_opera_toxicology_endpoints <- tmp[["anno_opera_toxicology_endpoints"]]
        }
        
        if(nrow(tmp[["anno_saagar"]]) > 0){
            topn_saagar <- tmp[["anno_saagar"]][, c("alert_name", "descript_num")]
            topn_saagar_split <- split(topn_saagar, topn_saagar$alert_name)
            topn_saagar_split_num_of_chemicals <- unlist(lapply(topn_saagar_split, function(x) paste0("Avg. ", sum(x$descript_num)/nrow(x), "/", num_chemicals)))
            topn_saagar_split <- unlist(lapply(topn_saagar_split, function(x) sum(unlist(lapply(x$descript_num, function(y) y/num_chemicals)))/nrow(x)))
            topn_saagar <- data.frame(annotation=names(topn_saagar_split), percentage=topn_saagar_split, num_of_chemicals=topn_saagar_split_num_of_chemicals, stringsAsFactors=FALSE)
        }
        
        if(nrow(tmp[["anno_t3db"]]) > 0){
            topn_t3db <- tmp[["anno_t3db"]][, c("annotation_name", "annotation_num")]
            topn_t3db_split <- split(topn_t3db, topn_t3db$annotation_name)
            topn_t3db_split_num_of_chemicals <- unlist(lapply(topn_t3db_split, function(x) paste0("Avg. ", sum(x$annotation_num)/nrow(x), "/", num_chemicals)))
            topn_t3db_split <- unlist(lapply(topn_t3db_split, function(x) sum(unlist(lapply(x$annotation_num, function(y) y/num_chemicals)))/nrow(x)))
            topn_t3db <- data.frame(annotation=names(topn_t3db_split), percentage=topn_t3db_split, num_of_chemicals=topn_t3db_split_num_of_chemicals, stringsAsFactors=FALSE)
        }
        
        if(nrow(tmp[["anno_toxrefdb_neoplastic"]]) > 0){
            topn_toxrefdb_neoplastic <- tmp[["anno_toxrefdb_neoplastic"]][, c("annotation_name", "annotation_num")]
            topn_toxrefdb_neoplastic_split <- split(topn_toxrefdb_neoplastic, topn_toxrefdb_neoplastic$annotation_name)
            topn_toxrefdb_neoplastic_split_num_of_chemicals <- unlist(lapply(topn_toxrefdb_neoplastic_split, function(x) paste0("Avg. ", sum(x$annotation_num)/nrow(x), "/", num_chemicals)))
            topn_toxrefdb_neoplastic_split <- unlist(lapply(topn_toxrefdb_neoplastic_split, function(x) sum(unlist(x$annotation_num))/nrow(x)))
            topn_toxrefdb_neoplastic <- data.frame(annotation=names(topn_toxrefdb_neoplastic_split), percentage=topn_toxrefdb_neoplastic_split, num_of_chemicals=topn_toxrefdb_neoplastic_split_num_of_chemicals, stringsAsFactors=FALSE)
        }
        if(nrow(tmp[["anno_toxrefdb_nonneoplastic"]]) > 0){
            topn_toxrefdb_nonneoplastic <- tmp[["anno_toxrefdb_nonneoplastic"]][, c("annotation_name", "annotation_num")]
            topn_toxrefdb_nonneoplastic_split <- split(topn_toxrefdb_nonneoplastic, topn_toxrefdb_nonneoplastic$annotation_name)
            topn_toxrefdb_nonneoplastic_split_num_of_chemicals <- unlist(lapply(topn_toxrefdb_nonneoplastic_split, function(x) paste0("Avg. ", sum(x$annotation_num)/nrow(x), "/", num_chemicals)))
            topn_toxrefdb_nonneoplastic_split <- unlist(lapply(topn_toxrefdb_nonneoplastic_split, function(x) sum(unlist(x$annotation_num))/nrow(x)))
            topn_toxrefdb_nonneoplastic <- data.frame(annotation=names(topn_toxrefdb_nonneoplastic_split), percentage=topn_toxrefdb_nonneoplastic_split, num_of_chemicals=topn_toxrefdb_nonneoplastic_split_num_of_chemicals, stringsAsFactors=FALSE)
        }
        
        if(nrow(tmp[["anno_pubchem"]]) > 0){
            topn_pubchem <- tmp[["anno_pubchem"]][, c("annotation_name", "annotation_num")]
            topn_pubchem_split <- split(topn_pubchem, topn_pubchem$annotation_name)
            topn_pubchem_split_num_of_chemicals <- unlist(lapply(topn_pubchem_split, function(x) paste0("Avg. ", sum(x$annotation_num)/nrow(x), "/", num_chemicals)))
            topn_pubchem_split <- unlist(lapply(topn_pubchem_split, function(x) sum(unlist(x$annotation_num))/nrow(x)))
            topn_pubchem <- data.frame(annotation=names(topn_pubchem_split), percentage=topn_pubchem_split, num_of_chemicals=topn_pubchem_split_num_of_chemicals, stringsAsFactors=FALSE)
        }
        
        if(nrow(tmp[["anno_pubchem_props"]]) > 0){
            topn_pubchem_props <- tmp[["anno_pubchem_props"]][, c("annotation_name", "annotation_num")]
            topn_pubchem_props_split <- split(topn_pubchem_props, topn_pubchem_props$annotation_name)
            topn_pubchem_props_split_num_of_chemicals <- unlist(lapply(topn_pubchem_props_split, function(x) paste0("Avg. ", sum(x$annotation_num)/nrow(x), "/", num_chemicals)))
            topn_pubchem_props_split <- unlist(lapply(topn_pubchem_props_split, function(x) sum(unlist(x$annotation_num))/nrow(x)))
            topn_pubchem_props <- data.frame(annotation=names(topn_pubchem_props_split), percentage=topn_pubchem_props_split, num_of_chemicals=topn_pubchem_props_split_num_of_chemicals, stringsAsFactors=FALSE)
        }
        
        if(nrow(tmp[["anno_gras"]]) > 0){
            topn_gras <- tmp[["anno_gras"]][, c("annotation_name", "annotation_num")]
            topn_gras_split <- split(topn_gras, topn_gras$annotation_name)
            topn_gras_split_num_of_chemicals <- unlist(lapply(topn_gras_split, function(x) paste0("Avg. ", sum(x$annotation_num)/nrow(x), "/", num_chemicals)))
            topn_gras_split <- unlist(lapply(topn_gras_split, function(x) sum(unlist(x$annotation_num))/nrow(x)))
            topn_gras <- data.frame(annotation=names(topn_gras_split), percentage=topn_gras_split, num_of_chemicals=topn_gras_split_num_of_chemicals, stringsAsFactors=FALSE)
        }
        
        if(nrow(tmp[["anno_foodb_enzymes"]]) > 0){
            topn_foodb_enzymes <- tmp[["anno_foodb_enzymes"]][, c("annotation_name", "annotation_num")]
            topn_foodb_enzymes_split <- split(topn_foodb_enzymes, topn_foodb_enzymes$annotation_name)
            topn_foodb_enzymes_split_num_of_chemicals <- unlist(lapply(topn_foodb_enzymes_split, function(x) paste0("Avg. ", sum(x$annotation_num)/nrow(x), "/", num_chemicals)))
            topn_foodb_enzymes_split <- unlist(lapply(topn_foodb_enzymes_split, function(x) sum(unlist(x$annotation_num))/nrow(x)))
            topn_foodb_enzymes <- data.frame(annotation=names(topn_foodb_enzymes_split), percentage=topn_foodb_enzymes_split, num_of_chemicals=topn_foodb_enzymes_split_num_of_chemicals, stringsAsFactors=FALSE)
        }
        if(nrow(tmp[["anno_foodb_flavors"]]) > 0){
            topn_foodb_flavors <- tmp[["anno_foodb_flavors"]][, c("annotation_name", "annotation_num")]
            topn_foodb_flavors_split <- split(topn_foodb_flavors, topn_foodb_flavors$annotation_name)
            topn_foodb_flavors_split_num_of_chemicals <- unlist(lapply(topn_foodb_flavors_split, function(x) paste0("Avg. ", sum(x$annotation_num)/nrow(x), "/", num_chemicals)))
            topn_foodb_flavors_split <- unlist(lapply(topn_foodb_flavors_split, function(x) sum(unlist(x$annotation_num))/nrow(x)))
            topn_foodb_flavors <- data.frame(annotation=names(topn_foodb_flavors_split), percentage=topn_foodb_flavors_split, num_of_chemicals=topn_foodb_flavors_split_num_of_chemicals, stringsAsFactors=FALSE)
        }
        if(nrow(tmp[["anno_foodb_foodcontent"]]) > 0){
            topn_foodb_foodcontent <- tmp[["anno_foodb_foodcontent"]][, c("annotation_name", "annotation_num")]
            topn_foodb_foodcontent_split <- split(topn_foodb_foodcontent, topn_foodb_foodcontent$annotation_name)
            topn_foodb_foodcontent_split_num_of_chemicals <- unlist(lapply(topn_foodb_foodcontent_split, function(x) paste0("Avg. ", sum(x$annotation_num)/nrow(x), "/", num_chemicals)))
            topn_foodb_foodcontent_split <- unlist(lapply(topn_foodb_foodcontent_split, function(x) sum(unlist(x$annotation_num))/nrow(x)))
            topn_foodb_foodcontent <- data.frame(annotation=names(topn_foodb_foodcontent_split), percentage=topn_foodb_foodcontent_split, num_of_chemicals=topn_foodb_foodcontent_split_num_of_chemicals, stringsAsFactors=FALSE)
        }
        if(nrow(tmp[["anno_foodb_healtheffects"]]) > 0){
            topn_foodb_healtheffects <- tmp[["anno_foodb_healtheffects"]][, c("annotation_name", "annotation_num")]
            topn_foodb_healtheffects_split <- split(topn_foodb_healtheffects, topn_foodb_healtheffects$annotation_name)
            topn_foodb_healtheffects_split_num_of_chemicals <- unlist(lapply(topn_foodb_healtheffects_split, function(x) paste0("Avg. ", sum(x$annotation_num)/nrow(x), "/", num_chemicals)))
            topn_foodb_healtheffects_split <- unlist(lapply(topn_foodb_healtheffects_split, function(x) sum(unlist(x$annotation_num))/nrow(x)))
            topn_foodb_healtheffects <- data.frame(annotation=names(topn_foodb_healtheffects_split), percentage=topn_foodb_healtheffects_split, num_of_chemicals=topn_foodb_healtheffects_split_num_of_chemicals, stringsAsFactors=FALSE)
        }
        if(nrow(tmp[["anno_foodb_ontology"]]) > 0){
            topn_foodb_ontology <- tmp[["anno_foodb_ontology"]][, c("annotation_name", "annotation_num")]
            topn_foodb_ontology_split <- split(topn_foodb_ontology, topn_foodb_ontology$annotation_name)
            topn_foodb_ontology_split_num_of_chemicals <- unlist(lapply(topn_foodb_ontology_split, function(x) paste0("Avg. ", sum(x$annotation_num)/nrow(x), "/", num_chemicals)))
            topn_foodb_ontology_split <- unlist(lapply(topn_foodb_ontology_split, function(x) sum(unlist(x$annotation_num))/nrow(x)))
            topn_foodb_ontology <- data.frame(annotation=names(topn_foodb_ontology_split), percentage=topn_foodb_ontology_split, num_of_chemicals=topn_foodb_ontology_split_num_of_chemicals, stringsAsFactors=FALSE)
        }
        
        if(nrow(tmp[["anno_epa_props"]]) > 0){
            
            topn_epa_props <- tmp[["anno_epa_props"]][, c("annotation_name", "annotation_num")]
            topn_epa_props_split <- split(topn_epa_props, topn_epa_props$annotation_name)
            topn_epa_props_split_num_of_chemicals <- unlist(lapply(topn_epa_props_split, function(x) paste0("Avg. ", sum(x$annotation_num)/nrow(x), "/", num_chemicals)))
            topn_epa_props_split <- unlist(lapply(topn_epa_props_split, function(x) sum(unlist(x$annotation_num))/nrow(x)))
            topn_epa_props <- data.frame(annotation=names(topn_epa_props_split), percentage=topn_epa_props_split, num_of_chemicals=topn_epa_props_split_num_of_chemicals, stringsAsFactors=FALSE)
        }
        
        annotation_datasets <- list(
            "plotout_admet_binr"=topn_admet_binr,
            "plotout_chembl"=topn_chembl,
            "plotout_cpd"=topn_cpd,
            "plotout_ctd_bioprocess"=topn_ctd_bioprocess,
            "plotout_ctd_cellcomp"=topn_ctd_cellcomp,
            "plotout_ctd_diseases"=topn_ctd_diseases,
            "plotout_ctd_genes"=topn_ctd_diseases,
            "plotout_ctd_molfunct"=topn_ctd_molfunct,
            "plotout_ctd_phenotypes"=topn_ctd_phenotypes,
            "plotout_drugbank_atccodes"=topn_drugbank_atccodes,
            "plotout_drugbank_carriers"=topn_drugbank_carriers,
            "plotout_drugbank_enzymes"=topn_drugbank_enzymes,
            "plotout_drugbank_targets"=topn_drugbank_targets,
            "plotout_drugbank_transporters"=topn_drugbank_transporters,
            "plotout_hmdb_biospecimenlocations"=topn_hmdb_biospecimenlocations,
            "plotout_hmdb_cellularlocations"=topn_hmdb_cellularlocations,
            "plotout_hmdb_diseases"=topn_hmdb_diseases,
            "plotout_hmdb_genes"=topn_hmdb_genes,
            "plotout_hmdb_tissuelocations"=topn_hmdb_tissuelocations,
            "plotout_ice"=topn_ice,
            "plotout_invitrodb"=topn_invitrodb,
            "plotout_leadscope"=topn_leadscope,
            "plotout_ochem"=topn_ochem,
            "plotout_t3db"=topn_t3db,
            "plotout_saagar"=topn_saagar,
            "plotout_toxrefdb_neoplastic"=topn_toxrefdb_neoplastic,
            "plotout_toxrefdb_nonneoplastic"=topn_toxrefdb_nonneoplastic,
            "plotout_pubchem"=topn_pubchem,
            "plotout_pubchem_props"=topn_pubchem_props,
            "plotout_gras"=topn_gras,
            "plotout_foodb_enzymes"=topn_foodb_enzymes,
            "plotout_foodb_flavors"=topn_foodb_flavors,
            "plotout_foodb_foodcontent"=topn_foodb_foodcontent,
            "plotout_foodb_healtheffects"=topn_foodb_healtheffects,
            "plotout_foodb_ontology"=topn_foodb_ontology,
            "plotout_epa_props"=topn_epa_props
        )
        
        final_plots <- lapply(names(annotation_datasets), function(dataset) {
            
            if(dataset %like% "opera"){
                # Handle Opera differently for now:
                if(nrow(topn_opera_adme) > 0){
                    colnames(topn_opera_adme) <- c("Chemical Name", "DTXSID", "Model", "Value", "Category", "Description")
                    return(topn_opera_adme)
                } else {
                    return(data.frame(no_data=c("")))
                }
                
                if(nrow(topn_opera_environmental_fate) > 0){
                    colnames(topn_opera_environmental_fate) <- c("Chemical Name", "DTXSID", "Model", "Value", "Category", "Description")
                    return(topn_opera_environmental_fate)
                } else {
                    return(data.frame(no_data=c("")))
                }
                
                if(nrow(topn_opera_physicochem) > 0){
                    colnames(topn_opera_physicochem) <- c("Chemical Name", "DTXSID", "Model", "Value", "Category", "Description")
                    return(topn_opera_physicochem)
                } else {
                    return(data.frame(no_data=c("")))
                }
                
                if(nrow(topn_opera_structural_properties) > 0){
                    colnames(topn_opera_structural_properties) <- c("Chemical Name", "DTXSID", "Model", "Value", "Category", "Description")
                    return(topn_opera_structural_properties)
                } else {
                    return(data.frame(no_data=c("")))
                }
                
                if(nrow(topn_opera_toxicology_endpoints) > 0){
                    colnames(topn_opera_toxicology_endpoints) <- c("Chemical Name", "DTXSID", "Model", "Value", "Category", "Description")
                    return(topn_opera_toxicology_endpoints)
                } else {
                    return(data.frame(no_data=c("")))
                }
            }
            
            if(nrow(annotation_datasets[[dataset]]) > 0){
                return(generate_plotly(
                    df=annotation_datasets[[dataset]],
                    title=dataset,
                    mode="annotations"
                ))
            } else {
                return(generate_plotly(
                    df=NULL,
                    title=NULL,
                    mode="blank"
                ))
            }
        })
        
        names(final_plots) <- names(annotation_datasets)
    } else if(mode == "chemical"){
        
        tmp <- chemical_annotations(selected=as.integer(epa_ids$epa_id), input=input, selected_node=selected_node, dtxsids=dtxsids, switch_mode=type)
        
        topn_admet_binr <- tmp[['anno_admet_binr']]
        topn_admet_catg <- tmp[['anno_admet_catg']]
        topn_admet_cont <- tmp[['anno_admet_cont']]
        
        topn_chembl <- data.frame(matrix(ncol=2, nrow=0))
        topn_cpd <- data.frame(matrix(ncol=2, nrow=0))
        topn_ctd_bioprocess <- data.frame(matrix(ncol=2, nrow=0))
        topn_ctd_cellcomp <- data.frame(matrix(ncol=2, nrow=0))
        topn_ctd_diseases <- data.frame(matrix(ncol=2, nrow=0))
        topn_ctd_genes <- data.frame(matrix(ncol=2, nrow=0))
        topn_ctd_molfunct <- data.frame(matrix(ncol=2, nrow=0))
        topn_ctd_phenotypes <- data.frame(matrix(ncol=2, nrow=0))
        topn_drugbank_atccodes <- data.frame(matrix(ncol=2, nrow=0))
        topn_drugbank_carriers <- data.frame(matrix(ncol=2, nrow=0))
        topn_drugbank_enzymes <- data.frame(matrix(ncol=2, nrow=0))
        topn_drugbank_targets <- data.frame(matrix(ncol=2, nrow=0))
        topn_drugbank_transporters <- data.frame(matrix(ncol=2, nrow=0))
        topn_hmdb_biospecimenlocations <- data.frame(matrix(ncol=2, nrow=0))
        topn_hmdb_cellularlocations <- data.frame(matrix(ncol=2, nrow=0))
        topn_hmdb_diseases <- data.frame(matrix(ncol=2, nrow=0))
        topn_hmdb_genes <- data.frame(matrix(ncol=2, nrow=0))
        topn_hmdb_tissuelocations <- data.frame(matrix(ncol=2, nrow=0))
        topn_ice <- data.frame(matrix(ncol=2, nrow=0))
        topn_invitrodb <- data.frame(matrix(ncol=2, nrow=0))
        topn_leadscope <- data.frame(matrix(ncol=2, nrow=0))
        topn_ochem <- data.frame(matrix(ncol=2, nrow=0))
        topn_opera_adme <- data.frame(matrix(ncol=6, nrow=0))
        topn_opera_environmental_fate <- data.frame(matrix(ncol=6, nrow=0))
        topn_opera_physicochem <- data.frame(matrix(ncol=6, nrow=0))
        topn_opera_structural_properties <- data.frame(matrix(ncol=6, nrow=0))
        topn_opera_toxicology_endpoints <- data.frame(matrix(ncol=6, nrow=0))
        topn_saagar <- data.frame(matrix(ncol=2, nrow=0))
        topn_t3db <- data.frame(matrix(ncol=2, nrow=0))
        topn_toxrefdb_neoplastic <- data.frame(matrix(ncol=2, nrow=0))
        topn_toxrefdb_nonneoplastic <- data.frame(matrix(ncol=2, nrow=0))
        topn_pubchem <- data.frame(matrix(ncol=2, nrow=0))
        topn_pubchem_props <- data.frame(matrix(ncol=2, nrow=0))
        topn_gras <- data.frame(matrix(ncol=2, nrow=0))
        topn_foodb_enzymes <- data.frame(matrix(ncol=2, nrow=0))
        topn_foodb_flavors <- data.frame(matrix(ncol=2, nrow=0))
        topn_foodb_foodcontent <- data.frame(matrix(ncol=2, nrow=0))
        topn_foodb_healtheffects <- data.frame(matrix(ncol=2, nrow=0))
        topn_foodb_ontology <- data.frame(matrix(ncol=2, nrow=0))
        topn_epa_props <- data.frame(matrix(ncol=2, nrow=0))
        
        # if(nrow(tmp[["anno_admet_binr"]]) > 0){
        #     topn_admet_binr <- tmp[["anno_admet_binr"]]
        # }
        if(nrow(tmp[["anno_chembl"]]) > 0){
            topn_chembl <- tmp[["anno_chembl"]]
        }
        if(nrow(tmp[["anno_cpd"]]) > 0){
            topn_cpd <- tmp[["anno_cpd"]]
        }
        if(nrow(tmp[["anno_ctd_bioprocess"]]) > 0){
            topn_ctd_bioprocess <- tmp[["anno_ctd_bioprocess"]]
        }
        if(nrow(tmp[["anno_ctd_cellcomp"]]) > 0){
            topn_ctd_cellcomp <- tmp[["anno_ctd_cellcomp"]]
        }
        if(nrow(tmp[["anno_ctd_diseases"]]) > 0){
            topn_ctd_diseases <- tmp[["anno_ctd_diseases"]]
        }
        if(nrow(tmp[["anno_ctd_genes"]]) > 0){
            topn_ctd_genes <- tmp[["anno_ctd_genes"]]
        }
        if(nrow(tmp[["anno_ctd_molfunct"]]) > 0){
            topn_ctd_molfunct <- tmp[["anno_ctd_molfunct"]]
        }
        if(nrow(tmp[["anno_ctd_phenotypes"]]) > 0){
            topn_ctd_phenotypes <- tmp[["anno_ctd_phenotypes"]]
        }
        if(nrow(tmp[["anno_drugbank_atccodes"]]) > 0){
            topn_drugbank_atccodes <- tmp[["anno_drugbank_atccodes"]]
        }
        if(nrow(tmp[["anno_drugbank_carriers"]]) > 0){
            topn_drugbank_carriers <- tmp[["anno_drugbank_carriers"]]
        }
        if(nrow(tmp[["anno_drugbank_enzymes"]]) > 0){
            topn_drugbank_enzymes <- tmp[["anno_drugbank_enzymes"]]
        }
        if(nrow(tmp[["anno_drugbank_targets"]]) > 0){
            topn_drugbank_targets <- tmp[["anno_drugbank_targets"]]
        }
        if(nrow(tmp[["anno_drugbank_transporters"]]) > 0){
            topn_drugbank_transporters <- tmp[["anno_drugbank_transporters"]]
        }
        if(nrow(tmp[["anno_hmdb_biospecimenlocations"]]) > 0){
            topn_hmdb_biospecimenlocations <- tmp[["anno_hmdb_biospecimenlocations"]]
        }
        if(nrow(tmp[["anno_hmdb_cellularlocations"]]) > 0){
            topn_hmdb_cellularlocations <- tmp[["anno_hmdb_cellularlocations"]]
        }
        if(nrow(tmp[["anno_hmdb_diseases"]]) > 0){
            topn_hmdb_diseases <- tmp[["anno_hmdb_diseases"]]
        }
        if(nrow(tmp[["anno_hmdb_genes"]]) > 0){
            topn_hmdb_genes <- tmp[["anno_hmdb_genes"]]
        }
        if(nrow(tmp[["anno_hmdb_tissuelocations"]]) > 0){
            topn_hmdb_tissuelocations <- tmp[["anno_hmdb_tissuelocations"]]
        }
        if(nrow(tmp[["anno_ice"]]) > 0){
            topn_ice <- tmp[["anno_ice"]]
        }
        if(nrow(tmp[["anno_invitrodb"]]) > 0){
            topn_invitrodb <- tmp[["anno_invitrodb"]]
        }
        if(nrow(tmp[["anno_leadscope"]]) > 0){
            topn_leadscope <- tmp[["anno_leadscope"]]
        }
        if(nrow(tmp[["anno_ochem"]]) > 0){
            topn_ochem <- tmp[["anno_ochem"]]
        }
        if(nrow(tmp[["anno_opera_adme"]]) > 0){
            topn_opera_adme <- tmp[["anno_opera_adme"]]
        }
        if(nrow(tmp[["anno_opera_environmental_fate"]]) > 0){
            topn_opera_environmental_fate <- tmp[["anno_opera_environmental_fate"]]
        }
        if(nrow(tmp[["anno_opera_physicochem"]]) > 0){
            topn_opera_physicochem <- tmp[["anno_opera_physicochem"]]
        }
        if(nrow(tmp[["anno_opera_structural_properties"]]) > 0){
            topn_opera_structural_properties <- tmp[["anno_opera_structural_properties"]]
        }
        if(nrow(tmp[["anno_opera_toxicology_endpoints"]]) > 0){
            topn_opera_toxicology_endpoints <- tmp[["anno_opera_toxicology_endpoints"]]
        }
        if(nrow(tmp[["anno_saagar"]]) > 0){
            topn_saagar <- tmp[["anno_saagar"]]
        }
        if(nrow(tmp[["anno_t3db"]]) > 0){
            topn_t3db <- tmp[["anno_t3db"]]
        }
        if(nrow(tmp[["anno_toxrefdb_neoplastic"]]) > 0){
            topn_toxrefdb_neoplastic <- tmp[["anno_toxrefdb_neoplastic"]]
        }
        if(nrow(tmp[["anno_toxrefdb_nonneoplastic"]]) > 0){
            topn_toxrefdb_nonneoplastic <- tmp[["anno_toxrefdb_nonneoplastic"]]
        }
        
        if(nrow(tmp[["anno_pubchem"]]) > 0){
            topn_pubchem <- tmp[["anno_pubchem"]]
        }
        
        if(nrow(tmp[["anno_pubchem_props"]]) > 0){
            topn_pubchem_props <- tmp[["anno_pubchem_props"]]
        }
        
        if(nrow(tmp[["anno_gras"]]) > 0){
            topn_gras <- tmp[["anno_gras"]]
        }
        
        if(nrow(tmp[["anno_foodb_enzymes"]]) > 0){
            topn_foodb_enzymes <- tmp[["anno_foodb_enzymes"]]
        }
        if(nrow(tmp[["anno_foodb_flavors"]]) > 0){
            topn_foodb_flavors <- tmp[["anno_foodb_flavors"]]
        }
        if(nrow(tmp[["anno_foodb_foodcontent"]]) > 0){
            topn_foodb_foodcontent <- tmp[["anno_foodb_foodcontent"]]
        }
        if(nrow(tmp[["anno_foodb_healtheffects"]]) > 0){
            topn_foodb_healtheffects <- tmp[["anno_foodb_healtheffects"]]
        }
        if(nrow(tmp[["anno_foodb_ontology"]]) > 0){
            topn_foodb_ontology <- tmp[["anno_foodb_ontology"]]
        }
        
        if(nrow(tmp[["anno_epa_props"]]) > 0){
            topn_epa_props <- tmp[["anno_epa_props"]]
        }
        

        # if(nrow(topn_admet_binr) > 0){
        #     colnames(topn_admet_binr) <- c("Chemical Name", "CASRN", "DTXSID", "Model", "Interpretation", "Description")
        # }
        if(nrow(topn_chembl) > 0){
            # topn_chembl_df <- as.data.frame(topn_chembl)
            # alerts <- setdiff(colnames(topn_chembl_df), c("preferred_name", "casrn", "dsstox_substance_id"))
            # topn_chembl_df[, alerts] <- sapply(topn_chembl_df[, alerts], as.numeric)
            # topn_chembl_mat <- data.matrix(topn_chembl_df)
            # rownames(topn_chembl_mat) <- unlist(lapply(topn_chembl_df$preferred_name, function(x){
            #     name <- x
            #     if(nchar(name) > 30){
            #         name <- paste0(substr(name, 1, 30), "...")
            #     }
            #     return(name)
            # }))
            # plot_chembl <- plot_ly(x=colnames(topn_chembl_mat[, alerts]), y=rownames(topn_chembl_mat), z=topn_chembl_mat[, alerts], type="heatmap") %>%
            #     layout(margin=list(l=120))
            # 
            # output$dt_chembl_plot <- renderPlotly(plot_chembl)
            colnames(topn_chembl) <- c("Chemical Name", "CASRN", "DTXSID", colnames(topn_chembl)[seq(4, ncol(topn_chembl))])
        }

        if(nrow(topn_cpd) > 0){
            colnames(topn_cpd) <- c("Chemical Name", "CASRN", "DTXSID", "General Category", "Product Family", "Product Type")
        }

        if(nrow(topn_ctd_bioprocess) > 0){
            colnames(topn_ctd_bioprocess) <- c("Chemical Name", "CASRN", "DTXSID", "GO Term", "Corrected p-value")
        }

        if(nrow(topn_ctd_cellcomp) > 0){
            colnames(topn_ctd_cellcomp) <- c("Chemical Name", "CASRN", "DTXSID", "GO Term", "Corrected p-value")
        }

        if(nrow(topn_ctd_diseases) > 0){
            colnames(topn_ctd_diseases) <- c("Chemical Name", "CASRN", "DTXSID", "Disease Name", "Disease ID", "Direct Evidence", "Inference Gene Symbol", "Inference Score", "OMIM IDs", "PubMed IDs")
        }

        if(nrow(topn_ctd_genes) > 0){
            colnames(topn_ctd_genes) <- c("Chemical Name", "CASRN", "DTXSID", "Gene Symbol", "Gene ID", "Gene Forms", "Organism", "Interaction", "Interaction Actions", "PubMed IDs")
        }

        if(nrow(topn_ctd_molfunct) > 0){
            colnames(topn_ctd_molfunct) <- c("Chemical Name", "CASRN", "DTXSID", "GO Term", "Corrected p-value")
        }

        if(nrow(topn_ctd_phenotypes) > 0){
            colnames(topn_ctd_phenotypes) <- c("Chemical Name", "CASRN", "DTXSID", "Phenotype Name", "Phenotype ID", "Comentioned Terms", "Organism", "Interaction", "Interaction Actions", "Anatomy Terms", "Inference Gene Symbols", "PubMed IDs")
        }

        if(nrow(topn_drugbank_atccodes) > 0){
            colnames(topn_drugbank_atccodes) <- c("Chemical Name", "CASRN", "DTXSID", "ATC Code", "ATC Annotation", "ATC Level")
        }

        if(nrow(topn_drugbank_carriers) > 0){
            colnames(topn_drugbank_carriers) <- c("Chemical Name", "CASRN", "DTXSID", "Carrier")
        }

        if(nrow(topn_drugbank_enzymes) > 0){
            colnames(topn_drugbank_enzymes) <- c("Chemical Name", "CASRN", "DTXSID", "Enzyme")
        }

        if(nrow(topn_drugbank_targets) > 0){
            colnames(topn_drugbank_targets) <- c("Chemical Name", "CASRN", "DTXSID", "Target")
        }

        if(nrow(topn_drugbank_transporters) > 0){
            colnames(topn_drugbank_transporters) <- c("Chemical Name", "CASRN", "DTXSID", "Transporter")
        }

        if(nrow(topn_hmdb_biospecimenlocations) > 0){
            colnames(topn_hmdb_biospecimenlocations) <- c("Chemical Name", "CASRN", "DTXSID", "Biospecimen Location")
        }

        if(nrow(topn_hmdb_cellularlocations) > 0){
            colnames(topn_hmdb_cellularlocations) <- c("Chemical Name", "CASRN", "DTXSID", "Cellular Location")
        }

        if(nrow(topn_hmdb_diseases) > 0){
            colnames(topn_hmdb_diseases) <- c("Chemical Name", "CASRN", "DTXSID", "Disease Name")
        }

        if(nrow(topn_hmdb_genes) > 0){
            colnames(topn_hmdb_genes) <- c("Chemical Name", "CASRN", "DTXSID", "Gene")
        }

        if(nrow(topn_hmdb_tissuelocations) > 0){
            colnames(topn_hmdb_tissuelocations) <- c("Chemical Name", "CASRN", "DTXSID", "Location")
        }

        if(nrow(topn_ice) > 0){
            colnames(topn_ice) <- c("Assay", "Endpoint", "Substance Type", "CASRN", "QSAR Ready ID", "Value", "Unit", "Species", "Receptor Species", "Route", "Sex", "Strain", "Life Stage", "Tissue", "Lesion", "Location", "Assay Source", "In Vitro Assay Format", "Reference", "Reference URL", "DTXSID", "Chemical Name", "PubMed ID")
        }

        if(nrow(topn_invitrodb) > 0){
            colnames(topn_invitrodb) <- c("Chemical Name", "CASRN", "DTXSID", "Assay Name", "Assay Endpoint", "Endpoint Value")
        }

        if(nrow(topn_leadscope) > 0){
            colnames(topn_leadscope) <- c("Chemical Name", "CASRN", "DTXSID", "Model", "Interpretation", "Description")
        }

        if(nrow(topn_ochem) > 0){
            colnames(topn_ochem) <- c("Chemical Name", "CASRN", "DTXSID", "Alert", "Value")
        }

        if(nrow(topn_opera_adme) > 0){
            colnames(topn_opera_adme) <- c("Chemical Name", "CASRN", "DTXSID", "Model", "Value", "Model Category", "Description")
        }

        if(nrow(topn_opera_environmental_fate) > 0){
            colnames(topn_opera_environmental_fate) <- c("Chemical Name", "CASRN", "DTXSID", "Model", "Value", "Model Category", "Description")
        }

        if(nrow(topn_opera_physicochem) > 0){
            colnames(topn_opera_physicochem) <- c("Chemical Name", "CASRN", "DTXSID", "Model", "Value", "Model Category", "Description")
        }

        if(nrow(topn_opera_structural_properties) > 0){
            colnames(topn_opera_structural_properties) <- c("Chemical Name", "CASRN", "DTXSID", "Model", "Value", "Model Category", "Description")
        }

        if(nrow(topn_opera_toxicology_endpoints) > 0){
            colnames(topn_opera_toxicology_endpoints) <- c("Chemical Name", "DTXSID", "Model", "Value", "Model Category", "Description")
        }

        if(nrow(topn_saagar) > 0){
            colnames(topn_saagar) <- c("Chemical Name", "CASRN", "DTXSID", "Alert", "Value")
        }

        if(nrow(topn_t3db) > 0){
            colnames(topn_t3db) <- c("Chemical Name", "CASRN", "DTXSID", "Target ID", "Target Name", "UniProt ID")
        }

        if(nrow(topn_toxrefdb_neoplastic) > 0){
            colnames(topn_toxrefdb_neoplastic) <- c("Chemical Name", "CASRN", "DTXSID", "Annotation")
        }

        if(nrow(topn_toxrefdb_nonneoplastic) > 0){
            colnames(topn_toxrefdb_nonneoplastic) <- c("Chemical Name", "CASRN", "DTXSID", "Annotation")
        }
        
        if(nrow(topn_pubchem) > 0){
            colnames(topn_pubchem) <- c("Chemical Name", "CASRN", "DTXSID", "Bioassay Name", "Source Name", "Source ID", "Substance Type", "Outcome Type", "Project Category", "Bioassay Group", "Bioassay Types", "Protein Accessions", "AID", "SID", "SID Group", "Activity Outcome", "Activity Name", "Activity Qualifier", "Activity Value", "Protein Accession", "Gene ID", "Target Tax ID", "PMID")
        }
        
        if(nrow(topn_pubchem_props) > 0){
            colnames(topn_pubchem_props) <- c("Chemical Name", "CASRN", "DTXSID", "PubChem Attribute", "Value") # TODO
        }
        
        if(nrow(topn_gras) > 0){
            colnames(topn_gras) <- c("Chemical Name", "CASRN", "DTXSID", "GRAS Substance", "SCOGS Report Number", "Year", "SCOGS Conclusion", "SCOGS Interpretation", "21 CFR Regulation", "NTIS Accession")
        }
        
        if(nrow(topn_foodb_enzymes) > 0){
            colnames(topn_foodb_enzymes) <- c("Chemical Name", "CASRN", "DTXSID", "Enzyme Name", "Gene Name", "Kingdom", "Superclass", "Class", "Subclass", "Citations")
        }
        if(nrow(topn_foodb_flavors) > 0){
            colnames(topn_foodb_flavors) <- c("Chemical Name", "CASRN", "DTXSID", "Flavor Name", "Flavor Group", "Category", "Kingdom", "Superclass", "Class", "Subclass", "Citations")
        }
        if(nrow(topn_foodb_foodcontent) > 0){
            colnames(topn_foodb_foodcontent) <- c("Chemical Name", "CASRN", "DTXSID", "Name", "Scientific Name", "Food Group", "Food Subgroup", "Food Type", "Description", "Kingdom", "Superclass", "Class", "Subclass")
        }
        if(nrow(topn_foodb_healtheffects) > 0){
            colnames(topn_foodb_healtheffects) <- c("Chemical Name", "CASRN", "DTXSID", "Health Effect Name", "Kingdom", "Superclass", "Class", "Subclass", "Description")
        }
        if(nrow(topn_foodb_ontology) > 0){
            colnames(topn_foodb_ontology) <- c("Chemical Name", "CASRN", "DTXSID", "Term", "Definition", "Level", "Kingdom", "Superclass", "Class", "Subclass")
        }
        if(nrow(topn_epa_props) > 0){
            colnames(topn_epa_props) <- c("Chemical Name", "CASRN", "DTXSID", "Attribute", "Value")
        }
        
        annotation_datasets <- list(
            "dt_admet_binr"=topn_admet_binr,
            "dt_admet_catg"=topn_admet_catg,
            "dt_admet_cont"=topn_admet_cont,
            "dt_chembl"=topn_chembl,
            "dt_cpd"=topn_cpd,
            "dt_ctd_bioprocess"=topn_ctd_bioprocess,
            "dt_ctd_cellcomp"=topn_ctd_cellcomp,
            "dt_ctd_diseases"=topn_ctd_diseases,
            "dt_ctd_genes"=topn_ctd_diseases,
            "dt_ctd_molfunct"=topn_ctd_molfunct,
            "dt_ctd_phenotypes"=topn_ctd_phenotypes,
            "dt_drugbank_atccodes"=topn_drugbank_atccodes,
            "dt_drugbank_carriers"=topn_drugbank_carriers,
            "dt_drugbank_enzymes"=topn_drugbank_enzymes,
            "dt_drugbank_targets"=topn_drugbank_targets,
            "dt_drugbank_transporters"=topn_drugbank_transporters,
            "dt_hmdb_biospecimenlocations"=topn_hmdb_biospecimenlocations,
            "dt_hmdb_cellularlocations"=topn_hmdb_cellularlocations,
            "dt_hmdb_diseases"=topn_hmdb_diseases,
            "dt_hmdb_genes"=topn_hmdb_genes,
            "dt_hmdb_tissuelocations"=topn_hmdb_tissuelocations,
            "dt_ice"=topn_ice,
            "dt_invitrodb"=topn_invitrodb,
            "dt_leadscope"=topn_leadscope,
            "dt_ochem"=topn_ochem,
            "dt_t3db"=topn_t3db,
            "dt_saagar"=topn_saagar,
            "dt_toxrefdb_neoplastic"=topn_toxrefdb_neoplastic,
            "dt_toxrefdb_nonneoplastic"=topn_toxrefdb_nonneoplastic,
            "dt_pubchem"=topn_pubchem,
            "dt_pubchem_props"=topn_pubchem_props,
            "dt_gras"=topn_gras,
            "dt_foodb_enzymes"=topn_foodb_enzymes,
            "dt_foodb_flavors"=topn_foodb_flavors,
            "dt_foodb_foodcontent"=topn_foodb_foodcontent,
            "dt_foodb_healtheffects"=topn_foodb_healtheffects,
            "dt_foodb_ontology"=topn_foodb_ontology,
            "dt_epa_props"=topn_epa_props
        )
        
        final_plots <- annotation_datasets
    }
    return(final_plots)
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

# Load cluster annotations
cluster_annotations <- function(selected, clust_l3, input, reactive_clusters, reactive_selected_cluster, dtxsids, switch_mode=""){
    anno_admet_binr <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_admet_binr) <- c("model_name", "positives")
    anno_chembl <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_chembl) <- c("alert_name", "alert_num")
    anno_cpd <- data.frame(matrix(ncol=2, nrow=0))
    
    anno_ctd_bioprocess <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_ctd_bioprocess) <- c("term_name", "term_num")
    anno_ctd_cellcomp <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_ctd_cellcomp) <- c("term_name", "term_num")
    anno_ctd_diseases <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_ctd_diseases) <- c("term_name", "term_num")
    anno_ctd_genes <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_ctd_genes) <- c("term_name", "term_num")
    anno_ctd_molfunct <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_ctd_molfunct) <- c("term_name", "term_num")
    anno_ctd_phenotypes <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_ctd_phenotypes) <- c("term_name", "term_num")
    
    anno_drugbank_atccodes <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_drugbank_atccodes) <- c("annotation_name", "annotation_num")
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
    colnames(anno_leadscope) <- c("model_name", "id_count", "id_mean")
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
    
    anno_pubchem <- data.frame(matrix(ncol=2, nrow=0))
    anno_pubchem_props <- data.frame(matrix(ncol=2, nrow=0))
    
    
    anno_saagar <- data.frame(matrix(ncol=2, nrow=0))
    colnames(anno_saagar) <- c("descript_name", "descript_num")
    anno_t3db <- data.frame(matrix(ncol=2, nrow=0))
    anno_toxrefdb_neoplastic <- data.frame(matrix(ncol=2, nrow=0))
    anno_toxrefdb_nonneoplastic <- data.frame(matrix(ncol=2, nrow=0))
    
    anno_gras <- data.frame(matrix(ncol=2, nrow=0))
    
    anno_foodb_enzymes <- data.frame(matrix(ncol=2, nrow=0))
    anno_foodb_flavors <- data.frame(matrix(ncol=2, nrow=0))
    anno_foodb_foodcontent <- data.frame(matrix(ncol=2, nrow=0))
    anno_foodb_healtheffects <- data.frame(matrix(ncol=2, nrow=0))
    anno_foodb_ontology <- data.frame(matrix(ncol=2, nrow=0))
    
    anno_epa_props <- data.frame(matrix(ncol=2, nrow=0))
    
    # ADMET
    if(input[[paste0("switch__admet_binr", switch_mode)]] == TRUE){
        anno_admet_binr <- run_query(paste0("
                        SELECT DISTINCT
                            mv.*,
                            apm.description
                        FROM
                            admet_binary_saagar_clusters mv,
                            admet_predictive_models apm
                        WHERE
                            mv.clust_l3 = $1
                        AND mv.amodel_id = apm.amodel_id
                    "), args=list(clust_l3))
    }
    
    # ChEMBL
    if(input[[paste0("switch__chembl", switch_mode)]] == TRUE){
        anno_chembl <- run_query(paste0("
                        SELECT DISTINCT
                            mv.*
                        FROM
                            chembl_structural_alerts_saagar_clusters mv
                        WHERE
                            mv.clust_l3 = $1
                    "), args=list(clust_l3))
    }
    
    # CPD
    if(input[[paste0("switch__cpd", switch_mode)]] == TRUE){
        anno_cpd <- run_query(paste0("
                        SELECT DISTINCT
                            mv.epa_id,
                            mv.prod_type
                        FROM
                            mv_cpd_categories_saagar mv
                        WHERE
                            mv.clust_l3 = $1
                    "), args=list(clust_l3))
        if(nrow(anno_cpd) > 0){
            anno_cpd <- split(anno_cpd, anno_cpd$prod_type)
            anno_cpd <- unlist(lapply(anno_cpd, function(x) length(unique(x$epa_id))))
            anno_cpd <- data.frame(product_name=names(anno_cpd), product_num=anno_cpd)
        }
    }
    
    # CTD
    if(input[[paste0("switch__ctd_bioprocess", switch_mode)]] == TRUE){
        anno_ctd_bioprocess <- run_query(paste0("
                        SELECT DISTINCT
                            mv.epa_id,
                            mv.go_term_name
                        FROM
                            ctd_biological_process_to_base_chemical_saagar_clusters mv
                        WHERE
                            mv.clust_l3 = $1
                    "), args=list(clust_l3))
        if(nrow(anno_ctd_bioprocess) > 0){
            anno_ctd_bioprocess <- split(anno_ctd_bioprocess, anno_ctd_bioprocess$go_term_name)
            anno_ctd_bioprocess <- unlist(lapply(anno_ctd_bioprocess, function(x) length(unique(x$epa_id))))
            anno_ctd_bioprocess <- data.frame(term_name=names(anno_ctd_bioprocess), term_num=anno_ctd_bioprocess)
        }
    }
    if(input[[paste0("switch__ctd_cellcomp", switch_mode)]] == TRUE){
        anno_ctd_cellcomp <- run_query(paste0("
                        SELECT DISTINCT
                            mv.epa_id,
                            mv.go_term_name
                        FROM
                            ctd_cellular_component_to_base_chemical_saagar_clusters mv
                        WHERE
                            mv.clust_l3 = $1
                    "), args=list(clust_l3))
        if(nrow(anno_ctd_cellcomp) > 0){
            anno_ctd_cellcomp <- split(anno_ctd_cellcomp, anno_ctd_cellcomp$go_term_name)
            anno_ctd_cellcomp <- unlist(lapply(anno_ctd_cellcomp, function(x) length(unique(x$epa_id))))
            anno_ctd_cellcomp <- data.frame(term_name=names(anno_ctd_cellcomp), term_num=anno_ctd_cellcomp)
        }
    }
    if(input[[paste0("switch__ctd_diseases", switch_mode)]] == TRUE){
        anno_ctd_diseases <- run_query(paste0("
                        SELECT DISTINCT
                            mv.epa_id,
                            mv.disease_name
                        FROM
                            ctd_diseases_to_base_chemical_saagar_clusters mv
                        WHERE
                            mv.clust_l3 = $1
                    "), args=list(clust_l3))
        if(nrow(anno_ctd_diseases) > 0){
            anno_ctd_diseases <- split(anno_ctd_diseases, anno_ctd_diseases$disease_name)
            anno_ctd_diseases <- unlist(lapply(anno_ctd_diseases, function(x) length(unique(x$epa_id))))
            anno_ctd_diseases <- data.frame(term_name=names(anno_ctd_diseases), term_num=anno_ctd_diseases)
        }
    }
    if(input[[paste0("switch__ctd_genes", switch_mode)]] == TRUE){
        anno_ctd_genes <- run_query(paste0("
                        SELECT DISTINCT
                            mv.epa_id,
                            mv.gene_symbol
                        FROM
                            ctd_genes_to_base_chemical_saagar_clusters mv
                        WHERE
                            mv.clust_l3 = $1
                    "), args=list(clust_l3))
        if(nrow(anno_ctd_genes) > 0){
            anno_ctd_genes <- split(anno_ctd_genes, anno_ctd_genes$gene_symbol)
            anno_ctd_genes <- unlist(lapply(anno_ctd_genes, function(x) length(unique(x$epa_id))))
            anno_ctd_genes <- data.frame(term_name=names(anno_ctd_genes), term_num=anno_ctd_genes)
        }
    }
    if(input[[paste0("switch__ctd_molfunct", switch_mode)]] == TRUE){
        anno_ctd_molfunct <- run_query(paste0("
                        SELECT DISTINCT
                            mv.epa_id,
                            mv.go_term_name
                        FROM
                            ctd_molecular_function_to_base_chemical_saagar_clusters mv
                        WHERE
                            mv.clust_l3 = $1
                    "), args=list(clust_l3))
        if(nrow(anno_ctd_molfunct) > 0){
            anno_ctd_molfunct <- split(anno_ctd_molfunct, anno_ctd_molfunct$go_term_name)
            anno_ctd_molfunct <- unlist(lapply(anno_ctd_molfunct, function(x) length(unique(x$epa_id))))
            anno_ctd_molfunct <- data.frame(term_name=names(anno_ctd_molfunct), term_num=anno_ctd_molfunct)
        }
    }
    if(input[[paste0("switch__ctd_phenotypes", switch_mode)]] == TRUE){
        anno_ctd_phenotypes <- run_query(paste0("
                        SELECT DISTINCT
                            mv.epa_id,
                            mv.phenotype_name
                        FROM
                            ctd_phenotypes_to_base_chemical_saagar_clusters mv
                        WHERE
                            mv.clust_l3 = $1
                    "), args=list(clust_l3))
        if(nrow(anno_ctd_phenotypes) > 0){
            anno_ctd_phenotypes <- split(anno_ctd_phenotypes, anno_ctd_phenotypes$phenotype_name)
            anno_ctd_phenotypes <- unlist(lapply(anno_ctd_phenotypes, function(x) length(unique(x$epa_id))))
            anno_ctd_phenotypes <- data.frame(term_name=names(anno_ctd_phenotypes), term_num=anno_ctd_phenotypes)
        }
    }
    
    # DrugBank
    if(input[[paste0("switch__drugbank_atccodes", switch_mode)]] == TRUE){
        anno_drugbank_atccodes <- run_query(paste0("
                        SELECT DISTINCT
                            mv.epa_id,
                            mv.annotation
                        FROM
                            base_chemical_drugbank_atccodes_saagar_clusters mv
                        WHERE
                            mv.clust_l3 = $1
                    "), args=list(clust_l3))
        if(nrow(anno_drugbank_atccodes) > 0){
            anno_drugbank_atccodes <- split(anno_drugbank_atccodes, anno_drugbank_atccodes$annotation)
            anno_drugbank_atccodes <- unlist(lapply(anno_drugbank_atccodes, function(x) length(unique(x$epa_id))))
            anno_drugbank_atccodes <- data.frame(annotation_name=names(anno_drugbank_atccodes), annotation_num=anno_drugbank_atccodes)
        }
    }
    
    if(input[[paste0("switch__drugbank_carriers", switch_mode)]] == TRUE){
        anno_drugbank_carriers <- run_query(paste0("
                        SELECT DISTINCT
                            mv.epa_id,
                            mv.annotation
                        FROM
                            base_chemical_drugbank_carriers_saagar_clusters mv
                        WHERE
                            mv.clust_l3 = $1
                    "), args=list(clust_l3))
        if(nrow(anno_drugbank_carriers) > 0){
            anno_drugbank_carriers <- split(anno_drugbank_carriers, anno_drugbank_carriers$annotation)
            anno_drugbank_carriers <- unlist(lapply(anno_drugbank_carriers, function(x) length(unique(x$epa_id))))
            anno_drugbank_carriers <- data.frame(annotation_name=names(anno_drugbank_carriers), annotation_num=anno_cpd)
        }
    }
    if(input[[paste0("switch__drugbank_enzymes", switch_mode)]] == TRUE){
        anno_drugbank_enzymes <- run_query(paste0("
                        SELECT DISTINCT
                            mv.epa_id,
                            mv.annotation
                        FROM
                            base_chemical_drugbank_enzymes_saagar_clusters mv
                        WHERE
                            mv.clust_l3 = $1
                    "), args=list(clust_l3))
        if(nrow(anno_drugbank_enzymes) > 0){
            anno_drugbank_enzymes <- split(anno_drugbank_enzymes, anno_drugbank_enzymes$annotation)
            anno_drugbank_enzymes <- unlist(lapply(anno_drugbank_enzymes, function(x) length(unique(x$epa_id))))
            anno_drugbank_enzymes <- data.frame(annotation_name=names(anno_drugbank_enzymes), annotation_num=anno_cpd)
        }
    }
    if(input[[paste0("switch__drugbank_targets", switch_mode)]] == TRUE){
        anno_drugbank_targets <- run_query(paste0("
                        SELECT DISTINCT
                            mv.epa_id,
                            mv.annotation
                        FROM
                            base_chemical_drugbank_targets_saagar_clusters mv
                        WHERE
                            mv.clust_l3 = $1
                    "), args=list(clust_l3))
        if(nrow(anno_drugbank_targets) > 0){
            anno_drugbank_targets <- split(anno_drugbank_targets, anno_drugbank_targets$annotation)
            anno_drugbank_targets <- unlist(lapply(anno_drugbank_targets, function(x) length(unique(x$epa_id))))
            anno_drugbank_targets <- data.frame(annotation_name=names(anno_drugbank_targets), annotation_num=anno_drugbank_targets)
        }
    }
    if(input[[paste0("switch__drugbank_transporters", switch_mode)]] == TRUE){
        anno_drugbank_transporters <- run_query(paste0("
                        SELECT DISTINCT
                            mv.epa_id,
                            mv.annotation
                        FROM
                            base_chemical_drugbank_transporters_saagar_clusters mv
                        WHERE
                            mv.clust_l3 = $1
                    "), args=list(clust_l3))
        if(nrow(anno_drugbank_transporters) > 0){
            anno_drugbank_transporters <- split(anno_drugbank_transporters, anno_drugbank_transporters$annotation)
            anno_drugbank_transporters <- unlist(lapply(anno_drugbank_transporters, function(x) length(unique(x$epa_id))))
            anno_drugbank_transporters <- data.frame(annotation_name=names(anno_drugbank_transporters), annotation_num=anno_cpd)
        }
    }
    
    # HMDB
    if(input[[paste0("switch__hmdb_biospecimenlocations", switch_mode)]] == TRUE){
        anno_hmdb_biospecimenlocations <- run_query(paste0("
                        SELECT DISTINCT
                            mv.annotation,
                            mv.epa_id
                        FROM
                            base_chemical_hmdb_biospecimen_locations_saagar_clusters mv
                        WHERE
                            mv.clust_l3 = $1
                    "), args=list(clust_l3))
        if(nrow(anno_hmdb_biospecimenlocations) > 0){
            anno_hmdb_biospecimenlocations <- split(anno_hmdb_biospecimenlocations, anno_hmdb_biospecimenlocations$annotation)
            anno_hmdb_biospecimenlocations <- unlist(lapply(anno_hmdb_biospecimenlocations, function(x) length(unique(x$epa_id))))
            anno_hmdb_biospecimenlocations <- data.frame(annotation_name=names(anno_hmdb_biospecimenlocations), annotation_num=anno_hmdb_biospecimenlocations)
        }
    }
    
    if(input[[paste0("switch__hmdb_cellularlocations", switch_mode)]] == TRUE){
        anno_hmdb_cellularlocations <- run_query(paste0("
                        SELECT DISTINCT
                            mv.annotation,
                            mv.epa_id
                        FROM
                            base_chemical_hmdb_cellular_locations_saagar_clusters mv
                        WHERE
                            mv.clust_l3 = $1
                    "), args=list(clust_l3))
        if(nrow(anno_hmdb_cellularlocations) > 0){
            anno_hmdb_cellularlocations <- split(anno_hmdb_cellularlocations, anno_hmdb_cellularlocations$annotation)
            anno_hmdb_cellularlocations <- unlist(lapply(anno_hmdb_cellularlocations, function(x) length(unique(x$epa_id))))
            anno_hmdb_cellularlocations <- data.frame(annotation_name=names(anno_hmdb_cellularlocations), annotation_num=anno_hmdb_cellularlocations)
        }
    }
    
    if(input[[paste0("switch__hmdb_diseases", switch_mode)]] == TRUE){
        anno_hmdb_diseases <- run_query(paste0("
                        SELECT DISTINCT
                            mv.annotation,
                            mv.epa_id
                        FROM
                            base_chemical_hmdb_diseases_saagar_clusters mv
                        WHERE
                            mv.clust_l3 = $1
                    "), args=list(clust_l3))
        if(nrow(anno_hmdb_diseases) > 0){
            anno_hmdb_diseases <- split(anno_hmdb_diseases, anno_hmdb_diseases$annotation)
            anno_hmdb_diseases <- unlist(lapply(anno_hmdb_diseases, function(x) length(unique(x$epa_id))))
            anno_hmdb_diseases <- data.frame(annotation_name=names(anno_hmdb_diseases), annotation_num=anno_hmdb_diseases)
        }
    }
    
    if(input[[paste0("switch__hmdb_genes", switch_mode)]] == TRUE){
        anno_hmdb_genes <- run_query(paste0("
                        SELECT DISTINCT
                            mv.annotation,
                            mv.epa_id
                        FROM
                            base_chemical_hmdb_genes_saagar_clusters mv
                        WHERE
                            mv.clust_l3 = $1
                    "), args=list(clust_l3))
        if(nrow(anno_hmdb_genes) > 0){
            anno_hmdb_genes <- split(anno_hmdb_genes, anno_hmdb_genes$annotation)
            anno_hmdb_genes <- unlist(lapply(anno_hmdb_genes, function(x) length(unique(x$epa_id))))
            anno_hmdb_genes <- data.frame(annotation_name=names(anno_hmdb_genes), annotation_num=anno_hmdb_genes)
        }
    }
    
    if(input[[paste0("switch__hmdb_tissuelocations", switch_mode)]] == TRUE){
        anno_hmdb_tissuelocations <- run_query(paste0("
                        SELECT DISTINCT
                            mv.annotation,
                            mv.epa_id
                        FROM
                            base_chemical_hmdb_tissue_locations_saagar_clusters mv
                        WHERE
                            mv.clust_l3 = $1
                    "), args=list(clust_l3))
        if(nrow(anno_hmdb_tissuelocations) > 0){
            anno_hmdb_tissuelocations <- split(anno_hmdb_tissuelocations, anno_hmdb_tissuelocations$annotation)
            anno_hmdb_tissuelocations <- unlist(lapply(anno_hmdb_tissuelocations, function(x) length(unique(x$epa_id))))
            anno_hmdb_tissuelocations <- data.frame(annotation_name=names(anno_hmdb_tissuelocations), annotation_num=anno_hmdb_tissuelocations)
        }
    }
    
    
    
    
    # ICE
    if(input[[paste0("switch__ice", switch_mode)]] == TRUE){
        anno_ice <- query_ice(chemicals=dtxsids)
        anno_ice <- aggregate_ice_assays(anno_ice)
    }
    
    # InVitroDB
    if(input[[paste0("switch__invitrodb", switch_mode)]] == TRUE){
        anno_invitrodb <- run_query(paste0("
                        SELECT DISTINCT
                            mv.assay_endpoint_name,
                            mv.smi_id --this table doesn't track epa_ids
                        FROM
                            invitrodb_saagar_clusters mv
                        WHERE
                            mv.clust_l3 = $1
                        AND mv.hit_call = 1
                    "), args=list(clust_l3))
        if(nrow(anno_invitrodb) > 0){
            anno_invitrodb <- split(anno_invitrodb, anno_invitrodb$assay_endpoint_name)
            anno_invitrodb <- data.frame(assay_endpoint_name=names(anno_invitrodb), assay_num=unlist(lapply(anno_invitrodb, function(x) nrow(x))))
        }
    }
    
    # Leadscope
    if(input[[paste0("switch__leadscope", switch_mode)]] == TRUE){
        anno_leadscope <- run_query(paste0("
                        SELECT DISTINCT
                            mv.model_name,
                            mv.id_count,
                            mv.id_mean
                        FROM
                            base_chemical_saagar_cluster_statistics mv
                        WHERE
                            mv.clust_l3 = $1
                    "), args=list(clust_l3))
    }
    
    # OChem
    if(input[[paste0("switch__ochem", switch_mode)]] == TRUE){
        anno_ochem <- run_query(paste0("
                        SELECT DISTINCT
                            mv.alert_name,
                            mv.alert_num
                        FROM
                            ochem_structural_alerts_saagar_clusters mv
                        WHERE
                            mv.clust_l3 = $1
                    "), args=list(clust_l3))
    }
    
    # Opera
    if(input[[paste0("switch__opera_adme", switch_mode)]] == TRUE){
        anno_opera_adme <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bc.dsstox_substance_id,
                            mv.model_name,
                            mv.model_value,
                            opm.category,
                            opm.model_description
                        FROM
                            base_chemical_opera_adme_saagar_clusters mv,
                            base_chemicals bc,
                            base_chemical_compounds bcc,
                            opera_predictive_models opm
                        WHERE
                            mv.clust_l3 = $1
                        AND bc.epa_id = mv.epa_id
                        AND bc.epa_id = bcc.epa_id
                        AND mv.omodel_id = opm.omodel_id
                    "), args=list(clust_l3))
    }
    
    if(input[[paste0("switch__opera_environmental_fate", switch_mode)]] == TRUE){
        anno_opera_environmental_fate <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bc.dsstox_substance_id,
                            mv.model_name,
                            mv.model_value,
                            opm.category,
                            opm.model_description
                        FROM
                            base_chemical_opera_environmental_fate_saagar_clusters mv,
                            base_chemicals bc,
                            base_chemical_compounds bcc,
                            opera_predictive_models opm
                        WHERE
                            mv.clust_l3 = $1
                        AND bc.epa_id = mv.epa_id
                        AND bc.epa_id = bcc.epa_id
                        AND mv.omodel_id = opm.omodel_id
                    "), args=list(clust_l3))
    }
    
    if(input[[paste0("switch__opera_physicochem", switch_mode)]] == TRUE){
        anno_opera_physicochem <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bc.dsstox_substance_id,
                            mv.model_name,
                            mv.model_value,
                            opm.category,
                            opm.model_description
                        FROM
                            base_chemical_opera_physicochem_saagar_clusters mv,
                            base_chemicals bc,
                            base_chemical_compounds bcc,
                            opera_predictive_models opm
                        WHERE
                            mv.clust_l3 = $1
                        AND bc.epa_id = mv.epa_id
                        AND bc.epa_id = bcc.epa_id
                        AND mv.omodel_id = opm.omodel_id
                    "), args=list(clust_l3))
    }
    
    if(input[[paste0("switch__opera_structural_properties", switch_mode)]] == TRUE){
        anno_opera_structural_properties <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bc.dsstox_substance_id,
                            mv.model_name,
                            mv.model_value,
                            opm.category,
                            opm.model_description
                        FROM
                            base_chemical_opera_structural_properties_saagar_clusters mv,
                            base_chemicals bc,
                            base_chemical_compounds bcc,
                            opera_predictive_models opm
                        WHERE
                            mv.clust_l3 = $1
                        AND bc.epa_id = mv.epa_id
                        AND bc.epa_id = bcc.epa_id
                        AND mv.omodel_id = opm.omodel_id
                    "), args=list(clust_l3))
    }
    
    if(input[[paste0("switch__opera_toxicology_endpoints", switch_mode)]] == TRUE){
        anno_opera_toxicology_endpoints <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bc.dsstox_substance_id,
                            mv.model_name,
                            mv.model_value,
                            opm.category,
                            opm.model_description
                        FROM
                            base_chemical_opera_toxicology_endpoints_saagar_clusters mv,
                            base_chemicals bc,
                            base_chemical_compounds bcc,
                            opera_predictive_models opm
                        WHERE
                            mv.clust_l3 = $1
                        AND bc.epa_id = mv.epa_id
                        AND bc.epa_id = bcc.epa_id
                        AND mv.omodel_id = opm.omodel_id
                    "), args=list(clust_l3))
    }
    
    # Saagar
    if(input[[paste0("switch__saagar", switch_mode)]] == TRUE){
        anno_saagar <- run_query(paste0("
                        SELECT DISTINCT
                            mv.alert_name,
                            mv.descript_num
                        FROM
                            saagar_descriptors_saagar_clusters mv
                        WHERE
                           mv.clust_l3 = $1
                    "), args=list(clust_l3))
    }
    
    # T3DB
    if(input[[paste0("switch__t3db", switch_mode)]] == TRUE){
        anno_t3db <- run_query(paste0("
                        SELECT DISTINCT
                            bcs.epa_id,
                            tct.target_name
                        FROM
                            t3db_chemicals_to_targets tct,
                            t3db_to_base_chemicals tbc,
                            t3db_curated_chemicals tcc,
                            base_chemical_compounds bcc,
                            base_chemicals bc,
                            base_chemical_to_smiles bcs
                        WHERE
                            tbc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
                        AND tbc.epa_id = bcs.epa_id
                        AND tbc.epa_id = bc.epa_id
                        AND tbc.t3db_id = tcc.t3db_id
                        AND tcc.t3_id = tct.t3_id
                    "), args=as.list(selected))
        
        if(nrow(anno_t3db) > 0){
            anno_t3db <- split(anno_t3db, anno_t3db$target_name)
            anno_t3db <- unlist(lapply(anno_t3db, function(x) length(unique(x$epa_id))))
            anno_t3db <- data.frame(annotation_name=names(anno_t3db), annotation_num=anno_t3db)
        }
    }
        
    # ToxRefDB
    if(input[[paste0("switch__toxrefdb_neoplastic", switch_mode)]] == TRUE){
        anno_toxrefdb_neoplastic <- run_query(paste0("
                        SELECT DISTINCT
                            mv.epa_id,
                            CONCAT(mv.disease_type, '_',
                            mv.study_type, '_',
                            mv.species, '_',
                            mv.sex, '_',
                            mv.life_stage, '_',
                            mv.endpoint_type, '_',
                            mv.endpoint_target) AS annotation
                        FROM
                            mv_toxrefdb_cancer_saagar mv
                        WHERE
                            mv.clust_l3 = $1
                    "), args=list(clust_l3))
        if(nrow(anno_toxrefdb_neoplastic) > 0){
            anno_toxrefdb_neoplastic <- split(anno_toxrefdb_neoplastic, anno_toxrefdb_neoplastic$annotation)
            anno_toxrefdb_neoplastic <- unlist(lapply(anno_toxrefdb_neoplastic, function(x) length(unique(x$epa_id))))
            anno_toxrefdb_neoplastic <- data.frame(annotation_name=names(anno_toxrefdb_neoplastic), annotation_num=anno_toxrefdb_neoplastic)
        }
    }
    
    if(input[[paste0("switch__toxrefdb_nonneoplastic", switch_mode)]] == TRUE){
        anno_toxrefdb_nonneoplastic <- run_query(paste0("
                        SELECT DISTINCT
                            mv.epa_id,
                            CONCAT(mv.disease_type, '_',
                            mv.study_type, '_',
                            mv.species, '_',
                            mv.sex, '_',
                            mv.life_stage, '_',
                            mv.endpoint_type, '_',
                            mv.endpoint_target) AS annotation
                        FROM
                            mv_toxrefdb_nonneoplastic_saagar mv
                        WHERE
                            mv.clust_l3 = $1
                    "), args=list(clust_l3))
        if(nrow(anno_toxrefdb_nonneoplastic) > 0){
            anno_toxrefdb_nonneoplastic <- split(anno_toxrefdb_nonneoplastic, anno_toxrefdb_nonneoplastic$annotation)
            anno_toxrefdb_nonneoplastic <- unlist(lapply(anno_toxrefdb_nonneoplastic, function(x) length(unique(x$epa_id))))
            anno_toxrefdb_nonneoplastic <- data.frame(annotation_name=names(anno_toxrefdb_nonneoplastic), annotation_num=anno_toxrefdb_nonneoplastic)
        }
    }
    
    # PubChem
    if(input[[paste0("switch__pubchem", switch_mode)]] == TRUE){
        anno_pubchem <- run_query(paste0("
            SELECT
                bcs.epa_id,
                pb.bioassay_group AS annotation
            FROM
                base_chemical_to_pubchem_cid bcpc,
                pubchem_bioactivity_mapping pbm,
                pubchem_bioassays pb,
                base_chemical_to_smiles bcs
            WHERE
                bcpc.pubchem_cid = pbm.cid
            AND bcpc.epa_id = bcs.epa_id
            AND pbm.aid = pb.aid
            AND bcpc.epa_id IN (
                ", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), "
            )
        "), args=as.list(selected))
        if(nrow(anno_pubchem) > 0){
            anno_pubchem <- split(anno_pubchem, anno_pubchem$annotation)
            anno_pubchem <- unlist(lapply(anno_pubchem, function(x) length(unique(x$epa_id))))
            anno_pubchem <- data.frame(annotation_name=names(anno_pubchem), annotation_num=anno_pubchem)
        }
    }
    
    if(input[[paste0("switch__pubchem_props", switch_mode)]] == TRUE){
        anno_pubchem_props <- run_query(paste0("
            SELECT DISTINCT
                bcs.epa_id,
                pba.pubchem_attribute AS annotation
            FROM
                base_chemical_to_pubchem_cid bcpc,
                pubchem_attributes pba,
                base_chemical_to_smiles bcs
            WHERE
                bcpc.pubchem_cid = pba.cid
            AND bcpc.epa_id = bcs.epa_id
            AND bcpc.epa_id IN (
                ", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), "
            )
        "), args=as.list(selected))
        if(nrow(anno_pubchem_props) > 0){
            anno_pubchem_props <- split(anno_pubchem_props, anno_pubchem_props$annotation)
            anno_pubchem_props <- unlist(lapply(anno_pubchem_props, function(x) length(unique(x$epa_id))))
            anno_pubchem_props <- data.frame(annotation_name=names(anno_pubchem_props), annotation_num=anno_pubchem_props)
        }
    }
    
    # GRAS
    if(input[[paste0("switch__gras", switch_mode)]] == TRUE){
        anno_gras <- run_query(paste0("
            SELECT
                bcs.epa_id,
                bcg.scogs_report_number AS annotation
            FROM
                base_chemical_to_gras bcg,
                base_chemical_to_smiles bcs
            WHERE
                bcg.epa_id = bcs.epa_id
            AND bcg.epa_id IN (
                ", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), "
            )
        "), args=as.list(selected))
        if(nrow(anno_gras) > 0){
            anno_gras <- split(anno_gras, anno_gras$annotation)
            anno_gras <- unlist(lapply(anno_gras, function(x) length(unique(x$epa_id))))
            anno_gras <- data.frame(annotation_name=names(anno_gras), annotation_num=anno_gras)
        }
    }
    
    # FooDB
    if(input[[paste0("switch__foodb_enzymes", switch_mode)]] == TRUE){
        anno_foodb_enzymes <- run_query(paste0("
            SELECT
                bcs.epa_id,
                fce.enzyme_name AS annotation
            FROM
                base_chemical_to_foodb_compound bcf,
                foodb_compound_enzymes fce,
                base_chemical_to_smiles bcs
            WHERE
                bcf.epa_id = bcs.epa_id
            AND bcf.foodb_compound_id = fce.compound_id
            AND bcf.epa_id IN (
                ", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), "
            )
        "), args=as.list(selected))
        if(nrow(anno_foodb_enzymes) > 0){
            anno_foodb_enzymes <- split(anno_foodb_enzymes, anno_foodb_enzymes$annotation)
            anno_foodb_enzymes <- unlist(lapply(anno_foodb_enzymes, function(x) length(unique(x$epa_id))))
            anno_foodb_enzymes <- data.frame(annotation_name=names(anno_foodb_enzymes), annotation_num=anno_foodb_enzymes)
        }
    }
    
    if(input[[paste0("switch__foodb_flavors", switch_mode)]] == TRUE){
        anno_foodb_flavors <- run_query(paste0("
            SELECT
                bcs.epa_id,
                fcf.flavor_name AS annotation
            FROM
                base_chemical_to_foodb_compound bcf,
                foodb_compound_flavors fcf,
                base_chemical_to_smiles bcs
            WHERE
                bcf.epa_id = bcs.epa_id
            AND bcf.foodb_compound_id = fcf.compound_id
            AND bcf.epa_id IN (
                ", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), "
            )
        "), args=as.list(selected))
        if(nrow(anno_foodb_flavors) > 0){
            anno_foodb_flavors <- split(anno_foodb_flavors, anno_foodb_flavors$annotation)
            anno_foodb_flavors <- unlist(lapply(anno_foodb_flavors, function(x) length(unique(x$epa_id))))
            anno_foodb_flavors <- data.frame(annotation_name=names(anno_foodb_flavors), annotation_num=anno_foodb_flavors)
        }
    }
    
    if(input[[paste0("switch__foodb_foodcontent", switch_mode)]] == TRUE){
        anno_foodb_foodcontent <- run_query(paste0("
            SELECT
                bcs.epa_id,
                fcfc.name AS annotation
            FROM
                base_chemical_to_foodb_compound bcf,
                foodb_compound_food_content fcfc,
                base_chemical_to_smiles bcs
            WHERE
                bcf.epa_id = bcs.epa_id
            AND bcf.foodb_compound_id = fcfc.compound_id
            AND bcf.epa_id IN (
                ", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), "
            )
        "), args=as.list(selected))
        if(nrow(anno_foodb_foodcontent) > 0){
            anno_foodb_foodcontent <- split(anno_foodb_foodcontent, anno_foodb_foodcontent$annotation)
            anno_foodb_foodcontent <- unlist(lapply(anno_foodb_foodcontent, function(x) length(unique(x$epa_id))))
            anno_foodb_foodcontent <- data.frame(annotation_name=names(anno_foodb_foodcontent), annotation_num=anno_foodb_foodcontent)
        }
    }
    
    if(input[[paste0("switch__foodb_healtheffects", switch_mode)]] == TRUE){
        anno_foodb_healtheffects <- run_query(paste0("
            SELECT
                bcs.epa_id,
                fche.health_effect_name AS annotation
            FROM
                base_chemical_to_foodb_compound bcf,
                foodb_compound_health_effects fche,
                base_chemical_to_smiles bcs
            WHERE
                bcf.epa_id = bcs.epa_id
            AND bcf.foodb_compound_id = fche.compound_id
            AND bcf.epa_id IN (
                ", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), "
            )
        "), args=as.list(selected))
        if(nrow(anno_foodb_healtheffects) > 0){
            anno_foodb_healtheffects <- split(anno_foodb_healtheffects, anno_foodb_healtheffects$annotation)
            anno_foodb_healtheffects <- unlist(lapply(anno_foodb_healtheffects, function(x) length(unique(x$epa_id))))
            anno_foodb_healtheffects <- data.frame(annotation_name=names(anno_foodb_healtheffects), annotation_num=anno_foodb_healtheffects)
        }
    }
    
    if(input[[paste0("switch__foodb_ontology", switch_mode)]] == TRUE){
        anno_foodb_ontology <- run_query(paste0("
            SELECT
                bcs.epa_id,
                fcot.term AS annotation
            FROM
                base_chemical_to_foodb_compound bcf,
                foodb_compound_ontology_terms fcot,
                base_chemical_to_smiles bcs
            WHERE
                bcf.epa_id = bcs.epa_id
            AND bcf.foodb_compound_id = fcot.compound_id
            AND bcf.epa_id IN (
                ", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), "
            )
        "), args=as.list(selected))
        if(nrow(anno_foodb_ontology) > 0){
            anno_foodb_ontology <- split(anno_foodb_ontology, anno_foodb_ontology$annotation)
            anno_foodb_ontology <- unlist(lapply(anno_foodb_ontology, function(x) length(unique(x$epa_id))))
            anno_foodb_ontology <- data.frame(annotation_name=names(anno_foodb_ontology), annotation_num=anno_foodb_ontology)
        }
    }
    
    if(input[[paste0("switch__epa_props", switch_mode)]] == TRUE){
        anno_epa_props <- run_query(paste0("
            SELECT
                bcs.epa_id,
                ep.epa_attribute AS annotation,
                bc.dsstox_substance_id
            FROM
                epa_properties ep,
                base_chemical_to_smiles bcs,
                base_chemicals bc
            WHERE
                ep.epa_id = bcs.epa_id
            AND ep.epa_id = bc.epa_id
            AND ep.epa_id IN (
                ", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), "
            )
        "), args=as.list(selected))
        if(nrow(anno_epa_props) > 0){
            anno_epa_props <- split(anno_epa_props, anno_epa_props$annotation)
            anno_epa_props <- unlist(lapply(anno_epa_props, function(x) length(unique(x$epa_id))))
            anno_epa_props <- data.frame(annotation_name=names(anno_epa_props), annotation_num=anno_epa_props)
        }
    }
    
    return(list(
        anno_admet_binr=anno_admet_binr,
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

# Load chemical annotations for a given cluster
chemical_annotations <- function(selected, input, selected_node, dtxsids, switch_mode=""){
    
    #print(input)
    
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
    # if(input[[paste0("switch__admet_binr", switch_mode)]] == TRUE){
    #     anno_admet_binr_plot <- run_query(paste0("
    #                     SELECT DISTINCT
    #                         bc.dsstox_substance_id,
    #                         abcp.model_results
    #                     FROM
    #                         base_chemical_to_smiles bcs,
    #                         base_chemicals bc,
    #                         admet_base_chemical_predictions_interpretation_binary_v11_json abcp
    #                     WHERE
    #                         bcs.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
    #                     AND bcs.epa_id = bc.epa_id
    #                     AND bcs.smi_id = abcp.smi_id
    #                 "), args=as.list(selected))
    #     anno_admet_binr_plot_rows <- lapply(anno_admet_binr_plot$model_results, function(x){
    #         fromJSON(x)
    #     })
    #     anno_admet_binr_plot_rows <- rbindlist(anno_admet_binr_plot_rows) %>% mutate_if(is.character, as.numeric)
    #     anno_admet_binr_plot_rows <- cbind(anno_admet_binr_plot$dsstox_substance_id, anno_admet_binr_plot_rows)
    #     colnames(anno_admet_binr_plot_rows)[1] <- "dsstox_substance_id"
    #     #anno_admet_binr_hm <- annotation_heatmap(anno_admet_binr_plot_rows)
    # }
    # 
    # if(input[[paste0("switch__admet_catg", switch_mode)]] == TRUE){
    #     anno_admet_catg_plot <- run_query(paste0("
    #                     SELECT DISTINCT
    #                         bc.dsstox_substance_id,
    #                         abcp.model_results
    #                     FROM
    #                         base_chemical_to_smiles bcs,
    #                         base_chemicals bc,
    #                         admet_base_chemical_predictions_interpretation_catego_v11_json abcp
    #                     WHERE
    #                         bcs.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
    #                     AND bcs.epa_id = bc.epa_id
    #                     AND bcs.smi_id = abcp.smi_id
    #                 "), args=as.list(selected))
    #     anno_admet_catg_plot_rows <- lapply(anno_admet_catg_plot$model_results, function(x){
    #         fromJSON(x)
    #     })
    #     anno_admet_catg_plot_rows <- rbindlist(anno_admet_catg_plot_rows) %>% mutate_if(is.character, as.numeric)
    #     anno_admet_catg_plot_rows <- cbind(anno_admet_catg_plot$dsstox_substance_id, anno_admet_catg_plot_rows)
    #     colnames(anno_admet_catg_plot_rows)[1] <- "dsstox_substance_id"
    #     
    #     #n_electr messes up scale
    #     #anno_admet_catg_plot_rows$N_Electr <- NULL
    #     
    #     #anno_admet_catg_hm <- annotation_heatmap(anno_admet_catg_plot_rows)
    # }
    # 
    # if(input[[paste0("switch__admet_cont", switch_mode)]] == TRUE){
    #     anno_admet_cont_plot <- run_query(paste0("
    #                     SELECT DISTINCT
    #                         bc.dsstox_substance_id,
    #                         abcp.model_results
    #                     FROM
    #                         base_chemical_to_smiles bcs,
    #                         base_chemicals bc,
    #                         admet_base_chemical_predictions_interpretation_contin_v11_json abcp
    #                     WHERE
    #                         bcs.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
    #                     AND bcs.epa_id = bc.epa_id
    #                     AND bcs.smi_id = abcp.smi_id
    #                 "), args=as.list(selected))
    #     anno_admet_cont_plot_rows <- lapply(anno_admet_cont_plot$model_results, function(x){
    #         fromJSON(x)
    #     })
    #     anno_admet_cont_plot_rows <- rbindlist(anno_admet_cont_plot_rows) %>% mutate_if(is.character, as.numeric)
    #     anno_admet_cont_plot_rows <- cbind(anno_admet_cont_plot$dsstox_substance_id, anno_admet_cont_plot_rows)
    #     colnames(anno_admet_cont_plot_rows)[1] <- "dsstox_substance_id"
    #     anno_admet_cont_plot_rows <- anno_admet_cont_plot_rows %>% select_if(~ !any(is.na(.)))
    #     #anno_admet_cont_hm <- annotation_heatmap(anno_admet_cont_plot_rows)
    # }
    
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
                iaea.assay_endpoint_value
            FROM
                invitrodb_summarized_results isr,
                invitrodb_annotation ia,
                invitrodb_assay_endpoint_annotation iaea,
                base_chemicals bc,
                base_chemical_to_smiles bcs,
                base_chemical_compounds bcc
            WHERE
                bc.epa_id = bcs.epa_id
            AND bc.epa_id IN (", paste0(lapply(seq_len(length(selected)), function(x) paste0("$", x)), collapse=", "), ")
            AND bc.epa_id = bcc.epa_id
            AND bcs.smi_id = isr.smi_id
            AND isr.invitrodb_attributes = ia.invitrodb_attributes
            AND isr.invitrodb_attributes = 'aenm'
            AND isr.invitrodb_values = iaea.aenm
        "), args=as.list(selected))
    }
    
    # Leadscope
    if(input[[paste0("switch__leadscope", switch_mode)]] == TRUE){
        anno_leadscope <- run_query(paste0("
                        SELECT DISTINCT
                            bcc.preferred_name,
                            bcc.casrn,
                            bc.dsstox_substance_id,
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
                            fcfc.subclass
                        FROM
                            base_chemical_to_foodb_compound bcfc,
                            foodb_compound_food_content fcfc,
                            base_chemicals bc,
                            base_chemical_compounds bcc
                        WHERE
                            bcfc.epa_id = bc.epa_id
                        AND bcfc.foodb_compound_id = fcfc.compound_id
                        AND bcfc.epa_id = bcc.epa_id
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