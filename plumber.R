#
# This is a Plumber API. You can run the API by clicking
# the 'Run API' button above.
#
# Find out more about building APIs with Plumber here:
#
#    https://www.rplumber.io/
#

library(plumber)
library(bit64)
library(data.table)
library(DBI)
library(dplyr)
library(DT)
library(heatmap3)
library(heatmaply)
library(htmlwidgets)
library(openxlsx)
library(plotly)
library(pool)
library(reshape2)
library(rsvg)
library(shiny)
library(shinycssloaders)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyjs)
library(shinyWidgets)
library(sortable)
library(stringi)
library(stringr)
library(svglite)
library(RPostgres)
library(visNetwork)
# Source other code files
source("modules/enrichment.R")
source("modules/file.R")
source("modules/ice.R")
source("modules/leadscope.R")
source("modules/metabolite.R")
source("modules/plot.R")
source("modules/report.R")
source("modules/postgres.R")

# Redefine %like% operator to be case insensitive (from https://stackoverflow.com/questions/41425699/how-to-get-the-like-operator-to-be-case-insensitive)
`%like%` <- function (x, pattern) { 
    stringi::stri_detect_regex(x, pattern, case_insensitive=TRUE)
}

#* @apiTitle ChemBioTox API
#* @apiDescription This API is for interfacing with the ChemBioTox (CBT) database, a comprehensive chemical database built on Postgres that contains property and toxicological annotations for over 1 million chemicals of interest to the EPA. Data within CBT is sourced from various publicly available databases (i.e., CTD), calculated using bioinformatics tools (i.e., RDKit), or produced from proprietary technology from Instem and Simulations Plus (not available in the public version of the API: you must log in using NIH credentials to access these endpoints). <br><br> Datasources and tools used in the development of the ChemBioTox (CBT) database include:<br> <ul> <li> ChEMBL v3.4 - Zdrazil B, et al. The ChEMBL Database in 2023: a drug discovery platform spanning multiple bioactivity data types and time periods. Nucleic Acids Research, Volume 52, Issue D1, 5 January 2024, Pages D1180–D1192, <a href="https://doi.org/10.1093/nar/gkad1004">https://doi.org/10.1093/nar/gkad1004</a> </li><li> Chemical and Products Database (CPDat) - Dionisio K, Philips K, Price P, et al. The Chemical and Products Database, a resource for exposure-relevant data on chemicals in consumer products. Sci Data 5, 180125 (2018). <a href="https://doi.org/10.1038/sdata.2018.125"> https://doi.org/10.1038/sdata.2018.125</a> </li><li> Comparative Toxicogenomics Database (CTD) - Davis AP, Wiegers TC, Johnson RJ, Sciaky D, Wiegers J, Mattingly CJ. Comparative Toxicogenomics Database (CTD): update 2023. Nucleic Acids Res. 2022 Sep 28. <a href="https://ctdbase.org/">https://ctdbase.org/</a> </li><li> DrugBank Online - Knox C, et al. DrugBank 6.0: the DrugBank Knowledgebase for 2024. Nucleic Acids Res. 2024 Jan 5;52(D1):D1265-D1275. <a href="https://doi.org/10.1093/nar/gkad976">https://doi.org/10.1093/nar/gkad976</a>, <a href="https://go.drugbank.com/">https://go.drugbank.com/</a> </li><li> EPA CompTox Chemicals Dashboard v2.4.1 - Williams AJ, Grulke CM, Edwards J, et al. The CompTox Chemistry Dashboard: a community data resource for environmental chemistry. J Cheminform 9, 61 (2017). <a href="https://doi.org/10.1186/s13321-017-0247-6">https://doi.org/10.1186/s13321-017-0247-6</a>, <a href="https://comptox.epa.gov/dashboard/">https://comptox.epa.gov/dashboard/</a> </li><li> FooDB v1.0 - https://foodb.ca</li><li> Generally Recognized As Safe (GRAS) Substances SCOGS Database - US FDA <a href="https://www.cfsanappsexternal.fda.gov/scripts/fdcc/?set=SCOGS">https://www.cfsanappsexternal.fda.gov/scripts/fdcc/?set=SCOGS</a>, <a href="https://www.fda.gov/food/generally-recognized-safe-gras/gras-substances-scogs-database">https://www.fda.gov/food/generally-recognized-safe-gras/gras-substances-scogs-database</a> </li><li> Human Metabolome Database (HMDB) - Wishart DS, Guo A, Oler E, Wang F, Anjum A, Peters H, Dizon R, Sayeeda Z, Tian S, Lee BL, Berjanskii M, Mah R, Yamamoto M, Jovel J, Torres-Calzada C, Hiebert-Giesbrecht M, Lui VW, Varshavi D, Varshavi D, Allen D, Arndt D, Khetarpal N, Sivakumaran A, Harford K, Sanford S, Yee K, Cao X, Budinski Z, Liigand J, Zhang L, Zheng J, Mandal R, Karu N, Dambrova M, Schiöth HB, Greiner R, Gautam V. HMDB 5.0: the Human Metabolome Database for 2022. Nucleic Acids Res. 2022 Jan 7;50(D1):D622-D631. <a href="https://doi.org/10.1093/nar/gkab1062">https://doi.org/10.1093/nar/gkab1062</a>. PMID: 34986597; PMCID: PMC8728138.</li><li> InVitroDB v3.5 - <a href="https://www.epa.gov/comptox-tools/exploring-toxcast-data">https://www.epa.gov/comptox-tools/exploring-toxcast-data</a> </li><li> OCHEM - Sushko I, et al. Online chemical modeling environment (OCHEM): web platform for data storage, model development and publishing of chemical information. J Comput Aided Mol Des. 2011; 25(6):533-54. <a href="https://doi.org/10.1007/s10822-011-9440-2">https://doi.org/10.1007/s10822-011-9440-2</a>, <a href="https://ochem.eu/home/show.do">https://ochem.eu/home/show.do</a> </li><li> PubChem - Kim S, Chen J, Cheng T, et al. PubChem 2023 update. Nucleic Acids Res. 2023;51(D1):D1373–D1380. <a href="https://doi.org/10.1093/nar/gkac956">https://doi.org/10.1093/nar/gkac956</a>, <a href="https://pubchem.ncbi.nlm.nih.gov/ ">https://pubchem.ncbi.nlm.nih.gov/</a> </li><li> RDKit - <a href="https://github.com/rdkit/rdkit">https://github.com/rdkit/rdkit</a> </li><li> Toxin and Toxin Target Database (T3DB) - Wishart D, et al. T3DB: the toxic exposome database. Nucleic Acids Res. 2015 Jan;43(Database issue):D928-34. <a href="ttps://doi.org/10.1093/nar/gku1004">https://doi.org/10.1093/nar/gku1004</a>; Lim E, et al. T3DB: a comprehensively annotated database of common toxins and their targets. Nucleic Acids Res. 2010 Jan 38(Database issue):D781-6. <a href="https://doi.org/10.1093/nar/gkp934">https://doi.org/10.1093/nar/gkp934</a>; <a href="http://www.t3db.ca/">http://www.t3db.ca/</a> </li><li> ToxRefDB - Martin MT and Judson R. ToxRefDB - Release user-friendly web-based tool for mining ToxRefDB. U.S. Environmental Protection Agency, Washington, DC, 2010. <a href="https://cfpub.epa.gov/si/si_public_record_report.cfm?Lab=NCCT&dirEntryId=227139">https://cfpub.epa.gov/si/si_public_record_report.cfm?Lab=NCCT&dirEntryId=227139</a>, <a href="https://www.epa.gov/comptox-tools/exploring-toxcast-data">https://www.epa.gov/comptox-tools/exploring-toxcast-data</a> </ul>


# #* Get a list of all chemicals supported by ChemBioTox
# #* @get /lists/base_chemicals
# #* @tag "Lists"
# function() {
#     # Get internal IDs for given chemical
#     chemicals <- run_query(paste0("
#         SELECT DISTINCT
#             bc.dsstox_substance_id,
#             bcc.dsstox_compound_id,
#             bcc.preferred_name,
#             bcc.casrn
#         FROM
#             base_chemicals bc
#         LEFT JOIN base_chemical_compounds bcc
#             ON bc.epa_id = bcc.epa_id
#     "))
# }

#* Given a chemical name, return the DSSTox substance ID.
#* @param name chemical name
#* @get /name2dtxsid
#* @tag "Identifier Lookup"
function(name = "") {
    # Get internal IDs for given chemical
    internal_id <- run_query(paste0("
        SELECT DISTINCT
            bc.dsstox_substance_id
        FROM
            base_chemicals bc,
            base_chemical_compounds bcc
        WHERE
            UPPER(bcc.preferred_name) = UPPER($1)
            AND bc.epa_id = bcc.epa_id
    "), args=as.list(name))
}

#* Given a CASRN, return the DSSTox substance ID.
#* @param name chemical name
#* @get /casrn2dtxsid
#* @tag "Identifier Lookup"
function(casrn = "") {
    # Get internal IDs for given chemical
    internal_id <- run_query(paste0("
        SELECT DISTINCT
            bc.dsstox_substance_id
        FROM
            base_chemicals bc,
            base_chemical_compounds bcc
        WHERE
            bcc.casrn = $1
            AND bc.epa_id = bcc.epa_id
    "), args=as.list(casrn))
}

#* Given a SMILES string, return the DSSTox substance ID. If the given SMILES is not in the database, will return the closest match using rdkit fingerprint-based structural similarity.
#* @param smiles SMILES string
#* @get /smiles2dtxsid
#* @tag "Identifier Lookup"
function(smiles = "") {
    # Get internal IDs for given chemical
    internal_id <- run_query(paste0("
        SELECT DISTINCT
            bc.dsstox_substance_id,
            sim.similarity
        FROM
            base_chemicals bc,
            base_chemical_to_smiles bcs,
            get_rdkit_fp_neighbors($1) sim
        WHERE
            sim.smi_id = bcs.smi_id
            AND bcs.epa_id = bc.epa_id
            AND sim.similarity >= 0.5
        ORDER BY sim.similarity DESC
        LIMIT 1
    "), args=as.list(smiles))
}

#* Given a DSSTox substance ID, return the canonical SMILES string.
#* @param dtxsid DSSTox substance ID
#* @get /dtxsid2smiles
#* @tag "Identifier Lookup"
function(dtxsid = "") {
    # Get internal IDs for given chemical
    internal_id <- run_query(paste0("
        SELECT DISTINCT
            css.canonical_smiles
        FROM
            base_chemicals bc,
            base_chemical_to_smiles bcs,
            smiles_to_canonical_smiles scs,
            canonical_smiles_strings css
        WHERE
            bc.dsstox_substance_id = $1
        AND bcs.epa_id = bc.epa_id
        AND bcs.smi_id = scs.smi_id
        AND scs.csm_id = css.csm_id
    "), args=as.list(dtxsid))
}

#* Given a DSSTox substance ID, return the corresponding PubChem compound ID.
#* @param dtxsid DSSTox substance ID
#* @get /cid
#* @tag "Identifier Lookup"
function(dtxsid = "") {
    # Get internal IDs for given chemical
    internal_id <- run_query(paste0("
        SELECT DISTINCT
            bpc.pubchem_cid
        FROM
            base_chemicals bc,
            base_chemical_to_pubchem_cid bpc
        WHERE
            bc.dsstox_substance_id = $1
        AND bc.epa_id = bpc.epa_id
    "), args=as.list(dtxsid))
}


#* Given a DSSTox substance ID, return a chemical binary vector fingerprint generated with RDKit. Supported types are: "atompair", "avalon", "maccskeys", "morgan", "rdkit", "rdklayer", "rdkpattern", "topologicaltorsion", default is "rdkit".
#* @param dtxsid DSSTox substance ID
#* @param mode Fingerprint mode. Supported modes are: "atompair", "avalon", "maccskeys", "morgan", "rdkit", "rdklayer", "rdkpattern", "topologicaltorsion"
#* @get /fingerprint
#* @tag "Identifier Lookup"
function(dtxsid="", mode="rdkit") {
    
    if(!(mode %in% c("atompair", "avalon", "maccskeys", "morgan", "rdkit", "rdklayer", "rdkpattern", "topologicaltorsion"))){
        return(paste0("Error: unsupported fingerprint type: ", mode))
    }
    
    # Get internal IDs for given chemical
    fp <- run_query(paste0("
        SELECT DISTINCT
            fp.bv
        FROM
            base_chemicals bc,
            base_chemical_to_smiles bcs,
            base_chemical_fingerprint_", mode, " fp
        WHERE
            bc.dsstox_substance_id = $1
        AND bcs.epa_id = bc.epa_id
        AND bcs.smi_id = fp.smi_id
    "), args=as.list(dtxsid))
}


#* Given a DSSTox substance ID, annotate using all datasources available in CBT.
#* @param dtxsid DSSTox substance ID
#* @get /annotate
#* @tag "Summary"
function(dtxsid = "") {
    # Get internal IDs for given chemical
    internal_id <- run_query(paste0("
        SELECT DISTINCT
            bc.epa_id,
            bc.dsstox_substance_id
        FROM
            base_chemicals bc
        WHERE
            bc.dsstox_substance_id = $1
    "), args=as.list(dtxsid))
    
    if(nrow(internal_id) > 0){
        
        input_all_annotations <- list(
            "switch__admet_binr"=FALSE,
            "switch__admet_catg"=FALSE,
            "switch__admet_cont"=FALSE,
            "switch__chembl"=TRUE,
            "switch__cpd"=TRUE,
            "switch__ctd_bioprocess"=TRUE,
            "switch__ctd_cellcomp"=TRUE,
            "switch__ctd_diseases"=TRUE,
            "switch__ctd_genes"=TRUE,
            "switch__ctd_molfunct"=TRUE,
            "switch__ctd_phenotypes"=TRUE,
            "switch__drugbank_atccodes"=TRUE,
            "switch__drugbank_carriers"=TRUE,
            "switch__drugbank_enzymes"=TRUE,
            "switch__drugbank_targets"=TRUE,
            "switch__drugbank_transporters"=TRUE,
            "switch__hmdb_biospecimenlocations"=TRUE,
            "switch__hmdb_cellularlocations"=TRUE,
            "switch__hmdb_diseases"=TRUE,
            "switch__hmdb_genes"=TRUE,
            "switch__hmdb_tissuelocations"=TRUE,
            "switch__ice"=TRUE,
            "switch__invitrodb"=TRUE,
            "switch__leadscope"=FALSE,
            "switch__ochem"=TRUE,
            "switch__opera_adme"=TRUE,
            "switch__opera_environmental_fate"=TRUE,
            "switch__opera_physicochem"=TRUE,
            "switch__opera_structural_properties"=TRUE,
            "switch__opera_toxicology_endpoints"=TRUE,
            "switch__saagar"=TRUE,
            "switch__t3db"=TRUE,
            "switch__toxrefdb_neoplastic"=TRUE,
            "switch__toxrefdb_nonneoplastic"=TRUE,
            "switch__pubchem"=TRUE,
            "switch__pubchem_props"=TRUE,
            "switch__gras"=TRUE,
            "switch__foodb_enzymes"=TRUE,
            "switch__foodb_flavors"=TRUE,
            "switch__foodb_foodcontent"=TRUE,
            "switch__foodb_healtheffects"=TRUE,
            "switch__foodb_ontology"=TRUE,
            "switch__epa_props"=TRUE
        )
        annotations <- chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode="")
    }
}

# #* Given a DSSTox substance ID, return chemical properties from GenRA.
# #* @param dtxsid DSSTox substance ID
# #* @get /genra_properties
# #* @tag "Chemical Properties"
# function(dtxsid="") {
#     # Get internal IDs for given chemical
#     prop <- run_query(paste0("
#         SELECT DISTINCT
#             fp.bv
#         FROM
#             base_chemicals bc,
#             base_chemical_to_smiles bcs,
#             base_chemical_fingerprint_", mode, " fp
#         WHERE
#             bc.dsstox_substance_id = $1
#         AND bcs.epa_id = bc.epa_id
#         AND bcs.smi_id = fp.smi_id
#     "), args=as.list(dtxsid))
# }

#* Given a DSSTox substance ID, annotate using the Chemical and Products Database (CPDat) available in CBT. Contains information about a chemical's usage in commercially-available products.
#* @param dtxsid DSSTox substance ID
#* @get /cpdat
#* @tag "Commercial and Industrial Usage"
#* @tag "Environmental Fate and Exposure Annotations"
function(dtxsid = "") {
    # Get internal IDs for given chemical
    internal_id <- run_query(paste0("
        SELECT DISTINCT
            bc.epa_id,
            bc.dsstox_substance_id
        FROM
            base_chemicals bc
        WHERE
            bc.dsstox_substance_id = $1
    "), args=as.list(dtxsid))
    
    if(nrow(internal_id) > 0){
        
        input_all_annotations <- list(
            "switch__admet_binr"=FALSE,
            "switch__admet_catg"=FALSE,
            "switch__admet_cont"=FALSE,
            "switch__chembl"=FALSE,
            "switch__cpd"=TRUE,
            "switch__ctd_bioprocess"=FALSE,
            "switch__ctd_cellcomp"=FALSE,
            "switch__ctd_diseases"=FALSE,
            "switch__ctd_genes"=FALSE,
            "switch__ctd_molfunct"=FALSE,
            "switch__ctd_phenotypes"=FALSE,
            "switch__drugbank_atccodes"=FALSE,
            "switch__drugbank_carriers"=FALSE,
            "switch__drugbank_enzymes"=FALSE,
            "switch__drugbank_targets"=FALSE,
            "switch__drugbank_transporters"=FALSE,
            "switch__hmdb_biospecimenlocations"=FALSE,
            "switch__hmdb_cellularlocations"=FALSE,
            "switch__hmdb_diseases"=FALSE,
            "switch__hmdb_genes"=FALSE,
            "switch__hmdb_tissuelocations"=FALSE,
            "switch__ice"=FALSE,
            "switch__invitrodb"=FALSE,
            "switch__leadscope"=FALSE,
            "switch__ochem"=FALSE,
            "switch__opera_adme"=FALSE,
            "switch__opera_environmental_fate"=FALSE,
            "switch__opera_physicochem"=FALSE,
            "switch__opera_structural_properties"=FALSE,
            "switch__opera_toxicology_endpoints"=FALSE,
            "switch__saagar"=FALSE,
            "switch__t3db"=FALSE,
            "switch__toxrefdb_neoplastic"=FALSE,
            "switch__toxrefdb_nonneoplastic"=FALSE,
            "switch__pubchem"=FALSE,
            "switch__pubchem_props"=FALSE,
            "switch__gras"=FALSE,
            "switch__foodb_enzymes"=FALSE,
            "switch__foodb_flavors"=FALSE,
            "switch__foodb_foodcontent"=FALSE,
            "switch__foodb_healtheffects"=FALSE,
            "switch__foodb_ontology"=FALSE,
            "switch__epa_props"=FALSE
        )
        annotations <- chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode="")[["anno_cpd"]]
    }
}


#* Given a DSSTox substance ID, return a structural image representing the corresponding chemical in SVG format generated by RDKit.
#* @param dtxsid DSSTox substance ID
#* @get /image/structure
#* @tag "Visualization"
function(dtxsid = "") {
    # Get internal IDs for given chemical
    img <- run_query(paste0("
        SELECT DISTINCT
            CONCAT(bci.svg1, bci.svg2) AS image
        FROM
            base_chemicals bc,
            base_chemical_to_smiles bcs,
            base_chemical_images bci
        WHERE
            bc.dsstox_substance_id = $1
        AND bc.epa_id = bcs.epa_id
        AND bcs.smi_id = bci.smi_id
    "), args=as.list(dtxsid))
    return(img)
}



#* Given a DSSTox substance ID, annotate using the FooDB available in CBT. Contains information about a chemical's usage in food products and flavoring.
#* @param dtxsid DSSTox substance ID
#* @get /foodb
#* @tag "Commercial and Industrial Usage"
#* @tag "Environmental Fate and Exposure Annotations"
function(dtxsid = "") {
    # Get internal IDs for given chemical
    internal_id <- run_query(paste0("
        SELECT DISTINCT
            bc.epa_id,
            bc.dsstox_substance_id
        FROM
            base_chemicals bc
        WHERE
            bc.dsstox_substance_id = $1
    "), args=as.list(dtxsid))
    
    if(nrow(internal_id) > 0){
        
        input_all_annotations <- list(
            "switch__admet_binr"=FALSE,
            "switch__admet_catg"=FALSE,
            "switch__admet_cont"=FALSE,
            "switch__chembl"=FALSE,
            "switch__cpd"=FALSE,
            "switch__ctd_bioprocess"=FALSE,
            "switch__ctd_cellcomp"=FALSE,
            "switch__ctd_diseases"=FALSE,
            "switch__ctd_genes"=FALSE,
            "switch__ctd_molfunct"=FALSE,
            "switch__ctd_phenotypes"=FALSE,
            "switch__drugbank_atccodes"=FALSE,
            "switch__drugbank_carriers"=FALSE,
            "switch__drugbank_enzymes"=FALSE,
            "switch__drugbank_targets"=FALSE,
            "switch__drugbank_transporters"=FALSE,
            "switch__hmdb_biospecimenlocations"=FALSE,
            "switch__hmdb_cellularlocations"=FALSE,
            "switch__hmdb_diseases"=FALSE,
            "switch__hmdb_genes"=FALSE,
            "switch__hmdb_tissuelocations"=FALSE,
            "switch__ice"=FALSE,
            "switch__invitrodb"=FALSE,
            "switch__leadscope"=FALSE,
            "switch__ochem"=FALSE,
            "switch__opera_adme"=FALSE,
            "switch__opera_environmental_fate"=FALSE,
            "switch__opera_physicochem"=FALSE,
            "switch__opera_structural_properties"=FALSE,
            "switch__opera_toxicology_endpoints"=FALSE,
            "switch__saagar"=FALSE,
            "switch__t3db"=FALSE,
            "switch__toxrefdb_neoplastic"=FALSE,
            "switch__toxrefdb_nonneoplastic"=FALSE,
            "switch__pubchem"=FALSE,
            "switch__pubchem_props"=FALSE,
            "switch__gras"=FALSE,
            "switch__foodb_enzymes"=TRUE,
            "switch__foodb_flavors"=TRUE,
            "switch__foodb_foodcontent"=TRUE,
            "switch__foodb_healtheffects"=TRUE,
            "switch__foodb_ontology"=TRUE,
            "switch__epa_props"=FALSE
        )
        annotations <- chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode="")[c("anno_foodb_enzymes", "anno_foodb_flavors", "anno_foodb_foodcontent", "anno_foodb_healtheffects", "anno_foodb_ontology")]
    }
}

#* Given a DSSTox substance ID, annotate using the Generally Recognized As Safe (GRAS) database available in CBT. Contains SCOGS reports for chemicals.
#* @param dtxsid DSSTox substance ID
#* @get /gras
#* @tag "Safety Studies"
function(dtxsid = "") {
    # Get internal IDs for given chemical
    internal_id <- run_query(paste0("
        SELECT DISTINCT
            bc.epa_id,
            bc.dsstox_substance_id
        FROM
            base_chemicals bc
        WHERE
            bc.dsstox_substance_id = $1
    "), args=as.list(dtxsid))
    
    if(nrow(internal_id) > 0){
        
        input_all_annotations <- list(
            "switch__admet_binr"=FALSE,
            "switch__admet_catg"=FALSE,
            "switch__admet_cont"=FALSE,
            "switch__chembl"=FALSE,
            "switch__cpd"=FALSE,
            "switch__ctd_bioprocess"=FALSE,
            "switch__ctd_cellcomp"=FALSE,
            "switch__ctd_diseases"=FALSE,
            "switch__ctd_genes"=FALSE,
            "switch__ctd_molfunct"=FALSE,
            "switch__ctd_phenotypes"=FALSE,
            "switch__drugbank_atccodes"=FALSE,
            "switch__drugbank_carriers"=FALSE,
            "switch__drugbank_enzymes"=FALSE,
            "switch__drugbank_targets"=FALSE,
            "switch__drugbank_transporters"=FALSE,
            "switch__hmdb_biospecimenlocations"=FALSE,
            "switch__hmdb_cellularlocations"=FALSE,
            "switch__hmdb_diseases"=FALSE,
            "switch__hmdb_genes"=FALSE,
            "switch__hmdb_tissuelocations"=FALSE,
            "switch__ice"=FALSE,
            "switch__invitrodb"=FALSE,
            "switch__leadscope"=FALSE,
            "switch__ochem"=FALSE,
            "switch__opera_adme"=FALSE,
            "switch__opera_environmental_fate"=FALSE,
            "switch__opera_physicochem"=FALSE,
            "switch__opera_structural_properties"=FALSE,
            "switch__opera_toxicology_endpoints"=FALSE,
            "switch__saagar"=FALSE,
            "switch__t3db"=FALSE,
            "switch__toxrefdb_neoplastic"=FALSE,
            "switch__toxrefdb_nonneoplastic"=FALSE,
            "switch__pubchem"=FALSE,
            "switch__pubchem_props"=FALSE,
            "switch__gras"=TRUE,
            "switch__foodb_enzymes"=FALSE,
            "switch__foodb_flavors"=FALSE,
            "switch__foodb_foodcontent"=FALSE,
            "switch__foodb_healtheffects"=FALSE,
            "switch__foodb_ontology"=FALSE,
            "switch__epa_props"=FALSE
        )
        annotations <- chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode="")[["anno_gras"]]
    }
}

#* Given a DSSTox substance ID, annotate using the Toxin and Toxin-Target Database (T3DB) available in CBT.
#* @param dtxsid DSSTox substance ID
#* @get /t3db
#* @tag "Target Annotations"
function(dtxsid = "") {
    # Get internal IDs for given chemical
    internal_id <- run_query(paste0("
        SELECT DISTINCT
            bc.epa_id,
            bc.dsstox_substance_id
        FROM
            base_chemicals bc
        WHERE
            bc.dsstox_substance_id = $1
    "), args=as.list(dtxsid))
    
    if(nrow(internal_id) > 0){
        
        input_all_annotations <- list(
            "switch__admet_binr"=FALSE,
            "switch__admet_catg"=FALSE,
            "switch__admet_cont"=FALSE,
            "switch__chembl"=FALSE,
            "switch__cpd"=FALSE,
            "switch__ctd_bioprocess"=FALSE,
            "switch__ctd_cellcomp"=FALSE,
            "switch__ctd_diseases"=FALSE,
            "switch__ctd_genes"=FALSE,
            "switch__ctd_molfunct"=FALSE,
            "switch__ctd_phenotypes"=FALSE,
            "switch__drugbank_atccodes"=FALSE,
            "switch__drugbank_carriers"=FALSE,
            "switch__drugbank_enzymes"=FALSE,
            "switch__drugbank_targets"=FALSE,
            "switch__drugbank_transporters"=FALSE,
            "switch__hmdb_biospecimenlocations"=FALSE,
            "switch__hmdb_cellularlocations"=FALSE,
            "switch__hmdb_diseases"=FALSE,
            "switch__hmdb_genes"=FALSE,
            "switch__hmdb_tissuelocations"=FALSE,
            "switch__ice"=FALSE,
            "switch__invitrodb"=FALSE,
            "switch__leadscope"=FALSE,
            "switch__ochem"=FALSE,
            "switch__opera_adme"=FALSE,
            "switch__opera_environmental_fate"=FALSE,
            "switch__opera_physicochem"=FALSE,
            "switch__opera_structural_properties"=FALSE,
            "switch__opera_toxicology_endpoints"=FALSE,
            "switch__saagar"=FALSE,
            "switch__t3db"=TRUE,
            "switch__toxrefdb_neoplastic"=FALSE,
            "switch__toxrefdb_nonneoplastic"=FALSE,
            "switch__pubchem"=FALSE,
            "switch__pubchem_props"=FALSE,
            "switch__gras"=FALSE,
            "switch__foodb_enzymes"=FALSE,
            "switch__foodb_flavors"=FALSE,
            "switch__foodb_foodcontent"=FALSE,
            "switch__foodb_healtheffects"=FALSE,
            "switch__foodb_ontology"=FALSE,
            "switch__epa_props"=FALSE
        )
        annotations <- chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode="")[["anno_t3db"]]
    }
}

#* Given a DSSTox substance ID, annotate using the InVitroDB data available in CBT.
#* @param dtxsid DSSTox substance ID
#* @get /invitrodb
#* @tag "Chemical Properties"
#* @tag "In Vitro Annotations"
function(dtxsid = "") {
    # Get internal IDs for given chemical
    internal_id <- run_query(paste0("
        SELECT DISTINCT
            bc.epa_id,
            bc.dsstox_substance_id
        FROM
            base_chemicals bc
        WHERE
            bc.dsstox_substance_id = $1
    "), args=as.list(dtxsid))
    
    if(nrow(internal_id) > 0){
        
        input_all_annotations <- list(
            "switch__admet_binr"=FALSE,
            "switch__admet_catg"=FALSE,
            "switch__admet_cont"=FALSE,
            "switch__chembl"=FALSE,
            "switch__cpd"=FALSE,
            "switch__ctd_bioprocess"=FALSE,
            "switch__ctd_cellcomp"=FALSE,
            "switch__ctd_diseases"=FALSE,
            "switch__ctd_genes"=FALSE,
            "switch__ctd_molfunct"=FALSE,
            "switch__ctd_phenotypes"=FALSE,
            "switch__drugbank_atccodes"=FALSE,
            "switch__drugbank_carriers"=FALSE,
            "switch__drugbank_enzymes"=FALSE,
            "switch__drugbank_targets"=FALSE,
            "switch__drugbank_transporters"=FALSE,
            "switch__hmdb_biospecimenlocations"=FALSE,
            "switch__hmdb_cellularlocations"=FALSE,
            "switch__hmdb_diseases"=FALSE,
            "switch__hmdb_genes"=FALSE,
            "switch__hmdb_tissuelocations"=FALSE,
            "switch__ice"=FALSE,
            "switch__invitrodb"=TRUE,
            "switch__leadscope"=FALSE,
            "switch__ochem"=FALSE,
            "switch__opera_adme"=FALSE,
            "switch__opera_environmental_fate"=FALSE,
            "switch__opera_physicochem"=FALSE,
            "switch__opera_structural_properties"=FALSE,
            "switch__opera_toxicology_endpoints"=FALSE,
            "switch__saagar"=FALSE,
            "switch__t3db"=FALSE,
            "switch__toxrefdb_neoplastic"=FALSE,
            "switch__toxrefdb_nonneoplastic"=FALSE,
            "switch__pubchem"=FALSE,
            "switch__pubchem_props"=FALSE,
            "switch__gras"=FALSE,
            "switch__foodb_enzymes"=FALSE,
            "switch__foodb_flavors"=FALSE,
            "switch__foodb_foodcontent"=FALSE,
            "switch__foodb_healtheffects"=FALSE,
            "switch__foodb_ontology"=FALSE,
            "switch__epa_props"=FALSE
        )
        annotations <- chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode="")[["anno_invitrodb"]]
    }
}

#* Given a DSSTox substance ID, find neoplasticity data for the corresponding chemical in ToxRefDB data available in CBT.
#* @param dtxsid DSSTox substance ID
#* @get /toxrefdb/neoplasticity
#* @tag "Disease Annotations"
function(dtxsid = "") {
    # Get internal IDs for given chemical
    internal_id <- run_query(paste0("
        SELECT DISTINCT
            bc.epa_id,
            bc.dsstox_substance_id
        FROM
            base_chemicals bc
        WHERE
            bc.dsstox_substance_id = $1
    "), args=as.list(dtxsid))
    
    if(nrow(internal_id) > 0){
        
        input_all_annotations <- list(
            "switch__admet_binr"=FALSE,
            "switch__admet_catg"=FALSE,
            "switch__admet_cont"=FALSE,
            "switch__chembl"=FALSE,
            "switch__cpd"=FALSE,
            "switch__ctd_bioprocess"=FALSE,
            "switch__ctd_cellcomp"=FALSE,
            "switch__ctd_diseases"=FALSE,
            "switch__ctd_genes"=FALSE,
            "switch__ctd_molfunct"=FALSE,
            "switch__ctd_phenotypes"=FALSE,
            "switch__drugbank_atccodes"=FALSE,
            "switch__drugbank_carriers"=FALSE,
            "switch__drugbank_enzymes"=FALSE,
            "switch__drugbank_targets"=FALSE,
            "switch__drugbank_transporters"=FALSE,
            "switch__hmdb_biospecimenlocations"=FALSE,
            "switch__hmdb_cellularlocations"=FALSE,
            "switch__hmdb_diseases"=FALSE,
            "switch__hmdb_genes"=FALSE,
            "switch__hmdb_tissuelocations"=FALSE,
            "switch__ice"=FALSE,
            "switch__invitrodb"=FALSE,
            "switch__leadscope"=FALSE,
            "switch__ochem"=FALSE,
            "switch__opera_adme"=FALSE,
            "switch__opera_environmental_fate"=FALSE,
            "switch__opera_physicochem"=FALSE,
            "switch__opera_structural_properties"=FALSE,
            "switch__opera_toxicology_endpoints"=FALSE,
            "switch__saagar"=FALSE,
            "switch__t3db"=FALSE,
            "switch__toxrefdb_neoplastic"=TRUE,
            "switch__toxrefdb_nonneoplastic"=TRUE,
            "switch__pubchem"=FALSE,
            "switch__pubchem_props"=FALSE,
            "switch__gras"=FALSE,
            "switch__foodb_enzymes"=FALSE,
            "switch__foodb_flavors"=FALSE,
            "switch__foodb_foodcontent"=FALSE,
            "switch__foodb_healtheffects"=FALSE,
            "switch__foodb_ontology"=FALSE,
            "switch__epa_props"=FALSE
        )
        annotations <- chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode="")[c("anno_toxrefdb_neoplastic", "anno_toxrefdb_nonneoplastic")]
    }
}

#* Given a DSSTox substance ID, annotate using the Comparative Toxicogenomics Database (CTD) available in CBT. Contains annotations related to a chemical's biological processes, cellular components, molecular functions, associated diseases and phenotypes, and affected genes.
#* @param dtxsid DSSTox substance ID
#* @get /ctd
#* @tag "Gene Annotations"
#* @tag "Disease Annotations"
#* @tag "Biological Process Annotations"
function(dtxsid = "") {
    # Get internal IDs for given chemical
    internal_id <- run_query(paste0("
        SELECT DISTINCT
            bc.epa_id,
            bc.dsstox_substance_id
        FROM
            base_chemicals bc
        WHERE
            bc.dsstox_substance_id = $1
    "), args=as.list(dtxsid))
    
    if(nrow(internal_id) > 0){
        
        input_all_annotations <- list(
            "switch__admet_binr"=FALSE,
            "switch__admet_catg"=FALSE,
            "switch__admet_cont"=FALSE,
            "switch__chembl"=FALSE,
            "switch__cpd"=FALSE,
            "switch__ctd_bioprocess"=TRUE,
            "switch__ctd_cellcomp"=TRUE,
            "switch__ctd_diseases"=TRUE,
            "switch__ctd_genes"=TRUE,
            "switch__ctd_molfunct"=TRUE,
            "switch__ctd_phenotypes"=TRUE,
            "switch__drugbank_atccodes"=FALSE,
            "switch__drugbank_carriers"=FALSE,
            "switch__drugbank_enzymes"=FALSE,
            "switch__drugbank_targets"=FALSE,
            "switch__drugbank_transporters"=FALSE,
            "switch__hmdb_biospecimenlocations"=FALSE,
            "switch__hmdb_cellularlocations"=FALSE,
            "switch__hmdb_diseases"=FALSE,
            "switch__hmdb_genes"=FALSE,
            "switch__hmdb_tissuelocations"=FALSE,
            "switch__ice"=FALSE,
            "switch__invitrodb"=FALSE,
            "switch__leadscope"=FALSE,
            "switch__ochem"=FALSE,
            "switch__opera_adme"=FALSE,
            "switch__opera_environmental_fate"=FALSE,
            "switch__opera_physicochem"=FALSE,
            "switch__opera_structural_properties"=FALSE,
            "switch__opera_toxicology_endpoints"=FALSE,
            "switch__saagar"=FALSE,
            "switch__t3db"=FALSE,
            "switch__toxrefdb_neoplastic"=FALSE,
            "switch__toxrefdb_nonneoplastic"=FALSE,
            "switch__pubchem"=FALSE,
            "switch__pubchem_props"=FALSE,
            "switch__gras"=FALSE,
            "switch__foodb_enzymes"=FALSE,
            "switch__foodb_flavors"=FALSE,
            "switch__foodb_foodcontent"=FALSE,
            "switch__foodb_healtheffects"=FALSE,
            "switch__foodb_ontology"=FALSE,
            "switch__epa_props"=FALSE
        )
        annotations <- as.list(chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode=""))[c("anno_ctd_bioprocess", "anno_ctd_cellcomp", "anno_ctd_diseases", "anno_ctd_genes", "anno_ctd_molfunct", "anno_ctd_phenotypes")]
    }
}

#* Given a DSSTox substance ID, annotate using the Human Metabolome Database (HMDB) available in CBT. Contains annotations related to a chemical's biospecimen locations, cellular locations, tissue locations, genes, and diseases.
#* @param dtxsid DSSTox substance ID
#* @get /hmdb
#* @tag "Gene Annotations"
#* @tag "Disease Annotations"
#* @tag "Biological Process Annotations"
function(dtxsid = "") {
    # Get internal IDs for given chemical
    internal_id <- run_query(paste0("
        SELECT DISTINCT
            bc.epa_id,
            bc.dsstox_substance_id
        FROM
            base_chemicals bc
        WHERE
            bc.dsstox_substance_id = $1
    "), args=as.list(dtxsid))
    
    if(nrow(internal_id) > 0){
        
        input_all_annotations <- list(
            "switch__admet_binr"=FALSE,
            "switch__admet_catg"=FALSE,
            "switch__admet_cont"=FALSE,
            "switch__chembl"=FALSE,
            "switch__cpd"=FALSE,
            "switch__ctd_bioprocess"=FALSE,
            "switch__ctd_cellcomp"=FALSE,
            "switch__ctd_diseases"=FALSE,
            "switch__ctd_genes"=FALSE,
            "switch__ctd_molfunct"=FALSE,
            "switch__ctd_phenotypes"=FALSE,
            "switch__drugbank_atccodes"=FALSE,
            "switch__drugbank_carriers"=FALSE,
            "switch__drugbank_enzymes"=FALSE,
            "switch__drugbank_targets"=FALSE,
            "switch__drugbank_transporters"=FALSE,
            "switch__hmdb_biospecimenlocations"=TRUE,
            "switch__hmdb_cellularlocations"=TRUE,
            "switch__hmdb_diseases"=TRUE,
            "switch__hmdb_genes"=TRUE,
            "switch__hmdb_tissuelocations"=TRUE,
            "switch__ice"=FALSE,
            "switch__invitrodb"=FALSE,
            "switch__leadscope"=FALSE,
            "switch__ochem"=FALSE,
            "switch__opera_adme"=FALSE,
            "switch__opera_environmental_fate"=FALSE,
            "switch__opera_physicochem"=FALSE,
            "switch__opera_structural_properties"=FALSE,
            "switch__opera_toxicology_endpoints"=FALSE,
            "switch__saagar"=FALSE,
            "switch__t3db"=FALSE,
            "switch__toxrefdb_neoplastic"=FALSE,
            "switch__toxrefdb_nonneoplastic"=FALSE,
            "switch__pubchem"=FALSE,
            "switch__pubchem_props"=FALSE,
            "switch__gras"=FALSE,
            "switch__foodb_enzymes"=FALSE,
            "switch__foodb_flavors"=FALSE,
            "switch__foodb_foodcontent"=FALSE,
            "switch__foodb_healtheffects"=FALSE,
            "switch__foodb_ontology"=FALSE,
            "switch__epa_props"=FALSE
        )
        annotations <- as.list(chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode=""))[c("anno_hmdb_biospecimenlocations", "anno_hmdb_cellularlocations", "anno_hmdb_diseases", "anno_hmdb_genes", "anno_hmdb_tissuelocations")]
    }
}

#* Given a DSSTox substance ID, annotate using results from PASS predictive models available in CBT.
#* @param dtxsid DSSTox substance ID
#* @get /pass
#* @tag "Predictive Models"
function(dtxsid = "") {
    pass_pred <- run_query(paste0("
        SELECT DISTINCT
            ppv.pass_prediction_results::text
        FROM
            base_chemicals bc,
            base_chemical_compounds bcc,
            pass_chemical_predictions_vector ppv
        WHERE
            bc.dsstox_substance_id = $1
        AND bc.epa_id = bcc.epa_id
        AND bcc.dsstox_compound_id = ppv.dsstox_compound_id
    "), args=as.list(dtxsid))
    
    pass_models <- run_query(paste0("
        SELECT DISTINCT
            ppm.position,
            ppm.pass_name
        FROM
            pass_predictive_models ppm
        ORDER BY ppm.position ASC
    "))
    
    pass_pred <- gsub("\\[", "", pass_pred)
    pass_pred <- gsub("\\]", "", pass_pred)
    pass_pred <- unlist(str_split(pass_pred, ","))
    pass_pred <- as.numeric(pass_pred)
    names(pass_pred) <- pass_models$pass_name
    as.list(pass_pred)
}

#* Given a DSSTox substance ID, annotate using results from CPDat's QSUR predictive models available in CBT.
#* @param dtxsid DSSTox substance ID
#* @get /cpdat/qsur
#* @tag "Predictive Models"
function(dtxsid = "") {
    qsur_pred <- run_query(paste0("
        SELECT DISTINCT
            cqd.harmonized_function,
            cqd.probability
        FROM
            cpd_qsur_data cqd
        WHERE
            cqd.dtxsid = $1
    "), args=as.list(dtxsid))
    return(qsur_pred)
}


# #* Given a DSSTox substance ID, annotate using predictive models from Instem (Leadscope) available in CBT.
# #* @param dtxsid DSSTox substance ID
# #* @param positives only return positive models?
# #* @get /leadscope
# #* @tag "Predictive Models"
# #* @tag "Environmental Fate and Exposure Annotations"
# function(dtxsid = "", positives = FALSE) {
#     # Get internal IDs for given chemical
#     internal_id <- run_query(paste0("
#         SELECT DISTINCT
#             bc.epa_id,
#             bc.dsstox_substance_id
#         FROM
#             base_chemicals bc
#         WHERE
#             bc.dsstox_substance_id = $1
#     "), args=as.list(dtxsid))
#     
#     if(nrow(internal_id) > 0){
#         
#         input_all_annotations <- list(
#             "switch__admet_binr"=FALSE,
#             "switch__chembl"=FALSE,
#             "switch__cpd"=FALSE,
#             "switch__ctd_bioprocess"=FALSE,
#             "switch__ctd_cellcomp"=FALSE,
#             "switch__ctd_diseases"=FALSE,
#             "switch__ctd_genes"=FALSE,
#             "switch__ctd_molfunct"=FALSE,
#             "switch__ctd_phenotypes"=FALSE,
#             "switch__drugbank_atccodes"=FALSE,
#             "switch__drugbank_carriers"=FALSE,
#             "switch__drugbank_enzymes"=FALSE,
#             "switch__drugbank_targets"=FALSE,
#             "switch__drugbank_transporters"=FALSE,
#             "switch__hmdb_biospecimenlocations"=FALSE,
#             "switch__hmdb_cellularlocations"=FALSE,
#             "switch__hmdb_diseases"=FALSE,
#             "switch__hmdb_genes"=FALSE,
#             "switch__hmdb_tissuelocations"=FALSE,
#             "switch__ice"=FALSE,
#             "switch__invitrodb"=FALSE,
#             "switch__leadscope"=TRUE,
#             "switch__ochem"=FALSE,
#             "switch__opera_adme"=FALSE,
#             "switch__opera_environmental_fate"=FALSE,
#             "switch__opera_physicochem"=FALSE,
#             "switch__opera_structural_properties"=FALSE,
#             "switch__opera_toxicology_endpoints"=FALSE,
#             "switch__saagar"=FALSE,
#             "switch__t3db"=FALSE,
#             "switch__toxrefdb_neoplastic"=FALSE,
#             "switch__toxrefdb_nonneoplastic"=FALSE,
#             "switch__pubchem"=FALSE,
#             "switch__pubchem_props"=FALSE,
#             "switch__gras"=FALSE,
#             "switch__foodb_enzymes"=FALSE,
#             "switch__foodb_flavors"=FALSE,
#             "switch__foodb_foodcontent"=FALSE,
#             "switch__foodb_healtheffects"=FALSE,
#             "switch__foodb_ontology"=FALSE,
#             "switch__epa_props"=FALSE
#         )
#         annotations <- chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode="")
#         
#         annotations <- annotations[["anno_leadscope"]]
#         
#         if(as.booltype(positives) == TRUE){
#             annotations <- annotations %>% filter(interpretation == 1)
#             #print(annotations %>% filter(interpretation == 1))
#         }
#         
#         return(annotations)
#         
#     }
# }

# #* Given a DSSTox substance ID, annotate using predictive models from ADMET Predictor available in CBT.
# #* @param dtxsid DSSTox substance ID
# #* @param positives only return positive models?
# #* @get /admet
# #* @tag "Predictive Models"
# #* @tag "Environmental Fate and Exposure Annotations"
# #* @tag "Transporter and Pathway Annotations"
# function(dtxsid = "", positives = FALSE) {
#     # Get internal IDs for given chemical
#     internal_id <- run_query(paste0("
#         SELECT DISTINCT
#             bc.epa_id,
#             bc.dsstox_substance_id
#         FROM
#             base_chemicals bc
#         WHERE
#             bc.dsstox_substance_id = $1
#     "), args=as.list(dtxsid))
#     
#     if(nrow(internal_id) > 0){
#         
#         input_all_annotations <- list(
#             "switch__admet_binr"=TRUE,
#             "switch__chembl"=FALSE,
#             "switch__cpd"=FALSE,
#             "switch__ctd_bioprocess"=FALSE,
#             "switch__ctd_cellcomp"=FALSE,
#             "switch__ctd_diseases"=FALSE,
#             "switch__ctd_genes"=FALSE,
#             "switch__ctd_molfunct"=FALSE,
#             "switch__ctd_phenotypes"=FALSE,
#             "switch__drugbank_atccodes"=FALSE,
#             "switch__drugbank_carriers"=FALSE,
#             "switch__drugbank_enzymes"=FALSE,
#             "switch__drugbank_targets"=FALSE,
#             "switch__drugbank_transporters"=FALSE,
#             "switch__hmdb_biospecimenlocations"=FALSE,
#             "switch__hmdb_cellularlocations"=FALSE,
#             "switch__hmdb_diseases"=FALSE,
#             "switch__hmdb_genes"=FALSE,
#             "switch__hmdb_tissuelocations"=FALSE,
#             "switch__ice"=FALSE,
#             "switch__invitrodb"=FALSE,
#             "switch__leadscope"=FALSE,
#             "switch__ochem"=FALSE,
#             "switch__opera_adme"=FALSE,
#             "switch__opera_environmental_fate"=FALSE,
#             "switch__opera_physicochem"=FALSE,
#             "switch__opera_structural_properties"=FALSE,
#             "switch__opera_toxicology_endpoints"=FALSE,
#             "switch__saagar"=FALSE,
#             "switch__t3db"=FALSE,
#             "switch__toxrefdb_neoplastic"=FALSE,
#             "switch__toxrefdb_nonneoplastic"=FALSE,
#             "switch__pubchem"=FALSE,
#             "switch__pubchem_props"=FALSE,
#             "switch__gras"=FALSE,
#             "switch__foodb_enzymes"=FALSE,
#             "switch__foodb_flavors"=FALSE,
#             "switch__foodb_foodcontent"=FALSE,
#             "switch__foodb_healtheffects"=FALSE,
#             "switch__foodb_ontology"=FALSE,
#             "switch__epa_props"=FALSE
#         )
#         annotations <- chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode="")
#         
#         annotations <- annotations[["anno_admet_binr"]]
#         
#         if(as.booltype(positives) == TRUE){
#             annotations <- annotations %>% filter(interpretation == 1)
#             #print(annotations %>% filter(interpretation == 1))
#         }
#         
#         return(annotations)
#         
#     }
# }

#* Given a DSSTox substance ID, annotate using the DrugBank database available in CBT. Contains information regarding a chemical's usage in drugs and their targets.
#* @param dtxsid DSSTox substance ID
#* @get /drugbank
#* @tag "Drug Annotations"
#* @tag "Commercial and Industrial Usage"
function(dtxsid = "") {
    # Get internal IDs for given chemical
    internal_id <- run_query(paste0("
        SELECT DISTINCT
            bc.epa_id,
            bc.dsstox_substance_id
        FROM
            base_chemicals bc
        WHERE
            bc.dsstox_substance_id = $1
    "), args=as.list(dtxsid))
    
    if(nrow(internal_id) > 0){
        
        input_all_annotations <- list(
            "switch__admet_binr"=FALSE,
            "switch__admet_catg"=FALSE,
            "switch__admet_cont"=FALSE,
            "switch__chembl"=FALSE,
            "switch__cpd"=FALSE,
            "switch__ctd_bioprocess"=FALSE,
            "switch__ctd_cellcomp"=FALSE,
            "switch__ctd_diseases"=FALSE,
            "switch__ctd_genes"=FALSE,
            "switch__ctd_molfunct"=FALSE,
            "switch__ctd_phenotypes"=FALSE,
            "switch__drugbank_atccodes"=TRUE,
            "switch__drugbank_carriers"=TRUE,
            "switch__drugbank_enzymes"=TRUE,
            "switch__drugbank_targets"=TRUE,
            "switch__drugbank_transporters"=TRUE,
            "switch__hmdb_biospecimenlocations"=FALSE,
            "switch__hmdb_cellularlocations"=FALSE,
            "switch__hmdb_diseases"=FALSE,
            "switch__hmdb_genes"=FALSE,
            "switch__hmdb_tissuelocations"=FALSE,
            "switch__ice"=FALSE,
            "switch__invitrodb"=FALSE,
            "switch__leadscope"=FALSE,
            "switch__ochem"=FALSE,
            "switch__opera_adme"=FALSE,
            "switch__opera_environmental_fate"=FALSE,
            "switch__opera_physicochem"=FALSE,
            "switch__opera_structural_properties"=FALSE,
            "switch__opera_toxicology_endpoints"=FALSE,
            "switch__saagar"=FALSE,
            "switch__t3db"=FALSE,
            "switch__toxrefdb_neoplastic"=FALSE,
            "switch__toxrefdb_nonneoplastic"=FALSE,
            "switch__pubchem"=FALSE,
            "switch__pubchem_props"=FALSE,
            "switch__gras"=FALSE,
            "switch__foodb_enzymes"=FALSE,
            "switch__foodb_flavors"=FALSE,
            "switch__foodb_foodcontent"=FALSE,
            "switch__foodb_healtheffects"=FALSE,
            "switch__foodb_ontology"=FALSE,
            "switch__epa_props"=FALSE
        )
        annotations <- chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode="")[c("anno_drugbank_atccodes", "anno_drugbank_carriers", "anno_drugbank_enzymes", "anno_drugbank_targets", "anno_drugbank_transporters")]
        
        return(annotations)
        
    }
}

# #* Given a DSSTox substance ID, annotate using the predicted metabolites from ADMET Predictor available in CBT.
# #* @param dtxsid DSSTox substance ID
# #* @param enzyme cyp or ugt or "both"
# #* @param maxlevel maximum level of metabolism
# #* @get /metabolites
# #* @tag "Predictive Models"
# #* @tag "Transporter and Pathway Annotations"
# function(dtxsid = "", enzyme = "", maxlevel = 1) {
#     # Get internal IDs for given chemical
#     internal_id <- run_query(paste0("
#         SELECT DISTINCT
#             bc.epa_id,
#             bc.dsstox_substance_id
#         FROM
#             base_chemicals bc
#         WHERE
#             bc.dsstox_substance_id = $1
#     "), args=as.list(dtxsid))
#     
#     if(nrow(internal_id) > 0){
#         metabolites_full <- list()
#         if(enzyme == "both"){
#             metabolites_full_cyp <- lapply(seq_len(maxlevel), function(x){
#                 metabolites <- metabolites_query("cyp", x, as.integer(internal_id$epa_id))
#             })
#             metabolites_full_cyp <- rbindlist(metabolites_full_cyp)
#             metabolites_full_ugt <- lapply(seq_len(maxlevel), function(x){
#                 metabolites <- metabolites_query("ugt", x, as.integer(internal_id$epa_id))
#             })
#             metabolites_full_ugt <- rbindlist(metabolites_full_ugt)
#             metabolites_full <- list(metabolites_full_cyp, metabolites_full_ugt)
#         } else {
#             metabolites_full <- lapply(seq_len(maxlevel), function(x){
#                 metabolites <- metabolites_query(enzyme, x, as.integer(internal_id$epa_id))
#             })
#         }
#         metabolites_full <- rbindlist(metabolites_full, fill=TRUE)
#         return(metabolites_full)
#     }
# }

#* Given a SMILES string, check for (OChem, ChEMBL, Saagar) structural alerts and notable substructures.
#* @param smiles SMILES string
#* @get /alerts
#* @tag "Substructure Search"
function(smiles = "") {
    # OChem
    alerts_ochem <- run_query(paste0("
        SELECT DISTINCT
            osa.smarts AS alert,
            'ochem' AS source
        FROM
            ochem_structural_alerts_mols osa
        WHERE
            mol_from_smiles($1)@> osa.m
    "), args=as.list(smiles))
    
    # ChEMBL
    alerts_chembl <- run_query(paste0("
        SELECT DISTINCT
            csa.alert_name AS alert,
            'chembl' AS source
        FROM
            chembl_structural_alerts csa
        WHERE
            mol_from_smiles($1)@> mol_from_smarts(csa.smarts)
    "), args=as.list(smiles))
    
    # Saagar
    alerts_saagar <- run_query(paste0("
        SELECT DISTINCT
            ssa.alert_name AS alert,
            'saagar' AS source
        FROM
            saagar_structural_alerts_mols ssa
        WHERE
            mol_from_smiles($1)@> ssa.m
    "), args=as.list(smiles))
    
    alerts <- rbindlist(list(alerts_ochem, alerts_chembl, alerts_saagar))
    return(alerts)
}


#* Given a DSSTox substance ID, annotate with chemical properties as detailed in PubChem that are available in CBT.
#* @param dtxsid DSSTox substance ID
#* @get /pubchem/properties
#* @tag "Chemical Properties"
function(dtxsid = "") {
    props <- run_query(paste0("
        SELECT DISTINCT
            pa.pubchem_attribute,
            pa.pubchem_value
        FROM
            base_chemicals bc,
            base_chemical_to_pubchem_cid bpc,
            pubchem_attributes pa
        WHERE
            bc.dsstox_substance_id = $1
        AND bc.epa_id = bpc.epa_id
        AND bpc.pubchem_cid = pa.cid
    "), args=as.list(dtxsid))
    return(props)
}

#* Given a DSSTox substance ID, annotate with chemical properties as detailed in the EPA CompTox Chemicals Dashboard that are available in CBT.
#* @param dtxsid DSSTox substance ID
#* @get /epa/properties
#* @tag "Chemical Properties"
function(dtxsid = "") {
    props <- run_query(paste0("
        SELECT DISTINCT
            ep.epa_attribute,
            ep.epa_attribute_value AS epa_value
        FROM
            base_chemicals bc,
            epa_properties ep
        WHERE
            bc.dsstox_substance_id = $1
        AND bc.epa_id = ep.epa_id
    "), args=as.list(dtxsid))
    return(props)
}


#* Given a DSSTox substance ID, return bioassays in PubChem that are available in CBT.
#* @param dtxsid DSSTox substance ID
#* @get /pubchem/bioassays
#* @tag "Biological Process Annotations"
#* @tag "Assays"
function(dtxsid = "") {
    assays <- run_query(paste0("
        SELECT DISTINCT
            pb.bioassay_name,
            pb.deposit_date,
            pb.modify_date,
            pb.source_name,
            pb.source_id,
            pb.substance_type,
            pb.outcome_type,
            pb.project_category,
            pb.bioassay_group,
            pb.bioassay_types,
            pbm.activity_outcome,
            pbm.activity_name,
            pbm.activity_qualifier,
            pbm.activity_value,
            pbm.protein_accession
        FROM
            base_chemicals bc,
            base_chemical_to_pubchem_cid bpc,
            pubchem_bioactivity_mapping pbm,
            pubchem_bioassays pb
        WHERE
            bc.dsstox_substance_id = $1
        AND bc.epa_id = bpc.epa_id
        AND bpc.pubchem_cid = pbm.cid
        AND pbm.aid = pb.aid
    "), args=as.list(dtxsid))
    return(assays)
}

#* Given a DSSTox substance ID, return test articles from the NTP that are available in CBT.
#* @param dtxsid DSSTox substance ID
#* @get /ntp/test_articles
#* @tag "Assays"
function(dtxsid = "") {
    assays <- run_query(paste0("
        SELECT DISTINCT
            bcta.testarticle,
            bcta.study_area,
            bcta.study_number,
            nta.web_url
        FROM
            base_chemicals bc,
            base_chemical_test_articles bcta,
            ntp_test_articles nta
        WHERE
            bc.dsstox_substance_id = $1
        AND bc.epa_id = bcta.epa_id
        AND bcta.ntpid = nta.ntpid
    "), args=as.list(dtxsid))
    return(assays)
}


#* Given a DSSTox substance ID, annotate with a chemical's upper 95th percentile SEEM3 exposure available in CBT.
#* @param dtxsid DSSTox substance ID
#* @get /seem3
#* @tag "Environmental Fate and Exposure Annotations"
function(dtxsid = "") {
    # Get internal IDs for given chemical
    internal_id <- run_query(paste0("
        SELECT DISTINCT
            bc.epa_id,
            bc.dsstox_substance_id
        FROM
            base_chemicals bc
        WHERE
            bc.dsstox_substance_id = $1
    "), args=as.list(dtxsid))
    
    if(nrow(internal_id) > 0){
        
        # SEEM3
        exp_seem3 <- run_query(paste0("
            SELECT DISTINCT
                bcs.data_value AS annotations,
                'seem3 u95' AS source
            FROM
                base_chemical_seem3 bcs
            WHERE
                bcs.epa_id::text = $1
            AND bcs.data_name = 'seem3u95'
        "), args=as.list(as.character(internal_id$epa_id)))
        
        
        exp <- rbindlist(list(exp_seem3))
        return(exp)
    }
}

#* Given a DSSTox substance ID, annotate with a chemical's Cmax concentration data.
#* @param dtxsid DSSTox substance ID
#* @get /cmax
#* @tag "Environmental Fate and Exposure Annotations"
function(dtxsid = "") {
    cmax <- run_query(paste0("
        SELECT DISTINCT
            bcc.mw,
            bcc.fu,
            bcc.clint,
            bcc.pka_accept,
            bcc.pka_donor,
            bcc.logp,
            bcc.hl,
            bcc.cmax_calc
        FROM
            base_chemicals bc,
            base_chemical_cmax bcc
        WHERE
            bc.dsstox_substance_id = $1
        AND bc.epa_id = bcc.epa_id
    "), args=as.list(dtxsid))
    return(cmax)
}

#* Given a DSSTox substance ID, annotate with a chemical's Cmax concentration data.
#* @param dtxsid DSSTox substance ID
#* @get /superfund
#* @tag "Environmental Fate and Exposure Annotations"
function(dtxsid = "") {
    sfund <- run_query(paste0("
        SELECT DISTINCT
            sdd.rmedia_desc,
            sdd.cc_max_conc_value_nmbr,
            sdd.cc_cont_of_cncrn_flag,
            sdd.rconc_units_desc,
            ssi.site_name,
            ssi.city,
            ssi.zipcode,
            ssi.latitude,
            ssi.longitude
        FROM
            superfund_to_base_chemicals sbc,
            superfund_data_dsstox sdd,
            superfund_site_information ssi
        WHERE
            sbc.dsstox_substance_id = $1
        AND sbc.epasfund_id = sdd.epasfund_id
        AND sdd.sfundsite_id = ssi.sfundsite_id
    "), args=as.list(dtxsid))
    return(sfund)
}



# #* Plot a histogram
# #* @serializer png
# #* @get /plot
# function() {
#     rand <- rnorm(100)
#     hist(rand)
# }
# 
# #* Return the sum of two numbers
# #* @param a The first number to add
# #* @param b The second number to add
# #* @post /sum
# function(a, b) {
#     as.numeric(a) + as.numeric(b)
# }
# 
# # Programmatically alter your API
# #* @plumber
# function(pr) {
#     pr %>%
#         # Overwrite the default serializer to return unboxed JSON
#         pr_set_serializer(serializer_unboxed_json())
# }
