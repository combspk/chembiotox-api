#
# This is a Plumber API. You can run the API by clicking
# the 'Run API' button above.
#
# Find out more about building APIs with Plumber here:
#
#    https://www.rplumber.io/
#
library(connectapi)
library(plumber)
library(bit64)
library(data.table)
library(DBI)
library(dplyr)
library(httr)
library(jsonlite)
library(openssl)
library(openxlsx)
library(plotly)
library(pool)
library(reshape2)
library(rsvg)
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
#* @apiDescription This is an API developed using Plumber for R designed for interfacing with the ChemBioTox (CBT) database. CBT is a comprehensive chemical database built on Postgres that contains property and toxicological annotations for over 1 million chemicals of interest to the EPA. Data within CBT is sourced from various publicly available databases (i.e., CTD), calculated using bioinformatics tools (i.e., RDKit), or produced from proprietary technology from Instem and Simulations Plus (not available in the public version of the API: you must log in using NIH credentials to access these endpoints). <br><br> Datasources and tools used in the development of the ChemBioTox (CBT) database include:<br> <ul> <li> ChEMBL v3.4 - Zdrazil B, et al. The ChEMBL Database in 2023: a drug discovery platform spanning multiple bioactivity data types and time periods. Nucleic Acids Research, Volume 52, Issue D1, 5 January 2024, Pages D1180–D1192, <a href="https://doi.org/10.1093/nar/gkad1004">https://doi.org/10.1093/nar/gkad1004</a> </li><li> Chemical and Products Database (CPDat) - Dionisio K, Philips K, Price P, et al. The Chemical and Products Database, a resource for exposure-relevant data on chemicals in consumer products. Sci Data 5, 180125 (2018). <a href="https://doi.org/10.1038/sdata.2018.125"> https://doi.org/10.1038/sdata.2018.125</a> </li><li> Comparative Toxicogenomics Database (CTD) - Davis AP, Wiegers TC, Johnson RJ, Sciaky D, Wiegers J, Mattingly CJ. Comparative Toxicogenomics Database (CTD): update 2023. Nucleic Acids Res. 2022 Sep 28. <a href="https://ctdbase.org/">https://ctdbase.org/</a> </li><li> DrugBank Online - Knox C, et al. DrugBank 6.0: the DrugBank Knowledgebase for 2024. Nucleic Acids Res. 2024 Jan 5;52(D1):D1265-D1275. <a href="https://doi.org/10.1093/nar/gkad976">https://doi.org/10.1093/nar/gkad976</a>, <a href="https://go.drugbank.com/">https://go.drugbank.com/</a> </li><li> EPA CompTox Chemicals Dashboard v2.4.1 - Williams AJ, Grulke CM, Edwards J, et al. The CompTox Chemistry Dashboard: a community data resource for environmental chemistry. J Cheminform 9, 61 (2017). <a href="https://doi.org/10.1186/s13321-017-0247-6">https://doi.org/10.1186/s13321-017-0247-6</a>, <a href="https://comptox.epa.gov/dashboard/">https://comptox.epa.gov/dashboard/</a> </li><li> FooDB v1.0 - https://foodb.ca</li><li> Generally Recognized As Safe (GRAS) Substances SCOGS Database - US FDA <a href="https://www.cfsanappsexternal.fda.gov/scripts/fdcc/?set=SCOGS">https://www.cfsanappsexternal.fda.gov/scripts/fdcc/?set=SCOGS</a>, <a href="https://www.fda.gov/food/generally-recognized-safe-gras/gras-substances-scogs-database">https://www.fda.gov/food/generally-recognized-safe-gras/gras-substances-scogs-database</a> </li><li> Human Metabolome Database (HMDB) - Wishart DS, Guo A, Oler E, Wang F, Anjum A, Peters H, Dizon R, Sayeeda Z, Tian S, Lee BL, Berjanskii M, Mah R, Yamamoto M, Jovel J, Torres-Calzada C, Hiebert-Giesbrecht M, Lui VW, Varshavi D, Varshavi D, Allen D, Arndt D, Khetarpal N, Sivakumaran A, Harford K, Sanford S, Yee K, Cao X, Budinski Z, Liigand J, Zhang L, Zheng J, Mandal R, Karu N, Dambrova M, Schiöth HB, Greiner R, Gautam V. HMDB 5.0: the Human Metabolome Database for 2022. Nucleic Acids Res. 2022 Jan 7;50(D1):D622-D631. <a href="https://doi.org/10.1093/nar/gkab1062">https://doi.org/10.1093/nar/gkab1062</a>. PMID: 34986597; PMCID: PMC8728138.</li><li> InVitroDB v3.5 - <a href="https://www.epa.gov/comptox-tools/exploring-toxcast-data">https://www.epa.gov/comptox-tools/exploring-toxcast-data</a> </li><li> OCHEM - Sushko I, et al. Online chemical modeling environment (OCHEM): web platform for data storage, model development and publishing of chemical information. J Comput Aided Mol Des. 2011; 25(6):533-54. <a href="https://doi.org/10.1007/s10822-011-9440-2">https://doi.org/10.1007/s10822-011-9440-2</a>, <a href="https://ochem.eu/home/show.do">https://ochem.eu/home/show.do</a> </li><li> PubChem - Kim S, Chen J, Cheng T, et al. PubChem 2023 update. Nucleic Acids Res. 2023;51(D1):D1373–D1380. <a href="https://doi.org/10.1093/nar/gkac956">https://doi.org/10.1093/nar/gkac956</a>, <a href="https://pubchem.ncbi.nlm.nih.gov/ ">https://pubchem.ncbi.nlm.nih.gov/</a> </li><li> RDKit - <a href="https://github.com/rdkit/rdkit">https://github.com/rdkit/rdkit</a> </li><li> Toxin and Toxin Target Database (T3DB) - Wishart D, et al. T3DB: the toxic exposome database. Nucleic Acids Res. 2015 Jan;43(Database issue):D928-34. <a href="ttps://doi.org/10.1093/nar/gku1004">https://doi.org/10.1093/nar/gku1004</a>; Lim E, et al. T3DB: a comprehensively annotated database of common toxins and their targets. Nucleic Acids Res. 2010 Jan 38(Database issue):D781-6. <a href="https://doi.org/10.1093/nar/gkp934">https://doi.org/10.1093/nar/gkp934</a>; <a href="http://www.t3db.ca/">http://www.t3db.ca/</a> </li><li> ToxRefDB - Martin MT and Judson R. ToxRefDB - Release user-friendly web-based tool for mining ToxRefDB. U.S. Environmental Protection Agency, Washington, DC, 2010. <a href="https://cfpub.epa.gov/si/si_public_record_report.cfm?Lab=NCCT&dirEntryId=227139">https://cfpub.epa.gov/si/si_public_record_report.cfm?Lab=NCCT&dirEntryId=227139</a>, <a href="https://www.epa.gov/comptox-tools/exploring-toxcast-data">https://www.epa.gov/comptox-tools/exploring-toxcast-data</a> </li> <li>Generalized Read-Across (GenRA) - Patlewicz G and Shah I. Towards systematic read-across using Generalised Read-Across (GenRA). Computational Toxicology 2023, Vol 25. <a href="https://doi.org/10.1016/j.comtox.2022.100258">https://doi.org/10.1016/j.comtox.2022.100258</a></li></ul><br>The API's GitHub repository can be accessed <a href="https://github.com/combspk/chembiotox-api">here</a>.
#* @apiVersion 0.1.0


# Authentication Check - from https://docs.posit.co/connect/user/plumber/
# Returns a list containing "user" and "groups" information 
# populated by incoming request data.
getUserMetadata <- function(req) {
    rawUserData <- req[["HTTP_RSTUDIO_CONNECT_CREDENTIALS"]]
    if (!is.null(rawUserData)) {
        jsonlite::fromJSON(rawUserData)
    } else {
        list()
    }
}

#* Get a list of all chemicals supported by ChemBioTox
#* @get /lists/chemicals/chembiotox
#* @tag "Lists"
function(res, req) {
    # Get internal IDs for given chemical
    chemicals <- run_query(paste0("
        SELECT DISTINCT
            bc.dsstox_substance_id,
            bcc.dsstox_compound_id,
            bcc.preferred_name,
            bcc.casrn,
            bcc.inchi,
            bcc.inchikey
        FROM
            base_chemicals bc
        LEFT JOIN base_chemical_compounds bcc
            ON bc.epa_id = bcc.epa_id
    "))
    fwrite(chemicals, "lists/lists_chemicals_chembiotox.csv", sep=",")
    return(FALSE)
    # include_file("lists/lists_chemicals_chembiotox.csv", res, "text/csv")
}

#* Get a list of all chemicals supported by DrugMatrix
#* @get /lists/chemicals/drugmatrix
#* @tag "Lists"
function(res, req) {
    # # Get internal IDs for given chemical
    # chemicals <- run_query(paste0("
    #     SELECT
    #         dbc.dsstox_substance_id,
    #         dbc.compound,
    #         dbc.chemical_name
    #     FROM
    #         drugmatrix_to_base_chemicals dbc
    #     LEFT JOIN base_chemicals bc
    #         ON dbc.dsstox_substance_id = bc.dsstox_substance_id
    #     LEFT JOIN base_chemical_compounds bcc
    #         ON bc.epa_id = bcc.epa_id
    # "))
    # fwrite(chemicals, "lists/lists_chemicals_drugmatrix.csv", sep=",")
    
    include_file("lists/lists_chemicals_drugmatrix.csv", res, "text/csv")
}

#* Get a list of all chemicals supported by DrugBank
#* @get /lists/chemicals/drugbank
#* @tag "Lists"
function(res, req) {
    # Get internal IDs for given chemical
    # chemicals <- run_query(paste0("
    #     SELECT
    #         dbc.dsstox_substance_id,
    #         dbc.drugbank_id,
    #         bcc.preferred_name
    #     FROM
    #         drugbank_to_base_chemicals dbc
    #     LEFT JOIN base_chemical_compounds bcc
    #         ON dbc.epa_id = bcc.epa_id
    # "))
    # fwrite(chemicals, "lists/lists_chemicals_drugbank.csv", sep=",")
    
    include_file("lists/lists_chemicals_drugbank.csv", res, "text/csv")
}

#* Get a list of all chemicals supported by the T3DB
#* @get /lists/chemicals/t3db
#* @tag "Lists"
function(res, req) {
    # Get internal IDs for given chemical
    # chemicals <- run_query(paste0("
    #     SELECT
    #         ttbc.dsstox_substance_id,
    #         tcc.chemical_name,
    #         tcc.casrn,
    #         tcc.smiles,
    #         tcc.inchi,
    #         tcc.inchikey,
    #         tcc.pubchem_id,
    #         tcc.chemspider_id,
    #         tcc.kegg_id
    #     FROM
    #         t3db_curated_chemicals tcc
    #     LEFT JOIN t3db_to_base_chemicals ttbc
    #         ON tcc.t3db_id = ttbc.t3db_id
    # "))
    # fwrite(chemicals, "lists/lists_chemicals_t3db.csv", sep=",")
    
    include_file("lists/lists_chemicals_t3db.csv", res, "text/csv")
}

#* Get a list of all chemicals supported by the T3DB
#* @get /lists/chemicals/t3db
#* @tag "Lists"
function(res, req) {
    # Get internal IDs for given chemical
    # chemicals <- run_query(paste0("
    #     SELECT
    #         ttbc.dsstox_substance_id,
    #         tcc.chemical_name,
    #         tcc.casrn,
    #         tcc.smiles,
    #         tcc.inchi,
    #         tcc.inchikey,
    #         tcc.pubchem_id,
    #         tcc.chemspider_id,
    #         tcc.kegg_id
    #     FROM
    #         t3db_curated_chemicals tcc
    #     LEFT JOIN t3db_to_base_chemicals ttbc
    #         ON tcc.t3db_id = ttbc.t3db_id
    # "))
    # fwrite(chemicals, "lists/lists_chemicals_t3db.csv", sep=",")
    
    include_file("lists/lists_chemicals_t3db.csv", res, "text/csv")
}

#* Get a list of all chemicals supported by Tox21
#* @get /lists/chemicals/tox21
#* @tag "Lists"
function(res, req) {
    # Get internal IDs for given chemical
    # chemicals <- run_query(paste0("
    #     SELECT DISTINCT
    #         bc.dsstox_substance_id,
    #         bcc.dsstox_compound_id,
    #         bcc.preferred_name,
    #         bcc.casrn,
    #         bcc.inchi,
    #         bcc.inchikey
    #     FROM
    #         base_chemicals bc
    #     LEFT JOIN base_chemical_compounds bcc
    #         ON bc.epa_id = bcc.epa_id
    #     WHERE bc.dsstox_substance_id IN (SELECT DISTINCT dsstox_substance_id FROM tox21_to_base_chemicals)
    # "))
    # fwrite(chemicals, "lists/lists_chemicals_tox21.csv", sep=",")
    include_file("lists/lists_chemicals_tox21.csv", res, "text/csv")
}

#* Get a list of all chemicals supported by PubChem
#* @get /lists/chemicals/pubchem
#* @tag "Lists"
function(res, req) {
    # Get internal IDs for given chemical
    # chemicals <- run_query(paste0("
    #     SELECT DISTINCT
    #         bcpc.pubchem_cid,
    #         bc.dsstox_substance_id,
    #         bcc.dsstox_compound_id,
    #         bcc.preferred_name,
    #         bcc.casrn,
    #         bcc.inchi,
    #         bcc.inchikey
    #     FROM
    #         base_chemical_to_pubchem_cid bcpc,
    #         base_chemicals bc
    #     LEFT JOIN base_chemical_compounds bcc
    #         ON bc.epa_id = bcc.epa_id
    #     WHERE bc.epa_id = bcpc.epa_id
    #
    # "))
    # fwrite(chemicals, "lists/lists_chemicals_pubchem.csv", sep=",")
    include_file("lists/lists_chemicals_pubchem.csv", res, "text/csv")
}

#* Given a DSSTox substance ID, return vendors that sell the chemical
#* @param dtxsid DSSTox substance ID
#* @get /availability/vendors
#* @tag "Availability"
function(res, req, dtxsid = "", pages=1) {
    # Get internal IDs for given chemical
    internal_id <- run_query(paste0("
         SELECT DISTINCT
            bc.dsstox_substance_id,
            bcc.preferred_name,
            bcc.casrn,
            pca.source_name,
            pca.source_record_url AS url
        FROM
            base_chemical_to_pubchem_cid bpc,
            pubchem_chemical_availability pca,
            base_chemicals bc,
            base_chemical_compounds bcc
        WHERE
            bc.dsstox_substance_id = $1
        AND bc.epa_id = bcc.epa_id
        AND bc.epa_id = bpc.epa_id
        AND bpc.pubchem_cid = pca.pubchem_cid
        LIMIT $2
    "), args=list(dtxsid, as.numeric(pages)*10))
}


#* Given a chemical name, return the DSSTox substance ID.
#* @param name chemical name
#* @get /name2dtxsid
#* @tag "Identifier Lookup"
function(res, req, name = "") {
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
#* @param casrn chemical name
#* @get /casrn2dtxsid
#* @tag "Identifier Lookup"
function(res, req, casrn = "") {
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

#* Given a chemical name, return the DSSTox substance ID. Takes into account synonyms when mapping.
#* @param name chemical name
#* @get /synonym2dtxsid
#* @tag "Identifier Lookup"
function(res, req, name = "") {
    # Get internal IDs for given chemical
    internal_id <- run_query(paste0("
        SELECT DISTINCT
            bc.dsstox_substance_id
        FROM
            base_chemicals bc,
            base_chemical_to_pubchem_synonyms bpcs
        WHERE
            UPPER(bpcs.synonym) = UPPER($1)
        AND bc.epa_id = bpcs.epa_id
        
    "), args=as.list(name))
}

#* Given a SMILES string, return the DSSTox substance ID. If the given SMILES is not in the database, will return the closest match using rdkit fingerprint-based structural similarity.
#* @param smiles SMILES string
#* @get /smiles2dtxsid
#* @tag "Identifier Lookup"
function(res, req, smiles = "") {
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
function(res, req, dtxsid = "") {
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
function(res, req, dtxsid = "") {
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



#* Given a DSSTox substance ID and fingerprint type, return a list of structurally similar chemicals. Supported fingerprint types are: "atompair", "avalon", "maccskeys", "morgan", "rdkit", "rdklayered", "rdkpattern", "topologicaltorsion"
#* @param smiles SMILES string
#* @param fp Fingerprint type. Supported fingerprint types are: "atompair", "avalon", "maccskeys", "morgan", "rdkit", "rdklayered", "rdkpattern", "topologicaltorsion". Default "morgan"
#* @param threshold Similarity threshold for chemicals to return. Only chemicals with at least this level of structural similarity will be returned. Default 0.5.
#* @get /similarity/structural
#* @tag "Similarity Lookup"
function(res, req, smiles="", fp="morgan", threshold=0.5, n=10, exact=TRUE) {
    
    exact <- as.logical(exact) # have to do this because Swagger converts boolean to string for some reason
    
    fp_types <- c("atompair", "avalon", "maccskeys", "morgan", "rdkit", "rdklayered", "rdkpattern", "topologicaltorsion")
    if(!(fp %in% fp_types)){
        return("Error: invalid fingerprint type provided.")
    }
    
    # Get internal IDs for given chemical
    internal_id <- run_query(paste0("
        SELECT DISTINCT
            bc.dsstox_substance_id,
            bcc.preferred_name,
            sim.similarity
        FROM
            base_chemicals bc,
            base_chemical_compounds bcc,
            base_chemical_to_smiles bcs,
            get_", fp, "_fp_neighbors($1) sim
        WHERE
            sim.smi_id = bcs.smi_id
            AND bc.epa_id = bcc.epa_id
            AND bcs.epa_id = bc.epa_id
            AND sim.similarity >= $2
            ", if(exact==FALSE) paste0("AND sim.similarity < 1"), "
        ORDER BY sim.similarity DESC
        LIMIT $3
        
    "), args=list(smiles, threshold, n))
}


#* Given a DSSTox substance ID, return a chemical binary vector fingerprint generated with RDKit. Supported types are: "atompair", "avalon", "maccskeys", "morgan", "rdkit", "rdklayer", "rdkpattern", "topologicaltorsion", default is "rdkit".
#* @param dtxsid DSSTox substance ID
#* @param mode Fingerprint mode. Supported modes are: "atompair", "avalon", "maccskeys", "morgan", "rdkit", "rdklayer", "rdkpattern", "topologicaltorsion"
#* @get /fingerprint
#* @tag "Identifier Lookup"
function(res, req, dtxsid="", mode="rdkit") {
    
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

#* Given a DSSTox substance ID, return chemical properties from GenRA. NOTE: This endpoint is still in development and responses may not be correct.
#* @param dtxsid DSSTox substance ID
#* @get /genra/properties
#* @tag "Chemical Properties"
function(res, req, dtxsid="", pages=1) {
    # Get internal IDs for given chemical
    prop <- run_query(paste0("
        SELECT DISTINCT
            bach.dsstox_substance_id,
            bcco.dsstox_compound_id,
            bcco.preferred_name,
            geca.genra_catagory,
            geca.genra_catagory_name,
            bcgr.genra_result
        FROM
            base_chemicals bach,
            base_chemical_compounds bcco,
            base_chemical_genra_results bcgr,
            genra_catagories geca
        WHERE
            bcco.epa_id = bach.epa_id
        AND bcgr.epa_id = bach.epa_id
        AND bcgr.cmpd_id = bcco.cmpd_id
        AND bcgr.gcat_id = geca.gcat_id
        AND bach.dsstox_substance_id = $1
        LIMIT $2
    "), args=list(dtxsid, as.numeric(pages)*10))
}

#* Given a DSSTox substance ID, annotate using the Chemical and Products Database (CPDat) available in CBT. Contains information about a chemical's usage in commercially-available products.
#* @param dtxsid DSSTox substance ID
#* @get /cpdat
#* @tag "Commercial and Industrial Usage"
#* @tag "Environmental Fate and Exposure Annotations"
function(res, req, dtxsid = "", pages=1) {
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
            "switch__cpd"
        )
        annotations <- chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode="")[["anno_cpd"]]
        if(nrow(annotations) > as.numeric(pages)*10){
            annotations <- annotations[seq_len(as.numeric(pages)*10), ]
        }
        return(annotations)
    }
}


#* Given a DSSTox substance ID, annotate using QSUR model predictions from the Chemical and Products Database (CPDat) available in CBT.
#* @param dtxsid DSSTox substance ID
#* @get /cpdat/qsur
#* @tag "Commercial and Industrial Usage"
#* @tag "Environmental Fate and Exposure Annotations"
function(res, req, dtxsid = "", pages=1) {
    anno_cpd_qsur <- run_query(paste0("
        SELECT DISTINCT
            bcc.preferred_name,
            bcc.casrn,
            bc.dsstox_substance_id,
            cqm.qsur_model,
            cqm.qsur_description,
            cqm.qsur_model_release,
            bcq.qsur_values,
            bcd.domain

        FROM
            base_chemicals bc,
            base_chemical_compounds bcc,
            base_chemical_qsur_predictions_2024 bcq,
            base_chemical_qsur_domain_2024 bcd,
            cpd_qsur_models cqm
        WHERE
            bc.dsstox_substance_id = $1
        AND bc.epa_id = bcc.epa_id
        AND bc.epa_id = bcq.epa_id
        AND bcq.epa_id = bcd.epa_id
        AND bcq.qsur_id = bcd.qsur_id
        AND bcq.qsur_id = cqm.qsur_id
        LIMIT $2
    "), args=list(dtxsid, as.numeric(pages)*10))
    return(anno_cpd_qsur)
}


#* Given a DSSTox substance ID, return a structural image representing the corresponding chemical in SVG format generated by RDKit.
#* @param dtxsid DSSTox substance ID
#* @get /image/structure
#* @tag "Visualization"
function(res, req, dtxsid = "") {
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
function(res, req, dtxsid = "", pages=1) {
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
            "switch__foodb_enzymes",
            "switch__foodb_flavors",
            "switch__foodb_foodcontent",
            "switch__foodb_healtheffects",
            "switch__foodb_ontology"
        )
        annotations <- chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode="")[c("anno_foodb_enzymes", "anno_foodb_flavors", "anno_foodb_foodcontent", "anno_foodb_healtheffects", "anno_foodb_ontology")]
        annotations <- lapply(annotations, function(x){
            if(nrow(x) > as.numeric(pages)*10){
                x <- x[seq_len(as.numeric(pages)*10), ]
            }
            return(x)
        })
        return(annotations)
    }
}

#* Given a DSSTox substance ID, annotate using the FooDB available in CBT. Contains information about a chemical's associated enzymes in food products.
#* @param dtxsid DSSTox substance ID
#* @get /foodb/enzymes
#* @tag "Commercial and Industrial Usage"
#* @tag "Environmental Fate and Exposure Annotations"
function(res, req, dtxsid = "", pages=1) {
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
            "switch__foodb_enzymes"
        )
        annotations <- chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode="")[[c("anno_foodb_enzymes")]]
        if(nrow(annotations) > as.numeric(pages)*10){
            annotations <- annotations[seq_len(as.numeric(pages)*10), ]
        }
        return(annotations)
    }
}

#* Given a DSSTox substance ID, annotate using the FooDB available in CBT. Contains information about a chemical's usage as a flavor in food products.
#* @param dtxsid DSSTox substance ID
#* @get /foodb/flavors
#* @tag "Commercial and Industrial Usage"
#* @tag "Environmental Fate and Exposure Annotations"
function(res, req, dtxsid = "", pages=1) {
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
            "switch__foodb_flavors"
        )
        annotations <- chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode="")[[c("anno_foodb_flavors")]]
        if(nrow(annotations) > as.numeric(pages)*10){
            annotations <- annotations[seq_len(as.numeric(pages)*10), ]
        }
        return(annotations)
    }
}

#* Given a DSSTox substance ID, annotate using the FooDB available in CBT. Contains information about a chemical's usage as an ingredient in food products.
#* @param dtxsid DSSTox substance ID
#* @get /foodb/content
#* @tag "Commercial and Industrial Usage"
#* @tag "Environmental Fate and Exposure Annotations"
function(res, req, dtxsid = "", pages=1) {
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
            "switch__foodb_foodcontent"
        )
        annotations <- chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode="")[[c("anno_foodb_foodcontent")]]
        if(nrow(annotations) > as.numeric(pages)*10){
            annotations <- annotations[seq_len(as.numeric(pages)*10), ]
        }
        return(annotations)
    }
}

#* Given a DSSTox substance ID, annotate using the FooDB available in CBT. Contains information about a chemical's health effects when used in food products.
#* @param dtxsid DSSTox substance ID
#* @get /foodb/effects
#* @tag "Commercial and Industrial Usage"
#* @tag "Environmental Fate and Exposure Annotations"
function(res, req, dtxsid = "", pages=1) {
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
            "switch__foodb_healtheffects"
        )
        annotations <- chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode="")[[c("anno_foodb_healtheffects")]]
        if(nrow(annotations) > as.numeric(pages)*10){
            annotations <- annotations[seq_len(as.numeric(pages)*10), ]
        }
        return(annotations)
    }
}


#* Given a DSSTox substance ID, annotate using the FooDB available in CBT. Contains information about a chemical's ontological relationships when used in food products.
#* @param dtxsid DSSTox substance ID
#* @get /foodb/ontology
#* @tag "Commercial and Industrial Usage"
#* @tag "Environmental Fate and Exposure Annotations"
function(res, req, dtxsid = "", pages=1) {
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
            "switch__foodb_ontology"
        )
        annotations <- chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode="")[[c("anno_foodb_ontology")]]
        if(nrow(annotations) > as.numeric(pages)*10){
            annotations <- annotations[seq_len(as.numeric(pages)*10), ]
        }
        return(annotations)
    }
}



#* Given a DSSTox substance ID, annotate using the Generally Recognized As Safe (GRAS) database available in CBT. Contains SCOGS reports for chemicals.
#* @param dtxsid DSSTox substance ID
#* @get /gras
#* @tag "Safety Studies"
function(res, req, dtxsid = "", pages=1) {
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
            "switch__gras"
        )
        annotations <- chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode="")[["anno_gras"]]
        if(nrow(annotations) > as.numeric(pages)*10){
            annotations <- annotations[seq_len(as.numeric(pages)*10), ]
        }
        return(annotations)
    }
}

#* Given a DSSTox substance ID, annotate using the Toxin and Toxin-Target Database (T3DB) available in CBT.
#* @param dtxsid DSSTox substance ID
#* @get /t3db
#* @tag "Target Annotations"
function(res, req, dtxsid = "", pages=1) {
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
            "switch__t3db"
        )
        annotations <- chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode="")[["anno_t3db"]]
        if(nrow(annotations) > as.numeric(pages)*10){
            annotations <- annotations[seq_len(as.numeric(pages)*10), ]
        }
        return(annotations)
    }
}

#* Given a DSSTox substance ID, annotate using the InVitroDB data available in CBT.
#* @param dtxsid DSSTox substance ID
#* @get /invitrodb
#* @tag "Chemical Properties"
#* @tag "In Vitro Annotations"
function(res, req, dtxsid = "", pages=1) {
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
            "switch__invitrodb"
        )
        annotations <- chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode="")[["anno_invitrodb"]]
        if(nrow(annotations) > as.numeric(pages)*10){
            annotations <- annotations[seq_len(as.numeric(pages)*10), ]
        }
        return(annotations)
    }
}

#* Given a DSSTox substance ID, find neoplasticity data for the corresponding chemical in ToxRefDB data available in CBT.
#* @param dtxsid DSSTox substance ID
#* @get /toxrefdb/neoplasticity
#* @tag "Disease Annotations"
function(res, req, dtxsid = "", pages=1) {
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
            "switch__toxrefdb_neoplastic",
            "switch__toxrefdb_nonneoplastic"
        )
        annotations <- chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode="")[c("anno_toxrefdb_neoplastic", "anno_toxrefdb_nonneoplastic")]
        
        annotations <- lapply(annotations, function(x){
            if(nrow(x) > as.numeric(pages)*10){
                x <- x[seq_len(as.numeric(pages)*10), ]
            }
            return(x)
        })
        return(annotations)
    }
}

#* Given a DSSTox substance ID, find studies and points of departure for the corresponding chemical in ToxRefDB data available in CBT.
#* @param dtxsid DSSTox substance ID
#* @get /toxrefdb/studies
#* @tag "TODO"
function(res, req, dtxsid = "", pages=1) {
    
    ep <- run_query(paste0("
        SELECT DISTINCT
            tdbp.pod_type,
            tdbp.qualifier,
            tdbp.pod_value,
            tdbp.pod_unit,
            tdbp.dose_level,
            tdbp.max_dose_level,
            tdbp.staggered_dosing,
            tdbs.study_citation,
            tdbs.study_year,
            tdbs.study_type_guideline,
            tdbs.species,
            tdbs.strain_group,
            tdbs.admin_route,
            tdbs.substance_purity,
            tdbt.sex,
            tdbt.generation,
            tdbt.dose_period,
            tdbt.dose_duration,
            tdbt.dose_duration_unit,
            tdbt.n,
            tdbte.life_stage,
            tdbte.effect_desc_free,
            tdbte.target_site,
            tdbte.direction,
            tdbte.effect_comment,
            tdbe.effect_desc,
            tdbe.cancer_related,
            tdbep.endpoint_category,
            tdbep.endpoint_type,
            tdbep.endpoint_target,
            tdbo.tested_status,
            tdbo.reported_status
        FROM
            base_chemicals bc,
            toxrefdb_chemical tdbc,
            toxrefdb_study tdbs,
            toxrefdb_tg tdbt,
            toxrefdb_tg_effect tdbte,
            toxrefdb_effect tdbe,
            toxrefdb_endpoint tdbep,
            toxrefdb_pod tdbp,
            toxrefdb_obs tdbo
        WHERE
            bc.dsstox_substance_id = $1
        AND bc.dsstox_substance_id = tdbc.dsstox_substance_id
        AND tdbc.chemical_id = tdbs.chemical_id

        AND tdbc.chemical_id = tdbp.chemical_id
        AND tdbp.study_id = tdbt.study_id

        AND tdbs.study_id = tdbt.study_id
        AND tdbt.tg_id = tdbte.tg_id
        AND tdbte.effect_id = tdbe.effect_id
        AND tdbe.endpoint_id = tdbep.endpoint_id
        
        AND tdbo.study_id = tdbt.study_id
        AND tdbo.endpoint_id = tdbe.endpoint_id
        AND tdbo.tested_status = true
        AND tdbo.reported_status = true
        LIMIT $2

    "), args=list(dtxsid, as.numeric(pages)*10))
    return(ep)
}

#* Given a DSSTox substance ID, find studies and points of departure for the corresponding chemical in ToxRefDB data available in CBT.
#* @param dtxsid DSSTox substance ID
#* @get /toxrefdb/studies/llm
#* @tag "TODO"
function(res, req, dtxsid = "", pages=1) {
    
    ep <- run_query(paste0("
        SELECT DISTINCT
            tdbs.study_citation,
            tdbe.effect_desc,
            tdbep.endpoint_target
        FROM
            base_chemicals bc,
            toxrefdb_chemical tdbc,
            toxrefdb_study tdbs,
            toxrefdb_tg tdbt,
            toxrefdb_tg_effect tdbte,
            toxrefdb_effect tdbe,
            toxrefdb_endpoint tdbep,
            toxrefdb_pod tdbp,
            toxrefdb_obs tdbo
        WHERE
            bc.dsstox_substance_id = $1
        AND bc.dsstox_substance_id = tdbc.dsstox_substance_id
        AND tdbc.chemical_id = tdbs.chemical_id

        AND tdbc.chemical_id = tdbp.chemical_id
        AND tdbp.study_id = tdbt.study_id

        AND tdbs.study_id = tdbt.study_id
        AND tdbt.tg_id = tdbte.tg_id
        AND tdbte.effect_id = tdbe.effect_id
        AND tdbe.endpoint_id = tdbep.endpoint_id
        
        AND tdbo.study_id = tdbt.study_id
        AND tdbo.endpoint_id = tdbe.endpoint_id
        AND tdbo.tested_status = true
        AND tdbo.reported_status = true
        LIMIT $2

    "), args=list(dtxsid, as.numeric(pages)*10))
    return(ep)
}


#* Given a DSSTox substance ID, annotate using the Comparative Toxicogenomics Database (CTD) available in CBT. Contains annotations related to a chemical's biological processes, cellular components, molecular functions, associated diseases and phenotypes, and affected genes.
#* @param dtxsid DSSTox substance ID
#* @get /ctd
#* @tag "Gene Annotations"
#* @tag "Disease Annotations"
#* @tag "Biological Process Annotations"
function(res, req, dtxsid = "", pages=1) {
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
            "switch__ctd_bioprocess",
            "switch__ctd_cellcomp",
            "switch__ctd_diseases",
            "switch__ctd_genes",
            "switch__ctd_molfunct",
            "switch__ctd_phenotypes"
        )
        annotations <- as.list(chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode=""))[c("anno_ctd_bioprocess", "anno_ctd_cellcomp", "anno_ctd_diseases", "anno_ctd_genes", "anno_ctd_molfunct", "anno_ctd_phenotypes")]
        
        annotations <- lapply(annotations, function(x){
            if(nrow(x) > as.numeric(pages)*10){
                x <- x[seq_len(as.numeric(pages)*10), ]
            }
            return(x)
        })
        return(annotations)
    }
}

#* Given a DSSTox substance ID, annotate using the Comparative Toxicogenomics Database (CTD) available in CBT. Contains annotations related to a chemical's biological processes.
#* @param dtxsid DSSTox substance ID
#* @get /ctd/bp
#* @tag "Biological Process Annotations"
function(res, req, dtxsid = "", pages=1) {
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
            "switch__ctd_bioprocess"
        )
        annotations <- as.list(chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode=""))[[c("anno_ctd_bioprocess")]]
        if(nrow(annotations) > as.numeric(pages)*10){
            annotations <- annotations[seq_len(as.numeric(pages)*10), ]
        }
        return(annotations)
    }
}


#* Given a DSSTox substance ID, annotate using the Comparative Toxicogenomics Database (CTD) available in CBT. Contains annotations related to a chemical's cellular components.
#* @param dtxsid DSSTox substance ID
#* @get /ctd/cc
#* @tag "Chemical Properties"
function(res, req, dtxsid = "", pages=1) {
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
            "switch__ctd_cellcomp"
        )
        annotations <- as.list(chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode=""))[[c("anno_ctd_cellcomp")]]
        if(nrow(annotations) > as.numeric(pages)*10){
            annotations <- annotations[seq_len(as.numeric(pages)*10), ]
        }
        return(annotations)
    }
}

#* Given a DSSTox substance ID, annotate using the Comparative Toxicogenomics Database (CTD) available in CBT. Contains annotations related to a chemical's associated diseases.
#* @param dtxsid DSSTox substance ID
#* @get /ctd/diseases
#* @tag "Disease Annotations"
function(res, req, dtxsid = "", pages=1) {
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
            "switch__ctd_diseases"
        )
        annotations <- as.list(chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode=""))[[c("anno_ctd_diseases")]]
        if(nrow(annotations) > as.numeric(pages)*10){
            annotations <- annotations[seq_len(as.numeric(pages)*10), ]
        }
        return(annotations)
    }
}


#* Given a DSSTox substance ID, annotate using the Comparative Toxicogenomics Database (CTD) available in CBT. Contains annotations related to a chemical's associated genes.
#* @param dtxsid DSSTox substance ID
#* @get /ctd/genes
#* @tag "Gene Annotations"
function(res, req, dtxsid = "", pages=1) {
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
            "switch__ctd_genes"
        )
        annotations <- as.list(chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode=""))[[c("anno_ctd_genes")]]
        if(nrow(annotations) > as.numeric(pages)*10){
            annotations <- annotations[seq_len(as.numeric(pages)*10), ]
        }
        return(annotations)
    }
}



#* Given a DSSTox substance ID, annotate using the Comparative Toxicogenomics Database (CTD) available in CBT. Contains annotations related to a chemical's associated genes. Reduced data version designed for use with LLM processing.
#* @param dtxsid DSSTox substance ID
#* @get /ctd/genes/llm
#* @tag "Gene Annotations"
function(res, req, dtxsid = "", pages=1) {
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
            "switch__ctd_genes"
        )
        annotations <- as.list(chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode=""))[[c("anno_ctd_genes")]]
        if(nrow(annotations) > as.numeric(pages)*10){
            annotations <- annotations[seq_len(as.numeric(pages)*10), ]
        }
        
        annotations <- annotations$gene_symbol
        return(annotations)
    }
    
    return(annotations)
}


#* Given a DSSTox substance ID, annotate using the Comparative Toxicogenomics Database (CTD) available in CBT. Contains annotations related to a chemical's molecular functions.
#* @param dtxsid DSSTox substance ID
#* @get /ctd/mf
#* @tag "Chemical Properties"
function(res, req, dtxsid = "", pages=1) {
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
            "switch__ctd_molfunct"
        )
        annotations <- as.list(chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode=""))[[c("anno_ctd_molfunct")]]
        if(nrow(annotations) > as.numeric(pages)*10){
            annotations <- annotations[seq_len(as.numeric(pages)*10), ]
        }
        return(annotations)
    }
}

#* Given a DSSTox substance ID, annotate using the Comparative Toxicogenomics Database (CTD) available in CBT. Contains annotations related to a chemical's associated phenotypes.
#* @param dtxsid DSSTox substance ID
#* @get /ctd/phenotypes
#* @tag "Gene Annotations"
function(res, req, dtxsid = "", pages=1) {
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
            "switch__ctd_phenotypes"
        )
        annotations <- as.list(chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode=""))[[c("anno_ctd_phenotypes")]]
        if(nrow(annotations) > as.numeric(pages)*10){
            annotations <- annotations[seq_len(as.numeric(pages)*10), ]
        }
        return(annotations)
    }
}

#* Given a DSSTox substance ID, annotate using the Human Metabolome Database (HMDB) available in CBT. Contains annotations related to a chemical's biospecimen locations.
#* @param dtxsid DSSTox substance ID
#* @get /hmdb/locations/biospecimen
#* @tag "Biological Process Annotations"
#* @tag "Environmental Fate and Exposure Annotations"
function(res, req, dtxsid = "", pages=1) {
    hmdb_results <- run_query_hmdb(paste0("
        SELECT DISTINCT
            eb.dsstox_substance_id,
            hbl.id_value AS location
        FROM
            epa_base eb,
            hmdb_to_dsstox hd,
            hmdb_biological_locations hbl
        WHERE
            eb.dsstox_substance_id = $1
        AND eb.epa_id = hd.epa_id
        AND hd.hmdb_id = hbl.hmdb_id
        AND hbl.id_type = 'biospecimen_locations'
        LIMIT $2
    "), args=list(dtxsid, as.numeric(pages)*10))
    return(hmdb_results)
}


#* Given a DSSTox substance ID, annotate using the Human Metabolome Database (HMDB) available in CBT. Contains annotations related to a chemical's cellular locations.
#* @param dtxsid DSSTox substance ID
#* @get /hmdb/locations/cell
#* @tag "Biological Process Annotations"
#* @tag "Environmental Fate and Exposure Annotations"
function(res, req, dtxsid = "", pages=1) {
    hmdb_results <- run_query_hmdb(paste0("
        SELECT DISTINCT
            eb.dsstox_substance_id,
            hbl.id_value AS location
        FROM
            epa_base eb,
            hmdb_to_dsstox hd,
            hmdb_biological_locations hbl
        WHERE
            eb.dsstox_substance_id = $1
        AND eb.epa_id = hd.epa_id
        AND hd.hmdb_id = hbl.hmdb_id
        AND hbl.id_type = 'cellular_locations'
        LIMIT $2
    "), args=list(dtxsid, as.numeric(pages)*10))
    return(hmdb_results)
}

#* Given a DSSTox substance ID, annotate using the Human Metabolome Database (HMDB) available in CBT. Contains annotations related to a chemical's tissue locations.
#* @param dtxsid DSSTox substance ID
#* @get /hmdb/locations/tissue
#* @tag "Biological Process Annotations"
#* @tag "Environmental Fate and Exposure Annotations"
function(res, req, dtxsid = "", pages=1) {
    hmdb_results <- run_query_hmdb(paste0("
        SELECT DISTINCT
            eb.dsstox_substance_id,
            hbl.id_value AS location
        FROM
            epa_base eb,
            hmdb_to_dsstox hd,
            hmdb_biological_locations hbl
        WHERE
            eb.dsstox_substance_id = $1
        AND eb.epa_id = hd.epa_id
        AND hd.hmdb_id = hbl.hmdb_id
        AND hbl.id_type = 'tissue_locations'
        LIMIT $2
    "), args=list(dtxsid, as.numeric(pages)*10))
    return(hmdb_results)
}

#* Given a DSSTox substance ID, annotate using the Human Metabolome Database (HMDB) available in CBT. Contains annotations related to a chemical's associated diseases.
#* @param dtxsid DSSTox substance ID
#* @get /hmdb/diseases
#* @tag "Disease Annotations"
function(res, req, dtxsid = "", pages=1) {
    hmdb_results <- run_query_hmdb(paste0("
        SELECT DISTINCT
            eb.dsstox_substance_id,
            hbd.id_value AS disease_name
        FROM
            epa_base eb,
            hmdb_to_dsstox hd,
            hmdb_diseases hbd
        WHERE
            eb.dsstox_substance_id = $1
        AND eb.epa_id = hd.epa_id
        AND hd.hmdb_id = hbd.hmdb_id
        AND hbd.id_type = 'disease_name'
        LIMIT $2
    "), args=list(dtxsid, as.numeric(pages)*10))
    return(hmdb_results)
}

#* Given a DSSTox substance ID, annotate using the Human Metabolome Database (HMDB) available in CBT. Contains annotations related to a chemical's associated genes.
#* @param dtxsid DSSTox substance ID
#* @get /hmdb/genes
#* @tag "Gene Annotations"
function(res, req, dtxsid = "", pages=1) {
    hmdb_results <- run_query_hmdb(paste0("
        SELECT DISTINCT
            eb.dsstox_substance_id,
            hbp.id_value AS gene_name
        FROM
            epa_base eb,
            hmdb_to_dsstox hd,
            hmdb_proteins hbp
        WHERE
            eb.dsstox_substance_id = $1
        AND eb.epa_id = hd.epa_id
        AND hd.hmdb_id = hbp.hmdb_id
        AND hbp.id_type = 'gene name'
        LIMIT $2
    "), args=list(dtxsid, as.numeric(pages)*10))
    return(hmdb_results)
}

#* Given a DSSTox substance ID, annotate using the Human Metabolome Database (HMDB) available in CBT. Contains annotations related to a chemical's associated proteins.
#* @param dtxsid DSSTox substance ID
#* @get /hmdb/proteins
#* @tag "Gene Annotations"
function(res, req, dtxsid = "", pages=1) {
    hmdb_results <- run_query_hmdb(paste0("
        SELECT DISTINCT
            eb.dsstox_substance_id,
            hbp.id_value AS protein_name
        FROM
            epa_base eb,
            hmdb_to_dsstox hd,
            hmdb_proteins hbp
        WHERE
            eb.dsstox_substance_id = $1
        AND eb.epa_id = hd.epa_id
        AND hd.hmdb_id = hbp.hmdb_id
        AND hbp.id_type = 'protein name'
        LIMIT $2
    "), args=list(dtxsid, as.numeric(pages)*10))
    return(hmdb_results)
}


#* Given a DSSTox substance ID, annotate using results from PASS predictive models available in CBT.
#* @param dtxsid DSSTox substance ID
#* @get /models/pass
#* @tag "Predictive Models"
function(res, req, dtxsid = "", pages=1) {
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
        LIMIT $1
    "), args=as.list(as.numeric(pages)*10))
    
    pass_pred <- gsub("\\[", "", pass_pred)
    pass_pred <- gsub("\\]", "", pass_pred)
    pass_pred <- unlist(str_split(pass_pred, ","))
    pass_pred <- as.numeric(pass_pred)
    names(pass_pred) <- pass_models$pass_name
    as.list(pass_pred)
}

#* Given a DSSTox substance ID, annotate using results from CPDat's QSUR predictive models available in CBT.
#* @param dtxsid DSSTox substance ID
#* @get /models/cpdat/qsur
#* @tag "Predictive Models"
function(res, req, dtxsid = "", pages=1) {
    qsur_pred <- run_query(paste0("
        SELECT DISTINCT
            cqd.harmonized_function,
            cqd.probability
        FROM
            cpd_qsur_data cqd
        WHERE
            cqd.dtxsid = $1
        LIMIT $2
    "), args=list(dtxsid, as.numeric(pages)*10))
    return(qsur_pred)
}


#* Given a DSSTox substance ID, annotate using results with Tox21 pathway endpoints available in CBT.
#* @param dtxsid DSSTox substance ID
#* @get /tox21
#* @tag "Pathways"
function(res, req, dtxsid = "", pages=1) {
    tox21 <- run_query(paste0("
        SELECT DISTINCT
            tbc.endpoint_name
        FROM
            tox21_to_base_chemicals tbc
        WHERE
            tbc.dsstox_substance_id = $1
        LIMIT $2
    "), args=list(dtxsid, as.numeric(pages)*10))
    return(tox21)
}

#* Given a DSSTox substance ID, annotate using results with Tox21 assay predictive models available in CBT.
#* @param dtxsid DSSTox substance ID
#* @get /models/tox21
#* @tag "Predictive Models"
function(res, req, dtxsid = "", pages=1) {
    tox21 <- run_query(paste0("
        SELECT DISTINCT
            tbc.activity_score,
            tbc.assay_model,
            tbc.ad,
            tbc.tc
        FROM
            tox21_base_chemical_assay_predictions tbc
        WHERE
            tbc.dsstox_substance_id = $1
        LIMIT $2
    "), args=list(dtxsid, as.numeric(pages)*10))
    return(tox21)
}

#* Given a DSSTox substance ID, annotate using the DrugBank database available in CBT. Contains information regarding a chemical's usage in drugs and their targets.
#* @param dtxsid DSSTox substance ID
#* @get /drugbank
#* @tag "Drug Annotations"
#* @tag "Commercial and Industrial Usage"
function(res, req, dtxsid = "", pages=1) {
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
            "switch__drugbank_atccodes",
            "switch__drugbank_carriers",
            "switch__drugbank_enzymes",
            "switch__drugbank_targets",
            "switch__drugbank_transporters"
        )
        annotations <- chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode="")[c("anno_drugbank_atccodes", "anno_drugbank_carriers", "anno_drugbank_enzymes", "anno_drugbank_targets", "anno_drugbank_transporters")]
        
        annotations <- lapply(annotations, function(x){
            if(nrow(x) > as.numeric(pages)*10){
                x <- x[seq_len(as.numeric(pages)*10), ]
            }
            return(x)
        })
        return(annotations)
        
    }
}


#* Given a DSSTox substance ID, annotate using the DrugBank database available in CBT. Contains a chemical's ATC codes.
#* @param dtxsid DSSTox substance ID
#* @get /drugbank/atc_codes
#* @tag "Drug Annotations"
function(res, req, dtxsid = "", pages=1) {
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
            "switch__drugbank_atccodes"
        )
        annotations <- chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode="")[[c("anno_drugbank_atccodes")]]
        if(nrow(annotations) > as.numeric(pages)*10){
            annotations <- annotations[seq_len(as.numeric(pages)*10), ]
        }
        return(annotations)
    }
}

#* Given a DSSTox substance ID, annotate using the DrugBank database available in CBT. Contains information about a chemical's carriers.
#* @param dtxsid DSSTox substance ID
#* @get /drugbank/carriers
#* @tag "Drug Annotations"
#* @tag "Biological Process Annotations"
function(res, req, dtxsid = "", pages=1) {
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
            "switch__drugbank_carriers"
        )
        annotations <- chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode="")[[c("anno_drugbank_carriers")]]
        if(nrow(annotations) > as.numeric(pages)*10){
            annotations <- annotations[seq_len(as.numeric(pages)*10), ]
        }
        return(annotations)
    }
}

#* Given a DSSTox substance ID, annotate using the DrugBank database available in CBT. Contains information about a chemical's associated enzymes.
#* @param dtxsid DSSTox substance ID
#* @get /drugbank/enzymes
#* @tag "Drug Annotations"
#* @tag "Biological Process Annotations"
function(res, req, dtxsid = "", pages=1) {
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
            "switch__drugbank_enzymes"
        )
        annotations <- chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode="")[[c("anno_drugbank_enzymes")]]
        if(nrow(annotations) > as.numeric(pages)*10){
            annotations <- annotations[seq_len(as.numeric(pages)*10), ]
        }
        return(annotations)
    }
}

#* Given a DSSTox substance ID, annotate using the DrugBank database available in CBT. Contains information about a chemical's associated targets.
#* @param dtxsid DSSTox substance ID
#* @get /drugbank/targets
#* @tag "Drug Annotations"
#* @tag "Biological Process Annotations"
function(res, req, dtxsid = "", pages=1) {
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
            "switch__drugbank_targets"
        )
        annotations <- chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode="")[[c("anno_drugbank_targets")]]
        if(nrow(annotations) > as.numeric(pages)*10){
            annotations <- annotations[seq_len(as.numeric(pages)*10), ]
        }
        return(annotations)
    }
}

#* Given a DSSTox substance ID, annotate using the DrugBank database available in CBT. Contains information about a chemical's associated transporters.
#* @param dtxsid DSSTox substance ID
#* @get /drugbank/transporters
#* @tag "Drug Annotations"
#* @tag "Biological Process Annotations"
function(res, req, dtxsid = "", pages=1) {
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
            "switch__drugbank_transporters"
        )
        annotations <- chemical_annotations(selected=as.integer(internal_id$epa_id), input=input_all_annotations, selected_node=1, dtxsids=internal_id$dsstox_substance_id, switch_mode="")[[c("anno_drugbank_transporters")]]
        if(nrow(annotations) > as.numeric(pages)*10){
            annotations <- annotations[seq_len(as.numeric(pages)*10), ]
        }
        return(annotations)
    }
}

#* Given a SMILES string, check for OChem structural alerts and notable substructures.
#* @param smiles SMILES string
#* @get /alerts/ochem
#* @tag "Substructure Search"
function(res, req, smiles = "", pages=1) {
    # OChem
    alerts_ochem <- run_query(paste0("
        SELECT DISTINCT
            osa.smarts AS alert,
            'ochem' AS source,
            ost.response AS description
        FROM
            ochem_structural_alerts_mols osa,
            ochem_structural_alerts_toxpipe_annotated ost
        WHERE
            mol_from_smiles($1)@> osa.m
        AND osa.ochem_id = ost.ochem_id
        LIMIT $2
    "), args=list(smiles, as.numeric(pages)*10))
    return(alerts_ochem)
}

#* Given a SMILES string, check for ChEMBL structural alerts and notable substructures.
#* @param smiles SMILES string
#* @get /alerts/chembl
#* @tag "Substructure Search"
function(res, req, smiles = "", pages=1) {
    # ChEMBL
    alerts_chembl <- run_query(paste0("
        SELECT DISTINCT
            csa.alert_name AS alert,
            'chembl' AS source,
            cst.response AS description
        FROM
            chembl_structural_alerts csa,
            chembl_structural_alerts_toxpipe_annotated cst
        WHERE
            mol_from_smiles($1)@> mol_from_smarts(csa.smarts)
        AND csa.chemblsa_id = cst.chemblsa_id
        LIMIT $2
    "), args=list(smiles, as.numeric(pages)*10))
    return(alerts_chembl)
}

#* Given a SMILES string, check for Saagar structural alerts and notable substructures.
#* @param smiles SMILES string
#* @get /alerts/saagar
#* @tag "Substructure Search"
function(res, req, smiles = "", pages=1) {
    # Saagar
    alerts_saagar <- run_query(paste0("
        SELECT DISTINCT
            ssa.alert_name AS alert,
            'saagar' AS source,
            sst.response AS description
        FROM
            saagar_structural_alerts_mols ssa,
            saagar_structural_alerts_toxpipe_annotated sst
        WHERE
            mol_from_smiles($1)@> ssa.m
        AND ssa.saagarsa_id = sst.saagarsa_id
        LIMIT $2
    "), args=list(smiles, as.numeric(pages)*10))
    return(alerts_saagar)
}

#* Given a DssTox Substance ID, check for ToxPrint structural alerts and notable substructures.
#* @param dtxsid DSSTox Substance ID
#* @get /alerts/toxprint
#* @tag "Substructure Search"
function(res, req, dtxsid = "", pages=1) {
    
    # ToxPrint
    alerts_toxprint <- run_query(paste0("
        SELECT DISTINCT
            ta.toxprint_alert,
            'toxprint' AS source
        FROM
            base_chemical_toxprint_positives_2024 mv,
            toxprint_alerts ta,
            base_chemical_to_smiles bcs,
            base_chemicals bc
        WHERE
            bc.epa_id = mv.epa_id
        AND bc.dsstox_substance_id = $1
        AND mv.toxpr_id = ta.toxpr_id
        AND mv.toxprint_values = 1
        LIMIT $2
    "), args=list(dtxsid, as.numeric(pages)*10))
    
    return(alerts_toxprint)
}

#* Given a DSSTox substance ID, annotate with chemical properties as detailed in PubChem that are available in CBT.
#* @param dtxsid DSSTox substance ID
#* @get /pubchem/properties
#* @tag "Chemical Properties"
#* @tag "PubChem"
function(res, req, dtxsid = "", pages=1) {
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
        LIMIT $2
    "), args=list(dtxsid, as.numeric(pages)*10))
    return(props)
}

#* Given a DSSTox substance ID, annotate with synonyms from PubChem.
#* @param dtxsid DSSTox substance ID
#* @get /pubchem/properties/synonyms
#* @tag "PubChem"
function(res, req, dtxsid = "", pages=1) {
    props <- run_query(paste0("
        SELECT DISTINCT
            bpcs.synonym
        FROM
            base_chemicals bc,
            base_chemical_to_pubchem_synonyms bpcs
        WHERE
            bc.dsstox_substance_id = $1
        AND bc.epa_id = bpcs.epa_id
        LIMIT $2
    "), args=list(dtxsid, as.numeric(pages)*10))
    
    if(nrow(props) > 0){
        props <- props$synonym
    }
    return(props)
}

#* Given a DSSTox substance ID, annotate with its mass from PubChem.
#* @param dtxsid DSSTox substance ID
#* @get /pubchem/properties/mass
#* @tag "PubChem"
function(res, req, dtxsid = "") {
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
    
    if(nrow(props) > 0){
        props <- props[props$pubchem_attribute == "exactmass", "pubchem_value"]
        props <- paste0(props, " Da")
    }
    return(props)
}

#* Given a DSSTox substance ID, annotate with its molecular formula from PubChem.
#* @param dtxsid DSSTox substance ID
#* @get /pubchem/properties/formula
#* @tag "PubChem"
function(res, req, dtxsid = "") {
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
    
    if(nrow(props) > 0){
        props <- props[props$pubchem_attribute == "mf", "pubchem_value"]
    }
    return(props)
}

#* Given a DSSTox substance ID, annotate with its molecular weight from PubChem.
#* @param dtxsid DSSTox substance ID
#* @get /pubchem/properties/weight
#* @tag "PubChem"
function(res, req, dtxsid = "") {
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
    
    if(nrow(props) > 0){
        props <- props[props$pubchem_attribute == "mw", "pubchem_value"]
        props <- paste0(props, " g/mol")
    }
    return(props)
}

#* Given a DSSTox substance ID, annotate with its XLogP from PubChem.
#* @param dtxsid DSSTox substance ID
#* @get /pubchem/properties/xlogp
#* @tag "PubChem"
function(res, req, dtxsid = "") {
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
    
    if(nrow(props) > 0){
        props <- props[props$pubchem_attribute == "xlogp", "pubchem_value"]
    }
    return(props)
}



#* Given a DSSTox substance ID, annotate with chemical properties as detailed in the EPA CompTox Chemicals Dashboard that are available in CBT.
#* @param dtxsid DSSTox substance ID
#* @get /epa/properties
#* @tag "Chemical Properties"
function(res, req, dtxsid = "", pages=1) {
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
        LIMIT $2
    "), args=list(dtxsid, as.numeric(pages)*10))
    return(props)
}


#* Given a DSSTox substance ID, return bioassays in PubChem that are available in CBT.
#* @param dtxsid DSSTox substance ID
#* @get /pubchem/bioassays
#* @tag "Biological Process Annotations"
#* @tag "Assays"
function(res, req, dtxsid = "", pages=1) {
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
        LIMIT $2
    "), args=list(dtxsid, as.numeric(pages)*10))
    return(assays)
}

#* Given a DSSTox substance ID, return test articles from the NTP that are available in CBT.
#* @param dtxsid DSSTox substance ID
#* @get /ntp/test_articles
#* @tag "Assays"
function(res, req, dtxsid = "", pages=1) {
    assays <- run_query(paste0("
        SELECT DISTINCT
            bcta.testarticle,
            bcta.study_area,
            bcta.study_number,
            nta.web_url
        FROM
            base_chemicals bc,
            base_chemical_test_articles bcta,
            ntp_test_article_urls nta
        WHERE
            bc.dsstox_substance_id = $1
        AND bc.epa_id = bcta.epa_id
        AND bcta.ntpid = nta.ntpid
        LIMIT $2
    "), args=list(dtxsid, as.numeric(pages)*10))
    return(assays)
}


#* Given a DSSTox substance ID, annotate with a chemical's upper 95th percentile SEEM3 exposure available in CBT.
#* @param dtxsid DSSTox substance ID
#* @get /seem3
#* @tag "Environmental Fate and Exposure Annotations"
function(res, req, dtxsid = "", pages=1) {
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
            LIMIT $2
        "), args=list(as.character(internal_id$epa_id), as.numeric(pages)*10))
        
        
        exp <- rbindlist(list(exp_seem3))
        return(exp)
    }
}

#* Given a DSSTox substance ID, annotate with a chemical's Cmax concentration data.
#* @param dtxsid DSSTox substance ID
#* @get /cmax
#* @tag "Environmental Fate and Exposure Annotations"
function(res, req, dtxsid = "", pages=1) {
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
        LIMIT $2
    "), args=list(dtxsid, as.numeric(pages)*10))
    return(cmax)
}

#* Given a DSSTox substance ID, annotate with a chemical's Cmax concentration data.
#* @param dtxsid DSSTox substance ID
#* @get /superfund
#* @tag "Environmental Fate and Exposure Annotations"
function(res, req, dtxsid = "", pages=1) {
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
        LIMIT $2
    "), args=list(dtxsid, as.numeric(pages)*10))
    return(sfund)
}

query_mesh <- function(q="", n=10, exact=FALSE){
    res <- tryCatch({
        match <- "contains"
        if(exact == TRUE) match="exact" 
        
        body_content <- list(
            label=q,
            match=match,
            year="current",
            limit=n
        )
        res <- httr::POST(url="https://id.nlm.nih.gov/mesh/lookup/descriptor", body=body_content, encode="json")
        json <- fromJSON(httr::content(res, "text"))
        return(json)
        
    }, error=function(cond){
        print(cond)
        return(NULL)
    })
}



#* Given a disease name, check for associated chemicals
#* @param name disease name
#* @param n num of disease name matches to return
#* @param exact match exact name?
#* @param n_chem number of chemicals per disease name to return
#* @get /ctd/diseases/chemicals
#* @tag "Disease Annotations"
function(res, req, name = "", n=10, exact=FALSE, n_chem=10) {
    
    # Get MeSH IDs for input disease
    mesh_ids <- query_mesh(q=name, n=n, exact=as.booltype(exact))
    mesh_ids$mesh_id <- unlist(lapply(mesh_ids$resource, function(x){
        tmp <- unlist(str_split(x, "/"))
        return(
            paste0("MESH:", tmp[length(tmp)])
        )
    }))
    
    if(length(mesh_ids) < 1){
        return(data.frame())
    }
    
    disease2chemical <- run_query(paste0("
    SELECT DISTINCT
        bcc.preferred_name,
        ccd.disease_name,
        ccd.direct_evidence,
        ccd.inference_score
    FROM
        ctd_chemicals_to_diseases ccd,
        ctd_to_base_chemicals cbc,
        base_chemical_compounds bcc
    WHERE 
        ccd.disease_id IN (", paste0(lapply(seq_len(length(mesh_ids$mesh_id)), function(x) paste0("$", x)), collapse=", "), ")
    AND ccd.ctd_id = cbc.ctd_id
    AND cbc.epa_id = bcc.epa_id
    ORDER BY ccd.direct_evidence ASC, ccd.inference_score DESC
    "), args=as.list(mesh_ids$mesh_id))
    
    disease2chemical <- split(disease2chemical, disease2chemical$disease_name)
    
    disease2chemical <- lapply(disease2chemical, function(x) {
        x[seq_len(n_chem), ]
    })
    
    return(disease2chemical)
}


#* Given a gene name, check for associated chemicals
#* @param name gene name
#* @param n num of chemicals to return
#* @param exact match exact name?
#* @get /ctd/genes/chemicals
#* @tag "Gene Annotations"
function(res, req, name = "", n=20, exact=FALSE) {
    gene2chemical <- run_query(paste0("
    SELECT DISTINCT
        bcc.preferred_name,
        ccg.gene_symbol,
        ccg.interaction
    FROM
        ctd_chemicals_to_genes ccg,
        ctd_to_base_chemicals cbc,
        base_chemical_compounds bcc
    WHERE 
    ", if(as.booltype(exact) == FALSE) paste0("UPPER(ccg.gene_symbol) LIKE UPPER(CONCAT('%', $1::text, '%'))") else paste0("UPPER(ccg.gene_symbol) = UPPER($1::text)"), "
    AND ccg.ctd_id = cbc.ctd_id
    AND cbc.epa_id = bcc.epa_id
    LIMIT $2
    "), args=list(name, n))
    
    gene2chemical <- split(gene2chemical, gene2chemical$gene_symbol)
    
    gene2chemical <- lapply(gene2chemical, function(x){
        split(x, x$interaction)
    })
    
    return(gene2chemical)
}



#* Given a disease name, check for associated chemicals
#* @param name disease name
#* @param n num of chemicals to return
#* @param exact match exact name?
#* @get /hmdb/diseases/chemicals
#* @tag "Disease Annotations"
function(res, req, name = "", n=20, exact=FALSE) {
    disease2chemical <- run_query(paste0("
    SELECT DISTINCT
        bcc.preferred_name,
        hd.annotation AS disease_name
    FROM
        base_chemical_hmdb_diseases_saagar_clusters hd,
        base_chemical_compounds bcc
    WHERE 
    ", if(as.booltype(exact) == FALSE) paste0("UPPER(hd.annotation) LIKE UPPER(CONCAT('%', $1::text, '%'))") else paste0("UPPER(hd.annotation) = UPPER($1::text)"), "
    AND hd.epa_id = bcc.epa_id
    LIMIT $2
    "), args=list(name, n))
    
    disease2chemical <- split(disease2chemical, disease2chemical$disease_name)
    
    return(disease2chemical)
}


#* Given a gene name, check for associated chemicals
#* @param name gene name
#* @param n num of chemicals to return
#* @param exact match exact name?
#* @get /hmdb/genes/chemicals
#* @tag "Gene Annotations"
function(res, req, name = "", n=20, exact=FALSE) {
    gene2chemical <- run_query(paste0("
    SELECT DISTINCT
        bcc.preferred_name,
        hg.annotation AS gene_symbol
    FROM
        base_chemical_hmdb_genes_saagar_clusters hg,
        base_chemical_compounds bcc
    WHERE 
    ", if(as.booltype(exact) == FALSE) paste0("UPPER(hg.annotation) LIKE UPPER(CONCAT('%', $1::text, '%'))") else paste0("UPPER(hg.annotation) = UPPER($1::text)"), "
    AND hg.epa_id = bcc.epa_id
    LIMIT $2
    "), args=list(name, n))
    
    gene2chemical$interaction <- apply(gene2chemical, 1, function(x){
        return(paste0("affects expression of ", x["gene_symbol"]))
    })
    
    gene2chemical <- split(gene2chemical, gene2chemical$gene_symbol)
    
    gene2chemical <- lapply(gene2chemical, function(x){
        split(x, x$interaction)
    })
    
    return(gene2chemical)
}

#* Given a gene name, check for associated chemicals
#* @param name gene name
#* @param n num of chemicals to return
#* @param exact match exact name?
#* @get /drugbank/genes/chemicals
#* @tag "Gene Annotations"
function(res, req, name = "", n=20, exact=FALSE) {
    gene2chemical <- run_query(paste0("
    SELECT DISTINCT
        bcc.preferred_name,
        dcg.genename AS gene_symbol,
        dcg.specific_function AS interaction
    FROM
        drugbank_chemicals_to_genes dcg,
        drugbank_to_base_chemicals dbc,
        drugbank_curated_chemicals dcc,
        base_chemical_compounds bcc
    WHERE 
    ", if(as.booltype(exact) == FALSE) paste0("UPPER(dcg.genename) LIKE UPPER(CONCAT('%', $1::text, '%'))") else paste0("UPPER(dcg.genename) = UPPER($1::text)"), "
    AND dbc.drugbank_id = dcc.drugbank_id
    AND dcc.db_id = dcg.db_id
    AND dbc.epa_id = bcc.epa_id
    LIMIT $2
    "), args=list(name, n))
    
    gene2chemical <- split(gene2chemical, gene2chemical$gene_symbol)
    
    gene2chemical <- lapply(gene2chemical, function(x){
        split(x, x$interaction)
    })
    
    return(gene2chemical)
}


#* Given a DTXSID, check for associated genes
#* @param dtxsid Chemical DSSTox Substance ID
#* @get /drugbank/genes
#* @tag "Gene Annotations"
function(res, req, dtxsid = "", pages=1) {
    chemical2gene <- run_query(paste0("
    SELECT DISTINCT
        bc.dsstox_substance_id,
        bcc.preferred_name,
        dcg.genename AS gene_symbol,
        dcg.specific_function AS interaction
    FROM
        drugbank_chemicals_to_genes dcg,
        drugbank_to_base_chemicals dbc,
        drugbank_curated_chemicals dcc,
        base_chemical_compounds bcc,
        base_chemicals bc
    WHERE 
        bc.dsstox_substance_id = $1
    AND bc.epa_id = bcc.epa_id
    AND dbc.drugbank_id = dcc.drugbank_id
    AND dcc.db_id = dcg.db_id
    AND dbc.epa_id = bcc.epa_id
    LIMIT $2
    "), args=list(dtxsid, as.numeric(pages)*10))

    
    return(chemical2gene)
}