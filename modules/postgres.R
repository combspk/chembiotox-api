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
library(sortable)
library(stringi)
library(stringr)
library(svglite)
library(RPostgres)

# Load variables from config.yml & assign to global variables
config_vars <- config::get("db") # Do not attach the config pkg directly, as per: https://github.com/rstudio/config
DB_HOST <- config_vars$host
DB_PORT <- config_vars$port
DB_USER <- config_vars$user
DB_PASS <- config_vars$pass

# Define database driver
pgdrv <- RPostgres::Postgres()

# Define reusable functions

# Redefine %like% operator to be case insensitive (from https://stackoverflow.com/questions/41425699/how-to-get-the-like-operator-to-be-case-insensitive)
`%like%` <- function (x, pattern) { 
    stringi::stri_detect_regex(x, pattern, case_insensitive=TRUE)
}

pool_main <- dbPool(
    drv=pgdrv,
    dbname="chembiotox_v2",
    host=DB_HOST,
    port=DB_PORT,
    user=DB_USER,
    password=DB_PASS
)

pool_vec <- dbPool(
    drv=pgdrv,
    dbname="chembiotox_vectors_v1",
    host=DB_HOST,
    port=DB_PORT,
    user=DB_USER,
    password=DB_PASS
)

pool_chembl <- dbPool(
    drv=pgdrv,
    dbname="chembl_34",
    host=DB_HOST,
    port=DB_PORT,
    user=DB_USER,
    password=DB_PASS
)

run_query <- function(query="", args=list()){
    df <- data.frame()
    try_query <- tryCatch({ # Attempt to query db
        # Checkout connection from db pool
        conn <- localCheckout(pool_main)
        
        res <- dbSendQuery(conn, query)
        if(length(args) > 0){
            dbBind(res, args)
        }
        df <- dbFetch(res)
        dbClearResult(res)
    }, error=function(cond){
        print("Error in query")
        print(cond)
    })
    
    resp <- data.frame()
    if(nrow(df) > 0){
        resp <- df
    }
    return(resp)
}

run_query_vec <- function(query="", args=list()){
    df <- data.frame()
    try_query <- tryCatch({ # Attempt to query db
        # Checkout connection from db pool
        conn <- localCheckout(pool_vec)
        
        res <- dbSendQuery(conn, query)
        if(length(args) > 0){
            dbBind(res, args)
        }
        df <- dbFetch(res)
        dbClearResult(res)
    }, error=function(cond){
        print("Error in query")
        print(cond)
    })
    
    resp <- data.frame()
    if(nrow(df) > 0){
        resp <- df
    }
    return(resp)
}

run_query_chembl <- function(query="", args=list()){
    df <- data.frame()
    try_query <- tryCatch({ # Attempt to query db
        # Checkout connection from db pool
        conn <- localCheckout(pool_chembl)
        
        res <- dbSendQuery(conn, query)
        if(length(args) > 0){
            dbBind(res, args)
        }
        df <- dbFetch(res)
        dbClearResult(res)
    }, error=function(cond){
        print("Error in query")
        print(cond)
    })
    
    resp <- data.frame()
    if(nrow(df) > 0){
        resp <- df
    }
    return(resp)
}

pool_metabolites <- dbPool(
    drv=pgdrv,
    #dbname="2023_leadscope_results",
    dbname="chembiotox_v2",
    host=DB_HOST,
    port=DB_PORT,
    user=DB_USER,
    password=DB_PASS
)

run_query_metabolites <- function(query, args=list()){
    df <- data.frame()
    try_query <- tryCatch({ # Attempt to query db
        # Checkout connection from db pool
        conn <- localCheckout(pool_metabolites)
        
        res <- dbSendQuery(conn, query)
        if(length(args) > 0){
            dbBind(res, args)
        }
        df <- dbFetch(res)
        dbClearResult(res)
    }, error=function(cond){
        print("Error in query")
        print(cond)
    })
    
    resp <- data.frame()
    if(nrow(df) > 0){
        resp <- df
    }
    return(resp)
}