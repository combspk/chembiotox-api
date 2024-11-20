metabolites_query <- function(epa_id, enzyme="both", level=3){
    
    if(level < 1 | level > 3){
        return("Error: level must be 1, 2, or 3.")
    }
    
    if(!enzyme %in% c("cyp", "ugt", "both")){
        return("Error: unknown enzyme specified.")
    }
    
    metabolite_identifier <- "metalev"
    if(enzyme == "ugt") {
        metabolite_identifier <- "ugtlevel"
    }

    res <- run_query(paste0("
        SELECT
            '", enzyme, "' AS enzyme,
            ct.", metabolite_identifier, level, "_id AS metabolite_id,
            ", level, " AS level,
            mlcs.epa_id,
            
            -- Use different previous levels as the parent depending on the current level of metabolism
            ", if(level == 1){paste0("mlcs.base_smi AS parent_smi_id, 0 AS parent_level")}, "
            ", if(level == 2){paste0("mlcs.level1_smi AS parent_smi_id, mlcs.level1_id AS parent_level")}, "
            ", if(level == 3){paste0("mlcs.level2_smi AS parent_smi_id, mlcs.level2_id AS parent_level")}, "
            ,
            mlcs.", if(enzyme=="cyp"){paste0("level", level, "_smi")} else {paste0("ugt", level, "_smi")}, " AS child_smi_id,
            mlcs.", if(enzyme=="cyp"){paste0("level", level, "_id")} else {paste0("ugt_id")}, " AS child_level,
            ss.smiles,
            ct.enzymes,
            ct.mass_list,
        
            -- Only include enzyme info for CYP*-driven metabolism
            ", if(enzyme=="cyp"){ paste0("
            ct.pctyield,
            ct.\"CYP2C19\",
            ct.\"CYP2D6\",
            ct.\"CYP2C9\",
            ct.\"CYP3A4\",
            ct.\"CYP1A2\",
            ct.sum,
            ")}, "
        
            mol_to_svg(mol_from_smiles(ss.smiles::cstring))::text AS svg
        FROM
            crosstab(format(
                'SELECT
                    mlcs.", metabolite_identifier, level, "_id,
                    m1ci.", if(enzyme=="ugt") paste0("ugt"), "level", level, "_attribute AS level", level, "_attribute,
                    m1ci.", if(enzyme=="ugt") paste0("ugt"), "level", level, "_value AS level", level, "_value
                FROM
                    metabolite_level", level, "_cyp_to_", if(enzyme=="ugt") paste0("ugt_to_"), "smiles mlcs,
                    metabolite_level", level, "_cyp_", if(enzyme=="ugt") paste0("to_ugt_"), "information m1ci
                WHERE
                    m1ci.", metabolite_identifier, level, "_id = mlcs.", metabolite_identifier, level, "_id
                AND mlcs.epa_id = %s
                ORDER BY 1', $1::int),
                $$VALUES(\'enzymes\'::text), (\'mass_list\'::text)", if(enzyme=="cyp") paste0(", (\'pctyield\'::text), (\'CYP2C19\'::text), (\'CYP2D6\'::text), (\'CYP2C9\'::text), (\'CYP3A4\'::text), (\'CYP1A2\'::text), (\'sum\'::text)"), "$$
            ) AS ct (\"", metabolite_identifier, level, "_id\" int, \"enzymes\" text, \"mass_list\" text", if(enzyme=="cyp") paste0(", \"pctyield\" text, \"CYP2C19\" text, \"CYP2D6\" text, \"CYP2C9\" text, \"CYP3A4\" text, \"CYP1A2\" text, \"sum\" text"), "),
            metabolite_level", level, "_cyp_to_", if(enzyme=="ugt") paste0("ugt_to_"), "smiles mlcs,
            smiles_strings ss
        WHERE
            mlcs.", metabolite_identifier, level, "_id = ct.", metabolite_identifier, level, "_id
        AND mlcs.level", level, "_smi = ss.smi_id
    "), args=list(epa_id))
    return(res)
}


metabolites_network <- function(df){
    df$parent_smi_id <- as.double(df$parent_smi_id)
    df$child_smi_id <- as.double(df$child_smi_id)
    # Remove anything without a percent yield (these metabolites are very unlikely to be created in real life)
    # Just do this for CYP* for now, as we don't track %yield for UGT*
    # df <- df[!(is.na(df$pctyield) & df$enzyme == "cyp"), ]
    
    df_removed <- df[(is.na(df$pctyield) & df$enzyme == "cyp"), "child_smi_id"]
    
    # Remove anything without a parent node
    # df <- df[df$parent_smi_id %in% df$child_smi_id | df$parent_smi_id == as.integer64(0), ]
    
    # 
    # #TODO - fix this
    # # df <- do.call(rbind, apply(df, 1, function(x){
    # #     if(x["parent_smi_id"] %in% df$child_smi_id == TRUE){
    # #         return(x)
    # #     }
    # #     return(NULL)
    # # }))
    # 
    # # print(df[as.double(df$child_smi_id) == 1510447765, ])
    # # 
    # # print(df[as.double(df$child_smi_id) == 574726927, ])
    # # 
    # # print(as.integer64(574726927) %in% df$child_smi_id)
    # 
    # print(df$parent_smi_id %in% df$child_smi_id)
    
    nodes <- unique(data.frame(
        id=df$child_smi_id,
        parent_id=df$parent_smi_id,
        child_id=df$child_smi_id,
        label=paste0(df$child_smi_id),
        shape="image",
        color=unlist(lapply(df$enzyme, function(enz){
            if(!is.na(enz)){
                if(enz == "cyp"){
                    return("blue")
                } else {
                    return("orange")
                }
            }
            return("black")
        })),
        image=paste0("data:image/svg+xml;charset=utf-8,", URLencode(df$svg, reserved=TRUE)),
        size=50,
        enzyme=df$enzyme,
        group=df$enzyme,
        level=df$level,
        title=paste0("
            <b>Metabolite ID</b>: ", df$child_smi_id, "<br>
            <b>Enzyme</b>: ", df$enzyme, "<br>
            <b>Level</b>: ", df$level, "<br>
            <b>Parent ID</b>: ", df$parent_smi_id, "
        "),
        stringsAsFactors=FALSE
    ))
    
    # TODO: handle this better - this is just getting the first instance if we have multiple nodes with the same SMI_ID - parts of the same chemical entry
    nodes <- nodes[!duplicated(nodes[, c("id")]), ]
    
    edges <- data.frame(
        from=df$parent_smi_id,
        to=df$child_smi_id,
        label=df$pctyield,
        color=unlist(lapply(df$enzyme, function(enz){
            if(!is.na(enz)){
                if(enz == "cyp"){
                    return("blue")
                } else {
                    return("orange")
                }
            }
            return("black")
        })),
        
        title=paste0("
            <b>Child Metabolite ID</b>: ", df$child_smi_id, "<br>
            <b>Parent Metabolite ID</b>: ", df$parent_smi_id, "<br>
            <b>Metabolism Level</b>: ", df$enzyme, "; level ", df$level, "<br>
            <b>Mass List</b>: ", df$mass_list, "<br>
            <b>% yield</b>: ", df$pctyield, "<br>
            <b>Enzymes</b>: ", df$enzymes, "<br>
            <b>- CYP2C19</b>: ", df$CYP2C19, "<br>
            <b>- CYP2D6</b>: ", df$CYP2D6, "<br>
            <b>- CYP2C9</b>: ", df$CYP2C9, "<br>
            <b>- CYP3A4</b>: ", df$CYP3A4, "<br>
            <b>- CYP1A2</b>: ", df$CYP1A2, "<br>
            <b>Sum</b>: ", df$sum, "<br>
        "),
        
        id=paste0(df$parent_smi_id, "__", df$child_smi_id)
    )
    
    edges <- edges[!duplicated(edges[, c("from", "to")]), ]
    nodes <- nodes[!(nodes$parent_id %in% df_removed), ]
    nodes <- nodes[!(nodes$child_id %in% df_removed), ]
    nodes <- nodes[(nodes$parent_id %in% nodes$child_id) | nodes$parent_id == 0, ]
    
    
    #print(df[df$parent_smi_id == 574726945 | df$child_smi_id == 574726945,])
    
     
    metabolites_graph <- visNetwork(nodes=nodes, edges=edges, height="1000px", width="100%", main="Select one or more metabolites") %>%
        visNodes(borderWidth=2, fixed=FALSE, shapeProperties=list(useBorderWithImage=TRUE)) %>%
        visEdges(physics=FALSE, smooth=FALSE, width=2, arrows=c("to")) %>%
        visHierarchicalLayout(levelSeparation=500) %>%
        visIgraphLayout() %>%
        visOptions(collapse=list(
            enabled=TRUE,
            keepCoord=TRUE
        ), nodesIdSelection=TRUE, selectedBy="enzyme", highlightNearest=list(enabled=TRUE, algorithm="hierarchical")) %>%
        visInteraction(multiselect=TRUE, navigationButtons=TRUE) %>%
        visPhysics(stabilization=FALSE, enabled=FALSE) %>%
        visClusteringByGroup(groups=c("cyp", "ugt")) %>%
        visEvents(
            selectEdge='function(properties) {Shiny.setInputValue("selectEdge_metabolite", properties, {priority:"event"});}',
            selectNode='function(properties) {Shiny.setInputValue("selectNode_metabolite", properties, {priority:"event"});}'
        )
    
    return(metabolites_graph)
}


load_metabolite_to_workspace <- function(smi_ids, metabolites){
    
    # TODO: handle this better - this is just getting the first instance if we have multiple nodes with the same SMI_ID - parts of the same chemical entry
    metabolites <- metabolites[!duplicated(metabolites[, c("child_smi_id")]), ]
    
    # Load metabolite leadscope model results
    model_results <- get_metabolite_leadscope_models(smi_ids)
    metabolite_cards <- apply(metabolites, 1, function(x){
        # Make metabolite data card
        div(
            box(
                fluidRow(
                    column(7,
                        HTML(
                            gsub("width='500", "width='250",
                                gsub("height='500", "height='200",
                                    x["svg"]
                                )
                            )
                        )
                    ),
                    column(5,
                        HTML(paste0("<p><b>ID</b>: ", x["child_smi_id"], "</p>")),
                        HTML(paste0("<p><b>SMILES</b>: ", x["smiles"], "</p>")),
                        HTML(paste0("<p><b>Enzyme</b>: ", x["enzyme"], "</p>"))
                    )
                ),
                title=x["child_smi_id"],
                footer=actionButton(inputId=paste0("actionButton__remove_", x["child_smi_id"]), label="Remove chemical"),
                width=4,
                solidHeader=TRUE,
                collapsible=TRUE,
                collapsed=FALSE
            )
        )
    })
    names(metabolite_cards) <- as.double(metabolites$child_smi_id)
    return(list(cards=metabolite_cards, models=model_results))
}

generate_metabolite_leadscope_graph <- function(models){
    p <- plot_ly(
        models[[1]],
        x=~model_name,
        y=~interpretation,
        name=~smi_id[1],
        type="bar"
    )
    if(length(models) > 1){
        for(i in seq(2, length(models))){
            p <- p %>% add_trace(
                models[[i]],
                x=models[[i]]$model_name,
                y=models[[i]]$interpretation,
                name=models[[i]]$smi_id[1],
                type="bar"
            )
        }
    }
    return(p)
}