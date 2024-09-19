generate_plotly <- function(df, title="", mode="annotations"){
    if(mode == "blank"){
        return(plotly_empty(type="bar") %>% config(displayModeBar=FALSE) %>% layout(title=list(text="No results found for the given cluster.", yref="paper", y=0.5)))
    }
    p <- plot_ly(
        df,
        y=df$annotation,
        x=df$percentage,
        text=df$num_of_chemicals,
        type="bar",
        name="summary",
        height=1000,
        orientation="h"
    ) %>% layout(xaxis=list(title="Ratio of Chemicals with Annotation", range=c(0, 1)), yaxis = list(title=title, categoryorder = "total ascending"))
    return(p)
}