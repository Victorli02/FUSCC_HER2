
heatmap_map <- ComplexHeatmap::pheatmap(choose_matrix,
                                         annotation_col =  group_list,
                                         column_split = group_list[,1],
                                         annotation_colors = list( "subtype" = subtype_color),
                                         show_rownames = T, show_colnames = F, 
                                         annotation_legend = T, cluster_cols = F, cluster_rows = F,
                                         border_color = NA, 
                                         legend = T, 
                                         color = mycolors,
                                         breaks= mybreaks)

