library(Seurat)
library(ggplot2)
library(openxlsx)
library(dplyr)


res_all_xlsx = function(so_int, name){
    require(openxlsx)
    require(stringr)
    
    DefaultAssay(so_int) <- "integrated"

    n_dir = paste(name, "_DEG-ALL", sep="")
    if (!dir.exists(n_dir)){
        dir.create(n_dir)
    }

    ind = grep("integrated_", colnames(so_int@meta.data))
    res = colnames(so_int@meta.data)[ind]

    for(r in res) {
        Idents(so_int) = r
        mk = FindAllMarkers(so_int, only.pos = TRUE, min.pct = 0.25,
            logfc.threshold = 0.25)
        mk = mk[, c(7,2,1,3:6)]
        wb = createWorkbook()
        name_res = gsub("integrated_snn_res.0.","",r)

        for(c in levels(mk$cluster)) {
            clust = mk[mk$cluster == c, ]
            clust = clust[order(clust$avg_log2FC, decreasing = T),]
            addWorksheet(wb, paste("Clust.", c, sep=""))
            writeData(wb, sheet=paste("Clust.",c, sep=""),
            x=clust)
        }
        saveWorkbook(wb, file=paste(n_dir, "/", name,"_0.", name_res," .xlsx",sep=""))
    }
}


##### heatmap from deg_xlsx

int = readRDS("RDS/int_set_rpca.rds")
markers = read.xlsx("markers_list.xlsx")

markers %>% group_by(cluster) %>% top_n(n= 10, wt = avg_log2FC) -> top10

pdf("hm_int_res_05_v3.pdf")
	DoHeatmap(int, features = top10$gene) + NoLegend()+
	labs(title = "int_set_res_") + theme(text = element_text(size = 6))
dev.off()



# heatmap_xlsx = function(so, wb, name){
	
	# cluster = names(wb)
	# cluster_list = c()
	# genes_list = c()
	
	# markers = 
	# write.xlsx(markers, "top10_int_set.xlsx")
						
	# png(paste(name,".png", sep = ""))
		# DoHeatmap(so, features = markers$topgene) + NoLegend()+
		# labs(title = name)
    # dev.off()
	# }
	
# markers = function(wb, name){
	# df_clust = read.xlsx(wb, names(wb)[1])
	# for (c in names(wb)[2: length(names(wb))]){
		# df_clust = rbind(df_clust, read.xlsx(wb, sheet = c))
		# }
	# write.xlsx(df_clust, paste(name, ".xlsx", sep=""))
	# }
		
		
# Ditterplot
library(dittoSeq)

pdf("int_freq_tnt.pdf")
dittoBarPlot(int, var="integrated_snn_res.0.5", group.by = "tnt")
dev.off()
