library(iTALK)
library(dplyr)
library(Seurat)
library(SeuratDisk)
library(Matrix)

# system("DISPLAY='' R --vanilla")

# so_combined = readRDS("S.combined-536-547-552-NT.rds")
# names(so_combined@meta.data)


# print("umap")
# p = DimPlot(so_combined, group.by = "integrated_snn_res.0.1" , label = TRUE)
# png('UMAP_res_01.png')
# print(p)
# dev.off()


# count_mtx = so_combined@assays$RNA@counts
# class(count_mtx)
# dim(count_mtx)
# transpose à faire sur le cluster car très long
# print("transpose")
# count_mtx = Transpose(count_mtx)

# count_mtx = as.matrix(count_mtx)
# dim(count_mtx)
# cell_type = so_combined$integrated_snn_res.1
# cell_type = paste("Cluster_",as.character(cell_type), sep="")


# compare_group = so_combined$treatment

# count_mtx = cbind(count_mtx, cell_type)

# count_mtx = cbind(count_mtx, compare_group)
# saveRDS(count_mtx, "36_47_52_count_mtx.rds")


# count_mtx = readRDS('36_47_52_count_mtx.rds')
# count_mtx = as.data.frame(count_mtx)

# on récupère la colonne comparaison
# compare_group = count_mtx[, "compare_group"]

# count_mtx = count_mtx


# for (i in 1:(dim(count_mtx)[2]-2)){
    # count_mtx[,i] = as.numeric(count_mtx[, i])
# }

# saveRDS(count_mtx, "num_count_mtx.rds")
# so = readRDS("S.combined-536-547-552-NT.rds")
# cell_type = so$integrated_snn_res.0.1
# count_mtx = readRDS('num_count_mtx.rds')
# count_mtx[,"cell_type"] = cell_type
count_mtx = readRDS("num_count_mtx.rds")
# highly_exprs_genes<-rawParse(count_mtx[,1:(dim(count_mtx)[2]-1)],top_genes=50,stats='mean')
# saveRDS(highly_exprs_genes, "highly_exprs_genes.rds")
highly_exprs_genes = readRDS("highly_exprs_genes.rds")
# find the ligand-receptor pairs from highly expressed genes
comm_list<-c('growth factor','other','cytokine','checkpoint')
cell_col<-structure(c('#e9f35f','#e87474','#52c63b','#a53bc6', '#3B54C6'),names=unique(count_mtx$cell_type))
# cell_col

print("on est là")
res<-NULL
for(comm_type in comm_list[1:2]){
    res_cat<-FindLR(highly_exprs_genes,datatype='mean count',comm_type=comm_type)
    res_cat<-res_cat[order(res_cat$cell_from_mean_exprs*res_cat$cell_to_mean_exprs,decreasing=T),]
    # plot by ligand category
    # overall network plot
    # NetView(res_cat,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
    # top 20 ligand-receptor pairs
	pdf(paste(comm_type,".pdf"))
    
	LRPlot(res_cat[1:20,],datatype='mean count',cell_col=cell_col,link.arr.lwd=res_cat$cell_from_mean_exprs[1:20],link.arr.width=res_cat$cell_to_mean_exprs[1:20])
    title(comm_type)
	
	
    dev.off()
    res<-rbind(res,res_cat)
}
# res<-res[order(res$cell_from_mean_exprs*res$cell_to_mean_exprs,decreasing=T),][1:20,]
# NetView(res,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
# LRPlot(res[1:20,],datatype='mean count',cell_col=cell_col,link.arr.lwd=res$cell_from_mean_exprs[1:20],link.arr.width=res$cell_to_mean_exprs[1:20])

#  Ajout du groupe traitement (colonne "compare_group")
# count_mtx$compare_group = compare_group

# find DEGenes of regulatory T cells and NK cells between these 2 groups
# deg_Tgd<-DEG(count_mtx %>% filter(cell_type=='Tgd Delta1'),method='Wilcox',contrast=c("S63","NT"))
# deg_NKr<-DEG(count_mtx %>% filter(cell_type=='NK au repos'),method='Wilcox',contrast=c("S63","NT"))
# deg_NKp<-DEG(count_mtx %>% filter(cell_type=='NK prolif.'),method='Wilcox',contrast=c("S63","NT"))
# deg_NKr_NKp = rbind(deg_NKr, deg_NKp)

# saveRDS(deg_Tgd, "deg_Tgd.rds")
# saveRDS(deg_NKr_NKp, "deg_NKr_NKp.rds")

# deg_Tgd = readRDS("deg_Tgd.rds")
# deg_NKr_NKp = readRDS('deg_NKr_NKp.rds')
# par(mfrow=c(1,1))
# res<-NULL
# for(comm_type in comm_list){
    # res_cat<-FindLR(deg_Tgd,deg_NKr_NKp,datatype='DEG',comm_type=comm_type)
    # res_cat<-res_cat[order(res_cat$cell_from_logFC*res_cat$cell_to_logFC,decreasing=T),]
    # plot by ligand category
    # if(nrow(res_cat)==0){
        # next
    # }else if(nrow(res_cat>=20)){
        # LRPlot(res_cat[1:20,],datatype='DEG',cell_col=cell_col,link.arr.lwd=res_cat$cell_from_logFC[1:20],link.arr.width=res_cat$cell_to_logFC[1:20])

    # }else{
        # LRPlot(res_cat,datatype='DEG',cell_col=cell_col,link.arr.lwd=res_cat$cell_from_logFC,link.arr.width=res_cat$cell_to_logFC)

    # }
    # NetView(res_cat,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)

    # title(comm_type)

    # res<-rbind(res,res_cat)
# }
# if(is.null(res)){
    # print('No significant pairs found')
# }else if(nrow(res)>=20){
    # res<-res[order(res$cell_from_logFC*res$cell_to_logFC,decreasing=T),][1:20,]
    # NetView(res,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
    # LRPlot(res[1:20,],datatype='DEG',cell_col=cell_col,link.arr.lwd=res$cell_from_logFC[1:20],link.arr.width=res$cell_to_logFC[1:20])
# }else{
    # NetView(res,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
    # LRPlot(res,datatype='DEG',cell_col=cell_col,link.arr.lwd=res$cell_from_logFC,link.arr.width=res$cell_to_logFC)
# }

# head(res)
# Sabine mail 24/06/22
# iTalk sur NT ou S63 seul
# NT res 0.1 /:/  NK=6 et Treg =7
# S63 res 0.2 : "NK ou Tgd"=2 / "NK ou Tgd cytotox"=5 / CD4 Treg=3 / CD8 Treg=6

# so_NT = readRDS('../../RDS/NT_intFromScratch.rds')
# so_S63 = readRDS('../../RDS/S63_intFromScratch.rds')

# subsetting NK / Treg
# names(so_NT@meta.data)
# sub_NT = subset(so_NT, subset = integrated_snn_res.0.1 == 6 | integrated_snn_res.0.1 == 7)

# clust_nb = sub_NT$integrated_snn_res.0.1
# func_1 = function(x){
    # return( if(x==6) "NK" else "Treg" )
# } 
# clust_id = unlist(lapply(clust_nb, func_1))
# sub_NT$clust_id = clust_id


# sub_S63 = subset(so_S63, subset = integrated_snn_res.0.2 == 2 | 
                    # integrated_snn_res.0.2 == 5 | 
                    # integrated_snn_res.0.2 == 3 | 
                    # integrated_snn_res.0.2 == 6 )

# clust_nb = sub_S63$integrated_snn_res.0.2
# func_2 = function(x){
    # return( if(x==6) "CD8_Treg" else if(x==2) "NK/Tgd" else if(x==5) "NK/Tgd_cyto" else  "CD4_Treg")
# } 

# clust_id = unlist(lapply(clust_nb, func_2))
# sub_S63$clust_id = clust_id

# count_mtx_NT = sub_NT@assays$RNA@counts
# count_mtx_S63 = sub_S63@assays$RNA@counts


# count_mtx_NT = Transpose(count_mtx_NT)
# count_mtx_S63 = Transpose(count_mtx_S63)

# saveRDS(count_mtx_NT, "../../RDS/iTalk_count_mtx_NT.rds")
# saveRDS(count_mtx_S63, "../../RDS/iTalk_count_mtx_S63.rds")

# count_mtx_NT = as.matrix(count_mtx_NT)
# cell_type = as.character(sub_NT$clust_id)
# count_mtx_NT = cbind(count_mtx_NT, cell_type)

# count_mtx_S63 = as.matrix(count_mtx_S63)
# cell_type = as.character(sub_S63$clust_id)
# count_mtx_S63 = cbind(count_mtx_S63, cell_type)

# char_to_num = function(mtx){
    # for (i in 1:(dim(mtx)[2]-1)){
        # mtx[,i]= as.numeric(mtx[,i])
    # }
    # return(mtx)
# }

# count_mtx_NT = as.data.frame(count_mtx_NT)
# count_mtx_NT = char_to_num(count_mtx_NT)
# count_mtx_S63 = as.data.frame(count_mtx_S63)
# count_mtx_S63 = char_to_num(count_mtx_S63)


# iTalk_LRPlot = function(mtx){

    # highly_exprs_genes<-rawParse(mtx,top_genes=50,stats='mean')

    # find the ligand-receptor pairs from highly expressed genes
    # comm_list<-c('growth factor','other','cytokine','checkpoint')
    # cell_col<-structure(c('#e9f35f','#e87474','#52c63b', '#3f4ae0'),names=unique(mtx$cell_type))
    # cell_col
    # name_mtx = deparse(substitute(mtx))

    # res<-NULL
    # for(comm_type in comm_list){
        # res_cat<-FindLR(highly_exprs_genes,datatype='mean count',comm_type=comm_type)
        # res_cat<-res_cat[order(res_cat$cell_from_mean_exprs*res_cat$cell_to_mean_exprs,decreasing=T),]
        #plot by ligand category
        #overall network plot
        # NetView(res_cat,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
        #top 20 ligand-receptor pairs
        # png(paste(name_mtx, "_",comm_type,".png", sep=""))
        # LRPlot(res_cat[1:20,],datatype='mean count',cell_col=cell_col,link.arr.lwd=res_cat$cell_from_mean_exprs[1:20],link.arr.width=res_cat$cell_to_mean_exprs[1:20])
        # title(comm_type)
        # dev.off()
        # res<-rbind(res,res_cat)
    # }
# }
# iTalk_LRPlot(count_mtx_NT)
# iTalk_LRPlot(count_mtx_S63)
