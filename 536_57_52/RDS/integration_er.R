library(Seurat)

list_rds=list.files()[grep("QC-", list.files())]


list_so = c()
list_names = c()
for (rds in list_rds){
	name_rds = gsub("-", "_", rds, fixed=T)
	name_rds = gsub(".rds", "", name_rds)
	assign(name_rds, readRDS(rds))
	temp = eval(parse(text= name_rds))
	c_names = colnames(temp)
	c_names = paste(c_names, "_", name_rds, sep="")
	temp = RenameCells(temp, new.names = c_names)
	list_so = append(list_so, temp)
	list_names = append(list_names, name_rds)
	}
names(list_so) = list_names

wu = readRDS("Wu_epithelial.rds")

##### subseting by subtype
wu_er = subset(wu, subset= subtype == "ER+")

##### split by orig.ident
wu_er_split = SplitObject(wu_er, split.by="orig.ident")

list_so_er = append(list_so, wu_er_split)

saveRDS(list_so_er, "list_so_er.rds")

list_so_er = lapply(X= list_so_er, FUN = function(x){
					x = NormalizeData(x)
					x = FindVariableFeatures(x)
					})

features = SelectIntegrationFeatures(list_so_er)

list_so_er = lapply(list_so_er, FUN = function(x){
					x = ScaleData(x, features = features)
					x = RunPCA(x, features = features)
					})
					
anchors_rpca_er = FindIntegrationAnchors(object.list = list_so_er,
			anchor.features = features, reduction = "rpca")
			
saveRDS(anchors_rpca_er, "anchors_rpca_er.rds")

int_set = IntegrateData(anchorset = anchors_rpca_er)

saveRDS(int_set, "int_set_rpca_er.rds")