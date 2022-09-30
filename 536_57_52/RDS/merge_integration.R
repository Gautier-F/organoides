library(Seurat)

list_rds=list.files()[grep("QC-", list.files())]


# creer la liste d'objet

# list_so = c()
# list_names = c()
# for (rds in list_rds){
	# name_rds = gsub("-", "_", rds, fixed=T)
	# name_rds = gsub(".rds", "", name_rds)
	# assign(name_rds, readRDS(rds))
	# temp = eval(parse(text= name_rds))
	# c_names = colnames(temp)
	# c_names = paste(c_names, "_", name_rds, sep="")
	# temp = RenameCells(temp, new.names = c_names)
	# list_so = append(list_so, temp)
	# list_names = append(list_names, name_rds)
	# }
# names(list_so) = list_names
# print("list_so done")
# wu_so = readRDS("Wu_epithelial.rds")
# wu_list  = SplitObject(wu_so, split.by="orig.ident")

# list_so = append(list_so, wu_list)

# saveRDS(list_so, "list_all_so.rds")

# list_so = readRDS("list_all_so.rds")

#integration (Tips for integrating large datasets, vignette seurat )

# features = SelectIntegrationFeatures(list_so)
# list_so = lapply(X= list_so, FUN = function(x){
					# x = NormalizeData(x)
					# x = FindVariableFeatures(x)
					# })


# features = SelectIntegrationFeatures(list_so)

# list_so = lapply(list_so, FUN = function(x){
					# x = ScaleData(x, features = features)
					# x = RunPCA(x, features = features)
					# })

# anchors_rpca = FindIntegrationAnchors(object.list = list_so,
			# anchor.features = features, reduction = "rpca")

# anchors_cca = FindIntegrationAnchors(object.list = list_so, anchor.features = features)
# print("anchors done")
# print("save anchors")
# saveRDS(anchors_rpca, "anchors_rpca.rds")
# saveRDS(anchors_cca, "anchors_cca.rds")

anchors = readRDS("anchors_rpca.rds")

int_set = IntegrateData(anchorset = anchors)
# int_set <- ScaleData(int_set, verbose = FALSE)
# int_set <- RunPCA(int_set, npcs = 30, verbose = FALSE)
# int_set <- RunUMAP(int_set, reduction = "pca", dims = 1:30)
# int_set <- FindNeighbors(int_set, reduction = "pca", dims = 1:30)
# int_set <- FindClusters(int_set, resolution = c( 0.5)
 # saveRDS(int_set, "int_set_rpca.rds")					


