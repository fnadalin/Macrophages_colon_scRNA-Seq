# Created on 2019/08/29

OBJECT <- "object.Robj"

library("Seurat")

object <- eval(parse(text=load(OBJECT)))
object <- UpdateSeuratObject(object)

