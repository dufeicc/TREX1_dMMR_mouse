#GSEA
rm(list = ls())
options(stringsAsFactors = F)library(clusterProfiler)
library(org.Mm.eg.db)
library(GseaVis)
library(stringr)
library(GSEABase)
library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(ggsci)
library(fs)

env <- rlang::new_environment()
env$pdi <- "./_i" %>% dir_create() %>% path_real()
env$pdo <- "./_o" %>% dir_create() %>% path_real()
env$pdt <- "./tmp" %>% dir_create() %>% path_real()

# all_mt <- read.delim("./_i/total_count.genename_mod.csv",sep=",",rownames=14)
D_M <- read.delim("./_i/DKO_MLH1.genename.xls") %>% distinct(GeneName,.keep_all = TRUE) %>% filter(Biotype=="protein_coding")
T_W <- read.delim("./_i/TREX1_WT.genename.xls") %>% distinct(GeneName,.keep_all = TRUE) %>% filter(Biotype=="protein_coding")

D_M_mt <- D_M %>% dplyr:::select(GeneName,Log2FoldChange) %>%
  arrange(desc(Log2FoldChange))%>%column_to_rownames(var="GeneName") %>%dplyr::rename("DKO_vs_MLH1"="Log2FoldChange")
D_M_list <-D_M_mt$DKO_vs_MLH1  
names(D_M_list) <- rownames(D_M_mt)  

T_W_mt <- T_W %>% dplyr:::select(GeneName,Log2FoldChange) %>%
  arrange(desc(Log2FoldChange))%>%column_to_rownames(var="GeneName") %>%dplyr::rename("TREX1_vs_WT"="Log2FoldChange")
T_W_list <-T_W_mt$TREX1_vs_WT  
names(T_W_list) <- rownames(T_W_mt)  

all_glist <- list(T_W_list, D_M_list)

# GSEA
lapply(1:2, function(x){
  ego3 <- gseGO(geneList     = all_glist[[x]],   
                OrgDb        = org.Mm.eg.db,    
                ont          = "BP",            
                keyType      = "SYMBOL",         
                minGSSize    = 10,             
                maxGSSize    = 500,            
                pvalueCutoff = 1,                
                verbose      = TRUE)          
  return(ego3)
}) -> m_gsea_list
tmp <- m_gsea_list[[2]]@result$Description %>% unique()
t <- tmp[tmp %>% str_detect("interferon")]

result_gsea<- m_gsea_list[[2]]@result

sel_path<-c("cellular response to type I interferon",
            "cellular response to interferon-gamma",
            "cellular response to interferon-beta",
            "cellular response to interferon-alpha",
            "regulation of type I interferon production",
            "regulation of type I interferon-mediated signaling pathway",
            "type I interferon production",
            "type I interferon signaling pathway",
            "response to interferon-alpha",
            "response to interferon-beta",
            "response to interferon-gamma",
            "response to type I interferon")
sel_id <- result_gsea %>% filter(Description %in%  sel_path )%>%dplyr::select(ID,Description)

col <- pal_npg("nrc")(10) 

for(i in 1:length(sel_path)){
p<- GSEAmultiGP(gsea_list = m_gsea_list,
            curve.col=col[c(4,5)],
            geneSetID = sel_id[i,"ID"],
            addPval = T,
            pvalX = 0.99,pvalY = 0.99,
            legend.position = c(0.85, 0.65),
            exp_name = c("TREX1_vs_WT","DKO_vs_MLH1"))
pdf(paste0(env$pdo,"/", sel_id[i,"Description"],"_TvsW_DvsM.pdf"),width = 6.5,height =6)
print(p)
dev.off()
} 
for(i in 1:length(sel_path)){
p_tmp <- gseaNb(object = m_gsea_list[[2]],
                curveCol=col[5],
       geneSetID ="cellular response to interferon-beta", #sel_id[i,"ID"],
       subPlot = 2)
pdf(paste0(env$pdo,"/TvsW/", sel_id[i,"Description"],"_TvsW.pdf"),width = 6.5,height =6)
print(p_tmp)
dev.off()
}
p4<- dotplotGsea(data =m_gsea_list[1],
            topn= 100,
            order.by ="NES",
            add.seg =T,
            str.width=100,
            line.col ='orange',
            line.type = 1)
pdf(paste0(env$pdo,"/top100_TvsW_GSEA_GO.pdf"),width =15,height =40)
print(p4)
dev.off()
