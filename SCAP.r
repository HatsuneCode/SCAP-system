## 1. read counts ##
counts = read.table('raw_counts_All_Samples.txt', sep = '\t', header = T, row.names = 1)
counts = counts[rowSums(counts) > 0,]
## check gene symbols
RNAseq.checkDupRow = function(expr, method = 'mean') {
  gene  = sub('^ENS.*?_', '', rownames(expr))
  dgene = unique(gene[duplicated(gene)])
  if (length(dgene)) {
    expr1 = expr[!gene %in% dgene,]
    expr2 = do.call(rbind, lapply(dgene, function(g) {
      e = expr[gene == g, , drop = F]
      if (method == 'mean')
        t(setNames(data.frame(colMeans(e)), g))
    }))
    expr = rbind(expr1, expr2)
    rm(expr1, expr2)
  }
  rm(gene, dgene)
  rownames(expr) = sub('^ENS.*?_', '', rownames(expr))
  expr
}
counts = RNAseq.checkDupRow(counts)
counts = round(counts)

## 2. DEG identified by DESeq2 ##
RNAseq.DESeq2 = function(expr, pos = NULL, neg = NULL, name = NULL, exp_cut = 10) {
  suppressMessages(library(DESeq2))
  ## expr split
  exprP = if (length(pos)) expr[, colnames(expr) %in% pos, drop = F] else 
    expr[, !colnames(expr) %in% neg, drop = F]
  exprN = if (length(neg)) expr[, colnames(expr) %in% neg, drop = F] else 
    expr[, !colnames(expr) %in% pos, drop = F]
  ## expr type
  type = paste(paste(colnames(exprP), collapse = ','), 'vs', paste(colnames(exprN), collapse = ',') )
  if (!length(name)) name = type
  message('DEG: ', type)
  ## condition control ~ treatment
  condition = factor( c(rep('Neg', ncol(exprN)), rep('Pos', ncol(exprP))), c('Neg', 'Pos') )
  ## counts extract
  expr  = cbind(exprN, exprP)
  expr  = expr[rowSums(expr) > 0, , drop = F]
  if (length(exp_cut))
    expr = expr[apply(expr, 1, function(i) !any(i < exp_cut) ), ]
  ##
  exprP = expr[, condition == 'Pos', drop = F]
  exprN = expr[, condition == 'Neg', drop = F]
  ## meta information
  meta = data.frame(row.names = colnames(expr), condition)
  ## DESeq2 run
  dds = DESeqDataSetFromMatrix(countData = expr, colData = meta, design = ~ condition)
  dds = DESeq(dds)
  dds = data.frame(results(dds), check.names = F)
  ## output
  data.frame(p_val = dds$pvalue, avg_log2FC = dds$log2FoldChange, 
             pct.1 = apply(exprP, 1, function(i) sum(i > 0)/ncol(exprP) ),
             pct.2 = apply(exprN, 1, function(i) sum(i > 0)/ncol(exprN) ),
             p_val_adj = dds$padj, gene = rownames(dds), 
             average = rowMeans(expr), median = apply(expr, 1, median), 
             posAvg = rowMeans(exprP), posMed = apply(exprP, 1, median),
             negAvg = rowMeans(exprN), negMed = apply(exprN, 1, median),
             type = name, 
             upDown = factor(ifelse(dds$log2FoldChange > 0, 'Up', 'Down'), c('Up', 'Down')), 
             row.names = NULL )
}
group = sub('\\..*', '', colnames(counts))
DEG = do.call(rbind, lapply(unique(group), function(pos)
  do.call(rbind, lapply('Ctrl', function(neg) {
    if (pos != neg) 
      RNAseq.DESeq2(counts, 
                    colnames(counts)[group == pos], 
                    colnames(counts)[group == neg], 
                    paste(pos, 'vs', neg) )
  })) ))
write.table(DEG, '2.DEG.bulk.xls', sep = '\t', row.names = F)
## stat DEG ##
fc = 1
pv = .05
deg = DEG[which(DEG$p_val < pv), ]
deg = deg[abs(deg$avg_log2FC) > fc,]
deg = deg[grep('vs Ctrl', deg$type), ]
pd = table(deg$type, deg$upDown)
pd[,2] = -pd[,2]
p  = Heatmap(pd, name = 'nDEG', rect_gp = gpar(col = 'grey90'), 
             #split = sub(' .*', '', rownames(pd)),
             cluster_rows = F, cluster_columns = F,
             row_names_side = 'left', column_names_side = 'top', column_names_rot = 0,
             column_title = paste('DEG: Log2FC', fc, 'P', pv), 
             column_title_gp = gpar(fontsize = 17, fontfamily = 'serif'),
             cell_fun = function(j, i, x, y, w, h, col) grid.text(abs(pd[i,j]), x, y) )
pdf('2.statDEG.pdf', w = 4, h = 2.5); draw(p); dev.off()

## 3. PCA ##
## Normalize counts by DESeq2
RNAseq.Normalize = function(expr, log2 = T, method = 'DESeq2') {
  if (method == 'DESeq2') {
    suppressMessages(library(DESeq2))
    dds  = DESeqDataSetFromMatrix(
      countData = expr,
      colData = data.frame(row.names = colnames(expr), 
                           samples = factor(colnames(expr))),
      design = ~ samples)
    dds  = estimateSizeFactors(dds)
    expr = counts(dds, normalize = T)
  }
  if (log2) expr = log2(expr + 1)
  expr
}
counts = RNAseq.Normalize(counts)
## Run PCA
PCA = function(expr, n.varGene = 3e+3, ...) {
  if (length(n.varGene))
    expr = expr[RNAseq.varGene(expr, n.varGene = n.varGene, ...),]
  pca = data.frame(prcomp(t(expr), ...)$x, check.rows = F)
  pca$sample = rownames(pca)
  pca$sample = factor(pca$sample, pca$sample)
  pca
}
RNAseq.varGene = function(expr, n.varGene = 3000, outFig = NULL) {
  suppressMessages(library(Seurat))
  suppressMessages(library(ggplot2))
  genes = grep('_', rownames(expr), value = T)
  if (length(genes))
    genes = data.frame(Raw = genes, Pro = gsub('_', '-', genes))
  obj = suppressWarnings(CreateSeuratObject(expr))
  obj = FindVariableFeatures(obj, nfeatures = n.varGene, verbose = F)
  p   = VariableFeaturePlot(obj)
  if (length(outFig)) ggsave(outFig, p, w = 6, h = 5)
  var = suppressWarnings(VariableFeatures(obj))
  rm(obj)
  if (length(genes)) {
    comm = intersect(var, genes$Pro)
    if (length(comm)) var[match(comm, var)] = genes$Raw[match(comm, genes$Pro)]
  }
  var
}
pca = PCA(counts, outFig = '3.PCA.nVar.png')
pca$group = sub('\\..*', '', pca$sample)
library(ggrepel)
p = ggplot(pca, aes(PC1, PC2)) + 
  geom_point(aes(color = group)) +
  geom_text_repel(label = pca$sample, family = 'serif') +
  labs(title = 'PCA') +
  theme_classic() + 
  theme(plot.title = element_text(hjust = .5), text = element_text(size = 14, family = 'serif'))
ggsave('3.PCA.pdf', p, w = 6, h = 5)

## 4. Sample correlation ##
library(ComplexHeatmap)
library(circlize)
cntV = counts[RNAseq.varGene(counts, outFig = '4.Cor.nVar.png'),]
cor = cor(cntV, method = 'spearman')
for (i in 1:nrow(cor)) cor[i,i] = NA
sp  = sub('\\..*', '', colnames(cor))
p = Heatmap(cor, name = 'Spearman',
            cluster_rows = F, cluster_columns = F,
            col = colorRamp2(c(.8, 1), c('white', 'red')),
            na_col = 'white', rect_gp = gpar(col = 'grey90'),
            split = sp, column_split = sp,
            show_row_names = F, show_column_names = F,
            row_dend_side = 'right', column_dend_side = 'bottom',
            cell_fun = function(j, i, x, y, w, h, col) if (!is.na(cor[i,j])) 
              grid.text(round(cor[i,j], 2), x, y) )
pdf('4.Cor.pdf', w = 8, h = 7); draw(p); dev.off()

## 5. GSEA ##
fGSEA = function(gene, sig = NULL, scoreType = 'std', minSize = 2, maxSize = 500, type = NULL, species = 'Homo sapiens') {
  # species: Homo sapiens / Mus musculus
  suppressMessages(library(fgsea))
  suppressMessages(library(msigdbr))
  if ('KEGG' %in% type) {
    kegg  = msigdbr(species, 'C2', 'KEGG')
    keggn = setNames(lapply(unique(kegg$gs_name), function(i)
      unique(as.character(kegg$gene_symbol)[kegg$gs_name == i])), unique(kegg$gs_name))
    rm(kegg)
    sig = keggn
  }
  if ('GOBP' %in% type) {
    bp   = msigdbr(species, 'C5', 'BP')
    bpn  = setNames(lapply(unique(bp$gs_name), function(i)
      unique(as.character(bp$gene_symbol)[bp$gs_name == i] )), unique(bp$gs_name))
    rm(bp)
    sig = bpn
  }
  if ('GOCC' %in% type) {
    cc   = msigdbr(species, 'C5', 'CC')
    ccn  = setNames(lapply(unique(cc$gs_name), function(i)
      unique(as.character(cc$gene_symbol)[cc$gs_name == i] )), unique(cc$gs_name))
    rm(cc)
    sig = ccn
  }
  if ('GOMF' %in% type) {
    mf   = msigdbr(species, 'C5', 'MF')
    mfn  = setNames(lapply(unique(mf$gs_name), function(i)
      unique(as.character(mf$gene_symbol)[mf$gs_name == i] )), unique(mf$gs_name))
    rm(mf)
    sig = mfn
  }
  set.seed(1)
  ## run gsea analysis by fgsea
  gsea = fgsea(sig, gene, minSize = minSize, maxSize = maxSize, scoreType = scoreType)
  gsea$gene = unlist(lapply(gsea$leadingEdge, function(i) paste(i, collapse = ', ') ))
  data.frame(gsea[, c('pathway', 'NES', 'ES', 'pval', 'padj', 'gene')])
}
GSEA = do.call(rbind, lapply(c('KEGG', 'GOBP', 'GOCC', 'GOMF'), function(i) {
    gsea = do.call(rbind, lapply(unique(DEG$type), function(g) {
      tmp = DEG[DEG$type == g,]
      gene = sort(setNames(tmp$avg_log2FC, tmp$gene), T)
      df = fGSEA(gene, type = i, species = 'Homo sapiens')
      df$group = g
      df
    }))
    gsea$type = i
    gsea
}))
write.table(GSEA, '5.GSEA.xls', sep = '\t', row.names = F)

## 6. Plot target pathways ##
tgs = c('GOBP_CELLULAR_RESPONSE_TO_VIRUS', 
        'GOBP_INFLAMMATORY_RESPONSE_TO_ANTIGENIC_STIMULUS', 
        'GOBP_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY', 
        'GOBP_REGULATION_OF_DEFENSE_RESPONSE_TO_VIRUS', 
        'GOMF_CYTOKINE_ACTIVITY', 
        'GOBP_I_KAPPAB_KINASE_NF_KAPPAB_SIGNALING', 
        'GOBP_NEGATIVE_REGULATION_OF_INTRACELLULAR_SIGNAL_TRANSDUCTION', 
        'GOBP_REGULATION_OF_RESPONSE_TO_CYTOKINE_STIMULUS', 
        'GOBP_POSITIVE_REGULATION_OF_IMMUNE_SYSTEM_PROCESS', 
        'GOBP_INNATE_IMMUNE_RESPONSE', 
        'GOBP_PATTERN_RECOGNITION_RECEPTOR_SIGNALING_PATHWAY', 
        'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION', 
        'GOBP_PHAGOCYTOSIS', 
        'GOBP_CYTOKINE_MEDIATED_SIGNALING_PATHWAY', 
        'GOBP_IMMUNE_EFFECTOR_PROCESS' )
gsea = GSEA[GSEA$pathway %in% tgs,]
## plot p-value heatmap
library(reshape2)
pd = acast(gsea, pathway ~ group, value.var = 'pval', fill = 1)
pd = -log10(pd)
library(ComplexHeatmap); library(circlize)
p = Heatmap(pd, name = '-log10(P-value)', rect_gp = gpar(col = 'grey90'),
            row_names_max_width = max_text_width(rownames(pd)),
            show_row_dend = F, row_names_side = 'left', column_names_side = 'top',
            column_order = c('SCAP vs Ctrl', 'AnaPIE vs Ctrl', 'Linear vs Ctrl'),
            col = colorRamp2(c(0, -log10(.001)), c('white', 'red')) );p
pdf('6.GSEA.target.pdf', w = 7.5, h = 4.5); draw(p); dev.off()
