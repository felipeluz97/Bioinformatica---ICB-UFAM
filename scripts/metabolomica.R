# --- Preparação: instalar / carregar pacotes -----------------------------
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
pkgs_cran <- c("tidyverse", "limma")
pkgs_bioc  <- c("xcms", "MSnbase", "CAMERA", "MetaboAnalystR")

for(p in pkgs_cran) if(!requireNamespace(p,quietly=TRUE)) install.packages(p)
for(p in pkgs_bioc) if(!requireNamespace(p,quietly=TRUE)) BiocManager::install(p)

library(tidyverse)
library(xcms)        # processamento LC-MS (peak detection, align)
library(MSnbase)     # leitura mzML,mzXML
library(CAMERA)      # anotação isotopologues/adducts
library(MetaboAnalystR) # análises estatísticas / funcionais
library(limma)       # alternativa para DE (linear models)

# --- 1) Ler dados brutos LC-MS (mzML) e detectar picos --------------------
# Ajuste o diretório e padrões conforme seus arquivos
raw_dir <- "data/mzML_files/"
fnames <- list.files(raw_dir, pattern = "\\.(mzML|mzXML)$", full.names = TRUE)
raw <- readMSData(fnames, pdata = NULL, mode = "onDisk")  # MSnbase

# Reproducible typical xcms pipeline (centWave peak picking)
cwp <- CentWaveParam(ppm=15, peakwidth = c(5,30), snthr = 10)
xdata <- findChromPeaks(raw, param = cwp)

# Grouping / retention time correction / fill peaks
pdp <- PeakDensityParam(sampleGroups = rep(1, length(fnames)), bw = 5)
xdata <- groupChromPeak(xdata, param = pdp)
# Rt correction & re-group (you may prefer obiwarp)
xdata <- adjustRtime(xdata, param = ObiwarpParam(binSize=0.6))
xdata <- groupChromPeak(xdata, param = pdp)
xdata <- fillChromPeaks(xdata)

# Convert xcms results to a feature table: features x samples (intensity)
featureValues <- featureValues(xdata, method = "medret") # matrix: features x samples
finfo <- featureDefinitions(xdata) # info de m/z, rt, etc.

# --- 2) Anotação com CAMERA (opcional) -----------------------------------
xs <- as(xdata, "xcmsSet") # algumas versões/compatibilidades pedem conversão
an <- xsAnnotate(xs)
an <- groupFWHM(an)
an <- findIsotopes(an)
an <- findAdducts(an)
annot <- getPeaklist(an) # tabela com anotações (adducts/isotopes)

# --- 3) Preparar arquivos para MetaboAnalystR (salvar CSVs) ---------------
# MetaboAnalyst espera: feature table (rows = features, cols = samples) e metadata
feat_tbl <- as.data.frame(featureValues) %>%
  rownames_to_column("FeatureID") %>%
  bind_cols(finfo %>% select(mz = mz, rt = rt)) # anexar mz/rt (ajuste nomes conforme finfo)

# Se tiver planilha de metadados com colunas: SampleID, Group (ex: Control/Case), etc.
# metadata <- read.csv("data/sample_metadata.csv") 

write.csv(feat_tbl, "metabo_feature_table.csv", row.names = FALSE)
# --- Crie um metadata.csv com pelo menos: SampleID, Group
# write.csv(metadata, "metabo_metadata.csv", row.names = FALSE)

# --- 4) Workflow MetaboAnalystR: carregar e checar dados -------------------
# Crie um diretório de trabalho para salvar outputs do MetaboAnalyst
dir.create("MA_results", showWarnings = FALSE)
setwd("MA_results")

# Inicializar objeto mSet (ex.: 'pktable' para peak table untargeted; tipo 'stat' para análises estatísticas)
# Tipo de dados: "pktable" (peak table), formato: "rowu" (features nas linhas, amostras nas colunas)
mSet <- InitDataObjects("pktable", "stat", FALSE) 

# Read.TextData lê: (mSetObj, filePath, format, lbl.type)
# format = "rowu" = rows = features, first column feature names; lbl.type = 1 se metadados com labels
mSet <- Read.TextData(mSetObj = mSet, filePath = "../metabo_feature_table.csv", format = "rowu", lbl.type = 1)

# Se tiver metadata separado: use Read.PhenoData or Read.TextData com lbl.type adequado.
# Opcional: checar integridade
mSet <- SanityCheckData(mSet)

# --- 5) Tratamento de missing, filtragem, normalização ---------------------
# Substituir zeros/NA por um pequeno valor (ReplaceMin), filtrar features pouco informativas
mSet <- ReplaceMin(mSet)                 # substitui zeros/NA por min detectável
mSet <- FilterVariable(mSet, "none", "F", 0.5) # Ex.: filtrar features com >50% missing (ajuste se desejar)

# Normalização (várias opções: "NULL", "SumNorm", "MedianNorm", "ProbNorm", etc.)
# Normalization function aceita: sampleNorm (row-wise), transNorm, scaleNorm
mSet <- Normalization(mSet, "PQN", "LogNorm", "Pareto", ratio=FALSE, ref=NULL) 
# Exemplo: PQN para amostras, log transform, pareto scaling — ajuste conforme seu dado

# --- 6) Análises exploratórias: PCA, PLS-DA -------------------------------
mSet <- PCA.Anal(mSet)          # PCA (cria arquivos de plot & output)
mSet <- PLSR.Anal(mSet)         # PLS-DA / PLSR (use OPLSR.Anal para OPLS-DA se disponível)

# Plots gerados por MetaboAnalystR serão salvos em mSet$imgSet — exibir no R:
# Ver score plot
print(mSet$imgSet$stat$PCA.score)  # ou ver o ficheiro PNG gerado no diretório

# Clustering heatmap (todas as features, ou as top variáveis)
mSet <- Heatmap.Pairwise(mSet, method="euclidean", clust="ward", scale="row")
# ou
mSet <- HeatmapDendro(mSet, "euclidean", "ward", scale="row")


mSet <- VolcanoPlot(mSet, fc.thr = 2, p.thr = 0.05, 
                    useFDR = TRUE, paircomp = NULL)

# --- 7) Testes diferenciais com MetaboAnalystR (t-test / ANOVA) -------------
# t-test (para 2 grupos)
mSet <- Ttests.Anal(mSet, nonpar = FALSE, threshp = 0.05, paired = FALSE, equal.var = TRUE)
# ANOVA para >2 agrupamentos:
mSet <- ANOVA.Anal(mSet, "fisher", thresh = 0.05)  # exemplo, veja documentação para parâmetros

# Os resultados estão gravados em arquivos CSV no diretório de resultados e em:
# mSet$dataSet$norm -> matriz normalizada
# mSet$analSet$tt$raw.pvalues etc (dependendo da versão)

# --- 8) Correção por múltiplos testes (FDR) - exemplo com p.adjust ----------------
# Extrair tabela de t-test gerada (substitua path/nome conforme a versão do MetaboAnalyst)
tt_file <- list.files(pattern = "t_test.csv", recursive = TRUE, full.names = TRUE)
if(length(tt_file)>0){
  tt <- read.csv(tt_file[1], stringsAsFactors = FALSE)
  tt$adj.p <- p.adjust(tt$p.value, method = "BH")
  write.csv(tt, "t_test_withFDR.csv", row.names = FALSE)
}

# --- 9) Alternativa (recomendada) — usar limma para DE (mais controle) -----
# Supondo que você tenha:
# - matrix 'norm_mat' com features (rows) x samples (cols)
# - metadata dataframe com SampleID, Group
norm_mat <- as.matrix(mSet$dataSet$norm) # verifique o caminho exato no seu mSet
# Carregar metadata (exemplo)
# metadata <- read.csv("../metabo_metadata.csv", stringsAsFactors = FALSE)
# Garanta que colnames(norm_mat) correspondem a metadata$SampleID

# Construir design model (exemplo 2 grupos: Control vs Case)
metadata$Group <- factor(metadata$Group)
design <- model.matrix(~0 + metadata$Group)
colnames(design) <- levels(metadata$Group)

# Limma pipeline
fit <- lmFit(norm_mat, design)
cont.matrix <- makeContrasts(Case_vs_Control = Case - Control, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
res <- topTable(fit2, number = Inf, adjust.method = "BH", sort.by = "P")
write.csv(res, "limma_results_all_features.csv", row.names = TRUE)

# --- 10) Análise funcional / vias (mummichog / MetPA) -----------------------
# Se houver m/z identificados (ou usou annotações CAMERA/HMDB), MetaboAnalystR tem módulos de Pathway Analysis
# Exemplo: se tiver IDs HMDB ou m/z -> use MS Peaks to Pathway (mummichog)
# mSet <- PerformPSEA(mSet, method = "mummichog", pval.cutoff = 0.05)  # ajuste conforme versão
# Consulte vignettes oficiais para parâmetros detalhados.

# ----------------- FIM -----------------
cat("Workflow concluído. Cheque o diretório 'MA_results' para plots e tabelas.\n")

# --- 7b) Selecionar metabólitos diferenciais ---------------------------
# Suponha que você rodou o t-test ou limma
tt_res <- read.csv("t_test_withFDR.csv", stringsAsFactors = FALSE)

# Filtro de significativos (ajuste thresholds conforme seu estudo)
sig_metabs <- tt_res %>%
  filter(adj.p < 0.05 & abs(log2(FC)) > 1)

# Salvar lista de IDs (depende do que você tem: HMDB, KEGG, m/z)
write.csv(sig_metabs, "sig_metabolites.csv", row.names = FALSE)

# --- 11) Enriquecimento de vias (MetaboAnalystR) -----------------------
# Se você tiver IDs de metabolitos (ex: HMDB), use MSEA:
mSet <- InitDataObjects("mset", "msetora", FALSE)

# Ler lista de IDs significativos
mSet <- Setup.MapData(mSet, "sig_metabolites.csv")

# Definir organismo-alvo (ex.: "hsa" para humano, "mmu" para mouse)
mSet <- CrossReferencing(mSet, "name")     # ou "hmdb", "kegg", "pubchem"
mSet <- CreateMappingResultTable(mSet)

# Enriquecimento via Over-Representation Analysis
mSet <- SetMetabolomeFilter(mSet, F)       # usa todos os metabólitos
mSet <- SetCurrentMsetLib(mSet, "smpdb_pathway", 2)  # biblioteca de vias
mSet <- CalculateOraScore(mSet, "rbc")     # "rbc" = referência por background completo

# Resultados (tabela de vias enriquecidas)
enrich_table <- mSet$analSet$ora.mat
write.csv(enrich_table, "pathway_enrichment_results.csv")

# --- 11b) Se for untargeted (apenas picos m/z): usar mummichog ----------
# Partindo do objeto mSet de estatística
mSet <- PerformPSEA(mSet, method = "mummichog", permNum = 100, pval.cutoff = 0.05)

# Resultados
pathway_res <- mSet$analSet$mummi.resTable
write.csv(pathway_res, "mummichog_results.csv")

