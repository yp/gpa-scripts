########################################################################
#
#    Gene/protein analyses
#    Copyright (C) 2015  Yuri Pirola, Raffaella Rizzi
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
########################################################################

##############
##
## This is the main script.
## To ensure reproducibility we distribute with the script
## the composition of the KEGG pathway as of Jan 22, 2015 (downloaded using
## package GAGE).
## Furthermore, the script avoids to repeat time-consuming computations by
## saving intermediate results.
## If you want to recompute all the results from scratch (and use the new
## composition of KEGG pathways), please set the following variable to TRUE
recompute.from.scratch <- FALSE
#recompute.from.scratch <- TRUE

### IMPORTANT
## Input data are not provided with the script.
## CEL files are available at GEO and must be placed (gzip-ed) in the
## subdirectory inputs/cel.
## Protein data are available as supplementary material.
## (More information will be provided when available.)


## Update packages (if needed)
#source("http://bioconductor.org/biocLite.R")
#biocLite(ask=F)
#update.packages(ask=F)



########################################################################
##
## Initialize constants and dependencies
##
options(stringsAsFactors=FALSE)
dir.create("results/", showWarnings=FALSE)
source("func-utils.R")
source("analyze-uarray_piano-utils.R")

dependencies(c("dplyr", "reshape2"))
dependencies("piano", "install.from.bioconductor")
dependencies(c("affy", "plier"), "install.from.bioconductor")
source("piano-modified-funct.R")




## Microarray setup
ctype <- c("normal", "transformed")
cmod <- c("minusFSK", "plusFSK")

contrasts.plus_minus <- data.frame(condition1=c(paste(ctype[1], cmod[2], sep="_"),
                                       paste(ctype[2], cmod[2], sep="_")),
                                   condition2=c(paste(ctype[1], cmod[1], sep="_"),
                                       paste(ctype[2], cmod[1], sep="_")))
contrasts.plus_minus <- base::transform(contrasts.plus_minus,
                                        text=paste(condition1, condition2, sep=" - "))
contrasts <- contrasts.plus_minus

fnames <- sort(apply(expand.grid(c("N", "T"),
                                 c("mF", "pF"),
                                 1:3),
                     1, paste, collapse=""))

uarray.setup <- data.frame(ctype=rep(ctype, each=6),
                           cmod=rep(cmod, each=3))
rownames(uarray.setup) <- fnames
rm(fnames)

## Gene-set enrichment analyses
tests.t <- c("mean", "median", "sum", "maxmean", "gsea", "page")
tests.p <- c("reporter", "wilcoxon", "tailStrength")
signif.method <- c(rep("geneSampling", each=5),
                   rep("nullDist", each=3),
                   "geneSampling")
which.values <- c(rep("t-values", each=length(tests.t)),
                  rep("p-values", each=length(tests.p)))

test.configs <- data.frame(
    gene.set.stat=c(tests.t, tests.p),
    signif.method=signif.method,
    which.values=which.values)
test.configs.t <- test.configs[test.configs$which.values == "t-values",]
test.configs.p <- test.configs[test.configs$which.values == "p-values",]

test.configs

diff.expr.cutoff <- 0.05

cel.data.dir <- "inputs/cel"
abund.data.filename <- "inputs/input-protein-abundance.csv"

print("Setup (PLEASE CHECK!)")
print(uarray.setup)
print(contrasts)



########################################################################
##
## Load uarray expression data (and related data)
##
input.data <- new.env(parent=emptyenv())
compute.with.caching(
    input.data, "input.data", "rdata-uarray-piano_inputdata.RData",
    recompute.from.scratch,
    compute.annotation.and.expression.from.cel,
    cel.data.dir, uarray.setup)
annotation <- input.data$input.data$annotation
array.data <- input.data$input.data$array.data
rm(input.data)

## Perform QC (using the custom function, see file 'piano-modified-funct.R'!!!)
#runQC(array.data, nuseRle=F, rnaDeg=T, hist=T, boxplot=T, pca=T,
#      plot.fnt=print.ggplot2.to.file, prefix="piano_uarray_QC", suffix="")
gc()
close.all.plots()



########################################################################
##
## Perform "classical" differential expression analysis using package
## 'piano' (based on the same approach of 'limma')
##

env.diffexp.res <- new.env(parent=emptyenv())
compute.with.caching(
    env.diffexp.res, "diffexp.res",
    "rdata-uaa-env.diffexp.res.RData",
    recompute.from.scratch,
    diffExp,
    array.data,
    contrasts=contrasts$text,
    significance=diff.expr.cutoff,
    plot=FALSE)
diffexp.res <- env.diffexp.res$diffexp.res
rm(env.diffexp.res)

## Select columns
diffexp.res.proc <- lapply(diffexp.res$resTable,
                           process.piano.diffexp.result)

## Reshape separate tables to a single wide table (for output)
diffexp.res.wide <- reshape.piano.diffexp.results(annotation, diffexp.res.proc)

## Select tcids with minimum median p-value along contrasts
selected.tcids <- (    diffexp.res.wide
                   %>% dplyr::select(egid, tcid, ends_with(" pVal"))
                   %>% melt(1:2)
                   %>% tbl_df()
                   %>% dplyr::select(-variable)
                   %>% dplyr::group_by(egid, tcid)
                   %>% dplyr::summarize(pval=median(value))
                   %>% dplyr::group_by(egid)
                   %>% dplyr::arrange(pval)
                   %>% dplyr::slice(1)
                   %>% dplyr::ungroup()
                   %>% dplyr::select(tcid)
                   )
selected.tcids.annotation <- annotation %>% inner_join(selected.tcids)


diffexp.res.wide <- (    diffexp.res.wide
                     %>% dplyr::inner_join(selected.tcids))

rm(diffexp.res, diffexp.res.proc)


## Output microarray logFC table
write.csv(diffexp.res.wide[order(apply(diffexp.res.wide[, paste(LETTERS[seq_along(contrasts$text)], "pVal")], 1, min)),],
          file="results/out-piano-ua-logfc.csv",
          row.names=FALSE, na="")
gc()

## Reformat expression and diff-expression tables
exp.table <- (    as.data.frame(array.data$dataNorm)
              %>% add_rownames("tcid")
              %>% tbl_df()
              %>% dplyr::mutate(tcid=as.numeric(tcid))
              %>% dplyr::inner_join(selected.tcids)
              %>% dplyr::inner_join(annotation %>% dplyr::select(tcid, egid), .)
              %>% tbl_df()
              %>% dplyr::select(-tcid))
exp.table

exp.table.single <-
    (    uarray.setup
     %>% add_rownames("sample")
     %>% dplyr::mutate(desc=sub("[0-9]$", "", chartr("pm", "+-", sample)))
     %>% dplyr::mutate(desc=sub("-F$", "", desc))
     %>% dplyr::select(sample, desc)
     %>% tbl_df()
     %>% inner_join(melt(exp.table,
                         id.vars="egid", variable.name="sample") %>% tbl_df())
     %>% tbl_df()
     %>% dplyr::select(-sample)
     %>% dplyr::group_by(egid, desc)
     %>% dplyr::summarise_each(funs(mean))
     %>% dcast(egid ~ desc)
     %>% tbl_df()
     )
exp.table.single %>% tbl_df()
diffexp.res.wide %>% tbl_df()




###########################
##
## Protein data
##

prot.input.data <- new.env(parent=emptyenv())
compute.with.caching(
    prot.input.data, "prot.input.data", "rdata-prota-piano_inputdata.RData",
    recompute.from.scratch,
    compute.annotation.and.expression.from.prot.abund,
    abund.data.filename)
prot.abund.headers <- prot.input.data$prot.input.data$headers
prot.abund <- prot.input.data$prot.input.data$prot.abund
rm(prot.input.data)

## Prepare data in the formats used by piano for diff expression analysis
prot.abund.mat <- prot.abund[, rownames(uarray.setup)]
rownames(prot.abund.mat) <- prot.abund$`Master spot no.`
prot.annot <- data.frame(geneName=prot.abund$`Gene Name`,
                         chromosome=prot.abund$chromosome,
                         start=prot.abund$start,
                         egid=prot.abund$egid,
                         msid=prot.abund$`Master spot no.`,
                         row.names=prot.abund$`Master spot no.`)

prot.array.data <- loadMAdata(setup=uarray.setup,
                              dataNorm=prot.abund.mat,
                              annotation=prot.annot)

#runQC(prot.array.data, nuseRle=F, rnaDeg=F, hist=T, boxplot=T, pca=T,
#      plot.fnt=print.ggplot2.to.file, prefix="piano_prota_QC", suffix="")
gc()
close.all.plots()

## Perform the diff expression analysis
env.diffexp.res <- new.env(parent=emptyenv())
compute.with.caching(
    env.diffexp.res, "diffexp.res",
    "rdata-prota-env.diffexp.res.RData",
    recompute.from.scratch,
    diffExp,
    prot.array.data,
    contrasts=contrasts$text,
    significance=diff.expr.cutoff,
    plot=FALSE)
prot.diffexp.res <- env.diffexp.res$diffexp.res
rm(env.diffexp.res)

prot.diffexp.res$resTable[[1]] %>% tbl_df()

## Select columns
prot.diffexp.res.proc <- lapply(prot.diffexp.res$resTable,
                                process.piano.diffexp.result)

## Reshape separate tables to a single wide table (for output)
prot.annot.as.tcid <- prot.annot %>% dplyr::rename(tcid=msid)
prot.diffexp.res.wide <- reshape.piano.diffexp.results(prot.annot.as.tcid, prot.diffexp.res.proc)


## Select tcids with minimum median p-value along contrasts
prot.selected.tcids <- (    prot.diffexp.res.wide
                        %>% dplyr::select(egid, tcid, ends_with(" pVal"))
                        %>% melt(1:2)
                        %>% tbl_df()
                        %>% dplyr::select(-variable)
                        %>% dplyr::group_by(egid, tcid)
                        %>% dplyr::summarize(pval=median(value))
                        %>% dplyr::group_by(egid)
                        %>% dplyr::arrange(pval)
                        %>% dplyr::slice(1)
                        %>% dplyr::ungroup()
                        %>% dplyr::select(tcid)
                        )
prot.selected.tcids.annotation <- prot.annot.as.tcid %>% inner_join(prot.selected.tcids)

prot.diffexp.res.wide <- (    prot.diffexp.res.wide
                          %>% dplyr::inner_join(prot.selected.tcids))

rm(prot.diffexp.res, prot.diffexp.res.proc)

## Output protein logFC table
write.csv(prot.diffexp.res.wide[order(apply(prot.diffexp.res.wide[, paste(LETTERS[seq_along(contrasts$text)], "pVal")], 1, min)),],
          file="results/out-piano-prot-logfc.csv",
          row.names=FALSE, na="")

## Reformat expression and diff-expression tables
prot.exp.table <- (    as.data.frame(prot.array.data$dataNorm)
                   %>% add_rownames("tcid")
                   %>% tbl_df()
                   %>% dplyr::mutate(tcid=as.numeric(tcid))
                   %>% dplyr::inner_join(prot.selected.tcids)
                   %>% dplyr::inner_join(prot.annot.as.tcid %>% dplyr::select(tcid, egid), .)
                   %>% tbl_df()
                   %>% dplyr::select(-tcid))
prot.exp.table

mean.no.na <- function(...) mean(..., na.rm=TRUE)

prot.exp.table.single <-
    (    uarray.setup
     %>% add_rownames("sample")
     %>% dplyr::mutate(desc=sub("[0-9]$", "", chartr("pm", "+-", sample)))
     %>% dplyr::mutate(desc=sub("-F$", "", desc))
     %>% dplyr::select(sample, desc)
     %>% tbl_df()
     %>% inner_join(melt(prot.exp.table,
                         id.vars="egid", variable.name="sample") %>% tbl_df())
     %>% tbl_df()
     %>% dplyr::select(-sample)
     %>% dplyr::group_by(egid, desc)
     %>% dplyr::summarise_each(funs(mean.no.na))
     %>% dcast(egid ~ desc)
     %>% tbl_df()
     )
prot.exp.table.single %>% tbl_df()
prot.diffexp.res.wide %>% tbl_df()



########################################################################
##
##
##
##   G E N E   S E T   A N A L Y S I S
##   (using package 'piano')
##
##


########################################################################
##
## Load gene sets from KEGG (using those stored in package 'gage')
##
piano.gsets <-
    load.gene.sets.from.kegg.using.gage(species="mmu",
                                        recompute.from.scratch=recompute.from.scratch)
print(piano.gsets)


env.kegg.sets <- new.env(parent=emptyenv())
compute.with.caching(
    env.kegg.sets, "kg.gsets",
    s.join("rdata-kegg-gsets-mmu.RData"),
    recompute.from.scratch,
    kegg.gsets,
    species="mmu",
    id.type="entrez")
kg.gsets <- env.kegg.sets$kg.gsets

## Annotate pathways with their type
gset.classif <- data.frame(Pathway=names(kg.gsets$kg.sets),
                           Type="")
gset.classif$Type[kg.gsets$sig.idx] <- "Signaling"
gset.classif$Type[kg.gsets$met.idx] <- "Metabolic"
gset.classif$Type[kg.gsets$dise.idx] <- "Disease"
gset.classif <- gset.classif %>% tbl_df()

## Map pathways with their genes
gset.genes <- (    piano.gsets$gsc
               %>% unlist()
               %>% setNames(NULL)
               %>% unique())

gset.gene.mapping <-
    lapply(names(piano.gsets$gsc), function(x) data_frame(Pathway=x, `Gene ID`=piano.gsets$gsc[[x]]))
gset.gene.mapping <-
    Reduce(function(...) merge(..., all=TRUE), gset.gene.mapping)
gset.gene.mapping <- gset.gene.mapping %>% tbl_df()
gset.gene.mapping <- gset.classif %>% dplyr::inner_join(gset.gene.mapping) %>% tbl_df()
gset.gene.mapping$`Gene ID` <- as.numeric(gset.gene.mapping$`Gene ID`)

tmp.exp.table.single <- exp.table.single
colnames(tmp.exp.table.single)[-1] <- paste("Transcript", colnames(tmp.exp.table.single)[-1])
tmp.prot.exp.table.single <- prot.exp.table.single
colnames(tmp.prot.exp.table.single)[-1] <- paste("Protein", colnames(tmp.prot.exp.table.single)[-1])

gset.gene.mapping$`Gene ID` <- as.integer(gset.gene.mapping$`Gene ID`)
prot.selected.tcids.annotation$`egid` <- as.integer(prot.selected.tcids.annotation$`egid`)
selected.tcids.annotation$`egid` <- as.integer(selected.tcids.annotation$`egid`)
tmp.prot.exp.table.single$`egid` <- as.integer(tmp.prot.exp.table.single$`egid`)

gset.gene.mapping <- ( (    gset.gene.mapping
                      %>% inner_join(selected.tcids.annotation
                                     %>% dplyr::select(`Gene ID`=egid, geneName, `Transcript ID`=tcid))
                        %>% inner_join(tmp.exp.table.single, by=c("Gene ID"="egid")))
                      %>% full_join(    gset.gene.mapping
                                    %>% inner_join(prot.selected.tcids.annotation
                                                   %>% dplyr::select(`Gene ID`=egid, geneName, `Protein Spot ID`=tcid))
                                    %>% inner_join(tmp.prot.exp.table.single, by=c("Gene ID"="egid"))
                                    )
                      %>% arrange(Pathway, `Gene ID`) %>% tbl_df()
                      )
write.csv(gset.gene.mapping, "results/geneset-mapping.csv", row.names=FALSE, na="")

gene.to.pathways <- gset.gene.mapping %>% dplyr::select(egid=`Gene ID`, Pathway) %>% group_by(egid) %>% summarise(Pathways=paste(Pathway, collapse="; "))

#################################
##
## Plot expression levels
##
min.expr.value <- min(exp.table.single[,-1], na.rm=TRUE) ## strict
max.expr.value <- max(exp.table.single[,-1], na.rm=TRUE) ## strict

## ...of all genes (requires a large amount of memory!!)
plot.expression.levels(exp.table.single, diffexp.res.wide,
                       selected.tcids.annotation,
                       gene.to.pathways,
                       prefix.fname="results/ua-genes_all-",
                       contrasts=LETTERS[seq_along(contrasts$text)],
                       min.value=min.expr.value,
                       max.value=max.expr.value)


## ...of genes with abs(logFC)>= 1 (for some contrast)
plot.expression.levels(exp.table.single %>%
                           inner_join(    diffexp.res.wide
                                      %>% dplyr::select(egid, ends_with("logFC"))
                                      %>% melt(id.vars="egid")
                                      %>% tbl_df()
                                      %>% dplyr::filter(abs(value) >= 1)
                                      %>% dplyr::select(egid)
                                      %>% dplyr::group_by(egid)
                                      %>% dplyr::slice(1)
                                      %>% dplyr::ungroup()
                                      ),
                       diffexp.res.wide,
                       selected.tcids.annotation,
                       gene.to.pathways,
                       prefix.fname="results/ua-genes_abs_logfc_gt_1-",
                       contrasts=LETTERS[seq_along(contrasts$text)],
                       min.value=min.expr.value,
                       max.value=max.expr.value)

for (i in seq_along(contrasts$text)) {
    write.csv((    diffexp.res.wide[abs(diffexp.res.wide[,paste(LETTERS[i], "logFC")]) >= 1,]
               %>% dplyr::select(geneName)
               %>% tbl_df()
               ), file=s.join("results/genes-", LETTERS[i], "-abs_logfc_gt_1.csv"),
              row.names=FALSE, na="")
}


## Protein data
prot.min.expr.value <- min(prot.exp.table.single[,-1], na.rm=TRUE)   ## strict
prot.max.expr.value <- max(prot.exp.table.single[,-1], na.rm=TRUE)   ## strict

## ...of all proteins
plot.expression.levels(prot.exp.table.single, prot.diffexp.res.wide,
                       prot.selected.tcids.annotation,
                       gene.to.pathways,
                       prefix.fname="results/prot-genes_all-",
                       contrasts=LETTERS[seq_along(contrasts$text)],
                       min.value=prot.min.expr.value,
                       max.value=prot.max.expr.value)
##
## end plot expression levels
##
#################################




########################################################################
##
## Perform the gene-set analyses
##
run.gene.set.analyses.all.contrasts <- function(
    contrasts,
    test.configs,
    piano.gsets,
    .data,
    sep=" ", egid.col="egid",
    t.values.col="tVal", p.values.col="pVal", log.fc.col="logFC",
    ...) {
    gsa.res <- list()
    for (i in seq_along(contrasts$text)) {
        gsa.res[[contrasts$text[i]]] <-
            run.gene.set.analyses(test.configs, piano.gsets,
                                  .data, prefix.col=LETTERS[i], ...)
    }
    gsa.res
}

##  analyses on t-values on ua
gsa.res.uaa.t <- new.env(parent=emptyenv())
compute.with.caching(
    gsa.res.uaa.t, "res", "rdata-piano-gsa-results-uaa-t.RData",
    recompute.from.scratch,
    run.gene.set.analyses.all.contrasts,
    contrasts,
    test.configs.t,
    piano.gsets,
    diffexp.res.wide,
    sep=" ",
    egid.col="egid",
    t.values.col="tVal",
    p.values.col="pVal",
    log.fc.col="logFC")
gsa.res.uaa.t <- gsa.res.uaa.t$res

##  analyses on p-values on ua
gsa.res.uaa.p <- new.env(parent=emptyenv())
compute.with.caching(
    gsa.res.uaa.p, "res", "rdata-piano-gsa-results-uaa-p.RData",
    recompute.from.scratch,
    run.gene.set.analyses.all.contrasts,
    contrasts,
    test.configs.p,
    piano.gsets,
    diffexp.res.wide,
    sep=" ",
    egid.col="egid",
    t.values.col="tVal",
    p.values.col="pVal",
    log.fc.col="logFC")
gsa.res.uaa.p <- gsa.res.uaa.p$res

##  analyses on t-values on prot
gsa.res.prot.t <- new.env(parent=emptyenv())
compute.with.caching(
    gsa.res.prot.t, "res", "rdata-piano-gsa-results-prot-t.RData",
    recompute.from.scratch,
    run.gene.set.analyses.all.contrasts,
    contrasts,
    test.configs.t,
    piano.gsets,
    prot.diffexp.res.wide,
    sep=" ",
    egid.col="egid",
    t.values.col="tVal",
    p.values.col="pVal",
    log.fc.col="logFC")
gsa.res.prot.t <- gsa.res.prot.t$res

##  analyses on p-values on prot
gsa.res.prot.p <- new.env(parent=emptyenv())
compute.with.caching(
    gsa.res.prot.p, "res", "rdata-piano-gsa-results-prot-p.RData",
    recompute.from.scratch,
    run.gene.set.analyses.all.contrasts,
    contrasts,
    test.configs.p,
    piano.gsets,
    prot.diffexp.res.wide,
    sep=" ",
    egid.col="egid",
    t.values.col="tVal",
    p.values.col="pVal",
    log.fc.col="logFC")
gsa.res.prot.p <- gsa.res.prot.p$res

warnings()

##
## end of basic gene-set analyses
##
########################################################################

## Merge results
merge.results <- function(res1, res2) {
    stopifnot(names(res1)==names(res2))
    dependencies("magrittr")
    names(res1) %>%
        lapply( . %>% { c(res1[[.]], res2[[.]]) } ) %>%
            setNames(names(res1))
}

dependencies("magrittr")
gsa.res.uaa <- merge.results(gsa.res.uaa.t, gsa.res.uaa.p) %>%
    lapply(. %>% {setNames(., paste("UA", names(.)))})
gsa.res.prot <- merge.results(gsa.res.prot.t, gsa.res.prot.p) %>%
    lapply(. %>% { setNames(., paste("PROT", names(.))) })
gsa.res <- merge.results(gsa.res.uaa, gsa.res.prot)

stopifnot(all(contrasts$text %in% names(gsa.res.uaa)))
stopifnot(all(contrasts$text %in% names(gsa.res.prot)))
stopifnot(all(contrasts$text %in% names(gsa.res)))


## Function for reformatting PIANO results
enlarge.all.results <- function(gsa.res, original.gsets) {
    dependencies("magrittr")

    .enlarge.results <- function(this.res, common.gsets) {
        fields.data <- data.frame(
            field=c("nGenesTot", "nGenesUp", "nGenesDn",
                "statDistinctDir", "statDistinctDirUp", "statDistinctDirDn", "statNonDirectional",
                "statMixedDirUp", "statMixedDirDn",
                "pDistinctDirUp", "pDistinctDirDn", "pNonDirectional",
                "pMixedDirUp", "pMixedDirDn",
                "pAdjDistinctDirUp", "pAdjDistinctDirDn", "pAdjNonDirectional",
                "pAdjMixedDirUp", "pAdjMixedDirDn"),
            default.value=rep(c(0, NA, NA), times=c(3, 6, 10)))
        
        .enlarge.matrix <- function(orig.mat, new.nrow, match.rows,
                                   default.value=NA) {
            new.mat <- matrix(data=default.value,
                              nrow=new.nrow, ncol=ncol(orig.mat),
                              dimnames=dimnames(orig.mat))
            new.mat[match.rows,] <- orig.mat
            new.mat
        }

        match.gsc <- match(names(this.res$gsc), names(common.gsets))
        new.res <- this.res
        for (i in 1:nrow(fields.data)) {
            new.res[[fields.data$field[i]]] <-
                .enlarge.matrix(this.res[[fields.data$field[i]]],
                                length(common.gsets), match.gsc,
                                default.value=fields.data$default.value[i])
        }
        new.res$gsc <- common.gsets
        new.res
    }

    common.gsets <- (    names(unlist(gsa.res, recursive=FALSE))
                     %>% lapply(. %>% { names(unlist(gsa.res, recursive=FALSE)[[.]]$gsc) })
                     %>% unlist
                     %>% unique
                     %>% { original.gsets$gsc[names(original.gsets$gsc) %in% .] }
                     )

    return (    gsa.res
            %>% lapply(. %>% { lapply(., .enlarge.results,
                                      common.gsets=common.gsets) } ) )
}

new.gsa.res <- enlarge.all.results(gsa.res, piano.gsets)

## Funtion for computing the consensus rankings/heatmaps for each contrast
compute.all.consensus.heatmaps <- function(gsa.res, contrasts, cutoff=10) {
    stopifnot(all(contrasts$text %in% names(gsa.res)))
    ch <- list()
    for (i in seq_along(contrasts$text)) {
        contrast <- contrasts$text[i]
        ch[[contrast]] <-
            consensusHeatmap(
                gsa.res[[contrast]],
                cutoff=cutoff,
                method="Borda",
                adjusted=TRUE,
                columnnames="abbr"
            )
    }
    invisible(ch)
}

## Compute the consensus rankings (cutoff > number of pathways, thus all pathways are included in the rankings)
ch.uaa <- compute.all.consensus.heatmaps(gsa.res.uaa, contrasts, cutoff=1000)
ch.prot <- compute.all.consensus.heatmaps(gsa.res.prot, contrasts, cutoff=1000)
ch.all <- compute.all.consensus.heatmaps(new.gsa.res, contrasts, cutoff=1000)

## Reformat the rankings (from long to wide)
res.df <- NULL
for (i in seq_along(contrasts$text)) {
    contrast <- contrasts$text[i]
    for (matn in c("rankMat", "pMat")) {
        rm.df <- as.data.frame(ch.uaa[[contrast]][[matn]])
        colnames(rm.df) <- c("Dist(dn)","Mix(dn)","Nondir","Mix(up)","Dist(up)")
        colnames(rm.df) <- paste(LETTERS[i], matn, "UA", colnames(rm.df))
        rm.df$`Gene Set` <- rownames(rm.df)
        rm.df <- rm.df[,c("Gene Set", names(rm.df)[names(rm.df) != "Gene Set"])]
        res.df <- if (is.null(res.df)) rm.df else full_join(res.df, rm.df)

        rm.df <- as.data.frame(ch.prot[[contrast]][[matn]])
        colnames(rm.df) <- c("Dist(dn)","Mix(dn)","Nondir","Mix(up)","Dist(up)")
        colnames(rm.df) <- paste(LETTERS[i], matn, "PR", colnames(rm.df))
        rm.df$`Gene Set` <- rownames(rm.df)
        res.df <- full_join(res.df, rm.df)

        rm.df <- as.data.frame(ch.all[[contrast]][[matn]])
        colnames(rm.df) <- c("Dist(dn)","Mix(dn)","Nondir","Mix(up)","Dist(up)")
        colnames(rm.df) <- paste(LETTERS[i], matn, "ALL", colnames(rm.df))
        rm.df$`Gene Set` <- rownames(rm.df)
        res.df <- full_join(res.df, rm.df)
    }
    stopifnot(all(ch.uaa[[contrast]]$nGenesMat[,1] == ch.uaa[[contrast]]$nGenesMat[,3]))
    stopifnot(all(ch.uaa[[contrast]]$nGenesMat[,1] == ch.uaa[[contrast]]$nGenesMat[,5]))
    rm.df <- as.data.frame(ch.uaa[[contrast]]$nGenesMat[,c(1,2,4)])
    colnames(rm.df) <- c("n", "nDn", "nUp")
    colnames(rm.df) <- paste(LETTERS[i], "genes", "UA", colnames(rm.df))
    rm.df$`Gene Set` <- rownames(rm.df)
    res.df <- if (is.null(res.df)) rm.df else full_join(res.df, rm.df)
    
    stopifnot(all(ch.prot[[contrast]]$nGenesMat[,1] == ch.prot[[contrast]]$nGenesMat[,3]))
    stopifnot(all(ch.prot[[contrast]]$nGenesMat[,1] == ch.prot[[contrast]]$nGenesMat[,5]))
    rm.df <- as.data.frame(ch.prot[[contrast]]$nGenesMat[,c(1,2,4)])
    colnames(rm.df) <- c("n", "nDn", "nUp")
    colnames(rm.df) <- paste(LETTERS[i], "genes", "PR", colnames(rm.df))
    rm.df$`Gene Set` <- rownames(rm.df)
    res.df <- full_join(res.df, rm.df)

    stopifnot(all(ch.all[[contrast]]$nGenesMat[,1] == ch.all[[contrast]]$nGenesMat[,3]))
    stopifnot(all(ch.all[[contrast]]$nGenesMat[,1] == ch.all[[contrast]]$nGenesMat[,5]))
    rm.df <- as.data.frame(ch.all[[contrast]]$nGenesMat[,c(1,2,4)])
    colnames(rm.df) <- c("n", "nDn", "nUp")
    colnames(rm.df) <- paste(LETTERS[i], "genes", "ALL", colnames(rm.df))
    rm.df$`Gene Set` <- rownames(rm.df)
    res.df <- full_join(res.df, rm.df)
}
res.df.small <- res.df %>% dplyr::select(-contains(" Mix"), -contains(" pMat "), -contains(" genes ALL ")) %>% tbl_df()
colnames(res.df.small) <- sub("rankMat", "rank", colnames(res.df.small))
write.csv(res.df.small, "results/piano-consensus-geneset_table.csv", row.names=FALSE, na="")

## Function for printing the consensus heatmaps
plot.rich.consensus.heatmaps <- function(ch.res, contrasts, gset.classif,
                                         prefix.fname="piano-consensus-",
                                         print.ranks=TRUE, max.rank=NULL,
                                         add.signif=TRUE, cutoff=10) {
    dependencies("dplyr")
    dependencies("reshape2")
    dependencies("RColorBrewer")
    dependencies("ggplot2")
    dependencies("grid")
    stopifnot(all(contrasts$text %in% names(ch.res)))

    .join.or.assign <- function(t1, t2) {
        if (is.null(t1))
            t2
        else
            inner_join(t1, t2)
    }

    tot.rank.mat <- NULL
    for (i in seq_along(contrasts$text)) {
        rank.mat <- ch.res[[i]]$rankMat
        rank.mat <- rank.mat[,c(1,3,5)]
        colnames(rank.mat) <- paste(LETTERS[i], c("Down","Mixed","Up"))
        rank.mat <- rank.mat %>% as.data.frame() %>% add_rownames("Pathway") %>% tbl_df()
        tot.rank.mat <- .join.or.assign(tot.rank.mat, rank.mat)
    }

    rank.mat <- tot.rank.mat

    rank.mask <- apply(rank.mat %>% dplyr::select(-Pathway) %>% as.matrix(),
                       1, function(x) any(x<=cutoff))
    rank.mat <- rank.mat[rank.mask,]

    rank.mat$Pathway <- sub("^mmu[0-9]+ ", "", rank.mat$Pathway, perl=TRUE)
    rank.coef <- rep(c(1, 1, -1), nrow(contrasts))
    rank.coef2 <- rep(c(1.5, 0.05, 0.5), nrow(contrasts))
    rank.coef2 <- c(1.5, 0.05, 0.5, 0.5, 0.02, 0.2)
    row.order <- order(rowSums(t(t(apply(-apply(t(t(rank.mat %>% dplyr::select(-Pathway) %>% as.matrix()) * rank.coef), 2, rank, ties.method="min"), 2,
                                         rank, ties.method="min")-1)*rank.coef2)))
    rank.mat.long <- melt(rank.mat, id.vars=c("Pathway"))
    rank.mat.long$Pathway <- factor(rank.mat.long$Pathway, levels=rank.mat$Pathway[row.order], ordered=TRUE)
    rank.mat.long$Contrast <- substr(rank.mat.long$variable, start=1, stop=1)
    rank.mat.long$Contrast <- factor(rank.mat.long$Contrast, levels=LETTERS[1:length(contrasts$text)], ordered=TRUE)
    levels(rank.mat.long$Contrast) <- contrasts$abbrev
    rank.mat.long$variable <- substr(rank.mat.long$variable, start=3, stop=10000)

    real.max.rank <- max(rank.mat.long$value, na.rm=TRUE)
    if (is.null(max.rank))
        max.rank <- length(unique(as.character(gset.classif$Pathway)))
    real.max.rank <- max.rank

    pl.fname <- s.join(prefix.fname, ".pdf")


    fill.breaks <- c(1) #, 10)
    if (real.max.rank >=50)
        fill.breaks <- c(fill.breaks, seq(from=50, to=real.max.rank, by=50))
    fill.breaks <- c(fill.breaks, real.max.rank)
    p <- (  ggplot(rank.mat.long, aes(x=variable, y=Pathway, fill=value))
          + theme_bw()
          + facet_grid(.~Contrast)
          + geom_tile()
          + xlab("")
          + ylab("Pathway enrichment")
          + theme(
              panel.margin = unit(0, "pt"),
              panel.grid = element_blank(),
              panel.border = element_blank(),
              axis.ticks = element_blank(),
              axis.ticks.margin = unit(-3, "pt"),
              axis.text.x = element_text(hjust=1, vjust=0.5, angle=90, size=8),
              axis.text.y = element_text(hjust=0, size=rel(0.75)),
              axis.title.y = element_text(hjust=0.5, vjust=0.5, angle=90, size=10, face="bold"),
              strip.background = element_blank(),
              strip.text = element_text(face="bold", size=9),
              legend.position = c(1.4, 0.865),
              complete=FALSE
          )
          + guides(fill = guide_colourbar(
                       title.position = "left",
                       title.hjust = 0.5,
                       title.vjust = 0.5,
                       title.theme=element_text(angle=90, size=8, face="bold"),
                       label.hjust=0,
                       label.theme=element_text(size=7, angle=0),
                       barheight=8, barwidth=0.6,
                       ticks=FALSE, reverse=TRUE, nbin=max.rank+1
                   )
          )
          + scale_fill_gradientn(
              name="Pathway ranking",
              colours = c("#67000d", "#ef3b2c", "#fc9272", "#fee0d2", "#ffffff"), ##fcae91
              values = c(1, 10, 11, 60, max.rank),
              breaks = fill.breaks,
              labels = fill.breaks,
              limits = c(1, max.rank),
              rescaler = function(x, ...) x, oob = identity
          )
          + scale_x_discrete(expand=c(0,0))
          + scale_y_discrete(expand=c(0,0))
          + coord_fixed(ratio=1)
          )
    if (print.ranks) p <- p+geom_text(aes(label=value), size=2)

    ggsave(filename=pl.fname, p)
    invisible(NULL)
}

contrasts$abbrev <- c("NF/N", "TF/T")
heatmap.max.rank <- 200

plot.rich.consensus.heatmaps(ch.uaa, contrasts, gset.classif,
                             prefix.fname="results/piano-consensus-1_uaa",
                             print.ranks=TRUE, add.signif=FALSE, cutoff=10, max.rank=heatmap.max.rank)
plot.rich.consensus.heatmaps(ch.prot, contrasts, gset.classif,
                             prefix.fname="results/piano-consensus-2_prot",
                             print.ranks=TRUE, add.signif=FALSE, cutoff=5, max.rank=heatmap.max.rank)
plot.rich.consensus.heatmaps(ch.all, contrasts, gset.classif,
                             prefix.fname="results/piano-consensus-3_prot_uaa",
                             print.ranks=TRUE, add.signif=FALSE, cutoff=10, max.rank=heatmap.max.rank)

message("Completed")
