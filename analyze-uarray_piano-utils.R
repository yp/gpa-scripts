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
## This script provides some convenience functions to use for performing
## the gene/protein analyses
source("func-utils.R")


## Convenience function called by the modified function "runQC" of piano package
## for saving the plots to file
print.ggplot2.to.file <- function(p, plot.name="plot", extension="pdf",
                                  prefix="", suffix="", word.sep="-", ...) {
    prefix.t <- if (nchar(prefix)>0) s.join(prefix, word.sep) else prefix;
    suffix.t <- if (nchar(suffix)>0) s.join(word.sep, suffix) else suffix;
    ggsave(filename=s.join(prefix.t, plot.name, suffix.t, ".", extension), plot=p, ...);
    invisible(NULL);
}




## Read uarray raw expression level from CEL files and perform normalization using the
## specialized package 'oligo', load the annotation from biomart, and then pass all these
## data to package 'piano'
compute.annotation.and.expression.from.cel <- function(data.dir,
                                                       uarray.setup) {
    dependencies(c("oligo", "piano"), "install.from.bioconductor")
    fnames <- rownames(uarray.setup)
    affyRaw <- read.celfiles(s.join(data.dir, "/", fnames, ".CEL.gz"),
                             sampleNames=fnames)
    eset <- rma(affyRaw,
                target = "core",
                normalize = TRUE,
                background = TRUE)
    basic.expr.levels <- exprs(eset)

    ## Import dependencies
    dependencies("dplyr")
    annotation <- enrich.from.biomart(data.frame(tcid=sort(as.numeric(rownames(basic.expr.levels)))),
                                      key.field="tcid", biomart.key.field="affy_mogene_1_0_st_v1",
                                      attributes=c("entrezgene", "mgi_symbol", "chromosome_name", "start_position", "strand"),
                                      attributes.names=c("egid", "geneName", "chromosome", "start_position", "strand"))
    annotation <-
        data.frame(
            annotation
            %>% filter(!is.na(egid))
            %>% mutate(start=strand*start_position)
            %>% dplyr::select(geneName, chromosome, start, egid, tcid)
            %>% arrange(egid)
            %>% group_by(tcid)
            %>% dplyr::slice(1)
        )
    rownames(annotation) <- annotation$tcid
    detach("package:pd.mogene.1.0.st.v1", unload=TRUE)
    detach("package:oligo", unload=TRUE)
    myArrayData <- loadMAdata(setup=uarray.setup, datadir=data.dir,
                              dataRaw=basic.expr.levels,
                              dataNorm=basic.expr.levels, annotation=annotation)

    myArrayData$dataRaw <- ReadAffy(celfile.path=data.dir)
    colnames(exprs(myArrayData$dataRaw)) <- gsub("\\.CEL(\\.gz)?", "",
                                                 colnames(exprs(myArrayData$dataRaw)),
                                                 ignore.case = TRUE)
    exprs(myArrayData$dataRaw) <- exprs(myArrayData$dataRaw)[,fnames]
    stopifnot(colnames(exprs(myArrayData$dataRaw)) %in% fnames)
    stopifnot(fnames %in% colnames(exprs(myArrayData$dataRaw)))
    list(annotation=annotation,
         array.data=myArrayData)
}




## Post-process the table resulting from differential expression analysis (on uarray data)
## performed with function 'piano::diffExp'
process.piano.diffexp.result <- function(res.table) {
    dependencies("dplyr")
    stat.values <- data_frame(tcid=as.numeric(res.table$ProbesetID),
                              pVal=res.table$adj.P.Val,
                              tVal=res.table$t,
                              logFC=res.table$logFC)
    stat.values
}




## Reshape the tables resulting from differential expression analysis (on uarray data)
## performed with function 'piano::diffExp' to a 'wide' single table
reshape.piano.diffexp.results <- function(annotation, diffexp.res.proc) {
    diffexp.res.wide <- data_frame(tcid=diffexp.res.proc[[1]]$tcid)
    for (i in seq_along(diffexp.res.proc)) {
        tmp <- diffexp.res.proc[[i]]
        colnames(tmp)[-1] <- paste(LETTERS[i], colnames(tmp)[-1])
        diffexp.res.wide <- diffexp.res.wide %>% inner_join(tmp) %>% tbl_df()
        rm(tmp)
    }
    annotation %>% inner_join(diffexp.res.wide)
}




## Reshape the tables resulting from differential expression analysis (on uarray data)
## performed with function 'piano::diffExp' to a 'wide' single table
annotate.piano.diffexp.wide.results <- function(diffexp.res.wide,
                                                diff.expr.cutoff,
                                                contrasts) {
    ## Generate all possible combinations of contrasts
    contrast.combs <- list()
    for (i in seq_along(contrasts$text)) {
        contrast.combs <- c(contrast.combs,
                            lapply(apply(gtools::combinations(length(contrasts$text), i, seq_along(contrasts$text)),
                                         1, list), function(x) x[[1]]))
    }
    names(contrast.combs) <- sapply(contrast.combs, function(x) paste(LETTERS[x], collapse=""))

    ## Mark as significant (with a '*') genes that have p-value less than diff.expr.cutoff
    ## for, at least, one contrast
    tmp <- apply(diffexp.res.wide[, paste(LETTERS[seq_along(contrasts$text)], "pVal")], 1, min) <= diff.expr.cutoff
    diffexp.res.wide$`Signif.` <- FALSE
    diffexp.res.wide$`Signif.`[tmp] <- TRUE
    rm(tmp)

    ## Select columns with p-values and elements with significant p-values
    p.vals <- diffexp.res.wide[,match(paste(LETTERS[seq_along(contrasts$text)], "pVal"), colnames(diffexp.res.wide))]
    p.vals <- p.vals <= diff.expr.cutoff

    ## Mark significant genes for a specific combination (-> 'yes.set')
    compute.uniquely.elements <- function(x, yes.set, no.set=-yes.set) {
        not.any <- function(x) !any(x)
        apply(x[, yes.set, drop=FALSE], 1, all) & apply(x[, no.set, drop=FALSE], 1, not.any)
    }
    tmp <- lapply(contrast.combs, compute.uniquely.elements, x=p.vals)
    names(tmp) <- paste("Uniquely in", names(tmp))
    diffexp.res.wide <- cbind(diffexp.res.wide, data.frame(tmp, check.names=FALSE))
    diffexp.res.wide
}




## Prepare gene sets in the format used by package 'piano' starting from the gene sets
## provided by package 'gage' (that is, the ones from KEGG)
load.gene.sets.from.kegg.using.gage <- function(species="mmu",
                                                recompute.from.scratch=FALSE) {
    dependencies(c("gage", "piano"), "install.from.bioconductor")
    env.kegg.sets <- new.env(parent=emptyenv())
    compute.with.caching(
        env.kegg.sets, "kg.gsets",
        s.join("rdata-kegg-gsets-", species, ".RData"),
        recompute.from.scratch,
        kegg.gsets,
        species=species,
        id.type="entrez")
    kg.gsets <- env.kegg.sets$kg.gsets
    #### UNCOMMENT ONLY ONE OF THE FOLLOWING !!!
    ### keep all gene sets
    ## gsets <- kg.gsets$kg.sets
    ### keep only SIGNALLING and METHABOLISM
    gsets <- kg.gsets$kg.sets[kg.gsets$sigmet.idx]

    piano.gsets.df <- do.call("rbind", lapply(names(gsets),
                                              function (x) { data.frame(genes=gsets[[x]], set=x) }))
    piano.gsets <- loadGSC(as.matrix(piano.gsets.df))

    detach("package:gage", unload=TRUE)
    piano.gsets
}




compute.annotation.and.expression.from.prot.abund <- function(data.filename) {
    dependencies("dplyr")
    dependencies("rentrez")

    prot.abund <- read.csv(data.filename, check.names=FALSE)
    stopifnot(!any(duplicated(prot.abund$`Master spot no.`)))

    enrich.prot.abund <- enrich.from.biomart(prot.abund,
                                             key.field="Gene Name", biomart.key.field="mgi_symbol",
                                             attributes=c("entrezgene", "chromosome_name", "start_position", "strand"),
                                             attributes.names=c("egid", "chromosome", "start", "strand"))

    gene.names <- enrich.prot.abund$`Gene Name`[is.na(enrich.prot.abund$egid)]
    res <- data.frame(egid=c(), `Gene Name`=c(), check.names=FALSE)
    for (gn in gene.names) {
        ressearch <- entrez_search(db="gene", s.join("\"Mouse\"[Organism] AND \"", gn, "\"[Gene]"))
        this.res <- data.frame(egid=ressearch$ids, `Gene Name`=gn, check.names=FALSE)
        res <- rbind(res, this.res)
    }
    prot.na <- (    enrich.prot.abund[is.na(enrich.prot.abund$egid), ]
                %>% dplyr::select(-`egid`) ) %>% inner_join(res)
    prot.na <- prot.na[,colnames(enrich.prot.abund)]
    enrich.prot.abund <- rbind(enrich.prot.abund[!is.na(enrich.prot.abund$egid), ], prot.na)
    stopifnot(!any(is.na(enrich.prot.abund$egid)))

    dup.master.spot.no <- duplicated.values(enrich.prot.abund$`Master spot no.`)
    enrich.prot.abund <- subset(enrich.prot.abund,
                                `Master spot no.` %ni% dup.master.spot.no |
                                (`Master spot no.` %in% dup.master.spot.no & `chromosome` %in% c(1:22, "X", "Y")))
    dup.master.spot.no <- duplicated.values(enrich.prot.abund$`Master spot no.`)
    enrich.prot.abund.dup <- subset(enrich.prot.abund, `Master spot no.` %in% dup.master.spot.no)
    enrich.prot.abund.ndup <- subset(enrich.prot.abund, `Master spot no.` %ni% dup.master.spot.no)
    message("No. of ambiguous gene names: ", length(unique(enrich.prot.abund.dup$`Gene Name`)))
    message("--> ", paste(unique(enrich.prot.abund.dup$`Gene Name`), collapse=", "))


    message("Finding (hopefully) correct Entrez gene IDs for these genes by querying the Entrez web server..")
    entrez.res <- entrez_summary(db="gene", id=unique(enrich.prot.abund.dup$egid))
    new.entrez.res <- NULL
    for (egid in names(entrez.res)) {
        if (("name" %in% names(entrez.res[[egid]])) && ("status" %in% names(entrez.res[[egid]]))) {
            this.entrez.res <- data.frame(`egid`=egid,
                                          `Entrez Gene Name`=entrez.res[[egid]][["name"]],
                                          `Status`=entrez.res[[egid]][["status"]],
                             check.names=FALSE)
            new.entrez.res <- if (is.null(new.entrez.res)) this.entrez.res else rbind(new.entrez.res, this.entrez.res)
        }
    }
    entrez.res <- new.entrez.res
    enrich.prot.abund.dup.filt <- inner_join(enrich.prot.abund.dup, entrez.res)
    enrich.prot.abund.dup.filt <- subset(enrich.prot.abund.dup.filt,
                                         `Gene Name`==`Entrez Gene Name` & Status==0)

    enrich.prot.abund <- rbind(enrich.prot.abund.ndup,
                               subset(enrich.prot.abund.dup.filt,
                                      select=-c(`Entrez Gene Name`, `Status`)))

    message("Are there duplicated master spot numbers? ", any(duplicated(enrich.prot.abund$`Master spot no.`)))

    stopifnot(all(prot.abund$`Master spot no.` %in% enrich.prot.abund$`Master spot no.`))
    stopifnot(!any(duplicated(enrich.prot.abund$`Master spot no.`)))

    enrich.prot.abund$start <- enrich.prot.abund$strand * enrich.prot.abund$start
    enrich.prot.abund$strand <- NULL

    headers <- c("Master spot no.", "egid", "Gene Name",  "Protein Name", "chromosome", "start")
    list(headers=headers,
         prot.abund=enrich.prot.abund[, c(headers, colnames(enrich.prot.abund)[colnames(enrich.prot.abund) %ni% headers])])
}




## Perform various gene set anayses as specified in the test configs
run.gene.set.analyses <- function(test.configs,
                                  piano.gsets,
                                  diffexp.df,
                                  egid.col="egid",
                                  prefix.col="",
                                  sep=" ",
                                  t.values.col="tVal",
                                  p.values.col="pVal",
                                  log.fc.col="logFC",
                                  ncpus=NULL, nPerm=NULL,
                                  ...) {
    dependencies("magrittr")
    dependencies("piano", "install.from.bioconductor")

    stopifnot(all(test.configs$which.values %in% c("t-values", "p-values")))

    .to.df <- function(.df, .names, .col) {
        .df[,.col] %>% data.frame(row.names=as.data.frame(.df)[,.names]) %>% setNames(.col)
    }

    if (prefix.col != "") {
        t.values.col <- s.join(prefix.col, sep, t.values.col)
        p.values.col <- s.join(prefix.col, sep, p.values.col)
        log.fc.col <- s.join(prefix.col, sep, log.fc.col)
    }

    diffexp.df <- diffexp.df[order(diffexp.df[,p.values.col]),]
    diffexp.df <- diffexp.df[!duplicated(diffexp.df[,egid.col]),]

    if (is.null(ncpus)) {
        dependencies("parallel")
        ncpus <- detectCores(logical=TRUE)
        message("Number of CPUs not given. Setting to ", ncpus)
    }
    ncpus <- as.numeric(ncpus)
    if (is.null(nPerm)) {
        nPerm <- 10000
        message("Number of permutations not given. Setting to ", nPerm)
    }
    if ((nPerm %% ncpus) != 0) {
        message("No. of permutations (", nPerm,
                ") is not a multiple of the no. of CPUs (", ncpus,
                "). Setting the no. of permutations to ", (nPerm + (ncpus - (nPerm %% ncpus))))
        nPerm <- nPerm + (ncpus - (nPerm %% ncpus))
    }

    gsa.res <- list()
    for (i in 1:nrow(test.configs)) {
        message("Performing test: ", paste(test.configs[i,], collapse=", "))
        if (test.configs$which.values[i] == "t-values") {
            geneLevelStats <- .to.df(diffexp.df, egid.col, t.values.col)
            directions <- NULL
        } else { # --> (test.configs$which.values[i] == "p-values")
            geneLevelStats <- .to.df(diffexp.df, egid.col, p.values.col)
            directions <- .to.df(diffexp.df, egid.col, log.fc.col)
        }
        gsa.res[[test.configs$gene.set.stat[i]]] <-
            runGSA(geneLevelStats=geneLevelStats, directions=directions,
                   geneSetStat=test.configs$gene.set.stat[i],
                   signifMethod=test.configs$signif.method[i],
                   gsc=piano.gsets, ..., ncpus=ncpus, nPerm=nPerm)
    }
    gsa.res
}



plot.expression.levels <- function(expr.lvs, p.values, annotation=NULL,
                                   annotation.gset=NULL,
                                   prefix.fname, contrasts,
                                   annotate.pvals= FALSE,
                                   min.value=NULL, max.value=NULL) {
    dependencies(c("bigmemory", "NMF", "dplyr", "RColorBrewer"))
    stopifnot(all(paste(contrasts, "pVal") %in% colnames(p.values)))
    stopifnot(all(expr.lvs$egid %in% p.values$egid))

    signif.annot <- function(pvals, lims=c(0.001, 0.01, 0.05, 1)) {
        my.lims <- unique(c(sort(lims), 1))
        s.join("<=", sapply(pvals, function(x) min(my.lims[my.lims >= x])))
    }

    expr.mat <- expr.lvs %>% dplyr::select(-egid) %>% as.matrix()
    rownames(expr.mat) <- expr.lvs$egid

    annots <- p.values[match(rownames(expr.mat), p.values$egid), paste(contrasts, "pVal")]

    pval.lims <- c(0.01, 0.05, 1)
    pval.levels <- s.join("<=", pval.lims)
    annots <- lapply((annots %>% as.list()), signif.annot, lims=pval.lims) %>% data.frame(check.names=FALSE)

    pval.colors <- setNames(c("#99000d", "#ef3b2c", "#fff5f0"),
                            pval.levels)

    pval.colors.list <- list()
    for (pval.name in colnames(annots)) {
        pval.colors.list[[pval.name]] <- pval.colors[names(pval.colors) %in% annots[,pval.name]]
    }

    if (is.null(max.value)) {
        max.value <- max(expr.mat)
    }
    if (is.null(min.value)) {
        min.value <- min(expr.mat)
    }
    plot.to.file(s.join(prefix.fname, "heatmap-expr.png"), my.width=12, my.height=16)
    pal.name <- "RdYlBu"
    n.pal.colors <- 11
    heat.colors <- rev(colorRampPalette(c(brewer.pal(n.pal.colors, pal.name)[1],
                                          brewer.pal(n.pal.colors, pal.name),
                                          brewer.pal(n.pal.colors, pal.name)[n.pal.colors]))(50))
    heat.breaks <- compute.aheatmap.breaks(min.val=min.value, max.val=max.value, nstep=50)
    aheatmap(expr.mat,
             Colv=FALSE, labRow=NA, Rowv=FALSE,
             breaks=heat.breaks,
             color=heat.colors,
             annRow=if (annotate.pvals) annots else NULL,
             annColors=pval.colors.list,
             verbose=1,
             fontsize=7, cexCol=1
             )
    close.all.plots()
    expr.mat.df <- expr.mat %>% as.data.frame() %>% add_rownames("egid") %>% tbl_df() %>% mutate(egid=as.numeric(egid))
    if (!is.null(annotation)) {
        prev.nrows <- nrow(expr.mat.df)
        expr.mat.df <- annotation %>% inner_join(expr.mat.df)
        stopifnot(prev.nrows == nrow(expr.mat.df))
    }
    if (!is.null(annotation.gset)) {
        prev.nrows <- nrow(expr.mat.df)
        saved.egids <- expr.mat.df$egid
        expr.mat.df <- expr.mat.df %>% left_join(annotation.gset)
        expr.mat.df <- expr.mat.df[match(saved.egids, expr.mat.df$egid),]
        stopifnot(prev.nrows == nrow(expr.mat.df))
    }
    write.csv(expr.mat.df,
              file=s.join(prefix.fname, "heatmap-expr.csv"),
              na="", row.names=FALSE)
}
