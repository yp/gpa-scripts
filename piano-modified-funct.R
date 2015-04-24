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
## This script provides an alternative definition of function piano::runQC
## in order to address some graphics issues
my.print.plot <- function(p, ...) {
    if (!is.null(p)) {
        dev.new();
        print(p, ...);
    }
    invisible(NULL);
}

my.piano.runQC <- function (arrayData,
                            rnaDeg = TRUE,
                            nuseRle = TRUE,
                            hist = TRUE, 
                            boxplot = TRUE,
                            pca = TRUE,
                            plot.fnt = my.print.plot,
                            ...
                            ) 
{
    dependencies(c("plyr", "dplyr", "ggplot2", "reshape2"))
    dependencies(c("affyPLM", "affy"), "install.from.bioconductor")
    my.setup <- data.frame(Sample=rownames(arrayData$setup),
                           Condition=do.call("paste", arrayData$setup),
                           stringsAsFactors=TRUE)
    .runRnaDeg <- function(arrayData, plot.fnt, ...) {
        .verb("Generating RNA degradation plot...")
        rna.deg.obj <- AffyRNAdeg(arrayData$dataRaw, log.it = TRUE)
        mns <- rna.deg.obj$means.by.number
        sds <- rna.deg.obj$ses
        mn <- mns[, 1]
        mns <- sweep(mns, 1, mn)
        mns <- mns/(sds)
        mns <- sweep(mns, 1, 1:(dim(mns)[1]), "+")
        x <- as.data.frame(t(mns))
        colnames(x) <- rna.deg.obj$sample.names
        x$X <- 0:(nrow(x)-1)
        x <- melt(x, id.var="X", variable.name="Sample")
        x <- inner_join(x, my.setup)
        rdp <- (  ggplot(x, aes(x=X, y=value, group=Sample, colour=Condition))
                + geom_line()
                + scale_x_continuous("5' <-----> 3'\n Probe Number")
                + scale_y_continuous("Mean Intensity : shifted and scaled")
                + theme(legend.position = "bottom"))
        plot.fnt(rdp, plot.name="dna-deg", ...)
        .verb("...done")
    }
    .runNuseRle <- function(arrayData, plot.fnt, ...) {
        .verb("Generating NUSE and RLE plot...")
        Pset <- fitPLM(arrayData$dataRaw)
        plot.fnt(NULL, ...)
        par(mar = c(10, 4, 4, 2))
        RLE(Pset, main = "RLE (Relative log expression)", 
            las = 2)
        plot.fnt(NULL, ...)
        par(mar = c(10, 4, 4, 2))
        NUSE(Pset, main = "NUSE (Normalized unscaled standard error)", 
             las = 2)
        .verb("...done")
    }
    .runHist <- function(arrayData, plot.fnt,
                         width = par("din")[1], height = par("din")[2],
                         fixed.size = FALSE, ...) {
        expr.raw <- data.frame()
        if ("dataRaw" %in% attributes(arrayData)$names) {
            .verb("Generating raw data distribution plot...")
            expr.raw <- melt(log2(exprs(arrayData$dataRaw)))
            expr.raw$Var1 <- NULL
            colnames(expr.raw)[colnames(expr.raw) == "Var2"] <- "Sample"
            expr.raw$Status <- factor("Raw", levels=c("Raw", "Normalized"))
            .verb("...done")
        }
        .verb("Generating normalized data distribution plot...")
        expr.norm <- melt(arrayData$dataNorm, id.vars=c(), variable.name="Sample")
        expr.norm$Status <- factor("Normalized", levels=c("Raw", "Normalized"))
        expr.tot <- inner_join(my.setup, rbind(expr.raw, expr.norm))
        rm(expr.raw, expr.norm)
        gc()

        dp <- (  ggplot(expr.tot, aes(x=value, group=Sample, colour=Condition))
               + geom_density()
               + facet_grid(~Status)
               + scale_x_continuous(expression(log[2]~intensity)))
        width <- if (!fixed.size) 1.5*width else width;
        plot.fnt(dp, plot.name="expr-val-hist", width=width, height=height, ...)
        rm(expr.tot)
        .verb("...done")
    }
    .runBoxplot <- function(arrayData, plot.fnt,
                            width = par("din")[1], height = par("din")[2],
                            fixed.size = FALSE, ...) {
        dd.raw <- data.frame()
        if ("dataRaw" %in% attributes(arrayData)$names) {
            .verb("Generating raw data boxplot...")
            dd.raw <- cbind(Status="Raw",
                            adply(log2(exprs(arrayData$dataRaw)), 2, quantile, c(0, .25, .5, .75, 1)))
            .verb("...done")
        }
        .verb("Generating normalized data boxplot...")
        dd.norm <- cbind(Status="Normalized",
                         adply(as.matrix(arrayData$dataNorm), 2, quantile, c(0, .25, .5, .75, 1), na.rm=TRUE))
        dd <- rbind(dd.raw, dd.norm)
        .verb("Generating normalized data boxplot...")
        dd$Status <- factor(dd$Status, levels=c("Raw", "Normalized"))
        colnames(dd)[colnames(dd) == "X1"] <- "Sample"
        dd <- inner_join(dd, my.setup)
        bp <- (  ggplot(dd)
               + geom_boxplot(aes(x=Sample,
                                  ymin=`0%`, lower = `25%`, middle = `50%`, upper = `75%`, ymax = `100%`,
                                  fill=Condition),
                              stat="identity")
               + scale_x_discrete("Samples") + scale_y_continuous(expression(log[2]~intensity))
               + facet_grid(~Status)
               + theme(legend.position = "bottom"))
        width <- if (!fixed.size) 1.5*width else width;
        plot.fnt(bp, plot.name="expr-val-boxpl", width=width, height=height, ...)
        .verb("...done")
    }
    .runPca <- function(arrayData, plot.fnt, ...) {
        dependencies("gridExtra")
        .verb("Generating PCA...")
        dataForPca <- arrayData$dataNorm
        dataForPca <- dataForPca - rowMeans(dataForPca, na.rm=TRUE)
        dataPrcomp <- prcomp(na.omit(dataForPca))
        dd <- cbind(data.frame(Sample=rownames(dataPrcomp$rotation)),
                    dataPrcomp$rotation[,1:3])
        dd <- inner_join(dd, my.setup)
        p1 <- (  ggplot(dd, aes(x=PC1, y=PC2, size=PC3, colour=Condition))
               + geom_point()
               + scale_size(range=c(3, 9))
               + geom_text(aes(label=Sample), colour="black", size=4)
               )

        dd <- as.data.frame(t(summary(dataPrcomp)$importance))
        dd$Dimension <- rownames(dd)
        p2 <- (  ggplot(dd, aes(x=factor(Dimension, levels=Dimension), y=100*`Proportion of Variance`))
               + geom_bar(stat="identity", width=.6)
               + geom_text(aes(y=3+100*`Proportion of Variance`,
                               label=sprintf("%1.1f%%", 100*`Proportion of Variance`)),
                           size=4)
               + scale_x_discrete("") + scale_y_continuous("Proportion of variance")
               )
        ptot <- arrangeGrob(p1, p2, ncol=1, heights=c(0.4, 0.2))
        plot.fnt(ptot, plot.name="pca", ...)
        .verb("...done")
    }
    .verb <- message
    if (class(arrayData) != "ArrayData") {
        stop("argument arrayData is not of class ArrayData")
    }
    if (rnaDeg) {
        if ("dataRaw" %in% attributes(arrayData)$names) {
            .runRnaDeg(arrayData, plot.fnt = plot.fnt, ...)
        } else {
            warning("can not run rnaDeg: argument arrayData does not contain dataRaw")
        }
    }
    if (nuseRle) {
        if ("dataRaw" %in% attributes(arrayData)$names) {
            .runNuseRle(arrayData, plot.fnt = plot.fnt, ...)
        } else {
            warning("can not run nuseRle: argument arrayData does not contain dataRaw")
        }
    }
    if (hist) {
        .runHist(arrayData, plot.fnt = plot.fnt, ...)
    }
    if (boxplot) {
        .runBoxplot(arrayData, plot.fnt = plot.fnt, ...)
    }
    if (pca) {
        .runPca(arrayData, plot.fnt = plot.fnt, ...)
    }
}






unlockBinding("runQC", env=as.environment("package:piano"))
assign("runQC", my.piano.runQC, envir=as.environment("package:piano"))
lockBinding("runQC", env=as.environment("package:piano"))
