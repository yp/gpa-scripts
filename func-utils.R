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
## This script provides some convenience general-purpose functions


`%ni%` <- Negate(`%in%`)




plot.to.file <- function(filename,
                         my.width=8, my.height=8,
                         my.ps=7) {
    if (dev.cur() > 1)
        dev.off()
    png(filename=filename,
        width=100*my.width, height=100*my.height,
        res=200, pointsize=my.ps, type="cairo-png", antialias="subpixel")
    invisible()
}




duplicated.values <- function(x, collapse=FALSE) {
    if (collapse)
        unique(x[duplicated(x)])
    else
        x[duplicated(x)]
}




install.from.bioconductor <- function(lib) {
    source("http://bioconductor.org/biocLite.R")
    biocLite(lib)
}




load.dependent.library <- function(lib, inst.func=install.packages) {
    loaded <- require(lib, character.only=TRUE)
    if (loaded) {
        message("'", lib, "' loaded.")
    } else {
        message("Trying to install '", lib,"'... with ", inst.func)
        eval(call(inst.func, lib))
        loaded <- require(lib, character.only=TRUE)
    }
    loaded
}





dependencies <- function(libs, inst.func="install.packages") {
    libs.df <- data.frame(lib=libs, inst.func=inst.func)
    apply(libs.df, 1,
          function(x) { load.dependent.library(x[1], x[2]); NULL})
    invisible()
}




enrich.from.biomart <- function(df, key.field, biomart.key.field,
                                attributes, attributes.names=attributes,
                                biomart.dataset="mmusculus_gene_ensembl") {
    dependencies("dplyr")
    dependencies("biomaRt", "install.from.bioconductor")

    ensembl <- useMart("ensembl", dataset=biomart.dataset)
    res <- getBM(attributes=c(biomart.key.field, attributes),
                 filters = biomart.key.field,
                 values = df[, key.field],
                 mart = ensembl)
    names(res) <- c(key.field, attributes.names)
    left_join(df, res)
}




compute.with.caching <- function(envir, data.name, file.name, recompute.from.scratch, FUN, ...) {
    dependencies(c("digest"))
    detailed.input.digest <- lapply(list(...), digest, algo="sha256")
    input.digest <- digest(object=list(...), algo="sha256")
    succesfully.loaded <- FALSE
    q.file.name <- file.path('rdata', file.name)
    if (!recompute.from.scratch && file.exists(q.file.name)) {
        message("Loading '", data.name, "' from file '", q.file.name, "'")
        load(q.file.name, envir=envir)
        if (!exists("__hash_sha256__", envir=envir, inherits=FALSE)) {
            message("No hash found in file '", q.file.name, "'. Recomputing the data...")
        } else {
            if (envir[["__hash_sha256__"]] != input.digest) {
                message("Read hash value '", envir[["__hash_sha256__"]], "' does not match ",
                        "the input hash value '", input.digest, "'. Recomputing the data...");
                if (exists("__detailed_hash_sha256__", envir=envir, inherits=FALSE)) {
                    message("Detailed input hash: ", detailed.input.digest)
                    message("Detailed read hash: ", envir[["__detailed_hash_sha256__"]])
                }
            } else {
                message("Read hash value and input hash values match!")
                succesfully.loaded <- TRUE
            }                                  
        }
    }
    if (!succesfully.loaded) {
        message("Recomputing '", data.name, "' and saving to file '", q.file.name, "'")
        envir[[data.name]] <- FUN(...)
        envir[["__hash_sha256__"]] <- input.digest
        envir[["__detailed_hash_sha256__"]] <- detailed.input.digest
        save(list=c(data.name, "__detailed_hash_sha256__", "__hash_sha256__"), envir=envir, file=q.file.name)
    }
    invisible()
}




close.all.plots <- function() {
    while (dev.cur() != 1) {
        dev.off()
    }
}




s.join <- function(...) paste(..., sep="")



compute.aheatmap.breaks <- function(mat=NULL, min.val=NULL, max.val=NULL, nstep=100) {
    .round.to.mult <- function(num, mult=0.5, round.fn=floor) {
        round.fn(num/mult)*mult
    }
    lims <- 10^(-8:8)
    lims <- sort(c(lims, lims*2, lims*2.5, lims*5))
    if (is.null(mat)) {
        min.value <- min.val
        max.value <- max.val
    } else {
        min.value <- min(mat, na.rm=TRUE) ## strict
        max.value <- max(mat, na.rm=TRUE) ## strict
    }
    print(min.value)
    print(max.value)
    break.step <- min(lims[lims >= (max.value-min.value)/6])
    min.value <- .round.to.mult(min.value, break.step, floor)
    max.value <- .round.to.mult(max.value, break.step, ceiling)
    seq(min.value, max.value, length.out=nstep+1)
}
