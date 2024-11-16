library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v75)
library(tidyverse)
BiocManager::install("biovizBase")


# what is a fragment file? How is it generated?
# https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments
frag.file <- read.delim('data/10k_pbmc_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz', header = F, nrows = 10)
head(frag.file)


# 1. Read in data -----------------

counts <- Read10X_h5('data/10k_pbmc_ATACv2_nextgem_Chromium_Controller_filtered_peak_bc_matrix.h5')
counts[1:10,1:10]


CreateChromatinAssayTest <- function(
    counts,
    data,
    min.cells = 0,
    min.features = 0,
    max.cells = NULL,
    ranges = NULL,
    motifs = NULL,
    fragments = NULL,
    genome = NULL,
    annotation = NULL,
    bias = NULL,
    positionEnrichment = NULL,
    sep = c("-", "-"),
    validate.fragments = TRUE,
    verbose = TRUE,
    ...
) {
  if (missing(x = counts) && missing(x = data)) {
    stop("Must provide either 'counts' or 'data'")
  } else if (!missing(x = counts) && !missing(x = data)) {
    stop("Either 'counts' or 'data' must be missing; both cannot be provided")
  } else if (!missing(x = counts)) {
    data.use <- counts
  } else {
    data.use <- data
  }
  if (!is.null(x = ranges)) {
    if (length(x = ranges) != nrow(x = data.use)) {
      stop("Length of supplied genomic ranges does not match number
           of rows in matrix")
    }
  } else {
    ranges <- StringToGRanges(regions = rownames(x = data.use), sep = sep)
  }
  if (!isDisjoint(x = ranges)) {
    warning("Overlapping ranges supplied. Ranges should be non-overlapping.")
  }
  if (!is.null(x = annotation) & !inherits(x = annotation, what = "GRanges")) {
    stop("Annotation must be a GRanges object.")
  }
  # remove low-count cells
  ncount.cell <- colSums(x = data.use > 0)
  data.use <- data.use[, ncount.cell >= min.features]
  
  if (ncol(x = data.use) == 0) {
    stop("No cells retained due to minimum feature cutoff supplied")
  }
  
  ncell.feature <- rowSums(x = data.use > 0)
  if (!is.null(x = max.cells)) {
    if (is(object = max.cells, class2 = "character")) {
      percent.cutoff <- as.numeric(
        x = gsub(pattern = "q", replacement = "", x = max.cells)
      )
      max.cells <- (percent.cutoff / 100) * ncol(x = data.use)
    }
  } else {
    max.cells <- ncol(x = data.use)
  }
  features.keep <- (ncell.feature >= min.cells) & (ncell.feature <= max.cells)
  if (sum(features.keep) == 0) {
    stop("No features retained due to minimum cell cutoff supplied")
  }
  data.use <- data.use[features.keep, ]
  ranges <- ranges[features.keep, ]
  # re-assign row names of matrix so that it's a known granges transformation
  new.rownames <- GRangesToString(grange = ranges, sep = c("-", "-"))
  rownames(x = data.use) <- new.rownames
  if (!missing(x = counts)) {
    seurat.assay <- CreateAssayObject(
      counts = data.use,
      data = data,
      min.cells = -1,
      min.features = -1 # min cell/feature filtering already done
    )
  } else {
    seurat.assay <- CreateAssayObject(
      counts = counts,
      data = data.use,
      min.cells = min.cells,
      min.features = min.features
    )
  }
  if (inherits(x = fragments, what = "list")) {
    # check each object in the list is a fragment object
    # fragment list usually supplied when doing object merge,
    # so don't validate cells here, we can assume that was done in
    # individual object creation
    obj.class <- sapply(
      X = fragments, FUN = function(x) inherits(x = x, what = "Fragment")
    )
    if (!all(obj.class)) {
      stop("All objects in fragments list must be Fragment-class objects")
    }
    frags <- lapply(
      X = fragments,
      FUN = AssignFragCellnames,
      cellnames = colnames(x = seurat.assay)
    )
    # subset to cells in the assay
    frags <- lapply(
      X = fragments,
      FUN = subset,
      cells = colnames(x = seurat.assay)
    )
  } else if (inherits(x = fragments, what = "Fragment")) {
    # single Fragment object supplied
    frags <- AssignFragCellnames(
      fragments = fragments, cellnames = colnames(x = seurat.assay)
    )
    # subset to cells in the assay
    frags <- subset(x = frags, cells = colnames(x = seurat.assay))
  } else {
    # path to fragment file supplied, create fragment object
    frags <- list()
    if (!is.null(x = fragments)) {
      if (nchar(x = fragments) > 0) {
        cells <- colnames(x = seurat.assay)
        names(x = cells) <- cells
        frags[[1]] <- CreateFragmentObject(
          path = fragments,
          cells = cells,
          validate.fragments = validate.fragments,
          verbose = verbose,
          ...
        )
      }
    }
  }
  
  if (!is.null(x = motifs)) {
    # pre-computed motif object, make sure features are formatted the same
    # as peak matrix and subset features
    if (!inherits(x = motifs, what = "Motif")) {
      stop("Provided motif object is not a Motif-class object")
    }
    if (!(all(rownames(x = motifs) == rownames(x = seurat.assay)))) {
      # rownames don't match
      motif.mat <- GetMotifData(object = motifs)
      motif.granges <- StringToGRanges(
        regions = rownames(x = motifs), sep = sep
      )
      rownames(x = motif.mat) <- GRangesToString(grange = motif.granges)
      # subset
      if (!all(rownames(x = seurat.assay) %in% rownames(x = motif.mat))) {
        warning("Some peak regions missing from supplied motif object. ",
                "Motif information will not be added")
        motifs <- NULL
      }
      motif.mat <- motif.mat[rownames(x = seurat.assay), ]
      motifs <- SetMotifData(
        object = motifs, slot = "data", new.data = motif.mat
      )
    }
  }
  chrom.assay <- as.ChromatinAssay(
    x = seurat.assay,
    ranges = ranges,
    seqinfo = genome,
    motifs = motifs,
    fragments = frags,
    annotation = annotation,
    bias = bias,
    positionEnrichment = positionEnrichment
  )
  return(chrom.assay)
}
chrom_assay <- CreateChromatinAssayTest(
  counts = counts,
  sep = c(":", "-"),
  fragments = "data/10k_pbmc_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

str(chrom_assay)

metadata <- read.csv(file = 'data/10k_pbmc_ATACv2_nextgem_Chromium_Controller_singlecell.csv', header = T, row.names = 1)
View(metadata)


# create a seurat Object
pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  meta.data = metadata,
  assay = 'ATAC'
)

str(pbmc)


# ....Adding Gene Annotation -------------------

pbmc@assays$ATAC@annotation
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)


# change to UCSC style since the data was mapped to hg19
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))


# add the gene information to the object
Annotation(pbmc) <- annotations
pbmc@assays$ATAC@annotation



# 2. Computing QC ---------------------

# compute nucleosome signal score per cell
pbmc <- NucleosomeSignal(pbmc)

# compute TSS enrichment score per cell
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100

View(pbmc@meta.data)

# ....Visualizing QC --------------------

colnames(pbmc@meta.data)
a1 <- DensityScatter(pbmc, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
a2 <- DensityScatter(pbmc, x = 'nucleosome_signal', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
a1
a1 | a2

VlnPlot(object = pbmc, 
        features = c('nCount_ATAC', 'nFeature_ATAC', 'TSS.enrichment', 'nucleosome_signal', 'blacklist_ratio', 'pct_reads_in_peaks'),
        pt.size = 0.1,
        ncol = 6)


# ....Filtering poor quality cells --------------------
# look at above plot to choose following values
pbmc <- subset(x = pbmc,
               subset = nCount_ATAC > 3000 &
                 nCount_ATAC < 40000 &
                 pct_reads_in_peaks > 15 & 
                 blacklist_ratio < 0.05 &
                 nucleosome_signal < 4 &
                 TSS.enrichment > 1)



# 3. Normalization and linear dimensional reduction ------------------
pbmc <- RunTFIDF(pbmc) # normalization
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0') # selecting top features
pbmc <- RunSVD(pbmc) # dimensionality reduction

DepthCor(pbmc)

# 4. Non-linear dimensional reduction and Clustering -------------------
# Get rid of first dimension. Chose by looking at the above plot
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, algorithm = 3)

DimPlot(object = pbmc, label = TRUE) + NoLegend()
