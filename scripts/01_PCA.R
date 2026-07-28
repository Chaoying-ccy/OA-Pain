
# 01_PCA.R
#
# Description:
# This script performs Principal Component Analysis (PCA)
# separately for three monocyte subsets:
#
#   1. Classical monocytes (CLAS)
#   2. Intermediate monocytes (INT)
#   3. Non-classical monocytes (NONCLAS)
#
# PCA is performed using normalized miRNA expression data.
# Samples are coloured according to their pain group.
#
# The three PCA plots are combined into one figure and exported
# as a PDF file.
#
# NOTE:
# The input datasets are not included in the GitHub repository
# due to data access restrictions.
#
# Expected project structure:
#
# OA-Pain/
# ├── data/
# │   ├── miRNA Data Sample Key.csv
# │   └── miRNA_normalized_counts.csv
# ├── scripts/
# │   └── 01_PCA.R
# └── results/
#
# Required R packages:
#   ggplot2
#   dplyr
#   patchwork

# 1. Clear R environment

rm(list = ls())

# 2. Load required packages

library(ggplot2)
library(dplyr)
library(patchwork)

# 3. Define input and output paths
# All paths are relative to the OA-Pain project root.
# No personal/local absolute paths are used.
meta_file <-
  "data/miRNA Data Sample Key.csv"

counts_file <-
  "data/miRNA_normalized_counts.csv"

output_dir <-
  "results/PCA"


# Create output directory if it does not already exist

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# 4. Read input data

# Read sample metadata

meta <-
  read.csv(
    meta_file,
    check.names = FALSE
  )


# Read normalized miRNA expression matrix
# The first column is used as the miRNA identifier.

counts <-
  read.csv(
    counts_file,
    row.names = 1,
    check.names = FALSE
  )

# 5. Clean column names
# Remove accidental leading/trailing spaces from column names.

colnames(meta) <-
  trimws(
    colnames(meta)
  )

colnames(counts) <-
  trimws(
    colnames(counts)
  )


# 6. Identify metadata columns#
# Automatically identify columns corresponding to:
#   Sample ID
#   Cell type
#   Pain group
#

sample_col <-
  grep(
    "Sample",
    colnames(meta),
    value = TRUE
  )[1]

cell_col <-
  grep(
    "Cell|Type",
    colnames(meta),
    value = TRUE
  )[1]

group_col <-
  grep(
    "Group|Pain",
    colnames(meta),
    value = TRUE
  )[1]


# Rename these columns to standardised names

colnames(meta)[
  colnames(meta) == sample_col
] <- "SampleID"

colnames(meta)[
  colnames(meta) == cell_col
] <- "CellType"

colnames(meta)[
  colnames(meta) == group_col
] <- "Group"


# Remove accidental spaces from sample IDs

meta$SampleID <-
  trimws(
    meta$SampleID
  )


# ============================================================
# 7. Match metadata and expression samples
# ============================================================
#
# Only samples present in both the metadata and expression
# matrix are retained.
#

common <-
  intersect(
    meta$SampleID,
    colnames(counts)
  )


# Stop the analysis if no matching samples are found.

if (length(common) == 0) {
  stop(
    "No matching sample IDs were found between metadata and expression data."
  )
}


# Reorder metadata according to matched sample IDs

meta <-
  meta[
    match(
      common,
      meta$SampleID
    ),
    ,
    drop = FALSE
  ]


# Keep matched samples in expression matrix

counts <-
  counts[
    ,
    common,
    drop = FALSE
  ]


cat(
  "Matched samples:",
  length(common),
  "\n"
)


# 8. Prepare expression matrix

expr <-
  as.matrix(
    counts
  )


# Ensure expression values are numeric

mode(expr) <-
  "numeric"


# Remove miRNAs containing missing values

expr <-
  expr[
    complete.cases(expr),
    ,
    drop = FALSE
  ]


# Remove miRNAs with zero variance across all samples
# because they do not contribute to PCA.

expr <-
  expr[
    apply(
      expr,
      1,
      var
    ) > 0,
    ,
    drop = FALSE
  ]


cat(
  "miRNAs retained for PCA:",
  nrow(expr),
  "\n"
)


# ============================================================
# 9. Define plot colours
# ============================================================

cols <-
  c(
    "#4DAF4A",
    "#377EB8",
    "#E41A1C"
  )


# 10. PCA plotting function
# PCA is performed independently within each monocyte subset.
# Input:
#   ct = cell type (CLAS, INT or NONCLAS)
#
# Output:
#   ggplot PCA object
#

make_pca <-
  function(ct) {
    
    # Select metadata for the specified cell type

    meta_sub <-
      meta %>%
      filter(
        CellType == ct
      )
    
    
    # Check that samples exist for this cell type
    
    if (nrow(meta_sub) < 2) {
      stop(
        paste(
          "Insufficient samples for cell type:",
          ct
        )
      )
    }
    
    
    # --------------------------------------------------------
    # Select corresponding expression data
    # --------------------------------------------------------
    
    expr_sub <-
      expr[
        ,
        meta_sub$SampleID,
        drop = FALSE
      ]
    

    # Remove zero-variance miRNAs within this cell type
    
    expr_sub <-
      expr_sub[
        apply(
          expr_sub,
          1,
          var
        ) > 0,
        ,
        drop = FALSE
      ]
    
    
    # --------------------------------------------------------
    # Perform PCA
    # --------------------------------------------------------
    #
    # Samples are treated as observations.
    # miRNAs are treated as variables.
    #
    
    pca <-
      prcomp(
        t(expr_sub),
        center = TRUE,
        scale. = TRUE
      )
    

    # Calculate percentage of variance explained

    variance_explained <-
      round(
        100 *
          pca$sdev^2 /
          sum(
            pca$sdev^2
          ),
        1
      )
    

    # Create PCA plotting dataframe

    
    df <-
      data.frame(
        Group =
          meta_sub$Group,
        
        PC1 =
          pca$x[, 1],
        
        PC2 =
          pca$x[, 2]
      )
    

    # Define colours according to groups present

    
    group_levels <-
      sort(
        unique(
          df$Group
        )
      )
    
    group_colours <-
      setNames(
        cols[
          seq_along(
            group_levels
          )
        ],
        group_levels
      )
    

    # Generate PCA plot

    
    p <-
      ggplot(
        df,
        aes(
          x = PC1,
          y = PC2,
          colour = Group
        )
      ) +
      
      geom_point(
        size = 3.5,
        alpha = 0.85
      ) +
      
      stat_ellipse(
        aes(
          group = Group
        ),
        type = "t",
        linewidth = 0.9
      ) +
      
      scale_colour_manual(
        values =
          group_colours
      ) +
      
      theme_bw(
        base_size = 15
      ) +
      
      labs(
        title =
          ct,
        
        x =
          paste0(
            "PC1 (",
            variance_explained[1],
            "%)"
          ),
        
        y =
          paste0(
            "PC2 (",
            variance_explained[2],
            "%)"
          )
      ) +
      
      theme(
        plot.title =
          element_text(
            face = "bold",
            hjust = 0.5
          ),
        
        legend.position =
          "right"
      )
    
    
    return(p)
  }


# ============================================================
# 11. Generate PCA plots
# ============================================================

p1 <-
  make_pca(
    "CLAS"
  )

p2 <-
  make_pca(
    "INT"
  )

p3 <-
  make_pca(
    "NONCLAS"
  )


# 12. Combine the three PCA plots
# The legends are collected into a single legend on the right.

final <-
  (
    p1 +
      p2 +
      p3
  ) +
  
  plot_layout(
    guides = "collect"
  ) +
  
  plot_annotation(
    title =
      "PCA of miRNA Expression by Cell Type"
  ) &
  
  theme(
    legend.position =
      "right",
    
    plot.title =
      element_text(
        size = 18,
        face = "bold",
        hjust = 0.5
      )
  )


# 13. Export PCA figure
ggsave(
  filename =
    file.path(
      output_dir,
      "PCA_All_CellTypes.pdf"
    ),
  
  plot =
    final,
  
  width =
    18,
  
  height =
    6
)


# 14. Finish

cat(
  "\nPCA analysis completed successfully.\n"
)

cat(
  "Output:",
  file.path(
    output_dir,
    "PCA_All_CellTypes.pdf"
  ),
  "\n"
)