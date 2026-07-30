# 02_DE_Volcano.R
# Differential Expression Analysis and Volcano Plots

# Description:
# This script performs differential expression (DE) analysis
# of miRNA expression between High Pain and Low Pain groups
# within three monocyte subsets:
#
#   1. Classical monocytes (CLAS)
#   2. Intermediate monocytes (INT)
#   3. Non-classical monocytes (NONCLAS)
#
# Differential expression analysis is performed using limma.
#
# Comparison:
#   High Pain vs Low Pain
#
# Therefore:
#   Positive logFC = higher expression in High Pain
#   Negative logFC = higher expression in Low Pain
#
# Exploratory candidate threshold:
#   P.Value < 0.05
#   |logFC| > 0.5
# NOTE:
# The input datasets are not included in this repository
# due to data access restrictions.

#
# Required R packages:
#   limma
#   dplyr
#   ggplot2
#   ggrepel

# 1. Clear R environment

rm(list = ls())

# 2. Load required packages

library(limma)
library(dplyr)
library(ggplot2)
library(ggrepel)


# 3. Define input and output paths
# All paths are relative to the OA-Pain project root.

counts_file <-
  "data/miRNA_normalized_counts.csv"

meta_file <-
  "data/miRNA Data Sample Key.csv"

de_dir <-
  "results/DE"


# Create output directory if necessary

dir.create(
  de_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# 4. Read input data

counts <-
  read.csv(
    counts_file,
    check.names = FALSE
  )

meta <-
  read.csv(
    meta_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )


# 5. Prepare expression matrix
# Use the first column as miRNA identifiers

rownames(counts) <-
  trimws(
    counts[, 1]
  )

counts <-
  counts[, -1]


# Convert expression data to numeric matrix

expr <-
  as.matrix(
    counts
  )

mode(expr) <-
  "numeric"


# Clean sample names

colnames(expr) <-
  trimws(
    colnames(expr)
  )


# 6. Prepare metadata

# Clean metadata column names

colnames(meta) <-
  trimws(
    colnames(meta)
  )


# Automatically identify relevant metadata columns

sample_col <-
  grep(
    "Sample",
    colnames(meta),
    value = TRUE
  )[1]

cell_col <-
  grep(
    "Cell",
    colnames(meta),
    value = TRUE
  )[1]

group_col <-
  grep(
    "Group",
    colnames(meta),
    value = TRUE
  )[1]


# Rename columns to standardised names

colnames(meta)[
  colnames(meta) == sample_col
] <- "SampleID"

colnames(meta)[
  colnames(meta) == cell_col
] <- "CellType"

colnames(meta)[
  colnames(meta) == group_col
] <- "Group"


# Clean metadata values

meta$SampleID <-
  trimws(
    meta$SampleID
  )

meta$CellType <-
  trimws(
    meta$CellType
  )

meta$Group <-
  trimws(
    meta$Group
  )


# 7. Match metadata with expression data
common <-
  intersect(
    colnames(expr),
    meta$SampleID
  )


# Stop if no samples match

if (length(common) == 0) {
  
  stop(
    "No matching sample IDs were found between metadata and expression data."
  )
  
}


# Keep matched samples only

expr <-
  expr[
    ,
    common,
    drop = FALSE
  ]


# Reorder metadata to match expression matrix

meta <-
  meta[
    match(
      common,
      meta$SampleID
    ),
    ,
    drop = FALSE
  ]


# Confirm correct sample order

stopifnot(
  all(
    colnames(expr) ==
      meta$SampleID
  )
)


cat(
  "\nMatched samples:",
  ncol(expr),
  "\n"
)

# 8. Differential expression function

run_DE <-
  function(cell) {
    
    cat(
      "\nRunning DE analysis:",
      cell,
      "\n"
    )
    
    
    # --------------------------------------------------------
    # Select cell type and remove Pain-Free Controls
    # --------------------------------------------------------
    
    meta_sub <-
      meta %>%
      filter(
        CellType == cell,
        Group != "Pain-Free Control"
      )
    
    
    # --------------------------------------------------------
    # Extract corresponding expression data
    # --------------------------------------------------------
    
    expr_sub <-
      expr[
        ,
        meta_sub$SampleID,
        drop = FALSE
      ]
    
    
    # --------------------------------------------------------
    # Remove miRNAs with zero variance
    # --------------------------------------------------------
    
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
    # Define pain groups
    # --------------------------------------------------------
    #
    # The original metadata label:
    #
    #   "High 17-HDHA - Low Pain"
    #
    # is defined as LowPain.
    #
    # All remaining non-control samples are defined as HighPain.
    #
    
    group <-
      ifelse(
        trimws(
          meta_sub$Group
        ) ==
          "High 17-HDHA - Low Pain",
        
        "LowPain",
        
        "HighPain"
      )
    
    
    # LowPain is the reference group.
    #
    # Therefore the limma coefficient represents:
    #
    #   HighPain - LowPain
    #
    # Positive logFC = High Pain higher
    # Negative logFC = Low Pain higher
    
    group <-
      factor(
        group,
        levels =
          c(
            "LowPain",
            "HighPain"
          )
      )
    
    
    cat(
      "\nGroup distribution:\n"
    )
    
    print(
      table(group)
    )
    
    
    # Ensure both comparison groups are present
    
    if (length(unique(group)) < 2) {
      
      stop(
        paste(
          "Both HighPain and LowPain samples are required for",
          cell
        )
      )
      
    }
    
    
    # --------------------------------------------------------
    # Create design matrix
    # --------------------------------------------------------
    #
    # This section follows the standard limma workflow for
    # differential expression analysis using a design matrix,
    # lmFit(), eBayes(), and topTable().
    # Reference:
    # https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf
    #
    # Modifications were made for the present dataset by comparing
    # High Pain and Low Pain samples separately within each monocyte subset.
    
    design <-
      model.matrix(
        ~ group
      )
    
    
    # --------------------------------------------------------
    # Fit linear model using limma
    # --------------------------------------------------------
    
    fit <-
      lmFit(
        expr_sub,
        design
      )
    
    fit <-
      eBayes(
        fit
      )
    
    
    # --------------------------------------------------------
    # Extract DE results
    # --------------------------------------------------------
    
    res <-
      topTable(
        fit,
        coef = 2,
        number = Inf,
        sort.by = "P"
      )
    
    
    # Add miRNA identifiers
    
    res <-
      data.frame(
        miRNA =
          rownames(res),
        
        res,
        
        row.names = NULL
      )
    
    
    # Keep relevant columns
    
    res <-
      res[
        ,
        c(
          "miRNA",
          "logFC",
          "AveExpr",
          "P.Value",
          "adj.P.Val"
        )
      ]
    
    
    # 9. Candidate classification
    #
    # Exploratory threshold:
    #
    #   P.Value < 0.05
    #   |logFC| > 0.5
    #
    
    res$status <-
      "NotSig"
    
    
    # Higher expression in High Pain
    
    res$status[
      res$P.Value < 0.05 &
        res$logFC > 0.5
    ] <-
      "HighPain"
    
    
    # Higher expression in Low Pain
    
    res$status[
      res$P.Value < 0.05 &
        res$logFC < -0.5
    ] <-
      "LowPain"
    
    
    # --------------------------------------------------------
    # Print DE summary
    # --------------------------------------------------------
    
    cat(
      "\nCandidate summary:",
      cell,
      "\n"
    )
    
    print(
      table(
        res$status
      )
    )
    
    
    # 10. Save DE result
    
    # Save complete DE table
    
    write.csv(
      res,
      
      file.path(
        de_dir,
        paste0(
          cell,
          "_DE_all.csv"
        )
      ),
      
      row.names = FALSE
    )
    
    
    # Save candidate miRNAs only
    
    write.csv(
      subset(
        res,
        status != "NotSig"
      ),
      
      file.path(
        de_dir,
        paste0(
          cell,
          "_candidates.csv"
        )
      ),
      
      row.names = FALSE
    )
    
    
    return(res)
  }


# 11. Volcano plot function

draw_volcano <-
  function(
    res,
    cell
  ) {
    

    # Select top 10 candidate miRNAs
    # run_DE() already sorts the table by P.Value.
    # Therefore head(10) selects the ten candidates with
    # the smallest nominal P-values.
    
    top <-
      res %>%
      filter(
        status != "NotSig"
      ) %>%
      head(10)
    
    # Generate volcano plot
    #
    # This section follows a standard volcano plot approach, plotting
    # log fold change against -log10(P.Value) and highlighting candidate
    # miRNAs using project-specific thresholds.
    # Reference:
    # https://ggplot2.tidyverse.org/reference/geom_point.html
    #
    # The miRNA labels are added using ggrepel to reduce label overlap.
    # Reference:
    # https://ggrepel.slowkow.com/reference/geom_text_repel.html
    #
    # Modifications were made to colour candidates by pain-group direction
    # and label the top candidate miRNAs for each monocyte subset.
    
    p <-
      ggplot(
        res,
        aes(
          x = logFC,
          y = -log10(P.Value)
        )
      ) +
      
      geom_point(
        aes(
          color = status
        ),
        size = 2
      ) +
      
      scale_color_manual(
        values =
          c(
            HighPain = "#619CFF",
            LowPain = "#F8766D",
            NotSig = "grey75"
          )
      ) +
      
      # logFC thresholds
      
      geom_vline(
        xintercept =
          c(
            -0.5,
            0.5
          ),
        linetype = 2
      ) +
      
      # Nominal P-value threshold
      
      geom_hline(
        yintercept =
          -log10(
            0.05
          ),
        linetype = 2
      ) +
      
      # Label top candidate miRNAs
      
      geom_text_repel(
        data = top,
        aes(
          label = miRNA
        ),
        size = 4,
        max.overlaps = Inf
      ) +
      
      theme_bw(
        base_size = 18
      ) +
      
      labs(
        title =
          paste0(
            cell,
            ": High Pain vs Low Pain"
          ),
        
        x =
          "logFC (High Pain - Low Pain)",
        
        y =
          expression(
            -log[10](P)
          ),
        
        colour =
          "Expression"
      )
    

    
    ggsave(
      filename =
        file.path(
          de_dir,
          paste0(
            "Volcano_",
            cell,
            ".pdf"
          )
        ),
      
      plot = p,
      
      width = 10,
      
      height = 8
    )
    
  }


# 13. Run DE analysis

CLAS <-
  run_DE(
    "CLAS"
  )

INT <-
  run_DE(
    "INT"
  )

NONCLAS <-
  run_DE(
    "NONCLAS"
  )


# 14. Generate volcano plots

draw_volcano(
  CLAS,
  "CLAS"
)

draw_volcano(
  INT,
  "INT"
)

draw_volcano(
  NONCLAS,
  "NONCLAS"
)


# 15. Finish

cat(
  "\nDE analysis completed successfully.\n"
)

cat(
  "\nResults saved to:",
  de_dir,
  "\n"
)
