# Function to prune pathways based on shared genes
prune_pathways <- function(tmp) {
  # Get all genes and pathways
  gen <- unique(unlist(lapply(tmp$intersection, function(c) {
    strsplit(c, ",")[[1]]
  })))
  path <- tmp$term_id # Pathways ordered by p-value

  # Create a matrix
  mat <- array(
    data = 0,
    dim = c(length(gen), length(path)),
    dimnames = list(gen, path)
  )

  # Fill the matrix with gene-pathway information
  for (j in 1:ncol(mat)) {
    p.gen <- strsplit(
      tmp$intersection[which(tmp$term_id == colnames(mat)[j])],
      ","
    )[[1]]
    mat[p.gen, j] <- 1
  }

  # Store the number of genes belonging to each pathway
  g.path <- colSums(mat)

  # Initialize empty list of pathways to include
  in.path <- c()

  # Continue while there are pathways remaining
  while (nrow(mat) > 0 & ncol(mat) > 0) {
    # Add the strongest remaining pathway (first pathway)
    in.path <- c(in.path, path[1])

    # Drop all genes belonging to this pathway
    mat <- mat[-which(mat[, path[1]] == 1), , drop = FALSE]

    # Drop empty and depleted pathways
    ii <- colSums(mat)
    ii <- ii / g.path
    mat <- mat[, which(ii > 0.5), drop = FALSE]

    # Reduce the list of pathways to look at
    path <- path[which(path %in% colnames(mat))]
    g.path <- g.path[which(names(g.path) %in% colnames(mat))]
  }

  # Report only the pathways included in in.path
  tmp <- tmp[tmp$term_id %in% in.path, ]

  return(tmp) # Return the pruned results
}
