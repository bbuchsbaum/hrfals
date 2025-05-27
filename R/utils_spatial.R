#' Build a Voxel Graph Laplacian
#'
#' Constructs an unnormalised graph Laplacian matrix for the voxels
#' contained in a 3D mask volume. The function accepts objects that can
#' be coerced to an array (e.g. a `neuroim2::NeuroVol`). Non-zero
#' elements of the volume are treated as included voxels.
#'
#' @param volume 3D array or `neuroim2::NeuroVol` mask defining the
#'   spatial domain.
#' @param connectivity Neighbourhood definition. Either 6, 18 or 26.
#' @return List with sparse Laplacian matrix `L` (`dgCMatrix`) and the
#'   numeric degree vector `degree`.
#' @export
build_voxel_laplacian <- function(volume, connectivity = 6) {
  vol <- as.array(volume)
  if (length(dim(vol)) != 3) {
    stop("'volume' must be 3D")
  }
  mask <- vol != 0
  dims <- dim(mask)
  coords <- which(mask, arr.ind = TRUE)
  v <- nrow(coords)
  if (v == 0) {
    stop("mask contains no voxels")
  }
  index_map <- array(0L, dims)
  index_map[matrix(t(coords), ncol = 3)] <- seq_len(v)

  offsets <- expand.grid(dx = -1:1, dy = -1:1, dz = -1:1)
  offsets <- offsets[!(offsets$dx == 0 & offsets$dy == 0 & offsets$dz == 0), ]
  if (connectivity == 6) {
    offsets <- offsets[rowSums(abs(offsets)) == 1, ]
  } else if (connectivity == 18) {
    offsets <- offsets[rowSums(abs(offsets)) %in% c(1, 2), ]
  } else if (connectivity == 26) {
    offsets <- offsets
  } else {
    stop("connectivity must be one of 6, 18, 26")
  }

  nnz <- 0
  rows <- integer(0)
  cols <- integer(0)

  for (i in seq_len(nrow(offsets))) {
    off <- as.integer(offsets[i, ])
    nbr_coords <- sweep(coords, 2, off, "+")
    inside <- nbr_coords[, 1] >= 1 & nbr_coords[, 1] <= dims[1] &
              nbr_coords[, 2] >= 1 & nbr_coords[, 2] <= dims[2] &
              nbr_coords[, 3] >= 1 & nbr_coords[, 3] <= dims[3]
    if (!any(inside)) next
    nbr_idx <- index_map[matrix(t(nbr_coords[inside, , drop = FALSE]), ncol = 3)]
    valid <- nbr_idx > 0L
    if (!any(valid)) next
    from <- which(inside)[valid]
    to <- nbr_idx[valid]
    rows <- c(rows, from)
    cols <- c(cols, to)
    nnz <- nnz + length(to)
  }

  if (nnz == 0) {
    Adj <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                                dims = c(v, v))
  } else {
    Adj <- Matrix::sparseMatrix(i = c(rows, cols), j = c(cols, rows),
                                x = 1, dims = c(v, v))
  }
  degree <- Matrix::rowSums(Adj)
  L <- Matrix::Diagonal(x = degree) - Adj
  list(L = L, degree = as.numeric(degree))
}

