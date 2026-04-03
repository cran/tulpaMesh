## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  dev = "svglite",
  fig.ext = "svg"
)
library(tulpaMesh)

## ----sphere-basic-------------------------------------------------------------
globe <- tulpa_mesh_sphere(subdivisions = 3)
globe

## ----sphere-radius------------------------------------------------------------
radii <- sqrt(rowSums(globe$vertices^2))
range(radii)

## ----sphere-fem---------------------------------------------------------------
fem <- fem_matrices(globe, lumped = TRUE)

# Total surface area should approximate 4*pi for unit sphere
cat("Total area:", sum(fem$va), "\n")
cat("4*pi:      ", 4 * pi, "\n")
cat("Error:     ", abs(sum(fem$va) - 4 * pi) / (4 * pi) * 100, "%\n")

## ----sphere-proj--------------------------------------------------------------
obs <- cbind(lon = c(0, 90, -45, 170), lat = c(0, 45, -30, 60))
fem <- fem_matrices(globe, obs_coords = obs)
dim(fem$A)
range(Matrix::rowSums(fem$A))  # row sums = 1

## ----sphere-earth-------------------------------------------------------------
earth <- tulpa_mesh_sphere(subdivisions = 3, radius = 6371)
sqrt(sum(earth$vertices[1, ]^2))  # 6371 km

## ----euler--------------------------------------------------------------------
globe$n_vertices - globe$n_edges + globe$n_triangles

## ----mesh1d-------------------------------------------------------------------
# Monthly observations over 5 years
times <- seq(2020, 2025, by = 1/12)
m1d <- tulpa_mesh_1d(times)
m1d

## ----mesh1d-range-------------------------------------------------------------
range(times)       # data range
range(m1d$knots)   # mesh range (extended)

## ----mesh1d-fem---------------------------------------------------------------
# No extension for cleaner demonstration
m <- tulpa_mesh_1d(seq(0, 1, by = 0.1), n_extend = 0)

# Symmetric, positive definite mass
Matrix::isSymmetric(m$C)
all(Matrix::diag(m$C) > 0)

# Stiffness row sums = 0
max(abs(Matrix::rowSums(m$G)))

# Total mass = domain length
sum(m$C)

## ----mesh1d-irregular---------------------------------------------------------
# Denser sampling in summer, sparser in winter
times_irr <- c(
  seq(2020, 2020.25, by = 1/52),      # weekly Jan-Mar
  seq(2020.25, 2020.75, by = 1/365),   # daily Apr-Sep
  seq(2020.75, 2021, by = 1/52)        # weekly Oct-Dec
)
m_irr <- tulpa_mesh_1d(times_irr, n_extend = 0)
cat(m_irr$n, "knots from", length(times_irr), "unique time points\n")

## ----graph--------------------------------------------------------------------
# Simple river network: main channel + tributary
edges <- list(
  cbind(x = seq(0, 10, by = 0.5), y = rep(0, 21)),     # main channel
  cbind(x = c(5, 5, 5, 5), y = c(0, 2, 4, 6))          # tributary
)

g <- tulpa_mesh_graph(edges)
g

## ----graph-junctions----------------------------------------------------------
cat("Junctions (degree > 2):", sum(g$degree > 2), "\n")
cat("Endpoints (degree = 1):", sum(g$degree == 1), "\n")

## ----graph-fem----------------------------------------------------------------
Matrix::isSymmetric(g$C)
max(abs(Matrix::rowSums(g$G)))  # row sums ~ 0

## ----graph-refine-------------------------------------------------------------
g_fine <- tulpa_mesh_graph(edges, max_edge = 0.3)
cat("Coarse:", g$n_vertices, "vertices\n")
cat("Fine:  ", g_fine$n_vertices, "vertices\n")

## ----graph-sf, eval = requireNamespace("sf", quietly = TRUE)------------------
library(sf)
line1 <- st_linestring(cbind(c(0, 5, 10), c(0, 3, 0)))
line2 <- st_linestring(cbind(c(5, 5), c(3, 8)))
g_sf <- tulpa_mesh_graph(st_sfc(line1, line2), max_edge = 1)
g_sf

## ----nonstationary------------------------------------------------------------
set.seed(42)
mesh <- tulpa_mesh(cbind(runif(50), runif(50)), max_edge = 0.15)
n <- mesh$n_vertices

# Range decreases from left to right
kappa <- sqrt(8) / (2 - mesh$vertices[, 1])  # shorter range on the right
tau <- rep(1, n)

ns <- fem_matrices_nonstationary(mesh, kappa, tau)
names(ns)

## ----ns-constant--------------------------------------------------------------
ns_const <- fem_matrices_nonstationary(mesh, rep(2, n), rep(3, n))
fem <- fem_matrices(mesh)

max(abs(ns_const$Ck - 4 * fem$C))  # kappa^2 = 4
max(abs(ns_const$Ct - 9 * fem$C))  # tau^2 = 9

## ----p2-----------------------------------------------------------------------
set.seed(42)
mesh <- tulpa_mesh(cbind(runif(30), runif(30)))
p2 <- fem_matrices_p2(mesh)

cat("P1 nodes:", mesh$n_vertices, "\n")
cat("P2 nodes:", p2$n_mesh, "(", p2$n_vertices, "vertices +",
    p2$n_midpoints, "midpoints)\n")

## ----p2-area------------------------------------------------------------------
fem_p1 <- fem_matrices(mesh)
cat("P1 total area:", sum(fem_p1$C), "\n")
cat("P2 total area:", sum(p2$C), "\n")

## ----parallel-----------------------------------------------------------------
set.seed(42)
mesh_large <- tulpa_mesh(cbind(runif(500), runif(500)), max_edge = 0.03)
cat(mesh_large$n_triangles, "triangles\n")

fem_seq <- fem_matrices(mesh_large, parallel = FALSE)
fem_par <- fem_matrices(mesh_large, parallel = TRUE)

# Results are identical
max(abs(fem_seq$C - fem_par$C))

