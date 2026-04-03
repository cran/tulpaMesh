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

## ----basic-mesh---------------------------------------------------------------
set.seed(42)
coords <- cbind(x = runif(100), y = runif(100))
mesh <- tulpa_mesh(coords)
mesh

## ----plot-basic---------------------------------------------------------------
plot(mesh, vertex_col = "steelblue", main = "Basic mesh")

## ----refined-mesh-------------------------------------------------------------
mesh_fine <- tulpa_mesh(coords, max_edge = 0.08)
mesh_fine
plot(mesh_fine, main = "Refined mesh (max_edge = 0.08)")

## ----fem----------------------------------------------------------------------
fem <- fem_matrices(mesh_fine, obs_coords = coords)
dim(fem$C)
dim(fem$A)

# Verify key properties
all(Matrix::diag(fem$C) > 0)        # positive diagonal
max(abs(Matrix::rowSums(fem$G)))     # row sums ~ 0
range(Matrix::rowSums(fem$A))        # row sums = 1

## ----lumped-------------------------------------------------------------------
fem_l <- fem_matrices(mesh_fine, obs_coords = coords, lumped = TRUE)
Matrix::isDiagonal(fem_l$C0)

## ----formula------------------------------------------------------------------
df <- data.frame(lon = runif(50), lat = runif(50), y = rnorm(50))
mesh_f <- tulpa_mesh(~ lon + lat, data = df)
mesh_f

## ----quality------------------------------------------------------------------
mesh_summary(mesh_fine)

## ----quality-plot-------------------------------------------------------------
plot(mesh_fine, color = "quality", main = "Colored by minimum angle")

## ----ruppert------------------------------------------------------------------
mesh_r <- tulpa_mesh(coords, min_angle = 25, max_edge = 0.15)
mesh_summary(mesh_r)

