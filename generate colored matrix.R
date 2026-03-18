# colored matrices for the paper

source("plot_colored_matrix.R")

# use these colors for color-blind accessibility
col_plus  <- "#009E73"
col_minus <- "#D55E00"
col_maybe <- "#F0E442"

# define tau and eta
K <- 6
cells <- expand.grid(x = 1:K, y = 1:K)
diag_cells <- subset(cells, x == y)
ul_cells <- subset(cells, y > x)

cairo_pdf("figures/eta_tau.pdf", width = 4, height = 4)
  plot_colored_matrix(
    K = 6,
    sets = list(
      list(x = diag_cells$x, y = diag_cells$y, col = col_maybe),
      list(x = ul_cells$x, y = ul_cells$y, col = col_plus)
    ),
    main = NULL,
    xlab = expression(italic(Y)^{(0)}),
    ylab = expression(italic(Y)^{(1)}),
  )
dev.off()


# an example of \Delta_j
K <- 6
cells <- expand.grid(x = 1:K, y = 1:K)
pos.mass <- subset(cells, (x < 3 & y >= 3))
neg.mass <- subset(cells, (x >= 3 & y < 3))

cairo_pdf("figures/Delta3.pdf", width = 4, height = 4)
plot_colored_matrix(
  K = 6,
  sets = list(
    list(x = pos.mass$x, y = pos.mass$y, col = col_plus),
    list(x = neg.mass$x, y = neg.mass$y, col = col_minus)
  ),
  main = NULL,
  xlab = expression(italic(Y)^{(0)}),
  ylab = expression(italic(Y)^{(1)}),
)
dev.off()


# upper bound for tau
K <- 6
cells <- expand.grid(x = 1:K, y = 1:K)
pos.mass <- cells 
neg.mass <- subset(cells, (x >= 3 & y < 3))
mass.2 <- subset(cells, (x < 3 & y >= 3))

cairo_pdf("figures/upper_tau.pdf", width = 4, height = 4)
plot_colored_matrix(
  K = 6,
  sets = list(
    list(x = cells$x, y = cells$y, col = col_plus),
    list(x = neg.mass$x, y = neg.mass$y, col = "white")
  ),
  marks = list(x = mass.2$x, y = mass.2$y),
  main = NULL,
  xlab = expression(italic(Y)^{(0)}),
  ylab = expression(italic(Y)^{(1)}),
)
dev.off()


# lower bound for tau
K <- 6
cells <- expand.grid(x = 1:K, y = 1:K)
pos.mass <- subset(cells, (x < 3 & y >= 3))
neg.mass <- subset(cells, (x >= 3 & y < 3))
add.mass <- subset(cells, (x == 3 & y >= 3 ))
no.mass <- subset(cells, (x == 3 & y < 3 ))

cairo_pdf("figures/lower_tau.pdf", width = 4, height = 4)
plot_colored_matrix(
  K = 6,
  sets = list(
    list(x = pos.mass$x, y = pos.mass$y, col = col_plus),
    list(x = neg.mass$x, y = neg.mass$y, col = col_minus),
    list(x = add.mass$x, y = add.mass$y, col = col_plus),
    list(x = no.mass$x, y = no.mass$y, col = "white")
  ),
  main = NULL,
  xlab = expression(italic(Y)^{(0)}),
  ylab = expression(italic(Y)^{(1)}),
)
dev.off()


# lower bounds - all possible j's

cairo_pdf("figures/lower_tau_all.pdf", width = 6, height = 9)
par(mfrow = c(3, 2))
for (j in 1:6) {
K <- 6
cells <- expand.grid(x = 1:K, y = 1:K)
pos.mass <- subset(cells, (x < j & y >= j))
neg.mass <- subset(cells, (x >= j & y < j))
add.mass <- subset(cells, (x == j & y >= j ))
no.mass <- subset(cells, (x == j & y < j ))

plot_colored_matrix(
  K = 6,
  sets = list(
    list(x = pos.mass$x, y = pos.mass$y, col = col_plus),
    list(x = neg.mass$x, y = neg.mass$y, col = col_minus),
    list(x = add.mass$x, y = add.mass$y, col = col_plus),
    list(x = no.mass$x, y = no.mass$y, col = "white")
  ),
  main = bquote(italic(j) == .(j)),
  xlab = expression(italic(Y)^{(0)}),
  ylab = expression(italic(Y)^{(1)}),
)
}
dev.off()
par(mfrow = c(1, 1))


# upper bounds - all possible j's
cairo_pdf("figures/upper_tau_all.pdf", width = 6, height = 9)
par(mfrow = c(3, 2))
for (j in 1:6) {
  K <- 6
  cells <- expand.grid(x = 1:K, y = 1:K)
  pos.mass <- cells 
  neg.mass <- subset(cells, (x >= j & y < j))
  mass.2 <- subset(cells, (x < j & y >= j))
  
  if (j==1) {
  plot_colored_matrix(
    K = 6,
    sets = list(
      list(x = cells$x, y = cells$y, col = col_plus),
      list(x = neg.mass$x, y = neg.mass$y, col = "white")
    ),
    main = bquote(italic(j) == .(j)),
    xlab = expression(italic(Y)^{(0)}),
    ylab = expression(italic(Y)^{(1)}),
  )
  } else {
  
  plot_colored_matrix(
    K = 6,
    sets = list(
      list(x = cells$x, y = cells$y, col = col_plus),
      list(x = neg.mass$x, y = neg.mass$y, col = "white")
    ),
    marks = list(x = mass.2$x, y = mass.2$y),
    main = bquote(italic(j) == .(j)),
    xlab = expression(italic(Y)^{(0)}),
    ylab = expression(italic(Y)^{(1)}),
  )
  }
}
dev.off()
par(mfrow = c(1, 1))


# example of matrix attaining the lower bound for tau
K <- 6
cells <- expand.grid(x = 1:K, y = 1:K)
cond <- (cells$x>=4 & cells$y<=2) |
  (cells$y>=cells$x & cells$x>=4 & cells$y>=4) |
  (cells$y>=cells$x & cells$x<=2 & cells$y<=2)
mass.0 <- cells[cond, ]
mass.1 <- cells[!cond, ]

cairo_pdf("figures/lower_tau_example.pdf", width = 4, height = 4)
plot_colored_matrix(
  K = 6,
  #sets = list(list(x = mass.1$x, y = mass.1$y, col = col_plus)),
  marks = mass.0,
  mark_label = "0",
  main = NULL,
  xlab = expression(italic(Y)^{(0)}),
  ylab = expression(italic(Y)^{(1)}),
)
dev.off()


# example of matrix attaining the upper bound for tau
K <- 6
cells <- expand.grid(x = 1:K, y = 1:K)
cond <- (cells$x<=2 & cells$y>=3) |
  (cells$x==2 & cells$y==1) |
  (cells$y<cells$x & cells$x>=4 & cells$y>=3)
mass.0 <- cells[cond, ]
mass.1 <- cells[!cond, ]

cairo_pdf("figures/upper_tau_example.pdf", width = 4, height = 4)
plot_colored_matrix(
  K = 6,
  #sets = list(list(x = mass.1$x, y = mass.1$y, col = col_plus)),
  marks = mass.0,
  mark_label = "0",
  main = NULL,
  xlab = expression(italic(Y)^{(0)}),
  ylab = expression(italic(Y)^{(1)}),
)
dev.off()


