# Exercice 1


# Question 2

# Fonction qui simule selon g
rgen_g <- function(n) {
  g1 <- runif(n, min = (-pi / 2), max = (pi / 2))
  g2 <- runif(n, min = -1, max = 1)
  return(cbind(g1, g2))
}

# Fonction qui calcule psi(x,y)
psi <- function(x, y) {
  result1 <- abs(sin((2 / pi) * (x**2) - (pi / 4))) + (4 * cos(x)**2) + y**4
  result2 <- exp(-2 * (x + abs(y)))
  return(result1 * result2)
}

# Variables utilisees dans la suite de l'exercice
n <- 10000
m <- (sqrt(2) * pi) * exp(pi)
c_g <- Sys.time() # Temps de calcul pour question 6
g <- rgen_g(n)
c_g <- Sys.time() - c_g
c_u <- Sys.time() # Temps de calcul pour question 6
uni <- runif(n, 0, 1)
ratio <- psi(g[, 1], g[, 2]) / (m / (2 * pi))
c_u <- Sys.time() - c_u
val_ratio <- ratio

# Fonction qui simule selon f par la methode du rejet
rgen_f <- function(g_1, g_2, uni, ratio) {
  n_l <- 0
  while (uni > ratio) {
    uni <- runif(1, 0, 1)
    g <- rgen_g(1)
    g_1 <- g[1]
    g_2 <- g[2]
    ratio <- psi(g_1, g_2) / (m / (2 * pi))
    val_ratio <- c(val_ratio, ratio)
    n_l <- n_l + 1
  }
  assign("val_ratio", val_ratio, envir = .GlobalEnv)
  return(cbind(g_1, g_2, n_l))
}


# Question 3

c_f <- Sys.time() # Temps de calcul pour question 6
z <- Vectorize(rgen_f)(g[, 1], g[, 2], uni, ratio)
c_f <- Sys.time() - c_f
z <- t(z) # Pour une visualisation plus agreable
n_l <- sum(z[, 3]) / n
n_t <- length(val_ratio)


# Question 4b)

c_b_hat <- Sys.time() # Temps de calcul pour question 6
# Evaluation de bn chapeau
wk <- 2 * pi * psi(g[, 1], g[, 2])
wk_hzk <- rep(1, n)
bn_hat <- sum(wk_hzk) / sum(wk)
c_b_hat <- Sys.time() - c_b_hat
# Evaluation de l'intervalle de confiance asymptotique
sigma_hat_wk <- (sum((wk - sum(wk) / n)**2)) / (n - 1)
ecart <- ((qnorm(0.975) * sqrt(sigma_hat_wk) * (bn_hat**2)) / sqrt(n))
ic_95_a_inf <- (n / sum(wk)) - ecart
ic_95_a_sup <- (n / sum(wk)) + ecart


# Question 4c)

# Estimation du biais de bn chapeau par une methode de bootstrap
k <- 1000
bn_bootstrap <- matrix(sample(wk, n * k, replace = TRUE), ncol = n)
bn_bootstrap <- n / rowSums(bn_bootstrap)
estim_biais <- mean(bn_bootstrap) - bn_hat


# Question 5b)

c_a_hat <- Sys.time() # Temps de calcul pour question 6
# Estimation de an chapeau
h_xi <- 1 / (2 * pi * psi(z[, 1], z[, 2]))
an_hat <- (sum(h_xi)) / n
c_a_hat <- Sys.time() - c_a_hat
# Evaluation de l'intervalle de confiance asymptotique
sigma_hat_hx <- sqrt(sum((h_xi + sum(h_xi) / n)**2) / (n - 1))
ic2_95_a_inf <- an_hat - ((qnorm(0.975) * sigma_hat_hx) / sqrt(n))
ic2_95_a_sup <- an_hat + ((qnorm(0.975) * sigma_hat_hx) / sqrt(n))


# Question 6

# Estimation du cout de calcul des estimateurs
c_b <- c_g + c_b_hat
c_a <- c_f + c_a_hat + c_u
# Estimation des variances des estimateurs
sigma_hat_b <- var(bn_bootstrap)
sigma_hat_a <- (sigma_hat_hx^2) / n
# Efficacite relative de bn chapeau par rapport a an chapeau
r_b_a <- (as.numeric(c_b) * sigma_hat_b) / (as.numeric(c_a) * sigma_hat_a)


# Question 7b)

# Fonction qui a x et n associe fXn(x) chapeau
fn_hat_x <- function(x, n) {
  x_i <- matrix(runif(n * length(x)), ncol = n)
  # Et non simplement runif(n) afin que la fonction supporte un vecteur
  # en entree et qu'elle puisse ainsi nous permettre le tracer 
  # de sa courbe dans les lignes ci-dessous.
  cx <- abs(sin((2 * (x^2) / pi) + (pi / 4))) + 4 * (cos(x)^2)
  hx_x_i <- (cx + x_i^4) * exp(-2 * x_i)
  resultat <- 2 * an_hat * exp(-2 * x) * rowSums(hx_x_i) / n
  return(resultat)
}

titre <- "Comparaison de la densite marginale empirique de f a l'estimateur"
axex <- "x"
axey <- "f_X(x)"
t <- seq(-pi / 2, pi / 2, 0.1)
hist(z[, 1], freq = FALSE, ylim = c(0, 1), xlab = axex, ylab = axey, main = titre)
lines(t, fn_hat_x(t, n), col = "red")


# Question 10

# Liste des coordonnees (x,y) qui forment la densite empirique associee a z[,1]
fx_emp <- density(z[, 1], from = -1.6, to = 1.6, n = 321)
fx_emp <- cbind(fx_emp[["x"]], fx_emp[["y"]])

# Fonction qui a x associe fX(x)
# Calcule empiriquement avec l'echantillon z
densite_emp <- function(x) {
  resultat <- fx_emp[which(round(fx_emp[, 1], 2) == round(x, 2)), 2]
  return(cbind(resultat))
}

c_w_hat <- Sys.time() # Temps de calcul pour question 11
# Estimation de wn(-1) chapeau
x <- -1
wk_emp <- (sapply(z[, 1], densite_emp))
uk_x <- (psi(x, z[, 2]) * wk_emp) / psi(z[, 1], z[, 2])
w_hat_x <- sum(uk_x) / n
c_w_hat <- Sys.time() - c_w_hat
# Evaluation de l'intervalle de confiance asymptotique de fX(-1)
sigma_hat_uk <- var(uk_x)
ic_95_fx_inf <- w_hat_x - qnorm(0.975) * sqrt(sigma_hat_uk / n)
ic_95_fx_sup <- w_hat_x + qnorm(0.975) * sqrt(sigma_hat_uk / n)

# Question 11
c_fx_hat <- Sys.time() # Temps de calcul
# Estimation de fXn(-1) chapeau
x_i <- runif(n)
cmoins1 <- abs(sin((2 * (x^2) / pi) + (pi / 4))) + 4 * (cos(x)^2)
hmoins1_x <- (cmoins1 + x_i^4) * exp(-2 * x_i)
fn_hat_moins1 <- 2 * an_hat * exp(-2 * x) * sum(hmoins1_x) / n
c_fx_hat <- Sys.time() - c_fx_hat
# Estimation des variance de fXn(-1) chapeau et wn(-1) chapeau
uk_x_bootstrap <- matrix(sample(uk_x, n * k, replace = TRUE), ncol = n)
w_x_bootstrap <- rowSums(uk_x_bootstrap) / n
sigma_hat_fx <- var(hmoins1_x) * ((4 * (an_hat^2) * exp(4))) / n
sigma_hat_w <- var(w_x_bootstrap)
# Estimation du cout de calcul de estimateurs
c_w_moins1 <- as.numeric(c_w_hat + c_f + c_u)
c_fx_moins1 <- as.numeric(c_a + c_fx_hat)
# Efficacite relative de wn(-1) chapeau par rapport a fXn(-1) chapeau
r_fx_w_moins1 <- (c_w_moins1 * sigma_hat_w) / (c_fx_moins1 * sigma_hat_fx)


####################


# Exercice 2


# Question 1

# Fonction qui simule selon une loi normale multivariee
rmvnorm <- function(n, mu, sigma) {
  normale <- matrix(rnorm(length(mu) * n), nrow = length(mu))
  return(t(chol(sigma)) %*% normale + mu)
  # chol() donne la matrice triangulaire superieur de Cholesky
  # --> On transpose
}

# Variables utilisees dans la suite de l'exercice
mu <- c(0.1, 0, 0.1)
sigma <- matrix(c(0.047, 0, 0.0117, 0, 0.047, 0, 0.0117, 0, 0.047), nrow = 3)
n <- 10000
x <- rmvnorm(n, mu, sigma)
x <- t(x) # Pour une visualisation plus agreable


# Question 2b)

# Application de h a notre echantillon x
hx <- rowSums(exp(-x)) / 3
hx[hx > 3] <- 3
delta_barre <- mean(hx)
mse_barre <- var(hx) / n


# Question 3b)

# Application de h(A) a notre echantillon x
ax <- 2 * mu - t(x)
ax <- t(ax) # Pour une visualisation plus agreable
hax <- rowSums(exp(-ax)) / 3
hax[hax > 3] <- 3
delta_hat <- sum(hx + hax) / (2 * n)
mse_hat <- var((hx + hax) / 2) / n
r1 <- mse_hat / mse_barre
print(r1)


# Question 4 a)

# Test 1
h0x1 <- rowSums(x) / 3
rho1 <- cor(hx, h0x1)
# Test 2
h0x2 <- rowSums((x - colSums(x) / 3)^2) / 2
rho2 <- cor(hx, h0x2)
# Test 3
h0x3 <- rowSums(x^2) / 3
rho3 <- cor(hx, h0x3)

# Question 4 b)

m <- mean(mu)
# Estimation de b* avec la methode de burn-in
l <- 1000
x_tilde <- t(rmvnorm(l, mu, sigma))
h0x_tilde <- rowSums(x_tilde) / 3
hx_tilde <- rowSums(exp(-x_tilde)) / 3
hx_tilde[hx_tilde > 3] <- 3
b_star_hat <- cov(h0x_tilde, hx_tilde) / var(h0x_tilde)
delta_hat_b_star <- mean(hx - b_star_hat * (h0x1 - m))
mse_hat_b <- var(hx - b_star_hat * (h0x1 - m)) / n
r2 <- mse_hat_b / mse_barre
print(r2)


####################


# Exercice 3


# Question 1

# Fonction qui simule S pour un Y donne
simul_s <- function(y) {
  return(sum(log(rgamma(y, m, theta) + 1)))
}

# Variables utilisees dans la suite de l'exercice
n <- 10000
p <- 0.2
m <- 2
theta <- 2

c_mc <- Sys.time() # Temps de calcul pour question 2b)
geo <- rgeom(n, p) + 1 # Sur R, rgeom est a valeur dans N et non N*
tirage_s <- sapply(geo, simul_s) # echantillon de n tirages de S
c_mc <- Sys.time() - c_mc
delta_mc <- (sum(tirage_s)) / n

# estimateur sans biais --> mse(delta_mc)=var(delta_mc)
mse_mc <- var(tirage_s) / n

# Histogramme pour se donner une idee du tirage obtenu
bords <- seq(0, round(max(tirage_s)) + 1)
titre <- "Histogramme de 10 000 tirages de S"
axex <- "Tirages de S"
axey <- "Frequence d'apparition"
hist(tirage_s, breaks = bords, freq = F, xlab = axex, ylab = axey, main = titre)


# Question 2b)

c_strat <- Sys.time()
vect_proba <- dgeom(0:13, p) # Car la geometrique commence en 0 sur R
vect_proba[15] <- 1 - sum(vect_proba)
allocation <- round(n * vect_proba)
allocation[10] <- allocation[10] + 1 # Afin d'avoir "somme des allocations = n"
# Simulation de Y sachant qu'elle appartient a la derniere strate
unif <- runif(allocation[15])
y_15 <- qgeom((pgeom(13, p) + unif * (1 - pgeom(13, p))), p) + 1
# On evalue F en 13 = 14-1 car la geometrique commence a 0 sur R.
# On ajoute 1 pour la meme raison.
y_2 <- c(rep(1:14, allocation[1:14]), y_15)
tirage_s_2 <- sapply(y_2, simul_s)
c_strat <- Sys.time() - c_strat
delta_strat <- (sum(tirage_s_2)) / n

# estimateur sans biais --> mse(delta_strat) = var(delta_strat)
mse_strat <- var(vect_proba) / n
# Efficacite relative de la methode
eff_relative <- (as.numeric(c_mc) * mse_mc) / (as.numeric(c_strat) * mse_strat)