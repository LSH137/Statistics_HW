# ==== ==== ==== 1 ==== ===== =====
cat("-- problem 1 --\n")
# --------------- (a) ----------------
TRIN <- function(n, p1, p2, x, y){
  z = n - x - y
  p3 = 1 - p1 - p2
  n_fac = factorial(n)
  x_fac = factorial(x)
  y_fac = factorial(y)
  z_fac = factorial(z)
  return((n_fac * p1^x * p2^y * p3^z) / (x_fac * y_fac * z_fac))
}

n = 32
p1 = 1/8
p2 = 1/4
p3 = 1 - p1 - p2

answer = 0
for(x in 0:n){
  for(y in 0:(n-x)){
      z = n - x - y
      #cat("x = ", x, "y = ", y, "z = ", z, "\n")
      answer = answer + TRIN(n, p1, p2, x, y) * ((y-2*x)*(y+2*z))
  }
}
cat("(a) : ", answer, "\n")

# --------------- (b) ----------------
marginal_X <- function(x){
  answer = 0
  for(y in 0:(n-x)){
    answer = answer + TRIN(n, p1, p2, x, y)
  }
  
  return(answer)
}

marginal_Y <- function(y){
  answer = 0
  for(x in 0:(n-y)){
    answer = answer + TRIN(n, p1, p2, x, y)
  }
  return(answer)
}

ux2 = 0
for(x in 0:n){
  ux2 = ux2 + (x^2) * marginal_X(x)
}

ux4 = 0
for(x in 0:n){
  ux4 = ux4 + (x^4) * marginal_X(x)
}

uy2 = 0
for(y in 0:n){
  uy2 = uy2 + (y^2) * marginal_Y(y)
}

uy4 = 0
for(y in 0:n){
  uy4 = uy4 + (y^4) * marginal_Y(y)
}

cov_X2Y2 = 0
for(x in 0:n){
  for(y in 0:(n-x)){
    cov_X2Y2 = cov_X2Y2 + TRIN(n, p1, p2, x, y) * ((x^2-ux2)*(y^2-uy2))
  }
}

var_x2 = ux4 - ux2^2
var_y2 = uy4 - uy2^2

sigma_x2 = sqrt(var_x2)
sigma_y2 = sqrt(var_y2)

corr_X2Y2 = cov_X2Y2 / (sigma_x2 * sigma_y2)
cat("(b) : ", corr_X2Y2, "\n")

# ==== ==== ==== ==== ==== ==== 2 ==== ==== ==== ==== ===== =====
cat("-- problem 2 --\n")
n_sample = 5000
max_x = 15
max_y = 15
dx = max_x / n_sample
dy = max_y / n_sample

X <- seq(0, max_x, dx)
Y <- seq(0, max_y, dy)

JointPDF <- function(x, y){
  return(4 * x^3 * exp((-1)*x*(y+2)))
}

GammaDistriution <- function(a, t, x){
  return(1/(gamma(a) * t^a) * x^(a-1)*exp(-1*x / t))
} 

marginal_X_2 <- function(x){
  return(GammaDistriution(3, 1/2, x))
}


marginal_Y_2 <- function(y){
  result = 0
  for(x in X){
    result = result + JointPDF(x, y) * dx
  }
  return(result)
}

conditional_Y_2 <- function(y, x){
  return(GammaDistriution(1, 1/x, y))
}

p_lst_x = marginal_X_2(X) * dx
plot(X, p_lst_x, type="l") # maximum x is reasonably large 
cat("sum of probability = ", sum(p_lst_x), "\n") # maximum x is reasonably large 

sum_x = 0
sum_y = 0
Ex = 0
Ey = 0
sample_result_list_x <- vector()
sample_result_list_y <- vector()
for(i in 0:n_sample){
  sample_result_x = sample(X, 1, replace = T, prob = p_lst_x)
  
  p_lst_y = conditional_Y_2(Y, sample_result_x) * dy
  sample_result_y = sample(Y, 1, replace = T, prob = p_lst_y)
  
  sum_x = sum_x + sample_result_x
  sum_y = sum_y + sample_result_y
  
  sample_result_list_x <- c(sample_result_list_x, sample_result_x)
  sample_result_list_y <- c(sample_result_list_y, sample_result_y)
}

average_x = sum_x / n_sample
average_y = sum_y / n_sample

diff2_sum_x = 0
diff2_sum_y = 0
for(i in 1:n_sample){
  diff2_sum_x = diff2_sum_x + (sample_result_list_x[i] - average_x)^2
  diff2_sum_y = diff2_sum_y + (sample_result_list_y[i] - average_y)^2
}

var_x = diff2_sum_x / (n_sample - 1)
var_y = diff2_sum_y / (n_sample - 1)

cov_XY = 0
sample_list <- list()
for(i in 1:n_sample){
  x = sample_result_list_x[i]
  y = sample_result_list_y[i]
  cov_XY = cov_XY + (1/(n_sample - 1)) * (x - average_x)*(y - average_y)
}
corr = cov_XY / (sqrt(var_x) * sqrt(var_y))

cat("cov_xy = ", cov_XY, "\n")
cat("- - - -")



# print result
cat("average of x = ", average_x, "\n")
cat("average of y = ", average_y, "\n")
cat("variance of x = ", var_x, "\n")
cat("variance of y = ", var_y, "\n")
cat("correlation of x and y = ", corr, "\n")
