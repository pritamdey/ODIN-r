external_file_path = "E:/Outlier Detection for Multi-Network Data/01. Code/R Code/github_upload/"
add_code_path = ""
add_data_path = "data"

code_path = paste0(external_file_path, add_code_path)
setwd(code_path)

data_path = paste0(external_file_path, add_data_path)

library(readr)

source(paste0(code_path, "logistic_fit.R"))
source(paste0(code_path, "helper_functions.R"))

#X
xmeta = read_csv(paste0(data_path,"LobeHemiMatrix.csv"))
X = make_X(xmeta)

#A
A = as.matrix(read_csv(paste0(data_path,"A_example.csv"), col_names = FALSE))
N=ncol(A)

#Algorithm Parameters
lam=0.001
tol=5e-6

#Model Fitting and influence calculation
clean_fit = mm_logistic(A,X,lam,tol,maxiter=10000,print_out = T)
influence_clean = influence_measure(clean_fit,A,X,lam)


#Outlier Detection

# First influence measure threshhold
med = median(influence_clean[[1]])
quartiles = quantile(influence_clean[[1]], c(0.25, 0.75))
iqr = quartiles[2] - quartiles[1]
threshold1 = med + 1.5 * iqr

# Second influence measure threshold
threshold2 = qchisq(0.999, 23)

# Indices of the outliers
out1 = which(influence_clean[[1]] > threshold1)
out2 = which(influence_clean[[2]] > threshold2)
out = union(out1, out2)

# Inliers
ins = setdiff(1:N, out)

#Influence Plot
plot(influence_clean[[1]],xlab = "Observation No.",ylab="IM(i)",main="")
abline(h=threshold1,col="red")

plot(influence_clean[[2]],xlab = "Observation No.",ylab="d_beta(i)",main="")
abline(h=threshold2,col="red")
