#!/mnt/software/R/3.5.1/bin/Rscript

#SBATCH --mail-user=dariush.mollet@students.unibe.ch
#SBATCH --mail-type=ALL
#SBATCH --job-name="lasso"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=72:00:00
#SBATCH --mem=48G
#SBATCH --output=good.out
#SBATCH --error=bad.error

library(glmnet)
library(plotrix)

options(warn=1)

#get files
files = sort(list.files(path="stats_all", pattern="*.stats_all", full.names=TRUE))
#files = sort(list.files(path="stats", pattern="*.stats", full.names=TRUE))
all.list = lapply(files, read.csv, sep=" ")

r_rates = read.csv("repl_rates.csv", sep=",", header=FALSE, col.names=c("strain", "rate"))
r_rates = r_rates[order(r_rates$strain),]

#all in one nice matrix
feats = c("rank", "dir", "occ")
stacked.list = lapply(all.list, stack, select=feats)
raw_data = do.call("cbind", stacked.list)

data = raw_data[, -which(names(raw_data) %in% c("ind"))]

data = scale(data)

data = t(as.matrix(data))

#features
genes = all.list[[1]]$gene
features = c()
for (i in 1:length(feats)){
    for (j in 1:length(genes)){
        features = rbind(features, paste(genes[j],"_", feats[i], sep=""))
    }
}

colnames(data) = features
rownames(data) = r_rates[,1]

#lasso multiple runs different alphas:

alpha_seq = seq(0,1,0.01)
results = c()
for (i in alpha_seq){
    
    print(i)
    
    lasso_reg = cv.glmnet(data, r_rates[,2], type.measure="mse", alpha=i, lambda=1:500, nfolds=length(r_rates[,1]))

    min_lambda = lasso_reg$lambda.min

    lasso_pred = predict(lasso_reg, s = min_lambda, newx = data) # Use best lambda to predict test data
    MSE = mean((lasso_pred - r_rates[,2])^2)

    #results
    lasso_coef = coef(lasso_reg, s = min_lambda)
    n_params = length(which(lasso_coef != 0))

    results = rbind(results,cbind(i,MSE,round(n_params)))
}
colnames(results) = c("alpha", "MSE", "number of params")

optimal_alpha = alpha_seq[which(results[,2] == min(results[-1,2]))]
optimal_alpha

png("results/lasso.png")
twoord.plot(alpha_seq, results[,2], alpha_seq, results[,3], main="Lasso optimization", lcol="black", rcol="red", xlab="alpha", ylab="MSE", rylab="number of parameters", lpch=20, rpch=20, type="o")

#abline(v=optimal_alpha, col="blue", lty=5, lwd=1)                                                                                                              
#lasso regression - single run
lasso_reg = cv.glmnet(data, r_rates[,2], type.measure="mse", alpha=optimal_alpha, lambda=1:500, nfolds=length(r_rates[,1]))

min_lambda = lasso_reg$lambda.min

lasso_coef = coef(lasso_reg, s = min_lambda)

lasso_pred = predict(lasso_reg, s = min_lambda, newx = data) # Use best lambda to predict test data
MSE = mean((lasso_pred - r_rates[,2])^2)

MSE
lasso_pred
r_rates
rownames(lasso_coef)[which(lasso_coef != 0)]


print("ok done")
stop()

# plot
pdf("results/lasso.pdf")
plot(lasso_reg, xvar = "lambda")
lasso_reg = glmnet(data, r_rates[,2])
plot(lasso_reg, xvar = "lambda")

