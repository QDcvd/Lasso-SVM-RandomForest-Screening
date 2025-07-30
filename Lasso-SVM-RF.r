# -*- coding: utf-8 -*-
#
#  Created on Tue Apr 14 14:53:46 2025
#
#  @author: 真弓快车
#
#  该脚本的作用是联合使用LASSO-SVM-RF进行核心基因筛选、重要性排序与诊断模型构建
#
library(data.table)
library(dplyr)
setwd('./example')
data <- read.csv('crucialGene_1.csv', header = T)

# LASSO出图
dir.create('Lasso')
library(glmnet)
# 读取输入数据
# expression_data <- data[, -length(colnames(data))]
expression_data <- read.csv('crucialGene_1.csv', header = TRUE, check.names = FALSE, row.names = 1)
#模型输入
expression_data_t <- t(expression_data)
feature_matrix <- as.matrix(expression_data_t)[, 1:15]

# feature_matrix <- t(as.matrix(expression_data)[, 2:47])

# sample_labels <- data$Group
sample_labels <- gsub("(.*)\\_(.*)", "\\2", rownames(expression_data_t))[1:47]
table(sample_labels)
# 构建LASSO回归模型
lasso_model <- glmnet(feature_matrix, sample_labels, family = 'binomial', alpha = 1)
print(lasso_model$df)

# 图1 LASSO系数路径图
pdf(file = "./Lasso/1_lasso_coefficient_paths.pdf", width=6, height=5.5)
plot(lasso_model)
# 展示系数随log(lambda)变化路径
# 纵轴: 特征系数值
# 每条曲线代表一个特征基因在不同惩罚强度下的系数变化
# 顶部数字表示当前lambda下的非零特征数量
dev.off()
# 交叉验证确定最优模型
set.seed(114514)
# 10折交叉验证
cv_model <- cv.glmnet(feature_matrix, sample_labels, family="binomial", alpha=1, nfolds=10)
# 图2 交叉验证误差曲线
pdf("./Lasso/2_cross_validation_error.pdf", width=6, height=5)
plot(cv_model)  # 展示交叉验证误差随惩罚强度变化
# 横轴: log(lambda) - 惩罚强度
# 纵轴: 二项偏差（分类误差）
# 红点: 不同lambda下平均误差
# 误差线: ±1个标准差
# 垂直线:
#   - 左: lambda.min (最小误差对应的lambda)
#   - 右: lambda.1se (最小误差1个标准误内的最简模型)
dev.off()

# 提取最优模型系数
optimal_lambda <- cv_model$lambda.min # 1个标准误规则下的lambda
model_coef <- coef(cv_model, s=optimal_lambda) # 提取最优lambda下的系数
# model_coef <- coef(cv_model, s=0.05)
# 筛选重要的生物标志物
selected_features <- rownames(model_coef)[which(model_coef != 0)]
selected_features <- selected_features[-1] # 移除截距项(第一个元素)
# 准备特征重要性数据
feature_importance <- as.matrix(model_coef)[selected_features, ]
sorted_idx <- order(abs(feature_importance), decreasing = FALSE)
sorted_importance <- feature_importance[sorted_idx] # 排序后的系数

# 图3 带方向的特征重要性
pdf("./Lasso/3_feature_importance_directional.pdf", width=7, height=5)
barplot(sorted_importance, 
        horiz=TRUE, 
        las=1,
        col=ifelse(sorted_importance > 0,"tomato","steelblue"),
        xlab="Coefficient Value",
        cex.names=0.7,
        main="Biomarker Impact Direction")

dev.off()

# 图4 特征绝对重要性
pdf("./Lasso/4_feature_importance_absolute.pdf", width=7, height=5)
barplot(abs(sorted_importance), horiz=TRUE, las=1,
        col = 'tomato',
        xlab = "Absolute Coefficient Value",
        main = "Gene Importance Ranking",
        cex.names=0.7
        )

dev.off()

# 保存结果
write.table(data.frame(Biomarker=selected_features, Coefficient=feature_importance), "./Lasso/Lasso_selected_genes.txt", sep="\t", quote=FALSE, row.names=FALSE)


# SVM-RFE 出图
# SVM
dir.create("SVM-RFE")
# 设置工作目录
# 加载必要的R包
library(limma)
library(e1071)
library(caret)
library(ggplot2)
library(pROC)
# 加载自定义SVM-RFE函数
source("./SVM-RFE/msvmRFE.R")

### 数据准备阶段 ---------------------------------------------------
# 读取基因表达数据文件
#expression_data <- data[,-length(colnames(data))]
# 获取样本分组信息
#sample_groups <- data$Group
expression_data <- read.csv('crucialGene_1.csv', header = TRUE, check.names = FALSE, row.names = 1)
#模型输入
expression_data_t <- t(expression_data)
#feature_matrix <- as.matrix(expression_data_t)[, 1:16]

# feature_matrix <- t(as.matrix(expression_data)[, 2:47])

# sample_labels <- data$Group
sample_labels <- gsub("(.*)\\_(.*)", "\\2", rownames(expression_data_t))[1:47]
table(sample_labels)
# 创建分析数据集
analysis_data <- data.frame(
                    Group = factor(sample_labels),
                    expression_data_t)
### SVM-RFE特征选择 -----------------------------------------------
# 执行初始SVM-RFE
svmRFE(analysis_data, k = 10, halve.above = 50)
# 设置5折交叉验证
n_folds <- 10
num_samples <- nrow(analysis_data)
# 创建分层交叉验证折叠
set.seed(123)
fold_index <- rep(1:n_folds, length.out = num_samples)
fold_index <- sample(fold_index)
fold_list <- lapply(1:n_folds, function(x) which(fold_index == x))
# 执行交叉验证
cv_results <- lapply(fold_list, function(test_indices) {svmRFE.wrap(test.fold = test_indices,
    X = analysis_data,
    k = 10,
    halve.above = 50)})

### 特征排序与结果输出 --------------------------------------------
# 获取特征重要性排序
feature_ranking <- WriteFeatures(cv_results, analysis_data, save = FALSE)
# 保存特征排序结果
write.table(feature_ranking, file = "SVM-RFE/SVM_Feature_Ranking.txt", sep ="\t",
            quote = FALSE, 
            row.names = FALSE)
cat("SVM-RFE selected", nrow(feature_ranking), "features\n")


### 交叉验证性能评估 --------------------------------------------
# 确定评估的特征数量范围
max_features <- ifelse(ncol(analysis_data) - 1 > 30, 30, ncol(analysis_data) - 1)
# 每完成10%输出一次进度
performance_list <- lapply(1:max_features, function(n_features) {
    if(n_features %% round(max_features/10) == 0) {
        cat(sprintf(">> 进度: %.0f%% [%s]\n",
        100*n_features/max_features,
        Sys.time()))
        }
        FeatSweep.wrap(n_features, cv_results, analysis_data)
        })
# 提取误差率
error_rates <- sapply(performance_list, function(x) {
    if(is.null(x)) return(NA) else x$error})

# 计算准确率
accuracy_rates <- 1 - error_rates
# 计算随机基线
baseline_error <- min(prop.table(table(analysis_data$Group)))
baseline_accuracy <- 1 - baseline_error
### 性能可视化 ----------------------------------------------------
# 图1: 交叉验证误差曲线
dev.new()
pdf("SVM-RFE/1_SVM_Error_Rate.pdf", width = 7, height = 5)
PlotErrors(error_rates, no.info = baseline_error,
            xlab = "Number of Features", ylab = "Classification Error Rate")
dev.off()
# 图2: 交叉验证准确率曲线
pdf("SVM-RFE/2_SVM_Accuracy.pdf", width = 7, height = 5)
Plotaccuracy(accuracy_rates, no.info = baseline_accuracy,
            xlab = "Number of Features", ylab = "Classification Accuracy")
dev.off()
# 图3: 特征数量与性能关系图
pdf("SVM-RFE/3_SVM_Performance_Comparison.pdf", width = 8, height = 6)
plot(1:max_features, accuracy_rates, type = "b", col = "dodgerblue", lwd = 1,
        xlab = "Number of Features", ylab = "Performance Metric",
        main = "SVM-RFE Feature Optimization",
        ylim = c(0, 1))
lines(1:max_features, error_rates, type = "b", col = "firebrick", lwd = 2)
abline(h = baseline_accuracy, col = "gray50", lty = 2)
abline(h = baseline_error, col = "gray50", lty = 2)
legend("topright", legend = c("Accuracy", "Error Rate", "Random Baseline"),
        col = c("dodgerblue", "firebrick", "gray50"), lwd = 1, lty = c(1, 1, 2), cex = 0.5)
dev.off()
# 图4: 特征重要性条形图
top_features <- head(feature_ranking, 15)
pdf("SVM-RFE/4_SVM_Feature_Importance.pdf", width = 10, height = 6)
ggplot(top_features, aes(x = reorder(FeatureName, AvgRank),
                        y = AvgRank, fill = AvgRank)) + 
                        geom_bar(stat = "identity") + 
                        scale_fill_gradient(low = "gray80", high = "darkred") +
                        coord_flip() + 
                        labs(title = "Top Features Selected by SVM-RFE",
                            x = "Gene Name", y = "Average Rank (Lower = More Important)") +
                            theme_minimal(base_size = 12) +
                            theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()
# 确定最佳特征数量
optimal_feature_count <- which.min(error_rates)
selected_genes <- feature_ranking[1:optimal_feature_count, "FeatureName", drop = FALSE]
# 保存筛选的特征基因
write.table(selected_genes, file = "SVM-RFE/SVM_Selected_Genes.txt", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)


# 随机森林出图
### RF
dir.create("RF")
# 设置工作目录
# 加载必要的包
library(randomForest)
library(pROC)
library(ggplot2)
library(caret)
### 数据准备阶段 ---------------------------------------------------
# 读取基因表达数据
#expression_data <- data[,-length(colnames(data))]
# 提取样本标签
#sample_labels <- data$Group
# 创建分析数据集


# expression_data <- data[, -length(colnames(data))]
expression_data <- read.csv('crucialGene_1.csv', header = TRUE, check.names = FALSE, row.names = 1)
#模型输入
expression_data_t <- t(expression_data)
#feature_matrix <- as.matrix(expression_data_t)[, 1:16]

# feature_matrix <- t(as.matrix(expression_data)[, 2:47])

# sample_labels <- data$Group
sample_labels <- gsub("(.*)\\_(.*)", "\\2", rownames(expression_data_t))[1:47]
table(sample_labels)

analysis_data <- data.frame(
  Status = factor(sample_labels),
  expression_data_t
)

### 随机森林模型构建 ----------------------------------------------
# 设置随机种子保证结果可复现
set.seed(123)
# 构建随机森林模型
rf_model <- randomForest(Status ~ .,data = analysis_data,ntree = 500,importance = TRUE,proximity = TRUE)

# 图1: 随机森林误差率随树数量变化
pdf("RF/1_RF_Error_Rate.pdf", width=6, height=5)
plot(rf_model, main="Random Forest Error Rate", lwd=2)

# 展示不同类别和总体误差随树数量增加的变化趋势
# 黑线: 总体误差
# 其他颜色线: 各类别误差
dev.off()
### 特征重要性分析 ------------------------------------------------
# 提取特征重要性
feature_importance <- importance(rf_model, type=1) # 使用MeanDecreaseAccuracy
feature_importance <- as.data.frame(feature_importance)
feature_importance$Gene <- rownames(feature_importance)

# 按重要性排序
feature_importance <- feature_importance[order(-feature_importance$MeanDecreaseAccuracy), ]
selected_features <- feature_importance$Gene[feature_importance$MeanDecreaseAccuracy > 0]  # 筛选正贡献特征
# 图2: 特征重要性条形图（前30个基因）
top_features <- head(feature_importance, 30)
pdf("RF/2_RF_Feature_Importance.pdf", width=10, height=8)
ggplot(top_features, aes(x=reorder(Gene, MeanDecreaseAccuracy),
                        y=MeanDecreaseAccuracy,
                        fill=MeanDecreaseAccuracy)) + geom_bar(stat="identity") + scale_fill_gradient(low="skyblue", high="firebrick") + coord_flip() + labs(title="Features by Random Forest", x="Gene", y="Mean Decrease in Accuracy") + theme_minimal() + theme(legend.position="none",
                        plot.title=element_text(hjust=0.5, size=14, face="bold"))

dev.off()

# 保存重要特征
write.table(top_features, "RF/RF_selected_genes.txt",
            sep="\t", quote=FALSE, row.names=FALSE)

rm(list = ls())
###取交集与可视化
# 安装并加载UpSetR包
# install.packages("UpSetR")
library(UpSetR)
Lasso <- read.table("Lasso/Lasso_selected_genes.txt",header = T)
SVM <- read.table('SVM-RFE/SVM_Selected_Genes.txt',header = F)
RF <- read.table('RF/RF_selected_genes.txt',header = T)

# 提取各算法选中的基因列表
gene_lists <- list(
            Lasso = Lasso$Biomarker,
            SVM = SVM$V1,
            RF = RF$Gene)

# 创建交集矩阵
upset_data <- fromList(gene_lists)
pdf("1_Lasso_SVM_RF_gene_intersection_upset.pdf", width=8.07, height=5.64)
upset(upset_data, sets = names(gene_lists),
        order.by = "freq",
        sets.bar.color = rainbow(3),
        main.bar.color ='black',
        matrix.color = rainbow(7)[7],
        text.scale = 1.2,
        mainbar.y.label = "Gene number intersected",
        sets.x.label = "Gene number selected")

dev.off()

library(VennDiagram)
venn.plot = venn.diagram(gene_lists,filename = NULL,fill = rainbow(length(gene_lists)),scaled = FALSE,cat.pos = c(-1, 0,1),cat.col = rainbow(length(gene_lists)), cat.cex = 1.2)

pdf("2_Lasso_SVM_RF_gene_intersection_venn_diagram.pdf", width = 5.49, height = 3.26)
grid.draw(venn.plot)
dev.off()

core_genes <- Reduce(intersect, gene_lists)
write.table(core_genes, "3_Lasso_SVM_RF_core_genes.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
core_genes
