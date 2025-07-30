# -*- coding: utf-8 -*-
#
#  Created on Tue Apr 14 14:53:46 2025
#
#  @author: 真弓快车
#
#  该脚本的作用是使用LASSO进行核心基因筛选、重要性排序与诊断模型构建
#

# 设置工作目录
# setwd("E:/2025-6-15-ObesityWGCNA-machinelearning/MachineLearningForTarget")
rm(list = ls())
dir.create('Result')
setwd('Result')
# 用于Lasso回归分析
library(glmnet)
# 输入数据
input_data <- "../dataSet/2025-7-15-dataSet/crucialGene_1.csv"
# 读取输入数据
expression_data <- read.csv(input_data, header = TRUE, check.names = FALSE, row.names = 1)
# 准备模型输入
expression_data_t <- t(expression_data)
feature_matrix <- as.matrix(expression_data_t)[, 1:12]
# 从行名提取样本标签
sample_labels <- gsub("(.*)\\_(.*)", "\\2", rownames(expression_data_t))


# 构建Lasso回归模型
# alpha=1 表示Lasso回归
lasso_model <- glmnet(feature_matrix, sample_labels, family = "binomial", alpha = 1)

# Lasso系数路径图
pdf(file = "1_lasso_coefficient_paths.pdf", width = 6, height = 5.5)
plot(lasso_model)
dev.off()

# 交叉验证误差曲线
#set.seed(123)
# 十折交叉验证
cv_model <- cv.glmnet(feature_matrix, sample_labels, family = "binomial", alpha = 1, nfolds = 10)

pdf("2_cross_validation_error.pdf", width = 6, height = 5)
plot(cv_model)
dev.off()


# 特征基因相关系数排序
optimal_lambda <- cv_model$lambda.min
model_coef <- coef(cv_model, s = optimal_lambda)
# 筛选重要生物标志物
selected_features <- rownames(model_coef)[which(model_coef != 0)]
selected_features <- selected_features[-1]
# 准备特征重要性数据
feature_importance <- as.matrix(model_coef)[selected_features, ]
sorted_idx <- order(abs(feature_importance), decreasing = FALSE)
sorted_importance <- feature_importance[sorted_idx]
# 带方向特征的重要性
pdf("3_feature_importance_directional.pdf", width = 7, height = 5)

barplot(sorted_importance, horiz = TRUE, las = 1, col = ifelse(sorted_importance > 0, "tomato", "steelblue"), xlab = "Coefficient Value", main = "Biomarker Impact Direction")
dev.off()

# 特征基因重要性排序
pdf("4_feature_importance_absolute.pdf", width = 7, height = 5)

barplot(abs(sorted_importance), horiz = TRUE, las = 1, col = "tomato", xlab = "Absolute Coefficient Value", main = "Gene Importance Ranking")
dev.off()
# 保存结果
write.table(data.frame(Biomarker = selected_features, Coefficient = feature_importance), "selected_biomarkers.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# 单基因ROC曲线
# 定义图形的颜色
bioCol <- rainbow(length(selected_features), s = 0.9, v = 0.9)
# 提取所有基因表达数据 (前面只用了前50个基因，这里用全部181个)
expression_data <- t(read.csv(input_data, header = TRUE, check.names = FALSE, row.names = 1))
expression_data <- expression_data[, 1:12]
# 提取样品的分组信息
#y <- gsub("(.*)\\_(.*)\\_(.*)", "\\3", rownames(expression_data))
y <- gsub("(.*)\\_(.*)", "\\2", rownames(expression_data))
#y <- ifelse(y == "Lean", 0, 1)
print(table(y))
# 对选中的基因绘制ROC曲线
aucText <- c()
k <- 0
# 生成PDF文件
pdf(file = "5_ROC.genes.pdf", width = 9, height = 9)
# 循环绘制每个基因的ROC曲线
#feature_test <- "ANGPT2"
selected_features

for (x in as.vector(selected_features)) {
    #k <- k + 1
    k <- k + 1
    # 绘制ROC曲线
    roc1 <- pROC::roc(y, as.numeric(t(expression_data)[x, ]))
    if (k == 1) {
        plot(roc1, print.auc = F, col = bioCol[k], legacy.axes = T, main = "", lwd = 3)
        aucText <- c(aucText, paste0(x, ", AUC=", sprintf("%.3f", roc1$auc[1])))
    } else {
        plot(roc1, print.auc = F, col = bioCol[k], legacy.axes = T, main = "", lwd = 3, add = TRUE)
        aucText <- c(aucText, paste0(x, ", AUC=", sprintf("%.3f", roc1$auc[1])))
    }
}
# 添加图例
legend("bottomright", aucText, lwd = 3, bty = "n", cex = 1.3, col = bioCol[1:k], inset = c(0.05, 0))
# title(main= feature_test)
dev.off()



# 多联合模型ROC曲线
# 提取被选基因的表达数据和样本标签
selected_data <- data.frame(expression_data[, selected_features[1:4], drop = FALSE], Status = factor(sample_labels))
# 转换分组为0/1变量
#selected_data$Status <- ifelse(selected_data$Status == "DN", 1, 0)
library(pROC)
# 构建包含所有选定基因的模型
full_model <- glm(Status ~ ., data = selected_data, family = binomial())
# 计算综合预测概率
selected_data$pred_prob <- predict(full_model, type = "response")
# 绘制ROC曲线
roc_obj <- roc(Status ~ pred_prob, data = selected_data)
# 生成PDF
pdf("6_logistical_curve.pdf", width = 6.62, height = 5.25)
# 绘制多基因联合ROC曲线
plot(roc_obj, print.auc = TRUE)
# AUC值表示多基因联合预测能力
# 比较单个基因AUC和联合模型AUC可评估基因组合效果
dev.off()

