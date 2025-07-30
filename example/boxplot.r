# 加载必要的包
library(tidyverse)

# 读取CSV文件
data <- read.csv("crucialGene_1.csv", header = TRUE, row.names = 1)

# 转置数据：行变为样本，列变为基因
transposed_data <- as.data.frame(t(data))

# 添加样本分组列（从行名中提取"_"后的数字）
transposed_data$Group <- ifelse(grepl("_1$", rownames(transposed_data)), "Obese", "Normal")

# 筛选目标基因
target_genes <- c("HPD", "MSC", "ANGPT2", "NGFR")
plot_data <- transposed_data[, c(target_genes, "Group")]

# 转换为长格式（便于ggplot绘图）
long_data <- plot_data %>%
  pivot_longer(
    cols = all_of(target_genes),
    names_to = "Gene",
    values_to = "Expression"
  )

# 绘制箱线图
ggplot(long_data, aes(x = Gene, y = Expression, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Normal" = "skyblue", "Obese" = "salmon")) +
  labs(
    title = "Gene Expression in Normal vs Obese Groups",
    x = "Gene Symbol",
    y = "Expression Level",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
