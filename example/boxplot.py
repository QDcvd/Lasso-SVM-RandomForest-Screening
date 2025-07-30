import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import stats

# 读取CSV文件
df = pd.read_csv('crucialGene_1.csv')

# 设置基因符号为索引
df.set_index('symbol', inplace=True)

# 转置数据框，使基因成为列，样本成为行
df_t = df.T.reset_index()
df_t.rename(columns={'index': 'sample'}, inplace=True)

# 从样本名称中提取分组信息（0=正常组，1=肥胖组）
df_t['group'] = df_t['sample'].apply(lambda x: 'Obese' if x.endswith('_1') else 'Normal')

# 选择目标基因
target_genes = ['HPD', 'NCEH1', 'ANGPT2', 'ALDH1A3']

# 设置绘图风格
# sns.set(style="whitegrid")
plt.figure(figsize=(20, 5))

# 对每个基因进行循环处理
for i, gene in enumerate(target_genes, 1):
    plt.subplot(1, 4, i)
    
    # 提取当前基因的数据
    gene_data = df_t[['sample', 'group', gene]].dropna()
    
    # 分离两组数据
    normal_vals = gene_data[gene_data['group'] == 'Normal'][gene].values
    obese_vals = gene_data[gene_data['group'] == 'Obese'][gene].values
    
    # 进行t检验
    t_stat, p_val = stats.ttest_ind(obese_vals, normal_vals, equal_var=False)
    
    # 确定显著性标记
    if p_val < 0.001:
        sig_symbol = '***'
    elif p_val < 0.01:
        sig_symbol = '**'
    elif p_val < 0.05:
        sig_symbol = '*'
    else:
        sig_symbol = 'ns'
    
    # 绘制箱型图
    sns.boxplot(x='group', y=gene, data=gene_data, order=['Normal', 'Obese'], 
                palette=['#1f77b4', '#ff7f0e'])
    
    # 添加标题和标签
    plt.title(f'{gene}', fontsize=14)
    plt.xlabel('')
    plt.ylabel('Expression Level', fontsize=12)
    
    # 添加p值和显著性标记
    plt.text(0.5, np.max(gene_data[gene]) * 0.95, 
            #  f'p = {p_val:.4f} {sig_symbol}',
            f'p = {p_val}',
            ha='center', fontsize=6, 
            bbox=dict(facecolor='white', alpha=0.8))

# 调整布局
plt.tight_layout(pad=3.0)
# plt.suptitle('Gene Expression in Normal vs Obese Groups', fontsize=16, y=1.02)
plt.savefig('gene_expression_boxplots.png', dpi=300, bbox_inches='tight')
plt.show()