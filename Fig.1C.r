library(ggplot2)
library(ggpubr)
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\01CAF.data")
#set color
blue   <- "#1b90be"
red ="#e64b35"
grey   <- "#8693ab"
darkred   <- "#F2042C"
lightred  <- "#FF7FBF"
#tmp$riskScore=log2(tmp$riskScore+1)
#----------------------------IntegratedScore-------------------------
tmp <- read.csv("CS.bikini.csv", row.names = NULL, check.names = F, header = T, stringsAsFactors = F)
head(tmp)
table(tmp$CancerType)

ggplot(data = tmp,aes(x = CancerType, #分组列名
                      y = IntegratedScore, #连续变量列名
                      fill = CancerType))+ #按分组填充颜???
  #scale_fill_manual(values = c(red,blue,  lightred, darkred)) + #用自定义颜色填充
  geom_violin(alpha = 0.4, position = position_dodge(width = .75),
              size = 0.8, color="black") +
  geom_boxplot(notch = TRUE, outlier.size = -1,
               color="black", lwd=0.8, alpha = 0.7) +
  geom_point(shape = 21, size=3, # 点的性状和大???
             position = position_jitterdodge(), # 让点散开
             color="black", 
             alpha = 1) +
  theme_classic() +
  ylab("Integrated Score") +
  xlab("Pan-cancer cohort") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5,vjust=1, size = 36),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 36),
        axis.text = element_text(size = 36)) +
# 如果不要组间比较就注释掉下面这行
#  stat_compare_means(comparisons = my_comparisons,label = "p.signif")
#+stat_compare_means(method = "kruskal.test", label.y = min(tmp$IntegratedScore))
ggsave("bikini.CancerType.risk.IntegratedScore.pdf", width = 36, height = 12)


