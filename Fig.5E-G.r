library(ggplot2)
library(ggpubr)
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\19pan-cancer.msi.tmb")
#set color
blue   <- "#1b90be"
red ="#e64b35"
green   <- "#5FC1C2"
darkred   <- "#F2042C"
lightred  <- "#FF7FBF"
#tmp$Rank30=log2(tmp$Rank30+1)
#----------------------------COAD-------------------------
tmp <- read.csv("MSI.COAD.csv", row.names = NULL, check.names = F, header = T, stringsAsFactors = F)
head(tmp)
table(tmp$MSIstatus)
# 设置组间对比，排列组�?
my_comparisons <- list( c("MSS","MSI-L"), 
                        c("MSI-L", "MSI-H"), 
                        c("MSS", "MSI-H")
                        #, 
                        #c("Pink-TGFB","Grey-TGFB"),
                        #c("Pink-TGFB","Blue-TGFB"),
                        #c("Grey-TGFB","Blue-TGFB")
                        )
ggplot(data = tmp,aes(x = MSIstatus, #分组列名
                      y = Rank30, #连续变量列名
                      fill = MSIstatus))+ #按分组填充颜�?
  scale_fill_manual(values = c( blue, green,red,darkred)) + #用自定义颜色填充
  geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
              size = 0.8, color="black") +
  geom_boxplot(notch = TRUE, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7) +
  geom_point(shape = 21, size=2, # 点的性状和大�?
             position = position_jitterdodge(), # 让点散开
             color="black", alpha = 1) +
  theme_classic() + 
  ylab("TEL score") +
  xlab("COAD") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5,vjust=1, size = 15),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
# 如果不要组间比较就注释掉下面这行
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") 
#+stat_compare_means(method = "kruskal.test", label.y = min(tmp$Rank30))
p
ggsave("bikini.MSIstatus.Rank30.COAD.pdf", width = 4, height = 4)

#----------------------------READ-------------------------
tmp <- read.csv("MSI.READ.csv", row.names = NULL, check.names = F, header = T, stringsAsFactors = F)
head(tmp)
table(tmp$MSIstatus)
# 设置组间对比，排列组�?
my_comparisons <- list( c("MSS","MSI-L"), 
                        c("MSI-L", "MSI-H"), 
                        c("MSS", "MSI-H")
                        #, 
                        #c("Pink-TGFB","Grey-TGFB"),
                        #c("Pink-TGFB","Blue-TGFB"),
                        #c("Grey-TGFB","Blue-TGFB")
)
ggplot(data = tmp,aes(x = MSIstatus, #分组列名
                      y = Rank30, #连续变量列名
                      fill = MSIstatus))+ #按分组填充颜�?
  scale_fill_manual(values = c( blue, green,red,darkred)) + #用自定义颜色填充
  geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
              size = 0.8, color="black") +
  geom_boxplot(notch = TRUE, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7) +
  geom_point(shape = 21, size=2, # 点的性状和大�?
             position = position_jitterdodge(), # 让点散开
             color="black", alpha = 1) +
  theme_classic() + 
  ylab("TEL score") +
  xlab("READ") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5,vjust=1, size = 15),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  # 如果不要组间比较就注释掉下面这行
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") 
#+stat_compare_means(method = "kruskal.test", label.y = min(tmp$Rank30))
p
ggsave("bikini.MSIstatus.Rank30.READ.pdf", width = 4, height = 4)


#----------------------------COADREAD-------------------------
tmp <- read.csv("MSI.COADREAD.csv", row.names = NULL, check.names = F, header = T, stringsAsFactors = F)
head(tmp)
table(tmp$MSIstatus)
# 设置组间对比，排列组�?
my_comparisons <- list( c("MSS","MSI-L"), 
                        c("MSI-L", "MSI-H"), 
                        c("MSS", "MSI-H")
                        #, 
                        #c("Pink-TGFB","Grey-TGFB"),
                        #c("Pink-TGFB","Blue-TGFB"),
                        #c("Grey-TGFB","Blue-TGFB")
)
ggplot(data = tmp,aes(x = MSIstatus, #分组列名
                      y = Rank30, #连续变量列名
                      fill = MSIstatus))+ #按分组填充颜�?
  scale_fill_manual(values = c( blue, green,red,darkred)) + #用自定义颜色填充
  geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
              size = 0.8, color="black") +
  geom_boxplot(notch = TRUE, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7) +
  geom_point(shape = 21, size=2, # 点的性状和大�?
             position = position_jitterdodge(), # 让点散开
             color="black", alpha = 1) +
  theme_classic() + 
  ylab("TEL score") +
  xlab("COADREAD") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5,vjust=1, size = 15),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  # 如果不要组间比较就注释掉下面这行
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") 
#+stat_compare_means(method = "kruskal.test", label.y = min(tmp$Rank30))
p
ggsave("bikini.MSIstatus.Rank30.COADREAD.pdf", width = 4, height = 4)


#----------------------------STAD-------------------------
tmp <- read.csv("MSI.STAD.csv", row.names = NULL, check.names = F, header = T, stringsAsFactors = F)
head(tmp)
table(tmp$MSIstatus)
# 设置组间对比，排列组�?
my_comparisons <- list( c("MSS","MSI-L"), 
                        c("MSI-L", "MSI-H"), 
                        c("MSS", "MSI-H")
                        #, 
                        #c("Pink-TGFB","Grey-TGFB"),
                        #c("Pink-TGFB","Blue-TGFB"),
                        #c("Grey-TGFB","Blue-TGFB")
)
ggplot(data = tmp,aes(x = MSIstatus, #分组列名
                      y = Rank30, #连续变量列名
                      fill = MSIstatus))+ #按分组填充颜�?
  scale_fill_manual(values = c( blue, green,red,darkred)) + #用自定义颜色填充
  geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
              size = 0.8, color="black") +
  geom_boxplot(notch = TRUE, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7) +
  geom_point(shape = 21, size=2, # 点的性状和大�?
             position = position_jitterdodge(), # 让点散开
             color="black", alpha = 1) +
  theme_classic() + 
  ylab("TEL score") +
  xlab("STAD") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5,vjust=1, size = 15),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  # 如果不要组间比较就注释掉下面这行
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") 
#+stat_compare_means(method = "kruskal.test", label.y = min(tmp$Rank30))
p
ggsave("bikini.MSIstatus.Rank30.STAD.pdf", width = 4, height = 4)


