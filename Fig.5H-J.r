library(ggplot2)
library(ggpubr)
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\20.3Bikini")
#set color
blue   <- "#1b90be"
red ="#e64b35"
grey   <- "#8693ab"
darkred   <- "#F2042C"
lightred  <- "#FF7FBF"
#tmp$riskScore=log2(tmp$riskScore+1)
# #----------------------------Kim-------------------------
# tmp <- read.csv("response.risk.Kim.csv", row.names = NULL, check.names = F, header = T, stringsAsFactors = F)
# head(tmp)
# table(tmp$Response)
# # 设置组间对比，排列组???
# my_comparisons <- list( c("Response","Non-response")
#                         #, 
#                         #c("Red-TGFB", "Grey-TGFB"), 
#                         #c("Red-TGFB", "Blue-TGFB"), 
#                         #c("Pink-TGFB","Grey-TGFB"),
#                         #c("Pink-TGFB","Blue-TGFB"),
#                         #c("Grey-TGFB","Blue-TGFB")
#                         )
# ggplot(data = tmp,aes(x = Response, #分组列名
#                       y = Kim, #连续变量列名
#                       fill = Response))+ #按分组填充颜???
#   scale_fill_manual(values = c(red,blue,  lightred, darkred)) + #用自定义颜色填充
#   geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
#               size = 0.8, color="black") +
#   geom_boxplot(notch = TRUE, outlier.size = -1, 
#                color="black", lwd=0.8, alpha = 0.7) +
#   geom_point(shape = 21, size=2, # 点的性状和大???
#              position = position_jitterdodge(), # 让点散开
#              color="black", alpha = 1) +
#   theme_classic() + 
#   ylab("Expression") +
#   xlab("Kim's cohort") +
#   theme(axis.text.x = element_text(angle = 0, hjust = 0.5,vjust=1, size = 15),
#         #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
#         axis.ticks = element_line(size=0.2, color="black"),
#         axis.ticks.length = unit(0.2, "cm"),
#         legend.position = "none",
#         axis.title = element_text(size = 15),
#         axis.text = element_text(size = 15)) +
# # 如果不要组间比较就注释掉下面这行
#   stat_compare_means(comparisons = my_comparisons,label = "p.signif") 
# #+stat_compare_means(method = "kruskal.test", label.y = min(tmp$Kim))
# ggsave("bikini.response.risk.Kim.pdf", width = 4, height = 4)


#----------------------------IMvigor210-------------------------
tmp <- read.csv("response.risk.IMvigor210.csv", row.names = NULL, check.names = F, header = T, stringsAsFactors = F)
head(tmp)
table(tmp$Response)
# 设置组间对比，排列组???
my_comparisons <- list( c("Response","Non-response")
                        #, 
                        #c("Red-TGFB", "Grey-TGFB"), 
                        #c("Red-TGFB", "Blue-TGFB"), 
                        #c("Pink-TGFB","Grey-TGFB"),
                        #c("Pink-TGFB","Blue-TGFB"),
                        #c("Grey-TGFB","Blue-TGFB")
)
ggplot(data = tmp,aes(x = Response, #分组列名
                      y = MIR100HG, #连续变量列名
                      fill = Response))+ #按分组填充颜???
  scale_fill_manual(values = c(red,blue,  lightred, darkred)) + #用自定义颜色填充
  geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
              size = 0.8, color="black") +
  geom_boxplot(notch = TRUE, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7) +
  geom_point(shape = 21, size=2, # 点的性状和大???
             position = position_jitterdodge(), # 让点散开
             color="black", alpha = 1) +
  theme_classic() + 
  ylab("MIR100HG Expression") +
  xlab("IMvigor210's cohort") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5,vjust=1, size = 15),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  # 如果不要组间比较就注释掉下面这行
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") 
#+stat_compare_means(method = "kruskal.test", label.y = min(tmp$IMvigor210))
ggsave("bikini.response.risk.IMvigor210.pdf", width = 4, height = 4)


# #----------------------------Braun-------------------------
# tmp <- read.csv("response.risk.Braun.csv", row.names = NULL, check.names = F, header = T, stringsAsFactors = F)
# head(tmp)
# table(tmp$Response)
# # 设置组间对比，排列组???
# my_comparisons <- list( c("Response","Non-response")
#                         #, 
#                         #c("Red-TGFB", "Grey-TGFB"), 
#                         #c("Red-TGFB", "Blue-TGFB"), 
#                         #c("Pink-TGFB","Grey-TGFB"),
#                         #c("Pink-TGFB","Blue-TGFB"),
#                         #c("Grey-TGFB","Blue-TGFB")
# )
# ggplot(data = tmp,aes(x = Response, #分组列名
#                       y = Braun, #连续变量列名
#                       fill = Response))+ #按分组填充颜???
#   scale_fill_manual(values = c(red,blue,  lightred, darkred)) + #用自定义颜色填充
#   geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
#               size = 0.8, color="black") +
#   geom_boxplot(notch = TRUE, outlier.size = -1, 
#                color="black", lwd=0.8, alpha = 0.7) +
#   geom_point(shape = 21, size=2, # 点的性状和大???
#              position = position_jitterdodge(), # 让点散开
#              color="black", alpha = 1) +
#   theme_classic() + 
#   ylab("Expression") +
#   xlab("Braun's cohort") +
#   theme(axis.text.x = element_text(angle = 0, hjust = 0.5,vjust=1, size = 15),
#         #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
#         axis.ticks = element_line(size=0.2, color="black"),
#         axis.ticks.length = unit(0.2, "cm"),
#         legend.position = "none",
#         axis.title = element_text(size = 15),
#         axis.text = element_text(size = 15)) +
#   # 如果不要组间比较就注释掉下面这行
#   stat_compare_means(comparisons = my_comparisons,label = "p.signif") 
# #+stat_compare_means(method = "kruskal.test", label.y = min(tmp$Braun))
# ggsave("bikini.response.risk.Braun.pdf", width = 4, height = 4)


# #----------------------------Liu-------------------------
# tmp <- read.csv("response.risk.Liu.csv", row.names = NULL, check.names = F, header = T, stringsAsFactors = F)
# head(tmp)
# table(tmp$Response)
# # 设置组间对比，排列组???
# my_comparisons <- list( c("Response","Non-response")
#                         #, 
#                         #c("Red-TGFB", "Grey-TGFB"), 
#                         #c("Red-TGFB", "Blue-TGFB"), 
#                         #c("Pink-TGFB","Grey-TGFB"),
#                         #c("Pink-TGFB","Blue-TGFB"),
#                         #c("Grey-TGFB","Blue-TGFB")
# )
# ggplot(data = tmp,aes(x = Response, #分组列名
#                       y = Liu, #连续变量列名
#                       fill = Response))+ #按分组填充颜???
#   scale_fill_manual(values = c(red,blue,  lightred, darkred)) + #用自定义颜色填充
#   geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
#               size = 0.8, color="black") +
#   geom_boxplot(notch = TRUE, outlier.size = -1, 
#                color="black", lwd=0.8, alpha = 0.7) +
#   geom_point(shape = 21, size=2, # 点的性状和大???
#              position = position_jitterdodge(), # 让点散开
#              color="black", alpha = 1) +
#   theme_classic() + 
#   ylab("Expression") +
#   xlab("Liu's cohort") +
#   theme(axis.text.x = element_text(angle = 0, hjust = 0.5,vjust=1, size = 15),
#         #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
#         axis.ticks = element_line(size=0.2, color="black"),
#         axis.ticks.length = unit(0.2, "cm"),
#         legend.position = "none",
#         axis.title = element_text(size = 15),
#         axis.text = element_text(size = 15)) +
#   # 如果不要组间比较就注释掉下面这行
#   stat_compare_means(comparisons = my_comparisons,label = "p.signif") 
# #+stat_compare_means(method = "kruskal.test", label.y = min(tmp$Liu))
# ggsave("bikini.response.risk.Liu.pdf", width = 4, height = 4)


# #----------------------------Miao-------------------------
# tmp <- read.csv("response.risk.Miao.csv", row.names = NULL, check.names = F, header = T, stringsAsFactors = F)
# head(tmp)
# table(tmp$Response)
# # 设置组间对比，排列组???
# my_comparisons <- list( c("Response","Non-response")
#                         #, 
#                         #c("Red-TGFB", "Grey-TGFB"), 
#                         #c("Red-TGFB", "Blue-TGFB"), 
#                         #c("Pink-TGFB","Grey-TGFB"),
#                         #c("Pink-TGFB","Blue-TGFB"),
#                         #c("Grey-TGFB","Blue-TGFB")
# )
# ggplot(data = tmp,aes(x = Response, #分组列名
#                       y = Miao, #连续变量列名
#                       fill = Response))+ #按分组填充颜???
#   scale_fill_manual(values = c(red,blue,  lightred, darkred)) + #用自定义颜色填充
#   geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
#               size = 0.8, color="black") +
#   geom_boxplot(notch = TRUE, outlier.size = -1, 
#                color="black", lwd=0.8, alpha = 0.7) +
#   geom_point(shape = 21, size=2, # 点的性状和大???
#              position = position_jitterdodge(), # 让点散开
#              color="black", alpha = 1) +
#   theme_classic() + 
#   ylab("Expression") +
#   xlab("Miao's cohort") +
#   theme(axis.text.x = element_text(angle = 0, hjust = 0.5,vjust=1, size = 15),
#         #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
#         axis.ticks = element_line(size=0.2, color="black"),
#         axis.ticks.length = unit(0.2, "cm"),
#         legend.position = "none",
#         axis.title = element_text(size = 15),
#         axis.text = element_text(size = 15)) +
#   # 如果不要组间比较就注释掉下面这行
#   stat_compare_means(comparisons = my_comparisons,label = "p.signif") 
# #+stat_compare_means(method = "kruskal.test", label.y = min(tmp$Miao))
# ggsave("bikini.response.risk.Miao.pdf", width = 4, height = 4)


#----------------------------Gide-------------------------
tmp <- read.csv("response.risk.Gide.csv", row.names = NULL, check.names = F, header = T, stringsAsFactors = F)
head(tmp)
table(tmp$Response)
# 设置组间对比，排列组???
my_comparisons <- list( c("Response","Non-response")
                        #, 
                        #c("Red-TGFB", "Grey-TGFB"), 
                        #c("Red-TGFB", "Blue-TGFB"), 
                        #c("Pink-TGFB","Grey-TGFB"),
                        #c("Pink-TGFB","Blue-TGFB"),
                        #c("Grey-TGFB","Blue-TGFB")
)
ggplot(data = tmp,aes(x = Response, #分组列名
                      y = MSC.AS1, #连续变量列名
                      fill = Response))+ #按分组填充颜???
  scale_fill_manual(values = c(red,blue,  lightred, darkred)) + #用自定义颜色填充
  geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
              size = 0.8, color="black") +
  geom_boxplot(notch = TRUE, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7) +
  geom_point(shape = 21, size=2, # 点的性状和大???
             position = position_jitterdodge(), # 让点散开
             color="black", alpha = 1) +
  theme_classic() + 
  ylab("MSC-AS1 Expression") +
  xlab("Gide's cohort") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5,vjust=1, size = 15),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  # 如果不要组间比较就注释掉下面这行
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") 
#+stat_compare_means(method = "kruskal.test", label.y = min(tmp$Gide))
ggsave("bikini.response.risk.Gide.pdf", width = 4, height = 4)


#----------------------------Nathanson-------------------------
tmp <- read.csv("response.risk.Nathanson.csv", row.names = NULL, check.names = F, header = T, stringsAsFactors = F)
head(tmp)
table(tmp$Response)
# 设置组间对比，排列组???
my_comparisons <- list( c("Response","Non-response")
                        #, 
                        #c("Red-TGFB", "Grey-TGFB"), 
                        #c("Red-TGFB", "Blue-TGFB"), 
                        #c("Pink-TGFB","Grey-TGFB"),
                        #c("Pink-TGFB","Blue-TGFB"),
                        #c("Grey-TGFB","Blue-TGFB")
)
ggplot(data = tmp,aes(x = Response, #分组列名
                      y = DNM3OS, #连续变量列名
                      fill = Response))+ #按分组填充颜???
  scale_fill_manual(values = c(red,blue,  lightred, darkred)) + #用自定义颜色填充
  geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
              size = 0.8, color="black") +
  geom_boxplot(notch = TRUE, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7) +
  geom_point(shape = 21, size=2, # 点的性状和大???
             position = position_jitterdodge(), # 让点散开
             color="black", alpha = 1) +
  theme_classic() + 
  ylab("DNM3OS Expression") +
  xlab("Nathanson's cohort") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5,vjust=1, size = 15),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  # 如果不要组间比较就注释掉下面这行
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") 
#+stat_compare_means(method = "kruskal.test", label.y = min(tmp$Nathanson))
ggsave("bikini.response.risk.Nathanson.pdf", width = 4, height = 4)


#----------------------------GSE91061-------------------------
tmp <- read.csv("response.risk.GSE91061.csv", row.names = NULL, check.names = F, header = T, stringsAsFactors = F)
head(tmp)
table(tmp$Response)
# 设置组间对比，排列组???
my_comparisons <- list( c("Response","Non-response")
                        #, 
                        #c("Red-TGFB", "Grey-TGFB"), 
                        #c("Red-TGFB", "Blue-TGFB"), 
                        #c("Pink-TGFB","Grey-TGFB"),
                        #c("Pink-TGFB","Blue-TGFB"),
                        #c("Grey-TGFB","Blue-TGFB")
)
ggplot(data = tmp,aes(x = Response, #分组列名
                      y = MAGI2.AS3, #连续变量列名
                      fill = Response))+ #按分组填充颜???
  scale_fill_manual(values = c(red,blue,  lightred, darkred)) + #用自定义颜色填充
  geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
              size = 0.8, color="black") +
  geom_boxplot(notch = TRUE, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7) +
  geom_point(shape = 21, size=2, # 点的性状和大???
             position = position_jitterdodge(), # 让点散开
             color="black", alpha = 1) +
  theme_classic() + 
  ylab("MAGI2-AS3 Expression") +
  xlab("GSE91061's cohort") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5,vjust=1, size = 15),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  # 如果不要组间比较就注释掉下面这行
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") 
#+stat_compare_means(method = "kruskal.test", label.y = min(tmp$GSE91061))
ggsave("bikini.response.risk.GSE91061.pdf", width = 4, height = 4)


# #----------------------------GSE78220-------------------------
# tmp <- read.csv("response.risk.GSE78220.csv", row.names = NULL, check.names = F, header = T, stringsAsFactors = F)
# head(tmp)
# table(tmp$Response)
# # 设置组间对比，排列组???
# my_comparisons <- list( c("Response","Non-response")
#                         #, 
#                         #c("Red-TGFB", "Grey-TGFB"), 
#                         #c("Red-TGFB", "Blue-TGFB"), 
#                         #c("Pink-TGFB","Grey-TGFB"),
#                         #c("Pink-TGFB","Blue-TGFB"),
#                         #c("Grey-TGFB","Blue-TGFB")
# )
# ggplot(data = tmp,aes(x = Response, #分组列名
#                       y = GSE78220, #连续变量列名
#                       fill = Response))+ #按分组填充颜???
#   scale_fill_manual(values = c(red,blue,  lightred, darkred)) + #用自定义颜色填充
#   geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
#               size = 0.8, color="black") +
#   geom_boxplot(notch = TRUE, outlier.size = -1, 
#                color="black", lwd=0.8, alpha = 0.7) +
#   geom_point(shape = 21, size=2, # 点的性状和大???
#              position = position_jitterdodge(), # 让点散开
#              color="black", alpha = 1) +
#   theme_classic() + 
#   ylab("Expression") +
#   xlab("GSE78220's cohort") +
#   theme(axis.text.x = element_text(angle = 0, hjust = 0.5,vjust=1, size = 15),
#         #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
#         axis.ticks = element_line(size=0.2, color="black"),
#         axis.ticks.length = unit(0.2, "cm"),
#         legend.position = "none",
#         axis.title = element_text(size = 15),
#         axis.text = element_text(size = 15)) +
#   # 如果不要组间比较就注释掉下面这行
#   stat_compare_means(comparisons = my_comparisons,label = "p.signif") 
# #+stat_compare_means(method = "kruskal.test", label.y = min(tmp$GSE78220))
# ggsave("bikini.response.risk.GSE78220.pdf", width = 4, height = 4)
# 
# 
# #----------------------------GSE67501-------------------------
# tmp <- read.csv("response.risk.GSE67501.csv", row.names = NULL, check.names = F, header = T, stringsAsFactors = F)
# head(tmp)
# table(tmp$Response)
# # 设置组间对比，排列组???
# my_comparisons <- list( c("Response","Non-response")
#                         #, 
#                         #c("Red-TGFB", "Grey-TGFB"), 
#                         #c("Red-TGFB", "Blue-TGFB"), 
#                         #c("Pink-TGFB","Grey-TGFB"),
#                         #c("Pink-TGFB","Blue-TGFB"),
#                         #c("Grey-TGFB","Blue-TGFB")
# )
# ggplot(data = tmp,aes(x = Response, #分组列名
#                       y = GSE67501, #连续变量列名
#                       fill = Response))+ #按分组填充颜???
#   scale_fill_manual(values = c(red,blue,  lightred, darkred)) + #用自定义颜色填充
#   geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
#               size = 0.8, color="black") +
#   geom_boxplot(notch = TRUE, outlier.size = -1, 
#                color="black", lwd=0.8, alpha = 0.7) +
#   geom_point(shape = 21, size=2, # 点的性状和大???
#              position = position_jitterdodge(), # 让点散开
#              color="black", alpha = 1) +
#   theme_classic() + 
#   ylab("Expression") +
#   xlab("GSE67501's cohort") +
#   theme(axis.text.x = element_text(angle = 0, hjust = 0.5,vjust=1, size = 15),
#         #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
#         axis.ticks = element_line(size=0.2, color="black"),
#         axis.ticks.length = unit(0.2, "cm"),
#         legend.position = "none",
#         axis.title = element_text(size = 15),
#         axis.text = element_text(size = 15)) +
#   # 如果不要组间比较就注释掉下面这行
#   stat_compare_means(comparisons = my_comparisons,label = "p.signif") 
# #+stat_compare_means(method = "kruskal.test", label.y = min(tmp$GSE67501))
# ggsave("bikini.response.risk.GSE67501.pdf", width = 4, height = 4)


