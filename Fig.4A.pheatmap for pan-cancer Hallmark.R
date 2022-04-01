library(pheatmap)
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\14pan-cancer.hallmark\\05plot")
test=read.table("heatmap.pan-cancer.correlation.Hallmark.txt", header=T, sep="\t", check.names=F, row.names=1)
data=as.matrix(test)
pdf("heatmap_hallmark.pdf", width=8, height=8)
pheatmap(test,cluster_cols = F,cluster_rows = F,
         color = colorRampPalette(c(rep("#3c5488",2), "white", rep("#e64b35",2)))(100),
         border_color = "black",
         treeheight_row = 10, treeheight_col = 50,
         #display_numbers = TRUE,number_color = "blue",
         #display_numbers = matrix(ifelse(test > 0.5, "*", ""), nrow(test)),
         #cellwidth = 12, cellheight = 12,
         main = "Pan-cancer correlation between CAFDL signature and Hallmark gene sets",
         #gaps_row = c(11, 21, 32),
         )
dev.off()


#½Ì³Ì£ºhttps://www.jianshu.com/p/1c55ea64ff3f