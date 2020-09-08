# Takes the 200 most important genes from random forest, input them in david to get the conversion ID 
# and original gene name. Then p-value (form t-test) associated to original gene name is used in pathfindR

small.df = ex2[,1:110]
f <- factor(c(rep("Tumor",80), rep("Normal",30)))
library("genefilter")
ttest <- rowttests(small.df,f)
ttest

# from random forest output, it takes the most 200 important genes
important.genes = sort(rf2$importance, decreasing = T)[1:200]
probe.names = rownames(rf2$importance)

top200 = probe.names[order(rf2$importance, decreasing = T)[1:200]]
imp = read.table("out_david2.txt", sep = "\t")

l = read.table("list_removed.txt", sep = "\t") # list file coming from david conversion
row.names(l) = 1:196


# T-test performed to associate a p-value to each gene
ttest$Name = as.character(rownames(ttest))

ttest = ttest[imp$V1,]  
ttest$ID = l
ttest= ttest[,c("ID","p.value")]

ID = regmatches(imp$V4, gregexpr("(?<=\\().*?(?=\\))", imp$V4, perl=T)) # to filter out not relevant IDs
ID = as.factor(unlist(ID))

write(ID,"list_to_remove.txt", ncolumns = 1)
write.table(ID,"list_to_remove.txt", sep = "\t")
l = read.table("list_removed.txt", sep = "\t")
row.names(l) = 1:196
imp$ID = l
imp = imp[,c("ID", "p.val")]

# Write the file composed of IDs and P-values associated with them to be used in pathffindR
write.csv("pathfindr", imp) 
View(imp)
