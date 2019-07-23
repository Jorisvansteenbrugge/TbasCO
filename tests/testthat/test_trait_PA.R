context("Test Trait presence/absence")
input <-
"Gene;Sample1;Annotation;Bin
x;1;K1;Bin1
x;1;K2;Bin1
x;1;K3;Bin1
x;1;K4;Bin1
x;1;K1;Bin2
x;1;K2;Bin2
x;1;K3;Bin2
x;1;K2;Bin3
x;1;K3;Bin3
x;1;K2;Bin4
x;1;K3;Bin4
x;1;K4;Bin4
x;1;K2;Bin4
"
modules.input <- "Module;Annotation
M1;K1
M1;K2
M1;K3
M1;K4
M2;K2
M2;K3
M2;K4
M2;K5"


data <- read.csv2(text=input,stringsAsFactors = F)
modules <- read.csv2(text=modules.input, stringsAsFactors = F)

annotation.db <- Create.Module.groups(modules)

expect_identical(as.character(annotation.db$module.dict$M1), c('K1', 'K2', 'K3', 'K4') )

# Test annotation presence absence
anno_pa <- Get_annotation_presence_absence(data, c('Bin1','Bin2' ,'Bin3','Bin4' ),
                                           annotation.db )
res <- matrix(nrow=5, ncol = 4, c(T,T,F,F,T,T,T,T,T,T,T,T,T,F,F,T,F,F,F,F), byrow = T)
colnames(res) <- data$Bin %>% unique %>% sort
rownames(res) <- annotation.db$`all annotations in a module`

expect(identical(res,anno_pa))


features <- Get_matrix_features(data,annotation.db)

res <- matrix(nrow = 2, ncol = 4, c(T,T,F,T,T,F,F,T), byrow = T)
colnames(res) <- data$Bin %>% unique %>% sort
rownames(res) <- c("M1" , "M2")

expect(identical(features$trait_presence_absence, res))
