
rownames(other_cohort_clinical_data) <- other_cohort_clinical_data$Sample.ID

HER2_list_clinical   <- other_cohort_clinical_data[which(other_cohort_clinical_data$HER2_Status == 'positive'),]
HER2_list   <- rownames(HER2_list_clinical)

other_cohort_HER2_mRNA <- other_cohort_mRNA[,which(colnames(other_cohort_mRNA) %in% HER2_list)]


# combat -------------------------------------------------------

PCA_data <- cbind(FUSCC_HER2_RNA,other_cohort_HER2_mRNA)
pca_group <- matrix(data = NA,ncol = 2,nrow  = ncol(PCA_data))
pca_group[,1]  <- colnames(PCA_data)
pca_group[,2]  <- c(rep('FUSCC_HER2',ncol(FUSCC_HER2_RNA)),rep('other_cohort',ncol(other_cohort_HER2_mRNA)))
colnames(pca_group) <- c('sample','group')
pca_group <- as.data.frame(pca_group)

library(sva)
modcombat = model.matrix(~1, data = pca_group)
batch = pca_group$group
combat_edata = ComBat(dat=PCA_data, batch=batch, mod=modcombat,   
                      par.prior=TRUE, prior.plots=TRUE)

FUSCC_HER2_combat_RNA <- combat_edata[differential_genes,1:ncol(FUSCC_HER2_RNA)]
other_cohort_combat_RNA <- combat_edata[differential_genes,ncol(FUSCC_HER2_RNA)+1:ncol(combat_edata)]


# RF -------------------------------------------------------

library(dplyr)
library(tidyr)
library(stringr)
library(stringi)
library(pheatmap)
library(export)
library(RColorBrewer)
library(caTools)
library(randomForest)
library(janitor)
library(e1071)
library(caret)

RNA_RF <- as.data.frame(t(FUSCC_HER2_RNA))
RNA_RF$patient_code <- rownames(RNA_RF)
RF_data <- merge(RNA_RF,HER2_NMF,by='patient_code')

set.seed(1234)
sample_data <- sample.split(RF_data[,'patient_code'], SplitRatio = 2/3)
gam_train <- RF_data[sample_data == TRUE, ] %>% select(-1)
gam_train <- clean_names(gam_train)
gam_test <- RF_data[sample_data == FALSE, ] %>% select(-1)
gam_test <- clean_names(gam_test)

trControl <- trainControl(method = "cv",
                          number = 10,
                          search = "grid")

# select mtry
tuneGrid <- expand.grid(.mtry = c(1:10))
rf_mtry <- train(nmf ~ .,          
                 data = gam_train,  
                 method = "rf",
                 metric = "Accuracy",
                 tuneGrid = tuneGrid,
                 trControl = trControl,
                 importance = TRUE,
                 nodesize = 10,
                 ntree = 300)
print(rf_mtry)
best_mtry <- rf_mtry$bestTune$mtry 

# select maxnodes
store_maxnode <- list()
tuneGrid <- expand.grid(.mtry = best_mtry)

for (maxnodes in 10:20) {              
  set.seed(1234)
  rf_maxnode <- train(nmf ~ .,
                      data = gam_train,
                      method = "rf",
                      metric = "Accuracy",
                      tuneGrid = tuneGrid,
                      trControl = trControl,
                      importance = TRUE,
                      nodesize = 14,
                      maxnodes = maxnodes,
                      ntree = 300)
  current_iteration <- toString(maxnodes)
  store_maxnode[[current_iteration]] <- rf_maxnode
  print(maxnodes)
}
results_mtry <- resamples(store_maxnode)
summary(results_mtry)  # maxnodes 

# select ntrees
store_maxtrees <- list()
for (ntree in seq(250, 2000, by = 50)) {
  set.seed(5678)
  rf_maxtrees <- train(nmf ~ .,
                       data = gam_train,
                       method = "rf",
                       metric = "Accuracy",
                       tuneGrid = tuneGrid,
                       trControl = trControl,
                       importance = TRUE,
                       nodesize = 14,
                       maxnodes = maxnodes,  
                       ntree = ntree)
  key <- toString(ntree)
  store_maxtrees[[key]] <- rf_maxtrees
  print(ntree)
}
results_tree <- resamples(store_maxtrees)
summary(results_tree)  
results_tree_df <- results_tree$values %>% 
  select(-1) %>% 
  t() %>% 
  data.frame()
results_tree_df$MEDIAN <- apply(results_tree_df, 1, median)  
results_tree_df$MEAN <- apply(results_tree_df[, 1:10], 1, mean)

# random forest establishment
gam_train$nmf <- as.factor(gam_train$nmf)
rf_fit <- randomForest(nmf ~ ., data = gam_train, importance = T, proximity = T, 
                       mtry = best_mtry, maxnodes = maxnodes, ntree = tree_num)  

plot(rf_fit)
varImpPlot(rf_fit)
rf_fit_imp <- varImp(rf_fit)
mdsplot <- MDSplot(rf_fit, fac = gam_train$nmf)

# random forest validation
rf_prediction <- predict(rf_fit, gam_test)
gam_test$nmf <- as.factor(gam_test$nmf)
rf_pred_cm <- confusionMatrix(rf_prediction, gam_test$nmf)
rf_pred_cmat <- rf_pred_cm$table %>% 
  data.frame() %>% 
  pivot_wider(id_cols = Prediction,
              names_from = Reference,
              values_from = Freq) %>% 
  data.frame() %>% 
  `rownames<-`(.$Prediction) %>% 
  select(-1)

rf_pred_cmat2 <- apply(rf_pred_cmat, 1, function(x){x/sum(x)}) %>% t()

pheatmap(rf_pred_cmat2, annotation_legend = F,
         cluster_rows = F, cluster_cols = F,
         show_rownames = T, show_colnames = T,
         legend = T, display_numbers = T, border_color = NA,
         fontsize_row = 12, fontsize_col = 12, fontsize_number = 12, number_color = 'white',
         color = rev(colorRampPalette(c('DarkRed', 'white', 'darkblue'))(100)))




# other_cohort HER2 subtype -------------------------------------------------------------------

other_cohort_RNA <- t(other_cohort_RNA)
other_cohort_RNA <- as.data.frame(other_cohort_RNA)
other_cohort_RNA <- clean_names(other_cohort_RNA)

other_cohort_RNA_NMF <- predict(rf_fit, other_cohort_RNA)
other_cohort_RNA_NMF <- as.data.frame(other_cohort_RNA_NMF)