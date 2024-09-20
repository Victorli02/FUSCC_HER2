
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

