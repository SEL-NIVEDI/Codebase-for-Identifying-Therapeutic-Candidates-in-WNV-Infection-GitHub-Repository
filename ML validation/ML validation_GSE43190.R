####################----------------svmRadial-------------######################

# Load necessary libraries
library(caret)
library(randomForest)
library(S4Vectors)

# Read expression data
filepath <- "E:/VARSHA RAMESH/Transcriptomics/West Nile Virus/GSE136342/Expression Data_GSE43190.csv"
expression_data <- read.csv(filepath, header = TRUE)
rownames(expression_data) <- expression_data[, 1]
expression_data <- expression_data[ , -1]  # Removing the first column with gene names

# Define class labels (e.g., 8 infected, 8 control samples)
class <- factor(rep(c("control", "infected"), c(22, 17)))

# Set the hub genes
hub_genes <- c("IL15RA", "ADM", "CCL2", "IRF1", "GBP1", "BATF2", "TRIM21", "TLR7", "AIM2", "AXL")
data <- expression_data[hub_genes, ]

# Split percentage for test data (e.g., 40%)
split_ratio <- 0.4

set.seed(123)
# Separate the samples by class
control_samples <- which(class == "control")
infected_samples <- which(class == "infected")

# Determine the number of test samples for each class
n_control_test <- ceiling(length(control_samples) * split_ratio)
n_infected_test <- ceiling(length(infected_samples) * split_ratio)

# Randomly sample indices for test data from both classes
control_test_indices <- sample(control_samples, n_control_test, replace = FALSE)
infected_test_indices <- sample(infected_samples, n_infected_test, replace = FALSE)

# Combine test indices from both classes
test_indices <- c(control_test_indices, infected_test_indices)

# Split data into training and testing datasets
data.train <- data[ , -test_indices]  # Columns not in 'test_indices' for training
data.test <- data[ , test_indices]    # Columns in 'test_indices' for testing

# Transpose the data to match samples as rows and genes as columns
train_data <- as.data.frame(t(data.train))
test_data <- as.data.frame(t(data.test))

# Add class labels to train and test data
train_data$condition <- class[-test_indices]  # Class labels for training data
test_data$condition <- class[test_indices]    # Class labels for testing data

# Train Random Forest model
set.seed(1234)
fit.svm <- train(condition ~ ., data = train_data, method = "svmRadial", 
                 trControl = trainControl(method = "repeatedcv", number = 5, repeats = 5),
                 tuneLength = 10)

# Print model summary
print(fit.svm)

# Make predictions on the test set
pred.svm <- predict(fit.svm, newdata = test_data)

# Evaluate model performance
confusionMatrix(pred.svm, test_data$condition)
cm_results <- confusionMatrix(pred.svm, test_data$condition)

write.table(cm_results$table, file = "confusion_matrix_output.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

write.table(cm_results$overall, file = "confusion_matrix_overall.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

write.table(cm_results$byClass, file = "confusion_matrix_byClass.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

##---------------result-------------------##
#Confusion Matrix and Statistics

#Reference
#Prediction control infected
#control        7        1
#infected       2        6

#Accuracy : 0.8125          
#95% CI : (0.5435, 0.9595)
#No Information Rate : 0.5625          
#P-Value [Acc > NIR] : 0.03511         

#              Kappa : 0.625           

#Mcnemar's Test P-Value : 1.00000         
                                          
#            Sensitivity : 0.7778          
#            Specificity : 0.8571          
#         Pos Pred Value : 0.8750          
#         Neg Pred Value : 0.7500          
#             Prevalence : 0.5625          
#         Detection Rate : 0.4375          
#   Detection Prevalence : 0.5000          
#      Balanced Accuracy : 0.8175          
                                          
#       'Positive' Class : control         

#######################--------------RF--------------###########################
set.seed(1234)
fit.rf <- train(condition ~ ., data = train_data, method = "rf", 
                 trControl = trainControl(method = "repeatedcv", number = 6, repeats = 6),
                 tuneLength = 10)

# Print model summary
print(fit.rf)

# Make predictions on the test set
pred.rf <- predict(fit.rf, newdata = test_data)

# Evaluate model performance
confusionMatrix(pred.rf, test_data$condition)
cm_results_rf <- confusionMatrix(pred.rf, test_data$condition)

###Results same as svmRadial.

              ######-------ROC curve----------#########

library(pROC)
actual_rf <- ifelse(test_data$condition == "control", 0, 1)
predicted_probs_rf <- predict(fit.rf, newdata = test_data, type = "prob")[, "control"]
roc_curve_rf <- roc(actual_rf, predicted_probs_rf)
plot(roc_curve_rf, col = "blue", lwd = 2, main = "ROC Curve for SVM Model on Hub Genes")
auc_value_rf <- auc(roc_curve_rf)
legend("bottomright", legend = paste("AUC =", round(auc_value_rf, 3)), col = "blue", lwd = 2, cex = 1.5)
print(paste("AUC Value:", round(auc_value_rf, 3)))


