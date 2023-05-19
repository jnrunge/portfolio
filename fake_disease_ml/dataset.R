# I want to try out some algorithms and I make a fake dataset to get it.
# this data set has a couple of genotyped loci and some additional data
# then each individual gets a binary trait based on the random data and
# I want to see how well we can uncover the logic behind the binary trait assignment
library(tidyverse)
library(e1071)
library(randomForest)
count_individuals = 2000
count_loci = 4

genotype_matrix <-
  matrix(
    rbinom(count_individuals * count_loci, 1, 0.5) + rbinom(count_individuals * count_loci, 1, 0.5),
    nrow = count_individuals,
    ncol = count_loci
  )
genotype_df <- as_tibble(genotype_matrix)
colnames(genotype_df) <- paste0("locus_", 1:count_loci)

additional_data <-
  tibble(
    ID = 1:count_individuals,
    weight = rnorm(count_individuals, mean = 80, sd = 5),
    age = sample(20:80, count_individuals, replace = TRUE)
  )

data_combined <- bind_cols(additional_data, genotype_df)



make_binary_trait = function(row) {
  # here are the rules how the binary trait is generated
  data_combined_slice <- slice(data_combined, row)
  #print(data_combined_slice)
  individual_risk <- 0
  individual_risk <-
    individual_risk + (pull(data_combined_slice, locus_1)) * 0.1
  individual_risk <-
    individual_risk + (pull(data_combined_slice, locus_2)) * 0.1
  individual_risk <-
    individual_risk + (pull(data_combined_slice, locus_3)) * 0.1
  individual_risk <-
    individual_risk + ((pull(data_combined_slice, locus_1) == 2) * (pull(data_combined_slice, locus_2) ==
                                                                      2)) * 0.4
  individual_risk <-
    individual_risk + ((pull(data_combined_slice, locus_1) == 2) * (pull(data_combined_slice, locus_2) ==
                                                                      2) * (pull(data_combined_slice, locus_3) == 2)) * 0.9 # multiple homozygous alleles needed
  if (individual_risk > 1) {
    individual_risk <- 1
  }
  phenotype <- rbinom(1, 1, individual_risk)
  return(phenotype)
}

add_affected_trait <- function(row) {
  data_combined_slice <- slice(data_combined, row)
  mean <-
    40 # the average in the healthy population of a trait, say VO2max
  mean <-
    mean - (pull(data_combined_slice, phenotype) == 1) * 5 # decrease in 5 for those with the phenotype
  mean <- mean - pull(data_combined_slice, weight) * 0.1
  mean <- mean - pull(data_combined_slice, age) * 0.1
  return(rnorm(1, mean, sd = 3))
}



data_combined$phenotype <-
  as.factor(unlist(lapply(
    1:nrow(data_combined), make_binary_trait
  )))
data_combined$affected_trait <-
  unlist(lapply(1:nrow(data_combined), add_affected_trait))
data_combined$weight = as.numeric(scale(data_combined$weight, center = TRUE, scale =
                                          TRUE))
data_combined$age = as.numeric(scale(data_combined$age, center = TRUE, scale =
                                       TRUE))
data_combined$affected_trait = as.numeric(scale(
  data_combined$affected_trait,
  center = TRUE,
  scale = TRUE
))

summary(data_combined$phenotype)


# plots showing how difficult it will be

library(ggplot2)
ggplot(data_combined, aes(age, fill = phenotype)) +
  geom_histogram(position = "fill")

ggplot(data_combined, aes(weight, fill = phenotype)) +
  geom_histogram(position = "fill")

ggplot(data_combined, aes(affected_trait)) +
  geom_histogram()

ggplot(data_combined, aes(affected_trait, fill = phenotype)) +
  geom_histogram(position = "fill")

ggplot(data_combined, aes(weight, affected_trait, color = phenotype)) +
  geom_point(shape = 1)

ggplot(data_combined, aes(age, affected_trait, color = phenotype)) +
  geom_point(shape = 1)

####

# simple logistic regression

# basically impossible to generate anything from it other than luckily combining the right loci with interactions like second example
summary(glm(
  paste0(
    "phenotype~weight+age+affected_trait",
    paste0("+locus_", 1:count_loci, collapse = "")
  ),
  data_combined,
  family = "binomial"
))
# but to find that one would have to do so many comparisons that we could not trust the results
summary(
  glm(
    "phenotype~weight+age+affected_trait+locus_1*locus_2*locus_3",
    data_combined,
    family = "binomial"
  )
)

# clustering



# unsupervised random forest

data_combined.urf <-
  randomForest(select(data_combined,-1,-(ncol(data_combined) - 1)), ntree =
                 count_individuals * 10)
MDSplot(data_combined.urf, data_combined$phenotype)

# classification random forest

data_combined.rf <-
  randomForest(
    phenotype ~ .,
    data = data_combined %>% select(-1),
    importance = TRUE,
    proximity = TRUE,
    ntree = count_individuals * 10,
    maxnodes = 500
  )
data_combined.rf$confusion
## Look at variable importance:
round(importance(data_combined.rf), 2)

# not great, not bad, lets try to tune the model

library(caret)
# split the data
ids <- createDataPartition(data_combined$phenotype, p = 0.8, list = F)
train_set <- data_combined[ids, ]
test_set <- data_combined[-ids, ]

#here we set the resampling method for hyperparameters tuning
#in this case we choose 10-fold cross validation
cn <- trainControl(method = "cv", number = 10)

grid <- expand.grid(mtry = 2:(ncol(train_set) - 1))

fit <-
  train(
    phenotype ~ .,
    data = train_set,
    
    method = "rf",
    trControl = cn,
    tuneGrid = grid,
    ntree = count_individuals * 10,
    maxnodes = 500
  )

#predict the output on the test set.
p <- predict(fit, test_set %>% select(-(ncol(test_set) - 1)))

confusionMatrix(p, test_set$phenotype)

# SVM
svm.model <- svm(phenotype ~ ., data = train_set %>% select(-1))
svm.model
svm.pred  <-
  predict(svm.model, select(test_set,-1,-(ncol(test_set) - 1)))
table(pred = svm.pred, true = test_set %>% pull(phenotype))
