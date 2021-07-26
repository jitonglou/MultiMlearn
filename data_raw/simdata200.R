#### Simulate data ####
rseed = 0
set.seed(rseed)
simdata200 = simulate.data(
   N = 200, p = 20, K = 4, J = 4,
   propensity_func = pi_true, main_func = mu_true, interaction_func = delta_true
)
usethis::use_data(simdata200) # https://community.rstudio.com/t/how-to-store-and-document-data-to-be-used-within-package/4359

group = 1
data_covariate = simdata200 %>%
  filter(cluster == group) %>%
  select(ID, cluster, reward, treatment, starts_with("feature"))
n_feature = ncol(data_covariate) - 4
feature_name = paste0("feature_",seq_len(n_feature))

table(simdata200$treatment)
table(data_covariate$treatment)

#### Estimate prognostic score ####
## tuning parameters
metric = "RMSE"
tunegrid = expand.grid(.mtry=seq(1,n_feature/3, 1), .ntree=seq(500,1000,500))
control = caret::trainControl(method="repeatedcv", number=10, repeats=3)
## train model
set.seed(rseed)
prognostic_fit = caret::train(
  reward~., data=data_covariate %>% select(-ID, -cluster, -treatment),
  method=rfcv2("Regression"), metric=metric,
  tuneGrid=tunegrid, trControl=control
)
print(prognostic_fit)
## estimators
pat_prog = predict(
  prognostic_fit, newdata=data_covariate %>% select(-ID, -cluster, -treatment)
)
caret::RMSE(pat_prog, prognostic_fit$trainingData$.outcome) # RMSE


#### Estimate propensity score ####
## tuning parameters
metric = "Kappa"
tunegrid = expand.grid(.mtry=seq(1,sqrt(n_feature), 1), .ntree=seq(500,1000,500))
control = caret::trainControl(method="repeatedcv", number=10, repeats=3)
## train model
set.seed(rseed)
propensity_fit = caret::train(
  treatment~., data=data_covariate %>% select(-ID, -cluster, -reward),
  method=rfcv2("Classification"), metric=metric,
  tuneGrid=tunegrid, trControl=control
)
print(propensity_fit)

## estimators
pat_prop = predict(
  propensity_fit, newdata=data_covariate %>% select(-ID, -cluster, -reward), type = "prob"
)
caret::confusionMatrix(
  predict(propensity_fit, data_covariate %>% select(-ID, -cluster, -reward)),
  data_covariate$treatment
) # confusion matrix

#### Create the final dataset used for Mlearning ####
data_feature = data.frame(
  data_covariate %>% mutate(reward_res = reward),
  pat_prop %>% select(-levels(simdata200$treatment)[1]), # remove propensity scores of the first treatment to avoid colinearity of propensity scores
  prog = pat_prog
)

## compute the Euclidean distance between features of subjects
dist_feature = data_feature %>%
  select(-ID, -cluster, -treatment, -reward, -reward_res) %>%
  dist() %>%
  as.matrix()
dim(dist_feature)

idx_data_feature = data.frame(ID=data_feature$ID, index=1:nrow(data_feature)) # index is used to select the row of the distance matrix

table(simdata200$cluster)
summary(data_feature)

#### [OLD] Learn ITR (OVO, svm) ####
ksvm.grid = expand.grid(C = 2^(-1:1)) # ksvm.grid = expand.grid(C = 2^(-1:1))
nfolds_obj = 3
nfolds = 3 # more folds yield more training samples
eps = 1e-16
g_func = function(x){abs(x)} # weights on the outcome in SVM
SNN = 1 # size of the matched set


set.seed(rseed)
ITR_OVO = ksvm_cv_wvf_obj_knownfold_v2(
  dat=data_feature, idx=idx_data_feature,
  trts=levels(data_feature$treatment), SNN=SNN, nfolds_obj=nfolds_obj,
  kernel="rbfdot", kpar="automatic", eps=eps, nfolds=nfolds, tuneGrid=ksvm.grid,
  delta=max(dist_feature), propensity=pat_prop, dist_mat=dist_feature
)

#### [NEW] Learn ITR (OVO, svm) ####
ksvm.grid = expand.grid(C = 2^(-1:1)) # ksvm.grid = expand.grid(C = 2^(-1:1))
nfolds_outer = 3
nfolds_inner = 3 # more folds yield more training samples
g_func = function(x){abs(x)} # weights on the outcome in SVM
max_size = 1 # size of the matched set

set.seed(rseed)
ITR_OVO = weighted.ksvm.2.cv(
  data=data_feature, idx=idx_data_feature,
  trts=levels(data_feature$treatment), max_size=max_size,
  delta=max(dist_feature), dist_mat=dist_feature, g_func=g_func,
  kernel="rbfdot", kpar="automatic",
  nfolds_outer=nfolds_outer, nfolds_inner=nfolds_inner, tuneGrid=ksvm.grid, propensity=pat_prop,
  foldid_outer=NULL
)

data=data_feature
idx=idx_data_feature
trts=levels(data_feature$treatment)
delta=max(dist_feature)
dist_mat=dist_feature
kernel="rbfdot"
kpar="automatic"
tuneGrid=ksvm.grid
propensity=pat_prop
foldid_outer=NULL


pi_vec = rep(NA, nrow(ITR_OVO$prediction))
for(i in 1:nrow(ITR_OVO$prediction)){
  pi_vec[i] = sub_prop[i,ITR_OVO$prediction$treatment[i]]
}
df_evf = data.frame(ITR_OVO$prediction, pi=pi_vec, iter=rseed)


#### Print results ####
## Confusion matrix of ITR and assigned treatment
caret::confusionMatrix(as.factor(df_evf$vote), df_evf$treatment)

## EVF of one-size-fits-all rules
df_evf %>%
  dplyr::group_by(treatment) %>%
  dplyr::summarize(evf=sum(reward/pi)/sum(1/pi)) %>%
  as.data.frame()

## Different versions of EVF
my_evf = function(reward, pi, x, y){
  flag = x==y
  return(sum((reward/pi)[flag])/sum((1/pi)[flag]))
}

summary_evf = df_evf %>%
  dplyr::left_join(
    simdata200 %>% select(ID, pi_true, contains("opt")),
    by = "ID") %>%
  dplyr::group_by(fold) %>%
  dplyr::summarize(
    evf_max = my_evf(reward, pi_true, treatment, treatment_opt),
    evf_ipw = my_evf(reward, pi, treatment, vote),
    evf_ipw_true = my_evf(reward, pi_true, treatment, vote),
    evf_ipw_opt = my_evf(reward, pi_opt, treatment_opt, vote),
    n=dplyr::n()
  ) %>%
  as.data.frame()

summary(summary_evf)
