library(R.matlab)
library(readr)
library(Matrix)


DATA_PATH = "Data_YL/"

groups = read.csv(file.path(DATA_PATH, "communities_aparc.csv"))
groups = groups[order(groups$id_ROI), ][, 3]
avg_matrix = sapply(unique(groups), FUN = function(x) as.numeric(groups == x)/sum(groups == x))

kinship = readMat(file.path(DATA_PATH, "kinship.mat"))
K_G = as(kinship$K[[1]], "TsparseMatrix")
id_subjects = as.character(kinship$K[[2]])
rownames(K_G) = id_subjects;
colnames(K_G) = id_subjects;

fun_connectomes = list()
fun_connectomes_averaged = list()
session = 1
for (i in seq(id_subjects)) {
  ts = as.matrix(read_csv(file.path(DATA_PATH, paste0('3T_HCP1200_aparc_ts2/', id_subjects[i], '_', session, '.csv')), col_names = F, col_types = cols()))
  fun_connectomes[[i]] = cor(ts)
  fun_connectomes_averaged[[i]] = cor(ts %*% avg_matrix)
}
names(fun_connectomes) = id_subjects
names(fun_connectomes_averaged) = id_subjects

str_connectomes = list()
str_connectomes_averaged = list()
for (i in seq(id_subjects)) {
  tryCatch({
    idx_ctx = c(1:3,5:35,36:38,40:70) # indices of valid cortical regions
    str_raw = read_delim(file.path(DATA_PATH, paste0('HCP1200_desikan_str/sub-',
                                                     id_subjects[i], '_ses-', session, '_run-1_dwi_Desikan_space-MNI152NLin6_res-1x1x1_connectome.csv')), col_names = F, col_types = cols(), delim = " ")
    str_raw_spmat = sparseMatrix(str_raw$X1, str_raw$X2, x = str_raw$X3, symmetric = TRUE)
    str_connectomes_full = as.matrix(str_raw_spmat)
    str_connectomes[[i]] = str_connectomes_full[idx_ctx, idx_ctx]
    str_connectomes_averaged[[i]] = t(avg_matrix) %*% str_connectomes_full[idx_ctx, idx_ctx] %*% avg_matrix
    diag(str_connectomes_averaged[[i]]) = 0
  },
  warning = function(cond){
    str_connectomes[[i]] = NULL
    str_connectomes_averaged[[i]] = NULL
  },
  error = function(e){
    str_connectomes[[i]] = NULL
    str_connectomes_averaged[[i]] = NULL
  })
}
names(str_connectomes) = id_subjects
names(str_connectomes_averaged) = id_subjects

common_subjects = names(str_connectomes)[!sapply(str_connectomes, is.null)]
common_subjects = setdiff(common_subjects, "180836")
K_G = K_G[common_subjects, common_subjects]
rownames(K_G) = as.character(common_subjects);
colnames(K_G) = as.character(common_subjects);

K_G = as.matrix(K_G)
order = hclust(as.dist(-K_G))$order
K_G = K_G[order, order]

common_subjects = colnames(K_G)
fun_connectomes = fun_connectomes[common_subjects]
fun_connectomes_averaged = fun_connectomes_averaged[common_subjects]
str_connectomes = str_connectomes[common_subjects]
str_connectomes_averaged = str_connectomes_averaged[common_subjects]

X = read.csv(file.path(DATA_PATH, "conf.csv"))
rownames(X) = X$subject
X = X[common_subjects, c("Age", "Age.2", "Sex", "FS_IntraCranial_Vol..1.3.", "FS_BrainSeg_Vol..1.3.")]
X = cbind(rep(1, nrow(X)), X)

data = list(
  fun_connectomes = fun_connectomes,
  fun_connectomes_averaged = fun_connectomes_averaged,
  str_connectomes = str_connectomes,
  str_connectomes_averaged = str_connectomes_averaged,
  K_G = K_G,
  X = X,
  groups = groups
)

saveRDS(data, file.path(DATA_PATH, "clean_data.rds"))
