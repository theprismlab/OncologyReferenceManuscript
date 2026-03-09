#' Calculates and returns univariate analysis results by correlating each column of Y with features sets from depmap.org.
#'
#' @param Y : Matrix n x m, make sure rownames are depmap_id's and columns are named. There can be NA's.
#' @param file : Please point out to the downloaded depmap_datasets.h5 file.
#' @param features : You can give a subset of available feature sets, but if left as NULL, it computes for everything.
#' @param n.X.min : Results with less thana given sample size are dropped, default is 100
#' @param q_val.max : Results with q-values less than q_val.max are dropped, default is 0.2
#' @param rank.max : Results with ranks (by q-value) greater than rank.max are dropped, default is 250
#'
#' @return Returns a data-table with each row corresponds to a particular (feature, feature_set, y) triplet. See linear_model for the other columns.
#' @export
#'
#' @examples
#'
univariate_biomarker_table <- function(Y, file, features = NULL, n.X.min = 250, v.X.min = 0.0025, q_val_max = .1, rank.max = 250){
  require(tidyverse)
  require(magrittr)
  require(rhdf5)
  require(WGCNA)
  
  
  if(!is.matrix(Y)){
    Y <- as.matrix(Y)
  }
  
  
  if(is.null(features)){
    features <-  substr(unique(dplyr::filter(h5ls(file), name == "mat")$group),2,100)
  }else{
    features <- intersect(features, substr(unique(dplyr::filter(h5ls(file), name == "mat")$group),2,100))
  }
  print(features)
  
  RESULTS <- list()
  for(feat in features){
    
    
    X <- read_dataset(file , feat)
    cl = intersect(rownames(X), rownames(Y))
    X <- X[cl, ]
    
    if((dim(X)[1] >= n.X.min) & (dim(X)[2] > 0)){
      
      print(paste0(feat, " - ", dim(X)[1] ,'x', dim(X)[2]))
      
      RESULTS[[feat]] <- linear_model(X = X, Y = Y[cl, , drop = FALSE],
                                             v.X.min = v.X.min, n.min = n.X.min) %>% 
        dplyr::filter(rank <= rank.max, q_val <= q_val_max) %>% 
        dplyr::rename(feature = x) %>%
        dplyr::mutate(feature_set = feat) 
    }
    print(feat)
  }
  
  RESULTS <- dplyr::bind_rows(RESULTS)
  return(RESULTS)
}



#' Fits a simple linear model by regressing y over each column of X.
#'
#' @param X : matrix of n by m, it can have NA's in it, columns should be named.
#' @param Y : matrix of n by p, rows are in the same order of with the rows of X.
#' @param v.th : Minimum variance for the columns of X, columns with smaller variances will be dropped.
#' @param n.min : Minimum number of finite pairs between columns of X and Y, column pairs not satisfying this condition will be dropped.
#' 
#'
#' @return A data frame with: x (corresponding column of X), y (corresponding column of Y), correlation_coeff (Pearson correlation), 
#'         regression_coeff (regression coefficient), p_val / q_val (homoskedastic p-value / q-value),
#'         n (number of non-na samples), rank (rank of the significance for each column of Y), 
#
#' @export
#'
#' @examples
#'
linear_model <- function(X, Y, v.X.min = 0.0025, n.min = 100) {
  require(tidyverse)
  require(magrittr)
  require(matrixStats)
  require(WGCNA)
  
  
  cor.table <- WGCNA::corAndPvalue(X,Y, use = "p")[c(1,2,5)] %>%
    reshape2::melt() %>% 
    tidyr::pivot_wider(names_from = L1, values_from = value) %>% 
    dplyr::filter(nObs >= n.min) %>% 
    dplyr::mutate(x = as.character(Var1), y= as.character(Var2)) %>% 
    dplyr::select(-Var1, -Var2) %>%
    dplyr::rename(correlation_coef = cor, p_val = p, n = nObs) 
  
  
  masks <- is.finite(Y)
  masks <- masks[, !duplicated(t(masks)), drop = F]
  colnames(masks) <- paste0("m", 1:dim(masks)[2])
  
  map <- apply(masks, 2, function(m) apply(is.finite(Y) == m, 2, all)) %>% 
    apply(1, which.max)
  map <- tibble(y = names(map),
                m = colnames(masks)[map])
  masks[masks == 0] = NA
  vX <- apply(masks, 2, function(m) colVars(X * m, na.rm = T))
  vX <- vX[,map$m]
  colnames(vX) <- map$y
  
  masks <- is.finite(X)
  masks <- masks[, !duplicated(t(masks)), drop = F]
  colnames(masks) <- paste0("m", 1:dim(masks)[2])
  map <- apply(masks, 2, function(m) apply(is.finite(X) == m, 2, all)) %>% 
    apply(1, which.max)
  map <- tibble(x = names(map),
                m = colnames(masks)[map])
  masks[masks == 0] = NA
  vY <- apply(masks, 2, function(m) colVars(Y * m, na.rm = T))
  vY <- vY[,map$m, drop = F]
  colnames(vY) <- map$x
  
  cor.table %<>% 
    dplyr::left_join(vX %>% 
                       reshape2::melt() %>% 
                       dplyr::mutate(x = as.character(Var1), y = as.character(Var2)) %>% 
                       dplyr::rename(var.x = value)  %>% 
                       dplyr::select(x,y, var.x)) %>% 
    dplyr::left_join(vY %>% 
                       reshape2::melt() %>% 
                       dplyr::mutate(y = as.character(Var1), x = as.character(Var2)) %>% 
                       dplyr::rename(var.y = value)  %>% 
                       dplyr::select(x,y, var.y)) %>% 
    dplyr::mutate(regression_coef = correlation_coef * sqrt(var.y) / sqrt(var.x)) %>% 
    dplyr::filter( n >= n.min, var.x >= v.X.min)
  
  if(nrow(cor.table) > 0 ){
    cor.table %<>% 
      dplyr::group_by(y) %>% 
      dplyr::mutate(q_val = p.adjust(p_val, method = "BH")) %>% 
      dplyr::arrange(q_val) %>% dplyr::mutate(rank = 1:n()) %>%
      dplyr::group_by(y, q_val) %>% dplyr::mutate(rank = min(rank, na.rm = T)) %>% 
      dplyr::ungroup() %>%
      dplyr::select(x,y,correlation_coef, regression_coef, q_val, rank, p_val, n, var.x, var.y)
    
  }
  
  return(cor.table) 
}



#' Exports individual datasets from depmap_datasets.h5 file. You can specify the row names either as ccle_names or depmap_ids.
#'
#' @param file : Please point out to the downloaded depmap_datasets.h5 file.
#' @param dataset : The dataset you want to export. You can check the available datasets with substr(setdiff(rhdf5::h5ls(file)$group, "/"),2,100)
#' @param rownames_depmap_ids : Default TRUE, you can get the rownames as ccle_names by switching to FALSE.
#'
#' @return The requested data matrix.
#' @export
#'
#' @examples
read_dataset <- function(file = '/data/biomarker/current/depmap_datasets_public.h5', dataset, rownames_depmap_ids = TRUE) {
  require(rhdf5)
  if(word(file, sep = fixed("://")) %in% c("s3", "http", "https")){
    s3 = TRUE
    print(paste0("Reading ", file, " from S3"))
  } else{
    s3 = FALSE
    print(paste0("Reading ", file, " from local"))
  }
  X <- h5read(file, name = paste0(dataset, "/mat"), s3 = s3)
  row_meta <- h5read(file, name = paste0(dataset, "/row_meta"), s3 = s3)
  column_meta <- h5read(file, name = paste0(dataset, "/column_meta"), s3 = s3)
  colnames(X) <- column_meta$column_name
  if(rownames_depmap_ids){
    rownames(X) <- row_meta$ModelID
  }else{
    rownames(X) <- row_meta$CCLEName
  }
  
  X <- X[rownames(X) != "NA", colnames(X) != "NA", drop = FALSE]
  X <- X[!duplicated(rownames(X)), !duplicated(colnames(X)), drop = FALSE]
  return(X)
}






#' Title
#'
#' @param Y 
#' @param file 
#' @param compound_annotations 
#' @param features 
#' @param rank.max 
#'
#' @returns
#' @export
#'
#' @examples
target_recovery <- function(Y, file, compound_annotations,  features = c("CRISPR", "RNAi",  "Expression", "Mutation", "CopyNumber",  "Fusion",  "Lineage",   "Repurposing.Primary"), rank.max = 100, q.max = 0.1, n.min = 250) {
  require(tidyverse)
  require(magrittr)
  
  # Compute the univariate biomarkers
  bm <- univariate_biomarker_table(Y, file , features = features, 
                                          rank.max = rank.max, q_val_max = q.max, n.X.min = n.min)
  
  if (nrow(bm) > 0) {
    bm %<>% dplyr::mutate(cn = word(y, 1,2, sep = fixed("::")))
  }
  
  
  # Mapping features into gene symbols
  feats <- bm %>%
    dplyr::distinct(feature) %>%
    dplyr::mutate(feature2 = word(feature, sep = fixed("."))) %>%
    dplyr::mutate(
      feature3 = word(feature2, sep = fixed("--")),
      feature4 = word(feature2, 2, sep = fixed("--"))
    ) %>%
    tidyr::pivot_longer(c("feature2", "feature3", "feature4"),
                        names_to = "dummy",
                        values_to = "GeneSymbolOfTargets") %>%
    dplyr::select(-dummy) %>%
    tidyr::drop_na() %>%
    dplyr::distinct()
  
  # Highlight the targets 
  bm <- compound_annotations %>%
    dplyr::distinct(cn, GeneSymbolOfTargets) %>%
    tidyr::separate_rows(GeneSymbolOfTargets, sep = ";") %>%
    tidyr::drop_na() %>% 
    dplyr::distinct() %>% 
    dplyr::inner_join(feats) %>% 
    dplyr::select(-GeneSymbolOfTargets) %>%
    dplyr::distinct() %>% 
    dplyr::mutate(is.target = TRUE) %>%
    dplyr::right_join(bm)
  
  
  feats <- bm %>%
    dplyr::distinct(feature, feature_set, is.target, cn)
  
  # Read the lineage matrix
  L <- read_dataset(file, "Lineage")
  cl = intersect(rownames(L), rownames(Y))
  L <- L[cl, ]
  
  
  res <- list()
  ix <- 1
  for (feat in  unique(feats$feature_set)) {
    Z <- read_dataset(file, feat)
    cl = intersect(rownames(Z), rownames(Y))
    Z <- Z[cl, ]
    
    # all the target features
    a <- feats %>%
      dplyr::filter(is.target) %>%
      dplyr::distinct(feature_set, feature)
    
    # features in feat
    b <- feats %>%
      dplyr::filter(feature_set == feat) %>%
      .$feature %>% unique() %>% as.character()
    
    # for each feature in a, compute the correlations with each feature of b in Z
    if (nrow(a) > 0) {
      cors <- list()
      kx <- 1
      for (ff in unique(a$feature_set)) {
        Z2 <- read_dataset(file, ff)
        cll = intersect(cl, rownames(Z2))
        Z2 <- Z2[cll, ]
        a2 <- dplyr::filter(a, feature_set == ff)$feature %>% unique %>% as.character()
        
        cors[[kx]] <- WGCNA::cor(Z2[, a2, drop = F], Z[cll, b, drop = F], use = "p") %>%
          reshape2::melt()  %>%
          dplyr::rename(
            target = Var1,
            feature = Var2,
            feature_target.cor = value
          ) %>%
          dplyr::mutate(
            feature_set_target = ff,
            target = as.character(target),
            feature = as.character(feature)
          ) %>%
          as_tibble()
        
        kx <- kx + 1
      }
      
      cors %<>%
        dplyr::bind_rows()
      
      tars2 <- bm %>%
        dplyr::filter(is.target) %>%
        dplyr::rename(
          target = feature,
          target.cor = correlation_coef,
          feature_set_target = feature_set
        ) %>%
        dplyr::left_join(cors) %>%
        dplyr::distinct(
          cn,
          target,
          feature,
          target.cor,
          feature_target.cor,
          feature_set_target
        ) %>%
        dplyr::mutate(feature_set = feat)
      
      temp <- bm %>%
        dplyr::filter(feature_set == feat) %>%
        dplyr::left_join(tars2) %>%
        dplyr::mutate(target.cor = ifelse(is.na(target.cor), 0 , target.cor)) %>%
        dplyr::mutate(
          target_part_cor_coef = (correlation_coef - feature_target.cor * target.cor) / sqrt(1 - target.cor^2) / sqrt(1 - feature_target.cor^2)
        ) %>%
        dplyr::group_by(cn, y, feature, feature_set) %>%
        dplyr::arrange(target_part_cor_coef^2 * (1 - target.cor^2), target) %>%
        dplyr::mutate(dix = 1:n()) %>%
        dplyr::ungroup() %>%
        dplyr::filter(dix == 1) %>%
        dplyr::select(-dix)
      
    } else{
      temp <- bm %>%
        dplyr::filter(feature_set == feat)
    }
    
    
    # lineage correlations for each b in Z
    cl <- intersect(rownames(Z), rownames(L))
    l.cor <- WGCNA::cor(L[cl, , drop = F], Z[cl, b, drop = F], use = "p") %>%
      reshape2::melt() %>%
      dplyr::rename(
        lineage = Var1,
        feature = Var2,
        feature_lineage.cor = value
      ) %>%
      dplyr::mutate(lineage = as.character(lineage),
                    feature = as.character(feature)) %>%
      as_tibble()
    
    tars3 <- bm %>%
      dplyr::filter(feature_set == "Lineage") %>%
      dplyr::rename(lineage = feature, lineage.cor = correlation_coef) %>%
      dplyr::left_join(l.cor) %>%
      dplyr::distinct(cn,
                      lineage,
                      feature,
                      lineage.cor,
                      feature_lineage.cor) %>%
      dplyr::mutate(feature_set = feat)
    
    
    
    temp <- temp %>%
      dplyr::left_join(tars3) %>%
      dplyr::mutate( lineage_part_cor_coef = (correlation_coef - feature_lineage.cor * lineage.cor) / sqrt(1 - lineage.cor^2) / sqrt(1 - feature_lineage.cor^2)) %>%
      dplyr::group_by(cn, y, feature, feature_set) %>%
      dplyr::arrange(lineage_part_cor_coef^2 * (1 - lineage.cor^2), lineage) %>%
      dplyr::mutate(dix = 1:n()) %>%
      dplyr::ungroup() %>%
      dplyr::filter(dix == 1) %>%
      dplyr::select(-dix) %>%
      dplyr::mutate(is.target = !is.na(is.target))
    res[[ix]] <- temp
    
    ix <- ix + 1
  }
  
  res <- dplyr::bind_rows(res) %>%
    dplyr::mutate(target_part_cor_coef = ifelse(is.na(target_part_cor_coef), correlation_coef, target_part_cor_coef),
                  target.cor = ifelse(is.na(target.cor), 0, target.cor),
                  lineage_part_cor_coef = ifelse(is.na(lineage_part_cor_coef), correlation_coef, lineage_part_cor_coef),
                  lineage.cor = ifelse(is.na(lineage.cor), 0, lineage.cor),
                  tf = (1 - target.cor^2)*target_part_cor_coef^2 / correlation_coef^2,
                  lf = (1 - lineage.cor^2)*lineage_part_cor_coef^2 / correlation_coef^2) %>%
    dplyr::mutate(status = ifelse(is.target, "Target", ifelse(tf < 0.5, "Target-Correlate", ifelse(lf < 0.5, "Lineage-Correlate", "Other")))) %>%
    dplyr::left_join(compound_annotations)
  
  return(res)
}



biomarker_suite_rf <- function(X, Y, biomarker_file, CompoundList, test_samples = NULL, bm_th = 0.05, bm_R = 10,  bm_R2 = 50, features =c("CRISPR", "RNAi", "Expression", "Mutation", "CopyNumber", "Fusion", "Lineage")){
  require(tidyverse)
  require(ranger)
  
  train <- setdiff(rownames(Y), test_samples) %>% intersect(rownames(X))
  test <- intersect(rownames(Y), test_samples) %>% intersect(rownames(X))
  
  bm.auc <- target_recovery(Y[train,, drop = F], biomarker_file, CompoundList, features = features) #, rc = F, q.max = 1, rank.max = 100000)
  
  targets <- CompoundList %>% 
    dplyr::distinct(cn, CompoundName,GeneSymbolOfTargets) %>% 
    tidyr::separate_rows(GeneSymbolOfTargets, sep = ";") %>% 
    dplyr::mutate(GeneSymbolOfTargets = trimws(GeneSymbolOfTargets)) %>% 
    dplyr::distinct() %>%
    tidyr::drop_na()
  
  
  
  targets <- tibble(cn = colnames(X)) %>%  
    dplyr::mutate(t1 =  word(word(cn, 2, sep = fixed("_")), sep = fixed(".")),
                  t2 = word(word(cn, 2, sep = fixed("_")), sep = fixed("--")),
                  t3 = word(word(cn, 2, sep = fixed("_")),-1, sep = fixed("--"))) %>% 
    tidyr::pivot_longer(c(t1,t2,t3), values_to = "GeneSymbolOfTargets", names_to = "d") %>% 
    dplyr::filter(GeneSymbolOfTargets != "X", !is.na(GeneSymbolOfTargets)) %>%
    dplyr::select(-d) %>% 
    dplyr::distinct() %>% 
    dplyr::rename(cn.feat = cn) %>% 
    dplyr::inner_join(targets)
  
  
  selected_features <- bm.auc %>% 
    dplyr::distinct(cn, CompoundName, feature_set, feature, rank, correlation_coef, q_val, status) %>%
    dplyr::filter((status == "Other") | (feature_set == "Lineage")) %>% 
    dplyr::filter(((rank <= bm_R) & (correlation_coef^2 >= bm_th)) | (rank <= 1)) 
  
  if(nrow(selected_features) > 0){
    selected_features <- selected_features %>%
      dplyr::group_by(cn) %>% 
      dplyr::arrange(q_val) %>%
      dplyr::mutate(rank_ = 1:n()) %>% 
      dplyr::group_by(cn, q_val) %>%
      dplyr::mutate(rank_ = min(rank_, na.rm = T)) %>% 
      dplyr::ungroup() %>% 
      dplyr::filter(rank_ <= bm_R2) %>% 
      dplyr::ungroup() 
  }
  
  
  
  selected_columns <- selected_features %>%
    dplyr::distinct(cn, CompoundName, feature_set, feature) %>% 
    dplyr::mutate(tar = ifelse(feature_set %in% c("CRISPR", "RNAi"), paste0(feature_set, "_"), 
                               ifelse(feature_set == "Expression",  "EXP_", 
                                      ifelse(feature_set == "CopyNumber",  "CN_" , 
                                             ifelse(feature_set == "Mutation",  "MUT_" ,
                                                    ifelse(feature_set == "Fusion",   "FUS_", 
                                                           ifelse(feature_set == "Lineage", "LIN_" , NA))))))) %>% 
    dplyr::filter(!is.na(tar)) %>% 
    dplyr::mutate(y = paste0(tar, feature)) 
  
  
  fit <- function(x,y){
    require(ranger)
    cl <- intersect(train, names(y))
    cl.test <- intersect(test, names(y))
    
    rf <- ranger::ranger(x = x[cl, , drop = F], y = y[cl] , importance = "impurity")
    
    pr <- predict(rf, data = x[union(cl, cl.test), , drop = F])
    
    y.hat <- tibble(depmap_id = union(cl, cl.test), 
                    y.hat = pr$predictions, y = y[union(cl, cl.test)],
                    type = ifelse(depmap_id %in% cl, "train", "test")) %>% 
      dplyr::left_join(tibble(depmap_id = cl, y.hat.oob = rf$predictions))
    
    
    imp <- tibble(var = names(rf$variable.importance), imp = rf$variable.importance) %>%
      dplyr::arrange(desc(imp))
    
    res <- tibble(mse.oob = mean((rf$predictions - y[cl])^2, na.rm = T),
                  var.y.train = var(y[cl], na.rm = T), 
                  r2.oob = mse.oob/var.y.train,
                  r.oob = cor(rf$predictions, y[cl], use = "p")[,1])
    
    if(!is.null(test)){
      res <- y.hat %>% 
        dplyr::filter(type == "test") %>% 
        dplyr::summarise(var.y.test = var(y, na.rm = T),
                         mse = mean((y - y.hat)^2, na.rm = T),
                         r2 = 1- mse/var.y.test,
                         r = cor(y, y.hat, use = "p")[,1]) %>%
        bind_cols(res)
    }
    
    
    return(list(res, y.hat, imp))
  }
  
  
  biomarker_table <- list(); prediction_table <- list(); importance_table <- list(); jx <- 1
  for(cmp in colnames(Y)){
    tars <- dplyr::filter(targets, cn == cmp)$GeneSymbolOfTargets %>% unique()
    extras <- dplyr::filter(selected_columns , cn == cmp)$y %>% intersect(colnames(X))
    
    y <- Y[, cmp]; y <- y[is.finite(y)]
    res <- list();    pred <- list(); imp <- list(); ix <- 1 
    
    if(length(tars) > 0){
      # fit a model for each target
      for(tar in tars){
        
        x <- X[names(y), unique(dplyr::filter(targets, cn == cmp, GeneSymbolOfTargets == tar)$cn.feat), drop = F]
        temp <- fit(x,y)
        
        res[[ix]] <-  temp[[1]] %>% 
          dplyr::mutate(model = tar)
        
        pred[[ix]] <- temp[[2]] %>% 
          dplyr::mutate(model = tar)
        
        imp[[ix]] <- temp[[3]] %>% 
          dplyr::mutate(model = tar)
        
        ix <- ix + 1
      }
      
      # fit all the targets together
      x <- X[names(y), unique(dplyr::filter(targets, cn == cmp, GeneSymbolOfTargets %in% tars)$cn.feat), drop = F]
      temp <- fit(x,y)
      
      res[[ix]] <-  temp[[1]] %>% 
        dplyr::mutate(model = "targets")
      
      pred[[ix]] <- temp[[2]] %>% 
        dplyr::mutate(model = "targets")
      
      imp[[ix]] <- temp[[3]] %>% 
        dplyr::mutate(model = "targets")
      
      ix <- ix + 1
    }
    
    cols <- unique(union(dplyr::filter(targets, cn == cmp, GeneSymbolOfTargets %in% tars)$cn.feat, extras))
    if(length(cols) > 0){
      # fit the extended model
      x <- X[names(y), cols, drop = F]
      
      temp <- fit(x,y)
      
      res[[ix]] <-  temp[[1]] %>% 
        dplyr::mutate(model = "extended")
      
      pred[[ix]] <- temp[[2]] %>% 
        dplyr::mutate(model = "extended")
      
      imp[[ix]] <- temp[[3]] %>% 
        dplyr::mutate(model = "extended")
      
      
    }
    
    # put them all together
    biomarker_table[[jx]] <- res %>% 
      dplyr::bind_rows() %>% 
      dplyr::mutate(cn = cmp)
    
    prediction_table[[jx]] <- pred %>% 
      dplyr::bind_rows() %>% 
      dplyr::mutate(cn = cmp)
    
    importance_table[[jx]] <- imp %>% 
      dplyr::bind_rows() %>% 
      dplyr::mutate(cn = cmp)
    
    print(paste0(cmp, " - ", jx))
    jx <- jx + 1
  }
  
  return(list(dplyr::bind_rows(biomarker_table) %>% 
                dplyr::left_join(CompoundList %>% dplyr::distinct(cn, CompoundName)), 
              dplyr::bind_rows(prediction_table) %>% 
                dplyr::left_join(CompoundList %>% dplyr::distinct(cn, CompoundName)), 
              dplyr::bind_rows(importance_table) %>% 
                dplyr::left_join(CompoundList %>% dplyr::distinct(cn, CompoundName)), 
              bm.auc))
}


biomarker_suite_rf_cv <- function(X, Y, biomarker_file, CompoundList, bm_th = 0.05, bm_R = 10, bm_R2 = 50, K = 10, seed = NULL){
  require(tidyverse)
  
  if(!is.null(seed)) set.seed(seed)
  
  cl <- intersect(rownames(X), rownames(Y)) %>% sample()
  RES <- list(); PRED <- list(); IMP <- list()
  for(k in 1:K){
    print(k)
    temp <- biomarker_suite_rf(X, Y, biomarker_file = biomarker_file, CompoundList = CompoundList, test_samples = cl[seq.int(k, length(cl), by = K)], bm_th = bm_th, bm_R = bm_R, bm_R2 = bm_R2)
    
    RES[[k]] <- temp[[1]] %>% 
      dplyr::mutate(K = k)
    
    PRED[[k]] <- temp[[2]] %>% 
      dplyr::mutate(K = k)
    
    IMP[[k]] <- temp[[3]] %>% 
      dplyr::mutate(K = k)
  }
  
  temp <- biomarker_suite_rf(X, Y, biomarker_file = biomarker_file, CompoundList = CompoundList, test_samples = NULL, bm_th = bm_th, bm_R = bm_R, bm_R2 = bm_R2)
  
  RES[[K + 1]] <- temp[[1]] %>% 
    dplyr::mutate(K = 0)
  
  PRED[[K + 1]] <- temp[[2]] %>% 
    dplyr::mutate(K = 0)
  
  IMP[[K + 1]] <- temp[[3]] %>% 
    dplyr::mutate(K = 0)
  
  
  return(list(model_performances = dplyr::bind_rows(RES), predictions = dplyr::bind_rows(PRED),  variable_importances = dplyr::bind_rows(IMP), univariate_biomarkers = temp[[4]]))
}





# Data processing -----


#' Fitting a dose response curve to the given dose and viability (FC) values.
#' This function fits 5 different dose-response functiosn to the given dose, viability pairs using dr4pl and drc packages and returns
#' the best one (lowest mse) among them.
#'
#' @param FC : Measured viability vector
#' @param dose : Dose vector corresponding to FC
#' @param UL_low : Lower limit for the upper asymptote
#' @param UL_up : Upper limit for the upper asympotote
#' @param slope_decreasing: Should the curve to be constrained to be decreasing or not.
#'
#' @return Returns a single row data-frame with following columns:
#'          fit_name : Name of the fitting method with the highest explained variance (lowest mse)
#'          lower_limit : Lower asmpytote for the selected fit
#'          upper_limit : Upper asymptote for the selected fit
#'          slope : The Hill slope for the selected fit
#'          inflection : inflection point of the selected fit (EC50)
#'          mse : mse of the selected fit
#'          mad : mad of the selected fit
#'          frac_var_explained : Explained variance for the best fit
#'          successful_fit: If any of the fitting methods has positive frac_var_explained
#'          auc_riemann : The average measured viability value across doses
#'          minimum_dose : Minimum measured dose
#'          maximum_dose : Maximum measured dose
#'          auc : auc of the log-dose vs viability curve normalized to the measured dose-range (takes values between 0 and 1)
#'          log2_ic50 : Log2 of IC50 for the fitted curve
#' @export
#'
#' @examples
get_best_fit <- function(FC, dose, UL_low=0.8, UL_up=1.01, slope_decreasing=TRUE, seed = NULL) {
  require(dr4pl)
  require(drc)
  require(tidyverse)
  require(magrittr)
  
  if(!is.null(seed)) set.seed(seed)
  
  # Fits a number of alternate models  to the DRC and chooses the best fit.
  
  # UL low is the lowerbound of UL we pass to the optimizer and UL_up is the upper bound of UL that we pass to the optimizer
  # fomat of output will be:-
  # results.df <- data.frame("fit_name"=character(),
  #                          "lower_limit"=double(),
  #                          "upper_limit"=double(),
  #                          "slope"=double(),
  #                          "inflection"=double(),
  #                          "mse"=double(), "mad" =double(),
  #                          "frac_var_explained"=double())
  
  dose = as.numeric(dose)
  FC = as.numeric(FC)
  slope_bound <- ifelse(slope_decreasing, 1e-5, Inf)  # bound the slopes by default unless passed another option
  riemann_auc <- mean(pmin(1,FC)) ## mean fold-change after rounding FC to 1.
  var_data = var(FC)
  
  minimum_dose = min(dose); maximum_dose = max(dose)
  
  results.df <- list(); ix = 1
  
  
  # FIT 1 ---
  drc_model <-  tryCatch(drc::drm(FC ~ dose, data= data.frame(FC = FC, dose = dose),
                                  fct=LL.4(names = c("slope", "Lower Limit", "Upper Limit", "ED50")),
                                  lowerl = c(-slope_bound,0.0, UL_low, -Inf),upperl = c(Inf,1.01,UL_up, Inf)),
                         error = function(e)
                         {return(list(convergence=FALSE, error=TRUE,fit=list(convergence=FALSE)))})
  # "slope" in drc package is -ve of slope in dr4pl package
  
  
  if (drc_model$fit$convergence & all(is.finite(unlist(drc_model$coefficients)))){
    mse_mad <- compute_mse_mad(FC, dose, as.numeric(drc_model$coefficients[[3]]), as.numeric(drc_model$coefficients[[2]]),
                               -as.numeric(drc_model$coefficients[[1]]), as.numeric(drc_model$coefficients[[4]]))
    # "slope" in drc package is -ve of slope in dr4pl package and so -ve sign needs to be put in here.
    
    results.df[[ix]] <- tibble( fit_name = "drc_drm_constrained",
                                lower_limit = as.numeric(drc_model$coefficients[[2]]),
                                upper_limit = as.numeric(drc_model$coefficients[[3]]),
                                slope = -as.numeric(drc_model$coefficients[[1]]),
                                inflection = as.numeric(drc_model$coefficients[[4]]),
                                mse = mse_mad$mse, mad = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
    ix = ix + 1
  }
  
  
  
  # FIT 2 ---
  drc_model <-  tryCatch(drc::drm(FC ~ dose, data= data.frame(FC = FC, dose = dose),
                                  fct=LL.4(names = c("slope", "Lower Limit", "Upper Limit", "ED50")),
                                  lowerl = c(-Inf,0.0, UL_low, -Inf),upperl = c(Inf,1.01,UL_up, Inf)),
                         error = function(e)
                         {return(list(convergence=FALSE, error=TRUE,fit=list(convergence=FALSE)))})
  # "slope" in drc package is -ve of slope in dr4pl package
  
  
  if (drc_model$fit$convergence & all(is.finite(unlist(drc_model$coefficients)))){
    if((!slope_decreasing) | (as.numeric(drc_model$coefficients[[1]]) > 0)){
      mse_mad <- compute_mse_mad(FC, dose, as.numeric(drc_model$coefficients[[3]]), as.numeric(drc_model$coefficients[[2]]),
                                 -as.numeric(drc_model$coefficients[[1]]), as.numeric(drc_model$coefficients[[4]]))
      # "slope" in drc package is -ve of slope in dr4pl package and so -ve sign needs to be put in here.
      
      results.df[[ix]] <- tibble( fit_name = "drc_drm_unconstrained",
                                  lower_limit = as.numeric(drc_model$coefficients[[2]]),
                                  upper_limit = as.numeric(drc_model$coefficients[[3]]),
                                  slope = -as.numeric(drc_model$coefficients[[1]]),
                                  inflection = as.numeric(drc_model$coefficients[[4]]),
                                  mse = mse_mad$mse, mad = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
      ix = ix + 1
    }
  }
  
  
  # FIT 3 ---
  dr4pl_initMan_optNM <- tryCatch(dr4pl(dose, FC,
                                        init.parm = dr4pl::dr4pl_theta(theta_1 = 1, theta_2 = 8*min(dose),
                                                                       theta_3= -3, theta_4 = 0.01),
                                        lowerl = c(UL_low, -Inf, -Inf, 0),
                                        upperl = c(UL_up, Inf, slope_bound, 1.01),
                                        method.optim="Nelder-Mead"),
                                  error= function(e){return(list(convergence=FALSE, error=TRUE))}
  )
  
  if (!dr4pl_initMan_optNM$convergence){
    if (!is.null(dr4pl_initMan_optNM$dr4pl.robust)) {
      dr4pl_initMan_optNM <- dr4pl_initMan_optNM$dr4pl.robust
    }
  }
  
  if (dr4pl_initMan_optNM$convergence & all(is.finite(unlist(dr4pl_initMan_optNM$parameters)))){
    mse_mad <- compute_mse_mad(FC, dose, dr4pl_initMan_optNM$parameters[[1]], dr4pl_initMan_optNM$parameters[[4]],
                               dr4pl_initMan_optNM$parameters[[3]], dr4pl_initMan_optNM$parameters[[2]])
    
    results.df[[ix]] <- tibble( fit_name = "dr4pl_initMan_constrained_optNM",
                                lower_limit = as.numeric(dr4pl_initMan_optNM$parameters[[4]]),
                                upper_limit = as.numeric(dr4pl_initMan_optNM$parameters[[1]]),
                                slope = as.numeric(dr4pl_initMan_optNM$parameters[[3]]),
                                inflection = as.numeric(dr4pl_initMan_optNM$parameters[[2]]),
                                mse = mse_mad$mse, mad = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
    ix = ix + 1
  }
  
  # FIT 4 ---
  dr4pl_unconstrained <- tryCatch(dr4pl(dose, FC,
                                        init.parm = dr4pl::dr4pl_theta(theta_1 = 1, theta_4 = 0.3),
                                        method.init = "logistic",
                                        lowerl = c(0.99, -Inf, -Inf, 0),
                                        upperl = c(1.01, Inf, Inf, 1.01)),
                                  error = function(e) {print(e); return(NA)})
  
  if (!all(is.na(dr4pl_unconstrained))) {
    if (!dr4pl_unconstrained$convergence) {
      dr4pl_unconstrained <- dr4pl_unconstrained$dr4pl.robust
    }
  }
  
  
  param <- tryCatch(dr4pl_unconstrained$parameters, error = function(e) return(NA))
  if (!all(is.na(param))){
    if(as.numeric(dr4pl_unconstrained$parameters[[3]])<slope_bound){ ### while slope bound is not passed to this last optimizer, we do not accept a solution not within the bound
      mse_mad <- compute_mse_mad(FC, dose, dr4pl_unconstrained$parameters[[1]], dr4pl_unconstrained$parameters[[4]],
                                 dr4pl_unconstrained$parameters[[3]], dr4pl_unconstrained$parameters[[2]])
      results.df[[ix]] <- tibble( fit_name = "dr4pl_initL_unconstrained",
                                  lower_limit = as.numeric(dr4pl_unconstrained$parameters[[4]]),
                                  upper_limit = as.numeric(dr4pl_unconstrained$parameters[[1]]),
                                  slope = as.numeric(dr4pl_unconstrained$parameters[[3]]),
                                  inflection = as.numeric(dr4pl_unconstrained$parameters[[2]]),
                                  mse = mse_mad$mse, mad = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
      ix = ix + 1
    }
  }
  
  # FIT 5 ---
  dr4pl_initL <- tryCatch(dr4pl(dose, FC,
                                init.parm = dr4pl::dr4pl_theta(theta_1 = 1, theta_4 = 0.005),
                                method.init = "logistic",
                                lowerl = c(UL_low, -Inf, -Inf, 0),
                                upperl = c(UL_up, Inf, slope_bound, 1.01)),
                          error= function(e){return(list(convergence=FALSE, error=TRUE))}
  )
  
  if (!dr4pl_initL$convergence){
    if (!is.null(dr4pl_initL$dr4pl.robust)) {
      dr4pl_initL <- dr4pl_initL$dr4pl.robust
    }
  }
  
  if (dr4pl_initL$convergence & all(is.finite(unlist(dr4pl_initL$parameters)))){
    mse_mad <- compute_mse_mad(FC,dose, dr4pl_initL$parameters[[1]], dr4pl_initL$parameters[[4]],
                               dr4pl_initL$parameters[[3]], dr4pl_initL$parameters[[2]])
    
    results.df[[ix]] <- tibble( fit_name = "dr4pl_initL_constrained",
                                lower_limit = as.numeric(dr4pl_initL$parameters[[4]]),
                                upper_limit = as.numeric(dr4pl_initL$parameters[[1]]),
                                slope = as.numeric(dr4pl_initL$parameters[[3]]),
                                inflection = as.numeric(dr4pl_initL$parameters[[2]]),
                                mse = mse_mad$mse, mad = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
    ix = ix + 1
  }
  
  
  # Choose the best fit among the successful fits---
  results.df <- dplyr::bind_rows(results.df)
  
  if (nrow(results.df)>0){
    results.df <- results.df %>%
      dplyr::arrange(desc(frac_var_explained)) %>%
      head(1) %>%
      dplyr::mutate(successful_fit = TRUE,
                    auc_riemann = as.numeric(riemann_auc),
                    minimum_dose = minimum_dose, maximum_dose = maximum_dose) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(auc = compute_auc(lower_limit, upper_limit, inflection, slope, minimum_dose, maximum_dose),
                    log2_auc = log2(auc),
                    log2_ic50 = compute_log_ic50(lower_limit, upper_limit, inflection, slope, minimum_dose, maximum_dose))
    
  }else{
    results.df  <- data.frame(successful_fit=FALSE,
                              auc_riemann = as.numeric(riemann_auc),
                              minimum_dose = minimum_dose, maximum_dose = maximum_dose, auc = NA, log2_ic50 = NA)
  }
  
  return (results.df)
}



#' Computing area under a 4 parameter log-logistic dose response curve
#'
#' @param LL : Lower asymptote
#' @param UL : Upper asymptote
#' @param inflection : inflection point (EC50)
#' @param slope : Hill slope ( > 0 for decreasing curves)
#' @param minimum_dose : Minimum dose
#' @param maximum_dose : Maximum dose
#'
#' @return auc value of the dose-response function (x: log2(dose), y: response) scaled to the dose range.
#' @export
#'
#' @examples
compute_auc <- function(LL, UL, inflection, slope, minimum_dose, maximum_dose) {
  f1 = function(x) pmax(pmin((UL + (LL - UL)/(1 + (2^x/inflection)^slope)), 1, na.rm = T), 0, na.rm = T)
  return(tryCatch(integrate(f1, log2(minimum_dose), log2(maximum_dose))$value/(log2(maximum_dose/minimum_dose)),
                  error = function(e) {print(e); NA}))
}


#' Computing area under a 4 parameter log-logistic dose response curve after tanh scaling
#'
#' @param LL : Lower asymptote
#' @param UL : Upper asymptote
#' @param inflection : inflection point (EC50)
#' @param slope : Hill slope ( > 0 for decreasing curves)
#' @param minimum_dose : Minimum dose
#' @param maximum_dose : Maximum dose
#' @param m : mean parameter for the transformation
#' @param s : spread parameter for the transformation
#'
#' @return auc value of the dose-response function (x: log2(dose), y: response) scaled to the dose range.
#' @export
#'
#' @examples
compute_scaled_auc <- function(LL, UL, inflection, slope, minimum_dose, maximum_dose, m, s) {
  
  f1 = function(x){
    y = (UL + (LL - UL)/(1 + (2^x/inflection)^slope))
    z = pmax(pmin((1/2 + tanh(atanh(0.5) * (y - m)/s)/2), 1, na.rm = T), 0, na.rm = T)
    return(z)
  } 
  return(tryCatch(integrate(f1, log2(minimum_dose), log2(maximum_dose))$value/(log2(maximum_dose/minimum_dose)),
                  error = function(e) {print(e); NA}))
}


#' Computing IC50 for a 4 parameter log-logistic dose response curve
#'
#' @param LL : Lower asymptote
#' @param UL : Upper asymptote
#' @param inflection : inflection point (EC50)
#' @param slope : Hill slope ( > 0 for decreasing curves)
#' @param minimum_dose : Minimum dose
#' @param maximum_dose : Maximum dose
#'
#' @return IC50 : The dose value where the curve intersects with y = 0.5, NA returned if they don't intersect in the given dose range.
#' @export
#'
#' @examples
compute_log_ic50 <- function(LL, UL, inflection, slope, minimum_dose, maximum_dose) {
  if((LL >= 0.5) | (UL <= 0.5)) {
    return(NA)
  } else {
    f1 = function(x) (UL + (LL - UL)/(1 + (2^x/inflection)^slope)- 0.5)
    return(tryCatch(uniroot(f1, c(log2(minimum_dose), log2(maximum_dose)))$root,  error = function(x) NA))
  }
}

#' Computing mean square and median absolute errors for for a 4 parameter log-logistic dose response curve and the corresponding (dose,viability) pairs.
#'
#' @param FC : Measured viability vector
#' @param dose : Dose vector corresponding to FC
#' @param UL : Upper asymptote
#' @param LL : Lower asymptote
#' @param slope : Hill slope
#' @param inflection : inflection point (EC50)
#'
#' @return List of mse and mad values.
#' @export
#'
#' @examples
compute_mse_mad <- function(FC, dose,  UL, LL,  slope, inflection) {
  FC.pred = UL  + (LL -UL )/(1 + (dose/inflection)^slope)
  residuals = FC - FC.pred
  return(list(mse = mean(residuals^2), mad = median(abs(residuals))))
}




