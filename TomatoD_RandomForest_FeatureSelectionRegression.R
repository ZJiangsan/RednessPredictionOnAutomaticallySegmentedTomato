

###  read the extracted HSI information of segmented tomato, do feature selection and random forest regression based on selected features


setwd("/Users/zhaojiangsan/Documents/RScript")
getwd()

## calculate NDVI of constructed tomato HSIs
##       for "file 18"
library(rhdf5)
library(baseline) # Collection of baseline correction algorithms

normalize = function(x) {(x-min(x))/(max(x)-min(x))}


tomato_spec_hsi='/Users/zhaojiangsan/Documents/RScript/tomato_hsi_mat' #out_hsi_reconstruction

files_spec= list.files(tomato_spec_hsi)#[-c(2,3,9,10)] #[c(1,8:12)]


wi=read.csv("wavelength_interval.csv") # here is the wavebands of specim IQ hyperspectral camera, 204 bands from 397 to 1003 nm
spec_color=wi$w_i
length(spec_color)


## extract NDVI of each date set
file_spec_ids = list.files(paste(tomato_spec_hsi, files_spec[2], sep = "/"))

## for tomato D
## need to adapt the id to treatments
#  the order of tomato in segmented image is manually registered here
id_2_od = data.frame(id_1_4 = c(1,11,5,8,10,7,6,12,3,4,9,2), 
                     id_5_8 = c(5,4,6,12,2,9,7,8,10,11,1,3),
                     id_9_12 = c(7,8,4,5,6,9,12,1,10,2,11,3),
                     od = c(1:12))

## extract the median values of the tomato overal spectral reflectance
hsi_spec_D_median = data.frame()
for (m in c(1:length(file_spec_ids))){
  print(paste("m", m, sep = ":"))
  tomato_spec_hsi_var = paste(tomato_spec_hsi, files_spec[2], file_spec_ids[m], sep = "/")
  hsi_spec_D_i <- data.frame(h5read(tomato_spec_hsi_var,"img"))
  hsi_spec_D_i[,1] = hsi_spec_D_i[,1]+1
  
  ## read the chemical information of each tomato
  id_2_cod = c(rep(0, dim(hsi_spec_D_i)[1]))
  
  for (i in c(1:12)){
    id_2_cod[hsi_spec_D_i[,1]==id_2_od[i,m]] = id_2_od[i,4]
  }
  
  ## change the original ID to color orders
  hsi_spec_D_i[,1] =  id_2_cod
  
  ## do shape correction: Asymmetric least squares (AsLs) baseline correction of logarithmic linearised reflectance log(1/R) 
  if (min(hsi_spec_D_i[,-1])<=0){ # in case something is divided by 0
    hsi_spec_D_i[,-1][hsi_spec_D_i[,-1]==min(hsi_spec_D_i[,-1])] =0.0001
  }
  
  hsi_spec_D_i_c= baseline.als(as.matrix(log(1/hsi_spec_D_i[,-1])), lambda = 3, p = 0.001, maxit = 20)$corrected # lambda controls smoothness, p for baseline?
  
  # get the median of each tomato
  hsi_spec_D_i_median = aggregate(hsi_spec_D_i_c, list(hsi_spec_D_i$X1), median) # median
  
  hsi_spec_D_i_median$Group.1 = c(((m-1)*12 +1): (m*12))
  
  size_id = 12/length(file_spec_ids)
  hsi_spec_D_i_median = data.frame(Group.2 = rep(((m-1)*size_id +1) : (m*size_id), each =length(file_spec_ids)), hsi_spec_D_i_median)
  print(paste("id", paste(hsi_spec_D_i_median$Group.1, collapse = ","), sep = ": "))
  hsi_spec_D_median= rbind(hsi_spec_D_median, hsi_spec_D_i_median)
}


## save the date in mat format
# H5close() # show what is inside the data
# h5createFile("tomato_color_reflectance_median.h5") # first create the file (only create the file in the first time)
# h5write(hsi_spec_D_median, file = "tomato_color_reflectance_median.h5", name="img") # write the data in
# shit_a <- h5read("tomato_color_reflectance_median.h5","img") # read the data

names(hsi_spec_D_median) = c("cd", "id", paste("spec", c(1:204), sep = "_"  ) ) # rename variables as wavebands

## save the data in csv

write.csv(hsi_spec_D_median, "hsi_spec_D_median_reflectance_cr_lambda1_p0.001.csv", row.names = F)
# write.csv(hsi_spec_D_median, "hsi_spec_D_median_reflectance.csv", row.names = F)

####
hsi_spec_D_median = read.csv("hsi_spec_D_median_reflectance_cr_lambda1_p0.001.csv")
hsi_spec_D_median[, 1:5]

## READ TOMATO CHEMICS
tomato_D_chemics = read.csv("tomato_D_chemics_NAI.csv") # tomato chemical properties measured


## combine chemical measurements and spectral information

tomato_D_spec_chemics = data.frame(tomato_D_chemics, hsi_spec_D_median)

tomato_D_spec_chemics[, 1:10]



###      a complete case to do feature selection, model training, prediction and plotting in one time
library(Metrics)
rsq <- function (x, y) cor(x, y) ^ 2

##
# load the library
library(mlbench)
library(caret)
library(plsVarSel)
library(pls)
library(ade4)
library(randomForest)
############################3
## remove NAs
tomato_D_spec_chemics_c = na.omit(tomato_D_spec_chemics)[-c(1, 6,7)]
tomato_D_spec_chemics_c[1:10]
names(tomato_D_spec_chemics_c)[1:4]

##  do repeated feature selection
subsets <- c(seq(1,100,1), seq(101,204,2)) # set length of variables to sample
set.seed(123)
seeds <- vector(mode = "list", length = 51) # set seed for each resampling iteration 50+1 = repeats *number +1
for(i in 1:50) seeds[[i]] <- sample.int(1000, length(subsets) + 1)  #repeats *number
seeds[[51]] <- sample.int(1000, 1) # the last one
##
library(doParallel) 
# ## check how many cores are there
detectCores(all.tests = FALSE, logical = TRUE)

tomato_D_var_rmse_r2_5s = data.frame()

tomato_D_pred = data.frame(ind = c(1:33))

for (i in c(1:4)){ # four labels to be predicted
  print(names(tomato_D_spec_chemics_c)[i])
  
  cl <- makeCluster(8)
  registerDoParallel(cl)
  ##
  
  #set.seed(1)
  system.time(RF_model_rfe_i <- rfe(tomato_D_spec_chemics_c[,-c(1:4)], tomato_D_spec_chemics_c[,i],
                                    sizes = subsets,
                                    rfeControl = rfeControl(functions = rfFuncs,
                                                            method = "repeatedcv",
                                                            number=10, # number of folds or rasampling iterations
                                                            repeats = 5, #the number of complete sets of folds to compute
                                                            
                                                            seeds = seeds,
                                                            verbose = TRUE)))
  stopCluster(cl)
  
  # RF_model_rfe_i$fit
  
  ## create a model with top five important variables
  for(x in c(1:length(predictors(RF_model_rfe_i)))){
    print(paste(x, "out of", length(predictors(RF_model_rfe_i)), sep = " "))
    # choose x  predictors to build the final model
    tomato_D_spec_chemics_c_i = data.frame(tomato_D_spec_chemics_c[i], tomato_D_spec_chemics_c[as.character(predictors(RF_model_rfe_i)[1:x])])
    # tomato_D_spec_chemics_c_i
    # use loov here
    RF_model_pred_v_i =c()
    for (j in c(1:dim(tomato_D_spec_chemics_c_i)[1])){
      set.seed(5433)
      if( x ==1) {
        tomato_D_spec_chemics_c_i_j = data.frame(tomato_D_spec_chemics_c_i[-j, -1])
        names(tomato_D_spec_chemics_c_i_j) = names(tomato_D_spec_chemics_c_i)[2]
      } else {
        tomato_D_spec_chemics_c_i_j = tomato_D_spec_chemics_c_i[-j, -1]
      }
      
      RF_model_final_i_j = randomForest(tomato_D_spec_chemics_c_i_j, tomato_D_spec_chemics_c_i[-j,1])
      # RF_model_final_i
      
      ## make predictions
      RF_model_pred_i_j=predict(RF_model_final_i_j, tomato_D_spec_chemics_c_i[j,-1])
      
      RF_model_pred_v_i = append(RF_model_pred_v_i, unname(RF_model_pred_i_j))
    }
    
    # save the predictions
    RF_model_pred_v_i_d = data.frame(RF_model_pred_v_i)
    names(RF_model_pred_v_i_d) =paste(names(tomato_D_spec_chemics_c_i)[1],x, sep = "_")
    
    tomato_D_pred = cbind(tomato_D_pred, RF_model_pred_v_i_d)
    
    rmse = rmse(tomato_D_spec_chemics_c_i[,1], RF_model_pred_v_i)
    
    
    # cord_min = min(c(tomato_D_spec_chemics_c_i[,1], RF_model_pred_v_i ))*1.0
    # cord_max = max(c(tomato_D_spec_chemics_c_i[,1], RF_model_pred_v_i))*1.0
    
    # x11()
    # plot(x= tomato_D_spec_chemics_c_i[,1], y =RF_model_pred_v_i,  xlim = c(cord_min, cord_max), ylim =  c(cord_min, cord_max))
    # abline(glm(tomato_D_spec_chemics_c_i[,1]~RF_model_pred_v_i), col="red")
    
    ###########################
    df_i=data.frame(pred=RF_model_pred_v_i, actual=tomato_D_spec_chemics_c_i[,1],ind=c(1:length(tomato_D_spec_chemics_c_i[,1])))
    


    # evaluate the quality of prediction model
    lm_i = lm(RF_model_pred_v_i ~ tomato_D_spec_chemics_c[,i])
    
    p_value_i = format(summary(lm_i)$coefficients[2,4], digits = 3) 
    r2 = round(summary(lm_i)$r.squared, digits = 2)
    
    
    ## save the usefult imformation
    tomato_D_var_rmse_r2_5s_i = data.frame(target = names(tomato_D_spec_chemics_c_i)[1], 
                                           r2= r2,
                                           rmse = rmse,
                                           # predictors = paste(predictors(RF_model_rfe_i)[1:5], collapse = ","),
                                           predictors = paste(names(tomato_D_spec_chemics_c_i)[-1], collapse = ","),
                                           no. = x)
    
    tomato_D_var_rmse_r2_5s = rbind(tomato_D_var_rmse_r2_5s, tomato_D_var_rmse_r2_5s_i)
  }
}



# tomato_D_var_rmse_r2_5s 
# tomato_D_var_rmse_r2_5s[tomato_D_var_rmse_r2_5s$r2 == max(tomato_D_var_rmse_r2_5s$r2),]

# tomato_D_pred 

##  extract the best model
tomato_D_var_r2_max = data.frame()
for (i in c(as.character(unique(tomato_D_var_rmse_r2_5s$target)))){
  t_i = tomato_D_var_rmse_r2_5s[tomato_D_var_rmse_r2_5s$target ==i,]
  t_i_max = t_i[t_i$r2 == max(t_i$r2),]
  # if there are multiple predictors lengths that lead to the same r2, then choose the shortest length
  if( dim(t_i_max)[1] >1) t_i_max_x = t_i_max[t_i_max$no. == min(t_i_max$no.),][1,] else t_i_max_x = t_i_max
  
  tomato_D_var_r2_max = rbind(tomato_D_var_r2_max, t_i_max_x)
}

tomato_D_var_r2_max


################################ if you want, can do furthre feature selection from the already selected features, can be skipped


tomato_D_imp_rmse_r2_5s = data.frame()

tomato_D_imp_pred = data.frame(ind = c(1:33))
for (i in c(as.character(unique(tomato_D_var_r2_max$target)))){
  print(i)
  
  ## first run, used the originally selected important variable
  pred_imp = unlist(strsplit(as.character(tomato_D_var_r2_max[tomato_D_var_r2_max$target == i,]$predictors), ","))
  
  # from the second run, use the newly selected important variables
  # pred_imp = unlist(strsplit(as.character(tomato_D_imp_r2_max[tomato_D_imp_r2_max$target == i,]$predictors), ","))
  
  
  if (length(pred_imp) <=2) { # if the number of important variables is less than 5, then skip feature selection
    ## put the same features in the new data frame
    tomato_D_imp_rmse_r2_5s_i = tomato_D_imp_r2_max[tomato_D_imp_r2_max$target == i,]
    # tomato_D_imp_rmse_r2_5s = rbind(tomato_D_imp_rmse_r2_5s, tomato_D_var_rmse_r2_5s_i) # first time
    tomato_D_imp_rmse_r2_5s = rbind(tomato_D_imp_rmse_r2_5s, tomato_D_imp_rmse_r2_5s_i) # second time
    
    next
  }
  
  print(pred_imp)
  tomato_D_spec_chemics_c_imp = data.frame(tomato_D_spec_chemics_c[i], tomato_D_spec_chemics_c[pred_imp])
  
  
  subsets <- c(seq(1,length(pred_imp) ,1))
  set.seed(123)
  seeds <- vector(mode = "list", length = 51) # set seed for each resampling iteration 50+1 = repeats *number +1
  for(i in 1:50) seeds[[i]] <- sample.int(1000, length(subsets) + 1)  #repeats *number
  seeds[[51]] <- sample.int(1000, 1) # the last one
  ##
  library(doParallel) 
  # ## check how many cores are there
  detectCores(all.tests = FALSE, logical = TRUE)
  
  cl <- makeCluster(8)
  registerDoParallel(cl)
  ##
  
  #set.seed(1)
  system.time(RF_model_rfe_i <- rfe(tomato_D_spec_chemics_c_imp[,-1], tomato_D_spec_chemics_c_imp[,1],
                                    sizes = subsets,
                                    rfeControl = rfeControl(functions = rfFuncs,
                                                            method = "repeatedcv",
                                                            number=10, # number of folds or rasampling iterations
                                                            repeats = 5, #the number of complete sets of folds to compute
                                                            seeds = seeds,
                                                            verbose = TRUE)))
  stopCluster(cl)
  
  # RF_model_rfe_i$fit
  
  ## create a model with top five important variables
  for(x in c(1:length(predictors(RF_model_rfe_i)))){
    print(paste(x, "out of", length(predictors(RF_model_rfe_i)), sep = " "))
    # choose x  predictors to build the final model
    tomato_D_spec_chemics_c_imp_x = data.frame(tomato_D_spec_chemics_c_imp[1], tomato_D_spec_chemics_c_imp[as.character(predictors(RF_model_rfe_i)[1:x])])
    # tomato_D_spec_chemics_c_i
    # use loov here
    RF_model_pred_v_i =c()
    for (j in c(1:dim(tomato_D_spec_chemics_c_imp_x)[1])){
      set.seed(5433)
      if( x ==1) {
        tomato_D_spec_chemics_c_imp_x_j = data.frame(tomato_D_spec_chemics_c_imp_x[-j, -1])
        names(tomato_D_spec_chemics_c_imp_x_j) = names(tomato_D_spec_chemics_c_imp_x)[2]
      } else {
        tomato_D_spec_chemics_c_imp_x_j = tomato_D_spec_chemics_c_imp_x[-j, -1]
      }
      
      RF_model_final_i_j = randomForest(tomato_D_spec_chemics_c_imp_x_j, tomato_D_spec_chemics_c_imp_x[-j,1])
      # RF_model_final_i
      
      ## make predictions
      RF_model_pred_i_j=predict(RF_model_final_i_j, tomato_D_spec_chemics_c_imp_x[j,-1])
      
      RF_model_pred_v_i = append(RF_model_pred_v_i, unname(RF_model_pred_i_j))
    }
    
    # save the predictions
    RF_model_pred_v_i_d = data.frame(RF_model_pred_v_i)
    names(RF_model_pred_v_i_d) =paste(names(tomato_D_spec_chemics_c_imp_x)[1],x, sep = "_")
    
    tomato_D_imp_pred = cbind(tomato_D_imp_pred, RF_model_pred_v_i_d)
    
    rmse = rmse(tomato_D_spec_chemics_c_imp_x[,1], RF_model_pred_v_i)
    
    
    ###########################
    df_i=data.frame(pred=RF_model_pred_v_i, actual=tomato_D_spec_chemics_c_imp_x[,1],ind=c(1:length(tomato_D_spec_chemics_c_imp_x[,1])))
    
    ###
    # r2 = round(rsq(tomato_D_spec_chemics_c_i[,1],RF_model_pred_v_i), 2)
    # evaluate the quality of prediction model
    lm_i = lm(RF_model_pred_v_i ~ tomato_D_spec_chemics_c_imp_x[,1])
    
    p_value_i = format(summary(lm_i)$coefficients[2,4], digits = 3) 
    r2 = round(summary(lm_i)$r.squared, digits = 2)
    
    
    ## save the usefult imformation
    tomato_D_imp_rmse_r2_5s_i = data.frame(target = names(tomato_D_spec_chemics_c_imp_x)[1], 
                                           r2= r2,
                                           rmse = rmse,
                                           # predictors = paste(predictors(RF_model_rfe_i)[1:5], collapse = ","),
                                           predictors = paste(names(tomato_D_spec_chemics_c_imp_x)[-1], collapse = ","),
                                           no. = x)
    
    tomato_D_imp_rmse_r2_5s = rbind(tomato_D_imp_rmse_r2_5s, tomato_D_imp_rmse_r2_5s_i)
    
  }
}


##  it is repeated for three times, selecting the newly selected important variable
##  extract the best model
tomato_D_imp_r2_max = data.frame()
for (i in c(as.character(unique(tomato_D_imp_rmse_r2_5s$target)))){
  t_i = tomato_D_imp_rmse_r2_5s[tomato_D_imp_rmse_r2_5s$target ==i,]
  t_i_max = t_i[t_i$r2 == max(t_i$r2),]
  # if there are multiple predictors lengths that lead to the same r2, then choose the shortest length
  if( dim(t_i_max)[1] >1) t_i_max_x = t_i_max[t_i_max$no. == min(t_i_max$no.),][1,] else t_i_max_x = t_i_max
  
  tomato_D_imp_r2_max = rbind(tomato_D_imp_r2_max, t_i_max_x)
}

tomato_D_imp_r2_max# where important features, R2 and P values are stored


## can be plotted
####   one example of plotting  "SSC"

####   1, SSC
# tomato_D_imp_r2_max = tomato_D_var_r2_max
for (i in c(as.character(unique(tomato_D_imp_rmse_r2_5s$target)))[1]){
  print(i)
  pred_imp = unlist(strsplit(as.character(tomato_D_imp_r2_max[tomato_D_imp_r2_max$target == i,]$predictors), ","))
  tomato_D_spec_chemics_c_i = data.frame(tomato_D_spec_chemics_c[i], tomato_D_spec_chemics_c[pred_imp])
  
  RF_model_pred_v_i =c()
  for (j in c(1:dim(tomato_D_spec_chemics_c_i)[1])){
    set.seed(5433)
    if( length(pred_imp) ==1) {
      tomato_D_spec_chemics_c_i_j = data.frame(tomato_D_spec_chemics_c_i[-j, -1])
      names(tomato_D_spec_chemics_c_i_j) = names(tomato_D_spec_chemics_c_i)[2]
    } else {
      tomato_D_spec_chemics_c_i_j = tomato_D_spec_chemics_c_i[-j, -1]
    }
    
    RF_model_final_i_j = randomForest(tomato_D_spec_chemics_c_i_j, tomato_D_spec_chemics_c_i[-j,1])
    # RF_model_final_i
    
    ## make predictions
    RF_model_pred_i_j=predict(RF_model_final_i_j, tomato_D_spec_chemics_c_i[j,-1])
    
    RF_model_pred_v_i = append(RF_model_pred_v_i, unname(RF_model_pred_i_j))
  }
  
  ###########################
  df_i=data.frame(pred=RF_model_pred_v_i, actual=tomato_D_spec_chemics_c_i[,1],ind=c(1:length(tomato_D_spec_chemics_c_i[,1])))
  
  library(reshape2)
  df_i_m <- melt(df_i, id.vars="ind")
  # head(df__im)
  
  size = 18
  #
  ## sepcify the ylab 
  y_lab_1 = as.character(names(tomato_D_spec_chemics_c_i)[1])
  y_lab_1 = if (y_lab_1 =="TTA") paste(y_lab_1, " (%)", sep = "") else y_lab_1
  y_lab_1 = if (y_lab_1 =="SSC") paste(y_lab_1, " (°Brix)", sep = "") else y_lab_1
  
  imp_pt_i_a_1 = ggplot(df_i_m, aes(ind,value, col=variable)) +
    geom_point(aes(shape=variable))+
    geom_line(aes(linetype=variable))+
    scale_x_continuous(limits = c(0, 34), breaks = c(seq(1, 33, 4)))+
    scale_y_continuous(limits = c(3, 6), breaks = c(seq(3, 6, 0.5)))+
    scale_color_manual(labels = c(paste(i, " (Prediction)", sep = ""), paste(i, " (Ground truth)", sep = "")), values = c("blue", "red")) +
    scale_shape_manual(labels = c(paste(i, " (Prediction)", sep = ""), paste(i, " (Ground truth)", sep = "")), values = c(16, 17)) +
    scale_linetype_manual(labels = c(paste(i, " (Prediction)", sep = ""), paste(i, " (Ground truth)", sep = "")), values=c("dotted", "solid")) +
    theme(axis.text.x = element_text(color = "black", angle = 0, vjust = 0.1,size=size),
          axis.text.y = element_text(color = "black",size=size),
          plot.title = element_text(colour = "black"), legend.title=element_blank(),
          legend.text = element_text(color = "black",size=size),
          axis.title=element_text(size=size),
          legend.position = c(0.7, 0.2),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.key = element_rect(fill = NA, color = NA),
          strip.background = element_rect(fill = "white", colour = "black"),
          strip.text.x = element_text(size = 8, colour = "black", angle = 0, face="bold"),
          strip.text = element_text(size = size, colour = "black", angle = 0, face="bold")) +
    ylab(y_lab_1)+ xlab("Individual tomato")
  
  tiff(paste(paste("tomato", as.character(names(tomato_D_spec_chemics_c_i)[1]), "prediction_shit_a", "imp", "FT_final", "median", sep = "_"), ".tiff", sep = ""), 
       res=300, compression = "lzw", height=4, width=6, units="in")
  print(imp_pt_i_a_1)
  dev.off()
  # 
  
  ###
  # r2 = round(rsq(tomato_D_spec_chemics_c_i[,1],RF_model_pred_v_i), 2)
  # evaluate the quality of prediction model
  lm_i = lm(RF_model_pred_v_i ~ tomato_D_spec_chemics_c[,i])
  
  p_value_i = format(summary(lm_i)$coefficients[2,4], digits = 3) 
  r2 = round(summary(lm_i)$r.squared, digits = 2)
  
  ## TRICK TO PLOT LEGNEND FOR TWO LINES
  modele = glm(pred ~ actual, data=df_i)
  coefs = data.frame(intercept=coef(modele)[1],slope=coef(modele)[2])
  coefs= rbind(list(0,1), coefs)
  regression=as.factor(c('1:1 line', 'Fitting line'))
  
  coefs=cbind(coefs, regression)
  coefs$regression = factor(coefs$regression, levels =c('Fitting line', '1:1 line'))
  
  imp_pt_i_b_1 <- ggplot()+
    geom_point(data=df_i,aes(x=actual,y=pred),size=1)+
    geom_abline(data=coefs, aes(intercept=intercept,slope=slope,linetype=regression), show.legend=TRUE) +
    # ggtitle(paste("RandomForest Regression R^2=", r2, sep=""))+
    scale_x_continuous(limits=c(3.5, 5.5),
                       breaks=seq(3.5, 5.5, 0.5), name=paste(y_lab_1, " (Ground truth)", sep = ""))+
    scale_y_continuous(limits=c(3.5, 5.5),
                       breaks=seq(3.5, 5.5, 0.5), name=paste(y_lab_1, " (Prediction)", sep = ""))+
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5,colour="black",size=size),
          axis.text.y = element_text(colour="black",size=size),
          axis.title=element_text(size=size),
          panel.background = element_rect(fill = 'white', colour = 'black',size=0.1),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
          legend.title=element_blank(),
          legend.text =element_text(colour = "black",size=size), 
          legend.key = element_rect(colour = "transparent", fill = "white"),
          legend.position=c(0.83, 0.13))+
    annotate(geom="text", x = 3.5 +(5.5 - 3.5)/4, y = 5.5 - (5.5 - 3.5)/5, hjust = 0.9, vjust = 0,
             label= deparse(bquote(atop(R^2 == .(r2), "P" == .(p_value_i)))), parse=TRUE, size = 6,
             color="black")+
    annotate(geom="text", x = 3.5*1.01, y = 5.5*0.99, hjust = 0.9, vjust = 0,
             label= "(a)", size = 6,
             color="black")
  
  tiff(paste(paste("tomato", as.character(names(tomato_D_spec_chemics_c_i)[1]), "prediction_shit_b", "imp", "FT_final", "median", sep = "_"), ".tiff", sep = ""), 
       res=300, compression = "lzw", height=4, width=6, units="in")
  print(imp_pt_i_b_1)
  dev.off()
}










































  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  