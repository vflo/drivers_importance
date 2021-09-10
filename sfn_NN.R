sfn_NN <- function( data, predictors, nam_target, weights = NULL,nn = NULL,
                    do_predict = TRUE,trainfrac = 0.75, method="nnet", 
                    seed = 1, hidden = NULL ){
  library(magrittr)
  scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)
  data %>% 
    ungroup() %>% 
    dplyr::select(one_of(nam_target),
                  one_of(predictors))%>% 
    # mutate_if(is.numeric, scale2) %>% 
    mutate_if(is.character, as.factor) %>% 
    mutate_if(is.factor, as.numeric)-> data
  data <- data[complete.cases(data),] 
  if(!is.null(weights)){weights <- data[[weights]]}
  forml  <- as.formula(  paste( nam_target, "~", paste( predictors, collapse=" + " ) ) )
  
  if (method == "nnet"| method == 'rf'){
    
    ###############################################################
    ## NNET
    ###############################################################
    
    if (method == "nnet"){
      
      require( nnet )
      require( caret )
      
      if (is.null(nn)){
        
        preprocessParams <- caret::preProcess( data, method = c('center',"scale") )
        
        
        max_neurons <- nrow(data)*trainfrac/(2*(length(predictors)+1))
        # https://stats.stackexchange.com/questions/181/how-to-choose-the-
        # number-of-hidden-layers-and-nodes-in-a-feedforward-neural-netw
        
        if (is.null(hidden)){
          tune_grid <- expand.grid( .decay = c(0.5, 0.1, 1e-3, 1e-5, 1e-7), 
                                    .size = seq(1,round(max_neurons),1)
          )
        } else {
          tune_grid <- expand.grid( .decay = c(0.5, 0.1, 1e-3, 1e-5, 1e-7), .size = c(hidden) )
        }
        
        
        particiones  <- 10
        repeticiones <- 5
        
        set.seed(seed)
        seeds <- vector(mode = "list", length = (particiones * repeticiones) + 1)
        for (i in 1:(particiones * repeticiones)) {
          seeds[[i]] <- sample.int(1000, nrow(tune_grid))
        }
        seeds[[(particiones * repeticiones) + 1]] <- sample.int(1000, 1)
        
        traincotrlParams <- caret::trainControl( method = "repeatedcv",
                                                 number = particiones, 
                                                 repeats = repeticiones,
                                                 verboseIter = FALSE,
                                                 seeds = seeds,
                                                 p = trainfrac ) ## take best of 5 repetitions of 
        ## training franction

        nn <- caret::train(
          y           = data[[nam_target]],
          x           = data %>%
            dplyr::select(one_of(predictors)) %>% as.data.frame(),
          # form        = forml,
          # data        = data, #training,
          weights     = weights,
          method      = "nnet",
          linout      = TRUE,
          tuneGrid    = tune_grid,
          preProcess  = c('center',"scale"),
          trControl   = traincotrlParams,
          trace       = FALSE,
          metric      = "RMSE",
          maximize    = FALSE,
          importance=TRUE
        )
        
      }
      
    }   
    
    ###############################################################
    ## RF
    ###############################################################
    
    if (method == 'rf'){
      
      require( randomForest )
      require( caret )
      
      if (is.null(nn)){
        
        preprocessParams <- caret::preProcess( data, method = c('center',"scale") )
        
        traincotrlParams <- caret::trainControl( method = "repeatedcv",
                                                 number = 5, 
                                                 repeats = 5,
                                                 verboseIter = FALSE, 
                                                 p = trainfrac ) ## take best of 5 repetitions of 
        ## training fraction
        
        # if (is.null(hidden)){
        #   tune_grid <- expand.grid( mtry = seq(1,round(max_neurons),1)
        #   )
        # } else {
        #   tune_grid <- expand.grid( .decay = c(0.1), .size = c(hidden) )
        # }
        # 
        mtry <- sqrt(ncol(data))
        set.seed(seed)
        nn <- caret::train(
          form        = forml,
          data        = data, #training,
          weights     = weights,
          method      = "rf",
          linout      = TRUE,
          tuneGrid    = expand.grid(.mtry=mtry),
          preProcess  = c('center',"scale"),#'range',
          trControl   = traincotrlParams,
          trace       = FALSE,
          importance  = TRUE
        )
        
      }
      
    }
    
    ###############################################################
    ## Do predicts
    ###############################################################
    
    if (do_predict){
      vals <- as.vector( predict( nn, data ) )  # try( predict( nn, newdata=testing ) )
    } else {
      vals <- rep( NA, nrow(data) )
    }
    
    
    ###############################################################
    ## Error message
    ###############################################################
  } else {
    
    rlang::abort("predict_nn(): No other training methods implemented than nnet or rf.")
    
  }
  
  r2 <- caret::R2(pred = vals, obs = data[[nam_target]], na.rm = TRUE)
  
  return( list( nn=nn$results[rownames(nn$bestTune),], nn_model = nn, vals=vals, hidden_best=nn$bestTune,r2 =r2 ) )
  
}
