# Functions for calculating I2 and R2 marginal and conditional from metafor models. 
# Functions written by Malgorzata (Losia) Lagisz, Shinichi Nakagawa and Wolfgang Viechtbauer 


### Functions for I2
I2 <- function(model, method = c("Wolfgang", "Shinichi")){
  warning("Make sure you have the observation (effec size) level random effect\n")
  ## evaluate choices
  method <- match.arg(method)
  
  # Wolfgang's method
  if(method == "Wolfgang"){
    W <- solve(model$V) 
    X <- model.matrix(model)
    P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
    I2_total <- sum(model$sigma2) / (sum(model$sigma2) + (model$k - model$p) / sum(diag(P)))
    I2_each  <- model$sigma2 / (sum(model$sigma2) + (model$k - model$p) / sum(diag(P)))
    names(I2_each) = paste0("I2_", model$s.names)
    
    # putting all together
    I2s <- c(I2_total = I2_total, I2_each)
    
    # or my way
  } else {
    # sigma2_v = typical sampling error variance
    sigma2_v <- sum(1/model$vi) * (model$k-1) / (sum(1/model$vi)^2 - sum((1/model$vi)^2)) 
    I2_total <- sum(model$sigma2) / (sum(model$sigma2) + sigma2_v) #s^2_t = total variance
    I2_each  <- model$sigma2 / (sum(model$sigma2) + sigma2_v)
    names(I2_each) = paste0("I2_", model$s.names)
    
    # putting all together
    I2s <- c(I2_total = I2_total, I2_each)
  }
  return(I2s)
}

### Functions for R2

R2 <- function(model){
  warning("Make sure you have the observation (effec size) level random effect as the last in the formula\n")
  
  # fixed effect variance
  fix <- var(as.numeric(as.vector(model$b) %*% t(as.matrix(model$X))))
  
  # marginal
  R2m <- fix / (fix + sum(model$sigma2))
  R2
  #Rm <- round(100*R2m, 3)
  
  # conditional
  R2c <- (fix + sum(model$sigma2) - model$sigma2[length(model$sigma2)]) / 
    (fix + sum(model$sigma2))
  
  R2s <- c(R2_marginal = R2m, R2_coditional = R2c)
  
  return(R2s)
}

### Functions for plotting obtained from: 
# https://github.com/JustinCally/SexualSelection
# which code originally came from: https://github.com/jaredlander/coefplot/blob/master/R/position.r

# Powers dodging, stacking and filling.
collidev <- function(data, height = NULL, name, strategy, check.height = TRUE) {
  # Determine height
  if (!is.null(height)) {
    # height set manually
    if (!(all(c("ymin", "ymax") %in% names(data)))) {
      data$ymin <- data$y - height / 2
      data$ymax <- data$y + height / 2
    }
  } else {
    if (!(all(c("ymin", "ymax") %in% names(data)))) {
      data$ymin <- data$y
      data$ymax <- data$y
    }
    
    # height determined from data, must be floating point constant
    heights <- unique(data$ymax - data$ymin)
    heights <- heights[!is.na(heights)]
    height <- heights[1]
  }
  
  # Reorder by x position, relying on stable sort to preserve existing
  # ordering, which may be by group or order.
  data <- data[order(data$ymin), ]
  
  # Check for overlap
  intervals <- as.numeric(t(unique(data[c("ymin", "ymax")])))
  intervals <- intervals[!is.na(intervals)]
  
  if (length(unique(intervals)) > 1 & any(diff(scale(intervals)) < -1e-6)) {
    warning(name, " requires non-overlapping y intervals", call. = FALSE)
    # This is where the algorithm from [L. Wilkinson. Dot plots.
    # The American Statistician, 1999.] should be used
  }
  
  if (!is.null(data$xmax)) {
    plyr::ddply(data, "ymin", strategy, height = height)
  } else if (!is.null(data$x)) {
    data$xmax <- data$x
    data <- plyr::ddply(data, "ymin", strategy, height = height)
    data$x <- data$xmax
    data
  } else {
    stop("Neither x nor xmax defined")
  }
}

# Dodge overlapping interval.
# Assumes that each set has the same horizontal position.
pos_dodgev <- function(df, height) {
  n <- length(unique(df$group))
  if (n == 1) return(df)
  
  if (!all(c("ymin", "ymax") %in% names(df))) {
    df$ymin <- df$y
    df$ymax <- df$y
  }
  
  d_height <- max(df$ymax - df$ymin)
  
  # df <- data.frame(n = c(2:5, 10, 26), div = c(4, 3, 2.666666,  2.5, 2.2, 2.1))
  # ggplot(df, aes(n, div)) + geom_point()
  
  # Have a new group index from 1 to number of groups.
  # This might be needed if the group numbers in this set don't include all of 1:n
  groupidy <- match(df$group, sort(unique(df$group)))
  
  # Find the center for each group, then use that to calculate xmin and xmax
  df$y <- df$y + height * ((groupidy - 0.5) / n - .5)
  df$ymin <- df$y - d_height / n / 2
  df$ymax <- df$y + d_height / n / 2
  df
}

position_dodgev <- function(height = NULL) {
  ggproto(NULL, PositionDodgeV, height = height)
}

PositionDodgeV <- ggproto(
  "PositionDodgeV", Position,
  required_aes = "y",
  height = NULL,
  setup_params = function(self, data) {
    if (is.null(data$ymin) && is.null(data$ymax) && is.null(self$height)) {
      warning("height not defined. Set with `position_dodgev(height = ?)`",
              call. = FALSE)
    }
    list(height = self$height)
  },
  
  compute_panel = function(data, params, scales) {
    collidev(data, params$height, "position_dodgev", pos_dodgev, check.height = FALSE)
  }
)


# extracting mean and CI from each metafor model
estimates.CI <- function(model){
  db.mf <- data.frame(model$b,row.names = 1:nrow(model$b))
  db.mf <- cbind(db.mf,model$ci.lb,model$ci.ub,row.names(model$b))
  names(db.mf) <- c("mean","lower","upper","estimate")
  return(db.mf[,c("estimate","mean","lower","upper")])
}