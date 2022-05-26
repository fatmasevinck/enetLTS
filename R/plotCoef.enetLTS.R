
plotCoef.enetLTS <- function(object,vers=c("reweighted","raw"),
                             colors=NULL,...){

  nam <- NULL

  if(is.null(colors)){
    colors <- list(bars="#0000AA",errorbars="red",
                   background="#BBBBEE",abline="#21A0D2",
                   scores="#0000AA",cutoffs="#00EEEE",
                   badouts="darkred", modouts="black")
  }
  family <- object$inputs$family
  vers <- match.arg(vers)

  coefficients <- object$coefficients
  raw.coefficients <- object$raw.coefficients


    if (family=="multinomial") {
      if (vers=="reweighted") {
        main <- "enetLTS coefficients for multinomial logistic regression"
        plotcoefs <- data.frame(Wavelength = rep(1:nrow(coefficients),ncol(coefficients)),
                                class = as.character(rep(1:ncol(coefficients),each=nrow(coefficients))),
                                coefficients = as.vector(coefficients))
        plot <- ggplot(data = plotcoefs, aes(x = Wavelength,
                                             y = coefficients,
                                             colour=class,
                                             linetype=class)) +
          geom_line() + geom_hline(yintercept=0, linetype="dashed",color="gray") +
          labs(title=paste(names(object$inputs$yy),main)) +
          theme(plot.title=element_text(size=rel(1.3)),
                axis.text.x=element_text(angle=90),
                axis.text.y=element_text(size = 13),
                axis.title.x=element_blank(),
                axis.title.y=element_blank()) +
          scale_x_continuous("", labels = as.character(1:nrow(coefficients)),
                             breaks = 1:nrow(coefficients))
        print(plot)
      } else if (vers=="raw"){
        main <- "enetLTS raw coefficients for multinomial logistic regression"
        raw.plotcoefs <- data.frame(Wavelength = rep(1:nrow(raw.coefficients),ncol(raw.coefficients)),
                                    class = as.character(rep(1:ncol(raw.coefficients),each=nrow(raw.coefficients))),
                                    raw.coefficients = as.vector(raw.coefficients))
        raw.plot <- ggplot(data = raw.plotcoefs, aes(x = Wavelength,
                                                     y = raw.coefficients,
                                                     colour = class,
                                                     linetype = class)) +
          geom_line() + geom_hline(yintercept=0, linetype="dashed",color="gray") +
          labs(title=paste(names(object$inputs$yy),main)) +
          theme(plot.title=element_text(size=rel(1.3)),
                axis.text.x=element_text(angle=90),
                axis.text.y=element_text(size = 13),
                axis.title.x=element_blank(),
                axis.title.y=element_blank()) +
          scale_x_continuous("", labels = as.character(1:nrow(raw.coefficients)), breaks = 1:nrow(raw.coefficients))
        print(raw.plot)
      }
    } else if (family=="binomial" | family=="gaussian") {
      # if coefficents is a vector
      if (vers=="reweighted") {
        if (family=="binomial") main <- "enetLTS coefficients for binary logistic regression"
        if (family=="gaussian") main <- "enetLTS coefficients for regression"
        plotcoefs <- data.frame(Wavelength = (1:length(coefficients)),
                                class = as.character(rep(1,length(coefficients))),
                                coefficients = as.vector(coefficients))
        plot <- ggplot(data = plotcoefs, aes(x = Wavelength,
                                             y = coefficients,
                                             colour = class,
                                             linetype = class)) +
          geom_line() + geom_hline(yintercept=0, linetype="dashed",color="gray") +
          labs(title=paste(names(object$inputs$yy),main)) +
          theme(plot.title=element_text(size=rel(1.3)),
                legend.position="none",
                axis.text.x=element_text(angle=90),
                axis.text.y=element_text(size = 13),
                axis.title.x=element_blank(),
                axis.title.y=element_blank())  +
          scale_x_continuous("", labels = as.character(1:length(coefficients)), breaks = 1:length(coefficients))
        print(plot)
      } else if (vers=="raw"){
        if (family=="binomial") main <- "enetLTS raw coefficients for binary logistic regression"
        if (family=="gaussian") main <- "enetLTS raw coefficients for regression"
        raw.plotcoefs <- data.frame(Wavelength=(1:length(raw.coefficients)),
                                    class = as.character(rep(1,length(raw.coefficients))),
                                    raw.coefficients = as.vector(raw.coefficients))
        plot <- ggplot(data = raw.plotcoefs, aes(x = Wavelength,
                                                 y = raw.coefficients,
                                                 colour = class,
                                                 linetype = class)) +
          geom_line() + geom_hline(yintercept=0, linetype="dashed",color="gray") +
          labs(title=paste(names(object$inputs$yy),main)) +
          theme(plot.title=element_text(size=rel(1)),
                legend.position="none",
                axis.text.x=element_text(angle=90),
                axis.text.y=element_text(size = 13),
                axis.title.x=element_blank(),
                axis.title.y=element_blank())  +
          scale_x_continuous("", labels = as.character(1:length(raw.coefficients)), breaks = 1:length(raw.coefficients))
        print(plot)
      }
    }


}
