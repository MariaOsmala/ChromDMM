


gradientDescent<- function(LambdaK,
                data,
                Z,
                hkm,
                hkm_lb,
                gradient,
                lambda_iter,
                lb,
                eta,
                nu,
                etah,
                nuh,
                verbose,
                MAX_GRAD_ITER=1000,
                reltol=1e-3, learning.rate=1e-3){

  

  
  #options(digits=20)
    iter <- 0
    i <- 0
    

    
    hkm_index <- vector(mode='integer', 1)
    hkm_lb_index <- vector(mode='integer', 1)
    gradient_index <- vector(mode='integer', 1)
    lb_index <- vector(mode='integer', 1)
    lambda_index <- vector(mode='integer', 1)
    params <- list(pi = Z, data = data, eta=eta, nu=nu, etah=etah, nuh=nuh,
                   hkm=hkm, hkm_lb=hkm_lb, gradient=gradient, lb=lb,hkm_index=hkm_index,
                   hkm_lb_index=hkm_lb_index,gradient_index=gradient_index,
                   lb_index= lb_index, lambda_index=lambda_index)
    #hkm_index <- vector(mode='integer', 1)
    #params <- list(pi = Z, data = data, eta=eta, nu=nu, etah=etah, nuh=nuh,
    #               hkm=hkm, hkm_index=hkm_index)
    
    cost <- neg_log_evidence_lambda_pi(LambdaK, lambda_iter, params)
    delta <- 1
    gradient.iter <- neg_log_derive_evidence_lambda_pi(LambdaK, lambda_iter, params)
    gradient.iter<- sqrt(sum(gradient.iter*gradient.iter))
    while(delta > reltol && i < MAX_GRAD_ITER){
      #print(i)
      i <- i + 1
      gr=neg_log_derive_evidence_lambda_pi(LambdaK, lambda_iter, params)
      LambdaK <- LambdaK - learning.rate * gr
      cval <- neg_log_evidence_lambda_pi(LambdaK, lambda_iter, params)
      cost <- append(cost, cval)
      gradient.iter=append(gradient.iter, sqrt( sum(gr*gr)) )
      delta <- abs(cost[i+1] - cost[i])
      if((cost[i+1] - cost[i]) > 0){
        print("The cost is increasing.  Try reducing learning rate.")
        return()
      }
      iter <- append(iter, i)
    }
    
    
    print(sprintf("Completed in %i iterations.", i))
  #   plot(cost)
  #   plot(gradient.iter)
  #   
  #   alpha_list =exp( do.call(rbind,lambda_iter))
  #   subset= seq(0,nrow(alpha_list), 500)
  #   subset[1]=1
  #   subset=c(subset, nrow(alpha_list))
  #   alpha.data=alpha_list[subset,]
  #   attr(alpha.data, "dimnames")=NULL
  #   alpha.data=alpha.data/rowSums(alpha.data)
  #   
  #   alpha_mk.data=reshape2::melt(alpha.data, varnames=c("iter", "bin"), value.name="prob")
  #   alpha_mk.data$iter=as.factor(alpha_mk.data$iter)
  #   
  #   
  #   axt=c(1,10,20,30,40,50)
  # ayt=subset
  #   ggplot(data=alpha_mk.data ,aes(x = bin, y = iter, height = prob)) +
  #     geom_density_ridges( stat = "identity", scale = 10)+theme(panel.background=element_blank(), panel.border = element_blank())
  #                                                               
  #   
    result <- list()
    result$par=LambdaK
    result$iter.nro=i
    result$value=cval
    result$gradient.norm=sqrt( sum(gr*gr))
    return(result)
   }