sim <- function(classSize, noYears, maxTime, backgroundROI,
                gammaParams, testType, probs, pThreshold, quarantineTime, catchRate, testDays,
                startingInfected, quarantineThreshold, teacherInfective, teacherInfectible,
                teacherTeacher, classMult, nonClass, weekend, maxIPeriod, percVacc, vaccEff,
                whoMask, maskEff,asymRate,asymDiscount){
  ## ----setup, include=FALSE----------------------------------------------------
  tb<-Sys.time()
  knitr::opts_chunk$set(echo = TRUE)
  library(tibble)
  library(knitr)
  library(formatR)
  library(dplyr)
  opts_chunk$set(tidy.opts=list(width.cutoff=60),tidy=TRUE)
  
  
  ## ----functions---------------------------------------------------------------
  gammaRate <- function(infected, year, tOrS){

    if(nrow(infected)==0){return (0)}
    inClass <- infected[which(infected["yearID"]==year),]
    outClass <- infected[!infected["studentID"] %in% inClass["studentID"],]
    
    if (dim(inClass[which(inClass["studentOrTeacher"]=="t"),])[1]!=0){
      inClass[which(inClass["studentOrTeacher"]=="t"),]["spreadMult"] <- sapply(inClass[which(inClass["studentOrTeacher"]=="t"),]["spreadMult"], function(a) as.character(as.double(a)*teacherInfective)) #Teachers n times infective as typical student
    }
  
    if (tOrS == "t"){
      if (dim(inClass)[1]!=0){
        inClass["spreadMult"] <- sapply(inClass["spreadMult"], function(a) as.character(as.double(a)*teacherInfectible)) #Teacher n times infectible
      }
      if (dim(outClass[which(outClass["studentOrTeacher"]=="t"),])[1]!=0){
        outClass[which(outClass["studentOrTeacher"]=="t"),]["spreadMult"] <- sapply(outClass[which(outClass["studentOrTeacher"]=="t"),]["spreadMult"], function(a) as.character(as.double(a)*teacherTeacher)) #Teachers n times likely to spread to each other
      }
    }
  
    
    inGamma <- apply(inClass, 1, function (a)
    dgamma(as.double(a["timeInState"][[1]]),shape=gammaParams["spread.shape"],scale=gammaParams["spread.scale"])*as.double(a["spreadMult"]) )*classMult
    outGamma <-apply(outClass, 1, function (a)
    dgamma(as.double(a["timeInState"][[1]]),shape=gammaParams["spread.shape"],scale=gammaParams["spread.scale"])*as.double(a["spreadMult"]) )*nonClass #Discount non class
    
    scaledGamma <- (sum(inGamma)+sum(outGamma))/meanGamma
    return(scaledGamma)
  }
  
  newState <- function(time, student){
    if(student["state"]=="S"){
      infected <- data.frame(as.matrix(subset(timeList[[time]],(state=="I")))) #Auto doesn't include student as already known to be in S
      studentNo <- gammaRate(infected, student["yearID"], student["studentOrTeacher"])
      if ((time%%7) %in% c(5,6)){
        studentNo <- studentNo*weekend[1]
        backgroundROI <- backgroundROI*weekend[2]
      }#Basic weekend discount
  
      lambda <- backgroundROI + as.double(student["catchMult"][[1]]) * catchRate *  studentNo/popSize
      coinFlipI <- rbinom(1,1,min(1,lambda))
      coinFlipQ <- rbinom(1,pThreshold,(1-probs[[testType]][["TN"]]) * as.numeric((time%%7) %in% testDays))
      if (coinFlipQ==pThreshold){
        student["state"]="Qs"
        student["timeInState"]=1
      }else if (coinFlipI==1){
        student["state"]="I"
        student["timeInState"]=1
      }else{
        student["timeInState"]<-as.double(student["timeInState"])+1
      }
    }
    
    else if(student["state"]=="I"){
      lambda <- diff(pgamma(c(as.double(student["timeInState"][[1]])-1,as.double(student["timeInState"][[1]])),shape=gammaParams["recover.shape"],scale=gammaParams["recover.scale"]))
      scaledInfectiousness <-  dgamma(as.double(student["timeInState"][[1]]),shape=gammaParams["spread.shape"],scale=gammaParams["spread.scale"]) / meanGamma
      coinFlipR <- rbinom(1,1,lambda)
      coinFlipQ <- rbinom(1,pThreshold,min(1,probs[[testType]][["TP"]]*scaledInfectiousness))
      if (coinFlipQ == pThreshold && (time%%7) %in% testDays){
        student["state"]="Qi"
        student["timeInState"]=1
      }else if ( coinFlipR == 1 || as.double(student["timeInState"])==maxIPeriod){
        student["state"]="R"
        student["timeInState"]=1
      }else{
        student["timeInState"]<-as.double(student["timeInState"])+1
      }
    }
    
    else if(student["state"] %in% c("Qs","Qi")){
      if (student["timeInState"]==quarantineTime & student["state"]=="Qi"){
        student["state"]="R"
        student["timeInState"]=1
      }else if (student["timeInState"]==quarantineTime & student["state"]=="Qs"){
        student["state"]="S"
        student["timeInState"]=1         
      }else{
        student["timeInState"]=as.double(student["timeInState"])+1
      }
    }
    
    else{
      student["timeInState"]<-as.double(student["timeInState"])+1
    }
    return (student)
  }
  
  
  ## ----Start variables---------------------------------------------------------
  meanGamma<-integrate(function(x)dgamma(x,shape=gammaParams["spread.shape"],scale=gammaParams["spread.scale"])^2,0,Inf)$value
  popSize <- (classSize+1)*noYears 
  studentID = c(1:popSize)
  yearID <- sort(rep(1:noYears,(classSize+1)))
  spreadMult <- rnorm(popSize, mean=1, sd=0.1)
  catchMult <- rnorm(popSize, mean=1, sd=0.1)
  state = rep("S", popSize)
  state[sample.int(popSize,startingInfected)]<-"I"
  timeInState = rep(1,popSize)
  studentOrTeacher <- rep(c(rep("s",classSize),"t"),noYears)
  vaccinated<-sample(popSize,floor(popSize*percVacc))
  spreadMult[vaccinated]<-spreadMult[vaccinated]*(1-vaccEff$spread)
  catchMult[vaccinated]<-catchMult[vaccinated]*(1-vaccEff$catch)
  spreadMult[which(studentOrTeacher %in% whoMask)] <- spreadMult[which(studentOrTeacher %in% whoMask)]*(1-maskEff)
  asym<-sample(popSize,round(popSize*asymRate))
  spreadMult[asym]<-spreadMult[asym]*asymDiscount
  school <- tibble(studentID, yearID, spreadMult, catchMult, state, timeInState, studentOrTeacher)
  multFrame <- tibble(studentID, spreadMult, catchMult)
  
  
  ## ----Run timesteps-----------------------------------------------------------
  timeList <- list(school)
  for (i in (2:maxTime)){
    tempYear <- as_tibble(t(apply(timeList[[i-1]],1,newState,time=i-1)))
    timeList <- append(timeList, list(tempYear))
  }

  
  ## ----Plotting----------------------------------------------------------------
  population <- data.frame("S"=popSize-startingInfected,"I"=startingInfected, "Q"=0, "R"=0)
  infectionEnd <- -1
  for (i in (2:maxTime)){
    Sno <- nrow(filter(timeList[[i]], state=="S"))
    Ino <- nrow(filter(timeList[[i]], state=="I"))
    Qno <- nrow(filter(timeList[[i]], state%in% c("Qs","Qi")))
    Rno <- popSize-Sno-Ino-Qno
    if (infectionEnd == -1 && Ino == 0){infectionEnd <- i}
    population <- add_row(population, "S"=Sno, "I"=Ino, "Q"=Qno, "R"=Rno)
  }
  
  ## ----Class Info--------------------------------------------------------------
  classVar <- c()
  classMax <- 0
  
  for (i in 1:noYears){
    tmp <- sapply(timeList,function(a)(length(which(a["yearID"]==i & a["state"]==c("Qs","Qi")) )))
    classVar <- c(classVar, var(tmp))
    classMax <- max(classMax, max(tmp))
  }

  
  ## ----Summary-----------------------------------------------------------------
  summaryInfo <- data.frame(
  "totalQ"=sum(population["Q"]),
  "avgQ"=round(sum(population["Q"])/popSize,digits=1),
  "finalI"=population[maxTime,"R"]+population[maxTime,"I"],
  "infectionEnd"=infectionEnd,
  "overThreshold"=length(population["Q"][population["Q"]>quarantineThreshold]),
  "meanVar"=mean(classVar),
  "meanInfperDay"=mean(population[2:maxTime,"R"]-population[1:(maxTime-1),"R"]),
  "percVaccInfected"=if(percVacc==0) NA else length(which(timeList[[100]][vaccinated,"state"]=="R" | timeList[[100]][vaccinated,"state"]=="I"))/length(vaccinated),
  "percUnVaccInfeced"=if(percVacc==1) NA else length(which(timeList[[100]][-vaccinated,"state"]=="R" | timeList[[100]][-vaccinated,"state"]=="I"))/(popSize-length(vaccinated)),
  "runTime"=as.double(Sys.time()-tb)
  ) 
  return(list("summary"=summaryInfo,"populations"=population))
}