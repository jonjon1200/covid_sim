sim <- function(classSize, noYears, subjectTeachers, maxTime, backgroundROI, gammaParams, testType, probs, quarantineTime, catchRate, testDays, startingInfected, quarantineThreshold, teacherInfective, teacherInfectible, teacherTeacher, nonClass, weekend, maxIPeriod){
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
    dgamma(as.double(a["timeInState"][[1]]),shape=gammaParams["spread.shape"],scale=gammaParams["spread.scale"])*as.double(a["spreadMult"]) )
    outGamma <-apply(outClass, 1, function (a)
    dgamma(as.double(a["timeInState"][[1]]),shape=gammaParams["spread.shape"],scale=gammaParams["spread.scale"])*as.double(a["spreadMult"]) )*nonClass #Discount non class
    
    scaledGamma <- (sum(inGamma)+sum(outGamma))/dgamma(gammaParams["spread.shape"]*gammaParams["spread.scale"],shape=gammaParams["spread.shape"],scale=gammaParams["spread.scale"])
    return(scaledGamma)
  }
  
  newState <- function(time, student){
    if(student["state"]=="S"){
      infected <- data.frame(as.matrix(subset(timeList[[time]],(state=="I")))) #Auto doesn't include student as already known to be in S
      studentNo <- gammaRate(infected, student["yearID"], student["studentOrTeacher"])
      if ((time%%7) %in% c(6,7)){
        studentNo <- studentNo*weekend[1]
        backgroundROI <- backgroundROI*weekend[2]
      }#Basic weekend discount
  
      lambda <- backgroundROI + as.double(student["catchMult"][[1]]) * catchRate *  studentNo/popSize
      coinFlipI <- rbinom(1,1,min(1,lambda))
      coinFlipQ <- rbinom(1,1,(1-probs[[testType]][["TN"]]) * as.numeric((time%%7) %in% testDays))
      if (coinFlipQ==1){
        student["state"]="Q"
        student["timeInState"]=1
      }else if (coinFlipI>0){
        student["state"]="I"
        student["timeInState"]=1
      }else{
        student["timeInState"]<-as.double(student["timeInState"])+1
      }
    }
    
    else if(student["state"]=="I"){
      lambda <- dgamma(as.double(student["timeInState"][[1]]),shape=gammaParams["recover.shape"],scale=gammaParams["recover.scale"])
      scaledInfectiousness <-  dgamma(as.double(student["timeInState"][[1]]),shape=gammaParams["spread.shape"],scale=gammaParams["spread.scale"]) / dgamma(gammaParams["spread.shape"]*gammaParams["spread.scale"],shape=gammaParams["spread.shape"],scale=gammaParams["spread.scale"])
      coinFlipI <- rbinom(1,1,min(1,lambda))
      coinFlipQ <- rbinom(1,1,probs[[testType]][["TP"]]*scaledInfectiousness)
      if (coinFlipQ == 1 && (time%%7) %in% testDays){
        student["state"]="Q"
        student["timeInState"]=1
      }else if ( coinFlipI == 1 || as.double(student["timeInState"])==maxIPeriod){
        student["state"]="R"
        student["timeInState"]=1
      }else{
        student["timeInState"]<-as.double(student["timeInState"])+1
      }
    }
    
    else if(student["state"]=="Q"){
      if (student["timeInState"]==quarantineTime){
        student["state"]="R"
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
  popSize <- (classSize+1)*noYears + subjectTeachers
  studentID = c(1:popSize)
  yearID <- sort(rep(1:noYears,(classSize+1)))
  spreadMult <- rnorm(popSize, mean=1, sd=0.1)
  catchMult <- rnorm(popSize, mean=1, sd=0.1)
  state = rep("S", popSize)
  state[sample.int(popSize,startingInfected)]<-"I"
  timeInState = rep(1,popSize)
  studentOrTeacher <- rep(c(rep("s",classSize),"t"),noYears)
  
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
    Qno <- nrow(filter(timeList[[i]], state=="Q"))
    Rno <- popSize-Sno-Ino-Qno
    if (infectionEnd == -1 && Ino == 0){infectionEnd <- i}
    population <- add_row(population, "S"=Sno, "I"=Ino, "Q"=Qno, "R"=Rno)
  }
  
  ## ----Class Info--------------------------------------------------------------
  classVar <- c()
  classMax <- 0
  for (i in 1:noYears){
    tmp <- sapply(timeList,function(a)(nrow(a[which(a["yearID"]==i & a["state"]=="Q"),] )))
      classVar <- c(classVar, var(tmp))
      classMax <- max(classMax, max(tmp))
  }
  
  
  ## ----Summary-----------------------------------------------------------------
  summaryInfo <- data.frame(
  "totalQ"=sum(population["Q"]),
  "avgQ"=round(sum(population["Q"])/popSize,digits=1),
  "finalI"=population[maxTime,"R"],
  "infectionEnd"=infectionEnd,
  "overThreshold"=length(population["Q"][population["Q"]>quarantineThreshold]),
  "meanVar"=mean(classVar),
  "runTime"=as.double(Sys.time()-tb)
  ) 
  return(list("summary"=summaryInfo,"populations"=population))
}