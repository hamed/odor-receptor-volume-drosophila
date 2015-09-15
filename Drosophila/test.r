#This simple code provided to test volume sensitivity of Olfactory Receptors.
#Ref article:Olfactory receptors are sensitive to molecular volume of odorants
#Authors: Majid Saberi, Hamed Seyed-allaei
#doi: http://dx.doi.org/10.1101/013516
#January 6, 2015



rm ( list=ls())

setwd("~/")                                #Path of files location setting

#Extract molecular volumes and responses from supplementary files
DoOR <- as.matrix( read.csv("proposed_odorants.csv") )
rownames( DoOR ) <- DoOR[, 1]
DoOR <- DoOR[, 2:33]
mode( DoOR ) <- "numeric"

odorants <- as.matrix( read.csv("odorants.csv") )
volume <- as.numeric( odorants[, 4] )
names( volume ) <- odorants[, 1]
###



wt.mean <- function (x, wt)                   #Weighted mean function
  {
    s = which(!is.na(x * wt))
    wt = wt[s]
    x = x[s]
    return(sum(wt * x)/sum(wt))
  }
wt.sd <- function (x, wt)                     #Weighted sd function
  {
    s = which(!is.na(x + wt))
    wt = wt[s]
    x = x[s]
    xbar = wt.mean(x, wt)
    return(sqrt( sum(wt * (x^2))/sum(wt) - xbar^2)  )
  }

gauss <- function (x, mu=0, sig=1)            #Gaussian function
  {
    return( exp( -(x - mu)^2 / (2 * sig^2) ) )
  }


plotting <- function(OR)                            
  {
    response = DoOR[, OR]                           #OR's response selection
    sig_g = sd(volume, na.rm=T)                      #SD of g(v)
    mu_g = mean(volume, na.rm=T)                     #MEAN of g(v)
    sig_h = wt.sd(volume, response)                  #SD of h(v), 10th formula
    mu_h =  wt.mean(volume, response)                #MEAN of h(v), 9th formual

    sig_f = sqrt( sig_g^2 * sig_h^2 / (sig_g^2 - sig_h^2) )    #SD of f(v), 11th formula
    mu_f = sig_f^2 * (mu_h/sig_h^2 - mu_g/sig_g^2)             #MEAN of f(v), 12th formula

    f = gauss(volume, mu=mu_f, sig=sig_f)            #Making f(v) as a gaussian function
    plot(volume, response, ylim=c(0,1), xlim=c(0,250))
    lines(volume,  f )
  }


#Everyone can see volume sensitivity of ORs: 
#Scatterplot shows response versus molecular volume
#Line represents f(v), volume dependancy of OR or volume profile of OR
#for example: 
plotting("Or47a")


#List of jak-knife significant simulated binding packet profiles of olfactory receptors 
colnames(DoOR)
