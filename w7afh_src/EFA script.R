### Notes ################################################################################

# This script reproduces the following analyses:
# 1) Hierarchical EFA + factor congruence coefficients
# 2) One-factor EFA based on factor correlations 
# 3) Correlations between four factors and other traits
# 4) Four-factor EFA and factor scoring (reported in Appendix)
# 5) Descriptive statistics and reliability (reported in Appendix)

# Important abbreviations:
#   USM = US Mturk sample
#   USS = US student sample
#   ITS = Italian student sample


### Load data and packages #############################################################

library(data.table)
library(psych)
library(haven)

# 128 items from NFSC, MAAS, IUS-12, URS
items=c(
  "IAL_absolutism_01", "IAL_absolutism_02", "IAL_absolutism_03", "IAL_absolutism_04", "IAL_absolutism_05", "IAL_absolutism_06", "IAL_absolutism_07", "IAL_absolutism_08", "IAL_absolutism_09", 
  "IAL_complexity_01", "IAL_complexity_02", "IAL_complexity_03", "IAL_complexity_04", "IAL_complexity_05", "IAL_complexity_06", "IAL_complexity_07", "IAL_complexity_08", "IAL_complexity_09", "IAL_complexity_10", 
  "IAL_discomfort_01", "IAL_discomfort_02", "IAL_discomfort_03", "IAL_discomfort_04", "IAL_discomfort_05", "IAL_discomfort_06", "IAL_discomfort_07", "IAL_discomfort_08", "IAL_discomfort_09", 
  "NFC_ambiguity_01", "NFC_ambiguity_02", "NFC_ambiguity_03", "NFC_ambiguity_04", "NFC_ambiguity_05", "NFC_ambiguity_06", "NFC_ambiguity_07", "NFC_ambiguity_08", "NFC_ambiguity_09", 
  "NFC_predict_01_R", "NFC_predict_03_R", "NFC_predict_04", "NFC_predict_05_R", "NFC_predict_06", "NFC_predict_07", "NFC_predict_08", "NFC_predict_09", 
  "NFC_order_01", "NFC_order_02_R", "NFC_order_03", "NFC_order_04", "NFC_order_05_R", "NFC_order_06", "NFC_order_07_R", "NFC_order_08", "NFC_order_09", 
  "NFC_order_10", "NFC_closedmind_01_R", "NFC_closedmind_02", "NFC_closedmind_03", "NFC_closedmind_04_R", "NFC_closedmind_05_R", "NFC_closedmind_06_R", "NFC_closedmind_07_R", "NFC_closedmind_08", 
  "NFC_decideQuickly_01", "NFC_decideQuickly_02", "NFC_decideQuickly_03", "NFC_decideQuickly_04", "NFC_decideQuickly_05", "NFC_decideQuickly_06", 
  "IU_inhibitory_01", "IU_inhibitory_02", "IU_inhibitory_03", "IU_inhibitory_04", "IU_inhibitory_05", 
  "IU_prospective_01", "IU_prospective_02", "IU_prospective_03", "IU_prospective_04", "IU_prospective_05", "IU_prospective_06", "IU_prospective_07", 
  "URS_change_01", "URS_change_02", "URS_change_03", "URS_change_04", "URS_change_05", "URS_change_06", "URS_change_07", "URS_change_08", "URS_change_09", "URS_change_10", "URS_change_11", "URS_change_12", "URS_change_13", "URS_change_14", "URS_change_15", "URS_change_16", 
  "URS_cognitive_01", "URS_cognitive_02", "URS_cognitive_03", "URS_cognitive_04", "URS_cognitive_05", "URS_cognitive_06", "URS_cognitive_07", "URS_cognitive_08", "URS_cognitive_09", "URS_cognitive_10", "URS_cognitive_11", "URS_cognitive_12", "URS_cognitive_13", "URS_cognitive_14", "URS_cognitive_15", "URS_cognitive_16", "URS_cognitive_17", 
  "URS_emotional_01", "URS_emotional_02", "URS_emotional_03", "URS_emotional_04", "URS_emotional_05", "URS_emotional_06", "URS_emotional_07", "URS_emotional_08", "URS_emotional_09", "URS_emotional_10", "URS_emotional_11", "URS_emotional_12", "URS_emotional_13", "URS_emotional_14"
)

# reading data
dt=data.table(read_sav(file.choose()))
usm = dt[sample=="mturk",(items),with=F] # US Mturk sample
uss = dt[sample=="student",(items),with=F] # US Student sample

dt2=data.table(read_sav(file.choose()))
its = dt2[,(items),with=F]


### Exploratory Factor Analysis #######################################################

library(psych)

# Parallel analysis to determine how many factors to extract

par(mfrow=c(1,1))
par.usm=fa.parallel(usm, fm='ml', fa='fa', show.legend=T, main = 'US Mturk') 
par.uss=fa.parallel(uss, fm='ml', fa='fa', show.legend=F, main = 'US student')
par.its=fa.parallel(its, fm='ml', fa='fa', show.legend=F, main = 'Italian student')

# BassAckward EFA

# We use bassAckward() to get the correlations between factors from successive solutions, so that we can create the hierarchical diagrams. The correlations are contained in the "bass.ack" result (See capture.output below). 
# In multi-factor solutions, the factors obtained from fa() and bassAckward() can be matched by their sequential order (i.e., first factor in fa() is the same as the first factor from bassAckward()).
# Note: BassAckwards is consistent with faCor but not factor.scores

bass.usm=bassAckward(usm, nfactors=7, fm='ml', cut = .45, lr=F, items=F, plot=T)
bass.uss=bassAckward(uss, nfactors=8, fm='ml', cut = .45, lr=F, items=F, plot=T)
bass.its=bassAckward(its, nfactors=9, fm='ml', cut = .45, lr=F, items=F, plot=T)

bassAckward.diagram(bass.its, lr=T, items=F, cut=.6)

capture.output(bass.usm[["bass.ack"]],file='bassAck USM.csv')
capture.output(bass.uss[["bass.ack"]],file='bassAck USS.csv')
capture.output(bass.its[["bass.ack"]],file='bassAck ITS.csv')


# EFA for factor interpretation and other 
#
# Here, we manually conduct the EFA for each level so that we can:
#   (1) interpret the factors in each solution
#   (2) compute the factor congruence coefficients
#   (3) compute the proportion of the variance accounted by each factor
#
# First, we will load the item labels to be appended to the fa results and create lists for storing the datasets and the results

labels=read.csv(file.choose(),row.names = 'names',stringsAsFactors = F) # load the item labels
dl=list(usm,uss,its) # list of datasets for EFA
for (i in c("res","load","pv")) assign(i, replicate(3, list())) # create new lists with 3 empty lists within
fcong=list()

# In the loop below, users define the datasets and the number of solutions to obtain, and the script returns the fa() results, the structure matrices with item labels, and factor congruence coefficients across samples.
# res stores complete fa results 
# load stores structure matrices
# pv stores percent of variance accounted by the factors in each solution
# fcong stores matrices with factor congruence coefficients

for (i in 1:9) { # i = number of solutions we want
  for(s in seq_along(dl)) { # for each sample in "dl"...
    
    id=c('M','U','I')[s] # letter identifier for each sample: M = Mturk, U = US students, I = Italian students
    
    fa1=fa(dl[[s]],nfactors=i,fm='ml') # conduct EFA (oblimin rotation, ML estimation)
    
    pv[[s]][[i]]=fa1$Vaccounted[2,] # proportion variance accounted for by each factor
    
    colnames(fa1$Structure)=sub('ML',id,colnames(fa1$Structure)) # rename the factor labels according to sample
    
    res[[s]][[i]]=fa1 # append fa results to list
    
    load[[s]][[i]]=merge(unclass(fa.sort(fa1$Structure)), labels, by='row.names', sort=F) # attach item labels to structure loadings and append to list
    
    write.table(rbind(load[[s]][[i]],''),
                file=paste0("load",s,".csv"), row.names=F, sep=',', append=T) # write matrix of structure loadings to file
  }
  
  fcong[[paste0(i,'_USM_USS')]]=fa.congruence(res[[1]][[i]],res[[2]][[i]],structure=T)  # append factor congruence coefficients for usm and uss
  fcong[[paste0(i,'_USM_ITS')]]=fa.congruence(res[[1]][[i]],res[[3]][[i]],structure=T)
  fcong[[paste0(i,'_USS_ITS')]]=fa.congruence(res[[2]][[i]],res[[3]][[i]],structure=T)
}

capture.output(fcong,file='fcong.csv')
capture.output(pv,file='pvaccounted.csv')

rm(fa1,fcong,load,pv,res,i,s,id) # cleaning afterwards



### One-factor EFA ##########################################################################
#
# Here, we obtain the correlations among the four factors (obliquely rotated) and use them to find the factor loadings on a general factor
# We also compute and save the factor scores for the 4 factors and the general factor

# 4-factor EFA
fa4.usm=fa(usm[,(items),with=F], nfactors=4, fm='ml') 
fa4.uss=fa(uss[,(items),with=F], nfactors=4, fm='ml')
fa4.its=fa(its[,(items),with=F], nfactors=4, fm='ml')

# 1-factor EFA using factor correlation matrics
fa(fa4.usm[["Phi"]], fm='ml')
fa(fa4.uss[["Phi"]], fm='ml')
fa(fa4.its[["Phi"]], fm='ml')

# saving factor scores
fnames=c("Fscores_UD","Fscores_PNC","Fscores_PS","Fscores_CT")
usm[,(fnames):=data.frame(factor.scores(.SD,fa4.usm)[["scores"]]),.SDcols=items]
uss[,(fnames):=data.frame(factor.scores(.SD,fa4.uss)[["scores"]]),.SDcols=items]
its[,(fnames):=data.frame(factor.scores(.SD,fa4.its)[["scores"]]),.SDcols=items]

# 1-factor EFA using the factor scores (we need these scores in order to obtain the general factor scores)
fa1.usm=fa(usm[,(fnames),with=F], nfactors=1, fm='ml')
fa1.uss=fa(uss[,(fnames),with=F], nfactors=1, fm='ml')
fa1.its=fa(its[,(fnames),with=F], nfactors=1, fm='ml')

# saving scores for the general factor
usm[,Fscores_GF:=factor.scores(.SD,fa1.usm)[["scores"]],.SDcols=fnames]
uss[,Fscores_GF:=factor.scores(.SD,fa1.uss)[["scores"]],.SDcols=fnames]
its[,Fscores_GF:=factor.scores(.SD,fa1.its)[["scores"]],.SDcols=fnames]

# correlations between the general factor and NFC, IA, IU

vars1=c("Fscores_GF", "IA", "NFC", "IU", "URS")
lowerCor(usm[,(vars1),with=F])
lowerCor(uss[,(vars1),with=F])
lowerCor(its[,(vars1),with=F])




### Correlations between the four factors and other traits #########################
# Note: This analysis was done using the combined US datasets (Mturk + students).

fa4=fa(dt[,(items),with=F], nfactors=4, fm='ml')
fnames=c("Fscores_UD","Fscores_PNC","Fscores_PS","Fscores_CT")
dt[,(fnames):=data.frame(factor.scores(.SD,fa4)[["scores"]]),.SDcols=items]

vars2=c(
  # these are the saved factor scores
  "Fscores_UD","Fscores_PNC","Fscores_PS","Fscores_CT", 
  # the Big Five and their aspects
  "BFA_Ac", "BFA_Ap", "BFA_Ci", "BFA_Co", "BFA_Ea", "BFA_Ee", "BFA_Nv", "BFA_Nw", "BFA_Oi", "BFA_Oo", 
  "BFA_A", "BFA_C", "BFA_E", "BFA_N", "BFA_O",
  # Sensation Seeking, Curiosity/Exploration, Need for Cognition, Naive Epistemic Beliefs, Preventive Coping, Anxious Coping
  "BSSS", "CEI", "NFCognition", "EBI", "TO_control", "TO_hypervig", 
  # Decision-making variables
  "DM_avoidant", "DM_certainty", "DM_hypervig", "DM_intuitive", "DM_procrast", "DM_spontaneous", "DM_rational")

lowerCor(dt[,(vars2),with=F])



### Measuring the Four factors ##########################################################################
#
# In the Appendix, we provide the structure matrix for the four-factor solution based on the combined US and Italian samples. 
# We also suggest items for scoring the four factors


# merging US and Italian data for the 128 items
us = dt[,(items)]
usits=rbindlist(list(us,its),fill=T)

# Four-factor EFA using the merged data
labels=read.csv(file.choose(),row.names = 'names',stringsAsFactors = F)
fa1=fa(usits,nfactors=4,fm='ml')
load=merge(unclass(fa.sort(fa1$Structure)), labels, by='row.names', sort=F)
write.csv(load,"USIT 4-factor EFA.csv")

# Items for scoring the four factors
# Note: these scores were not used in any of the analyses in the paper

toscore= list(
  meanscore_UD=c("URS_emotional_08", "IU_inhibitory_03", "IAL_discomfort_05", "URS_emotional_09", "IU_prospective_01", "URS_emotional_01", 
                "NFC_ambiguity_03", "IU_inhibitory_02", "URS_emotional_02", "IAL_discomfort_04", "IU_inhibitory_04"),
  meanscore_PS=c("NFC_order_10", "URS_cognitive_01", "URS_cognitive_17", "URS_cognitive_16", "IU_prospective_07", "NFC_order_09", "URS_cognitive_15", 
                 "URS_cognitive_09", "NFC_order_08", "URS_cognitive_10", "URS_cognitive_05", "NFC_order_01"),
  meanscore_CT=c("IAL_absolutism_03", "IAL_absolutism_01", "IAL_absolutism_05", "IAL_absolutism_09", "IAL_absolutism_04", "IAL_absolutism_06", 
                   "IAL_absolutism_02", "IAL_absolutism_08", "IAL_absolutism_07", "NFC_closedmind_02"), 
  meanscore_PNC=c("URS_change_01", "URS_change_16", "URS_change_03", "NFC_predict_03_R", "URS_change_15", "IAL_complexity_06", "URS_change_11", "IAL_complexity_02", 
                  "IAL_complexity_04", "URS_change_04", "IAL_complexity_10", "URS_change_09", "IAL_complexity_01")
#  meanscore_PNC_all=c("URS_change_03", "URS_change_16", "URS_change_01", "URS_change_11", "URS_change_15", "NFC_predict_03_R", "URS_change_09", "URS_change_12", "IAL_complexity_02", "URS_change_08", "URS_change_06", 
#               "IAL_complexity_01", "IAL_complexity_04", "IAL_complexity_10", "IAL_complexity_06", "IAL_complexity_05", "URS_change_04", "IAL_complexity_07", "NFC_closedmind_06_R", "IAL_complexity_08", "IAL_complexity_09", "NFC_closedmind_07_R"),
#  meanscore_Novelty=c("URS_change_03", "URS_change_16", "URS_change_01", "URS_change_11", "URS_change_15", "NFC_predict_03_R", "URS_change_09", 
#               "URS_change_12", "IAL_complexity_02", "URS_change_08", "URS_change_06"),
#  meanscore_Complexity=c("IAL_complexity_01", "IAL_complexity_04", "IAL_complexity_10", "IAL_complexity_06", "IAL_complexity_05", "URS_change_04", 
#                  "IAL_complexity_07", "NFC_closedmind_06_R", "IAL_complexity_08", "IAL_complexity_09", "NFC_closedmind_07_R")
)

loadings=sapply(toscore, function(x) fa.sort(fa(dt[, (x), with=F])$loadings[])) # item loadings on each factor
dt[,(names(toscore)):=lapply(toscore, function(x) rowMeans(.SD[,x,with=F],na.rm=T))] # compute the variables
dt[, meanscore_NovelComplex:=rowMeans(.SD,na.rm=T), .SDcols=c("meanscore_Complexity","meanscore_Novelty")] # average of Novelty + Complexity subfactors


### Descriptive Statistics and Reliability ######################################################################

# in order to use psych::describe make sure all the data is of data.frame class not data.table
usm = data.frame(dt[sample=="mturk"]) # US Mturk sample
uss = data.frame(dt[sample=="student"]) # US Student sample
its = data.frame(dt2)

# note: these are all the scored traits in the US data. For the Italian data, include only the IU, IA, NFC, URS variables.
vars=c("IU", "IU_inhibitory", "IU_prospective", "IAL", "IAL_absolutism", "IAL_complexity", "IAL_discomfort", "URS", "URS_change", "URS_cognitive", 
       "URS_emotional", "NFC", "NFC_order", "NFC_predict", "NFC_ambiguity", "NFC_closedmind", "NFC_decideQuickly"
       ,"BFA_Ac", "BFA_Ap", "BFA_Ci", "BFA_Co", "BFA_Ea", "BFA_Ee", "BFA_Nv", "BFA_Nw", "BFA_Oi", "BFA_Oo", "BFA_A", "BFA_C", "BFA_E", "BFA_N", "BFA_O", 
       "BSSS", "CEI", "NFCognition", "TO_control", "TO_hypervig", "DM_avoidant", "DM_certainty", "DM_hypervig", "DM_intuitive", "DM_procrast", "DM_spontaneous", "DM_rational"
      )
dstats=describe(usm[,(vars)])[,c("n", "mean", "sd", "min", "max", "range", "skew", "kurtosis")] # statistics we want

vars2=names(usm)[!(names(usm)) %in% vars] # select only the observed (unscored) variables in order to find alphas
alphas=sapply(vars, function(x) alpha(usm[, grep(x, vars2)],check.keys = T)$total$std.alpha) # Chronbach's alpha
dstatsAlpha=merge(descriptives,alphas,by="row.names") # merge dstats and alphas

# double-check if items are scored in the correct direction
# traits=unique(str_match(names(usm),".*(?=_[0-1])"))
# loadings=sapply(traits, function(x) fa.sort(fa(its[, grep(x, names(its)), with=F])$loadings[]))