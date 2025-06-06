### General notes ####################################################################################

# This script replicates the CFA analysis, which involves:
# 1) defining the bifactor model
# 2) generating 200 parcel allocations for each datset (sample)
# 3) running CFA for each allocation and pooling the results
# 4) finding the N (allocations) at which the estimates stabilize

### Load data and packages #############################################################################

library(haven) # for SPSS sav files
library(semTools) # for parcelAllocation

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
dt = data.frame(read_sav(file.choose()))
usm = dt[dt$sample=="mturk",c(items)] # US Mturk sample
uss = dt[dt$sample=="student",c(items)] # US Student sample

dt2 = data.frame(read_sav(file.choose())) # Italian sample
its = dt2[,c(items)]

### Bifactor model Notes ##################################################################################

#  The models below were designed to be tested with semTools::parcelAllocation, which requires a model with all the items and another model with the parcels. Both models needs to be defined prior to analysis.
#  IU_inhibitory has only 5 items, so one item was left out of the parcellation (labelled IU_inhibitory_05)
#  General factors are uncorrelated with specific factors as specified by the Fg ~~ 0*Fs syntax, but the specific factors are allowed to correlate among themselves

### Bifactor model with items ##############################################################################

mod.items='
IA_ =~ IAL_absolutism_01 + IAL_absolutism_02 + IAL_absolutism_03 + IAL_absolutism_04 + IAL_absolutism_05 + IAL_absolutism_06 + IAL_absolutism_07 + IAL_absolutism_08 + IAL_absolutism_09 + IAL_complexity_01 + IAL_complexity_02 + IAL_complexity_03 + IAL_complexity_04 + IAL_complexity_05 + IAL_complexity_06 + IAL_complexity_07 + IAL_complexity_08 + IAL_complexity_09 + IAL_complexity_10 + IAL_discomfort_01 + IAL_discomfort_02 + IAL_discomfort_03 + IAL_discomfort_04 + IAL_discomfort_05 + IAL_discomfort_06 + IAL_discomfort_07 + IAL_discomfort_08 + IAL_discomfort_09
IA_abs =~ IAL_absolutism_01 + IAL_absolutism_02 + IAL_absolutism_03 + IAL_absolutism_04 + IAL_absolutism_05 + IAL_absolutism_06 + IAL_absolutism_07 + IAL_absolutism_08 + IAL_absolutism_09
IA_com =~ IAL_complexity_03 + IAL_complexity_04 + IAL_complexity_05 + IAL_complexity_06 + IAL_complexity_07 + IAL_complexity_08 + IAL_complexity_09 + IAL_complexity_10
IA_dis =~ IAL_discomfort_01 + IAL_discomfort_02 + IAL_discomfort_03 + IAL_discomfort_04 + IAL_discomfort_05 + IAL_discomfort_06 + IAL_discomfort_07 + IAL_discomfort_08 + IAL_discomfort_09

NFC_ =~ NFC_ambiguity_01 + NFC_ambiguity_02 + NFC_ambiguity_03 + NFC_ambiguity_04 + NFC_ambiguity_05 + NFC_ambiguity_06 + NFC_ambiguity_07 + NFC_ambiguity_08 + NFC_ambiguity_09 + NFC_closedmind_01_R + NFC_closedmind_02 + NFC_closedmind_03 + NFC_closedmind_04_R + NFC_closedmind_05_R + NFC_closedmind_06_R + NFC_closedmind_07_R + NFC_closedmind_08 + NFC_order_01 + NFC_order_02_R + NFC_order_03 + NFC_order_04 + NFC_order_05_R + NFC_order_06 + NFC_order_07_R + NFC_order_08 + NFC_order_09 + NFC_order_10 + NFC_predict_01_R + NFC_predict_03_R + NFC_predict_04 + NFC_predict_05_R + NFC_predict_06 + NFC_predict_07 + NFC_predict_08 + NFC_predict_09 + NFC_decideQuickly_01 + NFC_decideQuickly_02 + NFC_decideQuickly_03 + NFC_decideQuickly_04 + NFC_decideQuickly_05 + NFC_decideQuickly_06
NFC_amb =~ NFC_ambiguity_01 + NFC_ambiguity_02 + NFC_ambiguity_03 + NFC_ambiguity_04 + NFC_ambiguity_05 + NFC_ambiguity_06 + NFC_ambiguity_07 + NFC_ambiguity_08 + NFC_ambiguity_09
NFC_pre =~ NFC_predict_01_R + NFC_predict_03_R + NFC_predict_04 + NFC_predict_05_R + NFC_predict_06 + NFC_predict_07 + NFC_predict_08 + NFC_predict_09
NFC_ord =~ NFC_order_01 + NFC_order_02_R + NFC_order_03 + NFC_order_04 + NFC_order_05_R + NFC_order_06 + NFC_order_07_R + NFC_order_08 + NFC_order_09 + NFC_order_10
NFC_clo =~ NFC_closedmind_01_R + NFC_closedmind_02 + NFC_closedmind_03 + NFC_closedmind_04_R + NFC_closedmind_05_R + NFC_closedmind_06_R + NFC_closedmind_07_R + NFC_closedmind_08
NFC_dec =~ NFC_decideQuickly_01 + NFC_decideQuickly_02 + NFC_decideQuickly_03 + NFC_decideQuickly_04 + NFC_decideQuickly_05 + NFC_decideQuickly_06

IU_ =~ IU_inhibitory_01 + IU_inhibitory_02 + IU_inhibitory_03 + IU_inhibitory_04 + IU_inhibitory_05 + IU_prospective_01 + IU_prospective_02 + IU_prospective_03 + IU_prospective_04 + IU_prospective_05 + IU_prospective_06 + IU_prospective_07
IU_inh =~ IU_inhibitory_01 + IU_inhibitory_02 + IU_inhibitory_03 + IU_inhibitory_04 + IU_inhibitory_05
IU_pro =~ IU_prospective_01 + IU_prospective_02 + IU_prospective_03 + IU_prospective_04 + IU_prospective_05 + IU_prospective_06 + IU_prospective_07

URS_ =~ URS_change_01 + URS_change_02 + URS_change_03 + URS_change_04 + URS_change_05 + URS_change_06 + URS_change_07 + URS_change_08 + URS_change_09 + URS_change_10 + URS_change_11 + URS_change_12 + URS_change_13 + URS_change_14 + URS_change_15 + URS_change_16 + URS_cognitive_01 + URS_cognitive_02 + URS_cognitive_03 + URS_cognitive_04 + URS_cognitive_05 + URS_cognitive_06 + URS_cognitive_07 + URS_cognitive_08 + URS_cognitive_09 + URS_cognitive_10 + URS_cognitive_11 + URS_cognitive_12 + URS_cognitive_13 + URS_cognitive_14 + URS_cognitive_15 + URS_cognitive_16 + URS_cognitive_17 + URS_emotional_01 + URS_emotional_02 + URS_emotional_03 + URS_emotional_04 + URS_emotional_05 + URS_emotional_06 + URS_emotional_07 + URS_emotional_08 + URS_emotional_09 + URS_emotional_10 + URS_emotional_11 + URS_emotional_12 + URS_emotional_13 + URS_emotional_14
URS_des =~ URS_change_01 + URS_change_02 + URS_change_03 + URS_change_04 + URS_change_05 + URS_change_06 + URS_change_07 + URS_change_08 + URS_change_09 + URS_change_10 + URS_change_11 + URS_change_12 + URS_change_13 + URS_change_14 + URS_change_15 + URS_change_16
URS_cog =~ URS_cognitive_01 + URS_cognitive_02 + URS_cognitive_03 + URS_cognitive_04 + URS_cognitive_05 + URS_cognitive_06 + URS_cognitive_07 + URS_cognitive_08 + URS_cognitive_09 + URS_cognitive_10 + URS_cognitive_11 + URS_cognitive_12 + URS_cognitive_13 + URS_cognitive_14 + URS_cognitive_15 + URS_cognitive_16 + URS_cognitive_17
URS_emo =~ URS_emotional_01 + URS_emotional_02 + URS_emotional_03 + URS_emotional_04 + URS_emotional_05 + URS_emotional_06 + URS_emotional_07 + URS_emotional_08 + URS_emotional_09 + URS_emotional_10 + URS_emotional_11 + URS_emotional_12 + URS_emotional_13 + URS_emotional_14

IA_ ~~ 0*IA_abs
IA_ ~~ 0*IA_com
IA_ ~~ 0*IA_dis
NFC_ ~~ 0*IA_abs
NFC_ ~~ 0*IA_com
NFC_ ~~ 0*IA_dis
IU_ ~~ 0*IA_abs
IU_ ~~ 0*IA_com
IU_ ~~ 0*IA_dis
URS_ ~~ 0*IA_abs
URS_ ~~ 0*IA_com
URS_ ~~ 0*IA_dis

NFC_ ~~ 0*NFC_amb
NFC_ ~~ 0*NFC_pre
NFC_ ~~ 0*NFC_ord
NFC_ ~~ 0*NFC_clo
NFC_ ~~ 0*NFC_dec
IA_ ~~ 0*NFC_amb
IA_ ~~ 0*NFC_pre
IA_ ~~ 0*NFC_ord
IA_ ~~ 0*NFC_clo
IA_ ~~ 0*NFC_dec
IU_ ~~ 0*NFC_amb
IU_ ~~ 0*NFC_pre
IU_ ~~ 0*NFC_ord
IU_ ~~ 0*NFC_clo
IU_ ~~ 0*NFC_dec
URS_ ~~ 0*NFC_amb
URS_ ~~ 0*NFC_pre
URS_ ~~ 0*NFC_ord
URS_ ~~ 0*NFC_clo
URS_ ~~ 0*NFC_dec

IU_ ~~ 0*IU_inh
IU_ ~~ 0*IU_pro
IA_ ~~ 0*IU_inh
IA_ ~~ 0*IU_pro
NFC_ ~~ 0*IU_inh
NFC_ ~~ 0*IU_pro
URS_ ~~ 0*IU_inh
URS_ ~~ 0*IU_pro

URS_ ~~ 0*URS_des
URS_ ~~ 0*URS_cog
URS_ ~~ 0*URS_emo
NFC_ ~~ 0*URS_des
NFC_ ~~ 0*URS_cog
NFC_ ~~ 0*URS_emo
IA_ ~~ 0*URS_des
IA_ ~~ 0*URS_cog
IA_ ~~ 0*URS_emo
IU_ ~~ 0*URS_des
IU_ ~~ 0*URS_cog
IU_ ~~ 0*URS_emo

IA_ ~~ r1*NFC_
IA_ ~~ r2*IU_
IA_ ~~ r3*URS_
NFC_ ~~ r4*IU_
NFC_ ~~ r5*URS_
IU_ ~~ r6*URS_

#constraints
r1 < 1; r1 > -1
r2 < 1; r2 > -1
r3 < 1; r3 > -1
r4 < 1; r4 > -1
r5 < 1; r5 > -1
r6 < 1; r6 > -1
'

### Bifactor model with parcels ##############################################################################

mod.par ='
IA_ =~ IA_abs1 + IA_abs2 + IA_abs3 + IA_com1 + IA_com2 + IA_com3 +  IA_dis1 + IA_dis2 + IA_dis3
IA_abs =~ IA_abs1 + IA_abs2 + IA_abs3
IA_com =~ IA_com1 + IA_com2 + IA_com3
IA_dis =~ IA_dis1 + IA_dis2 + IA_dis3

NFC_ =~ NFC_amb1 + NFC_amb2 + NFC_amb3 + NFC_pre1 + NFC_pre2 + NFC_pre3 + NFC_ord1 + NFC_ord2 + NFC_ord3 + NFC_clo1 + NFC_clo2 + NFC_clo3 + NFC_dec1 + NFC_dec2 + NFC_dec3
NFC_amb =~ NFC_amb1 + NFC_amb2 + NFC_amb3
NFC_pre =~ NFC_pre1 + NFC_pre2 + NFC_pre3
NFC_ord =~ NFC_ord1 + NFC_ord2 + NFC_ord3
NFC_clo =~ NFC_clo1 + NFC_clo2 + NFC_clo3
NFC_dec =~ NFC_dec1 + NFC_dec2 + NFC_dec3

IU_ =~ IU_inh1 + IU_inh2 + IU_inhibitory_05 + IU_pro1 + IU_pro2 + IU_pro3
IU_inh =~ IU_inh1 + IU_inh2 + IU_inhibitory_05
IU_pro =~ IU_pro1 + IU_pro2 + IU_pro3

URS_ =~ URS_des1 + URS_des2 + URS_des3 + URS_cog1 + URS_cog2 + URS_cog3 + URS_emo1 + URS_emo2 + URS_emo3
URS_des =~ URS_des1 + URS_des2 + URS_des3
URS_cog =~ URS_cog1 + URS_cog2 + URS_cog3
URS_emo =~ URS_emo1 + URS_emo2 + URS_emo3

IA_ ~~ 0*IA_abs
IA_ ~~ 0*IA_com
IA_ ~~ 0*IA_dis
NFC_ ~~ 0*IA_abs
NFC_ ~~ 0*IA_com
NFC_ ~~ 0*IA_dis
IU_ ~~ 0*IA_abs
IU_ ~~ 0*IA_com
IU_ ~~ 0*IA_dis
URS_ ~~ 0*IA_abs
URS_ ~~ 0*IA_com
URS_ ~~ 0*IA_dis

NFC_ ~~ 0*NFC_amb
NFC_ ~~ 0*NFC_pre
NFC_ ~~ 0*NFC_ord
NFC_ ~~ 0*NFC_clo
NFC_ ~~ 0*NFC_dec
IA_ ~~ 0*NFC_amb
IA_ ~~ 0*NFC_pre
IA_ ~~ 0*NFC_ord
IA_ ~~ 0*NFC_clo
IA_ ~~ 0*NFC_dec
IU_ ~~ 0*NFC_amb
IU_ ~~ 0*NFC_pre
IU_ ~~ 0*NFC_ord
IU_ ~~ 0*NFC_clo
IU_ ~~ 0*NFC_dec
URS_ ~~ 0*NFC_amb
URS_ ~~ 0*NFC_pre
URS_ ~~ 0*NFC_ord
URS_ ~~ 0*NFC_clo
URS_ ~~ 0*NFC_dec

IU_ ~~ 0*IU_inh
IU_ ~~ 0*IU_pro
IA_ ~~ 0*IU_inh
IA_ ~~ 0*IU_pro
NFC_ ~~ 0*IU_inh
NFC_ ~~ 0*IU_pro
URS_ ~~ 0*IU_inh
URS_ ~~ 0*IU_pro

URS_ ~~ 0*URS_des
URS_ ~~ 0*URS_cog
URS_ ~~ 0*URS_emo
NFC_ ~~ 0*URS_des
NFC_ ~~ 0*URS_cog
NFC_ ~~ 0*URS_emo
IA_ ~~ 0*URS_des
IA_ ~~ 0*URS_cog
IA_ ~~ 0*URS_emo
IU_ ~~ 0*URS_des
IU_ ~~ 0*URS_cog
IU_ ~~ 0*URS_emo

IA_ ~~ r1*NFC_
IA_ ~~ r2*IU_
IA_ ~~ r3*URS_
NFC_ ~~ r4*IU_
NFC_ ~~ r5*URS_
IU_ ~~ r6*URS_

#constraints
r1 < 1; r1 > -1
r2 < 1; r2 > -1
r3 < 1; r3 > -1
r4 < 1; r4 > -1
r5 < 1; r5 > -1
r6 < 1; r6 > -1
'
# parcel names
par.names=c("IA_abs1", "IA_abs2", "IA_abs3", "IA_com1", "IA_com2", "IA_com3", "IA_dis1", "IA_dis2", "IA_dis3",
            "NFC_amb1", "NFC_amb2", "NFC_amb3", "NFC_pre1", "NFC_pre2", "NFC_pre3", "NFC_ord1", "NFC_ord2", "NFC_ord3", "NFC_clo1", "NFC_clo2", "NFC_clo3", "NFC_dec1", "NFC_dec2", "NFC_dec3",
            "IU_inh1", "IU_inh2", "IU_pro1", "IU_pro2", "IU_pro3",
            "URS_des1", "URS_des2", "URS_des3", "URS_cog1", "URS_cog2", "URS_cog3", "URS_emo1", "URS_emo2", "URS_emo3")


### Generate parcellated data and do CFA ##############################################################################

# generate parcellated datasets
# must specify data and the number of allocations (nAlloc)
list1=parcelAllocation(mod.par, data=usm, par.names, mod.items, nAlloc=2, do.fit=F, std.lv=T)

# Multifit() conducts a CFA for each data.frame in a list (saved from parcelAllocation) and returns all of the results of interest on one row per data.frame
# The results consist of 21 columns with fit measures ("npar", "chisq", "df", "pvalue", "cfi", "tli", "rmsea", "rmsea.pvalue", "srmr") followed by the latent r between scales ("IA_NFC", "IA_URS", etc.)
# WARNING: running each CFA takes roughly 20 seconds (Intel i5-8350U CPU), so running it for 200 data.frames should take 1-1.5 hours. 

multifit=function(data) { 
  
  rows=c(131:136) # rows indexing the parameter estimates of interest (correlations among general factors)
  names=c("npar", "chisq", "df", "pvalue", "cfi", "tli", "rmsea", "rmsea.pvalue", "srmr",
          "IA_NFC", "IA_IU", "IA_URS", "NFC_IU", "NFC_URS", "IU_URS", "IA_NFC_se", "IA_IU_se", "IA_URS_se", "NFC_IU_se", "NFC_URS_se", "IU_URS_se")
  results=sapply(data, function(x){
    output=cfa(data=x, model=mod.par, std.lv=T) # save CFA output
    est=as.matrix(out@ParTable[["est"]][rows]) # save the parameter estimates of interest
    se=as.matrix(out@ParTable[["se"]][rows])
    fit=as.matrix(fitMeasures(output)[c(1,3,4,5,9,10,23,26,29)]) # save fit measures of interest
    bind=round(do.call("rbind", (list(fit,est,se))),3) # bind the data
  })
  results.t=t(results) # transpose the matrix
  colnames(results.t)=names # name columns
  results.t 
}

# example:
# fit.usm = multifit(list1)

# to find estimates of interest by row, you need to use a dataset with the parcellated items. est1 gives all the estimates
# fit1=cfa(data=usm, model=mod.par, std.lv=T)
# est1=t(rbind(fit1@ParTable[["lhs"]],fit1@ParTable[["op"]],fit1@ParTable[["rhs"]],fit1@ParTable[["est"]]))


### Average estimates and determine N of stability ###############################################################

# Average the CFA results produced by multifit(). This can also be done in Excel
res=read.table(file='clipboard', header=T)  # Copy to clipboard the table consisting of the fit measures and correlations
write.table(t(round(colMeans(abs(res)),3)), file='clipboard', sep='\t', row.names = F) # write to clipboard the averages

# Find at what N (number of unique CFA results) the correlations become stable
# Since the row order (of the CFA results) can produce different Ns, we run the analysis multiple times, each time randomizing the row order.
cors=read.table(file='clipboard', header=T) # Copy to clipboard only the columns with correlations

add=10
x=replicate(1000, for(N in seq(from=10, to=200, by=10)) {
  
  cors=cors[sample(nrow(cors)),] # randomize the row order
  
  d1=colMeans(cors[1:(N+add),])-colMeans(cors[1:N,]) # difference between the averages from the first 20 and first 10 results
  count1=as.numeric(length(d1[abs(d1)>.005])>0) # are any of the differences greater than .005?
  
  d2=colMeans(cors[1:(N+add*2),])-colMeans(cors[1:(N+add),]) # 30-20
  count2=as.numeric(length(d2[abs(d2)>.005])>0) # are any of the differences greater than .005?

  if(count1+count2==0) { # if there are no changes (.01) twice in a row
    return(N)
    break
  }
})
mean(unlist(x), na.rm=T) # the average N

