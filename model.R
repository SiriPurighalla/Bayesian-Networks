#This codes helps on the construction of a Gaussian Graphical model to study
#the correlation of hormonal contraceptives and cervicl cander diagnosis.

library( lavaan )
library(dagitty)

#Build the dag 

g <- dagitty('
 dag {
    Age [pos="0.312,-1.570"]
    Dx.Cancer [outcome,pos="-0.506,1.446"]
    Dx.HPV [pos="0.318,1.348"]
    First.sexual.intercourse [pos="-1.452,-1.647"]
    Hormonal.Contraceptives..years. [exposure,pos="0.044,-0.199"]
    IUD..years. [pos="-0.523,0.528"]
    Num.of.pregnancies [pos="-0.552,-0.142"]
    Number.of.sexual.partners [pos="-1.535,-0.327"]
    STDs [pos="-1.497,0.694"]
    Smokes [pos="0.354,-0.384"]
    Age -> First.sexual.intercourse
    Age -> Hormonal.Contraceptives..years.
    Age -> IUD..years.
    Age -> Num.of.pregnancies
    Age -> Number.of.sexual.partners
    Age -> Smokes
    Dx.HPV -> Dx.Cancer
    First.sexual.intercourse -> Hormonal.Contraceptives..years.
    First.sexual.intercourse -> IUD..years.
    First.sexual.intercourse -> Num.of.pregnancies
    First.sexual.intercourse -> Number.of.sexual.partners
    Hormonal.Contraceptives..years. -> Dx.HPV
    Hormonal.Contraceptives..years. -> Num.of.pregnancies
    IUD..years. -> Dx.Cancer
    IUD..years. -> Num.of.pregnancies
    Num.of.pregnancies -> Dx.HPV
    Number.of.sexual.partners -> Num.of.pregnancies
    Number.of.sexual.partners -> STDs
    STDs -> Dx.HPV
    Smokes -> Dx.HPV
  }
')

plot(g)

#Load the document
d <- new_df

#Local test 
localTests(g, d, type = 'cis')

#Covariance Matrix and local Test 
M <- lavCor(d)
r <- localTests(g, sample.cov=M, sample.nobs=nrow(d))

plotLocalTestResults( r, xlim=c(-1,1) )

#Correlation analysis 

info = c('Dx.Cancer ~ IUD..years. + Dx.HPV', 
         'Dx.HPV ~ Hormonal.Contraceptives..years. + STDs + Smokes + Num.of.pregnancies' , 
         'STDs ~ Number.of.sexual.partners' , 
         'Num.of.pregnancies ~ Age + Hormonal.Contraceptives..years. + Number.of.sexual.partners + IUD..years. + First.sexual.intercourse',
         'Smokes ~ Age' ,
         'First.sexual.intercourse ~ Age' ,
         'Number.of.sexual.partners ~ First.sexual.intercourse + Age')

res <- list()

for(d in info){
  mdl <-lm( d,data=data.frame(scale( new_df ) ) )
  c = coef(mdl)
  for(i in 2:length(c)){
    res = append(res, c)
  }
  print(length(c))
}

print(res)       

#Adjustment sets
adjustmentSets(g, 'Hormonal.Contraceptives..years.', 'Dx.Cancer'
               , type = "canonical", effect = "total")
