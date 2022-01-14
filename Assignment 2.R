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

plot(g,show.coefficients = TRUE)

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

adjustmentSets(g, 'Hormonal.Contraceptives..years.', 'Dx.Cancer', type = "all", effect = "total")

fit <- sem( toString(g,"lavaan"), sample.cov=M, sample.nobs=nrow(d) )

fg <- lavaanToGraph(fit, digits=2)
coordinates(fg) <- res
plot(fg, show.coefficients=TRUE)

install.packages('bnlearn')
library(bnlearn)


# Data pre processing for Structure Learning
updated_data <- new_df



# Converting Age into Age groups
s <- rep("<10", nrow(updated_data))
s[updated_data$Age >= 10 & updated_data$Age < 20] <- "10-20"
s[updated_data$Age >= 20 & updated_data$Age < 30] <- "20-30"
s[updated_data$Age >= 30 & updated_data$Age < 40] <- "30-40"
s[updated_data$Age >= 40 & updated_data$Age < 50] <- "40-50"
s[updated_data$Age >= 50 & updated_data$Age < 60] <- "50-60"
s[updated_data$Age >= 60 & updated_data$Age < 70] <- "60-70"
s[updated_data$Age >= 70 & updated_data$Age < 80] <- "70-80"
s[updated_data$Age >= 80] <- ">80"

updated_data$Age <- as.factor(s)

# Converting First Sexual Intercourse into caergorical
fsi <- rep("<10", nrow(updated_data))
fsi[updated_data$First.sexual.intercourse >= 10 & updated_data$First.sexual.intercourse < 20] <- "10-20"
fsi[updated_data$First.sexual.intercourse >= 20 & updated_data$First.sexual.intercourse < 30] <- "20-30"
fsi[updated_data$First.sexual.intercourse >= 30 & updated_data$First.sexual.intercourse < 40] <- "30-40"
fsi[updated_data$First.sexual.intercourse >= 40 & updated_data$First.sexual.intercourse < 50] <- "40-50"
updated_data$First.sexual.intercourse <- as.factor(fsi)

# Converting No. of Sexual Partners into categorical
nsp <- rep("<0", nrow(updated_data))
nsp[updated_data$Number.of.sexual.partners >= 0 & updated_data$Number.of.sexual.partners < 5] <- "0-5"
nsp[updated_data$Number.of.sexual.partners >= 5 & updated_data$Number.of.sexual.partners < 10] <- "5-10"
nsp[updated_data$Number.of.sexual.partners >= 10 & updated_data$Number.of.sexual.partners < 15] <- "10-15"
nsp[updated_data$Number.of.sexual.partners >= 15 & updated_data$Number.of.sexual.partners < 20] <- "15-20"
nsp[updated_data$Number.of.sexual.partners >= 20 & updated_data$Number.of.sexual.partners < 25] <- "20-25"
nsp[updated_data$Number.of.sexual.partners >= 25 & updated_data$Number.of.sexual.partners < 30] <- "25-30"

updated_data$Number.of.sexual.partners <- as.factor(nsp)

# Converting No. of Pregnancies into categorical
np <- rep("<0", nrow(updated_data))
np[updated_data$Num.of.pregnancies >=0 & updated_data$Num.of.pregnancies < 2] <- '0-2'
np[updated_data$Num.of.pregnancies >=2 & updated_data$Num.of.pregnancies < 4] <- '2-4'
np[updated_data$Num.of.pregnancies >=4 & updated_data$Num.of.pregnancies < 6] <- '4-6'
np[updated_data$Num.of.pregnancies >=6 & updated_data$Num.of.pregnancies < 8] <- '6-8'
np[updated_data$Num.of.pregnancies >=10 & updated_data$Num.of.pregnancies < 12] <- '8-10'

updated_data$Num.of.pregnancies <- as.factor(np)

# Converting Hormonal Contraceptive Years
hcy <- rep("<0", nrow(updated_data))
hcy[updated_data$Hormonal.Contraceptives..years. >=0 & updated_data$Hormonal.Contraceptives..years. < 5] <- '0-5'
hcy[updated_data$Hormonal.Contraceptives..years. >=5 & updated_data$Hormonal.Contraceptives..years. < 10] <- '5-10'
hcy[updated_data$Hormonal.Contraceptives..years. >=10 & updated_data$Hormonal.Contraceptives..years. < 15] <- '10-15'
hcy[updated_data$Hormonal.Contraceptives..years. >=15 & updated_data$Hormonal.Contraceptives..years. < 20] <- '15-20'
hcy[updated_data$Hormonal.Contraceptives..years. >=20 & updated_data$Hormonal.Contraceptives..years. < 25] <- '20-25'
hcy[updated_data$Hormonal.Contraceptives..years. >=25 & updated_data$Hormonal.Contraceptives..years. < 30] <- '25-30'
hcy[updated_data$Hormonal.Contraceptives..years. >=30 & updated_data$Hormonal.Contraceptives..years. < 35] <- '30-35'

updated_data$Hormonal.Contraceptives..years. <- as.factor(hcy)


# Converting IUD Years
iud <- rep("<0", nrow(updated_data))
iud[updated_data$IUD..years. >=0 & updated_data$IUD..years. < 5] <- '0-5'
iud[updated_data$IUD..years. >=5 & updated_data$IUD..years. < 10] <- '5-10'
iud[updated_data$IUD..years. >=10 & updated_data$IUD..years. < 15] <- '10-15'
iud[updated_data$IUD..years. >=15 & updated_data$IUD..years. < 20] <- '15-20'

updated_data$IUD..years. <- as.factor(iud)

# Converting Average Smoking years
smoke <- rep("<0", nrow(updated_data))
smoke[updated_data$Smokes >=0 & updated_data$Smokes < 200] <- '0-200'
smoke[updated_data$Smokes >=200 & updated_data$Smokes < 400] <- '200-400'
smoke[updated_data$Smokes >=400 & updated_data$Smokes < 600] <- '400-600'
smoke[updated_data$Smokes >=600 & updated_data$Smokes < 800] <- '600-800'
smoke[updated_data$Smokes >=800 & updated_data$Smokes < 1000] <- '800-1000'
smoke[updated_data$Smokes >=1000 & updated_data$Smokes < 1200] <- '1000-1200'
smoke[updated_data$Smokes >=1200 & updated_data$Smokes < 1400] <- '1200-1400'

updated_data$Smokes <- as.factor(smoke)

updated_data$STDs <- as.factor(updated_data$STDs)
updated_data$Dx.Cancer <- as.factor(updated_data$Dx.Cancer)
updated_data$Dx.HPV <- as.factor(updated_data$Dx.HPV)
 
# Applying PC Algorithm
pc_fit <- pc.stable(updated_data,alpha=0.05, debug = TRUE, test='mi')
plot(pc_fit)
pc_fit

for( v in setdiff(names(d),c("Hormonal.Contraceptives..years.","Dx.HPV")) ){
  print( v )
  print( ci.test("Hormonal.Contraceptives..years.", "Dx.Cancer", v, data=d, test="cor")$p.value )
}

ci.test("Hormonal.Contraceptives..years.", "Dx.Cancer","Age",data=updated_data, test="x2")
labels <- updated_data[,c("Age","Hormonal.Contraceptives..years.","IUD..years.","Dx.HPV","Dx.Cancer","Smokes")]
plot(pc.stable(labels, test = 'x2'))

# Applying HC algorithm
algo2 <- hc(updated_data, start = NULL, whitelist = NULL, blacklist = NULL, score = NULL,
            debug = TRUE, perturb = 1, max.iter = Inf, maxp = Inf, optimized = TRUE)
plot(algo2)
algo2

compare(g, pc_fit)

hamming(pc_fit, algo2)
shd(pc_fit, algo2)    

