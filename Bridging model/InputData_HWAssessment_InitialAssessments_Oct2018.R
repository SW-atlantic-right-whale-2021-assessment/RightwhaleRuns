#Data for Humpback Assessment Model - Initial Runs Oct 2018
#Sent to Grant Adams and John Best on 28 Oct 2018
##################################################################


######################
# INPUT DATA - ABUNDANCE AND TREND ESTIMATES
######################

#2005 Absoulute Abundance
Abs.Abundance.2005 <- data.frame(Year=c(2005), 
                                 N.obs=c(6404), 
                                 CV.obs=c(0.12))

#2008 Absolute Abundance
Abs.Abundance.2008 <- data.frame(Year=c(2008), 
                                 N.obs=c(14264), 
                                 CV.obs=c(0.084))

#2012 Absolute Abundance
Abs.Abundance.2012 <- data.frame(Year=c(2012), 
                                 N.obs=c(20389), 
                                 CV.obs=c(0.071))

#Index of Abundance from Branch 2011
Rel.Abundance.Branch=data.frame(
  Index = rep(1,3),
  Year=c(1982, 1986, 1997), 
  IA.obs=c(45, 259, 200), 
  CV.IA.obs=c(0.91, 0.59, 0.64))


#Index of Abundance from Pavanato et al. 2017, Ecol. Modelling, Table 3.
Rel.Abundance.Pavanato <- data.frame(
  Index = rep(1,3),
  Year=c(2008, 2011, 2015), 
  IA.obs=c(7689, 8652, 12123), 
  CV.IA.obs=c(0.078, 0.066, 0.07))

#Index of Abundance from Wedekin et al. 2017, Mar. Ecol. Prog. Series, Table 2.
Rel.Abundance.Wedekin=data.frame(
  Index = rep(1,6),
  Year=c(2002, 2003, 2004, 2005, 2008, 2011), 
  IA.obs=c(3026, 2999, 3763, 4113, 5399, 8832), 
  CV.IA.obs=c(0.132,0.131, 0.179,0.09, 0.136, 0.141))


################
# CATCH SERIES
################

# 
Modern.Empty <- data.frame(Year = seq(1830, 1900), Catch = rep(0, length(seq(1830, 1900))))

#Catch Series Modern (from Zerbini et al., 2011)

Core.Catches <- data.frame(Year=seq(1901,2017), 
                           Catch=c(0, 0, 0, 180, 288, 240, 1261, 1849, 3391, 6468, 5832, 2881, 999, 1155, 
                                   1697, 447, 121, 129, 111, 102, 9, 364, 133, 266, 254, 7, 0, 19, 51, 107, 
                                   18, 23, 132, 57, 48, 105, 242, 0, 2, 36, 13, 0, 4, 60, 238, 30, 35, 48, 
                                   83, 698, 45, 34, 140, 44, 96, 167, 61, 16, 15, 27, 13, 24, 12, 0, 52, 0,
                                   189, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                   0, 0 ,0))
Core.Catches2 <- rbind(Modern.Empty, Core.Catches)

sum(Core.Catches$Catch) # should be 31170

Fringe.Catches <- data.frame(Year=seq(1901, 2017),
                             Catch=c(0, 0, 0, 180, 288, 240, 1261, 1849, 3391, 6468, 5832, 2881, 999, 1155, 
                                     1697, 447, 121, 129, 111, 102, 9, 364, 133, 266, 254, 7, 0, 19, 55.5, 120, 
                                     18.5, 24, 151, 64, 148.5, 149, 274.5, 0, 2, 91.5, 13, 0, 4, 60, 238, 30.5, 35.5, 67, 
                                     211.5, 712, 102.5, 50.5, 155.5, 70, 137.5, 199.5, 77.5, 19, 18.5, 29, 13, 26,
                                     12, 0, 69, 0, 192, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                     0, 0, 0, 0, 0, 0))
Fringe.Catches2 <- rbind(Modern.Empty, Fringe.Catches)

sum(Fringe.Catches$Catch) # should be 31846.5

Overlap.Catches <- data.frame(Year=seq(1901, 2017),
                              Catch=c(0, 0, 0.1, 144, 232.7, 241.8, 1045.4, 1604.7, 2870.3, 5434.1, 4891.6, 2471.6, 974.1, 
                                      1054.2, 1395.9, 373.3, 116.1, 123.7, 112.7, 96.5, 7.2, 310.1, 116, 223, 219.5, 15.7, 
                                      0.4, 16.7, 42.4, 92, 14.7, 20.2, 114.1, 49.4, 68.4, 109.1, 212.5, 0, 1.6, 53.1, 10.4, 
                                      0, 3.2, 48, 190.4, 24.1, 30.3, 51.4, 116.3, 613.7, 84.4, 49.1, 123.7, 70.5, 94.3, 209.5, 
                                      61.4, 27.7, 39.7, 45.4, 131.5, 53.3, 12.2, 0, 133, 14.7, 225.5, 0, 0, 0, 0, 1.7, 0, 0, 0, 
                                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
Overlap.Catches2 <- rbind(Modern.Empty, Overlap.Catches)

sum(Overlap.Catches$Catch) # should be 27334.3

Falkland.Catches <- data.frame(Year=seq(1901, 2017),
                               Catch=c(0, 0, 0, 0, 0, 0, 0, 6, 66, 49, 12, 6, 5, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 
                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 36, 
                                       0, 4, 1, 22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

sum(Falkland.Catches$Catch) # should be 219
Falkland.Catches2 <- rbind(Modern.Empty, Falkland.Catches)


#Catch Series Pre-Modern (from Morais et al., 2017)
PreModern.Catch.Min <- data.frame(Year=seq(1830,1924), 
                                  Catch=c(rep(120, 10), #1830-1839
                                          rep(122.8, 10), #1840-1849
                                          rep(120, 10), #1850-1859
                                          rep(138.1, 10), #1860-1869
                                          rep(120, 24), #1870-1893
                                          168, #1894
                                          rep(120, 6), #1895-1900
                                          rep(271.38, 2), #1901-1902
                                          120, #1903
                                          271.38, #1904 (first year with pre-modern + modern)
                                          271.38, #1905 (pre-modern + modern)
                                          120, #1906
                                          120, #1907
                                          338.87, #1908
                                          190.02, #1909
                                          205.52, #1910
                                          rep(30, 14))) #1911-1924

sum(PreModern.Catch.Min$Catch)

PreModern.Catch.Max <- data.frame(Year=seq(1830,1924), 
                                  Catch=c(rep(400, 10), #1830-1839
                                          rep(402.8, 10), #1840-1849
                                          rep(400, 10), #1850-1859
                                          rep(418.1, 10), #1860-1869
                                          rep(400, 24), #1870-1893
                                          448, #1894
                                          rep(400, 6), #1895-1900
                                          rep(581.65, 2), #1901-1902
                                          400, #1903
                                          581.65, #1904 (first year with pre-modern + modern)
                                          581.65, #1905 (pre-modern + modern)
                                          400, #1906
                                          400, #1907
                                          406.65, #1908
                                          228.03, #1909
                                          246.62, #1910
                                          rep(50, 14))) #1911-1924

sum(PreModern.Catch.Max$Catch)

###########################
# STRUCK AND LOST RATES
###########################

# Pre-modern whaling period (from Smith and Reeves, 2010) 
# Multiplier for pre-modern whaling catches is N(1.71, 0.073) 

# Modern whaling between 1900 and 1939 (from Best, 2010)
# A distribution that has a 5% the probability that the SLR is greater than 13.9% 
# Multiplier for this period is a distribution for which there is a 5% probability of a value greater than 1.16 and a maximum value of 1.42

# Modern whaling between 1939 and 1945 (from Tonessen and Johnsen, 1982 and Best, 2010)
# Multiplier is U(1.25, 1.42)

# Modern whaling after 1945 (from Smith and Reeves, 2010)
# Multiplier is N(1.0185, 0.0028)


##########################
# GENETIC CONSTRAINT
##########################

# This is seen as a lower bound on the minimum population. It is calculated as the number of haplotypes for the population times 3 and serves as a conservative estimate of the number of whales in the population when it went through a bottleneck after severe exploitation. For now, lets use the number of haplotypes in Zerbini et al. (2011), which is 66, and multiply it by 3 to get the genetic constraint value (of 198). When the number of haplotypes is updated, we can re-run with the right number.





