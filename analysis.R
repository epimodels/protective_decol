global <- as.data.frame(read.csv("decol_simulations.csv",header=TRUE))

colonization <- glm(Colonization.Cases ~ Duration.of.Protection*CHG.Replacement.Effectiveness,data=global,family="poisson")
summary(colonization)

bacteremia <- glm(Bacteremia.Cases ~ Duration.of.Protection*CHG.Replacement.Effectiveness,data=global,family="poisson")
summary(bacteremia)

death <- glm(Bacteremia.Related.Deaths ~ Duration.of.Protection*CHG.Replacement.Effectiveness,data=global,family="poisson")
summary(death)