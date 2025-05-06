
### R code for the publication in JAMES
# Develope by Cuihai You


#############################
##Gelman-Rubin statistic calculation
chain1<-read_excel("****/para_chain.xlsx",1)
chain2<-read_excel("****/para_chain.xlsx",2)
chain3<-read_excel("****/para_chain.xlsx",3)
chain4<-read_excel("****/para_chain.xlsx",4)

chain1<-as.mcmc(chain1)
chain2<-as.mcmc(chain2)
chain3<-as.mcmc(chain3)
chain4<-as.mcmc(chain4)
y<-mcmc.list(chain1,chain2,chain3,chain4)

#根据est. 判断，你可以help一下，再看看这个功能，我这个代码不是很全了
gelman.diag(y, confidence = 0.95, transform=FALSE, autoburnin=TRUE,
            multivariate=TRUE)


############################
######Slope compare
library(readxl)
data<-read_excel("***/your_data.xlsx",1)
lm1<-lm(Simulated~Measured-1,data)
summary(lm1)
confint(lm1, level=0.95)## includes 1, indacating the slope has no difference with 1



#########################
### anova and multi-compare of parameters

data<-read_excel("*****/para_normal_dry_wet_linear.xlsx")
# ANOVA
fit <- aov(turnp~ Treat, data = data)
summary(fit)  # This gives the summary of the ANOVA 
# Post Hoc Test (Tukey's test for multiple comparisons)
posthoc <- TukeyHSD(fit)
print(posthoc)





######################
###PCA analysis
library("FactoMineR")
library("factoextra")

rm(list=ls())
data<-read_excel("*******/para_normal_dry_wet_linear.xlsx",1)##????????

decathlon2.active <- data[1:12, 2:9]
head(decathlon2.active[, 1:6], 4)
##：
PCA(decathlon2.active, scale.unit = TRUE, ncp = 5, graph = TRUE)
##
res.pca <- PCA(decathlon2.active, graph = TRUE)
print(res.pca)
##：
head(var$contrib, 4)
#picture of PCA
my_color <- c("#e41a1c", "#252525", "#377eb8")
p1<-fviz_pca_biplot(res.pca,arrowsize = 0.8,
                    col.ind = data$Treat, palette = my_color,  pointsize = 4, 
                    addEllipses = TRUE, label = "var",mean.point = FALSE,
                    col.var = "#252525", repel = TRUE,geom.var = c("arrow"))+#
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"),
        axis.title = element_text(size = 32), 
        axis.text = element_text(size = 28),
        legend.text = element_text(size = 16),
        legend.spacing.y = unit(3, 'cm'),
        legend.key.size= unit(1.2, 'cm'),
        legend.title = element_blank()) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2));p1
ggsave(filename = "*****/PCA_linear.jpeg", plot = p1, device = "jpg", width = 9, height = 8,  dpi = 300)

