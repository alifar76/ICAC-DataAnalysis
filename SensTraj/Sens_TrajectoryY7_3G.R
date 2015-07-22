#' Setwd
setwd("S:/RhoFED/ICAC/Studies/URECA/Prog/Derive/TrajectoriesY7")

#' hadleyverse Pakages
pacman::p_load(magrittr,dplyr,tidyr,broom,reshape2,lubridate,stringr,rio,devtools)
#' Figures Packages
pacman::p_load(lattice,latticeExtra,grid,RColorBrewer)
#' Models
pacman::p_load(lcmm)
#' Session Info
print(session_info(),locale = F)

#` Load Data
a <- import("Data_TrajectoryY7.csv")

#' Reshape
a1 <- a %>%
  select(studyid,asthma_agey7,all_aerosens_y2,all_aerosens_y3,all_aerosens_y5,all_aerosens_y7) %>%
  gather(variable,atopy,all_aerosens_y2:all_aerosens_y7) %>%
  mutate( year = extract_numeric(variable),
          studyid = as.factor(studyid)) %>%
  select(-variable) %>%
  arrange(studyid,year)

#' Models (3 groups)
m3<-hlme(atopy ~ year + I(year^2), random=~ year + I(year^2), mixture=~ year + I(year^2),
         subject='studyid', ng=3, idiag=F, nwg=T, data=a1)

#' Probabilities
t3 <- postprob(m3)

#' Create Data
dl3 <- data_frame(studyid=levels(a1$studyid),class=m3$pprob$class)
dl3 <- left_join(a1,dl3)
ml3 <- dl3 %>%
  filter(year==7) %>%
  mutate(asthma_agey7 = as.numeric(asthma_agey7)) %>%
  group_by(class) %>%
  summarize(m=mean(asthma_agey7)*100)

#' Create Data
class <- data_frame(studyid=dl3[dl3$year==7,'studyid'],
                    c3=m3$pprob[,2],
                    p3c1=m3$pprob[,3], p3c2=m3$pprob[,4], p3c3=m3$pprob[,5]) %>%
  mutate_each(funs(round(.,6)),p3c1,p3c2,p3c3) %>%
  mutate(g3 = car::recode(c3,"1='Late Onset'; 2='Minimal/None'; 3='Early Onset'")) %>%
  rowwise() %>%
  mutate(p3 = max(p3c1,p3c2,p3c3)) %>%
  select(studyid,c3,g3,p3,p3c1,p3c2,p3c3) %>%
  export("Sens_TrajectoryY7_3G.csv")

#' PDF
pdf("Sens_TrajectoryY7_3G.pdf",width=10,height=8)

with(dl3,
     xyplot(atopy~year,
            group=class,
            sub=paste("BIC",round(m3$BIC,0)),
            type="smooth", lwd=2, span=1,
            xlim=c(-0.2,7.5),
            ylab=list("% Atopic",cex=1.8),
            xlab=list("Year",cex=2),
            scales=list(y=list(alternating=3),
                        x=list(at=c(1,2,3,4,5,6,7),tck=c(1,0)),cex=1.5),
            panel = panel.superpose,
            panel.groups = function(x,y,...,subscripts,group.number) {
              panel.smoother(x,y,...)
              d <- tapply(y,x,mean)
              panel.text(x=1,y=d[4],adj=0,
                         label=paste0("G",group.number," ",round(t3[[1]][2,group.number],1),"% \n",
                                      "A  ",round(ml3[group.number,2],1),"%"),
                         col=trellis.par.get('superpose.line')$col[group.number],
                         cex=1.1)
            }))

dev.off()
