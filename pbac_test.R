library( knitr )
library( MASS )
library( ANTsR )
library( oro.nifti )
library( visreg )
library( ggplot2 )
library( boot )
library( candisc )
library( pheatmap )
library( ANTsR )
rdir="./RMI_Data"
rdir="./"
pbacc<-list.files(path = rdir, pattern = glob2rx("pbac*csv")  ,
full.names = T )
pbacTEcog<-read.csv(pbacc[1])
pbacTRcog<-read.csv(pbacc[2])
pbaci<-list.files(path = rdir, pattern = glob2rx("pbac*mha") ,
full.names = T )
pbacTEimg<-as.matrix( antsImageRead(pbaci[1], 2 ) )
pbacTRimg<-as.matrix( antsImageRead(pbaci[2], 2 ) )


# also pbac imaging data comes from this mask
mask<-antsImageRead( list.files(path = rdir, pattern=glob2rx("gmask_2mmb.nii.gz") , full.names=T ) , 3 )
# with anatomical labels
pbacaal<-antsImageRead( list.files(path = rdir,
pattern=glob2rx("pbac_aal.nii.gz"), full.names=T ) , 3 )
data("aal",package="ANTsR") # description of aal

#####################
inmask <-  mask > 0.5
mylabs<-sort( unique( pbacaal[ inmask  &  pbacaal > 0.5 &  pbacaal <
91   &  pbacaal != 51 &  pbacaal != 52 &  pbacaal != 53 &  pbacaal != 54 ] ) )
roimatrix<-matrix( rep( NA, length( mylabs ) * nrow( pbacTRimg ) ) ,
ncol=length(mylabs ) )
for ( i in 1:length(mylabs) ) {
# get vector for this label
  labelVec <- as.numeric( pbacaal[ inmask ] == mylabs[ i ] )
  roimatrix[   , i  ] <- pbacTRimg %*% ( labelVec / sum( labelVec ) )
  }
colnames( roimatrix ) <- aal$label_name[ mylabs ]
mydf<-data.frame( pbacTRcog, roimatrix )


## ----out.width='.7\\linewidth',dev='pdf',eval=TRUE,results='show',echo=TRUE----
pheatmap( cor( pbacTRcog ) , cluster_rows = F , cluster_cols =  F )


## ----out.width='.75\\linewidth',dev='pdf',eval=TRUE,results='show',echo=FALSE----
stars( pbacTRcog,
       labels = row.names(mydf), cex = 0.2, scale = TRUE, radius = FALSE, full = TRUE, flip.labels = FALSE,
       mar = c( 0, 0, 2, 0 ), main = "Brain Constellation Map of PBAC Cognition" )


pheatmap( cor( roimatrix ) , cluster_rows = F , cluster_cols =  F )


stars( roimatrix,
       labels = row.names(mydf), cex = 0.2, scale = TRUE, radius = FALSE, full = TRUE, flip.labels = FALSE,
       mar = c( 0, 0, 2, 0 ), main = "Brain Constellation Map of PBAC ROIs" )


## ----out.width='.7\\linewidth',dev='pdf',eval=TRUE,results='show',echo=TRUE----
pheatmap( cor( pbacTRcog , roimatrix ) , cluster_rows = F , cluster_cols =  F )


## ----out.width='.4\\linewidth',dev='jpeg',eval=TRUE,results='hide',echo=TRUE----
myform<-paste( colnames( roimatrix ), collapse='+'  )
myform<-as.formula( paste( "delay_free_adj~", myform , "+edu") )
mydf<-data.frame( pbacTRcog, roimatrix )
row.names(mydf)<- paste( c(1:nrow(pbacTRcog)),"_",as.character( pbacTRcog$mmse ),sep='')
mdl <- lm(  myform , data = mydf )
mdla<-stepAIC( mdl , direction =  c("forward" ) , k = 20 , steps= 20  )
ageregions<-gsub("_","",as.character(mdla$call$formula)[3])


## ----out.width='.35\\linewidth',dev='pdf',eval=TRUE,results='show',echo=TRUE----
visreg( mdla, xvar="Angular_L")
visreg( mdla, xvar="Frontal_Mid_R")
visreg( mdla, xvar="Temporal_Pole_Sup_L")


## ----out.width='.6\\linewidth',dev='pdf',echo=TRUE,eval=TRUE,results='hide'----
 coplot( delay_free_adj ~
   Angular_L + Frontal_Mid_R + Temporal_Pole_Sup_L | age ,
   data = mydf , panel = panel.smooth, rows = 1)


## ----out.width='.6\\linewidth',dev='pdf',echo=2,eval=TRUE,results='hide'----
nv<-5
m1=scale( as.matrix(pbacTRcog) ); # m1=m1-min(m1)
m2=scale( as.matrix(pbacTRimg) ); # m2=m2-min(m2)
mysccan<-sparseDecom2(
  inmatrix=list(m1 , m2),
  inmask=c( NA , mask ),
  smooth = 0, its=50,
  mycoption = 0, verbose=1,
  sparseness=c( -0.1, -0.05 ),
  nvecs=nv, perms=0, cthresh=c(0,250), ell1=1 )

sumcor=0
for ( ind in 1:ncol(pbacTRcog) ) {
#  myform<-paste( mytests , collapse="+" )
  lowmat = pbacTRimg %*% mysccan$eig2
  traindf<-data.frame(  pbacTRimg %*% mysccan$eig2,
#    cog = data.matrix(pbacTRcog) %*% as.numeric( mysccan$eig1[,ind] ) )
    cog = pbacTRcog[,ind] )
  testdf <-data.frame(  pbacTEimg %*% mysccan$eig2,
#    cog = data.matrix(pbacTEcog) %*% as.numeric( mysccan$eig1[,ind] ) )
    cog = pbacTEcog[,ind] )
  myform<-as.formula( paste( "gm~",myform) )
  predlm<-lm( cog ~ . , data=traindf )
  predcog<-predict( predlm , newdata=traindf )
  predcog2<-predict( predlm , newdata=testdf )
  print(  colnames( pbacTRcog )[ind] )
  print( paste("Train Correlation:",ind, cor.test( predcog, traindf$cog )$est  ) )
  print( paste("Test Correlation:",ind, cor.test( predcog2, testdf$cog )$est  ) )
  sumcor=sumcor+cor.test( predcog2, testdf$cog )$est
}
print(sumcor) # 7.57058  vs 6 ...
