##
## Code for "Networks of care for the modern adolescent"
##
## Simon R. White, Emma Soneson, OxWell Study Team, and Mina Fazel (2024). "Networks of care for the modern adolescent".
##
## Note: The data are not publicly available because of ethical and information governance restrictions
##



library("igraph")
library("RColorBrewer")

library("pscl")
library("sjPlot")
library("lmtest")
library("AICcmodavg")



DIRS <- list()

DIRS$USER.SSHKEY <- "~/.ssh" ## Location of ssh-key to decrypt data (must not be on secure drive!)

DIRS$MOUNT <- "OxWell" ## Mount location of secure drive
DIRS$WITHIN <- "networks_of_care" ## Project folder

DIRS$DATA      <<- file.path( DIRS$MOUNT, DIRS$WITHIN, "data" )
DIRS$PROJECT   <<- file.path( DIRS$MOUNT, DIRS$WITHIN, "scripts" )
DIRS$OUTPUTS   <<- file.path( DIRS$MOUNT, DIRS$WITHIN, "outputs" )


##
## The following line decrypts the data-key (using your ssh-rsa keyphrase)
KEY <- cyphr::data_key(path_data=DIRS$DATA,path_user=DIRS$USER.SSHKEY)

if( 0 ) {
    ## Test that encryption/decryption is working: my favourite game is?
    cyphr::decrypt_string(data=file.path(DIRS$DATA,"test.string-cyphr"), key=KEY )
}

## Using OxWell2023 Release 7 (2024-03-18)
DATA <- cyphr::decrypt_object(key=KEY, data=file.path(DIRS$DATA,"DATA2023.r07.limited.rds-cyphr") )
##

DIRTAG  <- "Final"
FILETAG <- "Final"
if( ! dir.exists(file.path(DIRS$OUTPUTS,DIRTAG)) ) {
    dir.create(file.path(DIRS$OUTPUTS,DIRTAG))
}
PNG.SCALER <- 1.25
PNG.BASE <- list(W=1200*PNG.SCALER,H=900*PNG.SCALER,R=200,SCALER=PNG.SCALER)
rm(PNG.SCALER)


if( 1 ) {
    ##
    ## Create subset for analysis
    ##

    ##
    SUBSET <- DATA[ (LASTPAGE.R6a>=34) & (FAB.FLAGCOMBINED=="PASS") & (YEARGROUP %in% sprintf("Y%02i",7:13)) ]
    ##
    ## * Subset to participants that respond to items on or beyond pages with service use items (ie page 34)
    ## * Subset to those who consent and engage with survey for 10+ minutes
    ## * Subset to Y07-Y13
    ##

    
    if( 1 ) {
        ## Sample summary table
        ##
        SUMMARY.DT <- DATA[ ,.N,by=.((!is.na(LASTPAGE.R6a))&(LASTPAGE.R6a>=34),(FAB.FLAGCOMBINED=="PASS"),(YEARGROUP %in% sprintf("Y%02i",7:13))) ][,Total:=sum(N)][,Flag:=sum(N),by=.(FAB.FLAGCOMBINED)][order(FAB.FLAGCOMBINED,YEARGROUP)]
        SUMMARY.DT[order(`FAB.FLAGCOMBINED`,`is.na`,`YEARGROUP`),]
        write.table(SUMMARY.DT,row.names=FALSE,quote=FALSE,sep="|",file=file.path(DIRS$OUTPUTS,DIRTAG,sprintf("%s.summary-sample.csv",FILETAG)))
    }

}

if( 1 ) {
    ##
    ## Clean and process subset
    ##
    
    SUBSET[,YEARGROUP:=droplevels(YEARGROUP)]
    
    SUBSET[,XID:=factor(1:.N)]

    MH.TRANSLATE <- rbindlist( list(
        data.table(old=c("No","Yes - more than a year ago"),new="No"),
        data.table(old=c("Yes - in the past 12 months"),new="Yes"),
        data.table(old=c("Prefer not to say", "NoResponse", "Skipped", "Incomplete(likely)","Incomplete(uncertain)"),new="Unknown")
    ) )
    SUBSET[ MH.TRANSLATE, on=.(X2160=old), MH.PROBLEM:=i.new ]
    if( 0 ) {
        SUBSET[,.N,by=.(X2160,MH.PROBLEM)]
    }
    rm(MH.TRANSLATE)
    
    ## 20240320
    ## Release 7 encode GD/GND (will be added to subsequent release)
    if ( 1 ) {
        SUBSET[,GENDER:=GENDER.broad]
        SUBSET[X1020=="Male",GENDER:="Boy"]
        SUBSET[X1020=="Female",GENDER:="Girl"]
        SUBSET[X1020=="Prefer not to say",GENDER:="Prefer not to say"]
        SUBSET[X1020=="Other" & X1021F=="NULL",GENDER:="Gender diverse"]
        SUBSET[X1020=="NoResponse",GENDER:="NoResponse"]

        Sex.TRANSLATE <- rbindlist( list(
            data.table(old="Boy",new="Boy"),
            data.table(old="Girl",new="Girl"),
            data.table(old="Gender diverse",new="GD/GND"),
            data.table(old="Prefer not to say",new="GD/GND"),
            data.table(old=c("Silly/malicious/not gender","Unsure"),new=NA_character_)
        ))
        SUBSET[ Sex.TRANSLATE, on=.(GENDER=old), Sex:=i.new ]
        Sex2.TRANSLATE <- rbindlist( list(
            data.table(old="Boy",new="Girl/Boy"),
            data.table(old="Girl",new="Girl/Boy"),
            data.table(old="Gender diverse",new="GD/GND"),
            data.table(old="Prefer not to say",new="GD/GND"),
            data.table(old=c("Silly/malicious/not gender","Unsure"),new=NA_character_)
        ))
        SUBSET[ Sex2.TRANSLATE, on=.(GENDER=old), Sex2:=i.new ]
    }
    
    SUBSET[,Yeargroup:=c(7,8,9,10,11,12,13)[YEARGROUP]]
    SUBSET[,YeargroupZ:=(c(7,8,9,10,11,12,13)-7)[YEARGROUP]]

    ## RCADS-11 cut-points
    ##
    SUBSET[Sex=="Girl",RCADS2:=cut(RCADS11.Total.MeanImputed,breaks=c(-Inf,14,Inf),labels=c("Normal","Clinical"),right=FALSE)]
    SUBSET[Sex=="Boy",RCADS2:=cut(RCADS11.Total.MeanImputed,breaks=c(-Inf,9,Inf),labels=c("Normal","Clinical"),right=FALSE)]
    SUBSET[Sex=="GD/GND",RCADS2:=cut(RCADS11.Total.MeanImputed,breaks=c(-Inf,14,Inf),labels=c("Normal","Clinical"),right=FALSE)]


    if( 1 ) {
        ## Sample summary table
        ##
        SUM <- list()
        SUM$A <- SUBSET[,.N,by=.(Sex,RCADS2,MH.PROBLEM)][,Prop:=paste0(format(100*N/sum(N),digits=2,nsmall=1),"%"),by=.(Sex,RCADS2)][order(Sex,RCADS2,MH.PROBLEM)][is.na(Sex),Sex:="NoResponse"][is.na(RCADS2),RCADS2:="Unable to impute"][]
        SUM$B <- SUBSET[,.N,by=.(Sex)][,Prop:=paste0(format(100*N/sum(N),digits=2,nsmall=1),"%")][order(Sex)][is.na(Sex),Sex:="NoResponse"][]
        SUM$C <- SUBSET[,.N,by=.(RCADS2)][,Prop:=paste0(format(100*N/sum(N),digits=2,nsmall=1),"%")][order(RCADS2)][is.na(RCADS2),RCADS2:="Unable to impute"][]
        SUM$D <- SUBSET[,.N,by=.(MH.PROBLEM)][,Prop:=paste0(format(100*N/sum(N),digits=2,nsmall=1),"%")][order(MH.PROBLEM)][]
        SUM$E <- data.table(Sex="Total",N=SUBSET[,.N])

        write.table(rbindlist(SUM[c("A","B","C","D","E")],use.names=TRUE,fill=TRUE),na="",
                    row.names=FALSE,quote=FALSE,sep="|",file=file.path(DIRS$OUTPUTS,DIRTAG,sprintf("%s.summary-breakdown.csv",FILETAG)))
        
    }

    MAP1 <- list("CARER"=c(shortlabel="Carer",label="Parent, step-parent or carer",item="X2190",wish="X2220a",section="informal"),
                 "SIBLING"=c(shortlabel="Sibling(s)",label="Sibling(s)",item="X2191",wish="X2220b",section="informal"),
                 "FAMILY"=c(shortlabel="Family",label="Someone else in your family",item="X2192",wish="X2220c",section="informal"),
                 "FRIPERSON"=c(shortlabel="Friends\nIn person",label="Friend(s), mainly known in person",item="X2193",wish="X2220d",section="informal"),
                 "FRIONLINE"=c(shortlabel="Friends\nOnline",label="Friend(s), mainly known online",item="X2194",wish="X2220e",section="informal"),
                 "ALTADULT"=c(shortlabel="Non-school\nAdult",label="An adult outside of school (at a sport club, another parent, family friend)",item="X2195",wish="X2220n",section="informal"),
                 ##
                 "SCHMH"=c(shortlabel="School\nMH",label="School nurse/counsellor/other pastoral staff at school",item="X2200",wish="X2220h",section="semiformal"),
                 "EMHP"=c(shortlabel="EMHP",label="Educational Mental Health Practitioner (EMHP)",item="X2201",wish="X2220i",section="semiformal"),
                 "SCHADULT"=c(shortlabel="School\nAdult",label="Another adult at school",item="X2202",wish=NA_character_,section="semiformal"),
                 "PEER"=c(shortlabel="Peer\nmentor",label="A peer mentor at school",item="X2203",wish="X2220j",section="semiformal"),
                 "OTHERSCH"=c(shortlabel="Other\nSch (F)",label="Other school services (please specify)",item="X2204",wish=NA_character_,section="freetext"),
                 ##
                 "GP"=c(shortlabel="GP",label="GP (family doctor)",item="X2210",wish="X2220f",section="formal"),
                 "SOCIAL"=c(shortlabel="Social\nWorker",label="Social worker",item="X2211",wish="X2220g",section="formal"),
                 "CAMHS"=c(shortlabel="CAMHS",label="CAMHS (NHS Child and Adolescent Mental Health Services)",item="X2212",wish="X2220k",section="formal"),
                 "THERAPIST"=c(shortlabel="Therapist",label="Private counsellor/therapist",item="X2213",wish="X2220l",section="formal"),

                 "CHARITY"=c(shortlabel="Charity",label="Support service given by a charity",item="X2214",wish="X2220m",section="semiformal"),
                 "HELPLINE"=c(shortlabel="Helpline",label="A telephone/text help-line",item="X2215",wish="X2220o",section="semiformal"),
                 "ONLINE"=c(shortlabel="Website",label="Website or online forum",item="X2216",wish="X2220p",section="semiformal"),
                 "ANON"=c(shortlabel="Anonymous\nonline",label="From an anonymous user on an online platform/chatroom/forum/server",item="X2217",wish="X2220q",section="semiformal"),
                 "OTHER"=c(shortlabel="Other\n(F)",label="Other services (please specify)",item="X2218",wish="X2220r",section="freetext")
                 )

    MAP <- lapply( MAP1, function(Z){
        c(Z,
          LABEL1=attr(SUBSET[[Z["item"]]],"meta")$text,
          LABEL2=if(is.na(Z["wish"])){NA_character_}else{attr(SUBSET[[Z["wish"]]],"meta")$text})
    })
    rm(MAP1)


    ## Mapping of script-terms to question (sub-)headers
    MAP.DT <- rbindlist(lapply(MAP,as.list))
    MAP.DT[,IDX:=1:.N]
    MAP.DT[,LABEL:=names(MAP)]

    if( 1 ) {   
        write.table(MAP.DT,row.names=FALSE,quote=FALSE,sep="|",file=file.path(DIRS$OUTPUTS,DIRTAG,sprintf("%s.label-mapping.csv",FILETAG)))
    }

    STATE3.TRANSLATE <- rbindlist( list(
        data.table(old=c("Currently being offered support"),new=c("Currently being offered support")),
        data.table(old=c("Previously been offered support"),new=c("Previously been offered support")),
        data.table(old=c("Not been offered support/been turned away", "Changed mind before getting the support", "NoResponse", "Skipped", "Incomplete(likely)", "Incomplete(uncertain)"),new=c("Other/NoResponse"))
    ))
    STATE2.TRANSLATE <- rbindlist( list(
        data.table(old=c("Currently being offered support","Previously been offered support"),new=c("Currently/Previously")),
        data.table(old=c("Not been offered support/been turned away", "Changed mind before getting the support", "NoResponse", "Skipped", "Incomplete(likely)", "Incomplete(uncertain)"),new=c("Other/NoResponse"))
        ))        
    HELP3.TRANSLATE <- rbindlist( list(
        data.table(old=c("Quite helpful", "Just about helpful enough", "Very helpful"),new=c("Helpful")),
        data.table(old=c("Not helpful at all", "Not helpful enough"),new=c("Not helpful")),
        data.table(old=c("NoResponse"),new=c("NoResponse"))
    ))


    for( lNAME in names(MAP) ) {
        cat( "=== ", lNAME, "===\n" )
        lITEM <- MAP[[lNAME]]["item"]

        SUBSET[ , (sprintf("SUP.%s.IND.all",lNAME)):=fifelse(get(sprintf("%s",lITEM))=="Ticked (Yes)",yes="YES",no="NO")]

        STATE2.TRANSLATE[,(sprintf("%sa",lITEM)):=old]
        SUBSET[ STATE2.TRANSLATE, on=sprintf("%sa",lITEM), (sprintf("SUP.%s.STATE2",lNAME)):=i.new ]
        SUBSET[ get(sprintf("SUP.%s.IND.all",lNAME))=="NO", (sprintf("SUP.%s.STATE2",lNAME)):="Other/NoResponse" ]        
        STATE2.TRANSLATE[,(sprintf("%sa",lITEM)):=NULL]

        STATE3.TRANSLATE[,(sprintf("%sa",lITEM)):=old]
        SUBSET[ STATE3.TRANSLATE, on=sprintf("%sa",lITEM), (sprintf("SUP.%s.STATE3",lNAME)):=i.new ]
        SUBSET[ get(sprintf("SUP.%s.IND.all",lNAME))=="NO", (sprintf("SUP.%s.STATE3",lNAME)):="Other/NoResponse" ]        
        STATE3.TRANSLATE[,(sprintf("%sa",lITEM)):=NULL]
        
        HELP3.TRANSLATE[,(sprintf("%sb",lITEM)):=old]
        SUBSET[ HELP3.TRANSLATE, on=sprintf("%sb",lITEM), (sprintf("SUP.%s.HELP3",lNAME)):=i.new ]
        SUBSET[ get(sprintf("SUP.%s.IND.all",lNAME))=="NO", (sprintf("SUP.%s.HELP3",lNAME)):="NoResponse" ]        
        HELP3.TRANSLATE[,(sprintf("%sb",lITEM)):=NULL]

        SUBSET[ , (sprintf("SUP.%s.IND.ever",lNAME)):=fifelse(get(sprintf("SUP.%s.IND.all",lNAME))=="YES" & get(sprintf("SUP.%s.STATE2",lNAME))=="Currently/Previously",
                                                              yes="YES",no="NO")]
        
    }


    ## Create subset of types excluding the freetext (noF)
    ##
    INDset <- list("oth"=MAP.DT[,LABEL])
    INDset[["noF"]] <- INDset[["oth"]][ !INDset[["oth"]] %in% c("OTHER","OTHERSCH") ]

    for( lSET in c("noF") ) { ## exclude freetext (noF)
        for( lSTATE in c("all","ever") ) { ## exclude responses that are not past year current/present, e.g. "did not attend" (ever)
            SUBSET[,(sprintf("SUP.TOTAL.IND.%s.%s",lSTATE,lSET)):=Reduce(`+`,lapply(mget(sprintf(sprintf("SUP.%%s.IND.%s",lSTATE),INDset[[lSET]])),function(Z){Z=="YES"})) ]
            SUBSET[,(sprintf("SUP.BANDS.IND.%s.%s",lSTATE,lSET)):=cut(get(sprintf("SUP.TOTAL.IND.%s.%s",lSTATE,lSET)),c(-0.5,0.5,1.5,2.5,Inf),c("0","1","2",sprintf("3-%i",length(INDset[[lSET]])))) ]
        }
    }

    
    AGG.SUP <- list()
    for( lSET in c("noF") ) {
        for( lSTATE in c("ever") ) {    
            AGG.SUP[[sprintf("SUP.BANDS.IND.%s.%s",lSTATE,lSET)]] <- dcast(data=SUBSET[!is.na(Sex),.N,by=c("Sex","YEARGROUP","MH.PROBLEM",sprintf("SUP.BANDS.IND.%s.%s",lSTATE,lSET))],
                                                                       formula=as.formula(sprintf("Sex + YEARGROUP + MH.PROBLEM ~ SUP.BANDS.IND.%s.%s",lSTATE,lSET)),
                                                                       fill=0, value.var="N" )
            AGG.SUP[[sprintf("SUP.BANDS.IND.%s.%s",lSTATE,lSET)]][MH.PROBLEM!="Unknown",ORDER:=order(MH.PROBLEM,Sex,YEARGROUP)]
        }
    }

    MELT.ID.VARS <- c("XID","Sex","Sex2",
                      "MH.PROBLEM","YEARGROUP","Yeargroup","YeargroupZ",
                      "RCADS2",
                      "SUP.TOTAL.IND.ever.noF")
    MELT.MEASURES <- list(
        SUP.IND.ever=sprintf("SUP.%s.IND.ever",MAP.DT$LABEL),
        SUP.STATE2=sprintf("SUP.%s.STATE2",MAP.DT$LABEL),
        SUP.STATE3=sprintf("SUP.%s.STATE3",MAP.DT$LABEL),
        SUP.HELP3=sprintf("SUP.%s.HELP3",MAP.DT$LABEL)
        )

    
    FULL.DT <- melt( SUBSET, id.vars=MELT.ID.VARS, variable.name="Support.Type", measure.vars=MELT.MEASURES )[order(XID)]

    FULL.DT[,Support:=MAP.DT$LABEL[as.numeric(Support.Type)]]
    FULL.DT[,Section:=MAP.DT$section[as.numeric(Support.Type)]]
   

    LINKS.DT <- list()
    for( lSET in c("noF") ) {
        for( lSTATE in c("ever") ) {
            cat(lSTATE,lSET,"\n")
            EXPR <- parse(text=sprintf("SUP.TOTAL.IND.%s.%s>1 & SUP.IND.%s=='YES'",lSTATE,lSET,lSTATE))
            LINKS.DT[[sprintf("%s.%s",lSTATE,lSET)]] <- FULL.DT[Support %in% INDset[[lSET]],][eval(EXPR),][, rbindlist(combn(Support,2,simplify=FALSE,FUN=function(Z){as.data.table(as.list(Z))})),by=.(XID,Sex,Sex2,LTCYN,EXCYN,SENYN,NeuroYN,NeuroNX,YEARGROUP,RCADS2,MH.PROBLEM)]
        }
    }
    
    
}


if( 1 ) {
    ##
    ## Summary tables of counts
    ##
    TABLE.HEADER <- c("0","1",sprintf("2-%i",length(INDset$noF)))
    SUM <- list()
    SUM$Y1 <- addmargins(FULL.DT[][,.N,by=.(XID,Sex,RCADS2,SUP.TOTAL.IND.ever.noF)][,xtabs( ~ 1 + cut(SUP.TOTAL.IND.ever.noF,breaks=c(-0.5,0.5,1.5,Inf),labels=TABLE.HEADER))],1)
    SUM$Y2 <- addmargins(FULL.DT[][,.N,by=.(XID,Sex,RCADS2,SUP.TOTAL.IND.ever.noF)][,xtabs( ~ Sex + cut(SUP.TOTAL.IND.ever.noF,breaks=c(-0.5,0.5,1.5,Inf),labels=TABLE.HEADER))],2)
    SUM$Y3 <- addmargins(FULL.DT[][,.N,by=.(XID,Sex,RCADS2,SUP.TOTAL.IND.ever.noF)][
       ,xtabs( ~ factor(interaction(RCADS2,Sex),
                        levels=c("Normal.Boy", "Clinical.Boy", "Normal.GD/GND", "Clinical.GD/GND","Normal.Girl", "Clinical.Girl"),
                        labels=c("Normal","ClinicalBoy","Normal","ClinicalGD/GND","Normal","ClinicalGirl")) + cut(SUP.TOTAL.IND.ever.noF,breaks=c(-0.5,0.5,1.5,Inf),labels=TABLE.HEADER))],2)
    SUM$Y4 <- addmargins(FULL.DT[][,.N,by=.(XID,MH.PROBLEM,SUP.TOTAL.IND.ever.noF)][,xtabs( ~ MH.PROBLEM + cut(SUP.TOTAL.IND.ever.noF,breaks=c(-0.5,0.5,1.5,Inf),labels=TABLE.HEADER))],2)
    SUM$Y5 <- addmargins(FULL.DT[][,.N,by=.(XID,YEARGROUP,Sex,RCADS2,SUP.TOTAL.IND.ever.noF)][,xtabs( ~ YEARGROUP + cut(SUP.TOTAL.IND.ever.noF,breaks=c(-0.5,0.5,1.5,Inf),labels=TABLE.HEADER))],2)

    
    
    PROP <- list()
    for( lLAB in names(SUM) ) {
        if( length(dim(SUM[[lLAB]]))==1 ) {
            B <- c("All",
                   format(SUM[[lLAB]],big.mark=",",trim=TRUE))            
            A <- c("",sprintf("[%4.1f%%]",round(100*SUM[[lLAB]][-4]/SUM[[lLAB]][4],1)),"-")
            PROP[[lLAB]] <- rbind(B,A)
        } else {
            B <- cbind(row.names(SUM[[lLAB]]),
                       format(SUM[[lLAB]],big.mark=",",trim=TRUE))            
            A <- cbind(rep("",NROW(SUM[[lLAB]])),
                       t(apply(round(100*prop.table(SUM[[lLAB]][,-4],1),1),1,FUN=sprintf,fmt="[%4.1f%%]")),
                       sprintf("(%4.1f%%)",round(100*prop.table(SUM[[lLAB]][,4,drop=FALSE],2),1)))
            PROP[[lLAB]] <- rbind(B,A)[order(c(seq_along(1:NROW(A)),seq_along(1:NROW(B)))),]
        }
    }
    ## do.call(rbind,PROP)
    write.table(do.call(rbind,PROP),
                row.names=FALSE,quote=FALSE,sep="|",file=file.path(DIRS$OUTPUTS,DIRTAG,sprintf("%s.summary-counts.csv",FILETAG)))


    ##
    ## Summary by category
    ## Note: not combined with above table since counts are no longer 0-18 in each row (0-X depending on category), too confusing to merge
    ##
    BOB <- dcast.data.table(data=FULL.DT[Support %in% INDset[["noF"]],.(COUNT=sum(SUP.IND.ever=="YES")),by=.(XID,RCADS2,Section)][,BAND:=cut(COUNT,breaks=c(-0.5,0.5,Inf),labels=c("0","1-18"))][,.N,by=.(RCADS2,Section,BAND)],formula=RCADS2 + Section ~ BAND, value.var="N" )[is.na(RCADS2),RCADS2:="Unable to impute RCADS"]
    BEN <- dcast.data.table(data=FULL.DT[Support %in% INDset[["noF"]],.(COUNT=sum(SUP.IND.ever=="YES")),by=.(XID,Section)][,BAND:=cut(COUNT,breaks=c(-0.5,0.5,Inf),labels=c("0","1-18"))][,.N,by=.(Section,BAND)],formula=Section ~ BAND, value.var="N" )[,RCADS2:=""]

    BILL <- rbindlist(list(BEN,BOB),use.names=TRUE)[,Total:=`0`+`1-18`][,.(RCADS2,Section,`0`,`1-18`,Total)]
    write.table(BILL,row.names=FALSE,quote=FALSE,sep="|",
                file=file.path(DIRS$OUTPUTS,DIRTAG,sprintf("%s.summary-counts-by-category.csv",FILETAG)))

    
}


if( 1 ) {
    ##
    ## Stacked histogram plots
    ##

    COLOURS.STACK <- gray.colors(4)
    TITLES <- c("ever"="only for those currently or previously offered support")    

    for( lSET in c("noF") ) {
        for( lSTATE in c("ever") ) {

            lAGG <- sprintf("SUP.BANDS.IND.%s.%s",lSTATE,lSET)
            
            png(file=file.path(DIRS$OUTPUTS,DIRTAG,sprintf("%s.SFigure1.%s.%s-barplots.png",FILETAG,lSTATE,lSET)),width=1.5*PNG.BASE$W,height=2*PNG.BASE$H,res=PNG.BASE$R)
            layout(matrix(1:2,nrow=2))
            par(oma=c(8,0,4,3.5))
            par(mar=c(0.5,4,0,0))
            GROUPINGS <- c("0","1","2",sprintf("3-%i",length(INDset[[lSET]])))
            BLOCK <- t(as.matrix(AGG.SUP[[lAGG]][!is.na(ORDER)][(ORDER),.SD,.SDcols=GROUPINGS]))
            BP <- barplot(height=BLOCK,col=COLOURS.STACK,las=2)
            text(x=BP,y=colSums(BLOCK),labels=format(colSums(BLOCK),big.mark=",",trim=TRUE),las=2,cex=0.45,font=2,xpd=NA,pos=3)
            mtext(text="Frequency",side=2,line=2.75)
            par(mar=c(0,4,0.5,0))
            
            BLOCK <- t(as.matrix(AGG.SUP[[lAGG]][!is.na(ORDER)][(ORDER),.SD/rowSums(.SD),.SDcols=c("0","1","2",sprintf("3-%i",length(INDset[[lSET]])))]))
            BP <- barplot(height=BLOCK,col=COLOURS.STACK,las=2)
            axis(side=1,at=BP,labels=AGG.SUP[[lAGG]][!is.na(ORDER)][(ORDER),sprintf("%s (%s) %s",Sex,YEARGROUP,MH.PROBLEM)],las=2,cex.axis=0.8)
            legend(x=grconvertX(x=1,from="ndc","user")-strwidth("x"),y=grconvertY(0.5,"ndc","user"),
                   legend=GROUPINGS,fill=COLOURS.STACK,bg="white",title="Support count",xpd=NA,xjust=1,yjust=0.5,cex=0.8)
            mtext(text="Proportion",side=2,line=2.75)
            text(labels="Accessing support for a mental health problem",xpd=NA,x=grconvertX(0.5,"ndc","user"),y=grconvertY(1,"ndc","user")-strheight("X",cex=1.5,font=2),cex=1.5,font=2)
            text(labels=sprintf("by gender, (yeargroup), and ever felt had a mental health problem in the past year %s",TITLES[lSTATE]),
                 xpd=NA,x=grconvertX(0.5,"ndc","user"),y=grconvertY(1,"ndc","user")-2.5*strheight("X",cex=1.5,font=2),cex=0.8,font=2)
            dev.off()
        }
    }

}



if( 1 ) {
    ##
    ## Hurdle modelling and plot
    ##

    SupportCountModel <- list()

    TITLES <- c("ever"="currently or previously offered support")

    FORMULAE <- list("A"=" SUP.TOTAL.IND.%s.noF ~ 1 + Sex + YeargroupZ + RCADS2 + MH.PROBLEM | 1 + Sex + YeargroupZ + RCADS2 + MH.PROBLEM",
                     "B"=" SUP.TOTAL.IND.%s.noF ~ 1 + Sex * YeargroupZ + RCADS2 * MH.PROBLEM | 1 + Sex * YeargroupZ + RCADS2 * MH.PROBLEM",
                     "C"=" SUP.TOTAL.IND.%s.noF ~ 1 + Sex * YeargroupZ + RCADS2 * MH.PROBLEM | 1 + Sex + MH.PROBLEM + RCADS2 + MH.PROBLEM",
                     "D"=" SUP.TOTAL.IND.%s.noF ~ 1 + Sex + YeargroupZ + RCADS2 + MH.PROBLEM | 1 + Sex * YeargroupZ + RCADS2 * MH.PROBLEM",
                     "E"=" SUP.TOTAL.IND.%s.noF ~ 1 + Sex * YeargroupZ + RCADS2 + MH.PROBLEM | 1 + Sex * YeargroupZ + RCADS2 + MH.PROBLEM",
                     "F"=" SUP.TOTAL.IND.%s.noF ~ 1 + Sex + YeargroupZ + RCADS2 * MH.PROBLEM | 1 + Sex + YeargroupZ + RCADS2 * MH.PROBLEM"
                     )

    for( lSUBSET in "ever" ) {
        FORMULA.SET <- names(FORMULAE)
        for( lFORMULA in FORMULA.SET ) {
            cat(lSUBSET,lFORMULA,"\n")
            lLABEL <- sprintf("%s-h%s",lSUBSET,lFORMULA)
            SupportCountModel[[lLABEL]] <-FULL.DT[MH.PROBLEM!="Unknown",.N,by=c("XID","Sex","YeargroupZ","RCADS2","MH.PROBLEM",sprintf("SUP.TOTAL.IND.%s.noF",lSUBSET))][,{
                hurdle(formula=as.formula(sprintf(FORMULAE[[lFORMULA]],lSUBSET)), data=.SD, dist="poisson", zero.dist="binomial" )
            }]
        }
    }
    

    ## AIC comparison
    ##
    aictab( SupportCountModel[grep("ever-",names(SupportCountModel))] )
    ##
    ## ==> best is model-B
    
    SJTAB.SIMPLE <- tab_model(SupportCountModel[grep("-hA",names(SupportCountModel))], show.aic=FALSE, transform=NULL, dv.labels=TITLES["ever"],
                              file=file.path(DIRS$OUTPUTS,DIRTAG,sprintf("%s.support-count-hurdle-simple.html",FILETAG)) )
    SJTAB.BEST <- tab_model( SupportCountModel[grep("-hB",names(SupportCountModel))], show.aic=FALSE, transform=NULL, dv.labels=TITLES["ever"],
                            file=file.path(DIRS$OUTPUTS,DIRTAG,sprintf("%s.support-count-hurdle-best.html",FILETAG)) )
    SJTAB.SIMPLE
    SJTAB.BEST



    if( 1 ) {
        
        BEST.MODEL <- "ever-hB"
        
        PRED <- CJ(Yeargroup=7:13,Sex=unique(FULL.DT$Sex),RCADS2=unique(FULL.DT$RCADS2),MH.PROBLEM=unique(FULL.DT$MH.PROBLEM)[c(1,2)])
        PRED[,YeargroupZ:=Yeargroup-7]
        PRED[,OUTever:=predict(object=SupportCountModel[[BEST.MODEL]], newdata=PRED, type="zero" )]
        PRED[,COUNTever:=predict(object=SupportCountModel[[BEST.MODEL]], newdata=PRED, type="count" )]
        PRED[,EXPECTEDever:=COUNTever/(1-exp(-1*COUNTever))]
        PRED[,RESPONSEever:=predict(object=SupportCountModel[[BEST.MODEL]], newdata=PRED, type="response" )]
        
        LTY.MH <- c("No"=1,"Yes"=2)
        COL.SEX <- c("Boy"="#2297E6","Girl"="#DF536B","GD/GND"="#000000")
        LWD.RCADS <- c("Normal"=1,"Clinical"=3)


        png(file=file.path(DIRS$OUTPUTS,DIRTAG,sprintf("%s.support-count-model-graphs.png",FILETAG)),width=1.5*PNG.BASE$W,height=PNG.BASE$H,res=PNG.BASE$R)
        layout(matrix(1:2,byrow=TRUE,nrow=1))
        par(oma=c(0,0,4,15),mar=c(4,4,0.5,0.5))
        plot(0,type="n",xlab="Year group",ylab="",xlim=c(7,13)-7,ylim=c(0,1),xaxt="n",las=1,yaxs="i")

        mtext(text="(a) Predicted probablity of accessing any support\n\u00A0",side=3,line=0.5)
        axis(side=1,at=(7:13)-7,labels=7:13)
        PRED[order(Yeargroup),{lines(OUTever~YeargroupZ,lty=LTY.MH[.BY[["MH.PROBLEM"]]],col=COL.SEX[.BY[["Sex"]]],lwd=LWD.RCADS[.BY[["RCADS2"]]])},
             by=.(Sex,RCADS2,MH.PROBLEM)]

        plot(0,type="n",xlab="Year group",ylab="",xlim=c(7,13)-7,ylim=c(1,3),xaxt="n",las=1,yaxs="i")

        mtext(text="(b) Expected conditional number of different\ntypes of support accessed",side=3,line=0.5)
        axis(side=1,at=(7:13)-7,labels=7:13)

        PRED[!is.na(Sex)&!is.na(RCADS2)&!is.na(MH.PROBLEM)][order(Yeargroup),{lines(COUNTever~YeargroupZ,lty=LTY.MH[.BY[["MH.PROBLEM"]]],col=COL.SEX[.BY[["Sex"]]],lwd=LWD.RCADS[.BY[["RCADS2"]]])},by=.(Sex,RCADS2,MH.PROBLEM)]

        text(x=grconvertX(0.5,"ndc","user"),y=grconvertY(1,"ndc","user")-1*par("cxy")[2],
             labels=sprintf("Hurdle model of the number of different types of support accessed (%s)",TITLES["ever"]),
             font=2,cex=1.25,xpd=NA)


        LEGEND <- expand.grid(RCADS=names(LWD.RCADS),MH=names(LTY.MH),SEX=names(COL.SEX),stringsAsFactors=FALSE)
        LEGEND$label <- apply(LEGEND,1,function(Z){paste(rev(Z),collapse=", ")})

        BILL <- na.omit(FULL.DT[,.N,by=.(XID,Sex,YeargroupZ,RCADS2,MH.PROBLEM)][,.N,by=.(Sex,MH.PROBLEM,RCADS2)])[ LEGEND, on=.(RCADS2=RCADS,MH.PROBLEM=MH,Sex=SEX) ][,sum(N)]
        BOB <- na.omit(FULL.DT[,.N,by=.(XID,Sex,YeargroupZ,RCADS2,MH.PROBLEM,SUP.TOTAL.IND.ever.noF)][,.N,by=.(Sex,MH.PROBLEM,RCADS2)])[ LEGEND, on=.(RCADS2=RCADS,MH.PROBLEM=MH,Sex=SEX) ][,sprintf("%s (%s)",label,format(N,trim=TRUE,big.mark=","))]        
        
        legend(x=mean(c(grconvertX(1,"ndc"),grconvertX(1,"nfc"))),y=grconvertY(0.5,"ndc"),seg.len=3,bg="white",title=sprintf("Gender, MH problem, RCADS11 (%s)",format(BILL,big.mark=",")),
               legend=BOB,col=COL.SEX[LEGEND$SEX],lty=LTY.MH[LEGEND$MH],lwd=LWD.RCADS[LEGEND$RCADS],xpd=NA,xjust=0.5,yjust=0.5,cex=0.8)

        dev.off()
        

    }
}
    

if( 1 ) {
    ##
    ## Network plots (via igraph)
    ##
    
    VERTEX.BREAKS <- c(0,0.4,0.5,0.6,0.7,0.8,1) 
    VERTEX.COLOURS <- c("#EEEEEE",brewer.pal(5,"YlGn"))
    VERTEX.SELECT.COLOURS <- c(brewer.pal(4,"Blues")[c(3,4)])
    
    Add.Gradient.Legend <- function(Title="",Breaks=c(0,1),Colours="blue",Position=c(0.475,0.525,5/8,7/8),Cex=0.8,Mar=if(horizontal){c(0,0,3,0)}else{c(0,3,0,0)},horizontal=FALSE) {
        subplot.par <- par(fig=Position,new=TRUE,mar=Mar)
        if( horizontal ) {
            plot(0,type="n",ylim=c(0,1),xlim=range(Breaks),axes=FALSE,xlab="",ylab="")
            rect(ybottom=0,ytop=1,xleft=Breaks[-length(Breaks)],xright=Breaks[-1],col=Colours)
            mtext(text=Breaks,side=3,at=Breaks,las=1,cex=Cex)
            mtext(text=Title,side=3,at=grconvertX(0.5,from="npc",to="user"),line=1,xpd=NA,cex=1.1*Cex)            
        } else {
            plot(0,type="n",xlim=c(0,1),ylim=range(Breaks),axes=FALSE,xlab="",ylab="")
            rect(xleft=0,xright=1,ybottom=Breaks[-length(Breaks)],ytop=Breaks[-1],col=Colours)
            mtext(text=Breaks,side=2,at=Breaks,las=2,cex=Cex)
            mtext(text=Title,side=3,at=0.5,line=1,xpd=NA,cex=1.1*Cex)
        }
    }

    mystar <- function(coords, v = NULL, params) {
        vertex.color <- params("vertex", "color")
        if (length(vertex.color) != 1 && !is.null(v)) {
            vertex.color <- vertex.color[v]
        }
        vertex.frame.color <- params("vertex", "frame.color")
        if (length(vertex.frame.color) != 1 && !is.null(v)) {
            vertex.frame.color <- vertex.frame.color[v]
        }        
        vertex.size <- 1 / 200 * params("vertex", "size")
        if (length(vertex.size) != 1 && !is.null(v)) {
            vertex.size <- vertex.size[v]
        }
        norays <- params("vertex", "norays")
        if (length(norays) != 1 && !is.null(v)) {
            norays <- norays[v]
        }
        
        mapply(coords[, 1], coords[, 2], vertex.color, vertex.frame.color, vertex.size, norays,
               FUN = function(x, y, bg, fg, size, nor) {
                   symbols(
                       x = x, y = y, bg = bg, fg=fg,
                       stars = matrix(c(size, (2*size / sqrt(3))), nrow = 1, ncol = nor * 2),
                       add = TRUE, inches = FALSE
                   )
               }
               )
    }

    add_shape("star",clip = shape_noclip, plot = mystar, parameters = list(vertex.norays = 6) ) 

    SHAPE.MAP <- data.table(section=c("informal", "semiformal", "formal"),
                            shape=c("circle","square","star"),
                            lty=1,
                            lambda=c(1,sqrt(pi),sqrt((2*pi)/3*sqrt(3))))
    MAP.DT[SHAPE.MAP,on=.(section),c("shape","lambda"):=.(i.shape,i.lambda)] ## NEED TO FIX RELATIVE AREA OF SHAPES
    MAP.DT[,lty:=fifelse(LABEL%in%c("OTHER","OTHERSCH"),2,1)]
    ## hexagon size a = 3*sqrt(3)*a^2/2
    ## circle size a = pi*a^2
    ## square size a = a^2
    
    MakeGraphElement <- function(Title, State="ever", Set="noF", Select="TRUE", EdgeDT, FullDT, Highlight=NULL, Lowlight=NULL, VertexMap=MAP.DT,
                                 vertex.selector=sprintf("SUP.IND.%s =='YES'",State),
                                 verbose=TRUE ) {
        if( verbose ) cat(Title,"\n")

        EdgeDT <- EdgeDT[eval(parse(text=Select))]
        FullDT <- FullDT[Support %in% INDset[[Set]],][eval(parse(text=Select))]
        VertexDT <- FullDT[Support %in% INDset[[Set]],][eval(parse(text=vertex.selector))]
        
        OUT <- list()
        OUT$FULLDT <- FullDT ## useful for making grids later in script, but very memory inefficient
        
        ## NOTE: edge.style must be for global and anchored separately
        if( 1 ) {
            if( !is.null(Highlight) | !is.null(Lowlight) ) {
                UPPER <- c(0.0100,0.0500,
                           0.1000,0.1250,0.200,
                           0.3000,0.6500,1.0000)
                BANDS.STYLE <- "anchored"
            } else {
                UPPER <- c(0.0050,0.0100,
                           0.0250,0.0500,0.0750,
                           0.1000,0.2000,1.0000)
                BANDS.STYLE <- "standard"

            }
            EDGE.STYLE.DT <- data.table(Upper=UPPER,Bands=BANDS.STYLE)
            EDGE.STYLE.DT[,color:=c("darkgray","darkgrey","darkgrey","black","black","black","black","black")]
            EDGE.STYLE.DT[,width:=c(0.00,0.75,1.50,1.50,3.00,5.00,7.00,0.00)]            
        } else {
            ## UPPER names are approximate quantiles across all edge weights (above are global/anchored specific)
            UPPER <- c("0100"=0.0010,"0250"=0.0050,"0500"=0.0100,"0750"=0.0250,"0900"=0.0500,"0950"=0.1000,"0975"=0.1500,"0999"=0.6500,"1000"=1.0000)
            EDGE.STYLE.DT <- data.table(Upper=UPPER,Bands="global")
            EDGE.STYLE.DT[,color:=c("darkgray","darkgray","darkgrey","darkgrey","black","black","black","black","black")]
            EDGE.STYLE.DT[,width:=c(0.00,0.25,0.50,1.00,1.00,3.00,5.00,7.00,0.00)]
        }
        EDGE.STYLE.DT[,Lower:=c(0,Upper[-1*(.N)])]
        EDGE.STYLE.DT[,BIN:=1:(.N)]
        
        VERTEX.MAX <- VertexDT[,.N,by=.(XID)][,.N]
        OUT$EDGEDT <- EdgeDT[,.(Title=Title,weight=.N,Proportion=.N/VERTEX.MAX),by=.(V1,V2)][,BIN:=findInterval(x=Proportion,vec=c(0,EDGE.STYLE.DT$Upper),all.inside=TRUE)][EDGE.STYLE.DT,on=.(BIN),nomatch=0]

        if( any( OUT$EDGEDT$BIN == 8 ) ) {
            cat(Title,"\n")
            stop("Edge intervals are HARD CODED, choice was made based on networks for paper, they may not be valid in all cases!")
        }
        
        OUT$EDGE.STYLE <- copy(EDGE.STYLE.DT)
        OUT$VERTEXDT <- data.table(V=INDset[[Set]],
                                   Title=Title,
                                   label.cex=0.90,
                                   label.color="black",
                                   Vprop.t=VertexDT[SUP.HELP3!="NoResponse",.(Prop=sum(SUP.HELP3=="Helpful")),by=.(Support)][match(INDset[[Set]],Support),Prop],
                                   Vprop.b=VertexDT[SUP.HELP3!="NoResponse",.(Prop=.N),by=.(Support)][match(INDset[[Set]],Support),Prop],
                                   Vcount=VertexDT[,.(N=.N),by=.(Support)][match(INDset[[Set]],Support),N],
                                   label=VertexMap[INDset[[Set]],on=.(LABEL),shortlabel],
                                   shape=VertexMap[INDset[[Set]],on=.(LABEL),shape]
                                   )

        OUT$VERTEXDT[,size:=10 + (10*sqrt(Vcount)/25)] ## function to vary node size

        OUT$VERTEXDT[,Vprop:=Vprop.t/Vprop.b]
        
        OUT$VERTEXDT[,Vprop5:=cut(x=Vprop,breaks=VERTEX.BREAKS,include.lowest=TRUE)]
        OUT$VERTEXDT[,color:=VERTEX.COLOURS[Vprop5]]
        if( !is.null(Highlight) ) {
            MATCH <- OUT$VERTEXDT[,match(Highlight,V)]
            if( any(is.na(MATCH)) ) warning("Some Highlight vertices do not match")
            OUT$VERTEXDT[na.omit(MATCH),color:=VERTEX.SELECT.COLOURS[2]]
        }
        if( !is.null(Lowlight) ) {
            MATCH <- OUT$VERTEXDT[,match(Lowlight,V)]
            if( any(is.na(MATCH)) ) warning("Some Lowlight vertices do not match")
            OUT$VERTEXDT[na.omit(MATCH),color:=VERTEX.SELECT.COLOURS[1]]
        }

        OUT$TITLE <- Title
        OUT$SUBTITLE1 <- sprintf("(n=%i, v=%i, e=%i)",
                                 FullDT[,.N,by=XID][,.N],
                                 VertexDT[,.N,by=XID][,.N],
                                 EdgeDT[,.N,by=XID][,.N]
                                 )
        OUT$SUBTITLE2 <- sprintf("(v/n=%2.1f%%, e/v=%2.1f%%)",
                                 round(100*VertexDT[,.N,by=XID][,.N]/FullDT[,.N,by=XID][,.N],1),
                                 round(100*EdgeDT[,.N,by=XID][,.N]/VertexDT[,.N,by=XID][,.N],1)
                                 )
        return(OUT)
    }

  
    
    GRAPH.LIST <- list()    
    for( lSET in c("noF") ) {
        for( lSTATE in c("ever") ) {
            LABEL <- sprintf("%s.%s",lSTATE,lSET)
            cat(LABEL,"\n")
            
            GRAPH.LIST[[LABEL]] <- list()
            GRAPH.LIST[[LABEL]][["ALL"]] <- MakeGraphElement(Title="All",State=lSTATE,Set=lSET,
                                                             EdgeDT=LINKS.DT[[LABEL]],
                                                             FullDT=FULL.DT
                                                             )

            for( lSEX in c("Girl","Boy","GD/GND") ) {
                GRAPH.LIST[[LABEL]][[sprintf("%s",lSEX)]] <- MakeGraphElement(Title=lSEX,State=lSTATE,Select=sprintf("Sex=='%s'",lSEX),
                                                                              EdgeDT=LINKS.DT[[LABEL]],
                                                                              FullDT=FULL.DT
                                                                              )
            }
            
            GRAPH.LIST[[LABEL]][["Rnorm"]] <- MakeGraphElement(Title="RCADS(Normal)",State=lSTATE,Select="RCADS2=='Normal'",
                                                               EdgeDT=LINKS.DT[[LABEL]],
                                                               FullDT=FULL.DT
                                                               )
            GRAPH.LIST[[LABEL]][["Rclinic"]] <- MakeGraphElement(Title="RCADS(Clinical)",State=lSTATE,Select="RCADS2=='Clinical'",
                                                                 EdgeDT=LINKS.DT[[LABEL]],
                                                                 FullDT=FULL.DT
                                                                 )
            

            for( lCOLUMN in c("CARER","SCHMH","CAMHS") ) {
                for( lHELP in c("Helpful","Not helpful") ) {
                    SELECT.DT <- FULL.DT[ Support==lCOLUMN & SUP.HELP3==lHELP & SUP.STATE2=="Currently/Previously", .(XID)]
                    GRAPH.LIST[[LABEL]][[sprintf("%s-%s",lCOLUMN,lHELP)]] <- MakeGraphElement(
                        EdgeDT=LINKS.DT[[LABEL]][ SELECT.DT, , on=.(XID), nomatch=NULL ],State=lSTATE,
                        FullDT=FULL.DT[ SELECT.DT, , on=.(XID), nomatch=NULL ],
                        Highlight=if(lHELP=="Helpful"){lCOLUMN}else{NULL},
                        Lowlight=if(lHELP=="Not helpful"){lCOLUMN}else{NULL},
                        Title=sprintf("%s anchored as %s",MAP.DT[LABEL==lCOLUMN,label],tolower(lHELP)),
                        Select="TRUE"
                    )
                }
            }
            
        }
    }    

    for( lLABEL in names(GRAPH.LIST) ) {
        for( lVAR in names(GRAPH.LIST[[lLABEL]]) ) {
            cat(lLABEL,lVAR,"\n")
            GRAPH.LIST[[lLABEL]][[lVAR]]$GRAPH <- graph_from_data_frame(d=GRAPH.LIST[[lLABEL]][[lVAR]]$EDGEDT,
                                                                        directed=FALSE,
                                                                        vertices=GRAPH.LIST[[lLABEL]][[lVAR]]$VERTEXDT)
        }
    }


    if( 0 ) {
        set.seed(1234)
        GRAPH.PLOTTING <- list(
            COORDa=norm_coords( layout_( graph=GRAPH.LIST[["ever.noF"]][["ALL"]]$GRAPH, layout=nicely() ) ),
            COORDb=norm_coords( layout_( graph=GRAPH.LIST[["ever.noF"]][["ALL"]]$GRAPH, layout=with_fr(weights=NULL) ) ),
            COORDc=norm_coords( layout_( graph=GRAPH.LIST[["ever.noF"]][["ALL"]]$GRAPH, layout=with_fr() ) ),
            COORDe=norm_coords( layout_( graph=GRAPH.LIST[["ever.noF"]][["ALL"]]$GRAPH, layout=with_fr(weights=sqrt(E(GRAPH.LIST[["ever.noF"]][["ALL"]]$GRAPH)$weight)) ) ),
            COORDf=norm_coords( layout_( graph=GRAPH.LIST[["ever.noF"]][["ALL"]]$GRAPH, layout=with_fr(weights=log(E(GRAPH.LIST[["ever.noF"]][["ALL"]]$GRAPH)$weight)) ) ),
            COORDd=norm_coords( layout_( graph=GRAPH.LIST[["ever.noF"]][["ALL"]]$GRAPH, layout=in_circle() ) ),
            XLIM=c(-1.0,1.0),
            YLIM=c(-1.1,1.1))

        ##
        ## Investigate various layouts (note use of random seed to ensure reproducibility)
        ##
        ## 
        ##
        
        saveRDS(GRAPH.PLOTTING,file=file.path(DIRS$OUTPUTS,DIRTAG,sprintf("%s.GRAPH.PLOTTING.rds",FILETAG)))
    } else {
        GRAPH.PLOTTING <- readRDS(file=file.path(DIRS$OUTPUTS,DIRTAG,sprintf("%s.GRAPH.PLOTTING.rds",FILETAG)))

    }
    
    COORD.SET <- c("COORDb") ## Best layout
    
    VAR.SETS <- list("p11-sex"=c("ALL"=1,"Girl"=2,"Boy"=3,"GD/GND"=4),
                     ##
                     "p12-rcads"=c("ALL"=1,"Rclinic"=2),
                     ##
                     "p13-carer+schmh+camhs"=c("CARER-Helpful"=1,"CARER-Not helpful"=2, "SCHMH-Helpful"=3,"SCHMH-Not helpful"=4, "CAMHS-Helpful"=5,"CAMHS-Not helpful"=6)
                     )

    LIMIT.SET <- list( list(E=0,V=0) ) ## Plot all nodes and all edges
    SINGLE <- FALSE ## multiple networks per plot
    for( lSTYLE in c("BoG","WoG") ) { ## Black-on-green (BoG) or white-on-green (WoG) for the node label text colour
        for( lLIMIT in 1:length(LIMIT.SET) ) {
            for( lCOORD in COORD.SET ) {
                for( lSUBSET in names(GRAPH.LIST) ) {
                    for( iIND in 1:length(VAR.SETS) ) {

                        E.LIMIT <- LIMIT.SET[[lLIMIT]]$E
                        V.LIMIT <- LIMIT.SET[[lLIMIT]]$V

                        cat(lCOORD,lSUBSET,names(VAR.SETS)[iIND],E.LIMIT,V.LIMIT,lSTYLE,"\n")

                        if( length(VAR.SETS[[iIND]])==2 ) {
                            GRID <- c(nrow=1,ncol=2)
                        } else if( length(VAR.SETS[[iIND]])==4 ) {
                            GRID <- c(nrow=2,ncol=2)
                        } else if( length(VAR.SETS[[iIND]])==6 ) {
                            GRID <- c(nrow=3,ncol=2)
                        } else {
                            stop("Unsupported length")
                        }
                        
                        if(!SINGLE) {
                            png(file=file.path(DIRS$OUTPUTS,DIRTAG,sprintf("%s.%s-%s-%s.graphsE%sV%s.%s.png",FILETAG,lSUBSET,names(VAR.SETS)[iIND],lCOORD,sprintf("%02i",E.LIMIT),sprintf("%02i",V.LIMIT),lSTYLE)),
                                width=GRID["ncol"]*PNG.BASE$W,height=GRID["nrow"]*PNG.BASE$H,res=PNG.BASE$R)
                        }
                        layout(matrix(1:prod(GRID),nrow=GRID["nrow"],ncol=GRID["ncol"],byrow=TRUE))
                        par(oma=c(0,0,0,0),mar=c(0,0,0,0),cex=0.66)
                        for( INDEX in sort(VAR.SETS[[iIND]]) ) { 

                            lVAR <- names( which( VAR.SETS[[iIND]] == INDEX ) )
                            E.TMP <- E(GRAPH.LIST[[lSUBSET]][[lVAR]]$GRAPH)$width
                            V.TMP <- V(GRAPH.LIST[[lSUBSET]][[lVAR]]$GRAPH)$color
                            C.TMP <- V(GRAPH.LIST[[lSUBSET]][[lVAR]]$GRAPH)$label.color

                            E(GRAPH.LIST[[lSUBSET]][[lVAR]]$GRAPH)$width <- ifelse(E(GRAPH.LIST[[lSUBSET]][[lVAR]]$GRAPH)$weight>E.LIMIT,
                                                                                   E(GRAPH.LIST[[lSUBSET]][[lVAR]]$GRAPH)$width,
                                                                                   0)
                            V(GRAPH.LIST[[lSUBSET]][[lVAR]]$GRAPH)$color <- ifelse(V(GRAPH.LIST[[lSUBSET]][[lVAR]]$GRAPH)$Vprop.b>V.LIMIT,
                                                                                   V(GRAPH.LIST[[lSUBSET]][[lVAR]]$GRAPH)$color,
                                                                                   "white")
                            V(GRAPH.LIST[[lSUBSET]][[lVAR]]$GRAPH)$frame.color <- ifelse(V(GRAPH.LIST[[lSUBSET]][[lVAR]]$GRAPH)$Vprop.b>V.LIMIT,
                                                                                         "black",
                                                                                         "white")

                            if( lSTYLE=="BoG" ) {
                                LABEL.COLOUR <- rep("black",length(V(GRAPH.LIST[[lSUBSET]][[lVAR]]$GRAPH)$color))
                            } else if( lSTYLE=="WoG" ){
                                LABEL.COLOUR <- c("#EEEEEE"="black", "#FFFFCC"="black", "#C2E699"="black", "#78C679"="black", "#31A354"="black", "#006837"="white")[V(GRAPH.LIST[[lSUBSET]][[lVAR]]$GRAPH)$color]
                            } else {
                                stop("Unknown style")
                            }
                            V(GRAPH.LIST[[lSUBSET]][[lVAR]]$GRAPH)$label.color <- LABEL.COLOUR
                            
                            PIG <- plot.igraph(x=GRAPH.LIST[[lSUBSET]][[lVAR]]$GRAPH, layout=GRAPH.PLOTTING[[lCOORD]],
                                               frame=FALSE, axes=FALSE, rescale=FALSE, xlim=GRAPH.PLOTTING$XLIM, ylim=GRAPH.PLOTTING$YLIM )

                            E(GRAPH.LIST[[lSUBSET]][[lVAR]]$GRAPH)$width <- E.TMP ## reset with saved values
                            V(GRAPH.LIST[[lSUBSET]][[lVAR]]$GRAPH)$color <- V.TMP ## reset with saved values
                            V(GRAPH.LIST[[lSUBSET]][[lVAR]]$GRAPH)$label.color <- C.TMP ## reset with saved values
                            V(GRAPH.LIST[[lSUBSET]][[lVAR]]$GRAPH)$frame.color <- "black"
                            
                            text(x=grconvertX(0,"nfc","user")+2*strwidth("X")-1.0*strwidth("X"),
                                 y=grconvertY(1,"nfc","user")-1.5*strheight("X"),
                                 labels=sprintf("(%s)",letters[INDEX]),cex=1.25,font=2,xpd=NA,adj=0)
                            text(x=grconvertX(0,"nfc","user")+2*strwidth("X")+2*strwidth("X"),
                                 y=grconvertY(1,"nfc","user")-1.5*strheight("X"),
                                 labels=sprintf("%s",GRAPH.LIST[[lSUBSET]][[lVAR]]$TITLE),cex=1.25,font=2,xpd=NA,adj=0)
                            text(x=grconvertX(0,"nfc","user")+2*strwidth("X")+2*strwidth("X"),
                                 y=grconvertY(1,"nfc","user")-2*strheight("X")-1.5*strheight("X"),
                                 labels=sprintf("%s",GRAPH.LIST[[lSUBSET]][[lVAR]]$SUBTITLE1),cex=1,font=1,xpd=NA,adj=0)
                            text(x=grconvertX(0,"nfc","user")+2*strwidth("X")+2*strwidth("X"),
                                 y=grconvertY(1,"nfc","user")-2*strheight("X")-3.25*strheight("X"),
                                 labels=sprintf("%s",GRAPH.LIST[[lSUBSET]][[lVAR]]$SUBTITLE2),cex=1,font=1,xpd=NA,adj=0)

                            text(x=grconvertX(0,"nfc","user")+2*strwidth("X")+2*strwidth("X"),
                                 y=grconvertY(1,"nfc","user")-2*strheight("X")-5.00*strheight("X"),
                                 labels=if(E.LIMIT>0){sprintf("Edges with fewer than %i individuals omitted",E.LIMIT)}else{""},cex=1,font=1,xpd=NA,adj=0)
                            text(x=grconvertX(0,"nfc","user")+2*strwidth("X")+2*strwidth("X"),
                                 y=grconvertY(1,"nfc","user")-2*strheight("X")-6.75*strheight("X"),
                                 labels=if(V.LIMIT>0){sprintf("Vertices with fewer than %i individuals omitted",V.LIMIT)}else{""},cex=1,font=1,xpd=NA,adj=0)
                        }

                        if( GRID["nrow"] > 1 ) {
                            Add.Gradient.Legend("Found Helpful\n(%)",100*VERTEX.BREAKS,VERTEX.COLOURS,
                                                Position=c((4/8)-0.04-0.015,(4/8)+0.04-0.015,1/(4*GRID["nrow"]),3/(4*GRID["nrow"])))
                        } else {
                            Add.Gradient.Legend("Found Helpful\n(%)",100*VERTEX.BREAKS,VERTEX.COLOURS,horizontal=TRUE,
                                                Position=c((4/8)-0.125-0.015,(4/8)+0.125-0.015,0,0.65/(4*GRID["nrow"])))
                        }
                        GRAPH.LIST[[lSUBSET]][[lVAR]]$EDGE.STYLE[,{
                            LB <- legend(x=grconvertX(0.5,from="ndc","user"),y=grconvertY(1-(2/((4*GRID["nrow"]))),from="ndc",to="user"),lwd=width,col=color,xpd=NA,bg="white",
                                         title="Co-occurring type (% of v)",
                                         xjust=0.5,yjust=0.5,bty="n",
                                         legend=sprintf("(%s,%s)",format(100*Lower,digits=0,nsmall=1,justify="right"),format(100*Upper,digits=0,nsmall=1,justify="right")))
                            text(x=LB$text$x,y=LB$text$y,labels=fifelse(BIN==min(BIN),"(omitted)",NA_character_),pos=2,xpd=NA)
                            text(x=LB$text$x,y=LB$text$y,labels=fifelse(BIN==max(BIN),"(no edges)",NA_character_),pos=2,xpd=NA)
                        }]                    

                        if(!SINGLE) dev.off()
                    }
                }
            }
        }
    }
}


if( 1 ) {
    ##
    ## Grid plots
    ##

    BOB <- FULL.DT[,{list(list(SUP.STATE3=unique(SUP.STATE3),SUP.HELP3=unique(SUP.HELP3)))}]
    CJ.DT <- do.call(what=CJ,args=c(INDset["noF"],BOB$V1))
    setnames(CJ.DT,new=c("Support","SUP.STATE3","SUP.HELP3"))
    SUP.STATE3.ORDER <- c("Currently being offered support","Previously been offered support","Other/NoResponse")
    SUP.HELP3.ORDER <- c("Helpful","Not helpful","NoResponse")
    PERPAGE <- 20

    for( THRESHOLD in c(10) ) { ## Censor counts <=10
        for( lGRAPH in c("ALL") ) {
            cat( THRESHOLD, lGRAPH, "\n" )
            
            COUNT.GRIDS <- GRAPH.LIST[["ever.noF"]][[lGRAPH]]$FULLDT[ SUP.IND.all=="YES",.(Count=.N), by=.(Support,SUP.STATE3,SUP.HELP3)][ CJ.DT,, on=.(Support,SUP.STATE3,SUP.HELP3) ]
            COUNT.GRIDS[is.na(Count),Count:=0]

            COUNT.GRIDS[,SUP.STATE3.ord:=match(SUP.STATE3,SUP.STATE3.ORDER)]
            COUNT.GRIDS[,SUP.HELP3.ord:=match(SUP.HELP3,SUP.HELP3.ORDER)]
            

            PAGE.MAP <- rep(1:ceiling(COUNT.GRIDS[,length(unique(Support))]/PERPAGE),each=PERPAGE)#[1:length(OFFSETS)]
            if( COUNT.GRIDS[,length(unique(Support))] < PERPAGE ) {
                PAGE.MAP[(COUNT.GRIDS[,length(unique(Support))]+1):length(PAGE.MAP)] <- NA
            }
            for( lPAGE in sort(na.omit(unique(PAGE.MAP))) ) {
                png(file=file.path(DIRS$OUTPUTS,DIRTAG,sprintf("%s.grids%i.%s-%i.png",FILETAG,THRESHOLD,lGRAPH,lPAGE)),width=2*PNG.BASE$W,height=3*PNG.BASE$H,res=PNG.BASE$R)
                SEQ <- seq(1,length(which(PAGE.MAP==lPAGE)))
                if(length(SEQ)<PERPAGE){
                    PAD <- rep(0,PERPAGE-length(SEQ))
                } else {
                    PAD <- numeric(0)
                }
                layout( matrix(c(SEQ,PAD),ncol=4,byrow=TRUE) )
                par(mar=c(10,8,3,1))
                BREAKS <- c(0,0.05,0.10,0.15,0.20,0.25,0.50,1)
                HEAT.COLOURS <- brewer.pal(length(BREAKS)-1,"Blues")
                CEX <- 1.5
                for( lIDX in sort(which(PAGE.MAP==lPAGE)) ) {
                    MAT <- COUNT.GRIDS[order(SUP.STATE3.ord,SUP.HELP3.ord)][ Support==INDset[["noF"]][lIDX] & !is.na(SUP.STATE3), matrix(Count,nrow=3,ncol=3,byrow=TRUE) ]
                    MAT.ROWS <- rowSums(MAT)
                    MAT.COLS <- colSums(MAT)
                    MAT.SUM <- sum(MAT)

                    image(z=(MAT)/sum(MAT),x=1:NROW(MAT),y=1:NCOL(MAT),xlim=c(0.5,NROW(MAT)+1.5),ylim=c(0.5,NCOL(MAT)+1.5),axes=FALSE,xlab="",ylab="",breaks=c(0,0.9,1),col=c("#FFFFFE","#FFFFFF"))
                    MAT.COLOUR <- VERTEX.COLOURS[findInterval(x=sum(MAT[1:2,1])/sum(MAT[1:2,1:2]),vec=VERTEX.BREAKS,all.inside=TRUE,left.open=TRUE)]
                    if( grepl("helpful",lGRAPH,ignore.case=TRUE) ) {
                        if( INDset[["noF"]][lIDX] == sub("(.*)-.*","\\1",lGRAPH) ) {
                            MAT.COLOUR <- VERTEX.SELECT.COLOURS[ if(sub(".*-(.*)","\\1",lGRAPH)=="Helpful"){2}else{1} ]
                        }
                    }
                    if(1) {
                        ## checking we match the network plots
                        if( GRAPH.LIST[["ever.noF"]][[lGRAPH]]$VERTEXDT[V==INDset[["noF"]][lIDX],color] != MAT.COLOUR ) {
                            cat(INDset[["noF"]][lIDX], "colour mismatch with network plot\n")
                        }
                    }
                    

                    MAT.CHECK <- any(MAT<THRESHOLD)
                    if( MAT.CHECK ) {
                        WHERE <- (MAT<THRESHOLD)
                        MAT[which(WHERE,arr.ind=TRUE)] <- NA
                        MAT.ROWS[apply(WHERE,1,any)] <- paste0((1+ceiling(MAT.ROWS[apply(WHERE,1,any)]/THRESHOLD))*THRESHOLD,"~")
                        MAT.COLS[apply(WHERE,2,any)] <- paste0((1+ceiling(MAT.COLS[apply(WHERE,2,any)]/THRESHOLD))*THRESHOLD,"~")
                        MAT.SUM <- paste0((1+ceiling(MAT.SUM/(2*THRESHOLD)))*2*THRESHOLD,"~")
                    }

                    axis(1,at=1:NROW(MAT),labels=lapply(strwrap(SUP.STATE3.ORDER,width=20,simplify=FALSE),function(Z){paste(Z,collapse="\n")}),las=2)
                    axis(2,at=1:NCOL(MAT),labels=lapply(strwrap(SUP.HELP3.ORDER,width=20,simplify=FALSE),function(Z){paste(Z,collapse="\n")}),las=2)
                    abline(v=seq(-0.5,NROW(MAT)+0.5,by=1),lty=2,xpd=FALSE)
                    abline(h=seq(-0.5,NCOL(MAT)+0.5,by=1),lty=2,xpd=FALSE)
                    abline(v=NROW(MAT)+0.5,lty=1,xpd=FALSE)
                    abline(h=NCOL(MAT)+0.5,lty=1,xpd=FALSE)            
                    rect(xleft=0.5,xright=2.5,ytop=2.5,ybottom=0.5,border=NA,col=MAT.COLOUR)
                    if( NROW(which(!is.na(MAT),arr.ind=TRUE)) > 0 ) {
                        text(x=which(!is.na(MAT),arr.ind=TRUE),labels=MAT[which(!is.na(MAT),arr.ind=TRUE)],cex=CEX)
                    }
                    if( MAT.CHECK ) {
                        text(x=which(is.na(MAT),arr.ind=TRUE),labels="~",cex=CEX)
                    }
                    text(x=1:NROW(MAT),y=NCOL(MAT)+1,labels=MAT.ROWS,font=2,cex=CEX)

                    
                    text(x=NROW(MAT)+1,y=1:NCOL(MAT),labels=MAT.COLS,font=2,cex=CEX)
                    text(x=NROW(MAT)+1,y=NCOL(MAT)+1,labels=MAT.SUM,font=2,cex=CEX)
                    text(x=grconvertX(0.5,"nfc","user"),y=grconvertY(1,"npc","user")+(3.5*strheight("X")),
                         labels=paste(strwrap(MAP.DT[LABEL==INDset[["noF"]][lIDX],label],40),collapse="\n"),cex=CEX,font=2,xpd=NA,adj=0.5)

                }
                if( any(is.na(PAGE.MAP)) & lPAGE==max(PAGE.MAP,na.rm=TRUE) ) {
                    Add.Gradient.Legend("Found Helpful\n(%)",100*VERTEX.BREAKS,VERTEX.COLOURS,Position=c(7/8-0.05,7/8+0.05,0/16,3/16))
                }
                dev.off()
            }
        }
    }
    
}
