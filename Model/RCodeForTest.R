if(a_3<0){result1<-paste0("ln(Methyl)=",a_2,"Xln(ChIP)",a_3,"; R^2=",result_r.squared,"; n=",n)}else{result1<-paste0("ln(Methyl)=",a_2,"Xln(ChIP)+",a_3,"; R^2=",result_r.squared,"; n=",n)}
p1=ggplot()+
     geom_hex(mapping=aes(x=mydata.GCH_1$Methyl,y=mydata.GCH_1$predicted),color="black")+
     labs(x = "Observed", y = "Model-predicted", title = "chr1: The Scatter Plot of mGpC(JY608-bg)",subtitle=result1)+
     theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=8))+
     expand_limits(x = 0)
p2=ggplot()+
     geom_hex(mapping=aes(x=mydata.GCH_1$Methyl,y=mydata.GCH_1$predicted),color="black",bins=60)+
     labs(x = "Observed", y = "Model-predicted", title = "chr1: The Scatter Plot of mGpC(JY608-bg)",subtitle=result1)+
     theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=8))+
     expand_limits(x = 0)
p3=ggplot()+
     geom_point(mapping=aes(x=mydata.GCH_1$Methyl,y=mydata.GCH_1$predicted),col=rgb(0, 0, 255, 40, maxColorValue=255))+
     labs(x = "Observed", y = "Model-predicted", title = "chr1: The Scatter Plot of mGpC(JY608-bg)",subtitle=result1)+
     theme(plot.title = element_text(hjust = 0.5,vjust = 1),plot.subtitle = element_text(hjust = 0.5,size=8))+
     expand_limits(x = 0)
     
ggplot2.multiplot(p1,p2,p3,p2,cols=2)		


Hexagonal heatmap of 2d bin counts


Transparency is set to 40%

It utilizes the shade of the color dots to indicate the density. 
