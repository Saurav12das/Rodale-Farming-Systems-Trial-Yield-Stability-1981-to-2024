# FST Farming Systems Trial (1981–2024)
# SPI-3 drought exposure and observed grain yields; descriptive figure supplement
# Generated from SPI3_crop_year.csv, SPI3_drought_yield_by_system.csv,
# and SPI3_wheat_drought_tillage.csv. Uses no imputed annual yields.
# Required CRAN packages: ggplot2, patchwork
# Run from the folder containing this script and three CSV files:
#   Rscript plot_SPI3_drought_FST.R
# Outputs: Figure_SPI3_Drought_Yields.{png,svg,tiff},
#          Figure_S_SPI3_Wheat_Tillage.{png,svg,tiff}
stopifnot(requireNamespace("ggplot2", quietly=TRUE),
          requireNamespace("patchwork",quietly=TRUE))

csvdir <- getwd()
req <- c("SPI3_crop_year.csv", "SPI3_drought_yield_by_system.csv",
         "SPI3_wheat_drought_tillage.csv")
missing <- req[!file.exists(file.path(csvdir,req))]
if(length(missing)) stop("Required files missing: ",paste(missing,collapse=", "))
spi <- read.csv(file.path(csvdir,req[1]),stringsAsFactors=FALSE)
yields <- read.csv(file.path(csvdir,req[2]),stringsAsFactors=FALSE)
till <- read.csv(file.path(csvdir,req[3]),stringsAsFactors=FALSE)
stopifnot(nrow(spi)==132L, all(yields$SPI3<=-1),
          all(till$SPI3<=-1))

library(ggplot2)
library(patchwork)
pal <- c(CNV="#475569",LEG="#CB8644",MNR="#15857A")
base <- theme_bw(base_size=10) + theme(
  panel.grid.minor=element_blank(),
  panel.grid.major.x=element_blank(),
  strip.background=element_blank(),
  strip.text=element_text(face="bold"),
  axis.text.x=element_text(size=8),
  plot.title=element_blank(),
  legend.position="bottom")

# One summer SPI-3 series shared by corn and soybean, and one spring series for wheat.
series <- rbind(
  transform(subset(spi,Crop=="corn"),Season="June–August (corn and soybean)"),
  transform(subset(spi,Crop=="wheat"),Season="March–May (wheat)")
)
series$Season <- factor(series$Season,
 levels=c("June–August (corn and soybean)","March–May (wheat)"))
p_spi <- ggplot(series,aes(x=Year,y=SPI3))+
  annotate("rect",xmin=2007.5,xmax=2024.5,ymin=-Inf,ymax=Inf,
           fill="#E9EFF5",alpha=0.7)+
  geom_hline(yintercept=-1.0,linetype="dashed",colour="#BE4C3B",linewidth=.45)+
  geom_hline(yintercept=0,colour="#D5DDE5",linewidth=.3)+
  geom_vline(xintercept=2008,linetype="dotted",colour="#7E8A98")+
  geom_line(colour="#758596",linewidth=.45)+
  geom_point(aes(colour=SPI3<=-1),size=1.45)+
  scale_colour_manual(values=c(`FALSE`="#A9B5C3",`TRUE`="#BE4C3B"),guide="none")+
  scale_x_continuous(breaks=c(1981,1990,2000,2008,2016,2024))+
  coord_cartesian(ylim=c(-3.25,2.5))+
  facet_wrap(~Season,nrow=1)+
  labs(x="Harvest year",y="SPI-3")+base

yields$System <- factor(yields$System,levels=c("CNV","LEG","MNR"))
yields$Crop <- factor(yields$Crop,levels=c("corn","soybean","wheat"),
                    labels=c("Corn","Soybean","Wheat"))
years <- unique(yields[,c("Year","Crop","SPI3")])
# Crop-specific discrete labels preserve actual missing crop–system years.
years <- years[order(years$Year),,drop=FALSE]
years$Year_SPI <- sprintf("%d\n(%.2f)",years$Year,years$SPI3)
yields$Year_SPI <- sprintf("%d\n(%.2f)",yields$Year,yields$SPI3)
yields$Year_SPI <- factor(yields$Year_SPI, levels=unique(years$Year_SPI))
p_yield <- ggplot(yields,aes(x=Year_SPI,y=Yield_Mg_ha,colour=System))+
  geom_point(size=2.8,position=position_dodge(width=.55))+
  scale_color_manual(values=pal)+
  facet_wrap(~Crop,nrow=1,scales="free")+
  labs(x="Harvest year (SPI-3 in parentheses)",y=expression(paste("Observed yield (Mg ",ha^{-1},")")),
       color=NULL)+base+
  theme(panel.grid.major.x=element_blank(),axis.text.x=element_text(size=8))

# Crop-wise observed drought yields below the time series.
main_fig <- (p_spi / p_yield)+plot_layout(heights=c(.9,1.3),guides="collect") &
            theme(legend.position="bottom")
ggsave("Figure_SPI3_Drought_Yields.png",main_fig,width=13.6,height=8.5,dpi=400,bg="white")
ggsave("Figure_SPI3_Drought_Yields.svg",main_fig,width=13.6,height=8.5,bg="white")
ggsave("Figure_SPI3_Drought_Yields.tiff",main_fig,width=13.6,height=8.5,
       dpi=450,compression="lzw",bg="white")

# Descriptive wheat tillage contrasts in drought years with both treatments observed.
till$System <- factor(till$System,levels=c("LEG","MNR"))
till$Year_label <- ifelse(till$Year==2012,"2012 (SPI-3 = −1.18)",
                       "2023 (SPI-3 = −1.53)")
till$Year_label <- factor(till$Year_label,
                        levels=c("2012 (SPI-3 = −1.18)","2023 (SPI-3 = −1.53)"))
p_till <- ggplot(till,aes(x=Year_label,y=RT_minus_FT_Mg_ha,colour=System))+
  geom_hline(yintercept=0,linetype="dashed",colour="#7E8A98")+
  geom_point(size=3.7,position=position_dodge(width=.36))+
  scale_colour_manual(values=pal[c("LEG","MNR")])+
  labs(x="Spring drought harvest year",y=expression(paste("RT − FT wheat yield (Mg ",ha^{-1},")")),
       colour=NULL)+base
for(e in c("png","svg","tiff")) {
  if(e=="tiff") ggsave("Figure_S_SPI3_Wheat_Tillage.tiff",p_till,width=7.5,
                      height=4.5,dpi=450,compression="lzw",bg="white")
  else ggsave(paste0("Figure_S_SPI3_Wheat_Tillage.",e),p_till,
              width=7.5,height=4.5,dpi=if(e=="png")400 else 72,bg="white")
}
cat("Figure export completed. Drought-event and wheat tillage plots are descriptive only.\n")
