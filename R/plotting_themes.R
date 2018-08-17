
theme_sticky <- function() {
  (ggthemes::theme_foundation()
   + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                       size = ggplot2::rel(1.2), hjust = 0.5),
                    text = ggplot2::element_text(),
                    panel.background = ggplot2::element_rect(colour = NA),
                    plot.background = ggplot2::element_rect(colour = NA),
                    panel.border = ggplot2::element_rect(colour = NA),
                    axis.title = ggplot2::element_text(face = "bold",size = ggplot2::rel(1)),
                    axis.title.y = ggplot2::element_text(angle=90,vjust =2),
                    axis.title.x = ggplot2::element_text(vjust = -0.2),
                    axis.text = ggplot2::element_text(), 
                    axis.line = ggplot2::element_line(colour="black"),
                    axis.ticks = ggplot2::element_line(),
                    panel.grid.major = ggplot2::element_line(colour="#f0f0f0"),
                    panel.grid.minor = ggplot2::element_blank(),
                    legend.key = ggplot2::element_rect(colour = NA),
                    legend.position = "right",
                    legend.direction = "vertical",
                    legend.key.size= ggplot2::unit(0.2, "cm"),
                    legend.spacing = ggplot2::unit(0, "cm"),
                    legend.title = ggplot2::element_text(face="italic"),
                    plot.margin=ggplot2::unit(c(10,5,5,5),"mm"),
                    strip.background=ggplot2::element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                    strip.text = ggplot2::element_text(face="bold")
   ))
}

theme_CN_comparison <- function() {
  
  (ggthemes::theme_foundation()
   + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                       size = ggplot2::rel(1.2), hjust = 0.5),
                    text = ggplot2::element_text(),
                    panel.background = ggplot2::element_rect(colour = NA),
                    plot.background = ggplot2::element_rect(colour = NA),
                    panel.border = ggplot2::element_rect(colour = NA),
                    axis.title = ggplot2::element_text(face = "bold",size = ggplot2::rel(1)),
                    axis.title.y = ggplot2::element_text(angle=90,vjust =2),
                    axis.title.x = ggplot2::element_text(vjust = -0.2),
                    axis.text = ggplot2::element_text(), 
                    axis.line = ggplot2::element_blank(),
                    axis.ticks = ggplot2::element_line(),
                    panel.grid.major = ggplot2::element_blank(),
                    panel.grid.minor = ggplot2::element_blank(),
                    plot.margin=ggplot2::unit(c(4,4,4,4),"mm"),
                    strip.background=ggplot2::element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                    strip.text = ggplot2::element_text(face="bold")
   ))
  
}

scale_fill_sticky <- function(values, ...){
  ggplot2::discrete_scale("fill","Publication",manual_pal(values = values), ...)
}

scale_colour_sticky <- function(values, ...){
  ggplot2::discrete_scale("colour","Publication",manual_pal(values = values), ...)
}
