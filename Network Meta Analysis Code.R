####Network Meta Analysis
library(readxl)
library(circlize)
library(export)
library(coda)
library(gemtc)
library(rjags)
library(ggplot2)
library(tidyr)
library(dplyr)
library(forestploter)
library(extrafont)
library(grid)
library(car)
library(meta)
library(metafor)
library(RColorBrewer)
library(gridExtra)
library(netmeta)
library(ggrepel)
library(ggforce)
library(scales)
library(tibble)

###Chord Figure
Chord_Frame <- read_excel(" ")
Chord_Description <- read_excel(" ")
circos.par(start.degree = 90, clock.wise = TRUE)
chordDiagram(Chord_Frame, 
             transparency = 0.7,
             grid.col = c("A"="#B71B33", "B"="#EB5A35", "C"="#E9AE00", "D"="#FFFF77", "E"="#98C51A", "F"="#00AA58", "G"="#2A7278", "H"="#6DC6D6", "I"="#0097CE", "J"="#28286E", "K"="#C5A900", "L"="#B96F1D", "M"="#D5006A", "N"="#A23E92", "O"="#9D9D9E"),
             order = c("E","F","G","H","I","J","K","L","M","N","O","A","B","C","D"),
             link.sort = TRUE,
             annotationTrack = c("name","grid")
)
plot <- recordPlot()
graph2office(x=plot,file=" ",type=c("ppt"), width = 9, height = 9)


###Results
##Survival Data with HR in K-M Curve
Network_Frame <- read_excel(" ", sheet = " ")
Network_Description <- read_excel(" ", sheet = " ")
Network_Frame$diff[Network_Frame$diff == "NA"] <- NA
Network_Frame$diff <- as.numeric(Network_Frame$diff)
Network_Frame$std.err[Network_Frame$std.err == "NA"] <- NA
Network_Frame$std.err <- as.numeric(Network_Frame$std.err)
network <- mtc.network(data.re = Network_Frame, treatments = Network_Description)
model <- mtc.model(network, likelihood ="binom", link = "cloglog", type = "consistency", linearModel = "random", dic = TRUE)
cat(model$code)
result <- mtc.run(model, n.adapt = 10000, n.iter = 50000, thin = 10)
gelman.diag(result)
par(mar=c(1,1,1,1))
gelman.plot(result)
summary(result)
plot(result)

##Survival Data with Events in Different Arms
Network_Frame <- read_excel(" ", sheet = " ")
Network_Description <- read_excel(" ", sheet = " ")
network <- mtc.network(Network_Frame, treatments = Network_Description)
model <- mtc.model(network, likelihood ="binom", link = "cloglog", type = "consistency", linearModel = "random", dic = TRUE, om.scale = 2, re.prior.sd = 2)
cat(model$code)
result <- mtc.run(model, n.adapt = 10000, n.iter = 50000, thin = 10)
gelman.diag(result)
par(mar=c(1,1,1,1))
gelman.plot(result)
summary(result)
plot(result)

##Efficacy Data with Events in Different Arms
Network_Frame <- read_excel(" ", sheet = " ")
Network_Description <- read_excel(" ", sheet = " ")
network <- mtc.network(Network_Frame, treatments = Network_Description)
model <- mtc.model(network, likelihood ="binom", link = "log", type = "consistency", linearModel = "random", dic = TRUE, om.scale = 2, re.prior.sd = 2)
cat(model$code)
result <- mtc.run(model, n.adapt = 10000, n.iter = 50000, thin = 10)
gelman.diag(result)
par(mar=c(1,1,1,1))
gelman.plot(result)
summary(result)
plot(result)

#Safety Data with Events in Different Arms
Network_Frame <- read_excel("D:/桌面/血液内科/BTK抑制剂网状meta/网状meta分析数据框架.xlsx", sheet = "SAE")
Network_Description <- read_excel("D:/桌面/血液内科/BTK抑制剂网状meta/网状meta分析数据描述.xlsx", sheet = "SAE")
network <- mtc.network(Network_Frame, treatments = Network_Description)
model <- mtc.model(network, likelihood ="binom", link = "log", type = "consistency", linearModel = "random", dic = TRUE, om.scale = 2, re.prior.sd = 2)
cat(model$code)
result <- mtc.run(model, n.adapt = 10000, n.iter = 50000, thin = 10)
gelman.diag(result)
par(mar=c(1,1,1,1))
gelman.plot(result)
summary(result)
plot(result)

#Safety Data Adjusted by Exposure Time
Network_Frame <- read_excel(" ", sheet = " ")
Network_Description <- read_excel(" ", sheet = " ")
network <- mtc.network(Network_Frame, treatments = Network_Description)
model <- mtc.model(network, likelihood ="poisson", link = "log", type = "consistency", linearModel = "random", dic = TRUE, om.scale = 2, re.prior.sd = 2)
cat(model$code)
result <- mtc.run(model, n.adapt = 10000, n.iter = 50000, thin = 10)
gelman.diag(result)
par(mar=c(1,1,1,1))
gelman.plot(result)
summary(result)
plot(result)
#Network Plot
networkplot <- plot(network,use.description = T, dynamic.edge.width = T, 
                    vertex.label.color = "black", vertex.label.dist=2.5, 
                    vertex.label.degree=pi/2, vertex.label.font.family="Arial", 
                    vertex.label.cex=0.8, vertex.color="red")

#Comparison Table
table <- relative.effect.table(result)
table <- round(exp(table),2)
print(table)

#SUCRA Data
ranks <- rank.probability(result,preferredDirection = 1)
print(ranks)
sucra(ranks)
print(rank.quantiles(ranks))

#SUCRA Plot
plotdata <- read_excel(" ", sheet = " ")
plot_data <- plot_data %>%
  arrange(desc(SUCRA)) %>% 
  mutate(Treatment = factor(Treatment, levels = rev(Treatment)))
sucraplot <- ggplot(plot_data) +
  geom_point(aes(x = SUCRA, y = Treatment), color = "black", size = 4) +
  geom_point(aes(x = X_base, y = Treatment), color = "red", size = 0) +
  geom_segment(aes(x = X_base, xend = SUCRA, y = Treatment, yend = Treatment), color = "black", size = 1) +
  geom_text(aes(x = SUCRA, y = Treatment, label = sprintf("%.5f", SUCRA)), color = "black", size = 3, hjust = -0.5) + 
  labs(x = " ", y = " ") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.00), limits = c(0, 1)) + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black"),   
    axis.title = element_text(size = 12),   
    axis.text = element_text(size = 10)    
  )
print(sucraplot)

##Forest Plot Comparing to Assigned Treatment Arm
forestdata <- read_excel(" ")
forestdata$` ` <- paste(rep(" ", 20), collapse = " ")
forestdata$`HR (95% CI)` <- ifelse(is.na(forestdata$HR), "", sprintf("%.2f (%.2f to %.2f)", forestdata$HR, forestdata$LowerCI, forestdata$UpperCI))
foresttheme <- forest_theme(base_size = 10, base_family = "Arial", ci_pch = 15,  ci_alpha = 1, ci_fill = NULL, ci_col = "black",
                            ci_lty = 1, ci_lwd = 1, legend_position = "none", refline_gp = gpar(lwd = 1, lty = "solid", col = "#AD002A"))


#Forest Plot for Difference Comparisons
result_anohe <- mtc.anohe(network,n.adapt=10000,n.iter=50000, likelihood ="binom", link = "cloglog", linearModel = "random")
summary_anohe <- summary(result_anohe)
summary_anohe
plot(summary_anohe, digit = 3)

#Publication Bias
library(readxl)
library(ggplot2)
biasdata <- read_excel(" ", sheet = " ")
t_result <- t.test(biasdata$`HR Difference`)
mean_value <- mean(biasdata$`HR Difference`)
se_value <- (t_result$conf.int[2] - t_result$conf.int[1]) / (2 * qt(0.975, df = t_result$parameter))
iasdata$SE_inv <- 1 / biasdata$SElnHR...17
biasdata$SND <- biasdata$`HR Difference` / biasdata$SElnHR...17
egger_model <- lm(SND ~ SE_inv, data = biasdata)
summary(egger_model)

#Funnel Plot
slope_lower <- (se_value - 0) / (t_result$conf.int[1] - 0)
slope_upper <- (se_value - 0) / (t_result$conf.int[2] - 0)
par(family = "Arial")
polygon_data <- data.frame(
  x = c(-0.7 / slope_lower, 0.7 / slope_lower, 0), 
  y = c(0.7, 0.7, 0)  
)
category_names <- c(" ", " ",  " ", " ", " ", " ", " ", " ", " "," ", " ", " ", " ", " ", " ", " ", " ", " ")
category_colors <- c(
  "#F4C1CA", "#D3A3C9", "#CA9A8E", "#FF8674", "#D0006F", "#D22630", "#653279",
  "#A4343A", "#FFC56E", "#F4DA40", "#97D700", "#009B77", "#6BA539", "#00ABAB",
  "#0092BC", "#008C95", "#004B87", "#B58500", "#9C4F01")
color_map <- setNames(category_colors, category_names)
funnelplot <- ggplot(biasdata, aes(x = biasdata$`HR Difference`, y = biasdata$SElnHR...17)) +
  geom_polygon(data = polygon_data, aes(x = x, y = y), fill = "white", alpha = 1) +
  geom_point(aes(fill = Comparison), shape = 21, color = "black", size = 3) +  
  geom_vline(xintercept = 0, linetype = "solid", color = "red") +   
  geom_abline(slope = slope_lower, intercept = 0, linetype = "dashed", color = "black") + 
  geom_abline(slope = -slope_lower, intercept = 0, linetype = "dashed", color = "black") +  
  scale_x_continuous(limits = c(-1.8, 1.8)) +  
  scale_y_reverse(limits = c(0.7, 0), expand = c(0, 0)) +  
  scale_fill_manual(values = color_map, name = "Comparison", limits = category_names) +
  theme_minimal(base_family = "Arial") + 
  theme(
    panel.grid = element_blank(),  
    axis.ticks = element_line(color = "black"),  
    axis.text = element_text(size = 12, face = "bold"), 
    axis.title = element_text(size = 12, face = "bold"), 
    axis.line = element_line(color = "black", size = 0.2),
    panel.background = element_rect(fill = "gray90", color = NA),  
    legend.position = "null"
  ) +
  guides(fill = guide_legend(keyheight = unit(0.5, "cm"),   
                             keywidth = unit(0.5, "cm"))) +
  annotate("text", x = -1.8, y = 0.01, label = "p = (Egger)", size = 5, hjust = 0, vjust = 1, fontface = "bold", family = "Arial", color = "black") +
  labs(x = "ln(HR) centered at comparison-specific pooled ln(HR)",
       y = "Standard Error of ln(HR)")
funnelplot
graph2office(x=funnelplot, file="PFS-funnelplot",type=c("PPT"), width = 12, height = 6)

#Inconsistency Model and Heterogeneity Test
mtc.nodesplit.comparisons(network)
result_consistency <- mtc.nodesplit(network,comparisons = mtc.nodesplit.comparisons(network),linearModel = "random", n.adapt=10000, n.iter=50000, likelihood ="binom", link = "log")
summary_consistency <- summary(result_consistency)
summary_consistency
plot(summary_consistency, digit = 3)
summary(result_consistency$)

##Similarity Test
Similarity_data <- read_excel(" ", sheet = " ")
cols_to_check <- colnames(Similarity_data)[
  !(colnames(Similarity_data) %in% c("Author", "Line", "DrugsA", "Drugs")) 
]
results <- data.frame(
  Variable = character(),
  Shapiro_P = numeric(),
  Test = character(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)
for (col in cols_to_check) {
  data_subset <- Similarity_data[, c("Drugs", col)]
  colnames(data_subset) <- c("Drugs", "Value")
  data_subset <- data_subset %>% filter(!is.na(Value)) 
  shapiro_p <- tryCatch(shapiro.test(data_subset$Value)$p.value, error = function(e) NA)
  if (!is.na(shapiro_p) && shapiro_p > 0.05) {
    aov_model <- aov(Value ~ Drugs, data = data_subset)
    anova_p <- summary(aov_model)[[1]][["Pr(>F)"]][1]
    test_used <- "ANOVA"
    final_p <- anova_p
  } else {
    kw_p <- kruskal.test(Value ~ Drugs, data = data_subset)$p.value
    test_used <- "Kruskal-Wallis"
    final_p <- kw_p
  }
  results <- rbind(results, data.frame(
    Variable = col,
    Shapiro_P = round(shapiro_p, 4),
    Test = test_used,
    P_value = round(final_p, 4)
  ))
}
print(results)

#Boxplot for Similarity Plot
Similarity_data <- read_excel(" ", sheet = " ")
boxplotdata <- Similarity_data[ ,c(" "," ")]
boxplotdata$  <- factor(boxplotdata$ , levels = rev(c(" ", " ", " ", " ", " ", " ", " ", " ",  " ", " ", " ", " ",  " ", " ", " ")))
boxplot <- ggplot(boxplotdata, aes(x = ` `, y = ` `)) + 
  geom_boxplot(fill = "lightblue",
               alpha = 0.6,                     
               outlier.color = "black",           
               outlier.shape = 16,
               lwd = 0.35,
               fatten = 1) + 
  coord_flip() +
  scale_y_continuous(breaks = seq( , by = 1)) +
  xlab("Treatment Group") +
  ylab(" ") +
  theme_minimal() +                           
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(size = 0.5, color = scales::alpha("gray", 0.2), linetype = "solid"),
        axis.text.y = element_text(family = "Arial", size = 6, face = "bold" , color = "black"),
        axis.text.x = element_text(family = "Arial", size = 6, face = "bold", margin = margin(b = 10), color = "black"),
        axis.title.y = element_text(family = "Arial", size = 8, face = "bold"),
        axis.title.x = element_text(family = "Arial", size = 8, face = "bold"),
        legend.position = "none")
boxplot

##Publication Bias (for classic meta analysis)
publicationbiasdata <- read_excel(" ", sheet = " ")
m.netmeta <- netmeta(TE = TE, seTE = seTE, treat1 = treatment1, treat2 = treatment2, studlab = study, data = publicationbiasdata, sm = "HR", random = T, common =  F, details.chkmultiarm = TRUE, sep.trts = "vs", tol.multiarm = 0.01)
summary(m.netmeta)
metabias(m.netmeta,pooled = "random", order = c("1","2","3","4","5","6","10","12","7","8","9","11","13","14","15"),method.bias = "Egger")
par(mar = c(5, 4, 4, 2))
funnel(m.netmeta, 
       order = c("1", "2","3","4","5","6","10","12","7","8","9","11","13","14","15"), 
       pooled = "random", 
       method.bias = "Egger" ,
       xlab = "Hazard Ratio centered at comparison-specific effect",
       ylab = "Standard Error",
       pos.tests = "topleft",
       pch = 21)
HRdata <- rma(yi=lnHR, sei=selnHR, method='DL', slab = paste(name))
category_names <- c(" ", " ",  " ", " ", " ", " ", " ", " ", " "," ", " ", " ", " ", " ", " ", " ", " ", " ", " ")
category_colors <- c(
  "#F4C1CA", "#D3A3C9", "#CA9A8E", "#FF8674", "#D0006F", "#D22630", "#653279",
  "#A4343A", "#FFC56E", "#F4DA40", "#97D700", "#009B77", "#6BA539", "#00ABAB",
  "#0092BC", "#008C95", "#004B87", "#B58500", "#9C4F01")
color_map <- setNames(category_colors, category_names)
names(category_colors) <- unique(HRdata$Category)
par(family = "Arial")
funnel(HRdata,  
       level=c(90, 95, 99),  
       shade=c("white","gray", "darkgray"), 
       pch = 21,  
       bg = color_map[publicationbiasdata$Comparison],
       digits = 2, 
       steps = 5, 
       xlab = "Hazard Ratio", ylab = "selnHR", 
       atransf=exp, 
       at=log(c(0.3, 0.6, 1, 1.6, 2.7, 4.5, 7.4, 12.5)),
       refline = 0, 
       refline2 =  , 
       back=c("lightgray"), 
       lty = 1, lty2 = 3, 
       col = color_map[publicationbiasdata$Comparison], 
       addtau2 = FALSE, 
       label = "all", 
       offset = 0.5, 
       cex = 0.7, 
       font = 2, 
       family = "Arial", 
       cex.axis = 1.2, 
       font.axis = 2, 
       font.main = 2,
       cex.lab = 1.2, 
       font.lab = 2) 
legend("topleft",  
       legend = c("0.1 > p > 0.05", "0.05 > p > 0.01", "p < 0.01"),  
       fill = c("white", "gray", "darkgray"), 
       cex = 0.7,  
       text.font = 2,  
       box.lwd = 0.5, 
       box.col = "black", 
       title = "P-value Ranges", 
       title.cex = 1,  
       title.font = 2,  
       xjust = 0, yjust = 1,) 
funnelplot <- recordPlot()
par(family = "Arial")
legend("bottom",
       legend = category_names,    
       fill = category_colors, 
       cex = 0.5,  
       text.font = 2, 
       box.lwd = 0.5,  
       box.col = "black",  
       title = "Treatment Comparisons",  
       title.cex = 1,  
       title.font = 2, 
       xjust = 0.5, yjust = 1, 
       ncol = 1) 
legendplot <- recordPlot() 
metabias(lnHR, selnHR, method.bias = "Egger")


#Bubble Plot (Used in concomitant analysis )
Network_Bubble <- read_excel(" ", sheet = " ")
Network_Bubble$HR_inv <- 1 / Network_Bubble$HR_PFS
Network_Bubble$Sig_Label <- with(Network_Bubble, ifelse(
  (LCI_RR_SAE > 1 | HCI_RR_SAE < 1) & (LCI_RR_SeriousAE > 1 | HCI_RR_SeriousAE < 1), "*#",
  ifelse((LCI_RR_SAE > 1 | HCI_RR_SAE < 1), "*",
         ifelse((LCI_RR_SeriousAE > 1 | HCI_RR_SeriousAE < 1), "#", "")
  )
))
group_colors <- c(
)
bunnelplot <- ggplot(
  Network_Bubble, aes(x =  , y =  , size = HR_inv, fill = Group)) +
  geom_point(color = "black", stroke = 1, shape = 21, alpha = 0.8) +
  geom_text_repel(
    aes(label = Group),
    size = 3,
    color = "black",
    family = "Arial",
    max.overlaps = Inf,
    box.padding = 0.2,
    point.padding = 0.2,
    segment.color = "grey60",
    segment.size = 0.3,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_size_continuous(
    range = c(4, 16),
    name = " ",
    guide = guide_legend(override.aes = list(fill = "grey50"))
  ) +
  geom_vline(xintercept = 1, linetype = "solid", color = "grey50") +
  geom_hline(yintercept = 1, linetype = "solid", color = "grey50") +
  geom_abline(slope = 1, intercept = 0, color = "#B12222", linetype = "dashed", size = 0.6) +
  coord_cartesian(xlim = c(0.6, 1.2), ylim = c(0.7, 1.7)) +
  labs(
    title = NULL,
    x = " ",
    y = " "
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Arial"),
    axis.text = element_text(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
bunnelplot

#Radar Plot (Used in Concomitant Analysis)
Network_Radar <- read_excel(" ", sheet = " ")
radar_long <- Network_Radar %>%
  pivot_longer(-Group, names_to = "Feature", values_to = "Value")
features <- unique(radar_long$Feature)
n_features <- length(features)
angle_map <- tibble(
  Feature = features,
  angle = seq(0, 2 * pi, length.out = n_features + 1)[- (n_features + 1)]
)
min_radius <- 0.5
radar_plot_data <- radar_long %>%
  mutate(Value_capped = pmin(pmax(Value, 0), 6)) %>%
  left_join(angle_map, by = "Feature") %>%
  mutate(
    radius = min_radius + Value_capped,
    x = radius * sin(angle),
    y = radius * cos(angle)
  )
circle_data <- tibble(
  x0 = 0,
  y0 = 0,
  r = min_radius + 0:6
)
labels_data <- tibble(
  label = 0:6,
  r = min_radius + 0:6
)
group_colors <- c(
  
)
groups_no_na <- radar_plot_data %>%
  group_by(Group) %>%
  filter(!any(is.na(Value_capped))) %>%
  ungroup() %>%
  distinct(Group) %>%
  pull(Group)
radar_plot_no_na <- radar_plot_data %>%
  filter(Group %in% groups_no_na)
radarplot <- ggplot() +
  geom_circle(data = circle_data, aes(x0 = x0, y0 = y0, r = r),
              inherit.aes = FALSE, color = "grey80", linetype = "dotted") +
  geom_polygon(data = radar_plot_no_na, aes(x = x, y = y, group = Group, color = Group ),
               fill = NA, size = 1, na.rm = TRUE , alpha = 0.5) +
  geom_point(data = radar_plot_data, aes(x = x, y = y, color = Group), size = 2, na.rm = TRUE) +
  geom_text(
    data = angle_map %>% mutate(x = (min_radius + 7) * sin(angle), y = (min_radius + 7) * cos(angle)),
    aes(x = x, y = y, label = Feature),
    inherit.aes = FALSE,
    size = 3,
    family = "Arial",
    fontface = "bold"
  ) +
  geom_text(
    data = labels_data,
    aes(x = 0, y = r, label = label),
    inherit.aes = FALSE,
    hjust = -0.2,
    size = 3,
    family = "Arial",
    fontface = "bold"
  ) +
  scale_color_manual(values = group_colors) +
  coord_fixed(xlim = c(-(min_radius + 7), min_radius + 7), ylim = c(-(min_radius + 7), min_radius + 7)) +
  theme_void() +
  theme(
    legend.position = "none",
    text = element_text(family = "Arial", face = "bold")
  ) +
  ggtitle(" ") +
  theme(plot.title = element_text(family = "Arial", face = "bold"))
radarplot

##Network Meta Regression
#Survival Data with HR in K-M Curve
covariatedata <- read_excel(" ", sheet = " ")
networkdata <- read_excel(" ", sheet = " ")
networkdata$diff[networkdata$diff == "NA"] <- NA
networkdata$diff <- as.numeric(networkdata$diff)
networkdata$std.err[networkdata$std.err == "NA"] <- NA
networkdata$std.err <- as.numeric(networkdata$std.err)
regressiondata <- covariatedata %>% select(study, `Years`)
regressiondata <- regressiondata[regressiondata$study %in% networkdata$study, ]
network <- mtc.network(data.re =networkdata, studies = regressiondata)
model <- mtc.model(network, likelihood ="binom", link = "cloglog", type = "regression", regressor = list(coefficient = "exchangeable", variable = "Years", control = "D"), linearModel = "random")
regressionresult <-  mtc.run(model, n.adapt = 10000, n.iter = 50000, thin = 10)
summary(regressionresult)

#Efficacy Data with Events in Different Arms
networkdata <- read_excel(" ", sheet = " ")
Network_Description <- read_excel(" ", sheet = " ")
regressiondata <- covariatedata %>% select(study, `Years`)
regressiondata <- regressiondata[regressiondata$study %in% networkdata$study, ]
network <- mtc.network(networkdata, studies = regressiondata)
model <- mtc.model(network, type = "regression", regressor = list(coefficient = "exchangeable", variable = "Years", control = "D"), linearModel = "random")
regressionresult <-  mtc.run(model, n.adapt = 10000, n.iter = 50000, thin = 10)
summary(regressionresult)

#Safety Data with Events in Different Arms
networkdata <- read_excel(" ", sheet = " ")
Network_Description <- read_excel(" ", sheet = " ")
regressiondata <- covariatedata %>% select(study, `Years`)
regressiondata <- regressiondata[regressiondata$study %in% networkdata$study, ]
network <- mtc.network(networkdata, studies = regressiondata)
model <- mtc.model(network, type = "regression", regressor = list(coefficient = "exchangeable", variable = "Years", control = "D"), linearModel = "random")
regressionresult <-  mtc.run(model, n.adapt = 10000, n.iter = 50000, thin = 10)
summary(regressionresult)

#Safety Data Adjusted by Exposure Time
networkdata <- read_excel(" ", sheet = " ")
Network_Description <- read_excel(" ", sheet = " ")
regressiondata <- covariatedata %>% select(study, `Years`)
regressiondata <- regressiondata[regressiondata$study %in% networkdata$study, ]
network <- mtc.network(networkdata, studies = regressiondata)
model <- mtc.model(network, , likelihood ="poisson", link = "log", type = "regression", regressor = list(coefficient = "exchangeable", variable = "Years", control = "D"), linearModel = "random")
regressionresult <-  mtc.run(model, n.adapt = 10000, n.iter = 50000, thin = 10)
summary(regressionresult)
samples_matrix <- as.matrix(regressionresult$samples)
main_effects <- grep("^d\\.", colnames(samples_matrix), value = TRUE)
results <- data.frame(
  comparison = main_effects,
  Q2.5 = NA,
  Q50 = NA,
  Q97.5 = NA
)
for (i in seq_along(main_effects)) {
  effect_samples <- samples_matrix[, main_effects[i]]
  quantiles_log <- quantile(effect_samples, probs = c(0.025, 0.5, 0.975))
  
  results$Q2.5[i] <- exp(quantiles_log[1])
  results$Q50[i]  <- exp(quantiles_log[2])
  results$Q97.5[i]<- exp(quantiles_log[3])
}
print(results)

##Data Cleaning for Subgroup Analysis
covariatedata <- read_excel(" ", sheet = " ")
line1_studies <- covariatedata %>% 
  filter(Line == 2) %>% 
  pull(study) %>% 
  unique()
file_path <- " "
sheets <- excel_sheets(file_path)
filtered_data_list <- map(set_names(sheets), function(sheet_name) {
  sheet_data <- read_excel(file_path, sheet = sheet_name)
  filtered <- sheet_data %>% filter(study %in% line1_studies)
  return(filtered)
})
output_path <- " "
write_xlsx(filtered_data_list, path = output_path)
frame_path <- " "
desc_path <- " "
sheet_names <- excel_sheets(frame_path)
filtered_desc_list <- map(set_names(sheet_names), function(sheet_name) {
  frame_data <- read_excel(frame_path, sheet = sheet_name)
  treatments <- unique(frame_data$treatment)
  desc_data <- read_excel(desc_path, sheet = sheet_name)
  filtered_desc <- desc_data %>% filter(id %in% treatments)
  return(filtered_desc)
})
output_path <- "D:/桌面/血液内科/BTK抑制剂网状meta/网状meta分析数据描述-2.xlsx"
write_xlsx(filtered_desc_list, path = output_path)