# Packages -------------------
library(tidyverse)
library(readxl)
library(see)
library(cowplot)
library(ggthemes)
library(vegan)
library(iNEXT)
library(mvabund)
library(RInSp)

# Data ------------------
data = read_excel("Mimagoniates_diet.xlsx", sheet = "geral")
str(data)

# Summary -------------
data$Sexo = as.factor(data$Sexo)
summary(data$Sexo)

data$CP_mm = data$CP * 10

resultados_CP <- data %>%
  group_by(Sexo) %>% # Agrupa pelos níveis do fator
  summarise(
    Minimo = min(CP_mm, na.rm = TRUE),
    Maximo = max(CP_mm, na.rm = TRUE),
    Media = mean(CP_mm, na.rm = TRUE),
    Desvio_Padrao = sd(CP_mm, na.rm = TRUE)
  )

resultados_peso <- data %>%
  group_by(Sexo) %>% 
  summarise(
    Minimo = min(Peso, na.rm = TRUE),
    Maximo = max(Peso, na.rm = TRUE),
    Media = mean(Peso, na.rm = TRUE),
    Desvio_Padrao = sd(Peso, na.rm = TRUE)
  )



# Sexual proportion -------------
tab_sexo <- table(data$Sexo)
tab_sexo

chisq.test(tab_sexo, p = c(0.5, 0.5))
binom.test(x = 37, n = 163, conf.level = 0.95)
binom.test(x = 126, n = 163, conf.level = 0.95)

ks.test(
  x = data$CP[data$Sexo == "M"],
  y = data$CP[data$Sexo == "F"]
)


# Size differences (sex) ---------------
t.test(CP_mm ~ Sexo, data = data) # t = 0.29, p = 0.76
t.test(Peso ~ Sexo, data = data) # t = 1.26, p = 0.20

p1 = data %>% 
  ggplot(aes(x = Sexo, y = CP_mm, fill = Sexo))+
  geom_jitter(width = 0.08, alpha = 0.5,
              shape = 21, color = "black",
              show.legend = FALSE)+
  geom_boxplot(width = 0.2,, alpha = 0.6,
               show.legend = FALSE)+
  geom_violinhalf(
    position = position_nudge(x = 0.15) ,
    alpha = 0.5,
    show.legend = FALSE
  )+
  scale_fill_manual(values = c("darkorchid", "cyan4")) +
  theme_tufte(base_size = 15)+
  theme(axis.line = element_line(colour = "grey50", linewidth = 0.5))+
  labs(x = "Sex", y = "Standard length (mm)")+
  scale_y_continuous(limits = c(10, 50), breaks = seq(10, 50, by = 10))
  

p2 = data %>% 
  ggplot(aes(x = Sexo, y = Peso, fill = Sexo))+
  geom_jitter(width = 0.08, alpha = 0.5,
              shape = 21, color = "black",
              show.legend = FALSE)+
  geom_boxplot(width = 0.2,, alpha = 0.6,
               show.legend = FALSE)+
  geom_violinhalf(
    position = position_nudge(x = 0.15) ,
    alpha = 0.5,
    show.legend = FALSE
  )+
  scale_fill_manual(values = c("darkorchid", "cyan4")) +
  theme_tufte(base_size = 15)+
  theme(axis.line = element_line(colour = "grey50", linewidth = 0.5))+
  labs(x = "Sex", y = "Weight (g)")+
  scale_y_continuous(limits = c(0, 1.6), breaks = seq(0, 1.5, by = .5))

p3 = plot_grid(p1, p2, labels = "AUTO", nrow = 1)
ggsave("p3.jpg", height =  5, width = 10)

# SR --------------------
gr_prop <- data %>%
  dplyr::count(GR) %>% 
  dplyr::mutate(prop = round(n/sum(n), 4)*100)

p2 = gr_prop %>% 
  ggplot(aes(x = 2, y = prop, fill = factor(GR)))+
  geom_bar(stat = "identity", color = "black",
           alpha = 0.8) +
  geom_text(aes(label = paste0(prop, "%")), color = "white",
            position = position_stack(vjust = .5), size = 6) +
  xlim(0, 2.5) +
  coord_polar(theta = "y", start = 0) +
  theme_void() +
  theme(legend.position = c(.5, .5),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15)) +
  labs(fill = "SR")+
  scale_fill_brewer(palette = "Set1")


# DIET --------------------
dieta = read_excel("Mimagoniates.xlsx", sheet = "ordem")
dieta$CP_mm <- dieta$CP * 10

str(dieta)
GPA = colSums(dieta[, 3:23]) / 124

diet_pa = decostand(dieta[,3:23], method = "pa")

freq_abs <- colSums(diet_pa)
freq_pct <- (freq_abs / 124) * 100

# iNEXT ----------------------------
datlist <- list()
com_inext <- diet_pa
#com_inext$ambiente <- env_com$ambiente
#com_inext <- arrange(com_inext, ambiente)

com_inext_mima <- decostand(diet_pa,
                            method = "pa")
# com_inext_rd <- decostand(com_inext[c(25:48),-21],
                          method = "pa")

datlist$Mimagoniates <-data.frame(t(com_inext_mima))
# datlist$Poças <-data.frame(t(com_inext_poca))

datlist

?iNEXT

result <- iNEXT(datlist,
                q = 0,
                datatype = "incidence_raw",
                endpoint = 200,
                se = TRUE, 
                nboot = 999)



fig2 = ggiNEXT(result, type = 1)+
  scale_colour_manual(values = c("darkorchid4")) +
  scale_fill_manual(values = c("darkorchid4")) +
  theme_tufte(base_size = 15)+
  theme(legend.position = "none", 
        axis.line = element_line(colour = "grey50", linewidth = 0.5))+
  labs(x = "Number of stomachs analyzed", 
       y = "Number of food items") 

ggsave("Figure_2.jpg", fig2)

# Multivariate GLM -----------------------

freq_occ <- colSums(diet_pa) / nrow(diet_pa)
keep_cols <- freq_occ >= 0.05
diet_pa_filtered <- diet_pa[, keep_cols]
ncol(diet_pa_filtered)
freq_occ[keep_cols]

abund <- mvabund(diet_pa_filtered)

mod <- manyglm(
  abund ~ CP_mm, 
  data = dieta,      
  family = "binomial"
)

summary(mod)

anova_res <- anova(mod)
anova_res

coef(mod)
coefs <- coef(mod)

slopes <- coefs["CP_mm", ]

slope_df <- data.frame(
  Item  = names(slopes),
  Slope = as.numeric(slopes)
)

nomes_novos <- c(
  "Detritus", 
  "Plant fragments", 
  "Fruits", 
  "Seeds", 
  "Fragments of terrestrial insects", 
  "Hymenoptera (adults)", 
  "Diptera (adults)", 
  "Diptera (larvae)", 
  "Araneae", 
  "Trombidiformes",
  "Coleoptera (adults)",
  "Ostracoda",
  "Collembola"
)

slope_df$item_new = nomes_novos

fig_3 = ggplot(slope_df, aes(x = reorder(item_new, Slope), y = Slope, fill = Slope > 0)) +
  geom_bar(stat = "identity", show.legend = F, alpha = 0.9,
           color = "black") +
  coord_flip() +
  scale_fill_manual(values = c("red", "steelblue"), 
                    labels = c("Negativo", "Positivo")) +
  labs(
    x = NULL,
    y = "Slope (Effect of Standard Length \n on Presence Probability)"
  ) +
  theme_tufte(base_size = 15)+
  theme(legend.position = "none", 
        axis.line = element_line(colour = "grey50", linewidth = 0.5))+
  scale_y_continuous(limits = c(-0.13, 0.155), 
                     breaks = seq(-0.1, 0.15, by = 0.05))

ggsave("Figure_3.jpg", fig_3)



plot(mod, which = 1:2, legend.pos = "none", main = "")

# Extrair resíduos e valores ajustados do modelo
residuals_mod <- residuals(mod, type = "response")
fitted_mod <- fitted(mod)

# Criar DataFrame com os dados necessários
diagnostics_df <- data.frame(
  Fitted = as.vector(fitted_mod),
  Residuals = as.vector(residuals_mod)
)

# Gráfico A: Resíduos vs Valores Ajustados
plot_A <- ggplot(diagnostics_df, aes(x = Fitted, y = Residuals)) +
  geom_point(alpha = 0.5, size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted Values", y = "Residuals") +
  theme_tufte(base_size = 15)+
  theme(axis.line = element_line(colour = "grey50", linewidth = 0.5))

plot_A

# Gráfico B: QQ Plot
# Combinar resíduos de todos os itens em um único vetor
residuals_combined <- as.vector(residuals_mod)

# Gerar valores teóricos e amostrais para o QQ plot
qq_vals <- qqnorm(residuals_combined, plot.it = FALSE)

# Criar o data frame para ggplot2
qq_data <- data.frame(
  Theoretical = qq_vals$x,
  Sample = qq_vals$y
)

# Criar o gráfico QQ com ggplot2
plot_B <- ggplot(qq_data, aes(x = Theoretical, y = Sample)) +
  geom_point(alpha = 0.5, size = 2.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_tufte(base_size = 15)+
  theme(axis.line = element_line(colour = "grey50", linewidth = 0.5))

# Exibir o gráfico
plot_B

# Combinar os dois gráficos lado a lado
Figure_4 = plot_grid(plot_A, plot_B, labels = "AUTO")

ggsave("Figure_4.jpg", Figure_4, height =  5, width = 10)
