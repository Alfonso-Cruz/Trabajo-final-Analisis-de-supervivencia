
# Apéndice

################################################################################
### Leer los datos
################################################################################

library(haven)
library(survival)
library(knitr)
library(kableExtra)
library(tidyverse)

datos <- read_dta("divorce3.dta")


################################################################################
### Creación de variables y selección de primeros matrimonios
################################################################################

datos$tiempo <- (datos$lower + datos$upper) / 2
datos$divorcio <- 1 - datos$censor

primeros_matrimonios <- subset(datos, marnum == 1)

table(datos$marnum)

prop.table(table(primeros_matrimonios$divorcio)) * 100

table(primeros_matrimonios$divorcio)

summary(primeros_matrimonios$tiempo)


################################################################################
### Estimación Kaplan - Meier
################################################################################


ajuste_km <- survfit(Surv(tiempo, divorcio) ~ 1, data = primeros_matrimonios)

plot(ajuste_km,
     xlab = "Duración del matrimonio en años",
     ylab = "Probabilidad de supervivencia",
     main = "Estimación de Kaplan-Meier",
     conf.int = TRUE,
     mark.time = TRUE)


################################################################################
### Log-Rank test
################################################################################

km_heblack <- survfit(Surv(tiempo, divorcio) ~ heblack, 
                      data = primeros_matrimonios)

km_sheblack <- survfit(Surv(tiempo, divorcio) ~ sheblack, 
                       data = primeros_matrimonios)

par(mfrow = c(1, 2))

plot(km_heblack,
     xlab = "Duración del matrimonio (años)",
     ylab = "Probabilidad de supervivencia",
     main = "Raza del esposo",
     col = c("black", "grey50"),
     lty = c(1, 2),
     lwd = 1.5,
     mark.time = FALSE)

legend("topright",
       legend = c("No afroamericano", "Afroamericano"),
       col = c("black", "grey50"),
       lty = c(1, 2),
       lwd = 1.5,
       bty = "n",
       cex = 0.75)

plot(km_sheblack,
     xlab = "Duración del matrimonio (años)",
     ylab = "Probabilidad de supervivencia",
     main = "Raza de la esposa",
     col = c("black", "grey50"),
     lty = c(1, 2),
     lwd = 1.5,
     mark.time = FALSE)

legend("topright",
       legend = c("No afroamericana", "Afroamericana"),
       col = c("black", "grey50"),
       lty = c(1, 2),
       lwd = 1.5,
       bty = "n",
       cex = 0.75)

par(mfrow = c(1, 1))

logrank_heblack <- survdiff(Surv(tiempo, divorcio) ~heblack, 
                            data = primeros_matrimonios)

logrank_sheblack <- survdiff(Surv(tiempo, divorcio) ~sheblack, 
                             data = primeros_matrimonios)

################################################################################
### Modelo de regresión de Cox
################################################################################
modelo_cox <- coxph(Surv(tiempo, divorcio) ~ age + agediff + hiseduc + 
                      hereduc + heblack + sheblack, 
                    data = primeros_matrimonios)

summary(modelo_cox)

prueba_residuos <- cox.zph(modelo_cox)


# Covariables dependientes del tiempo

aux <- primeros_matrimonios |> select(id, spans, starts_with("time"), 
                                      starts_with("numkid"))

aux <- aux |> pivot_longer(
  cols = -c(id,spans),names_to = c(".value", "periodo"),
  names_pattern = "(time|numkid)([0-9]+)") |> 
  mutate(periodo = as.integer(periodo)) |> 
  filter(periodo <=spans)

cambios <- aux |> group_by(id) |> arrange(periodo, .by_group = TRUE) |>
  mutate(y=lag(time), x=numkid) |> 
  filter(periodo > 1) |> ungroup() |> select(id, y, x)

base <- primeros_matrimonios |> select(id, tiempo, divorcio, age, agediff, 
                                       hiseduc, hereduc, heblack, sheblack)

tmerge1 <- tmerge(base,base,
                  id = id,
                  divorcio = event(tiempo, divorcio))

tmerge2 <- tmerge(tmerge1, cambios,
                  id = id,
                  numkid = tdc(y, x, init = 0))

modelo_cox_hijos <- coxph(Surv(tstart, tstop, divorcio) ~ numkid + age + 
                            agediff + hiseduc + hereduc + heblack + sheblack, 
                          data = tmerge2)

prueba_residuos_hijos <- cox.zph(modelo_cox_hijos)

tabla_cox_hijos <- data.frame(
  Variable = names(coef(modelo_cox_hijos)),
  Coeficiente = round(coef(modelo_cox_hijos), 3),
  "Hazard ratio" = round(exp(coef(modelo_cox_hijos)), 3),
  "IC 95%" = paste0(
    "[",
    round(exp(confint(modelo_cox_hijos))[, 1], 3),
    "; ",
    round(exp(confint(modelo_cox_hijos))[, 2], 3),
    "]"
  ),
  "p-valor" = round(summary(modelo_cox_hijos)$coefficients[, "Pr(>|z|)"], 4),
  check.names = FALSE
)

kable(tabla_cox_hijos,
      caption = "Modelo de Cox con el número de hijos 
      como covariable dependiente del tiempo",
      align = "lccccc") |>
  kable_styling(full_width = FALSE)

################################################################################
### Inclusión de todos los matrimonios
################################################################################

modelo_cox_total <- coxph(Surv(tiempo, divorcio) ~ age + agediff + hiseduc + 
                            hereduc + heblack + sheblack, data = datos)

summary(modelo_cox_total)

prueba_residuos_total <- cox.zph(modelo_cox_total)

prueba_residuos_total


modelo_cox_cluster <- coxph(Surv(tiempo, divorcio) ~ age + agediff + 
                              hiseduc + hereduc + heblack + sheblack + 
                              cluster(id), data = datos)

summary(modelo_cox_cluster)


modelo_frailty <- coxph(Surv(tiempo, divorcio) ~ age + agediff + 
                          hiseduc + hereduc+ heblack + sheblack + 
                          frailty(id), data = datos)

summary(modelo_frailty)


tabla_cox_total <- data.frame(
  Variable = names(coef(modelo_cox_total)),
  Coeficiente = round(coef(modelo_cox_total), 3),
  `Hazard ratio` = round(exp(coef(modelo_cox_total)), 3),
  `IC 95%` = paste0(
    "[",
    round(exp(confint(modelo_cox_total))[, 1], 3),
    "; ",
    round(exp(confint(modelo_cox_total))[, 2], 3),
    "]"
  ),
  `p-valor` = signif(summary(modelo_cox_total)$coefficients[, "Pr(>|z|)"], 3),
  check.names = FALSE
)

kable(tabla_cox_total,
      caption = "Modelo de Cox convencional para todos los matrimonios",
      align = "lccccc") %>%
  kable_styling(full_width = FALSE)


tabla_cox_cluster <- data.frame(
  Variable = names(coef(modelo_cox_cluster)),
  Coeficiente = round(coef(modelo_cox_cluster), 3),
  `Hazard ratio` = round(exp(coef(modelo_cox_cluster)), 3),
  `IC 95%` = paste0(
    "[",
    round(exp(confint(modelo_cox_cluster))[, 1], 3),
    "; ",
    round(exp(confint(modelo_cox_cluster))[, 2], 3),
    "]"
  ),
  `p-valor` = signif(summary(modelo_cox_cluster)$coefficients[, "Pr(>|z|)"], 3),
  check.names = FALSE
)

kable(tabla_cox_cluster,
      caption = "Modelo de Cox con errores estándar 
      robustos agrupados por individuo",
      align = "lccccc") %>%
  kable_styling(full_width = FALSE)


tabla_cox_frailty <- data.frame(
  Variable = names(coef(modelo_frailty)),
  Coeficiente = round(coef(modelo_frailty), 3),
  `Hazard ratio` = round(exp(coef(modelo_frailty)), 3),
  `IC 95%` = paste0(
    "[",
    round(exp(confint(modelo_frailty))[, 1], 3),
    "; ",
    round(exp(confint(modelo_frailty))[, 2], 3),
    "]"
  ),
  `p-valor` = signif(
    summary(modelo_frailty)$coefficients[1:6, "p"],
    3
  ),
  check.names = FALSE
)

kable(tabla_cox_frailty,
      caption = "Modelo de Cox con fragilidad compartida por individuo",
      align = "lccccc") %>%
  kable_styling(full_width = FALSE)
