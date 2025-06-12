library(dplyr)
library(tidyr)

# Data loading
data <- read.csv("https://raw.githubusercontent.com/KotaHirobe/hirobe_et_al_2025_LEE/refs/heads/main/250113data.csv")
exp <- read.csv("https://raw.githubusercontent.com/KotaHirobe/hirobe_et_al_2025_LEE/refs/heads/main/exp250113.csv")

result <- data %>%
  dplyr::group_by(No, anim) %>%
  dplyr::summarize(count = n(), .groups = "drop")


# Fox data analysis ####
fox_data <- data.frame(No = 1:16)

# Data merging
fox_data <- fox_data %>%
  left_join(
    result %>% filter(anim == "fox"),
    by = "No"
  )　%>%
  mutate(count = replace_na(count, 0))

print(fox_data)

fox_ana <- merge(fox_data, exp, by = "No")

# Geological data loading
log <- read.csv("https://raw.githubusercontent.com/KotaHirobe/hirobe_et_al_2025_LEE/refs/heads/main/log.csv")
salvage <- read.csv("https://raw.githubusercontent.com/KotaHirobe/hirobe_et_al_2025_LEE/refs/heads/main/salvage.csv")
geo <- rbind(log, salvage)
colnames(geo)[colnames(geo) == "name"] <- "No"

# Data merging
fox_glm <- merge(fox_ana, geo, by = "No")
fox_glm$area <- as.factor(fox_glm$area)

### Moran's I calculation
library(spdep)
library(sp)
coords <- data.frame(x = fox_glm$x, y = fox_glm$y)

knn <- knearneigh(coords, k = 15)
nb <- knn2nb(knn)
listw <- nb2listw(nb, style = "W") 

moran.test(fox_glm$count, listw = listw)


###
# Multicollinearity check
fox_glm_excluded <- subset(fox_glm, select = -c(anim, fid, x, y, No, count, shrub, vines, litter, moss, mouse, day, number_insect))
fox_glm_binary <- as.data.frame(lapply(fox_glm_excluded, function(col) {
  if(is.character(col) || is.factor(col)) {
    as.numeric(factor(col)) -1
  } else {
    col
  }
}))
print(fox_glm_binary)
cor(fox_glm_binary)



# VIF check
library(MASS)
library(performance)
fox_res <- glm.nb(count ~ 
                    area + 
                    family_insect +
                    ferns +
                    bare +
                    dwarf_bamboos +
                    hpa +
                    mouse_p_a +
                    hight +
                    offset(log(day)),
                    data = fox_glm)
check_collinearity(fox_res)

fox_res <- glm.nb(count ~ 
                    area + 
                    family_insect +
                    ferns +
                    dwarf_bamboos +
                    hight +
                    offset(log(day)),
                  data = fox_glm)
check_collinearity(fox_res)
summary(fox_res)


library(performance)
check_overdispersion(fox_res)
check_zeroinflation(fox_res)


# Best model selection
library(MuMIn)
options(na.action = "na.fail")  # dredgeのためにNA処理のオプションを変更
fox_models <- dredge(fox_res, fixed = ~ offset(log(day)))

print(fox_models)


fox_best_model <- get.models(fox_models, 1)[[1]]
summary(fox_best_model)


library(performance)
vif(fox_best_model)
check_overdispersion(fox_best_model)
check_zeroinflation(fox_best_model)

# Predict animal counts by area type
library(ggplot2)
library(ggsignif)
library(dplyr)

newdata <- data.frame(
  area = c("legacy", "plantation"),
  family_insect = mean(fox_glm$family_insect, na.rm = TRUE),
  hight = mean(fox_glm$hight, na.rm = TRUE),
  day = mean(fox_glm$day, na.rm = TRUE)
)

pred <- predict(fox_best_model, newdata = newdata, type = "link", se.fit = TRUE)
newdata$fit <- exp(pred$fit)
newdata$lwr <- exp(pred$fit - 1.96 * pred$se.fit)
newdata$upr <- exp(pred$fit + 1.96 * pred$se.fit)

ggplot(newdata, aes(x = area, y = fit)) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.1, linewidth = 1) +
  labs(x = "Habitat Type", y = "Predicted Count",
       title = "Red fox") +
  scale_x_discrete(labels = c("legacy" = "Unsalvaged", "plantation" = "Salvaged")) +
  ylim(0, max(newdata$upr) * 1.2) +  
  theme_classic(base_size = 16) +    
  geom_signif(
    comparisons = list(c("legacy", "plantation")),
    annotations = "**",
    y_position = max(newdata$upr) * 1.1,  
    tip_length = 0.01,
    textsize = 8                        
  )



# Raccoon dog data analysis ####
rdog_data <- data.frame(No = 1:16)


rdog_data <- rdog_data %>%
  left_join(
    result %>% filter(anim == "rdog"),
    by = "No"
  )　%>%
  mutate(count = replace_na(count, 0))


print(rdog_data)

# Data merging
rdog_ana <- merge(rdog_data, exp, by = "No")
rdog_glm <- merge(rdog_ana, geo, by = "No")
rdog_glm$area <- as.factor(rdog_glm$area)

###
# Moran's I calculation
library(spdep)
library(sp)
moran.test(rdog_glm$count, listw = listw)



#glm
library(MASS)
library(performance)
rdog_local <- glm.nb(count ~ 
                    area + 
                    family_insect +
                    ferns +
                    bare +
                    dwarf_bamboos +
                    hpa +
                    mouse_p_a +
                    hight +
                    offset(log(day)),
                  data = rdog_glm)
check_collinearity(rdog_local)

rdog_local <- glm.nb(count ~ 
                       area + 
                       family_insect +
                       ferns +
                       dwarf_bamboos +
                       hpa +
                       mouse_p_a +
                       hight +
                       offset(log(day)),
                     data = rdog_glm)
check_collinearity(rdog_local)

summary(rdog_local)


# Best model selection
library(MuMIn)
options(na.action = "na.fail") 
rdog_models <- dredge(rdog_local, fixed = ~ offset(log(day)))

print(rdog_models)


rdog_best_model <- get.models(rdog_models, 1)[[1]]
summary(rdog_best_model)


library(performance)
check_overdispersion(rdog_best_model)
check_zeroinflation(rdog_best_model)

# Predicted raccoon dog counts by area type
library(MASS)
library(ggplot2)
library(ggsignif)

rdog_newdata <- data.frame(area = c("legacy", "plantation"),
                           day = mean(rdog_glm$day))  

rdog_pred <- predict(rdog_best_model, rdog_newdata, type = "link", se.fit = TRUE)
rdog_newdata$fit <- exp(rdog_pred$fit)
rdog_newdata$lwr <- exp(rdog_pred$fit - 1.96 * rdog_pred$se.fit)
rdog_newdata$upr <- exp(rdog_pred$fit + 1.96 * rdog_pred$se.fit)

ggplot(rdog_newdata, aes(x = area, y = fit)) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.1, linewidth = 1) +
  labs(x = "Habitat Type", y = "Predicted Count",
       title = "Raccoon dog") +
  scale_x_discrete(labels = c("legacy" = "Unsalvaged", "plantation" = "Salvaged")) +
  ylim(0, max(rdog_newdata$upr) * 1.2) +  
  theme_classic(base_size = 16) +    
  geom_signif(
    comparisons = list(c("legacy", "plantation")),
    annotations = "*",  
    y_position = max(rdog_newdata$upr) * 1.1,
    tip_length = 0.01,
    textsize = 8
  )


# Raccoon data analysis####
raccoon_data <- data.frame(No = 1:16)


raccoon_data <- raccoon_data %>%
  left_join(
    result %>% filter(anim == "raccoon"),
    by = "No"
  )　%>%
  mutate(count = replace_na(count, 0))

print(raccoon_data)

# Data merging
raccoon_ana <- merge(raccoon_data, exp, by = "No")
raccoon_glm <- merge(raccoon_ana, geo, by = "No")
raccoon_glm$area <- as.factor(raccoon_glm$area)

### 
# Moran's I calculation
moran.test(raccoon_glm$count, listw = raccoon_listw)

#glm
library(MASS)
library(performance)
raccoon_local <- glm.nb(count ~ 
                       area + 
                       family_insect +
                       ferns +
                       bare +
                       dwarf_bamboos +
                       hpa +
                       mouse_p_a +
                       hight +
                       offset(log(day)),
                     data = raccoon_glm)
check_collinearity(raccoon_local)

# Best model selection
library(MuMIn)
options(na.action = "na.fail")  
raccoon_models <- dredge(raccoon_local, fixed = ~ offset(log(day)))

print(raccoon_models)

# Best model
raccoon_best_model <- get.models(raccoon_models, 1)[[1]]
summary(raccoon_best_model)


library(performance)

check_collinearity(raccoon_best_model)
check_overdispersion(raccoon_best_model)
check_zeroinflation(raccoon_best_model)

library(MASS)
library(ggplot2)
library(ggsignif)

raccoon_newdata <- data.frame(area = c("legacy", "plantation"),
                              family_insect = mean(raccoon_glm$family_insect, na.rm = TRUE),
                              hight = mean(raccoon_glm$hight, na.rm = TRUE),
                              day = mean(raccoon_glm$day, na.rm = TRUE))  

raccoon_pred <- predict(raccoon_best_model, raccoon_newdata, type = "link", se.fit = TRUE)
raccoon_newdata$fit <- exp(raccoon_pred$fit)
raccoon_newdata$lwr <- exp(raccoon_pred$fit - 1.96 * raccoon_pred$se.fit)
raccoon_newdata$upr <- exp(raccoon_pred$fit + 1.96 * raccoon_pred$se.fit)

ggplot(raccoon_newdata, aes(x = area, y = fit)) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.1, linewidth = 1) +
  labs(x = "Habitat Type", y = "Predicted Count",
       title = "Raccoon") +
  scale_x_discrete(labels = c("legacy" = "Unsalvaged", "plantation" = "Salvaged")) +
  ylim(0, max(raccoon_newdata$upr) * 1.2) +
  theme_classic(base_size = 16) +
  geom_signif(
    comparisons = list(c("legacy", "plantation")),
    annotations = "*",  # p = 0.01722 → 有意水準 0.05 未満でアスタリスク1つ
    y_position = max(raccoon_newdata$upr) * 1.1,
    tip_length = 0.01,
    textsize = 8
  )


# Local factor plotting ####
library(ggplot2)
library(ggplot2)

library(ggplot2)

ggplot(fox_glm, aes(x = area, y = weight_insect)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("legacy" = "Unsalvaged", "plantation" = "Salvaged")) +
  xlab("Area") +
  ylab("Weight of Insect") +
  theme_classic()

ggplot(fox_glm, aes(x = area, y = family_insect)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("legacy" = "Unsalvaged", "plantation" = "Salvaged")) +
  xlab("Area") +
  ylab("Insect families") +
  theme_classic()

ggplot(fox_glm, aes(x = area, y = ferns)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("legacy" = "Unsalvaged", "plantation" = "Salvaged")) +
  xlab("Area") +
  ylab("Cover of ferns") +
  theme_classic()

ggplot(fox_glm, aes(x = area, y = herbaceous_plants)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("legacy" = "Unsalvaged", "plantation" = "Salvaged")) +
  xlab("Area") +
  ylab("Cover of herbaceous plants") +
  theme_classic()

ggplot(fox_glm, aes(x = area, y = dwarf_bamboos)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("legacy" = "Unsalvaged", "plantation" = "Salvaged")) +
  xlab("Area") +
  ylab("Cover of dwarf bamboos") +
  theme_classic()

ggplot(fox_glm, aes(x = area, y = bare)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("legacy" = "Unsalvaged", "plantation" = "Salvaged")) +
  xlab("Area") +
  ylab("Cover of bare ground") +
  theme_classic()

ggplot(fox_glm, aes(x = area, y = hight)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("legacy" = "Unsalvaged", "plantation" = "Salvaged")) +
  xlab("Area") +
  ylab("Maximum vegetation height") +
  theme_classic()

ggplot(fox_glm, aes(x = area, y = hpa)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("legacy" = "Unsalvaged", "plantation" = "Salvaged")) +
  xlab("Area") +
  ylab("Volume of logjams") +
  theme_classic()

ggplot(subset(fox_glm, mouse_p_a == "presence"), aes(x = area)) +
  geom_bar(fill = "gray") +
  scale_x_discrete(labels = c("legacy" = "Unsalvaged", "plantation" = "Salvaged")) +
  xlab("Area") +
  ylab("Number of mouse presences") +
  theme_classic()
