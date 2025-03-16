#  Charger les bibliothèques
library(survival)  # package principal pour les analyses de survie
library(survminer) # package pour améliorer la visualisation des résultats de survie
library(dplyr)     # package pour la manipulation des données
library(ggplot2)  #  package pour des visualisations graphiques avancées

#  Charger le dataset Haberman avec le bon séparateur
data = read.csv("Haberman.csv", header=TRUE, sep=";")

#  Vérification des dimensions et structure du dataset
print(dim(data))  # affiche le nombre de lignes et colonnes du dataset
print(colnames(data))  # affiche les noms des colonnes
head(data)  # Aperçu des premières lignes

#  Vérification que le dataset contient bien 4 colonnes comme attendu
if (ncol(data) != 4) {
  stop(" Erreur : Le fichier ne contient pas 4 colonnes ! Vérifiez le séparateur.")
}

# Attribution de noms explicites aux colonnes
colnames(data) = c("Age", "Year", "Nodes", "Survival_Status")

#  Statistiques descriptives pour comprendre les données
summary(data)
print(table(data$Survival_Status))  

#  Transformer Survival_Status : 1 = survie (censuré), 2 = décès (événement)
# Dans l'analyse de survie, l'événement est généralement codé comme 1 (décès) et 0 (censuré)
data$Event = ifelse(data$Survival_Status == 2, 1, 0)

# Dans le dataset Haberman, on a :
# Survival_Status = 1 signifie que le patient a survécu 5 ans ou plus après la chirurgie
# Survival_Satus = 2 signifie que le patient est décédé dans les 5 ans suivant la chirurgie
# mais dans R, Event = 1 indique que l'événement d'intérêt et donc le décès s'est produit
# Event = 0 indique que l'observation est censurée c'est à dire que l'événement d'intérêt n'a pas été observé pendant la période de suivi
# donc si Survival_Status = 2 alors Event = 1; si Survival_Status = 1, alors Event = 0

#  Vérifier que les variables sont bien numériques pour éviter les erreurs car les fonctions d'analyse statistique exigent des données numériques correctement formatées
data$Age = as.numeric(data$Age)
data$Nodes = as.numeric(data$Nodes)
data$Year = as.numeric(data$Year)
data$Survival_Status = as.numeric(data$Survival_Status)

#  Création d'une variable catégorielle pour l'âge
data$Age_Group = cut(data$Age, breaks = c(20, 40, 60, 90), 
                      labels = c("20-40", "41-60", "61-90"))
# cela permet de comparer les courbes de survie entre différents tranches d'âge

#  Elimination des valeurs aberrantes ou extrêmes
data = data[data$Age > 20 & data$Age <= 90 & data$Nodes >= 0 & data$Nodes <= 30, ]
# donc ici, âge entre 20 et 90 ans; le nombre de ganglions entre 0 et 30

#  Définition de la variable de survie
# Dans le contexte de cette étude, "Year" représente l'année de suivi après chirurgie
# Surv(time, event) crée un objet pour l'analyse de survie où:
# - time est le temps jusqu'à l'événement
# - event indique si l'événement s'est produit (1) ou non (0)
surv_object = Surv(time = data$Year, event = data$Event)

#  Vérifier que l'objet de survie est bien construit
print(summary(surv_object))

#  Estimation de Kaplan-Meier
# La formule S(t) = ∏(t_i≤t) (1 - d_i/n_i) est utilisée où:
# - d_i est le nombre d'événements au temps t_i
# - n_i est le nombre de sujets à risque juste avant t_i
km_fit = survfit(surv_object ~ 1, data = data)
print(summary(km_fit))
# ~1 indique qu'on ne stratifie pas encore par une variable. Comme vu dans le cours, l'estimateur Kaplan Meier
# est une méthode non paramétrique qui calcule la probabilité de survie à chaque point de temps où un évenement se produit

#  Calcul du temps médian de survie
median_survival = summary(km_fit)$table["median"]
cat("Temps médian de survie:", median_survival, "années\n")

#  Visualisation de Kaplan-Meier
ggsurvplot(km_fit, 
           data = data, 
           pval = TRUE,    # affiche la p-value
           conf.int = TRUE, # ajoute les intervalles de confiance
           risk.table = TRUE, # ajoute une table du nombre de patients à risque
           title = "Courbe de survie Kaplan-Meier globale",
           xlab = "Année après chirurgie",
           ylab = "Probabilité de survie",
           break.time.by = 5,  # intervalle des graduations sur l'axe X
           palette = "lancet",  # palette de couleurs 
           ggtheme = theme_bw()) # thème du graph

#  Comparaison des groupes selon le nombre de ganglions
# Basée sur la littérature médicale, un seuil cliniquement significatif
data$Nodes_Group = ifelse(data$Nodes > 5, "Élevé (>5)", "Faible (≤5)")
data$Nodes_Group = factor(data$Nodes_Group)

#  Vérifier la répartition des groupes
print(table(data$Nodes_Group))

#  Courbes de survie pour les groupes de ganglions
km_fit_nodes = survfit(surv_object ~ Nodes_Group, data = data)

# Visualisation avancée avec ggsurvplot
ggsurvplot(km_fit_nodes, 
           data = data, 
           pval = TRUE,             # Affiche la p-valeur du test de Log-Rank
           conf.int = TRUE,         # Ajoute les intervalles de confiance
           risk.table = TRUE,       # Ajoute une table de risque
           risk.table.col = "strata", # Colorie la table de risque par groupe
           linetype = "strata",     # Lignes différentes par groupe
           title = "Comparaison de la survie par nombre de ganglions",
           xlab = "Année après chirurgie",
           ylab = "Probabilité de survie",
           palette = "jco",         # Palette de couleurs
           legend.labs = c("Faible (≤5)", "Élevé (>5)"), # Étiquettes de légende
           legend.title = "Nombre de ganglions",
           ggtheme = theme_bw())    # Thème ggplot2

#  Test du Log-Rank pour comparer les courbes de survie
# H0: Pas de différence entre les courbes de survie
# H1: Il existe une différence entre les courbes
log_rank_test = survdiff(surv_object ~ Nodes_Group, data = data)
print(log_rank_test)

# Calcul manuel du p-value pour interprétation
p_value = 1 - pchisq(log_rank_test$chisq, df = length(log_rank_test$n) - 1)
cat("P-value du test Log-Rank:", p_value, "\n")

#  Ajout d'une analyse par groupe d'âge
km_fit_age = survfit(surv_object ~ Age_Group, data = data)

ggsurvplot(km_fit_age, 
           data = data, 
           pval = TRUE, 
           conf.int = TRUE, 
           risk.table = TRUE,
           title = "Comparaison de la survie par groupe d'âge",
           xlab = "Année après chirurgie",
           ylab = "Probabilité de survie",
           palette = "npg",
           ggtheme = theme_bw())

#  Vérifier la colinéarité entre variables
correlation_matrix = cor(data[, c("Age", "Year", "Nodes")])
print("Matrice de corrélation:")
print(correlation_matrix)
# permet d'identifier d'eventuelles relations linéaires fortes qui pourraient affecter le modèle de Cox

#  Modèle de Cox (Cox Proportional Hazards Model)
# h(t) = h0(t) * exp(β1X1 + β2X2 + ... + βPXP)
# où h0(t) est le risque de base et les βi sont les coefficients estimés
cox_model = coxph(surv_object ~ Age + Nodes, data = data, 
                   control = coxph.control(iter.max = 1000))
# le modèle de Cox est semi paramétrique et permet d'estimer les hazards ratios

# Résumé du modèle avec interprétation des coefficients
cox_summary <- summary(cox_model)
print(cox_summary)

#  Visualisation des Hazard Ratios
ggforest(cox_model, data = data,
         main = "Hazard Ratios du modèle de Cox",
         fontsize = 0.8)

#  Test de l'hypothèse des risques proportionnels
# Vérifie l'hypothèse fondamentale du modèle de Cox
ph_test = cox.zph(cox_model)
print(ph_test)

# Visualisation des résidus de Schoenfeld
plot(ph_test)

#  Prédiction de survie pour différents profils de patients
# Créer des profils types
newdata = data.frame(
  Age = c(45, 65),
  Nodes = c(2, 10)
)

# Calcul des courbes de survie prédites
predicted_surv = survfit(cox_model, newdata = newdata)

# Visualisation des prédictions
ggsurvplot(predicted_surv, 
           data = newdata,
           conf.int = TRUE,
           palette = "lancet",
           title = "Survie prédite pour différents profils de patients",
           xlab = "Année après chirurgie",
           ylab = "Probabilité de survie",
           legend.labs = c("45 ans, 2 ganglions", "65 ans, 10 ganglions"),
           ggtheme = theme_bw())

#  Calcul du C-index pour évaluer la capacité prédictive du modèle
# Le C-index varie de 0.5 (pas mieux que le hasard) à 1 (discrimation parfaite)
c_index = cox_summary$concordance[1]  
cat("Indice de concordance (C-index):", c_index, "\n")
