load("Intermediate/transients.RData")
load("Intermediate/aci.RData")
load("Intermediate/lightflecks.RData")


transients_mean = subset(transients_mean, Genotype %in% c("col","npq4", "npq1", "rwt43","rca2","spsa"))
lf_data_mean = subset(lf_data_mean, Genotype %in% c("col","npq4", "npq1", "rwt43", "rca2","spsa"))
aci_inputs_means = subset(aci_inputs_means, Genotype %in% c("col","npq4", "npq1", "rwt43", "rca2","spsa"))

transients_se = subset(transients_se, Genotype %in% c("col","npq4", "npq1", "rwt43", "rca2","spsa"))
lf_data_se = subset(lf_data_se, Genotype %in% c("col","npq4","npq1", "rwt43","rca2","spsa"))
names(lf_data_se) = names(lf_data_mean)
aci_inputs_se = subset(aci_inputs_se, Genotype %in% c("col","npq4","npq1", "rwt43","rca2","spsa"))

aci_df_genotype = subset(aci_df_genotype, Genotype %in% c("col","npq4","npq1", "rwt43","rca2","spsa"))

# Convert data into lists -------------------------------------------------
l_transients_mean = dlply(transients_mean, c("Genotype","TransientType"), function(x) x)
l_transients_se = dlply(transients_se, c("Genotype","TransientType"), function(x) x)
l_lf_data_mean = dlply(lf_data_mean, c("Genotype","Amplitude"), function(x) x)
l_lf_data_se = dlply(lf_data_se, c("Genotype","Amplitude"), function(x) x)
l_aci_inputs_mean = dlply(aci_inputs_means, c("Genotype"), function(x) x)
l_aci_inputs_se = dlply(aci_inputs_se, c("Genotype"), function(x) x)


