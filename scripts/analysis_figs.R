library(ggplot2)
##Spearman correlations
#load data
#france
fr2p5k<-read.table('res/excl_filtered_sites/merged/fmt_france4_nucdivers_vcfPi_2500_watt.txt')
fr5k<-read.table('res/excl_filtered_sites/merged/fmt_france4_nucdivers_vcfPi_5000_watt.txt')
fr50k<-read.table('res/excl_filtered_sites/merged/fmt_france4_nucdivers_vcfPi_50000_watt.txt')
fr100k<-read.table('res/excl_filtered_sites/merged/fmt_france4_nucdivers_vcfPi_100000_watt.txt')
fr500k<-read.table('res/excl_filtered_sites/merged/fmt_france4_nucdivers_vcfPi_500000_watt.txt')
fr1000k<-read.table('res/excl_filtered_sites/merged/fmt_france4_nucdivers_vcfPi_1000000_watt.txt')
#germany
ger2p5k<-read.table('res/excl_filtered_sites/merged/fmt_germany8_nucdivers_vcfPi_2500_watt.txt')
ger5k<-read.table('res/excl_filtered_sites/merged/fmt_germany8_nucdivers_vcfPi_5000_watt.txt')
ger50k<-read.table('res/excl_filtered_sites/merged/fmt_germany8_nucdivers_vcfPi_50000_watt.txt')
ger100k<-read.table('res/excl_filtered_sites/merged/fmt_germany8_nucdivers_vcfPi_100000_watt.txt')
ger500k<-read.table('res/excl_filtered_sites/merged/fmt_germany8_nucdivers_vcfPi_500000_watt.txt')
ger1000k<-read.table('res/excl_filtered_sites/merged/fmt_germany8_nucdivers_vcfPi_1000000_watt.txt')
#gough island
gi2p5k<-read.table('res/excl_filtered_sites/merged/fmt_gough_nucdivers_vcfPi_2500_watt.txt')
gi5k<-read.table('res/excl_filtered_sites/merged/fmt_gough_nucdivers_vcfPi_5000_watt.txt')
gi50k<-read.table('res/excl_filtered_sites/merged/fmt_gough_nucdivers_vcfPi_50000_watt.txt')
gi100k<-read.table('res/excl_filtered_sites/merged/fmt_gough_nucdivers_vcfPi_100000_watt.txt')
gi500k<-read.table('res/excl_filtered_sites/merged/fmt_gough_nucdivers_vcfPi_500000_watt.txt')
gi1000k<-read.table('res/excl_filtered_sites/merged/fmt_gough_nucdivers_vcfPi_1000000_watt.txt')
#concatenate data
fr_cat_df<-cbind('winsize' = c(rep(2500, times = nrow(fr2p5k)), rep(5000, times = nrow(fr5k)), rep(50000, times = nrow(fr50k)), rep(100000, times = nrow(fr100k)), rep(500000, times = nrow(fr500k)), rep(1000000, times = nrow(fr1000k))),
                               rbind(fr2p5k, fr5k, fr50k, fr100k, fr500k, fr1000k))
fr_cat_df<-cbind('population' = rep('France', times = nrow(fr_cat_df)), fr_cat_df)
ger_cat_df<-cbind('winsize' = c(rep(2500, times = nrow(ger2p5k)), rep(5000, times = nrow(ger5k)), rep(50000, times = nrow(ger50k)), rep(100000, times = nrow(ger100k)), rep(500000, times = nrow(ger500k)), rep(1000000, times = nrow(ger1000k))),
                  rbind(ger2p5k, ger5k, ger50k, ger100k, ger500k, ger1000k))
ger_cat_df<-cbind('population' = rep('Germany', times = nrow(ger_cat_df)), ger_cat_df)
gi_cat_df<-cbind('winsize' = c(rep(2500, times = nrow(gi2p5k)), rep(5000, times = nrow(gi5k)), rep(50000, times = nrow(gi50k)), rep(100000, times = nrow(gi100k)), rep(500000, times = nrow(gi500k)), rep(1000000, times = nrow(gi1000k))),
                 rbind(gi2p5k, gi5k, gi50k, gi100k, gi500k, gi1000k))
gi_cat_df<-cbind('population' = rep('Gough Island', times = nrow(gi_cat_df)), gi_cat_df)
cat_df<-rbind(fr_cat_df, ger_cat_df, gi_cat_df)
#init results dataframe
cor.res<-data.frame('population' = c(rep('France', times = 6), rep('Germany', times = 6), rep('Gough Island', times = 6)),
                    'winsize' = rep(c(2500, 5000, 50000, 100000, 500000, 1000000), times = 3),
                    'pi_avg' = numeric(18), 'pi_sd' = numeric(18), 'watt_avg' = numeric(18), 'watt_sd' = numeric(18), 'Spearman_rho_tajima' = numeric(18), 'P_val_tajima' = numeric(18), 'Spearman_rho_watterson' = numeric(18), 'P_val_watterson' = numeric(18), 'Spearman_rho_tajVwat' = numeric(18), 'P_val_tajVwat' = numeric(18))
#loop over population and windowsize, writing results to res dataframe cor.res
res_row_counter<-0
for(pop in unique(cat_df$population)){
  for(win in unique(cat_df$winsize)){
    res_row_counter<-res_row_counter + 1
    sub_df<-subset(cat_df, population == pop & winsize == win)
    #pi avg and sd
    cor.res$pi_avg[res_row_counter]<-mean(sub_df$nuc_divers, na.rm = TRUE)
    cor.res$pi_sd[res_row_counter]<-sd(sub_df$nuc_divers, na.rm = TRUE)
    #watt avg and sd
    cor.res$watt_avg[res_row_counter]<-mean(sub_df$watt_theta/sub_df$n_nonTX_sites, na.rm = TRUE)
    cor.res$watt_sd[res_row_counter]<-sd(sub_df$watt_theta/sub_df$n_nonTX_sites, na.rm = TRUE)
    #pi correlation
    res<-cor.test(sub_df$cm.Mb, sub_df$nuc_divers, method = 'spearman')
    cor.res$Spearman_rho_tajima[res_row_counter]<-res$estimate
    cor.res$P_val_tajima[res_row_counter]<-res$p.value
    #watterson's theta correlation
    res<-cor.test(sub_df$cm.Mb, sub_df$watt_theta/sub_df$n_nonTX_sites, method = 'spearman')
    cor.res$Spearman_rho_watterson[res_row_counter]<-res$estimate
    cor.res$P_val_watterson[res_row_counter]<-res$p.value
    #test correlation between watterson and tajima
    res<-cor.test(sub_df$nuc_divers, sub_df$watt_theta/sub_df$n_nonTX_sites, method = 'spearman')
    cor.res$Spearman_rho_tajVwat[res_row_counter]<-res$estimate
    cor.res$P_val_tajVwat[res_row_counter]<-res$p.value
  }
}
#write concat table
write.table(cat_df, file = '/media/tb_disk/projects/mmdom_linked_selection/res/excl_filtered_sites/master_df.txt', row.names = FALSE, sep = '\t')
#upload and concatenate new master dataframes
ch1<-read.table('data/lm_data/new_master_df_chr1.txt', header = TRUE, sep = '\t')
ch2<-read.table('data/lm_data/new_master_df_chr2.txt', header = TRUE, sep = '\t')
ch3<-read.table('data/lm_data/new_master_df_chr3.txt', header = TRUE, sep = '\t')
ch4<-read.table('data/lm_data/new_master_df_chr4.txt', header = TRUE, sep = '\t')
ch5<-read.table('data/lm_data/new_master_df_chr5.txt', header = TRUE, sep = '\t')
ch6<-read.table('data/lm_data/new_master_df_chr6.txt', header = TRUE, sep = '\t')
ch7<-read.table('data/lm_data/new_master_df_chr7.txt', header = TRUE, sep = '\t')
ch8<-read.table('data/lm_data/new_master_df_chr8.txt', header = TRUE, sep = '\t')
ch9<-read.table('data/lm_data/new_master_df_chr9.txt', header = TRUE, sep = '\t')
ch10<-read.table('data/lm_data/new_master_df_chr10.txt', header = TRUE, sep = '\t')
ch11<-read.table('data/lm_data/new_master_df_chr11.txt', header = TRUE, sep = '\t')
ch12<-read.table('data/lm_data/new_master_df_chr12.txt', header = TRUE, sep = '\t')
ch13<-read.table('data/lm_data/new_master_df_chr13.txt', header = TRUE, sep = '\t')
ch14<-read.table('data/lm_data/new_master_df_chr14.txt', header = TRUE, sep = '\t')
ch15<-read.table('data/lm_data/new_master_df_chr15.txt', header = TRUE, sep = '\t')
ch16<-read.table('data/lm_data/new_master_df_chr16.txt', header = TRUE, sep = '\t')
ch17<-read.table('data/lm_data/new_master_df_chr17.txt', header = TRUE, sep = '\t')
ch18<-read.table('data/lm_data/new_master_df_chr18.txt', header = TRUE, sep = '\t')
ch19<-read.table('data/lm_data/new_master_df_chr19.txt', header = TRUE, sep = '\t')
#cat
full_master_df<-rbind(ch1,ch2,ch3,ch4,ch5,ch6,ch7,ch8,ch9,ch10,ch11,ch12,ch13,ch14,ch15,ch16,ch17,ch18,ch19)
#write
write.table(full_master_df, file = 'data/lm_data/full_master_df.txt', sep = '\t')
#for future analyses, just load full_master_df.txt
full_master_df<-read.table('data/lm_data/full_master_df.txt', header = TRUE)
gough_1Mb<-subset(full_master_df, population == 'Gough Island' & winsize == 1000000)
france_1Mb<-subset(full_master_df, population == 'France' & winsize == 1000000)
germany_1Mb<-subset(full_master_df, population == 'Germany' & winsize == 1000000)
gough_5Kb<-subset(full_master_df, population == 'Gough Island' & winsize == 5000)
france_5Kb<-subset(full_master_df, population == 'France' & winsize == 5000)
germany_5Kb<-subset(full_master_df, population == 'Germany' & winsize == 5000)
#examine correlation between CpG and proportion of tx sites at 1Mb windows
cor.test(gough_1Mb$CpG_count, gough_1Mb$prop_genic, method = 'spearman')
#examine correlation between prop tx and pi
cor.test(gough_1Mb$prop_genic, gough_1Mb$nuc_divers, method = 'spearman')
cor.test(france_1Mb$prop_genic, france_1Mb$nuc_divers, method = 'spearman')
cor.test(germany_1Mb$prop_genic, germany_1Mb$nuc_divers, method = 'spearman')
cor.test(gough_5Kb$prop_genic, gough_5Kb$nuc_divers, method = 'spearman')
cor.test(france_5Kb$prop_genic, france_5Kb$nuc_divers, method = 'spearman')
cor.test(germany_5Kb$prop_genic, germany_5Kb$nuc_divers, method = 'spearman')
#correlation between prop tx and recombination rate
cor.test(gough_1Mb$prop_genic, gough_1Mb$cm.Mb, method = 'spearman')
cor.test(germany_1Mb$prop_genic, germany_1Mb$cm.Mb)
#fit models
model_gough_1Mb<-lm(gough_1Mb$nuc_divers ~ gough_1Mb$cm.Mb + gough_1Mb$prop_div_JC + gough_1Mb$notx_CpG_count * gough_1Mb$prop_genic)
model_germany_1Mb<-lm(germany_1Mb$nuc_divers ~ germany_1Mb$cm.Mb + germany_1Mb$prop_div_JC + germany_1Mb$notx_CpG_count * germany_1Mb$prop_genic)
model_france_1Mb<-lm(france_1Mb$nuc_divers ~ france_1Mb$cm.Mb + france_1Mb$prop_div_JC + france_1Mb$notx_CpG_count * france_1Mb$prop_genic)
model_gough_5Kb<-lm(gough_5Kb$nuc_divers ~ gough_5Kb$cm.Mb + gough_5Kb$prop_div_JC + gough_5Kb$notx_CpG_count * gough_5Kb$prop_genic)
model_germany_5Kb<-lm(germany_5Kb$nuc_divers ~ germany_5Kb$cm.Mb + germany_5Kb$prop_div_JC + germany_5Kb$notx_CpG_count * germany_5Kb$prop_genic)
model_france_5Kb<-lm(france_5Kb$nuc_divers ~ france_5Kb$cm.Mb + france_5Kb$prop_div_JC + france_5Kb$notx_CpG_count * france_5Kb$prop_genic)
#prepare data for viz -- remove rows with NAs for fitted values
#1Mb windows
viz_df_1Mb<-data.frame('Pop' = c(as.character(gough_1Mb$population[-model_gough_1Mb$na.action]), as.character(germany_1Mb$population[-model_germany_1Mb$na.action]), as.character(france_1Mb$population[-model_france_1Mb$na.action])),
                      'Pi' = c(gough_1Mb$nuc_divers[-model_gough_1Mb$na.action], germany_1Mb$nuc_divers[-model_germany_1Mb$na.action], france_1Mb$nuc_divers[-model_france_1Mb$na.action]),
                      'RR' = c(gough_1Mb$cm.Mb[-model_gough_1Mb$na.action], germany_1Mb$cm.Mb[-model_germany_1Mb$na.action], france_1Mb$cm.Mb[-model_france_1Mb$na.action]),
                      'fitted' = c(model_gough_1Mb$fitted.values, model_germany_1Mb$fitted.values, model_france_1Mb$fitted.values))
#plots for lm's. 1Mb windows.
ggplot(viz_df_1Mb, aes(x = RR)) +
  geom_point(size = 3, alpha = 0.4, aes(y = Pi)) +
  geom_line(size = 0.8, color = 'steelblue3', aes(y = fitted)) +
  xlab('Recombination Rate (cM/Mb)') + ylab(expression(pi)) +
  facet_wrap(~Pop) +
  theme_classic() +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22), strip.text = element_text(size = 22)) +
  theme(panel.spacing = unit(1, "lines"), panel.background = element_rect(color = 'black')) +
  theme(plot.margin = margin(l = 9, t = 3, b = 3, r = 20))
#wilcox rank-sum tests for hightest 5% vs lowest 5% recombination rates -- all windows, all populations
#init res dataframe
summary_df<-data.frame('Population' = c(), 'winsize' = c(), 'upper_95_cmMb' = c(), 'lower_05_cmMb' = c(),
                       'Mean_Low_Pi' = c(), 'Mean_High_Pi' = c(), 'Pi_Wval' = c(), 'Pi_pval' = c(), 'Mean_Low_Watt' = c(),
                       'Mean_High_Watt' = c(), 'watt_Wval' = c(), 'watt_pval' = c())
#loop over populations and window sizes
for(pop in unique(full_master_df$population)){
  for(win in unique(full_master_df$winsize)){
    sub_df<-subset(full_master_df, population == pop & winsize == win)
    cutoffs<-quantile(sub_df$cm.Mb, c(0.05, 0.95), na.rm = TRUE)
    low_df<-subset(sub_df, cm.Mb < cutoffs[1])
    high_df<-subset(sub_df, cm.Mb > cutoffs[2])
    wilcox_res_pi<-wilcox.test(low_df$nuc_divers, high_df$nuc_divers)
    wilcox_res_watt<-wilcox.test(low_df$watt_theta/low_df$n_nonTX_sites, high_df$watt_theta/high_df$n_nonTX_sites)
    pi_W<-wilcox_res_pi[1]
    pi_pval<-wilcox_res_pi[3]
    watt_W<-wilcox_res_watt[1]
    watt_pval<-wilcox_res_watt[3]
    add_df<-data.frame('Population' = pop, 'winsize' = win, 'upper_95_cmMb' = cutoffs[2], 'lower_05_cmMb' = cutoffs[1],
                       'Mean_Low_Pi' = mean(low_df$nuc_divers, na.rm = TRUE), 'Mean_High_Pi' = mean(high_df$nuc_divers, na.rm = TRUE),
                       'Pi_Wval' = pi_W, 'Pi_pval' = pi_pval, 'Mean_Low_Watt' = mean(low_df$watt_theta/low_df$n_nonTX_sites, na.rm = TRUE),
                       'Mean_High_Watt' = mean(high_df$watt_theta/high_df$n_nonTX_sites, na.rm = TRUE), 'watt_Wval' = watt_W, 'watt_pval' = watt_pval)
    summary_df<-rbind(summary_df, add_df)
  }
}
#write out results
write.csv(summary_df, file = '/home/mikey/projects/mmdom_linked_selection/res/tables/wilcox_res.csv')
#boxplot for 1Mb window size
#prepare dataframe
dat_1Mb<-subset(full_master_df, winsize == 1000000)
cutoffs<-quantile(dat_1Mb$cm.Mb, probs = c(0.05, 0.95), na.rm = TRUE)
dat_1Mb_rr_outliers<-subset(dat_1Mb, cm.Mb <= cutoffs[1] | cm.Mb >= cutoffs[2])
percentile<-character(nrow(dat_1Mb_rr_outliers))
percentile[which(dat_1Mb_rr_outliers$cm.Mb <= cutoffs[1])]<-'Lower 5th'
percentile[which(dat_1Mb_rr_outliers$cm.Mb >= cutoffs[2])]<-'Upper 5th'
percentile<-as.factor(percentile)
dat_1Mb_rr_outliers<-cbind(dat_1Mb_rr_outliers, percentile)
#plot
#omit 3 windows with pi > 0.007
dat_1Mb_rr_outliers<-subset(dat_1Mb_rr_outliers, nuc_divers <0.007)
ggplot(dat_1Mb_rr_outliers) +
  geom_boxplot(alpha = 0.5, size = 0.8, width = 0.75, notch = TRUE, aes(x = percentile, y = nuc_divers, fill = percentile)) +
  xlab('Percentile') + ylab(expression(pi)) +
  facet_wrap(~population) +
  theme_classic() +
  theme(legend.position = 'none', axis.text = element_text(size = 11), axis.title = element_text(size = 12), strip.text = element_text(size = 12))
#correlation between window size and rho
cor_df<-data.frame('Population' = c(), 'winsize' = c(), 'avg_pi' = c(), 'sd_pi' = c(), 'avg_rr' = c(), 'sd_rr' = c(), 'rho' = c(), 'p_val' = c())
for(pop in unique(full_master_df$population)){
  for(win in unique(full_master_df$winsize)){
    win_df<-subset(full_master_df, population == pop & winsize == win)
    avg_pi<-mean(win_df$nuc_divers,  na.rm = TRUE)
    sd_pi<-sd(win_df$nuc_divers, na.rm = TRUE)
    avg_rr<-mean(win_df$cm.Mb, na.rm = TRUE)
    sd_rr<-sd(win_df$cm.Mb, na.rm = TRUE)
    win_cor<-cor.test(win_df$cm.Mb, win_df$nuc_divers, method = 'spearman')
    rho<-win_cor$estimate
    p_val<-win_cor$p.value
    new_row<-data.frame('Population' = pop, 'winsize' = win, 'avg_pi' = avg_pi, 'sd_pi' = sd_pi, 'avg_rr' = avg_rr, 'sd_rr' = sd_rr, 'rho' = rho, 'p_val' = p_val)
    cor_df<-rbind(cor_df, new_row)
  }
}
#plot
ggplot(cor_df) +
  geom_point(size = 3, alpha = 0.75, aes(x = winsize/1000, y = rho, color = Population)) +
  xlab('Window Size (Kb)') + ylab(expression(paste("Spearman's ", rho))) +
  theme_classic()
#check covariate correlations
#Pi v Divergence (Gough, 1Mb)
gough_dat_5Kb<-subset(full_master_df, population == 'Gough Island' & winsize == 5000)
cor.test(gough_dat_5Kb$nuc_divers, gough_dat_5Kb$prop_div_JC, method = 'spearman')
gough_dat_1Mb<-subset(dat_1Mb, population == 'Gough Island')
cor.test(gough_dat_1Mb$nuc_divers, gough_dat_1Mb$prop_div_JC, method = 'spearman')
#Pi v Divergence (Germany, 1Mb)
ger_dat_1Mb<-subset(dat_1Mb, population == 'Germany')
cor.test(ger_dat_1Mb$nuc_divers, ger_dat_1Mb$prop_div_JC, method = 'spearman')
#Pi v Divergence (France, 1Mb)
fr_dat_1Mb<-subset(dat_1Mb, population == 'France')
cor.test(fr_dat_1Mb$nuc_divers, fr_dat_1Mb$prop_div_JC, method = 'spearman')
#Correlation between nucleotide diversity and fine rr (Paigen et al chr1; Billings et al chr11)
#compute recombination rate from dense map data. get function.
rr_fixedwin_single_chrom<-function(map_df, rr_win_size){
  cm_mb_df<-data.frame('window' = c(), 'start' = c(), 'end' = c(), 'n_snps' = c(), 'cM/Mb' = c())
  #the 'temp_df' naming is relic from previous ver., in which unique chromosomes were iterated over.
  temp_df<-map_df
  endpos<-signif(max(temp_df$position_bp) + 5 * 10^(floor(log10(max(temp_df$position_bp))) - 1), digits = 1)
  win_start<-seq(3000000, endpos, by = rr_win_size)
  win_stop<-seq(3000000, endpos, by = rr_win_size) + (rr_win_size - 1)
  for(win in 1:length(win_start)) {
    win_df<-dplyr::filter(temp_df, position_bp >= win_start[win] & position_bp < win_stop[win])
    if(nrow(win_df) == 0) {
      win_add<-data.frame('window' = win, 'start' = win_start[win], 'end' = win_stop[win], 'n_snps' = 0, 'cm/Mb' = NA)
      cm_mb_df<-rbind(cm_mb_df, win_add)
    }
    else {
      cm_mb<-as.numeric(summary(lm(win_df$ave_cM ~ win_df$position_bp))[4]$coefficients[2]) * 1000000
      snp_num<-nrow(win_df)
      win_add<-data.frame('window' = win, 'start' = win_start[win], 'end' = win_stop[win], 'n_snps' = snp_num, 'cm/Mb' = cm_mb)
      cm_mb_df<-rbind(cm_mb_df, win_add)
    }
    #print(paste(win, ' of ', length(win_start), ' for chrom ', chrom))
  }
  return(cm_mb_df)
}
#function to append densemap rr to master df
# WILL NEED TO MAKE ADJUSTMENTS TO ACCOMODATAE FOR NEW HEADERS
append_densemap_rr<-function(master_df, rr_df){
  densemap_rr<-numeric(nrow(master_df))
  for(i in 1:nrow(master_df)){
    master_start<-master_df$start[i]
    master_end<-master_df$stop[i]
    rr_win<-which(rr_df$start <= master_start & rr_df$end >= master_end)
    densemap_rr[i]<-rr_df$cm.Mb[rr_win]
  }
  return_df<-cbind(master_df, densemap_rr)
}
#estimate rr for paigen map
chr1_map<-read.csv('data/maps/c1_cM.csv', header = TRUE)
#adjust some header names
names(chr1_map)[4]<-'position_bp'
names(chr1_map)[ncol(chr1_map)]<-'ave_cM'
#compute rr over 5Mb windows
chr1_paigen_rr_df<-rr_fixedwin_single_chrom(chr1_map, 5000000)
#merge dfs for chromosome 1
chr1_cox_data<-subset(full_master_df, chr == 'chr1')
chr1_cox_paigen<-append_densemap_rr(chr1_cox_data, chr1_paigen_rr_df)
#check correlation between maps
ggplot(chr1_cox_paigen) +
  geom_line(size = 1, aes(x = window, y = cm.Mb)) +
  geom_line(size = 1, color = 'blue', aes(x = window, y = densemap_rr)) +
  xlab('1Mb Window') + ylab('Recombination Rate (cM/Mb)') +
  theme_classic()
#examine correlations between dense map and nucleotide diversity at 5Kb and 1Mb windows
chr1_dense_res<-data.frame('population' = c(), 'winsize' = c(), 'rho' = c(), 'pval' = c())
for(pop in unique(chr1_cox_paigen$population)){
  for(win in unique(chr1_cox_paigen$winsize)){
    sub_df<-subset(chr1_cox_paigen, population == pop & winsize == win)
    res<-cor.test(sub_df$densemap_rr, sub_df$nuc_divers, method = 'spearman')
    new_row<-data.frame('population' = pop, 'winsize' = win, 'rho' = res$estimate, 'pval' = res$p.value)
    chr1_dense_res<-rbind(chr1_dense_res, new_row)
  }
}
write.csv(chr1_dense_res, file = 'res/tables/chr1_paigen_cors.csv')
#estimate rr for bililngs map
chr11_map<-read.csv('data/maps/c11_cM.csv', header = TRUE)
#change Mb position to bp position
chr11_map$Position.Mb.B37<-chr11_map$Position.Mb.B37 * 1000000
#adjust some col names
names(chr11_map)[3]<-'position_bp'
names(chr11_map)[ncol(chr11_map)]<-'ave_cM'
#compute rr over 5Mb windows
chr11_billings_rr_df<-rr_fixedwin_single_chrom(chr11_map, 5000000)
#merge dfs for chromosome 11
chr11_cox_data<-subset(full_master_df, chr == 'chr11')
chr11_cox_billings<-append_densemap_rr(chr11_cox_data, chr11_billings_rr_df)
#check correlation between maps
ggplot(chr11_cox_billings) +
  geom_line(size = 1, aes(x = window, y = cm.Mb)) +
  geom_line(size = 1, color = 'blue', aes(x = window, y = densemap_rr)) +
  xlab('1Mb Window') + ylab('Recombination Rate (cM/Mb)') +
  theme_classic()
cor.test(chr11_cox_billings$cm.Mb, chr11_cox_billings$densemap_rr, method = 'spearman')
#examine correlations between dense map and nucleotide diversity at 5Kb and 1Mb windows
chr11_dense_res<-data.frame('population' = c(), 'winsize' = c(), 'rho' = c(), 'pval' = c())
for(pop in unique(chr11_cox_billings$population)){
  for(win in unique(chr11_cox_billings$winsize)){
    sub_df<-subset(chr11_cox_billings, population == pop & winsize == win)
    res<-cor.test(sub_df$densemap_rr, sub_df$nuc_divers, method = 'spearman')
    new_row<-data.frame('population' = pop, 'winsize' = win, 'rho' = res$estimate, 'pval' = res$p.value)
    chr11_dense_res<-rbind(chr11_dense_res, new_row)
  }
}
write.csv(chr11_dense_res, file = 'res/tables/chr11_billings_cors.csv')
#Correlation between hotspot density and nucleotide diversity
hotspot_df<-read.csv('data/Smagulovaetal201_DSB_hotspot_list.csv', header = TRUE)
#merge master df and hotspot counts for all populations, 5kb and 1Mb window sizes
#5Kb
dat_5Kb<-subset(full_master_df, winsize == 5000)
DSB_hot_count<-numeric(nrow(dat_5Kb))
counter = 0
for(pop in unique(dat_5Kb$population)){
  for(chrom in unique(dat_5Kb$chr)){
    sub_df<-subset(dat_5Kb, chr == chrom & population == pop)
    sub_hot<-subset(hotspot_df, chromosome == chrom)
    for(win in 1:nrow(sub_df)){
      counter = counter + 1
      hotct<-0
      #starts in window, ends in window
      hotct<-hotct + length(which(sub_hot$start >= sub_df$start[win] & sub_hot$end <= sub_df$stop[win]))
      #starts in window, ends outside
      hotct<-hotct + length(which(sub_hot$start >= sub_df$start[win] & sub_hot$start <= sub_df$stop[win] & sub_hot$end > sub_df$stop[win]))
      #starts outide of window, ends inside
      hotct<-hotct + length(which(sub_hot$start < sub_df$start[win] & sub_hot$end >= sub_df$start[win] & sub_hot$end <= sub_df$stop[win]))
      #starts before window, ends after window
      hotct<-hotct + length(which(sub_hot$start < sub_df$start[win] & sub_hot$end > sub_df$stop[win]))
      DSB_hot_count[counter]<-hotct
    }
  }
}
dat_5Kb<-cbind(dat_5Kb, DSB_hot_count)
pi_dsb_wilcox<-data.frame('Population' = c(), 'W' = c(), 'pval' = c())
for(pop in unique(dat_5Kb$population)){
  sub_df<-subset(dat_5Kb, population == pop)
  res<-wilcox.test(sub_df$nuc_divers[sub_df$DSB_hot_count > 0], sub_df$nuc_divers[sub_df$DSB_hot_count == 0], alternative = 'greater')
  new_row<-data.frame('Population' = pop, 'W' = res$statistic, 'pval' = res$p.value)
  pi_dsb_wilcox<-rbind(pi_dsb_wilcox, new_row)
}
pi_dsb_wilcox
#1Mb
dat_1Mb<-subset(full_master_df, winsize = 1000000)
DSB_hot_count<-numeric(nrow(dat_1Mb))
counter = 0
for(pop in unique(dat_1Mb$population)){
  for(chrom in unique(dat_1Mb$chr)){
    sub_df<-subset(dat_1Mb, chr == chrom & population == pop)
    sub_hot<-subset(hotspot_df, chromosome == chrom)
    for(win in 1:nrow(sub_df)){
      counter = counter + 1
      hotct<-0
      #starts in window, ends in window
      hotct<-hotct + length(which(sub_hot$start >= sub_df$start[win] & sub_hot$end <= sub_df$stop[win]))
      #starts in window, ends outside
      hotct<-hotct + length(which(sub_hot$start >= sub_df$start[win] & sub_hot$start <= sub_df$stop[win] & sub_hot$end > sub_df$stop[win]))
      #starts outide of window, ends inside
      hotct<-hotct + length(which(sub_hot$start < sub_df$start[win] & sub_hot$end >= sub_df$start[win] & sub_hot$end <= sub_df$stop[win]))
      #starts before window, ends after window
      hotct<-hotct + length(which(sub_hot$start < sub_df$start[win] & sub_hot$end > sub_df$stop[win]))
      DSB_hot_count[counter]<-hotct
    }
  }
}
dat_1Mb<-cbind(dat_1Mb, DSB_hot_count)
pi_dsb_cor<-data.frame('Population' = c(), 'rho' = c(), 'pval'= c())
for(pop in unique(dat_1Mb$population)){
  sub_df<-subset(dat_1Mb, population == pop)
  res<-cor.test(sub_df$DSB_hot_count, sub_df$nuc_divers)
  new_row<-data.frame('Population' = pop, 'rho' = res$estimate, 'pval' = res$p.value)
  pi_dsb_cor<-rbind(pi_dsb_cor, new_row)
}
