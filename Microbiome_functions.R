if(!require("pacman")) {
  install.packages("pacman", repos = "http://cran.us.r-project.org")
}
pacman::p_load(
  BiocManager,
  tidyverse,
  ggpubr,
  ggsci,
  rstatix,
  phyloseq,
  vegan,
  zoo,
  mixOmics,
  lemon,
  microbiome,
  cowplot,
  Maaslin2
)


####################### COLOR PALETTES #########################################
stacked_bar.palette <- c(
  "#E64B35FF",
  "#4DBBD5FF",
  "#00A087FF",
  "#3C5488FF",
  "#F39B7FFF",
  "#8491B4FF",
  "#91D1C2FF",
  "#DC0000FF",
  "#7E6148FF",
  "#B09C85FF",
  "#E4E9B2",
  "#F9A620",
  "#054A29",
  "#52414C",
  "#D81E5B",
  "#331832",
  "#27474E",
  "#573D1C",
  "#404E4D",
  "#DAD4EF",
  "#E86A92",
  "#044389",
  "#6C4B5E",
  "#4E6E58",
  "#826AED",
  "#FF0054",
  "#9E0059",
  "#387D7A",
  "#395E66",
  "#1BE7FF"
)

######################### DATA WRANGLING #########################

taxonomy <- function (ps) {
  return(as.data.frame(tax_table(ps)))
}

meta_to_df <- function(ps, rownames=T) {
  if (rownames == T){
    return(as(sample_data(ps), "data.frame")) 
  } else {
    return(as(sample_data(ps), "data.frame") %>% 
             rownames_to_column("sample_id"))
  }
  
}

ps_to_feattab <- function(ps) {
  return(as.data.frame(ps@otu_table))
}

get_top_n <- function(ps, n, level = "species", agg=F) {
  if (level != "species" & agg == T) {
    ps <- ps %>% tax_fix %>%
      tax_agg(rank = level)
    ps <- ps_get(ps)
  }
  
  topn <- ps %>%
    transform(transform = "relative") %>%
    psmelt() %>%
    group_by(OTU) %>%
    summarise(Mean_abund = mean(Abundance)) %>%
    # remove those with zero abundance
    filter(Mean_abund > 0) %>% 
    slice_max(Mean_abund, n = n) %>%
    pull(OTU)
  return(topn)
}

get_top_n_group <- function(ps, n, level = "species",var,  
                            group=NULL, agg=F) {
  if (level != "species" & agg == T) {
    ps <- ps %>% tax_fix %>%
      tax_agg(rank = level)
    ps <- ps_get(ps)
  }
  topn <- ps %>%
    transform(transform = "relative") %>%
    psmelt() %>%
    filter({{ var }} == {{ group }}) %>% 
    group_by(OTU) %>%
    summarise(Mean_abund = mean(Abundance)) %>%
    # remove those with zero abundance
    filter(Mean_abund > 0) %>%
    slice_max(Mean_abund, n = n) %>%
    pull(OTU)
  
  return(topn)
}

filter_rownames <- function(df, filt_vector) {
  # wrapper script around filter to handle rownames
  df_filt <- df %>% 
    rownames_to_column(var="id") %>%
    filter(id %in% filt_vector) %>%
    column_to_rownames(var="id")
  return(df_filt)
}

common_cols <- function(list_df, n=length(list_df)) {
  col_list <- map(list_df,names)
  
  if (n < length(list_df)){
    n_cols <- table(unlist(col_list))
    return(names(n_cols[n_cols >= n]))
  }
  
  return(Reduce(intersect,col_list))
}


combine_counts <- function(counts_list, taxa_are_cols=T, prev_thresh=0.3) {
  if (taxa_are_cols == F){
    counts_list <- map(counts_list, function(x) t(x) %>% as.data.frame)
  }
  # purely presence absence
  shared_taxa <- common_cols(counts_list, n=round(length(counts_list) * prev_thresh))
  
  combined_counts <- bind_rows(counts_list) %>% 
    dplyr::select(all_of(shared_taxa))
  return(combined_counts)
}

combine_meta <- function(meta_list, col_list) {
  
  combined_meta <- bind_rows(meta_list) %>% 
    dplyr::select(all_of(col_list))
}

t_df <- function(x) {
  return(as.data.frame(t(x)))
}

ps_mutate <- function(ps, ...) {
  
  if(is(ps, "phyloseq")){
    df <- meta_to_df(ps)
    saved_rownames <- rownames(df)
    df <- dplyr::mutate(df, ...)
    rownames(df) <- saved_rownames
    ps@sam_data <- phyloseq::sample_data(df) 
  }
  
  else {
    stop("Not a phyloseq object")
  }
  
  return(ps)
}

#########################    IMPORTING DATA #########################
read_tab_delim_metag <- function(df, comment_char="#") {
  # read all tab delimited files using these params
  df_out <-
    read.table(
      df,
      row.names = 1,
      header = 1,
      sep = "\t",
      check.names = F,
      quote = "",
      comment.char = comment_char
    )
  return(df_out)
}

load_phylo <- function(asvtab, taxa, mapping, tree = NULL) {
  # convert to phyloseq and return list
  phylo_asv <- otu_table(asvtab, taxa_are_rows = T)
  
  phylo_tax <- tax_table(as.matrix(taxa))
  
  phylo_map <- sample_data(mapping)
  
  if (exists("tree")) {
    phylo_tree <- read_tree(tree)
    return(merge_phyloseq(phylo_asv, phylo_tax, phylo_tree, phylo_map))
  }
  else {
    return(merge_phyloseq(phylo_asv, phylo_tax, phylo_map))
  }
}


return_level <- function(level) {
  if (level %in% c("Subspecies", "subspecies", "t")) {
    return("t")
  } else if (level %in% c("Species", "species", "s")) {
    return("s")
  } else if (level %in% c("Genus", "genus", "g")) {
    return("g")
  } else if (level %in% c("Family", "family", "f")) {
    return("f")
  } else if (level %in% c("Order", "order", "o")) {
    return("o")
  } else if (level %in% c("Class", "class", "c")) {
    return("c")
  } else if (level %in% c("Phylum", "phylum", "p")) {
    return("p")
  } else if (level %in% c("Kingdom", "kingdom", "k", 
                          "Domain", "domain", "d")) {
    return("k")
  } else {
    stop("Not a valid taxonomic rank")
  }
}

select_rank <- function(merged_table, level, tax_format = "metaphlan") {
  level <- return_level(level)
  
  tax_table <-  merged_table %>%
    rownames_to_column("taxa") %>%
    dplyr::select(taxa)
  
  if (tax_format == "metaphlan") {
    tax_table_sep <- tax_table %>%
      rowwise() %>%
      separate(
        taxa,
        into = c(
          "Kingdom",
          "Phylum",
          "Class",
          "Order",
          "Family",
          "Genus",
          "Species",
          "Subspecies"
        ),
        remove = F,
        sep = "\\|"
      )
  } else if (tax_format == "gtdb") {
    tax_table_sep <- tax_table %>%
      rowwise() %>%
      separate(
        taxa,
        into = c(
          "Kingdom",
          "Phylum",
          "Class",
          "Order",
          "Family",
          "Genus",
          "Species",
          "Subspecies"
        ),
        remove = F,
        sep = ";"
      )
  }
  
  tax_table_fill <- t(zoo::na.locf(t(tax_table_sep))) %>%
    as.data.frame()
  
  tax_table_out <- tax_table_fill %>%
    rowwise() %>%
    separate(Subspecies, into = c("Level", "tmp"), sep = "_") %>%
    dplyr::select(-tmp) %>% 
    filter(Level == level) %>%
    column_to_rownames("taxa") %>%
    dplyr::select(-Level)
  
  merged_table_out <- merged_table %>%
    rownames_to_column("taxa") %>%
    dplyr::filter(taxa %in% rownames(tax_table_out)) %>%
    column_to_rownames("taxa")
  
  return_list <-
    list(counts = merged_table_out, tax_table = tax_table_out)
  return(return_list)
}

import_pseq_metag <- function(merged_table_path, metapath, level, 
                              tax_format="metaphlan", table_comment_char="#", 
                              meta_comment_char="#") {
  # read files
  
  merged_table <- read_tab_delim_metag(merged_table_path, comment_char = table_comment_char)
  metadata <- read_tab_delim_metag(metapath, comment_char = meta_comment_char)
  
  merged_table_rank <- select_rank(merged_table, level=level, tax_format=tax_format)
  
  out <- load_phylo(merged_table_rank$counts, merged_table_rank$tax_table, metadata)
  
  return(out)
}

######################### NORMALISATION & FILTERING ############################

transform <- function(ps, transform = "relative", offset=1) {
  if (length(is(ps)) == 1 && class(ps) == "phyloseq") {
    x <- ps_to_feattab(ps)
    mean_sum <- mean(colSums(x))
  }
  else {
    print("not a phyloseq object, exiting")
    stop()
  }
  
  if (transform %in% c("relative", "arcsin", "log10")) {
    if (transform == "relative") {
      ps_t <- t(t(x) / colSums(x))
      
    } else if (transform == "arcsin") {
      if (mean_sum == 1) {
        ps_t <- asin(sqrt(x))
      } else if (mean_sum == 100) {
        x <- x / 100
        ps_t <- asin(sqrt(x))
      } else {
        stop("not proportion data, exiting")
      } 
    } else if (transform == "log10"){
      if (max_x == 100 | max_x == 1) {
        warning("data are relative abundances", call. = F)
        ps_t <- log10(x + offset)
      }   else {
        ps_t <- log10(x + offset)
      }
    }
    otu_table(ps)@.Data <- as.matrix(ps_t)
    
    return(ps)
    
  } else {
    print("Not a valid transform, exiting")
    stop()
  }
  
}

filter_ps <- function(ps, abun, prev=NULL) {
  
  x <- ps_to_feattab(ps) %>% 
    as.matrix()
  
  
  if (is.null(prev)) {
    ps_filt <- x[rowSums(x)  >= abun, ]
  } else if (!is.null(prev)) {
    ps_filt <- x[rowSums(x >= abun) >= ncol(x) * prev, ]
  }
  
  otu_table(ps)@.Data <- ps_filt
  
  return(ps)
}

subset_samples_func <- function(func_profile, filter_statement) {
  
  keep <- func_profile$Metadata %>% 
    filter(eval(rlang::parse_expr(filter_statement))) %>% 
    rownames()
  
  func_profile_filt <- map(func_profile, function(x) filter_rownames(x, keep))
  
  return(func_profile_filt)
  
}


prune_samples_func  <- function(func_profile, keep) {
  
  func_profile_filt <- map(func_profile, function(x) filter_rownames(x, keep))
  return(func_profile_filt)
  
}
######################### ALPHA DIVERSITY #########################

# calculate shannon effective
Shannon.E <- function(x) {
  summed <- sum(x)
  shannon.e <-
    round(exp(-sum(x[x > 0] / summed * log(x[x > 0] / summed))), digits = 2)
  return(shannon.e)
}

# calculate richness
Richness <- function(x, detection = 1e-5) {
  observed <- sum(x > detection)
  return(observed)
}



# calculate all alpha diversity matrices and return dataframe
calc_alpha <- function(ps, ...) {
  mat_in <- ps_to_feattab(ps) %>%
    t()
  
  diversity <-
    setNames(data.frame(matrix(ncol = 2, nrow = nsamples(ps))),
             c("Richness", "Shannon.Effective"))
  rownames(diversity) <- rownames(meta_to_df(ps))
  
  diversity$Richness <- apply(mat_in, 1, Richness, ...)
  diversity$Shannon.Effective <- apply(mat_in, 1, Shannon.E)
  
  diversity_meta <- merge(diversity, meta_to_df(ps), by=0)
  
  return(diversity_meta)
}



######################### BETA DIVERSITY #########################

calc_betadiv <- function(ps, dist, ord_method = "NMDS") {
  if (ord_method %in% c("NMDS", "MDS", "PCoA")) {
    
    if (dist %in% c("bray", "jaccard")) {
      dist_mat <- distance(ps, dist)
      if (ord_method == "NMDS"){
        ord <- ordinate(ps, ord_method, dist_mat, trace=FALSE)
      } else {
        ord <- ordinate(ps, ord_method, dist_mat)
      }
      
      return_list <- list("Distance_Matrix" = dist_mat,
                          "Ordination" = ord)
      return(return_list)
    }
    else {
      print(
        "Distance metric not supported, supported metrics are; bray, jaccard"
      )
    }
  }
  else {
    print("Ordination method not supported, supported methods are: NMDS, MDS, PCoA")
    stop()
  }
  
}


######################### PLOTTING AND STATS #########################

# calculate colour palettes
calc_pal <- function(ps, group_variable) {
  # update function to accept nonps objects and add support for continuous palettes
  meta <- meta_to_df(ps)
  groups <- unique(meta[, group_variable])
  
  if (length(groups) < 10) {
    pal <- pal_npg()(length(groups))
  }
  else {
    print("Exceeded colour palette limit, use custom palette")
  }
}

# calculate adonis R2 and p-value
# better NA handling here
phyloseq_adonis <- function(ps, dist_matrix, group_variable, ...) {
  meta_df <- meta_to_df(ps)
  if (sum(is.na(meta_df[[group_variable]])) > 0) {
    print("metadata contains NAs, remove these samples with subset_samples
          before continuing")
    return (NULL)
  } else {
    # convert distance matix object to string
    dist_str <- deparse(substitute(dist_matrix))
    # define formula
    form <- as.formula(paste(dist_str, group_variable, sep="~"))
    ps_ad <- adonis2(form, data = meta_df, ...)
    return(ps_ad)
  }
}

# calculate betadisper
phyloseq_betadisper <- function(ps, dist_matrix, group_variable, ...) {
  meta_df <- meta_to_df(ps)
  if (sum(is.na(meta_df[[group_variable]])) > 0) {
    print("metadata contains NAs, remove these samples with subset_samples
          before continuing")
    return (NULL)
  } else {
    
    bd <- betadisper(dist_matrix, meta_df[[group_variable]])
    anova_res <- anova(bd)
    return(anova_res)
  }
}

# create statistical formula
xyform <- function (y_var, x_vars) {
  # y_var: a length-one character vector
  # x_vars: a character vector of object names
  as.formula(sprintf("%s ~ %s", y_var, paste(x_vars, collapse = " + ")))
}

round_any <- function(x, accuracy, f=ceiling) {
  f(x/ accuracy) * accuracy
}

# plot boxplot with stats
plot_boxplot <- function(df,
                         variable_col,
                         value_col,
                         comparisons_list,
                         fill_var = variable_col,
                         xlab = variable_col,
                         ylab = value_col,
                         p_title = NULL,
                         multiple_groups = FALSE,
                         cols = NULL,
                         group.order = NULL,
                         paired = FALSE,
                         ...) {
  # extend color palette with transparent value - required due to way we are
  # layering plot
  if (is.null(cols)) {
    cols <- pal_npg()(length(unique(df[, variable_col])))
  }
  cols <- c(cols, "transparent")
  
  if (!is.null(group.order)) {
    df[, variable_col] <-
      factor(df[, variable_col], levels = group.order)
  }
  
  formula <- xyform(value_col, variable_col)
  
  if (multiple_groups == TRUE) {
    if (paired == TRUE) {
      stat_variance <- df %>%
        friedman_test(formula)
      stat_test <- df %>%
        pairwise_wilcox_test(
          formula,
          comparisons = comparisons_list,
          p.adjust.method = "BH",
          paired = TRUE
        ) %>%
        add_significance() %>%
        add_xy_position(x = variable_col) 
    }
    else {
      stat_variance <- df %>%
        kruskal_test(formula)
      stat_test <- df %>%
        pairwise_wilcox_test(formula,
                             comparisons = comparisons_list,
                             p.adjust.method = "BH") %>%
        add_significance() %>%
        add_xy_position(x = variable_col) 
    }
  }
  else if (multiple_groups == FALSE) {
    if (paired == TRUE) {
      stat_test <- df %>%
        wilcox_test(formula, paired = TRUE) %>%
        add_significance() %>%
        add_xy_position(x = variable_col) 
    }
    else {
      stat_test <- df %>%
        wilcox_test(formula) %>%
        add_significance() %>%
        add_xy_position(x = variable_col) 
    }
    
  }
  print(stat_test)
  # aes string accepts strings as column names, this code plots boxplot and adds error bars
  plot <- ggplot(
    df,
    aes_string(
      x = variable_col,
      y = value_col,
      fill = variable_col,
      color = variable_col
    )
  ) +
    geom_boxplot(
      color = "black",
      alpha = 0.8,
      outlier.shape = 5,
      outlier.size = 1
    ) +
    geom_point(size = 1.5, position = position_jitterdodge()) +
    labs(x = xlab, y = ylab) +
    stat_boxplot(color = "black",
                 geom = "errorbar",
                 width = 0.2)
  # creates new 'finalised plot' and adds statistical significance, labels and adjusts theme and title
  final_plot <- plot +
    theme_classic2() +
    ggtitle(p_title) +
    theme(
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position = "None"
    ) +
    scale_fill_manual(values = cols) +
    scale_color_manual(values = cols) +
    rotate_x_text(angle = 45)
  
  if (dim(stat_test)[1] == 0) {
    plot_out <- final_plot
  }
  else {
    
    if (multiple_groups == T) {
      plot_out <- final_plot +
        stat_pvalue_manual(
          stat_test,
          label = "p.adj",
          inherit.aes = FALSE,
          bracket.nudge.y = (max(df[, value_col]) * 0.02), 
          ...
        )
    }
    else {
      plot_out <- final_plot +
        stat_pvalue_manual(
          stat_test,
          label = "p",
          inherit.aes = FALSE,
          bracket.nudge.y = (max(df[, value_col]) * 0.02),
          ...
        )
      
    }
  }
  
  return(plot_out)
}

plot_taxonomic_comp  <- function(ps, tax_level, var, ord=NULL, n_taxa=10, 
                                 per_group=F, groups=NULL, agg=F) {
  ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  if (!tax_level %in% ranks){
    stop("Provide a valid taxonomic rank")
  }
  var <- enquo(var)
  
  if (per_group == T){
    top_tax <- purrr::map(groups, function(x) ps %>% 
                     get_top_n_group(n = n_taxa, level = tax_level,var= !!var,
                                     group =x, agg=agg))
    top_tax <- c(unique(unlist(top_tax)), "other")
    # redefine n_taxa as in this case it reflects n across groups
    n_taxa <- length(top_tax)
  } else {
    top_tax <- c(get_top_n(ps, n=n_taxa, level = tax_level), "other") 
  }
  
  if (!is.null(ord)) {
    ps <- ps %>%
     ps_mutate(plot_var = factor(.data[[var]], levels = ord))
  }
  
  if (agg == T){
    melt_ps <- ps %>% 
      aggregate_taxa(level = tax_level) %>% 
      psmelt() 
  } else {
    melt_ps <- ps %>% 
      psmelt()
  }
  
  taxa_plot <- melt_ps %>% filter(OTU %in% top_tax)
  
  # group rest into other
  repl_cols <- c("OTU", ranks)
  
  other_plot <- melt_ps %>% filter(!OTU %in% top_tax) %>% 
    group_by(Sample) %>% 
    mutate(Abundance = sum(Abundance)) %>% 
    ungroup %>% 
    distinct(Sample, .keep_all = T) %>% 
    mutate_at(vars(repl_cols), ~ paste("other"))
  
  # bind dfs
  plot_df <- bind_rows(taxa_plot, other_plot)
  
  taxa_pal <- c(stacked_bar.palette[1:length(top_tax)-1], "#DCDCDC", "white")
  names(taxa_pal) <- plot_df %>% filter(OTU %in% top_tax) %>% 
    pull(.data[[ tax_level ]]) %>% 
    unique()
  
  comp_fig <- plot_df %>%
    ggplot(aes(fill=.data[[ tax_level ]], y=Abundance, x=Sample)) + 
    geom_bar(position="fill", stat="identity") + 
    scale_fill_manual(values = taxa_pal) + 
    facet_grid(. ~ factor(plot_var),
               scales = "free", space = "free") +
    theme_cowplot() +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    theme(
      panel.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 22),
      axis.text.y = element_text(size = 18),
      axis.ticks = element_blank()
    )
  return(comp_fig)
}

plot_beta_div <- function(ps,
                          beta_div,
                          group_variable,
                          add_ellipse = FALSE,
                          cols = NULL,
                          shape = NULL) {
  
  
  if (is.null(cols)) {
    cols <- calc_pal(ps, group_variable)
  }
  
  ordination <- beta_div$Ordination
  dist_matrix <- beta_div$Distance_Matrix#
  
  # significance
  ad <- phyloseq_adonis(ps, dist_matrix, group_variable)
  betadisp <- phyloseq_betadisper(ps, dist_matrix, group_variable)
  
  # check sample data dimensions >1 otherwise phyloseq
  # fails to colour samples due to bug:https://github.com/joey711/phyloseq/issues/541
  if (dim(sample_data(ps))[2] < 2) {
    # add repeat of first column as dummy
    sample_data(ps)[, 2] <- sample_data(ps)[, 1]
  }
  plot <- plot_ordination(ps, ordination , color = group_variable,  shape = shape)
  plot$layers[[1]] <- NULL
  
  plot_out <- plot + geom_point(size = 3, alpha = 0.75) +
    theme_bw() +
    scale_fill_manual(values = cols) +
    scale_color_manual(values = cols) +
    labs(caption = bquote(Adonis ~ R ^ 2 ~ .(round(ad$R2[1], 2)) ~
                            ~ p - value ~ .(ad$`Pr(>F)`[1])))
  if (add_ellipse == TRUE){
    plot_out <- plot_out +
      geom_polygon(stat = "ellipse", aes(fill = .data [[ group_variable ]] ), alpha = 0.3)
  }
  
  # ideally this needs to be added to the main plot eventually underneath adonis
  if (betadisp$`Pr(>F)`[[1]] < 0.05){
    warning("Group dispersion is not homogenous, interpret results carefully",
            call. = F)
  }
  
  return(plot_out)
}

maaslin2_tax <- function(tax_profile, out, 
                         fixed, abun_thresh=0, prev_thresh=0, taxonomy="metaphlan",...){
  
  mat_in <- ps_to_feattab(tax_profile) %>% 
    t_df()
  
  # this only works for species level
  if (taxonomy == "metaphlan"){
    colnames(mat_in) <- sapply(strsplit(colnames(mat_in), split= "\\|", fixed = TRUE), tail, 1L)
  } else if (taxonomy == "gtdb"){
    colnames(mat_in) <- sapply(strsplit(colnames(mat_in), split= ";", fixed = TRUE), tail, 1L)
  } else {
    stop("unsupported taxonomy")
  }
  
  
  metadata_in <- meta_to_df(tax_profile)
  
  Maaslin2(input_data = mat_in,
           input_metadata = metadata_in, 
           output = out, 
           normalization = "NONE",
           transform = "NONE",
           min_abundance = abun_thresh,
           min_prevalence = prev_thresh,
           fixed_effects = fixed,
           cores=6,
           ...
  )
}

plot_da <- function(maaslin_res, groups, cols = NULL, pthresh=0.25) {
  if (is.null(cols)) {
    cols <- pal_npg()(length(groups))
  }
  
  da_plot <- maaslin_res %>%
    mutate(enriched_in = ifelse(coef > 0, groups[2],
                                groups[1])) %>% 
    # species name
    mutate(species = str_extract(feature, "(?<=s__).*")) %>% 
    filter(qval <= pthresh)
  
  da_plot$species <- gsub("\\.", " ", da_plot$species)
  
  da_plot_sort <- da_plot %>%
    arrange(coef) %>%
    mutate(species = factor(species, levels = species))
  
  
  p <-
    ggplot(da_plot_sort,
           aes(x = species, y = coef, color = enriched_in)) +
    geom_point(size = 5) +
    labs(y = paste("coeffiecient", groups[1], "vs", groups[2], sep =
                     " "),
         x = "") +
    theme_bw() +
    theme(
      axis.text.x = element_text(color = "black", size = 14),
      axis.text.y = element_text(color = "black", size = 14),
      axis.title.y = element_text(size = 16),
      axis.title.x = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14),
      legend.position = "none"
    ) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dotted") +
    scale_color_manual(values = cols) 
  print(p)
}
