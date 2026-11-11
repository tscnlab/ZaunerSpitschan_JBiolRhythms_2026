# H02 sequence logic reproduced locally so upstream files remain read only.
assign_sequences <- function(x) {
  x <- x |> arrange(site, Id, local_date, source_utc_start, clock_bin) |>
    group_by(site, Id, local_date) |>
    mutate(elapsed_from_previous_seconds = as.numeric(difftime(source_utc_start, lag(source_utc_end), units = 'secs')),
      ar_start_reason = case_when(row_number() == 1L ~ 'participant_day_start',
        !one_to_one_elapsed_coordinate ~ 'non_one_to_one_wall_outcome',
        !coalesce(lag(one_to_one_elapsed_coordinate), TRUE) ~ 'after_non_one_to_one_wall_outcome',
        is.na(elapsed_from_previous_seconds) ~ 'missing_elapsed_coordinate',
        elapsed_from_previous_seconds != 0 ~ 'elapsed_discontinuity', TRUE ~ 'continuous'),
      AR_start = ar_start_reason != 'continuous', seq_no = cumsum(AR_start)) |>
    ungroup() |> mutate(true_elapsed_sequence_id = paste(participant_day_key, seq_no, sep = '::')) |> select(-seq_no)
  stopifnot(all(x$AR_start[!duplicated(x$participant_day_key)]), all(x$elapsed_from_previous_seconds[!x$AR_start] == 0), !anyNA(x$AR_start))
  x
}
quartile_cutoffs <- function(x) {
  p <- distinct(x, participant_key, msf_hours)
  stopifnot(!anyDuplicated(p$participant_key), all(is.finite(p$msf_hours)))
  unname(quantile(p$msf_hours, c(.25,.5,.75), type = 7))
}
apply_quartiles <- function(x, cutoffs) {
  x |> mutate(quartile = factor(cut(msf_hours, c(-Inf, cutoffs, Inf), labels = paste0('Q',1:4), right = TRUE), levels = paste0('Q',1:4)),
    direction = factor(if_else(quartile %in% c('Q1','Q2'), 'Morning','Evening'), levels = c('Morning','Evening')),
    extremity = factor(if_else(quartile %in% c('Q2','Q3'), 'Moderate','Extreme'), levels = c('Moderate','Extreme')),
    group_label = paste(extremity, tolower(as.character(direction))))
}
input_manifest <- function() {
  manifest <- read.csv(file.path(project_root,'input-data/manifest.csv'),check.names=FALSE)
  manifest$snapshot_path <- manifest$path
  manifest$path <- file.path(project_root,manifest$snapshot_path)
  stopifnot(!anyDuplicated(manifest$role),all(file.exists(manifest$path)),
    identical(unname(tools::md5sum(manifest$path)),manifest$md5))
  as_tibble(manifest)
}
input_paths <- function() {
  manifest <- input_manifest()
  setNames(manifest$path,manifest$role)
}
prepare_data <- function() {
  manifest <- input_manifest()
  paths <- setNames(manifest$path,manifest$role)
  save_table(manifest, 'inputs', 'results/provenance')
  ch <- readRDS(paths[['chronotype']]) |> transmute(site, Id, msf_hours = as.numeric(msf_sc)/3600, meq, meq_type = as.character(meq_type))
  diary <- readRDS(paths[['sleepdiaries']]) |> transmute(site, Id, local_date = as.Date(wake_wall), diary_present = TRUE,
    daytype = factor(as.character(daytype2), levels = c('a work day','a free day'), labels = c('Work','Free')))
  stopifnot(!anyDuplicated(ch[c('site','Id')]), !anyDuplicated(diary[c('site','Id','local_date')]))
  frames <- flows <- exclusions <- list()
  for (placement in c('glasses','chest')) {
    h <- readRDS(paths[[placement]])
    context <- readRDS(paths[[paste0(placement,'_context')]]) |>
      select(site,Id,position,local_date,clock_bin,photoperiod_hours)
    key <- c('site','Id','position','local_date','clock_bin')
    stopifnot(!anyDuplicated(h[key]),!anyDuplicated(context[key]),all(h$valid_medi_wall_minutes >= 15),all(h$bin_admissible),
      all(h$metric_value_lx >= 0),isTRUE(all.equal(h$response,log10(h$metric_value_lx+.1))))
    x <- h |> left_join(ch, by = c('site','Id'), relationship = 'many-to-one') |>
      left_join(diary, by = c('site','Id','local_date'), relationship = 'many-to-one') |>
      left_join(context, by = key, relationship = 'one-to-one')
    stopifnot(nrow(x) == nrow(h))
    a <- filter(x, is.finite(msf_hours)); b <- filter(a, !is.na(daytype))
    stopifnot(all(is.finite(b$photoperiod_hours)))
    flows[[placement]] <- bind_rows(mutate(count_sample(h),stage = 'Prepared sensor sample'),
      mutate(count_sample(a),stage = 'MSFsc available'),mutate(count_sample(b),stage = 'MSFsc and day type available')) |> mutate(placement = placement, .before = 1)
    exclusions[[placement]] <- x |> mutate(exclusion = case_when(!is.finite(msf_hours) ~ 'Missing MSFsc',
      is.na(diary_present) ~ 'No diary for waking date',is.na(daytype) ~ 'Missing work/free label',TRUE ~ 'Included')) |>
      distinct(participant_key,participant_day_key,site,local_date,exclusion) |> mutate(placement = placement)
    frames[[placement]] <- b |> mutate(clock = (clock_bin+15)/60,phase = wrap24(clock-msf_hours),phase_display = wrap12(clock-msf_hours)) |> assign_sequences()
  }
  cutoffs <- quartile_cutoffs(frames$glasses)
  frames <- lapply(frames, apply_quartiles, cutoffs = cutoffs)
  registry <- bind_rows(frames, .id = 'placement') |>
    select(placement,site,Id,participant_key,participant_day_key,local_date,clock_bin,clock,phase,phase_display,quartile,direction,extremity,group_label,daytype,msf_hours,photoperiod_hours,valid_medi_wall_minutes)
  save_table(bind_rows(flows), 'sample_flow')
  save_table(bind_rows(exclusions), 'day_inclusion_registry')
  save_table(bind_cols(tibble(percentile=c(25,50,75),msf_hours=cutoffs,clock_time=clock_label(cutoffs)),count_sample(frames$glasses)[rep(1,3),]), 'quartile_cutoffs')
  p <- distinct(frames$glasses,participant_key,quartile,direction,extremity,msf_hours,meq,meq_type)
  save_table(frames$glasses |> summarise(participants=n_distinct(participant_key),participant_days=n_distinct(participant_day_key),half_hours=n(),.by=c(quartile,meq_type)), 'quartile_meq_overlap')
  qcounts<-frames$glasses |> summarise(participant_days=n_distinct(participant_day_key),half_hours=n(),.by=quartile)
  save_table(p |> summarise(participants=n(),minimum=min(msf_hours),median=median(msf_hours),maximum=max(msf_hours),.by=c(quartile,direction,extremity)) |> left_join(qcounts,by='quartile'), 'quartile_distribution')
  save_table(registry |> summarise(participants=n_distinct(participant_key),participant_days=n_distinct(participant_day_key),half_hours=n(),.by=c(placement,quartile,daytype)), 'group_day_counts')
  save_table(registry, 'sample_registry', 'results/data')
  out <- list(frames=frames,registry=registry,cutoffs=cutoffs,manifest=manifest,flow=bind_rows(flows),participants=p)
  atomic_rds(out,'results/data/prepared.rds')
  versions <- tibble(package=c('R','mgcv','gratia','dplyr','tidyr','ggplot2','gt','knitr','renv','LightLogR','scales'),
    version=c(as.character(getRversion()),vapply(c('mgcv','gratia','dplyr','tidyr','ggplot2','gt','knitr','renv','LightLogR','scales'),function(p)as.character(packageVersion(p)),character(1))))
  save_table(versions,'versions','results/provenance')
  writeLines(c(capture.output(sessionInfo()),'',paste('OpenMP preload:',omp),'Command: Rscript --vanilla run.R [stage]','Quarto: R_PROFILE_USER=.Rprofile quarto render'), 'results/provenance/session-info.txt')
  out
}
fit_preprocess <- function(d, cutoffs = quartile_cutoffs(d)) {
  centers <- d |> distinct(site,local_date,photoperiod_hours) |> summarise(photo_center=mean(photoperiod_hours),.by=site)
  list(cutoffs=cutoffs, centers=centers)
}
apply_preprocess <- function(d, prep, axis) {
  d <- d |> select(-any_of(c('photo_center','photo'))) |> mutate(site=as.character(site)) |>
    left_join(prep$centers |> mutate(site=as.character(site)),by='site',relationship='many-to-one') |>
    mutate(photo=photoperiod_hours-photo_center, t=.data[[axis]]) |> apply_quartiles(prep$cutoffs) |> assign_sequences()
  stopifnot(all(is.finite(d$photo)))
  d |> mutate(site=factor(site),participant=factor(participant_key),participant_day=factor(participant_day_key))
}
