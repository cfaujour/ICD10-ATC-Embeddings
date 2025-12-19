# toy_data/01_dataset_toy_gen.R
# Toy dataset generator for co-occurrence / embedding demos.
# Output: toy_data/data/sequences.csv.gz and toy_data/data/code_topics.csv

suppressPackageStartupMessages({
  library(data.table)
})

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# ---- Output folder (relative to repo root) ----
out_dir <- file.path("toy_data", "data")
ensure_dir(out_dir)

# ------------------------- Concept vocabularies -------------------------

concepts <- list(
  diabetes = list(
    icd = c("ICD-E11","ICD-E119","ICD-E1165","ICD-E114","ICD-E118","ICD-E10","ICD-E109"),
    atc = c("ATC-A10BA02","ATC-A10BD02","ATC-A10BB12","ATC-A10BG03","ATC-A10BH01","ATC-A10BJ01")
  ),
  hypertension = list(
    icd = c("ICD-I10","ICD-I119","ICD-I129","ICD-I1310","ICD-I150","ICD-I152"),
    atc = c("ATC-C09AA02","ATC-C09CA01","ATC-C07AB02","ATC-C08CA01","ATC-C03CA01")
  ),
  airway = list(
    icd = c("ICD-J45","ICD-J459","ICD-J440","ICD-J441","ICD-J449"),
    atc = c("ATC-R03AC02","ATC-R03AK06","ATC-R03BA02","ATC-R03BB04","ATC-R03DC02")
  ),
  depression = list(
    icd = c("ICD-F320","ICD-F329","ICD-F330","ICD-F331","ICD-F339"),
    atc = c("ATC-N06AB03","ATC-N06AB10","ATC-N06AB04","ATC-N06AX16","ATC-N06AX11")
  ),
  ischemic_hd = list(
    icd = c("ICD-I209","ICD-I214","ICD-I2510","ICD-I252","ICD-I259"),
    atc = c("ATC-C01DA14","ATC-C07AB02","ATC-B01AC06","ATC-C01EB15","ATC-C10AA05")
  ),
  infections = list(
    icd = c("ICD-J029","ICD-J189","ICD-N390","ICD-J00","ICD-J069"),
    atc = c("ATC-J01CA04","ATC-J01CR02","ATC-J01FA10","ATC-J01MA12","ATC-J01XE01")
  )
)

background_icd <- c("ICD-Z000","ICD-R05","ICD-M545","ICD-K219","ICD-Z139","ICD-R519")
background_atc <- c("ATC-A02BC02","ATC-N02BE01","ATC-R05DA04","ATC-A03FA01")

# ------------------------- Generator -------------------------

draw_patient_profile <- function() {
  K <- sample(2:3, 1, prob = c(0.65, 0.35))
  picks <- sample(names(concepts), K, replace = FALSE)
  w <- rep(0, length(concepts)); names(w) <- names(concepts)
  w[picks] <- runif(K, 0.4, 0.8)
  bg <- max(0.05, 1 - sum(w))
  list(weights = w, bg = bg)
}

emit_visit_codes <- function(concept_name, t0,
                             base_icd_n = 2:3, base_atc_n = 1:3,
                             jitter_prob = 0.25) {
  if (concept_name == "background") {
    icd <- sample(background_icd, sample(0:2, 1), replace = TRUE)
    atc <- sample(background_atc, sample(0:2, 1), replace = TRUE)
  } else {
    icd_pool <- concepts[[concept_name]]$icd
    atc_pool <- concepts[[concept_name]]$atc
    icd <- sample(icd_pool, sample(base_icd_n, 1), replace = TRUE)
    atc <- sample(atc_pool, sample(base_atc_n, 1), replace = TRUE)
    if (runif(1) < 0.20) icd <- c(icd, sample(background_icd, 1))
    if (runif(1) < 0.20) atc <- c(atc, sample(background_atc, 1))
  }
  
  codes_vec <- c(icd, atc)
  onto_vec  <- substr(codes_vec, 1, 3)
  
  codes <- data.table(Event_code = codes_vec, onto = onto_vec, t = t0)
  
  if (runif(1) < jitter_prob && length(codes_vec) > 0) {
    jitter_t <- t0 + sample(c(-1L, 1L), 1)
    pick_n <- sample(1:min(2L, length(codes_vec)), 1)
    pick_idx <- sample(seq_along(codes_vec), pick_n)
    jcodes <- data.table(Event_code = codes_vec[pick_idx], onto = onto_vec[pick_idx], t = jitter_t)
    codes <- rbind(codes, jcodes, use.names = TRUE)
  }
  
  codes[]
}

generate_toy_sequences <- function(n_patients = 1500, mean_visits = 12, seed = 42) {
  set.seed(seed)
  
  code_topics <- unique(rbindlist(list(
    data.table(Event_code = unlist(lapply(concepts, `[[`, "icd")),
               onto = "ICD",
               topic = rep(names(concepts), lengths(lapply(concepts, function(x) x$icd)))),
    data.table(Event_code = unlist(lapply(concepts, `[[`, "atc")),
               onto = "ATC",
               topic = rep(names(concepts), lengths(lapply(concepts, function(x) x$atc)))),
    data.table(Event_code = background_icd, onto = "ICD", topic = "background"),
    data.table(Event_code = background_atc, onto = "ATC", topic = "background")
  ), use.names = TRUE))
  
  out_list <- vector("list", n_patients)
  
  for (pid in seq_len(n_patients)) {
    prof <- draw_patient_profile()
    visits <- max(1L, rpois(1, lambda = mean_visits))
    
    start_age_days <- sample(18:75, 1) * 365L
    gaps <- pmax(1L, as.integer(rexp(visits, rate = 1/90)))
    t_seq <- start_age_days + cumsum(gaps)
    
    rows <- vector("list", visits)
    
    for (v in seq_len(visits)) {
      t0 <- t_seq[v]
      
      w <- pmax(prof$weights, 1e-9)
      focus_n <- sample(1:2, 1, prob = c(0.65, 0.35))
      focus <- sample(names(concepts), focus_n, prob = w)
      
      bundle <- rbindlist(lapply(focus, function(cname) emit_visit_codes(cname, t0)), use.names = TRUE)
      
      if (runif(1) < prof$bg) {
        bundle <- rbind(bundle, emit_visit_codes("background", t0), use.names = TRUE)
      }
      
      rows[[v]] <- data.table(
        ID = sprintf("P%05d", pid),
        Event_code = bundle$Event_code,
        t = bundle$t,
        onto = bundle$onto
      )
    }
    
    out_list[[pid]] <- rbindlist(rows, use.names = TRUE)
  }
  
  events <- rbindlist(out_list, use.names = TRUE)
  setkeyv(events, c("ID", "t", "Event_code"))
  
  list(events = events[], code_topics = code_topics[])
}

# ------------------------- Run -------------------------

toy <- generate_toy_sequences(n_patients = 1500, mean_visits = 12, seed = 42)

fwrite(toy$events, file.path(out_dir, "sequences.csv.gz"))
fwrite(toy$code_topics, file.path(out_dir, "code_topics.csv"))

message("Toy dataset written to: ", out_dir)
