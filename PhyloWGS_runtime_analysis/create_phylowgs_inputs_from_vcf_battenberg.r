# ============================================================================
# Create PhyloWGS inputs (ssm_data.txt, cnv_data.txt) from VCF and Battenberg
# ============================================================================
#
# Reads:
#   - MuTect VCF (e.g. {sample}.mutect.vcf) via VariantAnnotation
#   - Battenberg copy-number file (e.g. {sample}.battenberg.txt)
#
# Writes:
#   - PhyloWGS_runtime_analysis/mutation_table/{sample}_ssm_data.txt
#   - PhyloWGS_runtime_analysis/copy_number_file/{sample}_cnv_data.txt
#
# Usage:
#   - Set SAMPLE_DIR to a directory with subdirs like P3-noXY containing
#     P3-noXY.mutect.vcf and P3-noXY.battenberg.txt; then run for all samples.
#   - Or set VCF_FILE and BATTENBERG_FILE (and SAMPLE_ID) for a single run.
#
# ============================================================================

suppressPackageStartupMessages({
    if (!requireNamespace("VariantAnnotation", quietly = TRUE)) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
        BiocManager::install("VariantAnnotation")
    }
    library(VariantAnnotation)
})

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------

# Root of PhyloWGS repo (script lives in PhyloWGS_runtime_analysis/)
PHYLOWGS_ROOT <- "/Users/khanhngocdinh/Documents/YanjieChen/GitHub/PhyloWGS"
setwd(PHYLOWGS_ROOT)

# DREAM-style data: one directory per sample with {sample}.mutect.vcf and {sample}.battenberg.txt
SAMPLE_DIR <- "/Volumes/My Passport for Mac/DECODE_project/data/DREAM_29-2026-02-22"
# Example single sample (override for one-off run):
# SAMPLE_DIR <- file.path(PHYLOWGS_ROOT, "data", "DREAM_29-2026-02-22", "P3-noXY")

# Output directories (under PhyloWGS)
OUT_MUTATION_TABLE <- file.path(PHYLOWGS_ROOT, "PhyloWGS_runtime_analysis", "mutation_table")
OUT_COPY_NUMBER    <- file.path(PHYLOWGS_ROOT, "PhyloWGS_runtime_analysis", "copy_number_file")

# Optional: single-file mode (set to non-NULL to use)
# Example for P3-noXY (uncomment and set paths to run one sample):
# VCF_FILE       <- file.path(PHYLOWGS_ROOT, "data", "DREAM_29-2026-02-22", "P3-noXY", "P3-noXY.mutect.vcf")
# BATTENBERG_FILE <- file.path(PHYLOWGS_ROOT, "data", "DREAM_29-2026-02-22", "P3-noXY", "P3-noXY.battenberg.txt")
# SAMPLE_ID      <- "P3-noXY"
VCF_FILE       <- NULL
BATTENBERG_FILE <- NULL
SAMPLE_ID      <- NULL

# PhyloWGS constants (README)
MU_R <- 0.999   # reference-allele fraction in reference population
MU_V <- 0.5     # reference-allele fraction in variant population (diploid)
READ_LENGTH    <- 150
DEFAULT_DP_CNV <- 500    # fallback mean depth for CNV a/d when not estimated from VCF

# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------

dir_create_if_not <- function(path) {
    if (!dir.exists(path)) {
        dir.create(path, recursive = TRUE, showWarnings = FALSE)
    }
}

# Normalize chromosome for matching (e.g. "chr1" -> "1", "X" -> "x")
norm_chr <- function(x) {
    x <- as.character(x)
    x <- sub("^chr", "", x, ignore.case = TRUE)
    x <- tolower(x)
    x
}

# ----------------------------------------------------------------------------
# VCF -> ssm_data.txt (PhyloWGS format)
# ----------------------------------------------------------------------------
# Columns: id, gene, a, d, mu_r, mu_v
# id: s0, s1, ... ; gene: e.g. chr_pos; a = ref reads, d = total (ref+alt or DP)
# ----------------------------------------------------------------------------
create_ssm_data_from_vcf <- function(vcf_file, sample_id, out_path) {
    if (!file.exists(vcf_file)) stop("VCF file not found: ", vcf_file)

    vcf <- readVcf(vcf_file)
    n <- nrow(vcf)

    if (n == 0) {
        message("No variants in VCF; writing empty ssm_data header only.")
        dir_create_if_not(dirname(out_path))
        writeLines("id\tgene\ta\td\tmu_r\tmu_v", out_path)
        return(invisible(data.frame()))
    }

    chr  <- as.character(seqnames(vcf))
    pos  <- start(ranges(vcf))
    gene <- paste0(norm_chr(chr), "_", pos)

    # Tumor sample: use last sample if multiple (MuTect often normal,tumor)
    snames <- colnames(geno(vcf)$GT)
    if (length(snames) == 0) stop("No genotype columns in VCF")
    tumor_col <- if (length(snames) >= 2 && any(tolower(snames) == "tumor")) {
        which(tolower(snames) == "tumor")[1]
    } else {
        length(snames)
    }

    ad <- geno(vcf)$AD
    dp <- geno(vcf)$DP
    if (is.null(ad)) stop("VCF has no AD (allelic depth) in FORMAT")

    # AD can be IntegerList per variant (ref, alt)
    ref <- integer(n)
    tot <- integer(n)
    for (i in seq_len(n)) {
        ad_i <- ad[i, tumor_col][[1]]
        if (is.null(ad_i) || length(ad_i) < 2) {
            ref[i] <- 0L
            tot[i] <- if (!is.null(dp)) as.integer(dp[i, tumor_col]) else 0L
            if (is.na(tot[i])) tot[i] <- 0L
        } else {
            ref[i] <- as.integer(ad_i[1])
            tot[i] <- as.integer(sum(ad_i, na.rm = TRUE))
        }
        if (tot[i] == 0 && !is.null(dp)) tot[i] <- as.integer(dp[i, tumor_col])
        if (is.na(tot[i]) || tot[i] < 1) tot[i] <- 1L
    }

    id <- paste0("s", seq(0, n - 1))
    a  <- as.character(ref)
    d  <- as.character(tot)

    out <- data.frame(
        id   = id,
        gene = gene,
        a    = a,
        d    = d,
        mu_r = MU_R,
        mu_v = MU_V,
        stringsAsFactors = FALSE
    )

    dir_create_if_not(dirname(out_path))
    write.table(out, out_path, sep = "\t", row.names = FALSE, quote = FALSE)
    message("Wrote ", n, " SSMs to ", out_path)

    # Return with chr/pos for CNV overlap later
    out$chromosome <- chr
    out$position   <- pos
    invisible(out)
}

# ----------------------------------------------------------------------------
# Battenberg file parsing (support DECODE-style and PhyloWGS-style)
# ----------------------------------------------------------------------------
# DECODE: tab-delimited, header chr, startpos, endpos, frac1_A, nMaj1_A, nMin1_A (and possibly subclonal cols)
# PhyloWGS Python: space/tab, no header, chrom=fields[2], start=fields[3], end=fields[4], major=fields[9], minor=fields[10], frac=fields[11]
# ----------------------------------------------------------------------------
read_battenberg <- function(battenberg_file, cellularity = 1.0) {
    if (!file.exists(battenberg_file)) stop("Battenberg file not found: ", battenberg_file)

    raw <- readLines(battenberg_file)
    raw <- raw[raw != ""]
    if (length(raw) < 2) {
        return(data.frame(
            chr = character(0), start = integer(0), end = integer(0),
            major_cn = integer(0), minor_cn = integer(0), cell_prev = numeric(0)
        ))
    }

    first <- read.table(text = raw[1], header = FALSE, nrows = 1, stringsAsFactors = FALSE)
    has_header <- length(first) >= 6 && any(grepl("chr|start|nMaj|nMin|frac", raw[1], ignore.case = TRUE))

    if (has_header) {
        bb <- read.table(battenberg_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        # DECODE-style: chr, startpos, endpos, frac1_A, nMaj1_A, nMin1_A
        chr_col  <- which(tolower(names(bb)) %in% c("chr", "chromosome", "chrom"))[1]
        start_col <- which(tolower(names(bb)) %in% c("startpos", "start"))[1]
        end_col   <- which(tolower(names(bb)) %in% c("endpos", "end"))[1]
        maj_col  <- which(grepl("nmaj|major", names(bb), ignore.case = TRUE))[1]
        min_col  <- which(grepl("nmin|minor", names(bb), ignore.case = TRUE))[1]
        frac_col <- which(grepl("frac|cellular|cell_prev|prev", names(bb), ignore.case = TRUE))[1]

        if (is.na(chr_col))  chr_col  <- 1
        if (is.na(start_col)) start_col <- 2
        if (is.na(end_col))   end_col   <- 3
        if (is.na(maj_col))  maj_col  <- which(tolower(names(bb)) == "nmaj1_a")[1]
        if (is.na(min_col))  min_col  <- which(tolower(names(bb)) == "nmin1_a")[1]
        if (is.na(frac_col)) frac_col <- which(tolower(names(bb)) == "frac1_a")[1]
        if (is.na(maj_col))  maj_col  <- 8L
        if (is.na(min_col))  min_col  <- 9L
        if (is.na(frac_col)) frac_col <- 10L

        chr   <- as.character(bb[[chr_col]])
        start <- as.integer(bb[[start_col]])
        end   <- as.integer(bb[[end_col]])
        major <- as.integer(round(as.numeric(bb[[maj_col]])))
        minor <- as.integer(round(as.numeric(bb[[min_col]])))
        if (!is.na(frac_col)) {
            cell_prev <- as.numeric(bb[[frac_col]]) * cellularity
        } else {
            cell_prev <- rep(cellularity, length(chr))
        }
    } else {
        # PhyloWGS-style: space/tab, no header; indices 2=chrom, 3=start, 4=end, 9=major, 10=minor, 11=frac
        bb <- read.table(battenberg_file, header = FALSE, sep = "", stringsAsFactors = FALSE)
        if (ncol(bb) < 11) {
            return(data.frame(
                chr = character(0), start = integer(0), end = integer(0),
                major_cn = integer(0), minor_cn = integer(0), cell_prev = numeric(0)
            ))
        }
        chr   <- as.character(bb[[2]])
        start <- as.integer(bb[[3]])
        end   <- as.integer(bb[[4]])
        major <- as.integer(bb[[9]])
        minor <- as.integer(bb[[10]])
        cell_prev <- as.numeric(bb[[11]]) * cellularity
    }

    data.frame(
        chr = norm_chr(chr),
        start = start,
        end = end,
        major_cn = major,
        minor_cn = minor,
        cell_prev = cell_prev,
        stringsAsFactors = FALSE
    )
}

# ----------------------------------------------------------------------------
# Build cnv_data.txt from Battenberg + SSM overlap
# ----------------------------------------------------------------------------
# PhyloWGS cnv_data: cnv, a, d, ssms, physical_cnvs
# ssms: "s2,1,2;s4,0,1" (SSM id, minor_cn, major_cn for that segment)
# physical_cnvs: "chrom=1,start=1234,end=5678,major_cn=2,minor_cn=1,cell_prev=0.8;..."
# ----------------------------------------------------------------------------
create_cnv_data_from_battenberg <- function(battenberg_file, ssm_df, sample_id, out_path,
                                            mean_dp_from_vcf = NULL, cellularity = 1.0) {
    segs <- read_battenberg(battenberg_file, cellularity = cellularity)
    dir_create_if_not(dirname(out_path))

    if (nrow(segs) == 0) {
        writeLines("cnv\ta\td\tssms\tphysical_cnvs", out_path)
        message("No Battenberg segments; wrote empty cnv_data to ", out_path)
        return(invisible(NULL))
    }

    # Estimate d (total reads) per segment if we have mean DP
    mean_dp <- mean_dp_from_vcf
    if (is.null(mean_dp) || !is.finite(mean_dp) || mean_dp < 1) mean_dp <- DEFAULT_DP_CNV
    segs$d_est <- pmax(1, round((segs$end - segs$start + 1) / READ_LENGTH * mean_dp))
    segs$a_est <- round(segs$d_est * segs$minor_cn / pmax(1, segs$major_cn + segs$minor_cn))

    # Assign SSMs to segments (by chr and position overlap)
    segs$ssms_str <- ""
    if (!is.null(ssm_df) && nrow(ssm_df) > 0 && "chromosome" %in% names(ssm_df)) {
        ssm_chr <- norm_chr(ssm_df$chromosome)
        ssm_pos <- ssm_df$position
        ssm_id  <- ssm_df$id
        for (i in seq_len(nrow(segs))) {
            in_seg <- which(ssm_chr == segs$chr[i] &
                            ssm_pos >= segs$start[i] &
                            ssm_pos <= segs$end[i])
            if (length(in_seg) > 0) {
                parts <- paste0(ssm_id[in_seg], ",",
                                segs$minor_cn[i], ",",
                                segs$major_cn[i])
                segs$ssms_str[i] <- paste(parts, collapse = ";")
            }
        }
    }

    # physical_cnvs string per row (one segment per row here)
    segs$physical_cnvs <- paste0(
        "chrom=", segs$chr, ",start=", segs$start, ",end=", segs$end,
        ",major_cn=", segs$major_cn, ",minor_cn=", segs$minor_cn,
        ",cell_prev=", round(segs$cell_prev, 6)
    )

    cnv_id <- paste0("c", seq(0, nrow(segs) - 1))
    out <- data.frame(
        cnv = cnv_id,
        a   = as.character(segs$a_est),
        d   = as.character(segs$d_est),
        ssms = segs$ssms_str,
        physical_cnvs = segs$physical_cnvs,
        stringsAsFactors = FALSE
    )

    write.table(out, out_path, sep = "\t", row.names = FALSE, quote = FALSE)
    message("Wrote ", nrow(segs), " CNV segments to ", out_path)
    invisible(out)
}

# ----------------------------------------------------------------------------
# Process one sample
# ----------------------------------------------------------------------------
process_one_sample <- function(sample_id, vcf_file, battenberg_file,
                               out_ssm_path, out_cnv_path,
                               cellularity = 1.0) {
    message("Processing sample: ", sample_id)

    ssm_df <- create_ssm_data_from_vcf(vcf_file, sample_id, out_ssm_path)
    mean_dp <- NULL
    if (nrow(ssm_df) > 0 && "d" %in% names(ssm_df)) {
        d_vals <- as.numeric(ssm_df$d)
        d_vals <- d_vals[is.finite(d_vals) & d_vals > 0]
        if (length(d_vals) > 0) mean_dp <- mean(d_vals)
    }
    create_cnv_data_from_battenberg(
        battenberg_file, ssm_df, sample_id, out_cnv_path,
        mean_dp_from_vcf = mean_dp, cellularity = cellularity
    )
    invisible(list(ssm = ssm_df))
}

# ----------------------------------------------------------------------------
# Discover samples (DREAM-style: subdirs with {name}.mutect.vcf and {name}.battenberg.txt)
# ----------------------------------------------------------------------------
get_samples_from_dir <- function(sample_dir) {
    if (!dir.exists(sample_dir)) return(character(0))
    subdirs <- list.dirs(sample_dir, full.names = FALSE, recursive = FALSE)
    subdirs <- subdirs[!subdirs %in% c("", ".", "..")]
    samples <- character(0)
    for (s in subdirs) {
        vcf <- file.path(sample_dir, s, paste0(s, ".mutect.vcf"))
        bat <- file.path(sample_dir, s, paste0(s, ".battenberg.txt"))
        if (file.exists(vcf) && file.exists(bat)) samples <- c(samples, s)
    }
    samples
}

# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
main <- function() {
    dir_create_if_not(OUT_MUTATION_TABLE)
    dir_create_if_not(OUT_COPY_NUMBER)

    if (!is.null(VCF_FILE) && !is.null(BATTENBERG_FILE)) {
        sid <- if (!is.null(SAMPLE_ID)) SAMPLE_ID else basename(dirname(VCF_FILE))
        process_one_sample(
            sample_id = sid,
            vcf_file = VCF_FILE,
            battenberg_file = BATTENBERG_FILE,
            out_ssm_path = file.path(OUT_MUTATION_TABLE, paste0(sid, "_ssm_data.txt")),
            out_cnv_path = file.path(OUT_COPY_NUMBER, paste0(sid, "_cnv_data.txt"))
        )
        return(invisible(NULL))
    }

    samples <- get_samples_from_dir(SAMPLE_DIR)
    if (length(samples) == 0) {
        stop("No samples found in ", SAMPLE_DIR,
             " (expect subdirs with {name}.mutect.vcf and {name}.battenberg.txt)")
    }

    message("Found ", length(samples), " sample(s): ", paste(samples, collapse = ", "))
    for (s in samples) {
        vcf_path <- file.path(SAMPLE_DIR, s, paste0(s, ".mutect.vcf"))
        bat_path <- file.path(SAMPLE_DIR, s, paste0(s, ".battenberg.txt"))
        tryCatch(
            process_one_sample(
                sample_id = s,
                vcf_file = vcf_path,
                battenberg_file = bat_path,
                out_ssm_path = file.path(OUT_MUTATION_TABLE, paste0(s, "_ssm_data.txt")),
                out_cnv_path = file.path(OUT_COPY_NUMBER, paste0(s, "_cnv_data.txt"))
            ),
            error = function(e) message("Failed ", s, ": ", conditionMessage(e))
        )
    }
}

main()
