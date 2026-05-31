# =============================================================================
# 15_NicheNet_analysis.R
#
# Current purpose:
#   Supplementary Figure S22 locked status / reproduction guard.
#
# Important:
#   This script intentionally does NOT rerun NicheNet.
#
# Rationale:
#   The manuscript-matched S22 result is locked with:
#     sender cells   = GlycoHigh hepatocytes, n = 5,589
#     receiver cells = tumor-derived T/NK + myeloid, n = 11,383
#     expressed ligands = 325
#     MIF rank = 37, AUPR = 0.091
#     SPP1 rank = 250, AUPR = 0.026
#
#   Recomputing S22 from the current seurat_final.rds median split would use the
#   current GlycoHigh hepatocyte count of 7,695 and may generate a different
#   NicheNet ranking. That would be less safe than preserving the already
#   manuscript-matched S22 output.
#
# Generates:
#   results/S22_locked_status.csv
#   results/S22_locked_focal_ligands.csv
#   results/S22_locked_qc.csv
#   results/S22_locked_file_presence.csv
#   results/S22_locked_reproduction_note.txt
#   results/S22_locked_manifest.csv
#
# Does NOT generate or overwrite:
#   Supplementary Figure S22 PDF/PNG
# =============================================================================

options(stringsAsFactors = FALSE)

out_dir <- "results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_csv <- function(x, filename) {
  path <- file.path(out_dir, filename)
  utils::write.csv(x, path, row.names = FALSE, quote = TRUE)
  message("Wrote: ", normalizePath(path, winslash = "/", mustWork = FALSE))
  invisible(path)
}

write_text <- function(lines, filename) {
  path <- file.path(out_dir, filename)
  writeLines(lines, con = path, useBytes = TRUE)
  message("Wrote: ", normalizePath(path, winslash = "/", mustWork = FALSE))
  invisible(path)
}

# -----------------------------------------------------------------------------
# Locked S22 manuscript values
# -----------------------------------------------------------------------------

locked <- list(
  script = "15_NicheNet_analysis.R",
  figure = "Supplementary Figure S22",
  analysis = "NicheNet ligand activity analysis",
  status = "LOCKED_STATUS_REPRODUCTION_GUARD",
  recompute_nichenet = "NO",
  sender_definition = "GlycoHigh hepatocytes",
  sender_n = 5589L,
  receiver_definition = "tumor-derived T/NK + myeloid",
  receiver_n = 11383L,
  expressed_ligands_n = 325L,
  ligand_activity_panel = "Top 30 ligand activity",
  mif_rank = 37L,
  mif_aupr = 0.091,
  spp1_rank = 250L,
  spp1_aupr = 0.026,
  current_object_warning = paste(
    "Do not regenerate S22 from the current seurat_final.rds median split;",
    "that current object uses GlycoHigh = 7,695 tumor-derived hepatocytes,",
    "whereas the locked manuscript S22 sender count is 5,589."
  )
)

status_df <- data.frame(
  field = c(
    "script",
    "figure",
    "analysis",
    "status",
    "recompute_nichenet",
    "sender_definition",
    "sender_n",
    "receiver_definition",
    "receiver_n",
    "expressed_ligands_n",
    "ligand_activity_panel",
    "current_object_warning"
  ),
  value = c(
    locked$script,
    locked$figure,
    locked$analysis,
    locked$status,
    locked$recompute_nichenet,
    locked$sender_definition,
    as.character(locked$sender_n),
    locked$receiver_definition,
    as.character(locked$receiver_n),
    as.character(locked$expressed_ligands_n),
    locked$ligand_activity_panel,
    locked$current_object_warning
  )
)

focal_ligands_df <- data.frame(
  ligand = c("MIF", "SPP1"),
  locked_rank = c(locked$mif_rank, locked$spp1_rank),
  locked_AUPR = c(locked$mif_aupr, locked$spp1_aupr),
  manuscript_status = c("LOCKED", "LOCKED"),
  note = c(
    "Manuscript-matched S22 ligand activity value; do not overwrite by current-object rerun.",
    "Manuscript-matched S22 ligand activity value; do not overwrite by current-object rerun."
  )
)

# -----------------------------------------------------------------------------
# QC/status checks
# These checks confirm that this script is behaving as a locked-status guard.
# They do not claim to fully rerun NicheNet.
# -----------------------------------------------------------------------------

qc_df <- data.frame(
  check_id = c(
    "mode_is_locked_status_guard",
    "nichenet_recompute_disabled",
    "sender_n_locked",
    "receiver_n_locked",
    "expressed_ligands_n_locked",
    "mif_rank_locked",
    "mif_aupr_locked",
    "spp1_rank_locked",
    "spp1_aupr_locked",
    "no_figure_overwrite"
  ),
  expected = c(
    "LOCKED_STATUS_REPRODUCTION_GUARD",
    "NO",
    "5589",
    "11383",
    "325",
    "37",
    "0.091",
    "250",
    "0.026",
    "Do not generate or overwrite S22 PDF/PNG"
  ),
  observed = c(
    locked$status,
    locked$recompute_nichenet,
    as.character(locked$sender_n),
    as.character(locked$receiver_n),
    as.character(locked$expressed_ligands_n),
    as.character(locked$mif_rank),
    sprintf("%.3f", locked$mif_aupr),
    as.character(locked$spp1_rank),
    sprintf("%.3f", locked$spp1_aupr),
    "Do not generate or overwrite S22 PDF/PNG"
  ),
  pass = c(
    identical(locked$status, "LOCKED_STATUS_REPRODUCTION_GUARD"),
    identical(locked$recompute_nichenet, "NO"),
    identical(locked$sender_n, 5589L),
    identical(locked$receiver_n, 11383L),
    identical(locked$expressed_ligands_n, 325L),
    identical(locked$mif_rank, 37L),
    isTRUE(all.equal(locked$mif_aupr, 0.091)),
    identical(locked$spp1_rank, 250L),
    isTRUE(all.equal(locked$spp1_aupr, 0.026)),
    TRUE
  )
)

if (!all(qc_df$pass)) {
  write_csv(qc_df, "S22_locked_qc.csv")
  stop("S22 locked-status QC failed. Do not proceed.")
}

# -----------------------------------------------------------------------------
# Optional file-presence report
# This is informational only. Missing files here do not fail the guard script,
# because this script is meant to document locked manuscript values rather than
# rebuild historical NicheNet output.
# -----------------------------------------------------------------------------

candidate_files <- c(
  "FigS22_NicheNet_ligand_activity_top30.pdf",
  "FigS22_NicheNet_ligand_activity_top30.png",
  file.path(out_dir, "FigS22_NicheNet_ligand_activity_top30.pdf"),
  file.path(out_dir, "FigS22_NicheNet_ligand_activity_top30.png"),
  file.path(out_dir, "S22_NicheNet_ligand_activities.csv"),
  file.path(out_dir, "S22_NicheNet_top30_ligand_activity.csv")
)

file_presence_df <- data.frame(
  file = candidate_files,
  present = file.exists(candidate_files),
  note = ifelse(
    file.exists(candidate_files),
    "Existing file detected; this script does not modify it.",
    "Not detected in this local checkout; this is informational only."
  )
)

# -----------------------------------------------------------------------------
# Reproduction note
# -----------------------------------------------------------------------------

note_lines <- c(
  "Supplementary Figure S22 locked status / reproduction guard",
  "==========================================================",
  "",
  "This script intentionally does not rerun NicheNet.",
  "",
  "Locked manuscript S22 values:",
  "  sender cells   = GlycoHigh hepatocytes, n = 5,589",
  "  receiver cells = tumor-derived T/NK + myeloid, n = 11,383",
  "  expressed ligands = 325",
  "  MIF rank = 37; MIF AUPR = 0.091",
  "  SPP1 rank = 250; SPP1 AUPR = 0.026",
  "  ligand activity panel = Top 30 ligand activity",
  "",
  "Reason for locked/status mode:",
  "  The current seurat_final.rds median split uses GlycoHigh = 7,695",
  "  tumor-derived hepatocytes, whereas the manuscript-matched S22",
  "  NicheNet sender count is 5,589. Recomputing from the current object",
  "  may produce a different ligand ranking and would risk overwriting",
  "  the manuscript-matched S22 result.",
  "",
  "For a full computational rerun of S22, the original S22-specific",
  "NicheNet input/output objects that produced sender n = 5,589 should",
  "be restored first. Until then, this script records and protects the",
  "locked manuscript values.",
  "",
  "This script does not create or overwrite any S22 PDF/PNG figure."
)

# -----------------------------------------------------------------------------
# Write outputs
# -----------------------------------------------------------------------------

status_path <- write_csv(status_df, "S22_locked_status.csv")
ligand_path <- write_csv(focal_ligands_df, "S22_locked_focal_ligands.csv")
qc_path <- write_csv(qc_df, "S22_locked_qc.csv")
presence_path <- write_csv(file_presence_df, "S22_locked_file_presence.csv")
note_path <- write_text(note_lines, "S22_locked_reproduction_note.txt")

manifest_df <- data.frame(
  output_file = c(
    basename(status_path),
    basename(ligand_path),
    basename(qc_path),
    basename(presence_path),
    basename(note_path)
  ),
  purpose = c(
    "Locked S22 status summary",
    "Locked focal ligand ranks and AUPR values",
    "QC confirming locked-status guard mode",
    "Informational check for existing local S22 files",
    "Human-readable reproduction note"
  ),
  overwrites_manuscript_figure = c("NO", "NO", "NO", "NO", "NO")
)

write_csv(manifest_df, "S22_locked_manifest.csv")

message("")
message("S22 locked/status guard complete.")
message("No NicheNet recomputation was performed.")
message("No S22 PDF/PNG figure was generated or overwritten.")
message("")
message("Locked manuscript values:")
message("  sender cells   = GlycoHigh hepatocytes, n = ", locked$sender_n)
message("  receiver cells = tumor-derived T/NK + myeloid, n = ", locked$receiver_n)
message("  expressed ligands = ", locked$expressed_ligands_n)
message("  MIF  rank = ", locked$mif_rank,  ", AUPR = ", locked$mif_aupr)
message("  SPP1 rank = ", locked$spp1_rank, ", AUPR = ", locked$spp1_aupr)
