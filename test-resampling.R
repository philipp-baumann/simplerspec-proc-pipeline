################################################################################
## Test if new resampling wrapper function `resample_spc()` has identical output
## as old one
################################################################################

library(here)
library(tidyverse)
library(vctrs)
library(simplerspec)

## Read spectra from different spectrometer models =============================

files_alpha <- dir(here("data", "test-data", "bruker-alpha"), full.names = TRUE)
spc_alpha_list <- read_opus_univ(fnames = files_alpha, extract = c("spc"))
spc_alpha <- gather_spc(data = spc_alpha_list)

## Test old vs. new resampling function ========================================

# "old" function:
# https://github.com/philipp-baumann/simplerspec/commit/9f70dade266397d96f5e9639c9706278c4796622

# Default order of `interpol` argument in `prospectr::resample()` changed: 
#   prospectr 0.1.0 on CRAN had "linear" interpolation as default,
#    now in 0.2.0 default is "spline"; causes breaking change if this argument 
#    in code of previous projects: 
# https://github.com/l-ramirez-lopez/prospectr/commit/48a33463f6c14264ee41302e46883c35572b6477#diff-e3dd90c3dc434794e91701bfc0cf53c2
# New `resampling_spc()` wrapper function
source(here("R", "resample-spc-column-in.R"))

# Bruker Alpha (mid-IR) --------------------------------------------------------

spc_alpha_rs_old <-
  spc_alpha %>%
  simplerspec::resample_spc(wn_lower = 500, wn_upper = 3996, wn_interval = 2)

spc_alpha_rs_new <-
  spc_alpha %>%
  resample_spc(wn_lower = 500, wn_upper = 3996, wn_interval = 2)

vctrs::vec_equal(spc_alpha_rs_old, spc_alpha_rs_new)
# Test if corresponding (list-)columns are equall
map2(.x = spc_alpha_rs_old, .y = spc_alpha_rs_new, ~ vec_equal(.x, .y))

# Test equality of first resampled spectrum
vec_equal(
  spc_alpha_rs_old$spc_rs[[1]],
  spc_alpha_rs_new$spc_rs[[1]]
)

# ASD spectra (vis--NIR) -------------------------------------------------------

files_asd <- dir(here("data", "test-data", "asd"), full.names = TRUE)
spc_asd <-
  read_asd_bin(fnames = files_asd) %>%
  simplerspec:::correct_join_offset()

spc_asd_rs_old <-
  spc_asd %>%
  simplerspec::resample_spc(
    x_unit = "wavelength",  wl_lower = 350, wl_upper = 2500, wl_interval = 1)

spc_asd_rs_new <-
  spc_asd %>%
  resample_spc(
    x_unit = "wavelength", wl_lower = 350, wl_upper = 2500, wl_interval = 1)

map2(.x = spc_alpha_rs_old, .y = spc_alpha_rs_new, ~ vec_equal(.x, .y))

vec_equal(
  spc_asd_rs_old$spc_rs[[1]],
  spc_asd_rs_new$spc_rs[[1]]
)

# Büchi spectra (NIR) ----------------------------------------------------------

source(here("R", "read-buechi-xml.R"))

file_buechi <- dir(here("data", "test-data", "buechi-xml"), full.names = TRUE)
spc_buechi_list <- read_buechi_xml(file = file_buechi)
spc_buechi <- gather_spc(spc_buechi_list)

spc_buechi_rs_old <-
  spc_buechi %>%
  simplerspec::resample_spc(wn_lower = 4000, wn_upper = 10000, wn_interval = 4)

spc_buechi_rs_new <-
  spc_buechi %>%
  resample_spc(wn_lower = 4000, wn_upper = 10000, wn_interval = 4)

map2(.x = spc_buechi_rs_old, .y = spc_buechi_rs_new, ~ vec_equal(.x, .y))

vec_equal(
  spc_buechi_rs_old$spc_rs[[1]],
  spc_buechi_rs_new$spc_rs[[1]]
)
