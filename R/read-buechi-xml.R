
read_buechi_xml <- function(file) {
  
  if (!file.exists(file)) {
    stop(glue::glue("File {file} does not exist"))
  }
  
  xml_in <- xml2::read_xml(x = file)
  xml_list <- xml2::as_list(xml_in)
  data <- xml_list[["CATDataSetImportExport"]]
  
  # Get the metadata XML elements ----------------------------------------------
  
  # keys/IDs
  meta_ids_dfs <- xml_elem_to_dfs(data, tag = "RE_SampleSpectrum")
  # measurement metadata
  meta_spc_dfs <- xml_elem_to_dfs(data, tag = "RE_Spectrum")
  # sample metadata
  meta_sample_dfs <- xml_elem_to_dfs(data, tag = "SM_Sample")
  
  ## Get spectra data
  spc_lst <- data[names(data) %in% "RE_SpectrumData"]
  
  # Decode base64 spectra ------------------------------------------------------
  
  # map number of points / wavelengths
  npoints_lst <- purrr::map(meta_spc_dfs, ~ as.integer(.x[["XLength"]]))
  # extract list of base64 character vectors
  spc_base64_lst <- 
    purrr::modify_depth(spc_lst, .depth = 1, "YVector") %>% purrr::flatten()
  # Decode base64 to raw vector and interpret raw/binary as double
  spc_vec_lst <- purrr::map2(
    .x = spc_base64_lst,
    .y = npoints_lst,
    .f = ~ readBin(base64enc::base64decode(.x), what = "double", n = .y)
  )
  
  # Generate data.tables of spectra --------------------------------------------
  
  # generate wavenumbers (revert vectors to have high to low wavenumbers)
  xstart_lst <- purrr::map(meta_spc_dfs, ~ as.integer(.x[["XStart"]]))
  xdelta_lst <- purrr::map(meta_spc_dfs, ~ as.integer(.x[["DataResolution"]]))
  xvalues_lst <- purrr::pmap(
    .l= list(
      ..1 = xstart_lst,
      ..2 = npoints_lst,
      ..3 = xdelta_lst
    ),
    .f = ~ rev(seq(from = ..1, to = ..1 + (..2 - 1L) * ..3, by = ..3))
  )
  
  # make list of 1 row data.table's
  spc_dt_lst <- purrr::map2(.x = spc_vec_lst, .y = xvalues_lst,
    ~ vec_as_dt(x = t(rev(.x)), colnames = as.character(.y)))
  
  # Join ID, sample, and measurement metadata ----------------------------------
  
  meta_sample_df <- 
    dplyr::bind_rows(meta_sample_dfs) %>% 
    dplyr::rename(GUID_sample = GUID) %>%
    tibble::add_column(sample_id = .[["Name"]], .before = 1L)
  
  meta_spc_dfs <- purrr::map(meta_spc_dfs, 
    ~ dplyr::rename(.x, GUID_spc = GUID, TimeStamp_spc = TimeStamp))
  
  # combine ID metadata (SampleSpectrum, Sample, Spectrum) with sample metadata
  # `list(sample_meta_df)` is recycled
  meta_dfs <- join_metadata(
    x = meta_ids_dfs, y = list(meta_sample_df), by = "SampleID")
  
  # combine metadata with measurement metadata
  meta_all_dfs <- join_metadata(
    x = meta_dfs, y = meta_spc_dfs, by = "SpectrumID")
  
  # Rearrange spectrum data for list output ------------------------------------
  
  # Generate element names for output list; each element is for one spectrum
  sample_spc_id_vec <- purrr::map_chr(meta_all_dfs, ~ .x[["SampleSpectrumID"]])
  sample_id_vec <- purrr::map_chr(meta_all_dfs, ~ .x[["sample_id"]])
  unique_id_vec <- paste0(sample_spc_id_vec, "_", sample_id_vec)
  
  # Append missing `file_id` to metadata ; for later `simplerspec::gather_spc()`
  meta_all_dfs <- purrr::map2(.x = meta_all_dfs, .y = unique_id_vec,
    ~ .x %>%
      tibble::add_column(unique_id = .y, .before = 1L) %>%
      tibble::add_column(file_id = .y, .after = 1L)    
  )
  
  data_out_lst <- purrr:::pmap(
    .l = list(
      ..1 = meta_all_dfs,
      ..2 = spc_dt_lst,
      ..3 = xvalues_lst
    ),
    .f = ~ list(
      "metadata" = ..1,
      "spc" = ..2,
      "wavenumbers" = ..3
    )
  ) %>% rlang::set_names(nm = unique_id_vec)
  
  return(data_out_lst)
}


# Helpers ----------------------------------------------------------------------

# Collect xml_elements to data frames
xml_elem_to_dfs <- function(data, tag, .depth = 1) {
  stopifnot(is.character(tag) && length(tag == 1))
  
  data_sel <- data[names(data) %in% tag]
  
  purrr::modify_depth(
    .x = data_sel, 
    .depth, 
    .f = ~ tibble::as_tibble(purrr::flatten(.x))
  )
}

# Mapping joins for metadata
join_metadata <- function(x, y, by) { 
  purrr::map2(.x = x, .y = y, ~ dplyr::inner_join(x = .x, y = .y, by))
}

vec_as_dt <- function(x, colnames) {
  dt <- data.table::data.table(x)
  data.table::setnames(dt, colnames)
}


## Notes for XML structure =====================================================

## Spectrum metadata (6 blocks)
# <RE_SampleSpectrum> :
#   <SampleSpectrumID>
#   <SampleID>        
#   <SpectrumID>  

## Measurement metadata (6 blocks)
# <RE_Spectrum> :
#   <SpectrumID>
#   <SpectrumType>
#   <SpectrumSourceType>
#   <GUID>
#   <SpectrumCharacteristicType>
#   <NumberOfScans>
#   <Comment>
#   <VersionNo>
#   <InstrumentResolution>
#   <DataResolution>
#   <InstrumentType>
#   <InstrumentSerialNo>
#   <MeasurementCellType>
#   <MeasurementCellSerialNo>
#   <MeasurementCellOptionType>
#   <TimeStamp>
#   <UserName>A
#   <SpectrumXAxisType>
#   <SpectrumYAxisType>
#   <XStart>
#   <InstrumentTemperature>
#   <SampleTemperature>
#   <XLength>
#   <GainFactor>
#   <Gain>

## Spectra y data (6 blocks)
# <RE_SpectrumData>
#   <SpectrumID>
#   <YVector>

## Sample metadata (2 blocks, one per sample; 3 spectra per sample)
# <SM_Sample>
#   <SampleID>
#   <Name>
#   <GUID>
#   <Comment>
#   <Origin>
#   <Supplier>
#   <ApplicationName>
#   <TimeStamp>