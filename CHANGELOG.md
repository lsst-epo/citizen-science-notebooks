# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Release Candidate 1]

## 0.5-rc.1
### Changed
- `Citizen_Science_Testing.ipynb` notebook has been updated to install the backend code via `pip install`
- Renamed the `Citizen_Science_TAP_Tutorial.ipynb` notebook to `TAP_Tutorial.ipynb` and moved it to the `experimental_notebooks` subfolder
- Added the `Citizen_Science_Tabular.ipynb` to the `experimental_notebooks` subfolder


## [Unreleased]

## 0.4.0 - 2023-7-26
### Changed
- Added installation of new rubin citsci pypi package

### Removed
- `rubin_citsci_core_pipeline.py` as it is bundled in the rubin citsci pypi package now
- `Citizen_Science_Install.ipynb` as it was installing dependencies now handled by the rubin citsci pypi package

## 0.3.0 - 2023-7-11
### Changed
- Pedagogy text updates by bnevin
- Testing notebook modified to make use of new Install notebook and `rubin_citsci_core_pipeline.py`

### Added
- Migrated majority of the code over from SDK notebook to the new `rubin_citsci_core_pipeline.py`

### Removed
- Removed the SDK notebook

## 0.2.0 - 2023-6-22
### Changed
- Removed the `__cit_sci_data_type` variable from the `Citizen_Science_Testing` notebook and its reference in the `Citizen_Science_SDK` notebook
- Renamed the `send_data()` funtion to `send_image_data()` as part of an effort to define single-use-case functions

## 0.1.0 - 2023-3-02
### Added
- This changelog

### Changed
- Added timestamp and email to error output
- Changed EDC endpoints to reference non-dev version of the RSP-data-exporter service

