# Model

This repository now includes a generalized calibration pipeline that supports multiple countries. The `multi_country_model.R` script reads parameter and target CSV files located under `data/<country>` and performs a simple calibration using synthetic inputs. The script produces summary tables and plots for the requested country.

## Running the calibration

```bash
Rscript multi_country_model.R --country rwanda --print-summary
```

Available example countries: `rwanda`, `kenya`, and `uganda`.

Outputs (plots and CSV summaries) are written to the `outputs/<country>` directory.
