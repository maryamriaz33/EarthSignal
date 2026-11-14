# Earth Signal

A trend-screening dashboard for gridded Earth science data, built for the NASA Space Apps 2026 challenge
"Be an Earth System Trend Detective". Pick a variable, a region and a period; Earth Signal shows what is changing,
where, how much per decade, and whether the change is statistically significant, including places where the same
variable rises in one region and falls in a neighbouring one.

## Status: read this first

| Item | State |
|---|---|
| Trend engine, significance control, opposite-trend detector, dashboard | Built and tested (36 automated tests) |
| Data in the box | **Synthetic demo data with planted trends. Not real observations.** The app shows a warning on every synthetic variable. |
| Real NASA data | Pipeline scripts included; **not run against real NASA files** in the build environment (no Earthdata access). Step-by-step in `docs/DATA.md`. |
| Scientific claims | Only claims about the method are supported (see `docs/VALIDATION.md`). Any finding about the real Earth needs real data first. |

## Quick start

```
pip install -r requirements.txt
streamlit run app.py          # opens the dashboard in demo mode
python -m pytest -q           # runs the tests
python pipeline/validate_engine.py   # regenerates docs/VALIDATION.md
```

## What the dashboard does

| Challenge asks | Earth Signal feature |
|---|---|
| Select variables | Variable catalog (6 variables with source, units, caveats) and a ranking tab that lists which variables show the clearest significant change in the chosen region |
| Visualise change over time | Trend map plus a per-cell time series with the fitted trend |
| What is changing, where, how much | Map of significant cells, slope per decade and percent per decade, region summary sentence, CSV export |
| Is it significant | Theil-Sen slope, Mann-Kendall test with serial-correlation correction, false-discovery-rate control across all cells |
| Opposite trends | Adjacent clusters of significant positive and negative cells, outlined on the map, with a year-shuffle control |

Extras: reliability grade per significant cell (Low, Medium, High), sensitivity of significance to the chosen area,
choice of correction and multiple-testing method.

## Project layout

```
app.py                     Streamlit dashboard
earthsignal/stats.py       Theil-Sen, Mann-Kendall, serial-correlation corrections, FDR, stability, signal-to-noise
earthsignal/clusters.py    Same-sign clusters (longitude wrap) and opposite-trend contacts
earthsignal/analysis.py    One-call analysis per variable, region ranking, area sensitivity, year-shuffle control, summary text
earthsignal/geo.py         Cell areas, region boxes, area-weighted means
earthsignal/catalog.py     Variable catalog (candidate NASA sources, to be verified)
earthsignal/io.py          Loads data/processed/<var>.nc or falls back to synthetic demo data
earthsignal/synthetic.py   Demo data generator with planted trends and known truth
pipeline/                  download_earthaccess.py, build_annual.py, validate_engine.py
tests/                     Statistics, pipeline, and data-conversion tests
docs/                      METHODS.md, DATA.md, VALIDATION.md, SUBMISSION.md
```

## Known limits

- About 20 to 25 annual values per cell give modest power. A total change of about three times the year-to-year noise is found about half the time (`docs/VALIDATION.md`).
- No analytic serial-correlation correction is exact at this length. False positives stay somewhat above nominal for strongly persistent variables.
- Neighbouring cells share noise, so the false-discovery procedure can be optimistic.
- Trends are associations in data. Irrigation, land-cover change, sensor changes and cloud gaps can produce them.
- `build_annual.py` does not apply per-product quality flags; use quality-screened inputs or add masking.
- The 6 variables and their NASA product names in `catalog.py` are candidates; check each on NASA Earthdata Search.
