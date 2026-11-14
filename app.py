"""Earth Signal: Streamlit dashboard.  Run:  streamlit run app.py"""
from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from earthsignal import analysis, catalog, clusters, geo, io
from earthsignal.geo import area_grid, box_mask

st.set_page_config(page_title="Earth Signal", layout="wide")

DECREASE, INCREASE, NS_GRAY = "#5e3c99", "#e66101", "#d9d9d9"   # colour-blind-safe diverging pair
SCALE = [[0, DECREASE], [0.25, "#b2abd2"], [0.5, "#f7f7f7"], [0.75, "#fdb863"], [1, INCREASE]]
GRADE_NAMES = {0: "Not significant", 1: "Low", 2: "Medium", 3: "High"}


# ---------------------------------------------------------------- data
@st.cache_data(show_spinner="Loading data...")
def get_data():
    return io.load_all()


@st.cache_data(show_spinner="Running trend tests...")
def get_result(var_id: str, y0: int, y1: int, alpha: float, method: str, correction: str, min_area: float, gap_deg: float):
    return analysis.run_variable(get_data()[var_id], y0, y1, alpha=alpha, method=method,
                                 correction=correction, min_area_km2=min_area, gap_deg=gap_deg)


@st.cache_data(show_spinner="Running year-shuffle controls...")
def get_null(var_id: str, y0: int, y1: int, alpha: float, method: str, correction: str, min_area: float, gap_deg: float, n_perm: int):
    return analysis.opposition_null(get_data()[var_id], n_perm=n_perm, seed=0, y0=y0, y1=y1, alpha=alpha,
                                    method=method, correction=correction, min_area_km2=min_area, gap_deg=gap_deg)


data = get_data()
ids = catalog.IDS

# ---------------------------------------------------------------- sidebar
sb = st.sidebar
sb.title("Earth Signal")
sb.caption("Which Earth system variables are changing, where, by how much, and is it significant?")
var_id = sb.selectbox("Variable", ids, format_func=catalog.label)
yrs_all = data[var_id]["years"].astype(int)
y0, y1 = sb.slider("Period (years)", int(yrs_all.min()), int(yrs_all.max()), (int(yrs_all.min()), int(yrs_all.max())))
region_name = sb.selectbox("Region", list(geo.REGIONS) + ["Custom box"])
if region_name == "Custom box":
    la0, la1 = sb.slider("Latitude", -90, 90, (-30, 30))
    lo0, lo1 = sb.slider("Longitude", -180, 180, (-60, 60))
    box = (la0, la1, lo0, lo1)
else:
    box = geo.REGIONS[region_name]

with sb.expander("Statistics settings", expanded=False):
    alpha = st.select_slider("False discovery rate (alpha)", [0.01, 0.05, 0.10, 0.20], value=0.10)
    method = st.radio("Multiple-testing method", ["bh", "by"], format_func=lambda m: {"bh": "Benjamini-Hochberg", "by": "Benjamini-Yekutieli (stricter)"}[m])
    correction = st.radio("Serial-correlation correction", ["ar1", "hr", "none"],
                          format_func=lambda c: {"ar1": "AR(1) effective sample size (default)", "hr": "Hamed-Rao", "none": "None (plain Mann-Kendall)"}[c])
    min_area = 1000.0 * st.slider("Minimum cluster area (thousand km2)", 50, 2000, 400, step=50)
    gap_deg = st.slider("Max gap for 'adjacent' clusters (degrees)", 0.0, 15.0, 8.0, step=1.0)
show_opp = sb.toggle("Outline opposite-trend contacts", value=True)

if y1 - y0 + 1 < 12:
    st.error("Choose at least 12 years: shorter records cannot support a meaningful trend test.")
    st.stop()

r = get_result(var_id, y0, y1, alpha, method, correction, min_area, gap_deg)
meta = catalog.CATALOG[var_id]
lat, lon, area = r["lat"], r["lon"], area_grid(r["lat"], r["lon"])
inbox = box_mask(lat, lon, box) & r["land"]

# ---------------------------------------------------------------- header
st.title("Earth Signal")
if data[var_id]["source"] == "synthetic":
    st.warning("**Demo mode: this variable is synthetic.** Trends shown are planted by the demo generator, not real observations. "
               "Place processed NASA files in `data/processed/` (see docs/DATA.md) to analyse real data.")
else:
    st.success("Showing processed observation data from data/processed/.")
st.markdown(f"**{meta['label']}** ({meta['units']} per decade) - {region_name}, {y0}-{y1}")

c = st.columns(4)
sig_in = r["sig"] & inbox
tot_area = area[inbox].sum() or np.nan
c[0].metric("Area-weighted mean change / decade", f"{geo.wmean(r['slope_dec'][inbox], area[inbox]):+.3g}")
c[1].metric("Area significantly increasing", f"{area[sig_in & (r['slope_dec'] > 0)].sum() / tot_area:.0%}")
c[2].metric("Area significantly decreasing", f"{area[sig_in & (r['slope_dec'] < 0)].sum() / tot_area:.0%}")
c[3].metric("Opposite-trend contact cells", int((r["contact"] & inbox).sum()))

tab_map, tab_cell, tab_rank, tab_opp, tab_method = st.tabs(["Trend map", "Cell inspector", "Variable ranking", "Opposite trends", "Method and limits"])

# ---------------------------------------------------------------- map
with tab_map:
    lim = float(np.nanpercentile(np.abs(r["slope_dec"][r["sig"]]), 98)) if r["sig"].any() else 1.0
    z_sig = np.where(r["sig"], r["slope_dec"], np.nan)
    z_ns = np.where(r["land"] & ~r["sig"], 1.0, np.nan)
    fig = go.Figure()
    fig.add_trace(go.Heatmap(x=lon, y=lat, z=z_ns, colorscale=[[0, NS_GRAY], [1, NS_GRAY]], showscale=False,
                             hoverinfo="skip", name="not significant"))
    fig.add_trace(go.Heatmap(x=lon, y=lat, z=z_sig, colorscale=SCALE, zmid=0, zmin=-lim, zmax=lim,
                             colorbar=dict(title=f"{meta['units']}<br>per decade"), name="significant",
                             hovertemplate="lat %{y:.1f}, lon %{x:.1f}<br>%{z:.3g} per decade<extra></extra>"))
    if show_opp and r["contact"].any():
        LON, LAT = np.meshgrid(lon, lat)
        fig.add_trace(go.Scatter(x=LON[r["contact"]], y=LAT[r["contact"]], mode="markers", name="opposite-trend contact",
                                 marker=dict(symbol="circle-open", size=11, color="black", line=dict(width=2))))
    for b in ([box] if region_name != "Global" else []):
        fig.add_shape(type="rect", x0=b[2], x1=b[3], y0=b[0], y1=b[1], line=dict(color="black", dash="dot"))
    fig.update_layout(height=500, margin=dict(l=10, r=10, t=10, b=10), plot_bgcolor="white",
                      xaxis=dict(range=[-180, 180], title="Longitude", tickvals=list(range(-180, 181, 60))),
                      yaxis=dict(range=[-90, 90], title="Latitude", tickvals=list(range(-90, 91, 30))),
                      legend=dict(orientation="h", y=-0.15))
    st.plotly_chart(fig, width="stretch")
    st.caption(f"Coloured cells are significant after {method.upper()} false-discovery-rate control at alpha={alpha:.2f} "
               f"(p* = {r['p_star']:.4f}). Grey cells have data but no significant trend. White cells are ocean or masked.")
    opp_in = int(((r["opp"]["contact_lat"] >= box[0]) & (r["opp"]["contact_lat"] <= box[1])).sum()) if len(r["opp"]) else 0
    st.info(analysis.summary_text(var_id, r, region_name, box, opp_in))

    tbl = pd.DataFrame({"lat": np.repeat(lat, lon.size), "lon": np.tile(lon, lat.size),
                        "slope_per_decade": r["slope_dec"].ravel(), "pct_per_decade": r["pct_dec"].ravel(),
                        "p_value": r["p"].ravel(), "significant": r["sig"].ravel(),
                        "reliability": [GRADE_NAMES[g] for g in r["grade"].ravel()],
                        "coverage": r["coverage"].ravel()})
    tbl = tbl[inbox.ravel() & np.isfinite(tbl["slope_per_decade"])]
    st.download_button("Download cell results (CSV)", tbl.to_csv(index=False), file_name=f"earth_signal_{var_id}.csv", mime="text/csv")

# ---------------------------------------------------------------- cell inspector
with tab_cell:
    cc = st.columns(2)
    zabs = np.where(r["sig"] & inbox, np.abs(r["z"]), -1.0)
    bi, bj = np.unravel_index(int(np.argmax(zabs)), zabs.shape) if zabs.max() > 0 else (int(np.argmin(abs(lat - (box[0] + box[1]) / 2))), 0)
    qlat = cc[0].number_input("Latitude", -90.0, 90.0, float(lat[bi]), step=1.0, key=f"qlat_{var_id}_{region_name}")
    qlon = cc[1].number_input("Longitude", -180.0, 180.0, float(lon[bj]), step=1.0, key=f"qlon_{var_id}_{region_name}")
    st.caption("Starts at the most significant cell in the selected region. Type any coordinates to inspect another cell.")
    i, j = int(np.argmin(abs(lat - qlat))), int(np.argmin(abs(lon - qlon)))
    d = data[var_id]
    sel = (d["years"] >= y0) & (d["years"] <= y1)
    yy, vv = d["years"][sel], d["values"][sel][:, i, j]
    if not np.isfinite(vv).any():
        st.warning(f"No data at {lat[i]:.1f}, {lon[j]:.1f} (ocean or masked). Pick another location.")
    else:
        slope_yr = r["slope_dec"][i, j] / 10
        ok = np.isfinite(vv)
        icpt = np.nanmedian(vv[ok] - slope_yr * (yy[ok] - yy.mean()))
        f2 = go.Figure()
        f2.add_trace(go.Scatter(x=yy[ok], y=vv[ok], mode="markers+lines", name="annual value", line=dict(color="#888", width=1)))
        f2.add_trace(go.Scatter(x=yy, y=icpt + slope_yr * (yy - yy.mean()), mode="lines", name="Theil-Sen trend",
                                line=dict(color=INCREASE if slope_yr > 0 else DECREASE, width=3)))
        f2.update_layout(height=340, margin=dict(l=10, r=10, t=30, b=10), title=f"Cell {lat[i]:.1f}, {lon[j]:.1f}",
                         yaxis_title=meta["units"], xaxis_title="Year")
        st.plotly_chart(f2, width="stretch")
        g = int(r["grade"][i, j])
        m1 = st.columns(5)
        m1[0].metric("Slope / decade", f"{r['slope_dec'][i, j]:+.3g}")
        m1[1].metric("Change / decade", f"{r['pct_dec'][i, j]:+.1f}%")
        m1[2].metric("p-value", f"{r['p'][i, j]:.4f}")
        m1[3].metric("Significant (FDR)", "yes" if r["sig"][i, j] else "no")
        m1[4].metric("Reliability", GRADE_NAMES[g])
        st.table(pd.DataFrame({
            "check": ["Significant after false-discovery control", "Sign holds after dropping first/last years",
                      "Data coverage at least 70%", "Total change at least 1x year-to-year noise"],
            "result": ["pass" if r["sig"][i, j] else "fail", "pass" if r["stable"][i, j] else "fail",
                       "pass" if r["coverage"][i, j] >= 0.7 else "fail", "pass" if r["snr"][i, j] >= 1 else "fail"],
            "value": [f"p={r['p'][i, j]:.4f} (threshold {r['p_star']:.4f})", "", f"{r['coverage'][i, j]:.0%}", f"{r['snr'][i, j]:.2f}"]}))
        st.caption(f"Serial-correlation variance inflation for this cell: {r['infl'][i, j]:.2f}x. "
                   "Reliability is a screening grade built from these checks, not a probability.")

# ---------------------------------------------------------------- ranking
with tab_rank:
    st.write(f"Which variables show the clearest significant change in **{region_name}**? (FDR re-applied inside the region.)")
    synth = [catalog.label(v) for v in ids if data[v]["source"] == "synthetic"]
    if synth:
        st.warning("Synthetic demo data in this ranking: " + ", ".join(synth) + ". Ranking real and synthetic variables together is not meaningful.")
    res_all = {v: get_result(v, y0, y1, alpha, method, correction, min_area, gap_deg) for v in ids}
    rk = analysis.region_ranking(res_all, box, alpha=alpha, method=method)
    if len(rk):
        show = rk.drop(columns=["id"]).copy()
        for col in ("area significant", "increasing", "decreasing"):
            show[col] = show[col].map("{:.0%}".format)
        show["mean slope/decade"] = show["mean slope/decade"].map("{:+.3g}".format)
        show["median SNR (sig. cells)"] = show["median SNR (sig. cells)"].map(lambda x: "-" if not np.isfinite(x) else f"{x:.1f}")
        st.dataframe(show, hide_index=True, width="stretch")
        st.caption("'mixed signs' = at least 3% of the region significantly increasing and 3% significantly decreasing, so the regional mean hides opposing change. "
                   "Ranking is by significant area share, not by importance.")
    else:
        st.warning("Region too small for ranking.")
    st.subheader("Does significance depend on the area you choose?")
    sens = analysis.aoi_sensitivity(res_all[var_id], box, alpha=alpha, method=method)
    if len(sens):
        s2 = sens.copy()
        s2["area significant"] = s2["area significant"].map("{:.0%}".format)
        s2["p* threshold"] = s2["p* threshold"].map("{:.4f}".format)
        st.dataframe(s2, hide_index=True, width="stretch")
        st.caption("FDR is re-applied inside nested boxes around the region centre. If the share changes a lot, "
                   "the result depends on the domain chosen, so report the domain with the result.")

# ---------------------------------------------------------------- opposite trends
with tab_opp:
    st.write("Adjacent clusters of significant cells with opposite sign for the same variable.")
    if len(r["opp"]):
        o = r["opp"].copy()
        for col in ("pos_slope_dec", "neg_slope_dec"):
            o[col] = o[col].round(3)
        o["pos_area_km2"] = o["pos_area_km2"].round(-3)
        o["neg_area_km2"] = o["neg_area_km2"].round(-3)
        st.dataframe(o, hide_index=True, width="stretch")
    else:
        st.info("No adjacent opposite-sign cluster pairs at these settings.")
    ct = clusters.cluster_table(r["cid"], r["csign"], r["slope_dec"], lat, lon)
    with st.expander("All significant clusters"):
        st.dataframe(ct, hide_index=True, width="stretch")
    st.subheader("Could this happen by chance?")
    n_perm = st.slider("Year-shuffle controls", 10, 50, 20, step=10)
    if st.button("Run controls"):
        out = get_null(var_id, y0, y1, alpha, method, correction, min_area, gap_deg, n_perm)
        f3 = go.Figure(go.Histogram(x=out["null"], name="shuffled years", marker_color="#999"))
        f3.add_vline(x=out["observed"], line_color=INCREASE, line_width=3, annotation_text="observed")
        f3.update_layout(height=300, margin=dict(l=10, r=10, t=30, b=10), xaxis_title="opposition contact cells", yaxis_title="runs")
        st.plotly_chart(f3, width="stretch")
        st.write(f"Observed {out['observed']} contact cells; shuffled controls reached at least that in "
                 f"{int(np.sum(out['null'] >= out['observed']))} of {out['n_perm']} runs (p = {out['p_value']:.3f}).")
        st.caption("Years are shuffled identically in every cell, which removes real trends but keeps spatial structure of the noise. "
                   "It also removes year-to-year persistence, so the control is slightly optimistic for persistent variables. "
                   f"With {out['n_perm']} runs the smallest possible p-value is {1 / (out['n_perm'] + 1):.3f}.")

# ---------------------------------------------------------------- method
with tab_method:
    st.markdown(f"""
**Source for this variable:** {meta['source']}
**Caveat:** {meta['note']}

**Pipeline.** Annual values per grid cell, then Theil-Sen slope, Mann-Kendall test with a serial-correlation correction,
then false-discovery-rate control across all cells, then labelling of same-sign clusters and detection of adjacent opposite-sign clusters.

**Reliability grade** (only for significant cells): High = passes all three extra checks, Medium = two, Low = one or none.
The checks are sign stability when the first or last years are dropped, data coverage of at least 70%, and total change at least equal to the year-to-year noise.

**Limits.**
- Records are about 20 to 25 annual values, so power is modest: a total change of about 3x the year-to-year noise is detected about half the time.
- No analytic correction for serial correlation is exact at this length; on simulated noise the default leaves false positives close to nominal,
  but not exactly (see docs/VALIDATION.md).
- Cells are spatially correlated, and the false-discovery procedure assumes weak dependence, so the realised false-discovery share can exceed alpha.
- Significance can change with the region chosen (see the sensitivity table in the ranking tab).
- Trends describe association in the data, not physical cause. Irrigation, land-cover change, sensor changes and cloud gaps can all produce trends.
""")
