"""Build a standalone, offline HTML dashboard: Earth_Signal_Dashboard.html

Run:  python build_dashboard.py
Then double-click Earth_Signal_Dashboard.html. No server, no internet, no IDE needed.

It runs the same trend engine as the Streamlit app (earthsignal/) for every variable at three
false-discovery-rate levels and embeds the results in one HTML file. If real processed files exist in
data/processed/ (and xarray is installed) they are used; otherwise the synthetic demo data is used and the
dashboard says so on screen.
"""
from __future__ import annotations

import base64
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from earthsignal import analysis, catalog, clusters, synthetic  # noqa: E402
from earthsignal.geo import REGIONS, area_grid, box_mask, wmean  # noqa: E402

ALPHAS = [0.05, 0.10, 0.20]
N_PERM = 50
METHOD, CORRECTION = "bh", "ar1"


def load_data() -> dict:
    try:
        from earthsignal import io  # needs xarray
        return io.load_all()
    except ImportError:
        lat, lon = synthetic.make_grid(3.0)
        years = np.arange(2003, 2026)
        return {v: synthetic.make_variable(v, years, lat, lon) for v in catalog.IDS}


def clean(o):
    """Make numpy / pandas / NaN values JSON safe."""
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, np.ndarray):
        return clean(o.tolist())
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else round(float(o), 6)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return o


def digits(a: np.ndarray) -> str:
    return "".join(str(int(x)) for x in a.ravel())


def encode_values(vals: np.ndarray) -> dict:
    fin = np.isfinite(vals)
    lo, hi = float(np.nanmin(vals)), float(np.nanmax(vals))
    scale = (hi - lo) / 60000 or 1.0
    q = np.where(fin, np.round((vals - lo) / scale) - 30000, -32768).astype("<i2")
    return {"b64": base64.b64encode(q.tobytes()).decode(), "lo": lo, "scale": scale, "shape": list(vals.shape)}


def main() -> None:
    data = load_data()
    ids = catalog.IDS
    ref = data[ids[0]]
    lat, lon = ref["lat"], ref["lon"]
    area = area_grid(lat, lon)
    years = [int(y) for y in ref["years"]]

    res = {a: {v: analysis.run_variable(data[v], alpha=a, method=METHOD, correction=CORRECTION) for v in ids} for a in ALPHAS}

    out = {
        "lat": lat, "lon": lon, "area": area[:, 0], "years": years,
        "regions": {k: list(b) for k, b in REGIONS.items()},
        "alphas": ALPHAS, "method": METHOD,
        "vars": {}, "stats": {}, "rank": {}, "sens": {}, "opp": {}, "clusters": {}, "grade": {}, "contact": {},
        "pstar": {}, "null": {},
    }

    for v in ids:
        d, r0 = data[v], res[0.10][v]
        meta = catalog.CATALOG[v]
        out["vars"][v] = {
            "label": meta["label"], "units": meta["units"], "short": meta["short"], "sphere": meta["sphere"],
            "source": meta["source"], "note": meta["note"], "synthetic": d["source"] == "synthetic",
            "land": digits(d["land"]),
            "slope": np.round(r0["slope_dec"], 5).ravel(), "pct": np.round(r0["pct_dec"], 3).ravel(),
            "p": np.round(r0["p"], 5).ravel(), "z": np.round(r0["z"], 2).ravel(),
            "cov": np.round(r0["coverage"], 2).ravel(), "snr": np.round(r0["snr"], 2).ravel(),
            "stable": digits(r0["stable"]), "infl": np.round(r0["infl"], 2).ravel(),
            "values": encode_values(d["values"]),
        }

    for a in ALPHAS:
        ak = f"{a:.2f}"
        out["stats"][ak], out["rank"][ak], out["sens"][ak] = {}, {}, {}
        out["opp"][ak], out["clusters"][ak], out["grade"][ak], out["contact"][ak], out["pstar"][ak] = {}, {}, {}, {}, {}
        for v in ids:
            r = res[a][v]
            out["grade"][ak][v] = digits(r["grade"])
            out["contact"][ak][v] = digits(r["contact"])
            out["pstar"][ak][v] = r["p_star"]
            o = r["opp"].copy()
            out["opp"][ak][v] = o.round(3).to_dict("records") if len(o) else []
            out["clusters"][ak][v] = clusters.cluster_table(r["cid"], r["csign"], r["slope_dec"], lat, lon).round(3).to_dict("records")
        for rn, box in REGIONS.items():
            out["stats"][ak][rn] = {}
            out["sens"][ak][rn] = {}
            rk = analysis.region_ranking(res[a], box, alpha=a, method=METHOD)
            out["rank"][ak][rn] = rk.to_dict("records") if len(rk) else []
            for v in ids:
                r = res[a][v]
                m = box_mask(lat, lon, box) & r["land"] & np.isfinite(r["p"])
                sig = r["sig"] & m
                tot = area[m].sum()
                opp_cells = int((r["contact"] & box_mask(lat, lon, box) & r["land"]).sum())
                opp_pairs = int(((r["opp"]["contact_lat"] >= box[0]) & (r["opp"]["contact_lat"] <= box[1])).sum()) if len(r["opp"]) else 0
                if m.sum() >= 5:
                    out["stats"][ak][rn][v] = {
                        "mean": wmean(r["slope_dec"][m], area[m]),
                        "pos": float(area[sig & (r["slope_dec"] > 0)].sum() / tot),
                        "neg": float(area[sig & (r["slope_dec"] < 0)].sum() / tot),
                        "contacts": opp_cells,
                        "summary": analysis.summary_text(v, r, rn, box, opp_pairs).replace(" Significance does not prove a physical cause.", ""),
                    }
                else:
                    out["stats"][ak][rn][v] = None
                sens = analysis.aoi_sensitivity(r, box, alpha=a, method=METHOD)
                out["sens"][ak][rn][v] = sens.to_dict("records") if len(sens) else []

    # year-shuffle controls (default alpha only, fixed seed)
    for v in ids:
        n = analysis.opposition_null(data[v], n_perm=N_PERM, seed=0, alpha=0.10, method=METHOD, correction=CORRECTION)
        out["null"][v] = {"observed": n["observed"], "null": n["null"], "p": n["p_value"], "n_perm": n["n_perm"]}

    out["any_synthetic"] = any(data[v]["source"] == "synthetic" for v in ids)

    tpl = (ROOT / "dashboard_template.html").read_text(encoding="utf-8")
    html = tpl.replace("/*__DATA__*/null", json.dumps(clean(out), separators=(",", ":")))
    target = ROOT / "Earth_Signal_Dashboard.html"
    target.write_text(html, encoding="utf-8")
    print(f"Wrote {target} ({target.stat().st_size / 1e6:.2f} MB)  synthetic={out['any_synthetic']}")


if __name__ == "__main__":
    main()
