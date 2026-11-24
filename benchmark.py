"""Gate 002 baseline scaffold. Recorded data only; no active intervention claims."""
import argparse
import json
from pathlib import Path
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

LAGS = 50
RIDGE = 1e-5
KINDS = ("linear_fir", "polynomial_fir", "wiener_fir")

def unpack(record):
    u, y = np.asarray(record.u, dtype=float).ravel(), np.asarray(record.y, dtype=float).ravel()
    if len(u) != len(y) or len(u) < LAGS or not np.isfinite(u).all() or not np.isfinite(y).all():
        raise ValueError("Expected finite aligned SISO time series of at least 50 samples")
    return u, y

def records(x):
    return list(x) if isinstance(x, (list, tuple)) else [x]

def design(u, kind, linear=None):
    x = sliding_window_view(u, LAGS)
    if kind == "polynomial_fir":
        x = np.column_stack((x, x*x, x*x*x))
    elif kind == "wiener_fir":
        z = x @ linear[1:] + linear[0]
        x = np.column_stack((z, z*z, z*z*z))
    return np.column_stack((np.ones(len(x)), x))

def fit(data, kind, mu, scale, linear=None):
    gram = rhs = None
    n = 0
    for u, y in data:
        # Chunked accumulation; never construct features crossing independent records.
        for start in range(0, len(u)-LAGS+1, 8192):
            stop = min(start+8192, len(u)-LAGS+1)
            X = design((u[start:stop+LAGS-1]-mu)/scale, kind, linear)
            target = y[start+LAGS-1:stop+LAGS-1]
            g, r = X.T @ X, X.T @ target
            gram = g if gram is None else gram+g
            rhs = r if rhs is None else rhs+r
            n += len(X)
    penalty = np.eye(len(rhs))*RIDGE*n
    penalty[0,0] = 0
    return np.linalg.solve(gram+penalty, rhs)

def evaluate(train, test):
    train = [unpack(r) for r in records(train)]
    all_u = np.concatenate([u for u,_ in train])
    mu, scale = float(all_u.mean()), float(all_u.std())
    if scale <= 1e-12:
        raise ValueError("Training input has no excitation")
    coeffs = {}
    for kind in KINDS:
        coeffs[kind] = fit(train, kind, mu, scale, coeffs.get("linear_fir"))
    report = []
    for record in records(test):
        u,y = unpack(record)
        start = max(LAGS-1, int(getattr(record,"state_initialization_window_length",0) or 0))
        if start >= len(y):
            raise ValueError("No test samples after warm-up")
        item = {"record":str(getattr(record,"name","test")), "evaluated_samples":len(y)-start,
                "warmup":start, "models":{}}
        for kind,w in coeffs.items():
            preds = []
            for i in range(0,len(u)-LAGS+1,8192):
                j = min(i+8192,len(u)-LAGS+1)
                X = design((u[i:j+LAGS-1]-mu)/scale,kind,coeffs["linear_fir"])
                preds.append(X @ w)
            pred = np.concatenate(preds)[start-(LAGS-1):]
            err = y[start:]-pred
            item["models"][kind] = {
                "rmse_native_units":float(np.sqrt(np.mean(err**2))),
                "mean_residual":float(err.mean()),
                "parameters":len(w)}
        report.append(item)
    return {"lags":LAGS,"ridge":RIDGE,"test_outputs_used_for_fitting":False,
            "physical_location_inference":"unresolved; predictive fit is not a structural certificate",
            "records":report}

def smoke():
    from types import SimpleNamespace
    rng = np.random.default_rng(19)
    def rec(n):
        u = rng.normal(size=n)
        y = np.convolve(u,[.6,.2,-.1],mode="full")[:n]
        return SimpleNamespace(u=u,y=y,name="synthetic",state_initialization_window_length=50)
    report = evaluate(rec(4000),rec(1500))
    assert report["records"][0]["models"]["linear_fir"]["rmse_native_units"] < .001
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("dataset",choices=["wh","silverbox","smoke"])
    p.add_argument("--output",type=Path)
    args=p.parse_args()
    if args.dataset=="smoke":
        report=smoke()
    else:
        import nonlinear_benchmarks as nb
        loader = nb.WienerHammerBenchMark if args.dataset=="wh" else nb.Silverbox
        train,test=loader()
        report=evaluate(train,test)
    report["dataset"]=args.dataset
    raw=json.dumps(report,indent=2,allow_nan=False)
    print(raw)
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(raw+"\n",encoding="utf-8")

if __name__=="__main__":
    main()
