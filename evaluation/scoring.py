import json
import os
import sys
import numpy as np
from sklearn.metrics import accuracy_score
import glob
import warnings
from scipy.stats import spearmanr
import pandas as pd
import argparse

required_id_cols = ['barcode', 'pxl_row_in_fullres', 'pxl_col_in_fullres']
protein_whitelist = ['synd','FOXP3','CD16','CD31','CXCL13','Ki67','OLIG2','CXCR5','HLA-A','PD-L1','PSD95','CD20','CD68','CD44','SMA','MSH6','CD23','GFAP','SYNA','Podoplanin','Vimentin','CD47','CD74','SIRP','Granzyme B','IDH1','MPO','CD45','CD21','FIBR','C-KIT','CD3e','TOX','PD-1','PDGFR','CD4','MAP2','CD8','MGMT','CD38','HLA-DR','CD14','ICOS','Granzyme K']
protein_whitelist = None

def _find_single_csv(dir):
    paths = glob.glob(os.path.join(dir, "*.csv"))
    if len(paths) != 1:
        raise FileNotFoundError(
            f"Expected exactly one CSV in {dir}"
        )
    return paths[0]

def _align_on_barcode(pred, gt):
    pred = pred.copy()
    gt = gt.copy()
    pred["barcode"] = pred["barcode"].astype(str)
    gt["barcode"] = gt["barcode"].astype(str)

    if pred['barcode'].duplicated().any():
        warnings.warn("Duplicate barcodes found in prediction; keeping first occurrence")
        pred = pred[~pred['barcode'].duplicated(keep='first')]
    if gt['barcode'].duplicated().any():
        warnings.warn("Duplicate barcodes found in ground truth; keeping first occurrence")
        gt = gt[~gt['barcode'].duplicated(keep='first')]

    common = pd.Index(gt["barcode"]).intersection(pred["barcode"])
    if common.empty:
        raise ValueError("No overlapping barcodes between prediction and ground truth")

    gt_aligned = gt.set_index("barcode").loc[common].reset_index()
    pred_aligned = pred.set_index('barcode').loc[common].reset_index()

    mismatch = ((pred_aligned['pxl_col_in_fullres'] != gt_aligned['pxl_col_in_fullres']) |
                (pred_aligned['pxl_row_in_fullres'] != gt_aligned['pxl_row_in_fullres']))
    if mismatch.any():
        warnings.warn(
            f"{int(mismatch.sum())} rows have mismatched arrow_col/arrow_row between prediction and ground truth"
        )

    return pred_aligned.reset_index(drop=True), gt_aligned.reset_index(drop=True)

def _infer_protein_cols(df):
    return [c for c in df.columns if c not in required_id_cols]

def _select_proteins(pred, gt):
    prot_pred = _infer_protein_cols(pred)
    prot_gt = _infer_protein_cols(gt)

    if protein_whitelist:
        proteins = [p for p in protein_whitelist if (p in prot_pred and p in prot_gt)]
        missing = [p for p in protein_whitelist if p not in proteins]

        if missing:
            warnings.warn(f"Missing proteins: {missing}")
    else:
        proteins = sorted(list(set(prot_pred).intersection(set(prot_gt))))
    
    if not proteins:
        raise ValueError("No overlapping proteins between prediction and ground truth")

    if protein_whitelist and (len(proteins) != len(protein_whitelist)):
        warnings.warn(
            f"Detected {len(proteins)} proteins for scoring, expected {len(protein_whitelist)}"
        )

    return proteins

def _write_scores(output_dir, results, name):
    os.makedirs(output_dir, exist_ok=True)

    scores = dict()
    per_prot_sorted = sorted(results["per_protein_scc"].items(), key=lambda kv: kv[0])
    
    scores['scc_mean'] = results['top10_mean_scc']

    for p, v in per_prot_sorted:
        scores["scc_"+p.lower()] = v 

    print("scores:", scores)
    with open(os.path.join(output_dir, f'{name}_scores.json'), 'w') as score_file:
        score_file.write(json.dumps(scores))

def score_submission(pred_csv, gt_csv):
    print("Start scoring")
    pred = pd.read_csv(pred_csv)
    pred.columns = [c.strip() for c in pred.columns]

    gt = pd.read_csv(gt_csv)
    gt.columns = [c.strip() for c in gt.columns]

    print("Start alignment")
    pred_aligned, gt_aligned = _align_on_barcode(pred, gt)
    print(pred_aligned.head(), gt_aligned.head())
    proteins = _select_proteins(pred_aligned, gt_aligned)

    print("Computing SCC")
    per_protein_scc = {}
    for p in proteins:
        scc = spearmanr(pred_aligned[p], gt_aligned[p])
        per_protein_scc[p] = None if np.isnan(scc[0]) else float(scc[0])
        print(p, scc)
    
    def _score_for_sort(v):
        return -np.inf if v is None else v
    
    sorted_prots_by_scc = sorted(per_protein_scc.items(), key=lambda kv: _score_for_sort(kv[1]), reverse=True)
    top10 = [kv[0] for kv in sorted_prots_by_scc[:10]]

    top10_vals = [per_protein_scc[p] for p in top10 if per_protein_scc[p] is not None]
    if not top10_vals:
        raise ValueError("All top-10 SCCs are NaN")
    top10_mean_scc = float(np.mean(top10_vals))
    print("TOP 10 SCC Mean:", top10_mean_scc)

    results = {
        'score': top10_mean_scc,
        "top10_mean_scc": top10_mean_scc,
        "num_proteins": len(proteins),
        "n_rows_scored": int(len(pred_aligned)),
        "top10_proteins": top10,
        "per_protein_scc": per_protein_scc
    }

    return results

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", type=str, required=True)
    parser.add_argument("--gt_path", type=str, required=True)
    parser.add_argument("--pred_path", type=str, required=True)
    args = parser.parse_args()

    pred_csv = args.pred_path
    gt_csv = args.gt_path 

    results = score_submission(pred_csv, gt_csv)
    _write_scores('.', results, args.name)
    print(f"Scoring done")

if __name__ == '__main__':
    main()