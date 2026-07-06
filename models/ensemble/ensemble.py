#!/usr/bin/env python3
"""
Consensus-based unsupervised ensemble strategy.

Pipeline:
1. auto-discover available team prediction CSVs and align them to the
   ground-truth file by `cell_id` or `barcode`
2. z-score each team prediction per cell across proteins
3. build z-scored pairwise per-protein SCC agreement matrices between teams
4. consensus-weighted mean: for every protein, softmax-weight each team by its
   mean cross-team SCC agreement, then average the team predictions
5. spatial smoothing: average each cell prediction over its k nearest spatial
   neighbours (KDTree on coordinates)
6. evaluate per-protein SCC, all-marker mean SCC, and top-10 mean SCC using the
   dataset's reference scoring.py

All input/output paths (prediction dir, ground-truth CSV, scoring script, and
output root) are supplied via command-line arguments; nothing is hardcoded.

Example:
    python run_zscore_consensus_weighted_mean_spatial_smooth.py \\
        --prediction-dir /path/to/team_predictions \\
        --truth-path /path/to/ground_truth.csv \\
        --scoring-script /path/to/scoring.py \\
        --output-root /path/to/outputs \\
        --dataset-name my_dataset \\
        --consensus-temperature 1.0 \\
        --spatial-k 5
"""

from __future__ import annotations

import argparse
import contextlib
import importlib.util
import io
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


METHOD_NAME = "zscore_consensus_weighted_mean_spatial_smooth"


@dataclass(frozen=True)
class DatasetConfig:
    name: str
    prediction_dir: Path
    truth_path: Path
    scoring_script_path: Path


ID_CANDIDATES = ["cell_id", "barcode"]
COORD_CANDIDATES = [
    ("x_centroid", "y_centroid"),
    ("pxl_row_in_fullres", "pxl_col_in_fullres"),
]
MIN_ID_OVERLAP_RATIO = 0.5


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=f"Run the {METHOD_NAME} unsupervised ensemble baseline."
    )
    parser.add_argument(
        "--prediction-dir",
        type=Path,
        required=True,
        help="Directory containing per-team prediction CSVs.",
    )
    parser.add_argument(
        "--truth-path",
        type=Path,
        required=True,
        help="Path to the ground-truth CSV.",
    )
    parser.add_argument(
        "--scoring-script",
        type=Path,
        required=True,
        help="Path to the reference scoring.py used for SCC evaluation.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        required=True,
        help="Directory for ensemble predictions and evaluation summaries.",
    )
    parser.add_argument(
        "--dataset-name",
        type=str,
        default=None,
        help="Name used for the output subdirectory. Defaults to the prediction-dir name.",
    )
    parser.add_argument(
        "--consensus-temperature",
        type=float,
        default=1.0,
        help="Softmax temperature for consensus-weighted averaging.",
    )
    parser.add_argument(
        "--spatial-k",
        type=int,
        default=5,
        help="k for spatial smoothing.",
    )
    return parser.parse_args()


def canonicalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [str(col).strip() for col in df.columns]
    return df


def infer_id_column(columns: list[str]) -> str:
    for col in ID_CANDIDATES:
        if col in columns:
            return col
    raise ValueError(f"Cannot infer ID column from columns: {columns[:10]} ...")


def infer_coord_columns(columns: list[str]) -> tuple[str, str] | None:
    for x_col, y_col in COORD_CANDIDATES:
        if x_col in columns and y_col in columns:
            return x_col, y_col
    return None


def extract_team_name(path: Path) -> str:
    import re

    match = re.search(r"(team\d+)", path.stem)
    if match:
        return match.group(1)
    return path.stem


def read_csv(path: Path) -> pd.DataFrame:
    return canonicalize_columns(pd.read_csv(path, low_memory=False))


def discover_team_paths(prediction_dir: Path) -> list[Path]:
    team_paths = sorted(prediction_dir.glob("*.csv"), key=lambda p: extract_team_name(p))
    if not team_paths:
        raise FileNotFoundError(f"No CSV files found in {prediction_dir}")
    return team_paths


def load_scoring_module(scoring_script_path: Path):
    spec = importlib.util.spec_from_file_location(
        f"scoring_{scoring_script_path.stem}", scoring_script_path
    )
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load scoring module from {scoring_script_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SCORING_MODULE_CACHE: dict[Path, object] = {}


def get_scoring_module(scoring_script_path: Path):
    scoring_script_path = scoring_script_path.resolve()
    if scoring_script_path not in SCORING_MODULE_CACHE:
        SCORING_MODULE_CACHE[scoring_script_path] = load_scoring_module(scoring_script_path)
    return SCORING_MODULE_CACHE[scoring_script_path]


def filter_team_dfs_by_id_overlap(
    truth_df: pd.DataFrame,
    team_paths: list[Path],
    team_dfs: list[pd.DataFrame],
    id_col: str,
) -> tuple[list[Path], list[pd.DataFrame]]:
    truth_ids = set(truth_df[id_col].astype(str))
    kept_paths: list[Path] = []
    kept_dfs: list[pd.DataFrame] = []

    for path, df in zip(team_paths, team_dfs):
        if id_col not in df.columns:
            print(f"[Skip] {path.name}: missing ID column `{id_col}`")
            continue
        pred_ids = set(df[id_col].astype(str))
        overlap = len(truth_ids & pred_ids)
        overlap_ratio = overlap / max(len(truth_ids), 1)
        if overlap_ratio < MIN_ID_OVERLAP_RATIO:
            print(
                f"[Skip] {path.name}: ID overlap {overlap}/{len(truth_ids)} "
                f"({overlap_ratio:.2%}) is too low for this dataset."
            )
            continue
        kept_paths.append(path)
        kept_dfs.append(df)

    if not kept_paths:
        raise ValueError("No valid team prediction files remain after ID-overlap filtering.")
    return kept_paths, kept_dfs


def prepare_truth_and_predictions(
    truth_path: Path, team_paths: list[Path]
) -> tuple[pd.DataFrame, list[pd.DataFrame], list[Path], str, list[str], list[str]]:
    truth_df = read_csv(truth_path)
    id_col = infer_id_column(list(truth_df.columns))
    coord_cols = infer_coord_columns(list(truth_df.columns))

    team_dfs = [read_csv(path) for path in team_paths]
    team_paths, team_dfs = filter_team_dfs_by_id_overlap(truth_df, team_paths, team_dfs, id_col)
    common_columns = set(truth_df.columns)
    for df in team_dfs:
        common_columns &= set(df.columns)

    excluded = {id_col}
    if coord_cols is not None:
        excluded.update(coord_cols)
    protein_cols = [col for col in truth_df.columns if col in common_columns and col not in excluded]
    if not protein_cols:
        raise ValueError(f"No common protein columns found for {truth_path}")

    base_columns = [id_col]
    if coord_cols is not None:
        base_columns.extend(coord_cols)
    ordered_columns = base_columns + protein_cols

    truth_df = truth_df[ordered_columns].copy()
    truth_df[id_col] = truth_df[id_col].astype(str)
    truth_df = truth_df.drop_duplicates(subset=id_col, keep="first").set_index(id_col)

    team_subsets = []
    for path, df in zip(team_paths, team_dfs):
        subset = df[[id_col] + protein_cols].copy()
        subset[id_col] = subset[id_col].astype(str)
        subset = subset.drop_duplicates(subset=id_col, keep="first").set_index(id_col)
        team_subsets.append((path, subset))

    common_ids = truth_df.index
    for _path, subset in team_subsets:
        common_ids = common_ids[common_ids.isin(subset.index)]
    if len(common_ids) == 0:
        raise ValueError("No shared IDs remain after intersecting truth with all kept team files.")
    if len(common_ids) < len(truth_df):
        print(
            f"[Info] Restricted evaluation to {len(common_ids)} shared IDs "
            f"(dropped {len(truth_df) - len(common_ids)} rows due to partial team coverage)."
        )

    truth_df = truth_df.loc[common_ids]
    aligned_team_dfs = []
    for _path, subset in team_subsets:
        aligned_team_dfs.append(subset.loc[common_ids].reset_index())

    truth_df = truth_df.reset_index()
    return truth_df, aligned_team_dfs, team_paths, id_col, list(coord_cols or []), protein_cols


def build_team_array(team_dfs: list[pd.DataFrame], protein_cols: list[str]) -> np.ndarray:
    arrays = []
    for df in team_dfs:
        arr = df[protein_cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
        arrays.append(arr)
    return np.stack(arrays, axis=0)


def zscore_team_array(team_array: np.ndarray) -> np.ndarray:
    mean = np.nanmean(team_array, axis=1, keepdims=True)
    std = np.nanstd(team_array, axis=1, keepdims=True)
    std = np.where((~np.isfinite(std)) | (std == 0), 1.0, std)
    return (team_array - mean) / std


def build_neighbor_index(coords: np.ndarray, k: int) -> np.ndarray:
    coords = np.asarray(coords, dtype=float)
    n = coords.shape[0]
    if n == 0:
        raise ValueError("Cannot build neighbor index for empty coordinates.")
    k_eff = min(max(int(k), 1), n)
    tree = cKDTree(coords)
    _, nn_idx = tree.query(coords, k=k_eff)
    if np.ndim(nn_idx) == 1:
        nn_idx = np.asarray(nn_idx)[:, None]
    return np.asarray(nn_idx, dtype=int)


def spatial_smooth_prediction(pred: np.ndarray, nn_idx: np.ndarray) -> np.ndarray:
    smoothed = np.nanmean(pred[nn_idx], axis=1)
    bad = np.isnan(smoothed).all(axis=1)
    if bad.any():
        smoothed[bad] = pred[bad]
    return smoothed


def to_scoring_frame(
    df: pd.DataFrame,
    id_col: str,
    coord_cols: list[str],
    protein_cols: list[str],
    scoring_module,
) -> pd.DataFrame:
    if len(coord_cols) != 2:
        raise ValueError("The reference scoring script requires exactly two coordinate columns.")
    required_id_cols = list(getattr(scoring_module, "required_id_cols", []))
    if len(required_id_cols) != 3:
        raise ValueError("Scoring script must define exactly three required_id_cols.")
    out = pd.DataFrame(
        {
            required_id_cols[0]: df[id_col].astype(str),
            required_id_cols[1]: df[coord_cols[0]],
            required_id_cols[2]: df[coord_cols[1]],
        }
    )
    for protein in protein_cols:
        out[protein] = df[protein]
    return out


def score_with_reference_scoring(pred_csv: Path, gt_csv: Path, scoring_module) -> dict:
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        results = scoring_module.score_submission(str(pred_csv), str(gt_csv))
    return results


def build_pairwise_protein_scc_matrices(
    team_input_paths: list[Path],
    protein_cols: list[str],
    scoring_module,
) -> dict[str, np.ndarray]:
    n_teams = len(team_input_paths)
    matrices = {protein: np.eye(n_teams, dtype=float) for protein in protein_cols}
    for i in range(n_teams):
        for j in range(i + 1, n_teams):
            results = score_with_reference_scoring(
                team_input_paths[i], team_input_paths[j], scoring_module
            )
            pair_scores = results["per_protein_scc"]
            for protein in protein_cols:
                score = pair_scores.get(protein)
                score = 0.0 if score is None or not np.isfinite(score) else float(score)
                matrices[protein][i, j] = score
                matrices[protein][j, i] = score
    return matrices


def softmax_weights(values: np.ndarray, temperature: float) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    if x.size == 0:
        return x
    shifted = x - np.nanmax(x)
    exp_x = np.exp(shifted / max(float(temperature), 1e-8))
    denom = np.sum(exp_x)
    if not np.isfinite(denom) or denom <= 0:
        return np.full_like(exp_x, 1.0 / len(exp_x))
    return exp_x / denom


def ensemble_consensus_weighted_mean(
    team_array: np.ndarray,
    team_names: list[str],
    protein_cols: list[str],
    pairwise_protein_scc: dict[str, np.ndarray],
    method_name: str,
    temperature: float,
) -> tuple[np.ndarray, pd.DataFrame]:
    n_teams, n_cells, n_proteins = team_array.shape
    output = np.full((n_cells, n_proteins), np.nan, dtype=float)
    weight_rows: list[dict[str, float | str]] = []

    for protein_idx, protein in enumerate(protein_cols):
        corr = pairwise_protein_scc[protein]
        team_scores = np.zeros(n_teams, dtype=float)
        for team_idx in range(n_teams):
            peer_scores = np.delete(corr[team_idx], team_idx)
            team_scores[team_idx] = float(np.mean(peer_scores)) if len(peer_scores) > 0 else 0.0
        weights = softmax_weights(team_scores, temperature=temperature)
        output[:, protein_idx] = np.tensordot(weights, team_array[:, :, protein_idx], axes=(0, 0))

        row = {
            "method": method_name,
            "protein": protein,
            "temperature": float(temperature),
        }
        for team_idx, team_name in enumerate(team_names):
            row[f"consensus_score_{team_name}"] = float(team_scores[team_idx])
            row[f"weight_{team_name}"] = float(weights[team_idx])
        weight_rows.append(row)

    return output, pd.DataFrame(weight_rows)


def prediction_frame(
    truth_df: pd.DataFrame,
    id_col: str,
    coord_cols: list[str],
    protein_cols: list[str],
    pred: np.ndarray,
) -> pd.DataFrame:
    out = truth_df[[id_col] + coord_cols].copy()
    out[protein_cols] = pred
    return out[[id_col] + coord_cols + protein_cols]


def evaluate_method(
    method_name: str,
    pred_path: Path,
    truth_path: Path,
    scoring_module,
) -> tuple[dict[str, float | str], pd.DataFrame]:
    results = score_with_reference_scoring(pred_path, truth_path, scoring_module)
    per_protein = results["per_protein_scc"]
    scc_values = [score for score in per_protein.values() if score is not None and np.isfinite(score)]
    if not scc_values:
        raise ValueError(f"No valid SCC values returned by scoring.py for {pred_path}")
    all_scc = float(np.mean(scc_values))
    top10_scc = float(results["top10_mean_scc"])
    metric_row = {
        "method": method_name,
        "all_markers_mean_scc": all_scc,
        "top10_mean_scc": top10_scc,
        "score": float(results["score"]),
        "n_rows_scored": int(results["n_rows_scored"]),
        "num_proteins": int(results["num_proteins"]),
        "top10_proteins": ",".join(results["top10_proteins"]),
    }
    per_marker = pd.DataFrame(
        [
            {
                "method": method_name,
                "protein": protein,
                "scc": np.nan if score is None else float(score),
            }
            for protein, score in per_protein.items()
        ]
    )
    return metric_row, per_marker


def run_dataset(
    config: DatasetConfig,
    output_root: Path,
    consensus_temperature: float,
    spatial_k: int,
) -> dict[str, float | str]:
    print(f"\n{'=' * 100}")
    print(f"Dataset: {config.name}")
    print(f"{'=' * 100}")

    team_paths = discover_team_paths(config.prediction_dir)
    discovered_team_names = [extract_team_name(path) for path in team_paths]
    print(f"Found {len(discovered_team_names)} team prediction files: {discovered_team_names}")
    scoring_module = get_scoring_module(config.scoring_script_path)

    truth_df, aligned_team_dfs, team_paths, id_col, coord_cols, protein_cols = prepare_truth_and_predictions(
        truth_path=config.truth_path,
        team_paths=team_paths,
    )
    team_names = [extract_team_name(path) for path in team_paths]
    print(f"Aligned {len(aligned_team_dfs)} teams on {len(truth_df)} rows and {len(protein_cols)} proteins.")

    raw_team_array = build_team_array(aligned_team_dfs, protein_cols)
    zscore_team_pred_array = zscore_team_array(raw_team_array)

    dataset_output_dir = output_root / config.name
    dataset_output_dir.mkdir(parents=True, exist_ok=True)
    scoring_input_dir = dataset_output_dir / "scoring_inputs"
    scoring_input_dir.mkdir(parents=True, exist_ok=True)

    truth_scoring_df = to_scoring_frame(truth_df, id_col, coord_cols, protein_cols, scoring_module)
    truth_scoring_path = scoring_input_dir / "ground_truth.csv"
    truth_scoring_df.to_csv(truth_scoring_path, index=False)

    coords = truth_df[coord_cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    nn_idx = build_neighbor_index(coords, k=spatial_k)

    # Build z-scored per-team scoring inputs, then z-scored pairwise SCC matrices.
    zscore_team_scoring_paths: list[Path] = []
    for team_idx, team_name in enumerate(team_names):
        zscore_team_df = prediction_frame(
            truth_df=truth_df,
            id_col=id_col,
            coord_cols=coord_cols,
            protein_cols=protein_cols,
            pred=zscore_team_pred_array[team_idx],
        )
        zscore_team_scoring_df = to_scoring_frame(
            zscore_team_df, id_col, coord_cols, protein_cols, scoring_module
        )
        zscore_team_path = scoring_input_dir / f"{team_name}_zscore.csv"
        zscore_team_scoring_df.to_csv(zscore_team_path, index=False)
        zscore_team_scoring_paths.append(zscore_team_path)

    zscore_pairwise_protein_scc = build_pairwise_protein_scc_matrices(
        zscore_team_scoring_paths, protein_cols, scoring_module
    )

    # zscore_consensus_weighted_mean_spatial_smooth:
    #   consensus-weighted mean on z-scored predictions, then spatial smoothing.
    pred, selection_df = ensemble_consensus_weighted_mean(
        team_array=zscore_team_pred_array,
        team_names=team_names,
        protein_cols=protein_cols,
        pairwise_protein_scc=zscore_pairwise_protein_scc,
        method_name=METHOD_NAME,
        temperature=consensus_temperature,
    )
    pred = spatial_smooth_prediction(pred, nn_idx)

    pred_df = prediction_frame(
        truth_df=truth_df,
        id_col=id_col,
        coord_cols=coord_cols,
        protein_cols=protein_cols,
        pred=pred,
    )
    pred_path = dataset_output_dir / f"{METHOD_NAME}_prediction.csv"
    pred_df.to_csv(pred_path, index=False)

    pred_scoring_df = to_scoring_frame(pred_df, id_col, coord_cols, protein_cols, scoring_module)
    pred_scoring_path = scoring_input_dir / f"{METHOD_NAME}_prediction.csv"
    pred_scoring_df.to_csv(pred_scoring_path, index=False)

    metric_row, per_marker_df = evaluate_method(
        METHOD_NAME, pred_scoring_path, truth_scoring_path, scoring_module
    )
    metric_row.update(
        {
            "dataset": config.name,
            "n_teams": len(team_names),
            "teams": ",".join(team_names),
            "n_rows": int(len(truth_df)),
            "n_proteins": len(protein_cols),
        }
    )

    if not selection_df.empty:
        selection_path = dataset_output_dir / f"{METHOD_NAME}_selection.csv"
        selection_df.to_csv(selection_path, index=False)
        metric_row["selection_path"] = str(selection_path)

    per_marker_path = dataset_output_dir / f"{METHOD_NAME}_per_marker_scc.csv"
    per_marker_df.to_csv(per_marker_path, index=False)

    print(
        f"{METHOD_NAME:<45} all_markers_mean_scc={metric_row['all_markers_mean_scc']:.4f}  "
        f"top10_mean_scc={metric_row['top10_mean_scc']:.4f}"
    )
    print(f"Saved prediction to {pred_path}")
    print(f"Saved per-marker SCC table to {per_marker_path}")
    return metric_row


def main() -> None:
    args = parse_args()
    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    prediction_dir = args.prediction_dir.resolve()
    dataset_name = args.dataset_name or prediction_dir.name
    config = DatasetConfig(
        name=dataset_name,
        prediction_dir=prediction_dir,
        truth_path=args.truth_path.resolve(),
        scoring_script_path=args.scoring_script.resolve(),
    )

    metric_row = run_dataset(
        config,
        output_root,
        consensus_temperature=args.consensus_temperature,
        spatial_k=args.spatial_k,
    )

    summary_path = output_root / config.name / f"{METHOD_NAME}_summary.csv"
    pd.DataFrame([metric_row]).to_csv(summary_path, index=False)
    print(f"\nSaved summary to {summary_path}")


if __name__ == "__main__":
    main()
