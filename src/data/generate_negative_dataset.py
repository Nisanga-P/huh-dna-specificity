import random
import pandas as pd
import argparse
from Bio.Align import PairwiseAligner

from pathlib import Path

"""
Performs pairwise similarity for each protein against all other possible proteins, picks a protein randomly from the bottom 20% dissimilar 
proteins and appends the 9-mer ori to generate a negative sample. If the selected 9-mer is still exactly the same, mutates one of the first 
7 bases randomly.

Saves the negative dataset in "data/processed/negative_datasets".

Example:
    uv run src/data/generate_negative_dataset.py --dataset Rep
    uv run src/data/generate_negative_dataset.py --dataset Relaxase
"""


POSITIVE_DATASET_PATH = Path("data/processed/final_datasets")
NEGATIVE_DATASET_PATH = Path("data/processed/negative_datasets")


def add_wrong_9mers(
    df: pd.DataFrame,
    *,
    nucleotide_col=None,
    protein_col=None,
    distant_fraction=0.20,
    seed=3242,
) -> pd.DataFrame:
    """Replace each row's 9-mer with one from a dissimilar protein row."""
    result = df.copy()
    cols = list(df.columns)

    nucleotide_col = nucleotide_col or cols[2]
    protein_col = protein_col or cols[3]

    rng = random.Random(seed)
    proteins = df[protein_col].fillna("").astype(str).tolist()
    original_9mers = df[nucleotide_col].tolist()
    n = len(df)

    if n < 2:
        raise ValueError("At least two rows are needed to select a donor.")

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5

    self_scores = [aligner.score(p, p) for p in proteins]
    similarities = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i + 1, n):
            score = aligner.score(proteins[i], proteins[j])
            denom = self_scores[i] + self_scores[j]
            similarity = 2 * score / denom if denom else 0.0
            similarities[i][j] = similarities[j][i] = similarity

    bases = "ACGT"
    replacements = []

    for i in range(n):
        candidates = [j for j in range(n) if j != i]
        candidates.sort(key=lambda j: similarities[i][j])

        k = max(1, round(len(candidates) * distant_fraction))
        donor = rng.choice(candidates[:k])
        selected_9mer = original_9mers[donor]

        if selected_9mer == original_9mers[i]:
            if not isinstance(selected_9mer, str) or len(selected_9mer) != 9:
                raise ValueError(
                    f"Row {df.index[i]!r} has a matching 9-mer that is not a "
                    "9-character string."
                )

            # Mutate one of the first seven bases to a different base.
            position = rng.randrange(7)
            current_base = selected_9mer[position].upper()
            alternatives = [base for base in bases if base != current_base]
            mutated_base = rng.choice(alternatives)

            selected_9mer = (
                selected_9mer[:position]
                + mutated_base
                + selected_9mer[position + 1:]
            )

        replacements.append(selected_9mer)

    result[nucleotide_col] = replacements
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dataset",
        required=True,
        choices=["Rep", "Relaxase"],
        help="Dataset type",
    )
    parser.add_argument("--seed", type=int, default=None)
    args = parser.parse_args()

    if args.dataset == "Rep":
        input_name = "final_rep_dataset.csv"
        output_name = "neg_dataset_rep.csv"
    else:
        input_name = "final_relax_dataset.csv"
        output_name = "neg_dataset_relax.csv"

    input_path = POSITIVE_DATASET_PATH / input_name
    output_path = NEGATIVE_DATASET_PATH / output_name

    df = pd.read_csv(input_path)
    negative_df = add_wrong_9mers(df, seed=args.seed)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    negative_df.to_csv(output_path, index=False)
    print(f"Saved negative dataset to {output_path}")


# can only be executed as a standalone program
if __name__ == "__main__":
    main()