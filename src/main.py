import argparse
import csv
import logging
from pathlib import Path

from src.calculator import RiskResult, calculate_all_risks
from src.data_loader import (
    load_contributions,
    load_samples,
    load_thresholds,
    load_traits,
    load_variants,
)


def save_results(
    results: list[RiskResult],
    output_path: Path,
) -> None:
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with output_path.open("w", encoding="utf-8-sig", newline="") as file:
        writer = csv.writer(file, delimiter=";")

        writer.writerow(
            [
                "sample_id",
                "trait_code",
                "score",
                "level",
                "conclusion",
                "valid_inputs",
                "total_inputs",
            ]
        )

        for result in results:
            writer.writerow(
                [
                    result.sample_id,
                    result.trait_code,
                    result.score if result.score is not None else "Недостаточно данных",
                    result.level,
                    result.conclusion,
                    result.valid_inputs,
                    result.total_inputs,
                ]
            )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Calculate genetic risk scores from an Excel file."
    )

    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/test_dataset.xlsx"),
        help="Path to the input Excel file",
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/results.csv"),
        help="Path to the output CSV file",
    )

    args = parser.parse_args()

    if args.input.resolve() == args.output.resolve():
        parser.error("Input and output paths must be different.")

    logging.basicConfig(
        level=logging.WARNING,
        format="%(levelname)s: %(message)s",
    )

    variants = load_variants(args.input)
    traits = load_traits(args.input, variants)

    contributions = load_contributions(args.input, variants, traits)

    thresholds = load_thresholds(args.input, traits)
    samples = load_samples(args.input, variants)

    results = calculate_all_risks(
        samples=samples,
        traits=traits,
        contributions=contributions,
        thresholds=thresholds,
    )

    save_results(results, args.output)

    print(f"Loaded variants: {len(variants)}")
    print(f"Loaded traits: {len(traits)}")
    print(f"Loaded contributions: {len(contributions)}")
    print(f"Loaded thresholds: {sum(len(items) for items in thresholds.values())}")
    print(f"Loaded samples: {len(samples)}")
    print(f"Calculated results: {len(results)}")
    print(f"Results saved to: {args.output}")


if __name__ == "__main__":
    main()
