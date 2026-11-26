import csv
import subprocess
import sys
from pathlib import Path


def test_full_pipeline(tmp_path):
    project_root = Path(__file__).resolve().parents[1]

    input_path = project_root / "data" / "test_dataset.xlsx"
    output_path = tmp_path / "results.csv"

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "src.main",
            "--input",
            str(input_path),
            "--output",
            str(output_path),
        ],
        cwd=project_root,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert output_path.exists()

    with output_path.open(encoding="utf-8-sig", newline="") as file:
        rows = list(csv.DictReader(file, delimiter=";"))

    assert len(rows) == 15

    actual = {
        (row["sample_id"], row["trait_code"]): (
            row["score"],
            row["level"],
            row["valid_inputs"],
        )
        for row in rows
    }

    expected = {
        ("S01", "OXIDATIVE_STRESS"): ("-0.38", "medium", "4"),
        ("S01", "LIPID_TG"): ("0.13", "low", "2"),
        ("S01", "VITAMIN_D"): ("-0.34", "medium", "3"),
        ("S02", "OXIDATIVE_STRESS"): ("-0.84", "high", "3"),
        ("S02", "LIPID_TG"): ("0.00", "low", "3"),
        ("S02", "VITAMIN_D"): ("-0.75", "medium", "3"),
        ("S03", "OXIDATIVE_STRESS"): ("-0.34", "medium", "3"),
        ("S03", "LIPID_TG"): ("-0.67", "medium", "3"),
        ("S03", "VITAMIN_D"): ("0.00", "low", "2"),
        ("S04", "OXIDATIVE_STRESS"): ("-0.50", "medium", "4"),
        ("S04", "LIPID_TG"): ("0.00", "low", "3"),
        ("S04", "VITAMIN_D"): ("Недостаточно данных", "", "1"),
        ("S05", "OXIDATIVE_STRESS"): ("0.00", "low", "4"),
        ("S05", "LIPID_TG"): ("0.17", "protective", "3"),
        ("S05", "VITAMIN_D"): ("0.00", "low", "3"),
    }

    assert actual == expected

    conclusions = {
        (row["sample_id"], row["trait_code"]): row["conclusion"] for row in rows
    }

    assert conclusions[("S01", "OXIDATIVE_STRESS")] == (
        "Умеренно повышенный генетический риск нарушений антиоксидантной защиты."
    )

    assert conclusions[("S04", "VITAMIN_D")] == ("Недостаточно данных")

    assert conclusions[("S05", "LIPID_TG")] == (
        "Выявлены генотипы с защитным эффектом "
        "в отношении повышения уровня триглицеридов."
    )

    assert "Conflicting contributions" in result.stderr
    assert "Invalid genotype" in result.stderr
    assert "Missing or ambiguous contribution" in result.stderr
