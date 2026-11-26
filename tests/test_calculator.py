from decimal import Decimal

import pytest

from src.calculator import (
    calculate_risk,
    round_score,
    select_threshold,
)
from src.data_loader import (
    Sample,
    Threshold,
    Trait,
)


@pytest.mark.parametrize(
    "score, expected",
    [
        ("0.33333", "0.34"),
        ("-0.33333", "-0.34"),
        ("0.125", "0.13"),
        ("-0.125", "-0.13"),
        ("-0.375", "-0.38"),
        ("0.50", "0.50"),
        ("-0.50", "-0.50"),
        ("0", "0.00"),
        ("0.001", "0.01"),
        ("-0.001", "-0.01"),
        ("1.5", "1.00"),
        ("-1.5", "-1.00"),
        ("1", "1.00"),
        ("-1", "-1.00"),
    ],
)
def test_round_score(score, expected):
    result = round_score(Decimal(score))

    assert result == Decimal(expected)
    assert result.as_tuple().exponent == -2


@pytest.fixture
def thresholds():
    return [
        Threshold(Decimal(-1), "high", "High risk"),
        Threshold(Decimal("-0.75"), "medium", "Medium risk"),
        Threshold(Decimal("-0.15"), "low", "Low risk"),
        Threshold(Decimal("0.15"), "protective", "Protective effect"),
    ]


@pytest.mark.parametrize(
    "score, expected_level",
    [
        ("-1", "high"),
        ("-0.76", "high"),
        ("-0.75", "medium"),
        ("-0.38", "medium"),
        ("-0.16", "medium"),
        ("-0.15", "low"),
        ("0", "low"),
        ("0.14", "low"),
        ("0.15", "protective"),
        ("1", "protective"),
    ],
)
def test_select_threshold(thresholds, score, expected_level):
    result = select_threshold(Decimal(score), thresholds)

    assert result.level == expected_level


@pytest.fixture
def trait():
    return Trait(
        code="TEST_TRAIT",
        name="Test trait",
        variant_ids=["rs1", "rs2", "rs3", "rs4"],
        min_inputs=3,
    )


@pytest.fixture
def contributions():
    return {
        ("TEST_TRAIT", "rs1", "C/T"): Decimal("-0.5"),
        ("TEST_TRAIT", "rs2", "C/C"): Decimal(0),
        ("TEST_TRAIT", "rs3", "T/T"): Decimal(-1),
        ("TEST_TRAIT", "rs4", "C/C"): Decimal(0),
    }


def test_calculate_risk(trait, contributions, thresholds):
    sample = Sample(
        sample_id="S01",
        sex="female",
        genotypes={
            "rs1": "C/T",
            "rs2": "C/C",
            "rs3": "T/T",
            "rs4": "C/C",
        },
    )

    result = calculate_risk(sample, trait, contributions, thresholds)

    assert result.sample_id == "S01"
    assert result.trait_code == "TEST_TRAIT"
    assert result.score == Decimal("-0.38")
    assert result.level == "medium"
    assert result.conclusion == "Medium risk"
    assert result.valid_inputs == 4
    assert result.total_inputs == 4


def test_calculate_risk_minimum_inputs(trait, contributions, thresholds):
    sample = Sample(
        sample_id="S02",
        sex="female",
        genotypes={
            "rs1": "C/T",
            "rs2": "C/C",
            "rs3": "T/T",
            "rs4": None,
        },
    )

    result = calculate_risk(sample, trait, contributions, thresholds)

    assert result.score == Decimal("-0.50")
    assert result.level == "medium"
    assert result.valid_inputs == 3
    assert result.total_inputs == 4


def test_calculate_risk_insufficient_data(trait, contributions, thresholds):
    sample = Sample(
        sample_id="S03",
        sex="female",
        genotypes={
            "rs1": "C/T",
            "rs2": "C/C",
            "rs3": None,
            "rs4": None,
        },
    )

    result = calculate_risk(sample, trait, contributions, thresholds)

    assert result.score is None
    assert result.level is None
    assert result.conclusion == "Недостаточно данных"
    assert result.valid_inputs == 2
    assert result.total_inputs == 4


def test_calculate_risk_missing_contribution(trait, contributions, thresholds, caplog):
    sample = Sample(
        sample_id="S04",
        sex="female",
        genotypes={
            "rs1": "T/T",
            "rs2": "C/C",
            "rs3": "T/T",
            "rs4": "C/C",
        },
    )

    result = calculate_risk(sample, trait, contributions, thresholds)

    assert result.score == Decimal("-0.34")
    assert result.level == "medium"
    assert result.valid_inputs == 3
    assert result.total_inputs == 4

    assert "Missing or ambiguous contribution" in caplog.text
