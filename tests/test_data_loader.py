from decimal import Decimal

import pytest
from openpyxl import Workbook

from src.data_loader import (
    Trait,
    Variant,
    load_contributions,
    load_samples,
    load_thresholds,
    normalize_genotype,
    parse_genotype,
    validate_genotype,
)


@pytest.mark.parametrize(
    "genotype, expected",
    [
        ("C/T", "C/T"),
        ("T/C", "C/T"),
        ("G/A", "A/G"),
        ("c/c", "C/C"),
        ("T", "T"),
        (" t / c ", "C/T"),
    ],
)
def test_normalize_genotype(genotype, expected):
    assert normalize_genotype(genotype) == expected


@pytest.fixture
def variants():
    return {
        "rs4880": Variant(
            gene="SOD2",
            rsid="rs4880",
            chromosome="6",
            alleles={"C", "T"},
        )
    }


@pytest.fixture
def traits():
    return {
        "OXIDATIVE_STRESS": Trait(
            code="OXIDATIVE_STRESS",
            name="Оксидативный стресс",
            variant_ids=["rs4880"],
            min_inputs=1,
        )
    }


def create_contributions_file(tmp_path, rows):
    workbook = Workbook()
    sheet = workbook.active
    sheet.title = "Вклады"

    sheet.append(
        [
            "Код признака",
            "Ген",
            "rsID",
            "Генотип",
            "Вклад",
        ]
    )

    for row in rows:
        sheet.append(row)

    file_path = tmp_path / "test_contributions.xlsx"

    workbook.save(file_path)
    workbook.close()

    return file_path


@pytest.mark.parametrize(
    "values, expected",
    [
        ([0], Decimal(0)),
        ([-0.5], Decimal("-0.5")),
        ([-0.5, -0.5], Decimal("-0.5")),
        ([-0.5, -0.4], None),
        ([-0.5, -0.4, -0.5], None),
    ],
)
def test_contribution_duplicates(tmp_path, variants, traits, values, expected):
    rows = [("OXIDATIVE_STRESS", "SOD2", "rs4880", "C/T", value) for value in values]

    file_path = create_contributions_file(tmp_path, rows)

    contributions = load_contributions(file_path, variants, traits)

    key = ("OXIDATIVE_STRESS", "rs4880", "C/T")

    if expected is None:
        assert key not in contributions
    else:
        assert contributions == {key: expected}


def test_contribution_genotype_normalization(tmp_path, variants, traits):
    rows = [
        ("OXIDATIVE_STRESS", "SOD2", "rs4880", "C/T", -0.5),
        ("OXIDATIVE_STRESS", "SOD2", "rs4880", "T/C", -0.5),
    ]

    file_path = create_contributions_file(tmp_path, rows)

    contributions = load_contributions(file_path, variants, traits)

    assert contributions == {("OXIDATIVE_STRESS", "rs4880", "C/T"): Decimal("-0.5")}


@pytest.mark.parametrize(
    "rsid, genotype, error_message",
    [
        ("rs4880", "A/G", "Invalid genotype"),
        ("rs4880", "C/C/T", "Invalid genotype ploidy"),
        ("rs4880", "C", "Invalid genotype ploidy"),
        ("rs999999", "C/T", "Unknown variant"),
    ],
)
def test_invalid_contributions(
    tmp_path, variants, traits, rsid, genotype, error_message
):
    rows = [("OXIDATIVE_STRESS", "SOD2", rsid, genotype, -0.5)]

    file_path = create_contributions_file(tmp_path, rows)

    with pytest.raises(ValueError, match=error_message):
        load_contributions(file_path, variants, traits)


def create_thresholds_file(tmp_path, rows):
    workbook = Workbook()
    sheet = workbook.active
    sheet.title = "Пороги"

    sheet.append(
        [
            "Код признака",
            "Нижняя граница",
            "Уровень",
            "Заключение",
        ]
    )

    for row in rows:
        sheet.append(row)

    file_path = tmp_path / "test_thresholds.xlsx"

    workbook.save(file_path)
    workbook.close()

    return file_path


def test_threshold_sorting(tmp_path, traits):
    rows = [
        ("OXIDATIVE_STRESS", 0.15, "protective", "Protective effect"),
        ("OXIDATIVE_STRESS", -0.75, "medium", "Medium risk"),
        ("OXIDATIVE_STRESS", -1, "high", "High risk"),
        ("OXIDATIVE_STRESS", -0.15, "low", "Low risk"),
    ]

    file_path = create_thresholds_file(tmp_path, rows)

    thresholds = load_thresholds(file_path, traits)
    items = thresholds["OXIDATIVE_STRESS"]

    assert [item.lower_bound for item in items] == [
        Decimal(-1),
        Decimal("-0.75"),
        Decimal("-0.15"),
        Decimal("0.15"),
    ]

    assert [item.level for item in items] == [
        "high",
        "medium",
        "low",
        "protective",
    ]

    assert items[0].conclusion == "High risk"


@pytest.mark.parametrize(
    "bounds, error_message",
    [
        ([-1, -1], "Duplicate thresholds"),
        ([-0.75, -0.15], "Thresholds must start at -1"),
        ([-1.1], "Threshold out of range"),
        ([-1, 1.1], "Threshold out of range"),
        ([], "Missing thresholds"),
    ],
)
def test_invalid_thresholds(tmp_path, traits, bounds, error_message):
    rows = [("OXIDATIVE_STRESS", bound, "level", "Conclusion") for bound in bounds]

    file_path = create_thresholds_file(tmp_path, rows)

    with pytest.raises(ValueError, match=error_message):
        load_thresholds(file_path, traits)


def test_threshold_unknown_trait(tmp_path, traits):
    rows = [
        ("UNKNOWN_TRAIT", -1, "high", "High risk"),
    ]

    file_path = create_thresholds_file(tmp_path, rows)

    with pytest.raises(ValueError, match="Unknown trait"):
        load_thresholds(file_path, traits)


@pytest.mark.parametrize(
    "chromosome, sex, genotype, expected",
    [
        ("6", "female", "C/C", True),
        ("6", "male", "C/T", True),
        ("6", "male", "A/G", False),
        ("6", "male", "T", False),
        ("6", "female", "C/C/T", False),
        ("X", "male", "T", True),
        ("X", "male", "C/T", False),
        ("X", "female", "C/T", True),
        ("X", "female", "T", False),
        ("X", "male", "A", False),
    ],
)
def test_validate_genotype(chromosome, sex, genotype, expected):
    variant = Variant(
        gene="TEST",
        rsid="rs123",
        chromosome=chromosome,
        alleles={"C", "T"},
    )

    assert validate_genotype(genotype, variant, sex) == expected


@pytest.mark.parametrize(
    "value, expected",
    [
        ("T/C", "C/T"),
        (" t / c ", "C/T"),
        ("C/C", "C/C"),
        ("A/G", None),
        ("ПРОВЕРИТЬ", None),
        ("", None),
        (None, None),
    ],
)
def test_parse_genotype(variants, value, expected):
    result = parse_genotype(
        value=value,
        variant=variants["rs4880"],
        sex="female",
        sample_id="TEST_SAMPLE",
    )

    assert result == expected


def test_invalid_genotype_warning(variants, caplog):
    result = parse_genotype(
        value="A/G",
        variant=variants["rs4880"],
        sex="female",
        sample_id="TEST_SAMPLE",
    )

    assert result is None
    assert "Invalid genotype" in caplog.text
    assert "TEST_SAMPLE" in caplog.text
    assert "rs4880" in caplog.text


def test_load_samples(tmp_path, variants):
    variants["rs1050828"] = Variant(
        gene="G6PD",
        rsid="rs1050828",
        chromosome="X",
        alleles={"C", "T"},
    )

    workbook = Workbook()
    sheet = workbook.active
    sheet.title = "Образцы"

    sheet.append(
        [
            "Образец",
            "Пол",
            "G6PD rs1050828",
            "SOD2 rs4880",
        ]
    )

    sheet.append(["S01", "female", "C/T", "T/C"])
    sheet.append(["S02", "male", "T", "A/G"])
    sheet.append(["S03", "male", "C/T", None])

    file_path = tmp_path / "test_samples.xlsx"

    workbook.save(file_path)
    workbook.close()

    samples = load_samples(file_path, variants)

    assert len(samples) == 3

    assert samples["S01"].genotypes == {
        "rs4880": "C/T",
        "rs1050828": "C/T",
    }

    assert samples["S02"].genotypes == {
        "rs4880": None,
        "rs1050828": "T",
    }

    assert samples["S03"].genotypes == {
        "rs4880": None,
        "rs1050828": None,
    }
