import logging
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from pathlib import Path

from openpyxl import load_workbook

logger = logging.getLogger(__name__)

ContributionKey = tuple[str, str, str]


@dataclass
class Variant:
    gene: str
    rsid: str
    chromosome: str
    alleles: set[str]


@dataclass
class Trait:
    code: str
    name: str
    variant_ids: list[str]
    min_inputs: int


@dataclass
class Threshold:
    lower_bound: Decimal
    level: str
    conclusion: str


@dataclass
class Sample:
    sample_id: str
    sex: str
    genotypes: dict[str, str | None]


def normalize_genotype(genotype: str) -> str:
    alleles = genotype.strip().upper().split("/")
    alleles = [allele.strip() for allele in alleles]

    return "/".join(sorted(alleles))


def validate_genotype(
    genotype: str,
    variant: Variant,
    sex: str,
) -> bool:
    alleles = genotype.split("/")

    if not set(alleles).issubset(variant.alleles):
        return False

    expected_ploidy = 1 if variant.chromosome == "X" and sex == "male" else 2

    return len(alleles) == expected_ploidy


def parse_genotype(
    value: object,
    variant: Variant,
    sex: str,
    sample_id: str,
) -> str | None:
    if value is None:
        return None

    genotype = normalize_genotype(str(value))

    if not validate_genotype(genotype, variant, sex):
        logger.warning(
            "Invalid genotype for sample %s, variant %s (%s): %r",
            sample_id,
            variant.rsid,
            variant.gene,
            value,
        )

        return None

    return genotype


def read_sheet(
    file_path: Path,
    sheet_name: str,
    max_col: int,
) -> tuple[tuple, list[tuple]]:
    workbook = load_workbook(
        file_path,
        read_only=True,
        data_only=True,
    )

    try:
        if sheet_name not in workbook.sheetnames:
            raise ValueError(f"Missing worksheet: {sheet_name}")

        sheet = workbook[sheet_name]

        headers = next(
            sheet.iter_rows(
                min_row=1,
                max_row=1,
                max_col=max_col,
                values_only=True,
            ),
            (),
        )

        rows = [
            row
            for row in sheet.iter_rows(
                min_row=2,
                max_col=max_col,
                values_only=True,
            )
            if any(value is not None for value in row)
        ]

        return headers, rows

    finally:
        workbook.close()


def validate_headers(
    headers: tuple,
    expected_headers: tuple[str, ...],
    sheet_name: str,
) -> None:
    if headers != expected_headers:
        raise ValueError(
            f"Unexpected headers in '{sheet_name}': "
            f"{headers}. Expected: {expected_headers}"
        )


def load_variants(file_path: Path) -> dict[str, Variant]:
    headers, rows = read_sheet(
        file_path,
        sheet_name="Варианты",
        max_col=4,
    )

    validate_headers(
        headers,
        ("Ген", "rsID", "Хромосома", "Аллели"),
        "Варианты",
    )

    variants: dict[str, Variant] = {}

    for row in rows:
        if any(value is None for value in row):
            raise ValueError(f"Incomplete variant record: {row}")

        gene, rsid, chromosome, alleles = row

        variant = Variant(
            gene=str(gene).strip(),
            rsid=str(rsid).strip(),
            chromosome=str(chromosome).strip().upper(),
            alleles={allele.strip().upper() for allele in str(alleles).split("/")},
        )

        if not variant.alleles.issubset({"A", "C", "G", "T"}):
            raise ValueError(f"Invalid alleles for {variant.rsid}: {variant.alleles}")

        if variant.rsid in variants:
            raise ValueError(f"Duplicate rsID: {variant.rsid}")

        variants[variant.rsid] = variant

    return variants


def load_traits(
    file_path: Path,
    variants: dict[str, Variant],
) -> dict[str, Trait]:
    headers, rows = read_sheet(
        file_path,
        sheet_name="Признаки",
        max_col=4,
    )

    validate_headers(
        headers,
        ("Код признака", "Название", "Варианты", "Минимум входов"),
        "Признаки",
    )

    traits: dict[str, Trait] = {}

    for row in rows:
        if any(value is None for value in row):
            raise ValueError(f"Incomplete trait record: {row}")

        code, name, variant_ids, min_inputs = row

        if type(min_inputs) is not int:
            raise ValueError(f"Invalid minimum inputs for {code}: {min_inputs}")

        trait = Trait(
            code=str(code).strip(),
            name=str(name).strip(),
            variant_ids=[rsid.strip() for rsid in str(variant_ids).split(",")],
            min_inputs=min_inputs,
        )

        if trait.code in traits:
            raise ValueError(f"Duplicate trait: {trait.code}")

        if len(trait.variant_ids) != len(set(trait.variant_ids)):
            raise ValueError(f"Duplicate variant IDs in {trait.code}")

        unknown_ids = set(trait.variant_ids) - variants.keys()

        if unknown_ids:
            raise ValueError(f"Unknown variants in {trait.code}: {unknown_ids}")

        if not 1 <= trait.min_inputs <= len(trait.variant_ids):
            raise ValueError(f"Invalid minimum inputs for {trait.code}")

        traits[trait.code] = trait

    return traits


def load_contributions(
    file_path: Path,
    variants: dict[str, Variant],
    traits: dict[str, Trait],
) -> dict[ContributionKey, Decimal]:
    headers, rows = read_sheet(
        file_path,
        sheet_name="Вклады",
        max_col=5,
    )

    validate_headers(
        headers,
        ("Код признака", "Ген", "rsID", "Генотип", "Вклад"),
        "Вклады",
    )

    contributions: dict[ContributionKey, Decimal] = {}
    conflicts: set[ContributionKey] = set()

    for row in rows:
        if any(value is None for value in row):
            raise ValueError(f"Incomplete contribution record: {row}")

        trait_code, gene, rsid, genotype, value = row

        trait_code = str(trait_code).strip()
        gene = str(gene).strip()
        rsid = str(rsid).strip()
        genotype = normalize_genotype(str(genotype))

        if trait_code not in traits:
            raise ValueError(f"Unknown trait: {trait_code}")

        if rsid not in variants:
            raise ValueError(f"Unknown variant: {rsid}")

        if rsid not in traits[trait_code].variant_ids:
            raise ValueError(f"Variant {rsid} is not linked to {trait_code}")

        variant = variants[rsid]

        if gene != variant.gene:
            raise ValueError(f"Gene mismatch for {rsid}")

        alleles = genotype.split("/")

        if not set(alleles).issubset(variant.alleles):
            raise ValueError(f"Invalid genotype {genotype} for {rsid}")

        allowed_ploidy = {1, 2} if variant.chromosome == "X" else {2}

        if len(alleles) not in allowed_ploidy:
            raise ValueError(f"Invalid genotype ploidy for {rsid}: {genotype}")

        try:
            contribution = Decimal(str(value))
        except InvalidOperation as exc:
            raise ValueError(f"Invalid contribution: {value}") from exc

        if not contribution.is_finite():
            raise ValueError(f"Non-finite contribution: {value}")

        key = (trait_code, rsid, genotype)

        if key in conflicts:
            continue

        if key in contributions:
            if contributions[key] != contribution:
                logger.warning(
                    "Conflicting contributions for %s: %s and %s",
                    key,
                    contributions[key],
                    contribution,
                )

                del contributions[key]
                conflicts.add(key)

            continue

        contributions[key] = contribution

    return contributions


def load_thresholds(
    file_path: Path,
    traits: dict[str, Trait],
) -> dict[str, list[Threshold]]:
    headers, rows = read_sheet(
        file_path,
        sheet_name="Пороги",
        max_col=4,
    )

    validate_headers(
        headers,
        ("Код признака", "Нижняя граница", "Уровень", "Заключение"),
        "Пороги",
    )

    thresholds: dict[str, list[Threshold]] = {code: [] for code in traits}

    for row in rows:
        if any(value is None for value in row):
            raise ValueError(f"Incomplete threshold record: {row}")

        trait_code, lower_bound, level, conclusion = row

        trait_code = str(trait_code).strip()

        if trait_code not in traits:
            raise ValueError(f"Unknown trait: {trait_code}")

        try:
            lower_bound = Decimal(str(lower_bound))
        except InvalidOperation as exc:
            raise ValueError(f"Invalid threshold: {lower_bound}") from exc

        if not lower_bound.is_finite():
            raise ValueError(f"Non-finite threshold: {lower_bound}")

        if not Decimal(-1) <= lower_bound <= Decimal(1):
            raise ValueError(f"Threshold out of range: {lower_bound}")

        threshold = Threshold(
            lower_bound=lower_bound,
            level=str(level).strip(),
            conclusion=str(conclusion).strip(),
        )

        if not threshold.level or not threshold.conclusion:
            raise ValueError(f"Empty level or conclusion for {trait_code}")

        thresholds[trait_code].append(threshold)

    for trait_code, items in thresholds.items():
        items.sort(key=lambda item: item.lower_bound)

        if not items:
            raise ValueError(f"Missing thresholds for {trait_code}")

        if items[0].lower_bound != Decimal(-1):
            raise ValueError(f"Thresholds must start at -1 for {trait_code}")

        bounds = {item.lower_bound for item in items}

        if len(bounds) != len(items):
            raise ValueError(f"Duplicate thresholds for {trait_code}")

    return thresholds


def load_samples(
    file_path: Path,
    variants: dict[str, Variant],
) -> dict[str, Sample]:
    expected_columns = 2 + len(variants)

    headers, rows = read_sheet(
        file_path,
        sheet_name="Образцы",
        max_col=expected_columns,
    )

    validate_headers(
        headers[:2],
        ("Образец", "Пол"),
        "Образцы",
    )

    variant_columns: dict[int, str] = {}

    for index, header in enumerate(headers[2:], start=2):
        parts = str(header).strip().split()

        if len(parts) != 2:
            raise ValueError(f"Invalid genotype column: {header}")

        gene, rsid = parts

        if rsid not in variants:
            raise ValueError(f"Unknown variant: {rsid}")

        if gene != variants[rsid].gene:
            raise ValueError(f"Gene mismatch for {rsid}")

        if rsid in variant_columns.values():
            raise ValueError(f"Duplicate genotype column: {rsid}")

        variant_columns[index] = rsid

    missing_variants = variants.keys() - set(variant_columns.values())

    if missing_variants:
        raise ValueError(f"Missing genotype columns: {missing_variants}")

    samples: dict[str, Sample] = {}

    for row in rows:
        sample_id, sex = row[:2]

        if sample_id is None or sex is None:
            raise ValueError(f"Incomplete sample metadata: {row[:2]}")

        sample_id = str(sample_id).strip()
        sex = str(sex).strip().lower()

        if not sample_id:
            raise ValueError("Empty sample ID")

        if sample_id in samples:
            raise ValueError(f"Duplicate sample ID: {sample_id}")

        if sex not in {"male", "female"}:
            raise ValueError(f"Invalid sex for sample {sample_id}: {sex}")

        genotypes: dict[str, str | None] = {}

        for index, rsid in variant_columns.items():
            genotypes[rsid] = parse_genotype(
                value=row[index],
                variant=variants[rsid],
                sex=sex,
                sample_id=sample_id,
            )

        samples[sample_id] = Sample(
            sample_id=sample_id,
            sex=sex,
            genotypes=genotypes,
        )

    return samples
