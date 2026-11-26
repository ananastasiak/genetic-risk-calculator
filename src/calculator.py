import logging
from dataclasses import dataclass
from decimal import ROUND_UP, Decimal

from src.data_loader import (
    ContributionKey,
    Sample,
    Threshold,
    Trait,
)

logger = logging.getLogger(__name__)


@dataclass
class RiskResult:
    sample_id: str
    trait_code: str
    score: Decimal | None
    level: str | None
    conclusion: str
    valid_inputs: int
    total_inputs: int


def round_score(score: Decimal) -> Decimal:
    bounded_score = max(
        Decimal(-1),
        min(Decimal(1), score),
    )

    return bounded_score.quantize(
        Decimal("0.01"),
        rounding=ROUND_UP,
    )


def select_threshold(
    score: Decimal,
    thresholds: list[Threshold],
) -> Threshold:
    return max(
        (threshold for threshold in thresholds if threshold.lower_bound <= score),
        key=lambda threshold: threshold.lower_bound,
    )


def calculate_risk(
    sample: Sample,
    trait: Trait,
    contributions: dict[ContributionKey, Decimal],
    thresholds: list[Threshold],
) -> RiskResult:
    valid_contributions: list[Decimal] = []

    for rsid in trait.variant_ids:
        genotype = sample.genotypes.get(rsid)

        if genotype is None:
            continue

        key = (trait.code, rsid, genotype)
        contribution = contributions.get(key)

        if contribution is None:
            logger.warning(
                "Missing or ambiguous contribution for sample %s: %s",
                sample.sample_id,
                key,
            )
            continue

        valid_contributions.append(contribution)

    valid_inputs = len(valid_contributions)
    total_inputs = len(trait.variant_ids)

    if valid_inputs < trait.min_inputs:
        return RiskResult(
            sample_id=sample.sample_id,
            trait_code=trait.code,
            score=None,
            level=None,
            conclusion="Недостаточно данных",
            valid_inputs=valid_inputs,
            total_inputs=total_inputs,
        )

    mean_score = sum(valid_contributions, Decimal(0)) / valid_inputs

    score = round_score(mean_score)
    threshold = select_threshold(score, thresholds)

    return RiskResult(
        sample_id=sample.sample_id,
        trait_code=trait.code,
        score=score,
        level=threshold.level,
        conclusion=threshold.conclusion,
        valid_inputs=valid_inputs,
        total_inputs=total_inputs,
    )


def calculate_all_risks(
    samples: dict[str, Sample],
    traits: dict[str, Trait],
    contributions: dict[ContributionKey, Decimal],
    thresholds: dict[str, list[Threshold]],
) -> list[RiskResult]:
    results: list[RiskResult] = []

    for sample in samples.values():
        for trait in traits.values():
            result = calculate_risk(
                sample=sample,
                trait=trait,
                contributions=contributions,
                thresholds=thresholds[trait.code],
            )

            results.append(result)

    return results
