# Causal DAG Specification

## Research Question

What is the causal effect of habitual weekday sleep timing on glycemic health among U.S. adults?

## Target Population

U.S. adults aged >=18 years participating in the NHANES Mobile Examination Center sample.

Development cohort:
- NHANES 2015-2016

Temporal test cohort:
- NHANES 2017-2018

## Primary Exposure

Habitual weekday sleep midpoint.

Sleep midpoint is reconstructed from:

- SLQ300: usual weekday sleep time
- SLQ310: usual weekday wake time

The midpoint is represented as a circular time variable using:

- midpoint_sin
- midpoint_cos

Sleep duration is retained as a separate variable.

## Primary Outcome

HbA1c:

- LBXGH

HbA1c is treated as a continuous glycemic-health outcome.

## Core Causal Assumption

The primary estimand is the total causal effect of habitual weekday sleep timing on HbA1c.

Conceptually:

Sleep timing -> HbA1c

Sleep duration is treated separately because sleep timing and sleep duration may jointly influence glycemic health.

## Candidate Confounders

The primary adjustment set will be restricted to variables considered plausibly causally prior to habitual sleep timing and HbA1c.

Candidate variables:

- RIDAGEYR: age
- RIAGENDR: sex
- RIDRETH3: race/ethnicity
- DMDEDUC2: educational attainment
- INDFMPIR: family income-to-poverty ratio
- DMDMARTL: marital status

These variables are candidates rather than automatically required adjustment variables.

## Potential Effect Modifiers

The following variables may be evaluated for heterogeneous effects:

- age
- sex
- socioeconomic status

Interaction/effect-heterogeneity analyses will be conducted only after the primary causal effect has been estimated.

## Potential Mediators

The following variables may plausibly lie downstream of sleep timing:

- physical activity
- BMI
- waist circumference
- blood pressure
- depressive symptoms
- metabolic disease

These variables will NOT automatically be included in the primary adjustment set.

Adjusting for a mediator would estimate a controlled/direct effect rather than the total effect and could introduce overadjustment bias.

## Variables Requiring Special Caution

### Physical activity

Potential structure:

Sleep timing -> physical activity -> HbA1c

Therefore physical activity may be a mediator.

### BMI

Potential structure:

Sleep timing -> adiposity -> HbA1c

Therefore BMI may be a mediator.

### Blood pressure

Potential structure:

Sleep timing -> cardiometabolic dysfunction -> blood pressure

Therefore blood pressure is not automatically a confounder.

### Depression

Potential structure:

Sleep timing -> depressive symptoms -> HbA1c

Depressive symptoms may also share common causes with sleep timing.

Because of this ambiguity, DPQ_TOTAL will not be required for the primary analysis.

### Diabetes diagnosis

DIQ010 will not automatically be included in the primary adjustment set.

Diabetes diagnosis may be downstream of metabolic dysfunction and may therefore be inappropriate to condition on when HbA1c is the outcome.

## Primary DAG

A simplified causal structure is:

Age
    |
    v
Sleep timing ------> HbA1c
    ^
    |
Socioeconomic factors

Sex, race/ethnicity, education, and socioeconomic status are treated as candidate common causes of sleep timing and glycemic health.

## Primary Adjustment Set

The initial minimal adjustment set will be:

- age
- sex
- race/ethnicity
- education
- income-to-poverty ratio

Marital status will be considered as a sensitivity covariate rather than automatically included.

The final adjustment set must be justified using the DAG and not selected by statistical significance.

## Sleep Duration

Sleep duration is not treated as a conventional confounder by default.

It is closely related to sleep timing and may represent a second dimension of sleep behavior.

Primary analysis:

Sleep midpoint -> HbA1c

with sleep duration included as a prespecified adjustment/modification variable where appropriate.

Secondary analysis:

Sleep midpoint x sleep duration

to evaluate whether the relationship between timing and glycemic health differs by duration.

## Total Effect vs Direct Effect

The primary analysis targets the total effect of sleep timing.

Therefore variables that may lie downstream of sleep timing should not be included simply because they improve predictive performance.

Predictive importance is not sufficient justification for causal adjustment.

## Selection Bias

The analysis will restrict the cohort to adults with valid MEC examination status and available exposure/outcome measurements.

Potential selection effects caused by missing exposure, outcome, or covariate information will be evaluated.

Complete-case analysis will not be assumed to be unbiased.

## Measurement Bias

Sleep timing is self-reported.

Therefore the primary exposure represents habitual reported sleep timing rather than directly measured biological circadian phase.

The study must not claim to measure DLMO or biological circadian alignment.

## Circadian Interpretation

The study does not directly measure circadian phase.

Therefore:

"behavioral sleep timing"

or

"behavioral sleep-circadian timing"

is preferred over:

"biological circadian phase"

or

"measured circadian misalignment."

## Temporal Validation

NHANES 2015-2016:

- development
- preprocessing decisions
- model development
- causal estimation

NHANES 2017-2018:

- frozen temporal test
- no hyperparameter tuning
- no model selection based on test performance

The temporal test evaluates whether the estimated relationship/generalization is stable under a later survey cycle.

## Primary Estimand

The primary estimand is the causal effect of habitual weekday sleep timing on continuous HbA1c.

Because sleep timing is circular, estimation should respect the 24-hour structure.

## Sensitivity Analyses

Planned sensitivity analyses include:

1. Sleep timing represented using circular sine/cosine components.
2. Alternative sleep-duration definitions.
3. Exclusion of implausible sleep durations.
4. Alternative covariate specifications.
5. Bootstrap uncertainty intervals.
6. Robustness to missing-data assumptions.
7. Temporal validation in NHANES 2017-2018.
8. Heterogeneous effects by age, sex, and socioeconomic status.

## Prohibited Shortcuts

The analysis will not:

- select confounders solely by p-value;
- adjust for every available NHANES variable;
- use 2017-2018 for model tuning;
- define circadian misalignment solely as sleeping after midnight;
- claim biological circadian phase from self-reported sleep timing;
- optimize exclusively for predictive performance;
- interpret association as causation without causal assumptions.
