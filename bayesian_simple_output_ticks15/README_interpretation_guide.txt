Simplified Bayesian analysis guide
==================================

Models fitted:
1. Scratch_bouts: negative-binomial model.
2. Ticks_remaining: beta-binomial model with 15 starting ticks.
3. Repeated Guinea_Pig measurements receive a random intercept.

How to read predicted_means_all_datasets.csv:
- expected_mean is the estimated mean outcome for each Infestation level.
- lower_95 and upper_95 are the 95% credible interval.
- For Ticks_remaining, values are ticks remaining out of 15.

How to read pairwise_differences_all_datasets.csv:
- contrast is level_2 minus level_1.
- A positive difference means level_2 has a higher expected outcome than level_1.
- A negative difference means level_2 has a lower expected outcome than level_1.
- prob_greater_0 is the posterior probability that the difference is positive.
- prob_less_0 is the posterior probability that the difference is negative.
- prob_practically_small is the probability that the difference is within the ROPE.

Main assumptions:
- Animals are independent experimental units.
- Infestation is coded correctly and ordered as intended.
- Scratch_bouts is a non-negative count.
- Ticks_remaining is a count between 0 and 15.
- Guinea_Pig identifies repeated measurements correctly.
- The posterior predictive check plots should look compatible with the observed data.
- MCMC diagnostics should show no divergences and Rhat near 1.00.

If diagnostics are poor:
- Increase adapt_delta to 0.99 and rerun.
- Increase iter to 6000 or 8000 if ESS is low.
- Check for data entry issues, outliers, or accidental duplicate rows.
- Consider whether a zero-inflated model is needed if there are many zeros.
