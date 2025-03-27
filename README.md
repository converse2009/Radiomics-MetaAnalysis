Radiomics Meta-Analysis Code
Description
This repository contains the R scripts used for performing meta-analyses on radiomics studies. The analysis includes sensitivity, specificity, area under the curve (AUC), and other performance metrics from various studies. Additionally, the repository includes complementary analysis steps, such as Hartung-Knapp-Sidik-Jonkman (HKSJ) method.

Contents
Meta-analysis with subgroups: Code for performing meta-analysis of sensitivity and specificity using metaprop.

AUC Calculation: Bootstrapping method for AUC estimation.

SROC (Summary Receiver Operating Characteristic): Code for plotting the sROC curve.

Deeks' Funnel Plot: Code to check publication bias using Deeks’ funnel plot method.

Leave-one-out analysis: A method to assess the influence of each study on the overall results.

Threshold effect analysis: Spearman correlation test to explore the threshold effect.

Meta-regression: Regression analysis to explore covariates that may influence the diagnostic accuracy.
